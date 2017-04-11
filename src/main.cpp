static char help[] = "Driver for heat";

#include "externVars.h"  // must be included or many linker errors will occur
#include "dendro.h"
#include "mpi.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <functional>

#include "petscksp.h"

#include "timeInfo.h"
#include "hcurvedata.h"

#include "massMatrix.h"
#include "parabolic.h"
#include "nonlinear.h"
#include "stiffnessMatrix.h"
#include "VecIO.h"
#include "rhs.h"
#include "DendroIO.h"

#include "TalyMat.h"
#include "TalyVec.h"
#include "HTEquation.h"


// defined below
void setScalarByFunction(ot::DA* da, Vec vec, int dof, std::function<double(double,double,double,int)> f);
int setScalarByFunction(DM da, int Ns, Vec vec, int dof, std::function<double(double,double,double,int)> f);

static double gSize[3] = {1, 1, 1};

float uniform() {
  return float(rand()) / RAND_MAX; // [0,1)
}

double gaussian(double mean = 0.5, double std_deviation = 0.1) {
  static double t = 0;
  double x1, x2, r;

  // reuse previous calculations
  if (t) {
    const double tmp = t;
    t = 0;
    return mean + std_deviation * tmp;
  }

  // pick randomly a point inside the unit disk
  do {
    x1 = 2 * uniform() - 1;
    x2 = 2 * uniform() - 1;
    r = x1 * x1 + x2 * x2;
  } while (r >= 1);

  // Box-Muller transform
  r = sqrt(-2.0 * log(r) / r);

  // save for next call
  t = r * x2;

  // only use one of the coordinates of a bivariate distribution
  return mean + std_deviation * r * x1;
}

std::vector<double> genPoints(int n_pts, double mean = 0.5, double dev = 0.1)
{
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  std::vector<double> pts;
  pts.reserve(n_pts*3);
  /*for (int i = 0; i < n_pts * 3; i++) {
    //pts.push_back(gaussian(mean, dev));
    pts.push_back(uniform());
  }*/

  if (mpi_rank == 0) {
    int n_pts_1d = 16;
    for (int i = 1; i < n_pts_1d; i++) {
      double x = ((double)i) / n_pts_1d;
      for (int j = 1; j < n_pts_1d; j++) {
        double y = ((double)j) / n_pts_1d;
        for (int k = 1; k < n_pts_1d; k++) {
          double z = ((double)k) / n_pts_1d;
          pts.push_back(x);
          pts.push_back(y);
          pts.push_back(z);
        }
      }
    }
  }

/*
  for (int i = 0; i < 10*10*10 * 3; i++) {
    pts.push_back(gaussian(mean, dev));
    //pts.push_back(uniform());
  }*/

/*
  int n_start = 8;
  int n_end = 24;
  for (int i = n_start; i < n_end; i++) {
    double x = ((double)i) / n_pts_1d;
    for (int j = n_start; j < n_end; j++) {
      double y = ((double)j) / n_pts_1d;
      for (int k = n_start; k < n_end; k++) {
        double z = ((double)k) / n_pts_1d;
        pts.push_back(x);
        pts.push_back(y);
        pts.push_back(z);
      }
    }
  }*/

  return pts;
}

std::vector<ot::TreeNode> buildOct()
{
  int n_pts = 8*8*8;
  std::vector<double> pts = genPoints(n_pts);
  double gSize[3] = {1.0, 1.0, 1.0};
  int dim = 3;
  int maxOctDepth = 8;
  int maxPtsPerOctant = 1;
  bool incCorner = true;

  /*std::vector<ot::TreeNode> linOct, balOct, newLinOct;

  std::cout << "before points2octree" << std::endl;
  ot::points2Octree(pts, gSize, linOct, dim, maxOctDepth, maxPtsPerOctant, MPI_COMM_WORLD);
  std::cout << "after points2octree" << std::endl;
  par::sampleSort<ot::TreeNode>(linOct, newLinOct, MPI_COMM_WORLD);
  
  MPI_Barrier(MPI_COMM_WORLD);

  std::cout << "after sampleSort" << std::endl;
  ot::balanceOctree (newLinOct, balOct, dim, maxOctDepth, incCorner, MPI_COMM_WORLD);
  std::cout << "made octree" << std::endl;

  return balOct;*/

  std::vector<ot::TreeNode> out;
  createRegularOctree(out, 4, 3, 8, MPI_COMM_WORLD);
  return out;
}

int main(int argc, char **argv)
{       
  PetscInitialize(&argc, &argv, "ht.opt", help);
  _InitializeHcurve(3);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int Ns = 16;
  int dof = 1;

  char problemName[PETSC_MAX_PATH_LEN];
  char filename[PETSC_MAX_PATH_LEN];

  double t0 = 0.0;
  double dt = 0.001;
  double t1 = 0.001;

  // Initial conditions
  Vec initialTemperature; 
  Vec rho;        // density - elemental scalar

  timeInfo ti;

  PetscBool mf = PETSC_FALSE;
  PetscOptionsGetBool(NULL, 0, "-mfree", &mf, 0);
  bool mfree = (mf == PETSC_TRUE);

  bool useOctree = true;
  {
    PetscBool temp = PETSC_TRUE;
    PetscOptionsGetBool(NULL, NULL, "-use_octree", &temp, NULL);
    useOctree = (temp == PETSC_TRUE);
  }

  // get Ns
  CHKERRQ ( PetscOptionsGetInt(NULL, 0,"-Ns",&Ns,0) );
  CHKERRQ ( PetscOptionsGetScalar(NULL, 0,"-t0",&t0,0) );
  CHKERRQ ( PetscOptionsGetScalar(NULL, 0,"-t1",&t1,0) );
  CHKERRQ ( PetscOptionsGetScalar(NULL, 0,"-dt",&dt,0) );
  CHKERRQ ( PetscOptionsGetString(NULL, PETSC_NULL, "-pn",problemName,PETSC_MAX_PATH_LEN-1,PETSC_NULL));

  // Time info for timestepping
  ti.start = t0;
  ti.stop  = t1;
  ti.step  = dt;

  ot::DA* octDA = NULL;
  DM petscDA = NULL;

  if (!useOctree) {
    // create DA
    CHKERRQ ( DMDACreate3d ( PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, 
                      Ns+1, Ns+1, Ns+1, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                      1, 1, 0, 0, 0, &petscDA) );

    // create vectors
    CHKERRQ( DMCreateGlobalVector(petscDA, &rho) );
    CHKERRQ( DMCreateGlobalVector(petscDA, &initialTemperature) );
  } else {
    // use octree
    auto oct = buildOct();  // TODO pass in Ns
    octDA = new ot::DA(oct, PETSC_COMM_WORLD, PETSC_COMM_WORLD, 0.3);

    // create vectors 
    octDA->createVector(rho, true, false, 1 /*dof*/);
    octDA->createVector(initialTemperature, false, false, 1 /*dof*/);
  }

  std::cout << "Created DA and vectors" << std::endl;
  MPI_Barrier(PETSC_COMM_WORLD);

  // set up equation
  auto matType = (octDA ? feMat::OCT : feMat::PETSC);
  auto vecType = (octDA ? feVec::OCT : feVec::PETSC);
  massMatrix* Mass = new massMatrix(matType); // Mass Matrix
  stiffnessMatrix* Stiffness = new stiffnessMatrix(matType); // Stiffness matrix
  auto talyMat = new TalyMatrix<HTEquation, HTNodeData>(matType);  // mass + stiffness matrix
  forceVector* Force = new forceVector(vecType);  // force term
  auto talyVec = new TalyVector<HTEquation, HTNodeData>(vecType);

  // Set initial conditions
  CHKERRQ( VecSet ( initialTemperature, 0.0) ); 
  CHKERRQ( VecSet ( rho, 1.0 ) );

  if (petscDA) {
    setScalarByFunction(petscDA, Ns, initialTemperature, dof, [](double x, double y, double z, int dof) {
      return sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z);
    });
  }

  if (octDA) {
    setScalarByFunction(octDA, initialTemperature, dof, [](double x, double y, double z, int dof) {
      return sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z);
    });

    double norm_inf, norm_l2;
    VecNorm(initialTemperature, NORM_INFINITY, &norm_inf);
    VecNorm(initialTemperature, NORM_2, &norm_l2);
    if (!rank)
      std::cout << "norm_inf: " << norm_inf << ", norm_l2: " << norm_l2 << std::endl;
  }

  // print IC
  if (petscDA) {
    write_vector("ic", initialTemperature, dof, petscDA);
  }

  if (octDA) {
    octree2VTK(octDA, initialTemperature, dof, "ic");
  }

  unsigned int numSteps = (unsigned int)(ceil(( ti.stop - ti.start)/ti.step));
  std::cout << "Numsteps is " << numSteps << std::endl;

  // Setup Matrices and Force Vector ...
  Mass->setProblemDimensions(1.0, 1.0, 1.0);
  Mass->setDof(dof);

  Stiffness->setProblemDimensions(1.0, 1.0, 1.0);
  Stiffness->setDof(dof);
  Stiffness->setNuVec(rho);

  talyMat->setProblemDimensions(1.0, 1.0, 1.0);
  talyMat->setDof(dof);

  Force->setProblemDimensions(1.0, 1.0, 1.0);
  Force->setDof(dof);

  talyVec->setProblemDimensions(1.0, 1.0, 1.0);
  talyVec->setDof(dof);

  if (petscDA) {
    Mass->setDA(petscDA);
    Stiffness->setDA(petscDA);
    Force->setDA(petscDA);
    talyMat->setDA(petscDA);
    talyVec->setDA(petscDA);
  }
  if (octDA) {
    Mass->setDA(octDA);
    Stiffness->setDA(octDA);
    Force->setDA(octDA);
    talyMat->setDA(octDA);
    talyVec->setDA(octDA);
  }

  // time stepper ...
  parabolic *ts = new parabolic; 

  //ts->setMassMatrix(Mass);
  //ts->setStiffnessMatrix(Stiffness);
  //ts->setForceVector(Force);

  ts->setTalyMatrix(talyMat);
  ts->setTalyVector(talyVec);
  ts->setDof(dof);

  ts->setTimeFrames(1);

  ts->setInitialTemperature(initialTemperature);
  if (petscDA) {
    ts->setDAForMonitor(petscDA);
  }
  if (octDA) {
    ts->setDAForMonitor(octDA);
  }

  ts->setTimeInfo(&ti);
  ts->setAdjoint(false); // set if adjoint or forward
  ts->setMatrixFree(mfree);
  ts->setBoundaryCondition([] (double x, double y, double z) -> Boundary {
    static const double eps = 1e-6;

    Boundary b;
    if (fabs(x) < eps || fabs(y) < eps || fabs(z) < eps || fabs(x - 1.0) < eps || fabs(y - 1.0) < eps || fabs(z - 1.0) < eps) {
      b.addDirechlet(0, 0.0);
    }
    return b;
  });

  if (!rank)
    std::cout <<"Initializing parabolic"<< std::endl;

  double itime = MPI_Wtime();
  ts->init(); // initialize IMPORTANT 
  if (!rank)
    std::cout <<"Starting parabolic Solve"<< std::endl;
  double stime = MPI_Wtime();
  ts->solve();// solve 
  double etime = MPI_Wtime();
  if (!rank)
    std::cout <<"Done parabolic"<< std::endl;
  if (!rank) {
		std::cout << "Total time for init is " << stime - itime << std::endl;
    std::cout << "Total time for solve is " << etime - stime << std::endl;
  }

  PetscFinalize();
}

// set IC for PETSc DA
int setScalarByFunction(DM da, int Ns, Vec vec, int dof, std::function<double(double,double,double,int)> f) {
  int x, y, z, m, n, p;
  int mx,my,mz, xne, yne, zne;

  CHKERRQ( DMDAGetCorners(da, &x, &y, &z, &m, &n, &p) ); 
  CHKERRQ( DMDAGetInfo(da,0, &mx, &my, &mz, 0,0,0,0,0,0,0,0,0) ); 

  if (x+m == mx) {
    xne=m-1;
  } else {
    xne=m;
  }
  if (y+n == my) {
    yne=n-1;
  } else {
    yne=n;
  }
  if (z+p == mz) {
    zne=p-1;
  } else {
    zne=p;
  }

  double hx = 1.0/((double)Ns);

  // allocate for temporary buffers ...
  unsigned int elemSize = Ns*Ns*Ns;  // number of elements
  unsigned int nodeSize = (Ns+1)*(Ns+1)*(Ns+1);  // number of nodes

  // Set Elemental material properties
  PetscScalar ***vec_data = NULL;
  CHKERRQ(DMDAVecGetArray(da, vec, &vec_data));

  // loop through all nodes ...
  for (int k=z; k < z+p; k++) {
    for (int j=y; j < y+n; j++) {
      for (int i=x; i < x+m; i++) {
        double coords[3] = { hx*i, hx*j, hx*k };
        for (int d = 0; d < dof; d++) {
          double val = f(coords[0], coords[1], coords[2], d);
          vec_data[k][j][i*dof + d] = val;
        }
      } // end i
    } // end j
  } // end k

  CHKERRQ( DMDAVecRestoreArray ( da, vec, &vec_data ) );
}

void setScalarByFunction(ot::DA* da, Vec vec, int dof, std::function<double(double,double,double,int)> f) {

  PetscScalar *_vec = NULL; 
  da->vecGetBuffer(vec, _vec, false, false, false, dof);
  
  //da->ReadFromGhostsBegin<PetscScalar>(_vec, dof);
  //da->ReadFromGhostsEnd<PetscScalar>(_vec);
		
  unsigned int maxD = da->getMaxDepth();
  unsigned int lev;
  double hx, hy, hz;
  Point pt, parPt;

  double xFac = gSize[0]/((double)(1<<(maxD-1)));
  double yFac = gSize[1]/((double)(1<<(maxD-1)));
  double zFac = gSize[2]/((double)(1<<(maxD-1)));
  double xx[8], yy[8], zz[8];
  unsigned int idx[8];

  for ( da->init<ot::DA_FLAGS::ALL>(); da->curr() < da->end<ot::DA_FLAGS::ALL>(); da->next<ot::DA_FLAGS::ALL>() ) { 
    // set the value
    lev = da->getLevel(da->curr());
    hx = xFac*(1<<(maxD - lev));
    hy = yFac*(1<<(maxD - lev));
    hz = zFac*(1<<(maxD - lev));

    pt = da->getCurrentOffset();
		
    //! get the correct coordinates of the nodes ...
    // unsigned int chNum = da->getChildNumber();
    unsigned char hangingMask = da->getHangingNodeIndex(da->curr());
    
    // if hanging, use parents, else mine. 
    xx[0] = pt.x()*xFac; yy[0] = pt.y()*yFac; zz[0] = pt.z()*zFac;
    xx[1] = pt.x()*xFac+hx; yy[1] = pt.y()*yFac; zz[1] = pt.z()*zFac;
    xx[2] = pt.x()*xFac; yy[2] = pt.y()*yFac+hy; zz[2] = pt.z()*zFac;
    xx[3] = pt.x()*xFac+hx; yy[3] = pt.y()*yFac+hy; zz[3] = pt.z()*zFac;
    
    xx[4] = pt.x()*xFac; yy[4] = pt.y()*yFac; zz[4] = pt.z()*zFac+hz;
    xx[5] = pt.x()*xFac+hx; yy[5] = pt.y()*yFac; zz[5] = pt.z()*zFac+hz;
    xx[6] = pt.x()*xFac; yy[6] = pt.y()*yFac+hy; zz[6] = pt.z()*zFac+hz;
    xx[7] = pt.x()*xFac+hx; yy[7] = pt.y()*yFac+hy; zz[7] = pt.z()*zFac +hz;
    
    da->getNodeIndices(idx);
    for (int i=0; i<8; ++i) {
      if ( ! ( hangingMask & ( 1u << i ) ) ) {
        for (int d = 0; d < dof; d++) {
          _vec[idx[i]*dof+d] =  f(xx[i], yy[i], zz[i], d); //! use correct coordinate
        }
      }
    }
    //  std::cout << "Setting value: " << idx[0] << " to " << f(xx,yy,zz) << std::endl;
  }

  da->vecRestoreBuffer(vec,  _vec, false, false, false,  dof);
}
