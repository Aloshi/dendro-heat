static char help[] = "Driver for heat";

#include "externVars.h"  // must be included or many linker errors will occur
#include "dendro.h"
#include "mpi.h"
#include <iostream>
#include <fstream>
#include <vector>

#include "petscksp.h"

//#include "oda.h"
//#include "TreeNode.h"
//#include "parUtils.h"
//#include "omg.h"

#include "timeInfo.h"
#include "hcurvedata.h"
#include "femUtils.h"
//#include "feMatrix.h"
//#include "feVector.h"
//#include "femUtils.h"
//#include "timeStepper.h"
//#include "newmark.h"
//#include "elasStiffness.h"
//#include "elasMass.h"
//#include "raleighDamping.h"
//#include "cardiacDynamic.h"

#include "massMatrix.h"
#include "parabolic.h"
#include "stiffnessMatrix.h"
#include "VecIO.h"
#include "rhs.h"
#include "DendroIO.h"

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
  std::vector<double> pts;
  pts.reserve(n_pts*3);
  for (int i = 0; i < n_pts * 3; i++) {
    pts.push_back(gaussian(mean, dev));
  }
  return pts;
}

std::vector<ot::TreeNode> buildOct()
{
  int n_pts = 20*20*20;
  std::vector<double> pts = genPoints(n_pts);
  double gSize[3] = {1.0, 1.0, 1.0};
  int dim = 3;
  int maxOctDepth = 8;
  int maxPtsPerOctant = 1;
  bool incCorner = true;

  std::vector<ot::TreeNode> linOct, balOct, newLinOct;

  ot::points2Octree(pts, gSize, linOct, dim, maxOctDepth, maxPtsPerOctant, MPI_COMM_WORLD);
  par::sampleSort<ot::TreeNode>(linOct, newLinOct, MPI_COMM_WORLD);
  ot::balanceOctree (newLinOct, balOct, dim, maxOctDepth, incCorner, MPI_COMM_WORLD);

  return balOct;
}

int main(int argc, char **argv)
{       
  PetscInitialize(&argc, &argv, "ht.opt", help);
  _InitializeHcurve(3);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int Ns = 32;
  unsigned int dof = 1;

  char problemName[PETSC_MAX_PATH_LEN];
  char filename[PETSC_MAX_PATH_LEN];

  double t0 = 0.0;
  double dt = 0.001;
  double t1 = 0.01;

  //DM  da;         // Underlying scalar DA - for scalar properties
  Vec rho;        // density - elemental scalar

  // Initial conditions
  Vec initialTemperature; 
  Vec elementalTemp;  // hack since we can't loop over nodes, gets interpolated onto nodes later

  timeInfo ti;

  PetscBool mf = PETSC_FALSE;
  bool mfree = false;
  PetscOptionsGetBool(0, "-mfree", &mf, 0);
  
  if (mf == PETSC_TRUE) {
    mfree = true;
  } else 
    mfree = false;

  // get Ns
  CHKERRQ ( PetscOptionsGetInt(0,"-Ns",&Ns,0) );
  CHKERRQ ( PetscOptionsGetScalar(0,"-t0",&t0,0) );
  CHKERRQ ( PetscOptionsGetScalar(0,"-t1",&t1,0) );
  CHKERRQ ( PetscOptionsGetScalar(0,"-dt",&dt,0) );
  CHKERRQ ( PetscOptionsGetString(PETSC_NULL, "-pn",problemName,PETSC_MAX_PATH_LEN-1,PETSC_NULL));

  // Time info for timestepping
  ti.start = t0;
  ti.stop  = t1;
  ti.step  = dt;

  /*if (!rank) {
    std::cout << "Grid size is " << Ns+1 << " and NT is " << (int)ceil(1.0/dt) << std::endl;
  }*/

  

  // create DA
  /*CHKERRQ ( DMDACreate3d ( PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, 
                    Ns+1, Ns+1, Ns+1, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                    1, 1, 0, 0, 0, &da) );*/
  auto oct = buildOct();
  auto octCopy = oct;
  std::cout << "Created oct tree\n";
  std::cout << "oct.size(): " << oct.size() << "\n";
  ot::DA da (oct, PETSC_COMM_WORLD, PETSC_COMM_WORLD, 0.3);
  std::cout << "Created DA\n";

  MPI_Barrier(PETSC_COMM_WORLD);
  massMatrix* Mass = new massMatrix(feMat::OCT); // Mass Matrix
  stiffnessMatrix* Stiffness = new stiffnessMatrix(feMat::OCT); // Stiffness matrix
  forceVector* Force = new forceVector(feVec::OCT);  // force term

  // create vectors 
  da.createVector(rho, true, false, 1);  // args??
  da.createVector(initialTemperature, false, true, 1);
  da.createVector(elementalTemp, true, true, 1);
  //CHKERRQ( DMCreateGlobalVector(da, &rho) );
  //CHKERRQ( DMCreateGlobalVector(da, &initialTemperature) );

  // Set initial conditions
  CHKERRQ( VecSet ( initialTemperature, 0.0) ); 
  CHKERRQ( VecSet ( rho, 1.0 ) );
  CHKERRQ( VecZeroEntries(elementalTemp) );

  PetscScalar *elemInitTempArray;
  da.vecGetBuffer(elementalTemp, elemInitTempArray, true, true, false, dof);

  int maxD = da.getMaxDepth();
  double hx = 1.0 / ((double)(1 << (maxD-1)));
  for ( da.init<ot::DA_FLAGS::ALL>(), da.init<ot::DA_FLAGS::WRITABLE>(); da.curr() < da.end<ot::DA_FLAGS::ALL>(); da.next<ot::DA_FLAGS::ALL>()) {
    unsigned int i = da.curr();
    unsigned int lev = da.getLevel(i);
    unsigned int half = ((1 << (maxD - lev))) / 2;

    Point pt = da.getCurrentOffset();
    double coords[3] = { (pt.x() + half) * hx, (pt.y() + half) * hx, (pt.z() + half) * hx };
    std::cout << "elem " << i << ": " << coords[0] << ", " << coords[1] << ", " << coords[2] << "\n";
    elemInitTempArray[i] = sin(M_PI * coords[0]) * sin(M_PI * coords[1]) * sin(M_PI * coords[2]);
  }

  da.vecRestoreBuffer(elementalTemp, elemInitTempArray, true, true, false, dof);

  elementToNode(da, elementalTemp, initialTemperature, dof);

  // print IC
{
  PetscScalar* data;
  //da.vecGetBuffer(initialTemperature, data, false, false, true, dof);
  da.vecGetBuffer(elementalTemp, data, true, false, true, dof);
  octree2VTK(da, rank, data, "ic.plt");
  //da.vecRestoreBuffer(initialTemperature, data, false, false, true, dof);
  da.vecRestoreBuffer(elementalTemp, data, true, false, true, dof);
}


  /*int x, y, z, m, n, p;
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

  double acx,acy,acz;
  double hx = 1.0/((double)Ns);

  // SET MATERIAL PROPERTIES ...

  // @todo - Write routines to read/write in Parallel
  // allocate for temporary buffers ...
  unsigned int elemSize = Ns*Ns*Ns;  // number of elements
  std::cout << "Elem size is " << elemSize << std::endl;
  unsigned int nodeSize = (Ns+1)*(Ns+1)*(Ns+1);  // number of nodes

  // Set Elemental material properties
  PetscScalar ***initialTemperatureArray, ***rhoArray;

  CHKERRQ(DMDAVecGetArray(da, initialTemperature, &initialTemperatureArray));
  CHKERRQ(DMDAVecGetArray(da, rho, &rhoArray));

  std::cout << "Setting initial guess." << std::endl;

  // loop through all nodes ...
  for (int k=z; k < z+p; k++) {
    for (int j=y; j < y+n; j++) {
      for (int i=x; i < x+m; i++) {
        double coords[3] = { hx*i, hx*j, hx*k };
        double ic = sin(M_PI * coords[0]) * sin(M_PI * coords[1]) * sin(M_PI * coords[2]);
        //std::cout << "ic at " << coords[0] << ", " << coords[1] << ", " << coords[2] << ": " << ic << "\n";

        initialTemperatureArray[k][j][i] = ic;
        rhoArray[k][j][i] = 1.0;
      } // end i
    } // end j
  } // end k

  std::cout << "Finished initial conditions loop." << std::endl;

  CHKERRQ( DMDAVecRestoreArray ( da, initialTemperature, &initialTemperatureArray ) );
  CHKERRQ( DMDAVecRestoreArray ( da, rho, &rhoArray ) );
  */

  //write_vector("ic.plt", initialTemperature, da);
 
  // DONE - SET MATERIAL PROPERTIES ...

  unsigned int numSteps = (unsigned int)(ceil(( ti.stop - ti.start)/ti.step));
  std::cout << "Numsteps is " << numSteps << std::endl;

  // Setup Matrices and Force Vector ...
  Mass->setProblemDimensions(1.0, 1.0, 1.0);
  Mass->setDA(&da);
  Mass->setDof(dof);

  Stiffness->setProblemDimensions(1.0, 1.0, 1.0);
  Stiffness->setDA(&da);
  Stiffness->setDof(dof);
  Stiffness->setNuVec(rho);

  Force->setProblemDimensions(1.0, 1.0, 1.0);
  Force->setDA(&da);
  Force->setDof(dof);

  // time stepper ...
  parabolic *ts = new parabolic; 

  ts->setMassMatrix(Mass);
  ts->setStiffnessMatrix(Stiffness);
  ts->setForceVector(Force);
  ts->setTimeFrames(1);

  ts->setInitialTemperature(initialTemperature);

  ts->setTimeInfo(&ti);
  ts->setAdjoint(false); // set if adjoint or forward
  //ts->useMatrixFree(mfree);

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

