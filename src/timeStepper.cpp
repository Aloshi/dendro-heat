#include "timeStepper.h"

//static int _internal_ierr;
// #define iC(fun) {_internal_ierr = fun; CHKERRQ(_internal_ierr);}

timeStepper::timeStepper()
{
  // MASS, STIFFNESS
  m_Mass = NULL;
  m_Stiffness = NULL;
  m_TalyMat = NULL;
  m_TalyVec = NULL;

  // Set initial displacement and velocity to null
  m_vecInitialSolution = NULL;
  m_vecInitialVelocity = NULL; 

  // Matrix 
  m_matJacobian = NULL;

  m_vecRHS = NULL;
  m_vecSolution = NULL;

  // Linear Solver
  m_ksp = NULL;

  // Adjoint flag false
  m_bIsAdjoint = false;

  m_uiDof = 0;
}
timeStepper::~timeStepper()
{
}

/**
 *	@brief  set the initial displacement for the problem
 * @param  PETSC Vec InitialDisplacement
 * @return bool true if successful, false otherwise
 **/
int timeStepper::setInitialTemperature(Vec initialTemperature)
{
  m_vecInitialSolution = initialTemperature;
  return(0);
}
/**
 *	@brief set the initial velocity for the problem
 * @param PETSC Vec InitialVelocity
 * @return bool true if successful, false otherwise
 **/
int timeStepper::setInitialVelocity(Vec InitialVelocity)
{
  m_vecInitialVelocity = InitialVelocity;
  return(0);
}
/**
 *	@brief This function sets the Mass Matrix
 * @param Mass operator
 * @return bool true if successful, false otherwise
 **/
int timeStepper::setMassMatrix(feMat* Mass)
{
  m_Mass = Mass;
  return(0);
}

/**
 *	@brief This function sets the Stiffness Matrix
 * @param Stiffness operator
 * @return bool true if successful, false otherwise
 **/
int timeStepper::setStiffnessMatrix(feMat* Stiffness)
{
  m_Stiffness = Stiffness;
  return(0);
}

int timeStepper::setTalyMatrix(feMat* taly)
{
  m_TalyMat = taly;
  return(0);
}

int timeStepper::setTalyVector(feVec* taly) {
  m_TalyVec = taly;
  return 0;
}

/**
 *	@brief This function sets the Force vector
 * @param Force vector
 * @return bool true if successful, false otherwise
 **/
int timeStepper::setForceVector(feVec* Force)
{
  m_Force = Force;
  return(0);
}

/**
 *	@brief This function sets all the time parameters
 * @param StartTime double starting time of the timestepper
 * @param StopTime  double stopping time of the timestepper
 * @param TimeStep  double timestep used in the simulation
 * @param TimeStepRatio double timestep ratio used for reaction...not used now
 * @return bool true if successful, false otherwise
 **/
int timeStepper::setTimeParameters(double StartTime, double StopTime, double TimeStep, double TimeStepRatio)
{
  m_dStartTime = StartTime;
  m_dStopTime  = StopTime;
  m_dTimeStep = TimeStep;
  m_dTimeStepRatio = TimeStepRatio;

  return(0);
}
/**
 *	@brief This function sets the time Info
 * @param *ti, pointer to ti
 * @return bool true if successful, false otherwise
 **/
int timeStepper::setTimeInfo(timeInfo *ti)
{
  m_ti = ti;
  return(0);
}
/**
 *	@brief This function sets the solve direction for the adjoint
 * @param flag, bool true implies adjoint is on, by default adjoint is false
 * @return bool, true if successuful, false otherwise
 **/

int timeStepper::setAdjoint(bool flag)
{
  m_bIsAdjoint= flag;
  return(0);
}

void getBoundaries(ot::DA* da, const double* problemSize, int ndof,
                   std::vector<PetscInt>* rows_out, std::vector<PetscScalar>* values_out,
                   const std::function<Boundary(double, double, double)>& f)
{
  std::map<unsigned int, Boundary> bdy;

  unsigned int maxD = da->getMaxDepth();
  unsigned int lev;
  double hx, hy, hz, dist, half;
  Point pt;

  double xFac = problemSize[0]/((double)(1<<(maxD-1))); 
  double yFac = problemSize[1]/((double)(1<<(maxD-1)));
  double zFac = problemSize[2]/((double)(1<<(maxD-1)));
  double xx[8], yy[8], zz[8];
  unsigned int idx[8];

  for ( da->init<ot::DA_FLAGS::WRITABLE>(); da->curr() < da->end<ot::DA_FLAGS::WRITABLE>(); da->next<ot::DA_FLAGS::WRITABLE>() ) {
    // set the value
    lev = da->getLevel(da->curr());
    hx = xFac * (1 << (maxD - lev));
    hy = yFac * (1 << (maxD - lev));
    hz = zFac * (1 << (maxD - lev));

    half = hx;
    if (half > hy)
      half = hy;
    if (half > hz)
      half = hz;

    half = half * 0.5;

    pt = da->getCurrentOffset();

    //! get the correct coordinates of the nodes ...
    unsigned char hangingMask = da->getHangingNodeIndex(da->curr());

    xx[0] = pt.x() * xFac;
    yy[0] = pt.y() * yFac;
    zz[0] = pt.z() * zFac;
    xx[1] = pt.x() * xFac + hx;
    yy[1] = pt.y() * yFac;
    zz[1] = pt.z() * zFac;
    xx[2] = pt.x() * xFac;
    yy[2] = pt.y() * yFac + hy;
    zz[2] = pt.z() * zFac;
    xx[3] = pt.x() * xFac + hx;
    yy[3] = pt.y() * yFac + hy;
    zz[3] = pt.z() * zFac;

    xx[4] = pt.x() * xFac;
    yy[4] = pt.y() * yFac;
    zz[4] = pt.z() * zFac + hz;
    xx[5] = pt.x() * xFac + hx;
    yy[5] = pt.y() * yFac;
    zz[5] = pt.z() * zFac + hz;
    xx[6] = pt.x() * xFac;
    yy[6] = pt.y() * yFac + hy;
    zz[6] = pt.z() * zFac + hz;
    xx[7] = pt.x() * xFac + hx;
    yy[7] = pt.y() * yFac + hy;
    zz[7] = pt.z() * zFac + hz;

    da->getNodeIndices(idx);
    for (int i = 0; i < 8; ++i) {
      if (!(hangingMask & (1u << i))) {
        auto boundary = f(xx[i], yy[i], zz[i]);
        if (!boundary.empty()) {
          bdy[idx[i]] = boundary;
        }
      }
    }
  } // loop over elements

  if (!da->computedLocalToGlobal())
    da->computeLocalToGlobalMappings();

  assert(da->computedLocalToGlobal());

  rows_out->clear();
  values_out->clear();

  // preallocate some memory
  rows_out->reserve(bdy.size());
  values_out->reserve(bdy.size());

  const DendroIntL* localToGlobal = da->getLocalToGlobalMap();
  for (const auto& it : bdy) {
    for (const auto& direchlet : it.second.direchlets) {
      rows_out->push_back(ndof * localToGlobal[it.first] + direchlet.first);
      values_out->push_back(direchlet.second);
    }
  }
}

void timeStepper::setBoundaryCondition(const std::function<Boundary(double, double, double)>& f)
{
  m_boundaryCondition = f;
}

void timeStepper::updateBoundaries(ot::DA* da)
{
  if (m_boundaryCondition)
    getBoundaries(da, getProblemSize(), getDof(), &m_boundaryRows, &m_boundaryValues, m_boundaryCondition);
}

void timeStepper::updateBoundaries(DM da)
{
  assert(false);
  // not implemented
}

void timeStepper::applyMatBoundaryConditions(DM da, Mat mat)
{
  assert(false);
  // not implemented
}

void timeStepper::applyMatBoundaryConditions(ot::DA* da, Mat mat)
{
  int errCode = MatZeroRows(mat, m_boundaryRows.size(), m_boundaryRows.data(), 1.0, NULL, NULL);
  assert(errCode == 0);
}

void timeStepper::applyMatVecBoundaryConditions(ot::DA* da, Vec _in, Vec _out)
{
  PetscScalar* in;
  PetscScalar* out;

  da->vecGetBuffer(_in, in, false, false, true, m_uiDof);
  da->vecGetBuffer(_out, out, false, false, false, m_uiDof);

  PetscInt low, high;
  int ierr = VecGetOwnershipRange(_out, &low, &high);
  assert(ierr == 0);

  for (unsigned int i = 0; i < m_boundaryRows.size(); i++) {
    PetscInt row = m_boundaryRows[i];
    if (row >= low && row < high) {
      out[row] = in[row];
    }
  }

  da->vecRestoreBuffer(_in, in, false, false, true, m_uiDof);
  da->vecRestoreBuffer(_out, out, false, false, false, m_uiDof);
}

void timeStepper::applyVecBoundaryConditions(DM da, Vec rhs)
{
  assert(false);
  // not implemented
}

void timeStepper::applyVecBoundaryConditions(ot::DA* da, Vec rhs)
{
  int errCode = VecSetValues(rhs, m_boundaryRows.size(), m_boundaryRows.data(), m_boundaryValues.data(), INSERT_VALUES);
  assert(errCode == 0);

  errCode = VecAssemblyBegin(rhs);
  assert(errCode == 0);
  
  errCode = VecAssemblyEnd(rhs);
  assert(errCode == 0);
}

/// Jacobian matmult, setRhs will be in the derived class

