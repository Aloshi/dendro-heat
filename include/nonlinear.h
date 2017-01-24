#pragma once

#include "timeStepper.h"
#include <sstream>
#include "VecIO.h"
#include "rhs.h"
#include "DendroIO.h"


class nonlinear : public timeStepper
{
 public:
  nonlinear() : m_da(NULL), m_octDA(NULL), m_matrixFree(false) {}

  virtual int init();

  virtual int solve();
  /**
	*	@brief The Jacobian matmult routine by matrix free method using the given stiffness, damping and mass matrices
	*  @param m_matJacobian PETSC Matrix, this is of type matrix shell
	*  @param In  PETSC Vec, input vector
	*  @param Out PETSC Vec, Output vector
	**/
  virtual void jacobianMatMult(Vec In, Vec Out) override;
  virtual void jacobianGetDiagonal(Vec diag) {};

  virtual void mgjacobianMatMult(DM da, Vec In, Vec Out);
  /**
	*	@brief Sets the right hand side vector given the current solution
	*  @param m_vecCurrentSolution PETSC Vec, current solution
	*  @param m_vecNextRHS PETSC Vec, next right hand side
	**/
  virtual bool setRHS();

  virtual bool setRHSFunction(Vec In, Vec Out);

    // set the number of levels for multigrid
  void setLevels(int levels)
	 {
		m_inlevels = levels;
	 }

  // set time frames for monitor
  int setTimeFrames(int mon){
	 m_iMon = mon;
  }

  // monitor routine which saves solution
  int monitor();
  std::vector<Vec> getSolution()
	 {
		return m_solVector;
	 }

  inline void setDAForMonitor(DM da) {
    m_da = da;
  }

  inline void setDAForMonitor(ot::DA* da) {
    m_octDA = da;
  }

  inline void setMatrixFree(bool matrixFree) {
    m_matrixFree = matrixFree;
  }

  static PetscErrorCode FormJacobian(SNES snes, Vec sol, Mat jac, Mat /*precond_matrix*/, void* ctx) {
    nonlinear* nl = (nonlinear*) ctx;
    if (!nl->m_matrixFree) {
      MatZeroEntries(jac);
      nl->m_TalyMat->GetAssembledMatrix_new(&jac, 0, sol);
    }
    return 0;
  }
  static PetscErrorCode FormFunction(SNES snes, Vec in, Vec out, void* ctx) {
    nonlinear* nl = (nonlinear*) ctx;

    VecZeroEntries(out);
    nl->m_TalyVec->addVec_new(in, out, 1.0);
    return 0;
  }

 private:
  int m_inlevels;
  int m_iMon;
  std::vector<Vec> m_solVector;
  DM m_da;
  ot::DA* m_octDA;
  bool m_matrixFree;
  SNES m_snes;
};


/**
 *	@brief The initialization function where the Jacobian matrix,
 *         KSP context to be used at every time step are created
 * @return bool true if successful, false otherwise
 * Matrix of shell type is created which does only a matvec.
 * This requires the mass matrix, stiffness matrix, damping matrix to do the matvec			 
 **/
int nonlinear::init()
{

  int matsize;

  int ierr;
  // Allocate memory for working vectors
  ierr = VecDuplicate(m_vecInitialSolution,&m_vecSolution); CHKERRQ(ierr);
  ierr = VecDuplicate(m_vecInitialSolution,&m_vecRHS); CHKERRQ(ierr);

  ierr = VecGetLocalSize(m_vecInitialSolution,&matsize); CHKERRQ(ierr);


  // FOR NON-MATRIX-FREE
  if (!m_matrixFree) {
    if (m_da) {  // petsc
      ierr = DMCreateMatrix(m_da, &m_matJacobian); CHKERRQ(ierr);
    } else if (m_octDA) {  // octree
      m_octDA->createMatrix(m_matJacobian, MATAIJ, 1);
    } else {
      assert(false);
    }
  } else {
    // matrix-free (same for both petsc/octree)
    ierr = MatCreateShell(PETSC_COMM_WORLD,matsize,matsize,PETSC_DETERMINE,PETSC_DETERMINE,this,&m_matJacobian); CHKERRQ(ierr);
    ierr = MatShellSetOperation(m_matJacobian,MATOP_MULT,(void(*)(void))(MatMult)); CHKERRQ(ierr);

    /*ierr = MatCreate(PETSC_COMM_WORLD, &m_matJacobian); CHKERRQ(ierr);
    ierr = MatSetSizes(m_matJacobian, matsize, matsize, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
    ierr = MatSetFromOptions(m_matJacobian); CHKERRQ(ierr);
    ierr = MatSetUp(m_matJacobian); CHKERRQ(ierr);*/
  }

  // Create a KSP context to solve  @ every timestep
  /*ierr = KSPCreate(PETSC_COMM_WORLD,&m_ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(m_ksp,m_matJacobian,m_matJacobian); CHKERRQ(ierr);
  //ierr = KSPSetType(m_ksp,KSPCG); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(m_ksp); CHKERRQ(ierr);*/
  ierr = SNESCreate(PETSC_COMM_WORLD, &m_snes); CHKERRQ(ierr);
  ierr = SNESSetFunction(m_snes, m_vecRHS, FormFunction, this); CHKERRQ(ierr);  // TODO does m_vecRHS need to be initiaized before we call this?
  ierr = SNESSetJacobian(m_snes, m_matJacobian, m_matJacobian, FormJacobian, this); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(m_snes); CHKERRQ(ierr);
  

  return(0);
}

/**
 *	@brief Solve the timestepping problem
 * @return bool, return true if successful , false otherwise
 **/
int nonlinear::solve()
{
  // Time info

  //  m_ti.start = m_dStartTime;
  //  m_ti.stop = m_dStopTime;
  //  m_ti.step = m_dTimeStep;
  int ierr;	
  unsigned int NT = (int)(ceil( m_ti->stop - m_ti->start)/m_ti->step);

  // Set initial conditions
  ierr = VecCopy(m_vecInitialSolution,m_vecSolution); CHKERRQ(ierr);
  // Time stepping

  //if(!m_bIsAdjoint) // FORWARD PROBLEM
  {
    m_ti->current = m_ti->start;
    m_ti->currentstep = 0;
    //while(m_ti->current < m_ti->stop)
    //{
      m_ti->currentstep++;
      m_ti->current += m_ti->step;

      // Get the Right hand side of the ksp solve using the current solution
      //setRHS();


//#ifdef __DEBUG__
      double norm;
      ierr = VecNorm(m_vecSolution,NORM_INFINITY,&norm); CHKERRQ(ierr);
      std::cout << "norm of the solution @ " << m_ti->current  << " = " << norm << std::endl;
//#endif

      // assemble matrix if this isn't matrix-free
      /*if (!m_matrixFree) {
        MatZeroEntries(m_matJacobian);
        m_TalyMat->GetAssembledMatrix_new(&m_matJacobian, 0, m_vecSolution);
      }*/

      // for matrix-free, implicitly "assembled" since m_matJacobian is a shell matrix

      // Solve the ksp using the current rhs and non-zero initial guess
      /*ierr = KSPSetInitialGuessNonzero(m_ksp,PETSC_TRUE); CHKERRQ(ierr);
      ierr = KSPSolve(m_ksp,m_vecRHS,m_vecSolution); CHKERRQ(ierr);*/

      ierr = SNESSolve(m_snes, PETSC_NULL, m_vecSolution); CHKERRQ(ierr);

		if(m_iMon > 0){
		  monitor();
		}

#ifdef __DEBUG__
      ierr = VecNorm(m_vecSolution,NORM_INFINITY,&norm); CHKERRQ(ierr);
      std::cout << "norm of the solution @ " << m_ti->current  << " = " << norm << std::endl;
#endif
    //}
  }

  return(0);

}

/**
 *	@brief Jacobian matrix-vec product done at every timestep
 * @param In, PETSC Vec which is the current solution
 * @param Out, PETSC Vec which is Out = J*in, J is the Jacobian (here J = Mass + dt*Stiffness)
 **/
void nonlinear::jacobianMatMult(Vec In, Vec Out)
{
  assert(m_matrixFree);
  VecZeroEntries(Out); /* Clear to zeros*/

  if (m_TalyMat) {
    m_TalyMat->MatVec_new(In, Out);
    //m_TalyMat->GetAssembledMatrix_new(&M, 0, In);
  } else {
    m_Mass->MatVec(In, Out);
    m_Stiffness->MatVec(In, Out, -m_ti->step); /* -dt factor for stiffness*/
    assert(false);
  }
}

void nonlinear::mgjacobianMatMult(DM da, Vec In, Vec Out){
  assert(false);
}
/**
 *	@brief Set right hand side used in stepping at every time step
 * @return true if successful, false otherwise
 **/

bool nonlinear::setRHS()
{
  VecZeroEntries(m_vecRHS);

  if (m_TalyVec) {
    m_TalyVec->addVec_new(m_vecSolution, m_vecRHS, 1.0);
  } else {
    ((forceVector*)m_Force)->setPrevTS(m_vecSolution);
    m_Force->addVec(m_vecRHS, 1.0/*1.0 / m_ti->step*/);
  }
  return true;
}

bool nonlinear::setRHSFunction(Vec In, Vec Out)
{
  return true;
}

int nonlinear::monitor()
{
  if(fmod(m_ti->currentstep,(double)(m_iMon)) < 0.0001)
	 {
		std::cout << "current step in the monitor " << m_ti->currentstep << std::endl;
		double norm;
		int ierr;
		Vec tempSol;
		ierr = VecDuplicate(m_vecSolution,&tempSol);  CHKERRQ(ierr);
#ifdef __DEBUG__
		//		VecNorm(m_vecSolution,NORM_INFINITY,&norm);
		//		PetscPrintf(0,"solution norm b4 push back %f\n",norm);
#endif
		ierr = VecCopy(m_vecSolution,tempSol); CHKERRQ(ierr);
		m_solVector.push_back(tempSol);

    std::stringstream ss;
    ss << "timestep_" << m_ti->currentstep;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (m_da) {
      write_vector(ss.str().c_str(), m_vecSolution, m_da);
    } else {
      std::string asdf = ss.str();
      octree2VTK(m_octDA, m_vecSolution, asdf);
    }

    /*PetscViewer view;
    PetscViewerCreate(PETSC_COMM_WORLD, &view);
    PetscViewerPushFormat(view, PETSC_VIEWER_ASCII_MATLAB); 
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, ss.str().c_str(), &view);
    VecView(m_vecSolution, view);
    PetscViewerDestroy(&view);*/

#ifdef __DEBUG__
		//		VecNorm(m_solVector[m_solVector.size()-1],NORM_INFINITY,&norm);
		//		PetscPrintf(0,"solution norm after push back %f\n",norm);
#endif
	 }
  return(0);
}

