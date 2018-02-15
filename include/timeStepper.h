#ifndef __TIME_STEPPER_H_
#define __TIME_STEPPER_H_

#include <vector>
#include <string>

#include "petscksp.h"
#include "petscdmda.h"
#include "petscts.h"
#include "feMat.h"
#include "feVec.h"
#include "timeInfo.h"
//#include "stsdamgHeader.h"
//#include "rpHeader.h"

// Boundary conditions stuff
typedef std::pair<int /*dof*/, PetscScalar /*value*/> Direchlet;

struct Boundary {
  std::vector<Direchlet> direchlets;

  inline void addDirechlet(int dof, PetscScalar value) {
    direchlets.push_back(Direchlet(dof, value));
  }
  inline bool empty() const { return direchlets.empty(); }
};

class timeStepper {
  public:
  // typedefs and enums
  enum matType {STIFFNESS, MASS, DAMPING };
  
  timeStepper();
  virtual ~timeStepper();

  // operations ...

  virtual int init() = 0;
  virtual int solve() = 0;

  int setMassMatrix(feMat* Mass);

  int setStiffnessMatrix(feMat* Stiffness);

  int setTalyMatrix(feMat* mat);

  int setForceVector(feVec* Force);

  int setTalyVector(feVec* vec);

  void setBoundaryCondition(const std::function<Boundary(double, double, double)>& f);
  
  int setAdjoint(bool flag);

  int setComputeEigenvaluesFlag(bool flag);

  int setTimeParameters(double StartTime, double StopTime, double TimeStep, double TimeStepRatio);

  int setTimeInfo(timeInfo *ti);

  inline int setDof(int dof) { m_uiDof = dof; }
  inline int getDof() const { return m_uiDof; }

  timeInfo* getTimeInfo() {
    return m_ti;
  }
  
  int setInitialTemperature(Vec initTemp);

  int setInitialVelocity(Vec InitialVelocity);

  /**
	*	@brief The Jacobian Matmult operation done a matrix free fashion
	*  @param _jac PETSC Matrix which is of shell type used in the time stepping
	*  @param _in  PETSC Vector which is the input vector
	*  @param _out PETSC Vector which is the output vector _jac*_in
	*  @return bool true if successful, false otherwise
	*
	*  See feMatrix.h for similar implementation
	**/

  virtual void  jacobianMatMult(Vec _in, Vec _out)= 0;
  virtual void jacobianGetDiagonal(Vec diag) = 0;

  virtual void  mgjacobianMatMult(DM _da, Vec _in, Vec _out)= 0;

  virtual bool  setRHSFunction(Vec _in, Vec _out) = 0;
  /**
	*	@brief The set right hand side function using the stiffness etc
	*  @param _currentSolutio PETSC Vector current solution
	*  @param _nextRHS    PETSC Vector next Right hand side
	*  @return bool true if succesfful, false otherwise
	*  See feMatrix.h for similar implementation
	**/

  virtual bool setRHS()=0;
  /*
	  {
		return asLeaf().setRHS(_currentSolution, _nextRHS);
	  }
   
   T& asLeaf() {return static_cast<T&>(*this);}
  */

  static void MatMult(Mat M, Vec In, Vec Out){
	 timeStepper *contxt;
	 MatShellGetContext(M,(void**)&contxt);

	 contxt->jacobianMatMult(In,Out);
  }

  static void MatGetDiagonal(Mat M, Vec diag){
	 timeStepper *contxt;
	 MatShellGetContext(M, (void**)&contxt);

	 contxt->jacobianGetDiagonal(diag);
  }

  static PetscErrorCode RhsFunction(TS ts, PetscReal t, Vec u, Vec F, void *ctxt){

	 timeStepper *context = (timeStepper*)ctxt;
	 context->setRHSFunction(u,F);

	 return(0);
  }

/*  static void MGMatMult(Mat M, Vec In, Vec Out){
	 stsDMMG *contxt;
	 MatShellGetContext(M,(void**)&contxt);

	 DM da = (DM)(((stsDMMG)contxt)->dm);

	 ((timeStepper*)(((stsDMMG)contxt)->user))->mgjacobianMatMult(da,In,Out);
  }*/

  
/*  static PetscErrorCode CreateJacobian(stsDMMG dmmg, Mat *J){
	 DM da = (DM)(dmmg->dm);
	 int m,n,xm,ym,zm;
	 int ierr;
	 DALocalInfo info;

	 ierr = DAGetLocalInfo(da,&info); CHKERRQ(ierr);
	 ierr = DAGetCorners(da,0,0,0,&xm,&ym,&zm); CHKERRQ(ierr);
	 m=n=xm*ym*zm;

	 std::cout << "size @ level = "<< m << std::endl;
	 ierr = MatCreateShell(PETSC_COMM_WORLD,m*info.dof,n*info.dof,PETSC_DETERMINE,PETSC_DETERMINE,dmmg,J); CHKERRQ(ierr);
	 ierr = MatShellSetOperation(*J,MATOP_MULT,(void(*)(void))MGMatMult); CHKERRQ(ierr);
								  
	 return(0);
  }

  static PetscErrorCode ComputeJacobian(stsDMMG dmmg, Mat A, Mat B){

	 return(0);
  }

  static PetscErrorCode ComputeRHS(stsDMMG dmmg, Vec b){

	 // int ierr;
	 // int size;

	 VecCopy(((timeStepper*)dmmg->user)->m_vecRHS,b);
	 return(0);
  }

  static PetscErrorCode FormInitialGuess(stsDMMG dmmg, Vec X){

	 VecCopy(((timeStepper*)dmmg->user)->m_vecSolution,X);

	 double norm;
	 VecNorm(X,NORM_INFINITY,&norm);
	 std::cout << " InitialGuess norm " << norm << std::endl;
	 return(0);
  }*/

  /**
	*	@brief return the force vector (this is for static vector)
	*  @return feVec*
	**/
  feVec* getForce()
	 {
		return m_Force;
	 }
  
  feMat* getMass()
	 {
		return m_Mass;
	 }

  feMat* getStiffness()
	 {
		return m_Stiffness;
	 }

  double getEMax()
	 {
		return m_demax;
	 }
  double getEMin()
	 {
		return m_demin;
	 }

  void updateBoundaries(ot::DA* da);
  void updateBoundaries(DM da);

 protected:
  inline const double* getProblemSize() {
    static const double problemSize[3] = { 1.0, 1.0, 1.0};
    return problemSize;
  }

  void applyMatBoundaryConditions(DM da, Mat mat);
  void applyMatBoundaryConditions(ot::DA* da, Mat mat);
  void applyMatVecBoundaryConditions(ot::DA* da, Vec in, Vec out);
  void applyVecBoundaryConditions(DM da, Vec rhs);
  void applyVecBoundaryConditions(ot::DA* da, Vec rhs);

  feMat* m_Mass;

  feMat* m_Stiffness;

  feMat* m_TalyMat;

  feVec* m_Force;

  feVec* m_TalyVec;
  
  // Time stepper parameters
  double			m_dTimeStep;
  double			m_dStartTime;
  double			m_dStopTime;

  double			m_dTimeStepRatio;

  // Initial Conditions set to null
  Vec				m_vecInitialSolution;
  Vec				m_vecInitialVelocity; 

  // Matrix 
  Mat				m_matJacobian;

  // Working Vectors
  Vec				m_vecRHS;
  Vec				m_vecSolution;
  Vec      		m_vecVelocity;
  Vec   			m_vecAccn;

  // Linear Solver
  KSP				m_ksp;

  // stsDMMG Multigrid
  //stsDMMG      *m_dmmg;
  
  // Time info
  timeInfo     *m_ti;
  
  // solve direction default is false
  bool				m_bIsAdjoint;

  // Flag to compute eigenvalues of the forward solver
  bool            m_bComputeEigenvalues;

  // Eigenvalues
  double          m_demax;
  double          m_demin;

  int m_uiDof;

  // for generating
  std::function<Boundary(double, double, double)> m_boundaryCondition;

  // applied to matrix/vector
  std::vector<PetscInt> m_boundaryRows;
  std::vector<PetscScalar> m_boundaryValues;
};

//#include "timeStepper.cpp"

#endif

/*
template <typename T>
class Newmark : public timeStepper<T> {

};

template <typename T>
class forwardEuler : public timeStepper<T> {

};

template <typename T>
class backwardEuler : public timeStepper<T> {

};

*/
