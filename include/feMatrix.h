#ifndef __FE_MATRIX_H_
#define __FE_MATRIX_H_

#include <string>
#include "odaUtils.h"
#include "feMat.h"
#include "timeInfo.h"

template <typename T>
class feMatrix : public feMat {
  public:  

    enum stdElemType {
      ST_0,ST_1,ST_2,ST_3,ST_4,ST_5,ST_6,ST_7
    };

    enum exhaustiveElemType {
      //Order: 654321
      //YZ ZX Z XY Y X
      ET_N = 0,
      ET_Y = 2,
      ET_X = 1,
      ET_XY = 3,
      ET_Z = 8,
      ET_ZY = 10,
      ET_ZX = 9,
      ET_ZXY = 11,
      ET_XY_XY = 7,
      ET_XY_ZXY = 15,
      ET_YZ_ZY = 42,
      ET_YZ_ZXY = 43,
      ET_YZ_XY_ZXY = 47,
      ET_ZX_ZX = 25,
      ET_ZX_ZXY = 27,
      ET_ZX_XY_ZXY = 31,
      ET_ZX_YZ_ZXY = 59,
      ET_ZX_YZ_XY_ZXY = 63
    };

  feMatrix();  
  feMatrix(daType da);
  ~feMatrix();

  // operations.
  void setStencil(void* stencil);

  void setName(std::string name);

  /**
   * 	@brief		The matrix-vector multiplication routine that is used by
   * 				matrix-free methods. 
   * 	@param		_in	PETSc Vec which is the input vector with whom the 
   * 				product is to be calculated.
   * 	@param		_out PETSc Vec, the output of M*_in
   * 	@return		bool true if successful, false otherwise.
   * 
   *  The matrix-vector multiplication routine that is used by matrix-free 
   * 	methods. The product is directly calculated from the elemental matrices,
   *  which are computed by the ElementalMatVec() function. Use the Assemble()
   *  function for matrix based methods.
   **/ 
  virtual bool MatVec(Vec _in, Vec _out, double scale=1.0);

  virtual bool MatGetDiagonal(Vec _diag, double scale=1.0);

  virtual bool GetAssembledMatrix(Mat *J, MatType mtype);

  /**
   * 	@brief		The elemental matrix-vector multiplication routine that is used
   *				by matrix-free methods. 
   * 	@param		_in	PETSc Vec which is the input vector with whom the 
   * 				product is to be calculated.
   * 	@param		_out PETSc Vec, the output of M*_in
   * 	@return		bool true if successful, false otherwise.
   *  @todo		Might have to change _in and _out to std. C arrays for speed.
   *
   *  The implementation for this function shall be in derived classes, based on
   * 	the problem formulation. Look at MassMatrix and StiffnessMatrix for standard
   * 	implementations. 
   **/ 

  inline bool GetElementalMatrix(int i, int j, int k, PetscScalar *mat) {
    return asLeaf().GetElementalMatrix(i, j, k, mat);  
  }

  inline bool ElementalMatVec(int i, int j, int k, PetscScalar ***in, PetscScalar ***out, double scale) {
    return asLeaf().ElementalMatVec(i,j,k,in,out,scale);  
  }

  inline bool ElementalMatGetDiagonal(int i, int j, int k, PetscScalar ***diag, double scale) {
    return asLeaf().ElementalMatGetDiagonal(i,j,k,diag,scale);  
  }

  /**
   * 	@brief		The elemental matrix-vector multiplication routine that is used
   *				by matrix-free methods. 
   * 	@param		_in	PETSc Vec which is the input vector with whom the 
   * 				product is to be calculated.
   * 	@param		_out PETSc Vec, the output of M*_in
   * 	@return		bool true if successful, false otherwise.
   *  @todo		Might have to change _in and _out to std. C arrays for speed.
   *
   *  The implementation for this function shall be in derived classes, based on
   * 	the problem formulation. Look at MassMatrix and StiffnessMatrix for standard
   * 	implementations. 
   **/ 
  inline bool GetElementalMatrix(unsigned int index, std::vector<ot::MatRecord>& records) {
    return asLeaf().GetElementalMatrix(index, records);  
  }

  inline bool ElementalMatVec(unsigned int index, PetscScalar *in, PetscScalar *out, double scale) {
    return asLeaf().ElementalMatVec(index, in, out, scale);  
  }

  inline bool ElementalMatGetDiagonal(unsigned int index, PetscScalar *diag, double scale) {
    return asLeaf().ElementalMatGetDiagonal(index, diag, scale);  
  }

  // PetscErrorCode matVec(Vec in, Vec out, timeInfo info);

  /**
 * @brief		Allows static polymorphism access to the derived class. Using the Barton Nackman trick.
 * 
 **/

  T& asLeaf() { return static_cast<T&>(*this);}  

  bool initStencils() {
    return asLeaf().initStencils();
  }

  bool preMatVec() {
    return asLeaf().preMatVec();
  }

  bool postMatVec() {
    return asLeaf().postMatVec();
  }

  void setDof(unsigned int dof) { m_uiDof = dof; }
  unsigned int getDof() { return m_uiDof; }

  void setTimeInfo(timeInfo *t) { m_time =t; }
  timeInfo* getTimeInfo() { return m_time; }

  void initOctLut();

	inline int getEtype(unsigned char hnMask, unsigned char cNum) {
		unsigned char type=0;
		if(hnMask) {
	    switch(cNum) {
	      case 0:
					type = ot::getElemType<0>(hnMask);
					break;
	      case 1:
					type = ot::getElemType<1>(hnMask);
					break;
	      case 2:
					type = ot::getElemType<2>(hnMask);
					break;
	      case 3:
					type = ot::getElemType<3>(hnMask);
					break;
	      case 4:
					type = ot::getElemType<4>(hnMask);
					break;
	      case 5:
					type = ot::getElemType<5>(hnMask);
					break;
	      case 6:
					type = ot::getElemType<6>(hnMask);
					break;
	      case 7:
					type = ot::getElemType<7>(hnMask);
					break;
	      default:
					assert(false);
	    }
		}
		return type;
  }
	
  inline PetscErrorCode alignElementAndVertices(ot::DA * da, stdElemType & sType, unsigned int* indices);
  inline PetscErrorCode mapVtxAndFlagsToOrientation(int childNum, unsigned int* indices, unsigned char & mask);
  inline PetscErrorCode reOrderIndices(unsigned char eType, unsigned int* indices);

protected:
  void *          	m_stencil;

  std::string     	m_strMatrixType;

  timeInfo		*m_time;

  unsigned int		m_uiDof;

  // Octree specific stuff ...
  unsigned char **	m_ucpLut;
};


template <typename T>
feMatrix<T>::feMatrix() {
	m_daType = PETSC;
	m_DA    = NULL;
	m_octDA   = NULL;
	m_stencil = NULL;
	m_uiDof = 1;
	m_ucpLut  = NULL;

	// initialize the stencils ...
	initStencils();
}

template <typename T>
feMatrix<T>::feMatrix(daType da) {
#ifdef __DEBUG__
	assert ( ( da == PETSC ) || ( da == OCT ) );
#endif
	m_daType = da;
	m_DA    = NULL;
	m_octDA   = NULL;
	m_stencil = NULL;
	m_ucpLut  = NULL;

	// initialize the stencils ...
	initStencils();
	if (da == OCT)
		initOctLut();
}

template <typename T>
void feMatrix<T>::initOctLut() {
	//Note: It is not symmetric.
	unsigned char tmp[8][8]={
		{0,1,2,3,4,5,6,7},
		{1,3,0,2,5,7,4,6},
		{2,0,3,1,6,4,7,5},
		{3,2,1,0,7,6,5,4},
		{4,5,0,1,6,7,2,3},
		{5,7,1,3,4,6,0,2},
		{6,4,2,0,7,5,3,1},
		{7,6,3,2,5,4,1,0}
	};

	//Is Stored in  ROW_MAJOR Format.  
	typedef unsigned char* charPtr;
	m_ucpLut = new charPtr[8];
	for (int i=0;i<8;i++) {
		m_ucpLut[i] = new unsigned char[8]; 
		for (int j=0;j<8;j++) {
			m_ucpLut[i][j] = tmp[i][j];
		}
	}
}

template <typename T>
feMatrix<T>::~feMatrix() {
}


#undef __FUNCT__
#define __FUNCT__ "feMatrix_MatGetDiagonal"
template <typename T>
bool feMatrix<T>::MatGetDiagonal(Vec _diag, double scale){
	PetscFunctionBegin;
#ifdef __DEBUG__
	assert ( ( m_daType == PETSC ) || ( m_daType == OCT ) );
#endif

	int ierr;

	// PetscScalar zero=0.0;

	if (m_daType == PETSC) {

		PetscInt x,y,z,m,n,p;
		PetscInt mx,my,mz;
		int xne,yne,zne;

		PetscScalar ***diag;
		Vec diagLocal;

		/* Get all corners*/
		if (m_DA == NULL)
			std::cerr << "Da is null" << std::endl;
		ierr = DMDAGetCorners(m_DA, &x, &y, &z, &m, &n, &p); CHKERRQ(ierr); 
		/* Get Info*/
		ierr = DMDAGetInfo(m_DA,0, &mx, &my, &mz, 0,0,0,0,0,0,0,0,0); CHKERRQ(ierr); 

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

		ierr = DMGetLocalVector(m_DA, &diagLocal); CHKERRQ(ierr);
		ierr = VecZeroEntries(diagLocal);

		// ierr = DMGlobalToLocalBegin(m_DA, _diag, INSERT_VALUES, diagLocal); CHKERRQ(ierr);
		// ierr = DMGlobalToLocalEnd(m_DA, _diag, INSERT_VALUES, diagLocal); CHKERRQ(ierr);


		ierr = DMDAVecGetArray(m_DA, diagLocal, &diag);

		// Any derived class initializations ...
		preMatVec();

		// loop through all elements ...
		for (int k=z; k<z+zne; k++) {
			for (int j=y; j<y+yne; j++) {
				for (int i=x; i<x+xne; i++) {
					ElementalMatGetDiagonal(i, j, k, diag, scale);
				} // end i
			} // end j
		} // end k

		postMatVec();

		ierr = DMDAVecRestoreArray(m_DA, diagLocal, &diag); CHKERRQ(ierr);  


		ierr = DMLocalToGlobalBegin(m_DA, diagLocal, ADD_VALUES, _diag); CHKERRQ(ierr);  
		ierr = DMLocalToGlobalEnd(m_DA, diagLocal, ADD_VALUES, _diag); CHKERRQ(ierr);  

		ierr = DMRestoreLocalVector(m_DA, &diagLocal); CHKERRQ(ierr);  


	} else {
		// loop for octree DA.
		PetscScalar *diag=NULL;

		// get Buffers ...
		//Nodal,Non-Ghosted,Read,1 dof, Get in array and get ghosts during computation
		m_octDA->vecGetBuffer(_diag, diag, false, false, false, m_uiDof);

		preMatVec();

		// loop through all elements ...
		for ( m_octDA->init<ot::DA_FLAGS::ALL>(); m_octDA->curr() < m_octDA->end<ot::DA_FLAGS::ALL>(); m_octDA->next<ot::DA_FLAGS::ALL>() ) {
			ElementalMatGetDiagonal( m_octDA->curr(), diag, scale); 
		}//end 

		postMatVec();

		// Restore Vectors ..
		m_octDA->vecRestoreBuffer(_diag, diag, false, false, false, m_uiDof);
	}

	PetscFunctionReturn(0);
}



/**
* 	@brief		The matrix-vector multiplication routine that is used by
* 				matrix-free methods. 
* 	@param		_in	PETSc Vec which is the input vector with whom the 
* 				product is to be calculated.
* 	@param		_out PETSc Vec, the output of M*_in
* 	@return		bool true if successful, false otherwise.
* 
*  The matrix-vector multiplication routine that is used by matrix-free 
* 	methods. The product is directly calculated from the elemental matrices,
*  which are computed by the ElementalMatrix() function. Use the Assemble()
*  function for matrix based methods.
**/
#undef __FUNCT__
#define __FUNCT__ "feMatrix_MatVec"
template <typename T>
bool feMatrix<T>::MatVec(Vec _in, Vec _out, double scale){
	PetscFunctionBegin;

#ifdef __DEBUG__
	assert ( ( m_daType == PETSC ) || ( m_daType == OCT ) );
#endif

	int ierr;
	// PetscScalar zero=0.0;

	if (m_daType == PETSC) {

		PetscInt x,y,z,m,n,p;
		PetscInt mx,my,mz;
		int xne,yne,zne;

		PetscScalar ***in, ***out;
		Vec inlocal, outlocal;

		/* Get all corners*/
		if (m_DA == NULL)
			std::cerr << "Da is null" << std::endl;
		ierr = DMDAGetCorners(m_DA, &x, &y, &z, &m, &n, &p); CHKERRQ(ierr); 
		/* Get Info*/
		ierr = DMDAGetInfo(m_DA,0, &mx, &my, &mz, 0,0,0,0,0,0,0,0,0); CHKERRQ(ierr); 

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

		// std::cout << x << "," << y << "," << z << " + " << xne <<","<<yne<<","<<zne<<std::endl;

		// Get the local vector so that the ghost nodes can be accessed
		ierr = DMGetLocalVector(m_DA, &inlocal); CHKERRQ(ierr);
		ierr = DMGetLocalVector(m_DA, &outlocal); CHKERRQ(ierr);
		// ierr = VecDuplicate(inlocal, &outlocal); CHKERRQ(ierr);

		ierr = DMGlobalToLocalBegin(m_DA, _in, INSERT_VALUES, inlocal); CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(m_DA, _in, INSERT_VALUES, inlocal); CHKERRQ(ierr);
		// ierr = DMGlobalToLocalBegin(m_DA, _out, INSERT_VALUES, outlocal); CHKERRQ(ierr);
		// ierr = DMGlobalToLocalEnd(m_DA, _out, INSERT_VALUES, outlocal); CHKERRQ(ierr);

		ierr = VecZeroEntries(outlocal);

		ierr = DMDAVecGetArray(m_DA, inlocal, &in);
		ierr = DMDAVecGetArray(m_DA, outlocal, &out);

		// Any derived class initializations ...
		preMatVec();

		// loop through all elements ...
		for (int k=z; k<z+zne; k++) {
			for (int j=y; j<y+yne; j++) {
				for (int i=x; i<x+xne; i++) {
					// std::cout << i <<"," << j << "," << k << std::endl;
					ElementalMatVec(i, j, k, in, out, scale);
				} // end i
			} // end j
		} // end k

		postMatVec();

		ierr = DMDAVecRestoreArray(m_DA, inlocal, &in); CHKERRQ(ierr);  
		ierr = DMDAVecRestoreArray(m_DA, outlocal, &out); CHKERRQ(ierr);  

		ierr = DMLocalToGlobalBegin(m_DA, outlocal, ADD_VALUES, _out); CHKERRQ(ierr);  
		ierr = DMLocalToGlobalEnd(m_DA, outlocal, ADD_VALUES,_out); CHKERRQ(ierr);  

		ierr = DMRestoreLocalVector(m_DA, &inlocal); CHKERRQ(ierr);  
		ierr = DMRestoreLocalVector(m_DA, &outlocal); CHKERRQ(ierr);  
		// ierr = VecDestroy(outlocal); CHKERRQ(ierr);  

	} else {
		// loop for octree DA.


		PetscScalar *out=NULL;
		PetscScalar *in=NULL; 

		// get Buffers ...
		//Nodal,Non-Ghosted,Read,1 dof, Get in array and get ghosts during computation
		m_octDA->vecGetBuffer(_in,   in, false, false, true,  m_uiDof);
		m_octDA->vecGetBuffer(_out, out, false, false, false, m_uiDof);

		// start comm for in ...
		//m_octDA->updateGhostsBegin<PetscScalar>(in, false, m_uiDof);
		// m_octDA->ReadFromGhostsBegin<PetscScalar>(in, false, m_uiDof);
		m_octDA->ReadFromGhostsBegin<PetscScalar>(in, m_uiDof);
		preMatVec();

		// Independent loop, loop through the nodes this processor owns..
		for ( m_octDA->init<ot::DA_FLAGS::INDEPENDENT>(), m_octDA->init<ot::DA_FLAGS::WRITABLE>(); m_octDA->curr() < m_octDA->end<ot::DA_FLAGS::INDEPENDENT>(); m_octDA->next<ot::DA_FLAGS::INDEPENDENT>() ) {
			ElementalMatVec( m_octDA->curr(), in, out, scale); 
		}//end INDEPENDENT

		// Wait for communication to end.
		//m_octDA->updateGhostsEnd<PetscScalar>(in);
		m_octDA->ReadFromGhostsEnd<PetscScalar>(in);

		// Dependent loop ...
		for ( m_octDA->init<ot::DA_FLAGS::DEPENDENT>(), m_octDA->init<ot::DA_FLAGS::WRITABLE>(); m_octDA->curr() < m_octDA->end<ot::DA_FLAGS::DEPENDENT>(); m_octDA->next<ot::DA_FLAGS::DEPENDENT>() ) {
			ElementalMatVec( m_octDA->curr(), in, out, scale); 
		}//end DEPENDENT

		postMatVec();

		// Restore Vectors ...
		m_octDA->vecRestoreBuffer(_in,   in, false, false, true,  m_uiDof);
		m_octDA->vecRestoreBuffer(_out, out, false, false, false, m_uiDof);

	}

	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "feMatrix_MatAssemble"
template <typename T>
bool feMatrix<T>::GetAssembledMatrix(Mat *J, MatType mtype) {
	PetscFunctionBegin;

	int ierr;

#ifdef __DEBUG__
	assert ( ( m_daType == PETSC ) || ( m_daType == OCT ) );
#endif
	if (m_daType == PETSC) {
		Mat K;

		// Petsc Part ..
		unsigned int elemMatSize = m_uiDof*8;

		PetscScalar* Ke = new PetscScalar[elemMatSize*elemMatSize];
		MatStencil *idx = new MatStencil[elemMatSize];

		PetscInt x,y,z,m,n,p;
		PetscInt mx,my,mz;
		int xne,yne,zne;

		/* Get all corners*/
		if (m_DA == NULL)
			std::cerr << "Da is null" << std::endl;
		ierr = DMDAGetCorners(m_DA, &x, &y, &z, &m, &n, &p); CHKERRQ(ierr); 
		/* Get Info*/
		ierr = DMDAGetInfo(m_DA,0, &mx, &my, &mz, 0,0,0,0,0,0,0,0,0); CHKERRQ(ierr); 

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

		// Get the matrix from the DA ...
		// DAGetMatrix(m_DA, mtype, &J); 
		//DMGetMatrix(m_DA, MATAIJ, &K);
		DMCreateMatrix(m_DA, &K);

		MatZeroEntries(K);

		preMatVec();
		// loop through all elements ...
		for (int k=z; k<z+zne; k++) {
			for (int j=y; j<y+yne; j++) {
				for (int i=x; i<x+xne; i++) {
					int idxMap[8][3]={
						{k, j, i},
						{k,j,i+1},
						{k,j+1,i},
						{k,j+1,i+1},
						{k+1, j, i},
						{k+1,j,i+1},
						{k+1,j+1,i},
						{k+1,j+1,i+1}               
					};
					for (unsigned int q=0; q<8; q++) {
						for (unsigned int dof=0; dof<m_uiDof; dof++) {
							idx[m_uiDof*q + dof].i = idxMap[q][2];
							idx[m_uiDof*q + dof].j = idxMap[q][1];
							idx[m_uiDof*q + dof].k = idxMap[q][0];
							idx[m_uiDof*q + dof].c = dof;
						}
					}
					GetElementalMatrix(i, j, k, Ke);
					// Set Values 
					// @check if rows/cols need to be interchanged.
					MatSetValuesStencil(K, elemMatSize, idx, elemMatSize, idx, Ke, ADD_VALUES);
				} // end i
			} // end j
		} // end k
		postMatVec();

		MatAssemblyBegin (K, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd   (K, MAT_FINAL_ASSEMBLY);

		*J = K;

		delete [] Ke;
		delete [] idx;

	} else {
		if(!(m_octDA->computedLocalToGlobal())) {
			m_octDA->computeLocalToGlobalMappings();
		}
		// Octree part ...
		char matType[30];
		PetscBool typeFound;
		PetscOptionsGetString(PETSC_NULL, "-fullJacMatType", matType, 30, &typeFound);
		if(!typeFound) {
			std::cout<<"I need a MatType for the full Jacobian matrix!"<<std::endl;
			MPI_Finalize();
			exit(0);		
		} 
		m_octDA->createMatrix(*J, matType, 1);
		MatZeroEntries(*J);
		std::vector<ot::MatRecord> records;

		preMatVec();

		for(m_octDA->init<ot::DA_FLAGS::WRITABLE>(); m_octDA->curr() < m_octDA->end<ot::DA_FLAGS::WRITABLE>();	m_octDA->next<ot::DA_FLAGS::WRITABLE>()) {
			GetElementalMatrix(m_octDA->curr(), records);	
			if(records.size() > 500) {
				m_octDA->setValuesInMatrix(*J, records, 1, ADD_VALUES);
			}
		}//end writable
		m_octDA->setValuesInMatrix(*J, records, 1, ADD_VALUES);

		postMatVec();

		MatAssemblyBegin(*J, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(*J, MAT_FINAL_ASSEMBLY);
	}


	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "alignElementAndVertices"
template <typename T>
PetscErrorCode feMatrix<T>::alignElementAndVertices(ot::DA * da, stdElemType & sType, unsigned int* indices) {
	PetscFunctionBegin;

	sType = ST_0;
	da->getNodeIndices(indices); 

	// not required ....
	// int rank;
	// MPI_Comm_rank(da->getComm(), &rank);

	if (da->isHanging(da->curr())) {

		int childNum = da->getChildNumber();
		Point pt = da->getCurrentOffset();   

		unsigned char hangingMask = da->getHangingNodeIndex(da->curr());    

		//Change HangingMask and indices based on childNum
		mapVtxAndFlagsToOrientation(childNum, indices, hangingMask);    

		unsigned char eType = ((126 & hangingMask)>>1);

		reOrderIndices(eType, indices);
	}//end if hangingElem.
	PetscFunctionReturn(0);
}//end function.

#undef __FUNCT__
#define __FUNCT__ "mapVtxAndFlagsToOrientation"
template <typename T>
PetscErrorCode feMatrix<T>::mapVtxAndFlagsToOrientation(int childNum, unsigned int* indices, unsigned char & mask) {
	PetscFunctionBegin;
	unsigned int tmp[8];
	unsigned char tmpFlags = 0;
	for (int i=0;i<8;i++) {
		tmp[i] = indices[m_ucpLut[childNum][i]];
		tmpFlags = ( tmpFlags | ( ( (1<<(m_ucpLut[childNum][i])) & mask ) ? (1<<i) : 0 ) );
	}
	for (int i=0;i<8;i++) {
		indices[i] = tmp[i];
	}
	mask = tmpFlags;
	PetscFunctionReturn(0);
}//end function

#undef __FUNCT__
#define __FUNCT__ "reOrderIndices"
template <typename T>
PetscErrorCode feMatrix<T>::reOrderIndices(unsigned char eType, unsigned int* indices) {
#ifdef __DEBUG_1
	std::cout << "Entering " << __func__ << std::endl;
#endif
	PetscFunctionBegin;
	unsigned int tmp;
	switch (eType) {
	case  ET_N: 
		break;
	case  ET_Y:
		break;
	case  ET_X:
		//Swap 1 & 2, Swap 5 & 6
		tmp = indices[1];
		indices[1] = indices[2];
		indices[2] = tmp;
		tmp = indices[5];
		indices[5] = indices[6];
		indices[6] = tmp;
		break;
	case  ET_XY:
		break;
	case  ET_Z:
		//Swap 2 & 4, Swap 3 & 5
		tmp = indices[2];
		indices[2] = indices[4];
		indices[4] = tmp;
		tmp = indices[3];
		indices[3] = indices[5];
		indices[5] = tmp;
		break;
	case  ET_ZY:
		//Swap 1 & 4, Swap 3 & 6
		tmp = indices[1];
		indices[1] = indices[4];
		indices[4] = tmp;
		tmp = indices[3];
		indices[3] = indices[6];
		indices[6] = tmp;
		break;
	case  ET_ZX:
		//Swap 2 & 4, Swap 3 & 5
		tmp = indices[2];
		indices[2] = indices[4];
		indices[4] = tmp;
		tmp = indices[3];
		indices[3] = indices[5];
		indices[5] = tmp;
		break;
	case  ET_ZXY:
		break;
	case  ET_XY_XY:
		break;
	case  ET_XY_ZXY:
		break;
	case  ET_YZ_ZY:
		//Swap 1 & 4, Swap 3 & 6
		tmp = indices[1];
		indices[1] = indices[4];
		indices[4] = tmp;
		tmp = indices[3];
		indices[3] = indices[6];
		indices[6] = tmp;
		break;
	case  ET_YZ_ZXY:
		//Swap 1 & 4, Swap 3 & 6
		tmp = indices[1];
		indices[1] = indices[4];
		indices[4] = tmp;
		tmp = indices[3];
		indices[3] = indices[6];
		indices[6] = tmp;
		break;
	case  ET_YZ_XY_ZXY:
		break;
	case  ET_ZX_ZX:
		//Swap 2 & 4, Swap 3 & 5
		tmp = indices[2];
		indices[2] = indices[4];
		indices[4] = tmp;
		tmp = indices[3];
		indices[3] = indices[5];
		indices[5] = tmp;
		break;
	case  ET_ZX_ZXY:
		//Swap 2 & 4, Swap 3 & 5
		tmp = indices[2];
		indices[2] = indices[4];
		indices[4] = tmp;
		tmp = indices[3];
		indices[3] = indices[5];
		indices[5] = tmp;
		break;
	case  ET_ZX_XY_ZXY:
		//Swap 1 & 2, Swap 5 & 6
		tmp = indices[1];
		indices[1] = indices[2];
		indices[2] = tmp;
		tmp = indices[5];
		indices[5] = indices[6];
		indices[6] = tmp;
		break;
	case  ET_ZX_YZ_ZXY:
		//Swap 2 & 4, Swap 3 & 5
		tmp = indices[2];
		indices[2] = indices[4];
		indices[4] = tmp;
		tmp = indices[3];
		indices[3] = indices[5];
		indices[5] = tmp;
		break;
	case  ET_ZX_YZ_XY_ZXY:
		break;
	default:
		std::cout<<"in reOrder Etype: "<< (int) eType << std::endl;
		assert(false);
	}
#ifdef __DEBUG_1
	std::cout << "Leaving " << __func__ << std::endl;
#endif
	PetscFunctionReturn(0);
}



#endif
