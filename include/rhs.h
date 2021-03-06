#pragma once

#include <iostream>
#include <assert.h>
#include "feVector.h"

class forceVector : public feVector<forceVector>
{
  public:
  forceVector(daType da);

  inline bool ElementalAddVec(unsigned int index, PetscScalar *in, double scale);
  inline bool ElementalAddVec(PetscScalar*, PetscScalar*, PetscScalar* coords, double scale) {
    assert(false);
  }

  inline bool ElementalAddVec(int i, int j, int k, PetscScalar ***in, double scale);
    inline bool initStencils();

    bool preAddVec();
    bool postAddVec();

    inline void setPrevTS(Vec ts) {
      m_prevTSVec = ts;
    }

   private:
    Vec m_prevTSVec;
    void* m_prevTSArray;

    double 		m_dHx;
    double xFac, yFac, zFac;
    unsigned int maxD;
};


forceVector::forceVector(daType da) {
#ifdef __DEBUG__
  assert ( ( da == PETSC ) || ( da == OCT ) );
#endif
  m_daType = da;
  m_DA 		= NULL;
  m_octDA 	= NULL;
  m_stencil	= NULL;
  m_prevTSArray = NULL;

  // initialize the stencils ...
  initStencils();
  if (da == OCT)
    initOctLut();
}

bool forceVector::initStencils() {
  typedef int* int1Ptr;
  typedef int** int2Ptr;
 
  if (m_daType == PETSC) {
    int Bjk[8][8] =    {
      { 64, 32, 32, 16, 32, 16, 16,  8},
      { 32, 64, 16, 32, 16, 32,  8, 16},
      { 32, 16, 64, 32, 16,  8, 32, 16},
      { 16, 32, 32, 64,  8, 16, 16, 32},
      { 32, 16, 16,  8, 64, 32, 32, 16},
      { 16, 32,  8, 16, 32, 64, 16, 32},
      { 16,  8, 32, 16, 32, 16, 64, 32},
      {  8, 16, 16, 32, 16, 32, 32, 64}
    };

    int** Ajk = new int1Ptr[8];
    for (int j=0;j<8;j++) {
      Ajk[j] = new int[8];
      for (int k=0;k<8;k++) {
        Ajk[j][k] = Bjk[j][k];
      }//end k
    }//end j
    m_stencil = Ajk;
  } else {
    int Bijk[8][8][8] = {
      //Type-0:No Hanging
      {{ 64, 32, 32, 16, 32, 16, 16,  8},
        { 32, 64, 16, 32, 16, 32,  8, 16},
        { 32, 16, 64, 32, 16,  8, 32, 16},
        { 16, 32, 32, 64,  8, 16, 16, 32},
        { 32, 16, 16,  8, 64, 32, 32, 16},
        { 16, 32,  8, 16, 32, 64, 16, 32},
        { 16,  8, 32, 16, 32, 16, 64, 32},
        {  8, 16, 16, 32, 16, 32, 32, 64}},

      //Type-1: Y Hanging
      {{ 112,  40,  32,  32,  40,  20,  32,  16},
        {  40,  64,   8,  32,  16,  32,   8,  16},
        {  32,   8,  16,  16,   8,   4,  16,   8},
        {  32,  32,  16,  64,   8,  16,  16,  32},
        {  40,  16,   8,   8,  64,  32,  32,  16},
        {  20,  32,   4,  16,  32,  64,  16,  32},
        {  32,   8,  16,  16,  32,  16,  64,  32},
        {  16,  16,   8,  32,  16,  32,  32,  64}},

      //Type-2: X and Y Hanging
      {{ 168,  36,  36,  48,  48,  36,  36,  24},
        {  36,  16,   4,  16,   8,  16,   4,   8},
        {  36,   4,  16,  16,   8,   4,  16,   8},
        {  48,  16,  16,  64,   8,  16,  16,  32},
        {  48,   8,   8,   8,  64,  32,  32,  16},
        {  36,  16,   4,  16,  32,  64,  16,  32},
        {  36,   4,  16,  16,  32,  16,  64,  32},
        {  24,   8,   8,  32,  16,  32,  32,  64}},

      //Type-3: X and Y and Z Hanging
      {{ 232,  40,  40,  52,  40,  52,  52,  32},
        {  40,  16,   4,  16,   4,  16,   4,   8},
        {  40,   4,  16,  16,   4,   4,  16,   8},
        {  52,  16,  16,  64,   4,  16,  16,  32},
        {  40,   4,   4,   4,  16,  16,  16,   8},
        {  52,  16,   4,  16,  16,  64,  16,  32},
        {  52,   4,  16,  16,  16,  16,  64,  32},
        {  32,   8,   8,  32,   8,  32,  32,  64}},

      //Type-4:XY and X and Y Hanging
      {{ 196,  56,  56,  16,  50,  40,  40,  32},
        {  56,  28,  16,   8,  10,  20,   8,  16},
        {  56,  16,  28,   8,  10,   8,  20,  16},
        {  16,   8,   8,   4,   2,   4,   4,   8},
        {  50,  10,  10,   2,  64,  32,  32,  16},
        {  40,  20,   8,   4,  32,  64,  16,  32},
        {  40,   8,  20,   4,  32,  16,  64,  32},
        {  32,  16,  16,   8,  16,  32,  32,  64}},

      //Type-5:XY and X and Y and Z Hanging
      {{ 262,  61,  61,  17,  41,  56,  56,  40},
        {  61,  28,  16,   8,   5,  20,   8,  16},
        {  61,  16,  28,   8,   5,   8,  20,  16},
        {  17,   8,   8,   4,   1,   4,   4,   8},
        {  41,   5,   5,   1,  16,  16,  16,   8},
        {  56,  20,   8,   4,  16,  64,  16,  32},
        {  56,   8,  20,   4,  16,  16,  64,  32},
        {  40,  16,  16,   8,   8,  32,  32,  64}},

      //Type-6:XY and YZ and X and Y and Z Hanging
      {{ 294,  63,  84,  18,  63,  60,  18,  48},
        {  63,  28,  18,   8,   7,  20,   2,  16},
        {  84,  18,  42,   9,  18,  12,   9,  24},
        {  18,   8,   9,   4,   2,   4,   1,   8},
        {  63,   7,  18,   2,  28,  20,   8,  16},
        {  60,  20,  12,   4,  20,  64,   4,  32},
        {  18,   2,   9,   1,   8,   4,   4,   8},
        {  48,  16,  24,   8,  16,  32,   8,  64}},

      //Type-7: All 6 Hanging
      {{ 328,  87,  87,  19,  87,  19,  19,  56},
        {  87,  42,  21,   9,  21,   9,   3,  24},
        {  87,  21,  42,   9,  21,   3,   9,  24},
        {  19,   9,   9,   4,   3,   1,   1,   8},
        {  87,  21,  21,   3,  42,   9,   9,  24},
        {  19,   9,   3,   1,   9,   4,   1,   8},
        {  19,   3,   9,   1,   9,   1,   4,   8},
        {  56,  24,  24,   8,  24,   8,   8,  64}}
    };

    int ***Aijk = new int2Ptr[8];
    for (int i=0;i<8;i++) {
      Aijk[i] = new int1Ptr[8];
      for (int j=0;j<8;j++) {
        Aijk[i][j] = new int[8];
        for (int k=0;k<8;k++) {
          Aijk[i][j][k] = Bijk[i][j][k];
        }//end k
      }//end j
    }//end i
    m_stencil = Aijk;
  }
  return true;
}

bool forceVector::preAddVec() {
  if (m_daType == PETSC) {
    // compute Hx
    int ierr;
    PetscInt mx,my,mz;
    ierr = DMDAGetInfo(m_DA,0, &mx, &my, &mz, 0,0, 0,0,0,0,0,0,0); CHKERRQ(ierr); 

    m_dHx = m_dLx/(mx -1);
    m_dHx = m_dHx*m_dHx*m_dHx;
    m_dHx /= 1728.0;
    CHKERRQ(ierr);

    PetscScalar*** arr;
    ierr = DMDAVecGetArray(m_DA, m_prevTSVec, &arr);
    m_prevTSArray = arr;
    CHKERRQ(ierr);
  } else {
    maxD = m_octDA->getMaxDepth();
    
    // Get the  x,y,z factors 
    xFac = 1.0/((double)(1<<(maxD-1)));
    if (m_octDA->getDimension() > 1) {
      yFac = 1.0/((double)(1<<(maxD-1)));
      if (m_octDA->getDimension() > 2) {
        zFac = 1.0/((double)(1<<(maxD-1)));
      }
    }

    PetscScalar* arr;
    m_octDA->vecGetBuffer(m_prevTSVec, arr, false, true, false, m_uiDof);
    m_prevTSArray = arr;
  }

  return true;
}

/*
bool massMatrix::ElementalMatVec(unsigned int i, PetscScalar *in, PetscScalar *out, double scale) {
  unsigned int lev = m_octDA->getLevel(i);
  double hx = xFac*(1<<(maxD - lev));
  double hy = yFac*(1<<(maxD - lev));
  double hz = zFac*(1<<(maxD - lev));

  double fac = scale*hx*hy*hz/1728.0;
  
  stdElemType elemType;
  unsigned int idx[8];
  
  int ***Aijk = (int ***)m_stencil;
  
  alignElementAndVertices(m_octDA, elemType, idx);       
  

  for (int k = 0;k < 8;k++) {
    for (int j=0;j<8;j++) {
      out[m_uiDof*idx[k]] += fac*(Aijk[elemType][k][j])*in[m_uiDof*idx[j]];
	}//end for j
  }//end for k
  
  return true;
}*/

bool forceVector::ElementalAddVec(unsigned int i, PetscScalar *out, double scale) {
  unsigned int lev = m_octDA->getLevel(i);
  double hx = xFac*(1<<(maxD - lev));
  double hy = yFac*(1<<(maxD - lev));
  double hz = zFac*(1<<(maxD - lev));

  double fac = scale*hx*hy*hz/1728.0;
  
  stdElemType elemType;
  unsigned int idx[8];
  
  int ***Aijk = (int ***)m_stencil;
  
  alignElementAndVertices(m_octDA, elemType, idx);       
  
  PetscScalar* in = (PetscScalar*) m_prevTSArray;
  for (int k = 0;k < 8;k++) {
    for (int j=0;j<8;j++) {
      out[m_uiDof*idx[k]] += fac*(Aijk[elemType][k][j])*in[m_uiDof*idx[j]];
    }//end for j
  }//end for k
}

bool forceVector::ElementalAddVec(int i, int j, int k, PetscScalar ***out, double scale) {
  int dof= m_uiDof;
  int idx[8][3]={
    {k, j, dof*i},
    {k,j,dof*(i+1)},
    {k,j+1,dof*i},
    {k,j+1,dof*(i+1)},
    {k+1, j, dof*i},
    {k+1,j,dof*(i+1)},
    {k+1,j+1,dof*i},
    {k+1,j+1,dof*(i+1)}               
  };             

  double stencilScale =  m_dHx*scale;
  int **Ajk = (int **)m_stencil;
  PetscScalar*** in = (PetscScalar***) m_prevTSArray;
  
  for (int q = 0; q < 8; q++) {
    for (int r = 0; r < 8; r++) {
      double contrib = stencilScale*Ajk[q][r]*in[idx[r][0]][idx[r][1]][idx[r][2]];
      out[idx[q][0]][idx[q][1]][idx[q][2]] += contrib;
    }
  }

  return true;
  // std::cout << "Stencil scale is " << stencilScale << std::endl;
}

bool forceVector::postAddVec() {
  if (m_daType == PETSC) {
    PetscScalar*** arr = (PetscScalar***) m_prevTSArray;
    int ierr = DMDAVecRestoreArray(m_DA, m_prevTSVec, &arr);
    CHKERRQ(ierr);
  } else {
    PetscScalar* arr = (PetscScalar*) m_prevTSArray;
    int ierr = m_octDA->vecRestoreBuffer(m_prevTSVec, arr, false,true,false,m_uiDof);
    CHKERRQ(ierr);
  }
  m_prevTSArray = NULL;
  return true;
}

