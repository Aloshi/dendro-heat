#include "interp.h"

#define NODE_0 1u
#define NODE_1 2u
#define NODE_2 4u
#define NODE_3 8u
#define NODE_4 16u
#define NODE_5 32u
#define NODE_6 64u
#define NODE_7 128u


void interp_global_to_local(PetscScalar* glo, PetscScalar* /*__restrict*/ loc, ot::DA* da) {
	unsigned int idx[8];
	unsigned char hangingMask = da->getHangingNodeIndex(da->curr());
	unsigned int chNum = da->getChildNumber();
	da->getNodeIndices(idx);

  unsigned int m_uiDof = 1;
  
  //std::cout << "chNum: " << chNum << std::endl;

  switch (chNum) {

    case 0:
      // 0,7 are not hanging
		  for (size_t i = 0; i < m_uiDof; i++) {
				loc[i] = glo[m_uiDof*idx[0]+i];

        if ( hangingMask & NODE_1 ) 
          loc[m_uiDof + i] = 0.5 * ( glo[m_uiDof*idx[0]+i] + glo[m_uiDof*idx[1]+i] );
        else
          loc[m_uiDof + i] = glo[m_uiDof*idx[1]+i];
				
        if ( hangingMask & NODE_2 ) 
          loc[2*m_uiDof + i] = 0.5 * ( glo[m_uiDof*idx[0]+i] + glo[m_uiDof*idx[2]+i] );
        else
          loc[2*m_uiDof + i] = glo[m_uiDof*idx[2]+i];
        
        if ( hangingMask & NODE_3 ) 
          loc[3*m_uiDof + i] = 0.25 * ( glo[m_uiDof*idx[0]+i] + glo[m_uiDof*idx[1]+i] + glo[m_uiDof*idx[2]+i] + glo[m_uiDof*idx[3]+i]);
        else
          loc[3*m_uiDof + i] = glo[m_uiDof*idx[3]+i];
        
        if ( hangingMask & NODE_4 ) 
          loc[4*m_uiDof + i] = 0.5 * ( glo[m_uiDof*idx[0]+i] + glo[m_uiDof*idx[4]+i] );
        else
          loc[4*m_uiDof + i] = glo[m_uiDof*idx[4]+i];
        
        if ( hangingMask & NODE_5 ) 
          loc[5*m_uiDof + i] = 0.25 * ( glo[m_uiDof*idx[0]+i] + glo[m_uiDof*idx[1]+i] + glo[m_uiDof*idx[4]+i] + glo[m_uiDof*idx[5]+i]);
        else
          loc[5*m_uiDof + i] = glo[m_uiDof*idx[5]+i];
        
        if ( hangingMask & NODE_6 ) 
          loc[6*m_uiDof + i] = 0.25 * ( glo[m_uiDof*idx[0]+i] + glo[m_uiDof*idx[2]+i] + glo[m_uiDof*idx[4]+i] + glo[m_uiDof*idx[6]+i]);
        else
          loc[6*m_uiDof + i] = glo[m_uiDof*idx[6]+i];

        loc[7*m_uiDof + i] = glo[m_uiDof*idx[7]+i];
      }
      break;
    case 1:
      // 1,6 are not hanging
		  for (size_t i = 0; i < m_uiDof; i++) {
        
        if ( hangingMask & NODE_0 ) 
          loc[i] = 0.5 * ( glo[m_uiDof*idx[0]+i] + glo[m_uiDof*idx[1]+i] );
        else
          loc[i] = glo[m_uiDof*idx[0]+i] ;
          
        loc[m_uiDof + i] = glo[m_uiDof*idx[1]+i];
				
        if ( hangingMask & NODE_2 ) 
          loc[2*m_uiDof + i] = 0.25 * ( glo[m_uiDof*idx[0]+i] + glo[m_uiDof*idx[1]+i] + glo[m_uiDof*idx[2]+i] + glo[m_uiDof*idx[3]+i]);
        else
          loc[2*m_uiDof + i] = glo[m_uiDof*idx[2]+i];
        
        if ( hangingMask & NODE_3 ) 
          loc[3*m_uiDof + i] = 0.5 * ( glo[m_uiDof*idx[1]+i] + glo[m_uiDof*idx[3]+i] );
        else
          loc[3*m_uiDof + i] = glo[m_uiDof*idx[3]+i];
        
        if ( hangingMask & NODE_4 ) 
          loc[4*m_uiDof + i] = 0.25 * ( glo[m_uiDof*idx[0]+i] + glo[m_uiDof*idx[1]+i] + glo[m_uiDof*idx[4]+i] + glo[m_uiDof*idx[5]+i]);
        else
          loc[4*m_uiDof + i] = glo[m_uiDof*idx[4]+i];
        
        if ( hangingMask & NODE_5 ) 
          loc[5*m_uiDof + i] = 0.5 * ( glo[m_uiDof*idx[1]+i] + glo[m_uiDof*idx[5]+i] );
        else
          loc[5*m_uiDof + i] = glo[m_uiDof*idx[5]+i];
        
        loc[6*m_uiDof + i] = glo[m_uiDof*idx[6]+i];
        
        if ( hangingMask & NODE_7 ) 
          loc[7*m_uiDof + i] = 0.25 * ( glo[m_uiDof*idx[1]+i] + glo[m_uiDof*idx[3]+i] + glo[m_uiDof*idx[5]+i] + glo[m_uiDof*idx[7]+i]);
        else
          loc[7*m_uiDof + i] = glo[m_uiDof*idx[7]+i];
      }
      break;
    case 2:
      // 2,5 are not hanging
		  for (size_t i = 0; i < m_uiDof; i++) {
        
        if ( hangingMask & NODE_0 ) 
          loc[i] = 0.5 * ( glo[m_uiDof*idx[0]+i] + glo[m_uiDof*idx[2]+i] );
        else
          loc[i] = glo[m_uiDof*idx[0]+i] ;
        
        if ( hangingMask & NODE_1 ) 
          loc[1*m_uiDof + i] = 0.25 * ( glo[m_uiDof*idx[0]+i] + glo[m_uiDof*idx[1]+i] + glo[m_uiDof*idx[2]+i] + glo[m_uiDof*idx[3]+i]);
        else
          loc[1*m_uiDof + i] = glo[m_uiDof*idx[1]+i];
				
        loc[2*m_uiDof + i] = glo[m_uiDof*idx[2]+i];
        
        if ( hangingMask & NODE_3 ) 
          loc[3*m_uiDof + i] = 0.5 * ( glo[m_uiDof*idx[2]+i] + glo[m_uiDof*idx[3]+i] );
        else
          loc[3*m_uiDof + i] = glo[m_uiDof*idx[3]+i];
        
        if ( hangingMask & NODE_4 ) 
          loc[4*m_uiDof + i] = 0.25 * ( glo[m_uiDof*idx[0]+i] + glo[m_uiDof*idx[2]+i] + glo[m_uiDof*idx[4]+i] + glo[m_uiDof*idx[6]+i]);
        else
          loc[4*m_uiDof + i] = glo[m_uiDof*idx[4]+i];
  
        loc[5*m_uiDof + i] = glo[m_uiDof*idx[5]+i];
        
        if ( hangingMask & NODE_6 ) 
          loc[6*m_uiDof + i] = 0.5 * ( glo[m_uiDof*idx[2]+i] + glo[m_uiDof*idx[6]+i] );
        else
          loc[6*m_uiDof + i] = glo[m_uiDof*idx[6]+i];
        
        if ( hangingMask & NODE_7 ) 
          loc[7*m_uiDof + i] = 0.25 * ( glo[m_uiDof*idx[2]+i] + glo[m_uiDof*idx[3]+i] + glo[m_uiDof*idx[6]+i] + glo[m_uiDof*idx[7]+i]);
        else
          loc[7*m_uiDof + i] = glo[m_uiDof*idx[7]+i];
      }
      break;
    case 3:
      // 3,4 are not hanging
		  for (size_t i = 0; i < m_uiDof; i++) {
        if ( hangingMask & NODE_0 ) 
          loc[i] = 0.25 * ( glo[m_uiDof*idx[0]+i] + glo[m_uiDof*idx[1]+i] + glo[m_uiDof*idx[2]+i] + glo[m_uiDof*idx[3]+i]);
        else
          loc[i] = glo[m_uiDof*idx[0]+i];
        
        if ( hangingMask & NODE_1 ) 
          loc[m_uiDof + i] = 0.5 * ( glo[m_uiDof*idx[1]+i] + glo[m_uiDof*idx[3]+i] );
        else
          loc[m_uiDof + i] = glo[m_uiDof*idx[1]+i];
				
        if ( hangingMask & NODE_2 ) 
          loc[2*m_uiDof + i] = 0.5 * ( glo[m_uiDof*idx[2]+i] + glo[m_uiDof*idx[3]+i] );
        else
          loc[2*m_uiDof + i] = glo[m_uiDof*idx[2]+i];
          
        loc[3*m_uiDof + i] = glo[m_uiDof*idx[3]+i];
        
        loc[4*m_uiDof + i] = glo[m_uiDof*idx[4]+i];
        
        if ( hangingMask & NODE_5 ) 
          loc[5*m_uiDof + i] = 0.25 * ( glo[m_uiDof*idx[1]+i] + glo[m_uiDof*idx[3]+i] + glo[m_uiDof*idx[5]+i] + glo[m_uiDof*idx[7]+i]);
        else
          loc[5*m_uiDof + i] = glo[m_uiDof*idx[5]+i];
        
        if ( hangingMask & NODE_6 ) 
          loc[6*m_uiDof + i] = 0.25 * ( glo[m_uiDof*idx[2]+i] + glo[m_uiDof*idx[3]+i] + glo[m_uiDof*idx[6]+i] + glo[m_uiDof*idx[7]+i]);
        else
          loc[6*m_uiDof + i] = glo[m_uiDof*idx[6]+i];
        
        if ( hangingMask & NODE_7 ) 
          loc[7*m_uiDof + i] = 0.5 * ( glo[m_uiDof*idx[3]+i] + glo[m_uiDof*idx[7]+i] );
        else
          loc[7*m_uiDof + i] = glo[m_uiDof*idx[7]+i];
      }
      break;
    case 4:
		  // 4,3 are not hanging
      for (size_t i = 0; i < m_uiDof; i++) {
        if ( hangingMask & NODE_0 ) 
          loc[i] = 0.5 * ( glo[m_uiDof*idx[0]+i] + glo[m_uiDof*idx[4]+i] );
        else
          loc[i] = glo[m_uiDof*idx[0]+i];
				
        if ( hangingMask & NODE_1 ) 
          loc[1*m_uiDof + i] = 0.25 * ( glo[m_uiDof*idx[0]+i] + glo[m_uiDof*idx[1]+i] + glo[m_uiDof*idx[4]+i] + glo[m_uiDof*idx[5]+i]);
        else
          loc[1*m_uiDof + i] = glo[m_uiDof*idx[1]+i];
        
        if ( hangingMask & NODE_2 ) 
          loc[2*m_uiDof + i] = 0.25 * ( glo[m_uiDof*idx[0]+i] + glo[m_uiDof*idx[2]+i] + glo[m_uiDof*idx[4]+i] + glo[m_uiDof*idx[6]+i]);
        else
          loc[2*m_uiDof + i] = glo[m_uiDof*idx[2]+i];
          
        loc[3*m_uiDof + i] = glo[m_uiDof*idx[3]+i];
          
        loc[4*m_uiDof + i] = glo[m_uiDof*idx[4]+i];
        
        if ( hangingMask & NODE_5 ) 
          loc[5*m_uiDof + i] = 0.5 * ( glo[m_uiDof*idx[4]+i] + glo[m_uiDof*idx[5]+i] );
        else
          loc[5*m_uiDof + i] = glo[m_uiDof*idx[5]+i];
        
        if ( hangingMask & NODE_6 ) 
          loc[6*m_uiDof + i] = 0.5 * ( glo[m_uiDof*idx[4]+i] + glo[m_uiDof*idx[6]+i] );
        else
          loc[6*m_uiDof + i] = glo[m_uiDof*idx[6]+i];
        
        if ( hangingMask & NODE_7 ) 
          loc[7*m_uiDof + i] = 0.25 * ( glo[m_uiDof*idx[4]+i] + glo[m_uiDof*idx[5]+i] + glo[m_uiDof*idx[6]+i] + glo[m_uiDof*idx[7]+i]);
        else
          loc[7*m_uiDof + i] = glo[m_uiDof*idx[7]+i];
      }
      break;
    case 5:
      // 5,2 are not hanging
      for (size_t i = 0; i < m_uiDof; i++) {
        if ( hangingMask & NODE_0 ) 
          loc[i] = 0.25 * ( glo[m_uiDof*idx[0]+i] + glo[m_uiDof*idx[1]+i] + glo[m_uiDof*idx[4]+i] + glo[m_uiDof*idx[5]+i]);
        else
          loc[i] = glo[m_uiDof*idx[0]+i];
        
        if ( hangingMask & NODE_1 ) 
          loc[m_uiDof + i] = 0.5 * ( glo[m_uiDof*idx[1]+i] + glo[m_uiDof*idx[5]+i] );
        else
          loc[m_uiDof + i] = glo[m_uiDof*idx[1]+i];
        
        loc[2*m_uiDof + i] = glo[m_uiDof*idx[2]+i];
        
        if ( hangingMask & NODE_3 ) 
          loc[3*m_uiDof + i] = 0.25 * ( glo[m_uiDof*idx[1]+i] + glo[m_uiDof*idx[3]+i] + glo[m_uiDof*idx[5]+i] + glo[m_uiDof*idx[7]+i]);
        else
          loc[3*m_uiDof + i] = glo[m_uiDof*idx[3]+i];
          
        if ( hangingMask & NODE_4 ) 
          loc[4*m_uiDof + i] = 0.5 * ( glo[m_uiDof*idx[4]+i] + glo[m_uiDof*idx[5]+i] );
        else
          loc[4*m_uiDof + i] = glo[m_uiDof*idx[4]+i];
          
        loc[5*m_uiDof + i] = glo[m_uiDof*idx[5]+i];
        
        if ( hangingMask & NODE_6 ) 
          loc[6*m_uiDof + i] = 0.25 * ( glo[m_uiDof*idx[4]+i] + glo[m_uiDof*idx[5]+i] + glo[m_uiDof*idx[6]+i] + glo[m_uiDof*idx[7]+i]);
        else
          loc[6*m_uiDof + i] = glo[m_uiDof*idx[6]+i];
        
        if ( hangingMask & NODE_7 ) 
          loc[7*m_uiDof + i] = 0.5 * ( glo[m_uiDof*idx[5]+i] + glo[m_uiDof*idx[7]+i] );
        else
          loc[7*m_uiDof + i] = glo[m_uiDof*idx[7]+i];
      }
      break;
    case 6:
      // 6,1 are not hanging
      for (size_t i = 0; i < m_uiDof; i++) {
        if ( hangingMask & NODE_0 ) 
          loc[i] = 0.25 * ( glo[m_uiDof*idx[0]+i] + glo[m_uiDof*idx[1]+i] + glo[m_uiDof*idx[4]+i] + glo[m_uiDof*idx[5]+i]);
        else
          loc[i] = glo[m_uiDof*idx[0]+i];
      
        loc[m_uiDof + i] = glo[m_uiDof*idx[1]+i];

        if ( hangingMask & NODE_2 ) 
          loc[2*m_uiDof + i] = 0.5 * ( glo[m_uiDof*idx[2]+i] + glo[m_uiDof*idx[6]+i] );
        else
          loc[2*m_uiDof + i] = glo[m_uiDof*idx[2]+i];

        if ( hangingMask & NODE_3 ) 
          loc[3*m_uiDof + i] = 0.25 * ( glo[m_uiDof*idx[2]+i] + glo[m_uiDof*idx[3]+i] + glo[m_uiDof*idx[6]+i] + glo[m_uiDof*idx[7]+i]);
        else
          loc[3*m_uiDof + i] = glo[m_uiDof*idx[3]+i];

        if ( hangingMask & NODE_4 ) 
          loc[4*m_uiDof + i] = 0.5 * ( glo[m_uiDof*idx[4]+i] + glo[m_uiDof*idx[6]+i] );
        else
          loc[4*m_uiDof + i] = glo[m_uiDof*idx[4]+i];

        if ( hangingMask & NODE_5 ) 
          loc[5*m_uiDof + i] = 0.25 * ( glo[m_uiDof*idx[4]+i] + glo[m_uiDof*idx[5]+i] + glo[m_uiDof*idx[6]+i] + glo[m_uiDof*idx[7]+i]);
        else
          loc[5*m_uiDof + i] = glo[m_uiDof*idx[5]+i];
          
        loc[6*m_uiDof + i] = glo[m_uiDof*idx[6]+i];
        
        if ( hangingMask & NODE_7 ) 
          loc[7*m_uiDof + i] = 0.5 * ( glo[m_uiDof*idx[6]+i] + glo[m_uiDof*idx[7]+i] );
        else
          loc[7*m_uiDof + i] = glo[m_uiDof*idx[7]+i];
      
      }
      break;
    case 7:
      // 7,0 are not hanging
      for (size_t i = 0; i < m_uiDof; i++) {
        loc[i] = glo[m_uiDof*idx[0]+i];
        
        if ( hangingMask & NODE_1 ) 
          loc[1*m_uiDof + i] = 0.25 * ( glo[m_uiDof*idx[1]+i] + glo[m_uiDof*idx[3]+i] + glo[m_uiDof*idx[5]+i] + glo[m_uiDof*idx[7]+i]);
        else
          loc[1*m_uiDof + i] = glo[m_uiDof*idx[1]+i];
        
        if ( hangingMask & NODE_2 ) 
          loc[2*m_uiDof + i] = 0.25 * ( glo[m_uiDof*idx[2]+i] + glo[m_uiDof*idx[3]+i] + glo[m_uiDof*idx[6]+i] + glo[m_uiDof*idx[7]+i]);
        else
          loc[2*m_uiDof + i] = glo[m_uiDof*idx[2]+i];
        
        if ( hangingMask & NODE_3 ) 
          loc[3*m_uiDof + i] = 0.5 * ( glo[m_uiDof*idx[3]+i] + glo[m_uiDof*idx[7]+i] );
        else
          loc[3*m_uiDof + i] = glo[m_uiDof*idx[3]+i];
        
        if ( hangingMask & NODE_4 ) 
          loc[4*m_uiDof + i] = 0.25 * ( glo[m_uiDof*idx[4]+i] + glo[m_uiDof*idx[5]+i] + glo[m_uiDof*idx[6]+i] + glo[m_uiDof*idx[7]+i]);
        else
          loc[4*m_uiDof + i] = glo[m_uiDof*idx[4]+i];

        if ( hangingMask & NODE_5 ) 
          loc[5*m_uiDof + i] = 0.5 * ( glo[m_uiDof*idx[5]+i] + glo[m_uiDof*idx[7]+i] );
        else
          loc[5*m_uiDof + i] = glo[m_uiDof*idx[5]+i];

        if ( hangingMask & NODE_6 ) 
          loc[6*m_uiDof + i] = 0.5 * ( glo[m_uiDof*idx[6]+i] + glo[m_uiDof*idx[7]+i] );
        else
          loc[6*m_uiDof + i] = glo[m_uiDof*idx[6]+i];
        
        loc[7*m_uiDof + i] = glo[m_uiDof*idx[7]+i];
      }
      break;
    default:
			std::cout<<"in glo_to_loc: incorrect child num = " << chNum << std::endl;
			assert(false);
      break;
  } // switch


}

void interp_local_to_global(PetscScalar* /*__restrict*/ loc, PetscScalar* glo, ot::DA* da) {
  unsigned int idx[8];
	unsigned char hangingMask = da->getHangingNodeIndex(da->curr());
	unsigned int chNum = da->getChildNumber();
	da->getNodeIndices(idx);

  unsigned int m_uiDof = 1;

  switch (chNum) {
    case 0:
      // 0,7 are not hanging
      for (size_t i = 0; i < m_uiDof; i++) {
        glo[m_uiDof*idx[0]+i] += loc[i];
        if ( hangingMask & NODE_1 ) {
          glo[m_uiDof*idx[1]+i] += 0.5*loc[m_uiDof + i];
          glo[m_uiDof*idx[0]+i] += 0.5*loc[m_uiDof + i];
        } else {
          glo[m_uiDof*idx[1]+i] += loc[m_uiDof + i];
        }
        if ( hangingMask & NODE_2 ) {
          glo[m_uiDof*idx[2]+i] += 0.5*loc[2*m_uiDof + i];
          glo[m_uiDof*idx[0]+i] += 0.5*loc[2*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[2]+i] += loc[2*m_uiDof + i];
        }
        if ( hangingMask & NODE_3 ) {
          glo[m_uiDof*idx[0]+i] += 0.25*loc[3*m_uiDof + i];
          glo[m_uiDof*idx[1]+i] += 0.25*loc[3*m_uiDof + i];
          glo[m_uiDof*idx[2]+i] += 0.25*loc[3*m_uiDof + i];
          glo[m_uiDof*idx[3]+i] += 0.25*loc[3*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[3]+i] += loc[3*m_uiDof + i];
        }
        if ( hangingMask & NODE_4 ) {
          glo[m_uiDof*idx[4]+i] += 0.5*loc[4*m_uiDof + i];
          glo[m_uiDof*idx[0]+i] += 0.5*loc[4*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[4]+i] += loc[4*m_uiDof + i];
        }
        if ( hangingMask & NODE_5 ) {
          glo[m_uiDof*idx[0]+i] += 0.25*loc[5*m_uiDof + i];
          glo[m_uiDof*idx[1]+i] += 0.25*loc[5*m_uiDof + i];
          glo[m_uiDof*idx[4]+i] += 0.25*loc[5*m_uiDof + i];
          glo[m_uiDof*idx[5]+i] += 0.25*loc[5*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[5]+i] += loc[5*m_uiDof + i];
        }
        if ( hangingMask & NODE_6 ) {
          glo[m_uiDof*idx[0]+i] += 0.25*loc[6*m_uiDof + i];
          glo[m_uiDof*idx[2]+i] += 0.25*loc[6*m_uiDof + i];
          glo[m_uiDof*idx[4]+i] += 0.25*loc[6*m_uiDof + i];
          glo[m_uiDof*idx[6]+i] += 0.25*loc[6*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[6]+i] += loc[6*m_uiDof + i];
        }
        glo[m_uiDof*idx[7]+i] += loc[7*m_uiDof + i];
      }
      break;
    case 1:
		  for (size_t i = 0; i < m_uiDof; i++) {
        if ( hangingMask & NODE_0 ) {
          glo[m_uiDof*idx[0]+i] += 0.5*loc[i];
          glo[m_uiDof*idx[1]+i] += 0.5*loc[i];
        } else {
          glo[m_uiDof*idx[0]+i] += loc[i];
        }

        glo[m_uiDof*idx[1]+i] += loc[m_uiDof + i];

        if ( hangingMask & NODE_2 ) {
          glo[m_uiDof*idx[0]+i] += 0.25*loc[2*m_uiDof + i];
          glo[m_uiDof*idx[1]+i] += 0.25*loc[2*m_uiDof + i];
          glo[m_uiDof*idx[2]+i] += 0.25*loc[2*m_uiDof + i];
          glo[m_uiDof*idx[3]+i] += 0.25*loc[2*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[2]+i] += loc[2*m_uiDof + i];
        }

        if ( hangingMask & NODE_3 ) {
          glo[m_uiDof*idx[1]+i] += 0.5*loc[3*m_uiDof + i];
          glo[m_uiDof*idx[3]+i] += 0.5*loc[3*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[3]+i] += loc[3*m_uiDof + i];
        }

        if ( hangingMask & NODE_4 ) {
          glo[m_uiDof*idx[0]+i] += 0.25*loc[4*m_uiDof + i];
          glo[m_uiDof*idx[1]+i] += 0.25*loc[4*m_uiDof + i];
          glo[m_uiDof*idx[4]+i] += 0.25*loc[4*m_uiDof + i];
          glo[m_uiDof*idx[5]+i] += 0.25*loc[4*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[4]+i] += loc[4*m_uiDof + i];
        }

        if ( hangingMask & NODE_5 ) {
          glo[m_uiDof*idx[1]+i] += 0.5*loc[5*m_uiDof + i];
          glo[m_uiDof*idx[5]+i] += 0.5*loc[5*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[5]+i] += loc[5*m_uiDof + i];
        }
          
        glo[m_uiDof*idx[6]+i] += loc[6*m_uiDof + i];
        
        if ( hangingMask & NODE_7 ) {
          glo[m_uiDof*idx[1]+i] += 0.25*loc[7*m_uiDof + i];
          glo[m_uiDof*idx[3]+i] += 0.25*loc[7*m_uiDof + i];
          glo[m_uiDof*idx[5]+i] += 0.25*loc[7*m_uiDof + i];
          glo[m_uiDof*idx[7]+i] += 0.25*loc[7*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[7]+i] += loc[7*m_uiDof + i];
        }
      }
      break;
    case 2:
      // 2,5 are not hanging
		  for (size_t i = 0; i < m_uiDof; i++) {
        if ( hangingMask & NODE_0 ) {
          glo[m_uiDof*idx[0]+i] += 0.5*loc[i];
          glo[m_uiDof*idx[2]+i] += 0.5*loc[i];
        } else {
          glo[m_uiDof*idx[0]+i] += loc[i];
        }

        if ( hangingMask & NODE_1 ) {
          glo[m_uiDof*idx[0]+i] += 0.25*loc[m_uiDof + i];
          glo[m_uiDof*idx[1]+i] += 0.25*loc[m_uiDof + i];
          glo[m_uiDof*idx[2]+i] += 0.25*loc[m_uiDof + i];
          glo[m_uiDof*idx[3]+i] += 0.25*loc[m_uiDof + i];
        } else {
          glo[m_uiDof*idx[1]+i] += loc[m_uiDof + i];
        }

        glo[m_uiDof*idx[2]+i] += loc[2*m_uiDof + i];

        if ( hangingMask & NODE_3 ) {
          glo[m_uiDof*idx[2]+i] += 0.5*loc[3*m_uiDof + i];
          glo[m_uiDof*idx[3]+i] += 0.5*loc[3*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[3]+i] += loc[3*m_uiDof + i];
        }

        if ( hangingMask & NODE_4 ) {
          glo[m_uiDof*idx[0]+i] += 0.25*loc[4*m_uiDof + i];
          glo[m_uiDof*idx[2]+i] += 0.25*loc[4*m_uiDof + i];
          glo[m_uiDof*idx[4]+i] += 0.25*loc[4*m_uiDof + i];
          glo[m_uiDof*idx[6]+i] += 0.25*loc[4*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[4]+i] += loc[4*m_uiDof + i];
        }

        glo[m_uiDof*idx[5]+i] += loc[5*m_uiDof + i];
          
        if ( hangingMask & NODE_6 ) {
          glo[m_uiDof*idx[2]+i] += 0.5*loc[6*m_uiDof + i];
          glo[m_uiDof*idx[6]+i] += 0.5*loc[6*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[6]+i] += loc[6*m_uiDof + i];
        }
        
        if ( hangingMask & NODE_7 ) {
          glo[m_uiDof*idx[2]+i] += 0.25*loc[7*m_uiDof + i];
          glo[m_uiDof*idx[3]+i] += 0.25*loc[7*m_uiDof + i];
          glo[m_uiDof*idx[6]+i] += 0.25*loc[7*m_uiDof + i];
          glo[m_uiDof*idx[7]+i] += 0.25*loc[7*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[7]+i] += loc[7*m_uiDof + i];
        }
      }
      break;
    case 3:
      // 3,4 are not hanging
		  for (size_t i = 0; i < m_uiDof; i++) {
        if ( hangingMask & NODE_0 ) {
          glo[m_uiDof*idx[0]+i] += 0.25*loc[i];
          glo[m_uiDof*idx[1]+i] += 0.25*loc[i];
          glo[m_uiDof*idx[2]+i] += 0.25*loc[i];
          glo[m_uiDof*idx[3]+i] += 0.25*loc[i];
        } else {
          glo[m_uiDof*idx[0]+i] += loc[i];
        }

        if ( hangingMask & NODE_1 ) {
          glo[m_uiDof*idx[1]+i] += 0.5*loc[m_uiDof + i];
          glo[m_uiDof*idx[3]+i] += 0.5*loc[m_uiDof + i];
        } else {
          glo[m_uiDof*idx[1]+i] += loc[m_uiDof + i];
        }

        if ( hangingMask & NODE_2 ) {
          glo[m_uiDof*idx[2]+i] += 0.5*loc[2*m_uiDof + i];
          glo[m_uiDof*idx[3]+i] += 0.5*loc[2*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[2]+i] += loc[2*m_uiDof + i];
        }

        glo[m_uiDof*idx[3]+i] += loc[3*m_uiDof + i];
        glo[m_uiDof*idx[4]+i] += loc[4*m_uiDof + i];

        if ( hangingMask & NODE_5 ) {
          glo[m_uiDof*idx[1]+i] += 0.25*loc[5*m_uiDof + i];
          glo[m_uiDof*idx[3]+i] += 0.25*loc[5*m_uiDof + i];
          glo[m_uiDof*idx[5]+i] += 0.25*loc[5*m_uiDof + i];
          glo[m_uiDof*idx[7]+i] += 0.25*loc[5*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[5]+i] += loc[5*m_uiDof + i];
        }
        if ( hangingMask & NODE_6 ) {
          glo[m_uiDof*idx[2]+i] += 0.25*loc[6*m_uiDof + i];
          glo[m_uiDof*idx[3]+i] += 0.25*loc[6*m_uiDof + i];
          glo[m_uiDof*idx[6]+i] += 0.25*loc[6*m_uiDof + i];
          glo[m_uiDof*idx[7]+i] += 0.25*loc[6*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[6]+i] += loc[6*m_uiDof + i];
        }
        if ( hangingMask & NODE_7 ) {
          glo[m_uiDof*idx[3]+i] += 0.5*loc[7*m_uiDof + i];
          glo[m_uiDof*idx[7]+i] += 0.5*loc[7*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[7]+i] += loc[7*m_uiDof + i];
        }
      }
      break;
    case 4:
		  // 4,3 are not hanging
      for (size_t i = 0; i < m_uiDof; i++) {
        if ( hangingMask & NODE_0 ) {
          glo[m_uiDof*idx[0]+i] += 0.5*loc[i];
          glo[m_uiDof*idx[4]+i] += 0.5*loc[i];
        } else {
          glo[m_uiDof*idx[0]+i] += loc[i];
        }
        if ( hangingMask & NODE_1 ) {
          glo[m_uiDof*idx[0]+i] += 0.25*loc[m_uiDof + i];
          glo[m_uiDof*idx[1]+i] += 0.25*loc[m_uiDof + i];
          glo[m_uiDof*idx[4]+i] += 0.25*loc[m_uiDof + i];
          glo[m_uiDof*idx[5]+i] += 0.25*loc[m_uiDof + i];
        } else {
          glo[m_uiDof*idx[1]+i] += loc[m_uiDof + i];
        }

        if ( hangingMask & NODE_2 ) {
          glo[m_uiDof*idx[0]+i] += 0.25*loc[2*m_uiDof + i];
          glo[m_uiDof*idx[2]+i] += 0.25*loc[2*m_uiDof + i];
          glo[m_uiDof*idx[4]+i] += 0.25*loc[2*m_uiDof + i];
          glo[m_uiDof*idx[6]+i] += 0.25*loc[2*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[2]+i] += loc[2*m_uiDof + i];
        }

        glo[m_uiDof*idx[3]+i] += loc[3*m_uiDof + i];
        glo[m_uiDof*idx[4]+i] += loc[4*m_uiDof + i];
        
        if ( hangingMask & NODE_5 ) {
          glo[m_uiDof*idx[4]+i] += 0.5*loc[5*m_uiDof + i];
          glo[m_uiDof*idx[5]+i] += 0.5*loc[5*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[5]+i] += loc[5*m_uiDof + i];
        }
        if ( hangingMask & NODE_6 ) {
          glo[m_uiDof*idx[4]+i] += 0.5*loc[6*m_uiDof + i];
          glo[m_uiDof*idx[6]+i] += 0.5*loc[6*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[6]+i] += loc[6*m_uiDof + i];
        }
        if ( hangingMask & NODE_7 ) {
          glo[m_uiDof*idx[4]+i] += 0.25*loc[7*m_uiDof + i];
          glo[m_uiDof*idx[5]+i] += 0.25*loc[7*m_uiDof + i];
          glo[m_uiDof*idx[6]+i] += 0.25*loc[7*m_uiDof + i];
          glo[m_uiDof*idx[7]+i] += 0.25*loc[7*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[7]+i] += loc[7*m_uiDof + i];
        }
      }
      break;
    case 5:
      // 5,2 are not hanging
      for (size_t i = 0; i < m_uiDof; i++) {
        if ( hangingMask & NODE_0 ) {
          glo[m_uiDof*idx[0]+i] += 0.25*loc[i];
          glo[m_uiDof*idx[1]+i] += 0.25*loc[i];
          glo[m_uiDof*idx[4]+i] += 0.25*loc[i];
          glo[m_uiDof*idx[5]+i] += 0.25*loc[i];
        } else {
          glo[m_uiDof*idx[0]+i] += loc[i];
        }
        if ( hangingMask & NODE_1 ) {
          glo[m_uiDof*idx[1]+i] += 0.5*loc[m_uiDof + i];
          glo[m_uiDof*idx[5]+i] += 0.5*loc[m_uiDof + i];
        } else {
          glo[m_uiDof*idx[1]+i] += loc[m_uiDof + i];
        }
        glo[m_uiDof*idx[2]+i] += loc[2*m_uiDof + i];
        if ( hangingMask & NODE_3 ) {
          glo[m_uiDof*idx[1]+i] += 0.25*loc[3*m_uiDof + i];
          glo[m_uiDof*idx[3]+i] += 0.25*loc[3*m_uiDof + i];
          glo[m_uiDof*idx[5]+i] += 0.25*loc[3*m_uiDof + i];
          glo[m_uiDof*idx[7]+i] += 0.25*loc[3*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[3]+i] += loc[3*m_uiDof + i];
        }
        if ( hangingMask & NODE_4 ) {
          glo[m_uiDof*idx[4]+i] += 0.5*loc[4*m_uiDof + i];
          glo[m_uiDof*idx[5]+i] += 0.5*loc[4*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[4]+i] += loc[4*m_uiDof + i];
        }
        glo[m_uiDof*idx[5]+i] += loc[5*m_uiDof + i];
        if ( hangingMask & NODE_6 ) {
          glo[m_uiDof*idx[4]+i] += 0.25*loc[6*m_uiDof + i];
          glo[m_uiDof*idx[5]+i] += 0.25*loc[6*m_uiDof + i];
          glo[m_uiDof*idx[6]+i] += 0.25*loc[6*m_uiDof + i];
          glo[m_uiDof*idx[7]+i] += 0.25*loc[6*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[6]+i] += loc[6*m_uiDof + i];
        }
        if ( hangingMask & NODE_7 ) {
          glo[m_uiDof*idx[5]+i] += 0.5*loc[7*m_uiDof + i];
          glo[m_uiDof*idx[7]+i] += 0.5*loc[7*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[7]+i] += loc[7*m_uiDof + i];
        }
      }
      break;
    case 6:
      // 6,1 are not hanging
      for (size_t i = 0; i < m_uiDof; i++) {
        if ( hangingMask & NODE_0 ) {
          glo[m_uiDof*idx[0]+i] += 0.25*loc[i];
          glo[m_uiDof*idx[1]+i] += 0.25*loc[i];
          glo[m_uiDof*idx[4]+i] += 0.25*loc[i];
          glo[m_uiDof*idx[5]+i] += 0.25*loc[i];
        } else {
          glo[m_uiDof*idx[0]+i] += loc[i];
        }
        glo[m_uiDof*idx[1]+i] += loc[m_uiDof + i];
        if ( hangingMask & NODE_2 ) {
          glo[m_uiDof*idx[2]+i] += 0.5*loc[2*m_uiDof + i];
          glo[m_uiDof*idx[6]+i] += 0.5*loc[2*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[2]+i] += loc[2*m_uiDof + i];
        }
        if ( hangingMask & NODE_3 ) {
          glo[m_uiDof*idx[2]+i] += 0.25*loc[3*m_uiDof + i];
          glo[m_uiDof*idx[3]+i] += 0.25*loc[3*m_uiDof + i];
          glo[m_uiDof*idx[6]+i] += 0.25*loc[3*m_uiDof + i];
          glo[m_uiDof*idx[7]+i] += 0.25*loc[3*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[3]+i] += loc[3*m_uiDof + i];
        }
        if ( hangingMask & NODE_4 ) {
          glo[m_uiDof*idx[4]+i] += 0.5*loc[4*m_uiDof + i];
          glo[m_uiDof*idx[6]+i] += 0.5*loc[4*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[4]+i] += loc[4*m_uiDof + i];
        }
        if ( hangingMask & NODE_5 ) {
          glo[m_uiDof*idx[4]+i] += 0.25*loc[5*m_uiDof + i];
          glo[m_uiDof*idx[5]+i] += 0.25*loc[5*m_uiDof + i];
          glo[m_uiDof*idx[6]+i] += 0.25*loc[5*m_uiDof + i];
          glo[m_uiDof*idx[7]+i] += 0.25*loc[5*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[5]+i] += loc[5*m_uiDof + i];
        }
        glo[m_uiDof*idx[6]+i] += loc[6*m_uiDof + i];
        if ( hangingMask & NODE_7 ) {
          glo[m_uiDof*idx[6]+i] += 0.5*loc[7*m_uiDof + i];
          glo[m_uiDof*idx[7]+i] += 0.5*loc[7*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[7]+i] += loc[7*m_uiDof + i];
        }
      }
      break;
    case 7:
      // 7,0 are not hanging
      for (size_t i = 0; i < m_uiDof; i++) {
        glo[m_uiDof*idx[0]+i] += loc[i];
        if ( hangingMask & NODE_1 ) {
          glo[m_uiDof*idx[1]+i] += 0.25*loc[m_uiDof + i];
          glo[m_uiDof*idx[3]+i] += 0.25*loc[m_uiDof + i];
          glo[m_uiDof*idx[5]+i] += 0.25*loc[m_uiDof + i];
          glo[m_uiDof*idx[7]+i] += 0.25*loc[m_uiDof + i];
        } else {
          glo[m_uiDof*idx[1]+i] += loc[m_uiDof + i];
        }
        if ( hangingMask & NODE_2 ) {
          glo[m_uiDof*idx[2]+i] += 0.25*loc[2*m_uiDof + i];
          glo[m_uiDof*idx[3]+i] += 0.25*loc[2*m_uiDof + i];
          glo[m_uiDof*idx[6]+i] += 0.25*loc[2*m_uiDof + i];
          glo[m_uiDof*idx[7]+i] += 0.25*loc[2*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[2]+i] += loc[2*m_uiDof + i];
        }
        if ( hangingMask & NODE_3 ) {
          glo[m_uiDof*idx[3]+i] += 0.5*loc[3*m_uiDof + i];
          glo[m_uiDof*idx[7]+i] += 0.5*loc[3*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[3]+i] += loc[3*m_uiDof + i];
        }

        if ( hangingMask & NODE_4 ) {
          glo[m_uiDof*idx[4]+i] += 0.25*loc[4*m_uiDof + i];
          glo[m_uiDof*idx[5]+i] += 0.25*loc[4*m_uiDof + i];
          glo[m_uiDof*idx[6]+i] += 0.25*loc[4*m_uiDof + i];
          glo[m_uiDof*idx[7]+i] += 0.25*loc[4*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[4]+i] += loc[4*m_uiDof + i];
        }
        if ( hangingMask & NODE_5 ) {
          glo[m_uiDof*idx[5]+i] += 0.5*loc[5*m_uiDof + i];
          glo[m_uiDof*idx[7]+i] += 0.5*loc[5*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[5]+i] += loc[5*m_uiDof + i];
        }
        if ( hangingMask & NODE_6 ) {
          glo[m_uiDof*idx[6]+i] += 0.5*loc[6*m_uiDof + i];
          glo[m_uiDof*idx[7]+i] += 0.5*loc[6*m_uiDof + i];
        } else {
          glo[m_uiDof*idx[6]+i] += loc[6*m_uiDof + i];
        }
        glo[m_uiDof*idx[7]+i] += loc[7*m_uiDof + i];
      }
    break;
    default:
			std::cout<<"in loc_to_glo: incorrect child num = " << chNum << std::endl;
			assert(false);
      break;
  } // switch chNum
} // loc_to_glo



void interp_local_to_global_still_local(PetscScalar* /*__restrict*/ loc, PetscScalar* loc_new, ot::DA* da) {
  //unsigned int idx[8];
	unsigned char hangingMask = da->getHangingNodeIndex(da->curr());
	unsigned int chNum = da->getChildNumber();
  //unsigned int idx[8];
	//da->getNodeIndices(idx);

  unsigned int m_uiDof = 1;

  switch (chNum) {
    case 0:
      // 0,7 are not hanging
      for (size_t i = 0; i < m_uiDof; i++) {
        loc_new[m_uiDof*0+i] += loc[i];
        if ( hangingMask & NODE_1 ) {
          loc_new[m_uiDof*1+i] += 0.5*loc[m_uiDof + i];
          loc_new[m_uiDof*0+i] += 0.5*loc[m_uiDof + i];
        } else {
          loc_new[m_uiDof*1+i] += loc[m_uiDof + i];
        }
        if ( hangingMask & NODE_2 ) {
          loc_new[m_uiDof*2+i] += 0.5*loc[2*m_uiDof + i];
          loc_new[m_uiDof*0+i] += 0.5*loc[2*m_uiDof + i];
        } else {
          loc_new[m_uiDof*2+i] += loc[2*m_uiDof + i];
        }
        if ( hangingMask & NODE_3 ) {
          loc_new[m_uiDof*0+i] += 0.25*loc[3*m_uiDof + i];
          loc_new[m_uiDof*1+i] += 0.25*loc[3*m_uiDof + i];
          loc_new[m_uiDof*2+i] += 0.25*loc[3*m_uiDof + i];
          loc_new[m_uiDof*3+i] += 0.25*loc[3*m_uiDof + i];
        } else {
          loc_new[m_uiDof*3+i] += loc[3*m_uiDof + i];
        }
        if ( hangingMask & NODE_4 ) {
          loc_new[m_uiDof*4+i] += 0.5*loc[4*m_uiDof + i];
          loc_new[m_uiDof*0+i] += 0.5*loc[4*m_uiDof + i];
        } else {
          loc_new[m_uiDof*4+i] += loc[4*m_uiDof + i];
        }
        if ( hangingMask & NODE_5 ) {
          loc_new[m_uiDof*0+i] += 0.25*loc[5*m_uiDof + i];
          loc_new[m_uiDof*1+i] += 0.25*loc[5*m_uiDof + i];
          loc_new[m_uiDof*4+i] += 0.25*loc[5*m_uiDof + i];
          loc_new[m_uiDof*5+i] += 0.25*loc[5*m_uiDof + i];
        } else {
          loc_new[m_uiDof*5+i] += loc[5*m_uiDof + i];
        }
        if ( hangingMask & NODE_6 ) {
          loc_new[m_uiDof*0+i] += 0.25*loc[6*m_uiDof + i];
          loc_new[m_uiDof*2+i] += 0.25*loc[6*m_uiDof + i];
          loc_new[m_uiDof*4+i] += 0.25*loc[6*m_uiDof + i];
          loc_new[m_uiDof*6+i] += 0.25*loc[6*m_uiDof + i];
        } else {
          loc_new[m_uiDof*6+i] += loc[6*m_uiDof + i];
        }
        loc_new[m_uiDof*7+i] += loc[7*m_uiDof + i];
      }
      break;
    case 1:
		  for (size_t i = 0; i < m_uiDof; i++) {
        if ( hangingMask & NODE_0 ) {
          loc_new[m_uiDof*0+i] += 0.5*loc[i];
          loc_new[m_uiDof*1+i] += 0.5*loc[i];
        } else {
          loc_new[m_uiDof*0+i] += loc[i];
        }

        loc_new[m_uiDof*1+i] += loc[m_uiDof + i];

        if ( hangingMask & NODE_2 ) {
          loc_new[m_uiDof*0+i] += 0.25*loc[2*m_uiDof + i];
          loc_new[m_uiDof*1+i] += 0.25*loc[2*m_uiDof + i];
          loc_new[m_uiDof*2+i] += 0.25*loc[2*m_uiDof + i];
          loc_new[m_uiDof*3+i] += 0.25*loc[2*m_uiDof + i];
        } else {
          loc_new[m_uiDof*2+i] += loc[2*m_uiDof + i];
        }

        if ( hangingMask & NODE_3 ) {
          loc_new[m_uiDof*1+i] += 0.5*loc[3*m_uiDof + i];
          loc_new[m_uiDof*3+i] += 0.5*loc[3*m_uiDof + i];
        } else {
          loc_new[m_uiDof*3+i] += loc[3*m_uiDof + i];
        }

        if ( hangingMask & NODE_4 ) {
          loc_new[m_uiDof*0+i] += 0.25*loc[4*m_uiDof + i];
          loc_new[m_uiDof*1+i] += 0.25*loc[4*m_uiDof + i];
          loc_new[m_uiDof*4+i] += 0.25*loc[4*m_uiDof + i];
          loc_new[m_uiDof*5+i] += 0.25*loc[4*m_uiDof + i];
        } else {
          loc_new[m_uiDof*4+i] += loc[4*m_uiDof + i];
        }

        if ( hangingMask & NODE_5 ) {
          loc_new[m_uiDof*1+i] += 0.5*loc[5*m_uiDof + i];
          loc_new[m_uiDof*5+i] += 0.5*loc[5*m_uiDof + i];
        } else {
          loc_new[m_uiDof*5+i] += loc[5*m_uiDof + i];
        }
          
        loc_new[m_uiDof*6+i] += loc[6*m_uiDof + i];
        
        if ( hangingMask & NODE_7 ) {
          loc_new[m_uiDof*1+i] += 0.25*loc[7*m_uiDof + i];
          loc_new[m_uiDof*3+i] += 0.25*loc[7*m_uiDof + i];
          loc_new[m_uiDof*5+i] += 0.25*loc[7*m_uiDof + i];
          loc_new[m_uiDof*7+i] += 0.25*loc[7*m_uiDof + i];
        } else {
          loc_new[m_uiDof*7+i] += loc[7*m_uiDof + i];
        }
      }
      break;
    case 2:
      // 2,5 are not hanging
		  for (size_t i = 0; i < m_uiDof; i++) {
        if ( hangingMask & NODE_0 ) {
          loc_new[m_uiDof*0+i] += 0.5*loc[i];
          loc_new[m_uiDof*2+i] += 0.5*loc[i];
        } else {
          loc_new[m_uiDof*0+i] += loc[i];
        }

        if ( hangingMask & NODE_1 ) {
          loc_new[m_uiDof*0+i] += 0.25*loc[m_uiDof + i];
          loc_new[m_uiDof*1+i] += 0.25*loc[m_uiDof + i];
          loc_new[m_uiDof*2+i] += 0.25*loc[m_uiDof + i];
          loc_new[m_uiDof*3+i] += 0.25*loc[m_uiDof + i];
        } else {
          loc_new[m_uiDof*1+i] += loc[m_uiDof + i];
        }

        loc_new[m_uiDof*2+i] += loc[2*m_uiDof + i];

        if ( hangingMask & NODE_3 ) {
          loc_new[m_uiDof*2+i] += 0.5*loc[3*m_uiDof + i];
          loc_new[m_uiDof*3+i] += 0.5*loc[3*m_uiDof + i];
        } else {
          loc_new[m_uiDof*3+i] += loc[3*m_uiDof + i];
        }

        if ( hangingMask & NODE_4 ) {
          loc_new[m_uiDof*0+i] += 0.25*loc[4*m_uiDof + i];
          loc_new[m_uiDof*2+i] += 0.25*loc[4*m_uiDof + i];
          loc_new[m_uiDof*4+i] += 0.25*loc[4*m_uiDof + i];
          loc_new[m_uiDof*6+i] += 0.25*loc[4*m_uiDof + i];
        } else {
          loc_new[m_uiDof*4+i] += loc[4*m_uiDof + i];
        }

        loc_new[m_uiDof*5+i] += loc[5*m_uiDof + i];
          
        if ( hangingMask & NODE_6 ) {
          loc_new[m_uiDof*2+i] += 0.5*loc[6*m_uiDof + i];
          loc_new[m_uiDof*6+i] += 0.5*loc[6*m_uiDof + i];
        } else {
          loc_new[m_uiDof*6+i] += loc[6*m_uiDof + i];
        }
        
        if ( hangingMask & NODE_7 ) {
          loc_new[m_uiDof*2+i] += 0.25*loc[7*m_uiDof + i];
          loc_new[m_uiDof*3+i] += 0.25*loc[7*m_uiDof + i];
          loc_new[m_uiDof*6+i] += 0.25*loc[7*m_uiDof + i];
          loc_new[m_uiDof*7+i] += 0.25*loc[7*m_uiDof + i];
        } else {
          loc_new[m_uiDof*7+i] += loc[7*m_uiDof + i];
        }
      }
      break;
    case 3:
      // 3,4 are not hanging
		  for (size_t i = 0; i < m_uiDof; i++) {
        if ( hangingMask & NODE_0 ) {
          loc_new[m_uiDof*0+i] += 0.25*loc[i];
          loc_new[m_uiDof*1+i] += 0.25*loc[i];
          loc_new[m_uiDof*2+i] += 0.25*loc[i];
          loc_new[m_uiDof*3+i] += 0.25*loc[i];
        } else {
          loc_new[m_uiDof*0+i] += loc[i];
        }

        if ( hangingMask & NODE_1 ) {
          loc_new[m_uiDof*1+i] += 0.5*loc[m_uiDof + i];
          loc_new[m_uiDof*3+i] += 0.5*loc[m_uiDof + i];
        } else {
          loc_new[m_uiDof*1+i] += loc[m_uiDof + i];
        }

        if ( hangingMask & NODE_2 ) {
          loc_new[m_uiDof*2+i] += 0.5*loc[2*m_uiDof + i];
          loc_new[m_uiDof*3+i] += 0.5*loc[2*m_uiDof + i];
        } else {
          loc_new[m_uiDof*2+i] += loc[2*m_uiDof + i];
        }

        loc_new[m_uiDof*3+i] += loc[3*m_uiDof + i];
        loc_new[m_uiDof*4+i] += loc[4*m_uiDof + i];

        if ( hangingMask & NODE_5 ) {
          loc_new[m_uiDof*1+i] += 0.25*loc[5*m_uiDof + i];
          loc_new[m_uiDof*3+i] += 0.25*loc[5*m_uiDof + i];
          loc_new[m_uiDof*5+i] += 0.25*loc[5*m_uiDof + i];
          loc_new[m_uiDof*7+i] += 0.25*loc[5*m_uiDof + i];
        } else {
          loc_new[m_uiDof*5+i] += loc[5*m_uiDof + i];
        }
        if ( hangingMask & NODE_6 ) {
          loc_new[m_uiDof*2+i] += 0.25*loc[6*m_uiDof + i];
          loc_new[m_uiDof*3+i] += 0.25*loc[6*m_uiDof + i];
          loc_new[m_uiDof*6+i] += 0.25*loc[6*m_uiDof + i];
          loc_new[m_uiDof*7+i] += 0.25*loc[6*m_uiDof + i];
        } else {
          loc_new[m_uiDof*6+i] += loc[6*m_uiDof + i];
        }
        if ( hangingMask & NODE_7 ) {
          loc_new[m_uiDof*3+i] += 0.5*loc[7*m_uiDof + i];
          loc_new[m_uiDof*7+i] += 0.5*loc[7*m_uiDof + i];
        } else {
          loc_new[m_uiDof*7+i] += loc[7*m_uiDof + i];
        }
      }
      break;
    case 4:
		  // 4,3 are not hanging
      for (size_t i = 0; i < m_uiDof; i++) {
        if ( hangingMask & NODE_0 ) {
          loc_new[m_uiDof*0+i] += 0.5*loc[i];
          loc_new[m_uiDof*4+i] += 0.5*loc[i];
        } else {
          loc_new[m_uiDof*0+i] += loc[i];
        }
        if ( hangingMask & NODE_1 ) {
          loc_new[m_uiDof*0+i] += 0.25*loc[m_uiDof + i];
          loc_new[m_uiDof*1+i] += 0.25*loc[m_uiDof + i];
          loc_new[m_uiDof*4+i] += 0.25*loc[m_uiDof + i];
          loc_new[m_uiDof*5+i] += 0.25*loc[m_uiDof + i];
        } else {
          loc_new[m_uiDof*1+i] += loc[m_uiDof + i];
        }

        if ( hangingMask & NODE_2 ) {
          loc_new[m_uiDof*0+i] += 0.25*loc[2*m_uiDof + i];
          loc_new[m_uiDof*2+i] += 0.25*loc[2*m_uiDof + i];
          loc_new[m_uiDof*4+i] += 0.25*loc[2*m_uiDof + i];
          loc_new[m_uiDof*6+i] += 0.25*loc[2*m_uiDof + i];
        } else {
          loc_new[m_uiDof*2+i] += loc[2*m_uiDof + i];
        }

        loc_new[m_uiDof*3+i] += loc[3*m_uiDof + i];
        loc_new[m_uiDof*4+i] += loc[4*m_uiDof + i];
        
        if ( hangingMask & NODE_5 ) {
          loc_new[m_uiDof*4+i] += 0.5*loc[5*m_uiDof + i];
          loc_new[m_uiDof*5+i] += 0.5*loc[5*m_uiDof + i];
        } else {
          loc_new[m_uiDof*5+i] += loc[5*m_uiDof + i];
        }
        if ( hangingMask & NODE_6 ) {
          loc_new[m_uiDof*4+i] += 0.5*loc[6*m_uiDof + i];
          loc_new[m_uiDof*6+i] += 0.5*loc[6*m_uiDof + i];
        } else {
          loc_new[m_uiDof*6+i] += loc[6*m_uiDof + i];
        }
        if ( hangingMask & NODE_7 ) {
          loc_new[m_uiDof*4+i] += 0.25*loc[7*m_uiDof + i];
          loc_new[m_uiDof*5+i] += 0.25*loc[7*m_uiDof + i];
          loc_new[m_uiDof*6+i] += 0.25*loc[7*m_uiDof + i];
          loc_new[m_uiDof*7+i] += 0.25*loc[7*m_uiDof + i];
        } else {
          loc_new[m_uiDof*7+i] += loc[7*m_uiDof + i];
        }
      }
      break;
    case 5:
      // 5,2 are not hanging
      for (size_t i = 0; i < m_uiDof; i++) {
        if ( hangingMask & NODE_0 ) {
          loc_new[m_uiDof*0+i] += 0.25*loc[i];
          loc_new[m_uiDof*1+i] += 0.25*loc[i];
          loc_new[m_uiDof*4+i] += 0.25*loc[i];
          loc_new[m_uiDof*5+i] += 0.25*loc[i];
        } else {
          loc_new[m_uiDof*0+i] += loc[i];
        }
        if ( hangingMask & NODE_1 ) {
          loc_new[m_uiDof*1+i] += 0.5*loc[m_uiDof + i];
          loc_new[m_uiDof*5+i] += 0.5*loc[m_uiDof + i];
        } else {
          loc_new[m_uiDof*1+i] += loc[m_uiDof + i];
        }
        loc_new[m_uiDof*2+i] += loc[2*m_uiDof + i];
        if ( hangingMask & NODE_3 ) {
          loc_new[m_uiDof*1+i] += 0.25*loc[3*m_uiDof + i];
          loc_new[m_uiDof*3+i] += 0.25*loc[3*m_uiDof + i];
          loc_new[m_uiDof*5+i] += 0.25*loc[3*m_uiDof + i];
          loc_new[m_uiDof*7+i] += 0.25*loc[3*m_uiDof + i];
        } else {
          loc_new[m_uiDof*3+i] += loc[3*m_uiDof + i];
        }
        if ( hangingMask & NODE_4 ) {
          loc_new[m_uiDof*4+i] += 0.5*loc[4*m_uiDof + i];
          loc_new[m_uiDof*5+i] += 0.5*loc[4*m_uiDof + i];
        } else {
          loc_new[m_uiDof*4+i] += loc[4*m_uiDof + i];
        }
        loc_new[m_uiDof*5+i] += loc[5*m_uiDof + i];
        if ( hangingMask & NODE_6 ) {
          loc_new[m_uiDof*4+i] += 0.25*loc[6*m_uiDof + i];
          loc_new[m_uiDof*5+i] += 0.25*loc[6*m_uiDof + i];
          loc_new[m_uiDof*6+i] += 0.25*loc[6*m_uiDof + i];
          loc_new[m_uiDof*7+i] += 0.25*loc[6*m_uiDof + i];
        } else {
          loc_new[m_uiDof*6+i] += loc[6*m_uiDof + i];
        }
        if ( hangingMask & NODE_7 ) {
          loc_new[m_uiDof*5+i] += 0.5*loc[7*m_uiDof + i];
          loc_new[m_uiDof*7+i] += 0.5*loc[7*m_uiDof + i];
        } else {
          loc_new[m_uiDof*7+i] += loc[7*m_uiDof + i];
        }
      }
      break;
    case 6:
      // 6,1 are not hanging
      for (size_t i = 0; i < m_uiDof; i++) {
        if ( hangingMask & NODE_0 ) {
          loc_new[m_uiDof*0+i] += 0.25*loc[i];
          loc_new[m_uiDof*1+i] += 0.25*loc[i];
          loc_new[m_uiDof*4+i] += 0.25*loc[i];
          loc_new[m_uiDof*5+i] += 0.25*loc[i];
        } else {
          loc_new[m_uiDof*0+i] += loc[i];
        }
        loc_new[m_uiDof*1+i] += loc[m_uiDof + i];
        if ( hangingMask & NODE_2 ) {
          loc_new[m_uiDof*2+i] += 0.5*loc[2*m_uiDof + i];
          loc_new[m_uiDof*6+i] += 0.5*loc[2*m_uiDof + i];
        } else {
          loc_new[m_uiDof*2+i] += loc[2*m_uiDof + i];
        }
        if ( hangingMask & NODE_3 ) {
          loc_new[m_uiDof*2+i] += 0.25*loc[3*m_uiDof + i];
          loc_new[m_uiDof*3+i] += 0.25*loc[3*m_uiDof + i];
          loc_new[m_uiDof*6+i] += 0.25*loc[3*m_uiDof + i];
          loc_new[m_uiDof*7+i] += 0.25*loc[3*m_uiDof + i];
        } else {
          loc_new[m_uiDof*3+i] += loc[3*m_uiDof + i];
        }
        if ( hangingMask & NODE_4 ) {
          loc_new[m_uiDof*4+i] += 0.5*loc[4*m_uiDof + i];
          loc_new[m_uiDof*6+i] += 0.5*loc[4*m_uiDof + i];
        } else {
          loc_new[m_uiDof*4+i] += loc[4*m_uiDof + i];
        }
        if ( hangingMask & NODE_5 ) {
          loc_new[m_uiDof*4+i] += 0.25*loc[5*m_uiDof + i];
          loc_new[m_uiDof*5+i] += 0.25*loc[5*m_uiDof + i];
          loc_new[m_uiDof*6+i] += 0.25*loc[5*m_uiDof + i];
          loc_new[m_uiDof*7+i] += 0.25*loc[5*m_uiDof + i];
        } else {
          loc_new[m_uiDof*5+i] += loc[5*m_uiDof + i];
        }
        loc_new[m_uiDof*6+i] += loc[6*m_uiDof + i];
        if ( hangingMask & NODE_7 ) {
          loc_new[m_uiDof*6+i] += 0.5*loc[7*m_uiDof + i];
          loc_new[m_uiDof*7+i] += 0.5*loc[7*m_uiDof + i];
        } else {
          loc_new[m_uiDof*7+i] += loc[7*m_uiDof + i];
        }
      }
      break;
    case 7:
      // 7,0 are not hanging
      for (size_t i = 0; i < m_uiDof; i++) {
        loc_new[m_uiDof*0+i] += loc[i];
        if ( hangingMask & NODE_1 ) {
          loc_new[m_uiDof*1+i] += 0.25*loc[m_uiDof + i];
          loc_new[m_uiDof*3+i] += 0.25*loc[m_uiDof + i];
          loc_new[m_uiDof*5+i] += 0.25*loc[m_uiDof + i];
          loc_new[m_uiDof*7+i] += 0.25*loc[m_uiDof + i];
        } else {
          loc_new[m_uiDof*1+i] += loc[m_uiDof + i];
        }
        if ( hangingMask & NODE_2 ) {
          loc_new[m_uiDof*2+i] += 0.25*loc[2*m_uiDof + i];
          loc_new[m_uiDof*3+i] += 0.25*loc[2*m_uiDof + i];
          loc_new[m_uiDof*6+i] += 0.25*loc[2*m_uiDof + i];
          loc_new[m_uiDof*7+i] += 0.25*loc[2*m_uiDof + i];
        } else {
          loc_new[m_uiDof*2+i] += loc[2*m_uiDof + i];
        }
        if ( hangingMask & NODE_3 ) {
          loc_new[m_uiDof*3+i] += 0.5*loc[3*m_uiDof + i];
          loc_new[m_uiDof*7+i] += 0.5*loc[3*m_uiDof + i];
        } else {
          loc_new[m_uiDof*3+i] += loc[3*m_uiDof + i];
        }

        if ( hangingMask & NODE_4 ) {
          loc_new[m_uiDof*4+i] += 0.25*loc[4*m_uiDof + i];
          loc_new[m_uiDof*5+i] += 0.25*loc[4*m_uiDof + i];
          loc_new[m_uiDof*6+i] += 0.25*loc[4*m_uiDof + i];
          loc_new[m_uiDof*7+i] += 0.25*loc[4*m_uiDof + i];
        } else {
          loc_new[m_uiDof*4+i] += loc[4*m_uiDof + i];
        }
        if ( hangingMask & NODE_5 ) {
          loc_new[m_uiDof*5+i] += 0.5*loc[5*m_uiDof + i];
          loc_new[m_uiDof*7+i] += 0.5*loc[5*m_uiDof + i];
        } else {
          loc_new[m_uiDof*5+i] += loc[5*m_uiDof + i];
        }
        if ( hangingMask & NODE_6 ) {
          loc_new[m_uiDof*6+i] += 0.5*loc[6*m_uiDof + i];
          loc_new[m_uiDof*7+i] += 0.5*loc[6*m_uiDof + i];
        } else {
          loc_new[m_uiDof*6+i] += loc[6*m_uiDof + i];
        }
        loc_new[m_uiDof*7+i] += loc[7*m_uiDof + i];
      }
    break;
    default:
			std::cout<<"in loc_to_glo: incorrect child num = " << chNum << std::endl;
			assert(false);
      break;
  } // switch chNum
} // loc_to_glo


void interp_local_to_global_matrix(PetscScalar* Ke, std::vector<ot::MatRecord>& out, ot::DA* da) {
  unsigned int m_uiDof = 1;
  unsigned int m = m_uiDof * 8;  // matrix size in one dimension

  unsigned int idx[8];
  da->getNodeIndices(idx);

  // get continuous columns of Ke (row major -> column major)
  PetscScalar* Ke_cols = new PetscScalar[m*m];
  for (unsigned int i = 0; i < m; i++) {
    for (unsigned int j = 0; j < m; j++) {
      Ke_cols[i*m + j] = Ke[i + j*m];
    }
  }

  PetscScalar* interp = new PetscScalar[m];
  for (unsigned int i = 0; i < m; i++) {  // column
    // TODO change interp function to be = and not += so we dont have to zero out first
    for (unsigned int j = 0; j < m; j++) {
      interp[j] = 0.0;
    }

    interp_local_to_global_still_local(Ke_cols + m*i, interp, da);

    // convert interp to m ot::MatRecords
    for (unsigned int j = 0; j < m; j++) {  // row
      ot::MatRecord mr;
      mr.rowIdx = idx[j / m_uiDof];
      mr.colIdx = idx[i / m_uiDof];
      mr.rowDim = j % m_uiDof;
      mr.colDim = i % m_uiDof;
      mr.val = interp[j];
      out.push_back(mr);
    }
  }

  delete[] Ke_cols;
  delete[] interp;
}

