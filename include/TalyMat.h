#pragma once

#include "feMatrix.h"

// TalyFEM includes
#include <Grid/grid_types/grid.h>
#include <Grid/gridfield.h>
#include <Grid/elem.h>
#include <Grid/elem_types/elem3dhexahedral.h>
#include <Grid/femelm.h>
#include <DataStructure/zeroarray.h>
#include <DataStructure/zeromatrix.h>

template <typename Equation, typename NodeData>
class TalyMatrix : public feMatrix< TalyMatrix<Equation, NodeData> >
{
  typedef feMatrix< TalyMatrix<Equation, NodeData> > Parent;

 public:
  TalyMatrix(feMat::daType da)
    : Parent(da), taly_grid_(), taly_elem_(NULL), taly_fe_(&taly_grid_, TALYFEMLIB::BASIS_ALL)
  {
    // 8 nodes, 1 element
    taly_grid_.redimArrays(8, 1);
    for (int i = 0; i < 8; i++) {
      taly_grid_.node_array_[i] = new TALYFEMLIB::NODE();
    }

    taly_elem_ = new TALYFEMLIB::ELEM3dHexahedral();
    taly_grid_.elm_array_[0] = taly_elem_;

    int node_id_array[8] = {
      0, 1, 2, 3, 4, 5, 6, 7
    };
    taly_elem_->redim(8, node_id_array);

    taly_gf_.redimGrid(&taly_grid_);
    taly_gf_.redimNodeData();

    taly_eq_.p_data_ = &taly_gf_;
    taly_eq_.p_grid_ = &taly_grid_;

    Ae_.redim(8, 8);
    be_.redim(8);
  }

  inline bool initStencils() {}

  inline bool ElementalMatVec(int i, int j, int k, PetscScalar ***in, PetscScalar ***out, double scale) {
    assert(false);
  }

  inline bool ElementalMatVec(unsigned int index, PetscScalar *in, PetscScalar *out, double scale) {
    assert(false);
  }

  inline bool ElementalMatVec(PetscScalar* in_local, PetscScalar* out_local, PetscScalar* coords, double scale)
  {
    const int ndof = 1;  // TODO

    // update node coordinates and values
    for (unsigned int i = 0; i < 8; i++) {
      taly_grid_.node_array_[i]->setCoor(coords[i*3], coords[i*3+1], coords[i*3+2]);
      //std::cout << "Node [" << i << "] coords: " << coords[i*3] << ", " << coords[i*3+1] << ", " << coords[i*3+2] << "\n";

      for (unsigned int dof = 0; dof < ndof; dof++) {
        taly_gf_.GetNodeData(i).value(dof) = in_local[i*ndof+dof];
      }
    }
    //std::cout << "\n";

    Ae_.fill(0.0);
    be_.fill(0.0);

    // sum integrands over all gauss points
    taly_fe_.refill(taly_elem_, TALYFEMLIB::BASIS_LINEAR, 0);
    while (taly_fe_.next_itg_pt()) {
      taly_eq_.Integrands(taly_fe_, Ae_, be_);
    }

    /*for (unsigned int y = 0; y < 8; y++) {
      for (unsigned int x = 0; x < 8; x++) {
        std::cout << "Ae(" << x << ", " << y << ") = " << Ae_(x, y) << "  ";
      }
      std::cout << "\n";
    }*/

    // copy Ae/be to dendro structures
    for (int k = 0; k < 8; k++) {
      for (int j=0; j<8; j++) {
        out_local[ndof * k] += Ae_(k, j) * in_local[ndof * j];  // TODO * scale?
      }
    }
    /*for (int i = 0; i < 8; i++) {
      std::cout << "in_local[" << i << "] = " << in_local[i] << "   (coords: " << coords[i*3] << ", " << coords[i*3+1] << ", " << coords[i*3+2] << ")\n";
    }
    std::cout << "\n";*/

    // TODO be is ignored
    return true;
  }

  bool preMatVec() {}
  bool postMatVec() {}

  // void Integrands(const FEMElm& fe, ZEROMATRIX<double>& Ae, ZEROMATRIX<double>& be);

 private:
  TALYFEMLIB::GRID taly_grid_;
  TALYFEMLIB::GridField<NodeData> taly_gf_;
  Equation taly_eq_;  // not fully initialized, just use to access Integrands

  TALYFEMLIB::ELEM* taly_elem_;
  TALYFEMLIB::FEMElm taly_fe_;
  TALYFEMLIB::ZeroMatrix<double> Ae_;
  TALYFEMLIB::ZEROARRAY<double> be_;
};
