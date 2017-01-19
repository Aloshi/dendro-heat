/*
  Copyright 2014-2016 Baskar Ganapathysubramanian

  This file is part of TALYFem.

  TALYFem is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 2.1 of the
  License, or (at your option) any later version.

  TALYFem is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TALYFem.  If not, see <http://www.gnu.org/licenses/>.
*/
// --- end license text --- //
#ifndef INCLUDE_HTEQUATION_H_
#define INCLUDE_HTEQUATION_H_

#include <FEM/cequation.h>
#include "HTNodeData.h"

class HTEquation : public TALYFEMLIB::CEquation<HTNodeData> {
 public:

  explicit HTEquation()
      : TALYFEMLIB::CEquation<HTNodeData>(false, TALYFEMLIB::kAssembleGaussPoints) {
    K_ = 1.0;
    dt_ = 0.001;  // TODO
  }

  virtual void Solve(double dt, double t) {}

  /**
   * Fills the Ae and be structures with data for a single Gauss point.
   *
   * This uses a Backward Euler scheme to discretize the heat equation.
   *
   * @param fe the element we are assembling a Gauss point from
   * @param Ae the element matrix to put data in
   * @param be the element vector to put data in
   */
  void Integrands(const TALYFEMLIB::FEMElm& fe, TALYFEMLIB::ZeroMatrix<double>& Ae,
                          TALYFEMLIB::ZEROARRAY<double>& be) {
    const int n_dimensions = fe.nsd();  // # of dimensions: 1D, 2D, or 3D
    const int n_basis_functions = fe.nbf();  // # of basis functions
    const double detJxW = fe.detJxW();  // (determinant of J) cross W
    double k_val = K_;  // thermal diffusivity in heat equation

    // ValueFEM works for hermite because HTNodeData is ordered so that the
    // DU_PRE values come right after the U_PRE value.
    // (U_PRE, DU_PRE_1, DU_PRE_2, etc.)
    const double u_pre_curr = p_data_->valueFEM(fe, U);

    // in order to assemble the gauss point, we loop over each pair of basis
    // functions and calculate the individual contributions to the 'Ae' matrix
    // and 'be' vector.
    for (int a = 0; a < n_basis_functions; a++) {
      for (int b = 0; b < n_basis_functions; b++) {
        double M = fe.N(a) * fe.N(b) * detJxW;
        double N = 0;
        for (int k = 0; k < n_dimensions; k++) {
          N += k_val * fe.dN(a, k) * fe.dN(b, k) * detJxW;
        }
        // Add term to the A element matrix
        Ae(a, b) += M / dt_ + N;
      }

      // Add term to the b element vector
      be(a) += fe.N(a) / dt_ * u_pre_curr * detJxW;
    }
  }

  void Integrands_Ae(const TALYFEMLIB::FEMElm& fe, TALYFEMLIB::ZeroMatrix<double>& Ae) {
    const int n_dimensions = fe.nsd();  // # of dimensions: 1D, 2D, or 3D
    const int n_basis_functions = fe.nbf();  // # of basis functions
    const double detJxW = fe.detJxW();  // (determinant of J) cross W
    double k_val = K_;  // thermal diffusivity in heat equation

    // in order to assemble the gauss point, we loop over each pair of basis
    // functions and calculate the individual contributions to the 'Ae' matrix
    // and 'be' vector.
    for (int a = 0; a < n_basis_functions; a++) {
      for (int b = 0; b < n_basis_functions; b++) {
        double M = fe.N(a) * fe.N(b) * detJxW;
        double N = 0;
        for (int k = 0; k < n_dimensions; k++) {
          N += k_val * fe.dN(a, k) * fe.dN(b, k) * detJxW;
        }
        // Add term to the A element matrix
        Ae(a, b) += M / dt_ + N;
      }
    }
  }


  void Integrands_be(const TALYFEMLIB::FEMElm& fe, TALYFEMLIB::ZEROARRAY<double>& be) {
    const int n_basis_functions = fe.nbf();  // # of basis functions
    const double detJxW = fe.detJxW();  // (determinant of J) cross W

    const double u_pre_curr = p_data_->valueFEM(fe, U);

    for (int a = 0; a < n_basis_functions; a++) {
      be(a) += fe.N(a) / dt_ * u_pre_curr * detJxW;
    }
  }

 private:
  double K_;
};

#endif  // INCLUDE_HTEQUATION_H_
