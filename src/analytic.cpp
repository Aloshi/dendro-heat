#include "analytic.h"

#include <HTEquation.h>
#include <Grid/elem-types.h>
#include <interp.h>
#include <feMatrix.h>

double calc_u_analytical(double x, double y, double z, double ts) {
  return sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*exp((-ts));
}

void addAnalyticalSolution(ot::DA* da, const double* problemSize, Vec u_vec, Vec* output_vec, int ndof, double ts)
{
  unsigned int maxD = da->getMaxDepth();
  unsigned int lev;
  double hx, hy, hz, dist, half;
  Point pt;

  double xFac = problemSize[0]/((double)(1<<(maxD-1))); 
  double yFac = problemSize[1]/((double)(1<<(maxD-1)));
  double zFac = problemSize[2]/((double)(1<<(maxD-1)));
  double xx[8], yy[8], zz[8];
  unsigned int idx[8];

  if (!da->computedLocalToGlobal())
    da->computeLocalToGlobalMappings();
  assert(da->computedLocalToGlobal());

  da->createVector(*output_vec, false, false, ndof+1);

  PetscScalar* output_vec_data;
  da->vecGetBuffer(*output_vec, output_vec_data, false, false, false, ndof+1);

  PetscScalar* u_vec_data;
  da->vecGetBuffer(u_vec, u_vec_data, false, false, false, ndof);

  double err = 0;

  const DendroIntL* localToGlobal = da->getLocalToGlobalMap();
  //for ( da->init<ot::DA_FLAGS::WRITABLE>(); da->curr() < da->end<ot::DA_FLAGS::WRITABLE>(); da->next<ot::DA_FLAGS::WRITABLE>() ) {
  for ( da->init<ot::DA_FLAGS::ALL>(); da->curr() < da->end<ot::DA_FLAGS::ALL>(); da->next<ot::DA_FLAGS::ALL>() ) {
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

        // no need to map here, vecGetBuffer has taken care of it
        //PetscInt row = localToGlobal[idx[i]];
        double u_analytical = calc_u_analytical(xx[i], yy[i], zz[i], ts);
        output_vec_data[idx[i]*2] = u_vec_data[idx[i]];
        output_vec_data[idx[i]*2 + 1] = u_analytical;
      }
    }
  } // loop over elements

  da->vecRestoreBuffer(*output_vec, output_vec_data, false, false, false, ndof+1);
  da->vecRestoreBuffer(u_vec, u_vec_data, false, false, false, ndof);
}

double calc_l2_error(ot::DA* da, const double* problemSize, Vec u_vec, int ndof, double ts)
{
  using namespace TALYFEMLIB;
  GRID grid;
  grid.redimArrays(8, 1);
  for (int i = 0; i < 8; i++) {
    grid.node_array_[i] = new NODE();
  }

  ELEM* elem = new ELEM3dHexahedral();
  grid.elm_array_[0] = elem;

  static const int node_id_array[8] = {
    0, 1, 2, 3, 4, 5, 6, 7
  };
  elem->redim(8, node_id_array);

  GridField<HTNodeData> gf;
  gf.redimGrid(&grid);
  gf.redimNodeData();
  
  FEMElm fe(&grid, BASIS_FIRST_DERIVATIVE | BASIS_POSITION);

  double coords[8*3];
  double local_in_dendro[8 * ndof];
  double local_in_taly[8 * ndof];

  const unsigned int maxD = da->getMaxDepth();
  const double xFac = problemSize[0]/((double)(1<<(maxD-1))); 
  const double yFac = problemSize[1]/((double)(1<<(maxD-1)));
  const double zFac = problemSize[2]/((double)(1<<(maxD-1)));

  PetscScalar* u_vec_data;
  da->vecGetBuffer(u_vec, u_vec_data, false, false, true, ndof);
  da->ReadFromGhostsBegin<PetscScalar>(u_vec_data, ndof);
  da->ReadFromGhostsEnd<PetscScalar>(u_vec_data);

  double l2_error = 0.0;
  double a_norm = 0.0, c_norm = 0.0;
  //for ( da->init<ot::DA_FLAGS::ALL>(); da->curr() < da->end<ot::DA_FLAGS::ALL>(); da->next<ot::DA_FLAGS::ALL>() ) {
  for ( da->init<ot::DA_FLAGS::WRITABLE>(); da->curr() < da->end<ot::DA_FLAGS::WRITABLE>(); da->next<ot::DA_FLAGS::WRITABLE>() ) {
    int lev = da->getLevel(da->curr());
    Point h(xFac*(1<<(maxD - lev)), yFac*(1<<(maxD - lev)), zFac*(1<<(maxD - lev)));
    Point pt = da->getCurrentOffset();
    pt.x() *= xFac;
    pt.y() *= yFac;
    pt.z() *= zFac;

    // build node coordinates (fill coords)
    build_taly_coordinates(coords, pt, h);

    // interpolate (in -> node_data_temp)
    interp_global_to_local(u_vec_data, local_in_dendro, da, ndof);

    // map from dendro order to taly order
    dendro_to_taly(local_in_taly, local_in_dendro, sizeof(local_in_taly[0])*ndof);

    // update node coordinates and values
    for (unsigned int i = 0; i < 8; i++) {
      grid.node_array_[i]->setCoor(coords[i*3], coords[i*3+1], coords[i*3+2]);

      for (unsigned int dof = 0; dof < ndof; dof++) {
        gf.GetNodeData(i).value(dof) = local_in_taly[i*ndof+dof];
      }
    }

    fe.refill(elem, BASIS_LINEAR, 0);
    while (fe.next_itg_pt()) {
      double val_c = gf.valueFEM(fe, 0);
      double val_a = calc_u_analytical(fe.position().x(), fe.position().y(), fe.position().z(), ts);
      c_norm += val_c * val_c * fe.detJxW();
      a_norm += val_a * val_a * fe.detJxW();
      l2_error += (val_c - val_a) * (val_c - val_a) * fe.detJxW();
    }
  }

  da->vecRestoreBuffer(u_vec, u_vec_data, false, false, true, ndof);

  double all_err, all_a_norm, all_c_norm;
  MPI_Allreduce(&l2_error, &all_err, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&a_norm, &all_a_norm, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&c_norm, &all_c_norm, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  /*int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank)
    std::cout << "a norm: " << all_a_norm << ", c norm: " << all_c_norm << "\n";*/

  return sqrt(all_err);
}

