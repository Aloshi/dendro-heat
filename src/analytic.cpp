#include "analytic.h"

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

        // no need to map here, vecGetBuffer has taken care of it
        //PetscInt row = localToGlobal[idx[i]];
        double u_analytical = calc_u_analytical(xx[i], yy[i], zz[i], ts);
        output_vec_data[idx[i]*2] = u_vec_data[idx[i]];
        output_vec_data[idx[i]*2 + 1] = u_analytical;
        double err_r = fabs(u_vec_data[idx[i]] - u_analytical) / fabs(u_analytical);
        if (err_r>err && fabs(u_analytical)>1e-6)
          err = err_r;
      }
    }
  } // loop over elements

  //std::cout<<" maximum error: "<<err<<std::endl;
  da->vecRestoreBuffer(*output_vec, output_vec_data, false, false, false, ndof+1);
  da->vecRestoreBuffer(u_vec, u_vec_data, false, false, false, ndof);
}
