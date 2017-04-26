#pragma once

#include <cmath>
#include "dendro.h"
#include "oda.h"

void addAnalyticalSolution(ot::DA* da, const double* problemSize, Vec u_vec, Vec* output_vec_out, int ndof, double ts);

double calc_l2_error(ot::DA* da, const double* problemSize, Vec u_vec, int ndof, double ts);
