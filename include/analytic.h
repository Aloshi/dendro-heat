#pragma once

#include <cmath>
#include "dendro.h"
#include "oda.h"

void addAnalyticalSolution(ot::DA* da, const double* problemSize, Vec u_vec, Vec* output_vec_out, int ndof, double ts);
