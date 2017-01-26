#pragma once

#include "parUtils.h"
#include "octUtils.h"
#include "dendro.h"

void interp_global_to_local(PetscScalar* glo, PetscScalar* /* __restrict*/ loc, ot::DA* m_octDA, int dof);
void interp_local_to_global(PetscScalar* /*__restrict*/ loc, PetscScalar* glo, ot::DA* da, int dof);
void interp_local_to_global_matrix(PetscScalar* Ke, std::vector<ot::MatRecord>& out, ot::DA* da, int dof);
