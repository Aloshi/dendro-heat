#pragma once

#include "parUtils.h"
#include "octUtils.h"
#include "dendro.h"

void interp_global_to_local(PetscScalar* glo, PetscScalar* /* __restrict*/ loc, ot::DA* m_octDA);
void interp_local_to_global(PetscScalar* /*__restrict*/ loc, PetscScalar* glo, ot::DA* da);
