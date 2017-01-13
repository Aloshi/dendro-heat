#pragma once

#include "petscdmda.h"
#include "petscsys.h"

int write_vector(const char* file_prefix, Vec v, DM da);
