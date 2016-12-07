#pragma once

#include "petscdmda.h"
#include "petscsys.h"

int write_vector(const char* fileName, Vec v, DM da);
