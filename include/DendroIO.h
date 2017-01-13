#pragma once

#include <vector>
#include <string>
#include "TreeNode.h"
#include "oda.h"

void octree2PLT(ot::DA* da, Vec u, std::string file_prefix);
void octree2VTK(ot::DA* da, Vec u, std::string file_prefix);
