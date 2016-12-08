#pragma once

#include <vector>
#include <string>
#include "TreeNode.h"
#include "oda.h"

//void octree2VTK(const std::vector<ot::TreeNode>& nodes, unsigned int rank, const double* u, std::string vtk_file_name);
void octree2VTK(ot::DA& da, unsigned int rank, const double* u, std::string file_name);
