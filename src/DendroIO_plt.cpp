#include "DendroIO.h"
#include "mpi.h"
#include "TecplotIO_ascii.h"
#include "interp.h"

#include <fstream>
#include <sstream>

void octree2PLT(ot::DA* da, Vec vec, int ndof, std::string file_name)
{
  unsigned int nsd = 3;
  ElemType elem_type = kElem3dHexahedral;
  unsigned int nodes_per_elem = 8;
  unsigned int n_elements = da->getElementSize();
  unsigned int n_nodes = n_elements * nodes_per_elem;

{
  int size;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  assert(size == 1);
}

  PetscScalar* u = NULL;
  da->vecGetBuffer(vec, u, false, false, true, ndof);

  // TODO - is this necessary?
  da->ReadFromGhostsBegin<PetscScalar>(u, ndof);
  da->ReadFromGhostsEnd<PetscScalar>(u);

  TecplotWriterASCII w;
  w.open(file_name.c_str(), false);

  TecplotHeader header;
  header.title = file_name;
  header.variables = { "x", "y", "z" };

  // add variables based on ndof
  for (int i = 0; i < ndof; i++) {
    header.variables.push_back(std::string("u") + std::to_string(i));
  }

  w.write_header(header);

  TecplotZone zone;
  zone.num_nodes = n_nodes;
  zone.num_elements = n_elements;
  zone.format = kFiniteElementPoint;
  zone.elem_type = elem_type;
  w.write_zone(zone);

  double* node_data_loc = new double[8*ndof];

  int maxD = da->getMaxDepth();
  double hx = 1.0 / ((double)(1 << (maxD-1)));  // TODO do hy/hz to support varied problem sizes
  for ( da->init<ot::DA_FLAGS::ALL>(), da->init<ot::DA_FLAGS::WRITABLE>(); da->curr() < da->end<ot::DA_FLAGS::ALL>(); da->next<ot::DA_FLAGS::ALL>()) {
    unsigned int i = da->curr();
    Point pt = da->getCurrentOffset();  // TODO this is wrong, need to multiply by xFac

    unsigned int lev = da->getLevel(i);
    unsigned int width = ((1 << (maxD - lev)));
    double coords[8][3] = {
      {pt.x(), pt.y(), pt.z()},
      {pt.x()+width, pt.y(), pt.z()},
      {pt.x()+width, pt.y()+width, pt.z()},
      {pt.x(), pt.y()+width, pt.z()},

      {pt.x(), pt.y(), pt.z()+width},
      {pt.x()+width, pt.y(), pt.z()+width},
      {pt.x()+width, pt.y()+width, pt.z()+width},
      {pt.x(), pt.y()+width, pt.z()+width},
    };


    unsigned int node_idx[8];
    da->getNodeIndices(node_idx); 
    interp_global_to_local(u, node_data_loc, da, ndof);

    std::vector<double> data(header.variables.size());
    for (int j = 0; j < 8; j++) {
      double p[3] = { coords[j][0] * hx, coords[j][1] * hx, coords[j][2] * hx };
      //double val = node_data_loc[j];
      for (int d = 0; d < nsd; d++) {
        data[d] = p[d];
      }
      for (int d = 0; d < ndof; d++) {
        data[nsd + d] = node_data_loc[j*ndof+d];
      }

      w.write_node(data.data());
    }
  }

  delete[] node_data_loc;

  for (unsigned int i = 0; i < n_elements; i++) {
    PhysicalNodeID conn[8];
    for (int node_idx = 0; node_idx < nodes_per_elem; node_idx++) {
      conn[node_idx] = 1 + i * nodes_per_elem + node_idx;
    }
    w.write_elem(conn);
  }

  da->vecRestoreBuffer(vec, u, false, false, true, ndof);
}

