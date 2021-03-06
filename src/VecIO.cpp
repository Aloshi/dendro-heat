#include "VecIO.h"

#include "mpi.h"
#include <assert.h>
#include <sstream>
#include <iostream>
#include "TecplotIO_ascii.h"

// get node in an element
int get_node_id(int i, int j, int k, int nx, int ny, int nz, int node_idx) {
  int idx[8][3]={
    {k, j, i},
    {k+1, j, i},
    {k+1, j+1, i},
    {k, j+1, i},

    {k, j, i+1},
    {k+1, j, i+1},
    {k+1, j+1, i+1},
    {k, j+1, i+1},
  };

  return idx[node_idx][2] + idx[node_idx][1] * nx + idx[node_idx][0] * nx * ny;
}

// 3D only
int write_vector(const char* file_prefix, Vec vec, int ndof, DM da)
{
  int mpi_size;
  MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
  if (mpi_size != 1) {
    std::cerr << "Cannot save vector (" << file_prefix << ".plt) on more than one process yet\n";
    return 0;
  }
  
  int x, y, z, m, n, p;
  int mx,my,mz, xne, yne, zne;

  CHKERRQ( DMDAGetCorners(da, &x, &y, &z, &m, &n, &p) ); 
  CHKERRQ( DMDAGetInfo(da,0, &mx, &my, &mz, 0,0,0,0,0,0,0,0,0) ); 

  assert(mx == my && my == mz);
  int Ns = mx - 1;
  double h = 1.0 / ((double)Ns);

  if (x+m == mx) {
    xne=m-1;
  } else {
    xne=m;
  }
  if (y+n == my) {
    yne=n-1;
  } else {
    yne=n;
  }
  if (z+p == mz) {
    zne=p-1;
  } else {
    zne=p;
  }

  unsigned int n_elements = Ns*Ns*Ns;
  unsigned int n_nodes = (Ns+1)*(Ns+1)*(Ns+1);
  unsigned int nsd = 3;
  unsigned int nodes_per_elem = 8;
  ElemType elem_type = kElem3dHexahedral;

  std::vector<PhysicalNodeID> connectivity;
  connectivity.resize(n_elements * nodes_per_elem);

  PetscScalar ***array;
  CHKERRQ(DMDAVecGetArray(da, vec, &array));

  TecplotWriterASCII w;
  std::stringstream ss;
  ss << file_prefix << ".plt";
  w.open(ss.str().c_str(), false);

  TecplotHeader header;
  header.title = file_prefix;
  header.variables = { "x", "y", "z" };
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

  std::vector<double> node_data(header.variables.size());
  unsigned int last_node_id = 0;
  for (int k=z; k < z+p; k++) {
    for (int j=y; j < y+n; j++) {
      for (int i=x; i < x+m; i++) {
        node_data[0] = h*i;
        node_data[1] = h*j;
        node_data[2] = h*k;
        for (int d = 0; d < ndof; d++) {
          node_data[nsd+d] = array[k][j][i*ndof+d];
        }
        //double data[4] = { h*i, h*j, h*k, array[k][j][i] };
        w.write_node(node_data.data());
        assert(get_node_id(i, j, k, Ns+1, Ns+1, Ns+1, 0) == last_node_id++);
      } // end i
    } // end j
  } // end k

  CHKERRQ( DMDAVecRestoreArray ( da, vec, &array ) );

  for (int k=z; k < zne; k++) {
    for (int j=y; j < yne; j++) {
      for (int i=x; i < xne; i++) {
        PhysicalNodeID conn[8];
        for (int node_idx = 0; node_idx < nodes_per_elem; node_idx++) {
          conn[node_idx] = get_node_id(i, j, k, Ns+1, Ns+1, Ns+1, node_idx) + 1;
        }
        w.write_elem(conn);
      } // end i
    } // end j
  } // end k

  return 0;
}
