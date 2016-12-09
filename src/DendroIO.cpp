#include "DendroIO.h"
#include "mpi.h"
#include "TecplotIO_ascii.h"
#include "feMatrix.h"  // haaaack

#include <fstream>
#include <sstream>

#define VTK_HEXAHEDRON 12

class Aligner : public feMatrix<Aligner> {
  public:
    inline bool ElementalMatVec(int i, int j, int k, PetscScalar ***in, PetscScalar ***out, double scale) { return true; }
    inline bool ElementalMatVec(unsigned int idx, PetscScalar *in, PetscScalar *out, double scale) { return true; }
    inline bool ElementalMatVec(PetscScalar* in_local, PetscScalar* out_local, PetscScalar* coords, double scale) {
      assert(false);
    }

    inline bool GetElementalMatrix(int i, int j, int k, PetscScalar *mat) { return true; }
    inline bool GetElementalMatrix(unsigned int idx, std::vector<ot::MatRecord> &records) { return true; }

    inline bool ElementalMatGetDiagonal(int i, int j, int k, PetscScalar ***diag, double scale) { return true; }
    inline bool ElementalMatGetDiagonal(unsigned int idx, PetscScalar *diag, double scale) { return true; }

    inline bool initStencils() { initOctLut(); return true; }
};

void octree2VTK(ot::DA& da, unsigned int rank, const double* u, std::string file_name)
{
  unsigned int nsd = 3;
  unsigned int ndof = 1;
  ElemType elem_type = kElem3dHexahedral;
  unsigned int nodes_per_elem = 8;
  unsigned int n_elements = da.getElementSize();
  unsigned int n_nodes = n_elements * nodes_per_elem;

{
  int size;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  assert(size == 1);
}

  TecplotWriterASCII w;
  w.open(file_name.c_str(), false);

  TecplotHeader header;
  header.title = file_name;
  header.variables = { "x", "y", "z", "u" };
  w.write_header(header);

  TecplotZone zone;
  zone.num_nodes = n_nodes;
  zone.num_elements = n_elements;
  zone.format = kFiniteElementPoint;
  zone.elem_type = elem_type;
  w.write_zone(zone);

  int maxD = da.getMaxDepth();
  double hx = 1.0 / ((double)(1 << (maxD-1)));
  for ( da.init<ot::DA_FLAGS::ALL>(), da.init<ot::DA_FLAGS::WRITABLE>(); da.curr() < da.end<ot::DA_FLAGS::ALL>(); da.next<ot::DA_FLAGS::ALL>()) {
    unsigned int i = da.curr();
    Point pt = da.getCurrentOffset();

    unsigned int lev = da.getLevel(i);
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
{
    Aligner hack;
    Aligner::stdElemType elemType;
    hack.alignElementAndVertices(&da, elemType, node_idx);
}

    for (int j = 0; j < 8; j++) {
      double p[3] = { coords[j][0] * hx, coords[j][1] * hx, coords[j][2] * hx };
      double val = u[node_idx[j]];
      //double val = u[i];
      double data[4] = { p[0], p[1], p[2], val };
      w.write_node(data);
    }
  }

  for (unsigned int i = 0; i < n_elements; i++) {
    PhysicalNodeID conn[8];
    for (int node_idx = 0; node_idx < nodes_per_elem; node_idx++) {
      conn[node_idx] = 1 + i * nodes_per_elem + node_idx;
    }
    w.write_elem(conn);
  }
}

/*void octree2VTK(ot::DA& da, unsigned int rank, const double* u, std::string vtk_file_name)
{
    static const int NUM_FUNC_VALUES = 1;
    static const int U = 0;

    if (!rank) std::cout << "writing mesh to VTK file: " << vtk_file_name << std::endl;
    std::ostringstream convert;

#ifdef  HILBERT_ORDERING
    convert << vtk_file_name << "_H_" <<rank << ".vtk";
#else
    convert << vtk_file_name << "_M_" << rank << ".vtk";
#endif

    vtk_file_name = convert.str();

    std::ofstream myfile;
    myfile.open(vtk_file_name.c_str());

    myfile << "# vtk DataFile Version 2.0" << std::endl;
    myfile << "DENDRO OCTREES" << std::endl;
    myfile << "ASCII" << std::endl;
    myfile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    int dim = 3;

    int unit_points = 1 << dim;  // 2^dim, nodes per element
    int num_cells = nodes.size();  // TODO
    int num_verticies = num_cells * unit_points;

    if (!rank) std::cout << "  writing " << num_verticies << " vertices...\n";

    myfile << "POINTS " << num_verticies << " float" << std::endl;

    if (dim == 3) {
        unsigned int len;
        unsigned int xl, yl, zl;
        int num_data_field = 2;

        int maxD = da.getMaxDepth();
        double hx = 1.0 / ((double)(1 << (maxD-1)));
        for ( da.init<ot::DA_FLAGS::ALL>(), da.init<ot::DA_FLAGS::WRITABLE>(); da.curr() < da.end<ot::DA_FLAGS::ALL>(); da.next<ot::DA_FLAGS::ALL>()) {
          unsigned int i = da.curr();
          Point pt = da.getCurrentOffset();

          unsigned int lev = da.getLevel(i);
          unsigned int w = ((1 << (maxD - lev)));
          double coords[8][3] = {
            {pt.x(), pt.y(), pt.z()},
            {pt.x()+w, pt.y(), pt.z()},
            {pt.x()+w, pt.y()+w, pt.z()},
            {pt.x(), pt.y()+w, pt.z()},

            {pt.x(), pt.y(), pt.z()+1},
            {pt.x()+w, pt.y(), pt.z()+1},
            {pt.x()+w, pt.y()+w, pt.z()+1},
            {pt.x(), pt.y()+w, pt.z()+1},
          };

          for (int i = 0; i < 8; i++) {
            double p[3] = { coords[i][0] * hx, coords[i][1] * hx, coords[i][2] * hx };
            myfile << p[0] << " " << p[1] << " " << p[2] << std::endl;
            std::cout << "elem " << i << ": " << coords[0] << ", " << coords[1] << ", " << coords[2] << "\n";
          }
        }

        int num_cells_elements = num_cells * unit_points + num_cells;
        myfile << "CELLS " << num_cells << " " << num_cells_elements << std::endl;

        for (int i = 0; i < num_cells; i++) {
            myfile << unit_points << " ";
            for (int j = 0; j < unit_points; j++) {
                myfile << (i * unit_points + j) << " ";
            }
            myfile << std::endl;
        }

        myfile << "CELL_TYPES " << num_cells << std::endl;
        for (int i = 0; i < num_cells; i++) {
            myfile << VTK_HEXAHEDRON << std::endl;
        }

        myfile<<std::endl;

        myfile<< "POINT_DATA "<<num_verticies<<std::endl;


        myfile <<"SCALARS U FLOAT"<<std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;
        for(unsigned int j=0;j<num_verticies;j++) {
            myfile << u[(j) * (NUM_FUNC_VALUES) + U] << std::endl;
        }
        myfile<<std::endl;


        myfile<<"CELL_DATA "<<num_cells<<std::endl;
        //myfile<<"POINT_DATA "<<(num_cells*unit_points)<<std::endl;

        //myfile << "FIELD OCTREE_DATA " << num_data_field << std::endl;

        //myfile << "cell_level 1 " << num_cells << " int" << std::endl;
        myfile << "SCALARS cell_level FLOAT"<< std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;

        for (int i = 0; i < nodes.size(); i++)
            myfile << nodes[i].getLevel() <<std::endl;

        myfile << std::endl;
        //myfile << std::endl;

        //myfile << "mpi_rank 1 " << num_cells << " int" << std::endl;
        myfile << "SCALARS mpi_rank FLOAT"<< std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;
        for (int i = 0; i < nodes.size(); i++)
            myfile <<rank << std::endl;

        myfile << std::endl;
    }

    myfile.close();
}*/
