#include "DendroIO.h"
#include "interp.h"

extern double gSize[3];

void octree2VTK(ot::DA* da, Vec vec, int dof, std::string file_prefix) {
  int rank, size;
  char fname[256];

  
	MPI_Comm_rank(da->getComm(), &rank);
	MPI_Comm_size(da->getComm(), &size);

  std::stringstream ss;
  sprintf(fname, "%s_%05d.vtk", file_prefix.c_str(), rank);

  //if ( !rank ) std::cout << "Writing to VTK file: " << fname << std::endl;

  std::ofstream out;
  out.open( fname );
  
  
  out << "# vtk DataFile Version 2.0" << std::endl;
  out << "DENDRO OCTREES" << std::endl;
  out << "ASCII" << std::endl;
  out << "DATASET UNSTRUCTURED_GRID" << std::endl;

  int dim = 3;

  int unit_points = 1 << dim;
  int num_vertices = da->getElementSize() * (unit_points);
  int num_cells = da->getElementSize();

  out << "POINTS " << num_vertices << " float" << std::endl;

  { // dim = 3
    
    unsigned int len; //!  ??
    unsigned int xl, yl, zl;  //! ??

    int num_data_field = 2; // rank and data 
    int num_cells_elements = num_cells * unit_points + num_cells;

    PetscScalar *_vec=NULL; 

    da->vecGetBuffer(vec,   _vec, false, false, true,  dof);

    da->ReadFromGhostsBegin<PetscScalar>(_vec, dof);
    da->ReadFromGhostsEnd<PetscScalar>(_vec);

    unsigned int maxD = da->getMaxDepth();
    unsigned int lev;
    double hx, hy, hz;
    Point pt;

    double xFac = gSize[0]/((double)(1<<(maxD-1)));
    double yFac = gSize[1]/((double)(1<<(maxD-1)));
    double zFac = gSize[2]/((double)(1<<(maxD-1)));
    double xx, yy, zz;
    unsigned int idx[8];

    for ( da->init<ot::DA_FLAGS::WRITABLE>(); da->curr() < da->end<ot::DA_FLAGS::WRITABLE>(); da->next<ot::DA_FLAGS::WRITABLE>() ) { 
      // set the value
      lev = da->getLevel(da->curr());
      hx = xFac*(1<<(maxD - lev));
      hy = yFac*(1<<(maxD - lev));
      hz = zFac*(1<<(maxD - lev));

      pt = da->getCurrentOffset();

      xx = pt.x()*xFac; yy = pt.y()*yFac; zz = pt.z()*zFac;
			
      out << pt.x()*xFac << " " <<  pt.y()*yFac << " " << pt.z()*zFac << std::endl;
      out << pt.x()*xFac + hx << " " <<  pt.y()*yFac << " " << pt.z()*zFac << std::endl;
      out << pt.x()*xFac + hx << " " <<  pt.y()*yFac + hy << " " << pt.z()*zFac << std::endl;
      out << pt.x()*xFac << " " <<  pt.y()*yFac + hy << " " << pt.z()*zFac << std::endl;

      out << pt.x()*xFac << " " <<  pt.y()*yFac << " " << pt.z()*zFac + hz<< std::endl;
      out << pt.x()*xFac + hx << " " <<  pt.y()*yFac << " " << pt.z()*zFac + hz << std::endl;
      out << pt.x()*xFac + hx << " " <<  pt.y()*yFac + hy << " " << pt.z()*zFac + hz << std::endl;
      out << pt.x()*xFac << " " <<  pt.y()*yFac + hy << " " << pt.z()*zFac + hz << std::endl;
    }

    
    
    out << "CELLS " << da->getElementSize() << " " << num_cells_elements << std::endl;

    for (int i = 0; i < num_cells; i++) {
      out << unit_points << " ";
      for (int j = 0; j < unit_points; j++) {
        out << (i * unit_points + j) << " ";
      }
      out << std::endl;
    }

    out << "CELL_TYPES " << num_cells << std::endl;
    for (int i = 0; i < num_cells; i++) {
      out << VTK_HEXAHEDRON << std::endl;
    }

    //myfile<<"CELL_DATA "<<num_cells<<std::endl;
    
    out << std::endl;
    out << "POINT_DATA " << num_vertices  << std::endl;

    // this is inefficient - we will interpolate all data dof times
    PetscScalar* local = new PetscScalar[8*dof];
    for (int d = 0; d < dof; d++) {
      // name fields data0, data1, etc.
      out << "SCALARS data" << d << " float 1" << std::endl;
      out << "LOOKUP_TABLE default" << std::endl;
      
      for ( da->init<ot::DA_FLAGS::WRITABLE>(); da->curr() < da->end<ot::DA_FLAGS::WRITABLE>(); da->next<ot::DA_FLAGS::WRITABLE>() ) { 
        da->getNodeIndices(idx);
        interp_global_to_local(_vec, local, da, dof);

          out << local[0*dof+d] << " ";
          out << local[1*dof+d] << " ";
          out << local[3*dof+d] << " ";
          out << local[2*dof+d] << " ";
          out << local[4*dof+d] << " ";
          out << local[5*dof+d] << " ";
          out << local[7*dof+d] << " ";
          out << local[6*dof+d] << " ";
      }
    }
    delete[] local;

    out << std::endl;

/*

    out << "FIELD OCTREE_DATA " << num_data_field << std::endl;

    out << "cell_level 1 " << num_cells << " int" << std::endl;

    for ( da->init<ot::DA_FLAGS::ALL>(); da->curr() < da->end<ot::DA_FLAGS::ALL>(); da->next<ot::DA_FLAGS::ALL>() ) { 
      int l = da->getLevel(da->curr()); 
      out << l << " ";
    }

    out << std::endl;

    out << "mpi_rank 1 " << num_cells << " int" << std::endl;
    for (int i = 0; i < num_cells; i++)
      out << rank << " ";

    out << std::endl;
*/
    da->vecRestoreBuffer(vec,  _vec, false, false, true,  dof);
  }

  out.close();
}







void octree2VTK_ghost(ot::DA* da, Vec vec, int dof, std::string file_prefix) {
  int rank, size;
  char fname[256];

  
	MPI_Comm_rank(da->getComm(), &rank);
	MPI_Comm_size(da->getComm(), &size);

  std::stringstream ss;
  sprintf(fname, "%s_%05d.vtk", file_prefix.c_str(), rank);

  //if ( !rank ) std::cout << "Writing to VTK file: " << fname << std::endl;

  std::ofstream out;
  out.open( fname );
  
  
  out << "# vtk DataFile Version 2.0" << std::endl;
  out << "DENDRO OCTREES" << std::endl;
  out << "ASCII" << std::endl;
  out << "DATASET UNSTRUCTURED_GRID" << std::endl;

  int dim = 3;

  int unit_points = 1 << dim;
  int num_vertices = da->getElementSize() * (unit_points);
  int num_cells = da->getElementSize();

  out << "POINTS " << num_vertices << " float" << std::endl;

  { // dim = 3
    
    unsigned int len; //!  ??
    unsigned int xl, yl, zl;  //! ??

    int num_data_field = 2; // rank and data 
    int num_cells_elements = num_cells * unit_points + num_cells;

    PetscScalar *_vec=NULL; 

    da->vecGetBuffer(vec,   _vec, false, false, true,  dof);

    da->ReadFromGhostsBegin<PetscScalar>(_vec, dof);
    da->ReadFromGhostsEnd<PetscScalar>(_vec);

    unsigned int maxD = da->getMaxDepth();
    unsigned int lev;
    double hx, hy, hz;
    Point pt;

    double xFac = gSize[0]/((double)(1<<(maxD-1)));
    double yFac = gSize[1]/((double)(1<<(maxD-1)));
    double zFac = gSize[2]/((double)(1<<(maxD-1)));
    double xx, yy, zz;
    unsigned int idx[8];

    for ( da->init<ot::DA_FLAGS::ALL>(); da->curr() < da->end<ot::DA_FLAGS::ALL>(); da->next<ot::DA_FLAGS::ALL>() ) { 
      // set the value
      lev = da->getLevel(da->curr());
      hx = xFac*(1<<(maxD - lev));
      hy = yFac*(1<<(maxD - lev));
      hz = zFac*(1<<(maxD - lev));

      pt = da->getCurrentOffset();

      xx = pt.x()*xFac; yy = pt.y()*yFac; zz = pt.z()*zFac;
			
      out << pt.x()*xFac << " " <<  pt.y()*yFac << " " << pt.z()*zFac << std::endl;
      out << pt.x()*xFac + hx << " " <<  pt.y()*yFac << " " << pt.z()*zFac << std::endl;
      out << pt.x()*xFac + hx << " " <<  pt.y()*yFac + hy << " " << pt.z()*zFac << std::endl;
      out << pt.x()*xFac << " " <<  pt.y()*yFac + hy << " " << pt.z()*zFac << std::endl;

      out << pt.x()*xFac << " " <<  pt.y()*yFac << " " << pt.z()*zFac + hz<< std::endl;
      out << pt.x()*xFac + hx << " " <<  pt.y()*yFac << " " << pt.z()*zFac + hz << std::endl;
      out << pt.x()*xFac + hx << " " <<  pt.y()*yFac + hy << " " << pt.z()*zFac + hz << std::endl;
      out << pt.x()*xFac << " " <<  pt.y()*yFac + hy << " " << pt.z()*zFac + hz << std::endl;
    }

    
    
    out << "CELLS " << da->getElementSize() << " " << num_cells_elements << std::endl;

    for (int i = 0; i < num_cells; i++) {
      out << unit_points << " ";
      for (int j = 0; j < unit_points; j++) {
        out << (i * unit_points + j) << " ";
      }
      out << std::endl;
    }

    out << "CELL_TYPES " << num_cells << std::endl;
    for (int i = 0; i < num_cells; i++) {
      out << VTK_HEXAHEDRON << std::endl;
    }

    //myfile<<"CELL_DATA "<<num_cells<<std::endl;
    
    out << std::endl;
    out << "POINT_DATA " << num_vertices  << std::endl;

    // this is inefficient - we will interpolate all data dof times
    PetscScalar* local = new PetscScalar[8*dof];
    for (int d = 0; d < dof; d++) {
      // name fields data0, data1, etc.
      out << "SCALARS data" << d << " float 1" << std::endl;
      out << "LOOKUP_TABLE default" << std::endl;
      
      for ( da->init<ot::DA_FLAGS::ALL>(); da->curr() < da->end<ot::DA_FLAGS::ALL>(); da->next<ot::DA_FLAGS::ALL>() ) { 
        da->getNodeIndices(idx);
        interp_global_to_local(_vec, local, da, dof);

          out << local[0*dof+d] << " ";
          out << local[1*dof+d] << " ";
          out << local[3*dof+d] << " ";
          out << local[2*dof+d] << " ";
          out << local[4*dof+d] << " ";
          out << local[5*dof+d] << " ";
          out << local[7*dof+d] << " ";
          out << local[6*dof+d] << " ";
      }
    }
    delete[] local;

    out << std::endl;

/*

    out << "FIELD OCTREE_DATA " << num_data_field << std::endl;

    out << "cell_level 1 " << num_cells << " int" << std::endl;

    for ( da->init<ot::DA_FLAGS::ALL>(); da->curr() < da->end<ot::DA_FLAGS::ALL>(); da->next<ot::DA_FLAGS::ALL>() ) { 
      int l = da->getLevel(da->curr()); 
      out << l << " ";
    }

    out << std::endl;

    out << "mpi_rank 1 " << num_cells << " int" << std::endl;
    for (int i = 0; i < num_cells; i++)
      out << rank << " ";

    out << std::endl;
*/
    da->vecRestoreBuffer(vec,  _vec, false, false, true,  dof);
  }

  out.close();
}
