cmake_minimum_required(VERSION 2.8)
project(dendro)

set (PETSC_DIR /home/lofquist/LIBS/petsc-3.5.4)
set (PETSC_ARCH arch-linux2-cxx-debug)
set (DENDRO_DIR /home/lofquist/LIBS/Dendro4)
set (TALYFEM_DIR /home/lofquist/talyfemlib)
set (LIBCONFIG_DIR /home/lofquist/LIBS/libconfig-1.5/install)

# version number
set (Dendro_VERSION_MAJOR 4)
set (Dendro_VERSION_MINOR 0)

#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${dendro_SOURCE_DIR}/lib)
#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${dendro_SOURCE_DIR}/lib)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${dendro_SOURCE_DIR}/bin)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

#set (PAPI_LIB /home/milinda/Programs/papi-5.4.3/src)
#set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -l/home/milinda/Research/Utah/DENDRO/Libs/petsc/lib/libpapi.a")

message (STATUS ${OpenMP_CXX_FLAGS} ${OpenMP_EXE_LINKER_FLAGS})

if(OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(MPI_COMPILE_FLAGS)
  set(COMPILE_FLAGS "${COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set(LINK_FLAGS "${LINK_FLAGS} ${MPI_LINK_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")
set(LINK_FLAGS "${LINK_FLAGS} -rdynamic")

# Options

##------
include_directories(${PROJECT_BINARY_DIR}
                    include
                    ${DENDRO_DIR}/build
                    ${DENDRO_DIR}/include
                    ${DENDRO_DIR}/include/fem
                    ${DENDRO_DIR}/include/oda
                    ${DENDRO_DIR}/include/omg
                    ${DENDRO_DIR}/include/test
                    ${DENDRO_DIR}/examples/include
                    $ENV{PETSC_DIR}/include
                    $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include
                    ${TALYFEM_DIR}
                    ${LIBCONFIG_DIR}/include
                    ${MPI_INCLUDE_PATH}
                    )

link_directories($ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
                 $ENV{PAPI_LIB}
                 ${DENDRO_DIR}/build 
                 ${TALYFEM_DIR}
                 ${LIBCONFIG_DIR}/lib
        )

# Sources

## executables

##---------------------------------------------------------------------------------------
##  Helper APPS
##---------------------------------------------------------------------------------------
add_executable(dendroHeat src/main.cpp src/timeStepper.cpp src/TecplotIO_ascii.cpp src/VecIO.cpp src/DendroIO_vtk.cpp src/DendroIO_plt.cpp src/interp.cpp)
target_link_libraries(dendroHeat dendro dendroDA petsc talyfem parmetis config++ ${MPI_LIBRARIES} m)
