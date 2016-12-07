# Set this to the path to your Dendro4 directory.
DENDRO_DIR = $(HOME)/Dendro4
# PETSC_DIR = $(HOME)/LIBS/petsc-3.5.4

include ${PETSC_DIR}/conf/variables

# List all the files in src/ that you want to compile.
SOURCES = main.cpp timeStepper.cpp

# Name of the executable you want to build.
EXECUTABLE = dendroHeat

CXX = g++
CXX_FLAGS += -std=c++11 -fopenmp

INCLUDES = -I./include -I$(DENDRO_DIR)/include -I$(DENDRO_DIR)/include/oda -I$(DENDRO_DIR)/include/fem -I$(DENDRO_DIR)/include/omg -I$(DENDRO_DIR)/build $(PETSC_CC_INCLUDES)
LIBS = $(PETSC_LIB) -L$(DENDRO_DIR)/build -ldendroDA -ldendro

# Generate the list of object files by replacing suffixes in SOURCES (e.g. main.cpp -> main.o).
OBJECTS = $(SOURCES:.cpp=.o)

# This tells the computer where to look for files.
VPATH = src

# Our executable depends on all the object files.
$(EXECUTABLE): $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LIBS)

# When running 'make' with no target specified, make defaults to the 'all' target.
# By default, clean and build our program.
all:
	$(MAKE) clean
	$(MAKE) $(EXECUTABLE)

# Remove all '.o' files and the executable, forcing make to rebuild everything from scratch next time.
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

### Implicit rules.
# We use this to say how to convert .cpp files to .o files, without having to name .cpp and .o files by hand.

.SUFFIXES: .cpp

.cpp.o:
	$(CXX) -c $(CXX_FLAGS) $(INCLUDES) $<
