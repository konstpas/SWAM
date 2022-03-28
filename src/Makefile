# Compiler
# COMP_NAME = mpif90
COMP_NAME = gfortran

# Compiler options
COMP_OPT = -O2 -Wall

# Linker
LINK_NAME = g++

# Amrex
# AMREX_LIBRARY_HOME ?= /home/konstantinos/amrex_local_mpi
AMREX_LIBRARY_HOME ?= /home/konstantinos/amrex_local
LIBDIR := $(AMREX_LIBRARY_HOME)/lib
INCDIR := $(AMREX_LIBRARY_HOME)/include
COMP_FLAGS := -I$(INCDIR) 
LINK_FLAGS := -L$(LIBDIR) -lamrex -lgfortran
# MPI_COMPILE_FLAGS = $(shell $(COMP_NAME) --showme:compile)
# MPI_LINK_FLAGS = $(shell $(COMP_NAME) --showme:link)

# Files                                                                                              
SRCf90 = $(wildcard *.f90)
SRCF90 = $(wildcard *.F90)
OBJS = $(patsubst %.f90, %.o, $(SRCf90))
OBJS += $(patsubst %.F90, %.o, $(SRCF90))

# Name of the executable                                                                             
EXECUTABLE = swam

# Compile                                                                                 
all: $(EXECUTABLE)

# Dependencies
include .depend

# Compiling rules
%.o: %.f90
	$(COMP_NAME) $(COMP_OPT) -fopenmp $(COMP_FLAGS) $(MPI_COMPILE_FLAGS) -c $<

# Link                                                                                               
$(EXECUTABLE): $(OBJS)
	$(LINK_NAME) -fopenmp $^ -o $@ $(LINK_FLAGS) $(MPI_LINK_FLAGS)

# Generate dependencies
depend .depend:
	makedepf90 *.f90 >.depend

# Clean
clean:
	rm *.mod *.o .depend
