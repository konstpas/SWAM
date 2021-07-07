# Compiler
COMP_NAME = gfortran

# Linker
LINK_NAME = g++

# Amrex
AMREX_LIBRARY_HOME ?= /home/federico/amrex_local
LIBDIR := $(AMREX_LIBRARY_HOME)/lib
INCDIR := $(AMREX_LIBRARY_HOME)/include
CFLAGS := -I$(INCDIR) 
LFLAGS := -L$(LIBDIR) -lamrex -lgfortran

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
	$(COMP_NAME) $(CFLAGS) -c $<

# Link                                                                                               
$(EXECUTABLE): $(OBJS)
	$(LINK_NAME) -fopenmp $^ -o $@ $(LFLAGS)

# Generate dependencies
depend .depend:
	makedepf90 *.f90 >.depend

# Clean
clean:
	rm *.mod *.o .depend
