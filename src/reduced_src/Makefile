# Compiler
COMP_NAME = gfortran

# Linker
LINK_NAME = g++

# Amrex library
AMREX_LIBRARY_HOME ?= /home/x_fedlu/amrex
AMREX_LIBDIR := $(AMREX_LIBRARY_HOME)/lib
AMREX_INCDIR := $(AMREX_LIBRARY_HOME)/include
CFLAGS := -I$(AMREX_INCDIR) 
LFLAGS := -lgfortran -L$(AMREX_LIBDIR) -lamrex

# FLAP library
FLAP_LIBRARY_HOME ?= /home/x_fedlu/FLAP
FLAP_LIBDIR = $(FLAP_LIBRARY_HOME)/static
FLAP_INCDIR := $(FLAP_LIBRARY_HOME)/static/mod
CFLAGS += -I$(FLAP_INCDIR)
LFLAGS += -L$(FLAP_LIBDIR) -lflap

# Files                                                                                             
SRCF90 = $(wildcard *.F90)
OBJS += $(patsubst %.F90, %.o, $(SRCF90))

# Name of the executable                                                                             
EXECUTABLE = swam

# Compile                                                                                 
all: $(EXECUTABLE)

# Dependencies
include .depend

# Compiling rules
%.o: %.F90
	$(COMP_NAME) $(CFLAGS) -c $<

# Link                                                                                               
$(EXECUTABLE): $(OBJS)
	$(LINK_NAME) -fopenmp $^ -o $@ $(LFLAGS)

# Generate dependencies
depend .depend:
	makedepf90 *.F90 >.depend

# Clean
clean:
	rm *.mod *.o .depend
