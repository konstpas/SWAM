# Compiler
CC = gfortran

# Amrex
AMREX_LIBRARY_HOME ?= /home/federico/amrex_local
LIBDIR := $(AMREX_LIBRARY_HOME)/lib
INCDIR := $(AMREX_LIBRARY_HOME)/include
COMPILE_CPP_FLAGS ?= $(shell awk '/Cflags:/ {$$1=$$2=""; print $$0}' $(LIBDIR)/pkgconfig/amrex.pc)
COMPILE_LIB_FLAGS ?= $(shell awk '/Libs:/ {$$1=$$2=""; print $$0}' $(LIBDIR)/pkgconfig/amrex.pc)
CFLAGS := -I$(INCDIR) $(COMPILE_CPP_FLAGS)
LFLAGS := -L$(LIBDIR)

# Files                                                                                              
SRC = amr_data_mod.F90
OBJS = $(patsubst %.F90, %.o, $(SRC))

# Name of the executable                                                         
EXECUTABLE = swam

# Compile                                                                                            
all: $(EXECUTABLE)

%.o: %.F90
	$(CC) $(CFLAGS) -c $<

# Link                                                                                               
$(EXECUTABLE): $(OBJS)
	$(CC) $^ -o $@ $(LFLAGS)
