AMREX_HOME ?= ../../../amrex

DEBUG     = TRUE
DEBUG     = FALSE

USE_MPI   = FALSE
USE_OMP   = FALSE
USE_MPI3  = FALSE

MEM_PROFILE  = FALSE
TINY_PROFILE = FALSE

COMP      = gnu

DIM       = 3

USE_F_INTERFACES = TRUE

Bpack   := ./Make.package 
Blocs   := . 

include $(Bpack)
include $(AMREX_HOME)/Tools/GNUMake/Make.defs

Pdirs 	:= Base Boundary AmrCore Particle F_Interfaces/Base F_Interfaces/AmrCore F_Interfaces/Particle
Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

$(info $(Ppack))

include $(Ppack)


include $(AMREX_HOME)/Tools/GNUMake/Make.rules



