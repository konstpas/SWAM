module bc_module

  ! -----------------------------------------------------------------
  ! This module is used to specify the boundary conditions on the
  ! entire heat condution domain.
  ! -----------------------------------------------------------------
  
  use amrex_base_module

  implicit none

  ! Homogeneous Neumann boundary condition (hoextrapcc stands implies that the
  ! cell-center of the ghost cells are filled with a copy of the closest point
  ! inside the domain)
  ! In the following the second argument of lo_bc and of hi_bc represents the
  ! number of components of the multifab to which the boundary conditions
  ! should be applied
  integer, save :: lo_bc(amrex_spacedim,1) = amrex_bc_hoextrapcc
  integer, save :: hi_bc(amrex_spacedim,1) = amrex_bc_hoextrapcc
  
  
  
end module bc_module
