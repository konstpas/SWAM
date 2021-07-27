module amr_data_module
  
  ! -----------------------------------------------------------------
  ! This module is used to declare, allocate and destroy the global
  ! variables employed in the AMReX calculations
  ! -----------------------------------------------------------------
  
  use iso_c_binding
  use amrex_amr_module
  
  implicit none

  private

  ! ------------------------------------------------------------------
  ! Public variables
  ! ------------------------------------------------------------------
  ! Time step (one for each level)
  public :: dt
  ! Flux registers
  public :: flux_reg
  ! Indexes used to label material and background
  public :: idomain_new
  public :: idomain_old
  ! Boundary conditions of the entire simulation box
  public :: lo_bc
  public :: hi_bc
  ! Bottom of the melt pool (2D)
  public :: melt_pos
  ! Melt velocity (2D)
  public :: melt_vel 
  ! Enthalpy
  public :: phi_new
  public :: phi_old
  ! Grid parameters for the 2D grid used in shallow water
  public :: surf_dx
  public :: surf_ind
  public :: surf_pos
  public :: surf_xlo
  ! Temperature
  public :: temp
  ! Time
  public :: t_new
  public :: t_old
  ! Parameters used by AMReX during subcycling
  public :: stepno
  public :: nsubsteps

  ! ------------------------------------------------------------------
  ! Public subroutines
  ! ------------------------------------------------------------------
  public :: amr_data_init
  public :: amr_data_finalize

  
  ! ------------------------------------------------------------------
  ! Declare public variables
  ! ------------------------------------------------------------------
  integer, allocatable, save :: lo_bc(:,:)
  integer, allocatable, save :: hi_bc(:,:)
  integer, save  :: surf_ind(2,2)
  integer, allocatable, save :: stepno(:)
  integer, allocatable, save :: nsubsteps(:)
  real(amrex_real), allocatable, save :: dt(:)
  real(amrex_real), allocatable, save :: melt_pos(:,:)
  real(amrex_real), allocatable, save :: melt_vel(:,:,:)
  real(amrex_real), allocatable, save :: t_new(:)
  real(amrex_real), allocatable, save :: t_old(:)
  real(amrex_real), allocatable, save :: surf_pos(:,:)
  real(amrex_real), save :: surf_xlo(2)
  real(amrex_real), save :: surf_dx(2)
  type(amrex_fluxregister), allocatable, save :: flux_reg(:)
  type(amrex_imultifab), allocatable, save :: idomain_new(:)
  type(amrex_imultifab), allocatable, save :: idomain_old(:)
  type(amrex_multifab), allocatable, save :: phi_new(:)
  type(amrex_multifab), allocatable, save :: phi_old(:)
  type(amrex_multifab), allocatable, save :: temp(:)

  
contains

  ! ------------------------------------------------------------------
  ! Subroutine to allocate and initialize the global variables.
  ! Note: this subroutine initializes all the variables that are not
  ! Multifabs or flux registers. Multifabs and flux registers are 
  ! initialized whenever a new level is created and this is done with
  ! some specific subroutines given in the regrid module.
  ! ------------------------------------------------------------------ 
  subroutine amr_data_init()
    
    use read_input_module, only : surf_pos_init
    
    integer ::  lev
    integer ::  lo_x
    integer ::  hi_x
    integer ::  lo_z
    integer ::  hi_z

    lo_x = amrex_geom(amrex_max_level)%domain%lo(1)
    hi_x = amrex_geom(amrex_max_level)%domain%hi(1)
    lo_z = amrex_geom(amrex_max_level)%domain%lo(3)
    hi_z = amrex_geom(amrex_max_level)%domain%hi(3)
    
    ! Allocate
    allocate(dt(0:amrex_max_level))
    allocate(flux_reg(0:amrex_max_level))
    allocate(idomain_new(0:amrex_max_level))
    allocate(idomain_old(0:amrex_max_level))
    allocate(lo_bc(amrex_spacedim,1)) ! The second argument is the number of components
    allocate(hi_bc(amrex_spacedim,1))
    allocate(melt_pos(lo_x:hi_x, lo_z:hi_z))
    allocate(melt_vel(lo_x:hi_x+1,lo_z:hi_z+1, 1:amrex_spacedim-1))
    allocate(phi_new(0:amrex_max_level))
    allocate(phi_old(0:amrex_max_level))
    allocate(surf_pos(lo_x:hi_x,lo_z:hi_z))
    allocate(temp(0:amrex_max_level))
    allocate(t_new(0:amrex_max_level))
    allocate(t_old(0:amrex_max_level))
    allocate(stepno(0:amrex_max_level))
    allocate(nsubsteps(0:amrex_max_level))

    ! Initialize
    dt = 1.0_amrex_real
    ! Homogeneous Neumann boundary condition (foextrap implies that the ghost
    ! cells are filled with a copy of the closest point inside the domain)
    lo_bc = amrex_bc_foextrap
    hi_bc = amrex_bc_foextrap
    ! It is assumed that there is no melting pool at the beginning of the
    ! simulation (melt_pos = surf_pos)
    melt_pos = surf_pos_init
    melt_vel = 0.0_amrex_real
    surf_dx(1) = amrex_geom(amrex_max_level)%dx(1)
    surf_dx(2) = amrex_geom(amrex_max_level)%dx(3)
    surf_ind(1,1) = lo_x
    surf_ind(1,2) = hi_x
    surf_ind(2,1) = lo_z
    surf_ind(2,2) = hi_z 
    surf_pos = surf_pos_init 
    surf_xlo(1) = amrex_problo(1) 
    surf_xlo(2) = amrex_problo(3)
    t_new = 0.0_amrex_real
    t_old = -1.0_amrex_real
    stepno = 0 
    nsubsteps(0) = 1
    do lev = 1, amrex_max_level
       nsubsteps(lev) = amrex_ref_ratio(lev-1)
    end do
    
  end subroutine amr_data_init

  ! ------------------------------------------------------------------
  ! Subroutine to free the memory
  ! ------------------------------------------------------------------
  subroutine amr_data_finalize()

    use read_input_module, only : deallocate_input
    
    integer :: lev
    
    deallocate(dt)
    deallocate(lo_bc)
    deallocate(hi_bc)
    deallocate(melt_pos)
    deallocate(melt_vel)
    deallocate(surf_pos)
    deallocate(t_new)
    deallocate(t_old)
    deallocate(stepno)
    deallocate(nsubsteps)
    
    do lev = 0, amrex_max_level

       call amrex_imultifab_destroy(idomain_new(lev))
       call amrex_imultifab_destroy(idomain_old(lev))
       call amrex_multifab_destroy(phi_new(lev))
       call amrex_multifab_destroy(phi_old(lev))
       call amrex_multifab_destroy(temp(lev))
      
    end do
    
    do lev = 1, amrex_max_level
       call amrex_fluxregister_destroy(flux_reg(lev))
    end do

    call deallocate_input
    
  end subroutine amr_data_finalize
  
end module amr_data_module
