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
  ! Public subroutines
  ! ------------------------------------------------------------------
  public :: amr_data_init
  public :: amr_data_finalize
  
  
  ! ------------------------------------------------------------------
  ! Public variables
  ! ------------------------------------------------------------------
  ! Time step 
  public :: dt
  ! Flux registers
  public :: flux_reg
  ! Enthalpy in the heat transfer domain
  public :: phi_new
  public :: phi_old
  ! Temperature in the heat transfer domain
  public :: temp
  ! Indexes used to label material and background
  public :: idomain
  ! Boundary conditions of the heat transfer domain
  public :: lo_bc
  public :: hi_bc
  ! Bottom of the melt pool
  public :: melt_pos
  ! Top of the melt pool
  public :: melt_top
  ! Melt velocity
  public :: melt_vel
  ! Maximum melt velocity (absolute value)
  public :: max_melt_vel
  ! Tabulated heat flux from plasma
  public :: plasma_flux_time_mesh
  public :: plasma_flux_surf_mesh
  public :: plasma_flux_side_time_mesh
  public :: plasma_flux_side_surf_mesh
  public :: plasma_flux_table
  public :: plasma_flux_side_table
  ! Resolution of the shallow water grid
  public :: surf_dx
  ! Index of the shallow water grid
  public :: surf_ind
  ! Lower corner of the shallow water domain
  public :: surf_xlo
  ! Position of the free surface
  public :: surf_pos 
  ! Position of the free surface in the heat transfer domain
  ! public :: surf_pos_grid
  ! Temperature on the free surface
  public :: surf_temperature
  ! Enthalpy on the free surface
  public :: surf_enthalpy
  ! Evaporation flux on the free surface
  public :: surf_evap_flux
  ! Thermionic current on the free surface
  public :: surf_current
  ! Time
  public :: t_new
  public :: t_old
  ! Parameters used by amrex during subcycling
  public :: stepno
  public :: nsubsteps
  ! Variable used to check when to regrid
  public :: last_regrid_step

  ! ------------------------------------------------------------------
  ! Declare public variables
  ! ------------------------------------------------------------------
  integer, allocatable, save :: hi_bc(:,:)
  integer, allocatable, save :: lo_bc(:,:)
  integer, allocatable, save :: last_regrid_step(:)
  integer, allocatable, save :: nsubsteps(:)
  integer, allocatable, save :: stepno(:)
  integer, save  :: surf_ind(1,2)
  real(amrex_real), allocatable, save :: dt(:)
  real(amrex_real), allocatable, save :: melt_pos(:)
  real(amrex_real), allocatable, save :: melt_top(:)
  real(amrex_real), save :: max_melt_vel
  real(amrex_real), allocatable, save :: melt_vel(:,:)
  real(amrex_real), allocatable, save :: plasma_flux_time_mesh(:)
  real(amrex_real), allocatable, save :: plasma_flux_surf_mesh(:)
  real(amrex_real), allocatable, save :: plasma_flux_table(:,:)
  real(amrex_real), allocatable, save :: plasma_flux_side_time_mesh(:)
  real(amrex_real), allocatable, save :: plasma_flux_side_surf_mesh(:)
  real(amrex_real), allocatable, save :: plasma_flux_side_table(:,:)
  real(amrex_real), allocatable, save :: surf_current(:)
  real(amrex_real), allocatable, save :: surf_enthalpy(:)
  real(amrex_real), allocatable, save :: surf_evap_flux(:)
  real(amrex_real), allocatable, save :: surf_pos(:)
  real(amrex_real), allocatable, save :: surf_pos_grid(:)
  real(amrex_real), allocatable, save :: t_new(:)
  real(amrex_real), allocatable, save :: t_old(:)
  real(amrex_real), save :: surf_dx(1)
  real(amrex_real), allocatable, save :: surf_temperature(:)
  real(amrex_real), save :: surf_xlo(1)
  type(amrex_fluxregister), allocatable, save :: flux_reg(:)
  type(amrex_multifab), allocatable, save :: idomain(:)
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
    
    use read_input_module, only : sw_surf_pos_init
                                 
    integer ::  lev
    integer ::  lo_x
    integer ::  hi_x

    lo_x = amrex_geom(amrex_max_level)%domain%lo(1)
    hi_x = amrex_geom(amrex_max_level)%domain%hi(1)
    
    ! Allocate
    allocate(dt(0:amrex_max_level))
    allocate(flux_reg(0:amrex_max_level))
    allocate(idomain(0:amrex_max_level))
    allocate(lo_bc(amrex_spacedim,1)) ! The second argument is the number of components
    allocate(hi_bc(amrex_spacedim,1))
    allocate(last_regrid_step(0:amrex_max_level))
    allocate(melt_pos(lo_x:hi_x))
    allocate(melt_top(lo_x:hi_x))
    allocate(melt_vel(lo_x:hi_x+1, 1:amrex_spacedim-1))
    allocate(phi_new(0:amrex_max_level))
    allocate(phi_old(0:amrex_max_level))
    allocate(surf_current(lo_x:hi_x))
    allocate(surf_enthalpy(lo_x:hi_x))
    allocate(surf_evap_flux(lo_x:hi_x))
    allocate(surf_pos(lo_x:hi_x))
    allocate(surf_pos_grid(lo_x:hi_x))
    allocate(surf_temperature(lo_x:hi_x))
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
    last_regrid_step = 0
    ! It is assumed that there is no melting pool at the beginning of the
    ! simulation (melt_pos = surf_pos)
    melt_pos = sw_surf_pos_init
    melt_top = sw_surf_pos_init
    max_melt_vel = 0.0_amrex_real
    melt_vel = 0.0_amrex_real
    surf_dx(1) = amrex_geom(amrex_max_level)%dx(1)
    surf_ind(1,1) = lo_x
    surf_ind(1,2) = hi_x
    surf_current = 0.0_amrex_real
    surf_enthalpy = 0.0_amrex_real
    surf_evap_flux = 0.0_amrex_real
    surf_pos = sw_surf_pos_init
    surf_pos_grid = sw_surf_pos_init
    surf_temperature = 0.0_amrex_real
    surf_xlo(1) = amrex_problo(1)
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

    ! Deallocate public variables
    deallocate(dt)
    deallocate(lo_bc)
    deallocate(hi_bc)
    deallocate(last_regrid_step)
    deallocate(melt_pos)
    deallocate(melt_vel)
    deallocate(surf_current)
    deallocate(surf_enthalpy)
    deallocate(surf_evap_flux)
    deallocate(surf_pos)
    deallocate(surf_pos_grid)
    deallocate(surf_temperature)
    deallocate(t_new)
    deallocate(t_old)
    deallocate(stepno)
    deallocate(nsubsteps)

    ! Deallocate public variables allocated only for certain geometries
    if (allocated(plasma_flux_time_mesh)) deallocate(plasma_flux_time_mesh)
    if (allocated(plasma_flux_surf_mesh)) deallocate(plasma_flux_surf_mesh)
    if (allocated(plasma_flux_table)) deallocate(plasma_flux_table)
    if (allocated(plasma_flux_side_time_mesh)) deallocate(plasma_flux_side_time_mesh)
    if (allocated(plasma_flux_side_surf_mesh)) deallocate(plasma_flux_side_surf_mesh)
    if (allocated(plasma_flux_side_table)) deallocate(plasma_flux_side_table)
    
    ! Destroy multifabs
    do lev = 0, amrex_max_level

       call amrex_multifab_destroy(idomain(lev))
       call amrex_multifab_destroy(phi_new(lev))
       call amrex_multifab_destroy(phi_old(lev))
       call amrex_multifab_destroy(temp(lev))
    
    end do

    ! Destroy flux registers
    do lev = 1, amrex_max_level
       call amrex_fluxregister_destroy(flux_reg(lev))
    end do

    ! Deallocate variables defined in the input module
    call deallocate_input
    
  end subroutine amr_data_finalize
  
end module amr_data_module
