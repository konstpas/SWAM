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
  ! Tabulated plasma heat flux
  public :: plasma_flux_time_mesh
  public :: plasma_flux_surf_mesh
  public :: plasma_flux_table
  ! Tabulated plasma heat flux on second exposed surface
  public :: plasma_flux_side_time_mesh
  public :: plasma_flux_side_surf_mesh
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
  ! Free surface deformation (dh/dt), used for regridding
  public :: surf_deformation
  ! Time
  public :: t_new
  public :: t_old
  ! Parameters used by amrex during subcycling
  public :: stepno
  public :: nsubsteps
  ! Variable used to check when to regrid
  public :: last_regrid_step

  ! Electrostatic potential
  public :: psi
  
  ! ------------------------------------------------------------------
  ! Declare public variables
  ! ------------------------------------------------------------------

  ! Variables related to time
  integer, allocatable, save :: last_regrid_step(:)  
  integer, allocatable, save :: nsubsteps(:)
  integer, allocatable, save :: stepno(:)
  real(amrex_real), allocatable, save :: dt(:)
  real(amrex_real), allocatable, save :: t_new(:)
  real(amrex_real), allocatable, save :: t_old(:)

  ! Variables related to the boundary conditions
  integer, allocatable, save :: hi_bc(:,:)
  integer, allocatable, save :: lo_bc(:,:)

  ! Variables related to the free surface and melt pool
  integer, save  :: surf_ind(1,2)
  real(amrex_real), allocatable, save :: melt_pos(:)
  real(amrex_real), save :: max_melt_vel
  real(amrex_real), allocatable, save :: melt_top(:)
  real(amrex_real), allocatable, save :: melt_vel(:,:)
  real(amrex_real), allocatable, save :: surf_current(:)
  real(amrex_real), allocatable, save :: surf_deformation(:)
  real(amrex_real), allocatable, save :: surf_enthalpy(:)
  real(amrex_real), allocatable, save :: surf_evap_flux(:)
  real(amrex_real), allocatable, save :: surf_pos(:)
  real(amrex_real), allocatable, save :: surf_pos_grid(:)
  real(amrex_real), save :: surf_dx(1)
  real(amrex_real), allocatable, save :: surf_temperature(:)
  real(amrex_real), save :: surf_xlo(1)
  
  ! Variables related to the heat flux read from files
  real(amrex_real), allocatable, save :: plasma_flux_time_mesh(:)
  real(amrex_real), allocatable, save :: plasma_flux_surf_mesh(:)
  real(amrex_real), allocatable, save :: plasma_flux_table(:,:)
  real(amrex_real), allocatable, save :: plasma_flux_side_time_mesh(:)
  real(amrex_real), allocatable, save :: plasma_flux_side_surf_mesh(:)
  real(amrex_real), allocatable, save :: plasma_flux_side_table(:,:)

  ! Multifabs
  type(amrex_multifab), allocatable, save :: idomain(:)
  type(amrex_multifab), allocatable, save :: phi_new(:)
  type(amrex_multifab), allocatable, save :: phi_old(:)
  type(amrex_multifab), allocatable, save :: temp(:)

  ! Flux registers
  type(amrex_fluxregister), allocatable, save :: flux_reg(:)

  ! Electrostatic potential
  type(amrex_multifab), allocatable, save :: psi(:)
contains

  ! ------------------------------------------------------------------
  ! Subroutine to allocate and initialize the global variables.
  ! Note: this subroutine initializes all the variables that are not
  ! Multifabs or flux registers.
  ! ------------------------------------------------------------------ 
  subroutine amr_data_init()
    
    ! Allocate
    ! The variables for the plasma heat flux are allocated, only if
    ! necessary, in the heat_transfer_flux module
    call allocate_time_variables
    call allocate_bc_variables
    call allocate_free_surface_variables
    call allocate_multifabs
    call allocate_flux_registers
    
    ! Initialize
    ! The variables for the plasma heat flux are initialized in the
    ! heat_transfer_flux module. Multifabs and flux registers are 
    ! initialized whenever a new level is created in the regrid
    ! module.
    call init_time_variables
    call init_bc_variables
    call init_free_surface_variables
    
  end subroutine amr_data_init


  subroutine allocate_time_variables()
    
    allocate(dt(0:amrex_max_level))
    allocate(last_regrid_step(0:amrex_max_level))
    allocate(t_new(0:amrex_max_level))
    allocate(t_old(0:amrex_max_level))
    allocate(stepno(0:amrex_max_level))
    allocate(nsubsteps(0:amrex_max_level))
    
  end subroutine allocate_time_variables


  subroutine allocate_bc_variables()

    ! In the following, the second argument is the number of components
    allocate(lo_bc(amrex_spacedim,1))
    allocate(hi_bc(amrex_spacedim,1))
        
  end subroutine allocate_bc_variables


  subroutine allocate_free_surface_variables()

    integer ::  lo_x
    integer ::  hi_x

    lo_x = amrex_geom(amrex_max_level)%domain%lo(1)
    hi_x = amrex_geom(amrex_max_level)%domain%hi(1)

    allocate(melt_pos(lo_x:hi_x))
    allocate(melt_top(lo_x:hi_x))
    allocate(melt_vel(lo_x:hi_x+1, 1:amrex_spacedim-1))
    allocate(surf_current(lo_x:hi_x))
    allocate(surf_deformation(lo_x:hi_x))  
    allocate(surf_enthalpy(lo_x:hi_x))
    allocate(surf_evap_flux(lo_x:hi_x))
    allocate(surf_pos(lo_x:hi_x))
    allocate(surf_pos_grid(lo_x:hi_x))
    allocate(surf_temperature(lo_x:hi_x))
  
  end subroutine allocate_free_surface_variables
  

  subroutine allocate_multifabs()

    allocate(idomain(0:amrex_max_level))
    allocate(phi_new(0:amrex_max_level))
    allocate(phi_old(0:amrex_max_level))
    allocate(temp(0:amrex_max_level))
    allocate(psi(0:amrex_max_level))
    
  end subroutine allocate_multifabs


  subroutine allocate_flux_registers()

    allocate(flux_reg(0:amrex_max_level))
    
  end subroutine allocate_flux_registers


  subroutine init_time_variables()
                    
    integer ::  lev
    
    dt = 1.0_amrex_real
    last_regrid_step = 0
    t_new = 0.0_amrex_real
    t_old = -1.0_amrex_real
    stepno = 0
    nsubsteps(0) = 1
    do lev = 1, amrex_max_level
       nsubsteps(lev) = amrex_ref_ratio(lev-1)
    end do
    
  end subroutine init_time_variables

  
  subroutine init_bc_variables()
      
    ! Homogeneous Neumann boundary condition (foextrap implies that the ghost
    ! cells are filled with a copy of the closest point inside the domain)
    lo_bc = amrex_bc_foextrap
    hi_bc = amrex_bc_foextrap
    
  end subroutine init_bc_variables


  subroutine init_free_surface_variables()
      
    use read_input_module, only : sw_surf_pos_init
                                 
    integer ::  lo_x
    integer ::  hi_x

    lo_x = amrex_geom(amrex_max_level)%domain%lo(1)
    hi_x = amrex_geom(amrex_max_level)%domain%hi(1)
    
    melt_pos = sw_surf_pos_init
    melt_top = sw_surf_pos_init
    max_melt_vel = 0.0_amrex_real
    melt_vel = 0.0_amrex_real
    surf_dx(1) = amrex_geom(amrex_max_level)%dx(1)
    surf_ind(1,1) = lo_x
    surf_ind(1,2) = hi_x
    surf_current = 0.0_amrex_real
    surf_deformation = 0.0_amrex_real
    surf_enthalpy = 0.0_amrex_real
    surf_evap_flux = 0.0_amrex_real
    surf_pos = sw_surf_pos_init
    surf_pos_grid = sw_surf_pos_init
    surf_temperature = 0.0_amrex_real
    surf_xlo(1) = amrex_problo(1)
    
  end subroutine init_free_surface_variables
  
  ! ------------------------------------------------------------------
  ! Subroutine to free the memory
  ! ------------------------------------------------------------------
  subroutine amr_data_finalize()

    use read_input_module, only : deallocate_input
    
    call deallocate_time_variables
    call deallocate_bc_variables
    call deallocate_free_surface_variables
    call deallocate_multifabs
    call deallocate_flux_registers
    call deallocate_plasma_flux_variables

    ! Additionally, also free the memory related to the
    ! variable defined when reading the input file
    call deallocate_input
    
  end subroutine amr_data_finalize


  subroutine deallocate_time_variables()
    
    deallocate(dt)
    deallocate(last_regrid_step)
    deallocate(t_new)
    deallocate(t_old)
    deallocate(stepno)
    deallocate(nsubsteps)
    
  end subroutine deallocate_time_variables


  subroutine deallocate_bc_variables()

    deallocate(lo_bc)
    deallocate(hi_bc)
        
  end subroutine deallocate_bc_variables


  subroutine deallocate_free_surface_variables()

    deallocate(melt_pos)
    deallocate(melt_top)
    deallocate(melt_vel)
    deallocate(surf_current)
    deallocate(surf_deformation)  
    deallocate(surf_enthalpy)
    deallocate(surf_evap_flux)
    deallocate(surf_pos)
    deallocate(surf_pos_grid)
    deallocate(surf_temperature)
  
  end subroutine deallocate_free_surface_variables
  

  subroutine deallocate_multifabs()

    integer :: lev
     
    do lev = 0, amrex_max_level
           
       call amrex_multifab_destroy(idomain(lev))
       call amrex_multifab_destroy(phi_new(lev))
       call amrex_multifab_destroy(phi_old(lev))
       call amrex_multifab_destroy(temp(lev))
       call amrex_multifab_destroy(psi(lev))
    
    end do
    
  end subroutine deallocate_multifabs


  subroutine deallocate_flux_registers()

    integer :: lev

    do lev = 1, amrex_max_level
       call amrex_fluxregister_destroy(flux_reg(lev))
    end do
    
  end subroutine deallocate_flux_registers

  
  subroutine deallocate_plasma_flux_variables()

    if (allocated(plasma_flux_time_mesh)) deallocate(plasma_flux_time_mesh)
    if (allocated(plasma_flux_surf_mesh)) deallocate(plasma_flux_surf_mesh)
    if (allocated(plasma_flux_table)) deallocate(plasma_flux_table)
    if (allocated(plasma_flux_side_time_mesh)) deallocate(plasma_flux_side_time_mesh)
    if (allocated(plasma_flux_side_surf_mesh)) deallocate(plasma_flux_side_surf_mesh)
    if (allocated(plasma_flux_side_table)) deallocate(plasma_flux_side_table)
    
  end subroutine deallocate_plasma_flux_variables
  
end module amr_data_module
