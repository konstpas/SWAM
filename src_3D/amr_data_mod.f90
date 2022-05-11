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
  ! Time step (one for each level)
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
  ! Boundary conditions of the entire simulation box
  public :: lo_bc
  public :: hi_bc
  ! Bottom of the melt pool
  public :: melt_pos
  ! Top of the melt pool
  public :: melt_top
  ! Melt velocity
  public :: melt_vel 
  ! Maximum melt velocity in x direction(absolute value)
  public :: max_melt_vel_x
  ! Maximum melt velocity in z direction(absolute value)
  public :: max_melt_vel_z
  ! Tabulated plasma heat flux
  public :: plasma_flux_time_mesh
  public :: plasma_flux_surf_x_mesh
  public :: plasma_flux_surf_z_mesh
  public :: heat_flux_table
  ! Tabulated plasma heat flux on second exposed surface
  public :: plasma_side_flux_time_mesh
  public :: plasma_side_flux_surf_y_mesh
  public :: plasma_side_flux_surf_z_mesh
  public :: heat_side_flux_table
  ! Resolution of the shallow water grid
  public :: surf_dx 
  ! Index of the shallow water grid
  public :: surf_ind 
  ! Lower corner of the shallow water domain
  public :: surf_xlo 
  ! Position of the free surface
  public :: surf_pos 
  ! Temperature on the free surface
  public :: surf_temperature 
  ! Orientation of the surface normal
  public :: surf_normal
  ! Enthalpy on the free surface
  public :: surf_enthalpy 
  ! Evaporation flux on the free surface
  public :: surf_evap_flux 
  ! Variable to store the solution of the geoclaw solver
  public :: qnew 
  ! Thermionic current on the free surface
  public :: surf_current 
  public :: domain_top ! Position of the top of the simulation domain
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
  ! Variable used to store the total energy poured into the system
  public :: Qplasma
  ! Variables used to store the total energy lost to the cooling pipe, thermionic emission, radiation and vaporization
  public :: Qpipe
  public :: Qtherm
  public :: Qrad
  public :: Qvap

  
  ! ------------------------------------------------------------------
  ! Declare public variables
  ! ------------------------------------------------------------------

  ! Variables related to time
  integer, allocatable, save :: stepno(:)
  integer, allocatable, save :: last_regrid_step(:)
  integer, allocatable, save :: nsubsteps(:)
  real(amrex_real), allocatable, save :: dt(:)
  real(amrex_real), allocatable, save :: t_new(:)
  real(amrex_real), allocatable, save :: t_old(:)

  ! Variables related to boundary conditions
  integer, allocatable, save :: lo_bc(:,:)
  integer, allocatable, save :: hi_bc(:,:)
  real(amrex_real), save :: domain_top

  ! Variables related to the free surface and melt pool
  integer, save  :: surf_ind(2,2)
  real(amrex_real), allocatable, save :: melt_pos(:,:)
  real(amrex_real), allocatable, save :: melt_top(:,:)
  real(amrex_real), allocatable, save :: melt_vel(:,:,:)
  real(amrex_real), save :: max_melt_vel_x
  real(amrex_real), save :: max_melt_vel_z
  real(amrex_real), allocatable, save :: qnew(:,:,:)
  real(amrex_real), allocatable, save :: surf_evap_flux(:,:)
  real(amrex_real), allocatable, save :: surf_pos(:,:)
  real(amrex_real), allocatable, save :: surf_temperature(:,:)
  real(amrex_real), allocatable, save :: surf_normal(:,:)
  real(amrex_real), allocatable, save :: surf_enthalpy(:,:)
  real(amrex_real), save :: surf_xlo(2)
  real(amrex_real), save :: surf_dx(2)
  real(amrex_real), allocatable, save :: surf_current(:,:)
  real(amrex_real), allocatable, save :: surf_deformation(:,:)

 
  ! Variables related to the heat flux read from files
  real(amrex_real), allocatable, dimension(:), save :: plasma_flux_time_mesh
  real(amrex_real), allocatable, dimension(:), save :: plasma_flux_surf_x_mesh
  real(amrex_real), allocatable, dimension(:), save :: plasma_flux_surf_z_mesh
  real(amrex_real), allocatable, dimension(:,:,:), save :: heat_flux_table
  real(amrex_real), allocatable, dimension(:), save :: plasma_side_flux_time_mesh
  real(amrex_real), allocatable, dimension(:), save :: plasma_side_flux_surf_y_mesh
  real(amrex_real), allocatable, dimension(:), save :: plasma_side_flux_surf_z_mesh
  real(amrex_real), allocatable, dimension(:,:,:), save :: heat_side_flux_table
  
  ! Multifabs
  type(amrex_multifab), allocatable, save :: idomain(:)
  type(amrex_multifab), allocatable, save :: phi_new(:)
  type(amrex_multifab), allocatable, save :: phi_old(:)
  type(amrex_multifab), allocatable, save :: temp(:)

  ! Flux registers
  type(amrex_fluxregister), allocatable, save :: flux_reg(:)

  ! Variables that keep track of the energy balance
  real(amrex_real), save :: Qplasma
  real(amrex_real), save :: Qpipe
  real(amrex_real), save :: Qtherm
  real(amrex_real), save :: Qrad
  real(amrex_real), save :: Qvap

  
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
    integer ::  lo_z
    integer ::  hi_z

    lo_x = amrex_geom(amrex_max_level)%domain%lo(1)
    hi_x = amrex_geom(amrex_max_level)%domain%hi(1)
    lo_z = amrex_geom(amrex_max_level)%domain%lo(3)
    hi_z = amrex_geom(amrex_max_level)%domain%hi(3)


    allocate(qnew(1:3,lo_x-3:hi_x+3,lo_z-3:hi_z+3))
    allocate(surf_evap_flux(lo_x:hi_x,lo_z:hi_z))
    allocate(surf_pos(lo_x:hi_x,lo_z:hi_z))
    allocate(surf_temperature(lo_x:hi_x,lo_z:hi_z))
    allocate(surf_normal(lo_x:hi_x,lo_z:hi_z))
    allocate(surf_enthalpy(lo_x:hi_x,lo_z:hi_z))
    allocate(melt_pos(lo_x:hi_x, lo_z:hi_z))
    allocate(melt_top(lo_x:hi_x, lo_z:hi_z))
    allocate(melt_vel(lo_x:hi_x+1,lo_z:hi_z+1, 1:amrex_spacedim-1))
    allocate(surf_current(lo_x:hi_x, lo_z:hi_z))
    allocate(surf_deformation(lo_x:hi_x,lo_z:hi_z))  

  end subroutine allocate_free_surface_variables


  subroutine allocate_multifabs()

    allocate(idomain(0:amrex_max_level))
    allocate(phi_new(0:amrex_max_level))
    allocate(phi_old(0:amrex_max_level))
    allocate(temp(0:amrex_max_level))
    
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
      
    use read_input_module, only : sw_surf_pos_init, &
                                  sw_read_free_surface_file
    use read_free_surface_module, only :  get_free_surface
                                 
    integer ::  lo_x
    integer ::  hi_x        
    integer ::  lo_z
    integer ::  hi_z

    lo_x = amrex_geom(amrex_max_level)%domain%lo(1)
    hi_x = amrex_geom(amrex_max_level)%domain%hi(1)
    lo_z = amrex_geom(amrex_max_level)%domain%lo(3)
    hi_z = amrex_geom(amrex_max_level)%domain%hi(3)

    surf_dx(1) = amrex_geom(amrex_max_level)%dx(1)
    surf_dx(2) = amrex_geom(amrex_max_level)%dx(3) 
    surf_xlo(1) = amrex_problo(1) 
    surf_xlo(2) = amrex_problo(3)  
    surf_ind(1,1) = lo_x
    surf_ind(1,2) = hi_x
    surf_ind(2,1) = lo_z
    surf_ind(2,2) = hi_z 
    if(sw_read_free_surface_file) then
      call get_free_surface(surf_ind, surf_xlo, surf_dx, melt_top, melt_pos, surf_pos)
    else
      melt_pos = sw_surf_pos_init
      melt_top = sw_surf_pos_init
      surf_pos = sw_surf_pos_init
    end if
    melt_vel = 0.0_amrex_real
    max_melt_vel_x = 0.0_amrex_real
    max_melt_vel_z = 0.0_amrex_real
    qnew(1,:,:) = 0.0_amrex_real
    qnew(2,:,:) = 0.0_amrex_real
    qnew(3,:,:) = 0.0_amrex_real
    surf_evap_flux = 0.0_amrex_real
    surf_temperature = 0.0_amrex_real
    surf_normal = 0.0_amrex_real
    surf_enthalpy = 0.0_amrex_real
    domain_top = amrex_probhi(2)
    surf_current = 0.0_amrex_real
    surf_deformation = 0.0_amrex_real
    
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
    deallocate(surf_temperature)
    deallocate(surf_normal)
  
  end subroutine deallocate_free_surface_variables
  

  subroutine deallocate_multifabs()

    integer :: lev
     
    do lev = 0, amrex_max_level
           
       call amrex_multifab_destroy(idomain(lev))
       call amrex_multifab_destroy(phi_new(lev))
       call amrex_multifab_destroy(phi_old(lev))
       call amrex_multifab_destroy(temp(lev))
    
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
    if (allocated(plasma_flux_surf_x_mesh)) deallocate(plasma_flux_surf_x_mesh)
    if (allocated(plasma_flux_surf_z_mesh)) deallocate(plasma_flux_surf_z_mesh)
    if (allocated(heat_flux_table)) deallocate(heat_flux_table)
    if (allocated(plasma_side_flux_time_mesh)) deallocate(plasma_side_flux_time_mesh)
    if (allocated(plasma_side_flux_surf_y_mesh)) deallocate(plasma_side_flux_surf_y_mesh)
    if (allocated(plasma_side_flux_surf_z_mesh)) deallocate(plasma_side_flux_surf_z_mesh)
    if (allocated(heat_side_flux_table)) deallocate(heat_side_flux_table)
    
  end subroutine deallocate_plasma_flux_variables

  
end module amr_data_module
