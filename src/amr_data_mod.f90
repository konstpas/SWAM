module amr_data_module
  
  ! -----------------------------------------------------------------
  ! This module is used to declare, allocate and destroy the global
  ! variables employed in the rest of the code
  ! -----------------------------------------------------------------
  
  use iso_c_binding
  use amrex_amr_module
  use amrex_fort_module, only : rt => amrex_real
  
  implicit none

  private

  ! ------------------------------------------------------------------
  ! Variables for the physical solution
  ! ------------------------------------------------------------------
  public :: t_new, t_old ! Time
  public :: dt ! Timestep
  public :: phi_new, phi_old ! Enthalpy
  public :: temp ! Temperature
  public :: idomain_new, idomain_old ! Indexes used to distinguish between material and background
  public :: surf_ind, surf_xlo, surf_dx ! 2D surface grid parameters 
  public :: melt_pos, surf_pos, melt_vel  ! Melt position, free surface position and melt velocity

  ! ------------------------------------------------------------------
  ! Variables for the AMReX calculations
  ! ------------------------------------------------------------------
  public :: flux_reg
  public :: stepno
  public :: nsubsteps

  ! ------------------------------------------------------------------
  ! Public subroutines
  ! ------------------------------------------------------------------
  public :: amr_data_init, amr_data_finalize
  
  ! ------------------------------------------------------------------
  ! Declare public variables
  ! ------------------------------------------------------------------
  integer  :: surf_ind(2,2) = 0 ! fluid domain index bounds (x,z) (lo,hi)
  integer, allocatable, save :: stepno(:)
  integer, allocatable, save :: nsubsteps(:)
  
  real(rt), allocatable :: t_new(:)
  real(rt), allocatable :: t_old(:)
  real(rt), allocatable, save :: dt(:)
  real(rt), allocatable :: surf_pos(:,:), melt_pos(:,:), melt_vel(:,:,:)
  real(rt) :: surf_xlo(2), surf_dx(2) 

  type(amrex_fluxregister), allocatable :: flux_reg(:)
  
  type(amrex_imultifab), allocatable :: idomain_new(:)
  type(amrex_imultifab), allocatable :: idomain_old(:)
  
  type(amrex_multifab), allocatable :: phi_new(:)
  type(amrex_multifab), allocatable :: phi_old(:)
  type(amrex_multifab), allocatable :: temp(:)


  ! --------------------------------------------------------
  ! Declare private constants
  ! --------------------------------------------------------
  integer, private, parameter :: ncomp = 1, nghost = 0
  
  
contains

  ! -----------------------------------------------
  ! Subroutine to allocate and initialize variables
  ! -----------------------------------------------
  subroutine amr_data_init()

    use read_input_module
    
    integer :: lo_x, hi_x, lo_z, hi_z, lev
    
    lo_x = amrex_geom(amrex_max_level)%domain%lo(1)
    hi_x = amrex_geom(amrex_max_level)%domain%hi(1)
    lo_z = amrex_geom(amrex_max_level)%domain%lo(3)
    hi_z = amrex_geom(amrex_max_level)%domain%hi(3)

    ! Allocate variables for physical solution
    allocate(t_new(0:amrex_max_level))
    allocate(t_old(0:amrex_max_level))
    allocate(dt(0:amrex_max_level))
    allocate(phi_new(0:amrex_max_level))
    allocate(phi_old(0:amrex_max_level))
    allocate(temp(0:amrex_max_level))
    allocate(idomain_new(0:amrex_max_level))
    allocate(idomain_old(0:amrex_max_level))
    allocate(surf_pos(lo_x:hi_x, lo_z:hi_z))
    allocate(melt_pos(lo_x:hi_x, lo_z:hi_z))	
    allocate(melt_vel(lo_x:hi_x+1, lo_z:hi_z+1, 1:amrex_spacedim-1))

    ! Allocate variables for AMReX calculations
    allocate(flux_reg(0:amrex_max_level))
    allocate(stepno(0:amrex_max_level))
    allocate(nsubsteps(0:amrex_max_level))

    ! Initialize variables for physical solution (CHECK IF THIS IS REALLY NECESSARY)
    ! The multifabs are initialized when the levels are created.
    ! This is done in the routines called my_make_new_level_...
    t_new = 0.0_rt
    t_old = -1.0e100_rt
    dt = huge(1._rt)
    surf_pos = surf_pos_init 	
    melt_pos = surf_pos_init 
    melt_vel = 0.0_rt 	 
    surf_xlo(1) = amrex_problo(1) 
    surf_xlo(2) = amrex_problo(3)
    surf_dx(1) = amrex_geom(amrex_max_level)%dx(1)
    surf_dx(2) = amrex_geom(amrex_max_level)%dx(3)
    surf_ind(1,1) = lo_x
    surf_ind(1,2) = hi_x
    surf_ind(2,1) = lo_z
    surf_ind(2,2) = hi_z
    
    ! Initialize variables for AMReX calculations (CHECK IF THIS IS REALLY NECESSARY)
    ! The flux_registers are initialized when the levels are created
    ! This is done in the routines called my_make_new_level_...
    stepno = 0
    nsubsteps(0) = 1
    do lev = 1, amrex_max_level
       nsubsteps(lev) = amrex_ref_ratio(lev-1)
    end do        
    
    
  end subroutine amr_data_init


  ! ------------------------------------------------------------------
  ! Subroutine to free multifab variables for physical solution
  ! ------------------------------------------------------------------
  subroutine amr_data_finalize()
    
    integer :: lev
    
    do lev = 0, amrex_max_level
       call amrex_multifab_destroy(phi_new(lev))
       call amrex_multifab_destroy(phi_old(lev))
       call amrex_multifab_destroy(temp(lev))
       call amrex_imultifab_destroy(idomain_new(lev))
       call amrex_imultifab_destroy(idomain_old(lev))
    end do
    
    do lev = 1, amrex_max_level
       call amrex_fluxregister_destroy(flux_reg(lev))
    end do
    
  end subroutine amr_data_finalize

  
end module amr_data_module
