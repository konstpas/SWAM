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
  public :: amr_data_init, amr_data_finalize


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

  
contains

  ! ------------------------------------------------------------------
  ! Subroutine to allocate variables for physical solution
  ! ------------------------------------------------------------------ 
  subroutine amr_data_init()
    
    use bc_module, only : lo_bc, hi_bc
    use read_input_module, only : surf_pos_init
    
    integer ::  lev
    
    allocate(t_new(0:amrex_max_level))
    allocate(t_old(0:amrex_max_level))
    allocate(phi_new(0:amrex_max_level))
    allocate(phi_old(0:amrex_max_level))
    allocate(temp(0:amrex_max_level))
    allocate(idomain_new(0:amrex_max_level))
    allocate(idomain_old(0:amrex_max_level))
    allocate(flux_reg(0:amrex_max_level))

    t_old = -1.0e100_rt ! This should be moved to where the variables are initialized
    t_new = 0.0_rt ! This should be moved to where the variables are initialized

    if (.not. amrex_is_all_periodic()) then
       lo_bc = amrex_bc_foextrap
       hi_bc = amrex_bc_foextrap
    end if

    allocate(stepno(0:amrex_max_level))
    stepno = 0

    allocate(nsubsteps(0:amrex_max_level))
    nsubsteps(0) = 1
    do lev = 1, amrex_max_level
       nsubsteps(lev) = amrex_ref_ratio(lev-1)
    end do

    allocate(dt(0:amrex_max_level))
    dt = huge(1._rt)

	
    ! 2D domain initialized 
    surf_xlo(1) = amrex_problo(1) 
    surf_xlo(2) = amrex_problo(3) ! 2nd dimension in fluid is z direction, 3rd dimension in 3d. 
    
    surf_dx(1) = amrex_geom(amrex_max_level)%dx(1)!
    surf_dx(2) = amrex_geom(amrex_max_level)%dx(3)!
    
    ! Initiate surface position array (y-position as function of x and z)
    allocate(surf_pos(amrex_geom(amrex_max_level)%domain%lo(1):amrex_geom(amrex_max_level)%domain%hi(1), &
         amrex_geom(amrex_max_level)%domain%lo(3):amrex_geom(amrex_max_level)%domain%hi(3)))
    ! Initiate melt interface position array (y-position as function of x and z)
    allocate(melt_pos(amrex_geom(amrex_max_level)%domain%lo(1):amrex_geom(amrex_max_level)%domain%hi(1), &
         amrex_geom(amrex_max_level)%domain%lo(3):amrex_geom(amrex_max_level)%domain%hi(3)))	
    ! Initiate melt velocity array (2D velocity, index 1 x, index 2 z)
    allocate(melt_vel(amrex_geom(amrex_max_level)%domain%lo(1):amrex_geom(amrex_max_level)%domain%hi(1)+1, &
         amrex_geom(amrex_max_level)%domain%lo(3):amrex_geom(amrex_max_level)%domain%hi(3)+1, & 
         1:amrex_spacedim-1))	
    
    surf_pos = surf_pos_init  ! Initial value of surface array 	
    melt_pos = surf_pos_init     ! Initial value of melt position array (no melting) 
    melt_vel = 0 		! Initial value of melt velocity array 
    
    surf_ind(1,1) = amrex_geom(amrex_max_level)%domain%lo(1)
    surf_ind(1,2) = amrex_geom(amrex_max_level)%domain%hi(1)
    surf_ind(2,1) = amrex_geom(amrex_max_level)%domain%lo(3)
    surf_ind(2,2) = amrex_geom(amrex_max_level)%domain%hi(3) 
    
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
