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
  ! Variables pertaining to physical solution
  ! ------------------------------------------------------------------
  public :: t_new, t_old ! Time 
  public :: phi_new, phi_old ! Enthalpy
  public :: temp ! Temperature
  public :: idomain_new, idomain_old ! Indexes used to distinguish between material and background
  public :: surf_ind, surf_xlo, surf_dx ! 2D surface grid parameters 
  public :: melt_pos, surf_pos, melt_vel  ! Melt position, free surface position and melt velocity
  ! ------------------------------------------------------------------
  ! Variables and subroutines pertaining to AMReX calculations
  ! ------------------------------------------------------------------
  public :: flux_reg 
  public :: amr_data_init, amr_data_finalize


  integer  :: surf_ind(2,2) = 0 ! fluid domain index bounds (x,z) (lo,hi)
  
  real(rt), allocatable :: t_new(:)
  real(rt), allocatable :: t_old(:)
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
    
    allocate(t_new(0:amrex_max_level))
    t_new = 0.0_rt ! This should be moved to where the variables are initialized
    allocate(t_old(0:amrex_max_level))
    t_old = -1.0e100_rt ! This should be moved to where the variables are initialized
    allocate(phi_new(0:amrex_max_level))
    allocate(phi_old(0:amrex_max_level))
    allocate(temp(0:amrex_max_level))
    allocate(idomain_new(0:amrex_max_level))
    allocate(idomain_old(0:amrex_max_level))
    allocate(flux_reg(0:amrex_max_level))
    
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
