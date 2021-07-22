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

  ! ------------------------------------------------------------------
  ! Subroutine to create a new AMReX level from scratch
  ! ------------------------------------------------------------------
  subroutine my_make_new_level_from_scratch (lev, time, pba, pdm) bind(c)
    
    use initdomain_module, only : init_phi
    use read_input_module, only : tempinit, surf_pos_init, do_reflux
    use material_properties_module, only : get_temp 

    ! Input and output variables
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm

    ! Local variables 
    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: phi(:,:,:,:), ptemp(:,:,:,:)

    ! Assign pointers for box array and distribution mapping
    ba = pba
    dm = pdm

    ! Set time
    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real

    ! Clear previous definition of the level
    call my_clear_level(lev)

    ! Build multifabs 
    call amrex_multifab_build(phi_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_old(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(temp(lev), ba, dm, ncomp, nghost)
    call amrex_imultifab_build(idomain_new(lev), ba, dm, ncomp, nghost)
    call amrex_imultifab_build(idomain_old(lev), ba, dm, ncomp, nghost)

    ! Build flux registers
    if (lev > 0 .and. do_reflux) then
       call amrex_fluxregister_build(flux_reg(lev), ba, dm, &
                                     amrex_ref_ratio(lev-1), lev, ncomp)
    end if

    ! Fill multifabs
    call amrex_mfiter_build(mfi, phi_new(lev))	
    do while (mfi%next())
       
       bx = mfi%tilebox()
       phi => phi_new(lev)%dataptr(mfi)
       ptemp => temp(lev)%dataptr(mfi)

       ! Initialize enthalpy
       call init_phi(lev, t_new(lev), bx%lo, bx%hi, tempinit, &
                     phi, lbound(phi), ubound(phi), &
                     amrex_geom(lev)%dx, amrex_problo)

       ! Initialize temperature (from the enthalpy)
       call get_temp(bx%lo, bx%hi, & 
       	             lbound(phi),   ubound(phi),   phi, &
                     lbound(ptemp), ubound(ptemp), ptemp)

       ! Initialize the idomain multifab (to be added)
       
    end do
    call amrex_mfiter_destroy(mfi)
    

  end subroutine my_make_new_level_from_scratch


  ! ------------------------------------------------------------------
  ! Subroutine to create a new AMReX level from a coarser level
  ! ------------------------------------------------------------------  
  subroutine my_make_new_level_from_coarse (lev, time, pba, pdm) bind(c)

    use fillpatch_module, only : fillcoarsepatch
    use read_input_module, only : surf_pos_init, do_reflux
    !use domain_module, only : get_idomain
    use material_properties_module, only : get_temp

    ! Input and output variables
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm

    ! Local variables
    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: phi(:,:,:,:), ptemp(:,:,:,:)
    !integer, contiguous, pointer :: idom(:,:,:,:) 

    ! Assign pointers for box array and distribution mapping
    ba = pba
    dm = pdm

    ! Clear previous definitions of the level
    call my_clear_level(lev)

    ! Set time
    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real

    ! Build multifabs
    call amrex_multifab_build(phi_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_old(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(temp(lev), ba, dm, ncomp, nghost)
    call amrex_imultifab_build(idomain_new(lev), ba, dm, ncomp, nghost)
    call amrex_imultifab_build(idomain_old(lev), ba, dm, ncomp, nghost)

    ! Build flux registers
    if (lev > 0 .and. do_reflux) then
       call amrex_fluxregister_build(flux_reg(lev), ba, dm, amrex_ref_ratio(lev-1), lev, ncomp)
    end if

    ! Fill the enthalpy multifab from the coarse level 
    call fillcoarsepatch(lev, time, phi_new(lev))
    
    ! Fill other multifabs
    call amrex_mfiter_build(mfi, temp(lev))
    do while (mfi%next())
       
       bx = mfi%tilebox()
       phi => phi_new(lev)%dataptr(mfi)
       ptemp => temp(lev)%dataptr(mfi)

       ! Initialize temperature (from the enthalpy) 
       call get_temp(bx%lo, bx%hi, & 
                     lbound(phi),   ubound(phi),   phi, &
                     lbound(ptemp), ubound(ptemp), ptemp)

       ! Initialize the idomain multifab (maybe add to fillpatch)
       !idom => idomain_new(lev)%dataptr(mfi)
       ! call get_idomain(bx%lo, bx%hi, & 
       ! 		lbound(idom), ubound(idom), idom)
       
    end do    
    call amrex_mfiter_destroy(mfi)
   
  end subroutine my_make_new_level_from_coarse

  
  ! ------------------------------------------------------------------
  ! Subroutine to copy a new AMReX level from a coarser level
  ! ------------------------------------------------------------------  
  subroutine my_remake_level (lev, time, pba, pdm) bind(c)

    use fillpatch_module, only : fillpatch
    use read_input_module, only : surf_pos_init, do_reflux
    use material_properties_module, only : get_temp     

    ! Input and output variables
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm

    ! Local variables
    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm
    type(amrex_multifab) :: new_phi_new
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: phi(:,:,:,:), ptemp(:,:,:,:)
    
    ! Assign pointers for box array and distribution mapping
    ba = pba
    dm = pdm
	
    ! Clear previous definitions of the level
    call my_clear_level(lev)

    ! Set time
    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real

    ! Build multifabs
    call amrex_multifab_build(new_phi_new, ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_old(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(temp(lev), ba, dm, ncomp, nghost)
    call amrex_imultifab_build(idomain_new(lev), ba, dm, ncomp, nghost)
    call amrex_imultifab_build(idomain_old(lev), ba, dm, ncomp, nghost)

    ! Build flux registers
    if (lev > 0 .and. do_reflux) then
       call amrex_fluxregister_build(flux_reg(lev), ba, dm, amrex_ref_ratio(lev-1), lev, ncomp)
    end if

    ! Make a copy of the content of new_phi_new into phi_new
    call fillpatch(lev, time, new_phi_new)
    call phi_new(lev)%copy(new_phi_new, 1, 1, ncomp, 0)
    call amrex_multifab_destroy(new_phi_new)
    
    ! Fill other multifabs
    call amrex_mfiter_build(mfi, temp(lev))
    do while (mfi%next())
       
       bx = mfi%tilebox()
       phi => phi_new(lev)%dataptr(mfi)
       ptemp => temp(lev)%dataptr(mfi)
       call get_temp(bx%lo, bx%hi, & 
       	             lbound(phi),   ubound(phi),   phi, &
                     lbound(ptemp), ubound(ptemp), ptemp)
       
    end do
    call amrex_mfiter_destroy(mfi)

    
  end subroutine my_remake_level

  ! ------------------------------------------------------------------
  ! Subroutine to clear a level
  ! ------------------------------------------------------------------  
  subroutine my_clear_level(lev) bind(c)
    
    integer, intent(in), value :: lev
    
    call amrex_multifab_destroy(phi_new(lev))
    call amrex_multifab_destroy(phi_old(lev))
    call amrex_multifab_destroy(temp(lev))
    call amrex_imultifab_destroy(idomain_new(lev))
    call amrex_imultifab_destroy(idomain_old(lev))
    call amrex_fluxregister_destroy(flux_reg(lev))
    
  end subroutine my_clear_level

  ! ------------------------------------------------------------------
  ! Subroutine used to construct the error estimate employed to 
  ! the grid points where it is necessary to regrid
  ! ------------------------------------------------------------------ 
  subroutine my_error_estimate(lev, cp, t, settag, cleartag) bind(c)

    use tagging_module, only : tag_phi_error
    use read_input_module,  only : surfdist 

    ! Input and output variables 
    integer, intent(in), value :: lev
    type(c_ptr), intent(in), value :: cp
    real(amrex_real), intent(in), value :: t
    character(kind=c_char), intent(in), value :: settag, cleartag

    ! Local variables
    type(amrex_geometry) :: geom  	! geometry at level
    real(amrex_real), allocatable, save :: phierr(:)
    type(amrex_parmparse) :: pp
    type(amrex_tagboxarray) :: tag
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: phi(:,:,:,:)
    character(kind=c_char), contiguous, pointer :: tag(:,:,:,:)

    ! Set pointers
    tag = cp

    ! Geometry
    geom = amrex_geom(lev) 


    
    !$omp parallel private(mfi, bx, phi, tag)
    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.true.)
    do while(mfi%next())
       
       bx = mfi%tilebox()
       phi => phi_new(lev)%dataptr(mfi)
       tag => tag%dataptr(mfi)
       call tag_phi_error(lev, t, bx%lo, bx%hi, &
                          geom%get_physical_location(bx%lo), &
                          geom%dx, surfdist(lev+1), & 
                          phi, lbound(phi), ubound(phi), &
                          tag, lbound(tag), ubound(tag), &
                          settag, cleartag)
       
    end do
    call amrex_mfiter_destroy(mfi)
    !$omp end parallel

  end subroutine my_error_estimate

  
end module amr_data_module
