module regridding_module
  
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
  ! Public subroutines
  ! ------------------------------------------------------------------
  public :: my_make_new_level_from_scratch
  public :: my_make_new_level_from_coarse
  public :: my_remake_level
  public :: my_clear_level
  public :: my_error_estimate
  
  ! --------------------------------------------------------
  ! Declare private constants
  ! --------------------------------------------------------
  integer, private, parameter :: ncomp = 1, nghost = 0
  
  
contains

  ! ------------------------------------------------------------------
  ! Subroutine to create a new AMReX level from scratch
  ! ------------------------------------------------------------------
  subroutine my_make_new_level_from_scratch (lev, time, pba, pdm) bind(c)

    use amr_data_module, only : phi_new, phi_old, temp, idomain_new, &
                                idomain_old, t_new, t_old, flux_reg
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

    use amr_data_module, only : phi_new, phi_old, temp, idomain_new, &
                                idomain_old, t_new, t_old, flux_reg
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

    use amr_data_module, only : phi_new, phi_old, temp, idomain_new, &
                                idomain_old, t_new, t_old, flux_reg
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

    use amr_data_module, only : phi_new, phi_old, temp, idomain_new, &
                                idomain_old, t_new, t_old, flux_reg
    
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

    use amr_data_module, only : phi_new
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
    character(kind=c_char), contiguous, pointer :: tagarr(:,:,:,:)

    ! Set pointers
    tag = cp

    ! Geometry
    geom = amrex_geom(lev) 
    
    !$omp parallel private(mfi, bx, phiarr, tagarr)
    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.true.)
    do while(mfi%next())

       bx = mfi%tilebox()
       phi => phi_new(lev)%dataptr(mfi)
       tagarr => tag%dataptr(mfi)
       call tag_phi_error(lev, t, bx%lo, bx%hi, &
                          geom%get_physical_location(bx%lo), &
                          geom%dx, surfdist(lev+1), & 
                          phi, lbound(phi), ubound(phi), &
                          tagarr, lbound(tagarr), ubound(tagarr), &
                          settag, cleartag)
       
    end do
    call amrex_mfiter_destroy(mfi)
    !$omp end parallel

  end subroutine my_error_estimate

  
end module regridding_module
