module my_amr_module

  use iso_c_binding
  use amrex_amr_module
  use amrex_fort_module, only : rt => amrex_real

  use amr_data_module
  implicit none
  
  integer, private, parameter :: ncomp = 1, nghost = 0
  
contains

  subroutine my_amr_init ()
    
    use bc_module, only : lo_bc, hi_bc
    use read_input_module
    use material_properties_module
    
    type(amrex_parmparse) :: pp
    integer :: ilev
    
    if (.not.amrex_amrcore_initialized()) call amrex_amrcore_init()
    
    call amrex_init_virtual_functions (my_make_new_level_from_scratch, &
         &                             my_make_new_level_from_coarse,  &
         &                             my_remake_level,                &
         &                             my_clear_level,                 &
         &                             my_error_estimate)

    call read_input_file
    ! ! some default parameters
    ! allocate(character(len=3)::check_file)
    ! check_file = "chk"
    ! allocate(character(len=3)::plot_file)
    ! plot_file = "plt"
    ! allocate(character(len=0)::restart)

    ! ! Read parameters
    ! call amrex_parmparse_build(pp)
    ! call pp%query("max_step", max_step)
    ! call pp%query("stop_time", stop_time)
    ! call amrex_parmparse_destroy(pp)
    
    ! ! Parameters amr.*
    ! call amrex_parmparse_build(pp, "amr")
    ! call pp%query("regrid_int", regrid_int)
    ! call pp%query("check_int", check_int)
    ! call pp%query("plot_int", plot_int)
    ! call pp%query("check_file", check_file)
    ! call pp%query("plot_file", plot_file)
    ! call pp%query("restart", restart)
    ! call pp%query("max_grid_size_2d", max_grid_size_2d)
    ! call amrex_parmparse_destroy(pp)
   
    ! ! Domain parameters 
    ! call amrex_parmparse_build(pp, "domain")
    ! call pp%query("meltvel", meltvel)  
    ! call pp%query("tempinit", tempinit) 
    ! call pp%query("surf_pos", surf_pos_init)   
    ! call pp%getarr("surfdist", surfdist)
    ! call pp%query("flux_peak", flux_peak) 
    ! call pp%getarr("flux_pos", flux_pos)
    ! call pp%getarr("flux_width", flux_width)
    ! call pp%query("exp_time", exp_time) 
    ! call amrex_parmparse_destroy(pp)

    ! ! Material parameters
    ! call amrex_parmparse_build(pp, "material")
    ! call pp%query("material", material)
    ! call pp%query("phiT_max_T", phiT_table_max_T)
    ! call pp%query("phiT_n_points", phiT_table_n_points)
    ! call amrex_parmparse_destroy(pp)
    ! call init_mat_prop ! Initialize material properties 
    
    ! ! Parameters myamr.*
    ! call amrex_parmparse_build(pp, "myamr")
    ! call pp%query("v", verbose)
    ! call pp%query("verbose", verbose)
    ! call pp%query("cfl", cfl)
    ! call pp%query("do_reflux", do_reflux)
    ! call amrex_parmparse_destroy(pp)

    call init_mat_prop
    
    if (.not. amrex_is_all_periodic()) then
       lo_bc = amrex_bc_foextrap
       hi_bc = amrex_bc_foextrap
    end if

    allocate(stepno(0:amrex_max_level))
    stepno = 0

    allocate(nsubsteps(0:amrex_max_level))
    nsubsteps(0) = 1
    do ilev = 1, amrex_max_level
       nsubsteps(ilev) = amrex_ref_ratio(ilev-1)
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
	
    call amr_data_init()
  end subroutine my_amr_init


  subroutine my_amr_finalize ()
    call amr_data_finalize()
  end subroutine my_amr_finalize

  ! Make a new level from scratch and put the data in phi_new.
  ! Note that phi_old contains no valid data after this.
  subroutine my_make_new_level_from_scratch (lev, time, pba, pdm) bind(c)
    use initdomain_module, only : init_phi
    use read_input_module, only : tempinit, surf_pos_init, do_reflux
    use material_properties_module, only : get_temp 
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm

    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: phi(:,:,:,:), ptemp(:,:,:,:)
    
    ba = pba
    dm = pdm

    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real

    call my_clear_level(lev)

    call amrex_multifab_build(phi_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_old(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(temp(lev), ba, dm, ncomp, nghost)
    call amrex_imultifab_build(idomain_new(lev), ba, dm, ncomp, nghost)
    call amrex_imultifab_build(idomain_old(lev), ba, dm, ncomp, nghost)


   if (lev > 0 .and. do_reflux) then
      call amrex_fluxregister_build(flux_reg(lev), ba, dm, amrex_ref_ratio(lev-1), lev, ncomp)
   end if

    call amrex_mfiter_build(mfi, phi_new(lev))
	
    do while (mfi%next())
       bx = mfi%tilebox()
       phi => phi_new(lev)%dataptr(mfi)
       ptemp => temp(lev)%dataptr(mfi)
       call init_phi(lev, t_new(lev), bx%lo, bx%hi, tempinit, phi, lbound(phi), ubound(phi), &
            amrex_geom(lev)%dx, amrex_problo)
       call get_temp(bx%lo, bx%hi, & 
       	lbound(phi),   ubound(phi),   phi, &
       	lbound(ptemp), ubound(ptemp), ptemp)
       
       ! call get domain integers 	
    end do

    call amrex_mfiter_destroy(mfi)
    

  end subroutine my_make_new_level_from_scratch

  ! Make a new level from coarse level and put the data in phi_new.
  ! Note tha phi_old contains no valid data after this.
  subroutine my_make_new_level_from_coarse (lev, time, pba, pdm) bind(c)
    use fillpatch_module, only : fillcoarsepatch
    use read_input_module, only : surf_pos_init, do_reflux
    use domain_module, only : get_idomain
    use material_properties_module, only : get_temp 
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm

    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: phi(:,:,:,:), ptemp(:,:,:,:)
    integer, contiguous, pointer :: idom(:,:,:,:) 

    ba = pba
    dm = pdm

    call my_clear_level(lev)

    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real

    call amrex_multifab_build(phi_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_old(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(temp(lev), ba, dm, ncomp, nghost)
    call amrex_imultifab_build(idomain_new(lev), ba, dm, ncomp, nghost)
    call amrex_imultifab_build(idomain_old(lev), ba, dm, ncomp, nghost)
    
    if (lev > 0 .and. do_reflux) then
       call amrex_fluxregister_build(flux_reg(lev), ba, dm, amrex_ref_ratio(lev-1), lev, ncomp)
    end if

     
    call fillcoarsepatch(lev, time, phi_new(lev))
    
    ! Fill temperature data 
    call amrex_mfiter_build(mfi, temp(lev))
    ! parallelize 
    do while (mfi%next())
       bx = mfi%tilebox()
       phi => phi_new(lev)%dataptr(mfi)
       ptemp => temp(lev)%dataptr(mfi)
       idom => idomain_new(lev)%dataptr(mfi)
       call get_temp(bx%lo, bx%hi, & 
       	lbound(phi),   ubound(phi),   phi, &
       	lbound(ptemp), ubound(ptemp), ptemp)  
       call get_idomain(bx%lo, bx%hi, & 
       		lbound(idom), ubound(idom), idom) 	
       ! call get domain integers 
    end do    
    
    call amrex_mfiter_destroy(mfi)
   
   
    
  end subroutine my_make_new_level_from_coarse

  ! Remake a level from current and coarse levels and put the data in phi_new.
  ! Note that phi_old contains no valid data after this.
  subroutine my_remake_level (lev, time, pba, pdm) bind(c)
    use fillpatch_module, only : fillpatch
    use read_input_module, only : surf_pos_init, do_reflux
    use material_properties_module, only : get_temp     
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm
    
    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm
    type(amrex_multifab) :: new_phi_new
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: phi(:,:,:,:), ptemp(:,:,:,:)
    

    ba = pba
    dm = pdm
	
    call amrex_multifab_build(new_phi_new, ba, dm, ncomp, 0)
    call fillpatch(lev, time, new_phi_new)

    call my_clear_level(lev)

    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real

    call amrex_multifab_build(phi_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_old(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(temp(lev), ba, dm, ncomp, nghost)
    call amrex_imultifab_build(idomain_new(lev), ba, dm, ncomp, nghost)
    call amrex_imultifab_build(idomain_old(lev), ba, dm, ncomp, nghost)
    
    if (lev > 0 .and. do_reflux) then
       call amrex_fluxregister_build(flux_reg(lev), ba, dm, amrex_ref_ratio(lev-1), lev, ncomp)
    end if

    call phi_new(lev)%copy(new_phi_new, 1, 1, ncomp, 0)
    
    call amrex_multifab_destroy(new_phi_new)
    
    ! Fill temperature data 
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

  subroutine my_clear_level (lev) bind(c)
    integer, intent(in), value :: lev
    call amrex_multifab_destroy(phi_new(lev))
    call amrex_multifab_destroy(phi_old(lev))
    call amrex_multifab_destroy(temp(lev))
    call amrex_imultifab_destroy(idomain_new(lev))
    call amrex_imultifab_destroy(idomain_old(lev))
    call amrex_fluxregister_destroy(flux_reg(lev))
  end subroutine my_clear_level

  subroutine my_error_estimate (lev, cp, t, settag, cleartag) bind(c)

    use tagging_module, only : tag_phi_error
    use read_input_module,  only : surfdist 

    integer, intent(in), value :: lev
    type(c_ptr), intent(in), value :: cp
    real(amrex_real), intent(in), value :: t
    character(kind=c_char), intent(in), value :: settag, cleartag
    
    type(amrex_geometry) :: geom  	! geometry at level
    real(amrex_real), allocatable, save :: phierr(:)
    type(amrex_parmparse) :: pp
    type(amrex_tagboxarray) :: tag
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: phiarr(:,:,:,:)
    character(kind=c_char), contiguous, pointer :: tagarr(:,:,:,:)

    ! ! PART OF THE INPUT IS READ HERE!, REMOVE THIS!
    ! if (.not.allocated(phierr)) then
    !    call amrex_parmparse_build(pp, "myamr")
    !    call pp%getarr("phierr", phierr)
    !    call amrex_parmparse_destroy(pp)
    ! end if

    tag = cp
    geom = amrex_geom(lev) 

	
    !$omp parallel private(mfi, bx, phiarr, tagarr)
    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.true.)
    do while(mfi%next())
       bx = mfi%tilebox()
       phiarr => phi_new(lev)%dataptr(mfi)
       tagarr => tag%dataptr(mfi)
       ! call tag_phi_error(lev, t, bx%lo, bx%hi, &
       !      geom%get_physical_location(bx%lo), geom%dx, surfdist(lev+1), & 
       !      phiarr, lbound(phiarr), ubound(phiarr), &
       !      tagarr, lbound(tagarr), ubound(tagarr), &
       !      phierr(lev+1), settag, cleartag)  ! +1 because level starts with 0, but phierr starts with 1
       call tag_phi_error(lev, t, bx%lo, bx%hi, &
            geom%get_physical_location(bx%lo), geom%dx, surfdist(lev+1), & 
            phiarr, lbound(phiarr), ubound(phiarr), &
            tagarr, lbound(tagarr), ubound(tagarr), &
            settag, cleartag)
    end do
    call amrex_mfiter_destroy(mfi)
    !$omp end parallel

  end subroutine my_error_estimate
  
end module my_amr_module
