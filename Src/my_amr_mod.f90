module my_amr_module

  use iso_c_binding
  use amrex_amr_module
  use amrex_fort_module, only : rt => amrex_real

  use amr_data_module
  implicit none

  ! runtime parameters
  integer :: verbose    = 0

  integer :: max_step   = huge(1)
  integer :: regrid_int = 2
  integer :: check_int  = -1
  integer :: plot_int   = -1
  integer :: max_grid_size_2d   = 16
  !
  logical :: do_reflux  = .true.
  !
  real(rt) :: stop_time  = huge(1._rt)
  real(rt) :: cfl        = 0.7_rt
  !
  character(len=:), allocatable, save :: check_file
  character(len=:), allocatable, save :: plot_file
  character(len=:), allocatable, save :: restart

  integer, allocatable, save :: stepno(:)
  integer, allocatable, save :: nsubsteps(:)

  real(rt), allocatable, save :: dt(:)

  integer, private, parameter :: ncomp = 1, nghost = 0
  
  
contains

  subroutine my_amr_init ()
    use bc_module, only : lo_bc, hi_bc
    use domain_module 
    use material_properties_module, only : init_mat_prop, material 
    type(amrex_parmparse) :: pp
    integer :: ilev
    
    if (.not.amrex_amrcore_initialized()) call amrex_amrcore_init()
    
    call amrex_init_virtual_functions (my_make_new_level_from_scratch, &
         &                             my_make_new_level_from_coarse,  &
         &                             my_remake_level,                &
         &                             my_clear_level,                 &
         &                             my_error_estimate)

    ! some default parameters
    allocate(character(len=3)::check_file)
    check_file = "chk"
    allocate(character(len=3)::plot_file)
    plot_file = "plt"
    allocate(character(len=0)::restart)

    ! Read parameters
    call amrex_parmparse_build(pp)
    call pp%query("max_step", max_step)
    call pp%query("stop_time", stop_time)
    call amrex_parmparse_destroy(pp)
    
    ! Parameters amr.*
    call amrex_parmparse_build(pp, "amr")
    call pp%query("regrid_int", regrid_int)
    call pp%query("check_int", check_int)
    call pp%query("plot_int", plot_int)
    call pp%query("check_file", check_file)
    call pp%query("plot_file", plot_file)
    call pp%query("restart", restart)
    call pp%query("max_grid_size_2d", max_grid_size_2d)
    call amrex_parmparse_destroy(pp)
    
    ! Domain parameters 
    call amrex_parmparse_build(pp, "domain")
    call pp%query("meltvel", meltvel)  
    call pp%query("tempinit", tempinit) 
    call pp%get("surf_pos", surf_pos_init)   
    call pp%getarr("surfdist", surfdist)  
    call pp%get("material", material)  
    call pp%query("flux_peak", flux_peak) 
    call pp%getarr("flux_pos", flux_pos)
    call pp%getarr("flux_width", flux_width)
    call pp%query("exp_time", exp_time) 
    call amrex_parmparse_destroy(pp)
    call init_mat_prop ! Initialize material properties 
    
    
    ! Parameters myamr.*
    call amrex_parmparse_build(pp, "myamr")
    call pp%query("v", verbose)
    call pp%query("verbose", verbose)
    call pp%query("cfl", cfl)
    call pp%query("do_reflux", do_reflux)
    call amrex_parmparse_destroy(pp)

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
	
	surf_dx(1) = amrex_geom(0)%dx(1)/(2**amrex_max_level)
	surf_dx(2) = amrex_geom(0)%dx(3)/(2**amrex_max_level)
	
	
	
    call amr_data_init()
  end subroutine my_amr_init


  subroutine my_amr_finalize ()
    call amr_data_finalize()
  end subroutine my_amr_finalize

  ! Make a new level from scratch and put the data in phi_new.
  ! Note that phi_old contains no valid data after this.
  subroutine my_make_new_level_from_scratch (lev, time, pba, pdm) bind(c)
    use initdomain_module, only : init_phi
    use domain_module, only : tempinit, surf_pos_init
    use material_properties_module, only : get_temp 
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm

    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: phi(:,:,:,:), ptemp(:,:,:,:)
    
    ! Transfer ba, dm, mf, bx 
    type(amrex_boxarray) :: trba
    type(amrex_distromap) :: trdm
    type(amrex_box) :: trbx
    logical :: firstbox = .true. 

    ba = pba
    dm = pdm

    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real

    call my_clear_level(lev)
    call amrex_multifab_destroy(surface) 
    
    call amrex_multifab_build(phi_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_old(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(temp(lev), ba, dm, ncomp, nghost)
    ! Build multifab for surface position for ba where each box has only one y-component 
    
    

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
       	  	! index space covered at highest level 
    			if (lev.eq.amrex_max_level) then 
       	  	 if (firstbox) then 
  			surf_ind(1,1) = bx%lo(1)
  			surf_ind(2,1) = bx%lo(3)
  			surf_ind(1,2) = bx%hi(1)
  			surf_ind(2,2) = bx%hi(3)
  			firstbox = .false. 
    			 else 
  			surf_ind(1,1) = min(surf_ind(1,1),bx%lo(1))
  			surf_ind(2,1) = min(surf_ind(2,1),bx%lo(3))
  			surf_ind(1,2) = max(surf_ind(1,2),bx%hi(1))
  			surf_ind(2,2) = max(surf_ind(2,2),bx%hi(3))
    			 end if 
    			end if 
    end do


    call amrex_mfiter_destroy(mfi)
    
    if (lev.eq.amrex_max_level) then 



		    trbx%lo(1) = surf_ind(1,1)
		    trbx%lo(2) = 0 
		    trbx%lo(3) = surf_ind(2,1)
		    trbx%hi(1) = surf_ind(1,2)
		    trbx%hi(2) = 0 
		    trbx%hi(3) = surf_ind(2,2)
			
			! Intermediate multifab
			call amrex_boxarray_build(trba,trbx)
			
			call trba%maxSize(max_grid_size_2d)
						
			call amrex_distromap_build(trdm,trba)
			call amrex_multifab_build(surface, trba, trdm, 1, 0 )
			call surface%setval(surf_pos_init)

			call amrex_distromap_destroy(trdm) 
			call amrex_boxarray_destroy(trba) 

    end if 

  end subroutine my_make_new_level_from_scratch

  ! Make a new level from coarse level and put the data in phi_new.
  ! Note tha phi_old contains no valid data after this.
  subroutine my_make_new_level_from_coarse (lev, time, pba, pdm) bind(c)
    use fillpatch_module, only : fillcoarsepatch
    use domain_module, only : surf_pos_init
    use material_properties_module, only : get_temp 
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm

    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: phi(:,:,:,:), ptemp(:,:,:,:)

    type(amrex_multifab) :: trmf 
    type(amrex_boxarray) :: trba
    type(amrex_distromap) :: trdm
    type(amrex_box) :: trbx
    integer :: tr_ind(2,2)=0 ! transfer fluid domain index bounds 


    ba = pba
    dm = pdm

    call my_clear_level(lev)

    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real

    call amrex_multifab_build(phi_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_old(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(temp(lev), ba, dm, ncomp, nghost)
    ! Build surface multifab 
    
    if (lev > 0 .and. do_reflux) then
       call amrex_fluxregister_build(flux_reg(lev), ba, dm, amrex_ref_ratio(lev-1), lev, ncomp)
    end if

     
    call fillcoarsepatch(lev, time, phi_new(lev))
    ! Interpolate surface data 
    
    ! Fill temperature data 
    call amrex_mfiter_build(mfi, temp(lev))
    
    do while (mfi%next())
       bx = mfi%tilebox()
       phi => phi_new(lev)%dataptr(mfi)
       ptemp => temp(lev)%dataptr(mfi)
       call get_temp(bx%lo, bx%hi, & 
       	lbound(phi),   ubound(phi),   phi, &
       	lbound(ptemp), ubound(ptemp), ptemp)  
       	
    			if (lev.eq.amrex_max_level) then 
  			tr_ind(1,1) = min(surf_ind(1,1),bx%lo(1))
  			tr_ind(2,1) = min(surf_ind(2,1),bx%lo(3))
  			tr_ind(1,2) = max(surf_ind(1,2),bx%hi(1))
  			tr_ind(2,2) = max(surf_ind(2,2),bx%hi(3))
    			end if  
    end do    
    
    call amrex_mfiter_destroy(mfi)
   
   
	! 1. define transfer surface multifab defined on low-hi box bounds in x,z coordinates at highest level 
	! 2. Get data based on existing surf multifab (new multifab has different boxarray) (as by fillpatch) 
	! 3. destroy current surface multifab 
	! 4. make new surface multifab on level (with updated boxarray) 
	! 5. Transfer data from transfer to surface multifab     
if (lev.eq.amrex_max_level) then 

		if ((tr_ind(1,1).lt.surf_ind(1,1)).or.&
		    (tr_ind(2,1).lt.surf_ind(2,1)).or.&
		    (tr_ind(1,2).gt.surf_ind(1,2)).or.&
		    (tr_ind(2,2).gt.surf_ind(2,2))) then 
		    
		    ! New index bounds of surface domain 
		    trbx%lo(1) = min(tr_ind(1,1),surf_ind(1,1))
		    trbx%lo(2) = 0 
		    trbx%lo(3) = min(tr_ind(2,1),surf_ind(2,1))
		    trbx%hi(1) = max(tr_ind(1,2),surf_ind(1,2))
		    trbx%hi(2) = 0 
		    trbx%hi(3) = max(tr_ind(2,2),surf_ind(2,2))

			
			! Intermediate multifab
			call amrex_boxarray_build(trba,trbx)
			call trba%maxSize(max_grid_size_2d)		
			call amrex_distromap_build(trdm,trba)
			call amrex_multifab_build(trmf, trba, trdm, 1, 0 )
			call trmf%setval(surf_pos_init)
			
			! copy data from surface to transfermf 
			! Parallel_copy performs on regions of intersection (here, intersection is previous surface domain)
			call trmf%parallel_copy(surface,amrex_geom(lev))
			
			! Remake surface multifab on new region 
			call amrex_multifab_destroy(surface) 
			call amrex_multifab_build(surface, trba, trdm, 1, 0 )
			
			! swap transfer and surface multifabs
			call amrex_multifab_swap(trmf, surface)
			
			call amrex_multifab_destroy(trmf) 
			call amrex_distromap_destroy(trdm) 
			call amrex_boxarray_destroy(trba) 

		end if
    
    
    end if 
    
  end subroutine my_make_new_level_from_coarse

  ! Remake a level from current and coarse levels and put the data in phi_new.
  ! Note that phi_old contains no valid data after this.
  subroutine my_remake_level (lev, time, pba, pdm) bind(c)
    use fillpatch_module, only : fillpatch
    use domain_module, only : surf_pos_init
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

    type(amrex_multifab) :: trmf 
    type(amrex_boxarray) :: trba
    type(amrex_distromap) :: trdm
    type(amrex_box) :: trbx
    integer :: tr_ind(2,2)=0 ! transfer fluid domain index bounds 
    

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
    			if (lev.eq.amrex_max_level) then 
  			tr_ind(1,1) = min(surf_ind(1,1),bx%lo(1))
  			tr_ind(2,1) = min(surf_ind(2,1),bx%lo(3))
  			tr_ind(1,2) = max(surf_ind(1,2),bx%hi(1))
  			tr_ind(2,2) = max(surf_ind(2,2),bx%hi(3))
    			end if  
     end do    
    call amrex_mfiter_destroy(mfi)


	! 1. define transfer surface multifab defined on low-hi box bounds in x,z coordinates at highest level 
	! 2. Get data based on existing surf multifab (new multifab has different boxarray) (as by fillpatch) 
	! 3. destroy current surface multifab 
	! 4. make new surface multifab on level (with updated boxarray) 
	! 5. Transfer data from transfer to surface multifab 
if (lev.eq.amrex_max_level) then 
		if ((tr_ind(1,1).lt.surf_ind(1,1)).or.&
		    (tr_ind(2,1).lt.surf_ind(2,1)).or.&
		    (tr_ind(1,2).gt.surf_ind(1,2)).or.&
		    (tr_ind(2,2).gt.surf_ind(2,2))) then 
		    
		    ! New index bounds of surface domain 
		    trbx%lo(1) = min(tr_ind(1,1),surf_ind(1,1))
		    trbx%lo(2) = 0 
		    trbx%lo(3) = min(tr_ind(2,1),surf_ind(2,1))
		    trbx%hi(1) = max(tr_ind(1,2),surf_ind(1,2))
		    trbx%hi(2) = 0 
		    trbx%hi(3) = max(tr_ind(2,2),surf_ind(2,2))

			
			! Intermediate multifab
			call amrex_boxarray_build(trba,trbx)
			call trba%maxSize(max_grid_size_2d)		
			call amrex_distromap_build(trdm,trba)
			call amrex_multifab_build(trmf, trba, trdm, 1, 0 )
			call trmf%setval(surf_pos_init)
			
			! copy data from surface to transfermf 
			! Parallel_copy performs on regions of intersection (here, intersection is previous surface domain)
			call trmf%parallel_copy(surface,amrex_geom(lev))
			
			! Remake surface multifab on new region 
			call amrex_multifab_destroy(surface) 
			call amrex_multifab_build(surface, trba, trdm, 1, 0 )
			
			! swap transfer and surface multifabs
			call amrex_multifab_swap(trmf, surface)
			
			call amrex_multifab_destroy(trmf) 
			call amrex_distromap_destroy(trdm) 
			call amrex_boxarray_destroy(trba) 

		end if
    
    
    end if 
    
  end subroutine my_remake_level

  subroutine my_clear_level (lev) bind(c)
    integer, intent(in), value :: lev
    call amrex_multifab_destroy(phi_new(lev))
    call amrex_multifab_destroy(phi_old(lev))
    call amrex_multifab_destroy(temp(lev))
    ! destroy surface multifabs 
    call amrex_fluxregister_destroy(flux_reg(lev))
  end subroutine my_clear_level

  subroutine my_error_estimate (lev, cp, t, settag, cleartag) bind(c)
    use tagging_module, only : tag_phi_error
    use domain_module,  only : surfdist 
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

    if (.not.allocated(phierr)) then
       call amrex_parmparse_build(pp, "myamr")
       call pp%getarr("phierr", phierr)
       call amrex_parmparse_destroy(pp)
    end if

    tag = cp
    geom = amrex_geom(lev) 
    
	
    !$omp parallel private(mfi, bx, phiarr, tagarr)
    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.true.)
    do while(mfi%next())
       bx = mfi%tilebox()
       phiarr => phi_new(lev)%dataptr(mfi)
       tagarr => tag%dataptr(mfi)
       call tag_phi_error(lev, t, bx%lo, bx%hi, &
            geom%get_physical_location(bx%lo), geom%dx, surfdist(lev+1), & 
            phiarr, lbound(phiarr), ubound(phiarr), &
            tagarr, lbound(tagarr), ubound(tagarr), &
            phierr(lev+1), settag, cleartag)  ! +1 because level starts with 0, but phierr starts with 1
    end do
    call amrex_mfiter_destroy(mfi)
    !$omp end parallel

  end subroutine my_error_estimate
  
end module my_amr_module
