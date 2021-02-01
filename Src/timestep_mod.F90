module timestep_module

  use amrex_amr_module

  implicit none
  private

  public :: timestep 

contains


  recursive subroutine timestep (lev, time, substep)
    use my_amr_module, only : regrid_int, stepno, nsubsteps, dt, do_reflux, verbose !! input
    use amr_data_module, only : t_old, t_new, phi_old, phi_new, temp, surface, surf_ind, flux_reg, loop_test ! (remove loop_test) !! transients of the solution 
    use averagedown_module, only : averagedownto
    use fillpatch_module, only : fillpatch
    integer, intent(in) :: lev, substep
    real(amrex_real), intent(in) :: time
    
    integer, allocatable, save :: last_regrid_step(:)
    integer :: k, old_finest_level, finest_level, fine_substep
    integer :: redistribute_ngrow
    
    !integer, intent(in) :: step, nsub
    integer :: step, nsub
    integer, parameter :: ngrow = 1    			! number of ghost points in each spatial direction 
    
    integer :: ncomp, idim, lowbound(2), hbound(2), i,j 	! components, ranges, and indexes
    logical :: nodal(3)   					! logical for flux multifabs 
    type(amrex_multifab) :: phiborder, tempborder 		! multifabs on mfi owned tilebox, with ghost points 
    type(amrex_multifab) :: surface_pos			! Surface position multifab test 
    
    type(amrex_mfiter) :: mfi					! mfi iterator 
    type(amrex_box) :: bx, tbx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pin, pout, ptempin, ptemp ! input, output pointers 
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: psurf, ptransient  
    
    
    type(amrex_multifab) :: fluxes(amrex_spacedim)    
    
    type(amrex_multifab) :: testmf, transientmf 
    type(amrex_box)      :: testbox 
    type(amrex_boxarray) :: testba
    type(amrex_distromap):: testdm 
    real(amrex_real) :: mfval = 2.0, mfval2 = 1.0  
    integer :: ba_maxSize_test  
    
    ! Currently face velocity and face fluxes are declared as arrays in the increment subroutine 
    ! Could be declared here as amrex_fab or amreX_multifab, could make it more neat with indexing, using nodalize
    ! and all components of a velocity field contained in single object. 
    ! Also necessary if output of those quantities is wanted on the amrex standard form.  
    
    
    
    ! Regridding 
    !---------------------------------------------------------
    if (regrid_int .gt. 0) then
       if (.not.allocated(last_regrid_step)) then
          allocate(last_regrid_step(0:amrex_max_level))
          last_regrid_step = 0
       end if

       ! regrid doesn't change the base level, so we don't regrid on amrex_max_level
       if (lev .lt. amrex_max_level .and. stepno(lev) .gt. last_regrid_step(lev)) then
          if (mod(stepno(lev), regrid_int) .eq. 0) then

             old_finest_level = amrex_get_finest_level()
             call amrex_regrid(lev, time) ! the finest level may change during regrid
             finest_level = amrex_get_finest_level()

             do k = lev, finest_level
                last_regrid_step(k) = stepno(k)
             end do

             do k = old_finest_level+1, finest_level
                dt(k) = dt(k-1) / amrex_ref_ratio(k-1)
             end do
             
          end if
       end if
    end if
    !---------------------------------------------------------
    
    
    ! Prepare to advance solution 
    stepno(lev) = stepno(lev)+1
    ! We need to update t_old(lev) and t_new(lev) before advance is called because of fillpatch.
    ! Fillpatch fills multifab for incrementing solution, using data from given level and level below, and interpolates ghost points 
    t_old(lev) = time
    t_new(lev) = time + dt(lev)
    ! swap phi_new(lev) and phi_old(lev) so they are consistent with t_new(lev) and t_old(lev)
    call amrex_multifab_swap(phi_old(lev), phi_new(lev))


    ! Increment solution at current level and substep 
    !--------------------------------------------------------

    ncomp = phi_new(lev)%ncomp()

    if (do_reflux) then
       do idim = 1, amrex_spacedim
          nodal = .false.
          nodal(idim) = .true.
          !call amrex_multifab_build(fluxes(idim), phi_new(lev)%ba, phi_new(lev)%dm, ncomp, 0, nodal)
       end do
    end if

    call amrex_multifab_build(phiborder, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, ngrow) 
    call amrex_multifab_build(tempborder, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, ngrow)
    !call amrex_multifab_build()
    
    !testbox = phi_new(lev)%ba%get_box!(1) 
    !!print *, phi_new(lev)%ba%num_pts()
    !print *, phi_new(lev)%ba%get_box()
    !print *, 'testbox lo hi is ', testbox%lo, testbox%hi 
    !print *, '' 
    !call amrex_print(phi_new(lev)%ba)


    call fillpatch(lev, time, phiborder)



		! Procedure as in my_remake_level 
	! 1. define transient surface multifab defined on low-hi box bounds in x,z coordinates at highest level 
	! 2. Get data based on existing surf multifab (new multifab has different boxarray) (as by fillpatch) 
	! 3. destroy current surface multifab 
	! 4. make new surface multifab on level (with updated boxarray) 
	! 5. Transfer data from transient to surface multifab 

	! Then how to store surface data for removed points ?
	! Perhaps never remove points, only add. 
	! In that case, keep old boxarray and proceed past step 1 only if bigger index space is encompassed. 
	! 

!--------------Testing-----------------------------------------
	! Multifab 'surface' defined in amr_data_mod 
	! 
	!print *, surface%ba%get_box(0)  !phi_new(lev)%ba%get_box!(1)
	!pause 
	
	testbox%lo(1) = 6 
	testbox%lo(2) = 0  
	testbox%lo(3) = 6 - loop_test 
	
	testbox%hi(1) = 10! + loop_test
	testbox%hi(2) = 0  
	testbox%hi(3) = 10 + loop_test 
	
	call amrex_boxarray_build(testba,testbox)
	ba_maxSize_test = 4 
	call testba%maxSize(ba_maxSize_test)
	!testba%maxSize(2)
	
	call amrex_distromap_build(testdm,testba) 
	
	! print *, 'testbox indexes are lo, hi ', testbox%lo, testbox%hi 
	!call amrex_print(testba)
	print *, 'loop is ', loop_test 
	if (loop_test.eq.1) then ! second loop, test multifab expansion idea  
		call amrex_multifab_build(surface, testba, testdm, 1, 0 )
		surf_ind(1,1) = testbox%lo(1)
		surf_ind(2,1) = testbox%lo(3)
		surf_ind(1,2) = testbox%hi(1)
		surf_ind(2,2) = testbox%hi(3)
		call surface%setval(mfval)  
	else 
		! 1. Find if new ba encompasses larger index space 
		! 2. Create transient mf 
		! 3. Transfer data from surface to transient mf 
		! 4. Remake surface multifab
		! 5. Transfer data back (multifab swap)
		
		! if transient box not contained within surface mf ba
		if ((testbox%lo(1).lt.surf_ind(1,1)).or.&
		    (testbox%lo(3).lt.surf_ind(2,1)).or.&
		    (testbox%hi(1).gt.surf_ind(1,2)).or.&
		    (testbox%hi(3).gt.surf_ind(2,2))) then 
		    
		    ! New bounds 
		    testbox%lo(1) = min(testbox%lo(1),surf_ind(1,1))
		    testbox%lo(3) = min(testbox%lo(3),surf_ind(2,1))
		    testbox%hi(1) = max(testbox%hi(1),surf_ind(1,2))
		    testbox%hi(3) = max(testbox%hi(3),surf_ind(2,2))
		  	surf_ind(1,1) = testbox%lo(1)
			surf_ind(2,1) = testbox%lo(3)
			surf_ind(1,2) = testbox%hi(1)
			surf_ind(2,2) = testbox%hi(3)
			
			! Intermediate multifab
			call amrex_boxarray_build(testba,testbox)
			ba_maxSize_test = 4 
			call testba%maxSize(ba_maxSize_test)
			call amrex_distromap_build(testdm,testba)
			call amrex_multifab_build(transientmf, testba, testdm, 1, 0 )
			call transientmf%setval(mfval2)
			
			! copy data from surface to transientmf 
			! Parallel_copy performs on regions of intersection (here, intersection is previous surface domain)
			call transientmf%parallel_copy(surface,amrex_geom(lev))
			
			! Remake surface multifab on new region 
			call amrex_multifab_destroy(surface) 
			call amrex_multifab_build(surface, testba, testdm, 1, 0 )
			
			! swap transient and surface multifabs
			call amrex_multifab_swap(transientmf, surface)
			

		end if 
		
	end if 
	
	!call amrex_print(testba)
	print *, 'sum of surface is', surface%sum() 
	
    
    
    loop_test = loop_test +1 ! for testing only 
    pause 	
! -------------------------------------------------------------







!$omp parallel private(mfi,bx,tbx,pin,pout,ptemp)
    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.true.) ! Tiling splits validbox into several tile boxes 
    								! could be useful depending on parallelization approach 
    do while(mfi%next())
    bx = mfi%tilebox()   
   
    ! if lev eq amrex_max level (2d fluid domain corresponds to )
    ! find lowest and highest box indexes in x,z plane.  
    ! compare to stored indexes of surface multifab 
   
       pin     => phiborder%dataptr(mfi)
       pout    => phi_new(lev)%dataptr(mfi)
       ptempin => tempborder%dataptr(mfi)
       ptemp   => temp(lev)%dataptr(mfi)
       
       ! Increment solution at given mfi tilebox 
	call increment(time, bx%lo, bx%hi, &
		pin,     lbound(pin),     ubound(pin),     &
		pout,    lbound(pout),    ubound(pout),    &
		ptempin, lbound(ptempin), ubound(ptempin), &
		ptemp,   lbound(ptemp),   ubound(ptemp),   &
		amrex_geom(lev), dt(lev))  

    end do ! while(mfi%next())
    call amrex_mfiter_destroy(mfi)
!$omp end parallel

    call amrex_multifab_destroy(phiborder)


!---------------------------------------------------------

    ! Recursive call to propagate every substep of every level and make the appropriate 
    ! level synchronization via averaging and flux divergence 
    if (lev .lt. amrex_get_finest_level()) then
       do fine_substep = 1, nsubsteps(lev+1)
          call timestep(lev+1, time+(fine_substep-1)*dt(lev+1), fine_substep) 
       end do

       ! Update coarse solution at coarse-fine interface via dPhidt = -div(+F) where F is stored flux (flux_reg)
       if (do_reflux) then
          !call flux_reg(lev+1)%reflux(phi_new(lev), 1.0_amrex_real) 
       end if

       ! A problem occurs when mesh size in region containing free interface is changed 
       ! Finer mesh resolves interface somewhere withing coarse grid point
       call averagedownto(lev)  ! set covered coarse cells to be the average of fine
        
    end if
  
  end subroutine timestep





		! Subroutine takes input enthalpy and indexes (with ghost points interpolated from nearby multifab at level and below) 
		! Updates enthalpy by divergence of enthalpy flux 

  subroutine increment(time, lo, hi, &
  			uin,    ui_lo, ui_hi, &
  			uout,   uo_lo, uo_hi, & 
  			tempin, ti_lo, ti_hi, & 
  			temp,   t_lo , t_hi , & 
  			geom, dt)
  			
    use domain_module 	
    use material_properties_module, only : get_temp  		
    integer, intent(in) :: lo(3), hi(3)  		! bounds of current tile box
    real(amrex_real), intent(in) :: dt, time		! grid size, sub time step, and time 
    type(amrex_geometry), intent(in) :: geom  	! geometry at level
    integer, intent(in) :: ui_lo(3), ui_hi(3)		! bounds of input enthalpy box 
    integer, intent(in) :: uo_lo(3), uo_hi(3)		! bounds of output enthalpy box  
    integer, intent(in) :: ti_lo(3), ti_hi(3)		! bounds of temperature box  
    integer, intent(in) :: t_lo (3), t_hi (3)		! bounds of temperature box    
    real(amrex_real) :: dx(3) 			! Grid size 

    real(amrex_real), intent(in   ) :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! Enthalpy in 
    real(amrex_real), intent(inout) :: uout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3)) ! Enthalpy out 
    real(amrex_real), intent(inout) :: tempin (ti_lo(1):ti_hi(1),ti_lo(2):ti_hi(2),ti_lo(3):ti_hi(3)) ! Temperature on enthalpy in-box
    real(amrex_real), intent(inout) :: temp   (t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) ! Temperature out 

    
    real(amrex_real) :: uface (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! face velocity x direction (nodal)
    real(amrex_real) :: fluxx (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! flux x direction (nodal)
    real(amrex_real) :: fluxy (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! flux y direction (nodal)
    
    real(amrex_real) :: qbound(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))	! Volumetric heating (boundary)
    real(amrex_real) :: qheat (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))	! Volumetric heating
    
    logical :: xfluxflag(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! surface flag for x-nodes 
    logical :: yfluxflag(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! surface flag for y-nodes   
    integer :: i,j,k  
#if AMREX_SPACEDIM == 3 
    real(amrex_real) :: wface (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! face velocity z direction (nodal)	 
    real(amrex_real) :: fluxz (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! flux z direction (nodal)
    logical :: zfluxflag(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! surface flag for z-nodes
#endif        

   
    dx = geom%dx(1:3) ! grid width at level 
    
    

   	call get_temp(ti_lo, ti_hi,             &	! tilebox indexes 
  		      ui_lo, ui_hi, uin, &	! Output enthalpy indexes and data array
  		      ti_lo, ti_hi, tempin)	! Temperature indexes and data array    
    
    
    
        ! Subroutine assigns logical arrays denoting free interface boundary 
  	call surface_tag(time, geom%get_physical_location(ui_lo), dx, lo, hi, &
  			xfluxflag, yfluxflag, ui_lo, ui_hi &
#if AMREX_SPACEDIM == 3 
  			, zfluxflag & 
#endif  
			)	
  	
  	! Subroutine assigns tangential (x,z) velocity on grid edges in whole domain  
  	call get_face_velocity(time, geom%get_physical_location(lo), dx, lo, hi, uface, &
#if AMREX_SPACEDIM == 3  
				wface, &
#endif 
  				ui_lo, ui_hi ) 	
 
 		
  				
  				 	

  	! Subroutine assigns enthalpy flux on grid edges in whole domain 
  	call create_face_flux(time, geom%get_physical_location(ui_lo), dx, lo, hi, & 			! domain module 
  				uin, uface, xfluxflag, yfluxflag, 	& 
  				fluxx, fluxy,  			& 
#if AMREX_SPACEDIM == 3 
				wface, zfluxflag, fluxz, 		&
#endif 	
				        ui_lo, ui_hi, 			&
  				tempin, ti_lo, ti_hi)
  	
  	
  	
  	! Zero flux across surface boundary. Volumetric heat deposition in first internal cell constitutes absorbed boundary flux. 
  	! Incorporates all absorption and cooling terms 			
  	call get_bound_heat(time, geom%get_physical_location(lo), dx, lo, hi, ui_lo, ui_hi, yfluxflag, qbound) 	! domain module 			
  	!call volume_heating(time, geom%get_physical_location(ui_lo), dx, qheat) 
  	


  	
  	! Increment enthalpy 
  	do   i = lo(1),hi(1)
  	 do  j = lo(2),hi(2) 
  	  do k = lo(3),hi(3)
  	  uout(i,j,k) = uin(i,j,k) &
  	     - dt/dx(1)      * (fluxx(i+1,j  ,k  )-fluxx(i,j,k))	&		! flux divergence x-direction 
  	     - dt/dx(2)      * (fluxy(i  ,j+1,k  )-fluxy(i,j,k))	& 		! flux divergence y-direction 
#if AMREX_SPACEDIM == 3 
  	     - dt/dx(3)      * (fluxz(i  ,j  ,k+1)-fluxz(i,j,k))	&		! flux divergence z-direction
#endif 
  	     + dt*qbound(i,j,k)	
  	  end do    						! 'boundary volumetric' source
  	 end do 
  	end do 
  	
	
  	! Find temperature 
  	call get_temp(lo, hi,             &	! tilebox indexes 
  		      uo_lo, uo_hi, uout, &	! Output enthalpy indexes and data array
  		      t_lo , t_hi , temp)	! Temperature indexes and data array 
 
  	
  	
  
  end subroutine increment 







end module timestep_module
