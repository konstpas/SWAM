module timestep_module

  use amrex_amr_module

  implicit none
  private

  public :: timestep 

contains


  recursive subroutine timestep (lev, time, substep)
    use my_amr_module, only : regrid_int, stepno, nsubsteps, dt, do_reflux, verbose !! input
    use amr_data_module, only : t_old, t_new, phi_old, phi_new, flux_reg !! transients of the solution 
    use averagedown_module, only : averagedownto
    use fillpatch_module, only : fillpatch
    integer, intent(in) :: lev, substep
    real(amrex_real), intent(in) :: time
    
    integer, allocatable, save :: last_regrid_step(:)
    integer :: k, old_finest_level, finest_level, fine_substep
    integer :: redistribute_ngrow
    
    !integer, intent(in) :: step, nsub
    integer :: step, nsub
    integer, parameter :: ngrow = 1 ! 3 			! number of ghost points in each spatial direction 
    
    integer :: ncomp, idim, lowbound(2), hbound(2), i,j 	! components, ranges, and indexes
    logical :: nodal(3)   					! logical for flux multifabs 
    type(amrex_multifab) :: phiborder 			! multifab on mfi owned tilebox, with ghost points 
    type(amrex_mfiter) :: mfi					! mfi iterator 
    type(amrex_box) :: bx, tbx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pin,pout ! input, output pointers 
    type(amrex_multifab) :: fluxes(amrex_spacedim)    
    
    
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

    call amrex_multifab_build(phiborder, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, ngrow) ! check that ngrow is indeed number of ghost points


    call fillpatch(lev, time, phiborder)

!$omp parallel private(mfi,bx,tbx,pin,pout)
    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.true.)
    do while(mfi%next())
    bx = mfi%tilebox()
   
       pin  => phiborder%dataptr(mfi)
       pout => phi_new(lev)%dataptr(mfi)
       
       ! Increment solution at given mfi tilebox 
	call increment(time, bx%lo, bx%hi, &
		pin, lbound(pin), ubound(pin), &
		pout,lbound(pout),ubound(pout), &
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







  subroutine increment(time, lo, hi, &
  			uin, ui_lo, ui_hi, &
  			uout, uo_lo, uo_hi, & 
  			geom, dt)
  			
    use domain_module 			
    integer, intent(in) :: lo(2), hi(2)  		! bounds of current tile box
    real(amrex_real), intent(in) :: dt, time		! grid size, sub time step, and time 
    type(amrex_geometry), intent(in) :: geom  	! geometry at level
    integer, intent(in) :: ui_lo(2), ui_hi(2)		! bounds of input tilebox 
    integer, intent(in) :: uo_lo(2), uo_hi(2)		! bounds of output tilebox 
    real(amrex_real), intent(in   ) :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2)) ! 
    real(amrex_real), intent(inout) :: uout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2)) ! 
    real(amrex_real) :: uface(ui_lo(1):ui_hi(1)+1,ui_lo(2):ui_hi(2)  ) ! face velocity x direction (nodal)
    real(amrex_real) :: fluxx(ui_lo(1):ui_hi(1)+1,ui_lo(2):ui_hi(2)  ) ! flux x direction (nodal)
    real(amrex_real) :: fluxy(ui_lo(1):ui_hi(1)  ,ui_lo(2):ui_hi(2)+1) ! flux y direction (nodal)
    real(amrex_real) :: qbound(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2))	! Volumetric heating (boundary)
    real(amrex_real) :: qheat(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2))	! Volumetric heating
    real(amrex_real) :: dx(2) 
    logical :: xfluxflag(ui_lo(1):ui_hi(1)+1,ui_lo(2):ui_hi(2)  ) 	! surface flag for x-nodes 
    logical :: yfluxflag(ui_lo(1):ui_hi(1)  ,ui_lo(2):ui_hi(2)+1) 	! surface flag for y-nodes   
    
    ! uface is face velocity in whole domain, extracted from SW equations. 
    
    integer :: i,j 
   
    dx = geom%dx(1:amrex_spacedim) ! grid width at level 
    
  	call surface_tag(time, geom%get_physical_location(ui_lo), dx, ui_lo, ui_hi, xfluxflag, yfluxflag)
  	call get_face_velocity(time, geom%get_physical_location(ui_lo), dx, uface, ui_lo, ui_hi ) 	! domain module 
  	call create_face_flux(time, geom%get_physical_location(ui_lo), dx, & 			! domain module 
  				uin, uface, ui_lo, ui_hi, xfluxflag, yfluxflag, & 
  				fluxx, fluxy)
  	! Zero flux across surface boundary. Volumetric heat deposition in first internal cell constitutes absorbed boundary flux. 
  	! Incorporates all absorption and cooling terms 			
  	call get_bound_heat(time, geom%get_physical_location(ui_lo), dx, ui_lo, ui_hi, yfluxflag, qbound) 	! domain module 			
  	!call volume_heating(time, geom%get_physical_location(ui_lo), dx, qheat) 
  	


  	
  	
  	do i = lo(1),hi(1)
  	 do j = lo(2),hi(2) 
  	  uout(i,j) = uin(i,j) &
  	     - dt/dx(1)      * (fluxx(i+1,j  )-fluxx(i,j))	&		! flux divergence x-direction 
  	     - dt/dx(2)      * (fluxy(i  ,j+1)-fluxy(i,j))	& 		! flux divergence y-direction 
  	     + dt*qbound(i,j)							! 'boundary volumetric' source
  	 end do 
  	end do 
  	
  
  
  end subroutine increment 





end module timestep_module
