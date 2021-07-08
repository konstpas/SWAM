module timestep_module

  use amrex_amr_module

  implicit none
  private

  public :: timestep 

contains


  recursive subroutine timestep (lev, time, substep)

    use my_amr_module, only : regrid_int, stepno, nsubsteps, dt, do_reflux, verbose
    use amr_data_module, only : t_old, t_new, phi_old, phi_new, temp, flux_reg  
    use averagedown_module, only : averagedownto

    integer, intent(in) :: lev, substep
    real(amrex_real), intent(in) :: time    
    integer, allocatable, save :: last_regrid_step(:)
    integer :: k, old_finest_level, finest_level, fine_substep    
    
    ! Regridding 
    if (regrid_int .gt. 0) then
       
       if (.not.allocated(last_regrid_step)) then

          allocate(last_regrid_step(0:amrex_max_level))
          last_regrid_step = 0
          
       end if

      
       if (lev .lt. amrex_max_level .and. stepno(lev) .gt. last_regrid_step(lev)) then

          if (mod(stepno(lev), regrid_int) .eq. 0) then

             old_finest_level = amrex_get_finest_level()
             call amrex_regrid(lev, time)
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
    
    
    ! Advance solution
    stepno(lev) = stepno(lev)+1
    ! We need to update t_old(lev) and t_new(lev) before advance is called because of fillpatch.
    ! Fillpatch fills multifab for incrementing solution, using data from given level and level below, and interpolates ghost points 
    t_old(lev) = time
    t_new(lev) = time + dt(lev)
    call amrex_multifab_swap(phi_old(lev), phi_new(lev))
    call advance(lev, time, dt(lev))


    ! Propagate solution and synchronize levels
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
       ! so that surface temperature is interpolated into vacuum region and diffuses on the wrong side of free interface boundary 
       call averagedownto(lev)  ! set covered coarse cells to be the average of fine
        
    end if
  
  end subroutine timestep


  subroutine advance(lev, time, dt)

    use my_amr_module, only : do_reflux, verbose
    use amr_data_module, only : phi_new, temp, surf_ind, flux_reg  
    use fillpatch_module, only : fillpatch
    use domain_module, only : get_melt_pos, reset_melt_pos 
    use shallow_water_module, only : increment_SW
    use heat_transfer_module, only: increment_enthalpy
    
    integer, intent(in) :: lev
    real(amrex_real), intent(in) :: time, dt        
    integer, parameter :: ngrow = 1    			! number of ghost points in each spatial direction 
    integer :: ncomp, idim, i,j 	! components, ranges, and indexes
    logical :: nodal(3)   					! logical for flux multifabs 
    type(amrex_multifab) :: phiborder, tempborder 		! multifabs on mfi owned tilebox, with ghost points 
    type(amrex_mfiter) :: mfi					! mfi iterator 
    type(amrex_box) :: bx, tbx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pin, pout, ptempin, ptemp ! input, output pointers 
    type(amrex_multifab) :: fluxes(amrex_spacedim)    
    
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
    call fillpatch(lev, time, phiborder)
    
    ! Propagate SW equations (only at max level)
    if (lev.eq.amrex_max_level) then 
       call increment_SW(time, amrex_geom(lev), dt)
    end if

    ! RE-distribute energy from 'lost and gained' domain points how


    ! Set melt interface position array equal to free interface position array 
    ! Since melt layer may span several tile boxes in y-direction (in mfiterator below), we cannot reset within each loop 
    ! Therefore we reset melt position after solving SW, and before propagating temperature 
    ! Melt position is then found after heat has been propagated 
    call reset_melt_pos() 
   
    

!$omp parallel private(mfi,bx,tbx,pin,pout,ptemp)
    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.true.) ! Tiling splits validbox into several tile boxes 
    								! could be useful depending on parallelization approach 
    do while(mfi%next())

       bx = mfi%tilebox()   
   
       pin     => phiborder%dataptr(mfi)
       pout    => phi_new(lev)%dataptr(mfi)
       ptempin => tempborder%dataptr(mfi)
       ptemp   => temp(lev)%dataptr(mfi)
       
       ! Increment solution at given mfi tilebox 
       call increment_enthalpy(time, bx%lo, bx%hi, &
                               pin, lbound(pin),     ubound(pin),     &
                               pout,    lbound(pout),    ubound(pout),    &
                               ptempin, lbound(ptempin), ubound(ptempin), &
                               ptemp,   lbound(ptemp),   ubound(ptemp),   &
                               amrex_geom(lev), dt)  
		

	! Find melt interface y position 
       if (lev.eq.amrex_max_level) then
          call get_melt_pos(bx%lo, bx%hi,                      	&
                          ptemp, lbound(ptemp), ubound(ptemp), 	&
			  amrex_geom(lev))
	end if  

    end do

    call amrex_mfiter_destroy(mfi)
!$omp end parallel

    call amrex_multifab_destroy(phiborder)
    
 end subroutine advance


end module timestep_module
