module simulation_module

  ! -----------------------------------------------------------------
  ! This module is used to run the simulation, i.e. to coordinate
  ! all the necessary calculations involved in the heat transfer,
  ! phase change and melt motion
  ! -----------------------------------------------------------------
  
  use amrex_amr_module

  implicit none
  private

  public :: run_simulation

contains

  ! run the simulation
  subroutine run_simulation()

    use amrex_amr_module 
    use amr_data_module, only : t_new, stepno, dt
    use read_input_module, only : max_step, stop_time, plot_int
    use plotfile_module, only: writeplotfile, write2dplotfile
    use energy_module, only: sum_enthalpy
   
    integer :: last_plot_file_step, step, lev, substep
    real(amrex_real) :: cur_time, total_enthalpy

    ! Initialize time
    cur_time = t_new(0)

    ! Initialize counter for output
    last_plot_file_step = 0;

    ! Output initial configuration
    if (plot_int .gt. 0) call writeplotfile

    ! Loop over time-steps
    do step = stepno(0), max_step-1

       ! Print timestep information on screen (start)
       if (amrex_parallel_ioprocessor()) then
          print *, ""
          print *, "STEP", step+1, "starts ..."
       end if

       ! Compute time step size
       call compute_dt()    
		 
       ! Advance all levels of one time step  
       lev = 0
       substep = 1
       call advance_one_timestep(lev, cur_time, substep)

       ! Get total enthalpy in the domain
       call sum_enthalpy(total_enthalpy)
 	
       ! Update time on all levels
       cur_time = cur_time + dt(0)
       do lev = 0, amrex_get_finest_level()
          t_new(lev) = cur_time
       end do

       ! Print timestep information on screen (end)
       if (amrex_parallel_ioprocessor()) then
          print *, 'Enthalpy is', total_enthalpy 
          print *, '' 
          print *, "STEP", step+1, "end. TIME =", cur_time, "DT =", dt(0)
       end if

       ! Print output to file
       if (plot_int .gt. 0 .and. mod(step+1,plot_int) .eq. 0) then
          last_plot_file_step = step+1
          call writeplotfile
          call write2dplotfile
       end if

       ! Stopping criteria
       if (cur_time .ge. stop_time) exit

    end do

    ! NOT SURE WHY THIS PART IS HERE
    ! if (plot_int .gt. 0 .and. stepno(0) .gt. last_plot_file_step) then
    !    call writeplotfile()
    !    call write2dplotfile() 
    ! end if
    
  end subroutine run_simulation


    subroutine compute_dt()
    
    use amr_data_module, only : t_new, dt, nsubsteps
    use read_input_module, only : stop_time, dt_change_max
    
    integer :: lev, nlevs, n_factor
    real(amrex_real) :: dt_0, eps
    real(amrex_real), allocatable :: dt_tmp(:)

    nlevs = amrex_get_numlevels()
    
    allocate(dt_tmp(0:nlevs-1))

    ! Estimate timestep in order to avoid instabilities
    do lev = 0, nlevs-1
       call est_timestep(lev, dt_tmp(lev))
    end do
    call amrex_parallel_reduce_min(dt_tmp, nlevs)

    ! Ensure that the timestep does not change too much
    ! from one timestep to the next
    dt_0 = dt_tmp(0)
    n_factor = 1
    do lev = 0, nlevs-1
       dt_tmp(lev) = min(dt_tmp(lev), dt_change_max*dt(lev))
       n_factor = n_factor * nsubsteps(lev)
       dt_0 = min(dt_0, n_factor*dt_tmp(lev))
    end do

    ! Compute timestep for all levels
    dt(0) = dt_0
    do lev = 1, nlevs-1
       dt(lev) = dt(lev-1) / nsubsteps(lev)
    end do

    
  end subroutine compute_dt


  subroutine est_timestep(lev, dt)

    use read_input_module, only : cfl
    use material_properties_module, only : max_diffus 

    ! Input and output variables 
    integer, intent(in) :: lev
    real(amrex_real), intent(out) :: dt

    ! Local variables 
    real(amrex_real) :: dxsqr

    ! NOTE: There are two stability criteria to take into
    ! account, the von Neumann stability and the CFL
    ! condition. The CFL is not yet implemented
    
    ! Von Neumann stability criterion 
    dxsqr= (1/amrex_geom(lev)%dx(1)**2 + &
            1/amrex_geom(lev)%dx(2)**2 + &
            1/amrex_geom(lev)%dx(3)**2) 
    dt = 0.5/(dxsqr*max_diffus)
    dt = dt * cfl

  end subroutine est_timestep

  

  ! Subroutine used to advance of one time step. Note that this subroutine is
  ! recursive, which implies that it calls itself
  recursive subroutine advance_one_timestep(lev, time, substep)

    use read_input_module, only : regrid_int, do_reflux, verbose
    use amr_data_module, only : t_old, t_new, phi_old, phi_new, temp, flux_reg, stepno, nsubsteps, dt  
    use regrid_module, only : averagedownto

    ! Input and output variables
    integer, intent(in) :: lev, substep
    real(amrex_real), intent(in) :: time

    ! Local variables
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
    call advance(lev, time, dt(lev), substep)


    ! Propagate solution and synchronize levels
    if (lev .lt. amrex_get_finest_level()) then
       
       do fine_substep = 1, nsubsteps(lev+1)
          call advance_one_timestep(lev+1, time+(fine_substep-1)*dt(lev+1), fine_substep) 
       end do

       ! Update coarse solution at coarse-fine interface via dPhidt = -div(+F) where F is stored flux (flux_reg)
       if (do_reflux) then
          call flux_reg(lev+1)%reflux(phi_new(lev), 1.0_amrex_real) 
       end if
       
       ! A problem occurs when mesh size in region containing free interface is changed 
       ! Finer mesh resolves interface somewhere withing coarse grid point
       ! so that surface temperature is interpolated into vacuum region and diffuses on the wrong side of free interface boundary 
       call averagedownto(lev)  ! set covered coarse cells to be the average of fine
        
    end if
  
  end subroutine advance_one_timestep


  subroutine advance(lev, time, dt, substep)

    use read_input_module, only : do_reflux, verbose
    use amr_data_module, only : phi_new, temp, surf_ind, flux_reg  
    use regrid_module, only : fillpatch
    use heat_transfer_module, only : get_melt_pos, reset_melt_pos 
    use shallow_water_module, only : increment_SW
    use heat_transfer_module, only: increment_enthalpy
    
    integer, intent(in) :: lev, substep
    real(amrex_real), intent(in) :: time, dt        
    integer, parameter :: ngrow = 1    			! number of ghost points in each spatial direction 
    integer :: ncomp, idim, i,j 	! components, ranges, and indexes
    logical :: nodal(3)   					! logical for flux multifabs 
    type(amrex_multifab) :: phiborder, tempborder 		! multifabs on mfi owned tilebox, with ghost points 
    type(amrex_mfiter) :: mfi					! mfi iterator 
    type(amrex_box) :: bx, tbx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pin, pout, ptempin, ptemp, pfx, pfy, pfz, pf, pfab ! input, output pointers
    type(amrex_fab) :: flux(amrex_spacedim)
    type(amrex_multifab) :: fluxes(amrex_spacedim)    
    
    ncomp = phi_new(lev)%ncomp()

    if (do_reflux) then
       do idim = 1, amrex_spacedim
          nodal = .false.
          nodal(idim) = .true.
          call amrex_multifab_build(fluxes(idim), phi_new(lev)%ba, phi_new(lev)%dm, ncomp, 0, nodal)
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
   
    

    !$omp parallel private(mfi,bx,tbx,pin,pout,ptemp,ptempin,pfx,pfy,pfz,pf,pfab,flux)

    do idim = 1, amrex_spacedim
       call flux(idim)%reset_omp_private()
    end do
    
    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.false.) ! Tiling splits validbox into several tile boxes 
    								! could be useful depending on parallelization approach 
    do while(mfi%next())

       bx = mfi%validbox()   
   
       pin     => phiborder%dataptr(mfi)
       pout    => phi_new(lev)%dataptr(mfi)
       ptempin => tempborder%dataptr(mfi)
       ptemp   => temp(lev)%dataptr(mfi)

       do idim = 1, amrex_spacedim
          tbx = bx
          call tbx%nodalize(idim)
          call flux(idim)%resize(tbx,ncomp)
          call tbx%grow(substep)
       end do

       pfx => flux(1)%dataptr()
       pfy => flux(2)%dataptr()
       pfz => flux(3)%dataptr()
       
       ! Increment solution at given mfi tilebox
       call increment_enthalpy(time, bx%lo, bx%hi, &
                               pin, lbound(pin),     ubound(pin),     &
                               pout,    lbound(pout),    ubound(pout),    &
                               ptempin, lbound(ptempin), ubound(ptempin), &
                               ptemp,   lbound(ptemp),   ubound(ptemp),   &
                               pfx, lbound(pfx), ubound(pfx), &
                               pfy, lbound(pfy), ubound(pfy), &
                               pfz, lbound(pfz), ubound(pfz), &
                               amrex_geom(lev), dt)

       if (do_reflux) then

          do idim = 1, amrex_spacedim

             pf => fluxes(idim)%dataptr(mfi)
             pfab => flux(idim)%dataptr()
             tbx = mfi%nodaltilebox(idim)
             pf       (tbx%lo(1):tbx%hi(1), tbx%lo(2):tbx%hi(2), tbx%lo(3):tbx%hi(3), :)  = & 
                  pfab(tbx%lo(1):tbx%hi(1), tbx%lo(2):tbx%hi(2), tbx%lo(3):tbx%hi(3), :) 
              
          end do
          
       end if
       
       ! Find melt interface y position 
       if (lev.eq.amrex_max_level) then
          call get_melt_pos(bx%lo, bx%hi,                      	&
                          ptemp, lbound(ptemp), ubound(ptemp), 	&
        		  amrex_geom(lev))
        end if  

    end do

    call amrex_mfiter_destroy(mfi)
    do idim = 1, amrex_spacedim
       call amrex_fab_destroy(flux(idim))
    end do
    
!$omp end parallel

    if (do_reflux) then

       ! The fluxes have already been scaled by dt and area in the increment_enthalpy subroutine
       if (lev > 0) then
          call flux_reg(lev)%fineadd(fluxes, 1.0_amrex_real)
       end if

       if (lev < amrex_get_finest_level()) then
          call flux_reg(lev+1)%crseinit(fluxes, -1.0_amrex_real)
       end if

       do idim = 1, amrex_spacedim
          call amrex_multifab_destroy(fluxes(idim))
       end do
       
    end if
       
    call amrex_multifab_destroy(phiborder)
    
 end subroutine advance


end module simulation_module
