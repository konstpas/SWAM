module simulation_module

  ! -----------------------------------------------------------------
  ! This module is used to run the simulation, i.e. to coordinate
  ! all the necessary calculations involved in the heat transfer,
  ! phase change and melt motion
  ! -----------------------------------------------------------------
  
  use amrex_amr_module

  implicit none
  
  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: run_simulation

contains

  ! -----------------------------------------------------------------
  ! Subroutine used to run the simulation
  ! -----------------------------------------------------------------
  subroutine run_simulation()

    use amr_data_module, only : t_new, &
                                stepno, &
                                dt
    use read_input_module, only : time_step_max, &
                                  time_max, &
                                  plot_int, &
                                  heat_solver, &
                                  do_reflux, &
                                  heat_solve
    use plotfile_module, only: writeplotfile

    ! Local variables
    integer :: last_plot_file_step
    integer :: step
    integer :: lev
    integer :: substep
    real(amrex_real) :: time

    ! Initialize time
    time = t_new(0)

    ! Initialize counter for output
    last_plot_file_step = 0;

    ! Automatically disable reflux if the heat equation is not solved explicitly
    if (heat_solver.ne."explicit" .or. .not.heat_solve) then
       do_reflux = .false.
    end if
    
    ! Output initial configuration
    if (plot_int .gt. 0) call writeplotfile

    ! Loop over time-steps
    do step = stepno(0), time_step_max - 1

       ! Print timestep information on screen (start)
       if (amrex_parallel_ioprocessor()) then
          print *, ""
          print *, "STEP", step+1, "starts ..."
       end if

       ! Compute time step size
       call compute_dt
		 
       ! Advance all levels of one time step (with or withouth subcycling)
       if (heat_solve) then
          if (heat_solver.eq."explicit") then
             lev = 0
             substep = 1
             call advance_one_timestep_subcycling(lev, time, substep)
          else if (heat_solver.eq."implicit") then
             call advance_one_timestep(time)
          else
             STOP "Unknown heat solver prescribed in input"
          end if
       else
          call advance_one_timestep(time)
       end if
       
       ! Update time on all levels
       time = time + dt(0)
       do lev = 0, amrex_max_level
          t_new(lev) = time
       end do

       ! Print timestep information on screen (end)
       if (amrex_parallel_ioprocessor()) then
          print *, "STEP", step+1, "end. TIME =", time, "DT =", dt(0)
       end if

       ! Print output to file
       if (plot_int .gt. 0 .and. mod(step+1,plot_int) .eq. 0) then
          last_plot_file_step = step+1
          call writeplotfile
       end if

       ! Stopping criteria
       if (time .ge. time_max) exit

    end do
    
  end subroutine run_simulation

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the timestep for each level
  ! -----------------------------------------------------------------
  subroutine compute_dt()
    
    use amr_data_module, only : dt, &
                                nsubsteps
    use read_input_module, only : dt_change_max, &
                                  time_dt, &
                                  heat_solver

    ! Local variables
    integer :: lev
    integer :: nlevs
    integer ::n_factor
    real(amrex_real) :: dt_0
    real(amrex_real), allocatable :: dt_tmp(:)

    nlevs = amrex_get_numlevels()
    
    allocate(dt_tmp(0:nlevs-1))

    ! Estimate timestep in order to avoid instabilities
    do lev = 0, nlevs-1
       call check_timestep_stability(lev, dt_tmp(lev))
    end do
    call amrex_parallel_reduce_min(dt_tmp, nlevs)

    ! Ensure that the timestep does not change too much
    ! from one timestep to the next
    dt_0 = dt_tmp(0)
    n_factor = 1
    do lev = 0, nlevs-1
       ! Limit the change of the timestep between two successive timesteps
       dt_tmp(lev) = min(dt_tmp(lev), dt_change_max*dt(lev))
       n_factor = n_factor * nsubsteps(lev)
       dt_0 = min(dt_0, n_factor*dt_tmp(lev))
       ! Cap the timestep to the value prescribed in input
       dt_0 = min(dt_0, time_dt)
    end do

    ! Compute timestep for all levels
    dt(0) = dt_0
    do lev = 1, nlevs-1
       if (heat_solver.eq."explicit") then
          dt(lev) = dt(lev-1) / nsubsteps(lev)
       else
          dt(lev) = dt(lev-1)
       end if
    end do
    
  end subroutine compute_dt


  ! -----------------------------------------------------------------
  ! Subroutine used to estimate the timestep based on the stability
  ! criteria for the solver of the heat equation
  ! -----------------------------------------------------------------
  subroutine check_timestep_stability(lev, dt)

    use amr_data_module, only : max_melt_vel
    use read_input_module, only : cfl, &
                                  heat_solver, &
                                  heat_solve, &
                                  time_dt
    use material_properties_module, only : max_diffus 

    ! Input and output variables 
    integer, intent(in) :: lev
    real(amrex_real), intent(out) :: dt

    ! Local variables 
    real(amrex_real) :: dxsqr
    real(amrex_real) :: dx
    real(amrex_real) :: dy

    ! Grid resolutions
    dx = amrex_geom(lev)%dx(1)
    dy = amrex_geom(lev)%dx(2)
    
    ! Von Neumann stability condition (enforced)
    if (heat_solver.eq."explicit" .and. heat_solve) then
       dxsqr= (1.0/dx**2 + 1.0/dy**2) 
       dt = 0.5/(dxsqr*max_diffus)
       dt = dt * cfl
    else
       dt = time_dt
    end if

    ! CFL stability condition (only checked)
    if (max_melt_vel*dt.ge.dx) then
       print *, "CFL stability condition is not satisfied"
    end if
    
  end subroutine check_timestep_stability

  
  ! -----------------------------------------------------------------
  ! Subroutine used to advance the simulation of one timestep on
  ! one level with the fully explicit method with subcycling.
  ! Note that the subroutine is recursive, i.e. it calls itself.
  ! -----------------------------------------------------------------
  recursive subroutine advance_one_timestep_subcycling(lev, time, substep)

    use read_input_module, only : do_reflux
    use amr_data_module, only : t_old, &
                                t_new, &
                                phi_new, &
                                flux_reg, &
                                stepno, &
                                nsubsteps, &
                                dt  
    use regrid_module, only : averagedownto

    ! Input and output variables
    integer, intent(in) :: lev
    integer, intent(in) :: substep
    real(amrex_real), intent(in) :: time

    ! Local variables
    integer :: fine_substep    
    
    ! Regridding
    call check_regridding(lev, time)    
    
    ! Advance solution
    stepno(lev) = stepno(lev)+1
    t_old(lev) = time
    t_new(lev) = time + dt(lev)
    call advance_one_level_subcycling(lev, time, dt(lev), substep)

    ! Propagate solution and synchronize levels
    if (lev .lt. amrex_get_finest_level()) then
       
       do fine_substep = 1, nsubsteps(lev+1)
          call advance_one_timestep_subcycling(lev+1, time+(fine_substep-1)*dt(lev+1), fine_substep) 
       end do

       ! Update coarse solution at coarse-fine interface via dPhidt = -div(+F) where F is stored flux (flux_reg)
       if (do_reflux) then
          call flux_reg(lev+1)%reflux(phi_new(lev), 1.0_amrex_real) 
       end if
       
       ! Define the solution at the coarser level to be the average of the solution at the finer level
       call averagedownto(lev)
        
    end if
  
  end subroutine advance_one_timestep_subcycling

  ! -----------------------------------------------------------------
  ! Subroutine used to check if it is time to regrid
  ! -----------------------------------------------------------------
  subroutine check_regridding(lev, time)

    use read_input_module, only : regrid_int
    use amr_data_module, only : last_regrid_step, stepno, dt

    ! Input and output variables
    integer, intent(in) :: lev
    real(amrex_real), intent(in) :: time
    
    ! Local variables
    integer :: ilev
    integer :: finest_level
    integer :: old_finest_level
    
    if (regrid_int .gt. 0) then
       
       if (lev .lt. amrex_max_level .and. stepno(lev) .gt. last_regrid_step(lev)) then

          if (mod(stepno(lev), regrid_int) .eq. 0) then

             old_finest_level = amrex_get_finest_level()
             call amrex_regrid(lev, time)
             finest_level = amrex_get_finest_level()

             do ilev = lev, finest_level
                last_regrid_step(ilev) = stepno(ilev)
             end do

             do ilev = old_finest_level+1, finest_level
                dt(ilev) = dt(ilev-1) / amrex_ref_ratio(ilev-1)
             end do
             
          end if
          
       end if
       
    end if
    
  end subroutine check_regridding

  
  ! -----------------------------------------------------------------
  ! Subroutine used to advance the shallow water solver and the
  ! heat equation solver of one time step at a given level. Only
  ! used if the fully explicit solver with subcycling is used.
  ! -----------------------------------------------------------------
  subroutine advance_one_level_subcycling(lev, time, dt, substep)

    use read_input_module, only : sw_solve, &
                                  heat_solve
    use heat_transfer_module, only : advance_heat_solver_explicit_level
    use shallow_water_module, only : advance_SW

    ! Input and output variables 
    integer, intent(in) :: lev
    integer, intent(in) :: substep
    real(amrex_real), intent(in) :: time
    real(amrex_real), intent(in) :: dt

    ! Advance shallow-water equations (only at max level)
    if (sw_solve .and. lev.eq.amrex_max_level) then
       call advance_SW(time)    
    end if 

    ! Advance heat equation
    if(heat_solve) then
       call advance_heat_solver_explicit_level(lev, time, dt, substep)
    end if

    
  end subroutine advance_one_level_subcycling

  ! -----------------------------------------------------------------
  ! Subroutine used to advance the simulation of one timestep with
  ! the implicit method for the heat equation.
  ! -----------------------------------------------------------------
  subroutine advance_one_timestep(time)

    use amr_data_module, only : t_old, t_new, stepno, dt
    use regrid_module, only : averagedown

    ! Input and output variables
    real(amrex_real), intent(in) :: time
    
    ! Local variables
    integer :: ilev
     
    ! Regridding
    do ilev = 0,amrex_max_level
       call check_regridding(ilev,time)
    end do
    
    ! Advance solution
    do ilev = 0, amrex_max_level
       t_old(ilev) = time
       t_new(ilev) = time + dt(0)
       stepno(ilev) = stepno(ilev) + 1
    end do

    call advance_all_levels(time,dt(0)) 

    call averagedown
    
  end subroutine advance_one_timestep
  
  ! -----------------------------------------------------------------
  ! Subroutine used to advance the shallow water solver and the
  ! heat equation solver of one time step for all levels. Only
  ! used if the implicit solver for the heat equation is employed.
  ! -----------------------------------------------------------------
  subroutine advance_all_levels(time, dt)

    use read_input_module, only : sw_solve, &
                                  heat_solve
    use material_properties_module, only : get_temp
    use shallow_water_module, only : advance_SW
    use heat_transfer_module, only : advance_heat_solver_implicit
    
    ! Input and output variables
    real(amrex_real), intent(in) :: dt
    real(amrex_real), intent(in) :: time

    ! Advance SW equations (only at max level)
    if (sw_solve) then
       call advance_SW(time)    
    end if

    ! Advance heat equation
    if(heat_solve) then
       call advance_heat_solver_implicit(time, dt)
    end if

  end subroutine advance_all_levels

  
end module simulation_module
