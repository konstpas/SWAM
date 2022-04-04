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

    use amr_data_module, only : t_new, stepno, dt
    use read_input_module, only : max_step, stop_time, plot_int, heat_solver, do_reflux, solve_heat
    use plotfile_module, only: writeplotfile, write2dplotfile
    use energy_module, only: sum_enthalpy

    ! Local variables
    integer :: last_plot_file_step
    integer :: step
    integer :: lev
    integer :: substep
    real(amrex_real) :: cur_time
    real(amrex_real) :: total_enthalpy
    
    ! Initialize time
    cur_time = t_new(0)

    ! Initialize counter for output
    last_plot_file_step = 0;

    ! Automatically disable reflux if the heat equation is not solved explicitly
    if (heat_solver.ne."explicit" .or. .not.solve_heat) then
      do_reflux = .false.
    end if

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
       call compute_dt
		 
       ! Advance all levels of one time step
       if (heat_solver.eq."explicit" .and. solve_heat) then
          lev = 0
          substep = 1
          call advance_one_timestep_subcycling(lev, cur_time, substep)
       else
          call advance_one_timestep(cur_time)
       end if
      
       ! Get total enthalpy in the domain
       if(solve_heat) then
         call sum_enthalpy(total_enthalpy)
       else
         total_enthalpy = 0.0
       end if
 	
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

  end subroutine run_simulation

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the timestep for each level
  ! -----------------------------------------------------------------
  subroutine compute_dt()
    
    use amr_data_module, only : dt, nsubsteps
    use read_input_module, only : dt_change_max, in_dt, heat_solver

    ! Local variables
    integer :: lev
    integer :: nlevs
    integer :: n_factor
    real(amrex_real) :: dt_0
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
       ! Limit the change of the timestep between two successive timesteps
       dt_tmp(lev) = min(dt_tmp(lev), dt_change_max*dt(lev))
       n_factor = n_factor * nsubsteps(lev)
       dt_0 = min(dt_0, n_factor*dt_tmp(lev))
       ! Cap the timestep to the value prescribed in input
       dt_0 = min(dt_0, in_dt)
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
  ! criteria for the solver of the energy equation
  ! -----------------------------------------------------------------
  subroutine est_timestep(lev, dt)

    use read_input_module, only : cfl, heat_solver, solve_heat
    use material_properties_module, only : max_diffus 

    ! Input and output variables 
    integer, intent(in) :: lev
    real(amrex_real), intent(out) :: dt

    ! Local variables 
    real(amrex_real) :: dxsqr

    if (heat_solver.eq."explicit" .and. solve_heat) then

      ! NOTE: There are two stability criteria to take into
      ! account, the von Neumann stability and the CFL
      ! condition. The CFL is not yet implemented
      
      ! Von Neumann stability criterion 
      dxsqr= (1/amrex_geom(lev)%dx(1)**2 + &
               1/amrex_geom(lev)%dx(2)**2 + &
               1/amrex_geom(lev)%dx(3)**2) 
      dt = 0.5/(dxsqr*max_diffus)
      dt = dt * cfl

    else

       ! Implicit diffusion, explicit advection

       ! NOTE: The CFL condition is not yet implemented
       dt = 1.0e100_amrex_real
      
    end if
  end subroutine est_timestep

  
  ! -----------------------------------------------------------------
  ! Subroutine used to advance the simulation of one timestep. Note
  ! that the subroutine is recursive, i.e. it calls itself
  ! -----------------------------------------------------------------
  recursive subroutine advance_one_timestep_subcycling(lev, time, substep)

    use read_input_module, only : do_reflux
    use amr_data_module, only : t_old, t_new, phi_new, flux_reg, stepno, nsubsteps, dt  
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

   use read_input_module, only : solve_sw
   use heat_transfer_module, only : advance_heat_solver_explicit_level
   use domain_module, only : reset_melt_pos
   use shallow_water_module, only : advance_SW

   ! Input and output variables 
   integer, intent(in) :: lev
   integer, intent(in) :: substep
   real(amrex_real), intent(in) :: time
   real(amrex_real), intent(in) :: dt

   ! Advance shallow-water equations (only at max level)
   if (solve_sw .and. lev.eq.amrex_max_level) then
      call advance_SW(time)    
   end if 
   
   ! Set melt interface position array equal to free interface position array 
   ! Since melt layer may span several tile boxes in y-direction (in mfiterator below), we cannot reset within each loop 
   ! Therefore we reset melt position after solving SW, and before propagating temperature 
   ! Melt position is then found after heat has been propagated 
   call reset_melt_pos 

   ! Advance heat equation
   call advance_heat_solver_explicit_level(lev, time, dt, substep)

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

    use read_input_module, only : solve_sw, solve_heat
    use domain_module, only : reset_melt_pos
    use material_properties_module, only : get_temp
    use shallow_water_module, only : advance_SW
    use heat_transfer_module, only : advance_heat_solver_implicit
   
    ! Input and output variables
    real(amrex_real), intent(in) :: dt
    real(amrex_real), intent(in) :: time

    ! Propagate SW equations (only at max level)
    if (solve_sw) then
       call advance_SW(time)    
    end if

    if(solve_heat) then
      ! Set melt interface position array equal to free interface position array 
      ! Since melt layer may span several tile boxes in y-direction (in mfiterator below), we cannot reset within each loop 
      ! Therefore we reset melt position after solving SW, and before propagating temperature 
      ! Melt position is then found after heat has been propagated 
      call reset_melt_pos 
      
      ! Advance heat equation
      call advance_heat_solver_implicit(time, dt)
   end if

 end subroutine advance_all_levels

   
end module simulation_module
