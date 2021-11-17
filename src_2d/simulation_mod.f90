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
    use read_input_module, only : max_step, stop_time, plot_int, heat_solver
    use plotfile_module, only: writeplotfile, write1dplotfile
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
       if (heat_solver.eq."explicit") then
          lev = 0
          substep = 1
          call advance_one_timestep_explicit(lev, cur_time, substep)
       else
          call advance_one_timestep_implicit(cur_time)
       end if
          
       ! Get total enthalpy in the domain
       call sum_enthalpy(total_enthalpy)
 	
       ! Update time on all levels
       cur_time = cur_time + dt(0)
       do lev = 0, amrex_max_level
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
          call write1dplotfile
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
    integer ::n_factor
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
  ! criteria for the solver of the heat equation
  ! -----------------------------------------------------------------
  subroutine est_timestep(lev, dt)

    use read_input_module, only : cfl, heat_solver
    use material_properties_module, only : max_diffus 

    ! Input and output variables 
    integer, intent(in) :: lev
    real(amrex_real), intent(out) :: dt

    ! Local variables 
    real(amrex_real) :: dxsqr

    if (heat_solver.eq."explicit") then

       ! Fully explicit solver
       
       ! NOTE: There are two stability criteria to take into
       ! account, the von Neumann stability and the CFL
       ! condition. Here we only consider the von Neumann
       ! stability condition which is usually more
       ! stringent
       
       ! Von Neumann stability criterion 
       dxsqr= (1/amrex_geom(lev)%dx(1)**2 + &
            1/amrex_geom(lev)%dx(2)**2) 
       dt = 0.5/(dxsqr*max_diffus)
       dt = dt * cfl

    else

       ! Implicit diffusion, explicit advection

       ! NOTE: The CFL condition is not yet implemented
       dt = 1.0e100_amrex_real
       
    end if
       
  end subroutine est_timestep

  
  ! -----------------------------------------------------------------
  ! Subroutine used to advance the simulation of one timestep on
  ! one level with the fully explicit method with subcycling.
  ! Note that the subroutine is recursive, i.e. it calls itself.
  ! -----------------------------------------------------------------
  recursive subroutine advance_one_timestep_explicit(lev, time, substep)

    use read_input_module, only : regrid_int, do_reflux
    use amr_data_module, only : t_old, t_new, phi_old, phi_new, flux_reg, stepno, nsubsteps, dt  
    use regrid_module, only : averagedownto

    ! Input and output variables
    integer, intent(in) :: lev
    integer, intent(in) :: substep
    real(amrex_real), intent(in) :: time

    ! Local variables
    integer, allocatable, save :: last_regrid_step(:)
    integer :: finest_level
    integer :: fine_substep    
    integer :: k
    integer :: old_finest_level
    
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
    t_old(lev) = time
    t_new(lev) = time + dt(lev)
    call amrex_multifab_swap(phi_old(lev), phi_new(lev))
    call advance_one_level_explicit(lev, time, dt(lev), substep)

    ! Propagate solution and synchronize levels
    if (lev .lt. amrex_get_finest_level()) then
       
       do fine_substep = 1, nsubsteps(lev+1)
          call advance_one_timestep_explicit(lev+1, time+(fine_substep-1)*dt(lev+1), fine_substep) 
       end do

       ! Update coarse solution at coarse-fine interface via dPhidt = -div(+F) where F is stored flux (flux_reg)
       if (do_reflux) then
          call flux_reg(lev+1)%reflux(phi_new(lev), 1.0_amrex_real) 
       end if
       
       ! Define the solution at the coarser level to be the average of the solution at the finer level
       call averagedownto(lev)
        
    end if
  
  end subroutine advance_one_timestep_explicit


  ! -----------------------------------------------------------------
  ! Subroutine used to advance the shallow water solver and the
  ! heat equation solver of one time step at a given level. Only
  ! used if the fully explicit solver with subcycling is used.
  ! -----------------------------------------------------------------
  subroutine advance_one_level_explicit(lev, time, dt, substep)

    use read_input_module, only : do_reflux, solve_sw, temp_fs
    use amr_data_module, only : phi_new, temp, idomain, flux_reg  
    use regrid_module, only : fillpatch
    use heat_transfer_module, only : get_idomain, get_melt_pos, reset_melt_pos, increment_enthalpy
    use heat_transfer_fixT_module, only : increment_enthalpy_fixT
    use material_properties_module, only : get_temp
    use shallow_water_module, only : increment_SW

    ! Input and output variables 
    integer, intent(in) :: lev
    integer, intent(in) :: substep
    real(amrex_real), intent(in) :: time
    real(amrex_real), intent(in) :: dt

    ! Local variables
    integer, parameter :: nghost = 1 ! number of ghost points in each spatial direction 
    integer :: ncomp
    integer :: idim
    logical :: nodal(2) ! logical for flux multifabs 
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pin
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pout
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptempin
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptemp
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pfx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pfy
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pf
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pfab
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidin
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidout
    type(amrex_multifab) :: phi_tmp ! Enthalpy multifab with ghost points
    type(amrex_multifab) :: temp_tmp ! Temperature multifab with ghost points
    type(amrex_multifab) :: idomain_tmp ! Idomain multifab to distinguish material and background
    type(amrex_mfiter) :: mfi ! Multifab iterator
    type(amrex_box) :: bx
    type(amrex_box) :: tbx
    type(amrex_fab) :: flux(amrex_spacedim)
    type(amrex_multifab) :: fluxes(amrex_spacedim)
    type(amrex_geometry) :: geom

    ! Get geometry
    geom = amrex_geom(lev)
    
    ! Get number of components
    ncomp = phi_new(lev)%ncomp()

    ! Initialize fluxes 
    if (do_reflux) then
       do idim = 1, amrex_spacedim
          nodal = .false.
          nodal(idim) = .true.
          call amrex_multifab_build(fluxes(idim), phi_new(lev)%ba, phi_new(lev)%dm, ncomp, 0, nodal)
       end do
    end if

    ! Build enthalpy and temperature multifabs with ghost points
    call amrex_multifab_build(phi_tmp, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, nghost) 
    call amrex_multifab_build(temp_tmp, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, nghost)
    call amrex_multifab_build(idomain_tmp, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, nghost)
    
    ! Fill enthalpy multifab
    call fillpatch(lev, time, phi_tmp)

    ! Swap idomain solution before computing the surface deformation
    call amrex_multifab_swap(idomain_tmp, idomain(lev))
    
    ! Propagate SW equations (only at max level)
    if (solve_sw .and. lev.eq.amrex_max_level) then
       
       call amrex_mfiter_build(mfi, idomain_tmp, tiling=.false.)  
    						 
       do while(mfi%next())

          ! Box
          bx = mfi%validbox()   
          
          ! Pointers
          ptemp   => temp(lev)%dataptr(mfi)
          pidin => idomain_tmp%dataptr(mfi)
          
          ! Increment shallow water solution
          call increment_SW(bx%lo, bx%hi, &
                            ptemp, lbound(ptemp), ubound(ptemp), &
                            pidin, lbound(pidin), ubound(pidin), dt)
       
       end do

       ! Clean memory
       call amrex_mfiter_destroy(mfi)    
    
    end if 
    
    ! Set melt interface position array equal to free interface position array 
    ! Since melt layer may span several tile boxes in y-direction (in mfiterator below), we cannot reset within each loop 
    ! Therefore we reset melt position after solving SW, and before propagating temperature 
    ! Melt position is then found after heat has been propagated 
    call reset_melt_pos 
       
    ! Increment heat solver on all levels
    do idim = 1, amrex_spacedim
       call flux(idim)%reset_omp_private()
    end do
    
    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.false.)  
    						 
    do while(mfi%next())

       ! Box
       bx = mfi%validbox()   

       ! Pointers
       pin     => phi_tmp%dataptr(mfi)
       pout    => phi_new(lev)%dataptr(mfi)
       ptempin => temp_tmp%dataptr(mfi)
       ptemp   => temp(lev)%dataptr(mfi)
       do idim = 1, amrex_spacedim
          tbx = bx
          call tbx%nodalize(idim)
          call flux(idim)%resize(tbx,ncomp)
          call tbx%grow(substep)
       end do
       pfx => flux(1)%dataptr()
       pfy => flux(2)%dataptr()
       pidin => idomain_tmp%dataptr(mfi)
       pidout => idomain(lev)%dataptr(mfi)

       ! Get temperature corresponding to the enthalpy
       call get_temp(lbound(ptempin), ubound(ptempin), &
                     pin, lbound(pin), ubound(pin), &
                     ptempin, lbound(ptempin), ubound(ptempin),.true.)
       
       ! Get configuration of the system after the deformation
       call get_idomain(geom%get_physical_location(bx%lo), geom%dx, &
                        bx%lo, bx%hi, &
                        pidout, lbound(pidout), ubound(pidout), &
                        ptempin, lbound(ptempin), ubound(ptempin))
          
       ! Increment enthalpy at given box depending on the condition of the free surface
       if (temp_fs.gt.0) then
          call increment_enthalpy_fixT(bx%lo, bx%hi, &
                                       pin, lbound(pin),     ubound(pin),     &
                                       pout,    lbound(pout),    ubound(pout),    &
                                       ptempin, lbound(ptempin), ubound(ptempin), &
                                       ptemp,   lbound(ptemp),   ubound(ptemp),   &
                                       pfx, lbound(pfx), ubound(pfx), &
                                       pfy, lbound(pfy), ubound(pfy), &
                                       pidin, lbound(pidin), ubound(pidin), &
                                       pidout, lbound(pidout), ubound(pidout), &
                                       geom, dt)
       else
          call increment_enthalpy(time, bx%lo, bx%hi, &
                                  pin, lbound(pin),     ubound(pin),     &
                                  pout,    lbound(pout),    ubound(pout),    &
                                  ptempin, lbound(ptempin), ubound(ptempin), &
                                  ptemp,   lbound(ptemp),   ubound(ptemp),   &
                                  pfx, lbound(pfx), ubound(pfx), &
                                  pfy, lbound(pfy), ubound(pfy), &
                                  pidin, lbound(pidin), ubound(pidin), &
                                  pidout, lbound(pidout), ubound(pidout), &
                                  geom, dt)
       end if
       
       ! Update pointers for flux registers
       if (do_reflux) then

          do idim = 1, amrex_spacedim

             pf => fluxes(idim)%dataptr(mfi)
             pfab => flux(idim)%dataptr()
             tbx = mfi%nodaltilebox(idim)
             pf(tbx%lo(1):tbx%hi(1), tbx%lo(2):tbx%hi(2), tbx%lo(3):tbx%hi(3), :)  = & 
                  pfab(tbx%lo(1):tbx%hi(1), tbx%lo(2):tbx%hi(2), tbx%lo(3):tbx%hi(3), :) 
              
          end do
          
       end if
       
       ! Find melt interface y position 
       if (lev.eq.amrex_max_level) then
          call get_melt_pos(bx%lo, bx%hi, &
                            pidout, lbound(pidout), ubound(pidout), &
                            geom)
       end if

    end do

    ! Clean memory
    call amrex_mfiter_destroy(mfi)
    do idim = 1, amrex_spacedim
       call amrex_fab_destroy(flux(idim))
    end do
    
    ! Update flux registers (fluxes have already been scaled by dt and area in the increment_enthalpy subroutine)
    if (do_reflux) then

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

    ! Clean memory
    call amrex_multifab_destroy(phi_tmp)
    call amrex_multifab_destroy(temp_tmp)
    call amrex_multifab_destroy(idomain_tmp)

  end subroutine advance_one_level_explicit


  ! -----------------------------------------------------------------
  ! Subroutine used to advance the simulation of one timestep with
  ! the implicit method for the heat equation.
  ! -----------------------------------------------------------------
  subroutine advance_one_timestep_implicit(time)

    use read_input_module, only : regrid_int, in_dt, heat_solver
    use amr_data_module, only : phi_old, phi_new, t_old, t_new, stepno, dt
    use regrid_module, only : averagedown

    ! Input and output variables
    real(amrex_real), intent(in) :: time
    
    ! Local variables
    integer :: ilev
    integer :: jlev
    integer, allocatable, save :: last_regrid_step(:)
    integer :: finest_level
    integer :: fine_substep    
    integer :: old_finest_level
    
    ! Regridding 
    if (regrid_int .gt. 0) then
       
       if (.not.allocated(last_regrid_step)) then

          allocate(last_regrid_step(0:amrex_max_level))
          last_regrid_step = 0
          
       end if

       do ilev = 0,amrex_max_level
          
          if (ilev .lt. amrex_max_level .and. stepno(ilev) .gt. last_regrid_step(ilev)) then

             if (mod(stepno(ilev), regrid_int) .eq. 0) then

                old_finest_level = amrex_get_finest_level()
                call amrex_regrid(ilev, time)
                finest_level = amrex_get_finest_level()
                
                do jlev = ilev, finest_level
                   last_regrid_step(jlev) = stepno(jlev)
                end do
                
                do jlev = old_finest_level+1, finest_level
                   dt(jlev) = dt(jlev-1)
                end do
                
             end if
             
          end if
       
       end do
    end if
    
    ! Advance solution
    do ilev = 0, amrex_max_level
       t_old(ilev) = time
       t_new(ilev) = time + dt(0)
       stepno(ilev) = stepno(ilev) + 1
       call amrex_multifab_swap(phi_old(ilev), phi_new(ilev))       
    end do

    call advance_all_levels_implicit(time,dt(0)) 

    call averagedown
    
  end subroutine advance_one_timestep_implicit

  ! -----------------------------------------------------------------
  ! Subroutine used to advance the shallow water solver and the
  ! heat equation solver of one time step for all levels. Only
  ! used if the implicit solver for the heat equation is employed.
  ! -----------------------------------------------------------------
  subroutine advance_all_levels_implicit(time, dt)

    use read_input_module, only : solve_sw, ls_agglomeration, &
                                  ls_consolidation, ls_max_coarsening_level, &
                                  ls_linop_maxorder, ls_bottom_solver, &
                                  ls_bottom_verbose, ls_max_fmg_iter, &
                                  ls_max_iter, ls_verbose

    use amr_data_module, only : phi_new, temp, idomain
    use regrid_module, only : fillpatch
    use heat_transfer_module, only : get_idomain, get_melt_pos, reset_melt_pos, increment_enthalpy
    use material_properties_module, only : get_temp
    use shallow_water_module, only : increment_SW
    use amrex_linear_solver_module
    
    ! Input and output variables
    real(amrex_real), intent(in) :: dt
    real(amrex_real), intent(in) :: time
  
    ! Local variables
    logical :: nodal(2)
    integer, parameter :: nghost = 1 ! number of ghost points in each spatial direction 
    integer :: ncomp
    integer :: idim
    integer :: ilev
    real(amrex_real) :: err ! Error from the solution of the linear system of equations
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pin
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pout
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptempin
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptemp
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidin
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidout
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pac
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pbc
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: prhs
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pu
    type(amrex_multifab) :: phi_tmp(0:amrex_max_level) ! Enthalpy multifab with ghost points
    type(amrex_multifab) :: temp_tmp(0:amrex_max_level) ! Temperature multifab with ghost points
    type(amrex_multifab) :: idomain_tmp(0:amrex_max_level) ! Idomain multifab to distinguish material and background
    type(amrex_multifab) :: ls_rhs(0:amrex_max_level) ! multifab for the right hand side of the linear system of equations
    type(amrex_multifab) :: ls_alpha(0:amrex_max_level) ! multifab for the first term on the left hand side of the linear system of equations
    type(amrex_multifab) :: ls_beta(amrex_spacedim,0:amrex_max_level) ! multifab for the second term on the left hand side of the linear system of equations
    type(amrex_multifab) :: nullmf ! Null multifab used to define the boundary condition at the coarse-fine interface
    type(amrex_mfiter) :: mfi ! Multifab iterator
    type(amrex_box) :: bx
    type(amrex_boxarray) :: ba(0:amrex_max_level) ! Box array
    type(amrex_distromap):: dm(0:amrex_max_level) ! Distribution mapping
    type(amrex_abeclaplacian) :: abeclap ! Object for the solution of the linear systems of equations
    type(amrex_multigrid) :: multigrid ! 
    
    ! Build necessary multifabs on all levels
    do ilev = 0, amrex_max_level

       ! Enthalpy and temperature multifabs with ghost points
       ncomp = phi_new(ilev)%ncomp()
       call amrex_multifab_build(phi_tmp(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, ncomp, nghost) 
       call amrex_multifab_build(temp_tmp(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, ncomp, nghost)

       ! Temporary idomain to store the deformation computed via shallow water
       call amrex_multifab_build(idomain_tmp(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, ncomp, nghost)

       ! Multifabs for the linear solver
       call amrex_multifab_build(ls_rhs(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, ncomp, 0)
       call amrex_multifab_build(ls_alpha(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, ncomp, 0)
       do idim = 1, amrex_spacedim
          nodal = .false.
          nodal(idim) = .true.
          call amrex_multifab_build(ls_beta(idim,ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, ncomp, 0, nodal)
       end do
 
       ! Get boxarray and distribution mapping on level
       ba(ilev) = phi_new(ilev)%ba
       dm(ilev) = phi_new(ilev)%dm
       
       ! Fill enthalpy multifab
       call fillpatch(ilev, time, phi_tmp(ilev))
       
       ! Swap idomain solution before computing the surface deformation
       call amrex_multifab_swap(idomain_tmp(ilev), idomain(ilev))
       
    end do
    

    ! Propagate SW equations (only at max level)
    if (solve_sw) then
       
       call amrex_mfiter_build(mfi, idomain_tmp(amrex_max_level), tiling=.false.)  
    						 
       do while(mfi%next())

          ! Box
          bx = mfi%validbox()   
          
          ! Pointers
          ptemp   => temp(amrex_max_level)%dataptr(mfi)
          pidin => idomain_tmp(amrex_max_level)%dataptr(mfi)
          
          ! Increment shallow water solution
          call increment_SW(bx%lo, bx%hi, &
                            ptemp, lbound(ptemp), ubound(ptemp), &
                            pidin, lbound(pidin), ubound(pidin), dt)
       
       end do

       ! Clean memory
       call amrex_mfiter_destroy(mfi)    
    
    end if 
    
    ! Set melt interface position array equal to free interface position array 
    ! Since melt layer may span several tile boxes in y-direction (in mfiterator below), we cannot reset within each loop 
    ! Therefore we reset melt position after solving SW, and before propagating temperature 
    ! Melt position is then found after heat has been propagated 
    call reset_melt_pos 
       
    ! Prepare all levels for heat equation solution
    do ilev = 0, amrex_max_level 
    
       call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.false.)  
    						 
       do while(mfi%next())

          ! Box
          bx = mfi%validbox()   
          
          ! Pointers
          pin     => phi_tmp(ilev)%dataptr(mfi)
          pout    => phi_new(ilev)%dataptr(mfi)
          ptempin => temp_tmp(ilev)%dataptr(mfi)
          pidin   => idomain_tmp(ilev)%dataptr(mfi)
          pidout  => idomain(ilev)%dataptr(mfi)
          pac     => ls_alpha(ilev)%dataptr(mfi)
          prhs    => ls_rhs(ilev)%dataptr(mfi)
          
          ! Get temperature corresponding to the enthalpy
          call get_temp(lbound(ptempin), ubound(ptempin), &
                        pin, lbound(pin), ubound(pin), &
                        ptempin, lbound(ptempin), ubound(ptempin), .true.)

          ! Get configuration of the system after the deformation
          call get_idomain(amrex_geom(ilev)%get_physical_location(bx%lo), amrex_geom(ilev)%dx, &
                           bx%lo, bx%hi, &
                           pidout, lbound(pidout), ubound(pidout), &
                           ptempin, lbound(ptempin), ubound(ptempin))

          ! Update enthalpy according to deformation (still missing)
          
          ! Get alpha matrix for the linear solver (first term on left hand side)
          call get_alpha(bx%lo, bx%hi, &
                         ptempin, lbound(ptempin), ubound(ptempin), &
                         pac, lbound(pac), ubound(pac))

          ! Get beta matrix for the linear solver (second term on left hand side)
          do idim = 1,amrex_spacedim
             pbc => ls_beta(idim,ilev)%dataptr(mfi)
             call get_beta(bx%lo, bx%hi, idim, &
                           pidout, lbound(pidout), ubound(pidout), &
                           ptempin, lbound(ptempin), ubound(ptempin), &
                           pbc, lbound(pbc), ubound(pbc))
          end do

          ! Get right hand for the linear solver
          call get_rhs(bx%lo, bx%hi, time, dt, &
                       amrex_geom(ilev)%get_physical_location(bx%lo), amrex_geom(ilev)%dx, &
                       pidout, lbound(pidout), ubound(pidout), &
                       ptempin, lbound(ptempin), ubound(ptempin), &
                       prhs, lbound(prhs), ubound(prhs))


       end do

       call amrex_mfiter_destroy(mfi)
       
    end do

    ! -----------------------------------------------------------------------------------
    ! This should go in a dedicated subroutine
    ! ----------------------------------------------------------------------------------- 
    ! Solve linear system of equations
    call amrex_abeclaplacian_build(abeclap, amrex_geom, ba, dm, &
                                   metric_term=.false., agglomeration=ls_agglomeration, consolidation=ls_consolidation, &
                                   max_coarsening_level=ls_max_coarsening_level)

    call abeclap%set_maxorder(ls_linop_maxorder)

    ! This is set up to have homogeneous Neumann BC
    call abeclap%set_domain_bc([amrex_lo_neumann, amrex_lo_neumann], &
                               [amrex_lo_neumann, amrex_lo_neumann])

    do ilev = 0, amrex_max_level
       ! for problem with pure homogeneous Neumann BC, we could pass an empty multifab
       call abeclap % set_level_bc(ilev, nullmf)
    end do

    call abeclap%set_scalars(1.0_amrex_real, dt) ! A and B scalars
    do ilev = 0, amrex_max_level
       call abeclap%set_acoeffs(ilev, ls_alpha(ilev))
       call abeclap%set_bcoeffs(ilev, ls_beta(:,ilev))
    end do

    call amrex_multigrid_build(multigrid, abeclap)
    call multigrid%set_verbose(ls_verbose)
    call multigrid%set_bottom_verbose(ls_bottom_verbose)
    call multigrid%set_max_iter(ls_max_iter)
    call multigrid%set_max_fmg_iter(ls_max_fmg_iter)
    call multigrid%set_bottom_solver(ls_bottom_solver)
    
    err = multigrid%solve(temp_tmp, ls_rhs, 1.e-10_amrex_real, 0.0_amrex_real)
    
    call amrex_abeclaplacian_destroy(abeclap)
    call amrex_multigrid_destroy(multigrid)
    ! -----------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------
    
    ! Find melt interface y position
    call amrex_mfiter_build(mfi, idomain(amrex_max_level), tiling=.false.)
    do while(mfi%next())

       ! Box
       bx = mfi%validbox()   
       
       ! Pointers
       pidin => idomain(amrex_max_level)%dataptr(mfi)
       
       call get_melt_pos(bx%lo, bx%hi, &
                         pidin, lbound(pidin), ubound(pidin), &
                         amrex_geom(amrex_max_level))
       
    end do
    call amrex_mfiter_destroy(mfi)    

    ! Map solution to temperature and temperature to enthalpy (not really efficient, for testing purposes only)
    do ilev = 0, amrex_max_level
       
       call amrex_mfiter_build(mfi, temp_tmp(ilev), tiling=.false.)
       do while(mfi%next())

          ! Box
          bx = mfi%validbox()   
          
          ! Pointers
          ptemp   => temp(ilev)%dataptr(mfi)
          pu      => phi_new(ilev)%dataptr(mfi)
          ptempin => temp_tmp(ilev)%dataptr(mfi)
          
          ! Map solution on temperature
          call get_temp_ls(bx%lo, bx%hi, &
                           ptemp, lbound(ptemp), ubound(ptemp), &
                           ptempin, lbound(ptempin), ubound(ptempin))

          ! Map temperature on enthalpy
          call get_temp(bx%lo, bx%hi, &
                        ptemp, lbound(ptemp), ubound(ptemp), &
                        pu, lbound(pu), ubound(pu), .false.)
                    
       end do
       call amrex_mfiter_destroy(mfi)
       
    end do

        
    ! Clean memory
    do ilev = 0, amrex_max_level
       call amrex_multifab_destroy(phi_tmp(ilev))
       call amrex_multifab_destroy(temp_tmp(ilev))
       call amrex_multifab_destroy(idomain_tmp(ilev))
       call amrex_multifab_destroy(ls_rhs(ilev))
       call amrex_multifab_destroy(ls_alpha(ilev))
       do idim = 1, amrex_spacedim
          call amrex_multifab_destroy(ls_beta(idim,ilev))
       end do
       call amrex_distromap_destroy(dm)
    end do
    
  end subroutine advance_all_levels_implicit
  
  ! -----------------------------------------------------------------
  ! Subroutine used to fill the alpha matrix for the linear solver
  ! -----------------------------------------------------------------  
  subroutine get_alpha(lo, hi, &
                       temp, t_lo, t_hi, &
                       alpha, a_lo, a_hi)
  				
    use material_properties_module, only: get_mass_density, get_heat_capacity

    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)  
    integer, intent(in) :: t_lo(2), t_hi(2)
    integer, intent(in) :: a_lo(2), a_hi(2)    
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2))
    real(amrex_real), intent(out) :: alpha(a_lo(1):a_hi(1),a_lo(2):a_hi(2))

    ! Local variables
    integer :: i,j
    real(amrex_real) :: cp
    real(amrex_real) :: rho
    
    ! Fill alpha matrix
    do i = lo(1), hi(1)
       do j = lo(2), hi(2)          
          call get_mass_density(temp(i,j), rho)
          call get_heat_capacity(temp(i,j), cp)
          alpha(i,j) = rho*cp
       end do
    end do

  end subroutine get_alpha

  ! -----------------------------------------------------------------
  ! Subroutine used to fill the beta matrix for the linear solver
  ! -----------------------------------------------------------------  
  subroutine get_beta(lo, hi, dim, &
                      idom, id_lo, id_hi, &
                      temp, t_lo, t_hi, &
                      beta, b_lo, b_hi)
  				
    use material_properties_module, only: get_conductivity

    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: dim
    integer, intent(in) :: t_lo(2), t_hi(2)
    integer, intent(in) :: id_lo(2), id_hi(2)
    integer, intent(in) :: b_lo(2), b_hi(2)    
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2))
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2))
    real(amrex_real), intent(out) :: beta(b_lo(1):b_hi(1),b_lo(2):b_hi(2))

    ! Local variables
    integer :: i,j
    real(amrex_real) :: temp_face
    
    if (dim == 1) then ! x-direction
       do i = lo(1), hi(1)+1
          do j = lo(2), hi(2)
             
             if (nint(idom(i-1,j)).eq.0 .or. nint(idom(i,j)).eq.0) then
                 beta(i,j) = 0_amrex_real ! Suppress flux at the free surface
              else
                temp_face = (temp(i,j) + temp(i-1,j))/2_amrex_real
                call get_conductivity(temp_face, beta(i,j))
             end if

          end do
       end do

    else if (dim == 2) then ! y-direction
       do i = lo(1), hi(1)
          do j = lo(2), hi(2)+1
             
             if (nint(idom(i,j-1)).eq.0 .or. nint(idom(i,j)).eq.0) then
                 beta(i,j) = 0_amrex_real ! Suppress flux at the free surface
             else
                temp_face = (temp(i,j) + temp(i,j-1))/2_amrex_real
                call get_conductivity(temp_face, beta(i,j))
             end if
             
          end do
       end do

          
    end if
    
  end subroutine get_beta

  ! -----------------------------------------------------------------
  ! Subroutine used to get the right hand side of the linear solver
  ! -----------------------------------------------------------------  
  subroutine get_rhs(lo, hi, time, dt, &
                     lo_phys, dx, &
                     idom, id_lo, id_hi, &
                     temp, t_lo, t_hi, &
                     rhs, r_lo, r_hi)
  				
    use material_properties_module, only: get_mass_density, get_heat_capacity
    use heat_flux_module, only: get_boundary_heat_flux
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)  
    integer, intent(in) :: t_lo(2), t_hi(2)
    integer, intent(in) :: id_lo(2), id_hi(2)
    integer, intent(in) :: r_lo(2), r_hi(2)    
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2))
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2))
    real(amrex_real), intent(out) :: rhs(r_lo(1):r_hi(1),r_lo(2):r_hi(2))
    real(amrex_real), intent(in) :: lo_phys(2)
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(in) :: time
    real(amrex_real), intent(in) :: dt
    
    ! Local variables
    integer :: i,j
    real(amrex_real) :: rho
    real(amrex_real) :: cp
    
    ! Get the boundary heat flux
    call get_boundary_heat_flux(time, lo_phys, &
                                dx, lo, hi, &
                                idom, id_lo, id_hi, &
                                temp, t_lo, t_hi, rhs)
    
    ! Fill rhs of linear problem
    do i = lo(1), hi(1)
       do j = lo(2), hi(2)
          call get_mass_density(temp(i,j), rho)
          call get_heat_capacity(temp(i,j), cp)
          rhs(i,j) = rhs(i,j)*dt + temp(i,j)*rho*cp
       end do
    end do

  end subroutine get_rhs


  ! -----------------------------------------------------------------
  ! Subroutine used to ...
  ! -----------------------------------------------------------------  
  subroutine get_temp_ls(lo, hi, &
                         temp, t_lo, t_hi, &
                         sol, s_lo, s_hi)
  				
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)  
    integer, intent(in) :: t_lo(2), t_hi(2)
    integer, intent(in) :: s_lo(2), s_hi(2)
    real(amrex_real), intent(out) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2))
    real(amrex_real), intent(in) :: sol(s_lo(1):s_hi(1),s_lo(2):s_hi(2))
    
    ! Local variables
    integer :: i,j
        
    ! Map solution to temperature
    do i = lo(1), hi(1)
       do j = lo(2), hi(2)          
          temp(i,j) = sol(i,j)
       end do
    end do
    
  end subroutine get_temp_ls

  
end module simulation_module
