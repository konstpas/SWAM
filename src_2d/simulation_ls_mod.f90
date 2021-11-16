module simulation_ls_module

  ! -----------------------------------------------------------------
  ! This module is used to run the simulation, WORK IN PROGRESS...
  ! -----------------------------------------------------------------
  
  use amrex_amr_module

  implicit none
  
  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: run_simulation_ls

contains

  ! -----------------------------------------------------------------
  ! Subroutine used to run the simulation
  ! -----------------------------------------------------------------
  subroutine run_simulation_ls()

    use read_input_module, only : max_step, stop_time, plot_int, ls_dt
    use plotfile_module, only: writeplotfile, write1dplotfile
    use energy_module, only: sum_enthalpy

    ! Local variables
    integer :: last_plot_file_step
    integer :: step
    real(amrex_real) :: cur_time
    real(amrex_real) :: total_enthalpy

    ! Initialize time
    cur_time = 0_amrex_real

    ! Initialize counter for output
    last_plot_file_step = 0;

    ! Output initial configuration
    if (plot_int .gt. 0) call writeplotfile

    ! Loop over time-steps
    do step = 0, max_step-1

       ! Print timestep information on screen (start)
       if (amrex_parallel_ioprocessor()) then
          print *, ""
          print *, "STEP", step+1, "starts ..."
       end if
		 
       ! Advance of one timestep
       call advance_one_timestep(cur_time)

       ! Get total enthalpy in the domain
       call sum_enthalpy(total_enthalpy)
 	
       ! Update time
       cur_time = cur_time + ls_dt

       ! Print timestep information on screen (end)
       if (amrex_parallel_ioprocessor()) then
          print *, 'Enthalpy is', total_enthalpy 
          print *, '' 
          print *, "STEP", step+1, "end. TIME =", cur_time, "DT =", ls_dt
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
    
  end subroutine run_simulation_ls
  
  ! -----------------------------------------------------------------
  ! Subroutine used to advance the simulation of one timestep
  ! -----------------------------------------------------------------
  subroutine advance_one_timestep(time)

    !use read_input_module, only : regrid_int
    use amr_data_module, only : phi_old, phi_new

    ! Input and output variables
    real(amrex_real), intent(in) :: time

    ! Local variables
    integer :: ilev
    ! integer, allocatable, save :: last_regrid_step(:)
    ! integer :: finest_level
    ! integer :: fine_substep    
    ! integer :: k
    ! integer :: old_finest_level
    
    ! Regridding 
    ! if (regrid_int .gt. 0) then
       
    !    if (.not.allocated(last_regrid_step)) then

    !       allocate(last_regrid_step(0:amrex_max_level))
    !       last_regrid_step = 0
          
    !    end if

      
    !    if (lev .lt. amrex_max_level .and. stepno(lev) .gt. last_regrid_step(lev)) then

    !       if (mod(stepno(lev), regrid_int) .eq. 0) then

    !          old_finest_level = amrex_get_finest_level()
    !          call amrex_regrid(lev, time)
    !          finest_level = amrex_get_finest_level()

    !          do k = lev, finest_level
    !             last_regrid_step(k) = stepno(k)
    !          end do

    !          do k = old_finest_level+1, finest_level
    !             dt(k) = dt(k-1) / amrex_ref_ratio(k-1)
    !          end do
             
    !       end if
          
    !    end if
       
    ! end if
    
    
    ! Advance solution
    do ilev = 0, amrex_max_level
       call amrex_multifab_swap(phi_old(ilev), phi_new(ilev))       
    end do

    call advance(time) 
    
  end subroutine advance_one_timestep


  ! -----------------------------------------------------------------
  ! Subroutine used to advance the shallow water solver and the
  ! heat equation solver of one time step at a given level
  ! -----------------------------------------------------------------
  subroutine advance(time)

    use read_input_module, only : solve_sw, ls_dt, ls_agglomeration, &
                                  ls_consolidation, ls_max_coarsening_level, &
                                  ls_linop_maxorder, ls_bottom_solver, &
                                  ls_bottom_verbose, ls_max_fmg_iter, &
                                  ls_max_iter, ls_verbose

    use amr_data_module, only : phi_new, temp, idomain, ls_solution, ls_rhs, ls_acoef, ls_bcoef
    use regrid_module, only : fillpatch
    use heat_transfer_module, only : get_idomain, get_melt_pos, reset_melt_pos, increment_enthalpy
    use material_properties_module, only : get_temp
    use shallow_water_module, only : increment_SW
    use amrex_linear_solver_module
    
    ! Input and output variables 
    real(amrex_real), intent(in) :: time

    ! Local variables
    integer, parameter :: nghost = 1 ! number of ghost points in each spatial direction 
    integer :: ncomp
    integer :: idim
    integer :: ilev
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pin
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pout
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptempin
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptemp
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidin
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidout
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pac
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pbc
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: prhs
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: psol
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pu
    type(amrex_multifab) :: phi_tmp(0:amrex_max_level) ! Enthalpy multifab with ghost points
    type(amrex_multifab) :: temp_tmp(0:amrex_max_level) ! Temperature multifab with ghost points
    type(amrex_multifab) :: idomain_tmp(0:amrex_max_level) ! Idomain multifab to distinguish material and background
    type(amrex_mfiter) :: mfi ! Multifab iterator
    type(amrex_box) :: bx
    type(amrex_abeclaplacian) :: abeclap
    type(amrex_boxarray) :: ba(0:amrex_max_level)
    type(amrex_distromap):: dm(0:amrex_max_level)
    type(amrex_multifab) :: nullmf
    type(amrex_multigrid) :: multigrid
    real(amrex_real) :: err
    
    ! ...
    do ilev = 0, amrex_max_level

       ! Build enthalpy and temperature multifabs with ghost points
       ncomp = phi_new(ilev)%ncomp()
       call amrex_multifab_build(phi_tmp(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, ncomp, nghost) 
       call amrex_multifab_build(temp_tmp(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, ncomp, nghost)
       call amrex_multifab_build(idomain_tmp(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, ncomp, nghost)
       
       ! Fill enthalpy multifab
       call fillpatch(ilev, time, phi_tmp(ilev))
       
       ! Swap idomain solution before computing the surface deformation
       call amrex_multifab_swap(idomain_tmp(ilev), idomain(ilev))

       ba(ilev) = phi_new(ilev)%ba
       dm(ilev) = phi_new(ilev)%dm
       
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
                            pidin, lbound(pidin), ubound(pidin), ls_dt)
       
       end do

       ! Clean memory
       call amrex_mfiter_destroy(mfi)    
    
    end if 
    
    ! Set melt interface position array equal to free interface position array 
    ! Since melt layer may span several tile boxes in y-direction (in mfiterator below), we cannot reset within each loop 
    ! Therefore we reset melt position after solving SW, and before propagating temperature 
    ! Melt position is then found after heat has been propagated 
    call reset_melt_pos 
       

    ! Get idomain after the deformation
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
          pac     => ls_acoef(ilev)%dataptr(mfi)
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
          
          ! Get matrix of coefficients (alpha)
          call get_alpha(bx%lo, bx%hi, &
                         ptempin, lbound(ptempin), ubound(ptempin), &
                         pac, lbound(pac), ubound(pac))

          ! Get matrix of coefficients for the gradient term (beta)
          do idim = 1,amrex_spacedim

             pbc => ls_bcoef(idim,ilev)%dataptr(mfi)
             call get_beta(bx%lo, bx%hi, idim, &
                           pidout, lbound(pidout), ubound(pidout), &
                           ptempin, lbound(ptempin), ubound(ptempin), &
                           pbc, lbound(pbc), ubound(pbc))
             
          end do

          ! Get right hand side of the equation
          call get_rhs(bx%lo, bx%hi, time, ls_dt, &
                       amrex_geom(ilev)%get_physical_location(bx%lo), amrex_geom(ilev)%dx, &
                       pidout, lbound(pidout), ubound(pidout), &
                       ptempin, lbound(ptempin), ubound(ptempin), &
                       prhs, lbound(prhs), ubound(prhs))



       end do

       call amrex_mfiter_destroy(mfi)
       
    end do

    ! -----------------------------------------------------------------------------------
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

    call abeclap%set_scalars(1.0_amrex_real, ls_dt) ! A and B scalars
    do ilev = 0, amrex_max_level
       call abeclap%set_acoeffs(ilev, ls_acoef(ilev))
       call abeclap%set_bcoeffs(ilev, ls_bcoef(:,ilev))
    end do

    call amrex_multigrid_build(multigrid, abeclap)
    call multigrid%set_verbose(ls_verbose)
    call multigrid%set_bottom_verbose(ls_bottom_verbose)
    call multigrid%set_max_iter(ls_max_iter)
    call multigrid%set_max_fmg_iter(ls_max_fmg_iter)
    call multigrid%set_bottom_solver(ls_bottom_solver)
    
    err = multigrid%solve(ls_solution, ls_rhs, 1.e-10_amrex_real, 0.0_amrex_real)
    
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
       
       call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.false.)
       do while(mfi%next())

          ! Box
          bx = mfi%validbox()   
          
          ! Pointers
          ptemp   => temp(ilev)%dataptr(mfi)
          pu      => phi_new(ilev)%dataptr(mfi)
          psol    => ls_solution(ilev)%dataptr(mfi)
          
          ! Map solution on temperature
          call get_temp_ls(bx%lo, bx%hi, &
                           ptemp, lbound(ptemp), ubound(ptemp), &
                           psol, lbound(psol), ubound(psol))

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
       call amrex_distromap_destroy(dm)
    end do
    
  end subroutine advance

  
  ! -----------------------------------------------------------------
  ! Subroutine used to ...
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
    
    ! Fill matrix of coefficients
    do i = lo(1), hi(1)
       do j = lo(2), hi(2)          
          
          call get_mass_density(temp(i,j), rho)
          call get_heat_capacity(temp(i,j), cp)
          alpha(i,j) = rho*cp
          
       end do
    end do

  end subroutine get_alpha

  ! -----------------------------------------------------------------
  ! Subroutine used to ...
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
                beta(i,j) = 0_amrex_real
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
                beta(i,j) = 0_amrex_real
             else
                temp_face = (temp(i,j) + temp(i,j-1))/2_amrex_real
                call get_conductivity(temp_face, beta(i,j))
             end if
             
          end do
       end do

          
    end if
    
  end subroutine get_beta

  ! -----------------------------------------------------------------
  ! Subroutine used to ...
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
    
    ! Get the boundary heat flux
    call get_boundary_heat_flux(time, lo_phys, &
                                dx, lo, hi, &
                                idom, id_lo, id_hi, &
                                temp, t_lo, t_hi, rhs)
    
    ! Fill rhs of linear problem
    do i = lo(1), hi(1)
       do j = lo(2), hi(2)          
          rhs(i,j) = rhs(i,j)*dt + temp(i,j)          
       end do
    end do

  end subroutine get_rhs


  ! -----------------------------------------------------------------
  ! Subroutine used to ...
  ! -----------------------------------------------------------------  
  subroutine get_temp_ls(lo, hi, &
                         temp, t_lo, t_hi, &
                         sol, s_lo, s_hi)
  				
    use material_properties_module, only: get_mass_density, get_heat_capacity
    use heat_flux_module, only: get_boundary_heat_flux
    
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

  
end module simulation_ls_module
