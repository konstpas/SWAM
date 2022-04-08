module heat_transfer_implicit_module
  
  ! -----------------------------------------------------------------
  ! This module contains the implicit solver for the heat transfer
  ! equations
  ! -----------------------------------------------------------------
  
  use amrex_amr_module

  implicit none

  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: advance_heat_solver_implicit
  
contains

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy at a new time step for
  ! all levels via an implicit update.
  ! -----------------------------------------------------------------
  subroutine advance_heat_solver_implicit(time, dt)
    
    ! Input and output variables
    real(amrex_real), intent(in) :: dt
    real(amrex_real), intent(in) :: time
    
    ! Local variables
    integer :: ilev
    type(amrex_multifab) :: phi_tmp(0:amrex_max_level) ! Enthalpy multifab with ghost points
    type(amrex_multifab) :: temp_tmp(0:amrex_max_level) ! Temperature multifab with ghost points
    type(amrex_multifab) :: idomain_tmp(0:amrex_max_level) ! Idomain multifab to distinguish material and background

    ! Initialize temporary multifabs
    call init_tmp_multifab(time, phi_tmp, temp_tmp, idomain_tmp)

    ! Advance the conductive part of the heat equation
    call advance_conduction(time, dt, phi_tmp, temp_tmp, idomain_tmp)

    ! Advance the advective part of the heat equation
    call advance_advection(dt, temp_tmp)

    ! Update position of the melt bottom
    call get_melt_pool_bottom   
    
    ! Clean memory
    do ilev = 0, amrex_max_level
       call amrex_multifab_destroy(phi_tmp(ilev))
       call amrex_multifab_destroy(temp_tmp(ilev))
       call amrex_multifab_destroy(idomain_tmp(ilev))
    end do

    
  end subroutine advance_heat_solver_implicit

  ! -----------------------------------------------------------------
  ! Subroutine used to initialize the multifabs used in the
  ! update of the full heat equation
  ! -----------------------------------------------------------------
  subroutine init_tmp_multifab(time, phi_tmp, temp_tmp, idomain_tmp)
    
    use amr_data_module, only : phi_new, &
                                phi_old, &
                                idomain
    use regrid_module, only : fillpatch
 
    ! Input and output variables
    real(amrex_real), intent(in) :: time
    type(amrex_multifab), intent(out) :: phi_tmp(0:amrex_max_level)
    type(amrex_multifab), intent(out) :: temp_tmp(0:amrex_max_level)
    type(amrex_multifab), intent(out) :: idomain_tmp(0:amrex_max_level)
    
    ! Local variables
    integer, parameter :: nghost = 1 ! number of ghost points in each spatial direction
    integer :: ncomp
    integer :: ilev
    type(amrex_boxarray) :: ba(0:amrex_max_level) ! Box array
    type(amrex_distromap):: dm(0:amrex_max_level) ! Distribution mapping
    type(amrex_geometry) :: geom

    do ilev = 0, amrex_max_level
       
       ! Swap old and new enthalpies
       call amrex_multifab_swap(phi_old(ilev), phi_new(ilev))
       
       ! Get boxarray and distribution mapping on level
       ba(ilev) = phi_new(ilev)%ba
       dm(ilev) = phi_new(ilev)%dm

       ! Number of components
       ncomp = phi_new(ilev)%ncomp()

       ! Geometry
       geom = amrex_geom(ilev)
       
       ! Enthalpy and temperature multifabs with ghost points
       call amrex_multifab_build(phi_tmp(ilev), ba(ilev), dm(ilev), ncomp, nghost) 
       call amrex_multifab_build(temp_tmp(ilev), ba(ilev), dm(ilev), ncomp, nghost)
       call fillpatch(ilev, time, phi_tmp(ilev))

       ! Build temporary idomain multifab to store the domain configuration before SW deformation
       call amrex_multifab_build(idomain_tmp(ilev), ba(ilev), dm(ilev), ncomp, nghost)
       call amrex_multifab_swap(idomain_tmp(ilev), idomain(ilev))
       
    end do

    ! Clean memory
    do ilev = 0, amrex_max_level 
       call amrex_distromap_destroy(dm(ilev))
    end do
    
  end subroutine init_tmp_multifab
  
  ! -----------------------------------------------------------------
  ! Subroutine used to advance of one time step the conduction part
  ! of the heat transfer equation
  ! -----------------------------------------------------------------
  subroutine advance_conduction(time, dt, phi_tmp, temp_tmp, idomain_tmp)
    
    ! Input and output variables
    real(amrex_real), intent(in) :: dt
    real(amrex_real), intent(in) :: time
    type(amrex_multifab), intent(inout) :: phi_tmp(0:amrex_max_level)
    type(amrex_multifab), intent(inout) :: temp_tmp(0:amrex_max_level)
    type(amrex_multifab), intent(inout) :: idomain_tmp(0:amrex_max_level) 
    
    ! Local variables
    integer :: idim
    integer :: ilev
    type(amrex_multifab) :: rhs(0:amrex_max_level) ! multifab for the right hand side of the linear system of equations
    type(amrex_multifab) :: alpha(0:amrex_max_level) ! multifab for the first term on the left hand side of the linear system of equations
    type(amrex_multifab) :: beta(amrex_spacedim,0:amrex_max_level) ! multifab for the second term on the left hand side of the linear system of equations
    type(amrex_boxarray) :: ba(0:amrex_max_level) ! Box array
    type(amrex_distromap):: dm(0:amrex_max_level) ! Distribution mapping

    ! Initialize temporary multifabs
    call conduction_init_tmp_multifab(ba, dm, alpha, beta, rhs)

    ! Predictor step on all levels
    call conduction_predict(time, dt, phi_tmp, temp_tmp, idomain_tmp, alpha, &
                            beta, rhs)    
    
    ! New temperature via implicit update
    call conduction_get_temperature(ba, dm, dt, alpha, beta, rhs, temp_tmp)
    
    ! Corrector step on all levels
    call conduction_correct(phi_tmp, temp_tmp, alpha)

    ! Synchronize idomains
    call conduction_synch_idomain(temp_tmp)

    ! Clean memory
    do ilev = 0, amrex_max_level
       call amrex_multifab_destroy(rhs(ilev))
       call amrex_multifab_destroy(alpha(ilev))
       do idim = 1, amrex_spacedim
          call amrex_multifab_destroy(beta(idim,ilev))
       end do
       call amrex_distromap_destroy(dm(ilev))
    end do
   
  end subroutine advance_conduction

  ! -----------------------------------------------------------------
  ! Subroutine used to initialize the multifabs used in the
  ! update of the conductive part of the heat equation
  ! -----------------------------------------------------------------
  subroutine conduction_init_tmp_multifab(ba, dm, alpha, beta, rhs)
    
    use amr_data_module, only : phi_new

    ! Input and output variables
    type(amrex_multifab), intent(out) :: rhs(0:amrex_max_level)
    type(amrex_multifab), intent(out) :: alpha(0:amrex_max_level)
    type(amrex_multifab), intent(out) :: beta(amrex_spacedim,0:amrex_max_level)
    type(amrex_boxarray), intent(out) :: ba(0:amrex_max_level)
    type(amrex_distromap), intent(out) :: dm(0:amrex_max_level)
    
    ! Local variables
    logical :: nodal(2)
    integer :: idim
    integer :: ilev
    integer :: ncomp
    type(amrex_geometry) :: geom

    ! Initialize temporary multifabs
    do ilev = 0, amrex_max_level
       
       ! Get boxarray and distribution mapping
       ba(ilev) = phi_new(ilev)%ba
       dm(ilev) = phi_new(ilev)%dm

       ! Number of components
       ncomp = phi_new(ilev)%ncomp()

       ! Geometry
       geom = amrex_geom(ilev)

       ! Multifabs for the linear solver
       call amrex_multifab_build(rhs(ilev), ba(ilev), dm(ilev), ncomp, 0)
       call amrex_multifab_build(alpha(ilev), ba(ilev), dm(ilev), ncomp, 0)
       do idim = 1, amrex_spacedim
          nodal = .false.
          nodal(idim) = .true.
          call amrex_multifab_build(beta(idim,ilev), ba(ilev), dm(ilev), ncomp, 0, nodal)
       end do

    end do
    
  end subroutine conduction_init_tmp_multifab


  ! -----------------------------------------------------------------
  ! Subroutine used to compute the predictor step of the conductive
  ! part of the heat transfer equation
  ! -----------------------------------------------------------------
  subroutine conduction_predict(time, dt, phi_tmp, temp_tmp, &
                                idomain_tmp, alpha, beta, rhs)
    
    use amr_data_module, only : phi_new
    
    ! Input and output variables
    real(amrex_real), intent(in) :: dt
    real(amrex_real), intent(in) :: time
    type(amrex_multifab), intent(inout) :: phi_tmp(0:amrex_max_level)
    type(amrex_multifab), intent(inout) :: temp_tmp(0:amrex_max_level)
    type(amrex_multifab), intent(inout) :: idomain_tmp(0:amrex_max_level) 
    type(amrex_multifab), intent(inout) :: rhs(0:amrex_max_level)
    type(amrex_multifab), intent(inout) :: alpha(0:amrex_max_level)
    type(amrex_multifab), intent(inout) :: beta(amrex_spacedim,0:amrex_max_level)
    
    ! Local variables
    integer :: ilev
    type(amrex_geometry) :: geom
    type(amrex_mfiter) :: mfi

    ! Predictor step on all levels
    do ilev = 0, amrex_max_level

       ! Geometry
       geom = amrex_geom(ilev)
       
       !$omp parallel private(mfi)
       call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.false.)
       do while(mfi%next())
          call conduction_predict_box(ilev, time, dt, mfi, &
                                      geom, phi_tmp(ilev), temp_tmp(ilev), &
                                      idomain_tmp(ilev), alpha(ilev), beta(:,ilev), &
                                      rhs(ilev))

       end do
       call amrex_mfiter_destroy(mfi)
       !$omp end parallel

    end do
 
  end subroutine conduction_predict
  
  ! -----------------------------------------------------------------
  ! Subroutine used to perform the prediction step for the implicit
  ! enthalpy update at a given box on a given level 
  ! -----------------------------------------------------------------
  subroutine conduction_predict_box(lev, time, dt, mfi, &
                                    geom, phi_tmp, temp_tmp, &
                                    idomain_tmp, alpha, beta, &
                                    rhs)

    use read_input_module, only : temp_fs
    use amr_data_module, only : phi_new, temp, idomain
    use heat_transfer_domain_module, only : get_idomain, &
                                            revaluate_heat_domain
    use material_properties_module, only : get_temp
    
    ! Input and output variables
    integer, intent(in) :: lev
    real(amrex_real), intent(in) :: dt
    real(amrex_real), intent(in) :: time
    type(amrex_mfiter), intent(in) :: mfi
    type(amrex_geometry), intent(in) :: geom
    type(amrex_multifab), intent(inout) :: phi_tmp
    type(amrex_multifab), intent(inout) :: temp_tmp
    type(amrex_multifab), intent(inout) :: idomain_tmp
    type(amrex_multifab), intent(inout) :: alpha
    type(amrex_multifab), intent(inout) :: beta(amrex_spacedim)
    type(amrex_multifab), intent(inout) :: rhs

    ! Local variables
    integer :: idim
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pin
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pout
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptempin
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptemp
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidin
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidout
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pac
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pbc
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: prhs
    type(amrex_box) :: bx
    
    ! Box
    bx = mfi%validbox()   
          
    ! Pointers
    pin     => phi_tmp%dataptr(mfi)
    pout    => phi_new(lev)%dataptr(mfi)
    ptemp   => temp(lev)%dataptr(mfi)
    ptempin => temp_tmp%dataptr(mfi)
    pidin   => idomain_tmp%dataptr(mfi)
    pidout  => idomain(lev)%dataptr(mfi)
    pac     => alpha%dataptr(mfi)
    prhs    => rhs%dataptr(mfi)
    
    ! Get temperature corresponding to the enthalpy (with ghost points)
    call get_temp(lbound(ptempin), ubound(ptempin), &
                  pin, lbound(pin), ubound(pin), &
                  ptempin, lbound(ptempin), ubound(ptempin), .true.)
    
    ! Get configuration of the system after the deformation
    call get_idomain(geom%get_physical_location(bx%lo), geom%dx, &
                     bx%lo, bx%hi, &
                     pidout, lbound(pidout), ubound(pidout), &
                     ptempin, lbound(ptempin), ubound(ptempin))
    
    ! Re-evaluate heat transfer domain after the deformation
    call revaluate_heat_domain(geom%get_physical_location(bx%lo), geom%dx, &
                               bx%lo, bx%hi, &
                               pidin, lbound(pidin), ubound(pidin), &
                               pidout, lbound(pidout), ubound(pidout), &
                               pin, lbound(pin), ubound(pin), &
                               ptempin, lbound(ptempin), ubound(ptempin))
    
    ! Additional call to update the temperature without the ghost points
    ! the temperature with ghost points is already updated in the
    ! subroutine revaluate_heat_domain
    call get_temp(lbound(ptemp), ubound(ptemp), &
                  pin, lbound(pin), ubound(pin), &
                  ptemp, lbound(ptemp), ubound(ptemp), .true.)
    
    ! Get alpha matrix for the linear solver (first term on left hand side)
    call get_alpha(bx%lo, bx%hi, &
                   pin, lbound(pin), ubound(pin), &
                   ptempin, lbound(ptempin), ubound(ptempin), &                       
                   pac, lbound(pac), ubound(pac))
    

    ! Get beta matrix for the linear solver (second term on left hand side)
    do idim = 1,amrex_spacedim
       pbc => beta(idim)%dataptr(mfi)
       call get_beta(bx%lo, bx%hi, idim, &
                     pidout, lbound(pidout), ubound(pidout), &
                     ptempin, lbound(ptempin), ubound(ptempin), &
                     pbc, lbound(pbc), ubound(pbc))
    end do
    
    ! Get right hand for the linear solver
    if (temp_fs.gt.0) then
       call get_rhs_fixT(bx%lo, bx%hi, &
                         ptempin, lbound(ptempin), ubound(ptempin), &
                         pac, lbound(pac), ubound(pac), &                      
                         prhs, lbound(prhs), ubound(prhs))
    else
       call get_rhs(bx%lo, bx%hi, time, dt, &
                    geom%get_physical_location(bx%lo), geom%dx, &
                    pidout, lbound(pidout), ubound(pidout), &
                    ptempin, lbound(ptempin), ubound(ptempin), &
                    pac, lbound(pac), ubound(pac), &                      
                    prhs, lbound(prhs), ubound(prhs), lev)
    end if
    
  end subroutine conduction_predict_box
  
  
  ! -----------------------------------------------------------------------------------
  ! Subroutine used to solve the linear system of equations for the implict update
  ! of the temperature
  ! ----------------------------------------------------------------------------------- 
  subroutine conduction_get_temperature(ba, dm, dt, alpha, beta, rhs, sol)

    use read_input_module, only : ls_accuracy, &
                                  ls_composite_solve, &
                                  ls_agglomeration, &
                                  ls_consolidation, &
                                  ls_max_coarsening_level, &
                                  ls_linop_maxorder, &
                                  ls_bottom_solver, &
                                  ls_bottom_verbose, &
                                  ls_max_fmg_iter, &
                                  ls_max_iter, &
                                  ls_verbose
    use amrex_linear_solver_module

        
    ! Input and output variables
    type(amrex_multifab), intent(inout) :: sol(0:amrex_max_level)
    type(amrex_multifab), intent(in) :: rhs(0:amrex_max_level)
    type(amrex_multifab), intent(in) :: alpha(0:amrex_max_level)
    type(amrex_multifab), intent(in) :: beta(amrex_spacedim,0:amrex_max_level)
    type(amrex_boxarray), intent(in) :: ba(0:amrex_max_level)
    type(amrex_distromap), intent(in) :: dm(0:amrex_max_level)
    real(amrex_real), intent(in) :: dt
    
    ! Local variables
    integer ilev
    real(amrex_real) :: err ! Error from the solution of the linear system of equations
    type(amrex_multifab) :: nullmf ! Null multifab used to define the boundary condition at the coarse-fine interface
    type(amrex_abeclaplacian) :: abeclap ! Object for the solution of the linear systems of equations
    type(amrex_multigrid) :: multigrid ! Object for the solution of the linear systems of equations


    ! Solve linear system of equations
    if (ls_composite_solve) then
       
       call amrex_abeclaplacian_build(abeclap, amrex_geom, ba, dm, &
                                      metric_term=.false., agglomeration=ls_agglomeration, &
                                      consolidation=ls_consolidation, &
                                      max_coarsening_level=ls_max_coarsening_level)

       call abeclap%set_maxorder(ls_linop_maxorder)
       
       ! This is set up to have homogeneous Neumann BC
       call abeclap%set_domain_bc([amrex_lo_neumann, amrex_lo_neumann], &
                                  [amrex_lo_neumann, amrex_lo_neumann])
       
       ! Coarse-fine interface  (for problem with pure homogeneous Neumann BC we can pass an empty multifab)
       do ilev = 0, amrex_max_level
          call abeclap%set_level_bc(ilev, nullmf)
       end do

       ! Set A and B scalar
       call abeclap%set_scalars(1.0_amrex_real, dt)
       
       ! Set alpha and beta matrices
       do ilev = 0, amrex_max_level
          call abeclap%set_acoeffs(ilev, alpha(ilev))
          call abeclap%set_bcoeffs(ilev, beta(:,ilev))
       end do
       
       ! Build multigrid solver
       call amrex_multigrid_build(multigrid, abeclap)
       call multigrid%set_verbose(ls_verbose)
       call multigrid%set_bottom_verbose(ls_bottom_verbose)
       call multigrid%set_max_iter(ls_max_iter)
       call multigrid%set_max_fmg_iter(ls_max_fmg_iter)
       call multigrid%set_bottom_solver(ls_bottom_solver)
       
       ! Solve the linear system
       err = multigrid%solve(sol, rhs, ls_accuracy, 0.0_amrex_real)
       
       ! Clean memory
       call amrex_abeclaplacian_destroy(abeclap)
       call amrex_multigrid_destroy(multigrid)

    else

       do ilev = 0, amrex_max_level

          call amrex_abeclaplacian_build(abeclap, [amrex_geom(ilev)], [ba(ilev)], [dm(ilev)], &
                                         metric_term=.false., agglomeration=ls_agglomeration, &
                                         consolidation=ls_consolidation, &
                                        max_coarsening_level=ls_max_coarsening_level)
          
          call abeclap%set_maxorder(ls_linop_maxorder)

          ! This is set up to have homogeneous Neumann BC
          call abeclap%set_domain_bc([amrex_lo_neumann, amrex_lo_neumann], &
                                     [amrex_lo_neumann, amrex_lo_neumann])
       
          ! Coarse-fine interface  (for problem with pure homogeneous Neumann BC we can pass an empty multifab)
          if (ilev > 0) then
             ! use coarse level data to set up bc at corase/fine boundary
             call abeclap % set_coarse_fine_bc(sol(ilev-1), amrex_ref_ratio(ilev-1))
          end if
          call abeclap%set_level_bc(0, nullmf)

          ! Set A and B scalar
          call abeclap%set_scalars(1.0_amrex_real, dt)
          
          ! Set alpha and beta matrices
          call abeclap%set_acoeffs(0, alpha(ilev))
          call abeclap%set_bcoeffs(0, beta(:,ilev))

          ! Build multigrid solver
          call amrex_multigrid_build(multigrid, abeclap)
          call multigrid%set_verbose(ls_verbose)
          call multigrid%set_bottom_verbose(ls_bottom_verbose)
          call multigrid%set_max_iter(ls_max_iter)
          call multigrid%set_max_fmg_iter(ls_max_fmg_iter)
          call multigrid%set_bottom_solver(ls_bottom_solver)
        
          ! Solve the linear system
          err = multigrid%solve([sol(ilev)], [rhs(ilev)], ls_accuracy, 0.0_amrex_real)
       
          ! Clean memory
          call amrex_abeclaplacian_destroy(abeclap)
          call amrex_multigrid_destroy(multigrid)
          
       end do
       
    end if
  end subroutine conduction_get_temperature

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the correction step of the conductive
  ! part of the heat transfer equation
  ! -----------------------------------------------------------------
  subroutine conduction_correct(phi_tmp, temp_tmp, alpha)
    
    use amr_data_module, only : phi_new
    
    ! Input and output variables
    type(amrex_multifab), intent(inout) :: phi_tmp(0:amrex_max_level)
    type(amrex_multifab), intent(inout) :: temp_tmp(0:amrex_max_level)
    type(amrex_multifab), intent(inout) :: alpha(0:amrex_max_level)
    
    ! Local variables
    integer :: ilev
    type(amrex_geometry) :: geom
    type(amrex_mfiter) :: mfi

    ! Corrector step on all levels
    do ilev = 0, amrex_max_level

       ! Geometry
       geom = amrex_geom(ilev)
       
       !$omp parallel private(mfi)
       call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.false.)
       do while(mfi%next())
          call conduction_correct_box(ilev, mfi, phi_tmp(ilev), &
                                      temp_tmp(ilev), alpha(ilev))
          
       end do
       call amrex_mfiter_destroy(mfi)
       !$omp end parallel
       
       call temp_tmp(ilev)%fill_boundary(geom)
       
    end do
    
  end subroutine conduction_correct
  
  ! -----------------------------------------------------------------
  ! Subroutine used to perform the correction step for the implicit
  ! enthalpy update at a given box on a given level 
  ! -----------------------------------------------------------------
  subroutine conduction_correct_box(lev, mfi, phi_tmp, &
                                    temp_tmp, alpha)

    use amr_data_module, only : phi_new, temp, idomain
    use heat_transfer_domain_module, only : get_melt_pos, get_idomain
    use material_properties_module, only : get_temp
    use read_input_module, only : temp_fs
    
    ! Input and output variables
    integer, intent(in) :: lev
    type(amrex_mfiter), intent(in) :: mfi
    type(amrex_multifab), intent(in) :: phi_tmp
    type(amrex_multifab), intent(inout) :: temp_tmp
    type(amrex_multifab), intent(inout) :: alpha

    ! Local variables
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pin
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pout
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptempin
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptemp
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pid
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pac
    type(amrex_box) :: bx
    
    ! Box
    bx = mfi%validbox()   
    
    ! Pointers
    pin     => phi_tmp%dataptr(mfi)
    pout    => phi_new(lev)%dataptr(mfi)
    ptempin => temp_tmp%dataptr(mfi)
    ptemp   => temp(lev)%dataptr(mfi)
    pid     => idomain(lev)%dataptr(mfi)
    pac => alpha%dataptr(mfi)
    
    ! Get corrected enthalpy
    if (temp_fs.gt.0) then
      STOP 'The implicit solver cannot be used to solve problems with a fixed temperature on the free surface'
    else
       call get_enthalpy_implicit(bx%lo, bx%hi, &
                                  pin, lbound(pin), ubound(pin), &
                                  pout, lbound(pout), ubound(pout), &
                                  ptemp, lbound(ptemp), ubound(ptemp), &
                                  ptempin, lbound(ptempin), ubound(ptempin), &
                                  pac, lbound(pac), ubound(pac))                     
    end if

  end subroutine conduction_correct_box

  ! -----------------------------------------------------------------
  ! Subroutine used to synchronize the idomains in the conduction 
  ! part of the heat transfer equation
  ! -----------------------------------------------------------------
  subroutine conduction_synch_idomain(temp_tmp)
    
    use amr_data_module, only : phi_new, &
                                idomain
    use heat_transfer_domain_module, only : get_idomain
        
    ! Input and output variables
    type(amrex_multifab), intent(inout) :: temp_tmp(0:amrex_max_level)
    
    ! Local variables
    integer :: ilev
    type(amrex_geometry) :: geom
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidom
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptemp_tmp
    
    do ilev = 0, amrex_max_level

       ! Geometry
       geom = amrex_geom(ilev)
       
       !$omp parallel private(mfi, bx, pidom, ptemp_tmp)
       call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.false.)
       do while(mfi%next())
          bx = mfi%validbox()
          pidom  => idomain(ilev)%dataptr(mfi)
          ptemp_tmp   => temp_tmp(ilev)%dataptr(mfi)
          call get_idomain(geom%get_physical_location(bx%lo), geom%dx, &
                           bx%lo, bx%hi, &
                           pidom, lbound(pidom), ubound(pidom), &
                           ptemp_tmp, lbound(ptemp_tmp), ubound(ptemp_tmp))
          
       end do
       call amrex_mfiter_destroy(mfi)   
       !$omp end parallel
       
    end do
    
  end subroutine conduction_synch_idomain
  
  ! -----------------------------------------------------------------
  ! Subroutine used to fill the alpha matrix for the linear solver
  ! -----------------------------------------------------------------  
  subroutine get_alpha(lo, hi, &
                       u_old, uo_lo, uo_hi, &
                       temp, t_lo, t_hi, &
                       alpha, a_lo, a_hi)

    use material_properties_module, only: get_mass_density, get_heat_capacity, &
                                          enth_at_melt, temp_melt, latent_heat
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: uo_lo(2), uo_hi(2)
    integer, intent(in) :: t_lo(2), t_hi(2)
    integer, intent(in) :: a_lo(2), a_hi(2)
    real(amrex_real), intent(in) :: u_old(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2))  
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2))
    real(amrex_real), intent(out) :: alpha(a_lo(1):a_hi(1),a_lo(2):a_hi(2))

    ! Local variables
    integer :: i,j
    real(amrex_real) :: cp
    real(amrex_real) :: cps
    real(amrex_real) :: cpl
    real(amrex_real) :: rho
    real(amrex_real) :: rhos
    real(amrex_real) :: rhol
    real(amrex_real) :: lf

    ! Fill alpha matrix
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (temp(i,j).ne.temp_melt) then
             call get_mass_density(temp(i,j), rho)
             call get_heat_capacity(temp(i,j), cp)
             alpha(i,j) = rho*cp
          else
             ! Compute liquid fraction
             lf = (u_old(i,j) - enth_at_melt)/latent_heat
             if (lf .gt. 1.0) lf = 1.0
             if (lf .lt. 0.0) lf = 0.0
             ! Get solid and liquid properties close to the melt transition
             call get_mass_density(temp_melt - 1e-10_amrex_real, rhos)
             call get_mass_density(temp_melt, rhol)
             call get_heat_capacity(temp_melt - 1e-10_amrex_real, cps)
             call get_heat_capacity(temp_melt, cpl)
             ! Properties of the liquid at the melting point
             cp = (1-lf)*cps + lf*cpl
             rho = (1-lf)*rhos + lf*rhol
             alpha(i,j) = rho*cp
          end if
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
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             
             if(nint(idom(i-1,j)).eq.0 .or. nint(idom(i,j)).eq.0 .or. & 
                nint(idom(i-1,j)).eq.-1 .or. nint(idom(i,j)).eq.-1) then
                 beta(i,j) = 0_amrex_real ! Suppress flux at the free surface and cooling pipe surface
              else
                temp_face = (temp(i,j) + temp(i-1,j))/2_amrex_real
                call get_conductivity(temp_face, beta(i,j))
             end if

          end do
       end do

    else if (dim == 2) then ! y-direction
       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
        
             if(nint(idom(i,j-1)).eq.0 .or. nint(idom(i,j)).eq.0 .or. &
                nint(idom(i,j-1)).eq.-1 .or. nint(idom(i,j)).eq.-1) then
                 beta(i,j) = 0_amrex_real ! Suppress flux at the free surface and cooling pipe surface
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
                     alpha, a_lo, a_hi, &
                     rhs, r_lo, r_hi, lev)

    use heat_transfer_flux_module, only: get_boundary_heat_flux
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: id_lo(2), id_hi(2)
    integer, intent(in) :: t_lo(2), t_hi(2)
    integer, intent(in) :: r_lo(2), r_hi(2)
    integer, intent(in) :: a_lo(2), a_hi(2)
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2))
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2))
    real(amrex_real), intent(in) :: alpha(a_lo(1):a_hi(1),a_lo(2):a_hi(2))   
    real(amrex_real), intent(out) :: rhs(r_lo(1):r_hi(1),r_lo(2):r_hi(2))
    real(amrex_real), intent(in) :: lo_phys(2)
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(in) :: time
    real(amrex_real), intent(in) :: dt
    integer, intent(in) :: lev
    
    ! Local variables
    integer :: i,j

    ! Get the boundary heat flux
    call get_boundary_heat_flux(time, lo_phys, &
                                dx, lo, hi, &
                                idom, id_lo, id_hi, &
                                temp, t_lo, t_hi, lev, rhs)
    
    ! Fill rhs of linear problem
    do j = lo(2), hi(2)          
       do i = lo(1), hi(1)
          rhs(i,j) = rhs(i,j)*dt + temp(i,j)*alpha(i,j)
       end do
    end do

  end subroutine get_rhs

  
  ! -----------------------------------------------------------------
  ! Subroutine used to get the right hand side of the linear solver.
  ! Case with fixed temperature at the free surface
  ! -----------------------------------------------------------------  
  subroutine get_rhs_fixT(lo, hi, &
                          temp, t_lo, t_hi, &
                          alpha, a_lo, a_hi, &
                          rhs, r_lo, r_hi)
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: t_lo(2), t_hi(2)
    integer, intent(in) :: r_lo(2), r_hi(2)
    integer, intent(in) :: a_lo(2), a_hi(2)
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2))
    real(amrex_real), intent(in) :: alpha(a_lo(1):a_hi(1),a_lo(2):a_hi(2))   
    real(amrex_real), intent(out) :: rhs(r_lo(1):r_hi(1),r_lo(2):r_hi(2))
    
    ! Local variables
    integer :: i,j

    ! Fill rhs of linear problem
    do j = lo(2), hi(2)          
       do i = lo(1), hi(1)
          rhs(i,j) = temp(i,j)*alpha(i,j)
       end do
    end do

  end subroutine get_rhs_fixT

  ! -----------------------------------------------------------------
  ! Subroutine used to get the corrected enthalpy
  ! -----------------------------------------------------------------  
  subroutine get_enthalpy_implicit(lo, hi, &
                                   u_old, uo_lo, uo_hi, &
                                   u_new, un_lo, un_hi, &
                                   temp_old, to_lo, to_hi, &
                                   temp_new, tn_lo, tn_hi, &
                                   alpha, a_lo, a_hi)

    use material_properties_module, only : get_temp, &
                                           temp_melt, &
                                           enth_at_melt, &
                                           latent_heat
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)  
    integer, intent(in) :: uo_lo(2), uo_hi(2)
    integer, intent(in) :: un_lo(2), un_hi(2)
    integer, intent(in) :: to_lo(2), to_hi(2)
    integer, intent(in) :: tn_lo(2), tn_hi(2)
    integer, intent(in) :: a_lo(2), a_hi(2)
    real(amrex_real), intent(in) :: u_old(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2))
    real(amrex_real), intent(out) :: u_new(un_lo(1):un_hi(1),un_lo(2):un_hi(2))
    real(amrex_real), intent(inout) :: temp_old(to_lo(1):to_hi(1),to_lo(2):to_hi(2))
    real(amrex_real), intent(inout) :: temp_new(tn_lo(1):tn_hi(1),tn_lo(2):tn_hi(2))
    real(amrex_real), intent(in) :: alpha(a_lo(1):a_hi(1),a_lo(2):a_hi(2))
    
    ! Local variables
    integer :: i,j
    real(amrex_real) :: ubase
    
    ! General case, when no phase transition
    ! (will be overwritten for phase transitioning control volumes.)
    call get_temp(un_lo, un_hi, &
                  u_new, un_lo, un_hi, &
                  temp_new, tn_lo, tn_hi, .false.)
    
    ! Get new enthalpy for points undergoing transition 
    ! (might overwrite enthalpy found previously.)
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (temp_old(i,j).le.temp_melt .and. temp_new(i,j).gt.temp_melt) then
             ubase = max(u_old(i,j), enth_at_melt)
             u_new(i,j) = ubase + alpha(i,j)*(temp_new(i,j) - temp_melt)
          elseif (temp_old(i,j).ge.temp_melt .and. temp_new(i,j).lt.temp_melt) then
             ubase = min(u_old(i,j), enth_at_melt+latent_heat)
             u_new(i,j) = ubase + alpha(i,j)*(temp_new(i,j) - temp_melt)
          end if
       end do
    end do
    
    ! Update temperature acordingly (with ghost points)
    call get_temp(lo, hi, &
         u_new, un_lo, un_hi, &
         temp_new, tn_lo, tn_hi, .true.)
    
    ! Update temperature acordingly (without ghost points)
    call get_temp(lo, hi, &
         u_new, un_lo, un_hi, &
         temp_old, to_lo, to_hi, .true.)
      
  end subroutine get_enthalpy_implicit
  
  ! -----------------------------------------------------------------
  ! Subroutine used to advance of one time step the advective part
  ! of the heat transfer equation
  ! -----------------------------------------------------------------
  subroutine advance_advection(dt, temp_tmp)
    
    use amr_data_module, only : phi_new, &
                                temp, &
                                idomain

    ! Input and output variables
    real(amrex_real), intent(in) :: dt
    type(amrex_multifab), intent(inout) :: temp_tmp(0:amrex_max_level)
    
    ! Local variables
    integer :: ilev
    type(amrex_geometry) :: geom
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidom
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptemp_tmp
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptemp
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pout
   
    do ilev = 0, amrex_max_level
 
       ! Geometry
       geom = amrex_geom(ilev)
       
       !$omp parallel private(mfi, bx, pidom, ptemp, ptemp_tmp, pout)
       call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.false.)
       do while(mfi%next())
          bx = mfi%validbox()
          pidom  => idomain(ilev)%dataptr(mfi)
          ptemp   => temp(ilev)%dataptr(mfi)
          ptemp_tmp   => temp_tmp(ilev)%dataptr(mfi)
          pout    => phi_new(ilev)%dataptr(mfi)
          call advection_box(geom%get_physical_location(bx%lo), bx%lo, &
                             bx%hi, dt, geom%dx(1:3), &
                             lbound(pidom), ubound(pidom), pidom, &
                             pout, lbound(pout), ubound(pout), &
                             ptemp_tmp, lbound(ptemp_tmp), ubound(ptemp_tmp), &
                             ptemp, lbound(ptemp), ubound(ptemp))
          
       end do
       call amrex_mfiter_destroy(mfi)
       !$omp end parallel   
       
    end do
    
  end subroutine advance_advection
  
  ! -----------------------------------------------------------------  
  ! Subroutine used for an explicit update of the temperature
  ! equation only with heat advection and no heat conduction.
  ! -----------------------------------------------------------------  
  subroutine advection_box(lo_phys, lo, hi, dt, dx, &
                           id_lo, id_hi, idom, &
                           u_new, un_lo, un_hi, &
                           temp_old, to_lo, to_hi, &
                           temp, t_lo, t_hi)

    use material_properties_module, only : get_temp, get_enthalpy
    use heat_transfer_domain_module, only: get_face_velocity
    
    ! Input and output variables
    integer, intent(in) :: lo(2)
    integer, intent(in) :: hi(2)
    integer, intent(in) :: t_lo(2)
    integer, intent(in) :: t_hi(2)
    integer, intent(in) :: to_lo(2)
    integer, intent(in) :: to_hi(2)
    integer, intent(in) :: id_lo(2)
    integer, intent(in) :: id_hi(2)
    integer, intent(in) :: un_lo(2)
    integer, intent(in) :: un_hi(2)
    real(amrex_real), intent(in) :: lo_phys(2)
    real(amrex_real), intent(in) :: dt
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2))
    real(amrex_real), intent(inout) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2))
    real(amrex_real), intent(in) :: temp_old(to_lo(1):to_hi(1),to_lo(2):to_hi(2))
    real(amrex_real), intent(inout) :: u_new(un_lo(1):un_hi(1),un_lo(2):un_hi(2))
    
    ! Local variables
    integer :: i,j
    integer :: ux_lo(2), ux_hi(2)
    real(amrex_real) :: ux(lo(1):hi(1)+1,lo(2):hi(2))
    real(amrex_real) :: vx_l, vx_r
    
    ux_lo = lo
    ux_hi(1) = hi(1)+1
    ux_hi(2) = hi(2)
    
    ! Construct 2D melt velocity profile from the 2D shallow water solution
    call get_face_velocity(lo_phys, lo, hi, dx, &
                           ux, ux_lo, ux_hi, &
                           idom, id_lo, id_hi)

    ! Update temperature profile
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          vx_l = ux(i,j)
          vx_r = ux(i+1,j)

          if (nint(idom(i,j)).ge.2) then
             if ((vx_l.gt.0 .and. nint(idom(i-1,j)).ge.2) .or. &
                  (vx_r.lt.0 .and. nint(idom(i+1,j)).ge.2)) then
                   temp(i,j) = temp_old(i,j) & 
                               - dt/dx(1) * ( (vx_l+ABS(vx_l))*(temp_old(i,j)-temp_old(i-1,j)) &
                               + (vx_r-ABS(vx_r))*(temp_old(i+1,j)-temp_old(i,j)) )/2.0_amrex_real
             end if
          end if
          
          if(temp(i,j).ne.temp_old(i,j)) then
             call get_enthalpy(temp(i,j),u_new(i,j))
          end if
          
       end do
    end do
    
  end subroutine advection_box


  ! -----------------------------------------------------------------
  ! Subroutine used to update the position of the bottom of the
  ! melt pool
  ! -----------------------------------------------------------------
  subroutine get_melt_pool_bottom()
    
    use amr_data_module, only : phi_new,&
                                idomain
    use heat_transfer_domain_module, only : get_melt_pos

    
    ! Local variables
    integer :: ilev
    type(amrex_geometry) :: geom
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidom
    
    ! Find the melt position
    ilev = amrex_max_level
    
    ! Geometry
    geom = amrex_geom(ilev)
    
    ! Loop through all the boxes in the level
    !$omp parallel private(mfi, bx, pidom)
    call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.false.)
    do while(mfi%next())
       bx = mfi%validbox()
       pidom  => idomain(ilev)%dataptr(mfi)
       call get_melt_pos(bx%lo, bx%hi, &
                         pidom, lbound(pidom), ubound(pidom), &
                         geom)
       
    end do
    call amrex_mfiter_destroy(mfi)
    !$omp end parallel   

    
  end subroutine get_melt_pool_bottom
  
end module heat_transfer_implicit_module
