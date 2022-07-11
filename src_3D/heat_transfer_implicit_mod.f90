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
    public :: advance_heat_solver_implicit_level
    
  contains



  ! -----------------------------------------------------------------
  ! Subroutine used to advance the heat solver on all levels
  ! -----------------------------------------------------------------
  subroutine advance_heat_solver_implicit(time, dt)

   ! Input and output variables
   real(amrex_real), intent(in) :: dt
   real(amrex_real), intent(in) :: time

   ! Local variables
   integer :: ilev
   type(amrex_multifab) :: temp_lm1(0:amrex_max_level) ! Temperature multifab at lev-1
   

   ! Advance heat solver
   do ilev = 0, amrex_max_level
      call advance_heat_solver_implicit_level(ilev, time, dt, 1, temp_lm1)
   end do

   ! Clean memory
   do ilev = 0, amrex_max_level
      call amrex_multifab_destroy(temp_lm1(ilev))
   end do
   
 end subroutine advance_heat_solver_implicit  

  
  ! -----------------------------------------------------------------
  ! Subroutine used to advance the heat solver on one level
  ! -----------------------------------------------------------------
  subroutine advance_heat_solver_implicit_level(lev, time, dt, substep, &
                                                temp_tmp_lm1)

    use heat_transfer_domain_module, only : reset_melt_pos, &
                                            update_melt_pos
    use heat_transfer_explicit_module, only : init_tmp_multifab, &
                                              free_memory, &
                                              update_flux_registers
    
    ! Input and output variables
    integer, intent(in) :: lev
    integer, intent(in) :: substep
    real(amrex_real), intent(in) :: dt
    real(amrex_real), intent(in) :: time
    type(amrex_multifab), intent(inout) :: temp_tmp_lm1(0:amrex_max_level)
    
    ! Local variables
    type(amrex_multifab) :: phi_tmp  ! Enthalpy multifab with ghost points
    type(amrex_multifab) :: temp_tmp ! Temperature multifab with ghost points
    type(amrex_multifab) :: idomain_tmp ! Idomain multifab to distinguish material and background
    type(amrex_multifab) :: fluxes(amrex_spacedim)

    ! Initialize temporary multifabs
    call init_tmp_multifab(lev, time, phi_tmp, temp_tmp, idomain_tmp, fluxes)

    ! Reset melt position (to properly treat resolidification)
    if (lev.eq.amrex_max_level) call reset_melt_pos
    
    ! Advance the conductive part of the heat equation
    call advance_conduction(lev, time, dt, phi_tmp, temp_tmp, &
                            temp_tmp_lm1, idomain_tmp)

    ! Advance the advective part of the heat equation
    call advance_advection(lev, dt, substep, temp_tmp, phi_tmp, fluxes)

    ! Update flux registers 
    call update_flux_registers(lev, fluxes)
   
    ! Update melt position
    call update_melt_pos(lev)

    ! Set background to avoind errors at the call of averagedown in simulation module
    call synch_background(lev)

   ! Clean memory
   call free_memory(phi_tmp, temp_tmp, idomain_tmp, fluxes)

  end subroutine advance_heat_solver_implicit_level

  ! -----------------------------------------------------------------
  ! Subroutine used to advance of one time step the conduction part
  ! of the heat transfer equation
  ! -----------------------------------------------------------------
  subroutine advance_conduction(lev, time, dt, phi_tmp, temp_tmp, &
                                temp_tmp_lm1, idomain_tmp)
    
    ! Input and output variables
    integer, intent(in) :: lev
    real(amrex_real), intent(in) :: dt
    real(amrex_real), intent(in) :: time
    type(amrex_multifab), intent(inout) :: temp_tmp_lm1(0:amrex_max_level)
    type(amrex_multifab), intent(inout) :: phi_tmp
    type(amrex_multifab), intent(inout) :: temp_tmp
    type(amrex_multifab), intent(inout) :: idomain_tmp 
    
    ! Local variables
    type(amrex_multifab) :: rhs ! multifab for the right hand side of the linear system of equations
    type(amrex_multifab) :: alpha ! multifab for the first term on the left hand side of the linear system of equations
    type(amrex_multifab) :: beta(amrex_spacedim) ! multifab for the second term on the left hand side of the linear system of equations
    type(amrex_boxarray) :: ba ! Box array
    type(amrex_distromap):: dm ! Distribution mapping

    ! Initialize temporary multifabs
    call conduction_init_tmp_multifab(lev, ba, dm, alpha, beta, rhs)

    ! Predictor step on all levels
    call conduction_predict(lev, time, dt, phi_tmp, temp_tmp, &
                            idomain_tmp, alpha, beta, rhs)    
    
    ! New temperature via implicit update
    call conduction_get_temperature(lev, ba, dm, dt, alpha, beta, &
                                    rhs, temp_tmp_lm1(lev-1), temp_tmp)

    ! Call temperature multifab at level-1
    call update_temperature_lm1(lev, temp_tmp, temp_tmp_lm1(lev))

    ! Corrector step on all levels
    call conduction_correct(lev, phi_tmp, temp_tmp, alpha)

    ! Synchronize idomains
    call conduction_synch_idomain(lev, phi_tmp, temp_tmp)

    ! Clean memory
    call conduction_free_memory(rhs, alpha, beta, dm)
    
  end subroutine advance_conduction


  ! -----------------------------------------------------------------
  ! Subroutine used to advance the advective part of the heat
  ! transfer equation on one level
  ! -----------------------------------------------------------------
  subroutine advance_advection(lev, dt, substep, temp_tmp, phi_tmp, fluxes)
    
    use amr_data_module, only : phi_new
    
    ! Input and output variables
    integer, intent(in) :: lev
    integer, intent(in) :: substep
    real(amrex_real), intent(in) :: dt
    type(amrex_multifab), intent(inout) :: phi_tmp
    type(amrex_multifab), intent(inout) :: temp_tmp
    type(amrex_multifab), intent(inout) :: fluxes(amrex_spacedim)
    
    ! Local variables
    integer :: ncomp
    integer :: idim
    type(amrex_fab) :: flux(amrex_spacedim)
    type(amrex_geometry) :: geom
    type(amrex_mfiter) :: mfi     

    ! Get geometry
    geom = amrex_geom(lev)

    ! Get number of components
    ncomp = phi_new(lev)%ncomp()

    ! Update heat equation
    !$omp parallel private(mfi, flux)    
    do idim = 1, amrex_spacedim
       call flux(idim)%reset_omp_private()
    end do

    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.false.)
    do while(mfi%next())

       call advection_box(lev, dt, substep, mfi, geom, ncomp, &
                          phi_tmp, temp_tmp, flux, fluxes)
    end do
    call amrex_mfiter_destroy(mfi)


    ! Clean memory
    do idim = 1, amrex_spacedim
      call amrex_fab_destroy(flux(idim))
   end do
    !$omp end parallel

       
  end subroutine advance_advection


  ! -----------------------------------------------------------------
  ! Subroutine used update the temperature solution at one level
  ! below (needed for the implicit update of the conduction)
  ! -----------------------------------------------------------------
  subroutine update_temperature_lm1(lev, temp_tmp, temp_tmp_lm1)

   use heat_transfer_domain_module, only : reset_melt_pos
   use amr_data_module, only : phi_new

   ! Input and output variables
   integer, intent(in) :: lev
   type(amrex_multifab), intent(in) :: temp_tmp    
   type(amrex_multifab), intent(out) :: temp_tmp_lm1

   ! Local variables
   integer, parameter :: nghost = 1 ! number of ghost points in each spatial direction
   integer :: ncomp
   type(amrex_boxarray) :: ba ! Box array
   type(amrex_distromap):: dm ! Distribution mapping
   type(amrex_geometry) :: geom

   ncomp = phi_new(lev)%ncomp()
   ba = phi_new(lev)%ba
   dm = phi_new(lev)%dm
   geom = amrex_geom(lev)
      
   ! Destroy previous definitions of the lm1 multifab
   if (lev > 0) call amrex_multifab_destroy(temp_tmp_lm1)

   ! Rebuild the multifab to store the temperature data of the current level
   call amrex_multifab_build(temp_tmp_lm1, ba, dm, ncomp, nghost)

   ! Copy temperature of the current level to the lm1 multifab
   call temp_tmp_lm1%copy(temp_tmp, 1, 1, ncomp, nghost) 
   call temp_tmp_lm1%fill_boundary(geom)
   
 end subroutine update_temperature_lm1


   ! -----------------------------------------------------------------
   ! Subroutine used to initialize the multifabs used in the
   ! update of the conductive part of the heat equation
   ! -----------------------------------------------------------------
   subroutine conduction_init_tmp_multifab(lev, ba, dm, alpha, beta, rhs)
      
      use amr_data_module, only : phi_new

      ! Input and output variables
      integer, intent(in) :: lev
      type(amrex_multifab), intent(out) :: rhs
      type(amrex_multifab), intent(out) :: alpha
      type(amrex_multifab), intent(out) :: beta(amrex_spacedim)
      type(amrex_boxarray), intent(out) :: ba
      type(amrex_distromap), intent(out) :: dm
      
      ! Local variables
      logical :: nodal(3)
      integer :: idim
      integer :: ncomp
      type(amrex_geometry) :: geom
      
      ! Get boxarray and distribution mapping
      ba = phi_new(lev)%ba
      dm = phi_new(lev)%dm

      ! Number of components
      ncomp = phi_new(lev)%ncomp()

      ! Geometry
      geom = amrex_geom(lev)

      ! Multifabs for the linear solver
      call amrex_multifab_build(rhs, ba, dm, ncomp, 0)
      call amrex_multifab_build(alpha, ba, dm, ncomp, 0)
      do idim = 1, amrex_spacedim
         nodal = .false.
         nodal(idim) = .true.
         call amrex_multifab_build(beta(idim), ba, dm, ncomp, 0, nodal)
      end do
      
   end subroutine conduction_init_tmp_multifab

  
    ! -----------------------------------------------------------------
    ! Subroutine used to compute the predictor step of the conductive
    ! part of the heat transfer equation
    ! -----------------------------------------------------------------
    subroutine conduction_predict(lev, time, dt, phi_tmp, temp_tmp, &
                                    idomain_tmp, alpha, beta, rhs)
        
        use amr_data_module, only : phi_new

        ! Input and output variables
        integer, intent(in) :: lev
        real(amrex_real), intent(in) :: dt
        real(amrex_real), intent(in) :: time
        type(amrex_multifab), intent(inout) :: phi_tmp
        type(amrex_multifab), intent(inout) :: temp_tmp
        type(amrex_multifab), intent(inout) :: idomain_tmp 
        type(amrex_multifab), intent(inout) :: rhs
        type(amrex_multifab), intent(inout) :: alpha
        type(amrex_multifab), intent(inout) :: beta(amrex_spacedim)

        ! Local variables
        type(amrex_geometry) :: geom
        type(amrex_mfiter) :: mfi

        ! Geometry
        geom = amrex_geom(lev)
        
        !$omp parallel private(mfi)
        call amrex_mfiter_build(mfi, phi_new(lev), tiling=.false.)
        do while(mfi%next())
            call conduction_predict_box(lev, time, dt, mfi, &
                                        geom, phi_tmp, temp_tmp, &
                                        idomain_tmp, alpha, beta, &
                                        rhs)

        end do
        call amrex_mfiter_destroy(mfi)
        !$omp end parallel

    
    end subroutine conduction_predict

    ! -----------------------------------------------------------------
    ! Subroutine used to perform the prediction step for the implicit
    ! enthalpy update at a given box on a given level 
    ! -----------------------------------------------------------------
    subroutine conduction_predict_box(lev, time, dt, mfi, &
                                      geom, phi_tmp, temp_tmp, &
                                      idomain_tmp, alpha, beta, &
                                      rhs)
  
      use read_input_module, only : heat_temp_surf
      use amr_data_module, only : phi_new, & 
                                  temp, & 
                                  idomain, &
                                  Qplasma, &
                                  Qpipe, &
                                  Qtherm, &
                                  Qvap, &
                                  Qrad
      use heat_transfer_domain_module, only : get_idomain, revaluate_heat_domain
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
      real(amrex_real) :: Qpipe_box, Qplasma_box, Qtherm_box, Qrad_box, Qvap_box
      
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
      if (heat_temp_surf.gt.0) then   
         STOP 'The implicit solver cannot be used to solve problems with a fixed temperature on the free surface'
      else
         call get_rhs(lev, bx%lo, bx%hi, time, dt, &
                      geom%get_physical_location(bx%lo), geom%dx, &
                      pidout, lbound(pidout), ubound(pidout), &
                      ptempin, lbound(ptempin), ubound(ptempin), &
                      pac, lbound(pac), ubound(pac), &                      
                      prhs, lbound(prhs), ubound(prhs), &
                      Qpipe_box, Qtherm_box, Qvap_box, Qrad_box, Qplasma_box)
      end if
      !$omp critical
      Qplasma = Qplasma + Qplasma_box
      Qpipe = Qpipe + Qpipe_box
      Qtherm = Qtherm + Qtherm_box
      Qvap = Qvap + Qvap_box
      Qrad = Qrad + Qrad_box
      !$omp end critical
      
    end subroutine conduction_predict_box

        ! -----------------------------------------------------------------------------------
    ! Subroutine used to solve the linear system of equations for the implict update
    ! of the temperature
    ! ----------------------------------------------------------------------------------- 
    subroutine conduction_get_temperature(lev, ba, dm, dt, alpha, beta, rhs, sol_lm1, sol)
  
    use read_input_module, only : num_accuracy, &
                                  num_agglomeration, &
                                  num_consolidation, &
                                  num_max_coarsening_level, &
                                  num_linop_maxorder, &
                                  num_bottom_solver, &
                                  num_bottom_verbose, &
                                  num_max_fmg_iter, &
                                  num_max_iter, &
                                  num_verbose
      use amrex_linear_solver_module
  
          
      ! Input and output variables
      integer, intent(in) :: lev
      type(amrex_multifab), intent(inout) :: sol
      type(amrex_multifab), intent(in) :: rhs
      type(amrex_multifab), intent(in) :: alpha
      type(amrex_multifab), intent(in) :: beta(amrex_spacedim)
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_multifab), intent(in) :: sol_lm1    
      real(amrex_real), intent(in) :: dt
      
      ! Local variables
      real(amrex_real) :: err ! Error from the solution of the linear system of equations
      type(amrex_multifab) :: nullmf ! Null multifab used to define the boundary condition at the coarse-fine interface
      type(amrex_abeclaplacian) :: abeclap ! Object for the solution of the linear systems of equations
      type(amrex_multigrid) :: multigrid ! Object for the solution of the linear systems of equations
  
  
      ! Solve linear system of equations
         
      call amrex_abeclaplacian_build(abeclap, [amrex_geom(lev)], [ba], [dm], &
                                       metric_term=.false., agglomeration=num_agglomeration, &
                                       consolidation=num_consolidation, &
                                       max_coarsening_level=num_max_coarsening_level)

      call abeclap%set_maxorder(num_linop_maxorder)
      
      ! This is set up to have homogeneous Neumann BC
      call abeclap%set_domain_bc([amrex_lo_neumann, amrex_lo_neumann, amrex_lo_neumann], &
                                 [amrex_lo_neumann, amrex_lo_neumann, amrex_lo_neumann])
      
      ! Coarse-fine interface  (for problem with pure homogeneous Neumann BC we can pass an empty multifab)
      if (lev > 0) then
         ! use coarse level data to set up bc at corase/fine boundary
         call abeclap % set_coarse_fine_bc(sol_lm1, amrex_ref_ratio(lev-1))
      end if
      call abeclap%set_level_bc(0, nullmf)

      ! Set A and B scalar
      call abeclap%set_scalars(1.0_amrex_real, dt)
      
      ! Set alpha and beta matrices
      call abeclap%set_acoeffs(0, alpha)
      call abeclap%set_bcoeffs(0, beta)
      
      ! Build multigrid solver
      call amrex_multigrid_build(multigrid, abeclap)
      call multigrid%set_verbose(num_verbose)
      call multigrid%set_bottom_verbose(num_bottom_verbose)
      call multigrid%set_max_iter(num_max_iter)
      call multigrid%set_max_fmg_iter(num_max_fmg_iter)
      call multigrid%set_bottom_solver(num_bottom_solver)
      
      ! Solve the linear system
      err = multigrid%solve([sol], [rhs], num_accuracy, 0.0_amrex_real)
      
      ! Clean memory
      call amrex_abeclaplacian_destroy(abeclap)
      call amrex_multigrid_destroy(multigrid) 
  
            
    end subroutine conduction_get_temperature

    ! -----------------------------------------------------------------
    ! Subroutine used to compute the correction step of the conductive
    ! part of the heat transfer equation
    ! -----------------------------------------------------------------
    subroutine conduction_correct(lev, phi_tmp, temp_tmp, alpha)
    
        use amr_data_module, only : phi_new
        
        ! Input and output variables
        integer, intent(in) :: lev
        type(amrex_multifab), intent(inout) :: phi_tmp
        type(amrex_multifab), intent(inout) :: temp_tmp
        type(amrex_multifab), intent(inout) :: alpha
        
        ! Local variables
        type(amrex_geometry) :: geom
        type(amrex_mfiter) :: mfi
    
    
         ! Geometry
         geom = amrex_geom(lev)
         
         !$omp parallel private(mfi)
         call amrex_mfiter_build(mfi, phi_new(lev), tiling=.false.)
         do while(mfi%next())
            call conduction_correct_box(lev, mfi, phi_tmp, &
                                       temp_tmp, alpha)
            
         end do
         call amrex_mfiter_destroy(mfi)
         !$omp end parallel
         
         ! call temp_tmp%fill_boundary(geom)
           
        
    end subroutine conduction_correct  
  
    ! -----------------------------------------------------------------
    ! Subroutine used to perform the correction step for the implicit
    ! enthalpy update at a given box on a given level 
    ! -----------------------------------------------------------------
    subroutine conduction_correct_box(lev, mfi, phi_tmp, &
                                                temp_tmp, alpha)
  
      use amr_data_module, only : phi_new, temp, idomain
      use read_input_module, only : heat_temp_surf
      
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
      if (heat_temp_surf.gt.0) then
        STOP 'The implicit solver cannot be used to solve problems with a fixed temperature on the free surface'
      else
         call get_enthalpy_conduction(bx%lo, bx%hi, &
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
    subroutine conduction_synch_idomain(lev, phi_tmp, temp_tmp)
        
      use amr_data_module, only : phi_new, &
                                 temp, &
                                 idomain
      use heat_transfer_domain_module, only : get_idomain
         
      ! Input and output variables
      integer, intent(in) :: lev
      type(amrex_multifab), intent(inout) :: phi_tmp
      type(amrex_multifab), intent(inout) :: temp_tmp
        
      ! Local variables
      integer :: ncomp
      type(amrex_geometry) :: geom
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidom
      real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptemp_tmp
      

      ! Geometry
      geom = amrex_geom(lev)
   
      ! Number of components
      ncomp = phi_new(lev)%ncomp()

      ! Synchronize enthalpy multifab with ghost points
      call phi_tmp%copy(phi_new(lev), 1, 1, ncomp, 0) ! The last 0 is the number of ghost points
      call phi_tmp%fill_boundary(geom)

      ! Synchronize temperature multifab with ghost points
      call temp_tmp%copy(temp(lev), 1, 1, ncomp, 0)
      call temp_tmp%fill_boundary(geom)

      !$omp parallel private(mfi, bx, pidom, ptemp_tmp)
      call amrex_mfiter_build(mfi, phi_new(lev), tiling=.false.)
      do while(mfi%next())
         bx = mfi%validbox()
         pidom  => idomain(lev)%dataptr(mfi)
         ptemp_tmp   => temp_tmp%dataptr(mfi)
         call get_idomain(geom%get_physical_location(bx%lo), geom%dx, &
                           bx%lo, bx%hi, &
                           pidom, lbound(pidom), ubound(pidom), &
                           ptemp_tmp, lbound(ptemp_tmp), ubound(ptemp_tmp))
         
      end do
      call amrex_mfiter_destroy(mfi)   
      !$omp end parallel
      
        
    end subroutine conduction_synch_idomain    
    
  
  ! -----------------------------------------------------------------
  ! Subroutine used to advance of one time step the conduction part
  ! of the heat transfer equation
  ! -----------------------------------------------------------------
    subroutine conduction_free_memory(rhs, alpha, beta, dm)
    
      ! Input and output variables
      type(amrex_multifab) :: rhs
      type(amrex_multifab) :: alpha
      type(amrex_multifab) :: beta(amrex_spacedim) 
      type(amrex_distromap):: dm
      
      ! Local variables
      integer :: idim
  
      call amrex_multifab_destroy(rhs)
      call amrex_multifab_destroy(alpha)
      do idim = 1, amrex_spacedim
         call amrex_multifab_destroy(beta(idim))
      end do
      call amrex_distromap_destroy(dm)
      
    end subroutine conduction_free_memory

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
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: uo_lo(3), uo_hi(3)
      integer, intent(in) :: t_lo(3), t_hi(3)
      integer, intent(in) :: a_lo(3), a_hi(3)
      real(amrex_real), intent(in) :: u_old(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3))  
      real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
      real(amrex_real), intent(out) :: alpha(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))
  
      ! Local variables
      integer :: i,j,k
      real(amrex_real) :: cp
      real(amrex_real) :: cps
      real(amrex_real) :: cpl
      real(amrex_real) :: rho
      real(amrex_real) :: rhos
      real(amrex_real) :: rhol
      real(amrex_real) :: lf
  
      ! Fill alpha matrix
      do i = lo(1), hi(1)
         do j = lo(2), hi(2)
             do k = lo(3), hi(3)
                if (temp(i,j,k).ne.temp_melt) then
                   call get_mass_density(temp(i,j,k), rho)
                   call get_heat_capacity(temp(i,j,k), cp)
                   alpha(i,j,k) = rho*cp
                else
                   ! Compute liquid fraction
                   lf = (u_old(i,j,k) - enth_at_melt)/latent_heat
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
                   alpha(i,j,k) = rho*cp
                end if
          end do
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
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: dim
      integer, intent(in) :: t_lo(3), t_hi(3)
      integer, intent(in) :: id_lo(3), id_hi(3)
      integer, intent(in) :: b_lo(3), b_hi(3)    
      real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2),id_lo(3):id_hi(3))
      real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
      real(amrex_real), intent(out) :: beta(b_lo(1):b_hi(1),b_lo(2):b_hi(2),b_lo(3):b_hi(3))
  
      ! Local variables
      integer :: i,j,k
      real(amrex_real) :: temp_face
  
      if (dim == 1) then ! x-direction
         do i = lo(1), hi(1)+1
            do j = lo(2), hi(2)
               do k = lo(3), hi(3)
               
                  if (nint(idom(i-1,j,k)).eq.0 .or. nint(idom(i,j,k)).eq.0 .or. &
                   nint(idom(i-1,j,k)).eq.-1 .or. nint(idom(i,j,k)).eq.-1) then
                        beta(i,j,k) = 0_amrex_real ! Suppress flux at the free surface
                  else
                     temp_face = (temp(i,j,k) + temp(i-1,j,k))/2_amrex_real
                     call get_conductivity(temp_face, beta(i,j,k))
                  end if
 
               end do 
            end do
         end do
  
      else if (dim == 2) then ! y-direction
         do i = lo(1), hi(1)
            do j = lo(2), hi(2)+1
               do k = lo(3), hi(3)
               
                  if (nint(idom(i,j-1,k)).eq.0 .or. nint(idom(i,j,k)).eq.0 .or. &
                  nint(idom(i,j-1,k)).eq.-1 .or. nint(idom(i,j,k)).eq.-1) then
                        beta(i,j,k) = 0_amrex_real ! Suppress flux at the free surface and surface of cooling pipe
                  else
                     temp_face = (temp(i,j,k) + temp(i,j-1,k))/2_amrex_real
                     call get_conductivity(temp_face, beta(i,j,k))
                  end if
 
               end do
            end do
         end do
 
 
       else if (dim == 3) then ! z-direction
          do i = lo(1), hi(1)
             do j = lo(2), hi(2)
                do k = lo(3), hi(3)+1
                
                   if (nint(idom(i,j,k-1)).eq.0 .or. nint(idom(i,j,k)).eq.0 .or. &
                    nint(idom(i,j,k-1)).eq.-1 .or. nint(idom(i,j,k)).eq.-1) then
                         beta(i,j,k) = 0_amrex_real ! Suppress flux at the free surface and surface of cooling pipe
                   else
                      temp_face = (temp(i,j,k) + temp(i,j,k-1))/2_amrex_real
                      call get_conductivity(temp_face, beta(i,j,k))
                   end if
  
                end do
             end do
          end do
  
            
      end if
  
    end subroutine get_beta
  
    ! -----------------------------------------------------------------
    ! Subroutine used to get the right hand side of the linear solver
    ! -----------------------------------------------------------------  
    subroutine get_rhs(lev, lo, hi, time, dt, &
                       lo_phys, dx, &
                       idom, id_lo, id_hi, &
                       temp, t_lo, t_hi, &
                       alpha, a_lo, a_hi, &
                       rhs, r_lo, r_hi, &
                       Qpipe_box, Qtherm_box, Qvap_box, Qrad_box, Qplasma_box)
  
      use heat_flux_module, only: get_boundary_heat_flux
      
      ! Input and output variables
      integer, intent(in) :: lev
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: id_lo(3), id_hi(3)
      integer, intent(in) :: t_lo(3), t_hi(3)
      integer, intent(in) :: r_lo(3), r_hi(3)
      integer, intent(in) :: a_lo(3), a_hi(3)
      real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2),id_lo(3):id_hi(3))
      real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
      real(amrex_real), intent(in) :: alpha(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))   
      real(amrex_real), intent(out) :: rhs(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
      real(amrex_real), intent(in) :: lo_phys(3)
      real(amrex_real), intent(in) :: dx(3)
      real(amrex_real), intent(in) :: time
      real(amrex_real), intent(in) :: dt
      real(amrex_real), intent(out) :: Qplasma_box
      real(amrex_real), intent(out) :: Qpipe_box
      real(amrex_real), intent(out) :: Qtherm_box
      real(amrex_real), intent(out) :: Qrad_box
      real(amrex_real), intent(out) :: Qvap_box
      
      ! Local variables
      integer :: i,j,k
  
      ! Get the boundary heat flux
      call get_boundary_heat_flux(time, lo_phys, &
                                  dx, lo, hi, &
                                  idom, id_lo, id_hi, &
                                  temp, t_lo, t_hi, lev, dt, &
                                  Qpipe_box, Qtherm_box, Qvap_box, Qrad_box, &
                                  Qplasma_box, rhs)
      
      ! Fill rhs of linear problem
      do i = lo(1), hi(1)
         do j = lo(2), hi(2)
            do k = lo(3), hi(3)          
                rhs(i,j,k) = rhs(i,j,k)*dt + temp(i,j,k)*alpha(i,j,k)
            end do 
         end do
      end do
  
    end subroutine get_rhs  
    
    ! -----------------------------------------------------------------
    ! Subroutine used to get the corrected enthalpy
    ! -----------------------------------------------------------------  
    subroutine get_enthalpy_conduction(lo, hi, &
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
       integer, intent(in) :: lo(3), hi(3)  
       integer, intent(in) :: uo_lo(3), uo_hi(3)
       integer, intent(in) :: un_lo(3), un_hi(3)
       integer, intent(in) :: to_lo(3), to_hi(3)
       integer, intent(in) :: tn_lo(3), tn_hi(3)
       integer, intent(in) :: a_lo(3), a_hi(3)
       real(amrex_real), intent(in) :: u_old(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3))
       real(amrex_real), intent(out) :: u_new(un_lo(1):un_hi(1),un_lo(2):un_hi(2),un_lo(3):un_hi(3))
       real(amrex_real), intent(inout) :: temp_old(to_lo(1):to_hi(1),to_lo(2):to_hi(2),to_lo(3):to_hi(3))
       real(amrex_real), intent(inout) :: temp_new(tn_lo(1):tn_hi(1),tn_lo(2):tn_hi(2),tn_lo(3):tn_hi(3))
       real(amrex_real), intent(in) :: alpha(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))
      
       ! Local variables
       integer :: i,j,k
       real(amrex_real) :: ubase
 
       ! General case, when no phase transition
       ! (will be overwritten for phase transitioning control volumes.)
       call get_temp(un_lo, un_hi, &
                     u_new, un_lo, un_hi, &
                     temp_new, tn_lo, tn_hi, .false.)
                     
       ! Get new enthalpy for points undergoing transition 
       ! (might overwrite enthalpy found previously.)
       do i = lo(1), hi(1)
          do j = lo(2), hi(2)
             do k = lo(3), hi(3)
               
                ! Avoid bugs caused by numerical accuracy at ghost points
                if (temp_new(i,j,k).lt.10.0) then
                  temp_new(i,j,k) = 0
                  u_new(i,j,k) = 0
                end if

                if (temp_old(i,j,k).le.temp_melt .and. temp_new(i,j,k).gt.temp_melt) then
                   ubase = max(u_old(i,j,k), enth_at_melt)
                   u_new(i,j,k) = ubase + alpha(i,j,k)*(temp_new(i,j,k) - temp_melt)
                elseif (temp_old(i,j,k).ge.temp_melt .and. temp_new(i,j,k).lt.temp_melt) then
                   ubase = min(u_old(i,j,k), enth_at_melt+latent_heat)
                   u_new(i,j,k) = ubase + alpha(i,j,k)*(temp_new(i,j,k) - temp_melt)
                end if
             end do
          end do
       end do                  
                                  
       ! Update temperature acordingly (without ghost points)
       call get_temp(lo, hi, &
                     u_new, un_lo, un_hi, &
                     temp_old, to_lo, to_hi, .true.)
 
    end subroutine get_enthalpy_conduction
 
  ! -----------------------------------------------------------------  
  ! Subroutine used for an explicit update of the advective part
  ! of the heat equation on a box of a given level
  ! -----------------------------------------------------------------  
  subroutine advection_box(lev, dt, substep, mfi, geom, ncomp, &
                           phi_tmp, temp_tmp, flux, fluxes)

    use amr_data_module, only : phi_new, &
                                temp, &
                                idomain
    use read_input_module, only : heat_reflux
    use material_properties_module, only : get_temp
    
    ! Input and output variables
    integer, intent(in) :: lev
    integer, intent(in) :: substep
    integer, intent(in) :: ncomp
    real(amrex_real), intent(in) :: dt
    type(amrex_mfiter), intent(in) :: mfi
    type(amrex_geometry), intent(in) :: geom
    type(amrex_multifab), intent(inout) :: phi_tmp
    type(amrex_multifab), intent(inout) :: temp_tmp
    type(amrex_multifab), intent(inout) :: fluxes(amrex_spacedim)
    type(amrex_fab), intent(inout) :: flux(amrex_spacedim)
     
    ! Local variables
    integer :: idim
    type(amrex_box) :: bx
    type(amrex_box) :: tbx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidom
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pin
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptempin
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptemp
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pout
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pfx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pfy
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pfz
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pf
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pfab


        
    ! Box
    bx = mfi%validbox()

    ! Pointer
    pin   => phi_tmp%dataptr(mfi)
    pout    => phi_new(lev)%dataptr(mfi)
    ptempin   => temp_tmp%dataptr(mfi)
    ptemp => temp(lev)%dataptr(mfi)
    pidom  => idomain(lev)%dataptr(mfi)
    do idim = 1, amrex_spacedim
       tbx = bx
       call tbx%nodalize(idim)
       call flux(idim)%resize(tbx,ncomp)
       call tbx%grow(substep)
    end do
    pfx => flux(1)%dataptr()
    pfy => flux(2)%dataptr()
    pfz => flux(3)%dataptr()


    ! Update the enthalpy
    call get_enthalpy_advection(bx%lo, bx%hi, &
                                pin, lbound(pin), ubound(pin), &
                                pout, lbound(pout), ubound(pout), &
                                ptempin, lbound(ptempin), ubound(ptempin), &
                                pfx, lbound(pfx), ubound(pfx), &
                                pfy, lbound(pfy), ubound(pfy), &
                                pfz, lbound(pfz), ubound(pfz), &
                                pidom, lbound(pidom), ubound(pidom), &
                                geom, dt)
    
    ! Update the temperature
    call get_temp(bx%lo, bx%hi, &
                  pout, lbound(pout), ubound(pout), &
                  ptemp, lbound(ptemp), ubound(ptemp), .true.) 


    ! Update pointers for flux registers
    if (heat_reflux) then
       
       do idim = 1, amrex_spacedim
          pf => fluxes(idim)%dataptr(mfi)
          pfab => flux(idim)%dataptr()
          tbx = mfi%nodaltilebox(idim)
          pf(tbx%lo(1):tbx%hi(1), tbx%lo(2):tbx%hi(2), tbx%lo(3):tbx%hi(3), :)  = & 
               pfab(tbx%lo(1):tbx%hi(1), tbx%lo(2):tbx%hi(2), tbx%lo(3):tbx%hi(3), :) 
              
       end do
       
    end if
    
  end subroutine advection_box


  ! -----------------------------------------------------------------  
  ! Subroutine used for an explicit update of the advective part
  ! of the heat equation on a box of a given level
  ! -----------------------------------------------------------------  
  subroutine get_enthalpy_advection(lo, hi, &
                                    u_old,  uo_lo, uo_hi, &
                                    u_new, un_lo, un_hi, &
                                    temp, t_lo, t_hi, &
                                    flxx, fx_lo, fx_hi, &
                                    flxy, fy_lo, fy_hi, &
                                    flxz, fz_lo, fz_hi, &
                                    idom, id_lo, id_hi, &
                                    geom, dt)
    
    use read_input_module, only : num_subcycling
    use heat_transfer_explicit_module, only : get_volumetric_heat_source, &
                                              get_face_flux
 
    ! Input and output variables
    integer, intent(in) :: lo(3), hi(3) 
    integer, intent(in) :: uo_lo(3), uo_hi(3)
    integer, intent(in) :: un_lo(3), un_hi(3)
    integer, intent(in) :: t_lo(3), t_hi(3)  
    integer, intent(in) :: fx_lo(3), fx_hi(3)
    integer, intent(in) :: fy_lo(3), fy_hi(3)
    integer, intent(in) :: fz_lo(3), fz_hi(3)
    integer, intent(in) :: id_lo(3), id_hi(3)
    real(amrex_real), intent(in) :: dt
    real(amrex_real), intent(in) :: u_old(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3))  
    real(amrex_real), intent(inout) :: u_new(un_lo(1):un_hi(1),un_lo(2):un_hi(2),un_lo(3):un_hi(3))
    real(amrex_real), intent(inout) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(amrex_real), intent(out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3))
    real(amrex_real), intent(out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3))
    real(amrex_real), intent(out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3))
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2),id_lo(3):id_hi(3))
    type(amrex_geometry), intent(in) :: geom 
    
    ! Local variables
    integer :: i,j,k
    real(amrex_real) :: dx(3)
    real(amrex_real) :: lo_phys(3) 
    real(amrex_real) :: qvol(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
     
    ! Get grid size
    dx = geom%dx(1:3) ! grid width at level 

    ! Get physical location of the lowest corner of the tile box
    lo_phys = geom%get_physical_location(lo)
    
    ! Get enthalpy flux (only advective component)
    call get_face_flux(dx, lo_phys, lo, hi, &
                       u_old, uo_lo, uo_hi, &
                       flxx, fx_lo, fx_hi, &
                       flxy, fy_lo, fy_hi, &
                       flxz, fz_lo, fz_hi, &
                       temp, t_lo, t_hi, &
                       idom, id_lo, id_hi, &
                       .false., .true.)

    ! Volumetric sources (we only need the terms with the velocity divergence)
    call get_volumetric_heat_source(dx, lo_phys, lo, hi, &
                                    u_old, uo_lo, uo_hi, &
                                    idom, id_lo, id_hi, &
                                    qvol)
    
    ! Compute enthalpy at the new timestep
    do k = lo(3),hi(3)
      do  j = lo(2),hi(2)
         do i = lo(1),hi(1)
            u_new(i,j,k) = u_old(i,j,k) &
                           - dt/dx(1) * (flxx(i+1,j,k) - flxx(i,j,k)) &  
                           - dt/dx(2) * (flxy(i,j+1,k) - flxy(i,j,k)) & 
                           - dt/dx(3) * (flxz(i,j,k+1) - flxz(i,j,k)) & 
                           + dt*qvol(i,j,k)
         end do
      end do
   end do
  	
    ! Scale the fluxes for the flux registers
    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
         do i = lo(1), hi(1) + 1
            flxx(i,j,k) = flxx(i,j,k) * dx(2)*dx(3)
            if (num_subcycling) flxx(i,j,k) = flxx(i,j,k) * dt 
         end do
      end do
    end do

    do k = lo(3), hi(3)
      do j = lo(2), hi(2) + 1
         do i = lo(1), hi(1)
            flxy(i,j,k) = flxy(i,j,k) * dx(1)*dx(3)
            if (num_subcycling) flxy(i,j,k) = flxy(i,j,k) * dt 
         end do
      end do
   end do

   do k = lo(3), hi(3) + 1
      do j = lo(2), hi(2) 
         do i = lo(1), hi(1)
            flxz(i,j,k) = flxz(i,j,k) * dx(2)*dx(3)
            if (num_subcycling) flxz(i,j,k) = flxz(i,j,k) * dt 
         end do
      end do
   end do

  end subroutine get_enthalpy_advection


  ! -----------------------------------------------------------------
  ! Subroutine used to set the background temperature to zero and the
  ! background enthalpy to the one of the closest free-surface cell
  ! -----------------------------------------------------------------
  subroutine synch_background(lev)
    
   use amr_data_module, only : phi_new, &
                               idomain, &
                               temp    
   use heat_transfer_domain_module, only : set_backgroung_domain
       
   ! Input and output variables
   integer, intent(in) :: lev
   
   ! Local variables
   type(amrex_geometry) :: geom
   type(amrex_mfiter) :: mfi
   type(amrex_box) :: bx
   real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidom
   real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptemp
   real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: penth
   
   ! Geometry
   geom = amrex_geom(lev)
   
   !$omp parallel private(mfi, bx, pidom, ptemp, penth)
   call amrex_mfiter_build(mfi, phi_new(lev), tiling=.false.)
   do while(mfi%next())
      bx = mfi%validbox()
      pidom  => idomain(lev)%dataptr(mfi)
      ptemp   => temp(lev)%dataptr(mfi)
      penth   => phi_new(lev)%dataptr(mfi)

      call set_backgroung_domain(geom%get_physical_location(bx%lo), geom%dx, &
                                 bx%lo, bx%hi, &
                                 pidom, lbound(pidom), ubound(pidom), &
                                 penth, lbound(penth), ubound(penth), &
                                 ptemp, lbound(ptemp), ubound(ptemp))
      
   end do
   call amrex_mfiter_destroy(mfi)   
   !$omp end parallel
      
  end subroutine synch_background

end module heat_transfer_implicit_module
  
 
