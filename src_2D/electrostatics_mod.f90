module electrostatics_module
  
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
  public :: advance_electrostatics
  public :: advance_electrostatics_level
  
contains

  
  
  ! -----------------------------------------------------------------
  ! Subroutine used to advance the heat solver on all levels
  ! -----------------------------------------------------------------
  subroutine advance_electrostatics()

    ! Local variables
    integer :: ilev

    ! Advance heat solver
    do ilev = 0, amrex_max_level
       call advance_electrostatics_level(ilev)
    end do
    
  end subroutine advance_electrostatics

  
  ! -----------------------------------------------------------------
  ! Subroutine used to advance the heat solver on one level
  ! -----------------------------------------------------------------
  subroutine advance_electrostatics_level(lev)
    
    ! Input and output variables
    integer, intent(in) :: lev
    
    ! Advance the conductive part of the heat equation
    call compute_potential(lev)
    
  end subroutine advance_electrostatics_level
  
  ! -----------------------------------------------------------------
  ! Subroutine used to advance of one time step the conduction part
  ! of the heat transfer equation
  ! -----------------------------------------------------------------
  subroutine compute_potential(lev)
    
    ! Input and output variables
    integer, intent(in) :: lev
    
    ! Local variables
    type(amrex_multifab) :: rhs ! multifab for the right hand side of the linear system of equations
    type(amrex_multifab) :: alpha ! multifab for the first term on the left hand side of the linear system of equations
    type(amrex_multifab) :: beta(amrex_spacedim) ! multifab for the second term on the left hand side of the linear system of equations
    type(amrex_boxarray) :: ba ! Box array
    type(amrex_distromap):: dm ! Distribution mapping

    ! Initialize temporary multifabs
    call init_tmp_multifab(lev, ba, dm, alpha, beta, rhs)

    ! Predictor step on all levels
    call setup(lev, alpha, beta, rhs)    
    
    ! New temperature via implicit update
    call get_potential(lev, ba, dm, alpha, beta, rhs)
    
    ! Clean memory
    call free_memory(rhs, alpha, beta, dm)
    
  end subroutine compute_potential

  ! -----------------------------------------------------------------
  ! Subroutine used to initialize the multifabs used in the
  ! update of the conductive part of the heat equation
  ! -----------------------------------------------------------------
  subroutine init_tmp_multifab(lev, ba, dm, alpha, beta, rhs)
    
    use amr_data_module, only : phi_new

    ! Input and output variables
    integer, intent(in) :: lev
    type(amrex_multifab), intent(out) :: rhs
    type(amrex_multifab), intent(out) :: alpha
    type(amrex_multifab), intent(out) :: beta(amrex_spacedim)
    type(amrex_boxarray), intent(out) :: ba
    type(amrex_distromap), intent(out) :: dm
    
    ! Local variables
    logical :: nodal(2)
    integer :: ncomp
    integer :: idim
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
       
  end subroutine init_tmp_multifab
  
  ! -----------------------------------------------------------------
  ! Subroutine used to compute the predictor step of the conductive
  ! part of the heat transfer equation
  ! -----------------------------------------------------------------
  subroutine setup(lev, alpha, beta, rhs)
    
    use amr_data_module, only : phi_new
    
    ! Input and output variables
    integer, intent(in) :: lev
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
       call setup_box(lev, mfi, geom, alpha, beta, rhs)
    end do
    call amrex_mfiter_destroy(mfi)
    !$omp end parallel

  end subroutine setup
  
  ! -----------------------------------------------------------------
  ! Subroutine used to perform the prediction step for the implicit
  ! enthalpy update at a given box on a given level 
  ! -----------------------------------------------------------------
  subroutine setup_box(lev, mfi, geom, alpha, beta, rhs)

    use amr_data_module, only : temp, idomain
    
    ! Input and output variables
    integer, intent(in) :: lev
    type(amrex_mfiter), intent(in) :: mfi
    type(amrex_geometry), intent(in) :: geom
    type(amrex_multifab), intent(inout) :: alpha
    type(amrex_multifab), intent(inout) :: beta(amrex_spacedim)
    type(amrex_multifab), intent(inout) :: rhs

    ! Local variables
    integer :: idim
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptemp
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidout
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pac
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pbc
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: prhs
    type(amrex_box) :: bx
    
    ! Box
    bx = mfi%validbox()   
          
    ! Pointers
    ptemp   => temp(lev)%dataptr(mfi)
    pidout  => idomain(lev)%dataptr(mfi)
    pac     => alpha%dataptr(mfi)
    prhs    => rhs%dataptr(mfi)
    
    ! Get alpha matrix for the linear solver (first term on left hand side)
    call get_alpha(bx%lo, bx%hi, pac, lbound(pac), ubound(pac))

    
    ! Get beta matrix for the linear solver (second term on left hand side)
    do idim = 1,amrex_spacedim
       pbc => beta(idim)%dataptr(mfi)
       call get_beta(bx%lo, bx%hi, idim, &
                     pidout, lbound(pidout), ubound(pidout), &
                     ptemp, lbound(ptemp), ubound(ptemp), &
                     pbc, lbound(pbc), ubound(pbc))
    end do
    
    ! Get right hand for the linear solver
    call get_rhs(bx%lo, bx%hi, &
                 geom%get_physical_location(bx%lo), geom%dx, &
                 pidout, lbound(pidout), ubound(pidout), &
                 prhs, lbound(prhs), ubound(prhs))
        
  end subroutine setup_box
  
  
  ! -----------------------------------------------------------------------------------
  ! Subroutine used to solve the linear system of equations for the implict update
  ! of the temperature
  ! ----------------------------------------------------------------------------------- 
  subroutine get_potential(lev, ba, dm, alpha, beta, rhs)

    use amr_data_module, only : psi
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
    type(amrex_multifab), intent(in) :: rhs
    type(amrex_multifab), intent(in) :: alpha
    type(amrex_multifab), intent(in) :: beta(amrex_spacedim)
    type(amrex_boxarray), intent(in) :: ba
    type(amrex_distromap), intent(in) :: dm
    
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
    call abeclap%set_domain_bc([amrex_lo_dirichlet, amrex_lo_dirichlet], &
                               [amrex_lo_dirichlet, amrex_lo_dirichlet])
    
    ! Coarse-fine interface  (for problem with pure homogeneous Neumann BC we can pass an empty multifab)
    if (lev > 0) then
       ! use coarse level data to set up bc at corase/fine boundary
       call abeclap % set_coarse_fine_bc(psi(lev-1), amrex_ref_ratio(lev-1))
    end if
    call abeclap%set_level_bc(0, nullmf)
    
    ! Set A and B scalar
    call abeclap%set_scalars(0.0_amrex_real, 1.0_amrex_real)
    
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
    err = multigrid%solve([psi(lev)], [rhs], num_accuracy, 0.0_amrex_real)
    
    ! Clean memory
    call amrex_abeclaplacian_destroy(abeclap)
    call amrex_multigrid_destroy(multigrid)
    

  end subroutine get_potential

  ! -----------------------------------------------------------------
  ! Subroutine used to advance of one time step the conduction part
  ! of the heat transfer equation
  ! -----------------------------------------------------------------
  subroutine free_memory(rhs, alpha, beta, dm)
    
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
    
  end subroutine free_memory

  
  ! -----------------------------------------------------------------
  ! Subroutine used to fill the alpha matrix for the linear solver
  ! -----------------------------------------------------------------  
  subroutine get_alpha(lo, hi, alpha, a_lo, a_hi)

    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: a_lo(2), a_hi(2)
    real(amrex_real), intent(out) :: alpha(a_lo(1):a_hi(1),a_lo(2):a_hi(2))

    ! Local variables
    integer :: i,j

    ! Fill alpha matrix
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          alpha(i,j) = 0.0_amrex_real
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
  				
    use material_properties_module, only: get_electrical_resistivity, temp_melt

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
                ! temp_face = (temp(i,j) + temp(i-1,j))/2_amrex_real
                 ! call get_conductivity(temp_face, beta(i,j))
                 call get_electrical_resistivity(temp_melt, beta(i,j))
                 beta(i,j) = 1.0/beta(i,j)
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
                ! temp_face = (temp(i,j) + temp(i,j-1))/2_amrex_real
                ! call get_conductivity(temp_face, beta(i,j))
                call get_electrical_resistivity(temp_melt, beta(i,j))
                beta(i,j) = 1.0/beta(i,j)
             end if
             
          end do
       end do

          
    end if

  end subroutine get_beta

  ! -----------------------------------------------------------------
  ! Subroutine used to get the right hand side of the linear solver
  ! -----------------------------------------------------------------  
  subroutine get_rhs(lo, hi, &
                     lo_phys, dx, &
                     idom, id_lo, id_hi, &
                     rhs, r_lo, r_hi)

    use amr_data_module, only : surf_current
    use heat_transfer_domain_module, only : interp_to_max_lev
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: id_lo(2), id_hi(2)
    integer, intent(in) :: r_lo(2), r_hi(2)
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2))
    real(amrex_real), intent(out) :: rhs(r_lo(1):r_hi(1),r_lo(2):r_hi(2))
    real(amrex_real), intent(in) :: lo_phys(2)
    real(amrex_real), intent(in) :: dx(2)
    
    ! Local variables
    integer :: i,j
    integer :: xind
    
    ! Fill rhs of linear problem
    do j = lo(2), hi(2)          
       do i = lo(1), hi(1)
          if(nint(idom(i,j)).ne.0 .and. nint(idom(i,j+1)).eq.0) then
             call interp_to_max_lev(lo, lo_phys, dx, i, xind)
             rhs(i,j) = surf_current(xind)/dx(2) 
          else
             rhs(i,j) = 0.0_amrex_real
          end if
       end do
    end do

  end subroutine get_rhs

    
end module electrostatics_module
