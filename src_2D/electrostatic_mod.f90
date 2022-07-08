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
    subroutine advance_electrostatics(time)
  
      ! Input and output variables
      real(amrex_real), intent(in) :: time
      ! Local variables
      integer :: ilev
  
      ! Advance heat solver
      do ilev = 0, amrex_max_level
         call advance_electrostatics_level(ilev, time)
      end do
      
    end subroutine advance_electrostatics
  
    
    ! -----------------------------------------------------------------
    ! Subroutine used to advance the heat solver on one level
    ! -----------------------------------------------------------------
    subroutine advance_electrostatics_level(lev, time)
      
      ! Input and output variables
      integer, intent(in) :: lev
      real(amrex_real), intent(in) :: time
      
      ! Advance the conductive part of the heat equation
      call compute_potential(lev, time)
      
    end subroutine advance_electrostatics_level
    
    ! -----------------------------------------------------------------
    ! Subroutine used to advance of one time step the conduction part
    ! of the heat transfer equation
    ! -----------------------------------------------------------------
    subroutine compute_potential(lev, time)
      
      ! Input and output variables
      integer, intent(in) :: lev
      real(amrex_real), intent(in) :: time
      
      ! Local variables
      type(amrex_multifab) :: rhs ! multifab for the right hand side of the linear system of equations
      type(amrex_multifab) :: alpha ! multifab for the first term on the left hand side of the linear system of equations
      type(amrex_multifab) :: beta(amrex_spacedim) ! multifab for the second term on the left hand side of the linear system of equations
      type(amrex_multifab) :: temp_tmp
      type(amrex_multifab) :: phi_tmp
      type(amrex_boxarray) :: ba ! Box array
      type(amrex_distromap):: dm ! Distribution mapping
  
      ! Initialize temporary multifabs
      call init_tmp_multifab(lev, time, ba, dm, phi_tmp, temp_tmp, alpha, beta, rhs)
  
      ! Predictor step on all levels
      call setup(lev, phi_tmp, temp_tmp, alpha, beta, rhs)    
      
      ! New temperature via implicit update
      call get_potential(lev, ba, dm, alpha, beta, rhs)

      ! Calculate the averaged current inside the melt pool - only on the highest level
      call get_current (temp_tmp, lev)
      
      ! Clean memory
      call free_memory(phi_tmp, temp_tmp, rhs, alpha, beta, dm)
      
    end subroutine compute_potential
  
    ! -----------------------------------------------------------------
    ! Subroutine used to initialize the multifabs used in the
    ! update of the conductive part of the heat equation
    ! -----------------------------------------------------------------
    subroutine init_tmp_multifab(lev, time, ba, dm, phi_tmp, temp_tmp, alpha, beta, rhs)
      
      use amr_data_module, only : phi_new
      use regrid_module, only : fillpatch
  
      ! Input and output variables
      integer, intent(in) :: lev
      real(amrex_real), intent(in) :: time
      type(amrex_boxarray), intent(out) :: ba
      type(amrex_distromap), intent(out) :: dm
      type(amrex_multifab), intent(out) :: phi_tmp
      type(amrex_multifab), intent(out) :: temp_tmp
      type(amrex_multifab), intent(out) :: alpha
      type(amrex_multifab), intent(out) :: beta(amrex_spacedim)
      type(amrex_multifab), intent(out) :: rhs
      
      ! Local variables
      logical :: nodal(2)
      integer :: ncomp
      integer :: nghost
      integer :: idim
      type(amrex_geometry) :: geom
         
      ! Get boxarray and distribution mapping
      ba = phi_new(lev)%ba
      dm = phi_new(lev)%dm
      
      ! Number of components
      ncomp = phi_new(lev)%ncomp()

      ! Number of ghost points
      nghost = 1
      
      ! Geometry
      geom = amrex_geom(lev)

      ! Build and fill in the enthalpy multifab with ghost points which will
      ! be used to create the temperature multifab with ghost points
      call amrex_multifab_build(phi_tmp, ba, dm, ncomp, nghost) 
      call phi_tmp%copy(phi_new(lev), 1, 1, ncomp, 0) ! The last 0 is the number of ghost points
      call phi_tmp%fill_boundary(geom)
      call fillpatch(lev, time, phi_tmp)

      ! Build the temperature multifab with ghost points
      call amrex_multifab_build(temp_tmp, ba, dm, ncomp, nghost) 
      
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
    subroutine setup(lev, phi_tmp, temp_tmp, alpha, beta, rhs)
      
      use amr_data_module, only : phi_new
      
      ! Input and output variables
      integer, intent(in) :: lev
      type(amrex_multifab), intent(inout) :: phi_tmp
      type(amrex_multifab), intent(inout) :: temp_tmp
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
         call setup_box(lev, mfi, geom, phi_tmp, temp_tmp, alpha, beta, rhs)
      end do
      call amrex_mfiter_destroy(mfi)
      !$omp end parallel
  
    end subroutine setup
    
    ! -----------------------------------------------------------------
    ! Subroutine used to perform the prediction step for the implicit
    ! enthalpy update at a given box on a given level 
    ! -----------------------------------------------------------------
    subroutine setup_box(lev, mfi, geom, phi_tmp, temp_tmp, alpha, beta, rhs)
  
      use amr_data_module, only : idomain
      use material_properties_module, only : get_temp
      ! use heat_transfer_domain_module, only : set_physical_boundary
      
      ! Input and output variables
      integer, intent(in) :: lev
      type(amrex_mfiter), intent(in) :: mfi
      type(amrex_geometry), intent(in) :: geom
      type(amrex_multifab), intent(inout) :: phi_tmp
      type(amrex_multifab), intent(inout) :: temp_tmp
      type(amrex_multifab), intent(inout) :: alpha
      type(amrex_multifab), intent(inout) :: beta(amrex_spacedim)
      type(amrex_multifab), intent(inout) :: rhs
  
      ! Local variables
      integer :: idim
      real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: penth_tmp
      real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptemp_tmp
      real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidout
      real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pac
      real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pbc
      real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: prhs
      ! real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: p
      type(amrex_box) :: bx
      
      ! Box
      bx = mfi%validbox()   
            
      ! Pointers
      ptemp_tmp   => temp_tmp%dataptr(mfi)
      penth_tmp   => phi_tmp%dataptr(mfi)
      pidout  => idomain(lev)%dataptr(mfi)
      pac     => alpha%dataptr(mfi)
      prhs    => rhs%dataptr(mfi)
      
      ! Fill in the temperature multifab with ghost points
      call get_temp(lbound(ptemp_tmp), ubound(ptemp_tmp), &
                    penth_tmp, lbound(penth_tmp), ubound(penth_tmp), &
                    ptemp_tmp, lbound(ptemp_tmp), ubound(ptemp_tmp),.true.)


      ! ! Check if this box contains domain ghost points, if so fill the in
      ! p => mf%dataptr(mfi)
      ! if (.not. geom%domain%contains(p)) then ! part of this box is outside the domain
      !    call set_physical_boundary(geom, ptemp_tmp, lbound(ptemp_tmp), ubound(ptemp_tmp))
      ! end if

      ! Get alpha matrix for the linear solver (first term on left hand side)
      call get_alpha(bx%lo, bx%hi, pac, lbound(pac), ubound(pac))
  
      
      ! Get beta matrix for the linear solver (second term on left hand side)
      do idim = 1,amrex_spacedim
         pbc => beta(idim)%dataptr(mfi)
         call get_beta(bx%lo, bx%hi, idim, &
                       pidout, lbound(pidout), ubound(pidout), &
                       ptemp_tmp, lbound(ptemp_tmp), ubound(ptemp_tmp), &
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
    ! of the auxillary potential
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
      call abeclap%set_domain_bc([amrex_lo_dirichlet, amrex_lo_neumann], &
                                 [amrex_lo_neumann, amrex_lo_neumann])
      ! call abeclap%set_domain_bc([amrex_lo_neumann, amrex_lo_dirichlet], &
      !                            [amrex_lo_neumann, amrex_lo_neumann])
      
      ! Coarse-fine interface  (for problem with pure homogeneous Neumann BC we can pass an empty multifab)
      if (lev > 0) then
         ! use coarse level data to set up bc at corase/fine boundary
         call abeclap % set_coarse_fine_bc(psi(lev-1), amrex_ref_ratio(lev-1))
      end if
      !call abeclap%set_level_bc(0, psi(lev))
      call abeclap%set_level_bc(0, nullmf)
      
      ! Set A and B scalar
      call abeclap%set_scalars(1.0_amrex_real, -1.0_amrex_real)
      
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

    ! -----------------------------------------------------------------------------------
    ! Subroutine used to get the averaged current within the melt pool
    ! ----------------------------------------------------------------------------------- 
    subroutine get_current(temp_tmp, lev)

      use amr_data_module, only : psi, &
                                  idomain, &
                                  surf_accum_current_x, &
                                  surf_accum_current_y, &
                                  surf_ind
      
      ! Input and output variables
      type(amrex_multifab), intent(in) :: temp_tmp
      integer, intent(in) :: lev

      ! Local variables
      real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptemp_tmp
      real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ppsi
      real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pid
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(amrex_real) :: surf_box_current_x(surf_ind(1,1):surf_ind(1,2),1:2)
      real(amrex_real) :: surf_box_current_y(surf_ind(1,1):surf_ind(1,2),1:2)
      type(amrex_geometry) :: geom

      geom = amrex_geom(lev)

      surf_accum_current_x = 0.0
      surf_accum_current_y = 0.0
      !$omp parallel private(mfi, bx, pid, ptemp_tmp, ppsi, surf_box_current_x, surf_box_current_y)
      ! Loop through the boxes on the maximum level
      call amrex_mfiter_build(mfi, idomain(lev), tiling=.false.)
      do while(mfi%next())
         surf_box_current_x = 0.0
         surf_box_current_y = 0.0
      
         ! Box
         bx = mfi%validbox()   
         
         ! Pointers
         ptemp_tmp => temp_tmp%dataptr(mfi)
         pid   => idomain(lev)%dataptr(mfi)
         ppsi   => psi(lev)%dataptr(mfi)
         
         ! Terms that depend on the temperature
         call get_current_box(bx%lo, bx%hi, geom%dx, &
                              ptemp_tmp, lbound(ptemp_tmp), ubound(ptemp_tmp), &
                              pid, lbound(pid), ubound(pid), &
                              ppsi, lbound(ppsi), ubound(ppsi), &
                              surf_box_current_x(bx%lo(1):bx%hi(1),:), &
                              surf_box_current_y(bx%lo(1):bx%hi(1),:))
        
         !$omp critical
         surf_accum_current_x = surf_accum_current_x + surf_box_current_x
         surf_accum_current_y = surf_accum_current_y + surf_box_current_y
         !$omp end critical
      
      end do
      !$omp end parallel   
      call amrex_mfiter_destroy(mfi)

    end subroutine get_current

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the current within one amrex box
  ! -----------------------------------------------------------------
  subroutine get_current_box(lo, hi, dx, &
                             temp, temp_lo, temp_hi, &
                             idom, id_lo, id_hi, &
                             psi, psi_lo, psi_hi, &
                             surf_box_current_x, surf_box_current_y)
    
    use amr_data_module, only : surf_current
    use material_properties_module, only : get_electrical_resistivity
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2) ! Bounds of the current tile box
    real(amrex_real), intent(in) :: dx(2)
    integer, intent(in) :: temp_lo(2), temp_hi(2) ! Bounds of temperature box
    integer, intent(in) :: id_lo(2), id_hi(2) ! Bounds of the idomain box
    integer, intent(in) :: psi_lo(2), psi_hi(2) ! Bounds of electrostatic potential box
    real(amrex_real), intent(in) :: temp(temp_lo(1):temp_hi(1),temp_lo(2):temp_hi(2))
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2))
    real(amrex_real), intent(in) :: psi(psi_lo(1):psi_hi(1),psi_lo(2):psi_hi(2))
    real(amrex_real), intent(out) :: surf_box_current_x(lo(1):hi(1),1:2)
    real(amrex_real), intent(out) :: surf_box_current_y(lo(1):hi(1),1:2)

    ! Local variables
    integer :: i,j
    real(amrex_real) :: rho_e
    real(amrex_real) :: temp_face
    real(amrex_real) :: Jx_face
    real(amrex_real) :: Jy_face
   !  integer :: molten_cells

    surf_box_current_x = 0.0
    surf_box_current_y = 0.0

    ! Loop through the box
    ! x - component
    do j = lo(2), hi(2)
      do i = lo(1),hi(1)
         ! Catch case of motlen elements with no neighboring material cells
         if(nint(idom(i-1,j)).le.0 .or. nint(idom(i,j)).le.0) then

           Jx_face =  0.0_amrex_real
         !   molten_cells = min(nint((surf_pos(i-1)-melt_pos(i-1))/dx(2)), nint((surf_pos(i)-melt_pos(i))/dx(2)))
         !    ! molten_cells = nint((surf_pos(i)-melt_pos(i))/dx(2))
         !   surf_box_current_x(i) = surf_box_current_x(i) +  1.0/molten_cells*Jx_face
         elseif(nint(idom(i,j-1)).ge.2 .or. nint(idom(i,j)).ge.2) then
           temp_face = (temp(i-1,j)+temp(i,j))/2
           call get_electrical_resistivity(temp_face, rho_e)
           Jx_face =  - 1/(rho_e*dx(1))*(psi(i,j)-psi(i-1,j))
           
           ! The number of molten cells
         !   molten_cells = min(nint((surf_pos(i-1)-melt_pos(i-1))/dx(2)), nint((surf_pos(i)-melt_pos(i))/dx(2)))
         !    ! molten_cells = nint((surf_pos(i)-melt_pos(i))/dx(2))
           surf_box_current_x(i,1) = surf_box_current_x(i,1) + Jx_face
           surf_box_current_x(i,2) = surf_box_current_x(i,2) + 1.0
         end if
      end do
   end do   

    ! y - component
    do j = lo(2), hi(2)
       do i = lo(1),hi(1)
          if(nint(idom(i,j-1)).ge.2 .and. nint(idom(i,j)).eq.0) then

            Jy_face =  - surf_current(i)
            ! molten_cells = nint((surf_pos(i)-melt_pos(i))/dx(2))+1
            ! surf_box_current_y(i) = surf_box_current_y(i) +  1.0/molten_cells*Jy_face
            surf_box_current_y(i,1) = surf_box_current_y(i,1) +  Jy_face
            surf_box_current_y(i,2) = surf_box_current_y(i,2) +  1.0
          elseif(nint(idom(i,j)).ge.2 .and. nint(idom(i,j-1)).ge.1) then

            temp_face = (temp(i,j-1)+temp(i,j))/2
            call get_electrical_resistivity(temp_face, rho_e)
            Jy_face =  - 1/(rho_e*dx(2))*(psi(i,j)-psi(i,j-1))
            
            ! molten_cells = nint((surf_pos(i)-melt_pos(i))/dx(2))+1
            ! surf_box_current_y(i) = surf_box_current_y(i) + 1.0/molten_cells*Jy_face
            surf_box_current_y(i,1) = surf_box_current_y(i,1) + Jy_face
            surf_box_current_y(i,2) = surf_box_current_y(i,2) + 1.0
          end if
       end do
    end do   
    
  end subroutine get_current_box    

    ! -----------------------------------------------------------------
    ! Subroutine used to advance of one time step the conduction part
    ! of the heat transfer equation
    ! -----------------------------------------------------------------
    subroutine free_memory(phi_tmp, temp_tmp, rhs, alpha, beta, dm)
      
      ! Input and output variables
      type(amrex_multifab) :: phi_tmp
      type(amrex_multifab) :: temp_tmp
      type(amrex_multifab) :: rhs
      type(amrex_multifab) :: alpha
      type(amrex_multifab) :: beta(amrex_spacedim) 
      type(amrex_distromap):: dm
      
      ! Local variables
      integer :: idim
  
      call amrex_multifab_destroy(phi_tmp)
      call amrex_multifab_destroy(temp_tmp)
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

                   call get_electrical_resistivity(temp_melt, beta(i,j))
                   beta(i,j) = 1E-16/beta(i,j)

                  !  beta(i,j) = 0_amrex_real ! Suppress flux at the free surface and cooling pipe surface
                else
                  !  temp_face = (temp(i,j) + temp(i-1,j))/2_amrex_real
                  !  if (temp_face.lt.2900.0) temp_face = 2900.0
                  !  call get_electrical_resistivity(temp_face, beta(i,j))

                  if (temp(i,j).lt.300) then
                     temp_face = temp_melt
                  else
                     temp_face = temp(i,j)
                  end if
                  call get_electrical_resistivity(temp_face, beta(i,j))

                  !  call get_electrical_resistivity(temp_melt, beta(i,j))

                   beta(i,j) = 1.0/beta(i,j)
               end if
               
            end do
         end do
  
      else if (dim == 2) then ! y-direction
         do j = lo(2), hi(2)+1
            do i = lo(1), hi(1)
          
               if(nint(idom(i,j-1)).eq.0 .or. nint(idom(i,j)).eq.0 .or. &
                    nint(idom(i,j-1)).eq.-1 .or. nint(idom(i,j)).eq.-1) then

               ! if(nint(idom(i,j)).ne.0 .and. nint(idom(i,j+1)).eq.0) then

                  call get_electrical_resistivity(temp_melt, beta(i,j))
                  beta(i,j) = 1E-16/beta(i,j)

                  ! beta(i,j) = 0.0
                  
               else

                  if (temp(i,j).lt.300.0) then
                     temp_face = temp_melt
                  else
                     temp_face = temp(i,j)
                  end if
                  call get_electrical_resistivity(temp_face, beta(i,j))

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
      ! use read_input_module, only : sw_current
      
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
      type(amrex_geometry) :: geom

      geom = amrex_geom(0)
      
      ! if (geom%domain%hi(2).lt.j+1) then
      !    rhs(i,j) = 100.0
      ! end if
      ! Fill rhs of linear problem
      do j = lo(2), hi(2)          
         do i = lo(1), hi(1)
               ! rhs(i,j) = 0.0_amrex_real
            if(nint(idom(i,j)).ne.0 .and. nint(idom(i,j+1)).eq.0) then
               call interp_to_max_lev(lo, lo_phys, dx, i, xind)
               rhs(i,j) = -surf_current(xind)/dx(2)
            else
               rhs(i,j) = 0.0_amrex_real
            end if
            ! if (geom%domain%hi(2).lt.j+1) then
            !    rhs(i,j) = -1E11
            ! else
            !    rhs(i,j) = 0.0_amrex_real
            ! end if
            ! if (geom%domain%hi(2).lt.j+1) then
            !    call interp_to_max_lev(lo, lo_phys, dx, i, xind)
            !    rhs(i,j) = -surf_current(xind)/dx(2)
            ! else
            !    rhs(i,j) = 0.0_amrex_real
            ! end if
         end do
      end do
  
    end subroutine get_rhs
  
      
  end module electrostatics_module