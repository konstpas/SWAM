module heat_transfer_module
  
  ! -----------------------------------------------------------------
  ! This module is used to perform all the calculations relative
  ! to the heat transfer part of the code.
  ! -----------------------------------------------------------------
  
  use amrex_amr_module
  
  implicit none

  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: advance_heat_solver_explicit_level
  public :: advance_heat_solver_implicit
  
contains

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy at a new time step for
  ! a given level via an explicit update
  ! -----------------------------------------------------------------
  subroutine advance_heat_solver_explicit_level(lev, time, dt, substep)

    use read_input_module, only : do_reflux, temp_fs
    use amr_data_module, only : phi_new, phi_old, temp, idomain, flux_reg
    use domain_module, only : get_idomain, get_melt_pos, reset_melt_pos
    use regrid_module, only : fillpatch
    use material_properties_module, only : get_temp
    
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

    ! Swap old and new enthalpies before updating solution
    call amrex_multifab_swap(phi_old(lev), phi_new(lev))
    
    ! Build enthalpy and temperature multifabs with ghost points
    call amrex_multifab_build(phi_tmp, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, nghost) 
    call amrex_multifab_build(temp_tmp, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, nghost)
    call fillpatch(lev, time, phi_tmp)

    ! Build temporary idomain multifab to store the domain configuration before SW deformation
    call amrex_multifab_build(idomain_tmp, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, nghost)
    call amrex_multifab_swap(idomain_tmp, idomain(lev))

    ! Initialize fluxes  
    if (do_reflux) then
       do idim = 1, amrex_spacedim
          nodal = .false.
          nodal(idim) = .true.
          call amrex_multifab_build(fluxes(idim), phi_new(lev)%ba, phi_new(lev)%dm, ncomp, 0, nodal)
       end do       
    end if
    do idim = 1, amrex_spacedim
       call flux(idim)%reset_omp_private()
    end do
       
    ! Advance heat solver on all boxes on the level
    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.false.)  
    						 
    do while(mfi%next())

       ! Box
       bx = mfi%validbox()   

       ! Pointers
       pin     => phi_tmp%dataptr(mfi)
       pout    => phi_new(lev)%dataptr(mfi)
       ptempin => temp_tmp%dataptr(mfi)
       ptemp   => temp(lev)%dataptr(mfi)
       pidin => idomain_tmp%dataptr(mfi)
       pidout => idomain(lev)%dataptr(mfi)   
       do idim = 1, amrex_spacedim
          tbx = bx
          call tbx%nodalize(idim)
          call flux(idim)%resize(tbx,ncomp)
          call tbx%grow(substep)
       end do
       pfx => flux(1)%dataptr()
       pfy => flux(2)%dataptr()

       
       ! Get temperature corresponding to the enthalpy (with ghost points)
       call get_temp(lbound(ptempin), ubound(ptempin), &
                     pin, lbound(pin), ubound(pin), &
                     ptempin, lbound(ptempin), ubound(ptempin),.true.)
       
       ! Get configuration of the system after the deformation
       call get_idomain(geom%get_physical_location(bx%lo), geom%dx, &
                        bx%lo, bx%hi, &
                        pidout, lbound(pidout), ubound(pidout), &
                        ptempin, lbound(ptempin), ubound(ptempin))

       ! Re-evaluate heat transfer domain after the deformation
       call revaluate_heat_domain(bx%lo, bx%hi, &
                                  pidin, lbound(pidin), ubound(pidin), &
                                  pidout, lbound(pidout), ubound(pidout), &
                                  pin, lbound(pin), ubound(pin))

       
       ! Increment enthalpy at given box depending on the condition of the free surface
       if (temp_fs.gt.0) then
          call get_enthalpy_explicit_fixT(bx%lo, bx%hi, &
                                          pin, lbound(pin),     ubound(pin),     &
                                          pout,    lbound(pout),    ubound(pout),    &
                                          ptempin, lbound(ptempin), ubound(ptempin), &
                                          pfx, lbound(pfx), ubound(pfx), &
                                          pfy, lbound(pfy), ubound(pfy), &
                                          pidout, lbound(pidout), ubound(pidout), &
                                          geom, dt)
       else
          call get_enthalpy_explicit(time, bx%lo, bx%hi, &
                                     pin, lbound(pin),     ubound(pin),     &
                                     pout,    lbound(pout),    ubound(pout),    &
                                     ptempin, lbound(ptempin), ubound(ptempin), &
                                     pfx, lbound(pfx), ubound(pfx), &
                                     pfy, lbound(pfy), ubound(pfy), &
                                     pidout, lbound(pidout), ubound(pidout), &
                                     geom, dt)
       end if

       ! Get temperature corresponding to the new enthalpy
       call get_temp(bx%lo, bx%hi, &
                     pout, lbound(pout), ubound(pout), &
                     ptemp, lbound(ptemp), ubound(ptemp), .true.) 
       
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

    ! Update flux registers (fluxes have already been scaled by dt and area in the advance_heat_solver_box subroutine)
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
    call amrex_mfiter_destroy(mfi)
    do idim = 1, amrex_spacedim
       call amrex_fab_destroy(flux(idim))
    end do
    call amrex_multifab_destroy(phi_tmp)
    call amrex_multifab_destroy(temp_tmp)
    call amrex_multifab_destroy(idomain_tmp)
    
  end subroutine advance_heat_solver_explicit_level

  
  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy at a new time step for
  ! a given box of a certain level via an explicit update
  ! -----------------------------------------------------------------
  subroutine get_enthalpy_explicit(time, lo, hi, &
                                   u_old,  uo_lo, uo_hi, &
                                   u_new, un_lo, un_hi, &
                                   temp, t_lo, t_hi, &
                                   flxx, fx_lo, fx_hi, &
                                   flxy, fy_lo, fy_hi, &
                                   idom, id_lo, id_hi, &
                                   geom, dt)

    use heat_flux_module, only: get_boundary_heat_flux
    use material_properties_module, only : get_temp
 
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2) ! bounds of current tile box
    integer, intent(in) :: uo_lo(2), uo_hi(2) ! bounds of input enthalpy box 
    integer, intent(in) :: un_lo(2), un_hi(2) ! bounds of output enthalpy box  
    integer, intent(in) :: t_lo(2), t_hi(2) ! bounds of the temperature box  
    integer, intent(in) :: fx_lo(2), fx_hi(2) ! bounds of the enthalpy flux along x
    integer, intent(in) :: fy_lo(2), fy_hi(2) ! bounds of the enthalpy flux along y
    integer, intent(in) :: id_lo(2), id_hi(2) ! bounds of the idomain box
    real(amrex_real), intent(in) :: dt ! time step
    real(amrex_real), intent(in) :: time ! time
    real(amrex_real), intent(in) :: u_old(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2)) ! Input enthalpy 
    real(amrex_real), intent(inout) :: u_new(un_lo(1):un_hi(1),un_lo(2):un_hi(2)) ! Output enthalpy
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2)) ! Temperature
    real(amrex_real), intent(out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2)) ! flux along the x direction  			
    real(amrex_real), intent(out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2)) ! flux along the y direction
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2)) ! idomain
    type(amrex_geometry), intent(in) :: geom ! geometry
    
    !Local variables
    integer :: i,j
    real(amrex_real) :: dx(2) ! Grid size
    real(amrex_real) :: lo_phys(2) ! Physical location of the lowest corner of the tile box
    real(amrex_real) :: qbound(lo(1):hi(1),lo(2):hi(2)) ! Volumetric heating (boundary)

    ! Get grid size
    dx = geom%dx(1:2) ! grid width at level 

    ! Get physical location of the lowest corner of the tile box
    lo_phys = geom%get_physical_location(lo)
    
    ! Get enthalpy flux 
    call create_face_flux(dx, lo, hi, &
                          u_old, uo_lo, uo_hi, &
                          flxx, fx_lo, fx_hi, &
                          flxy, fy_lo, fy_hi, &
                          temp, t_lo, t_hi, &
                          idom, id_lo, id_hi)
  				  	
    ! Prescribe external heat flux on the free surface
    call get_boundary_heat_flux(time, lo_phys, &
                                dx, lo, hi, &
                                idom, id_lo, id_hi, &
                                temp, t_lo, t_hi, qbound)

    ! Compute enthalpy at the new timestep
    do i = lo(1),hi(1)
       do  j = lo(2),hi(2)
          u_new(i,j) = u_old(i,j) &
               - dt/dx(1) * (flxx(i+1,j) - flxx(i,j)) & ! flux divergence x-direction 
               - dt/dx(2) * (flxy(i,j+1) - flxy(i,j)) & ! flux divergence y-direction 
               + dt*qbound(i,j) ! 'boundary volumetric' source
       end do
    end do
    
  	
    ! Scale the fluxes for the flux registers
    do j = lo(2), hi(2)
       do i = lo(1), hi(1) + 1
          flxx(i,j) = flxx(i,j) * (dt * dx(2))
       end do
    end do

    do j = lo(2), hi(2) + 1
       do i = lo(1), hi(1)
          flxy(i,j) = flxy(i,j) * (dt * dx(1))
       end do
    end do
    
  end subroutine get_enthalpy_explicit


  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy at a new time step for
  ! a given box on a certain level via an explicit update.
  ! Case with fixed temperature at the free surface
  ! -----------------------------------------------------------------
  subroutine get_enthalpy_explicit_fixT(lo, hi, &
                                        u_old,  uo_lo, uo_hi, &
                                        u_new, un_lo, un_hi, &
                                        temp, t_lo, t_hi, &
                                        flxx, fx_lo, fx_hi, &
                                        flxy, fy_lo, fy_hi, &
                                        idom, id_lo, id_hi, &
                                        geom, dt)
    
    use material_properties_module, only : get_temp, get_enthalpy
    use read_input_module, only : temp_fs
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2) ! bounds of current tile box
    integer, intent(in) :: uo_lo(2), uo_hi(2) ! bounds of input enthalpy box 
    integer, intent(in) :: un_lo(2), un_hi(2) ! bounds of output enthalpy box  
    integer, intent(in) :: t_lo(2), t_hi(2) ! bounds of the temperature box  
    integer, intent(in) :: fx_lo(2), fx_hi(2) ! bounds of the enthalpy flux along x
    integer, intent(in) :: fy_lo(2), fy_hi(2) ! bounds of the enthalpy flux along y
    integer, intent(in) :: id_lo(2), id_hi(2) ! bounds of the idomain box
    real(amrex_real), intent(in) :: dt ! time step
    real(amrex_real), intent(in) :: u_old(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2)) ! Input enthalpy
    real(amrex_real), intent(inout) :: u_new(un_lo(1):un_hi(1),un_lo(2):un_hi(2)) ! Output enthalpy
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2)) ! Temperature
    real(amrex_real), intent(out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2)) ! flux along the x direction  			
    real(amrex_real), intent(out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2)) ! flux along the y direction
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2)) ! idomain
    type(amrex_geometry), intent(in) :: geom ! geometry
    
    !Local variables
    integer :: i,j
    real(amrex_real) :: u_fs ! Enthalpy of the points on the free surface
    real(amrex_real) :: dx(2) ! Grid size
    real(amrex_real) :: lo_phys(2) ! Physical location of the lowest corner of the tile box

    ! Enthalpy at the free surface
    call get_enthalpy(temp_fs, u_fs)
    
    ! Get grid size
    dx = geom%dx(1:2) ! grid width at level 

    ! Get physical location of the lowest corner of the tile box
    lo_phys = geom%get_physical_location(lo)
    
    ! Get enthalpy flux 
    call create_face_flux_fixT(dx, lo, hi, &
                               flxx, fx_lo, fx_hi, &
                               flxy, fy_lo, fy_hi, &
                               temp, t_lo, t_hi)

    ! Compute enthalpy at the new timestep
    do i = lo(1),hi(1)
       do  j = lo(2),hi(2)
          
          if (nint(idom(i,j)).ne.0 .and. nint(idom(i,j+1)).eq.0) then
             u_new(i,j) = u_fs ! Impose temperature on the free surface
          else if(nint(idom(i,j)).eq.0) then
             u_new(i,j) = u_fs ! Set background temperature equal to the free surface temperature
          else
             ! Update enthalpy according to heat equation
             u_new(i,j) = u_old(i,j) &
                          - dt/dx(1) * (flxx(i+1,j) - flxx(i,j)) & ! flux divergence x-direction 
                          - dt/dx(2) * (flxy(i,j+1) - flxy(i,j)) ! flux divergence y-direction
          end if 

       end do
    end do
    
  	
    ! Scale the fluxes for the flux registers
    do j = lo(2), hi(2)
       do i = lo(1), hi(1) + 1
          flxx(i,j) = flxx(i,j) * (dt * dx(2))
       end do
    end do

    do j = lo(2), hi(2) + 1
       do i = lo(1), hi(1)
          flxy(i,j) = flxy(i,j) * (dt * dx(1))
       end do
    end do
    
  end subroutine get_enthalpy_explicit_fixT


  ! -----------------------------------------------------------------
  ! Subroutine used to re-evaluate the heat equation domain
  ! -----------------------------------------------------------------  
  subroutine revaluate_heat_domain(lo, hi, &
                                   idom_old, ido_lo, ido_hi, &
                                   idom_new, idn_lo, idn_hi, &
                                   u_in, u_lo, u_hi)


    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2) ! bounds of current tile box
    integer, intent(in) :: u_lo(2), u_hi(2) ! bounds of input enthalpy box 
    integer, intent(in) :: ido_lo(2), ido_hi(2) ! bounds of the input idomain box
    integer, intent(in) :: idn_lo(2), idn_hi(2) ! bounds of the output idomain box
    real(amrex_real), intent(inout) :: u_in(u_lo(1):u_hi(1),u_lo(2):u_hi(2)) ! Input enthalpy 
    real(amrex_real), intent(in) :: idom_old(ido_lo(1):ido_hi(1),ido_lo(2):ido_hi(2))
    real(amrex_real), intent(in) :: idom_new(idn_lo(1):idn_hi(1),idn_lo(2):idn_hi(2))
    
    !Local variables
    integer :: i,j
    
    do i = lo(1)-1,hi(1)+1
       do  j = lo(2)-1,hi(2)+1

          ! Points added to the domain
          if (nint(idom_old(i,j)).eq.0 .and. nint(idom_new(i,j)).ne.0) then
             u_in(i,j) = u_in(i,j-1)
          ! Points removed from the domain
          else if (nint(idom_new(i,j)).eq.0) then
             u_in(i,j) = 0_amrex_real
          end if
          
       end do
    end do
    
  end subroutine revaluate_heat_domain
  
  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy fluxes on the edges of
  ! the grid
  ! -----------------------------------------------------------------  
  subroutine create_face_flux(dx, lo, hi, &
                              u_old, uo_lo, uo_hi, &
                              flxx, fx_lo, fx_hi, &
                              flxy, fy_lo, fy_hi, &
                              temp, t_lo, t_hi, &
                              idom, id_lo, id_hi)
  				
    use material_properties_module, only: get_conductivity

    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)  
    integer, intent(in) :: uo_lo(2), uo_hi(2)
    integer, intent(in) :: t_lo(2), t_hi(2)
    integer, intent(in) :: fx_lo(2), fx_hi(2)
    integer, intent(in) :: fy_lo(2), fy_hi(2)
    integer, intent(in) :: id_lo(2), id_hi(2)
    real(amrex_real), intent(in) :: dx(2)    
    real(amrex_real), intent(in) :: u_old(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2))
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2)) 
    real(amrex_real), intent(out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))
    real(amrex_real), intent(out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2))
    real(amrex_real), intent(in) :: temp (t_lo(1):t_hi(1),t_lo(2):t_hi(2))

    ! Local variables
    integer :: i,j
    real(amrex_real) :: ktherm
    real(amrex_real) :: temp_face
    real(amrex_real) :: vx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))

    ! Construct 3D melt velocity profile from the 2D shallow water solution
    call get_face_velocity(lo, hi, &
                           vx, fx_lo, fx_hi, &
                           temp, t_lo, t_hi)
    
    ! Flux along the x direction
    do i = lo(1), hi(1)+1
       do j = lo(2), hi(2)          

          if (nint(idom(i-1,j)).eq.0 .or. nint(idom(i,j)).eq.0) then

             ! Suppress flux at the free surface
             flxx(i,j) = 0_amrex_real

          else
             
             ! Advective component
             if (vx(i,j) > 0_amrex_real) then 
                flxx(i,j)  = u_old(i-1,j)*vx(i,j)
             else 
                flxx(i,j)  = u_old(i,j)*vx(i,j)
             end if
             
             ! Diffusive component
             temp_face = (temp(i,j) + temp(i-1,j))/2_amrex_real
             call get_conductivity(temp_face, ktherm)
             flxx(i,j) = flxx(i,j) - ktherm*(temp(i,j)-temp(i-1,j))/dx(1)
             
          end if
          
       end do
    end do
    
    ! Flux along the y direction
    do i = lo(1), hi(1)
       do j = lo(2), hi(2)+1

          if (nint(idom(i,j-1)).eq.0 .or. nint(idom(i,j)).eq.0) then

             ! Suppress flux at the free surface
             flxy(i,j) = 0_amrex_real

          else
             
             ! Diffusive component (there is no advection in the y direction)
             temp_face = (temp(i,j) + temp(i,j-1))/2_amrex_real
             call get_conductivity(temp_face, ktherm)
             flxy(i,j) = -ktherm*(temp(i,j)-temp(i,j-1))/dx(2)
             
          end if
          
       end do
    end do
 
  end subroutine create_face_flux
  

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy fluxes on the edges of
  ! the grid. Case with fixed temperature at the free surface
  ! -----------------------------------------------------------------  
  subroutine create_face_flux_fixT(dx, lo, hi, &
                                   flxx, fx_lo, fx_hi, &
                                   flxy, fy_lo, fy_hi, &
                                   temp, t_lo, t_hi)
  				
    use material_properties_module, only: get_conductivity

    ! Input and output variables 
    integer, intent(in) :: lo(2), hi(2)  
    integer, intent(in) :: t_lo(2), t_hi(2)
    integer, intent(in) :: fx_lo(2), fx_hi(2)
    integer, intent(in) :: fy_lo(2), fy_hi(2)
    real(amrex_real), intent(in) :: dx(2) 
    real(amrex_real), intent(out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))
    real(amrex_real), intent(out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2))
    real(amrex_real), intent(in) :: temp (t_lo(1):t_hi(1),t_lo(2):t_hi(2))

    ! Local variables
    integer :: i,j
    real(amrex_real) :: ktherm
    real(amrex_real) :: temp_face

 
    
    ! Flux along the x direction
    do i = lo(1), hi(1)+1
       do j = lo(2), hi(2)          

          ! Diffusive component
          temp_face = (temp(i,j) + temp(i-1,j))/2_amrex_real
          call get_conductivity(temp_face, ktherm)
          flxx(i,j) = -ktherm*(temp(i,j)-temp(i-1,j))/dx(1)

       end do
    end do
    
    ! Flux along the y direction
    do i = lo(1), hi(1)
       do j = lo(2), hi(2)+1
          
          ! Diffusive component (there is no advection in the y direction)
          temp_face = (temp(i,j) + temp(i,j-1))/2_amrex_real
          call get_conductivity(temp_face, ktherm)
          flxy(i,j) = -ktherm*(temp(i,j)-temp(i,j-1))/dx(2)
          
       end do
    end do
 
  end subroutine create_face_flux_fixT
  
  ! -----------------------------------------------------------------
  ! Subroutine used to the velocity on the faces of each grid cell.
  ! This subroutine translates to 3D the 2D velocity field obtained
  ! from the solution of the shallow water equations
  ! -----------------------------------------------------------------  
  subroutine get_face_velocity(lo, hi, &
                               vx, vx_lo, vx_hi, &
                               temp, t_lo, t_hi)

    use material_properties_module, only : temp_melt
    use read_input_module, only : meltvel
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: t_lo(2), t_hi(2)
    integer, intent(in) :: vx_lo(2), vx_hi(2)
    real(amrex_real), intent(in)  :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2))
    real(amrex_real), intent(out) :: vx(vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2))

    ! Local variables
    integer :: i,j
    
    ! THIS MUST BE UPDATED 
    do i = lo(1), hi(1)+1
       do j = lo(2), hi(2) 
          if (temp(i,j).gt.temp_melt) then
             vx(i,j) = meltvel
          else
             vx(i,j) = 0_amrex_real
          end if
       end do   
    end do

    
  end subroutine get_face_velocity

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy at a new time step for
  ! all levels via an implicit update.
  ! -----------------------------------------------------------------
  subroutine advance_heat_solver_implicit(time, dt)

    use amr_data_module, only : phi_new, phi_old, temp, idomain
    use regrid_module, only : fillpatch
    use domain_module, only : get_idomain, get_melt_pos, reset_melt_pos
    use material_properties_module, only : get_temp

    ! Input and output variables
    real(amrex_real), intent(in) :: dt
    real(amrex_real), intent(in) :: time
  
    ! Local variables
    logical :: nodal(2)
    integer, parameter :: nghost = 1 ! number of ghost points in each spatial direction 
    integer :: ncomp
    integer :: idim
    integer :: ilev
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pin
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pout
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptempin
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptemp
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pfx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pfy
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidin
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidout
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pac
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pbc
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: prhs
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ppstar
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptempstar
    type(amrex_multifab) :: phi_tmp(0:amrex_max_level) ! Enthalpy multifab with ghost points
    type(amrex_multifab) :: temp_tmp(0:amrex_max_level) ! Temperature multifab with ghost points
    type(amrex_multifab) :: idomain_tmp(0:amrex_max_level) ! Idomain multifab to distinguish material and background
    type(amrex_multifab) :: rhs(0:amrex_max_level) ! multifab for the right hand side of the linear system of equations
    type(amrex_multifab) :: alpha(0:amrex_max_level) ! multifab for the first term on the left hand side of the linear system of equations
    type(amrex_multifab) :: beta(amrex_spacedim,0:amrex_max_level) ! multifab for the second term on the left hand side of the linear system of equations
    type(amrex_multifab) :: phi_star(0:amrex_max_level) ! Enthalpy multifab for predictor-corrector procedure
    type(amrex_multifab) :: temp_star(0:amrex_max_level) ! Temperature multifab for predictor-corrector procedure
    type(amrex_mfiter) :: mfi ! Multifab iterator
    type(amrex_box) :: bx
    type(amrex_box) :: tbx
    type(amrex_boxarray) :: ba(0:amrex_max_level) ! Box array
    type(amrex_distromap):: dm(0:amrex_max_level) ! Distribution mapping
    type(amrex_fab) :: flux(amrex_spacedim)

        
    ! Loop through all the levels
    do ilev = 0, amrex_max_level

       ! Swap old and new enthalpies before updating solution
       call amrex_multifab_swap(phi_old(ilev), phi_new(ilev))
       
       ! Get boxarray and distribution mapping on level
       ba(ilev) = phi_new(ilev)%ba
       dm(ilev) = phi_new(ilev)%dm

       ! Number of components
       ncomp = phi_new(ilev)%ncomp()
       
       ! Enthalpy and temperature multifabs with ghost points
       call amrex_multifab_build(phi_tmp(ilev), ba(ilev), dm(ilev), ncomp, nghost) 
       call amrex_multifab_build(temp_tmp(ilev), ba(ilev), dm(ilev), ncomp, nghost)
       call fillpatch(ilev, time, phi_tmp(ilev))

       ! Enthalpy and temperature multifabs for predictor-corrector procedure
       call amrex_multifab_build(phi_star(ilev), ba(ilev), dm(ilev), ncomp, 0)
       call amrex_multifab_build(temp_star(ilev), ba(ilev), dm(ilev), ncomp, 0)
       
       ! Build temporary idomain multifab to store the domain configuration before SW deformation
       call amrex_multifab_build(idomain_tmp(ilev), ba(ilev), dm(ilev), ncomp, nghost)
       call amrex_multifab_swap(idomain_tmp(ilev), idomain(ilev))

       ! Initialize fluxes (used in the predictor step)
       do idim = 1, amrex_spacedim
          call flux(idim)%reset_omp_private()
       end do
    
       ! Multifabs for the linear solver
       call amrex_multifab_build(rhs(ilev), ba(ilev), dm(ilev), ncomp, 0)
       call amrex_multifab_build(alpha(ilev), ba(ilev), dm(ilev), ncomp, 0)
       do idim = 1, amrex_spacedim
          nodal = .false.
          nodal(idim) = .true.
          call amrex_multifab_build(beta(idim,ilev), ba(ilev), dm(ilev), ncomp, 0, nodal)
       end do
 
       ! Loop through all the boxes in the level
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
          pac     => alpha(ilev)%dataptr(mfi)
          prhs    => rhs(ilev)%dataptr(mfi)
          do idim = 1, amrex_spacedim
             tbx = bx
             call tbx%nodalize(idim)
             call flux(idim)%resize(tbx,ncomp)
             call tbx%grow(1)
          end do
          pfx => flux(1)%dataptr()
          pfy => flux(2)%dataptr()
          ppstar => phi_star(ilev)%dataptr(mfi)
          ptempstar => temp_star(ilev)%dataptr(mfi)
          
          ! Get temperature corresponding to the enthalpy (with ghost points)
          call get_temp(lbound(ptempin), ubound(ptempin), &
                        pin, lbound(pin), ubound(pin), &
                        ptempin, lbound(ptempin), ubound(ptempin), .true.)

          ! Get configuration of the system after the deformation
          call get_idomain(amrex_geom(ilev)%get_physical_location(bx%lo), amrex_geom(ilev)%dx, &
                           bx%lo, bx%hi, &
                           pidout, lbound(pidout), ubound(pidout), &
                           ptempin, lbound(ptempin), ubound(ptempin))

          ! Re-evaluate heat transfer domain after the deformation
          call revaluate_heat_domain(bx%lo, bx%hi, &
                                     pidin, lbound(pidin), ubound(pidin), &
                                     pidout, lbound(pidout), ubound(pidout), &
                                     pin, lbound(pin), ubound(pin))

          ! Get temperature corresponding to the enthalpy (prediction step)
          call get_temp(lbound(ptempin), ubound(ptempin), &
                        pin, lbound(pin), ubound(pin), &
                        ptempin, lbound(ptempin), ubound(ptempin), .true.)
          
          ! Explicit update of the enthalpy (prediction step)
          call get_enthalpy_explicit(time, bx%lo, bx%hi, &
                                     pin, lbound(pin),     ubound(pin),     &
                                     ppstar, lbound(ppstar), ubound(ppstar), &
                                     ptempin, lbound(ptempin), ubound(ptempin), &
                                     pfx, lbound(pfx), ubound(pfx), &
                                     pfy, lbound(pfy), ubound(pfy), &
                                     pidout, lbound(pidout), ubound(pidout), &
                                     amrex_geom(ilev), dt)
          
          ! Get temperature corresponding to the enthalpy (prediction step)
          call get_temp(lbound(ptempstar), ubound(ptempstar), &
                        ppstar, lbound(ppstar), ubound(ppstar), &
                        ptempstar, lbound(ptempstar), ubound(ptempstar), .true.)
          
          ! Get alpha matrix for the linear solver (first term on left hand side)
          call get_alpha(bx%lo, bx%hi, &
                         pin, lbound(pin), ubound(pin), &
                         ppstar, lbound(ppstar), ubound(ppstar), &
                         ptempin, lbound(ptempin), ubound(ptempin), &
                         ptempstar, lbound(ptempstar), ubound(ptempstar), &                         
                         pac, lbound(pac), ubound(pac))

          ! Get beta matrix for the linear solver (second term on left hand side)
          do idim = 1,amrex_spacedim
             pbc => beta(idim,ilev)%dataptr(mfi)
             call get_beta(bx%lo, bx%hi, idim, &
                           pidout, lbound(pidout), ubound(pidout), &
                           ptempin, lbound(ptempin), ubound(ptempin), &
                           pbc, lbound(pbc), ubound(pbc))
          end do

          ! Get right hand for the linear solver
          call get_rhs(bx%lo, bx%hi, time, dt, &
                       amrex_geom(ilev)%get_physical_location(bx%lo), amrex_geom(ilev)%dx, &
                       pidout, lbound(pidout), ubound(pidout), &
                       pin, lbound(pin), ubound(pin), &
                       ppstar, lbound(ppstar), ubound(ppstar), &
                       ptempin, lbound(ptempin), ubound(ptempin), &
                       ptempstar, lbound(ptempstar), ubound(ptempstar), &                      
                       prhs, lbound(prhs), ubound(prhs))

       end do

       call amrex_mfiter_destroy(mfi)

       
    end do
    

    ! Compute new temperature via implicit update
    call get_temperature_implicit(ba, dm, dt, alpha, beta, rhs, temp_tmp)

    
    ! ! Re-compute enthalpy and temperature on all levels (correction step)
    ! do ilev = 0, amrex_max_level
 
    !    ! Loop through all the boxes in the level
    !    call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.false.)  
    						 
    !    do while(mfi%next())

    !       ! Box
    !       bx = mfi%validbox()   
          
    !       ! Pointers
    !       pin     => phi_tmp(ilev)%dataptr(mfi)
    !       pout    => phi_new(ilev)%dataptr(mfi)
    !       ptempin => temp_tmp(ilev)%dataptr(mfi)
    !       ptemp   => temp(ilev)%dataptr(mfi)
    !       ppstar => phi_star(ilev)%dataptr(mfi)
    !       ptempstar => temp_star(ilev)%dataptr(mfi)

    !       ! Get corrected enthalpy
    !       call get_enthalpy_implicit(bx%lo, bx%hi, &
    !                                  pin, lbound(pin), ubound(pin), &
    !                                  ppstar, lbound(ppstar), ubound(ppstar), &
    !                                  pout, lbound(pout), ubound(pout), &
    !                                  ptempin, lbound(ptempin), ubound(ptempin), &
    !                                  ptempstar, lbound(ptempstar), ubound(ptempstar), &
    !                                  ptempin, lbound(ptempin), ubound(ptempin))
          
    !       ! Get corrected temperature 
    !       call get_temp(lbound(ptemp), ubound(ptemp), &
    !                     pout, lbound(pout), ubound(pout), &
    !                     ptemp, lbound(ptemp), ubound(ptemp), .true.)


    !    end do

    !    call amrex_mfiter_destroy(mfi)

    ! end do
    
    ! ! Find melt interface y position
    ! call amrex_mfiter_build(mfi, idomain(amrex_max_level), tiling=.false.)
    ! do while(mfi%next())

    !    ! Box
    !    bx = mfi%validbox()   
       
    !    ! Pointers
    !    pidin => idomain(amrex_max_level)%dataptr(mfi)
       
    !    call get_melt_pos(bx%lo, bx%hi, &
    !                      pidin, lbound(pidin), ubound(pidin), &
    !                      amrex_geom(amrex_max_level))
       
    ! end do
    ! call amrex_mfiter_destroy(mfi)    
        
    ! ! Clean memory
    ! do ilev = 0, amrex_max_level
    !    call amrex_multifab_destroy(phi_tmp(ilev))
    !    call amrex_multifab_destroy(temp_tmp(ilev))
    !    call amrex_multifab_destroy(idomain_tmp(ilev))
    !    call amrex_multifab_destroy(rhs(ilev))
    !    call amrex_multifab_destroy(alpha(ilev))
    !    do idim = 1, amrex_spacedim
    !       call amrex_multifab_destroy(beta(idim,ilev))
    !    end do
    !    call amrex_multifab_destroy(phi_star)
    !    call amrex_multifab_destroy(temp_star)
    !    call amrex_distromap_destroy(dm)
    ! end do

    
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
          pout      => phi_new(ilev)%dataptr(mfi)
          ptempin => temp_tmp(ilev)%dataptr(mfi)
          
          ! Map solution on temperature
          call get_temp_ls(bx%lo, bx%hi, &
                           ptemp, lbound(ptemp), ubound(ptemp), &
                           ptempin, lbound(ptempin), ubound(ptempin))

          ! Map temperature on enthalpy
          call get_temp(bx%lo, bx%hi, &
                        ptemp, lbound(ptemp), ubound(ptemp), &
                        pout, lbound(pout), ubound(pout), .false.)
                    
       end do
       call amrex_mfiter_destroy(mfi)
       
    end do

  end subroutine advance_heat_solver_implicit

  
  ! -----------------------------------------------------------------------------------
  ! ....
  ! ----------------------------------------------------------------------------------- 
  subroutine get_temperature_implicit(ba, dm, dt, alpha, beta, rhs, sol)

    use read_input_module, only : ls_agglomeration, &
                                  ls_consolidation, ls_max_coarsening_level, &
                                  ls_linop_maxorder, ls_bottom_solver, &
                                  ls_bottom_verbose, ls_max_fmg_iter, &
                                  ls_max_iter, ls_verbose
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
    call amrex_abeclaplacian_build(abeclap, amrex_geom, ba, dm, &
                                   metric_term=.false., agglomeration=ls_agglomeration, consolidation=ls_consolidation, &
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
    err = multigrid%solve(sol, rhs, 1.e-10_amrex_real, 0.0_amrex_real)

    ! Clean memory
    call amrex_abeclaplacian_destroy(abeclap)
    call amrex_multigrid_destroy(multigrid)
    
  end subroutine get_temperature_implicit

    
  ! -----------------------------------------------------------------
  ! Subroutine used to fill the alpha matrix for the linear solver
  ! -----------------------------------------------------------------  
  subroutine get_alpha(lo, hi, &
                       u_old, uo_lo, uo_hi, &
                       u_star, us_lo, us_hi, &
                       temp, t_lo, t_hi, &
                       temp_star, ts_lo, ts_hi, &
                       alpha, a_lo, a_hi)

    use material_properties_module, only: get_mass_density, get_heat_capacity
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)  
    integer, intent(in) :: uo_lo(2), uo_hi(2)
    integer, intent(in) :: us_lo(2), us_hi(2)
    integer, intent(in) :: t_lo(2), t_hi(2)
    integer, intent(in) :: ts_lo(2), ts_hi(2)
    integer, intent(in) :: a_lo(2), a_hi(2)
    real(amrex_real), intent(in) :: u_old(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2))
    real(amrex_real), intent(in) :: u_star(us_lo(1):us_hi(1),us_lo(2):us_hi(2))    
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2))
    real(amrex_real), intent(in) :: temp_star(ts_lo(1):ts_hi(1),ts_lo(2):ts_hi(2))
    real(amrex_real), intent(out) :: alpha(a_lo(1):a_hi(1),a_lo(2):a_hi(2))

    ! Local variables
    integer :: i,j
    real(amrex_real) :: cp
    real(amrex_real) :: rho
       
    ! Fill alpha matrix
    do i = lo(1), hi(1)
       do j = lo(2), hi(2)          
          !alpha(i,j) = (u_star(i,j) - u_old(i,j))/(temp_star(i,j) - temp(i,j))
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
                     u_old, uo_lo, uo_hi, &
                     u_star, us_lo, us_hi, &
                     temp, t_lo, t_hi, &
                     temp_star, ts_lo, ts_hi, &
                     rhs, r_lo, r_hi)

    use material_properties_module, only: get_mass_density, get_heat_capacity
    use heat_flux_module, only: get_boundary_heat_flux
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: id_lo(2), id_hi(2)
    integer, intent(in) :: uo_lo(2), uo_hi(2)
    integer, intent(in) :: us_lo(2), us_hi(2)
    integer, intent(in) :: t_lo(2), t_hi(2)
    integer, intent(in) :: ts_lo(2), ts_hi(2)
    integer, intent(in) :: r_lo(2), r_hi(2)    
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2))
    real(amrex_real), intent(in) :: u_old(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2))
    real(amrex_real), intent(in) :: u_star(us_lo(1):us_hi(1),us_lo(2):us_hi(2))    
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2))
    real(amrex_real), intent(in) :: temp_star(ts_lo(1):ts_hi(1),ts_lo(2):ts_hi(2))
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
          ! rhs(i,j) = rhs(i,j)*dt + &
          !       temp(i,j)*(u_star(i,j) - u_old(i,j))/(temp_star(i,j) - temp(i,j))
          call get_mass_density(temp(i,j), rho)
          call get_heat_capacity(temp(i,j), cp)
          rhs(i,j) = rhs(i,j)*dt + temp(i,j)*rho*cp
       end do
    end do

  end subroutine get_rhs
 

  ! ! -----------------------------------------------------------------
  ! ! Subroutine used to get the corrected enthalpy
  ! ! -----------------------------------------------------------------  
  ! subroutine get_enthalpy_implicit(lo, hi, &
  !                                  u_old, uo_lo, uo_hi, &
  !                                  u_star, us_lo, us_hi, &
  !                                  u_new, un_lo, un_hi, &
  !                                  temp_old, to_lo, to_hi, &
  !                                  temp_star, ts_lo, ts_hi, &
  !                                  temp_new, tn_lo, tn_hi)
  				
  !   ! Input and output variables
  !   integer, intent(in) :: lo(2), hi(2)  
  !   integer, intent(in) :: uo_lo(2), uo_hi(2)
  !   integer, intent(in) :: us_lo(2), us_hi(2)
  !   integer, intent(in) :: un_lo(2), un_hi(2)
  !   integer, intent(in) :: to_lo(2), to_hi(2)
  !   integer, intent(in) :: ts_lo(2), ts_hi(2)
  !   integer, intent(in) :: tn_lo(2), tn_hi(2)
  !   real(amrex_real), intent(in) :: u_old(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2))
  !   real(amrex_real), intent(in) :: u_star(us_lo(1):us_hi(1),us_lo(2):us_hi(2))
  !   real(amrex_real), intent(out) :: u_new(un_lo(1):un_hi(1),un_lo(2):un_hi(2))
  !   real(amrex_real), intent(in) :: temp_old(to_lo(1):to_hi(1),to_lo(2):to_hi(2))
  !   real(amrex_real), intent(in) :: temp_star(ts_lo(1):ts_hi(1),ts_lo(2):ts_hi(2))
  !   real(amrex_real), intent(in) :: temp_new(tn_lo(1):tn_hi(1),tn_lo(2):tn_hi(2))

  !   ! Local variables
  !   integer :: i,j
    
  !   ! Fill alpha matrix
  !   do i = lo(1), hi(1)
  !      do j = lo(2), hi(2)          
  !         u_new(i,j) = u_old(i,j) + (temp_new(i,j) - temp_old(i,j)) * &
  !                      (u_star(i,j) - u_old(i,j))/(temp_star(i,j) - temp_old(i,j))
  !      end do
  !   end do

  ! end subroutine get_enthalpy_implicit


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
  
end module heat_transfer_module
