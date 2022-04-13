module heat_transfer_explicit_no_subcycling
  
  ! -----------------------------------------------------------------
  ! This module contains the explicit solver for the heat transfer
  ! equations
  ! -----------------------------------------------------------------
  
  use amrex_amr_module
  
  implicit none

  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: advance_heat_solver_explicit
 
contains

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy at a new time step for
  ! a given level via an explicit update.
  ! -----------------------------------------------------------------
  subroutine advance_heat_solver_explicit(time, dt)

    use heat_transfer_domain_module, only : reset_melt_pos

    ! Input and output variables
    real(amrex_real), intent(in) :: dt
    real(amrex_real), intent(in) :: time
    
    ! Local variables
    type(amrex_multifab) :: phi_tmp(0:amrex_max_level) ! Enthalpy multifab with ghost points
    type(amrex_multifab) :: temp_tmp(0:amrex_max_level) ! Temperature multifab with ghost points
    type(amrex_multifab) :: idomain_tmp(0:amrex_max_level) ! Idomain multifab to distinguish material and background
    type(amrex_multifab) :: fluxes(0:amrex_max_level, amrex_spacedim)

    ! Initialize temporary multifabs and fluxes
    call init_tmp_multifab(time, phi_tmp, temp_tmp, idomain_tmp, fluxes)

    ! Reset melt position (to properly treat resolidification)
    call reset_melt_pos
    
    ! Advance heat transfer equation of one timestep
    call advance(time, dt, phi_tmp, temp_tmp, idomain_tmp, fluxes)

    ! Update flux registers 
    call update_flux_registers(fluxes)

    ! Synchronize idomains
    call synch_idomain(temp_tmp)
    
    ! Update melt position
    call update_melt_pos(amrex_max_level)
      
    ! Clean memory
    call free_memory(phi_tmp, temp_tmp, idomain_tmp, fluxes)
    
  end subroutine advance_heat_solver_explicit


  ! -----------------------------------------------------------------
  ! Subroutine used to initialize the multifabs and fluxes used in the
  ! update of the heat equation
  ! -----------------------------------------------------------------
  subroutine init_tmp_multifab(time, phi_tmp, temp_tmp, idomain_tmp, &
                               fluxes)

    use read_input_module, only : heat_reflux
    use amr_data_module, only : phi_new, &
                                phi_old, &
                                idomain 
    use regrid_module, only : fillpatch
    
    ! Input and output variables
    real(amrex_real), intent(in) :: time
    type(amrex_multifab), intent(out) :: phi_tmp(0:amrex_max_level)
    type(amrex_multifab), intent(out) :: temp_tmp(0:amrex_max_level)
    type(amrex_multifab), intent(out) :: idomain_tmp(0:amrex_max_level)
    type(amrex_multifab), intent(out) :: fluxes(0:amrex_max_level, amrex_spacedim)
    
    ! Local variables
    integer, parameter :: nghost = 1 ! number of ghost points in each spatial direction 
    integer :: ncomp
    integer :: ilev
    integer :: idim
    logical :: nodal(2) 
    type(amrex_geometry) :: geom
    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm

    do ilev = 0, amrex_max_level
    
       ! Get geometry
       geom = amrex_geom(ilev)
       
       ! Get number of components
       ncomp = phi_new(ilev)%ncomp()
       
       ! Get boxarray and distribution mapping 
       ba = phi_new(ilev)%ba
       dm = phi_new(ilev)%dm
       
       ! Swap old and new enthalpies before updating solution
       call amrex_multifab_swap(phi_old(ilev), phi_new(ilev))
       
       ! Build enthalpy and temperature multifabs with ghost points
       call amrex_multifab_build(phi_tmp(ilev), ba, dm, ncomp, nghost) 
       call amrex_multifab_build(temp_tmp(ilev), ba, dm, ncomp, nghost)
       call fillpatch(ilev, time, phi_tmp(ilev))
       
       ! Build temporary idomain multifab to store the domain configuration before SW deformation
       call amrex_multifab_build(idomain_tmp(ilev), ba, dm, ncomp, nghost)
       call amrex_multifab_swap(idomain_tmp(ilev), idomain(ilev))
       
       ! Initialize fluxes  
       if (heat_reflux) then
          do idim = 1, amrex_spacedim
             nodal = .false.
             nodal(idim) = .true.
             call amrex_multifab_build(fluxes(ilev, idim), ba, dm, ncomp, 0, nodal)
          end do
       end if
       
       ! Clean memory
       call amrex_distromap_destroy(dm)
       
    end do
    
  end subroutine init_tmp_multifab

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy at a new time step for
  ! a given level via an explicit update. No flux register
  ! synchronization is performed
  ! -----------------------------------------------------------------
  subroutine advance(time, dt, phi_tmp, temp_tmp, idomain_tmp, fluxes)

    use amr_data_module, only : phi_new
    
    ! Input and output variables
    real(amrex_real), intent(in) :: time
    real(amrex_real), intent(in) :: dt
    type(amrex_multifab), intent(inout) :: phi_tmp(0:amrex_max_level)
    type(amrex_multifab), intent(inout) :: temp_tmp(0:amrex_max_level)
    type(amrex_multifab), intent(inout) :: idomain_tmp(0:amrex_max_level)
    type(amrex_multifab), intent(inout) :: fluxes(0:amrex_max_level, amrex_spacedim)
    
    ! Local variables
    integer :: ncomp
    integer :: idim
    integer :: ilev
    type(amrex_fab) :: flux(amrex_spacedim)
    type(amrex_geometry) :: geom
    type(amrex_mfiter) :: mfi
    
    do ilev = 0, amrex_max_level
       
       ! Get geometry
       geom = amrex_geom(ilev)
       
       ! Get number of components
       ncomp = phi_new(ilev)%ncomp()
       
       !$omp parallel private(mfi, flux)
       
       do idim = 1, amrex_spacedim
          call flux(idim)%reset_omp_private()
       end do
       
       call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.false.)  
       do while(mfi%next())
          call advance_box(ilev, time, dt, mfi, &
                           geom, ncomp, phi_tmp(ilev), &
                           temp_tmp(ilev), idomain_tmp(ilev), &
                           flux, fluxes(ilev,:))
       end do
       call amrex_mfiter_destroy(mfi)
       
       !$omp end parallel
       
       ! Clean memory
       do idim = 1, amrex_spacedim
          call amrex_fab_destroy(flux(idim))
       end do

    end do
       
  end subroutine advance

  ! -----------------------------------------------------------------
  ! Subroutine used to update the flux registers after the
  ! heat transfer equation has been advanced in time
  ! -----------------------------------------------------------------
  subroutine update_flux_registers(fluxes)

    use read_input_module, only : heat_reflux
    use amr_data_module, only : flux_reg
    
    ! Input and output variables 
    type(amrex_multifab), intent(inout) :: fluxes(0:amrex_max_level, amrex_spacedim)

    ! Local variables
    integer :: ilev
    
    ! Update flux registers (fluxes have already been scaled by dt and area
    ! in the advance subroutine)
    if (heat_reflux) then

       do ilev = 0, amrex_max_level

          if (ilev > 0) then
             call flux_reg(ilev)%fineadd(fluxes(ilev,:), 1.0_amrex_real)
          end if
          
          if (ilev < amrex_get_finest_level()) then
             call flux_reg(ilev+1)%crseinit(fluxes(ilev,:), -1.0_amrex_real)
          end if

       end do
       
    end if
    
  end subroutine update_flux_registers

  ! -----------------------------------------------------------------
  ! Subroutine used to synchronize the idomains on all levels
  ! -----------------------------------------------------------------
  subroutine synch_idomain(temp_tmp)
    
    use amr_data_module, only : phi_new, &
                                idomain, &
                                temp    
    use heat_transfer_domain_module, only : get_idomain
        
    ! Input and output variables
    type(amrex_multifab), intent(inout) :: temp_tmp(0:amrex_max_level)
    
    ! Local variables
    integer :: ilev
    integer :: ncomp
    type(amrex_geometry) :: geom
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidom
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptemp_tmp
    
    do ilev = 0, amrex_max_level

       ! Geometry
       geom = amrex_geom(ilev)

       ! Get number of components
       ncomp = phi_new(ilev)%ncomp()
       
       ! Synchronize temperature with ghost points
       call temp_tmp(ilev)%copy(temp(ilev), 1, 1, ncomp, 1) ! The last 1 is the number of ghost points
       call temp_tmp(ilev)%fill_boundary(geom)
       
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
    
  end subroutine synch_idomain

  ! -----------------------------------------------------------------
  ! Subroutine used to free memory on all levels
  ! -----------------------------------------------------------------
  subroutine free_memory(phi_tmp, temp_tmp, idomain_tmp, fluxes)

    use read_input_module, only: heat_reflux
    
    ! Input and output variables
    type(amrex_multifab), intent(inout) :: phi_tmp(0:amrex_max_level)
    type(amrex_multifab), intent(inout) :: temp_tmp(0:amrex_max_level)
    type(amrex_multifab), intent(inout) :: idomain_tmp(0:amrex_max_level)
    type(amrex_multifab), intent(inout) :: fluxes(0:amrex_max_level, amrex_spacedim)
    
    ! Local variables
    integer :: idim
    integer :: ilev

    do ilev = 0, amrex_max_level
       
       call amrex_multifab_destroy(phi_tmp(ilev))
       call amrex_multifab_destroy(temp_tmp(ilev))
       call amrex_multifab_destroy(idomain_tmp(ilev))
       if (heat_reflux) then
          do idim = 1, amrex_spacedim
             call amrex_multifab_destroy(fluxes(ilev, idim))
          end do
       end if
       
    end do
       
  end subroutine free_memory
  
  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy at a new time step for
  ! a given box on a given level via an explicit update
  ! -----------------------------------------------------------------
  subroutine advance_box(lev, time, dt, mfi, &
                         geom, ncomp, phi_tmp, temp_tmp, &
                         idomain_tmp, flux, fluxes)
    
    use amr_data_module, only : phi_new, &
                                temp, &
                                idomain
    use read_input_module, only : heat_reflux, &
                                  heat_temp_surf
    use heat_transfer_domain_module, only : get_idomain, &
                                            get_melt_pos, &
                                            revaluate_heat_domain
    use material_properties_module, only : get_temp
    
    ! Input and output variables
    integer, intent(in) :: lev
    integer, intent(in) :: ncomp
    real(amrex_real), intent(in) :: dt
    real(amrex_real), intent(in) :: time
    type(amrex_mfiter), intent(in) :: mfi
    type(amrex_geometry), intent(in) :: geom
    type(amrex_multifab), intent(inout) :: phi_tmp
    type(amrex_multifab), intent(inout) :: temp_tmp
    type(amrex_multifab), intent(inout) :: idomain_tmp
    type(amrex_fab), intent(inout) :: flux(amrex_spacedim)
    type(amrex_multifab), intent(inout) :: fluxes(amrex_spacedim)

    ! Local variables
    integer :: idim
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
    type(amrex_box) :: bx
    type(amrex_box) :: tbx
    
    ! Box
    bx = mfi%validbox()   
    
    ! Pointers
    pin     => phi_tmp%dataptr(mfi)
    pout    => phi_new(lev)%dataptr(mfi)
    ptempin => temp_tmp%dataptr(mfi)
    ptemp   => temp(lev)%dataptr(mfi)
    pidin   => idomain_tmp%dataptr(mfi)
    pidout  => idomain(lev)%dataptr(mfi)   
    do idim = 1, amrex_spacedim
       tbx = bx
       call tbx%nodalize(idim)
       call flux(idim)%resize(tbx,ncomp)
       call tbx%grow(lev+1)
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
    call revaluate_heat_domain(geom%get_physical_location(bx%lo), geom%dx, &
                               bx%lo, bx%hi, &
                               pidin, lbound(pidin), ubound(pidin), &
                               pidout, lbound(pidout), ubound(pidout), &
                               pin, lbound(pin), ubound(pin), &
                               ptempin, lbound(ptempin), ubound(ptempin))
    
    ! Increment enthalpy at given box depending on the condition of the free surface
    if (heat_temp_surf.gt.0) then
       call get_enthalpy_fixT(bx%lo, bx%hi, &
                              pin, lbound(pin), ubound(pin),     &
                              pout,    lbound(pout), ubound(pout),    &
                              ptempin, lbound(ptempin), ubound(ptempin), &
                              pfx, lbound(pfx), ubound(pfx), &
                              pfy, lbound(pfy), ubound(pfy), &
                              pidout, lbound(pidout), ubound(pidout), &
                              geom, dt)
    else
       call get_enthalpy(time, bx%lo, bx%hi, &
                         pin, lbound(pin),     ubound(pin),     &
                         pout,    lbound(pout),    ubound(pout),    &
                         ptempin, lbound(ptempin), ubound(ptempin), &
                         pfx, lbound(pfx), ubound(pfx), &
                         pfy, lbound(pfy), ubound(pfy), &
                         pidout, lbound(pidout), ubound(pidout), &
                         geom, dt, lev)
    end if

    ! Get temperature corresponding to the new enthalpy
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
    
  end subroutine advance_box


  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy at a new time step for
  ! a given box of a certain level via an explicit update
  ! -----------------------------------------------------------------
  subroutine get_enthalpy(time, lo, hi, &
                          u_old,  uo_lo, uo_hi, &
                          u_new, un_lo, un_hi, &
                          temp, t_lo, t_hi, &
                          flxx, fx_lo, fx_hi, &
                          flxy, fy_lo, fy_hi, &
                          idom, id_lo, id_hi, &
                          geom, dt, lev)
    
    use heat_transfer_flux_module, only: get_boundary_heat_flux
 
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
    integer, intent(in) :: lev
    
    ! Local variables
    integer :: i,j
    real(amrex_real) :: dx(2) ! Grid size
    real(amrex_real) :: lo_phys(2) ! Physical location of the lowest corner of the tile box
    real(amrex_real) :: qbound(lo(1):hi(1),lo(2):hi(2)) ! Volumetric heating (boundary)
    real(amrex_real) :: qvol(lo(1):hi(1),lo(2):hi(2)) ! Volumetric heating 
 
    ! Get grid size
    dx = geom%dx(1:2) ! grid width at level 

    ! Get physical location of the lowest corner of the tile box
    lo_phys = geom%get_physical_location(lo)
                  
    ! Get enthalpy flux 
    call create_face_flux(dx, lo_phys, lo, hi, &
                          u_old, uo_lo, uo_hi, &
                          flxx, fx_lo, fx_hi, &
                          flxy, fy_lo, fy_hi, &
                          temp, t_lo, t_hi, &
                          idom, id_lo, id_hi)
  				  	
    ! Prescribe external heat flux on the free surface
    call get_boundary_heat_flux(time, lo_phys, &
                                dx, lo, hi, &
                                idom, id_lo, id_hi, &
                                temp, t_lo, t_hi, lev, qbound)

    ! Volumetric sources
    call get_volumetric_heat_source(dx, lo_phys, lo, hi, &
                                    u_old, uo_lo, uo_hi, &
                                    idom, id_lo, id_hi, &
                                    qvol)
    
    ! Compute enthalpy at the new timestep
    do  j = lo(2),hi(2)
       do i = lo(1),hi(1)
          u_new(i,j) = u_old(i,j) &
                       - dt/dx(1) * (flxx(i+1,j) - flxx(i,j)) & ! flux divergence x-direction 
                       - dt/dx(2) * (flxy(i,j+1) - flxy(i,j)) & ! flux divergence y-direction 
                       + dt*qbound(i,j) & ! 'boundary volumetric' source
                       + dt*qvol(i,j)  ! volumetric source
       end do
    end do
    
  	
    ! Scale the fluxes for the flux registers
    do j = lo(2), hi(2)
       do i = lo(1), hi(1) + 1
          ! flxx(i,j) = flxx(i,j) * (dt * dx(2))
          flxx(i,j) = flxx(i,j) * (dx(2))
       end do
    end do

    do j = lo(2), hi(2) + 1
       do i = lo(1), hi(1)
          ! flxy(i,j) = flxy(i,j) * (dt * dx(1))
          flxy(i,j) = flxy(i,j) * (dx(1))
       end do
    end do
    
  end subroutine get_enthalpy



  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy at a new time step for
  ! a given box on a certain level via an explicit update.
  ! Case with fixed temperature at the free surface
  ! -----------------------------------------------------------------
  subroutine get_enthalpy_fixT(lo, hi, &
                                        u_old,  uo_lo, uo_hi, &
                                        u_new, un_lo, un_hi, &
                                        temp, t_lo, t_hi, &
                                        flxx, fx_lo, fx_hi, &
                                        flxy, fy_lo, fy_hi, &
                                        idom, id_lo, id_hi, &
                                        geom, dt)
    
    use material_properties_module, only : get_enthalpy
    use read_input_module, only : heat_temp_surf

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
    call get_enthalpy(heat_temp_surf, u_fs)
    
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
    do  j = lo(2),hi(2)
       do i = lo(1),hi(1)
           
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
    
  end subroutine get_enthalpy_fixT


  
  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy fluxes on the edges of
  ! the grid
  ! -----------------------------------------------------------------  
  subroutine create_face_flux(dx, lo_phys, lo, hi, &
                              u_old, uo_lo, uo_hi, &
                              flxx, fx_lo, fx_hi, &
                              flxy, fy_lo, fy_hi, &
                              temp, t_lo, t_hi, &
                              idom, id_lo, id_hi)
  				
    use material_properties_module, only: get_conductivity
    use heat_transfer_domain_module, only: get_face_velocity
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)  
    integer, intent(in) :: uo_lo(2), uo_hi(2)
    integer, intent(in) :: t_lo(2), t_hi(2)
    integer, intent(in) :: fx_lo(2), fx_hi(2)
    integer, intent(in) :: fy_lo(2), fy_hi(2)
    integer, intent(in) :: id_lo(2), id_hi(2)
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(in) :: lo_phys(2)
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
    call get_face_velocity(lo_phys, lo, hi, dx, &
                           vx, fx_lo, fx_hi, &
                           idom, id_lo, id_hi)
    
    ! Flux along the x direction
    do j = lo(2), hi(2)          
       do i = lo(1), hi(1)+1

          if (nint(idom(i-1,j)).eq.0 .or. nint(idom(i,j)).eq.0) then

             ! Suppress flux at the free surface
             flxx(i,j) = 0.0_amrex_real

          else if (nint(idom(i-1,j)).eq.-1 .or. nint(idom(i,j)).eq.-1) then

             ! Suppress flux at the surface of the cooling pipe
             flxx(i,j) = 0.0_amrex_real
             
          else
             
             ! Advective component
             if (vx(i,j).gt.0.0_amrex_real) then
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
    do j = lo(2), hi(2)+1
       do i = lo(1), hi(1)

          if (nint(idom(i,j-1)).eq.0 .or. nint(idom(i,j)).eq.0) then

             ! Suppress flux at the free surface
             flxy(i,j) = 0.0_amrex_real

          else if (nint(idom(i,j-1)).eq.-1 .or. nint(idom(i,j)).eq.-1) then

             ! Suppress flux at the surface of the cooling pipe
             flxy(i,j) = 0.0_amrex_real
             
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
  ! the grid. Case with fixed temperature at the free surface without
  ! convection
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
    do j = lo(2), hi(2)          
       do i = lo(1), hi(1)+1

          ! Diffusive component
          temp_face = (temp(i,j) + temp(i-1,j))/2_amrex_real
          call get_conductivity(temp_face, ktherm)
          flxx(i,j) = -ktherm*(temp(i,j)-temp(i-1,j))/dx(1)

       end do
    end do
    
    ! Flux along the y direction
    do j = lo(2), hi(2)+1
       do i = lo(1), hi(1)
       
          ! Diffusive component (there is no advection in the y direction)
          temp_face = (temp(i,j) + temp(i,j-1))/2_amrex_real
          call get_conductivity(temp_face, ktherm)
          flxy(i,j) = -ktherm*(temp(i,j)-temp(i,j-1))/dx(2)
          
       end do
    end do
 
  end subroutine create_face_flux_fixT

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the volumetric enthalpy source terms
  ! -----------------------------------------------------------------  
  subroutine get_volumetric_heat_source(dx, lo_phys, lo, hi, &
                                        u_old, uo_lo, uo_hi, &
                                        idom, id_lo, id_hi, &
                                        qvol)
  				
    use material_properties_module, only: get_conductivity
    use heat_transfer_domain_module, only: get_face_velocity
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)  
    integer, intent(in) :: uo_lo(2), uo_hi(2)
    integer, intent(in) :: id_lo(2), id_hi(2)
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(in) :: lo_phys(2)
    real(amrex_real), intent(in) :: u_old(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2))
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2))
    real(amrex_real), intent(out) :: qvol(lo(1):hi(1),lo(2):hi(2))

    ! Local variables
    integer :: i,j
    integer :: fx_lo(2), fx_hi(2)
    real(amrex_real) :: vx(lo(1):hi(1)+1,lo(2):hi(2))


    ! Construct 3D melt velocity profile from the 2D shallow water solution
    fx_lo = lo
    fx_hi = hi
    fx_hi(1) = fx_hi(1) + 1
    call get_face_velocity(lo_phys, lo, hi, dx, &
                           vx, fx_lo, fx_hi, &
                           idom, id_lo, id_hi)
    
    ! Volumetric heat source terms
    do j = lo(2), hi(2)          
       do i = lo(1), hi(1)
          qvol(i,j) = u_old(i,j) * (vx(i+1,j) - vx(i,j))/dx(1)
       end do
    end do
 
  end subroutine get_volumetric_heat_source

  
  ! -----------------------------------------------------------------
  ! Subroutine used to update the position of the bottom of the
  ! melt pool
  ! -----------------------------------------------------------------
  subroutine update_melt_pos(lev)
    
    use amr_data_module, only : phi_new,&
                                idomain
    use heat_transfer_domain_module, only : get_melt_pos

    ! Input variables
    integer, intent(in) :: lev
    
    ! Local variables
    type(amrex_geometry) :: geom
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidom
    
    if (lev.eq.amrex_max_level) then
       
       ! Geometry
       geom = amrex_geom(lev)
       
       ! Loop through all the boxes in the level
       !$omp parallel private(mfi, bx, pidom)
       call amrex_mfiter_build(mfi, phi_new(lev), tiling=.false.)
       do while(mfi%next())
          bx = mfi%validbox()
          pidom  => idomain(lev)%dataptr(mfi)
          call get_melt_pos(bx%lo, bx%hi, &
                             pidom, lbound(pidom), ubound(pidom), &
                             geom)
          
       end do
       call amrex_mfiter_destroy(mfi)
       !$omp end parallel   

    end if
    
  end subroutine update_melt_pos

  
end module heat_transfer_explicit_no_subcycling