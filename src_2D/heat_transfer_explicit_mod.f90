module heat_transfer_explicit_module
  
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
  public :: advance_heat_solver_explicit_level
 
contains

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy at a new time step for
  ! a given level via an explicit update.
  ! -----------------------------------------------------------------
  subroutine advance_heat_solver_explicit_level(lev, time, dt, substep)
    
    ! Input and output variables 
    integer, intent(in) :: lev
    integer, intent(in) :: substep
    real(amrex_real), intent(in) :: time
    real(amrex_real), intent(in) :: dt
    
    ! Local variables
    type(amrex_multifab) :: phi_tmp ! Enthalpy multifab with ghost points
    type(amrex_multifab) :: temp_tmp ! Temperature multifab with ghost points
    type(amrex_multifab) :: phi_tmp2 ! Enthalpy multifab with 2 ghost points
    type(amrex_multifab) :: temp_tmp2 ! Temperature multifab with 2 ghost points
    type(amrex_multifab) :: idomain_tmp ! Idomain multifab to distinguish material and background
    type(amrex_multifab) :: fluxes(amrex_spacedim)

    ! Initialize temporary multifabs and fluxes
    call init_tmp_multifab(lev, time, phi_tmp, temp_tmp, phi_tmp2, temp_tmp2, &
                           idomain_tmp, fluxes)

    ! Advance heat transfer equation of one timestep
    call advance(lev, time, dt, substep, phi_tmp, temp_tmp, &
                 idomain_tmp, fluxes, phi_tmp2, temp_tmp2)

    ! Update flux registers (fluxes have already been scaled by dt and area
    ! in the advance_heat_solver_box subroutine)
    call update_flux_registers(lev, fluxes)
    
    ! Clean memory
    call amrex_multifab_destroy(phi_tmp)
    call amrex_multifab_destroy(temp_tmp)
    call amrex_multifab_destroy(phi_tmp2)
    call amrex_multifab_destroy(temp_tmp2)
    call amrex_multifab_destroy(idomain_tmp)
    
  end subroutine advance_heat_solver_explicit_level


  ! -----------------------------------------------------------------
  ! Subroutine used to initialize the multifabs and fluxes used in the
  ! update of the heat equation
  ! -----------------------------------------------------------------
  subroutine init_tmp_multifab(lev, time, phi_tmp, temp_tmp, phi_tmp2, &
                               temp_tmp2, idomain_tmp, fluxes)

    use read_input_module, only : do_reflux
    use amr_data_module, only : phi_new, &
                                phi_old, &
                                temp, &
                                idomain 
    use regrid_module, only : fillpatch
    
    ! Input and output variables 
    integer, intent(in) :: lev
    real(amrex_real), intent(in) :: time
    type(amrex_multifab), intent(out) :: phi_tmp
    type(amrex_multifab), intent(out) :: temp_tmp
    type(amrex_multifab), intent(out) :: phi_tmp2
    type(amrex_multifab), intent(out) :: temp_tmp2
    type(amrex_multifab), intent(out) :: idomain_tmp
    type(amrex_multifab), intent(out) :: fluxes(amrex_spacedim)
    
    ! Local variables
    integer, parameter :: nghost = 1 ! number of ghost points in each spatial direction 
    integer :: ncomp
    integer :: srccmp
    integer :: dstcmp
    integer :: idim
    logical :: nodal(2) 
    type(amrex_geometry) :: geom
    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm
    
    ! Get geometry
    geom = amrex_geom(lev)
    
    ! Get number of components
    ncomp = phi_new(lev)%ncomp()
    
    ! Get boxarray and distribution mapping 
    ba = phi_new(lev)%ba
    dm = phi_new(lev)%dm
    
    ! Swap old and new enthalpies before updating solution
    call amrex_multifab_swap(phi_old(lev), phi_new(lev))
    
    ! Build enthalpy and temperature multifabs with ghost points
    call amrex_multifab_build(phi_tmp, ba, dm, ncomp, nghost) 
    call amrex_multifab_build(temp_tmp, ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_tmp2, ba, dm, ncomp, nghost+1) 
    call amrex_multifab_build(temp_tmp2, ba, dm, ncomp, nghost+1)
    call fillpatch(lev, time, phi_tmp)

    ! Build temporary idomain multifab to store the domain configuration before SW deformation
    call amrex_multifab_build(idomain_tmp, ba, dm, ncomp, nghost)
    call amrex_multifab_swap(idomain_tmp, idomain(lev))

    ! Fill in multifabs with two ghost points
    srccmp = 1
    dstcmp = 1
    call temp_tmp2%copy(temp(lev), srccmp, dstcmp, ncomp, nghost+1)
    call temp_tmp2%fill_boundary(geom)
    call phi_tmp2%copy(phi_old(lev), srccmp, dstcmp, ncomp, nghost+1)
    call phi_tmp2%fill_boundary(geom)

    ! Initialize fluxes  
    if (do_reflux) then
       do idim = 1, amrex_spacedim
          nodal = .false.
          nodal(idim) = .true.
          call amrex_multifab_build(fluxes(idim), ba, dm, ncomp, 0, nodal)
       end do       
    end if

    ! Clean memory
    call amrex_distromap_destroy(dm)
    
  end subroutine init_tmp_multifab

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy at a new time step for
  ! a given level via an explicit update. No flux register
  ! synchronization is performed
  ! -----------------------------------------------------------------
  subroutine advance(lev, time, dt, substep, phi_tmp, temp_tmp, &
                     idomain_tmp, fluxes, phi_tmp2, temp_tmp2)

    use amr_data_module, only : phi_new
    
    ! Input and output variables 
    integer, intent(in) :: lev
    integer, intent(in) :: substep
    real(amrex_real), intent(in) :: time
    real(amrex_real), intent(in) :: dt
    type(amrex_multifab), intent(inout) :: phi_tmp
    type(amrex_multifab), intent(inout) :: temp_tmp
    type(amrex_multifab), intent(inout) :: phi_tmp2
    type(amrex_multifab), intent(inout) :: temp_tmp2
    type(amrex_multifab), intent(inout) :: idomain_tmp
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
    
    !$omp parallel private(mfi, flux)

    do idim = 1, amrex_spacedim
      call flux(idim)%reset_omp_private()
    end do

    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.false.)  
    do while(mfi%next())
       call advance_box(lev, time, dt, substep, mfi, &
                        geom, ncomp, phi_tmp, temp_tmp, &
                        idomain_tmp, flux, fluxes, phi_tmp2, temp_tmp2)
    end do
    call amrex_mfiter_destroy(mfi)
    
    !$omp end parallel

    ! Clean memory
    do idim = 1, amrex_spacedim
      call amrex_fab_destroy(flux(idim))
   end do
   
  end subroutine advance

  ! -----------------------------------------------------------------
  ! Subroutine used to update the flux registers after the
  ! heat transfer equation has been advanced in time
  ! -----------------------------------------------------------------
  subroutine update_flux_registers(lev, fluxes)

    use read_input_module, only : do_reflux
    use amr_data_module, only : flux_reg
    
    ! Input and output variables 
    integer, intent(in) :: lev
    type(amrex_multifab), intent(inout) :: fluxes(amrex_spacedim)

    ! Local variables
    integer :: idim
    
    ! Update flux registers (fluxes have already been scaled by dt and area
    ! in the advance subroutine)
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
    
  end subroutine update_flux_registers
  
  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy at a new time step for
  ! a given box on a given level via an explicit update
  ! -----------------------------------------------------------------
  subroutine advance_box(lev, time, dt, substep, mfi, &
                                              geom, ncomp, phi_tmp, temp_tmp, &
                                              idomain_tmp, flux, fluxes, &
                                              phi_tmp2, temp_tmp2)

    use amr_data_module, only : phi_new, &
                                temp, &
                                idomain
    use read_input_module, only : do_reflux, &
                                  temp_fs
    use heat_transfer_domain_module, only : get_idomain, &
                                            get_melt_pos, &
                                            revaluate_heat_domain
    use material_properties_module, only : get_temp
    
    ! Input and output variables
    integer, intent(in) :: lev
    integer, intent(in) :: ncomp
    integer, intent(in) :: substep
    real(amrex_real), intent(in) :: dt
    real(amrex_real), intent(in) :: time
    type(amrex_mfiter), intent(in) :: mfi
    type(amrex_geometry), intent(in) :: geom
    type(amrex_multifab), intent(inout) :: phi_tmp
    type(amrex_multifab), intent(inout) :: temp_tmp
    type(amrex_multifab), intent(in) :: phi_tmp2
    type(amrex_multifab), intent(in) :: temp_tmp2
    type(amrex_multifab), intent(inout) :: idomain_tmp
    type(amrex_fab), intent(inout) :: flux(amrex_spacedim)
    type(amrex_multifab), intent(inout) :: fluxes(amrex_spacedim)

    ! Local variables
    integer :: idim
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pin
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pin2
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pout
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptempin
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptempin2
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
    pin2     => phi_tmp2%dataptr(mfi)
    pout    => phi_new(lev)%dataptr(mfi)
    ptempin => temp_tmp%dataptr(mfi)
    ptempin2 => temp_tmp2%dataptr(mfi)
    ptemp   => temp(lev)%dataptr(mfi)
    pidin   => idomain_tmp%dataptr(mfi)
    pidout  => idomain(lev)%dataptr(mfi)   
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
    call revaluate_heat_domain(geom%get_physical_location(bx%lo), geom%dx, &
                               bx%lo, bx%hi, &
                               pidin, lbound(pidin), ubound(pidin), &
                               pidout, lbound(pidout), ubound(pidout), &
                               pin, lbound(pin), ubound(pin), &
                               ptempin, lbound(ptempin), ubound(ptempin), &
                               pin2, lbound(pin2), ubound(pin2), &
                               ptempin2, lbound(ptempin2), ubound(ptempin2))
    
    ! Increment enthalpy at given box depending on the condition of the free surface
    if (temp_fs.gt.0) then
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
          flxx(i,j) = flxx(i,j) * (dt * dx(2))
       end do
    end do

    do j = lo(2), hi(2) + 1
       do i = lo(1), hi(1)
          flxy(i,j) = flxy(i,j) * (dt * dx(1))
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
  
end module heat_transfer_explicit_module
