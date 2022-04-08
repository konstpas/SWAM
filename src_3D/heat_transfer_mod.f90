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
 
     use read_input_module, only : do_reflux
     use amr_data_module, only : phi_new, phi_old, temp, idomain, flux_reg
     use regrid_module, only : fillpatch
     
     ! Input and output variables 
     integer, intent(in) :: lev
     integer, intent(in) :: substep
     real(amrex_real), intent(in) :: time
     real(amrex_real), intent(in) :: dt
     
     ! Local variables
     integer, parameter :: nghost = 1 ! number of ghost points in each spatial direction 
     integer :: ncomp, srccmp, dstcmp
     integer :: idim
     logical :: nodal(3) ! logical for flux multifabs 
     type(amrex_multifab) :: phi_tmp ! Enthalpy multifab with ghost points
     type(amrex_multifab) :: temp_tmp ! Temperature multifab with ghost points
     type(amrex_multifab) :: phi_tmp2 ! Enthalpy multifab with 2 ghost points
     type(amrex_multifab) :: temp_tmp2 ! Temperature multifab with 2 ghost points
     type(amrex_multifab) :: idomain_tmp ! Idomain multifab to distinguish material and background
     type(amrex_fab) :: flux(amrex_spacedim)
     type(amrex_multifab) :: fluxes(amrex_spacedim)
     type(amrex_geometry) :: geom
     type(amrex_mfiter) :: mfi
     
     ! Get geometry
     geom = amrex_geom(lev)
     
     ! Get number of components
     ncomp = phi_new(lev)%ncomp()
 
     ! Swap old and new enthalpies before updating solution
     call amrex_multifab_swap(phi_old(lev), phi_new(lev))
     
     ! Build enthalpy and temperature multifabs with ghost points
     call amrex_multifab_build(phi_tmp, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, nghost) 
     call amrex_multifab_build(temp_tmp, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, nghost)
     call amrex_multifab_build(phi_tmp2, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, nghost+1) 
     call amrex_multifab_build(temp_tmp2, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, nghost+1)
     call fillpatch(lev, time, phi_tmp)
 
     ! Build temporary idomain multifab to store the domain configuration before SW deformation
     call amrex_multifab_build(idomain_tmp, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, nghost)
     call amrex_multifab_swap(idomain_tmp, idomain(lev))

     ! Fill in multifabs with two ghose points
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
           call amrex_multifab_build(fluxes(idim), phi_new(lev)%ba, phi_new(lev)%dm, ncomp, 0, nodal)
        end do       
     end if

    !$OMP parallel private(mfi, flux)
     do idim = 1, amrex_spacedim
        call flux(idim)%reset_omp_private()
     end do
        
     ! Advance heat solver on all boxes on the level
     call amrex_mfiter_build(mfi, phi_new(lev), tiling=.false.)
     do while(mfi%next())
        call advance_heat_solver_explicit_box(lev, time, dt, substep, mfi, &
                                              geom, ncomp, phi_tmp, temp_tmp, &
                                              idomain_tmp, flux, fluxes, phi_tmp2, temp_tmp2)
     end do
     call amrex_mfiter_destroy(mfi)
    !$OMP end parallel

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
   ! a given box on a given level via an explicit update
   ! -----------------------------------------------------------------
   subroutine advance_heat_solver_explicit_box(lev, time, dt, substep, mfi, &
                                               geom, ncomp, phi_tmp, temp_tmp, &
                                               idomain_tmp, flux, fluxes, phi_tmp2, temp_tmp2)
 
     use amr_data_module, only : phi_new, &
                                 temp, &
                                 idomain, &
                                 Qplasma, &
                                 Qpipe, &
                                 Qtherm, &
                                 Qvap, &
                                 Qrad
     use read_input_module, only : do_reflux, temp_fs
     use domain_module, only : get_idomain, get_melt_pos, revaluate_heat_domain
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
     real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pfz
     real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pf
     real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pfab
     real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidin
     real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidout
     type(amrex_box) :: bx
     type(amrex_box) :: tbx
     real(amrex_real) :: Qpipe_box, Qplasma_box, Qtherm_box, Qrad_box, Qvap_box
     
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
     pfz => flux(3)%dataptr()
     
     
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
     call revaluate_heat_domain(lev, geom%get_physical_location(bx%lo), geom%dx, &
                                bx%lo, bx%hi, &
                                pidin, lbound(pidin), ubound(pidin), &
                                pidout, lbound(pidout), ubound(pidout), &
                                pin, lbound(pin), ubound(pin), &
                                ptempin, lbound(ptempin), ubound(ptempin), &
                                pin2, lbound(pin2), ubound(pin2), &
                                ptempin2, lbound(ptempin2), ubound(ptempin2))
     
     
     ! Increment enthalpy at given box depending on the condition of the free surface
     if (temp_fs.gt.0) then
        call get_enthalpy_explicit_fixT(bx%lo, bx%hi, &
                                        pin, lbound(pin),     ubound(pin),     &
                                        pout,    lbound(pout),    ubound(pout),    &
                                        ptempin, lbound(ptempin), ubound(ptempin), &
                                        pfx, lbound(pfx), ubound(pfx), &
                                        pfy, lbound(pfy), ubound(pfy), &
                                        pfz, lbound(pfz), ubound(pfz), &
                                        pidout, lbound(pidout), ubound(pidout), &
                                        geom, dt)
     else
        call get_enthalpy_explicit(time, bx%lo, bx%hi, &
                                   pin, lbound(pin),     ubound(pin),     &
                                   pout,    lbound(pout),    ubound(pout),    &
                                   ptempin, lbound(ptempin), ubound(ptempin), &
                                   pfx, lbound(pfx), ubound(pfx), &
                                   pfy, lbound(pfy), ubound(pfy), &
                                   pfz, lbound(pfz), ubound(pfz), &
                                   pidout, lbound(pidout), ubound(pidout), &
                                   geom, lev, dt, &
                                   Qpipe_box, Qtherm_box, Qvap_box, Qrad_box, Qplasma_box)
     end if
     !$omp critical
     Qplasma = Qplasma + Qplasma_box
     Qpipe = Qpipe + Qpipe_box
     Qtherm = Qtherm + Qtherm_box
     Qvap = Qvap + Qvap_box
     Qrad = Qrad + Qrad_box
     !$omp end critical
 
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
     
   end subroutine advance_heat_solver_explicit_box
   
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
                                    flxz, fz_lo, fz_hi, &
                                    idom, id_lo, id_hi, &
                                    geom, lev, dt, &
                                    Qpipe_box, Qtherm_box, Qvap_box, Qrad_box, Qplasma_box)

      use heat_flux_module, only: get_boundary_heat_flux

      ! Input and output variables
      integer, intent(in) :: lo(3), hi(3) ! bounds of current tile box
      integer, intent(in) :: uo_lo(3), uo_hi(3) ! bounds of input enthalpy box 
      integer, intent(in) :: un_lo(3), un_hi(3) ! bounds of output enthalpy box   
      integer, intent(in) :: t_lo (3), t_hi (3) ! bounds of the temperature box
      integer, intent(in) :: fx_lo(3), fx_hi(3) ! bounds of the enthalpy flux along x
      integer, intent(in) :: fy_lo(3), fy_hi(3) ! bounds of the enthalpy flux along y
      integer, intent(in) :: fz_lo(3), fz_hi(3) ! bounds of the enthalpy flux along z
      integer, intent(in) :: id_lo(3), id_hi(3) ! bounds of the output idomain box
      real(amrex_real), intent(in) :: dt ! time step
      real(amrex_real), intent(in) :: time ! time
      integer, intent(in) :: lev
      real(amrex_real), intent(inout) :: u_old (uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3)) ! Input enthalpy 
      real(amrex_real), intent(inout) :: u_new(un_lo(1):un_hi(1),un_lo(2):un_hi(2),un_lo(3):un_hi(3)) ! Output enthalpy
      real(amrex_real), intent(inout) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) ! Temperature
      real(amrex_real), intent(out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3)) ! flux along the x direction  			
      real(amrex_real), intent(out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3)) ! flux along the y direction
      real(amrex_real), intent(out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3)) ! flux along the z direction
      real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2),id_lo(3):id_hi(3))
      type(amrex_geometry), intent(in) :: geom ! geometry
      real(amrex_real), intent(out) :: Qplasma_box
      real(amrex_real), intent(out) :: Qpipe_box
      real(amrex_real), intent(out) :: Qtherm_box
      real(amrex_real), intent(out) :: Qrad_box
      real(amrex_real), intent(out) :: Qvap_box

      !Local variables
      integer :: i,j,k
      real(amrex_real) :: dx(3) ! Grid size
      real(amrex_real) :: lo_phys(3) ! Physical location of the lowest corner of the tile box
      real(amrex_real) :: qbound(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) ! Volumetric heating (boundary)

      ! Get grid size
      dx = geom%dx(1:3) ! grid width at level 

      ! Get physical location of the lowest corner of the tile box
      lo_phys = geom%get_physical_location(lo)

      ! Get enthalpy flux 
      call create_face_flux(dx, lo_phys, lo, hi, &
                            u_old, uo_lo, uo_hi, &
                            flxx, fx_lo, fx_hi, &
                            flxy, fy_lo, fy_hi, &
                            flxz, fz_lo, fz_hi, &
                            temp, t_lo, t_hi, &
                            idom, id_lo, id_hi)

      ! Prescribe external heat flux on the free surface
      call get_boundary_heat_flux(time, lo_phys, &
                                  dx, lo, hi, &
                                  idom, id_lo, id_hi, &
                                  temp, t_lo, t_hi, lev, dt, &
                                  Qpipe_box, Qtherm_box, Qvap_box, Qrad_box, &
                                  Qplasma_box, qbound)

      ! Compute enthalpy at the new timestep
      do   i = lo(1),hi(1)
         do  j = lo(2),hi(2) 
            do k = lo(3),hi(3)
               u_new(i,j,k) = u_old(i,j,k) &
               - dt/dx(1) * (flxx(i+1,j,k) - flxx(i,j,k)) &	! flux divergence x-direction 
               - dt/dx(2) * (flxy(i,j+1,k) - flxy(i,j,k)) &	! flux divergence y-direction 
               - dt/dx(3) * (flxz(i,j,k+1) - flxz(i,j,k)) &	! flux divergence z-direction
               + dt*qbound(i,j,k) ! 'boundary volumetric' source
            end do
         end do
      end do

      ! Scale the fluxes for the flux registers
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1) + 1
               flxx(i,j,k) = flxx(i,j,k) * (dt * dx(2)*dx(3))
            end do
         end do
      end do

      do k = lo(3), hi(3)
         do j = lo(2), hi(2) + 1
            do i = lo(1), hi(1)
               flxy(i,j,k) = flxy(i,j,k) * (dt * dx(1)*dx(3))
            end do
         end do
      end do

      do k = lo(3), hi(3) + 1
         do j = lo(2), hi(2)
               do i = lo(1), hi(1)
               flxz(i,j,k) = flxz(i,j,k) * (dt * dx(1)*dx(2))
               end do
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
                                         flxz, fz_lo, fz_hi, &
                                         idom, id_lo, id_hi, &
                                         geom, dt)

      use material_properties_module, only : get_enthalpy
      use read_input_module, only : temp_fs

      ! Input and output variables
      integer, intent(in) :: lo(3), hi(3) ! bounds of current tile box
      integer, intent(in) :: uo_lo(3), uo_hi(3) ! bounds of input enthalpy box 
      integer, intent(in) :: un_lo(3), un_hi(3) ! bounds of output enthalpy box  
      integer, intent(in) :: t_lo (3), t_hi (3) ! bounds of temperature box
      integer, intent(in) :: fx_lo(3), fx_hi(3) ! bounds of the enthalpy flux along x
      integer, intent(in) :: fy_lo(3), fy_hi(3) ! bounds of the enthalpy flux along y
      integer, intent(in) :: fz_lo(3), fz_hi(3) ! bounds of the enthalpy flux along z
      integer, intent(in) :: id_lo(3), id_hi(3) ! bounds of the output idomain box
      real(amrex_real), intent(in) :: dt ! time step
      real(amrex_real), intent(inout) :: u_old (uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3)) ! Input enthalpy 
      real(amrex_real), intent(inout) :: u_new(un_lo(1):un_hi(1),un_lo(2):un_hi(2),un_lo(3):un_hi(3)) ! Output enthalpy
      real(amrex_real), intent(inout) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) ! Temperature
      real(amrex_real), intent(out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3)) ! flux along the x direction  			
      real(amrex_real), intent(out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3)) ! flux along the y direction
      real(amrex_real), intent(out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3)) ! flux along the z direction
      real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2),id_lo(3):id_hi(3))
      type(amrex_geometry), intent(in) :: geom ! geometry

      !Local variables
      integer :: i,j,k
      real(amrex_real) :: u_fs ! Enthalpy of the points on the free surface
      real(amrex_real) :: dx(3) ! Grid size
      real(amrex_real) :: lo_phys(3) ! Physical location of the lowest corner of the tile box

      ! Enthalpy at the free surface
      call get_enthalpy(temp_fs, u_fs)

      ! Get grid size
      dx = geom%dx(1:3) ! grid width at level 

      ! Get physical location of the lowest corner of the tile box
      lo_phys = geom%get_physical_location(lo)

      ! Get enthalpy flux 
      call create_face_flux_fixT(dx, lo, hi, &
                                 flxx, fx_lo, fx_hi, &
                                 flxy, fy_lo, fy_hi, &
                                 flxz, fz_lo, fz_hi, &
                                 temp, t_lo, t_hi)

      ! Compute enthalpy at the new timestep
      do   i = lo(1),hi(1)
         do  j = lo(2),hi(2) 
            do k = lo(3),hi(3)

               if (nint(idom(i,j,k)).ne.0 .and. nint(idom(i,j+1,k)).eq.0) then
                  u_new(i,j,k) = u_fs ! Impose temperature on the free surface
               else if(nint(idom(i,j,k)).eq.0) then
                  u_new(i,j,k) = u_fs ! Set background temperature equal to the free surface temperature
               else
                  ! Update enthalpy according to heat equation
                  u_new(i,j,k) = u_old(i,j,k) &
                  - dt/dx(1) * (flxx(i+1,j,k) - flxx(i,j,k)) &	! flux divergence x-direction 
                  - dt/dx(2) * (flxy(i,j+1,k) - flxy(i,j,k)) &	! flux divergence y-direction 
                  - dt/dx(3) * (flxz(i,j,k+1) - flxz(i,j,k))       ! flux divergence z-direction
               end if

            end do
         end do
      end do

      ! Scale the fluxes for the flux registers
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1) + 1
               flxx(i,j,k) = flxx(i,j,k) * (dt * dx(2)*dx(3))
            end do
         end do
      end do

      do k = lo(3), hi(3)
         do j = lo(2), hi(2) + 1
            do i = lo(1), hi(1)
               flxy(i,j,k) = flxy(i,j,k) * (dt * dx(1)*dx(3))
            end do
         end do
      end do

      do k = lo(3), hi(3) + 1
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               flxz(i,j,k) = flxz(i,j,k) * (dt * dx(1)*dx(2))
            end do
         end do
      end do

   end subroutine get_enthalpy_explicit_fixT
 
 
  ! -----------------------------------------------------------------
  ! Subroutine used to the enthalpy fluxes on the edges of the grid
  ! -----------------------------------------------------------------  
  subroutine create_face_flux(dx, lo_phys, lo, hi, &
                              u_old, uo_lo, uo_hi, &
                              flxx, fx_lo, fx_hi, &
                              flxy, fy_lo, fy_hi, &
                              flxz, fz_lo, fz_hi, &
                              temp, t_lo, t_hi, &
                              idom, id_lo, id_hi)
  				
    use material_properties_module, only: get_conductivity

    ! Input and output variables
    integer, intent(in) :: lo(3), hi(3)  
    integer, intent(in) :: uo_lo(3), uo_hi(3)
    integer, intent(in) :: t_lo(3), t_hi(3)
    integer, intent(in) :: fx_lo(3), fx_hi(3)
    integer, intent(in) :: fy_lo(3), fy_hi(3)
    integer, intent(in) :: fz_lo(3), fz_hi(3)
    integer, intent(in) :: id_lo(3), id_hi(3)
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(in) :: lo_phys(3)
    real(amrex_real), intent(in) :: u_old(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3))
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2),id_lo(3):id_hi(3))
    real(amrex_real), intent(out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3))
    real(amrex_real), intent(out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3))
    real(amrex_real), intent(out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3))
    real(amrex_real), intent(in) :: temp (t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))

    ! Local variables
    integer :: i,j,k 
    real(amrex_real) :: ktherm
    real(amrex_real) :: temp_face
    real(amrex_real) :: vx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3))
    real(amrex_real) :: vz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3))

    ! Construct 3D melt velocity profile from the 2D shallow water solution
    call get_face_velocity(lo_phys, lo, hi, dx, &
                           vx, fx_lo, fx_hi, &
                           vz, fz_lo, fz_hi, &
                           idom, id_lo, id_hi)
    
    ! Flux along the x direction
    do i = lo(1), hi(1)+1
       do j = lo(2), hi(2)
          do k = lo(3), hi(3)

             if (nint(idom(i-1,j,k)).eq.0 .or. nint(idom(i,j,k)).eq.0 .or. &
              nint(idom(i-1,j,k)).eq.-1 .or. nint(idom(i,j,k)).eq.-1) then

                ! Suppress flux at the free surface and surface of cooling pipe
                flxx(i,j,k) = 0_amrex_real

             else
                
                ! Advective component
                if (nint(idom(i-1,j,k)).gt.2 .and. nint(idom(i,j,k)).gt.2) then
                   
                   if (vx(i,j,k).gt.0.0_amrex_real) then 
                      flxx(i,j,k)  = u_old(i-1,j,k)*vx(i,j,k)
                   else 
                      flxx(i,j,k)  = u_old(i,j,k)*vx(i,j,k)
                   end if

                else
                   flxx(i,j,k) = 0.0_amrex_real
                end if
                
                ! Diffusive component
                temp_face = (temp(i,j,k) + temp(i-1,j,k))/2_amrex_real
                call get_conductivity(temp_face, ktherm)
                flxx(i,j,k) = flxx(i,j,k) - ktherm*(temp(i,j,k)-temp(i-1,j,k))/dx(1)
                                
                end if
                
          end do
       end do
    end do

    ! Flux along the y direction
    do i = lo(1), hi(1)
       do j = lo(2), hi(2)+1
          do k = lo(3), hi(3)

             if (nint(idom(i,j-1,k)).eq.0 .or. nint(idom(i,j,k)).eq.0 .or. &
              nint(idom(i,j-1,k)).eq.-1 .or. nint(idom(i,j,k)).eq.-1) then

                ! Suppress flux at the free surface and surface of cooling pipe
                flxy(i,j,k) = 0_amrex_real

             else
                
                ! Diffusive component (there is no advection in the y direction)
                temp_face = (temp(i,j,k) + temp(i,j-1,k))/2_amrex_real
                call get_conductivity(temp_face, ktherm)
                flxy(i,j,k) = -ktherm*(temp(i,j,k)-temp(i,j-1,k))/dx(2)
                               
             end if
             
          end do
       end do
    end do

    ! Flux along the z direction
    do i = lo(1), hi(1)
       do j = lo(2), hi(2)
          do k = lo(3), hi(3)+1

             if (nint(idom(i,j,k-1)).eq.0 .or. nint(idom(i,j,k)).eq.0 .or. &
              nint(idom(i,j,k-1)).eq.-1 .or. nint(idom(i,j,k)).eq.-1) then

                ! Suppress flux at the free surface and surface of cooling pipe
                flxz(i,j,k) = 0_amrex_real

             else
                
                ! Advective component
                if (nint(idom(i,j,k-1)).gt.2 .and. nint(idom(i,j,k)).gt.2) then
                   
                   if (vz(i,j,k) > 0_amrex_real) then 
                      flxz(i,j,k)  = u_old(i,j,k-1)*vz(i,j,k)
                   else 
                      flxz(i,j,k)  = u_old(i,j,k)*vz(i,j,k)
                   end if

                else
                   flxz(i,j,k) = 0.0_amrex_real
                end if
                
                ! Diffusive component
                temp_face = (temp(i,j,k) + temp(i,j,k-1))/2_amrex_real
                call get_conductivity(temp_face, ktherm)
                flxz(i,j,k) = flxz(i,j,k) - ktherm*(temp(i,j,k)-temp(i,j,k-1))/dx(3)
                
             end if
             
          end do
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
                              flxz, fz_lo, fz_hi, &
                              temp, t_lo, t_hi)
  				
    use material_properties_module, only: get_conductivity

    ! Input and output variables
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: t_lo(3), t_hi(3)
    integer, intent(in) :: fx_lo(3), fx_hi(3)
    integer, intent(in) :: fy_lo(3), fy_hi(3)
    integer, intent(in) :: fz_lo(3), fz_hi(3)
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3))
    real(amrex_real), intent(out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3))
    real(amrex_real), intent(out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3))
    real(amrex_real), intent(in) :: temp (t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))

    ! Local variables
    integer :: i,j,k 
    real(amrex_real) :: ktherm
    real(amrex_real) :: temp_face
    
    ! Flux along the x direction
    do i = lo(1), hi(1)+1
       do j = lo(2), hi(2)
          do k = lo(3), hi(3)
                
             ! Diffusive component
             temp_face = (temp(i,j,k) + temp(i-1,j,k))/2_amrex_real
             call get_conductivity(temp_face, ktherm)
             flxx(i,j,k) = -ktherm*(temp(i,j,k)-temp(i-1,j,k))/dx(1)
                                
          end do
       end do
    end do

    ! Flux along the y direction
    do i = lo(1), hi(1)
       do j = lo(2), hi(2)+1
          do k = lo(3), hi(3)
             
             ! Diffusive component (there is no advection in the y direction)
             temp_face = (temp(i,j,k) + temp(i,j-1,k))/2_amrex_real
             call get_conductivity(temp_face, ktherm)
             flxy(i,j,k) = -ktherm*(temp(i,j,k)-temp(i,j-1,k))/dx(2)
             
          end do
       end do
    end do

    ! Flux along the z direction
    do i = lo(1), hi(1)
       do j = lo(2), hi(2)
          do k = lo(3), hi(3)+1
                
             ! Diffusive component
             temp_face = (temp(i,j,k) + temp(i,j,k-1))/2_amrex_real
             call get_conductivity(temp_face, ktherm)
             flxz(i,j,k) = -ktherm*(temp(i,j,k)-temp(i,j,k-1))/dx(3)
                
          end do
       end do
    end do
    
    
   end subroutine create_face_flux_fixT
   
   ! -----------------------------------------------------------------
   ! Subroutine used to the velocity on the faces of each grid cell.
   ! This subroutine translates to 3D the 2D velocity field obtained
   ! from the solution of the shallow water equations
   ! -----------------------------------------------------------------  
   subroutine get_face_velocity(lo_phys, lo, hi, dx, &
                                vx, vx_lo, vx_hi, &
                                vz, vz_lo, vz_hi, &
                                idom, id_lo, id_hi)

    use amr_data_module, only : surf_dx, &
                                surf_xlo, &
                                surf_ind, &
                                melt_vel
    
    ! Input and output variables
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: id_lo(3), id_hi(3)
    integer, intent(in) :: vx_lo(3), vx_hi(3)
    integer, intent(in) :: vz_lo(3), vz_hi(3)
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(in) :: lo_phys(3)
    real(amrex_real), intent(in)  :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2),id_lo(3):id_hi(3))
    real(amrex_real), intent(out) :: vx(vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2),vx_lo(3):vx_hi(3))
    real(amrex_real), intent(out) :: vz(vz_lo(1):vz_hi(1),vz_lo(2):vz_hi(2),vz_lo(3):vz_hi(3))

    ! Local variables
    integer :: i,j,k
    integer :: xind, zind
    real(amrex_real) :: xpos, zpos
    
    ! Assign x-component of the velocity 
    do i = lo(1), hi(1)+1
       do j = lo(2), hi(2)
          do k = lo(3), hi(3)
             
             if (nint(idom(i,j,k)).ge.2 .and. nint(idom(i-1,j,k)).ge.2) then

                ! Interpolation similar to the one used to find the free
                ! surface position in the heat transfer domain. See
                ! subroutine get_surf_pos in the domain module
                xpos = lo_phys(1) + (0.5 + i-lo(1))*dx(1) 
                zpos = lo_phys(3) + (0.5 + k-lo(3))*dx(3)
                xind = nint((xpos - surf_dx(1)/2 - surf_xlo(1))/surf_dx(1)) 
                zind = nint((zpos - surf_dx(2)/2 - surf_xlo(2))/surf_dx(2))
                if (xind.lt.surf_ind(1,1)) xind = surf_ind(1,1)
                if (xind.gt.surf_ind(1,2)) xind = surf_ind(1,2) 
                if (zind.lt.surf_ind(2,1)) zind = surf_ind(2,1)
                if (zind.gt.surf_ind(2,2)) zind = surf_ind(2,2) 
                vx(i,j,k) = melt_vel(xind,zind,1)
               !  vx(i,j,k) = melt_vel(i,k,1)
                
             else
                vx(i,j,k) = 0.0_amrex_real
             end if
             
          end do
       end do   
    end do

    ! Assign z-component of the velocity
    do i = lo(1), hi(1)
       do j = lo(2), hi(2)
          do k = lo(3), hi(3)+1
             
             if (nint(idom(i,j,k)).ge.2 .and. nint(idom(i,j,k-1)).ge.2) then

                ! Interpolation similar to the one used to find the free
                ! surface position in the heat transfer domain. See
                ! subroutine get_surf_pos in the domain module
                xpos = lo_phys(1) + (0.5 + i-lo(1))*dx(1) 
                zpos = lo_phys(3) + (0.5 + k-lo(3))*dx(3)
                xind = nint((xpos - surf_dx(1)/2 - surf_xlo(1))/surf_dx(1)) 
                zind = nint((zpos - surf_dx(2)/2 - surf_xlo(2))/surf_dx(2))
                if (xind.lt.surf_ind(1,1)) xind = surf_ind(1,1)
                if (xind.gt.surf_ind(1,2)) xind = surf_ind(1,2) 
                if (zind.lt.surf_ind(2,1)) zind = surf_ind(2,1)
                if (zind.gt.surf_ind(2,2)) zind = surf_ind(2,2) 
                vz(i,j,k) = melt_vel(xind,zind,2)
               !  vz(i,j,k) = melt_vel(i,k,2)
                
             else
                vz(i,j,k) = 0.0_amrex_real
             end if

          end do
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
     use domain_module, only : get_melt_pos, get_idomain
 
     ! Input and output variables
     real(amrex_real), intent(in) :: dt
     real(amrex_real), intent(in) :: time
   
     ! Local variables
     logical :: nodal(3)
     integer, parameter :: nghost = 1 ! number of ghost points in each spatial direction 
     integer :: ncomp, srccmp, dstcmp
     integer :: idim
     integer :: ilev
     type(amrex_multifab) :: phi_tmp(0:amrex_max_level) ! Enthalpy multifab with ghost points
     type(amrex_multifab) :: temp_tmp(0:amrex_max_level) ! Temperature multifab with ghost points
     type(amrex_multifab) :: phi_tmp2(0:amrex_max_level) ! Enthalpy multifab with two ghost points
     type(amrex_multifab) :: temp_tmp2(0:amrex_max_level) ! Temperature multifab with two ghost points
     type(amrex_multifab) :: idomain_tmp(0:amrex_max_level) ! Idomain multifab to distinguish material and background
     type(amrex_multifab) :: rhs(0:amrex_max_level) ! multifab for the right hand side of the linear system of equations
     type(amrex_multifab) :: alpha(0:amrex_max_level) ! multifab for the first term on the left hand side of the linear system of equations
     type(amrex_multifab) :: beta(amrex_spacedim,0:amrex_max_level) ! multifab for the second term on the left hand side of the linear system of equations
     type(amrex_boxarray) :: ba(0:amrex_max_level) ! Box array
     type(amrex_distromap):: dm(0:amrex_max_level) ! Distribution mapping
     type(amrex_fab) :: flux(amrex_spacedim)
     type(amrex_geometry) :: geom
     type(amrex_mfiter) :: mfi ! Multifab iterator
     type(amrex_box) :: bx
     real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidom
     real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptemp_tmp
     real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptemp
     real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pout
         
     ! Loop through all the levels
     do ilev = 0, amrex_max_level
 
        ! Swap old and new enthalpies before updating solution
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
        call amrex_multifab_build(phi_tmp2(ilev), ba(ilev), dm(ilev), ncomp, nghost+1) 
        call amrex_multifab_build(temp_tmp2(ilev), ba(ilev), dm(ilev), ncomp, nghost+1)
        call fillpatch(ilev, time, phi_tmp(ilev))

       ! Fill in the temperature and enthalpy ultifabs with 2 ghost points from previous solution  
        srccmp = 1
        dstcmp = 1
        call temp_tmp2(ilev)%copy(temp(ilev), srccmp, dstcmp, ncomp, nghost+1)
        call temp_tmp2(ilev)%fill_boundary(geom)
        call phi_tmp2(ilev)%copy(phi_old(ilev), srccmp, dstcmp, ncomp, nghost+1)
        call phi_tmp2(ilev)%fill_boundary(geom)
       
        ! Build temporary idomain multifab to store the domain configuration before SW deformation
        call amrex_multifab_build(idomain_tmp(ilev), ba(ilev), dm(ilev), ncomp, nghost)
        call amrex_multifab_swap(idomain_tmp(ilev), idomain(ilev))
 
     
        ! Multifabs for the linear solver
        call amrex_multifab_build(rhs(ilev), ba(ilev), dm(ilev), ncomp, 0)
        call amrex_multifab_build(alpha(ilev), ba(ilev), dm(ilev), ncomp, 0)
        do idim = 1, amrex_spacedim
           nodal = .false.
           nodal(idim) = .true.
           call amrex_multifab_build(beta(idim,ilev), ba(ilev), dm(ilev), ncomp, 0, nodal)
        end do

      !$omp parallel private(mfi, flux)
        ! Initialize fluxes (used in the predictor step)
        do idim = 1, amrex_spacedim
         call flux(idim)%reset_omp_private()
        end do
        ! Loop through all the boxes in the level
        call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.false.)
        do while(mfi%next())
           call predict_heat_solver_implicit_box(ilev, time, dt, mfi, &
                                                 geom, ncomp, phi_tmp(ilev), temp_tmp(ilev), &
                                                 idomain_tmp(ilev), flux, alpha(ilev), beta(:,ilev), &
                                                 rhs(ilev), phi_tmp2(ilev), temp_tmp2(ilev))
 
        end do
        call amrex_mfiter_destroy(mfi)
      !$omp end parallel
 
     end do
     
 
     ! Compute new temperature via implicit update
     call get_temperature_implicit(ba, dm, dt, alpha, beta, rhs, temp_tmp)
     
     ! Re-compute enthalpy and temperature on all levels (correction step)
     do ilev = 0, amrex_max_level
 
        ! Geometry
        geom = amrex_geom(ilev)
        
       !$omp parallel private(mfi)
        ! Loop through all the boxes in the level
        call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.false.)
        do while(mfi%next())
           call correct_heat_solver_implicit_box(ilev, mfi, &
                                                 phi_tmp(ilev), temp_tmp(ilev), &
                                                 alpha(ilev))
           
        end do
        call amrex_mfiter_destroy(mfi)
       !$omp end parallel
        call temp_tmp(ilev)%fill_boundary(geom) ! Fill in the solution to the ghost points
        
     end do

    ! Update the idomains accordingly
     do ilev = 0, amrex_max_level

      ! Geometry
      geom = amrex_geom(ilev)
      
     !$omp parallel private(mfi, bx, pidom, ptemp_tmp)
      ! Loop through all the boxes in the level
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

     ! Explicit update of the temperature due to heat advection
     do ilev = 0, amrex_max_level

         ! Geometry
         geom = amrex_geom(ilev)

         ! Loop through all the boxes in the level
        !$omp parallel private(mfi, bx, pidom, ptemp, ptemp_tmp, pout)
         call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.false.)
         do while(mfi%next())
            bx = mfi%validbox()
            pidom  => idomain(ilev)%dataptr(mfi)
            ptemp   => temp(ilev)%dataptr(mfi)
            ptemp_tmp   => temp_tmp(ilev)%dataptr(mfi)
            pout    => phi_new(ilev)%dataptr(mfi)
            call solve_heat_advection_box(geom%get_physical_location(bx%lo), bx%lo, bx%hi, dt, geom%dx(1:3), &
                                          lbound(pidom), ubound(pidom), pidom, &
                                          pout, lbound(pout), ubound(pout), &
                                          lbound(ptemp_tmp), ubound(ptemp_tmp), ptemp_tmp, &
                                          lbound(ptemp), ubound(ptemp), ptemp)
            
         end do
         call amrex_mfiter_destroy(mfi)
        !$omp end parallel   

     end do        

      ! Clean memory
      do ilev = 0, amrex_max_level
         call amrex_multifab_destroy(phi_tmp(ilev))
         call amrex_multifab_destroy(temp_tmp(ilev))
         call amrex_multifab_destroy(phi_tmp2(ilev))
         call amrex_multifab_destroy(temp_tmp2(ilev))
         call amrex_multifab_destroy(idomain_tmp(ilev))
         call amrex_multifab_destroy(rhs(ilev))
         call amrex_multifab_destroy(alpha(ilev))
         do idim = 1, amrex_spacedim
            call amrex_multifab_destroy(beta(idim,ilev))
         end do
         call amrex_distromap_destroy(dm)
      end do
 
     
   end subroutine advance_heat_solver_implicit
 
 
   ! -----------------------------------------------------------------
   ! Subroutine used to perform the prediction step for the implicit
   ! enthalpy update at a given box on a given level 
   ! -----------------------------------------------------------------
   subroutine predict_heat_solver_implicit_box(lev, time, dt, mfi, &
                                               geom, ncomp, phi_tmp, temp_tmp, &
                                               idomain_tmp, flux, alpha, beta, &
                                               rhs,  phi_tmp2, temp_tmp2)
 
     use read_input_module, only : temp_fs
     use amr_data_module, only : phi_new, & 
                                 temp, & 
                                 idomain, &
                                 Qplasma, &
                                 Qpipe, &
                                 Qtherm, &
                                 Qvap, &
                                 Qrad
     use domain_module, only : get_idomain, revaluate_heat_domain
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
     type(amrex_multifab), intent(inout) :: phi_tmp2
     type(amrex_multifab), intent(inout) :: temp_tmp2
     type(amrex_multifab), intent(inout) :: idomain_tmp
     type(amrex_multifab), intent(inout) :: alpha
     type(amrex_multifab), intent(inout) :: beta(amrex_spacedim)
     type(amrex_multifab), intent(inout) :: rhs
     type(amrex_fab), intent(inout) :: flux(amrex_spacedim)
 
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
     real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pfz
     real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidin
     real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidout
     real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pac
     real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pbc
     real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: prhs
     type(amrex_box) :: tbx
     type(amrex_box) :: bx
     real(amrex_real) :: Qpipe_box, Qplasma_box, Qtherm_box, Qrad_box, Qvap_box
     
     ! Box
     bx = mfi%validbox()   
           
     ! Pointers
     pin     => phi_tmp%dataptr(mfi)
     pin2     => phi_tmp2%dataptr(mfi)
     pout    => phi_new(lev)%dataptr(mfi)
     ptemp   => temp(lev)%dataptr(mfi)
     ptempin => temp_tmp%dataptr(mfi)
     ptempin2 => temp_tmp2%dataptr(mfi)
     pidin   => idomain_tmp%dataptr(mfi)
     pidout  => idomain(lev)%dataptr(mfi)
     pac     => alpha%dataptr(mfi)
     prhs    => rhs%dataptr(mfi)
     do idim = 1, amrex_spacedim
        tbx = bx
        call tbx%nodalize(idim)
        call flux(idim)%resize(tbx,ncomp)
        call tbx%grow(1)
     end do
     pfx => flux(1)%dataptr()
     pfy => flux(2)%dataptr()
     pfz => flux(3)%dataptr()
     
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
     call revaluate_heat_domain(lev, geom%get_physical_location(bx%lo), geom%dx, &
                                bx%lo, bx%hi, &
                                pidin, lbound(pidin), ubound(pidin), &
                                pidout, lbound(pidout), ubound(pidout), &
                                pin, lbound(pin), ubound(pin), &
                                ptempin, lbound(ptempin), ubound(ptempin), &
                                pin2, lbound(pin2), ubound(pin2), &
                                ptempin2, lbound(ptempin2), ubound(ptempin2))
     
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
     
   end subroutine predict_heat_solver_implicit_box
 
 
   ! -----------------------------------------------------------------
   ! Subroutine used to perform the correction step for the implicit
   ! enthalpy update at a given box on a given level 
   ! -----------------------------------------------------------------
   subroutine correct_heat_solver_implicit_box(lev, mfi, phi_tmp, &
                                               temp_tmp, alpha)
 
     use amr_data_module, only : phi_new, temp, idomain
     use domain_module, only : get_melt_pos
     use material_properties_module, only : get_temp
     use read_input_module, only : temp_fs
     
     ! Input and output variables
     integer, intent(in) :: lev
     type(amrex_mfiter), intent(in) :: mfi
     type(amrex_multifab), intent(inout) :: phi_tmp
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
        call get_enthalpy_implicit_fixT(bx%lo, bx%hi, &
                                        pin, lbound(pin), ubound(pin), &
                                        pout, lbound(pout), ubound(pout), &
                                        ptemp, lbound(ptemp), ubound(ptemp), &
                                        ptempin, lbound(ptempin), ubound(ptempin), &
                                        pid, lbound(pid), ubound(pid), &
                                        pac, lbound(pac), ubound(pac))
     else
        call get_enthalpy_implicit(bx%lo, bx%hi, &
                                   pin, lbound(pin), ubound(pin), &
                                   pout, lbound(pout), ubound(pout), &
                                   ptemp, lbound(ptemp), ubound(ptemp), &
                                   ptempin, lbound(ptempin), ubound(ptempin), &
                                   pac, lbound(pac), ubound(pac))
     end if
 
     
   end subroutine correct_heat_solver_implicit_box
   
   
   ! -----------------------------------------------------------------------------------
   ! Subroutine used to solve the linear system of equations for the implict update
   ! of the temperature
   ! ----------------------------------------------------------------------------------- 
   subroutine get_temperature_implicit(ba, dm, dt, alpha, beta, rhs, sol)
 
     use read_input_module, only : ls_composite_solve, ls_agglomeration, &
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
     if (ls_composite_solve) then
        
        call amrex_abeclaplacian_build(abeclap, amrex_geom, ba, dm, &
                                       metric_term=.false., agglomeration=ls_agglomeration, consolidation=ls_consolidation, &
                                       max_coarsening_level=ls_max_coarsening_level)
 
        call abeclap%set_maxorder(ls_linop_maxorder)
        
        ! This is set up to have homogeneous Neumann BC
        call abeclap%set_domain_bc([amrex_lo_neumann, amrex_lo_neumann, amrex_lo_neumann], &
                                   [amrex_lo_neumann, amrex_lo_neumann, amrex_lo_neumann])
        
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
        err = multigrid%solve(sol, rhs, 1.e-5_amrex_real, 0.0_amrex_real)
        
        ! Clean memory
        call amrex_abeclaplacian_destroy(abeclap)
        call amrex_multigrid_destroy(multigrid)
 
     else
 
        do ilev = 0, amrex_max_level
 
           call amrex_abeclaplacian_build(abeclap, [amrex_geom(ilev)], [ba(ilev)], [dm(ilev)], &
                                         metric_term=.false., agglomeration=ls_agglomeration, consolidation=ls_consolidation, &
                                         max_coarsening_level=ls_max_coarsening_level)
           
           call abeclap%set_maxorder(ls_linop_maxorder)
 
           ! This is set up to have homogeneous Neumann BC
           call abeclap%set_domain_bc([amrex_lo_neumann, amrex_lo_neumann, amrex_lo_neumann], &
                                      [amrex_lo_neumann, amrex_lo_neumann, amrex_lo_neumann])
        
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
           err = multigrid%solve([sol(ilev)], [rhs(ilev)], 1.e-5_amrex_real, 0.0_amrex_real)
        
           ! Clean memory
           call amrex_abeclaplacian_destroy(abeclap)
           call amrex_multigrid_destroy(multigrid)
           
        end do
        
     end if
   end subroutine get_temperature_implicit
 
     
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
   ! Subroutine used to get the right hand side of the linear solver.
   ! Case with fixed temperature at the free surface
   ! -----------------------------------------------------------------  
   subroutine get_rhs_fixT(lo, hi, &
                           temp, t_lo, t_hi, &
                           alpha, a_lo, a_hi, &
                           rhs, r_lo, r_hi)
     
     ! Input and output variables
     integer, intent(in) :: lo(3), hi(3)
     integer, intent(in) :: t_lo(3), t_hi(3)
     integer, intent(in) :: r_lo(3), r_hi(3)
     integer, intent(in) :: a_lo(3), a_hi(3)
     real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
     real(amrex_real), intent(in) :: alpha(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))   
     real(amrex_real), intent(out) :: rhs(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
     
     ! Local variables
     integer :: i,j,k
 
     ! Fill rhs of linear problem
     do i = lo(1), hi(1)
        do j = lo(2), hi(2)   
           do k = lo(3), hi(3)       
              rhs(i,j,k) = temp(i,j,k)*alpha(i,j,k)
           end do
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
   ! Subroutine used to get the corrected enthalpy. Case with fixed
   ! temperature at the free surface
   ! -----------------------------------------------------------------  
   subroutine get_enthalpy_implicit_fixT(lo, hi, &
                                         u_old, uo_lo, uo_hi, &
                                         u_new, un_lo, un_hi, &
                                         temp_old, to_lo, to_hi, &
                                         temp_new, tn_lo, tn_hi, &
                                         idom, id_lo, id_hi, &
                                         alpha, a_lo, a_hi)
 
     use material_properties_module, only : get_enthalpy, &
                                            get_temp, &
                                            temp_melt, &
                                            enth_at_melt, &
                                            latent_heat
     use read_input_module, only : temp_fs
     
     ! Input and output variables
     integer, intent(in) :: lo(3), hi(3)  
     integer, intent(in) :: uo_lo(3), uo_hi(3)
     integer, intent(in) :: un_lo(3), un_hi(3)
     integer, intent(in) :: to_lo(3), to_hi(3)
     integer, intent(in) :: tn_lo(3), tn_hi(3)
     integer, intent(in) :: id_lo(3), id_hi(3)
     integer, intent(in) :: a_lo(3), a_hi(3)
     real(amrex_real), intent(in) :: u_old(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3))
     real(amrex_real), intent(out) :: u_new(un_lo(1):un_hi(1),un_lo(2):un_hi(2),un_lo(3):un_hi(3))
     real(amrex_real), intent(inout) :: temp_old(to_lo(1):to_hi(1),to_lo(2):to_hi(2),to_lo(3):to_hi(3))
     real(amrex_real), intent(inout) :: temp_new(tn_lo(1):tn_hi(1),tn_lo(2):tn_hi(2),tn_lo(3):tn_hi(3))
     real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2),id_lo(3):id_hi(3))
     real(amrex_real), intent(in) :: alpha(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))
     
     ! Local variables
     integer :: i,j,k
     real(amrex_real) :: u_fs
     real(amrex_real) :: ubase
     real(amrex_real) :: u_tmp(un_lo(1):un_hi(1),un_lo(2):un_hi(2),un_lo(3):un_hi(3))
 
     ! Enthalpy at the free surface
     call get_enthalpy(temp_fs, u_fs)

    ! Get enthalpy in the general case, no transition of free surface
     call get_temp(un_lo, un_hi, &
                   u_tmp, un_lo, un_hi, &
                   temp_new, tn_lo, tn_hi, .false.)

     ! Get new enthalpy
     do i = lo(1), hi(1)
        do j = lo(2), hi(2)
           do k = lo(3), hi(3)
              if (nint(idom(i,j,k)).ne.0 .and. nint(idom(i,j+1,k)).eq.0) then
                 u_new(i,j,k) = u_fs ! Impose temperature on the free surface
              else if(nint(idom(i,j,k)).eq.0) then
                 u_new(i,j,k) = u_fs ! Set background temperature equal to the free surface temperature
              else
                 if (temp_old(i,j,k).le.temp_melt .and. temp_new(i,j,k).gt.temp_melt) then
                    ubase = max(u_old(i,j,k), enth_at_melt)
                    u_new(i,j,k) = ubase + alpha(i,j,k)*(temp_new(i,j,k) - temp_melt)
                 elseif (temp_old(i,j,k).ge.temp_melt .and. temp_new(i,j,k).lt.temp_melt) then
                    ubase = min(u_old(i,j,k), enth_at_melt+latent_heat)
                    u_new(i,j,k) = ubase + alpha(i,j,k)*(temp_new(i,j,k) - temp_melt)
                 else 
                    u_new(i,j,k) = u_tmp(i,j,k)
                 end if
            end if
           end do
        end do
     end do
 
   end subroutine get_enthalpy_implicit_fixT

  ! -----------------------------------------------------------------  
  ! Subroutine used for an explicit update of the temperature
  ! equation only with heat advection and no heat conduction.
  ! -----------------------------------------------------------------  
  subroutine solve_heat_advection_box(lo_phys, lo, hi, dt, dx, &
                                      id_lo, id_hi, idom, &
                                      u_new, un_lo, un_hi, &
                                      to_lo, to_hi, temp_old, &
                                      t_lo, t_hi, temp)

    use material_properties_module, only : get_temp, get_enthalpy

    ! Input and output variables
    integer, intent(in) :: lo(3)
    integer, intent(in) :: hi(3)
    integer, intent(in) :: t_lo(3)
    integer, intent(in) :: t_hi(3)
    integer, intent(in) :: to_lo(3)
    integer, intent(in) :: to_hi(3)
    integer, intent(in) :: id_lo(3)
    integer, intent(in) :: id_hi(3)
    integer, intent(in) :: un_lo(3)
    integer, intent(in) :: un_hi(3)
    real(amrex_real), intent(in) :: lo_phys(3)
    real(amrex_real), intent(in) :: dt
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(inout) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2),id_lo(3):id_hi(3))
    real(amrex_real), intent(inout) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(amrex_real), intent(in) :: temp_old(to_lo(1):to_hi(1),to_lo(2):to_hi(2),to_lo(3):to_hi(3))
    real(amrex_real), intent(inout) :: u_new(un_lo(1):un_hi(1),un_lo(2):un_hi(2),un_lo(3):un_hi(3))
    
    ! Local variables
    integer :: i,j,k
    integer :: ux_lo(3), ux_hi(3), uz_lo(3), uz_hi(3)
    real(amrex_real) :: ux(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3))
    real(amrex_real) :: uz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1)
    real(amrex_real) :: vx_l, vx_r, vz_l, vz_r
    real(amrex_real) :: partial_x
    real(amrex_real) :: partial_z

    
    ux_lo = lo
    ux_hi = hi
    ux_hi(1) = hi(1)+1
    uz_lo = lo
    uz_hi = hi
    uz_hi(3) = hi(3)+1
    
    ! Construct 2D melt velocity profile from the 2D shallow water solution
    call get_face_velocity(lo_phys, lo, hi, dx, &
                           ux, ux_lo, ux_hi, &
                           uz, uz_lo, uz_hi, &
                           idom, id_lo, id_hi)

    ! Update temperature profile
    do i = lo(1), hi(1)
       do j = lo(2), hi(2)
         do k = lo(3), hi(3)
            vx_l = ux(i,j,k)
            vx_r = ux(i+1,j,k)
            vz_l = uz(i,j,k)
            vz_r = uz(i,j,k+1)

            if (nint(idom(i,j,k)).ge.2) then
               partial_x = 0.0_amrex_real
               partial_z = 0.0_amrex_real
               ! Partial update based on motion in the x direction
               if (vx_l.gt.0 .and. nint(idom(i-1,j,k)).ge.2) then
                  partial_x = - dt/dx(1) * (vx_l+ABS(vx_l))*(temp_old(i,j,k)-temp_old(i-1,j,k))/2.0_amrex_real
               end if
               if (vx_r.lt.0 .and. nint(idom(i+1,j,k)).ge.2) then
                  partial_x = partial_x - dt/dx(1) * (vx_r-ABS(vx_r))*(temp_old(i+1,j,k)-temp_old(i,j,k))/2.0_amrex_real
               end if
               ! Partial update based on motion in the z direction
               if ((vz_l.gt.0 .and. nint(idom(i,j,k-1)).ge.2)) then
                  partial_z =  - dt/dx(3) * (vz_l+ABS(vz_l))*(temp_old(i,j,k)-temp_old(i,j,k-1))/2.0_amrex_real
               end if
               if (vz_r.lt.0 .and. nint(idom(i,j,k+1)).ge.2) then
                  partial_z = partial_z - dt/dx(3) * (vz_r-ABS(vz_r))*(temp_old(i,j,k+1)-temp_old(i,j,k))/2.0_amrex_real
               end if
               temp(i,j,k) = temp_old(i,j,k) + partial_x + partial_z
            end if
            if(temp(i,j,k).ne.temp_old(i,j,k)) then
               call get_enthalpy(temp(i,j,k),u_new(i,j,k))
            end if
            
         end do
       end do
    end do
    
  end subroutine solve_heat_advection_box  
   
 end module heat_transfer_module
 