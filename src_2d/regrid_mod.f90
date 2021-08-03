module regrid_module

  ! -----------------------------------------------------------------
  ! This module is used to perform the regridding and to deal with
  ! the different mesh levels generated by amrex
  ! -----------------------------------------------------------------
  
  use iso_c_binding
  use amrex_amr_module
  use amr_data_module, only : idomain_new, &
                              idomain_old, &
                              flux_reg, &
                              phi_new, &
                              phi_old, &
                              temp, &
                              t_new, &
                              t_old
  
  implicit none

  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: averagedown
  public :: averagedownto
  public :: fillpatch
  public :: my_clear_level
  public :: my_error_estimate
  public :: my_make_new_level_from_scratch
  public :: my_make_new_level_from_coarse
  public :: my_remake_level
  
  ! -----------------------------------------------------------------
  ! Declare private variables shared by all subroutines
  ! -----------------------------------------------------------------  
  integer, parameter :: ncomp = 1
  integer, parameter :: nghost = 0
  
contains

  ! NOTE: the bind(c) command is necessary for all those functions
  ! that have to be called directly from the underlying c++ code.
  
  ! -----------------------------------------------------------------
  ! Subroutine used to create a new level and fill it with the
  ! initialization values specified in input
  ! -----------------------------------------------------------------
  subroutine my_make_new_level_from_scratch(lev, time, pba, pdm) bind(c)

    use read_input_module, only : tempinit, do_reflux
    use material_properties_module, only : get_temp
    use heat_transfer_module, only : get_idomain

    ! Input and output variables
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm

    ! Local variables
    real(amrex_real), contiguous, pointer :: pid(:,:,:,:)
    real(amrex_real), contiguous, pointer :: phi(:,:,:,:)
    real(amrex_real), contiguous, pointer :: ptemp(:,:,:,:)
    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    type(amrex_geometry) :: geom
    
    ! Pointers for box array and distribution mapping
    ba = pba
    dm = pdm

    ! Geometry
    geom = amrex_geom(lev)
    
    ! Time
    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real

    ! Clear any previous definition of the level
    call my_clear_level(lev)

    ! Build the multifabs
    call amrex_multifab_build(phi_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_old(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(temp(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(idomain_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(idomain_old(lev), ba, dm, ncomp, nghost)

    ! Build the flux registers
    if (lev > 0 .and. do_reflux) then
       call amrex_fluxregister_build(flux_reg(lev), ba, dm, &
                                     amrex_ref_ratio(lev-1), lev, ncomp)
    end if

    ! Fill the multifabs with the initialization values
    call amrex_mfiter_build(mfi, phi_new(lev))
    do while (mfi%next())
       
       bx = mfi%tilebox()
       phi => phi_new(lev)%dataptr(mfi)
       ptemp => temp(lev)%dataptr(mfi)
       pid => idomain_new(lev)%dataptr(mfi)
       
       ! Enthalpy
       call init_phi(bx%lo, bx%hi, tempinit, &
                      phi, lbound(phi), ubound(phi))
       ! Temperature
       call get_temp(bx%lo, bx%hi, &
                     phi, lbound(phi), ubound(phi), &
                     ptemp, lbound(ptemp), ubound(ptemp))
       
       ! Integer domain to distinguish material and background
       call get_idomain(geom%get_physical_location(bx%lo), geom%dx, &
                        bx%lo, bx%hi, &
                        pid, lbound(pid), ubound(pid))
       
       
    end do
    call amrex_mfiter_destroy(mfi)

  end subroutine my_make_new_level_from_scratch

  ! -----------------------------------------------------------------
  ! Subroutine used to initialize the enthalpy multifab
  ! -----------------------------------------------------------------
  subroutine init_phi(lo, hi, tempinit, &
                      phi, phi_lo, phi_hi)

    use material_properties_module, only : get_enthalpy

    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: phi_lo(2), phi_hi(2)
    real(amrex_real), intent(inout) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2))
    real(amrex_real), intent(in) :: tempinit

    ! Local variables
    integer          :: i,j
    real(amrex_real) :: enth_init

    ! Get enthalpy consistent with the initialization temperature
    call get_enthalpy(tempinit,enth_init)

    !$omp parallel do private(i,j) collapse(2)
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          phi(i,j) = enth_init
       end do
    end do
    !$omp end parallel do
    
    
  end subroutine init_phi

  
  ! -----------------------------------------------------------------
  ! Subroutine used to create a new level and fill it with the values
  ! obtained from a coarser level
  ! -----------------------------------------------------------------
  subroutine my_make_new_level_from_coarse(lev, time, pba, pdm) bind(c)

    use read_input_module, only : do_reflux
    use material_properties_module, only : get_temp
    use heat_transfer_module, only : get_idomain

    ! Input and output variables    
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm

    ! Local variables
    real(amrex_real), contiguous, pointer :: pid(:,:,:,:)
    real(amrex_real), contiguous, pointer :: phi(:,:,:,:)
    real(amrex_real), contiguous, pointer :: ptemp(:,:,:,:)
    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    type(amrex_geometry) :: geom
    
    ! Pointers for box array and distribution mapping
    ba = pba
    dm = pdm

    ! Geometry
    geom = amrex_geom(lev)
    
    ! Time
    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real

    ! Clear any previous definition of the level
    call my_clear_level(lev)

    ! Build the multifabs
    call amrex_multifab_build(phi_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_old(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(temp(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(idomain_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(idomain_old(lev), ba, dm, ncomp, nghost)

    ! Build the flux registers
    if (lev > 0 .and. do_reflux) then
       call amrex_fluxregister_build(flux_reg(lev), ba, dm, &
                                     amrex_ref_ratio(lev-1), lev, ncomp)
    end if
    
    ! Fill enthalpy multifab from coarser level
    call fillcoarsepatch(lev, time, phi_new(lev))
    
    ! Fill the temperature and idomain multifabs
    call amrex_mfiter_build(mfi, temp(lev))    
    do while (mfi%next())

       bx = mfi%tilebox()
       phi => phi_new(lev)%dataptr(mfi)
       ptemp => temp(lev)%dataptr(mfi)
       pid => idomain_new(lev)%dataptr(mfi)
       
       ! Temperature
       call get_temp(bx%lo, bx%hi, & 
                     phi, lbound(phi), ubound(phi), &
                     ptemp, lbound(ptemp), ubound(ptemp))
       
       ! Integer domain to distinguish material and background
       call get_idomain(geom%get_physical_location(bx%lo), geom%dx, &
                        bx%lo, bx%hi, &
                        pid, lbound(pid), ubound(pid))
       
    end do    
    call amrex_mfiter_destroy(mfi)   
    
  end subroutine my_make_new_level_from_coarse


  ! -----------------------------------------------------------------
  ! Subroutine used to fill the enthalpy multifab at a given level
  ! when that level is created from a coarser one
  ! ----------------------------------------------------------------- 
  subroutine fillcoarsepatch(lev, time, phi)

    use amr_data_module, only : t_old, t_new, phi_old, phi_new, lo_bc, hi_bc

    ! Input and output variables
    integer, intent(in) :: lev
    real(amrex_real), intent(in) :: time
    type(amrex_multifab), intent(inout) :: phi

    
    call amrex_fillcoarsepatch(phi, t_old(lev-1), phi_old(lev-1),  &
                               t_new(lev-1), phi_new(lev-1),  &
                               amrex_geom(lev-1),    fill_physbc,  &
                               amrex_geom(lev  ),    fill_physbc,  &
                               time, ncomp, ncomp, ncomp, &
                               amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
                               lo_bc, hi_bc)
    
  end subroutine fillcoarsepatch

  ! -----------------------------------------------------------------
  ! Subroutine used to apply the boundary conditions
  ! ----------------------------------------------------------------- 
  subroutine fill_physbc(pmf, scomp, ncomp, time, pgeom) bind(c)

    use amrex_filcc_module, only : amrex_filcc
    use amr_data_module, only : lo_bc, hi_bc

    ! Input and output variables
    type(c_ptr), value :: pmf, pgeom
    integer(c_int), value :: scomp, ncomp
    real(amrex_real), value :: time

    ! Local variables
    type(amrex_geometry) :: geom
    type(amrex_multifab) :: mf
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: p
    integer :: plo(4), phi(4)

    ! Pointers to the geometry and multifab to which the boundary
    ! conditions must be applied
    geom = pgeom
    mf = pmf
	
	
    !$omp parallel private(mfi,p,plo,phi)
    call amrex_mfiter_build(mfi, mf, tiling=.false.)
    do while(mfi%next())

       p => mf%dataptr(mfi)

       if (.not. geom%domain%contains(p)) then ! part of this box is outside the domain

          plo = lbound(p)
          phi = ubound(p)
          call amrex_filcc(p, plo, phi,         & 
                           geom%domain%lo, geom%domain%hi,  & 
                           geom%dx,                         & 
                           geom%get_physical_location(plo), & 
                           lo_bc, hi_bc)                      
       
       end if
       
    end do
    call amrex_mfiter_destroy(mfi)   
    !$omp end parallel

  end subroutine fill_physbc 
  
  ! -----------------------------------------------------------------
  ! Subroutine used to remake an existing level and fill it with the
  ! values obtained from a coarser level
  ! -----------------------------------------------------------------
  subroutine my_remake_level(lev, time, pba, pdm) bind(c)

    use read_input_module, only : do_reflux
    use material_properties_module, only : get_temp
    use heat_transfer_module, only : get_idomain

    ! Input and output variables    
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm

    ! Local variables
    real(amrex_real), contiguous, pointer :: pid(:,:,:,:)
    real(amrex_real), contiguous, pointer :: phi(:,:,:,:)
    real(amrex_real), contiguous, pointer :: ptemp(:,:,:,:)
    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm
    type(amrex_mfiter) :: mfi
    type(amrex_multifab) :: new_phi_new
    type(amrex_box) :: bx
    type(amrex_geometry) :: geom
    
    ! Pointers for box array and distribution mapping
    ba = pba
    dm = pdm

    ! Geometry
    geom = amrex_geom(lev)
    
    ! Create a copy of phi_new and fill it with fillpatch
    call amrex_multifab_build(new_phi_new, ba, dm, ncomp, 0)
    call fillpatch(lev, time, new_phi_new)

    ! Time
    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real

    ! Clear previous definition of the level
    call my_clear_level(lev)

    ! Build the multifabs
    call amrex_multifab_build(phi_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_old(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(temp(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(idomain_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(idomain_old(lev), ba, dm, ncomp, nghost)

    ! Build the flux registers
    if (lev > 0 .and. do_reflux) then
       call amrex_fluxregister_build(flux_reg(lev), ba, dm, &
                                     amrex_ref_ratio(lev-1), lev, ncomp)
    end if

    ! Copy the content of new_phi_new into phi_new
    call phi_new(lev)%copy(new_phi_new, 1, 1, ncomp, 0)
    call amrex_multifab_destroy(new_phi_new)
    
    ! Fill temperature and idomain
    call amrex_mfiter_build(mfi, temp(lev))
    do while (mfi%next())
       
       bx = mfi%tilebox()
       phi => phi_new(lev)%dataptr(mfi)
       ptemp => temp(lev)%dataptr(mfi)
       pid => idomain_new(lev)%dataptr(mfi)

       ! Temperature
       call get_temp(bx%lo, bx%hi, & 
                     phi, lbound(phi), ubound(phi), &
                     ptemp, lbound(ptemp), ubound(ptemp))

       ! Integer domain to distinguish material and background
       call get_idomain(geom%get_physical_location(bx%lo), geom%dx, &
                        bx%lo, bx%hi, &
                        pid, lbound(pid), ubound(pid))
       
     end do    
    call amrex_mfiter_destroy(mfi)

    
  end subroutine my_remake_level

  ! -----------------------------------------------------------------
  ! Subroutine used to fill the enthalpy multifab at a given level
  ! ----------------------------------------------------------------- 
  subroutine fillpatch(lev, time, phi)

    use amr_data_module, only : t_old, t_new, phi_old, phi_new, lo_bc, hi_bc
    
    ! Input and output variables
    integer, intent(in) :: lev
    real(amrex_real), intent(in) :: time
    type(amrex_multifab), intent(inout) :: phi
    
    
    if (lev .eq. 0) then

       call amrex_fillpatch(phi, t_old(lev), phi_old(lev), &
                            t_new(lev), phi_new(lev), &
                            amrex_geom(lev), fill_physbc , &
                            time, ncomp, ncomp, ncomp)
       
    else
       
       call amrex_fillpatch(phi, t_old(lev-1), phi_old(lev-1), &
                            t_new(lev-1), phi_new(lev-1), &
                            amrex_geom(lev-1), fill_physbc   , &
                            t_old(lev  ), phi_old(lev  ), &
                            t_new(lev  ), phi_new(lev  ), &
                            amrex_geom(lev  ), fill_physbc   , &
                            time, ncomp, ncomp, ncomp, &
                            amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
                            lo_bc, hi_bc)

    end if
    
  end subroutine fillpatch

  
  ! -----------------------------------------------------------------
  ! Subroutine used to clear a given level by destroying all the
  ! multifabs defined on that level
  ! -----------------------------------------------------------------
  subroutine my_clear_level(lev) bind(c)
    
    integer, intent(in), value :: lev
    
    call amrex_multifab_destroy(phi_new(lev))
    call amrex_multifab_destroy(phi_old(lev))
    call amrex_multifab_destroy(temp(lev))
    call amrex_multifab_destroy(idomain_new(lev))
    call amrex_multifab_destroy(idomain_old(lev))
    call amrex_fluxregister_destroy(flux_reg(lev))
    
  end subroutine my_clear_level

  ! -----------------------------------------------------------------
  ! Subroutine used to construct the error estimate used to
  ! identify the grid points that need regridding
  ! -----------------------------------------------------------------
  subroutine my_error_estimate(lev, cp, t, settag, cleartag) bind(c)

    use read_input_module,  only : surfdist 

    ! Input and output variables
    character(kind=c_char), intent(in), value :: cleartag
    character(kind=c_char), intent(in), value :: settag
    integer, intent(in), value :: lev
    type(c_ptr), intent(in), value :: cp
    real(amrex_real), intent(in), value :: t

    ! Local variables
    character(kind=c_char), contiguous, pointer :: tagarr(:,:,:,:)
    real(amrex_real), contiguous, pointer :: phiarr(:,:,:,:)
    type(amrex_geometry) :: geom
    type(amrex_tagboxarray) :: tag
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    
    ! ??? 
    tag = cp

    ! Geometry
    geom = amrex_geom(lev) 

    ! Get error estimate
    !$omp parallel private(mfi, bx, phiarr, tagarr)
    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.false.)
    do while(mfi%next())
       
       bx = mfi%validbox()
       phiarr => phi_new(lev)%dataptr(mfi)
       tagarr => tag%dataptr(mfi)
       
       call tag_phi_error(bx%lo, bx%hi, &
                          geom%get_physical_location(bx%lo), &
                          geom%dx, surfdist(lev+1), & 
                          phiarr, lbound(phiarr), ubound(phiarr), &
                          tagarr, lbound(tagarr), ubound(tagarr), &
                          settag)
       
    end do
    call amrex_mfiter_destroy(mfi)
    !$omp end parallel

  end subroutine my_error_estimate


  ! -----------------------------------------------------------------
  ! Subroutine used to tag the grid points that need regridding
  ! -----------------------------------------------------------------  
  subroutine tag_phi_error(lo, hi, xlo, dx, surfdist, &
                           phi, philo, phihi, &
                           tag, taglo, taghi, &
                           settag)
    
    use heat_transfer_module, only : get_surf_pos   
    use material_properties_module, only : enth_at_melt

    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: philo(3), phihi(3) ! WHY DOES THIS HAVE 3 ELEMENTS?!
    integer, intent(in) :: taglo(3), taghi(3)
    real(amrex_real), intent(in) :: dx(2)  
    real(amrex_real), intent(in) :: phi(philo(1):phihi(1),philo(2):phihi(2))
    real(amrex_real), intent(in) :: xlo(2)
    real(amrex_real), intent(in) :: surfdist
    character(kind=c_char), intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2))
    character(kind=c_char), intent(in) :: settag
    
    ! Local variables
    integer :: i,j
    real(amrex_real) :: surfpos(lo(1):hi(1)) 
    real(amrex_real) :: ydist


    ! Get position of the free surface
    call get_surf_pos(xlo, dx, lo, hi, surfpos)

    ! Loop through the domain
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          
          ! Regrid based on the distance from the free surface
          ydist = abs(xlo(2) + (j-lo(2))*dx(2) - surfpos(i) ) 
          if (ydist .le. surfdist) then 
             tag(i,j) = settag
          endif
          
          ! Regrid based on the state: all points belonging to the melt pool
          ! should be highly resolved 
          if (phi(i,j).ge.enth_at_melt) then
             tag(i,j) = settag  
          endif
          
       enddo
    enddo
    
  end subroutine tag_phi_error

  ! ------------------------------------------------------------------
  ! Subroutine to average down to the lowest level
  ! ------------------------------------------------------------------ 
  subroutine averagedown()

    use amr_data_module, only : phi_new
    
    integer :: lev, finest_level
    finest_level = amrex_get_finest_level()
    do lev = finest_level-1, 0, -1
       call amrex_average_down(phi_new(lev+1), phi_new(lev), &
                               amrex_geom(lev+1), amrex_geom(lev), &
                               1, 1, amrex_ref_ratio(lev))
    end do
    
  end subroutine averagedown

  ! ------------------------------------------------------------------
  ! Subroutine to average one level down
  ! ------------------------------------------------------------------ 
  subroutine averagedownto(clev)

    use amr_data_module, only : phi_new
    
    integer, intent(in) :: clev

    call amrex_average_down(phi_new(clev+1), phi_new(clev), &
                            amrex_geom(clev+1), amrex_geom(clev), &
                            1, 1, amrex_ref_ratio(clev))
    
  end subroutine averagedownto
  
end module regrid_module
