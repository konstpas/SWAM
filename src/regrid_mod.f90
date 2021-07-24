module regrid_module

  use iso_c_binding
  use amrex_amr_module

  use amr_data_module
  implicit none
  
  integer, private, parameter :: ncomp = 1, nghost = 0
  
contains


  ! Make a new level from scratch and put the data in phi_new.
  ! Note that phi_old contains no valid data after this.
  subroutine my_make_new_level_from_scratch (lev, time, pba, pdm) bind(c)

   
    use read_input_module, only : tempinit, surf_pos_init, do_reflux
    use material_properties_module, only : get_temp 
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm

    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: phi(:,:,:,:), ptemp(:,:,:,:)
    
    ba = pba
    dm = pdm

    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real

    call my_clear_level(lev)

    call amrex_multifab_build(phi_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_old(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(temp(lev), ba, dm, ncomp, nghost)
    call amrex_imultifab_build(idomain_new(lev), ba, dm, ncomp, nghost)
    call amrex_imultifab_build(idomain_old(lev), ba, dm, ncomp, nghost)


   if (lev > 0 .and. do_reflux) then
      call amrex_fluxregister_build(flux_reg(lev), ba, dm, amrex_ref_ratio(lev-1), lev, ncomp)
   end if

    call amrex_mfiter_build(mfi, phi_new(lev))
	
    do while (mfi%next())
       bx = mfi%tilebox()
       phi => phi_new(lev)%dataptr(mfi)
       ptemp => temp(lev)%dataptr(mfi)
       call init_phi(lev, t_new(lev), bx%lo, bx%hi, tempinit, phi, lbound(phi), ubound(phi), &
            amrex_geom(lev)%dx, amrex_problo)
       call get_temp(bx%lo, bx%hi, & 
       	lbound(phi),   ubound(phi),   phi, &
       	lbound(ptemp), ubound(ptemp), ptemp)
       
       ! call get domain integers 	
    end do

    call amrex_mfiter_destroy(mfi)
    

  end subroutine my_make_new_level_from_scratch

  ! Make a new level from coarse level and put the data in phi_new.
  ! Note tha phi_old contains no valid data after this.
  subroutine my_make_new_level_from_coarse (lev, time, pba, pdm) bind(c)
    use read_input_module, only : surf_pos_init, do_reflux
    use heat_transfer_module, only : get_idomain
    use material_properties_module, only : get_temp 
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm

    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: phi(:,:,:,:), ptemp(:,:,:,:)
    integer, contiguous, pointer :: idom(:,:,:,:) 

    ba = pba
    dm = pdm

    call my_clear_level(lev)

    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real

    call amrex_multifab_build(phi_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_old(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(temp(lev), ba, dm, ncomp, nghost)
    call amrex_imultifab_build(idomain_new(lev), ba, dm, ncomp, nghost)
    call amrex_imultifab_build(idomain_old(lev), ba, dm, ncomp, nghost)
    
    if (lev > 0 .and. do_reflux) then
       call amrex_fluxregister_build(flux_reg(lev), ba, dm, amrex_ref_ratio(lev-1), lev, ncomp)
    end if

     
    call fillcoarsepatch(lev, time, phi_new(lev))
    
    ! Fill temperature data 
    call amrex_mfiter_build(mfi, temp(lev))
    ! parallelize 
    do while (mfi%next())
       bx = mfi%tilebox()
       phi => phi_new(lev)%dataptr(mfi)
       ptemp => temp(lev)%dataptr(mfi)
       idom => idomain_new(lev)%dataptr(mfi)
       call get_temp(bx%lo, bx%hi, & 
       	lbound(phi),   ubound(phi),   phi, &
       	lbound(ptemp), ubound(ptemp), ptemp)  
       call get_idomain(bx%lo, bx%hi, & 
       		lbound(idom), ubound(idom), idom) 	
       ! call get domain integers 
    end do    
    
    call amrex_mfiter_destroy(mfi)
   
   
    
  end subroutine my_make_new_level_from_coarse

  ! Remake a level from current and coarse levels and put the data in phi_new.
  ! Note that phi_old contains no valid data after this.
  subroutine my_remake_level (lev, time, pba, pdm) bind(c)

    use read_input_module, only : surf_pos_init, do_reflux
    use material_properties_module, only : get_temp     
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm
    
    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm
    type(amrex_multifab) :: new_phi_new
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: phi(:,:,:,:), ptemp(:,:,:,:)
    

    ba = pba
    dm = pdm
	
    call amrex_multifab_build(new_phi_new, ba, dm, ncomp, 0)
    call fillpatch(lev, time, new_phi_new)

    call my_clear_level(lev)

    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real

    call amrex_multifab_build(phi_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_old(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(temp(lev), ba, dm, ncomp, nghost)
    call amrex_imultifab_build(idomain_new(lev), ba, dm, ncomp, nghost)
    call amrex_imultifab_build(idomain_old(lev), ba, dm, ncomp, nghost)
    
    if (lev > 0 .and. do_reflux) then
       call amrex_fluxregister_build(flux_reg(lev), ba, dm, amrex_ref_ratio(lev-1), lev, ncomp)
    end if

    call phi_new(lev)%copy(new_phi_new, 1, 1, ncomp, 0)
    
    call amrex_multifab_destroy(new_phi_new)
    
    ! Fill temperature data 
    call amrex_mfiter_build(mfi, temp(lev))
     do while (mfi%next())
       bx = mfi%tilebox()
       phi => phi_new(lev)%dataptr(mfi)
       ptemp => temp(lev)%dataptr(mfi)
       call get_temp(bx%lo, bx%hi, & 
       	lbound(phi),   ubound(phi),   phi, &
       	lbound(ptemp), ubound(ptemp), ptemp)   
     end do    
    call amrex_mfiter_destroy(mfi)

    
  end subroutine my_remake_level

  subroutine my_clear_level (lev) bind(c)
    integer, intent(in), value :: lev
    call amrex_multifab_destroy(phi_new(lev))
    call amrex_multifab_destroy(phi_old(lev))
    call amrex_multifab_destroy(temp(lev))
    call amrex_imultifab_destroy(idomain_new(lev))
    call amrex_imultifab_destroy(idomain_old(lev))
    call amrex_fluxregister_destroy(flux_reg(lev))
  end subroutine my_clear_level

  subroutine my_error_estimate (lev, cp, t, settag, cleartag) bind(c)

    use read_input_module,  only : surfdist 

    integer, intent(in), value :: lev
    type(c_ptr), intent(in), value :: cp
    real(amrex_real), intent(in), value :: t
    character(kind=c_char), intent(in), value :: settag, cleartag
    
    type(amrex_geometry) :: geom  	! geometry at level
    real(amrex_real), allocatable, save :: phierr(:)
    type(amrex_parmparse) :: pp
    type(amrex_tagboxarray) :: tag
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: phiarr(:,:,:,:)
    character(kind=c_char), contiguous, pointer :: tagarr(:,:,:,:)

    ! ! PART OF THE INPUT IS READ HERE!, REMOVE THIS!
    ! if (.not.allocated(phierr)) then
    !    call amrex_parmparse_build(pp, "myamr")
    !    call pp%getarr("phierr", phierr)
    !    call amrex_parmparse_destroy(pp)
    ! end if

    tag = cp
    geom = amrex_geom(lev) 

	
    !$omp parallel private(mfi, bx, phiarr, tagarr)
    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.true.)
    do while(mfi%next())
       bx = mfi%tilebox()
       phiarr => phi_new(lev)%dataptr(mfi)
       tagarr => tag%dataptr(mfi)
       ! call tag_phi_error(lev, t, bx%lo, bx%hi, &
       !      geom%get_physical_location(bx%lo), geom%dx, surfdist(lev+1), & 
       !      phiarr, lbound(phiarr), ubound(phiarr), &
       !      tagarr, lbound(tagarr), ubound(tagarr), &
       !      phierr(lev+1), settag, cleartag)  ! +1 because level starts with 0, but phierr starts with 1
       call tag_phi_error(lev, t, bx%lo, bx%hi, &
            geom%get_physical_location(bx%lo), geom%dx, surfdist(lev+1), & 
            phiarr, lbound(phiarr), ubound(phiarr), &
            tagarr, lbound(tagarr), ubound(tagarr), &
            settag, cleartag)
    end do
    call amrex_mfiter_destroy(mfi)
    !$omp end parallel

  end subroutine my_error_estimate


  subroutine init_phi(level, time, lo, hi, tempinit, phi, phi_lo, phi_hi, &
       dx, prob_lo)

    use amrex_fort_module, only : amrex_spacedim, amrex_real
    use material_properties_module, only : get_enthalpy
    
    integer, intent(in) :: level, lo(3), hi(3), phi_lo(3), phi_hi(3)
    real(amrex_real), intent(in) :: time, tempinit
    real(amrex_real), intent(inout) :: phi(phi_lo(1):phi_hi(1), &
         &                                 phi_lo(2):phi_hi(2), &
         &                                 phi_lo(3):phi_hi(3))
    real(amrex_real), intent(in) :: dx(3), prob_lo(3)
    
    integer          :: i,j,k
    real(amrex_real) :: x,y,z,r2
    real(amrex_real) :: enth_init

    call get_enthalpy(tempinit,enth_init)
    
    !$omp parallel do private(i,j,k,x,y,z,r2) collapse(2)
    ! do k=lo(3),hi(3)
    !    do j=lo(2),hi(2)
    !       z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
    !       y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
    !       do i=lo(1),hi(1)
    !          x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)
             
    !          if ( amrex_spacedim .eq. 2) then
    !             r2 = ((x-0.5d0)**2 + (y-0.75d0)**2) / 0.01d0
    !             !phi(i,j,k) = tempinit! 1.d0 !+ exp(-r2)
    !             phi(i,j,k) = enth_init
    !          else
    !             r2 = ((x-0.5d0)**2 + (y-0.75d0)**2 + (z-0.5d0)**2) / 0.01d0
    !             !phi(i,j,k) = tempinit! + exp(-r2)
    !             phi(i,j,k) = enth_init 
    !          end if
    !       end do
    !    end do
    ! end do
    !$omp end parallel do

    !$omp parallel do private(i,j,k,x,y,z,r2) collapse(2)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             phi(i,j,k) = enth_init
          end do
       end do
    end do
    !$omp end parallel do
    
    
  end subroutine init_phi
  
  
  subroutine tag_phi_error(level, time, lo, hi, xlo, dx, &
                           surfdist, phi, philo, phihi, tag, &
                           taglo, taghi, settag, cleartag)
    
    use heat_transfer_module, only : get_surf_pos   
    use material_properties_module, only : enth_at_melt  
    integer               , intent(in   ) :: level, lo(3), hi(3), philo(4), phihi(4), taglo(4), taghi(4)
    real(amrex_real)      , intent(in   ) :: phi(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
    real(amrex_real)      , intent(in   ) :: xlo(3), dx(3)  ! Surface y-position (to be multifab on DIM-1) 
    real(amrex_real)      , intent(in   ) :: surfdist ! distance from interface to refine grid
    character(kind=c_char), intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    !real(amrex_real)      , intent(in   ) :: time, phierr
    real(amrex_real)      , intent(in   ) :: time
    character(kind=c_char), intent(in   ) :: settag, cleartag
    real(amrex_real)                      :: surfpos(lo(1):hi(1),lo(3):hi(3))  ! Array surface position 

    integer :: i,j,k

    real(amrex_real) :: phierr_this
    real(amrex_real) :: ydist			! distance from free interface



    call get_surf_pos(xlo, dx, lo, hi, surfpos)
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
            
             ydist = abs(xlo(2) + (j-lo(2))*dx(2) - surfpos(i,k) ) 
             if (ydist .le. surfdist) then ! Free surface interface highly resolved 
                tag(i,j,k) = settag
             endif 
             if (phi(i,j,k).ge.enth_at_melt) then ! all molten elements highly resolved 
             	tag(i,j,k) = settag  
             endif 
             
          enddo
       enddo
    enddo
    
    
    
  end subroutine tag_phi_error


    ! Fill phi with data from phi_old and phi_new of current level and one level below.
  subroutine fillpatch (lev, time, phi)

    use amr_data_module, only : t_old, t_new, phi_old, phi_new, lo_bc, hi_bc
    
    ! Input and output variables
    integer, intent(in) :: lev
    real(amrex_real), intent(in) :: time
    type(amrex_multifab), intent(inout) :: phi

    ! Local variables
    integer, parameter :: src_comp=1, dst_comp=1, num_comp=1  ! for this test code


    if (lev .eq. 0) then

       call amrex_fillpatch(phi, t_old(lev), phi_old(lev), &
                            t_new(lev), phi_new(lev), &
                            amrex_geom(lev), fill_physbc , &
                            time, src_comp, dst_comp, num_comp)
       
    else
       
       call amrex_fillpatch(phi, t_old(lev-1), phi_old(lev-1), &
                            t_new(lev-1), phi_new(lev-1), &
                            amrex_geom(lev-1), fill_physbc   , &
                            t_old(lev  ), phi_old(lev  ), &
                            t_new(lev  ), phi_new(lev  ), &
                            amrex_geom(lev  ), fill_physbc   , &
                            time, src_comp, dst_comp, num_comp, &
                            amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
                            lo_bc, hi_bc)

    end if
    
  end subroutine fillpatch


  
  ! This is called when a new level is created from a coarser one.
  subroutine fillcoarsepatch (lev, time, phi)

    use amr_data_module, only : t_old, t_new, phi_old, phi_new, lo_bc, hi_bc

    ! Input and output variables
    integer, intent(in) :: lev
    real(amrex_real), intent(in) :: time
    type(amrex_multifab), intent(inout) :: phi

    ! Local variables
    integer, parameter :: src_comp=1, dst_comp=1, num_comp=1
    
    
    call amrex_fillcoarsepatch(phi, t_old(lev-1), phi_old(lev-1),  &
                               t_new(lev-1), phi_new(lev-1),  &
                               amrex_geom(lev-1),    fill_physbc,  &
                               amrex_geom(lev  ),    fill_physbc,  &
                               time, src_comp, dst_comp, num_comp, &
                               amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
                               lo_bc, hi_bc)
    
  end subroutine fillcoarsepatch


  ! Apply the boundary conditions (this needs to be passed to
  ! the underlying c++ code, see the arguments of amrex_fillpatch,
  ! and this is why the bind(c) is necessary) 
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
    
    !$omp end parallel

 
  end subroutine fill_physbc 




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
