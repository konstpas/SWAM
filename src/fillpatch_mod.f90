module fillpatch_module

  ! -----------------------------------------------------------------
  ! This module is used to fill the multifabs used in the heat
  ! solver (both the points inside the domain and the boundary
  ! conditions are specified). NOTE: This module was adapted from
  ! the AMReX tutorial with very little modification, but the
  ! following functionalities are not fully understood:
  !  1) what is the difference between fillpatch and fillcoarsepatch?
  ! -----------------------------------------------------------------
  
  use iso_c_binding
  use amrex_amr_module

  implicit none

  private

  public :: fillpatch, fillcoarsepatch


contains

  ! Fill phi with data from phi_old and phi_new of current level and one level below.
  subroutine fillpatch (lev, time, phi)

    use amr_data_module, only : t_old, t_new, phi_old, phi_new
    use bc_module, only : lo_bc, hi_bc

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

    use amr_data_module, only : t_old, t_new, phi_old, phi_new
    use bc_module, only : lo_bc, hi_bc

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
    use bc_module, only : lo_bc, hi_bc

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
  

end module fillpatch_module
