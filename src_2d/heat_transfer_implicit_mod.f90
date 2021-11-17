module heat_transfer_implicit_module
  
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
  public :: increment_enthalpy_implicit

contains

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy at a new time step
  ! -----------------------------------------------------------------
  ! subroutine increment_enthalpy_implicit(time, lo, hi, &
  !                                        u_old,  uo_lo, uo_hi, &
  !                                        u_new, un_lo, un_hi, &
  !                                        temp_old, to_lo, to_hi, &
  !                                        temp_new, tn_lo, tn_hi , &
  !                                        flxx, fx_lo, fx_hi, &
  !                                        flxy, fy_lo, fy_hi, &
  !                                        idom_old, ido_lo, ido_hi, &
  !                                        idom_new, idn_lo, idn_hi, &
  !                                        geom, dt)


  !   ! ! Re-evaluate domain
  !   ! do i = lo(1)-1,hi(1)+1
  !   !    do  j = lo(2)-1,hi(2)+1

  !   !       ! Points added to the domain
  !   !       if (nint(idom_old(i,j)).eq.0 .and. nint(idom_new(i,j)).ne.0) then
  !   !          u_old(i,j) = u_old(i,j-1)
  !   !       ! Points removed from the domain
  !   !       else if (nint(idom_new(i,j)).eq.0) then
  !   !          u_old(i,j) = 0_amrex_real
  !   !       end if
          
  !   !    end do
  !   ! end do
    
    
  ! end subroutine increment_enthalpy_implicit
  subroutine increment_enthalpy_implicit()
  end subroutine increment_enthalpy_implicit
  

  
end module heat_transfer_implicit_module
