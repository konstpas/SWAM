module averagedown_module

  ! -----------------------------------------------------------------
  ! This module is used to perform level averaging , one of the steps
  ! necessary to synchronize the solution on different levels. The
  ! other step necessary for the synchronization, which involves
  ! using flux registers, is performed in the heat equation module 
  ! -----------------------------------------------------------------
  
  use amrex_amr_module
  use amr_data_module, only : phi_new

  implicit none
  
  private
  
  public :: averagedown, averagedownto

contains

  ! ------------------------------------------------------------------
  ! Subroutine to average down to the lowest level
  ! ------------------------------------------------------------------ 
  subroutine averagedown()

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
    
    integer, intent(in) :: clev

    call amrex_average_down(phi_new(clev+1), phi_new(clev), &
                            amrex_geom(clev+1), amrex_geom(clev), &
                            1, 1, amrex_ref_ratio(clev))
    
  end subroutine averagedownto

end module averagedown_module
