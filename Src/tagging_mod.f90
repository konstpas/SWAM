module tagging_module

  use iso_c_binding

  use amrex_fort_module

  implicit none

  private

  public :: tag_phi_error

contains

       
  subroutine tag_phi_error (level, time, lo, hi, xlo, dx, surfdist, phi, philo, phihi, tag, taglo, taghi, phierr, &
       settag, cleartag)  
    use domain_module, only : get_surf_pos    
    integer               , intent(in   ) :: level, lo(3), hi(3), philo(4), phihi(4), taglo(4), taghi(4)
    real(amrex_real)      , intent(in   ) :: phi(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
    real(amrex_real)      , intent(in   ) :: xlo(3), dx(3)  ! Surface y-position (to be multifab on DIM-1) 
    real(amrex_real)      , intent(in   ) :: surfdist ! distance from interface to refine grid
    character(kind=c_char), intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(amrex_real)      , intent(in   ) :: time, phierr
    character(kind=c_char), intent(in   ) :: settag, cleartag
    real(amrex_real)                      :: surfpos(lo(1):hi(1))  ! Array surface position 

    integer :: i,j,k

    real(amrex_real) :: phierr_this
    real(amrex_real) :: ydist			! distance from free interface



    call get_surf_pos(xlo, dx, lo, hi, surfpos)
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
            
             ydist = abs(xlo(2) + (j-lo(2))*dx(2) - surfpos(i) ) 
             if (ydist .le. surfdist) then 
                tag(i,j,k) = settag
             endif 
             
          enddo
       enddo
    enddo
    
    
    
  end subroutine tag_phi_error

end module tagging_module

