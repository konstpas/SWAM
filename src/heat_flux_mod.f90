module heat_flux_module
  
  ! -----------------------------------------------------------------
  ! This module is used to compute the heat flux imposed on the free
  ! surface
  ! -----------------------------------------------------------------
  
  use amrex_amr_module
  
  implicit none

  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: get_boundary_heat_flux

contains
  
  ! -----------------------------------------------------------------
  ! Subroutine used to prescribe the boundary heating on the free
  ! surface. Note that the incoming heat flux is assigned to the
  ! first node under the free surface in the form of a volumetric
  ! heating (after an appropriate dimensionality correction)
  ! -----------------------------------------------------------------   
  subroutine get_boundary_heat_flux(time, xlo, &
                                    dx, lo, hi, &
                                    idom, id_lo, id_hi, &
                                    qb) 

    use read_input_module, only : flux_type

    ! Input and output variables
    integer, intent(in) :: lo(3), hi(3)  
    integer, intent(in) :: id_lo(3), id_hi(3)
    real(amrex_real), intent(in) :: time
    real(amrex_real), intent(in) :: xlo(3)
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2),id_lo(3):id_hi(3))
    real(amrex_real), intent(out) :: qb(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

    if (flux_type.eq.'Gaussian') then 
       call gaussian_heat_flux(time, xlo, dx, lo, hi, &
                               idom, id_lo, id_hi, &
                               qb)
    else
       STOP 'Unknown flux type'
    end if

       
  end subroutine get_boundary_heat_flux
  

  ! -----------------------------------------------------------------
  ! Subroutine used to prescribe a gaussian heat flux active for
  ! time <= time_exposure
  ! -----------------------------------------------------------------   
  subroutine gaussian_heat_flux(time, xlo, &
                                dx, lo, hi, &
                                idom, id_lo, id_hi, &
                                qb) 

    use read_input_module, only : flux_params

    ! Input and output variables
    integer, intent(in) :: lo(3), hi(3)  
    integer, intent(in) :: id_lo(3), id_hi(3)
    real(amrex_real), intent(in) :: time
    real(amrex_real), intent(in) :: xlo(3)
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2),id_lo(3):id_hi(3))
    real(amrex_real), intent(out) :: qb(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

    ! Local variables
    real(amrex_real) :: xpos, zpos
    integer :: i,j,k
    
    ! Prescribe boundary heating
    do   i = lo(1), hi(1) 
       do  j = lo(2), hi(2)
          do k = lo(3), hi(3)
             
             qb(i,j,k) = 0_amrex_real
             
             if (time.ge.flux_params(1) .and. time.lt.flux_params(2)) then

                if(nint(idom(i,j,k)).ne.0 .and. nint(idom(i,j+1,k)).eq.0) then 
                
                   xpos = xlo(1) + (i-lo(1))*dx(1)
                   zpos = xlo(3) + (k-lo(3))*dx(3)
                   
                   ! Note: the term /dx(2) converts a surface heat flux [W/m^2]
                   ! into a volumetric heat flux [W/m^3]
                   qb(i,j,k) = flux_params(3) &
                        *EXP(-((xpos-flux_params(4))**2)/(flux_params(5)**2) &
                             -((zpos-flux_params(6))**2)/(flux_params(7)**2))/dx(2)

             end if
             
          end if

       end do
    end do
    end do
     
  end subroutine gaussian_heat_flux  
  
end module heat_flux_module
