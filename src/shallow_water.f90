module shallow_water_module  

  ! -----------------------------------------------------------------
  ! This module is used to perform all the calculations relative
  ! to the shallow water part of the code.
  ! -----------------------------------------------------------------
  
  use amrex_amr_module
  use amr_data_module, only : surf_ind, &
                              surf_xlo, &
                              surf_dx, &
                              surf_pos, &
                              melt_pos, &
                              melt_vel

  implicit none 

  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: increment_SW
  
  contains 
 
  
  ! -----------------------------------------------------------------
  ! Subroutine used to advance the shallow water equations in time
  ! -----------------------------------------------------------------
  subroutine increment_SW(time, geom, dt)


    type(amrex_geometry), intent(in) :: geom  	! geometry at highest level  
    real(amrex_real), intent(in) :: dt, time		! sub time step, and time 
    real(amrex_real) :: melt_height(surf_ind(1,1):surf_ind(1,2),surf_ind(2,1):surf_ind(2,2))
    real(amrex_real) :: height_flux(surf_ind(1,1):surf_ind(1,2)+1,surf_ind(2,1):surf_ind(2,2)+1,2) 
    !real(amrex_real) :: xfluxp, xfluxm, zfluxp, zfluxm 
    integer :: i,k 
    
 
    ! Momentum continuity equation (not solved for now, only prescribed) 
    melt_vel = 0.05 

    
    ! Mass (column height) continuity equation 
    height_flux = 0. 
  
  
  
    ! Find 'height flux' 
    melt_height = surf_pos-melt_pos 
  
  
    ! X flux 
    do  i = surf_ind(1,1),surf_ind(1,2)+1
       do k = surf_ind(2,1),surf_ind(2,2)
          if (i.eq.surf_ind(1,1)) then ! low boundary             
             if (melt_vel(i,k,1) > 0_amrex_real) then 
                height_flux(i,k,1) = 0 !no influx
             else 
                height_flux(i,k,1) = melt_height(i,k)*melt_vel(i,k,1)
             end if
          elseif (i.eq.surf_ind(1,2)+1) then ! high boundary 
             if (melt_vel(i,k,1) > 0_amrex_real) then 
                height_flux(i,k,1) = melt_height(i-1,k)*melt_vel(i,k,1)
             else 
                height_flux(i,k,1) = 0 !no influx
             end if
          else 	
             if (melt_vel(i,k,1) > 0_amrex_real) then 
                height_flux(i,k,1) = melt_height(i-1,k)*melt_vel(i,k,1)
             else 
                height_flux(i,k,1) = melt_height(i,k)*melt_vel(i,k,1)
             end if
          end if
       end do
    end do
    
    ! Z flux 
    do i = surf_ind(1,1),surf_ind(1,2)
       do k = surf_ind(2,1),surf_ind(2,2)+1
          if     (k.eq.surf_ind(2,1)  ) then 
             if (melt_vel(i,k,2) > 0_amrex_real) then 
                height_flux(i,k,2) = 0! no influx
             else 
                height_flux(i,k,2) = melt_height(i,k)*melt_vel(i,k,2)
             end if
          elseif (k.eq.surf_ind(2,2)+1) then 
             if (melt_vel(i,k,2) > 0_amrex_real) then 
                height_flux(i,k,2) = melt_height(i,k-1)*melt_vel(i,k,2)
             else 
                height_flux(i,k,2) = 0 !no influx  
             end if
          else 
             if (melt_vel(i,k,2) > 0_amrex_real) then 
                height_flux(i,k,2) = melt_height(i,k-1)*melt_vel(i,k,2)
             else 
                height_flux(i,k,2) = melt_height(i,k)*melt_vel(i,k,2)
             end if
          end if
       end do
    end do
  	
		
    do  i = surf_ind(1,1),surf_ind(1,2)
       do k = surf_ind(2,1),surf_ind(2,2)
          surf_pos(i,k) = surf_pos(i,k) & 
               - dt/surf_dx(1) * (height_flux(i+1,k  ,1) - height_flux(i,k,1)) & 
               - dt/surf_dx(2) * (height_flux(i  ,k+1,2) - height_flux(i,k,2))		
       end do
    end do


  end subroutine increment_SW
 

end module shallow_water_module 
