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
                              surf_temperature, &
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
  subroutine increment_SW(lo, hi, &
                          temp, temp_lo, temp_hi, &
                          idom, id_lo, id_hi, dt)
    
    !use read_input_module, only : meltvel

    ! Input and output variables
    integer, intent(in) :: lo(3), hi(3) ! Bounds of the current tile box
    integer, intent(in) :: temp_lo(3), temp_hi(3) ! Bounds of temperature box
    integer, intent(in) :: id_lo(3), id_hi(3) ! Bounds of the idomain box
    real(amrex_real), intent(in) :: dt
    real(amrex_real), intent(in) :: temp(temp_lo(1):temp_hi(1),temp_lo(2):temp_hi(2),temp_lo(3):temp_hi(3))
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2),id_lo(3):id_hi(3))

    ! Local variables
    real(amrex_real) :: xdot_vap
    integer :: i,j,k    

    ! Update melt thickness (only evaporation erosion)
    do j = lo(2), hi(2)
       do i = lo(1),hi(1)
          do k = lo(3),hi(3)
             if(nint(idom(i,j,k)).ne.0 .and. nint(idom(i,j+1,k)).eq.0) then
                call get_evaporation_flux(temp(i,j,k), xdot_vap)
                surf_pos(i,k) = surf_pos(i,k) - xdot_vap*dt
                surf_temperature(i,k) = temp(i,j,k)
             end if
          end do
       end do
    end do   
    
    ! ! Momentum continuity equation (not solved for now, only prescribed) 
    ! melt_vel = meltvel

    
    ! ! Mass (column height) continuity equation 
    ! height_flux = 0. 
  
    ! ! Find 'height flux' 
    ! melt_height = surf_pos-melt_pos 
  
  
    ! ! X flux 
    ! do  i = surf_ind(1,1),surf_ind(1,2)+1
    !    do k = surf_ind(2,1),surf_ind(2,2)
    !       if (i.eq.surf_ind(1,1)) then ! low boundary             
    !          if (melt_vel(i,k,1) > 0_amrex_real) then 
    !             height_flux(i,k,1) = 0 !no influx
    !          else 
    !             height_flux(i,k,1) = melt_height(i,k)*melt_vel(i,k,1)
    !          end if
    !       elseif (i.eq.surf_ind(1,2)+1) then ! high boundary 
    !          if (melt_vel(i,k,1) > 0_amrex_real) then 
    !             height_flux(i,k,1) = melt_height(i-1,k)*melt_vel(i,k,1)
    !          else 
    !             height_flux(i,k,1) = 0 !no influx
    !          end if
    !       else
    !          if (melt_vel(i,k,1) > 0_amrex_real) then 
    !             height_flux(i,k,1) = melt_height(i-1,k)*melt_vel(i,k,1)
    !          else 
    !             height_flux(i,k,1) = melt_height(i,k)*melt_vel(i,k,1)
    !          end if
    !       end if
    !    end do
    ! end do
    
    ! ! Z flux 
    ! do i = surf_ind(1,1),surf_ind(1,2)
    !    do k = surf_ind(2,1),surf_ind(2,2)+1
    !       if     (k.eq.surf_ind(2,1)  ) then 
    !          if (melt_vel(i,k,2) > 0_amrex_real) then 
    !             height_flux(i,k,2) = 0! no influx
    !          else 
    !             height_flux(i,k,2) = melt_height(i,k)*melt_vel(i,k,2)
    !          end if
    !       elseif (k.eq.surf_ind(2,2)+1) then 
    !          if (melt_vel(i,k,2) > 0_amrex_real) then 
    !             height_flux(i,k,2) = melt_height(i,k-1)*melt_vel(i,k,2)
    !          else 
    !             height_flux(i,k,2) = 0 !no influx  
    !          end if
    !       else 
    !          if (melt_vel(i,k,2) > 0_amrex_real) then 
    !             height_flux(i,k,2) = melt_height(i,k-1)*melt_vel(i,k,2)
    !          else 
    !             height_flux(i,k,2) = melt_height(i,k)*melt_vel(i,k,2)
    !          end if
    !       end if
    !    end do
    ! end do
  	
		
    ! do  i = surf_ind(1,1),surf_ind(1,2)
    !    do k = surf_ind(2,1),surf_ind(2,2)
    !       surf_pos(i,k) = surf_pos(i,k) & 
    !            - dt/surf_dx(1) * (height_flux(i+1,k  ,1) - height_flux(i,k,1)) &
    !            - dt/surf_dx(2) * (height_flux(i  ,k+1,2) - height_flux(i,k,2))
    !    end do
    ! end do


  end subroutine increment_SW
 

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the surface erosion due to evaporation
  ! -----------------------------------------------------------------
  subroutine get_evaporation_flux(Ts, xdot_vap)

    use material_properties_module, only : get_vapor_pressure, &
                                          get_atomic_mass, &
                                          get_mass_density
    
    ! Input and output variables
    real(amrex_real), intent(in) :: Ts
    real(amrex_real), intent(out) :: xdot_vap
    
    ! Local variables
    real(amrex_real) :: kb = 1.38064852E-23 ! Boltzmann constant [m^2*kg/(s^2*K)]
    real(amrex_real) :: pv ! Vapor pressure
    real(amrex_real) :: m_A ! Atomic mass [g/mol]
    real(amrex_real) :: pi = 3.1415927
    real(amrex_real) :: Na = 6.02214076E23 ! Avogadro's number [mol^-1]  
    real(amrex_real) :: rho_m
    
    call get_atomic_mass(m_A)
    call get_vapor_pressure(Ts, pv)
    call get_mass_density(Ts, rho_m)
    
    ! Conversion from g/mol to kg
    m_A = m_A*1E-3/Na
    
    ! Evaporation flux
    xdot_vap = pv*sqrt(m_A/(2*pi*kb*Ts))/rho_m
    
  end subroutine get_evaporation_flux
  
end module shallow_water_module 
