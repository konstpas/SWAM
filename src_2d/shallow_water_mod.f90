module shallow_water_module  

  ! -----------------------------------------------------------------
  ! This module is used to perform all the calculations relative
  ! to the shallow water part of the code.
  ! -----------------------------------------------------------------
  
  use amrex_amr_module

  implicit none 

  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: advance_SW
  
contains 

  ! -----------------------------------------------------------------
  ! Subroutine used to advance the shallow water equations in time
  ! for the entire maximum level  
  ! -----------------------------------------------------------------
  subroutine advance_SW()

    use amr_data_module, only : idomain, temp, dt
    use read_input_module, only : solve_sw_momentum
    
    ! Local variables
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pid
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptemp
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx

    ! Compute terms that depend on the temperature. Since the temperature
    ! is a multifab, it must be accessed in blocks according to the rules
    ! defined by amrex
    call amrex_mfiter_build(mfi, idomain(amrex_max_level), tiling=.false.)
    do while(mfi%next())

       ! Box
       bx = mfi%validbox()   
       
       ! Pointers
       ptemp => temp(amrex_max_level)%dataptr(mfi)
       pid   => idomain(amrex_max_level)%dataptr(mfi)
       
       ! Terms that depend on the temperature
       call compute_SW_temperature_terms(bx%lo, bx%hi, &
                                         ptemp, lbound(ptemp), ubound(ptemp), &
                                         pid, lbound(pid), ubound(pid))
       
    end do
    call amrex_mfiter_destroy(mfi)    

    ! Advance shallow water equations in time
    if (solve_sw_momentum) then
       STOP "Solver for the shallow water momentum equations is not implemented in the 2D version of the code"
    else
       call advance_SW_fixed_velocity(dt(amrex_max_level))
    end if
    
  end subroutine advance_SW
  

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the temperature dependent terms
  ! in the shallow water equation
  ! -----------------------------------------------------------------
  subroutine compute_SW_temperature_terms(lo, hi, &
                                          temp, temp_lo, temp_hi, &
                                          idom, id_lo, id_hi)
    
    use amr_data_module, only : surf_temperature, surf_evap_flux
    use read_input_module, only : cooling_vaporization
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2) ! Bounds of the current tile box
    integer, intent(in) :: temp_lo(2), temp_hi(2) ! Bounds of temperature box
    integer, intent(in) :: id_lo(2), id_hi(2) ! Bounds of the idomain box
    real(amrex_real), intent(in) :: temp(temp_lo(1):temp_hi(1),temp_lo(2):temp_hi(2))
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2))

    ! Local variables
    integer :: i,j

    ! Loop through the box
    do j = lo(2), hi(2)
       do i = lo(1),hi(1)
          if(nint(idom(i,j)).ne.0 .and. nint(idom(i,j+1)).eq.0) then
             
             ! Evaporation flux
             if (cooling_vaporization) then
                call get_evaporation_flux(temp(i,j), surf_evap_flux(i))
             end if
             
             ! Surface temperature
             surf_temperature(i) = temp(i,j)
             
          end if
       end do
    end do   
    
  end subroutine compute_SW_temperature_terms

  
  ! -----------------------------------------------------------------
  ! Subroutine used to compute the surface erosion due to evaporation
  ! -----------------------------------------------------------------
  subroutine get_evaporation_flux(Ts, evap_flux)

    use material_properties_module, only : get_vapor_pressure, &
                                           get_atomic_mass, &
                                           get_mass_density
    
    ! Input and output variables
    real(amrex_real), intent(in) :: Ts
    real(amrex_real), intent(out) :: evap_flux
    
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
    evap_flux = pv*sqrt(m_A/(2*pi*kb*Ts))/rho_m
    
  end subroutine get_evaporation_flux


  ! -----------------------------------------------------------------
  ! Subroutine used to advance the shallow water equations in time
  ! without solving the momentum equation
  ! -----------------------------------------------------------------
  subroutine advance_SW_fixed_velocity(dt)
    
    use amr_data_module, only : surf_ind, &
                                surf_dx, &
                                surf_evap_flux, &
                                surf_pos, &
                                melt_pos, &
                                melt_vel
    
    use read_input_module, only : meltvel
    
    ! Input and output variables
    real(amrex_real), intent(in) :: dt
    
    ! Local variables
    integer :: i
    real(amrex_real) :: melt_height(surf_ind(1,1):surf_ind(1,2))
    real(amrex_real) :: height_flux(surf_ind(1,1):surf_ind(1,2)+1,1)
    
    ! Momentum continuity equation (not solved for now, only prescribed) 
    melt_vel = meltvel


    ! Mass (column height) continuity equation 
    height_flux = 0. 
    
    ! Find 'height flux' 
    melt_height = surf_pos-melt_pos 
    
    
    ! X flux 
    do  i = surf_ind(1,1),surf_ind(1,2)+1
       
       if (i.eq.surf_ind(1,1)) then ! low boundary             
          if (melt_vel(i,1) > 0_amrex_real) then 
             height_flux(i,1) = 0 !no influx
          else 
             height_flux(i,1) = melt_height(i)*melt_vel(i,1)
          end if
       elseif (i.eq.surf_ind(1,2)+1) then ! high boundary 
          if (melt_vel(i,1) > 0_amrex_real) then 
             height_flux(i,1) = melt_height(i-1)*melt_vel(i,1)
          else 
             height_flux(i,1) = 0 !no influx
          end if
       else
          if (melt_vel(i,1) > 0_amrex_real) then 
             height_flux(i,1) = melt_height(i-1)*melt_vel(i,1)
          else 
             height_flux(i,1) = melt_height(i)*melt_vel(i,1)
          end if
       end if
    end do

    do  i = surf_ind(1,1),surf_ind(1,2)
       surf_pos(i) = surf_pos(i) & 
                     - dt/surf_dx(1) * (height_flux(i+1,1) - height_flux(i,1)) &
                     - dt*surf_evap_flux(i)
    end do
    
  end subroutine advance_SW_fixed_velocity
  
end module shallow_water_module
