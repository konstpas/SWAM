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

    use amr_data_module, only : idomain, temp, phi_new, dt
    use read_input_module, only : solve_sw_momentum
    
    ! Local variables
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pid
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptemp
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: penth
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
       penth => phi_new(amrex_max_level)%dataptr(mfi)
       pid   => idomain(amrex_max_level)%dataptr(mfi)
       
       ! Terms that depend on the temperature
       call compute_SW_temperature_terms(bx%lo, bx%hi, &
                                         ptemp, lbound(ptemp), ubound(ptemp), &
                                         pid, lbound(pid), ubound(pid), &
                                         penth, lbound(penth), ubound(penth))
       
    end do
    call amrex_mfiter_destroy(mfi)    

    ! Advance shallow water equations in time
    if (solve_sw_momentum) then
       call advance_SW_explicit(dt(amrex_max_level))
    else
       call advance_SW_fixed_velocity(dt(amrex_max_level))
    end if
    
  end subroutine advance_SW
  

  ! -----------------------------------------------------------------
  ! Subroutine used to advance the shallow water equations in time
  ! without solving the momentum equation
  ! -----------------------------------------------------------------
  subroutine advance_SW_fixed_velocity(dt)

    use amr_data_module, only : melt_vel
    use read_input_module, only : fixed_melt_velocity
    
    ! Input and output variables
    real(amrex_real), intent(in) :: dt

    ! Advance the column height equation
    call advance_SW_explicit_height(dt)

    ! Momentum continuity equation (not solved for now, only prescribed) 
    melt_vel = fixed_melt_velocity

  end subroutine advance_SW_fixed_velocity


  ! -----------------------------------------------------------------
  ! Subroutine used to advance the shallow water equations in time
  ! with the explicit solver
  ! -----------------------------------------------------------------
  subroutine advance_SW_explicit(dt)
    
    ! Input and output variables
    real(amrex_real), intent(in) :: dt

    ! Advance the column height equation
    call advance_SW_explicit_height(dt)

    ! Advance the momentum equation 
    call advance_SW_explicit_momentum(dt)

  end subroutine advance_SW_explicit

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the temperature dependent terms
  ! in the shallow water equations
  ! -----------------------------------------------------------------
  subroutine compute_SW_temperature_terms(lo, hi, &
                                          temp, temp_lo, temp_hi, &
                                          idom, id_lo, id_hi, &
                                          enth, enth_lo, enth_hi)
    
    use amr_data_module, only : surf_temperature, surf_evap_flux, surf_enthalpy
    use read_input_module, only : cooling_vaporization
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2) ! Bounds of the current tile box
    integer, intent(in) :: temp_lo(2), temp_hi(2) ! Bounds of temperature box
    integer, intent(in) :: id_lo(2), id_hi(2) ! Bounds of the idomain box
    integer, intent(in) :: enth_lo(2), enth_hi(2) ! Bounds of enthalpy box
    real(amrex_real), intent(in) :: temp(temp_lo(1):temp_hi(1),temp_lo(2):temp_hi(2))
    real(amrex_real), intent(in) :: enth(enth_lo(1):enth_hi(1),enth_lo(2):enth_hi(2))
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

             ! Surface enthalpy
             surf_enthalpy(i) = enth(i,j)
             
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
  ! Subroutine used to advance the shallow water equation for the
  ! column height with an explicit upwind scheme
  ! -----------------------------------------------------------------
  subroutine advance_SW_explicit_height(dt)
    
    use amr_data_module, only : surf_ind, &
                                surf_dx, &
                                surf_evap_flux, &
                                surf_pos, &
                                melt_pos, &
                                melt_vel
    
    
    ! Input and output variables
    real(amrex_real), intent(in) :: dt
    
    ! Local variables
    integer :: i
    real(amrex_real) :: melt_height(surf_ind(1,1):surf_ind(1,2))
    real(amrex_real) :: height_flux(surf_ind(1,1):surf_ind(1,2)+1,1)
   
    ! Initialize the height fluxes
    height_flux = 0. 
    
    ! Compute old column height (input from heat solver)
    melt_height = surf_pos - melt_pos 
    
    ! Height fluxes along the x direction
    do  i = surf_ind(1,1),surf_ind(1,2)+1

       ! Boundary condition (low boundary)
       if (i.eq.surf_ind(1,1)) then
          
          if (melt_vel(i,1).gt.0.0_amrex_real) then 
             height_flux(i,1) = 0 !no influx
          else 
             height_flux(i,1) = melt_height(i)*melt_vel(i,1)
          end if
          
       ! Boundary condition (high boundary)   
       elseif (i.eq.surf_ind(1,2)+1) then
          
          if (melt_vel(i,1).gt.0.0_amrex_real) then 
             height_flux(i,1) = melt_height(i-1)*melt_vel(i,1)
          else 
             height_flux(i,1) = 0 !no influx
          end if
          
       ! Points inside the domain   
       else
          
          if (melt_vel(i,1).gt.0.0_amrex_real) then 
             height_flux(i,1) = melt_height(i-1)*melt_vel(i,1)
          else 
             height_flux(i,1) = melt_height(i)*melt_vel(i,1)
          end if
          
       end if
    end do

    ! Update the column height equation
    do  i = surf_ind(1,1),surf_ind(1,2)
       
       surf_pos(i) = surf_pos(i) & 
                     - dt/surf_dx(1) * (height_flux(i+1,1) - height_flux(i,1)) &
                     - dt*surf_evap_flux(i)
       
    end do
    
  end subroutine advance_SW_explicit_height


  ! -----------------------------------------------------------------
  ! Subroutine used to advance the shallow water equation for the
  ! momentum with an explicit scheme
  ! -----------------------------------------------------------------
  subroutine advance_SW_explicit_momentum(dt)
    
    use amr_data_module, only : surf_ind, &
                                surf_dx, &
                                surf_pos, &
                                surf_temperature, &
                                melt_pos, &
                                melt_vel, &
                                J_th
    
    use read_input_module, only : sw_magnetic

    use material_properties_module, only : get_mass_density, &
                                           get_viscosity
    
    ! Input and output variables
    real(amrex_real), intent(in) :: dt
    
    ! Local variables
    integer :: i
    real(amrex_real) :: melt_height(surf_ind(1,1):surf_ind(1,2))
    real(amrex_real) :: adv_term(surf_ind(1,1)+1:surf_ind(1,2))
    real(amrex_real) :: src_term(surf_ind(1,1)+1:surf_ind(1,2))
    real(amrex_real) :: abs_vel_l, abs_vel_r
    real(amrex_real) :: vel_l_half, vel_r_half
    real(amrex_real) :: abs_vel
    real(amrex_real) :: hh
    real(amrex_real) :: temp_face
    real(amrex_real) :: visc
    real(amrex_real) :: rho
    real(amrex_real) :: max_vel
    real(amrex_real) :: J_face
    
    ! Initialize advective and source terms
    adv_term = 0.0_amrex_real
    src_term = 0.0_amrex_real
    
    ! Compute column height (defined on staggered grid)
    melt_height = surf_pos - melt_pos 

    ! Advective term
    do i = surf_ind(1,1)+1,surf_ind(1,2)

       ! vel_l_half = (melt_vel(i-1,1)+melt_vel(i,1))/2.0
       ! vel_r_half = (melt_vel(i+1,1)+melt_vel(i,1))/2.0
       ! adv_term(i) = (vel_l_half + ABS(vel_l_half))/2.0_amrex_real * &
       !               (melt_vel(i,1) - melt_vel(i-1,1))/surf_dx(1) + &
       !               (vel_r_half - ABS(vel_r_half))/2.0_amrex_real * &
       !               (melt_vel(i+1,1) - melt_vel(i,1))/surf_dx(1)

       abs_vel = abs(melt_vel(i,1))
       adv_term(i) = (melt_vel(i,1) + abs_vel)/2.0_amrex_real * &
                     (melt_vel(i,1) - melt_vel(i-1,1))/surf_dx(1) + &
                     (melt_vel(i,1) - abs_vel)/2.0_amrex_real * &
                     (melt_vel(i+1,1) - melt_vel(i,1))/surf_dx(1)
       
    end do

    ! Source term
    do i = surf_ind(1,1)+1,surf_ind(1,2)

       ! Melt thickness
       hh = (melt_height(i) + melt_height(i-1))/2.0_amrex_real
       temp_face = (surf_temperature(i) + surf_temperature(i-1))/2.0_amrex_real
       
       ! Compute source terms only for grid points with a finite melt thickness
       if (hh.gt.0.0_amrex_real) then

          ! Material properties
          call get_viscosity(temp_face,visc)
          call get_mass_density(temp_face,rho)
          
          J_face = (J_th(i-1)+J_th(i))/2
          !  J_face = J_face*surf_dx(1)**2 ! Dimension correction
          ! Update source term
          src_term(i) =  sw_magnetic*J_face/4 & ! Lorentz force
                         + visc * (melt_vel(i+1,1)-2*melt_vel(i,1)+melt_vel(i-1,1))/surf_dx(1)**2
          
          ! Fix dimensionality
          src_term(i) = src_term(i)/rho

       end if
       
    end do

    max_vel = 0.0
    ! Update momentum equation
    do  i = surf_ind(1,1)+1,surf_ind(1,2)
      hh = (melt_height(i) + melt_height(i-1))/2.0_amrex_real
      temp_face = (surf_temperature(i) + surf_temperature(i-1))/2.0_amrex_real
      call get_viscosity(temp_face,visc)
      call get_mass_density(temp_face,rho)
      if (hh.gt.0.0_amrex_real) then
         !melt_vel(i,1) = (melt_vel(i,1) + dt * (src_term(i) - adv_term(i)))/(1+3*visc*dt/(rho*hh**2))
         melt_vel(i,1) = (melt_vel(i,1) + dt * (src_term(i) - adv_term(i)))
      else
          melt_vel(i,1) = 0.0_amrex_real
      endif
       if (ABS(melt_vel(i,1)).gt.ABS(max_vel)) then
         max_vel = melt_vel(i,1)
       end if
    end do

    ! Apply boundary conditions
    melt_vel(surf_ind(1,1),1) = melt_vel(surf_ind(1,1)+1,1)
    melt_vel(surf_ind(1,2)+1,1) = melt_vel(surf_ind(1,2),1)

    if (max_vel*dt.ge.surf_dx(1)) then
      write(*,*) 'CFL not satisfied'
    end if
  end subroutine advance_SW_explicit_momentum
  
  
end module shallow_water_module
