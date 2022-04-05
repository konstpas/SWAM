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
  public :: init_melt_pos
  
contains 

  ! -----------------------------------------------------------------
  ! Subroutine used to advance the shallow water equations in time
  ! for the entire maximum level  
  ! -----------------------------------------------------------------
  subroutine advance_SW()

    use amr_data_module, only : idomain, temp, phi_new, dt
    use read_input_module, only : solve_sw_momentum, solve_heat
    
    ! Local variables
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pid
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: ptemp
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: penth
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx

    if(solve_heat) then
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
   else
      call compute_SW_temperature_terms_decoupled
   end if 

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
  ! Subroutine used to compute the temperature dependent terms
  ! in the shallow water equations
  ! -----------------------------------------------------------------
  subroutine compute_SW_temperature_terms_decoupled()
    
      use amr_data_module, only : surf_temperature, surf_evap_flux, surf_enthalpy, J_th
      use material_properties_module, only : get_enthalpy
      
      ! Local variables
      real(amrex_real) :: enth, temp

      temp = 4000.0
      surf_evap_flux = 0.0
      surf_temperature = temp
      call get_enthalpy(temp, enth)
      surf_enthalpy = enth
      j_th = 2e6;
    
  end subroutine compute_SW_temperature_terms_decoupled
  
  
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
                              !   melt_top, &
                                melt_vel!, &
                              !   surf_pos_grid
    
    
    ! Input and output variables
    real(amrex_real), intent(in) :: dt
    
    ! Local variables
    integer :: i
    real(amrex_real) :: melt_height_old(surf_ind(1,1):surf_ind(1,2))
    real(amrex_real) :: melt_height(surf_ind(1,1):surf_ind(1,2))
    real(amrex_real) :: height_flux(surf_ind(1,1):surf_ind(1,2)+1,1)
   
    ! Initialize the height fluxes
    height_flux = 0. 
    
    ! Compute old column height (input from heat solver)
    melt_height_old = surf_pos - melt_pos 
   !  melt_height_old = surf_pos_grid - melt_pos 

    ! Height fluxes along the x direction
    do  i = surf_ind(1,1),surf_ind(1,2)+1

       ! Boundary condition (low boundary)
       if (i.eq.surf_ind(1,1)) then
          
          if (melt_vel(i,1).gt.0.0_amrex_real) then 
             height_flux(i,1) = 0 !no influx
          else 
             height_flux(i,1) = melt_height_old(i)*melt_vel(i,1)
          end if
          
       ! Boundary condition (high boundary)   
       elseif (i.eq.surf_ind(1,2)+1) then
          
          if (melt_vel(i,1).gt.0.0_amrex_real) then 
             height_flux(i,1) = melt_height_old(i-1)*melt_vel(i,1)
          else 
             height_flux(i,1) = 0 !no influx
          end if
          
       ! Points inside the domain   
       else

          ! No outflow from solid points
          if (melt_height_old(i).eq.0) then
             if (melt_vel(i,1).lt.0) melt_vel(i,1) = 0.0_amrex_real
             if (melt_vel(i+1,1).gt.0) melt_vel(i+1,1) = 0.0_amrex_real
          end if

          ! Compute height fluxes
          if (melt_vel(i,1).gt.0.0_amrex_real) then 
             height_flux(i,1) = melt_height_old(i-1)*melt_vel(i,1)
          else 
             height_flux(i,1) = melt_height_old(i)*melt_vel(i,1)
          end if
          
       end if
    end do

    ! Update the column height equation
    do  i = surf_ind(1,1),surf_ind(1,2)
       
       surf_pos(i) = surf_pos(i) & 
                     - dt/surf_dx(1) * (height_flux(i+1,1) - height_flux(i,1)) &
                     - dt*surf_evap_flux(i)
       
    end do
    melt_height = surf_pos - melt_pos

    ! Update the velocity of newly molten surface elements 
    do  i = surf_ind(1,1),surf_ind(1,2)
       
       if (melt_height_old(i).eq.0.0_amrex_real .and. melt_height(i).gt.0.0_amrex_real) then
      ! if (surf_pos_grid(i).lt.surf_pos(i) .and. melt_pos(i).eq.melt_top(i)) then

          ! It is still possible that both conditions are satisfied, in this case it could
          ! become problematic
          if (melt_vel(i,1).gt.0.0_amrex_real) then
             melt_vel(i+1,1)  = melt_vel(i,1)
          else if (melt_vel(i+1,1).lt.0.0_amrex_real) then
             melt_vel(i,1) = melt_vel(i+1,1)
          end if
          
       end if
       
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
                              !   surf_pos_grid, &
                                surf_temperature, &
                                melt_pos, &
                                melt_vel, &
                                J_th
    
    use read_input_module, only : sw_magnetic, sw_h_cap

    use material_properties_module, only : get_mass_density, &
                                           get_viscosity
    
    ! Input and output variables
    real(amrex_real), intent(in) :: dt
    
    ! Local variables
    integer :: i
    real(amrex_real) :: melt_height(surf_ind(1,1):surf_ind(1,2))
    real(amrex_real) :: adv_term(surf_ind(1,1)+1:surf_ind(1,2))
    real(amrex_real) :: src_term(surf_ind(1,1)+1:surf_ind(1,2))
    real(amrex_real) :: abs_vel
    real(amrex_real) :: hh
    real(amrex_real) :: temp_face
    real(amrex_real) :: visc
    real(amrex_real) :: rho
    real(amrex_real) :: max_vel
    real(amrex_real) :: J_face
    real(amrex_real) :: laplacian_term
    integer :: adv_flag1, adv_flag2
    
    ! Initialize advective and source terms
    adv_term = 0.0_amrex_real
    src_term = 0.0_amrex_real
    
    ! Compute column height (defined on staggered grid)
    melt_height = surf_pos - melt_pos 

    ! Advective term
    do i = surf_ind(1,1)+1,surf_ind(1,2)

       abs_vel = abs(melt_vel(i,1))
       adv_flag1 = 0
       adv_flag2 = 0
       if (melt_vel(i-1,1).ne.0) adv_flag1 = 1
       if (melt_vel(i+1,1).ne.0) adv_flag2 = 1
       adv_term(i) = adv_flag1*(melt_vel(i,1) + abs_vel)/2.0_amrex_real * &
                     (melt_vel(i,1) - melt_vel(i-1,1))/surf_dx(1) + &
                     adv_flag2*(melt_vel(i,1) - abs_vel)/2.0_amrex_real * &
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
          
          ! Calculate the laplacian term only in points inside the melt pool - not on its edges
          if (melt_vel(i-1, 1).ne.0.0 .and. melt_vel(i+1, 1).ne.0.0) then
            laplacian_term = visc * (melt_vel(i+1,1)-2*melt_vel(i,1)+melt_vel(i-1,1))/surf_dx(1)**2
          else
            laplacian_term = 0.0
          end if
          ! Update source term
          src_term(i) =  sw_magnetic*J_face/4 & ! Lorentz force
                         + laplacian_term
          
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
         if (hh.lt.sw_h_cap) hh = sw_h_cap
         melt_vel(i,1) = (melt_vel(i,1) + dt * (src_term(i) - adv_term(i)))/(1+3*visc*dt/(rho*hh**2))
         ! melt_vel(i,1) = (melt_vel(i,1) + dt * (src_term(i) - adv_term(i)))/(1+3*visc*dt/(rho*(hh + sw_h_cap)**2))
         ! melt_vel(i,1) = (melt_vel(i,1) + dt * (src_term(i) - adv_term(i)))
      else
          melt_vel(i,1) = 0.0_amrex_real
      endif
     
      if (ABS(melt_vel(i,1)).gt.ABS(max_vel)) then
         max_vel = melt_vel(i,1)
      end if
      
    end do
    write(*,*) max_vel

    ! Apply boundary conditions
    melt_vel(surf_ind(1,1),1) = melt_vel(surf_ind(1,1)+1,1)
    melt_vel(surf_ind(1,2)+1,1) = melt_vel(surf_ind(1,2),1)

    if (abs(max_vel)*dt.ge.surf_dx(1)) then
      write(*,*) 'CFL not satisfied'
    end if
   
  end subroutine advance_SW_explicit_momentum

  ! -----------------------------------------------------------------
  ! Subroutine used to prescribe a melt pool when the heat response
  ! in not calculated.
  ! -----------------------------------------------------------------
  subroutine init_melt_pos()
    use amr_data_module, only : melt_pos, &
                                surf_pos, &
                                surf_dx, &
                                surf_ind

      ! Local variables
      real(amrex_real) :: x_phys
      integer :: i

      do i = surf_ind(1,1), surf_ind(1,2)
         x_phys = (i+0.5)*surf_dx(1)
         melt_pos(i) = surf_pos(i) - 200e-6 * EXP(-((x_phys-2e-2)**2)/(0.005**2))
         if (surf_pos(i)-melt_pos(i).lt.3e-6) then
            melt_pos(i) = surf_pos(i)
         end if
      end do

  end subroutine init_melt_pos

  
  
end module shallow_water_module
