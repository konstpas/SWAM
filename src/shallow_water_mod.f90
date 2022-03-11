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
     use read_input_module, only : solve_sw_momentum, sw_solver
     
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
         if (sw_solver.eq.'geoclaw') then
            call advance_SW_geoclaw(dt(amrex_max_level))
         else if (sw_solver.eq.'explicit') then
            call advance_SW_explicit_height(dt(amrex_max_level))
            call advance_SW_explicit_momentum(dt(amrex_max_level))
         else
            STOP 'Unknown shallow water solver'
         end if

     else
        call advance_SW_fixed_velocity(dt(amrex_max_level))
     end if
     
   end subroutine advance_SW
   
 
   ! -----------------------------------------------------------------
   ! Subroutine used to compute the temperatuer dependent terms in
   ! the shallow water equations
   ! -----------------------------------------------------------------
   subroutine compute_SW_temperature_terms(lo, hi, &
                                           temp, temp_lo, temp_hi, &
                                           idom, id_lo, id_hi, &
                                           enth, enth_lo, enth_hi)
     
     use amr_data_module, only : surf_temperature, surf_evap_flux, surf_enthalpy
     use read_input_module, only : cooling_vaporization
     
     ! Input and output variables
     integer, intent(in) :: lo(3), hi(3) ! Bounds of the current tile box
     integer, intent(in) :: temp_lo(3), temp_hi(3) ! Bounds of temperature box
     integer, intent(in) :: enth_lo(3), enth_hi(3) ! Bounds of enthalpy box
     integer, intent(in) :: id_lo(3), id_hi(3) ! Bounds of the idomain box
     real(amrex_real), intent(in) :: temp(temp_lo(1):temp_hi(1),temp_lo(2):temp_hi(2),temp_lo(3):temp_hi(3))
     real(amrex_real), intent(in) :: enth(enth_lo(1):enth_hi(1),enth_lo(2):enth_hi(2),enth_lo(3):enth_hi(3))
     real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2),id_lo(3):id_hi(3))
 
     ! Local variables
     integer :: i,j,k
 
     ! Update melt thickness (only evaporation erosion)
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1),hi(1)
              if(nint(idom(i,j,k)).ne.0 .and. nint(idom(i,j+1,k)).eq.0) then

                 ! Evaporation flux
                 if (cooling_vaporization) then
                    call get_evaporation_flux(temp(i,j,k), surf_evap_flux(i,k))
                 end if

                 ! Surface temperature
                 surf_temperature(i,k) = temp(i,j,k)
                 surf_enthalpy(i,k) = enth(i,j,k)
              end if
           end do
        end do 
     end do  
     
   end subroutine compute_SW_temperature_terms
 
   
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
   

   ! -----------------------------------------------------------------
   ! Subroutine used to advance the shallow water equations in time
   ! with a prescribed velocity
   ! -----------------------------------------------------------------
   subroutine advance_SW_fixed_velocity(dt)
     
     use amr_data_module, only : surf_ind, &
                                 surf_dx, &
                                 surf_evap_flux, &
                                 surf_pos, &
                                 melt_pos, &
                                 melt_vel
 
     use read_input_module, only : fixed_melt_velocity
     
     ! Input and output variables
     real(amrex_real), intent(in) :: dt
     
     ! Local variables
     integer :: i,k
     real(amrex_real) :: melt_height(surf_ind(1,1):surf_ind(1,2),surf_ind(2,1):surf_ind(2,2))
     real(amrex_real) :: height_flux(surf_ind(1,1):surf_ind(1,2)+1,surf_ind(2,1):surf_ind(2,2)+1,2) 
         
     ! Momentum continuity equation (not solved, only prescribed) 
     melt_vel(surf_ind(1,1):surf_ind(1,2)+1,surf_ind(2,1):surf_ind(2,2)+1,1) = fixed_melt_velocity(1)
     melt_vel(surf_ind(1,1):surf_ind(1,2)+1,surf_ind(2,1):surf_ind(2,2)+1,2) = fixed_melt_velocity(2)
     
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
                           - dt/surf_dx(1) * (height_flux(i+1,k,1) - height_flux(i,k,1)) &
                           - dt/surf_dx(2) * (height_flux(i,k+1,2) - height_flux(i,k,2)) &
                           - dt*surf_evap_flux(i,k)
        end do
     end do
     
   end subroutine advance_SW_fixed_velocity
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
      
      use read_input_module, only : sw_Bx, sw_Bz, sw_h_cap
  
      use material_properties_module, only : get_mass_density, &
                                             get_viscosity
      
      ! Input and output variables
      real(amrex_real), intent(in) :: dt
      
      ! Local variables
      integer :: i, j
      real(amrex_real) :: melt_height(surf_ind(1,1):surf_ind(1,2),surf_ind(2,1):surf_ind(2,2))
      real(amrex_real) :: adv_term_x(surf_ind(1,1):surf_ind(1,2),surf_ind(2,1):surf_ind(2,2))
      real(amrex_real) :: adv_term_z(surf_ind(1,1):surf_ind(1,2),surf_ind(2,1):surf_ind(2,2))
      real(amrex_real) :: src_term_x(surf_ind(1,1):surf_ind(1,2),surf_ind(2,1):surf_ind(2,2))
      real(amrex_real) :: src_term_z(surf_ind(1,1):surf_ind(1,2),surf_ind(2,1):surf_ind(2,2))
      real(amrex_real) :: abs_vel_x, abs_vel_z
      real(amrex_real) :: hh
      real(amrex_real) :: temp_face
      real(amrex_real) :: visc
      real(amrex_real) :: rho
      real(amrex_real) :: max_vel_x, max_vel_z
      real(amrex_real) :: J_face
      real(amrex_real) :: laplacian_term
      
      ! Initialize advective and source terms
      adv_term_x = 0.0_amrex_real
      adv_term_z = 0.0_amrex_real
      src_term_x = 0.0_amrex_real
      src_term_z = 0.0_amrex_real
      
      ! Compute column height (defined on staggered grid)
      melt_height = surf_pos - melt_pos 
  
      ! X component of advective acceleration
      do i = surf_ind(1,1)+1,surf_ind(1,2)
         do j = surf_ind(2,1),surf_ind(2,2)
  
            abs_vel_x = abs(melt_vel(i,j,1))
            abs_vel_z = abs(melt_vel(i,j,2)) 

            if (melt_vel(i-1,j,1).ne.0) then
               adv_term_x(i,j) = (melt_vel(i,j,1) + abs_vel_x)/2.0_amrex_real * &
                                 (melt_vel(i,j,1) - melt_vel(i-1,j,1))/surf_dx(1)
            end if
            if (melt_vel(i+1,j,1).ne.0) then
               adv_term_x(i,j) = adv_term_x(i,j) + (melt_vel(i,j,1) - abs_vel_x)/2.0_amrex_real * &
                                 (melt_vel(i+1,j,1) - melt_vel(i,j,1))/surf_dx(1) 
            end if
            if (melt_vel(i,j-1,1).ne.0 .and. j.gt.surf_ind(2,1)) then
               adv_term_x(i,j) = adv_term_x(i,j) + (melt_vel(i,j,2) + abs_vel_z)/2.0_amrex_real * &
                                 (melt_vel(i,j,1) - melt_vel(i,j-1,1))/surf_dx(2)
            end if
            if (melt_vel(i,j+1,1).ne.0 .and. j.lt.surf_ind(2,2)) then
               adv_term_x(i,j) = adv_term_x(i,j) + (melt_vel(i,j,2) - abs_vel_z)/2.0_amrex_real * &
                                 (melt_vel(i,j+1,1) - melt_vel(i,j,1))/surf_dx(2)
            end if   
         end do
      end do

      ! Z component of advective acceleration
      do i = surf_ind(1,1),surf_ind(1,2)
         do j = surf_ind(2,1)+1,surf_ind(2,2)
            
            abs_vel_x = abs(melt_vel(i,j,1))
            abs_vel_z = abs(melt_vel(i,j,2))

            if (melt_vel(i-1,j,2).ne.0 .and. i.gt.surf_ind(1,1)) then
               adv_term_z(i,j) = (melt_vel(i,j,1) + abs_vel_x)/2.0_amrex_real * &
                                 (melt_vel(i,j,2) - melt_vel(i-1,j,2))/surf_dx(1)
            end if
            if (melt_vel(i+1,j,2).ne.0 .and. i.lt.surf_ind(1,2)) then
               adv_term_z(i,j) = adv_term_z(i,j) + (melt_vel(i,j,1) - abs_vel_x)/2.0_amrex_real * &
                                 (melt_vel(i+1,j,2) - melt_vel(i,j,2))/surf_dx(1)
            end if
            if (melt_vel(i,j-1,2).ne.0) then
               adv_term_z(i,j) = adv_term_z(i,j) + (melt_vel(i,j,2) + abs_vel_z)/2.0_amrex_real * &
                                 (melt_vel(i,j,2) - melt_vel(i,j-1,2))/surf_dx(2)
            end if
            if (melt_vel(i,j+1,2).ne.0) then 
               adv_term_z(i,j) = adv_term_z(i,j) + (melt_vel(i,j,2) - abs_vel_z)/2.0_amrex_real * &
                                 (melt_vel(i,j+1,2) - melt_vel(i,j,2))/surf_dx(2) 
            end if
         end do
      end do
  
      ! Source term in the x direction
      do i = surf_ind(1,1)+1, surf_ind(1,2)
         do j = surf_ind(2,1), surf_ind(2,2)
  
            ! Melt thickness - interpolation in x direction
            hh = (melt_height(i,j) + melt_height(i-1,j))/2.0_amrex_real
            temp_face = (surf_temperature(i,j) + surf_temperature(i-1,j))/2.0_amrex_real
            
            ! Compute source terms only for grid points with a finite melt thickness
            if (hh.gt.0.0_amrex_real) then
   
               ! Material properties
               call get_viscosity(temp_face,visc) ! Should you use different viscosity for direvatives in the z direction
               call get_mass_density(temp_face,rho)
               
               J_face = (J_th(i-1,j)+J_th(i,j))/2
               
               ! Calculate the laplacian term only in points inside the melt pool - not on its edges
               ! Partial update from second derivative in direction of x
               laplacian_term = 0.0
               if (melt_vel(i-1, j, 1).ne.0.0 .and. melt_vel(i+1, j, 1).ne.0.0) then
                  laplacian_term = visc * (melt_vel(i+1,j,1)-2*melt_vel(i,j,1)+melt_vel(i-1,j,1)/surf_dx(1)**2)
               end if
               ! Partial update from second derivative in direction of z
               if (j.gt.surf_ind(2,1) .and. j.lt.surf_ind(2,2) .and. melt_vel(i,j-1,1).ne.0.0 .and. melt_vel(i,j+1,1).ne.0.0) then
                  laplacian_term  = laplacian_term + visc*(melt_vel(i,j+1,1)-2*melt_vel(i,j,1)+melt_vel(i,j-1,1)/surf_dx(2)**2)
               end if
               ! Update source term for accelerationin the x direction
               src_term_x(i,j) =  sw_Bz*J_face/4 & ! Lorentz force
                                 + laplacian_term
               ! Fix dimensionality
               src_term_x(i,j) = src_term_x(i,j)/rho
            end if
         end do
      end do

      ! Source term in the z direction (it assumed that jxb is allways alligned with the x axis)
      do i = surf_ind(1,1),surf_ind(1,2)
         do j = surf_ind(2,1)+1,surf_ind(2,2)
            ! Melt thickness - interpolation in z direction
            hh = (melt_height(i,j) + melt_height(i,j-1))/2.0_amrex_real
            temp_face = (surf_temperature(i,j) + surf_temperature(i,j-1))/2.0_amrex_real
            
            ! Compute source terms only for grid points with a finite melt thickness
            if (hh.gt.0.0_amrex_real) then
   
               ! Material properties
               call get_viscosity(temp_face,visc)
               call get_mass_density(temp_face,rho)
               
               J_face = (J_th(i,j)+J_th(i,j-1))/2
               
               ! Calculate the laplacian term only in points inside the melt pool - not on its edges
               ! Partial update from second derivative in direction of x
               laplacian_term = 0.0
               if (melt_vel(i-1,j,2).ne.0.0 .and. melt_vel(i+1,j,2).ne.0.0 .and. i.gt.surf_ind(1,1) .and. i.lt.surf_ind(1,2)) then
                  laplacian_term = visc * (melt_vel(i+1,j,2)-2*melt_vel(i,j,2)+melt_vel(i-1,j,2))/surf_dx(1)**2
               end if
               ! Partial update from second derivative in direction of z
               if (melt_vel(i,j-1,2).ne.0.0 .and. melt_vel(i,j+1,2).ne.0.0) then
                  laplacian_term  = laplacian_term + visc*(melt_vel(i,j+1,2)-2*melt_vel(i,j,2)+melt_vel(i,j-1,2))/surf_dx(2)**2
               end if
               ! Update source term for accelerationin the z direction
               src_term_z(i,j) =  -sw_Bx*J_face/4 & ! Lorentz force
                                 + laplacian_term
               
               ! Fix dimensionality
               src_term_z(i,j) = src_term_z(i,j)/rho
   
            end if
         
         end do
      end do
  
      max_vel_x = 0.0
      max_vel_z = 0.0
      ! Update momentum equation
      do  i = surf_ind(1,1),surf_ind(1,2)
         do  j = surf_ind(2,1),surf_ind(2,2)
         
            ! Update of x component of velocity
            if (i.gt.surf_ind(1,1)) then
               hh = (melt_height(i,j) + melt_height(i-1,j))/2.0_amrex_real
               temp_face = (surf_temperature(i,j) + surf_temperature(i-1,j))/2.0_amrex_real
               
               call get_viscosity(temp_face,visc)
               call get_mass_density(temp_face,rho)
               
               if (hh.gt.0.0_amrex_real) then
                  if (hh.lt.sw_h_cap) hh = sw_h_cap
                  melt_vel(i,j,1) = (melt_vel(i,j,1) + dt * (src_term_x(i,j) - adv_term_x(i,j)))/(1+3*visc*dt/(rho*hh**2))
               else
                  melt_vel(i,j,1) = 0.0_amrex_real
               endif

               if (ABS(melt_vel(i,j,1)).gt.ABS(max_vel_x)) then
                  max_vel_x = melt_vel(i,j,1)
               end if
            ! No penetration boundary condition
            else
               melt_vel(i,j,1) = 0
            end if

            ! Update of z component of velocity
            if (j.gt.surf_ind(2,1)) then
               hh = (melt_height(i,j) + melt_height(i,j-1))/2.0_amrex_real
               temp_face = (surf_temperature(i,j) + surf_temperature(i,j-1))/2.0_amrex_real
               
               call get_viscosity(temp_face,visc)
               call get_mass_density(temp_face,rho)
               
               if (hh.gt.0.0_amrex_real) then
                  if (hh.lt.sw_h_cap) hh = sw_h_cap
                  melt_vel(i,j,2) = (melt_vel(i,j,2) + dt * (src_term_z(i,j) - adv_term_z(i,j)))/(1+3*visc*dt/(rho*hh**2))
               else
                  melt_vel(i,j,1) = 0.0_amrex_real
               endif
               
               if (ABS(melt_vel(i,j,2)).gt.ABS(max_vel_z)) then
                  max_vel_z = melt_vel(i,j,2)
               end if
            ! No penetration boundary condition
            else
               melt_vel(i,j,2) = 0
            end if
         
         end do
      end do
   !  write(*,*) max_vel_x
   !  write(*,*) max_vel_z

      ! Apply boundary conditions - no penetration at the edge
      do j = surf_ind(2,1), surf_ind(2,2)
         melt_vel(surf_ind(1,2)+1,j,1) = 0
      end do
      do i = surf_ind(1,1), surf_ind(1,2)
         melt_vel(i,surf_ind(2,2)+1,2) = 0
      end do

      if (abs(max_vel_x)*dt.ge.surf_dx(1)) then
         STOP 'CFL not satisfied in x direction'
      end if
      if (abs(max_vel_z)*dt.ge.surf_dx(2)) then
         STOP 'CFL not satisfied in z direction'
      end if
    end subroutine advance_SW_explicit_momentum


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
      integer :: i, j
      real(amrex_real) :: melt_height_old(surf_ind(1,1):surf_ind(1,2),surf_ind(2,1):surf_ind(2,2))
      real(amrex_real) :: melt_height(surf_ind(1,1):surf_ind(1,2),surf_ind(2,1):surf_ind(2,2))
      real(amrex_real) :: height_flux(surf_ind(1,1):surf_ind(1,2)+1,surf_ind(2,1):surf_ind(2,2)+1,2)
     
      ! Initialize the height fluxes
      height_flux = 0. 
      
      ! Compute old column height (input from heat solver)
      melt_height_old = surf_pos - melt_pos 
  
      ! Height fluxes along the x direction
      do  i = surf_ind(1,1),surf_ind(1,2)+1
         do j = surf_ind(2,1),surf_ind(2,2)
  
            ! Boundary condition (low boundary)
            if (i.eq.surf_ind(1,1)) then
               
               if (melt_vel(i,j,1).gt.0.0_amrex_real) then 
                  height_flux(i,j,1) = 0 !no influx
               else 
                  height_flux(i,j,1) = melt_height_old(i,j)*melt_vel(i,j,1)
               end if
               
            ! Boundary condition (high boundary)   
            elseif (i.eq.surf_ind(1,2)+1) then
               
               if (melt_vel(i,j,1).gt.0.0_amrex_real) then 
                  height_flux(i,j,1) = melt_height_old(i-1,j)*melt_vel(i,j,1)
               else 
                  height_flux(i,j,1) = 0 !no influx
               end if
               
            ! Points inside the domain   
            else
   
               ! No outflow from solid points
               if (melt_height_old(i,j).eq.0) then
                  if (melt_vel(i,j,1).lt.0) melt_vel(i,j,1) = 0.0_amrex_real
                  if (melt_vel(i+1,j,1).gt.0) melt_vel(i+1,j,1) = 0.0_amrex_real
               end if
   
               ! Compute height fluxes
               if (melt_vel(i,j,1).gt.0.0_amrex_real) then 
                  height_flux(i,j,1) = melt_height_old(i-1,j)*melt_vel(i,j,1)
               else 
                  height_flux(i,j,1) = melt_height_old(i,j)*melt_vel(i,j,1)
               end if
               
            end if
         end do
      end do

      ! Height fluxes along the z direction
      do  i = surf_ind(1,1),surf_ind(1,2)
         do j = surf_ind(2,1),surf_ind(2,2)+1
  
            ! Boundary condition (low boundary)
            if (j.eq.surf_ind(2,1)) then
               
               if (melt_vel(i,j,2).gt.0.0_amrex_real) then 
                  height_flux(i,j,2) = 0 !no influx
               else 
                  height_flux(i,j,2) = melt_height_old(i,j)*melt_vel(i,j,2)
               end if
               
            ! Boundary condition (high boundary)   
            elseif (j.eq.surf_ind(2,2)+1) then
               
               if (melt_vel(i,j,2).gt.0.0_amrex_real) then 
                  height_flux(i,j,2) = melt_height_old(i,j-1)*melt_vel(i,j,2)
               else 
                  height_flux(i,j,2) = 0 !no influx
               end if
               
            ! Points inside the domain   
            else
   
               ! No outflow from solid points
               if (melt_height_old(i,j).eq.0) then
                  if (melt_vel(i,j,2).lt.0) melt_vel(i,j,2) = 0.0_amrex_real
                  if (melt_vel(i,j+1,2).gt.0) melt_vel(i,j+1,2) = 0.0_amrex_real
               end if
   
               ! Compute height fluxes
               if (melt_vel(i,j,2).gt.0.0_amrex_real) then 
                  height_flux(i,j,2) = melt_height_old(i,j-1)*melt_vel(i,j,2)
               else 
                  height_flux(i,j,2) = melt_height_old(i,j)*melt_vel(i,j,2)
               end if
               
            end if
         end do
      end do
  
      ! Update the column height equation
      do  i = surf_ind(1,1),surf_ind(1,2)
         do j = surf_ind(2,1), surf_ind(2,2)
         
            surf_pos(i,j) = surf_pos(i,j) & 
                        - dt/surf_dx(1) * (height_flux(i+1,j,1) - height_flux(i,j,1)) &
                        - dt/surf_dx(2) * (height_flux(i,j+1,2) - height_flux(i,j,2)) &
                        - dt*surf_evap_flux(i,j)

         end do         
      end do
      melt_height = surf_pos - melt_pos
  
      ! Update the velocity of newly molten surface elements 
      do  i = surf_ind(1,1),surf_ind(1,2)
         do  j = surf_ind(2,1),surf_ind(2,2)
         
            ! Fluid parcles upwind in the x-direction will be used to initiliaze x-component of
            ! velocity, while upwind parcels in the z-direction are used to initialize the z component.
            ! If a newly added column has no upwind fluid in the x or z direction the corresponding velocity
            ! is initialized to zero.
            if (melt_height_old(i,j).eq.0.0_amrex_real .and. melt_height(i,j).gt.0.0_amrex_real) then 

               ! It is still possible that both conditions are satisfied, in this case it could
               ! become problematic.
               ! x component of velocity
               if (melt_vel(i,j,1).gt.0.0_amrex_real) then
                  melt_vel(i+1,j,1)  = melt_vel(i,j,1)
               else if (melt_vel(i+1,j,1).lt.0.0_amrex_real) then
                  melt_vel(i,j,1) = melt_vel(i+1,j,1)
               end if
               ! z component of velocity
               if (melt_vel(i,j,2).gt.0.0_amrex_real) then
                  melt_vel(i,j+1,2)  = melt_vel(i,j,2)
               else if (melt_vel(i,j+1,2).lt.0.0_amrex_real) then
                  melt_vel(i,j,2) = melt_vel(i,j+1,2)
               end if

            end if
         end do
      end do
      
    end subroutine advance_SW_explicit_height

   ! -----------------------------------------------------------------
   ! Subroutine used to advance the shallow water equations adapted
   ! from geoclaw
   ! -----------------------------------------------------------------
   subroutine advance_SW_geoclaw(dt)

     use amr_data_module, only : surf_ind, &
                                 surf_dx, &
                                 surf_evap_flux, &
                                 surf_pos, &
                                 surf_temperature, &
                                 melt_pos, &
                                 melt_vel, &
                                 qnew
     
     use material_properties_module, only : get_viscosity, get_mass_density
     use read_input_module, only : cooling_vaporization, &
                                   sw_iter, &
                                   sw_gravity, &
                                   sw_drytol, &
                                   sw_jxb

     ! Input and output variables
     real(amrex_real), intent(in) :: dt
     
     ! Local variables
     integer :: i,j,m,mw,ibc,jbc
     real(amrex_real) :: hL,hR,huL,huR,hvL,hvR,bL,bR,pL,pR,uL,vL,phiL,uR,vR,phiR,sL,sR,uhat,chat,sRoe1,sRoe2,sE1,sE2
     real(amrex_real) :: wall(1:3), sw(1:3), fw(1:3,1:3)
     real(amrex_real) :: hstar,s1m,s2m,rare1,rare2,hstartest
     real(amrex_real) , allocatable :: q1d(:,:),s(:,:),qadd(:,:),amdq(:,:),apdq(:,:),fadd(:,:),aux1d(:,:),fwave(:,:,:)
     real(amrex_real) , dimension(1:3,surf_ind(1,1)-3:surf_ind(1,2)+3,surf_ind(2,1)-3:surf_ind(2,2)+3) :: qold
     real(amrex_real) , dimension(surf_ind(1,1)-3:surf_ind(1,2)+3,surf_ind(2,1)-3:surf_ind(2,2)+3) :: hold,uold,vold,aux
     real(amrex_real) , dimension(1:3,surf_ind(1,1):surf_ind(1,2),surf_ind(2,1):surf_ind(2,2)) :: Srce
     real(amrex_real) :: dtdx,dtdy,dx,dy
     real(amrex_real) :: visc,rho

     ! Define grid and time-step
     dx = surf_dx(1)
     dy = surf_dx(2)
     dtdx = dt/surf_dx(1) 
     dtdy = dt/surf_dx(2)
     
     ! Bathymetry
     aux(surf_ind(1,1):surf_ind(1,2),surf_ind(2,1):surf_ind(2,2)) = melt_pos
         
     ! Fill qnew with results from heat solver
     do j = surf_ind(2,1),surf_ind(2,2)
        do i = surf_ind(1,1),surf_ind(1,2)
          qnew(1,i,j) = surf_pos(i,j) - melt_pos(i,j) 
          if (qnew(1,i,j).gt.sw_drytol) then
             qnew(2,i,j) = melt_vel(i,j,1)*qnew(1,i,j)
             qnew(3,i,j) = melt_vel(i,j,2)*qnew(1,i,j)
          else
             qnew(2,i,j) = 0.0_amrex_real
             qnew(3,i,j) = 0.0_amrex_real
          end if
       end do
    end do
     
     ! Apply boundary conditions 
     ! solid wall (assumes 2'nd component is velocity or momentum in x): left
     do j = surf_ind(2,1) - 3,surf_ind(2,2) + 3
        do ibc = 1,3
           do m=1,3
              qnew(m,surf_ind(1,1)-ibc,j) = qnew(m,surf_ind(1,1)+ibc,j)
           end do
        end do
     end do
     do j = surf_ind(2,1) - 3,surf_ind(2,2) + 3
        do ibc = 1,3
           qnew(2,surf_ind(1,1)-ibc,j) =-qnew(2,surf_ind(1,1)+ibc,j)
           aux(surf_ind(1,1)-ibc,j) = aux(surf_ind(1,1),j)
        end do
     end do
     ! solid wall (assumes 2'nd component is velocity or momentum in x): right
     do j = surf_ind(2,1) - 3,surf_ind(2,2) + 3
        do ibc = 1,3
           do m=1,3
              qnew(m,surf_ind(1,2)+ibc,j) = qnew(m,surf_ind(1,2)-ibc,j)
           end do
        end do
     end do
     do j = surf_ind(2,1) - 3,surf_ind(2,2) + 3
        do ibc = 1,3
           qnew(2,surf_ind(1,2)+ibc,j) =-qnew(2,surf_ind(1,2)-ibc,j)
           aux(surf_ind(1,2)+ibc,j) = aux(surf_ind(1,2),j)
        end do
     end do
     ! solid wall (assumes 3'rd component is velocity or momentum in y): bottom
     do jbc = 1,3
        do i = surf_ind(1,1) - 3,surf_ind(1,2) + 3
           do m=1,3
              qnew(m,i,surf_ind(2,1)-jbc) = qnew(m,i,surf_ind(2,1)+jbc)
           end do
        end do
     end do
     do jbc = 1,3
        do i = surf_ind(1,1) - 3,surf_ind(1,2) + 3
           qnew(3,i,surf_ind(2,1)-jbc) =-qnew(3,i,surf_ind(2,1)+jbc)
           aux(i,surf_ind(2,1)-jbc) = aux(i,surf_ind(2,1))
        end do
     end do
     ! solid wall (assumes 3'rd component is velocity or momentum in y): top
     do jbc = 1,3
        do i = surf_ind(1,1) - 3,surf_ind(1,2) + 3
           do m=1,3
              qnew(m,i,surf_ind(2,2)+jbc) = qnew(m,i,surf_ind(2,2)-jbc)
           end do
        end do
     end do
     do jbc = 1,3
        do i = surf_ind(1,1) - 3,surf_ind(1,2) + 3
           qnew(3,i,surf_ind(2,2)+jbc) =-qnew(3,i,surf_ind(2,2)-jbc)
           aux(i,surf_ind(2,2)+jbc) = aux(i,surf_ind(2,2))
        end do
     end do
     
     ! Initialize solution (fill old solution with new values)
     qold = qnew
     
     do j = surf_ind(2,1) - 3,surf_ind(2,2) + 3
        do i = surf_ind(1,1) - 3,surf_ind(1,2) + 3       
           if (qnew(1,i,j) .gt. sw_drytol) then
              hold(i,j) = qnew(1,i,j)
              uold(i,j) = qnew(2,i,j) / qnew(1,i,j)
              vold(i,j) = qnew(3,i,j) / qnew(1,i,j)
           else
              hold(i,j) = 0.  
              uold(i,j) = 0.  
              vold(i,j) = 0.  
           end if
        end do
     end do
     
     ! Source terms
     
     Srce = 0.
     
     do j = surf_ind(2,1),surf_ind(2,2)
        do i = surf_ind(1,1),surf_ind(1,2)
           
           if ( hold(i,j) .gt. sw_drytol) then
              call get_viscosity(surf_temperature(i,j),visc)
              call get_mass_density(surf_temperature(i,j),rho)
              ! Source terms for the continuity equation
              if (cooling_vaporization) then
                 Srce(1,i,j) = -surf_evap_flux(i,j)
              end if
              ! Source terms for the momentum equation along x
              Srce(2,i,j) = sw_jxb(1) &
                            - 3.*visc*uold(i,j)/(hold(i,j)**2) &
                            + visc * ( (uold(i+1,j) + uold(i-1,j) - 2.*uold(i,j))/(dx**2) & 
                            + (uold(i,j+1) + uold(i,j-1) - 2.*uold(i,j))/(dy**2))
              Srce(2,i,j) = Srce(2,i,j)/rho
              ! Source terms for the momentum equation along z
              Srce(3,i,j) = sw_jxb(2) &
                            - 3.*visc*vold(i,j)/(hold(i,j)**2) &
                            + visc * ( (vold(i+1,j) + vold(i-1,j) - 2.*vold(i,j))/(dx**2) & 
                            + (vold(i,j+1) + vold(i,j-1) - 2.*vold(i,j))/(dy**2))
              Srce(3,i,j) = Srce(3,i,j)/rho
           else
              Srce(1,i,j) = 0. 
              Srce(2,i,j) = 0. 
              Srce(3,i,j) = 0. 
           end if
           
        end do
     end do
     
     ! Update with fluxes along x direction
     
     ! Allocate
     allocate(  q1d(1:3    ,surf_ind(1,1)-3:surf_ind(1,2)+3)) 
     allocate(aux1d(1:3    ,surf_ind(1,1)-3:surf_ind(1,2)+3)) 
     allocate( qadd(1:3    ,surf_ind(1,1)-3:surf_ind(1,2)+3)) 
     allocate( amdq(1:3    ,surf_ind(1,1)-3:surf_ind(1,2)+3)) 
     allocate( apdq(1:3    ,surf_ind(1,1)-3:surf_ind(1,2)+3)) 
     allocate( fadd(1:3    ,surf_ind(1,1)-3:surf_ind(1,2)+3)) 
     allocate(    s(1:3    ,surf_ind(1,1)-3:surf_ind(1,2)+3)) 
     allocate(fwave(1:3,1:3,surf_ind(1,1)-3:surf_ind(1,2)+3)) 
     
     
     ! Flux2 for x sweeps
     do j = surf_ind(2,1) - 3,surf_ind(2,2) + 3
        
        do i = surf_ind(1,1) - 3,surf_ind(1,2) + 3
           do m=1,3 
              q1d(m,i) = qold(m,i,j)
              aux1d(1,i) = aux(i,j)
           end do
        end do
        
        
        !rpn2
        do i = surf_ind(1,1) - 2,surf_ind(1,2) + 3
           
           do m=1,3
              if (q1d(1,i-1).lt.0.) then
                 q1d(m,i-1) = 0.
              end if
              if (q1d(1,i).lt.0.) then
                 q1d(m,i) = 0.
              end if
           end do
           
           do mw=1,3
              s(mw,i)=0.
              fwave(1,mw,i)=0.
              fwave(2,mw,i)=0.
              fwave(3,mw,i)=0.
           end do
           
           
           if ( (q1d(1,i-1) .lt. sw_drytol) .and. (q1d(1,i) .lt. sw_drytol) ) goto 30

           hL = q1d(1,i-1)
           hR = q1d(1,i)
           huL = q1d(2,i-1)
           huR = q1d(2,i)
           hvL = q1d(3,i-1)
           hvR = q1d(3,i)
           bL = aux1d(1,i-1)
           bR = aux1d(1,i)
           pL = 0.
           pR = 0.
           
           ! check for wet/dry boundary
           if (hR .gt. sw_drytol) then
              uR=huR/hR
              vR=hvR/hR
              phiR = 0.5*sw_gravity*hR**2 + huR**2/hR
           else
              hR = 0.
              huR = 0.
              hvR = 0.
              uR = 0.
              vR = 0.
              phiR = 0.
           end if
           
           if (hL .gt. sw_drytol) then
              uL=huL/hL
              vL=hvL/hL
              phiL = 0.5*sw_gravity*hL**2 + huL**2/hL
           else
              hL=0.
              huL=0.
              hvL=0.
              uL=0.
              vL=0.
              phiL = 0.
           end if
           
           wall(1:3)=1.
           
           if (hR .lt. sw_drytol) then
              call riemanntype(hL,hL,uL,-uL,1,sw_drytol,sw_gravity,hstar,s1m,s2m,rare1,rare2)
              hstartest=max(hL,hstar)
              if ( (hstartest+bL) .lt. bR) then  !right state should become ghost values that mirror left for wall problem
                 wall(2)=0.
                 wall(3)=0.
                 hR=hL
                 huR=-huL
                 bR=bL
                 phiR=phiL
                 uR=-uL
                 vR=vL
              elseif ( (hL+bL) .lt. bR) then
                 bR=hL+bL
              end if
              
           elseif (hL .lt. sw_drytol) then  ! right surface is lower than left topo
              call riemanntype(hR,hR,-uR,uR,1,sw_drytol,sw_gravity,hstar,s1m,s2m,rare1,rare2)
              hstartest=max(hR,hstar)
              if ( (hstartest+bR) .lt. bL) then !left state should become ghost values that mirror right
                 wall(1)=0.
                 wall(2)=0.
                 hL=hR
                 huL=-huR
                 bL=bR
                 phiL=phiR
                 uL=-uR
                 vL=vR
              elseif ( (hR+bR) .lt. bL) then
                 bL=hR+bR
              end if
              
           end if
           
           !determine wave speeds
           sL = uL - sqrt(sw_gravity*hL) ! 1 wave speed of left state
           sR = uR + sqrt(sw_gravity*hR) ! 2 wave speed of right state
           
           uhat=(sqrt(sw_gravity*hL)*uL + sqrt(sw_gravity*hR)*uR)/(sqrt(sw_gravity*hR)+sqrt(sw_gravity*hL)) ! Roe average
           chat=sqrt(sw_gravity*0.5*(hR+hL)) ! Roe average
           sRoe1=uhat-chat ! Roe wave speed 1 wave
           sRoe2=uhat+chat ! Roe wave speed 2 wave
           
           sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
           sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave
           
           call riemann_aug_JCP(sw_iter,3,3,hL,hR,huL,huR,hvL,hvR,bL,bR,uL,uR,vL,&
                vR,phiL,phiR,pL,pR,sE1,sE2,sw_drytol,sw_gravity,1.d0,sw,fw)
           
           ! eliminate ghost fluxes for wall
           do mw=1,3
              sw(mw)=sw(mw)*wall(mw)
              fw(1,mw)=fw(1,mw)*wall(mw)
              fw(2,mw)=fw(2,mw)*wall(mw)
              fw(3,mw)=fw(3,mw)*wall(mw)
           end do
           
           do mw=1,3
              s(mw,i)=sw(mw)
              fwave(1,mw,i)=fw(1,mw)
              fwave(2,mw,i)=fw(2,mw)
              fwave(3,mw,i)=fw(3,mw)
           end do
           
30         continue
        end do ! i
        
        amdq = 0.
        apdq = 0.
        do i=surf_ind(1,1)-2,surf_ind(1,2)+3 
           do  mw=1,3
              if (s(mw,i) .lt. 0.) then
                 amdq(1:3,i) = amdq(1:3,i) + fwave(1:3,mw,i)
              elseif (s(mw,i) .gt. 0.) then
                 apdq(1:3,i)  = apdq(1:3,i) + fwave(1:3,mw,i)
              else
                 amdq(1:3,i) = amdq(1:3,i) + 0.5*fwave(1:3,mw,i)
                 apdq(1:3,i) = apdq(1:3,i) + 0.5*fwave(1:3,mw,i)
              end if
           end do
        end do
        
        
        qadd = 0.
        do i = surf_ind(1,1),surf_ind(1,2) 
           do m = 1,3
              qadd(m,i-1) = qadd(m,i-1) - dtdx*amdq(m,i)
           end do
           do m = 1,3
              qadd(m,i) = qadd(m,i) - dtdx*apdq(m,i)
           end do
        end do
        
        
        ! second order correction terms , with a flux limiter as specified by mthlim.
        fadd = 0. 
        
        call fwave_limiter(surf_ind(1,1),surf_ind(1,2),fwave,s,dtdx,fadd) ! optional
        
        do i=surf_ind(1,1),surf_ind(1,2) 
           do m=1,3
              qnew(m,i,j) = qnew(m,i,j) + qadd(m,i) - dtdx*(fadd(m,i+1) - fadd(m,i))
           end do
        end do
        
     end do
     
     deallocate(  q1d)
     deallocate(aux1d)
     deallocate( qadd)
     deallocate( amdq)
     deallocate( apdq)
     deallocate( fadd)
     deallocate(    s)
     deallocate(fwave)
     
     ! Update with fluxes along y direction
     
     allocate(  q1d(1:3    ,surf_ind(2,1)-3:surf_ind(2,2)+3))
     allocate(aux1d(1:3    ,surf_ind(2,1)-3:surf_ind(2,2)+3))
     allocate( qadd(1:3    ,surf_ind(2,1)-3:surf_ind(2,2)+3))
     allocate( amdq(1:3    ,surf_ind(2,1)-3:surf_ind(2,2)+3))
     allocate( apdq(1:3    ,surf_ind(2,1)-3:surf_ind(2,2)+3))
     allocate( fadd(1:3    ,surf_ind(2,1)-3:surf_ind(2,2)+3))
     allocate(    s(1:3    ,surf_ind(2,1)-3:surf_ind(2,2)+3))
     allocate(fwave(1:3,1:3,surf_ind(2,1)-3:surf_ind(2,2)+3))
     
     
     do j = surf_ind(1,1) - 3,surf_ind(1,2) + 3
        
        do i = surf_ind(2,1) - 3,surf_ind(2,2) + 3
           do m=1,3
              q1d(m,i) = qold(m,j,i)
              aux1d(1,i) = aux(j,i)
           end do
        end do
        
        
        do i = surf_ind(2,1) - 2,surf_ind(2,2) + 3
           
           do m=1,3
              if (q1d(1,i-1).lt.0.) then
                 q1d(m,i-1) = 0.
              end if
              if (q1d(1,i).lt.0.) then
                 q1d(m,i) = 0.
              end if
           end do
           
           do mw=1,3
              s(mw,i)=0.
              fwave(1,mw,i)=0.
              fwave(2,mw,i)=0.
              fwave(3,mw,i)=0.
           end do
           
           
           if ( (q1d(1,i-1) .lt. sw_drytol) .and. (q1d(1,i) .lt. sw_drytol) ) goto 40
           
           hL = q1d(1,i-1)
           hR = q1d(1,i)
           huL = q1d(3,i-1)
           huR = q1d(3,i)
           hvL = q1d(2,i-1)
           hvR = q1d(2,i)
           bL = aux1d(1,i-1)
           bR = aux1d(1,i)
           pL = 0.
           pR = 0.
           
           ! check for wet/dry boundary
           if (hR .gt. sw_drytol) then
              uR=huR/hR
              vR=hvR/hR
              phiR = 0.5*sw_gravity*hR**2 + huR**2/hR
           else
              hR = 0.
              huR = 0.
              hvR = 0.
              uR = 0.
              vR = 0.
              phiR = 0.
           end if
           
           if (hL .gt. sw_drytol) then
              uL=huL/hL
              vL=hvL/hL
              phiL = 0.5*sw_gravity*hL**2 + huL**2/hL
           else
              hL=0.
              huL=0.
              hvL=0.
              uL=0.
              vL=0.
              phiL = 0.
           end if
           
           
           wall(1:3)=1.
           
           if (hR .lt. sw_drytol) then
              call riemanntype(hL,hL,uL,-uL,1,sw_drytol,sw_gravity,hstar,s1m,s2m,rare1,rare2)
              hstartest=max(hL,hstar)
              if ( (hstartest+bL) .lt. bR) then  !right state should become ghost values that mirror left for wall problem
                 wall(2)=0.
                 wall(3)=0.
                 hR=hL
                 huR=-huL
                 bR=bL
                 phiR=phiL
                 uR=-uL
                 vR=vL
              elseif ( (hL+bL) .lt. bR) then
                 bR=hL+bL
              end if
              
              
           elseif (hL .lt. sw_drytol) then  ! right surface is lower than left topo
              call riemanntype(hR,hR,-uR,uR,1,sw_drytol,sw_gravity,hstar,s1m,s2m,rare1,rare2)
              hstartest=max(hR,hstar)
              if ( (hstartest+bR) .lt. bL) then !left state should become ghost values that mirror right
                 wall(1)=0.
                 wall(2)=0.
                 hL=hR
                 huL=-huR
                 bL=bR
                 phiL=phiR
                 uL=-uR
                 vL=vR
              elseif ( (hR+bR) .lt. bL) then
                 bL=hR+bR
              end if
           end if
           
           !determine wave speeds
           sL = uL - sqrt(sw_gravity*hL) ! 1 wave speed of left state
           sR = uR + sqrt(sw_gravity*hR) ! 2 wave speed of right state
           
           uhat=(sqrt(sw_gravity*hL)*uL + sqrt(sw_gravity*hR)*uR)/(sqrt(sw_gravity*hR)+sqrt(sw_gravity*hL)) ! Roe average
           chat=sqrt(sw_gravity*0.5*(hR+hL)) ! Roe average
           sRoe1=uhat-chat ! Roe wave speed 1 wave
           sRoe2=uhat+chat ! Roe wave speed 2 wave
           
           sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
           sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave
           
           call riemann_aug_JCP(sw_iter,3,3,hL,hR,huL,huR,hvL,hvR,bL,bR,uL,uR,vL,&
                vR,phiL,phiR,pL,pR,sE1,sE2,sw_drytol,sw_gravity,1.d0,sw,fw)
          
           ! eliminate ghost fluxes for wall
           do mw=1,3
              sw(mw)=sw(mw)*wall(mw)
              fw(1,mw)=fw(1,mw)*wall(mw)
              fw(2,mw)=fw(2,mw)*wall(mw)
              fw(3,mw)=fw(3,mw)*wall(mw)
           end do
           
           do mw=1,3
              s(mw,i)=sw(mw)
              fwave(1,mw,i)=fw(1,mw)
              fwave(2,mw,i)=fw(2,mw)
              fwave(3,mw,i)=fw(3,mw)
           end do
           
40         continue
        end do ! i
        
        amdq = 0.
        apdq = 0.
        do i=surf_ind(2,1)-2,surf_ind(2,2)+3
           do  mw=1,3
              if (s(mw,i) .lt. 0.) then
                 amdq(1:3,i) = amdq(1:3,i) + fwave(1:3,mw,i)
              elseif (s(mw,i) .gt. 0.) then
                 apdq(1:3,i)  = apdq(1:3,i) + fwave(1:3,mw,i)
              else
                 amdq(1:3,i) = amdq(1:3,i) + 0.5*fwave(1:3,mw,i)
                 apdq(1:3,i) = apdq(1:3,i) + 0.5*fwave(1:3,mw,i)
              end if
           end do
        end do
        
        qadd = 0.
        do i = surf_ind(2,1),surf_ind(2,2)
           do m = 1,3
              qadd(m,i-1) = qadd(m,i-1) - dtdy*amdq(m,i)
           end do
           do m = 1,3
              qadd(m,i) = qadd(m,i) - dtdy*apdq(m,i)
           end do
       end do
       
       
       ! second order correction terms , with a flux limiter as specified by mthlim.
       fadd = 0.
       
       call fwave_limiter(surf_ind(2,1),surf_ind(2,2),fwave,s,dtdy,fadd) ! optional
       
       do i=surf_ind(2,1),surf_ind(2,2)
          do m=1,3
             qnew(m,j,i) = qnew(m,j,i) + qadd(m,i) - dtdy*(fadd(m,i+1) - fadd(m,i))
          end do
       end do
       
    end do
    
    ! Update with source terms
    do j=surf_ind(2,1),surf_ind(2,2)
       do i=surf_ind(1,1),surf_ind(1,2)
          
          qnew(1,i,j) = qnew(1,i,j) + dt*(Srce(1,i,j))
          qnew(2,i,j) = qnew(2,i,j) + dt*(Srce(1,i,j)*uold(i,j) + Srce(2,i,j)*hold(i,j))
          qnew(3,i,j) = qnew(3,i,j) + dt*(Srce(1,i,j)*vold(i,j) + Srce(3,i,j)*hold(i,j))
          
       end do
    end do
    
    
    ! Apply boundary conditions to updated solution before passing results to heat solver
    
    ! solid wall (assumes 2'nd component is velocity or momentum in x): left
    do j = surf_ind(2,1) - 3,surf_ind(2,2) + 3
       do ibc = 1,3
          do m=1,3
             qnew(m,surf_ind(1,1)-ibc,j) = qnew(m,surf_ind(1,1)+ibc,j)
          end do
       end do
    end do
    do j = surf_ind(2,1) - 3,surf_ind(2,2) + 3
       do ibc = 1,3
          qnew(2,surf_ind(1,1)-ibc,j) =-qnew(2,surf_ind(1,1)+ibc,j)
          aux(surf_ind(1,1)-ibc,j) = aux(surf_ind(1,1),j)
       end do
    end do
    ! solid wall (assumes 2'nd component is velocity or momentum in x): right
    do j = surf_ind(2,1) - 3,surf_ind(2,2) + 3
       do ibc = 1,3
          do m=1,3
             qnew(m,surf_ind(1,2)+ibc,j) = qnew(m,surf_ind(1,2)-ibc,j)
          end do
       end do
    end do
    do j = surf_ind(2,1) - 3,surf_ind(2,2) + 3
       do ibc = 1,3
          qnew(2,surf_ind(1,2)+ibc,j) =-qnew(2,surf_ind(1,2)-ibc,j)
          aux(surf_ind(1,2)+ibc,j) = aux(surf_ind(1,2),j)
       end do
    end do
    ! solid wall (assumes 3'rd component is velocity or momentum in y): bottom
    do jbc = 1,3
       do i = surf_ind(1,1) - 3,surf_ind(1,2) + 3
          do m=1,3
             qnew(m,i,surf_ind(2,1)-jbc) = qnew(m,i,surf_ind(2,1)+jbc)
          end do
       end do
    end do
    do jbc = 1,3
       do i = surf_ind(1,1) - 3,surf_ind(1,2) + 3
          qnew(3,i,surf_ind(2,1)-jbc) =-qnew(3,i,surf_ind(2,1)+jbc)
          aux(i,surf_ind(2,1)-jbc) = aux(i,surf_ind(2,1))
       end do
    end do
    ! solid wall (assumes 3'rd component is velocity or momentum in y): top
    do jbc = 1,3
       do i = surf_ind(1,1) - 3,surf_ind(1,2) + 3
          do m=1,3
             qnew(m,i,surf_ind(2,2)+jbc) = qnew(m,i,surf_ind(2,2)-jbc)
          end do
       end do
    end do
    do jbc = 1,3
       do i = surf_ind(1,1) - 3,surf_ind(1,2) + 3
          qnew(3,i,surf_ind(2,2)+jbc) =-qnew(3,i,surf_ind(2,2)-jbc)
          aux(i,surf_ind(2,2)+jbc) = aux(i,surf_ind(2,2))
       end do
    end do

    ! Update surface position
    do j=surf_ind(2,1),surf_ind(2,2)
       do i=surf_ind(1,1),surf_ind(1,2)
         !  if (qnew(1,i,j).gt.sw_drytol) then
             surf_pos(i,j) = melt_pos(i,j) + qnew(1,i,j)
         !  end if
       end do
    end do

    ! Update melt velocity along x
    do j=surf_ind(2,1),surf_ind(2,2)
       do i=surf_ind(1,1),surf_ind(1,2)+1
          if (qnew(1,i,j).gt.sw_drytol) then
             melt_vel(i,j,1) = qnew(2,i,j)/qnew(1,i,j)
          end if
       end do
    end do

    ! Update melt velocity along z
    do j=surf_ind(2,1),surf_ind(2,2)+1
       do i=surf_ind(1,1),surf_ind(1,2)
          if (qnew(1,i,j).gt.sw_drytol) then
             melt_vel(i,j,2) = qnew(3,i,j)/qnew(1,i,j)
          end if
       end do
    end do
    
    
  end subroutine advance_SW_geoclaw


  subroutine riemanntype(hL,hR,uL,uR,maxiter,drytol,g,hm,s1m,s2m,rare1,rare2)
    real(amrex_real), intent(in) :: hL,hR,uL,uR,drytol,g
    integer, intent(in) :: maxiter
    real(amrex_real), intent(out) :: hm,s1m,s2m,rare1,rare2
    integer :: iter
    real(amrex_real) :: h_min,h_max,delu,F_min,F_max,um,h0,gL,gR,F0,dfdh,slope,u1m,u2m
    
    
    ! Riemann structure
    h_min=min(hR,hL)
    h_max=max(hR,hL)
    delu=uR-uL
    if (h_min .lt. drytol) then
       hm=0.
       um=0.
       s1m=uR+uL-2.*sqrt(g*hR)+2.*sqrt(g*hL)
       s2m=uR+uL-2.*sqrt(g*hR)+2.*sqrt(g*hL)
       if (hL .lt. 0.) then
          rare2=1.
          rare1=0.
       else
          rare1=1.
          rare2=0.
       end if
    else
       F_min= delu+2.*(sqrt(g*h_min)-sqrt(g*h_max))
       F_max= delu+(h_max-h_min)*(sqrt(0.5*g*(h_max+h_min)/(h_max*h_min)))
       if (F_min .gt. 0) then !2-rarefactions
          hm=(1./(16.*g))*max(0.,-delu+2.*(sqrt(g*hL)+sqrt(g*hR)))**2
          um=sign(1.d0,hm)*(uL+2.*(sqrt(g*hL)-sqrt(g*hm)))
          s1m=uL+2.*sqrt(g*hL)-3.*sqrt(g*hm)
          s2m=uR-2.*sqrt(g*hR)+3.*sqrt(g*hm)
          rare1=1.
          rare2=1.
       elseif (F_max.le.0) then !2 shocks
          !root finding using a Newton iteration on sqrt(h)===
          h0=h_max
          do iter=1,maxiter
             gL=sqrt(0.5*g*(1./h0 + 1./hL))
             gR=sqrt(0.5*g*(1./h0 + 1./hR))
             F0=delu+(h0-hL)*gL + (h0-hR)*gR
             dfdh=gL-g*(h0-hL)/(4.*(h0**2)*gL)+gR-g*(h0-hR)/(4.*(h0**2)*gR)
             slope=2.*sqrt(h0)*dfdh
             h0=(sqrt(h0)-F0/slope)**2
          end do
          hm=h0
          u1m=uL-(hm-hL)*sqrt((0.5*g)*(1./hm + 1./hL))
          u2m=uR+(hm-hR)*sqrt((0.5*g)*(1./hm + 1./hR))
          um=0.5*(u1m+u2m)
          s1m=u1m-sqrt(g*hm)
          s2m=u2m+sqrt(g*hm)
          rare1=0.
          rare2=0.
       else! %one shock one rarefaction
          h0=h_min
          do iter=1,maxiter
             F0=delu+2.*(sqrt(g*h0)-sqrt(g*h_max))+(h0-h_min)*sqrt(0.5*g*(1./h0+1./h_min))
             slope=(F_max-F0)/(h_max-h_min)
             h0=h0-F0/slope
          end do
          hm=h0
          if (hL .gt. hR) then
             um=uL+2.*sqrt(g*hL)-2.*sqrt(g*hm)
             s1m=uL+2.*sqrt(g*hL)-3.*sqrt(g*hm)
             s2m=uL+2.*sqrt(g*hL)-sqrt(g*hm)
             rare1=1.
             rare2=0.
          else
             s2m=uR-2.*sqrt(g*hR)+3.*sqrt(g*hm)
             s1m=uR-2.*sqrt(g*hR)+sqrt(g*hm)
             um=uR-2.*sqrt(g*hR)+2.*sqrt(g*hm)
             rare2=1.
             rare1=0.
          end if
       end if
    end if
    
 end subroutine riemanntype





 subroutine riemann_aug_JCP(maxiter,meqn,mwaves,hL,hR,huL,huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2,drytol,g,rho,sw,fw)

   integer, intent(in) :: maxiter,meqn,mwaves
   real(amrex_real), intent(in) :: hL,hR,huL,huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,drytol,g,rho
   real(amrex_real), intent(inout) :: sE1,sE2
   real(amrex_real), intent(out) :: sw(1:mwaves),fw(1:meqn,1:mwaves)
   real(amrex_real) :: A(1:3,1:3),r(1:3,1:3),lambda(1:3),del(1:3),beta(1:3)
   real(amrex_real) :: delh,delhu,delphi,delb,delnorm,rare1st,rare2st,sdelta,raremin,raremax,criticaltol,convergencetol,raretol
   real(amrex_real) :: criticaltol_2, hustar_interface,s1s2bar,s1s2tilde,hbar,hLstar,hRstar,s1m,s2m,hm
   real(amrex_real) :: huRstar,huLstar,uRstar,uLstar,hstarHLL,deldelh,deldelphi,delP,det1,det2,det3,determinant
   real(amrex_real) :: rare1,rare2,rarecorrector,rarecorrectortest,sonic
   integer :: mw,k,iter
   

   
   !determine del vectors
   delh = hR-hL
   delhu = huR-huL
   delphi = phiR-phiL
   delb = bR-bL
   delp = pR - pL
   delnorm = delh**2 + delphi**2
   
   call riemanntype(hL,hR,uL,uR,1,drytol,g,hm,s1m,s2m,rare1,rare2)
   
   lambda(1)= min(sE1,s2m)  ! Modified Einfeldt speed
   lambda(3)= max(sE2,s1m)  ! Modified Eindfeldt speed
   sE1=lambda(1)
   sE2=lambda(3)
   lambda(2) = 0.  ! ### Fix to avoid uninitialized value in loop on mw -- Correct?? ###
   

   hstarHLL = max((huL-huR+sE2*hR-sE1*hL)/(sE2-sE1),0.) ! middle state in an HLL solve
   
   !determine the middle entropy corrector wave------------------------
   rarecorrectortest=0.
   rarecorrector=0.
   if (rarecorrectortest .eq. 1.) then
      sdelta=lambda(3)-lambda(1)
      raremin = 0.5
      raremax = 0.9
      if ( (rare1 .eq. 1.) .and. (sE1*s1m .lt. 0.) ) then
         raremin=0.2
      end if
      if ( (rare2 .eq. 1.) .and. (sE2*s2m .lt. 0.) ) then
         raremin=0.2
      end if
      if ( (rare1 .eq. 1.) .or. (rare2 .eq. 1.) ) then
         !see which rarefaction is larger
         rare1st=3.*(sqrt(g*hL)-sqrt(g*hm))
         rare2st=3.*(sqrt(g*hR)-sqrt(g*hm))
         if ( (max(rare1st,rare2st) .gt. raremin*sdelta) .and. (max(rare1st,rare2st) .lt. raremax*sdelta) ) then
            rarecorrector=1.
            if (rare1st.gt.rare2st) then
               lambda(2)=s1m
            elseif (rare2st.gt.rare1st) then
               lambda(2)=s2m
            else
               lambda(2)=0.5*(s1m+s2m)
            end if
         end if
      end if
      if (hstarHLL .lt. (min(hL,hR)/5.) ) then
         rarecorrector=0.
      end if
   end if
   
   !     ## Is this correct 2-wave when rarecorrector == .true. ??
   do mw=1,mwaves
      r(1,mw)=1.
      r(2,mw)=lambda(mw)
      r(3,mw)=(lambda(mw))**2
   end do
   if (rarecorrector .eq. 0.) then
      lambda(2) = 0.5*(lambda(1)+lambda(3))
      r(1,2)=0.
      r(2,2)=0.
      r(3,2)=1.
   end if
   !     !---------------------------------------------------
   
   
   !determine the steady state wave -------------------
   !criticaltol = 1.d-6
   ! MODIFIED:
   criticaltol = max(drytol*g, 1e-6)
   criticaltol_2 = sqrt(criticaltol)
   deldelh = -delb
   deldelphi = -0.5 * (hR + hL) * (g * delb + delp / rho)
   
   !determine a few quanitites needed for steady state wave if iterated
   hLstar=hL
   hRstar=hR
   uLstar=uL
   uRstar=uR
   huLstar=uLstar*hLstar
   huRstar=uRstar*hRstar
   
   !iterate to better determine the steady state wave
   convergencetol=1e-6
   do iter=1,maxiter
      !determine steady state wave (this will be subtracted from the delta vectors)
      if ( (min(hLstar,hRstar) .lt. drytol) .and. (rarecorrector.eq.1.) ) then
         rarecorrector=0.
         hLstar=hL
         hRstar=hR
         uLstar=uL
         uRstar=uR
         huLstar=uLstar*hLstar
         huRstar=uRstar*hRstar
         lambda(2) = 0.5*(lambda(1)+lambda(3))
         r(1,2)=0.
         r(2,2)=0.
         r(3,2)=1.
      end if
      
      hbar =  max(0.5*(hLstar+hRstar),0.)
      s1s2bar = 0.25*(uLstar+uRstar)**2 - g*hbar
      s1s2tilde= max(0.,uLstar*uRstar) - g*hbar
      
      !find if sonic problem
      ! MODIFIED from 5.3.1 version
      sonic = 0.
      if (abs(s1s2bar) .le. criticaltol) then
         sonic = 1.
      elseif (s1s2bar*s1s2tilde .le. criticaltol**2) then
         sonic = 1.
      elseif (s1s2bar*sE1*sE2 .le. criticaltol**2) then
         sonic = 1.
      elseif (min(abs(sE1),abs(sE2)) .lt. criticaltol_2) then
         sonic = 1.
      elseif ( (sE1 .lt.  criticaltol_2) .and. (s1m .gt. -criticaltol_2) ) then
         sonic = 1.
      elseif ( (sE2 .gt. -criticaltol_2) .and. (s2m .lt. criticaltol_2) ) then
         sonic = 1.
      elseif ( (uL+sqrt(g*hL))*(uR+sqrt(g*hR)) .lt. 0.) then
         sonic = 1.
      elseif ((uL- sqrt(g*hL))*(uR- sqrt(g*hR)) .lt. 0.) then
         sonic = 1.
      end if

      !find jump in h, deldelh
      if (sonic.eq.1.) then
         deldelh =  -delb
      else
         deldelh = delb*g*hbar/s1s2bar
      end if
      !find bounds in case of critical state resonance, or negative states
      if ( (sE1.lt.-criticaltol) .and. (sE2.gt.criticaltol) ) then
         deldelh = min(deldelh,hstarHLL*(sE2-sE1)/sE2)
         deldelh = max(deldelh,hstarHLL*(sE2-sE1)/sE1)
      elseif (sE1.ge.criticaltol) then
         deldelh = min(deldelh,hstarHLL*(sE2-sE1)/sE1)
         deldelh = max(deldelh,-hL)
      elseif (sE2 .le. -criticaltol) then
         deldelh = min(deldelh,hR)
         deldelh = max(deldelh,hstarHLL*(sE2-sE1)/sE2)
      end if
      
      ! adjust deldelh for well-balancing of atmospheric pressure difference 
      if (g.ne.0.) then
         deldelh = deldelh - delp/(rho*g)
      end if
      !find jump in phi, deldelphi
      if (sonic.eq.1.) then
         deldelphi = -g*hbar*delb
      else
         deldelphi = -delb*g*hbar*s1s2tilde/s1s2bar
      end if
      !find bounds in case of critical state resonance, or negative states
      deldelphi=min(deldelphi,g*max(-hLstar*delb,-hRstar*delb))
      deldelphi=max(deldelphi,g*min(-hLstar*delb,-hRstar*delb))
      deldelphi = deldelphi - hbar * delp / rho
      
      del(1)=delh-deldelh
      del(2)=delhu
      del(3)=delphi-deldelphi

      !Determine determinant of eigenvector matrix========
      det1=r(1,1)*(r(2,2)*r(3,3)-r(2,3)*r(3,2))
      det2=r(1,2)*(r(2,1)*r(3,3)-r(2,3)*r(3,1))
      det3=r(1,3)*(r(2,1)*r(3,2)-r(2,2)*r(3,1))
      determinant=det1-det2+det3
      
      !solve for beta(k) using Cramers Rule=================
      do k=1,3
         do mw=1,3
            A(1,mw)=r(1,mw)
            A(2,mw)=r(2,mw)
            A(3,mw)=r(3,mw)
         end do
         A(1,k)=del(1)
         A(2,k)=del(2)
         A(3,k)=del(3)
         det1=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
         det2=A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
         det3=A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
         beta(k)=(det1-det2+det3)/determinant
      end do
      
      !exit if things aren't changing
      if (abs(del(1)**2+del(3)**2-delnorm).lt.convergencetol) exit
      
      delnorm = del(1)**2+del(3)**2
      !find new states qLstar and qRstar on either side of interface
      hLstar=hL
      hRstar=hR
      uLstar=uL
      uRstar=uR
      huLstar=uLstar*hLstar
      huRstar=uRstar*hRstar
      do mw=1,mwaves
         if (lambda(mw) .lt. 0.) then
            hLstar= hLstar + beta(mw)*r(1,mw)
            huLstar= huLstar + beta(mw)*r(2,mw)
         end if
      end do
      do mw=mwaves,1,-1
         if (lambda(mw).gt.0.) then
            hRstar= hRstar - beta(mw)*r(1,mw)
            huRstar= huRstar - beta(mw)*r(2,mw)
         end if
      end do
      if (hLstar.gt.drytol) then
         uLstar=huLstar/hLstar
      else
         hLstar=max(hLstar,0.)
         uLstar=0.
      end if
      if (hRstar.gt.drytol) then
         uRstar=huRstar/hRstar
      else
         hRstar=max(hRstar,0.)
         uRstar=0.
      end if
   end do ! end iteration on Riemann problem
   do mw=1,mwaves
      sw(mw)=lambda(mw)
      fw(1,mw)=beta(mw)*r(2,mw)
      fw(2,mw)=beta(mw)*r(3,mw)
      fw(3,mw)=beta(mw)*r(2,mw)
   end do
   !find transverse components (ie huv jumps).
   ! MODIFIED from 5.3.1 version
   fw(3,1)=fw(3,1)*vL
   fw(3,3)=fw(3,3)*vR
   fw(3,2)= 0.
   
   hustar_interface = huL + fw(1,1)   ! = huR - fw(1,3)
   if (hustar_interface .le. 0.) then
      fw(3,1) = fw(3,1) + (hR*uR*vR - hL*uL*vL - fw(3,1)- fw(3,3))
   else
      fw(3,3) = fw(3,3) + (hR*uR*vR - hL*uL*vL - fw(3,1)- fw(3,3))
   end if
    
 end subroutine riemann_aug_JCP
 



 
 subroutine fwave_limiter(i1,i2,wave,s,dtdx,fadd)

   integer, intent(in) :: i1,i2  
   real(amrex_real), intent(inout) :: wave(1:3,1:3,i1-3:i2+3)
   real(amrex_real), intent(in) :: s(1:3,i1-3:i2+3),dtdx
   real(amrex_real), intent(out) :: fadd(1:3,i1-3:i2+3)
   real(amrex_real) :: dotr(1:3),wnorm2,dotl,wlimiter,cqxx(1:3,i1-3:i2+3)
   integer :: m,mw,i 
   !%%%%%%%% limiter   
   cqxx = 0.
   dotr = 0.
   
   do i = i1-1,i2+1
      
      do mw=1,3
         
         ! Construct dot products
         wnorm2 = 0.
         dotl = dotr(mw)
         dotr(mw) = 0.
         do m=1,3
            wnorm2 = wnorm2 + wave(m,mw,i)**2
            dotr(mw) = dotr(mw) + wave(m,mw,i)*wave(m,mw,i+1)
         end do
         ! Skip this loop if it's on the boundary or the size of the wave is
         ! zero (but still want dot products to be initialized above)
         if (i .eq. i1-1) goto 110 
         if (i .eq. i2+1) goto 110
         if (wnorm2 .eq. 0.) goto 110
         ! Compute ratio of this wave's strength to upwind wave's strength
         if (s(mw,i) .gt. 0.) then
            !                wlimiter=max(0.,min((dotl/wnorm2),1.))  !minmod
            !                wlimiter = max(0.,min(1.,2.*(dotl/wnorm2)),min(2.,(dotl/wnorm2))) ! superbee
            !                wlimiter = ((dotl/wnorm2) + abs((dotl/wnorm2))) / (1. + abs((dotl/wnorm2))) !van Leer
                wlimiter=max(0.,min( 0.5*(1.+(dotl/wnorm2)),2.,2.*(dotl/wnorm2) )) ! monotinized centered
                !                wlimiter = (dotl/wnorm2)
             else
                !                wlimiter=max(0.,min((dotr(mw)/wnorm2),1.)) !minmod
                !                wlimiter = max(0.,min(1.,2.*(dotr(mw)/wnorm2)),min(2.,(dotr(mw)/wnorm2))) ! superbee
                !                wlimiter = ((dotr(mw)/wnorm2) + abs((dotr(mw)/wnorm2))) / (1. + abs((dotr(mw)/wnorm2))) !van Leer
                wlimiter=max(0.,min( 0.5*(1.+(dotr(mw)/wnorm2)),2.,2.*(dotr(mw)/wnorm2) )) ! monotinized centered     
                !                wlimiter = (dotl/wnorm2)
             end if
             
             ! Apply resulting limit
             do m=1,3
                wave(m,mw,i) = wlimiter * wave(m,mw,i)
             end do
110          continue
          end do
       end do
       !%%%%%%%%%%%%%%%%%   

   do i = i1-2,i2+3
      do m = 1,3
         cqxx(m,i) = 0.
      end do
      do mw=1,3
         do m=1,3
            ! second order corrections:
            cqxx(m,i) = cqxx(m,i) + sign(1.d0,s(mw,i))*(1.-abs(s(mw,i))*dtdx)*wave(m,mw,i)
         end do
      end do
      do m = 1,3
         fadd(m,i) = 0.5*cqxx(m,i)
      end do
   end do

  end subroutine fwave_limiter
   
 end module shallow_water_module
 
