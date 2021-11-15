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
   public :: debug_cooling_fluxes  
   
   ! ------------------------------------------------------------------
   ! Public variables
   ! ------------------------------------------------------------------
   public :: input_time_mesh
   public :: input_surf_mesh
   public :: heatflux_table
   real(amrex_real), allocatable, dimension(:), save :: input_time_mesh
   real(amrex_real), allocatable, dimension(:), save :: input_surf_mesh
   real(amrex_real), allocatable, dimension(:,:), save :: heatflux_table
   
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
                                     temp, t_lo, t_hi, qb) 
 
      use read_input_module, only : plasma_flux_type, &
                                    cooling_thermionic, &
                                    cooling_vaporization, &
                                    cooling_radiation
      

      ! Input and output variables
      integer, intent(in) :: lo(2), hi(2)  
      integer, intent(in) :: id_lo(2), id_hi(2)
      integer, intent(in) :: t_lo(2), t_hi(2)
      real(amrex_real), intent(in) :: time
      real(amrex_real), intent(in) :: xlo(2)
      real(amrex_real), intent(in) :: dx(2)
      real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2))
      real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2))
      real(amrex_real), intent(out) :: qb(lo(1):hi(1),lo(2):hi(2))

      ! Local variables
      integer :: i, j, k, l
      real(amrex_real) :: q_plasma
      real(amrex_real) :: q_vap, q_rad, q_therm
      real(amrex_real) :: xpos, ypos

      qb = 0.0
      q_plasma = 0.0
      q_rad = 0.0
      q_vap = 0.0
      q_therm = 0.0

      do i = lo(1), hi(1)
         do j = lo(2), hi(2)

            ! Assign fluxes only on free surface
            if(nint(idom(i,j)).ne.0 .and. nint(idom(i,j+1)).eq.0) then

               ! Location of the free surface
               xpos = xlo(1) + (i-lo(1))*dx(1)

               ! Plasma flux
               if (plasma_flux_type.eq.'Gaussian') then
                  call gaussian_heat_flux(time, xpos, q_plasma)
               elseif (plasma_flux_type.eq.'Uniform') then
                  call uniform_heat_flux(time, xpos, q_plasma)
               elseif (plasma_flux_type.eq.'Input_file') then
                  call get_input_heat_flux (time, xpos, q_plasma)
               else
                  STOP "Unknown plasma heat flux type"
               end if

               ! Thermionic cooling flux
               if (cooling_thermionic) then
                  call thermionic_cooling(temp(i,j), q_plasma, q_therm)
               end if

               ! Vaporization cooling flux
               if (cooling_vaporization) then
                  call vaporization_cooling(temp(i ,j), q_vap)
               end if

               ! Radiative cooling flux
               if (cooling_radiation) then
                  call radiation_cooling(temp(i,j), q_rad)
               end if

               ! Sum all flux contributions
               qb(i,j) = q_plasma - q_rad - q_vap - q_therm
               
               ! Note: the term /dx(2) converts a surface heat flux [W/m^2]
               ! into a volumetric heat flux [W/m^3]
               qb(i,j) = qb(i,j)/dx(2)
               ! if(qb(i,j).le.0) then
               !    write(*,*) 'Negative surface heat flux at cell ' 
               !    write(*,*) i, j
               !    write(*,*) 'Where the surface temperature is '
               !    write(*,*) temp(i,j)      
               !    do k = lo(1), hi(1)
               !       do l = lo(2), hi(2)
               !          write(*,*) k,l
               !          write(*,*) temp(k,l)   
               !       end do 
               !    end do
               !    STOP 'Negative heat flux.'
               ! end if      
               
               if(qb(i,j).ne.qb(i,j)) then
                  ypos = xlo(2) + (j-lo(2))*dx(2)
                  write(*,*) 'Nan heat flux at cell ' 
                  write(*,*) i, j
                  write(*,*) xpos, ypos
                  write(*,*) 'Where the surface temperature is '
                  write(*,*) temp(i,j)
                  do k = lo(1), hi(1)
                     do l = lo(2), hi(2)
                        write(*,*) k,l
                        write(*,*) temp(k,l)   
                     end do 
                  end do
                  STOP 'Nan heat flux.'
               end if             
            end if
            
         end do
      end do
        
   end subroutine get_boundary_heat_flux
   
 
   ! -----------------------------------------------------------------
   ! Subroutine used to prescribe a gaussian heat flux active for
   ! time <= time_exposure
   ! -----------------------------------------------------------------   
   subroutine gaussian_heat_flux(time, xpos, qb) 
 
     use read_input_module, only : plasma_flux_params
     
     ! Input and output variables
     real(amrex_real), intent(in) :: time
     real(amrex_real), intent(in) :: xpos
     real(amrex_real), intent(out) :: qb
     
     
     qb = 0_amrex_real
     
     if (time.ge.plasma_flux_params(1) .and. time.le.plasma_flux_params(2)) then
        qb = plasma_flux_params(3) &
             *EXP(-((xpos-plasma_flux_params(4))**2)/(plasma_flux_params(5)**2))
     end if
     
   end subroutine gaussian_heat_flux
   
   ! -----------------------------------------------------------------
   ! Subroutine used to prescribe a uniform heat flux active for
   ! time <= time_exposure
   ! -----------------------------------------------------------------   
   subroutine uniform_heat_flux(time, xpos, qb) 
 
      use read_input_module, only : plasma_flux_params
      
      ! Input and output variables
      real(amrex_real), intent(in) :: time
      real(amrex_real), intent(in) :: xpos
      real(amrex_real), intent(out) :: qb
      
      
      qb = 0_amrex_real
      
      if (time.ge.plasma_flux_params(1) .and. time.le.plasma_flux_params(2)) then
         if(xpos.ge.plasma_flux_params(4) .and. xpos.le.plasma_flux_params(5)) then
            qb = plasma_flux_params(3) 
         end if
      end if
      
    end subroutine uniform_heat_flux  
    
   ! -----------------------------------------------------------------
   ! Subroutine used to prescribe a heat flux defined in the
   ! an input file.
   ! -----------------------------------------------------------------   
    subroutine get_input_heat_flux(time, xpos, qb) 
      
      ! Input and output variables
      real(amrex_real), intent(in) :: time
      real(amrex_real), intent(in) :: xpos
      real(amrex_real), intent(out) :: qb

      ! Local variables 
      integer :: i_x, i_t, n, m
      real(amrex_real) :: x(2), t(2), val(4)
      real(amrex_real) :: tx_query(2)
      
      
      qb = 0_amrex_real
      
      n = size(input_time_mesh,1)
      m = size(input_surf_mesh,1)
      call locate(input_time_mesh, n, time, i_t)
      call locate(input_surf_mesh, m, xpos, i_x)

      ! If query point is outside the 4 bounds
      ! then take the closest corner
      if(i_t.eq.0 .and. i_x.eq.0) then
         qb = heatflux_table(1,1)
      elseif (i_t.eq.n .and. i_x.eq.m) then
         qb = heatflux_table(n,m)
      elseif (i_t.eq.0 .and. i_x.eq.m) then
         qb = heatflux_table(1,m)
      elseif (i_t.eq.n .and. i_x.eq.0) then
         qb = heatflux_table(n,1)
      ! If query point is outside 2 bounds
      ! then linear interpolation
      elseif (i_t.eq.0 .or. i_t.eq.n) then
         x(1) = input_surf_mesh(i_x)
         x(2) = input_surf_mesh(i_x+1)
         if(i_t.eq.0) i_t = 1
         val(1) = heatflux_table(i_t, i_x)
         val(2) = heatflux_table(i_t, i_x+1)
         call lin_intrp(x, val(1:2), xpos, qb)
      elseif (i_x.eq.0 .or. i_x.eq.m) then
         t(1) = input_time_mesh(i_t)
         t(2) = input_time_mesh(i_t+1)
         if(i_x.eq.0) i_x = 1
         val(1) = heatflux_table(i_t, i_x)
         val(2) = heatflux_table(i_t+1, i_x)
         call lin_intrp(t, val(1:2), time, qb)
      else
         t(1) = input_time_mesh(i_t)
         t(2) = input_time_mesh(i_t+1)
         x(1) = input_surf_mesh(i_x)
         x(2) = input_surf_mesh(i_x+1)
         val(1) = heatflux_table(i_t,i_x)
         val(2) = heatflux_table(i_t+1,i_x)
         val(3) = heatflux_table(i_t+1,i_x+1)
         val(4) = heatflux_table(i_t,i_x+1)
         tx_query(1) = time
         tx_query(2) = xpos
         call bilin_intrp(t, x, val, tx_query, qb)
      endif
      
    end subroutine get_input_heat_flux
   
   ! -----------------------------------------------------------------
   ! Subroutine used to find the surface cooling flux due to 
   ! radiation given the surface temperature
   ! -----------------------------------------------------------------   
   subroutine radiation_cooling(Ts, q_rad)
 
    use material_properties_module, only : get_emissivity
 
    ! Input and output variables
    real(amrex_real), intent(in) :: Ts     ! Temperature at the center of cells adjacent to the free surface [K]
    real(amrex_real), intent(out) :: q_rad ! Radiated power [W/m^2]
 
    ! Local variables 
    real(amrex_real) :: eps_t ! Emissivity
    real(amrex_real) :: sigma = 5.670374419E-8 ! Stefan-Boltzmann constant [W/(m^-2*K^-4)]

    call get_emissivity(Ts, eps_t)
    q_rad = sigma*eps_t*Ts**4
 
   end subroutine radiation_cooling
   
 
   ! -----------------------------------------------------------------
   ! Subroutine used to find the surface cooling flux due to 
   ! thermionic emission given the surface temperature temperature
   ! see E. Thorén et al. Plasma Phys. Control. Fusion 63 035021 (2021)
   ! -----------------------------------------------------------------   
   subroutine thermionic_cooling(Ts, q_plasma, q_therm)
 
     use material_properties_module, only : get_work_function, &
                                            get_Richardson
     
     use read_input_module, only : thermionic_alpha
 
     ! Input and output variables                                       
     real(amrex_real), intent(in) :: Ts        ! Temperature at the center of cells adjacent to the free surface [K]
     real(amrex_real), intent(in) :: q_plasma  ! Plasma heat flux [K]
     real(amrex_real), intent(out) :: q_therm  ! Flux of energy due to thermionic emission [W/m^2]
     
     ! Local variables
     real(amrex_real) :: kb = 1.38064852E-23 ! Boltzmann constant [m^2*kg/(s^2*K)]
     real(amrex_real) :: Jth_nom
     real(amrex_real) :: Jth_lim
     real(amrex_real) :: Jth
     real(amrex_real) :: Aeff
     real(amrex_real) :: Wf
     real(amrex_real) :: e = 1.60217662E-19
     real(amrex_real) :: pi = 3.1415927
     
     call get_work_function(Wf)
     call get_Richardson(Aeff)

     ! Nominal thermionic current from the Richardson-Dushman formula
     Jth_nom = Aeff*EXP(-Wf/(kb*Ts))*Ts**2

     ! Space-charge limited current (semi-empirical expression)
     Jth_lim = 1.51e4 * q_plasma**(1.0/3.0) * (SIN(thermionic_alpha/180*pi))**2

     ! Minimum between nominal and space-charge limited
     Jth = MIN(Jth_lim, Jth_nom)

     ! Heat flux
     q_therm = Jth/e*(Wf+2*kb*Ts)
     
   end subroutine thermionic_cooling
 
 
   ! -----------------------------------------------------------------
   ! Subroutine used to find the surface cooling flux due to 
   ! vaporization given the surface temperature temperature
   ! see E. Thorén et al. Plasma Phys. Control. Fusion 63 035021 (2021)
   ! -----------------------------------------------------------------   
    subroutine vaporization_cooling(Ts, q_vap)
 
      use material_properties_module, only : get_vapor_pressure, &
                                             get_enthalpy_of_vaporization, &
                                             get_atomic_mass
 
      ! Input and output variables
      real(amrex_real), intent(in) :: Ts
      real(amrex_real), intent(out) :: q_vap
      ! Local variables
      real(amrex_real) :: gm
      real(amrex_real) :: kb = 1.38064852E-23 ! Boltzmann constant [m^2*kg/(s^2*K)]
      real(amrex_real) :: h_vap ! Enthalpy of vaporization [kJ/mol]
      real(amrex_real) :: pv ! Vapor pressure
      real(amrex_real) :: m_A ! Atomic mass [g/mol]
      real(amrex_real) :: pi = 3.1415927
      real(amrex_real) :: Na = 6.02214076E23 ! Avogadro's number [mol^-1]  
      
      call get_atomic_mass(m_A)
      call get_vapor_pressure(Ts, pv)
      call get_enthalpy_of_vaporization(Ts, h_vap)
      
      ! Conversion from g/mol to kg
      m_A = m_A*1E-3/Na
      
      ! Conversion from kJ/mol to J
      h_vap = h_vap*1E3/Na
      
      gm = pv*sqrt(m_A/(2*pi*kb*Ts))    
      q_vap = gm/m_A*(h_vap + 2*kb*Ts) 
      
    end subroutine vaporization_cooling

    ! -----------------------------------------------------------------
    ! Subroutine used to print the cooling fluxes as a function of
    ! the temperature to a file. Only useful for debugging purposes
    ! -----------------------------------------------------------------   
    subroutine debug_cooling_fluxes() 
 
      use read_input_module, only : cooling_debug, &
                                    thermionic_alpha

      integer :: i
      real(amrex_real) :: dT
      real(amrex_real) :: temp
      real(amrex_real) :: q_therm
      real(amrex_real) :: q_vap
      real(amrex_real) :: q_rad
      
      dT = (cooling_debug(3) - cooling_debug(2))/nint(cooling_debug(4))
      temp = cooling_debug(2)
      
      open (2, file = 'cooling_fluxes.dat', status = 'unknown')
      write(2, *) '# plasma flux and magnetic inclination: ', cooling_debug(5), thermionic_alpha
      write(2, *) '# Temperature[K], Thermionic flux [W/m^2], Radiative flux [W/m^2], Vaporization flux [W/m^2]'
      do i = 0,nint(cooling_debug(4)) 
         call thermionic_cooling(temp, cooling_debug(5), q_therm)
         call vaporization_cooling(temp, q_vap)
         call radiation_cooling(temp, q_rad)
         write(2,*) temp, q_therm, q_vap, q_rad
         temp = temp + dT
      end do
      close(2) 
      
    end subroutine debug_cooling_fluxes


   ! Given an array xx(1:n), and given a value x, returns a value j such that x is between
   ! xx(j) and xx(j+1). xx(1:n) must be monotonic, either increasing or decreasing. j=0
   ! or j=n is returned to indicate that x is out of range.
   subroutine locate(xx,n,x,j)

      ! Input and output variables

      integer, intent(in) :: n
      real(amrex_real), intent(in) :: xx(n)
      real(amrex_real), intent(in) :: x
      integer, intent(out) :: j

      ! Local variables

      integer jl,jm,ju

      jl=0 ! Initialize lower
      ju=n+1 ! upper limits.
      do while(ju-jl.gt.1)
         jm=(ju+jl)/2
         if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
            jl=jm
         else
            ju=jm
         endif
      end do
      if(x.eq.xx(1)) then
         j=1
      else if(x.eq.xx(n)) then
         j=n-1
      else
         j=jl
      endif
   end subroutine locate

   subroutine bilin_intrp(x, y, val, xy_query, val_query)
      ! Input and output variables
      real(amrex_real), intent(in) :: x(1:2)
      real(amrex_real), intent(in) :: y(1:2)
      real(amrex_real), intent(in) :: val(1:4)
      real(amrex_real), intent(in) :: xy_query(1:2)
      real(amrex_real), intent(out) :: val_query

      ! Local variables
      real(amrex_real) t, u

      t = (xy_query(1)-x(1))/(x(2)-x(1))
      u = (xy_query(2)-y(1))/(y(2)-y(1))

      val_query = (1-t)*(1-u)*val(1) + t*(1-u)*val(2) &
                  + t*u*val(3) + (1-t)*u*val(4)

   end subroutine bilin_intrp

   subroutine lin_intrp(x, val, x_query, val_query)
      ! Input and output variables
      real(amrex_real), intent(in) :: x(1:2)
      real(amrex_real), intent(in) :: val(1:2)
      real(amrex_real), intent(in) :: x_query
      real(amrex_real), intent(out) :: val_query

      ! Local variables
      real(amrex_real) t

      t = (x_query-x(1))/(x(2)-x(1))

      val_query = (1-t)*val(1) + t*val(2)

   end subroutine lin_intrp
    
 end module heat_flux_module
 
