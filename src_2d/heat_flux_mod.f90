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
                                     temp, t_lo, t_hi, qb) 
 
      use read_input_module, only : flux_type, &
                                    cool_flux

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
      integer :: i, j
      real(amrex_real) :: q_plasma
      real(amrex_real) :: q_vap, q_rad, q_therm
      real(amrex_real) :: xpos

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
               if (flux_type.eq.'Gaussian') then
                  call gaussian_heat_flux(time, xpos, q_plasma)
               end if

               ! Radiative cooling flux
               if (cool_flux.eq.'all' .or. cool_flux.eq.'radiation') then
                  call radiation_cooling(temp(i, j), q_rad)
               end if

               ! Thermionic cooling flux
               if (cool_flux.eq.'all' .or. cool_flux.eq.'thermionic') then
                  call thermionic_cooling(temp(i, j), q_therm)
               end if

               ! Vapor cooling flux
               if (cool_flux.eq.'all' .or. cool_flux.eq.'vaporization') then
                  call vaporization_cooling(temp(i ,j), q_vap)
               end if

               ! Sum all flux contributions
               qb(i,j) = q_plasma - q_rad - q_vap - q_therm
               
               ! Note: the term /dx(2) converts a surface heat flux [W/m^2]
               ! into a volumetric heat flux [W/m^3]
               qb(i,j) = qb(i,j)/dx(2)
               
            end if
            
         end do
      end do
        
   end subroutine get_boundary_heat_flux
   
 
   ! -----------------------------------------------------------------
   ! Subroutine used to prescribe a gaussian heat flux active for
   ! time <= time_exposure
   ! -----------------------------------------------------------------   
   subroutine gaussian_heat_flux(time, xpos, qb) 
 
     use read_input_module, only : flux_params
     
     ! Input and output variables
     real(amrex_real), intent(in) :: time
     real(amrex_real), intent(in) :: xpos
     real(amrex_real), intent(out) :: qb
     
     
     qb = 0_amrex_real
     
     if (time.ge.flux_params(1) .and. time.lt.flux_params(2)) then
        qb = flux_params(3) &
             *EXP(-((xpos-flux_params(4))**2)/(flux_params(5)**2))
     end if
     
   end subroutine gaussian_heat_flux
   

   ! -----------------------------------------------------------------
   ! Subroutine used to find the surface cooling flux due to 
   ! radiation given the surface temperature
   ! -----------------------------------------------------------------   
   subroutine radiation_cooling (Ts, q_rad)
 
    use material_properties_module, only : get_emissivity
 
    ! Input and output variables
    real(amrex_real), intent(in) :: Ts     ! Temperature at the center of cells adjacent to the free surface [K]
    real(amrex_real), intent(out) :: q_rad ! Radiated power [W/m^2]
 
    ! Local variables 
    real(amrex_real) :: eps_t ! Emissivity
    real(amrex_real) :: sigma = 5.670374419E-8 ! Stefan-Bolzmann constant [W/(m^-2*K^-4)]

    call get_emissivity(Ts, eps_t)
    q_rad = sigma*eps_t*Ts**4
 
   end subroutine radiation_cooling
   
 
   ! -----------------------------------------------------------------
   ! Subroutine used to find the surface cooling flux due to 
   ! thermionic emission given the surface temperature temperature
   ! see E. Thorén et al. Plasma Phys. Control. Fusion 63 035021 (2021)
   ! -----------------------------------------------------------------   
   subroutine thermionic_cooling (Ts, q_therm)
 
     use material_properties_module, only : get_work_function, &
                                            get_Richardson
     
     use read_input_module, only : u_the, &
                                   n_e, &
                                   alpha_B
 
     ! Input and output variables                                       
     real(amrex_real), intent(in) :: Ts        ! Temperature at the center of cells adjacent to the free surface [K]
     real(amrex_real), intent(out) :: q_therm  ! Flux of energy due to thermionic emission [W/m^2]
     
     ! Local variables
     real(amrex_real) :: kb = 1.38064852E-23 ! Bolzmann constant [m^2*kg/(s^2*K)]
     real(amrex_real) :: Jth_nom
     real(amrex_real) :: Jth_lim
     real(amrex_real) :: Jth
     real(amrex_real) :: Aeff
     real(amrex_real) :: Wf
     real(amrex_real) :: e = 1.60217662E-19
     real(amrex_real) :: pi = 3.1415927
     
     call get_work_function(Wf)
     call get_Richardson(Aeff)
     
     Jth_lim = 0.43*e*n_e*u_the*(SIN(alpha_B/180*pi))**2  
     Jth_nom = Aeff*EXP(-Wf/(kb*Ts))*Ts**2                
     Jth = MIN(Jth_lim, Jth_nom)
     
     q_therm = Jth/e*(Wf+2*kb*Ts)
     
   end subroutine thermionic_cooling
 
 
   ! -----------------------------------------------------------------
   ! Subroutine used to find the surface cooling flux due to 
   ! vaporization given the surface temperature temperature
   ! see E. Thorén et al. Plasma Phys. Control. Fusion 63 035021 (2021)
   ! -----------------------------------------------------------------   
    subroutine vaporization_cooling(Ts, q_vap)
 
      use material_properties_module, only : get_vapor_pressure, &
                                             get_enthalpy_of_vap, &
                                             get_m_A
 
      ! Input and output variables
      real(amrex_real), intent(in) :: Ts
      real(amrex_real), intent(out) :: q_vap
      ! Local variables
      real(amrex_real) :: gm
      real(amrex_real) :: kb = 1.38064852E-23 ! Bolzmann constant [m^2*kg/(s^2*K)]
      real(amrex_real) :: h_vap ! Enthalpy of vaporization [kJ/mol]
      real(amrex_real) :: pv ! Vapor pressure
      real(amrex_real) :: m_A ! Atomic mass [g/mol]
      real(amrex_real) :: pi = 3.1415927
      real(amrex_real) :: Na = 6.02214076E23 ! Avogadro's number [mol^-1]  
      
      call get_m_A(m_A)
      call get_vapor_pressure(Ts, pv)
      call get_enthalpy_of_vap(Ts, h_vap)
      
      ! Conversion from g/mol to kg
      m_A = m_A*1E-3/Na
      
      ! Conversion from kJ/mol to J
      h_vap = h_vap*1E3/Na
      
      gm = pv*sqrt(m_A/(2*pi*kb*Ts))    
      q_vap = gm/m_A*(h_vap + 2*kb*Ts) 
      
    end subroutine vaporization_cooling

 
 end module heat_flux_module
 
