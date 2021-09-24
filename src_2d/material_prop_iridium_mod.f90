module material_properties_iridium_module  

  use amrex_amr_module
  
  implicit none 

  ! -----------------------------------------------------------------
  ! This module is used to compute the material properties of
  ! iridium.
  ! -----------------------------------------------------------------
  
  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: get_Cp_iridium
  public :: get_ktherm_iridium
  public :: get_m_A_iridium
  public :: get_melting_point_iridium 
  public :: get_rho_iridium
  public :: get_electrical_resistivity_iridium
  public :: get_emissivity_iridium
  public :: get_enthalpy_of_vap_iridium
  public :: get_Richardson_iridium
  public :: get_surf_tension_iridium
  public :: get_thermelec_power_iridium
  public :: get_vapor_pressure_iridium
  public :: get_viscosity_iridium
  public :: get_work_function_iridium

contains 

  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the thermal conductivity
  ! -----------------------------------------------------------------
  subroutine get_ktherm_iridium(temp,ktherm)

    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: ktherm ! Thermal conductivity [W/mK] 
    
    if (temp.lt.300.0) then
       ktherm = 145.8297 
    else if(temp.lt.2719.0) then 
       ktherm = 169.1053 - 0.7527E6/temp**2 - 52.5153E-3*temp &
         + 9.3593E-6*temp**2
    else  
       ktherm = 18.0033 + 1.834E-2*temp
    end if
    
  end subroutine get_ktherm_iridium


  ! -----------------------------------------------------------------
  ! Subroutine used to computed the mass density
  ! -----------------------------------------------------------------
  subroutine get_rho_iridium(temp,rho)

    ! Input and output variables 
    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: rho    ! Mass density [kg/m3]

    ! Local variables
    real(amrex_real) :: temp_0 = 293.15
    
    if (temp.lt.300.0) then
       rho = 22.5592
    else if(temp.lt.2719.0) then 
       rho = 22.562 - 0.4068E-3*(temp-temp_0) - 9.4842E-8*(temp-temp_0)**2 - &
         7.4395E-12*(temp-temp_0)**3
    else
      rho = 19.5 - 0.85E-3*(temp-2719.0)
    end if
    ! Conversion from g/cm3 to kg/m3 
    rho = rho*1E3  
   
  end subroutine get_rho_iridium
  

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the heat capacity
  ! -----------------------------------------------------------------
  subroutine get_Cp_iridium(temp,Cp) 

    ! Input and output variables 
    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: Cp     ! Specific heat capacity [J/kgK]

    ! Local variables
    real(amrex_real) :: m_A ! Atomic mass [g/mol]
    
    if (temp.lt.300.0) then
       Cp = 25.041
    else if(temp.lt.2719.0) then 
       Cp = 25.1128 + 3.3293E-3*temp + 1.0441E-6*temp**2 - 0.1048E6/temp**2
    else 
       Cp = 44.210
    end if
    ! Conversion from J/(mol*K) to J/(kg*K) 
    call get_m_A_iridium(m_A)
    Cp = 1E3*(Cp/m_A)
    
  end subroutine get_Cp_iridium

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the electrical resistivity
  ! -----------------------------------------------------------------
  subroutine get_electrical_resistivity_iridium(temp,rho_e)

   ! Input and output variables 
   real(amrex_real), intent(in) :: temp      ! Temperature [K]
   real(amrex_real), intent(out) :: rho_e    ! Mass density [Ohm*m]

   ! Local variables
   real(amrex_real) :: temp_m = 2719.0
   
   if (temp.lt.300.0) then
      rho_e = 5.1995 
   else if(temp.lt.temp_m) then 
      rho_e = -0.290595 + 1.66688E-2*temp + 5.70656E-6*temp**2 - &
         8.95612E-10*temp**3
   else 
      rho_e = 60.5+1.316E-2*temp
   end if
   ! Conversion from uOhm*cm to Ohm*m 
   rho_e = rho_e*1E-8  
  
  end subroutine get_electrical_resistivity_iridium


  ! -----------------------------------------------------------------
  ! Subroutine used to computed the surface tension
  ! -----------------------------------------------------------------
  subroutine get_surf_tension_iridium(temp,sigma)

   real(amrex_real), intent(in) :: temp    ! Temperature [K]
   real(amrex_real), intent(out) :: sigma  ! Surface tension [N/m] 
   
   if (temp.lt.2719.0) then
      sigma = 2.262
   else 
      sigma = 2.262 - 0.28E-3*(temp-2719)
   end if
   
  end subroutine get_surf_tension_iridium


  ! -----------------------------------------------------------------
  ! Subroutine used to computed vapor pressure
  ! -----------------------------------------------------------------
  subroutine get_vapor_pressure_iridium(temp,pv)

   real(amrex_real), intent(in) :: temp    ! Temperature [K]
   real(amrex_real), intent(out) :: pv     ! Vapor pressure [Pa]
   
   pv = 10**(11.7592 - 30934.6/(temp-82.47))
   
  end subroutine get_vapor_pressure_iridium  


 ! -----------------------------------------------------------------
 ! Subroutine used to computed the viscosity
 ! -----------------------------------------------------------------
 subroutine get_viscosity_iridium(temp,mu)

   real(amrex_real), intent(in) :: temp   ! Temperature [K]
   real(amrex_real), intent(out) :: mu    ! Viscosity [N/m] 
   
   if (temp.lt.2719.0) then
      mu = 0.0085;
   else 
      mu = 0.59E-3*EXP(2.309*2719.0/temp)
   end if
   
 end subroutine get_viscosity_iridium  
  

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the atomic mass
  ! -----------------------------------------------------------------
  subroutine get_m_A_iridium(m_A)

    real(amrex_real), intent(out) :: m_A ! Atomic mass [g/mol]

    m_A = 192.217
    
  end subroutine get_m_A_iridium

  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the melting point properties
  ! -----------------------------------------------------------------
  subroutine get_melting_point_iridium(temp_melt, enth_fus, rho_melt)

    real(amrex_real), intent(out) :: temp_melt ! Temperature at melting [K] 
    real(amrex_real), intent(out) :: enth_fus ! Enthalpy of fusion [kJ/mol] 
    real(amrex_real), intent(out) :: rho_melt ! Density at melting [kg/m3] 

    temp_melt = 2719.0
    enth_fus = 39.37
    rho_melt = 19.5E3
    
  end subroutine get_melting_point_iridium

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the enthalpy of vaporization
  ! -----------------------------------------------------------------
  subroutine get_enthalpy_of_vap_iridium(enth_vap)

   real(amrex_real), intent(out) :: enth_vap ! Enthalpy of vaporization [kj/mol]

   enth_vap = 564.0
   
  end subroutine get_enthalpy_of_vap_iridium  


  ! -----------------------------------------------------------------
  ! Subroutine used to compute the work function
  ! -----------------------------------------------------------------
  subroutine get_work_function_iridium(Wf)

   real(amrex_real), intent(out) :: Wf ! Enthalpy of vaporization [J]

   Wf = 5.3
   ! Conversion from eV to J
   Wf = Wf*1.60218E-19
   
  end subroutine get_work_function_iridium 


  ! -----------------------------------------------------------------
  ! Subroutine used to compute Richardson constant
  ! -----------------------------------------------------------------
  subroutine get_Richardson_iridium(Aeff)

   real(amrex_real), intent(out) :: Aeff ! Richardson constant [A/m^2*K^2]

   Aeff = 100E4
   
  end subroutine get_Richardson_iridium


  ! -----------------------------------------------------------------
  ! Subroutine used to compute emissivity
  ! -----------------------------------------------------------------
  subroutine get_emissivity_iridium(temp, eps_t)

   real(amrex_real), intent(in) :: temp   ! Temperature [K]
   real(amrex_real), intent(out) :: eps_t ! Emissivity [dimensionless]

   if (temp.lt.300) then
      eps_t = 1.7802E-2 
   else if(temp.lt.2719) then
      eps_t = 1.7802E-2 + 1.2155E-4*(temp-300) - 1.5421E-8*(temp-300)**2
   else
      eps_t = 0.25
   endif
   
 end subroutine get_emissivity_iridium


  ! -----------------------------------------------------------------
  ! Subroutine used to compute the absolute thermoelectric power
  ! -----------------------------------------------------------------
 subroutine get_thermelec_power_iridium(temp, S)

   real(amrex_real), intent(in) :: temp   ! Temperature [K]
   real(amrex_real), intent(out) :: S     ! Thermoelectric power [V/K]

   if(temp.lt.300.0) then
      S = 0.86738
   elseif(temp.le.2719.0) then
      S = 0.86738 - 1.9990E-3*(temp-300) - 0.75620E-6*(temp-300)**2
   else
       S = 0 ! No available data, place-holder value
   endif

   ! Conversion from uV/K to V/K
   S = S*1E-6
   
 end subroutine get_thermelec_power_iridium

  
 
end module material_properties_iridium_module
