module material_properties_niobium_module  

  use amrex_amr_module
  
  implicit none 

  ! -----------------------------------------------------------------
  ! This module is used to compute the material properties of
  ! niobium.
  ! -----------------------------------------------------------------
  
  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: get_Cp_niobium
  public :: get_ktherm_niobium
  public :: get_m_A_niobium
  public :: get_melting_point_niobium 
  public :: get_rho_niobium
  public :: get_electrical_resistivity_niobium
  public :: get_emissivity_niobium
  public :: get_Richardson_niobium
  public :: get_surf_tension_niobium
  public :: get_viscosity_niobium
  public :: get_work_function_niobium
  public :: get_enthalpy_of_vap_niobium
  public :: get_vapor_pressure_niobium
  public :: get_thermelec_power_niobium

contains 

  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the thermal conductivity
  ! -----------------------------------------------------------------
  subroutine get_ktherm_niobium(temp,ktherm)

    real(amrex_real), intent(in) :: temp    ! Temperature [K]
    real(amrex_real), intent(out) :: ktherm ! Thermal conductivity [W/mK] 
    
    if (temp.lt.300.0) then
       ktherm = 54.0812 
    else if(temp.lt.2745.0) then 
      ktherm = 45.27025 + 0.2305E6/temp**2 + 21.5609E-3*temp &
         - 2.4274E-6*temp**2
    else 
      ktherm = 11.3096 + 22.6768E-3*temp - 1.5865E-6*temp**2
    end if
    
  end subroutine get_ktherm_niobium


  ! -----------------------------------------------------------------
  ! Subroutine used to computed the mass density
  ! -----------------------------------------------------------------
  subroutine get_rho_niobium(temp,rho)

    ! Input and output variables 
    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: rho    ! Mass density [kg/m3]

    ! Local variables
    real(amrex_real) :: temp_0 = 293.15
    real(amrex_real) :: temp_m = 2745.0
    
    if (temp.lt.300.0) then
       rho = 8.5805
    else if(temp.lt.temp_m) then 
      rho = 8.582 - 0.1978E-3*(temp-temp_0) + 3.2727E-9*(temp-temp_0)**2 & 
         - 8.2514E-12*(temp-temp_0)**3
    else 
       rho = 7.73 - 0.39E-3*(temp-temp_m)
    end if
    ! Conversion from g/cm3 to kg/m3 
    rho = rho*1E3  
   
  end subroutine get_rho_niobium
  

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the heat capacity
  ! -----------------------------------------------------------------
  subroutine get_Cp_niobium(temp,Cp) 

    ! Input and output variables 
    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: Cp     ! Specific heat capacity [J/kgK]

    ! Local variables
    real(amrex_real) :: m_A ! Atomic mass [g/mol]
    
    if (temp.lt.300.0) then
      Cp = 24.4121
    elseif(temp.lt.2745.0) then
      Cp = 21.4900 + 11.6459E-3*temp - 6.9824E-6*temp**2 &
         + 2.1021E-9*temp**3
    else
      Cp = 43.294
    end if
    ! Conversion from J/(mol*K) to J/(kg*K) 
    call get_m_A_niobium(m_A)
    Cp = 1E3*(Cp/m_A)
    
  end subroutine get_Cp_niobium
  

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the electrical resistivity
  ! -----------------------------------------------------------------
  subroutine get_electrical_resistivity_niobium(temp,rho_e)

    ! Input and output variables 
    real(amrex_real), intent(in) :: temp      ! Temperature [K]
    real(amrex_real), intent(out) :: rho_e    ! Mass density [Ohm*m]

    ! Local variables
    real(amrex_real) :: temp_m = 2745.0
    
    if (temp.lt.300.0) then
       rho_e = 15.639
    else if(temp.lt.temp_m) then 
       rho_e = 1.8095 + 4.887E-2*temp - 9.7589E-6*temp**2 &
         + 1.7339E-9*temp**3
    else 
      rho_e = 68.33 + 1.477E-2*temp
    end if
    ! Conversion from uOhm*cm to Ohm*m 
    rho_e = rho_e*1E-8  
   
  end subroutine get_electrical_resistivity_niobium


  ! -----------------------------------------------------------------
  ! Subroutine used to computed the surface tension
  ! -----------------------------------------------------------------
  subroutine get_surf_tension_niobium(temp,sigma)

    real(amrex_real), intent(in) :: temp    ! Temperature [K]
    real(amrex_real), intent(out) :: sigma  ! Surface tension [N/m] 
    
    if (temp.lt.2745.0) then
       sigma = 1.967
    else 
       sigma = 1.967 - 0.17E-3*(temp-2745.0)
    end if
    
  end subroutine get_surf_tension_niobium


  ! -----------------------------------------------------------------
  ! Subroutine used to computed vapor pressure
  ! -----------------------------------------------------------------
  subroutine get_vapor_pressure_niobium(temp,pv)

    real(amrex_real), intent(in) :: temp    ! Temperature [K]
    real(amrex_real), intent(out) :: pv     ! Vapor pressure [Pa]
    
    pv = 10**(12.0502-35189.3/(temp-21.93))
    
  end subroutine get_vapor_pressure_niobium  
  

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the viscosity
  ! -----------------------------------------------------------------
  subroutine get_viscosity_niobium(temp,mu)

    real(amrex_real), intent(in) :: temp    ! Temperature [K]
    real(amrex_real), intent(out) :: mu     ! Viscosity [N/m] 
    
    if (temp.lt.2745.0) then
       mu = 0.0050;
    else 
      mu = 1.05E-3*EXP(1.556*2745/temp)
    end if
    
  end subroutine get_viscosity_niobium  
  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the atomic mass
  ! -----------------------------------------------------------------
  subroutine get_m_A_niobium(m_A)

    real(amrex_real), intent(out) :: m_A ! Atomic mass [g/mol]

    m_A = 92.90638
    
  end subroutine get_m_A_niobium

  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the melting point properties
  ! -----------------------------------------------------------------
  subroutine get_melting_point_niobium(temp_melt, enth_fus, rho_melt)

    real(amrex_real), intent(out) :: temp_melt ! Temperature at melting [K] 
    real(amrex_real), intent(out) :: enth_fus   ! Enthalpy of fusion [kJ/mol]
    real(amrex_real), intent(out) :: rho_melt   ! Density at metling [Kg/m^3]

    temp_melt = 2745.0
    enth_fus = 32.98
    rho_melt = 7730
    
  end subroutine get_melting_point_niobium


  ! -----------------------------------------------------------------
  ! Subroutine used to computed the enthalpy of vaporization
  ! -----------------------------------------------------------------
  subroutine get_enthalpy_of_vap_niobium(enth_vap)

    real(amrex_real), intent(out) :: enth_vap ! Enthalpy of vaporization [kj/mol]

    enth_vap = 694.0
    
  end subroutine get_enthalpy_of_vap_niobium  


  ! -----------------------------------------------------------------
  ! Subroutine used to compute the work function
  ! -----------------------------------------------------------------
  subroutine get_work_function_niobium(Wf)

    real(amrex_real), intent(out) :: Wf ! Enthalpy of vaporization [J]

    Wf = 4.2
    
    ! Conversion from eV to J
    Wf = Wf*1.60218E-19
    
  end subroutine get_work_function_niobium    


  ! -----------------------------------------------------------------
  ! Subroutine used to compute Richardson constant
  ! -----------------------------------------------------------------
  subroutine get_Richardson_niobium(Aeff)

    real(amrex_real), intent(out) :: Aeff ! Richardson constant [A/m^2*K^2]

    Aeff = 50E4
    
  end subroutine get_Richardson_niobium


  ! -----------------------------------------------------------------
  ! Subroutine used to compute emissivity
  ! -----------------------------------------------------------------
  subroutine get_emissivity_niobium(temp, eps_t)

    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: eps_t ! Emissivity [dimensionless]

    if (temp.lt.300.0) then
      eps_t = 2.68362E-2
    else if(temp.lt.2745.0) then
      eps_t = 2.68362E-2 + 1.54653E-4*(temp-300) - 2.23163E-8*(temp-300)**2
    else
      eps_t = 0.29
    endif
    
  end subroutine get_emissivity_niobium      

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the absolute thermoelectric power
  ! -----------------------------------------------------------------
  subroutine get_thermelec_power_niobium(temp, S)

    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: S     ! Thermoelectric power [V/K]

    if(temp.lt.300.0) then
      S = -1.2462
    elseif(temp.le.1700.0) then
      S = 88.280 + 2.69689E6/temp**2 - 24.97805E3/temp - 155.698E-3*temp + &
         126.954E-6*temp**2 - 35.10907E-9*temp**3
    elseif(temp.le.2745.0) then
      S = 65.199852 - 50.81273E-3*temp + 8.796398E-6*temp**2
    else
      S = 0 ! No available data, place-holder value
    endif

    ! Convversion from uV/K to V/K
    S = S*1E-6
    
  end subroutine get_thermelec_power_niobium

  
end module material_properties_niobium_module
