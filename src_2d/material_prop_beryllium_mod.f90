module material_properties_beryllium_module  

  use amrex_amr_module
  
  implicit none 

  ! -----------------------------------------------------------------
  ! This module is used to compute the material properties of
  ! beryllium as described in P. Tolias, Nucl. Mater. Energy 13, 42
  ! (2017)
  ! -----------------------------------------------------------------
  
  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: get_Cp_beryllium
  public :: get_ktherm_beryllium
  public :: get_m_A_beryllium
  public :: get_melting_point_beryllium 
  public :: get_rho_beryllium
  public :: get_electrical_resistivity_beryllium
  public :: get_emissivity_beryllium
  public :: get_Richardson_beryllium
  public :: get_surf_tension_beryllium
  public :: get_viscosity_beryllium
  public :: get_work_function_beryllium
  public :: get_enthalpy_of_vap_beryllium
  public :: get_hcp_to_bcc_point_beryllium
  public :: get_vapor_pressure_beryllium
  public :: get_thermelec_power_beryllium

contains 

  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the thermal conductivity
  ! -----------------------------------------------------------------
  subroutine get_ktherm_beryllium(temp,ktherm)

    real(amrex_real), intent(in) :: temp    ! Temperature [K]
    real(amrex_real), intent(out) :: ktherm ! Thermal conductivity [W/mK] 
    
    if (temp.lt.300.0) then
       ktherm = 199.7338  
    else if(temp.lt.1560.0) then 
       ktherm = 148.8912 + 6.5407E6/temp**2 - 76.3780E-3*temp+ &
            12.0174E-6*temp**2
    else 
       ktherm = 76.44+35E-3*(temp-1560.0)
    end if
    
  end subroutine get_ktherm_beryllium


  ! -----------------------------------------------------------------
  ! Subroutine used to computed the mass density
  ! -----------------------------------------------------------------
  subroutine get_rho_beryllium(temp,rho)

    ! Input and output variables 
    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: rho    ! Mass density [kg/m3]

    ! Local variables
    real(amrex_real) :: temp_0 = 300.0
    real(amrex_real) :: temp_m = 1560.0
    
    if (temp.lt.300.0) then
       rho = 1.85 
    else if(temp.lt.temp_m) then 
      rho = 1.85 - 6.86479E-5*(temp-temp_0) - 4.1660E-8*(temp-temp_0)**2 &
           + 1.1354E-11*(temp-temp_0)**3
    else 
       rho = 1.690-0.116E-3*(temp-temp_m)
    end if
    ! Conversion from g/cm3 to kg/m3 
    rho = rho*1E3  
   
  end subroutine get_rho_beryllium
  

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the heat capacity
  ! -----------------------------------------------------------------
  subroutine get_Cp_beryllium(temp,Cp) 

    ! Input and output variables 
    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: Cp     ! Specific heat capacity [J/kgK]

    ! Local variables
    real(amrex_real) :: m_A ! Atomic mass [g/mol]
    
    if (temp.lt.300.0) then
       Cp = 16.5438629378 
    else if(temp.lt.1543.0) then 
       Cp = 21.5390 + 4.94572E-3*temp + 1.356324E-6*temp**2 - 0.594083E6/temp**2
    else if(temp.lt.1560) then
       Cp = 30.0 
    else 
       Cp = 25.4345+2.15008E-3*temp 
    end if
    ! Conversion from J/(mol*K) to J/(kg*K) 
    call get_m_A_beryllium(m_A)
    Cp = 1E3*(Cp/m_A)
    
  end subroutine get_Cp_beryllium
  

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the electrical resistivity
  ! -----------------------------------------------------------------
  subroutine get_electrical_resistivity_beryllium(temp,rho_e)

    ! Input and output variables 
    real(amrex_real), intent(in) :: temp      ! Temperature [K]
    real(amrex_real), intent(out) :: rho_e    ! Mass density [Ohm*m]

    ! Local variables
    real(amrex_real) :: temp_0 = 300.0
    real(amrex_real) :: temp_m = 1560.0
    
    if (temp.lt.300.0) then
       rho_e = 3.71 
    else if(temp.lt.temp_m) then 
       rho_e = 3.71 + 30.412E-3*(temp-temp_0) + 2.785E-6*(temp-temp_0)**2 &
            + 3.252E-9*(temp-temp_0)**3
    else 
       rho_e = 45
    end if
    ! Conversion from uOhm*cm to Ohm*m 
    rho_e = rho_e*1E-8  
   
  end subroutine get_electrical_resistivity_beryllium


  ! -----------------------------------------------------------------
  ! Subroutine used to computed the surface tension
  ! -----------------------------------------------------------------
  subroutine get_surf_tension_beryllium(temp,sigma)

    real(amrex_real), intent(in) :: temp    ! Temperature [K]
    real(amrex_real), intent(out) :: sigma  ! Surface tension [N/m] 
    
    if (temp.lt.1560.0) then
       sigma = 1.143
    else 
       sigma = 1.143 - 0.20E-3*(temp-1560.0)
    end if
    
  end subroutine get_surf_tension_beryllium


  ! -----------------------------------------------------------------
  ! Subroutine used to computed vapor pressure
  ! -----------------------------------------------------------------
  subroutine get_vapor_pressure_beryllium(temp,pv)

    real(amrex_real), intent(in) :: temp    ! Temperature [K]
    real(amrex_real), intent(out) :: pv     ! Vapor pressure [Pa]
    
    pv = 10**(10.2089-13696.6/(temp-124.63))
    
  end subroutine get_vapor_pressure_beryllium  
  

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the viscosity
  ! -----------------------------------------------------------------
  subroutine get_viscosity_beryllium(temp,mu)

    real(amrex_real), intent(in) :: temp ! Temperature [K]
    real(amrex_real), intent(out) :: mu  ! Viscosity [N/m] 
    
    if (temp.lt.1560.0) then
       mu = 0.005090697766920;
    else 
       mu = 0.1E-3*EXP(3.93*1560.0/temp)
    end if
    
  end subroutine get_viscosity_beryllium  
  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the atomic mass
  ! -----------------------------------------------------------------
  subroutine get_m_A_beryllium(m_A)

    real(amrex_real), intent(out) :: m_A ! Atomic mass [g/mol]

    m_A = 9.0121831
    
  end subroutine get_m_A_beryllium

  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the melting point properties
  ! -----------------------------------------------------------------
  subroutine get_melting_point_beryllium(temp_melt, enth_fus, rho_melt)

    real(amrex_real), intent(out) :: temp_melt ! Temperature at melting [K] 
    real(amrex_real), intent(out) :: enth_fus   ! Enthalpy of fusion [kJ/mol]
    real(amrex_real), intent(out) :: rho_melt   ! Density at metling [Kg/m^3]

    temp_melt = 1560.0
    enth_fus = 7.959
    rho_melt = 1690.0
    
  end subroutine get_melting_point_beryllium


  ! -----------------------------------------------------------------
  ! Subroutine used to computed the enthalpy of vaporization
  ! -----------------------------------------------------------------
  subroutine get_enthalpy_of_vap_beryllium(enth_vap)

    real(amrex_real), intent(out) :: enth_vap ! Enthalpy of vaporization [kj/mol]

    enth_vap = 324.0
    
  end subroutine get_enthalpy_of_vap_beryllium  


  ! -----------------------------------------------------------------
  ! Subroutine used to computed the hcp to bcc point properties 
  ! -----------------------------------------------------------------
  subroutine get_hcp_to_bcc_point_beryllium(enth_hcp2bcc, hcp2bcc_point, rho_hcp2bcc)

    real(amrex_real), intent(out) :: enth_hcp2bcc     ! Enthalpy of vaporization [kj/mol]
    real(amrex_real), intent(out) :: hcp2bcc_point    ! Temperature of phase transition [K]
    real(amrex_real), intent(out) :: rho_hcp2bcc      ! Density at solid phase transition point [kg/m^3]

    enth_hcp2bcc = 6.855
    hcp2bcc_point = 1543.0
    rho_hcp2bcc = 1722.1
    
  end subroutine get_hcp_to_bcc_point_beryllium  



  ! -----------------------------------------------------------------
  ! Subroutine used to compute the work function
  ! -----------------------------------------------------------------
  subroutine get_work_function_beryllium(Wf)

    real(amrex_real), intent(out) :: Wf ! Enthalpy of vaporization [J]

    Wf = 4.98
    ! Conversion from eV to Joule
    Wf = Wf*1.60218E-19
    
  end subroutine get_work_function_beryllium    


  ! -----------------------------------------------------------------
  ! Subroutine used to compute Richardson constant
  ! -----------------------------------------------------------------
  subroutine get_Richardson_beryllium(Aeff)

    real(amrex_real), intent(out) :: Aeff ! Richardson constant [A/m^2*K^2]

    Aeff = 120E4
    
  end subroutine get_Richardson_beryllium


  ! -----------------------------------------------------------------
  ! Subroutine used to compute emissivity
  ! -----------------------------------------------------------------
  subroutine get_emissivity_beryllium(temp, eps_t)

    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: eps_t ! Emissivity [dimensionless]

    if (temp.lt.300.0) then
       eps_t = 4.3865E-2
    else if(temp.lt.1560.0) then
       eps_t = 4.3865E-2 + 0.5728E-4*(temp-300) - &
            0.2184E-6*(temp-300)**2 + 0.52076E-9*(temp-300)**3
    else
       eps_t = 0.811
    endif
    
  end subroutine get_emissivity_beryllium      

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the absolute thermoelectric power
  ! -----------------------------------------------------------------
  subroutine get_thermelec_power_beryllium(temp, S)

    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: S     ! Thermoelectric power [V/K]

    if(temp.lt.300.0) then
       S = 6.75364
    elseif(temp.le.1560.0) then
       S = 6.75364 + 43.09365E-3*(temp-300) - 8.21233E-6*(temp-300)**2 &
            + 3.17939E-9*(temp-300)**3
    else
       S = 0 ! No available data, place-holder value
    endif

    ! Conversion from uV/K to V/K
    S = S*1E-6
    
  end subroutine get_thermelec_power_beryllium

  
end module material_properties_beryllium_module
