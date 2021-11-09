module material_properties_test_evaporation_module  

  use amrex_amr_module
  
  implicit none 

  ! -----------------------------------------------------------------
  ! This module is used to compute the material properties of
  ! the material used to test the evaporation routines.
  ! -----------------------------------------------------------------
  
  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: get_heat_capacity_test_evaporation
  public :: get_conductivity_test_evaporation
  public :: get_atomic_mass_test_evaporation
  public :: get_melting_point_test_evaporation 
  public :: get_mass_density_test_evaporation
  public :: get_electrical_resistivity_test_evaporation
  public :: get_emissivity_test_evaporation
  public :: get_Richardson_test_evaporation
  public :: get_surface_tension_test_evaporation
  public :: get_viscosity_test_evaporation
  public :: get_work_function_test_evaporation
  public :: get_enthalpy_of_vaporization_test_evaporation
  public :: get_vapor_pressure_test_evaporation
  public :: get_thermelectric_power_test_evaporation

contains 

  
  ! -----------------------------------------------------------------
  ! Subroutine used to compute the thermal conductivity
  ! -----------------------------------------------------------------
  subroutine get_conductivity_test_evaporation(temp,ktherm)

    real(amrex_real), intent(in) :: temp    ! Temperature [K]
    real(amrex_real), intent(out) :: ktherm ! Thermal conductivity [W/mK] 
    
    if (temp.lt.1560.0) then
       ktherm = 60 
    else 
       ktherm = 75
    end if
    
  end subroutine get_conductivity_test_evaporation


  ! -----------------------------------------------------------------
  ! Subroutine used to compute the mass density
  ! -----------------------------------------------------------------
  subroutine get_mass_density_test_evaporation(rho)

    ! Input and output variables 
    real(amrex_real), intent(out) :: rho    ! Mass density [kg/m3]
    
    rho = 1705.0
   
  end subroutine get_mass_density_test_evaporation
  

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the heat capacity
  ! -----------------------------------------------------------------
  subroutine get_heat_capacity_test_evaporation(temp,Cp) 

    ! Input and output variables 
    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: Cp     ! Specific heat capacity [J/kgK]

    ! Local variables
    real(amrex_real) :: m_A ! Atomic mass [g/mol]
    
    if (temp.lt.1560.0) then
       Cp = 30.0
    else 
       Cp = 30.0
    end if
    ! Conversion from J/(mol*K) to J/(kg*K) 
    call get_atomic_mass_test_evaporation(m_A)
    Cp = 1E3*(Cp/m_A)
    
  end subroutine get_heat_capacity_test_evaporation
  

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the electrical resistivity
  ! -----------------------------------------------------------------
  subroutine get_electrical_resistivity_test_evaporation(temp,rho_e)

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
   
  end subroutine get_electrical_resistivity_test_evaporation


  ! -----------------------------------------------------------------
  ! Subroutine used to compute the surface tension
  ! -----------------------------------------------------------------
  subroutine get_surface_tension_test_evaporation(temp,sigma)

    real(amrex_real), intent(in) :: temp    ! Temperature [K]
    real(amrex_real), intent(out) :: sigma  ! Surface tension [N/m] 
    
    if (temp.lt.1560.0) then
       sigma = 1.143
    else 
       sigma = 1.143 - 0.20E-3*(temp-1560.0)
    end if
    
  end subroutine get_surface_tension_test_evaporation


  ! -----------------------------------------------------------------
  ! Subroutine used to compute vapor pressure
  ! -----------------------------------------------------------------
  subroutine get_vapor_pressure_test_evaporation(temp,pv)

    real(amrex_real), intent(in) :: temp    ! Temperature [K]
    real(amrex_real), intent(out) :: pv     ! Vapor pressure [Pa]

    if(temp.lt.300.0) then
      pv = 0.0
    else
      pv = 10**(10.2089-13696.6/(temp-124.63))
    end if
    
  end subroutine get_vapor_pressure_test_evaporation  
  

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the viscosity
  ! -----------------------------------------------------------------
  subroutine get_viscosity_test_evaporation(temp,mu)

    real(amrex_real), intent(in) :: temp ! Temperature [K]
    real(amrex_real), intent(out) :: mu  ! Viscosity [Pa*s] 
    
    if (temp.lt.1560.0) then
       mu = 0.005090697766920
    else 
       mu = 0.1E-3*EXP(3.93*1560.0/temp)
    end if
    
  end subroutine get_viscosity_test_evaporation  
  
  ! -----------------------------------------------------------------
  ! Subroutine used to compute the atomic mass
  ! -----------------------------------------------------------------
  subroutine get_atomic_mass_test_evaporation(m_A)

    real(amrex_real), intent(out) :: m_A ! Atomic mass [g/mol]

    m_A = 9.0121831
    
  end subroutine get_atomic_mass_test_evaporation

  
  ! -----------------------------------------------------------------
  ! Subroutine used to compute the melting point properties
  ! -----------------------------------------------------------------
  subroutine get_melting_point_test_evaporation(temp_melt, enth_fus, rho_melt)

    real(amrex_real), intent(out) :: temp_melt ! Temperature at melting [K] 
    real(amrex_real), intent(out) :: enth_fus   ! Enthalpy of fusion [kJ/mol]
    real(amrex_real), intent(out) :: rho_melt   ! Density at metling [Kg/m^3]

    temp_melt = 1560.0
    enth_fus = 7.959 + 6.855 ! The latent heat of the hcp-to-bcc transition is included in the enthalpy of fusion
    rho_melt = 1705.0
    
  end subroutine get_melting_point_test_evaporation


  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy of vaporization
  ! -----------------------------------------------------------------
  subroutine get_enthalpy_of_vaporization_test_evaporation(enth_vap)

    real(amrex_real), intent(out) :: enth_vap ! Enthalpy of vaporization [kj/mol]

    enth_vap = 324.0
    
  end subroutine get_enthalpy_of_vaporization_test_evaporation  


  ! -----------------------------------------------------------------
  ! Subroutine used to compute the work function
  ! -----------------------------------------------------------------
  subroutine get_work_function_test_evaporation(Wf)

    real(amrex_real), intent(out) :: Wf ! Enthalpy of vaporization [J]

    Wf = 4.98
    ! Conversion from eV to Joule
    Wf = Wf*1.60218E-19
    
  end subroutine get_work_function_test_evaporation    


  ! -----------------------------------------------------------------
  ! Subroutine used to compute Richardson constant
  ! -----------------------------------------------------------------
  subroutine get_Richardson_test_evaporation(Aeff)
    
    real(amrex_real), intent(out) :: Aeff ! Richardson constant [A/m^2*K^2]

    Aeff = 120E4
    
  end subroutine get_Richardson_test_evaporation


  ! -----------------------------------------------------------------
  ! Subroutine used to compute emissivity
  ! -----------------------------------------------------------------
  subroutine get_emissivity_test_evaporation(eps_t)

    real(amrex_real), intent(out) :: eps_t ! Emissivity [dimensionless]

    eps_t = 0.0
    
  end subroutine get_emissivity_test_evaporation      

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the absolute thermoelectric power
  ! -----------------------------------------------------------------
  subroutine get_thermelectric_power_test_evaporation(temp, S)

    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: S     ! Thermoelectric power [V/K]
    
    if(temp.lt.300.0) then
       S = 6.75364
    elseif(temp.lt.1560.0) then
       S = 6.75364 + 43.09365E-3*(temp-300) - 8.21233E-6*(temp-300)**2 &
            + 3.17939E-9*(temp-300)**3
    else
       S = 0 ! No available data, place-holder value
    endif

    ! Conversion from uV/K to V/K
    S = S*1E-6
    
  end subroutine get_thermelectric_power_test_evaporation

  
end module material_properties_test_evaporation_module
