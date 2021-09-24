module material_properties_tungsten_module  

  use amrex_amr_module
  
  implicit none 

  ! -----------------------------------------------------------------
  ! This module is used to compute the material properties of
  ! tungsten as described in P. Tolias, Nucl. Mater. Energy 13, 42
  ! (2017)
  ! -----------------------------------------------------------------
  
  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: get_Cp_tungsten
  public :: get_ktherm_tungsten
  public :: get_m_A_tungsten
  public :: get_melting_point_tungsten 
  public :: get_rho_tungsten
  public :: get_electrical_resistivity_tungsten
  public :: get_emissivity_tungsten
  public :: get_enthalpy_of_vap_tungsten
  public :: get_Richardson_tungsten
  public :: get_surf_tension_tungsten
  public :: get_thermelec_power_tungsten
  public :: get_vapor_pressure_tungsten
  public :: get_viscosity_tungsten
  public :: get_work_function_tungsten

contains 

  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the thermal conductivity
  ! -----------------------------------------------------------------
  subroutine get_ktherm_tungsten(temp,ktherm)

    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: ktherm ! Thermal conductivity [W/mK] 
    
    if (temp.lt.300.0) then
       ktherm = 179.9041  
    else if(temp.lt.3695.0) then 
       ktherm = 149.441 - 45.466E-3*temp + 13.193E-6*temp**2 - &
            1.484E-9*temp**3 + 3.866E6/temp**2 
    else  
       ktherm = 66.6212 + 0.02086*(temp-3695.0) - &
            3.7585E-6*(temp-3695.0)**2
    end if
    
  end subroutine get_ktherm_tungsten


  ! -----------------------------------------------------------------
  ! Subroutine used to computed the mass density
  ! -----------------------------------------------------------------
  subroutine get_rho_tungsten(temp,rho)

    ! Input and output variables 
    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: rho    ! Mass density [kg/m3]

    ! Local variables
    real(amrex_real) :: temp_0 = 293.15
    
    if (temp.lt.300.0) then
       rho = 19.2482 
    else if(temp.lt.3695.0) then 
       rho = 19.25 - 2.66207E-4*(temp-temp_0)-3.0595E-9*(temp-temp_0)**2 - &
             9.5185E-12*(temp-temp_0)**3
    else
       rho = 16.267 - 7.679E-4*(temp-3695.0) - 8.091E-8*(temp-3695.0)**2 
    end if
    ! Conversion from g/cm3 to kg/m3 
    rho = rho*1E3  
   
  end subroutine get_rho_tungsten
  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the heat capacity
  ! -----------------------------------------------------------------
  subroutine get_Cp_tungsten(temp,Cp) 

    ! Input and output variables 
    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: Cp     ! Specific heat capacity [J/kgK]

    ! Local variables
    real(amrex_real) :: m_A ! Atomic mass [g/mol]
    
    if (temp.lt.300.0) then
       Cp = 24.1363 
    else if(temp.lt.3080.0) then 
       Cp = 21.868372 + 8.068661E-3*temp - 3.756196E-6*temp**2 + 1.075862E-9*temp**3 + 1.406637E4/(temp**2)
    else if(temp.le.3695.0) then 
       Cp = 2.022 + 1.315E-2*temp 
    else 
       Cp = 51.3 
    end if
    ! Conversion from J/(mol*K) to J/(kg*K) 
    call get_m_A_tungsten(m_A)
    Cp = 1E3*(Cp/m_A)
    
  end subroutine get_Cp_tungsten

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the electrical resistivity
  ! -----------------------------------------------------------------
  subroutine get_electrical_resistivity_tungsten(temp,rho_e)

   ! Input and output variables 
   real(amrex_real), intent(in) :: temp      ! Temperature [K]
   real(amrex_real), intent(out) :: rho_e    ! Mass density [Ohm*m]

   ! Local variables
   real(amrex_real) :: temp_m = 3695.0
   
   if (temp.lt.300.0) then
      rho_e = 5.4702 
   else if(temp.lt.temp_m) then 
      rho_e = -0.9680 + 1.9274E-2*temp + 7.8260E-6*temp**2 - 1.8517E-9*temp**3 &
               + 2.0790E-13*temp**4
   else 
      rho_e = 135 - 1.855E-3*(temp-temp_m) + 4.420E-6*(temp-temp_m)**2
   end if
   ! Conversion from uOhm*cm to Ohm*m 
   rho_e = rho_e*1E-8  
  
  end subroutine get_electrical_resistivity_tungsten


  ! -----------------------------------------------------------------
  ! Subroutine used to computed the surface tension
  ! -----------------------------------------------------------------
  subroutine get_surf_tension_tungsten(temp,sigma)

   real(amrex_real), intent(in) :: temp    ! Temperature [K]
   real(amrex_real), intent(out) :: sigma  ! Surface tension [N/m] 
   
   if (temp.lt.3695.0) then
      sigma = 2.48
   else 
      sigma = 2.48 - 0.31E-3*(temp-3695.0)
   end if
   
  end subroutine get_surf_tension_tungsten


  ! -----------------------------------------------------------------
  ! Subroutine used to computed vapor pressure
  ! -----------------------------------------------------------------
  subroutine get_vapor_pressure_tungsten(temp,pv)

   real(amrex_real), intent(in) :: temp    ! Temperature [K]
   real(amrex_real), intent(out) :: pv     ! Vapor pressure [Pa]
   
   pv = 10**(12.108 - 40387.965/(temp-141.51))
   
  end subroutine get_vapor_pressure_tungsten  


 ! -----------------------------------------------------------------
 ! Subroutine used to computed the viscosity
 ! -----------------------------------------------------------------
 subroutine get_viscosity_tungsten(temp,mu)

   real(amrex_real), intent(in) :: temp   ! Temperature [K]
   real(amrex_real), intent(out) :: mu    ! Viscosity [N/m] 
   
   if (temp.lt.3695.0) then
      mu = 0.0085;
   else 
      mu = 0.16E-3*EXP(3.9713*3695.0/temp)
   end if
   
 end subroutine get_viscosity_tungsten  
  

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the atomic mass
  ! -----------------------------------------------------------------
  subroutine get_m_A_tungsten(m_A)

    real(amrex_real), intent(out) :: m_A ! Atomic mass [g/mol]

    m_A = 183.84
    
  end subroutine get_m_A_tungsten

  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the melting point properties
  ! -----------------------------------------------------------------
  subroutine get_melting_point_tungsten(temp_melt, enth_fus, rho_melt)

    real(amrex_real), intent(out) :: temp_melt ! Temperature at melting [K] 
    real(amrex_real), intent(out) :: enth_fus ! Enthalpy of fusion [kJ/mol] 
    real(amrex_real), intent(out) :: rho_melt ! Density at melting [kg/m3] 

    temp_melt = 3695.0
    enth_fus = 52.3
    rho_melt = 17.1E3
    
  end subroutine get_melting_point_tungsten

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the enthalpy of vaporization
  ! -----------------------------------------------------------------
  subroutine get_enthalpy_of_vap_tungsten(temp, enth_vap)

   real(amrex_real), intent(in) :: temp      ! Temperature [K]
   real(amrex_real), intent(out) :: enth_vap ! Enthalpy of vaporization [kj/mol]

   enth_vap = -0.0137245*temp+859.0920
   
  end subroutine get_enthalpy_of_vap_tungsten  


  ! -----------------------------------------------------------------
  ! Subroutine used to compute the work function
  ! -----------------------------------------------------------------
  subroutine get_work_function_tungsten(Wf)

   real(amrex_real), intent(out) :: Wf ! Enthalpy of vaporization [J]

   Wf = 4.55
   ! Conversion from eV to J
   Wf = Wf*1.60218E-19
   
  end subroutine get_work_function_tungsten 


  ! -----------------------------------------------------------------
  ! Subroutine used to compute Richardson constant
  ! -----------------------------------------------------------------
  subroutine get_Richardson_tungsten(Aeff)

   real(amrex_real), intent(out) :: Aeff ! Richardson constant [A/m^2*K^2]

   Aeff = 60E4
   
  end subroutine get_Richardson_tungsten


  ! -----------------------------------------------------------------
  ! Subroutine used to compute emissivity
  ! -----------------------------------------------------------------
  subroutine get_emissivity_tungsten(temp, eps_t)

   real(amrex_real), intent(in) :: temp   ! Temperature [K]
   real(amrex_real), intent(out) :: eps_t ! Emissivity [dimensionless]

   if (temp.lt.300) then
      eps_t = 1.85076E-2 
   else if(temp.lt.3695) then
      eps_t = 1.85076E-2 + 1.6048E-4*(temp-300) - 1.762608E-8*(temp-300)**2
   else
      eps_t = 0.38
   endif
   
 end subroutine get_emissivity_tungsten


  ! -----------------------------------------------------------------
  ! Subroutine used to compute the absolute thermoelectric power
  ! -----------------------------------------------------------------
 subroutine get_thermelec_power_tungsten(temp, S)

   real(amrex_real), intent(in) :: temp   ! Temperature [K]
   real(amrex_real), intent(out) :: S     ! Thermoelectric power [V/K]

   if(temp.lt.300.0) then
      S = 1.67403
   elseif(temp.le.3695.0) then
      S = 1.67403 + 35.2054E-3*(temp-300) - 16.0719E-6*(temp-300)**2
   else
      S = 0 ! No available data, place-holder value
   endif

   ! Conversion from uV/K to V/K
   S = S*1E-6
   
 end subroutine get_thermelec_power_tungsten

  
 
end module material_properties_tungsten_module
