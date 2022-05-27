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
  public :: get_heat_capacity_iridium
  public :: get_conductivity_iridium
  public :: get_atomic_mass_iridium
  public :: get_melting_point_iridium 
  public :: get_mass_density_iridium
  public :: get_electrical_resistivity_iridium
  public :: get_emissivity_iridium
  public :: get_enthalpy_of_vaporization_iridium
  public :: get_richardson_constant_iridium
  public :: get_surface_tension_iridium
  public :: get_temp_deriv_surface_tension_iridium
  public :: get_thermelectric_power_iridium
  public :: get_vapor_pressure_iridium
  public :: get_viscosity_iridium
  public :: get_work_function_iridium

contains 

  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the thermal conductivity
  ! -----------------------------------------------------------------
  subroutine get_conductivity_iridium(temp,ktherm)

    ! Solid phase Iridium data (300-1900K) adopted from the synthetic / extrapolated dataset of Ho, Powell
    ! and Liley, Thermal conductivity of the elements, J. Phys. Chem. Ref. Data 3, 689–704 (1974).
    ! Solid phase Iridium data (2000-2719K) adopted from Cagran and Pottlacher, Physical properties & normal
    ! spectral emissivity of Ir up to 3500K, 16 th symposium on thermophysical properties, Boulder, USA (2006).
    ! Liquid phase Iridium data adopted from Cagran and Pottlacher, Physical properties & normal spectral
    ! emissivity of Ir up to 3500K, 16 th symposium on thermophysical properties, Boulder, USA (2006).
    ! Shomate fit to the solid data and linear fit to the liquid data.
    
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
    
  end subroutine get_conductivity_iridium


  ! -----------------------------------------------------------------
  ! Subroutine used to computed the mass density
  ! -----------------------------------------------------------------
  subroutine get_mass_density_iridium(temp,rho)

    ! Solid phase Iridium data adopted from J. W. Arblaster, Crystallographic properties of iridium, Platinum
    ! Metals Rev. 54, 93 (2010)
    ! Liquid phase Iridium data adopted from J. W. Arblaster, Selected Values for the Densities and Molar
    ! Volumes of the Liquid Platinum Group Metals, Johnson Matthey Technol. Rev. 61, 80–86 (2017)
    ! Taylor fit around room temperature to the solid data, Taylor fit around melting point to the liquid data.
    
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
   
  end subroutine get_mass_density_iridium
  

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the heat capacity
  ! -----------------------------------------------------------------
  subroutine get_heat_capacity_iridium(temp,Cp) 

    ! Solid phase Iridium data adopted from J. W. Arblaster, The thermodynamic properties of iridium on ITS-
    ! 90, Calphad 19, 365-372 (1995).
    ! Liquid phase Iridium data adopted from Cagran and Pottlacher, Physical properties & normal spectral
    ! emissivity of iridium up to 3500K, 16 th symposium on thermophysical properties, Boulder, USA (2006).
    ! Shomate fit to the solid data, while the liquid data are nearly constant.
    
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
    call get_atomic_mass_iridium(m_A)
    Cp = 1E3*(Cp/m_A)
    
  end subroutine get_heat_capacity_iridium

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the electrical resistivity
  ! -----------------------------------------------------------------
  subroutine get_electrical_resistivity_iridium(temp,rho_e)

    ! Solid and liquid phase Iridium fits adopted from J. W. Arblaster, Selected Electrical Resistivity Values for
    ! the Platinum Group of Metals II:Rhodium and Iridium, Johnson Matthey Technol. Rev. 60, 4–11 (2016).
    ! Cubic fit to the solid data, linear fit to the liquid data.
    
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
  subroutine get_surface_tension_iridium(temp,sigma)

    ! Liquid phase Iridium data and fits adopted from Paradis, Ishikawa and Okada, Thermophysical Properties
    ! of Platinum Group Metals in their Liquid Undercooled and Superheated Phases, Johnson Matthey Technol.
    ! Rev. 58, 124–136 (2014).
    ! Linear fit around the melting point.
    
    real(amrex_real), intent(in) :: temp    ! Temperature [K]
    real(amrex_real), intent(out) :: sigma  ! Surface tension [N/m] 
    
    if (temp.lt.2719.0) then
       sigma = 2.262
    else 
       sigma = 2.262 - 0.28E-3*(temp-2719)
    end if
    
  end subroutine get_surface_tension_iridium


  ! -----------------------------------------------------------------
  ! Subroutine used to computed the surface tension
  ! -----------------------------------------------------------------
  subroutine get_temp_deriv_surface_tension_iridium(temp,dsigma_dT)

   ! Liquid phase Iridium data and fits adopted from Paradis, Ishikawa and Okada, Thermophysical Properties
   ! of Platinum Group Metals in their Liquid Undercooled and Superheated Phases, Johnson Matthey Technol.
   ! Rev. 58, 124–136 (2014).
   ! Linear fit around the melting point.
   
   real(amrex_real), intent(in) :: temp    ! Temperature [K]
   real(amrex_real), intent(out) :: dsigma_dT  ! Surface tension [N/(K*m)] 
   
   if (temp.lt.2719.0) then
      dsigma_dT = 0.0
   else 
      dsigma_dT =  - 0.28E-3
   end if
   
 end subroutine get_temp_deriv_surface_tension_iridium


  ! -----------------------------------------------------------------
  ! Subroutine used to computed vapor pressure
  ! -----------------------------------------------------------------
  subroutine get_vapor_pressure_iridium(temp,pv)

    ! Iridium fitting parameters adopted from C. L. Yaws, The Yaws handbook of vapor pressure, Elsevier,
    ! Oxford, 2018.
    ! Standard Antoine fit.
    
    real(amrex_real), intent(in) :: temp    ! Temperature [K]
    real(amrex_real), intent(out) :: pv     ! Vapor pressure [Pa]
    
    pv = 10**(11.7592 - 30934.6/(temp-82.47))
    
  end subroutine get_vapor_pressure_iridium


 ! -----------------------------------------------------------------
 ! Subroutine used to computed the viscosity
 ! -----------------------------------------------------------------
 subroutine get_viscosity_iridium(temp,mu)

   ! Liquid phase Iridium data and fits adopted from Paradis, Ishikawa and Okada, Thermophysical Properties
   ! of Platinum Group Metals in their Liquid Undercooled and Superheated Phases, Johnson Matthey Technol.
   ! Rev. 58, 124–136 (2014).
   ! Arrhenius type fit around the melting point.
   
   real(amrex_real), intent(in) :: temp   ! Temperature [K]
   real(amrex_real), intent(out) :: mu    ! Viscosity [Pa*s] 
   
   if (temp.lt.2719.0) then
      mu = 0.0085;
   else 
      mu = 0.59E-3*EXP(2.309*2719.0/temp)
   end if
   
 end subroutine get_viscosity_iridium  
  

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the atomic mass
  ! -----------------------------------------------------------------
  subroutine get_atomic_mass_iridium(m_A)

    real(amrex_real), intent(out) :: m_A ! Atomic mass [g/mol]

    m_A = 192.217
    
  end subroutine get_atomic_mass_iridium

  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the melting point properties
  ! -----------------------------------------------------------------
  subroutine get_melting_point_iridium(temp_melt, enth_fus, rho_melt)

    ! The iridium enthalpy of fusion has been adopted from Cagran and Pottlacher, Physical properties & normal spectral emissivity
    ! of iridium up to 3500K, 16 th symposium on thermophysical properties, Boulder, USA (2006).
    
    real(amrex_real), intent(out) :: temp_melt ! Temperature at melting [K] 
    real(amrex_real), intent(out) :: enth_fus ! Enthalpy of fusion [kJ/mol] 
    real(amrex_real), intent(out) :: rho_melt ! Density at melting [kg/m3] 

    temp_melt = 2719.0
    enth_fus = 39.37
    rho_melt = 20.205E3
    
  end subroutine get_melting_point_iridium

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the enthalpy of vaporization
  ! -----------------------------------------------------------------
  subroutine get_enthalpy_of_vaporization_iridium(enth_vap)

    ! The iridium enthalpy of vaporization has been adopted from Zhang, Evans and Yang, Corrected values for boiling points and 
    ! enthalpies of vaporization of elements in handbooks, J. Chem. Eng. Data 56, 328-337 (2011).
    
    real(amrex_real), intent(out) :: enth_vap ! Enthalpy of vaporization [kj/mol]
    
    enth_vap = 564.0
   
  end subroutine get_enthalpy_of_vaporization_iridium  


  ! -----------------------------------------------------------------
  ! Subroutine used to compute the work function
  ! -----------------------------------------------------------------
  subroutine get_work_function_iridium(Wf)

    ! Very good agreement on the work function value between old, classical and contemporary literature
    ! see for instance R.G Wilson, Vacuum thermionic work functions of polycrystalline Nb, Mo, Ta, W, Re,
    ! Os and Ir, J. Appl. Phys. 37, 3170 (1966); V. S. Fomenko, The Handbook of thermionic properties,
    ! Plenum, New York (1966); H. B. Michaelson, The work function of the elements and its periodicity, J.
    ! Appl. Phys. 48, 4729 (1977); H. Kawano, Effective work functions for ionic and electronic emissions
    ! from mono- and polycrystalline surfaces, Prog. Surf. Sci. 83, 1-165 (2008).
    
    real(amrex_real), intent(out) :: Wf ! Enthalpy of vaporization [J]
    
    Wf = 5.3
    
    ! Conversion from eV to J
    Wf = Wf*1.60218E-19
    
  end subroutine get_work_function_iridium 


  ! -----------------------------------------------------------------
  ! Subroutine used to compute Richardson constant
  ! -----------------------------------------------------------------
  subroutine get_richardson_constant_iridium(Aeff)

    ! Relatively good agreement on the Richardson constant value, after correcting for small work function
    ! discrepancies see for instance V. S. Fomenko, The Handbook of thermionic properties, Plenum, New
    ! York (1966); D. L. Goldwater and W. E. Danforth, Thermionic emission constants of Iridium, Phys. Rev.
    ! 103, 871 (1956).
    
    real(amrex_real), intent(out) :: Aeff ! Richardson constant [A/m^2*K^2]
   
    Aeff = 100E4
   
  end subroutine get_richardson_constant_iridium


  ! -----------------------------------------------------------------
  ! Subroutine used to compute emissivity
  ! -----------------------------------------------------------------
  subroutine get_emissivity_iridium(temp, eps_t)
 
    ! Very few available measurements of the total hemispherical emissivity at elevated temperatures. The
    ! sparse spectral emissivity results refer to a very narrow wavelength range that does not cover the
    ! Wien maximum in a satisfactory manner rendering them useless.
    ! For solid iridium, a synthetic dataset has been constructed that features 4 measured values, 1
    ! extrapolated value and 1 interpolated value. The dataset was fitted with a quadratic polynomial.
    ! For liquid iridium, a constant value was assumed.
    
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
 subroutine get_thermelectric_power_iridium(temp, S)

   ! The discontinuity jump at the melting point has not been
   ! measured, thus it is impossible to justify our assumption that thermoelectric effects do not influence
   ! the replacement current.
   ! Solid iridium data (300-1400K) were adopted from Brodowsky, Albus, Chen, Xiao and Yin,
   ! Measurement and Calculation of the Absolute Thermoelectric Power of Rhodium and Iridium, Journal
   ! of Electronic Materials 40, 907 (2011). The measurements are between 100K and 1400K, but only the
   ! values above 300K were selected in order to avoid the low temperature maximum.
   ! There is a large extrapolation range (1400-2719K) but the behavior of the curve (in terms of values
   ! and slopes) remains reasonable up to the melting point.
   ! There are no liquid iridium measurements and any educated guess is impossible. Existing liquid metal
   ! measurements focused on low melting point metals, for instance alkaline earths are characterized by
   ! a positive small discontinuity but the Seebeck slope can either switch sign or not at the melting point,
   ! lathanides are characterized by a negative small discontinuity but the Seebeck slope can either switch
   ! sign or not at the melting point, palladium is characterized by a large negative discontinuity etc etc.
   
   real(amrex_real), intent(in) :: temp   ! Temperature [K]
   real(amrex_real), intent(out) :: S     ! Thermoelectric power [V/K]
   
   if(temp.lt.300.0) then
      S = 0.86738
   elseif(temp.lt.2719.0) then
      S = 0.86738 - 1.9990E-3*(temp-300) - 0.75620E-6*(temp-300)**2
   else
       S = 0 ! No available data, place-holder value
   endif

   ! Conversion from uV/K to V/K
   S = S*1E-6
   
 end subroutine get_thermelectric_power_iridium

  
 
end module material_properties_iridium_module
