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
  public :: get_heat_capacity_niobium
  public :: get_conductivity_niobium
  public :: get_atomic_mass_niobium
  public :: get_melting_point_niobium 
  public :: get_mass_density_niobium
  public :: get_electrical_resistivity_niobium
  public :: get_emissivity_niobium
  public :: get_richardson_constant_niobium
  public :: get_surface_tension_niobium
  public :: get_temp_deriv_surface_tension_niobium
  public :: get_viscosity_niobium
  public :: get_work_function_niobium
  public :: get_enthalpy_of_vaporization_niobium
  public :: get_vapor_pressure_niobium
  public :: get_thermelectric_power_niobium

contains 

  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the thermal conductivity
  ! -----------------------------------------------------------------
  subroutine get_conductivity_niobium(temp,ktherm)

    ! Solid phase niobium data (300-2000K) adopted from the synthetic dataset of Ho, Powell and Liley,
    ! Thermal conductivity of the elements, J. Phys. Chem. Ref. Data 3, 689–704 (1974).
    ! Solid phase niobium data (2000-2745K) generated from the electrical resistivity expression of Hüpf,
    ! Cagran, Lohöfer and Pottlacher, Electrical resistivity of high melting metals up into the liquid phase (V, Nb,
    ! Ta, Mo, W), J. Phys.: Conf. Ser 98, 062002 (2008). The construction was based on the Wiedemann-Franz
    ! law and a re-scaling of the high temperature data to ensure continuity at 2000K.
    ! Liquid phase niobium data generated from the electrical resistivity expression of Hüpf, Cagran, Lohöfer
    ! and Pottlacher, Electrical resistivity of high melting metals up into the liquid phase (V, Nb, Ta, Mo, W), J.
    ! Phys.: Conf. Ser 98, 062002 (2008). The construction was based on the Wiedemann-Franz law.
    ! Shomate fits to the solid and the liquid data.
    
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
    
  end subroutine get_conductivity_niobium


  ! -----------------------------------------------------------------
  ! Subroutine used to computed the mass density
  ! -----------------------------------------------------------------
  subroutine get_mass_density_niobium(temp,rho)

    ! Solid phase niobium data adopted from J. W. Arblaster, Selected Values of the Crystallographic Properties
    ! of Elements, ASM International, USA, 2018.
    ! Liquid phase niobium fit adopted from Paradis, Ishikawa, Lee, Holland-Moritz, Brillo, Rhim and Okada,
    ! Materials properties measurements and particle beam interactions studies using electrostatic levitation,
    ! Materials Science and Engineering R 76, 1–53 (2014).
    ! Taylor fit around room temperature to the solid data, Taylor fit around melting point to the liquid data.
    
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
   
  end subroutine get_mass_density_niobium
  

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the heat capacity
  ! -----------------------------------------------------------------
  subroutine get_heat_capacity_niobium(temp,Cp) 

    ! Solid phase niobium data adopted from J. W. Arblaster, The thermodynamic properties of niobium, J.
    ! Phase Equilib. Diffus. 38, 707-722 (2017).
    ! Liquid phase niobium data adopted from Wilthan, Cagran and Pottlacher, Combined DSC and Pulse-
    ! Heating Measurements of Electrical Resistivity and Enthalpy of Tungsten, Niobium and Titanium, Int. J.
    ! Thermophys. 26, 1017 (2005).
    ! Shomate fit to the solid data, while the liquid data are nearly constant.
    
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
    call get_atomic_mass_niobium(m_A)
    Cp = 1E3*(Cp/m_A)
    
  end subroutine get_heat_capacity_niobium
  

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the electrical resistivity
  ! -----------------------------------------------------------------
  subroutine get_electrical_resistivity_niobium(temp,rho_e)

    ! Solid phase niobium fits adopted from Maglic, Perovic, Vukovic and Zekovic, Specific Heat and Electrical
    ! Resistivity of Niobium Measured by Subsecond Calorimetric Technique, Int. J. Thermophys. 15, 963 (1994).
    ! Thermal expansion effects were then included, a new dataset was generated and curve fitted.
    ! Liquid phase niobium fits adopted from Hüpf, Cagran, Lohöfer and Pottlacher, Electrical resistivity of high
    ! melting metals up into the liquid phase (V, Nb, Ta, Mo, W), J. Phys.: Conf. Ser 98, 062002 (2008).
    ! Cubic fit to the solid data, linear fit to the liquid data.

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
  subroutine get_surface_tension_niobium(temp,sigma)

    ! Liquid phase niobium data and fits adopted from Ishikawa, Paradis, Okada and Watanabe, Viscosity
    ! measurements of molten refractory metals using an electrostatic levitator, Meas. Sci. Technol. 23, 025305
    ! (2012). Linear fit around the melting point
    
    real(amrex_real), intent(in) :: temp    ! Temperature [K]
    real(amrex_real), intent(out) :: sigma  ! Surface tension [N/m] 
    
    if (temp.lt.2745.0) then
       sigma = 1.967
    else 
       sigma = 1.967 - 0.17E-3*(temp-2745.0)
    end if
    
  end subroutine get_surface_tension_niobium

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the derivative of the surface tension
  ! with respect to the temperature.
  ! -----------------------------------------------------------------
  subroutine get_temp_deriv_surface_tension_niobium(temp,dsigma_dT)

    ! Liquid phase niobium data and fits adopted from Ishikawa, Paradis, Okada and Watanabe, Viscosity
    ! measurements of molten refractory metals using an electrostatic levitator, Meas. Sci. Technol. 23, 025305
    ! (2012). Linear fit around the melting point
    
    real(amrex_real), intent(in) :: temp    ! Temperature [K]
    real(amrex_real), intent(out) :: dsigma_dT  ! Surface tension [N/(K*m)] 
    
    if (temp.lt.2745.0) then
       dsigma_dT = 0.0
    else 
       dsigma_dT = - 0.17E-3
    end if
    
  end subroutine get_temp_deriv_surface_tension_niobium


  ! -----------------------------------------------------------------
  ! Subroutine used to computed vapor pressure
  ! -----------------------------------------------------------------
  subroutine get_vapor_pressure_niobium(temp,pv)

    ! Niobium fitting parameters adopted from C. L. Yaws, The Yaws handbook of vapor pressure,
    ! Elsevier, Oxford, 2018.
    ! Standard Antoine fit.
    
    real(amrex_real), intent(in) :: temp    ! Temperature [K]
    real(amrex_real), intent(out) :: pv     ! Vapor pressure [Pa]
    
    pv = 10**(12.0502-35189.3/(temp-21.93))
    
  end subroutine get_vapor_pressure_niobium  
  

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the viscosity
  ! -----------------------------------------------------------------
  subroutine get_viscosity_niobium(temp,mu)

    ! Liquid phase niobium data and fits adopted from Ishikawa, Paradis, Okada and Watanabe, Viscosity
    ! measurements of molten refractory metals using an electrostatic levitator, Meas. Sci. Technol. 23, 025305
    ! (2012).
    ! Arrhenius type fit around the melting point.
    
    real(amrex_real), intent(in) :: temp    ! Temperature [K]
    real(amrex_real), intent(out) :: mu     ! Viscosity [Pa*s] 
    
    if (temp.lt.2745.0) then
       mu = 0.0050;
    else 
      mu = 1.05E-3*EXP(1.556*2745.0/temp)
    end if
    
  end subroutine get_viscosity_niobium  
  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the atomic mass
  ! -----------------------------------------------------------------
  subroutine get_atomic_mass_niobium(m_A)

    real(amrex_real), intent(out) :: m_A ! Atomic mass [g/mol]

    m_A = 92.90638
    
  end subroutine get_atomic_mass_niobium

  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the melting point properties
  ! -----------------------------------------------------------------
  subroutine get_melting_point_niobium(temp_melt, enth_fus, rho_melt)

    ! The niobium enthalpy of fusion has been adopted from Wilthan, Cagran and Pottlacher,
    ! Combined DSC and Pulse-Heating Measurements of Electrical Resistivity and Enthalpy of Tungsten,
    ! Niobium and Titanium, Int. J. Thermophys. 26, 1017 (2005).
    
    real(amrex_real), intent(out) :: temp_melt ! Temperature at melting [K] 
    real(amrex_real), intent(out) :: enth_fus   ! Enthalpy of fusion [kJ/mol]
    real(amrex_real), intent(out) :: rho_melt   ! Density at metling [Kg/m^3]

    temp_melt = 2745.0
    enth_fus = 32.98
    rho_melt = 7862.5
    
  end subroutine get_melting_point_niobium


  ! -----------------------------------------------------------------
  ! Subroutine used to computed the enthalpy of vaporization
  ! -----------------------------------------------------------------
  subroutine get_enthalpy_of_vaporization_niobium(enth_vap)

    ! The niobium enthalpy of vaporization has been adopted from Zhang, Evans and Yang,
    ! Corrected values for boiling points and enthalpies of vaporization of elements in
    ! handbooks, J. Chem. Eng. Data 56, 328-337 (2011).
    
    real(amrex_real), intent(out) :: enth_vap ! Enthalpy of vaporization [kj/mol]
    
    enth_vap = 694.0
    
  end subroutine get_enthalpy_of_vaporization_niobium  


  ! -----------------------------------------------------------------
  ! Subroutine used to compute the work function
  ! -----------------------------------------------------------------
  subroutine get_work_function_niobium(Wf)

    ! Very good agreement on the work function value between old, classical and contemporary literature
    ! see for instance R.G Wilson, Vacuum thermionic work functions of polycrystalline Nb, Mo, Ta, W, Re,
    ! Os and Ir, J. Appl. Phys. 37, 3170 (1966); S. Trasatti, Electronegativity, Work Function and Heat of
    ! Adsorption of Hydrogen on Metals, Chim. Ind. (Milan) 53, 559 (1971); H. B. Michaelson, The work
    ! function of the elements and its periodicity, J. Appl. Phys. 48, 4729 (1977); H. Kawano, Effective work
    ! functions for ionic and electronic emissions from mono- and polycrystalline surfaces, Prog. Surf. Sci.
    ! 83, 1-165 (2008). The recommended values vary between 4.1 − 4.3eV.
    ! The only exception is V. S. Fomenko, The Handbook of thermionic properties, Plenum, New York
    ! (1966) who recommends the rather low value of 4.0 eV.
    
    real(amrex_real), intent(out) :: Wf ! Enthalpy of vaporization [J]
    
    Wf = 4.2
    
    ! Conversion from eV to J
    Wf = Wf*1.60218E-19
    
  end subroutine get_work_function_niobium    


  ! -----------------------------------------------------------------
  ! Subroutine used to compute Richardson constant
  ! -----------------------------------------------------------------
  subroutine get_richardson_constant_niobium(Aeff)

    ! Relatively good agreement on the Richardson constant value, see V. S. Fomenko, The Handbook of
    ! thermionic properties, Plenum, New York (1966); W. C. Niehaus and E. A. Coomes, Surface-barrier
    ! analysis for Nb, Ta and Ta-on-Nb from periodic deviations in the thermionic Schottky effect, Surf. Sci.
    ! 27, 256 (1971).
    
    real(amrex_real), intent(out) :: Aeff ! Richardson constant [A/m^2*K^2]

    Aeff = 50E4
    
  end subroutine get_richardson_constant_niobium


  ! -----------------------------------------------------------------
  ! Subroutine used to compute emissivity
  ! -----------------------------------------------------------------
  subroutine get_emissivity_niobium(temp, eps_t)

    ! Solid phase niobium data (300-1000K) adopted from S. Cheng, P. Cebe, L. Hanssen, D. Riffe, A. Sievers,
    ! Hemispherical Emissivity of V, Nb, Ta, Mo, and W from 300 to 1000 K, J. Opt. Soc. Am. 68, 351 (1987).
    ! The data were digitized from a figure and fitted with a cubic polynomial.
    ! Solid phase niobium data (1100-2745K) adopted from Maglic, Perovic, Vukovic and Zekovic, Specific
    ! Heat and Electrical Resistivity of Niobium Measured by Subsecond Calorimetric Technique, Int. J.
    ! Thermophys. 15, 963 (1994). The data were tabulated and uniformly re-scaled in order to ensure
    ! continuity with the extrapolated low temperature data at 1100K. The two datasets were then merged
    ! and fitted with a quadratic polynomial.
    ! Liquid phase niobium data adopted from Ishikawa, Okada, Paradis, Watanabe, Emissivity Measurements
    ! of Molten Metals with an Electrostatic Levitator, Int. J. Microgravity Sci. Appl. 34, 340305 (2017). The
    ! value is nearly constant.
    
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
  subroutine get_thermelectric_power_niobium(temp, S)

    ! The discontinuity jump at the melting point has not been
    ! measured, thus it is impossible to justify our assumption that thermoelectric effects do not influence
    ! the replacement current.
    ! Solid niobium data (300-1700K) were adopted from Brodowsky, Chen, Xiao and Yin, Calculation of the
    ! Thermoelectric Power of Vanadium, Niobium and Tantalum, Journal of Electronic Materials 40, 1984
    ! (2011). The original curve lies between 0K and 1700K, but only the values above 300K were selected
    ! due to the low temperature uncertainties.
    ! There is a large extrapolation range (1700-2745K) where the behavior of the fitting curve becomes
    ! unreasonable (due to the very low values acquired) in spite of remaining monotonic. As a result, an
    ! extrapolation was attempted based on continuity at 1700K, first derivative continuity at 1700K and
    ! reasonable value at the melting point.
    ! There are no liquid niobium measurements and any educated guess is impossible. Existing liquid metal
    ! measurements focused on low melting point metals, for instance alkaline earths are characterized by a
    ! positive small discontinuity but the Seebeck slope can either switch sign or not at the melting point,
    ! lathanides are characterized by a negative small discontinuity but the Seebeck slope can either switch
    ! sign or not at the melting point, palladium is characterized by a large negative discontinuity etc etc.
    ! Shomate type fit and quadratic fit for the solid phase.
    
    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: S     ! Thermoelectric power [V/K]

    if(temp.lt.300.0) then
      S = -1.2462
    elseif(temp.le.1700.0) then
      S = 88.280 + 2.69689E6/temp**2 - 24.97805E3/temp - 155.698E-3*temp + &
         126.954E-6*temp**2 - 35.10907E-9*temp**3
    elseif(temp.lt.2745.0) then
      S = 65.199852 - 50.81273E-3*temp + 8.796398E-6*temp**2
    else
      S = 0 ! No available data, place-holder value
    endif

    ! Convversion from uV/K to V/K
    S = S*1E-6
    
  end subroutine get_thermelectric_power_niobium

  
end module material_properties_niobium_module
