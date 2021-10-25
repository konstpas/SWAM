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
  public :: get_heat_capacity_tungsten
  public :: get_conductivity_tungsten
  public :: get_atomic_mass_tungsten
  public :: get_melting_point_tungsten 
  public :: get_mass_density_tungsten
  public :: get_electrical_resistivity_tungsten
  public :: get_emissivity_tungsten
  public :: get_enthalpy_of_vaporization_tungsten
  public :: get_Richardson_tungsten
  public :: get_surface_tension_tungsten
  public :: get_thermelectric_power_tungsten
  public :: get_vapor_pressure_tungsten
  public :: get_viscosity_tungsten
  public :: get_work_function_tungsten

contains 

  
  ! -----------------------------------------------------------------
  ! Subroutine used to compute the thermal conductivity
  ! -----------------------------------------------------------------
  subroutine get_conductivity_tungsten(temp,ktherm)

    ! Solid phase tungsten fit is partially adopted from J.G. Hust and A.B. Lankford, Thermal Conductivity of Aluminum,
    ! Copper, Iron and Tungsten for  Temperatures from 1 K to the Melting Point, U. S. Department of Commerce, Boulder,
    ! 1984, pp. 199–256 . NBS Internal Report 84-3007. Their original complicated fit was modified by digitizing the fitting 
    ! function with sampling steps of 50K within 300K-3700K and by least square fitting the emerging dataset to a Shomate type fit.
    ! Liquid phase tungsten data for the electrical resistivity adopted from U. Seydel and W. Fucke, Electrical resistivity of
    ! liquid Ti, V, Mo and W, J. Phys. F 10, L203 (1980). The electrical resistivity data are converted to thermal conductivity 
    ! data with the aid of the Wiedemann-Franz law for a nominal Lorenz number and then fitted to a quadratic polynomial around
    ! the melting point.
    
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
    
  end subroutine get_conductivity_tungsten


  ! -----------------------------------------------------------------
  ! Subroutine used to compute the mass density
  ! -----------------------------------------------------------------
  subroutine get_mass_density_tungsten(temp,rho)

    ! J. Thermophys. 18, 1269 (1997). The polynomial fit for the linear expansion coefficient is translated to a fit for
    ! the normalized linear dimension (after integration), which is translated to a fit for the specific volume (after cubing),
    ! which is translated to a fit for the mass density (after using rho = 19.25g/cm^3 at 300K). Overall, cubic polynomial fit
    ! around the room temperature.
    ! Liquid phase tungsten data adopted from E. Kaschnitz, G. Pottlacher and L. Windholz, High-pressure, high-temperature
    ! thermophysical measurements on tungsten, High Press. Res. 4, 558 (1990). Specific volume curve digitized in steps of 100K,
    ! least-square fitted to a quadratic polynomial and translated to a fit for the mass density (after using rho = 19.25g/cm^3 at 
    ! 300K). Overall, quadratic polynomial fit around the melting temperature.
    
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
   
  end subroutine get_mass_density_tungsten
  
  ! -----------------------------------------------------------------
  ! Subroutine used to compute the heat capacity
  ! -----------------------------------------------------------------
  subroutine get_heat_capacity_tungsten(temp,Cp) 

    ! Solid phase tungsten fit in the range 300 ≤ T_K ≤ 3080 adopted from G.K. White and M.L. Minges, Thermophysical properties of 
    ! some key solids: An update, Int. J. Thermophys. 18, 1269 (1997). Shomate-type fit to the data.
    ! Solid phase tungsten fit in the range 3080 ≤ T_K ≤ 3695 constructed from the specific enthalpy fit of B. Wilthan,
    ! C. Cagran & G. Pottlacher, Combined DSC and Pulse-Heating Measurements of Electrical Resistivity and Enthalpy of Tungsten,
    ! Niobium, and Titanium, Int. J. Thermophys. 26, 1017 (2005). 
    ! The original fit is quadratic leading to a linear fit after differentiation.
    ! Liquid phase tungsten data and fit adopted from B. Wilthan, C. Cagran & G. Pottlacher, Combined DSC and
    ! Pulse-Heating Measurements of  Electrical Resistivity and Enthalpy of Tungsten, Niobium, and Titanium, Int. J. Thermophys.
    ! 26, 1017 (2005). The liquid data are nearly constant.
    
    ! Input and output variables 
    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: Cp     ! Specific heat capacity [J/kgK]

    ! Local variables
    real(amrex_real) :: m_A ! Atomic mass [g/mol]
    
    if (temp.lt.300.0) then
       Cp = 24.1363 
    else if(temp.lt.3080.0) then 
       Cp = 21.868372 + 8.068661E-3*temp - 3.756196E-6*temp**2 + 1.075862E-9*temp**3 + 1.406637E4/(temp**2)
    else if(temp.lt.3695.0) then 
       Cp = 2.022 + 1.315E-2*temp 
    else 
       Cp = 51.3 
    end if
    ! Conversion from J/(mol*K) to J/(kg*K) 
    call get_atomic_mass_tungsten(m_A)
    Cp = 1E3*(Cp/m_A)
    
  end subroutine get_heat_capacity_tungsten

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the electrical resistivity
  ! -----------------------------------------------------------------
  subroutine get_electrical_resistivity_tungsten(temp,rho_e)
    
    ! Solid phase tungsten fit adopted from G.K. White and M.L. Minges, Thermophysical properties of some key solids: An update,
    ! Int. J. Thermophys. 18, 1269 (1997). Fourth-order polynomial fit.
    ! Liquid phase tungsten fit adopted from U. Seydel and W. Fucke, Electrical resistivity of liquid Ti, V, Mo and W, J. Phys.
    ! F 10, L203 (1980). Quadratic fit around the melting point.
    
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
  ! Subroutine used to compute the surface tension
  ! -----------------------------------------------------------------
  subroutine get_surface_tension_tungsten(temp,sigma)

    ! Liquid phase tungsten data and fits adopted P.-F. Paradis, T. Ishikawa, R. Fujii and S. Yoda, Physical properties of liquid
    ! and undercooled tungsten by levitation techniques, Appl. Phys. Lett. 86, 041901 (2005).
    ! Linear fit around the melting point.
    
    real(amrex_real), intent(in) :: temp    ! Temperature [K]
    real(amrex_real), intent(out) :: sigma  ! Surface tension [N/m] 
    
    if (temp.lt.3695.0) then
       sigma = 2.48
    else 
       sigma = 2.48 - 0.31E-3*(temp-3695.0)
    end if
    
  end subroutine get_surface_tension_tungsten


  ! -----------------------------------------------------------------
  ! Subroutine used to compute vapor pressure
  ! -----------------------------------------------------------------
  subroutine get_vapor_pressure_tungsten(temp,pv)

    ! Tungsten fitting parameters adopted from C. L. Yaws, The Yaws handbook of vapor pressure, Elsevier, Oxford, 2018.
    ! Standard Antoine fit.
    
    real(amrex_real), intent(in) :: temp    ! Temperature [K]
    real(amrex_real), intent(out) :: pv     ! Vapor pressure [Pa]

    if(temp.lt.300.0) then
      pv = 0.0
    else
      pv = 10**(12.108 - 40387.965/(temp-141.51))
    end if
    
  end subroutine get_vapor_pressure_tungsten


  ! -----------------------------------------------------------------
  ! Subroutine used to compute the viscosity
  ! -----------------------------------------------------------------
  subroutine get_viscosity_tungsten(temp,mu)

    ! Liquid phase tungsten data and fits adopted from T. Ishikawa, P.-F. Paradis, J.T. Okada, M.V. Kumar and Y. Watanabe, 
    ! Viscosity of molten Mo, Ta, Os, Re and W measured by electrostatic levitation, J. Chem. Thermodyn. 65, 1 (2013).
    ! Arrhenius type fit around the melting point.
    
    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: mu    ! Viscosity [Pa*s] 
    
    if (temp.lt.3695.0) then
       mu = 0.0085;
    else 
       mu = 0.16E-3*EXP(3.9713*3695.0/temp)
    end if
    
  end subroutine get_viscosity_tungsten
  

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the atomic mass
  ! -----------------------------------------------------------------
  subroutine get_atomic_mass_tungsten(m_A)

    real(amrex_real), intent(out) :: m_A ! Atomic mass [g/mol]

    m_A = 183.84
    
  end subroutine get_atomic_mass_tungsten

  
  ! -----------------------------------------------------------------
  ! Subroutine used to compute the melting point properties
  ! -----------------------------------------------------------------
  subroutine get_melting_point_tungsten(temp_melt, enth_fus, rho_melt)

    ! The tungsten enthalpy of fusion has been adopted from the recommendation of P. Tolias, Analytical expressions for thermophysical
    ! properties of solid and liquid tungsten relevant for fusion applications, Nuclear Materials and Energy 13, 42 (2017).
    
    real(amrex_real), intent(out) :: temp_melt ! Temperature at melting [K] 
    real(amrex_real), intent(out) :: enth_fus ! Enthalpy of fusion [kJ/mol] 
    real(amrex_real), intent(out) :: rho_melt ! Density at melting [kg/m3] 

    temp_melt = 3695.0
    enth_fus = 52.3
    rho_melt = 17.1E3
    
  end subroutine get_melting_point_tungsten

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy of vaporization
  ! -----------------------------------------------------------------
  subroutine get_enthalpy_of_vaporization_tungsten(temp, enth_vap)

    ! The tungsten enthalpy of vaporization at the room temperature of 298.15K has been adopted from Arblaster, Thermodynamic
    ! properties of tungsten, J. Phase Equilib. Diffus, 39, 891 (2018). The recommended value is Deltah_v = 855kJ/mol. 
    ! The tungsten enthalpy of vaporization at the normal boiling point of 6200K has been adopted from Zhang, Evans and Yang, 
    ! Corrected values for boiling points and enthalpies of vaporization of elements in handbooks, 
    ! J. Chem. Eng. Data 56, 328-337 (2011). The recommended value is Deltah_v = 774kJ/mol. 
    ! There is a temperature dependence of the enthalpy of vaporization that can be considered to be linear.
    
    real(amrex_real), intent(in) :: temp      ! Temperature [K]
    real(amrex_real), intent(out) :: enth_vap ! Enthalpy of vaporization [kj/mol]

    enth_vap = -0.0137245*temp+859.0920
   
  end subroutine get_enthalpy_of_vaporization_tungsten  


  ! -----------------------------------------------------------------
  ! Subroutine used to compute the work function
  ! -----------------------------------------------------------------
  subroutine get_work_function_tungsten(Wf)

    ! Very good agreement on the recommended work function value between old, classical and contemporary reviews see for
    ! instance V. S. Fomenko, The Handbook of thermionic properties, Plenum, New York (1966); S. Trasatti, Electronegativity,
    ! work function, and heat of adsorption of hydrogen on metals, J. Chem. Soc., Faraday Trans. 1, 68, 229 (1972);
    ! H. B. Michaelson, The work function of the elements and its periodicity, J. Appl. Phys. 48, 4729 (1977); 
    ! H. Kawano, Effective work functions for ionic and electronic emissions from mono- and polycrystalline surfaces,
    ! Prog. Surf. Sci. 83, 1-165 (2008).
    ! For the first very accurate measurements see M. H. Nichols, Average Thermionic Constants of Polycrystalline Tungsten Wires,
    ! Phys. Rev. 78, 158 (1950); B. J. Hopkins and G. C. Riviere, The work function of polycrystalline tungsten foil, 
    ! Proc. Phys. Soc., 81, 590 (1963); R.G Wilson, Vacuum thermionic work functions of polycrystalline Nb, Mo, Ta, W, Re, Os and
    ! Ir, J. Appl. Phys. 37, 3170 (1966);
    
    real(amrex_real), intent(out) :: Wf ! Enthalpy of vaporization [J]
    
    Wf = 4.55
    ! Conversion from eV to J
    Wf = Wf*1.60218E-19
    
  end subroutine get_work_function_tungsten


  ! -----------------------------------------------------------------
  ! Subroutine used to compute Richardson constant
  ! -----------------------------------------------------------------
  subroutine get_Richardson_tungsten(Aeff)

    ! Relatively good agreement on the Richardson constant value, after correcting for small work function discrepancies see for
    ! instance M. H. Nichols, The Thermionic Constants of Tungsten as a Function of Crystallographic Direction,
    ! Phys. Rev. 57, 297 (1940); M. H. Nichols, Average Thermionic Constants of Polycrystalline Tungsten Wires, Phys. Rev.
    ! 78, 158 (1950); V. S. Fomenko, The Handbook of thermionic properties, Plenum, New York (1966);
    ! The recommended value is Aeff = 60A/ cm^2 K^22 that is half the nominal Richardson value of 𝐴𝐴 nom = 120A/ cm 2 K 2 . 
    ! This is common to many refractory metals.
    
    real(amrex_real), intent(out) :: Aeff ! Richardson constant [A/m^2*K^2]
    
    Aeff = 60E4
   
  end subroutine get_Richardson_tungsten


  ! -----------------------------------------------------------------
  ! Subroutine used to compute emissivity
  ! -----------------------------------------------------------------
  subroutine get_emissivity_tungsten(temp, eps_t)

    ! For solid tungsten, a synthetic dataset was constructed featuring 15 values from 2000K up to 3400K
    ! (steps of 100K) adopted from Matsumoto, Cezairliyan and Basak, Hemispherical total emissivity of Niobium, Molybdenum and 
    ! Tungsten at High Temperatures Using a Combined Transient and Brief Steady State Technique, Int. J. Thermophys. 20, 943 (1999)
    ! and featuring 11 values from 300K up to 1300K (steps of 100K) adopted from the ITER material handbook. The dataset was
    ! fitted with a quadratic polynomial around the room temperature.
    ! For liquid tungsten, a constant value was considered by assuming that a small positive emissivity jump of 0.02 occurs at the
    ! melting point (similar to Niobium).
    
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
 subroutine get_thermelectric_power_tungsten(temp, S)

   ! The discontinuity jump at the melting point has not been measured, thus it is impossible to justify our assumption
   ! that thermoelectric effects do not influence the replacement current.
   ! Solid tungsten data (273-1600K) adopted from L. Abadlia, F. Gasser, K. Khalouk, M. Mayoufi, and J. G. Gasser, New experimental
   ! methodology, setup and LabView program for accurate absolute thermoelectric power and electrical resistivity measurements
   ! between 25 and 1600 K: Application to pure copper, platinum, tungsten, and nickel at very high temperatures, Rev. Sci. Instrum.
   ! 85, 095121 (2014). Note that these measurements agree well will those of R. Roberts, F. Righini & R. Compton, Absolute scale of
   ! thermoelectricity III, Philosophical Magazine Part B, 52, 1147-1163 (1985); N. Cusack and P. Kendall, The Absolute Scale of 
   ! Thermoelectric Power at High Temperature, Proc. Phys. Soc. 72, 898 (1958).
   ! Unfortunately, there is a very large extrapolation range (1600-3695K) and due to the presence of a maximum at the end of the
   ! experimental range, it is dangerous to extrapolate the behavior of the curve up to the melting point. In absence of high
   ! temperature data, this is our only resort for the time being.
   ! There are no liquid tungsten measurements and any educated guess is impossible. Existing liquid metal measurements focused 
   ! on low melting point metals, for instance alkaline earths are characterized by a positive small discontinuity but the Seebeck 
   ! slope can either switch sign or not at the melting point, lathanides are characterized by a negative small discontinuity but
   ! the Seebeck slope can either switch sign or not at the melting point, palladium is characterized by a large negative
   ! discontinuity etc etc. Quadratic fit for the solid phase.
   
   real(amrex_real), intent(in) :: temp   ! Temperature [K]
   real(amrex_real), intent(out) :: S     ! Thermoelectric power [V/K]

   if(temp.lt.300.0) then
      S = 1.67403
   elseif(temp.lt.3695.0) then
      S = 1.67403 + 35.2054E-3*(temp-300) - 16.0719E-6*(temp-300)**2
   else
      S = 0 ! No available data, place-holder value
   endif

   ! Conversion from uV/K to V/K
   S = S*1E-6
   
 end subroutine get_thermelectric_power_tungsten

  
 
end module material_properties_tungsten_module
