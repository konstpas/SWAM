module material_properties_beryllium_module  

  use amrex_amr_module
  
  implicit none 

  ! -----------------------------------------------------------------
  ! This module is used to compute the material properties of
  ! beryllium.
  ! -----------------------------------------------------------------
  
  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: get_heat_capacity_beryllium
  public :: get_conductivity_beryllium
  public :: get_atomic_mass_beryllium
  public :: get_melting_point_beryllium 
  public :: get_mass_density_beryllium
  public :: get_electrical_resistivity_beryllium
  public :: get_emissivity_beryllium
  public :: get_richardson_constant_beryllium
  public :: get_surface_tension_beryllium
  public :: get_temp_deriv_surface_tension_beryllium
  public :: get_viscosity_beryllium
  public :: get_work_function_beryllium
  public :: get_enthalpy_of_vaporization_beryllium
  ! public :: get_hcp_to_bcc_point_beryllium
  public :: get_vapor_pressure_beryllium
  public :: get_thermelectric_power_beryllium

contains 

  
  ! -----------------------------------------------------------------
  ! Subroutine used to compute the thermal conductivity
  ! -----------------------------------------------------------------
  subroutine get_conductivity_beryllium(temp,ktherm)

    ! Solid beryllium data (300-1400K) adopted from the synthetic dataset of Ho, Powell and Liley, Thermal
    ! conductivity of the elements, J. Phys. Chem. Ref. Data 3, 689–704 (1974).
    ! In absence of liquid beryllium data, the following procedure was adopted. (a) The high temperature
    ! Lorenz number was extracted from the recommended solid Be descriptions of electrical resistivity and
    ! thermal conductivity. (b) The Wiedemann-Franz law was employed with the extrapolated Lorenz number
    ! and the recommended liquid Be description of electrical resistivity. (c) The linear curve us expressed as a
    ! Taylor fit around the melting point. (d) The slope is slightly reduced to account for resistivity variations.
    ! Shomate fit to the solid data and linear fit to the liquid data.
    
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
    
  end subroutine get_conductivity_beryllium


  ! -----------------------------------------------------------------
  ! Subroutine used to compute the mass density
  ! -----------------------------------------------------------------
  subroutine get_mass_density_beryllium(temp,rho)

    ! Solid phase beryllium data adopted from Y. S. Touloukian, Thermophysical properties of matter – The TPRC
    ! Data Series - Vol 12 Thermal expansion - Metallic elements and alloys (Plenum, New York, 1975).
    ! The tabulated spatial expansion data were translated into volume expansion data (cubed) and then into
    ! mass density data (rho_m = 1.850g/cm 3 at 300K).
    ! Liquid phase beryllium fit adopted from Thermophysical Properties of Materials For Nuclear Engineering: A
    ! Tutorial and Collection of Data (IAEA, Vienna, 2008). Nearly the same fit was recommended in Iida and
    ! Guthrie, The thermophysical properties of metallic liquids Vol. 2 (Oxford University Press, Oxford, 2015).
    ! Taylor fit around room temperature to the solid data, Taylor fit around melting point to the liquid data.

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
   
  end subroutine get_mass_density_beryllium
  

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the heat capacity
  ! -----------------------------------------------------------------
  subroutine get_heat_capacity_beryllium(temp,Cp) 

    ! Solid alpha-phase, solid beta-phase & liquid phase beryllium fits adopted from J. W. Arblaster, Thermodynamic 
    ! Properties of Beryllium, Journal of Phase Equilibria and Diffusion 37, 581-591 (2016).
    ! Shomate fit to the solid α-phase data, constant value for the solid β-phase data and linear fit for the
    ! liquid phase data.
    
    ! Input and output variables 
    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: Cp     ! Specific heat capacity [J/kgK]

    ! Local variables
    real(amrex_real) :: m_A ! Atomic mass [g/mol]
    
    if (temp.lt.300.0) then
       Cp = 16.5438629378 
    ! else if(temp.lt.1543.0) then 
    !    Cp = 21.5390 + 4.94572E-3*temp + 1.356324E-6*temp**2 - 0.594083E6/temp**2
    ! else if(temp.lt.1560) then
       !    Cp = 30.0
    else if(temp.lt.1560.0) then 
       Cp = 21.5390 + 4.94572E-3*temp + 1.356324E-6*temp**2 - 0.594083E6/temp**2   
    else 
       Cp = 25.4345+2.15008E-3*temp 
    end if
    ! Conversion from J/(mol*K) to J/(kg*K) 
    call get_atomic_mass_beryllium(m_A)
    Cp = 1E3*(Cp/m_A)
    
  end subroutine get_heat_capacity_beryllium
  

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the electrical resistivity
  ! -----------------------------------------------------------------
  subroutine get_electrical_resistivity_beryllium(temp,rho_e)

    ! Solid beryllium data adopted from T. C. Chi, "Electrical resistivity of alkaline earth elements", Journal of
    ! Physical and Chemical Reference Data 8, 439 (1979). The same data were suggested in Gmelin, Beryllium,
    ! A3, Supplement.
    ! Liquid beryllium value at melting point adopted from Iida and Guthrie, The thermophysical properties of
    ! metallic liquids Vol. 2 (Oxford University Press, Oxford, 2015). In absence of other reliable data, it is
    ! assumed that the value remains constant based on the small resistivity changes for most liquid metals.
    ! Taylor fit around room temperature to the solid data, constant value for the liquid data.
    
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
  ! Subroutine used to compute the surface tension
  ! -----------------------------------------------------------------
  subroutine get_surface_tension_beryllium(temp,sigma)

    ! Liquid phase beryllium data adopted from Thermophysical Properties of Materials For Nuclear
    ! Engineering: A Tutorial and Collection of Data (IAEA, Vienna, 2008). Only two data points available, one
    ! also reported in Iida and Guthrie, The thermophysical properties of metallic liquids Vol. 2 (Oxford
    ! University Press, Oxford, 2015).
    ! Linear fit around the melting point. This fit leads to a prediction of 7275K for the critical point which is
    ! relatively close to the 8080K prediction of Fortov and collaborators as quoted by Apfelbaum, Estimate of
    ! Beryllium Critical Point on the Basis of Correspondence between the Critical and Zeno-Line Parameters, J.
    ! Phys. Chem. B 116, 14660−14666 (2012).
    
    real(amrex_real), intent(in) :: temp    ! Temperature [K]
    real(amrex_real), intent(out) :: sigma  ! Surface tension [N/m] 
    
    if (temp.lt.1560.0) then
       sigma = 1.143
    else 
       sigma = 1.143 - 0.20E-3*(temp-1560.0)
    end if
    
  end subroutine get_surface_tension_beryllium

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the derivative of the surface tension
  ! with respect to the temperature.
  ! -----------------------------------------------------------------
  subroutine get_temp_deriv_surface_tension_beryllium(temp,dsigma_dT)

   ! Liquid phase beryllium data adopted from Thermophysical Properties of Materials For Nuclear
   ! Engineering: A Tutorial and Collection of Data (IAEA, Vienna, 2008). Only two data points available, one
   ! also reported in Iida and Guthrie, The thermophysical properties of metallic liquids Vol. 2 (Oxford
   ! University Press, Oxford, 2015).
   ! Linear fit around the melting point. This fit leads to a prediction of 7275K for the critical point which is
   ! relatively close to the 8080K prediction of Fortov and collaborators as quoted by Apfelbaum, Estimate of
   ! Beryllium Critical Point on the Basis of Correspondence between the Critical and Zeno-Line Parameters, J.
   ! Phys. Chem. B 116, 14660−14666 (2012).
   
   real(amrex_real), intent(in) :: temp    ! Temperature [K]
   real(amrex_real), intent(out) :: dsigma_dT  ! Surface tension [N/(K*m)] 
   
   if (temp.lt.1560.0) then
       dsigma_dT = 0.0
   else 
       dsigma_dT = - 0.20E-3
   end if
   
 end subroutine get_temp_deriv_surface_tension_beryllium


  ! -----------------------------------------------------------------
  ! Subroutine used to compute vapor pressure
  ! -----------------------------------------------------------------
  subroutine get_vapor_pressure_beryllium(temp,pv)

    ! Beryllium fitting parameters adopted from C. L. Yaws, The Yaws handbook of vapor pressure,
    ! Elsevier, Oxford, 2018. Standard Antoine fit.
    ! A much more complicated fit with branches in the solid α-, solid β-, liquid phases is proposed by J.
    ! W. Arblaster, Thermodynamic Properties of Beryllium, Journal of Phase Equilibria and Diffusion 37,
    ! 581-591 (2016). Except from very low temperatures, where vaporization is anyways negligible, the
    ! two fits are strongly overlapping. Therefore, there is no reason to utilize Arblaster’s fit.
    
    real(amrex_real), intent(in) :: temp    ! Temperature [K]
    real(amrex_real), intent(out) :: pv     ! Vapor pressure [Pa]

    if(temp.lt.300.0) then
      pv = 0.0
    else
      pv = 10**(10.2089-13696.6/(temp-124.63))
    end if
    
  end subroutine get_vapor_pressure_beryllium  
  

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the viscosity
  ! -----------------------------------------------------------------
  subroutine get_viscosity_beryllium(temp,mu)

    ! Liquid phase beryllium data were extracted from a figure provided in Dombrowski et al, Atomic and
    ! Plasma-Material Interaction Data for Fusion vol .5 (IAEA, Vienna, 1995) with the aid of software. The data
    ! exhibited a very pronounced disagreement with both the Andrade formula (at the melting point) and the
    ! Fowler-Born-Green formula (at the entire liquid range). The disagreement exceeded an order of
    ! magnitude and the data were dismissed.
    ! The liquid phase beryllium fit was adopted from L. Battezzati and A. L. Greer, The viscosity of liquid metals
    ! and alloys, Acta Metall. 37, 1791-1802 (1989). It exhibits good agreement with both the Andrade and the
    ! Fowler–Born–Green formulas.
    ! Arrhenius type fit around the melting point.
    
    real(amrex_real), intent(in) :: temp ! Temperature [K]
    real(amrex_real), intent(out) :: mu  ! Viscosity [Pa*s] 
    
    if (temp.lt.1560.0) then
       mu = 0.005090697766920;
    else 
       mu = 0.1E-3*EXP(3.93*1560.0/temp)
    end if
    
  end subroutine get_viscosity_beryllium  
  
  ! -----------------------------------------------------------------
  ! Subroutine used to compute the atomic mass
  ! -----------------------------------------------------------------
  subroutine get_atomic_mass_beryllium(m_A)

    real(amrex_real), intent(out) :: m_A ! Atomic mass [g/mol]

    m_A = 9.0121831
    
  end subroutine get_atomic_mass_beryllium

  
  ! -----------------------------------------------------------------
  ! Subroutine used to compute the melting point properties
  ! -----------------------------------------------------------------
  subroutine get_melting_point_beryllium(temp_melt, enth_fus, rho_melt)

    ! The beryllium enthalpy of fusion has been adopted from J. W. Arblaster, Thermodynamic Properties of Beryllium, Journal of
    ! Phase Equilibria and Diffusion 37, 581-591 (2016). A rather close value is provided by H. Kleykamp, Selected thermal
    ! properties of beryllium and phase equilibria in beryllium systems relevant for nuclear fusion reactor blankets, 
    ! Journal of Nuclear Materials 294, 88-93 (2001).
    
    real(amrex_real), intent(out) :: temp_melt ! Temperature at melting [K] 
    real(amrex_real), intent(out) :: enth_fus   ! Enthalpy of fusion [kJ/mol]
    real(amrex_real), intent(out) :: rho_melt   ! Density at metling [Kg/m^3]

    temp_melt = 1560.0
    enth_fus = 7.959 + 6.855 ! The latent heat of the hcp-to-bcc transition is included in the enthalpy of fusion
    rho_melt = 1705.0
    
  end subroutine get_melting_point_beryllium


  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy of vaporization
  ! -----------------------------------------------------------------
  subroutine get_enthalpy_of_vaporization_beryllium(enth_vap)

    ! The beryllium enthalpy of vaporization has been adopted from J. W. Arblaster, Thermodynamic Properties of Beryllium, Journal
    ! of Phase Equilibria and Diffusion 37, 581-591 (2016). A very close value is provided by Zhang, Evans and Yang, Corrected values
    ! for boiling points and enthalpies of vaporization of elements in handbooks, J. Chem. Eng. Data 56, 328-337 (2011).
    
    real(amrex_real), intent(out) :: enth_vap ! Enthalpy of vaporization [kj/mol]

    enth_vap = 324.0
    
  end subroutine get_enthalpy_of_vaporization_beryllium  


  ! ! -----------------------------------------------------------------
  ! ! Subroutine used to compute the hcp to bcc point properties 
  ! ! -----------------------------------------------------------------
  ! subroutine get_hcp_to_bcc_point_beryllium(enth_hcp2bcc, hcp2bcc_point, rho_hcp2bcc)

  !   ! The beryllium enthalpy of hcp-to-bcc transformation has been adopted from J. W. Arblaster, Thermodynamic Properties
  !   ! of Beryllium, Journal of Phase Equilibria and Diffusion 37, 581-591 (2016). A rather close value is provided by H. Kleykamp,
  !   ! Selected thermal properties of beryllium and phase equilibria in beryllium systems relevant for nuclear fusion reactor
  !   ! blankets, Journal of Nuclear Materials 294, 88-93 (2001).
    
  !   real(amrex_real), intent(out) :: enth_hcp2bcc     ! Enthalpy of vaporization [kj/mol]
  !   real(amrex_real), intent(out) :: hcp2bcc_point    ! Temperature of phase transition [K]
  !   real(amrex_real), intent(out) :: rho_hcp2bcc      ! Density at solid phase transition point [kg/m^3]

  !   enth_hcp2bcc = 6.855
  !   hcp2bcc_point = 1543.0
  !   rho_hcp2bcc = 1722.1
    
  ! end subroutine get_hcp_to_bcc_point_beryllium  



  ! -----------------------------------------------------------------
  ! Subroutine used to compute the work function
  ! -----------------------------------------------------------------
  subroutine get_work_function_beryllium(Wf)

    ! Early measurements even by well-known specialists referred to oxidized samples, which led to very
    ! low work function values around 3.7eV, see for instance R. G. Wilson, Vacuum Thermionic Work
    ! Functions of Polycrystalline Be, Ti, Cr, Fe, Ni, Cu, Pt and Type 304 Stainless Steel, J. Appl. Phys. 37, 2261
    ! (1966); V. S. Fomenko, The Handbook of thermionic properties, Plenum, New York (1966); R. C. Jernert
    ! and C. B. Magee, Effect of Surface Oxidation on the Thermionic Work Function of Beryllium, Ox. Met.
    ! 2, 1 (1970). This is pointed out in Gmelin, Beryllium, A3, Supplement.
    ! Very good agreement on the work function value in the classical literature, see for instance S. Trasatti,
    ! Electronegativity, Work Function and Heat of Adsorption of Hydrogen on Metals, Chim. Ind. (Milan)
    ! 53, 559 (1971); H. B. Michaelson, The work function of the elements and its periodicity, J. Appl. Phys.
    ! 48, 4729 (1977).
    
    real(amrex_real), intent(out) :: Wf ! Enthalpy of vaporization [J]

    Wf = 4.98
    ! Conversion from eV to Joule
    Wf = Wf*1.60218E-19
    
  end subroutine get_work_function_beryllium    


  ! -----------------------------------------------------------------
  ! Subroutine used to compute Richardson constant
  ! -----------------------------------------------------------------
  subroutine get_richardson_constant_beryllium(Aeff)
    
    ! Lack of measurements of the effective Richardson constant of clean polycrystalline samples.
    ! The recommended value is the nominal Richardson value

    real(amrex_real), intent(out) :: Aeff ! Richardson constant [A/m^2*K^2]

    Aeff = 120E4
    
  end subroutine get_richardson_constant_beryllium


  ! -----------------------------------------------------------------
  ! Subroutine used to compute emissivity
  ! -----------------------------------------------------------------
  subroutine get_emissivity_beryllium(temp, eps_t)

    ! Very few available measurements of the total hemispherical emissivity at elevated temperatures.
    ! Solid beryllium data were adopted from Thermophysical Properties of Materials For Nuclear
    ! Engineering: A Tutorial and Collection of Data (IAEA, Vienna, 2008). The original data stem from G. E.
    ! Darwin and J. H. Buddery, Metallurgy of the rarer metals no.7: Beryllium (Butterworths Scientific
    ! Publications, London, 1960).
    ! These data agree reasonably well with total normal emissivity data provided by Y. S. Touloukian,
    ! Thermophysical properties of matter– The TPRC Data Series - Vol 7 Thermal radiative properties -
    ! Metallic elements and alloys (Plenum, New York, 1970). The same normal emissivity data are quoted
    ! by two beryllium-centered ASM handbooks.
    ! The solid dataset was fitted with a third order Taylor expansion around the room temperature.
    ! In absence of measurements, for liquid beryllium, an absence of discontinuity at the melting point and
    ! a constant value across the liquid phase were assumed.
    
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
  subroutine get_thermelectric_power_beryllium(temp, S)

    ! The discontinuity jump at the melting point has not been
    ! measured, thus it is impossible to justify our assumption that thermoelectric effects do not influence
    ! the replacement current.
    ! Solid beryllium data (600-1100K) were adopted from Lillie’s experiments (1955) that are quoted in
    ! different handbooks, see for instance E. Vidal et al., ASM Handbook: Beryllium and Beryllium Alloys
    ! (ASM International, Ohio, 2010), E. Vidal et al., Beryllium Chemistry & Processing (ASM International,
    ! Ohio, 2009) and Gmelin, Beryllium A3. The data refer to polycrystalline beryllium versus platinum and
    ! have been fitted with a linear expansion around the room temperature.
    ! Absolute solid platinum data (273-1600K) were adopted from R. Roberts, F. Righini & R. Compton,
    ! Absolute scale of thermoelectricity III, Philosophical Magazine Part B, 52, 1147-1163 (1985). The ninth
    ! order Taylor expansion proposed in L. Abadlia et al. "New experimental methodology, setup and
    ! LabView program for accurate absolute thermoelectric power and electrical resistivity measurements
    ! between 25 and 1600 K: Application to pure Cu, Pt, W and Ni at very high temperatures", Review of
    ! Scientific Instruments 85,095121 (2014) has been utilized.
    ! The difference of the fits has been digitized and refitted to a third order Taylor expansion around the
    ! room temperature. The absolute numbers might seem high, but they are similar to other alkanine
    ! earth metals, see e.g. Ba and Ca in the Landolt-Börnstein database.
    ! There are no liquid beryllium measurements and any educated guess is impossible.
    
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
    
  end subroutine get_thermelectric_power_beryllium

  
end module material_properties_beryllium_module
