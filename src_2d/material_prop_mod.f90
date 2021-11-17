module material_properties_module  

  ! -----------------------------------------------------------------
  ! This module is used to compute the material properties.
  ! -----------------------------------------------------------------
  
   use amrex_amr_module
   use read_input_module, only : material, &
                                 phiT_table_max_T, &
                                 phiT_table_n_points
   
   use material_properties_tungsten_module, only : get_conductivity_tungsten, &
                                                   get_mass_density_tungsten, &
                                                   get_heat_capacity_tungsten, &
                                                   get_atomic_mass_tungsten, &
                                                   get_melting_point_tungsten, &
                                                   get_electrical_resistivity_tungsten, &
                                                   get_emissivity_tungsten, &
                                                   get_enthalpy_of_vaporization_tungsten, &
                                                   get_Richardson_tungsten, &
                                                   get_surface_tension_tungsten, &
                                                   get_thermelectric_power_tungsten, &
                                                   get_vapor_pressure_tungsten, &
                                                   get_viscosity_tungsten, &
                                                   get_work_function_tungsten
 
 
   use material_properties_test_melting_one_phase_module, only : get_conductivity_test_melting_one_phase, &
                                                                 get_mass_density_test_melting_one_phase, &
                                                                 get_heat_capacity_test_melting_one_phase, &
                                                                 get_atomic_mass_test_melting_one_phase, &
                                                                 get_melting_point_test_melting_one_phase
 
   use material_properties_test_melting_two_phase_module, only : get_conductivity_test_melting_two_phase, &
                                                                 get_mass_density_test_melting_two_phase, &
                                                                 get_heat_capacity_test_melting_two_phase, &
                                                                 get_atomic_mass_test_melting_two_phase, &
                                                                 get_melting_point_test_melting_two_phase    
                                                      
   use material_properties_test_evaporation_module, only :    get_conductivity_test_evaporation, &
                                                              get_mass_density_test_evaporation, &
                                                              get_heat_capacity_test_evaporation, &
                                                              get_atomic_mass_test_evaporation, &
                                                              get_melting_point_test_evaporation, &
                                                              get_electrical_resistivity_test_evaporation, &
                                                              get_emissivity_test_evaporation, &
                                                              get_enthalpy_of_vaporization_test_evaporation, &
                                                              get_Richardson_test_evaporation, &
                                                              get_surface_tension_test_evaporation, &
                                                              get_thermelectric_power_test_evaporation, &
                                                              get_vapor_pressure_test_evaporation, &
                                                              get_viscosity_test_evaporation, &
                                                              get_work_function_test_evaporation                                                      
 
   use material_properties_beryllium_module, only :   get_conductivity_beryllium, &
                                                      get_mass_density_beryllium, &
                                                      get_heat_capacity_beryllium, &
                                                      get_atomic_mass_beryllium, &
                                                      get_melting_point_beryllium, &
                                                      get_electrical_resistivity_beryllium, &
                                                      get_emissivity_beryllium, &
                                                      get_enthalpy_of_vaporization_beryllium, &
                                                      ! get_hcp_to_bcc_point_beryllium, &
                                                      get_Richardson_beryllium, &
                                                      get_surface_tension_beryllium, &
                                                      get_thermelectric_power_beryllium, &
                                                      get_vapor_pressure_beryllium, &
                                                      get_viscosity_beryllium, &
                                                      get_work_function_beryllium

  use material_properties_iridium_module, only :  get_heat_capacity_iridium, &
                                                  get_electrical_resistivity_iridium,&
                                                  get_emissivity_iridium, &
                                                  get_enthalpy_of_vaporization_iridium, &
                                                  get_conductivity_iridium, &
                                                  get_atomic_mass_iridium, &
                                                  get_melting_point_iridium, &
                                                  get_mass_density_iridium, &
                                                  get_Richardson_iridium, &
                                                  get_surface_tension_iridium, &
                                                  get_thermelectric_power_iridium, &
                                                  get_vapor_pressure_iridium, &
                                                  get_viscosity_iridium, &
                                                  get_work_function_iridium

  use material_properties_niobium_module,  only : get_heat_capacity_niobium, &
                                                  get_electrical_resistivity_niobium, &
                                                  get_emissivity_niobium, &
                                                  get_enthalpy_of_vaporization_niobium, &
                                                  get_conductivity_niobium, &
                                                  get_atomic_mass_niobium, &
                                                  get_melting_point_niobium, &
                                                  get_mass_density_niobium, &
                                                  get_Richardson_niobium, &
                                                  get_surface_tension_niobium, &
                                                  get_thermelectric_power_niobium, &
                                                  get_vapor_pressure_niobium, &
                                                  get_viscosity_niobium, &
                                                  get_work_function_niobium
    
     
   implicit none 
 
   private
 
   ! -----------------------------------------------------------------
   ! Public variables
   ! -----------------------------------------------------------------
   ! Enthalpy at the onset of melting
   public :: enth_at_melt
   ! Temperature at the melting point
   public :: temp_melt
   ! Maximum diffusivity
   public :: max_diffus  
   
   ! -----------------------------------------------------------------
   ! Public subroutines
   ! -----------------------------------------------------------------  
   public :: init_mat_prop
   public :: get_conductivity
   public :: get_temp
   public :: get_enthalpy
   public :: get_emissivity
   public :: get_work_function
   public :: get_Richardson
   public :: get_vapor_pressure
   public :: get_atomic_mass
   public :: get_enthalpy_of_vaporization
   public :: get_mass_density
   public :: finalize_mat_prop
   
 
   ! -----------------------------------------------------------------
   ! Declare public variables
   ! -----------------------------------------------------------------
   real(amrex_real), save :: enth_at_melt
   real(amrex_real), save :: temp_melt
   real(amrex_real), save :: max_diffus
   
   ! -----------------------------------------------------------------
   ! Declare private variables shared by all subroutines
   ! -----------------------------------------------------------------
   real(amrex_real), save :: enth_fus
   real(amrex_real), save :: m_A
   real(amrex_real), save :: rho_melt   
   real(amrex_real), allocatable, save :: enth_table(:)
   real(amrex_real), allocatable, save :: temp_table(:)
 
 contains 
   
 
   ! ------------------------------------------------------------------
   ! Subroutine used to compute the thermal conductivity
   ! ------------------------------------------------------------------ 
   subroutine get_conductivity(temp,ktherm)
     
     real(amrex_real), intent(in) :: temp    ! Temperature [K]
     real(amrex_real), intent(out) :: ktherm ! Thermal conductivity [W/mK] 
          
     if (material.eq.'Tungsten') then 
        call get_conductivity_tungsten(temp, ktherm)
     elseif (material.eq.'Beryllium') then
        call get_conductivity_beryllium(temp, ktherm)
     elseif (material.eq.'Iridium') then
         call get_conductivity_iridium(temp, ktherm)
     elseif (material.eq.'Niobium') then
         call get_conductivity_niobium(temp, ktherm)
     elseif(material.eq.'Test_melting_one_phase') then
        call get_conductivity_test_melting_one_phase(temp, ktherm)
     elseif(material.eq.'Test_melting_two_phase') then
        call get_conductivity_test_melting_two_phase(temp, ktherm)
     elseif(material.eq.'Test_Evaporation') then
        call get_conductivity_test_evaporation(temp, ktherm)
     else
        STOP 'Unknown material'
     end if
 
   end subroutine get_conductivity
 
   
   ! ------------------------------------------------------------------
   ! Subroutine used to compute the mass density
   ! ------------------------------------------------------------------ 
   subroutine get_mass_density(temp,rho)
     
     real(amrex_real), intent(in) :: temp   ! Temperature [K]
     real(amrex_real), intent(out) :: rho   ! Mass density [kg/m3] 
        
     if (material.eq.'Tungsten') then 
        call get_mass_density_tungsten(temp,rho)
     elseif (material.eq.'Beryllium') then
        call get_mass_density_beryllium(temp, rho)
     elseif (material.eq.'Iridium') then
        call get_mass_density_iridium(temp, rho)
     elseif (material.eq.'Niobium') then
        call get_mass_density_niobium(temp, rho)
     elseif (material.eq.'Test_melting_one_phase') then
        call get_mass_density_test_melting_one_phase(temp, rho)
     elseif (material.eq.'Test_melting_two_phase') then
        call get_mass_density_test_melting_two_phase(rho)
      elseif (material.eq.'Test_Evaporation') then
         call get_mass_density_test_evaporation(rho)
     else
        STOP 'Unknown material'
     end if
     
   end subroutine get_mass_density
 
   
   ! ------------------------------------------------------------------
   ! Subroutine used to compute the heat capacity
   ! ------------------------------------------------------------------ 
   subroutine get_heat_capacity(temp,Cp)
     
     real(amrex_real), intent(in) :: temp  ! Temperature [K]
     real(amrex_real), intent(out) :: Cp   ! Specific heat capacity [J/kgK] 
        
     if (material.eq.'Tungsten') then 
        call get_heat_capacity_tungsten(temp,Cp)
     elseif (material.eq.'Beryllium') then
        call get_heat_capacity_beryllium(temp, Cp)
     elseif (material.eq.'Iridium') then
        call get_heat_capacity_iridium(temp, Cp)
     elseif (material.eq.'Niobium') then
        call get_heat_capacity_niobium(temp, Cp)
     elseif (material.eq.'Test_melting_one_phase') then
        call get_heat_capacity_test_melting_one_phase(temp, Cp)
     elseif (material.eq.'Test_melting_two_phase') then
        call get_heat_capacity_test_melting_two_phase(Cp)
     elseif (material.eq.'Test_Evaporation') then
        call get_heat_capacity_test_evaporation(temp, Cp)
     else
        STOP 'Unknown material'
     end if
     
   end subroutine get_heat_capacity
         
 
   ! ------------------------------------------------------------------
   ! Subroutine used to compute the atomic mass
   ! ------------------------------------------------------------------ 
   subroutine get_atomic_mass(m_A)

      real(amrex_real), intent(out) :: m_A
 
     if (material.eq.'Tungsten') then 
        call get_atomic_mass_tungsten(m_A)
     elseif (material.eq.'Beryllium') then
        call get_atomic_mass_beryllium(m_A)
     elseif (material.eq.'Iridium') then
        call get_atomic_mass_iridium(m_A)
     elseif (material.eq.'Niobium') then
        call get_atomic_mass_niobium(m_A)
     elseif (material.eq.'Test_melting_one_phase') then
        call get_atomic_mass_test_melting_one_phase(m_A)
     elseif (material.eq.'Test_melting_two_phase') then
        call get_atomic_mass_test_melting_two_phase(m_A)
     elseif (material.eq.'Test_Evaporation') then
        call get_atomic_mass_test_evaporation(m_A)
     else
        STOP 'Unknown material'
     end if
     
   end subroutine get_atomic_mass
   
 
   ! ------------------------------------------------------------------
   ! Subroutine used to compute the properties at the melting point
   ! ------------------------------------------------------------------ 
   subroutine get_melting_point()
 
     if (material.eq.'Tungsten') then 
        call get_melting_point_tungsten(temp_melt, enth_fus, rho_melt)
     elseif (material.eq.'Beryllium') then
         call get_melting_point_beryllium(temp_melt, enth_fus, rho_melt)
     elseif (material.eq.'Iridium') then
         call get_melting_point_iridium(temp_melt, enth_fus, rho_melt)
     elseif (material.eq.'Niobium') then
         call get_melting_point_niobium(temp_melt, enth_fus, rho_melt)
     elseif (material.eq.'Test_melting_one_phase') then
        call get_melting_point_test_melting_one_phase(temp_melt, enth_fus, rho_melt)
     elseif (material.eq.'Test_melting_two_phase') then
        call get_melting_point_test_melting_two_phase(temp_melt, enth_fus, rho_melt)
     elseif (material.eq.'Test_Evaporation') then
        call get_melting_point_test_evaporation(temp_melt, enth_fus, rho_melt)
     else
        STOP 'Unknown material'
     end if
     
   end subroutine get_melting_point
   
   ! ------------------------------------------------------------------
   ! Subroutine used to compute the electrical resistivity
   ! ------------------------------------------------------------------ 
   subroutine get_electrical_resistivity(temp, rho_e)
     
     real(amrex_real), intent(in) :: temp    ! Temperature [K]
     real(amrex_real), intent(out) :: rho_e  ! Electrical resistivity [Ohm*m] 

      if (material.eq.'Tungsten') then
         call get_electrical_resistivity_tungsten(temp, rho_e)
      elseif (material.eq.'Beryllium') then
         call get_electrical_resistivity_beryllium(temp, rho_e)
      elseif (material.eq.'Iridium') then
         call get_electrical_resistivity_iridium(temp, rho_e)
      elseif (material.eq.'Niobium') then
         call get_electrical_resistivity_niobium(temp, rho_e)
      elseif (material.eq.'Test_melting_one_phase' .or. material.eq.'Test_melting_two_phase') then
         rho_e = 0 ! resistivity not defined for test materials, placeholder value
      elseif (material.eq.'Test_Evaporation') then
         call get_electrical_resistivity_test_evaporation(temp, rho_e)
      else
         STOP 'Unknown material'
      endif

    end subroutine get_electrical_resistivity


   ! ------------------------------------------------------------------
   ! Subroutine used to compute the surface tension
   ! ------------------------------------------------------------------ 
   subroutine get_surface_tension(temp, sigma)
     
      real(amrex_real), intent(in) :: temp    ! Temperature [K]
      real(amrex_real), intent(out) :: sigma  ! Surface tension [N/m]

      if (material.eq.'Tungsten') then
         call get_surface_tension_tungsten(temp, sigma)
      elseif (material.eq.'Beryllium') then
         call get_surface_tension_beryllium(temp, sigma)
      elseif (material.eq.'Iridium') then
         call get_surface_tension_iridium(temp, sigma)
      elseif (material.eq.'Niobium') then
         call get_surface_tension_niobium(temp, sigma)
      elseif (material.eq.'Test_melting_one_phase' .or. material.eq.'Test_melting_two_phase') then
         sigma = 0 ! surface tension not defined for test materials, placeholder value
      elseif (material.eq.'Test_Evaporation') then
         call get_surface_tension_test_evaporation(temp, sigma)
      else
         STOP 'Unknown material'
      endif

    end subroutine get_surface_tension


   ! ------------------------------------------------------------------
   ! Subroutine used to compute the viscosity
   ! ------------------------------------------------------------------ 
   subroutine get_viscosity(temp, mu)
     
      real(amrex_real), intent(in) :: temp ! Temperature [K]
      real(amrex_real), intent(out) :: mu  ! Viscosity [Pa*m] 

      if (material.eq.'Tungsten') then
         call get_viscosity_tungsten(temp, mu)
      elseif (material.eq.'Beryllium') then
         call get_viscosity_beryllium(temp, mu)
      elseif (material.eq.'Iridium') then
         call get_viscosity_iridium(temp, mu)
      elseif (material.eq.'Niobium') then
         call get_viscosity_niobium(temp, mu)
      elseif (material.eq.'Test_melting_one_phase' .or. material.eq.'Test_melting_two_phase') then
         mu = 0 ! viscosity not defined for test materials, placeholder value
      elseif (material.eq.'Test_Evaporation') then
         call get_viscosity_test_evaporation(temp, mu)
      else
         STOP 'Unknown material'
      endif

    end subroutine get_viscosity


   ! ------------------------------------------------------------------
   ! Subroutine used to compute the vapor pressure
   ! ------------------------------------------------------------------ 
   subroutine get_vapor_pressure(temp, pv)
     
      real(amrex_real), intent(in) :: temp ! Temperature [K]
      real(amrex_real), intent(out) :: pv  ! Vapor pressure [Pa] 

      if (material.eq.'Tungsten') then
         call get_vapor_pressure_tungsten(temp, pv)
      elseif (material.eq.'Beryllium') then
         call get_vapor_pressure_beryllium(temp, pv)
      elseif (material.eq.'Iridium') then
         call get_vapor_pressure_iridium(temp, pv)
      elseif (material.eq.'Niobium') then
         call get_vapor_pressure_niobium(temp, pv)
      elseif (material.eq.'Test_melting_one_phase' .or. material.eq.'Test_melting_two_phase') then
         pv = 0 ! vapor pressure not defined for test materials, placeholder value
      elseif (material.eq.'Test_Evaporation') then
         call get_vapor_pressure_test_evaporation(temp, pv)
      else
         STOP 'Unknown material'
      endif

    end subroutine get_vapor_pressure


   ! ------------------------------------------------------------------
   ! Subroutine used to compute the enthalpy of vaporization
   ! ------------------------------------------------------------------ 
   subroutine get_enthalpy_of_vaporization(temp, enth_vap)
     
      real(amrex_real), intent(in) :: temp      ! Temperature [K]
      real(amrex_real), intent(out) :: enth_vap ! Enthalpy of vaporization [kj/mol] 

      if (material.eq.'Tungsten') then
         call get_enthalpy_of_vaporization_tungsten(temp, enth_vap)
      elseif (material.eq.'Beryllium') then
         call get_enthalpy_of_vaporization_beryllium(enth_vap)
      elseif (material.eq.'Iridium') then
         call get_enthalpy_of_vaporization_iridium(enth_vap)
      elseif (material.eq.'Niobium') then
         call get_enthalpy_of_vaporization_niobium(enth_vap)
      elseif (material.eq.'Test_melting_one_phase' .or. material.eq.'Test_melting_two_phase') then
         enth_vap = 0 ! enthalpy of vaporization not defined for test materials, placeholder value
      elseif (material.eq.'Test_Evaporation') then
         call get_enthalpy_of_vaporization_test_evaporation(enth_vap)
      else
         STOP 'Unknown material'
      endif
      
   end subroutine get_enthalpy_of_vaporization



   ! ------------------------------------------------------------------
   ! Subroutine used to compute the work function
   ! ------------------------------------------------------------------ 
   subroutine get_work_function(wf)
     
      real(amrex_real), intent(out) :: wf  ! Work function [J] 

      if (material.eq.'Tungsten') then
         call get_work_function_tungsten(wf)
      elseif (material.eq.'Beryllium') then
         call get_work_function_beryllium(wf)
      elseif (material.eq.'Iridium') then
         call get_work_function_iridium(wf)
      elseif (material.eq.'Niobium') then
         call get_work_function_niobium(wf)
      elseif (material.eq.'Test_melting_one_phase' .or. material.eq.'Test_melting_two_phase') then
         wf = 0 ! Work function not defined for test materials, placeholder value
      elseif (material.eq.'Test_Evaporation') then
         call get_work_function_test_evaporation(wf)
      else
         STOP 'Unknown material'
      endif

   end subroutine get_work_function



   ! ------------------------------------------------------------------
   ! Subroutine used to compute the Richardson constant
   ! ------------------------------------------------------------------ 
   subroutine get_Richardson(Aeff)
     
      real(amrex_real), intent(out) :: Aeff  ! Richardson constant [A/(m^2K^2)] 

      if (material.eq.'Tungsten') then
         call get_Richardson_tungsten(Aeff)
      elseif (material.eq.'Beryllium') then
         call get_Richardson_beryllium(Aeff)
      elseif (material.eq.'Iridium') then
         call get_Richardson_iridium(Aeff)
      elseif (material.eq.'Niobium') then
         call get_Richardson_niobium(Aeff)
      elseif (material.eq.'Test_melting_one_phase' .or. material.eq.'Test_melting_two_phase') then
         Aeff = 0 ! Richardson constant not defined for test materials, placeholder value
      elseif (material.eq.'Test_Evaporation') then
         call get_Richardson_test_evaporation(Aeff)
      else
         STOP 'Unknown material'
      endif
      
   end subroutine get_Richardson


   ! ------------------------------------------------------------------
   ! Subroutine used to compute the emissivity
   ! ------------------------------------------------------------------ 
   subroutine get_emissivity(temp, eps_t)
     
      real(amrex_real), intent(in) :: temp    ! Temperature [K]
      real(amrex_real), intent(out) :: eps_t  ! Emissivity [dimensionless] 

      if (material.eq.'Tungsten') then
         call get_emissivity_tungsten(temp, eps_t)
      elseif (material.eq.'Beryllium') then
         call get_emissivity_beryllium(temp, eps_t)
      elseif (material.eq.'Iridium') then
         call get_emissivity_iridium(temp, eps_t)
      elseif (material.eq.'Niobium') then
         call get_emissivity_niobium(temp, eps_t)
      elseif (material.eq.'Test_melting_one_phase' .or. material.eq.'Test_melting_two_phase') then
         eps_t = 0 ! Richardson constant not defined for test materials, placeholder value
      elseif (material.eq.'Test_Evaporation') then
         call get_emissivity_test_evaporation(eps_t)
      else
         STOP 'Unknown material'
      endif
      
   end subroutine get_emissivity
   
   

   ! ------------------------------------------------------------------
   ! Subroutine used to compute the absolute thermoelectric power
   ! ------------------------------------------------------------------ 
   subroutine get_thermelectric_power(temp, S)
     
      real(amrex_real), intent(in) :: temp   ! Temperature [K]
      real(amrex_real), intent(out) :: S     ! Absolute thermoelectric power [V/K] 

      if (material.eq.'Tungsten') then
         call get_thermelectric_power_tungsten(temp, S)
      elseif (material.eq.'Beryllium') then
         call get_thermelectric_power_beryllium(temp, S)
      elseif (material.eq.'Iridium') then
         call get_thermelectric_power_iridium(temp, S)
      elseif (material.eq.'Niobium') then
         call get_thermelectric_power_niobium(temp, S)
      elseif (material.eq.'Test_melting_one_phase' .or. material.eq.'Test_melting_two_phase') then
         S = 0 ! Thermoelectric power not defined for test materials, placeholder value
      elseif (material.eq.'Test_Evaporation') then
         call get_thermelectric_power_test_evaporation(temp, S)
      else
         STOP 'Unknown material'
      endif
      
   end subroutine get_thermelectric_power
   

  
  ! ------------------------------------------------------------------
  ! Subroutine used to compute the tables that relate the enthalpy
  ! and the temperature. Note that the temperature is in K and the
  ! enthalpy in J/m^3
  ! ------------------------------------------------------------------ 
  subroutine init_mat_prop()

    integer :: i
    integer :: imelt = 0
    logical :: isolid = .true.  ! for enthalpy table, true before phase transfer
    real(amrex_real) :: Cp
    real(amrex_real) :: diffus
    real(amrex_real) :: ktherm 
    real(amrex_real) :: phiT_table_dT  
    real(amrex_real) :: rho
    real(amrex_real) :: rhocp_i
    real(amrex_real) :: rhocp_im1
    real(amrex_real) :: rho_e
    real(amrex_real) :: mu
    real(amrex_real) :: sigma
    real(amrex_real) :: pv
    real(amrex_real) :: Wf
    real(amrex_real) :: Aeff
    real(amrex_real) :: eps_t
    real(amrex_real) :: S
    real(amrex_real) :: enth_vap

    ! Allocate the temperature and enthalpy tables
    allocate(temp_table(0:phiT_table_n_points))
    allocate(enth_table(0:phiT_table_n_points))     
    
    ! Initialize maximum diffusivity
    max_diffus = 0.
    
    ! Atomic mass and thermodynamic properties at the melting point
    call get_atomic_mass(m_A)
    call get_melting_point
    
    ! The enthalpy is phi = int_0^T rho(T')*Cp(T') dT' and trapezoidal integration is used
    ! to compute the integral. The calculation of the enthalpy includes the phase transfer
    ! discountinuity. 

    ! Table increment
    phiT_table_dT = phiT_table_max_T/phiT_table_n_points

    ! Properties at zero temperature
    temp_table(0) = 0_amrex_real
    call get_mass_density(temp_table(0),rho) 
    call get_heat_capacity(temp_table(0),Cp)  
    rhocp_i = rho*Cp ! product of density and heat capacity 
    enth_table(0) = 0_amrex_real

    ! Fill the enthalpy-temperature table
    do i = 1,phiT_table_n_points                 

       ! Update temperature
       temp_table(i) = temp_table(i-1) + phiT_table_dT 
       
       ! Check if the updated temperature falls above the melting point
       if (temp_table(i).ge.temp_melt) then

          ! Operations to perform for the first data point above melting temperature  
          if (isolid) then  

             isolid = .false.   
             imelt = i          
             temp_table(i) = temp_melt 
             phiT_table_dT = temp_melt - temp_table(i-1)
             rhocp_im1 = rhocp_i 
             call get_mass_density(temp_table(i),rho) 
             call get_heat_capacity(temp_table(i),Cp)  
             rhocp_i = rho*Cp  
             enth_table(i) = enth_table(i-1) + (rhocp_i+rhocp_im1)*phiT_table_dT/2_amrex_real   ! Enthalpy at melt onset 
             enth_at_melt = enth_table(i) 
             phiT_table_dT = (phiT_table_max_T - temp_melt)/(phiT_table_n_points-1-i)  ! New phiT_table_dT to match phiT_table_max_T !
             
          end if

          ! Solid-liquid phase transfer jump
          if(imelt.eq.i-1) then

             temp_table(i) = temp_melt
             enth_table(i) = enth_table(i-1) + 1E6*enth_fus*rho_melt/m_A ! [J/m3] = 1E6*[kJ/mol]*[kg/m3]/[g/mol]

          ! Compute enthalpy for all the state points above the melting point
          elseif(imelt.ne.i) then
             
             rhocp_im1 = rhocp_i 
             call get_mass_density(temp_table(i),rho) 
             call get_heat_capacity(temp_table(i),Cp)  
             rhocp_i = rho*Cp       
             enth_table(i) = enth_table(i-1) + (rhocp_i+rhocp_im1)*phiT_table_dT/2_amrex_real  ! Trapezoidal integration
             
          end if
          
       else

          ! Compute enthalpy for all the state points below the melting point
          rhocp_im1 = rhocp_i 
          call get_mass_density(temp_table(i),rho) 
          call get_heat_capacity(temp_table(i),Cp)  
          rhocp_i = rho*Cp 
          enth_table(i) = enth_table(i-1) + (rhocp_i+rhocp_im1)*phiT_table_dT/2_amrex_real  ! Trapezoidal integration
          
       end if

       ! Initialize the maximum diffusivity  
       call get_conductivity(temp_table(i),ktherm)
       diffus = ktherm/rhocp_i  
       if (diffus.gt.max_diffus) max_diffus=diffus
      
    end do
    
    ! Output employed material properties to file
    open (2, file = 'material_properties_'//TRIM(material)//'.dat', status = 'unknown')
    write(2, *) 'Material properties employed' 
    write(2, *) 'Temperature[K], Cp [J/kgK], rho [kg/m^3], k [W/mk], enthalpy [J/m^3]', &
               'Resistivity [Ohm m], Surface tension [N/m], Viscosity [sPa], Vapour pressure [Pa]', &
               'Work function [J], Richardson constant [A/(m^2K^2)], Emissivity [-], Abs. Thermoelectric power [V/K]', &
               'Enthalpy of vaporization [kJ/mol]' 
    do i = 0,phiT_table_n_points
       call get_heat_capacity(temp_table(i),Cp) 
       call get_conductivity(temp_table(i),ktherm) 
       call get_mass_density(temp_table(i),rho)
       call get_electrical_resistivity(temp_table(i), rho_e)
       call get_surface_tension(temp_table(i), sigma)
       call get_viscosity(temp_table(i), mu)
       call get_vapor_pressure(temp_table(i), pv)
       call get_work_function(Wf)
       call get_Richardson(Aeff)
       call get_emissivity(temp_table(i),  eps_t)
       call get_thermelectric_power(temp_table(i), S)
       call get_enthalpy_of_vaporization(temp_table(i), enth_vap)
       write(2,*) temp_table(i), Cp, rho, ktherm, enth_table(i), rho_e, sigma, mu, pv, Wf, Aeff, S, enth_vap
    end do
    close(2) 


  end subroutine init_mat_prop

  ! ------------------------------------------------------------------
  ! Subroutine used to deallocate the tables that relate the enthalpy
  ! and the temperature
  ! ------------------------------------------------------------------ 
  subroutine finalize_mat_prop()

    deallocate(temp_table)
    deallocate(enth_table)
    
  end subroutine finalize_mat_prop
  
 
  ! ------------------------------------------------------------------
  ! Subroutine used to obtain the temperature given the enthalpy
  ! ------------------------------------------------------------------ 
  subroutine get_temp(lo, hi, &
                      ui, uo_lo, uo_hi, & 
                      temp, t_lo , t_hi) 

    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: uo_lo(2), uo_hi(2)
    integer, intent(in) :: t_lo(2), t_hi(2)
    real(amrex_real), intent(in) :: ui (uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2))
    real(amrex_real), intent(out) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2))

    ! Local variables
    integer :: idx 
    integer :: i,j
    real(amrex_real) :: int_coeff
    ! real(amrex_real) :: Cp
    ! real(amrex_real) :: diffus 
    ! real(amrex_real) :: ktherm
    ! real(amrex_real) :: rho
     

    ! Obtain the temperature from linear interpolation of the enthalpy-temperature tables
    do i = lo(1),hi(1)
       do j = lo(2),hi(2)

          if(ui(i,j).ne.ui(i,j)) then
             STOP 'Enthalpy query in temperature table is nan.'
          end if

          if(ui(i,j).eq.ui(i,j)-1) then
             STOP 'Enthalpy query for temperature table is infinity.'
          end if

          if(ui(i,j).lt.0) then
            write(*,*) ui(i,j)
            STOP 'Enthalpy query for temperature table is negative.'
          end if

          ! Routine bisection assumes that first argument is an array 
          ! with its index starting from one. Since enth_table
          ! starts from 0 a +1 is added in the second argument.
          call bisection(enth_table, phiT_table_n_points+1, ui(i,j), idx)
          
          if (idx.ge.phiT_table_n_points) then
             STOP 'Temperature table exceeded' 
          end if
         
          int_coeff = (ui(i,j)-enth_table(idx-1))/ &
               (enth_table(idx)-enth_table(idx-1))
          temp(i,j) = temp_table(idx-1) + &
                      int_coeff*(temp_table(idx)-temp_table(idx-1))
          
          ! Update maximum diffusivity (consider only material grid points and not the background)
          ! if (temp(i,j).gt.0) then
          !    call get_conductivity(temp(i,j),ktherm)
          !    call get_mass_density(temp(i,j),rho) 
          !    call get_heat_capacity(temp(i,j),Cp)
          !    diffus = ktherm/(rho*Cp)
          !    if (diffus.gt.max_diffus) then
          !       max_diffus = diffus
          !    end if
          ! end if
          
       end do
    end do
    
  end subroutine get_temp


  ! ------------------------------------------------------------------
  ! Subroutine used to obtain the enthalpy given the temperature. It
  ! is only used during the initialization phase when the temperature
  ! passed in input should be translated into an enthalpy or for
  ! simulations with fixed temperature on the free surface in order
  ! to prescribe the temperature on the free surface
  ! ------------------------------------------------------------------ 
  subroutine get_enthalpy(temp,enth) 

    use read_input_module, only : phase_init
    
    ! Input and output variables
    real(amrex_real), intent(in) :: temp
    real(amrex_real), intent(out) :: enth

    ! Local variables
    integer :: idx
    real(amrex_real) :: int_coeff

    if(temp.ne.temp) then
      STOP 'Temperature query for enthalpy table is nan.'
    end if

    if(temp.eq.temp-1) then
      STOP 'Temperature query for enthalpy table is infinity.'
    end if

    do idx = 0,phiT_table_n_points 
       if (temp .le. temp_table(idx)) exit 
    end do
    
    if (idx.eq.phiT_table_n_points) then 
      STOP 'Temperature table exceeded'
    end if
    
    ! If the input temperature is the melting temperature the enthalpy is ambiguous and
    ! it should be specified if the system is to be considered solid or liquid
    if (temp .eq. temp_melt) then
       
       if (phase_init .eq. "solid") then
          enth = enth_table(idx)
       else if (phase_init .eq. "liquid") then
          enth = enth_table(idx+1)
       else
          STOP "get_enthalpy: For systems at the melting temperature the phase (liquid or solid) should be specified"
       end if
       
    else ! In all other cases, interpolate from table

       int_coeff = (temp-temp_table(idx-1))/(temp_table(idx)-temp_table(idx-1))
       enth = enth_table(idx-1) + int_coeff*(enth_table(idx)-enth_table(idx-1))
    
    end if
    
  end subroutine get_enthalpy
  
   ! -----------------------------------------------------------------
   ! Given an array xx(1:n), and given a value x, returns a value j such that x is between
   ! xx(j) and xx(j+1). xx(1:n) must be monotonic, either increasing or decreasing. j=0
   ! or j=n is returned to indicate that x is out of range.
   ! Taken from Numerical recipes in Fortran 77.
   ! -----------------------------------------------------------------
  subroutine bisection(xx,n,x,j)

   ! Input and output variables
   integer, intent(in) :: n
   real(amrex_real), intent(in) :: xx(n)
   real(amrex_real), intent(in) :: x
   integer, intent(out) :: j

   ! Local variables

   integer jl,jm,ju

   jl=0 ! Initialize lower
   ju=n+1 ! upper limits.
   do while(ju-jl.gt.1)
      jm=(ju+jl)/2
      if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
         jl=jm
      else
         ju=jm
      endif
   end do
   if(x.eq.xx(1)) then
      j=1
   else if(x.eq.xx(n)) then
      j=n-1
   else
      j=jl
   endif
end subroutine bisection

end module material_properties_module 
