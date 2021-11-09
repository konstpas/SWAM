module material_properties_test_melting_one_phase_module  

  use amrex_amr_module
  
  implicit none 

  ! -----------------------------------------------------------------
  ! This module is used to compute the material properties of
  ! the material used to test the melting module via the
  ! one region stephan problem
  ! -----------------------------------------------------------------
  
  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: get_heat_capacity_test_melting_one_phase
  public :: get_conductivity_test_melting_one_phase
  public :: get_atomic_mass_test_melting_one_phase
  public :: get_melting_point_test_melting_one_phase 
  public :: get_mass_density_test_melting_one_phase

contains 

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the thermal conductivity
  ! -----------------------------------------------------------------
  subroutine get_conductivity_test_melting_one_phase(temp,ktherm)

    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: ktherm ! Thermal conductivity [W/mK] 
    
    if (temp .ge. 0_amrex_real) then
       ktherm = 150.0
    else
       STOP "Negative temperature!"
    end if
    
  end subroutine get_conductivity_test_melting_one_phase
  
  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the mass density
  ! -----------------------------------------------------------------
  subroutine get_mass_density_test_melting_one_phase(temp,rho)
    
    ! Input and output variables 
    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: rho   ! Mass density [kg/m3]
    
    if (temp .ge. 0_amrex_real) then
       rho = 20.0E3
    else
       STOP "Negative temperature!"
    end if
    
  end subroutine get_mass_density_test_melting_one_phase
  
  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the heat capacity
  ! -----------------------------------------------------------------
  subroutine get_heat_capacity_test_melting_one_phase(temp,Cp) 
    
    ! Input and output variables 
    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: Cp    ! Specific heat capacity [J/kgK]
    
    ! Local variables
    real(amrex_real) :: m_A ! Atomic mass [g/mol]
    
    if (temp .ge. 0_amrex_real) then
       Cp = 20.0
    else
       STOP "Negative temperature!"
    end if
    
    ! Conversion from J/(mol*K) to J/(kg*K) 
    call get_atomic_mass_test_melting_one_phase(m_A)
    Cp = 1E3*(Cp/m_A)
    
  end subroutine get_heat_capacity_test_melting_one_phase
  
  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the atomic mass
  ! -----------------------------------------------------------------
  subroutine get_atomic_mass_test_melting_one_phase(m_A)
    
    real(amrex_real), intent(out) :: m_A ! Atomic mass [g/mol]
    
    m_A = 183.84
    
  end subroutine get_atomic_mass_test_melting_one_phase
  
  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the melting point properties
  ! -----------------------------------------------------------------
  subroutine get_melting_point_test_melting_one_phase(temp_melt, enth_fus, rho_melt)
    
    real(amrex_real), intent(out) :: temp_melt ! Temperature at melting [K] 
    real(amrex_real), intent(out) :: enth_fus ! Enthalpy of fusion [kJ/mol] 
    real(amrex_real), intent(out) :: rho_melt ! Density at melting [kg/m3] 
    
    temp_melt = 3695.0
    enth_fus = 52.3
    rho_melt = 20.0E3
    
  end subroutine get_melting_point_test_melting_one_phase
  
  
 
end module material_properties_test_melting_one_phase_module
