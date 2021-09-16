module material_properties_test_module

  use amrex_amr_module
  
  implicit none 

  ! -----------------------------------------------------------------
  ! This module is used to compute the material properties of a
  ! test material used for benchmarking
  ! -----------------------------------------------------------------
  
  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: get_Cp_test
  public :: get_ktherm_test
  public :: get_m_A_test
  public :: get_melting_point_test 
  public :: get_rho_test

contains 

  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the thermal conductivity
  ! -----------------------------------------------------------------
  subroutine get_ktherm_test(temp,ktherm)

    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: ktherm ! Thermal conductivity [W/mK] 
    
    ktherm = 149.441
    
  end subroutine get_ktherm_test


  ! -----------------------------------------------------------------
  ! Subroutine used to computed the mass density
  ! -----------------------------------------------------------------
  subroutine get_rho_test(temp,rho)

    ! Input and output variables 
    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: rho    ! Mass density [kg/m3]

    rho = 19.25*1E3
   
  end subroutine get_rho_test
  

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the heat capacity
  ! -----------------------------------------------------------------
  subroutine get_Cp_test(temp,Cp) 

    ! Input and output variables 
    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: Cp     ! Specific heat capacity [J/kgK]

    ! Local variables
    real(amrex_real) :: m_A ! Atomic mass [g/mol]
    
    Cp = 21.868372
    ! Conversion from J/(mol*K) to J/(kg*K) 
    call get_m_A_test(m_A)
    Cp = 1E3*(Cp/m_A)
    
  end subroutine get_Cp_test
  

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the atomic mass
  ! -----------------------------------------------------------------
  subroutine get_m_A_test(m_A)

    real(amrex_real), intent(out) :: m_A ! Atomic mass [g/mol]

    m_A = 183.84
    
  end subroutine get_m_A_test

  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the melting point properties
  ! -----------------------------------------------------------------
  subroutine get_melting_point_test(melt_point, enth_fus, rho_melt)

    real(amrex_real), intent(out) :: melt_point ! Temperature at melting [K] 
    real(amrex_real), intent(out) :: enth_fus ! Enthalpy of fusion [kJ/mol] 
    real(amrex_real), intent(out) :: rho_melt ! Density at melting [kg/m3] 

    melt_point = 3695.0
    enth_fus = 52.3
    rho_melt = 19.25E3
    
  end subroutine get_melting_point_test
  
 
end module material_properties_test_module
