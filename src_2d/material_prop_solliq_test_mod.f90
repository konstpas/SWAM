module material_properties_solliq_test_module  

  use amrex_amr_module
  
  implicit none 

  ! -----------------------------------------------------------------
  ! This module is used to compute the material properties of
  ! a test material 
  ! -----------------------------------------------------------------
  
  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: get_ktherm_solliq_test
  public :: get_m_A_solliq_test
  public :: get_melting_point_solliq_test
  public :: get_rho_solliq_test
  public :: get_Cp_solliq_test

contains 

  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the thermal conductivity
  ! -----------------------------------------------------------------
  subroutine get_ktherm_solliq_test(temp, ktherm)

    real(amrex_real), intent(in) :: temp ! Temperature [K]
    real(amrex_real), intent(out) :: ktherm ! Thermal conductivity [W/mK] 

    if (temp.lt.3695) then
        ktherm = 150.0
    elseif (temp.ge. 3695) then
        ktherm = 80.0
    endif
     
    
  end subroutine get_ktherm_solliq_test


  ! -----------------------------------------------------------------
  ! Subroutine used to computed the mass density
  ! -----------------------------------------------------------------
  subroutine get_rho_solliq_test(rho)

    ! Input and output variables 
    real(amrex_real), intent(out) :: rho    ! Mass density [kg/m3]
 
    rho = 20.00 
    ! Conversion from g/cm3 to kg/cm3 
    rho = rho*1E3  
   
  end subroutine get_rho_solliq_test
  

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the heat capacity
  ! -----------------------------------------------------------------
  subroutine get_Cp_solliq_test(Cp) 

    ! Input and output variables 
    real(amrex_real), intent(out) :: Cp     ! Specific heat capacity [J/kgK]

    ! Local variables
    real(amrex_real) :: m_A !Atomic mass [g/mol]
    
    Cp = 20.0
    ! Conversion from J/(mol*K) to J/(kg*K) 
    call get_m_A_solliq_test(m_A)
    Cp = 1E3*(Cp/m_A)
    
  end subroutine get_Cp_solliq_test
  

  ! -----------------------------------------------------------------
  ! Subroutine used to computed the atomic mass
  ! -----------------------------------------------------------------
  subroutine get_m_A_solliq_test(m_A)

    real(amrex_real), intent(out) :: m_A ! Atomic mass [g/mol]

    m_A = 183.84
    
  end subroutine get_m_A_solliq_test

  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the melting point properties
  ! -----------------------------------------------------------------
  subroutine get_melting_point_solliq_test(temp_melt, enth_fus, rho_melt)

    real(amrex_real), intent(out) :: temp_melt ! Temperature at melting [K] 
    real(amrex_real), intent(out) :: enth_fus ! Enthalpy of fusion [kJ/mol] 
    real(amrex_real), intent(out) :: rho_melt ! Density at melting [kg/m3] 

    temp_melt = 3695.0
    enth_fus = 52.3
    rho_melt = 20.0E3
    
  end subroutine get_melting_point_solliq_test
  
 
end module material_properties_solliq_test_module
