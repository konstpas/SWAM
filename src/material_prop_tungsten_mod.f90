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
    else if(temp.le.6000.0) then 
       ktherm = 66.6212 + 0.02086*(temp-3695.0) - &
            3.7585E-6*(temp-3695.0)**2
    else 
       ktherm = 95.7345 
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
    else if(temp.le.6000.0) then 
       rho = 16.267 - 7.679E-4*(temp-3695.0) - 8.091E-8*(temp-3695.0)**2 
    else 
       rho = 14.0671 
    end if
    ! Conversion from g/cm3 to kg/cm3 
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
       Cp = 21.868372 + 8.068661E-3*temp - 3.756196E-6*temp**2 + 1.075862E-9*temp**3 + 1.406637E4/temp**2 
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
  ! Subroutine used to computed the atomic mass
  ! -----------------------------------------------------------------
  subroutine get_m_A_tungsten(m_A)

    real(amrex_real), intent(out) :: m_A ! Atomic mass [g/mol]

    m_A = 183.84
    
  end subroutine get_m_A_tungsten

  
  ! -----------------------------------------------------------------
  ! Subroutine used to computed the melting point properties
  ! -----------------------------------------------------------------
  subroutine get_melting_point_tungsten(melt_point, enth_fus, rho_melt)

    real(amrex_real), intent(out) :: melt_point ! Temperature at melting [K] 
    real(amrex_real), intent(out) :: enth_fus ! Enthalpy of fusion [kJ/mol] 
    real(amrex_real), intent(out) :: rho_melt ! Density at melting [kg/m3] 

    melt_point = 3695.0
    enth_fus = 52.3
    rho_melt = 17.1E3
    
  end subroutine get_melting_point_tungsten
  
 
end module material_properties_tungsten_module
