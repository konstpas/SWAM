module material_properties_module  

  ! -----------------------------------------------------------------
  ! This module is used to compute the material properties.
  ! NOTE: As of July 24, 2021 the maximum diffusivity is set in
  ! init_mat_prop and never changed again. This should probably
  ! be changed
  ! -----------------------------------------------------------------
  
  use amrex_amr_module
  use read_input_module, only : material, &
                                phiT_table_max_T, &
                                phiT_table_n_points
  use material_properties_tungsten_module, only : get_ktherm_tungsten, &
                                                  get_rho_tungsten, &
                                                  get_Cp_tungsten, &
                                                  get_m_A_tungsten, &
                                                  get_melting_point_tungsten

  
  implicit none 

  private

  ! -----------------------------------------------------------------
  ! Public variables
  ! -----------------------------------------------------------------
  ! Enthalpy at the onset of melting
  public :: enth_at_melt
  ! Temperature at the melting point
  public :: melt_point
  ! Maximum diffusivity
  public :: max_diffus  
  
  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------  
  public :: init_mat_prop
  public :: get_ktherm
  public :: get_temp
  public :: get_enthalpy

  ! -----------------------------------------------------------------
  ! Declare public variables
  ! -----------------------------------------------------------------
  real(amrex_real), save :: enth_at_melt
  real(amrex_real), save :: melt_point
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
  subroutine get_ktherm(temp,ktherm)
    
    real(amrex_real), intent(in) :: temp    ! Temperature [K]
    real(amrex_real), intent(out) :: ktherm ! Thermal conductivity [W/mK] 
    	 
    if (material.eq.'Tungsten') then 
       call get_ktherm_tungsten(temp,ktherm)
    else
       STOP 'Unknown material'
    end if

  end subroutine get_ktherm

  
  ! ------------------------------------------------------------------
  ! Subroutine used to compute the mass density
  ! ------------------------------------------------------------------ 
  subroutine get_rho(temp,rho)
    
    real(amrex_real), intent(in) :: temp   ! Temperature [K]
    real(amrex_real), intent(out) :: rho   ! Mass density [kg/m3] 
  	 
    if (material.eq.'Tungsten') then 
       call get_rho_tungsten(temp,rho)
    else
       STOP 'Unknown material'
    end if
    
  end subroutine get_rho

  
  ! ------------------------------------------------------------------
  ! Subroutine used to compute the heat capacity
  ! ------------------------------------------------------------------ 
  subroutine get_Cp(temp,Cp)
    
    real(amrex_real), intent(in) :: temp  ! Temperature [K]
    real(amrex_real), intent(out) :: Cp   ! Specific heat capacity [J/kgK] 
  	 
    if (material.eq.'Tungsten') then 
       call get_Cp_tungsten(temp,Cp)
    else
       STOP 'Unknown material'
    end if
    
  end subroutine get_Cp
		

  ! ------------------------------------------------------------------
  ! Subroutine used to compute the atomic mass
  ! ------------------------------------------------------------------ 
  subroutine get_m_A()

    if (material.eq.'Tungsten') then 
       call get_m_A_tungsten(m_A)
    else
       STOP 'Unknown material'
    end if
    
  end subroutine get_m_A
  

  ! ------------------------------------------------------------------
  ! Subroutine used to compute the properties at the melting point
  ! ------------------------------------------------------------------ 
  subroutine get_melting_point()

    if (material.eq.'Tungsten') then 
       call get_melting_point_tungsten(melt_point, enth_fus, rho_melt)
    else
       STOP 'Unknown material'
    end if
    
  end subroutine get_melting_point

  
  
  ! ------------------------------------------------------------------
  ! Subroutine used to compute the tables that relate the enthalpy
  ! and the temperature. Note that the temperature is in K and the
  ! enthalpy in J/m^3
  ! ------------------------------------------------------------------ 
  subroutine init_mat_prop()

    use read_input_module, only : tempinit, temp_fs
    
    integer :: i
    integer :: imelt = 0
    logical :: isolid = .true.  ! for enthalpy table, true before phase transfer
    real(amrex_real) :: Cp
    real(amrex_real) :: ktherm 
    real(amrex_real) :: phiT_table_dT  
    real(amrex_real) :: rho
    real(amrex_real) :: rhocp_i
    real(amrex_real) :: rhocp_im1

    ! Allocate the temperature and enthalpy tables
    allocate(temp_table(0:phiT_table_n_points))
    allocate(enth_table(0:phiT_table_n_points))     
    
    ! Initialize maximum diffusivity
    max_diffus = 0.
    
    ! Atomic mass and thermodynamic properties at the melting point
    call get_m_A
    call get_melting_point
    
    ! The enthalpy is phi = int_0^T rho(T')*Cp(T') dT' and trapezoidal integration is used
    ! to compute the integral. The calculation of the enthalpy includes the phase transfer
    ! discountinuity. 

    ! Table increment
    phiT_table_dT = phiT_table_max_T/phiT_table_n_points

    ! Properties at zero temperature
    temp_table(0) = 0_amrex_real
    call get_rho(temp_table(0),rho) 
    call get_Cp(temp_table(0),Cp)  
    rhocp_i = rho*Cp ! product of density and heat capacity 
    enth_table(0) = 0_amrex_real

    ! Fill the enthalpy-temperature table
    do i = 1,phiT_table_n_points                 

       ! Update temperature
       temp_table(i) = temp_table(i-1) + phiT_table_dT 
       
       ! Check if the updated temperature falls above the melting point
       if (temp_table(i).ge.melt_point) then

          ! Operations to perform for the first data point above melting temperature  
          if (isolid) then  

             isolid = .false.   
             imelt = i          
             temp_table(i) = melt_point 
             phiT_table_dT = melt_point - temp_table(i-1)
             rhocp_im1 = rhocp_i 
             call get_rho(temp_table(i),rho) 
             call get_Cp(temp_table(i),Cp)  
             rhocp_i = rho*Cp  
             enth_table(i) = enth_table(i-1) + (rhocp_i+rhocp_im1)*phiT_table_dT/2_amrex_real   ! Enthalpy at melt onset 
             enth_at_melt = enth_table(i) 
             phiT_table_dT = (phiT_table_max_T - melt_point)/(phiT_table_n_points-1-i)  ! New phiT_table_dT to match phiT_table_max_T !
             
          end if

          ! Solid-liquid phase transfer jump
          if(imelt.eq.i-1) then

             temp_table(i) = melt_point
             enth_table(i) = enth_table(i-1) + 1E6*enth_fus*rho_melt/m_A ! [J/m3] = 1E6*[kJ/mol]*[kg/m3]/[g/mol]

          ! Compute enthalpy for all the state points above the melting point
          elseif(imelt.ne.i) then
             
             rhocp_im1 = rhocp_i 
             call get_rho(temp_table(i),rho) 
             call get_Cp(temp_table(i),Cp)  
             rhocp_i = rho*Cp       
             enth_table(i) = enth_table(i-1) + (rhocp_i+rhocp_im1)*phiT_table_dT/2_amrex_real  ! Trapezoidal integration
             
          end if
          
       else

          ! Compute enthalpy for all the state points below the melting point
          rhocp_im1 = rhocp_i 
          call get_rho(temp_table(i),rho) 
          call get_Cp(temp_table(i),Cp)  
          rhocp_i = rho*Cp 
          enth_table(i) = enth_table(i-1) + (rhocp_i+rhocp_im1)*phiT_table_dT/2_amrex_real  ! Trapezoidal integration
          
       end if
   
    end do

    ! Compute diffusivity corresponding to
    ! (a) the input temperature for calculations with imposed flux at the free surface
    ! (b) the free surface for calculations with imposed temperature at the free surface
    if (temp_fs.gt.tempinit) then
       call get_ktherm(temp_fs,ktherm)
       call get_rho(temp_fs,rho) 
       call get_Cp(temp_fs,Cp) 
    else
       call get_ktherm(tempinit,ktherm)
       call get_rho(tempinit,rho) 
       call get_Cp(tempinit,Cp)
    end if
    max_diffus = ktherm/(rho*Cp)
    
    ! Output employed material properties to file
    open (2, file = 'material_properties_'//TRIM(material)//'.dat', status = 'unknown')
    write(2,*) 'Material properties employed' 
    write(2,*) 'Temperature[K], Cp [J/kgK], rho [kg/m^3], k [W/mk], enthalpy [J/m^3]' 
    do i = 0,phiT_table_n_points
       call get_Cp(temp_table(i),Cp) 
       call get_ktherm(temp_table(i),ktherm) 
       call get_rho(temp_table(i),rho)
       write(2,*) temp_table(i), Cp, rho, ktherm, enth_table(i)
    end do
    close(2) 

  end subroutine init_mat_prop
  
 
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
    integer :: e_ind 
    integer :: i,j
    real(amrex_real) :: Cp
    real(amrex_real) :: diffus
    real(amrex_real) :: int_coeff 
    real(amrex_real) :: ktherm
    real(amrex_real) :: rho
     
    ! Obtain the temperature from linear interpolation of the enthalpy-temperature tables
    do i = lo(1),hi(1)
       do j = lo(2),hi(2)
          
          do e_ind = 0,phiT_table_n_points 
             if (ui(i,j) .le. enth_table(e_ind) ) exit 
          end do
          
          if (e_ind.eq.phiT_table_n_points) STOP 'Temperature table exceeded' 
          
          int_coeff = (ui(i,j)-enth_table(e_ind-1))/ &
               (enth_table(e_ind)-enth_table(e_ind-1))
          temp(i,j) = temp_table(e_ind-1) + &
               int_coeff*(temp_table(e_ind)-temp_table(e_ind-1))
          
          ! Update maximum diffusivity (consider only material grid points and not the background)
          if (temp(i,j).gt.0) then
             call get_ktherm(temp(i,j),ktherm)
             call get_rho(temp(i,j),rho) 
             call get_Cp(temp(i,j),Cp)
             diffus = ktherm/(rho*Cp)
             if (diffus.gt.max_diffus) then
                max_diffus = diffus
             end if
          end if
          
       end do
    end do
    
  end subroutine get_temp


  ! ------------------------------------------------------------------
  ! Subroutine used to obtain the enthalpy given the temperature. It
  ! is only used during the initialization phase when the temperature
  ! passed in input should be translated into an enthalpy
  ! ------------------------------------------------------------------ 
  subroutine get_enthalpy(temp,enth) 

    integer :: e_ind 
    real(amrex_real), intent(in) :: temp
    real(amrex_real), intent(out) :: enth
    real(amrex_real) :: int_coeff

    ! Obtain the enthalpy from linear interpolation of the enthalpy-temperature tables
    do e_ind = 0,phiT_table_n_points 
       if (temp .le. temp_table(e_ind) ) exit 
    end do
    
    if (e_ind.eq.phiT_table_n_points) STOP 'Temperature table exceeded'
    
    int_coeff = (temp-temp_table(e_ind-1))/(temp_table(e_ind)-temp_table(e_ind-1))
    enth = enth_table(e_ind-1) + int_coeff*(enth_table(e_ind)-enth_table(e_ind-1))
     
  end subroutine get_enthalpy


end module material_properties_module 
