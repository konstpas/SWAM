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
  public :: get_maxdiffus

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
  ! and the temperature
  ! ------------------------------------------------------------------ 
  subroutine init_mat_prop()

    integer :: i
    integer :: imelt
    logical :: isolid = .true.  ! for enthalpy table, true before phase transfer
    real(amrex_real) :: Cp
    real(amrex_real) :: diffus
    real(amrex_real) :: ktherm 
    real(amrex_real) :: phiT_table_dT 	! end of temperature interval (decided in input), and temperature increment size  
    real(amrex_real) :: rho
    real(amrex_real) :: rhocp_i
    real(amrex_real) :: rhocp_im1

    ! Initialize maximum diffusivity
    max_diffus = 0.
    
    ! Atomic mass and thermodynamic properties at the melting point
    call get_m_A
    call get_melting_point
    
    ! give end temperature point and data points in input file? 
    ! now called several times, move call to initdata later 
    allocate(temp_table(0:phiT_table_n_points))
    allocate(enth_table(0:phiT_table_n_points))     
    
    ! Trapezoidal integration for enthalpy vs temperature, including phase transfer discontinuity  
    ! Enthalpy is phi = int_0^T rho(T')*Cp(T') dT' 
    ! Perhaps use larger number of data points for integration than for table. 

    ! Table increment
    phiT_table_dT = phiT_table_max_T/phiT_table_n_points  ! Equal increment to end temperature

    ! Properties at zero temperature
    temp_table(0) = 0_amrex_real  ! Temperature table 
    call get_rho(temp_table(0),rho) 
    call get_Cp(temp_table(0),Cp)  
    rhocp_i = rho*Cp ! product of density and heat capacity at temperature 
    enth_table(0) = 0_amrex_real  ! Enthalpy table

    ! Fill phiT table
    do i = 1,phiT_table_n_points                 

       temp_table(i) = temp_table(i-1) + phiT_table_dT  ! Update temperature 
       
       ! Solid-liquid phase transfer in temp-enth table 
       if (temp_table(i).ge.melt_point) then
          
          if (isolid) then  ! First data point above melting temperature  

             isolid = .false.  ! 
             imelt = i         ! Index of melt data point 
             temp_table(i) = melt_point             ! Exactly pin-point melt temperature in table 
             phiT_table_dT = melt_point - temp_table(i-1)   ! 
             rhocp_im1 = rhocp_i 
             call get_rho(temp_table(i),rho) 
             call get_Cp(temp_table(i),Cp)  
             rhocp_i = rho*Cp ! product of density and heat capacity at temperature  
             enth_table(i) = enth_table(i-1) + (rhocp_i+rhocp_im1)*phiT_table_dT/2_amrex_real   ! Enthalpy at melt onset 
             enth_at_melt = enth_table(i) 
             phiT_table_dT = (phiT_table_max_T - melt_point)/(phiT_table_n_points-1-i)  ! New phiT_table_dT to match phiT_table_max_T !
             
          end if
          
          if(imelt.eq.i-1) then ! solid-liquid phase transfer jump 

             temp_table(i) = melt_point ! 
             enth_table(i) = enth_table(i-1) + 1E6*enth_fus*rho_melt/m_A ! [J/m3] = 1E6*[kJ/mol]*[kg/m3]/[g/mol]
             
          elseif(imelt.ne.i) then
             
             rhocp_im1 = rhocp_i ! Product of density and heat capacity at previous temperature 
             call get_rho(temp_table(i),rho) 
             call get_Cp(temp_table(i),Cp)  
             rhocp_i = rho*Cp ! product of density and heat capacity at temperature       
             enth_table(i) = enth_table(i-1) + (rhocp_i+rhocp_im1)*phiT_table_dT/2_amrex_real  ! Trapezoidal integration
             
          end if
          
       else
          
          rhocp_im1 = rhocp_i            ! Product of density and heat capacity at previous temperature 
          call get_rho(temp_table(i),rho) 
          call get_Cp(temp_table(i),Cp)  
          rhocp_i = rho*Cp ! product of density and heat capacity at temperature 
          enth_table(i) = enth_table(i-1) + (rhocp_i+rhocp_im1)*phiT_table_dT/2_amrex_real  ! Trapezoidal integration
          
       end if

       ! Diffusivity 
       call get_ktherm(temp_table(i),ktherm)
       diffus = ktherm/rhocp_i  
       if (diffus.gt.max_diffus) max_diffus=diffus ! Find maximum thermal diffusivity for time step determination
       ! Is continually re-evaluated to correspond to maximum within domain. This evaluation adds comp. load - perhaps 
       ! better to use max value at all times? (i.e. gain of slightly larger step does not outweigh extra load to find diff)
       
   
    end do
		
    ! Output employed material properties to file
    open (2, file = 'Employed_'//TRIM(material)//'_properties.txt', status = 'unknown') 
     write(2,*) 'Material properties employed' 
     write(2,*) 'Temperature, Cp [J/kgK], rho [kg/m3], k [W/mk]' 
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
                      uo_lo, uo_hi, phi, & 
                      t_lo , t_hi , temp) 

    ! Input and output variables
    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: uo_lo(3), uo_hi(3)
    integer         , intent(in   ) :: t_lo(3), t_hi(3)
    real(amrex_real), intent(in   ) :: phi (uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3))  ! Enthalpy     [j/m3]
    real(amrex_real), intent(  out) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))  ! Temperature  [K]

    ! Local variables
    real(amrex_real) :: rho, Cp  ! Mass density [kg/m3], Specific heat capacity [j/kgK]
    integer :: i,j,k
    integer :: e_ind 
    real(amrex_real) :: phi_melt_onset, phi_melt_finish, enth_fus
    real(amrex_real) :: alph 
    ! Tabulated temp as function of enthalpy 
    ! Enthalpy is phi = int_0^T rho(T')*Cp(T') dT' 	
    
	
    ! Enthalpy given as integral phi = int_0^T rho(T')*Cp(T') dT'. When initializing problem, perform integral and create table 
    ! Get temp interpolates temp(phi) for each point in the lo,hi box.    
    do i = lo(1),hi(1)
       do j = lo(2),hi(2)
          do k = lo(3),hi(3)
             !Linear interpolation of temperature by given enthalpy, from table 
             do e_ind = 0,phiT_table_n_points 
                if (phi(i,j,k) .le. enth_table(e_ind) ) exit 
             end do
             if (e_ind.eq.phiT_table_n_points) STOP 'Temperature table exceeded' 
             alph = (phi(i,j,k)-enth_table(e_ind-1))/(enth_table(e_ind)-enth_table(e_ind-1))
             temp(i,j,k) = temp_table(e_ind-1) + alph*(temp_table(e_ind)-temp_table(e_ind-1))
             !print *, i, j, k, phi(i,j,k), temp(i,j,k)
          end do
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
    real(amrex_real) :: alph 

    !Linear interpolation of temperature by given enthalpy, from table
    do e_ind = 0,phiT_table_n_points 
       if (temp .le. temp_table(e_ind) ) exit 
    end do
    if (e_ind.eq.phiT_table_n_points) STOP 'Temperature table exceeded' 
    alph = (temp-temp_table(e_ind-1))/(temp_table(e_ind)-temp_table(e_ind-1))
    enth = enth_table(e_ind-1) + alph*(enth_table(e_ind)-enth_table(e_ind-1))
     
  end subroutine get_enthalpy
	
	
	
	
  ! subroutine to find maximum diffusivity for variable timestep 
  ! If used, called at timestep_mod, subr. increment after second get_temp. diff has to be reset before call to timestep (main)
  subroutine get_maxdiffus(lo, hi, 		 &
       t_lo , t_hi , temp) 
    integer         , intent(in) :: lo(3), hi(3)				
    integer         , intent(in) :: t_lo(3), t_hi(3)			
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))  ! Temperature  [K] 
    integer :: i,j,k 			
    real(amrex_real) :: rho, Cp, ktherm, diffus 
    
    do   i = lo(1),hi(1)
       do  j = lo(2),hi(2) 
  	  do k = lo(3),hi(3)
            call get_ktherm(temp(i,j,k),ktherm)
            call get_rho(temp(i,j,k),rho) 
            call get_Cp(temp(i,j,k),Cp) 
            diffus = ktherm/(rho*Cp) 
            if (diffus.gt.max_diffus) max_diffus = diffus 
         end do
      end do
   end do

 end subroutine get_maxdiffus

	


end module material_properties_module 
