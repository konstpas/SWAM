module material_properties_module  
  use amrex_amr_module
  
	implicit none 

	! Parameters known only here such as material interface position 
	! String for chosen material 
	! 

	public :: init_mat_prop, get_ktherm, get_temp, get_maxdiffus
	public :: enth_at_melt, melt_point, max_diffus  

	real(amrex_real), allocatable	:: temp_table(:), enth_table(:)
	real(amrex_real) 		:: m_A, melt_point, enth_fus, rhomelt 
	real(amrex_real)		:: enth_at_melt ! enthalpy at onset of melting  
	real(amrex_real) 		:: max_diffus ! Maximum thermal diffusivity 
	character(len=:), allocatable	:: material		! material name
	integer 			:: ntemp      		! ntemp temperature points for table
	
	contains 
	
	!---------------------------------
	subroutine get_ktherm(temp,ktherm) 
	 real(amrex_real), intent(in   ) :: temp   ! Temperature [K]
	 real(amrex_real), intent(  out) :: ktherm ! Thermal conductivity [W/mK] 
	
	 ! Outputs thermal conductivity as a function of temperature 

	 ktherm  = 100_amrex_real ! Constant from input file 
	 
	 
	 
	! Tungsten thermal conductivity ---------------------------
	if (material.eq.'Tungsten') then 
	 ! Temperature dependent properties known 
	 ! Nuclear Materials and Energy 13 (2017) 42-57 
	 
	 !W/(m K)
	  if (temp.lt.300.0) then ! lower edge value 
	      	ktherm = 179.9041  
	  else if(temp.lt.3695.0) then 
	 	ktherm = 149.441 - 45.466E-3*temp + 13.193E-6*temp**2 - 1.484E-9*temp**3 + 3.866E6/temp**2 
	  else if(temp.le.6000.0) then 
	 	ktherm = 66.6212 + 0.02086*(temp-3695.0) - 3.7585E-6*(temp-3695.0)**2
	  else 
	  	ktherm = 95.7345 
	  end if
	 end if  
	 !---------------------------------------------------------

	end subroutine get_ktherm 
	!---------------------------------
	
	!---------------------------------
	subroutine get_rho(temp,rho) 
	 real(amrex_real), intent(in   ) :: temp   ! Temperature [K]
	 real(amrex_real), intent(  out) :: rho    ! Mass density [kg/m3] 
	 real(amrex_real) :: temp_0 = 293.15 
	 ! Outputs mass density as a function of temperature
	
	 rho  = 17e3 ! Constant from input file 
	 
	 
	 
	! Tungsten mass density ---------------------------
	if (material.eq.'Tungsten') then 
	 ! Temperature dependent properties known 
	 ! Nuclear Materials and Energy 13 (2017) 42-57 
	  if (temp.lt.300.0) then ! lower edge value 
	      	rho = 19.2482 
	  else if(temp.lt.3695.0) then 
	 	rho = 19.25 - 2.66207E-4*(temp-temp_0)-3.0595E-9*(temp-temp_0)**2 - 9.5185E-12*(temp-temp_0)**3
	  else if(temp.le.6000.0) then 
	 	rho = 16.267 - 7.679E-4*(temp-3695.0) - 8.091E-8*(temp-3695.0)**2 
	  else 
	  	rho = 14.0671 
	  end if
	  rho = rho*1e3  ! g/cm3 to kg/m3 
	 end if  
	 !---------------------------------------------------------
	
	end subroutine get_rho
	!---------------------------------
	
	!---------------------------------
	subroutine get_Cp(temp,Cp) 
	 real(amrex_real), intent(in   ) :: temp   ! Temperature [K]
	 real(amrex_real), intent(  out) :: Cp     ! Specific heat capacity [J/kgK] 
	
	 ! Tabulate Cp as function of temperature 
	
	 Cp = 270_amrex_real ! Constant from input file 
	 
	 
	 
	 ! Tungsten specific heat capacity ---------------------------
	 if (material.eq.'Tungsten') then 
	  ! Temperature dependent properties known 
	  ! Nuclear Materials and Energy 13 (2017) 42-57 
	  if (temp.lt.300.0) then ! lower edge value 
	      	Cp = 24.1363 
	  else if(temp.lt.3080.0) then 
	 	Cp = 21.868372 + 8.068661E-3*temp - 3.756196E-6*temp**2 + 1.075862E-9*temp**3 + 1.406637E4/temp**2 
	  else if(temp.le.3695.0) then 
	 	Cp = 2.022 + 1.315E-2*temp 
	  else 
	  	Cp = 51.3 
	  end if
	  Cp = 1E3*(Cp/m_A)! J/(mol K) to J/(kg K) 
	 end if  
	 !---------------------------------------------------------
	 
	 
	 
	
	end subroutine get_Cp 
	!---------------------------------
	
	
	
  
  	! Initializes table for temperature as function of enthalpy 
  	subroutine init_mat_prop()
  	real(amrex_real) :: end_temp, dtemp 	! end of temperature interval (decided in input), and temperature increment size  
  	real(amrex_real) :: rho, Cp, ktherm 
  	integer :: i, imelt 
  	real(amrex_real) :: rhocp_i, rhocp_im1, diffus   
  	logical :: isolid = .true.  ! for enthalpy table, true before phase transfer  
  	
  	max_diffus = 0. 
  	
  	! Reads material properties values 
  	! Values constant unless known temperature dependence for 
  	! material is written in this subroutine 
  	open (1, file = TRIM(material)//'.dat', status = 'old') 
  		read(1,*) !material properties 
  		read(1,*) ! Atomic mass [g/mol]
  		read(1,*) m_A
  		read(1,*) ! Melting point [K]
  		read(1,*) melt_point
  		read(1,*) !Enthalpy of fusion [kJ/mol], density at melting [kg/m3] 
  		read(1,*) enth_fus, rhomelt 
  		read(1,*) ! End temperature [K] and number of data points [], used for table of temperature as function of enthalpy 
  		read(1,*) end_temp, ntemp 
  	close(1) 
  	

	! Tungsten thermophysical properties as function of temperature 
	if (material=='Tungsten') then 
  	 ! give end temperature point and data points in input file? 
  	 ! now called several times, move call to initdata later 
	 allocate(temp_table(0:ntemp))
	 allocate(enth_table(0:ntemp)) 
	
	
	! Trapezoidal integration for enthalpy vs temperature, including phase transfer discontinuity  
	! Enthalpy is phi = int_0^T rho(T')*Cp(T') dT' 
	! Perhaps use larger number of data points for integration than for table. 
	 dtemp = end_temp/ntemp  ! Equal increment to end temperature 
	 temp_table(0) = 0_amrex_real  ! Temperature table 
	 
	 call get_rho(temp_table(0),rho) 
	 call get_Cp(temp_table(0),Cp)  
	 rhocp_i = rho*Cp ! product of density and heat capacity at temperature 
	 enth_table(0) = 0_amrex_real  ! Enthalpy table 
	 do i = 1,ntemp                ! Integration do loop 
	      temp_table(i) = temp_table(i-1) + dtemp  ! Update temperature 
	      
	       ! Solid-liquid phase transfer in temp-enth table 
	       if (temp_table(i).ge.melt_point) then
	       
	        	if (isolid) then  ! First data point above melting temperature  
	        		isolid = .false.  ! 
	       		imelt = i         ! Index of melt data point 
	        		temp_table(i) = melt_point             ! Exactly pin-point melt temperature in table 
	      			dtemp = melt_point - temp_table(i-1)   ! 
	      	  		 rhocp_im1 = rhocp_i 
	      	  		   call get_rho(temp_table(i),rho) 
	 			   call get_Cp(temp_table(i),Cp)  
	 			   rhocp_i = rho*Cp ! product of density and heat capacity at temperature  
	      	  		 enth_table(i) = enth_table(i-1) + (rhocp_i+rhocp_im1)*dtemp/2_amrex_real   ! Enthalpy at melt onset 
	      	  		 enth_at_melt = enth_table(i) 
	        		 dtemp = (end_temp - melt_point)/(ntemp-1-i)  ! New dtemp to match end_temp ! 
	       	end if 	
	       
	       	if(imelt.eq.i-1) then ! solid-liquid phase transfer jump 
	       		temp_table(i) = melt_point ! 
	       		enth_table(i) = enth_table(i-1) + 1E6*enth_fus*rhomelt/m_A ! j/m3 = 1E6*[kJ/mol]*[kg/m3]/[g/mol]
	       	elseif(imelt.ne.i) then  
	      			rhocp_im1 = rhocp_i            ! Product of density and heat capacity at previous temperature 
	      	  		   call get_rho(temp_table(i),rho) 
	 			   call get_Cp(temp_table(i),Cp)  
	 			   rhocp_i = rho*Cp ! product of density and heat capacity at temperature       
	      			enth_table(i) = enth_table(i-1) + (rhocp_i+rhocp_im1)*dtemp/2_amrex_real  ! Trapezoidal integration
	       	end if 
	       	
	       else 
	      		rhocp_im1 = rhocp_i            ! Product of density and heat capacity at previous temperature 
	      	  		   call get_rho(temp_table(i),rho) 
	 			   call get_Cp(temp_table(i),Cp)  
	 			   rhocp_i = rho*Cp ! product of density and heat capacity at temperature 
	      		enth_table(i) = enth_table(i-1) + (rhocp_i+rhocp_im1)*dtemp/2_amrex_real  ! Trapezoidal integration 	
	      end if 
	  
	  call get_ktherm(temp_table(i),ktherm)    
	  diffus = ktherm/rhocp_i  
	  if (diffus.gt.max_diffus) max_diffus=diffus ! Find maximum thermal diffusivity for time step determination
	  ! Is continually re-evaluated to correspond to maximum within domain. This evaluation adds comp. load - perhaps 
	  ! better to use max value at all times? (i.e. gain of slightly larger step does not outweigh extra load to find diff)
	    
	      
	 end do 
	! End of trapezoidal integration 
	
	
	! Output employed material properties in table 
	open (2, file = 'Employed_'//TRIM(material)//'_properties.txt', status = 'unknown') 
		write(2,*) 'Material properties employed' 
		write(2,*) 'Temperature, Cp [J/kgK], rho [kg/m3], k [W/mk]' 
		do i = 0,ntemp
		call get_Cp(temp_table(i),Cp) 
		call get_ktherm(temp_table(i),ktherm) 
		call get_rho(temp_table(i),rho)
		write(2,*) temp_table(i), Cp, rho, ktherm
		end do  
  	close(2) 
	
	
	
	
	end if 



  	end subroutine init_mat_prop
  
  
	
	
	subroutine get_temp(lo, hi, 		 &
				uo_lo, uo_hi, phi, & 
				t_lo , t_hi , temp) 
	 integer         , intent(in   ) :: lo(3), hi(3)
	 integer         , intent(in   ) :: uo_lo(3), uo_hi(3)
	 integer         , intent(in   ) :: t_lo(3), t_hi(3)
	 real(amrex_real), intent(in   ) :: phi (uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3))  ! Enthalpy     [j/m3]
	 real(amrex_real), intent(  out) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))  ! Temperature  [K] 	 
	 real(amrex_real) :: rho, Cp  ! Mass density [kg/m3], Specific heat capacity [j/kgK]
	 integer :: i,j,k
	 integer :: e_ind 
	 real(amrex_real) :: phi_melt_onset, phi_melt_finish, enth_fus
	 real(amrex_real) :: alph 
	 ! Tabulated temp as function of enthalpy 
	 ! Enthalpy is phi = int_0^T rho(T')*Cp(T') dT' 	
	
	
	! Enthalpy given as integral phi = int_0^T rho(T')*Cp(T') dT'. When initializing problem, perform integral and create table 
	! Get temp interpolates temp(phi) for each point in the lo,hi box. 
	
	
   	do   i = lo(1),hi(1)
  	 do  j = lo(2),hi(2) 
  	  do k = lo(3),hi(3)
  	  	! Linear interpolation of temperature by given enthalpy, from table 
  	  	do e_ind = 0,ntemp 
  	  		if (phi(i,j,k) .le. enth_table(e_ind) ) exit 
  	  	end do 
  	  		if (e_ind.eq.ntemp) STOP 'Temperature table exceeded' 
			alph = (phi(i,j,k)-enth_table(e_ind-1))/(enth_table(e_ind)-enth_table(e_ind-1))
			temp(i,j,k) = temp_table(e_ind-1) + alph*(temp_table(e_ind)-temp_table(e_ind-1))
  	  end do    						
  	 end do 
  	end do 
	
	end subroutine get_temp
	
	
	
	
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
