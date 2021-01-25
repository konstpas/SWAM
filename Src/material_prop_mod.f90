module material_properties_module  
  use amrex_amr_module
  
	implicit none 

	! Parameters known only here such as material interface position 
	! String for chosen material 
	! 

	public :: get_ktherm, get_temp 

	real(amrex_real), allocatable	:: temp_table(:), enth_table(:)
	real(amrex_real) 		:: m_A, melt_point, enth_fus 
	character(len=:), allocatable	:: material
	
	contains 
	
	subroutine get_ktherm(temp,ktherm) 
	 real(amrex_real), intent(in   ) :: temp   ! Temperature [K]
	 real(amrex_real), intent(  out) :: ktherm ! Thermal conductivity [W/mK] 
	
	 ! Tabulate ktherm as function of temperature 
	
	 ! test case below (tungsten)
	 ktherm  = 100_amrex_real 
	
	end subroutine get_ktherm 
	
  
  	! Initializes table for temperature as function of enthalpy 
  	subroutine init_mat_prop()
  	integer :: ntemp      		! ntemp temperature points for table
  	real(amrex_real) :: end_temp, dtemp 	! end of temperature interval (decided in input), and temperature increment size  
  	real(amrex_real) :: rho, Cp, rhomelt  
  	integer :: i, imelt 
  	real(amrex_real) :: rhocp_i, rhocp_im1  
  	logical :: isolid = .true.  ! for enthalpy table, true before phase transfer  
  	
  	! Reads material properties values 
  	! Values constant unless known temperature dependence for 
  	! material is written in this subroutine 
  	open (1, file = TRIM(material)//'.dat', status = 'old') 
  		read(1,*) !material properties 
  		read(1,*) ! Atomic mass [g/mol]
  		read(1,*) m_A
  		read(1,*) ! Melting point [K]
  		read(1,*) melt_point
  		read(1,*) ! Enthalpy of fusion [kJ/mol] 
  		read(1,*) enth_fus 
  		read(1,*) ! End temperature [K] and number of data points [], used for table of temperature as function of enthalpy 
  		read(1,*) end_temp, ntemp 
  	close(1) 
  	







	! Tungsten thermophysical properties as function of temperature 
	if (material=='Tungsten') then 
  	 ! give end temperature point and data points in input file? 
  	 ! now called several times, move call to initdata later 
	 allocate(temp_table(0:ntemp))
	 allocate(enth_table(0:ntemp)) 
	 print *, end_temp, ntemp 
	
	
	! Trapezoidal integration for enthalpy vs temperature, including phase transfer discontinuity  
	! Enthalpy is phi = int_0^T rho(T')*Cp(T') dT' 
	 dtemp = end_temp/ntemp  ! Equal increment to end temperature 
	 temp_table(0) = 0_amrex_real  ! Temperature table 
	 rhocp_i = 17e3*270_amrex_real ! product of density and heat capacity at temperature 
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
	      	  		 rhocp_i = 17e3*270_amrex_real ! Solid density/heat capacity product 
	      	  		 enth_table(i) = enth_table(i-1) + (rhocp_i+rhocp_im1)*dtemp/2_amrex_real   ! Enthalpy at melt onset 
	        		 dtemp = (end_temp - melt_point)/(ntemp-1-i)  ! New dtemp to match end_temp ! 
	       	end if 	
	       
	       	if(imelt.eq.i-1) then ! solid-liquid phase transfer jump 
	       		rhomelt = 17e3 ! Density of phase transfer (average of liquid and solid densities at melting point)
	       		temp_table(i) = melt_point ! 
	       		enth_table(i) = enth_table(i-1) + 1E6*enth_fus*rhomelt/m_A ! j/m3 = 1E6*[kJ/mol]*[kg/m3]/[g/mol]
	       	elseif(imelt.ne.i) then  
	      			rhocp_im1 = rhocp_i            ! Product of density and heat capacity at previous temperature 
	      			rhocp_i = 17e3*270_amrex_real  ! product of density and heat capacity at temperature      
	      			enth_table(i) = enth_table(i-1) + (rhocp_i+rhocp_im1)*dtemp/2_amrex_real  ! Trapezoidal integration
	       	end if 
	       	
	       else 
	      		rhocp_im1 = rhocp_i            ! Product of density and heat capacity at previous temperature 
	      		rhocp_i = 17e3*270_amrex_real  ! product of density and heat capacity at temperature      
	      		enth_table(i) = enth_table(i-1) + (rhocp_i+rhocp_im1)*dtemp/2_amrex_real  ! Trapezoidal integration 	
	      end if 
	      
	 end do 
	
	
	open (2, file = 'Enthalpy.dat', status = 'unknown') 
	do i = 1,ntemp
  		write(2,*) temp_table(i), enth_table(i)
  	end do 
  	close(2) 
	
	end if 
	
	
	 pause 


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
	 real(amrex_real) :: phi_melt_onset, phi_melt_finish, enth_fus
	 ! Tabulate temp as function of enthalpy 
	 ! Enthalpy is phi = int_0^T rho(T')*Cp(T') dT' 	
	
	!  test case below (tungsten)
	rho = 17e3           ! kg/m3 
	Cp  = 270_amrex_real ! j/kgK
	enth_fus = 284.4865E3 ! j/kg 
	phi_melt_onset = 3695_amrex_real*rho*Cp 
	phi_melt_finish = phi_melt_onset + enth_fus*rho  
	
	
	! Enthalpy of fusion is 52.3 kJ/mol, at 183.84 g/mol, 0.28448651 kJ/g 
	! 284.48651 J/g, 284.48651E3 j/kg
	
	! Enthalpy given as integral phi = int_0^T rho(T')*Cp(T') dT'. When initializing problem, perform integral and create table 
	! Get temp interpolates temp(phi) for each point in the lo,hi box. 
	
	call init_mat_prop 
	
	!temp = phi/(rho*Cp)
   	do   i = lo(1),hi(1)
  	 do  j = lo(2),hi(2) 
  	  do k = lo(3),hi(3)
  	  	if (phi(i,j,k).lt.phi_melt_onset) then 
		temp(i,j,k) = phi(i,j,k)/(rho*Cp) 
		elseif (phi(i,j,k).le.phi_melt_finish) then 
		temp(i,j,k) = phi_melt_onset/(rho*Cp) 
		else 
		temp(i,j,k) = (phi(i,j,k)-enth_fus*rho)/(rho*Cp) 
		end if 
  	  end do    						! 'boundary volumetric' source
  	 end do 
  	end do 
	
	end subroutine 


end module material_properties_module 
