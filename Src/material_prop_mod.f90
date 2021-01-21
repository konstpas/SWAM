module material_properties_module  
  use amrex_amr_module
  
	implicit none 

	! Parameters known only here such as material interface position 
	! String for chosen material 
	! 


	public :: get_ktherm, get_temp 


	contains 
	
	subroutine get_ktherm(temp,ktherm) 
	 real(amrex_real), intent(in   ) :: temp   ! Temperature [K]
	 real(amrex_real), intent(  out) :: ktherm ! Thermal conductivity [W/mK] 
	
	 ! Tabulate ktherm as function of temperature 
	
	 ! test case below (tungsten)
	 ktherm  = 100_amrex_real 
	
	end subroutine get_ktherm 
	
	
	
	
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
	 real(amrex_real) :: phi_melt, enth_fus
	 ! Tabulate temp as function of enthalpy 
	 ! Enthalpy is phi = int_0^T rho(T')*Cp(T') dT' 	
	
	!  test case below (tungsten)
	rho = 17e3           ! kg/m3 
	Cp  = 270_amrex_real ! j/kgK
	
	! Enthalpy of fusion is 52.3 kJ/mol, at 183.84 g/mol, 0.28448651 kJ/g 
	! 284.48651 J/g, 284.48651E3 j/kg
	
	
	!temp = phi/(rho*Cp)
   	do   i = lo(1),hi(1)
  	 do  j = lo(2),hi(2) 
  	  do k = lo(3),hi(3)
  	  	
		temp(i,j,k) = phi(i,j,k)/(rho*Cp) 
  	  end do    						! 'boundary volumetric' source
  	 end do 
  	end do 
	
	end subroutine 


end module material_properties_module 
