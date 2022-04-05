module energy_module 

  ! -----------------------------------------------------------------
  ! This module is used to keep track of the energy balance
  ! -----------------------------------------------------------------
  
   use amrex_amr_module
   
   implicit none
   private 

   ! ------------------------------------------------------------------
   ! Public subroutines
   ! ------------------------------------------------------------------
   public :: sum_enthalpy
   
contains 

  ! ------------------------------------------------------------------
  ! Subroutine to compute the total enthalpy in the domain
  ! ------------------------------------------------------------------
  subroutine sum_enthalpy(enth_tot) 

     use amr_data_module, only : phi_new

     real(amrex_real), intent(out) :: enth_tot
   
     ! phi_new contains the enthalpy per unit volume
     enth_tot = phi_new(0)%sum(1) * (amrex_geom(0)%dx(1)* & 
                                     amrex_geom(0)%dx(2)* & 
                                     amrex_geom(0)%dx(3))
     
   end subroutine sum_enthalpy


end module energy_module 
