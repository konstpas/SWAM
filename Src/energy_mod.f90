module energy_module 

   use amrex_amr_module
   
   implicit none
   private 
   
   public :: sum_enthalpy
   
contains 

   subroutine sum_enthalpy(enth_tot) 
   use amr_data_module, only : phi_new
   real(amrex_real), intent(out) :: enth_tot			! total enthalpy (output)
   
   
   ! integrate enthalpy in solid and liquid domain 
   enth_tot = phi_new(0)%sum(1) * amrex_geom(0)%dx(1)*amrex_geom(0)%dx(2)
   
   end subroutine sum_enthalpy


end module energy_module 
