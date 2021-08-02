program SWAM

  use amrex_amr_module 
  use init_module, only : run_init, run_finalize
  use simulation_module, only : run_simulation
  
  implicit none

  call run_init
  
  call run_simulation

  call run_finalize
 
end program SWAM         
