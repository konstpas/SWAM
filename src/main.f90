program SWAM

  use amrex_amr_module 
  use amr_data_module 
  use initdata_module 
  use my_amr_module 
  use timestep_module
  
  implicit none
  
  call my_amr_init  
  call initdata 

  call run_simulation
  
  call amr_data_finalize 
  call amrex_amrcore_finalize 
  call amrex_finalize 
 
end program SWAM         
