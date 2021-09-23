program SWAM

  use omp_lib
  use amrex_amr_module 
  use init_module, only : run_init, run_finalize
  use simulation_module, only : run_simulation
  
  implicit none

  real(amrex_real) :: time_start
  real(amrex_real) :: time_end

  ! Start timing
  time_start = omp_get_wtime()

  ! Initialize data
  call run_init

  ! Run simulation
  call run_simulation

  ! Finalize data
  call run_finalize

  ! Stop timing
  time_end = omp_get_wtime()
  
  ! Print end of simulation message on screen
  print *, "Simulation completed. Elasped time", time_end - time_start, " seconds"
  
end program SWAM         
