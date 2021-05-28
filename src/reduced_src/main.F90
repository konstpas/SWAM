program SWAM

  use amrex_amr_module
  use parser
  use input
  
  implicit none

  ! Amrex variables
  
  ! real(amrex_real) :: cur_time
  ! real(amrex_real) :: total_enthalpy     
  ! integer :: last_plot_file_step, step, lev, substep, finest_level
  
  ! Variables for command line parsing

  integer, parameter :: buffer = 256
  character(len=buffer) :: input_file
  
  ! Open and read file variables
  type(in_data) :: input_data
  
  ! Parse command line
  call parse_command_line(input_file,buffer)
  
  ! Read input file
  call read_input(input_file,input_data)

    
  ! ---- Amrex ----
  ! call amrex_init 
  ! call amrex_amrcore_init 

  ! call amrex_amrcore_finalize 
  ! call amrex_finalize 
  ! ---- Amrex ----
  
end program SWAM         
