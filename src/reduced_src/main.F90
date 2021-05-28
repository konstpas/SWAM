program SWAM

  use iso_fortran_env
  use amrex_amr_module
  use flap
  use strings
  
  implicit none

  ! Amrex variables
  
  ! real(amrex_real) :: cur_time
  ! real(amrex_real) :: total_enthalpy     
  ! integer :: last_plot_file_step, step, lev, substep, finest_level
  
  ! FLAP variables
  type(command_line_interface):: cli
  integer :: cli_error
  integer, parameter :: buffer = 256
  character(len=buffer) :: input_file
  
  ! Open and read file variables
  integer, parameter :: read_unit = 99, n_tokens_max = 256
  character(len=buffer) :: line
  character(len=buffer), dimension(n_tokens_max) :: tokens
  integer :: ios, n_lines, n_tokens, ii
  
  ! Data structure to store the input
  type input_data

     ! Simulation duration
     integer :: max_step
     real(amrex_real) :: stop_time

     ! Problem size and geometry
     real(amrex_real) :: prob_lo_x, prob_lo_y, prob_lo_z
     real(amrex_real) :: prob_hi_x, prob_hi_y, prob_hi_z
     integer :: n_cell_x, n_cell_y, n_cell_z

     ! Domain
     real(amrex_real) :: melt_vel
     real(amrex_real) :: temp_init
     real(amrex_real) :: surf_pos
     real(amrex_real) :: surf_dist_x, surf_dist_y, surf_dist_z
     character(len=buffer) :: material
     real(amrex_real) :: flux_peak
     real(amrex_real) :: flux_pos_x, flux_pos_z
     real(amrex_real) :: flux_width_x, flux_width_z
     real(amrex_real) :: exposure_time
     
     ! Refinement
     integer :: max_level
     integer :: ref_ratio
     integer :: blocking_factor
     integer :: max_grid_size
     integer :: max_grid_size_2d
     integer :: regrid_interval

     ! Plot
     character(len=buffer) :: plot_file
     integer :: plot_interval

     ! Verbosity and extra properties
     integer :: verbose
     real(amrex_real) :: cfl
     integer :: do_reflux
     real(amrex_real) :: phi_err_x, phi_err_y, phi_err_z
     
  end type input_data

  type(input_data) :: in
  
  ! ---- Command line parsing (move to dedicated file) ----

  ! Initializing Command Line Interface
  call cli%init(progname    = 'SWAM',      &
                version     = '1.0',       &
                authors     = 'Fingerborg',&                              
                license     = '---',       &                              
                description = '---')
  
  ! Setting Command Line Argumenst
  call cli%add(switch='--input',  & 
               switch_ab='-in',   &
               help='Input file', &
               required=.true.,   &
               act='store',       &
               error=cli_error)

  ! Parsing Command Line Interface
  call cli%parse(error=cli_error)
  if (cli_error/=0) stop

  ! Assign Command Line Arguments to variables
  call cli%get(switch='--input', val=input_file, error=cli_error) ; if (cli_error/=0) stop

  ! Sanity check on the input data
  if (len(trim(input_file)) == buffer) then
        
     print *, "Input file name is too long, maximum", buffer, "characters are allowed"
     close(read_unit)
     stop
        
  end if
  
  ! ---- Command line parsing (move to dedicated file) ----

  ! ---- Read input file (move to dedicated file) ----

  ! Assign default values
  in%max_step = 1000000
  in%stop_time = 1.0
  in%prob_lo_x = 0.0
  in%prob_lo_y = 0.0
  in%prob_lo_z = 0.0
  in%prob_hi_x = 0.02
  in%prob_hi_y = 0.025
  in%prob_hi_z = 0.02
  in%n_cell_x = 16
  in%n_cell_y = 64
  in%n_cell_z = 16
  in%melt_vel = 10.0
  in%temp_init = 1.0
  in%surf_pos = 0.020
  in%surf_dist_x = 0.001
  in%surf_dist_y = 0.0005
  in%surf_dist_z = 0.0001
  in%material = "Tungsten"
  in%flux_peak = 300e6
  in%flux_pos_x = 0.0
  in%flux_pos_z = 0.0
  in%flux_width_x = 1.01
  in%flux_width_z = 1.01
  in%exposure_time = 0.4  
  in%max_level = 1
  in%ref_ratio = 1
  in%blocking_factor = 8
  in%max_grid_size = 16
  in%max_grid_size_2D = 16
  in%regrid_interval = 2
  in%plot_file = "plt"
  in%plot_interval = 10
  in%verbose = 1
  in%cfl = 0.95
  in%do_reflux = 1
  in%phi_err_x = 1.01
  in%phi_err_y = 1.1
  in%phi_err_z = 1.2

  ! Open file
  open(unit=read_unit, file=input_file, iostat=ios)
  if ( ios /= 0 ) then
     print *, "Error opening the input file: ", input_file
     stop
  end if
  
  n_lines = 0
  do
     
     ! Read line
     read(read_unit, '(A)', iostat=ios) line
     if ( ios /= 0 .and. ios /= iostat_end) stop "Unexpected error while reading input file"
  
     ! Update line counter
     n_lines = n_lines + 1

     ! Check that the line fits the buffer
     if (len(trim(line)) == buffer) then
        
        print *, "Risk of buffer overflow while reading input line number: ", n_lines
        close(read_unit)
        stop
        
     end if

     ! Parse line
     call parse(line,' ', tokens, n_tokens)

     ! Read command (skip comments and empty lines)
     if ( tokens(1) /= "#" .or. n_tokens /= 0) then
        print *, "---"
     end if

     ! Stop at the end of file
     if (ios == iostat_end) exit

 
  end do

  ! Close file
  close(read_unit)

  
  ! ---- Read input file (move to dedicated file) ----
  
  ! ---- Amrex ----
  ! call amrex_init 
  ! call amrex_amrcore_init 

  ! call amrex_amrcore_finalize 
  ! call amrex_finalize 
  ! ---- Amrex ----
  
end program SWAM         
