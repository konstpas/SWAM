program SWAM

  use amrex_amr_module
  use flap
  use input
  
  implicit none

  ! Amrex variables
  
  ! real(amrex_real) :: cur_time
  ! real(amrex_real) :: total_enthalpy     
  ! integer :: last_plot_file_step, step, lev, substep, finest_level
  
  ! Variables for command line parsing
  type(command_line_interface):: cli
  integer :: cli_error
  integer, parameter :: buffer = 256
  character(len=buffer) :: input_file
  
  ! Open and read file variables
  integer, parameter :: read_unit = 99, n_tokens_max = 256
  character(len=buffer) :: line
  character(len=buffer), dimension(n_tokens_max) :: tokens
  integer :: ios, n_lines, n_tokens, ii
  type(in_data) :: input_data
  
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

  ! Read input file
  call read_input(input_file,input_data)

    
  ! ---- Amrex ----
  ! call amrex_init 
  ! call amrex_amrcore_init 

  ! call amrex_amrcore_finalize 
  ! call amrex_finalize 
  ! ---- Amrex ----
  
end program SWAM         
