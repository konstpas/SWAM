module parser

  use flap
  
  implicit none

  private

  public :: parse_command_line
  
contains 

  subroutine parse_command_line(input_file,buffer)

    integer, intent(in) :: buffer
    character(len=buffer), intent(inout) :: input_file
    type(command_line_interface):: cli
    integer :: cli_error

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
       stop
    end if
    
  end  subroutine parse_command_line

  
end module parser
