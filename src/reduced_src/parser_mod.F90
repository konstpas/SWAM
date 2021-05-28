module parser

  use flap
  
  implicit none

  private

  public :: parse_command_line
  
contains 

  subroutine parse_command_line(input_file,buffer)

    integer intent(in) :: buffer = 256
    character(len=buffer), intent(inout) :: input_file
    type(command_line_interface):: cli
    integer :: cli_error
  end  subroutine parse_command_line

  
end module input
