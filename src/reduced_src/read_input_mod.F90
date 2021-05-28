module input

  use iso_fortran_env
  use amrex_amr_module
  use strings
  
  implicit none

  private

  public :: in_data
  public :: read_input
  
  ! String buffer
  integer, parameter :: buffer = 256
  
  ! Data structure to store the input
  type :: in_data

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
     
  end type in_data
  
contains 

  subroutine read_input(input_file,input_data)

    character(len=*), intent(in) :: input_file
    type(in_data), intent(inout) :: input_data
    integer :: fid = 99
    character(len=buffer) :: line
    character(len=buffer), dimension(buffer) :: tokens
    integer :: ios, n_lines, n_tokens, ii

    ! Assign default values
    call assign_default(input_data)

    ! Open input file
    open(unit=fid, file=input_file, iostat=ios)
    if ( ios /= 0 ) then
       print *, "Error opening the input file: ", input_file
       stop
    end if

    ! Read content of input file
    n_lines = 0
    do
     
       ! Read line
       read(fid, '(A)', iostat=ios) line
       if ( ios /= 0 .and. ios /= iostat_end) stop "Unexpected error while reading input file"
  
       ! Update line counter
       n_lines = n_lines + 1

       ! Check that the line fits the buffer
       if (len(trim(line)) == buffer) then
        
          print *, "Risk of buffer overflow while reading input line number: ", n_lines
          close(fid)
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
  close(fid)
    
  end subroutine read_input


  subroutine assign_default(input_data)

    type(in_data), intent(inout) :: input_data
    
    ! Default values for input data
    input_data%max_step = 1000000
    input_data%stop_time = 1.0
    input_data%prob_lo_x = 0.0
    input_data%prob_lo_y = 0.0
    input_data%prob_lo_z = 0.0
    input_data%prob_hi_x = 0.02
    input_data%prob_hi_y = 0.025
    input_data%prob_hi_z = 0.02
    input_data%n_cell_x = 16
    input_data%n_cell_y = 64
    input_data%n_cell_z = 16
    input_data%melt_vel = 10.0
    input_data%temp_init = 1.0
    input_data%surf_pos = 0.020
    input_data%surf_dist_x = 0.001
    input_data%surf_dist_y = 0.0005
    input_data%surf_dist_z = 0.0001
    input_data%material = "Tungsten"
    input_data%flux_peak = 300e6
    input_data%flux_pos_x = 0.0
    input_data%flux_pos_z = 0.0
    input_data%flux_width_x = 1.01
    input_data%flux_width_z = 1.01
    input_data%exposure_time = 0.4  
    input_data%max_level = 1
    input_data%ref_ratio = 1
    input_data%blocking_factor = 8
    input_data%max_grid_size = 16
    input_data%max_grid_size_2D = 16
    input_data%regrid_interval = 2
    input_data%plot_file = "plt"
    input_data%plot_interval = 10
    input_data%verbose = 1
    input_data%cfl = 0.95
    input_data%do_reflux = 1
    input_data%phi_err_x = 1.01
    input_data%phi_err_y = 1.1
    input_data%phi_err_z = 1.2

  end subroutine assign_default

  
end module input
