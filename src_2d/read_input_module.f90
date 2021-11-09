module read_input_module

  ! -----------------------------------------------------------------
  ! This module is used to read the input file. When reading the
  ! input file it is important to realize that there are two type of
  ! variables:
  !
  !     1) Variables that are used by amrex to generate the mesh.
  !        These variables are labelled with "amr" and with 
  !        "geometry" in the input file. The value of such
  !        variables is read from the amrex initialization routines
  !        routines, so in the following module they are used only
  !        to ensure self-consistency of the input. Among these
  !        variables, the following must be specified:
  !            a) geometry.is_periodic
  !            b) geometry.coord_sys
  !            c) geometry.prob_lo
  !            d) geometry.prob_hi
  !            e) amr.n_cell
  !            f) amr.max_level
  !        The following can be omitted, in which case they will
  !        assume the default value:
  !            a) amr.v (default: 0)
  !            b) amr.ref_ratio (default: 2 for each level)
  !            c) amr.blocking_factor (default: 8)
  !            d) amr.max_grid_size (default: 32)
  !
  !     2) Variables that are used from the rest of the code. Such
  !        variables have a diverse set of labels and are read
  !        by the routine specified in this module. All these
  !        variables are optional, if not specified they will
  !        be assigned the values given in set_default_values
  ! -----------------------------------------------------------------
  
  use amrex_amr_module

  implicit none 

  private

  ! -----------------------------------------------------------------
  ! Public variables (read from input file)
  ! -----------------------------------------------------------------
  public :: cfl
  public :: check_file
  public :: check_int
  public :: cooling_debug
  public :: cooling_thermionic
  public :: cooling_vaporization
  public :: cooling_radiation
  public :: do_reflux
  public :: dt_change_max
  public :: plasma_flux_type
  public :: plasma_flux_params
  public :: plasma_input_file
  public :: material
  public :: max_grid_size_1d
  public :: max_step
  public :: meltvel
  public :: phase_init
  public :: phiT_table_max_T
  public :: phiT_table_n_points
  public :: plot_file
  public :: plot_int
  public :: regrid_int
  public :: restart
  public :: solve_sw
  public :: stop_time
  public :: surfdist
  public :: surf_pos_init
  public :: temp_fs
  public :: temp_init
  public :: thermionic_alpha
  public :: verbose

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: read_input_file, deallocate_input

  ! -----------------------------------------------------------------
  ! Default values of public variables
  ! -----------------------------------------------------------------
  character(len=:), allocatable, save :: check_file
  character(len=:), allocatable, save :: material
  character(len=:), allocatable, save :: plasma_flux_type
  character(len=:), allocatable, save :: plasma_input_file
  character(len=:), allocatable, save :: phase_init
  character(len=:), allocatable, save :: plot_file
  character(len=:), allocatable, save :: restart   
  integer, save :: check_int
  integer, save :: max_grid_size_1d
  integer, save :: max_step
  integer, save :: phiT_table_n_points
  integer, save :: plot_int
  integer, save :: regrid_int
  integer, save :: verbose
  logical, save :: cooling_thermionic
  logical, save :: cooling_vaporization
  logical, save :: cooling_radiation
  logical, save :: do_reflux
  logical, save :: solve_sw
  real(amrex_real), save :: cfl
  real(amrex_real), save :: dt_change_max  
  real(amrex_real), save :: meltvel
  real(amrex_real), save :: phiT_table_max_T
  real(amrex_real), save :: stop_time
  real(amrex_real), save :: surf_pos_init
  real(amrex_real), save :: temp_fs
  real(amrex_real), save :: temp_init
  real(amrex_real), save :: thermionic_alpha
  real(amrex_real), allocatable, save :: cooling_debug(:)
  real(amrex_real), allocatable, save :: surfdist(:)
  real(amrex_real), allocatable, save :: plasma_flux_params(:)
  
contains

  ! ------------------------------------------------------------------
  ! Subroutine used to read the input file. Note: part of the input
  ! file is also read by the amrex initialization routines invoked
  ! with amrex_init (see the initialization module)
  ! ------------------------------------------------------------------
  subroutine read_input_file()

    type(amrex_parmparse) :: pp

    ! Default parameters
    call set_default_values

    ! Parameters for the simulation length
    call amrex_parmparse_build(pp, "length")
    call pp%query("max_step", max_step)
    call pp%query("stop_time", stop_time)
    call amrex_parmparse_destroy(pp)

    ! Parameters for the grid control
    call amrex_parmparse_build(pp, "grid")
    call pp%query("regrid_int", regrid_int)
    call pp%query("max_grid_size_1d", max_grid_size_1d)
    call pp%getarr("surfdist", surfdist)
    call amrex_parmparse_destroy(pp)

    ! Parameters for the output
    call amrex_parmparse_build(pp, "output")
    call pp%query("check_int", check_int)
    call pp%query("plot_int", plot_int)
    call pp%query("check_file", check_file)
    call pp%query("plot_file", plot_file)
    call pp%query("verbose", verbose)
    call amrex_parmparse_destroy(pp)

    ! Parameters for the restart
    call amrex_parmparse_build(pp, "restart")
    call pp%query("restart", restart)
    call amrex_parmparse_destroy(pp)
   
    ! Parameters for the heat solver
    call amrex_parmparse_build(pp, "heat")
    call pp%query("surf_pos", surf_pos_init)   
    call pp%query("meltvel", meltvel)  
    call pp%query("temp_init", temp_init)
    call pp%query("phase_init", phase_init)
    call pp%query("plasma_flux_type", plasma_flux_type) 
    call pp%getarr("plasma_flux_params", plasma_flux_params)
    call pp%query("plasma_input_file", plasma_input_file)
    call pp%query("temp_free_surface", temp_fs)
    call pp%query("cooling_thermionic",cooling_thermionic)
    call pp%query("cooling_vaporization",cooling_vaporization)
    call pp%query("cooling_radiation",cooling_radiation)
    call pp%getarr("cooling_debug",cooling_debug)
    call pp%query("thermionic_alpha",thermionic_alpha)
    call amrex_parmparse_destroy(pp)

    ! Parameters for the shallow waters solver
    call amrex_parmparse_build(pp, "sw")
    call pp%query("solve", solve_sw)
    call amrex_parmparse_destroy(pp)

    ! Parameters for the material
    call amrex_parmparse_build(pp, "material")
    call pp%query("material", material)
    call pp%query("phiT_max_T", phiT_table_max_T)
    call pp%query("phiT_n_points", phiT_table_n_points)
    call amrex_parmparse_destroy(pp)
    
    ! Parameters for the numerics
    call amrex_parmparse_build(pp, "numerics")
    call pp%query("cfl", cfl)
    call pp%query("do_reflux", do_reflux)
    call pp%query("dt_change_max", dt_change_max)
    call amrex_parmparse_destroy(pp)
	    
  end subroutine read_input_file


  ! Default values for the input parameters
  subroutine set_default_values()

    integer :: i
        
    allocate(character(len=3)::check_file)
    allocate(character(len=8)::plasma_flux_type)
    allocate(character(len=25)::plasma_input_file)
    allocate(character(len=8)::material)
    allocate(character(len=9)::phase_init)
    allocate(character(len=3)::plot_file)
    allocate(character(len=0)::restart)
    allocate(cooling_debug(5))    
    allocate(surfdist(0:amrex_max_level))
    allocate(plasma_flux_params(100))

    cfl = 0.70
    check_file = "chk"
    check_int = -1
    cooling_debug(1) = 0
    cooling_debug(2) = 300
    cooling_debug(3) = 301
    cooling_debug(4) = 1
    cooling_debug(5) = 1e6
    cooling_thermionic = .true.
    cooling_vaporization = .true.
    cooling_radiation = .true.
    do_reflux = .true.
    dt_change_max = 1.1
    plasma_flux_params(1) = 0.0
    plasma_flux_params(2) = 1.0
    plasma_flux_params(3) = 300e6
    plasma_flux_params(4) = 0.0
    plasma_flux_params(5) = 0.01
    plasma_flux_type = "Gaussian"
    plasma_input_file = "Incoming_heat.txt"
    material = "Tungsten"
    max_grid_size_1d = 16
    max_step = 10000
    phase_init = "undefined"
    phiT_table_max_T = 10000.0
    phiT_table_n_points = 10000
    plot_file = "plt"
    plot_int = -1
    regrid_int = 2
    solve_sw = .true.
    stop_time = 1.0
    do i = 0, amrex_max_level
       surfdist(i) = 0.0
    end do
    surf_pos_init = 0.020
    thermionic_alpha = 90.0
    temp_fs = -1.0
    verbose = 0
    
  end subroutine set_default_values

  ! ------------------------------------------------------------------
  ! Subroutine used to free the memory related to the input variables
  ! ------------------------------------------------------------------  
  subroutine deallocate_input()
    
    deallocate(check_file)
    deallocate(material)
    deallocate(phase_init)
    deallocate(plasma_flux_params)
    deallocate(plasma_flux_type)
    deallocate(plasma_input_file)
    deallocate(plot_file)
    deallocate(restart)
    deallocate(surfdist)

  end subroutine deallocate_input
    
end module read_input_module
