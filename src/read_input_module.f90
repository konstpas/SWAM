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
  !            c) amr.blocking_facto (default: 8)
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

  ! Public variables (read from input file)
  public :: cfl
  public :: check_file
  public :: check_int
  public :: do_reflux
  public :: exp_time
  public :: flux_peak
  public :: flux_pos
  public :: flux_width
  public :: material
  public :: max_grid_size_2d
  public :: max_step
  public :: meltvel
  public :: phiT_table_max_T
  public :: phiT_table_n_points
  public :: plot_file
  public :: plot_int
  public :: regrid_int
  public :: restart
  public :: stop_time
  public :: surfdist
  public :: surf_pos_init
  public :: tempinit
  public :: verbose
  
  ! Public subroutines
  public :: read_input_file
  
  ! Default values of public variables
  character(len=:), allocatable, save :: check_file
  character(len=:), allocatable, save :: material
  character(len=:), allocatable, save :: plot_file
  character(len=:), allocatable, save :: restart
  integer :: check_int
  integer :: max_grid_size_2d
  integer :: max_step
  integer :: phiT_table_n_points
  integer :: plot_int
  integer :: regrid_int
  integer :: verbose
  logical :: do_reflux
  real(amrex_real) :: cfl
  real(amrex_real) :: exp_time
  real(amrex_real) :: flux_peak
  real(amrex_real) :: meltvel
  real(amrex_real) :: phiT_table_max_T
  real(amrex_real) :: stop_time
  real(amrex_real) :: surf_pos_init
  real(amrex_real) :: tempinit
  real(amrex_real), allocatable, save :: surfdist(:)
  real(amrex_real), allocatable, save :: flux_pos(:)
  real(amrex_real), allocatable, save :: flux_width(:)
  
contains

  ! Subroutine used to update the default value of the public variables 
  subroutine read_input_file()

    type(amrex_parmparse) :: pp
    integer :: ilev
    
    ! Default parameters
    call set_default_values
    
    ! Read parameters
    call amrex_parmparse_build(pp, "length")
    call pp%query("max_step", max_step)
    call pp%query("stop_time", stop_time)
    call amrex_parmparse_destroy(pp)
    
    ! Parameters amr.*
    call amrex_parmparse_build(pp, "amr")
    call pp%query("regrid_int", regrid_int)
    call pp%query("check_int", check_int)
    call pp%query("plot_int", plot_int)
    call pp%query("check_file", check_file)
    call pp%query("plot_file", plot_file)
    call pp%query("restart", restart)
    call pp%query("max_grid_size_2d", max_grid_size_2d)
    call amrex_parmparse_destroy(pp)
   
    ! Domain parameters 
    call amrex_parmparse_build(pp, "domain")
    call pp%query("meltvel", meltvel)  
    call pp%query("tempinit", tempinit) 
    call pp%query("surf_pos", surf_pos_init)   
    call pp%getarr("surfdist", surfdist)
    call pp%query("flux_peak", flux_peak) 
    call pp%getarr("flux_pos", flux_pos)
    call pp%getarr("flux_width", flux_width)
    call pp%query("exp_time", exp_time) 
    call amrex_parmparse_destroy(pp)

    ! Material parameters
    call amrex_parmparse_build(pp, "material")
    call pp%query("material", material)
    call pp%query("phiT_max_T", phiT_table_max_T)
    call pp%query("phiT_n_points", phiT_table_n_points)
    call amrex_parmparse_destroy(pp)
    
    ! Parameters myamr.*
    call amrex_parmparse_build(pp, "myamr")
    call pp%query("v", verbose)
    call pp%query("verbose", verbose)
    call pp%query("cfl", cfl)
    call pp%query("do_reflux", do_reflux)
    call amrex_parmparse_destroy(pp)
	
 
  end subroutine read_input_file


  ! Default values for the input parameters
  subroutine set_default_values()

    allocate(character(len=3)::check_file)
    allocate(character(len=8)::material)
    allocate(character(len=3)::plot_file)
    allocate(character(len=0)::restart)
    allocate(surfdist(2))
    allocate(flux_pos(2))
    allocate(flux_width(2))
    
    cfl = 0.70
    check_file = "chk"
    check_int = -1
    do_reflux = .true.
    exp_time = 1.0
    flux_peak = 300d6
    flux_pos = (/0.0,0.0/)
    flux_width = (/1.01,1.01/)
    material = "tungsten"
    max_grid_size_2d = 16
    max_step = 10000
    phiT_table_max_T = 10000.0
    phiT_table_n_points = 10000
    plot_file = "plt"
    plot_int = -1
    regrid_int = 2
    stop_time = 1.0
    surfdist = (/0.001,0.0005/)
    surf_pos_init = 0.020
    verbose = 0

    
  end subroutine set_default_values

  
end module read_input_module
