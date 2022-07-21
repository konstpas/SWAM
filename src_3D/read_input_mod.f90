module read_input_module

! -----------------------------------------------------------------
  ! This module is used to read the input file. When reading the
  ! input file it is important to realize that there are two type of
  ! variables:
  !
  ! 1) Variables that are used by amrex to generate the mesh.
  !    These variables are labelled with "amr" and with 
  !    "geometry" in the input file. The value of such
  !     variables is read from the amrex initialization routines
  !    so in the following module they are used only
  !    to ensure self-consistency of the input.
  !
  !    Among these variables, the following must be specified:
  !     a) geometry.prob_lo: Coordinates of the bottom-left corner of
  !                          the domain (x,y). Distances are in meters
  !     b) geometry.prob_hi: Coordinates of the top-right corner of
  !                          the domain (x,y). Distances are in meters
  !     c) amr.n_cell: Number of cells in each direction (x,y)
  !     d) amr.max_level: Maximum amrex level (starting from 0)
  !
  !    The following can be omitted, in which case they will
  !    assume the default value:
  !     a) amr.v (default: 0)
  !     b) amr.ref_ratio (default: 2 for each level)
  !     c) amr.blocking_factor (default: 8)
  !     d) amr.max_grid_size (default: 32)
  !
  ! 2) Variables that are used from the rest of the code. Such
  !    variables have a diverse set of labels and are read
  !    by the routine specified in this module. All these
  !    variables are optional, if not specified they will
  !    be assigned the values given in set_default_values
  ! -----------------------------------------------------------------
  
  use amrex_amr_module
  use amrex_linear_solver_module
  
  implicit none 

  private

  ! -----------------------------------------------------------------
  ! Public variables (read from input file)
  ! -----------------------------------------------------------------

  ! --- Variables that control the length of the simulation ---

  ! Time step [s]
  public :: time_dt
  ! Maximum number of time steps
  public :: time_step_max  
  ! Maximum simulated time [s]
  public :: time_max
  ! Maximum fractional change of the timestep between two
  ! successive timesteps
  public :: time_ddt_max

  ! --- Variables for the geometry of the heat transfer domain ---

  ! Select type of geometry to simulate (Slab, West or West_rectangular)
  public :: geom_name
  ! West geometry: Coordinates for the center of the cooling pipe 
  public :: geom_cool_pipe_cntr
  ! West geometry: Radius of the cooling pipe [m]
  public :: geom_cool_pipe_radius

  ! --- Variables that control the grid generated by amrex ---

  ! How often to check if regrid is necessary (in number of time steps)
  public :: regrid_int
  ! Distance from the free surface that triggers regridding [m].
  ! One value for each level larger than 0 must be provided
  public :: regrid_dist
  ! Free surface deformation that triggers regridding [m/s].
  ! One value for each level larger than 0 must be provided 
  public :: regrid_def 
  
  ! --- Variables for the heat solver ---

  ! Parameters to plot the cooling fluxes as a function of
  ! the temperature. Array with 5 elements: 1) on-off switch,
  ! with 1 = on and 0 = off, 2) minimum temperature
  ! 3) maximum temperature, 4) temperature resolution
  ! 5) incoming heat flux
  public :: heat_cooling_debug
  ! Activate or deactivate thermionic cooling (1 = on, 0 = off)
  public :: heat_cooling_thermionic
  ! Activate or deactivate thermionic cooling on the exposed
  ! side in leading edge geometries (1 = on, 0 = off)
  public :: heat_cooling_thermionic_side
  ! Activate or deactivate vaporization cooling (1 = on, 0 = off)
  public :: heat_cooling_vaporization
  ! Activate or deactivate radiative cooling (1 = on, 0 = off)
  public :: heat_cooling_radiation
  ! Species density to be used in the limited thermionic current formula [m^-3]
  public :: heat_density 
  ! Mass of ions used in the limited thermionic current formula [kg]
  public :: heat_ion_mass
  ! Sheath transmission coefficient
  public :: heat_sheath_coeff
  ! Solver for the heat equation (explicit or implicit)
  public :: heat_solver
  ! Initial phase of the sample (solid, liquid). Relevant only
  ! if the sample is initialized at the melting temperature
  public :: heat_phase_init
  ! Type of plasma heat flux (Gaussian, Uniform, Input_file)
  public :: heat_plasma_flux_type
  ! Parameters for the plasma heat flux. An array with 5 parameters
  ! whose content depends on the type of plasma heat flux
  !  a) Gaussian heat flux, array with 5 elements: 1) switching on,
  !     time [s], 2) switching off time [s], 3) peak magnitude [W/m^2],
  !     4) peak position [m], 5) standard deviation [m]
  !  b) Uniform heat flux, array with 5 elements: 1) switching on
  !     time [s], 2) switching off time [s], 3) magnitude [W/m^2],
  !     4) minimum x coordinate [m], 5) maximum x coordinate [m]
  public :: heat_plasma_flux_params
  ! Name of file containing the plasma flux (used only if the plasma
  ! heat flux is of type Input_file)
  public :: heat_plasma_flux_file
  ! Type of plasma heat flux on the second exposed side for the
  ! West geometry (Gaussian, Uniform, Input_file)
  public :: heat_plasma_flux_side_type
  ! Parameters for the plasma heat flux on the second exposed side
  ! for the West geometry. Analogous to heat_plasma_flux_parameters
  public :: heat_plasma_flux_side_params
  ! Name of file containing the plasma flux on the second exposed
  ! side for the West geometry (Used only if the plasma
  ! heat flux is of type Input_file)
  public :: heat_plasma_flux_side_file
  ! Activate or deactivate reflux in the explicit heat
  ! solver (1 = on, 0 = off)
  public :: heat_reflux
  ! Edge of the simulation domain along the x direction [m]
  public :: heat_sample_edge
  ! Solve the heat equation (1 = yes, 0 = no) 
  public :: heat_solve
  ! Impose temperature on the free surface (only used if positive) [K]
  public :: heat_temp_surf
  ! Initial uniform temperature of the sample [K]
  public :: heat_temp_init
  ! Flag to control if the parallel heat-flux should be scaled with the local
  ! surface normals calculated by the code 1=on, 0=off(default)
  public :: heat_local_surface_normals
  ! Ratio of the parallel heat flux that satisfies the optical approximation
  ! see Nuclear Materials and Energy 17 (2018) 194-199 E. Thoren et al.
  public :: heat_Foa
  ! Prescribed a temperature at the bottom of the sample (only used if positive) [K]
  public :: heat_temp_bottom

  ! --- Variables for the shallow water solver ---

  ! Capping height for the viscous drag in the momentum equation [m]
  public :: sw_captol
  ! Thermionic current on the free surface [A/m^2]. An array of 3
  ! elements containing: 1) the magnitude of the current,
  !  2) the switching on time, 3) the switching off time.
  ! This is only used if the heat equation is not solved,
  ! i.e. if heat.solve = 0
  public :: sw_current
  ! Dry tolerance for the melt pool [m]  (any height below dry
  ! tolerance is assumed to be zero)
  public :: sw_drytol
  ! Inclination of the magnetic field [degrees]
  public :: sw_magnetic_inclination
  ! The angle of attack of the magnetic field w.r.t. the exposed side in WEST geometry [degrees]
  public :: sw_side_magnetic_inclination
  ! Magnitude of the magnetic field [T]
  public :: sw_Bx
  public :: sw_Bz
  ! Magnitude of the acceleration due to gravity [m/s^2]
  public :: sw_gx
  public :: sw_gz
  ! Fixed melt velocity [m/s]. Only used if the momentum equation
  ! is not solve, i.e. sw_solve_momentum = 0
  public :: sw_melt_velocity
  ! Parameters used to generate a gaussian pool if the heat equation
  ! is not solved. An array with 3 elements containing: 1) The
  ! maximum depth of the pool [m], 2) the position of the center of the
  ! pool [m], 3) the width of the pool (as a standard deviation) [m]
  public :: sw_pool_params
  ! Solve the shallow water equation (1 = yes, 0 = no) 
  public :: sw_solve
  ! Initial position of the free surface [m] (in the heat domain,
  ! assumed perpendicular to the y-direction)
  public :: sw_surf_pos_init
  ! Solve the momentum equation (1 = yes, 0 = no). If the momentum
  ! equation is not solved, the velocity is assumed to be equal to
  ! the velocity prescribed with sw.sw_melt_velocity
  public :: sw_solve_momentum
  ! Solver for the shallow water equations (geoclaw or explicit)
  public :: sw_solver
  ! Number of maximum iterations in the calculation of the advective update.
  ! Used only if the sw_solver is geoclaw.
  public :: sw_iter
  ! Acceleration due to gravity. Used only if the sw_solver is geoclaw [m/s^2]
  public :: sw_gravity
  ! Include Marangoni flow when solving shallow water equation. Four elements,
  ! first for positive x-direction, second for negative x-direction, third positive
  ! z-direction and fourth for negative z-direction. If an element if on (1) the 
  ! marangoni contribution in that direction is taken into account. If an element is
  ! off (0) the Marangoni contribution in that direction is neglected. 
  public :: sw_marangoni
  ! Capping of the 1/h pre-factor to the Marangoni contribution at the momentum
  ! equation of the shallow water equations [m]
  public :: sw_marang_cap
  ! Flag to control whether the free surface is read from a file (1 = yes, 0 = no) 
  public :: sw_read_free_surface_file
  ! Name of file that containt the coordinates for the free surface
  ! (used only if the read_free_surface_file is on)
  public :: sw_free_surface_file
  ! The direction of the incoming heat-flux
  public :: sw_B_unity
  ! A prefactor to the derivative of the surface tension with respect to the temperature
  ! due to possible uncertainties at the quantity
  public :: sw_surf_tension_deriv_prefactor

  ! --- Variables for the material properties ---

  ! Type of material to simulated, only one can be selected (Tungsten,
  ! Beryllium, Iridium, Niobium)
  public :: material_name
  ! Maximum temperature used to compute the material properties [K]
  public :: material_maxT
  ! Number of temperatures at which the material properties are computed
  public :: material_nPoints


  ! --- Variables that control the numerics ---

  ! Safety factor used to multiply the minimum time step
  ! that ensures the stability of the heat solverw
  public :: num_cfl
  ! Relative error used in the linear solvers (only used if the implicit
  ! solver for the heat equation is employed)
  public :: num_accuracy
  ! Internal parameters for the linear solvers (only used if the implicit
  ! solver for the heat equation is employed)
  public :: num_composite_solve
  public :: num_verbose
  public :: num_bottom_verbose
  public :: num_max_iter
  public :: num_max_fmg_iter
  public :: num_bottom_solver
  public :: num_linop_maxorder
  public :: num_agglomeration
  public :: num_consolidation
  public :: num_max_coarsening_level
  ! Activate or deactivate subcycling in time (1 = on, 0 = off) 
  public :: num_subcycling

  
  ! --- Variables for the input/output ---

  ! Name of the restart files written by amrex
  public :: io_check_file
  ! Name of the plot files written by amrex
  public :: io_plot_file
  ! How often to write restart files (in time steps)
  public :: io_check_int
  ! How often to write output to files (in time steps)
  public :: io_plot_int
  ! Name of file use to restart the simulation
  public :: io_restart
  ! Verbosity of the code (1 = on, 0 = off)
  public :: io_verbose
  


  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: read_input_file, deallocate_input

  ! -----------------------------------------------------------------
  ! Default values of public variables
  ! -----------------------------------------------------------------

  ! Length of simulation
  integer, save :: time_step_max
  real(amrex_real), save :: time_dt
  real(amrex_real), save :: time_ddt_max
  real(amrex_real), save :: time_max

  ! Geometry
  character(len=:), allocatable, save :: geom_name
  real(amrex_real), allocatable, save :: geom_cool_pipe_cntr(:)
  real(amrex_real), save :: geom_cool_pipe_radius

  ! Grid
  integer, save :: regrid_int
  real(amrex_real), allocatable, save :: regrid_dist(:)
  real(amrex_real), allocatable, save :: regrid_def(:)

  ! Heat solver
  character(len=:), allocatable, save :: heat_solver
  character(len=:), allocatable, save :: heat_plasma_flux_type
  character(len=:), allocatable, save :: heat_plasma_flux_file
  character(len=:), allocatable, save :: heat_plasma_flux_side_type
  character(len=:), allocatable, save :: heat_plasma_flux_side_file
  character(len=:), allocatable, save :: heat_phase_init
  logical, save :: heat_cooling_thermionic
  logical, save :: heat_cooling_thermionic_side
  logical, save :: heat_cooling_vaporization
  logical, save :: heat_cooling_radiation
  logical, save :: heat_reflux
  logical, save :: heat_solve
  logical, save :: heat_local_surface_normals
  real(amrex_real), save :: heat_sample_edge
  real(amrex_real), save :: heat_temp_surf
  real(amrex_real), save :: heat_temp_bottom
  real(amrex_real), save :: heat_temp_init
  real(amrex_real), save :: heat_Foa
  real(amrex_real), save :: heat_density
  real(amrex_real), save :: heat_ion_mass
  real(amrex_real), save :: heat_sheath_coeff
  real(amrex_real), allocatable, save :: heat_cooling_debug(:)
  real(amrex_real), allocatable, save :: heat_plasma_flux_params(:)
  real(amrex_real), allocatable, save :: heat_plasma_flux_side_params(:)

  ! Shallow water solver
  logical, save :: sw_solve
  logical, save :: sw_solve_momentum
  logical, allocatable, save :: sw_marangoni(:)
  logical, save :: sw_read_free_surface_file
  real(amrex_real), save :: sw_captol
  real(amrex_real), save :: sw_drytol
  real(amrex_real), save :: sw_marang_cap
  real(amrex_real), save :: sw_magnetic_inclination
  real(amrex_real), save :: sw_side_magnetic_inclination
  real(amrex_real), save :: sw_Bx
  real(amrex_real), save :: sw_Bz
  real(amrex_real), save :: sw_gx
  real(amrex_real), save :: sw_gz
  real(amrex_real), save :: sw_gravity
  real(amrex_real), save :: sw_surf_tension_deriv_prefactor
  real(amrex_real), allocatable, save :: sw_melt_velocity(:)
  real(amrex_real), save :: sw_surf_pos_init
  real(amrex_real), allocatable, save :: sw_current(:)
  real(amrex_real), allocatable, save :: sw_pool_params(:)
  real(amrex_real), allocatable, save :: sw_B_unity(:)
  character(len=:), allocatable, save :: sw_solver
  character(len=:), allocatable, save :: sw_free_surface_file
  integer, save :: sw_iter

  ! Material properties
  character(len=:), allocatable, save :: material_name
  integer, save :: material_nPoints
  real(amrex_real), save :: material_maxT

  ! Numerics
  integer, save :: num_verbose
  integer, save :: num_bottom_verbose
  integer, save :: num_max_iter
  integer, save :: num_max_fmg_iter
  integer, save :: num_bottom_solver
  integer, save :: num_linop_maxorder
  logical, save :: num_composite_solve  
  logical, save :: num_agglomeration
  logical, save :: num_consolidation
  logical, save :: num_subcycling
  integer, save :: num_max_coarsening_level
  real(amrex_real), save :: num_cfl
  real(amrex_real), save :: num_accuracy

  ! Input/output
  character(len=:), allocatable, save :: io_check_file
  character(len=:), allocatable, save :: io_plot_file
  character(len=:), allocatable, save :: io_restart
  integer, save :: io_check_int
  integer, save :: io_plot_int
  integer, save :: io_verbose
  
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
    call amrex_parmparse_build(pp, "time")
    call pp%query("max_step", time_step_max)
    call pp%query("stop_time", time_max)
    call pp%query("dt", time_dt)
    call pp%query("dt_change_max", time_ddt_max)
    call amrex_parmparse_destroy(pp)

    ! Parameters for the geometry
    call amrex_parmparse_build(pp, "geometry")
    call pp%query("name", geom_name)
    call pp%queryarr("cool_pipe_cntr",geom_cool_pipe_cntr)
    call pp%query("cool_pipe_radius", geom_cool_pipe_radius)
    call amrex_parmparse_destroy(pp)

    ! Parameters for the grid
    call amrex_parmparse_build(pp, "regrid")
    call pp%query("interval", regrid_int)
    call pp%queryarr("distance", regrid_dist)
    call pp%queryarr("deformation", regrid_def)
    call amrex_parmparse_destroy(pp)

    ! Parameters for the heat solver
    call amrex_parmparse_build(pp, "heat")
    call pp%query("solve",heat_solve) 
    call pp%queryarr("sw_melt_velocity", sw_melt_velocity)  
    call pp%query("sample_edge", heat_sample_edge)   
    call pp%query("temp_init", heat_temp_init)
    call pp%query("local_surface_normals", heat_local_surface_normals)
    call pp%query("optical_ratio", heat_Foa)
    call pp%query("phase_init", heat_phase_init)
    call pp%query("plasma_flux_type", heat_plasma_flux_type) 
    call pp%query("plasma_flux_side_type", heat_plasma_flux_side_type) 
    call pp%queryarr("plasma_flux_side_params", heat_plasma_flux_side_params)
    call pp%queryarr("plasma_flux_params", heat_plasma_flux_params)
    call pp%query("plasma_flux_file", heat_plasma_flux_file)
    call pp%query("plasma_flux_side_file", heat_plasma_flux_side_file)
    call pp%query("temp_free_surface", heat_temp_surf)
    call pp%query("temp_bottom_surface", heat_temp_bottom)
    call pp%query("cooling_thermionic",heat_cooling_thermionic)
    call pp%query("side_cooling_thermionic",heat_cooling_thermionic_side)
    call pp%query("cooling_vaporization",heat_cooling_vaporization)
    call pp%query("cooling_radiation",heat_cooling_radiation)
    call pp%queryarr("cooling_debug",heat_cooling_debug)
    call pp%query("solver",heat_solver)
    call pp%query("reflux", heat_reflux)
    call pp%query("density", heat_density)
    call pp%query("ion_mass", heat_ion_mass)
    call pp%query("sheath_transmission_coeff", heat_sheath_coeff)
    call amrex_parmparse_destroy(pp)

    ! Parameters for the shallow water solver
    call amrex_parmparse_build(pp, "sw")
    call pp%queryarr("melt_velocity", sw_melt_velocity) 
    call pp%query("surf_pos_init", sw_surf_pos_init) 
    call pp%query("solve", sw_solve)
    call pp%query("solve_momentum", sw_solve_momentum)
    call pp%query("marangoni_x_positive", sw_marangoni(1))
    call pp%query("marangoni_x_negative", sw_marangoni(2))
    call pp%query("marangoni_z_positive", sw_marangoni(3))
    call pp%query("marangoni_z_negative", sw_marangoni(4))
    call pp%query("marangoni_cap", sw_marang_cap)
    call pp%query("read_free_surface_file", sw_read_free_surface_file)
    call pp%query("magnetic_inclination",sw_magnetic_inclination)
    call pp%query("side_magnetic_inclination",sw_side_magnetic_inclination)
    call pp%query("captol", sw_captol)
    call pp%query("drytol", sw_drytol)
    call pp%queryarr("current", sw_current)
    call pp%queryarr("pool_params", sw_pool_params) 
    call pp%queryarr("B_unity", sw_B_unity) 
    call pp%query("Bx", sw_Bx)
    call pp%query("Bz", sw_Bz)
    call pp%query("gx", sw_gx)
    call pp%query("gz", sw_gz)
    call pp%query("surf_tension_deriv_prefactor", sw_surf_tension_deriv_prefactor)
    call pp%query("solver",sw_solver)
    call pp%query("free_surface_file",sw_free_surface_file)
    call pp%query("geoclaw_iter", sw_iter)
    call pp%query("geoclaw_gravity", sw_gravity)
    call amrex_parmparse_destroy(pp)
    if( abs(sqrt(sw_B_unity(1)**2+sw_B_unity(3)**2+sw_B_unity(3)**2)-1.0).gt.1E-6) then
      print *, 'Magnetic field unity vector does not have a lengh of 1. Instead it has a length of ', &
        sqrt(sw_B_unity(1)**2+sw_B_unity(3)**2+sw_B_unity(3)**2)
      print *, 'Renormalizing'
      sw_B_unity = sw_B_unity/(sqrt(sw_B_unity(1)**2+sw_B_unity(3)**2+sw_B_unity(3)**2)) 
    end if
    
    ! Parameters for the numerics
    call amrex_parmparse_build(pp, "numerics")
    call pp%query("cfl", num_cfl)
    call pp%query("accuracy", num_accuracy)
    call pp%query("composite_solve", num_composite_solve)  
    call pp%query("verbose", num_verbose)
    call pp%query("bottom_verbose_ls", num_bottom_verbose)
    call pp%query("max_iter", num_max_iter)
    call pp%query("max_fmg_iter", num_max_fmg_iter)
    call pp%query("bottom_solver", num_bottom_solver)
    call pp%query("linop_maxorder", num_linop_maxorder)
    call pp%query("agglomeration", num_agglomeration)
    call pp%query("consolidation", num_consolidation)
    call pp%query("max_coarsening_level", num_max_coarsening_level)
    call pp%query("subcycling", num_subcycling)
    call amrex_parmparse_destroy(pp)

    ! Parameters for the material
    call amrex_parmparse_build(pp, "material")
    call pp%query("name", material_name)
    call pp%query("max_temperature", material_maxT)
    call pp%query("points", material_nPoints)
    call amrex_parmparse_destroy(pp)

    ! Parameters for the input and output
    call amrex_parmparse_build(pp, "io")
    call pp%query("check_int", io_check_int)
    call pp%query("plot_int", io_plot_int)
    call pp%query("check_file", io_check_file)
    call pp%query("plot_file", io_plot_file)
    call pp%query("verbose", io_verbose)
    call pp%query("restart", io_restart)
    call amrex_parmparse_destroy(pp)
    
  end subroutine read_input_file


  ! ------------------------------------------------------------------
  ! Subroutines used to set the default values for the input variables
  ! ------------------------------------------------------------------
  subroutine set_default_values()

    ! Allocate input variables that are not scalars
    call allocate_input
    
    ! Assign default values
    call set_default_time
    call set_default_geom
    call set_default_regrid
    call set_default_heat
    call set_default_sw
    call set_default_material
    call set_default_numerics
    call set_default_io
     
  end subroutine set_default_values


  subroutine set_default_time()
    
    time_ddt_max = 1.1
    time_dt = 1e-6
    time_step_max = 0
    time_max = 0.0
    
  end subroutine set_default_time


  subroutine set_default_geom()

    real(amrex_real) :: Lx
    real(amrex_real) :: Ly

    Lx = amrex_probhi(1) - amrex_problo(1)
    Ly = amrex_probhi(2) - amrex_problo(2)
    
    geom_name = "Slab"
    geom_cool_pipe_cntr(1) = 0.5*Lx
    geom_cool_pipe_cntr(2) = 0.5*Ly
    geom_cool_pipe_radius = 0.0

  end subroutine set_default_geom


  subroutine set_default_regrid()

    regrid_int = 1
    regrid_dist = 0.0_amrex_real
    regrid_def = 0.0_amrex_real
    
  end subroutine set_default_regrid


  subroutine set_default_heat()    
    
    heat_cooling_debug(1) = 0
    heat_cooling_debug(2) = 300
    heat_cooling_debug(3) = 301
    heat_cooling_debug(4) = 1
    heat_cooling_debug(5) = 1e6
    heat_cooling_thermionic = .true.
    heat_cooling_vaporization = .true.
    heat_cooling_radiation = .true.
    heat_reflux = .true.    
    heat_cooling_thermionic_side = .false.
    heat_local_surface_normals = .false.
    heat_solver = "explicit"
    heat_phase_init = "undefined"
    heat_plasma_flux_file = "plasma_flux.dat"
    heat_Foa = 0.0
    heat_plasma_flux_params(1) = 0.0
    heat_plasma_flux_params(2) = 1.0
    heat_plasma_flux_params(3) = 300E6
    heat_plasma_flux_params(4) = 0.0
    heat_plasma_flux_params(5) = 0.01
    heat_plasma_flux_params(6) = 0.0
    heat_plasma_flux_params(7) = 0.01
    heat_plasma_flux_side_params(1) = 0.0
    heat_plasma_flux_side_params(2) = 1.0
    heat_plasma_flux_side_params(3) = 300e6
    heat_plasma_flux_side_params(4) = 0.0
    heat_plasma_flux_side_params(5) = 0.01
    heat_plasma_flux_type = "Gaussian"
    heat_plasma_flux_side_file = "plasma_side_flux.dat"
    heat_plasma_flux_side_params = heat_plasma_flux_params
    heat_plasma_flux_side_type = heat_plasma_flux_type
    heat_solve = .true.
    heat_sample_edge = 0.020
    heat_temp_surf = -1.0
    heat_temp_bottom = -1.0
    heat_density = 2E19
    heat_ion_mass = 3.34358377E-27
    heat_sheath_coeff = 7
    
  end subroutine set_default_heat


  subroutine set_default_sw()

    real(amrex_real) :: Ly

    Ly = amrex_probhi(2) - amrex_problo(2)

    sw_captol = 0.0
    sw_current = 0.0
    sw_drytol = 0.0    
    sw_marang_cap = 0.0
    sw_magnetic_inclination = 90.0
    sw_side_magnetic_inclination = 90.0
    sw_Bx = 0.0
    sw_Bz = 0.0
    sw_gx = 0.0
    sw_gz = 0.0
    sw_melt_velocity = 0.0
    sw_B_unity(1) = 0.0
    sw_B_unity(2) = 1.0
    sw_B_unity(3) = 0.0
    sw_pool_params(1) = 0.0
    sw_pool_params(2) = 0.0
    sw_pool_params(3) = 1.0
    sw_solve = .true.
    sw_solve_momentum = .true.
    sw_marangoni(1) = .true.
    sw_marangoni(2) = .true.
    sw_marangoni(3) = .true.
    sw_marangoni(4) = .true.
    sw_surf_tension_deriv_prefactor = 1.0
    sw_read_free_surface_file = .false.
    sw_surf_pos_init = 0.5*Ly
    ! sw_surf_pos_init = 0.006*Ly
    sw_iter = 1000
    sw_gravity = 1.0_amrex_real
    sw_solver = "explicit"
    sw_free_surface_file = "free_surface.dat"
    
  end subroutine set_default_sw


  subroutine set_default_material()
    
    material_name = "Tungsten"
    material_maxT = 10000.0
    material_nPoints = 10000
    
  end subroutine set_default_material


  subroutine set_default_numerics()

    num_cfl = 0.95
    num_accuracy = 10e-10
    num_composite_solve = .true.
    num_verbose = 0
    num_bottom_verbose = 0
    num_max_iter = 100
    num_max_fmg_iter = 0
    num_bottom_solver = amrex_bottom_default
    num_linop_maxorder = 2
    num_agglomeration = .true.
    num_consolidation = .true.
    num_max_coarsening_level = 30
    num_subcycling = .false.
    
  end subroutine set_default_numerics

  subroutine set_default_io()

    io_check_file = "chk"
    io_check_int = -1
    io_plot_file = "plt"
    io_plot_int = -1
    io_verbose = 0
    
  end subroutine set_default_io

  ! ------------------------------------------------------------------
  ! Subroutines used to allocate memory relate to the input variables
  ! ------------------------------------------------------------------
  subroutine allocate_input()

    call allocate_geom_variables
    call allocate_regrid_variables
    call allocate_heat_variables
    call allocate_sw_variables
    call allocate_material_variables
    call allocate_io_variables
    
  end subroutine allocate_input


  subroutine allocate_geom_variables()

    allocate(character(len=4)::geom_name)
    allocate(geom_cool_pipe_cntr(2))
    
  end subroutine allocate_geom_variables


  subroutine allocate_regrid_variables()
    
    allocate(regrid_dist(1:amrex_max_level))
    allocate(regrid_def(1:amrex_max_level))
        
  end subroutine allocate_regrid_variables


  subroutine allocate_heat_variables()

    allocate(character(len=8)::heat_solver)
    allocate(character(len=8)::heat_plasma_flux_type)
    allocate(character(len=8)::heat_plasma_flux_side_type)
    allocate(character(len=9)::heat_phase_init)
    allocate(character(len=25)::heat_plasma_flux_file)
    allocate(character(len=25)::heat_plasma_flux_side_file)
    allocate(heat_cooling_debug(5))
    allocate(heat_plasma_flux_params(100))
    allocate(heat_plasma_flux_side_params(100))
        
  end subroutine allocate_heat_variables


  subroutine allocate_sw_variables()

    allocate(sw_current(3))
    allocate(sw_pool_params(3))
    allocate(sw_B_unity(3))
    allocate(sw_melt_velocity(1:2))
    allocate(character(len=8)::sw_solver)  
    allocate(character(len=25)::sw_free_surface_file)   
    allocate(sw_marangoni(4)) 
        
  end subroutine allocate_sw_variables


  subroutine allocate_material_variables()

    allocate(character(len=8)::material_name)
        
  end subroutine allocate_material_variables


  subroutine allocate_io_variables()

    allocate(character(len=3)::io_check_file)
    allocate(character(len=3)::io_plot_file)
    allocate(character(len=0)::io_restart)
    
  end subroutine allocate_io_variables


  ! ------------------------------------------------------------------
  ! Subroutine used to free the memory related to the input variables
  ! ------------------------------------------------------------------  
  subroutine deallocate_input()
    
    deallocate(io_check_file)
    deallocate(heat_solver)
    deallocate(sw_solver)
    deallocate(sw_free_surface_file)
    deallocate(material_name)
    deallocate(geom_name)
    deallocate(heat_phase_init)
    deallocate(heat_plasma_flux_params)
    deallocate(heat_plasma_flux_type)
    deallocate(heat_plasma_flux_side_type)
    deallocate(heat_plasma_flux_file)
    deallocate(heat_plasma_flux_side_file)
    deallocate(io_plot_file)
    deallocate(io_restart)
    deallocate(regrid_dist)
    deallocate(sw_current)
    deallocate(sw_pool_params)
    deallocate(sw_B_unity)

  end subroutine deallocate_input
    
end module read_input_module
