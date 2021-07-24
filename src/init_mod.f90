module init_module

  ! -----------------------------------------------------------------
  ! This module is used to initialize and finalize the simulation
  ! -----------------------------------------------------------------
  
  use iso_c_binding
  use amrex_amr_module
  use amr_data_module
  
  implicit none

  private

  ! ------------------------------------------------------------------
  ! Public subroutines
  ! ------------------------------------------------------------------
  public :: run_init
  public :: run_finalize

contains

  
  ! -----------------------------------------------------------------
  ! Subroutine to initialize a simulation (from scratch or from
  ! restart file)
  ! -----------------------------------------------------------------
  subroutine run_init()
    
    use read_input_module
    use material_properties_module
    use regrid_module

    ! Initialize amrex 
    call amrex_init 
    call amrex_amrcore_init

    ! Read input file
    call read_input_file

    ! Initialize tables with material properties
    call init_mat_prop

    ! Initialize amrex data used in the simulation
    call amr_data_init
    
    ! Initialize amrex functions to generate grid
    call amrex_init_virtual_functions(my_make_new_level_from_scratch, &
                                      my_make_new_level_from_coarse, &
                                      my_remake_level, &
                                      my_clear_level, &
                                      my_error_estimate)

    ! Start new simulation or read restart file
    if (len_trim(restart) .eq. 0) then
       call amrex_init_from_scratch(0.0_amrex_real)
       call averagedown()
    else
       call amrex_abort("init from checkpoint not implemented yet")
    end if
    
  end subroutine run_init

  ! -----------------------------------------------------------------
  ! Subroutine to finalize a simulation and clean the memory
  ! -----------------------------------------------------------------
  subroutine run_finalize()

    ! Free amrex data
    call amr_data_finalize

    ! Finalize amrex
    call amrex_amrcore_finalize 
    call amrex_finalize
  
  end subroutine run_finalize
  
  
end module init_module
