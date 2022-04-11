module init_module

  ! -----------------------------------------------------------------
  ! This module is used to initialize and finalize the simulation
  ! -----------------------------------------------------------------
  
  use iso_c_binding
  use amrex_amr_module
  
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

    use amr_data_module, only : amr_data_init
    use read_input_module, only : cooling_debug, &
                                  read_input_file, &
                                  plasma_flux_type, &
                                  restart, &
                                  geom_name, &
                                  plasma_flux_side_type, &
                                  solve_heat 
    use heat_transfer_read_flux_module, only : construct_plasma_flux_table
    use material_properties_module, only : init_mat_prop
    use heat_transfer_flux_module, only : debug_cooling_fluxes
    use shallow_water_module, only : init_melt_pos
    use regrid_module, only : averagedown, &
                              my_make_new_level_from_scratch, &
                              my_make_new_level_from_coarse, &
                              my_remake_level, &
                              my_clear_level, &
                              my_error_estimate

    ! Initialize amrex
    call amrex_init 
    call amrex_amrcore_init
    
    ! Read input file
    call read_input_file

    ! Initialize tables with material properties
    call init_mat_prop

    ! Print cooling fluxes to table if debug is on
    if (nint(cooling_debug(1)) .eq. 1) call debug_cooling_fluxes

    ! Read heat flux from file
    if (plasma_flux_type.eq.'Input_file') then
       call construct_plasma_flux_table(.false.)
    end if
    if (geom_name.eq.'West' .and. &
        plasma_flux_side_type.eq.'Input_file') then
       call construct_plasma_flux_table(.true.)
    end if
    
    ! Initialize amrex data used in the simulation
    call amr_data_init

    ! Initialize melt pool if thermal response is not simulated
    if (.not.solve_heat) call init_melt_pos

    ! Initialize amrex functions to generate grid
    call amrex_init_virtual_functions(my_make_new_level_from_scratch, &
                                      my_make_new_level_from_coarse, &
                                      my_remake_level, &
                                      my_clear_level, &
                                      my_error_estimate)
    
    ! Start new simulation or read restart file
    if (len_trim(restart) .eq. 0) then
       call amrex_init_from_scratch(0.0_amrex_real)
       call averagedown      
    else
       call amrex_abort("init from checkpoint not implemented yet")
    end if
    
  end subroutine run_init

  ! -----------------------------------------------------------------
  ! Subroutine to finalize a simulation and clean the memory
  ! -----------------------------------------------------------------
  subroutine run_finalize()
    
    use amr_data_module, only : amr_data_finalize
    use material_properties_module, only : finalize_mat_prop
    
    ! Free amrex data
    call amr_data_finalize

    ! Free tables with material properties
    call finalize_mat_prop
    
    ! Finalize amrex
    call amrex_amrcore_finalize 
    call amrex_finalize
  
  end subroutine run_finalize


  
  
end module init_module
