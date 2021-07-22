module my_amr_module

  use iso_c_binding
  use amrex_amr_module
  use amr_data_module
  
  implicit none
  
contains

  subroutine my_amr_init ()
    
    use bc_module, only : lo_bc, hi_bc
    use read_input_module
    use material_properties_module
    use regrid_module
    
    call amrex_init 
    call amrex_amrcore_init
    
    call amrex_init_virtual_functions (my_make_new_level_from_scratch, &
                                       my_make_new_level_from_coarse,  &
                                       my_remake_level,                &
                                       my_clear_level,                 &
                                       my_error_estimate)

    call read_input_file

    call init_mat_prop
	
    call amr_data_init
    
  end subroutine my_amr_init
  
end module my_amr_module
