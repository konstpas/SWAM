module my_amr_module

  use iso_c_binding
  use amrex_amr_module
  use amrex_fort_module, only : rt => amrex_real

  use amr_data_module
  implicit none
  
  !integer, private, parameter :: ncomp = 1, nghost = 0
  
contains

  subroutine my_amr_init()
    
    use bc_module, only : lo_bc, hi_bc
    use read_input_module, only : read_input_file
    use material_properties_module, only : init_mat_prop
    
    integer :: ilev

    call amrex_init 
    call amrex_amrcore_init 
    call amrex_init_virtual_functions(my_make_new_level_from_scratch, &
                                      my_make_new_level_from_coarse,  &
                                      my_remake_level,                &
                                      my_clear_level,                 &
                                      my_error_estimate)

    call read_input_file

    call init_mat_prop
 
    call amr_data_init
 
  end subroutine my_amr_init


  ! subroutine my_amr_finalize ()
  !   call amr_data_finalize()
  ! end subroutine my_amr_finalize

  
end module my_amr_module
