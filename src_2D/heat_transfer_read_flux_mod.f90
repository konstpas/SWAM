module heat_transfer_read_flux_module
  
   ! -----------------------------------------------------------------
   ! This module is used to read a file that descibes the heat flux 
   ! imposed on the free surface
   ! -----------------------------------------------------------------
   
  use amrex_amr_module
  
  implicit none
    
  private

  public :: construct_plasma_flux_table

contains

  ! -----------------------------------------------------------------
  ! Subroutine used to construct the table containing the heat flux
  ! read from file
  ! -----------------------------------------------------------------
  subroutine construct_plasma_flux_table(side_flag)
    
    use amr_data_module, only: plasma_flux_time_mesh, &
                               plasma_flux_surf_mesh, &
                               plasma_flux_table, &
                               plasma_flux_side_time_mesh, &
                               plasma_flux_side_surf_mesh, &
                               plasma_flux_side_table
                               
    use read_input_module, only: heat_plasma_flux_file, &
                                 heat_plasma_flux_side_file
    
    use read_files_module, only: get_mesh_dimensions, &
                                 read_heatflux_file
    ! Input and output variables
    ! side_flag = true --> construct the heat-flux matrix for the exposed side (only in WEST geometry).
    ! side_flag = false --> construct the matrix for the heat-flux to the free surface
    logical, intent(in) :: side_flag
    
    ! Local variables
    integer :: dims(1:2)
    
    if(side_flag) then
       
      call get_mesh_dimensions (heat_plasma_flux_side_file, 2, len(heat_plasma_flux_side_file), dims)
      
      allocate (plasma_flux_side_time_mesh(1:dims(1)) )
      allocate (plasma_flux_side_surf_mesh(1:dims(2)) )
      allocate (plasma_flux_side_table(1:dims(1),1:dims(2)) )
      
      call read_heatflux_file(heat_plasma_flux_side_file, plasma_flux_side_time_mesh, &
                              plasma_flux_side_surf_mesh, plasma_flux_side_table)
      
   else
      
      call get_mesh_dimensions (heat_plasma_flux_file, 2, len(heat_plasma_flux_file), dims)
      
      allocate (plasma_flux_time_mesh(1:dims(1)) )
      allocate (plasma_flux_surf_mesh(1:dims(2)) )
      allocate (plasma_flux_table(1:dims(1),1:dims(2)) )
      
      call read_heatflux_file(heat_plasma_flux_file, plasma_flux_time_mesh, &
                              plasma_flux_surf_mesh, plasma_flux_table)
      
    end if
    
    
  end subroutine construct_plasma_flux_table
      

end module heat_transfer_read_flux_module