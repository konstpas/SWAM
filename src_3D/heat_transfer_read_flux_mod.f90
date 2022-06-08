module read_heat_flux_module
  
   ! -----------------------------------------------------------------
   ! This module is used to read a file that descibes the heat flux 
   ! imposed on the free surface
   ! -----------------------------------------------------------------
   
    use amrex_amr_module

    implicit none

    private

    public :: construct_heat_flux_table

  contains

    ! -----------------------------------------------------------------
    ! Subroutine used to construct the table containing the heat flux
    ! read from file
    ! -----------------------------------------------------------------
    subroutine construct_heat_flux_table(side_flag)
      
      use amr_data_module, only: plasma_flux_time_mesh, &
                                 plasma_flux_surf_x_mesh, &
                                 plasma_flux_surf_z_mesh, &
                                 heat_flux_table, &
                                 plasma_side_flux_time_mesh, &
                                 plasma_side_flux_surf_y_mesh, &
                                 plasma_side_flux_surf_z_mesh, &
                                 heat_side_flux_table

      use read_input_module, only: heat_plasma_flux_file, &
                                   heat_plasma_flux_side_file

      use read_files_module, only: get_mesh_dimensions

      use read_files_module, only: get_mesh_dimensions, &
                                   read_heatflux_file
      
      implicit none

     ! Input and output variables
      logical, intent(in) :: side_flag  ! True if contructing heat-flux matrix for the exposed side (only in WEST geometry)
                                        ! false if constructing the matrix for the heat-flux to the free surface.
      ! Local variables
      integer :: dims(1:3)
      
      if(side_flag) then
         call get_mesh_dimensions (heat_plasma_flux_side_file, 3, len(heat_plasma_flux_side_file), dims)

        allocate (plasma_side_flux_time_mesh(1:dims(1)))
        allocate (plasma_side_flux_surf_y_mesh(1:dims(2)))
        allocate (plasma_side_flux_surf_z_mesh(1:dims(3)))
        allocate (heat_side_flux_table(1:dims(1),1:dims(2),1:dims(3)) )
        
        call read_heatflux_file(heat_plasma_flux_side_file, len(heat_plasma_flux_side_file), plasma_side_flux_time_mesh, &
                                plasma_side_flux_surf_y_mesh, plasma_side_flux_surf_z_mesh, heat_side_flux_table)
      else
        call get_mesh_dimensions (heat_plasma_flux_file, 3, len(heat_plasma_flux_file), dims)

        allocate (plasma_flux_time_mesh(1:dims(1)))
        allocate (plasma_flux_surf_x_mesh(1:dims(2)))
        allocate (plasma_flux_surf_z_mesh(1:dims(3)))
        allocate (heat_flux_table(1:dims(1),1:dims(2),1:dims(3)) )
     
        call read_heatflux_file(heat_plasma_flux_file, len(heat_plasma_flux_file), plasma_flux_time_mesh, &
                                plasma_flux_surf_x_mesh, plasma_flux_surf_z_mesh, heat_flux_table)
      end if
      
    end subroutine construct_heat_flux_table
    
end module read_heat_flux_module
