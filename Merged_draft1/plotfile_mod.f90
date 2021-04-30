module plotfile_module

  use amrex_amr_module
!!  use amr_data_module, only : pc
  use my_amr_module, only : plot_file, phi_new, temp, t_new, stepno

  implicit none

  private

  public :: writeplotfile, write2dplotfile 

  logical :: first_2d = .true. 
  
  
  
contains

  subroutine writeplotfile ()

    integer :: nlevs
    character(len=127) :: name
    character(len=16)  :: current_step
    type(amrex_string) :: varname(1)
    
    if      (stepno(0) .lt. 1000000) then
       write(current_step,fmt='(i5.5)') stepno(0)
    else if (stepno(0) .lt. 10000000) then
       write(current_step,fmt='(i6.6)') stepno(0)
    else if (stepno(0) .lt. 100000000) then
       write(current_step,fmt='(i7.7)') stepno(0)
    else if (stepno(0) .lt. 1000000000) then
       write(current_step,fmt='(i8.8)') stepno(0)
    else
       write(current_step,fmt='(i15.15)') stepno(0)
    end if
    name = trim(plot_file) // current_step

    nlevs = amrex_get_numlevels()

    !call amrex_string_build(varname(1), "phi")
    call amrex_string_build(varname(1), "Temperature")
    !call amrex_write_plotfile(name, nlevs, phi_new, varname, amrex_geom, &
    !     t_new(0), stepno, amrex_ref_ratio)
    call amrex_write_plotfile(name, nlevs, temp, varname, amrex_geom, &
         t_new(0), stepno, amrex_ref_ratio)

    
  end subroutine writeplotfile
  
  
  subroutine write2dplotfile ()
  use domain_module, only : integrate_surf 
  use amr_data_module, only: surf_pos, melt_vel, height_flux, surf_ind, surf_xlo, surf_dx 
    integer :: nlevs 
    integer :: statio = 0 
    integer :: i,k
    real(amrex_real) :: melt_vol 
    real(amrex_real) :: xpos, zpos, cell_meltvel_x, cell_meltvel_z, h_flux_x, h_flux_z, h_flux 
    character(len=127) :: name, name2dplt 
    character(len=16)  :: current_step
    type(amrex_string) :: varname(1)
    
    if      (stepno(0) .lt. 1000000) then
       write(current_step,fmt='(i5.5)') stepno(0)
    else if (stepno(0) .lt. 10000000) then
       write(current_step,fmt='(i6.6)') stepno(0)
    else if (stepno(0) .lt. 100000000) then
       write(current_step,fmt='(i7.7)') stepno(0)
    else if (stepno(0) .lt. 1000000000) then
       write(current_step,fmt='(i8.8)') stepno(0)
    else
       write(current_step,fmt='(i15.15)') stepno(0)
    end if
    name = trim(plot_file) // current_step

    nlevs = amrex_get_numlevels() 
    
    call integrate_surf(melt_vol) 
    
    if (first_2d) then 
    	first_2d = .false. 
  	open (1, file = 'surf_1D_time.dat', status = 'unknown') 
	
	write(1,*) 'Surface properties vs time' 
	write(1,*) 'Time [s]' , 'Integrated molten volume [mm3]'
	write(1,*) ''
	write(1,*) t_new(0), melt_vol  
	
  	close(1) 
    else 
  	open (1, file = 'surf_1D_time.dat', status = 'old', access = 'append') 
	
	write(1,*) t_new(0), melt_vol  
	
  	close(1) 
    end if 





 	! 2D surface plot file 
	name2dplt = 'surface_2D_'//trim(name)//'.csv'
	open (2, file = name2dplt, status = 'unknown') 
	write(2,*) 'X (m), ','Z (m), ' , 'Surf_Pos (m), ', 'Vel. x (m/s), ', 'Vel. z (m/s),', 'h_flux_x (m2/s),','h_flux_z (m2/s),','h_flux (m2/s)' 
		do i = surf_ind(1,1), surf_ind(1,2)
		do k = surf_ind(2,1), surf_ind(2,2)
		xpos = surf_xlo(1) + (i+0.5)*surf_dx(1)
		zpos = surf_xlo(2) + (k+0.5)*surf_dx(2) 
		cell_meltvel_x = 0.5*(melt_vel(i,k,1) + melt_vel(i+1,k,1)) 
		cell_meltvel_z = 0.5*(melt_vel(i,k,2) + melt_vel(i,k+1,2)) 
		h_flux_x = 0.5*(height_flux(i,k,1) + height_flux(i+1,k,1))
		h_flux_z = 0.5*(height_flux(i,k,2) + height_flux(i,k+1,2))
		h_flux = (h_flux_x**2 + h_flux_z**2)**0.5
		
	write(2,*) xpos,',', zpos,',', surf_pos(i,k),',', cell_meltvel_x,',', cell_meltvel_z,',',h_flux_x,',',h_flux_z,',',h_flux    
		end do 
		end do 
  	close(2) 


	!pause 
    
  end subroutine write2dplotfile
  
  
  
  
  
  
  
  
  
  
  
  
  

end module plotfile_module

