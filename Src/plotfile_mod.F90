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
    integer :: nlevs 
    integer :: statio = 0 
    real(amrex_real) :: melt_vol 
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


    
  end subroutine write2dplotfile
  
  
  
  
  
  
  
  
  
  
  
  
  

end module plotfile_module

