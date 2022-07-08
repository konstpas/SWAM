module plotfile_module

  ! -----------------------------------------------------------------
  ! This module is used to write the output files
  ! -----------------------------------------------------------------
  
  implicit none

  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: writeplotfile
    
contains


  ! -----------------------------------------------------------------
  ! Subroutine used to write plotfiles (containing the multifabs)
  ! and regular text files (with free surface coordinates and
  ! similar)
  ! -----------------------------------------------------------------
  subroutine writeplotfile()

    use amrex_amr_module
    use amr_data_module, only : phi_new, &
                                temp, &
                                idomain, &
                                psi, &
                                t_new, &
                                surf_pos, &
                                surf_temperature, &
                                surf_deformation, &
                                surf_current, &
                                melt_pos, &
                                melt_top, &
                                surf_ind, &
                                surf_dx, &
                                stepno, &
                                surf_accum_current_x, & 
                                surf_accum_current_y
    
    use read_input_module, only : io_plot_file
    
    integer :: i, ilev
    integer :: nlevs
    character(len=127) :: name
    character(len=16)  :: current_step
    real(amrex_real) :: xpos, max_temp, max_temp_lev
    type(amrex_string) :: varname(1)    

    ! Time step output
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

    ! Get number of levels
    nlevs = amrex_get_numlevels()

    ! Output enthalpy
    name = trim(io_plot_file) // "_enthalpy_" //current_step 
    call amrex_string_build(varname(1), "phi")
    call amrex_write_plotfile(name, nlevs, phi_new, &
                              varname, amrex_geom, &
                              t_new(0), stepno, &
                              amrex_ref_ratio)

    ! Output temperature
    name = trim(io_plot_file) // "_temperature_" //current_step 
    call amrex_string_build(varname(1), "Temperature")
    call amrex_write_plotfile(name, nlevs, temp, &
                              varname, amrex_geom, &
                              t_new(0), stepno, &
                              amrex_ref_ratio)

    ! Output flag to distinguish between material and background
    name = trim(io_plot_file) // "_idomain_" //current_step 
    call amrex_string_build(varname(1), "idomain")
    call amrex_write_plotfile(name, nlevs, idomain, &
                              varname, amrex_geom, &
                              t_new(0), stepno, &
                              amrex_ref_ratio)

    ! Output electrostatic potential
    name = trim(io_plot_file) // "_potential_" //current_step 
    call amrex_string_build(varname(1), "psi")
    call amrex_write_plotfile(name, nlevs, psi, &
                              varname, amrex_geom, &
                              t_new(0), stepno, &
                              amrex_ref_ratio)

    max_temp = 0.0
    max_temp_lev = 0.0
    do ilev = 0, nlevs-1
      max_temp_lev = temp(ilev)%max(1,0)
      if (max_temp_lev.gt.max_temp) max_temp = max_temp_lev
    end do

    ! Output melt thickness
    name = "melt_thickness_" //trim(current_step)//".dat"
    open(2, file = name, status = 'unknown', action = "write")
    write(2, *) 'x-coordinate     Free Surface     Melt Bottom     Melt top     Free Surface Temperature     Surface deformation', &
            '    Jx    Jy    Jth    Max temperature'
    do i=surf_ind(1,1), surf_ind(1,2)
      ! i starts from 0 so to output the x-coord at the center of the cell add 0.5
      xpos = (i+0.5)*surf_dx(1)
      if (i.eq.surf_ind(1,1)) then
         write(2, '(10(es13.6, 4x))') xpos, surf_pos(i), melt_pos(i), melt_top(i), surf_temperature(i), surf_deformation(i), &
         surf_accum_current_x(i,1)/surf_accum_current_x(i,2), surf_accum_current_y(i,1)/surf_accum_current_y(i,2), & 
         surf_current(i), max_temp
      else
         write(2, '(9(es13.6, 4x))') xpos, surf_pos(i), melt_pos(i), melt_top(i), surf_temperature(i), surf_deformation(i), &
         surf_accum_current_x(i,1)/surf_accum_current_x(i,2), surf_accum_current_y(i,1)/surf_accum_current_y(i,2), surf_current(i)
      endif
    end do
    close(2)
  end subroutine writeplotfile

 
end module plotfile_module

