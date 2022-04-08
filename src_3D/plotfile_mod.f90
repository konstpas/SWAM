module plotfile_module

  ! -----------------------------------------------------------------
  ! This module is used to write the output files
  ! -----------------------------------------------------------------
  
  use amrex_amr_module
  use amr_data_module, only : phi_new, &
                              temp, &
                              idomain, &
                              t_new, &
                              surf_pos, &
                              surf_temperature, &
                              melt_pos, &
                              melt_top, &
                              surf_ind, &
                              surf_dx, &
                              stepno, &
                              melt_vel, &
                              qnew, &
                              Qpipe, &
                              Qtherm, &
                              Qrad, &
                              Qvap, &
                              Qplasma
  
  use read_input_module, only : plot_file, solve_heat
  
  implicit none

  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: writeplotfile
  public :: write2dplotfile 
  public :: initialize_heatfluxes_file

  ! -----------------------------------------------------------------
  ! Declare private variables shared by all subroutines
  ! -----------------------------------------------------------------
  logical, save :: init_2d = .true. 
  
    
contains


  ! -----------------------------------------------------------------
  ! Subroutine used to output the multifabs containing, for example
  ! the enthalpy and the temperature in all domain at all levels
  ! -----------------------------------------------------------------
  subroutine writeplotfile()

    integer :: i, k, ilev
    integer :: nlevs
    character(len=127) :: name
    character(len=16)  :: current_step
    character(len=16)  :: dashfmt
    real(amrex_real) :: xpos, zpos
    real(amrex_real) :: max_temp, max_temp_lev
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

    if (solve_heat) then
      ! Get number of levels
      nlevs = amrex_get_numlevels()

      ! Enthalpy output
      name = trim(plot_file) // "_enthalpy_" //current_step 
      call amrex_string_build(varname(1), "phi")
      call amrex_write_plotfile(name, nlevs, phi_new, &
                                 varname, amrex_geom, &
                                 t_new(0), stepno, &
                                 amrex_ref_ratio)

      ! Temperature output
      name = trim(plot_file) // "_temperature_" //current_step 
      call amrex_string_build(varname(1), "Temperature")
      call amrex_write_plotfile(name, nlevs, temp, &
                                 varname, amrex_geom, &
                                 t_new(0), stepno, &
                                 amrex_ref_ratio)

      ! Output flag to distinguish between material and background
      name = trim(plot_file) // "_idomain_" //current_step 
      call amrex_string_build(varname(1), "idomain")
      call amrex_write_plotfile(name, nlevs, idomain, &
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
      write(2, *) 'x-coordinate  z-coordinate   Free Surface     Melt Bottom     Melt top     Free Surface Temperature', &
          '    Max temperat'
      dashfmt = '(6(es13.6, 4x))'
      do i=surf_ind(1,1), surf_ind(1,2)
          do k=surf_ind(2,1), surf_ind(2,2)
              ! i starts from 0 so to output the x-coord at the center of the cell add 0.5
              ! the same applies for k and the z-coord
              xpos = (i+0.5)*surf_dx(1)
              zpos = (k+0.5)*surf_dx(2)
              if(i.eq.surf_ind(1,1) .and. k.eq.surf_ind(2,1)) then
                write(2, '(7(es13.6, 4x))') xpos, zpos, surf_pos(i,k), melt_pos(i,k), melt_top(i,k), surf_temperature(i,k), max_temp
              else
              write(2, dashfmt) xpos, zpos, surf_pos(i,k), melt_pos(i,k), melt_top(i,k), surf_temperature(i,k)
            end if
          end do
      end do
      close(2)
    end if

    ! Output melt thickness
    name = "heat_fluxes.dat"! //trim(current_step)//".dat"
    open(2, file = name, status = 'old', action = "write", access = 'append')
    dashfmt = '(6(es13.6, 4x))'
    write(2, dashfmt) t_new(0), Qtherm, Qrad, Qvap, Qpipe, Qplasma
    close(2)

    ! Output shallow water variables
    name = "sw_" //trim(current_step)//".dat"
    open(2, file = name, status = 'unknown', action = "write")
    write(2, *) 'x-coordinate  z-coordinate   Free Surface     Melt Bottom     Melt top     '&
                'Melt vx     Melt vz     qnew1     qnew2     qnew3'
    dashfmt = '(11(es13.6, 4x))'
    do i=surf_ind(1,1), surf_ind(1,2)
        do k=surf_ind(2,1), surf_ind(2,2)
           ! i starts from 0 so to output the x-coord at the center of the cell add 0.5
           ! the same applies for k and the z-coord
           xpos = (i+0.5)*surf_dx(1)
           zpos = (k+0.5)*surf_dx(2)
           write(2, dashfmt) xpos, zpos, surf_pos(i,k), melt_pos(i,k), melt_top(i,k), &
                             melt_vel(i,k,1) ,melt_vel(i,k,2), qnew(1,i,k), qnew(2,i,k), &
                             qnew(3,i,k), surf_temperature(i,k)
        end do
    end do
    close(2)
    
  end subroutine writeplotfile

  
  ! -----------------------------------------------------------------
  ! Subroutine used to output the total volume of the melt pool
  ! -----------------------------------------------------------------  
  subroutine write2dplotfile()
    
    use domain_module, only : integrate_surf
    
    real(amrex_real) :: melt_vol 

    ! Get total volume of molten pool
    call integrate_surf(melt_vol) 

    ! Write output header
    if (init_2d) then 

       init_2d = .false. 
       open (1, file = 'surf_1D_time.dat', &
             status = 'unknown') 
       write(1,*) 'Surface properties vs time' 
       write(1,*) 'Time [s]' , 'Integrated molten volume [mm3]'
       write(1,*) ''
       
    end if

    ! Write output
    open (1, file = 'surf_1D_time.dat', &
          status = 'old', &
          access = 'append')    
    write(1,*) t_new(0), melt_vol  
    close(1)
    
  end subroutine write2dplotfile

  subroutine initialize_heatfluxes_file()

   open (1, file = 'heat_fluxes.dat', &
         status = 'unknown') 
   write(1, *) 'Time[s]          Thermionic[J]    Radiation[J]     Vaporization[J]  Cooling Pipe[J]  Plasma[J]'
   close(1)
  end subroutine

  
  
end module plotfile_module
