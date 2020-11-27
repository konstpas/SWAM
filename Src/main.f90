program SWAM


use amrex_amr_module 
use amr_data_module 
use initdata_module 
use my_amr_module
use compute_dt_module 
use plotfile_module 
use timestep_module

implicit none

    real(amrex_real) :: cur_time
    integer :: last_plot_file_step, step, lev, substep, finest_level

call amrex_init 
call amrex_amrcore_init 

call my_amr_init() 
call initdata() 


!!!.................................................

    !use my_amr_module, only : stepno, max_step, stop_time, dt, plot_int
    !use amr_data_module, only : phi_old, phi_new, t_new
    !use compute_dt_module, only : compute_dt
    !use plotfile_module, only : writeplotfile
    !real(amrex_real) :: cur_time
    !integer :: last_plot_file_step, step, lev, substep, finest_level

    cur_time = t_new(0)   !! current time 
    last_plot_file_step = 0;
    
    do step = stepno(0), max_step-1  !! main time loop 
       if (cur_time .ge. stop_time) exit

       if (amrex_parallel_ioprocessor()) then
          print *, ""
          print *, "STEP", step+1, "starts ..."
       end if

       call compute_dt()   !! calculate time step size 

       lev = 0
       substep = 1
       call timestep(lev, cur_time, substep)
 
       cur_time = cur_time + dt(0)

       if (amrex_parallel_ioprocessor()) then
          print *, "STEP", step+1, "end. TIME =", cur_time, "DT =", dt(0)
       end if

       ! sync up time to avoid roundoff errors
       finest_level = amrex_get_finest_level()
       do lev = 0, finest_level
          t_new(lev) = cur_time
       end do

       if (plot_int .gt. 0 .and. mod(step+1,plot_int) .eq. 0) then
          last_plot_file_step = step+1
          call writeplotfile()
       end if

       if (cur_time .ge. stop_time - 1.e-6_amrex_real*dt(0)) exit
    end do

    if (plot_int .gt. 0 .and. stepno(0) .gt. last_plot_file_step) then
       call writeplotfile()
    end if


!!!!-------------------------------------------



call amr_data_finalize 


call amrex_amrcore_finalize 
call amrex_finalize 





end program SWAM         
