program SWAM


use amrex_amr_module 
use amr_data_module 
use initdata_module 
use my_amr_module
use compute_dt_module 
use plotfile_module 
use timestep_module
use energy_module 
implicit none

    real(amrex_real) :: cur_time
    real(amrex_real) :: total_enthalpy     
    integer :: last_plot_file_step, step, lev, substep, finest_level

call amrex_init 
call amrex_amrcore_init 

call my_amr_init() ! Read input file and the initialize multifabs (subr. amr_data_init in amr data mod)
call initdata()    !Initialize multifab data from scratch (initdata mod) 


call run_simulation

call amr_data_finalize 


call amrex_amrcore_finalize 
call amrex_finalize 





end program SWAM         
