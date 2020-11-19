program SWAM

use mod_initdata
use mod_Readinput
use amrex_base_module 
use amrex_amr_module 

implicit none


print *, 'Compile test'
call initdata 
call Readinput
call amrex_init 
call amrex_amrcore_init 



call amrex_amrcore_finalize 
call amrex_finalize 





end program SWAM         
