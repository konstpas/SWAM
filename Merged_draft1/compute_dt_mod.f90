module compute_dt_module

  use amrex_amr_module

  implicit none
  private

  public :: compute_dt

contains

  subroutine compute_dt ()
    use my_amr_module, only : t_new, dt, stop_time, nsubsteps

    integer :: lev, nlevs, n_factor
    real(amrex_real) :: dt_0, eps
    real(amrex_real), allocatable :: dt_tmp(:)
    real(amrex_real), parameter :: change_max = 1.1_amrex_real

    nlevs = amrex_get_numlevels()
    
    allocate(dt_tmp(0:nlevs-1))
    do lev = 0, nlevs-1
       dt_tmp(lev) = est_timestep(lev, t_new(lev))
    end do
    call amrex_parallel_reduce_min(dt_tmp, nlevs)
 
    dt_0 = dt_tmp(0)
    n_factor = 1
    do lev = 0, nlevs-1
       dt_tmp(lev) = min(dt_tmp(lev), change_max*dt(lev))
       n_factor = n_factor * nsubsteps(lev)
       dt_0 = min(dt_0, n_factor*dt_tmp(lev))
    end do
    
    ! Limit dt's by the value of stop_time.
    eps = 1.e-3_amrex_real * dt_0
    if (t_new(0) + dt_0 .gt. stop_time - eps) then
       dt_0 = (stop_time - t_new(0))/10.  
    end if

    dt(0) = dt_0
    do lev = 1, nlevs-1
       dt(lev) = dt(lev-1) / nsubsteps(lev)
    end do
  end subroutine compute_dt


  function est_timestep (lev, time) result(dt)
    use my_amr_module, only : phi_new, temp, cfl !! input
    use amr_data_module, only : melt_vel, surf_dx 
    use material_properties_module, only : max_diffus 
    real(amrex_real) :: maxvel_x,maxvel_z, dt
    integer, intent(in) :: lev
    real(amrex_real), intent(in) :: time

    real(amrex_real) :: dt_est, dt_est_heat, dt_est_sw 
    real(amrex_real) :: max_vel_val(2) 
    real(amrex_real) :: dxsqr 
    type(amrex_box) :: bx
    type(amrex_fab) :: u
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, pointer :: p(:,:,:,:)

    dt_est_heat = huge(1._amrex_real)
    dt_est_sw = huge(1._amrex_real)

    ! Von Neumann stability analysis dictates diff*dt*(1/dx2 + 1/dy2 + 1/dz2) .le. 0.5 
    dxsqr= (1/amrex_geom(lev)%dx(1)**2 + 1/amrex_geom(lev)%dx(2)**2 + 1/amrex_geom(lev)%dx(3)**2) 
    dt_est_heat = 0.5/(dxsqr*max_diffus)
    
    max_vel_val = MAXVAL(MAXVAL(abs(melt_vel),1),1) ! Maximum velocity in each direction (magnitude)
    dt_est_sw = 0.5/( max_vel_val(1)/surf_dx(1) + max_vel_val(2)/surf_dx(2) )
        
        
    dt_est = MIN(dt_est_heat,dt_est_sw) 
    dt = dt_est * cfl

  end function est_timestep

end module compute_dt_module
