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
    use my_amr_module, only : phi_new, cfl !! input
    
    real(amrex_real) :: dt
    integer, intent(in) :: lev
    real(amrex_real), intent(in) :: time

    real(amrex_real) :: dt_est
    type(amrex_box) :: bx
    type(amrex_fab) :: u
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, pointer :: p(:,:,:,:)

    dt_est = huge(1._amrex_real)


    dt_est = min(dt_est,amrex_geom(lev)%dx(1)*amrex_geom(lev)%dx(1))
    dt_est = min(dt_est,amrex_geom(lev)%dx(2)*amrex_geom(lev)%dx(2))
    
    dt = dt_est * cfl

  end function est_timestep

end module compute_dt_module
