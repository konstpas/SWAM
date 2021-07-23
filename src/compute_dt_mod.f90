module compute_dt_module

  ! -----------------------------------------------------------------
  ! This module is used to compute the timestep employed in the
  ! solution of both the heat conduction equation and the shallow-
  ! water equation
  ! -----------------------------------------------------------------
  
  use amrex_amr_module

  implicit none
  private

  ! ------------------------------------------------------------------
  ! Public subroutines
  ! ------------------------------------------------------------------
  public :: compute_dt

contains

  subroutine compute_dt()
    
    use amr_data_module, only : t_new, dt, nsubsteps
    use read_input_module, only : stop_time
    
    integer :: lev, nlevs, n_factor
    real(amrex_real) :: dt_0, eps
    real(amrex_real), allocatable :: dt_tmp(:)
    real(amrex_real), parameter :: change_max = 1.1_amrex_real ! This should be moved to input

    nlevs = amrex_get_numlevels()
    
    allocate(dt_tmp(0:nlevs-1))

    ! Estimate timestep in order to avoid instabilities
    do lev = 0, nlevs-1
       call est_timestep(lev, dt_tmp(lev))
    end do
    call amrex_parallel_reduce_min(dt_tmp, nlevs)

    ! Ensure that the timestep does not change too much
    ! from one timestep to the next
    dt_0 = dt_tmp(0)
    n_factor = 1
    do lev = 0, nlevs-1
       dt_tmp(lev) = min(dt_tmp(lev), change_max*dt(lev))
       n_factor = n_factor * nsubsteps(lev)
       dt_0 = min(dt_0, n_factor*dt_tmp(lev))
    end do

    ! Compute timestep for all levels
    dt(0) = dt_0
    do lev = 1, nlevs-1
       dt(lev) = dt(lev-1) / nsubsteps(lev)
    end do

    
  end subroutine compute_dt


  subroutine est_timestep(lev, dt)

    use read_input_module, only : cfl
    use material_properties_module, only : max_diffus 

    ! Input and output variables 
    integer, intent(in) :: lev
    real(amrex_real), intent(out) :: dt

    ! Local variables 
    real(amrex_real) :: dxsqr

    ! NOTE: There are two stability criteria to take into
    ! account, the von Neumann stability and the CFL
    ! condition. The CFL is not yet implemented
    
    ! Von Neumann stability criterion 
    dxsqr= (1/amrex_geom(lev)%dx(1)**2 + &
            1/amrex_geom(lev)%dx(2)**2 + &
            1/amrex_geom(lev)%dx(3)**2) 
    dt = 0.5/(dxsqr*max_diffus)
    dt = dt * cfl

  end subroutine est_timestep

end module compute_dt_module
