module timestep_module

  use amrex_amr_module

  implicit none
  private

  public :: timestep 

contains


  recursive subroutine timestep (lev, time, substep)
    use my_amr_module, only : regrid_int, stepno, nsubsteps, dt, do_reflux !! input
    use amr_data_module, only : t_old, t_new, phi_old, phi_new, flux_reg !! transients of the solution 
    use averagedown_module, only : averagedownto
    integer, intent(in) :: lev, substep
    real(amrex_real), intent(in) :: time
    
    integer, allocatable, save :: last_regrid_step(:)
    integer :: k, old_finest_level, finest_level, fine_substep
    integer :: redistribute_ngrow
    
    if (regrid_int .gt. 0) then
       if (.not.allocated(last_regrid_step)) then
          allocate(last_regrid_step(0:amrex_max_level))
          last_regrid_step = 0
       end if

       ! regrid doesn't change the base level, so we don't regrid on amrex_max_level
       if (lev .lt. amrex_max_level .and. stepno(lev) .gt. last_regrid_step(lev)) then
          if (mod(stepno(lev), regrid_int) .eq. 0) then

             old_finest_level = amrex_get_finest_level()
             call amrex_regrid(lev, time) ! the finest level may change during regrid
             finest_level = amrex_get_finest_level()

             do k = lev, finest_level
                last_regrid_step(k) = stepno(k)
             end do

             do k = old_finest_level+1, finest_level
                dt(k) = dt(k-1) / amrex_ref_ratio(k-1)
             end do
             
          end if
       end if
    end if

    stepno(lev) = stepno(lev)+1

    ! We need to update t_old(lev) and t_new(lev) before advance is called because of fillpath.
    t_old(lev) = time
    t_new(lev) = time + dt(lev)
    ! swap phi_new(lev) and phi_old(lev) so they are consistent with t_new(lev) and t_old(lev)
    call amrex_multifab_swap(phi_old(lev), phi_new(lev))

    call advance(lev, time, dt(lev), stepno(lev), substep, nsubsteps(lev))

    if (lev .lt. amrex_get_finest_level()) then
       do fine_substep = 1, nsubsteps(lev+1)
          call timestep(lev+1, time+(fine_substep-1)*dt(lev+1), fine_substep)
       end do

       if (do_reflux) then
          !call flux_reg(lev+1)%reflux(phi_new(lev), 1.0_amrex_real)
       end if

       call averagedownto(lev)  !! set covered coarse cells to be the average of fine 
    end if
  
  end subroutine timestep















  ! update phi_new(lev)
  subroutine advance (lev, time, dt, step, substep, nsub)
    use my_amr_module, only : verbose, do_reflux
    use amr_data_module, only : phi_new, phi_old, flux_reg
    use fillpatch_module, only : fillpatch
    integer, intent(in) :: lev, step, substep, nsub
    real(amrex_real), intent(in) :: time, dt


integer, parameter :: ngrow = 1 ! 3 ! number of ghost points in each spatial direction 
    integer :: ncomp, idim, lowbound(2), hbound(2), i,j 
    logical :: nodal(3)
    type(amrex_multifab) :: phiborder
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx, tbx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pin,pout
    real(amrex_real) :: dx1,dx2,transf
    type(amrex_multifab) :: fluxes(amrex_spacedim)    
    
    ncomp = phi_new(lev)%ncomp()

    if (do_reflux) then
       do idim = 1, amrex_spacedim
          nodal = .false.
          nodal(idim) = .true.
          !call amrex_multifab_build(fluxes(idim), phi_new(lev)%ba, phi_new(lev)%dm, ncomp, 0, nodal)
       end do
    end if

    call amrex_multifab_build(phiborder, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, ngrow)


    call fillpatch(lev, time, phiborder)

!$omp parallel private(mfi,bx,tbx,pin,pout)
    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.true.)
    do while(mfi%next())
       bx = mfi%tilebox()
       pin  => phiborder%dataptr(mfi)
       pout => phi_new(lev)%dataptr(mfi)
       

	do i=lbound(pout,1),ubound(pout,1) ! three ghost points, reaching outside of pout.. 
		do j = lbound(pout,2),ubound(pout,2)
		!print*, amrex_geom(lev)%dx
		!transf = 1.2 !pin(i,j,1,1)
		!pout(i,j,1,1) = 1.! pout(i,j,1,1)
		pout(i,j,1,1) = pin(i,j,1,1) + & !! why the extra dimensions if spacedim = 2 etc. 
		dt/(amrex_geom(lev)%dx(1)**2) * (pin(i-1,j  ,1,1)-2*pin(i,j,1,1) + pin(i+1,j  ,1,1)) + & 
		dt/(amrex_geom(lev)%dx(2)**2) * (pin(i  ,j-1,1,1)-2*pin(i,j,1,1) + pin(i  ,j+1,1,1))
		end do 
	end do 
    end do ! while(mfi%next())
    call amrex_mfiter_destroy(mfi)
!$omp end parallel

call amrex_multifab_destroy(phiborder)


  end subroutine advance

end module timestep_module
