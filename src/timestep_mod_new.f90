module timestep_module_new

  use amrex_amr_module

  implicit none
  private

  public :: timestep 

contains


  recursive subroutine timestep (lev, time, substep)

    use my_amr_module, only : regrid_int, stepno, nsubsteps, dt, do_reflux, verbose
    use amr_data_module, only : t_old, t_new, phi_old, phi_new, temp, flux_reg  
    use averagedown_module, only : averagedownto

    integer, intent(in) :: lev, substep
    real(amrex_real), intent(in) :: time    
    integer, allocatable, save :: last_regrid_step(:)
    integer :: k, old_finest_level, finest_level, fine_substep    
    
    ! Regridding 
    if (regrid_int .gt. 0) then
       
       if (.not.allocated(last_regrid_step)) then

          allocate(last_regrid_step(0:amrex_max_level))
          last_regrid_step = 0
          
       end if

      
       if (lev .lt. amrex_max_level .and. stepno(lev) .gt. last_regrid_step(lev)) then

          if (mod(stepno(lev), regrid_int) .eq. 0) then

             old_finest_level = amrex_get_finest_level()
             call amrex_regrid(lev, time)
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
    
    
    ! Advance solution
    stepno(lev) = stepno(lev)+1
    ! We need to update t_old(lev) and t_new(lev) before advance is called because of fillpatch.
    ! Fillpatch fills multifab for incrementing solution, using data from given level and level below, and interpolates ghost points 
    t_old(lev) = time
    t_new(lev) = time + dt(lev)
    call amrex_multifab_swap(phi_old(lev), phi_new(lev))
    call advance(lev, time, dt(lev))


    ! Propagate solution and synchronize levels
    if (lev .lt. amrex_get_finest_level()) then
       
       do fine_substep = 1, nsubsteps(lev+1)
          call timestep(lev+1, time+(fine_substep-1)*dt(lev+1), fine_substep) 
       end do

       ! Update coarse solution at coarse-fine interface via dPhidt = -div(+F) where F is stored flux (flux_reg)
       if (do_reflux) then
          !call flux_reg(lev+1)%reflux(phi_new(lev), 1.0_amrex_real) 
       end if

       ! A problem occurs when mesh size in region containing free interface is changed 
       ! Finer mesh resolves interface somewhere withing coarse grid point
       ! so that surface temperature is interpolated into vacuum region and diffuses on the wrong side of free interface boundary 
       call averagedownto(lev)  ! set covered coarse cells to be the average of fine
        
    end if
  
  end subroutine timestep


  subroutine advance(lev, time, dt)

    use my_amr_module, only : do_reflux, verbose
    use amr_data_module, only : phi_new, temp, surf_ind, flux_reg  
    use fillpatch_module, only : fillpatch
    use domain_module, only : get_melt_pos, reset_melt_pos 
    use shallow_water_module, only : increment_SW

    integer, intent(in) :: lev
    real(amrex_real), intent(in) :: time, dt        
    integer, parameter :: ngrow = 1    			! number of ghost points in each spatial direction 
    integer :: ncomp, idim, i,j 	! components, ranges, and indexes
    logical :: nodal(3)   					! logical for flux multifabs 
    type(amrex_multifab) :: phiborder, tempborder 		! multifabs on mfi owned tilebox, with ghost points 
    type(amrex_mfiter) :: mfi					! mfi iterator 
    type(amrex_box) :: bx, tbx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pin, pout, ptempin, ptemp ! input, output pointers 
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
    call amrex_multifab_build(tempborder, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, ngrow)
    call fillpatch(lev, time, phiborder)
    
    ! Propagate SW equations (only at max level)
    if (lev.eq.amrex_max_level) then 
       call increment_SW(time, amrex_geom(lev), dt)
    end if

    ! RE-distribute energy from 'lost and gained' domain points how


    ! Set melt interface position array equal to free interface position array 
    ! Since melt layer may span several tile boxes in y-direction (in mfiterator below), we cannot reset within each loop 
    ! Therefore we reset melt position after solving SW, and before propagating temperature 
    ! Melt position is then found after heat has been propagated 
    call reset_melt_pos() 
   
    

!$omp parallel private(mfi,bx,tbx,pin,pout,ptemp)
    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.true.) ! Tiling splits validbox into several tile boxes 
    								! could be useful depending on parallelization approach 
    do while(mfi%next())
    bx = mfi%tilebox()   
   
   
       pin     => phiborder%dataptr(mfi)
       pout    => phi_new(lev)%dataptr(mfi)
       ptempin => tempborder%dataptr(mfi)
       ptemp   => temp(lev)%dataptr(mfi)
       
       ! Increment solution at given mfi tilebox 
       call increment(time, bx%lo, bx%hi, &
                      pin,     lbound(pin),     ubound(pin),     &
                      pout,    lbound(pout),    ubound(pout),    &
                      ptempin, lbound(ptempin), ubound(ptempin), &
                      ptemp,   lbound(ptemp),   ubound(ptemp),   &
                      amrex_geom(lev), dt)  
		

	! Find melt interface y position 
       if (lev.eq.amrex_max_level) then
          call get_melt_pos(bx%lo, bx%hi,                      	&
                          ptemp, lbound(ptemp), ubound(ptemp), 	&
			  amrex_geom(lev))
	end if  

    end do

    call amrex_mfiter_destroy(mfi)
!$omp end parallel

    call amrex_multifab_destroy(phiborder)
    
 end subroutine advance

  !----------------------------------------------------------------------------------
  ! All this subroutine should be moved to a dedicated module for the heat conduction

  ! Subroutine takes input enthalpy and indexes (with ghost points interpolated from nearby multifab at level and below) 
  ! Updates enthalpy by divergence of enthalpy flux 

  subroutine increment(time, lo, hi, &
  			uin,    ui_lo, ui_hi, &
  			uout,   uo_lo, uo_hi, & 
  			tempin, ti_lo, ti_hi, & 
  			temp,   t_lo , t_hi , & 
  			geom, dt)
  			
    use domain_module 	
    use material_properties_module, only : get_temp, get_maxdiffus  		
    integer, intent(in) :: lo(3), hi(3)  		! bounds of current tile box
    real(amrex_real), intent(in) :: dt, time		! sub time step, and time 
    type(amrex_geometry), intent(in) :: geom  	! geometry at level
    integer, intent(in) :: ui_lo(3), ui_hi(3)		! bounds of input enthalpy box 
    integer, intent(in) :: uo_lo(3), uo_hi(3)		! bounds of output enthalpy box  
    integer, intent(in) :: ti_lo(3), ti_hi(3)		! bounds of temperature box  
    integer, intent(in) :: t_lo (3), t_hi (3)		! bounds of temperature box    
    real(amrex_real) :: dx(3) 			! Grid size 

    real(amrex_real), intent(in   ) :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! Enthalpy in 
    real(amrex_real), intent(inout) :: uout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3)) ! Enthalpy out 
    real(amrex_real), intent(inout) :: tempin(ti_lo(1):ti_hi(1),ti_lo(2):ti_hi(2),ti_lo(3):ti_hi(3)) ! Temperature on enthalpy in-box
    real(amrex_real), intent(inout) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) ! Temperature out 

    
    ! real(amrex_real) :: uface (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! face velocity x direction (nodal)
    ! real(amrex_real) :: fluxx (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! flux x direction (nodal)
    ! real(amrex_real) :: fluxy (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! flux y direction (nodal)
    
    ! real(amrex_real) :: qbound(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))	! Volumetric heating (boundary)
    ! real(amrex_real) :: qheat (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))	! Volumetric heating
    
    ! logical :: xfluxflag(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! surface flag for x-nodes 
    ! logical :: yfluxflag(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! surface flag for y-nodes   
    ! integer :: i,j,k  
    ! real(amrex_real) :: wface (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! face velocity z direction (nodal)	 
    ! real(amrex_real) :: fluxz (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! flux z direction (nodal)
    ! logical :: zfluxflag(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! surface flag for z-nodes

   
    ! dx = geom%dx(1:3) ! grid width at level 
    
    

    !     call get_temp(ti_lo, ti_hi,             &	! tilebox indexes 
    !     	      ui_lo, ui_hi, uin, &	! Output enthalpy indexes and data array
    !     	      ti_lo, ti_hi, tempin)	! Temperature indexes and data array    
    
    
  
    !     ! Subroutine assigns logical arrays denoting free interface boundary 
    !     call surface_tag(time, geom%get_physical_location(ui_lo), dx, lo, hi, &
    !     		xfluxflag, yfluxflag, ui_lo, ui_hi &
    !     		, zfluxflag & 
    !     		)	
  	
    !     ! Subroutine assigns tangential (x,z) velocity on grid edges in whole domain  
    !     call get_face_velocity(time, geom%get_physical_location(lo), dx, lo, hi, uface, &
    !     			wface, &
    !     			ui_lo, ui_hi ) 	
 
 		
  				
  				 	

    !     ! Subroutine assigns enthalpy flux on grid edges in whole domain 
    !     call create_face_flux(time, geom%get_physical_location(ui_lo), dx, lo, hi, & 			! domain module 
    !     			uin, uface, xfluxflag, yfluxflag, 	& 
    !     			fluxx, fluxy,  			& 
    !     			wface, zfluxflag, fluxz, 		&
    !     			        ui_lo, ui_hi, 			&
    !     			tempin, ti_lo, ti_hi)
  	
  	
  	
    !     ! Zero flux across surface boundary. Volumetric heat deposition in first internal cell constitutes absorbed boundary flux. 
    !     ! Incorporates all absorption and cooling terms 			
    !     call get_bound_heat(time, geom%get_physical_location(lo), dx, lo, hi, ui_lo, ui_hi, yfluxflag, qbound) 	! domain module 			
    !     !call volume_heating(time, geom%get_physical_location(ui_lo), dx, qheat) 
  	


  	
    !     ! Increment enthalpy 
    !     do   i = lo(1),hi(1)
    !      do  j = lo(2),hi(2) 
    !       do k = lo(3),hi(3)
    !       uout(i,j,k) = uin(i,j,k) &
    !          - dt/dx(1)      * (fluxx(i+1,j  ,k  )-fluxx(i,j,k))	&		! flux divergence x-direction 
    !          - dt/dx(2)      * (fluxy(i  ,j+1,k  )-fluxy(i,j,k))	& 		! flux divergence y-direction 
    !          - dt/dx(3)      * (fluxz(i  ,j  ,k+1)-fluxz(i,j,k))	&		! flux divergence z-direction
    !          + dt*qbound(i,j,k)	
    !       end do    						! 'boundary volumetric' source
    !      end do 
    !     end do 
  	
	
    !     ! Find temperature 
    !     call get_temp(lo, hi,             &	! tilebox indexes 
    !     	      uo_lo, uo_hi, uout, &	! Output enthalpy indexes and data array
    !     	      t_lo , t_hi , temp)	! Temperature indexes and data array 
  		      
    !     ! find maximum diffusivity for time step determination 
    !     ! Not called noe, constant max possible diffusivity used.	      
    !     !call get_maxdiffus(lo, hi, & 
    !     !		    t_lo, t_hi, temp)
  	
  	
  
 end subroutine increment

  ! !----------------------------------------------------------------------------------
  ! ! All this subroutine should be moved to a dedicated module for the heat conduction

end module timestep_module_new
