module heat_transfer_module

  use amrex_amr_module
  
  implicit none

  private

  public :: increment_enthalpy

contains

  subroutine increment_enthalpy(time, lo, hi, &
                                uin,    ui_lo, ui_hi, &
                                uout,   uo_lo, uo_hi, & 
  			        tempin, ti_lo, ti_hi, & 
  			        temp,   t_lo , t_hi , &
                                flxx, fx_lo, fx_hi, &
                                flxy, fy_lo, fy_hi, &
                                flxz, fz_lo, fz_hi, &
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
    integer         , intent(in   ) :: fx_lo(3), fx_hi(3)				! bounds of input tilebox
    integer         , intent(in   ) :: fy_lo(3), fy_hi(3)				! bounds of input tilebox
    integer         , intent(in   ) :: fz_lo(3), fz_hi(3)				! bounds of input tilebox

    
    real(amrex_real) :: dx(3) 			! Grid size 

    real(amrex_real), intent(in   ) :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! Enthalpy in 
    real(amrex_real), intent(inout) :: uout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3)) ! Enthalpy out 
    real(amrex_real), intent(inout) :: tempin(ti_lo(1):ti_hi(1),ti_lo(2):ti_hi(2),ti_lo(3):ti_hi(3)) ! Temperature on enthalpy in-box
    real(amrex_real), intent(inout) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) ! Temperature out 
    real(amrex_real), intent(  out) :: flxx    (fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3)) ! flux x direction  			
    real(amrex_real), intent(  out) :: flxy    (fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3)) ! flux y direction
    real(amrex_real), intent(  out) :: flxz    (fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3)) ! flux z direction	

    real(amrex_real) :: qbound(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))	! Volumetric heating (boundary)
    real(amrex_real) :: qheat (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))	! Volumetric heating
    
    integer :: i,j,k  

    
    ! Grid width 
    dx = geom%dx(1:3) ! grid width at level 

    
    ! Get temperature
    call get_temp(ti_lo, ti_hi,             &	! tilebox indexes 
        	  ui_lo, ui_hi, uin, &	! Output enthalpy indexes and data array
                  ti_lo, ti_hi, tempin)	! Temperature indexes and data array    
    
  
    ! Subroutine assigns logical arrays denoting free interface boundary 
    ! call surface_tag(time, geom%get_physical_location(ui_lo), dx, lo, hi, &
    !                  xfluxflag, yfluxflag, ui_lo, ui_hi, zfluxflag)	
  	
    ! Subroutine assigns tangential (x,z) velocity on grid edges in whole domain  
    ! call get_face_velocity(time, geom%get_physical_location(lo), dx, lo, hi, &
    !                        uface, wface, ui_lo, ui_hi ) 	
 			 	

    ! Subroutine assigns enthalpy flux on grid edges in whole domain 
    call create_face_flux(time, geom%get_physical_location(ui_lo), dx, lo, hi, 	        &
                          uin, ui_lo, ui_hi,              &
                          flxx, fx_lo, fx_hi,             &
                          flxy, fy_lo, fy_hi,             &
                          flxz, fz_lo, fz_hi,             &
                          tempin, ti_lo, ti_hi)
  				  	
    ! Zero flux across surface boundary. Volumetric heat deposition in first internal cell constitutes absorbed boundary flux. 
    ! Incorporates all absorption and cooling terms 			
    call get_bound_heat(time, geom%get_physical_location(lo), dx, lo, hi, ui_lo, ui_hi, 1, qbound) 	! domain module 			
    !call volume_heating(time, geom%get_physical_location(ui_lo), dx, qheat) 
  	


  	
    ! Increment enthalpy 
    do   i = lo(1),hi(1)
       do  j = lo(2),hi(2) 
          do k = lo(3),hi(3)
             uout(i,j,k) = uin(i,j,k) &
                  - dt/dx(1)      * (flxx(i+1,j  ,k  )-flxx(i,j,k))	&		! flux divergence x-direction 
                  - dt/dx(2)      * (flxy(i  ,j+1,k  )-flxy(i,j,k))	& 		! flux divergence y-direction 
                  - dt/dx(3)      * (flxz(i  ,j  ,k+1)-flxz(i,j,k))	&		! flux divergence z-direction
                  + dt*qbound(i,j,k)	                                                ! 'boundary volumetric' source
          end do
       end do
    end do
  	
    ! Scale the fluxes for the flux registers
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1) + 1
             flxx(i,j,k) = flxx(i,j,k) * (dt * dx(2)*dx(3))
          end do
       end do
    end do

    do k = lo(3), hi(3)
       do j = lo(2), hi(2) + 1
          do i = lo(1), hi(1)
             flxy(i,j,k) = flxy(i,j,k) * (dt * dx(1)*dx(3))
          end do
       end do
    end do
    
    do k = lo(3), hi(3) + 1
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             flxz(i,j,k) = flxz(i,j,k) * (dt * dx(1)*dx(2))
          end do
       end do
    end do
    
    ! Find temperature 
    call get_temp(lo, hi,             &
                  uo_lo, uo_hi, uout, &
                  t_lo , t_hi , temp) 

    ! find maximum diffusivity for time step determination 
    ! Not called noe, constant max possible diffusivity used.	      
    !call get_maxdiffus(lo, hi, & 
    !		    t_lo, t_hi, temp)
  	
  	
    
  end subroutine increment_enthalpy

end module heat_transfer_module
