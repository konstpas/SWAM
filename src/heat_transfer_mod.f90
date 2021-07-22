module heat_transfer_module

  use amrex_amr_module
  
  implicit none

  private

  public :: increment_enthalpy
  public :: get_face_velocity, create_face_flux, surface_tag, get_bound_heat
  public :: get_surf_pos, get_melt_pos, reset_melt_pos, integrate_surf
  !public :: get_idomain
  
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

    real(amrex_real) :: uface (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! face velocity x direction (nodal)
    real(amrex_real) :: wface (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! face velocity z direction (nodal)
    logical :: xfluxflag(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! surface flag for x-nodes 
    logical :: yfluxflag(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! surface flag for y-nodes
    logical :: zfluxflag(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! surface flag for z-nodes
    
    
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
    call surface_tag(time, geom%get_physical_location(ui_lo), dx, lo, hi, &
                     xfluxflag, yfluxflag, ui_lo, ui_hi, zfluxflag)	
  	
    ! Subroutine assigns tangential (x,z) velocity on grid edges in whole domain  
    call get_face_velocity(time, geom%get_physical_location(lo), dx, lo, hi, &
                           uface, wface, ui_lo, ui_hi ) 	
 			 	

    ! Subroutine assigns enthalpy flux on grid edges in whole domain 
    call create_face_flux(time, geom%get_physical_location(ui_lo), dx, lo, hi, 	        &
                          uin, ui_lo, ui_hi,              &
                          flxx, fx_lo, fx_hi,             &
                          flxy, fy_lo, fy_hi,             &
                          flxz, fz_lo, fz_hi,             &
                          tempin, ti_lo, ti_hi,           &
                          uface, wface,                   &
                          xfluxflag, yfluxflag, zfluxflag)
  				  	
    ! Zero flux across surface boundary. Volumetric heat deposition in first internal cell constitutes absorbed boundary flux. 
    ! Incorporates all absorption and cooling terms 			
    call get_bound_heat(time, geom%get_physical_location(lo), dx, lo, hi, ui_lo, ui_hi, yfluxflag, qbound) 	! domain module 			
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


  subroutine get_face_velocity(time, xlo, dx, lo, hi, uface, & 
                               wface, ui_lo, ui_hi )
    
    real(amrex_real), intent(in) :: dx(3), time, xlo(3) 					! grid size, time, and lower corner physical location 
    integer, intent(in) :: lo(3), hi(3)						! bounds of input tilebox       
    integer, intent(in) :: ui_lo(3), ui_hi(3)					! bounds of input box 
    real(amrex_real), intent(inout) :: uface(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! face velocity 
    real(amrex_real), intent(inout) :: wface(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! face velocity 
    integer 		:: i,j,k, yhigh, ylow  							! indexes and bounds 
    real(amrex_real)	:: yhighpos(ui_lo(1):ui_hi(1),ui_lo(3):ui_hi(3)), ylowpos (ui_lo(1):ui_hi(1),ui_lo(3):ui_hi(3))
    real(amrex_real) :: ypos 
    
    uface = 0_amrex_real 
    wface = 0_amrex_real 
    
    ! Some test velocity allocated 
    !      yhighpos = 0.90_amrex_real
    !      ylowpos = 0.70_amrex_real
    
    !      do i = lo(1), hi(1)+1  ! +1 because of staggering 
    !      do k = lo(3), hi(3)+1  ! +1 because of staggering 
    
    !      	do j = lo(2), hi(2) 
    !      						! stagger		x   x   x  	
    !      	ypos = xlo(2) + (j-lo(2))*dx(2)  	! index backwards 	!---!---!---!
    ! for now		 j-1  j  j+1   
    !      	if ((ypos < yhighpos(i,k)).and.(ypos > ylowpos(i,k))) then 
    !      		uface(i,j,k) = -10_amrex_real 
    !      		wface(i,j,k) = -10_amrex_real 
    !      	end if 
    
    !      	end do 
    
    !      end do 
    !      end do 
    
  end subroutine get_face_velocity
  
  
  

  
  ! Create enthalpy flux on edges in whole domain
  subroutine create_face_flux(time, xlo, dx, lo, hi, 	        &
                                uin, ui_lo, ui_hi,              &
                                flxx, fx_lo, fx_hi,             &
                                flxy, fy_lo, fy_hi,             &
                                flxz, fz_lo, fz_hi,             &
                                temp, t_lo, t_hi,               &
                                uface, wface,                   &
                                xfluxflag, yfluxflag, zfluxflag)
  				
    use material_properties_module, only: get_ktherm
    
    real(amrex_real), intent(in   ) :: time, xlo(3), dx(3)				! time, lower corner physical location, and grid size
    integer         , intent(in   ) :: lo(3), hi(3)					! bounds of input tilebox	  
    integer         , intent(in   ) :: ui_lo(3), ui_hi(3)				! bounds of input tilebox
    integer         , intent(in   ) :: t_lo(3), t_hi(3)				! bounds of input tilebox
    integer         , intent(in   ) :: fx_lo(3), fx_hi(3)				! bounds of input tilebox
    integer         , intent(in   ) :: fy_lo(3), fy_hi(3)				! bounds of input tilebox
    integer         , intent(in   ) :: fz_lo(3), fz_hi(3)				! bounds of input tilebox
    
    real(amrex_real), intent(in   ) :: uin      (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! phi (enthalpy)		
    real(amrex_real), intent(  out) :: flxx    (fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3)) ! flux x direction  			
    real(amrex_real), intent(  out) :: flxy    (fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3)) ! flux y direction
    real(amrex_real), intent(  out) :: flxz    (fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3)) ! flux z direction				
    real(amrex_real), intent(in   ) :: temp     (t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))  ! temperature
    real(amrex_real), intent(in) :: uface (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! face velocity x direction (nodal)
    real(amrex_real), intent(in) :: wface (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! face velocity z direction (nodal)
    logical, intent(in) :: xfluxflag(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! surface flag for x-nodes 
    logical, intent(in) :: yfluxflag(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! surface flag for y-nodes 
    logical, intent(in) :: zfluxflag(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! surface flag for z-nodes   
    
    integer :: i,j,k 
    real(amrex_real) :: ktherm, temp_face
    
    
    ! create face flux
    ! problems may arise since 2d or 1d flow not incompressible if vertical component neglected. 
    ! Also need to treat advection on empty-adjacent cells once surface is variable 
    
    ! construct logical array for edges bounding domain (true if edge constitutes boundary)
    !call surface_tag(time, xlo, dx, ui_lo, ui_hi, xfluxflag, yfluxflag)
    
    ! stagger		x   x   x  	
    ! index backwards 	!---!---!---!
    ! 			 i-1  i  i+1   
    ! Nodal enthalpy flux from temperature gradient and velocity field 	
    ! do i = lo(1), hi(1)+1  
    !    do j = lo(2), hi(2)
    !       do k = lo(3), hi(3)
   
    !            ! x direction flux 	
    !            temp_face = (temp(i,j,k) + temp(i-1,j,k))/2_amrex_real ! temperature on face 
    !            call get_ktherm(temp_face,ktherm)
    !            fluxx(i,j,k) = 0_amrex_real
    !            fluxx(i,j,k) = fluxx(i,j,k) - ktherm*(temp(i,j,k)-temp(i-1,j  ,k))/dx(1)  ! x-velocity and temperature gradient 
               
    !         end do
    !      end do
    !   end do
    
    do i = lo(1), hi(1)+1
       do j = lo(2), hi(2)
          do k = lo(3), hi(3)

             if (uface(i,j,k) > 0_amrex_real) then 
                flxx(i,j,k)  = uin(i-1,j,k)*uface(i,j,k)
             else 
                flxx(i,j,k)  = uin(i  ,j,k)*uface(i,j,k)
             end if
             
             temp_face = (temp(i,j,k) + temp(i-1,j,k))/2_amrex_real
             call get_ktherm(temp_face, ktherm)
             flxx(i,j,k) = flxx(i,j,k) - ktherm*(temp(i,j,k)-temp(i-1,j,k))/dx(1)

             if(xfluxflag(i,j,k)) then 
                flxx(i,j,k) = 0_amrex_real 
             end if
  
          end do
       end do
    end do

    do i = lo(1), hi(1)
       do j = lo(2), hi(2)+1
          do k = lo(3), hi(3)
             
             temp_face = (temp(i,j,k) + temp(i,j-1,k))/2_amrex_real
             call get_ktherm(temp_face, ktherm)
             flxy(i,j,k) = -ktherm*(temp(i,j,k)-temp(i,j-1,k))/dx(2)

             if(yfluxflag(i,j,k)) then 
                flxy(i,j,k) = 0_amrex_real 
             end if
  
          end do
       end do
    end do

    do i = lo(1), hi(1)
       do j = lo(2), hi(2)
          do k = lo(3), hi(3)+1

             if (wface(i,j,k) > 0_amrex_real) then 
                flxz(i,j,k)  = uin(i,j,k-1)*wface(i,j,k)
             else 
                flxz(i,j,k)  = uin(i,j,k  )*wface(i,j,k)
             end if
             
             temp_face = (temp(i,j,k) + temp(i,j,k-1))/2_amrex_real
             call get_ktherm(temp_face, ktherm)
             flxz(i,j,k) = flxz(i,j,k) - ktherm*(temp(i,j,k)-temp(i,j,k-1))/dx(3)

             if(zfluxflag(i,j,k)) then 
                flxz(i,j,k) = 0_amrex_real 
             end if
  
          end do
       end do
    end do
    
    
  end subroutine create_face_flux

	
    subroutine surface_tag(time, xlo, dx, lo, hi, &
 		           xflux, yflux, ui_lo, ui_hi, & 
                           zflux)
 
      real(amrex_real), intent(in   ) :: time, xlo(3), dx(3)			! time, lower corner physical location, and grid size
      integer, intent(in   ) :: lo(3), hi(3)			      		! bounds of input tilebox  
      integer, intent(in   ) :: ui_lo(3), ui_hi(3)			      	! bounds of input box
      logical, intent(  out) :: xflux(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! surface flag for x-nodes 
      logical, intent(  out) :: yflux(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! surface flag for y-nodes 
      logical, intent(  out) :: zflux(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! surface flag for z-nodes
      real(amrex_real) :: surfpos(ui_lo(1):ui_hi(1),ui_lo(3):ui_hi(3)) ! 
      integer          :: surfind(ui_lo(1):ui_hi(1),ui_lo(3):ui_hi(3)) 
      integer :: i,j,k, jsurf
      
      ! Inefficient as of now, since it fills whole box for every called tile (every box has several tiles) 
      call get_surf_pos(xlo, dx, ui_lo, ui_hi, surfpos)  
      
      xflux = .false. 
      yflux = .false. 
      zflux = .false. 
      
      
      ! Preliminary flux application on surface 
      ! valid for single valued height as function of tangential coordinate 
      ! step-wise approximation
      
      
      ! Find surface index and flag surface edge for flux y-direction
      
      do i = lo(1)-1, hi(1)+1 ! Including ghost points 
         do k = lo(3)-1, hi(3)+1 ! Including ghost points (needs rewrite for 2D) 

            jsurf = ui_lo(2) + floor(  (surfpos(i, k)-xlo(2))/dx(2) ) 
            surfind(i,k) = jsurf 		! y-index of surface element
            
            if (jsurf >= lo(2)) then 	! if surface is contained within this box

               if (jsurf <= hi(2)+1) then 	! if surface is contained within this box. +1 because of staggering 
                  yflux(i, jsurf, k) = .true. 	! no diffusion across surface 
               end if
               
            end if
            
         end do
      end do



      ! Flag surface edge for x-direction flux 
      do i = lo(1),hi(1)+1  ! +1 because of grid staggering 
         do k = lo(3),hi(3)
            
            if	(surfind(i,k) .gt. surfind(i-1,k)) then
               
               do j = 	max(surfind(i-1,k),lo(2)), & ! max, min since interface may be outside box 
                    min(surfind(i,k)-1,hi(2)) ! surfind(i)-1 because surfind(i) edge is outside domain
 		  xflux(i,j,k) = .true. 
               end do
  
            elseif	(surfind(i,k) .lt. surfind(i-1,k)) then
               
 		do j = 	max(surfind(i,k),lo(2)), & ! max, min since interface may be outside box 
 				min(surfind(i-1,k)-1,hi(2)) ! -1 because surfind(i-1) edge is outside domain
 		  xflux(i,j,k) = .true. 
                end do
  
            end if
             
         end do
      end do
  

  
      ! Flag surface edge for z-direction flux 
      do i = lo(1),hi(1)
         do k = lo(3),hi(3)+1   ! +1 because of grid staggering
            
            if	(surfind(i,k) .gt. surfind(i,k-1)) then
               
 		do j = 	max(surfind(i,k-1),lo(2)), & 
 				min(surfind(i,k)-1,hi(2)) 
 		  zflux(i,j,k) = .true. 
                end do
  
             elseif	(surfind(i,k) .lt. surfind(i,k-1)) then
                
 		do j = 	max(surfind(i,k),lo(2)), & 
 				min(surfind(i,k-1)-1,hi(2)) 
 		  zflux(i,j,k) = .true. 
                end do
  
             end if
             
          end do
       end do
  

     end subroutine surface_tag
  
  
 
 
 
 
     subroutine get_bound_heat(time, xlo, dx, lo, hi, ui_lo, ui_hi, yflux, qb) 

       use amr_data_module, only : surf_pos
       use read_input_module, only : exp_time, flux_peak, flux_width, flux_pos
       
       real(amrex_real), intent(in   ) :: time, xlo(3), dx(3)			! time, lower corner physical location, and grid size
       integer         , intent(in   ) :: lo(3), hi(3)			      	! bounds of input tilebox  
       integer         , intent(in   ) :: ui_lo(3), ui_hi(3)			      	! bounds of input tilebox (ghost points)
       logical         , intent(in   ) :: yflux(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! surface flag for y-nodes
       !integer, intent(in) :: yflux ! DUMMY DECLARATION (NOT USED FOR NOW)
       real(amrex_real), intent(  out) :: qb(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))    ! Volumetric heating localized to boundary 	
       real(amrex_real) :: xpos, zpos, ypos
       integer :: i,j,k
    
       qb = 0. 
    
       ! Boundary heating 
       do   i = lo(1), hi(1) 
          do  j = lo(2), hi(2)
             do k = lo(3), hi(3)
                
                if (time.lt.exp_time) then
                   
                   if(yflux(i,j+1,k)) then 
                      
                      xpos = xlo(1) + (i-lo(1))*dx(1)
                      zpos = xlo(3) + (k-lo(3))*dx(3)
                      
                      qb(i,j,k) = flux_peak*EXP( 	&
                           -((xpos-flux_pos(1))**2)/(flux_width(1)**2)	&
                           -((zpos-flux_pos(2))**2)/(flux_width(2)**2))/dx(2)   
                   end if
                   
                end if
                
             end do
          end do
       end do

       ! ! Boundary heating 
       ! do   i = lo(1), hi(1) 
       !    do  k = lo(3), hi(3)
       !       do j = lo(2), hi(2)

       !          if (xlo(2).lt.surf_pos(i,k)) then

       !             if (time.lt.exp_time) then

       !                ypos = xlo(2) + (j+1-lo(2))*dx(2)
                   
       !                if(ypos.gt.surf_pos(i,k)) then 
                         
       !                   xpos = xlo(1) + (i-lo(1))*dx(1)
       !                   zpos = xlo(3) + (k-lo(3))*dx(3)
                         
       !                   qb(i,j,k) = flux_peak*EXP( 	&
       !                        -((xpos-flux_pos(1))**2)/(flux_width(1)**2)	&
       !                        -((zpos-flux_pos(2))**2)/(flux_width(2)**2))/dx(2)
       !                   !qb(i,j,k) = flux_peak
                         
       !                   exit
                         
       !                end if
                   
       !             end if

       !          else

       !             exit
                 
       !          end if
                
       !       end do
       !    end do
       ! end do
       
     end subroutine get_bound_heat
  
  
  
     ! Subroutine to interpolate surface position as given by the fluid solver 
     ! in order to construct the heat conduction free interface  
     subroutine get_surf_pos(xlo, dx, ui_lo, ui_hi, surfpos)

       use amr_data_module, only : surf_ind, surf_pos, surf_xlo, surf_dx  
       
       real(amrex_real), intent(in   ) :: xlo(3), dx(3)			                ! lower corner physical location, and grid size
       integer         , intent(in   ) :: ui_lo(3), ui_hi(3)			      	 ! bounds of input tilebox  
       real(amrex_real), intent(  out) :: surfpos(ui_lo(1):ui_hi(1),ui_lo(3):ui_hi(3))     ! Surface position with x and z coordinates, at grid size   
       integer :: i, k, xind, zind
       real(amrex_real) :: xpos, zpos, x_alpha, z_alpha, valzp, valzm, valxp, valxm   
  
	  
       do  i = ui_lo(1),ui_hi(1) 
          do k = ui_lo(3),ui_hi(3)
   
             xpos = xlo(1) + (0.5 + i-ui_lo(1))*dx(1) 
             zpos = xlo(3) + (0.5 + k-ui_lo(3))*dx(3)
   
             xind = ceiling( (xpos - surf_dx(1)/2 - surf_xlo(1))/surf_dx(1)  ) ! -surf_dx(1)/2, and ceiling ceiling since staggered 'backwards' on faces w.r.t values which are centered. 
             x_alpha = mod(xpos - surf_dx(1)/2 - surf_xlo(1), surf_dx(1))
             zind = ceiling( (zpos - surf_dx(2)/2 - surf_xlo(2))/surf_dx(2)  ) !
             z_alpha = mod(zpos - surf_dx(2)/2 - surf_xlo(2), surf_dx(2))
   
             if (xind.lt.surf_ind(1,1)) xind = surf_ind(1,1)
             if (xind.ge.surf_ind(1,2)) xind = surf_ind(1,2)-1 
             if (zind.lt.surf_ind(2,1)) zind = surf_ind(2,1)
             if (zind.ge.surf_ind(2,2)) zind = surf_ind(2,2)-1 
   
   
             valzm = surf_pos(xind,zind  ) + x_alpha * (surf_pos(xind+1,zind  )-surf_pos(xind,zind  )) ! interpolated value at zind	
             valzp = surf_pos(xind,zind+1) + x_alpha * (surf_pos(xind+1,zind+1)-surf_pos(xind,zind+1)) ! int value at zind+1 
   
             surfpos(i,k) = valzm + z_alpha*(valzp - valzm)  ! 2D linear interpolation 

          end do
       end do
 
 
 
     end subroutine get_surf_pos
 
 
 
     subroutine reset_melt_pos()	

       use amr_data_module, only : surf_pos, melt_pos, surf_ind    ! 2D array of melt position 
       integer :: i,k 
  
       do i =  surf_ind(1,1), surf_ind(1,2) 
          do k = surf_ind(2,1), surf_ind(2,2)
             melt_pos(i,k) = surf_pos(i,k)
          end do
       end do
 
     end subroutine reset_melt_pos
 
     subroutine get_melt_pos(lo, hi, temp, t_lo, t_hi, geom)
       
       use amr_data_module, only : surf_pos, melt_pos   ! 2D array of melt position 
       use material_properties_module, only : melt_point

       integer,              intent(in) :: lo(3), hi(3)  		! bounds of current tile box
       integer,              intent(in) :: t_lo(3), t_hi(3)		! bounds of temperature box    
       real(amrex_real),     intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) ! Temperature out 
       type(amrex_geometry), intent(in) :: geom  	! geometry at level
       integer :: i,j,k, it(1:3) 
       real(amrex_real) :: grid_pos(1:3)

       
       do i = lo(1), hi(1)  ! x-direction
          do k = lo(3), hi(3)  ! z-direction 	
             do j = lo(2), hi(2) 

		if (temp(i,j,k).gt.melt_point) then
     
                   it(1) = i
		   it(2) = j
		   it(3) = k 
 		   grid_pos = geom%get_physical_location(it)
                   if (grid_pos(2).lt.melt_pos(i,k)) then    
                      melt_pos(i,k) = grid_pos(2) 
                   end if
                   
		end if
	
             end do
          end do
       end do

     end subroutine get_melt_pos
 
	  
     subroutine integrate_surf(melt_vol)	

       use amr_data_module, only : surf_pos, melt_pos, surf_ind, surf_dx   ! 2D array of surface position, melt position. Index space for surface and grid size for surface 
       real(amrex_real), intent(out) :: melt_vol ! Integrated melt volume [mm3] 
       integer :: i,k 
 
       melt_vol = 0 
       
       do i =  surf_ind(1,1), surf_ind(1,2) 
          do k = surf_ind(2,1), surf_ind(2,2)
             melt_vol = melt_vol +  surf_pos(i,k) - melt_pos(i,k)
          end do
       end do
 
       melt_vol = melt_vol*surf_dx(1)*surf_dx(2)*1E9  


     end subroutine integrate_surf
  
end module heat_transfer_module
