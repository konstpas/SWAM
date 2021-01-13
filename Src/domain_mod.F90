module domain_module

  use amrex_amr_module

  implicit none
  private
  
  public :: tempinit, meltvel, surfdist 
  public :: get_face_velocity, create_face_flux, surface_tag, get_bound_heat, get_surf_pos
  
  
  real(amrex_real) :: tempinit, meltvel  
  real(amrex_real), allocatable :: surfdist(:)
  
  contains 
  
  
  
  
    subroutine get_face_velocity(time, xlo, dx, uface, & 
#if AMREX_SPACEDIM == 3  
				wface, &
#endif     
     				ui_lo, ui_hi )
     				
      real(amrex_real), intent(in) :: dx(3), time, xlo(3) 					! grid size, time, and lower corner physical location 
      integer, intent(in) :: ui_lo(3), ui_hi(3)						! bounds of input tilebox 
      real(amrex_real), intent(inout) :: uface(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! face velocity 
#if AMREX_SPACEDIM == 3  
      real(amrex_real), intent(inout) :: wface(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! face velocity 
#endif  
      integer 		:: i,j,k, yhigh, ylow  							! indexes and bounds 
      real(amrex_real)	:: yhighpos(ui_lo(1):ui_hi(1),ui_lo(3):ui_hi(3)), ylowpos (ui_lo(1):ui_hi(1),ui_lo(3):ui_hi(3))
      real(amrex_real) :: ypos 
      
      uface = 0_amrex_real 
#if AMREX_SPACEDIM == 3 
      wface = 0_amrex_real 
#endif  
      
      ! Some test velocity allocated 
      yhighpos = 0.90_amrex_real
      ylowpos = 0.70_amrex_real
      
      do i = ui_lo(1), ui_hi(1) 
      do k = ui_lo(3), ui_hi(3)

      	do j = ui_lo(2), ui_hi(2) 
      						! stagger		x   x   x  	
      	ypos = xlo(2) + (j-ui_lo(2))*dx(2)  	! index backwards 	!---!---!---!
      						! for now		 j-1  j  j+1   
      	if ((ypos < yhighpos(i,k)).and.(ypos > ylowpos(i,k))) then 
      		uface(i,j,k) = -10_amrex_real 
#if AMREX_SPACEDIM == 3 
      		wface(i,j,k) = -10_amrex_real 
#endif       		
      	end if 
      	
      	end do 
      	
      end do 
      end do 

  end subroutine get_face_velocity 
  
  
  
  
  
 
  
  ! Create enthalpy flux on edges in whole domain  
  subroutine create_face_flux(time, xlo, dx, &
  				uin, uface, ui_lo, ui_hi, &
  				xfluxflag,yfluxflag, &  
  				fluxx, fluxy &
#if AMREX_SPACEDIM == 3 
				, wface, zfluxflag, fluxz &
#endif 	
  				)
  				
  real(amrex_real), intent(in   ) :: time, xlo(3), dx(3)				! time, lower corner physical location, and grid size
  integer         , intent(in   ) :: ui_lo(3), ui_hi(3)				! bounds of input tilebox	
  real(amrex_real), intent(in   ) :: uin      (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3))  ! phi (temperature)
  real(amrex_real), intent(in   ) :: uface    (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! edge velocity x direction 
  logical         , intent(in   ) :: xfluxflag(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! surface flag for x-nodes 
  logical         , intent(in   ) :: yfluxflag(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! surface flag for y-nodes 		
  real(amrex_real), intent(  out) :: fluxx    (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! flux x direction  			
  real(amrex_real), intent(  out) :: fluxy    (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! flux y direction  	
#if AMREX_SPACEDIM == 3 
  real(amrex_real), intent(in   ) :: wface    (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! face velocity z direction 
  logical         , intent(in   ) :: zfluxflag(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! surface flag for z-nodes 
  real(amrex_real), intent(  out) :: fluxz    (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) ! flux z direction				
#endif   

  integer :: i,j,k 


  	! create face flux
  	! problems may arise since 2d or 1d flow not incompressible if vertical component neglected. 
  	! Also need to treat advection on empty-adjacent cells once surface is variable 


	! construct logical array for edges bounding domain (true if edge constitutes boundary)
	!call surface_tag(time, xlo, dx, ui_lo, ui_hi, xfluxflag, yfluxflag)

		! stagger		x   x   x  	
      		! index backwards 	!---!---!---!
      		! 			 i-1  i  i+1   
      	! Nodal enthalpy flux from temperature gradient and velocity field 	
	do   i = ui_lo(1), ui_hi(1) 
	 do  j = ui_lo(2), ui_hi(2) 
	  do k = ui_lo(3), ui_hi(3) 
	  
	 	if (uface(i,j,k) > 0_amrex_real) then 
			fluxx(i,j,k)  = uin(i-1,j,k)*uface(i,j,k)
		else 
			fluxx(i,j,k)  = uin(i  ,j,k)*uface(i,j,k)
		end if 
	   fluxx(i,j,k) = fluxx(i,j,k) - (uin(i,j,k)-uin(i-1,j  ,k))/dx(1)  ! x-velocity and temperature gradient 
	   fluxy(i,j,k) =              - (uin(i,j,k)-uin(i  ,j-1,k))/dx(2)  ! no y-velocity, temperature gradient  
	   
	   if(xfluxflag(i,j,k)) then ! true if center on either side of node is empty 
	    fluxx(i,j,k) = 0_amrex_real ! 
	   end if 
	   if(yfluxflag(i,j,k)) then ! true if surface node 
	    fluxy(i,j,k) = 0_amrex_real ! 
	   end if 
	   
#if AMREX_SPACEDIM == 3 
	 	if (wface(i,j,k) > 0_amrex_real) then 
			fluxz(i,j,k)  = uin(i,j,k-1)*wface(i,j,k)
		else 
			fluxz(i,j,k)  = uin(i,j,k  )*wface(i,j,k)
		end if 
	   fluxz(i,j,k) = fluxz(i,j,k) - (uin(i,j,k)-uin(i,j,k-1))/dx(3)  ! x-velocity and temperature gradient 
	   
	   if(zfluxflag(i,j,k)) then ! true if surface node 
	    fluxz(i,j,k) = 0_amrex_real ! 
	   end if 
#endif 	   
	   
	  end do  
	 end do
	end do  
		
  end subroutine create_face_flux 

	
 subroutine surface_tag(time, xlo, dx, ui_lo, ui_hi, xflux, yflux & 
#if AMREX_SPACEDIM == 3 
			, zflux)
#else  
 			)
#endif  
 
  real(amrex_real), intent(in   ) :: time, xlo(3), dx(3)			! time, lower corner physical location, and grid size
  integer, intent(in   ) :: ui_lo(3), ui_hi(3)			      	! bounds of input tilebox
  logical, intent(  out) :: xflux(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! surface flag for x-nodes 
  logical, intent(  out) :: yflux(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! surface flag for y-nodes 
#if AMREX_SPACEDIM == 3 
  logical, intent(  out) :: zflux(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! surface flag for z-nodes
#endif 
  real(amrex_real) :: surfpos(ui_lo(1):ui_hi(1),ui_lo(3):ui_hi(3)) 
  integer          :: surfind(ui_lo(1):ui_hi(1),ui_lo(3):ui_hi(3)) 
  integer :: i,j,k, jsurf
  
  
  
  
  
  call get_surf_pos(xlo, dx, ui_lo, ui_hi, surfpos)
											
  xflux = .false. 
  yflux = .false. 
#if AMREX_SPACEDIM == 3 
  zflux = .false. 
#endif   
  
  
  ! Preliminary flux application on surface 
  ! valid for single valued height as function of tangential coordinate 
  ! step-wise approximation
  
  
  ! Find surface index and flag surface edge for flux y-direction
   
  ! overlap between boxes here, since ui_lo,hi is used.  
  do i = ui_lo(1),ui_hi(1) 
  do k = ui_lo(3),ui_hi(3)
      jsurf = ui_lo(2) + floor(  (surfpos(i, k)-xlo(2))/dx(2) ) 
      surfind(i,k) = jsurf 		! y-index of surface element 
    if (jsurf >= ui_lo(2)) then 	! if surface is contained within this box
      	if (jsurf <= ui_hi(2)) then 	! if surface is contained within this box 
      	yflux(i, jsurf, k) = .true. 	! no diffusion across surface 
      	end if 
    end if   
  end do   
  end do 

  
  ! Flag surface edge for x-direction flux 
  do i = ui_lo(1)+1,ui_hi(1)
  do k = ui_lo(3)+1,ui_hi(3)-1 
 	if	(surfind(i,k) .gt. surfind(i-1,k)) then ! 		
 		do j = 	max(surfind(i-1,k),ui_lo(2)), & ! max, min since interface may be outside box 
 				min(surfind(i,k)-1,ui_hi(2)) ! surfind(i)-1 because surfind(i) edge is outside domain
 		  xflux(i,j,k) = .true. 
 		end do 
 	elseif	(surfind(i,k) .lt. surfind(i-1,k)) then !
 		do j = 	max(surfind(i,k),ui_lo(2)), & ! max, min since interface may be outside box 
 				min(surfind(i-1,k)-1,ui_hi(2)) ! -1 because surfind(i-1) edge is outside domain
 		  xflux(i,j,k) = .true. 
 		end do 
 	end if
  end do  
  end do 
  
#if AMREX_SPACEDIM == 3   
  ! Flag surface edge for z-direction flux 
  do i = ui_lo(1)+1,ui_hi(1)-1
  do k = ui_lo(3)+1,ui_hi(3) 
 	if	(surfind(i,k) .gt. surfind(i,k-1)) then ! 		
 		do j = 	max(surfind(i,k-1),ui_lo(2)), & ! 
 				min(surfind(i,k)-1,ui_hi(2)) !
 		  zflux(i,j,k) = .true. 
 		end do 
 	elseif	(surfind(i,k) .lt. surfind(i,k-1)) then !
 		do j = 	max(surfind(i,k),ui_lo(2)), & !
 				min(surfind(i,k-1)-1,ui_hi(2)) !
 		  zflux(i,j,k) = .true. 
 		end do 
 	end if
  end do  
  end do 
#endif   
  
  

 end subroutine surface_tag
  
  
 
 
 
 
 subroutine get_bound_heat(time, xlo, dx, lo, hi, ui_lo, ui_hi, yflux, qb) 
  real(amrex_real), intent(in   ) :: time, xlo(3), dx(3)			! time, lower corner physical location, and grid size
  integer         , intent(in   ) :: lo(3), hi(3)			      	! bounds of input tilebox  
  integer         , intent(in   ) :: ui_lo(3), ui_hi(3)			      	! bounds of input tilebox (ghost points)
  logical         , intent(in   ) :: yflux(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3)) 	! surface flag for y-nodes  
  real(amrex_real), intent(  out) :: qb(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))    ! Volumetric heating localized to boundary 	
  integer :: i,j,k
  
  
  qb = 0. 
  
  
  ! Test case for boundary heating 
   do   i = lo(1), hi(1) 
    do  j = lo(2), hi(2)
     do k = lo(3), hi(3) 
            !if ( (xlo(1)+(i-lo(1))*dx(1)).gt.0.5 ) then 
            !if ( (xlo(1)+(i-lo(1))*dx(1)).lt.1.5 ) then 
            
            if ( (xlo(3)+(k-lo(3))*dx(3)).gt.0.5 ) then 
            if ( (xlo(3)+(k-lo(3))*dx(3)).lt.1.5 ) then
            
      if(yflux(i,j+1,k)) then 
    qb(i,j,k) = 3_amrex_real/dx(2)   
      end if 
      
            end if 
            end if
      
            !end if 
            !end if
     end do        
    end do  
   end do  
  
 end subroutine get_bound_heat
  
  
  
 ! Subroutine to interpolate surface position as given by the fluid solver 
 ! in order to construct the heat conduction free interface  
 subroutine get_surf_pos(xlo, dx, ui_lo, ui_hi, surfpos)
  real(amrex_real), intent(in   ) :: xlo(3), dx(3)			                ! lower corner physical location, and grid size
  integer         , intent(in   ) :: ui_lo(3), ui_hi(3)			      	! bounds of input tilebox  
  real(amrex_real), intent(  out) :: surfpos(ui_lo(1):ui_hi(1),ui_lo(3):ui_hi(3))  ! Surface position with x and z coordinates 
  integer :: i,k
  
  
 ! Subroutine to interpolate surface position as given by the fluid solver 
 ! in order to construct the heat conduction free interface  
  do  i = ui_lo(1),ui_hi(1) 
   do k = ui_lo(3),ui_hi(3)
   surfpos(i,k) = 0.5_amrex_real + & 
   0.1_amrex_real*SIN(6.18_amrex_real*(xlo(1) + (i-ui_lo(1))*dx(1)) ) & 
                 *SIN(6.18_amrex_real*(xlo(3) + (k-ui_lo(3))*dx(3)) )
   end do 
  end do 
 
 
 
 end subroutine get_surf_pos 
 
 
  
  
end module domain_module 

