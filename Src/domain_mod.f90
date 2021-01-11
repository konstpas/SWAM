module domain_module

  use amrex_amr_module

  implicit none
  private
  
  public :: tempinit, meltvel, surfdist 
  public :: get_face_velocity, create_face_flux, surface_tag, get_bound_heat, get_surf_pos
  
  
  real(amrex_real) :: tempinit, meltvel  
  real(amrex_real), allocatable :: surfdist(:)
  
  contains 
  
    subroutine get_face_velocity(time, xlo, dx, uface, ui_lo, ui_hi ) ! to be in separate module, using input from fluid solver and heat conduction solver 
      real(amrex_real), intent(in) :: dx(2), time, xlo(2) 					! grid size, time, and lower corner physical location 
      integer, intent(in) :: ui_lo(2), ui_hi(2)						! bounds of input tilebox 
      real(amrex_real), intent(inout) :: uface(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2)) 	! face velocity 
      integer 		:: i,j, yhigh, ylow  							! indexes and bounds 
      real(amrex_real)	:: yhighpos(ui_lo(1):ui_hi(1)), ylowpos (ui_lo(1):ui_hi(1))
      real(amrex_real) :: ypos 
      
      uface = 0_amrex_real 
      !call sw_vel()
      yhighpos = 0.90_amrex_real
      ylowpos = 0.70_amrex_real
      
      do i = ui_lo(1), ui_hi(1) 
      !find yhigh,ylow (yhighpos(i)) 
      	do j = ui_lo(2), ui_hi(2) 
      						! stagger		x   x   x  	
      	ypos = xlo(2) + (j-ui_lo(2))*dx(2)  	! index backwards 	!---!---!---!
      						! for now		 j-1  j  j+1   
      	if ((ypos < yhighpos(i)).and.(ypos > ylowpos(i))) then 
      		uface(i,j) = -40_amrex_real 
      	end if 
      	
      	
      	end do 
      end do 
  	

  end subroutine get_face_velocity 
  
  
  
  
  
    
  
  
  
  ! Create enthalpy flux on all nodes 
  subroutine create_face_flux(time, xlo, dx, &
  				uin, uface, ui_lo, ui_hi, &
  				xfluxflag,yfluxflag, &  
  				fluxx, fluxy)
  				
  real(amrex_real), intent(in   ) :: time, xlo(2), dx(2)				! time, lower corner physical location, and grid size
  integer         , intent(in   ) :: ui_lo(2), ui_hi(2)				! bounds of input tilebox	
  real(amrex_real), intent(in   ) :: uin      (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2))  	! phi (temperature)
  real(amrex_real), intent(in   ) :: uface    (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2)) 	! face velocity 
  logical         , intent(in   ) :: xfluxflag(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2)) 	! surface flag for x-nodes 
  logical         , intent(in   ) :: yfluxflag(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2)) 	! surface flag for y-nodes   			
  real(amrex_real), intent(  out) :: fluxx    (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2)) 	! flux x direction  				
  real(amrex_real), intent(  out) :: fluxy    (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2)) 	! flux y direction  	

  integer :: i,j 


  	! create face flux
  	! problems may arise since 2d or 1d flow not incompressible if vertical component neglected. 
  	! Also need to treat advection on empty-adjacent cells once surface is variable 


	! construct logical array for edges bounding domain (true if edge constitutes boundary)
	!call surface_tag(time, xlo, dx, ui_lo, ui_hi, xfluxflag, yfluxflag)

		! stagger		x   x   x  	
      		! index backwards 	!---!---!---!
      		! 			 i-1  i  i+1   
      	! Nodal enthalpy flux from temperature gradient and velocity field 	
	do i = ui_lo(1),ui_hi(1) 
	 do j = ui_lo(2), ui_hi(2) 
	 	if (uface(i,j) > 0_amrex_real) then 
			fluxx(i,j)  = uin(i-1,j)*uface(i,j)
		else 
			fluxx(i,j)  = uin(i  ,j)*uface(i,j)
		end if 
	   fluxx(i,j) = fluxx(i,j) - (uin(i,j)-uin(i-1,j  ))/dx(1)  ! x-velocity and temperature gradient 
	   fluxy(i,j) =            - (uin(i,j)-uin(i  ,j-1))/dx(2)  ! no y-velocity, temperature gradient  
	   
	   if(xfluxflag(i,j)) then ! true if center on either side of node is empty 
	    fluxx(i,j) = 0_amrex_real ! 
	   end if 
	   if(yfluxflag(i,j)) then ! true if surface node 
	    fluxy(i,j) = 0_amrex_real ! 
	   end if 
	   
	 end do
	end do  

		! For boundaries. Flag all grid boxes as filled or empty. x-direction, if i or i-1 is empty, flux(i) = 0 and correspondingly for y-direction. 
		! y-direction add flux at boundary corresponding to external surface flux. 
	
	
		! Nodal flux at domain boundaries 
		
  end subroutine create_face_flux 

	
 subroutine surface_tag(time, xlo, dx, ui_lo, ui_hi, xflux, yflux)
  real(amrex_real), intent(in   ) :: time, xlo(2), dx(2)			! time, lower corner physical location, and grid size
  integer, intent(in   ) :: ui_lo(2), ui_hi(2)			      	! bounds of input tilebox
  logical, intent(  out) :: xflux(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2)) 	! surface flag for x-nodes 
  logical, intent(  out) :: yflux(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2)) 	! surface flag for y-nodes 
  
  real(amrex_real) :: surfpos(ui_lo(1):ui_hi(1)) 
  integer          :: surfind(ui_lo(1):ui_hi(1)) 
  integer :: i,j, jsurf
  
  
  
  
  
   call get_surf_pos(xlo, dx, ui_lo, ui_hi, surfpos)
											
  xflux = .false. 
  yflux = .false. 
  
  
  
  ! Preliminary flux application on surface 
  ! valid for single valued height as function of tangential coordinate 
  ! step-wise approximation
  
  
  ! Find surface index and flag surface edge for flux y-direction 
  do i = ui_lo(1),ui_hi(1) 
      jsurf = ui_lo(2) + floor(  (surfpos(i)-xlo(2))/dx(2) ) 
      surfind(i) = jsurf 		! y-index of surface element 
    if (jsurf >= ui_lo(2)) then 	! if surface is contained within this box
      	if (jsurf <= ui_hi(2)) then 	! if surface is contained within this box 
      	yflux(i, jsurf) = .true. 	! no diffusion across surface 
      	end if 
    end if   
  end do 
  ! Surface-index of grid point left-, and right-adjacent to box to be given as interpolation of surfpos

  
  ! Flag surface edge for x-direction flux 
  do i = ui_lo(1)+1,ui_hi(1)
 	if	(surfind(i) .gt. surfind(i-1)) then ! 
 				
 		do j = 	max(surfind(i-1),ui_lo(2)), & ! max, min since interface may be outside box 
 				min(surfind(i)-1,ui_hi(2)) ! surfind(i)-1 because surfind(i) edge is outside domain
 		  xflux(i,j) = .true. 
 		end do 
 		
 	elseif	(surfind(i) .lt. surfind(i-1)) then !
 	
 		do j = 	max(surfind(i),ui_lo(2)), & ! max, min since interface may be outside box 
 				min(surfind(i-1)-1,ui_hi(2)) ! -1 because surfind(i-1) edge is outside domain
 		  xflux(i,j) = .true. 
 		end do 
 		
 	end if
  end do  
  
  

 end subroutine surface_tag
  
  
 
 
 
 
 subroutine get_bound_heat(time, xlo, dx, lo, hi, ui_lo, ui_hi, yflux, qb) 
  real(amrex_real), intent(in   ) :: time, xlo(2), dx(2)			! time, lower corner physical location, and grid size
  integer         , intent(in   ) :: lo(2), hi(2)			      	! bounds of input tilebox  
  integer         , intent(in   ) :: ui_lo(2), ui_hi(2)			      	! bounds of input tilebox (ghost points)
  logical         , intent(in   ) :: yflux(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2)) 	! surface flag for y-nodes  
  real(amrex_real), intent(  out) :: qb(lo(1):hi(1),lo(2):hi(2))    ! Volumetric heating localized to boundary 	
  integer :: i,j
  
  qb = 0. 
   do i = lo(1), hi(1) 
    do j = lo(2), hi(2)
            if ( (xlo(1)+(i-lo(1))*dx(1)).gt.0.95 ) then 
            if ( (xlo(1)+(i-lo(1))*dx(1)).lt.1.05 ) then 
      if(yflux(i,j+1)) then 
    qb(i,j) = 3_amrex_real/dx(2)   
      end if 
            end if 
            end if
    end do  
   end do  
  
 end subroutine get_bound_heat
  
  
  
 ! Surface position interpolated to cell centers at interface  
 subroutine get_surf_pos(xlo, dx, ui_lo, ui_hi, surfpos)
  real(amrex_real), intent(in   ) :: xlo(2), dx(2)			! lower corner physical location, and grid size
  integer         , intent(in   ) :: ui_lo(2), ui_hi(2)			      	! bounds of input tilebox  
  real(amrex_real), intent(  out) :: surfpos(ui_lo(1):ui_hi(1))
  integer :: i
  
  do i = ui_lo(1),ui_hi(1) 
   surfpos(i) = 0.5_amrex_real + & 
   0.1_amrex_real*SIN(6.18_amrex_real*(xlo(1) + (i-ui_lo(1))*dx(1)) )
  end do 
 
 
 
 
 end subroutine get_surf_pos 
 
 
  
  
end module domain_module 

