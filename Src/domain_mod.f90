module domain_module

  use amrex_amr_module

  implicit none
  private
  
  public :: surfpos, tempinit, meltvel, surfdist ! description of surface alignment, constant scalar for now, multifab later. Need to interpolate between multifabs. 
  public :: get_face_velocity, create_face_flux, surface_tag, get_bound_heat
  
  
  real(amrex_real) :: surfpos, tempinit, meltvel     ! y surface position, for now constant given by input file
  real(amrex_real), allocatable :: surfdist(:)
  
  contains 
  
    subroutine get_face_velocity(time, xlo, dx, uface, ui_lo, ui_hi ) ! to be in separate module, using input from fluid solver and heat conduction solver 
      real(amrex_real), intent(in) :: dx(2), time, xlo(2) 					! grid size, time, and lower corner physical location 
      integer, intent(in) :: ui_lo(2), ui_hi(2)						! bounds of input tilebox 
      real(amrex_real), intent(inout) :: uface(ui_lo(1):ui_hi(1)+1,ui_lo(2):ui_hi(2)) 	! face velocity 
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
  real(amrex_real), intent(in   ) :: uin  (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2))  	! phi (temperature)
  real(amrex_real), intent(in   ) :: uface(ui_lo(1):ui_hi(1)+1,ui_lo(2):ui_hi(2)  ) 	! face velocity 
  logical         , intent(in   ) :: xfluxflag(ui_lo(1):ui_hi(1)+1,ui_lo(2):ui_hi(2)  ) 	! surface flag for x-nodes 
  logical         , intent(in   ) :: yfluxflag(ui_lo(1):ui_hi(1)  ,ui_lo(2):ui_hi(2)+1) 	! surface flag for y-nodes   			
  real(amrex_real), intent(  out) :: fluxx(ui_lo(1):ui_hi(1)+1,ui_lo(2):ui_hi(2)  ) 	! flux x direction  				
  real(amrex_real), intent(  out) :: fluxy(ui_lo(1):ui_hi(1)  ,ui_lo(2):ui_hi(2)+1) 	! flux y direction  	

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
  logical, intent(  out) :: xflux(ui_lo(1):ui_hi(1)+1,ui_lo(2):ui_hi(2)  ) 	! surface flag for x-nodes 
  logical, intent(  out) :: yflux(ui_lo(1):ui_hi(1)  ,ui_lo(2):ui_hi(2)+1) 	! surface flag for y-nodes 
  
  real(amrex_real) :: surf(ui_lo(1):ui_hi(1)) ! To be multifab with ghost points
  integer          :: surfind(ui_lo(1):ui_hi(1)) 
  integer :: i,j, jsurf
  

  do i = ui_lo(1),ui_hi(1) 
  surf(i) = surfpos! - & 
  !1_amrex_real*(i-ui_lo(1))*dx(1)!*(xlo(1) + (i-ui_lo(1))*dx(1) )
  
  		!(1+0.1*SIN(  (i-ui_lo(1))*dx(1) + xlo(1)  ))			! y-position of surface, array as given by fluid solver. 
  							! what form. Multifab? if so, need to interpolated between grids. 
  							! Need to 'know' surfpos at ui_lo(1)-1 
    
  end do 					
  
  !surf(ui_lo(1)  ) = surfpos + 0.1
  !surf(ui_lo(1)-1) = surfpos + 0.1
  surf(ui_lo(1)+1) = surfpos + 0.1
  
  !surf(ui_hi(1)  ) = surfpos + 0.1
  !surf(ui_hi(1)+1) = surfpos + 0.1
  surf(ui_hi(1)-1) = surfpos + 0.1
  
  !surf(ui_lo(1)+2) = surfpos + 0.075	
  !surf(ui_lo(1)+3) = surfpos + 0.05					
  !surf(ui_lo(1)+4) = surfpos + 0.025						
  					
  xflux = .false. 
  yflux = .false. 
  
  
  
  ! Preliminary flux application on surface 
  ! valid for single valued height as function of tangential coordinate 
  ! step-wise approximation
  
  
  ! Find surface index and flag surface edge for flux y-direction 
  do i = ui_lo(1),ui_hi(1) 
      jsurf = ui_lo(2) + floor(  (surf(i)-xlo(2))/dx(2) )
      surfind(i) = jsurf 		! y-index of surface element 
    !if (xlo(2) < surf(i)) then ! 
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
  
  
 
 
 
 
  subroutine get_bound_heat(time, xlo, dx, ui_lo, ui_hi, yflux, qb) 
  real(amrex_real), intent(in   ) :: time, xlo(2), dx(2)			! time, lower corner physical location, and grid size
  integer         , intent(in   ) :: ui_lo(2), ui_hi(2)			      	! bounds of input tilebox
  logical         , intent(in   ) :: yflux(ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2)+1) 	! surface flag for y-nodes  
  real(amrex_real), intent(  out) :: qb   (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2)  )    ! Volumetric heating localized to boundary 	
  integer :: i,j
  
  qb = 0. 
   do i = ui_lo(1), ui_hi(1) 
    do j = ui_lo(2), ui_hi(2)
    if(yflux(i,j)) then 
    qb(i,j-1) = 3_amrex_real/dx(2)   
    end if 
    end do  
   end do  
  
  end subroutine get_bound_heat
  
 
  
  
end module domain_module 

