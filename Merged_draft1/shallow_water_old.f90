module shallow_water_module  
  use amrex_amr_module
  use amr_data_module, only : surf_ind, surf_xlo, surf_dx, surf_pos, melt_pos, melt_vel, height_flux 

  implicit none 

  contains 
 
! surfpos initialization, for SW testing only (to be removed)  
 subroutine init_surfpos()
 use domain_module, only: flux_pos, flux_width
 real(amrex_real) :: xpos, zpos 
 integer :: i,k 
   do  i = surf_ind(1,1),surf_ind(1,2)
    do k = surf_ind(2,1),surf_ind(2,2)
    xpos = surf_xlo(1) + surf_dx(1)*i
    zpos = surf_xlo(2) + surf_dx(2)*k
     surf_pos(i,k) = surf_pos(i,k)*(1 + 0.1*EXP(	&
     			-((xpos-flux_pos(1))**2)/(flux_width(1)**2) & 
     			-((zpos-flux_pos(2))**2)/(flux_width(2)**2) & 	
     ))
    end do 
   end do 
    				
 end subroutine init_surfpos 
 
 
  
  subroutine increment_SW(time, geom, dt)
    type(amrex_geometry), intent(in) :: geom  	! geometry at highest level  
    real(amrex_real), intent(in) :: dt, time		! sub time step, and time 
    real(amrex_real) :: melt_height(surf_ind(1,1):surf_ind(1,2),surf_ind(2,1):surf_ind(2,2))
    integer :: i,k 
    
 
  ! surf_ind(:,:) are the x/z, begin/end (1/2, 1/2) indexes for the 2D surface domain 
  ! surf_xlo(:) are the x/z lowest coordinate 
  ! surf_dx(:) are the x/z grid sizes   
 
 
  ! Momentum continuity equation 
  ! Melt velocity constant for testing 
  ! To be given by momentum cont. 
  melt_vel = 0.05 
	
  
  
  height_flux = 0. 
  
  ! Mass (column height) continuity equation 
  
  
  ! Find 'height flux' 
  
  ! Surf_pos is initialized to a gaussian (init_SW above, called in main script). In final implementation, surf_pos is initialized 
  ! surf_pos = melt_pos = init_surf_pos, i.e. there is no melting and the surface is flat. 
  ! 
  ! Melt (liquid) column height is the distance between surface and liquid/solid interfaces 
  ! surf_pos, melt_pos are declared in amr_data_mod, and allocated in my_amr_init of my_amr_mod. 
  ! melt_pos is found every time step in the heat solver substep  
  
  melt_height = surf_pos-melt_pos 
  
  
  
  ! height_flux is declared in amr_data_mod, and allocated in my_amr_init of my_amr_mod.
  ! X flux 
  do  i = surf_ind(1,1),surf_ind(1,2)+1
   do k = surf_ind(2,1),surf_ind(2,2)
   	if (i.eq.surf_ind(1,1)) then ! low boundary
		if (melt_vel(i,k,1) > 0_amrex_real) then 
			height_flux(i,k,1) = 0 !no influx
		else 
			height_flux(i,k,1) = melt_height(i,k)*melt_vel(i,k,1)
		end if 
   	elseif (i.eq.surf_ind(1,2)+1) then ! high boundary 
   	
		if (melt_vel(i,k,1) > 0_amrex_real) then 
			height_flux(i,k,1) = melt_height(i-1,k)*melt_vel(i,k,1)
		else 
			height_flux(i,k,1) = 0 !no influx
		end if 
   	else 	
		if (melt_vel(i,k,1) > 0_amrex_real) then 
			height_flux(i,k,1) = melt_height(i-1,k)*melt_vel(i,k,1)
		else 
			height_flux(i,k,1) = melt_height(i,k)*melt_vel(i,k,1)
		end if 
	end if 
   end do 
  end do 
  ! Z flux 
  do  i = surf_ind(1,1),surf_ind(1,2)
   do k = surf_ind(2,1),surf_ind(2,2)+1
	if     (k.eq.surf_ind(2,1)  ) then 
		if (melt_vel(i,k,2) > 0_amrex_real) then 
		 height_flux(i,k,2) = 0! no influx
		else 
		 height_flux(i,k,2) = melt_height(i,k)*melt_vel(i,k,2)
		end if 
	elseif (k.eq.surf_ind(2,2)+1) then 
		if (melt_vel(i,k,2) > 0_amrex_real) then 
		 height_flux(i,k,2) = melt_height(i,k-1)*melt_vel(i,k,2)
		else 
		 height_flux(i,k,2) = 0 !no influx  
		end if 
	else 
		if (melt_vel(i,k,2) > 0_amrex_real) then 
		 height_flux(i,k,2) = melt_height(i,k-1)*melt_vel(i,k,2)
		else 
		 height_flux(i,k,2) = melt_height(i,k)*melt_vel(i,k,2)
		end if 
	end if 
   end do 
  end do 
  	
  ! Update surface position by divergence of 'column height flux' 			
   do  i = surf_ind(1,1),surf_ind(1,2)
    do k = surf_ind(2,1),surf_ind(2,2)
     surf_pos(i,k) = surf_pos(i,k) & 
     		     - dt/surf_dx(1) * (height_flux(i+1,k  ,1) - height_flux(i,k,1)) & 
     		     - dt/surf_dx(2) * (height_flux(i  ,k+1,2) - height_flux(i,k,2))		
    end do 
   end do 		
		
		
		
		
		
  
  end subroutine increment_SW
 





end module shallow_water_module 
