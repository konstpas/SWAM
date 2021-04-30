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

 surf_pos(surf_ind(1,1):surf_ind(1,2)/2,surf_ind(2,1):surf_ind(2,2)) = 0.2
 surf_pos(surf_ind(1,2)/2 + 1 :surf_ind(1,2),surf_ind(2,1):surf_ind(2,2)) = 0.

!   do  i = surf_ind(1,1),surf_ind(1,2)
!    do k = surf_ind(2,1),surf_ind(2,2)
!    xpos = surf_xlo(1) + surf_dx(1)*i
!    zpos = surf_xlo(2) + surf_dx(2)*k
!     surf_pos(i,k) = surf_pos(i,k)*(1 + 0.1*EXP(                    &
!                        -((xpos-flux_pos(1))**2)/(flux_width(1)**2) & 
!                        -((zpos-flux_pos(2))**2)/(flux_width(2)**2) & 
!     ))
       


!    end do 
!   end do 
 end subroutine init_surfpos 
 
 
  
  subroutine increment_SW(time, geom, dt)
    type(amrex_geometry), intent(in) :: geom ! geometry at highest level  
    real(amrex_real), intent(in) :: dt, time ! sub time step, and time 
    real(amrex_real) :: melt_height(surf_ind(1,1):surf_ind(1,2),surf_ind(2,1):surf_ind(2,2))
    integer :: i,k 
    real(amrex_real) :: h,hm,u,um,hu,hum,c,vp,v,vm,vmm,gradvm,gradvp,hv,hvm,up,umm,gradum,gradup
    real(amrex_real) :: dh,dhu,dhv,fqp,fqm,fvp,fvm 
    real(amrex_real), dimension(surf_ind(1,1):surf_ind(1,2)+1,surf_ind(2,1):surf_ind(2,2)) :: fh_x,fq_x,fv_x
    real(amrex_real), dimension(surf_ind(1,1):surf_ind(1,2),surf_ind(2,1):surf_ind(2,2)+1) :: fh_y,fq_y,fv_y
    real(amrex_real) :: g,small

 

  ! surf_ind(:,:) are the x/z, begin/end (1/2, 1/2) indexes for the 2D surface domain 
  ! surf_xlo(:) are the x/z lowest coordinate 
  ! surf_dx(:) are the x/z grid sizes   
 
 
  ! Momentum continuity equation 
  ! Melt velocity constant for testing 
  ! To be given by momentum cont. 
  
  ! Find 'height flux' 
  
  ! Surf_pos is initialized to a gaussian (init_SW above, called in main script). In final implementation, surf_pos is initialized 
  ! surf_pos = melt_pos = init_surf_pos, i.e. there is no melting and the surface is flat. 
  ! 
  ! Melt (liquid) column height is the distance between surface and liquid/solid interfaces 
  ! surf_pos, melt_pos are declared in amr_data_mod, and allocated in my_amr_init of my_amr_mod. 
  ! melt_pos is found every time step in the heat solver substep  
  small = 1.*1e-10 

  g = 1. 
  melt_height = surf_pos-melt_pos 
  
 !computing the fluxes at surfaces using the Rusanov riemann solver
 !check damb example at basilisk.fr or popinet 2011 
 !we can either use Kurganov riemann solver or rusanov to compute fluxes; check riemann.h at basilisk.fr for detailed description of Kurganov solver
 
  do  i = surf_ind(1,1),surf_ind(1,2) + 1
    do k = surf_ind(2,1),surf_ind(2,2)
      
 !computing the fluxes at the boundaries; using the neumann boundary condition at all the boundaries for h,u,v 
 !fh_x represent hu and fq_x represent huu+gravity term
 !fv_x represent x-dir momentum (hu) projected by v 
      if (i.eq.surf_ind(1,1)) then
        h   = melt_height(i,k)    
        hm  = h
        u   = 0.5*( melt_vel(i,k,1) + melt_vel(i+1,k,1) ) 
        um  = u
        hu  = h  * u
        hum = hm * um
       
        c = max( (abs(um)+(g*hm)**0.5) , (abs(u)+(g*h)**0.5) )

        fh_x(i,k) = 0.5*(hum + hu) - 0.5*c*(h-hm)
        fq_x(i,k) = 0.5*(hum**2 + 0.5*g*hm**2 + hu**2 + 0.5*g*h**2) - 0.5*c*(hu-hum)

        vp  = 0.5*( melt_vel(i+1,k,2) + melt_vel(i+1,k+1,2) )
        v   = 0.5*( melt_vel(i,k,2) + melt_vel(i,k+1,2) )
        vm  = v
        gradvm = 0.
        gradvp = vp - v

        if (fh_x(i,k).gt.0.) fv_x(i,k) = ( vm + 0.5*gradvm )*fh_x(i,k)
        if (fh_x(i,k).lt.0.) fv_x(i,k) = ( v  - 0.5*gradvp )*fh_x(i,k)
      elseif (i.eq.surf_ind(1,2)+1) then
        h   = melt_height(i,k)
        hm  = melt_height(i,k)
        u   = 0.5*( melt_vel(i,k,1) + melt_vel(i-1,k,1) )
        um  = 0.5*( melt_vel(i,k,1) + melt_vel(i-1,k,1) )
        hu  = h  * u
        hum = hm * um
     
        c = max( (abs(um)+(g*hm)**0.5) , (abs(u)+(g*h)**0.5) )

        fh_x(i,k) = 0.5*(hum + hu) - 0.5*c*(h-hm)
        fq_x(i,k) = 0.5*(hum**2 + 0.5*g*hm**2 + hu**2 + 0.5*g*h**2) - 0.5*c*(hu-hum)

        v   = 0.5*( melt_vel(i,k,2) + melt_vel(i,k+1,2) )
        vm  = 0.5*( melt_vel(i-1,k,2) + melt_vel(i-1,k+1,2) )
        vmm = 0.5*( melt_vel(i-2,k,2) + melt_vel(i-2,k+1,2) )
        gradvp  = 0.
        gradvm  = vm  - vmm
          
        if (fh_x(i,k).gt.0.) fv_x(i,k) = ( vm + 0.5*gradvm )*fh_x(i,k)
        if (fh_x(i,k).lt.0.) fv_x(i,k) = ( v  - 0.5*gradvp )*fh_x(i,k)
      else
 !computing the fluxes within the domain
       
        h   = melt_height(i,k)
        hm  = melt_height(i-1,k)
        u   = 0.5*( melt_vel(i,k,1) + melt_vel(i+1,k,1) )
        um  = 0.5*( melt_vel(i,k,1) + melt_vel(i-1,k,1) )
        hu  = h  * u
        hum = hm * um

        c = max( (abs(um)+(g*hm)**0.5) , (abs(u)+(g*h)**0.5) )

        fh_x(i,k) = 0.5*(hum + hu) - 0.5*c*(h-hm)
        fq_x(i,k) = 0.5*(hum**2 + 0.5*g*hm**2 + hu**2 + 0.5*g*h**2) - 0.5*c*(hu-hum)

        if (i.eq.surf_ind(1,1)+1) then
          vp  = 0.5*( melt_vel(i+1,k,2) + melt_vel(i+1,k+1,2) )
          v   = 0.5*( melt_vel(i,k,2) + melt_vel(i,k+1,2) )
          vm  = 0.5*( melt_vel(i-1,k,2) + melt_vel(i-1,k+1,2) )
          gradvp  = vp  - v
          gradvm  = 0.
        elseif (i.eq.surf_ind(1,2)) then
          v    = 0.5*( melt_vel(i,k,2) + melt_vel(i,k+1,2) )
          vm   = 0.5*( melt_vel(i-1,k,2) + melt_vel(i-1,k+1,2) )
          vmm  = 0.5*( melt_vel(i-2,k,2) + melt_vel(i-2,k+1,2) )
          gradvp  = 0.
          gradvm  = vm  - vmm
        else
          vp  = 0.5*( melt_vel(i+1,k,2) + melt_vel(i+1,k+1,2) )
          v   = 0.5*( melt_vel(i,k,2) + melt_vel(i,k+1,2) )
          vm  = 0.5*( melt_vel(i-1,k,2) + melt_vel(i-1,k+1,2) )
          vmm = 0.5*( melt_vel(i-2,k,2) + melt_vel(i-2,k+1,2) )
          gradvp  = vp  - v
          gradvm  = vm  - vmm
        endif

        if (fh_x(i,k).gt.0.) fv_x(i,k) = ( vm + 0.5*gradvm )*fh_x(i,k)
        if (fh_x(i,k).lt.0.) fv_x(i,k) = ( v  - 0.5*gradvp )*fh_x(i,k)
 

      endif



   enddo
  enddo  





  do  i = surf_ind(1,1),surf_ind(1,2) 
    do k = surf_ind(2,1),surf_ind(2,2) + 1

 ! computing the fluxes at the boundaries 
 ! fh_y represent hv
 ! fq_y represent (hvv + gravity-term)
 ! fv_y represent x-dir of momentum(hu) projected by v : huv

      if (k.eq.surf_ind(2,1)) then
        h   = melt_height(i,k)
        hm  = h
        v   = 0.5*( melt_vel(i,k,2) + melt_vel(i,k+1,2) )
        vm  = v
        hv  = h  * v
        hvm = hm * vm

        c = max( (abs(vm)+(g*hm)**0.5) , (abs(v)+(g*h)**0.5) )

        fh_y(i,k) = 0.5*(hvm + hv) - 0.5*c*(h-hm)
        fq_y(i,k) = 0.5*(hvm**2 + 0.5*g*hm**2 + hv**2 + 0.5*g*h**2) - 0.5*c*(hv-hvm)

        up  = 0.5*( melt_vel(i,k+1,1) + melt_vel(i+1,k+1,1) )
        u   = 0.5*( melt_vel(i,k  ,1) + melt_vel(i+1,k  ,1) )
        um  = u
        gradum = 0.
        gradup = up - u

        if (fh_y(i,k).gt.0.) fv_y(i,k) = ( um + 0.5*gradum )*fh_y(i,k)
        if (fh_y(i,k).lt.0.) fv_y(i,k) = ( u  - 0.5*gradup )*fh_y(i,k)
      elseif (k.eq.surf_ind(2,2)+1) then
        h   = melt_height(i,k)
        hm  = melt_height(i,k)
        v   = 0.5*( melt_vel(i,k,2) + melt_vel(i,k-1,2) )
        vm  = 0.5*( melt_vel(i,k,2) + melt_vel(i,k-1,2) )
        hv  = h  * v
        hvm = hm * vm

        c = max( (abs(vm)+(g*hm)**0.5) , (abs(v)+(g*h)**0.5) )

        fh_y(i,k) = 0.5*(hvm + hv) - 0.5*c*(h-hm)
        fq_y(i,k) = 0.5*(hvm**2 + 0.5*g*hm**2 + hv**2 + 0.5*g*h**2) - 0.5*c*(hv-hvm)

        u   = 0.5*( melt_vel(i,k  ,1) + melt_vel(i+1,k  ,1) )
        um  = 0.5*( melt_vel(i,k-1,1) + melt_vel(i+1,k-1,1) )
        umm = 0.5*( melt_vel(i,k-2,1) + melt_vel(i+1,k-2,1) )
        gradup  = 0.
        gradum  = um  - umm

        if (fh_y(i,k).gt.0.) fv_y(i,k) = ( um + 0.5*gradum )*fh_y(i,k)
        if (fh_y(i,k).lt.0.) fv_y(i,k) = ( u  - 0.5*gradup )*fh_y(i,k)
      else
 !computing the flues within the domain 
        h   = melt_height(i,k)
        hm  = melt_height(i,k-1)
        v   = 0.5*( melt_vel(i,k,2) + melt_vel(i,k+1,2) )
        vm  = 0.5*( melt_vel(i,k,2) + melt_vel(i,k-1,2) )
        hv  = h  * v
        hvm = hm * vm

        c = max( (abs(vm)+(g*hm)**0.5) , (abs(v)+(g*h)**0.5) )

        fh_y(i,k) = 0.5*(hvm + hv) - 0.5*c*(h-hm)
        fq_y(i,k) = 0.5*(hvm**2 + 0.5*g*hm**2 + hv**2 + 0.5*g*h**2) - 0.5*c*(hv-hvm)

        if (k.eq.surf_ind(2,1)+1) then
          up  = 0.5*( melt_vel(i,k+1,1) + melt_vel(i+1,k+1,1) )
          u   = 0.5*( melt_vel(i,k  ,1) + melt_vel(i+1,k  ,1) )
          um  = 0.5*( melt_vel(i,k-1,1) + melt_vel(i+1,k-1,1) )
          gradup  = up  - u
          gradum  = 0.
        elseif (i.eq.surf_ind(1,2)) then
          u   = 0.5*( melt_vel(i,k  ,1) + melt_vel(i+1,k  ,1) )
          um  = 0.5*( melt_vel(i,k-1,1) + melt_vel(i+1,k-1,1) )
          umm = 0.5*( melt_vel(i,k-2,1) + melt_vel(i+1,k-2,1) )
          gradup  = 0.
          gradum  = um  - umm
        else
          up  = 0.5*( melt_vel(i,k+1,1) + melt_vel(i+1,k+1,1) )
          u   = 0.5*( melt_vel(i,k  ,1) + melt_vel(i+1,k  ,1) )
          um  = 0.5*( melt_vel(i,k-1,1) + melt_vel(i+1,k-1,1) )
          umm = 0.5*( melt_vel(i,k-2,1) + melt_vel(i+1,k-2,1) )
          gradup  = up  - u
          gradum  = um  - umm
        endif

        if (fh_y(i,k).gt.0.) fv_y(i,k) = ( um + 0.5*gradum )*fh_y(i,k)
        if (fh_y(i,k).lt.0.) fv_y(i,k) = ( u  - 0.5*gradup )*fh_y(i,k)


      endif



   enddo
  enddo



 !Discretization of the saint-venant equations and insrting the fluex 
 ! dh/dt + d(hu)/dx + d(hv)/dy = 0
 !updating the values of melt_height 

  do  i = surf_ind(1,1),surf_ind(1,2) 
    do k = surf_ind(2,1),surf_ind(2,2) 
  
       dh = (fh_x(i+1,k) - fh_x(i,k))/surf_dx(1) + (fh_y(i,k+1) - fh_y(i,k))/surf_dx(2)
       melt_height(i,k) = melt_height(i,k) - dt*dh
       melt_height(i,k) = max(melt_height(i,k),0.)
   enddo
  enddo
 ! d(hu)/dt + d(huu+0.5gh^2)/dx + d(huv)/dy = 0
 ! updating the values of melt_vel

  do  i = surf_ind(1,1),surf_ind(1,2) + 1
    do k = surf_ind(2,1),surf_ind(2,2)

       if (i.eq.surf_ind(1,1)) then
         dhu = (fv_x(i,k+1) - fv_x(i,k))/surf_dx(2)
         height_flux(i,k,1) = height_flux(i,k,1) - dt*dhu
         melt_vel(i,k,1)    = height_flux(i,k,1) / ( melt_height(i,k) + small )
       elseif (i.eq.surf_ind(1,2)+1) then 
         dhu = (fv_x(i-1,k+1) - fv_x(i-1,k))/surf_dx(2)
         height_flux(i,k,1) = height_flux(i,k,1) - dt*dhu
         melt_vel(i,k,1)    = height_flux(i,k,1) / ( melt_height(i-1,k) + small )
       else
         fqp = 0.5*(fq_x(i+1,k) + fq_x(i,k)) 
         fqm = 0.5*(fq_x(i-1,k) + fq_x(i,k))
         fvp = 0.5*(fv_x(i,k+1) + fv_x(i-1,k+1))
         fvm = 0.5*(fv_x(i,k  ) + fv_x(i-1,k  ))
         dhu = (fqp-fqm)/surf_dx(1) + (fvp-fvm)/surf_dx(2)
         height_flux(i,k,1) = height_flux(i,k,1) - dt*dhu
         melt_vel(i,k,1) = 2.*height_flux(i,k,1) / (melt_height(i,k)+melt_height(i-1,k)+2.*small)
       endif

    enddo
  enddo

  !updating the value of melt_vel in spanwise direction
  !d(hv)/dt + d(huv)/dx + d(hvv+0.5gh^2)/dy = 0

  do  i = surf_ind(1,1),surf_ind(1,2) 
    do k = surf_ind(2,1),surf_ind(2,2) + 1

       if (k.eq.surf_ind(2,1)) then
         dhv = (fv_y(i+1,k) - fv_y(i,k))/surf_dx(1)
         height_flux(i,k,2) = height_flux(i,k,2) - dt*dhv
         melt_vel(i,k,2)    = height_flux(i,k,2) / ( melt_height(i,k) + small )
       elseif (k.eq.surf_ind(2,2)+1) then
         dhv = (fv_y(i+1,k-1) - fv_y(i,k-1))/surf_dx(1)
         height_flux(i,k,2) = height_flux(i,k,2) - dt*dhu
         melt_vel(i,k,2)    = height_flux(i,k,2) / ( melt_height(i,k-1) + small )
       else
         fqp = 0.5*(fq_y(i,k+1) + fq_y(i,k))
         fqm = 0.5*(fq_y(i,k-1) + fq_y(i,k))
         fvp = 0.5*(fv_y(i+1,k-1) + fv_y(i+1,k))
         fvm = 0.5*(fv_y(i  ,k-1) + fv_y(i  ,k))
         dhv = (fqp-fqm)/surf_dx(2) + (fvp-fvm)/surf_dx(1)
         height_flux(i,k,2) = height_flux(i,k,2) - dt*dhv
         melt_vel(i,k,2) = 2.*height_flux(i,k,2) / (melt_height(i,k)+melt_height(i,k-1)+2.*small)
       endif

    enddo
  enddo


   
  surf_pos = melt_height + melt_pos
 
  

  
  end subroutine increment_SW
 





end module shallow_water_module 
