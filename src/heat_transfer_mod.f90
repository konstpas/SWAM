module heat_transfer_module
  
  ! -----------------------------------------------------------------
  ! This module is used to perform all the calculations relative
  ! to the heat transfer part of the code.
  ! NOTE: As of July 24, 2021 there are several parts of this module
  ! that have to be updated. The main changes have to do with:
  ! 1) Implementing the get_idomain routine
  ! 2) Implementing the get_face_velocity routine
  ! 3) Implementing appropriate bounds for the velocities that
  !    do not depend on the enthalpy
  ! -----------------------------------------------------------------
  
  use amrex_amr_module
  
  implicit none

  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: get_idomain
  public :: get_melt_pos
  public :: get_surf_pos
  public :: increment_enthalpy
  public :: integrate_surf
  public :: reset_melt_pos
  
contains

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy at a new time step
  ! -----------------------------------------------------------------
  subroutine increment_enthalpy(time, lo, hi, &
                                u_old,  uo_lo, uo_hi, &
                                u_new, un_lo, un_hi, &
                                temp_old, to_lo, to_hi, &
                                temp_new, tn_lo , tn_hi , &
                                flxx, fx_lo, fx_hi, &
                                flxy, fy_lo, fy_hi, &
                                flxz, fz_lo, fz_hi, &
                                idom_old, ido_lo, ido_hi, &
                                idom_new, idn_lo, idn_hi, &
                                geom, dt)

    use material_properties_module, only : get_temp, get_maxdiffus, get_enthalpy
    use read_input_module, only : tempinit
    
    ! Input and output variables
    integer, intent(in) :: lo(3), hi(3) ! bounds of current tile box
    integer, intent(in) :: uo_lo(3), uo_hi(3) ! bounds of input enthalpy box 
    integer, intent(in) :: un_lo(3), un_hi(3) ! bounds of output enthalpy box  
    integer, intent(in) :: to_lo(3), to_hi(3) ! bounds of input temperature box  
    integer, intent(in) :: tn_lo (3), tn_hi (3) ! bounds of output temperature box
    integer, intent(in) :: fx_lo(3), fx_hi(3) ! bounds of the enthalpy flux along x
    integer, intent(in) :: fy_lo(3), fy_hi(3) ! bounds of the enthalpy flux along y
    integer, intent(in) :: fz_lo(3), fz_hi(3) ! bounds of the enthalpy flux along z
     integer, intent(in) :: ido_lo(3), ido_hi(3) ! bounds of the input idomain box
    integer, intent(in) :: idn_lo(3), idn_hi(3) ! bounds of the output idomain box
    real(amrex_real), intent(in) :: dt ! time step
    real(amrex_real), intent(in) :: time ! time
    real(amrex_real), intent(inout) :: u_old (uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3)) ! Input enthalpy 
    real(amrex_real), intent(inout) :: u_new(un_lo(1):un_hi(1),un_lo(2):un_hi(2),un_lo(3):un_hi(3)) ! Output enthalpy
    real(amrex_real), intent(inout) :: temp_old(to_lo(1):to_hi(1),to_lo(2):to_hi(2),to_lo(3):to_hi(3)) ! Input temperature
    real(amrex_real), intent(inout) :: temp_new(tn_lo(1):tn_hi(1),tn_lo(2):tn_hi(2),tn_lo(3):tn_hi(3)) ! Output temperature
    real(amrex_real), intent(out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3)) ! flux along the x direction  			
    real(amrex_real), intent(out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3)) ! flux along the y direction
    real(amrex_real), intent(out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3)) ! flux along the z direction
    real(amrex_real), intent(in) :: idom_old(ido_lo(1):ido_hi(1),ido_lo(2):ido_hi(2),ido_lo(3):ido_hi(3))
    real(amrex_real), intent(in) :: idom_new(idn_lo(1):idn_hi(1),idn_lo(2):idn_hi(2),idn_lo(3):idn_hi(3))
    type(amrex_geometry), intent(in) :: geom ! geometry
    
    !Local variables
    integer :: i,j,k
    real(amrex_real) :: dx(3) ! Grid size
    real(amrex_real) :: lo_phys(3) ! Physical location of the lowest corner of the tile box
    real(amrex_real) :: qbound(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) ! Volumetric heating (boundary)
    real(amrex_real) :: u_back

    ! Enthalpy of the background (this should be computed once and stored)
    call get_enthalpy(tempinit,u_back)

    ! Re-evaluate domain
    do i = lo(1)-1,hi(1)+1
       do  j = lo(2)-1,hi(2)+1
          do k = lo(3)-1,hi(3)+1
             
             ! Points added to the domain
             if (nint(idom_old(i,j,k)).eq.0 .and. nint(idom_new(i,j,k)).eq.1) then
                u_old(i,j,k) = u_old(i,j-1,k)
             ! Points removed from the domain
             else if (nint(idom_new(i,j,k)).eq.0) then
                u_old(i,j,k) = u_back
             end if
             
          end do
       end do
    end do

    ! Get grid size
    dx = geom%dx(1:3) ! grid width at level 

    ! Get physical location of the lowest corner of the tile box
    lo_phys = geom%get_physical_location(lo)
    
    ! Get temperature corresponding to the input enthalpy
    call get_temp(to_lo, to_hi, &
                  u_old, uo_lo, uo_hi, &
                  temp_old, to_lo, to_hi)
   
    ! Get enthalpy flux 
    call create_face_flux(dx, lo, hi, &
                          u_old, uo_lo, uo_hi, &
                          flxx, fx_lo, fx_hi, &
                          flxy, fy_lo, fy_hi, &
                          flxz, fz_lo, fz_hi, &
                          temp_old, to_lo, to_hi, &
                          idom_new, idn_lo, idn_hi)
  				  	
    ! Prescribe external heat flux on the free surface
    call get_bound_heat(time, lo_phys, &
                        dx, lo, hi, &
                        idom_new, idn_lo, idn_hi, &
                        qbound)

    ! Compute output enthalpy, i.e. compute enthalpy at the next timestep
    do   i = lo(1),hi(1)
       do  j = lo(2),hi(2) 
          do k = lo(3),hi(3)
             u_new(i,j,k) = u_old(i,j,k) &
                  - dt/dx(1) * (flxx(i+1,j,k) - flxx(i,j,k)) &	! flux divergence x-direction 
                  - dt/dx(2) * (flxy(i,j+1,k) - flxy(i,j,k)) &	! flux divergence y-direction 
                  - dt/dx(3) * (flxz(i,j,k+1) - flxz(i,j,k)) &	! flux divergence z-direction
                  + dt*qbound(i,j,k) ! 'boundary volumetric' source
             
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
    
    ! Get temperature corresponding to the output enthalpy
    call get_temp(lo, hi, &
                  u_new, un_lo, un_hi, &
                  temp_new, tn_lo , tn_hi) 

    ! THIS MUST BE UPDATED
    ! find maximum diffusivity for time step determination 
    ! Not called noe, constant max possible diffusivity used.	      
    !call get_maxdiffus(lo, hi, & 
    !		        temp_new, tn_lo, tn_hi)  	
    
  end subroutine increment_enthalpy


  ! -----------------------------------------------------------------
  ! Subroutine used to obtain the integer field used to distinguish
  ! between material and background
  ! -----------------------------------------------------------------
  subroutine get_idomain(xlo, dx, lo, hi, &
                         idom, id_lo, id_hi)

    ! Input and output variables
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: id_lo(3), id_hi(3)
    real(amrex_real), intent(in) :: dx(3)    
    real(amrex_real), intent(inout) :: idom(id_lo(1):id_hi(1), id_lo(2):id_hi(2),id_lo(3):id_hi(3))
    real(amrex_real), intent(in) :: xlo(3)
    
    ! Local variables
    integer :: i,j,k
    integer :: surf_ind_heat_domain
    real(amrex_real) :: surf_pos_heat_domain(id_lo(1):id_hi(1),id_lo(3):id_hi(3))

    ! Get location of the free surface
    call get_surf_pos(xlo-dx, dx, id_lo, id_hi, surf_pos_heat_domain)

    ! Set flags to distinguish between material and background
    do i = lo(1)-1,hi(1)+1
       do k = lo(3)-1,hi(3)+1
          
          surf_ind_heat_domain = id_lo(2) + &
               floor((surf_pos_heat_domain(i,k) - &
               xlo(2)+dx(2))/dx(2))
          
          do j = lo(2)-1,hi(2)+1
             if (j .le. surf_ind_heat_domain) then
                idom(i,j,k) = 1
             else
                idom(i,j,k) = 0
             end if
          end do
          
       end do
    end do
    
  end subroutine get_idomain

  ! -----------------------------------------------------------------
  ! Subroutine used to the enthalpy fluxes on the edges of the grid
  ! -----------------------------------------------------------------  
  subroutine create_face_flux(dx, lo, hi, &
                              u_old, uo_lo, uo_hi, &
                              flxx, fx_lo, fx_hi, &
                              flxy, fy_lo, fy_hi, &
                              flxz, fz_lo, fz_hi, &
                              temp, t_lo, t_hi, &
                              idom, id_lo, id_hi)
  				
    use material_properties_module, only: get_ktherm

    ! Input and output variables
    integer, intent(in) :: lo(3), hi(3)  
    integer, intent(in) :: uo_lo(3), uo_hi(3)
    integer, intent(in) :: t_lo(3), t_hi(3)
    integer, intent(in) :: fx_lo(3), fx_hi(3)
    integer, intent(in) :: fy_lo(3), fy_hi(3)
    integer, intent(in) :: fz_lo(3), fz_hi(3)
    integer, intent(in) :: id_lo(3), id_hi(3)
    real(amrex_real), intent(in) :: dx(3)    
    real(amrex_real), intent(in) :: u_old(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3))
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2),id_lo(3):id_hi(3))
    real(amrex_real), intent(out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3))
    real(amrex_real), intent(out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3))
    real(amrex_real), intent(out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3))
    real(amrex_real), intent(in) :: temp (t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))

    ! Local variables
    integer :: i,j,k 
    real(amrex_real) :: ktherm
    real(amrex_real) :: temp_face
    real(amrex_real) :: vx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3))
    real(amrex_real) :: vz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3))

    ! Construct 3D melt velocity profile from the 2D shallow water solution
    call get_face_velocity(lo, hi, &
                           vx, fx_lo, fx_hi, &
                           vz, fz_lo, fz_hi, &
                           temp, t_lo, t_hi)
    
    ! Flux along the x direction
    do i = lo(1), hi(1)+1
       do j = lo(2), hi(2)
          do k = lo(3), hi(3)

             if (nint(idom(i-1,j,k)).eq.0 .or. nint(idom(i,j,k)).eq.0) then

                ! Suppress flux at the free surface
                flxx(i,j,k) = 0_amrex_real

             else
                
                ! Advective component
                if (vx(i,j,k) > 0_amrex_real) then 
                   flxx(i,j,k)  = u_old(i-1,j,k)*vx(i,j,k)
                else 
                   flxx(i,j,k)  = u_old(i,j,k)*vx(i,j,k)
                end if
                
                ! Diffusive component
                temp_face = (temp(i,j,k) + temp(i-1,j,k))/2_amrex_real
                call get_ktherm(temp_face, ktherm)
                flxx(i,j,k) = flxx(i,j,k) - ktherm*(temp(i,j,k)-temp(i-1,j,k))/dx(1)
                                
                end if
                
          end do
       end do
    end do

    ! Flux along the y direction
    do i = lo(1), hi(1)
       do j = lo(2), hi(2)+1
          do k = lo(3), hi(3)

             if (nint(idom(i,j-1,k)).eq.0 .or. nint(idom(i,j,k)).eq.0) then

                ! Suppress flux at the free surface
                flxy(i,j,k) = 0_amrex_real

             else
                
                ! Diffusive component (there is no advection in the y direction)
                temp_face = (temp(i,j,k) + temp(i,j-1,k))/2_amrex_real
                call get_ktherm(temp_face, ktherm)
                flxy(i,j,k) = -ktherm*(temp(i,j,k)-temp(i,j-1,k))/dx(2)
                               
             end if
             
          end do
       end do
    end do

    ! Flux along the z direction
    do i = lo(1), hi(1)
       do j = lo(2), hi(2)
          do k = lo(3), hi(3)+1

             if (nint(idom(i,j,k-1)).eq.0 .or. nint(idom(i,j,k)).eq.0) then

                ! Suppress flux at the free surface
                flxz(i,j,k) = 0_amrex_real

             else
                
                ! Advective component
                if (vz(i,j,k) > 0_amrex_real) then 
                   flxz(i,j,k)  = u_old(i,j,k-1)*vz(i,j,k)
                else 
                   flxz(i,j,k)  = u_old(i,j,k)*vz(i,j,k)
                end if
                
                ! Diffusive component
                temp_face = (temp(i,j,k) + temp(i,j,k-1))/2_amrex_real
                call get_ktherm(temp_face, ktherm)
                flxz(i,j,k) = flxz(i,j,k) - ktherm*(temp(i,j,k)-temp(i,j,k-1))/dx(3)
                
             end if
             
          end do
       end do
    end do
    
    
  end subroutine create_face_flux
  
  
  ! -----------------------------------------------------------------
  ! Subroutine used to the velocity on the faces of each grid cell.
  ! This subroutine translates to 3D the 2D velocity field obtained
  ! from the solution of the shallow water equations
  ! -----------------------------------------------------------------  
  subroutine get_face_velocity(lo, hi, &
                               vx, vx_lo, vx_hi, &
                               vz, vz_lo, vz_hi, &
                               temp, t_lo, t_hi)

    use material_properties_module, only : melt_point
    use read_input_module, only : meltvel
    
    ! Input and output variables
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: t_lo(3), t_hi(3)
    integer, intent(in) :: vx_lo(3), vx_hi(3)
    integer, intent(in) :: vz_lo(3), vz_hi(3)
    real(amrex_real), intent(in)  :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(amrex_real), intent(out) :: vx(vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2),vx_lo(3):vx_hi(3))
    real(amrex_real), intent(out) :: vz(vz_lo(1):vz_hi(1),vz_lo(2):vz_hi(2),vz_lo(3):vz_hi(3))

    ! Local variables
    integer :: i,j,k
    
    ! THIS MUST BE UPDATED 
    do i = lo(1), hi(1)+1
       do j = lo(2), hi(2)
          do k = lo(3), hi(3)
             if (temp(i,j,k).gt.melt_point) then
                vx(i,j,k) = meltvel
             else
                vx(i,j,k) = 0_amrex_real
             end if
          end do
       end do   
    end do

    do i = lo(1), hi(1)
       do j = lo(2), hi(2)
          do k = lo(3), hi(3)+1
             if (temp(i,j,k).gt.melt_point) then
                vz(i,j,k) = meltvel
             else
                vz(i,j,k) = 0_amrex_real
             end if
          end do
       end do   
    end do
    
  end subroutine get_face_velocity
    

  ! ! -----------------------------------------------------------------
  ! ! Subroutine used to get the flags to supress the enthalpy fluxes
  ! ! at the free surface
  ! ! -----------------------------------------------------------------
  ! subroutine surface_tag(xlo, dx, lo, hi, &
  !                        flxx_flag, fx_lo, fx_hi, &
  !                        flxy_flag, fy_lo, fy_hi, & 
  !                        flxz_flag, fz_lo, fz_hi)
 
  !   ! Input and output variables
  !   integer, intent(in) :: lo(3), hi(3)
  !   integer, intent(in) :: fx_lo(3), fx_hi(3)
  !   integer, intent(in) :: fy_lo(3), fy_hi(3)
  !   integer, intent(in) :: fz_lo(3), fz_hi(3)
  !   logical, intent(out) :: flxx_flag(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3))
  !   logical, intent(out) :: flxy_flag(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3))
  !   logical, intent(out) :: flxz_flag(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3))
  !   real(amrex_real), intent(in) :: xlo(3)
  !   real(amrex_real), intent(in) :: dx(3)
    
  !   ! Local variables
  !   integer :: i,j,k
  !   integer :: surf_lo(3)
  !   integer :: surf_hi(3)
  !   integer :: surf_ind_heat_domain(lo(1)-1:hi(1)+1,lo(3)-1:hi(3)+1)
  !   real(amrex_real) :: surf_pos_heat_domain(lo(1)-1:hi(1)+1,lo(3)-1:hi(3)+1)
  !   real(amrex_real) :: surf_xlo(3)
  
  !   ! Get location of the free surface in the heat solver domain
  !   surf_lo = lo - 1
  !   surf_hi = hi + 1
  !   surf_xlo = xlo - dx
  !   call get_surf_pos(surf_xlo, dx, surf_lo, surf_hi, surf_pos_heat_domain)
    
  !   ! Initialize flags
  !   flxx_flag = .false. 
  !   flxy_flag = .false. 
  !   flxz_flag = .false. 
    
  !   ! Get indexes corresponding to the free surface position
  !   do i = lo(1)-1, hi(1)+1
  !      do k = lo(3)-1, hi(3)+1 
  !         surf_ind_heat_domain(i,k)  =  surf_lo(2) + &
  !                                       floor((surf_pos_heat_domain(i,k) - &
  !                                              surf_xlo(2))/dx(2))  
  !      end do
  !   end do
    
  !   ! Set flag to suppress flux in the x direction across the free surface
  !   do i = lo(1),hi(1)+1 
  !      do k = lo(3),hi(3)

    
  !         if (surf_ind_heat_domain(i,k) .gt. surf_ind_heat_domain(i-1,k)) then

  !            do j = max(surf_ind_heat_domain(i-1,k),lo(2)), & ! max, min since interface may be outside box 
  !                   min(surf_ind_heat_domain(i,k)-1,hi(2)) ! surf_ind_heat_domain(i)-1 because surf_ind_heat_domain(i) edge is outside domain
  !               flxx_flag(i,j,k) = .true.
  !            end do
             
  !         elseif (surf_ind_heat_domain(i,k) .lt. surf_ind_heat_domain(i-1,k)) then
        
  !            do j = max(surf_ind_heat_domain(i,k),lo(2)), & ! max, min since interface may be outside box 
  !                   min(surf_ind_heat_domain(i-1,k)-1,hi(2)) ! -1 because surf_ind_heat_domain(i-1) edge is outside domain
  !               flxx_flag(i,j,k) = .true.
  !            end do
             
  !         end if
          
  !      end do
  !   end do
    
  !   ! Set flag to suppress flux in the y direction across the free surface
  !   do i = lo(1), hi(1)
  !      do k = lo(3), hi(3)
          
  !         if (surf_ind_heat_domain(i,k) >= lo(2) .and. &
  !             surf_ind_heat_domain(i,k) <= hi(2)+1) then
  !            flxy_flag(i,surf_ind_heat_domain(i,k),k) = .true.
  !         end if
            
  !      end do
  !   end do
    
  !   ! Set flag to suppress flux in the x direction across the free surface
  !   do i = lo(1),hi(1)
  !      do k = lo(3),hi(3)+1
          
  !         if (surf_ind_heat_domain(i,k) .gt. surf_ind_heat_domain(i,k-1)) then
             
  !            do j = max(surf_ind_heat_domain(i,k-1),lo(2)), & 
  !                   min(surf_ind_heat_domain(i,k)-1,hi(2)) 
  !               flxz_flag(i,j,k) = .true. 
  !            end do
             
  !         elseif(surf_ind_heat_domain(i,k) .lt. surf_ind_heat_domain(i,k-1)) then
                
  !            do j = max(surf_ind_heat_domain(i,k),lo(2)), & 
  !                   min(surf_ind_heat_domain(i,k-1)-1,hi(2)) 
  !               flxz_flag(i,j,k) = .true. 
  !            end do
  
  !         end if
             
  !      end do
  !   end do

    
  ! end subroutine surface_tag

  
  ! -----------------------------------------------------------------
  ! Subroutine used to interpolate the free surface position as given
  ! by the fluid solver in order to construct the heat conduction
  ! free interface
  ! -----------------------------------------------------------------     
  subroutine get_surf_pos(xlo, dx, lo, hi, surf_pos_heat_domain)

    use amr_data_module, only : surf_ind, surf_pos, surf_xlo, surf_dx  

    ! Input and output variables
    integer, intent(in) :: lo(3), hi(3) 
    real(amrex_real), intent(in) :: xlo(3)
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(out) :: surf_pos_heat_domain(lo(1):hi(1),lo(3):hi(3))

    ! Local variables
    integer :: i, k
    integer :: xind, zind
    real(amrex_real) :: xpos, zpos
    real(amrex_real) :: x_alpha, z_alpha
    real(amrex_real) :: valzp
    real(amrex_real) :: valzm
    
    do  i = lo(1),hi(1)
       do k = lo(3),hi(3)
   
          xpos = xlo(1) + (0.5 + i-lo(1))*dx(1) 
          zpos = xlo(3) + (0.5 + k-lo(3))*dx(3)

          ! In what follows -surf_dx(1)/2  and ceiling are used since
          ! staggered 'backwards' on faces w.r.t values which are centered.  
          xind = ceiling((xpos - surf_dx(1)/2 - surf_xlo(1))/surf_dx(1)) 
          x_alpha = mod(xpos - surf_dx(1)/2 - surf_xlo(1), surf_dx(1))
          zind = ceiling((zpos - surf_dx(2)/2 - surf_xlo(2))/surf_dx(2))
          z_alpha = mod(zpos - surf_dx(2)/2 - surf_xlo(2), surf_dx(2))
          
          if (xind.lt.surf_ind(1,1)) xind = surf_ind(1,1)
          if (xind.ge.surf_ind(1,2)) xind = surf_ind(1,2)-1 
          if (zind.lt.surf_ind(2,1)) zind = surf_ind(2,1)
          if (zind.ge.surf_ind(2,2)) zind = surf_ind(2,2)-1 

          ! interpolated value at zind
          valzm = surf_pos(xind,zind  ) + &
               x_alpha * (surf_pos(xind+1,zind)-surf_pos(xind,zind))
          ! interpolated value at zind+1 
          valzp = surf_pos(xind,zind+1) + &
               x_alpha * (surf_pos(xind+1,zind+1)-surf_pos(xind,zind+1)) 

          ! 2D interpolation
          surf_pos_heat_domain(i,k) = valzm + z_alpha*(valzp - valzm)
          
       end do
    end do
    
  end subroutine get_surf_pos
 
 
  ! -----------------------------------------------------------------
  ! Subroutine used to reset the position of the bottom of the melt
  ! pool to the position of the free surface. This routine is
  ! necessary to avoid problems during the re-solidification phase
  ! -----------------------------------------------------------------      
  subroutine reset_melt_pos()
    
    use amr_data_module, only : surf_pos, melt_pos, surf_ind
    
    integer :: i,k 
  
    do i =  surf_ind(1,1), surf_ind(1,2) 
       do k = surf_ind(2,1), surf_ind(2,2)
          melt_pos(i,k) = surf_pos(i,k)
       end do
    end do
 
  end subroutine reset_melt_pos


  ! -----------------------------------------------------------------
  ! Subroutine used to get the position of the bottom of the melt
  ! pool
  ! -----------------------------------------------------------------
  subroutine get_melt_pos(lo, hi, temp, t_lo, t_hi, geom)
       
    use amr_data_module, only : melt_pos
    use material_properties_module, only : melt_point
       
    ! Input and output variables
    integer, intent(in) :: lo(3), hi(3) 
    integer, intent(in) :: t_lo(3), t_hi(3)    
    real(amrex_real),     intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) 
    type(amrex_geometry), intent(in) :: geom
    
    ! Local variables
    integer :: i,j,k
    integer :: it(1:3) 
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

  
  ! -----------------------------------------------------------------
  ! Subroutine used to prescribe the boundary heating on the free
  ! surface. Note that the incoming heat flux is assigned to the
  ! first node under the free surface in the form of a volumetric
  ! heating (after an appropriate dimensionality correction)
  ! -----------------------------------------------------------------   
  subroutine get_bound_heat(time, xlo, &
                            dx, lo, hi, &
                            idom, id_lo, id_hi, &
                            qb) 

    use read_input_module, only : flux_peak, flux_width, flux_pos, exp_time
    
    ! Input and output variables
    integer, intent(in) :: lo(3), hi(3)  
    integer, intent(in) :: id_lo(3), id_hi(3)
    real(amrex_real), intent(in) :: time
    real(amrex_real), intent(in) :: xlo(3)
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(out) :: qb(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2),id_lo(3):id_hi(3))

    ! Local variables
    real(amrex_real) :: xpos, zpos, ypos
    integer :: i,j,k
    integer :: surf_ind_heat_domain
    real(amrex_real) :: surf_pos_heat_domain(lo(1):hi(1),lo(3):hi(3))

    ! Get location of the free surface
    call get_surf_pos(xlo, dx, lo, hi, surf_pos_heat_domain)

    ! Initialize the heat flux
    qb = 0. 
    
    ! Prescribe boundary heating
    do   i = lo(1), hi(1) 
       do  j = lo(2), hi(2)
          do k = lo(3), hi(3)

             if (time.lt.exp_time) then

                ypos = xlo(2) + (j-lo(2))*dx(2)

                !if (ypos .lt. surf_pos_heat_domain(i,k) .and. ypos .ge. surf_pos_heat_domain(i,k)-dx(2)) then
                !if(ypos .le. surf_pos_heat_domain(i,k) .and. ypos .gt. surf_pos_heat_domain(i,k)-dx(2)) then
                if(nint(idom(i,j,k)).eq.1 .and. nint(idom(i,j+1,k)).eq.0) then

                   xpos = xlo(1) + (i-lo(1))*dx(1)
                   zpos = xlo(3) + (k-lo(3))*dx(3)
                   
                   ! Note: the term /dx(2) converts a surface heat flux [W/m^2]
                   ! into a volumetric heat flux [W/m^3]
                   qb(i,j,k) = flux_peak &
                               *EXP(-((xpos-flux_pos(1))**2)/(flux_width(1)**2) &
                                   -((zpos-flux_pos(2))**2)/(flux_width(2)**2)) &
                                   /dx(2)

                   
                end if
                
             end if
             
          end do
       end do
    end do

       
  end subroutine get_bound_heat
  

  ! -----------------------------------------------------------------
  ! Subroutine used to get the total volume of molten material
  ! -----------------------------------------------------------------
  subroutine integrate_surf(melt_vol)

    use amr_data_module, only : surf_pos, melt_pos, surf_ind, surf_dx

    ! Input and output variables
    real(amrex_real), intent(out) :: melt_vol ! Integrated melt volume [mm3]

    ! Local variables
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
