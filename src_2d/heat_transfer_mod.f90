module heat_transfer_module
  
  ! -----------------------------------------------------------------
  ! This module is used to perform all the calculations relative
  ! to the heat transfer part of the code.
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
                                temp_new, tn_lo, tn_hi , &
                                flxx, fx_lo, fx_hi, &
                                flxy, fy_lo, fy_hi, &
                                idom_old, ido_lo, ido_hi, &
                                idom_new, idn_lo, idn_hi, &
                                geom, dt)

    use heat_flux_module, only: get_boundary_heat_flux
    use material_properties_module, only : get_temp
 
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2) ! bounds of current tile box
    integer, intent(in) :: uo_lo(2), uo_hi(2) ! bounds of input enthalpy box 
    integer, intent(in) :: un_lo(2), un_hi(2) ! bounds of output enthalpy box  
    integer, intent(in) :: to_lo(2), to_hi(2) ! bounds of input temperature box  
    integer, intent(in) :: tn_lo (2), tn_hi (2) ! bounds of output temperature box
    integer, intent(in) :: fx_lo(2), fx_hi(2) ! bounds of the enthalpy flux along x
    integer, intent(in) :: fy_lo(2), fy_hi(2) ! bounds of the enthalpy flux along y
    integer, intent(in) :: ido_lo(2), ido_hi(2) ! bounds of the input idomain box
    integer, intent(in) :: idn_lo(2), idn_hi(2) ! bounds of the output idomain box
    real(amrex_real), intent(in) :: dt ! time step
    real(amrex_real), intent(in) :: time ! time
    real(amrex_real), intent(inout) :: u_old(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2)) ! Input enthalpy 
    real(amrex_real), intent(inout) :: u_new(un_lo(1):un_hi(1),un_lo(2):un_hi(2)) ! Output enthalpy
    real(amrex_real), intent(inout) :: temp_old(to_lo(1):to_hi(1),to_lo(2):to_hi(2)) ! Input temperature
    real(amrex_real), intent(inout) :: temp_new(tn_lo(1):tn_hi(1),tn_lo(2):tn_hi(2)) ! Output temperature
    real(amrex_real), intent(out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2)) ! flux along the x direction  			
    real(amrex_real), intent(out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2)) ! flux along the y direction
    real(amrex_real), intent(in) :: idom_old(ido_lo(1):ido_hi(1),ido_lo(2):ido_hi(2))
    real(amrex_real), intent(in) :: idom_new(idn_lo(1):idn_hi(1),idn_lo(2):idn_hi(2))
    type(amrex_geometry), intent(in) :: geom ! geometry
    
    ! Local variables
    integer :: i,j
    real(amrex_real) :: dx(2) ! Grid size
    real(amrex_real) :: lo_phys(2) ! Physical location of the lowest corner of the tile box
    real(amrex_real) :: qbound(lo(1):hi(1),lo(2):hi(2)) ! Volumetric heating (boundary)

    ! Re-evaluate domain
    do i = lo(1)-1,hi(1)+1
       do  j = lo(2)-1,hi(2)+1

          ! Points added to the domain
          if (nint(idom_old(i,j)).eq.0 .and. nint(idom_new(i,j)).ne.0) then
             u_old(i,j) = u_old(i,j-1)
          ! Points removed from the domain
          else if (nint(idom_new(i,j)).eq.0) then
             u_old(i,j) = 0_amrex_real
          end if
          
       end do
    end do
    
    ! Get grid size
    dx = geom%dx(1:2) ! grid width at level 

    ! Get physical location of the lowest corner of the tile box
    lo_phys = geom%get_physical_location(lo)

                  
    ! Get enthalpy flux 
    call create_face_flux(dx, lo, hi, &
                          u_old, uo_lo, uo_hi, &
                          flxx, fx_lo, fx_hi, &
                          flxy, fy_lo, fy_hi, &
                          temp_old, to_lo, to_hi, &
                          idom_new, idn_lo, idn_hi)
  				  	
    ! Prescribe external heat flux on the free surface
    call get_boundary_heat_flux(time, lo_phys, &
                                dx, lo, hi, &
                                idom_new, idn_lo, idn_hi, &
                                temp_old, to_lo, to_hi, qbound)

    ! Compute enthalpy at the new timestep
    do i = lo(1),hi(1)
       do  j = lo(2),hi(2)
          u_new(i,j) = u_old(i,j) &
               - dt/dx(1) * (flxx(i+1,j) - flxx(i,j)) & ! flux divergence x-direction 
               - dt/dx(2) * (flxy(i,j+1) - flxy(i,j)) & ! flux divergence y-direction 
               + dt*qbound(i,j) ! 'boundary volumetric' source
       end do
    end do
    
  	
    ! Scale the fluxes for the flux registers
    do j = lo(2), hi(2)
       do i = lo(1), hi(1) + 1
          flxx(i,j) = flxx(i,j) * (dt * dx(2))
       end do
    end do

    do j = lo(2), hi(2) + 1
       do i = lo(1), hi(1)
          flxy(i,j) = flxy(i,j) * (dt * dx(1))
       end do
    end do
    
    ! Get temperature corresponding to the output enthalpy
    call get_temp(lo, hi, &
                  u_new, un_lo, un_hi, &
                  temp_new, tn_lo , tn_hi)
    
  end subroutine increment_enthalpy


  ! -----------------------------------------------------------------
  ! Subroutine used to obtain the integer field used to distinguish
  ! between material and background
  ! -----------------------------------------------------------------
  subroutine get_idomain(xlo, dx, lo, hi, &
                         idom, id_lo, id_hi, &
                         temp, t_lo, t_hi)

    use material_properties_module, only : temp_melt   
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: id_lo(2), id_hi(2)
    integer, intent(in) :: t_lo(2), t_hi(2)
    real(amrex_real), intent(in) :: dx(2)    
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1), t_lo(2):t_hi(2))
    real(amrex_real), intent(inout) :: idom(id_lo(1):id_hi(1), id_lo(2):id_hi(2))
    real(amrex_real), intent(in) :: xlo(2)

    ! Local variables
    logical :: find_liquid
    integer :: i,j
    integer :: surf_ind_heat_domain
    real(amrex_real) :: surf_pos_heat_domain(id_lo(1):id_hi(1))

    ! Get location of the free surface
    call get_surf_pos(xlo-dx, dx, id_lo, id_hi, surf_pos_heat_domain)

    ! Check whether the liquid should be distinguished from the solid or not
    if (t_lo(1).eq.lo(1)-1 .and. t_hi(1).eq.hi(1)+1 .and. &
         t_lo(2).eq.lo(2)-1 .and. t_hi(2).eq.hi(2)+1) then
       find_liquid = .true.
    else
       find_liquid = .false.
    end if
    
    ! Set flags to distinguish between material and background
    do i = lo(1)-1, hi(1)+1

       surf_ind_heat_domain = id_lo(2) + &
                              floor((surf_pos_heat_domain(i) - &
                              xlo(2)+dx(2))/dx(2))
       
       do j = lo(2)-1, hi(2)+1

          if (j .le. surf_ind_heat_domain) then

             if (find_liquid) then
                if (temp(i,j).gt.temp_melt) then
                   idom(i,j) = 3
                else if (temp(i,j).eq.temp_melt) then
                   idom(i,j) = 2
                else if (temp(i,j).gt.0) then
                   idom(i,j) = 1
                else
                  STOP 'Temperature is 0 below surf_ind_heat_domain'
                end if
             else
                idom(i,j) = 1
             end if
             
          else
             idom(i,j) = 0
          end if
          
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
                              temp, t_lo, t_hi, &
                              idom, id_lo, id_hi)
  				
    use material_properties_module, only: get_conductivity

    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)  
    integer, intent(in) :: uo_lo(2), uo_hi(2)
    integer, intent(in) :: t_lo(2), t_hi(2)
    integer, intent(in) :: fx_lo(2), fx_hi(2)
    integer, intent(in) :: fy_lo(2), fy_hi(2)
    integer, intent(in) :: id_lo(2), id_hi(2)
    real(amrex_real), intent(in) :: dx(2)    
    real(amrex_real), intent(in) :: u_old(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2))
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2)) 
    real(amrex_real), intent(out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))
    real(amrex_real), intent(out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2))
    real(amrex_real), intent(in) :: temp (t_lo(1):t_hi(1),t_lo(2):t_hi(2))

    ! Local variables
    integer :: i,j
    real(amrex_real) :: ktherm
    real(amrex_real) :: temp_face
    real(amrex_real) :: vx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))

    ! Construct 3D melt velocity profile from the 2D shallow water solution
    call get_face_velocity(lo, hi, &
                           vx, fx_lo, fx_hi, &
                           temp, t_lo, t_hi)
    
    ! Flux along the x direction
    do i = lo(1), hi(1)+1
       do j = lo(2), hi(2)          

          if (nint(idom(i-1,j)).eq.0 .or. nint(idom(i,j)).eq.0) then

             ! Suppress flux at the free surface
             flxx(i,j) = 0_amrex_real

          else
             
             ! Advective component
             if (vx(i,j) > 0_amrex_real) then 
                flxx(i,j)  = u_old(i-1,j)*vx(i,j)
             else 
                flxx(i,j)  = u_old(i,j)*vx(i,j)
             end if
             
             ! Diffusive component
             temp_face = (temp(i,j) + temp(i-1,j))/2_amrex_real
             call get_conductivity(temp_face, ktherm)
             flxx(i,j) = flxx(i,j) - ktherm*(temp(i,j)-temp(i-1,j))/dx(1)
             
          end if
          
       end do
    end do
    
    ! Flux along the y direction
    do i = lo(1), hi(1)
       do j = lo(2), hi(2)+1

          if (nint(idom(i,j-1)).eq.0 .or. nint(idom(i,j)).eq.0) then

             ! Suppress flux at the free surface
             flxy(i,j) = 0_amrex_real

          else
             
             ! Diffusive component (there is no advection in the y direction)
             temp_face = (temp(i,j) + temp(i,j-1))/2_amrex_real
             call get_conductivity(temp_face, ktherm)
             flxy(i,j) = -ktherm*(temp(i,j)-temp(i,j-1))/dx(2)
             
          end if
          
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
                               temp, t_lo, t_hi)

    use material_properties_module, only : temp_melt
    use read_input_module, only : meltvel
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: t_lo(2), t_hi(2)
    integer, intent(in) :: vx_lo(2), vx_hi(2)
    real(amrex_real), intent(in)  :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2))
    real(amrex_real), intent(out) :: vx(vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2))

    ! Local variables
    integer :: i,j
    
    ! THIS MUST BE UPDATED 
    do i = lo(1), hi(1)+1
       do j = lo(2), hi(2) 
          if (temp(i,j).gt.temp_melt) then
             vx(i,j) = meltvel
          else
             vx(i,j) = 0_amrex_real
          end if
       end do   
    end do

    
  end subroutine get_face_velocity
  
  
  ! -----------------------------------------------------------------
  ! Subroutine used to interpolate the free surface position as given
  ! by the fluid solver in order to construct the heat conduction
  ! free interface. Note that the surface position is defined on
  ! the faces of the cells and not on the centers
  ! -----------------------------------------------------------------     
  subroutine get_surf_pos(xlo, dx, lo, hi, surf_pos_heat_domain)

    use amr_data_module, only : surf_ind, surf_pos, surf_xlo, surf_dx  

    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2) 
    real(amrex_real), intent(in) :: xlo(2)
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(out) :: surf_pos_heat_domain(lo(1):hi(1))

    ! Local variables
    integer :: i
    integer :: xind
    real(amrex_real) :: xpos
    
    do  i = lo(1),hi(1)
   
      
      ! In what follows -surf_dx(1)/2 is used since
      ! the surface is staggered 'backwards' on the faces with
      ! respect to the xlo and xpos values which are defined on the
      ! centers. 

       xpos = xlo(1) + (0.5+i-lo(1))*dx(1) 

       ! The nearest integer is taken to round of numerical
       ! errors since we know that dx(1) is n*surf_dx(1) where
       ! n is an integer which depends on the current level and
       ! the refinment ratio between levels.
       xind = nint((xpos - surf_dx(1)/2 - surf_xlo(1))/surf_dx(1))

       if (xind.lt.surf_ind(1,1)) xind = surf_ind(1,1)
       if (xind.ge.surf_ind(1,2)) xind = surf_ind(1,2)-1 
       
       surf_pos_heat_domain(i) = surf_pos(xind)
    end do
    
  end subroutine get_surf_pos
 
 
  ! -----------------------------------------------------------------
  ! Subroutine used to reset the position of the bottom of the melt
  ! pool to the position of the free surface. This routine is
  ! necessary to avoid problems during the re-solidification phase
  ! -----------------------------------------------------------------      
  subroutine reset_melt_pos()
    
    use amr_data_module, only : surf_pos, melt_pos, melt_top, surf_ind
    
    integer :: i
  
    do i =  surf_ind(1,1), surf_ind(1,2) 
       melt_pos(i) = surf_pos(i)
       melt_top(i) = surf_pos(i)
    end do
 
  end subroutine reset_melt_pos


  ! -----------------------------------------------------------------
  ! Subroutine used to get the position of the bottom of the melt
  ! pool
  ! -----------------------------------------------------------------
  subroutine get_melt_pos(lo, hi, idom, id_lo, id_hi, geom)
       
    use amr_data_module, only : melt_pos, melt_top
       
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2) 
    integer, intent(in) :: id_lo(2), id_hi(2)    
    real(amrex_real),     intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2)) 
    type(amrex_geometry), intent(in) :: geom
    
    ! Local variables
    integer :: i,j
    integer :: it(1:2) 
    real(amrex_real) :: grid_pos(1:2)
    
    
    do i = lo(1), hi(1)  ! x-direction
       do j = lo(2), hi(2) 
             
          if (nint(idom(i,j)).eq.3 .and. nint(idom(i,j-1)).ne.3) then
             
             it(1) = i
             it(2) = j
             grid_pos = geom%get_physical_location(it)
             melt_pos(i) = grid_pos(2) 
             
          elseif(nint(idom(i,j)).ne.3 .and. nint(idom(i,j-1)).eq.3) then

            it(1) = i
            it(2) = j
            grid_pos = geom%get_physical_location(it)
            melt_top(i) = grid_pos(2)
            ! if (nint(idom(i,j)).ne.0) write(*,*) &
            !  'WARNING: Melt top not at free surface. Caution in shallow waters solver.'

          end if
          
       end do   
    end do
    
  end subroutine get_melt_pos

  ! -----------------------------------------------------------------
  ! Subroutine used to get the total volume of molten material
  ! -----------------------------------------------------------------
  subroutine integrate_surf(melt_vol)

    use amr_data_module, only : melt_top, melt_pos, surf_ind, surf_dx

    ! Input and output variables
    real(amrex_real), intent(out) :: melt_vol ! Integrated melt area [mm2]

    ! Local variables
    integer :: i 
 
    melt_vol = 0 
    
    do i =  surf_ind(1,1), surf_ind(1,2) 
       melt_vol = melt_vol +  melt_top(i) - melt_pos(i)
    end do
    
    melt_vol = melt_vol*surf_dx(1)*1E6  

  end subroutine integrate_surf

  
end module heat_transfer_module
