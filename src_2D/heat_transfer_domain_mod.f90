module heat_transfer_domain_module
  
  ! -----------------------------------------------------------------
  ! This module is used to manipulate the computational domain used
  ! for the solution of the heat transfer
  ! -----------------------------------------------------------------
  
  use amrex_amr_module
  
  implicit none

  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: get_idomain
  public :: get_melt_pos
  public :: update_melt_pos
  public :: get_surf_pos
  public :: get_face_velocity
  public :: reset_melt_pos
  public :: revaluate_heat_domain
  public :: get_surf_deformation
  
contains
  
  
  ! -----------------------------------------------------------------
  ! Subroutine used to obtain the integer field used to distinguish
  ! between material and background
  ! -----------------------------------------------------------------
  subroutine get_idomain(xlo, dx, lo, hi, &
                         idom, id_lo, id_hi, &
                         temp, t_lo, t_hi)

    use read_input_module, only : geom_name 
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: id_lo(2), id_hi(2)
    integer, intent(in) :: t_lo(2), t_hi(2)
    real(amrex_real), intent(in) :: dx(2)    
    real(amrex_real), intent(inout) :: idom(id_lo(1):id_hi(1), id_lo(2):id_hi(2))
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1), t_lo(2):t_hi(2))
    real(amrex_real), intent(in) :: xlo(2)

    ! Define idomains depending on the geometry specified in input
    if (geom_name .eq. "Slab") then
       
       call get_slab_idomain(xlo, dx, lo, hi, &
                             idom, id_lo, id_hi, &
                             temp, t_lo, t_hi)
       
    elseif (geom_name .eq. "West") then
      
       call get_west_idomain(xlo, dx, lo, hi, &
                             idom, id_lo, id_hi, &
                             temp, t_lo, t_hi)
       
    elseif (geom_name .eq. "West_rectangular") then
       
       call get_west_rectangular_idomain(xlo, dx, lo, hi, &
                                         idom, id_lo, id_hi, &
                                         temp, t_lo, t_hi)

    else

       STOP 'Unknown geometry name'
       
   end if
   
  end subroutine get_idomain

  ! -----------------------------------------------------------------
  ! Subroutine used to obtain the integer field used to distinguish
  ! between material and background in a slab geometry
  ! -----------------------------------------------------------------
  subroutine get_slab_idomain(xlo, dx, lo, hi, &
                              idom, id_lo, id_hi, &
                              temp, t_lo, t_hi)

    use material_properties_module, only : temp_melt    
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: id_lo(2), id_hi(2)
    integer, intent(in) :: t_lo(2), t_hi(2)
    real(amrex_real), intent(in) :: dx(2)    
    real(amrex_real), intent(inout) :: idom(id_lo(1):id_hi(1), id_lo(2):id_hi(2))
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1), t_lo(2):t_hi(2))
    real(amrex_real), intent(in) :: xlo(2)

    ! Local variables
    logical :: find_liquid
    integer :: i,j
    integer :: surf_ind_hd
    real(amrex_real) :: surf_pos_hd(id_lo(1):id_hi(1))
 
    ! Get location of the free surface
    call get_surf_pos(xlo-dx, dx, id_lo, id_hi, surf_pos_hd)

    ! Do not distinguish between liquid and solid if the temperature
    ! multifab passed in input doesn't have ghost point. This should
    ! happen only in the initialization phase
    if (t_lo(1).eq.lo(1)-1 .and. t_hi(1).eq.hi(1)+1 .and. &
        t_lo(2).eq.lo(2)-1 .and. t_hi(2).eq.hi(2)+1) then
       find_liquid = .true.
    else
       find_liquid = .false.
    end if
    
    ! Set flags to distinguish between material and background
    do j = lo(2)-1,hi(2)+1
       do i = lo(1)-1, hi(1)+1

          ! Index of the free surface in the heat transfer domain
          surf_ind_hd = id_lo(2) + &
                                 floor((surf_pos_hd(i) - &
                                 xlo(2)+dx(2))/dx(2))
          
          if (j .le. surf_ind_hd) then
             
             if (find_liquid) then
                if (temp(i,j).gt.temp_melt) then
                   idom(i,j) = 3 ! Liquid
                else if (temp(i,j).eq.temp_melt) then
                   idom(i,j) = 2 ! Material at melting (mushy)
                else
                   idom(i,j) = 1 ! Solid
                end if
             else
                idom(i,j) = 1 ! Liquid or solid (no distinction is made)
             end if
             
          else
             idom(i,j) = 0 ! Background
          end if
          
       end do       
    end do
      

  end subroutine get_slab_idomain

  ! -----------------------------------------------------------------
  ! Subroutine used to obtain the integer field used to distinguish
  ! between material and background in WEST geometry with circular
  ! cooling pipe
  ! -----------------------------------------------------------------
  subroutine get_west_idomain(xlo, dx, lo, hi, &
                              idom, id_lo, id_hi, &
                              temp, t_lo, t_hi)
  
    use material_properties_module, only : temp_melt  
    use read_input_module, only : heat_sample_edge, &
                                  geom_cool_pipe_cntr, &
                                  geom_cool_pipe_radius  
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: id_lo(2), id_hi(2)
    integer, intent(in) :: t_lo(2), t_hi(2)
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(inout) :: idom(id_lo(1):id_hi(1), id_lo(2):id_hi(2))
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1), t_lo(2):t_hi(2))
    real(amrex_real), intent(in) :: xlo(2)
  
    ! Local variables
    logical :: find_liquid
    integer :: i,j
    integer :: surf_ind_hd
    real(amrex_real) :: surf_pos_hd(id_lo(1):id_hi(1))
    real(amrex_real) :: xpos, ypos
    logical :: pipe_flag
    
    ! Get location of the free surface
    call get_surf_pos(xlo-dx, dx, id_lo, id_hi, surf_pos_hd)
    
    ! Do not distinguish between liquid and solid if the temperature
    ! multifab passed in input doesn't have ghost point. This should
    ! happen only in the initialization phase
    if (t_lo(1).eq.lo(1)-1 .and. t_hi(1).eq.hi(1)+1 .and. &
        t_lo(2).eq.lo(2)-1 .and. t_hi(2).eq.hi(2)+1) then
       find_liquid = .true.
    else
       find_liquid = .false.
    end if
    
    ! Set flags to distinguish between material and background
    do j = lo(2)-1, hi(2)+1  
       do i = lo(1)-1, hi(1)+1

          ! Index of the free surface in the heat transfer domain
          surf_ind_hd = id_lo(2) + &
                                 floor((surf_pos_hd(i) - &
                                 xlo(2)+dx(2))/dx(2))
          
          ! x position of the i-th center in the current box
          xpos = xlo(1) + (0.5 + i - lo(1))*dx(1) 
          
       
          ! y position of the j-th center in the current box
          ypos = xlo(2) + (0.5+j-lo(2))*dx(2)

          ! Check if the grid point (i,j) falls inside the cooling pipe
          if (sqrt((xpos-geom_cool_pipe_cntr(1))**2 &
                   +(ypos-geom_cool_pipe_cntr(2))**2) .lt. geom_cool_pipe_radius) then
             pipe_flag = .true.
          else
             pipe_flag = .false.
          end if

          ! Set flags
          if (j .le. surf_ind_hd .and. &
               xpos .le. heat_sample_edge .and. &
               (.not.pipe_flag)) then             
             if (find_liquid) then
                if (temp(i,j).gt.temp_melt) then
                   idom(i,j) = 3 ! Liquid
                else if (temp(i,j).eq.temp_melt) then
                   idom(i,j) = 2 ! Material at melting (mushy)
                else
                   idom(i,j) = 1 ! Solid 
                end if
             else
                idom(i,j) = 1 ! Liquid or solid (no distinction is made)
             end if
            
          elseif (pipe_flag) then
             idom(i,j) = -1 ! Background inside the cooling pipe
          else
             idom(i,j) = 0 ! Background above the free surface
          end if
          
       end do
    end do
    
  end subroutine get_west_idomain


  ! -----------------------------------------------------------------
  ! Subroutine used to obtain the integer field used to distinguish
  ! between material and background in WEST geometry with rectangular
  ! cooling pipe
  ! -----------------------------------------------------------------
  subroutine get_west_rectangular_idomain(xlo, dx, lo, hi, &
                                          idom, id_lo, id_hi, &
                                          temp, t_lo, t_hi)
    
    use material_properties_module, only : temp_melt  
    use read_input_module, only : heat_sample_edge, &
         geom_cool_pipe_cntr, &
         geom_cool_pipe_radius  
    
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
    integer :: surf_ind_hd
    real(amrex_real) :: surf_pos_hd(id_lo(1):id_hi(1))
    real(amrex_real) :: xpos, ypos
    logical :: pipe_flag
    real(amrex_real), parameter :: pi = 3.14159265358979
    
    ! Get location of the free surface
    call get_surf_pos(xlo-dx, dx, id_lo, id_hi, surf_pos_hd)
    
    ! Do not distinguish between liquid and solid if the temperature
    ! multifab passed in input doesn't have ghost point. This should
    ! happen only in the initialization phase
    if (t_lo(1).eq.lo(1)-1 .and. t_hi(1).eq.hi(1)+1 .and. &
         t_lo(2).eq.lo(2)-1 .and. t_hi(2).eq.hi(2)+1) then
       find_liquid = .true.
    else
       find_liquid = .false.
    end if
    
    ! Set flags to distinguish between material and background
    do j = lo(2)-1, hi(2)+1  
       do i = lo(1)-1, hi(1)+1       

          ! Index of the free surface in the heat transfer domain
          surf_ind_hd = id_lo(2) + &
                                 floor((surf_pos_hd(i) - &
                                 xlo(2)+dx(2))/dx(2))

          ! x position of the i-th center in the current box
          xpos = xlo(1) + (0.5+i-lo(1))*dx(1) 

          ! y position of the j-th center in the current box
          ypos = xlo(2) + (0.5+j-lo(2))*dx(2)

          ! Check if the grid point (i,j) falls inside the cooling pipe
          if((xpos.lt.geom_cool_pipe_cntr(1)+(pi)*geom_cool_pipe_radius/4) .and. & 
               (xpos.gt.geom_cool_pipe_cntr(1)-(pi)*geom_cool_pipe_radius/4) .and. &
               (ypos.lt.geom_cool_pipe_cntr(2)+(pi)*geom_cool_pipe_radius/4) .and. &
               (ypos.gt.geom_cool_pipe_cntr(2)-(pi)*geom_cool_pipe_radius/4)) then
             pipe_flag = .true.
          else
             pipe_flag = .false.
          end if

          ! Set flags
          if (j .le. surf_ind_hd .and. &
               xpos .le. heat_sample_edge .and. &
               (.not.pipe_flag)) then
             if (find_liquid) then
                if (temp(i,j).gt.temp_melt) then
                   idom(i,j) = 3 ! Liquid
                else if (temp(i,j).eq.temp_melt) then
                   idom(i,j) = 2 ! Material at melting (mushy)
                else
                   idom(i,j) = 1 ! Solid
                end if
             else
                idom(i,j) = 1 ! Solid or liquid (no distinction is made)
             end if
             
          elseif (pipe_flag) then
             idom(i,j) = -1 ! Background inside the cooling pipe
          else
             idom(i,j) = 0 ! Background above the free surface
          end if
            
       end do       
    end do
    
  end subroutine get_west_rectangular_idomain
    

   
  ! -----------------------------------------------------------------
  ! Subroutine used to interpolate the free surface position as given
  ! by the fluid solver in order to construct the heat conduction
  ! free interface at a given box on a given level characterized
  ! by xlo (lower corner) and dx (grid resolution) 
  ! -----------------------------------------------------------------     
  subroutine get_surf_pos(xlo, dx, lo, hi, surf_pos_hd)

    
    use amr_data_module, only : surf_pos

    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2) 
    real(amrex_real), intent(in) :: xlo(2)
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(out) :: surf_pos_hd(lo(1):hi(1))

    ! Local variables
    integer :: i
    integer :: xind

    ! Note: xlo is the lowest corner of a given box and is defined on the faces.
    ! xpos refers to the x-coordinate of one cell center of a given box
    ! and is defined on the centers.
    
    do  i = lo(1),hi(1)

       ! Map to maximum level
       call interp_to_max_lev(lo, xlo, dx, i, xind)
      
       ! Position of the free surface in the heat transfer domain at
       ! a given level
       surf_pos_hd(i) = surf_pos(xind)
       
    end do
    
  end subroutine get_surf_pos

 
  ! -----------------------------------------------------------------
  ! Subroutine used to reset the position of the bottom of the melt
  ! pool to the position of the free surface. This routine is
  ! necessary to avoid problems during the re-solidification phase
  ! -----------------------------------------------------------------      
  subroutine reset_melt_pos()
    
    use amr_data_module, only : surf_pos, &
                                melt_pos, &
                                melt_top, &
                                surf_ind
    
    integer :: i
    
    do i =  surf_ind(1,1), surf_ind(1,2) 
       melt_pos(i) = surf_pos(i)
       melt_top(i) = surf_pos(i)
    end do
 
  end subroutine reset_melt_pos


  ! -----------------------------------------------------------------
  ! Subroutine used to update the position of the bottom of the
  ! melt pool for a given level
  ! -----------------------------------------------------------------
  subroutine update_melt_pos(lev)
    
    use amr_data_module, only : phi_new,&
                                idomain

    ! Input variables
    integer, intent(in) :: lev
    
    ! Local variables
    type(amrex_geometry) :: geom
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidom
    
    if (lev.eq.amrex_max_level) then
       
       ! Geometry
       geom = amrex_geom(lev)
       
       ! Loop through all the boxes in the level
       !$omp parallel private(mfi, bx, pidom)
       call amrex_mfiter_build(mfi, phi_new(lev), tiling=.false.)
       do while(mfi%next())
          bx = mfi%validbox()
          pidom  => idomain(lev)%dataptr(mfi)
          call get_melt_pos(bx%lo, bx%hi, &
                            pidom, lbound(pidom), ubound(pidom), &
                            geom)
          
       end do
       call amrex_mfiter_destroy(mfi)
       !$omp end parallel   

    end if
    
  end subroutine update_melt_pos

  ! -----------------------------------------------------------------
  ! Subroutine used to get the position of the bottom of the melt
  ! pool for a box of a given level
  ! -----------------------------------------------------------------
  subroutine get_melt_pos(lo, hi, idom, id_lo, id_hi, geom)
       
    use amr_data_module, only : melt_pos, &
                                melt_top
   
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2) 
    integer, intent(in) :: id_lo(2), id_hi(2)    
    real(amrex_real),     intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2)) 
    type(amrex_geometry), intent(in) :: geom
    
    ! Local variables
    integer :: i,j
    integer :: it(1:2) 
    real(amrex_real) :: grid_pos(1:2)

    do j = lo(2), hi(2)        
       do i = lo(1), hi(1)  
          
          ! Positon of the melt bottom
          if (nint(idom(i,j)).ge.2 .and. nint(idom(i,j-1)).lt.2) then
             
             it(1) = i
             it(2) = j
             grid_pos = geom%get_physical_location(it)
             melt_pos(i) = grid_pos(2) 

             
          ! Position of the melt top (is not necessarily equal to the free surface,
          ! due to cooling effects)
          elseif(nint(idom(i,j)).lt.2 .and. nint(idom(i,j-1)).ge.2) then

             it(1) = i
             it(2) = j
             grid_pos = geom%get_physical_location(it)
             melt_top(i) = grid_pos(2)
            
          end if

       end do   
    end do
    
  end subroutine get_melt_pos
  
  ! -----------------------------------------------------------------
  ! Subroutine used to re-evaluate the heat transfer domain
  ! -----------------------------------------------------------------  
  subroutine revaluate_heat_domain(xlo, dx, lo, hi, &
                                   idom_old, ido_lo, ido_hi, &
                                   idom_new, idn_lo, idn_hi, &
                                   u_in, u_lo, u_hi, &
                                   temp, t_lo, t_hi)

    use material_properties_module, only : temp_melt
    use amr_data_module, only : surf_temperature, &
                                surf_enthalpy
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2) ! bounds of current tile box
    integer, intent(in) :: u_lo(2), u_hi(2) ! bounds of input enthalpy box 
    integer, intent(in) :: ido_lo(2), ido_hi(2) ! bounds of the input idomain box
    integer, intent(in) :: idn_lo(2), idn_hi(2) ! bounds of the output idomain box
    integer, intent(in) :: t_lo(2), t_hi(2) ! bounds of the temperature box
    real(amrex_real), intent(in) :: xlo(2) ! Physical location of box boundaries
    real(amrex_real), intent(in) :: dx(2) ! Grid size
    real(amrex_real), intent(inout) :: u_in(u_lo(1):u_hi(1),u_lo(2):u_hi(2)) ! Input enthalpy 
    real(amrex_real), intent(in) :: idom_old(ido_lo(1):ido_hi(1),ido_lo(2):ido_hi(2))
    real(amrex_real), intent(inout) :: idom_new(idn_lo(1):idn_hi(1),idn_lo(2):idn_hi(2))
    real(amrex_real), intent(inout) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2))
     
    ! Local variables
    integer :: i,j
    integer :: xind

    do  j = lo(2)-1,hi(2)+1
       do i = lo(1)-1,hi(1)+1

          ! Points added to the domain
          if (nint(idom_old(i,j)).eq.0 .and. nint(idom_new(i,j)).ne.0) then

             ! Index of free surface element
             call interp_to_max_lev(lo, xlo, dx, i, xind)

             ! Get data from upwind column for points added on top of solid
             if (surf_temperature(xind).lt.temp_melt) then 
                call get_upwind_column_xind(lo, xlo, dx, i, xind)
             end if

             ! Update properties
             u_in(i,j) = surf_enthalpy(xind)
             temp(i,j) = surf_temperature(xind)
             if (temp(i,j).gt.temp_melt) then
                idom_new(i,j) = 3
             else if (temp(i,j).eq.temp_melt) then
                idom_new(i,j) = 2
             else
                ! Add warning here, this should not happen
                idom_new(i,j) = 1
             end if
             
          ! Points removed from the domain
          else if (nint(idom_new(i,j)).eq.0) then
             
             u_in(i,j) = 0.0_amrex_real
             temp(i,j) = 0.0_amrex_real
             
          end if
          
       end do
    end do
    
  end subroutine revaluate_heat_domain

  ! -----------------------------------------------------------------
  ! Subroutine used to identify the upwind column while
  ! revaluating the heat domain
  ! -----------------------------------------------------------------  
  subroutine get_upwind_column_xind(lo, xlo, dx, xind_lev, xind)

    use amr_data_module, only : surf_ind, &
                                melt_vel
    
    ! Input and output variables
    integer, intent(in) :: lo(2) ! Bounds of current tile box
    integer, intent(in) :: xind_lev ! Index on the current level
    integer, intent(out) :: xind ! Index on the maximum level
    real(amrex_real), intent(in) :: xlo(2) ! Physical location of box boundaries
    real(amrex_real), intent(in) :: dx(2) ! Grid size

    ! Local variables
    logical found_upwind
    
    ! Map current level to the maximum level
    call interp_to_max_lev(lo, xlo, dx, xind_lev, xind)
    
    ! Identify upwind column
    if (melt_vel(xind,1).gt.0_amrex_real) then
       xind = xind - 1
       found_upwind = .true.
    else if (melt_vel(xind+1,1).lt.0_amrex_real) then
       xind = xind + 1
       found_upwind = .true.
    end if

    if (.not.found_upwind) then
       xind = xind - 1
       print *, 'Upwind column could not be identified:'&
            ' Taking the column on the left as upwind'
    end if

    ! Ensure that the upwind colum falls inside the computational domain
    if (xind.lt.surf_ind(1,1)) xind = surf_ind(1,1)
    if (xind.gt.surf_ind(1,2)) xind = surf_ind(1,2)                 
    
    
  end subroutine get_upwind_column_xind


  ! -----------------------------------------------------------------
  ! Subroutine used to the velocity on the faces of each grid cell.
  ! This subroutine translates to 2D the 1D velocity field obtained
  ! from the solution of the shallow water equations
  ! -----------------------------------------------------------------  
  subroutine get_face_velocity(xlo, lo, hi, dx, &
                               vx, vx_lo, vx_hi, &
                               idom, id_lo, id_hi)

    use amr_data_module, only : melt_vel
                                
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: id_lo(2), id_hi(2)
    integer, intent(in) :: vx_lo(2), vx_hi(2)
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(in) :: xlo(2)
    real(amrex_real), intent(in)  :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2))
    real(amrex_real), intent(out) :: vx(vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2))

    ! Local variables
    integer :: i,j
    integer :: xind
    
    ! Assign x-component of the velocity
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)+1
          
          if (nint(idom(i,j)).ge.2 .and. nint(idom(i-1,j)).ge.2) then

             call interp_to_max_lev(lo, xlo, dx, i, xind)
             vx(i,j) = melt_vel(xind,1)
             
          else

             vx(i,j) = 0.0_amrex_real
             
          end if
          
       end do
    end do

    
  end subroutine get_face_velocity


  ! -----------------------------------------------------------------
  ! Subroutine used to interpolate the free surface deformation as
  ! given by the fluid solver at a given box on a given level
  ! characterized by xlo (lower corner) and dx (grid resolution) 
  ! -----------------------------------------------------------------     
  subroutine get_surf_deformation(xlo, dx, lo, hi, surf_def_hd)

    use amr_data_module, only : surf_deformation
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2) 
    real(amrex_real), intent(in) :: xlo(2)
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(out) :: surf_def_hd(lo(1):hi(1))

    ! Local variables
    integer :: i
    integer :: xind

    ! Note: xlo is the lowest corner of a given box and is defined on the faces.
    ! xpos refers to the x-coordinate of one cell center of a given box
    ! and is defined on the centers.
    
    do  i = lo(1),hi(1)

       ! Map to maximum level
       call interp_to_max_lev(lo, xlo, dx, i, xind)
      
       ! Deformation of the free surface in the heat transfer domain at
       ! a given level
       surf_def_hd(i) = surf_deformation(xind)
       
    end do
    
  end subroutine get_surf_deformation

  ! -----------------------------------------------------------------
  ! Subroutine used to interpolate spatial locations to the maximum
  ! amrex level where the free surface is defined
  ! -----------------------------------------------------------------  
  subroutine interp_to_max_lev(lo, xlo, dx, xind_lev, xind)

    use amr_data_module, only : surf_ind, &
                                surf_dx, &
                                surf_xlo
    
    ! Input and output variables
    integer, intent(in) :: lo(2) ! Bounds of current tile box
    integer, intent(in) :: xind_lev ! Index on the current level
    integer, intent(out) :: xind ! Index on the maximum level
    real(amrex_real), intent(in) :: xlo(2) ! Physical location of box boundaries
    real(amrex_real), intent(in) :: dx(2) ! Grid size

    ! Local variables
    real(amrex_real) :: xpos

    ! Note: xlo is the lowest corner of a given box and is defined on the faces.
    ! xpos refers to the x-coordinate of one cell center of a given box
    ! and is defined on the centers.
    
    ! Position referring to xind_lev
    xpos = xlo(1) + (0.5 + xind_lev - lo(1))*dx(1)

    ! Index of xpos at the maximum level
    xind = nint((xpos - surf_dx(1)/2 - surf_xlo(1))/surf_dx(1))
    if (xind.lt.surf_ind(1,1)) xind = surf_ind(1,1)
    if (xind.gt.surf_ind(1,2)) xind = surf_ind(1,2)
    
  end subroutine interp_to_max_lev

  
end module heat_transfer_domain_module
