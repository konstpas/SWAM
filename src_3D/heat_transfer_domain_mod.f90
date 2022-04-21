module heat_transfer_domain_module
  
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
  public :: update_melt_pos
  public :: get_surf_pos
  public :: get_face_velocity
  public :: reset_melt_pos
  public :: revaluate_heat_domain
  public :: get_local_highest_level
  public :: get_surf_deformation
  
contains
  
  
    ! -----------------------------------------------------------------
    ! Subroutine used to obtain the integer field used to distinguish
    ! between material and background
    ! -----------------------------------------------------------------
    subroutine get_idomain( xlo, dx, lo, hi, &
                            idom, id_lo, id_hi, &
                            temp, t_lo, t_hi)  
 
       use read_input_module, only : geom_name 
       
       ! Input and output variables
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: id_lo(3), id_hi(3)
      integer, intent(in) :: t_lo(3), t_hi(3)
      real(amrex_real), intent(in) :: dx(3)    
      real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1), t_lo(2):t_hi(2), t_lo(3):t_hi(3))
      real(amrex_real), intent(inout) :: idom(id_lo(1):id_hi(1), id_lo(2):id_hi(2), id_lo(3):id_hi(3))
      real(amrex_real), intent(in) :: xlo(3)

       if (geom_name .eq. "Slab") then
          call get_slab_idomain(xlo, dx, lo, hi, &
                               idom, id_lo, id_hi, &
                               temp, t_lo, t_hi)
       elseif (geom_name .eq. "West") then
          call get_west_idomain(xlo, dx, lo, hi, &
                               idom, id_lo, id_hi, &
                               temp, t_lo, t_hi)
       else
         STOP 'Unkown geometry, please select Slab or West'
       end if
    end subroutine get_idomain 
    
    ! -----------------------------------------------------------------
    ! Subroutine used to obtain the integer field used to distinguish
    ! between material and background in slab geometry
    ! -----------------------------------------------------------------
    subroutine get_slab_idomain(xlo, dx, lo, hi, &
                           idom, id_lo, id_hi, &
                           temp, t_lo, t_hi)
  
      use material_properties_module, only : temp_melt    
      
      ! Input and output variables
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: id_lo(3), id_hi(3)
      integer, intent(in) :: t_lo(3), t_hi(3)
      real(amrex_real), intent(in) :: dx(3)    
      real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1), t_lo(2):t_hi(2), t_lo(3):t_hi(3))
      real(amrex_real), intent(inout) :: idom(id_lo(1):id_hi(1), id_lo(2):id_hi(2), id_lo(3):id_hi(3))
      real(amrex_real), intent(in) :: xlo(3)
  
      ! Local variables
      logical :: find_liquid
      integer :: i,j,k
      integer :: surf_ind_heat_domain
      real(amrex_real) :: surf_pos_heat_domain(id_lo(1):id_hi(1),id_lo(3):id_hi(3))
  
      ! Get location of the free surface
      call get_surf_pos(xlo-dx, dx, id_lo, id_hi, surf_pos_heat_domain)
  
      ! Check whether the liquid should be distinguished from the solid or not
      if (t_lo(1).eq.lo(1)-1 .and. t_hi(1).eq.hi(1)+1 .and. &
            t_lo(2).eq.lo(2)-1 .and. t_hi(2).eq.hi(2)+1 .and. &
            t_lo(3).eq.lo(3)-1 .and. t_hi(3).eq.hi(3)+1) then
         find_liquid = .true.
      else
         find_liquid = .false.
      end if
         
      ! Set flags to distinguish between material and background
      do i = lo(1)-1,hi(1)+1
         do k = lo(3)-1,hi(3)+1
            
            surf_ind_heat_domain = id_lo(2) + &
                  floor((surf_pos_heat_domain(i,k) - &
                  xlo(2)+dx(2))/dx(2))
            
            do j = lo(2)-1, hi(2)+1

               if (j .le. surf_ind_heat_domain) then
                  
                  if (find_liquid) then
                     if (temp(i,j,k).gt.temp_melt) then
                        idom(i,j,k) = 3 ! Liquid
                     else if (temp(i,j,k).eq.temp_melt) then
                        idom(i,j,k) = 2 ! Liquid
                     else
                        idom(i,j,k) = 1 ! Solid
                     end if
                  else
                     idom(i,j,k) = 1 ! Liquid or solid (no distinction is made)
                  end if
                  
               else
                  idom(i,j,k) = 0 ! Background
               end if
            
            end do

         end do
      end do
  
    end subroutine get_slab_idomain

    ! -----------------------------------------------------------------
    ! Subroutine used to obtain the integer field used to distinguish
    ! between material and background
    ! -----------------------------------------------------------------
    subroutine get_west_idomain(xlo, dx, lo, hi, &
                                idom, id_lo, id_hi, &
                                temp, t_lo, t_hi)
  
      use material_properties_module, only : temp_melt 
      use read_input_module, only : heat_sample_edge, &
                                    geom_cool_pipe_cntr, &
                                    geom_cool_pipe_radius     
      
      ! Input and output variables
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: id_lo(3), id_hi(3)
      integer, intent(in) :: t_lo(3), t_hi(3)
      real(amrex_real), intent(in) :: dx(3)    
      real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1), t_lo(2):t_hi(2), t_lo(3):t_hi(3))
      real(amrex_real), intent(inout) :: idom(id_lo(1):id_hi(1), id_lo(2):id_hi(2), id_lo(3):id_hi(3))
      real(amrex_real), intent(in) :: xlo(3)
  
      ! Local variables
      logical :: find_liquid
      integer :: i,j,k
      integer :: surf_ind_heat_domain
      real(amrex_real) :: xpos, ypos
      real(amrex_real) :: surf_pos_heat_domain(id_lo(1):id_hi(1),id_lo(3):id_hi(3))
      logical :: pipe_flag
  
      ! Get location of the free surface
      call get_surf_pos(xlo-dx, dx, id_lo, id_hi, surf_pos_heat_domain)
  
      ! Check whether the liquid should be distinguished from the solid or not
      if (t_lo(1).eq.lo(1)-1 .and. t_hi(1).eq.hi(1)+1 .and. &
            t_lo(2).eq.lo(2)-1 .and. t_hi(2).eq.hi(2)+1 .and. &
            t_lo(3).eq.lo(3)-1 .and. t_hi(3).eq.hi(3)+1) then
         find_liquid = .true.
      else
         find_liquid = .false.
      end if
         
      ! Set flags to distinguish between material and background
      do i = lo(1)-1,hi(1)+1
         do k = lo(3)-1,hi(3)+1
            
            surf_ind_heat_domain = id_lo(2) + &
                  floor((surf_pos_heat_domain(i,k) - &
                  xlo(2)+dx(2))/dx(2))

            xpos = xlo(1) + (0.5+i-lo(1))*dx(1) 
            
            do j = lo(2)-1, hi(2)+1
               
               ypos = xlo(2) + (0.5+j-lo(2))*dx(2) 
               if (sqrt((xpos-geom_cool_pipe_cntr(1))**2+(ypos-geom_cool_pipe_cntr(2))**2).lt.geom_cool_pipe_radius) then
                  pipe_flag = .true.
               else
                  pipe_flag = .false.
               end if

               if (j .le. surf_ind_heat_domain .and. xpos .le. heat_sample_edge .and. (.not.pipe_flag)) then
               ! if (j .le. surf_ind_heat_domain .and. (.not.pipe_flag)) then

                  if (find_liquid) then
                     if (temp(i,j,k).gt.temp_melt) then
                        idom(i,j,k) = 3 ! Liquid
                     else if (temp(i,j,k).eq.temp_melt) then
                        idom(i,j,k) = 2 ! Liquid
                     else
                        idom(i,j,k) = 1 ! Solid
                     end if
                  else
                     idom(i,j,k) = 1 ! Liquid or solid (no distinction is made)
                  end if
               elseif (pipe_flag) then
                  idom(i,j,k) = -1   
               else
                  idom(i,j,k) = 0 ! Background
               end if
            
            end do

         end do
      end do
  
    end subroutine get_west_idomain      
    
  ! -----------------------------------------------------------------
  ! Subroutine used to interpolate the free surface position as given
  ! by the fluid solver in order to construct the heat conduction
  ! free interface. Note that the surface position is defined on
  ! the faces of the cells and not on the centers
  ! -----------------------------------------------------------------     
    subroutine get_surf_pos(xlo, dx, lo, hi, surf_pos_heat_domain)
    
      ! use amr_data_module, only : surf_ind, surf_pos, surf_xlo, surf_dx  
      
      ! ! Input and output variables
      ! integer, intent(in) :: lo(3), hi(3) 
      ! real(amrex_real), intent(in) :: xlo(3)
      ! real(amrex_real), intent(in) :: dx(3)
      ! real(amrex_real), intent(out) :: surf_pos_heat_domain(lo(1):hi(1),lo(3):hi(3))
      
      ! ! Local variables
      ! integer :: i, k
      ! integer :: xind, zind
      ! real(amrex_real) :: xpos, zpos
      ! real(amrex_real) :: x(1:2), z(1:2)
      ! real(amrex_real) :: t,u
  
      ! do  i = lo(1),hi(1)
      !    do k = lo(3),hi(3)
            
      !       xpos = xlo(1) + (0.5 + i-lo(1))*dx(1) 
      !       zpos = xlo(3) + (0.5 + k-lo(3))*dx(3)
             
      !       xind = floor(xpos/surf_dx(1))
      !       zind = floor(zpos/surf_dx(2))
      !       if (xind.lt.surf_ind(1,1)) xind = surf_ind(1,1)+1
      !       if (xind.ge.surf_ind(1,2)) xind = surf_ind(1,2) 
      !       if (zind.lt.surf_ind(2,1)) zind = surf_ind(2,1)+1
      !       if (zind.ge.surf_ind(2,2)) zind = surf_ind(2,2) 
  
      !       x(1) = surf_xlo(1)+(xind-0.5)*surf_dx(1)
      !       x(2) = surf_xlo(1)+(xind+0.5)*surf_dx(1)
      !       z(1) = surf_xlo(2)+(zind-0.5)*surf_dx(2)
      !       z(2) = surf_xlo(2)+(zind+0.5)*surf_dx(2)
  
      !       t = min(1.0, max(0.0, (xpos-x(1))/(x(2)-x(1))))
      !       u = min(1.0, max(0.0, (zpos-z(1))/(z(2)-z(1))))

      !       ! Round off numerical accuracies
      !       if (abs(t-0.5).lt.1e-6) t = 0.5
      !       if (abs(u-0.5).lt.1e-6) u = 0.5
      !       if (abs(t-0.0).lt.1e-6) t = 0.0_amrex_real
      !       if (abs(u-0.0).lt.1e-6) u = 0.0_amrex_real
      !       if (abs(t-1.0).lt.1e-6) t = 1
      !       if (abs(u-1.0).lt.1e-6) u = 1
  
      !       surf_pos_heat_domain(i,k) = (1-t)*(1-u)*surf_pos(xind-1,zind-1) + t*(1-u)*surf_pos(xind,zind-1) &
      !                   + t*u*surf_pos(xind,zind) + (1-t)*u*surf_pos(xind-1,zind)
            

      !    end do
      ! end do

      use amr_data_module, only : surf_pos 
    
      ! Input and output variables
      integer, intent(in) :: lo(3), hi(3) 
      real(amrex_real), intent(in) :: xlo(3)
      real(amrex_real), intent(in) :: dx(3)
      real(amrex_real), intent(out) :: surf_pos_heat_domain(lo(1):hi(1),lo(3):hi(3))
      
      ! Local variables
      integer :: i, k
      integer :: xind, zind
      
      do  i = lo(1),hi(1)
         do k = lo(3),hi(3)

            ! Map to maximum level            
            call interp_to_max_lev(lo, xlo, dx, i, k, xind, zind)
       
            ! Position of the free surface in the heat transfer domain at
            ! a given level            
            surf_pos_heat_domain(i,k) = surf_pos(xind, zind)
              
         end do
      end do
      
    end subroutine get_surf_pos
  
    
  ! -----------------------------------------------------------------
  ! Subroutine used to reset the position of the bottom of the melt
  ! pool to the position of the free surface. This routine is
  ! necessary to avoid problems during the re-solidification phase
  ! -----------------------------------------------------------------      
  subroutine reset_melt_pos()
    
    use amr_data_module, only : surf_pos, melt_pos, melt_top, surf_ind
    
    integer :: i,k 
    
    do i =  surf_ind(1,1), surf_ind(1,2) 
       do k = surf_ind(2,1), surf_ind(2,2)
          melt_pos(i,k) = surf_pos(i,k)
          melt_top(i,k) = surf_pos(i,k)
       end do
    end do
    
  end subroutine reset_melt_pos
  
  ! -----------------------------------------------------------------
  ! Subroutine used to update the position of the bottom of the
  ! melt pool
  ! -----------------------------------------------------------------
  subroutine update_melt_pos(lev)
    
   use amr_data_module, only : phi_new,&
                               idomain

   ! Input and output variables
   integer, intent(in) :: lev                            
   ! Local variables
   type(amrex_geometry) :: geom
   type(amrex_mfiter) :: mfi
   type(amrex_box) :: bx
   real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pidom
   
   ! Find the melt position
   if(lev .eq. amrex_max_level) then
   
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
  ! pool
  ! -----------------------------------------------------------------
  subroutine get_melt_pos(lo, hi, idom, id_lo, id_hi, geom)
    
    use amr_data_module, only : melt_pos, melt_top
    
    ! Input and output variables
    integer, intent(in) :: lo(3), hi(3) 
    integer, intent(in) :: id_lo(3), id_hi(3)    
    real(amrex_real),     intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2),id_lo(3):id_hi(3)) 
    type(amrex_geometry), intent(in) :: geom
    
    ! Local variables
    !logical :: check_warning
    integer :: i,j,k
    integer :: it(1:3) 
    real(amrex_real) :: grid_pos(1:3)
    
    !check_warning = .true.
    
    do i = lo(1), hi(1)  ! x-direction
       do k = lo(3), hi(3)  ! z-direction 	
          do j = lo(2), hi(2) 
             
             if (nint(idom(i,j,k)).ge.2 .and. nint(idom(i,j-1,k)).lt.2) then
                
                it(1) = i
                it(2) = j
                it(3) = k 
                grid_pos = geom%get_physical_location(it)
                melt_pos(i,k) = grid_pos(2) 
                
             else if (nint(idom(i,j,k)).lt.2 .and. nint(idom(i,j-1,k)).ge.2) then
                
                it(1) = i
                it(2) = j
                it(3) = k 
                grid_pos = geom%get_physical_location(it)
                melt_top(i,k) = grid_pos(2) 
                ! if (nint(idom(i,j,k)).ne.0 .and. check_warning) then
                !    print *, 'WARNING: Melt top not at free surface. Results from the shallow water solver should not be trusted.'
                !    check_warning = .false.
                ! end if
                
             end if
             
          end do
       end do
    end do
    
  end subroutine get_melt_pos
  
  ! -----------------------------------------------------------------
  ! Subroutine used to re-evaluate the heat equation domain
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
    integer, intent(in) :: lo(3), hi(3) ! bounds of current tile box
    integer, intent(in) :: u_lo(3), u_hi(3) ! bounds of input enthalpy box 
    integer, intent(in) :: ido_lo(3), ido_hi(3) ! bounds of the input idomain box
    integer, intent(in) :: idn_lo(3), idn_hi(3) ! bounds of the output idomain box
    integer, intent(in) :: t_lo(3), t_hi(3) ! bounds of the temperature box
    real(amrex_real), intent(in) :: xlo(3) ! Physical location of box boundaries
    real(amrex_real), intent(in) :: dx(3) ! Grid size
    real(amrex_real), intent(inout) :: u_in(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3)) ! Input enthalpy 
    real(amrex_real), intent(in) :: idom_old(ido_lo(1):ido_hi(1),ido_lo(2):ido_hi(2),ido_lo(3):ido_hi(3))
    real(amrex_real), intent(inout) :: idom_new(idn_lo(1):idn_hi(1),idn_lo(2):idn_hi(2),idn_lo(3):idn_hi(3))
    real(amrex_real), intent(inout) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
       
    !Local variables
    integer :: i,j,k
    integer :: xind, zind
    
    ! Re-evaluate domain
    do i = lo(1)-1,hi(1)+1
       do  j = lo(2)-1,hi(2)+1
          do k = lo(3)-1,hi(3)+1
             
             ! Points added to the domain
             if (nint(idom_old(i,j,k)).eq.0 .and. nint(idom_new(i,j,k)).ne.0) then

               ! Index of free surface element
               call interp_to_max_lev(lo, xlo, dx, i, k, xind, zind)

              ! Get data from upwind column for points added on top of solid
               if (surf_temperature(xind,zind).lt.temp_melt) then 
                  call get_upwind_column(lo, xlo, dx, i, k, xind, zind)
               end if

               ! Update properties
               u_in(i,j,k) = surf_enthalpy(xind,zind)
               temp(i,j,k) = surf_temperature(xind,zind)

               ! Update properties
               u_in(i,j,k) = surf_enthalpy(xind,zind)
               temp(i,j,k) = surf_temperature(xind,zind)
               if (temp(i,j,k).gt.temp_melt) then
                  idom_new(i,j,k) = 3
               else if (temp(i,j,k).eq.temp_melt) then
                  idom_new(i,j,k) = 2
               else
                  ! Add warning here, this should not happen
                  idom_new(i,j,k) = 1
               end if

             ! Points removed from the domain     
             else if (nint(idom_new(i,j,k)).eq.0) then
                u_in(i,j,k) = 0.0_amrex_real
                temp(i,j,k) = 0.0_amrex_real
             end if
             
          end do
       end do
    end do
    
  end subroutine revaluate_heat_domain


  ! -----------------------------------------------------------------
  ! Subroutine used to identify the upwind column while
  ! revaluating the heat domain
  ! -----------------------------------------------------------------  
  subroutine get_upwind_column(lo, xlo, dx, xind_lev, zind_lev, xind, zind)

   use amr_data_module, only : surf_ind, &
                               melt_vel
   
   ! Input and output variables
   integer, intent(in) :: lo(3) ! Bounds of current tile box
   integer, intent(in) :: xind_lev ! Index on the current level in the x direction
   integer, intent(in) :: zind_lev ! Index on the current level in the z direction
   integer, intent(out) :: xind ! Index on the maximum level in the x direction
   integer, intent(out) :: zind ! Index on the maximum level in the z direction
   real(amrex_real), intent(in) :: xlo(3) ! Physical location of box boundaries
   real(amrex_real), intent(in) :: dx(3) ! Grid size

   ! Local variables
   logical found_upwind
   
   ! Map current level to the maximum level
   call interp_to_max_lev(lo, xlo, dx, xind_lev, zind_lev, xind, zind)

   ! Ensure that the upwind colum falls inside the computational domain
   if (xind.lt.surf_ind(1,1)) xind = surf_ind(1,1)
   if (xind.gt.surf_ind(1,2)) xind = surf_ind(1,2)    

    ! Figure out which column is upwind. First try to find upwind in the direction of x, only if that doesn't work
    ! look in the z direction.
    if (melt_vel(xind,zind,1).gt.0_amrex_real) then
      xind = xind-1 
      found_upwind = .true.
    elseif (melt_vel(xind+1,zind,1).lt.0_amrex_real) then
      xind = xind+1
      found_upwind = .true.
    elseif (melt_vel(xind,zind,2).gt.0_amrex_real) then
       zind = zind-1 
       found_upwind = .true.
    elseif (melt_vel(xind,zind+1,2).lt.0_amrex_real) then
      zind = zind+1
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
    if (zind.lt.surf_ind(2,1)) zind = surf_ind(2,1)
    if (zind.gt.surf_ind(2,2)) zind = surf_ind(2,2)   
       
   
 end subroutine get_upwind_column
  
    

  subroutine get_local_highest_level(xlo, dx, lo, hi, lev)

   use read_input_module, only: regrid_dist

   ! Input and output variables
   real(amrex_real), intent(in) :: xlo(3)
   real(amrex_real), intent(in) :: dx(3)
   integer, intent(in) :: lo(3)
   integer, intent(in) :: hi(3)
   integer, intent(out) :: lev(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

   ! Local variables
   real(amrex_real) :: surf_pos_local(lo(1):hi(1), lo(3):hi(3))
   real(amrex_real) :: ypos
   integer :: i,j,k
   integer :: lev_test

   call get_surf_pos(xlo, dx, lo, hi, surf_pos_local)

   do j = lo(2), hi(2)
      ypos = xlo(2)+(j+0.5-lo(2))*dx(2)
      do i = lo(1), hi(1)
         do k = lo(3), hi(3)
            lev_test = 1
            do while(abs(ypos-surf_pos_local(i,k)).lt.regrid_dist(lev_test))
                lev_test = lev_test+1
            end do
            lev(i,j,k) = lev_test-1
         end do
      end do
   end do

  end subroutine get_local_highest_level
  

    ! -----------------------------------------------------------------
    ! Subroutine used to the velocity on the faces of each grid cell.
    ! This subroutine translates to 3D the 2D velocity field obtained
    ! from the solution of the shallow water equations
    ! -----------------------------------------------------------------  
    subroutine get_face_velocity(lo_phys, lo, hi, dx, &
                                 vx, vx_lo, vx_hi, &
                                 vz, vz_lo, vz_hi, &
                                 idom, id_lo, id_hi)
 
     use amr_data_module, only : surf_dx, &
                                 surf_xlo, &
                                 surf_ind, &
                                 melt_vel
     
     ! Input and output variables
     integer, intent(in) :: lo(3), hi(3)
     integer, intent(in) :: id_lo(3), id_hi(3)
     integer, intent(in) :: vx_lo(3), vx_hi(3)
     integer, intent(in) :: vz_lo(3), vz_hi(3)
     real(amrex_real), intent(in) :: dx(3)
     real(amrex_real), intent(in) :: lo_phys(3)
     real(amrex_real), intent(in)  :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2),id_lo(3):id_hi(3))
     real(amrex_real), intent(out) :: vx(vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2),vx_lo(3):vx_hi(3))
     real(amrex_real), intent(out) :: vz(vz_lo(1):vz_hi(1),vz_lo(2):vz_hi(2),vz_lo(3):vz_hi(3))
 
     ! Local variables
     integer :: i,j,k
     integer :: xind, zind
     real(amrex_real) :: xpos, zpos
     
     ! Assign x-component of the velocity 
     do i = lo(1), hi(1)+1
        do j = lo(2), hi(2)
           do k = lo(3), hi(3)
              
              if (nint(idom(i,j,k)).ge.2 .and. nint(idom(i-1,j,k)).ge.2) then
 
                 ! Interpolation similar to the one used to find the free
                 ! surface position in the heat transfer domain. See
                 ! subroutine get_surf_pos in the domain module
                 xpos = lo_phys(1) + (0.5 + i-lo(1))*dx(1) 
                 zpos = lo_phys(3) + (0.5 + k-lo(3))*dx(3)
                 xind = nint((xpos - surf_dx(1)/2 - surf_xlo(1))/surf_dx(1)) 
                 zind = nint((zpos - surf_dx(2)/2 - surf_xlo(2))/surf_dx(2))
                 if (xind.lt.surf_ind(1,1)) xind = surf_ind(1,1)
                 if (xind.gt.surf_ind(1,2)) xind = surf_ind(1,2) 
                 if (zind.lt.surf_ind(2,1)) zind = surf_ind(2,1)
                 if (zind.gt.surf_ind(2,2)) zind = surf_ind(2,2) 
                 vx(i,j,k) = melt_vel(xind,zind,1)
                !  vx(i,j,k) = melt_vel(i,k,1)
                 
              else
                 vx(i,j,k) = 0.0_amrex_real
              end if
              
           end do
        end do   
     end do
 
     ! Assign z-component of the velocity
     do i = lo(1), hi(1)
        do j = lo(2), hi(2)
           do k = lo(3), hi(3)+1
              
              if (nint(idom(i,j,k)).ge.2 .and. nint(idom(i,j,k-1)).ge.2) then
 
                 ! Interpolation similar to the one used to find the free
                 ! surface position in the heat transfer domain. See
                 ! subroutine get_surf_pos in the domain module
                 xpos = lo_phys(1) + (0.5 + i-lo(1))*dx(1) 
                 zpos = lo_phys(3) + (0.5 + k-lo(3))*dx(3)
                 xind = nint((xpos - surf_dx(1)/2 - surf_xlo(1))/surf_dx(1)) 
                 zind = nint((zpos - surf_dx(2)/2 - surf_xlo(2))/surf_dx(2))
                 if (xind.lt.surf_ind(1,1)) xind = surf_ind(1,1)
                 if (xind.gt.surf_ind(1,2)) xind = surf_ind(1,2) 
                 if (zind.lt.surf_ind(2,1)) zind = surf_ind(2,1)
                 if (zind.gt.surf_ind(2,2)) zind = surf_ind(2,2) 
                 vz(i,j,k) = melt_vel(xind,zind,2)
                !  vz(i,j,k) = melt_vel(i,k,2)
                 
              else
                 vz(i,j,k) = 0.0_amrex_real
              end if
 
           end do
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
      integer, intent(in) :: lo(3), hi(3) 
      real(amrex_real), intent(in) :: xlo(3)
      real(amrex_real), intent(in) :: dx(3)
      real(amrex_real), intent(out) :: surf_def_hd(lo(1):hi(1), lo(3):hi(3))
  
      ! Local variables
      integer :: i,k
      integer :: xind, zind
  
      ! Note: xlo is the lowest corner of a given box and is defined on the faces.
      ! xpos refers to the x-coordinate of one cell center of a given box
      ! and is defined on the centers.
      
      do  i = lo(1),hi(1)
         do k = lo(3), hi(3)
            ! Map to maximum level
            call interp_to_max_lev(lo, xlo, dx, i, k, xind, zind)
        
            ! Deformation of the free surface in the heat transfer domain at
            ! a given level
            surf_def_hd(i,k) = surf_deformation(xind,zind)
         end do
      end do
      
    end subroutine get_surf_deformation

  ! -----------------------------------------------------------------
  ! Subroutine used to interpolate spatial locations to the maximum
  ! amrex level where the free surface is defined
  ! -----------------------------------------------------------------  
    subroutine interp_to_max_lev(lo, xlo, dx, xind_lev, zind_lev, xind, zind)

      use amr_data_module, only : surf_ind, &
                                  surf_dx, &
                                  surf_xlo
      
      ! Input and output variables
      integer, intent(in) :: lo(3) ! Bounds of current tile box
      integer, intent(in) :: xind_lev ! Index on the current level in the x-direction
      integer, intent(in) :: zind_lev ! Index on the current level in the z-direction
      integer, intent(out) :: xind ! Index on the maximum level in the x-direction
      integer, intent(out) :: zind ! Index on the maximum level in the z-direction
      real(amrex_real), intent(in) :: xlo(3) ! Physical location of box boundaries
      real(amrex_real), intent(in) :: dx(3) ! Grid size
  
      ! Local variables
      real(amrex_real) :: xpos, zpos
  
      ! Note: xlo is the lowest corner of a given box and is defined on the faces.
      ! xpos refers to the x-coordinate of one cell center of a given box
      ! and is defined on the centers.
      
      ! Position referring to xind_lev and zind_lev
      xpos = xlo(1) + (0.5 + xind_lev - lo(1))*dx(1)
      zpos = xlo(3) + (0.5 + zind_lev - lo(3))*dx(3)
  
      ! Index of xpos at the maximum level
      xind = nint((xpos - surf_dx(1)/2 - surf_xlo(1))/surf_dx(1))
      if (xind.lt.surf_ind(1,1)) xind = surf_ind(1,1)
      if (xind.gt.surf_ind(1,2)) xind = surf_ind(1,2)

      ! Index of zpos at the maximum level
      zind = nint((zpos - surf_dx(2)/2 - surf_xlo(2))/surf_dx(2))
      if (zind.lt.surf_ind(2,1)) zind = surf_ind(2,1)
      if (zind.gt.surf_ind(2,2)) zind = surf_ind(2,2)
      
    end subroutine interp_to_max_lev
end module heat_transfer_domain_module



