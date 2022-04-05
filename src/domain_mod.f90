module domain_module
  
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
  public :: integrate_surf
  public :: reset_melt_pos
  public :: revaluate_heat_domain
  public :: get_local_highest_level
  
contains
  
  
    ! -----------------------------------------------------------------
    ! Subroutine used to obtain the integer field used to distinguish
    ! between material and background
    ! -----------------------------------------------------------------
    subroutine get_idomain( xlo, dx, lo, hi, &
                            idom, id_lo, id_hi, &
                            temp, t_lo, t_hi)  
 
       use read_input_module, only : geometry_name 
       
       ! Input and output variables
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: id_lo(3), id_hi(3)
      integer, intent(in) :: t_lo(3), t_hi(3)
      real(amrex_real), intent(in) :: dx(3)    
      real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1), t_lo(2):t_hi(2), t_lo(3):t_hi(3))
      real(amrex_real), intent(inout) :: idom(id_lo(1):id_hi(1), id_lo(2):id_hi(2), id_lo(3):id_hi(3))
      real(amrex_real), intent(in) :: xlo(3)

       if (geometry_name .eq. "Slab") then
          call get_slab_idomain(xlo, dx, lo, hi, &
                               idom, id_lo, id_hi, &
                               temp, t_lo, t_hi)
       elseif (geometry_name .eq. "West") then
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
      use read_input_module, only : sample_edge, &
                                    cool_pipe_cntr, &
                                    cool_pipe_radius     
      
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
               if (sqrt((xpos-cool_pipe_cntr(1))**2+(ypos-cool_pipe_cntr(2))**2).lt.cool_pipe_radius) then
                  pipe_flag = .true.
               else
                  pipe_flag = .false.
               end if

               if (j .le. surf_ind_heat_domain .and. xpos .le. sample_edge .and. (.not.pipe_flag)) then
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
      
      do  i = lo(1),hi(1)
         do k = lo(3),hi(3)
            
            xpos = xlo(1) + (0.5 + i-lo(1))*dx(1) 
            zpos = xlo(3) + (0.5 + k-lo(3))*dx(3)
            
            ! The nearest integer is taken to round of numerical
            ! errors since we know that dx(1) is n*surf_dx(1) where
            ! n is an integer which depends on the current level and
            ! the refinment ratio between levels. 
            xind = nint((xpos - surf_dx(1)/2 - surf_xlo(1))/surf_dx(1)) 
            zind = nint((zpos - surf_dx(2)/2 - surf_xlo(2))/surf_dx(2))
            
            if (xind.lt.surf_ind(1,1)) xind = surf_ind(1,1)
            if (xind.gt.surf_ind(1,2)) xind = surf_ind(1,2)
            if (zind.lt.surf_ind(2,1)) zind = surf_ind(2,1)
            if (zind.gt.surf_ind(2,2)) zind = surf_ind(2,2) 
            
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
  subroutine revaluate_heat_domain(lev, xlo, dx, lo, hi, &
                                   idom_old, ido_lo, ido_hi, &
                                   idom_new, idn_lo, idn_hi, &
                                   u_in, u_lo, u_hi, &
                                   temp, t_lo, t_hi, &
                                   u_in2, u2_lo, u2_hi, &
                                   temp2, t2_lo, t2_hi)

    use material_properties_module, only : temp_melt
    use amr_data_module, only : surf_ind, &
                                surf_temperature, &
                                surf_enthalpy, &
                                surf_dx, &
                                surf_xlo, &
                                melt_vel
    
    ! Input and output variables
    integer, intent(in) :: lev
    integer, intent(in) :: lo(3), hi(3) ! bounds of current tile box
    integer, intent(in) :: u_lo(3), u_hi(3) ! bounds of input enthalpy box 
    integer, intent(in) :: u2_lo(3), u2_hi(3) ! bounds of input enthalpy box with 2 ghost points
    integer, intent(in) :: ido_lo(3), ido_hi(3) ! bounds of the input idomain box
    integer, intent(in) :: idn_lo(3), idn_hi(3) ! bounds of the output idomain box
    integer, intent(in) :: t_lo(3), t_hi(3) ! bounds of the temperature box
    integer, intent(in) :: t2_lo(3), t2_hi(3) ! bounds of the temperature box with two ghost points
    real(amrex_real), intent(in) :: xlo(3) ! Physical location of box boundaries
    real(amrex_real), intent(in) :: dx(3) ! Grid size
    real(amrex_real), intent(inout) :: u_in(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3)) ! Input enthalpy 
    real(amrex_real), intent(in) :: idom_old(ido_lo(1):ido_hi(1),ido_lo(2):ido_hi(2),ido_lo(3):ido_hi(3))
    real(amrex_real), intent(inout) :: idom_new(idn_lo(1):idn_hi(1),idn_lo(2):idn_hi(2),idn_lo(3):idn_hi(3))
    real(amrex_real), intent(inout) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(amrex_real), intent(in) :: temp2(t2_lo(1):t2_hi(1),t2_lo(2):t2_hi(2),t2_lo(3):t2_hi(3))
    real(amrex_real), intent(in) :: u_in2(u2_lo(1):u2_hi(1),u2_lo(2):u2_hi(2),u2_lo(3):u2_hi(3)) ! Input enthalpy 
       
    !Local variables
    integer :: i,j,k
    real(amrex_real) :: xpos, zpos
    integer :: xind, zind
    
    ! Re-evaluate domain
    do i = lo(1)-1,hi(1)+1
       do  j = lo(2)-1,hi(2)+1
          do k = lo(3)-1,hi(3)+1
             
             ! Points added to the domain
             if (nint(idom_old(i,j,k)).eq.0 .and. nint(idom_new(i,j,k)).ne.0) then

               ! Points added on top of melt
               if (temp2(i,j-1,k).gt.temp_melt) then

                  ! Update properties
                  u_in(i,j,k) = u_in2(i,j-1,k)
                  temp(i,j,k) = temp2(i,j-1,k)
                  idom_new(i,j,k) = 3
  
               ! Points added on top of mushy
               elseif (temp2(i,j-1,k).eq.temp_melt) then
  
                 ! Update properties
                 u_in(i,j,k) = u_in2(i,j-1,k)
                 temp(i,j,k) = temp2(i,j-1,k)
                 idom_new(i,j,k) = 2
  
               ! Points added on top of solid (take upwind temperature)   
               else
                  ! write(*,*) 'Point added on top of solid' 
  
                  ! Index for the surface properties (temperature and enthalpy)
                  xpos = xlo(1) + (0.5 + i-lo(1))*dx(1)  
                  xind = nint((xpos - surf_dx(1)/2 - surf_xlo(1))/surf_dx(1))
                  if (xind.lt.surf_ind(1,1)) xind = surf_ind(1,1)
                  if (xind.gt.surf_ind(1,2)) xind = surf_ind(1,2)  
                  zpos = xlo(3) + (0.5 + k-lo(3))*dx(3)  
                  zind = nint((zpos - surf_dx(2)/2 - surf_xlo(2))/surf_dx(2)) 
                  if (zind.lt.surf_ind(2,1)) zind = surf_ind(2,1)
                  if (zind.gt.surf_ind(2,2)) zind = surf_ind(2,2) 
                  
                  ! Figure out which column is upwind. Firsttry to find upwind in the direction of x, only if that doesn't work
                  ! look in the z direction.
                  if (melt_vel(xind,zind,1).gt.0_amrex_real) then
                    ! Boundary condition   
                    if (xind.gt.surf_ind(1,1)) xind = xind-1 
                  elseif (melt_vel(xind+1,zind,1).lt.0_amrex_real) then
                    xind = xind+1
                  elseif (melt_vel(xind,zind,2).gt.0_amrex_real) then
                     if (zind.gt.surf_ind(1,1)) zind = zind-1 
                  elseif (melt_vel(xind,zind+1,2).lt.0_amrex_real) then
                    zind = zind+1
                  else
                     if (lev.eq.amrex_max_level) then
                        write(*,*) & 
                        'Wind not blowing towards this column, cell shouldnt be added. Taking the column on the left as upwind'
                     end if
                    ! Boundary condition   
                    if (xind.gt.surf_ind(1,1)) xind = xind-1 
                  end if 
  
                  ! Update properties
                  u_in(i,j,k) = surf_enthalpy(xind,zind)
                  temp(i,j,k) = surf_temperature(xind,zind)
                  if (temp(i,j,k).gt.temp_melt) then
                     idom_new(i,j,k) = 3
                  else if (temp(i,j,k).eq.temp_melt) then
                     idom_new(i,j,k) = 2
                  else
                     idom_new(i,j,k) = 1
                  end if
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
  ! Subroutine used to get the total volume of molten material
  ! -----------------------------------------------------------------
  subroutine integrate_surf(melt_vol)
    
    use amr_data_module, only : melt_top, melt_pos, surf_ind, surf_dx
    
    ! Input and output variables
    real(amrex_real), intent(out) :: melt_vol ! Integrated melt volume [mm3]
    
    ! Local variables
    integer :: i,k 
    
    melt_vol = 0 
    
    do i =  surf_ind(1,1), surf_ind(1,2) 
       do k = surf_ind(2,1), surf_ind(2,2)
          melt_vol = melt_vol +  melt_top(i,k) - melt_pos(i,k)
       end do
    end do
    
    melt_vol = melt_vol*surf_dx(1)*surf_dx(2)*1E9  
    
  end subroutine integrate_surf
  

  subroutine get_local_highest_level(xlo, dx, lo, hi, lev)

   use read_input_module, only: surfdist

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
            do while(abs(ypos-surf_pos_local(i,k)).lt.surfdist(lev_test))
                lev_test = lev_test+1
            end do
            lev(i,j,k) = lev_test-1
         end do
      end do
   end do

  end subroutine get_local_highest_level
  
end module domain_module



