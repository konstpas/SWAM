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
  
contains
  
  
  ! -----------------------------------------------------------------
  ! Subroutine used to obtain the integer field used to distinguish
  ! between material and background
  ! -----------------------------------------------------------------
  subroutine get_idomain(xlo, dx, lo, hi, &
                         idom, id_lo, id_hi, &
                         temp, t_lo, t_hi)

    use read_input_module, only : geometry_name 
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: id_lo(2), id_hi(2)
    integer, intent(in) :: t_lo(2), t_hi(2)
    real(amrex_real), intent(in) :: dx(2)    
    real(amrex_real), intent(inout) :: idom(id_lo(1):id_hi(1), id_lo(2):id_hi(2))
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1), t_lo(2):t_hi(2))
    real(amrex_real), intent(in) :: xlo(2)

    ! Define idomains depending on the geometry specified in input
    if (geometry_name .eq. "Slab") then
       
       call get_slab_idomain(xlo, dx, lo, hi, &
                             idom, id_lo, id_hi, &
                             temp, t_lo, t_hi)
       
    elseif (geometry_name .eq. "West") then
      
       call get_west_idomain(xlo, dx, lo, hi, &
                             idom, id_lo, id_hi, &
                             temp, t_lo, t_hi)
       
    elseif (geometry_name .eq. "West_rectangular") then
       
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
    logical :: inside_domain
    integer :: i,j
    integer :: surf_ind_heat_domain
    real(amrex_real) :: surf_pos_heat_domain(id_lo(1):id_hi(1))
    real(amrex_real) :: xpos_face
    real(amrex_real) :: ypos_face

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

       ! Index of the free surface in the heat transfer domain
       surf_ind_heat_domain = id_lo(2) + &
                              floor((surf_pos_heat_domain(i) - &
                              xlo(2) + dx(2))/dx(2))

       xpos_face = xlo(1) + (i-lo(1))*dx(1)
       
       do j = lo(2)-1, hi(2)+1

          ypos_face = xlo(2) + (j-lo(2))*dx(2)

          ! Initialization
          inside_domain = .true.
          if (xpos_face.lt.amrex_problo(1) .or. xpos_face.gt.amrex_probhi(1) .or. &
              ypos_face.lt.amrex_problo(2) .or. ypos_face.gt.amrex_probhi(2)) then
               inside_domain = .false.
          end if

          if (j .le. surf_ind_heat_domain .and. inside_domain) then

             if (find_liquid) then
                if (temp(i,j).gt.temp_melt) then
                   idom(i,j) = 3 ! Liquid
                else if (temp(i,j).eq.temp_melt) then
                   idom(i,j) = 2 ! Liquid
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
    ! between material and background in WEST geometry
    ! -----------------------------------------------------------------
    subroutine get_west_idomain(xlo, dx, lo, hi, &
                           idom, id_lo, id_hi, &
                           temp, t_lo, t_hi)
  
      use material_properties_module, only : temp_melt  
      use read_input_module, only : sample_edge, &
                                    cool_pipe_cntr, &
                                    cool_pipe_radius  
      
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
      real(amrex_real) :: xpos, ypos
      logical :: pipe_flag
  
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

         xpos = xlo(1) + (0.5+i-lo(1))*dx(1) 
         
         do j = lo(2)-1, hi(2)+1  
            
            ypos = xlo(2) + (0.5+j-lo(2))*dx(2) 
            if (sqrt((xpos-cool_pipe_cntr(1))**2+(ypos-cool_pipe_cntr(2))**2).lt.cool_pipe_radius) then
              pipe_flag = .true.
            else
              pipe_flag = .false.
            end if

            if (j .le. surf_ind_heat_domain .and. xpos .le. sample_edge .and. (.not.pipe_flag)) then
               if (find_liquid) then
                  if (temp(i,j).gt.temp_melt) then
                      idom(i,j) = 3
                  else if (temp(i,j).eq.temp_melt) then
                      idom(i,j) = 2
                  else
                      idom(i,j) = 1
                  end if
               else
                  idom(i,j) = 1
               end if
               
            elseif (pipe_flag) then
               idom(i,j) = -1
            else
               idom(i,j) = 0
            end if
            
         end do
         
      end do
  
    end subroutine get_west_idomain


    ! -----------------------------------------------------------------
    ! Subroutine used to obtain the integer field used to distinguish
    ! between material and background in WEST-like geometry with rectangular
    ! cross-section cooling pipes
    ! -----------------------------------------------------------------
    subroutine get_west_rectangular_idomain(xlo, dx, lo, hi, &
                           idom, id_lo, id_hi, &
                           temp, t_lo, t_hi)
  
      use material_properties_module, only : temp_melt  
      use read_input_module, only : sample_edge, &
                                    cool_pipe_cntr, &
                                    cool_pipe_radius  
      
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
      real(amrex_real) :: xpos, ypos
      logical :: pipe_flag
      real(amrex_real), parameter :: pi = 3.14159265358979
  
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

         xpos = xlo(1) + (0.5+i-lo(1))*dx(1) 
         
         do j = lo(2)-1, hi(2)+1  
            
            ypos = xlo(2) + (0.5+j-lo(2))*dx(2) 
            if((xpos.lt.cool_pipe_cntr(1)+(pi)*cool_pipe_radius/4) .and. & 
             (xpos.gt.cool_pipe_cntr(1)-(pi)*cool_pipe_radius/4) .and. &
             (ypos.lt.cool_pipe_cntr(2)+(pi)*cool_pipe_radius/4) .and. &
             (ypos.gt.cool_pipe_cntr(2)-(pi)*cool_pipe_radius/4)) then
              pipe_flag = .true.
            else
              pipe_flag = .false.
            end if

            if (j .le. surf_ind_heat_domain .and. xpos .le. sample_edge .and. (.not.pipe_flag)) then
               if (find_liquid) then
                  if (temp(i,j).gt.temp_melt) then
                      idom(i,j) = 3
                  else if (temp(i,j).eq.temp_melt) then
                      idom(i,j) = 2
                  else
                      idom(i,j) = 1
                  end if
               else
                  idom(i,j) = 1
               end if
               
            elseif (pipe_flag) then
               idom(i,j) = -1
            else
               idom(i,j) = 0
            end if
            
         end do
         
      end do
  
    end subroutine get_west_rectangular_idomain
    

  
  ! -----------------------------------------------------------------
  ! Subroutine used to interpolate the free surface position as given
  ! by the fluid solver in order to construct the heat conduction
  ! free interface. Note that the surface position is defined on
  ! the faces of the cells and not on the centers
  ! -----------------------------------------------------------------     
  subroutine get_surf_pos(xlo, dx, lo, hi, surf_pos_heat_domain)

    
    use amr_data_module, only : surf_ind, &
                                surf_pos, &
                                surf_xlo, &
                                surf_dx  

    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2) 
    real(amrex_real), intent(in) :: xlo(2)
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(out) :: surf_pos_heat_domain(lo(1):hi(1))

    ! Local variables
    integer :: i
    integer :: xind
    real(amrex_real) :: xpos
   !  real(amrex_real) :: frac_part
   !  real(amrex_real) :: num_tol
    
   !  num_tol = 1e-5
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
       if (xind.gt.surf_ind(1,2)) xind = surf_ind(1,2) 
      
       ! Check if you fall inbetween two grid points of the heighest level
      !  frac_part = (xpos - surf_dx(1)/2 - surf_xlo(1))/surf_dx(1) - xind
      !  if (frac_part.gt.0.5-num_tol .and. frac_part.lt.0.5+num_tol ) then
      !    surf_pos_heat_domain(i) = (surf_pos(xind)+surf_pos(xind+1))/2
      !  else
         ! surf_pos_heat_domain(i) = surf_pos(xind)
      !  end if

       surf_pos_heat_domain(i) = surf_pos(xind)
       
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
  ! Subroutine used to get the position of the bottom of the melt
  ! pool
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
        
    do i = lo(1), hi(1)  ! x-direction
       do j = lo(2), hi(2) 
             
          if (nint(idom(i,j)).ge.2 .and. nint(idom(i,j-1)).lt.2) then
             
             it(1) = i
             it(2) = j
             grid_pos = geom%get_physical_location(it)
             melt_pos(i) = grid_pos(2) 
             
          elseif(nint(idom(i,j)).lt.2 .and. nint(idom(i,j-1)).ge.2) then

            
             it(1) = i
             it(2) = j
             grid_pos = geom%get_physical_location(it)
             melt_top(i) = grid_pos(2)
            
          end if

          if (nint(idom(i,j)).eq.0 .and. nint(idom(i,j-1)).gt.0) then
             
            it(1) = i
            it(2) = j
            grid_pos = geom%get_physical_location(it)

          end if
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
    integer, intent(in) :: lo(2), hi(2) ! bounds of current tile box
    integer, intent(in) :: u_lo(2), u_hi(2) ! bounds of input enthalpy box 
    integer, intent(in) :: u2_lo(2), u2_hi(2) ! bounds of input enthalpy box with 2 ghost points
    integer, intent(in) :: ido_lo(2), ido_hi(2) ! bounds of the input idomain box
    integer, intent(in) :: idn_lo(2), idn_hi(2) ! bounds of the output idomain box
    integer, intent(in) :: t_lo(2), t_hi(2) ! bounds of the temperature box
    integer, intent(in) :: t2_lo(2), t2_hi(2) ! bounds of the temperature box with 2 ghost points
    real(amrex_real), intent(in) :: xlo(2) ! Physical location of box boundaries
    real(amrex_real), intent(in) :: dx(2) ! Grid size
    real(amrex_real), intent(inout) :: u_in(u_lo(1):u_hi(1),u_lo(2):u_hi(2)) ! Input enthalpy 
    real(amrex_real), intent(in) :: idom_old(ido_lo(1):ido_hi(1),ido_lo(2):ido_hi(2))
    real(amrex_real), intent(inout) :: idom_new(idn_lo(1):idn_hi(1),idn_lo(2):idn_hi(2))
    real(amrex_real), intent(inout) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2))
    real(amrex_real), intent(in) :: u_in2(u2_lo(1):u2_hi(1),u2_lo(2):u2_hi(2))
    real(amrex_real), intent(in) :: temp2(t2_lo(1):t2_hi(1),t2_lo(2):t2_hi(2))
     
    !Local variables
    integer :: i,j
    integer :: xind
    real(amrex_real) :: xpos

    ! NOTE: as it is now, the domain revaluation is correct only if the velocity is positive
    ! i.e. if the velocity goes in the direction of positive x
    
    do i = lo(1)-1,hi(1)+1
       do  j = lo(2)-1,hi(2)+1

          ! Points added to the domain
          if (nint(idom_old(i,j)).eq.0 .and. nint(idom_new(i,j)).ne.0) then

             ! Points added on top of melt
             if (temp2(i,j-1).gt.temp_melt) then

                ! Update properties
                u_in(i,j) = u_in2(i,j-1)
                temp(i,j) = temp2(i,j-1)
                idom_new(i,j) = 3

            ! ! Points added on top of mushy
             elseif (temp2(i,j-1).eq.temp_melt) then

               ! Update properties
               u_in(i,j) = u_in2(i,j-1)
               temp(i,j) = temp2(i,j-1)
               idom_new(i,j) = 2

             ! Points added on top of solid (take upwind temperature)   
             else

                ! Index for the surface properties (temperature and enthalpy)
                xpos = xlo(1) + (0.5 + i-lo(1))*dx(1)  
                xind = nint((xpos - surf_dx(1)/2 - surf_xlo(1))/surf_dx(1)) 
                if (xind.lt.surf_ind(1,1)) xind = surf_ind(1,1)
                if (xind.ge.surf_ind(1,2)) xind = surf_ind(1,2)-1 
                
                ! Figure out which column is upwind
                if (melt_vel(xind,1).gt.0_amrex_real) then
                  ! Boundary condition   
                  if (xind.gt.surf_ind(1,1)) xind = xind-1 
                elseif (melt_vel(xind+1,1).lt.0_amrex_real) then
                  xind = xind+1
                else
                  write(*,*) 'Wind not blowing towards this column, cell shouldnt be added. Taking the column on the left as upwind'
                  ! Boundary condition   
                  if (xind.gt.surf_ind(1,1)) xind = xind-1 
                end if 

                ! Update properties
                u_in(i,j) = surf_enthalpy(xind)
                temp(i,j) = surf_temperature(xind)
                if (temp(i,j).gt.temp_melt) then
                   idom_new(i,j) = 3
                else if (temp(i,j).eq.temp_melt) then
                   idom_new(i,j) = 2
                else
                   idom_new(i,j) = 1
                end if
                
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

  
end module domain_module
