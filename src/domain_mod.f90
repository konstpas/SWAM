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
  
  contains
    
    
    ! -----------------------------------------------------------------
    ! Subroutine used to obtain the integer field used to distinguish
    ! between material and background
    ! -----------------------------------------------------------------
    subroutine get_idomain(xlo, dx, lo, hi, &
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
  
    end subroutine get_idomain
  
    
    ! -----------------------------------------------------------------
    ! Subroutine used to interpolate the free surface position as given
    ! by the fluid solver in order to construct the heat conduction
    ! free interface. Note that the surface position is defined on
    ! the faces of the cells and not on the centers
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
            if (xind.ge.surf_ind(1,2)) xind = surf_ind(1,2)-1 
            if (zind.lt.surf_ind(2,1)) zind = surf_ind(2,1)
            if (zind.ge.surf_ind(2,2)) zind = surf_ind(2,2)-1 
  
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
               
               if (nint(idom(i,j,k)).eq.3 .and. nint(idom(i,j-1,k)).ne.3) then
                  
                  it(1) = i
                  it(2) = j
                  it(3) = k 
                  grid_pos = geom%get_physical_location(it)
                  melt_pos(i,k) = grid_pos(2) 
  
               else if (nint(idom(i,j,k)).eq.3 .and. nint(idom(i,j-1,k)).ne.3) then
  
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
    
    
  end module domain_module
  
