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
                         idom, id_lo, id_hi)
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: id_lo(2), id_hi(2)
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(inout) :: idom(id_lo(1):id_hi(1), id_lo(2):id_hi(2))
    real(amrex_real), intent(in) :: xlo(2)

    ! Local variables
    integer :: i,j
    integer :: surf_ind_heat_domain
    real(amrex_real) :: surf_pos_heat_domain(id_lo(1):id_hi(1))

    ! Get location of the free surface
    call get_surf_pos(xlo-dx, dx, id_lo, id_hi, surf_pos_heat_domain)

    ! Set flags to distinguish between material and background
    do i = lo(1)-1, hi(1)+1

       ! Free surface position in the heat solver domain
       surf_ind_heat_domain = id_lo(2) + &
                              floor((surf_pos_heat_domain(i) - &
                              xlo(2)+dx(2))/dx(2))

       ! Distinguish material and background
       do j = lo(2)-1, hi(2)+1

          if (j .le. surf_ind_heat_domain) then
             idom(i,j) = 1 ! Material
          else
             idom(i,j) = 0 ! Background
          end if
          
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
    real(amrex_real) :: x_alpha
    
    do  i = lo(1),hi(1)
   
       xpos = xlo(1) + (0.5 + i-lo(1))*dx(1) 
       
       ! In what follows -surf_dx(1)/2  and ceiling are used since
       ! the surface is staggered 'backwards' on the faces with
       ! respect to the xlo and xpos values which are defined on the
       ! centers. 
       xind = ceiling((xpos - surf_dx(1)/2 - surf_xlo(1))/surf_dx(1)) 
       x_alpha = mod(xpos - surf_dx(1)/2 - surf_xlo(1), surf_dx(1))
       
       if (xind.lt.surf_ind(1,1)) xind = surf_ind(1,1)
       if (xind.ge.surf_ind(1,2)) xind = surf_ind(1,2)-1 
       
       ! Interpolation
       surf_pos_heat_domain(i) = surf_pos(xind) + &
            x_alpha * (surf_pos(xind+1)-surf_pos(xind))
       
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
  subroutine get_melt_pos(lo, hi, &
                          idom, id_lo, id_hi, &
                          temp, t_lo, t_hi, &
                          geom)
       
    use amr_data_module, only : melt_pos, &
                                melt_top

    use material_properties_module, only : temp_melt
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2) 
    integer, intent(in) :: id_lo(2), id_hi(2)
    integer, intent(in) :: t_lo(2), t_hi(2)    
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2))
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2))
    type(amrex_geometry), intent(in) :: geom
    
    ! Local variables
    logical :: check_warning
    integer :: i,j
    integer :: it(1:2) 
    real(amrex_real) :: grid_pos(1:2)
    
    do i = lo(1), hi(1)  ! x-direction
       do j = lo(2), hi(2) 

          ! Bottom of the melt pool
          if (temp(i,j).gt.temp_melt .and. temp(i,j-1).le.temp_melt) then

             it(1) = i
             it(2) = j
             grid_pos = geom%get_physical_location(it)
             melt_pos(i) = grid_pos(2)

          ! Top of the melt pool (should coincide with the free surface)
          else if (temp(i,j).gt.temp_melt .and. temp(i,j+1).le.temp_melt) then

            it(1) = i
            it(2) = j
            grid_pos = geom%get_physical_location(it)
            melt_top(i) = grid_pos(2)

            ! Consistency check: melt top should coincide with the free
            ! surface, i.e. the location where the idomains switch from
            ! one to zero
            if (nint(idom(i,j+1)).ne.0 .and. check_warning) then
               print *,  'WARNING: Melt top not at free surface'
               check_warning = .false.
            end if
             
          end if

       end do   
    end do
    
  end subroutine get_melt_pos

  ! -----------------------------------------------------------------
  ! Subroutine used to re-evaluate the heat equation domain
  ! -----------------------------------------------------------------  
  subroutine revaluate_heat_domain(lo, hi, &
                                   idom_old, ido_lo, ido_hi, &
                                   idom_new, idn_lo, idn_hi, &
                                   u_in, u_lo, u_hi)

    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2) ! bounds of current tile box
    integer, intent(in) :: u_lo(2), u_hi(2) ! bounds of input enthalpy box 
    integer, intent(in) :: ido_lo(2), ido_hi(2) ! bounds of the input idomain box
    integer, intent(in) :: idn_lo(2), idn_hi(2) ! bounds of the output idomain box
    real(amrex_real), intent(inout) :: u_in(u_lo(1):u_hi(1),u_lo(2):u_hi(2)) ! Input enthalpy 
    real(amrex_real), intent(in) :: idom_old(ido_lo(1):ido_hi(1),ido_lo(2):ido_hi(2))
    real(amrex_real), intent(in) :: idom_new(idn_lo(1):idn_hi(1),idn_lo(2):idn_hi(2))
    
    !Local variables
    integer :: i,j
    
    do i = lo(1)-1,hi(1)+1
       do  j = lo(2)-1,hi(2)+1

          ! Points added to the domain
          if (nint(idom_old(i,j)).eq.0 .and. nint(idom_new(i,j)).ne.0) then
             u_in(i,j) = u_in(i,j-1)
          ! Points removed from the domain
          else if (nint(idom_new(i,j)).eq.0) then
             u_in(i,j) = 0.0_amrex_real
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
