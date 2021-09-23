module heat_transfer_fixT_module
  
  ! -----------------------------------------------------------------
  ! This module is used to perform all the calculations relative
  ! to the heat transfer part of the code. All the calculations are
  ! performed with a fixed temperature on the free surface and only
  ! the diffusive part of the heat transfer equation is solved.
  ! -----------------------------------------------------------------
  
  use amrex_amr_module
  
  implicit none

  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: increment_enthalpy_fixT

contains

  ! -----------------------------------------------------------------
  ! Subroutine used to compute the enthalpy at a new time step
  ! -----------------------------------------------------------------
  subroutine increment_enthalpy_fixT(lo, hi, &
                                     u_old,  uo_lo, uo_hi, &
                                     u_new, un_lo, un_hi, &
                                     temp_old, to_lo, to_hi, &
                                     temp_new, tn_lo, tn_hi , &
                                     flxx, fx_lo, fx_hi, &
                                     flxy, fy_lo, fy_hi, &
                                     idom_old, ido_lo, ido_hi, &
                                     idom_new, idn_lo, idn_hi, &
                                     geom, dt)
    
    use material_properties_module, only : get_temp, get_enthalpy
    use read_input_module, only : temp_fs
    
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
    real(amrex_real), intent(inout) :: u_old(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2)) ! Input enthalpy 
    real(amrex_real), intent(inout) :: u_new(un_lo(1):un_hi(1),un_lo(2):un_hi(2)) ! Output enthalpy
    real(amrex_real), intent(inout) :: temp_old(to_lo(1):to_hi(1),to_lo(2):to_hi(2)) ! Input temperature
    real(amrex_real), intent(inout) :: temp_new(tn_lo(1):tn_hi(1),tn_lo(2):tn_hi(2)) ! Output temperature
    real(amrex_real), intent(out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2)) ! flux along the x direction  			
    real(amrex_real), intent(out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2)) ! flux along the y direction
    real(amrex_real), intent(in) :: idom_old(ido_lo(1):ido_hi(1),ido_lo(2):ido_hi(2))
    real(amrex_real), intent(in) :: idom_new(idn_lo(1):idn_hi(1),idn_lo(2):idn_hi(2))
    type(amrex_geometry), intent(in) :: geom ! geometry
    
    !Local variables
    integer :: i,j
    real(amrex_real) :: u_fs ! Enthalpy of the points on the free surface
    real(amrex_real) :: dx(2) ! Grid size
    real(amrex_real) :: lo_phys(2) ! Physical location of the lowest corner of the tile box

    ! Enthalpy at the free surface
    call get_enthalpy(temp_fs, u_fs)
    
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
                          flxx, fx_lo, fx_hi, &
                          flxy, fy_lo, fy_hi, &
                          temp_old, to_lo, to_hi)

    ! Compute enthalpy at the new timestep
    do i = lo(1),hi(1)
       do  j = lo(2),hi(2)
          
          if (nint(idom_new(i,j)).ne.0 .and. nint(idom_new(i,j+1)).eq.0) then
             u_new(i,j) = u_fs ! Impose temperature on the free surface
          else if(nint(idom_new(i,j)).eq.0) then
             u_new(i,j) = u_fs ! Set background temperature equal to the free surface temperature
          else
             ! Update enthalpy according to heat equation
             u_new(i,j) = u_old(i,j) &
                          - dt/dx(1) * (flxx(i+1,j) - flxx(i,j)) & ! flux divergence x-direction 
                          - dt/dx(2) * (flxy(i,j+1) - flxy(i,j)) ! flux divergence y-direction
          end if 

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
    
  end subroutine increment_enthalpy_fixT

  ! -----------------------------------------------------------------
  ! Subroutine used to the enthalpy fluxes on the edges of the grid
  ! -----------------------------------------------------------------  
  subroutine create_face_flux(dx, lo, hi, &
                              flxx, fx_lo, fx_hi, &
                              flxy, fy_lo, fy_hi, &
                              temp, t_lo, t_hi)
  				
    use material_properties_module, only: get_ktherm

    ! Input and output variables 
    integer, intent(in) :: lo(2), hi(2)  
    integer, intent(in) :: t_lo(2), t_hi(2)
    integer, intent(in) :: fx_lo(2), fx_hi(2)
    integer, intent(in) :: fy_lo(2), fy_hi(2)
    real(amrex_real), intent(in) :: dx(2) 
    real(amrex_real), intent(out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))
    real(amrex_real), intent(out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2))
    real(amrex_real), intent(in) :: temp (t_lo(1):t_hi(1),t_lo(2):t_hi(2))

    ! Local variables
    integer :: i,j
    real(amrex_real) :: ktherm
    real(amrex_real) :: temp_face

 
    
    ! Flux along the x direction
    do i = lo(1), hi(1)+1
       do j = lo(2), hi(2)          

          ! Diffusive component
          temp_face = (temp(i,j) + temp(i-1,j))/2_amrex_real
          call get_ktherm(temp_face, ktherm)
          flxx(i,j) = -ktherm*(temp(i,j)-temp(i-1,j))/dx(1)

       end do
    end do
    
    ! Flux along the y direction
    do i = lo(1), hi(1)
       do j = lo(2), hi(2)+1
          
          ! Diffusive component (there is no advection in the y direction)
          temp_face = (temp(i,j) + temp(i,j-1))/2_amrex_real
          call get_ktherm(temp_face, ktherm)
          flxy(i,j) = -ktherm*(temp(i,j)-temp(i,j-1))/dx(2)
          
       end do
    end do
 
  end subroutine create_face_flux
  
end module heat_transfer_fixT_module
