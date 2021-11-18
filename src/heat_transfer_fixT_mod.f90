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
                                     temp_new, tn_lo , tn_hi , &
                                     flxx, fx_lo, fx_hi, &
                                     flxy, fy_lo, fy_hi, &
                                     flxz, fz_lo, fz_hi, &
                                     idom_old, ido_lo, ido_hi, &
                                     idom_new, idn_lo, idn_hi, &
                                     geom, dt)

    use material_properties_module, only : get_temp, get_enthalpy
    use read_input_module, only : temp_fs

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
    real(amrex_real) :: u_fs ! Enthalpy of the points on the free surface
    real(amrex_real) :: dx(3) ! Grid size
    real(amrex_real) :: lo_phys(3) ! Physical location of the lowest corner of the tile box

    ! Enthalpy at the free surface
    call get_enthalpy(temp_fs, u_fs)
    
    ! Re-evaluate domain
    do i = lo(1)-1,hi(1)+1
       do  j = lo(2)-1,hi(2)+1
          do k = lo(3)-1,hi(3)+1
             
             ! Points added to the domain
             if (nint(idom_old(i,j,k)).eq.0 .and. nint(idom_new(i,j,k)).ne.0) then
                u_old(i,j,k) = u_old(i,j-1,k)
             ! Points removed from the domain
             else if (nint(idom_new(i,j,k)).eq.0) then
                u_old(i,j,k) = 0_amrex_real
             end if

          end do
       end do
    end do

    ! Get grid size
    dx = geom%dx(1:3) ! grid width at level 

    ! Get physical location of the lowest corner of the tile box
    lo_phys = geom%get_physical_location(lo)
   
    ! Get enthalpy flux 
    call create_face_flux(dx, lo, hi, &
                          flxx, fx_lo, fx_hi, &
                          flxy, fy_lo, fy_hi, &
                          flxz, fz_lo, fz_hi, &
                          temp_old, to_lo, to_hi)
  				  	   
    ! Compute enthalpy at the new timestep
    do   i = lo(1),hi(1)
       do  j = lo(2),hi(2) 
          do k = lo(3),hi(3)

             if (nint(idom_new(i,j,k)).ne.0 .and. nint(idom_new(i,j+1,k)).eq.0) then
                u_new(i,j,k) = u_fs ! Impose temperature on the free surface
             else if(nint(idom_new(i,j,k)).eq.0) then
                u_new(i,j,k) = u_fs ! Set background temperature equal to the free surface temperature
             else
                ! Update enthalpy according to heat equation
                u_new(i,j,k) = u_old(i,j,k) &
                               - dt/dx(1) * (flxx(i+1,j,k) - flxx(i,j,k)) &	! flux divergence x-direction 
                               - dt/dx(2) * (flxy(i,j+1,k) - flxy(i,j,k)) &	! flux divergence y-direction 
                               - dt/dx(3) * (flxz(i,j,k+1) - flxz(i,j,k))       ! flux divergence z-direction
             end if
             
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
                  temp_new, tn_lo , tn_hi, .true.) 

  end subroutine increment_enthalpy_fixT


  ! -----------------------------------------------------------------
  ! Subroutine used to the enthalpy fluxes on the edges of the grid
  ! -----------------------------------------------------------------  
  subroutine create_face_flux(dx, lo, hi, &
                              flxx, fx_lo, fx_hi, &
                              flxy, fy_lo, fy_hi, &
                              flxz, fz_lo, fz_hi, &
                              temp, t_lo, t_hi)
  				
    use material_properties_module, only: get_conductivity

    ! Input and output variables
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: t_lo(3), t_hi(3)
    integer, intent(in) :: fx_lo(3), fx_hi(3)
    integer, intent(in) :: fy_lo(3), fy_hi(3)
    integer, intent(in) :: fz_lo(3), fz_hi(3)
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3))
    real(amrex_real), intent(out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3))
    real(amrex_real), intent(out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3))
    real(amrex_real), intent(in) :: temp (t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))

    ! Local variables
    integer :: i,j,k 
    real(amrex_real) :: ktherm
    real(amrex_real) :: temp_face
    
    ! Flux along the x direction
    do i = lo(1), hi(1)+1
       do j = lo(2), hi(2)
          do k = lo(3), hi(3)
                
             ! Diffusive component
             temp_face = (temp(i,j,k) + temp(i-1,j,k))/2_amrex_real
             call get_conductivity(temp_face, ktherm)
             flxx(i,j,k) = -ktherm*(temp(i,j,k)-temp(i-1,j,k))/dx(1)
                                
          end do
       end do
    end do

    ! Flux along the y direction
    do i = lo(1), hi(1)
       do j = lo(2), hi(2)+1
          do k = lo(3), hi(3)
             
             ! Diffusive component (there is no advection in the y direction)
             temp_face = (temp(i,j,k) + temp(i,j-1,k))/2_amrex_real
             call get_conductivity(temp_face, ktherm)
             flxy(i,j,k) = -ktherm*(temp(i,j,k)-temp(i,j-1,k))/dx(2)
             
          end do
       end do
    end do

    ! Flux along the z direction
    do i = lo(1), hi(1)
       do j = lo(2), hi(2)
          do k = lo(3), hi(3)+1
                
             ! Diffusive component
             temp_face = (temp(i,j,k) + temp(i,j,k-1))/2_amrex_real
             call get_conductivity(temp_face, ktherm)
             flxz(i,j,k) = -ktherm*(temp(i,j,k)-temp(i,j,k-1))/dx(3)
                
          end do
       end do
    end do
    
    
  end subroutine create_face_flux
  
end module heat_transfer_fixT_module
