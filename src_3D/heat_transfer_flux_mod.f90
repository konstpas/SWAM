module heat_flux_module
  
  ! -----------------------------------------------------------------
  ! This module is used to compute the heat flux imposed on the free
  ! surface
  ! -----------------------------------------------------------------
  
  use amrex_amr_module
  
  implicit none

  private

  ! -----------------------------------------------------------------
  ! Public subroutines
  ! -----------------------------------------------------------------
  public :: get_boundary_heat_flux
  public :: debug_cooling_fluxes
   
  
contains
  
  ! -----------------------------------------------------------------
  ! Subroutine used to prescribe the boundary heating on the free
  ! surface. Note that the incoming heat flux is assigned to the
  ! first node under the free surface in the form of a volumetric
  ! heating (after an appropriate dimensionality correction)
  ! -----------------------------------------------------------------   
  subroutine get_boundary_heat_flux(time, xlo, &
                                    dx, lo, hi, &
                                    idom, id_lo, id_hi, &
                                    temp, t_lo, t_hi, lev, dt, &
                                    Qpipe_box, Qtherm_box, Qvap_box, Qrad_box, &
                                    Qplasma_box, qb)

    use read_input_module, only : heat_plasma_flux_type, &
                                  heat_cooling_thermionic, &
                                  heat_cooling_thermionic_side, &
                                  heat_cooling_vaporization, &
                                  heat_cooling_radiation, &
                                  geom_name, &
                                  heat_plasma_flux_side_type, &
                                  heat_sample_edge

    use amr_data_module, only : surf_current

    use heat_transfer_domain_module, only : get_local_highest_level

    ! Input and output variables
    integer, intent(in) :: lo(3), hi(3)  
    integer, intent(in) :: id_lo(3), id_hi(3)
    integer, intent(in) :: t_lo(3), t_hi(3)
    real(amrex_real), intent(in) :: time
    real(amrex_real), intent(in) :: xlo(3)
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2),id_lo(3):id_hi(3))
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    integer, intent(in) :: lev
    real(amrex_real), intent(in) :: dt
    real(amrex_real), intent(out) :: Qplasma_box
    real(amrex_real), intent(out) :: Qpipe_box
    real(amrex_real), intent(out) :: Qtherm_box
    real(amrex_real), intent(out) :: Qrad_box
    real(amrex_real), intent(out) :: Qvap_box
    real(amrex_real), intent(out) :: qb(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

    ! Local variables
    integer :: i, j, k
    real(amrex_real) :: q_plasma
    real(amrex_real) :: q_vap, q_rad, q_therm, q_cool
    real(amrex_real) :: xpos, ypos, zpos
    logical :: side_flag
    integer :: local_max_level(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

    qb = 0.0
    q_plasma = 0.0
    q_rad = 0.0
    q_vap = 0.0
    q_therm = 0.0
    q_cool = 0.0
    local_max_level = 0.0
    Qpipe_box = 0.0
    Qplasma_box = 0.0
    Qtherm_box = 0.0
    Qrad_box = 0.0
    Qvap_box = 0.0

    call get_local_highest_level(xlo, dx, lo, hi, local_max_level)

    do i = lo(1), hi(1)
       do j = lo(2), hi(2)
          do k = lo(3), hi(3)
          
             ! Assign fluxes only on free surface
             if(nint(idom(i,j,k)).ne.0 .and. nint(idom(i,j+1,k)).eq.0) then
                
               side_flag = .false.
               q_plasma = 0.0
               q_therm = 0.0
               q_vap = 0.0
               q_rad = 0.0

                ! Location of the free surface
                xpos = xlo(1) + (i-lo(1))*dx(1)
                zpos = xlo(3) + (k-lo(3))*dx(3)
                
                ! Plasma flux
                if (heat_plasma_flux_type.eq.'Gaussian') then
                   call gaussian_heat_flux(time, xpos, zpos, side_flag, q_plasma)
                elseif (heat_plasma_flux_type.eq.'Gaussian_tube') then
                  call gaussian_tube_heat_flux(time, xpos, side_flag, q_plasma)
                elseif (heat_plasma_flux_type.eq.'Uniform') then
                  call uniform_heat_flux(time, xpos, zpos, side_flag, q_plasma)
                elseif (heat_plasma_flux_type.eq.'Input_file') then
                  call file_heat_flux(time, xpos, zpos, side_flag, q_plasma)
                else
                   STOP "Unknown plasma heat flux type"
                end if
                if (lev.eq.amrex_max_level) Qplasma_box = Qplasma_box + q_plasma*dx(1)*dx(3)*dt
                
                ! Thermionic cooling flux
                if (heat_cooling_thermionic) then
                  if(lev.eq.amrex_max_level) then
                     call thermionic_cooling(temp(i,j,k), q_plasma, q_therm, surf_current(i,k))
                  else
                     call thermionic_cooling(temp(i,j,k), q_plasma, q_therm)
                  end if
                  if (lev.eq.amrex_max_level) Qtherm_box = Qtherm_box + q_therm*dx(1)*dx(3)*dt
                end if
                
                ! Vaporization cooling flux
                if (heat_cooling_vaporization) then
                   call vaporization_cooling(temp(i,j,k), q_vap)
                   if (lev.eq.amrex_max_level) Qvap_box = Qvap_box+ q_vap*dx(1)*dx(3)*dt
                end if
                
                ! Radiative cooling flux
                if (heat_cooling_radiation) then
                   call radiation_cooling(temp(i,j,k), q_rad)
                   if (lev.eq.amrex_max_level) Qrad_box = Qrad_box +q_rad*dx(1)*dx(3)*dt
                end if
                
                ! Sum all flux contributions
                qb(i,j,k) = q_plasma - q_rad - q_vap - q_therm
                
                ! Note: the term /dx(2) converts a surface heat flux [W/m^2]
                ! into a volumetric heat flux [W/m^3]
                qb(i,j,k) = qb(i,j,k)/dx(2)           
            end if
               
             if(nint(idom(i,j,k)).ne.-1 .and. nint(idom(i,j+1,k)).eq.-1) then
                call active_cooling(temp(i,j,k), q_cool)
                qb(i,j,k) = qb(i,j,k) - q_cool/dx(2)
                if(lev.eq.local_max_level(i,j,k)) Qpipe_box = Qpipe_box + q_cool*dx(1)*dx(3)*dt
             end if
             if(nint(idom(i,j,k)).ne.-1 .and. nint(idom(i,j-1,k)).eq.-1) then
                call active_cooling(temp(i,j,k), q_cool)
                qb(i,j,k) = qb(i,j,k) - q_cool/dx(2)
                if(lev.eq.local_max_level(i,j,k)) Qpipe_box = Qpipe_box + q_cool*dx(1)*dx(3)*dt
             end if
             if(nint(idom(i,j,k)).ne.-1 .and. nint(idom(i+1,j,k)).eq.-1) then
                call active_cooling(temp(i,j,k), q_cool)
                qb(i,j,k) = qb(i,j,k) - q_cool/dx(1)
                if(lev.eq.local_max_level(i,j,k)) Qpipe_box = Qpipe_box + q_cool*dx(2)*dx(3)*dt
             end if
             if(nint(idom(i,j,k)).ne.-1 .and. nint(idom(i-1,j,k)).eq.-1) then
                call active_cooling(temp(i,j,k), q_cool)
                qb(i,j,k) = qb(i,j,k) - q_cool/dx(1)
                if(lev.eq.local_max_level(i,j,k)) Qpipe_box = Qpipe_box + q_cool*dx(2)*dx(3)*dt
             end if                   
             if(nint(idom(i,j,k)).ne.-1 .and. nint(idom(i,j,k+1)).eq.-1) then
                call active_cooling(temp(i,j,k), q_cool)
                qb(i,j,k) = qb(i,j,k) - q_cool/dx(3)
                if(lev.eq.local_max_level(i,j,k)) Qpipe_box = Qpipe_box + q_cool*dx(1)*dx(2)*dt
             end if
             if(nint(idom(i,j,k)).ne.-1 .and. nint(idom(i,j,k-1)).eq.-1) then
                call active_cooling(temp(i,j,k), q_cool)
                qb(i,j,k) = qb(i,j,k) - q_cool/dx(3)
                if(lev.eq.local_max_level(i,j,k)) Qpipe_box = Qpipe_box + q_cool*dx(1)*dx(2)*dt
             end if   


             zpos = xlo(3) + (1+k-lo(3))*dx(3)
             ypos = xlo(2) + (1+j-lo(2))*dx(2)
             if(geom_name .eq. "West" .and. nint(idom(i,j,k)).ne.0 .and. & 
               (xpos.ge.heat_sample_edge .or. nint(idom(i+1,j,k)).eq.0)) then           
                side_flag = .true.
                q_plasma = 0.0
                q_therm = 0.0
                q_vap = 0.0
                q_rad = 0.0

                ! Plasma flux
                if (heat_plasma_flux_side_type.eq.'Gaussian') then
                   call gaussian_heat_flux(time, ypos, zpos, side_flag, q_plasma)
                elseif (heat_plasma_flux_side_type.eq.'Gaussian_tube') then
                   call gaussian_tube_heat_flux(time, ypos, side_flag, q_plasma)
                elseif (heat_plasma_flux_side_type.eq.'Uniform') then
                  call uniform_heat_flux(time, ypos, zpos, side_flag, q_plasma)
                elseif (heat_plasma_flux_side_type.eq.'Input_file') then
                   call file_heat_flux (time, ypos, zpos, side_flag, q_plasma)
                else
                   STOP "Unknown plasma heat flux type"
                end if
                
                if (lev.eq.local_max_level(i,j,k)) Qplasma_box = Qplasma_box + q_plasma*dx(2)*dx(3)*dt
 
               !  ! Thermionic cooling flux
                if (heat_cooling_thermionic_side) then
                   call thermionic_cooling(temp(i,j,k), q_plasma, q_therm)
                   if (lev.eq.local_max_level(i,j,k)) Qtherm_box = Qtherm_box + q_therm*dx(2)*dx(3)*dt
                end if
 
                ! Vaporization cooling flux
                if (heat_cooling_vaporization) then
                   call vaporization_cooling(temp(i,j,k), q_vap)
                   if (lev.eq.local_max_level(i,j,k)) Qvap_box = Qvap_box + q_vap*dx(2)*dx(3)*dt
                end if
 
                ! Radiative cooling flux
                if (heat_cooling_radiation) then
                   call radiation_cooling(temp(i,j,k), q_rad)
                   if (lev.eq.local_max_level(i,j,k)) Qrad_box = Qrad_box + q_rad*dx(2)*dx(3)*dt
                end if
                qb(i,j,k) = qb(i,j,k) + (q_plasma-q_therm-q_vap-q_rad)/dx(1)
             end if
              

          end do
       end do
    end do
    
  end subroutine get_boundary_heat_flux
  

   ! -----------------------------------------------------------------
   ! Subroutine used to prescribe a gaussian heat flux active for
   ! time <= time_exposure
   ! -----------------------------------------------------------------   
   subroutine gaussian_heat_flux(time, xpos, zpos, side_flag, qb) 
 
     use read_input_module, only : heat_plasma_flux_params, &
                                   heat_plasma_flux_side_params
     
     ! Input and output variables
     real(amrex_real), intent(in) :: time
     real(amrex_real), intent(in) :: xpos
     real(amrex_real), intent(in) :: zpos
     logical, intent(in) :: side_flag
     real(amrex_real), intent(out) :: qb
     
     ! Local variables
     real(amrex_real), dimension(1:7) :: plasma_params

     qb = 0_amrex_real
     plasma_params = 0.0

     if (side_flag) then
      plasma_params = heat_plasma_flux_side_params
     else
      plasma_params = heat_plasma_flux_params
     end if
     
     if (time.ge.plasma_params(1) .and. time.le.plasma_params(2)) then
        qb = plasma_params(3) &
             *EXP(-((xpos-plasma_params(4))**2)/(plasma_params(5)**2) &
                  -((zpos-plasma_params(6))**2)/(plasma_params(7)**2))
     end if
     
   end subroutine gaussian_heat_flux


   ! -----------------------------------------------------------------
   ! Subroutine used to prescribe a heat flux with a gaussinan profile
   ! in the x-t plane and uniform profile in the z-t plane, active for
   ! time <= time_exposure
   ! -----------------------------------------------------------------   
   subroutine gaussian_tube_heat_flux(time, xpos, side_flag, qb) 
 
      use read_input_module, only : heat_plasma_flux_params, &
                                    heat_plasma_flux_side_params
      
      ! Input and output variables
      real(amrex_real), intent(in) :: time
      real(amrex_real), intent(in) :: xpos
      logical, intent(in) :: side_flag
      real(amrex_real), intent(out) :: qb
      
      ! Local variables
      real(amrex_real), dimension(1:5) :: plasma_params
      
      qb = 0_amrex_real
      plasma_params = 0.0
      
      if (side_flag) then
         plasma_params = heat_plasma_flux_side_params
      else
         plasma_params = heat_plasma_flux_params
      end if

      if (time.ge.plasma_params(1) .and. time.le.plasma_params(2)) then
         qb = plasma_params(3) &
              *EXP(-((xpos-plasma_params(4))**2)/(plasma_params(5)**2))
      end if
      
    end subroutine gaussian_tube_heat_flux   

   ! -----------------------------------------------------------------
   ! Subroutine used to prescribe a uniform heat flux active for
   ! time <= time_exposure
   ! -----------------------------------------------------------------   
   subroutine uniform_heat_flux(time, xpos, zpos, side_flag, qb) 
 
      use read_input_module, only : heat_plasma_flux_params, &
                                    heat_plasma_flux_side_params
      
      ! Input and output variables
      real(amrex_real), intent(in) :: time
      real(amrex_real), intent(in) :: xpos
      real(amrex_real), intent(in) :: zpos
      logical, intent(in) :: side_flag
      real(amrex_real), intent(out) :: qb
      
      ! Local variables
      real(amrex_real), dimension(1:7) :: plasma_params

      qb = 0_amrex_real
      plasma_params = 0.0
      
      if (side_flag) then
         plasma_params = heat_plasma_flux_side_params
        else
         plasma_params = heat_plasma_flux_params
      end if
      
      if (time.ge.plasma_params(1) .and. time.le.plasma_params(2)) then
         if(xpos.ge.plasma_params(4) .and. xpos.le.plasma_params(5) .and. &
          zpos.ge.plasma_params(6) .and. zpos.le.plasma_params(7)) then
            qb = plasma_params(3) 
         end if
      end if
      
    end subroutine uniform_heat_flux   
   
   ! -----------------------------------------------------------------
   ! Subroutine used to prescribe a heat flux defined in an
   ! input file.
   ! -----------------------------------------------------------------   
    subroutine file_heat_flux(time, xpos, zpos, side_flag, qb) 

      use amr_data_module, only : heat_flux_table, &
                                  plasma_flux_time_mesh, &
                                  plasma_flux_surf_x_mesh, &
                                  plasma_flux_surf_z_mesh, &
                                  heat_side_flux_table, &
                                  plasma_side_flux_time_mesh, &
                                  plasma_side_flux_surf_y_mesh, &
                                  plasma_side_flux_surf_z_mesh

      
      ! Input and output variables
      real(amrex_real), intent(in) :: time
      real(amrex_real), intent(in) :: xpos
      real(amrex_real), intent(in) :: zpos
      logical, intent(in) :: side_flag
      real(amrex_real), intent(out) :: qb

      ! Local variables 
      integer :: i_z, i_x, i_t, k, n, m
      real(amrex_real) :: z(2), x(2), t(2), val(8)
      real(amrex_real) :: txz_query(3)
      real(amrex_real), allocatable, dimension(:) :: spatial_mesh_1
      real(amrex_real), allocatable, dimension(:) :: spatial_mesh_2
      real(amrex_real), allocatable, dimension(:) :: temporal_mesh
      real(amrex_real), allocatable, dimension(:,:,:) :: heat_flux
      
      
      qb = 0_amrex_real
      
      if(side_flag) then
         n = size(plasma_side_flux_time_mesh,1)
         m = size(plasma_side_flux_surf_y_mesh,1)
         k = size(plasma_side_flux_surf_z_mesh,1)
         allocate (heat_flux(1:n,1:m,1:k))
         allocate (temporal_mesh(1:n))
         allocate (spatial_mesh_1(1:m))
         allocate (spatial_mesh_2(1:k))
         spatial_mesh_1 = plasma_side_flux_surf_y_mesh
         spatial_mesh_2 = plasma_side_flux_surf_z_mesh
         temporal_mesh = plasma_side_flux_time_mesh
         heat_flux = heat_side_flux_table
      else
         n = size(plasma_flux_time_mesh,1)
         m = size(plasma_flux_surf_x_mesh,1)
         k = size(plasma_flux_surf_z_mesh,1)
         allocate (heat_flux(1:n,1:m,1:k))
         allocate (temporal_mesh(1:n))
         allocate (spatial_mesh_1(1:m))
         allocate (spatial_mesh_2(1:k))
         spatial_mesh_1 = plasma_flux_surf_x_mesh
         spatial_mesh_2 = plasma_flux_surf_z_mesh
         temporal_mesh = plasma_flux_time_mesh
         heat_flux = heat_flux_table
      end if

      ! Find the maximum index i_t such that the time
      ! falls in-between temporal_mesh(i_t) and
      ! temporal_mesh(i_t+1). Similar for i_x and i_z      
      call bisection(temporal_mesh, n, time, i_t)
      call bisection(spatial_mesh_1, m, xpos, i_x)
      call bisection(spatial_mesh_2, k, zpos, i_z)

      ! Check if the query falls outside the domain
      ! defined by the heat flux mesh.
      if (i_t.eq.0 .or. i_x.eq.0 .or. i_z.eq.0 .or.&
          i_t.eq.n .or. i_x.eq.m .or. i_z.eq.k) then

         ! If query point is outside the 6 bounds
         ! then take the closest of the 8 corners
         if(i_t.eq.0 .and. i_x.eq.0 .and. i_z.eq.0) then
            qb = heat_flux(1,1,1)
         elseif(i_t.eq.0 .and. i_x.eq.0 .and. i_z.eq.k) then
            qb = heat_flux(1,1,k)
         elseif(i_t.eq.0 .and. i_x.eq.m .and. i_z.eq.0) then
            qb = heat_flux(1,m,1)
         elseif(i_t.eq.0 .and. i_x.eq.m .and. i_z.eq.k) then
            qb = heat_flux(1,m,k)
         elseif(i_t.eq.n .and. i_x.eq.0 .and. i_z.eq.0) then
            qb = heat_flux(n,1,1)
         elseif(i_t.eq.n .and. i_x.eq.0 .and. i_z.eq.k) then
            qb = heat_flux(n,1,k)
         elseif(i_t.eq.n .and. i_x.eq.m .and. i_z.eq.0) then
            qb = heat_flux(n,m,1)
         elseif(i_t.eq.n .and. i_x.eq.m .and. i_z.eq.k) then
            qb = heat_flux(n,m,k)

         ! If query point is outside 2 bounds then
         ! linear interpolation on one of the 12 edges
         elseif (i_t.eq.0 .and. i_x.eq.0) then
            z(1) = spatial_mesh_2(i_z)
            z(2) = spatial_mesh_2(i_z+1)
            val(1) = heat_flux(1, 1, i_z)
            val(2) = heat_flux(1, 1, i_z+1)
            call lin_intrp(z, val(1:2), zpos, qb)
         elseif (i_t.eq.0 .and. i_x.eq.m) then
            z(1) = spatial_mesh_2(i_z)
            z(2) = spatial_mesh_2(i_z+1)
            val(1) = heat_flux(1, m, i_z)
            val(2) = heat_flux(1, m, i_z+1)
            call lin_intrp(z, val(1:2), zpos, qb)
         elseif (i_t.eq.n .and. i_x.eq.0) then
            z(1) = spatial_mesh_2(i_z)
            z(2) = spatial_mesh_2(i_z+1)
            val(1) = heat_flux(n, 1, i_z)
            val(2) = heat_flux(n, 1, i_z+1)
            call lin_intrp(z, val(1:2), zpos, qb)
         elseif (i_t.eq.n .and. i_x.eq.m) then
            z(1) = spatial_mesh_2(i_z)
            z(2) = spatial_mesh_2(i_z+1)
            val(1) = heat_flux(n, m, i_z)
            val(2) = heat_flux(n, m, i_z+1)
            call lin_intrp(z, val(1:2), zpos, qb)


         elseif (i_t.eq.0 .and. i_z.eq.0) then
            x(1) = spatial_mesh_1(i_x)
            x(2) = spatial_mesh_1(i_x+1)
            val(1) = heat_flux(1, i_x, 1)
            val(2) = heat_flux(1, i_x+1, 1)
            call lin_intrp(x, val(1:2), xpos, qb)
         elseif (i_t.eq.0 .and. i_z.eq.k) then
            x(1) = spatial_mesh_1(i_x)
            x(2) = spatial_mesh_1(i_x+1)
            val(1) = heat_flux(1, i_x, k)
            val(2) = heat_flux(1, i_x+1, k)
            call lin_intrp(x, val(1:2), xpos, qb)
         elseif (i_t.eq.n .and. i_z.eq.0) then
            x(1) = spatial_mesh_1(i_x)
            x(2) = spatial_mesh_1(i_x+1)
            val(1) = heat_flux(n, i_x, 1)
            val(2) = heat_flux(n, i_x+1, 1)
            call lin_intrp(x, val(1:2), xpos, qb)
         elseif (i_t.eq.n .and. i_z.eq.k) then
            x(1) = spatial_mesh_1(i_x)
            x(2) = spatial_mesh_1(i_x+1)
            val(1) = heat_flux(n, i_x, k)
            val(2) = heat_flux(n, i_x+1, k)
            call lin_intrp(x, val(1:2), xpos, qb)


         elseif (i_x.eq.0 .and. i_z.eq.0) then
            t(1) = temporal_mesh(i_t)
            t(2) = temporal_mesh(i_t+1)
            val(1) = heat_flux(i_t, 1, 1)
            val(2) = heat_flux(i_t+1, 1, 1)
            call lin_intrp(t, val(1:2), time, qb)
         elseif (i_x.eq.0 .and. i_z.eq.k) then
            t(1) = temporal_mesh(i_t)
            t(2) = temporal_mesh(i_t+1)
            val(1) = heat_flux(i_t, 1, k)
            val(2) = heat_flux(i_t+1, 1, k)
            call lin_intrp(t, val(1:2), time, qb)   
         elseif (i_x.eq.m .and. i_z.eq.0) then
            t(1) = temporal_mesh(i_t)
            t(2) = temporal_mesh(i_t+1)
            val(1) = heat_flux(i_t, m, 1)
            val(2) = heat_flux(i_t+1, m, 1)
            call lin_intrp(t, val(1:2), time, qb) 
         elseif (i_x.eq.m .and. i_z.eq.k) then
            t(1) = temporal_mesh(i_t)
            t(2) = temporal_mesh(i_t+1)
            val(1) = heat_flux(i_t, m, k)
            val(2) = heat_flux(i_t+1, m, k)
            call lin_intrp(t, val(1:2), time, qb)
            
         ! If query point outside one of the domain bounds
         ! then bilinear interpolation in one the 6 surfaces.
         elseif(i_z.eq.0 .or. i_z.eq.k) then
            if (i_z.eq.0) i_z=1
            t(1) = temporal_mesh(i_t)
            t(2) = temporal_mesh(i_t+1)
            x(1) = spatial_mesh_1(i_x)
            x(2) = spatial_mesh_1(i_x+1)
            val(1) = heat_flux(i_t,i_x,i_z)
            val(2) = heat_flux(i_t+1,i_x,i_z)
            val(3) = heat_flux(i_t+1,i_x+1,i_z)
            val(4) = heat_flux(i_t,i_x+1,i_z)
            txz_query(1) = time
            txz_query(2) = xpos
            call bilin_intrp(t, x, val(1:4), txz_query(1:2), qb)

         elseif(i_x.eq.0 .or. i_x.eq.m) then
            if (i_x.eq.0) i_x=1
            t(1) = temporal_mesh(i_t)
            t(2) = temporal_mesh(i_t+1)
            z(1) = spatial_mesh_2(i_z)
            z(2) = spatial_mesh_2(i_z+1)
            val(1) = heat_flux(i_t,i_x,i_z)
            val(2) = heat_flux(i_t+1,i_x,i_z)
            val(3) = heat_flux(i_t+1,i_x,i_z+1)
            val(4) = heat_flux(i_t,i_x,i_z+1)
            txz_query(1) = time
            txz_query(2) = zpos
            call bilin_intrp(t, z, val(1:4), txz_query(1:2), qb)

         elseif(i_t.eq.0 .or. i_t.eq.n) then
            if (i_t.eq.0) i_t=1
            x(1) = spatial_mesh_1(i_x)
            x(2) = spatial_mesh_1(i_x+1)
            z(1) = spatial_mesh_2(i_z)
            z(2) = spatial_mesh_2(i_z+1)
            val(1) = heat_flux(i_t,i_x,i_z)
            val(2) = heat_flux(i_t,i_x+1,i_z)
            val(3) = heat_flux(i_t,i_x+1,i_z+1)
            val(4) = heat_flux(i_t,i_x,i_z+1)
            txz_query(1) = xpos
            txz_query(2) = zpos
            call bilin_intrp(x, z, val(1:4), txz_query(1:2), qb)
         endif
      
      ! In all other cases, trilinear interpolation   
      else
         t(1) = temporal_mesh(i_t)
         t(2) = temporal_mesh(i_t+1)
         x(1) = spatial_mesh_1(i_x)
         x(2) = spatial_mesh_1(i_x+1)
         z(1) = spatial_mesh_2(i_z)
         z(2) = spatial_mesh_2(i_z+1)
         val(1) = heat_flux(i_t,i_x,i_z)
         val(2) = heat_flux(i_t+1,i_x,i_z)
         val(3) = heat_flux(i_t+1,i_x+1,i_z)
         val(4) = heat_flux(i_t,i_x+1,i_z)
         val(5) = heat_flux(i_t,i_x,i_z+1)
         val(6) = heat_flux(i_t+1,i_x,i_z+1)
         val(7) = heat_flux(i_t+1,i_x+1,i_z+1)
         val(8) = heat_flux(i_t,i_x+1,i_z+1)
         txz_query(1) = time
         txz_query(2) = xpos
         txz_query(3) = zpos
         call trilin_intrp(t, x, z, val, txz_query, qb)
      endif

      
    end subroutine file_heat_flux    

   ! -----------------------------------------------------------------
   ! Subroutine used to find the surface cooling flux due to 
   ! radiation given the surface temperature
   ! -----------------------------------------------------------------   
   subroutine radiation_cooling(Ts, q_rad)
 
    use material_properties_module, only : get_emissivity
 
    ! Input and output variables
    real(amrex_real), intent(in) :: Ts     ! Temperature at the center of cells adjacent to the free surface [K]
    real(amrex_real), intent(out) :: q_rad ! Radiated power [W/m^2]
 
    ! Local variables 
    real(amrex_real) :: eps_t ! Emissivity
    real(amrex_real) :: sigma = 5.670374419E-8 ! Stefan-Boltzmann constant [W/(m^-2*K^-4)]

    call get_emissivity(Ts, eps_t)
    q_rad = sigma*eps_t*Ts**4
 
   end subroutine radiation_cooling


   ! -----------------------------------------------------------------
   ! Subroutine used to find the surface cooling flux due to 
   ! thermionic emission given the surface temperature temperature
   ! see E. Thorén et al. Plasma Phys. Control. Fusion 63 035021 (2021)
   ! -----------------------------------------------------------------   
   subroutine thermionic_cooling(Ts, q_plasma, q_therm, Jth)
 
     use material_properties_module, only : get_work_function, &
                                            get_richardson_constant
     
     use read_input_module, only : sw_magnetic_inclination
 
     ! Input and output variables                                       
     real(amrex_real), intent(in) :: Ts        ! Temperature at the center of cells adjacent to the free surface [K]
     real(amrex_real), intent(in) :: q_plasma  ! Plasma heat flux [K]
     real(amrex_real), intent(out) :: q_therm  ! Flux of energy due to thermionic emission [W/m^2]
     real(amrex_real), intent(out), optional :: Jth
     
     ! Local variables
     real(amrex_real) :: kb = 1.38064852E-23 ! Boltzmann constant [m^2*kg/(s^2*K)]
     real(amrex_real) :: Jth_nom
     real(amrex_real) :: Jth_lim
     real(amrex_real) :: Aeff
     real(amrex_real) :: Wf
     real(amrex_real) :: e = 1.60217662E-19
     real(amrex_real) :: pi = 3.1415927
     real(amrex_real) :: J
     
     call get_work_function(Wf)
     call get_richardson_constant(Aeff)

     ! Nominal thermionic current from the Richardson-Dushman formula
     Jth_nom = Aeff*EXP(-Wf/(kb*Ts))*Ts**2

     ! Space-charge limited current (semi-empirical expression)
     Jth_lim = 1.51e4 * q_plasma**(1.0/3.0) * (SIN(sw_magnetic_inclination/180*pi))**2

     ! Minimum between nominal and space-charge limited
     J = MIN(Jth_lim, Jth_nom)

     ! Heat flux
     q_therm = J/e*(Wf+2*kb*Ts)

     if (present(Jth)) Jth = J
     
   end subroutine thermionic_cooling
 
 
   ! -----------------------------------------------------------------
   ! Subroutine used to find the surface cooling flux due to 
   ! vaporization given the surface temperature temperature
   ! see E. Thorén et al. Plasma Phys. Control. Fusion 63 035021 (2021)
   ! -----------------------------------------------------------------   
    subroutine vaporization_cooling(Ts, q_vap)
 
      use material_properties_module, only : get_vapor_pressure, &
                                             get_enthalpy_of_vaporization, &
                                             get_atomic_mass
 
      ! Input and output variables
      real(amrex_real), intent(in) :: Ts
      real(amrex_real), intent(out) :: q_vap
      ! Local variables
      real(amrex_real) :: gm
      real(amrex_real) :: kb = 1.38064852E-23 ! Boltzmann constant [m^2*kg/(s^2*K)]
      real(amrex_real) :: h_vap ! Enthalpy of vaporization [kJ/mol]
      real(amrex_real) :: pv ! Vapor pressure
      real(amrex_real) :: m_A ! Atomic mass [g/mol]
      real(amrex_real) :: pi = 3.1415927
      real(amrex_real) :: Na = 6.02214076E23 ! Avogadro's number [mol^-1]  
      
      call get_atomic_mass(m_A)
      call get_vapor_pressure(Ts, pv)
      call get_enthalpy_of_vaporization(Ts, h_vap)
      
      ! Conversion from g/mol to kg
      m_A = m_A*1E-3/Na
      
      ! Conversion from kJ/mol to J
      h_vap = h_vap*1E3/Na
      
      gm = pv*sqrt(m_A/(2*pi*kb*Ts))    
      q_vap = gm/m_A*(h_vap + 2*kb*Ts) 
      
    end subroutine vaporization_cooling

    ! -----------------------------------------------------------------
    ! Subroutine used to prescribe the cooling flux from cooling pipes
    ! -----------------------------------------------------------------
    subroutine active_cooling(T, q_cool)
      ! Input and output variables
      real(amrex_real), intent(in) :: T
      real(amrex_real), intent(out) :: q_cool

      ! Local variables
      real(amrex_real) :: h   ! [W/K m^2]
      real(amrex_real) :: T0  ! [K]

      T0 = 343
      call get_convection_coeff(T, h)
      q_cool = h*(T-T0)
    end subroutine active_cooling

    ! -----------------------------------------------------------------
    ! Subroutine used to print the cooling fluxes as a function of
    ! the temperature to a file. Only useful for debugging purposes
    ! -----------------------------------------------------------------   
    subroutine debug_cooling_fluxes() 
      
      use read_input_module, only : heat_cooling_debug, &
                                    sw_magnetic_inclination
      
      integer :: i
      real(amrex_real) :: dT
      real(amrex_real) :: temp
      real(amrex_real) :: q_therm
      real(amrex_real) :: q_vap
      real(amrex_real) :: q_rad
      
      dT = (heat_cooling_debug(3) - heat_cooling_debug(2))/nint(heat_cooling_debug(4))
      temp = heat_cooling_debug(2)
      
      open (2, file = 'cooling_fluxes.dat', status = 'unknown')
      write(2, *) '# plasma flux and magnetic inclination: ', heat_cooling_debug(5), sw_magnetic_inclination
      write(2, *) '# Temperature[K], Thermionic flux [W/m^2], Radiative flux [W/m^2], Vaporization flux [W/m^2]'
      do i = 0,nint(heat_cooling_debug(4)) 
         call thermionic_cooling(temp, heat_cooling_debug(5), q_therm)
         call vaporization_cooling(temp, q_vap)
         call radiation_cooling(temp, q_rad)
         write(2,*) temp, q_therm, q_vap, q_rad
         temp = temp + dT
      end do
      close(2) 
      
    end subroutine debug_cooling_fluxes
    
   ! -----------------------------------------------------------------
   ! Given an array xx(1:n), and given a value x, returns a value j 
   ! such that x is between xx(j) and xx(j+1). xx(1:n) must be 
   ! monotonic, either increasing or decreasing. j=0 or j=n is 
   ! returned to indicate that x is out of range.
   ! -----------------------------------------------------------------
    subroutine bisection(xx,n,x,j)

      ! Input and output variables
      integer, intent(in) :: n
      real(amrex_real), intent(in) :: xx(n)
      real(amrex_real), intent(in) :: x
      integer, intent(out) :: j

      ! Local variables
      integer jl,jm,ju

      jl=0 ! Initialize lower
      ju=n+1 ! upper limits.
      do while(ju-jl.gt.1)
         jm=(ju+jl)/2
         if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm))) then ! eqv is used so that the subroutine can work with
                                                     ! decreasing and increasing functions.
            jl=jm ! If the guess overshoot the solution, set lower limit to current guess.
         else
            ju=jm ! If the guess undershoot the solution, set lower limit to current guess.
         endif
      end do
      if(x.eq.xx(1)) then
         j=1 ! Treament of edge-case (query point at beginning of values vector)
      else if(x.eq.xx(n)) then
         j=n-1 ! Treament of edge-case (query point at end of values vector)
      else
         j=jl
      endif
   end subroutine bisection

   ! -----------------------------------------------------------------
   ! Bilinear interpolation. Values "val" are ordered so that 
   ! val(1)->val(4) correspond to accesing the grid starting from the 
   ! bottom left corner and moving counter clockwise (x1,y1)->(x2,y1)
   ! (x2->y2)->(x1,y2).
   ! -----------------------------------------------------------------
   subroutine bilin_intrp(x, y, val, xy_query, val_query)
      ! Input and output variables
      real(amrex_real), intent(in) :: x(1:2)
      real(amrex_real), intent(in) :: y(1:2)
      real(amrex_real), intent(in) :: val(1:4)
      real(amrex_real), intent(in) :: xy_query(1:2)
      real(amrex_real), intent(out) :: val_query

      ! Local variables
      real(amrex_real) t, u

      t = (xy_query(1)-x(1))/(x(2)-x(1))
      u = (xy_query(2)-y(1))/(y(2)-y(1))

      val_query = (1-t)*(1-u)*val(1) + t*(1-u)*val(2) &
                  + t*u*val(3) + (1-t)*u*val(4)

   end subroutine bilin_intrp

   ! -----------------------------------------------------------------
   ! Linear interpolation. Values should be passed so that x1->val1
   ! x2->val2.
   ! -----------------------------------------------------------------
   subroutine lin_intrp(x, val, x_query, val_query)
      ! Input and output variables
      real(amrex_real), intent(in) :: x(1:2)
      real(amrex_real), intent(in) :: val(1:2)
      real(amrex_real), intent(in) :: x_query
      real(amrex_real), intent(out) :: val_query

      ! Local variables
      real(amrex_real) t

      t = (x_query-x(1))/(x(2)-x(1))

      val_query = (1-t)*val(1) + t*val(2)

   end subroutine lin_intrp

   ! -----------------------------------------------------------------
   ! Trilinear interpolation. Values "val" are ordered so that 
   ! val(1)->val(8) correspond to accesing the grid starting from the 
   ! bottom left corner, first moving counter clockwise and move in 
   ! the z-direction after completing one circle. (x1,y1,z1)->(x2,y1,z1)
   ! (x2->y2,z1)->(x1,y2,z1)->(x1,y1,z2)->(x2,y1,z2)->(x2->y2,z2)->(x1,y2,z2).
   ! -----------------------------------------------------------------
   subroutine trilin_intrp(x, y, z, val, xyz_query, val_query)
      ! Input and output variables
      real(amrex_real), intent(in) :: x(1:2)
      real(amrex_real), intent(in) :: y(1:2)
      real(amrex_real), intent(in) :: z(1:2)
      real(amrex_real), intent(in) :: val(1:8)
      real(amrex_real), intent(in) :: xyz_query(1:3)
      real(amrex_real), intent(out) :: val_query

      ! Local variables
      real(amrex_real) t, u, w

      t = (xyz_query(1)-x(1))/(x(2)-x(1))
      u = (xyz_query(2)-y(1))/(y(2)-y(1))
      w = (xyz_query(3)-z(1))/(z(2)-z(1))

      val_query = (1-t)*(1-u)*(1-w)*val(1) + t*(1-u)*(1-w)*val(2) &
                  + t*u*(1-w)*val(3) + (1-t)*u*(1-w)*val(4) &
                  + (1-t)*(1-u)*w*val(5) + t*(1-u)*w*val(6) &
                  + t*u*w*val(7) + (1-t)*u*w*val(8)

   end subroutine trilin_intrp


   ! -----------------------------------------------------------------
   ! Subroutine that returns the convection coefficient given a 
   ! temperature. The values re obtained by performing a seventh order
   ! polyonymal fit to 30 extracted values of figure 3.11 of "Impact 
   ! of geometry and shaping of the plasma facing components on hot spot
   ! generation in tokamak devices", Alex GROSJEAN"
   ! -----------------------------------------------------------------
   subroutine get_convection_coeff(T, h)
      ! Input and output variables
      real(amrex_real), intent(in) :: T
      real(amrex_real), intent(out) :: h

      ! Local variables
      real(amrex_real) :: T0, Tc
      real(amrex_real) :: p1, p2, p3, p4, p5, p6, p7, p8

      T0 = 273.15
      p1 = 1.61185E-10
      p2 = -1.5378146E-7
      p3 = 5.927767E-5
      p4 = -0.0118479
      p5 = 1.311232
      p6 = -79.753407
      p7 = 2599.511
      p8 = 3.45107679E4

      Tc = T-T0
      if (Tc.gt.294.0) then
         write(*,*) 'Temperature near pipe requires extrapolation when computing the convection coefficient'
         write(*,*) 'The temperature adjacent to pipe is (in Kelvin)'
         write(*,*) T
         write(*,*) 'The convection coefficient will be arbitratly assumed constant for T > 567 K'
         Tc = 294
      end if

      h = p1*Tc**7 + p2*Tc**6 + p3*Tc**5 + p4*Tc**4 + p5*Tc**3 + p6*Tc**2 + p7*Tc + p8
   end subroutine get_convection_coeff

end module heat_flux_module
