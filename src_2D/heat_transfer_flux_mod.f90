module heat_transfer_flux_module
  
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
                                    temp, t_lo, t_hi, lev, qb)
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)  
    integer, intent(in) :: id_lo(2), id_hi(2)
    integer, intent(in) :: t_lo(2), t_hi(2)
    integer, intent(in) :: lev
    real(amrex_real), intent(in) :: time
    real(amrex_real), intent(in) :: xlo(2)
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2))
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2))
    real(amrex_real), intent(out) :: qb(lo(1):hi(1),lo(2):hi(2))
    
    ! Local variables
    real(amrex_real) :: q_free(lo(1):hi(1),lo(2):hi(2))
    real(amrex_real) :: q_active(lo(1):hi(1),lo(2):hi(2))
    real(amrex_real) :: q_side(lo(1):hi(1),lo(2):hi(2))

    ! Initialize heat flux
    qb = 0.0

    ! Free surface contribution
    call get_free_surface_heat_flux(time, xlo, &
                                    dx, lo, hi, &
                                    idom, id_lo, id_hi, &
                                    temp, t_lo, t_hi, &
                                    lev, q_free)

    ! Active cooling contribution
    call get_active_cooling_heat_flux(dx, lo, hi, &
                                      idom, id_lo, id_hi, &
                                      temp, t_lo, t_hi, &
                                      q_active) 
    
    ! Second exposed side contribution
    call get_exposed_side_heat_flux(time, xlo, &
                                    dx, lo, hi, &
                                    idom, id_lo, id_hi, &
                                    temp, t_lo, t_hi, &
                                    q_side) 

    ! Sum all contributions
    qb = q_free + q_active + q_side
    
  end subroutine get_boundary_heat_flux
  

  ! -----------------------------------------------------------------
  ! Subroutine used to prescribe the boundary heating on the free
  ! surface. Note that the incoming heat flux is assigned to the
  ! first node under the free surface in the form of a volumetric
  ! heating (after an appropriate dimensionality correction)
  ! -----------------------------------------------------------------   
  subroutine get_free_surface_heat_flux(time, xlo, &
                                        dx, lo, hi, &
                                        idom, id_lo, id_hi, &
                                        temp, t_lo, t_hi, &
                                        lev, qb) 
    
    use read_input_module, only : plasma_flux_type, &
                                  cooling_thermionic, &
                                  cooling_vaporization, &
                                  cooling_radiation, &
                                  cooling_radiation
    use amr_data_module, only : surf_current
    
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)  
    integer, intent(in) :: id_lo(2), id_hi(2)
    integer, intent(in) :: t_lo(2), t_hi(2)
    real(amrex_real), intent(in) :: time
    real(amrex_real), intent(in) :: xlo(2)
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2))
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2))
    integer, intent(in) :: lev
    real(amrex_real), intent(out) :: qb(lo(1):hi(1),lo(2):hi(2))
    
    ! Local variables
    integer :: i, j
    real(amrex_real) :: q_plasma
    real(amrex_real) :: q_rad
    real(amrex_real) :: q_therm
    real(amrex_real) :: q_vap
    real(amrex_real) :: xpos, ypos
    logical :: side_flag

    side_flag = .false.

    ! Initialize heat flux
    qb = 0.0

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
       
          ! Locate the free surface 
          if(nint(idom(i,j)).ne.0 .and. nint(idom(i,j+1)).eq.0) then
             
             ! Initialize heat flux contributions
             q_plasma = 0.0
             q_therm = 0.0
             q_vap = 0.0
             q_rad = 0.0
             
             ! Location of the free surface
             xpos = xlo(1) + (i - lo(1))*dx(1)
             ypos = xlo(2) + (j - lo(2))*dx(2)
             
             ! Plasma flux
             if (plasma_flux_type.eq.'Gaussian') then
                call gaussian_heat_flux(time, xpos, side_flag, q_plasma)
             elseif (plasma_flux_type.eq.'Uniform') then
                call uniform_heat_flux(time, xpos, side_flag, q_plasma)
             elseif (plasma_flux_type.eq.'Input_file') then
                call file_heat_flux (time, xpos, side_flag, q_plasma)
             else
                STOP "Unknown plasma heat flux type"
             end if
             
             ! Thermionic cooling flux (at the maximum level, store the
             ! value of the thermionic current needed in the shallow water solver)
             if (cooling_thermionic) then
                if(lev.eq.amrex_max_level) then
                   call thermionic_cooling(temp(i,j), q_plasma, q_therm, surf_current(i))
                else
                   call thermionic_cooling(temp(i,j), q_plasma, q_therm)
                end if
             end if
             
             ! Vaporization cooling flux
             if (cooling_vaporization) then
                call vaporization_cooling(temp(i ,j), q_vap)
             end if
             
             ! Radiative cooling flux
             if (cooling_radiation) then
                call radiation_cooling(temp(i,j), q_rad)
             end if
             
             ! Sum all flux contributions
             qb(i,j) = q_plasma - q_rad - q_vap - q_therm
             
             ! Note: the term /dx(2) converts a surface heat flux [W/m^2]
             ! into a volumetric heat flux [W/m^3]
             qb(i,j) = qb(i,j)/dx(2)
             
             ! Check the validity of the heat flux
             if(qb(i,j).ne.qb(i,j)) then
                
                print *, 'Nan heat flux on free surface at cell ' 
                print *, i, j
                print *, xpos, ypos
                print *, 'Where the surface temperature is '
                print *, temp(i,j)
                STOP
                
             end if
             
          end if
                    
       end do
    end do
    
  end subroutine get_free_surface_heat_flux

  ! -----------------------------------------------------------------
  ! Subroutine used to prescribe the heat flux obtained from active
  ! cooling via a cooling pipe
  ! -----------------------------------------------------------------   
  subroutine get_active_cooling_heat_flux(dx, lo, hi, &
                                          idom, id_lo, id_hi, &
                                          temp, t_lo, t_hi, qb)
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)  
    integer, intent(in) :: id_lo(2), id_hi(2)
    integer, intent(in) :: t_lo(2), t_hi(2)
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2))
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2))
    real(amrex_real), intent(out) :: qb(lo(1):hi(1),lo(2):hi(2))
    
    ! Local variables
    integer :: i, j
    real(amrex_real) :: q_cool

    ! Initialize heat flux
    qb = 0.0   

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          ! Get contributions from the different sides of the cooling pipe
          
          if(nint(idom(i,j)).ne.-1 .and. nint(idom(i,j+1)).eq.-1) then
             call active_cooling(temp(i,j), q_cool)
             qb(i,j) = -q_cool/dx(2)
          end if

          if(nint(idom(i,j)).ne.-1 .and. nint(idom(i,j-1)).eq.-1) then
             call active_cooling(temp(i,j), q_cool)
             qb(i,j) = -q_cool/dx(2)
          end if

          if(nint(idom(i,j)).ne.-1 .and. nint(idom(i+1,j)).eq.-1) then
             call active_cooling(temp(i,j), q_cool)
             qb(i,j) = -q_cool/dx(1)
          end if

          if(nint(idom(i,j)).ne.-1 .and. nint(idom(i-1,j)).eq.-1) then
             call active_cooling(temp(i,j), q_cool)
             qb(i,j) = -q_cool/dx(1)
          end if

                    
       end do
    end do
    
  end subroutine get_active_cooling_heat_flux


  ! -----------------------------------------------------------------
  ! Subroutine used to prescribe the heat  flux resulting from 
  ! a second side exposed to the plasma flux
  ! -----------------------------------------------------------------   
  subroutine get_exposed_side_heat_flux(time, xlo, &
                                        dx, lo, hi, &
                                        idom, id_lo, id_hi, &
                                        temp, t_lo, t_hi, qb) 
    
    use read_input_module, only : plasma_flux_side_type, &
                                  cooling_thermionic, &
                                  cooling_vaporization, &
                                  cooling_radiation, &
                                  sample_edge, &
                                  geom_name, & 
                                  cooling_radiation
    
    
    ! Input and output variables
    integer, intent(in) :: lo(2), hi(2)  
    integer, intent(in) :: id_lo(2), id_hi(2)
    integer, intent(in) :: t_lo(2), t_hi(2)
    real(amrex_real), intent(in) :: time
    real(amrex_real), intent(in) :: xlo(2)
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(in) :: idom(id_lo(1):id_hi(1),id_lo(2):id_hi(2))
    real(amrex_real), intent(in) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2))
    real(amrex_real), intent(out) :: qb(lo(1):hi(1),lo(2):hi(2))
    
    ! Local variables
    integer :: i, j
    real(amrex_real) :: q_plasma
    real(amrex_real) :: q_rad
    real(amrex_real) :: q_therm
    real(amrex_real) :: q_vap
    real(amrex_real) :: xpos, ypos
    logical :: side_flag

    side_flag = .true.
    
    ! Initialize heat flux to zero
    qb = 0.0

    ! Return immediately if the geometry does not include a second exposed side
    if (geom_name .ne. "West" .and. geom_name .ne. "West_rectangular") return

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          
          if( nint(idom(i,j)).ne.0 .and. & 
              (xpos.ge.sample_edge .or. nint(idom(i+1,j)).eq.0)) then

             ! Initialize heat flux contributions
             q_plasma = 0.0
             q_therm = 0.0
             q_vap = 0.0
             q_rad = 0.0

             ! Location of the exposed surface
             xpos = xlo(1) + (1 + i - lo(1))*dx(1)
             ypos = xlo(2) + (j - lo(2))*dx(2)
             
             ! Plasma flux
             if (plasma_flux_side_type.eq.'Gaussian') then
                call gaussian_heat_flux(time, ypos, side_flag, q_plasma)
             elseif (plasma_flux_side_type.eq.'Uniform') then
                call uniform_heat_flux(time, ypos, side_flag, q_plasma)
             elseif (plasma_flux_side_type.eq.'Input_file') then
                call file_heat_flux (time, xpos, side_flag, q_plasma)
             else
                STOP "Unknown plasma heat flux type"
             end if
             
             ! Thermionic cooling flux
             if (cooling_thermionic) then
                call thermionic_cooling(temp(i,j), q_plasma, q_therm)
             end if
             
             ! Vaporization cooling flux
             if (cooling_vaporization) then
                call vaporization_cooling(temp(i ,j), q_vap)
             end if
             
             ! Radiative cooling flux
             if (cooling_radiation) then
                call radiation_cooling(temp(i,j), q_rad)
             end if


             ! Sum all flux contributions
             qb(i,j) = q_plasma - q_rad - q_vap - q_therm
             
             ! Note: the term /dx(2) converts a surface heat flux [W/m^2]
             ! into a volumetric heat flux [W/m^3]
             qb(i,j) = qb(i,j)/dx(1)
             
             ! Check the validity of the heat flux
             if(qb(i,j).ne.qb(i,j)) then
           
                print *, 'Nan heat flux on second exposed face at cell ' 
                print *, i, j
                print *, xpos, ypos
                print *, 'Where the surface temperature is '
                print *, temp(i,j)
                STOP
                
             end if

          end if
          
       end do
    end do
    
  end subroutine get_exposed_side_heat_flux
  
  ! -----------------------------------------------------------------
  ! Subroutine used to prescribe a gaussian heat flux active for
  ! time <= time_exposure
  ! -----------------------------------------------------------------   
  subroutine gaussian_heat_flux(time, xpos, side_flag, qb) 
    
    use read_input_module, only : plasma_flux_params, &
                                  plasma_flux_side_params
    
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
       plasma_params = plasma_flux_side_params
    else
       plasma_params = plasma_flux_params
    end if
    
    if (time.ge.plasma_params(1) .and. time.le.plasma_params(2)) then
       qb = plasma_params(3) &
            *EXP(-((xpos-plasma_params(4))**2)/(plasma_params(5)**2))
    end if
    
  end subroutine gaussian_heat_flux
  
  ! -----------------------------------------------------------------
  ! Subroutine used to prescribe a uniform heat flux active for
  ! time <= time_exposure
  ! -----------------------------------------------------------------   
  subroutine uniform_heat_flux(time, xpos, side_flag, qb) 
    
    use read_input_module, only : plasma_flux_params, &
                                  plasma_flux_side_params
    
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
       plasma_params = plasma_flux_side_params
    else
       plasma_params = plasma_flux_params
    end if
    
    if (time.ge.plasma_params(1) .and. time.le.plasma_params(2)) then
       if(xpos.ge.plasma_params(4) .and. xpos.le.plasma_params(5)) then
          qb = plasma_params(3) 
       end if
    end if
    
  end subroutine uniform_heat_flux
    
  ! -----------------------------------------------------------------
  ! Subroutine used to prescribe a heat flux defined in
  ! an input file.
  ! -----------------------------------------------------------------   
  subroutine file_heat_flux(time, xpos, side_flag, qb) 

    use amr_data_module, only : plasma_flux_surf_mesh, &
                                plasma_flux_time_mesh, &
                                plasma_flux_table, &
                                plasma_flux_side_surf_mesh, &
                                plasma_flux_side_time_mesh, &
                                plasma_flux_side_table
    
    
    ! Input and output variables
    real(amrex_real), intent(in) :: time
    real(amrex_real), intent(in) :: xpos
    logical, intent(in) :: side_flag
    real(amrex_real), intent(out) :: qb

    ! Local variables 
    integer :: i_x, i_t, n, m
    real(amrex_real) :: x(2), t(2), val(4)
    real(amrex_real) :: tx_query(2)
    real(amrex_real), allocatable :: spatial_mesh(:)
    real(amrex_real), allocatable :: temporal_mesh(:)
    real(amrex_real), allocatable :: heat_flux(:,:)
    

    ! Initialize heat flux
    qb = 0_amrex_real
    
    if(side_flag) then
       n = size(plasma_flux_side_time_mesh,1)
       m = size(plasma_flux_side_surf_mesh,1)
       allocate (heat_flux(1:n,1:m))
       allocate (temporal_mesh(1:n))
       allocate (spatial_mesh(1:m))
       spatial_mesh = plasma_flux_side_surf_mesh
       temporal_mesh = plasma_flux_side_time_mesh
       heat_flux = plasma_flux_side_table
    else
       n = size(plasma_flux_time_mesh,1)
       m = size(plasma_flux_surf_mesh,1)
       allocate (heat_flux(1:n,1:m))
       allocate (temporal_mesh(1:n))
       allocate (spatial_mesh(1:m))
       spatial_mesh = plasma_flux_surf_mesh
       temporal_mesh = plasma_flux_time_mesh
       heat_flux = plasma_flux_table
    end if
    
    ! Find the maximum index i_t such that the
    ! time falls in-between plasma_flux_time_mesh(i_t)
    ! and plasma_flux_time_mesh(i_t+1). Similar for i_x
    call bisection(temporal_mesh, n, time, i_t)
    call bisection(spatial_mesh, m, xpos, i_x)
    
    ! If query point is outside the 4 bounds
    ! then take the closest corner
    if(i_t.eq.0 .and. i_x.eq.0) then
       qb = heat_flux(1,1)
    elseif (i_t.eq.n .and. i_x.eq.m) then
       qb = heat_flux(n,m)
    elseif (i_t.eq.0 .and. i_x.eq.m) then
       qb = heat_flux(1,m)
    elseif (i_t.eq.n .and. i_x.eq.0) then
       qb = heat_flux(n,1)
       ! If query point is outside 2 bounds
       ! then linear interpolation
    elseif (i_t.eq.0 .or. i_t.eq.n) then
       x(1) = spatial_mesh(i_x)
       x(2) = spatial_mesh(i_x+1)
       if(i_t.eq.0) i_t = 1
       val(1) = heat_flux(i_t, i_x)
       val(2) = heat_flux(i_t, i_x+1)
       call lin_intrp(x, val(1:2), xpos, qb)
    elseif (i_x.eq.0 .or. i_x.eq.m) then
       t(1) = temporal_mesh(i_t)
       t(2) = temporal_mesh(i_t+1)
       if(i_x.eq.0) i_x = 1
       val(1) = heat_flux(i_t, i_x)
       val(2) = heat_flux(i_t+1, i_x)
       call lin_intrp(t, val(1:2), time, qb)
    else
       t(1) = temporal_mesh(i_t)
       t(2) = temporal_mesh(i_t+1)
       x(1) = spatial_mesh(i_x)
       x(2) = spatial_mesh(i_x+1)
       val(1) = heat_flux(i_t,i_x)
       val(2) = heat_flux(i_t+1,i_x)
       val(3) = heat_flux(i_t+1,i_x+1)
       val(4) = heat_flux(i_t,i_x+1)
       tx_query(1) = time
       tx_query(2) = xpos
       call bilin_intrp(t, x, val, tx_query, qb)
    endif
    
    ! Free memory
    deallocate(spatial_mesh)
    deallocate(temporal_mesh)
    deallocate(heat_flux)
    
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
                                           get_Richardson
    
    use read_input_module, only : sw_magnetic_inclination
    
    ! Input and output variables                                       
    real(amrex_real), intent(in) :: Ts        ! Temperature at the center of cells adjacent to the free surface [K]
    real(amrex_real), intent(in) :: q_plasma  ! Plasma heat flux [K]
    real(amrex_real), intent(out) :: q_therm  ! Flux of energy due to thermionic emission [W/m^2]
    real(amrex_real), intent(out), optional :: Jth ! Thermionic current [A/m^2]
    
     ! Local variables
    real(amrex_real) :: kb = 1.38064852E-23 ! Boltzmann constant [m^2*kg/(s^2*K)]
    real(amrex_real) :: Jth_nom
    real(amrex_real) :: Jth_lim
    real(amrex_real) :: J
    real(amrex_real) :: Aeff
    real(amrex_real) :: Wf
    real(amrex_real) :: e = 1.60217662E-19
    real(amrex_real) :: pi = 3.1415927
    
    call get_work_function(Wf)
    call get_Richardson(Aeff)
    
    ! Nominal thermionic current from the Richardson-Dushman formula
    Jth_nom = Aeff*EXP(-Wf/(kb*Ts))*Ts**2
    
    ! Space-charge limited current (semi-empirical expression)
    Jth_lim = 1.51e4 * q_plasma**(1.0/3.0) * (SIN(sw_magnetic_inclination/180*pi))**2
    
    ! Minimum between nominal and space-charge limited
    J = MIN(Jth_lim, Jth_nom)
    
    ! Heat flux
    q_therm = J/e*(Wf+2*kb*Ts)
    
    ! Thermionic current
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
    real(amrex_real), parameter :: pi = 3.1415927
    
    T0 = 343
    
    ! Convection coefficient
    call get_convection_coeff(T, h)
    
    ! Heuristic correction to take into account the overestimation of
    ! the heat flux due to piecewise approximation of the circular
    ! cross-section of the pipe
    h = h*2*pi/8

    ! Cooling flux
    q_cool = h*(T-T0)
    
  end subroutine active_cooling
  
  ! -----------------------------------------------------------------
  ! Subroutine that returns the convection coefficient given a 
  ! temperature. The values are obtained by performing a seventh order
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
       print *, "Warning: Convection coefficient assume constant for T > 567 K"
       Tc = 294
    end if
    
    h = p1*Tc**7 + p2*Tc**6 + p3*Tc**5 + p4*Tc**4 + &
         p5*Tc**3 + p6*Tc**2 + p7*Tc + p8
    
  end subroutine get_convection_coeff
    
  ! -----------------------------------------------------------------
  ! Subroutine used to print the cooling fluxes as a function of
  ! the temperature to a file. Only useful for debugging purposes
  ! -----------------------------------------------------------------   
  subroutine debug_cooling_fluxes() 
    
    use read_input_module, only : cooling_debug, &
                                  sw_magnetic_inclination
    
    integer :: i
    real(amrex_real) :: dT
    real(amrex_real) :: temp
    real(amrex_real) :: q_therm
    real(amrex_real) :: q_vap
    real(amrex_real) :: q_rad
    
    dT = (cooling_debug(3) - cooling_debug(2))/nint(cooling_debug(4))
    temp = cooling_debug(2)
    
    open (2, file = 'cooling_fluxes.dat', status = 'unknown')
    write(2, *) '# plasma flux and magnetic inclination: ', cooling_debug(5), sw_magnetic_inclination
    write(2, *) '# Temperature[K], Thermionic flux [W/m^2], Radiative flux [W/m^2], Vaporization flux [W/m^2]'
    do i = 0,nint(cooling_debug(4)) 
       call thermionic_cooling(temp, cooling_debug(5), q_therm)
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

    
end module heat_transfer_flux_module
 
