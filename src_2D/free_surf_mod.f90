module free_surface_module
  
    ! -----------------------------------------------------------------
    ! This module is used to read a file that descibes the heat flux 
    ! imposed on the free surface
    ! -----------------------------------------------------------------
    
    use amrex_amr_module

    use read_files_module

    implicit none

    private

    public :: get_free_surface
    public :: get_surface_normal
 
   contains
 
    ! -----------------------------------------------------------------
    ! Subroutine used to construct the table containing the heat flux
    ! read from file
    ! -----------------------------------------------------------------
    subroutine get_free_surface(surf_ind, surf_xlo, surf_dx, melt_top, melt_pos, surf_pos, surf_normal, surf_normal_undeformed)
    
        use read_input_module, only: sw_free_surface_file
        
        implicit none
        ! Input and output variables
        integer, intent(in) :: surf_ind(1:1,1:2)
        real(amrex_real), intent(in) :: surf_xlo(1)
        real(amrex_real), intent(in) :: surf_dx(1)
        real(amrex_real), intent(out) :: melt_top(surf_ind(1,1):surf_ind(1,2))
        real(amrex_real), intent(out) :: melt_pos(surf_ind(1,1):surf_ind(1,2))
        real(amrex_real), intent(out) :: surf_pos(surf_ind(1,1):surf_ind(1,2))
        real(amrex_real), intent(out) :: surf_normal(surf_ind(1,1):surf_ind(1,2),1:2)
        real(amrex_real), intent(out) :: surf_normal_undeformed(surf_ind(1,1):surf_ind(1,2),1:2)

        ! Local variables
        integer :: points(1:1)
        real(amrex_real), dimension(:), allocatable :: x_file
        real(amrex_real), dimension(:), allocatable :: y_file
        real(amrex_real) :: x_query
        real(amrex_real) :: x_adj(1:2)
        real(amrex_real) :: y_adj(1:2)
        real(amrex_real) :: y_answer
        integer :: i
        integer :: iq
        
        call get_mesh_dimensions (sw_free_surface_file, 1, len(sw_free_surface_file), points)

        allocate (x_file(1:points(1)))
        allocate (y_file(1:points(1)))
        
        call read_free_surface_file(sw_free_surface_file, len(sw_free_surface_file), points(1), x_file, y_file)

        do i = surf_ind(1,1), surf_ind(1,2)

            x_query = surf_xlo(1) + (i-surf_ind(1,1)+0.5)*surf_dx(1)
            call bisection(x_file, points(1), x_query, iq)
            ! Treament of edge case - make sure that you
            ! dont do extrapolation
            if(iq.eq.0) then 
                x_query = x_file(1) 
                iq = 1
            elseif(iq.eq.points(1)) then
                x_query = x_file(points(1))
                iq = points(1)-1
            end if
            x_adj(1) = x_file(iq)
            x_adj(2) = x_file(iq+1)

            y_adj(1) = y_file(iq)
            y_adj(2) = y_file(iq+1)

            call lin_intrp(x_adj, y_adj, x_query, y_answer)
            surf_pos(i) = y_answer
            melt_pos(i) = y_answer
            melt_top(i) = y_answer
        end do

        call get_surface_normal(surf_ind, surf_dx, surf_pos, surf_normal)
        call get_surface_normal(surf_ind, surf_dx, surf_pos, surf_normal_undeformed)
    
    end subroutine get_free_surface

    ! -----------------------------------------------------------------
    ! Subroutine used to calculate the local surface normals
    ! -----------------------------------------------------------------
    subroutine get_surface_normal(surf_ind, surf_dx, surf_pos, surf_normal)
        ! use amr_data_module, only : surf_normal, &
        !                             surf_pos, &
        !                             surf_ind, &
        !                             surf_dx

        ! Input and output variables
        integer, intent(in) :: surf_ind(1:1,1:2)
        real(amrex_real), intent(in) :: surf_dx(1)
        real(amrex_real), intent(in) :: surf_pos(surf_ind(1,1):surf_ind(1,2))
        real(amrex_real), intent(out) :: surf_normal(surf_ind(1,1):surf_ind(1,2),1:2)
        ! ! Local variables
        real(amrex_real) :: db1dx(surf_ind(1,1):surf_ind(1,2))
        integer :: i
        real(amrex_real) :: normalization_factor

        db1dx = 0.0

        ! Derivatives along x-direction
        ! For edge points use non-central difference
        i =  surf_ind(1,1)
        db1dx(i) = (surf_pos(i+1) - surf_pos(i))/surf_dx(1)
        i =  surf_ind(1,2)
        db1dx(i) = (surf_pos(i) - surf_pos(i-1))/surf_dx(1)
        ! For the rest of the points use central 
        do i = surf_ind(1,1)+1,surf_ind(1,2)-1
            db1dx(i) = (surf_pos(i+1) - surf_pos(i-1))/(2*surf_dx(1))
        end do

        ! Update the angle of the surface normals
        do i = surf_ind(1,1),surf_ind(1,2)
            normalization_factor = sqrt(1+db1dx(i)**2)
            surf_normal(i,1) = -db1dx(i)/normalization_factor
            surf_normal(i,2) = 1/normalization_factor
        end do

    end subroutine get_surface_normal


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
     
end module free_surface_module