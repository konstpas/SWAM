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
        integer, intent(in) :: surf_ind(2,2)
        real(amrex_real), intent(in) :: surf_xlo(2)
        real(amrex_real), intent(in) :: surf_dx(2)
        real(amrex_real), intent(out) :: melt_top(surf_ind(1,1):surf_ind(1,2), surf_ind(2,1):surf_ind(2,2))
        real(amrex_real), intent(out) :: melt_pos(surf_ind(1,1):surf_ind(1,2), surf_ind(2,1):surf_ind(2,2))
        real(amrex_real), intent(out) :: surf_pos(surf_ind(1,1):surf_ind(1,2), surf_ind(2,1):surf_ind(2,2))
        real(amrex_real), intent(out) :: surf_normal(surf_ind(1,1):surf_ind(1,2), surf_ind(2,1):surf_ind(2,2),1:3)
        real(amrex_real), intent(out) :: surf_normal_undeformed(surf_ind(1,1):surf_ind(1,2), surf_ind(2,1):surf_ind(2,2),1:3)
    
        ! Local variables
        integer :: points(1:2)
        real(amrex_real), dimension(:), allocatable :: x_file
        real(amrex_real), dimension(:), allocatable :: y_file
        real(amrex_real), dimension(:), allocatable :: z_file
        real(amrex_real) :: xz_query(1:2)
        real(amrex_real) :: x_adj(1:2)
        real(amrex_real) :: z_adj(1:2)
        real(amrex_real) :: y_adj(1:4)
        real(amrex_real) :: y_answer
        integer :: i, j
        integer :: iq, jq
       
        call get_mesh_dimensions (sw_free_surface_file, 2, len(sw_free_surface_file), points(1:2))

        allocate (x_file(1:points(1)))
        allocate (y_file(1:points(1)*points(2)))
        allocate (z_file(1:points(2)))
        
        call read_free_surface_file(sw_free_surface_file, len(sw_free_surface_file), points, x_file, y_file, z_file)

        do i = surf_ind(1,1), surf_ind(1,2)

            xz_query(1) = surf_xlo(1) + (i-surf_ind(1,1)+0.5)*surf_dx(1)
            call bisection(x_file, points(1), xz_query(1), iq)
            ! Treament of edge case - make sure that you
            ! dont do extrapolation
            if(iq.eq.0) then 
                xz_query(1) = x_file(1) 
                iq = 1
            elseif(iq.eq.points(1)) then
                xz_query(1) = x_file(points(1))
                iq = points(1)-1
            end if
            x_adj(1) = x_file(iq)
            x_adj(2) = x_file(iq+1)

            do j = surf_ind(2,1), surf_ind(2,2)

                xz_query(2) = surf_xlo(2) + (j-surf_ind(2,1)+0.5)*surf_dx(2)
                call bisection(z_file, points(2), xz_query(2), jq)
                ! Treament of edge case - make sure that you
                ! dont do extrapolation
                if(jq.eq.0) then 
                    xz_query(2) = z_file(1) 
                    jq = 1
                elseif(jq.eq.points(2)) then
                    xz_query(2) = z_file(points(2))
                    jq = points(2)-1
                end if
                z_adj(1) = z_file(jq)
                z_adj(2) = z_file(jq+1)

                y_adj(1) = y_file((iq-1)*points(2)+jq)
                y_adj(2) = y_file(iq*points(2)+jq)
                y_adj(3) = y_file(iq*points(2)+jq+1)
                y_adj(4) = y_file((iq-1)*points(2)+jq+1)

                call bilin_intrp(x_adj, z_adj, y_adj, xz_query, y_answer)
                surf_pos(i,j) = y_answer
                melt_pos(i,j) = y_answer
                melt_top(i,j) = y_answer
            end do
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
           integer, intent(in) :: surf_ind(2,2)
           real(amrex_real), intent(in) :: surf_dx(2)
           real(amrex_real), intent(in) :: surf_pos(surf_ind(1,1):surf_ind(1,2), surf_ind(2,1):surf_ind(2,2))
           real(amrex_real), intent(out) :: surf_normal(surf_ind(1,1):surf_ind(1,2), surf_ind(2,1):surf_ind(2,2),1:3)
         ! ! Local variables
         real(amrex_real) :: db1dx(surf_ind(1,1):surf_ind(1,2),surf_ind(2,1):surf_ind(2,2))
         real(amrex_real) :: db1dz(surf_ind(1,1):surf_ind(1,2),surf_ind(2,1):surf_ind(2,2))
         integer :: i,j
         real(amrex_real) :: normalization_factor
   
         db1dx = 0.0
         db1dz = 0.0
   
         ! Derivatives along x-direction
         ! For edge points use non-central difference
         do j = surf_ind(2,1),surf_ind(2,2)
            i =  surf_ind(1,1)
            db1dx(i,j) = (surf_pos(i+1,j) - surf_pos(i,j))/surf_dx(1)
            i =  surf_ind(1,2)
            db1dx(i,j) = (surf_pos(i,j) - surf_pos(i-1,j))/surf_dx(1)
         end do
         ! For the rest of the points use central 
         do i = surf_ind(1,1)+1,surf_ind(1,2)-1
            do j = surf_ind(2,1),surf_ind(2,2)
               db1dx(i,j) = (surf_pos(i+1,j) - surf_pos(i-1,j))/(2*surf_dx(1))
            end do
         end do
   
         ! Derivatives along z-direction
         ! For edge points use non-central difference
         do i = surf_ind(1,1),surf_ind(1,2)
            j = surf_ind(2,1)
            db1dz(i,j) = (surf_pos(i,j+1) - surf_pos(i,j))/surf_dx(2)
            j = surf_ind(2,2)
            db1dz(i,j) = (surf_pos(i,j) - surf_pos(i,j-1))/surf_dx(2)
         end do
         ! For the rest of the points use central 
         do i = surf_ind(1,1),surf_ind(1,2)
            do j = surf_ind(2,1)+1,surf_ind(2,2)-1
               db1dz(i,j) = (surf_pos(i,j+1) - surf_pos(i,j-1))/(2*surf_dx(2))
            end do
         end do
   
   
         ! Update the angle of the surface normals
         do i = surf_ind(1,1),surf_ind(1,2)
            do j = surf_ind(2,1),surf_ind(2,2)
               normalization_factor = sqrt(1+db1dx(i,j)**2+db1dz(i,j)**2)
               surf_normal(i,j,1) = -db1dx(i,j)/normalization_factor
               surf_normal(i,j,2) = 1/normalization_factor
               surf_normal(i,j,3) = -db1dz(i,j)/normalization_factor
            end do
         end do
   
      end subroutine get_surface_normal

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