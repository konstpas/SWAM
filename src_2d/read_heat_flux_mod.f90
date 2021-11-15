module read_heat_flux_module
  
   ! -----------------------------------------------------------------
   ! This module is used to read a file that descibes the heat flux 
   ! imposed on the free surface
   ! -----------------------------------------------------------------
   
    use amrex_amr_module

    implicit none

    private

    public :: get_mesh_dimensions
    public :: read_heatflux_file

    contains

    ! -----------------------------------------------------------------
    ! Subroutine used to get the dimensions of the cartesian mesh 
    ! described in the heat flux input file.
    ! -----------------------------------------------------------------
    subroutine get_mesh_dimensions(input_filename, dims)

        ! Input and output variable
        character ( len = 80 ), intent(in) :: input_filename
        integer, intent(out) :: dims(2)

        ! Local variables
        logical got_one
        character( len = 255 ) line
        integer ierror
        integer input_unit, input_status

        ierror = 0
        
        call get_unit (input_unit)

        open ( unit = input_unit, file = input_filename, status = 'old', &
        iostat = input_status )
    
        if ( input_status /= 0 ) then
            ierror = 1
            stop 'Could not open input file'
        end if
    
        got_one = .false.
        do while (.not.got_one) ! Go through the document until you find a non
                                ! empty line that is not a comment line

            read ( input_unit, '(a)', iostat = input_status ) line ! Reads one line
        
            if ( input_status /= 0 ) then
                stop 'No data lines given'
            end if
        
            if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
                cycle ! If the line starts with #, then it's a comment line, start over the do loop.
            end if  
        
            call s_to_i4vec ( line, 2, dims, ierror)  
            if ( ierror /= 0 ) then
                cycle
            end if
            got_one = .true. ! The dimensions of the mesh ("dims") have been found.
        end do
    
        close ( unit = input_unit )
    end subroutine get_mesh_dimensions

    ! -----------------------------------------------------------------
    ! Subroutine used to read the heat flux described in an input.
    ! It also returns arrays for the vectors that generated the heat
    ! flux mesh in the input file.
    ! -----------------------------------------------------------------
    subroutine read_heatflux_file(input_filename, tpoints, &
                                xpoints, heatflux)

        ! Input and output variables
        character ( len = 80 ), intent(in) :: input_filename 
        real(amrex_real), allocatable, dimension(:), intent(out) :: tpoints
        real(amrex_real), allocatable, dimension(:), intent(out)  :: xpoints
        real(amrex_real), allocatable, dimension(:,:), intent(out)  :: heatflux
        ! Local variables
        logical got_one
        character ( len = 255 ) line
        integer :: dims(2)
        integer ierror
        integer input_unit, input_status
        integer i

        ierror = 0
        
        call get_unit (input_unit)

        open ( unit = input_unit, file = input_filename, status = 'old', &
        iostat = input_status )
    
        if ( input_status /= 0 ) then
            ierror = 1
            stop 'Could not open input file'
        end if
    
    
        got_one = .false.
        do while (.not.got_one) ! Go through the document until you find a non
                                ! empty line that is not a comment line

            read ( input_unit, '(a)', iostat = input_status ) line ! Reads one line
        
            if ( input_status /= 0 ) then
                stop 'No data given.'
            end if
        
            if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
                cycle ! If the line starts with #, then it's a comment line, start over the do loop.
            end if  
        
            call s_to_i4vec ( line, 2, dims, ierror)  
            if ( ierror /= 0 ) then
                cycle
            end if
            got_one = .true. ! The dimensions of the mesh ("dims") have been found.
        end do

        allocate (heatflux(1:dims(1),1:dims(2)))
        allocate (tpoints(1:dims(1)))
        allocate (xpoints(1:dims(2)))

        call r8vec_data_read (input_unit, input_status, dims(1), tpoints) ! Reads dims(1) amount of lines and finds the 
                                                                          ! vector that describes the time coordinates
        call r8vec_data_read (input_unit, input_status, dims(2), xpoints) ! Reads dims(2) amount of lines and finds the 
                                                                          ! vector that describes the space coordinates

        do i = 1, dims(1) ! Fills the heatflux table, line by line (one line corresponds to the heat-flux at all locations
                          ! at a given time.)
            call r8vec_data_read (input_unit, input_status, dims(2), heatflux(i,1:dims(2)))
        end do

    
        close ( unit = input_unit )
    end subroutine read_heatflux_file

    ! -----------------------------------------------------------------
    ! Subroutine used to find an available i/o unit. Taken from 
    ! https://people.sc.fsu.edu/~jburkardt/f_src/table_io/table_io.html
    ! -----------------------------------------------------------------
    subroutine get_unit ( iunit )

    !*****************************************************************************80
    !
    !! GET_UNIT returns a free FORTRAN unit number.
    !
    !  Discussion:
    !
    !    A "free" FORTRAN unit number is a value between 1 and 99 which
    !    is not currently associated with an I/O device.  A free FORTRAN unit
    !    number is needed in order to open a file with the OPEN command.
    !
    !    If IUNIT = 0, then no free FORTRAN unit could be found, although
    !    all 99 units were checked (except for units 5, 6 and 9, which
    !    are commonly reserved for console I/O).
    !
    !    Otherwise, IUNIT is a value between 1 and 99, representing a
    !    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
    !    are special, and will never return those values.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    26 October 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Output, integer IUNIT, the free unit number.
    !
        implicit none

        integer i
        integer ios
        integer iunit
        logical lopen

        iunit = 0

        do i = 1, 99

            if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

                inquire ( unit = i, opened = lopen, iostat = ios )

                if ( ios == 0 ) then
                if ( .not. lopen ) then
                    iunit = i
                    return
                end if
                end if

            end if

        end do
    end subroutine get_unit

    ! -----------------------------------------------------------------
    ! Subroutine used  through a string "s" and returns an integer(kind=4)
    ! https://people.sc.fsu.edu/~jburkardt/f_src/table_io/table_io.html
    ! -----------------------------------------------------------------
    subroutine s_to_i4 ( s, ival, ierror, length )

    !*****************************************************************************80
    !
    !! S_TO_I4 reads an I4 from a string.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    28 June 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) S, a string to be examined.
    !
    !    Output, integer IVAL, the integer value read from the string.
    !    If the string is blank, then IVAL will be returned 0.
    !
    !    Output, integer IERROR, an error flag.
    !    0, no error.
    !    1, an error occurred.
    !
    !    Output, integer LENGTH, the number of characters of S
    !    used to make IVAL.
    !
        implicit none

        character c
        integer i
        integer ierror
        integer isgn
        integer istate
        integer ival
        integer length
        character ( len = * ) s

        ierror = 0
        istate = 0
        isgn = 1
        ival = 0

        do i = 1, len_trim ( s )

            c = s(i:i)
        !
        !  Haven't read anything.
        !
            if ( istate == 0 ) then

                if ( c == ' ' ) then

                else if ( c == '-' ) then
                    istate = 1
                    isgn = -1
                else if ( c == '+' ) then
                    istate = 1
                    isgn = + 1
                else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
                    istate = 2
                    ival = ichar ( c ) - ichar ( '0' )
                else
                    ierror = 1
                    return
                end if
        !
        !  Have read the sign, expecting digits.
        !
            else if ( istate == 1 ) then

                if ( c == ' ' ) then

                else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
                    istate = 2
                    ival = ichar ( c ) - ichar ( '0' )
                else
                    ierror = 1
                    return
                end if
        !
        !  Have read at least one digit, expecting more.
        !
            else if ( istate == 2 ) then

                if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
                    ival = 10 * ival + ichar ( c ) - ichar ( '0' )
                else
                    ival = isgn * ival
                    length = i - 1
                    return
                end if

            end if

        end do
    !
    !  If we read all the characters in the string, see if we're OK.
    !
        if ( istate == 2 ) then
            ival = isgn * ival
            length = len_trim ( s )
        else
            ierror = 1
            length = 0
        end if

    end subroutine s_to_i4
    
    ! -----------------------------------------------------------------
    ! Goes through a string "s" and returns "n" an integers(kind=4).
    ! https://people.sc.fsu.edu/~jburkardt/f_src/table_io/table_io.html
    ! -----------------------------------------------------------------
    subroutine s_to_i4vec ( s, n, ivec, ierror )

    !*****************************************************************************80
    !
    !! S_TO_I4VEC reads an I4VEC from a string.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    08 October 2003
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) S, the string to be read.
    !
    !    Input, integer N, the number of values expected.
    !
    !    Output, integer IVEC(N), the values read from the string.
    !
    !    Output, integer IERROR, error flag.
    !    0, no errors occurred.
    !    -K, could not read data for entries -K through N.
    !
        implicit none

        integer n

        integer i
        integer ierror
        integer ilo
        integer ivec(n)
        integer length
        character ( len = * ) s

        i = 0
        ierror = 0
        ilo = 1
        length = 0

        do while ( i < n )

            i = i + 1

            call s_to_i4 ( s(ilo:), ivec(i), ierror, length )

            if ( ierror /= 0 ) then
                ierror = -i
                exit
            end if

            ilo = ilo + length

        end do

    end subroutine s_to_i4vec
 
    ! -----------------------------------------------------------------
    ! Subroutine that takes an open file and reads "n" real numbers 
    ! from in. Taken from 
    ! https://people.sc.fsu.edu/~jburkardt/f_src/table_io/table_io.html
    ! -----------------------------------------------------------------             
    subroutine r8vec_data_read (input_unit, input_status, n, table )

    !*****************************************************************************80
    !
    !! R8VEC_DATA_READ reads data from an R8VEC file.
    !
    !  Discussion:
    !
    !    An R8VEC is a vector of R8 values.
    !
    !  Discussion:
    !
    !    The file may contain more than N points, but this routine will
    !    return after reading N of them.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    10 June 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
    !
    !    Input, integer N, the number of points.
    !
    !    Output, real ( kind = rk ) TABLE(N), the data.
    !
        implicit none

        integer n

        integer ierror
        integer input_unit
        integer input_status
        integer j
        integer length
        character ( len = 255 ) line
        real (amrex_real) table(n)
        real (amrex_real) x

        ierror = 0


        if ( input_status /= 0 ) then
            stop 'Empty file passed.'
        end if

        j = 0

        do while ( j < n )

            read ( input_unit, '(a)', iostat = input_status ) line

            if ( input_status /= 0 ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'R8VEC_DATA_READ - Fatal error!'
                write ( *, '(a)' ) '  Error while reading lines of data.'
                write ( *, '(a,i8)' ) '  Number of data lines read, J =         ', j
                write ( *, '(a,i8)' ) '  Number of data lines needed, N =       ', n
                stop 1
            end if

            if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
                cycle
            end if

            call s_to_r8 ( line, x, ierror, length )

            if ( ierror /= 0 ) then
                cycle
            end if

            j = j + 1

            table(j) = x

        end do

    end subroutine r8vec_data_read

    ! -----------------------------------------------------------------
    ! Subroutine that takes a string "s" and reads a real number 
    ! from in. Taken from 
    ! https://people.sc.fsu.edu/~jburkardt/f_src/table_io/table_io.html
    ! -----------------------------------------------------------------
    subroutine s_to_r8 ( s, dval, ierror, length )

    !*****************************************************************************80
    !
    !! S_TO_R8 reads an R8 from a string.
    !
    !  Discussion:
    !
    !    The routine will read as many characters as possible until it reaches
    !    the end of the string, or encounters a character which cannot be
    !    part of the number.
    !
    !    Legal input is:
    !
    !       1 blanks,
    !       2 '+' or '-' sign,
    !       2.5 blanks
    !       3 integer part,
    !       4 decimal point,
    !       5 fraction part,
    !       6 'E' or 'e' or 'D' or 'd', exponent marker,
    !       7 exponent sign,
    !       8 exponent integer part,
    !       9 exponent decimal point,
    !      10 exponent fraction part,
    !      11 blanks,
    !      12 final comma or semicolon,
    !
    !    with most quantities optional.
    !
    !  Example:
    !
    !    S                 DVAL
    !
    !    '1'               1.0
    !    '     1   '       1.0
    !    '1A'              1.0
    !    '12,34,56'        12.0
    !    '  34 7'          34.0
    !    '-1E2ABCD'        -100.0
    !    '-1X2ABCD'        -1.0
    !    ' 2E-1'           0.2
    !    '23.45'           23.45
    !    '-4.2E+2'         -420.0
    !    '17d2'            1700.0
    !    '-14e-2'         -0.14
    !    'e2'              100.0
    !    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    07 September 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) S, the string containing the
    !    data to be read.  Reading will begin at position 1 and
    !    terminate at the end of the string, or when no more
    !    characters can be read to form a legal real.  Blanks,
    !    commas, or other nonnumeric data will, in particular,
    !    cause the conversion to halt.
    !
    !    Output, real ( kind = rk ) DVAL, the value read from the string.
    !
    !    Output, integer IERROR, error flag.
    !    0, no errors occurred.
    !    1, 2, 6 or 7, the input number was garbled.  The
    !    value of IERROR is the last type of input successfully
    !    read.  For instance, 1 means initial blanks, 2 means
    !    a plus or minus sign, and so on.
    !
    !    Output, integer LENGTH, the number of characters read
    !    to form the number, including any terminating
    !    characters such as a trailing comma or blanks.
    !
        implicit none

        character c
        ! logical ch_eqi
        real (amrex_real) dval
        integer ierror
        integer ihave
        integer isgn
        integer iterm
        integer jbot
        integer jsgn
        integer jtop
        integer length
        integer nchar
        integer ndig
        real (amrex_real) rbot
        real (amrex_real) rexp
        real (amrex_real) rtop
        character ( len = * ) s

        nchar = len_trim ( s )

        ierror = 0
        dval = 0.0D+00
        length = -1
        isgn = 1
        rtop = 0
        rbot = 1
        jsgn = 1
        jtop = 0
        jbot = 1
        ihave = 1
        iterm = 0

        do

            length = length + 1

            if ( nchar < length+1 ) then
                exit
            end if

            c = s(length+1:length+1)
        !
        !  Blank character.
        !
            if ( c == ' ' ) then

                if ( ihave == 2 ) then

                else if ( ihave == 6 .or. ihave == 7 ) then
                iterm = 1
                else if ( 1 < ihave ) then
                ihave = 11
                end if
        !
        !  Comma.
        !
            else if ( c == ',' .or. c == ';' ) then

                if ( ihave /= 1 ) then
                iterm = 1
                ihave = 12
                length = length + 1
                end if
        !
        !  Minus sign.
        !
            else if ( c == '-' ) then

                if ( ihave == 1 ) then
                    ihave = 2
                    isgn = -1
                else if ( ihave == 6 ) then
                    ihave = 7
                    jsgn = -1
                else
                    iterm = 1
                end if
        !
        !  Plus sign.
        !
            else if ( c == '+' ) then

                if ( ihave == 1 ) then
                    ihave = 2
                else if ( ihave == 6 ) then
                    ihave = 7
                else
                    iterm = 1
                end if
        !
        !  Decimal point.
        !
            else if ( c == '.' ) then

                if ( ihave < 4 ) then
                    ihave = 4
                else if ( 6 <= ihave .and. ihave <= 8 ) then
                    ihave = 9
                else
                    iterm = 1
                end if
        !
        !  Scientific notation exponent marker.
        !
            else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

                if ( ihave < 6 ) then
                    ihave = 6
                else
                    iterm = 1
                end if
        !
        !  Digit.
        !
            else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

                if ( ihave <= 2 ) then
                    ihave = 3
                else if ( ihave == 4 ) then
                    ihave = 5
                else if ( ihave == 6 .or. ihave == 7 ) then
                    ihave = 8
                else if ( ihave == 9 ) then
                    ihave = 10
                end if

                call ch_to_digit ( c, ndig )

                if ( ihave == 3 ) then
                    rtop = 10.0D+00 * rtop + real ( ndig, kind = amrex_real)
                else if ( ihave == 5 ) then
                    rtop = 10.0D+00 * rtop + real ( ndig, kind = amrex_real)
                    rbot = 10.0D+00 * rbot
                else if ( ihave == 8 ) then
                    jtop = 10 * jtop + ndig
                else if ( ihave == 10 ) then
                    jtop = 10 * jtop + ndig
                    jbot = 10 * jbot
                end if
        !
        !  Anything else is regarded as a terminator.
        !
            else
                iterm = 1
            end if
        !
        !  If we haven't seen a terminator, and we haven't examined the
        !  entire string, go get the next character.
        !
            if ( iterm == 1 ) then
                exit
            end if

        end do
    !
    !  If we haven't seen a terminator, and we have examined the
    !  entire string, then we're done, and LENGTH is equal to NCHAR.
    !
        if ( iterm /= 1 .and. length + 1 == nchar ) then
            length = nchar
        end if
    !
    !  Number seems to have terminated.  Have we got a legal number?
    !  Not if we terminated in states 1, 2, 6 or 7!
    !
        if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
            ierror = ihave
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
            write ( *, '(a)' ) '  Illegal or nonnumeric input:'
            write ( *, '(a)' ) '    ' // trim ( s )
            return
        end if
    !
    !  Number seems OK.  Form it.
    !
        if ( jtop == 0 ) then
            rexp = 1.0D+00
        else
        if ( jbot == 1 ) then
            rexp = 10.0D+00 ** ( jsgn * jtop )
        else
            rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = amrex_real) &
            / real ( jbot, kind = amrex_real ) )
        end if
        end if

        dval = real ( isgn, kind = amrex_real) * rexp * rtop / rbot

        return
    end subroutine s_to_r8

    ! -----------------------------------------------------------------
    ! A functions that takes two characters and check if they are 
    ! equal (case insensitive). Taken from
    ! https://people.sc.fsu.edu/~jburkardt/f_src/table_io/table_io.html
    ! -----------------------------------------------------------------
    function ch_eqi ( c1, c2 )

    !*****************************************************************************80
    !
    !! CH_EQI is a case insensitive comparison of two characters for equality.
    !
    !  Example:
    !
    !    CH_EQI ( 'A', 'a' ) is .TRUE.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    28 July 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character C1, C2, the characters to compare.
    !
    !    Output, logical CH_EQI, the result of the comparison.
    !
        implicit none

        logical ch_eqi
        character c1
        character c1_cap
        character c2
        character c2_cap

        c1_cap = c1
        c2_cap = c2

        call ch_cap ( c1_cap )
        call ch_cap ( c2_cap )

        if ( c1_cap == c2_cap ) then
            ch_eqi = .true.
        else
            ch_eqi = .false.
        end if

    end function ch_eqi

    ! -----------------------------------------------------------------
    ! A subroutine that takes a character and returns an integer with
    ! the same value (if character is an integer). Taken from
    ! https://people.sc.fsu.edu/~jburkardt/f_src/table_io/table_io.html
    ! -----------------------------------------------------------------
    subroutine ch_to_digit ( c, digit )

    !*****************************************************************************80
    !
    !! CH_TO_DIGIT returns the integer value of a base 10 digit.
    !
    !  Example:
    !
    !     C   DIGIT
    !    ---  -----
    !    '0'    0
    !    '1'    1
    !    ...  ...
    !    '9'    9
    !    ' '    0
    !    'X'   -1
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    04 August 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character C, the decimal digit, '0' through '9' or blank
    !    are legal.
    !
    !    Output, integer DIGIT, the corresponding integer value.
    !    If C was 'illegal', then DIGIT is -1.
    !
        implicit none

        character c
        integer digit

        if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

            digit = ichar ( c ) - 48

        else if ( c == ' ' ) then

            digit = 0

        else

            digit = -1

        end if

    end subroutine ch_to_digit
    
    ! -----------------------------------------------------------------
    ! A subroutine that takes a character and capcapitalizes it. 
    ! Taken from
    ! https://people.sc.fsu.edu/~jburkardt/f_src/table_io/table_io.html
    ! -----------------------------------------------------------------
    subroutine ch_cap ( c )

    !*****************************************************************************80
    !
    !! CH_CAP capitalizes a single character.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    19 July 1998
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input/output, character C, the character to capitalize.
    !
        implicit none

        character c
        integer itemp

        itemp = ichar ( c )

        if ( 97 <= itemp .and. itemp <= 122 ) then
        c = char ( itemp - 32 )
        end if

    end subroutine ch_cap
end module read_heat_flux_module