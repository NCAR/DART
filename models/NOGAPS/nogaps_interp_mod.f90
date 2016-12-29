! This code may (or may not) be part of the NOGAPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

module nogaps_interp_mod

    use types_mod, only: r8,            &
                         MISSING_R8
    use utilities_mod, only: E_ERR,     &
                             error_handler

implicit none
private


! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"



    public :: compute_neighbors
    public :: NUM_NEIGHBORS
    public :: get_val_at_pressure

    ! Use bilinear interpolation
    integer, parameter :: NUM_NEIGHBORS        = 4
    integer, parameter :: NEIGHBOR_LOWER_LEFT  = 1
    integer, parameter :: NEIGHBOR_UPPER_LEFT  = 2
    integer, parameter :: NEIGHBOR_UPPER_RIGHT = 3
    integer, parameter :: NEIGHBOR_LOWER_RIGHT = 4

    integer, parameter :: OUT_OF_BOUNDS_ERROR = 13

    character(len=128) :: msgbuf    ! for error handler


!------------------------------------------------------------------
contains
!------------------------------------------------------------------


    subroutine compute_neighbors(target_lat, target_lon, lats, lons, nlats, &
                                 nlons, neighbor_i, neighbor_j,             &
                                 neighbor_weights)
        real(kind=r8), intent(in)  :: target_lat
        real(kind=r8), intent(in)  :: target_lon
        real(kind=r8), intent(in)  :: lats(nlats)
        real(kind=r8), intent(in)  :: lons(nlons)
        integer,       intent(in)  :: nlats
        integer,       intent(in)  :: nlons
        integer,       intent(out) :: neighbor_i(NUM_NEIGHBORS)
        integer,       intent(out) :: neighbor_j(NUM_NEIGHBORS)
        real(kind=r8), intent(out) :: neighbor_weights(NUM_NEIGHBORS)

        real(kind=r8) :: lat_frac, lon_frac

        integer :: i_left, i_right, j_up, j_down
        integer :: error_status

        call get_bounds(target_lat, lats, nlats, j_up, j_down, lat_frac)
        call get_bounds(target_lon, lons, nlons, i_left, i_right, lon_frac, &
                        is_cyclic=.true.)

        ! Ordering goes clockwise starting at lower left 
        neighbor_i(NEIGHBOR_LOWER_LEFT)  = i_left
        neighbor_j(NEIGHBOR_LOWER_LEFT)  = j_down
        neighbor_i(NEIGHBOR_UPPER_LEFT)  = i_left
        neighbor_j(NEIGHBOR_UPPER_LEFT)  = j_up
        neighbor_i(NEIGHBOR_UPPER_RIGHT) = i_right
        neighbor_j(NEIGHBOR_UPPER_RIGHT) = j_up
        neighbor_i(NEIGHBOR_LOWER_RIGHT) = i_right
        neighbor_j(NEIGHBOR_LOWER_RIGHT) = j_down

        call compute_weights(lon_frac, lat_frac, neighbor_weights)
    end subroutine compute_neighbors


!------------------------------------------------------------------


    subroutine get_bounds(point_val, coord_vals, num_coords, upper_bound, &
                          lower_bound, delta, is_cyclic)
        real(kind=r8),     intent(in)  :: point_val
        real(kind=r8),     intent(in)  :: coord_vals(num_coords)
        integer,           intent(in)  :: num_coords           
        integer,           intent(out) :: upper_bound
        integer,           intent(out) :: lower_bound
        real(kind=r8),     intent(out) :: delta
        logical, optional, intent(in)  :: is_cyclic

        real(kind=r8), allocatable :: temp_coords(:)
        real(kind=r8)              :: step_size
        logical                    :: use_cyclic

        integer :: coord_index

        if (present(is_cyclic)) then
            use_cyclic = is_cyclic
        else
            use_cyclic = .false.
        end if

        if (use_cyclic) then
            allocate(temp_coords(num_coords+2))
            step_size                   = coord_vals(2) - coord_vals(1)
            temp_coords(1)              = coord_vals(1) - step_size
            temp_coords(num_coords + 2) = coord_vals(num_coords) + step_size
            temp_coords(2:num_coords+1) = coord_vals(:)
        else
            allocate(temp_coords(num_coords))
            temp_coords(:) = coord_vals(:)
        end if

        upper_bound = -1
        lower_bound = -1
        do coord_index = 2, size(temp_coords)
            if (point_val <= temp_coords(coord_index)) then
                lower_bound = coord_index - 1
                upper_bound = coord_index
                delta       = (point_val - temp_coords(lower_bound)) / &
                              (temp_coords(upper_bound) -              &
                               temp_coords(lower_bound))
                
                if (use_cyclic) then
                    ! Handle edge cases
                    if (upper_bound == num_coords+2) then
                        ! Ran off the right edge
                        upper_bound = 1
                        lower_bound = num_coords
                    elseif (lower_bound == 1) then
                        ! Ran off the left edge
                        lower_bound = num_coords
                        upper_bound = 1
                    else
                        upper_bound = upper_bound - 1
                        lower_bound = lower_bound - 1
                    end if
                end if
                return
            end if
        end do

        write(msgbuf, "(A,F12.8)") 'Ran off the end of the loop ', point_val
        !print *, temp_coords(144:)
        call error_handler(E_ERR, 'get_bounds', msgbuf, "", "", "")
    end subroutine get_bounds


!------------------------------------------------------------------


    subroutine compute_weights(lon_frac, lat_frac, weights)
        real(kind=r8), intent(in)  :: lon_frac
        real(kind=r8), intent(in)  :: lat_frac
        real(kind=r8), intent(out) :: weights(NUM_NEIGHBORS)

        weights(NEIGHBOR_LOWER_LEFT)  = (1 - lon_frac) * (1 - lat_frac)
        weights(NEIGHBOR_UPPER_LEFT)  = (1 - lon_frac) * (lat_frac)
        weights(NEIGHBOR_UPPER_RIGHT) = (lon_frac)     * (lat_frac)
        weights(NEIGHBOR_LOWER_RIGHT) = (lon_frac)     * (1 - lat_frac)
    end subroutine compute_weights 


!------------------------------------------------------------------


    subroutine get_val_at_pressure(input_values, input_pressures, &
                                   interp_pressure, interp_value, &
                                   error_status)
        real(kind=r8), intent(in)  :: input_values(:)
        real(kind=r8), intent(in)  :: input_pressures(:)
        real(kind=r8), intent(in)  :: interp_pressure
        real(kind=r8), intent(out) :: interp_value
        integer,       intent(out) :: error_status

        integer :: input_count
        integer :: cur_level

        real(kind=r8) :: pressure_above, pressure_below
        real(kind=r8) :: value_above, value_below
        real(kind=r8) :: interp_weight

        ! First check to make sure this will actually work
        if (interp_pressure > maxval(input_pressures) .or. &
            interp_pressure < minval(input_pressures)) then
            print *, 'Not in appropriate pressure bounds.'
            print *, interp_pressure, maxval(input_pressures), minval(input_pressures)
            error_status = OUT_OF_BOUNDS_ERROR
            return
        end if

        ! Note that here, "above" and "below" are in a *geometric*
        ! sense, i.e. pressure_below > pressure_above.
        input_count = size(input_values)
        do cur_level = 2, input_count
            if (interp_pressure < input_pressures(cur_level)) then
                pressure_below = input_pressures(cur_level)     
                pressure_above = input_pressures(cur_level - 1)     
                value_below    = input_values(cur_level)     
                value_above    = input_values(cur_level - 1)     

                interp_weight = (pressure_below - interp_pressure) / &
                                (pressure_below - pressure_above)
                interp_value  = (value_above * (    interp_weight)) + &
                                (value_below * (1 - interp_weight))
                error_status = 0
                return
            end if    
        end do

        ! Should not get here!
        interp_value = MISSING_R8
        error_status = 99
    end subroutine get_val_at_pressure 

!===================================================================
! End of nogaps_interp_mod
!===================================================================

end module nogaps_interp_mod


! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

