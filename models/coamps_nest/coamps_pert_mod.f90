! This code may (or may not) be part of the GITM distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

!------------------------------
! MODULE:       coamps_pert_mod
! AUTHOR:       T. R. Whitcomb
!               Naval Research Laboratory
! DART VERSION: Jamaica
!
! Module for perturbing 
!------------------------------ 

module coamps_pert_mod

    use coamps_domain_mod, only : coamps_domain, get_domain_nest

    use coamps_statevec_mod, only : state_vector, state_iterator, &
                                    get_iterator, has_next, get_next, &
                                    get_domain

    use coamps_statevar_mod, only : state_variable, &
                                    get_var_substate, &
                                    get_pert_type,    &
                                    get_pert_magnitude, &
                                    PERT_TYPE_NOPERTS,    &
                                    PERT_TYPE_UNIFORM,    &
                                    PERT_TYPE_INDIVID,    &
                                    get_nest_number

    use coamps_nest_mod, only : coamps_nest, &
                                nest_index_2d_to_1d, &
                                get_nest_i_width, &
                                get_nesT_j_width

    use mpi_utilities_mod, only : my_task_id

    use random_seq_mod, only: random_seq_type, &
                              init_random_seq, &
                              random_gaussian

    use types_mod, only : r8

    implicit none

    private

    !------------------------------
    ! BEGIN PUBLIC INTERFACE
    !------------------------------

    public :: perturb_state

    !------------------------------
    ! END PUBLIC INTERFACE
    !------------------------------

    !------------------------------
    ! BEGIN EXTERNAL INTERFACES
    !------------------------------
    ! [none]
    !------------------------------
    ! END EXTERNAL INTERFACES
    !------------------------------

    !------------------------------
    ! BEGIN TYPES AND CONSTANTS 
    !------------------------------
    ! [none]
    !------------------------------
    ! END TYPES AND CONSTANTS 
    !------------------------------

    !------------------------------
    ! BEGIN MODULE VARIABLES
    !------------------------------

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

    logical, save :: module_initialized = .false.

    real(kind=r8), save :: ZERO_MEAN = real(0.0, kind=r8)

    integer, save :: x_bound_skip = 0
    integer, save :: y_bound_skip = 0

    type(random_seq_type), save :: random_sequence
    type(coamps_domain)         :: domain

    !------------------------------
    ! END MODULE VARIABLES
    !------------------------------

contains

    !------------------------------
    ! BEGIN PUBLIC ROUTINES
    !------------------------------

    ! Returns a perturbed state vector based on the perturbations to 
    ! component variables defined in the state vector definition
    !  PARAMETERS
    !   IN  state           Complete state vector to perturb
    !   IN  state_layout    State vector layout
    !   IN  i_bound_skip    How many points in i to skip along the boundary
    !   IN  j_bound_skip    How many points in j to skip along the boundary
    function perturb_state(state, state_layout, i_bound_skip, j_bound_skip)
        real(kind=r8), dimension(:), intent(in)  :: state
        type(state_vector)  ,        intent(in)  :: state_layout
        integer,                     intent(in)  :: i_bound_skip
        integer,                     intent(in)  :: j_bound_skip
        real(kind=r8), dimension(:), pointer     :: perturb_state

        type(state_iterator) :: iterator

        if (.not. module_initialized) then
            call initialize_module(i_bound_skip, j_bound_skip, &
                                   get_domain(state_layout))
        end if

        iterator = get_iterator(state_layout)
        do while (has_next(iterator))
            call perturb_var(get_next(iterator), state, perturb_state)
        end do
    end function perturb_state

    !------------------------------
    ! END PUBLIC ROUTINES
    !------------------------------

    !------------------------------
    ! BEGIN PRIVATE ROUTINES
    !------------------------------

    ! perturb_var
    ! -----------
    ! Given a state variable and the state vector, generates a perturbed
    ! portion of the state vector corresponding to the variable
    !  PARAMETERS
    !   IN  var             The state variable to perturb
    !   IN  state           The (full) state to perturb 
    ! INOUT perturb_state   The state including this variable's perturbations
    subroutine perturb_var(var, state, perturb_state)
        type(state_variable),                intent(in)    :: var
        real(kind=r8), dimension(:), target, intent(in)    :: state
        real(kind=r8), dimension(:), target, intent(inout) :: perturb_state

        real(kind=r8), dimension(:), pointer :: state_subsect
        real(kind=r8), dimension(:), pointer :: pert_state_subsect

        state_subsect      => get_var_substate(var, state)
        pert_state_subsect => get_var_substate(var, perturb_state)

        ! Need to reset the perturbation generator for each variable since it
        ! is using possibly a different perturbation magnitude
        select case (get_pert_type(var))
            case (PERT_TYPE_NOPERTS)
                call null_pert()
                pert_state_subsect = apply_perturbation(var, state_subsect)
            case (PERT_TYPE_UNIFORM)
                call uniform_pert(magnitude=real(get_pert_magnitude(var), kind=r8))
                pert_state_subsect = apply_perturbation(var, state_subsect)
            case (PERT_TYPE_INDIVID)
                call individual_pert(magnitude=real(get_pert_magnitude(var), kind=r8))
                pert_state_subsect = apply_perturbation(var, state_subsect)
        end select

    end subroutine perturb_var

    ! apply_perturbation
    ! -----------------------
    ! Perturbs the given vector with a Gaussian perturbation    
    ! Note that the perturbation uses a static real size - this is because
    ! I can't get the interface statement to recognize the "r8" type
    !  PARAMETERS
    !   IN  var             State variable to perturb
    !   IN  state           Vector to perturb
    function apply_perturbation(cur_var, state)
        type(state_variable),        intent(in)  :: cur_var
        real(kind=r8), dimension(:), intent(in)  :: state
        integer                                     :: pert_type
        real(kind=r8), dimension(:), pointer        :: apply_perturbation
        type(coamps_nest)                           :: var_nest
        real(kind=r8)                               :: pert

        ! Loop over the interior of the domain, away from the boundary
        integer :: ii, jj, nn
        integer :: interior_i_start, interior_i_end
        integer :: interior_j_start, interior_j_end

        var_nest = get_domain_nest(domain, get_nest_number(cur_var))
        interior_i_start = 1 + x_bound_skip
        interior_j_start = 1 + y_bound_skip
        interior_i_end   = get_nest_i_width(var_nest) - x_bound_skip
        interior_j_end   = get_nest_j_width(var_nest) - y_bound_skip

        pert_type = get_pert_type(cur_var)

        do jj = interior_j_start, interior_j_end 
            do ii = interior_i_start, interior_i_end
                nn = nest_index_2d_to_1d(var_nest, ii, jj)

                select case (pert_type)
                  case (PERT_TYPE_NOPERTS)
                    call null_pert()
                  case (PERT_TYPE_UNIFORM)
                    call uniform_pert()
                  case (PERT_TYPE_INDIVID)
                    call individual_pert()
                end select

                apply_perturbation(nn) = state(nn) + pert
            end do
        end do
        print *, "Total perturbation: "
        print *, size(apply_perturbation), size(state)
        print *, "", sum(abs(apply_perturbation(1:nn) - state(1:nn)))
    end function apply_perturbation

    ! initialize_perturbation
    ! -----------------------
    ! Set up module to generate random perturbations
    !  PARAMETERS
    !   [none]
    subroutine initialize_module(i_bound_skip, j_bound_skip, state_domain)
        integer,             intent(in) :: i_bound_skip
        integer,             intent(in) :: j_bound_skip
        type(coamps_domain), intent(in) :: state_domain

        integer               :: rand_seq_seed
        integer, parameter    :: DT_MILLISECOND = 8
        integer, dimension(8) :: datetime 

        ! The boundary points to skip is static 
        x_bound_skip = i_bound_skip
        y_bound_skip = j_bound_skip
    
        ! Seed the random number generator based on the number mmmP where
        ! mmm is the current millisecond component of the time and P is 
        ! the current process ID
        call DATE_AND_TIME(values=datetime)
        rand_seq_seed = 1E1 * datetime(DT_MILLISECOND) + my_task_id()
        
        call init_random_seq(r=random_sequence, seed=rand_seq_seed)

        ! Need the domain to get the proper nest sizes
        domain = state_domain

        module_initialized = .true.
        
    end subroutine initialize_module

    !---
    ! Perturbation Methods
    !---
    ! Emulate closure-like behavior so we can easily add new perturbation
    ! methods.  Each subroutine should be set up with a first (optional)
    ! argument called "value" that will return a perturbation.  A second
    ! optional argument gives the magnitude of the perturbation.  The
    ! presence of the "magnitude" argument is interpreted as a "reset" 
    ! command.

    ! null_pert
    ! ---------
    ! Always return zero
    subroutine null_pert(value, magnitude)
        real(kind=r8), optional, intent(out) :: value
        real(kind=r8), optional, intent(in)  :: magnitude

        if (present(value)) then
            value = real(0, kind=r8)
        end if
    end subroutine null_pert

    ! uniform_pert
    ! ------------
    ! Return the same value until reset with a (possibly different) 
    ! perturbation magnitude
    subroutine uniform_pert(value, magnitude)
        real(kind=r8), optional :: value
        real(kind=r8), optional :: magnitude

        real(kind=r8), save :: pert_value

        if (present(magnitude)) then
            pert_value = real(random_gaussian(random_sequence, ZERO_MEAN, magnitude))
        end if

        if (present(value)) then
            value = pert_value
        end if
    end subroutine uniform_pert

    ! individual_pert
    ! ---------------
    ! Always return a different value consistent with the last-supplied
    ! perturbation magnitude
    subroutine individual_pert(value, magnitude)
        real(kind=r8), optional, intent(out) :: value
        real(kind=r8), optional, intent(in)  :: magnitude

        real(kind=r8), save :: pert_magnitude

        if (present(magnitude)) then
            pert_magnitude = magnitude
        end if

        if (present(value)) then
            value = real(random_gaussian(random_sequence, ZERO_MEAN, pert_magnitude))
        end if
    end subroutine individual_pert

    !------------------------------
    ! END PRIVATE ROUTINES
    !------------------------------

end module coamps_pert_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
