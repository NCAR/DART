! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

!------------------------------
! MODULE:       coamps_statevar_mod
! AUTHOR:       T. R. Whitcomb
!               Naval Research Laboratory
! DART VERSION: Jamaica
!
! Module containing the data structure and routines for dealing with
! a COAMPS state variable
!------------------------------ 

module coamps_statevar_mod

    use coamps_domain_mod, only : coamps_domain, &
                                  get_domain_num_levels,   &
                                  get_domain_msigma,       &
                                  get_domain_wsigma,       &
                                  get_domain_nest

    use coamps_nest_mod, only : get_nest_size,       &
                                get_nest_i_width,    &
                                get_nest_j_width

    use coamps_util_mod, only : check_io_status,     &
                                uppercase,           &
                                lowercase
    use obs_kind_mod
    use types_mod,       only : r4, r8, missing_r8
    use utilities_mod,   only : E_ERR,         &
                                error_handler
    use location_mod,    only : VERTISUNDEF, VERTISSURFACE, VERTISLEVEL, &
                                VERTISPRESSURE, VERTISHEIGHT


    implicit none

    private

    !------------------------------
    ! BEGIN PUBLIC INTERFACE
    !------------------------------

    public :: state_variable
    public :: new_state_variable

    public :: read_state_variable
    public :: dump_state_variable

    public :: get_record_in_nest

    public :: operator(==)
    ! Query functions
    public :: gets_update
    public :: get_nest_number
    public :: get_state_begin
    public :: get_state_end
    public :: is_var_at_index
    public :: is_nonnegative
    public :: is_sigma_level
    public :: get_vert_loc
    public :: get_vert_type
    public :: get_vert_value
    public :: get_var_name
    public :: get_var_nlevs
    public :: get_var_kind
    public :: get_pert_magnitude
    public :: get_pert_type
    public :: get_mean_flag
    public :: get_io_flag
    public :: get_mass_level_flag
    public :: get_sigma_record
    public :: get_var_max_level
    public :: get_var_min_level
    public :: get_nc_varid
    public :: get_var_stagger

    public :: set_sigma_record
    public :: set_2d_flag
    public :: set_position
    public :: set_nc_varid
    public :: set_var_stagger
    public :: define_mean_var

    ! Array sectioning
    public :: get_var_substate

    ! Perturbation information
    public :: PERT_TYPE_NOPERTS
    public :: PERT_TYPE_UNIFORM
    public :: PERT_TYPE_INDIVID

    !------------------------------
    ! END PUBLIC INTERFACE
    !------------------------------

    !------------------------------
    ! BEGIN EXTERNAL INTERFACE
    !------------------------------

    ! This is the external routine that grabs where fields are in the
    ! restart file - it's generated via shell scripts and is therefore
    ! included as a separate procedure.  Pass in the constants so we
    ! don't need to repeat them inside the subroutine.
    interface
        subroutine get_name_info(DIM_TYPE_2D, DIM_TYPE_3D, DIM_TYPE_3DW,    &
                                 SINGLEIO, MULTIIO, var_name, var_dim_type, &
                                 var_record_num)
        integer,               intent(in)  :: DIM_TYPE_2D
        integer,               intent(in)  :: DIM_TYPE_3D
        integer,               intent(in)  :: DIM_TYPE_3DW
        integer,               intent(in)  :: SINGLEIO
        integer,               intent(in)  :: MULTIIO
        character(len=*),      intent(in)  :: var_name
        integer,               intent(out) :: var_dim_type
        integer, dimension(2), intent(out) :: var_record_num
        end subroutine get_name_info
    end interface

    ! Generic procedure to allow selection of an array subsection
    ! corresponding to a particular state variable - need 2 since there
    ! is potentially a state vector of type r8 or r4
    interface get_var_substate
        module procedure get_var_substate_r8, get_var_substate_r4
    end interface get_var_substate

    ! Check if two state variables are talking about the same thing
    interface operator(==)
        module procedure is_same_variable
    end interface

    !------------------------------
    ! END EXTERNAL INTERFACE
    !------------------------------
  
    !------------------------------
    ! BEGIN TYPES AND CONSTANTS
    !------------------------------

    ! Dimension information - w levels include the surface so there are 
    ! kka+1 levels there
    integer, parameter :: DIM_TYPE_2D  = 1
    integer, parameter :: DIM_TYPE_3D  = 2
    integer, parameter :: DIM_TYPE_3DW = 3

    ! Length of string variables
    integer, parameter :: VAR_NAME_LEN  = 10
    integer, parameter :: VAR_QTY_LEN  = 32
    integer, parameter :: PERT_TYPE_LEN = 7

    ! Perturbation types - don't perturb at all, perturb whole sigma
    ! level field uniformly, or perturb each grid point individually.
    integer,          parameter :: PERT_TYPE_NOPERTS = 0
    integer,          parameter :: PERT_TYPE_UNIFORM = 1
    integer,          parameter :: PERT_TYPE_INDIVID = 2
    character(len=*), parameter :: NOPERTS_NAME      = 'NOPERTS' 
    character(len=*), parameter :: UNIFORM_NAME      = 'UNIFORM'
    character(len=*), parameter :: INDIVID_NAME      = 'INDIVID'

    ! If the variable is defined on a mass level or a w level
    character, parameter :: MASS_LEVEL    = 'M'
    character, parameter :: W_LEVEL       = 'W'
    character, parameter :: ZHT_LEVEL     = 'Z'
    character, parameter :: SFC_LEVEL     = 'S'
    character, parameter :: PRES_LEVEL    = 'P'
    character, parameter :: UNDEF_LEVEL   = 'U'

    logical,   parameter :: IS_MASS_LEVEL = .true.
    logical,   parameter :: IS_W_LEVEL    = .false.
    logical,   parameter :: IS_NOT_SIG    = .true.

    ! Indices for the variable record array to take into account the
    ! difference between the numbering in the restart files written by
    ! a single I/O process and an MPI I/O process
    integer, parameter :: SINGLEIO = 1
    integer, parameter :: MULTIIO  = 2

    ! Whether or not to write the field back to the COAMPS restart file
    integer,          parameter :: FIELD_UPDATE_LEN  = 6
    character(len=*), parameter :: FIELD_UPDATE      = 'UPDATE'
    character(len=*), parameter :: FIELD_FREEZE      = 'FREEZE'
    logical,          parameter :: FLAG_UPDATE_FIELD = .true.
    logical,          parameter :: FLAG_FREEZE_FIELD = .false.

    ! Define the strings used to signify whether we'll check for
    ! positive definiteness or not and their corresponding flags
    integer,          parameter :: POS_DEF_LEN         = 8
    character(len=*), parameter :: FORCE_POSITIVE      = 'ISPOSDEF'
    character(len=*), parameter :: ALLOW_NEGATIVE      = 'NOPOSDEF'
    logical,          parameter :: FLAG_FORCE_POSITIVE = .true.
    logical,          parameter :: FLAG_ALLOW_NEGATIVE = .false.

    ! Define an entry for a single field in the DART state vector.
    ! This encompasses several pieces of data
    !  Where to find the field in the COAMPS restart file:
    !   var_name, dim_type, var_record, sigma_record
    ! How to perturb the variable if necessary:
    !   pert_mag, pert_type
    ! Whether or not we actually assimilate the field (e.g. mean
    ! fields are needed in here to do interpolation but should
    ! not have their assimilated alterations written back to the
    ! COAMPS restart file):  
    !   update_field
    ! Help find specific variables for interpolation:
    !   mean_field, mass_level
    ! Forbid negative values and set anything less than zero to
    ! zero (e.g. for mixing ratios):
    !   positive_definite
    type :: state_variable
        private  

        ! Bounds of this variable within the DART state vector
        integer                     :: state_begin
        integer                     :: state_end
        logical                     :: is_2d_variable

        integer                     :: nest_number          ! IN

        character(len=VAR_NAME_LEN) :: var_name             ! IN
        integer                     :: dim_type             ! IN
        integer, dimension(2)       :: var_record           ! IN
        integer                     :: sigma_record         ! IN
        real(kind=r8)               :: pert_mag             ! IN
        integer                     :: pert_type            ! IN
        logical                     :: update_field         ! IN
        integer                     :: var_kind             ! IN
        logical                     :: mean_field           ! IN
        logical                     :: io_flag              ! IN
        logical                     :: mass_level           ! IN
        logical                     :: positive_definite    ! IN
        integer                     :: nc_varid             ! IN
        character(len=1)            :: var_stagger          ! IN
        integer                     :: vert_type            ! IN 
        real(kind=r8)               :: vert_value           ! IN  
    end type state_variable
  
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
  
    !------------------------------
    ! END MODULE VARIABLES
    !------------------------------

contains

    !------------------------------
    ! BEGIN PUBLIC ROUTINES
    !------------------------------

    ! new_state_variable
    ! ------------------
    ! Given initial parameters for a state variable, consolidates them into a 
    ! "state_variable" structure and returns it.  The non-optional parameters 
    ! are what makes the state variable unique (i.e. only need to check these
    ! for equality tests)
    function new_state_variable(nest_num, var_kind, is_mean_var,            &
                                level_type, sigma_index,                    &
                                pert_mag, pert_type, update_var, is_posdef, &
                                is_twod, var_name, is_flat_file, vert_value)
        integer,                    intent(in)  :: nest_num
        integer,                    intent(in)  :: var_kind
        logical,                    intent(in)  :: is_mean_var
        character(len=*),           intent(in)  :: level_type
        integer,                    intent(in)  :: sigma_index
        real(kind=r8),    optional, intent(in)  :: pert_mag
        integer,          optional, intent(in)  :: pert_type
        logical,          optional, intent(in)  :: update_var
        logical,          optional, intent(in)  :: is_posdef
        logical,          optional, intent(in)  :: is_twod
        character(len=*), optional, intent(in)  :: var_name
        logical,          optional, intent(in)  :: is_flat_file
        real(kind=r8),    optional, intent(in)  :: vert_value
        type(state_variable)                    :: new_state_variable

        new_state_variable%vert_type    = vert_type_from_level_type(level_type)
        new_state_variable%nest_number  = nest_num
        new_state_variable%var_kind     = var_kind
        new_state_variable%mean_field   = is_mean_var
        new_state_variable%io_flag      = .not.is_mean_var
        new_state_variable%mass_level   = mass_flag_from_level_type(level_type)
        new_state_variable%sigma_record = sigma_index

        if (present(pert_mag)) then
            new_state_variable%pert_mag = pert_mag
        end if

        if (present(pert_type)) then
            new_state_variable%pert_type = pert_type
        end if

        if (present(update_var)) then
            new_state_variable%update_field = update_var
        end if

        if (present(is_posdef)) then
            new_state_variable%positive_definite = is_posdef
        end if

        if (present(is_twod)) then
            new_state_variable%is_2d_variable = is_twod
        else
            new_state_variable%is_2d_variable = .true.
        end if

        if (present(vert_value)) then
            new_state_variable%vert_value = vert_value
        else
            new_state_variable%vert_value = missing_r8
        end if

        ! Knowing the name means we can figure out where this variable is in
        ! the big COAMPS restart file
        if (present(var_name)) then
            new_state_variable%var_name = var_name

            if (present(is_flat_file)) then
              if(.not.is_flat_file) then
                call get_name_info(DIM_TYPE_2D, DIM_TYPE_3D, DIM_TYPE_3DW, &
                                   SINGLEIO, MULTIIO, var_name,            & 
                                   new_state_variable%dim_type,            &
                                   new_state_variable%var_record)
              end if
            end if

        end if
    end function new_state_variable

    ! read_state_variable
    ! -------------------
    ! Read in a line from an open state vector definition file and use that
    ! data to populate a state variable
    !  PARAMETERS
    ! INOUT var             State variable to populate
    !   IN  state_def_file  [OPEN] state vector definition file unit
    function read_state_variable(state_def_file, var, nest_number, is_flat_file) result(is_eof)
        type(state_variable), intent(out) :: var
        integer,              intent(in)  :: state_def_file
        integer,              intent(in)  :: nest_number
        logical,              intent(in)  :: is_flat_file
        logical                           :: is_eof

        character(len=VAR_NAME_LEN)      :: variable_name
        integer                          :: sigma_index
        character(len=PERT_TYPE_LEN)     :: pert_type_name
        real(kind=r8)                    :: pert_magnitude
        character                        :: level_type
        character(len=VAR_QTY_LEN)      :: var_kind_name
        character(len=FIELD_UPDATE_LEN)  :: update_type
        logical                          :: is_mean_field
        character(len=POS_DEF_LEN)       :: positive_type_name
        real(kind=r8)                    :: vert_value

        character(len=*), parameter :: routine = 'read_state_variable'
        integer                     :: io_status

        ! Format is:
        !  1. Nest number
        !  2. Variable name - from the model
        !  3. Sigma Level (as an index: 1 - num_levels)
        !  4. Perturbation type (UNIFORM, INDIVID, NOPERTS)
        !  5. Perturbation magnitude
        !  6. Mass/W level (should be either 'M' or 'W')
        !  7. Variable type (from the raw observation kinds)
        !  8. Update field in restart vector or no?
        !  9. Whether field is a mean or not (TRUE/FALSE)
        ! 10. Positive definiteness (ISPOSDEF, NOPOSDEF)
        read(unit=state_def_file, fmt=*, iostat=io_status, end=300)   &
             variable_name, sigma_index, pert_type_name, &
             pert_magnitude, level_type, var_kind_name, update_type,  &
             is_mean_field, positive_type_name

        call check_io_status(io_status, routine, source, revision, &
                             revdate, 'Reading state variable entry')

        !FIXME:  Level values should be read in from the state_def file
        ! Set the level values for 10-m wind and 2-m temperature
        select case (var_kind_num_from_name(var_kind_name))
          case (QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT)
             if(vert_type_from_level_type(level_type) == VERTISHEIGHT) then
               vert_value = 10.0_r8
             else
               vert_value = missing_r8
             end if
               
          case (QTY_TEMPERATURE)
             if(vert_type_from_level_type(level_type) == VERTISHEIGHT) then
               vert_value = 10.0_r8
             else
               vert_value = missing_r8
             end if

          case default
            vert_value = missing_r8
        end select

        ! Some fields need to be converted into values for internal storage
        var = new_state_variable(nest_number,                               &
                                 var_kind_num_from_name(var_kind_name),     &
                                 is_mean_field,                             &
                                 level_type,                                &
                                 sigma_index,                               &
                                 pert_magnitude,                            &
                                 pert_type_from_name(pert_type_name),       &
                                 update_flag_from_type(update_type),        &
                                 posdef_flag_from_name(positive_type_name), &
                                 var_name = variable_name,                  &
                                 vert_value = vert_value,                   &
                                 is_flat_file = is_flat_file)

        is_eof = .false.
        return

300     continue
        is_eof = .true.
    end function read_state_variable

    ! define_mean_var
    ! -------------------
    ! Definitions for mean state variables. 
    !  PARAMETERS
    ! IN variable_name   
    ! IN nest_number
    ! IN current sigma level
    function define_mean_var(variable_name, nest_number, cur_level) result(var)
      character(len=*),  intent(in) :: variable_name 
      integer,           intent(in) :: nest_number
      integer, optional, intent(in) :: cur_level
      type(state_variable)          :: var
       
      logical,       parameter :: is_mean_field  = .true.
      real(kind=r8), parameter :: pert_magnitude = 0.0
      integer                  :: cur_level_local
      character(len=1)         :: level_type
      character(len=30)        :: var_kind_name

      if(present(cur_level)) then
        cur_level_local = cur_level
      else 
        cur_level_local = 0
      end if

      select case (trim(uppercase(variable_name)))
        case ('EXBM')
          level_type      = 'M' 
          var_kind_name = 'QTY_EXNER_FUNCTION'
        case ('THBM')
          level_type      = 'M' 
          var_kind_name = 'QTY_POTENTIAL_TEMPERATURE'
        case ('EXBW')
          level_type      = 'W' 
          var_kind_name = 'QTY_EXNER_FUNCTION'
      end select

      var = new_state_variable(nest_number,                               &
                               var_kind_num_from_name(var_kind_name),     &
                               is_mean_field,                             &
                               level_type,                                &
                               cur_level_local,                           &
                               pert_magnitude,                            &
                               PERT_TYPE_NOPERTS,                         &
                               FLAG_FREEZE_FIELD,                         &
                               FLAG_ALLOW_NEGATIVE,                       &
                               var_name = variable_name)

    end function define_mean_var

    ! get_var_substate_r8
    ! ---------------------
    ! Return a pointer to the subsection of a kind(r8) vector corresponding 
    ! to a given state variable.
    !  PARAMETERS
    !   IN  vector          The vector to take the section from
    !   IN  var             The state variable to base the section on
    function get_var_substate_r8(var, vector) result(var_substate)
        type(state_variable),                 intent(in)  :: var
        real(kind=r8), dimension(:), target,  intent(in)  :: vector
        real(kind=r8), dimension(:), pointer              :: var_substate

        var_substate => vector(var%state_begin:var%state_end)
    end function get_var_substate_r8

    ! get_var_substate_r4
    ! -----------------------
    ! Return a pointer to the subsection of a kind(r4) vector 
    ! corresponding to a given state variable
    !  PARAMETERS
    !   IN  vector          The vector to take the section from
    !   IN  var             The state variable to base the section on
    function get_var_substate_r4(var, vector) result(var_substate)
        type(state_variable),                     intent(in)  :: var
        real(kind=r4), dimension(:), target,  intent(in)      :: vector
        real(kind=r4), dimension(:), pointer                  :: var_substate

        var_substate => vector(var%state_begin:var%state_end)
    end function get_var_substate_r4

    ! gets_update
    ! -----------
    ! Given a state variable, returns true if the field should be
    ! updated in the COAMPS restart file
    !  PARAMETERS
    !   IN  var             state variable to query
    function gets_update(var)
        type(state_variable), intent(in)  :: var
        logical                           :: gets_update

        gets_update = var%update_field
    end function gets_update

    ! get_nest_number
    ! ---------------
    ! Get the nest number for a particular state variable
    !  PARAMETERS
    !   IN  var             state variable to query
    function get_nest_number(var)
        type(state_variable), intent(in)  :: var
        integer                           :: get_nest_number

        get_nest_number = var%nest_number
    end function get_nest_number

    ! get_state_begin
    ! --------------
    ! Returns the begining index of the state variable
    !  PARAMETERS
    !   IN  var             State variable to query
    function get_state_begin(var)
        type(state_variable), intent(in)  :: var
        integer                           :: get_state_begin

        get_state_begin = var%state_begin
    end function get_state_begin

    ! get_state_end
    ! --------------
    ! Returns the end index of the state variable
    !  PARAMETERS
    !   IN  var             State variable to query
    function get_state_end(var)
        type(state_variable), intent(in)  :: var
        integer                           :: get_state_end

        get_state_end = var%state_end
    end function get_state_end

    ! is_var_at_index
    ! --------------
    ! Returns true if the given variable exists at the index
    !  PARAMETERS
    !   IN  var             State variable to query
    !   IN  indx            index to check
    function is_var_at_index(var, index_in)
        type(state_variable), intent(in)  :: var
        integer,              intent(in)  :: index_in
        logical                           :: is_var_at_index

        is_var_at_index = ( var%state_begin <= index_in .AND. var%state_end >= index_in)
    end function is_var_at_index

    ! is_nonnegative
    ! --------------
    ! Returns true if the given variable is restricted to positive values
    !  PARAMETERS
    !   IN  var             State variable to query
    function is_nonnegative(var)
        type(state_variable), intent(in)  :: var
        logical                           :: is_nonnegative

        is_nonnegative = var%positive_definite
    end function is_nonnegative

    ! get_vert_loc
    ! ------------
    ! Returns the vertical level this variable is defined on
    !  PARAMETERS
    !   IN  var             State variable to query
    function get_vert_loc(var, domain, kk_in)
        type(state_variable), intent(in)  :: var
        type(coamps_domain),  intent(in)  :: domain
        integer, optional,    intent(in)  :: kk_in
        real(kind=r8)                     :: get_vert_loc

        integer :: kk

        if(present(kk_in)) then
           kk = kk_in
        else
           kk = var%sigma_record
        end if

        ! No surface types for now
        if (var%mass_level) then
            get_vert_loc = get_domain_msigma(domain, kk) 
        else
            get_vert_loc = get_domain_wsigma(domain, kk) 
        end if
    end function get_vert_loc

    ! get_vert_type
    ! -------------
    ! Return the type of vertical level this variable's defined on
    !  PARAMETERS
    !   IN  var             State variable to query
    function get_vert_type(var)
        type(state_variable), intent(in)  :: var
        integer                           :: get_vert_type

        get_vert_type = var%vert_type 
    end function get_vert_type

    ! get_vert_value
    ! -------------
    ! Return the value of vertical level this variable's defined on
    !  PARAMETERS
    !   IN  var             State variable to query
    function get_vert_value(var)
        type(state_variable), intent(in)  :: var
        real(kind=r8)                     :: get_vert_value

        get_vert_value = var%vert_value 
    end function get_vert_value

    ! get_var_min_level
    ! ------------
    ! Returns the minimum vertical level of the variable
    !  PARAMETERS
    !   IN  var             state variable to query
    !   IN  domain          domain of the variable
    function get_var_min_level(var, domain) result(level)
        type(state_variable), intent(in)  :: var
        type(coamps_domain),  intent(in)  :: domain
        real(kind=r8)                     :: level 

        if (var%mass_level) then
          level = get_domain_msigma(domain, get_domain_num_levels(domain))
        else
          level = get_domain_wsigma(domain, get_domain_num_levels(domain) + 1)
        end if
    end function get_var_min_level

    ! get_var_max_level
    ! ------------
    ! Returns the maximum vertical level of the variable
    !  PARAMETERS
    !   IN  var             state variable to query
    !   IN  domain          domain of the variable
    function get_var_max_level(var, domain) result(level)
        type(state_variable), intent(in)  :: var
        type(coamps_domain),  intent(in)  :: domain
        real(kind=r8)                     :: level 

        if (var%mass_level) then
            level = get_domain_msigma(domain, 1)
        else
            level = get_domain_wsigma(domain, 1)
        end if
    end function get_var_max_level

    ! get_var_name
    ! ------------
    ! Returns the kind of this restart variable
    !  PARAMETERS
    !   IN  var             Restart variable to query
    function get_var_name(var) result(var_name)
        type(state_variable), intent(in)  :: var
        character(len=VAR_NAME_LEN)       :: var_name

        var_name = var%var_name 
    end function get_var_name

    ! get_mass_level_flag
    ! ------------
    ! Returns the mass level flag
    !  PARAMETERS
    !   IN  var             Restart variable to query
    function get_mass_level_flag(var) result(mass_level)
        type(state_variable), intent(in)  :: var
        logical                           :: mass_level

        mass_level = var%mass_level 
    end function get_mass_level_flag

    ! get_var_kind
    ! ------------
    ! Returns the kind of this restart variable
    !  PARAMETERS
    !   IN  var             Restart variable to query
    function get_var_kind(var)
        type(state_variable), intent(in)  :: var
        integer                           :: get_var_kind

        get_var_kind = var%var_kind 
    end function get_var_kind

    ! get_pert_magnitude
    ! ----------------------
    ! Return the perturbation magnitude for this variable
    !  PARAMETERS
    !   IN  var             State variable to query
    function get_pert_magnitude(var)
        type(state_variable), intent(in)  :: var
        real(kind=r8)                     :: get_pert_magnitude

        get_pert_magnitude = var%pert_mag
    end function get_pert_magnitude

    ! get_pert_type
    ! -------------
    ! Return the perturbation type for this variable
    !  PARAMETERS
    !   IN  var             State variable to query
    function get_pert_type(var)
        type(state_variable), intent(in)  :: var
        integer                           :: get_pert_type

        get_pert_type = var%pert_type
    end function get_pert_type

    ! get_sigma_record
    ! -------------
    ! Return the sigma_record
    !  PARAMETERS
    !   IN  var             State variable to query
    function get_sigma_record(var)
        type(state_variable), intent(in)  :: var
        integer                           :: get_sigma_record

        get_sigma_record = var%sigma_record
    end function get_sigma_record

    ! set_sigma_record
    ! -------------
    ! Return the sigma_record
    !  PARAMETERS
    !   IN  var             State variable to query
    subroutine set_sigma_record(var,sigma_record)
        type(state_variable), intent(inout)  :: var
        integer             , intent(in)     :: sigma_record

        var%sigma_record = sigma_record
    end subroutine set_sigma_record

    ! get_record_in_nest
    ! -------------------
    ! Queries the given variable and returns the record number of a
    ! particular sigma level within its particular nest in the large 
    ! COAMPS restart file
    !  PARAMETERS
    !   IN  var                 state variable to pull compute position of
    !   IN  domain              COAMPS domain that the restart var is on
    !   IN  used_single_io      true if restart file is written using
    !                           a single I/O processor
    !   IN  used_digital_filter true if COAMPS used a digital filter (this
    !                           is "ldigit" in the COAMPS namelist)
    function get_record_in_nest(var, domain, used_single_io, &
                                used_digital_filter)
        type(state_variable), intent(in)  :: var
        type(coamps_domain),  intent(in)  :: domain
        logical,              intent(in)  :: used_single_io
        logical, optional,    intent(in)  :: used_digital_filter
        integer                           :: get_record_in_nest

        ! How many total variables of each kind there are - use this for
        ! calculating offsets between dimensions
        integer :: N2D, N3D

        ! Temporary shorthand so we don't need to keep referencing the
        ! structure in the select statement
        integer :: kka, k, n

        ! How many of each kind of variable we have depends on if we are
        ! doing single or multiple process I/O
        if (used_single_io) then
            N2D = 101
            N3D = 69
        else
            N2D = 112

            ! This is the sum of 28 "normal" 3d variables and 28 variables
            ! that gets used only if the digital filter doesn't
            N3D = 56
            if (present(used_digital_filter)) then
                if (used_digital_filter) N3D = 28
            end if
        end if

        kka = get_domain_num_levels(domain)
        k   = var%sigma_record

        if (used_single_io) then
            n = var%var_record(SINGLEIO)
        else
            n = var%var_record(MULTIIO)
        end if

        ! Take into account the number of variables that occur before the 
        ! target variable and all their associated sigma levels, then all 
        ! the sigma levels of THIS variable that have occured before the 
        ! target level.
        select case (var%dim_type)
        case (DIM_TYPE_2D)
            get_record_in_nest = n
        case (DIM_TYPE_3D)
            get_record_in_nest = (N2D) + ((kka) * (n-1) + 1) + (k - 1)
        case (DIM_TYPE_3DW)
            get_record_in_nest = (N2D + (N3D * kka)) + &
                                 ((kka + 1) * (n - 1) + 1) + (k - 1)
        case default
            call error_handler(E_ERR,'get_record_in_nest',             &
                               'Unrecognized dimension type!', source, &
                               revision, revdate)
        end select
    end function get_record_in_nest

    ! set_position
    ! ------------
    ! Figure the points in the state vector that correspond to this state
    ! variable based on its corresponding nest size.  Store the limits in 
    ! the state variable structure and return the location of this 
    ! variable's last point
    !  PARAMETERS
    ! INOUT var                 State variable to calculate the region for
    !   IN  initial_position    Offset of this variable in the state vector
    !   IN  domain              The COAMPS domain connected to this state
    function set_position(var, initial_position, domain)
        type(state_variable), intent(inout) :: var
        integer,              intent(inout) :: initial_position
        type(coamps_domain),  intent(in)    :: domain
        integer                             :: set_position

        var%state_begin = initial_position + 1
        var%state_end   = (var%state_begin - 1) +                          &
                          get_nest_size(get_domain_nest(domain, var%nest_number)) &
                        * get_var_nlevs(var, domain)

        !print '(A,A,I2.2,1x,A,I3.3,A,1x,A,1x,I10,1x,A,1x,I10)',  &
        !       trim(var%var_name), '_g', var%nest_number,   &
        !      '(sigma_index = ', var%sigma_record ,')',     &
        !      'stretches from', var%state_begin, 'to', var%state_end

        set_position = var%state_end
    end function set_position

    ! set_nc_varid
    !  PARAMETERS
    !  INOUT var        State variable to modify
    !  IN    varid      netcdf variable id
    subroutine set_nc_varid(var, nc_varid)
        type(state_variable), intent(inout) :: var
        integer,              intent(in)    :: nc_varid

        var%nc_varid = nc_varid
    end subroutine set_nc_varid

    ! get_nc_varid
    !  PARAMETERS
    !  INOUT var        State variable to modify
    !  IN    varid      netcdf variable id
    function get_nc_varid(var)
        type(state_variable), intent(in) :: var
        integer                          :: get_nc_varid

        get_nc_varid = var%nc_varid
    end function get_nc_varid

    ! get_io_flag
    ! -------------------
    ! Gets flag indicating if var is to be io'd
    !  PARAMETERS
    !   IN     var               state variable
    function get_io_flag(var)
        type(state_variable), intent(in) :: var
        logical                          :: get_io_flag

        get_io_flag=var%io_flag
    end function get_io_flag

    ! get_mean_flag
    ! -------------------
    ! Gets flag indicating if var is mean field
    !  PARAMETERS
    !   IN     var               state variable
    function get_mean_flag(var)
        type(state_variable), intent(in) :: var
        logical                          :: get_mean_flag

        get_mean_flag=var%mean_field
    end function get_mean_flag

    ! set_var_stagger
    ! -------------------
    ! Sets the staggering properties of a variable
    subroutine set_var_stagger(var)
      type(state_variable), intent(inout) :: var
    
            if(var%is_2d_variable) then
              var%var_stagger = 'T'
            else
              select case(trim(uppercase(var%var_name)))
                case('UUWIND')
                  if(is_sigma_level(var)) then
                  var%var_stagger = 'U'
                  else
                    var%var_stagger = 'T'
                  end if
                case('VVWIND')
                  if(is_sigma_level(var)) then
                  var%var_stagger = 'V'
                  else
                    var%var_stagger = 'T'
                  end if
                case('WWWIND', 'EXBW')
                  var%var_stagger = 'W'
                case default
                  var%var_stagger = 'T'
              end select
            end if
    end subroutine set_var_stagger

    ! is_sigma_level
    ! -------------
    ! Return the true if the var is defined on a sigma level (W or M)
    !  PARAMETERS
    !   IN  var             State variable to query
    function is_sigma_level(var)
        type(state_variable), intent(in)  :: var
        logical                           :: is_sigma_level

        if(var%vert_type .eq. VERTISLEVEL) then
           is_sigma_level = .true.
        else
           is_sigma_level = .false.
        end if
    end function is_sigma_level

    ! get_var_stagger
    ! -------------------
    ! Given an iterator to a variable list returns the stagering of that variable
    !  PARAMETERS
    !   OUT variable_stagger   Staggering of field (U|V|W|T)
    !   IN  var                state_variable
    function get_var_stagger(var) result(var_stagger)
      type(state_variable), intent(in) :: var
      character(len=1)                 :: var_stagger
      var_stagger = var%var_stagger
    end function get_var_stagger

    ! set_twod_flag
    ! -------------------
    ! Sets a flag to indicate if this is a two-dim field
    !  PARAMETERS
    !   INOUT  var               state variable
    !   IN     twod_flag         flag to set to
    subroutine set_2d_flag(var, twod_flag)
        type(state_variable), intent(inout) :: var
        logical,              intent(in)    :: twod_flag
        var%is_2d_variable = twod_flag
    end subroutine set_2d_flag

    ! dump_state_variable
    ! -------------------
    ! Writes out the contents of a state variable
    !  PARAMETERS
    !   IN  var               state variable to dump 
    subroutine dump_state_variable(var)
        type(state_variable), intent(in) :: var

        write (*, '(A12,T15,A10)')   "VAR NAME:", trim(var%var_name)
        write (*, '(A12,T15,I9.9,1x,I9.9)')   "STATE BNDRY:", var%state_begin, var%state_end
        !write (*, '(A12,T15I8.8)') "DIM TYPE:", var%dim_type
        !write (*, '(A12,T15I8.8)') "VAR RCRD:", var%var_record
        write (*, '(A12,T15,I8.8)') "SGM RCRD:", var%sigma_record
        write (*, '(A12,T15,F8.6)') "PTRB PCT:", var%pert_mag
        write (*, '(A12,T15,I8.8)') "PTB TYPE:", var%pert_type
        write (*, '(A12,T15,L1)')   "UPDATE??:", var%update_field
        write (*, '(A12,T15,A32)')  "VAR TYPE:", get_name_for_quantity(var%var_kind)
        write (*, '(A12,T15,I8.8)') "VAR NMBR:", var%var_kind
        write (*, '(A12,T15,L1)')   "MASSLVL?:", var%mass_level
        write (*, '(A12,T15,L1)')   "MEANFLD?:", var%mean_field
        write (*, '(A12,T15,L1)')   "IS_2D_VAR?:", var%is_2d_variable
        write (*, '(A12,T15,L1)')   "IO_FLAG?:", var%io_flag
        !write (*, '(A12,T15,A1)')   "STAGGERING:", var%var_stagger
        write (*, *) "-----------------------------------------"

    end subroutine dump_state_variable
    
    !------------------------------
    ! END PUBLIC ROUTINES
    !------------------------------

    !------------------------------
    ! BEGIN PRIVATE ROUTINES
    !------------------------------

    ! get_var_nlevs
    ! -------------------
    ! Returns the number of vertical levels in a state_variable
    !  PARAMETERS
    !   IN  var               state_variable
    function get_var_nlevs(var, domain) result(nlevels)
        type(state_variable), intent(in) :: var
        type(coamps_domain),  intent(in) :: domain
        integer                          :: nlevels

        if(var%is_2d_variable) then 
          nlevels = 1
        else
          nlevels = get_domain_num_levels(domain)
          if(.not. var%mass_level) then
            nlevels = nlevels + 1
          end if
        end if
    end function get_var_nlevs

    ! var_kind_num_from_name
    ! ----------------------
    ! Converts a character string describing a the kind of variable to a
    ! numeric representation of that variable kind.
    function var_kind_num_from_name(var_kind_name)
        character(len=*), intent(in)  :: var_kind_name
        integer                       :: var_kind_num_from_name

        var_kind_num_from_name = get_index_for_quantity(var_kind_name)

        if (var_kind_num_from_name .lt. 0) then
            call error_handler(E_ERR, 'var_kind_num_from_name',           &
                               'Error in locating numeric equivalent ' // &
                               'to raw observation kind' // var_kind_name,& 
                               source, revision, revdate)
        end if
    end function var_kind_num_from_name

    ! pert_type_from_name
    ! -------------------
    ! Converts a character string describing a perturbation type to the
    ! internal numeric representation
    function pert_type_from_name(pert_name)
        character(len=*), intent(in)  :: pert_name
        integer                       :: pert_type_from_name

        select case (adjustl(pert_name))
        case (adjustl(NOPERTS_NAME))
            pert_type_from_name = PERT_TYPE_NOPERTS
        case (adjustl(UNIFORM_NAME))
            pert_type_from_name = PERT_TYPE_UNIFORM
        case (adjustl(INDIVID_NAME))
            pert_type_from_name = PERT_TYPE_INDIVID
        case default
            call error_handler(E_ERR, 'pert_type_from_name',           &
                               'Perturbation type ' // pert_name //    &
                               ' not found!', source, revision, revdate)
        end select
    end function pert_type_from_name

    ! update_flag_from_type
    ! ---------------------
    ! Convert character representation of whether the field should be
    ! written back into the COAMPS restart file (UPDATE) or skipped
    ! (FREEZE) to the internal logical representation
    function update_flag_from_type(update_string)
        character(len=*), intent(in)  :: update_string
        logical                       :: update_flag_from_type

        select case (adjustl(update_string))
        case (adjustl(FIELD_UPDATE))
            update_flag_from_type = FLAG_UPDATE_FIELD
        case (adjustl(FIELD_FREEZE))
            update_flag_from_type = FLAG_FREEZE_FIELD
        case default
            call error_handler(E_ERR, 'update_flag_from_type',         &
                               'Update type' // update_string //       &
                               ' not found!', source, revision, revdate)
        end select
    end function update_flag_from_type

    ! mass_flag_from_level_type
    ! -------------------------
    ! Convert character representation of [M]ass or [W] velocity level
    ! type to an internal logical type 
    function mass_flag_from_level_type(level_type_char)
        character, intent(in)  :: level_type_char
        logical                :: mass_flag_from_level_type

        select case (adjustl(level_type_char))
        case (adjustl(MASS_LEVEL))
            mass_flag_from_level_type = IS_MASS_LEVEL
        case (adjustl(W_LEVEL))
            mass_flag_from_level_type = IS_W_LEVEL
        case default
            mass_flag_from_level_type = IS_NOT_SIG
        end select
    end function mass_flag_from_level_type

    ! vert_type_from_level_type
    ! -------------------------
    ! Convert character representation of [M]ass or [W] velocity level
    ! type to an internal logical type 
    function vert_type_from_level_type(level_type_char)
        character, intent(in)  :: level_type_char
        integer                :: vert_type_from_level_type

        select case (adjustl(level_type_char))
        case (adjustl(MASS_LEVEL), adjustl(W_LEVEL))
          vert_type_from_level_type = VERTISLEVEL
        case (adjustl(ZHT_LEVEL))
          vert_type_from_level_type = VERTISHEIGHT
        case (adjustl(SFC_LEVEL))
          vert_type_from_level_type = VERTISSURFACE
        case (adjustl(PRES_LEVEL))
          vert_type_from_level_type = VERTISPRESSURE
        case (adjustl(UNDEF_LEVEL))
          vert_type_from_level_type = VERTISUNDEF
        case default
            call error_handler(E_ERR, 'vert_type_from_level_type', 'Level '//&
                               'type' // level_type_char // ' not found!',   &
                               source, revision, revdate)
        end select
    end function vert_type_from_level_type

    ! posdef_flag_from_name
    ! ---------------------
    ! Converts a character string describing whether the entire field
    ! is positive definite (ISPOSDEF) or not (NOPOSDEF) to the internal
    ! logical representation 
    function posdef_flag_from_name(posdef_string)
        character(len=*), intent(in)  :: posdef_string
        logical                       :: posdef_flag_from_name

        select case (adjustl(posdef_string))
        case (adjustl(ALLOW_NEGATIVE))
            posdef_flag_from_name = FLAG_ALLOW_NEGATIVE
        case (adjustl(FORCE_POSITIVE))
            posdef_flag_from_name = FLAG_FORCE_POSITIVE
        case default
            call error_handler(E_ERR, 'posdef_flag_from_name',          &
                               'Positive definite ' // posdef_string // &
                               ' not found!', source, revision, revdate)
        end select
    end function posdef_flag_from_name

    ! is_same_variable
    ! ----------------
    ! Return true if the variables are the "same" - i.e. they are describing
    ! identical entries.  
    function is_same_variable(var1, var2) result(is_same)
        type(state_variable), intent(in)  :: var1
        type(state_variable), intent(in)  :: var2
        logical                           :: is_same

        is_same = is_same_kind(var1, var2)       .and. &
                  is_same_mean(var1, var2)       .and. &
                  is_same_mass_level(var1, var2) .and. &
                  is_same_nest(var1, var2)       .and. &
                  is_same_sigma_index(var1, var2) .and. &
                  is_same_vert_type_and_value(var1, var2)
    end function is_same_variable
    
    ! is_same_vert_type_and_value
    ! ----------------
    ! Return true if the variables are the "same" - i.e. they are describing
    ! entries with the same vert_type (SIGMA, HEIGHT, PRESSURE, SURFACE, UNDEF)
    ! If vert_type is not SIGMA then variables need to have the same vert value.
    function is_same_vert_type_and_value(var1, var2) result(is_same)
        type(state_variable), intent(in)  :: var1
        type(state_variable), intent(in)  :: var2
        logical                           :: is_same
        real(kind=r8), parameter :: ERR_TOL = 0.001_r8
    
        is_same = var1%vert_type .eq. var2%vert_type
        if( (var1%vert_type .eq. VERTISLEVEL) .or. &
            (var1%vert_type .eq. VERTISUNDEF) .or. &
            .not. is_same) return

        is_same = (abs(var1%vert_value - var2%vert_value) <= ERR_TOL)  &
                  .and. is_same

    end function is_same_vert_type_and_value

    ! is_same_kind_and_nest
    ! ----------------
    ! Return true if the variables are the "same" - i.e. they are describing
    ! entries with the same kind and nest.  
    function is_same_kind_and_nest(var1, var2) result(is_same)
        type(state_variable), intent(in)  :: var1
        type(state_variable), intent(in)  :: var2
        logical                           :: is_same
    
        is_same = is_same_kind(var1, var2) .and.  is_same_nest(var1, var2)
    end function is_same_kind_and_nest
    
    ! is_same_kind
    ! -------------
    ! Returns true if the variables have the same kind
    function is_same_kind(var1, var2)
        type(state_variable), intent(in)  :: var1
        type(state_variable), intent(in)  :: var2
        logical                           :: is_same_kind

        is_same_kind = var1%var_kind .eq. var2%var_kind
    end function is_same_kind

    ! is_same_mean
    ! ------------
    ! Returns true if the variables have the same mean variable status
    function is_same_mean(var1, var2)
        type(state_variable), intent(in)  :: var1
        type(state_variable), intent(in)  :: var2
        logical                           :: is_same_mean

        is_same_mean = var1%mean_field .eqv. var2%mean_field
    end function is_same_mean

    ! is_same_mass_level
    ! ------------------
    ! Returns true if the variables are defined on the same type of vertical
    ! sigma level.
    function is_same_mass_level(var1, var2)
        type(state_variable), intent(in)  :: var1
        type(state_variable), intent(in)  :: var2
        logical                           :: is_same_mass_level

        is_same_mass_level = var1%mass_level .eqv. var2%mass_level
    end function is_same_mass_level

    ! is_same_sigma_index
    ! -------------------
    ! Returns true if the variables are defined on the same sigma level
    function is_same_sigma_index(var1, var2)
        type(state_variable), intent(in)  :: var1
        type(state_variable), intent(in)  :: var2
        logical                           :: is_same_sigma_index
    
        is_same_sigma_index = var1%sigma_record .eq. var2%sigma_record
    end function is_same_sigma_index 

    ! is_same_nest
    ! -------------------
    ! Returns true if the variables are defined on the same nest
    function is_same_nest(var1, var2)
        type(state_variable), intent(in)  :: var1
        type(state_variable), intent(in)  :: var2
        logical                           :: is_same_nest
    
        is_same_nest = var1%nest_number .eq. var2%nest_number
    end function is_same_nest

    !------------------------------
    ! END PRIVATE ROUTINES
    !------------------------------

end module coamps_statevar_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
