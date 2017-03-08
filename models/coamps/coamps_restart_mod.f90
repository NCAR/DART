! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

module coamps_restart_mod

!------------------------------
! MODULE:       coamps_restart_mod
! AUTHOR:       T. R. Whitcomb
!               Naval Research Laboratory
! DART VERSION: Jamaica
!
! Module containing the data structure and routines for dealing with
! COAMPS restart files and the dynamic DART state vector definition
!------------------------------ 

  use coamps_grid_mod, only : coamps_grid,         &
                              get_grid_dims,       &
                              get_grid_field_size, &
                              initialize_grid
  use coamps_util_mod, only : check_alloc_status,  &
                              check_io_status 

  use location_mod,    only : VERTISLEVEL,   &
                              VERTISSURFACE
  use obs_kind_mod
  use types_mod,       only : r8
  use utilities_mod,   only : do_output,     &
                              E_ERR,         &
                              E_MSG,         &
                              E_WARN,        &
                              error_handler, &
                              get_unit

  implicit none

  private

  !------------------------------
  ! BEGIN PUBLIC INTERFACE
  !------------------------------

  ! Initialization
  public :: initialize_restart_info

  ! Size information
  public :: get_num_vars

  ! Grid information
  public :: get_restart_grid

  ! Search functions
  public :: get_restart_index_by_properties

  ! Accessors by index
  public :: get_full_record_num_by_index
  public :: get_sigma_by_index
  public :: get_var_type_by_index
  public :: get_mean_flag_by_index
  public :: get_posdef_flag_by_index
  public :: get_vert_coord_by_index
  public :: get_pert_magnitude_by_index
  public :: get_pert_type_by_index
  public :: get_kind_by_index
  public :: get_update_flag_by_index
  public :: get_var_info_by_abs_index

  ! Numeric constants
  public :: PERT_TYPE_NOPERTS
  public :: PERT_TYPE_UNIFORM
  public :: PERT_TYPE_INDIVID

  ! Diagnostics
  public :: dump_restart_vars
  
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
    subroutine get_name_info(DIM_TYPE_2D, DIM_TYPE_3D, DIM_TYPE_3DW,&
                             SINGLEIO, MULTIIO, var_name,           &
                             var_dim_type, var_record_num)
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

  !------------------------------
  ! END EXTERNAL INTERFACE
  !------------------------------

  !------------------------------
  ! BEGIN TYPES AND CONSTANTS
  !------------------------------

  ! Dimension information - pretty straightforward except that 3Dw
  ! includes the surface, so there are kka+1 levels there
  integer, parameter :: DIM_TYPE_2D  = 1
  integer, parameter :: DIM_TYPE_3D  = 2
  integer, parameter :: DIM_TYPE_3DW = 3

  ! Length of some string variables
  integer, parameter :: VAR_NAME_LEN  = 5
  integer, parameter :: VAR_TYPE_LEN  = 32
  integer, parameter :: PERT_TYPE_LEN = 7

  ! Perturbation types - don't perturb at all, perturb whole sigma
  ! level field uniformly, or perturb each grid point individually.
  integer, parameter :: PERT_TYPE_NOPERTS = 0
  integer, parameter :: PERT_TYPE_UNIFORM = 1
  integer, parameter :: PERT_TYPE_INDIVID = 2
  character(len=PERT_TYPE_LEN), parameter :: NOPERTS_NAME = 'NOPERTS' 
  character(len=PERT_TYPE_LEN), parameter :: UNIFORM_NAME = 'UNIFORM'
  character(len=PERT_TYPE_LEN), parameter :: INDIVID_NAME = 'INDIVID'

  ! If the variable is defined on a mass level or a w level
  character, parameter :: MASS_LEVEL = 'M'
  character, parameter :: W_LEVEL    = 'W'

  ! Indices for the variable record array to take into account the
  ! difference between the numbering in the single I/O case and
  ! multiple process I/O case
  integer, parameter :: SINGLEIO = 1
  integer, parameter :: MULTIIO  = 2

  ! Whether or not to write the field back to the COAMPS restart file
  logical, parameter :: FLAG_UPDATE_FIELD  = .true.
  logical, parameter :: FLAG_FREEZE_FIELD  = .false.

  ! Define the strings used to signify whether we'll do the updates
  ! or not
  integer, parameter :: FIELD_UPDATE_LEN = 6
  character(len=FIELD_UPDATE_LEN), parameter :: FIELD_UPDATE='UPDATE'
  character(len=FIELD_UPDATE_LEN), parameter :: FIELD_FREEZE='FREEZE'

  ! Define the strings used to signify whether we'll check for
  ! positive definiteness or not and their corresponding flags
  integer, parameter :: POS_DEF_LEN = 8
  character(len=POS_DEF_LEN), parameter :: FORCE_POSITIVE='ISPOSDEF'
  character(len=POS_DEF_LEN), parameter :: ALLOW_NEGATIVE='NOPOSDEF'
  logical, parameter :: FLAG_FORCE_POSITIVE = .true.
  logical, parameter :: FLAG_ALLOW_NEGATIVE = .false.

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
  type :: restart_var
     character(len=VAR_NAME_LEN) :: var_name
     integer                     :: dim_type
     integer,dimension(2)        :: var_record
     integer                     :: sigma_record
     real(kind=r8)               :: pert_mag
     integer                     :: pert_type
     logical                     :: update_field
     integer                     :: var_type
     logical                     :: mean_field
     logical                     :: mass_level
     logical                     :: positive_definite
  end type restart_var

  ! Derived data type that contains the information used to read and
  ! store the field information for the state vector
  type :: coamps_restart
     character(len=10)                            :: cdtg
     type(coamps_grid)                            :: grid
     type(restart_var), dimension(:), allocatable :: restart_vars
  end type coamps_restart

  !------------------------------
  ! END TYPES AND CONSTANTS
  !------------------------------

  !------------------------------
  ! BEGIN MODULE VARIABLES
  !------------------------------

  type(coamps_restart), target :: restart_info

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

  ! initialize_restart_info
  ! -----------------------
  ! Populate the restart_info structure by reading in the grid data
  ! and reading the state vector definition file
  !  PARAMETERS
  !   IN  dtg                  base date-time group for this run
  !   IN  restart_var_filename name of state vector definition file
  subroutine initialize_restart_info(dtg,restart_var_filename)
    character(len=10), intent(in) :: dtg
    character(len=*), intent(in)  :: restart_var_filename

    integer           :: restart_var_unit

    character(len=*), parameter :: routine = 'initialize_restart_info'
    integer           :: io_status

    restart_info%cdtg = dtg

    call initialize_grid(restart_info%cdtg, restart_info%grid)

    restart_var_unit = get_unit()
    open(unit=restart_var_unit, file=restart_var_filename, &
         status='old', access='sequential', action='read', &
         form='formatted', iostat=io_status)
    call check_io_status(io_status, routine, source, revision, revdate, &
                         'Opening file ' // restart_var_filename)

    call read_restart_var_file(restart_var_unit)

    close(restart_var_unit)

  end subroutine initialize_restart_info

  ! get_num_vars
  ! ------------
  ! Returns how many fields there are in the state vector
  !  PARAMETERS
  !   OUT num_vars          number of fields in state vector
  subroutine get_num_vars(num_vars)
    integer, intent(out) :: num_vars
    num_vars = size(restart_info%restart_vars)
  end subroutine get_num_vars

  ! get_restart_grid
  ! ----------------
  ! Returns the grid structure associated with the state vector
  ! definition
  !  PARAMETERS
  !   OUT restart_grid      coamps_grid this state vector holds
  subroutine get_restart_grid(restart_grid)
    type(coamps_grid), intent(out) :: restart_grid
    restart_grid = restart_info%grid
  end subroutine get_restart_grid

  ! get_restart_index_by_properties
  ! -------------------------------
  ! Given the type of a variable, whether it's a mean field, and the
  ! type of level it is defined on plus the corresponding sigma index
  ! returns the position of that variable in the restart_vars array.
  !  PARAMETERS
  !   IN  var_type          integer form of variable type
  !   IN  is_mean           true if variable is a mean field
  !   IN  on_mass_level     true if variable is on a mass level
  !   IN  sig_index         sigma level index for variable
  !   OUT restart_index     position in the restart_vars array
  subroutine get_restart_index_by_properties(var_type, is_mean, &
                                             on_mass_level,     &
                                             sig_index,         &
                                             restart_index)
    integer, intent(in)  :: var_type
    logical, intent(in)  :: is_mean
    logical, intent(in)  :: on_mass_level
    integer, intent(in)  :: sig_index
    integer, intent(out) :: restart_index

    integer :: cur_index, total_vars
    integer :: cur_var_type
    integer :: cur_sigma
    logical :: cur_mean_flag
    logical :: cur_mass_flag

    character(len=128) :: message

    restart_index = -1

    ! Assume there is only one - if there are multiples, return the
    ! first one found
    call get_num_vars(total_vars)
    do cur_index = 1, total_vars
       cur_var_type  = restart_info%restart_vars(cur_index)%var_type
       cur_sigma     = restart_info%restart_vars(cur_index)%sigma_record
       cur_mean_flag = restart_info%restart_vars(cur_index)%mean_field
       cur_mass_flag = restart_info%restart_vars(cur_index)%mass_level

       if (cur_var_type .eq. var_type) then
          if (cur_sigma .eq. sig_index) then
             if (cur_mean_flag .eqv. is_mean) then
                if (cur_mass_flag .eqv. on_mass_level) then
                    restart_index = cur_index
                    exit
                end if
             end if
          end if
       end if
    end do

    ! The proper error handling in this case is simply to allow the
    ! negative restart index to stick around - let the caller handle
    ! dealing with this.  Since some variables might not be available
    ! at all levels, quitting would be a bad idea.
    if (restart_index .le. 0) then 
      write (message,*) "Could not find variable type ", var_type, &
                        "on (mass?",on_mass_level,") level",       &
                        sig_index
       call error_handler(E_WARN, "get_restart_index_by_properties",&
                         message, source, revision, revdate) 
    end if
  end subroutine get_restart_index_by_properties

  ! get_full_record_num_by_index
  ! ----------------------------
  ! Call get_full_record_num using a variable at the specified index
  ! in the state vector definition
  !  PARAMETERS
  !   IN  var_index         index of variable to interrogate
  !   IN  use_singleio          true if restart file is written using
  !                         a single I/O processor
  !   OUT recnum            variable's record number within the
  !                         restart file
  subroutine get_full_record_num_by_index(var_index, use_singleio, &
                                          recnum)
    integer, intent(in)  :: var_index
    logical, intent(in)  :: use_singleio
    integer, intent(out) :: recnum

    if (var_index > size(restart_info%restart_vars)) then
       call error_handler(E_ERR, 'get_full_record_num_by_index', &
                          'Variable index out of range', source, &
                          revision, revdate)
    else
       call get_full_record_num(restart_info%restart_vars(var_index),&
                                use_singleio, recnum)
    end if
  end subroutine get_full_record_num_by_index

  ! get_sigma_by_index
  ! ------------------
  ! Returns the sigma level of a field at the specified index in the
  ! state vector.
  !  PARAMETERS
  !   IN  var_index         variable index in state vector to query
  !   OUT sigma_index       variable's sigma level
  subroutine get_sigma_by_index(var_index, sigma_index)
    integer, intent(in)  :: var_index
    integer, intent(out) :: sigma_index

    sigma_index = restart_info%restart_vars(var_index)%sigma_record
  end subroutine get_sigma_by_index

  ! get_var_type_by_index
  ! ---------------------
  ! Return the variable type for the variable at the given index in
  ! the state vector definition
  !  PARAMETERS
  !   IN  restart_index     index of variable to interrogate
  !   OUT var_type          type of given variable
  subroutine get_var_type_by_index(restart_index, var_type)
    integer, intent(in)  :: restart_index
    integer, intent(out) :: var_type

    var_type = restart_info%restart_vars(restart_index)%var_type
  end subroutine get_var_type_by_index

  ! get_mean_flag_by_index
  ! ----------------------
  ! Returns if the field at the specified index in the state vector
  ! is a mean field
  !  PARAMETERS
  !   IN  var_index         variable index in state vector to query
  !   OUT mean_field        true if the field is a mean field
  subroutine get_mean_flag_by_index(var_index, mean_field)
    integer, intent(in)  :: var_index
    logical, intent(out) :: mean_field

    mean_field = restart_info%restart_vars(var_index)%mean_field
  end subroutine get_mean_flag_by_index

  ! get_posdef_flag_by_index
  ! ----------------------------
  ! Call get_posdef_flag usign a variable at the specified index in
  ! the state vector definition
  !  PARAMETERS
  !   IN  var_index         index of variable to interrogate
  !   OUT update            true if this variable gets updated
  subroutine get_posdef_flag_by_index(var_index, posdef)
    integer, intent(in)  :: var_index
    logical, intent(out) :: posdef

    if (var_index > size(restart_info%restart_vars)) then
       call error_handler(E_ERR, 'get_posdef_flag_by_index',      &
                          'Variable index out of range', source,  &
                          revision, revdate)
    else
       call get_posdef_flag(restart_info%restart_vars(var_index), &
                            posdef)
    end if
  end subroutine get_posdef_flag_by_index

  ! get_vert_coord_by_index
  ! -----------------------
  ! Call get_vert_coord using a variable at the specified index in 
  ! the state vector definition
  !  PARAMETERS
  !   IN  var_index         index of variable to interrogate
  !   OUT vert_type         variable's vertical coordinate type
  !   OUT vert_level        variable's vertical coordinate
  subroutine get_vert_coord_by_index(var_index,vert_type,vert_level)
    integer, intent(in)        :: var_index
    integer, intent(out)       :: vert_type
    real(kind=r8), intent(out) :: vert_level

    if (var_index > size(restart_info%restart_vars)) then
       call error_handler(E_ERR, 'get_vert_coord_by_index',      &
                          'Variable index out of range', source, &
                          revision, revdate)
    else
       call get_vert_coord(restart_info%restart_vars(var_index), & 
                           restart_info%grid, vert_type, vert_level)
    end if
  end subroutine get_vert_coord_by_index

  ! get_pert_magnitude_by_index
  ! ---------------------------
  ! Returns the perturbation magnitude of a field at the specified
  ! index in the state vector.
  !  PARAMETERS
  !   IN  var_index         variable index in state vector to query
  !   OUT pert_mag          perturbation magnitude
  subroutine get_pert_magnitude_by_index(var_index, pert_mag)
    integer, intent(in)        :: var_index
    real(kind=r8), intent(out) :: pert_mag

    pert_mag = restart_info%restart_vars(var_index)%pert_mag
  end subroutine get_pert_magnitude_by_index

  ! get_pert_type_by_index
  ! ----------------------
  ! Returns the perturbation type (no perturbation, entire field 
  ! perturbed by same value, or each point perturbed with a random
  ! value) of a field at the specified index in the state vector.
  !  PARAMETERS
  !   IN  var_index         variable index in state vector to query
  !   OUT pert_type         perturbation type (defined as integer
  !                         constants in this module)
  subroutine get_pert_type_by_index(var_index, pert_type)
    integer, intent(in)  :: var_index
    integer, intent(out) :: pert_type
    
    pert_type = restart_info%restart_vars(var_index)%pert_type
  end subroutine get_pert_type_by_index

  ! get_kind_by_index
  ! -----------------
  ! Returns the variable kind of a field at the specified index in
  ! the state vector
  !  PARAMETERS
  !   IN  var_index         variable index in state vector to query
  !   OUT var_kind          numeric variable kind for that variable
  subroutine get_kind_by_index(var_index, var_kind)
    integer, intent(in)  :: var_index
    integer, intent(out) :: var_kind
    
    var_kind = restart_info%restart_vars(var_index)%var_type
  end subroutine get_kind_by_index

  ! get_update_flag_by_index
  ! ------------------------
  ! Call get_update_flag using a variable at the specified index in
  ! the state vector definition.
  !  PARAMETERS
  !   IN  var_index         index of variable to interrogate
  !   OUT update            true if this variable gets updated
  subroutine get_update_flag_by_index(var_index, update)
    integer, intent(in)  :: var_index
    logical, intent(out) :: update

    if (var_index > size(restart_info%restart_vars)) then
       call error_handler(E_ERR, 'get_update_flag_by_index',      &
                          'Variable index out of range', source,  &
                          revision, revdate)
    else
       call get_update_flag(restart_info%restart_vars(var_index), &
                            update)
    end if
  end subroutine get_update_flag_by_index

  ! get_var_info_by_abs_index
  ! -------------------------
  ! Maps an entry in the long state vector to a particular field,
  ! then returns that field's name and sigma level index.
  !  PARAMETERS
  !   IN  abs_index         Arbitrary index of an element in the
  !                         state vector
  !   OUT var_name          Name of the field containing the element
  !   OUT sigma_index       Sigma index of the field containing the 
  !                         element
  subroutine get_var_info_by_abs_index(abs_index, var_name, level)
    integer, intent(in)                      :: abs_index
    character(len=VAR_NAME_LEN), intent(out) :: var_name
    integer, intent(out)                     :: level

    integer :: invar
    integer :: gridsize

    call get_grid_field_size(restart_info%grid, gridsize)

    invar = int(abs_index / gridsize) + 1

    var_name = restart_info%restart_vars(invar)%var_name
    level = restart_info%restart_vars(invar)%sigma_record
  end subroutine get_var_info_by_abs_index

  ! dump_restart_vars
  ! -----------------
  ! Write out all the components of the state vector
  !  PARAMETERS
  !   [none]
  subroutine dump_restart_vars()
    integer :: ii

    if (do_output()) then
       write (*,*) "*** Restart File Contents ***"
       write (*,*) "-----------------------------"
    end if
    
    do ii=1,size(restart_info%restart_vars)
       call dump_restart_var(restart_info%restart_vars(ii))
    end do
  end subroutine dump_restart_vars

  !------------------------------
  ! END PUBLIC ROUTINES
  !------------------------------

  !------------------------------
  ! BEGIN PRIVATE ROUTINES
  !------------------------------

  ! read_restart_var_file
  ! ---------------------
  ! Reads in the state vector definition file - the first line 
  ! is the number of variables to follow, then each line after that
  ! contains information on each entry within the state vector
  !  PARAMETERS
  !   IN  statevec_def_unit    Unit number for the state vector
  !                            definition file
  subroutine read_restart_var_file(statevec_def_unit)
    integer, intent(in) :: statevec_def_unit

    integer                         :: num_recs

    character(len=*), parameter     :: routine = 'read_restart_var_file'
    logical                         :: is_opened
    integer                         :: io_status
    integer                         :: alloc_status

    ! Things we read in from the file
    character(len=VAR_NAME_LEN)     :: name      ! Variable name
    integer                         :: index     ! Sigma index
    character(len=PERT_TYPE_LEN)    :: type_name ! text pert type
    real(kind=r8)                   :: pert_mag  ! perturbation 
                                                 ! magnitude
    character                       :: lvl_type  ! mass or w level
    character(len=VAR_TYPE_LEN)     :: var_type  ! variable type 
    character(len=FIELD_UPDATE_LEN) :: write_fld ! write into restart
    logical                         :: mean_flag ! Mean or not
    character(len=POS_DEF_LEN)      :: pos_def   ! Positive definite?
    
    ! Things that we convert what we read in the file to
    integer                         :: pert_type  ! Numeric pert type
    logical                         :: lvl_flag   ! true if mass level
    logical                         :: write_flag ! Write/Skip flag
    logical                         :: pos_flag   ! positive definite?
    integer                         :: var_type_num !numeric var type

    integer :: cur_rec
    
    ! Assert
    inquire(unit=statevec_def_unit, opened=is_opened)
    if (.not. is_opened) then
       call error_handler(E_ERR, 'read_restart_var_file',          &
                          'State vector definition file not open', &
                          source, revision, revdate)
    end if
    
    ! Prepare storage
    read(unit=statevec_def_unit,fmt='(I8)',iostat=io_status) num_recs
    call check_io_status(io_status, routine, source, revision, revdate, &
                         'Read number of records in restart.vars')
    allocate(restart_info%restart_vars(num_recs), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision,    &
                            revdate, 'restart_info%restart_vars')

    ! Format is:
    !  1. Variable name - from the model
    !  2. Sigma Level
    !  3. Perturbation type (UNIFORM, INDIVID, NOPERTS)
    !  4. Perturbation magnitude
    !  5. Mass/W level (should be either 'M' or 'W')
    !  6. Variable type (from the raw observation kinds)
    !  7. Update field in restart vector or no?
    !  8. Whether field is a mean or not (TRUE/FALSE)
    !  9. Positive definiteness (ISPOSDEF, NOPOSDEF)
    do cur_rec = 1, num_recs
       read(unit=statevec_def_unit, fmt=*, iostat=io_status)      &
            name, index, type_name, pert_mag, lvl_type, var_type, &
            write_fld, mean_flag, pos_def

       ! We could gracefully continue if we hit an error, but bail
       ! instead
       call check_io_status(io_status, routine, source, revision, &
                            revdate, 'Reading restart.vars entry')
       
       ! Convert the text strings to their internally stored values
       call pert_name_to_type(type_name, pert_type)
       call level_type_to_flag(lvl_type, lvl_flag)
       call update_field_to_flag(write_fld, write_flag)
       call pos_def_to_flag(pos_def, pos_flag)
       var_type_num = get_index_for_quantity(var_type)

       if (var_type_num .lt. 0) then
          call error_handler(E_ERR, 'read_restart_var_file',            &
                             'Error in locating numeric equivalent ' // &
                             'to raw observation kind' // var_type,     & 
                             source, revision, revdate)
       end if

       call populate_info(name, index, pert_type, pert_mag,    &
                          var_type_num, write_flag, mean_flag, &
                          lvl_flag, pos_flag,                  &
                          restart_info%restart_vars(cur_rec))
    end do

  end subroutine read_restart_var_file

  ! level_type_to_flag
  ! ------------------
  ! Convert character representation of [M]ass or [W] velocity level
  ! type to an internal logical type 
  !  PARAMETERS
  !   IN  level_type_char      level type as character (M/W)
  !   OUT is_mass_level        True if levels is [M]ass
  subroutine level_type_to_flag(level_type_char, is_mass_level)
    character, intent(in)  :: level_type_char
    logical, intent(out)   :: is_mass_level

    select case (adjustl(level_type_char))
       case (adjustl(MASS_LEVEL))
          is_mass_level = .true.
       case (adjustl(W_LEVEL))
          is_mass_level = .false.
       case default
          call error_handler(E_ERR, 'level_type_to_flag', 'Level ' //    &
                             'type' // level_type_char // ' not found!', &
                             source, revision, revdate)
       end select
  end subroutine level_type_to_flag

  ! update_field_to_flag
  ! --------------------
  ! Convert character representation of whether the field should be
  ! written back into the COAMPS restart file (UPDATE) or skipped
  ! (FREEZE) to the internal logical representation
  !  PARAMETERS
  !   IN  update_string     UPDATE or FREEZE
  !   OUT update_flag       consistent with module-defined flag
  subroutine update_field_to_flag(update_string, update_flag)
    character(len=FIELD_UPDATE_LEN), intent(in) :: update_string
    logical, intent(out)                        :: update_flag

    select case (adjustl(update_string))
    case (adjustl(FIELD_UPDATE))
       update_flag = FLAG_UPDATE_FIELD
    case (adjustl(FIELD_FREEZE))
       update_flag = FLAG_FREEZE_FIELD
    case default
       call error_handler(E_ERR, 'update_field_to_flag',          &
                          'Update type' // update_string //       &
                          ' not found!', source, revision, revdate)
    end select
  end subroutine update_field_to_flag

  ! pos_def_to_flag
  ! ---------------
  ! Converts a character string describing whether the entire field
  ! is positive definite (ISPOSDEF) or not (NOPOSDEF) to the internal
  ! logical representation 
  !  PARAMETERS
  !   IN  posdef_string     ISPOSDEF or NOPOSDEF
  !   OUT posdef_flag       consistent with module-defined flag
  subroutine pos_def_to_flag(posdef_string, posdef_flag)
    character(len=POS_DEF_LEN), intent(in) :: posdef_string
    logical, intent(out)                   :: posdef_flag

    select case (adjustl(posdef_string))
    case (adjustl(ALLOW_NEGATIVE))
       posdef_flag = FLAG_ALLOW_NEGATIVE
    case (adjustl(FORCE_POSITIVE))
       posdef_flag = FLAG_FORCE_POSITIVE
    case default
       call error_handler(E_ERR, 'pos_def_to_flag', &
                          'Positive definite ' // posdef_string // &
                          ' not found!', source, revision, revdate)
    end select
  end subroutine pos_def_to_flag

  ! pert_name_to_type
  ! -----------------
  ! Converts a character string describing a perturbation type to the
  ! internal numeric representation
  !  PARAMETERS
  !   IN  pert_name         NOPERTS, INDIVID, UNIFORM
  !   OUT pert_type         consistent with module-defined values
  subroutine pert_name_to_type(pert_name, pert_type)
    character(len=PERT_TYPE_LEN), intent(in) :: pert_name
    integer, intent(out)                     :: pert_type

    select case (pert_name)
    case (NOPERTS_NAME)
       pert_type = PERT_TYPE_NOPERTS
    case (UNIFORM_NAME)
       pert_type = PERT_TYPE_UNIFORM
    case (INDIVID_NAME)
       pert_type = PERT_TYPE_INDIVID
    case default
       call error_handler(E_ERR, 'pert_name_to_type',             &
                          'Perturbation type ' // pert_name //    &
                          ' not found!', source, revision, revdate)
    end select
  end subroutine pert_name_to_type

  ! populate_info
  ! -------------
  ! Populates a restart_var structure with the information read in
  ! from the state vector definition file.  The dimension type and 
  ! variable record in that dimension are looked up given the name
  ! specified.
  !  PARAMETERS
  !   IN  var_name          the variable name
  !   IN  sig_index         sigma level index for this variable
  !   IN  pert_type         perturbation type (after conversion)
  !   IN  pert_magnitude    perturbation magnitude
  !   IN  var_type          the variable type (see the list of raw
  !                         types in obs_kind_mod definition)
  !   IN  update_flag       True if field is to be written back to
  !                         COAMPS restart file
  !   IN  mean_flag         True if field is a mean field
  !   IN  level_flag        True if field is a mass level
  !   IN  pos_flag          True if field is restricted to be 
  !                         non-negative
  !   OUT cur_var           restart_var structure containing the 
  !                         supplied information as well as the 
  !                         dimension type/location in restart file
  subroutine populate_info(var_name, sig_index, pert_type,         &
                           pert_magnitude, var_type, update_flag,  &
                           mean_flag, level_flag, pos_flag, cur_var) 
    character(len=VAR_NAME_LEN), intent(in)   :: var_name
    integer, intent(in)                       :: sig_index
    integer, intent(in)                       :: pert_type
    real(kind=r8), intent(in)                 :: pert_magnitude
    integer, intent(in)                       :: var_type
    logical, intent(in)                       :: update_flag
    logical, intent(in)                       :: mean_flag
    logical, intent(in)                       :: level_flag
    logical, intent(in)                       :: pos_flag
    type(restart_var), intent(out)            :: cur_var

    integer               :: cur_dim_type
    integer, dimension(2) :: cur_var_record

    ! These values we already have
    cur_var%var_name = var_name
    cur_var%sigma_record = sig_index
    cur_var%pert_mag = pert_magnitude
    cur_var%pert_type = pert_type
    cur_var%var_type = var_type
    cur_var%update_field = update_flag
    cur_var%mean_field = mean_flag
    cur_var%mass_level = level_flag
    cur_var%positive_definite = pos_flag

    ! Match the name up with the corresponding dimension/record
    ! information
    call get_name_info(DIM_TYPE_2D, DIM_TYPE_3D, DIM_TYPE_3DW,    &
                       SINGLEIO, MULTIIO, var_name, cur_dim_type, &
                       cur_var_record)
    cur_var%dim_type = cur_dim_type
    cur_var%var_record = cur_var_record
  end subroutine populate_info

  ! get_full_record_num
  ! -------------------
  ! Queries the given variable and returns the record number of a
  ! particular sigma level in the large COAMPS restart file
  !  PARAMETERS
  !   IN  cur_var           restart_var to pull data from
  !   IN  use_singleio      true if restart file is written using
  !                         a single I/O processor
  !   OUT recnum            variable's record number within the
  !                         restart file
  subroutine get_full_record_num(cur_var, use_singleio, recnum)
    type(restart_var), intent(in) :: cur_var
    logical, intent(in)           :: use_singleio
    integer, intent(out)          :: recnum

    ! How many total variables of each kind there are - use this for
    ! calculating offsets between dimensions
    integer :: N1D, N3D

    ! Temporary shorthand so we don't need to keep referencing the
    ! structure in the select statement
    integer :: kka, k, n

    ! How many of each kind of variable we have depends on if we are
    ! doing single or multiple process I/O
    if (use_singleio) then
       N1D = 101
       N3D = 64
    else
       N1D = 112
       N3D = 56
    end if

    kka = restart_info%grid%sigm_lvls
    k   = cur_var%sigma_record
    
    if (use_singleio) then
       n = cur_var%var_record(SINGLEIO)
    else
       n = cur_var%var_record(MULTIIO)
    end if

    ! Calculate offsets - do this by taking into account the number
    ! of variables that occur before the target variable and all
    ! their associated sigma levels, then all the sigma levels of
    ! THIS variable that have occured before the target level.
    select case (cur_var%dim_type)
    case (DIM_TYPE_2D)
       recnum = n
    case (DIM_TYPE_3D)
       recnum = (N1D) + ((kka)*(n-1) + 1) + (k - 1)
    case (DIM_TYPE_3DW)
       recnum = (N1D + N3D*kka) + ((kka + 1)*(n - 1) + 1) + (k - 1)
    case default
       call error_handler(E_ERR,'get_full_record_num',            &
                          'Unrecognized dimension type!', source, &
                          revision, revdate)
    end select
  end subroutine get_full_record_num

  ! get_posdef_flag
  ! ---------------
  ! Queries the given variable and returns whether the variable is
  ! constrained to never be negative (e.g. mixing ratios)
  !  PARAMETERS
  !   IN  cur_var           restart_var to pull data from
  !   OUT posdef            true if field shouldn't be negative
  subroutine get_posdef_flag(cur_var, posdef)
    type(restart_var), intent(in) :: cur_var
    logical, intent(out)          :: posdef

    posdef = cur_var%positive_definite
  end subroutine get_posdef_flag

  ! get_vert_coord
  ! --------------
  ! Queries the given variable and returns whether the variable's
  ! vertical coordinate type and vertical coordinate location.
  !  PARAMETERS
  !   IN  cur_var           restart_var to pull data from
  !   IN  grid              coamps_grid to use for coordinate calcs
  !   OUT vert_type         vertical coordinate type - currently
  !                         either surface or model level
  !   OUT vert_level        0 at surface or the sigma level
  subroutine get_vert_coord(cur_var, grid, vert_type, vert_level)
    type(restart_var), intent(in) :: cur_var
    type(coamps_grid), intent(in) :: grid
    integer, intent(out)          :: vert_type
    real(kind=r8), intent(out)    :: vert_level

    integer :: sigma_index

    sigma_index = cur_var%sigma_record

    ! Handle three cases - surface variables, on vertical mass levels
    ! and on vertical w levels
    if (sigma_index .eq. 0) then
       vert_type = VERTISSURFACE
       vert_level = 0.0
    else
       vert_type = VERTISLEVEL
       if (cur_var%mass_level) then
          vert_level = grid%msigma(sigma_index)
       else
          vert_level = grid%wsigma(sigma_index)
       end if
    end if
  end subroutine get_vert_coord

  ! get_update_flag
  ! ---------------
  ! Queries the given variable and returns whether the variable is
  ! written to the COAMPS restart file or not.
  !  PARAMETERS
  !   IN  cur_var           restart_var to pull data from
  !   OUT update            true if variable is written back to the
  !                         restart file
  subroutine get_update_flag(cur_var, update)
    type(restart_var), intent(in) :: cur_var
    logical, intent(out)          :: update

    update = cur_var%update_field
  end subroutine get_update_flag

  ! dump_restart_var
  ! ----------------
  ! Writes out the contents of a restart variable
  !  PARAMETERS
  !   IN  var               restart_var to dump contents of
  subroutine dump_restart_var(var)
    type(restart_var), intent(in) :: var

    character(len=32) :: kind_name

    kind_name = get_name_for_quantity(var%var_type)

    if (do_output()) then
       write (*, '(A9 T15A5)')   "VAR NAME:",var%var_name
       write (*, '(A9 T15I8.8)') "DIM TYPE:",var%dim_type
       write (*, '(A9 T15I8.8)') "VAR RCRD:",var%var_record
       write (*, '(A9 T15I8.8)') "SGM RCRD:",var%sigma_record
       write (*, '(A9 T15F8.6)') "PTRB PCT:",var%pert_mag
       write (*, '(A9 T15I8.8)') "PTB TYPE:",var%pert_type
       write (*, '(A9 T15L1)')   "UPDATE??:",var%update_field
       write (*, '(A9 T15A32)')  "VAR TYPE:",kind_name
       write (*, '(A9 T15I8.8)') "VAR NMBR:",var%var_type
       write (*, '(A9 T15L1)')   "MASSLVL?:",var%mass_level
       write (*, '(A9 T15L1)')   "MEANFLD?:",var%mean_field
       write (*, *) "-----------------------------------------"
    end if
  end subroutine dump_restart_var

  !------------------------------
  ! END PRIVATE ROUTINES
  !------------------------------
end module coamps_restart_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
