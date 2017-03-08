! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

!------------------------------
! MODULE:       coamps_statevec_mod
! AUTHOR:       T. R. Whitcomb
!               Naval Research Laboratory
! DART VERSION: Jamaica
!
! Module containing the data structure and routines for dealing with
! a COAMPS state vector (a collection of COAMPS state variables)
!------------------------------ 

module coamps_statevec_mod

    use coamps_statevar_mod, only : state_variable, new_state_variable,    &
                                    read_state_variable, set_position,     &
                                    set_2d_flag, dump_state_variable,      &
                                    set_var_stagger,                       &
                                    get_mean_flag, define_mean_var,        &
                                    get_sigma_record, set_sigma_record,    &
                                    get_mass_level_flag,                   &
                                    operator(==)

    use coamps_domain_mod,   only : coamps_domain, get_domain_num_levels, &
                                    get_nest_count

    use coamps_nest_mod,     only : coamps_nest, get_nest_id

    use coamps_util_mod,     only : check_alloc_status,  &
                                    check_io_status

    use types_mod,           only :  r8

    use obs_kind_mod

    use utilities_mod,       only : do_output,     &
                                    E_WARN,        &
                                    E_ERR,         &
                                    error_handler, &
                                    get_unit

    implicit none

    private

    !------------------------------
    ! BEGIN PUBLIC INTERFACE
    !------------------------------

    public :: state_vector
    public :: initialize_state_vector

    public :: get_file_type

    ! Iterator structure and methods
    public :: state_iterator    
    public :: get_iterator
    public :: has_next
    public :: get_next

    ! Diagnostics
    public :: dump_state_vector

    ! Search functions
    public :: find_state_variable

    public :: get_num_fields
    public :: get_total_size
    public :: get_domain
    public :: get_var_by_index
    

    !------------------------------
    ! END PUBLIC INTERFACE
    !------------------------------

    !------------------------------
    ! BEGIN EXTERNAL INTERFACE
    !------------------------------
    ! [none]
    !------------------------------
    ! END EXTERNAL INTERFACE
    !------------------------------
  
    !------------------------------
    ! BEGIN TYPES AND CONSTANTS
    !------------------------------

    ! Derived data type that contains the information used to read and
    ! store the field information for the state vector
    type :: state_vector
        private

        logical                                     :: FLAT_FILE_IO
        integer                                     :: num_fields
        integer                                     :: cur_fld_cnt

        type(coamps_domain),                pointer :: domain
        type(state_variable), dimension(:), pointer :: vars
        integer                                     :: total_size

    end type state_vector
  
    ! Iterate through the restart vector variables
    type :: state_iterator
        private

        type(state_vector), pointer :: state
        integer                     :: cur_index
    end type state_iterator


    ! Flat and restart files i/o constants 
    integer,                      parameter :: FILE_TYPE_LEN       = 12
    character(len=FILE_TYPE_LEN), parameter :: FLAT_FILE_STRING    = 'FLAT_FILE'
    character(len=FILE_TYPE_LEN), parameter :: RESTART_FILE_STRING = 'RESTART_FILE'
    logical,                      parameter :: FLAG_FLAT           = .true.
    logical,                      parameter :: FLAG_RESTART        = .false.

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

    ! initialize_state_vector
    ! -----------------------
    ! Populate a COAMPS state vector structure by reading in the grid
    ! data and reading the state vector definition file
    ! INOUT state           COAMPS state vector
    ! INOUT state_list      COAMPS state vector 
    !   IN  filename        name of state vector definition file
    !   IN  domain          domain for this run
    subroutine initialize_state_vector(state, filename, domain, list_only)
        type(state_vector),          intent(inout) :: state
        character(len=*),            intent(in)    :: filename
        type(coamps_domain), target, intent(in)    :: domain
        logical, optional,           intent(in)    :: list_only

        nullify(state%domain)
        nullify(state%vars)

        state%domain => domain
        state%cur_fld_cnt = 0

        if(present(list_only)) then
          call read_state_vector(state, filename, list_only)
        else
          call read_state_vector(state, filename)
        end if
        call calculate_state_positions(state)

    end subroutine initialize_state_vector

    ! get_file_type
    ! ------------
    ! Returns logical specifying restart or flat files
    !  PARAMETERS
    !  OUT file_type  Flat files or restart files?
    function get_file_type(state) result(file_type)
      type(state_vector), intent(in) :: state
      logical                        :: file_type
      file_type = state%FLAT_FILE_IO
    end function get_file_type

    ! dump_state_vector
    ! -----------------
    ! Write out all the components of the state vector
    !  PARAMETERS
    !   IN  state           State vector to dump
    subroutine dump_state_vector(state)
        type(state_vector), intent(in) :: state

        type(state_iterator) :: iterator

        if (do_output()) then
            write (*,*) "*** State Vector Contents ***"
            write (*,*) "-----------------------------"

            iterator = get_iterator(state)    
            do while (has_next(iterator))
                call dump_state_variable(get_next(iterator))
            end do
        end if
    end subroutine dump_state_vector
    
    ! find_state_variable
    ! -------------------
    ! Search through the state vector for a variable matching the attributes
    ! supplied.  Return the first state variable found to match or a null 
    ! pointer if no match was found.
    !  PARAMETERS
    !   IN  state           COAMPS state vector to search
    !   IN  nest            COAMPS nest the variable is defined on 
    !   IN  var_kind        integer form of variable type
    !   IN  is_mean         true if variable is a mean field
    !   IN  level_type      level of variable to find
    !   IN  sigma_index     sigma level index for variable
    function find_state_variable(state, nest, var_kind, is_mean, &
                                 level_type, sigma_index, vert_value)
        type(state_vector),            intent(in)  :: state
        type(coamps_nest),             intent(in)  :: nest    
        integer,                       intent(in)  :: var_kind
        logical,                       intent(in)  :: is_mean
        character(len=*),              intent(in)  :: level_type
        integer,                       intent(in)  :: sigma_index
        real(kind=r8), optional,       intent(in)  :: vert_value
        type(state_variable), pointer              :: find_state_variable

        type(state_iterator)         :: iterator
        type(state_variable), target :: test_var
        type(state_variable)         :: compare_var

        character(len=128) :: message
        character(len=*), parameter :: routine = 'find_state_variable'

        if(present(vert_value)) then
          compare_var = new_state_variable(get_nest_id(nest), var_kind, &
                                           is_mean, level_type,         &
                                           sigma_index,                 &
                                           vert_value = vert_value)
        else
          compare_var = new_state_variable(get_nest_id(nest), var_kind, &
                                           is_mean, level_type,         &
                                           sigma_index)
        end if

        ! Default assumption is no match found
        nullify(find_state_variable)

        iterator = get_iterator(state)
        do while (has_next(iterator))
            test_var = get_next(iterator)

            ! Return the first match found
            if (compare_var == test_var) then
                find_state_variable => test_var
                return
            end if
        end do

        ! If we got here, the search failed
        write (message,'(3A,A1,A,I3,L1)') "Could not find  ",    &
                          trim(get_name_for_quantity(var_kind)), &
                          " on ", level_type,"-level.", sigma_index, is_mean
    end function find_state_variable

    ! get_num_fields
    ! ------------
    ! Returns how many fields there are in the state vector
    !  PARAMETERS
    !   OUT num_vars          number of fields in state vector
    function get_num_fields(state) result(num_fields)
      type(state_vector), intent(in) :: state
      integer                        :: num_fields
      num_fields = state%num_fields
    end function get_num_fields

    ! get_iterator
    ! ------------
    ! Return an iterator for the state vector that can be used to access
    ! each element in turn
    !  PARAMETERS
    !   IN  state           State vector to iterate over
    function get_iterator(state)
        type(state_vector), target, intent(in)  :: state
        type(state_iterator)                    :: get_iterator

        get_iterator%state     => state
        get_iterator%cur_index =  1
    end function get_iterator

    ! has_next
    ! --------
    ! Return true if the given iterator has another value to give 
    !  PARAMETERS
    !   IN  iterator        The iterator to check
    function has_next(iterator)
        type(state_iterator), intent(in)  :: iterator
        logical                           :: has_next

        has_next = (iterator%cur_index <= size(iterator%state%vars))
    end function has_next

    ! get_next
    ! --------
    ! Returns the current state variable based on the given iterator 
    ! **AND** advance to the next position
    !  PARAMETERS
    !   IN  iterator        The iterator to grab a value from
    function get_next(iterator)
        type(state_iterator),          intent(inout) :: iterator
        type(state_variable), pointer                :: get_next

        ! Return the current value and advance to the next
        get_next           => iterator%state%vars(iterator%cur_index)
        iterator%cur_index = iterator%cur_index + 1

    end function get_next

    ! get_total_size
    ! --------------
    ! Get the total length of the state vector based on the variable
    ! definition and the domain.
    !  PARAMETERS
    !   IN  state           State vector to query
    function get_total_size(state)
        type(state_vector), intent(in)  :: state
        integer                         :: get_total_size

        get_total_size = state%total_size
    end function get_total_size

    ! get_domain
    ! ----------
    ! Returns the domain for a given state vector definition
    !  PARAMETERS
    !   IN  state           State vector definition to query
    function get_domain(state)
        type(state_vector),  intent(in)  :: state
        type(coamps_domain)              :: get_domain

        get_domain = state%domain
    end function get_domain

    ! get_var_by_index
    ! ----------
    ! Returns the var for a given index in the state vector
    !  PARAMETERS
    !   IN  state           State vector definition to query
    function get_var_by_index(state, index_in)
        type(state_vector),   intent(in) :: state
        integer,              intent(in) :: index_in
        type(state_variable), pointer    :: get_var_by_index

        get_var_by_index => state%vars(index_in)
    end function get_var_by_index

    !------------------------------
    ! END PUBLIC ROUTINES
    !------------------------------

    !------------------------------
    ! BEGIN PRIVATE ROUTINES
    !------------------------------

    ! add_to_list
    ! -------------------
    ! Adds the state_variable to the list
    !  PARAMETERS
    !   IN  var       Integer to add to the list
    ! INOUT list      List to add the entry to
    subroutine add_to_list(var, list)
      type(state_variable),  intent(in)    :: var
      type(state_vector),    intent(inout) :: list
      integer                              :: cur_index

      list%cur_fld_cnt = list%cur_fld_cnt + 1
      list%vars(list%cur_fld_cnt) = var

    end subroutine add_to_list

    ! read_state_vector
    ! -----------------
    ! Reads in the state vector definition file - the first line 
    ! is the number of variables to follow, then each line after that
    ! contains information on each entry within the state vector
    !  PARAMETERS
    ! INOUT state                   The state vector to populate
    !   IN  state_definition_file   The definition file name to read from 
    subroutine read_state_vector(state, state_definition_file, list_only)
        type(state_vector), intent(inout) :: state
        character(len=*),   intent(in)    :: state_definition_file
        logical, optional,  intent(in)    :: list_only

        integer :: state_definition_unit

        type(state_iterator) :: vars_iterator
        type(state_variable) :: cur_var

        character(len=*), parameter     :: routine = 'read_state_vector'
        integer                         :: io_status
        integer                         :: alloc_status

        integer                         :: num_lvls
        integer                         :: num_flds
        integer                         :: num_nests
        integer                         :: num_recs

        integer                         :: cur_lvl
        integer                         :: cur_fld
        integer                         :: cur_nest

        character(len=FILE_TYPE_LEN)    :: file_type
        logical                         :: list_only_local

        ! variables to deal with special case mean fields
        integer, parameter                                       :: NUM_DEFAULT_VARS = 3
        integer, parameter                                       :: VAR_NAME_LEN     = 4
        character(len=VAR_NAME_LEN), dimension(NUM_DEFAULT_VARS), parameter :: &
                  default_vars = (/'THBM','EXBM','EXBW'/)

        num_nests = get_nest_count(state%domain)

        if(present(list_only)) then
          list_only_local = list_only
        else
          list_only_local = .false.
        end if

        state_definition_unit = get_unit()
        open(unit=state_definition_unit, file=state_definition_file, &
             status='old', access='sequential', action='read',       &
             form='formatted', iostat=io_status)
        call check_io_status(io_status, routine, source, revision, revdate, &
                             'Opening file ' // state_definition_file)

        read(unit=state_definition_unit,fmt='(A)',iostat=io_status) file_type
        call check_io_status(io_status, routine, source, revision, revdate, &
                             'Read file type in state.vars')
        state%FLAT_FILE_IO = file_type_to_flag(file_type) 

        ! Count the number of itemes in the state vector file
        num_recs = 0
        count_recs: do 
          if(read_state_variable(state_definition_unit, cur_var,  &
                                 0, state%FLAT_FILE_IO)) exit count_recs

          if(state%FLAT_FILE_IO .AND. .not.list_only_local) then
            num_recs = num_recs + get_sigma_record(cur_var)
          else
            num_recs = num_recs + 1
          end if
        end do count_recs

        ! Account for 3 mean-state variable if this is flatfile i/o
        if(state%FLAT_FILE_IO)  then
          if(list_only_local) then
            num_recs = num_nests * (num_recs + 3)
          else
            num_recs = num_nests * (num_recs                         &
                     + 3 * get_domain_num_levels(state%domain) + 1)
          end if
        end if

        ! Rewind to the begining of the file 
        rewind(state_definition_unit)
        read(unit=state_definition_unit,fmt='(A)',iostat=io_status) file_type

        allocate(state%vars(num_recs), stat=alloc_status)
        call check_alloc_status(alloc_status, routine, source, revision,    &
                                revdate, 'state%vars')

        state%num_fields = 0
        nest_input: do cur_nest = 1, num_nests
          read_list: do 

            if(read_state_variable(state_definition_unit, cur_var,  &
                                   cur_nest, state%FLAT_FILE_IO)) exit read_list

            if(get_mean_flag(cur_var)) cycle read_list

            state%num_fields = state%num_fields + 1

            if(list_only_local) then
              if(get_sigma_record(cur_var) > 1) call set_2d_flag(cur_var, .false.)
              call set_var_stagger(cur_var)
              call add_to_list(cur_var, state)
              cycle read_list
            end if
              
            if(state%FLAT_FILE_IO) then
              num_lvls = get_sigma_record(cur_var)
              levels: do cur_lvl=1, num_lvls
                call set_sigma_record(cur_var, cur_lvl)
                call add_to_list(cur_var, state)
              end do levels
            else
              call add_to_list(cur_var, state)
            end  if

          end do read_list

          ! Add mean fields to the list.  Only for flat_file i/o, otherwise
          ! variables are defined in state_vector text file.
          if(state%FLAT_FILE_IO) then
            fld_loop: do cur_fld = 1, NUM_DEFAULT_VARS

              cur_var = define_mean_var(default_vars(cur_fld), cur_nest)
              num_lvls = get_domain_num_levels(state%domain)
              if(.not.get_mass_level_flag(cur_var)) num_lvls = num_lvls + 1

              if(list_only_local) then
                call set_sigma_record(cur_var, num_lvls)
                call set_2d_flag(cur_var, .false.)
                call set_var_stagger(cur_var)
                call add_to_list(cur_var, state)
              else
                lvl_loop: do cur_lvl = 1, num_lvls
                  call set_sigma_record(cur_var, cur_lvl)
                  call add_to_list(cur_var, state)
                end do lvl_loop
              end if

            end do fld_loop
          end if

          ! Rewind the file for the next nest
          rewind(state_definition_unit)
          read(unit=state_definition_unit,fmt='(A)',iostat=io_status) file_type
        end do nest_input

        close(state_definition_unit)

    end subroutine read_state_vector

    ! calculate_state_positions
    ! -------------------------
    ! Iterate through the state vector definition and calculate the piece of
    ! the large state vector corresponding to each state variable.    
    !  PARAMETERS
    !   IN  state         state vector definition to iterate over 
    subroutine calculate_state_positions(state)
        type(state_vector), intent(inout) :: state

        type(state_iterator)          :: iterator
        type(state_variable), pointer :: cur_var
        integer                       :: cur_state_pos

        cur_state_pos = 0
        iterator = get_iterator(state)
        do while (has_next(iterator))
          cur_var => get_next(iterator) 
          cur_state_pos = set_position(cur_var, cur_state_pos, state%domain)
        end do

        state%total_size = cur_state_pos
    end subroutine calculate_state_positions

    ! file_type_to_flag
    ! ---------------
    ! Converts a character string describing whether the entire field
    ! is read from a flat file (FLAT) or a restart file (RESTART) to 
    ! the interna llogical representation 
    !  PARAMETERS
    !   IN  file_type_string     FLAT or RESTART
    !   OUT flat_flag       consistent with module-defined flag
    function file_type_to_flag(file_type_string) result(file_type_flag)
      character(len=FILE_TYPE_LEN), intent(in) :: file_type_string
      logical                                  :: file_type_flag
      select case (adjustl(file_type_string))
      case (adjustl(FLAT_FILE_STRING))
         file_type_flag = FLAG_FLAT
      case (adjustl(RESTART_FILE_STRING))
         file_type_flag = FLAG_RESTART
      case default
         call error_handler(E_ERR, 'file_type_to_flag', &
                            'FILE_TYPE ' // file_type_string // &
                            ' not found!', source, revision, revdate)
      end select
    end function file_type_to_flag
    
    !------------------------------
    ! END PRIVATE ROUTINES
    !------------------------------
end module coamps_statevec_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
