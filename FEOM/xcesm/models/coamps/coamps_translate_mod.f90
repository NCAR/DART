! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

module coamps_translate_mod

!------------------------------
! MODULE:       coamps_translate_mod
! AUTHOR:       T. R. Whitcomb
!               Naval Research Laboratory
! DART VERSION: Jamaica
!
! Module containing storage for a DART state vector, a COAMPS state
! vector (spanning multiple restart files in the case of multi-
! processor I/O), and routines to read/write them and convert
! one to another.
!
! The DART time structure is represented as days and seconds, and
! I've included several routines in here to deal with converting
! the time to hours/minutes/seconds with no respect to dates 
! (since we're restricting to a single date-time group).
!------------------------------ 

  use coamps_grid_mod,    only : check_ij_within_grid,         &
                                 coamps_grid,                  &
                                 decompose_domain,             &
                                 get_grid_field_size,          &
                                 get_grid_dims,                &
                                 get_grid_iminf,               &
                                 get_grid_imaxf,               &
                                 get_grid_jminf,               &
                                 get_grid_jmaxf,               &
                                 grid_ij_to_vector_index,      &
                                 initialize_decomposition
  use coamps_restart_mod, only : dump_restart_vars,            &
                                 get_num_vars,                 &
                                 get_full_record_num_by_index, &
                                 get_restart_grid,             &
                                 get_update_flag_by_index,     &
                                 get_var_info_by_abs_index,    &
                                 get_posdef_flag_by_index,     &
                                 initialize_restart_info
  use coamps_util_mod,    only : C_REAL,                       &
                                 check_alloc_status,           &
                                 check_dealloc_status,         &
                                 check_io_status,              &
                                 fix_for_platform
  use time_manager_mod,   only : get_date,                     &
                                 get_time,                     &
                                 GREGORIAN,                    &
                                 operator(-),                  &
                                 operator(+),                  &
                                 read_time,                    &
                                 set_calendar_type,            &
                                 set_date,                     &
                                 set_time,                     &
                                 time_type,                    &
                                 write_time
  use types_mod,          only : r8
  use utilities_mod,      only : E_ERR,                        &
                                 error_handler,                &
                                 file_exist,                   &
                                 get_unit

  implicit none

  private

  !------------------------------
  ! BEGIN PUBLIC INTERFACE
  !------------------------------

  public :: initialize_translator
  public :: finalize_translator

  ! COAMPS restart file tools
  public :: generate_restart_filenames
  public :: open_coamps_files
  public :: coamps_restart_read_all_fields
  public :: coamps_restart_write_all_fields

  ! DART restart file tools
  public :: open_dart_file
  public :: dart_read
  public :: dart_write

  ! Conversion tools
  public :: convert_dart_state_to_coamps
  public :: convert_coamps_state_to_dart
  public :: fix_negative_values

  ! Time handling tools
  public :: set_dart_current_time
  public :: get_dart_current_time
  public :: get_dart_target_time
  public :: write_pickup_file

  !FIX ONLY FOR TESTING
  public :: dart_state

  !------------------------------
  ! END PUBLIC INTERFACE
  !------------------------------

  !------------------------------
  ! BEGIN EXTERNAL INTERFACE
  !------------------------------
  !  [none]
  !------------------------------
  ! END EXTERNAL INTERFACE
  !------------------------------

  !------------------------------
  ! BEGIN TYPES AND CONSTANTS
  !------------------------------

  ! Maximum number of nests that could appear in the namelist - we
  ! only care about the first one though but if there's more numbers
  ! in the namelist than this we'll get an error.
  integer, parameter :: MAXNEST=10


  ! Conversion parameters
  integer, parameter :: MIN_TO_SEC = 60
  integer, parameter :: HR_TO_MIN  = 60
  integer, parameter :: DAY_TO_HR  = 24
  integer, parameter :: HR_TO_SEC  = HR_TO_MIN * MIN_TO_SEC
  integer, parameter :: DAY_TO_SEC = DAY_TO_HR * HR_TO_SEC

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

  ! Namelist containing the date time group and lead time information
  ! Also include information about the domain decomposition: number
  ! of processors in x and y, the number of I/O processors, and the
  ! number of halo boundary points.
  integer, dimension(3)        :: ktauf    = (/ 1,0,0 /)
  integer, dimension(3)        :: ktaust   = (/ 0,0,0 /)
  integer, dimension(MAXNEST)  :: ndxnam   = (/ 1,1,1,1,1,1,1,1,1,1 /)
  integer, dimension(MAXNEST)  :: ndynam   = (/ 1,1,1,1,1,1,1,1,1,1 /) 
  integer                      :: npr0nam  = 1
  integer                      :: nbnam    = 2
  character(len=10)            :: cdtg     = '1999083100'
  character(len=11)            :: nml_file = 'convert.nml'
  namelist /convert/ ktauf, ktaust, ndxnam, ndynam, npr0nam, nbnam,&
       & cdtg

  ! Contents of DART state vector
  type(time_type), dimension(2)            :: dart_time
  real(kind=r8), dimension(:), allocatable, target :: dart_state

  ! Contents of COAMPS state vector - this is not the full restart
  ! file - only the fields that we have defined in the dynamic
  ! state vector definition 
  real(kind=C_REAL), dimension(:), allocatable, target :: coamps_state
  
  ! COAMPS grid details - need this for information about the domain
  ! size and the domain decomposition routines
  type(coamps_grid) :: grid

  ! Arrays allow reading/writing multiple COAMPS restart files 
  integer                                      :: total_coamps_files
  character(len=36), dimension(:), allocatable :: coamps_file_names
  integer, dimension(:), allocatable           :: coamps_file_units
  
  ! DART restart file
  character(len=11), parameter :: DART_FILENAME = 'dart_vector'
  integer                      :: dart_unit              

  !------------------------------
  ! END MODULE VARIABLES
  !------------------------------

contains

  !------------------------------
  ! BEGIN PUBLIC ROUTINES
  !------------------------------

  ! initialize_translator
  ! ---------------------
  ! Module "constructor" that calls the internal initialization 
  ! routines in their proper order 
  !  PARAMETERS
  !   [none]
  subroutine initialize_translator()

    ! Get the data we can't read from pre-exiting COAMPS files
    call read_convert_namelist()

    ! Set up the state vector field definitions 
    call initialize_restart_info(cdtg, 'restart.vars')

    ! Need domain size so we know how much memory we need
    call get_restart_grid(grid)

    call allocate_state()
    call allocate_coamps_file_info()

    ! Figure out what points go with what processes and hence which
    ! restart files - currently only support a single nest
    call initialize_decomposition(grid,ndxnam(1)*ndynam(1))
    call decompose_domain(grid, npr0nam, ndxnam(1), ndynam(1), nbnam)
  end subroutine initialize_translator

  ! finalize_translator
  ! -------------------
  ! Closes files and frees allocated memory
  !  PARAMETERS
  !   [none] 
  subroutine finalize_translator()
    integer :: cur_file 

    character(len=*), parameter :: routine = 'finalize_translator'
    integer                     :: dealloc_status

    deallocate(dart_state,   stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source, &
                              revision, revdate, 'dart_state'  )
    deallocate(coamps_state, stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source, &
                              revision, revdate, 'coamps_state')

    do cur_file = 1, total_coamps_files
       close(coamps_file_units(cur_file))
    end do
    close(dart_unit)

    deallocate(coamps_file_names, stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source,       &
                              revision, revdate, 'coamps_file_names' )
    deallocate(coamps_file_units, stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source,       &
                              revision, revdate, 'coamps_file_units' )
  end subroutine finalize_translator
  
  ! generate_restart_filenames
  ! --------------------------
  ! Generate the COAMPS restart file names based on information from
  ! the convert namelist.
  !  PARAMETERS
  !   IN  writing_coamps    True if we are going to be writing the
  !                         COAMPS restart file, false if reading it
  subroutine generate_restart_filenames(writing_coamps)
    logical, intent(in) :: writing_coamps

    integer               :: coamps_file_index
    integer, dimension(3) :: working_time

    ! Assert that the restart file names are already allocated
    if (.not. allocated(coamps_file_names)) then
       call error_handler(E_ERR, 'generate_restart_filenames', &
                         'COAMPS file names unallocated',      &
                         source, revision, revdate)
    end if

    ! Need to decide what time we are going to use.
    ! If we are going from COAMPS -> DART, need the COAMPS final
    ! forecast time from the convert namelist.  If we are going 
    ! from DART -> COAMPS, use the initial forecast time from the
    ! DART restart file.
    if (writing_coamps) then
       working_time = ktaust
    else
       working_time = ktauf
    end if

    print *,"COAMPS restart files in use:"
    print *,"-------------------------------------------------"
    do coamps_file_index = 1,total_coamps_files
       if (npr0nam > 0) then
          write (coamps_file_names(coamps_file_index), 100) &
                cdtg, working_time
       else
          write (coamps_file_names(coamps_file_index), 200) &
                (coamps_file_index-1), cdtg, working_time
       end if
       write (*,*) coamps_file_index,                       &
             coamps_file_names(coamps_file_index)
    end do
    print *,"-------------------------------------------------"
100 format( 'restarta1p001'A10I3.3I2.2I2.2'.nest1')
200 format( 'restarta1p'I3.3A10I3.3I2.2I2.2)
  end subroutine generate_restart_filenames

  ! open_coamps_files
  ! -----------------
  ! Uses the list of COAMPS file names to populate the vector of
  ! opened COAMPS restart file units
  ! PARAMETERS
  !  IN  writing_coamps     True if we're opening the files for
  !                         write access
  subroutine open_coamps_files(writing_coamps)
    logical, intent(in) :: writing_coamps

    integer           :: cur_file_index, cur_file_unit

    character(len=*), parameter :: routine = 'open_coamps_files'
    integer           :: io_status

    character(len=5)  :: fileaction

    if (writing_coamps) then
       fileaction = 'write'
    else
       fileaction = 'read'
    end if

    do cur_file_index = 1, total_coamps_files
       cur_file_unit = get_unit()
       open( unit=cur_file_unit,                                &
             file=coamps_file_names(cur_file_index),            &
             status='old', access='direct', action=fileaction,  &
             form='unformatted', recl=C_REAL, iostat=io_status)
       call check_io_status(io_status, routine, source, revision, &
                            revdate, 'Opening ' //                &
                            coamps_file_names(cur_file_index))

       coamps_file_units(cur_file_index) = cur_file_unit
    end do
  end subroutine open_coamps_files

  ! coamps_restart_read_all_fields
  ! ------------------------------
  ! Wrapper for using coamps_process_all_fields to read COAMPS
  ! restart files
  !  PARAMETERS
  !   [none] 
  subroutine coamps_restart_read_all_fields()
    logical, parameter :: READING_COAMPS = .false.
    call coamps_process_all_fields(READING_COAMPS)
  end subroutine coamps_restart_read_all_fields

  ! coamps_restart_write_all_fields
  ! -------------------------------
  ! Wrapper for using coamps_process_all_fields to write COAMPS
  ! restart files
  !  PARAMETERS
  !   [none]
  subroutine coamps_restart_write_all_fields()
    logical, parameter :: WRITING_COAMPS = .true.
    call coamps_process_all_fields(WRITING_COAMPS)
  end subroutine coamps_restart_write_all_fields

  ! open_dart_file
  ! -----------------
  ! Opens the DART restart file for reading or writing
  !  PARAMETERS
  !   IN  writing_dart      True if we are opening the DART file
  !                         for write access
  subroutine open_dart_file(writing_dart)
    logical, intent(in) :: writing_dart

    character(len=5) :: fileaction
    character(len=7) :: filestatus

    character(len=*), parameter :: routine = 'open_dart_file'
    integer :: io_status

    if (writing_dart) then
       fileaction = 'write'
       filestatus = 'replace'
    else
       fileaction = 'read'
       filestatus = 'old'
    end if

    dart_unit = get_unit()
    open( unit=dart_unit, file=DART_FILENAME, status=filestatus,   & 
          action=fileaction, form='unformatted', iostat=io_status)
    call check_io_status(io_status, routine, source, revision, revdate, &
                         'Opening ' // DART_FILENAME)
  end subroutine open_dart_file

  ! convert_dart_state_to_coamps
  ! ----------------------------
  ! Wrapper for DART -> COAMPS convert_state
  !  PARAMETERS
  !   [none]
  subroutine convert_dart_state_to_coamps()
    logical, parameter :: CONVERT_TO_COAMPS = .false.
    call convert_state(CONVERT_TO_COAMPS)
  end subroutine convert_dart_state_to_coamps

  ! convert_coamps_state_to_dart
  ! ----------------------------
  ! Wrapper for COAMPS -> DART convert_state
  !  PARAMETERS
  !   [none]
  subroutine convert_coamps_state_to_dart()
    logical, parameter :: CONVERT_TO_DART  = .true.
    call convert_state(CONVERT_TO_DART)
  end subroutine convert_coamps_state_to_dart

  ! dart_read
  ! ---------
  ! Reads two times and the DART state vector from the DART restart
  ! file.  The dart_(read|write) routines are separate unlike the
  ! COAMPS routines since this involves no offset math at all.
  !  PARAMETERS
  !   IN  test_read         OPTIONAL: True if we just want to read
  !                                   the state vector, not the 
  !                                   times (useful for tests) 
  subroutine dart_read(test_read)
    logical, optional, intent(in) :: test_read

    logical :: read_times

    integer :: dart_time_days, dart_time_seconds

    ! Performing a test read means not reading the times
    if (present(test_read)) then
       read_times = .not. test_read
    else
       read_times = .true. 
    end if

    if (read_times) then
       ! dart_time(1) is the time DART wants the model to advance to
       dart_time(1) = read_time(dart_unit, 'unformatted')
       call get_time(dart_time(1), dart_time_seconds, dart_time_days)
       write (*,*) "DART target time read is ", dart_time_days, &
                   " days and ", dart_time_seconds, " seconds."

       ! dart_time(2) is the current time
       dart_time(2) = read_time(dart_unit, 'unformatted')
       call get_time(dart_time(2), dart_time_seconds, dart_time_days)
       write (*,*) "DART current time read is ", dart_time_days, &
                   " days and ", dart_time_seconds, " seconds."
    end if
  
    read(dart_unit) dart_state

    call print_dart_diagnostics()
  end subroutine dart_read

  ! dart_write
  ! ---------
  ! Writes the time and the DART state vector to the DART restart
  ! file.  The dart_(read|write) routines are separate unlike the
  ! COAMPS since these involve no offset math at all.
  !  PARAMETERS
  !   IN  test_write        OPTIONAL: True if we just want to write
  !                                   out the state vector, not the
  !                                   time (useful for tests)
  subroutine dart_write(test_write)
    logical, optional, intent(in) :: test_write
    
    logical :: write_times

    integer :: dart_time_days, dart_time_seconds

    ! Performing a test write means not writing the time
    if (present(test_write)) then
       write_times = .not. test_write
    else
       write_times = .true.
    end if

    rewind(dart_unit)

    if (write_times) then
       call write_time(dart_unit, dart_time(1), 'unformatted')
       call get_time(dart_time(1), dart_time_seconds, dart_time_days)
       write (*,*) "DART time written is ", dart_time_days, &
                   " days and ", dart_time_seconds, " seconds."
    end if

    write(dart_unit) dart_state

    call print_dart_diagnostics()
  end subroutine dart_write

  ! set_dart_current_time
  ! ---------------------
  ! When converting from COAMPS to DART, the forecast is finished
  ! and the current time is the final forecast time of the model
  ! run - take that hours/minutes/seconds time and convert it to 
  ! the DART days/seconds format.
  !  PARAMETERS
  !   [none]
  subroutine set_dart_current_time()
    integer :: tau_hour
    integer :: tau_minute
    integer :: tau_second
    integer :: dart_second
    integer :: dart_day

    ! Get this information from the namelist
    tau_hour   = ktauf(1)
    tau_minute = ktauf(2)
    tau_second = ktauf(3)

    call hms_to_sd(tau_hour, tau_minute, tau_second, dart_day,&
         & dart_second)

    ! When writing a DART file, only have one time entry to worry
    ! about.
    dart_time(1) = set_time(dart_second, dart_day)
  end subroutine set_dart_current_time

  ! get_dart_current_time
  ! ---------------------
  ! When converting from DART to COAMPS, we get *two* times - the
  ! current time and the target time that it wants the model to
  ! advance to - take the current time in DART days/seconds format
  ! and convert it to the COAMPS hour/minute/second format.
  !  PARAMETERS
  !   [none]
  subroutine get_dart_current_time()
    integer :: cur_hour
    integer :: cur_minute
    integer :: cur_second
    integer :: dart_days
    integer :: dart_seconds

    ! Note that the *second* entry in the array is the current time
    call get_time(dart_time(2), dart_seconds, dart_days)
    call sd_to_hms(dart_days, dart_seconds, cur_hour, cur_minute,&
         & cur_second) 

    ktaust(1) = cur_hour
    ktaust(2) = cur_minute
    ktaust(3) = cur_second
  end subroutine get_dart_current_time

  ! get_dart_target_time
  ! ---------------------
  ! When converting from DART to COAMPS, we get *two* times - the
  ! current time and the target time that it wants the model to
  ! advance to - take the target time in DART days/seconds format
  ! and convert it to the COAMPS hour/minute/second format
  !  PARAMETERS
  !   [none]
  subroutine get_dart_target_time()
    integer :: tgt_hour
    integer :: tgt_minute
    integer :: tgt_second
    integer :: dart_days
    integer :: dart_seconds

    ! Note that the *first* entry in the array is the target time
    call get_time(dart_time(1), dart_seconds, dart_days)
    call sd_to_hms(dart_days, dart_seconds, tgt_hour, tgt_minute,&
         & tgt_second) 
    ktauf(1) = tgt_hour
    ktauf(2) = tgt_minute
    ktauf(3) = tgt_second
  end subroutine get_dart_target_time

  ! write_pickup_file
  ! -----------------
  ! Writes a file to disk containing the values of the DART-based
  ! ktaust and ktauf - these will be picked up by our run script and
  ! placed in the namelist to advance the forecast
  ! PARAMETERS
  !  [none]
  subroutine write_pickup_file()

    character(len=*), parameter :: PICKUP_NAME = 'dart.kstart'
    integer                     :: pickup_unit

    character(len=*), parameter :: routine = 'write_pickup_file' 
    integer                     :: io_status
    
    pickup_unit = get_unit()
    
    open( unit=pickup_unit, file=PICKUP_NAME, status='replace',         &
          form='formatted', iostat=io_status)
    call check_io_status(io_status, routine, source, revision, revdate, &
                         'Opening file ' // PICKUP_NAME)
    
    write(pickup_unit, fmt=300) ktaust, ktauf
    
    close(pickup_unit)

300 format( 5(I3",")I3 )
  end subroutine write_pickup_file

  !------------------------------
  ! END PUBLIC ROUTINES
  !------------------------------

  !------------------------------
  ! BEGIN PRIVATE ROUTINES
  !------------------------------

  ! read_convert_namelist
  ! ---------------------
  ! Reads the conversion namelist
  !  PARAMETERS
  !   [none]
  subroutine read_convert_namelist()

    integer :: nml_unit

    character(len=*), parameter :: routine = 'read_convert_namelist'
    integer                     :: io_status

    nml_unit = get_unit()
    if (file_exist(nml_file)) then
       
       open(nml_unit, file=nml_file, status='old', iostat=io_status)
       call check_io_status(io_status, routine, source, revision, revdate, &
                            'Opening ' // nml_file)

       read(nml_unit, nml=convert, iostat=io_status)
       call check_io_status(io_status, routine, source, revision, revdate, &
                            'Reading ' // nml_file)

       close(nml_unit)
    else
       call error_handler(E_ERR, 'coamps_translate_mod', 'Convert name' // &
                          'list read failed - target file not found',      &
                          source, revision, revdate)
    end if
    write (*, nml=convert)
  end subroutine read_convert_namelist

  ! allocate_state
  ! --------------
  ! Allocate memory for the DART and COAMPS state vectors
  !  PARAMETERS
  !   [none]
  subroutine allocate_state()

    integer :: fieldsize, num_vars

    character(len=*), parameter :: routine = 'allocate_state'
    integer :: alloc_status

    ! Figure out how much we need
    call get_num_vars(num_vars)
    call get_grid_field_size(grid, fieldsize)
    
    ! DART
    allocate( dart_state(num_vars * fieldsize), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'dart_state')

    ! COAMPS
    allocate( coamps_state(num_vars * fieldsize), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'coamps_state')
  end subroutine allocate_state

  ! allocate_coamps_file_info
  ! -------------------------
  ! Allocate memory for the COAMPS restart file names and their 
  ! corresponding unit numbers - do this dynamically since the 
  ! number will (potentially) change from run to run based on the
  ! number of processors used.
  !  PARAMETERS
  !   [none]
  subroutine allocate_coamps_file_info()

    character(len=*), parameter :: routine = 'allocate_coamps_file_info'
    integer :: alloc_status
    
    if (npr0nam > 0) then
       ! We're using an I/O processor
       total_coamps_files = 1
    else
       ! We're doing parallel I/O - one file per process
       ! Hard-code a 1 since we're supporting a single nest
       total_coamps_files = ndxnam(1) * ndynam(1)
    end if

    allocate(coamps_file_names(total_coamps_files), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'coamps_file_names')
    
    allocate( coamps_file_units(total_coamps_files), stat=alloc_status )
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'coamps_file_units')
  end subroutine allocate_coamps_file_info

 ! coamps_process_all_fields
  ! ---------------------------------
  ! Reads/writes all fields in the DART restart vector to/from their
  ! proper places in the COAMPS restart file(s)
  !  PARAMETERS
  !   IN  write_coamps      True if we're writing to COAMPS files
  subroutine coamps_process_all_fields(write_coamps) 
    logical, intent(in) :: write_coamps

    integer :: cur_var        ! Which field in the state vector
    integer :: cur_point      ! Which point in the field
    integer :: fieldsize      ! Number of points on this grid
    integer :: num_vars       ! Fields in DART restart vector
    integer :: record         ! Field's record in COAMPS restart

    integer :: l_index        ! Bounds of DART restart vector subset
    integer :: u_index        ! for the field in question

    logical :: write_field    ! Allow us to not write a field even if
                              ! it's in the DART state vector
    logical :: use_singleio   ! True if using single process I/O

    ! For calculating loop limits and partial array indices
    call get_num_vars(num_vars)
    call get_grid_field_size(grid,fieldsize)

    if (npr0nam .gt. 0) then
       use_singleio = .true.
    else
       use_singleio = .false.
    end if

    do cur_var = 1, num_vars
       ! Figure out where this variable is both in the COAMPS
       ! restart file and the large state vector
       call get_full_record_num_by_index(cur_var,use_singleio,record)
       l_index = (cur_var - 1) * fieldsize + 1
       u_index = cur_var * fieldsize

       if (write_coamps) then

          ! Check if we should update this field - used if we need to
          ! have some numbers in the DART state vector for
          ! computations, but don't want to assimilate them
          ! (like mean profiles) 
          call get_update_flag_by_index(cur_var, write_field)

          if (write_field) then
             call coamps_restart_write_field(record, &
                   coamps_state(l_index:u_index)) 
          end if
       else
          call coamps_restart_read_field(record, &
                coamps_state(l_index:u_index)) 
       end if
    end do
  end subroutine coamps_process_all_fields

  ! coamps_restart_read_field
  ! -------------------------
  ! Reads a single field out of one or more COAMPS restart file(s)
  !  PARAMETERS
  !   IN  record            Field's record location in restart file
  !   OUT statevec_slice    Section of the state vector to read into
  subroutine coamps_restart_read_field(record, statevec_slice)
    integer, intent(in)                          :: record
    real(kind=C_REAL), dimension(:), intent(out) :: statevec_slice
    
    integer :: subgrid_index

    do subgrid_index = 1, size(coamps_file_units)
       call coamps_restart_read_single_file(subgrid_index,record,&
                                            statevec_slice)
    enddo
  end subroutine coamps_restart_read_field

  ! coamps_restart_write_field
  ! --------------------------
  ! Writes a single field into one or more COAMPS restart file(s)
  !  PARAMETERS
  !   IN  record            Field's record location in restart file
  !   IN  statevec_slice    Section of the state vector to write out
  subroutine coamps_restart_write_field(record, statevec_slice)
    integer, intent(in)                         :: record
    real(kind=C_REAL), dimension(:), intent(inout) :: statevec_slice

    integer :: subgrid_index

    do subgrid_index = 1,size(coamps_file_units)
       call coamps_restart_write_single_file(subgrid_index,record,&
                                             statevec_slice)
    enddo
  end subroutine coamps_restart_write_field

  ! coamps_restart_read_single_file
  ! -------------------------------
  ! Reads either a complete (in the case of single-processor I/O)
  ! or partial (for multi-proc I/O) field out of one restart file
  !  PARAMETERS
  !   IN  subgrid_index     decomposed domain index that we'll
  !                         read from
  !   IN  record            Field's record location in restart file
  !   OUT statevec_slice    Section of the state vector to read into
  subroutine coamps_restart_read_single_file(subgrid_index, record,&
                                             statevec_slice)
    integer, intent(in)                          :: subgrid_index
    integer, intent(in)                          :: record
    real(kind=C_REAL), dimension(:), intent(out) :: statevec_slice

    logical :: READING_COAMPS = .false.

    call coamps_restart_do_single_file(subgrid_index,record,        &
                                       READING_COAMPS,statevec_slice) 
  end subroutine coamps_restart_read_single_file

  ! coamps_restart_write_single_file
  ! --------------------------------
  ! Writes either a complete (in the case of single-processor I/O)
  ! or partial (for multi-proc I/O) field to one restart file
  ! N.B. The true intent of statevec_slice is intent(in), but since
  ! coamps_restart_do_single_file either reads or writes it can
  ! cause an error
  !  PARAMETERS
  !   IN  subgrid_index     decomposed domain index that we'll
  !                         write to
  !   IN  record            Field's record location in restart file
  !   IN  statevec_slice    Section of the state vector to write from
  subroutine coamps_restart_write_single_file(subgrid_index, record,&
                                              statevec_slice)
    integer, intent(in)                            :: subgrid_index
    integer, intent(in)                            :: record
    real(kind=C_REAL), dimension(:), intent(inout) :: statevec_slice

    logical :: WRITING_COAMPS = .true.

    call coamps_restart_do_single_file(subgrid_index, record, &
                                       WRITING_COAMPS,        &
                                       statevec_slice         )
  end subroutine coamps_restart_write_single_file


  ! coamps_restart_do_single_file
  ! ------------------------------
  ! Computes the position of field data within a COAMPS restart file,
  ! then reads or writes the field section to the file.
  !  PARAMETERS
  !   IN  subgrid_index     decomposed domain to read or write
  !   IN  field_record      where the field is in the COAMPS restart
  !                         file
  !   IN  writing_coamps    True if we're writing to the COAMPS 
  !                         restart file
  ! INOUT statevec_slice    The section of the state vector under
  !                         consideration
  subroutine coamps_restart_do_single_file(subgrid_index, &
              field_record, writing_coamps, statevec_slice) 
    integer, intent(in)                            :: subgrid_index
    integer, intent(in)                            :: field_record
    logical, intent(in)                            :: writing_coamps
    real(kind=C_REAL), dimension(:), intent(inout) :: statevec_slice

    ! Calculate/store the locations of the data to read and write
    ! in both the state vector slice and the COAMPS restart file
    ! (which is measured as an offset from a base field record).
    integer :: iextent
    integer :: jextent
    integer :: ijarea
    integer :: subfield_record
    integer :: slice_pos
    integer :: restart_offset
    integer :: restart_pos

    ! Index locations in the local subdomain - includes current
    ! locations (ii/jj) and the i/j boundaries based (possibly)
    ! on the domain decomposition 
    integer :: ii,jj     
    integer :: ilbound  
    integer :: iubound 
    integer :: jlbound
    integer :: jubound
    logical :: ij_within_grid

    ! I/O
    integer :: restart_unit

    ! Error handling
    character(len=*), parameter :: routine = 'coamps_restart_do_single_file'
    integer :: io_status

    restart_unit = coamps_file_units(subgrid_index) 

    ! Start trying to read/write at position field_record in the 
    ! file and at position 1 in the data array - doing it this way 
    ! to avoid any odd (jj-1)*isize + ii computations
    restart_offset  = 0

    ! The "subfield" boundaries for the single-processor I/O case
    ! is just the i/j limits for the entire field - multi-processor
    ! needs the domain decomposition data
    if (npr0nam .gt. 0) then
       ilbound = 1
       jlbound = 1
       call get_grid_dims(grid, iubound, jubound)
    else
       call get_grid_iminf(grid, subgrid_index, ilbound)
       call get_grid_imaxf(grid, subgrid_index, iubound)

       call get_grid_jminf(grid, subgrid_index, jlbound)
       call get_grid_jmaxf(grid, subgrid_index, jubound)
    end if

    ! Need to take into account that the size of the subfield may
    ! vary from file to file based on the domain decomposition
    iextent = iubound - ilbound + 1
    jextent = jubound - jlbound + 1
    ijarea  = iextent * jextent
    subfield_record  = ijarea * (field_record-1) + 1

    do jj = jlbound, jubound
       do ii =ilbound, iubound

          ! Don't read or write any points if they are in the halo
          ! region around the domain
          call check_ij_within_grid(grid, ii, jj, ij_within_grid)
          if ( .not. ij_within_grid) then
             restart_offset = restart_offset + 1
             cycle
          end if

          ! We've established that we're in the region, read or write
          ! this particular value in the DART state vector
          call grid_ij_to_vector_index(grid, ii, jj, slice_pos)
          restart_pos = subfield_record + restart_offset

          if (writing_coamps) then
             write(restart_unit, rec=restart_pos, iostat=io_status)&
                  & statevec_slice(slice_pos)
          else
             read(restart_unit, rec=restart_pos, iostat=io_status)&
                  & statevec_slice(slice_pos)
          end if
          call check_io_status(io_status, routine, source, revision, &
                               revdate, 'COAMPS restart file')

          restart_offset  = restart_offset + 1
       enddo
    enddo
  end subroutine coamps_restart_do_single_file

  ! convert_state
  ! -------------
  ! Converts a state vector to or from the DART numeric format to or
  ! from the COAMPS format.  The state vectors are stored as module
  ! variables so no need to pass them in as parameters 
  ! PARAMETERS
  !  IN  to_dart            True if we're converting to DART format
  subroutine convert_state(to_dart)
    logical, intent(in)                             :: to_dart
    
    ! If COAMPS used an I/O processor, then the restart file written
    ! out is big-endian.  This means that either the vector *read*
    ! from COAMPS or the vector to be written *to* COAMPS needs to
    ! be big-endian as well.
    logical :: used_io_proc

    if (npr0nam .gt. 0) then
      used_io_proc = .true.
    else
      used_io_proc = .false.
    end if

    if (to_dart) then
      ! big_endian -> little-endian
      if (used_io_proc) call change_coamps_endian()
      dart_state = real(coamps_state,kind=r8)
    else
      coamps_state = real(dart_state,kind=C_REAL)
      ! little-endian -> big_endian
      if (used_io_proc) call change_coamps_endian()
    endif
  end subroutine convert_state
  

  ! fix_negative_values
  ! -------------------
  ! Sets values in the state vector that are less than zero to zero
  ! if they have been defined to be always positive. - this includes
  ! things like mixing ratios, etc.
  !  PARAMETERS
  !   [none]
  subroutine fix_negative_values()

    integer :: num_vars
    integer :: fieldsize

    integer :: cur_var

    logical :: var_is_nonnegative

    integer :: var_lower_index
    integer :: var_upper_index

    real(kind=r8), dimension(:), pointer :: dart_substate


    call get_num_vars(num_vars)
    call get_grid_field_size(grid, fieldsize)

    do cur_var = 1, num_vars
        call get_posdef_flag_by_index(cur_var, var_is_nonnegative) 

        var_lower_index = (cur_var - 1) * fieldsize + 1
        var_upper_index = (cur_var    ) * fieldsize

        if (var_is_nonnegative) then
            dart_substate => dart_state(var_lower_index:var_upper_index)
            where (dart_substate < 0) dart_substate = 0     
        end if
    end do
  end subroutine fix_negative_values

  ! change_coamps_endian
  ! --------------------
  ! Switch the byte ordering of the COAMPS state vector - note
  ! that fix_for_platform just returns the state vector untouched
  ! if we happen to be on a big-endian platform already
  !  PARAMETERS
  !   [none]
  subroutine change_coamps_endian()
    call fix_for_platform(coamps_state, size(coamps_state), C_REAL)
  end subroutine change_coamps_endian

  ! sd_to_hms
  ! ---------
  ! Convert a time given in the DART (second, day) format into the
  ! COAMPS (hour, minute, second) format
   !  PARAMETERS
  !   IN  sd_day            The "day" bit of (second, day)
  !   IN  sd_sec            The "second" bit of (second, day)
  !   OUT hms_hour          The "hour" bit of (hour, minute, second)
  !   OUT hms_min           The "minute" bit of (hour,minute,second)
  !   OUT hms_sec           The "second" bit of (hour,minute,second)
 subroutine sd_to_hms(sd_day, sd_sec, hms_hour, hms_min, hms_sec)
    integer, intent(in)  :: sd_day
    integer, intent(in)  :: sd_sec
    integer, intent(out) :: hms_hour
    integer, intent(out) :: hms_min
    integer, intent(out) :: hms_sec

    integer :: residual_secs
    
    residual_secs = (sd_day * DAY_TO_SEC) + sd_sec

    ! Note integer division
    hms_hour      = residual_secs / HR_TO_SEC
    residual_secs = residual_secs - (hms_hour * HR_TO_SEC)
    hms_min       = residual_secs / MIN_TO_SEC
    residual_secs = residual_secs - (hms_min * MIN_TO_SEC)
    hms_sec       = residual_secs
  end subroutine sd_to_hms

  ! hms_to_sd
  ! ---------
  ! Convert a time given in the COAMPS (hour, minute, second) format
  ! into the DART (second, day) format
  !  PARAMETERS
  !   IN  hms_hour          The "hour" bit of (hour, minute, second)
  !   IN  hms_min           The "minute" bit of (hour,minute,second)
  !   IN  hms_sec           The "second" bit of (hour,minute,second)
  !   OUT sd_day            The "day" bit of (second, day)
  !   OUT sd_sec            The "second" bit of (second, day)
  subroutine hms_to_sd(hms_hour, hms_min, hms_sec, sd_day, sd_sec)
    integer, intent(in)  :: hms_hour
    integer, intent(in)  :: hms_min
    integer, intent(in)  :: hms_sec
    integer, intent(out) :: sd_day
    integer, intent(out) :: sd_sec

    integer :: residual_secs

    residual_secs = (hms_hour * HR_TO_SEC)  + &
                    (hms_min  * MIN_TO_SEC) + &
                     hms_sec

    ! Note integer division
    sd_day        = residual_secs / (DAY_TO_SEC)
    residual_secs = residual_secs - (sd_day * DAY_TO_SEC)
    sd_sec        = residual_secs
  end subroutine hms_to_sd

  ! print_dart_diagnostics
  ! ----------------------
  ! Print out information about the maximum and minimum values
  ! contained in the state vector and which variable(s) they 
  ! correspond to.
  !  PARAMETERS
  !   [none]
  subroutine print_dart_diagnostics()
    real(kind=r8)               :: max_value, min_value
    integer, dimension(1)       :: max_location, min_location
    integer                     :: level
    character(len=5)            :: var_name

    ! Get the values
    max_value = maxval(dart_state)
    min_value = minval(dart_state)

    ! Get where they are
    max_location = maxloc(dart_state)
    min_location = minloc(dart_state)

    print *, "Values in DART State Vector"
    call get_var_info_by_abs_index(max_location(1), var_name, level)
    write (*,*) "Maximum is ", max_value, "in variable ", var_name,&
                "on level ", level
    call get_var_info_by_abs_index(min_location(1), var_name, level)
    write (*,*) "Minimum is ", min_value, "in variable ", var_name,&
                "on level ", level
    write (*,*) "First entry in state vector is", dart_state(1)
  end subroutine print_dart_diagnostics

  !------------------------------
  ! END PRIVATE ROUTINES
  !------------------------------

end module coamps_translate_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
