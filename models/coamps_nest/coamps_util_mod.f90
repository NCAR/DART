! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

module coamps_util_mod
!------------------------------
! MODULE:       coamps_util_mod
! AUTHOR:       T. R. Whitcomb
!               Naval Research Laboratory
! DART VERSION: Jamaica
!
! Module with various COAMPS utility routines to handle things
! like error checking
!------------------------------ 
  
  use types_mod,     only : r4, r8
  use utilities_mod, only : E_ERR,         &
                            E_MSG,         &
                            error_handler, &
                            do_output,     &
                            get_unit,      &
                            timestamp,     &
                            register_module

  implicit none

  private

  !------------------------------
  ! BEGIN PUBLIC INTERFACE
  !------------------------------

  ! Error checking
  public :: check_io_status
  public :: check_alloc_status
  public :: check_dealloc_status

  ! Pretty-print
  public :: print_label_name
  public :: print_2d_real8_array

  ! Icky platform-specific details
  public :: fix_for_platform

  ! Working with flat files
  public :: generate_flat_file_name
  public :: read_flat_file
  public :: write_flat_file
  public :: read_datahd_file
  public :: write_datahd_file

  public :: dump_data_file

  ! Variable type
  public :: C_REAL

  ! Change string case
  public :: lowercase
  public :: uppercase

  public :: set_debug_level
  public :: trace_message
  public :: timestamp_message

  public :: DATAHD_LEN
  public :: DATAHD_NUM_NESTS

  !------------------------------
  ! END PUBLIC INTERFACE
  !------------------------------

  !------------------------------
  ! BEGIN EXTERNAL INTERFACES
  !------------------------------
  interface fix_for_platform
    module procedure fix_for_platform4, fix_for_platform8
  end interface fix_for_platform
  !------------------------------
  ! END EXTERNAL INTERFACES
  !------------------------------

  !------------------------------
  ! BEGIN TYPES AND CONSTANTS 
  !------------------------------

  ! Size (in bytes) of the real number written/read by COAMPS
  ! Changing this obviously requires a recompile of the DART
  ! interface, but it also requires a recompile of COAMPS so
  ! I'm not too concerned.
  integer, parameter :: C_REAL  = r8

  ! Error handling
  integer, parameter :: NUM_ERROR_TYPES         = 3

  integer, parameter :: ERROR_TYPE_ALLOCATION   = 1
  integer, parameter :: ERROR_TYPE_DEALLOCATION = 2
  integer, parameter :: ERROR_TYPE_IO           = 3

  ! Constants for pretty-printing
  character, parameter :: DIVCHAR     = '-'
  integer,   parameter :: LABEL_WIDTH = 45
  integer,   parameter :: LABEL_LEAD  = 2

  ! Set this to true if we're on a little-endian platform
  logical, parameter :: LITTLE_ENDIAN_PLATFORM = .true.

  integer, parameter :: DATAHD_LEN       = 2000
  integer, parameter :: DATAHD_NUM_NESTS = 11

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

  logical :: module_initialized = .false.
  integer :: debug_level = 0

  character(len=28), dimension(NUM_ERROR_TYPES) :: error_msgs

  !------------------------------
  ! END MODULE VARIABLES
  !------------------------------

contains

  !------------------------------
  ! BEGIN PUBLIC ROUTINES
  !------------------------------
  subroutine set_debug_level(debug)
    integer, intent(in) :: debug

    if (.not. module_initialized) call initialize_module()
    debug_level = debug
  end subroutine set_debug_level

  subroutine timestamp_message(msg, sync)

  character(len=*), intent(in) :: msg
  logical, intent(in), optional :: sync

  ! Write current time and message to stdout and log file.
  ! if sync is present and true, sync mpi jobs before printing time.

  if (debug_level <= 0) return

  !if (present(sync)) then
  !  if (sync) call task_sync()
  !endif

  if (do_output()) call timestamp(' '//trim(msg), pos='brief')  ! was debug

  end subroutine timestamp_message


  subroutine trace_message(msg, label, threshold)

  character(len=*), intent(in)           :: msg
  character(len=*), intent(in), optional :: label
  integer,          intent(in), optional :: threshold

  ! Write message to stdout and log file.
  integer :: t

  t = 0
  if (present(threshold)) t = threshold

  !if (debug_level <= t) return

  !if (.not. do_output()) return

  if (present(label)) then
     call error_handler(E_MSG,trim(label),trim(msg))
  else
     call error_handler(E_MSG,'coamps trace:',trim(msg))
  endif

  end subroutine trace_message

  ! lowercase
  ! ---------------
  ! Wraper to change to case of a string to lower case
  !  PARAMETERS
  !   IN  str  Character string to change case of
  function lowercase(str) result(lstr)
    character(len=*), intent(in)        :: str
    character(len=len(str))             :: lstr
    logical                             :: TO_LOWER = .true.
    lstr = str
    call change_case(lstr,TO_LOWER)
  end function lowercase
  
  ! uppercase
  ! ---------------
  ! Wraper to change to case of a string to upper case
  ! Sets a string to upper case
  !  PARAMETERS
  !   IN  str  Character string to change case of
  function uppercase(str) result(ustr)
    character(len=*), intent(in)       :: str
    character(len=len(str))            :: ustr
    logical                            :: TO_LOWER = .false.
    ustr = str
    call change_case(ustr,TO_LOWER)
  end function uppercase
                           
  ! check_io_status
  ! ---------------
  ! Wrapper for check status function for I/O errors
  !  PARAMETERS
  !   IN  io_status         Return value from I/O function
  !   IN  routine           Subroutine that called I/O function
  !   IN  source            Source file that called I/O function
  !   IN  revision          Revision of source file
  !   IN  revdate           Revision date of source file
  !   IN  context           OPTIONAL error message context
  subroutine check_io_status(io_status, routine, source, revision, &
                             revdate, context)
    integer, intent(in)                    :: io_status
    character(len=*), intent(in)           :: routine
    character(len=*), intent(in)           :: source
    character(len=*), intent(in)           :: revision
    character(len=*), intent(in)           :: revdate
    character(len=*), intent(in), optional :: context

    character(len=128) :: error_context
    
    if (.not. module_initialized) call initialize_module()

    if (present(context)) then
      error_context = context
    else
      error_context = ""
    end if

    call check_error_status(ERROR_TYPE_IO, io_status, routine,  &
                            error_context, source, revision,    &
                            revdate)
  end subroutine check_io_status

  ! check_alloc_status
  ! ------------------
  ! Wrapper for check status function for memory allocation errors 
  !  PARAMETERS
  !   IN  alloc_status      Return value from allocate()
  !   IN  routine           Subroutine that called allocate() 
  !   IN  source            Source file that called allocate()
  !   IN  revision          Revision of source file
  !   IN  revdate           Revision date of source file
  !   IN  context           OPTIONAL error message context
  subroutine check_alloc_status(alloc_status, routine, source, &
                                revision, revdate, context)
    integer, intent(in)                    :: alloc_status
    character(len=*), intent(in)           :: routine
    character(len=*), intent(in)           :: source
    character(len=*), intent(in)           :: revision
    character(len=*), intent(in)           :: revdate
    character(len=*), intent(in), optional :: context

    character(len=128) :: error_context
    
    if (.not. module_initialized) call initialize_module()

    if (present(context)) then
      error_context = context
    else
      error_context = ""
    end if

    call check_error_status(ERROR_TYPE_ALLOCATION, alloc_status,&
                            routine, error_context, source,     &
                            revision, revdate)
  end subroutine check_alloc_status

  ! check_dealloc_status
  ! ------------------
  ! Wrapper for check status function for memory deallocation errors 
  !  PARAMETERS
  !   IN  dealloc_status    Return value from deallocate()
  !   IN  routine           Subroutine that called deallocate() 
  !   IN  source            Source file that called deallocate()
  !   IN  revision          Revision of source file
  !   IN  revdate           Revision date of source file
  !   IN  context           OPTIONAL error message context
  subroutine check_dealloc_status(dealloc_status, routine, source, &
                                  revision, revdate, context)
    integer, intent(in)                    :: dealloc_status
    character(len=*), intent(in)           :: routine
    character(len=*), intent(in)           :: source
    character(len=*), intent(in)           :: revision
    character(len=*), intent(in)           :: revdate
    character(len=*), intent(in), optional :: context

    character(len=128) :: error_context
    
    if (.not. module_initialized) call initialize_module()

    if (present(context)) then
      error_context = context
    else
      error_context = ""
    end if

    call check_error_status(ERROR_TYPE_DEALLOCATION, dealloc_status,&
                            routine, error_context, source,         &
                            revision, revdate)
  end subroutine check_dealloc_status

  ! print_array_name
  ! ----------------
  ! Given a name, prints a pretty label 
  !  PARAMETERS
  !   IN  label_name        Name to print
  subroutine print_label_name(label_name)
    character(len=*), intent(in) :: label_name

    integer :: name_len
    integer :: left_to_print

    name_len = len(label_name)
    left_to_print = LABEL_WIDTH - LABEL_LEAD - name_len

    write (*,*) repeat(DIVCHAR,LABEL_LEAD) // label_name // &
                repeat(DIVCHAR,left_to_print)
  end subroutine print_label_name
  
  ! print_2d_real8_array
  ! --------------------
  ! Prints a two-dimensional array of type real(r8) as an array
  ! instead of the usual write-everything-out-at-once behavior
  !  PARAMETERS
  !   IN  array             The 2D array to print
  subroutine print_2d_real8_array(array)
    real(kind=r8), dimension(:,:) :: array

    integer :: cur_row
    integer :: cur_col
    integer :: num_rows
    integer :: num_cols

    character(len=120) :: format_str

    num_rows = size(array, dim=1)
    num_cols = size(array, dim=2)

    ! Assemble format string
    write (format_str,'(A1 I2 A5 A1)') '(', num_cols, 'D15.5', ')'

    do cur_row = 1, num_rows
      write (*,trim(format_str)) array(cur_row,:)
    end do
  end subroutine print_2d_real8_array

  ! fix_for_platform4
  ! ----------------
  ! Performs any platform-specific fixes necessary on an array of r4's
  ! Currently these are:
  !  - Little-endian platform:  COAMPS writes big-endian files, so
  !                             need to byteswap 
  !  PARAMETERS
  ! INOUT array             array containing values byteswap
  !   IN  array_size        length of the array
  !   IN  element_size      length (in bytes) of each array entry 
  subroutine fix_for_platform4(array, array_size)
    real(kind=r4),      dimension(:), intent(inout) :: array
    integer,                          intent(in)    :: array_size

    if (LITTLE_ENDIAN_PLATFORM) then
      call c_byteswap_array(array, array_size, r4)
    end if
  end subroutine fix_for_platform4

  ! fix_for_platform8
  ! ----------------
  ! Performs any platform-specific fixes necessary on an array of r8's
  ! Currently these are:
  !  - Little-endian platform:  COAMPS writes big-endian files, so
  !                             need to byteswap 
  !  PARAMETERS
  ! INOUT array             array containing values byteswap
  !   IN  array_size        length of the array
  !   IN  element_size      length (in bytes) of each array entry 
  subroutine fix_for_platform8(array, array_size)
    real(kind=r8),      dimension(:), intent(inout) :: array
    integer,                          intent(in)    :: array_size

    if (LITTLE_ENDIAN_PLATFORM) then
      call c_byteswap_array(array, array_size, r8)
    end if
  end subroutine fix_for_platform8

  ! generate_flat_file_name
  ! -----------------------
  ! Given field, level, and grid information, generate the properly
  ! formatted 64-character COAMPS flat file name.  Note that this
  ! does *not* generate any path information - it only returns the
  ! file name.
  !  PARAMETERS
  !   IN  var_name          the field the file contains
  !   IN  level_type        vertical level type (height/pressure/etc)
  !   IN  level1            lowest vertical level in the file
  !   IN  level2            highest vertical level in the file
  !                         (for files for a single level, level1 is 
  !                          that level and level2 is left to 0)
  !   IN  gridnum           nest number (only 1 supported for now)
  !   IN  aoflag            field type: (a)tmosphere or (o)cean
  !   IN  xpts              number of points in the x direction
  !   IN  ypts              number of points in the y direction
  !   IN  dtg               base date-time group
  !   IN  tau_hh            forecast lead time - hour component
  !   IN  tau_mm            forecast lead time - minute component
  !   IN  tau_ss            forecast lead time - second component
  !   IN  field_type        type of field (e.g. fcstfld, infofld)
  !   OUT file_name         COAMPS flat file name
  subroutine generate_flat_file_name(var_name, level_type, level1, level2,     &
                                     gridnum, aoflag, xpts, ypts, dtg, tau_hh, &
                                     tau_mm, tau_ss, field_type, file_name)
    character(len=6),  intent(in)  :: var_name
    character(len=3),  intent(in)  :: level_type
    integer,           intent(in)  :: level1
    integer,           intent(in)  :: level2
    integer,           intent(in)  :: gridnum
    character(len=1),  intent(in)  :: aoflag
    integer,           intent(in)  :: xpts
    integer,           intent(in)  :: ypts
    character(len=10), intent(in)  :: dtg
    integer,           intent(in)  :: tau_hh
    integer,           intent(in)  :: tau_mm
    integer,           intent(in)  :: tau_ss
    character(len=7),  intent(in)  :: field_type
    character(len=64), intent(out) :: file_name

    write(file_name, 100) lowercase(var_name), level_type, level1, level2, &
         & gridnum, aoflag, xpts, ypts, dtg, tau_hh, tau_mm,    &
         & tau_ss, field_type

100 format(A6,'_',A3,'_',I6.6,'_',I6.6,'_',I1,A1,I4.4,'x',I4.4,'_',A10,  &
           '_',I4.4,I2.2,I2.2,'_',A7)
  end subroutine generate_flat_file_name

  ! write_flat_file
  ! ----------------
  ! Given the unit number of an *open* COAMPS flat 
  ! file, read the file into an array. 
  !  PARAMETERS
  !   IN  flat_unit      unit number of an open flat file
  !   OUT flat_array     coamps_grid structure to be filled
  subroutine write_flat_file(flat_unit, flat_array)
    integer,                          intent(in)  :: flat_unit
    real(kind=r8),      dimension(:), intent(in)  :: flat_array

    character(len=*), parameter                   :: routine = 'write_flat_file'

    real(kind=r4),      dimension(:), allocatable :: flat_array_tmp
    integer                                       :: io_status
    integer                                       :: alloc_status
    integer                                       :: dealloc_status
    integer                                       :: field_size
    integer                                       :: i

    field_size=size(flat_array)
    allocate(flat_array_tmp(field_size), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'flat_array_tmp')

    ! COAMPS flat files are real(kind=r4)
    flat_array_tmp(:)=real(flat_array(:) , kind=r4)
    call fix_for_platform(flat_array_tmp, field_size)

    write(unit=flat_unit, rec=1, iostat=io_status) flat_array_tmp
    call check_io_status(io_status, routine, source, revision, &
                         revdate, 'writing flat file')

    deallocate(flat_array_tmp, stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'flat_array_tmp')
  end subroutine write_flat_file

  ! read_flat_file
  ! ----------------
  ! Given the unit number of an *open* COAMPS flat 
  ! file, read the file into an array. 
  !  PARAMETERS
  !   IN  flat_unit      unit number of an open flat file
  !   OUT flat_array     coamps_grid structure to be filled
  subroutine read_flat_file(flat_unit, flat_array)
    integer, intent(in)                            :: flat_unit
    real(kind=r8),     dimension(:), intent(inout) :: flat_array

    real(kind=r4),       dimension(:), allocatable :: flat_array_tmp
    character(len=*), parameter :: routine = 'read_flat_file'
    integer :: io_status, alloc_status, dealloc_status
    integer :: field_size

    field_size=size(flat_array)
    allocate(flat_array_tmp(field_size), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'flat_array_tmp')

    ! Read in the data - COAMPS writes flat files as r4's
    read(unit=flat_unit, rec=1, iostat=io_status) flat_array_tmp
    call check_io_status(io_status, routine, source, revision, &
                         revdate, 'Reading flat file')
    call fix_for_platform(flat_array_tmp, field_size)
    flat_array(:)=real(flat_array_tmp(:) , kind=r8)

    deallocate(flat_array_tmp, stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'flat_array_tmp')
  end subroutine read_flat_file

  ! read_datahd_file
  ! ----------------
  ! Given the dtg of a COAMPS domain information file, read the domain 
  ! parameters from the file into an array.
  !  PARAMETERS
  !   IN  dtg               COAMPS date-time-group for filename generation
  !   OUT datahd            array to be filled
  subroutine read_datahd_file(dtg, datahd)
    character(len=10),            intent(in) :: dtg
    real(kind=r8), dimension(:), intent(out) :: datahd

    character(len=64)                    :: datahd_filename
    integer                              :: datahd_unit

    ! Error checking
    logical :: is_opened
    character(len=*), parameter :: routine = 'read_datahd_file'
    integer :: io_status, alloc_status

    integer :: ii

    call generate_datahd_filename(dtg, datahd_filename)
    datahd_unit = get_unit()
    open(unit=datahd_unit, file=datahd_filename, status='old', &
         access='sequential', action='read', form='formatted', &
         iostat=io_status)
    call check_io_status(io_status, routine, source, revision, &
                         revdate, 'Opening datahd file')

    read(unit=datahd_unit, fmt='(5e13.6)', iostat=io_status)   &
        (datahd(ii), ii=1,size(datahd))
    call check_io_status(io_status, routine, source, revision, &
                         revdate, 'Reading datahd file')

    close(datahd_unit)
  end subroutine read_datahd_file

  ! write_datahd_file
  ! ----------------
  ! Given the dtg of a COAMPS domain information file, read the domain 
  ! parameters from the file into an array.
  !  PARAMETERS
  !   IN  dtg               COAMPS date-time-group for filename generation
  !   OUT datahd            array to be filled
  subroutine write_datahd_file(dtg, datahd)
    character(len=10),            intent(in) :: dtg
    real(kind=r8), dimension(:),  intent(in) :: datahd

    character(len=64)                    :: datahd_filename
    integer                              :: datahd_unit

    ! Error checking
    logical :: is_opened
    character(len=*), parameter :: routine = 'read_datahd_file'
    integer :: io_status, alloc_status

    integer :: ii

    call generate_datahd_filename(dtg, datahd_filename)
    datahd_unit = get_unit()
    open(unit=datahd_unit, file=datahd_filename, status='replace', &
         access='sequential', action='write', form='formatted', &
         iostat=io_status)
    call check_io_status(io_status, routine, source, revision, &
                         revdate, 'Opening datahd file')

    write(unit=datahd_unit, fmt='(5e13.6)', iostat=io_status)   &
        (datahd(ii), ii=1,size(datahd))
    call check_io_status(io_status, routine, source, revision, &
                         revdate, 'Reading datahd file')

    close(datahd_unit)
  end subroutine write_datahd_file

  subroutine dump_data_file(fname, dat_array)
    character(len=*), intent(in) :: fname
    real(kind=r8),    intent(in) :: dat_array(:)

    integer                     :: r4_length
    integer                     :: rec_length
    integer                     :: flat_unit

    character(len=*), parameter :: routine = 'dump_data_file'
    integer                     :: io_status

    inquire(IOLENGTH=r4_length) 0_r4

    flat_unit = get_unit()

    rec_length = r4_length * size(dat_array)

     open( unit=flat_unit,                                     &
          file='./'//trim(adjustl(fname)),                     &
          status='unknown', access='direct', action='write',   &
          form='unformatted', recl=rec_length, iostat=io_status)

    call check_io_status(io_status, routine, source, revision, &
                         revdate, 'Opening ./' // fname )  

    call write_flat_file(flat_unit, dat_array)

    close(flat_unit)

  end subroutine dump_data_file

  !------------------------------
  ! END PUBLIC ROUTINES
  !------------------------------

  !------------------------------
  ! BEGIN PRIVATE ROUTINES
  !------------------------------

  ! initialize_module
  ! -----------------
  ! Initialize model by populating error message structure
  !  PARAMETERS
  !   [none]
  subroutine initialize_module()
    ! Only need to do this once
    module_initialized = .true.

    call register_module(source, revision, revdate)

    ! Populate error message types
    error_msgs(ERROR_TYPE_IO)           = "I/O Error::"
    error_msgs(ERROR_TYPE_ALLOCATION)   = "Allocation Error::"
    error_msgs(ERROR_TYPE_DEALLOCATION) = "Deallocation Error::"
  end subroutine initialize_module

  ! check_error_status
  ! ------------------
  ! Checks the error status returned from various function calls
  ! Prints an error message if action is needed
  !  PARAMETERS
  !   IN  error_type        type of error (I/O, allocate, etc.)
  !   IN  error_status      Return value from function
  !   IN  routine           Subroutine that called function 
  !   IN  context           error message context
  !   IN  source            Source file that called function
  !   IN  revision          Revision of source file
  !   IN  revdate           Revision date of source file
  subroutine check_error_status(error_type, error_status, routine,&
                                context, source, revision, revdate)
    integer, intent(in)                    :: error_type
    integer, intent(in)                    :: error_status
    character(len=*), intent(in)           :: routine
    character(len=*), intent(in)           :: context
    character(len=*), intent(in)           :: source
    character(len=*), intent(in)           :: revision
    character(len=*), intent(in)           :: revdate

    character(len=128) :: error_message

    if (.not. module_initialized) call initialize_module

    ! Only need to do something if an error occured
    if (error_status .ne. 0) then
      error_message = trim(error_msgs(error_type)) // trim(context)
      call error_handler(E_ERR, routine, error_message, source,    &
                         revision, revdate)
    end if
  end subroutine check_error_status

  ! change_case
  ! ---------------
  ! Changes the case of a string.
  !  PARAMETERS
  !   IN  string     Character string to change case of
  !   IN  TO_LOWER   Change to lowercase or uppercase.
  subroutine change_case(string,TO_LOWER)
    character(len=*), intent(inout)  :: string
    logical, intent(in)              :: TO_LOWER

    integer                          :: idiff, i1, i
    integer                          :: ia, iz

    if(TO_LOWER) then
      idiff = ichar('a') - ichar('A')
      ia = ichar('A') ; iz = ichar('Z')
    else
      idiff = ichar('A') - ichar('a')
      ia = ichar('a') ; iz = ichar('z')
    end if

    do i=1,len(string)
      i1 = ichar(string(i:i))
      if( i1.ge.ia .and. i1.le.iz) then
        string(i:i) = char(ichar(string(i:i)) + idiff)
      end if
    end do
  end subroutine change_case

  ! generate_datahd_filename
  ! ------------------------
  ! Generates the COAMPS domain information file name for the first
  ! nest at a given date-time group.  This does *not* generate any
  ! path information - it only returns the file name.
  !  PARAMETERS
  !   IN  dtg               base date-time group for model run
  !   OUT datahd_filename   name of domain information file
  subroutine generate_datahd_filename(dtg, datahd_filename)
    character(len=10), intent(in)  :: dtg
    character(len=64), intent(out) :: datahd_filename

    ! The format of the datahd filename is fixed except for the date
    ! -time group: mimics a 2000x1 flat file with no level
    ! information or tau information
    call generate_flat_file_name( var_name   = 'datahd',      &
                                  level_type = 'sfc',         &
                                  level1     = 0,             &
                                  level2     = 0,             &
                                  gridnum    = 1,             &
                                  aoflag     = 'a',           &
                                  xpts       = 2000,          &
                                  ypts       = 1,             &
                                  dtg        = dtg,           &
                                  tau_hh     = 0,             &
                                  tau_mm     = 0,             &
                                  tau_ss     = 0,             &
                                  field_type = 'infofld',     &
                                  file_name  = datahd_filename )
  end subroutine generate_datahd_filename

  !------------------------------
  ! END PRIVATE ROUTINES
  !------------------------------

end module coamps_util_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
