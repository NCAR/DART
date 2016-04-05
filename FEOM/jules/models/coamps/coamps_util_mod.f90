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
  
  use types_mod,     only : r8
  use utilities_mod, only : E_ERR,         &
                            error_handler, &
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

  ! Variable type
  public :: C_REAL

  !------------------------------
  ! END PUBLIC INTERFACE
  !------------------------------

  !------------------------------
  ! BEGIN EXTERNAL INTERFACES
  !------------------------------
  ! (This code calls a C routine but leave the interface as
  !  implicit to allow calling for various sizes of REAL)
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
  integer, parameter :: C_REAL = 4

  ! Error handling
  integer, parameter :: NUM_ERROR_TYPES       = 3

  integer, parameter :: ERROR_TYPE_ALLOCATION   = 1
  integer, parameter :: ERROR_TYPE_DEALLOCATION = 2
  integer, parameter :: ERROR_TYPE_IO           = 3

  ! Constants for pretty-printing
  character, parameter :: DIVCHAR     = '-'
  integer,   parameter :: LABEL_WIDTH = 45
  integer,   parameter :: LABEL_LEAD  = 2

  ! Set this to true if we're on a little-endian platform
  logical, parameter :: LITTLE_ENDIAN_PLATFORM = .true.

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

  character(len=28), dimension(NUM_ERROR_TYPES) :: error_msgs

  !------------------------------
  ! END MODULE VARIABLES
  !------------------------------

contains

  !------------------------------
  ! BEGIN PUBLIC ROUTINES
  !------------------------------
  
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

  ! fix_for_platform
  ! ----------------
  ! Performs any platform-specific fixes necessary on an array of
  ! C_REAL reals.  Currently these are:
  !  - Little-endian platform:  COAMPS writes big-endian files, so
  !                             need to byteswap 
  !  PARAMETERS
  ! INOUT array             array containing values byteswap
  !   IN  array_size        length of the array
  !   IN  element_size      length (in bytes) of each array entry 
  subroutine fix_for_platform(array, array_size, element_size)
    real(kind=C_REAL), dimension(:), intent(inout) :: array
    integer,                         intent(in)    :: array_size
    integer,                         intent(in)    :: element_size

    if (LITTLE_ENDIAN_PLATFORM) then
      call c_byteswap_array(array, array_size, element_size)
    end if
  end subroutine fix_for_platform

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

  !------------------------------
  ! END PRIVATE ROUTINES
  !------------------------------

end module coamps_util_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
