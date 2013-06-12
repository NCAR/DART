! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

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
module coamps_translate_mod

  use coamps_domain_mod,    only : coamps_domain, decompose_domain,           &
                                   get_domain_nest, get_nest_count,           &
                                   get_num_subdomains, get_domain_num_levels, &
                                   get_domain_msigma, get_domain_wsigma,      &
                                   initialize_domain

  use coamps_statevec_mod,  only : state_vector, get_total_size,              &
                                   state_iterator, get_next, has_next,        &
                                   get_iterator, get_file_type,               &
                                   get_num_fields, initialize_state_vector,   &
                                   dump_state_vector 

  use coamps_statevar_mod,  only : state_variable, gets_update,               &
                                   get_nest_number, get_var_substate,         &
                                   get_record_in_nest, is_nonnegative,        &
                                   get_var_name, get_var_max_level,           &
                                   get_var_min_level, get_mean_flag,          &
                                   get_io_flag, get_vert_type, get_vert_value,&
                                   get_sigma_record, get_mass_level_flag

  use coamps_nest_mod,      only : coamps_nest, get_nest_i_width,             &
                                   get_nest_j_width, get_subnest_iminf,       &
                                   get_subnest_imaxf, get_subnest_jminf,      &
                                   get_subnest_jmaxf, in_this_nest,           &
                                   nest_index_2d_to_1d, get_subnest_i_width,  &
                                   nest_index_1d_to_3d,                       &
                                   get_subnest_j_width, make_nest_point,      &
                                   get_nest_size, get_terrain

  use coamps_intrinsic_mod, only : define_mean_theta,                         &
                                   define_mean_exner

  use coamps_util_mod,      only : C_REAL,                                    &
                                   check_alloc_status,                        &
                                   check_dealloc_status,                      &
                                   check_io_status,                           &
                                   fix_for_platform,                          &
                                   generate_flat_file_name,                   &
                                   read_flat_file,                            &
                                   write_flat_file

  use time_manager_mod,     only : get_date,                                  &
                                   get_time,                                  &
                                   GREGORIAN,                                 &
                                   operator(-),                               &
                                   operator(+),                               &
                                   read_time,                                 &
                                   set_calendar_type,                         &
                                   get_calendar_type,                         &
                                   decrement_date,                            &
                                   increment_date,                            &
                                   set_date,                                  &
                                   set_time,                                  &
                                   time_type,                                 &
                                   write_time

  use types_mod,            only : r4, r8

  use utilities_mod,        only : E_ERR,                                     &
                                   error_handler,                             &
                                   file_exist,                                &
                                   get_unit

  use location_mod,         only : VERTISUNDEF, VERTISSURFACE, VERTISLEVEL,   &
                                   VERTISPRESSURE, VERTISHEIGHT

  implicit none

  private

  !------------------------------
  ! BEGIN PUBLIC INTERFACE
  !------------------------------

  public :: initialize_translator
  public :: finalize_translator

  ! COAMPS restart file tools
  public :: generate_coamps_filenames
  public :: open_coamps_files

  public :: coamps_read_all_fields
  public :: coamps_write_all_fields

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

  public :: print_dart_diagnostics

  !FIX ONLY FOR TESTING
  public :: dart_state

  public :: previous_dtg

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
  integer, parameter :: MAXNEST=7

  integer, parameter :: ONE = 1


  ! Conversion parameters
  integer, parameter :: MIN_TO_SEC = 60
  integer, parameter :: HR_TO_MIN  = 60
  integer, parameter :: DAY_TO_HR  = 24
  integer, parameter :: HR_TO_SEC  = HR_TO_MIN * MIN_TO_SEC
  integer, parameter :: DAY_TO_SEC = DAY_TO_HR * HR_TO_SEC

  ! Index dimensions for the offset matrix
  integer, parameter :: OFFSET_FILE_DIM   = 1
  integer, parameter :: OFFSET_NEST_DIM   = 2

  integer, parameter :: DART_TARGET_TIME  = 1
  integer, parameter :: DART_CURRENT_TIME = 2

  integer, parameter :: COAMPS_HOUR       = 1
  integer, parameter :: COAMPS_MINUTE     = 2
  integer, parameter :: COAMPS_SECOND     = 3

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
  ! number of halo boundary points.  Include the "digital filter" option
  ! as well since we need that to figure byte offsets in the restart file
  integer, dimension(3,MAXNEST):: ktauf = &
  reshape( (/6,0,0,6,0,0,6,0,0,6,0,0,6,0,0,6,0,0,6,0,0/),(/3,MAXNEST/))
  integer, dimension(3)        :: ktaust   = (/ 0,0,0 /)
  integer                      :: icycle   = 6
  integer, dimension(MAXNEST)  :: ndxnam   = (/ 1,1,1,1,1,1,1 /)
  integer, dimension(MAXNEST)  :: ndynam   = (/ 1,1,1,1,1,1,1 /) 
  integer                      :: npr0nam  = 1
  integer                      :: nbnam    = 2
  integer                      :: nbdypt   = 7
  logical                      :: ldigit   = .false.
  character(len=10)            :: cdtg     = '1999083100'
  character(len=11)            :: nml_file = 'convert.nml'
  logical                      :: is_pmo   = .false.
  logical                      :: is_first = .false.
  character(len=180)           :: dsnrff   = './'

  namelist /convert/ ktauf, ktaust, ndxnam, ndynam, npr0nam, nbnam,&
       & cdtg, icycle, is_pmo, dsnrff, is_first, nbdypt, ldigit

  logical  :: FLAT_FILE_IO = .true.
  logical  :: is_dart_async

  ! True if the COAMPS model used an I/O processor when writing out the
  ! restart files. This means:
  ! .true.    =   1 restart file for each nest, regardless of domain decomp
  ! .false.   =   1 restart file for each processor, regardless of nests
  logical :: coamps_used_io_proc

  ! Contents of DART state vector
  type(time_type), dimension(2)                      :: dart_time
  real(kind=r8),   dimension(:), allocatable, target :: dart_state

  ! Contents of COAMPS state vector - this is not the full restart
  ! file - only the fields that we have defined in the dynamic
  ! state vector definition 
  real(kind=C_REAL), dimension(:), allocatable, target :: coamps_state

  ! Arrays allow reading/writing multiple COAMPS files 
  integer                                      :: total_coamps_files
  character(len=64), dimension(:), allocatable :: coamps_file_names
  integer, dimension(:), allocatable           :: coamps_file_units
  character(len=64) :: coamps_file_names_foo

  ! DART restart file
  character(len=11), parameter :: DART_FILENAME = 'dart_vector'
  integer                      :: dart_unit              

  character(len=10)  :: cdtgm1
  character(len=180) :: dsnrff1

  ! The COAMPS domain and the contents of the DART state vector
  type(coamps_domain) :: domain
  type(state_vector)  :: state_layout
  type(state_vector)  :: file_layout

  ! Offsets of each nest within each restart file - this is indexed by (file, nest)
  integer, dimension(:,:), pointer :: restart_nest_offsets

  logical, save :: module_initialized = .false.

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

    character(len=*), parameter :: STATE_VEC_DEF_FILE = 'state.vars'

    if (module_initialized) return
    call set_calendar_type(GREGORIAN)

    ! Get the data we can't read from pre-exiting COAMPS files
    call read_convert_namelist()
    if(is_pmo) then
      cdtgm1 = cdtg
    else
      if(is_first) then
        cdtgm1 = cdtg
      else
        cdtgm1 = previous_dtg(cdtg, icycle)
      end if
    end if

    if (npr0nam > 0) then
      coamps_used_io_proc = .true.
    else
      coamps_used_io_proc = .false. 
    end if

    call initialize_domain(cdtgm1, domain)

    ! Set up the state vector field definitions 
    call initialize_state_vector(state_layout, STATE_VEC_DEF_FILE, domain)
    call initialize_state_vector(file_layout,  STATE_VEC_DEF_FILE, domain, .true.)

    ! is the state contained in flat files or restart files?
    FLAT_FILE_IO = get_file_type(file_layout)

    call allocate_state()

    if (FLAT_FILE_IO) then
      call allocate_coamps_flat_file_info()
    else
      call allocate_coamps_restart_file_info()
      call decompose_domain(domain, npr0nam, ndxnam, ndynam, nbnam)
      call compute_offsets()
    end if

   ! Is this a flat_file io run with/without dart async control?
   is_dart_async = .not.( FLAT_FILE_IO .and. (.not.is_pmo) )

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

    if(.not.FLAT_FILE_IO) then
      deallocate(restart_nest_offsets, stat=dealloc_status)
      call check_dealloc_status(dealloc_status, routine, source,       &
                               revision, revdate, 'restart_nest_offsets' )
    end if
  end subroutine finalize_translator

  ! generate_file_names
  ! --------------------------
  ! Wraper to generate the a list of either COAMPS flat file names 
  ! or restart file names 
  !  PARAMETERS
  !   IN  writing_coamps    True if we are going to be writing the
  !                         COAMPS files, false if reading files
  subroutine generate_coamps_filenames(writing_coamps)
    logical, intent(in)         :: writing_coamps
    if(FLAT_FILE_IO) then
      call generate_all_flat_filenames(writing_coamps)
    else
      call generate_restart_filenames(writing_coamps)
    end if
  end subroutine generate_coamps_filenames

  ! open_coamps_files
  ! -----------------
  ! Wraper to open a list of either COAMPS flat files or restart files
  ! PARAMETERS
  !  IN  writing_coamps     True if we're opening the files for
  !                         write access
  subroutine open_coamps_files(writing_coamps)
    logical, intent(in) :: writing_coamps
    if(FLAT_FILE_IO) then
      call open_coamps_flat_files(writing_coamps)
    else
      call open_coamps_restart_files(writing_coamps)
    end if
  end subroutine open_coamps_files
  
  ! generate_all_flat_filenames
  ! --------------------------
  ! Generate the COAMPS restart file names based on information from
  ! the convert namelist.
  !  PARAMETERS
  !   IN  writing_coamps    True if we are going to be writing the
  !                         COAMPS restart file, false if reading it
  subroutine generate_all_flat_filenames(writing_coamps)
    logical, intent(in)         :: writing_coamps

    !  Variables needed to define the flat file name
    type(state_variable)        :: cur_var
    type(state_iterator)        :: iterator
    integer                     :: coamps_file_index
    integer                     :: nest_number
    integer                     :: level1
    integer                     :: level2
    integer, dimension(3)       :: working_time
    character(len=1)            :: aotype
    character(len=7)            :: field_type
    character(len=3)            :: level_type
    character(len=6)            :: var_name
    integer                     :: i_width
    integer                     :: j_width

    integer                     :: num_lvls

    aotype      = 'a'

    ! Assert that the restart file names are already allocated
    if (.not. allocated(coamps_file_names)) then
       call error_handler(E_ERR, 'generate_flat_filenames', &
                         'COAMPS file names unallocated',      &
                         source, revision, revdate)
    end if

    ! Need to decide the time and the type of file.
    ! If we are going from COAMPS -> DART, we need the COAMPS fcstfld files at the DA interval.  
    ! If we are going 
    ! from DART -> COAMPS, we need the time to be zero and the files to be analfld
    ! DART restart file.
    if (writing_coamps) then
       working_time = (/ 0,0,0 /)
       field_type = 'analfld'
       dsnrff1 = trim(dsnrff)//'analyses/'
    else
       if(is_pmo) then
         working_time = ktauf(:,ONE)
       else
         working_time = (/icycle,0,0/)
       end if
       field_type = 'fcstfld'
       dsnrff1 = dsnrff
    end if

    print *,"COAMPS flat files in use:"
    print *,"-------------------------------------------------"

    coamps_file_index = 0
    iterator = get_iterator(file_layout)
    flat_file_loop: do while (has_next(iterator))
      cur_var = get_next(iterator)
      if( .not. get_io_flag(cur_var)) cycle flat_file_loop

      nest_number = get_nest_number(cur_var)

      var_name = get_var_name(cur_var)

      select case (get_vert_type(cur_var))
       case(VERTISLEVEL)
         level_type = 'sig'
         level1     = int(get_var_max_level(cur_var, domain))
         level2     = int(get_var_min_level(cur_var, domain))
       case(VERTISHEIGHT)
         level_type = 'zht'
         level1     = get_vert_value(cur_var)
         level2     = 0.0_r8
       case(VERTISPRESSURE)
         level_type = 'pre'
         level1     = get_vert_value(cur_var)
         level2     = 0.0_r8
       case(VERTISSURFACE, VERTISUNDEF)
         level_type = 'sfc'
         level1     = 0.0_r8
         level2     = 0.0_r8
      end select

      i_width = get_nest_i_width(get_domain_nest(domain, nest_number))
      j_width = get_nest_j_width(get_domain_nest(domain, nest_number))

      coamps_file_index = coamps_file_index + 1
      call generate_flat_file_name(                                   &
           var_name, level_type, level1, level2, nest_number, aotype, &
           i_width, j_width, cdtg, working_time(COAMPS_HOUR),         &
           working_time(COAMPS_MINUTE), working_time(COAMPS_SECOND),  &
           field_type, coamps_file_names(coamps_file_index) )

      write (*,'(I2.2,2x,A)') coamps_file_index,                      &
                  trim(dsnrff1)//coamps_file_names(coamps_file_index)

    end do flat_file_loop
    print *,"-------------------------------------------------"

    !FIXME Test pressure level output
!    if(.not.writing_coamps) then

!      iterator = get_iterator(file_layout)
!      flat_file_loop_plev: do while (has_next(iterator))
!        cur_var = get_next(iterator)
!        if( get_mean_flag(cur_var)) cycle flat_file_loop_plev
!        nest_number = get_nest_number(cur_var)

!        i_width = get_nest_i_width(get_domain_nest(domain, nest_number))
!        j_width = get_nest_j_width(get_domain_nest(domain, nest_number))

!        var_name = get_var_name(cur_var)
!        level1   = 500
!        level2   = 0
!        level_type = 'pre'

!        call generate_flat_file_name(                                   &
!             var_name, level_type, level1, level2, nest_number, aotype, &
!             i_width, j_width, cdtg, working_time(COAMPS_HOUR),         &
!             working_time(COAMPS_MINUTE), working_time(COAMPS_SECOND),  &
!             field_type, coamps_file_names_foo )

!        write (*,'(4x,A)') trim(dsnrff1)//coamps_file_names_foo

!      end do flat_file_loop_plev
!    print *,"-------------------------------------------------"

!    end if

  end subroutine generate_all_flat_filenames

  ! open_coamps_flat_files
  ! -----------------
  ! Uses the list of COAMPS flat file names to populate the vector of
  ! opened COAMPS flat file units
  ! PARAMETERS
  !  IN  writing_coamps     True if we're opening the files for
  !                         write access
  subroutine open_coamps_flat_files(writing_coamps)
    logical, intent(in) :: writing_coamps

    type(state_variable)        :: cur_var
    type(state_iterator)        :: iterator
    integer                     :: cur_file_index
    integer                     :: cur_file_unit
    integer                     :: rec_len
    integer                     :: r4_length

    character(len=*), parameter :: routine = 'open_coamps_flat_files'
    integer                     :: io_status

    character(len=5)            :: fileaction
    character(len=7)            :: filestatus

    if (writing_coamps) then
       fileaction = 'write'
       filestatus = 'replace'
    else
       fileaction = 'read'
       filestatus = 'old'
    end if

    cur_file_index = 0
    iterator = get_iterator(file_layout)
    flat_file_loop: do while (has_next(iterator))
      cur_var = get_next(iterator)

      if( .not. get_io_flag(cur_var) ) cycle flat_file_loop
      cur_file_index = cur_file_index + 1

      !COAMPS flatfiles are stored as single precesion
      inquire(IOLENGTH=r4_length) 0_r4
      rec_len = get_nest_size(get_domain_nest(domain, get_nest_number(cur_var))) &
                 * get_sigma_record(cur_var) * r4_length

      cur_file_unit = get_unit()
      open( unit=cur_file_unit,                                     &
            file=trim(dsnrff1)//coamps_file_names(cur_file_index),  &
            status=filestatus, access='direct', action=fileaction,  &
            form='unformatted', recl=rec_len, iostat=io_status)
      call check_io_status(io_status, routine, source, revision,    &
                           revdate, 'Opening ' //                   &
                           coamps_file_names(cur_file_index))

      coamps_file_units(cur_file_index) = cur_file_unit
    end do flat_file_loop
  end subroutine open_coamps_flat_files
  
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
       working_time = ktauf(:,ONE)
    end if

    print *,"COAMPS restart files in use:"
    print *,"-------------------------------------------------"
    do coamps_file_index = 1,total_coamps_files
       if (coamps_used_io_proc) then
          write (coamps_file_names(coamps_file_index), 100) &
                 cdtg, working_time, coamps_file_index
       else
          write (coamps_file_names(coamps_file_index), 200) &
                (coamps_file_index-1), cdtg, working_time
       end if
       write(*,*) coamps_file_index, coamps_file_names(coamps_file_index)
    end do
    print *,"-------------------------------------------------"
100 format( 'restarta1p001',A10,I3.3,I2.2,I2.2,'.nest',I1)
200 format( 'restarta1p',I3.3,A10,I3.3,I2.2,I2.2)
  end subroutine generate_restart_filenames

  ! open_coamps_restart_files
  ! -----------------
  ! Uses the list of COAMPS file names to populate the vector of
  ! opened COAMPS restart file units
  ! PARAMETERS
  !  IN  writing_coamps     True if we're opening the files for
  !                         write access
  subroutine open_coamps_restart_files(writing_coamps)
    logical, intent(in) :: writing_coamps

    integer           :: cur_file_index, cur_file_unit

    character(len=*), parameter :: routine = 'open_coamps_restart_files'
    integer           :: io_status

    character(len=5)  :: fileaction

    if (writing_coamps) then
       fileaction = 'write'
    else
       fileaction = 'read'
    end if

    do cur_file_index = 1, total_coamps_files
       cur_file_unit = get_unit()
       open( unit=cur_file_unit,                                  &
             file=coamps_file_names(cur_file_index),              &
             status='old', access='direct', action=fileaction,    &
             form='unformatted', recl=C_REAL, iostat=io_status)
       call check_io_status(io_status, routine, source, revision, &
                            revdate, 'Opening ' //                &
                            coamps_file_names(cur_file_index))

       coamps_file_units(cur_file_index) = cur_file_unit
    end do
  end subroutine open_coamps_restart_files

  ! coamps_read_all_fields
  ! ------------------------------
  ! Wrapper for using coamps_process_all_restart_fields to read COAMPS
  ! restart files
  !  PARAMETERS
  !   [none] 
  subroutine coamps_read_all_fields()
    logical, parameter :: READING_COAMPS = .false.
    if(FLAT_FILE_IO) then
      call coamps_process_all_flat_files(READING_COAMPS)
    else
      call coamps_process_all_restart_fields(READING_COAMPS)
    end if
    call coamps_process_default_vars()
  end subroutine coamps_read_all_fields

  ! coamps_write_all_fields
  ! -------------------------------
  ! Wrapper for using coamps_process_all_restart_fields to write COAMPS
  ! restart files
  !  PARAMETERS
  !   [none]
  subroutine coamps_write_all_fields()
    logical, parameter :: WRITING_COAMPS = .true.
    if(FLAT_FILE_IO) then
      call coamps_process_all_flat_files(WRITING_COAMPS)
    else
      call coamps_process_all_restart_fields(WRITING_COAMPS)
    end if
  end subroutine coamps_write_all_fields

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
    open( unit=dart_unit, file=DART_FILENAME, status=filestatus,        & 
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
       ! DART will output a restart file containing the current time
       ! and the time it wants the model to integrate to
       dart_time(DART_TARGET_TIME) = read_time(dart_unit, 'unformatted')
       call get_time(dart_time(DART_TARGET_TIME), dart_time_seconds, dart_time_days)
       write (*,*) "DART target time read is ", dart_time_days, &
                   " days and ", dart_time_seconds, " seconds."

      if(is_dart_async) then
       dart_time(DART_CURRENT_TIME) = read_time(dart_unit, 'unformatted')
       call get_time(dart_time(DART_CURRENT_TIME), dart_time_seconds, dart_time_days)
       write (*,*) "DART current time read is ", dart_time_days, &
                   " days and ", dart_time_seconds, " seconds."
      end if
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
       call write_time(dart_unit, dart_time(DART_CURRENT_TIME), 'unformatted')
       call get_time(dart_time(DART_CURRENT_TIME), dart_time_seconds, dart_time_days)
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
    integer :: ccyy, mm, dd, hh, dt_sec
    type(time_type) :: t0, dt

    if(FLAT_FILE_IO) then
      ! convert coamps time to dart time
      read(cdtg,100) ccyy,mm,dd,hh
      t0 = set_date(ccyy, mm, dd, hh)
      if(is_pmo) then
        dt_sec = ktauf(COAMPS_HOUR,ONE)*3600 + ktauf(COAMPS_MINUTE,ONE)*60 + ktauf(COAMPS_SECOND,ONE)
      else
        dt_sec = icycle*3600
      end if
      dt = set_time(dt_sec)
      dart_time(DART_CURRENT_TIME) = t0 + dt
    else
      ! Get this information from the namelist
      tau_hour   = ktauf(COAMPS_HOUR,ONE)
      tau_minute = ktauf(COAMPS_MINUTE,ONE)
      tau_second = ktauf(COAMPS_SECOND,ONE)

      call hms_to_sd(tau_hour, tau_minute, tau_second, dart_day,&
           & dart_second)

      !  When writing a DART file, only have one time entry to worry
      ! about.
      dart_time(DART_CURRENT_TIME) = set_time(dart_second, dart_day)
    end if
100 format((I4.4),3(I2.2))
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
    integer :: dart_days
    integer :: dart_seconds
    integer :: ccyy, mm, dd, hh
    type(time_type) :: t0, dt

    if(FLAT_FILE_IO) then
      read(cdtg,100) ccyy,mm,dd,hh

      t0 = set_date(ccyy, mm, dd, hh)
      dt = dart_time(DART_CURRENT_TIME) - t0

      call get_time(dt, dart_seconds, dart_days)
    else
      call get_time(dart_time(DART_CURRENT_TIME), dart_seconds, dart_days)
    end if

    call sd_to_hms(dart_days, dart_seconds, ktaust(COAMPS_HOUR), &
                   ktaust(COAMPS_MINUTE), ktaust(COAMPS_SECOND)) 

100 format((I4.4),3(I2.2))
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
    integer :: dart_days
    integer :: dart_seconds
    integer :: ccyy, mm, dd, hh
    type(time_type) :: t0, dt

    if(FLAT_FILE_IO) then
      read(cdtg,100) ccyy,mm,dd,hh

      t0 = set_date(ccyy, mm, dd, hh)
      dt = dart_time(DART_TARGET_TIME) - t0

      call get_time(dt, dart_seconds, dart_days)
    else
      call get_time(dart_time(DART_TARGET_TIME), dart_seconds, dart_days)
    end if
      call sd_to_hms(dart_days, dart_seconds, ktauf(COAMPS_HOUR,ONE), &
                     ktauf(COAMPS_MINUTE,ONE), ktauf(COAMPS_SECOND,ONE)) 
100 format((I4.4),3(I2.2))
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

    if(.not.is_dart_async) return

    pickup_unit = get_unit()

    open( unit=pickup_unit, file=PICKUP_NAME, status='replace',         &
          form='formatted', iostat=io_status)
    call check_io_status(io_status, routine, source, revision, revdate, &
                         'Opening file ' // PICKUP_NAME)

    write(pickup_unit, fmt=300) ktaust, ktauf(:,ONE)

    close(pickup_unit)

300 format( 5(I3,","),I3 )
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
                            'Reading ' // nml_file // ' convert')
       close(nml_unit)
    else
       call error_handler(E_ERR, 'coamps_translate_mod', 'Convert name' // &
                          'namelist read failed - target file not found',  &
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

    character(len=*), parameter :: routine = 'allocate_state'
    integer :: alloc_status

    ! DART
    allocate(dart_state(get_total_size(state_layout)), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'dart_state')

    ! COAMPS
    allocate(coamps_state(get_total_size(state_layout)), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'coamps_state')
  end subroutine allocate_state

  ! allocate_coamps_restart_file_info
  ! -------------------------
  ! Allocate memory for the COAMPS restart file names and their 
  ! corresponding unit numbers - do this dynamically since the 
  ! number will (potentially) change from run to run based on the
  ! number of processors used.
  !  PARAMETERS
  !   [none]
  subroutine allocate_coamps_restart_file_info()

    character(len=*), parameter :: routine = 'allocate_coamps_restart_file_info'
    integer :: alloc_status

    if (coamps_used_io_proc) then
       ! I/O processor uses 1 restart file per nest
       total_coamps_files = get_nest_count(domain)
    else
       ! Not using I/O processor means one restart file per process
       total_coamps_files = get_num_subdomains(domain)
    end if

    allocate(coamps_file_names(total_coamps_files), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'coamps_file_names')

    allocate( coamps_file_units(total_coamps_files), stat=alloc_status )
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'coamps_file_units')
  end subroutine allocate_coamps_restart_file_info

  ! allocate_coamps_flat_file_info
  ! -------------------------
  ! Allocate memory for the COAMPS flat file names and their 
  ! corresponding unit numbers - do this dynamically since the 
  ! number will (potentially) change from run to run based on the
  ! dimension of the state.
  !  PARAMETERS
  !   [none]
  subroutine allocate_coamps_flat_file_info()

    character(len=*), parameter :: routine = 'allocate_coamps_flat_file_info'
    integer :: alloc_status
    
    total_coamps_files = get_num_fields(file_layout) 

    allocate(coamps_file_names(total_coamps_files), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'coamps_file_names')
    
    allocate( coamps_file_units(total_coamps_files), stat=alloc_status )
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'coamps_file_units')
  end subroutine allocate_coamps_flat_file_info

 ! coamps_process_all_flat_files
  ! ---------------------------------
  ! Reads/writes all fields in the DART restart vector to/from their
  ! proper places in the COAMPS restart file(s)
  !  PARAMETERS
  !   IN  write_coamps      True if we're writing to COAMPS files
  subroutine coamps_process_all_flat_files(write_coamps) 
    logical, intent(in)       :: write_coamps

    type(state_variable)      :: cur_var
    type(state_iterator)      :: iterator

    integer :: cur_file       ! Which file is being read

    logical :: write_field    ! Allow us to not write a field even if
                              ! it's in the DART state vector
    real(kind=C_REAL), dimension(:), pointer :: var_state

    cur_file = 0
    iterator = get_iterator(file_layout)
    flat_file_loop: do while (has_next(iterator))

      cur_var  = get_next(iterator)
      if( .not. get_io_flag(cur_var)) cycle flat_file_loop

      var_state => get_var_substate(cur_var, coamps_state)

      cur_file = cur_file + 1

      if(write_coamps) then
        write_field = gets_update(cur_var)
        if(write_field) then
          call write_flat_file(coamps_file_units(cur_file), var_state)

          !FIXME Write pressure level data here
          ! Need to write out data on specified pressure levels here.

        end if

      else
        call read_flat_file(coamps_file_units(cur_file), var_state)
      end if
    end do flat_file_loop

  end subroutine coamps_process_all_flat_files

  !coamps_process_default_vars
  ! ---------------------------------
  ! adds variables to state vector that are
  ! defined internally
  ! 
  subroutine coamps_process_default_vars()
    type(state_variable)        :: cur_var
    type(state_iterator)        :: iterator
    real(kind=C_REAL), dimension(:), pointer :: var_state


    iterator = get_iterator(file_layout)
    mean_fld_loop: do while (has_next(iterator))
      cur_var   = get_next(iterator)
      if(.not. get_mean_flag(cur_var)) cycle mean_fld_loop

      var_state => get_var_substate(cur_var, coamps_state)
      call calculate_mean_var(cur_var, var_state)
    end do mean_fld_loop

  end subroutine coamps_process_default_vars

  ! compute_offsets
  ! ---------------
  ! Computes the offsets for each nest in use within each restart file
  !  PARAMETERS
  !   [none]
  subroutine compute_offsets()
        
    integer :: nests_per_file
    integer :: cur_file_num, cur_nest_num

    integer :: total_offset

    character(len=*), parameter :: routine = 'compute_offsets'
    integer                     :: alloc_status

    if (coamps_used_io_proc) then
      nests_per_file = 1
    else
      nests_per_file = get_nest_count(domain)
    end if

    allocate(restart_nest_offsets(total_coamps_files, nests_per_file), &
             stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision,   &
                            revdate, 'restart_nest_offsets')
        
    if (coamps_used_io_proc) then
       ! Each nest has its own restart file, so the portion of the file
       ! corresponding to a particular nest is the *whole file*
       restart_nest_offsets(:,:) = 0
    else
       do cur_file_num = 1, total_coamps_files

          total_offset = 0
          do cur_nest_num = 1, nests_per_file
             restart_nest_offsets(cur_file_num, cur_nest_num) = &
             total_offset
             total_offset = total_offset + &
             compute_nest_size_in_file(cur_nest_num, cur_file_num)
          end do
       end do
    end if

    write (*,*) restart_nest_offsets
  end subroutine compute_offsets

  ! coamps_process_all_restart_fields
  ! ---------------------------------
  ! Reads/writes all fields in the DART restart vector to/from their
  ! proper places in the COAMPS restart file(s)
  !  PARAMETERS
  !   IN  write_coamps      True if we're writing to COAMPS files
  subroutine coamps_process_all_restart_fields(write_coamps) 
    logical, intent(in) :: write_coamps

    type(state_iterator) :: iterator
    type(state_variable) :: var_to_process

    iterator = get_iterator(state_layout)
    do while (has_next(iterator))

       var_to_process = get_next(iterator)

       if (write_coamps) then
         if (gets_update(var_to_process)) then
           call coamps_restart_write_field(var_to_process)
         end if
       else
         call coamps_restart_read_field(var_to_process)
       end if
    end do
  end subroutine coamps_process_all_restart_fields

  ! coamps_restart_read_field
  ! -------------------------
  ! Reads a single field out of one or more COAMPS restart file(s)
  !  PARAMETERS
  !   IN  var_to_read     The restart variable to read from the file(s)
  subroutine coamps_restart_read_field(var_to_read)
     type(state_variable), intent(in) :: var_to_read

    integer :: cur_file_num

    do cur_file_num = 1, get_files_per_nest()
      call coamps_restart_read_single_file(var_to_read, cur_file_num)
    end do
  end subroutine coamps_restart_read_field

  ! coamps_restart_write_field
  ! --------------------------
  ! Writes a single field into one or more COAMPS restart file(s)
  !  PARAMETERS
  !   IN  var_to_write    The restart variable to write to the file(s)
  subroutine coamps_restart_write_field(var_to_write)
    type(state_variable), intent(in) :: var_to_write

    integer :: cur_file_num

    do cur_file_num = 1, get_files_per_nest()
      call coamps_restart_write_single_file(var_to_write, cur_file_num)
    enddo
  end subroutine coamps_restart_write_field

  ! coamps_restart_read_single_file
  ! -------------------------------
  ! Reads either a complete (in the case of single-processor I/O)
  ! or partial (for multi-proc I/O) field out of one restart file
  !  PARAMETERS
  !   IN  var_to_read     The restart variable to read from the file     
  !   IN  file_num        The index of the file to read in the units array
  subroutine coamps_restart_read_single_file(var_to_read, file_num)
    type(state_variable), intent(in) :: var_to_read
    integer,              intent(in) :: file_num

    logical :: READING_COAMPS = .false.

    call coamps_restart_do_single_file(var_to_read, file_num, &
                                       READING_COAMPS)
  end subroutine coamps_restart_read_single_file

  ! coamps_restart_write_single_file
  ! --------------------------------
  ! Writes either a complete (in the case of single-processor I/O)
  ! or partial (for multi-proc I/O) field to one restart file
  !  PARAMETERS
  !   IN  var_to_write    The restart variable to write to the file   
  !   IN  file_num        The index of the file unit to write to
  subroutine coamps_restart_write_single_file(var_to_write, file_num)
    type(state_variable), intent(in) :: var_to_write
    integer,              intent(in) :: file_num

    logical :: WRITING_COAMPS = .true.

    call coamps_restart_do_single_file(var_to_write, file_num, &
                                       WRITING_COAMPS)
  end subroutine coamps_restart_write_single_file

  ! coamps_restart_do_single_file
  ! ------------------------------
  ! Computes the position of field data within a COAMPS restart file,
  ! then reads or writes the field section to the file.
  !  PARAMETERS
  !   IN  var_to_process  The restart variable to read or write     
  !   IN  file_num        The index of the file to read from
  !                           I/O processor:    This is the nest number
  !                           No I/O processor: This is the subnest number 
  !   IN  writing_coamps  True if we're writing to the COAMPS restart file
  subroutine coamps_restart_do_single_file(var_to_process, file_index, &
                                           writing_coamps)
    type(state_variable), intent(in) :: var_to_process
    integer,              intent(in) :: file_index
    logical,              intent(in) :: writing_coamps

    ! Index locations in the local subdomain - includes current
    ! locations (ii/jj) and the i/j boundaries based (possibly)
    ! on the domain decomposition 
    integer :: ii,jj     
    integer :: ilbound  
    integer :: iubound 
    integer :: jlbound
    integer :: jubound
    integer :: subnest_index
    integer :: nest_in_file

    ! COAMPS restart file - which one to use and where we are in it
    integer :: restart_unit
    integer :: restart_pos
    integer :: field_nest_record

    type(coamps_nest)                        :: var_nest
    real(kind=C_REAL), dimension(:), pointer :: var_state

    ! Error handling
    character(len=*), parameter :: routine ='coamps_restart_do_single_file'
    logical                     :: file_is_open
    integer                     :: io_status

    ! Recall that for single-process I/O restart files, each nest is
    ! in its own restart file. 
    if (coamps_used_io_proc) then
      restart_unit = coamps_file_units(get_nest_number(var_to_process))
      print *, restart_unit, &
               coamps_file_names(get_nest_number(var_to_process)) 
    else
      restart_unit = coamps_file_units(file_index) 
      print *, restart_unit, coamps_file_names(file_index) 
    end if

    ! Assert
    inquire(unit=restart_unit, opened=file_is_open) 
    if (.not. file_is_open) then
            call error_handler(E_ERR, 'do_single_file', 'COAMPS restart file'&
                               // ' is not open!', source, revision, revdate)
    end if 

    var_nest  =  get_domain_nest(domain, get_nest_number(var_to_process))
    var_state => get_var_substate(var_to_process, coamps_state)

    ! The "subfield" boundaries for the single-processor I/O case
    ! is just the i/j limits for the entire field - multi-processor
    ! needs the domain decomposition data.  
    if (coamps_used_io_proc) then
            ilbound = 1
            iubound = get_nest_i_width(var_nest)
            jlbound = 1
            jubound = get_nest_j_width(var_nest)

            nest_in_file = 1
    else
            subnest_index = file_index
            ilbound = get_subnest_iminf(var_nest, subnest_index)
            iubound = get_subnest_imaxf(var_nest, subnest_index)
            jlbound = get_subnest_jminf(var_nest, subnest_index)
            jubound = get_subnest_jmaxf(var_nest, subnest_index)

      nest_in_file = get_nest_number(var_to_process)
    end if

    ! The stored record is in terms of fields - convert that (using the
    ! number of points in a field) to units of C_REAL.  Also need to take
    ! into account that this nest is possibly buried somewhere in the
    ! middle of the file
    field_nest_record = get_record_in_nest(var_to_process, domain,    &
                                               coamps_used_io_proc, ldigit)
    print *, "Record in nest is", field_nest_record
    restart_pos = (iubound - ilbound + 1) * (jubound - jlbound + 1) * &
                  (field_nest_record - 1) + 1 +                       &
                  restart_nest_offsets(file_index, nest_in_file)

    open(unit=76, file='read.dump', form='formatted', access='append')
    write (76,*) "Bounds:", ilbound, iubound, jlbound, jubound
    do jj = jlbound, jubound
      do ii =ilbound, iubound

        ! If the point is outside the domain, just skip it and move on
        if ( .not. in_this_nest(make_nest_point(var_nest,ii,jj))) then
          print *, "skipping point", ii, jj
          restart_pos = restart_pos + 1
          cycle
        end if

        if (writing_coamps) then
          write(restart_unit, rec=restart_pos, iostat=io_status)&
          var_state(nest_index_2d_to_1d(var_nest, ii, jj))
        else
          read(restart_unit, rec=restart_pos, iostat=io_status)&
               var_state(nest_index_2d_to_1d(var_nest, ii, jj))
        end if
        call check_io_status(io_status, routine, source, revision, &
                             revdate, 'COAMPS restart file')

        write (76, *) restart_pos, ii, jj, &
                      nest_index_2d_to_1d(var_nest, ii, jj), &
                      var_state(nest_index_2d_to_1d(var_nest, ii,jj))

        restart_pos  = restart_pos + 1
      enddo
    enddo
    close(76)
  end subroutine coamps_restart_do_single_file

  ! calculate_mean_var
  ! -------------
  ! Converts a state vector to or from the DART numeric format to or
  ! from the COAMPS format.  The state vectors are stored as module
  ! variables so no need to pass them in as parameters 
  subroutine calculate_mean_var(var, state_slice)
    type(state_variable)                           :: var   
    real(kind=C_REAL), dimension(:), intent(inout) :: state_slice

    real(kind=r8), allocatable, dimension(:,:)     :: mean_var
    real(kind=r8)                                  :: ztop

    integer                                        :: nest_size
    integer                                        :: num_levels

    integer                                        :: alloc_status
    integer                                        :: dealloc_status
    character(len=*),              parameter       :: routine = 'calculate_mean_var'

    nest_size  = get_nest_size(get_domain_nest(domain,get_nest_number(var)))
    num_levels = get_domain_num_levels(domain)  

    if(.not.get_mass_level_flag(var)) num_levels = num_levels + 1

    allocate(mean_var(nest_size, num_levels), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'mean_var--'//trim(get_var_name(var)))

    ztop = get_domain_wsigma(domain, 1)

    select case (trim(get_var_name(var)))
      case ('EXBM')
        call define_mean_exner(mean_var, nest_size, num_levels, ztop,  &
             reshape(get_terrain(get_domain_nest(domain, get_nest_number(var))), (/ nest_size /)),      &
             get_domain_msigma(domain))

      case ('THBM')
        call define_mean_theta(mean_var, nest_size, num_levels, ztop,  &
             reshape(get_terrain(get_domain_nest(domain, get_nest_number(var))), (/ nest_size /)),      &
             get_domain_msigma(domain))

      case ('EXBW')
        call define_mean_exner(mean_var, nest_size, num_levels, ztop,  &
             reshape(get_terrain(get_domain_nest(domain, get_nest_number(var))), (/ nest_size /)),      &
             get_domain_wsigma(domain))

    end select

    state_slice = reshape(real(mean_var, kind=C_REAL), (/ nest_size * num_levels /))

    deallocate(mean_var, stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'mean_var--'//trim(get_var_name(var)))

  end subroutine calculate_mean_var

  ! convert_state
  ! -------------
  ! Converts a state vector to or from the DART numeric format to or
  ! from the COAMPS format.  The state vectors are stored as module
  ! variables so no need to pass them in as parameters 
  ! PARAMETERS
  !  IN  to_dart            True if we're converting to DART format
  subroutine convert_state(to_dart)
    logical, intent(in) :: to_dart

    ! If COAMPS used an I/O processor, then the restart file written
    ! out and read in should be big-endian.

    if (to_dart) then
      ! big_endian -> little-endian
      if (coamps_used_io_proc .and. .not. FLAT_FILE_IO) call change_coamps_endian()
      !if (coamps_used_io_proc) call change_coamps_endian()
      dart_state = real(coamps_state,kind=r8)
    else
      coamps_state = real(dart_state,kind=C_REAL)
      ! little-endian -> big_endian
      !if (coamps_used_io_proc) call change_coamps_endian()
      if (coamps_used_io_proc .and. .not. FLAT_FILE_IO) call change_coamps_endian()
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

    real(kind=r8), dimension(:), pointer :: state_subsect
    type(state_iterator) :: var_iterator

    type(state_variable) :: cur_var

    var_iterator = get_iterator(state_layout)
    do while (has_next(var_iterator))
       cur_var = get_next(var_iterator)

       if (is_nonnegative(cur_var)) then
         state_subsect => get_var_substate(cur_var, dart_state)
         where (state_subsect < 0) state_subsect = 0     
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
    call fix_for_platform(coamps_state, size(coamps_state))
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

  ! previous_dtg
  ! ---------
  ! Decrements the current date time group by one DA cycle
  !  PARAMETERS
  !   IN  cdtg         current date time group
  function previous_dtg(dtg, icycle_loc) result(dtgm1)
    character(len=10), intent(in) :: dtg  
    integer,           intent(in) :: icycle_loc
    character(len=10)             :: dtgm1
    type(time_type)               :: t0, tm1, dt
    integer                       :: ccyy, mm, dd, hh, nn, ss

    ! convert coamps time to dart time
    read(dtg,100) ccyy,mm,dd,hh
    t0 = set_date(ccyy, mm, dd, hh)

    ! decrement the time by one DA cycle
    ss = 3600*icycle_loc
    dt = set_time(ss) 
    tm1 = t0 - dt
  
    ! convert decremented time back to ccyymmddhh form
    call get_date(tm1, ccyy, mm, dd, hh, nn, ss)
    write(dtgm1,100) ccyy,mm,dd,hh

100 format((I4.4),3(I2.2))
    end function previous_dtg

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
    integer :: ii

    real(kind=r8), dimension(:), pointer :: var_state
    integer,       dimension(3)          :: max_ijk, min_ijk
    integer                              :: nest_number
    type(state_iterator)                 :: iterator
    type(state_variable)                 :: cur_var
    type(coamps_nest)                    :: cur_nest

    print *, "Values in DART State Vector"

    iterator  = get_iterator(file_layout)
    do while(has_next(iterator))

      cur_var = get_next(iterator)
      var_state => get_var_substate(cur_var, dart_state)
      max_value = maxval(var_state)
      min_value = minval(var_state)

      max_location = maxloc(var_state)
      min_location = minloc(var_state)

      nest_number = get_nest_number(cur_var)
      cur_nest = get_domain_nest(domain, nest_number)

      max_ijk = nest_index_1d_to_3d(cur_nest, max_location(1))
      min_ijk = nest_index_1d_to_3d(cur_nest, min_location(1))

      write (*,1602) get_var_name(cur_var), 'max/min on nest', nest_number, ':', &
                     max_value, ' (',max_ijk,')',' / ', &
                     min_value, ' (',min_ijk,')'
      !write (*,1600) get_var_name(cur_var), "maximum is", max_value, "at index ", max_ijk
      !write (*,1601) "minimum is", min_value, "at index ", min_ijk
    end do
    write (*,*) "First entry in state vector is", dart_state(1)

1600 format(A6,1x,A10,1x,F14.5,1x,A8,3(1x,I3))
1601 format(7x,A10,1x,F14.5,1x,A8,3(1x,I3))
1602 format(A6,1x,A,1x,I1,1x,A1,1x,F18.5,1x,A,3(1x,I3),2A,F18.5,1x,A,3(1x,I3),A)
    

    !open(76,file='state.dump', form='formatted')
    !do ii=1,size(dart_state)
    !  write (76,*) dart_state(ii), coamps_state(ii)
    !end do
    !close(76)
  end subroutine print_dart_diagnostics
    
  ! compute_nest_size_in_file
  ! -------------------------
  ! Computes the size of a given nest within a given multiprocessor restart
  ! file - this will vary between files depending on the domain decomp.
  ! The answer is given in units of C_REAL
  ! ***WARNING*** This assumes that the size of REALs and INTEGERs are
  ! the same within the restart file (since we're using C_REAL as the
  ! unformatted file record length)
  !  PARAMETERS
  !   IN  nest            The nest to compute offsets from   
  function compute_nest_size_in_file(nest_num, file_num)
    integer, intent(in)  :: nest_num
    integer, intent(in)  :: file_num
    integer              :: compute_nest_size_in_file

    type(coamps_nest) :: nest
    integer           :: subnest_num

    integer :: num_levels
    integer :: num_subnests

    integer :: i_width
    integer :: j_width
    integer :: ij_area
    integer :: bc_size

    integer :: lateral_bc_points

    ! Get this particular nest and define which subdomain we're
    ! processing
    nest        = get_domain_nest(domain, nest_num)
    subnest_num = file_num

    ! These are constant between all nests
    num_levels   = get_domain_num_levels(domain) 
    num_subnests = get_num_subdomains(domain)

    ! Get the size of this subdomain
    i_width = get_subnest_i_width(nest, subnest_num)
    j_width = get_subnest_j_width(nest, subnest_num)
    ij_area = i_width * j_width
    bc_size = ij_area
                            
    ! Boundary points for this nest are stored only on the root process
    if (file_num .eq. 1) then
        lateral_bc_points = (2 * get_nest_i_width(nest) * nbdypt) + &
                            (2 * get_nest_j_width(nest) * nbdypt) - &
                            (4 * nbdypt                 * nbdypt) 
    else
        lateral_bc_points = 1
    end if

    ! This is based on counting in write_restart.F in COAMPS amlib
    compute_nest_size_in_file =  112 *  ij_area                           &
                               + 28  * (ij_area * (num_levels    ))       &
                               + 20  * (ij_area * (num_levels + 1))       &
                               + 25  * (bc_size * (num_levels    ))       &
                               + 2   * (bc_size * (num_levels + 1))       &
                               + 6   *  lateral_bc_points                 &
                               + 4   *  bc_size                           &
                               + 6   * (lateral_bc_points * num_subnests) &
                               + 3   *  num_subnests                      &
                               + 3   *  bc_size                           &
                               + 1   * (i_width * 2)                      &
                               + 1   * (j_width * 2)                      &
                               + 1   *  ij_area          

    ! Not using the digital filter requires having a few more entries
    if (.not. ldigit) then
        compute_nest_size_in_file =  compute_nest_size_in_file  & 
                                   + 28 * (ij_area * num_levels)
    end if

    write (*,*) "Calculated nest size for nest", nest_num, "in file", &
                file_num
    write (*,*) "     ", ij_area, num_levels, lateral_bc_points
    write (*,*) "     ", compute_nest_size_in_file, &
                         compute_nest_size_in_file*C_REAL
  end function compute_nest_size_in_file    

  ! get_files_per_nest
  ! ------------------
  ! Returns the number of files required to obtain an entire nest's worth
  ! of values - this is based on the domain decomposition and the use of
  ! a single I/O processor
  !  PARAMETERS
  !   [none]
  function get_files_per_nest()
    integer :: get_files_per_nest

    if (coamps_used_io_proc) then
      get_files_per_nest = 1
    else
      get_files_per_nest = get_num_subdomains(domain)
    end if
  end function get_files_per_nest
    
  !------------------------------
  ! END PRIVATE ROUTINES
  !------------------------------

end module coamps_translate_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
