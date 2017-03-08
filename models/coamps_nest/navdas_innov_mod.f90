! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

module navdas_innov_mod

!------------------------------
! MODULE:       navdas_innov_mod
! AUTHOR:       P. A. Reinecke
!               Naval Research Laboratory
!
! Module containing routines to process a navdas generate ascii
! innovation file
!------------------------------ 

  use location_mod,    only : location_type,    &
                              set_location,     &
                              VERTISUNDEF,      &
                              VERTISSURFACE,    &
                              VERTISLEVEL,      &
                              VERTISPRESSURE,   &
                              VERTISHEIGHT 

  use time_manager_mod,only : time_type,                                   &
                              operator(==),                                &
                              increment_time,                              &
                              decrement_time,                              &
                              set_date,                                    &
                              time_manager_init,                           &
                              set_calendar_type

  use obs_def_mod,     only : obs_def_type,                                &
                              init_obs_def

  use obs_sequence_mod,only : obs_sequence_type,                           &
                              obs_type,                                    &
                              insert_obs_in_seq,                           &
                              static_init_obs_sequence,                    &
                              destroy_obs_sequence,                        &
                              init_obs,                                    &
                              set_qc,                                      &
                              destroy_obs,                                 &
                              set_obs_def,                                 &
                              set_obs_values

  use obs_kind_mod,    only : get_index_for_type_of_obs,                          &
                              get_quantity_for_type_of_obs,                       &
                              QTY_VORTEX_LAT, QTY_VORTEX_LON

  use obs_err_mod,     only : rawin_temp_error,                            &
                              rawin_wind_error,                            &
                              acars_wind_error,                            &
                              acars_temp_error,                            &
                              fixed_marine_pres_error,                     &
                              land_pres_error

  use utilities_mod,   only : E_ERR,                                       &
                              E_MSG,                                       &
                              error_handler,                               &
                              file_exist,                                  &
                              get_unit,                                    &
                              find_namelist_in_file,                       &
                              check_namelist_read

  use coamps_util_mod, only : check_alloc_status,                          &
                              check_dealloc_status,                        &
                              check_io_status,                             &
                              uppercase

  use types_mod,       only : r8

  use coamps_intrinsic_mod, only : compute_altimeter

  implicit none

  private

  !------------------------------
  ! BEGIN PUBLIC INTERFACE
  !------------------------------

  public :: init_navdas_innov_mod,        &
            terminate_navdas_innov_mod,   &
            open_innov_file,              &
            read_innov_header,            &
            read_innov_data,              &
            read_ngt_file,                &
            open_ngt_file,                &
            get_max_obs

  public :: seq, num_copies, num_qc, obs_seq_in_name, innov_file_name
  
  !------------------------------
  ! END PUBLIC INTERFACE
  !------------------------------

  !------------------------------
  ! BEGIN EXTERNAL INTERFACE
  !------------------------------

  !------------------------------
  ! END EXTERNAL INTERFACE
  !------------------------------

  !------------------------------
  ! BEGIN TYPES AND CONSTANTS
  !------------------------------

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

  character(len=100) :: format_header(5)
  character(len=100) :: format_data

  integer :: max_num_instruments
  integer :: max_num_vars
  integer :: max_num_obs

  integer, allocatable :: instrument_key(:)

  integer, parameter :: instrument_list_len=30, &
                        variable_list_len=10
  character(len=variable_list_len),   allocatable :: variable_list(:)
  character(len=instrument_list_len), allocatable :: instrument_list(:)

  logical :: is_initialized=.false.
  logical :: ngt_exists=.false.

  integer, parameter :: num_copies=1, num_qc=1

  integer, parameter :: SECONDS_PER_HOUR=3600

  type(time_type)         :: time_base
  type(obs_def_type)      :: obs_def
  type(obs_sequence_type) :: seq
  type(obs_type)          :: obs

  character(len=129) :: innov_file_name   = 'innov.out'
  character(len=129) :: ngt_file_name     = 'ngt.out'
  character(len=129) :: obs_seq_in_name   = 'obs_seq.out'
  integer            :: obs_window = 1800
  namelist /navdas_innov_nml/ innov_file_name, obs_seq_in_name, obs_window, ngt_file_name

  integer :: innov_unit, ngt_unit

  integer :: alloc_status, dealloc_status, io_status

  real(kind=r8), parameter :: CONVERT_MB_TO_PA=100.0
  !------------------------------
  ! END MODULE VARIABLES
  !------------------------------

contains

  !------------------------------
  ! BEGIN PUBLIC ROUTINES
  !------------------------------
  
  ! init_navdas_innov_mod
  ! -------------------
  ! initializes this module.
  !  PARAMETERS [NONE]
  subroutine init_navdas_innov_mod()
    integer :: nml_unit

    is_initialized = .true.

    call find_namelist_in_file('input.nml', 'navdas_innov_nml', nml_unit)
    read (nml_unit,nml=navdas_innov_nml,iostat=io_status)
    call check_namelist_read(nml_unit, io_status, 'navdas_innov_nml')

    call set_calendar_type('GREGORIAN')

    call static_init_obs_sequence()
      
    call init_obs(obs, num_copies, num_qc)

    call set_innov_file_format()
    call set_variable_list()
    call set_instrument_list()
  end subroutine init_navdas_innov_mod

  ! terminate_navdas_innov_mod
  ! -------------------
  ! terminates this module.
  !  PARAMETERS [NONE]
  subroutine terminate_navdas_innov_mod()
    character(len=*), parameter :: routine='terminate_navdas_innov_mod'
    is_initialized = .false.

    close(innov_unit)
    close(ngt_unit)

    call destroy_obs_sequence(seq)
    call destroy_obs(obs)

    deallocate( variable_list, stat=dealloc_status )
    call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'variable_list')

    deallocate( instrument_list, stat=dealloc_status )
    call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'instrument_list')

    deallocate( instrument_key, stat=dealloc_status )
    call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'instrument_key')
  end subroutine terminate_navdas_innov_mod

  ! open_innov_file
  ! -------------------
  ! Opens an innovation file
  subroutine open_innov_file()
    character(len=*), parameter :: routine='open_innov_file'
    if(.not.is_initialized) call init_navdas_innov_mod()
    innov_unit = get_unit()
    if (file_exist(trim(innov_file_name))) then
      open(innov_unit, file=trim(innov_file_name), status='old', iostat=io_status)
      call check_io_status(io_status, routine, source, revision, revdate,   &
                           'Opening ' // trim(innov_file_name))
    else
      call error_handler(E_ERR, 'innov_to_obs_seq',  &
                         trim(innov_file_name) // ' not found', &
                         source, revision, revdate)
      call terminate_navdas_innov_mod()
      stop
    end if
  end subroutine open_innov_file

  ! open_ngt_file
  ! -------------------
  ! Opens an ngt file
  subroutine open_ngt_file()
    character(len=*), parameter :: routine='open_ngt_file'
    if(.not.is_initialized) call init_navdas_innov_mod()
    ngt_unit = get_unit()
    if (file_exist(trim(ngt_file_name))) then
      open(ngt_unit, file=trim(ngt_file_name), status='old', iostat=io_status)
      call check_io_status(io_status, routine, source, revision, revdate,   &
                           'Opening ' // trim(ngt_file_name))
      ngt_exists = .true.
    else
      call error_handler(E_MSG, 'innov_to_obs_seq',  &
                         trim(ngt_file_name) // ' not found', &
                         source, revision, revdate)
      ngt_exists = .false.
    end if
  end subroutine open_ngt_file

  ! read_innov_data
  ! -------------------
  ! reads the full innovatios file and fill the module level obs_seq variable
  !  PARAMETERS  [NONE]
  subroutine read_innov_data()
    integer :: n
    logical :: is_last, is_ob_defined

    if(.not.is_initialized) call init_navdas_innov_mod()
    do n=1,max_num_obs
      call read_innov_line(is_ob_defined,is_last)
      if(is_ob_defined) call insert_obs_in_seq(seq, obs)
    end do
  end subroutine 

  ! read_innov_header
  ! -------------------
  ! reads the header of the innovation file sets the module level variables
  ! max_num_obs and time_base .
  !  PARAMETERS  [NONE]
  subroutine read_innov_header()
    character(len=*), parameter :: routine='read_innov_header'

    integer, parameter :: num_in=20, str_len=20
    integer      :: igrid, iref, jref, im, jm, lm, int_in(num_in), tau, &
                   yyyy, mm, dd, hh
    real(kind=r8) :: reflat, reflon, stdlt1, stdlt2, stdlon, delx, dely
    character(len=str_len) :: str(num_in)
    character(len=10)      :: cdtg_bg

    real(kind=r8), allocatable :: pressure(:)

    if(.not.is_initialized) call init_navdas_innov_mod()

    read(innov_unit,format_header(1)) &
       str(1),igrid,str(2),iref,str(3),jref,str(4),im,str(5),jm &
      ,str(6),lm,str(7),reflat,str(8),reflon,str(9),stdlt1,str(10) &
      ,stdlt2,str(11),stdlon,str(12),delx,str(13),dely

    allocate(pressure(lm), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'pressure')

    read(innov_unit,format_header(2)) str(1), pressure
    read(innov_unit,format_header(3)) str(1),int_in(1),str(2),int_in(2)
    read(innov_unit,format_header(4)) str(1),max_num_obs,str(2),cdtg_bg,str(4),tau
    read(innov_unit,format_header(5)) str(1:19)

    read(cdtg_bg(1:4),'(I)')  yyyy
    read(cdtg_bg(5:6),'(I)')  mm
    read(cdtg_bg(7:8),'(I)')  dd
    read(cdtg_bg(9:10),'(I)') hh

    time_base = set_date(yyyy, mm, dd, hh)
    time_base = increment_time(time_base,tau*SECONDS_PER_HOUR) 

    deallocate(pressure,stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'pressure')
    return
  end subroutine
  ! get_max_obs
  ! -------------------
  ! returns the number of observations in the innovation file
  function get_max_obs() result(max_obs)
    integer :: max_obs
    max_obs=max_num_obs+100
    return
  end function get_max_obs

  !------------------------------
  ! END PUBLIC ROUTINES
  !------------------------------

  !------------------------------
  ! BEGIN PRIVATE ROUTINES
  !------------------------------

  ! read_ngt_file
  ! -------------------
  subroutine read_ngt_file()

    type(time_type)      :: time_ob
    type(obs_def_type)   :: obs_def
    type(location_type)  :: ob_loc

    integer :: nobs, yy, mm, dd, hh
    integer :: i, ierr 

    character(len=1) :: NS, EW

    real(kind=r8) :: ob_lat, ob_lon, ob_err, ob_qc
    real(kind=r8) :: max_wind

    character(len=13), parameter :: ngt_header_fmt='(I2,1x,4(I2))' 
    character(len=28), parameter :: ngt_data_fmt='(F3.1,A1,1x,F4.1,A1,1x,F3.0)'
    !character(len=32), parameter :: ngt_data_fmt='(F3.1,A1,1x,F4.1,A1,5x,A2,1x,A1)'

    if(.not.ngt_exists) return

    do while(.true.)
      read(ngt_unit, ngt_header_fmt,iostat=ierr) nobs, yy, mm, dd, hh
      if(ierr.ne.0) return

      time_ob = set_date(2000+yy, mm, dd, hh)
      if (time_ob == time_base) exit

      do i=1,nobs
        read(ngt_unit, ngt_data_fmt) ob_lat, NS, ob_lon, EW, max_wind 
      end do

    end do

    do i=1,nobs
      read(ngt_unit, ngt_data_fmt) ob_lat, NS, ob_lon, EW, max_wind 
      if(NS == 'S') ob_lat = -ob_lat
      if(EW == 'W') ob_lon = 360.0_r8 - ob_lon

      ob_loc  = set_location(ob_lon, ob_lat, 0.0_r8, VERTISSURFACE)

      if(max_wind <= 35.0_r8) then
        ob_err = 1.0_r8
      else if (max_wind <= 50.0_r8 .and. max_wind > 35.0_r8) then
        ob_err = 0.5_r8
      else
        ob_err = 0.25_r8
      end if
 
      ob_err = ob_err**2
      ob_qc  = 0.0_r8

      ! Set the vortex lat location
      call init_obs_def(obs_def, ob_loc, get_index_for_type_of_obs('VORTEX_LAT'),time_ob, ob_err) 
      call set_obs_def(obs, obs_def) 
      call set_obs_values(obs, (/ob_lat/))
      call set_qc(obs, (/ob_qc/))

      call insert_obs_in_seq(seq, obs)

      ! Set the vortex lon location
      call init_obs_def(obs_def, ob_loc, get_index_for_type_of_obs('VORTEX_LON'),time_ob, ob_err) 
      call set_obs_def(obs, obs_def) 
      call set_obs_values(obs, (/ob_lon/))
      call set_qc(obs, (/ob_qc/))

      call insert_obs_in_seq(seq, obs)
    end do

    return
  end subroutine read_ngt_file

  ! read_innov_line
  ! -------------------
  ! reads a single line an open innovation file assigns the data
  ! to the module level obs_type structure.
  !  PARAMETERS 
  !   OUT  is_last_line  	true if an EOF signal was reached			
  subroutine read_innov_line(is_ob_defined, is_last_line)
    logical, intent(out) :: is_last_line
    logical, intent(out) :: is_ob_defined

    type(time_type)      :: time_ob
    type(obs_def_type)   :: obs_def
    type(location_type)  :: ob_loc

    character(len=32) :: ob_type

    real(kind=r8)      :: ob_value(1), ob_bk, t_bk, iv, ob_err, etc, &
                         ob_lat, ob_lon, ob_lev, q_bk
    integer           :: nn, vtype, itype, nvp, ob_qc(1), ob_dt, idp
    integer           :: vert_level, ob_kind_indx, ob_type_indx
    character(len=11) :: ob_database
    character(len=17) :: ob_platform

    is_last_line=.false.
    is_ob_defined=.false.
    if(.not.is_initialized) call init_navdas_innov_mod()

    read(innov_unit,format_data) nn,ob_value,ob_bk,t_bk,iv,ob_err,etc,         &
                                 ob_lat,ob_lon,ob_lev,vtype,itype,nvp,ob_qc,   &
                                 ob_dt,ob_platform,ob_database,idp,q_bk

    !All navdas qc flags less than zero are good!
    ob_qc=max(ob_qc(1),0)

    ! ob_err is std, need to give variance 

 
    ! deal with special cases for surface based observations.
    select case (itype) 
      case(1, 2) ! Land observations (coastal, manual, automated) and Land observations (metar)

        if(vtype /= 1) then ! Obs err very large for all other variables
          is_ob_defined = .false. 
          return
        end if

        if(vtype==1) then 
          ob_err   = land_pres_error(ob_lev)*CONVERT_MB_TO_PA
          ob_value = ob_lev*CONVERT_MB_TO_PA ! For surface pressure
          if(etc.gt.1) ob_value = compute_altimeter(ob_value(1), etc) ! For altimeter
        end if
        ob_lev = etc

      case(10  ) ! Surface obs from ships, fixed and mobile, drifting buoys
        if(vtype /= 1 .and. vtype /= 3 .and. vtype /= 4) then 
          is_ob_defined = .false. 
          return
        end if

        if(vtype==1) then 
          ob_err   = fixed_marine_pres_error(ob_lev)*CONVERT_MB_TO_PA
          ob_value = ob_lev*CONVERT_MB_TO_PA ! For surface pressure
          if(etc.gt.1) ob_value = compute_altimeter(ob_value(1), etc) ! For altimeter
        end if

        if(vtype == 3 .or. vtype == 4) then
          ob_lev = 10.0_r8 ! Assume surface winds are at 10 m.
        else
          ob_lev = etc
        end if

      case(60,61,70,71,72,73,80,81,82,83,84,85,86,87) ! ssmi winds, wind sat, and scat winds
        ob_lev = 10.0_r8 ! sattelite winds are supposed to be 10-m winds

      case(250, 251) ! ssmi and nrl windsat TPW
        ob_lev = 0.0_r8 ! Does this matter?

    end select

    ob_err=ob_err**2

    ! Since dart currently doesn't support a finite-time window
    ! around the observations, fudge it here.
    if(abs(ob_dt) .le. obs_window) ob_dt=0

    call set_ob_type(vtype, itype, ob_platform, ob_type)
    if(ob_dt .lt. 0) then
      time_ob = decrement_time(time_base,abs(ob_dt)) 
    else
      time_ob = increment_time(time_base,abs(ob_dt)) 
    end if

    vert_level = get_vert_level_type(itype, vtype) 
    
    ! navdas innov vector has pressure in MB as vertical coord.  Dart expects PA.
    if(vert_level == VERTISPRESSURE) ob_lev = ob_lev*CONVERT_MB_TO_PA

    ob_loc  = set_location(ob_lon, ob_lat, ob_lev, vert_level)
    ob_type_indx = get_index_for_type_of_obs(ob_type) 

    if(ob_type_indx <= 0) return

    ob_kind_indx=get_quantity_for_type_of_obs(ob_type_indx)

    call init_obs_def(obs_def, ob_loc, ob_type_indx, time_ob, ob_err) 
    call set_obs_def(obs, obs_def) 
    call set_obs_values(obs, ob_value)
    call set_qc(obs, real(ob_qc,kind=r8))
  
    is_ob_defined=.true.
    return
  end subroutine read_innov_line

  ! set_ob_type
  ! -------------------
  ! sets the name of the ob_type for a observation
  !  PARAMETERS 
  !   IN  vtype			variable key
  !   IN  itype			instrument key
  !   OUT ob_type		ob_type character string	
  subroutine set_ob_type(vtype, itype, ob_platform, ob_type) 
    integer, intent(in) :: vtype, itype
    character(len=*), intent(in)  :: ob_platform 
    character(len=*), intent(out) :: ob_type

    character(len=16) :: instrument_type 
    character(len=8)  :: variable_type 
    logical           :: is_drop

    is_drop = .false.
    if(uppercase(trim(adjustl(ob_platform(1:5)))) .eq. 'RECO') is_drop = .true. 

    call get_ob_variable_from_indx(vtype, variable_type)

    if(is_drop) then
      instrument_type = 'DROPSONDE'
    else
      call get_instrument_from_indx(itype, instrument_type)
    end if

    select case (itype)
      case(1, 2, 10)
        if(vtype==1) variable_type = 'P'
      case(60, 61) ! Use WS for wind speed obs instead of FF found in .inc file
        if(vtype==8) variable_type = 'WS'
    end select

    ob_type=trim(adjustl(variable_type))//' '//trim(adjustl(instrument_type))

    call replace_white_space(ob_type,'_')
    ob_type=uppercase(ob_type)
    return

  end subroutine set_ob_type
  
  ! get_ob_variable_from_indx
  ! -------------------
  ! sets variable label according to integer key  
  !  PARAMETERS 
  !   IN  vtype			variable key
  !   OUT variable		variable character string	
  subroutine get_ob_variable_from_indx(vtype,variable)
    integer, intent(in) :: vtype
    character(len=*), intent(out) :: variable
    variable=variable_list(vtype)
    return
  end subroutine get_ob_variable_from_indx
  
  ! get_instrument_from_indx
  ! -------------------
  ! sets variable label according to integer key  
  !  PARAMETERS 
  !   IN  itype			instrument key
  !   OUT instrument	instrument character string	
  subroutine get_instrument_from_indx(itype,instrument)
    integer, intent(in) :: itype
    character(len=*), intent(out) :: instrument
    integer i
    do i=1,max_num_instruments
      if(itype == instrument_key(i)) exit
    end do
    instrument=instrument_list(i)
    return
  end subroutine get_instrument_from_indx
  ! set_innov_file_format
  ! -------------------
  ! sets the string format descriptors provided from navdas include file
  !  PARAMETERS  [NONE]
  subroutine set_innov_file_format()
    include 'xiv_format_default.inc'
    format_header(1) = fmt_hdr1
    format_header(2) = fmt_hdr2
    format_header(3) = fmt_hdr3
    format_header(4) = fmt_hdr4
    format_header(5) = fmt_hdr5
    format_data      = fmt_xiv
    return
  end subroutine set_innov_file_format

  ! set_variable_list
  ! -------------------
  ! sets the list of variables that may occur in innovation file
  ! include file is provided from navdas code
  !  PARAMETERS [NONE]
  subroutine set_variable_list()
    include 'var_list_coamps.inc'
    character(len=*), parameter :: routine='set_variable_list'
    integer :: ii
    max_num_vars=mx_nm_var
    allocate( variable_list(max_num_vars), stat=alloc_status )
    call check_alloc_status(alloc_status, routine, source, revision, &
                           revdate, 'variable_list')
    do ii=1,max_num_vars
      variable_list(ii) = c_label_var(ii)
    end do
  end subroutine set_variable_list
    
  ! set_instrument_list
  ! -------------------
  ! sets the list of instruments that may occur the innovation file
  ! include file is provided from navdas code
  !  PARAMETERS [NONE]
  subroutine set_instrument_list()
    include 'instru_list_coamps.inc'
    character(len=*), parameter :: routine='set_instrument_list'
    character(len=instrument_list_len) :: instrument_tmp
    integer :: ii,i_len
    integer :: i_under, i_dash, i_var
    max_num_instruments = mx_nm_instru
    allocate( instrument_list(max_num_instruments), stat=alloc_status )
    call check_alloc_status(alloc_status, routine, source, revision, &
                           revdate, 'instrument_list')
    allocate( instrument_key(max_num_instruments), stat=alloc_status )
    call check_alloc_status(alloc_status, routine, source, revision, &
                           revdate, 'instrument_key')
    do ii=1,max_num_instruments 
      instrument_tmp = c_label(ii)

      ! Remove TPPW from end of names
      i_var=index(trim(instrument_tmp),'TPPW')-1
      if(i_var>0) instrument_tmp = instrument_tmp(1:i_var)

      ! cuttoff instrument names at the first '-' or '_'
      i_under=index(trim(instrument_tmp),'_')
      i_dash=index(trim(instrument_tmp),'-')
      if(i_dash==0) i_dash=instrument_list_len
      if(i_under==0) i_under=instrument_list_len
      i_len = min(i_under,i_dash)-1
      if(i_len.lt.1) then
        instrument_list(ii) = instrument_tmp
      else
        instrument_list(ii) = instrument_tmp(1:i_len)
      end if
      instrument_key(ii) = nc(ii)
    end do
  end subroutine set_instrument_list
  
  ! get_vert_level_type
  ! -------------------
  ! set the level type of the observations
  ! sets the list of instruments that may occur the innovation file
  ! include file is provided from navdas code
  !  PARAMETERS 
  !   IN  itype	 			key for instrument type
  !   OUT get_vert_level	vertical level key
  integer function get_vert_level_type(itype, vtype) 
    integer, intent(in)  :: itype
    integer, intent(in), optional :: vtype
    select case (itype)
      case(1,2,10) ! sfc land, metar, ship
        if(present(vtype)) then
          select case (vtype)
            case(2,3,4) 
              get_vert_level_type=VERTISHEIGHT
            case default
              get_vert_level_type=VERTISSURFACE
            end select
        else
          get_vert_level_type=VERTISSURFACE
        end if

      case(60,61,70,71,72,73,80,81,82,83,84,85,86,87) ! ssmi winds, wind sat, and scat winds
        get_vert_level_type=VERTISHEIGHT
      case(250, 251) ! ssmi TPW, nrl windsat TPW
        get_vert_level_type=VERTISUNDEF
      case default
        get_vert_level_type=VERTISPRESSURE
    end select
    return
  end function get_vert_level_type

  ! replace_white_space
  ! -------------------
  ! replaces the white space in a string with str_fill
  subroutine replace_white_space(str_in,str_fill)
    character(len=*), intent(inout) :: str_in
    character(len=*), intent(in)    :: str_fill
    integer :: ibeg
    str_in=adjustl(str_in)
    do while(.true.)
      ibeg=index(trim(str_in),' ')
      if(ibeg == 0) exit
      str_in(ibeg:ibeg) = trim(str_fill)
    end do
    return
  end subroutine replace_white_space

  !------------------------------
  ! END PRIVATE ROUTINES
  !------------------------------
end module navdas_innov_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
