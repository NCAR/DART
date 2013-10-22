! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

program update_datahd_file

  use coamps_util_mod,      only : check_io_status,              &
                                   read_datahd_file,             &
                                   write_datahd_file,            &
                                   DATAHD_LEN,                   &
                                   DATAHD_NUM_NESTS

  use utilities_mod,        only : E_ERR,                        &
                                   error_handler,                &
                                   file_exist,                   &
                                   get_unit

  use time_manager_mod,     only : get_date, set_date, set_time, &
                                   operator(-),                  &
                                   GREGORIAN,                    &
                                   set_calendar_type

  use types_mod,            only : r8

  implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

  character(len=*), parameter :: routine = 'update_datahd_file'

  integer :: nn, np
  integer :: ccyy, mm, dd, hh, mn, ss
  integer :: nml_unit, io_status

  real(kind=r8), dimension(DATAHD_LEN) :: datahd

  character(len=17), parameter :: nml_file = 'update_datahd.nml'
  character(len=10)            :: dtg_new
  character(len=10)            :: dtg_old
  character(len=80)            :: dsnrff='./'

  integer, dimension(3,DATAHD_NUM_NESTS) :: idelay
  character(len=10)                      :: cdtg
  integer                                :: icycle=6
  namelist /update_datahd/ idelay, cdtg, icycle

 call set_calendar_type(GREGORIAN)

  idelay(:,:) = 0

  nml_unit = get_unit()
  
  print *,file_exist(nml_file)
  if (file_exist(nml_file)) then
    open(nml_unit, file=nml_file, status='old', iostat=io_status)
    call check_io_status(io_status, routine, source, revision, revdate,   &
                         'Opening ' // nml_file)

      read(nml_unit, nml=update_datahd, iostat=io_status)
      call check_io_status(io_status, routine, source, revision, revdate, &
                          'Reading ' // nml_file // ' compute_mean')
    close(nml_unit)
  else
    call error_handler(E_ERR, routine, 'namelist read failed', source, revision, revdate)
  end if

  dtg_new = cdtg

  read(dtg_new,'((I4.4),3(I2.2))') ccyy,mm,dd,hh
  call get_date(set_date(ccyy, mm, dd, hh)-set_time(3600*icycle), ccyy, mm, dd, hh, mn, ss)
  write(dtg_old,'((I4.4),3(I2.2))') ccyy,mm,dd,hh

  print *,dtg_new
  print *,dtg_old

  call read_datahd_file(dtg_old, datahd)

  print *,'nnest: ',datahd(11)
  do nn=1,DATAHD_NUM_NESTS
    np=30+(nn-1)*30
    if(any(idelay(:,nn) > 0)) datahd(np+25) = 1
    print *,np,datahd(np+25)
  end do

  !call write_datahd_file(dtg_new, datahd)

end program update_datahd_file

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
