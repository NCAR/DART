! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

PROGRAM driver

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! 1d WRF incorporates original WRF PBL routines for idealized runs.
! It can be initialized with idealized profiles of wind, temperature,
! and mixing ratio or with profiles derived from 3d WRF runs.
! All input is driven by parameters specified in wrf1d_namelist.input
! file. This code is compatible with DART software used for Data
! Assimilatiuon studies.

! The code currently has a problem when compiled with ifc
! as it crashes when calculating MINVAL of an array for unknown reason.
! It runs OK when compiled with pgf90


! The code is available as is but had been checked quite thoroughlly.
! It is of course free to use and modify but I would appreciate
! if you let me know your interest, acknowledge the authorship,
! and, even better if you contribute
! to the development and make the modifications available.

! File README gives more details on input parameters in the namelist.
! File PBLcolumn.pdf gives some general information on the model
! and its coupling with DART.

! 1d WRF was originally developed by Mariusz Pagowski 
! NOAA/GSD/CIRA, Boulder, CO, Mariusz.Pagowski@noaa.gov
! Josh Hacker,NCAR/RAL, hacker@ucar.edu
! contributed reading routines to real data and 3dWRF output and 
! interface to DART and made the code more flexible

! Based on WRFV2.1.2

  USE module_wrf
  USE utilities_mod,    only : initialize_utilities 
  USE time_manager_mod, only : time_type, set_time, get_time, &
                             increment_time, print_time, set_date, &
                             get_date, print_date,&
                             set_calendar_type, GREGORIAN, julian_day, &
                             operator(==), operator(<=), &
                             operator(-), operator(+)

  IMPLICIT NONE

  CHARACTER(len=120) :: &
       namelistfile='wrf1d_namelist.input',&
       outlogfile='wrf1d_log.out',         &
       driverfile='driver_log.out' ! gets its own unit
  INTEGER  :: unit_nml=151,unit_log=152
  LOGICAL :: is_it_there = .FALSE.
  INTEGER :: dart_days, dart_seconds
  INTEGER :: iyr,imo,idy,ihr,imm,iss,seconds_in_day

  INTEGER :: wrf_rnd_seed,itime ! equivalent to itimestep in module_wrf.F

  LOGICAL :: allocate_wrf = .TRUE.
  TYPE(time_type) :: initialization_time, time_step, real_time, &
                     time_into_forecast
  integer, parameter :: calendar_type = GREGORIAN

  CALL initialize_utilities('driver',driverfile)

  INQUIRE ( FILE = namelistfile , EXIST = is_it_there )

  IF ( is_it_there ) THEN
     OPEN ( FILE   = namelistfile     , UNIT   =  unit_nml        ,&
          & STATUS = 'OLD'            , FORM   = 'FORMATTED'      ,&
          & ACTION = 'READ'           , ACCESS = 'SEQUENTIAL'     )
  ELSE
     PRINT '(A)','Could not find the namelist: ',namelistfile
     STOP 'No_namelist_found'
  ENDIF

  OPEN ( FILE   = outlogfile        , UNIT   =  unit_log        ,&
       & STATUS = 'REPLACE'         , FORM   = 'FORMATTED'      ,&
       & ACTION = 'WRITE'           , ACCESS = 'SEQUENTIAL'     )

  CALL do_namelist_wrf1d(unit_nml,unit_log)

  CLOSE(unit_log)
  
  CALL static_init_wrf(allocate_wrf)
  
  wrf_rnd_seed = rnd_seed_val
  allocate_wrf=.FALSE.

! initialize some timing stuff
  time_step = set_time(int(dt), 0)
  call set_calendar_type(calendar_type)
  call init_time(initialization_time)
  call print_time(time_step)

  CALL init_wrf(wrf_rnd_seed)
  

! now this is done in init_wrf
!  IF (init_f) THEN 
!     OPEN(ncunit,file=out_f_file)
!  ENDIF

  real_time = initialization_time
  DO itime=1,ntime
     CALL output_wrf_profiles()
!     IF (init_f) THEN
        time_into_forecast = real_time - initialization_time
        call get_time(time_into_forecast,dart_seconds,dart_days)
!        dart_days=INT(REAL(REAL(itime-1)*dt+start_seconds)/86400.)
!        dart_seconds=NINT(REAL(itime-1)*dt-REAL(dart_days*86400))+&
!             &start_seconds
        call get_date(real_time,iyr,imo,idy,ihr,imm,iss)
        real_time = increment_time(real_time,int(dt),0)

!     ELSE
!        dart_days=0
!        dart_seconds=(itime-1)*dt
!     ENDIF
     IF ( forecast_length == 0 ) stop '0 time steps'
!     IF ( init_f ) THEN
       seconds_in_day = iss + 60*imm + 3600 * ihr
       CALL wrf(dart_seconds,dart_days,julian_day(iyr,imo,idy),seconds_in_day)
!     ELSE
!       iyr = 1970
!       imo = 1
!       idy = 0
!       seconds_in_day = 0
!       CALL wrf(dart_seconds,dart_days,julian_day(iyr,imo,idy),seconds_in_day)
!       CALL wrf(dart_seconds,dart_days)
!     ENDIF
  ENDDO
  
  CONTAINS
subroutine init_time(time)
!------------------------------------------------------------------
!
! Gets the initial time for a state from the model. Where should this info
! come from in the most general case?

type(time_type), intent(out) :: time

integer             :: start_seconds, start_hour, start_month, start_day
character(len=129)  :: errstring

! based on namelist
select case ( init_f_type )
   case ('WRF')
      start_hour = start_hour_f+int(start_forecast/3600)
      if ( start_hour < 24 ) then
      time = set_date(start_year_f, start_month_f,&
                      start_day_f, start_hour_f+int(start_forecast/3600), &
                      start_minute_f,0)
      else
      start_month = start_month_f
      start_day   = start_day_f+1
      start_hour = start_hour_f+int(start_forecast/3600)-24
      if ((start_month == 5 .and. start_day == 32) .or. &
          (start_month == 6 .and. start_day == 31)) then
        start_month = start_month + 1
        start_day = 1
      endif
        
      time = set_date(start_year_f, start_month,&
                      start_day, start_hour, &
                      start_minute_f,0)

      endif
   case ('OBS' , 'SFC' )
      time = set_date(start_year_f, start_month_f,&
                      start_day_f, start_hour_f, &
                      start_minute_f,0)
   case default
      write(*,*) "Do not know how to initialize ",init_f_type
      stop 'init_time'
end select

end subroutine init_time
  
END PROGRAM driver
