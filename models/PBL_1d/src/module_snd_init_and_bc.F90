MODULE module_snd_init_and_bc
!
! DART $Id$
!
  USE types_mod,             only: r8
  USE time_manager_mod,      only: time_type, GREGORIAN, &
                                   set_calendar_type, print_time, &
                                   print_date, set_date, set_time, &
                                   operator(-), operator(+), &
                                   operator(<=), operator(>=), &
                                   operator(/=)
  USE module_model_constants
  USE module_namelist,       only: start_hour_f, interval_f, &
                                   forecast_length, start_year_f, &
                                   start_month_f, start_day_f, &
                                   start_minute_f, interval_flux, &
                                   interval_soil, interval_smos, &
                                   init_f_file, isltyp_ref, &
                                   ivgtyp_ref, lu_index_ref, vegfra_ref, &
                                   lat_ref, lon_ref, julday_ref, mminlu_ref,&
                                   init_flux_file, init_soil_file, &
                                   rnd_init, start_forecast, init_smos_file, &
                                   force_uvg, init_f_type
  USE module_interpolations, only: linear

  IMPLICIT NONE

  private

  INTEGER, PARAMETER    :: npvars = 6   ! number of physical variables in profs
  INTEGER, PARAMETER    :: nptvars = 2  ! number of timing variables in profs
  INTEGER, PARAMETER    :: nfluxvars = 3! number of fluxes (SW, LW)
  INTEGER, PARAMETER    :: nftvars = 1  ! number of flux timing vars
  INTEGER, PARAMETER    :: nsoilvars = 5! number of soil vars &
                                        ! (depth, tslb,smois) W and E
  INTEGER, PARAMETER    :: nsltvars = 1 ! number of soil timing vars
  INTEGER, PARAMETER    :: nsmosvars = 6! PRECIP, SPD, DIR, T, P, e
  INTEGER, PARAMETER    :: nsmtvars = 1 ! smos timing variables

  INTEGER, PARAMETER    :: calendar_type = GREGORIAN
  INTEGER, DIMENSION(3) :: EPOCH_DATE = (/1970,1,1/)
  INTEGER, DIMENSION(nptvars) :: ptid  ! timing ids
  INTEGER, DIMENSION(npvars)  :: pid  ! variable ids
  INTEGER, DIMENSION(2) :: pdimid, pdimlen, lenp, stp

  INTEGER, DIMENSION(nfluxvars) :: flid  ! flux variable ids
  INTEGER, DIMENSION(nftvars)  :: fltid  ! flux timing ids
  INTEGER, DIMENSION(1) :: fldimid, fldimlen, lenfl, stfl

  INTEGER, DIMENSION(nsoilvars) :: slid  ! soil variable ids
  INTEGER, DIMENSION(nsltvars)  :: sltid  ! soil timing ids
  INTEGER, DIMENSION(2) :: sldimid, sldimlen, lensl, stsl

  INTEGER, DIMENSION(nsmosvars) :: smid  ! smos variable ids
  INTEGER, DIMENSION(nsmtvars)  :: smtid  ! smos timing ids
  INTEGER, DIMENSION(1) :: smdimid, smdimlen, lensm, stsm

  INTEGER, DIMENSION(2) :: stter, lenter ! terrain

  REAL,    DIMENSION(npvars)     :: missingVals ! associated netCDF missing flags 
  REAL,    DIMENSION(nfluxvars)  :: missingFluxVals ! associated netCDF missing flags 
  REAL,    DIMENSION(nsoilvars)  :: missingSoilVals ! associated netCDF missing flags 
  REAL,    DIMENSION(nsmosvars)  :: missingSMOSVals ! associated netCDF missing flags 
  INTEGER               :: nrecords, nflux_file, nsoil_file, nsmos_file


  INTEGER, DIMENSION(2)               :: control_index

  REAL                                :: control_w

 
  public                :: snd_f_dims, snd_init_and_bc

CONTAINS

  SUBROUTINE snd_f_dims(ncid, ncid_soil, ncid_flux, ncid_smos, &
                        nz, nt, nt_flux, nt_smos, &
                        ns, nt_soil, &
                        mminlu,julday,&
                        lat, lon, cor, terrain)
   IMPLICIT NONE
   INCLUDE 'netcdf.inc'

! gets all static data from the input file or *_ref, including dimensions

   INTEGER, INTENT(inout):: ncid, ncid_flux, ncid_soil, ncid_smos
   INTEGER, INTENT(out)  :: nz, ns, nt, nt_flux, nt_soil, nt_smos
   INTEGER, INTENT(out)  :: julday
   REAL, INTENT(out)     :: lat,lon,cor, terrain
   CHARACTER(len=4), INTENT(out) :: mminlu

! local
   INTEGER              :: ierr

! some error checking on the namelist choices
   IF ( .not. force_uvg .and. trim(init_f_type) .ne. 'SFC' ) THEN
   IF (mod(forecast_length,interval_f) /= 0) THEN
      print*,'Please choose forecast_length and hour_inteval such that'
      print*,'a profile is available at both init and end times'
      stop 'snd_init_and_bc'
   ENDIF 
   ENDIF

! some timing info
   nt = NINT(1+(REAL(forecast_length)) / REAL(interval_f))
   nt_flux = NINT(1+(REAL(forecast_length)) / REAL(interval_flux))
   nt_soil = NINT(1+(REAL(forecast_length)) / REAL(interval_soil))
   nt_smos = NINT(1+(REAL(forecast_length)) / REAL(interval_smos))

! open forcing file
   ierr = nf_open(init_f_file, 0, ncid)
   IF ( ierr /= NF_NOERR ) THEN
      PRINT*,"Problem opening forcing file, aborting!"
      STOP
   ENDIF

   ierr=nf_inq_dimid(ncid,'time',pdimid(1))
   ierr=nf_inq_dimlen(ncid,pdimid(1),pdimlen(1))
   ierr=nf_inq_dimid(ncid,'record',pdimid(2))
   ierr=nf_inq_dimlen(ncid,pdimid(2),pdimlen(2))
   nz = pdimlen(1)
   nrecords = pdimlen(2)

! variables
   ierr=nf_inq_varid(ncid,'alt',pid(1))
   ierr=nf_inq_varid(ncid,'pres',pid(2))
   ierr=nf_inq_varid(ncid,'tdry',pid(3))
   ierr=nf_inq_varid(ncid,'qv',pid(4))
   ierr=nf_inq_varid(ncid,'u_wind',pid(5))
   ierr=nf_inq_varid(ncid,'v_wind',pid(6))

   ierr=nf_inq_varid(ncid,'base_time',ptid(1))
   ierr=nf_inq_varid(ncid,'time_offset',ptid(2))

! terrain from first altitude  

   stter = (/1,1/)
   lenter = (/1,1/)
   ierr=nf_get_vara_double(ncid,pid(1),stter,lenter,terrain)

! lat and lon from namelist (probably want to make this part of the sounding
! file at some point
   lat=lat_ref
   lon=lon_ref
   cor=2.*eomeg*SIN(lat)

!   ierr = nf_close(ncid)

! fill for now
   julday = julday_ref
   mminlu = mminlu_ref

!-----------------------------------------------------
! radiative fluxes
!-----------------------------------------------------

! open forcing file
   ierr = nf_open(init_flux_file, 0, ncid_flux)
   IF ( ierr /= NF_NOERR ) THEN
      PRINT*,"Problem opening flux file, aborting!"
      STOP
   ENDIF

   ierr=nf_inq_dimid(ncid_flux,'time',fldimid(1))
   ierr=nf_inq_dimlen(ncid_flux,fldimid(1),fldimlen(1))
   nflux_file = fldimlen(1)

! variable ids
   ierr=nf_inq_varid(ncid_flux,'down_short_hemisp',flid(1))
   ierr=nf_inq_varid(ncid_flux,'down_long_hemisp_shaded',flid(2))
   ierr=nf_inq_varid(ncid_flux,'up_long_hemisp',flid(3))

   ierr=nf_inq_varid(ncid_flux,'time',fltid(1))

!-----------------------------------------------------
! soil
!-----------------------------------------------------

! open forcing file
   ierr = nf_open(init_soil_file, 0, ncid_soil)
   IF ( ierr /= NF_NOERR ) THEN
      PRINT*,"Problem opening soil file, aborting!"
      STOP
   ENDIF

   ierr=nf_inq_dimid(ncid_soil,'time',sldimid(1))
   ierr=nf_inq_dimlen(ncid_soil,sldimid(1),sldimlen(1))
   nsoil_file = sldimlen(1)

   ierr=nf_inq_dimid(ncid_soil,'depth',sldimid(2))
   ierr=nf_inq_dimlen(ncid_soil,sldimid(2),sldimlen(2))
   ns = sldimlen(2)
 
! variable ids

   ierr=nf_inq_varid(ncid_soil,'depth',slid(1))
   ierr=nf_inq_varid(ncid_soil,'tsoil_W',slid(2))
   ierr=nf_inq_varid(ncid_soil,'tsoil_E',slid(3))
   ierr=nf_inq_varid(ncid_soil,'watcont_W',slid(4))
   ierr=nf_inq_varid(ncid_soil,'watcont_E',slid(5))

   ierr=nf_inq_varid(ncid_soil,'time',sltid(1))

!-----------------------------------------------------
! surface meteorology (smos), including precip and more if necessary
!-----------------------------------------------------

! open forcing file
   ierr = nf_open(init_smos_file, 0, ncid_smos)
   IF ( ierr /= NF_NOERR ) THEN
      PRINT*,"Problem opening smos file, aborting!"
      STOP
   ENDIF

   ierr=nf_inq_dimid(ncid_smos,'time',smdimid(1))
   ierr=nf_inq_dimlen(ncid_smos,smdimid(1),smdimlen(1))
   nsmos_file = smdimlen(1)

! variable ids - could add screen values here

   ierr=nf_inq_varid(ncid_smos,'precip',smid(1))
   ierr=nf_inq_varid(ncid_smos,'wspd',smid(2))
   ierr=nf_inq_varid(ncid_smos,'wdir',smid(3))
   ierr=nf_inq_varid(ncid_smos,'temp',smid(4))
   ierr=nf_inq_varid(ncid_smos,'bar_pres',smid(5))
   ierr=nf_inq_varid(ncid_smos,'vap_pres',smid(6))

   ierr=nf_inq_varid(ncid_smos,'time',smtid(1))

   ierr=nf_close(ncid)
   ierr=nf_close(ncid_soil)
   ierr=nf_close(ncid_flux)
   ierr=nf_close(ncid_smos)

   END SUBROUTINE snd_f_dims

!**********************************************

  SUBROUTINE snd_init_and_bc(ncid, ncid_flux, ncid_soil, ncid_smos, &
       nz, ns, nt, nt_flux, nt_soil, nt_smos, &
       z_agl,t,u,v,q,p,&
       th2,t2, tsk,&
       u10, v10, q2, precip, &
       glw,glw_up,gsw,qsfc,z_s,tslb,smois,tmn,&
       vegfra,isltyp,lu_index,ivgtyp,&
       times,times_flux,times_soil,times_smos,idum)

    USE module_nr_procedures

    IMPLICIT NONE
    INCLUDE 'netcdf.inc'

! Subroutine selects a random forecast from observation input files
! For profiles and fluxes, there must not be any missing through the length
! of the simulation.  For soil, only need initialization.

! arguments

    INTEGER, INTENT(inout) :: ncid, ncid_flux, ncid_soil, ncid_smos 
    INTEGER, INTENT(in)  :: nz, ns, nt, nt_flux, nt_soil, nt_smos
    INTEGER, INTENT(inout) :: idum
    REAL, DIMENSION(:), INTENT(out) :: th2,t2,u10,v10,q2, &
                                       precip, &
                                       tsk,glw,glw_up,gsw,qsfc,&
                                       tmn,vegfra,times,times_flux, &
                                       times_soil, times_smos
    INTEGER, DIMENSION(:), INTENT(out) :: isltyp,ivgtyp,lu_index

    REAL, DIMENSION(:),   INTENT(out) :: z_s
    REAL, DIMENSION(:,:), INTENT(out) :: t,u,v,q,p,z_agl,tslb,smois
    
! local
    INTEGER :: i,k,kkl,kkr,imem, itran, itranflux,itransoil, itransmos
    INTEGER :: ierr, timeid
    INTEGER, DIMENSION(3) :: vec, lenvec

    INTEGER :: nmix, imix
    REAL :: rtran
    
   
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: work, work_soil
    REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: work_flux, work_smos
    INTEGER, ALLOCATABLE, DIMENSION(:):: time_snd, time_list, days_list, &
                                         secs_list, secs_smos,days_smos
    INTEGER                           :: timediff
    REAL, ALLOCATABLE, DIMENSION(:)   :: time_flux, time_soil, time_smos
    REAL, ALLOCATABLE, DIMENSION(:,:) :: offset_snd
    REAL, ALLOCATABLE, DIMENSION(:,:) :: invrho,p_r,p_l,p_u,p_d
    REAL, ALLOCATABLE, DIMENSION(:,:) :: z

    LOGICAL                           :: got_all_soundings 

    INTEGER                                    :: requested_snd_index,&
                                                  requested_flux_index, &
                                                  requested_soil_index, &
                                                  requested_smos_index, &
                                                  req_flux_secs, &
                                                  req_flux_days, &
                                                  req_soil_secs, &
                                                  req_soil_days, &
                                                  req_smos_secs, &
                                                  req_smos_days
    REAL                                       :: missingVal
    TYPE(time_type)                            :: requested_snd_time
    TYPE(time_type), ALLOCATABLE, DIMENSION(:) :: snd_times, smos_times
    TYPE(time_type)                            :: time_tolerance
    TYPE(time_type)                            :: epoch_time
    TYPE(time_type)                            :: req_flux_time
    TYPE(time_type)                            :: req_soil_time
    TYPE(time_type)                            :: req_smos_time

! set up calendar and time tolerance
   call set_calendar_type(calendar_type)
   time_tolerance = set_time(600,0)                   ! 10 mins either way
   epoch_time     = set_date(EPOCH_DATE(1),EPOCH_DATE(2),EPOCH_DATE(3))

! open forcing file
   ierr = nf_open(init_f_file, 0, ncid)
   IF ( ierr /= NF_NOERR ) THEN
      PRINT*,"Problem opening forcing file, aborting!"
      STOP
   ENDIF
! open soil file
   ierr = nf_open(init_soil_file, 0, ncid_soil)
   IF ( ierr /= NF_NOERR ) THEN
      PRINT*,"Problem opening soil file, aborting!"
      STOP
   ENDIF
! open flux file
   ierr = nf_open(init_flux_file, 0, ncid_flux)
   IF ( ierr /= NF_NOERR ) THEN
      PRINT*,"Problem opening flux file, aborting!"
      STOP
   ENDIF
! open smos file
   ierr = nf_open(init_smos_file, 0, ncid_smos)
   IF ( ierr /= NF_NOERR ) THEN
      PRINT*,"Problem opening smos file, aborting!"
      STOP
   ENDIF



! get the times
   ALLOCATE(time_list(nrecords), days_list(nrecords), secs_list(nrecords))
   ALLOCATE(snd_times(nrecords), smos_times(nsmos_file))
   ALLOCATE(time_flux(nflux_file), time_soil(nsoil_file))
   ALLOCATE(time_smos(nsmos_file), days_smos(nsmos_file), secs_smos(nsmos_file))

   ierr= nf_get_vara_int(ncid,ptid(1),(/1/),(/nrecords/),time_list(:))
   ierr= nf_get_att_double(ncid, ptid(2), "missing_value", missingVal)

   WHERE (time_list < 0)
     time_list = missingVal
   ELSE WHERE
     secs_list = mod(time_list,86400)
     days_list = (time_list - secs_list) / 86400 
   END WHERE

   DO i = 1, nrecords
     IF (time_list(i) /= missingVal) snd_times(i) = epoch_time + &
                                     set_time(secs_list(i),days_list(i))
   ENDDO

   ierr= nf_get_vara_double(ncid_smos,smtid(1),(/1/),(/nsmos_file/),time_smos)

   WHERE (time_smos < 0)
     time_smos = missingVal
   ELSE WHERE
     secs_smos = mod(int(time_smos),86400)
     days_smos = (int(time_smos) - secs_smos) / 86400 
   END WHERE

   DO i = 1, nsmos_file
     IF (int(time_smos(i)) /= missingVal) smos_times(i) = epoch_time + &
                                     set_time(secs_smos(i),days_smos(i))
   ENDDO

   ierr= nf_get_vara_double(ncid_flux,fltid(1),(/1/),(/nflux_file/),time_flux)
   ierr= nf_get_vara_double(ncid_soil,sltid(1),(/1/),(/nsoil_file/),time_soil)

! if the year is negative we get a random one valid at the
! same time of day, otherwise the specified sounding

   IF ( rnd_init /= 2 ) THEN
print*, 'CANNOT CHOOSE RANDOM OBS YET!'
stop 'snd_init_and_bc'
      nmix = 2
   ELSE
      nmix = 1
      requested_snd_index = -9999
      requested_flux_index = -9999
      requested_soil_index = -9999
      requested_smos_index = -9999
      requested_snd_time = set_date(start_year_f, start_month_f, &
                                    start_day_f, start_hour_f,   &
                                    start_minute_f, 0)

      !------------------------------------------------------
      IF ( trim(init_f_type) /= 'SFC' ) THEN ! key off sounding

      DO i = 1, nrecords
         IF ( snd_times(i) >= (requested_snd_time-time_tolerance) .and. &
              snd_times(i) <= (requested_snd_time+time_tolerance) ) THEN
            requested_snd_index = i 
         ENDIF
      ENDDO
      
      IF ( requested_snd_index < 0 ) THEN
         print*
         print*,'Could not find requested input sounding: '
         call print_date(requested_snd_time)
         if ( trim(init_f_type) /= 'SFC' ) stop 'module_snd_init_and_bc'
      ENDIF
      IF ( .not. force_uvg ) THEN
      do i = 0, nt-2
        timediff = time_list(requested_snd_index+i+1) - &
            time_list(requested_snd_index+i)
        IF (timediff > interval_f + 600) THEN
         print*
         print*,'The interval in the sounding file does not agree with ',&
                'your namelist interval_f: ',interval_f
         print*,timediff
         print*,'This could be caused by either a missing sounding or ', &
                'an incorrect namelist value'
         if ( trim(init_f_type) /= 'SFC' ) stop 'module_snd_init_and_bc'
         ENDIF
      enddo
      ENDIF ! init_type = SFC

      IF ( time_list(requested_snd_index) == missingVal ) THEN
         print*
         print*,'Requested input sounding is missing in the input file: '
         call print_date(requested_snd_time)
         if ( trim(init_f_type) /= 'SFC' ) stop 'module_snd_init_and_bc'
      ENDIF
      IF ( .not. force_uvg ) THEN
      IF ( minval(time_list(requested_snd_index:requested_snd_index+nt-1)) &
                  < 0 ) THEN
         print*
         print*,'A missing sounding prevents a simulation of ',&
                 forecast_length,' seconds'
         if ( trim(init_f_type) /= 'SFC' ) stop 'module_snd_init_and_bc'
      ENDIF
      ENDIF

      print*,'Getting requested sounding valid at date: '
      call print_date(requested_snd_time)
      print*
      print*,'If you are running DART, the correct time ',&
             ' for the obs_sequence is near'
      call print_time(requested_snd_time)

      call find_matching_index(int(time_smos),time_list(requested_snd_index), &
                           requested_smos_index, &
                           req_smos_secs, req_smos_days)

      IF ( requested_smos_index < 0 ) THEN
         print*
         print*,'Could not find requested input surface met: '
         CALL print_date(requested_snd_time)
         stop 'module_snd_init_and_bc'
      ENDIF

      !------------------------------------------------------
      ELSE ! key off smos (only care about sfc obs)

      DO i = 1, nsmos_file
         IF ( smos_times(i) >= (requested_snd_time-time_tolerance) .and. &
              smos_times(i) <= (requested_snd_time+time_tolerance) ) THEN
            requested_smos_index = i 
         ENDIF
      ENDDO
      IF ( requested_smos_index < 0 ) THEN
         print*
         print*,'Could not find requested input smos data: '
         call print_date(requested_snd_time)
         stop 'module_snd_init_and_bc'
      ENDIF

      req_smos_secs = secs_smos(requested_smos_index)
      req_smos_days = days_smos(requested_smos_index)

      ENDIF
!  if we are running to get the obs out only (init_f_type == 'SFC') then
!  use the smos times directly and ignore everything else
   

      IF (time_smos(requested_smos_index+1) - time_smos(requested_smos_index) &
          /= interval_smos ) THEN
         print*
         print*,'The interval in the smos file does not agree with ',&
                'your namelist interval_f_soil: ',interval_smos
         print*,time_smos(requested_smos_index+1)-time_smos(requested_smos_index)
         stop 'module_snd_init_and_bc'
      ENDIF
 
      req_smos_time = epoch_time + set_time(req_smos_secs,req_smos_days)
      print*,'Using corresponding SMOS date: '
      call print_date(req_smos_time)

   ! and the specific flux to go with it
      IF ( trim(init_f_type)  /= 'SFC' ) THEN

      CALL find_matching_index(INT(time_flux),time_list(requested_snd_index), &
                           requested_flux_index, &
                           req_flux_secs, req_flux_days)
 
      ELSE

      CALL find_matching_index(INT(time_flux), &
                           INT(time_smos(requested_smos_index)), &
                           requested_flux_index, &
                           req_flux_secs, req_flux_days)

      ENDIF

      IF ( requested_flux_index < 0 ) THEN
         print*
         print*,'Could not find requested input fluxes: '
         call print_date(requested_snd_time)
         stop 'module_snd_init_and_bc'
      ENDIF
      IF (time_flux(requested_flux_index+1) - time_flux(requested_flux_index) &
          /= interval_flux ) THEN
         print*
         print*,'The interval in the flux file does not agree with ',&
                'your namelist interval_f_flux: ',interval_flux
         print*,time_flux(requested_flux_index+1)-time_flux(requested_flux_index)
         stop 'module_snd_init_and_bc'
      ENDIF
 
      req_flux_time = epoch_time + set_time(req_flux_secs,req_flux_days)
      print*,'Using corresponding FLUX date: '
      call print_date(req_flux_time)

   ! and the specific soil info to go with it
      IF ( trim(init_f_type)  /= 'SFC' ) THEN

      CALL find_matching_index(INT(time_soil),time_list(requested_snd_index), &
                           requested_soil_index, &
                           req_soil_secs, req_soil_days)

      ELSE

      CALL find_matching_index(INT(time_soil),  &
                           INT(time_smos(requested_smos_index)), &
                           requested_soil_index, &
                           req_soil_secs, req_soil_days)

      ENDIF

      IF ( requested_soil_index < 0 ) THEN
         print*
         print*,'Could not find requested input soil state: '
         call print_date(requested_snd_time)
         stop 'module_snd_init_and_bc'
      ENDIF
      IF (time_soil(requested_soil_index+1) - time_soil(requested_soil_index) &
          /= interval_soil ) THEN
         print*
         print*,'The interval in the soil file does not agree with ',&
                'your namelist interval_f_soil: ',interval_soil
         print*,time_soil(requested_soil_index+1)-time_soil(requested_soil_index)
         stop 'module_snd_init_and_bc'
      ENDIF
 
      req_soil_time = epoch_time + set_time(req_soil_secs,req_soil_days)
      print*,'Using corresponding SOIL date: '
      call print_date(req_soil_time)

   ENDIF ! want a specific sounding

! init the arrays because they will not be completely filled
   t = -9999.
   u = -9999.
   v = -9999.
   q = -9999.
   p = -9999.
   th2 = -9999.
   t2 = -9999.
   tsk = -9999.
   u10 = -9999.
   v10 = -9999.
   q2 = -9999.
   precip = -9999.
   z_agl = -9999.
   glw = -9999.
   glw_up = -9999.
   gsw = -9999.
   tmn = -9999.
   qsfc = -9999.
   tslb = -9999.
   smois = -9999.
   times = -9999.
   times_flux = -9999.
   times_soil = -9999.
   times_smos = -9999.
   vegfra = -9999.
   lu_index=-9999
   isltyp = -9999
   ivgtyp = -9999

! allocations of work arrays
   ALLOCATE(work(npvars,nmix,nz,nt))
   ALLOCATE(work_flux(nfluxvars,nmix,nt_flux))
   ALLOCATE(work_soil(nsoilvars,nmix,ns,nt_soil))
   ALLOCATE(work_smos(nsmosvars,nmix,nt_smos))
   ALLOCATE(time_snd(nt),offset_snd(nz,nt))
   ALLOCATE(z(nz,nt))

   lenp = (/nz,nt/)
   lenfl = (/nt_flux/)
   lensl = (/ns,nt_soil/)
   lensm = (/nt_smos/)
   got_all_soundings = .FALSE.
print*,'got them all? ',got_all_soundings
   DO WHILE ( .not. got_all_soundings ) 

      got_all_soundings = .true.  ! default is good
      control_index = -9999
      DO imix = 1, nmix

! choose a random record and get the data
         IF ( nmix == 1 ) THEN
            itran = requested_snd_index
            itranflux = requested_flux_index
            itransoil = requested_soil_index
            itransmos = requested_smos_index
         ELSE
            rtran = ran1(idum)*(nrecords-nt) + 1.
            itran = AINT(rtran)
            CALL find_matching_index(INT(time_flux),time_list(itran), &
                           itranflux)
            CALL find_matching_index(INT(time_soil),time_list(itran), &
                           itransoil)

         ENDIF
         control_index(imix) = itran
         rtran = ran1(idum)
         control_w = rtran

         !profiles
         stp = (/1,itran/)

         DO i = 1, npvars
            ierr=nf_get_vara_double(ncid,pid(i),stp,lenp,work(i,imix,:,:))
            ierr=nf_get_att_double(ncid, pid(i), "missing_value", missingVals(i))
         ENDDO

         !fluxes
         stfl = (/itranflux/)

         DO i = 1, nfluxvars
            ierr=nf_get_vara_double(ncid_flux,flid(i),stfl,lenfl, &
                 work_flux(i,imix,:))
            ierr=nf_get_att_double(ncid_flux, flid(i), "missing_value", missingFluxVals(i))
         ENDDO

         !soil
         stsl = (/1,itransoil/)

         ! depth is a special case
         ierr=nf_get_vara_double(ncid_soil,slid(1),(/1/),(/ns/), &
                 work_soil(1,imix,:,1))

         DO i = 2, nsoilvars
            ierr=nf_get_vara_double(ncid_soil,slid(i),stsl,lensl, &
                 work_soil(i,imix,:,:))
            ierr=nf_get_att_double(ncid_soil, slid(i), "missing_value", missingsoilVals(i))
         ENDDO

         !smos
         stsm = (/itransmos/)

         DO i = 1, nsmosvars
            ierr=nf_get_vara_double(ncid_smos,smid(i),stsm,lensm, &
                 work_smos(i,imix,:))
            ierr=nf_get_att_double(ncid_smos, smid(i), "missing_value", missingSMOSVals(i))
         ENDDO

      END DO ! imix

      ! for now, be satisfied only if we have no missing values
      DO i = 1, npvars
         IF ( minval(work(i,1:nmix,:,1)) == missingVals(i) ) THEN
            got_all_soundings = .false.
         ENDIF
      ENDDO
      IF ( .not. force_uvg ) THEN
      DO i = 1, npvars
         IF ( minval(work(i,1:nmix,:,2:)) == missingVals(i) ) THEN
            got_all_soundings = .false.
         ENDIF
      ENDDO
      ENDIF
      DO i = 1, nfluxvars
         IF ( minval(work_flux(i,1:nmix,:)) == missingFluxVals(i) ) THEN
            got_all_soundings = .false.
         ENDIF
      ENDDO
      DO i = 2, nsoilvars ! no depth, and only care about first time and layer here
         IF ( minval(work_soil(i,1:nmix,1,1)) == missingSoilVals(i) ) THEN
            got_all_soundings = .false.
         ENDIF
      ENDDO
      DO i = 1, nsmosvars 
         IF ( minval(work_smos(i,1:nmix,:)) == missingSMOSVals(i) ) THEN
            got_all_soundings = .false.
         ENDIF
      ENDDO
      IF ( .not. got_all_soundings .and. nmix == 1 ) THEN 
         ! fail in this case because it is specified
         print*
         print*,'Missing data prevents a simulation of ',&
                 forecast_length,' seconds'
        
         if ( trim(init_f_type) == 'SFC' ) then
           got_all_soundings = .true.
         else
           stop 'module_snd_init_and_bc'
         endif
      ELSEIF ( .not. got_all_soundings ) THEN
         print*,'Trying again to get full forcing profiles'
      ENDIF
      
   END DO  ! still searching for soundings
   PRINT*,"Getting profile from times ",control_index
   PRINT*,"with alpha ", control_w
   PRINT*,"mixing ", nmix

! do grid differencing and weighing of profiles

   SELECT CASE ( nmix )
      CASE (1)
         z(:,1:nt) = work(1,1,:,:)
         p(:,1:nt) = work(2,1,:,:)
         t(:,1:nt) = work(3,1,:,:)
         q(:,1:nt) = work(4,1,:,:)
         u(:,1:nt) = work(5,1,:,:)
         v(:,1:nt) = work(6,1,:,:)
         gsw(1:nt_flux) = work_flux(1,1,:)
         glw(1:nt_flux) = work_flux(2,1,:)
         glw_up(1:nt_flux) = work_flux(3,1,:)
         tslb(:,1:nt_soil) = 0.5*(work_soil(2,1,:,:)+work_soil(3,1,:,:))
         smois(:,1:nt_soil) = 0.5*(work_soil(4,1,:,:)+work_soil(5,1,:,:))
         precip(1:nt_smos) = work_smos(1,1,:)
         t2(1:nt_smos) = work_smos(4,1,:) + SVPT0 ! C to K
! don't fill th2 for now - creates problems later because the sounding
! is AGL, and the first height is 0.0.  Probably want to replace
         th2(1:nt_smos) = t2(1:nt_smos) * &
                          (p1000mb/(work_smos(5,1,1:nt_smos)*1000.0_r8))**rcp 
         DO i = 1,nt_smos
            call get_uv(work_smos(2,1,i),work_smos(3,1,i),u10(i),v10(i))
            q2(i) = get_qv(work_smos(5,1,i),work_smos(6,1,i))
         ENDDO

      CASE (2)
         CALL blend_snd(z,p,t,q,u,v,work,gsw,glw,glw_up,work_flux,&
                        tslb, smois, work_soil, &
                        precip, work_smos, &
                        control_w,nt,nt_flux,nt_soil,nt_smos,nz,ns)
         
      CASE DEFAULT
   END SELECT

! move z to z_agl and set first level to screen
   DO i = 1, nt
     z_agl(:,i) = z(:,i) - z(1,i)
     z_agl(1,i) = 2.0_r8
     if ( z_agl(2,i) == 2.0_r8 ) z_agl(2,i) = 3.0_r8
   ENDDO

! store soil depth
   z_s = work_soil(1,1,:,1)

   ierr=nf_close(ncid)
   ierr=nf_close(ncid_soil)
   ierr=nf_close(ncid_flux)
   ierr=nf_close(ncid_smos)

   DEALLOCATE(work)
   DEALLOCATE(time_snd,offset_snd,time_list,days_list,secs_list)
   DEALLOCATE(work_flux)
   DEALLOCATE(time_flux)
   DEALLOCATE(work_soil)
   DEALLOCATE(time_soil)
   DEALLOCATE(work_smos)
   DEALLOCATE(time_smos)

! finally, forecast times  --- these are matched to the profile!
! The soil may be offset.  Start_forecast = 0 but put here for consistency

   DO i=1,nt
      times(i)=(i-1)*interval_f+start_forecast 
   ENDDO
   DO i=1,nt_flux
      times_flux(i)=(i-1)*interval_flux+start_forecast 
   ENDDO
   DO i=1,nt_soil
      times_soil(i)=(i-1)*interval_soil+start_forecast 
   ENDDO
   DO i=1,nt_smos
      times_smos(i)=(i-1)*interval_smos+start_forecast 
   ENDDO

!   ierr = nf_close(ncid)
!   ierr = nf_close(ncid_flux)
!   ierr = nf_close(ncid_soil)

! and any unit changes we need to make
   IF ( maxval(p) < 2000.0 ) p = p * 100.0   ! mb to Pa
   IF ( maxval(t) < 200.0 )  t = t + 273.16  ! C to K
   IF ( maxval(tslb) < 200.0 ) THEN
     WHERE (tslb > -500) tslb = tslb + 273.16
   ENDIF

  END SUBROUTINE snd_init_and_bc

!-------------------------------------------------------------------
  SUBROUTINE blend_snd(z,p,t,q,u,v,work,&
                       gsw,glw,glw_up,work_flux, &
                       precip, work_smos, &
                       tslb, smois, work_soil, &
                       w,nt,nt_flux,nt_soil,nt_smos,nz,ns)
  
  IMPLICIT NONE
!NOTE - this needs work to handle soil profiles better (i.e. missing
! should be dependent on levels, etc), and probably the sounding selection 
! is not the way to go.  Maybe should go the climo mean+pert route.
!
! Subroutine finds an average profile of f.  A simple approach:
! 1.  Pick one member based on which has the lowest top height at t=0
! 2.  Interpolate to the heights of those members, and average accordingly

  INTEGER,                  INTENT(IN)         :: nt, nz, ns, &
                                                  nt_flux, nt_soil, nt_smos
  REAL,                     INTENT(IN)         :: w
  REAL, DIMENSION(nz,nt),   INTENT(OUT)        :: z,p,t,q,u,v
  REAL, DIMENSION(nt_flux), INTENT(OUT)        :: gsw,glw,glw_up
  REAL, DIMENSION(ns,nt_soil), INTENT(OUT)     :: tslb,smois
  REAL, DIMENSION(nt_smos), INTENT(OUT)        :: precip
  REAL, DIMENSION(npvars,2,nz,nt), INTENT(IN)  :: work
  REAL, DIMENSION(nfluxvars,2,nt_flux), INTENT(IN)  :: work_flux
  REAL, DIMENSION(nsoilvars,2,ns,nt_soil), INTENT(IN)  :: work_soil
  REAL, DIMENSION(nsmosvars,2,nt_smos), INTENT(IN)  :: work_smos

  INTEGER                                      :: iz, it, asnd, bsnd
  REAL, DIMENSION(nz)                          :: zvals

   DO it = 1, nt

   asnd = 1
   bsnd = 2
   IF ( work(1,1,nz,it) > work(1,2,nz,it) ) THEN
      asnd = 2
      bsnd = 1
   ENDIF

     zvals = work(1,asnd,:,it)

     DO iz = 1, nz

        z(iz,it) = work(1,asnd,iz,it) 
        p(iz,it) = w*work(2,asnd,iz,it) + &
           (1.0-w)*linear(work(1,bsnd,:,it),work(2,bsnd,:,it),nz,z(iz,it)) 
        t(iz,it) = w*work(3,asnd,iz,it) + &
           (1.0-w)*linear(work(1,bsnd,:,it),work(3,bsnd,:,it),nz,z(iz,it))
        q(iz,it) = w*work(4,asnd,iz,it) + &
           (1.0-w)*linear(work(1,bsnd,:,it),work(4,bsnd,:,it),nz,z(iz,it))
        u(iz,it) = w*work(5,asnd,iz,it) + &
           (1.0-w)*linear(work(1,bsnd,:,it),work(5,bsnd,:,it),nz,z(iz,it))
        v(iz,it) = w*work(6,asnd,iz,it) + &
           (1.0-w)*linear(work(1,bsnd,:,it),work(6,bsnd,:,it),nz,z(iz,it))

     ENDDO !iz
   ENDDO !it

   DO it = 1, nt_flux
     gsw(it) = w*work_flux(1,asnd,it) + (1.0-w)*work_flux(1,bsnd,it)
     glw(it) = w*work_flux(2,asnd,it) + (1.0-w)*work_flux(2,bsnd,it)
     glw_up(it) = w*work_flux(3,asnd,it) + (1.0-w)*work_flux(3,bsnd,it)
   ENDDO !it

   DO it = 1, nt_soil
     DO iz = 1, ns
        WHERE ( work_soil(2,asnd,:,it) > -500 .and. &
                work_soil(3,asnd,:,it) > -500 .and. &
                work_soil(2,bsnd,:,it) > -500 .and. &
                work_soil(3,bsnd,:,it) > -500 ) &
        tslb(:,it) = w*0.5*(work_soil(2,asnd,:,it) + &
                            work_soil(3,asnd,:,it))  + &
                     (1.0-w)*0.5*(work_soil(2,bsnd,:,it) + &
                                  work_soil(3,bsnd,:,it))

        WHERE ( work_soil(4,asnd,:,it) > -500 .and. &
                work_soil(5,asnd,:,it) > -500 .and. &
                work_soil(4,bsnd,:,it) > -500 .and. &
                work_soil(5,bsnd,:,it) > -500 ) &
         smois(:,it) = w*0.5*(work_soil(4,asnd,:,it) + &
                             work_soil(5,asnd,:,it))  + &
                     (1.0-w)*0.5*(work_soil(4,bsnd,:,it) + &
                                  work_soil(5,bsnd,:,it))
     ENDDO !iz (soil)

   ENDDO !it

   DO it = 1, nt_smos
     precip(it) = w*work_smos(1,asnd,it) + (1.0-w)*work_smos(1,bsnd,it)
   ENDDO !it
  return

  END SUBROUTINE blend_snd

!--------------------------------------------------------------------
  SUBROUTINE find_matching_index(tlist, t_to_match, &
                           matching_index, &
                           matching_secs, matching_days )

      INTEGER, DIMENSION(:), INTENT(in)           :: tlist
      INTEGER,               INTENT(in)           :: t_to_match
      INTEGER,               INTENT(out)          :: matching_index
      INTEGER, OPTIONAL,     INTENT(out)          :: matching_secs, &
                                                     matching_days

      INTEGER                                     :: minlocs(1)

      matching_index = -9999

      minlocs = minloc((dble(tlist)-dble(t_to_match))**2)

      if ( minval((dble(tlist)-dble(t_to_match))**2) > 3600**2 ) then
         print*,'Matching time is greater than 1h different!'
         stop 'find_matching_index:module_snd_and_bc'
      endif

      matching_index = minlocs(1)

      IF ( matching_index < 0 ) THEN
         return
      ENDIF

      IF (PRESENT(matching_secs)) &
         matching_secs = mod(tlist(matching_index),86400)
      IF (PRESENT(matching_days)) &
         matching_days = (tlist(matching_index) &
                          - matching_secs) / 86400 

  END SUBROUTINE find_matching_index

!------------------------------------------------------------------------
  SUBROUTINE get_uv(spd,dir,u,v)

  IMPLICIT NONE

  REAL, intent(in)     :: spd,dir
  REAL, intent(out)    :: u,v

  u = -spd * dsin(dir*DEGRAD) 
  v = -spd * dcos(dir*DEGRAD) 

  END SUBROUTINE get_uv

!------------------------------------------------------------------------
  REAL FUNCTION get_qv(p,e)

! returns water vapor mixing ratio
! p, e should be in the same units.  output in g/g
  IMPLICIT NONE

  REAL, intent(in)     :: p,e

  get_qv = ( r_d / r_v ) * e / (p - e)

  END FUNCTION get_qv

END MODULE module_snd_init_and_bc

