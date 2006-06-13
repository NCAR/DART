MODULE module_uvg_force

  USE time_manager_mod,      only: time_type, GREGORIAN, &
                                   set_calendar_type, print_time, &
                                   print_date, set_date, set_time, &
                                   operator(-), operator(+), &
                                   operator(<=), operator(>=), &
                                   operator(/=)
  USE module_model_constants
  USE module_namelist,       only: start_hour_f, &
                                   forecast_length, start_year_f, &
                                   start_month_f, start_day_f, &
                                   start_minute_f, interval_uvg, &
                                   uvg_file, start_forecast, &
                                   rnd_force
  USE module_interpolations, only: linear

  IMPLICIT NONE

  private

  INTEGER, PARAMETER    :: npvars = 4   ! number of physical variables in profs
  INTEGER, PARAMETER    :: nptvars = 2  ! number of timing variables in profs

  INTEGER, PARAMETER    :: calendar_type = GREGORIAN
  INTEGER, DIMENSION(3) :: EPOCH_DATE = (/1970,1,1/)
  INTEGER, DIMENSION(nptvars) :: ptid  ! timing ids
  INTEGER, DIMENSION(npvars)  :: pid  ! variable ids
  INTEGER, DIMENSION(2) :: pdimid, pdimlen, lenp, stp

  REAL,    DIMENSION(npvars)     :: missingVals ! associated netCDF missing flags 
  INTEGER               :: nrecords

  INTEGER, DIMENSION(2)               :: control_index

  REAL                                :: control_w

 
  public                :: uvg_f_dims, uvg_f_bc

CONTAINS

  SUBROUTINE uvg_f_dims(ncid, nz, nt)

   IMPLICIT NONE
   INCLUDE 'netcdf.inc'

! gets all static data from the input file or *_ref, including dimensions

   INTEGER, INTENT(inout):: ncid
   INTEGER, INTENT(out)  :: nz, nt

! local
   INTEGER              :: ierr

! some timing info
   nt = NINT(1+(REAL(forecast_length)) / REAL(interval_uvg))

! open forcing file
   ierr = nf_open(uvg_file, 0, ncid)
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
   ierr=nf_inq_varid(ncid,'u_wind',pid(3))
   ierr=nf_inq_varid(ncid,'v_wind',pid(4))

   ierr=nf_inq_varid(ncid,'base_time',ptid(1))
   ierr=nf_inq_varid(ncid,'time_offset',ptid(2))

! close 
   ierr = nf_close(ncid)

   END SUBROUTINE uvg_f_dims

!**********************************************

  SUBROUTINE uvg_f_bc(ncid, nz, nt, terrain, nz_grid, z_grid, &
                      u_g_f,v_g_f, nsplinetimes, splinetimes, init_f_type, &
                      idum, ims,ime,jms,jme,kms,kme)

    USE module_nr_procedures

    IMPLICIT NONE
    INCLUDE 'netcdf.inc'

! Subroutine gets uvg forcing, and associated p and z profiles, and puts
! them on the column grid

! arguments

    INTEGER, INTENT(inout) :: ncid
    INTEGER, INTENT(in)  :: nz, nt, nz_grid, nsplinetimes
    REAL,    INTENT(in)  :: terrain
    REAL, DIMENSION(:,:,:), INTENT(in):: z_grid
    INTEGER, INTENT(inout) :: idum
    REAL, DIMENSION(:), INTENT(in) :: splinetimes
    REAL, DIMENSION(:,:), INTENT(out) :: u_g_f,v_g_f
    CHARACTER(len=*), INTENT(in) :: init_f_type
    INTEGER, INTENT(in) :: ims,ime,jms,jme,kms,kme
                                      
! local
    INTEGER :: i,k,kkl,kkr,imem, itran
    INTEGER :: ierr, timeid
    INTEGER, DIMENSION(3) :: vec, lenvec

    INTEGER :: nmix, imix
    REAL :: rtran
    
   
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: work
    INTEGER, ALLOCATABLE, DIMENSION(:):: time_snd, time_list, days_list, &
                                         secs_list
    REAL, ALLOCATABLE, DIMENSION(:)   :: times
    REAL, ALLOCATABLE, DIMENSION(:,:) :: z, z_agl, p, u, v, utmp, vtmp

    LOGICAL                           :: got_all_soundings

    INTEGER                                    :: requested_snd_index
    REAL                                       :: missingVal
    TYPE(time_type)                            :: requested_snd_time
    TYPE(time_type), ALLOCATABLE, DIMENSION(:) :: snd_times
    TYPE(time_type)                            :: time_tolerance
    TYPE(time_type)                            :: epoch_time
    TYPE(time_type)                            :: offset_snd

! set up calendar and time tolerance
   call set_calendar_type(calendar_type)
   time_tolerance = set_time(600,0)                   ! 10 mins either way
   epoch_time     = set_date(EPOCH_DATE(1),EPOCH_DATE(2),EPOCH_DATE(3))

   ierr = nf_open(uvg_file, 0, ncid)
   IF ( ierr /= NF_NOERR ) THEN
      PRINT*,"Problem opening forcing file, aborting!"
      STOP
   ENDIF


! get the times
   ALLOCATE(time_list(nrecords), days_list(nrecords), secs_list(nrecords))
   ALLOCATE(snd_times(nrecords), times(nt))

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

! if the year is negative we get a random one valid at the
! same time of day, otherwise the specified sounding
!NOTE - if we are initializing from a sounding, the time will be 30 mins
! off, so account for this.

   offset_snd = set_time(0,0)
   IF ( init_f_type == 'OBS' ) then
      offset_snd = set_time(1800,0) 
   ELSE
      offset_snd = set_time(start_forecast,0)
   ENDIF

   IF ( rnd_force ) THEN
      nmix = 2
   ELSE
      nmix = 1
      requested_snd_index = -9999
      requested_snd_time = set_date(start_year_f, start_month_f, &
                                    start_day_f, start_hour_f,   &
                                    start_minute_f, 0)
      DO i = 1, nrecords
         IF ( snd_times(i)-offset_snd >= (requested_snd_time-time_tolerance) .and. &
              snd_times(i)-offset_snd <= (requested_snd_time+time_tolerance) ) THEN
            requested_snd_index = i 
         ENDIF
      ENDDO
      
      IF ( requested_snd_index < 0 ) THEN
         print*
         print*,'Could not find requested input uvg forcing: '
         call print_date(requested_snd_time)
         stop 'uvg_force_bc'
      ENDIF
      do i = 0, nt-2
        IF (time_list(requested_snd_index+i+1) - &
            time_list(requested_snd_index+i) &
            > interval_uvg + 600 .or. &
            time_list(requested_snd_index+i+1) - &
            time_list(requested_snd_index+i) &
            < interval_uvg - 600 ) THEN
         print*
         print*,'The interval in the sounding file does not agree with ',&
                'your namelist interval_uvg: ',interval_uvg
         print*,time_list(requested_snd_index+1)-time_list(requested_snd_index)
         print*,'This could be caused by either a missing sounding or the ', &
                'an incorrect namelist value'
         stop 'module_snd_init_and_bc'
         ENDIF
      enddo
      IF ( time_list(requested_snd_index) == missingVal ) THEN
         print*
         print*,'Requested input sounding is missing in the input file: '
         call print_date(requested_snd_time)
         stop 'module_snd_init_and_bc'
      ENDIF
      IF ( minval(time_list(requested_snd_index:requested_snd_index+nt-1)) &
                  < 0 ) THEN
         print*
         print*,'A missing sounding prevents a simulation of ',&
                 forecast_length,' seconds'
         stop 'module_snd_init_and_bc'
      ENDIF

      print*,'Getting requested uvg forcing valid at date: '
      call print_date(requested_snd_time+offset_snd)

   ENDIF ! want a specific sounding

! allocations of work arrays
   ALLOCATE(work(npvars,nmix,nz,nt))
   ALLOCATE(time_snd(nt))
   ALLOCATE(z(nz,nt), z_agl(nz,nt), u(nz,nt), v(nz,nt), p(nz,nt))

   lenp = (/nz,nt/)
   got_all_soundings = .FALSE.
   DO WHILE ( .not. got_all_soundings ) 

      got_all_soundings = .true.  ! default is good
      control_index = -9999
      DO imix = 1, nmix

! choose a random record and get the data
! NEED A WAY TO GET THE SAME DATE HERE FOR RND
         IF ( nmix == 1 ) THEN
            itran = requested_snd_index
         ELSE
            rtran = ran1(idum)*(nrecords-nt) + 1.
            itran = AINT(rtran)
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

      END DO ! imix

      ! for now, be satisfied only if we have no missing values
      DO i = 1, npvars
         IF ( minval(work(i,1:nmix,:,:)) == missingVals(i) ) THEN
            got_all_soundings = .false.
         ENDIF
      ENDDO
      IF ( .not. got_all_soundings .and. nmix == 1 ) THEN 
         ! fail in this case because it is specified
         print*
         print*,'Missing data prevents a simulation of ',&
                 forecast_length,' seconds'
         stop 'uvg_f_bc'
      ELSEIF ( .not. got_all_soundings ) THEN
         print*,'Trying again to get full forcing profiles'
      ENDIF
      
   END DO  ! still searching for soundings
   PRINT*,"Getting uvg forcing from times ",control_index
   PRINT*,"with alpha ", control_w

! do grid differencing and weighing of profiles

   SELECT CASE ( nmix )
      CASE (1)
         z(:,1:nt) = work(1,1,:,:)
         p(:,1:nt) = work(2,1,:,:)
         u(:,1:nt) = work(3,1,:,:)
         v(:,1:nt) = work(4,1,:,:)

      CASE (2)
         z(:,1:nt) = work(1,1,:,:)*control_w + work(1,2,:,:)*(1.0-control_w)
         p(:,1:nt) = work(2,1,:,:)*control_w + work(2,2,:,:)*(1.0-control_w)
         u(:,1:nt) = work(3,1,:,:)*control_w + work(3,2,:,:)*(1.0-control_w)
         v(:,1:nt) = work(4,1,:,:)*control_w + work(4,2,:,:)*(1.0-control_w)
         
      CASE DEFAULT
   END SELECT

   ierr = nf_close(ncid)

! move z to z_agl 
   DO i = 1, nt
     z_agl(:,i) = z(:,i) - terrain
   ENDDO

   DEALLOCATE(work)
   DEALLOCATE(time_snd,time_list,days_list,secs_list)

! finally, forcing times  --- these are matched to the profile!

   DO i=1,nt
      times(i)=(i-1)*interval_uvg+start_forecast
   ENDDO

! just go ahead and interpolate to model grid and spline times here
   ALLOCATE(utmp(nz_grid,nt),vtmp(nz_grid,nt))

   DO i = 1,nt
      DO k = 1, nz_grid
         utmp(k,i) = linear(z_agl(:,i),u(:,i),nz,z_grid(1,k,1))
         vtmp(k,i) = linear(z_agl(:,i),v(:,i),nz,z_grid(1,k,1))
      ENDDO
   ENDDO

   DO k = 1, nz_grid
      DO i = 1, nsplinetimes
        u_g_f(k,i) = linear(times,utmp(k,:),nt,splinetimes(i))
        v_g_f(k,i) = linear(times,vtmp(k,:),nt,splinetimes(i))
      ENDDO
   ENDDO

   DEALLOCATE(u,v,p,z,utmp,vtmp)

   RETURN

  END SUBROUTINE uvg_f_bc

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

! p, e should be in the same units.  output in g/g
  IMPLICIT NONE

  REAL, intent(in)     :: p,e

  get_qv = ( r_d / r_v ) * e / (p - e)

  END FUNCTION get_qv

END MODULE module_uvg_force

