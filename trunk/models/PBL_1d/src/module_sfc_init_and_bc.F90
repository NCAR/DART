MODULE module_sfc_init_and_bc

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
                                   start_minute_f, interval_sfc,&
                                   sfc_file, &
                                   rnd_init, start_forecast

  USE module_interpolations, only: linear

  IMPLICIT NONE

  private

  INTEGER, PARAMETER    :: nsfcvars = 5   ! number of physical variables in profs
  INTEGER, PARAMETER    :: nftvars = 1  ! number of timing variables in fluxes

  INTEGER, PARAMETER    :: calendar_type = GREGORIAN
  INTEGER, DIMENSION(3) :: EPOCH_DATE = (/1970,1,1/)
  INTEGER, DIMENSION(nsfcvars) :: flid  ! sfc variable ids
  INTEGER, DIMENSION(nftvars)  :: fltid  ! sfc timing ids
  INTEGER, DIMENSION(1) :: fldimid, fldimlen, lenfl, stfl


  REAL,    DIMENSION(nsfcvars)     :: missingSfcVals ! associated netCDF missing flags 

  INTEGER               :: nsfc_file


  INTEGER, DIMENSION(2)               :: control_index

  REAL                                :: control_w

 
  public                :: sfc_f_dims, sfc_init_and_bc

CONTAINS

  SUBROUTINE sfc_f_dims(ncid_sfc, nt_sfc)

   IMPLICIT NONE
   INCLUDE 'netcdf.inc'

! gets all static data from the input file or *_ref, including dimensions

   INTEGER, INTENT(inout):: ncid_sfc
   INTEGER, INTENT(out)  :: nt_sfc

! local
   INTEGER              :: ierr

! some error checking on the namelist choices

! some timing info
   nt_sfc = NINT(1+(REAL(forecast_length)) / REAL(interval_sfc))

!-----------------------------------------------------
! surface data for sfc forcing option
!-----------------------------------------------------

! open forcing file
   ierr = nf_open(sfc_file, 0, ncid_sfc)
   IF ( ierr /= NF_NOERR ) THEN
      PRINT*,"Problem opening sfc file, aborting!"
      STOP
   ENDIF

   ierr=nf_inq_dimid(ncid_sfc,'time',fldimid(1))
   ierr=nf_inq_dimlen(ncid_sfc,fldimid(1),fldimlen(1))
   nsfc_file = fldimlen(1)

! variable ids
   ierr=nf_inq_varid(ncid_sfc,'ts',flid(1))
   ierr=nf_inq_varid(ncid_sfc,'qvs',flid(2))
   ierr=nf_inq_varid(ncid_sfc,'ustar',flid(3))
   ierr=nf_inq_varid(ncid_sfc,'hflux',flid(4))
   ierr=nf_inq_varid(ncid_sfc,'qvflux',flid(5))

   ierr=nf_inq_varid(ncid_sfc,'time',fltid(1))

   END SUBROUTINE sfc_f_dims

!**********************************************

   SUBROUTINE sfc_init_and_bc(ncid_sfc, nt_sfc,&
        &ts,qvs,ustar,hflux,qvflux,times_sfc,idum)

    USE module_nr_procedures

    IMPLICIT NONE
    INCLUDE 'netcdf.inc'

! Subroutine selects a random forecast from observation input files
! there must not be any missing sfcs,ts,qvs (there are not)
! need input for all times when forecast to be issued


    INTEGER, INTENT(inout) :: ncid_sfc
    INTEGER, INTENT(in)  :: nt_sfc
    INTEGER, INTENT(inout) :: idum
    REAL, DIMENSION(:), INTENT(out) :: qvs,ts,&
         &ustar,hflux,qvflux, times_sfc
    
! local
    INTEGER :: i,k,kkl,kkr,imem, itran, itransfc,itransoil, itransmos
    INTEGER :: ierr, timeid
    INTEGER, DIMENSION(3) :: vec, lenvec

    INTEGER :: nmix, imix
    REAL :: rtran
    
    REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: work_sfc
    INTEGER, ALLOCATABLE, DIMENSION(:):: days_list, &
                                         secs_list
    INTEGER                           :: timediff
    REAL, ALLOCATABLE, DIMENSION(:)   :: time_sfc
    REAL, ALLOCATABLE, DIMENSION(:,:) :: offset_sfc

    LOGICAL                           :: got_all_sfces = .false.

    TYPE(time_type)                            :: requested_sfc_time

    INTEGER                                    :: &
                                                  requested_sfc_index,&
                                                  req_sfc_secs, &
                                                  req_sfc_days

    REAL                                       :: missingVal
    TYPE(time_type), ALLOCATABLE, DIMENSION(:) :: sfc_times
    TYPE(time_type)                            :: time_tolerance
    TYPE(time_type)                            :: epoch_time
    TYPE(time_type)                            :: req_sfc_time

! set up calendar and time tolerance
   call set_calendar_type(calendar_type)
   time_tolerance = set_time(600,0)                   ! 10 mins either way
   epoch_time     = set_date(EPOCH_DATE(1),EPOCH_DATE(2),EPOCH_DATE(3))

   ALLOCATE(days_list(nsfc_file), secs_list(nsfc_file),sfc_times(nsfc_file))

   ALLOCATE(time_sfc(nsfc_file))

   ierr= nf_get_vara_double(ncid_sfc,fltid(1),(/1/),(/nsfc_file/),time_sfc)
   ierr= nf_get_att_double(ncid_sfc,flid(1), "missing_value", missingVal)


! get the times

   WHERE (time_sfc < 0)
     time_sfc = missingVal
   ELSE WHERE
     secs_list = mod(time_sfc,86400.0)
     days_list = (time_sfc - secs_list) / 86400 
   END WHERE

   DO i = 1, nsfc_file
     IF (time_sfc(i) /= missingVal) sfc_times(i) = epoch_time + &
                                     set_time(secs_list(i),days_list(i))
   ENDDO


! if the year is negative we get a random one valid at the
! same time of day, otherwise the specified sfces

   IF ( rnd_init ) THEN
      nmix = 2
   ELSE
      nmix = 1
      requested_sfc_index = -9999
      requested_sfc_time = set_date(start_year_f, start_month_f, &
                                    start_day_f, start_hour_f,   &
                                    start_minute_f, 0)
      DO i = 1, nsfc_file
         IF ( sfc_times(i) >= (requested_sfc_time-time_tolerance) .AND. &
              sfc_times(i) <= (requested_sfc_time+time_tolerance) ) THEN
            requested_sfc_index = i 
         ENDIF
      ENDDO
      
      IF ( requested_sfc_index < 0 ) THEN
         print*
         print*,'Could not find requested input sfces: '
         call print_date(requested_sfc_time)
         stop 'module_sfc_init_and_bc'
      ENDIF
      IF ( requested_sfc_index < 0 ) THEN
         print*
         print*,'Could not find requested input sfces: '
         call print_date(requested_sfc_time)
         stop 'module_sfc_init_and_bc'
      ENDIF
      IF (time_sfc(requested_sfc_index+1) - time_sfc(requested_sfc_index) &
          /= interval_sfc ) THEN
         print*
         print*,'The interval in the sfc file does not agree with ',&
                'your namelist interval_f_sfc: ',interval_sfc
         print*,time_sfc(requested_sfc_index+1)-time_sfc(requested_sfc_index)
         stop 'module_sfc_init_and_bc'
      ENDIF
 
      req_sfc_time = epoch_time + set_time(req_sfc_secs,req_sfc_days)
      print*,'Using corresponding SFC date: '
      call print_date(requested_sfc_time)

   ENDIF ! want a specific sfc

! init the arrays because they will not be completely filled
   ts = -9999.
   qvs = -9999.
   ustar = -9999.
   hflux = -9999.
   qvflux = -9999.
   times_sfc = -9999.

! allocations of work arrays
   ALLOCATE(work_sfc(nsfcvars,nmix,nt_sfc))
   lenfl = (/nt_sfc/)

   DO WHILE ( .not. got_all_sfces ) 

      got_all_sfces = .true.  ! default is good
      control_index = -9999
      DO imix = 1, nmix

! choose a random record and get the data
         IF ( nmix == 1 ) THEN
            itransfc = requested_sfc_index
         ELSE
            rtran = ran1(idum)*(nsfc_file-nt_sfc) + 1.
            itran = AINT(rtran)

         ENDIF
         control_index(imix) = itran
         rtran = ran1(idum)
         control_w = rtran

         !sfces
         stfl = (/itransfc/)

         DO i = 1, nsfcvars
            ierr=nf_get_vara_double(ncid_sfc,flid(i),stfl,lenfl, &
                 work_sfc(i,imix,:))
            ierr=nf_get_att_double(ncid_sfc, flid(i), "missing_value", missingSfcVals(i))
         ENDDO

      END DO ! imix

      ! for now, be satisfied only if we have no missing values
      DO i = 1, nsfcvars
         IF ( minval(work_sfc(i,1:nmix,:)) == missingSfcVals(i) ) THEN
            got_all_sfces = .false.
         ENDIF
      ENDDO

      IF ( .not. got_all_sfces .and. nmix == 1 ) THEN 
         ! fail in this case because it is specified
         print*
         print*,'Missing data prevents a simulation of ',&
                 forecast_length,' seconds'
         stop 'module_sfc_init_and_bc'
      ELSEIF ( .not. got_all_sfces ) THEN
         print*,'Trying again to get full sfces'
      ENDIF
      
   END DO  ! still searching for sfces
   PRINT*,"Getting profile from times ",control_index
   PRINT*,"with alpha ", control_w

! do grid differencing and weighing of profiles

   SELECT CASE ( nmix )
      CASE (1)
         ts(1:nt_sfc) = work_sfc(1,1,:)
         qvs(1:nt_sfc) = work_sfc(2,1,:)
         ustar(1:nt_sfc) = work_sfc(3,1,:)
         hflux(1:nt_sfc) = work_sfc(4,1,:)
         qvflux(1:nt_sfc) = work_sfc(5,1,:)
      CASE (2)
print*, 'CANNOT CHOOSE RANDOM OBS YET!'
stop 'sfc_init_and_bc'
         CALL blend_sfc(ts,qvs,ustar,hflux,qvflux,&
              &work_sfc,control_w,nt_sfc)
      CASE DEFAULT
   END SELECT


   DEALLOCATE(work_sfc)
   DEALLOCATE(time_sfc)

! finally, forecast times  --- these are matched to the profile!
! The soil may be offset.

   DO i=1,nt_sfc
      times_sfc(i)=(i-1)*interval_sfc 
   ENDDO

!   ierr = nf_close(ncid_sfc)


 END SUBROUTINE sfc_init_and_bc

!-------------------------------------------------------------------
  SUBROUTINE blend_sfc(&
       &ts,qvs,ustar,hflux,qvflux,&
       &work_sfc,w,nt_sfc)
  
  IMPLICIT NONE
!
! Subroutine finds an average sfc of f.  A simple approach:

  INTEGER,                  INTENT(IN)         :: &
                                                  nt_sfc
  REAL,                     INTENT(IN)         :: w
  REAL, DIMENSION(nt_sfc), INTENT(OUT)        :: ts,qvs,&
       &ustar,hflux,qvflux
  REAL, DIMENSION(nsfcvars,2,nt_sfc), INTENT(IN)  :: work_sfc

  INTEGER                                      :: it, asfc, bsfc

   asfc = 1
   bsfc = 2

   DO it = 1, nt_sfc
     ts(it) = w*work_sfc(1,asfc,it) + (1.0-w)*work_sfc(1,bsfc,it)
     qvs(it) = w*work_sfc(2,asfc,it) + (1.0-w)*work_sfc(2,bsfc,it)
     ustar(it) = w*work_sfc(3,asfc,it) + (1.0-w)*work_sfc(3,bsfc,it)
     hflux(it) = w*work_sfc(4,asfc,it) + (1.0-w)*work_sfc(4,bsfc,it)
     qvflux(it) = w*work_sfc(5,asfc,it) + (1.0-w)*work_sfc(5,bsfc,it)
   ENDDO !it

   RETURN

  END SUBROUTINE blend_sfc

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
         stop 'find_matching_index:module_sfc_and_bc'
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

END MODULE module_sfc_init_and_bc
