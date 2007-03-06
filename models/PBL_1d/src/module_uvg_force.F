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
                                   rnd_force, eofs_file_forc, n_eo,&
                                   force_f_type, rnd_init, scales
  USE module_interpolations, only: linear, seval, spline

  IMPLICIT NONE

  private

  INTEGER, PARAMETER    :: npvars = 4   ! number of physical variables in profs
  INTEGER, PARAMETER    :: nptvars = 2  ! number of timing variables in profs

  INTEGER, PARAMETER    :: calendar_type = GREGORIAN
  INTEGER, DIMENSION(3) :: EPOCH_DATE = (/1970,1,1/)
  INTEGER, DIMENSION(nptvars) :: ptid  ! timing ids
  INTEGER, DIMENSION(npvars)  :: pid  ! variable ids
  INTEGER, DIMENSION(2) :: pdimid, pdimlen, lenp, stp

! for WRF uvg
  INTEGER                           :: wrf_ind, wrf_end_ind, nt_wrfin, nt
  INTEGER, PARAMETER                :: num_prof_dims = 3
  INTEGER, PARAMETER                :: num_prof_vars = 7
  INTEGER, DIMENSION(num_prof_dims) :: pdimid_wrf,pdimlen_wrf,lenp_wrf,stp_wrf
  INTEGER, DIMENSION(num_prof_dims) :: pdimid_stag_wrf,pdimlen_stag_wrf
  INTEGER, DIMENSION(num_prof_vars) :: pid_wrf

! for WRF uvg eofs
  INTEGER, PARAMETER                   :: num_eo_sfc_vars = 2
  INTEGER, PARAMETER                   :: num_prof_eo_dims = 3
  INTEGER, PARAMETER                   :: num_sfc_eo_dims = 2
  INTEGER, DIMENSION(6)                :: pid_wrf_eo
  INTEGER, DIMENSION(num_prof_eo_dims) :: pdimid_wrf_eo, pdimlen_wrf_eo
  INTEGER, DIMENSION(num_prof_eo_dims) :: pdimid_stag_wrf_eo, pdimlen_stag_wrf_eo
  INTEGER, DIMENSION(num_prof_eo_dims) :: lenp_eo, stp_eo

  INTEGER, DIMENSION(num_sfc_eo_dims)  :: sfcdimid_wrf_eo, sfcdimlen_wrf_eo
  INTEGER, DIMENSION(num_sfc_eo_dims)  :: lensfc_eo, stsfc_eo
  INTEGER, DIMENSION(num_sfc_eo_dims)  :: lenev_eo, stev_eo

  REAL,    DIMENSION(npvars)     :: missingVals ! associated netCDF missing flags 
  INTEGER               :: nrecords

  INTEGER, DIMENSION(2)               :: control_index

  REAL                                :: control_w

! I have to define the relevant things for WRF eofs and WRF uvg 
  public :: uvg_ruc_f_dims, uvg_wrf_f_dims, uvg_ruc_f_bc, uvg_wrf_f_bc

CONTAINS

  SUBROUTINE uvg_ruc_f_dims(ncid, nz, nt)

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
       PRINT*,"Problem opening forcing file uvg (0), aborting!"
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

   END SUBROUTINE uvg_ruc_f_dims
!**********************************************

  SUBROUTINE uvg_wrf_f_dims(ncid,nz,nz_stag,nrecords,nt)

!DO nz, nz_stag, nt_wrfin, nrecords come from the uvg_wrf file: available vertical
!levels, number of hours for each run, and number of runs

!DO nt comes from the namelist parameters, it is the number of times we want to run, the
!length of the column run

   IMPLICIT NONE
   INCLUDE 'netcdf.inc'

   INTEGER, INTENT(inout):: ncid
   INTEGER, INTENT(out)  :: nz, nz_stag, nrecords, nt

! local
   INTEGER              :: ierr

! In init I used interval_f, here I use interval_uvg
! some timing info
    wrf_ind = NINT(1+(REAL(start_forecast/interval_uvg)))
    nt = NINT(1+(REAL(forecast_length)) / REAL(interval_uvg))
    wrf_end_ind = wrf_ind + nt - 1

! nt is the number of times that we have the WRF uvg variables in the uvg WRF file,
! according to the forecast length that I asked and to the interval at which
! those variables were written into the file
! the column run needs WRF input starting at wrf_ind and ending at wrf_end_ind

! the times in the column variables run from 1 to nt
! the times in the WRF file run from wrf_ind to wrf_end_ind

! open forcing file
   ierr = nf_open(uvg_file, 0, ncid)
   IF ( ierr /= NF_NOERR ) THEN
      PRINT*,"Problem opening forcing file uvg (1), aborting!"
      STOP
   ENDIF

   ierr=nf_inq_dimid(ncid,'Times',pdimid_wrf(2))
   ierr=nf_inq_dimlen(ncid,pdimid_wrf(2),pdimlen_wrf(2))

! nt_wrfin is the number of hours for each WRF run in the WRF uvg file, the
! length of the WRF time series in each run
   nt_wrfin = pdimlen_wrf(2)

!   abort if not enough times in the file 
    IF ( wrf_end_ind > nt_wrfin ) THEN
       print*,"Not enough times in the wrf forcing file for this forecast"
       print*,nt_wrfin - wrf_ind + 1, &
              " forecast times available from ",start_hour_f + &
                start_forecast/interval_uvg,"Z"
       stop
    ENDIF
   
   ierr=nf_inq_dimid(ncid,'z_amsl',pdimid_wrf(1))
   ierr=nf_inq_dimlen(ncid,pdimid_wrf(1),pdimlen_wrf(1))
   ierr=nf_inq_dimid(ncid,'record',pdimid_wrf(3))
   ierr=nf_inq_dimlen(ncid,pdimid_wrf(3),pdimlen_wrf(3))

   ! fill staggered vars too
   ierr=nf_inq_dimid(ncid,'z_amsl_stag',pdimid_stag_wrf(1))
   ierr=nf_inq_dimlen(ncid,pdimid_stag_wrf(1),pdimlen_stag_wrf(1))
   ierr=nf_inq_dimid(ncid,'Times',pdimid_stag_wrf(2))
   ierr=nf_inq_dimlen(ncid,pdimid_stag_wrf(2),pdimlen_stag_wrf(2))
   ierr=nf_inq_dimid(ncid,'record',pdimid_stag_wrf(3))
   ierr=nf_inq_dimlen(ncid,pdimid_stag_wrf(3),pdimlen_stag_wrf(3))
   
   nz = pdimlen_wrf(1)
   nz_stag = pdimlen_stag_wrf(1)
   nrecords = pdimlen_wrf(3)

! variables
   ierr=nf_inq_varid(ncid,'U_G',pid_wrf(1))
   ierr=nf_inq_varid(ncid,'V_G',pid_wrf(2))
   ierr=nf_inq_varid(ncid,'Z',pid_wrf(3))
   ierr=nf_inq_varid(ncid,'inityear',pid_wrf(4))
   ierr=nf_inq_varid(ncid,'initmonth',pid_wrf(5))
   ierr=nf_inq_varid(ncid,'initday',pid_wrf(6))
   ierr=nf_inq_varid(ncid,'inithour',pid_wrf(7))

! close 
   ierr = nf_close(ncid)

   END SUBROUTINE uvg_wrf_f_dims
!**********************************************

  SUBROUTINE uvg_eof_f_dims(ncid_eofs_forc)

! As of now this subroutine is only for WRF eofs, not eofs subroutine for RUC  

   IMPLICIT NONE
   INCLUDE 'netcdf.inc'

   INTEGER, INTENT(inout):: ncid_eofs_forc

! local
   INTEGER              :: ierr

! open eofs forcing file
   ierr = nf_open(eofs_file_forc, 0, ncid_eofs_forc)
   
   IF ( ierr /= NF_NOERR ) THEN
      PRINT*,"Problem opening eofs uvg forcing file, ",eofs_file_forc
      PRINT*,"aborting!"
      STOP
   ENDIF

  ! profile dimensions
   ierr=nf_inq_dimid(ncid_eofs_forc,'N',pdimid_wrf_eo(1))
   ierr=nf_inq_dimlen(ncid_eofs_forc,pdimid_wrf_eo(1),pdimlen_wrf_eo(1))
   ierr=nf_inq_dimid(ncid_eofs_forc,'z_amsl',pdimid_wrf_eo(2))
   ierr=nf_inq_dimlen(ncid_eofs_forc,pdimid_wrf_eo(2),pdimlen_wrf_eo(2))
   ierr=nf_inq_dimid(ncid_eofs_forc,'Times',pdimid_wrf_eo(3))
   ierr=nf_inq_dimlen(ncid_eofs_forc,pdimid_wrf_eo(3),pdimlen_wrf_eo(3))

  ! staggered profile dimensions
   ierr=nf_inq_dimid(ncid_eofs_forc,'N',pdimid_stag_wrf_eo(1))
   ierr=nf_inq_dimlen(ncid_eofs_forc,pdimid_stag_wrf_eo(1),pdimlen_stag_wrf_eo(1))
   ierr=nf_inq_dimid(ncid_eofs_forc,'z_amsl_stag',pdimid_stag_wrf_eo(2))
   ierr=nf_inq_dimlen(ncid_eofs_forc,pdimid_stag_wrf_eo(2),pdimlen_stag_wrf_eo(2))
   ierr=nf_inq_dimid(ncid_eofs_forc,'Times',pdimid_stag_wrf_eo(3))
   ierr=nf_inq_dimlen(ncid_eofs_forc,pdimid_stag_wrf_eo(3),pdimlen_stag_wrf_eo(3))

  ! surface dimensions
   ierr=nf_inq_dimid(ncid_eofs_forc,'N',sfcdimid_wrf_eo(1))
   ierr=nf_inq_dimlen(ncid_eofs_forc,sfcdimid_wrf_eo(1),sfcdimlen_wrf_eo(1))
   ierr=nf_inq_dimid(ncid_eofs_forc,'Times',sfcdimid_wrf_eo(2))
   ierr=nf_inq_dimlen(ncid_eofs_forc,sfcdimid_wrf_eo(2),sfcdimlen_wrf_eo(2))

! variables
   ierr=nf_inq_varid(ncid_eofs_forc,'U_G',pid_wrf_eo(1))
   ierr=nf_inq_varid(ncid_eofs_forc,'V_G',pid_wrf_eo(2))
   ierr=nf_inq_varid(ncid_eofs_forc,'Z',pid_wrf_eo(3))
   ierr=nf_inq_varid(ncid_eofs_forc,'GLW_G',pid_wrf_eo(4))
   ierr=nf_inq_varid(ncid_eofs_forc,'GSW_G',pid_wrf_eo(5))
   ierr=nf_inq_varid(ncid_eofs_forc,'EIGENVALS',pid_wrf_eo(6))

! close 
   ierr = nf_close(ncid_eofs_forc)

   END SUBROUTINE uvg_eof_f_dims

!**********************************************
  SUBROUTINE uvg_ruc_f_bc(ncid, nz, nt, terrain, nz_grid, z_grid, &
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
       PRINT*,"Problem opening forcing file uvg (2), aborting!"
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

   IF ( rnd_force==1 ) THEN
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

  END SUBROUTINE uvg_ruc_f_bc

!**********************************************************************************
! DO I am adding here more arguments
 SUBROUTINE uvg_wrf_f_bc(ncid, nz, nt, terrain, nz_grid, z_grid, &
                   u_g_f,v_g_f, nsplinetimes, splinetimes, &
                   glw_init, gsw_init, &
                   glw_f, gsw_f, &
                   ntimes_f_flux,nsplinetimes_flux, splinetimes_flux, &
                   init_f_type, &
                   idum, ims,ime,jms,jme,kms,kme,ncid_eofs_forc, &
                   itran1,itran2,control_w,gasdo,nz_stag,nrecords)

    USE module_nr_procedures

    IMPLICIT NONE
    INCLUDE 'netcdf.inc'

! Subroutine gets uvg forcing from WRF , and puts
! them on the column grid

! arguments

!DO I add here to pass from init routine 
     INTEGER, INTENT(in) :: itran1, itran2
     REAL, INTENT(in) :: control_w
     REAL, DIMENSION(:), INTENT(in):: gasdo

    INTEGER, INTENT(inout) :: ncid, ncid_eofs_forc
    INTEGER, INTENT(in)  :: nz, nt, nz_grid 
    INTEGER, INTENT(in)  :: ntimes_f_flux, nsplinetimes, nsplinetimes_flux
! DO I added those
    INTEGER, INTENT(in)  :: nz_stag, nrecords
    
    REAL,    INTENT(in)  :: terrain
    REAL, DIMENSION(:,:,:), INTENT(in):: z_grid
    INTEGER, INTENT(inout) :: idum
    REAL, DIMENSION(nsplinetimes_flux), INTENT(in) :: splinetimes_flux
    REAL, DIMENSION(nsplinetimes), INTENT(in) :: splinetimes
    REAL, DIMENSION(:,:), INTENT(out) :: u_g_f,v_g_f
    REAL, DIMENSION(ntimes_f_flux), INTENT(in) :: gsw_init,glw_init
    REAL, DIMENSION(:), INTENT(out) :: gsw_f,glw_f
    CHARACTER(len=*), INTENT(in) :: init_f_type
    INTEGER, INTENT(in) :: ims,ime,jms,jme,kms,kme

! DO I add here some WRF stuff
    INTEGER :: nt_wrf_uvg
! local
    INTEGER :: i,k,kkl,kkr,imem, itran
    INTEGER :: ierr, timeid, iwrf_err
    INTEGER, DIMENSION(3) :: vec, lenvec

    INTEGER :: nmix, imix
    REAL :: rtran
    
! to do spline interp of fluxes
    REAL, DIMENSION(ntimes_f_flux) :: aatmp_flux,bbtmp_flux,cctmp_flux

! Arrays to work with the WRF variables   
    INTEGER, DIMENSION(:), ALLOCATABLE:: wrf_year, wrf_month, wrf_day, wrf_hour
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ug, vg, zg
    REAL, ALLOCATABLE, DIMENSION(:    ) :: glw,gsw
    REAL, ALLOCATABLE, DIMENSION(:)   :: times
    REAL, ALLOCATABLE, DIMENSION(:,:) :: z, z_agl, p, u, v, utmp, vtmp
    REAL, ALLOCATABLE, DIMENSION(:,:) :: ug1

!DO. Add arrays for eofs of profiles
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: un, vn, zn 
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: surfn

! now declare the eigenvalues
    REAL, ALLOCATABLE, DIMENSION(:,:) :: eval
!for the eofs based perturbations
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: pun, pvn, pzn
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: psur
    REAL, ALLOCATABLE, DIMENSION(:,:) :: put, pvt, pzt
    REAL, ALLOCATABLE, DIMENSION(:,:) :: psurt
    REAL, ALLOCATABLE, DIMENSION(:,:) :: auxfac
    REAL, ALLOCATABLE, DIMENSION(:,:) :: damped_scale
    INTEGER :: ivar,iz,ieo, fac, v1, v2, s1, s2, rsq
    REAL :: facs

    INTEGER                                    :: requested_f_index
    REAL                                       :: missingVal
    TYPE(time_type)      :: requested_wrf_time

    INTEGER::i1,j1,it


   IF(rnd_force==3)THEN

     call uvg_eof_f_dims(ncid_eofs_forc)

     ! open file again
    ierr = nf_open(eofs_file_forc, 0, ncid_eofs_forc)

! profiles
     ALLOCATE(un(n_eo,pdimlen_wrf_eo(2),pdimlen_wrf_eo(3)))
     ALLOCATE(vn(n_eo,pdimlen_wrf_eo(2),pdimlen_wrf_eo(3)))
     ALLOCATE(zn(n_eo,pdimlen_stag_wrf_eo(2),pdimlen_wrf_eo(3)))
! surface
     ALLOCATE(surfn(num_eo_sfc_vars,n_eo,sfcdimlen_wrf_eo(2)))
! the eigenvalues
     ALLOCATE(eval(n_eo,pdimlen_wrf_eo(3)))
! auxiliar working array
     ALLOCATE(auxfac(n_eo,pdimlen_wrf_eo(3)))
! allocate the perturbations for each mode
     ALLOCATE(pun(n_eo,nz,nt_wrfin),pvn(n_eo,nz,nt_wrfin),pzn(n_eo,nz_stag,nt_wrfin))
     ALLOCATE(psur(num_eo_sfc_vars,n_eo,sfcdimlen_wrf_eo(2)))

    ENDIF ! eof_init (rnd_force == 3)

  ALLOCATE(utmp(nz_grid,nt),vtmp(nz_grid,nt))

! set up calendar and time tolerance
! for now this is only for information
     call set_calendar_type(calendar_type)
     
     requested_wrf_time = set_date(start_year_f, start_month_f, &
                                    start_day_f, start_hour_f,   &
                                    start_minute_f, 0)


! open forcing file
   ierr = nf_open(uvg_file, 0, ncid)
   IF ( ierr /= NF_NOERR ) THEN
       PRINT*,"Problem opening forcing file uvg (3), aborting!"
      STOP
   ENDIF

! DO I will copy into here some stuff from wrf_init_and_bc to deal with different
! forcings using WRF data

! use requested sounding or random?
! DO. rnd_init changed to integer. 1=undated random mixture, 2=dated perfect mode, 3=dated eofs.
    IF ( rnd_force /= 1) THEN
! DO this option means that we chose dated initialization    
       requested_f_index = -9999
       ALLOCATE(wrf_year(nrecords), wrf_month(nrecords), &
                wrf_day(nrecords), wrf_hour(nrecords))

       ierr=nf_get_vara_int(ncid,pid_wrf(4),(/1,1/),(/1,nrecords/),wrf_year)
       ierr=nf_get_vara_int(ncid,pid_wrf(5),(/1,1/),(/1,nrecords/),wrf_month)
       ierr=nf_get_vara_int(ncid,pid_wrf(6),(/1,1/),(/1,nrecords/),wrf_day)
       ierr=nf_get_vara_int(ncid,pid_wrf(7),(/1,1/),(/1,nrecords/),wrf_hour)

       DO i = 1, nrecords
          IF ( wrf_year(i)  == start_year_f .and. &
               wrf_month(i) == start_month_f .and. &
               wrf_day(i)   == start_day_f .and. &
               wrf_hour(i)  == start_hour_f ) &
               requested_f_index = i
       ENDDO

       IF ( requested_f_index < 0 ) THEN
          print*,"Could not find requested WRF date: "
          print*,start_year_f, start_month_f, start_day_f, start_hour_f
          print*,'in file ',uvg_file
          stop 'uvg_wrf_f_bc'
       ENDIF
          
       print*,'Getting requested WRF uvg forcing valid at date: '
       call print_date(requested_wrf_time)
       print*
       print*,'If you are running DART, the correct time ',&
              ' for the obs_sequence is near'
       call print_time(requested_wrf_time)

    ENDIF

   IF ( rnd_force /= 2) THEN
    nmix=2
   ELSE
    nmix=1
   ENDIF

! The dimensions come from uvg_wrf_f_dims
! I could avoid the variable nt_wrfin and just used the pdimlen_wrf I got from
! the files, but it it seems shorter and more convenient to use nt_wrfin

   ALLOCATE(ug(nz,nt_wrfin,nmix),vg(nz,nt_wrfin,nmix),zg(nz_stag,nt_wrfin,nmix))
   ALLOCATE(glw(nt),gsw(nt))
  
! The nt dimension here is the length of the time series I want to run with the
! column, it was calculated using the namelist
   ALLOCATE(z(nz_stag,nt), z_agl(nz,nt), u(nz,nt), v(nz,nt), p(nz,nt))

! profiles perturbations, dimensions have to be compatible with wrf profiles
! dimensions

! I calculate the perturbations for all the times I have in the files, then I
! will take the ones I need for the series I run in the column. For now I leave
! like this


 IF(rnd_force==3)THEN
   ALLOCATE(put(nz,nt_wrfin),pvt(nz,nt_wrfin),pzt(nz_stag,nt_wrfin))
   ALLOCATE(damped_scale(nz_stag,nt_wrfin))
   ALLOCATE(psurt(num_eo_sfc_vars,nt_wrfin))
   put=0
   pvt=0
   pzt=0
   psurt=0
 ENDIF  

! I have to read the U,V,Z from the uvg WRF file
  ierr=nf_get_vara_double(ncid,pid_wrf(1),(/1,1,itran1/),(/nz,nt_wrfin,1/),ug(:,:,1))
  ierr=nf_get_vara_double(ncid,pid_wrf(2),(/1,1,itran1/),(/nz,nt_wrfin,1/),vg(:,:,1))
  ierr=nf_get_vara_double(ncid,pid_wrf(3),(/1,1,itran1/),(/nz_stag,nt_wrfin,1/),zg(:,:,1))

! the base state glw and gsw come on input
  glw = glw_init
  gsw = gsw_init

  IF ( rnd_force == 1 ) THEN
    ierr=nf_get_vara_double(ncid,pid_wrf(1),(/1,1,itran2/),(/nz,nt_wrfin,1/),ug(:,:,2))
    ierr=nf_get_vara_double(ncid,pid_wrf(2),(/1,1,itran2/),(/nz,nt_wrfin,1/),vg(:,:,2))
    ierr=nf_get_vara_double(ncid,pid_wrf(3),(/1,1,itran2/),(/nz_stag,nt_wrfin,1/),zg(:,:,2))

    ! don't do anything for u_g and v_g because it is already consistent
  ENDIF

 IF(rnd_force==3)THEN ! add eofs
! Check dimensions compatibility between WRF forcing variables and eofs

! Times
   IF (pdimlen_wrf(2) /= pdimlen_wrf_eo(3)) THEN
      PRINT*,"Times in WRF forc file and in eofs forc file differ,aborting"
      STOP
   ENDIF

! Vertical levels, z_amsl
  IF (pdimlen_wrf(1) /= pdimlen_wrf_eo(2)) THEN
     PRINT*,"# of vertical levels z_amsl in WRF forc and eofs forc differ, aborting"
     STOP
  ENDIF

! Vertical staggered levels, z_amsl_stag
  IF (pdimlen_stag_wrf(2) /= pdimlen_stag_wrf_eo(3)) THEN
     PRINT*,"# of vertical levels z_amsl_stag in WRF forc and eofs forc differ, aborting"
     STOP
  ENDIF

! eofs of profile
      stp_eo = (/1,1,1/)
      lenp_eo = (/n_eo,pdimlen_wrf_eo(2),pdimlen_wrf_eo(3)/)
      ierr=nf_get_vara_double(ncid_eofs_forc,pid_wrf_eo(1),stp_eo,lenp_eo,un)
      ierr=nf_get_vara_double(ncid_eofs_forc,pid_wrf_eo(2),stp_eo,lenp_eo,vn)

! z eofs are in staggered vertical levels      
      lenp_eo = (/n_eo,pdimlen_stag_wrf_eo(2),pdimlen_stag_wrf_eo(3)/)
      ierr=nf_get_vara_double(ncid_eofs_forc,pid_wrf_eo(3),stp_eo,lenp_eo,zn)
      
! eofs of gls and gsw
      stsfc_eo = (/1,1/)
      lensfc_eo = (/n_eo,sfcdimlen_wrf_eo(2)/)
      ierr=nf_get_vara_double(ncid_eofs_forc,pid_wrf_eo(4),stsfc_eo,lensfc_eo,surfn(1,:,:))
      ierr=nf_get_vara_double(ncid_eofs_forc,pid_wrf_eo(5),stsfc_eo,lensfc_eo,surfn(2,:,:))

! read eigenvalues
      stev_eo = (/1,1/)
      lenev_eo = (/n_eo,pdimlen_wrf_eo(3)/)
      ierr=nf_get_vara_double(ncid_eofs_forc,pid_wrf_eo(6),stev_eo,lenev_eo,eval)

! Do eofs perturbations stuff

! auxfac is the weight of the perturbation
! I brought gasdo from init 
    do ieo=1,n_eo
      auxfac(ieo,:)=sqrt(eval(ieo,:))*gasdo(ieo)
    enddo

! calculate first the term for each eof
! profiles
    do iz=1,nz
      pun(:,iz,:)=un(:,iz,:)*auxfac(:,:)
      pvn(:,iz,:)=vn(:,iz,:)*auxfac(:,:)
    enddo
     
! z is in staggered vertical levels
! I am avoiding the first level because it has only missing values. Easiest
! though not smartest way to do it    

    do iz=2,nz_stag
      pzn(:,iz,:)=zn(:,iz,:)*auxfac(:,:)
    enddo

! surface
    do ivar=1,num_eo_sfc_vars
      psur(ivar,:,:)=surfn(ivar,:,:)*auxfac(:,:)
    enddo

    ! damp the scales so that forcing spread is 0 at top
    do k = 1,nz_stag
      do i = 1,nt_wrfin
        damped_scale(k,i) = scales * cos(dble(k)/dble(nz_stag)*acos(-1.0)/2)**2
      enddo
    enddo
    damped_scale = 1.0 ! replace with constant

! calculate the perturbation over all eofs
    do ieo=1,n_eo
! profiles perturbations of time and vertical level
     put(:,:)=put(:,:)+pun(ieo,:,:)*damped_scale(1:nz,:)
     pvt(:,:)=pvt(:,:)+pvn(ieo,:,:)*damped_scale(1:nz,:)
     pzt(2:nz_stag,:)=pzt(2:nz_stag,:)+pzn(ieo,2:nz_stag,:)*damped_scale
     psurt(:,:)=psurt(:,:)+ psur(:,ieo,:)*scales 
    enddo

  ENDIF

! Here I use from 1:nt for the z,u,v and wrf_ind, wrf_end_ind for the
! zg,ug,vg,and pzt, put, pvt
! DO here I have z in staggered levels, z_agl will be in non-staggered levels

   IF ( rnd_force==1 )THEN
!         z(:,1:nt) = zg(:,wrf_ind:wrf_end_ind,1)*control_w + &
!                     zg(:,wrf_ind:wrf_end_ind,2)*(1.0-control_w)
         z(:,1:nt) = MAX(0.,&
                     zg(:,wrf_ind:wrf_end_ind,1)*control_w + &
                     zg(:,wrf_ind:wrf_end_ind,2)*(1.0-control_w)&
                     -terrain)
     
      DO k=1,nz
       z_agl(k,1:nt)=.5*(&
             (zg(k,wrf_ind:wrf_end_ind,1)+zg(k+1,wrf_ind:wrf_end_ind,1))&
             *control_w & 
             +(zg(k,wrf_ind:wrf_end_ind,2)+zg(k+1,wrf_ind:wrf_end_ind,2))&
             *(1.0-control_w))-terrain
      ENDDO       
                     
         u(:,1:nt) = ug(:,wrf_ind:wrf_end_ind,1)*control_w + &
                     ug(:,wrf_ind:wrf_end_ind,2)*(1.0-control_w)
         v(:,1:nt) = vg(:,wrf_ind:wrf_end_ind,1)*control_w + &
                     vg(:,wrf_ind:wrf_end_ind,2)*(1.0-control_w)
   ELSEIF ( rnd_force==2 )THEN
         z(:,1:nt) = MAX(0.,zg(:,wrf_ind:wrf_end_ind,1)-terrain)

      DO k=1,nz
       z_agl(k,1:nt)=.5*(&
             zg(k,wrf_ind:wrf_end_ind,1)+zg(k+1,wrf_ind:wrf_end_ind,1))&
             -terrain
      ENDDO      
         
         u(:,1:nt) = ug(:,wrf_ind:wrf_end_ind,1)
         v(:,1:nt) = vg(:,wrf_ind:wrf_end_ind,1)

   ELSEIF ( rnd_force==3 )THEN

         z(:,1:nt) = MAX(0.,zg(:,wrf_ind:wrf_end_ind,1)+pzt(:,wrf_ind:wrf_end_ind)&
                     -terrain)
                     
       DO k=1,nz
       z_agl(k,1:nt)=.5*(&
              zg(k,wrf_ind:wrf_end_ind,1)+zg(k+1,wrf_ind:wrf_end_ind,1)&
             +pzt(k,wrf_ind:wrf_end_ind)+pzt(k+1,wrf_ind:wrf_end_ind))&
             -terrain
             
       ENDDO

       u(:,1:nt) = ug(:,wrf_ind:wrf_end_ind,1)+put(:,wrf_ind:wrf_end_ind)
       v(:,1:nt) = vg(:,wrf_ind:wrf_end_ind,1)+pvt(:,wrf_ind:wrf_end_ind)         
       glw(1:nt) = glw(1:nt)+psurt(1,wrf_ind:wrf_end_ind)
       gsw(1:nt) = gsw(1:nt)+psurt(2,wrf_ind:wrf_end_ind)

       WHERE(gsw < 0.0) gsw = 0.0
   ENDIF    


   ierr = nf_close(ncid)
   ierr = nf_close(ncid_eofs_forc)

! finally, forcing times  --- these are matched to the profile!

   ALLOCATE(times(nt))

   DO i=1,nt
      times(i)=(i-1)*interval_uvg+start_forecast
   ENDDO

! just go ahead and interpolate to model grid and spline times here

   DO i = 1,nt
      DO k = 1, nz_grid

!         utmp(k,i) = linear(z_agl(:,i),u(:,i),nz,z_grid(1,k,1))
!         vtmp(k,i) = linear(z_agl(:,i),v(:,i),nz,z_grid(1,k,1))

         utmp(k,i) = linear(z(:,i),u(:,i),nz,z_grid(1,k,1))
         vtmp(k,i) = linear(z(:,i),v(:,i),nz,z_grid(1,k,1))

      ENDDO
   ENDDO

   DO k = 1, nz_grid
      DO i = 1, nsplinetimes
        u_g_f(k,i) = linear(times,utmp(k,:),nt,splinetimes(i))
        v_g_f(k,i) = linear(times,vtmp(k,:),nt,splinetimes(i))
      ENDDO
   ENDDO

   CALL spline(ntimes_f_flux,times,glw,aatmp_flux,bbtmp_flux,cctmp_flux)
   DO i = 1, nsplinetimes_flux
        glw_f(i) = seval(ntimes_f_flux,splinetimes_flux(i),times,glw,&
               &aatmp_flux,bbtmp_flux,cctmp_flux)
        gsw_f(i) = seval(ntimes_f_flux,splinetimes_flux(i),times,gsw,&
               &aatmp_flux,bbtmp_flux,cctmp_flux)
   ENDDO

   DEALLOCATE(u,v,p,z,utmp,vtmp)

 IF(rnd_force==3)&
   DEALLOCATE(un,vn,zn,pun,pvn,pzn,auxfac,put,pvt,pzt,eval,glw,gsw)

   RETURN

  END SUBROUTINE uvg_wrf_f_bc

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

