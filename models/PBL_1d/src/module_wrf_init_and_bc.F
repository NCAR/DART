MODULE module_wrf_init_and_bc

  USE time_manager_mod,      only: time_type, GREGORIAN, &
                                   set_calendar_type, print_time, &
                                   print_date, set_date, set_time
  USE module_model_constants
  USE module_namelist,       only: start_hour_f, interval_f, &
                                   start_forecast, forecast_length, &
                                   start_year_f, start_month_f, & 
                                   start_day_f, start_minute_f, &
                                   init_f_file, isltyp_ref, &
                                   ivgtyp_ref, lu_index_ref, vegfra_ref, &
                                   rnd_init, rnd_force, t_advection, &
                                   qv_advection

  IMPLICIT NONE
  private
  INTEGER                ::  midx, midy, wrf_ind, wrf_end_ind, nt_wrfin
  INTEGER                :: nrecords
  INTEGER, DIMENSION(24) :: fid  !profiles, screen, surface, soil, z_agl
  INTEGER, DIMENSION(5) :: sdimid, sdimlen, lens, sts
  INTEGER, DIMENSION(5) :: sfcdimid, sfcdimlen, lensfc, stsfc
  INTEGER, DIMENSION(6) :: sldimid, sldimlen, lensl, stsl
  INTEGER, DIMENSION(5) :: pdimid, pdimlen, lenp, stp
  INTEGER, DIMENSION(5) :: pdimid_stag, pdimlen_stag, &
                           lenp_x_stag, lenp_y_stag, lenp_z_stag
  INTEGER, DIMENSION(4) :: stprecip, lenprecip
  INTEGER, DIMENSION(4) :: stter, lenter

  INTEGER, DIMENSION(2)               :: control_index

  INTEGER, PARAMETER    :: calendar_type = GREGORIAN

  REAL                                :: control_w
 
  public           :: wrf_f_dims, wrf_init_and_bc

CONTAINS

  SUBROUTINE wrf_f_dims(ncid, nz, ns, nt, dx,cent_lat,cent_lon, &
                        truelat1,truelat2,mminlu,julday,&
                        sw_corner_lon,sw_corner_lat,projcode,&
                        lat, lon, cor)
     IMPLICIT NONE
     INCLUDE 'netcdf.inc'

! gets all static data from the input file, including dimensions

     INTEGER, INTENT(inout):: ncid
     INTEGER, INTENT(out)  :: nz, ns, nt
     INTEGER, INTENT(out)  :: julday,projcode
     REAL, INTENT(out)     :: lat,lon,cor,dx
     REAL, INTENT(out)     :: cent_lat,cent_lon,truelat1,truelat2,&
                              sw_corner_lon,sw_corner_lat
     CHARACTER(len=4), INTENT(out) :: mminlu

! local
     INTEGER              :: ierr
     INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: isltyp2d,ivgtyp2d,lu_index2d
     REAL, ALLOCATABLE, DIMENSION(:,:,:) :: lat2d,lon2d,vegfra2d

! some timing info
!    wrf_ind = NINT(1+(REAL(start_forecast/3600.0 - start_hour_f))/&
!         &REAL(interval_f/3600.))
    wrf_ind = NINT(1+(REAL(start_forecast/interval_f)))
    nt = NINT(1+(REAL(forecast_length)) / REAL(interval_f))
    wrf_end_ind = wrf_ind + nt - 1

! open forcing file
     ierr = nf_open(init_f_file, 0, ncid)
     IF ( ierr /= NF_NOERR ) THEN
        PRINT*,"Problem opening forcing file, aborting!"
        STOP
     ENDIF

! get global attributes
    ierr= nf_get_att_int(ncid, NF_GLOBAL,'MAP_PROJ',projcode)
    ierr= nf_get_att_double(ncid, NF_GLOBAL,'DX',dx)
    ierr= nf_get_att_double(ncid, NF_GLOBAL,'CEN_LAT',cent_lat)
    ierr= nf_get_att_double(ncid, NF_GLOBAL,'CEN_LON',cent_lon)
    ierr= nf_get_att_double(ncid, NF_GLOBAL,'TRUELAT1',truelat1)
    ierr= nf_get_att_double(ncid, NF_GLOBAL,'TRUELAT2',truelat2)
    ierr= nf_get_att_text(ncid, NF_GLOBAL,'MMINLU',mminlu)
    ierr= nf_get_att_int(ncid, NF_GLOBAL,'JULDAY',julday)
    ierr= nf_get_att_double(ncid, NF_GLOBAL,'SW_LON',sw_corner_lon)
    ierr= nf_get_att_double(ncid, NF_GLOBAL,'SW_LAT',sw_corner_lat)

    ierr=nf_inq_dimid(ncid,'x',pdimid(1))
    ierr=nf_inq_dimlen(ncid,pdimid(1),pdimlen(1))
!    PRINT *,'pdimlen(1) =',pdimlen(1)
    ierr=nf_inq_dimid(ncid,'y',pdimid(2))
    ierr=nf_inq_dimlen(ncid,pdimid(2),pdimlen(2))
!    PRINT *,'pdimlen(2) =',pdimlen(2)
    ierr=nf_inq_dimid(ncid,'z_amsl',pdimid(3))
    ierr=nf_inq_dimlen(ncid,pdimid(3),pdimlen(3))
!    PRINT *,'pdimlen(3) =',pdimlen(3)
    ierr=nf_inq_dimid(ncid,'x_stag',pdimid_stag(1))
    ierr=nf_inq_dimlen(ncid,pdimid_stag(1),pdimlen_stag(1))
!    PRINT *,'pdimlen_stag(1) =',pdimlen_stag(1)
    ierr=nf_inq_dimid(ncid,'y_stag',pdimid_stag(2))
    ierr=nf_inq_dimlen(ncid,pdimid_stag(2),pdimlen_stag(2))
!    PRINT *,'pdimlen_stag(2) =',pdimlen_stag(2)
    ierr=nf_inq_dimid(ncid,'z_amsl_stag',pdimid_stag(3))
    ierr=nf_inq_dimlen(ncid,pdimid_stag(3),pdimlen_stag(3))
!    PRINT *,'pdimlen_stag(3) =',pdimlen_stag(3)

    ierr=nf_inq_dimid(ncid,'Times',pdimid(4))
    ierr=nf_inq_dimlen(ncid,pdimid(4),pdimlen(4))
    nt_wrfin = pdimlen(4)

    ierr=nf_inq_dimid(ncid,'record',pdimid(5))
    ierr=nf_inq_dimlen(ncid,pdimid(5),pdimlen(5))
    nrecords = pdimlen(5)

!   abort if not enough times in the file 
    IF ( wrf_end_ind > nt_wrfin ) THEN
       print*,"Not enough times in the wrf input file for this forecast"
       print*,nt_wrfin - wrf_ind + 1, &
              " forecast times available from ",start_hour_f + &
                start_forecast/interval_f,"Z"
       stop
    ENDIF

! screen info
    ierr=nf_inq_dimid(ncid,'variable_scr',sdimid(1))
    ierr=nf_inq_dimlen(ncid,sdimid(1),sdimlen(1))
    ierr=nf_inq_dimid(ncid,'Times',sdimid(2))
    ierr=nf_inq_dimlen(ncid,sdimid(2),sdimlen(2))
    ierr=nf_inq_dimid(ncid,'record',sdimid(3))
    ierr=nf_inq_dimlen(ncid,sdimid(3),sdimlen(3))
    lens = (/sdimlen(1),pdimlen(1),pdimlen(2),sdimlen(2),1/)

! surface info
    ierr=nf_inq_dimid(ncid,'variable_sfc',sfcdimid(1))
    ierr=nf_inq_dimlen(ncid,sfcdimid(1),sfcdimlen(1))
    ierr=nf_inq_dimid(ncid,'Times',sfcdimid(2))
    ierr=nf_inq_dimlen(ncid,sfcdimid(2),sfcdimlen(2))
    ierr=nf_inq_dimid(ncid,'record',sfcdimid(3))
    ierr=nf_inq_dimlen(ncid,sfcdimid(3),sfcdimlen(3))
    lensfc = (/sfcdimlen(1),pdimlen(1),pdimlen(2),sfcdimlen(2),1/)

! soil info
    ierr=nf_inq_dimid(ncid,'variable_soil',sldimid(1))
    ierr=nf_inq_dimlen(ncid,sldimid(1),sldimlen(1))
    ierr=nf_inq_dimid(ncid,'soil_levels',sldimid(2))
    ierr=nf_inq_dimlen(ncid,sldimid(2),sldimlen(2))
    ierr=nf_inq_dimid(ncid,'Times',sldimid(3))
    ierr=nf_inq_dimlen(ncid,sldimid(3),sldimlen(3))
    ierr=nf_inq_dimid(ncid,'record',sldimid(4))
    ierr=nf_inq_dimlen(ncid,sldimid(4),sldimlen(4))
    lensl = (/sldimlen(1),pdimlen(1),pdimlen(2),sldimlen(2),sldimlen(3),&
         &1/)

! variables
    ierr=nf_inq_varid(ncid,'U',fid(1))
    ierr=nf_inq_varid(ncid,'V',fid(2))
    ierr=nf_inq_varid(ncid,'T',fid(3))
    ierr=nf_inq_varid(ncid,'Q',fid(4))
    ierr=nf_inq_varid(ncid,'P',fid(5))
    ierr=nf_inq_varid(ncid,'Z',fid(6))
    ierr=nf_inq_varid(ncid,'screen',fid(7))
    ierr=nf_inq_varid(ncid,'surface',fid(8))
    ierr=nf_inq_varid(ncid,'soil',fid(9))
    ierr=nf_inq_varid(ncid,'lats',fid(10))
    ierr=nf_inq_varid(ncid,'lons',fid(11))
    ierr=nf_inq_varid(ncid,'terrain',fid(12))
    ierr=nf_inq_varid(ncid,'inityear',fid(13))
    ierr=nf_inq_varid(ncid,'initmonth',fid(14))
    ierr=nf_inq_varid(ncid,'initday',fid(15))
    ierr=nf_inq_varid(ncid,'inithour',fid(16))

    ierr=nf_inq_varid(ncid,'MU',fid(17))
    ierr=nf_inq_varid(ncid,'MUB',fid(18))
    ierr=nf_inq_varid(ncid,'MU0',fid(19))
    ierr=nf_inq_varid(ncid,'ZNU',fid(20))
    ierr=nf_inq_varid(ncid,'ZNW',fid(21))
    ierr=nf_inq_varid(ncid,'P_TOP',fid(22))
    ierr=nf_inq_varid(ncid,'MAPFAC_M',fid(23))
    ierr=nf_inq_varid(ncid,'PRECIP',fid(24))

    ns = sldimlen(2)
    nz = pdimlen(3)

! lat, lon, terrain  

    lenter = (/pdimlen(1),pdimlen(2),sdimlen(2),1/)
    stter = (/1,1,1,1/)

    ALLOCATE(lat2d(pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(lon2d(pdimlen(1),pdimlen(2),sdimlen(2)))

    ierr=nf_get_vara_double(ncid,fid(10),stter,lenter,lat2d)
    ierr=nf_get_vara_double(ncid,fid(11),stter,lenter,lon2d)

    midx=(1+pdimlen(1))/2
    midy=(1+pdimlen(2))/2

    lat=lat2d(midx,midy,1)*DEGRAD
    lon=lon2d(midx,midy,1)*DEGRAD
    cor=2.*eomeg*SIN(lat)

    DEALLOCATE(lat2d)
    DEALLOCATE(lon2d)

    ierr = nf_close(ncid)

  END SUBROUTINE wrf_f_dims

!**********************************************

  SUBROUTINE wrf_init_and_bc(ncid,nz, ns, nt, &
       z_agl,z_agl_stag,t,u,v,q,p, &
       t_ups_x, t_ups_y, &
       q_ups_x, q_ups_y, &
       tau_u, tau_v,&
       th2,t2,tsk,u10,v10,     &
       q2,glw,gsw,qsfc,tslb,smois,tmn,&
       precip, &
       vegfra,isltyp,lu_index,ivgtyp, terrain, dx, &
       times,times_flux,times_soil,times_smos,idum)

    USE module_nr_procedures

    IMPLICIT NONE
    INCLUDE 'netcdf.inc'

! Subroutine selects a random forecast from file wrf_cntl.nc

! arguments

    INTEGER, INTENT(inout) :: ncid
    INTEGER, INTENT(in)  :: nz, ns, nt
    INTEGER, INTENT(inout) :: idum
    REAL, DIMENSION(:), INTENT(out) :: &
         th2,t2,tsk,u10,v10,q2,glw,gsw,qsfc, precip, &
         tmn,vegfra,times,times_flux, times_soil, times_smos
    INTEGER, DIMENSION(:), INTENT(out) :: isltyp,ivgtyp,lu_index
    REAL, INTENT(out)                  :: terrain
    REAL, INTENT( in)                  :: dx

    REAL, DIMENSION(:,:), INTENT(out) :: t,u,v,q,p,tslb,smois
    REAL, DIMENSION(:,:), INTENT(out)   :: tau_u, tau_v
    REAL, DIMENSION(:,:), INTENT(out) :: t_ups_x, t_ups_y
    REAL, DIMENSION(:,:), INTENT(out) :: q_ups_x, q_ups_y
    
! local
    INTEGER :: i,k,kkl,kkr,imem, itran1, itran2
    INTEGER :: ierr, nmix, imix
    INTEGER, DIMENSION(3) :: vec, lenvec

    INTEGER                           :: requested_f_index
    INTEGER, DIMENSION(:), ALLOCATABLE:: wrf_year, wrf_month, wrf_day, wrf_hour
    REAL :: rtran, rtran2
    REAL :: t_l,t_m,t_r,dt_dz,dp
    REAL, DIMENSION(:,:) :: z_agl,z_agl_stag
    
   
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: u3d1,v3d1,t3d1,q3d1,p3d1, &
         &u3d2, v3d2, t3d2, q3d2, p3d2,z3d1,z3d2
    REAL, ALLOCATABLE, DIMENSION(:,:) :: invrho,p_r,p_l,p_u,p_d
    REAL, ALLOCATABLE, DIMENSION(:,:) :: z_l,z_m,z_r,z_d,z_u


    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: screen2d1, screen2d2,&
         &surf2d1,surf2d2
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) ::  soil2d1,soil2d2

    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: mu2d1,mu2d2,mub2d1,mub2d2,&
         &mu02d1,mu02d2,mapfac_m2d, ter2d,precip2d1, precip2d2
    REAL, ALLOCATABLE, DIMENSION(:,:) :: znu,znw
    REAL, ALLOCATABLE, DIMENSION(:) :: ptop
    REAL, ALLOCATABLE, DIMENSION(:,:) :: dpn

    REAL, ALLOCATABLE, DIMENSION(:) :: dnw,dn,fnp,fnm,nh_term
    REAL :: cof1,cof2,cf1,cf2,cf3,cfn,cfn1
    REAL :: missingVal
    TYPE(time_type)      :: requested_wrf_time

! for now this is only for information
     call set_calendar_type(calendar_type)
     requested_wrf_time = set_date(start_year_f, start_month_f, &
                                    start_day_f, start_hour_f,   &
                                    start_minute_f, 0)

! open forcing file
     ierr = nf_open(init_f_file, 0, ncid)
     IF ( ierr /= NF_NOERR ) THEN
        PRINT*,"Problem opening forcing file, aborting!"
        STOP
     ENDIF

! which index from the WRF file do we want to start from?
    IF (MOD(start_forecast,interval_f) /= 0) THEN
       PRINT *,'Cannot start a forecast at ',start_forecast,'Z'
       PRINT *,'Stopping'
       STOP 'S/R wrf_init_and_bc'
    ENDIF

! use requested sounding or random?
    IF ( rnd_init ) THEN
       nmix = 2
    ELSE
       requested_f_index = -9999
       ALLOCATE(wrf_year(nrecords), wrf_month(nrecords), &
                wrf_day(nrecords), wrf_hour(nrecords))
       nmix = 1
       ierr=nf_inq_varid(ncid,'inityear',fid(13))
       ierr=nf_get_vara_int(ncid,fid(13),(/1,1/),(/1,nrecords/),wrf_year)
       ierr=nf_inq_varid(ncid,'initmonth',fid(14))
       ierr=nf_get_vara_int(ncid,fid(14),(/1,1/),(/1,nrecords/),wrf_month)
       ierr=nf_inq_varid(ncid,'initday',fid(15))
       ierr=nf_get_vara_int(ncid,fid(15),(/1,1/),(/1,nrecords/),wrf_day)
       ierr=nf_inq_varid(ncid,'inithour',fid(16))
       ierr=nf_get_vara_int(ncid,fid(16),(/1,1/),(/1,nrecords/),wrf_hour)

       DO i = 1, nrecords
          IF ( wrf_year(i)  == start_year_f .and. &
               wrf_month(i) == start_month_f .and. &
               wrf_day(i)   == start_day_f .and. &
               wrf_hour(i)  == start_hour_f ) &
               requested_f_index = i
       ENDDO
 
       IF ( requested_f_index < 0 ) THEN
          print*,"Could not find requested sounding date: "
          print*,start_year_f, start_month_f, start_day_f, start_hour_f
          print*,'in file ',init_f_file
          stop 'wrf_init_and_bc'
       ENDIF
          
       print*,'Getting requested WRF forcing valid at date: '
       call print_date(requested_wrf_time)
       print*
       print*,'If you are running DART, the correct time ',&
              ' for the obs_sequence is near'
       call print_time(requested_wrf_time)

    ENDIF
! init the arrays because they will not be completely filled
    t = -9999.
    u = -9999.
    v = -9999.
    q = -9999.
    p = -9999.
    t_ups_x = -9999.
    t_ups_y = -9999.
    q_ups_x = -9999.
    q_ups_y = -9999.
    tau_u = -9999.
    tau_v = -9999.
    th2 = -9999.
    t2 = -9999.
    tsk = -9999.
    u10 = -9999.
    v10 = -9999.
    q2 = -9999.
    precip = -9999.
    glw = -9999.
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

!    PRINT *,'pdimlen(4) =',pdimlen(4)
!    PRINT *,'pdimlen(5) =',pdimlen(5)
    lenp = (/pdimlen(1),pdimlen(2),pdimlen(3),pdimlen(4),1/)
    lenp_x_stag = (/pdimlen_stag(1),pdimlen(2),pdimlen(3),pdimlen(4),1/)
    lenp_y_stag = (/pdimlen(1),pdimlen_stag(2),pdimlen(3),pdimlen(4),1/)
    lenp_z_stag = (/pdimlen(1),pdimlen(2),pdimlen_stag(3),pdimlen(4),1/)

! allocate work arrays
    ALLOCATE(u3d1(pdimlen_stag(1),pdimlen(2),pdimlen(3),pdimlen(4)))
    ALLOCATE(u3d2(pdimlen_stag(1),pdimlen(2),pdimlen(3),pdimlen(4)))
    ALLOCATE(v3d1(pdimlen(1),pdimlen_stag(2),pdimlen(3),pdimlen(4)))
    ALLOCATE(v3d2(pdimlen(1),pdimlen_stag(2),pdimlen(3),pdimlen(4)))
    ALLOCATE(t3d1(pdimlen(1),pdimlen(2),pdimlen(3),pdimlen(4)))
    ALLOCATE(t3d2(pdimlen(1),pdimlen(2),pdimlen(3),pdimlen(4)))
    ALLOCATE(q3d2(pdimlen(1),pdimlen(2),pdimlen(3),pdimlen(4)))
    ALLOCATE(q3d1(pdimlen(1),pdimlen(2),pdimlen(3),pdimlen(4)))
    ALLOCATE(p3d2(pdimlen(1),pdimlen(2),pdimlen(3),pdimlen(4)))
    ALLOCATE(p3d1(pdimlen(1),pdimlen(2),pdimlen(3),pdimlen(4)))
    ALLOCATE(z3d1(pdimlen(1),pdimlen(2),pdimlen_stag(3),pdimlen(4)))
    ALLOCATE(z3d2(pdimlen(1),pdimlen(2),pdimlen_stag(3),pdimlen(4)))
    ALLOCATE(invrho(pdimlen(3),pdimlen(4)))
    ALLOCATE(p_r(pdimlen(3),pdimlen(4)))
    ALLOCATE(p_l(pdimlen(3),pdimlen(4)))
    ALLOCATE(p_u(pdimlen(3),pdimlen(4)))
    ALLOCATE(p_d(pdimlen(3),pdimlen(4)))
    ALLOCATE(z_l(pdimlen(3),pdimlen(4)))
    ALLOCATE(z_m(pdimlen(3),pdimlen(4)))
    ALLOCATE(z_r(pdimlen(3),pdimlen(4)))
    ALLOCATE(z_d(pdimlen(3),pdimlen(4)))
    ALLOCATE(z_u(pdimlen(3),pdimlen(4)))
    ALLOCATE(screen2d1(sdimlen(1),pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(screen2d2(sdimlen(1),pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(surf2d1(sfcdimlen(1),pdimlen(1),pdimlen(2),sfcdimlen(2)))
    ALLOCATE(surf2d2(sfcdimlen(1),pdimlen(1),pdimlen(2),sfcdimlen(2)))
    ALLOCATE(ter2d(pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(soil2d1(sldimlen(1),pdimlen(1),pdimlen(2),sldimlen(2),&
         &sldimlen(3)))
    ALLOCATE(soil2d2(sldimlen(1),pdimlen(1),pdimlen(2),sldimlen(2),&
         &sldimlen(3)))
    ALLOCATE(mu2d1(pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(mu2d2(pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(mub2d1(pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(mub2d2(pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(mu02d1(pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(mu02d2(pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(mapfac_m2d(pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(precip2d1(pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(precip2d2(pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(znu(pdimlen(3),sdimlen(2)))
    ALLOCATE(znw(pdimlen_stag(3),sdimlen(2)))
    ALLOCATE(dpn(pdimlen(3),sdimlen(2)))
    ALLOCATE(dnw(pdimlen(3)))
    ALLOCATE(dn(pdimlen(3)))
    ALLOCATE(fnp(pdimlen(3)))
    ALLOCATE(fnm(pdimlen(3)))
    ALLOCATE(nh_term(sdimlen(2)))

    ALLOCATE(ptop(sdimlen(2)))

! choose a random record and get the data
    IF ( rnd_init ) THEN
       rtran = ran1(idum)*(nrecords-1) + 1.
       itran1 = AINT(rtran)
       rtran = ran1(idum)*(nrecords-1) + 1.
       itran2 = AINT(rtran)
       rtran = ran1(idum)
    ELSE
       itran1 = requested_f_index
       itran2 = -9999
       rtran = 1.0d0
    ENDIF

    control_index(1) = itran1
    control_index(2) = itran2
    control_w = rtran
    PRINT*,"Getting profile from ensembles ",control_index
    PRINT*,"with alpha ", control_w
    PRINT*,"at time index ", wrf_ind

! extra sfc vars and precip

    stter = (/1,1,1,itran1/)
    ierr=nf_get_vara_double(ncid,fid(17),stter,lenter,mu2d1)
    ierr=nf_get_vara_double(ncid,fid(18),stter,lenter,mub2d1)
    ierr=nf_get_vara_double(ncid,fid(19),stter,lenter,mu02d1)
    ierr=nf_get_vara_double(ncid,fid(23),stter,lenter,mapfac_m2d)
    ierr=nf_get_vara_double(ncid,fid(24),stter,lenter,precip2d1)

    IF ( nmix > 1 ) THEN
       stter = (/1,1,1,itran2/)
       ierr=nf_get_vara_double(ncid,fid(17),stter,lenter,mu2d2)
       ierr=nf_get_vara_double(ncid,fid(18),stter,lenter,mub2d2)
       ierr=nf_get_vara_double(ncid,fid(19),stter,lenter,mu02d2)
       ierr=nf_get_vara_double(ncid,fid(24),stter,lenter,precip2d2)
    ENDIF


! 1D height information
    lenvec = (/pdimlen(3),sdimlen(2),1/)
    vec = (/1,1,itran1/)
    ierr=nf_get_vara_double(ncid,fid(20),vec,lenvec,znu)
    lenvec = (/pdimlen_stag(3),sdimlen(2),1/)
    ierr=nf_get_vara_double(ncid,fid(21),vec,lenvec,znw)
    lenvec = (/1,sdimlen(2),1/)
    ierr=nf_get_vara_double(ncid,fid(22),vec,lenvec,ptop)

! 3d vars

    stp = (/1,1,1,1,itran1/)
    ierr=nf_get_vara_double(ncid,fid(1),stp,lenp_x_stag,u3d1)
    ierr=nf_get_vara_double(ncid,fid(2),stp,lenp_y_stag,v3d1)
    ierr=nf_get_vara_double(ncid,fid(3),stp,lenp,t3d1)
    ierr=nf_get_vara_double(ncid,fid(4),stp,lenp,q3d1)
    ierr=nf_get_vara_double(ncid,fid(5),stp,lenp,p3d1)
    ierr=nf_get_vara_double(ncid,fid(6),stp,lenp_z_stag,z3d1)

    IF ( nmix > 1 ) THEN
       stp = (/1,1,1,1,itran2/)
       ierr=nf_get_vara_double(ncid,fid(1),stp,lenp_x_stag,u3d2)
       ierr=nf_get_vara_double(ncid,fid(2),stp,lenp_y_stag,v3d2)
       ierr=nf_get_vara_double(ncid,fid(3),stp,lenp,t3d2)
       ierr=nf_get_vara_double(ncid,fid(4),stp,lenp,q3d2)
       ierr=nf_get_vara_double(ncid,fid(5),stp,lenp,p3d2)
       ierr=nf_get_vara_double(ncid,fid(6),stp,lenp_z_stag,z3d2)
    ENDIF

    ierr=nf_get_vara_double(ncid,fid(12),stter,lenter,ter2d)
    terrain = ter2d(2,2,1)

! do grid differencing and weighing of profiles

    t(:,1:nt) = &
         t3d1(midx,midy,:,wrf_ind:wrf_end_ind)*control_w + &
         t3d2(midx,midy,:,wrf_ind:wrf_end_ind)*(1.0-control_w)

    q(:,1:nt) = &
         q3d1(midx,midy,:,wrf_ind:wrf_end_ind)*control_w + &
         q3d2(midx,midy,:,wrf_ind:wrf_end_ind)*(1.0-control_w)

    p(:,1:nt) = 1.e2*(&
         p3d1(midx,midy,:,wrf_ind:wrf_end_ind)*control_w + &
         p3d2(midx,midy,:,wrf_ind:wrf_end_ind)*(1.0-control_w))

    u(:,1:nt) = &
         .5*(u3d1(midx,midy,:,wrf_ind:wrf_end_ind)+ &
         u3d1(midx+1,midy,:,wrf_ind:wrf_end_ind))*&
         control_w + &
         .5*(u3d2(midx,midy,:,wrf_ind:wrf_end_ind)+&
         u3d2(midx+1,midy,:,wrf_ind:wrf_end_ind))*&
         (1.0-control_w)

    v(:,1:nt) = &
         .5*(v3d1(midx,midy,:,wrf_ind:wrf_end_ind)+&
         v3d1(midx,midy+1,:,wrf_ind:wrf_end_ind))*&
         control_w + &
         .5*(v3d2(midx,midy,:,wrf_ind:wrf_end_ind)+&
         v3d2(midx,midy+1,:,wrf_ind:wrf_end_ind))*&
         (1.0-control_w)

    z_agl_stag(:,1:nt)=MAX(0.,&
         z3d1(midx,midy,:,wrf_ind:wrf_end_ind)*control_w + &
         z3d2(midx,midy,:,wrf_ind:wrf_end_ind)*(1.0-control_w) -&
         ter2d(midx,midy,1))

    DO k=1,nz
       z_agl(k,1:nt)=.5*(&
            (z3d1(midx,midy,k,wrf_ind:wrf_end_ind)+&
            z3d1(midx,midy,k+1,wrf_ind:wrf_end_ind))*control_w + &
            (z3d2(midx,midy,k,wrf_ind:wrf_end_ind)+&
            z3d2(midx,midy,k+1,wrf_ind:wrf_end_ind))*&
            (1.0-control_w))-ter2d(midx,midy,1)
    ENDDO

! upstream values and tau, as defined in Gahn et al. 1999 equation 2
    IF ( t_advection ) then
    t_ups_x(:,1:nt) = &
         t3d1(midx-1,midy,:,wrf_ind:wrf_end_ind)*control_w + &
         t3d2(midx-1,midy,:,wrf_ind:wrf_end_ind)*(1.0-control_w)
    t_ups_y(:,1:nt) = &
         t3d1(midx,midy-1,:,wrf_ind:wrf_end_ind)*control_w + &
         t3d2(midx,midy-1,:,wrf_ind:wrf_end_ind)*(1.0-control_w)
    ENDIF
    IF ( qv_advection ) then
    q_ups_x(:,1:nt) = &
         q3d1(midx-1,midy,:,wrf_ind:wrf_end_ind)*control_w + &
         q3d2(midx-1,midy,:,wrf_ind:wrf_end_ind)*(1.0-control_w)
    q_ups_y(:,1:nt) = &
         q3d1(midx,midy-1,:,wrf_ind:wrf_end_ind)*control_w + &
         q3d2(midx,midy-1,:,wrf_ind:wrf_end_ind)*(1.0-control_w)
    ENDIF
    IF ( t_advection .or. qv_advection ) THEN
    tau_u(:,1:nt) = &
        u3d1(midx,midy,:,wrf_ind:wrf_end_ind)*control_w + &
        u3d2(midx,midy,:,wrf_ind:wrf_end_ind)*(1.0-control_w)
    tau_u = dx / tau_u
    tau_v(:,1:nt) = &
        v3d1(midx,midy,:,wrf_ind:wrf_end_ind)*control_w + &
        v3d2(midx,midy,:,wrf_ind:wrf_end_ind)*(1.0-control_w)
    tau_v = dx / tau_v
    ENDIF

! screen 
    sts = (/1,1,1,1,itran1/)
    ierr=nf_get_vara_double(ncid,fid(7),sts,lens,screen2d1)

    IF ( nmix > 1 ) THEN
       sts = (/1,1,1,1,itran2/)
       ierr=nf_get_vara_double(ncid,fid(7),sts,lens,screen2d2)
    ENDIF

! do grid differencing and weighing of profiles

    t2(1:nt) = &
         &screen2d1(1,midx,midy,wrf_ind:wrf_end_ind)*control_w + &
         &screen2d2(1,midx,midy,wrf_ind:wrf_end_ind)*(1.0-control_w)
    th2 = t2 ! this is a kluge, will want to use P to correct

    u10(1:nt) = &
         &screen2d1(2,midx,midy,wrf_ind:wrf_end_ind)*control_w + &
         screen2d2(2,midx,midy,wrf_ind:wrf_end_ind)*(1.0-control_w)

    v10(1:nt) = &
         &screen2d1(3,midx,midy,wrf_ind:wrf_end_ind)*control_w + &
         screen2d2(3,midx,midy,wrf_ind:wrf_end_ind)*(1.0-control_w)

    q2(1:nt) = &
         &screen2d1(4,midx,midy,wrf_ind:wrf_end_ind)*control_w + &
         screen2d2(4,midx,midy,wrf_ind:wrf_end_ind)*(1.0-control_w)
    
! set all t0 screen values back to missing if they are 0.0
    if ( th2(1) == 0.0 ) th2(1) = -9999
    if ( t2(1) == 0.0 )  t2(1)  = -9999
    if ( u10(1) == 0.0 ) u10(1) = -9999
    if ( v10(1) == 0.0 ) v10(1) = -9999
    if ( q2(1) == 0.0 )  q2(1)  = -9999

! precip is a special case because we need to change from total
! to hourly - the first is zero
    precip(1) = 0.0d0
    DO i = wrf_ind+1,wrf_end_ind
       precip(i-wrf_ind+1) = &
       (precip2d1(midx,midy,i)-precip2d1(midx,midy,i-1))*control_w + &
       (precip2d2(midx,midy,i)-precip2d2(midx,midy,i-1))*(1.0d0-control_w)
    ENDDO

! surface
    stsfc = (/1,1,1,1,itran1/)
    ierr  = nf_get_vara_double(ncid,fid(8),stsfc,lensfc,surf2d1)

    IF ( nmix > 1 ) THEN
       stsfc = (/1,1,1,1,itran2/)
       ierr  = nf_get_vara_double(ncid,fid(8),stsfc,lensfc,surf2d2)
       ierr  = nf_get_att_double(ncid, fid(8), "_FillValue", missingVal)
    ENDIF

! do grid differencing and weighing of profiles
! required surface variables
    tsk(1:nt) = &
         &surf2d1(1,midx,midy,wrf_ind:wrf_end_ind)*control_w + &
         &surf2d2(1,midx,midy,wrf_ind:wrf_end_ind)*(1.0-control_w)

    glw(1:nt) = &
         &surf2d1(2,midx,midy,wrf_ind:wrf_end_ind)*control_w + &
         &surf2d2(2,midx,midy,wrf_ind:wrf_end_ind)*(1.0-control_w)

    gsw(1:nt) = &
         &surf2d1(3,midx,midy,wrf_ind:wrf_end_ind)*control_w + &
         &surf2d2(3,midx,midy,wrf_ind:wrf_end_ind)*(1.0-control_w)
    
    
! check to see if these are available or else default to *_ref values
    IF ( minval(surf2d1(7,midx,midy,wrf_ind:wrf_end_ind)) > missingVal ) then
       qsfc(1:nt) = &
            surf2d1(7,midx,midy,wrf_ind:wrf_end_ind)*control_w + &
            surf2d2(7,midx,midy,wrf_ind:wrf_end_ind)*(1.0-control_w)
    ELSE
       print*,'WARNING: not assigning qsfc'
    ENDIF

! 
    IF ( minval(surf2d1(8,midx,midy,wrf_ind:wrf_end_ind)) > missingVal ) then
       vegfra(1:nt) = &
         surf2d1(8,midx,midy,wrf_ind:wrf_end_ind)*control_w + &
         surf2d2(8,midx,midy,wrf_ind:wrf_end_ind)*(1.0-control_w)
    ELSE
       print*,'WARNING: using vegfra_ref'
       vegfra(1:nt) = vegfra_ref
    ENDIF

    IF ( minval(surf2d1(9,midx,midy,wrf_ind:wrf_end_ind)) > missingVal ) then
       isltyp(1:nt) = NINT(surf2d1(9,midx,midy,wrf_ind:wrf_end_ind))
    ELSE
       print*,'WARNING: using isltyp_ref'
       isltyp(1:nt) = isltyp_ref
    ENDIF

    IF ( minval(surf2d1(10,midx,midy,wrf_ind:wrf_end_ind)) > missingVal ) then
       ivgtyp(1:nt) = NINT(surf2d1(10,midx,midy,wrf_ind:wrf_end_ind))
    ELSE
       print*,'WARNING: using ivgtyp_ref'
       ivgtyp(1:nt) = ivgtyp_ref
    ENDIF

    IF ( minval(surf2d1(11,midx,midy,wrf_ind:wrf_end_ind)) > missingVal ) then
       lu_index(1:nt) = NINT(surf2d1(11,midx,midy,wrf_ind:wrf_end_ind))
    ELSE
       print*,'WARNING: using lu_index_ref'
       lu_index(1:nt) = lu_index_ref
    ENDIF

!soil
! do grid differencing and weighing of profiles
    stsl = (/1,1,1,1,1,itran1/)
    ierr=nf_get_vara_double(ncid,fid(9),stsl,lensl,soil2d1)

    IF ( nmix > 1 ) THEN
       stsl = (/1,1,1,1,1,itran2/)
       ierr=nf_get_vara_double(ncid,fid(9),stsl,lensl,soil2d2)
    ENDIF

    tslb(:,1:nt) = &
         &soil2d1(1,midx,midy,:,wrf_ind:wrf_end_ind)*control_w + &
         &soil2d2(1,midx,midy,:,wrf_ind:wrf_end_ind)*(1.0-control_w)

    smois(:,1:nt) = &
         &soil2d1(2,midx,midy,:,wrf_ind:wrf_end_ind)*control_w + &
         &soil2d2(2,midx,midy,:,wrf_ind:wrf_end_ind)*(1.0-control_w)

    tmn(1:nt) = &
            &tslb(ns,1:sldimlen(2)-wrf_ind+1)

! finally, forecast times 

    DO i=1,nt
       times(i)=(i-1)*interval_f+start_forecast
    ENDDO
    times_flux = times
    times_soil = times
    times_smos = times

    ierr = nf_close(ncid)

    IF ( .not. rnd_init ) DEALLOCATE(wrf_year, wrf_month, wrf_day, wrf_hour)
    DEALLOCATE(u3d1)
    DEALLOCATE(u3d2)
    DEALLOCATE(v3d1)
    DEALLOCATE(v3d2)
    DEALLOCATE(t3d1)
    DEALLOCATE(t3d2)
    DEALLOCATE(q3d2)
    DEALLOCATE(q3d1)
    DEALLOCATE(p3d2)
    DEALLOCATE(p3d1)
    DEALLOCATE(z3d1)
    DEALLOCATE(z3d2)
    DEALLOCATE(invrho)
    DEALLOCATE(p_l)
    DEALLOCATE(p_r)
    DEALLOCATE(p_u)
    DEALLOCATE(p_d)
    DEALLOCATE(z_m)
    DEALLOCATE(z_l)
    DEALLOCATE(z_r)
    DEALLOCATE(z_d)
    DEALLOCATE(z_u)
    DEALLOCATE(screen2d1)
    DEALLOCATE(screen2d2)
    DEALLOCATE(surf2d1)
    DEALLOCATE(surf2d2)
    DEALLOCATE(soil2d1)
    DEALLOCATE(soil2d2)
    DEALLOCATE(ter2d)
    DEALLOCATE(mu2d1)
    DEALLOCATE(mu2d2)
    DEALLOCATE(mub2d1)
    DEALLOCATE(mub2d2)
    DEALLOCATE(mu02d1)
    DEALLOCATE(mu02d2)
    DEALLOCATE(mapfac_m2d)
    DEALLOCATE(precip2d1)
    DEALLOCATE(precip2d2)
 
  END SUBROUTINE wrf_init_and_bc

!
END MODULE module_wrf_init_and_bc

