MODULE module_jh_procedures

  IMPLICIT NONE

CONTAINS

  SUBROUTINE genctrl(bamex,nz,numsoil,ntimes,Nforecasttimes, &
       &Dx,Cent_lat,Cent_lon,Truelat1,Truelat2,Mminlu,Julday,&
       &Sw_corner_lon,Sw_corner_lat,Projcode,&
       &Z_agl,Z_agl_stag,T,U,V,Q,P,&
       &Th2,Tsk,U10,V10,     &
       &Q2,Glw,Gsw,Qsfc,Tslb,Smois,Tmn,&
       &Vegfra,Isltyp,Lu_index,Ivgtyp,&
       &Lat,Lon,Cor,Times,Idum,&
       &start_hour_wrf, hour_interval_wrf, start_forecast,&
       &forecast_length,&
       &Control_index1,Control_index2,Control_w,infilename)

    USE module_nr_procedures
    USE module_model_constants

    INCLUDE 'netcdf.inc'

! Subroutine selects a random forecast from file wrf_cntl.nc

! arguments

    LOGICAL :: bamex
    INTEGER, INTENT(in) :: nz,numsoil,ntimes
    INTEGER, INTENT(out) :: Nforecasttimes
    INTEGER, INTENT(inout) :: idum
    INTEGER, INTENT(out) :: control_index1, control_index2
    REAL, INTENT(out) ::  control_w
    REAL, DIMENSION(:), INTENT(out) :: &
         &th2,tsk,u10,v10,q2,glw,gsw,qsfc,&
         &tmn,vegfra,times
    INTEGER, DIMENSION(:), INTENT(out) :: isltyp,ivgtyp,lu_index

    REAL, DIMENSION(:,:), INTENT(out) :: t,u,v,q,p,tslb,smois
    
    INTEGER, INTENT(in) :: start_hour_wrf, hour_interval_wrf, &
         &start_forecast, forecast_length
    CHARACTER(len=120), INTENT(in) :: infilename
    REAL, INTENT(out) :: lat,lon,cor,dx
    REAL, INTENT(out) :: cent_lat,cent_lon,truelat1,truelat2,&
         &sw_corner_lon,sw_corner_lat
    INTEGER :: julday,projcode
    CHARACTER(len=4) :: mminlu

! local
    INTEGER :: i,k,kkl,kkr,imem, nrecords, itran1, itran2
    INTEGER :: ierr, ncid, timeid
    INTEGER, DIMENSION(23) :: fid  !profiles, screen, surface, soil, z_agl
    INTEGER, DIMENSION(5) :: sdimid, sdimlen, lens, sts
    INTEGER, DIMENSION(5) :: sfcdimid, sfcdimlen, lensfc, stsfc
    INTEGER, DIMENSION(6) :: sldimid, sldimlen, lensl, stsl
    INTEGER, DIMENSION(5) :: pdimid, pdimlen, lenp, stp
    INTEGER, DIMENSION(5) :: pdimid_stag, pdimlen_stag, &
         &lenp_x_stag, lenp_y_stag, lenp_z_stag
    INTEGER, DIMENSION(4) :: stter, lenter
    INTEGER, DIMENSION(3) :: vec, lenvec

    INTEGER :: wrf_ind, wrf_times,midx,midy
    REAL :: rtran, rtran2
    REAL :: t_l,t_m,t_r,dt_dz,dp
    REAL, DIMENSION(:,:) :: z_agl,z_agl_stag
    
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: lat2d,lon2d,ter2d,vegfra2d


    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: isltyp2d,ivgtyp2d,lu_index2d
    
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: u3d1,v3d1,t3d1,q3d1,p3d1, &
         &u3d2, v3d2, t3d2, q3d2, p3d2,z3d1,z3d2
    REAL, ALLOCATABLE, DIMENSION(:,:) :: invrho,p_r,p_l,p_u,p_d
    REAL, ALLOCATABLE, DIMENSION(:,:) :: z_l,z_m,z_r,z_d,z_u


    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: screen2d1, screen2d2,&
         &surf2d1,surf2d2
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) ::  soil2d1,soil2d2

    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: mu2d1,mu2d2,mub2d1,mub2d2,&
         &mu02d1,mu02d2,mapfac_m2d
    REAL, ALLOCATABLE, DIMENSION(:,:) :: znu,znw
    REAL, ALLOCATABLE, DIMENSION(:) :: ptop
    REAL, ALLOCATABLE, DIMENSION(:,:) :: dpn

    REAL, ALLOCATABLE, DIMENSION(:) :: dnw,dn,fnp,fnm,nh_term
    REAL :: cof1,cof2,cf1,cf2,cf3,cfn,cfn1

! which index from the WRF file do we want to start from?
    IF (MOD(start_forecast,hour_interval_wrf) /= 0) THEN
       PRINT *,'Cannot start a forecast at ',start_forecast,'Z'
       PRINT *,'Stopping'
       STOP
    ENDIF

    wrf_ind = NINT(1+(REAL(start_forecast - start_hour_wrf))/&
         &REAL(hour_interval_wrf))

! init the arrays because they will not be completely filled
    t = -9999.
    u = -9999.
    v = -9999.
    q = -9999.
    p = -9999.
    th2 = -9999.
    tsk = -9999.
    u10 = -9999.
    v10 = -9999.
    q2 = -9999.
    glw = -9999.
    gsw = -9999.
    tmn = -9999.
    qsfc = -9999.
    tslb = -9999.
    smois = -9999.
    times = -9999.
    vegfra = -9999.
    lu_index=-9999
    isltyp = -9999
    ivgtyp = -9999

! open forecast file
    ierr = nf_open(infilename, 0, ncid)
    IF ( ierr /= NF_NOERR ) THEN
       PRINT*,"Problem opening wrf_cntl.nc, aborting!"
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
    IF (pdimlen(4) < ntimes) THEN
       PRINT *,'WARNING'
       PRINT *,'Number of output times available from WRF = ',pdimlen(4)
       PRINT *,'is smaller then expected from the namelist ntimes_gen =',&
            &ntimes
    ENDIF

    IF (forecast_length + start_forecast > &
        &(pdimlen(4)-1)*hour_interval_wrf+hour_interval_wrf+start_hour_wrf) THEN
       PRINT *,'Number of output times available from WRF = ',pdimlen(4)
       PRINT *,'Number of output times required to issue ', &
            &forecast_length,'h forecast = ',&
            &NINT(REAL(forecast_length + start_forecast)/&
            &REAL(hour_interval_wrf))+1
       PRINT *,'Stopping'
       STOP
    ENDIF

!    PRINT *,'pdimlen(4) =',pdimlen(4)
    ierr=nf_inq_dimid(ncid,'record',pdimid(5))
    ierr=nf_inq_dimlen(ncid,pdimid(5),pdimlen(5))
!    PRINT *,'pdimlen(5) =',pdimlen(5)
    lenp = (/pdimlen(1),pdimlen(2),pdimlen(3),pdimlen(4),1/)
    lenp_x_stag = (/pdimlen_stag(1),pdimlen(2),pdimlen(3),pdimlen(4),1/)
    lenp_y_stag = (/pdimlen(1),pdimlen_stag(2),pdimlen(3),pdimlen(4),1/)
    lenp_z_stag = (/pdimlen(1),pdimlen(2),pdimlen_stag(3),pdimlen(4),1/)

    midx=(1+pdimlen(1))/2
    midy=(1+pdimlen(2))/2


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
    ALLOCATE(soil2d1(sldimlen(1),pdimlen(1),pdimlen(2),sldimlen(2),&
         &sldimlen(3)))
    ALLOCATE(soil2d2(sldimlen(1),pdimlen(1),pdimlen(2),sldimlen(2),&
         &sldimlen(3)))
    ALLOCATE(lat2d(pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(lon2d(pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(ter2d(pdimlen(1),pdimlen(2),sdimlen(2)))

    ALLOCATE(mu2d1(pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(mu2d2(pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(mub2d1(pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(mub2d2(pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(mu02d1(pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(mu02d2(pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(mapfac_m2d(pdimlen(1),pdimlen(2),sdimlen(2)))
    ALLOCATE(znu(pdimlen(3),sdimlen(2)))
    ALLOCATE(znw(pdimlen_stag(3),sdimlen(2)))
    ALLOCATE(dpn(pdimlen(3),sdimlen(2)))
    ALLOCATE(dnw(pdimlen(3)))
    ALLOCATE(dn(pdimlen(3)))
    ALLOCATE(fnp(pdimlen(3)))
    ALLOCATE(fnm(pdimlen(3)))
    ALLOCATE(nh_term(sdimlen(2)))

    ALLOCATE(ptop(sdimlen(2)))

    IF ( pdimlen(3) /= nz ) THEN
       PRINT*,"wrf_cntl.nc dimensions do not match yours!"
       PRINT *,'nz = ',nz,' pdimlen(3) = ',pdimlen(3)
       STOP
    ENDIF

! choose a random record and get the data
    nrecords = sdimlen(3)
    rtran = ran1(idum)*(nrecords-1) + 1.
    itran1 = AINT(rtran)
    control_index1 = itran1
    rtran = ran1(idum)*(nrecords-1) + 1.
    itran2 = AINT(rtran)
    control_index2 = itran2
    rtran = ran1(idum)
    control_w = rtran
    PRINT*,"Getting profile from ensembles ",control_index1,control_index2
    PRINT*,"with alpha ", control_w
    PRINT*,"at time index ", wrf_ind

! temporary assign begin
    
!    control_w=1.

! temporary assign end

! lat, lon, terrain  

    lenter = (/pdimlen(1),pdimlen(2),sdimlen(2),1/)
    stter = (/1,1,1,itran1/)
    ierr=nf_get_vara_double(ncid,fid(12),stter,lenter,ter2d)
    ierr=nf_get_vara_double(ncid,fid(10),stter,lenter,lat2d)
    ierr=nf_get_vara_double(ncid,fid(11),stter,lenter,lon2d)

    lat=lat2d(midx,midy,1)*DEGRAD
    lon=lon2d(midx,midy,1)*DEGRAD
    cor=2.*eomeg*SIN(lat)

    
! new vars

    stter = (/1,1,1,itran1/)
    ierr=nf_get_vara_double(ncid,fid(17),stter,lenter,mu2d1)
    ierr=nf_get_vara_double(ncid,fid(18),stter,lenter,mub2d1)
    ierr=nf_get_vara_double(ncid,fid(19),stter,lenter,mu02d1)
    ierr=nf_get_vara_double(ncid,fid(23),stter,lenter,mapfac_m2d)

    stter = (/1,1,1,itran2/)
    ierr=nf_get_vara_double(ncid,fid(17),stter,lenter,mu2d2)
    ierr=nf_get_vara_double(ncid,fid(18),stter,lenter,mub2d2)
    ierr=nf_get_vara_double(ncid,fid(19),stter,lenter,mu02d2)


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

    stp = (/1,1,1,1,itran2/)
    ierr=nf_get_vara_double(ncid,fid(1),stp,lenp_x_stag,u3d2)
    ierr=nf_get_vara_double(ncid,fid(2),stp,lenp_y_stag,v3d2)
    ierr=nf_get_vara_double(ncid,fid(3),stp,lenp,t3d2)
    ierr=nf_get_vara_double(ncid,fid(4),stp,lenp,q3d2)
    ierr=nf_get_vara_double(ncid,fid(5),stp,lenp,p3d2)
    ierr=nf_get_vara_double(ncid,fid(6),stp,lenp_z_stag,z3d2)



! do grid differencing and weighing of profiles

    nforecasttimes=pdimlen(4)-wrf_ind+1

    t(:,1:nforecasttimes) = &
         &t3d1(midx,midy,:,wrf_ind:)*control_w + &
         &t3d2(midx,midy,:,wrf_ind:)*(1.0-control_w)

    q(:,1:nforecasttimes) = &
         &q3d1(midx,midy,:,wrf_ind:)*control_w + &
         &q3d2(midx,midy,:,wrf_ind:)*(1.0-control_w)

    p(:,1:nforecasttimes) = 1.e2*(&
         &p3d1(midx,midy,:,wrf_ind:)*control_w + &
         &p3d2(midx,midy,:,wrf_ind:)*(1.0-control_w))

    u(:,1:nforecasttimes) = &
         &.5*(u3d1(midx,midy,:,wrf_ind:)+u3d1(midx+1,midy,:,wrf_ind:))*&
         &control_w + &
         &.5*(u3d2(midx,midy,:,wrf_ind:)+u3d2(midx+1,midy,:,wrf_ind:))*&
         &(1.0-control_w)

    v(:,1:nforecasttimes) = &
         &.5*(v3d1(midx,midy,:,wrf_ind:)+v3d1(midx,midy+1,:,wrf_ind:))*&
         &control_w + &
         &.5*(v3d2(midx,midy,:,wrf_ind:)+v3d2(midx,midy+1,:,wrf_ind:))*&
         &(1.0-control_w)

    z_agl_stag(:,1:nforecasttimes)=MAX(0.,&
         &z3d1(midx,midy,:,wrf_ind:)*control_w + &
         &z3d2(midx,midy,:,wrf_ind:)*(1.0-control_w) -&
         &ter2d(midx,midy,1))

!    DO k=2,nz+1
!       WRITE(350,'(i5,f15.7)')k-1,z_agl_stag(k,nforecasttimes)
!   ENDDO
!    stop

    DO k=1,nz
       z_agl(k,1:nforecasttimes)=.5*(&
            &(z3d1(midx,midy,k,wrf_ind:)+&
            &z3d1(midx,midy,k+1,wrf_ind:))*control_w + &
            &(z3d2(midx,midy,k,wrf_ind:)+&
            &z3d2(midx,midy,k+1,wrf_ind:))*&
            &(1.0-control_w))-ter2d(midx,midy,1)
    ENDDO

! screen 
    sts = (/1,1,1,1,itran1/)
    ierr=nf_get_vara_double(ncid,fid(7),sts,lens,screen2d1)
    sts = (/1,1,1,1,itran2/)
    ierr=nf_get_vara_double(ncid,fid(7),sts,lens,screen2d2)

! do grid differencing and weighing of profiles

    th2(1:sdimlen(2)-wrf_ind+1) = &
         &screen2d1(1,midx,midy,wrf_ind:)*control_w + &
         &screen2d2(1,midx,midy,wrf_ind:)*(1.0-control_w)

    u10(1:sdimlen(2)-wrf_ind+1) = &
         &screen2d1(2,midx,midy,wrf_ind:)*control_w + &
         screen2d2(2,midx,midy,wrf_ind:)*(1.0-control_w)

    v10(1:sdimlen(2)-wrf_ind+1) = &
         &screen2d1(3,midx,midy,wrf_ind:)*control_w + &
         screen2d2(3,midx,midy,wrf_ind:)*(1.0-control_w)

    q2(1:sdimlen(2)-wrf_ind+1) = &
         &screen2d1(4,midx,midy,wrf_ind:)*control_w + &
         screen2d2(4,midx,midy,wrf_ind:)*(1.0-control_w)
    
! surface
    stsfc = (/1,1,1,1,itran1/)
    ierr=nf_get_vara_double(ncid,fid(8),stsfc,lensfc,surf2d1)
    stsfc = (/1,1,1,1,itran2/)
    ierr=nf_get_vara_double(ncid,fid(8),stsfc,lensfc,surf2d2)

! do grid differencing and weighing of profiles

    tsk(1:sfcdimlen(2)-wrf_ind+1) = &
         &surf2d1(1,midx,midy,wrf_ind:)*control_w + &
         &surf2d2(1,midx,midy,wrf_ind:)*(1.0-control_w)

    glw(1:sfcdimlen(2)-wrf_ind+1) = &
         &surf2d1(2,midx,midy,wrf_ind:)*control_w + &
         &surf2d2(2,midx,midy,wrf_ind:)*(1.0-control_w)

    gsw(1:sfcdimlen(2)-wrf_ind+1) = &
         &surf2d1(3,midx,midy,wrf_ind:)*control_w + &
         &surf2d2(3,midx,midy,wrf_ind:)*(1.0-control_w)

    IF (.NOT. bamex) THEN
       tmn(1:sfcdimlen(2)-wrf_ind+1) = &
            &surf2d1(4,midx,midy,wrf_ind:)*control_w + &
            &surf2d2(4,midx,midy,wrf_ind:)*(1.0-control_w)

       qsfc(1:sfcdimlen(2)-wrf_ind+1) = &
            &surf2d1(7,midx,midy,wrf_ind:)*control_w + &
            &surf2d2(7,midx,midy,wrf_ind:)*(1.0-control_w)
    ENDIF


    vegfra(1:sfcdimlen(2)-wrf_ind+1) = &
         &surf2d1(8,midx,midy,wrf_ind:)*control_w + &
         &surf2d2(8,midx,midy,wrf_ind:)*(1.0-control_w)

    isltyp(1:sfcdimlen(2)-wrf_ind+1) = &
         &NINT(surf2d1(9,midx,midy,wrf_ind:))

    ivgtyp(1:sfcdimlen(2)-wrf_ind+1) = &
         &NINT(surf2d1(10,midx,midy,wrf_ind:))

    lu_index(1:sfcdimlen(2)-wrf_ind+1) = &
         &NINT(surf2d1(11,midx,midy,wrf_ind:))


!soil
! do grid differencing and weighing of profiles
    stsl = (/1,1,1,1,1,itran1/)
    ierr=nf_get_vara_double(ncid,fid(9),stsl,lensl,soil2d1)
    stsl = (/1,1,1,1,1,itran2/)
    ierr=nf_get_vara_double(ncid,fid(9),stsl,lensl,soil2d2)

    tslb(:,1:sldimlen(3)-wrf_ind+1) = &
         &soil2d1(1,midx,midy,:,wrf_ind:)*control_w + &
         &soil2d2(1,midx,midy,:,wrf_ind:)*(1.0-control_w)

    smois(:,1:sldimlen(3)-wrf_ind+1) = &
         &soil2d1(2,midx,midy,:,wrf_ind:)*control_w + &
         &soil2d2(2,midx,midy,:,wrf_ind:)*(1.0-control_w)


    IF (bamex) THEN
       tmn(1:sfcdimlen(2)-wrf_ind+1) = &
            &tslb(numsoil,1:sldimlen(2)-wrf_ind+1)
    ENDIF


! finally, forecast times (3 hourly)

    DO i=1,nforecasttimes
       times(i)=(i-1)*hour_interval_wrf+start_forecast
    ENDDO

    ierr = nf_close(ncid)

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
    DEALLOCATE(lat2d)
    DEALLOCATE(lon2d)
    DEALLOCATE(ter2d)

  END SUBROUTINE genctrl

END MODULE module_jh_procedures

