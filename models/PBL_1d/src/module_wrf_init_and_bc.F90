MODULE module_wrf_init_and_bc

  USE time_manager_mod,      only: time_type, GREGORIAN, &
                                   set_calendar_type, print_time, &
                                   print_date, set_date, set_time
  USE module_model_constants
! *** need to add more here **** ???
!11/29/06 DO. Added n_eo
  USE module_namelist,       only: start_hour_f, interval_f, &
                                   start_forecast, forecast_length, &
                                   start_year_f, start_month_f, & 
                                   start_day_f, start_minute_f, &
                                   init_f_file, isltyp_ref, &
                                   ivgtyp_ref, lu_index_ref, vegfra_ref, &
                                   rnd_init, rnd_force, eofs_file_init,&
                                   eofs_file_forc, n_eo, scales,&
                                   t_advection, qv_advection

  IMPLICIT NONE
  private
  INTEGER                ::  midx, midy, wrf_ind, wrf_end_ind, nt_wrfin,nt
  INTEGER                :: nrecords
  INTEGER, PARAMETER     :: num_screen_vars = 4
  INTEGER, PARAMETER     :: num_sfc_vars = 12
  INTEGER, PARAMETER     :: num_soil_vars = 2

  INTEGER, PARAMETER     :: num_screen_dims = 4
  INTEGER, PARAMETER     :: num_sfc_dims = 4
  INTEGER, PARAMETER     :: num_soil_dims = 5
  INTEGER, PARAMETER     :: num_prof_dims = 5
  INTEGER, DIMENSION(38) :: fid  !profiles, screen, surface, soil, z_agl
  INTEGER, DIMENSION(num_screen_dims) :: sdimid, sdimlen, lens, sts
  INTEGER, DIMENSION(num_sfc_dims)    :: sfcdimid, sfcdimlen, lensfc, stsfc
  INTEGER, DIMENSION(num_soil_dims)   :: sldimid, sldimlen, lensl, stsl 
  INTEGER, DIMENSION(num_prof_dims)   :: pdimid, pdimlen, lenp, stp
  INTEGER, DIMENSION(num_prof_dims)   :: pdimid_stag, pdimlen_stag, &
                           lenp_x_stag, lenp_y_stag, lenp_z_stag
  INTEGER, DIMENSION(num_sfc_dims) :: stter, lenter

  INTEGER, DIMENSION(2)               :: control_index

  INTEGER, PARAMETER    :: calendar_type = GREGORIAN

  REAL                                :: control_w

!  DO. Add arrays for eofs
  INTEGER, PARAMETER :: num_prof_eo_dims = 3  ! Times, z, N
  INTEGER, PARAMETER :: num_soil_eo_dims = 3
  INTEGER, PARAMETER :: num_screen_eo_dims = 2 ! TImes, N
  INTEGER, PARAMETER :: num_sfc_eo_dims = 2

  INTEGER, DIMENSION(24) :: fid_eo  
  INTEGER, DIMENSION(num_prof_eo_dims)   :: pdimid_eo, pdimlen_eo
  INTEGER, DIMENSION(num_prof_eo_dims)   :: pdimid_stag_eo, pdimlen_stag_eo
  INTEGER, DIMENSION(num_soil_eo_dims)   :: sldimid_eo, sldimlen_eo
  INTEGER, DIMENSION(num_screen_eo_dims) :: sdimid_eo, sdimlen_eo
  INTEGER, DIMENSION(num_sfc_eo_dims)    :: sfcdimid_eo, sfcdimlen_eo

  INTEGER, DIMENSION(num_screen_eo_dims) :: lens_eo, sts_eo 
  INTEGER, DIMENSION(num_sfc_eo_dims)    :: lensfc_eo, stsfc_eo 
  INTEGER, DIMENSION(num_soil_eo_dims)   :: lensl_eo, stsl_eo, lenp_eo, stp_eo
  INTEGER, DIMENSION(num_sfc_eo_dims)    :: lenev_eo, stev_eo
 
!  DO. Add a subroutine to read from eofs_file, eofs_dims
  public           :: wrf_f_dims, wrf_init_and_bc, eofs_dims

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
! DO the following and vegfra2d never used
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
        PRINT*,"Problem opening forcing file in wrf_init_and_bc, aborting!"
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

! collent all information about profiles
    ierr=nf_inq_dimid(ncid,'x',pdimid(1))
    ierr=nf_inq_dimlen(ncid,pdimid(1),pdimlen(1))
    ierr=nf_inq_dimid(ncid,'y',pdimid(2))
    ierr=nf_inq_dimlen(ncid,pdimid(2),pdimlen(2))
    ierr=nf_inq_dimid(ncid,'z_amsl',pdimid(3))
    ierr=nf_inq_dimlen(ncid,pdimid(3),pdimlen(3))
    nz = pdimlen(3)
    ierr=nf_inq_dimid(ncid,'Times',pdimid(4))
    ierr=nf_inq_dimlen(ncid,pdimid(4),pdimlen(4))
    nt_wrfin = pdimlen(4)
    ierr=nf_inq_dimid(ncid,'record',pdimid(5))
    ierr=nf_inq_dimlen(ncid,pdimid(5),pdimlen(5))
    nrecords = pdimlen(5)

! staggered profiles
    ierr=nf_inq_dimid(ncid,'x_stag',pdimid_stag(1))
    ierr=nf_inq_dimlen(ncid,pdimid_stag(1),pdimlen_stag(1))
    ierr=nf_inq_dimid(ncid,'y_stag',pdimid_stag(2))
    ierr=nf_inq_dimlen(ncid,pdimid_stag(2),pdimlen_stag(2))
    ierr=nf_inq_dimid(ncid,'z_amsl_stag',pdimid_stag(3))
    ierr=nf_inq_dimlen(ncid,pdimid_stag(3),pdimlen_stag(3))
    ierr=nf_inq_dimid(ncid,'Times',pdimid_stag(4))
    ierr=nf_inq_dimlen(ncid,pdimid_stag(4),pdimlen_stag(4))
    ierr=nf_inq_dimid(ncid,'record',pdimid_stag(5))
    ierr=nf_inq_dimlen(ncid,pdimid_stag(5),pdimlen_stag(5))

    ! collect all information about screen
    ierr=nf_inq_dimid(ncid,'x',sdimid(1))
    ierr=nf_inq_dimlen(ncid,sdimid(1),sdimlen(1))
    ierr=nf_inq_dimid(ncid,'y',sdimid(2))
    ierr=nf_inq_dimlen(ncid,sdimid(2),sdimlen(2))
    ierr=nf_inq_dimid(ncid,'Times',sdimid(3))
    ierr=nf_inq_dimlen(ncid,sdimid(3),sdimlen(3))
    ierr=nf_inq_dimid(ncid,'record',sdimid(4))
    ierr=nf_inq_dimlen(ncid,sdimid(4),sdimlen(4))

    ! collect all information about skin (surface)
    ierr=nf_inq_dimid(ncid,'x',sfcdimid(1))
    ierr=nf_inq_dimlen(ncid,sfcdimid(1),sfcdimlen(1))
    ierr=nf_inq_dimid(ncid,'y',sfcdimid(2))
    ierr=nf_inq_dimlen(ncid,sfcdimid(2),sfcdimlen(2))
    ierr=nf_inq_dimid(ncid,'Times',sfcdimid(3))
    ierr=nf_inq_dimlen(ncid,sfcdimid(3),sfcdimlen(3))
    ierr=nf_inq_dimid(ncid,'record',sfcdimid(4))
    ierr=nf_inq_dimlen(ncid,sfcdimid(4),sfcdimlen(4))

    ! collect all information about soil
    ierr=nf_inq_dimid(ncid,'x',sldimid(1))
    ierr=nf_inq_dimlen(ncid,sldimid(1),sldimlen(1))
    ierr=nf_inq_dimid(ncid,'y',sldimid(2))
    ierr=nf_inq_dimlen(ncid,sldimid(2),sldimlen(2))
    ierr=nf_inq_dimid(ncid,'soil_levels',sldimid(3))
    ierr=nf_inq_dimlen(ncid,sldimid(3),sldimlen(3))
    ns = sldimlen(3)
    ierr=nf_inq_dimid(ncid,'Times',sldimid(4))
    ierr=nf_inq_dimlen(ncid,sldimid(4),sldimlen(4))
    ierr=nf_inq_dimid(ncid,'record',sldimid(5))
    ierr=nf_inq_dimlen(ncid,sldimid(5),sldimlen(5))

!   abort if not enough times in the file 
    IF ( wrf_end_ind > nt_wrfin ) THEN
       print*,"Not enough times in the wrf input file for this forecast"
       print*,nt_wrfin - wrf_ind + 1, &
              " forecast times available from ",start_hour_f + &
                start_forecast/interval_f,"Z"
       stop
    ENDIF

! only want a single record each time we read
    lens = (/sdimlen(1:num_screen_dims-1),1/)
    lensfc = (/sfcdimlen(1:num_sfc_dims-1),1/)
    lensl = (/sldimlen(1:num_soil_dims-1),1/)

! variables
    ierr=nf_inq_varid(ncid,'U',fid(1))
    ierr=nf_inq_varid(ncid,'V',fid(2))
    ierr=nf_inq_varid(ncid,'T',fid(3))
    ierr=nf_inq_varid(ncid,'Q',fid(4))
    ierr=nf_inq_varid(ncid,'P',fid(5))
    ierr=nf_inq_varid(ncid,'Z',fid(6))
    ierr=nf_inq_varid(ncid,'T2',fid(7))
    ierr=nf_inq_varid(ncid,'Q2',fid(8))
    ierr=nf_inq_varid(ncid,'U10',fid(9))
    ierr=nf_inq_varid(ncid,'V10',fid(10))
    ierr=nf_inq_varid(ncid,'TSK',fid(11))
    ierr=nf_inq_varid(ncid,'GLW',fid(12))
    ierr=nf_inq_varid(ncid,'GSW',fid(13))
    ierr=nf_inq_varid(ncid,'TMN',fid(14))
    ierr=nf_inq_varid(ncid,'HFX',fid(15))
    ierr=nf_inq_varid(ncid,'QFX',fid(16))
    ierr=nf_inq_varid(ncid,'QSFC',fid(17))
    ierr=nf_inq_varid(ncid,'VEGFRA',fid(18))
    ierr=nf_inq_varid(ncid,'ISLTYP',fid(19))
    ierr=nf_inq_varid(ncid,'IVGTYP',fid(20))
    ierr=nf_inq_varid(ncid,'LU_INDEX',fid(21))
    ierr=nf_inq_varid(ncid,'TSLB',fid(22))
    ierr=nf_inq_varid(ncid,'SMOIS',fid(23))
    ierr=nf_inq_varid(ncid,'PRECIP',fid(24))
    ierr=nf_inq_varid(ncid,'lats',fid(25))
    ierr=nf_inq_varid(ncid,'lons',fid(26))
    ierr=nf_inq_varid(ncid,'terrain',fid(27))
    ierr=nf_inq_varid(ncid,'inityear',fid(28))
    ierr=nf_inq_varid(ncid,'initmonth',fid(29))
    ierr=nf_inq_varid(ncid,'initday',fid(30))
    ierr=nf_inq_varid(ncid,'inithour',fid(31))

    ierr=nf_inq_varid(ncid,'MU',fid(32))
    ierr=nf_inq_varid(ncid,'MUB',fid(33))
    ierr=nf_inq_varid(ncid,'MU0',fid(34))
    ierr=nf_inq_varid(ncid,'ZNU',fid(35))
    ierr=nf_inq_varid(ncid,'ZNW',fid(36))
    ierr=nf_inq_varid(ncid,'P_TOP',fid(37))
    ierr=nf_inq_varid(ncid,'MAPFAC_M',fid(38))

! lat, lon, terrain  

    lenter = (/sfcdimlen(1:num_sfc_dims-1),1/)
    stter = (/1,1,1,1/)

    ALLOCATE(lat2d(sfcdimlen(1),sfcdimlen(2),sfcdimlen(3)))
    ALLOCATE(lon2d(sfcdimlen(1),sfcdimlen(2),sfcdimlen(3)))

    ierr=nf_get_vara_double(ncid,fid(25),stter,lenter,lat2d)
    ierr=nf_get_vara_double(ncid,fid(26),stter,lenter,lon2d)

    midx=(1+sfcdimlen(1))/2
    midy=(1+sfcdimlen(2))/2

    lat=lat2d(midx,midy,1)*DEGRAD
    lon=lon2d(midx,midy,1)*DEGRAD
    cor=2.*eomeg*SIN(lat)

    DEALLOCATE(lat2d)
    DEALLOCATE(lon2d)

    ierr = nf_close(ncid)

  END SUBROUTINE wrf_f_dims

!**********************************************

!  DO. open eofs file and get dimensions
  SUBROUTINE eofs_dims(ncid_eofs_init)

     IMPLICIT NONE
     INCLUDE 'netcdf.inc'

! gets all static data from the input file, including dimensions

     INTEGER, INTENT(inout):: ncid_eofs_init
! local
     INTEGER              :: ierr
     REAL :: missingVal

! open eofs initialization file
     ierr = nf_open(eofs_file_init, 0, ncid_eofs_init)
     IF ( ierr /= NF_NOERR ) THEN
        PRINT*,"Problem opening eofs_init file, aborting!"
        STOP
     ENDIF

! get dimensions for eof profiles
    ierr=nf_inq_dimid(ncid_eofs_init,'N',pdimid_eo(1))
    ierr=nf_inq_dimlen(ncid_eofs_init,pdimid_eo(1),pdimlen_eo(1))
    ierr=nf_inq_dimid(ncid_eofs_init,'z_amsl',pdimid_eo(2))
    ierr=nf_inq_dimlen(ncid_eofs_init,pdimid_eo(2),pdimlen_eo(2))
    ierr=nf_inq_dimid(ncid_eofs_init,'Times',pdimid_eo(3))
    ierr=nf_inq_dimlen(ncid_eofs_init,pdimid_eo(3),pdimlen_eo(3))

! get dimensions for eof staggered profiles
    ierr=nf_inq_dimid(ncid_eofs_init,'N',pdimid_stag_eo(1))
    ierr=nf_inq_dimlen(ncid_eofs_init,pdimid_stag_eo(1),pdimlen_stag_eo(1))
    ierr=nf_inq_dimid(ncid_eofs_init,'z_amsl_stag',pdimid_stag_eo(2))
    ierr=nf_inq_dimlen(ncid_eofs_init,pdimid_stag_eo(2),pdimlen_stag_eo(2))
    ierr=nf_inq_dimid(ncid_eofs_init,'Times',pdimid_stag_eo(3))
    ierr=nf_inq_dimlen(ncid_eofs_init,pdimid_stag_eo(3),pdimlen_stag_eo(3))

! get dimensions for eof soil
    ierr=nf_inq_dimid(ncid_eofs_init,'N',sldimid_eo(1))
    ierr=nf_inq_dimlen(ncid_eofs_init,sldimid_eo(1),sldimlen_eo(1))
    ierr=nf_inq_dimid(ncid_eofs_init,'soil_levels',sldimid_eo(2))
    ierr=nf_inq_dimlen(ncid_eofs_init,sldimid_eo(2),sldimlen_eo(2))
    ierr=nf_inq_dimid(ncid_eofs_init,'Times',sldimid_eo(3))
    ierr=nf_inq_dimlen(ncid_eofs_init,sldimid_eo(3),sldimlen_eo(3))

! get dimensions for eof screen
    ierr=nf_inq_dimid(ncid_eofs_init,'N',sdimid_eo(1))
    ierr=nf_inq_dimlen(ncid_eofs_init,sdimid_eo(1),sdimlen_eo(1))
    ierr=nf_inq_dimid(ncid_eofs_init,'Times',sdimid_eo(2))
    ierr=nf_inq_dimlen(ncid_eofs_init,sdimid_eo(2),sdimlen_eo(2))

! get dimensions for eof surface
    ierr=nf_inq_dimid(ncid_eofs_init,'N',sfcdimid_eo(1))
    ierr=nf_inq_dimlen(ncid_eofs_init,sfcdimid_eo(1),sfcdimlen_eo(1))
    ierr=nf_inq_dimid(ncid_eofs_init,'Times',sfcdimid_eo(2))
    ierr=nf_inq_dimlen(ncid_eofs_init,sfcdimid_eo(2),sfcdimlen_eo(2))

! variables for init
    ierr=nf_inq_varid(ncid_eofs_init,'U',fid_eo(1))
    ierr=nf_inq_varid(ncid_eofs_init,'V',fid_eo(2))
    ierr=nf_inq_varid(ncid_eofs_init,'T',fid_eo(3))
    ierr=nf_inq_varid(ncid_eofs_init,'Q',fid_eo(4))
    ierr=nf_inq_varid(ncid_eofs_init,'P',fid_eo(5))
    ierr=nf_inq_varid(ncid_eofs_init,'Z',fid_eo(6))
    ierr=nf_inq_varid(ncid_eofs_init,'T2',fid_eo(7))
    ierr=nf_inq_varid(ncid_eofs_init,'Q2',fid_eo(8))
    ierr=nf_inq_varid(ncid_eofs_init,'U10',fid_eo(9))
    ierr=nf_inq_varid(ncid_eofs_init,'V10',fid_eo(10))
    ierr=nf_inq_varid(ncid_eofs_init,'TSK',fid_eo(11))
    ierr=nf_inq_varid(ncid_eofs_init,'GLW',fid_eo(12))
    ierr=nf_inq_varid(ncid_eofs_init,'GSW',fid_eo(13))
    ierr=nf_inq_varid(ncid_eofs_init,'TMN',fid_eo(14))
    ierr=nf_inq_varid(ncid_eofs_init,'HFX',fid_eo(15))
    ierr=nf_inq_varid(ncid_eofs_init,'QFX',fid_eo(16))
    ierr=nf_inq_varid(ncid_eofs_init,'QSFC',fid_eo(17))
    ierr=nf_inq_varid(ncid_eofs_init,'VEGFRA',fid_eo(18))
    ierr=nf_inq_varid(ncid_eofs_init,'ISLTYP',fid_eo(19))
    ierr=nf_inq_varid(ncid_eofs_init,'IVGTYP',fid_eo(20))
    ierr=nf_inq_varid(ncid_eofs_init,'LU_INDEX',fid_eo(21))
    ierr=nf_inq_varid(ncid_eofs_init,'TSLB',fid_eo(22))
    ierr=nf_inq_varid(ncid_eofs_init,'SMOIS',fid_eo(23))
    ierr=nf_inq_varid(ncid_eofs_init,'EIGENVALS',fid_eo(24))
   
   ierr = nf_close(ncid_eofs_init)
 
  END SUBROUTINE eofs_dims

!**********************************************


  SUBROUTINE wrf_init_and_bc(ncid,nz, ns, nt, &
       z_agl,z_agl_stag,t,u,v,q,p,&
       t_ups_x, t_ups_y, &
       q_ups_x, q_ups_y, &
       tau_u, tau_v,&
       th2,t2,tsk,u10,v10,     &
       q2,glw,gsw,qsfc,tslb,smois,tmn,&
       precip, &
       vegfra,isltyp,lu_index,ivgtyp, terrain, dx, &
       times,times_flux,times_soil,times_smos,idum, &
       itran1, itran2, control_w, gasdo)

    USE module_nr_procedures

    IMPLICIT NONE
    INCLUDE 'netcdf.inc'


! arguments
! DO I add here to pass to uvg routine
    INTEGER, INTENT(out) :: itran1, itran2

! I add here to have the ncid of eofs
    INTEGER :: ncid_eofs_init

    INTEGER, INTENT(inout) :: ncid 
    INTEGER, INTENT(in)  :: nz, ns, nt
    INTEGER, INTENT(inout) :: idum
    REAL, DIMENSION(:), INTENT(out) :: &
         th2,t2,tsk,u10,v10,q2,glw,gsw,qsfc, precip, &
         tmn,vegfra,times,times_flux, times_soil, times_smos
    INTEGER, DIMENSION(:), INTENT(out) :: isltyp,ivgtyp,lu_index
    REAL, INTENT(out)                  :: terrain
    REAL, INTENT( in)                  :: dx

!DO I add here to pass to uvg routine *****
    REAL, INTENT(out) :: control_w
    REAL, DIMENSION(:), INTENT(out) :: gasdo

    REAL, DIMENSION(:,:), INTENT(out) :: t,u,v,q,p,tslb,smois
    REAL, DIMENSION(:,:), INTENT(out)   :: tau_u, tau_v
    REAL, DIMENSION(:,:), INTENT(out) :: t_ups_x, t_ups_y
    REAL, DIMENSION(:,:), INTENT(out) :: q_ups_x, q_ups_y
    
! local
    INTEGER :: i,k,kkl,kkr,imem
    INTEGER :: ierr, nmix, imix
    INTEGER, DIMENSION(3) :: vec, lenvec

    INTEGER                           :: requested_f_index
    INTEGER, DIMENSION(:), ALLOCATABLE:: wrf_year, wrf_month, wrf_day, wrf_hour
    REAL :: rtran, rtran2
    REAL :: t_l,t_m,t_r,dt_dz,dp
    REAL, DIMENSION(:,:) :: z_agl,z_agl_stag
    
   
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: u3d1,v3d1,t3d1,q3d1,p3d1, &
          u3d2, v3d2, t3d2, q3d2, p3d2,z3d1,z3d2
    REAL, ALLOCATABLE, DIMENSION(:,:) :: invrho,p_r,p_l,p_u,p_d
    REAL, ALLOCATABLE, DIMENSION(:,:) :: z_l,z_m,z_r,z_d,z_u


    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: screen2d1, screen2d2,&
          surf2d1,surf2d2
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) ::  soil2d1,soil2d2

    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: mu2d1,mu2d2,mub2d1,mub2d2,&
          mu02d1,mu02d2,mapfac_m2d, ter2d, precip2d1, precip2d2
!    REAL, ALLOCATABLE, DIMENSION(:,:) :: znu,znw
!    REAL, ALLOCATABLE, DIMENSION(:) :: ptop
!    REAL, ALLOCATABLE, DIMENSION(:,:) :: dpn

    REAL, ALLOCATABLE, DIMENSION(:) :: dnw,dn,fnp,fnm,nh_term
    REAL :: cof1,cof2,cf1,cf2,cf3,cfn,cfn1
    REAL :: missingVal, missingVal_gsw_eo
    TYPE(time_type)      :: requested_wrf_time

! DO. add arrays for eofs of profiles
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: un, vn, tn, qn, pn, zn, ugn, vgn
! for screen variables one of the index denotes variable (t2, u10, v10, q2), for surf vars too 
! (tsk, glw, gsw, qsfc)
    REAL,  ALLOCATABLE, DIMENSION(:,:,:) :: screenn, surfn
! for soil there are vert levels too, one of the index is var (tslb, smois)
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) ::  soiln

! now declare the eigenvalues
    REAL, ALLOCATABLE, DIMENSION(:,:) :: eval


! for the eof based perturbations
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: pu, pv, pt, pq, pp, pz
    REAL, ALLOCATABLE, DIMENSION(:,:) :: put, pvt, ptt, pqt, ppt, pzt
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: pscr, psur
    REAL, ALLOCATABLE, DIMENSION(:,:) :: pscrt, psurt
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: psl
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: pslt
    REAL, ALLOCATABLE, DIMENSION(:,:) :: auxfac
    REAL, ALLOCATABLE, DIMENSION(:,:) :: damped_scale
    INTEGER :: ivar,iz,islev,ieo

! for now this is only for information
     call set_calendar_type(calendar_type)
     requested_wrf_time = set_date(start_year_f, start_month_f, &
                                    start_day_f, start_hour_f,   &
                                    start_minute_f, 0)

! open forcing file
     ierr = nf_open(init_f_file, 0, ncid)
     IF ( ierr /= NF_NOERR ) THEN
        PRINT*,"Problem opening forcing file in wrf_init_and_bc (2), aborting!"
        STOP
     ENDIF

! which index from the WRF file do we want to start from?
    IF (MOD(start_forecast,interval_f) /= 0) THEN
       PRINT *,'Cannot start a forecast at ',start_forecast,' seconds'
       PRINT *,'Stopping'
       STOP 'S/R wrf_init_and_bc'
    ENDIF

! use requested sounding or random?
! DO. rnd_init changed to integer. 1=undated random mixture, 2=dated perfect mode, 3=dated eofs.
    IF ( rnd_init /= 1) THEN
! DO this option means that we chose dated initialization    
       requested_f_index = -9999
       ALLOCATE(wrf_year(nrecords), wrf_month(nrecords), &
                wrf_day(nrecords), wrf_hour(nrecords))
!DO eliminated nmix = 1
       ierr=nf_get_vara_int(ncid,fid(28),(/1,1/),(/1,nrecords/),wrf_year)
       ierr=nf_get_vara_int(ncid,fid(29),(/1,1/),(/1,nrecords/),wrf_month)
       ierr=nf_get_vara_int(ncid,fid(30),(/1,1/),(/1,nrecords/),wrf_day)
       ierr=nf_get_vara_int(ncid,fid(31),(/1,1/),(/1,nrecords/),wrf_hour)

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

    lenp = (/pdimlen(1:num_prof_dims-1),1/)
    ! need to take individual staggered dimensions
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
    ALLOCATE(screen2d1(num_screen_vars,sdimlen(1),sdimlen(2),sdimlen(3)))
    ALLOCATE(screen2d2(num_screen_vars,sdimlen(1),sdimlen(2),sdimlen(3)))
    ALLOCATE(surf2d1(num_sfc_vars,sfcdimlen(1),sfcdimlen(2),sfcdimlen(3)))
    ALLOCATE(surf2d2(num_sfc_vars,sfcdimlen(1),sfcdimlen(2),sfcdimlen(3)))
    ALLOCATE(ter2d(sfcdimlen(1),sfcdimlen(2),sfcdimlen(3)))
    ALLOCATE(soil2d1(num_soil_vars,sldimlen(1),sldimlen(2),sldimlen(3),sldimlen(4)))
    ALLOCATE(soil2d2(num_soil_vars,sldimlen(1),sldimlen(2),sldimlen(3),sldimlen(4)))
    ALLOCATE(mu2d1(sfcdimlen(1),sfcdimlen(2),sfcdimlen(3)))
    ALLOCATE(mu2d2(sfcdimlen(1),sfcdimlen(2),sfcdimlen(3)))
    ALLOCATE(mub2d1(sfcdimlen(1),sfcdimlen(2),sfcdimlen(3)))
    ALLOCATE(mub2d2(sfcdimlen(1),sfcdimlen(2),sfcdimlen(3)))
    ALLOCATE(mu02d1(sfcdimlen(1),sfcdimlen(2),sfcdimlen(3)))
    ALLOCATE(mu02d2(sfcdimlen(1),sfcdimlen(2),sfcdimlen(3)))
    ALLOCATE(mapfac_m2d(sfcdimlen(1),sfcdimlen(2),sfcdimlen(3)))
    ALLOCATE(precip2d1(sfcdimlen(1),sfcdimlen(2),sfcdimlen(3)))
    ALLOCATE(precip2d2(sfcdimlen(1),sfcdimlen(2),sfcdimlen(3)))
!    ALLOCATE(znu(pdimlen(3),pdimlen(4)))
!    ALLOCATE(znw(pdimlen_stag(3),pdimlen(4)))
!    ALLOCATE(dpn(pdimlen(3),pdimlen(4)))
!    ALLOCATE(dnw(pdimlen(3)))
!    ALLOCATE(dn(pdimlen(3)))
!    ALLOCATE(fnp(pdimlen(3)))
!    ALLOCATE(fnm(pdimlen(3)))
!    ALLOCATE(nh_term(pdimlen(4)))

!    ALLOCATE(ptop(pdimlen(4)))

! DO Allocate here the arrays for the perturbations. Dimensions have to be
! compatible with dimensions from WRF profiles, screen, surface and soil
! variables. I can allocate them even if I am
! not going to use it, not going to read eof file, because the dimensions are
! known to me anyway.

! profiles perturbations
! the dimensions are: # of vertical levels, time
     ALLOCATE(damped_scale(pdimlen_stag(3),pdimlen(4)))
     ALLOCATE(put(pdimlen(3),pdimlen(4)))
     ALLOCATE(pvt(pdimlen(3),pdimlen(4)))
     ALLOCATE(ptt(pdimlen(3),pdimlen(4)))
     ALLOCATE(pqt(pdimlen(3),pdimlen(4)))
     ALLOCATE(ppt(pdimlen(3),pdimlen(4)))
     ALLOCATE(pzt(pdimlen_stag(3),pdimlen(4)))

! screen variables
! the dimensions are: variable, time
     ALLOCATE(pscrt(num_screen_vars,sdimlen(3))) 

! surface variables
! the dimensions are: variable, time
     ALLOCATE(psurt(num_sfc_vars,sfcdimlen(3))) !surface

! soil variables
! the dimensions are: variable, # of soil levels, time
     ALLOCATE(pslt(num_soil_vars,sldimlen(3),sldimlen(4))) !soil

! precipitation - not used yet

! DO Zero the total perturbations arrays. This way they will not contribute in
! the other cases (perfect model and random profile)

    put=0
    pvt=0
    ptt=0
    pqt=0
    ppt=0
    pzt=0
    pscrt=0
    psurt=0
    pslt=0

    IF ( rnd_init == 3 ) THEN
      call eofs_dims(ncid_eofs_init)

! open eofs initialization file
    ierr = nf_open(eofs_file_init, 0, ncid_eofs_init)

    IF ( ierr /= NF_NOERR ) THEN
        PRINT*,"Problem opening eofs_init file, aborting!"
        STOP
    ENDIF

! Check dimensions compatibility between WRF initialization variables and eofs

! Times
     IF (pdimlen(4).ne.pdimlen_eo(3)) THEN
        PRINT*,"Times in WRF init file and in eofs init file differ,&
        aborting"
        STOP
     ENDIF

! Vertical levels, z_amsl
     IF (pdimlen(3).ne.pdimlen_eo(2)) THEN
        PRINT*,"# of vertical levels z_amsl in WRF and eofs differ, aborting"
        STOP
     ENDIF
! Vertical levels, staggered, z_amsl_stag
     IF (pdimlen_stag(3).ne.pdimlen_stag_eo(2)) THEN
        PRINT*,"# of staggered vertical levels z_amsl in WRF and eofs differ, aborting"
        STOP
     ENDIF
! Soil levels, soil_levels
     IF (sldimlen(3).ne.sldimlen_eo(2)) THEN
        PRINT*,"# of soil levels soil_levels in WRF and eofs differ, aborting"
        PRINT*,"sldimlen(3),sldimlen_eo(2)",sldimlen(3),sldimlen_eo(2)
        STOP
     ENDIF

! Here allocate working arrays    
! profiles
! 3rd dimension is time, 2nd is vert levels, 1st is number of eofs
! I WILL USE THE VARIABLES THAT DEFINE LEVELS, VARIABLES, TIME FOR WRF, BECAUSE
! IF THERE IS IMCOMPATIBILITY I SHOULD HAVE STOPPED ANYWAY
! I allocate according to n_eo, the # of eofs that I will use, this way
! I can work on the whole array without risk of undefined value, etc.

! profiles, 1st dimension is number of eofs, 2nd is vertical levels, 3rd is
! times
     ALLOCATE(un(n_eo,pdimlen(3),pdimlen(4)))
     ALLOCATE(vn(n_eo,pdimlen(3),pdimlen(4)))
     ALLOCATE(tn(n_eo,pdimlen(3),pdimlen(4)))
     ALLOCATE(qn(n_eo,pdimlen(3),pdimlen(4)))
     ALLOCATE(pn(n_eo,pdimlen(3),pdimlen(4)))
! zn is in staggered vertical levels     
     ALLOCATE(zn(n_eo,pdimlen_stag(3),pdimlen(4)))

! 1st dimension is variable, 2nd is number of eofs, 3rd is time
     ALLOCATE(screenn(num_screen_vars,n_eo,sdimlen(3))) !screen
     ALLOCATE(surfn(num_sfc_vars,n_eo,sfcdimlen(3))) !surface
! 1st dimension is variable, 2nd is number of eofs, 3rd is soil levels, 4th is time
     ALLOCATE(soiln(num_soil_vars,n_eo,sldimlen(3),sldimlen(4))) !soil
! eigenvalues
! 1st dimension is number of eofs, 2nd is time
     ALLOCATE(eval(n_eo,pdimlen(4))) !for init
! auxiliar working array
     ALLOCATE(auxfac(n_eo,pdimlen(4))) 

! to allocate the perturbations
! profiles
! for each mode
! the dimensions are: # of eofs, # of vertical levels, time
     ALLOCATE(pu(n_eo,pdimlen(3),pdimlen(4)))
     ALLOCATE(pv(n_eo,pdimlen(3),pdimlen(4)))
     ALLOCATE(pt(n_eo,pdimlen(3),pdimlen(4)))
     ALLOCATE(pq(n_eo,pdimlen(3),pdimlen(4)))
     ALLOCATE(pp(n_eo,pdimlen(3),pdimlen(4)))
! pz is in staggered vertical levels      
     ALLOCATE(pz(n_eo,pdimlen_stag(3),pdimlen(4)))

! screen variables, for each mode
! the dimensions are: variable, # of eofs, time
     ALLOCATE(pscr(num_screen_vars,n_eo,sdimlen(3))) !screen
     ALLOCATE(psur(num_sfc_vars,n_eo,sfcdimlen(3))) !surface

! soil variables, for each mode
! the dimensions are: variable, # of eofs, # of soil levels, time
     ALLOCATE(psl(num_soil_vars,n_eo,sldimlen(3),sldimlen(4))) !soil

   ENDIF ! eof_init (rnd_init == 3)

!  DO. Changed for 3 cases
! rnd_init=1 is random climo, rnd_init=2 is for defined profile (perfect model)
! rnd_init=3 is for eofs perturbation
! rnd_init=1 ===> random climo ===> control_w between 0 and 1
! rnd_init=2 ===> perfect model ===> control_w=1, and eofs perturbations stay zero
! rnd_init=3 ===> eofs perturb ===> control_w=1, and eofs pert different from zero
    IF ( rnd_init == 1) THEN
! this is used for random climo       
       rtran = ran1(idum)*(nrecords-1) + 1.
       itran1 = AINT(rtran)
       rtran = ran1(idum)*(nrecords-1) + 1.
       itran2 = AINT(rtran)
       rtran = ran1(idum)
       control_w = rtran
    ELSE 
! this is used for perfect model or eofs perturb
       itran1 = requested_f_index
       itran2 = -9999
       rtran = 1.0d0
       control_w = rtran
    ENDIF

    control_index(1) = itran1
    control_index(2) = itran2
    control_w = rtran
    PRINT*,"Getting profile from ensembles ",control_index
    PRINT*,"with alpha ", control_w
    PRINT*,"at time index ", wrf_ind
! extra sfc vars and precip

    stsfc = (/1,1,1,itran1/)
    ierr=nf_get_vara_double(ncid,fid(32),stsfc,lensfc,mu2d1)
    ierr=nf_get_vara_double(ncid,fid(33),stsfc,lensfc,mub2d1)
    ierr=nf_get_vara_double(ncid,fid(34),stsfc,lensfc,mu02d1)
    ierr=nf_get_vara_double(ncid,fid(38),stsfc,lensfc,mapfac_m2d)
    ierr=nf_get_vara_double(ncid,fid(24),stsfc,lensfc,precip2d1)

    IF ( rnd_init == 1 ) THEN
       stter = (/1,1,1,itran2/)
       ierr=nf_get_vara_double(ncid,fid(32),stsfc,lensfc,mu2d2)
       ierr=nf_get_vara_double(ncid,fid(33),stsfc,lensfc,mub2d2)
       ierr=nf_get_vara_double(ncid,fid(34),stsfc,lensfc,mu02d2)
       ierr=nf_get_vara_double(ncid,fid(24),stsfc,lensfc,precip2d2)
    ENDIF

! 1D height information
!    lenvec = (/pdimlen(3),pdimlen(4),1/)
!    vec = (/1,1,itran1/)
!    ierr=nf_get_vara_double(ncid,fid(35),vec,lenvec,znu)
!    lenvec = (/pdimlen_stag(3),pdimlen(4),1/)
!    ierr=nf_get_vara_double(ncid,fid(36),vec,lenvec,znw)
!    lenvec = (/1,pdimlen(4),1/)
!    ierr=nf_get_vara_double(ncid,fid(37),vec,lenvec,ptop)

! 3d vars

    stp = (/1,1,1,1,itran1/)
    ierr=nf_get_vara_double(ncid,fid(1),stp,lenp_x_stag,u3d1)
    ierr=nf_get_vara_double(ncid,fid(2),stp,lenp_y_stag,v3d1)
    ierr=nf_get_vara_double(ncid,fid(3),stp,lenp,t3d1)
    ierr=nf_get_vara_double(ncid,fid(4),stp,lenp,q3d1)
    ierr=nf_get_vara_double(ncid,fid(5),stp,lenp,p3d1)
    ierr=nf_get_vara_double(ncid,fid(6),stp,lenp_z_stag,z3d1)
 
    IF ( rnd_init == 1 ) THEN
       stp = (/1,1,1,1,itran2/)
       ierr=nf_get_vara_double(ncid,fid(1),stp,lenp_x_stag,u3d2)
       ierr=nf_get_vara_double(ncid,fid(2),stp,lenp_y_stag,v3d2)
       ierr=nf_get_vara_double(ncid,fid(3),stp,lenp,t3d2)
       ierr=nf_get_vara_double(ncid,fid(4),stp,lenp,q3d2)
       ierr=nf_get_vara_double(ncid,fid(5),stp,lenp,p3d2)
       ierr=nf_get_vara_double(ncid,fid(6),stp,lenp_z_stag,z3d2)
    ENDIF
 
    stter = (/1,1,1,itran1/)
    ierr=nf_get_vara_double(ncid,fid(27),stter,lenter,ter2d)
    terrain = ter2d(midx,midy,1)

! DO. read eofs
    IF ( rnd_init == 3 ) THEN

! I only read the n_eo eofs that I am going to use
! For initialization

! eofs of profiles
      stp_eo = (/1,1,1/)
      lenp_eo = (/n_eo,pdimlen(3:4)/)
      ierr=nf_get_vara_double(ncid_eofs_init,fid_eo(1),stp_eo,lenp_eo,un)
      ierr=nf_get_vara_double(ncid_eofs_init,fid_eo(2),stp_eo,lenp_eo,vn)
      ierr=nf_get_vara_double(ncid_eofs_init,fid_eo(3),stp_eo,lenp_eo,tn)
      ierr=nf_get_vara_double(ncid_eofs_init,fid_eo(4),stp_eo,lenp_eo,qn)
      ierr=nf_get_vara_double(ncid_eofs_init,fid_eo(5),stp_eo,lenp_eo,pn)
! z eofs are in staggered vertical levels      
      lenp_eo = (/n_eo,pdimlen_stag(3),pdimlen(4)/)
      ierr=nf_get_vara_double(ncid_eofs_init,fid_eo(6),stp_eo,lenp_eo,zn)

! eofs of screen variables
! missing values: there are no screen values for the analyses
      sts_eo = (/1,1/)
      lens_eo = (/n_eo,sdimlen(3)/)
      ierr=nf_get_vara_double(ncid_eofs_init,fid_eo(7),sts_eo,lens_eo,screenn(1,:,:))  !T2
      ierr=nf_get_vara_double(ncid_eofs_init,fid_eo(8),sts_eo,lens_eo,screenn(4,:,:))  !Q2 
      ierr=nf_get_vara_double(ncid_eofs_init,fid_eo(9),sts_eo,lens_eo,screenn(2,:,:))  !U10
      ierr=nf_get_vara_double(ncid_eofs_init,fid_eo(10),sts_eo,lens_eo,screenn(3,:,:)) !V10

! eofs of surface variables
      stsfc_eo = (/1,1/)
      lensfc_eo = (/n_eo,sfcdimlen(3)/)
      ierr=nf_get_vara_double(ncid_eofs_init,fid_eo(11),sts_eo,lensfc_eo,surfn(1,:,:)) !TSK
      ierr=nf_get_vara_double(ncid_eofs_init,fid_eo(12),sts_eo,lensfc_eo,surfn(2,:,:)) !GLW
      ierr=nf_get_vara_double(ncid_eofs_init,fid_eo(13),sts_eo,lensfc_eo,surfn(3,:,:)) !GSW
      ierr=nf_get_vara_double(ncid_eofs_init,fid_eo(14),sts_eo,lensfc_eo,surfn(7,:,:)) !TMN

! eofs of soil variables
      stsl_eo = (/1,1,1/)
      lensl_eo = (/n_eo,sldimlen(3:4)/)
      ierr=nf_get_vara_double(ncid_eofs_init,fid_eo(22),stsl_eo,lensl_eo,soiln(1,:,:,:)) !TSLB
      ierr=nf_get_vara_double(ncid_eofs_init,fid_eo(23),stsl_eo,lensl_eo,soiln(2,:,:,:)) !SMOIS

! read eigenvalues
      stev_eo = (/1,1/)
      lenev_eo = (/n_eo,sfcdimlen(3)/)
      ierr=nf_get_vara_double(ncid_eofs_init,fid_eo(24),stev_eo,lenev_eo,eval)

! get missing values for gsw, replace with zeros
      ierr  = nf_get_att_double(ncid, fid(13), "_FillValue", missingVal_gsw_eo)
      WHERE(surfn(3,:,:) == missingVal_gsw_eo) surfn(3,:,:) = 0.0

! the following lines are needed to get a gaussian distributed random number 
     DO ieo=1,n_eo
        gasdo(ieo)=gasdev(idum)
! auxfac is the weight of the perturbation
        auxfac(ieo,:)=sqrt(eval(ieo,:))*gasdo(ieo)
     ENDDO

! calculate first the term for each eof
! profiles
    do iz=1,nz
     pu(:,iz,:)=un(:,iz,:)*auxfac(:,:)
     pv(:,iz,:)=vn(:,iz,:)*auxfac(:,:)
     pt(:,iz,:)=tn(:,iz,:)*auxfac(:,:)
     pq(:,iz,:)=qn(:,iz,:)*auxfac(:,:)
     pp(:,iz,:)=pn(:,iz,:)*auxfac(:,:)
    enddo
     
! z is in staggered vertical levels
! I am avoiding the first level because it has only missing values. Easiest
! though not smartest way to do it    
    do iz=2,nz+1
     pz(:,iz,:)=zn(:,iz,:)*auxfac(:,:)
    enddo
     
!  PRINT*,'zn init',zn
     
! screen 
    do ivar=1,num_screen_vars
     pscr(ivar,:,:)=screenn(ivar,:,:)*auxfac(:,:)
    enddo

! surface
    do ivar=1,num_sfc_vars
     psur(ivar,:,:)=surfn(ivar,:,:)*auxfac(:,:)
    enddo

! soil
    do ivar=1,2
     do islev=1,ns
      psl(ivar,:,islev,:)=soiln(ivar,:,islev,:)*auxfac(:,:)
     enddo
    enddo

! damp the scales so that forcing spread is 0 at top
    do k = 1,pdimlen_stag(3)
      do i = 1, pdimlen(4)
        damped_scale(k,i) = scales * cos(dble(k)/dble(pdimlen_stag(3))*acos(-1.0)/2)**2
      enddo
    enddo
    damped_scale = scales ! replace with constant

! calculate the perturbation over all eofs
    do ieo=1,n_eo
! profiles perturbations of time and vertical level
     put(:,:)=put(:,:)+pu(ieo,:,:)*damped_scale(1:pdimlen(3),:)
     pvt(:,:)=pvt(:,:)+pv(ieo,:,:)*damped_scale(1:pdimlen(3),:)
     ptt(:,:)=ptt(:,:)+pt(ieo,:,:)*damped_scale(1:pdimlen(3),:)
     pqt(:,:)=pqt(:,:)+pq(ieo,:,:)*damped_scale(1:pdimlen(3),:)
     ppt(:,:)=ppt(:,:)+pp(ieo,:,:)*damped_scale(1:pdimlen(3),:)
     pzt(:,:)=pzt(:,:)+pz(ieo,:,:)*damped_scale
! screen and surface perturbatios of variable and time     
     pscrt(:,:)=pscrt(:,:)+ pscr(:,ieo,:)*scales
     psurt(:,:)=psurt(:,:)+ psur(:,ieo,:)*scales
! soil perturbations of variable, time and soil level
     pslt(:,:,:)=pslt(:,:,:)+psl(:,ieo,:,:)*scales

    enddo

   ENDIF
   
! do grid differencing and weighing of profiles

! if perfect model: control_w=1, perturbations are 0 and only the first term is kept.
! if random climo: control_w is between 0 and 1, perturbations are 0, the two first terms are
! kept.
! if eofs perturbations: control_w=1, first and third terms are kept.

    t(:,1:nt) = &
         t3d1(midx,midy,:,wrf_ind:wrf_end_ind)*control_w + &
         t3d2(midx,midy,:,wrf_ind:wrf_end_ind)*(1.0-control_w) + &
         ptt(:,wrf_ind:wrf_end_ind)

    q(:,1:nt) = &
         q3d1(midx,midy,:,wrf_ind:wrf_end_ind)*control_w + &
         q3d2(midx,midy,:,wrf_ind:wrf_end_ind)*(1.0-control_w) + &
         pqt(:,wrf_ind:wrf_end_ind)

    p(:,1:nt) = 1.e2*(&
         p3d1(midx,midy,:,wrf_ind:wrf_end_ind)*control_w + &
         p3d2(midx,midy,:,wrf_ind:wrf_end_ind)*(1.0-control_w)) + &
         ppt(:,wrf_ind:wrf_end_ind)

    u(:,1:nt) = &
         .5*(u3d1(midx,midy,:,wrf_ind:wrf_end_ind)+ &
         u3d1(midx+1,midy,:,wrf_ind:wrf_end_ind))*&
         control_w + &
         .5*(u3d2(midx,midy,:,wrf_ind:wrf_end_ind)+&
         u3d2(midx+1,midy,:,wrf_ind:wrf_end_ind))*&
         (1.0-control_w) + &
         put(:,wrf_ind:wrf_end_ind)

    v(:,1:nt) = &
         .5*(v3d1(midx,midy,:,wrf_ind:wrf_end_ind)+&
         v3d1(midx,midy+1,:,wrf_ind:wrf_end_ind))*&
         control_w + &
         .5*(v3d2(midx,midy,:,wrf_ind:wrf_end_ind)+&
         v3d2(midx,midy+1,:,wrf_ind:wrf_end_ind))*&
         (1.0-control_w) + &
         pvt(:,wrf_ind:wrf_end_ind)

    z_agl_stag(:,1:nt)=MAX(0.,&
         z3d1(midx,midy,:,wrf_ind:wrf_end_ind)*control_w + &
         z3d2(midx,midy,:,wrf_ind:wrf_end_ind)*(1.0-control_w)+&
         pzt(:,wrf_ind:wrf_end_ind) &
         - ter2d(midx,midy,1))


! I will do it from k=2, and for k=1 I will not perturb
! Since pzt is zero anyway, I can start from k=1 anyway
    DO k=2,nz
       z_agl(k,1:nt)=.5*(&
            (z3d1(midx,midy,k,wrf_ind:wrf_end_ind)+&
            z3d1(midx,midy,k+1,wrf_ind:wrf_end_ind))*control_w + &
            (z3d2(midx,midy,k,wrf_ind:wrf_end_ind)+&
            z3d2(midx,midy,k+1,wrf_ind:wrf_end_ind))*&
            (1.0-control_w) + &
            pzt(k,wrf_ind:wrf_end_ind)+ &
            pzt(k+1,wrf_ind:wrf_end_ind)) &
             -ter2d(midx,midy,1)
    ENDDO

        k=1
           z_agl(k,1:nt)=.5*(&
            (z3d1(midx,midy,k,wrf_ind:wrf_end_ind)+&
            z3d1(midx,midy,k+1,wrf_ind:wrf_end_ind))*control_w + &
            (z3d2(midx,midy,k,wrf_ind:wrf_end_ind)+&
            z3d2(midx,midy,k+1,wrf_ind:wrf_end_ind))*&
            (1.0-control_w) + &
            pzt(k,wrf_ind:wrf_end_ind)+ &
            pzt(k+1,wrf_ind:wrf_end_ind)) &
             -ter2d(midx,midy,1)

!*** I haven't done perturbations to the advection part yet****
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
    sts = (/1,1,1,itran1/)
    ierr=nf_get_vara_double(ncid,fid(7),sts,lens,screen2d1(1,:,:,:))
    ierr=nf_get_vara_double(ncid,fid(9),sts,lens,screen2d1(2,:,:,:))
    ierr=nf_get_vara_double(ncid,fid(10),sts,lens,screen2d1(3,:,:,:))
    ierr=nf_get_vara_double(ncid,fid(8),sts,lens,screen2d1(4,:,:,:))

    IF ( rnd_init == 1 ) THEN
       sts = (/1,1,1,itran2/)
    ierr=nf_get_vara_double(ncid,fid(7),sts,lens,screen2d2(1,:,:,:))
    ierr=nf_get_vara_double(ncid,fid(9),sts,lens,screen2d2(2,:,:,:))
    ierr=nf_get_vara_double(ncid,fid(10),sts,lens,screen2d2(3,:,:,:))
    ierr=nf_get_vara_double(ncid,fid(8),sts,lens,screen2d2(4,:,:,:))
    ENDIF

! do grid differencing and weighing of screen variables

    t2(1:nt) = &
         &screen2d1(1,midx,midy,wrf_ind:wrf_end_ind)*control_w + &
         &screen2d2(1,midx,midy,wrf_ind:wrf_end_ind)*(1.0-control_w) + &
         pscrt(1,wrf_ind:wrf_end_ind)
    th2 = t2 ! this is a kluge, will want to use P to correct

    u10(1:nt) = &
         &screen2d1(2,midx,midy,wrf_ind:wrf_end_ind)*control_w + &
         screen2d2(2,midx,midy,wrf_ind:wrf_end_ind)*(1.0-control_w) + &
         pscrt(2,wrf_ind:wrf_end_ind)

    v10(1:nt) = &
         &screen2d1(3,midx,midy,wrf_ind:wrf_end_ind)*control_w + &
         screen2d2(3,midx,midy,wrf_ind:wrf_end_ind)*(1.0-control_w) + &
         pscrt(3,wrf_ind:wrf_end_ind)

    q2(1:nt) = &
         &screen2d1(4,midx,midy,wrf_ind:wrf_end_ind)*control_w + &
         screen2d2(4,midx,midy,wrf_ind:wrf_end_ind)*(1.0-control_w) + &
         pscrt(4,wrf_ind:wrf_end_ind)

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
    stsfc = (/1,1,1,itran1/)
    
    ierr  = nf_get_vara_double(ncid,fid(11),stsfc,lensfc,surf2d1(1,:,:,:))!TSK
    ierr  = nf_get_vara_double(ncid,fid(12),stsfc,lensfc,surf2d1(2,:,:,:))!GLW
    ierr  = nf_get_vara_double(ncid,fid(13),stsfc,lensfc,surf2d1(3,:,:,:))!GSW
    ierr  = nf_get_vara_double(ncid,fid(17),stsfc,lensfc,surf2d1(7,:,:,:))!QSFC
    ierr  = nf_get_vara_double(ncid,fid(18),stsfc,lensfc,surf2d1(8,:,:,:))!VEGFRA
    ierr  = nf_get_vara_double(ncid,fid(19),stsfc,lensfc,surf2d1(9,:,:,:))!ISLTYP
    ierr  = nf_get_vara_double(ncid,fid(20),stsfc,lensfc,surf2d1(10,:,:,:))!IVGTYP
    ierr  = nf_get_vara_double(ncid,fid(21),stsfc,lensfc,surf2d1(11,:,:,:))!LU_INDEX

    IF ( rnd_init == 1 ) THEN
       stsfc = (/1,1,1,itran2/)
    ierr  = nf_get_vara_double(ncid,fid(11),stsfc,lensfc,surf2d2(1,:,:,:))
    ierr  = nf_get_vara_double(ncid,fid(12),stsfc,lensfc,surf2d2(2,:,:,:))
    ierr  = nf_get_vara_double(ncid,fid(13),stsfc,lensfc,surf2d2(3,:,:,:))
    ierr  = nf_get_vara_double(ncid,fid(17),stsfc,lensfc,surf2d2(7,:,:,:))
    ierr  = nf_get_vara_double(ncid,fid(18),stsfc,lensfc,surf2d2(8,:,:,:))
    ierr  = nf_get_vara_double(ncid,fid(19),stsfc,lensfc,surf2d2(9,:,:,:))
    ierr  = nf_get_vara_double(ncid,fid(20),stsfc,lensfc,surf2d2(10,:,:,:))
    ierr  = nf_get_vara_double(ncid,fid(21),stsfc,lensfc,surf2d2(11,:,:,:))

    ENDIF

    ierr  = nf_get_att_double(ncid, fid(17), "_FillValue", missingVal)

! do grid differencing and weighing of surface variables
! required surface variables
    tsk(1:nt) = &
         &surf2d1(1,midx,midy,wrf_ind:wrf_end_ind)*control_w + &
         &surf2d2(1,midx,midy,wrf_ind:wrf_end_ind)*(1.0-control_w) + &
         psurt(1,wrf_ind:wrf_end_ind)

    glw(1:nt) = &
         &surf2d1(2,midx,midy,wrf_ind:wrf_end_ind)*control_w + &
         &surf2d2(2,midx,midy,wrf_ind:wrf_end_ind)*(1.0-control_w) + &
         psurt(2,wrf_ind:wrf_end_ind)

    gsw(1:nt) = &
         &surf2d1(3,midx,midy,wrf_ind:wrf_end_ind)*control_w + &
         &surf2d2(3,midx,midy,wrf_ind:wrf_end_ind)*(1.0-control_w) + &
         psurt(3,wrf_ind:wrf_end_ind)

! need to keep gsw >= 0
    WHERE(gsw < 0.0) gsw = 0.0

! check to see if these are available or else default to *_ref values
    IF ( minval(surf2d1(7,midx,midy,wrf_ind:wrf_end_ind)) > missingVal ) then
       qsfc(1:nt) = &
            surf2d1(7,midx,midy,wrf_ind:wrf_end_ind)*control_w + &
            surf2d2(7,midx,midy,wrf_ind:wrf_end_ind)*(1.0-control_w) + &
            psurt(7,wrf_ind:wrf_end_ind)
    ELSE
       print*,'WARNING: not assigning qsfc'
    ENDIF

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
! do grid differencing and weighing of soil variables
    stsl = (/1,1,1,1,itran1/)
    ierr=nf_get_vara_double(ncid,fid(22),stsl,lensl,soil2d1(1,:,:,:,:))
    ierr=nf_get_vara_double(ncid,fid(23),stsl,lensl,soil2d1(2,:,:,:,:))

    IF ( rnd_init == 1 ) THEN
       stsl = (/1,1,1,1,itran2/)
       ierr=nf_get_vara_double(ncid,fid(22),stsl,lensl,soil2d2(1,:,:,:,:))
       ierr=nf_get_vara_double(ncid,fid(23),stsl,lensl,soil2d2(2,:,:,:,:))
    ENDIF

    tslb(:,1:nt) = &
         &soil2d1(1,midx,midy,:,wrf_ind:wrf_end_ind)*control_w + &
         &soil2d2(1,midx,midy,:,wrf_ind:wrf_end_ind)*(1.0-control_w) + &
         pslt(1,:,wrf_ind:wrf_end_ind)

    smois(:,1:nt) = &
         &soil2d1(2,midx,midy,:,wrf_ind:wrf_end_ind)*control_w + &
         &soil2d2(2,midx,midy,:,wrf_ind:wrf_end_ind)*(1.0-control_w) + &
         pslt(2,:,wrf_ind:wrf_end_ind)

! tslb is already perturbed before so no need to perturb it in this line
    tmn(1:nt) = &
            &tslb(ns,1:pdimlen(6)-wrf_ind+1)

! finally, forecast times 

    DO i=1,nt
       times(i)=(i-1)*interval_f+start_forecast
    ENDDO
    times_flux = times
    times_soil = times
    times_smos = times

    ierr = nf_close(ncid)
    ierr = nf_close(ncid_eofs_init)

    IF ( rnd_init /=  1 ) DEALLOCATE(wrf_year, wrf_month, wrf_day, wrf_hour)
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
    DEALLOCATE(put)
    DEALLOCATE(pvt)
    DEALLOCATE(ptt)
    DEALLOCATE(pqt)
    DEALLOCATE(ppt)
    DEALLOCATE(pzt)
    DEALLOCATE(pscrt)
    DEALLOCATE(psurt)
    DEALLOCATE(pslt)

    IF ( rnd_init == 3 ) THEN
     DEALLOCATE(un)
     DEALLOCATE(vn)
     DEALLOCATE(tn)
     DEALLOCATE(qn)
     DEALLOCATE(pn)
     DEALLOCATE(zn)
     DEALLOCATE(screenn)
     DEALLOCATE(surfn)
     DEALLOCATE(soiln)
     DEALLOCATE(eval)
     DEALLOCATE(pu)
     DEALLOCATE(pv)
     DEALLOCATE(pt)
     DEALLOCATE(pq)
     DEALLOCATE(pp)
     DEALLOCATE(pz)
     DEALLOCATE(pscr)
     DEALLOCATE(psur)
     DEALLOCATE(psl)
     DEALLOCATE(auxfac)
    ENDIF  
 
  END SUBROUTINE wrf_init_and_bc

END MODULE module_wrf_init_and_bc

  
