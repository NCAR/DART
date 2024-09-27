! Copyright 2019 University Corporation for Atmospheric Research and 
! Colorado Department of Public Health and Environment.
!
! Licensed under the Apache License, Version 2.0 (the "License"); you may not use 
! this file except in compliance with the License. You may obtain a copy of the 
! License at      http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed
! under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR 
! CONDITIONS OF ANY KIND, either express or implied. See the License for the 
! specific language governing permissions and limitations under the License.
!
! Development of this code utilized the RMACC Summit supercomputer, which is 
! supported by the National Science Foundation (awards ACI-1532235 and ACI-1532236),
! the University of Colorado Boulder, and Colorado State University.
! The Summit supercomputer is a joint effort of the University of Colorado Boulder
! and Colorado State University.

program iasi_co_profile_ascii_to_obs

! RAWR - observations in the format provided by the distribution center.
! RETR - observations in retrieval space format (ppbv) with the QOR subtractions.
! QOR  - observations in phase space format with QOR transformation.
! CPSR - observations in phase space format with the CPSR trasformations. 
!
!=============================================
! IASI CO retrieval obs
!=============================================
!
  use apm_cpsr_mod, only           : cpsr_calculation, &
                                     mat_prd, &
                                     mat_tri_prd, &
                                     vec_to_mat, &
                                     diag_inv_sqrt, &
                                     lh_mat_vec_prd, &
                                     rh_vec_mat_prd, &
                                     mat_transpose, &
                                     diag_vec
  
  use apm_mapping_mod, only        : w3fb06, &
                                     w3fb07, &
                                     w3fb08, &
                                     w3fb09, &
                                     w3fb11, &
                                     w3fb12, &
                                     w3fb13, &
                                     w3fb14

  use apm_model_fields_vertloc_mod, only : vertical_locate, &
                                           get_model_profile, &
                                           get_DART_diag_data, &
                                           handle_err, &
                                           interp_hori_vert, &
                                           interp_to_obs
  
  use apm_time_code_mod, only          : calc_greg_sec

  use apm_upper_bdy_mod, only      : get_upper_bdy_fld, &
                                     get_MOZART_INT_DATA, &
                                     get_MOZART_REAL_DATA, &
                                     wrf_dart_ubval_interp, &
                                     apm_get_exo_coldens, &
                                     apm_get_upvals, &
                                     apm_interpolate

use    utilities_mod, only : timestamp, 		&
                             register_module, 		&
                             open_file, 		&
                             close_file, 		&
                             initialize_utilities, 	&
                             find_namelist_in_file,  	&
                             check_namelist_read,    	&
                             error_handler, 		&
                             E_ERR,			& 
                             E_WARN,			& 
                             E_MSG, 			&
                             E_DBG

use obs_sequence_mod, only : obs_sequence_type, 	&
                             interactive_obs, 		&
                             write_obs_seq, 		&
                             interactive_obs_sequence,  &
                             static_init_obs_sequence,  &
                             init_obs_sequence,         &
                             init_obs,                  &
                             set_obs_values,            &
                             set_obs_def,               &
                             set_qc,                    &
                             set_qc_meta_data,          &
                             set_copy_meta_data,        &
                             insert_obs_in_seq,         &
                             obs_type

use obs_def_mod, only      : set_obs_def_location,      &
                             set_obs_def_time,          &
                             set_obs_def_key,           &
                             set_obs_def_error_variance,&
                             obs_def_type,              &
                             set_obs_def_type_of_obs

use obs_def_IASI_CO_PROFILE_mod, only : set_obs_def_iasi_co_profile

use  assim_model_mod, only : static_init_assim_model

use location_mod, only : location_type,                 &
                         set_location

use time_manager_mod, only : set_date, 			&
                             set_calendar_type, 	&
                             time_type, 		&
                             get_time

use obs_kind_mod, only   : QTY_CO,                      &
                           IASI_CO_PROFILE,         &
                           get_type_of_obs_from_menu

use random_seq_mod, only : random_seq_type,             &
                           init_random_seq,             &
                           random_uniform

use sort_mod, only       : index_sort

implicit none
!
! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'iasi_co_profile_ascii_to_obs.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''
!
! add variables AFA
type(obs_sequence_type) :: seq
type(obs_type)          :: obs
type(obs_type)          :: obs_old
type(obs_def_type)      :: obs_def
type(location_type)     :: obs_location
type(time_type)         :: obs_time
integer                 :: obs_kind
integer                 :: obs_key
!
integer,parameter       :: fileid=88
integer,parameter       :: max_num_obs=1000000
integer,parameter       :: ias_dim=19, ias_dimp=20
integer,parameter       :: num_copies=1, num_qc=1
integer,parameter       :: lwrk=5*ias_dim
!
real                    :: qc_del,lat_mean
integer                 :: nlon_qc,nlat_qc
real*8                  :: dlon_qc,dlat_qc
type (random_seq_type)  :: inc_ran_seq
!
integer                 :: year, month, day, hour
integer                 :: year1, month1, day1, hour1, minute, second 
integer                 :: year_lst, month_lst, day_lst, hour_lst, minute_lst, second_lst 
integer                 :: iunit, io, icopy, calendar_type
integer                 :: qc_count, ios
integer                 :: nlvls, nlvlsp, index_qc, klvls
integer                 :: lon_qc, lat_qc
integer                 :: i, j, k, l, kk, ik, ikk, k1, k2, kstr 
integer                 :: line_count, index, nlev, nlevp, prs_idx
integer                 :: seconds, days, which_vert, old_ob
integer                 :: kmax,itrm
integer,dimension(max_num_obs)             :: qc_iasi, qc_thinning
integer,dimension(12)                      :: days_in_month=(/ &
                                           31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31  /)
integer,dimension(1000)                    :: index_20
!
real*8                          :: eps_tol=1.e-3,log10e
real*8                          :: dofs, sdof, co_tot_col, co_tot_err
real*8                          :: latitude, longitude, level
real*8                          :: co_psurf, err, co_error, co_prior
real                            :: bin_beg_sec, bin_end_sec
real                            :: sec, lat, lon, nlevels
real                            :: pi ,rad2deg, re, wt, corr_err, fac, fac_obs_error
real                            :: ln_10, xg_sec_avg, co_log_max, co_log_min, co_min
real                            :: prs_loc
real                            :: lat_min,lat_max,lon_min,lon_max,lat_mn
integer                         :: irot, nlvls_fix, nqc_obs
real*8, dimension(1000)         :: unif
real*8, dimension(num_qc)       :: co_qc
real*8, dimension(ias_dim)      :: co_avgker
real*8, dimension(ias_dimp)     :: co_press
real*8, dimension(num_copies)   :: co_vmr
real,dimension(ias_dimp)        :: ias_prs
real,dimension(ias_dim)         :: x_r, x_p, x_p_col, air_col, ret_x_r,ret_x_p,raw_x_r, &
                                   raw_x_p, err2_rs_r, raw_err, ret_err
real,dimension(ias_dim)         :: xcomp, xcomperr, xapr
real,dimension(ias_dim,ias_dim) :: temp_mat, avgker, avg_k, adj_avg_k
real,dimension(ias_dim,ias_dim) :: raw_cov, ret_cov, cov_a, cov_r, cov_m, cov_use
!
double precision,dimension(ias_dim) :: raw_adj_x_r, ret_adj_x_r, adj_x_p 
!
character*129           :: qc_meta_data='IASI CO QC index'
character*129           :: file_name='iasi_co_profile_obs_seq'
character*2             :: chr_month, chr_day, chr_hour
character*4             :: chr_year
character*129           :: filedir, filename, fileout, copy_meta_data, filen
character*129           :: transform_typ
character*129           :: IASI_CO_retrieval_type
character*129           :: IASI_O3_retrieval_type
!
! SUPER OBBING ARRAYS
integer,allocatable,dimension(:,:)       :: xg_count
integer,allocatable,dimension(:,:,:)     :: xg
integer,allocatable,dimension(:,:)       :: xg_nlvls
real,allocatable,dimension(:,:)          :: xg_lon,xg_lat,xg_twt,xg_dof
real,allocatable,dimension(:,:,:)        :: xg_sec,xg_raw_err,xg_ret_err
real,allocatable,dimension(:,:,:)        :: xg_raw_adj_x_r,xg_raw_adj_x_p,xg_raw_x_r,xg_raw_x_p
real,allocatable,dimension(:,:,:)        :: xg_ret_adj_x_r,xg_ret_adj_x_p,xg_ret_x_r,xg_ret_x_p
real,allocatable,dimension(:,:,:)        :: xg_norm,xg_nint
real,allocatable,dimension(:,:,:,:)      :: xg_avg_k,xg_raw_cov,xg_ret_cov
real,allocatable,dimension(:,:,:)        :: xg_prs,xg_prs_norm
!
! QOR/CPSR variables
integer                                        :: info,nlvls_trc,qstatus
integer                                        :: cpsr_co_trunc_lim, cpsr_o3_trunc_lim
double precision,dimension(lwrk)               :: wrk
double precision,allocatable,dimension(:)      :: ZV,SV_cov
double precision,allocatable,dimension(:)      :: rr_x_r,rr_x_p
double precision,allocatable,dimension(:)      :: rs_x_r,rs_x_p
double precision,allocatable,dimension(:,:)    :: Z,ZL,ZR,SV,U_cov,V_cov,UT_cov,VT_cov
double precision,allocatable,dimension(:,:)    :: rr_avg_k,rr_cov
double precision,allocatable,dimension(:,:)    :: rs_avg_k,rs_cov
!
logical                 :: use_log_co
logical                 :: use_log_o3
logical                 :: use_cpsr_co_trunc
logical                 :: use_cpsr_o3_trunc
!
namelist /create_iasi_obs_nml/filedir,filename,fileout,year,month,day,hour, &
         bin_beg_sec,bin_end_sec,IASI_CO_retrieval_type,IASI_O3_retrieval_type, &
         fac_obs_error,use_log_co,use_log_o3,use_cpsr_co_trunc,cpsr_co_trunc_lim, &
         use_cpsr_o3_trunc,cpsr_o3_trunc_lim, &
         lon_min,lon_max,lat_min,lat_max
!
! Set constants
log10e=log10(exp(1.0))
pi=4.*atan(1.)
rad2deg=360./(2.*pi)
re=6371000.
corr_err=.15
corr_err=1.
year_lst=-9999
month_lst=-9999
day_lst=-9999
hour_lst=-9999
minute_lst=-9999
second_lst=-9999 
fac=1.0
lon_min=-9999
lon_max=-9999
lat_min=-9999
lat_max=-9999
!
! Record the current time, date, etc. to the logfile
call initialize_utilities(source)
call register_module(source,revision,revdate)
!
! Initialize the assim_model module, need this to get model
! state meta data for locations of identity observations
!call static_init_assim_model()
!
! Initialize the obs_sequence module
call static_init_obs_sequence()
!
call find_namelist_in_file("input.nml", "create_iasi_obs_nml", iunit)
read(iunit, nml = create_iasi_obs_nml, iostat = io)
call check_namelist_read(iunit, io, "create_iasi_obs_nml")
!
! Record the namelist values used for the run ...
call error_handler(E_MSG,'init_create_iasi_obs','create_iasi_obs_nml values are',' ',' ',' ')
write(     *     , nml=create_iasi_obs_nml)
IASI_CO_retrieval_type='RAWR'
IASI_O3_retrieval_type='RAWR'
!
! Initialize an obs_sequence structure
call init_obs_sequence(seq, num_copies, num_qc, max_num_obs)
!
! Initialize the obs variable
call init_obs(obs, num_copies, num_qc)
!
do icopy =1, num_copies
   if (icopy == 1) then
       copy_meta_data='IASI CO observation'
   else
       copy_meta_data='Truth'
   endif
   call set_copy_meta_data(seq, icopy, copy_meta_data)
enddo

call set_qc_meta_data(seq, 1, qc_meta_data)

qc_iasi(:)=100
qc_thinning(:)=100
!
! assign obs error scale factor
fac=fac_obs_error
!
! define qc arrays (add 1 degree to each side of the qc box)
! qc-del is the size of the qc_bax in km
!
! 39 degress is the central latitude for Colorado:
! nlat_qc=cos(39) * 2. * 3.14 * 6371 km /del_x_km
! nlon_qc=2. *3.14 *6371 km /del_y_km
! dlat_qc=180./nlat_qc
! dlon_qc=360./nlon_qc
!
nqc_obs=40
qc_del=30.
lat_min=lat_min-1.
lat_max=lat_max+1.
lon_min=lon_min-1.
lon_max=lon_max+1.
lat_mean=(lat_min+lat_max)/2.
nlon_qc=(lon_max-lon_min)/360.*2.*pi*re/1000./qc_del+1
nlat_qc=cos(lat_mean/rad2deg)*(lat_max-lat_min)/360.*2.*pi*re/1000./qc_del+1
dlon_qc=(lon_max-lon_min)/(nlon_qc-1)
dlat_qc=(lat_max-lat_min)/(nlat_qc-1)
!
allocate (xg_count(nlon_qc,nlat_qc))
allocate (xg(nlon_qc,nlat_qc,500))
allocate (xg_nlvls(nlon_qc,nlat_qc))
allocate (xg_lon(nlon_qc,nlat_qc),xg_lat(nlon_qc,nlat_qc), &
xg_twt(nlon_qc,nlat_qc),xg_dof(nlon_qc,nlat_qc))
allocate (xg_sec(nlon_qc,nlat_qc,ias_dim),xg_raw_err(nlon_qc,nlat_qc,ias_dim), &
xg_ret_err(nlon_qc,nlat_qc,ias_dim))
allocate (xg_raw_adj_x_r(nlon_qc,nlat_qc,ias_dim),xg_raw_adj_x_p(nlon_qc,nlat_qc,ias_dim), &
xg_raw_x_r(nlon_qc,nlat_qc,ias_dim),xg_raw_x_p(nlon_qc,nlat_qc,ias_dim))
allocate (xg_ret_adj_x_r(nlon_qc,nlat_qc,ias_dim),xg_ret_adj_x_p(nlon_qc,nlat_qc,ias_dim), &
xg_ret_x_r(nlon_qc,nlat_qc,ias_dim), xg_ret_x_p(nlon_qc,nlat_qc,ias_dim))
allocate (xg_norm(nlon_qc,nlat_qc,ias_dim),xg_nint(nlon_qc,nlat_qc,ias_dim))
allocate (xg_avg_k(nlon_qc,nlat_qc,ias_dim,ias_dim), &
xg_raw_cov(nlon_qc,nlat_qc,ias_dim,ias_dim),xg_ret_cov(nlon_qc,nlat_qc,ias_dim,ias_dim))
allocate (xg_prs(nlon_qc,nlat_qc,ias_dimp),xg_prs_norm(nlon_qc,nlat_qc,ias_dimp))
!
!-------------------------------------------------------
! Read IASI obs
!-------------------------------------------------------
!
! Set dates and initialize qc_count
  calendar_type=3                          !Gregorian
  call set_calendar_type(calendar_type)
!
! Perhaps make this a time loop for later runs
  write(chr_year,'(i4.4)') year
  write(chr_month,'(i2.2)') month
  write(chr_day,'(i2.2)') day
  write(chr_hour,'(i2.2)') hour

  if ( mod(year,4) == 0 ) then
       days_in_month(2) = days_in_month(2) + 1
  endif
  if ( mod(year,100) == 0 ) then
       days_in_month(2) = days_in_month(2) - 1
  endif
  if ( mod(year,400) == 0 ) then
       days_in_month(2) = days_in_month(2) + 1
  endif
  if(hour.gt.24) then
     print *, 'APM 1: hour error ',hour
     stop
  endif
!
! Open IASI binary file
  filen=chr_year//chr_month//chr_day//chr_hour//'.dat'
  write(6,*)'opening ',TRIM(filedir)//TRIM(filen)

! Read IASI file 1
  index_qc=0
  line_count = 0
  open(fileid,file=TRIM(filedir)//TRIM(filen),                     &
       form='formatted', status='old',  &
       iostat=ios)
!
! Error Check
  if (ios /=0) then
      write(6,*) 'no iasi file for the day ', day
      go to 999
  endif

! Read IASI
! lon = [-180,180]; lat = [-90,90]
  read(fileid,*,iostat=ios) transform_typ, sec, lat, lon, nlevels, dofs 
  if(lon.lt.0.) lon=lon+360.
!  print *, 'trans_typ, sec, lat, lon, nlevels, dofs ',trim(transform_typ),sec,lat,lon,nlevels,dofs
!
! Error Check
  if (ios /=0) then
      write(6,*) 'no data in file ', TRIM(filen)
      go to 999
  endif
  nlvls=nint(nlevels)
  nlvlsp=nlvls+1
!
!-------------------------------------------------------
! MAIN LOOP FOR IASI OBS
!-------------------------------------------------------
  do while(ios == 0)
       ! Read IASI variables
       read(fileid,*) ias_prs(1:nlvlsp)
       read(fileid,*) x_r(1:nlvls)
       read(fileid,*) x_p(1:nlvls)
       read(fileid,*) x_p_col(1:nlvls)
       read(fileid,*) temp_mat(1:nlvls,1:nlvls)
       do i=1,nlvls
          avg_k(i,1:nlvls)=temp_mat(1:nlvls,i)
       enddo
       read(fileid,*) temp_mat(1:nlvls,1:nlvls)
       do i=1,nlvls
          cov_a(i,1:nlvls)=temp_mat(1:nlvls,i)
       enddo
       read(fileid,*) temp_mat(1:nlvls,1:nlvls)
       do i=1,nlvls
          cov_r(i,1:nlvls)=temp_mat(1:nlvls,i)
       enddo
       read(fileid,*) temp_mat(1:nlvls,1:nlvls)
       do i=1,nlvls
          cov_m(i,1:nlvls)=temp_mat(1:nlvls,i)
       enddo
       read(fileid,*) co_tot_col,co_tot_err
!
       if(lon.ge.lon_min .and. lon.le.lon_max .and. &
       lat.ge.lat_min .and. lat.le.lat_max) then       
          index_qc = index_qc + 1
          qc_iasi(index_qc)=0
!
!-------------------------------------------------------
! Bin to nlat_qc x nlon_qc
!-------------------------------------------------------
! find lon_qc, lat_qc
          lon_qc=nint((lon - lon_min)/dlon_qc)+1
          lat_qc=nint((lat - lat_min)/dlat_qc)+1
!
          xg_count(lon_qc,lat_qc)=xg_count(lon_qc,lat_qc)+1
          xg(lon_qc,lat_qc,xg_count(lon_qc,lat_qc))=index_qc
       endif
!
!read next data point
! lon = [-180,180]
       read(fileid,*,iostat=ios) transform_typ, sec, lat, lon, nlevels, dofs 
       if(lon.lt.0.) lon=lon+360.
!       print *, 'trans_typ, sec, lat, lon, nlevels, dofs ',trim(transform_typ),sec,lat,lon,nlevels,dofs
       nlvls=nint(nlevels)
       nlvlsp=nlvls+1
  enddo !ios

9999   continue

  close(fileid)
!
! Now do the thinning (draw nqc_obs obs)
  call init_random_seq(inc_ran_seq)
  do i=1,nlon_qc
     do j=1,nlat_qc
        if (xg_count(i,j)>nqc_obs) then
           do ik=1,xg_count(i,j)
              unif(ik)=random_uniform(inc_ran_seq)
           enddo
           call index_sort(unif,index_20,xg_count(i,j))
           do ik=1,nqc_obs
              index=xg(i,j,index_20(ik))
              qc_thinning(index)=0
           enddo
        else
           do k=1,xg_count(i,j)
              index=xg(i,j,k)
              qc_thinning(index)=0
           enddo
        endif !xg_count
     enddo !j
  enddo !i
!
!===================================================================================
!
! Read IASI file AGAIN
  index_qc=0
  xg_lon(:,:)=0.
  xg_lat(:,:)=0.
  xg_twt(:,:)=0.
  xg_dof(:,:)=0.
  xg_prs(:,:,:)=0.          
  xg_raw_x_r(:,:,:)=0.          
  xg_raw_x_p(:,:,:)=0.          
  xg_raw_err(:,:,:)=0.          
  xg_raw_adj_x_r(:,:,:)=0.
  xg_raw_adj_x_p(:,:,:)=0.
  xg_raw_cov(:,:,:,:)=0.          
  xg_ret_x_r(:,:,:)=0.          
  xg_ret_x_p(:,:,:)=0.          
  xg_ret_err(:,:,:)=0.          
  xg_ret_adj_x_r(:,:,:)=0.
  xg_ret_adj_x_p(:,:,:)=0.
  xg_ret_cov(:,:,:,:)=0.          
  xg_avg_k(:,:,:,:)=0.          
  xg_norm(:,:,:)=0.
  xg_prs_norm(:,:,:)=0.
  xg_nint(:,:,:)=0.
  xg_sec(:,:,:)=0.
  qc_count=0
!
! NOTE NOTE NOTE Check if it should be BIG_ENDIAN
  open(fileid,file=TRIM(filedir)//TRIM(filen),                     &
       form='formatted', status='old',   &
       iostat=ios)

! Error Check
  if (ios /=0) then
      write(6,*) 'no iasi file for the day ', day
      go to 999
  endif
!
! Read IASI
! lon = [-180,180]
  read(fileid,*,iostat=ios) transform_typ, sec, lat, lon, nlevels, dofs 
  if(lon.lt.0.) lon=lon+360.
!  print *, 'trans_typ, sec, lat, lon, nlevels, dofs ',trim(transform_typ),sec,lat,lon,nlevels,dofs
  nlvls=nint(nlevels)
  nlvlsp=nlvls+1
!
! Error Check
  if (ios /=0) then
      write(6,*) 'no data on file ', TRIM(filen)
      go to 999
  endif
!
!-------------------------------------------------------
! MAIN LOOP FOR IASI OBS (ppbv)
!-------------------------------------------------------
!
  do while(ios == 0)
     index_qc=index_qc+1
     read(fileid,*) ias_prs(1:nlvlsp)
     read(fileid,*) x_r(1:nlvls)
     read(fileid,*) x_p(1:nlvls)
     read(fileid,*) x_p_col(1:nlvls)
       read(fileid,*) temp_mat(1:nlvls,1:nlvls)
       do i=1,nlvls
          avg_k(i,1:nlvls)=temp_mat(1:nlvls,i)
       enddo
       read(fileid,*) temp_mat(1:nlvls,1:nlvls)
       do i=1,nlvls
          cov_a(i,1:nlvls)=temp_mat(1:nlvls,i)
       enddo
       read(fileid,*) temp_mat(1:nlvls,1:nlvls)
       do i=1,nlvls
          cov_r(i,1:nlvls)=temp_mat(1:nlvls,i)
       enddo
       read(fileid,*) temp_mat(1:nlvls,1:nlvls)
       do i=1,nlvls
          cov_m(i,1:nlvls)=temp_mat(1:nlvls,i)
       enddo
       read(fileid,*) co_tot_col,co_tot_err
!        print *, 'nlvls ',nlvls
!        print *, 'x_r ',x_r(1:nlvls)
!        print *, 'x_p ',x_p(1:nlvls)
!        do k=1,nlvls
!           print *, 'k, avg_k ',k,avg_k(k,1:nlvls)
!           print *, 'k, cov_r ',k,cov_r(k,1:nlvls)
!        enddo
!     endif
!
! calculate the air column for column to vmr conversion
     do i=1,nlvls
        air_col(i)=x_p_col(i)/x_p(i)
     enddo
!
! convert averaging kernel from column to vmr units
     do i=1,nlvls
        do j=1,nlvls
           avg_k(i,j)=avg_k(i,j)/air_col(i)*air_col(j)
        enddo
     enddo       
!
     if ( (qc_iasi(index_qc)==0).and.(qc_thinning(index_qc)==0) ) then
        co_qc(1)=0
     else
        co_qc(1)=100
     endif
!
     if ( co_qc(1) == 0 )  then
!
! calculate bin indexes
        lon_qc=nint((lon - lon_min)/dlon_qc)+1
        lat_qc=nint((lat - lat_min)/dlat_qc)+1
!
! Assign cov_use
        cov_use(:,:)=cov_r(:,:)
!
! Calculate averaging kernel complement: (I-A) (RAWR form)
        raw_adj_x_r(:)=0.
        ret_adj_x_r(:)=0.
        adj_x_p(:)=0.
        adj_avg_k(:,:)=0.
        do i=1,nlvls
           do j=1,nlvls
              adj_avg_k(i,j)=-1.*avg_k(i,j)
           enddo
           adj_avg_k(i,i)=adj_avg_k(i,i)+1.
        enddo
!
! Calcuate the prior term: (I-A) x_p (RAWR form)
        call lh_mat_vec_prd(dble(adj_avg_k(1:nlvls,1:nlvls)),dble(x_p(1:nlvls)),adj_x_p(1:nlvls),nlvls)
!
! Calculate the QOR term: x_r - (I-A) x_p (RAWR form)
        raw_adj_x_r(1:nlvls)=x_r(1:nlvls)-adj_x_p(1:nlvls)
!
! Calculate the QOR term: x_r - (I-A) x_p (RETR form)
        ret_adj_x_r(1:nlvls)=x_r(1:nlvls)-adj_x_p(1:nlvls)
!
! Calculate RAWR retrieval and prior 
        do i=1,nlvls
           raw_x_r(i)=x_r(i)
           raw_x_p(i)=x_p(i)
        enddo
!
! Calculate RETR retrieval and prior 
        do i=1,nlvls
           ret_x_r(i)=x_r(i)
           ret_x_p(i)=x_p(i)
        enddo
!
! Calculate RAWR errors
        do i=1,nlvls
           do j=1,nlvls
              raw_cov(i,j)=cov_use(i,j)
           enddo
        enddo
!
! Calculate RETR errors
        do i=1,nlvls
           do j=1,nlvls
!              ret_cov(i,j)=cov_use(i,j)/raw_x_r(i)/raw_x_r(j)
              ret_cov(i,j)=cov_use(i,j)
           enddo
        enddo
!
! Calculate errors for RAWR case
        do j=1,nlvls
           raw_err(j)=sqrt(raw_cov(j,j))
        enddo
!
! Calculate errors for RETR case
        do j=1,nlvls
           ret_err(j)=sqrt(ret_cov(j,j))
        enddo
!
! Calculate superobs 
        kstr=ias_dim-nlvls+1
        wt=cos(lat/rad2deg)
        xg_twt(lon_qc,lat_qc)=xg_twt(lon_qc,lat_qc)+wt
        xg_lon(lon_qc,lat_qc)=xg_lon(lon_qc,lat_qc)+lon*wt
        xg_lat(lon_qc,lat_qc)=xg_lat(lon_qc,lat_qc)+lat*wt
        xg_dof(lon_qc,lat_qc)=xg_dof(lon_qc,lat_qc)+dofs*wt
        do i=kstr,ias_dim
           xg_norm(lon_qc,lat_qc,i)=xg_norm(lon_qc,lat_qc,i)+wt
           xg_prs_norm(lon_qc,lat_qc,i)=xg_prs_norm(lon_qc,lat_qc,i)+wt
           xg_nint(lon_qc,lat_qc,i)=xg_nint(lon_qc,lat_qc,i)+1
           if(hour.eq.24 .and. sec.le.10800) then
              xg_sec(lon_qc,lat_qc,i)=xg_sec(lon_qc,lat_qc,i)+(86400+sec)*wt
           else
              xg_sec(lon_qc,lat_qc,i)=xg_sec(lon_qc,lat_qc,i)+sec*wt
           endif
           xg_prs(lon_qc,lat_qc,i)=xg_prs(lon_qc,lat_qc,i)+ias_prs(i-kstr+1)*wt
           xg_raw_x_r(lon_qc,lat_qc,i)=xg_raw_x_r(lon_qc,lat_qc,i)+raw_x_r(i-kstr+1)*wt
           xg_raw_x_p(lon_qc,lat_qc,i)=xg_raw_x_p(lon_qc,lat_qc,i)+raw_x_p(i-kstr+1)*wt
           xg_raw_err(lon_qc,lat_qc,i)=xg_raw_err(lon_qc,lat_qc,i)+raw_err(i-kstr+1)*wt
           xg_raw_adj_x_r(lon_qc,lat_qc,i)=xg_raw_adj_x_r(lon_qc,lat_qc,i)+raw_adj_x_r(i-kstr+1)*wt
!           xg_raw_adj_x_p(lon_qc,lat_qc,i)=xg_raw_adj_x_p(lon_qc,lat_qc,i)+adj_x_p(i-kstr+1)*wt
           xg_ret_x_r(lon_qc,lat_qc,i)=xg_ret_x_r(lon_qc,lat_qc,i)+ret_x_r(i-kstr+1)*wt
           xg_ret_x_p(lon_qc,lat_qc,i)=xg_ret_x_p(lon_qc,lat_qc,i)+ret_x_p(i-kstr+1)*wt
           xg_ret_err(lon_qc,lat_qc,i)=xg_ret_err(lon_qc,lat_qc,i)+ret_err(i-kstr+1)*wt
           xg_ret_adj_x_r(lon_qc,lat_qc,i)=xg_ret_adj_x_r(lon_qc,lat_qc,i)+ret_adj_x_r(i-kstr+1)*wt
!           xg_ret_adj_x_p(lon_qc,lat_qc,i)=xg_ret_adj_x_p(lon_qc,lat_qc,i)+adj_x_p(i-kstr+1)*wt
           do j=kstr,ias_dim
              xg_ret_cov(lon_qc,lat_qc,i,j)=xg_ret_cov(lon_qc,lat_qc,i,j)+ret_cov(i-kstr+1,j-kstr+1)*wt
              xg_raw_cov(lon_qc,lat_qc,i,j)=xg_raw_cov(lon_qc,lat_qc,i,j)+raw_cov(i-kstr+1,j-kstr+1)*wt
              xg_avg_k(lon_qc,lat_qc,i,j)=xg_avg_k(lon_qc,lat_qc,i,j)+avg_k(i-kstr+1,j-kstr+1)*wt
           enddo
        enddo
        xg_prs_norm(lon_qc,lat_qc,ias_dimp)=xg_prs_norm(lon_qc,lat_qc,ias_dimp)+wt
        xg_prs(lon_qc,lat_qc,ias_dimp)=xg_prs(lon_qc,lat_qc,ias_dimp)+ias_prs(nlvlsp)*wt
     endif    ! co_qc(1)
!
! read next data point
! lon = [-180,180]
     read(fileid,*,iostat=ios) transform_typ, sec, lat, lon, nlevels, dofs 
     if(lon.lt.0) lon=lon+360.
     nlvls=nint(nlevels)
     nlvlsp=nlvls+1
  enddo    !ios
!
! Calculate number of vertical levels and averages
  qc_count=0
  do i=1,nlon_qc  
     do j=1,nlat_qc
        if(xg_twt(i,j).eq.0) cycle
        xg_lon(i,j)=xg_lon(i,j)/xg_twt(i,j)
        xg_lat(i,j)=xg_lat(i,j)/xg_twt(i,j)
        xg_dof(i,j)=xg_dof(i,j)/xg_twt(i,j)
        do k=1,ias_dim
           if(xg_norm(i,j,k).eq.0) cycle
           xg_sec(i,j,k)=xg_sec(i,j,k)/xg_norm(i,j,k)
           if(xg_sec(i,j,k).ge.86400) xg_sec(i,j,k)=xg_sec(i,j,k)-86400
           xg_prs(i,j,k)=xg_prs(i,j,k)/real(xg_norm(i,j,k))
           xg_raw_x_r(i,j,k)=xg_raw_x_r(i,j,k)/real(xg_norm(i,j,k))
           xg_raw_x_p(i,j,k)=xg_raw_x_p(i,j,k)/real(xg_norm(i,j,k))
           xg_raw_err(i,j,k)=xg_raw_err(i,j,k)/sqrt(real(xg_norm(i,j,k)))
           xg_raw_adj_x_r(i,j,k)=xg_raw_adj_x_r(i,j,k)/real(xg_norm(i,j,k))
!           xg_raw_adj_x_p(i,j,k)=xg_raw_adj_x_p(i,j,k)/real(xg_norm(i,j,k))
           xg_ret_x_r(i,j,k)=xg_ret_x_r(i,j,k)/real(xg_norm(i,j,k))
           xg_ret_x_p(i,j,k)=xg_ret_x_p(i,j,k)/real(xg_norm(i,j,k))
           xg_ret_err(i,j,k)=xg_ret_err(i,j,k)/sqrt(real(xg_norm(i,j,k)))
           xg_ret_adj_x_r(i,j,k)=xg_ret_adj_x_r(i,j,k)/real(xg_norm(i,j,k))
!           xg_ret_adj_x_p(i,j,k)=xg_ret_adj_x_p(i,j,k)/real(xg_norm(i,j,k))
           do l=1,ias_dim
              if(xg_norm(i,j,l).eq.0) cycle
              xg_raw_cov(i,j,k,l)=xg_raw_cov(i,j,k,l)/real(xg_norm(i,j,k))
              xg_ret_cov(i,j,k,l)=xg_ret_cov(i,j,k,l)/real(xg_norm(i,j,k))
              xg_avg_k(i,j,k,l)=xg_avg_k(i,j,k,l)/real(xg_norm(i,j,k))
           enddo
        enddo
        xg_prs(i,j,ias_dimp)=xg_prs(i,j,ias_dimp)/real(xg_prs_norm(i,j,ias_dimp))
!
! Get number of vertical levels
        klvls=ias_dim
        do k=1,ias_dim
          if(xg_norm(i,j,k).eq.0) then
             klvls=klvls-1
          endif
        enddo
        xg_nlvls(i,j)=klvls
        nlvls=xg_nlvls(i,j)
        kstr=ias_dim-xg_nlvls(i,j)+1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! QOR CODE BLOCK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        if(trim(IASI_CO_retrieval_type) .eq. 'QOR') then
!
! Calculate SVD of raw_cov (Z=U_xxx * SV_xxx * VT_xxx)
           allocate(Z(nlvls,nlvls),SV_cov(nlvls),SV(nlvls,nlvls))
           allocate(U_cov(nlvls,nlvls),UT_cov(nlvls,nlvls),V_cov(nlvls,nlvls),VT_cov(nlvls,nlvls))
           allocate(rs_avg_k(nlvls,nlvls),rs_cov(nlvls,nlvls),rs_x_r(nlvls),rs_x_p(nlvls))       
           allocate(rr_avg_k(nlvls,nlvls),rr_cov(nlvls,nlvls),rr_x_r(nlvls),rr_x_p(nlvls))       
           allocate(ZL(nlvls,nlvls),ZR(nlvls,nlvls),ZV(nlvls))
!
           Z(1:nlvls,1:nlvls)=dble(xg_raw_cov(i,j,kstr:ias_dim,kstr:ias_dim))
           call dgesvd('A','A',nlvls,nlvls,Z,nlvls,SV_cov,U_cov,nlvls,VT_cov,nlvls,wrk,lwrk,info)
           nlvls_trc=0
           do k=1,nlvls
              if(SV_cov(k).ge.eps_tol) then
                 nlvls_trc=k
              else
                 SV_cov(k)=0
                 U_cov(:,k)=0. 
                 VT_cov(k,:)=0.
              endif 
           enddo
!              print *,'nlvls_trc ',nlvls_trc
!              print *, 'SV ',SV_cov(:)
!
! Scale the singular vectors
           do k=1,nlvls_trc
              U_cov(:,k)=U_cov(:,k)/sqrt(SV_cov(k))
           enddo
!              print *, 'nlvls_trc ',nlvls_trc
!              print *, 'SV ',SV_cov(:)
!        
           call mat_transpose(U_cov,UT_cov,nlvls,nlvls)
           call mat_transpose(VT_cov,V_cov,nlvls,nlvls)
           call vec_to_mat(SV_cov,SV,nlvls)
!              do k=1,nlvls
!                print *, 'U ',k,(U_cov(k,l),l=1,nlvls)
!              enddo
!              do k=1,nlvls
!                print *, 'UT ',k,(UT_cov(k,l),l=1,nlvls)
!              enddo
!
! Rotate terms in the forward operator
           ZL(1:nlvls,1:nlvls)=dble(xg_avg_k(i,j,kstr:ias_dim,kstr:ias_dim))
           call mat_prd(UT_cov(1:nlvls,1:nlvls),ZL(1:nlvls,1:nlvls), &
           rs_avg_k(1:nlvls,1:nlvls),nlvls,nlvls,nlvls,nlvls)
!              do k=1,nlvls
!                print *, 'UT ',k,(UT_cov(k,l),l=1,nlvls)
!              enddo
!              do k=kstr,ias_dim
!                print *, 'xg_avg_k ',k,(xg_avg_k(i,j,k,l),l=kstr,ias_dim)
!              enddo
!              do k=1,nlvls
!                print *, 'rs_avg_k ',k,(rs_avg_k(k,l),l=1,nlvls)
!              enddo
           ZL(1:nlvls,1:nlvls)=dble(xg_raw_cov(i,j,kstr:ias_dim,kstr:ias_dim))
           call mat_tri_prd(UT_cov(1:nlvls,1:nlvls),ZL(1:nlvls,1:nlvls),U_cov(1:nlvls,1:nlvls), &
           rs_cov(1:nlvls,1:nlvls),nlvls,nlvls,nlvls,nlvls,nlvls,nlvls)
!              do k=1,nlvls
!                print *, 'UT ',k,(UT_cov(k,l),l=1,nlvls)
!              enddo
!              do k=kstr,ias_dim
!                print *, 'xg_raw_cov ',k,(xg_raw_cov(i,j,k,l),l=kstr,ias_dim)
!              enddo
!              do k=1,nlvls
!                print *, 'rs_cov ',k,(rs_cov(k,l),l=1,nlvls)
!              enddo
           ZV(1:nlvls)=dble(xg_raw_adj_x_r(i,j,kstr:ias_dim))
           call lh_mat_vec_prd(UT_cov(1:nlvls,1:nlvls),ZV(1:nlvls),rs_x_r(1:nlvls),nlvls)
!              do k=1,nlvls
!                print *, 'UT ',k,(UT_cov(k,l),l=1,nlvls)
!              enddo
!              print *, 'xg_adj_x_r ',(xg_raw_adj_x_r(i,j,l),l=kstr,ias_dim)
!              print *, 'rs_x_r ',(rs_x_r(l),l=1,nlvls)
           ZV(1:nlvls)=dble(xg_raw_adj_x_p(i,j,kstr:ias_dim))
           call lh_mat_vec_prd(UT_cov(1:nlvls,1:nlvls),ZV(1:nlvls),rs_x_p(1:nlvls),nlvls)
!              do k=1,nlvls
!                 print *, 'UT ',k,(UT_cov(k,l),l=1,nlvls)
!              enddo
!              print *, 'xg_adj_x_p ',(xg_raw_adj_x_p(i,j,l),l=kstr,ias_dim)
!              print *, 'rs_x_p ',(rs_x_p(l),l=1,nlvls)
!
! Get new errors (check if err2_rs_r < 0 the qstatus=1)
           qstatus=0.0
           do k=1,nlvls
              err2_rs_r(k)=sqrt(rs_cov(k,k))
           enddo
        endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! CPSR CODE BLOCK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate SVD of avg_k (Z=U_xxx * SV_xxx * VT_xxx) - FIRST ROTATION
        if(trim(IASI_CO_retrieval_type) .eq. 'CPSR') then
           allocate(Z(nlvls,nlvls),SV_cov(nlvls),SV(nlvls,nlvls))
           allocate(U_cov(nlvls,nlvls),UT_cov(nlvls,nlvls),V_cov(nlvls,nlvls),VT_cov(nlvls,nlvls))
           allocate(rs_avg_k(nlvls,nlvls),rs_cov(nlvls,nlvls),rs_x_r(nlvls),rs_x_p(nlvls))       
           allocate(rr_avg_k(nlvls,nlvls),rr_cov(nlvls,nlvls),rr_x_r(nlvls),rr_x_p(nlvls))       
           allocate(ZL(nlvls,nlvls),ZR(nlvls,nlvls),ZV(nlvls))
           Z(1:nlvls,1:nlvls)=dble(xg_avg_k(i,j,kstr:ias_dim,kstr:ias_dim))
           call dgesvd('A','A',nlvls,nlvls,Z,nlvls,SV_cov,U_cov,nlvls,VT_cov,nlvls,wrk,lwrk,info)
           nlvls_trc=0
           sdof=0.
!
! APM: the phase space truncation should not be done here
! because it impacts the projection of the compressed averaging kernel
! onto the left on-zero singular vectors of the error covariance.
           do k=1,nlvls
              if(SV_cov(k).ge.eps_tol) then
                 nlvls_trc=k
                 sdof=sdof+SV_cov(k)
              else
                 SV_cov(k)=0
                 U_cov(:,k)=0. 
                 VT_cov(k,:)=0.
              endif 
           enddo
!              print *,'nlvls_trc ',nlvls_trc
!              print *, 'SV ',SV_cov(:)
           call mat_transpose(U_cov,UT_cov,nlvls,nlvls)
           call mat_transpose(VT_cov,V_cov,nlvls,nlvls)
           call vec_to_mat(SV_cov,SV,nlvls)
!
! Rotate terms in the forward operator
           ZL(1:nlvls,1:nlvls)=dble(xg_avg_k(i,j,kstr:ias_dim,kstr:ias_dim))
! averaging kernel
           call mat_prd(UT_cov(1:nlvls,1:nlvls),ZL(1:nlvls,1:nlvls), &
           rr_avg_k(1:nlvls,1:nlvls),nlvls,nlvls,nlvls,nlvls)
!              do k=1,nlvls
!                print *, 'UT ',k,(UT_cov(k,l),l=1,nlvls)
!              enddo
!              do k=kstr,ias_dim
!                print *, 'xg_avg_k ',k,(xg_avg_k(i,j,k,l),l=kstr,ias_dim)
!              enddo
!              do k=1,nlvls
!                print *, 'rr_avg_k ',k,(rr_avg_k(k,l),l=1,nlvls)
!              enddo

! retrieval error covariance
           ZL(1:nlvls,1:nlvls)=dble(xg_raw_cov(i,j,kstr:ias_dim,kstr:ias_dim))
           call mat_tri_prd(UT_cov(1:nlvls,1:nlvls),ZL(1:nlvls,1:nlvls),U_cov(1:nlvls,1:nlvls), &
           rr_cov(1:nlvls,1:nlvls),nlvls,nlvls,nlvls,nlvls,nlvls,nlvls)
!              do k=1,nlvls
!                print *, 'UT ',k,(UT_cov(k,l),l=1,nlvls)
!              enddo
!              do k=kstr,ias_dim
!                print *, 'xg_raw_cov ',k,(xg_raw_cov(i,j,k,l),l=kstr,ias_dim)
!              enddo
!              do k=1,nlvls
!                print *, 'rr_cov ',k,(rr_cov(k,l),l=1,nlvls)
!              enddo
! adjusted retrieval
           ZV(1:nlvls)=dble(xg_raw_adj_x_r(i,j,kstr:ias_dim))
           call lh_mat_vec_prd(UT_cov(1:nlvls,1:nlvls),ZV(1:nlvls),rr_x_r(1:nlvls),nlvls)
!              do k=1,nlvls
!                print *, 'UT ',k,(UT_cov(k,l),l=1,nlvls)
!              enddo
!              print *, 'xg_adj_x_r ',(xg_raw_adj_x_r(i,j,l),l=kstr,ias_dim)
!              print *, 'rr_x_r ',(rr_x_r(l),l=1,nlvls)
           ZV(1:nlvls)=dble(xg_raw_adj_x_p(i,j,kstr:ias_dim))
           call lh_mat_vec_prd(UT_cov(1:nlvls,1:nlvls),ZV(1:nlvls),rr_x_p(1:nlvls),nlvls)
!              do k=1,nlvls
!                 print *, 'UT ',k,(UT_cov(k,l),l=1,nlvls)
!              enddo
!              print *, 'xg_adj_x_p ',(xg_raw_adj_x_p(i,j,l),l=kstr,ias_dim)
!              print *, 'rr_x_p ',(rr_x_p(l),l=1,nlvls)
!
! Calculate SVD of rr_cov (Z=U_xxx * SV_xxx * VT_xxx) - SECOND ROTATION
           Z(1:nlvls,1:nlvls)=rr_cov(1:nlvls,1:nlvls)
           call dgesvd('A','A',nlvls,nlvls,Z,nlvls,SV_cov,U_cov,nlvls,VT_cov,nlvls,wrk,lwrk,info)
           do k=nlvls_trc+1,nlvls
              SV_cov(k)=0
              U_cov(:,k)=0. 
              VT_cov(k,:)=0.
           enddo
!
! Scale the singular vectors
           do k=1,nlvls_trc
              U_cov(:,k)=U_cov(:,k)/sqrt(SV_cov(k))
           enddo
!              print *, 'nlvls_trc ',nlvls_trc
!              print *, 'SV ',SV_cov(:)
!          
           call mat_transpose(U_cov,UT_cov,nlvls,nlvls)
           call mat_transpose(VT_cov,V_cov,nlvls,nlvls)
           call vec_to_mat(SV_cov,SV,nlvls)
!              do k=1,nlvls
!                print *, 'U ',k,(U_cov(k,l),l=1,nlvls)
!              enddo
!              do k=1,nlvls
!                print *, 'UT ',k,(UT_cov(k,l),l=1,nlvls)
!              enddo
!
! Rotate terms in the forward operator
           ZL(:,:)=0.
           do k=1,nlvls
              do kk=1,nlvls
                 ZL(k,kk)=rr_avg_k(k,kk)
              enddo
           enddo
           call mat_prd(UT_cov(1:nlvls,1:nlvls),ZL(1:nlvls,1:nlvls), &
           rs_avg_k(1:nlvls,1:nlvls),nlvls,nlvls,nlvls,nlvls)
!              do k=1,nlvls
!                print *, 'UT ',k,(UT_cov(k,l),l=1,nlvls)
!              enddo
!              do k=1,nlvls
!                print *, 'rr_avg_k ',k,(rr_avg_k(k,l),l=1,nlvls)
!              enddo
!              do k=1,nlvls
!                print *, 'rs_avg_k ',k,(rs_avg_k(k,l),l=1,nlvls)
!              enddo
           ZL(:,:)=0.
           do k=1,nlvls
              do kk=1,nlvls
                 ZL(k,kk)=rr_cov(k,kk)
              enddo
           enddo
           call mat_tri_prd(UT_cov(1:nlvls,1:nlvls),ZL(1:nlvls,1:nlvls),U_cov(1:nlvls,1:nlvls), &
           rs_cov(1:nlvls,1:nlvls),nlvls,nlvls,nlvls,nlvls,nlvls,nlvls)
!              do k=1,nlvls
!                print *, 'UT ',k,(UT_cov(k,l),l=1,nlvls)
!              enddo
!              do k=1,nlvls
!                print *, 'rr_cov ',k,(rr_cov(k,l),l=1,nlvls)
!              enddo
!              do k=1,nlvls
!                print *, 'rs_cov ',k,(rs_cov(k,l),l=1,nlvls)
!              enddo
           ZV(:)=0.
           do k=1,nlvls    
              ZV(k)=rr_x_r(k)
           enddo
           call lh_mat_vec_prd(UT_cov(1:nlvls,1:nlvls),ZV(1:nlvls),rs_x_r(1:nlvls),nlvls)
!              do k=1,nlvls
!                print *, 'UT ',k,(UT_cov(k,l),l=1,nlvls)
!              enddo
!              print *, 'rr_x_r ',(rr_x_r(l),l=1,nlvls)
!              print *, 'rs_x_r ',(rs_x_r(l),l=1,nlvls)
           ZV(:)=0.
           do k=1,nlvls
              ZV(k)=rr_x_p(k)
           enddo
           call lh_mat_vec_prd(UT_cov(1:nlvls,1:nlvls),ZV(1:nlvls),rs_x_p(1:nlvls),nlvls)
!              do k=1,nlvls
!                 print *, 'UT ',k,(UT_cov(k,l),l=1,nlvls)
!              enddo
!              print *, 'rr_x_p ',(rr_x_p(l),l=1,nlvls)
!              print *, 'rs_x_p ',(rs_x_p(l),l=1,nlvls)
!
! Get new errors (check if err2_rs_r < 0 the qstatus=1)
           qstatus=0.0
           do k=1,nlvls
              err2_rs_r(k)=sqrt(rs_cov(k,k))
           enddo
        endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! APM make assignments to Ave's scaled variables
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Set vertical for levels (irot=0) or modes (irot=1)
        if(trim(IASI_CO_retrieval_type).eq.'RAWR' .or. &
        trim(IASI_CO_retrieval_type).eq.'RETR') then
           irot=0
           nlvls_fix=xg_nlvls(i,j)
        else
           irot=1
           nlvls_fix=nlvls_trc
        endif
!
! Truncate the number of CPSR modes
        if(use_cpsr_co_trunc .and. nlvls_fix.gt.cpsr_co_trunc_lim) then
           print *, 'APM: change limit ',nlvls_fix, cpsr_co_trunc_lim
           nlvls_fix=cpsr_co_trunc_lim
        endif
!        
        do k=1,nlvls_fix
           qc_count=qc_count+1
!
! RAWR
           if(trim(IASI_CO_retrieval_type) .eq. 'RAWR') then
              xcomp(k)=xg_raw_adj_x_r(i,j,k+kstr-1)
              xcomperr(k)=fac*xg_raw_err(i,j,k+kstr-1)
!              xapr(k)=xg_raw_adj_x_p(i,j,k+kstr-1)
              xapr(k)=0.
              do l=1,xg_nlvls(i,j)
                 avgker(k,l)=xg_avg_k(i,j,k+kstr-1,l+kstr-1)
              enddo
           endif
!
! RETR
           if(trim(IASI_CO_retrieval_type) .eq. 'RETR') then
              xcomp(k)=xg_ret_adj_x_r(i,j,k+kstr-1)
              xcomperr(k)=fac*xg_ret_err(i,j,k+kstr-1)
!              xapr(k)=xg_ret_adj_x_p(i,j,k+kstr-1)
              xapr(k)=0.
              do l=1,xg_nlvls(i,j)
                 avgker(k,l)=xg_avg_k(i,j,k+kstr-1,l+kstr-1)
              enddo
           endif
!
! RAWR QOR
           if(trim(IASI_CO_retrieval_type) .eq. 'QOR') then
              xcomp(k)=rs_x_r(k)
              xcomperr(k)=fac*err2_rs_r(k)
              xapr(k)=0.
              do l=1,xg_nlvls(i,j)
                 avgker(k,l)=rs_avg_k(k,l)
              enddo
           endif
!
! RAWR CPSR
           if(trim(IASI_CO_retrieval_type) .eq. 'CPSR') then
              xcomp(k)=rs_x_r(k)
              xcomperr(k)=fac*err2_rs_r(k)
              xapr(k)=0.
              do l=1,xg_nlvls(i,j)
                 avgker(k,l)=rs_avg_k(k,l)
              enddo
           endif
!
! Calculate average seconds
           xg_sec_avg=0.
           do l=1,xg_nlvls(i,j)
              xg_sec_avg=xg_sec_avg+xg_sec(i,j,l+kstr-1)/xg_nlvls(i,j)
           enddo
           if(irot.eq.0) xg_sec_avg=xg_sec(i,j,k+kstr-1)
!
!--------------------------------------------------------
! assign obs variables for obs_sequence
!--------------------------------------------------------
!
! location
           latitude=xg_lat(i,j) 
           longitude=xg_lon(i,j)
!
! time (get time from sec IASI variable)
           hour1 = int(xg_sec_avg/3600d0)
           if(hour1.gt.24) then
              print *, 'APM 2: hour error ',hour1,xg_sec_avg
              stop
           endif
           minute = int( (xg_sec_avg-hour1*3600d0)/60d0)
           second = int(xg_sec_avg - hour1*3600d0 - minute*60d0)
           if(hour1.gt.24.or.hour.gt.24) then
              print *, 'APM 3: hour error ',hour1,hour
              stop
           endif
           if ( hour == 24 ) then
              if (xg_sec_avg < 3.00*3600d0) then
                 day1 = day+1
                 if (day1 > days_in_month(month)) then
                    day1 = 1
                    if (month < 12) then
                       month1 = month + 1
                       year1 = year
                    else
                       month1 = 1
                       year1  = year+1
                    endif
                 else
                    month1 = month
                    year1 = year
                 endif
              else
                 day1 = day
                 month1 = month
                 year1 = year
              endif
           else
              day1 = day
              month1 = month
              year1 = year
           endif
           old_ob=0
           if(year1.lt.year_lst) then
              old_ob=1
           else if(year1.eq.year_lst .and. month1.lt.month_lst) then
              old_ob=1
           else if(year1.eq.year_lst .and. month1.eq.month_lst .and. &
           day1.lt.day_lst) then
              old_ob=1
           else if(year1.eq.year_lst .and. month1.eq.month_lst .and. &
           day1.eq.day_lst .and. hour1.lt.hour_lst) then
              old_ob=1
           else if(year1.eq.year_lst .and. month1.eq.month_lst .and. &
           day1.eq.day_lst .and. hour1.eq.hour_lst .and. minute.lt.minute_lst) then
              old_ob=1
           else if(year1.eq.year_lst .and. month1.eq.month_lst .and. &
           day1.eq.day_lst .and. hour1.eq.hour_lst .and. &
           minute.eq.minute_lst .and. second.lt.second_lst) then
              old_ob=1
           else if(year1.gt.year_lst .or. month1.gt.month_lst .or. &
           day1.gt.day_lst .or. hour1.gt.hour_lst .or. &
           minute.gt.minute_lst .or. second.gt.second_lst) then
              year_lst=year1
              month_lst=month1
              day_lst=day1
              hour_lst=hour1
              minute_lst=minute
              second_lst=second
           endif
           obs_time=set_date(year1,month1,day1,hour1,minute,second)
           call get_time(obs_time, seconds, days)
!
!--------------------------------------------------------
! Loop through the ias_dim levels for now
! Use each mixing ratio as a separate obs
!--------------------------------------------------------
!
           if(irot.eq.1) then
              level=1
              which_vert=-2         ! undefined
           elseif(irot.eq.0) then  
              level=(xg_prs(i,j,k+kstr-1)+xg_prs(i,j,k+kstr))/2*100.
              which_vert=2          ! pressure surface
           endif
           obs_kind = IASI_CO_PROFILE
           obs_location=set_location(longitude, latitude, level, which_vert)
           co_psurf=xg_prs(i,j,kstr)*100.
           co_avgker(1:xg_nlvls(i,j))=avgker(k,1:xg_nlvls(i,j))
           co_press(1:xg_nlvls(i,j)+1)=xg_prs(i,j,kstr:ias_dimp)*100.
           co_prior=xapr(k)
           co_vmr(1)=xcomp(k)
           err = xcomperr(k)
           co_error=err*err
!
           call set_obs_def_type_of_obs(obs_def, obs_kind)
           call set_obs_def_location(obs_def, obs_location)
           call set_obs_def_time(obs_def, obs_time)
           call set_obs_def_error_variance(obs_def, co_error)
           call set_obs_def_iasi_co_profile(qc_count, co_avgker, co_press, co_prior, co_psurf, &
           xg_nlvls(i,j), xg_nlvls(i,j)+1)
           call set_obs_def_key(obs_def, qc_count)
           call set_obs_values(obs, co_vmr, 1)
           call set_qc(obs, co_qc, num_qc)
           call set_obs_def(obs, obs_def)
!
           if ( qc_count == 1 .or. old_ob.eq.1) then
              call insert_obs_in_seq(seq, obs)
           else
              call insert_obs_in_seq(seq, obs, obs_old )
           endif
           obs_old=obs
        enddo
        if(trim(IASI_CO_retrieval_type).eq.'QOR' .or. &
        trim(IASI_CO_retrieval_type).eq.'CPSR') then 
           deallocate(Z,SV_cov,SV)
           deallocate(U_cov,UT_cov,V_cov,VT_cov)
           deallocate(rs_avg_k,rs_cov,rs_x_r,rs_x_p)       
           deallocate(rr_avg_k,rr_cov,rr_x_r,rr_x_p)       
           deallocate(ZL,ZR,ZV)

        endif
     enddo
  enddo   
  deallocate (xg_count,xg,xg_nlvls)
  deallocate (xg_lon,xg_lat,xg_twt,xg_dof)
  deallocate (xg_sec,xg_raw_err,xg_ret_err)
  deallocate (xg_raw_adj_x_r,xg_raw_adj_x_p,xg_raw_x_r,xg_raw_x_p)
  deallocate (xg_ret_adj_x_r,xg_ret_adj_x_p,xg_ret_x_r, xg_ret_x_p)
  deallocate (xg_norm,xg_nint)
  deallocate (xg_avg_k,xg_raw_cov,xg_ret_cov)
  deallocate (xg_prs,xg_prs_norm)
!
!----------------------------------------------------------------------
! Write the sequence to a file
!----------------------------------------------------------------------
 call write_obs_seq(seq, fileout)

999 continue
 close(fileid)

!-----------------------------------------------------------------------------
! Clean up
!-----------------------------------------------------------------------------
 call timestamp(string1=source,string2=revision,string3=revdate,pos='end')

end program iasi_co_profile_ascii_to_obs
