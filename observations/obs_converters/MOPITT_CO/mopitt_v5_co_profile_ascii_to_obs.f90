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
program mopitt_v5_co_profile_ascii_to_obs
!
!=============================================
! MOPITT CO retrieval obs
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

   use obs_def_mopitt_v5_co_profile_mod, only : set_obs_def_mopitt_v5_co_profile

   use  assim_model_mod, only : static_init_assim_model

   use location_mod, only : location_type,                 &
                         set_location

   use time_manager_mod, only : set_date, 			&
                             set_calendar_type, 	&
                             time_type, 		&
                             get_time

   use obs_kind_mod, only   : QTY_CO,                      &
                           MOPITT_V5_CO_PROFILE,         &
                           get_type_of_obs_from_menu

   use random_seq_mod, only : random_seq_type,             &
                           init_random_seq,             &
                           random_uniform

   use sort_mod, only       : index_sort
   implicit none
!
! version controlled file description for error handling, do not edit
   character(len=*), parameter :: source   = 'mopitt_v5_co_profile_ascii_to_obs.f90'
   character(len=*), parameter :: revision = ''
   character(len=*), parameter :: revdate  = ''
!
   integer,parameter                 :: num_copies=1, num_qc=1
   integer,parameter                 :: max_num_obs=1000000
   real*8,dimension(num_qc)          :: mopitt_qc
   real*8,dimension(num_copies)      :: obs_val
   type(obs_sequence_type)           :: seq
   type(obs_type)                    :: obs
   type(obs_type)                    :: obs_old
   type(obs_def_type)                :: obs_def
   type(location_type)               :: obs_location
   type(time_type)                   :: obs_time
                                   
   integer                           :: obs_kind
   integer                           :: iunit,io,ios,icopy,old_ob
   integer                           :: calendar_type,qc_count,obs_cnt
   integer                           :: fileid,i_min,j_min,klay
   integer                           :: nlay_obs,nlaym_obs,nlev_obs   
   integer                           :: seconds,days,which_vert
   integer                           :: seconds_last,days_last
   integer                           :: nx_model,ny_model,nz_model
   integer                           :: i,j,k,ilv,flg
   integer                           :: sum_reject,sum_accept,sum_total
   integer                           :: yyyy_obs,mm_obs,mn_obs,dy_obs,hh_obs,ss_obs
   integer                           :: year,month,day,hour,sec,fid
   integer                           :: kstr,nlev_obs_stg,nlay_obs_stg
!
   integer,dimension(12)             :: days_in_month=(/ &
        31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31  /)
!
   real                              :: bin_beg_sec,bin_end_sec,dofs,psfc_obs
   real                              :: lon_min,lon_max,lat_min,lat_max
   real                              :: lon_obs,lat_obs
   real                              :: fac_obs_error,fac_err
   real                              :: pi,rad2deg,re,level_crit
   real                              :: du2molpm2,log10e
   character*180                     :: filedir,filename,fileout,cmd
   character*180                     :: copy_meta_data,obser_type
   character*180                     :: qc_meta_data='MOPITT CO QC index'
   character*180                     :: chr_year,chr_month,chr_day,chr_hour
   character*180                     :: file_name='mopitt_co_profile_obs_seq'
   character*180                     :: path_model,file_model,file_in
   logical                           :: use_log_co
!
! Species-specific variables
   real                              :: prs_loc
   real                              :: co_sfc,co_sfc_err
   real                              :: prior_sfc,prior_sfc_err
   real                              :: co_col_amt,co_col_amt_err
   real,allocatable,dimension(:)     :: prs_lev,co_laym,co_err_laym
   real,allocatable,dimension(:)     :: prior_laym,prior_err_laym
   real,allocatable,dimension(:)     :: prs_obs,co_obs,prior_obs
   real,allocatable,dimension(:,:)   :: avgk_lay,cov_lay
   real,allocatable,dimension(:,:)   :: avgk_obs,cov_obs

   real*8                            :: lat_obs_r8,lon_obs_r8,level,psfc_obs_r8
   real*8                            :: obs_err_var
   real*8,allocatable,dimension(:)   :: avgk_obs_r8,prs_obs_r8,prior_obs_r8
   real,allocatable,dimension(:,:)   :: lon,lat,psfc_fld
   real,allocatable,dimension(:,:,:) :: prs_prt,prs_bas,prs_fld
   real,allocatable,dimension(:,:,:) :: tmp_prt,tmp_fld,vtmp_fld
   real,allocatable,dimension(:,:,:) :: co_fld,qmr_fld
!
   namelist /create_mopitt_obs_nml/filedir,filename,fileout,year,month,day,hour, &
        bin_beg_sec,bin_end_sec,fac_obs_error,use_log_co,lon_min,lon_max,lat_min, &
        lat_max,path_model,file_model,nx_model,ny_model,nz_model
!
! Set constants
   log10e=log10(exp(1.0))
   pi=4.*atan(1.)
   rad2deg=360./(2.*pi)
   re=6371000.
   fac_err=1.0
   days_last=-9999.
   seconds_last=-9999.
   level_crit=50000.
   sum_reject=0
   sum_accept=0
   sum_total=0
   qc_count=0
!
! Record the current time, date, etc. to the logfile
   call initialize_utilities(source)
   call register_module(source,revision,revdate)
!
! Initialize the obs_sequence module
   call static_init_obs_sequence()
!
   call find_namelist_in_file("input.nml", "create_mopitt_obs_nml", iunit)
   read(iunit, nml = create_mopitt_obs_nml, iostat = io)
   call check_namelist_read(iunit, io, "create_mopitt_obs_nml")
   write(chr_year,'(i4.4)') year
   write(chr_month,'(i2.2)') month
   write(chr_day,'(i2.2)') day
!
! Record the namelist values used for the run ...
   call error_handler(E_MSG,'init_create_mopitt_obs','create_mopitt_obs_nml values are',' ',' ',' ')
   write(     *     , nml=create_mopitt_obs_nml)
!
! Initialize an obs_sequence structure
   call init_obs_sequence(seq, num_copies, num_qc, max_num_obs)
!
! Initialize the obs variable
   call init_obs(obs, num_copies, num_qc)
!
   do icopy =1, num_copies
      if (icopy == 1) then
         copy_meta_data='MOPITT CO observation'
      else
         copy_meta_data='Truth'
      endif
      call set_copy_meta_data(seq, icopy, copy_meta_data)
   enddo
   call set_qc_meta_data(seq, 1, qc_meta_data)
!
!-------------------------------------------------------
! Read MOPITT obs
!-------------------------------------------------------
!
! Set dates and initialize qc_count
   calendar_type=3                          !Gregorian
   call set_calendar_type(calendar_type)
!
! Read model data
   allocate(lon(nx_model,ny_model))
   allocate(lat(nx_model,ny_model))
   allocate(psfc_fld(nx_model,ny_model))
   allocate(prs_prt(nx_model,ny_model,nz_model))
   allocate(prs_bas(nx_model,ny_model,nz_model))
   allocate(prs_fld(nx_model,ny_model,nz_model))
   allocate(tmp_prt(nx_model,ny_model,nz_model))
   allocate(tmp_fld(nx_model,ny_model,nz_model))
   allocate(qmr_fld(nx_model,ny_model,nz_model))
   allocate(co_fld(nx_model,ny_model,nz_model))
   file_in=trim(path_model)//'/'//trim(file_model)
   call get_DART_diag_data(trim(file_in),'XLONG',lon,nx_model,ny_model,1,1)
   call get_DART_diag_data(trim(file_in),'XLAT',lat,nx_model,ny_model,1,1)
   call get_DART_diag_data(trim(file_in),'PSFC',psfc_fld,nx_model,ny_model,1,1)   
   call get_DART_diag_data(trim(file_in),'P',prs_prt,nx_model,ny_model,nz_model,1)
   call get_DART_diag_data(trim(file_in),'PB',prs_bas,nx_model,ny_model,nz_model,1)
   call get_DART_diag_data(trim(file_in),'T',tmp_prt,nx_model,ny_model,nz_model,1)
   call get_DART_diag_data(trim(file_in),'QVAPOR',qmr_fld,nx_model,ny_model,nz_model,1)
   call get_DART_diag_data(file_in,'co',co_fld,nx_model,ny_model,nz_model,1)
   prs_fld(:,:,:)=prs_bas(:,:,:)+prs_prt(:,:,:)
   tmp_fld(:,:,:)=300.+tmp_prt(:,:,:)
   co_fld(:,:,:)=co_fld(:,:,:)*1.e-6

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
! Open MOPITT CO binary file
   fid=100   
   write(6,*)'opening ',TRIM(filedir)//TRIM(filename)
   open(fileid,file=TRIM(filedir)//TRIM(filename),form='formatted', &
   status='old',iostat=ios)

! Read MOPITT CO
   read(fileid,*,iostat=ios) obser_type,obs_cnt,i_min,j_min
   do while (ios == 0)
      sum_total=sum_total+1
      read(fileid,*) yyyy_obs,mn_obs,dy_obs,hh_obs,mm_obs,ss_obs
      sec=hh_obs*3600 + mm_obs*60 + ss_obs
      read(fileid,*) lat_obs,lon_obs
      if(lon_obs.lt.0) lon_obs=lon_obs+360.
      read(fileid,*) nlay_obs,nlev_obs
      nlaym_obs=nlay_obs-1
!
      allocate(prs_lev(nlay_obs))
      allocate(avgk_lay(nlay_obs,nlay_obs))
      allocate(co_laym(nlaym_obs))
      allocate(co_err_laym(nlaym_obs))
      allocate(prior_laym(nlaym_obs))
      allocate(prior_err_laym(nlaym_obs))
      allocate(cov_lay(nlay_obs,nlay_obs))
!      
      read(fileid,*) dofs
      read(fileid,*) psfc_obs
      read(fileid,*) prs_lev(1:nlev_obs-1)
      do k=1,nlay_obs
         read(fileid,*) avgk_lay(k,1:nlay_obs)
      enddo
      read(fileid,*) co_sfc
      read(fileid,*) co_laym(1:nlay_obs-1)
      read(fileid,*) co_sfc_err
      read(fileid,*) co_err_laym(1:nlay_obs-1)
      read(fileid,*) prior_sfc
      read(fileid,*) prior_laym(1:nlay_obs-1)
      read(fileid,*) prior_sfc_err
      read(fileid,*) prior_err_laym(1:nlay_obs-1)
      do k=1,nlay_obs
         read(fileid,*) cov_lay(k,1:nlay_obs)
      enddo
      read(fileid,*) co_col_amt
      read(fileid,*) co_col_amt_err
!
! Process observations
      if(psfc_obs.gt.prs_lev(1)) then
         kstr=1
      else
         do ilv=1,nlay_obs-1
            if(psfc_obs.le.prs_lev(ilv) .and. &
               psfc_obs.gt.prs_lev(ilv+1)) then
               kstr=ilv+1
               exit
            endif
         enddo
      endif
!
! Place observations on collapsed vertical grid
      nlev_obs_stg=nlev_obs-kstr+1
      nlay_obs_stg=nlev_obs_stg-1
!
      allocate(prs_obs(nlev_obs_stg))
      allocate(co_obs(nlay_obs_stg))
      allocate(prior_obs(nlay_obs_stg))
      allocate(avgk_obs(nlay_obs_stg,nlay_obs_stg))
      allocate(cov_obs(nlay_obs_stg,nlay_obs_stg))
      allocate(avgk_obs_r8(nlay_obs_stg))
      allocate(prs_obs_r8(nlev_obs_stg))
      allocate(prior_obs_r8(nlay_obs_stg))
!
      prs_obs(1)=psfc_obs
      prs_obs(2:nlev_obs_stg)=prs_lev(kstr:nlev_obs-1)
      co_obs(1)=co_sfc
      co_obs(2:nlay_obs_stg)=co_laym(kstr:nlaym_obs)
      prior_obs(1)=prior_sfc
      prior_obs(2:nlay_obs_stg)=prior_laym(kstr:nlaym_obs)
      avgk_obs(1:nlay_obs_stg,1:nlay_obs_stg)=avgk_lay(kstr:nlay_obs,kstr:nlay_obs)
      cov_obs(1:nlay_obs_stg,1:nlay_obs_stg)=cov_lay(kstr:nlay_obs,kstr:nlay_obs)
!
! QA/QC checks      
      flg=0
      do k=1,nlay_obs_stg
         if(isnan(cov_obs(k,k)) .or. cov_obs(k,k).lt.0.) then
            flg=1
            print *, 'mopitt error variance is NaN or negative ',cov_obs(k,k)
            exit
         endif
      enddo
      if(flg.eq.1) then
         deallocate(prs_lev)
         deallocate(avgk_lay)
         deallocate(co_laym)
         deallocate(co_err_laym)
         deallocate(prior_laym)
         deallocate(prior_err_laym)
         deallocate(cov_lay)
         deallocate(prs_obs)
         deallocate(co_obs)
         deallocate(prior_obs)
         deallocate(avgk_obs)
         deallocate(cov_obs)
         deallocate(avgk_obs_r8)
         deallocate(prs_obs_r8)
         deallocate(prior_obs_r8)
         read(fileid,*,iostat=ios) obser_type, obs_cnt
         cycle
      endif
!
      prs_obs(:)=prs_obs(:)*100.
      prs_obs_r8(:)=prs_obs(:)/100.
      prior_obs_r8(:)=prior_obs(:)
      psfc_obs_r8=psfc_obs
!      
      lon_obs_r8=lon_obs
      lat_obs_r8=lat_obs
!
!--------------------------------------------------------
! Find model CO profile corresponding to the observation
! kend is the MOPITT index for the top of the model.      
!--------------------------------------------------------
!      call get_model_profile(prf_model,nz_model, &
!      prs_obs,prs_fld(i_min,j_min,:),tmp_fld(i_min,j_min,:), &
!      qmr_fld(i_min,j_min,:),co_fld(i_min,j_min,:),psfc_fld(i_min,j_min), &
!      nlev_obs,avgk_obs,prior_obs,kend)
!
! Loop through vertical grid (MOPITT CO is bottom to top)
      do ilv=1,nlay_obs_stg
         avgk_obs_r8(1:nlay_obs_stg)=avgk_obs(ilv,1:nlay_obs_stg)
!
         if(prior_obs_r8(ilv).lt.0. .or. co_obs(ilv).lt.0.) then
            print *, 'APM: MOPITT Retrieval or Prior is negative. Layer: ',ilv
            print *, 'APM: Prior ',prior_obs_r8(ilv),'Avgk: ',avgk_obs_r8(ilv)
            cycle
         endif
!
!--------------------------------------------------------
! Find vertical location
!--------------------------------------------------------
!         call vertical_locate(prs_loc,prs_obs,nlev_obs,prf_model,nlay_obs,kend)
!         level=prs_loc
!
! Process accepted observations
         sum_accept=sum_accept+1
         qc_count=qc_count+1
!
! Obs value is the total vertical column      
         obs_val(:)=co_obs(ilv)
         obs_err_var=(fac_obs_error*fac_err*sqrt(cov_obs(ilv,ilv)))**2.
         mopitt_qc(:)=0      
         obs_time=set_date(yyyy_obs,mn_obs,dy_obs,hh_obs,mm_obs,ss_obs)
         call get_time(obs_time, seconds, days)
!
!         which_vert=-2      ! undefined
!         which_vert=-1      ! surface
!         which_vert=1       ! level
         which_vert=2       ! pressure surface
!
         obs_kind = MOPITT_V5_CO_PROFILE
! (0 <= lon_obs <= 360); (-90 <= lat_obs <= 90)
         klay=ilv
         level=prs_obs(ilv)   
         obs_location=set_location(lon_obs_r8, lat_obs_r8, level, which_vert)
!
         call set_obs_def_type_of_obs(obs_def, obs_kind)
         call set_obs_def_location(obs_def, obs_location)
         call set_obs_def_time(obs_def, obs_time)
         call set_obs_def_error_variance(obs_def, obs_err_var)
         call set_obs_def_mopitt_v5_co_profile(qc_count, prs_obs_r8, avgk_obs_r8, prior_obs_r8, nlay_obs_stg,klay)
         call set_obs_def_key(obs_def, qc_count)
         call set_obs_values(obs, obs_val, 1)
         call set_qc(obs, mopitt_qc, num_qc)
         call set_obs_def(obs, obs_def)
!
         old_ob=0
         if(days.lt.days_last) then
            old_ob=1
         elseif(days.eq.days_last .and. seconds.lt.seconds_last) then
            old_ob=1
         endif
         if(old_ob.eq.0) then
            days_last=days
            seconds_last=seconds
         endif
         if ( qc_count == 1 .or. old_ob.eq.1) then
            call insert_obs_in_seq(seq, obs)
         else
            call insert_obs_in_seq(seq, obs, obs_old )
         endif
         obs_old=obs
      enddo
!
! Read next data point
      deallocate(prs_lev)
      deallocate(avgk_lay)
      deallocate(co_laym)
      deallocate(co_err_laym)
      deallocate(prior_laym)
      deallocate(prior_err_laym)
      deallocate(cov_lay)
      deallocate(prs_obs)
      deallocate(co_obs)
      deallocate(prior_obs)
      deallocate(avgk_obs)
      deallocate(cov_obs)
      deallocate(avgk_obs_r8)
      deallocate(prs_obs_r8)
      deallocate(prior_obs_r8)
      read(fileid,*,iostat=ios) obser_type, obs_cnt
   enddo !ios
!
!----------------------------------------------------------------------
! Write the sequence to a file
!----------------------------------------------------------------------
   deallocate(lon)
   deallocate(lat)
   deallocate(psfc_fld)
   deallocate(prs_prt)
   deallocate(prs_bas)
   deallocate(prs_fld)
   deallocate(tmp_prt)
   deallocate(tmp_fld)
   deallocate(qmr_fld)
   deallocate(co_fld)
!
   print *, 'total obs ',sum_total
   print *, 'accepted ',sum_accept
   print *, 'rejected ',sum_reject
   call timestamp(string1=source,string2=revision,string3=revdate,pos='end')
   call write_obs_seq(seq, trim(fileout))
   close(fileid)
!
! Remove obs_seq if empty
   cmd='rm -rf '//trim(fileout)
   if(qc_count.eq.0) then
      call execute_command_line(trim(cmd))
   endif   
!
end program
