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
!
program mls_o3_profile_ascii_to_obs
!
!=============================================
! MLS O3 column obs
!=============================================
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

  use utilities_mod, only          : timestamp,                  &
                                      register_module,            &
                                      open_file,                  &
                                      close_file,                 &
                                      initialize_utilities,       &
                                      find_namelist_in_file,      &
                                      check_namelist_read,        &
                                      error_handler,              &
                                      E_ERR,                      & 
                                      E_WARN,                     & 
                                      E_MSG,                      &
                                      E_DBG

   use obs_sequence_mod, only       : obs_sequence_type,          &
                                      interactive_obs,            &
                                      write_obs_seq,              &
                                      interactive_obs_sequence,   &
                                      static_init_obs_sequence,   &
                                      init_obs_sequence,          &
                                      init_obs,                   &
                                      set_obs_values,             &
                                      set_obs_def,                &
                                      set_qc,                     &
                                      set_qc_meta_data,           &
                                      set_copy_meta_data,         &
                                      insert_obs_in_seq,          &
                                      obs_type

   use obs_def_mod, only            : set_obs_def_location,       &
                                      set_obs_def_time,           &
                                      set_obs_def_key,            &
                                      set_obs_def_error_variance, &
                                      obs_def_type,               &
                                      set_obs_def_type_of_obs

   use obs_def_mls_o3_profile_mod, only     : set_obs_def_mls_o3_profile

   use assim_model_mod, only        : static_init_assim_model

   use location_mod, only           : location_type,              &
                                      set_location

   use time_manager_mod, only       : set_date,                   &
                                      set_calendar_type,          &
                                      time_type,                  &
                                      get_time

   use obs_kind_mod, only           : QTY_O3,                     &
                                      MLS_O3_PROFILE,              &
                                      get_type_of_obs_from_menu

   use random_seq_mod, only         : random_seq_type,            &
                                      init_random_seq,            &
                                      random_uniform

   use sort_mod, only               : index_sort
   implicit none
!
! version controlled file description for error handling, do not edit
   character(len=*), parameter     :: source   = 'mls_o3_profile_ascii_to_obs.f90'
   character(len=*), parameter     :: revision = ''
   character(len=*), parameter     :: revdate  = ''
!
   integer,parameter               :: num_copies=1, num_qc=1
   integer,parameter               :: max_num_obs=1000000
   type(obs_sequence_type)         :: seq
   type(obs_type)                  :: obs
   type(obs_type)                  :: obs_old
   type(obs_def_type)              :: obs_def
   type(location_type)             :: obs_location
   type(time_type)                 :: obs_time
!
   integer                         :: obs_kind
   integer                         :: obs_key
   integer                         :: year,month,day,hour,sec
   integer                         :: iunit,io,ios,icopy,old_ob
   integer                         :: calendar_type,qc_count
   integer                         :: line_count,fileid,nlevels
   integer                         :: obs_id,yr_obs,mn_obs,dy_obs,hh_obs,mm_obs,ss_obs
   integer                         :: nlay_obs,nlev_obs,ret_nlay,tru_nlay,ilv
   integer                         :: seconds,days,which_vert
   integer                         :: seconds_last,days_last
   integer                         :: nx_model,ny_model,nz_model
   integer                         :: reject,k,l,kk,klev,kend,ilay
   integer                         :: i_min,j_min,icol,reject_ak
   integer                         :: sum_reject,sum_accept,sum_total
   integer                         :: obs_accept,obs_o3_reten_freq,obs_hno3_reten_freq
!
   integer,dimension(12)           :: days_in_month=(/ &
                                      31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31  /)
!
   real                            :: bin_beg_sec,bin_end_sec
   real                            :: lon_min,lon_max,lat_min,lat_max
   real                            :: fac_obs_error,fac_err
   real                            :: pi,rad2deg,re,level_crit
   real                            :: x_observ,y_observ,dofs
   real                            :: prs_loc,obs_sum
   real*8                          :: obs_err_var,level
!
   real*8,dimension(num_qc)        :: mls_qc
   real*8,dimension(num_copies)    :: obs_val
!
   character*129                   :: filedir,filename,fileout
   character*129                   :: copy_meta_data
   character*129                   :: qc_meta_data='MLS O3 QC index'
   character*129                   :: chr_year,chr_month,chr_day
   character*129                   :: file_name='mls_o3_profile_obs_seq'
   character*129                   :: data_type,cmd
   character*129                   :: path_model,file_model,file_in
!
   logical                         :: use_log_o3,use_log_hno3
!
! Species-specific variables
   real                              :: o3_col_obs,o3_col_obs_err,prior_err
   real                              :: lat_obs,lon_obs,dofs_obs,trop_prs
   real*8                            :: lat_obs_r8,lon_obs_r8
   real,allocatable,dimension(:,:)   :: avgk_obs,cov_obs,cov_total
   real,allocatable,dimension(:)     :: prior_obs,o3_obs,o3_obs_err
   real*8,allocatable,dimension(:)   :: avgk_obs_r8
   real*8,allocatable,dimension(:,:) :: cov_obs_r8
   real*8,allocatable,dimension(:)   :: prior_obs_r8,o3_obs_r8
   real,allocatable,dimension(:)     :: prs_obs,ret_lay_obs,tru_lay_obs
   real*8,allocatable,dimension(:)   :: prs_obs_r8
   real,allocatable,dimension(:,:)   :: lon,lat
   real,allocatable,dimension(:,:,:) :: prs_prt,prs_bas,prs_fld
   real,allocatable,dimension(:,:,:) :: tmp_prt,tmp_fld,vtmp_fld
   real,allocatable,dimension(:,:,:) :: o3_fld,qmr_fld
!
   namelist /create_mls_obs_nml/filedir,filename,fileout, &
   bin_beg_sec,bin_end_sec,fac_obs_error,use_log_o3,use_log_hno3, &
   lon_min,lon_max,lat_min,lat_max,path_model,file_model,nx_model,ny_model, &
   nz_model,obs_o3_reten_freq,obs_hno3_reten_freq
!
! Set constants
   pi=4.*atan(1.)
   rad2deg=360./(2.*pi)
   re=6371000.
   days_last=-9999.
   seconds_last=-9999.
   level_crit=500.
   sum_reject=0
   sum_accept=0
   sum_total=0
   obs_accept=0
   fac_err=0.1
!
! Record the current time, date, etc. to the logfile
   call initialize_utilities(source)
   call register_module(source,revision,revdate)
!
! Initialize the obs_sequence module
   call static_init_obs_sequence()
!
   call find_namelist_in_file("input.nml", "create_mls_obs_nml", iunit)
   read(iunit, nml = create_mls_obs_nml, iostat = io)
   call check_namelist_read(iunit, io, "create_mls_obs_nml")
!
! Record the namelist values used for the run ...
   call error_handler(E_MSG,'init_create_mls_obs','create_mls_obs_nml values are',' ',' ',' ')
   write(     *     , nml=create_mls_obs_nml)
!
! Initialize an obs_sequence structure
   call init_obs_sequence(seq, num_copies, num_qc, max_num_obs)
!
! Initialize the obs variable
   call init_obs(obs, num_copies, num_qc)
!
   do icopy =1, num_copies
      if (icopy == 1) then
         copy_meta_data='MLS O3 observation'
      else
         copy_meta_data='Truth'
      endif
      call set_copy_meta_data(seq, icopy, copy_meta_data)
   enddo
   call set_qc_meta_data(seq, 1, qc_meta_data)
!
!-------------------------------------------------------
! Read MLS O3 data
!-------------------------------------------------------
!
! Set dates and initialize qc_count
   qc_count=0
   calendar_type=3                          !Gregorian
   call set_calendar_type(calendar_type)
!
! Read model data
!   allocate(lon(nx_model,ny_model))
!   allocate(lat(nx_model,ny_model))
!   allocate(prs_prt(nx_model,ny_model,nz_model))
!   allocate(prs_bas(nx_model,ny_model,nz_model))
!   allocate(prs_fld(nx_model,ny_model,nz_model))
!   allocate(tmp_prt(nx_model,ny_model,nz_model))
!   allocate(tmp_fld(nx_model,ny_model,nz_model))
!   allocate(qmr_fld(nx_model,ny_model,nz_model))
!   allocate(o3_fld(nx_model,ny_model,nz_model))
!   file_in=trim(path_model)//'/'//trim(file_model)
!   call get_DART_diag_data(trim(file_in),'XLONG',lon,nx_model,ny_model,1,1)
!   call get_DART_diag_data(trim(file_in),'XLAT',lat,nx_model,ny_model,1,1)
!   call get_DART_diag_data(trim(file_in),'P',prs_prt,nx_model,ny_model,nz_model,1)
!   call get_DART_diag_data(trim(file_in),'PB',prs_bas,nx_model,ny_model,nz_model,1)
!   call get_DART_diag_data(trim(file_in),'T',tmp_prt,nx_model,ny_model,nz_model,1)
!   call get_DART_diag_data(trim(file_in),'QVAPOR',qmr_fld,nx_model,ny_model,nz_model,1)
!   call get_DART_diag_data(file_in,'o3',o3_fld,nx_model,ny_model,nz_model,1)
!   prs_fld(:,:,:)=prs_bas(:,:,:)+prs_prt(:,:,:)
!   tmp_fld(:,:,:)=300.+tmp_prt(:,:,:)
!   o3_fld(:,:,:)=o3_fld(:,:,:)*1.e-6
!
! Open MLS O3 binary file
   fileid=100
   write(6,*)'opening ',TRIM(filedir)//TRIM(filename)
   open(fileid,file=TRIM(filedir)//TRIM(filename), &
   form='formatted', status='old', iostat=ios)
!
! Read MLS O3
   read(fileid,*,iostat=ios) data_type, obs_id, i_min, j_min
   do while (ios == 0)
      read(fileid,*,iostat=ios) yr_obs, mn_obs, &
      dy_obs, hh_obs, mm_obs, ss_obs
      read(fileid,*,iostat=ios) lat_obs,lon_obs
      if(lon_obs.lt.0.) lon_obs=lon_obs+360.
      read(fileid,*,iostat=ios) nlay_obs,nlev_obs
      read(fileid,*,iostat=ios) ret_nlay,tru_nlay
!
! Check layer variables
      reject=0      
      if(nlay_obs.ne.ret_nlay .or. nlay_obs.ne.tru_nlay .or. &
      ret_nlay.ne.tru_nlay) reject=1
!
      allocate(prs_obs(nlay_obs))
      allocate(o3_obs(nlay_obs))
      allocate(o3_obs_err(nlay_obs))
      allocate(prior_obs(nlay_obs))
      allocate(ret_lay_obs(ret_nlay))
      allocate(tru_lay_obs(tru_nlay))
      allocate(avgk_obs(ret_nlay,tru_nlay))
      allocate(cov_obs(ret_nlay,tru_nlay))
      allocate(prior_obs_r8(nlay_obs))
      allocate(prs_obs_r8(nlay_obs))
      allocate(avgk_obs_r8(tru_nlay))
!
! MLS prs_obs, ret_lay_obs, and tru_lay_obs are all the same
      read(fileid,*,iostat=ios) prs_obs(1:nlay_obs)
      read(fileid,*,iostat=ios) ret_lay_obs(1:ret_nlay)
      read(fileid,*,iostat=ios) tru_lay_obs(1:tru_nlay)
      read(fileid,*,iostat=ios) o3_obs(1:nlay_obs)
      read(fileid,*,iostat=ios) o3_obs_err(1:nlay_obs)
      cov_obs(:,:)=0.
      do k=1,nlay_obs
         cov_obs(k,k)=o3_obs_err(k)*o3_obs_err(k)
      enddo
      read(fileid,*,iostat=ios) prior_obs(1:nlay_obs)
      read(fileid,*,iostat=ios) prior_err
      do k=1,ret_nlay
         read(fileid,*,iostat=ios) (avgk_obs(k,l),l=1,tru_nlay)
      enddo
!      read(fileid,*,iostat=ios) o3_col_obs
!      read(fileid,*,iostat=ios) o3_col_obs_err
      if(reject.eq.1) then
         deallocate(prs_obs) 
         deallocate(o3_obs) 
         deallocate(o3_obs_err) 
         deallocate(prior_obs)
         deallocate(ret_lay_obs)
         deallocate(tru_lay_obs)
         deallocate(avgk_obs)
         deallocate(cov_obs)
         deallocate(prior_obs_r8)
         deallocate(prs_obs_r8) 
         deallocate(avgk_obs_r8)
         read(fileid,*,iostat=ios) data_type, obs_id, i_min, j_min
         cycle
      endif
      print *, 'MLS O3: Completed data read'      
!
! Obs thinning test
      obs_accept=obs_accept+1
      if(obs_accept/obs_o3_reten_freq*obs_o3_reten_freq.eq.obs_accept) then
         prs_obs(:)=prs_obs(:)*100.
         prs_obs_r8(:)=prs_obs(:)
         prior_obs_r8(:)=prior_obs(:)
         lon_obs_r8=lon_obs
         lat_obs_r8=lat_obs
!      
! Loop through the profile retrievals
         do k=1,nlay_obs
            ilay=k
            if(all(avgk_obs(ilay,1:nlay_obs).eq.0.)) then
               cycle
            endif
            if(o3_obs(k).lt.0. .or. prior_obs(k).lt.0. .or. &
            isnan(o3_obs(k)) .or. isnan(prior_obs(k)) .or. &
            any(isnan(avgk_obs(k,:))) .or. all(avgk_obs(k,:).eq.0)) then
               print *, 'o3, prior, cov < 0 ',k,nlay_obs,o3_obs(k),prior_obs(k),cov_obs(k,k)
               cycle
            endif
!
! Set data for writing obs_sequence file
            sum_accept=sum_accept+1
            qc_count=qc_count+1
            avgk_obs_r8(1:nlay_obs)=avgk_obs(ilay,1:nlay_obs)         
            obs_val(:)=o3_obs(ilay)
            obs_err_var=(fac_obs_error*fac_err*sqrt(cov_obs(ilay,ilay)))**2
            mls_qc(:)=0
            obs_time=set_date(yr_obs,mn_obs,dy_obs,hh_obs,mm_obs,ss_obs)
            call get_time(obs_time, seconds, days)
!
!            which_vert=-2      ! undefined
!            which_vert=-1      ! surface
!            which_vert=1       ! level
            which_vert=2       ! pressure
!
            obs_kind = MLS_O3_PROFILE
! (0 <= lon_obs <= 360); (-90 <= lat_obs <= 90)
            klev=ilay
            kend=ilay
            level=prs_obs(ilay)
            obs_location=set_location(lon_obs_r8, lat_obs_r8, level, which_vert)
!
            call set_obs_def_type_of_obs(obs_def, obs_kind)
            call set_obs_def_location(obs_def, obs_location)
            call set_obs_def_time(obs_def, obs_time)
            call set_obs_def_error_variance(obs_def, obs_err_var)
            call set_obs_def_mls_o3_profile(qc_count, prs_obs_r8(1:nlay_obs), &
            avgk_obs_r8(1:nlay_obs), prior_obs_r8(1:nlay_obs), nlay_obs, klev, kend)
            call set_obs_def_key(obs_def, qc_count)
            call set_obs_values(obs, obs_val, 1)
            call set_qc(obs, mls_qc, num_qc)
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
      endif
      deallocate(prs_obs) 
      deallocate(o3_obs) 
      deallocate(o3_obs_err) 
      deallocate(prior_obs)
      deallocate(ret_lay_obs)
      deallocate(tru_lay_obs)
      deallocate(avgk_obs)
      deallocate(cov_obs)
      deallocate(prior_obs_r8)
      deallocate(prs_obs_r8) 
      deallocate(avgk_obs_r8)
      read(fileid,*,iostat=ios) data_type, obs_id, i_min, j_min
   enddo
!
!----------------------------------------------------------------------
! Write the sequence to a file
!----------------------------------------------------------------------
!   deallocate(lon)
!   deallocate(lat)
!   deallocate(prs_prt)
!   deallocate(prs_bas)
!   deallocate(prs_fld)
!   deallocate(tmp_prt)
!   deallocate(tmp_fld)
!   deallocate(qmr_fld)
!   deallocate(o3_fld)
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
end program mls_o3_profile_ascii_to_obs
