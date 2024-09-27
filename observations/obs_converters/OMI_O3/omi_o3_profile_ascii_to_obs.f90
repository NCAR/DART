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
program omi_o3_profile_ascii_to_obs
!
!=============================================
! OMI O3 profile obs
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

   use obs_def_omi_o3_profile_mod, only     : set_obs_def_omi_o3_profile

   use assim_model_mod, only        : static_init_assim_model

   use location_mod, only           : location_type,              &
                                      set_location

   use time_manager_mod, only       : set_date,                   &
                                      set_calendar_type,          &
                                      time_type,                  &
                                      get_time

   use obs_kind_mod, only           : QTY_O3,                     &
                                      OMI_O3_PROFILE,              &
                                      get_type_of_obs_from_menu

   use random_seq_mod, only         : random_seq_type,            &
                                      init_random_seq,            &
                                      random_uniform

   use sort_mod, only               : index_sort
   implicit none
!
! version controlled file description for error handling, do not edit
   character(len=*), parameter     :: source   = 'omi_o3_profile_ascii_to_obs.f90'
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
   integer                         :: nlay_obs,nlev_obs,ndim_cov,ilv
   integer                         :: nlay_adj,nlev_adj
   integer                         :: seconds,days,which_vert
   integer                         :: seconds_last,days_last
   integer                         :: nx_model,ny_model,nz_model
   integer                         :: reject,i,j,k,kk,icnt,klev,kend
   integer                         :: i_min,j_min,flg
   integer                         :: sum_reject,sum_accept,sum_total
   integer                         :: obs_o3_reten_freq,obs_no2_reten_freq,obs_so2_reten_freq, &
                                      obs_hcho_reten_freq
!
   integer,dimension(12)           :: days_in_month=(/ &
                                      31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31  /)
!
   real                            :: bin_beg_sec,bin_end_sec
   real                            :: lon_min,lon_max,lat_min,lat_max
   real                            :: fac_obs_error,fac,fac_err
   real                            :: pi,rad2deg,re,level_crit
   real                            :: x_observ,y_observ,dofs_obs
   real*8                          :: obs_err_var,level
   real                            :: prs_trop, obs_sum
   real                            :: omi_sum_fl,omi_sum_pr,omi_frac
   real                            :: du2molpm2
!
   real*8,dimension(num_qc)        :: omi_qc
   real*8,dimension(num_copies)    :: obs_val
!
   character*129                   :: filedir,filename,fileout
   character*129                   :: copy_meta_data
   character*129                   :: qc_meta_data='OMI O3 QC index'
   character*129                   :: chr_year,chr_month,chr_day
   character*129                   :: file_name='omi_o3_profile_obs_seq'
   character*129                   :: data_type,cmd
   character*129                   :: path_model,file_model,file_in
!
   logical                         :: use_log_o3,use_log_no2,use_log_so2,use_log_hcho
!
! Species-specific variables
   real                            :: prs_loc,prs_mdl_top
   real,allocatable,dimension(:)   :: prs_obs,o3_obs,prior_obs,prior_err_obs
   real,allocatable,dimension(:)   :: cov_obs_tmp,cov_prior_obs_tmp
   real,allocatable,dimension(:,:) :: avgk_obs,cov_obs,cov_prior_obs
   real*8,allocatable,dimension(:) :: avgk_obs_r8,prs_obs_r8,prior_obs_r8
   real,allocatable,dimension(:)   :: prf_model
   real                            :: lat_obs,lon_obs
   real*8                          :: lat_obs_r8,lon_obs_r8
   real,allocatable,dimension(:,:)     :: lon,lat,psfc_fld
   real,allocatable,dimension(:,:,:)   :: prs_prt,prs_bas,prs_fld
   real,allocatable,dimension(:,:,:)   :: tmp_prt,tmp_fld,vtmp_fld
   real,allocatable,dimension(:,:,:)   :: o3_fld,qmr_fld
!
   namelist /create_omi_obs_nml/filedir,filename,fileout,year,month,day,hour, &
   bin_beg_sec,bin_end_sec,fac_obs_error,use_log_o3,use_log_no2,use_log_so2,use_log_hcho, &
   lon_min,lon_max,lat_min,lat_max,path_model,file_model,nx_model,ny_model,nz_model, &
   obs_o3_reten_freq,obs_no2_reten_freq,obs_so2_reten_freq,obs_hcho_reten_freq
!
! Set constants
   du2molpm2=4.4615e-4
   pi=4.*atan(1.)
   rad2deg=360./(2.*pi)
   re=6371000.
   fac_err=0.3
   days_last=-9999.
   seconds_last=-9999.
   prs_trop=-9999.
   level_crit=50000.
   sum_reject=0
   sum_accept=0
   sum_total=0
   prs_mdl_top=15000. ! (Pa)
!
! Record the current time, date, etc. to the logfile
   call initialize_utilities(source)
   call register_module(source,revision,revdate)
!
! Initialize the obs_sequence module
   call static_init_obs_sequence()
!
   call find_namelist_in_file("input.nml", "create_omi_obs_nml", iunit)
   read(iunit, nml = create_omi_obs_nml, iostat = io)
   call check_namelist_read(iunit, io, "create_omi_obs_nml")
   write(chr_year,'(i4.4)') year
   write(chr_month,'(i2.2)') month
   write(chr_day,'(i2.2)') day
!
! Record the namelist values used for the run ...
   call error_handler(E_MSG,'init_create_omi_obs','create_omi_obs_nml values are',' ',' ',' ')
   write(     *     , nml=create_omi_obs_nml)
!
! Initialize an obs_sequence structure
   call init_obs_sequence(seq, num_copies, num_qc, max_num_obs)
!
! Initialize the obs variable
   call init_obs(obs, num_copies, num_qc)
!
   do icopy =1, num_copies
      if (icopy == 1) then
         copy_meta_data='OMI O3 observation'
      else
         copy_meta_data='Truth'
      endif
      call set_copy_meta_data(seq, icopy, copy_meta_data)
   enddo
   call set_qc_meta_data(seq, 1, qc_meta_data)
!
!-------------------------------------------------------
! Read OMI O3 data
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
!   allocate(psfc_fld(nx_model,ny_model))
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
!   call get_DART_diag_data(trim(file_in),'PSFC',psfc_fld,nx_model,ny_model,1,1)   
!   call get_DART_diag_data(trim(file_in),'P',prs_prt,nx_model,ny_model,nz_model,1)
!   call get_DART_diag_data(trim(file_in),'PB',prs_bas,nx_model,ny_model,nz_model,1)
!   call get_DART_diag_data(trim(file_in),'T',tmp_prt,nx_model,ny_model,nz_model,1)
!   call get_DART_diag_data(trim(file_in),'QVAPOR',qmr_fld,nx_model,ny_model,nz_model,1)
!   call get_DART_diag_data(file_in,'o3',o3_fld,nx_model,ny_model,nz_model,1)
!   prs_fld(:,:,:)=prs_bas(:,:,:)+prs_prt(:,:,:)
!   tmp_fld(:,:,:)=300.+tmp_prt(:,:,:)
!   o3_fld(:,:,:)=o3_fld(:,:,:)*1.e-6
!
! Open OMI O3 binary file
   fileid=100
   write(6,*)'opening ',TRIM(filedir)//TRIM(filename)
   open(unit=fileid,file=TRIM(filedir)//TRIM(filename),form='formatted', status='old', iostat=ios)
!
! Read OMI O3
   line_count = 0
   read(fileid,*,iostat=ios) data_type, obs_id, i_min, j_min
   do while (ios == 0)
      sum_total=sum_total+1
      read(fileid,*,iostat=ios) yr_obs, mn_obs, &
      dy_obs, hh_obs, mm_obs, ss_obs
      read(fileid,*,iostat=ios) lat_obs,lon_obs
      if(lon_obs.lt.0.) lon_obs=lon_obs+360.
      read(fileid,*,iostat=ios) nlay_obs,nlev_obs,ndim_cov
      read(fileid,*,iostat=ios) dofs_obs
      allocate(prs_obs(nlev_obs))
      allocate(o3_obs(nlay_obs))
      allocate(prior_obs(nlay_obs))
      allocate(prior_err_obs(nlay_obs))
      allocate(avgk_obs(nlay_obs,nlay_obs))
      allocate(cov_obs(nlay_obs,nlay_obs))
      allocate(cov_prior_obs(nlay_obs,nlay_obs))
      allocate(cov_obs_tmp(ndim_cov))
      allocate(cov_prior_obs_tmp(ndim_cov))
      allocate(prs_obs_r8(nlev_obs))
      allocate(prior_obs_r8(nlay_obs))
      allocate(avgk_obs_r8(nlay_obs))
      allocate(prf_model(nlay_obs))
      read(fileid,*,iostat=ios) prs_obs(1:nlev_obs)
      read(fileid,*,iostat=ios) o3_obs(1:nlay_obs)
      read(fileid,*,iostat=ios) prior_obs(1:nlay_obs)
      read(fileid,*,iostat=ios) prior_err_obs(1:nlay_obs)
      do k=1,nlay_obs
         read(fileid,*,iostat=ios) avgk_obs(1:nlay_obs,k)
      enddo
      read(fileid,*,iostat=ios) cov_obs_tmp(1:ndim_cov)
      read(fileid,*,iostat=ios) cov_prior_obs_tmp(1:ndim_cov)
!
      icnt=0
      do i=1,nlay_obs
         do j=1,i
            icnt=icnt+1
            cov_obs(i,j)=cov_obs_tmp(icnt)
            cov_prior_obs(i,j)=cov_prior_obs_tmp(icnt)
            if(i.ne.j) then
               cov_obs(j,i)=cov_obs(i,j)
               cov_prior_obs(j,i)=cov_obs(i,j)
            endif
         enddo
      enddo
!
! Adjust OMI data for low surface pressures (missing levels near surface)
      nlev_adj=nlev_obs
      nlay_adj=nlay_obs
      if(any(prs_obs(:).lt.0.)) then      
         do k=1,nlev_obs
            kk=nlev_obs-k+1
            if(prs_obs(kk).lt.0.) then
               nlev_adj=nlev_adj-1
               nlay_adj=nlay_adj-1
            else
               exit
            endif
         enddo
      endif
      flg=0
      do i=1,nlay_adj
         if(isnan(cov_obs(i,i)) .or. cov_obs(i,i).lt.0.) then
            flg=1
!            print *, 'omi_variance is NaN or negative ',cov_obs(i,i)
            exit
         endif
      enddo
      if(flg.eq.1) then
         deallocate(prs_obs,o3_obs,prior_obs)
         deallocate(prior_err_obs,avgk_obs,cov_obs)
         deallocate(cov_prior_obs,cov_obs_tmp)
         deallocate(cov_prior_obs_tmp,prs_obs_r8)
         deallocate(prior_obs_r8,avgk_obs_r8) 
         deallocate(prf_model) 
         read(fileid,*,iostat=ios) data_type, obs_id, i_min, j_min
         cycle
      endif

      prs_obs(1:nlev_adj)=prs_obs(1:nlev_adj)*100.
      prs_obs_r8(1:nlev_adj)=prs_obs(1:nlev_adj)/100.
      prior_obs_r8(1:nlay_adj)=prior_obs(1:nlay_adj)
      
      lon_obs_r8=lon_obs
      lat_obs_r8=lat_obs
!
!--------------------------------------------------------
! Find model O3 profile corresponding to the observation
! kend is the OMI index for the top of the model.      
!--------------------------------------------------------
!
! Loop through vertical grid (OMI O3 is top to bottom)
      reject=0
      do ilv=1,nlay_adj
         avgk_obs_r8(1:nlay_adj)=avgk_obs(ilv,1:nlay_adj)
!
! Check whether observation is above model vertical grid
         if(prs_obs(ilv).le.prs_mdl_top) cycle
!
! Check whether observation is below the surface based on OMI          
         if(prior_obs_r8(ilv).lt.0. .or. o3_obs(ilv).lt.0.) then
            print *, 'APM: OMI Retrieval or Prior is negative. Layer may be below surface. Layer: ',ilv
            print *, 'APM: Prior ',prior_obs_r8(ilv),'Avgk: ',avgk_obs_r8(ilv)
            cycle
         endif
!
! Process observations
!
! Obs value is the total vertical column      
         obs_val(:)=o3_obs(ilv)
!
! Assign observation error
         if(prior_obs_r8(ilv).gt.0. .and. o3_obs(ilv).gt.0.) then
!            obs_err_var=(fac_obs_error*fac_err*o3_obs(ilv))**2.
            obs_err_var=(fac_obs_error*o3_obs(ilv)*sqrt(cov_obs(ilv,ilv)))**2.
         else
            obs_err_var=0.
            cycle
         endif
         sum_accept=sum_accept+1
         qc_count=qc_count+1
         
         omi_qc(:)=0      
         obs_time=set_date(yr_obs,mn_obs,dy_obs,hh_obs,mm_obs,ss_obs)
         call get_time(obs_time, seconds, days)
!
!         which_vert=-2      ! undefined
!         which_vert=-1      ! surface
!         which_vert=1       ! level
         which_vert=2       ! pressure surface
!
         obs_kind = OMI_O3_PROFILE
! (0 <= lon_obs <= 360); (-90 <= lat_obs <= 90)
         klev=ilv
         kend=ilv
         level=prs_obs(ilv)
         obs_location=set_location(lon_obs_r8, lat_obs_r8, level, which_vert)
!
         call set_obs_def_type_of_obs(obs_def, obs_kind)
         call set_obs_def_location(obs_def, obs_location)
         call set_obs_def_time(obs_def, obs_time)
         call set_obs_def_error_variance(obs_def, obs_err_var)
         call set_obs_def_omi_o3_profile(qc_count, prs_obs_r8(1:nlay_adj+1), &
         avgk_obs_r8(1:nlay_adj), prior_obs_r8(1:nlay_adj), nlay_adj, klev, kend)
         call set_obs_def_key(obs_def, qc_count)
         call set_obs_values(obs, obs_val, 1)
         call set_qc(obs, omi_qc, num_qc)
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
      deallocate(prs_obs,o3_obs,prior_obs)
      deallocate(prior_err_obs,avgk_obs,cov_obs)
      deallocate(cov_prior_obs,cov_obs_tmp)
      deallocate(cov_prior_obs_tmp,prs_obs_r8)
      deallocate(prior_obs_r8,avgk_obs_r8) 
      deallocate(prf_model)
      read(fileid,*,iostat=ios) data_type, obs_id, i_min, j_min
   enddo   
!
!----------------------------------------------------------------------
! Write the sequence to a file
!----------------------------------------------------------------------
!   deallocate(lon)
!   deallocate(lat)
!   deallocate(psfc_fld)
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
end program omi_o3_profile_ascii_to_obs
