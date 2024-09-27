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
program omi_no2_total_col_ascii_to_obs
!
!=============================================
! OMI NO2 total column obs
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

 use        types_mod, only       : r8  
 use time_manager_mod, only       : set_date,                   &
                                    set_calendar_type,          &
                                    time_type,                  &
                                    get_time
 use    utilities_mod, only       : timestamp,                  &
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
 use     location_mod, only       : location_type,              &
                                    set_location

 use obs_sequence_mod, only       : obs_sequence_type,          &
                                    write_obs_seq,              &
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
 use      obs_def_mod, only       : set_obs_def_location,       &
                                    set_obs_def_time,           &
                                    set_obs_def_key,            &
                                    set_obs_def_error_variance, &
                                    obs_def_type,               &
                                    set_obs_def_type_of_obs
 use obs_def_omi_no2_total_col_mod, only : set_obs_def_omi_no2_total_col
 use     obs_kind_mod, only       : OMI_NO2_TOTAL_COL

 implicit none
!
! version controlled file description for error handling, do not edit
   character(len=*), parameter     :: source   = 'omi_no2_total_col_ascii_to_obs.f90'
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
   integer                         :: nlay_obs,nlev_obs,ilv
   integer                         :: seconds,days,which_vert
   integer                         :: seconds_last,days_last
   integer                         :: nx_model,ny_model,nz_model
   integer                         :: reject,k,kend
   integer                         :: i_min,j_min
   integer                         :: sum_reject,sum_accept,sum_total
!
   integer,dimension(12)           :: days_in_month=(/ &
                                      31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31  /)
!
   real                            :: bin_beg_sec,bin_end_sec
   real                            :: lon_min,lon_max,lat_min,lat_max
   real                            :: fac_obs_error,fac_err
   real                            :: pi,rad2deg,re,level_crit
   real                            :: x_observ,y_observ,dofs
   real(r8)                        :: obs_err_var,level
!
   real(r8),dimension(num_qc)      :: omi_qc
   real(r8),dimension(num_copies)  :: obs_val
!
   character(len=129)              :: filedir,filename,fileout
   character(len=129)              :: copy_meta_data
   character(len=129)              :: qc_meta_data='OMI NO2 QC index'
   character(len=129)              :: chr_year,chr_month,chr_day
   character(len=129)              :: file_name='omi_no2_total_col_obs_seq'
   character(len=129)              :: data_type,cmd
   character(len=129)              :: path_model,file_model,file_in
!
   logical                         :: use_log_o3,use_log_no2,use_log_so2,use_log_hcho
!
! Species-specific variables
   real                            :: amfstrat,amfstrat_clr,amfstrat_cld
   real                            :: amftrop,amftrop_clr,amftrop_cld
   real                            :: cld_frac,cld_prs,cld_rad_frac
   real                            :: col_amt,col_amt_err
   real                            :: col_amt_trop,col_amt_trop_err
   real                            :: slnt_col_amt,slnt_col_amt_err
   real                            :: prs_trop,zenang,obs_sum
   real(r8)                        :: prs_trop_r8
   real                            :: prs_loc
   real                            :: lat_obs,lon_obs
   real(r8)                        :: lat_obs_r8,lon_obs_r8
   real,allocatable,dimension(:)   :: scat_wt,prs_obs
   real(r8),allocatable,dimension(:) :: scat_wt_r8,prs_obs_r8
   real,allocatable,dimension(:)   :: prf_model
   real,allocatable,dimension(:,:)     :: lon,lat,psfc_fld
   real,allocatable,dimension(:,:,:)   :: prs_prt,prs_bas,prs_fld
   real,allocatable,dimension(:,:,:)   :: tmp_prt,tmp_fld,vtmp_fld
   real,allocatable,dimension(:,:,:)   :: no2_fld,qmr_fld
!
   namelist /create_omi_obs_nml/filedir,filename,fileout,year,month,day,hour, &
   bin_beg_sec,bin_end_sec,fac_obs_error,use_log_o3,use_log_no2,use_log_so2, &
   use_log_hcho,lon_min,lon_max,lat_min,lat_max, &
   path_model,file_model,nx_model,ny_model,nz_model
!
! Set constants
   pi=4.*atan(1.)
   rad2deg=360./(2.*pi)
   re=6371000.
   days_last=-9999.
   seconds_last=-9999.
   level_crit=50000.
   sum_reject=0
   sum_accept=0
   sum_total=0
   fac_err=1.
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
         copy_meta_data='OMI NO2 observation'
      else
         copy_meta_data='Truth'
      endif
      call set_copy_meta_data(seq, icopy, copy_meta_data)
   enddo
   call set_qc_meta_data(seq, 1, qc_meta_data)
!
!-------------------------------------------------------
! Read OMI NO2 data
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
!   allocate(no2_fld(nx_model,ny_model,nz_model))
!   file_in=trim(path_model)//'/'//trim(file_model)
!   call get_DART_diag_data(trim(file_in),'XLONG',lon,nx_model,ny_model,1,1)
!   call get_DART_diag_data(trim(file_in),'XLAT',lat,nx_model,ny_model,1,1)
!   call get_DART_diag_data(trim(file_in),'PSFC',psfc_fld,nx_model,ny_model,1,1)   
!   call get_DART_diag_data(trim(file_in),'P',prs_prt,nx_model,ny_model,nz_model,1)
!   call get_DART_diag_data(trim(file_in),'PB',prs_bas,nx_model,ny_model,nz_model,1)
!   call get_DART_diag_data(trim(file_in),'T',tmp_prt,nx_model,ny_model,nz_model,1)
!   call get_DART_diag_data(trim(file_in),'QVAPOR',qmr_fld,nx_model,ny_model,nz_model,1)
!   call get_DART_diag_data(file_in,'no2',no2_fld,nx_model,ny_model,nz_model,1)
!   prs_fld(:,:,:)=prs_bas(:,:,:)+prs_prt(:,:,:)
!   tmp_fld(:,:,:)=300.+tmp_prt(:,:,:)
!   no2_fld(:,:,:)=no2_fld(:,:,:)*1.e-6
!
! Open TROPOMI NO2 binary file
   fileid=100
   write(6,*)'opening ',TRIM(TRIM(filedir)//TRIM(filename))
   open(unit=fileid,file=TRIM(TRIM(filedir)//TRIM(filename)), &
   form='formatted', status='old', iostat=ios)
!
! Read OMI NO2
   line_count = 0
   read(fileid,*,iostat=ios) data_type, obs_id, i_min, j_min
!   print *, trim(data_type), obs_id, i_min, j_min
   do while (ios == 0)
      sum_total=sum_total+1
      read(fileid,*,iostat=ios) yr_obs, mn_obs, &
      dy_obs, hh_obs, mm_obs, ss_obs
      read(fileid,*,iostat=ios) lat_obs,lon_obs
      if(lon_obs.lt.0.) lon_obs=lon_obs+360.
      read(fileid,*,iostat=ios) nlay_obs,nlev_obs
      read(fileid,*,iostat=ios) amfstrat,amfstrat_clr,amfstrat_cld
      read(fileid,*,iostat=ios) amftrop,amftrop_clr,amftrop_cld
      read(fileid,*,iostat=ios) cld_frac,cld_prs,cld_rad_frac
      read(fileid,*,iostat=ios) col_amt,col_amt_err
      read(fileid,*,iostat=ios) col_amt_trop,col_amt_trop_err
      read(fileid,*,iostat=ios) slnt_col_amt,slnt_col_amt_err
      read(fileid,*,iostat=ios) prs_trop,zenang
      allocate(prs_obs(nlev_obs))
      allocate(scat_wt(nlay_obs))
      allocate(prs_obs_r8(nlev_obs))
      allocate(scat_wt_r8(nlay_obs))
      allocate(prf_model(nlay_obs))
      read(fileid,*,iostat=ios) scat_wt(1:nlay_obs)
      read(fileid,*,iostat=ios) prs_obs(1:nlev_obs)
      scat_wt_r8(:)=scat_wt(:)
      prs_obs(:)=prs_obs(:)*100.
      prs_obs_r8(:)=prs_obs(:)/100.
      prs_trop=prs_trop*100.
      prs_trop_r8=prs_trop/100.
      lon_obs_r8=lon_obs
      lat_obs_r8=lat_obs
!
! OMI vertical grid is bottom to top
      kend=0
      do k=1,nlay_obs
         if((prs_obs(k)+prs_obs(k+1))/2..lt.prs_trop) then
            kend=k
            exit
         endif
      enddo
!
! Process accepted observations
      sum_accept=sum_accept+1
!
! Set data for writing obs_sequence file
      qc_count=qc_count+1
!
! Obs value is the tropospheric slant column
! scd = vcd * amf      
      obs_val(:)=col_amt_trop*amftrop
      obs_err_var=(fac_obs_error*fac_err*col_amt_trop_err*amftrop)**2.
      omi_qc(:)=0
      obs_time=set_date(yr_obs,mn_obs,dy_obs,hh_obs,mm_obs,ss_obs)
      call get_time(obs_time, seconds, days)
!
      which_vert=-2      ! undefined
!      which_vert=-1      ! surface
!      which_vert=1       ! level
!      which_vert=2       ! pressure surface
!
      obs_kind = OMI_NO2_TOTAL_COL
! (0 <= lon_obs <= 360); (-90 <= lat_obs <= 90) 
      level=0.
      obs_location=set_location(lon_obs_r8, lat_obs_r8, level, which_vert)
!
      call set_obs_def_type_of_obs(obs_def, obs_kind)
      call set_obs_def_location(obs_def, obs_location)
      call set_obs_def_time(obs_def, obs_time)
      call set_obs_def_error_variance(obs_def, obs_err_var)
      call set_obs_def_omi_no2_total_col(qc_count, prs_obs_r8, scat_wt_r8, prs_trop_r8, kend, nlay_obs)
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
!      print *, 'APM: ',qc_count,days,seconds
      if ( qc_count == 1 .or. old_ob.eq.1) then
         call insert_obs_in_seq(seq, obs)
      else
         call insert_obs_in_seq(seq, obs, obs_old )
      endif
      obs_old=obs
      deallocate(scat_wt)
      deallocate(prs_obs) 
      deallocate(scat_wt_r8)
      deallocate(prs_obs_r8) 
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
!   deallocate(no2_fld)
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
end program omi_no2_total_col_ascii_to_obs
