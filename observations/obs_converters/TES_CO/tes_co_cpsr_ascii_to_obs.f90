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
program tes_co_cpsr_ascii_to_obs
!
!=============================================
! TES CO column obs
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

   use obs_def_tes_co_cpsr_mod, only     : set_obs_def_tes_co_cpsr

   use assim_model_mod, only        : static_init_assim_model

   use location_mod, only           : location_type,              &
                                      set_location

   use time_manager_mod, only       : set_date,                   &
                                      set_calendar_type,          &
                                      time_type,                  &
                                      get_time

   use obs_kind_mod, only           : QTY_CO,                     &
                                      TES_CO_CPSR,              &
                                      get_type_of_obs_from_menu

   use random_seq_mod, only         : random_seq_type,            &
                                      init_random_seq,            &
                                      random_uniform

   use sort_mod, only               : index_sort
   implicit none
!
! version controlled file description for error handling, do not edit
   character(len=*), parameter     :: source   = 'tes_co_cpsr_ascii_to_obs.f90'
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
   integer                         :: reject,k,l,kk,klev,kend,ilay
   integer                         :: i_min,j_min,reject_ak
   integer                         :: sum_reject,sum_accept,sum_total
   integer                         :: obs_accept,obs_o3_reten_freq,obs_co_reten_freq
   integer                         :: obs_co2_reten_freq,obs_ch4_reten_freq,obs_nh3_reten_freq
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
   real*8,dimension(num_qc)        :: tes_qc
   real*8,dimension(num_copies)    :: obs_val
!
   character*129                   :: filedir,filename,fileout
   character*129                   :: copy_meta_data
   character*129                   :: qc_meta_data='TES CO QC index'
   character*129                   :: chr_year,chr_month,chr_day
   character*129                   :: file_name='tes_co_cpsr_obs_seq'
   character*129                   :: data_type,cmd
   character*129                   :: path_model,file_model,file_in
!
   logical                         :: use_log_o3,use_log_co,use_log_co2,use_log_ch4,use_log_nh3
!
! Species-specific variables
   real                              :: co_trop_col_obs,co_trop_col_obs_prior,co_trop_col_obs_err
   real                              :: co_total_col_obs,co_total_col_obs_prior,co_total_col_obs_err
   real                              :: lat_obs,lon_obs,dofs_obs,trop_prs
   real*8                            :: lat_obs_r8,lon_obs_r8
   real,allocatable,dimension(:,:)   :: avgk_obs,cov_obs,cov_total
   real,allocatable,dimension(:)     :: avgk_diag_obs,co_obs,co_obs_prior,co_obs_err
   real*8,allocatable,dimension(:)   :: avgk_obs_r8
   real*8,allocatable,dimension(:,:) :: cov_obs_r8
   real*8,allocatable,dimension(:)   :: prior_obs_r8,co_obs_r8
   real,allocatable,dimension(:)     :: prs_obs
   real*8,allocatable,dimension(:)   :: prs_obs_r8
   real,allocatable,dimension(:)     :: prf_locl,prf_full
   real                              :: trop_sum,strat_sum
   real,allocatable,dimension(:,:)   :: lon,lat
   real,allocatable,dimension(:,:,:) :: prs_prt,prs_bas,prs_fld
   real,allocatable,dimension(:,:,:) :: tmp_prt,tmp_fld,vtmp_fld
   real,allocatable,dimension(:,:,:) :: co_fld,qmr_fld

   real                              :: trop_index,cov_trop,wt1_sum,wt2_sum
   real                              :: avgk_term, sum_co_obs
   real*8                            :: prior_total_r8
   real*8,allocatable,dimension(:)   :: avgk_total_obs_r8
!
! CPSR arrays
   real,allocatable,dimension(:)     :: co_shift,prior_shift
   real,allocatable,dimension(:)     :: co_cpsr,prior_cpsr
   real,allocatable,dimension(:,:)   :: avgk_shift,avgk_cpsr,cov_shift
   integer                           :: irow,irow_sh,icol,icol_sh
   integer                           :: nmodes,imds
!
   namelist /create_tes_obs_nml/filedir,filename,fileout, &
   bin_beg_sec,bin_end_sec,fac_obs_error,use_log_o3,use_log_co,use_log_co2,use_log_ch4, &
   use_log_nh3,lon_min,lon_max,lat_min,lat_max,path_model,file_model,nx_model,ny_model, &
   nz_model,obs_o3_reten_freq,obs_co_reten_freq,obs_co2_reten_freq,obs_ch4_reten_freq, &
   obs_nh3_reten_freq
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
   fac_err=1.0
!
! Record the current time, date, etc. to the logfile
   call initialize_utilities(source)
   call register_module(source,revision,revdate)
!
! Initialize the obs_sequence module
   call static_init_obs_sequence()
!
   call find_namelist_in_file("input.nml", "create_tes_obs_nml", iunit)
   read(iunit, nml = create_tes_obs_nml, iostat = io)
   call check_namelist_read(iunit, io, "create_tes_obs_nml")
!
! Record the namelist values used for the run ...
   call error_handler(E_MSG,'init_create_tes_obs','create_tes_obs_nml values are',' ',' ',' ')
   write(     *     , nml=create_tes_obs_nml)
!
! Initialize an obs_sequence structure
   call init_obs_sequence(seq, num_copies, num_qc, max_num_obs)
!
! Initialize the obs variable
   call init_obs(obs, num_copies, num_qc)
!
   do icopy =1, num_copies
      if (icopy == 1) then
         copy_meta_data='TES CO observation'
      else
         copy_meta_data='Truth'
      endif
      call set_copy_meta_data(seq, icopy, copy_meta_data)
   enddo
   call set_qc_meta_data(seq, 1, qc_meta_data)
!
!-------------------------------------------------------
! Read TES CO data
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
!   allocate(co_fld(nx_model,ny_model,nz_model))
!   file_in=trim(path_model)//'/'//trim(file_model)
!   call get_DART_diag_data(trim(file_in),'XLONG',lon,nx_model,ny_model,1,1)
!   call get_DART_diag_data(trim(file_in),'XLAT',lat,nx_model,ny_model,1,1)
!   call get_DART_diag_data(trim(file_in),'P',prs_prt,nx_model,ny_model,nz_model,1)
!   call get_DART_diag_data(trim(file_in),'PB',prs_bas,nx_model,ny_model,nz_model,1)
!   call get_DART_diag_data(trim(file_in),'T',tmp_prt,nx_model,ny_model,nz_model,1)
!   call get_DART_diag_data(trim(file_in),'QVAPOR',qmr_fld,nx_model,ny_model,nz_model,1)
!   call get_DART_diag_data(file_in,'co',co_fld,nx_model,ny_model,nz_model,1)
!   prs_fld(:,:,:)=prs_bas(:,:,:)+prs_prt(:,:,:)
!   tmp_fld(:,:,:)=300.+tmp_prt(:,:,:)
!   co_fld(:,:,:)=co_fld(:,:,:)*1.e-6
!
! Open TES CO binary file
   fileid=100
   write(6,*)'opening ',TRIM(filedir)//TRIM(filename)
   open(fileid,file=TRIM(filedir)//TRIM(filename),                     &
   form='formatted', status='old', iostat=ios)
!
! Read TES CO
   read(fileid,*,iostat=ios) data_type, obs_id, i_min, j_min
   do while (ios == 0)
      read(fileid,*,iostat=ios) yr_obs, mn_obs, &
      dy_obs, hh_obs, mm_obs, ss_obs
      read(fileid,*,iostat=ios) lat_obs,lon_obs
      if(lon_obs.lt.0.) lon_obs=lon_obs+360.
      read(fileid,*,iostat=ios) nlay_obs,nlev_obs
      allocate(prs_obs(nlay_obs))
      allocate(co_obs(nlay_obs))
      allocate(co_obs_prior(nlay_obs))
      allocate(co_obs_err(nlay_obs))
      allocate(avgk_total_obs_r8(nlay_obs))
      allocate(avgk_obs(nlay_obs,nlay_obs))
      allocate(avgk_diag_obs(nlay_obs))
      allocate(cov_obs(nlay_obs,nlay_obs))
      allocate(cov_total(nlay_obs,nlay_obs))
      allocate(prs_obs_r8(nlay_obs))
      allocate(prior_obs_r8(nlay_obs))
      allocate(avgk_obs_r8(nlay_obs))
      allocate(prf_locl(nlay_obs))
      allocate(prf_full(nlay_obs))
      read(fileid,*,iostat=ios) prs_obs(1:nlay_obs)
      read(fileid,*,iostat=ios) trop_prs
      read(fileid,*,iostat=ios) co_obs(1:nlay_obs)
      read(fileid,*,iostat=ios) co_obs_prior(1:nlay_obs)
      read(fileid,*,iostat=ios) co_obs_err(1:nlay_obs)
      read(fileid,*,iostat=ios) dofs_obs
      do k=1,nlay_obs
         read(fileid,*,iostat=ios) (avgk_obs(k,l),l=1,nlay_obs)
      enddo
      read(fileid,*,iostat=ios) avgk_diag_obs(1:nlay_obs)
      do k=1,nlay_obs
         read(fileid,*,iostat=ios) (cov_obs(k,l),l=1,nlay_obs)
      enddo
      do k=1,nlay_obs
         read(fileid,*,iostat=ios) (cov_total(k,l),l=1,nlay_obs)
      enddo
      read(fileid,*,iostat=ios) co_total_col_obs
      read(fileid,*,iostat=ios) co_total_col_obs_prior
      read(fileid,*,iostat=ios) co_total_col_obs_err

      print *, 'TES CO: Completed data read'
!
! Find first level above the ground
      klev=0
      do ilay=1,nlay_obs
         if(prs_obs(ilay).lt.0.) then
            klev=klev+1
         else
            exit
         endif
      enddo
      kend=nlay_obs-klev
!
      prs_obs(klev+1:nlay_obs)=prs_obs(klev+1:nlay_obs)*100.
      prs_obs_r8(1:kend)=prs_obs(klev+1:nlay_obs)
      prior_obs_r8(1:kend)=co_obs_prior(klev+1:nlay_obs)
      lon_obs_r8=lon_obs
      lat_obs_r8=lat_obs
!
! Obs thinning test
      obs_accept=obs_accept+1
      if(obs_accept/obs_co_reten_freq*obs_co_reten_freq.eq.obs_accept) then
!
! Check whether avgk row is missing.
         do ilay=klev+1,nlay_obs
            reject_ak=1
            do icol=klev+1,nlay_obs
               if(avgk_obs(ilay,icol).ne.-999 .or. cov_total(ilay,icol).ne.-999) then
                  reject_ak=0
                  exit
               endif  
            enddo
         enddo
         if(reject_ak.eq.1) cycle
!
! Apply CPSR transform
         allocate(co_shift(kend))
         allocate(prior_shift(kend))
         allocate(avgk_shift(kend,kend))
         allocate(cov_shift(kend,kend))
         allocate(co_cpsr(kend))
         allocate(prior_cpsr(kend))
         allocate(avgk_cpsr(kend,kend))
         do irow=klev+1,nlay_obs
            irow_sh=irow-klev
            co_shift(irow_sh)=log(co_obs(irow))
            prior_shift(irow_sh)=log(co_obs_prior(irow))
            do icol=klev+1,nlay_obs
               icol_sh=icol-klev
               avgk_shift(irow_sh,icol_sh)=avgk_obs(irow,icol)
               cov_shift(irow_sh,icol_sh)=fac_obs_error**2 * fac_err**2 * &
               cov_total(irow,icol)
            enddo
         enddo
!
! Calculate CPSRs for this retrieval profile
         call cpsr_calculation(nmodes,kend,co_cpsr,avgk_cpsr,prior_cpsr,co_shift,prior_shift, &
         avgk_shift,cov_shift,sum_accept)
!
! Loop through the dominant modes (OMI O3 is top to bottom)
         do imds=1,nmodes
            avgk_obs_r8(1:kend)=avgk_cpsr(imds,1:kend)
!
! Process accepted observations
            sum_accept=sum_accept+1
            qc_count=qc_count+1
!
! Obs is a CPSR (the compressed and rotated obs errors are unity)
            obs_val(:)=co_cpsr(imds)
            obs_err_var=1.
            tes_qc(:)=0
            obs_time=set_date(yr_obs,mn_obs,dy_obs,hh_obs,mm_obs,ss_obs)
            call get_time(obs_time, seconds, days)
!
            which_vert=-2      ! undefined
!            which_vert=-1      ! surface
!            which_vert=1       ! level
!            which_vert=2       ! pressure
!
            obs_kind = TES_CO_CPSR
! (0 <= lon_obs <= 360); (-90 <= lat_obs <= 90)
            level=real(imds)
            obs_location=set_location(lon_obs_r8, lat_obs_r8, level, which_vert)
!
            call set_obs_def_type_of_obs(obs_def, obs_kind)
            call set_obs_def_location(obs_def, obs_location)
            call set_obs_def_time(obs_def, obs_time)
            call set_obs_def_error_variance(obs_def, obs_err_var)
            call set_obs_def_tes_co_cpsr(qc_count, prs_obs_r8(1:kend), &
            avgk_obs_r8(1:kend), prior_obs_r8(1:kend), kend, nmodes, nlay_obs)
            call set_obs_def_key(obs_def, qc_count)
            call set_obs_values(obs, obs_val, 1)
            call set_qc(obs, tes_qc, num_qc)
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
         deallocate(co_shift)
         deallocate(prior_shift)
         deallocate(avgk_shift)
         deallocate(cov_shift)
         deallocate(co_cpsr)
         deallocate(prior_cpsr)
         deallocate(avgk_cpsr)
      endif
      deallocate(prs_obs) 
      deallocate(co_obs) 
      deallocate(co_obs_prior)
      deallocate(co_obs_err)
      deallocate(avgk_total_obs_r8)
      deallocate(avgk_obs)
      deallocate(avgk_diag_obs)
      deallocate(cov_obs)
      deallocate(cov_total)
      deallocate(prs_obs_r8) 
      deallocate(avgk_obs_r8)
      deallocate(prior_obs_r8)
      deallocate(prf_locl) 
      deallocate(prf_full)
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
!   deallocate(co_fld)
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
end program tes_co_cpsr_ascii_to_obs
