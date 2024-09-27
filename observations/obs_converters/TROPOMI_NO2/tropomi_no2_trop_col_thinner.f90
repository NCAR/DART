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
program tropomi_no2_thinner
   use    apm_stats_utilities, only : median_code,                &
                                      mode_code,                  &
                                      vertical_column
   
   use    utilities_mod, only       : register_module,            &
                                      initialize_utilities,       &
                                      find_namelist_in_file,      &
                                      check_namelist_read,        &
                                      open_file,                  &
                                      close_file                 
   implicit none
!
! Version controlled file description for error handling, do not edit
   character(len=*), parameter     :: source   = 'tropomi_no2_thinner.f90'
   character(len=*), parameter     :: revision = ''
   character(len=*), parameter     :: revdate  = ''
!
   integer,parameter               :: max_num_obs=200
   integer                         :: iunit,io,fileid,ios
   integer                         :: icnt
   integer                         :: i,j
!
! TROPOMI NO2 variable declarations
   integer                         :: obs_id,i_min,j_min
   integer                         :: yr_obs,mn_obs,dy_obs,hh_obs,mm_obs,ss_obs
   real                            :: lat_obs,lon_obs
   integer                         :: nlay_obs,nlev_obs,ilv
   real                            :: amfstrat,amfstrat_clr,amfstrat_cld
   real                            :: amftrop,amftrop_clr,amftrop_cld
   real                            :: cld_frac,cld_prs,cld_rad_frac
   real                            :: col_amt,col_amt_err
   real                            :: col_amt_trop,col_amt_trop_err
   real                            :: slnt_col_amt,slnt_col_amt_err
   real                            :: prs_trop,zenang
   real,allocatable,dimension(:)   :: scat_wt,prs_obs
!
! Namelist variable declarations
   character(len=129)              :: filedir,filename,fileout
   integer                         :: year,month,day,hour
   real                            :: bin_beg_sec,bin_end_sec
   real                            :: fac_obs_error
   logical                         :: use_log_o3,use_log_no2,use_log_so2, &
                                      use_log_hcho,use_log_co,use_log_ch4
   real                            :: lon_min,lon_max,lat_min,lat_max
   character(len=129)              :: path_model,file_model,data_type
   integer                         :: nx_model,ny_model,nz_model
   integer                         :: obs_o3_reten_freq,obs_no2_reten_freq, &
                                      obs_so2_reten_freq,obs_hcho_reten_freq, &
                                      obs_co_reten_freq,obs_ch4_reten_freq
!  
   type observation_data
      integer                          :: icnt
      real                             :: mean_col
      real                             :: mean_err
      real                             :: stdv_col
      real                             :: stdv_err
      real                             :: median_col      
      real                             :: mode_col      
      character(len=129)               :: data_type
      integer                          :: obs_id
      integer,allocatable,dimension(:) :: yr_obs
      integer,allocatable,dimension(:) :: mn_obs
      integer,allocatable,dimension(:) :: dy_obs
      integer,allocatable,dimension(:) :: hh_obs
      integer,allocatable,dimension(:) :: mm_obs
      integer,allocatable,dimension(:) :: ss_obs
      real,allocatable,dimension(:)    :: lat_obs
      real,allocatable,dimension(:)    :: lon_obs
      integer,allocatable,dimension(:) :: nlay_obs
      integer,allocatable,dimension(:) :: nlev_obs
      real,allocatable,dimension(:,:)  :: prs_obs
      real,allocatable,dimension(:,:)  :: avgk_obs
      real,allocatable,dimension(:)    :: no2_trop_col
      real,allocatable,dimension(:)    :: no2_trop_col_err
      real,allocatable,dimension(:)    :: amf_trop_col
      integer,allocatable,dimension(:) :: trop_indx
   end type observation_data
!
   type(observation_data),allocatable,dimension(:,:) :: tropomi_data
   real                                :: dist
   real,allocatable,dimension(:,:)     :: dist_min
   real,allocatable,dimension(:,:)     :: icnt_min
!
   namelist /create_tropomi_obs_nml/filedir,filename,fileout, &
   bin_beg_sec,bin_end_sec,fac_obs_error,use_log_co,use_log_o3,use_log_no2,use_log_so2, &
   use_log_ch4,use_log_hcho,lon_min,lon_max,lat_min,lat_max,path_model,file_model,nx_model, &
   ny_model,nz_model,obs_co_reten_freq,obs_o3_reten_freq,obs_no2_reten_freq,obs_so2_reten_freq, &
   obs_ch4_reten_freq,obs_hcho_reten_freq
!
! Record the current time, date, etc. to the logfile
   call initialize_utilities(source)
   call register_module(source,revision,revdate)
!
   call find_namelist_in_file("input.nml", "create_tropomi_obs_nml", iunit)
   read(iunit, nml = create_tropomi_obs_nml, iostat = io)
   call check_namelist_read(iunit, io, "create_tropomi_obs_nml")
!
! Read TROPOMI NO2
   allocate(tropomi_data(nx_model,ny_model))
   allocate(dist_min(nx_model,ny_model))
   allocate(icnt_min(nx_model,ny_model))
!
   tropomi_data(:,:)%icnt=0
   fileid=100
   write(6,*)'opening ',TRIM(TRIM(filedir)//TRIM(filename))
   open(unit=fileid,file=TRIM(TRIM(filedir)//TRIM(filename)), &
   form='formatted', status='old', iostat=ios)
   read(fileid,*,iostat=ios) data_type, obs_id, &
   i_min, j_min
   tropomi_data(i_min,j_min)%data_type=data_type
   tropomi_data(i_min,j_min)%obs_id=obs_id
   do while (ios == 0)
      tropomi_data(i_min,j_min)%icnt=tropomi_data(i_min,j_min)%icnt+1
      icnt=tropomi_data(i_min,j_min)%icnt
      if(tropomi_data(i_min,j_min)%icnt.eq.1) then
         allocate(tropomi_data(i_min,j_min)%yr_obs(max_num_obs))
         allocate(tropomi_data(i_min,j_min)%mn_obs(max_num_obs))
         allocate(tropomi_data(i_min,j_min)%dy_obs(max_num_obs))
         allocate(tropomi_data(i_min,j_min)%hh_obs(max_num_obs))
         allocate(tropomi_data(i_min,j_min)%mm_obs(max_num_obs))
         allocate(tropomi_data(i_min,j_min)%ss_obs(max_num_obs))
         allocate(tropomi_data(i_min,j_min)%lat_obs(max_num_obs))
         allocate(tropomi_data(i_min,j_min)%lon_obs(max_num_obs))
         allocate(tropomi_data(i_min,j_min)%nlay_obs(max_num_obs))
         allocate(tropomi_data(i_min,j_min)%nlev_obs(max_num_obs))
         allocate(tropomi_data(i_min,j_min)%no2_trop_col(max_num_obs))
         allocate(tropomi_data(i_min,j_min)%no2_trop_col_err(max_num_obs))
      endif
!
      read(fileid,*,iostat=ios) &
         tropomi_data(i_min,j_min)%yr_obs(icnt), &
         tropomi_data(i_min,j_min)%mn_obs(icnt), &
         tropomi_data(i_min,j_min)%dy_obs(icnt), &
         tropomi_data(i_min,j_min)%hh_obs(icnt), &
         tropomi_data(i_min,j_min)%mm_obs(icnt), &
         tropomi_data(i_min,j_min)%ss_obs(icnt)
      read(fileid,*,iostat=ios) &
         tropomi_data(i_min,j_min)%lat_obs(icnt), &
         tropomi_data(i_min,j_min)%lon_obs(icnt)
      read(fileid,*,iostat=ios) &
         tropomi_data(i_min,j_min)%nlay_obs(icnt), &
         tropomi_data(i_min,j_min)%nlev_obs(icnt)
         nlay_obs=tropomi_data(i_min,j_min)%nlay_obs(icnt)
         nlev_obs=tropomi_data(i_min,j_min)%nlev_obs(icnt)
      if(tropomi_data(i_min,j_min)%icnt.eq.1) then
         allocate(tropomi_data(i_min,j_min)%prs_obs(max_num_obs,nlay_obs))
         allocate(tropomi_data(i_min,j_min)%avgk_obs(max_num_obs,nlay_obs))
      endif
      read(fileid,*,iostat=ios) &
         tropomi_data(i_min,j_min)%prs_obs(icnt,1:nlay_obs)
      read(fileid,*,iostat=ios) &
         tropomi_data(i_min,j_min)%avgk_obs(icnt,1:nlay_obs)
      read(fileid,*,iostat=ios) &
         tropomi_data(i_min,j_min)%no2_trop_col(icnt), &
         tropomi_data(i_min,j_min)%no2_trop_col_err(icnt)
      read(fileid,*,iostat=ios) &
         tropomi_data(i_min,j_min)%amf_trop_col(icnt)
      read(fileid,*,iostat=ios) &
         tropomi_data(i_min,j_min)%trop_indx(icnt)      
      read(fileid,*,iostat=ios) data_type, obs_id, &
         i_min, j_min
      tropomi_data(i_min,j_min)%data_type=data_type
      tropomi_data(i_min,j_min)%obs_id=obs_id
   enddo
!
! Find the median
   tropomi_data(:,:)%median_col=0.   
   do i=1,nx_model
      do j=1,ny_model
         if(tropomi_data(i,j)%icnt.gt.1) then
            call median_code(tropomi_data(i,j)%median_col,tropomi_data(i,j)%no2_trop_col(1:tropomi_data(i,j)%icnt), &
            tropomi_data(i,j)%icnt)
         endif
      enddo
   enddo
!
! Find the mode
   tropomi_data(:,:)%mode_col=0.   
   do i=1,nx_model
      do j=1,ny_model
         if(tropomi_data(i,j)%icnt.gt.20) then
            call mode_code(tropomi_data(i,j)%mode_col,tropomi_data(i,j)%no2_trop_col(1:tropomi_data(i,j)%icnt), &
            tropomi_data(i,j)%icnt)
         endif
      enddo
   enddo
!
! Calculate mean
   do i=1,nx_model
      do j=1,ny_model
         if(tropomi_data(i,j)%icnt.gt.0) then
            tropomi_data(i,j)%mean_col=0.
            tropomi_data(i,j)%mean_err=0.
            do icnt=1,tropomi_data(i,j)%icnt
               tropomi_data(i,j)%mean_col=tropomi_data(i,j)%mean_col+tropomi_data(i,j)%no2_trop_col(icnt)
               tropomi_data(i,j)%mean_err=tropomi_data(i,j)%mean_err+tropomi_data(i,j)%no2_trop_col_err(icnt)
            enddo
            tropomi_data(i,j)%mean_col=tropomi_data(i,j)%mean_col/float(tropomi_data(i,j)%icnt)
            tropomi_data(i,j)%mean_err=tropomi_data(i,j)%mean_err/float(tropomi_data(i,j)%icnt)
         endif
      enddo
   enddo
!
! Calculate standard deviation   
   do i=1,nx_model
      do j=1,ny_model
         if(tropomi_data(i,j)%icnt.gt.1) then
            tropomi_data(i,j)%stdv_col=0.
            tropomi_data(i,j)%stdv_err=0.
            do icnt=1,tropomi_data(i,j)%icnt
               tropomi_data(i,j)%stdv_col=tropomi_data(i,j)%stdv_col+(tropomi_data(i,j)%no2_trop_col(icnt)- &
               tropomi_data(i,j)%mean_col)**2.
               tropomi_data(i,j)%stdv_err=tropomi_data(i,j)%stdv_err+(tropomi_data(i,j)%no2_trop_col_err(icnt)- &
               tropomi_data(i,j)%mean_err)**2.
            enddo
            tropomi_data(i,j)%stdv_col=sqrt(tropomi_data(i,j)%stdv_col/float(tropomi_data(i,j)%icnt-1))
            tropomi_data(i,j)%stdv_err=sqrt(tropomi_data(i,j)%stdv_err/float(tropomi_data(i,j)%icnt-1))
         endif
      enddo
   enddo
!
! Find the mean's closest neighbor based on trop_col
   dist_min(:,:)=1.e20
   icnt_min(:,:)=0
   rewind(fileid)
   do i=1,nx_model
      do j=1,ny_model
         dist=0
         if(tropomi_data(i,j)%icnt.eq.1) then
            icnt_min(:,:)=1
         else if(tropomi_data(i,j)%icnt.gt.1) then
            do icnt=1,tropomi_data(i,j)%icnt
               dist=abs(tropomi_data(i,j)%mean_col-tropomi_data(i,j)%no2_trop_col(icnt))
               if(dist.lt.dist_min(i,j)) then
                  dist_min(i,j)=dist
                  icnt_min(i,j)=icnt
               endif
            enddo
         endif
!
! Write data to file
         if(tropomi_data(i,j)%icnt.gt.0) then
            icnt=icnt_min(i,j)
            write(fileid,*,iostat=ios) &
               trim(tropomi_data(i,j)%data_type), &
               tropomi_data(i,j)%obs_id, i, j
            write(fileid,*,iostat=ios) &
               tropomi_data(i,j)%yr_obs(icnt), &
               tropomi_data(i,j)%mn_obs(icnt), &
               tropomi_data(i,j)%dy_obs(icnt), &
               tropomi_data(i,j)%hh_obs(icnt), &
               tropomi_data(i,j)%mm_obs(icnt), &
               tropomi_data(i,j)%ss_obs(icnt)
            write(fileid,*,iostat=ios) &
               tropomi_data(i,j)%lat_obs(icnt), &
               tropomi_data(i,j)%lon_obs(icnt)
            write(fileid,*,iostat=ios) &
               tropomi_data(i,j)%nlay_obs(icnt), &
               tropomi_data(i,j)%nlev_obs(icnt)
            write(fileid,*,iostat=ios) &
               tropomi_data(i,j)%prs_obs(icnt,1:nlay_obs)
            write(fileid,*,iostat=ios) &
               tropomi_data(i,j)%avgk_obs(icnt,1:nlay_obs)
            write(fileid,*,iostat=ios) &
               tropomi_data(i,j)%no2_trop_col(icnt), &
               tropomi_data(i,j)%no2_trop_col_err(icnt)
            write(fileid,*,iostat=ios) &
               tropomi_data(i,j)%amf_trop_col(icnt)
            write(fileid,*,iostat=ios) &
               tropomi_data(i,j)%trop_indx(icnt)
         endif
      enddo
   enddo
!               
! Deallocate structure arrays
   do i=1,nx_model
      do j=1,ny_model
         if(tropomi_data(i,j)%icnt.gt.0) then
            deallocate(tropomi_data(i,j)%yr_obs)
            deallocate(tropomi_data(i,j)%mn_obs)
            deallocate(tropomi_data(i,j)%dy_obs)
            deallocate(tropomi_data(i,j)%hh_obs)
            deallocate(tropomi_data(i,j)%mm_obs)
            deallocate(tropomi_data(i,j)%ss_obs)
            deallocate(tropomi_data(i,j)%lat_obs)
            deallocate(tropomi_data(i,j)%lon_obs)
            deallocate(tropomi_data(i,j)%nlay_obs)
            deallocate(tropomi_data(i,j)%nlev_obs)
            deallocate(tropomi_data(i,j)%prs_obs)
            deallocate(tropomi_data(i,j)%avgk_obs)
            deallocate(tropomi_data(i,j)%no2_trop_col)
            deallocate(tropomi_data(i,j)%no2_trop_col_err)
            deallocate(tropomi_data(i,j)%amf_trop_col)
            deallocate(tropomi_data(i,j)%trop_indx)
         endif
      enddo
   enddo
   deallocate (tropomi_data)
   deallocate (dist_min)
   deallocate (icnt_min)   
end program tropomi_no2_thinner
