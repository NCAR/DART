!
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
program modis_aod_thinner_matlab
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
   character(len=*), parameter     :: source   = 'modis_aod_thinning.f90'
   character(len=*), parameter     :: revision = ''
   character(len=*), parameter     :: revdate  = ''
!
   integer,parameter               :: max_num_obs=200
   integer                         :: iunit,io,fileid,ios
   integer                         :: icnt
   integer                         :: i,j
!
! MODIS AOD variable declarations
   integer                         :: obs_id,i_min,j_min
   integer                         :: yr_obs,mn_obs,dy_obs,hh_obs,mm_obs,ss_obs
   real                            :: lat_obs,lon_obs
   integer                         :: nlay_obs,nlev_obs,ilv
   real                            :: taudb, taudb_best, taudb_err
   real                            :: prs_trop,zenang
   real,allocatable,dimension(:)   :: scat_wt,prs_obs
!
! Namelist variable declarations
   character(len=129)              :: filedir,filename,fileout
   integer                         :: year,month,day,hour
   real                            :: bin_beg_sec,bin_end_sec
   real                            :: fac_obs_error
   logical                         :: use_log_aod
   real                            :: lon_min,lon_max,lat_min,lat_max
   character(len=129)              :: path_model,file_model,data_type
   integer                         :: nx_model,ny_model,nz_model
   integer                         :: obs_aod_reten_freq
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
      real,allocatable,dimension(:)    :: taudb
      real,allocatable,dimension(:)    :: taudb_best
      real,allocatable,dimension(:)    :: taudb_err
   end type observation_data
!
   type(observation_data),allocatable,dimension(:,:) :: modis_data
   real                                :: dist
   real,allocatable,dimension(:,:)     :: dist_min
   real,allocatable,dimension(:,:)     :: icnt_min
!
   namelist /create_modis_obs_nml/filedir,filename,fileout,year,month,day,hour, &
   bin_beg_sec,bin_end_sec,fac_obs_error,use_log_aod, &
   lon_min,lon_max,lat_min,lat_max,path_model,file_model,nx_model, &
   ny_model,nz_model,obs_aod_reten_freq
!
! Record the current time, date, etc. to the logfile
   call initialize_utilities(source)
   call register_module(source,revision,revdate)
!
   call find_namelist_in_file("input.nml", "create_modis_obs_nml", iunit)
   read(iunit, nml = create_modis_obs_nml, iostat = io)
   call check_namelist_read(iunit, io, "create_modis_obs_nml")
!
! Read MODIS AOD
   allocate(modis_data(nx_model,ny_model))
   allocate(dist_min(nx_model,ny_model))
   allocate(icnt_min(nx_model,ny_model))
!
   modis_data(:,:)%icnt=0
   fileid=100
   write(6,*)'opening ',TRIM(TRIM(filedir)//TRIM(filename))
   open(unit=fileid,file=TRIM(TRIM(filedir)//TRIM(filename)), &
   form='formatted', status='old', iostat=ios)
   read(fileid,*,iostat=ios) data_type, obs_id, &
   i_min, j_min
   modis_data(i_min,j_min)%data_type=data_type
   modis_data(i_min,j_min)%obs_id=obs_id
   do while (ios == 0)
      modis_data(i_min,j_min)%icnt=modis_data(i_min,j_min)%icnt+1
      icnt=modis_data(i_min,j_min)%icnt
      if(modis_data(i_min,j_min)%icnt.eq.1) then
         allocate(modis_data(i_min,j_min)%yr_obs(max_num_obs))
         allocate(modis_data(i_min,j_min)%mn_obs(max_num_obs))
         allocate(modis_data(i_min,j_min)%dy_obs(max_num_obs))
         allocate(modis_data(i_min,j_min)%hh_obs(max_num_obs))
         allocate(modis_data(i_min,j_min)%mm_obs(max_num_obs))
         allocate(modis_data(i_min,j_min)%ss_obs(max_num_obs))
         allocate(modis_data(i_min,j_min)%lat_obs(max_num_obs))
         allocate(modis_data(i_min,j_min)%lon_obs(max_num_obs))
         allocate(modis_data(i_min,j_min)%taudb(max_num_obs))
         allocate(modis_data(i_min,j_min)%taudb_best(max_num_obs))
         allocate(modis_data(i_min,j_min)%taudb_err(max_num_obs))

      endif
!
      read(fileid,*,iostat=ios) &
         modis_data(i_min,j_min)%yr_obs(icnt), &
         modis_data(i_min,j_min)%mn_obs(icnt), &
         modis_data(i_min,j_min)%dy_obs(icnt), &
         modis_data(i_min,j_min)%hh_obs(icnt), &
         modis_data(i_min,j_min)%mm_obs(icnt), &
         modis_data(i_min,j_min)%ss_obs(icnt)
      read(fileid,*,iostat=ios) &
         modis_data(i_min,j_min)%lat_obs(icnt), &
         modis_data(i_min,j_min)%lon_obs(icnt)
      read(fileid,*,iostat=ios) &
         modis_data(i_min,j_min)%taudb(icnt), &
         modis_data(i_min,j_min)%taudb_best(icnt), &
         modis_data(i_min,j_min)%taudb_err(icnt)
      read(fileid,*,iostat=ios) data_type, obs_id, &
         i_min, j_min
      modis_data(i_min,j_min)%data_type=data_type
      modis_data(i_min,j_min)%obs_id=obs_id
   enddo
!
! Find the median
   modis_data(:,:)%median_col=0.   
   do i=1,nx_model
      do j=1,ny_model
         if(modis_data(i,j)%icnt.gt.1) then
            call median_code(modis_data(i,j)%median_col,modis_data(i,j)%taudb(1:modis_data(i,j)%icnt), &
            modis_data(i,j)%icnt)
         endif
      enddo
   enddo
!
! Find the mode
   modis_data(:,:)%mode_col=0.   
   do i=1,nx_model
      do j=1,ny_model
         if(modis_data(i,j)%icnt.gt.20) then
            call mode_code(modis_data(i,j)%mode_col,modis_data(i,j)%taudb(1:modis_data(i,j)%icnt), &
            modis_data(i,j)%icnt)
         endif
      enddo
   enddo
!
! Calculate mean
   do i=1,nx_model
      do j=1,ny_model
         if(modis_data(i,j)%icnt.gt.0) then
            modis_data(i,j)%mean_col=0.
            modis_data(i,j)%mean_err=0.
            do icnt=1,modis_data(i,j)%icnt
               modis_data(i,j)%mean_col=modis_data(i,j)%mean_col+modis_data(i,j)%taudb(icnt)
               modis_data(i,j)%mean_err=modis_data(i,j)%mean_err+modis_data(i,j)%taudb_err(icnt)
            enddo
            modis_data(i,j)%mean_col=modis_data(i,j)%mean_col/float(modis_data(i,j)%icnt)
            modis_data(i,j)%mean_err=modis_data(i,j)%mean_err/float(modis_data(i,j)%icnt)
         endif
      enddo
   enddo
!
! Calculate standard deviation   
   do i=1,nx_model
      do j=1,ny_model
         if(modis_data(i,j)%icnt.gt.1) then
            modis_data(i,j)%stdv_col=0.
            modis_data(i,j)%stdv_err=0.
            do icnt=1,modis_data(i,j)%icnt
               modis_data(i,j)%stdv_col=modis_data(i,j)%stdv_col+(modis_data(i,j)%taudb(icnt)- &
               modis_data(i,j)%mean_col)**2.
               modis_data(i,j)%stdv_err=modis_data(i,j)%stdv_err+(modis_data(i,j)%taudb_err(icnt)- &
               modis_data(i,j)%mean_err)**2.
            enddo
            modis_data(i,j)%stdv_col=sqrt(modis_data(i,j)%stdv_col/float(modis_data(i,j)%icnt-1))
            modis_data(i,j)%stdv_err=sqrt(modis_data(i,j)%stdv_err/float(modis_data(i,j)%icnt-1))
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
         if(modis_data(i,j)%icnt.eq.1) then
            icnt_min(:,:)=1
         else if(modis_data(i,j)%icnt.gt.1) then
            do icnt=1,modis_data(i,j)%icnt
               dist=abs(modis_data(i,j)%mean_col-modis_data(i,j)%taudb(icnt))
               if(dist.lt.dist_min(i,j)) then
                  dist_min(i,j)=dist
                  icnt_min(i,j)=icnt
               endif
            enddo
         endif
!
! Write data to file
         if(modis_data(i,j)%icnt.gt.0) then
            icnt=icnt_min(i,j)
            write(fileid,*,iostat=ios) &
               trim(modis_data(i,j)%data_type), &
               modis_data(i,j)%obs_id, i, j
            write(fileid,*,iostat=ios) &
               modis_data(i,j)%yr_obs(icnt), &
               modis_data(i,j)%mn_obs(icnt), &
               modis_data(i,j)%dy_obs(icnt), &
               modis_data(i,j)%hh_obs(icnt), &
               modis_data(i,j)%mm_obs(icnt), &
               modis_data(i,j)%ss_obs(icnt)
            write(fileid,*,iostat=ios) &
               modis_data(i,j)%lat_obs(icnt), &
               modis_data(i,j)%lon_obs(icnt)
            write(fileid,*,iostat=ios) &
               modis_data(i,j)%taudb(icnt), &
               modis_data(i,j)%taudb_best(icnt), &
               modis_data(i,j)%taudb_err(icnt)
         endif
      enddo
   enddo
!               
! Deallocate structure arrays
   do i=1,nx_model
      do j=1,ny_model
         if(modis_data(i,j)%icnt.gt.0) then
            deallocate(modis_data(i,j)%yr_obs)
            deallocate(modis_data(i,j)%mn_obs)
            deallocate(modis_data(i,j)%dy_obs)
            deallocate(modis_data(i,j)%hh_obs)
            deallocate(modis_data(i,j)%mm_obs)
            deallocate(modis_data(i,j)%ss_obs)
            deallocate(modis_data(i,j)%lat_obs)
            deallocate(modis_data(i,j)%lon_obs)
            deallocate(modis_data(i,j)%taudb)
            deallocate(modis_data(i,j)%taudb_best)
            deallocate(modis_data(i,j)%taudb_err)
            endif
      enddo
   enddo
   deallocate (modis_data)
   deallocate (dist_min)
   deallocate (icnt_min)   
end program modis_aod_thinner_matlab
