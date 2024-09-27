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
program mls_o3_profile_thinner
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
   character(len=*), parameter      :: source   = 'mls_o3_thinning.f90'
   character(len=*), parameter      :: revision = ''
   character(len=*), parameter      :: revdate  = ''
!
   integer,parameter                :: max_num_obs=200
   integer                          :: iunit,io,fileid,ios
   integer                          :: icnt
   integer                          :: i,j,k
!
! MLS O3 variable declarations
   integer                           :: obs_id,i_min,j_min
   integer                           :: yr_obs,mn_obs,dy_obs,hh_obs,mm_obs,ss_obs
   real                              :: lat_obs,lon_obs
   integer                           :: nlay_obs,nlev_obs,ilv
   integer                           :: ret_nlay,tru_nlay
   real                              :: dofs
   real,allocatable,dimension(:,:)   :: prs_obs
   real,allocatable,dimension(:,:)   :: ret_lay_obs
   real,allocatable,dimension(:,:)   :: tru_lay_obs
   real,allocatable,dimension(:,:)   :: o3_obs
   real,allocatable,dimension(:,:)   :: o3_obs_err
   real,allocatable,dimension(:,:)   :: prior_obs
   real,allocatable,dimension(:)     :: prior_err
   real,allocatable,dimension(:,:,:) :: avgk_obs
   real,allocatable,dimension(:)     :: o3_col_obs
   real,allocatable,dimension(:)     :: o3_col_obs_err
!
! Namelist variable declarations
   character(len=129)               :: filedir,filename,fileout
   integer                          :: year,month,day,hour
   real                             :: bin_beg_sec,bin_end_sec
   real                             :: fac_obs_error
   logical                          :: use_log_o3,use_log_hno3
   real                             :: lon_min,lon_max,lat_min,lat_max
   character(len=129)               :: path_model,file_model,data_type
   integer                          :: nx_model,ny_model,nz_model
   integer                          :: obs_o3_reten_freq,obs_hno3_reten_freq
!  
   type observation_data
      integer                           :: icnt
      real                              :: mean_col
      real                              :: stdv_col
      real                              :: median_col      
      real                              :: mode_col      
      real,allocatable,dimension(:)     :: vert_col
      character(len=129)                :: data_type
      integer                           :: obs_id
      integer,allocatable,dimension(:)  :: yr_obs
      integer,allocatable,dimension(:)  :: mn_obs
      integer,allocatable,dimension(:)  :: dy_obs
      integer,allocatable,dimension(:)  :: hh_obs
      integer,allocatable,dimension(:)  :: mm_obs
      integer,allocatable,dimension(:)  :: ss_obs
      real,allocatable,dimension(:)     :: lat_obs
      real,allocatable,dimension(:)     :: lon_obs
      integer,allocatable,dimension(:)  :: nlay_obs
      integer,allocatable,dimension(:)  :: nlev_obs
      integer,allocatable,dimension(:)  :: ret_nlay
      integer,allocatable,dimension(:)  :: tru_nlay
      real,allocatable,dimension(:)     :: dofs
      real,allocatable,dimension(:,:)   :: prs_obs
      real,allocatable,dimension(:,:)   :: ret_lay_obs
      real,allocatable,dimension(:,:)   :: tru_lay_obs
      real,allocatable,dimension(:,:)   :: o3_obs
      real,allocatable,dimension(:,:)   :: o3_obs_err
      real,allocatable,dimension(:,:)   :: prior_obs
      real,allocatable,dimension(:)     :: prior_err
      real,allocatable,dimension(:,:,:) :: avgk_obs
      real,allocatable,dimension(:)     :: o3_col_obs
      real,allocatable,dimension(:)     :: o3_col_obs_err
   end type observation_data
!
   type(observation_data),allocatable,dimension(:,:) :: mls_data
   real                                :: dist
   real,allocatable,dimension(:,:)     :: dist_min
   real,allocatable,dimension(:,:)     :: icnt_min
!
   namelist /create_mls_obs_nml/filedir,filename,fileout,year,month,day,hour, &
   bin_beg_sec,bin_end_sec,fac_obs_error,use_log_o3,use_log_hno3, &
   lon_min,lon_max,lat_min,lat_max,path_model,file_model,nx_model, &
   ny_model,nz_model,obs_o3_reten_freq,obs_hno3_reten_freq
!
! Record the current time, date, etc. to the logfile
   call initialize_utilities(source)
   call register_module(source,revision,revdate)
!
   call find_namelist_in_file("input.nml", "create_mls_obs_nml", iunit)
   read(iunit, nml = create_mls_obs_nml, iostat = io)
   call check_namelist_read(iunit, io, "create_mls_obs_nml")
!
! Read MLS O3
   allocate(mls_data(nx_model,ny_model))
   allocate(dist_min(nx_model,ny_model))
   allocate(icnt_min(nx_model,ny_model))
!
   mls_data(:,:)%icnt=0
   fileid=100
   write(6,*)'opening ',TRIM(TRIM(filedir)//TRIM(filename))
   open(unit=fileid,file=TRIM(TRIM(filedir)//TRIM(filename)), &
   form='formatted', status='old', iostat=ios)
   read(fileid,*,iostat=ios) data_type, obs_id, &
   i_min, j_min
   mls_data(i_min,j_min)%data_type=data_type
   mls_data(i_min,j_min)%obs_id=obs_id
   do while (ios == 0)
      mls_data(i_min,j_min)%icnt=mls_data(i_min,j_min)%icnt+1
      icnt=mls_data(i_min,j_min)%icnt
      if(mls_data(i_min,j_min)%icnt.eq.1) then
         allocate(mls_data(i_min,j_min)%vert_col(max_num_obs))
         allocate(mls_data(i_min,j_min)%yr_obs(max_num_obs))
         allocate(mls_data(i_min,j_min)%mn_obs(max_num_obs))
         allocate(mls_data(i_min,j_min)%dy_obs(max_num_obs))
         allocate(mls_data(i_min,j_min)%hh_obs(max_num_obs))
         allocate(mls_data(i_min,j_min)%mm_obs(max_num_obs))
         allocate(mls_data(i_min,j_min)%ss_obs(max_num_obs))
         allocate(mls_data(i_min,j_min)%lat_obs(max_num_obs))
         allocate(mls_data(i_min,j_min)%lon_obs(max_num_obs))
         allocate(mls_data(i_min,j_min)%nlay_obs(max_num_obs))
         allocate(mls_data(i_min,j_min)%nlev_obs(max_num_obs))
         allocate(mls_data(i_min,j_min)%ret_nlay(max_num_obs))
         allocate(mls_data(i_min,j_min)%tru_nlay(max_num_obs))
      endif
!
      read(fileid,*,iostat=ios) &
         mls_data(i_min,j_min)%yr_obs(icnt), &
         mls_data(i_min,j_min)%mn_obs(icnt), &
         mls_data(i_min,j_min)%dy_obs(icnt), &
         mls_data(i_min,j_min)%hh_obs(icnt), &
         mls_data(i_min,j_min)%mm_obs(icnt), &
         mls_data(i_min,j_min)%ss_obs(icnt)      
      read(fileid,*,iostat=ios) &
         mls_data(i_min,j_min)%lat_obs(icnt), &
         mls_data(i_min,j_min)%lon_obs(icnt)
      read(fileid,*,iostat=ios) &
         mls_data(i_min,j_min)%nlay_obs(icnt), &
         mls_data(i_min,j_min)%nlev_obs(icnt)
      read(fileid,*,iostat=ios) &
         mls_data(i_min,j_min)%ret_nlay(icnt), &
         mls_data(i_min,j_min)%tru_nlay(icnt)
         nlay_obs=mls_data(i_min,j_min)%nlay_obs(icnt)
         nlev_obs=mls_data(i_min,j_min)%nlev_obs(icnt)
         ret_nlay=mls_data(i_min,j_min)%ret_nlay(icnt)
         tru_nlay=mls_data(i_min,j_min)%tru_nlay(icnt)
      if(mls_data(i_min,j_min)%icnt.eq.1) then
         allocate(mls_data(i_min,j_min)%prs_obs(max_num_obs,nlay_obs))
         allocate(mls_data(i_min,j_min)%ret_lay_obs(max_num_obs,ret_nlay))
         allocate(mls_data(i_min,j_min)%tru_lay_obs(max_num_obs,tru_nlay))
         allocate(mls_data(i_min,j_min)%o3_obs(max_num_obs,nlay_obs))
         allocate(mls_data(i_min,j_min)%o3_obs_err(max_num_obs,nlay_obs))
         allocate(mls_data(i_min,j_min)%prior_obs(max_num_obs,nlay_obs))
         allocate(mls_data(i_min,j_min)%prior_err(max_num_obs))
         allocate(mls_data(i_min,j_min)%avgk_obs(max_num_obs,ret_nlay,tru_nlay))
         allocate(mls_data(i_min,j_min)%o3_col_obs(max_num_obs))
         allocate(mls_data(i_min,j_min)%o3_col_obs_err(max_num_obs))
      endif
      read(fileid,*,iostat=ios) &
         mls_data(i_min,j_min)%prs_obs(icnt,1:nlay_obs)
      read(fileid,*,iostat=ios) &
         mls_data(i_min,j_min)%ret_lay_obs(icnt,1:ret_nlay)
      read(fileid,*,iostat=ios) &
         mls_data(i_min,j_min)%tru_lay_obs(icnt,1:tru_nlay)
      read(fileid,*,iostat=ios) &
         mls_data(i_min,j_min)%o3_obs(icnt,1:nlay_obs)
      read(fileid,*,iostat=ios) &
         mls_data(i_min,j_min)%o3_obs_err(icnt,1:nlay_obs)
      read(fileid,*,iostat=ios) &
         mls_data(i_min,j_min)%prior_obs(icnt,1:nlay_obs)
      read(fileid,*,iostat=ios) &
         mls_data(i_min,j_min)%prior_err(icnt)
      do k=1,ret_nlay
         read(fileid,*,iostat=ios) &
         mls_data(i_min,j_min)%avgk_obs(icnt,k,1:tru_nlay)
      enddo   
      read(fileid,*,iostat=ios) &
         mls_data(i_min,j_min)%o3_col_obs(icnt)
      read(fileid,*,iostat=ios) &
         mls_data(i_min,j_min)%o3_col_obs_err(icnt)
      read(fileid,*,iostat=ios) data_type, obs_id, &
         i_min, j_min
      mls_data(i_min,j_min)%data_type=data_type
      mls_data(i_min,j_min)%obs_id=obs_id
   enddo
!
! Calculate vertical column   
   do i=1,nx_model
      do j=1,ny_model
         do icnt=1,mls_data(i,j)%icnt
            call vertical_column(mls_data(i,j)%vert_col(icnt), &
            mls_data(i_min,j_min)%o3_obs(icnt,1:nlay_obs), &
            mls_data(i_min,j_min)%prs_obs(icnt,1:nlay_obs),nlay_obs)
         enddo
      enddo
   enddo
!
! Find the median
   mls_data(:,:)%median_col=0.   
   do i=1,nx_model
      do j=1,ny_model
         if(mls_data(i,j)%icnt.gt.1) then
            call median_code(mls_data(i,j)%median_col,mls_data(i,j)%vert_col(1:mls_data(i,j)%icnt), &
            mls_data(i,j)%icnt)
         endif
      enddo
   enddo
!
! Find the mode
   mls_data(:,:)%mode_col=0.   
   do i=1,nx_model
      do j=1,ny_model
         if(mls_data(i,j)%icnt.gt.20) then
            call mode_code(mls_data(i,j)%mode_col,mls_data(i,j)%vert_col(1:mls_data(i,j)%icnt), &
            mls_data(i,j)%icnt)
         endif
      enddo
   enddo
!
! Calculate mean
   do i=1,nx_model
      do j=1,ny_model
         if(mls_data(i,j)%icnt.gt.0) then
            mls_data(i,j)%mean_col=0.
            do icnt=1,mls_data(i,j)%icnt
               mls_data(i,j)%mean_col=mls_data(i,j)%mean_col+mls_data(i,j)%vert_col(icnt)
            enddo
            mls_data(i,j)%mean_col=mls_data(i,j)%mean_col/float(mls_data(i,j)%icnt)
         endif
      enddo
   enddo
!
! Calculate standard deviation   
   do i=1,nx_model
      do j=1,ny_model
         if(mls_data(i,j)%icnt.gt.1) then
            mls_data(i,j)%stdv_col=0.
            do icnt=1,mls_data(i,j)%icnt
               mls_data(i,j)%stdv_col=mls_data(i,j)%stdv_col+(mls_data(i,j)%vert_col(icnt)- &
               mls_data(i,j)%mean_col)**2.
            enddo
            mls_data(i,j)%stdv_col=sqrt(mls_data(i,j)%stdv_col/float(mls_data(i,j)%icnt-1))
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
         if(mls_data(i,j)%icnt.eq.1) then
            icnt_min(:,:)=1
         else if(mls_data(i,j)%icnt.gt.1) then
            do icnt=1,mls_data(i,j)%icnt
               dist=abs(mls_data(i,j)%mean_col-mls_data(i,j)%vert_col(icnt))
               if(dist.lt.dist_min(i,j)) then
                  dist_min(i,j)=dist
                  icnt_min(i,j)=icnt
               endif
            enddo
         endif
!
! Write data to file
         if(mls_data(i,j)%icnt.gt.0) then
            icnt=icnt_min(i,j)
            nlay_obs=mls_data(i,j)%nlay_obs(icnt)
            ret_nlay=mls_data(i,j)%ret_nlay(icnt)
            tru_nlay=mls_data(i,j)%tru_nlay(icnt)
            write(fileid,*,iostat=ios) &
               trim(mls_data(i,j)%data_type), &
               mls_data(i,j)%obs_id, i, j
            write(fileid,*,iostat=ios) &
               mls_data(i,j)%yr_obs(icnt), &
               mls_data(i,j)%mn_obs(icnt), &
               mls_data(i,j)%dy_obs(icnt), &
               mls_data(i,j)%hh_obs(icnt), &
               mls_data(i,j)%mm_obs(icnt), &
               mls_data(i,j)%ss_obs(icnt)
            write(fileid,*,iostat=ios) &
               mls_data(i,j)%lat_obs(icnt), &
               mls_data(i,j)%lon_obs(icnt)
            write(fileid,*,iostat=ios) &
               mls_data(i,j)%nlay_obs(icnt), &
               mls_data(i,j)%nlev_obs(icnt)
            write(fileid,*,iostat=ios) &
               mls_data(i,j)%ret_nlay(icnt), &
               mls_data(i,j)%tru_nlay(icnt)
            write(fileid,*,iostat=ios) &
               mls_data(i,j)%prs_obs(icnt,1:nlay_obs)
            write(fileid,*,iostat=ios) &
               mls_data(i,j)%ret_lay_obs(icnt,1:ret_nlay)
            write(fileid,*,iostat=ios) &
               mls_data(i,j)%tru_lay_obs(icnt,1:tru_nlay)
            write(fileid,*,iostat=ios) &
               mls_data(i,j)%o3_obs(icnt,1:nlay_obs)
            write(fileid,*,iostat=ios) &
               mls_data(i,j)%o3_obs_err(icnt,1:nlay_obs)
            write(fileid,*,iostat=ios) &
               mls_data(i,j)%prior_obs(icnt,1:nlay_obs)
            write(fileid,*,iostat=ios) &
               mls_data(i,j)%prior_err(icnt)
            do k=1,ret_nlay
               write(fileid,*,iostat=ios) &
               mls_data(i,j)%avgk_obs(icnt,k,1:tru_nlay)
            enddo
            write(fileid,*,iostat=ios) &
               mls_data(i,j)%o3_col_obs(icnt)
            write(fileid,*,iostat=ios) &
               mls_data(i,j)%o3_col_obs_err(icnt)
         endif
      enddo
   enddo
!               
! Deallocate structure arrays
   do i=1,nx_model
      do j=1,ny_model
         if(mls_data(i,j)%icnt.gt.0) then
            deallocate(mls_data(i,j)%yr_obs)
            deallocate(mls_data(i,j)%mn_obs)
            deallocate(mls_data(i,j)%dy_obs)
            deallocate(mls_data(i,j)%hh_obs)
            deallocate(mls_data(i,j)%mm_obs)
            deallocate(mls_data(i,j)%ss_obs)
            deallocate(mls_data(i,j)%lat_obs)
            deallocate(mls_data(i,j)%lon_obs)
            deallocate(mls_data(i,j)%nlay_obs)
            deallocate(mls_data(i,j)%nlev_obs)
            deallocate(mls_data(i,j)%ret_nlay)
            deallocate(mls_data(i,j)%tru_nlay)
            deallocate(mls_data(i,j)%o3_obs)
            deallocate(mls_data(i,j)%o3_obs_err)
            deallocate(mls_data(i,j)%prior_obs)
            deallocate(mls_data(i,j)%prior_err)
            deallocate(mls_data(i,j)%avgk_obs)
            deallocate(mls_data(i,j)%o3_col_obs)
            deallocate(mls_data(i,j)%o3_col_obs_err)
         endif
      enddo
   enddo
   deallocate (mls_data)
   deallocate (dist_min)
   deallocate (icnt_min)   
end program mls_o3_profile_thinner
