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
program mopitt_v8_co_profile_thinner
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
   character(len=*), parameter      :: source   = 'mopitt_v8_co_thinner.f90'
   character(len=*), parameter      :: revision = ''
   character(len=*), parameter      :: revdate  = ''
!
   integer,parameter                :: max_num_obs=200
   integer                          :: iunit,io,fileid,ios
   integer                          :: icnt
   integer                          :: i,j,k
!
! MOPITT CO variable declarations
   integer                          :: obs_id,i_min,j_min
   integer                          :: yr_obs,mn_obs,dy_obs,hh_obs,mm_obs,ss_obs
   real                             :: lat_obs,lon_obs
   integer                          :: nlay_obs,nlev_obs,ilv
   real                             :: dofs
   real,allocatable,dimension(:)    :: prs_lay,prs_lev
   real,allocatable,dimension(:)    :: co_lay, co_lay_err
   real,allocatable,dimension(:)    :: co_prior_lay,co_prior_lay_err
   real,allocatable,dimension(:,:)  :: avgk_lay,cov_s,cov_r,cov_m
   real                             :: co_col,co_col_err,co_col_prior
!
! Namelist variable declarations
   character(len=129)               :: filedir,filename,fileout
   integer                          :: year,month,day,hour
   real                             :: bin_beg_sec,bin_end_sec
   real                             :: fac_obs_error
   logical                          :: use_log_co,use_log_o3
   real                             :: lon_min,lon_max,lat_min,lat_max
   character(len=129)               :: path_model,file_model,data_type
   integer                          :: nx_model,ny_model,nz_model
   integer                          :: obs_co_reten_freq,obs_o3_reten_freq
!  
   type observation_data
      integer                           :: icnt
      real                              :: mean_col
      real                              :: stdv_col
      real                              :: median_col      
      real                              :: mode_col      
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
      real,allocatable,dimension(:)     :: dofs
      real,allocatable,dimension(:)     :: prs_sfc
      real,allocatable,dimension(:,:)   :: prs_lay
      real,allocatable,dimension(:,:)   :: prs_lev
      real,allocatable,dimension(:,:,:) :: avgk_lay
      real,allocatable,dimension(:,:)   :: co_lay
      real,allocatable,dimension(:,:)   :: co_lay_err
      real,allocatable,dimension(:,:)   :: co_prior_lay
      real,allocatable,dimension(:,:)   :: co_prior_lay_err
      real,allocatable,dimension(:,:,:) :: cov_s
      real,allocatable,dimension(:,:,:) :: cov_r
      real,allocatable,dimension(:,:,:) :: cov_m
      real,allocatable,dimension(:)     :: co_col
      real,allocatable,dimension(:)     :: co_col_err
      real,allocatable,dimension(:)     :: co_col_prior
   end type observation_data
!
   type(observation_data),allocatable,dimension(:,:) :: mopitt_data
   real                                :: dist
   real,allocatable,dimension(:,:)     :: dist_min
   real,allocatable,dimension(:,:)     :: icnt_min
!
   namelist /create_mopitt_obs_nml/filedir,filename,fileout,year,month,day,hour, &
   bin_beg_sec,bin_end_sec,fac_obs_error,use_log_co,use_log_o3, &
   lon_min,lon_max,lat_min,lat_max,path_model,file_model,nx_model, &
   ny_model,nz_model,obs_co_reten_freq,obs_o3_reten_freq
!
! Record the current time, date, etc. to the logfile
   call initialize_utilities(source)
   call register_module(source,revision,revdate)
!
   call find_namelist_in_file("input.nml", "create_mopitt_obs_nml", iunit)
   read(iunit, nml = create_mopitt_obs_nml, iostat = io)
   call check_namelist_read(iunit, io, "create_mopitt_obs_nml")
!
! Read MOPITT CO
   allocate(mopitt_data(nx_model,ny_model))
   allocate(dist_min(nx_model,ny_model))
   allocate(icnt_min(nx_model,ny_model))
!
   mopitt_data(:,:)%icnt=0
   fileid=100
   write(6,*)'opening ',TRIM(TRIM(filedir)//TRIM(filename))
   open(unit=fileid,file=TRIM(TRIM(filedir)//TRIM(filename)), &
   form='formatted', status='old', iostat=ios)
   read(fileid,*,iostat=ios) data_type, obs_id, i_min, j_min
   mopitt_data(i_min,j_min)%data_type=data_type
   mopitt_data(i_min,j_min)%obs_id=obs_id
   do while (ios == 0)
      mopitt_data(i_min,j_min)%icnt=mopitt_data(i_min,j_min)%icnt+1
      icnt=mopitt_data(i_min,j_min)%icnt
      if(mopitt_data(i_min,j_min)%icnt.eq.1) then
         allocate(mopitt_data(i_min,j_min)%yr_obs(max_num_obs))
         allocate(mopitt_data(i_min,j_min)%mn_obs(max_num_obs))
         allocate(mopitt_data(i_min,j_min)%dy_obs(max_num_obs))
         allocate(mopitt_data(i_min,j_min)%hh_obs(max_num_obs))
         allocate(mopitt_data(i_min,j_min)%mm_obs(max_num_obs))
         allocate(mopitt_data(i_min,j_min)%ss_obs(max_num_obs))
         allocate(mopitt_data(i_min,j_min)%lat_obs(max_num_obs))
         allocate(mopitt_data(i_min,j_min)%lon_obs(max_num_obs))
         allocate(mopitt_data(i_min,j_min)%nlay_obs(max_num_obs))
         allocate(mopitt_data(i_min,j_min)%nlev_obs(max_num_obs))
         allocate(mopitt_data(i_min,j_min)%dofs(max_num_obs))
         allocate(mopitt_data(i_min,j_min)%prs_sfc(max_num_obs))
         allocate(mopitt_data(i_min,j_min)%co_col(max_num_obs))
         allocate(mopitt_data(i_min,j_min)%co_col_err(max_num_obs))
         allocate(mopitt_data(i_min,j_min)%co_col_prior(max_num_obs))
      endif
!
      read(fileid,*,iostat=ios) &
         mopitt_data(i_min,j_min)%yr_obs(icnt), &
         mopitt_data(i_min,j_min)%mn_obs(icnt), &
         mopitt_data(i_min,j_min)%dy_obs(icnt), &
         mopitt_data(i_min,j_min)%hh_obs(icnt), &
         mopitt_data(i_min,j_min)%mm_obs(icnt), &
         mopitt_data(i_min,j_min)%ss_obs(icnt)      
      read(fileid,*,iostat=ios) &
         mopitt_data(i_min,j_min)%lat_obs(icnt), &
         mopitt_data(i_min,j_min)%lon_obs(icnt)
      read(fileid,*,iostat=ios) &
         mopitt_data(i_min,j_min)%nlay_obs(icnt), &
         mopitt_data(i_min,j_min)%nlev_obs(icnt)
         nlay_obs=mopitt_data(i_min,j_min)%nlay_obs(icnt)
         nlev_obs=mopitt_data(i_min,j_min)%nlev_obs(icnt)
      read(fileid,*,iostat=ios) &
         mopitt_data(i_min,j_min)%dofs(icnt)
      read(fileid,*,iostat=ios) &
         mopitt_data(i_min,j_min)%prs_sfc(icnt)
      if(mopitt_data(i_min,j_min)%icnt.eq.1) then
         allocate(mopitt_data(i_min,j_min)%prs_lay(max_num_obs,nlay_obs))
         allocate(mopitt_data(i_min,j_min)%prs_lev(max_num_obs,nlev_obs))
         allocate(mopitt_data(i_min,j_min)%avgk_lay(max_num_obs,nlay_obs,nlay_obs))
         allocate(mopitt_data(i_min,j_min)%co_lay(max_num_obs,nlay_obs))
         allocate(mopitt_data(i_min,j_min)%co_lay_err(max_num_obs,nlay_obs))
         allocate(mopitt_data(i_min,j_min)%co_prior_lay(max_num_obs,nlay_obs))
         allocate(mopitt_data(i_min,j_min)%co_prior_lay_err(max_num_obs,nlay_obs))
         allocate(mopitt_data(i_min,j_min)%cov_s(max_num_obs,nlay_obs,nlay_obs))
         allocate(mopitt_data(i_min,j_min)%cov_r(max_num_obs,nlay_obs,nlay_obs))
         allocate(mopitt_data(i_min,j_min)%cov_m(max_num_obs,nlay_obs,nlay_obs))
      endif
      read(fileid,*,iostat=ios) &
         mopitt_data(i_min,j_min)%prs_lay(icnt,1:nlay_obs)
      read(fileid,*,iostat=ios) &
         mopitt_data(i_min,j_min)%prs_lev(icnt,1:nlev_obs)
      do k=1,nlay_obs
         read(fileid,*,iostat=ios) &
         mopitt_data(i_min,j_min)%avgk_lay(icnt,1:nlay_obs,k)
      enddo   
      read(fileid,*,iostat=ios) &
         mopitt_data(i_min,j_min)%co_lay(icnt,1:nlay_obs)
      read(fileid,*,iostat=ios) &
         mopitt_data(i_min,j_min)%co_lay_err(icnt,1:nlay_obs)
      read(fileid,*,iostat=ios) &
         mopitt_data(i_min,j_min)%co_prior_lay(icnt,1:nlay_obs)
      read(fileid,*,iostat=ios) &
         mopitt_data(i_min,j_min)%co_prior_lay_err(icnt,1:nlay_obs)
      do k=1,nlay_obs
         read(fileid,*,iostat=ios) &
         mopitt_data(i_min,j_min)%cov_s(icnt,1:nlay_obs,k)
      enddo   
      do k=1,nlay_obs
         read(fileid,*,iostat=ios) &
         mopitt_data(i_min,j_min)%cov_r(icnt,1:nlay_obs,k)
      enddo   
      do k=1,nlay_obs
         read(fileid,*,iostat=ios) &
         mopitt_data(i_min,j_min)%cov_m(icnt,1:nlay_obs,k)
      enddo   
      read(fileid,*,iostat=ios) &
         mopitt_data(i_min,j_min)%co_col(icnt)
      read(fileid,*,iostat=ios) &
         mopitt_data(i_min,j_min)%co_col_err(icnt)
      read(fileid,*,iostat=ios) &
         mopitt_data(i_min,j_min)%co_col_prior(icnt)
      
      read(fileid,*,iostat=ios) data_type, obs_id, i_min, j_min
      mopitt_data(i_min,j_min)%data_type=data_type
      mopitt_data(i_min,j_min)%obs_id=obs_id
   enddo
!
! Calculate vertical column (not needed MOPITT has vertical column)
!   do i=1,nx_model
!      do j=1,ny_model
!         do icnt=1,mopitt_data(i,j)%icnt
!            call vertical_column(mopitt_data(i,j)%vert_col(icnt), &
!            mopitt_data(i_min,j_min)%co_lay(icnt,1:nlay_obs), &
!            mopitt_data(i_min,j_min)%prs_lev(icnt,1:nlev_obs),nlay_obs)
!         enddo
!      enddo
!   enddo
!
! Find the median
   mopitt_data(:,:)%median_col=0.   
   do i=1,nx_model
      do j=1,ny_model
         if(mopitt_data(i,j)%icnt.gt.1) then
            call median_code(mopitt_data(i,j)%median_col,mopitt_data(i,j)%co_col(1:mopitt_data(i,j)%icnt), &
            mopitt_data(i,j)%icnt)
         endif
      enddo
   enddo
!
! Find the mode
   mopitt_data(:,:)%mode_col=0.   
   do i=1,nx_model
      do j=1,ny_model
         if(mopitt_data(i,j)%icnt.gt.20) then
            call mode_code(mopitt_data(i,j)%mode_col,mopitt_data(i,j)%co_col(1:mopitt_data(i,j)%icnt), &
            mopitt_data(i,j)%icnt)
         endif
      enddo
   enddo
!
! Calculate mean
   do i=1,nx_model
      do j=1,ny_model
         if(mopitt_data(i,j)%icnt.gt.0) then
            mopitt_data(i,j)%mean_col=0.
            do icnt=1,mopitt_data(i,j)%icnt
               mopitt_data(i,j)%mean_col=mopitt_data(i,j)%mean_col+mopitt_data(i,j)%co_col(icnt)
            enddo
            mopitt_data(i,j)%mean_col=mopitt_data(i,j)%mean_col/float(mopitt_data(i,j)%icnt)
         endif
      enddo
   enddo
!
! Calculate standard deviation   
   do i=1,nx_model
      do j=1,ny_model
         if(mopitt_data(i,j)%icnt.gt.1) then
            mopitt_data(i,j)%stdv_col=0.
            do icnt=1,mopitt_data(i,j)%icnt
               mopitt_data(i,j)%stdv_col=mopitt_data(i,j)%stdv_col+(mopitt_data(i,j)%co_col(icnt)- &
               mopitt_data(i,j)%mean_col)**2.
            enddo
            mopitt_data(i,j)%stdv_col=sqrt(mopitt_data(i,j)%stdv_col/float(mopitt_data(i,j)%icnt-1))
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
         if(mopitt_data(i,j)%icnt.eq.1) then
            icnt_min(:,:)=1
         else if(mopitt_data(i,j)%icnt.gt.1) then
            do icnt=1,mopitt_data(i,j)%icnt
               dist=abs(mopitt_data(i,j)%mean_col-mopitt_data(i,j)%co_col(icnt))
               if(dist.lt.dist_min(i,j)) then
                  dist_min(i,j)=dist
                  icnt_min(i,j)=icnt
               endif
            enddo
         endif
!
! Write data to file
         if(mopitt_data(i,j)%icnt.gt.0) then
            icnt=icnt_min(i,j)
            nlay_obs=mopitt_data(i,j)%nlay_obs(icnt)
            nlev_obs=mopitt_data(i,j)%nlev_obs(icnt)
            write(fileid,*,iostat=ios) &
               trim(mopitt_data(i,j)%data_type), &
               mopitt_data(i,j)%obs_id, i, j
            write(fileid,*,iostat=ios) &
               mopitt_data(i,j)%yr_obs(icnt), &
               mopitt_data(i,j)%mn_obs(icnt), &
               mopitt_data(i,j)%dy_obs(icnt), &
               mopitt_data(i,j)%hh_obs(icnt), &
               mopitt_data(i,j)%mm_obs(icnt), &
               mopitt_data(i,j)%ss_obs(icnt)
            write(fileid,*,iostat=ios) &
               mopitt_data(i,j)%lat_obs(icnt), &
               mopitt_data(i,j)%lon_obs(icnt)
            write(fileid,*,iostat=ios) &
               mopitt_data(i,j)%nlay_obs(icnt), &
               mopitt_data(i,j)%nlev_obs(icnt)
            write(fileid,*,iostat=ios) &
               mopitt_data(i,j)%dofs(icnt)
            write(fileid,*,iostat=ios) &
               mopitt_data(i,j)%prs_sfc(icnt)
            write(fileid,*,iostat=ios) &
               mopitt_data(i,j)%prs_lay(icnt,1:nlay_obs)
            write(fileid,*,iostat=ios) &
               mopitt_data(i,j)%prs_lev(icnt,1:nlev_obs)
            do k=1,nlay_obs
               write(fileid,*,iostat=ios) &
               mopitt_data(i,j)%avgk_lay(icnt,1:nlay_obs,k)
            enddo
            
            do ilv=1,nlay_obs
               if(mopitt_data(i,j)%co_lay(icnt,ilv).eq.0.) then
                  print *, 'obs_id,i,j ',mopitt_data(i,j)%obs_id,i,j
                  print *, 'yr, mo, dy, mn, sec ', &
                  mopitt_data(i,j)%yr_obs(icnt), &
                  mopitt_data(i,j)%mn_obs(icnt), &
                  mopitt_data(i,j)%dy_obs(icnt), &
                  mopitt_data(i,j)%hh_obs(icnt), &
                  mopitt_data(i,j)%mm_obs(icnt), &
                  mopitt_data(i,j)%ss_obs(icnt)
                  print *, 'co_lay ',mopitt_data(i,j)%co_lay(icnt,1:nlay_obs)
                  stop
               endif
            enddo
               
            write(fileid,*,iostat=ios) &
               mopitt_data(i,j)%co_lay(icnt,1:nlay_obs)
            write(fileid,*,iostat=ios) &
               mopitt_data(i,j)%co_lay_err(icnt,1:nlay_obs)
            write(fileid,*,iostat=ios) &
               mopitt_data(i,j)%co_prior_lay(icnt,1:nlay_obs)
            write(fileid,*,iostat=ios) &
               mopitt_data(i,j)%co_prior_lay_err(icnt,1:nlay_obs)
            do k=1,nlay_obs
               write(fileid,*,iostat=ios) &
               mopitt_data(i,j)%cov_s(icnt,1:nlay_obs,k)
            enddo
            do k=1,nlay_obs
               write(fileid,*,iostat=ios) &
               mopitt_data(i,j)%cov_r(icnt,1:nlay_obs,k)
            enddo
            do k=1,nlay_obs
               write(fileid,*,iostat=ios) &
               mopitt_data(i,j)%cov_m(icnt,1:nlay_obs,k)
            enddo
            write(fileid,*,iostat=ios) &
               mopitt_data(i,j)%co_col(icnt)
            write(fileid,*,iostat=ios) &
               mopitt_data(i,j)%co_col_err(icnt)
            write(fileid,*,iostat=ios) &
               mopitt_data(i,j)%co_col_prior(icnt)
         endif
      enddo
   enddo
!               
! Deallocate structure arrays
   do i=1,nx_model
      do j=1,ny_model
         if(mopitt_data(i,j)%icnt.gt.0) then
            deallocate(mopitt_data(i,j)%yr_obs)
            deallocate(mopitt_data(i,j)%mn_obs)
            deallocate(mopitt_data(i,j)%dy_obs)
            deallocate(mopitt_data(i,j)%hh_obs)
            deallocate(mopitt_data(i,j)%mm_obs)
            deallocate(mopitt_data(i,j)%ss_obs)
            deallocate(mopitt_data(i,j)%lat_obs)
            deallocate(mopitt_data(i,j)%lon_obs)
            deallocate(mopitt_data(i,j)%nlay_obs)
            deallocate(mopitt_data(i,j)%nlev_obs)
            deallocate(mopitt_data(i,j)%dofs)
            deallocate(mopitt_data(i,j)%prs_sfc)
            deallocate(mopitt_data(i,j)%co_col)
            deallocate(mopitt_data(i,j)%co_col_err)
            deallocate(mopitt_data(i,j)%co_col_prior)
            deallocate(mopitt_data(i,j)%prs_lay)
            deallocate(mopitt_data(i,j)%prs_lev)
            deallocate(mopitt_data(i,j)%avgk_lay)
            deallocate(mopitt_data(i,j)%co_lay)
            deallocate(mopitt_data(i,j)%co_lay_err)
            deallocate(mopitt_data(i,j)%co_prior_lay)
            deallocate(mopitt_data(i,j)%co_prior_lay_err)
            deallocate(mopitt_data(i,j)%cov_s)
            deallocate(mopitt_data(i,j)%cov_r)
            deallocate(mopitt_data(i,j)%cov_m)
         endif
      enddo
   enddo
   deallocate (mopitt_data)
   deallocate (dist_min)
   deallocate (icnt_min)   
end program mopitt_v8_co_profile_thinner
