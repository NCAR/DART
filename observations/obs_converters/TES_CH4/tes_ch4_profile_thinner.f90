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
program tes_ch4_profile_thinner
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
   character(len=*), parameter      :: source   = 'tes_ch4_thinner.f90'
   character(len=*), parameter      :: revision = ''
   character(len=*), parameter      :: revdate  = ''
!
   integer,parameter                :: max_num_obs=200
   integer                          :: iunit,io,fileid,ios
   integer                          :: icnt
   integer                          :: i,j,k
!
! TES CH4 variable declarations
   integer                          :: obs_id,i_min,j_min
   integer                          :: yr_obs,mn_obs,dy_obs,hh_obs,mm_obs,ss_obs
   real                             :: lat_obs,lon_obs
   integer                          :: nlay_obs,nlev_obs,ndim_obs,ilv
   real                             :: dofs
   real,allocatable,dimension(:)    :: prs_lev
   real,allocatable,dimension(:)    :: ch4_lay
   real,allocatable,dimension(:)    :: ch4_prior_lay
   real,allocatable,dimension(:)    :: ch4_prior_err_lay
   real,allocatable,dimension(:,:)  :: avgk_lay
   real,allocatable,dimension(:)    :: cov_lay
   real,allocatable,dimension(:)    :: cov_prior_lay
!
! Namelist variable declarations
   character(len=129)               :: filedir,filename,fileout
   integer                          :: year,month,day,hour
   real                             :: bin_beg_sec,bin_end_sec
   real                             :: fac_obs_error
   logical                          :: use_log_co,use_log_o3,use_log_co2,use_log_ch4,use_log_nh3
   real                             :: lon_min,lon_max,lat_min,lat_max
   character(len=129)               :: path_model,file_model,data_type
   integer                          :: nx_model,ny_model,nz_model
   integer                          :: obs_co_reten_freq,obs_o3_reten_freq,obs_co2_reten_freq, &
                                       obs_ch4_reten_freq,obs_nh3_reten_freq
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
      real,allocatable,dimension(:,:)   :: prs_obs
      real,allocatable,dimension(:)     :: trop_prs
      real,allocatable,dimension(:,:)   :: ch4_obs
      real,allocatable,dimension(:,:)   :: ch4_obs_prior
      real,allocatable,dimension(:,:)   :: ch4_obs_err
      real,allocatable,dimension(:)     :: dofs
      real,allocatable,dimension(:,:,:) :: avgk_obs
      real,allocatable,dimension(:,:)   :: avgk_diag_obs
      real,allocatable,dimension(:,:,:) :: cov_obs
      real,allocatable,dimension(:,:,:) :: cov_total
      real,allocatable,dimension(:)     :: ch4_total_col
      real,allocatable,dimension(:)     :: ch4_total_col_prior
      real,allocatable,dimension(:)     :: ch4_total_col_err
   end type observation_data
!
   type(observation_data),allocatable,dimension(:,:) :: tes_data
   real                                :: dist
   real,allocatable,dimension(:,:)     :: dist_min
   real,allocatable,dimension(:,:)     :: icnt_min
!
   namelist /create_tes_obs_nml/filedir,filename,fileout, &
   bin_beg_sec,bin_end_sec,fac_obs_error,use_log_co,use_log_o3,use_log_co2,use_log_ch4, &
   use_log_nh3,lon_min,lon_max,lat_min,lat_max,path_model,file_model,nx_model,ny_model, &
   nz_model,obs_co_reten_freq,obs_o3_reten_freq,obs_co2_reten_freq,obs_ch4_reten_freq, &
   obs_nh3_reten_freq
!
! Record the current time, date, etc. to the logfile
   call initialize_utilities(source)
   call register_module(source,revision,revdate)
!
   call find_namelist_in_file("input.nml", "create_tes_obs_nml", iunit)
   read(iunit, nml = create_tes_obs_nml, iostat = io)
   call check_namelist_read(iunit, io, "create_tes_obs_nml")
!
! Read TES CH4
   allocate(tes_data(nx_model,ny_model))
   allocate(dist_min(nx_model,ny_model))
   allocate(icnt_min(nx_model,ny_model))
!
   tes_data(:,:)%icnt=0
   fileid=100
   write(6,*)'opening ',TRIM(TRIM(filedir)//TRIM(filename))
   open(unit=fileid,file=TRIM(TRIM(filedir)//TRIM(filename)), &
   form='formatted', status='old', iostat=ios)
   read(fileid,*,iostat=ios) data_type, obs_id, &
   i_min, j_min
   tes_data(i_min,j_min)%data_type=data_type
   tes_data(i_min,j_min)%obs_id=obs_id
   do while (ios == 0)
      tes_data(i_min,j_min)%icnt=tes_data(i_min,j_min)%icnt+1
      icnt=tes_data(i_min,j_min)%icnt
      if(tes_data(i_min,j_min)%icnt.eq.1) then
         allocate(tes_data(i_min,j_min)%vert_col(max_num_obs))
         allocate(tes_data(i_min,j_min)%yr_obs(max_num_obs))
         allocate(tes_data(i_min,j_min)%mn_obs(max_num_obs))
         allocate(tes_data(i_min,j_min)%dy_obs(max_num_obs))
         allocate(tes_data(i_min,j_min)%hh_obs(max_num_obs))
         allocate(tes_data(i_min,j_min)%mm_obs(max_num_obs))
         allocate(tes_data(i_min,j_min)%ss_obs(max_num_obs))
         allocate(tes_data(i_min,j_min)%lat_obs(max_num_obs))
         allocate(tes_data(i_min,j_min)%lon_obs(max_num_obs))
         allocate(tes_data(i_min,j_min)%nlay_obs(max_num_obs))
         allocate(tes_data(i_min,j_min)%nlev_obs(max_num_obs))
         allocate(tes_data(i_min,j_min)%trop_prs(max_num_obs))
         allocate(tes_data(i_min,j_min)%dofs(max_num_obs))
         allocate(tes_data(i_min,j_min)%ch4_total_col(max_num_obs))
         allocate(tes_data(i_min,j_min)%ch4_total_col_prior(max_num_obs))
         allocate(tes_data(i_min,j_min)%ch4_total_col_err(max_num_obs))
      endif
!
      read(fileid,*,iostat=ios) &
         tes_data(i_min,j_min)%yr_obs(icnt), &
         tes_data(i_min,j_min)%mn_obs(icnt), &
         tes_data(i_min,j_min)%dy_obs(icnt), &
         tes_data(i_min,j_min)%hh_obs(icnt), &
         tes_data(i_min,j_min)%mm_obs(icnt), &
         tes_data(i_min,j_min)%ss_obs(icnt)      
      read(fileid,*,iostat=ios) &
         tes_data(i_min,j_min)%lat_obs(icnt), &
         tes_data(i_min,j_min)%lon_obs(icnt)
      read(fileid,*,iostat=ios) &
         tes_data(i_min,j_min)%nlay_obs(icnt), &
         tes_data(i_min,j_min)%nlev_obs(icnt)
         nlay_obs=tes_data(i_min,j_min)%nlay_obs(icnt)
         nlev_obs=tes_data(i_min,j_min)%nlev_obs(icnt)
      if(tes_data(i_min,j_min)%icnt.eq.1) then
         allocate(tes_data(i_min,j_min)%prs_obs(max_num_obs,nlay_obs))
         allocate(tes_data(i_min,j_min)%ch4_obs(max_num_obs,nlay_obs))
         allocate(tes_data(i_min,j_min)%ch4_obs_prior(max_num_obs,nlay_obs))
         allocate(tes_data(i_min,j_min)%ch4_obs_err(max_num_obs,nlay_obs))
         allocate(tes_data(i_min,j_min)%avgk_obs(max_num_obs,nlay_obs,nlay_obs))
         allocate(tes_data(i_min,j_min)%avgk_diag_obs(max_num_obs,nlay_obs))
         allocate(tes_data(i_min,j_min)%cov_obs(max_num_obs,nlay_obs,nlay_obs))
         allocate(tes_data(i_min,j_min)%cov_total(max_num_obs,nlay_obs,nlay_obs))
      endif
      read(fileid,*,iostat=ios) &
         tes_data(i_min,j_min)%prs_obs(icnt,1:nlay_obs)
      read(fileid,*,iostat=ios) &
         tes_data(i_min,j_min)%trop_prs(icnt)
      read(fileid,*,iostat=ios) &
         tes_data(i_min,j_min)%ch4_obs(icnt,1:nlay_obs)
      read(fileid,*,iostat=ios) &
         tes_data(i_min,j_min)%ch4_obs_prior(icnt,1:nlay_obs)
      read(fileid,*,iostat=ios) &
         tes_data(i_min,j_min)%ch4_obs_err(icnt,1:nlay_obs)
      read(fileid,*,iostat=ios) &
         tes_data(i_min,j_min)%dofs(icnt)
      do k=1,nlay_obs
         read(fileid,*,iostat=ios) &
         tes_data(i_min,j_min)%avgk_obs(icnt,1:nlay_obs,k)
      enddo   
      read(fileid,*,iostat=ios) &
         tes_data(i_min,j_min)%avgk_diag_obs(icnt,1:nlay_obs)
      do k=1,nlay_obs
         read(fileid,*,iostat=ios) &
         tes_data(i_min,j_min)%cov_obs(icnt,1:nlay_obs,k)
      enddo   
      do k=1,nlay_obs
         read(fileid,*,iostat=ios) &
         tes_data(i_min,j_min)%cov_total(icnt,1:nlay_obs,k)
      enddo   
      read(fileid,*,iostat=ios) &
         tes_data(i_min,j_min)%ch4_total_col(icnt)
      read(fileid,*,iostat=ios) &
         tes_data(i_min,j_min)%ch4_total_col_prior(icnt)
      read(fileid,*,iostat=ios) &
         tes_data(i_min,j_min)%ch4_total_col_err(icnt)
      read(fileid,*,iostat=ios) data_type, obs_id, &
         i_min, j_min
      tes_data(i_min,j_min)%data_type=data_type
      tes_data(i_min,j_min)%obs_id=obs_id
   enddo
!
! Calculate vertical column   
   do i=1,nx_model
      do j=1,ny_model
         do icnt=1,tes_data(i,j)%icnt
            call vertical_column(tes_data(i,j)%vert_col(icnt), &
            tes_data(i_min,j_min)%ch4_obs(icnt,1:nlay_obs), &
            tes_data(i_min,j_min)%prs_obs(icnt,1:nlay_obs),nlay_obs)
         enddo
      enddo
   enddo
!
! Find the median
   tes_data(:,:)%median_col=0.   
   do i=1,nx_model
      do j=1,ny_model
         if(tes_data(i,j)%icnt.gt.1) then
            call median_code(tes_data(i,j)%median_col,tes_data(i,j)%vert_col(1:tes_data(i,j)%icnt), &
            tes_data(i,j)%icnt)
         endif
      enddo
   enddo
!
! Find the mode
   tes_data(:,:)%mode_col=0.   
   do i=1,nx_model
      do j=1,ny_model
         if(tes_data(i,j)%icnt.gt.20) then
            call mode_code(tes_data(i,j)%mode_col,tes_data(i,j)%vert_col(1:tes_data(i,j)%icnt), &
            tes_data(i,j)%icnt)
         endif
      enddo
   enddo
!
! Calculate mean
   do i=1,nx_model
      do j=1,ny_model
         if(tes_data(i,j)%icnt.gt.0) then
            tes_data(i,j)%mean_col=0.
            do icnt=1,tes_data(i,j)%icnt
               tes_data(i,j)%mean_col=tes_data(i,j)%mean_col+tes_data(i,j)%vert_col(icnt)
            enddo
            tes_data(i,j)%mean_col=tes_data(i,j)%mean_col/float(tes_data(i,j)%icnt)
         endif
      enddo
   enddo
!
! Calculate standard deviation   
   do i=1,nx_model
      do j=1,ny_model
         if(tes_data(i,j)%icnt.gt.1) then
            tes_data(i,j)%stdv_col=0.
            do icnt=1,tes_data(i,j)%icnt
               tes_data(i,j)%stdv_col=tes_data(i,j)%stdv_col+(tes_data(i,j)%vert_col(icnt)- &
               tes_data(i,j)%mean_col)**2.
            enddo
            tes_data(i,j)%stdv_col=sqrt(tes_data(i,j)%stdv_col/float(tes_data(i,j)%icnt-1))
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
         if(tes_data(i,j)%icnt.eq.1) then
            icnt_min(:,:)=1
         else if(tes_data(i,j)%icnt.gt.1) then
            do icnt=1,tes_data(i,j)%icnt
               dist=abs(tes_data(i,j)%mean_col-tes_data(i,j)%vert_col(icnt))
               if(dist.lt.dist_min(i,j)) then
                  dist_min(i,j)=dist
                  icnt_min(i,j)=icnt
               endif
            enddo
         endif
!
! Write data to file
         if(tes_data(i,j)%icnt.gt.0) then
            icnt=icnt_min(i,j)
            write(fileid,*,iostat=ios) &
               trim(tes_data(i,j)%data_type), &
               tes_data(i,j)%obs_id, i, j
            write(fileid,*,iostat=ios) &
               tes_data(i,j)%yr_obs(icnt), &
               tes_data(i,j)%mn_obs(icnt), &
               tes_data(i,j)%dy_obs(icnt), &
               tes_data(i,j)%hh_obs(icnt), &
               tes_data(i,j)%mm_obs(icnt), &
               tes_data(i,j)%ss_obs(icnt)
            write(fileid,*,iostat=ios) &
               tes_data(i,j)%lat_obs(icnt), &
               tes_data(i,j)%lon_obs(icnt)
            write(fileid,*,iostat=ios) &
               tes_data(i,j)%nlay_obs(icnt), &
               tes_data(i,j)%nlev_obs(icnt)
            write(fileid,*,iostat=ios) &
               tes_data(i,j)%prs_obs(icnt,1:nlay_obs)
            write(fileid,*,iostat=ios) &
               tes_data(i,j)%trop_prs(icnt)
            write(fileid,*,iostat=ios) &
               tes_data(i,j)%ch4_obs(icnt,1:nlay_obs)
            write(fileid,*,iostat=ios) &
               tes_data(i,j)%ch4_obs_prior(icnt,1:nlay_obs)
            write(fileid,*,iostat=ios) &
               tes_data(i,j)%ch4_obs_err(icnt,1:nlay_obs)
            write(fileid,*,iostat=ios) &
               tes_data(i,j)%dofs(icnt)
            do k=1,nlay_obs
               write(fileid,*,iostat=ios) &
               tes_data(i,j)%avgk_obs(icnt,1:nlay_obs,k)
            enddo
            write(fileid,*,iostat=ios) &
               tes_data(i,j)%avgk_diag_obs(icnt,1:nlay_obs)
            do k=1,nlay_obs
               write(fileid,*,iostat=ios) &
               tes_data(i,j)%cov_obs(icnt,1:nlay_obs,k)
            enddo
            do k=1,nlay_obs
               write(fileid,*,iostat=ios) &
               tes_data(i,j)%cov_total(icnt,1:nlay_obs,k)
            enddo
            write(fileid,*,iostat=ios) &
               tes_data(i,j)%ch4_total_col(icnt)
            write(fileid,*,iostat=ios) &
               tes_data(i,j)%ch4_total_col_prior(icnt)
            write(fileid,*,iostat=ios) &
               tes_data(i,j)%ch4_total_col_err(icnt)
         endif
      enddo
   enddo
!               
! Deallocate structure arrays
   do i=1,nx_model
      do j=1,ny_model
         if(tes_data(i,j)%icnt.gt.0) then
            deallocate(tes_data(i,j)%yr_obs)
            deallocate(tes_data(i,j)%mn_obs)
            deallocate(tes_data(i,j)%dy_obs)
            deallocate(tes_data(i,j)%hh_obs)
            deallocate(tes_data(i,j)%mm_obs)
            deallocate(tes_data(i,j)%ss_obs)
            deallocate(tes_data(i,j)%lat_obs)
            deallocate(tes_data(i,j)%lon_obs)
            deallocate(tes_data(i,j)%nlay_obs)
            deallocate(tes_data(i,j)%nlev_obs)
            deallocate(tes_data(i,j)%prs_obs)
            deallocate(tes_data(i,j)%trop_prs)
            deallocate(tes_data(i,j)%ch4_obs)
            deallocate(tes_data(i,j)%ch4_obs_prior)
            deallocate(tes_data(i,j)%ch4_obs_err)
            deallocate(tes_data(i,j)%dofs)
            deallocate(tes_data(i,j)%avgk_obs)
            deallocate(tes_data(i,j)%avgk_diag_obs)
            deallocate(tes_data(i,j)%cov_obs)
            deallocate(tes_data(i,j)%cov_total)
            deallocate(tes_data(i,j)%ch4_total_col)
            deallocate(tes_data(i,j)%ch4_total_col_prior)
            deallocate(tes_data(i,j)%ch4_total_col_err)
         endif
      enddo
   enddo
   deallocate (tes_data)
   deallocate (dist_min)
   deallocate (icnt_min)   
end program tes_ch4_profile_thinner
