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
program gome2a_no2_thinner
   use    utilities_mod, only     : register_module,            &
                                    initialize_utilities,       &
                                    find_namelist_in_file,      &
                                    check_namelist_read,        &
                                    open_file,                  &
                                    close_file                 
   implicit none
!
! Version controlled file description for error handling, do not edit
   character(len=*), parameter     :: source   = 'gome2a_no2_thinning.f90'
   character(len=*), parameter     :: revision = ''
   character(len=*), parameter     :: revdate  = ''
!
   integer,parameter               :: max_num_obs=200
   integer                         :: iunit,io,fileid,ios
   integer                         :: icnt
   integer                         :: i,j
!
! GOME2A NO2 variable declarations
   integer                         :: obs_id,i_min,j_min
   integer                         :: yr_obs,mn_obs,dy_obs,hh_obs,mm_obs,ss_obs
   real                            :: lat_obs,lon_obs
   integer                         :: nlay_obs,nlev_obs,ilv
   integer                         :: trop_index
   real                            :: no2_trop_col_obs,no2_trop_col_obs_err
   real                            :: no2_total_col_obs,no2_total_col_obs_err
   real                            :: amf_trop_col
   real                            :: amf_total_col
   real                            :: no2_slnt_col,no2_slnt_col_err
   real                            :: o3_slnt_col,o3_slnt_col_err
   real,allocatable,dimension(:)   :: avgk_obs,prs_obs
!
! Namelist variable declarations
   character(len=129)              :: filedir,filename,fileout
   integer                         :: year,month,day,hour
   real                            :: bin_beg_sec,bin_end_sec
   real                            :: fac_obs_error
   logical                         :: use_log_no2
   real                            :: lon_min,lon_max,lat_min,lat_max
   character(len=129)              :: path_model,file_model,data_type
   integer                         :: nx_model,ny_model,nz_model
   integer                         :: obs_no2_reten_freq
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
      integer,allocatable,dimension(:) :: trop_index
      real,allocatable,dimension(:)    :: no2_trop_col_obs
      real,allocatable,dimension(:)    :: no2_trop_col_obs_err
      real,allocatable,dimension(:)    :: no2_total_col_obs
      real,allocatable,dimension(:)    :: no2_total_col_obs_err
      real,allocatable,dimension(:)    :: amf_trop_col
      real,allocatable,dimension(:)    :: amf_total_col
      real,allocatable,dimension(:)    :: no2_slnt_col
      real,allocatable,dimension(:)    :: no2_slnt_col_err
      real,allocatable,dimension(:)    :: o3_slnt_col
      real,allocatable,dimension(:)    :: o3_slnt_col_err
      real,allocatable,dimension(:,:)  :: avgk_obs
      real,allocatable,dimension(:,:)  :: prs_obs
   end type observation_data
!
   type(observation_data),allocatable,dimension(:,:) :: gome2a_data
   real                                :: dist
   real,allocatable,dimension(:,:)     :: dist_min
   real,allocatable,dimension(:,:)     :: icnt_min
!
   namelist /create_gome2a_obs_nml/filedir,filename,fileout,year,month,day,hour, &
   bin_beg_sec,bin_end_sec,fac_obs_error,use_log_no2, &
   lon_min,lon_max,lat_min,lat_max,path_model,file_model,nx_model, &
   ny_model,nz_model,obs_no2_reten_freq
!
! Record the current time, date, etc. to the logfile
   call initialize_utilities(source)
   call register_module(source,revision,revdate)
!
   call find_namelist_in_file("input.nml", "create_gome2a_obs_nml", iunit)
   read(iunit, nml = create_gome2a_obs_nml, iostat = io)
   call check_namelist_read(iunit, io, "create_gome2a_obs_nml")
!
! Read GOME2A NO2
   allocate(gome2a_data(nx_model,ny_model))
   allocate(dist_min(nx_model,ny_model))
   allocate(icnt_min(nx_model,ny_model))
!
   gome2a_data(:,:)%icnt=0
   fileid=100
   write(6,*)'opening ',TRIM(TRIM(filedir)//TRIM(filename))
   open(unit=fileid,file=TRIM(TRIM(filedir)//TRIM(filename)), &
   form='formatted', status='old', iostat=ios)
   read(fileid,*,iostat=ios) data_type, obs_id, &
   i_min, j_min
   gome2a_data(i_min,j_min)%data_type=data_type
   gome2a_data(i_min,j_min)%obs_id=obs_id
   do while (ios == 0)
      gome2a_data(i_min,j_min)%icnt=gome2a_data(i_min,j_min)%icnt+1
      icnt=gome2a_data(i_min,j_min)%icnt
      if(gome2a_data(i_min,j_min)%icnt.eq.1) then
         allocate(gome2a_data(i_min,j_min)%yr_obs(max_num_obs))
         allocate(gome2a_data(i_min,j_min)%mn_obs(max_num_obs))
         allocate(gome2a_data(i_min,j_min)%dy_obs(max_num_obs))
         allocate(gome2a_data(i_min,j_min)%hh_obs(max_num_obs))
         allocate(gome2a_data(i_min,j_min)%mm_obs(max_num_obs))
         allocate(gome2a_data(i_min,j_min)%ss_obs(max_num_obs))
         allocate(gome2a_data(i_min,j_min)%lat_obs(max_num_obs))
         allocate(gome2a_data(i_min,j_min)%lon_obs(max_num_obs))
         allocate(gome2a_data(i_min,j_min)%nlay_obs(max_num_obs))
         allocate(gome2a_data(i_min,j_min)%nlev_obs(max_num_obs))
         allocate(gome2a_data(i_min,j_min)%trop_index(max_num_obs))
         allocate(gome2a_data(i_min,j_min)%no2_trop_col_obs(max_num_obs))
         allocate(gome2a_data(i_min,j_min)%no2_trop_col_obs_err(max_num_obs))
         allocate(gome2a_data(i_min,j_min)%no2_total_col_obs(max_num_obs))
         allocate(gome2a_data(i_min,j_min)%no2_total_col_obs_err(max_num_obs))
         allocate(gome2a_data(i_min,j_min)%amf_trop_col(max_num_obs))
         allocate(gome2a_data(i_min,j_min)%amf_total_col(max_num_obs))
         allocate(gome2a_data(i_min,j_min)%no2_slnt_col(max_num_obs))
         allocate(gome2a_data(i_min,j_min)%no2_slnt_col_err(max_num_obs))
         allocate(gome2a_data(i_min,j_min)%o3_slnt_col(max_num_obs))
         allocate(gome2a_data(i_min,j_min)%o3_slnt_col_err(max_num_obs))
      endif
!
      read(fileid,*,iostat=ios) &
         gome2a_data(i_min,j_min)%yr_obs(icnt), &
         gome2a_data(i_min,j_min)%mn_obs(icnt), &
         gome2a_data(i_min,j_min)%dy_obs(icnt), &
         gome2a_data(i_min,j_min)%hh_obs(icnt), &
         gome2a_data(i_min,j_min)%mm_obs(icnt), &
         gome2a_data(i_min,j_min)%ss_obs(icnt)
      read(fileid,*,iostat=ios) &
         gome2a_data(i_min,j_min)%lat_obs(icnt), &
         gome2a_data(i_min,j_min)%lon_obs(icnt)
      read(fileid,*,iostat=ios) &
         gome2a_data(i_min,j_min)%nlay_obs(icnt), &
         gome2a_data(i_min,j_min)%nlev_obs(icnt)
         nlay_obs=gome2a_data(i_min,j_min)%nlay_obs(icnt)
         nlev_obs=gome2a_data(i_min,j_min)%nlev_obs(icnt)
      read(fileid,*,iostat=ios) &
         gome2a_data(i_min,j_min)%trop_index(icnt)
      if(gome2a_data(i_min,j_min)%icnt.eq.1) then
         allocate(gome2a_data(i_min,j_min)%avgk_obs(max_num_obs,nlay_obs))
         allocate(gome2a_data(i_min,j_min)%prs_obs(max_num_obs,nlev_obs))
      endif
      read(fileid,*,iostat=ios) &
         gome2a_data(i_min,j_min)%prs_obs(icnt,1:nlev_obs)
      read(fileid,*,iostat=ios) &
         gome2a_data(i_min,j_min)%avgk_obs(icnt,1:nlay_obs)
      read(fileid,*,iostat=ios) &
         gome2a_data(i_min,j_min)%no2_trop_col_obs(icnt), &
         gome2a_data(i_min,j_min)%no2_trop_col_obs_err(icnt)
      read(fileid,*,iostat=ios) &
         gome2a_data(i_min,j_min)%no2_total_col_obs(icnt), &
         gome2a_data(i_min,j_min)%no2_total_col_obs_err(icnt)
      read(fileid,*,iostat=ios) &
         gome2a_data(i_min,j_min)%amf_trop_col(icnt)
      read(fileid,*,iostat=ios) &
         gome2a_data(i_min,j_min)%amf_total_col(icnt)
      read(fileid,*,iostat=ios) &
         gome2a_data(i_min,j_min)%no2_slnt_col(icnt), &
         gome2a_data(i_min,j_min)%no2_slnt_col_err(icnt)
      read(fileid,*,iostat=ios) &
         gome2a_data(i_min,j_min)%o3_slnt_col(icnt), &
         gome2a_data(i_min,j_min)%o3_slnt_col_err(icnt)
      read(fileid,*,iostat=ios) data_type, obs_id, &
         i_min, j_min
      gome2a_data(i_min,j_min)%data_type=data_type
      gome2a_data(i_min,j_min)%obs_id=obs_id
   enddo
!
! Find the median
   gome2a_data(:,:)%median_col=0.   
   do i=1,nx_model
      do j=1,ny_model
         if(gome2a_data(i,j)%icnt.gt.1) then
            call median_code(gome2a_data(i,j)%median_col,gome2a_data(i,j)%no2_trop_col_obs(1:gome2a_data(i,j)%icnt), &
            gome2a_data(i,j)%icnt)
         endif
      enddo
   enddo
!
! Find the mode
   gome2a_data(:,:)%mode_col=0.   
   do i=1,nx_model
      do j=1,ny_model
         if(gome2a_data(i,j)%icnt.gt.20) then
            call mode_code(gome2a_data(i,j)%mode_col,gome2a_data(i,j)%no2_trop_col_obs(1:gome2a_data(i,j)%icnt), &
            gome2a_data(i,j)%icnt)
         endif
      enddo
   enddo
!
! Calculate mean
   do i=1,nx_model
      do j=1,ny_model
         if(gome2a_data(i,j)%icnt.gt.0) then
            gome2a_data(i,j)%mean_col=0.
            gome2a_data(i,j)%mean_err=0.
            do icnt=1,gome2a_data(i,j)%icnt
               gome2a_data(i,j)%mean_col=gome2a_data(i,j)%mean_col+gome2a_data(i,j)%no2_trop_col_obs(icnt)
               gome2a_data(i,j)%mean_err=gome2a_data(i,j)%mean_err+gome2a_data(i,j)%no2_trop_col_obs_err(icnt)
            enddo
            gome2a_data(i,j)%mean_col=gome2a_data(i,j)%mean_col/float(gome2a_data(i,j)%icnt)
            gome2a_data(i,j)%mean_err=gome2a_data(i,j)%mean_err/float(gome2a_data(i,j)%icnt)
         endif
      enddo
   enddo
!
! Calculate standard deviation   
   do i=1,nx_model
      do j=1,ny_model
         if(gome2a_data(i,j)%icnt.gt.1) then
            gome2a_data(i,j)%stdv_col=0.
            gome2a_data(i,j)%stdv_err=0.
            do icnt=1,gome2a_data(i,j)%icnt
               gome2a_data(i,j)%stdv_col=gome2a_data(i,j)%stdv_col+(gome2a_data(i,j)%no2_trop_col_obs(icnt)- &
               gome2a_data(i,j)%mean_col)**2.
               gome2a_data(i,j)%stdv_err=gome2a_data(i,j)%stdv_err+(gome2a_data(i,j)%no2_trop_col_obs_err(icnt)- &
               gome2a_data(i,j)%mean_err)**2.
            enddo
            gome2a_data(i,j)%stdv_col=sqrt(gome2a_data(i,j)%stdv_col/float(gome2a_data(i,j)%icnt-1))
            gome2a_data(i,j)%stdv_err=sqrt(gome2a_data(i,j)%stdv_err/float(gome2a_data(i,j)%icnt-1))
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
         if(gome2a_data(i,j)%icnt.eq.1) then
            icnt_min(:,:)=1
         else if(gome2a_data(i,j)%icnt.gt.1) then
            do icnt=1,gome2a_data(i,j)%icnt
               dist=abs(gome2a_data(i,j)%mean_col-gome2a_data(i,j)%no2_trop_col_obs(icnt))
               if(dist.lt.dist_min(i,j)) then
                  dist_min(i,j)=dist
                  icnt_min(i,j)=icnt
               endif
            enddo
         endif
!
! Write data to file
         if(gome2a_data(i,j)%icnt.gt.0) then
            icnt=icnt_min(i,j)
            write(fileid,*,iostat=ios) &
               trim(gome2a_data(i,j)%data_type), &
               gome2a_data(i,j)%obs_id, i, j
            write(fileid,*,iostat=ios) &
               gome2a_data(i,j)%yr_obs(icnt), &
               gome2a_data(i,j)%mn_obs(icnt), &
               gome2a_data(i,j)%dy_obs(icnt), &
               gome2a_data(i,j)%hh_obs(icnt), &
               gome2a_data(i,j)%mm_obs(icnt), &
               gome2a_data(i,j)%ss_obs(icnt)
            write(fileid,*,iostat=ios) &
               gome2a_data(i,j)%lat_obs(icnt), &
               gome2a_data(i,j)%lon_obs(icnt)
            write(fileid,*,iostat=ios) &
               gome2a_data(i,j)%nlay_obs(icnt), &
               gome2a_data(i,j)%nlev_obs(icnt)
            write(fileid,*,iostat=ios) &
               gome2a_data(i,j)%trop_index(icnt)
            write(fileid,*,iostat=ios) &
               gome2a_data(i,j)%prs_obs(icnt,1:nlev_obs)
            write(fileid,*,iostat=ios) &
               gome2a_data(i,j)%avgk_obs(icnt,1:nlay_obs)
            write(fileid,*,iostat=ios) &
               gome2a_data(i,j)%no2_trop_col_obs(icnt), &
               gome2a_data(i,j)%no2_trop_col_obs_err(icnt)
            write(fileid,*,iostat=ios) &
               gome2a_data(i,j)%no2_total_col_obs(icnt), &
               gome2a_data(i,j)%no2_total_col_obs_err(icnt)
            write(fileid,*,iostat=ios) &
               gome2a_data(i,j)%amf_trop_col(icnt)
            write(fileid,*,iostat=ios) &
               gome2a_data(i,j)%amf_total_col(icnt)
            write(fileid,*,iostat=ios) &
               gome2a_data(i,j)%no2_slnt_col(icnt), &
               gome2a_data(i,j)%no2_slnt_col_err(icnt)
            write(fileid,*,iostat=ios) &
               gome2a_data(i,j)%o3_slnt_col(icnt), &
               gome2a_data(i,j)%o3_slnt_col_err(icnt)
         endif
      enddo
   enddo
!               
! Deallocate structure arrays
   do i=1,nx_model
      do j=1,ny_model
         if(gome2a_data(i,j)%icnt.gt.0) then
            deallocate(gome2a_data(i,j)%yr_obs)
            deallocate(gome2a_data(i,j)%mn_obs)
            deallocate(gome2a_data(i,j)%dy_obs)
            deallocate(gome2a_data(i,j)%hh_obs)
            deallocate(gome2a_data(i,j)%mm_obs)
            deallocate(gome2a_data(i,j)%ss_obs)
            deallocate(gome2a_data(i,j)%lat_obs)
            deallocate(gome2a_data(i,j)%lon_obs)
            deallocate(gome2a_data(i,j)%nlay_obs)
            deallocate(gome2a_data(i,j)%nlev_obs)
            deallocate(gome2a_data(i,j)%trop_index)
            deallocate(gome2a_data(i,j)%avgk_obs)
            deallocate(gome2a_data(i,j)%prs_obs)
            deallocate(gome2a_data(i,j)%no2_trop_col_obs)
            deallocate(gome2a_data(i,j)%no2_trop_col_obs_err)
            deallocate(gome2a_data(i,j)%no2_total_col_obs)
            deallocate(gome2a_data(i,j)%no2_total_col_obs_err)
            deallocate(gome2a_data(i,j)%amf_trop_col)
            deallocate(gome2a_data(i,j)%amf_total_col)
            deallocate(gome2a_data(i,j)%no2_slnt_col)
            deallocate(gome2a_data(i,j)%no2_slnt_col_err)
            deallocate(gome2a_data(i,j)%o3_slnt_col)
            deallocate(gome2a_data(i,j)%o3_slnt_col_err)
         endif
      enddo
   enddo
   deallocate (gome2a_data)
   deallocate (dist_min)
   deallocate (icnt_min)   
end program gome2a_no2_thinner
!
subroutine median_code(fld_median,fld,nfld)
   implicit none
   integer                 :: ipt,iptt,nfld
   real                    :: fld_median,temp
   real, dimension(nfld)   :: fld
   real, dimension(1000)   :: temp_arr
!
   temp_arr(1:nfld)=fld(1:nfld)
   do ipt=2,nfld
      if(temp_arr(ipt) .lt. temp_arr(ipt-1)) then
         iptt=ipt
         do while (temp_arr(iptt) .lt. temp_arr(iptt-1))
            temp=temp_arr(iptt-1)
            temp_arr(iptt-1)=temp_arr(iptt)
            temp_arr(iptt)=temp
            iptt=iptt-1
            if(iptt .eq. 1) exit
         enddo
      endif
   enddo
   if(nfld/2*2 .eq.nfld) then
      fld_median=(temp_arr(nfld/2)+temp_arr(nfld/2+1))/2.
   else
      fld_median=temp_arr(nfld/2+1)
   endif
   return
end subroutine median_code
!
subroutine mode_code(fld_mode,fld,nfld)
   implicit none
   integer                  :: ipt,iptt,nfld
   integer                  :: cnt_max,cnt_idx
   integer, dimension(nfld) :: cnt
   real                     :: fld_mode,temp
   real                     :: bin_tol
   real, dimension(nfld)    :: fld
   real, dimension(1000)    :: temp_arr
!
   bin_tol=.1   
   temp_arr(1:nfld)=fld(1:nfld)
   do ipt=2,nfld
      if(temp_arr(ipt) .lt. temp_arr(ipt-1)) then
         iptt=ipt
         do while (temp_arr(iptt) .lt. temp_arr(iptt-1))
            temp=temp_arr(iptt-1)
            temp_arr(iptt-1)=temp_arr(iptt)
            temp_arr(iptt)=temp
            iptt=iptt-1
            if(iptt .eq. 1) exit
         enddo
      endif
   enddo
!
   cnt(:)=1
   do ipt=1,nfld
      do iptt=1,nfld
         if(ipt.eq.iptt) cycle
         if(temp_arr(ipt).eq.temp_arr(iptt)) cnt(ipt)=cnt(ipt)+1
      enddo
   enddo   
   cnt_max=cnt(1)
   cnt_idx=1
   do ipt=2,nfld
      if(cnt(ipt).gt.cnt_max) then
         cnt_max=cnt(ipt)
         cnt_idx=ipt
      endif
   enddo
!
   fld_mode=0.    
   if(cnt_max.gt.1) then
      fld_mode=temp_arr(cnt_idx)
   endif
   return
end subroutine mode_code
          
     
