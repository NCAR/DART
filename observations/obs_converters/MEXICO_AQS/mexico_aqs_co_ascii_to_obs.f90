
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

program mexico_aqs_co_ascii_to_obs

!=============================================
! MEXICO AQS SURFACE AQ obs
! Based from create_obs_sequence.f90
!=============================================
!
      use    utilities_mod, only :    timestamp,                 &
                                      register_module,           &
                                      initialize_utilities,      &
                                      find_namelist_in_file,     &
                                      check_namelist_read,       &
                                      error_handler,             &
                                      E_ERR,                     & 
                                      E_WARN,                    & 
                                      E_MSG,                     &
                                      E_DBG
!
      use obs_sequence_mod, only :    obs_sequence_type,         &
                                      write_obs_seq,             &
                                      static_init_obs_sequence,  &
                                      init_obs_sequence,         &
                                      init_obs,                  &
                                      set_obs_values,            &
                                      set_obs_def,               &
                                      set_qc,                    &
                                      set_qc_meta_data,          &
                                      set_copy_meta_data,        &
                                      insert_obs_in_seq,         &
                                      obs_type
!
      use obs_def_mod, only :         obs_def_type,              &
                                      set_obs_def_time,          &
                                      set_obs_def_type_of_obs,   &
                                      set_obs_def_error_variance,&
                                      set_obs_def_location,      &
                                      set_obs_def_key
!
      use assim_model_mod, only :     static_init_assim_model
!
      use location_mod, only :        location_type,             &
                                      set_location
!
      use time_manager_mod, only :    set_date,                  &
                                      set_calendar_type,         &
                                      time_type,                 &
                                      get_time
!
      use obs_kind_mod, only :        MEXICO_AQS_CO,                 &
                                      get_type_of_obs_from_menu
!
      use random_seq_mod, only :      random_seq_type,           &
                                      init_random_seq,           &
                                      random_uniform
!
      use sort_mod, only :            index_sort

      implicit none

! version controlled file description for error handling, do not edit
      character(len=*), parameter :: source   = 'mexico_aqs_co_ascii_to_obs.f90'
      character(len=*), parameter :: revision = ''
      character(len=*), parameter :: revdate  = ''
!
      type(obs_sequence_type)      :: seq
      type(obs_type)               :: obs
      type(obs_type)               :: obs_old
      type(obs_def_type)           :: obs_def
      type(location_type)          :: obs_location
      type(time_type)              :: obs_time
      integer,parameter            :: max_num_obs=2000000
      integer,parameter            :: indx_max=2000000
      integer,parameter            :: num_copies=1
      integer,parameter            :: num_qc=1
      integer,parameter            :: fileid=88
      integer,parameter            :: maxkinds=20
      integer                      :: icopy,iunit,io,indx,ierr,nndx
      integer                      :: year0, month0, day0, hour0, min0, sec0
      integer                      :: year1, month1, day1, hour1, min1, sec1
      integer                      :: calendar_type
      integer                      :: qc_count,qstatus
      integer                      :: seconds, days, which_vert
      integer                      :: obs_kind, obs_key
      integer                      :: days_in_month(12) =(/ &
                                      31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31  /)
      integer                      :: beg_year,beg_mon,beg_day,beg_hour,beg_min,beg_sec 
      integer                      :: end_year,end_mon,end_day,end_hour,end_min,end_sec
      integer                      :: anal_greg_sec,beg_greg_sec,end_greg_sec
      integer                      :: save_greg_sec,calc_greg_sec
      integer                      :: data_greg_sec_temp
      integer,dimension(indx_max)  :: year_temp,month_temp,day_temp,hour_temp
      real                         :: lat_mn,lat_mx,lon_mn,lon_mx
      real*8                       :: latitude,longitude,level
      real*8                       :: ob_err_var
      real                         :: lat_temp,lon_temp
      real                         :: lat,lon,obs_val,obs_err
      real*8,dimension(num_qc)     :: obs_qc
      real*8,dimension(num_copies) :: obs_val_out
      real                         :: co_log_max, co_log_min
      character(len=2)             :: chr_month, chr_day
      character(len=4)             :: chr_year
      character(len=20)            :: dmy
      character(len=20)            :: var1,var2,var3,var4
      character(len=20),dimension(:),allocatable :: var5
      character(len=129)           :: copy_meta_data
      character(len=129)           :: qc_meta_data='MEXICO_AQS QC index'
      character(len=129)           :: file_name='mexico_aqs_obs_seq'
      character(len=180)           :: file_in_stations,file_in_data
      logical                      :: use_log_co,use_log_o3,use_log_nox,use_log_so2,use_log_pm10,use_log_pm25

      real                         :: pi,fac,fac_err,fac_obs_error,ran1,ran2,zfac
!
      real                         :: correction_old      
      integer                      :: iobs,nobs
      character(len=150),dimension(maxkinds) :: obs_list
      character(len=150)           :: correction_filename,path_filein
      logical                      :: does_file_exist
!      
      integer                                    :: iflg,ipt,ista,istaa,num_stations
      real,dimension(:),allocatable              :: sta_lon,sta_lat
      character(len=10),dimension(:),allocatable :: sta_ids,sta_ids_data
      character(len=200)                         :: record
      character(len=20),dimension(:),allocatable :: record_array                      
      real,dimension(:,:),allocatable :: obs_val_temp
      integer                                    :: i,m
!
      namelist /create_mexico_aqs_obs_nml/beg_year,beg_mon,beg_day, &
      beg_hour,beg_min,beg_sec,end_year,end_mon,end_day,end_hour,end_min,end_sec, &
      fac_obs_error,file_in_stations,file_in_data,lat_mn,lat_mx,lon_mn,lon_mx,use_log_co,use_log_o3,use_log_nox, &
      use_log_so2,use_log_pm10,use_log_pm25
!
      namelist /bias_correct_nml/path_filein,does_file_exist,correction_filename,nobs,obs_list

!============================================================
!obs sequence extra variables
!============================================================

      pi=4.*atan(1.0)
      fac=1.0
      fac_err=0.3
      fac_obs_error=1.0
      obs_qc(1)=0.

      save_greg_sec=-9999                                                 
! Record the current time, date, etc. to the logfile                                       
      call initialize_utilities(source)
      call register_module(source,revision,revdate)
!
! Initialize the obs_sequence module
      call static_init_obs_sequence()
!
! Initialize an obs_sequence structure
      call init_obs_sequence(seq, num_copies, num_qc, max_num_obs)
      calendar_type=3
      call set_calendar_type(calendar_type)
!
! Initialize the obs variable
      call init_obs(obs, num_copies, num_qc)
!
      do icopy =1, num_copies
         if (icopy == 1) then
             copy_meta_data='MEXICO_AQS observation'
         else
             copy_meta_data='Truth'
         endif
         call set_copy_meta_data(seq, icopy, copy_meta_data)
      enddo
      call set_qc_meta_data(seq, 1, qc_meta_data)
!
! READ MEXICO AQS OBSERVATIONS
! Read the namelist entry
      call find_namelist_in_file("create_mexico_aqs_obs_nml.nl", "create_mexico_aqs_obs_nml", iunit)
      read(iunit, nml = create_mexico_aqs_obs_nml, iostat = io)
      call check_namelist_read(iunit, io, "create_mexico_aqs_obs_nml")
!
      min0=0
      sec0=0
      print *, 'beg_year        ',beg_year
      print *, 'beg_mon         ',beg_mon
      print *, 'beg_day         ',beg_day
      print *, 'beg_hour        ',beg_hour
      print *, 'beg_min         ',beg_min
      print *, 'beg_sec         ',beg_sec
      print *, 'end_year        ',end_year
      print *, 'end_mon         ',end_mon
      print *, 'end_day         ',end_day
      print *, 'end_min         ',end_min
      print *, 'end_sec         ',end_sec
      print *, 'file_in_stations ',trim(file_in_stations)
      print *, 'file_in_data ',trim(file_in_data)
      print *, 'lat_mn          ',lat_mn
      print *, 'lat_mx          ',lat_mx
      print *, 'lon_mn          ',lon_mn
      print *, 'lon_mx          ',lon_mx
      print *, ' '
!
! Read the namelist entry
      call find_namelist_in_file("bias_correct_nml", "bias_correct_nml", iunit)
      read(iunit, nml = bias_correct_nml, iostat = io)
      call check_namelist_read(iunit, io, "bias_correct_nml")
      print *, 'path_filein ',trim(path_filein)
      print *, 'does_file_exit ',does_file_exist
      print *, 'correction_filename ',trim(correction_filename)
      print *, 'nobs ',nobs
      do iobs=1,nobs
         print *, 'obs_list ',iobs,trim(obs_list(iobs))
      enddo
!
! Determine bias correction
      correction_old=0.
      if(does_file_exist) then
         iunit=101
         open(unit=iunit,file=trim(correction_filename),form='unformatted', &
         status='old',action='READ')
         rewind(iunit)
         correction_old=-999.
         do iobs=1,nobs
            read(iunit) correction_old
            if(trim(obs_list(iobs)).eq.'MEXICO_AQS_CO') then
               exit
            endif
         enddo
         if(correction_old.eq.-999.) then
            print *, 'APM: MEXICO_AQS_CO - Error assigning bias correction'
            stop
         endif
         close(iunit)
      endif
!
! Fix leap year number of days
      days_in_month(2)=28
      if(beg_year/4*4.eq.0) days_in_month(2)=29
      beg_greg_sec=calc_greg_sec(beg_year,beg_mon,beg_day,beg_hour,beg_min,beg_sec,days_in_month)
      days_in_month(2)=28
      if(end_year/4*4.eq.0) days_in_month(2)=29
      end_greg_sec=calc_greg_sec(end_year,end_mon,end_day,end_hour,end_min,end_sec,days_in_month)
!
! Read the station id file
      iunit=101
      open(unit=iunit,file=trim(file_in_stations))
      read(iunit,*) var1
      read(var1(1:2),*) num_stations
!      print *, 'Number of Stations ',num_stations
!
      allocate(obs_val_temp(indx_max,num_stations))
      allocate(sta_ids_data(num_stations))
      allocate(sta_ids(num_stations))
      allocate(sta_lon(num_stations))
      allocate(sta_lat(num_stations))
      do ista=1,num_stations
         read(iunit,'(a3, 1x, f10.6, 1x, f10.6)') sta_ids(ista),sta_lon(ista),sta_lat(ista)
!         print *, 'Station: ID, LON, LAT ',trim(sta_ids(ista)),sta_lon(ista),sta_lat(ista)
      enddo
      close (iunit)
!
! Read data file
      iunit=101
      open(unit=iunit,file=trim(file_in_data))
! Need to adjust the format for the number of stations (assume 32
      read(iunit,'(10x,32(1x,a3))') sta_ids_data(1:num_stations)
!      do ista=1,num_stations
!         print *, 'APM: Data station ids ',ista,trim(sta_ids_data(ista))
!      enddo
!
! Read the data records
      indx=0
      do
         indx=indx+1
         read(iunit,'(a200)',iostat=ierr) record
         nndx=count(transfer(record,'a',len(record))==",")
         nndx=nndx+1
         allocate(record_array(nndx))
         m=1
         do i=1,nndx-1
            ipt=index(record(m:),',')
            record_array(i)=adjustl(record(m:m+ipt-2))
            m=m+ipt
         enddo
         record_array(nndx)=adjustl(record(m:))
!
         if(len_trim(record_array(2)).eq.1) then         
            read(record_array(2)(1:1),'(i1)') hour_temp(indx)
         else
            read(record_array(2)(1:2),'(i2)') hour_temp(indx)
         endif            
         read(record_array(1)(1:2),'(i2)') day_temp(indx)
         read(record_array(1)(4:5),'(i2)') month_temp(indx)
         read(record_array(1)(7:10),'(i4)') year_temp(indx)
         do ista=1,32
            if(record_array(ista+2).eq.'-99') then
               read(record_array(ista+2),'(f4.0)') obs_val_temp(indx,ista)
            else
               read(record_array(ista+2),'(f5.2)') obs_val_temp(indx,ista)
            endif
         enddo
         deallocate(record_array)
         if(ierr.lt.0) exit
      enddo
      if(ierr.gt.0) then
         print *,'APM: Read error '
         stop
      endif
      if(ierr.le.0) then
         print *,'APM: EOF reached '
      endif
      close(iunit)
      nndx=indx
!
! Put data in obs_seq_file
      qc_count=0
      do indx=1,nndx
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Need to fix this to change to next day.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
         if(hour_temp(indx).eq.24) cycle
!         
         data_greg_sec_temp=calc_greg_sec(year_temp(indx),month_temp(indx),day_temp(indx), &
         hour_temp(indx),0,0,days_in_month)
         do ista=1,num_stations
            if (obs_val_temp(indx,ista) .le. 0.) cycle
            if (beg_greg_sec .gt. data_greg_sec_temp .or. end_greg_sec .lt. data_greg_sec_temp) cycle
!
! Find lon and lat for this station
            iflg=0
            do istaa=1,num_stations
               if(trim(sta_ids_data(ista)).eq.trim(sta_ids(istaa))) then
                  lon_temp=sta_lon(istaa)
                  lat_temp=sta_lat(istaa)
                  iflg=1
                  exit
               endif
            enddo
            if(iflg.eq.0) then
               print *, 'APM: ERROR - MEXICO AQS CO failed to find lon and lat for station '
               stop
            endif
            if(lat_mn .gt. lat_temp .or. lat_mx .lt. lat_temp .or. &
            lon_mn .gt. lon_temp .or. lon_mx .lt. lon_temp) cycle
!
! Keep this observation and place in obs_sequence file               
            qstatus=0
            qc_count=qc_count+1
            if(data_greg_sec_temp.gt.save_greg_sec) save_greg_sec=data_greg_sec_temp
!
! location
            level=1
            latitude=lat_temp 
            if (lon_temp<0) then
               longitude=lon_temp+360
            else
               longitude=lon_temp
            endif
!
! time 
            year1=year_temp(indx)
            month1=month_temp(indx)
            day1=day_temp(indx)
            hour1=hour_temp(indx)
            min1=0
            sec1=0
!
! correction_old if for bias correction if used           
            obs_val_out(1:num_copies)=obs_val_temp(indx,ista)*fac - correction_old
            obs_val=obs_val_temp(indx,ista)*fac - correction_old
            obs_err=(obs_val_temp(indx,ista)*fac - correction_old)*fac_err*fac_obs_error
            ob_err_var = obs_err*obs_err
!
!            print *, year1,month1,day1,hour1,min1
            obs_time=set_date(year1,month1,day1,hour1,min1,sec1)
            call get_time(obs_time, seconds, days)
!        
!            which_vert           = -1       ! vert is surface
            which_vert           = 1       ! vert is level
            obs_location         = set_location(longitude, latitude, level, which_vert)
            obs_kind             = MEXICO_AQS_CO
!
            call set_obs_def_type_of_obs(obs_def, obs_kind)
            call set_obs_def_location(obs_def, obs_location)
            call set_obs_def_time(obs_def, obs_time)
            call set_obs_def_error_variance(obs_def, ob_err_var)
            call set_obs_def_key(obs_def, qc_count)
            call set_obs_values(obs, obs_val_out, 1)
            call set_qc(obs, obs_qc, num_qc)
            call set_obs_def(obs, obs_def)
            if ( qc_count == 1 .or. data_greg_sec_temp.lt.save_greg_sec) then
               call insert_obs_in_seq(seq, obs)
            else
               call insert_obs_in_seq(seq, obs, obs_old )
            endif
            obs_old=obs
         enddo
      enddo
!
!----------------------------------------------------------------------
! Write the sequence to a file
!----------------------------------------------------------------------
      call write_obs_seq(seq, file_name)
      close(fileid)
      deallocate(obs_val_temp)
      deallocate(sta_ids_data)
      deallocate(sta_ids)
      deallocate(sta_lon)
      deallocate(sta_lat)

!-----------------------------------------------------------------------------
! Clean up
!-----------------------------------------------------------------------------
      call timestamp(string1=source,string2=revision,string3=revdate,pos='end')
   end program mexico_aqs_co_ascii_to_obs
!
   integer function calc_greg_sec(year,month,day,hour,minute,sec,days_in_month)
      implicit none
      integer                  :: i,j,k,year,month,day,hour,minute,sec
      integer, dimension(12)   :: days_in_month
!
! assume time goes from 00:00:00 to 23:59:59  
      calc_greg_sec=0
      do i=1,month-1
         calc_greg_sec=calc_greg_sec+days_in_month(i)*24*60*60
      enddo
      do i=1,day-1
         calc_greg_sec=calc_greg_sec+24*60*60
      enddo
      do i=1,hour
         calc_greg_sec=calc_greg_sec+60*60
      enddo
      do i=1,minute
         calc_greg_sec=calc_greg_sec+60
      enddo
      calc_greg_sec=calc_greg_sec+sec
   end function calc_greg_sec

