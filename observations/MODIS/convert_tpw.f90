! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_tpw

! convert MODIS observations of Total Precipitable Water into
! DART observation sequence files.

use        types_mod, only : r8, metadatalength
use time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, set_time,&
                             increment_time, get_time, set_date, operator(-),  &
                             print_date
use    utilities_mod, only : initialize_utilities, find_namelist_in_file,    &
                             check_namelist_read, nmlfileunit, do_nml_file,   &
                             get_next_filename, error_handler, E_ERR, E_MSG, &
                             nc_check, find_textfile_dims, do_nml_term
use     location_mod, only : VERTISSURFACE, set_location
use obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq,       &
                             static_init_obs_sequence, init_obs, destroy_obs, &
                             write_obs_seq, init_obs_sequence, get_num_obs,   &
                             insert_obs_in_seq, destroy_obs_sequence,         &
                             set_copy_meta_data, set_qc_meta_data, set_qc,    &
                             set_obs_values, set_obs_def, insert_obs_in_seq
use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_kind, &
                             set_obs_def_error_variance, set_obs_def_location, &
                             set_obs_def_key
use     obs_kind_mod, only : MODIS_TOTAL_PRECIPITABLE_WATER

use           netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
 "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer, parameter ::   num_copies = 1,   &   ! number of copies in sequence
                        num_qc     = 1,   &   ! number of QC entries
                        num_new_obs= 5000     ! number of observations 

character(len = 129) :: ObsBase = '/home/hliu/data'
character (len=metadatalength) :: meta_data
character (len=129) :: msgstring, next_infile
character (len=80)  :: name
character (len=19)  :: datestr
character (len=6)   :: subset

type(obs_def_type)      :: obs_def
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs 
type(time_type)         :: time_obs

!--------------------------------------
      integer  :: i, ibin
      integer  :: year, month, day, hour, minute, tot_days, ierr, io
      integer  :: iyear, imonth, iday, ihour, imin, isec
      integer leng, nrec, nlev, idd, max_rec
      parameter (leng=9, max_rec= 900000)

      real(r8) lon, lat, tpw, time, obserr, seconds
      real(r8) prof(leng), prof0(leng), prof_ave(leng)
      real(r8) prof_all(max_rec, leng), prof_in(max_rec, leng), &
               lon_in(max_rec), lat_in(max_rec)
      real(r8) missing 
      real(r8) :: obs_val(1), qc_val(1), hght
      logical :: first_obs

      integer  iqc, ipc , itype, num, numq, j, k, nr, iuf
      real(r8)  lati, latk, loni, lonk
      character infile*100
      character(len = 8 ):: obsdate

      integer  :: nloc, navg, kk, nn, obs_num, oday, osec
      real(r8) :: nu, oerr, qc, dlat, dlon
      real(r8) :: bin_beg, bin_end
      integer  :: bin_start, bin_interval, bin_half_width
      character(len=128) :: out_file

      namelist /convert_tpw_nml/year, month, day, tot_days, bin_start, bin_interval, bin_half_width, ObsBase

! Read namelist ...

      open (10, file='input.nml')
      ierr = 1

      do while(ierr /= 0)
      read(10, nml = convert_tpw_nml, iostat = io, end = 11)
      enddo
   11 continue
      close(10)

      missing = -9999.9
!--------------------------------------
      do 3333 idd =day, day+ tot_days
!--------------------------------------

! initialize some values

call set_calendar_type(GREGORIAN)
call initialize_utilities()

call static_init_obs_sequence()
call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

  call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)

  do k = 1, num_copies
    meta_data = 'MODIS observation'
    call set_copy_meta_data(obs_seq, k, meta_data)
  end do
  do k = 1, num_qc
    meta_data = 'MODIS QC'
    call set_qc_meta_data(obs_seq, k, meta_data)
  end do


    write(obsdate, '(i4.4, i2.2, i2.2)') year, month, idd

!infile=trim(adjustl(ObsBase))//'/sinlaku/AQUA_MODIS_TPW_'//obsdate//'_Sinlaku.txt'
!out_file = 'tpw_AQUA_sinlaku_obs_seq.'//obsdate

 infile=trim(adjustl(ObsBase))//'/sinlaku/TERRA_MODIS_TPW_'//obsdate//'_Sinlaku.txt'
 out_file = 'tpw_TERRA_sinlaku_obs_seq.'//obsdate

        print*, 'input file= ', infile

        iuf = 90
        open(iuf,file=infile,form='formatted')

      prof_all(:,:) = missing
      nrec = 0

!   read in all observation records

 obsloop: do 
    read(iuf,888, end=200) lat, lon, tpw, iyear, imonth, iday, ihour, imin, seconds
 888 format(f11.6, f13.5, f10.4, 4x, i4,4i3, f7.3)
    nrec = nrec + 1
     if(lon < 0.0 ) lon = lon + 360.0
    prof(1) = lat
    prof(2) = lon
    prof(3) = tpw
    prof(4) = iyear
    prof(5) = imonth
    prof(6) = iday
    prof(7) = ihour
    prof(8) = imin
    prof(9) = seconds

    do k=1, leng
    prof_all(nrec, k) = prof(k)
    end do

 end do obsloop
 200  continue
      close(iuf)

      print*, 'nrec= ', nrec
      if(max_rec .le. nrec) then
      print*, 'max_rec too small'
      stop
      endif

 first_obs = .true.

! ---------------------------------------------------
!   select observations for specific time interval
! ---------------------------------------------------
      navg = 0
      obs_num = 0
      ibin = 0
      bin_beg = 0.0
      bin_end = 0.0

!    time bin loop
      binloop:  do while(bin_end < 25.0)

      nloc = 0

!  read in all of the observations in the selected bin info into arrays

      ibin = ibin + 1
      bin_beg = bin_start + (ibin-1)*bin_interval - bin_half_width
      bin_end = bin_start + (ibin-1)*bin_interval + bin_half_width
      print*, 'bin= ' , bin_beg, bin_end

      do 6000 nr =1, nrec

       do k=1, leng
       prof0(k) = prof_all(nr, k)
       end do

      time  = prof0(7) + prof0(8)/60.0 + prof0(9)/3600.0

      if(time .le. bin_beg .or. time .gt. bin_end) go to 600

       nloc = nloc + 1
!     print*, 'nloc=', nloc, lat, lon, time
       do k=1, leng
       prof_in(nloc, k) = prof0(k)
       lat_in(nloc)  = prof0(1)
       lon_in(nloc)  = prof0(2)
       end do

 600 continue
 6000 continue

       print*, 'nloc= ', nloc
!    do the average of the nearby observations
      do 900 kk = 1, nloc ! loop over all observations

      latk = prof_in(kk,1)
      lonk = prof_in(kk,2)

      if ( latk .ne. missing ) then
        nu = 1.0
        do k=1, leng
        prof_ave(k)  = prof_in(kk, k)
        end do
!       print*, 'count0= ', kk, latk, lonk
!
!   use the time of the first observation for the averaged one.
      iyear   = prof_ave(4)
      imonth  = prof_ave(5)
      iday    = prof_ave(6)
      ihour   = prof_ave(7)
      imin    = prof_ave(8)
      isec    = prof_ave(9)

   time_obs = set_date(iyear, imonth, iday, ihour, imin, isec)
   call get_time(time_obs,  osec, oday)

      do 80 nn = (kk+1), nloc

      lati = prof_in(nn,1)
      loni = prof_in(nn,2)
        if ( lati /= missing ) then

          dlat = abs(latk-lati)
          dlon = abs(lonk-loni)

!   select the nearby observations to average; lat/lon in degree.
!         if ( dlat < 0.5 .and. dlon < 0.5 ) then
          if ( dlat < 0.6 .and. dlon < 0.6 ) then
!          print*, 'count = ', nn, lati, loni
              nu = nu + 1.0
              do k=1, leng
              prof_ave(k)  = prof_ave(k) + prof_in(nn, k)
              end do 

            prof_in(nn, :) = missing     !  lat
          end if

        end if
      80 continue

      if ( nu > 0.0 ) then
        navg = navg + 1
        do k=1, leng
        prof_ave(k) = prof_ave(k) / nu
        end do

!   ----------------------
      lat   = prof_ave(1)
      lon   = prof_ave(2)
      tpw   = prof_ave(3)
      time  = prof_ave(7) + prof_ave(8)/60.0 + prof_ave(9)/3600.0
!       print*, 'avg= ', navg, nu, lat, lon, tpw

      if(lon < 0.0 ) lon = lon + 360.0

!   select the obs in ATLANTIC domain -----------------------------------
!       if(lat .lt.   5.0 .or. lat .gt. 35.0 )  go to 500
!       if(lon .lt. 250.0 .or. lon.gt. 330.0 )  go to 500

!    For W Pacific domain
        if(lat .lt.   5.0 .or. lat .gt. 45.0 )  go to 500
        if(lon .lt. 95.0 .or. lon.gt. 160.0 )  go to 500

     qc   = 0.0_r8 
     hght = 0.0_r8
     oerr = 0.5     !  cm
     call set_obs_def_location(obs_def,set_location(lon,lat,hght,VERTISSURFACE))
     call set_obs_def_kind(obs_def, MODIS_TOTAL_PRECIPITABLE_WATER)
     call set_obs_def_time(obs_def, set_time(osec, oday))
     call set_obs_def_error_variance(obs_def, oerr * oerr)
     call set_obs_def_key(obs_def, obs_num)
     call set_obs_def(obs, obs_def)

     obs_val(1) = tpw
     call set_obs_values(obs, obs_val)
     qc_val(1)  = qc
     call set_qc(obs, qc_val)

     if (first_obs) then
        call insert_obs_in_seq(obs_seq, obs)
        first_obs = .false.
     else
        call insert_obs_in_seq(obs_seq, obs, prev_obs)
     endif
     obs_num = obs_num+1
     prev_obs = obs

  500  continue
!----------------------------
      end if     ! nu> loop
      end if     ! if ( latk
 900 continue    ! nloc loop
! ------------------
      end do binloop
!-------------------

   call write_obs_seq(obs_seq, out_file)

   call destroy_obs(obs)

 3333 continue

      stop 
      end

! <next few lines under version control, do not edit>
! $URL$
! $Revision$
! $Date$

