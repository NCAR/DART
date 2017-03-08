! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

PROGRAM littler_tf_dart

use        types_mod, only : r8, DEG2RAD, RAD2DEG, MISSING_I, MISSING_R8
use    utilities_mod, only : open_file, close_file, file_exist, &
                             get_unit, initialize_utilities, &
                             register_module, error_handler, E_ERR, E_MSG, &
                             finalize_utilities, logfileunit
use obs_sequence_mod, only : obs_type, obs_sequence_type, init_obs_sequence, &
                             insert_obs_in_seq, write_obs_seq, read_obs_seq, &
                             set_qc, set_qc_meta_data, set_obs_values, set_copy_meta_data, &
                             assignment(=), get_obs_time_range, &
                             init_obs, static_init_obs_sequence, get_num_qc, set_obs_def, &
                             get_num_obs, get_max_num_obs, get_obs_values, &
                             get_time_range_keys, get_obs_from_key, get_obs_def, get_qc
use      obs_def_mod, only : copy_obs_def, obs_def_type, &
                             get_obs_def_time, get_obs_def_location, &
                             get_obs_def_error_variance, get_obs_def_type_of_obs, &
                             set_obs_def_type_of_obs, set_obs_def_location, set_obs_def_time, &
                             set_obs_def_error_variance
use     obs_kind_mod, only : RADIOSONDE_U_WIND_COMPONENT, &
                             RADIOSONDE_V_WIND_COMPONENT, &
                             RADIOSONDE_SURFACE_PRESSURE, &
                             RADIOSONDE_TEMPERATURE, &
                             RADIOSONDE_SPECIFIC_HUMIDITY, &
                             max_defined_types_of_obs, get_name_for_type_of_obs, map_type_of_obs_table
use     location_mod, only : location_type, get_location, query_location, set_location
use time_manager_mod, only : time_type, get_time, set_calendar_type, GREGORIAN, get_date, &
                             set_date, operator(/=), set_time

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type(obs_sequence_type) :: dart_seq
type(obs_type)          :: obs, prev_obs
type(obs_def_type)      :: obs_def
type(location_type)     :: location
type(time_type)         :: time, stime, ftime

INTEGER                 :: iunit, end_of_file

integer                 :: key_bounds(2), obs_kind, which_vert, pwhich_vert
real(r8)                :: obs_value(1), qc(1), erms
real(r8), dimension(3)  :: loc, sloc, floc
real(r8)                :: lon, lat, vloc, pvloc

integer           :: is, ie, iobs, k
integer           :: num_obs, num_copies, num_qc, max_num_obs, dart_seq_num_obs, num_obs_in_set
integer, allocatable :: keys(:)

character(len = 32)  :: obs_name, obs_no_support(max_defined_types_of_obs)
character(len = 129) :: msgstring

integer              :: n_no_support

character(len = 129) :: dart_file_name     = 'obs_seq.out', &
                        littler_file_name = 'little-r.dat'

integer              :: calendar_type      = GREGORIAN

character *20  :: date_char
CHARACTER *120 :: rpt_format
CHARACTER *120 :: meas_format
CHARACTER *120 :: end_format
character *40  :: tst_id, tst_name, tst_pltfrm, tst_src
real(r8)       :: tst_ter, tst_slp, tst_xlat, tst_xlon, tst_Psfc
integer        :: kx, iseq_num, tst_sut, tst_julian
logical        :: tst_sound, tst_bogus, tst_discard
integer        :: i1,i2,i3,i4,i5,i6,i7,i8,i9,i10
integer        :: i11,i12,i13,i14,i15,i16,i17,i18
real(r8)       :: f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12

integer, parameter :: n_wind_pres = 26, n_pres = 15
real(r8)           :: wind_error(n_wind_pres)
real(r8)           :: temp_error(n_pres), specific_humidity_error(n_pres)

! Observational errors are contained in a separate file 'obserr.txt'

CHARACTER (LEN=80) :: filein      != 'obserr.txt'
CHARACTER (LEN=80) :: platform    != 'RAOBS'
CHARACTER (LEN=80) :: sensor_name != 'WIND SENSOR ERRORS'

integer :: wind_pressures(n_wind_pres)
data wind_pressures / 10,   20,   30,   40,   50,  100,  150, &
                     200,  250,  300,  350,  400,  450,  500, &
                     550,  600,  650,  700,  750,  800,  850, &
                     900,  950, 1000, 1050, 1100/

integer :: pressure_levels(n_pres)
data pressure_levels / 1000,  850,  700,  500,  400, &
                        300,  250,  200,  150,  100, &
                         70,   50,   30,   20,   10/


real(r8), allocatable :: z(:),t(:),td(:),spd(:),dir(:), &
                         uu(:), vv(:), p(:), cld(:), ciel(:)
integer, allocatable  :: zpp_qc(:), tt_qc(:), td_qc(:), &
                         zuu_qc(:), zvv_qc(:), p_qc(:), &
                         spd_qc(:), dir_qc(:), cld_qc(:), ciel_qc(:)

logical :: littler_to_dart, out_of_range

!------------------------------------------------------------------------------

rpt_format =  ' ( 2f20.5 , 4a40 , f20.5 , 5i10 , 3L10 , ' &
                  // ' 2i10 , a20 ,  13( f13.5 , i7 ) ) '
meas_format =  ' ( 10( f13.5 , i7 ) ) '
end_format = ' ( 3 ( i7 ) ) '

call initialize_utilities('littler_tf_dart')
call register_module(source, revision, revdate)

call set_calendar_type(calendar_type)

write(*,*) 'littler to DART (.true./T) or DART to littler (.false./F)?'

read(*,*) littler_to_dart

tst_id = 'dart_id'
tst_name = 'dart_name'
tst_pltfrm = 'dart_pltfrm'
tst_src = 'dart_src'
tst_ter = MISSING_R8
tst_bogus = .false.
tst_discard = .false.
tst_slp = MISSING_R8
i1 = 0
i2 = 0
iseq_num = 0
i3 = 0
tst_sut = MISSING_I
tst_julian = MISSING_I
i6 = 0
i7 = 0
i8 = 0
i9 = 0
i10 = 0
i11 = 0
i12 = 0
i13 = 0
i14 = 0
i15 = 0
i16 = 0
i17 = 0
i18 = 0
f1 = MISSING_R8
f2 = MISSING_R8
f3 = MISSING_R8
f4 = MISSING_R8
f5 = MISSING_R8
f6 = MISSING_R8
f7 = MISSING_R8
f8 = MISSING_R8
f9 = MISSING_R8
f10 = MISSING_R8
f11 = MISSING_R8
f12 = MISSING_R8
wind_error(:) = -1.0_R8
temp_error(:) = -1.0_R8
specific_humidity_error(:) = -1.0_R8
n_no_support = 0

call static_init_obs_sequence()

num_copies = 1
num_qc = 1

! must use read_obs_seq_header

call init_obs(obs, num_copies, num_qc)

! -------------------------------------------------------------------
! Initialize the counter:

num_obs = 0

if(.not. littler_to_dart) then

   ! Conversion from DART to littler format.
   ! littler format needs soundings

   call read_obs_seq(dart_file_name, 0, 0, 0, dart_seq)

   dart_seq_num_obs = get_num_obs(dart_seq)
   stime = set_time(0, 0)
   ftime = set_time(0, 200000)
   call get_obs_time_range(dart_seq, stime, ftime, key_bounds, num_obs_in_set, &
                           out_of_range)
   if(num_obs_in_set /= dart_seq_num_obs) then
      write(msgstring, *) 'Did not get all obs. Got ',num_obs_in_set, &
           '. In file: ',dart_seq_num_obs
      call error_handler(E_ERR,'littler_tf_dart', &
           msgstring,source,revision,revdate)
   endif
   allocate(keys(dart_seq_num_obs))

   call get_time_range_keys(dart_seq, key_bounds, dart_seq_num_obs, keys)

   ! This should be get_first_obs(), get_last_obs()

   print*,'First obs key: ',keys(1)
   print*,'Last  obs key: ',keys(dart_seq_num_obs)
   print*,'Number of obs: ',dart_seq_num_obs

   iunit = get_unit()
   open(unit=iunit,file=littler_file_name,status='new')

   is = 1
   ie = 1
   do while (ie <= dart_seq_num_obs+1)

      call get_obs_from_key(dart_seq, keys(is), obs)
      call get_obs_def(obs, obs_def)

      sloc = get_location(get_obs_def_location(obs_def))
      stime = get_obs_def_time(obs_def)

      if(ie <= dart_seq_num_obs) then
         call get_obs_from_key(dart_seq, keys(ie), obs)
         call get_obs_def(obs, obs_def)
      endif

      floc = get_location(get_obs_def_location(obs_def))
      ftime = get_obs_def_time(obs_def)

      ! Write to littler file when a full sounding has been read.

      if ( ftime /= stime .or. sloc(1) /= floc(1) .or. &
           sloc(2) /= floc(2) .or. ie == dart_seq_num_obs+1) then

         tst_Psfc = MISSING_R8

         kx = 1
         call get_obs_from_key(dart_seq, keys(is), obs)
         call get_obs_def(obs, obs_def)
         location = get_obs_def_location(obs_def)
         loc = get_location(location)
         pwhich_vert = nint(query_location(location,'which_vert'))
         pvloc = loc(3)

         ! Count the number of pressure levels.

         do iobs = is+1, ie-1
            call get_obs_from_key(dart_seq, keys(iobs), obs)
            call get_obs_def(obs, obs_def)
            location = get_obs_def_location(obs_def)
            loc = get_location(location)
            which_vert = nint(query_location(location,'which_vert'))

!!$            if (which_vert == -1) then
!!$               obs_kind =  map_type_of_obs_table(get_obs_def_type_of_obs(obs_def))
!!$               if (obs_kind == RADIOSONDE_SURFACE_PRESSURE) then
!!$                  call get_obs_values(obs, obs_value, 1)
!!$                  tst_Psfc = obs_value(1)
!!$               endif
!!$            endif

            if (loc(3) /= pvloc .and. which_vert == pwhich_vert) kx = kx + 1
            pvloc = loc(3)
            pwhich_vert = which_vert
         enddo

!!$         if (kx > 1) then
         obs_kind = map_type_of_obs_table(get_obs_def_type_of_obs(obs_def))
         obs_name = get_name_for_type_of_obs(obs_kind)
         if ( obs_name(1:10) == 'RADIOSONDE' ) then
            tst_sound = .true.
            tst_pltfrm = 'FM-35 TEMP'
!!$         else
!!$            tst_sound = .false.
         endif

         call set_str_date(date_char, stime)

         ! In littler format, longitude is defined within [-180,180].

         tst_xlat = loc(2)
         tst_xlon = loc(1)
         if (tst_xlon > 180.0_r8) tst_xlon = tst_xlon - 360.0_r8

         ! Now the header can be written.

         WRITE ( UNIT = iunit , ERR = 19 , FMT = rpt_format ) &
              tst_xlat, tst_xlon, tst_id , tst_name, &
              tst_pltfrm, tst_src, tst_ter, kx, i1, i2, iseq_num, i3, &
              tst_sound, tst_bogus, tst_discard, tst_sut, tst_julian, date_char, &
              tst_slp, i6, f1, i7, f2, i8, f3, i9, tst_Psfc, i10, &
              f5, i11, f6, i12, f7, i13, f8, i14, f9, i15, &
              f10, i16, f11, i17, f12, i18

         allocate(p(kx),p_qc(kx), z(kx),zpp_qc(kx), t(kx),tt_qc(kx), td(kx),td_qc(kx), &
              spd(kx),spd_qc(kx), dir(kx),dir_qc(kx), &
              uu(kx),zuu_qc(kx), vv(kx),zvv_qc(kx), cld(kx),cld_qc(kx), &
              ciel(kx),ciel_qc(kx))

         p(:) = MISSING_R8
         z(:) = MISSING_R8
         t(:) = MISSING_R8
         td(:) = MISSING_R8
         spd(:) = MISSING_R8
         dir(:) = MISSING_R8
         uu(:) = MISSING_R8
         vv(:) = MISSING_R8
         cld(:) = MISSING_R8
         ciel(:) = MISSING_R8

         p_qc(:) = 0
         zpp_qc(:) = 0
         tt_qc(:) = 0
         td_qc(:) = 0
         spd_qc(:) = 0
         dir_qc(:) = 0
         zuu_qc(:) = 0
         zvv_qc(:) = 0
         cld_qc(:) = 0
         ciel_qc(:) = 0

         ! Write data at each (pressure) level.

         k = 1
         call get_obs_from_key(dart_seq, keys(is), obs)
         call get_obs_def(obs, obs_def)
         location = get_obs_def_location(obs_def)
         loc = get_location(location)
         pwhich_vert = nint(query_location(location,'which_vert'))
         pvloc = loc(3)
         do iobs = is, ie-1
            call get_obs_from_key(dart_seq, keys(iobs), obs)
            call get_obs_def(obs, obs_def)
            location = get_obs_def_location(obs_def)
            loc = get_location(location)
            which_vert = nint(query_location(location,'which_vert'))

            if (loc(3) /= pvloc .and. which_vert == pwhich_vert) k = k + 1
            pvloc = loc(3)
            pwhich_vert = which_vert

            if(which_vert == 2) then
               p(k) = loc(3)
!--------------------------------------------------------------------------
!  Flaging all data above 100 hPa.
!--------------------------------------------------------------------------
!!$               if(p(k) < 10000.0_r8) then
!!$                  p_qc(k) = MISSING_I
!!$                  zpp_qc(k) = MISSING_I
!!$                  tt_qc(k) = MISSING_I
!!$                  td_qc(k) = MISSING_I
!!$                  spd_qc(k) = MISSING_I
!!$                  dir_qc(k) = MISSING_I
!!$                  zuu_qc(k) = MISSING_I
!!$                  zvv_qc(k) = MISSING_I
!!$                  cld_qc(k) = MISSING_I
!!$                  ciel_qc(k) = MISSING_I
!!$               endif
            elseif(which_vert == 3) then
               z(k) = loc(3)
            endif

            obs_kind = map_type_of_obs_table(get_obs_def_type_of_obs(obs_def))

            call get_obs_values(obs, obs_value, 1)
            if(obs_kind == RADIOSONDE_U_WIND_COMPONENT) then
               num_obs = num_obs + 1
               uu(k) = obs_value(1)
            elseif(obs_kind == RADIOSONDE_V_WIND_COMPONENT) then
               num_obs = num_obs + 1
               vv(k) = obs_value(1)
               if(p(k) /= MISSING_R8) then
                  call StoreObsErr(get_obs_def_error_variance(obs_def), &
                       nint(0.01_r8 * p(k)), wind_pressures, n_wind_pres, wind_error)
               endif
            elseif(obs_kind == RADIOSONDE_TEMPERATURE) then
               num_obs = num_obs + 1
               t(k) = obs_value(1)
               if(p(k) /= MISSING_R8) then
                  call StoreObsErr(get_obs_def_error_variance(obs_def), &
                       nint(0.01_r8 * p(k)), pressure_levels, n_pres, temp_error)
               endif
!!$            elseif(obs_kind == RADIOSONDE_SPECIFIC_HUMIDITY) then
!!$               num_obs = num_obs + 1
!!$               t(k) = obs_value(1)
!!$               if(p(k) /= MISSING_R8) then
!!$                  call StoreObsErr(1.0e6_r8*get_obs_def_error_variance(obs_def), &
!!$                       nint(0.01_r8 * p(k)), pressure_levels, n_pres, specific_humidity_error)
!!$               endif
!!$            elseif(obs_kind == RADIOSONDE_SURFACE_PRESSURE) then
!!$               num_obs = num_obs + 1
            else
               if (n_no_support == 0) then
                  n_no_support = n_no_support + 1
                  obs_no_support(n_no_support) = get_name_for_type_of_obs(obs_kind)
               elseif (obs_no_support(n_no_support) /= get_name_for_type_of_obs(obs_kind)) then
                  n_no_support = n_no_support + 1
                  obs_no_support(n_no_support) = get_name_for_type_of_obs(obs_kind)
               endif
            endif
         enddo

         do k = 1 , kx

            ! Convert zonal and meridional wind components to
            ! wind speed and direction.

            if (uu(k) /= missing_r8 .and. vv(k) /= missing_r8) then

               spd(k) = sqrt(uu(k)*uu(k) + vv(k)*vv(k))

               if (spd(k) /= 0.0_r8) then

                  if (vv(k) == 0.0_r8) then

                     if (uu(k) > 0.0_r8) dir(k) = 270.0_r8
                     if (uu(k) < 0.0_r8) dir(k) = 90.0_r8

                  else

                     dir(k) = atan(uu(k)/vv(k))*RAD2DEG

                     if (vv(k) >= 0.0_r8) then
                        dir(k) = dir(k) + 180.0_r8
                     else
                        if (uu(k) >= 0.0_r8) dir(k) = dir(k) + 360.0_r8
                     endif

                  endif

               else

                  dir(k) = 0.0_r8

               endif

            endif

            WRITE ( UNIT = iunit , ERR = 20 , FMT = meas_format ) &
                 p(k),0, z(k),zpp_qc(k), t(k),tt_qc(k),td(k),td_qc(k), &
                 spd(k),zvv_qc(k), dir(k),zuu_qc(k), &
                 uu(k),zuu_qc(k), vv(k),zvv_qc(k), MISSING_R8,0, MISSING_R8,0

         enddo

         ! Write "footer" of the sounding.

         WRITE ( UNIT = iunit , ERR = 21 , FMT = meas_format ) &
              -777777.,0, -777777.,0,float(kx),0, &
              MISSING_R8,0, MISSING_R8,0, MISSING_R8,0, &
              MISSING_R8,0, MISSING_R8,0, MISSING_R8,0, MISSING_R8,0

         WRITE ( UNIT = iunit , ERR = 21 , FMT = end_format )  kx, 0, 0

         deallocate(p,p_qc, z,zpp_qc, t,tt_qc, td,td_qc, &
              spd,spd_qc, dir,dir_qc, &
              uu,zuu_qc, vv,zvv_qc, cld,cld_qc, &
              ciel,ciel_qc)

         is = ie

      endif

      ie = ie + 1

   enddo

   deallocate(keys)

   close(iunit)

   ! Write out the obserr.txt.

   open(unit=iunit,file='obserr.txt',status='new')

   write(iunit, FMT = ' ( 7f6.1 , "      BOGUS   WIND SENSOR ERRORS" ) ') wind_error(1:7)
   write(iunit, FMT = ' ( 7f6.1 , "        :" ) ') wind_error(8:14)
   write(iunit, FMT = ' ( 7f6.1 , "        :" ) ') wind_error(15:21)
   write(iunit, FMT = ' ( 7f6.1 , "        :" ) ') wind_error(22:26),wind_error(26),wind_error(26)
   write(iunit, FMT = ' ( 7f6.1 , "      RAOBS" ) ') wind_error(1:7)
   write(iunit, FMT = ' ( 7f6.1 , "        :" ) ') wind_error(8:14)
   write(iunit, FMT = ' ( 7f6.1 , "        :" ) ') wind_error(15:21)
   write(iunit, FMT = ' ( 7f6.1 , "        :" ) ') wind_error(22:26),wind_error(26),wind_error(26)
   write(iunit, FMT = ' ( 5f6.1 , "      BOGUS   TEMP SENSOR ERRORS" ) ') temp_error(1:5)
   write(iunit, FMT = ' ( 5f6.1 , "        :" ) ') temp_error(6:10)
   write(iunit, FMT = ' ( 5f6.1 , "        :" ) ') temp_error(11:15)
   write(iunit, FMT = ' ( 5f6.1 , "      RAOBS" ) ') temp_error(1:5)
   write(iunit, FMT = ' ( 5f6.1 , "        :" ) ') temp_error(6:10)
   write(iunit, FMT = ' ( 5f6.1 , "        :" ) ') temp_error(11:15)
   write(iunit, FMT = ' ( 5f6.2 , "      BOGUS   SPECIFIC HUMIDITY SENSOR ERRORS" ) ') specific_humidity_error(1:5)
   write(iunit, FMT = ' ( 5f6.2 , "        :" ) ') specific_humidity_error(6:10)
   write(iunit, FMT = ' ( 5f6.2 , "        :" ) ') specific_humidity_error(11:15)
   write(iunit, FMT = ' ( 5f6.2 , "      RAOBS" ) ') specific_humidity_error(1:5)
   write(iunit, FMT = ' ( 5f6.2 , "        :" ) ') specific_humidity_error(6:10)
   write(iunit, FMT = ' ( 5f6.2 , "        :" ) ') specific_humidity_error(11:15)

   close(iunit)

   if (n_no_support /= 0) then

      write(*,*) 'Observations not supported yet:'
      do k = 1 , n_no_support
         write(*,*) k,'. ',trim(adjustl(obs_no_support(k)))
      enddo

   endif

else

   ! Conversion from littler to DART format.

   ! Read in the observational errors from file 'obserr.txt'.

   filein = 'obserr.txt'
   platform = 'RAOBS'
   sensor_name = 'WIND SENSOR ERRORS'
   call READ_OBSERR (filein, platform, sensor_name, wind_error, n_wind_pres)
   sensor_name = 'TEMP SENSOR ERRORS'
   call READ_OBSERR (filein, platform, sensor_name, temp_error, n_pres)

   max_num_obs = 500000

   ! Initialize an obs_sequence structure
   call init_obs_sequence(dart_seq, num_copies, num_qc, max_num_obs)

   call set_copy_meta_data(dart_seq, 1, 'littler observations')
   call set_qc_meta_data(  dart_seq, 1, 'littler QC')

   iunit = open_file(littler_file_name, action = 'read')

   end_of_file=0

   ! Read littler until end of file.

   do while (end_of_file == 0)

      ! Read header of sounding.

      READ ( UNIT=iunit,ERR=19,FMT=rpt_format,iostat=end_of_file,end=25 ) &
           tst_xlat, tst_xlon, tst_id, tst_name, &
           tst_pltfrm, tst_src, tst_ter, kx, i1, i2, iseq_num, i3, &
           tst_sound, tst_bogus, tst_discard, tst_sut, tst_julian, date_char, &
           tst_slp, i6, f1,i7, f2,i8, f3,i9, tst_Psfc,i10, f5,i11, &
           f6,i12, f7,i13, f8,i14, f9,i15, f10,i16, f11,i17, f12,i18

      allocate(p(kx),p_qc(kx), z(kx),zpp_qc(kx), t(kx),tt_qc(kx), td(kx),td_qc(kx), &
           spd(kx),spd_qc(kx), dir(kx),dir_qc(kx), &
           uu(kx),zuu_qc(kx), vv(kx),zvv_qc(kx), cld(kx),cld_qc(kx), &
           ciel(kx),ciel_qc(kx))

      if (.not. tst_sound ) kx = 1

      ! Read each level of the sounding.

      do k = 1 , kx
         READ (UNIT=iunit,ERR=20,FMT=meas_format,iostat=end_of_file,end=25) &
              p(k),p_qc(k), z(k),zpp_qc(k), t(k),tt_qc(k), &
              td(k),td_qc(k), &
              spd(k),spd_qc(k), dir(k),dir_qc(k), &
              uu(k),zuu_qc(k), vv(k),zvv_qc(k), cld(k),cld_qc(k), &
              ciel(k),ciel_qc(k)

         if (p(k) /= missing_r8) then

            vloc = p(k)
            which_vert = 2

!!$         elseif (z(k)/= missing_r8) then
!!$
!!$            vloc = p(k)
!!$            which_vert = 3
!!$
         else

            call error_handler(E_ERR,'littler_tf_dart', &
                 'No vertical coordinate.', source, revision, revdate)

         endif

         ! In DART format, longitude is defined within [0, 360].

         if(tst_xlon < 0.0_r8) tst_xlon = tst_xlon + 360.0_r8
         lon = tst_xlon
         lat = tst_xlat

         location = set_location(lon, lat, vloc, which_vert)
         call set_obs_def_location(obs_def, location)

         call set_dart_time(date_char, time)
         call set_obs_def_time(obs_def, time)

         ! Insert obs into DART obs sequence.

         ! Observational errors are interpolated from data
         ! in obserr.txt

         if (t(k) /= missing_r8) then

            if (tst_pltfrm == 'FM-35 TEMP') then

               num_obs = num_obs + 1

               call set_obs_def_type_of_obs(obs_def, RADIOSONDE_TEMPERATURE)

               erms = intplog (p(k)*0.01_r8, pressure_levels, temp_error)

               call set_obs_def_error_variance(obs_def, erms*erms)

               call set_obs_def(obs, obs_def)

               obs_value(1) = t(k)
               call set_obs_values(obs, obs_value)
               qc(1) = tt_qc(k)   ! translate this from little-r to dart numbers?
               call set_qc(obs, qc)

               if(num_obs == 1) then
                  call insert_obs_in_seq(dart_seq, obs)
               else
                  call insert_obs_in_seq(dart_seq, obs, prev_obs)
               endif

               prev_obs = obs

            endif

         endif

         ! Littler may carry wind speed and direction:
         ! convert to zonal and meridional wind components.

         if (uu(k) == missing_r8 .and. vv(k) == missing_r8 .and. &
              spd(k) /= missing_r8 .and. dir(k) /= missing_r8) then

            dir(k) = dir(k)*DEG2RAD
            uu(k) = -spd(k)*SIN(dir(k))
            vv(k) = -spd(k)*COS(dir(k))

         endif

         if (uu(k) /= missing_r8) then

            if (tst_pltfrm == 'FM-35 TEMP') then

               num_obs = num_obs + 1

               call set_obs_def_type_of_obs(obs_def, RADIOSONDE_U_WIND_COMPONENT)

               erms = intplin (p(k)*0.01_r8, wind_pressures, wind_error)

               call set_obs_def_error_variance(obs_def, erms*erms)

               call set_obs_def(obs, obs_def)

               obs_value(1) = uu(k)
               call set_obs_values(obs, obs_value)
               qc(1) = zuu_qc(k)   ! translate this from little-r to dart numbers?
               call set_qc(obs, qc)

               if(num_obs == 1) then
                  call insert_obs_in_seq(dart_seq, obs)
               else
                  call insert_obs_in_seq(dart_seq, obs, prev_obs)
               endif

               prev_obs = obs

            endif

         endif

         if (vv(k) /= missing_r8) then

            if (tst_pltfrm == 'FM-35 TEMP') then

               num_obs = num_obs + 1

               call set_obs_def_type_of_obs(obs_def, RADIOSONDE_V_WIND_COMPONENT)

               erms = intplin (p(k)*0.01_r8, wind_pressures, wind_error)

               call set_obs_def_error_variance(obs_def, erms*erms)

               call set_obs_def(obs, obs_def)

               obs_value(1) = vv(k)
               call set_obs_values(obs, obs_value, 1)
               qc(1) = zvv_qc(k)   ! translate this from little-r to dart numbers?
               call set_qc(obs, qc)

               if(num_obs == 1) then
                  call insert_obs_in_seq(dart_seq, obs)
               else
                  call insert_obs_in_seq(dart_seq, obs, prev_obs)
               endif

               prev_obs = obs

            endif

         endif

      enddo

      deallocate(p,p_qc, z,zpp_qc, t,tt_qc, td,td_qc, &
           spd,spd_qc, dir,dir_qc, &
           uu,zuu_qc, vv,zvv_qc, cld,cld_qc, &
           ciel,ciel_qc)

      READ ( UNIT = iunit , ERR = 21 , FMT = meas_format ) &
           f1,i1, f2,i2, f3,i3, f4,i4, f5,i5, f6,i6, &
           f7,i7, f8,i8, f9,i9, f10,i10

      READ ( UNIT = iunit , ERR = 21 , FMT = end_format )  kx, i1, i2

   enddo

25 print*,"FOUND END"

   close(iunit)
 
! Write out the sequence
   call write_obs_seq(dart_seq, dart_file_name)

endif

write(unit=*, fmt='(5x,a,i7)') &
     'Total number of observations:  ', num_obs

goto 26

! These are for error messages when regarding littler files.

19 continue
call error_handler(E_ERR,'littler_tf_dart', &
     'Error when reading or writing header of sounding.', source, revision, revdate)
20 continue
call error_handler(E_ERR,'littler_tf_dart', &
     'Error when reading or writing sounding.', source, revision, revdate)
21 continue
call error_handler(E_ERR,'littler_tf_dart', &
     'Error when reading or writing footer of the sounding.', source, revision, revdate)

26 continue

call error_handler(E_MSG, 'littler_tf_dart', 'FINISHED littler_tf_dart.')
call error_handler(E_MSG, 'littler_tf_dart', 'Finished successfully.',&
                   source,revision,revdate)
call finalize_utilities()


contains

!#######################################################

subroutine set_str_date(timestring, dart_time)

implicit none

type(time_type),   intent(in)  :: dart_time
character(len=20), intent(out) :: timestring

integer           :: year, month, day, hour, minute, second

character(len=4)  :: ch_year
character(len=2)  :: ch_month, ch_day, ch_hour, ch_minute, ch_second

call get_date(dart_time, year, month, day, hour, minute, second)

write(ch_year,'(i4)') year
write(ch_month,'(i2.2)') month
write(ch_day,'(i2.2)') day
write(ch_hour,'(i2.2)') hour
write(ch_minute,'(i2.2)') minute
write(ch_second,'(i2.2)') second

timestring(1:6)   = "      "
timestring(7:10)  = ch_year
timestring(11:12) = ch_month
timestring(13:14) = ch_day
timestring(15:16) = ch_hour
timestring(17:18) = ch_minute
timestring(19:20) = ch_second

end subroutine set_str_date
 
!#######################################################

subroutine set_dart_time(tstring, dart_time)

implicit none

type(time_type),   intent(out) :: dart_time
character(len=20), intent(in)  :: tstring

integer           :: year, month, day, hour, minute, second

read(tstring(7:10),'(i4)') year
read(tstring(11:12),'(i2)') month
read(tstring(13:14),'(i2)') day
read(tstring(15:16),'(i2)') hour
read(tstring(17:18),'(i2)') minute
read(tstring(19:20),'(i2)') second

dart_time = set_date(year, month, day, hour, minute, second)

end subroutine set_dart_time

!#######################################################

Subroutine StoreObsErr(obs_err_var, pres, plevel, nlev, obs_err_std)

implicit none

integer,  intent(in)    :: nlev, pres
real(r8), intent(in)    :: obs_err_var
integer,  intent(in)    :: plevel(nlev)
real(r8), intent(inout) :: obs_err_std(nlev)

integer :: level_index


level_index = GetClosestLevel(pres, plevel, nlev)

if ( plevel(level_index) == pres ) then
   if ( obs_err_std(level_index) == -1.0_R8 ) then
      obs_err_std(level_index) = sqrt(obs_err_var)
   elseif (obs_err_std(level_index) /= sqrt(obs_err_var)) then
      print*,'Various observation errors at same level: ', &
           obs_err_std(level_index), sqrt(obs_err_var), ' level = ',pres
   endif
endif

end Subroutine StoreObsErr

!#######################################################

Function GetClosestLevel(ilev, vlev, nlev) result (level_index)

implicit none

integer,  intent(in) :: nlev, ilev
integer,  intent(in) :: vlev(nlev)

integer                  :: level_index, a(1)
integer, dimension(nlev) :: dx

dx = abs(ilev - vlev)
a = minloc(dx)
level_index = a(1)

end Function GetClosestLevel

!#######################################################

SUBROUTINE READ_OBSERR (filein, platform, sensor_name, err, nlevels)
!------------------------------------------------------------------------------!
!
! Read observational error on pressure levels (in hPa)
!
!------------------------------------------------------------------------------!
  IMPLICIT NONE
!------------------------------------------------------------------------------!

  CHARACTER (LEN=80), intent(in) :: filein      != 'obserr.txt'
  CHARACTER (LEN=80), intent(in) :: platform    != 'RAOBS'
  CHARACTER (LEN=80), intent(in) :: sensor_name != 'WIND SENSOR ERRORS'
  INTEGER,            intent(in) :: nlevels

  REAL(r8), intent(out) :: err(nlevels)

  character(len=129)  :: msgstring
  INTEGER             :: io_error, i, iunit, ncol
  CHARACTER (LEN=80)  :: fmt_err                    != '(5(1X,F5.1))'
  CHARACTER (LEN=80)  :: line1, line2, line3, line4
  LOGICAL             :: sensor_found, found

!------------------------------------------------------------------------------!

! 1.  OPEN INPUT FILE
! ===================

  iunit = get_unit()
  OPEN (UNIT = iunit , FILE = filein , FORM = 'FORMATTED'  , &
       ACTION = 'READ' , STATUS = 'OLD', IOSTAT = io_error)

  IF (io_error /= 0) THEN
     write(msgstring, *) 'Unable to open input observational error file ',TRIM (filein)
     call error_handler(E_ERR,'READ_OBSERR', &
          msgstring, source, revision, revdate)
  ENDIF


! 2.  READ DATA
! =============

! Read file until sensor_name is found
  found          = .FALSE.
  sensor_found   = .FALSE.
  io_error= 0

  DO WHILE (io_error == 0.)

     READ (UNIT = iunit, IOSTAT = io_error, FMT = '(A)') line1
     !  Exit when error or at end of file

     IF (io_error /= 0 .OR. (line1 (1:2) == '*.')) EXIT

     IF (TRIM (sensor (line1)) == 'BOGUS') THEN

        IF (TRIM (obstype (line1)) == TRIM (sensor_name)) sensor_found = .TRUE.

     ENDIF

     IF (sensor_found .and. (TRIM (sensor (line1)) == TRIM (platform))) THEN

        found = .TRUE.
        READ (UNIT = iunit, IOSTAT = io_error, FMT = '(A)') line2
        READ (UNIT = iunit, IOSTAT = io_error, FMT = '(A)') line3
        IF (TRIM (sensor_name) == 'WIND SENSOR ERRORS') &
          READ (UNIT = iunit, IOSTAT = io_error, FMT = '(A)') line4
        EXIT

     ELSE

        !  If obstype is not found, keep on reading

        CYCLE

     ENDIF
        
  ENDDO

  IF (found) THEN

     !  Sensor_name has been found, Error at mandatory pressure levels follow
     !  Break down data upon obs type

     SELECT CASE (TRIM (sensor_name))
     CASE ('WIND SENSOR ERRORS')
        ncol = 7
        fmt_err = '(7(1X,F5.1))'
     CASE ('RH SENSOR ERRORS')
        ncol = 5
        fmt_err = '(5(1X,F5.2))'
     CASE DEFAULT
        ncol = 5
        fmt_err = '(5(1X,F5.1))'
     END SELECT

     READ (line1, fmt_err) (err (i), i =        1,   ncol)
     READ (line2, fmt_err) (err (i), i =   ncol+1, 2*ncol)
     READ (line3, fmt_err) (err (i), i = 2*ncol+1, 3*ncol)
     !  Winds are given over 4 lines, any other data are given over 3 lines
     IF (TRIM (sensor_name) == 'WIND SENSOR ERRORS') &
     READ (line4, fmt_err) (err (i), i = 3*ncol+1, min(4*ncol,nlevels))

  ELSE

     WRITE (*, FMT = '(/,A,A,A,/)') TRIM (sensor_name), &
          ' were not found in file ', TRIM (filein)

  ENDIF


! 3.  CLOSE INPUT FILE
! ====================

  CLOSE (UNIT = iunit)

END SUBROUTINE READ_OBSERR

!#######################################################

FUNCTION obstype (line) RESULT (f_obstype)
!------------------------------------------------------------------------------!

  ! Read in a line the string present after keyword 'BOGUS'

!------------------------------------------------------------------------------!
  IMPLICIT NONE
  CHARACTER (LEN= 80), intent(in) :: line

  CHARACTER (LEN= 80) :: f_obstype
  INTEGER :: b,c
!------------------------------------------------------------------------------!

!  Find keyword bogus

  DO c = 1, LEN_TRIM (line) - 5
     IF (line (c:c+4) == 'BOGUS') EXIT
  ENDDO

  c = c + 5

!  Skip blank until next word

  DO b = c, LEN_TRIM (line)
     IF (line (b:b) /= ' ') EXIT
  ENDDO

!  String follows

  f_obstype = TRIM (line (b:LEN_TRIM (line)))

END FUNCTION obstype

!#######################################################

FUNCTION sensor (line) RESULT (f_sensor)

  ! Read first in a string after numbers

!------------------------------------------------------------------------------!
  IMPLICIT NONE
  CHARACTER (LEN= 80), intent(in) :: line

  CHARACTER (LEN= 80) :: f_sensor
  INTEGER :: b,c
!------------------------------------------------------------------------------!

!  Find the first non-blank, non point and non-number character

  DO c = 1, LEN_TRIM (line)
     IF (((iachar (line(c:c)) /= 32)  .AND. &
          (iachar (line(c:c)) /= 46)) .AND. &
          ((iachar (line(c:c)) .LT. 48)  .OR.  &   
          (iachar (line(c:c)) .GT. 57)))      &
          EXIT
  ENDDO

  DO b = c, LEN_TRIM (line)
     IF (line (b:b) == ' ') EXIT
  ENDDO

  f_sensor = line (c:b-1)

END FUNCTION sensor

!#######################################################

FUNCTION intplin (x,xx,yy) RESULT (val)
!------------------------------------------------------------------------------!
  IMPLICIT NONE

  INTEGER,  DIMENSION (:), intent(in) :: xx
  REAL(r8), DIMENSION (:), intent(in) :: yy
  REAL(r8),                intent(in) :: x

  REAL(r8) :: val

  INTEGER :: n,m,jl
!------------------------------------------------------------------------------!

  n = size (xx)
  m = size (yy)

  IF (n /= m) THEN
     call error_handler(E_ERR,'intplin', &
          'arrays xx and yy must have same size', source,revision,revdate)
  ENDIF

  jl = locate (x,xx)

  IF (jl .LE. 0) THEN
     if (yy (1) >= 0.0_r8) then
        val = yy (1)
     else
        call error_handler(E_ERR,'intplin', &
             'bad value in yy(1)', source,revision,revdate)
     endif
  ELSE IF (jl .GE. n) THEN    
     if (yy (n) >= 0.0_r8) then
        val = yy (n)
     else
        call error_handler(E_ERR,'intplin', &
             'bad value in yy(n)', source,revision,revdate)
     endif
  ELSE
     if (yy (jl) >= 0.0_r8 .and. yy (jl+1) >= 0.0_r8) then
        val = (xx (jl+1) - x) * yy (jl) + (x - xx (jl)) * yy (jl+1)
        val = val / (xx (jl+1) - xx (jl))
     else
        call error_handler(E_ERR,'intplin', &
             'bad value in yy(jl) or yy(jl+1)', source,revision,revdate)
     endif
  ENDIF

END FUNCTION intplin

!#######################################################

FUNCTION intplog (x,xx,yy) RESULT (val)
!------------------------------------------------------------------------------!
  IMPLICIT NONE

  INTEGER,  DIMENSION (:), intent(in) :: xx
  REAL(r8), DIMENSION (:), intent(in) :: yy
  REAL(r8),                intent(in) :: x

  REAL(r8) :: val

  INTEGER :: n,m,jl
!------------------------------------------------------------------------------!

  n = size (xx)
  m = size (yy)

  IF (n /= m) THEN
     call error_handler(E_ERR,'intplog', &
          'arrays xx and yy must have same size', source,revision,revdate)
  ENDIF

  jl = locate (x,xx)

  IF (jl .LE. 0) THEN    
     if (yy (1) >= 0.0_r8) then
        val = yy (1)
     else
        call error_handler(E_ERR,'intplog', &
             'bad value in yy(1)', source,revision,revdate)
     endif
  ELSE IF (jl .GE. n) THEN    
     if (yy (n) >= 0.0_r8) then
        val = yy (n)
     else
        call error_handler(E_ERR,'intplog', &
             'bad value in yy(n)', source,revision,revdate)
     endif
  ELSE
     if (yy (jl) >= 0.0_r8 .and. yy (jl+1) >= 0.0_r8) then
        val = log (real(xx(jl+1)) / x) * yy (jl) + log (x / real(xx(jl))) * yy (jl+1)
        val = val / log (real(xx(jl+1)) / real(xx(jl)))
     else
        call error_handler(E_ERR,'intplog', &
             'bad value in yy(jl) or yy(jl+1)', source,revision,revdate)
     endif
  ENDIF

END FUNCTION intplog

FUNCTION locate (x,xx) RESULT (level_index)
!------------------------------------------------------------------------------!
  IMPLICIT NONE

  INTEGER, DIMENSION (:), intent(in) :: xx
  REAL(r8),               intent(in) :: x

  INTEGER :: level_index

  INTEGER :: n,jl,jm,ju
  LOGICAL :: ascnd
!------------------------------------------------------------------------------!

  n = size (xx)
  ascnd = (xx (n) >= xx (1))   ! True if ascending order, false otherwise
  jl = 0                       ! Initialize lower limit
  ju = n+1                     ! Initialize upper limit

  DO

     IF (ju-jl <= 1) EXIT      ! Repeat until this condition is satisfied

     jm = (ju+jl) / 2.         ! Compute a mid point

     IF (ascnd .EQV. (x >= xx (jm))) THEN
        jl = jm               ! Replace mid point by lower limit
     ELSE
        ju = jm               ! Replace mid point by upper limit
     ENDIF

  ENDDO

  IF (x == xx (1)) THEN      ! Set the output, being carefull with endpoints
     level_index = 1
  ELSE IF (x == xx (n)) THEN
     level_index = n-1 
  ELSE
     level_index = jl
  ENDIF

END FUNCTION LOCATE

END PROGRAM littler_tf_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
