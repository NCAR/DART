! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

PROGRAM littler_tf_dart

! <next three lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$

use        types_mod, only : r8, DEG2RAD, RAD2DEG, MISSING_I, MISSING_R8
use    utilities_mod, only : open_file, check_nml_error, close_file, file_exist, &
                             get_unit, initialize_utilities, &
                             finalize_utilities, register_module, error_handler, E_ERR, &
                             logfileunit
use obs_sequence_mod, only : obs_type, obs_sequence_type, init_obs_sequence, &
                             insert_obs_in_seq, write_obs_seq, read_obs_seq, &
                             set_qc, set_obs_values, set_copy_meta_data, &
                             assignment(=), &
                             init_obs, static_init_obs_sequence, get_num_qc, set_obs_def, &
                             get_num_obs, get_max_num_obs, get_obs_values, &
                             get_obs_from_key, get_obs_def, get_qc
use      obs_def_mod, only : copy_obs_def, obs_def_type, &
                             get_obs_def_time, get_obs_def_location, &
                             get_obs_def_error_variance, get_obs_def_kind, &
                             set_obs_def_kind, set_obs_def_location, set_obs_def_time, &
                             set_obs_def_error_variance
use     obs_kind_mod, only : KIND_U, KIND_V, KIND_T, &
                             get_obs_kind, set_obs_kind
use     location_mod, only : location_type, get_location, query_location, set_location
use time_manager_mod, only : time_type, get_time, set_calendar_type, GREGORIAN, get_date, &
                             set_date, operator(/=)

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(obs_sequence_type) :: dart_seq
type(obs_type)          :: obs, prev_obs
type(obs_def_type)      :: obs_def
type(location_type)     :: location
type(time_type)         :: time, stime, ftime

INTEGER                 :: iunit, endfile

integer                 :: kind, which_vert
real(r8)                :: obs_value(1), qc(1)
real(r8), dimension(3)  :: loc, sloc, floc
real(r8)                :: lon, lat, vloc, pvloc

integer           :: is, ie, iobs, k
integer           :: num_obs, num_copies, num_qc, max_num_obs, dart_seq_num_obs

character(len = 129) :: dart_file_name     = 'obs_seq.out', &
                        littler_file_name = 'little-r.dat', &
                        errstring
integer              :: calendar_type      = GREGORIAN

character *20  :: date_char
CHARACTER *120 :: rpt_format
CHARACTER *120 :: meas_format
CHARACTER *120 :: end_format
character *40  :: tst_id, tst_name, tst_pltfrm, tst_src
real(r8)       :: tst_ter, tst_slp, tst_xlat, tst_xlon
integer        :: kx, iseq_num, tst_sut, tst_julian
logical        :: tst_sound, tst_bogus, tst_discard
integer        :: i1,i2,i3,i4,i5,i6,i7,i8,i9,i10
integer        :: i11,i12,i13,i14,i15,i16,i17,i18
real(r8)       :: f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12

real(r8), allocatable :: z(:),t(:),td(:),spd(:),dir(:), &
                         uu(:), vv(:), p(:), cld(:), ciel(:)
integer, allocatable  :: zpp_qc(:), tt_qc(:), td_qc(:), &
                         zuu_qc(:), zvv_qc(:), p_qc(:), &
                         spd_qc(:), dir_qc(:), cld_qc(:), ciel_qc(:)

logical :: littler_to_dart

!------------------------------------------------------------------------------

rpt_format =  ' ( 2f20.5 , 4a40 , f20.5 , 5i10 , 3L10 , ' &
                  // ' 2i10 , a20 ,  13( f13.5 , i7 ) ) '
meas_format =  ' ( 10( f13.5 , i7 ) ) '
end_format = ' ( 3 ( i7 ) ) '

call initialize_utilities
call register_module(source, revision, revdate)
write(logfileunit,*)'STARTING littler_tf_dart ...'

call set_calendar_type(calendar_type)

read(5,*) littler_to_dart

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

call static_init_obs_sequence()

num_copies = 1
num_qc = 0
call init_obs(obs, num_copies, num_qc)

! -------------------------------------------------------------------
! Initialize the counter:

num_obs = 0

if(.not. littler_to_dart) then

   call read_obs_seq(dart_file_name, 0, 0, 0, dart_seq)

   dart_seq_num_obs = get_num_obs(dart_seq)

   iunit = get_unit()

   open(unit=iunit,file=littler_file_name,status='new')

   is = 1
   ie = 1
   do while (ie <= dart_seq_num_obs+1)

      call get_obs_from_key(dart_seq, is, obs)
      call get_obs_def(obs, obs_def)

      sloc = get_location(get_obs_def_location(obs_def))
      stime = get_obs_def_time(obs_def)

      if(ie <= dart_seq_num_obs) then
         call get_obs_from_key(dart_seq, ie, obs)
         call get_obs_def(obs, obs_def)
      endif

      floc = get_location(get_obs_def_location(obs_def))
      ftime = get_obs_def_time(obs_def)

      if ( ftime /= stime .or. sloc(1) /= floc(1) .or. &
           sloc(2) /= floc(2) .or. ie == dart_seq_num_obs+1) then

         kx = 1
         call get_obs_from_key(dart_seq, is, obs)
         call get_obs_def(obs, obs_def)
         loc = get_location(get_obs_def_location(obs_def))
         pvloc = loc(3)
         do iobs = is+1, ie-1
            call get_obs_from_key(dart_seq, iobs, obs)
            call get_obs_def(obs, obs_def)
            loc = get_location(get_obs_def_location(obs_def))
            vloc = loc(3)
            if (vloc /= pvloc) kx = kx + 1
            pvloc = vloc
         enddo

!!$         if (kx > 1) then
            tst_sound = .true.
            tst_pltfrm = 'FM-35 TEMP'
!!$         else
!!$            tst_sound = .false.
!!$         endif

         call set_str_date(date_char, stime)

         tst_xlat = loc(2)
         tst_xlon = loc(1)
         if (tst_xlon > 180.0_r8) tst_xlon = tst_xlon - 360.0_r8

         WRITE ( UNIT = iunit , ERR = 19 , FMT = rpt_format ) &
              tst_xlat, tst_xlon, tst_id , tst_name, &
              tst_pltfrm, tst_src, tst_ter, kx, i1, i2, iseq_num, i3, &
              tst_sound, tst_bogus, tst_discard, tst_sut, tst_julian, date_char, &
              tst_slp, i6, f1, i7, f2, i8, f3, i9, f4, i10, &
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

         k = 1
         call get_obs_from_key(dart_seq, is, obs)
         call get_obs_def(obs, obs_def)
         loc = get_location(get_obs_def_location(obs_def))
         pvloc = loc(3)
         do iobs = is, ie-1
            call get_obs_from_key(dart_seq, iobs, obs)
            call get_obs_def(obs, obs_def)
            location = get_obs_def_location(obs_def)
            loc = get_location(location)
            if (loc(3) /= pvloc) k = k + 1
            pvloc = loc(3)
            which_vert = nint(query_location(location,'which_vert'))
            
            if(which_vert == 2) then
               p(k) = loc(3)
            elseif(which_vert == 3) then
               z(k) = loc(3)
            endif

            kind = get_obs_kind(get_obs_def_kind(obs_def))

            call get_obs_values(obs, obs_value, 1)
            if(kind == KIND_U) then
               num_obs = num_obs + 1
               uu(k) = obs_value(1)
            elseif(kind == KIND_V) then
               num_obs = num_obs + 1
               vv(k) = obs_value(1)
            elseif(kind == KIND_T) then
               num_obs = num_obs + 1
               t(k) = obs_value(1)
            else
               write(errstring,*)'Observation type ', &
                    kind,' not implemented yet.'
               call error_handler(E_ERR,'littler_tf_dart',errstring, &
                    source, revision, revdate)
            endif
         enddo

         do k = 1 , kx

            if (uu(k) /= missing_r8 .and. vv(k) /= missing_r8) then

               spd(k) = sqrt(uu(k)*uu(k) + vv(k)*vv(k))

               if (spd(k) /= 0.0_r8) then

                  if (vv(k) == 0.0_r8) then

                     if (uu(k) > 0.0_r8) dir(k) = 270.0_r8
                     if (uu(k) < 0.0_r8) dir(k) = 90.0_r8

                  else

                     dir(k) = atan(uu(k)/vv(k))*RAD2DEG

                     if (uu(k) <= 0.0_r8 .and. vv(k) >= 0.0_r8) dir(k) = dir(k) + 180.0_r8
                     if (uu(k) >= 0.0_r8 .and. vv(k) >= 0.0_r8) dir(k) = dir(k) + 180.0_r8
                     if (uu(k) >= 0.0_r8 .and. vv(k) <= 0.0_r8) dir(k) = dir(k) + 360.0_r8

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

         WRITE ( UNIT = iunit , ERR = 21 , FMT = meas_format ) &
              -777777.,0, -777777.,0,float(kx),0, &
              MISSING_R8,0, MISSING_R8,0, MISSING_R8,0, &
              MISSING_R8,0, MISSING_R8,0, MISSING_R8,0, MISSING_R8,0

         WRITE ( UNIT = iunit , ERR = 19 , FMT = end_format )  kx, 0, 0

         deallocate(p,p_qc, z,zpp_qc, t,tt_qc, td,td_qc, &
              spd,spd_qc, dir,dir_qc, &
              uu,zuu_qc, vv,zvv_qc, cld,cld_qc, &
              ciel,ciel_qc)

         is = ie

      endif

      ie = ie + 1

   enddo

   close(iunit)

else

   max_num_obs = 500000

! Initialize an obs_sequence structure
   call init_obs_sequence(dart_seq, num_copies, num_qc, max_num_obs)

   call set_copy_meta_data(dart_seq, 1, 'observations')

   iunit = open_file(littler_file_name, action = 'read')

   endfile=0
   do while (endfile .eq. 0)
      READ ( UNIT=iunit,ERR=19,FMT=rpt_format,iostat=endfile,end=25 ) &
           tst_xlat, tst_xlon, tst_id, tst_name, &
           tst_pltfrm, tst_src, tst_ter, kx, i1, i2, iseq_num, i3, &
           tst_sound, tst_bogus, tst_discard, tst_sut, tst_julian, date_char, &
           tst_slp, i6, f1,i7, f2,i8, f3,i9, f4,i10, f5,i11, &
           f6,i12, f7,i13, f8,i14, f9,i15, f10,i16, f11,i17, f12,i18

      print*,'kx = ',kx

      allocate(p(kx),p_qc(kx), z(kx),zpp_qc(kx), t(kx),tt_qc(kx), td(kx),td_qc(kx), &
           spd(kx),spd_qc(kx), dir(kx),dir_qc(kx), &
           uu(kx),zuu_qc(kx), vv(kx),zvv_qc(kx), cld(kx),cld_qc(kx), &
           ciel(kx),ciel_qc(kx))

      if (.not. tst_sound ) kx = 1

      do k = 1 , kx
         READ (UNIT=iunit,ERR=20,FMT=meas_format,iostat=endfile,end=25) &
              p(k),p_qc(k), z(k),zpp_qc(k), t(k),tt_qc(k), &
              td(k),td_qc(k), &
              spd(k),spd_qc(k), dir(k),dir_qc(k), &
              uu(k),zuu_qc(k), vv(k),zvv_qc(k), cld(k),cld_qc(k), &
              ciel(k),ciel_qc(k)

         if (p(k) /= missing_r8) then

            vloc = p(k)
            which_vert = 2

         elseif (z(k)/= missing_r8) then

            vloc = p(k)
            which_vert = 3

         else

            call error_handler(E_ERR,'littler_tf_dart', &
                 'No vertical coordinate.', source, revision, revdate)

         endif

         if(tst_xlon < 0.0_r8) tst_xlon = tst_xlon + 360.0_r8
         lon = tst_xlon
         lat = tst_xlat

         location = set_location(lon, lat, vloc, which_vert)
         call set_obs_def_location(obs_def, location)

         call set_dart_time(date_char, time)
         call set_obs_def_time(obs_def, time)

         if (t(k) /= missing_r8) then

            num_obs = num_obs + 1

            call set_obs_def_kind(obs_def, set_obs_kind(KIND_T))

            call set_obs_def_error_variance(obs_def, 1.0_r8)

            call set_obs_def(obs, obs_def)

            obs_value(1) = t(k)
            call set_obs_values(obs, obs_value)
            qc(1) = tt_qc(k)
!!$            call set_qc(obs, qc)

            if(num_obs == 1) then
               call insert_obs_in_seq(dart_seq, obs)
            else
               call insert_obs_in_seq(dart_seq, obs, prev_obs)
            endif

            prev_obs = obs

         endif

         if (uu(k) == missing_r8 .and. vv(k) == missing_r8 .and. &
              spd(k) /= missing_r8 .and. dir(k) /= missing_r8) then

            dir(k) = dir(k)*DEG2RAD
            uu(k) = -spd(k)*SIN(dir(k))
            vv(k) = -spd(k)*COS(dir(k))

         endif

         if (uu(k) /= missing_r8) then

            num_obs = num_obs + 1

            call set_obs_def_kind(obs_def, set_obs_kind(KIND_U))

            call set_obs_def_error_variance(obs_def, 1.0_r8)

            call set_obs_def(obs, obs_def)

            obs_value(1) = uu(k)
            call set_obs_values(obs, obs_value)
            qc(1) = zuu_qc(k)
!!$            call set_qc(obs, qc)

            if(num_obs == 1) then
               call insert_obs_in_seq(dart_seq, obs)
            else
               call insert_obs_in_seq(dart_seq, obs, prev_obs)
            endif

            prev_obs = obs

         endif

         if (vv(k) /= missing_r8) then

            num_obs = num_obs + 1

            call set_obs_def_kind(obs_def, set_obs_kind(KIND_V))

            call set_obs_def_error_variance(obs_def, 1.0_r8)

            call set_obs_def(obs, obs_def)

            obs_value(1) = vv(k)
            call set_obs_values(obs, obs_value, 1)
            qc(1) = zvv_qc(k)
!!$            call set_qc(obs, qc)

            if(num_obs == 1) then
               call insert_obs_in_seq(dart_seq, obs)
            else
               call insert_obs_in_seq(dart_seq, obs, prev_obs)
            endif

            prev_obs = obs

         endif

      enddo

      deallocate(p,p_qc, z,zpp_qc, t,tt_qc, td,td_qc, &
           spd,spd_qc, dir,dir_qc, &
           uu,zuu_qc, vv,zvv_qc, cld,cld_qc, &
           ciel,ciel_qc)

      READ ( UNIT = iunit , ERR = 21 , FMT = meas_format ) &
           f1,i1, f2,i2, f3,i3, f4,i4, f5,i5, f6,i6, &
           f7,i7, f8,i8, f9,i9, f10,i10

      READ ( UNIT = iunit , ERR = 19 , FMT = end_format )  kx, i1, i2

   enddo

25 print*,"FOUND END"

   close(iunit)
 
! Write out the sequence
   call write_obs_seq(dart_seq, dart_file_name)

endif

write(unit=*, fmt='(5x,a,i7)') &
     'Total number of observations:  ', num_obs

goto 26

19 continue
print *,' err=19 '
stop 19
20 continue
print *,' err=20 '
stop 20
21 continue
print *,' err=21 '
stop 21

26 continue

write(logfileunit,*)'FINISHED littler_tf_dart.'
write(logfileunit,*)

call finalize_utilities ! closes the log file.

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
 
END PROGRAM littler_tf_dart
