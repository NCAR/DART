! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module quikscat_JPL_mod

use types_mod,        only : r4, r8, digits12, deg2rad, rad2deg

use obs_def_mod,      only : obs_def_type, get_obs_def_time, read_obs_def, &
                             write_obs_def, destroy_obs_def, interactive_obs_def, &
                             copy_obs_def, set_obs_def_time, set_obs_def_type_of_obs, &
                             set_obs_def_error_variance, set_obs_def_location

use time_manager_mod, only : time_type, get_date, set_date, get_time, set_time, &
                             set_calendar_type, GREGORIAN, print_date, print_time, &
                             operator(==), operator(>), operator(<), operator(>=), &
                             operator(/=), operator(+), operator(-)

use    utilities_mod, only : get_unit, open_file, close_file, file_exist, &
                             register_module, error_handler, &
                             E_ERR, E_MSG, is_longitude_between

use     location_mod, only : location_type, set_location, VERTISHEIGHT

use     obs_kind_mod, only : get_index_for_type_of_obs, &
                             QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT

use obs_kind_mod,     only : QKSWND_U_WIND_COMPONENT, QKSWND_V_WIND_COMPONENT

use obs_sequence_mod, only : init_obs_sequence, init_obs, insert_obs_in_seq, &
                             set_obs_values, set_qc, obs_sequence_type, obs_type, &
                             copy_obs, set_copy_meta_data, set_qc_meta_data, set_obs_def, &
                             get_first_obs, get_last_obs, get_obs_def

implicit none
private

public :: real_obs_sequence, read_qscat2b, orbit_type, create_output_filename

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.
character(len=128) :: msgstring

logical :: DEBUG = .false.

! For the 25.0km resolution WVC ... MAX_ROWS = 1624 MAX_WVC =  76
! For the 12.5km resolution WVC ... MAX_ROWS = 3248 MAX_WVC = 152

integer :: MAX_ROWS, MAX_WVC, MAX_AMBIG
parameter (MAX_ROWS  = 1624)
parameter (MAX_WVC   = 76)
parameter (MAX_AMBIG = 4)

type orbit_type
   private
   character(len=32)                   :: solution

   type(time_type),dimension(MAX_ROWS) :: row_time
   real(r4),       dimension(MAX_ROWS) :: wvc_row

 ! real(r4), dimension(MAX_WVC,MAX_ROWS) :: wvc_index
 ! real(r4), dimension(MAX_WVC,MAX_ROWS) :: num_in_fore
 ! real(r4), dimension(MAX_WVC,MAX_ROWS) :: num_in_aft
 ! real(r4), dimension(MAX_WVC,MAX_ROWS) :: num_out_fore
 ! real(r4), dimension(MAX_WVC,MAX_ROWS) :: num_out_aft
 ! real(r4), dimension(MAX_WVC,MAX_ROWS) :: num_ambigs
 ! real(r4), dimension(MAX_WVC,MAX_ROWS) :: wvc_selection
 ! real(r4), dimension(MAX_WVC,MAX_ROWS) :: nof_rain_index
   integer,  dimension(MAX_WVC,MAX_ROWS) :: wvc_index
   integer,  dimension(MAX_WVC,MAX_ROWS) :: num_in_fore
   integer,  dimension(MAX_WVC,MAX_ROWS) :: num_in_aft
   integer,  dimension(MAX_WVC,MAX_ROWS) :: num_out_fore
   integer,  dimension(MAX_WVC,MAX_ROWS) :: num_out_aft
   integer,  dimension(MAX_WVC,MAX_ROWS) :: num_ambigs
   integer,  dimension(MAX_WVC,MAX_ROWS) :: wvc_selection
   integer,  dimension(MAX_WVC,MAX_ROWS) :: nof_rain_index

   real(r4), dimension(MAX_WVC,MAX_ROWS) :: wvc_lat, wvc_lon
   real(r4), dimension(MAX_WVC,MAX_ROWS) :: wvc_quality_flag
   real(r4), dimension(MAX_WVC,MAX_ROWS) :: atten_corr
   real(r4), dimension(MAX_WVC,MAX_ROWS) :: mp_rain_probability
   real(r4), dimension(MAX_WVC,MAX_ROWS) :: srad_rain_rate
   
   real(r4), dimension(MAX_AMBIG,MAX_WVC,MAX_ROWS) :: wind_speed
   real(r4), dimension(MAX_AMBIG,MAX_WVC,MAX_ROWS) :: wind_dir
   real(r4), dimension(MAX_AMBIG,MAX_WVC,MAX_ROWS) :: wind_speed_err
   real(r4), dimension(MAX_AMBIG,MAX_WVC,MAX_ROWS) :: wind_dir_err
   real(r4), dimension(MAX_AMBIG,MAX_WVC,MAX_ROWS) :: max_likelihood_est

   real(r4), dimension(MAX_WVC,MAX_ROWS) :: model_speed          ! NWP Wind Vector
   real(r4), dimension(MAX_WVC,MAX_ROWS) :: model_dir            ! NWP Wind Vector
   real(r4), dimension(MAX_WVC,MAX_ROWS) :: wind_speed_selection ! DRE Wind Vector
   real(r4), dimension(MAX_WVC,MAX_ROWS) :: wind_dir_selection   ! DRE Wind Vector
   
end type orbit_type

contains


subroutine initialize_module
!-------------------------------------------------
call register_module(source, revision, revdate)

call set_calendar_type(GREGORIAN)

module_initialized = .true.

end subroutine initialize_module


subroutine create_output_filename(l2bname, ofname)
!-------------------------------------------------
! The L2b filenames have a very long extension that
! records when the data was published - not very interesting
! for our purposes. replace with something DART-y.
character(len=*), intent(IN)  :: l2bname
character(len=*), intent(OUT) :: ofname

integer :: i, strlen, extstart

strlen = len_trim(l2bname)

dotloop : do i = strlen,1,-1
   if (l2bname(i:i) == '.' ) then
      extstart = i
      exit dotloop
   endif
enddo dotloop

ofname = l2bname(1:extstart)//'obs_seq_out'

end subroutine create_output_filename



function real_obs_sequence ( orbit, lon1, lon2, lat1, lat2, &
                             row_thin, col_thin )
!------------------------------------------------------------------------------
!  this function is to prepare data to DART sequence format
!
type(orbit_type), intent(in) :: orbit
real(r8), intent(in) :: lon1, lon2, lat1, lat2
integer,  intent(in) :: row_thin, col_thin

integer :: max_num=MAX_WVC*MAX_ROWS*2

type(obs_sequence_type) :: real_obs_sequence
type(obs_def_type)      :: obs_def
type(obs_type)          :: obs, prev_obs

integer :: i, irow, iwvc, iamb
integer, PARAMETER :: num_copies = 1
!integer, PARAMETER :: num_qc = 3
integer, PARAMETER :: num_qc = 1
integer :: days, seconds
integer :: obs_num
integer :: which_vert, uobstype, vobstype

real(r4) :: speed, dir
real(r8) :: lon, lat, vloc
real(r8) :: u_obs, v_obs, u_var, v_var
real(r8) :: aqc(num_qc)
real(r8) :: sintheta, costheta, dirvar, speedvar

type(time_type) :: time, pre_time

if ( .not. module_initialized ) call initialize_module

! Initialize an obs_sequence 
call init_obs_sequence(real_obs_sequence, num_copies, num_qc, max_num)

! set meta data of obs_seq
do i = 1, num_copies
   call set_copy_meta_data(real_obs_sequence, i, 'observation')
end do

call set_qc_meta_data(real_obs_sequence, 1, 'QC flag - wvc_quality_flag')
!call set_qc_meta_data(real_obs_sequence, 2, 'IMUDH - mp_rain_probability')
!call set_qc_meta_data(real_obs_sequence, 3, 'NOF Index - nof_rain_index')

! Initialize the obs variables
call init_obs(     obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

! assign each observation the correct observation type
uobstype = get_index_for_type_of_obs('QKSWND_U_WIND_COMPONENT')
vobstype = get_index_for_type_of_obs('QKSWND_V_WIND_COMPONENT')
if(( uobstype < 1) .or. ( vobstype < 1)) then
   msgstring = 'unknown observation type [QKSWND_U_WIND_COMPONENT] ... dying ...'
   call error_handler(E_ERR,'real_obs_sequence',msgstring,source,revision,revdate)
endif

!  loop over all observations within the file
!------------------------------------------------------------------------------

obs_num    = 0
vloc       = 10  ! all L2b data is 10m height
which_vert = VERTISHEIGHT

!rowloop:  do irow=400,403
rowloop:  do irow=1,MAX_ROWS

   ! if we're going to subset rows, cycle here
   if (row_thin > 0) then
      if (modulo(irow, row_thin) /= 1) cycle rowloop
   endif

   time = orbit%row_time(irow)
   call get_time(time, seconds, days)

   wvcloop:  do iwvc=1,MAX_WVC

      ! if we're going to subset columns, cycle here
      if (col_thin > 0) then
         if (modulo(iwvc, col_thin) /= 1) cycle wvcloop
      endif

      ! no ambiguities means no retrieval

      if (orbit%num_ambigs(iwvc,irow) < 1) cycle wvcloop

    !   only use data flagged as 'best'
    !   this should be namelist controlled at some point
    !
      if ( orbit%wvc_quality_flag(iwvc,irow) /= 0 ) cycle wvcloop

      aqc(1) = orbit%wvc_quality_flag(   iwvc,irow)
!     aqc(2) = orbit%mp_rain_probability(iwvc,irow)
!     aqc(3) = orbit%nof_rain_index(     iwvc,irow)

      lat  = orbit%wvc_lat(iwvc,irow) ! valid range [-90.00,  90.00]
      lon  = orbit%wvc_lon(iwvc,irow) ! valid range [  0.00, 359.99]

      ! verify the location is not outside valid limits
      if((lon > 360.0_r8) .or. (lon <   0.0_r8) .or.  &
         (lat >  90.0_r8) .or. (lat < -90.0_r8)) then
         write(*,*)'WARNING : invalid location.  wvc,row,lon,lat = ', iwvc,irow,lon,lat
         cycle wvcloop
      endif

      ! reject observations outside the bounding box (allowing wrapping)
      if(( lat < lat1) .or. ( lat > lat2 ) .or. &
         (.not. is_longitude_between(lon, lon1, lon2))) cycle wvcloop

      ! QuikSCAT uses the oceanographic/flow convention ... 
      ! 0.0 is TOWARD the north - in direct contradiction to 
      ! atmospheric convention. Convert by adding 180 modulo 360 
   
      iamb  = orbit%wvc_selection(     iwvc,irow)

      ! using the 'selected wind' ... aka highest-ranked ambiguity,

!     speed =     orbit%wind_speed(iamb,iwvc,irow)
!     dir   = mod(orbit%wind_dir(  iamb,iwvc,irow)+180.0_r4, 360.0_r4)

      ! Use this block for the DIRTH selection

      speed =     orbit%wind_speed_selection(iwvc,irow)
      dir   = mod(orbit%wind_dir_selection(  iwvc,irow) + 180.0_r4, 360.0_r4)

      if ( speed < 1.0_r4 ) cycle wvcloop ! everything unreliable 

      ! The requirements for QuikSCAT were 2 m/s speed (rms) over 3-20m/s
      ! or 10% of the speed from 20-30 m/s and 20 degrees (rms) 3-30 m/s

      dirvar = (20.0_r4*deg2rad)**2   ! 20 degree (rms) by spec
      speedvar = max(2.0_r4, speed*0.1_r4)**2

      sintheta = sin(dir*deg2rad)
      costheta = cos(dir*deg2rad)

      ! converting the speed and direction and uncertainties to U,V
      ! is a bit tedious for the uncertainties. Really could/should
      ! have a native type for speed and direction and a forward
      ! operator that takes U,V and calculates speed/dir ...(trivial)
      ! left for a future upgrade TJH

      ! U = R*sin(theta)
      ! V = R*cos(theta)
      ! theta ~ theta0 + etheta ... etheta ~ (0,dirvar)
      ! R     ~     R0 + espeed ... espeed ~ (0,speedvar)
      ! U ~ sin(theta0+etheta)(R0 + espeed)
      ! expand via Taylor ...
      ! U ~ sin(theta0)R0 + R0*cos(theta0)*etheta + sin(theta0)*espeed + ...
      ! The first term is the mean. 
      ! The variance of a constant times a random variable is the 
      ! constant-squared times the variance of the random variable.
      ! The random variables in the higher order terms is etheta and espeed. 

      u_obs = speed * sintheta
      v_obs = speed * costheta

      u_var = ((speed*costheta)**2)*dirvar + (sintheta**2)*speedvar
      v_var = ((speed*sintheta)**2)*dirvar + (costheta**2)*speedvar
      
      ! verify the location is not outside valid limits
      if((lon > 360.0_r8) .or. (lon <   0.0_r8) .or.  &
         (lat >  90.0_r8) .or. (lat < -90.0_r8)) then
         write(*,*) 'invalid location.  lon,lat = ', lon, lat
         cycle wvcloop
      endif

      obs_num = obs_num + 1

      ! create the obs_def for this observation, add to sequence
      !---------------------------------------------------------------
   
      call real_obs(num_copies, num_qc, obs, lon, lat, vloc, u_obs, &
                 u_var, aqc, uobstype, which_vert, seconds, days)
   
      if(obs_num == 1) then ! for the first observation 
         call insert_obs_in_seq(real_obs_sequence, obs)
         call copy_obs(prev_obs, obs)
         pre_time = time
      else  !  not the first observation 
         if(time == pre_time) then  ! same time as previous observation
            call insert_obs_in_seq(real_obs_sequence, obs, prev_obs)
            call copy_obs(prev_obs, obs)
         else  ! not the same time
            call insert_obs_in_seq(real_obs_sequence, obs)
            call copy_obs(prev_obs, obs)
            pre_time = time
         endif
      endif

      ! The (paired) V observation must be at the same time,
      ! and cannot be first, so it is easier.
      call real_obs(num_copies, num_qc, obs, lon, lat, vloc, v_obs, &
                 v_var, aqc, vobstype, which_vert, seconds, days)
      call insert_obs_in_seq(real_obs_sequence, obs, prev_obs)
      call copy_obs(prev_obs, obs)

      obs_num = obs_num + 1

      if ( DEBUG ) then
         write(*,*)irow, iwvc, aqc, lat, lon, &
            sqrt(u_obs**2+v_obs**2), atan(u_obs/v_obs)*rad2deg
      endif

   enddo wvcloop
enddo rowloop

! Print a little summary
print*, 'obs used = ', obs_num, ' obs skipped = ', 2*MAX_WVC*MAX_ROWS - obs_num

if ( get_first_obs(real_obs_sequence, obs) ) then
   call get_obs_def(obs, obs_def)
   pre_time = get_obs_def_time(obs_def)
   call print_time(pre_time,' first time in sequence is ')
   call print_date(pre_time,' first date in sequence is ')
endif
if( get_last_obs(real_obs_sequence, obs)) then
   call get_obs_def(obs, obs_def)
   time = get_obs_def_time(obs_def)
   call print_time(time,' last  time in sequence is ')
   call print_date(time,' last  date in sequence is ')
endif
print*, ''

end function real_obs_sequence



subroutine real_obs(num_copies, num_qc, obs, lon, lat, vloc, obs_value, &
                      var2, aqc, obs_kind, which_vert, seconds, days)
!------------------------------------------------------------------------------
integer,        intent(in)    :: num_copies, num_qc
type(obs_type), intent(inout) :: obs
real(r8),       intent(in)    :: lon, lat, vloc, obs_value, var2, aqc(:)
integer,        intent(in)    :: obs_kind, which_vert, seconds, days

integer            :: i
real(r8)           :: obs_value01(1)
type(obs_def_type) :: obsdef0

if ( .not. module_initialized ) call initialize_module

! Does real initialization of an observation type

call real_obs_def(obsdef0, lon, lat, vloc, &
                    var2, obs_kind, which_vert, seconds, days)
call set_obs_def(obs, obsdef0)

do i = 1, num_copies
   obs_value01(1) = obs_value
   call set_obs_values(obs, obs_value01(1:1) )
end do

call set_qc(obs, aqc)

end subroutine real_obs



subroutine real_obs_def(obs_def, lon, lat, vloc, &
                        var2, obs_kind, which_vert, seconds, days)
!----------------------------------------------------------------------
type(obs_def_type), intent(inout) :: obs_def
real(r8),intent(in) :: lon, lat, vloc, var2
integer, intent(in) :: obs_kind, which_vert, seconds, days

type(location_type) :: loc0

if ( .not. module_initialized ) call initialize_module

! set obs location
loc0 = set_location(lon, lat, vloc, which_vert )
call set_obs_def_location(obs_def, loc0)

! set obs kind
call set_obs_def_type_of_obs(obs_def, obs_kind)

call set_obs_def_time(obs_def, set_time(seconds, days) )
call set_obs_def_error_variance(obs_def, var2)

end subroutine real_obs_def



subroutine read_qscat2b(l2b_file, orbit)
!====================================================================
!
!     Description:
!
!       This file contains 3 subroutines in order to read the 
!       QuikSCAT Level 2B data in Hierarchical Data Format (HDF).  
!       The subroutines are as follows.
!
!       1. read_attrib_byname():  a subroutine to read the name 
!                                 and value(s) of a global attribute
!                                 referenced by its name.
!       
!       2. read_timetags():  a subroutine to read the timetag info
!                            contained in the HDF VDATA
!
!       3. extract_sds():  a subroutine to read the contents of an
!                          SDS from an HDF file
!
!     NOTES:
!     1. Please refer all questions concerning this program and
!        QuikSCAT data obtained from the JPL PO.DAAC to
!        qscat@podaac.jpl.nasa.gov.
!
!     2. The HDF library must be installed before this program will 
!        work properly.  The HDF library and further information 
!        about HDF may be obtained from the National Center for 
!        Supercomputing Applications (NCSA) at http://hdf.ncsa.uiuc.edu.
!
!     3. The L2B data are read in their entirety.  Examples of reading
!        the QuikSCAT data by slabs can be found in read_qscat1b.f
!        read_qscat2a.f.
!
!====================================================================

character(len=128), intent(in)  :: l2b_file
type(orbit_type),   intent(out) :: orbit

!     Set Parameters

integer :: DFACC_RDONLY
parameter (DFACC_RDONLY = 1)

!     Define Variables

! character(len=1)  :: product(8)
character(len=21) :: TimeTags(MAX_ROWS)

real(r4), dimension(MAX_WVC,MAX_ROWS) :: datmat

integer :: sd_id, retn, sfstart, sfend
integer :: irow, iwvc, iamb
! integer :: ntype, nval

integer  :: year, doy, hh, mm
real(r4) :: ss
type(time_type) :: basetime, offset

if ( .not. module_initialized ) call initialize_module

call error_handler(E_MSG,'read_qscat2b','FILENAME: '//trim(l2b_file))

if ( .not. file_exist(l2b_file)) then
   call error_handler(E_ERR,'read_qscat2b', &
           'QuikSCAT L2B file does not exist.', source, revision, revdate)
endif

! Open the HDF input file and initiate the SD interface
sd_id = sfstart(l2b_file,DFACC_RDONLY)

! Make sure that the file is a QuikSCAT Level 2B file
! I've never been able to make this work with F90 ... crazy games are afoot.
!----------------------------------------------------------------------
! call read_attrib_byname(sd_id,'ShortName',ntype,nval,product)
! write(*,*) product 
! write(*,*) product .ne. 'QSCATL2B'
!if ( product .ne. 'QSCATL2B') then
!  print *,'The input file is not a QuikSCAT Level 2B file'
!  print *,'*** Aborting program ***'
!  stop
!endif

! Read the timetag info contained in the HDF VDATA
call read_timetags(l2b_file, TimeTags)

!     Read each SDS in its entirety.  For an example of reading
!     the QuikSCAT SDS data in slabs, please refer to read_qscat2a.f.

irow=1
call extract_sds(sd_id,'wvc_row',             irow, MAX_ROWS, orbit%wvc_row)
call extract_sds(sd_id,'wvc_lat',             irow, MAX_ROWS, orbit%wvc_lat)
call extract_sds(sd_id,'wvc_lon',             irow, MAX_ROWS, orbit%wvc_lon)
call extract_sds(sd_id,'atten_corr',          irow, MAX_ROWS, orbit%atten_corr)
call extract_sds(sd_id,'model_speed',         irow, MAX_ROWS, orbit%model_speed)
call extract_sds(sd_id,'model_dir',           irow, MAX_ROWS, orbit%model_dir)
call extract_sds(sd_id,'wind_speed',          irow, MAX_ROWS, orbit%wind_speed)
call extract_sds(sd_id,'wind_dir',            irow, MAX_ROWS, orbit%wind_dir)
call extract_sds(sd_id,'wind_speed_err',      irow, MAX_ROWS, orbit%wind_speed_err)
call extract_sds(sd_id,'wind_dir_err',        irow, MAX_ROWS, orbit%wind_dir_err)
call extract_sds(sd_id,'max_likelihood_est',  irow, MAX_ROWS, orbit%max_likelihood_est)
call extract_sds(sd_id,'wind_speed_selection',irow, MAX_ROWS, orbit%wind_speed_selection)
call extract_sds(sd_id,'wind_dir_selection',  irow, MAX_ROWS, orbit%wind_dir_selection)
call extract_sds(sd_id,'mp_rain_probability', irow, MAX_ROWS, orbit%mp_rain_probability)
call extract_sds(sd_id,'srad_rain_rate',      irow, MAX_ROWS, orbit%srad_rain_rate)
call extract_sds(sd_id,'wvc_index',           irow, MAX_ROWS, datmat)
                  orbit%wvc_index = nint(datmat)
call extract_sds(sd_id,'num_in_fore',         irow, MAX_ROWS, datmat)
                  orbit%num_in_fore = nint(datmat)
call extract_sds(sd_id,'num_in_aft',          irow, MAX_ROWS, datmat)
                  orbit%num_in_aft = nint(datmat)
call extract_sds(sd_id,'num_out_fore',        irow, MAX_ROWS, datmat)
                  orbit%num_out_fore = nint(datmat)
call extract_sds(sd_id,'num_out_aft',         irow, MAX_ROWS, datmat)
                  orbit%num_out_aft = nint(datmat)
call extract_sds(sd_id,'wvc_quality_flag',    irow, MAX_ROWS, datmat)
                  orbit%wvc_quality_flag = int(datmat)
call extract_sds(sd_id,'num_ambigs',          irow, MAX_ROWS, datmat)
                  orbit%num_ambigs = nint(datmat)
call extract_sds(sd_id,'wvc_selection',       irow, MAX_ROWS, datmat)
                  orbit%wvc_selection = nint(datmat)
call extract_sds(sd_id,'nof_rain_index',      irow, MAX_ROWS, datmat)
                  orbit%nof_rain_index = nint(datmat)

! Convert time tag to a dart timetype.
!00000000011111111112
!12345678901234567890
!2008-001T01:16:24.584

do irow = 1,MAX_ROWS
   read(TimeTags(irow),'(      i4)')year
   read(TimeTags(irow),'( 5x,  i3)')doy
   read(TimeTags(irow),'( 9x,  i2)')hh
   read(TimeTags(irow),'(12x,  i2)')mm
   read(TimeTags(irow),'(15x,f7.3)')ss

   basetime = set_date(year,1,1,0,0,0) ! start of the (right) year
   offset   = set_time((hh*60 + mm)*60 + nint(ss), doy-1)
   orbit%row_time(irow) = basetime + offset 

!  if (DEBUG) call print_date(orbit%row_time(irow), TimeTags(irow))
enddo

! Print results to screen
if ( DEBUG ) then
RECORDLOOP : do irow=400,410  ! or 1,MAX_ROWS

   write(*,*) ' '
   write(*,*) 'TIME: ', TimeTags(irow)
   write(*,'(''WVC ROW: '',f5.0)') orbit%wvc_row(irow)

   write(*,105)
 105    format('WVC#','  WVC_Qual',        &
               '  WVC Latitude/Longitude', &
               '  Selected Wind Vector',   &
               ' NWP Wind Vector',         &
               '  Num/Sel ambig',          &
               '  DRE Wind Vector',        &
               ' MUDH Prob',               &
               ' NOF Index',               &
               ' SRR')

   WVCLOOP : do iwvc = 1,MAX_WVC
      if (orbit%num_ambigs(iwvc,irow).gt.0) then
         iamb=orbit%wvc_selection(iwvc,irow)
         write(*,110)  orbit%wvc_index(iwvc,irow), &    ! WVC#
            int(orbit%wvc_quality_flag(iwvc,irow)), &   ! WVC_Qual
            orbit%wvc_lat(iwvc,irow), orbit%wvc_lon(iwvc,irow), &
            orbit%wind_speed(iamb,iwvc,irow), orbit%wind_dir(iamb,iwvc,irow), & ! Selected Wind Vector
            orbit%model_speed(iwvc,irow), orbit%model_dir(iwvc,irow), &         ! NWP Wind Vector
            orbit%num_ambigs(iwvc,irow),  orbit%wvc_selection(iwvc,irow), &      ! Num/Sel ambig
            orbit%wind_speed_selection(iwvc,irow), orbit%wind_dir_selection(iwvc,irow), & ! DRE Wind Vector
            orbit%mp_rain_probability(iwvc,irow), &   ! MUDH Prob
            orbit%nof_rain_index(iwvc,irow), &        ! NOF Index
            orbit%srad_rain_rate(iwvc,irow)           ! SRR
      endif
   enddo WVCLOOP
enddo RECORDLOOP
endif

 110  format(i3,4x,"0X",z4.4,8x,f6.2,3x,f6.2,6x,f6.2,3x,f6.2, &
             4x,f6.2,3x,f6.2,7x,i2,3x,i2,3x,f6.2,3x,f6.2, &
             3x,f6.2,3x,i6,2x,f6.2)

retn = sfend(sd_id)
end subroutine read_qscat2b

!====================================================================
!    READ_ATTRIB_BYNAME:  a subroutine to read the name and
!                         value(s) of a global attribute
!                         referenced by its name.
!    
!    5/14/1998 R.S. Dunbar
!====================================================================

      subroutine read_attrib_byname(sd_id, in_attr_name, &
                            num_type, n_values, fvalues)
      
      character*(*) fvalues(*)

      integer :: sd_id,num_type,n_values
      integer :: attr_index,count,retn,n,oldn
      integer :: sffattr,sfgainfo,sfrattr
      integer :: ival, i

      integer :: MAX_NC_NAME
      parameter (MAX_NC_NAME=256)

      character*(*) in_attr_name
      character attr_name*(MAX_NC_NAME),attr_data*512
      character*(MAX_NC_NAME) values(20)
      character :: cr
 
      if ( .not. module_initialized ) call initialize_module

!     Find the attribute assigned to in_attr_name
      attr_index = sffattr(sd_id,in_attr_name)

!     Get information about the  file attribute
      retn = sfgainfo(sd_id,attr_index,attr_name,num_type,count)

!     Read the attribute data
      retn = sfrattr(sd_id,attr_index,attr_data)

      cr = char(10)
      ival = 0
      oldn = 1
 5    continue

!     QuikSCAT attributes have atleast three lines: 
!     metadata type, array size and metadata contents
!     Use "blank spaces" to identify the end of a line

      n = index(attr_data(oldn:(count-1)),cr)

!     Read all of the metadata lines
      if (n .eq. 0) then
         ival=ival+1
         values(ival) = attr_data(oldn:(count-1))
         goto 99
      else
         ival=ival+1
         values(ival) = attr_data(oldn:(oldn+n-2))
      endif
      oldn=n+oldn
      goto 5

 99   continue
      n_values = ival - 2
      do i=1,n_values
         fvalues(i) = values(i+2)
      enddo
      return
      end subroutine read_attrib_byname

!====================================================================
!    READ_TIMETAGS:  a subroutine to read the timetag info
!                    contained in the HDF VDATA
!    
!    5/1998 R.S. Dunbar
!
!    Revisions:
!    7/1999 Code adapted to read timetags in their entirety.
!           Commenter were also added.  K.L. Perry
!====================================================================
      subroutine read_timetags(filename,timetags)

      character(len=80) :: filename
      character(len=21) :: timetags(*)
      character(len=60) :: fields
      character vdata_name*30
      integer :: file_id,vdata_ref,vdata_id
      integer :: n_records,interlace,vdata_size
      integer :: hopen,vsfgid,vsfatch,vsfinq,vsfread,hclose
      integer :: retn

      integer :: DFACC_RDONLY,FULL_INTERLACE
      parameter(DFACC_RDONLY=1)
      parameter(FULL_INTERLACE=0)

      if ( .not. module_initialized ) call initialize_module

!     Open the HDF file
      file_id = hopen(filename,DFACC_RDONLY,0)

!     Initialize the VS interface
      call vfstart(file_id)

!     Get the reference number for the first vdata in the file
      vdata_ref = -1
      vdata_ref = vsfgid(file_id,vdata_ref)

!     Attach to the vdata for reading if it is found, otherwise 
!     exit the program.
      if (vdata_ref.eq.0) then
         print *,'No Timetags were found in the HDF VDATA'
         print *,'*** Aborting program ***'
         stop
      endif

      vdata_id = vsfatch(file_id,vdata_ref,'r')

!     Get n_records
      retn=vsfinq(vdata_id,n_records,interlace,fields,vdata_size,vdata_name)

!     Read the timetags
      retn = vsfread(vdata_id,timetags,n_records,FULL_INTERLACE)

!     Terminate access to the vdata and to the VS interface, 
!     then close the HDF file.

!     retn =  vsfdtch(vdata_id)
      call vsfdtch(vdata_id)
      call    vfend(file_id)
      retn = hclose(file_id)

      return
      end subroutine read_timetags



!====================================================================
!    EXTRACT_SDS:  a subroutine to read the contents of an
!                  SDS from an HDF file
!    
!    5/12/1998 R.S. Dunbar
!
!    Revisions:
!    7/1999   Code adapted to read input in bytes as well as ints 
!             and floats.  Comments were also added.  K.L. Perry
!
!    3/2000   Corrected code for 8-bit unsigned integers.  "buffer"
!             was used instead of "buffer2".  K.L. Perry
!
!    5/2000   Changed MAX_BUF_SIZE from 1000000 to 10000000.
!             Corrected code for 32-bit unsigned integers.  Created
!             "buffer3" array of int*4 to correctly read in uint32.
!             K.L. Perry
!
!====================================================================
      subroutine extract_sds(sd_id,in_var,irec,slab_size,out_var)

      integer :: MAX_BUF_SIZE
      parameter (MAX_BUF_SIZE=10000000)
      integer :: sd_id,sds_index,sds_id,retn
      integer :: rank,dim_sizes(3),data_type,nattrs,num_type
      integer :: edge(3),stride(3),start(3),irec,slab_size
      real(digits12) :: cal, cal_err, off, off_err
      integer :: iprod,i,itmp

      character*(*) in_var
      character(len=256) :: name
      integer :: sfn2index,sfselect,sfginfo,sfrdata,sfgcal,sfendacc

      integer*2 :: buffer(MAX_BUF_SIZE)
      byte buffer2(MAX_BUF_SIZE)
      integer*4 :: buffer3(MAX_BUF_SIZE)
      real(r4) :: out_var(*)

      if ( .not. module_initialized ) call initialize_module

!     Search for the index of "in_var"
      sds_index = sfn2index(sd_id, in_var)

!     Select data set corresponding to the returned index
      sds_id = sfselect(sd_id,sds_index)
      retn = sfginfo(sds_id,name,rank,dim_sizes,data_type,nattrs)

      do i=1,rank
         edge(i)=dim_sizes(i)
         start(i)=0
         stride(i)=1
      enddo
      edge(rank)=slab_size
      start(rank)=irec-1

      iprod=1
      do i=1,rank
         iprod=iprod*edge(i)
      enddo

!     Get the calibration and offset values of input
      retn = sfgcal(sds_id,cal,cal_err,off,off_err,num_type)

!     Read Arrays which are not float32 or int8 or uint8 or uint32
      if ((data_type.ne. 5).and.(data_type.ne.20).and. &
          (data_type.ne.21).and.(data_type.ne.25)) then

!     Read the data set into the "buffer" array
         retn=sfrdata(sds_id,start,stride,edge,buffer)

!     Calibrate the output
         do i=1,iprod
!     Correct for 16-bit unsigned integers
            if ((data_type.eq.23).and.(buffer(i).lt.0)) then
               out_var(i)=buffer(i)+65536.0

!     No correction needed for signed or positive unsigned integers
            else
               out_var(i)=buffer(i)
            endif

            out_var(i)=out_var(i)*cal
         enddo

!     Read uint32 arrays.
      else if (data_type.eq.25) then

!     Read the data set into the "buffer3" uint32 array
         retn=sfrdata(sds_id,start,stride,edge,buffer3)

!     Calibrate the output
         do i=1,iprod
!     Correct for 32-bit unsigned integers
            if ((data_type.eq.25).and.(buffer3(i).lt.0)) then
               out_var(i)=buffer3(i)+4294967296.0
            else
               out_var(i)=buffer3(i)
            endif
            out_var(i)=out_var(i)*cal
         enddo

!     Read int8 and uint8 arrays. 
      else if ((data_type.eq.20).or.(data_type.eq.21)) then

!     Read the data set into the "buffer2" byte-array
         retn=sfrdata(sds_id,start,stride,edge,buffer2)

!     Calibrate the output
         do i=1,iprod

!     Correct for 8-bit unsigned integers
            itmp=buffer2(i)
            if ((data_type.eq.21).and.(buffer2(i).lt.0)) then
               itmp=itmp+256
               out_var(i)=itmp
               
!     No correction needed for signed or positive unsigned integers
            else
               out_var(i)=itmp
            endif

            out_var(i)=out_var(i)*cal
         enddo

      else
!     Read float32 arrays directly into the "out_var" array
         retn=sfrdata(sds_id,start,stride,edge,out_var)

!     Calibrate the output
         do i=1,iprod
            out_var(i)=out_var(i)*cal
         enddo
      endif

!     Terminate access to the array data set.
      retn = sfendacc(sds_id)
      end subroutine extract_sds

end module quikscat_JPL_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
