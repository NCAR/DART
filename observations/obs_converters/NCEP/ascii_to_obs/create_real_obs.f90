! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program create_real_obs

use types_mod,        only : r8, deg2rad, PI
use obs_sequence_mod, only : obs_sequence_type, write_obs_seq, &
                             static_init_obs_sequence, destroy_obs_sequence 
use     real_obs_mod, only : real_obs_sequence
use    utilities_mod, only : initialize_utilities, register_module,            &
                             do_output, nmlfileunit, do_nml_file, do_nml_term, &
                             error_handler, finalize_utilities, E_ERR, E_MSG,  &
                             find_namelist_in_file, check_namelist_read

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type(obs_sequence_type) :: seq

character(len = 256) :: output_name
character(len = 8 ) :: obsdate
integer :: iunit, io, ii, day1, kkk, kbeg, kend

character(len = 2) :: obstime(4), hour1
data obstime/'06','12','18','24'/

real(r8) :: bin_beg(5), bin_end(5)
data bin_beg/ 3.001_r8,  9.001_r8, 15.001_r8, 21.001_r8,  3.001_r8/
data bin_end/ 9.000_r8, 15.000_r8, 21.000_r8, 27.000_r8, 27.000_r8/

! ----------------------------------------------------------------------
! Declare namelist parameters
! ----------------------------------------------------------------------
        
integer :: year = 2003, month =1, day =1, tot_days = 31
integer :: max_num = 800000, select_obs = 0
character(len = 256) :: ObsBase = 'temp_obs.'
logical :: ADPUPA = .false., AIRCAR = .false., AIRCFT = .false., &
           SATWND = .false., SATEMP = .false., SFCSHP = .false., &
           ADPSFC = .false.

logical :: obs_U  = .false., obs_V  = .false., obs_T  = .false. , &
           obs_PS = .false., obs_QV = .false., daily_file = .true., & 
           obs_time = .true., obs_Z = .false.

real(r8) :: lon1 =   0.0_r8,  &   !  lower longitude bound
            lon2 = 360.0_r8,  &   !  upper longitude bound 
            lat1 = -90.0_r8,  &   !  lower latitude bound
            lat2 =  90.0_r8       !  upper latitude bound

logical  :: include_specific_humidity = .true.,  &
            include_relative_humidity = .false., &
            include_dewpoint          = .false., &
            include_surface_pressure  = .false.

namelist /ncepobs_nml/ year, month, day, tot_days, max_num, select_obs,  &
        ObsBase, ADPUPA, AIRCAR, AIRCFT, SATEMP, SFCSHP, ADPSFC, SATWND, &
        obs_U, obs_V, obs_T, obs_PS, obs_QV, obs_Z, daily_file, lon1, lon2, & 
        lat1, lat2, obs_time, include_specific_humidity, &
        include_relative_humidity, include_dewpoint, include_surface_pressure

! ----------------------------------------------------------------------
! Select observation types using NCEP categories (when select_obs /= 0).
!  ADPUPA: upper-air reports (mostly radiosonde plus few dropsonde, PIBAL)
!  AIRCFT: Conv. (AIREP, PIREP) and ASDAR aircraft reports
!  AIRCAR: ACARS sircraft reports
!  SATEMP: ATOVS retrived temperature
!  SFCSHP: SURFACE MARINE reports
!  ADPSFC: SURFACE LAND SYNOPTIC STATION reports
!  SATWND: Satellite derived wind reports
! ----------------------------------------------------------------------
! Select variables of U, V, T, QV, PS using the logicals:
!  obs_U   obs_V   obs_PS   obs_T   obs_QV    obs_Z
! ----------------------------------------------------------------------

! start of executable program code

call initialize_utilities('create_real_obs')
call register_module(source,revision,revdate)

! Initialize the obs_sequence module ...

call static_init_obs_sequence()

! Read the namelist entry
call find_namelist_in_file("input.nml", "ncepobs_nml", iunit)
read(iunit, nml = ncepobs_nml, iostat = io)
call check_namelist_read(iunit, io, "ncepobs_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=ncepobs_nml)
if (do_nml_term()) write(     *     , nml=ncepobs_nml)

! Loop through the days interested.

do ii = 1, tot_days
 
  day1 = day + ii - 1

  ! define observation filename
  write(obsdate, '(i4.4,i2.2,i2.2)') year, month, day1

  ! set the obs sequence of the day (daily or 6 hourly)
  if(daily_file) then
    kbeg = 5
    kend = 5
    output_name = 'obs_seq'//obsdate
  else
    kbeg = 1
    kend = 4
  endif

  do kkk = kbeg, kend

    if (daily_file) then
      hour1 = ''
    else
      hour1 = obstime(kkk)
    end if

    seq = real_obs_sequence(year, month, day1, hour1, max_num, select_obs, &
         ObsBase, ADPUPA, AIRCAR, AIRCFT, SATEMP, SFCSHP, ADPSFC, SATWND, &
         obs_U, obs_V, obs_T, obs_PS, obs_QV, obs_Z, include_specific_humidity, &
         include_relative_humidity, include_dewpoint, include_surface_pressure, &
         bin_beg(kkk), bin_end(kkk), lon1, lon2, lat1, lat2, obs_time)

    ! output the daily sequence to a file
    if(.not. daily_file) output_name = 'obs_seq'//obsdate//obstime(kkk)
    call write_obs_seq(seq, output_name)

    ! release the memory of the seq.
    call destroy_obs_sequence(seq)

  enddo

enddo

call error_handler(E_MSG,'create_real_obs','Finished successfully.',source,revision,revdate)
call finalize_utilities()

end program create_real_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
