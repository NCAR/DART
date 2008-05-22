! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program create_ocean_obs

! <next few lines under version control, do not edit>
! $URL: http://subversion.ucar.edu/DAReS/DART/trunk/ncep_obs/create_ocean_obs.f90 $
! $Id: create_ocean_obs.f90 3304 2008-04-09 22:20:02Z nancy $
! $Revision: 3304 $
! $Date: 2008-04-09 16:20:02 -0600 (Wed, 09 Apr 2008) $

! Initial program to read the raw ocean observations and insert them
! into an observation sequence. To make things easy ... we will mandate
! an assimilatin interval of 1 day - so all observations will be
! redefined to occur at NOON on the day they were observed.

use types_mod,        only : r8, deg2rad, PI
use obs_sequence_mod, only : obs_sequence_type, write_obs_seq, &
                             static_init_obs_sequence, destroy_obs_sequence 
use    ocean_obs_mod, only : real_obs_sequence
use    utilities_mod, only : initialize_utilities, register_module, &
                             do_output, logfileunit, &
                             error_handler, timestamp, E_ERR, E_MSG, &
                             find_namelist_in_file, check_namelist_read

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL: http://subversion.ucar.edu/DAReS/DART/trunk/ncep_obs/create_ocean_obs.f90 $", &
   revision = "$Revision: 3304 $", &
   revdate  = "$Date: 2008-04-09 16:20:02 -0600 (Wed, 09 Apr 2008) $"

type(obs_sequence_type) :: seq

character(len = 129) :: output_name
character(len = 8 ) :: obsdate
integer :: iunit, io, ii, day1, kkk, kbeg, kend

! ----------------------------------------------------------------------
! Declare namelist parameters
! ----------------------------------------------------------------------
        
integer :: year = 1996, month =1, day =1, tot_days = 31
integer :: max_num = 800000
character(len = 129) :: fname = 'raw_ocean_obs.txt'

logical :: ADCP = .false., DRIFTERS = .false., GLIDERS = .false., &
           TMI = .false., SSH = .false.

real(r8) :: lon1 =   0.0_r8,  &   !  lower longitude bound
            lon2 = 360.0_r8,  &   !  upper longitude bound 
            lat1 = -90.0_r8,  &   !  lower latitude bound
            lat2 =  90.0_r8       !  upper latitude bound

namelist /ocean_obs_nml/ year, month, day, tot_days, max_num, &
        fname, ADCP, DRIFTERS, GLIDERS, TMI, SSH, &
        lon1, lon2, lat1, lat2

! ----------------------------------------------------------------------
!  ADPUPA: upper-air reports (mostly radiosonde plus few dropsonde, PIBAL)
!  AIRCFT: Conv. (AIREP, PIREP) and ASDAR aircraft reports
!  AIRCAR: ACARS sircraft reports
!  SATEMP: ATOVS retrived temperature
!  SFCSHP: SURFACE MARINE reports
!  ADPSFC: SURFACE LAND SYNOPTIC STATION reports
!  SATWND: Satellite derived wind reports
! ----------------------------------------------------------------------

! start of executable program code

call initialize_utilities('create_ocean_obs')
call register_module(source,revision,revdate)

! Initialize the obs_sequence module ...

call static_init_obs_sequence()

! Read the namelist entry
call find_namelist_in_file("input.nml", "ocean_obs_nml", iunit)
read(iunit, nml = ocean_obs_nml, iostat = io)
call check_namelist_read(iunit, io, "ocean_obs_nml")

! Record the namelist values used for the run ...
call error_handler(E_MSG,'create_ocean_obs','ocean_obs_nml values are',' ',' ',' ')
if (do_output()) write(logfileunit, nml=ocean_obs_nml)
if (do_output()) write(     *     , nml=ocean_obs_nml)

lon1 = min(max(lon1,0.0_r8),360.0_r8)
lon2 = min(max(lon2,0.0_r8),360.0_r8)
if ( lon1 > lon2 ) lon2 = lon2 + 360.0_r8

  ! define observation filename
  write(obsdate, '(i4.4,i2.2)') year, month

  seq = real_obs_sequence(fname, year, month, day1, max_num, &
        ADCP, DRIFTERS, GLIDERS, TMI, SSH, &
        lon1, lon2, lat1, lat2)

  output_name = 'obs_seq'//obsdate
  call write_obs_seq(seq, output_name)

  ! release the memory of the seq.
  call destroy_obs_sequence(seq)

call timestamp(source,revision,revdate,'end') ! close the log file.

end program create_ocean_obs

