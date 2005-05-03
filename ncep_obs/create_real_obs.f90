! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program create_real_obs

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

use obs_sequence_mod, only : obs_sequence_type, write_obs_seq, &
                             static_init_obs_sequence, destroy_obs_sequence 
use     real_obs_mod, only : real_obs_sequence
use    utilities_mod, only : get_unit, open_file, close_file, file_exist, &
                             check_nml_error, initialize_utilities, register_module, &
                             error_handler, timestamp, E_MSG, logfileunit

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(obs_sequence_type) :: seq

character(len = 129) :: output_name
character(len = 8 ) :: obsdate
integer :: iunit, ierr, io, ii, day1

! ----------------------------------------------------------------------
! Declare namelist parameters
! ----------------------------------------------------------------------
        
integer :: year = 2003, month =1, day =1, tot_days = 31
integer :: max_num = 800000, select_obs = 0
character(len = 129) :: ObsBase = '/project/dart/home/hliu/ncepobs/April05/temp_obs.'
logical :: ADDUPA = .false., AIRCAR = .false., AIRCFT = .false., SATWND = .false.
logical :: obs_U  = .false., obs_V  = .false., obs_T  = .false. , &
           obs_PS = .false., obs_QV = .false.

namelist /ncepobs_nml/year, month, day, tot_days, max_num, select_obs, ObsBase, &
        ADDUPA, AIRCAR, AIRCFT, SATWND, obs_U, obs_V, obs_T, obs_PS, obs_QV

! ----------------------------------------------------------------------
! Select observation types using NCEP categories (when select_obs \= 0).
!
!  ADPUPA: upper-air reports (mostly radiosonde plus few dropsonde, PIBAL)
!  AIRCFT: Conv. (AIREP, PIREP) and ASDAR aircraft reports
!  AIRCAR: ACARS sircraft reports
!  SATWND: Satellite derived wind reports
! ----------------------------------------------------------------------

!  Select variables of U, V, T, QV, PS using the logicals:
!  obs_U   obs_V   obs_PS   obs_T   obs_QV  

call initialize_utilities('create_real_obs')
call register_module(source,revision,revdate)
call error_handler(E_MSG,'create_real_obs','STARTING',source,revision,revdate)

! Initialize the obs_sequence module ...

call static_init_obs_sequence()

! Read namelist ...

if(file_exist('input.nml')) then
  iunit = open_file('input.nml', action = 'read')
  ierr = 1

  do while(ierr /= 0)
     read(iunit, nml = ncepobs_nml, iostat = io, end = 11)
     ierr = check_nml_error(io, 'ncepobs_nml')
   enddo
11 continue
   call close_file(iunit)
endif
write(logfileunit,nml=ncepobs_nml)

! Loop through the days interested.

do ii = 1, tot_days
 
   day1 = day + ii -1

   ! define observation filename

   write(obsdate, '(i4.4,i2.2,i2.2)') year, month, day1

   ! set the obs sequence of the day

   seq = real_obs_sequence(year, month, day1, max_num, select_obs, ObsBase, &
        ADDUPA, AIRCAR, AIRCFT, SATWND, obs_U, obs_V, obs_T, obs_PS, obs_QV)

   ! output the daily sequence to a file

   output_name = 'obs_seq'//obsdate
   call write_obs_seq(seq, output_name)

   ! release the memory of the seq.

   call destroy_obs_sequence(seq)

enddo

call timestamp(source,revision,revdate,'end') ! close the log file.

end program create_real_obs
