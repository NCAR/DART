! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program convert_airs_L2

! <next few lines under version control, do not edit>
! $URL: https://subversion.ucar.edu/DAReS/DART/trunk/observations/AIRS/convert_airs_L2.f90 $
! $Id: convert_airs_L2.f90 3809 2009-04-13 16:21:33Z nancy $
! $Revision: 3809 $
! $Date: 2009-04-13 10:21:33 -0600 (Mon, 13 Apr 2009) $

! Initial version of a program to read the AIRS retrievals for temperature
! and humidity. 

use types_mod,        only : r8, deg2rad, PI
use obs_sequence_mod, only : obs_sequence_type, write_obs_seq, &
                             static_init_obs_sequence, destroy_obs_sequence
use     airs_JPL_mod, only : airs_ret_rdr, airs_granule_type
use    utilities_mod, only : initialize_utilities, register_module, &
                             error_handler, timestamp, E_ERR, E_MSG, &
                             find_namelist_in_file, check_namelist_read, &
                             do_nml_file, do_nml_term, &
                             logfileunit, nmlfileunit

use airs_obs_mod   ! FIXME: need to add ,only :   

implicit none

! ----------------------------------------------------------------------
! Declare local parameters
! ----------------------------------------------------------------------

character(len=256)      :: datafile(1), output_name, dartfile, string1
type(airs_granule_type) :: granule
type(obs_sequence_type) :: seq

integer :: io, iunit

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL: https://subversion.ucar.edu/DAReS/DART/trunk/observations/AIRS/convert_airs_L2.f90 $", &
   revision = "$Revision: 3809 $", &
   revdate  = "$Date: 2009-04-13 10:21:33 -0600 (Mon, 13 Apr 2009) $"

! ----------------------------------------------------------------------
! Declare namelist parameters
! ----------------------------------------------------------------------
        
character(len=128) ::   l2_file = ''
character(len=128) ::   datadir = '.'
character(len=128) :: outputdir = '.'

real(r8) :: lon1 =   0.0_r8,  &   !  lower longitude bound
            lon2 = 360.0_r8,  &   !  upper longitude bound 
            lat1 = -90.0_r8,  &   !  lower latitude bound
            lat2 =  90.0_r8       !  upper latitude bound

namelist /convert_airs_L2_nml/ l2_file, datadir, outputdir, &
                           lon1, lon2, lat1, lat2

! ----------------------------------------------------------------------
! start of executable program code
! ----------------------------------------------------------------------

call initialize_utilities('convert_airs_L2')
call register_module(source,revision,revdate)

! Initialize the obs_sequence module ...

call static_init_obs_sequence()

!----------------------------------------------------------------------
! Read the namelist
!----------------------------------------------------------------------

call find_namelist_in_file('input.nml', 'convert_airs_L2_nml', iunit)
read(iunit, nml = convert_airs_L2_nml, iostat = io)
call check_namelist_read(iunit, io, 'convert_airs_L2_nml')

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=convert_airs_L2_nml)
if (do_nml_term()) write(    *      , nml=convert_airs_L2_nml)

call create_output_filename(l2_file, output_name)
datafile(1) = trim(datadir) // '/' // trim(l2_file)
dartfile = trim(outputdir) // '/' // trim(output_name)
!dartfile = trim(outputdir) // '/test.out'

call airs_ret_rdr(datafile, granule)   ! read from HDF file into a structure
seq = real_obs_sequence(granule, lon1, lon2, lat1, lat2) ! convert structure to a sequence
call write_obs_seq(seq, dartfile)
call destroy_obs_sequence(seq)       ! release the memory of the seq
call timestamp(source,revision,revdate,'end') ! close the log file

end program convert_airs_L2

