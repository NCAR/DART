! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program convert_gmi_L1

! Program to read GPM/GMI radiances

use types_mod,        only : r8, deg2rad, PI
use obs_sequence_mod, only : obs_sequence_type, write_obs_seq, &
                             static_init_obs_sequence, destroy_obs_sequence, &
                             init_obs_sequence
use       gmi_L1_mod, only : gmi_l1c_granule_type, gmi_load_l1c_granule, &
                             make_obs_sequence, MAXCHANS
use    utilities_mod, only : initialize_utilities, register_module, &
                             error_handler, finalize_utilities, E_ERR, E_MSG, &
                             find_namelist_in_file, check_namelist_read, &
                             do_nml_file, do_nml_term, set_filename_list, &
                             nmlfileunit

implicit none

! ----------------------------------------------------------------------
! Declare local parameters
! ----------------------------------------------------------------------

integer                 :: filecount
type(gmi_l1c_granule_type) :: granule
type(obs_sequence_type) :: seq

integer :: io, iunit, index, nchans

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'convert_gmi_L1.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

! ----------------------------------------------------------------------
! Declare namelist parameters
! ----------------------------------------------------------------------

integer, parameter :: MAXFILES = 512

character(len=256) :: l1_files(MAXFILES) = ''
character(len=256) :: l1_file_list       = ''
character(len=256) :: outputfile         = ''

real(r8) :: lon1 =   0.0_r8,  &   !  lower longitude bound
            lon2 = 360.0_r8,  &   !  upper longitude bound
            lat1 = -90.0_r8,  &   !  lower latitude bound
            lat2 =  90.0_r8       !  upper latitude bound

integer  :: cross_track_thin = 0
integer  :: along_track_thin = 0

integer  :: i

logical  :: use_channels(MAXCHANS)

! specify exactly which channel numbers to convert
integer  :: channel_list(MAXCHANS) = -1

namelist /convert_gmi_L1_nml/ l1_files, l1_file_list, &
                              outputfile, &
                              lon1, lon2, lat1, lat2, &
                              channel_list, &
                              cross_track_thin, along_track_thin

! ----------------------------------------------------------------------
! start of executable program code
! ----------------------------------------------------------------------

call initialize_utilities('convert_gmi_L1')
call register_module(source,revision,revdate)

! Initialize the obs_sequence module ...

call static_init_obs_sequence()

!----------------------------------------------------------------------
! Read the namelist
!----------------------------------------------------------------------

call find_namelist_in_file('input.nml', 'convert_gmi_L1_nml', iunit)
read(iunit, nml = convert_gmi_L1_nml, iostat = io)
call check_namelist_read(iunit, io, 'convert_gmi_L1_nml')

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=convert_gmi_L1_nml)
if (do_nml_term()) write(    *      , nml=convert_gmi_L1_nml)


! when this routine returns, the l1_files variable will have
! all the filenames, regardless of which way they were specified.
filecount = set_filename_list(l1_files, l1_file_list, "convert_gmi_L1")

! If no channels are specified in the namelist, 
! all the channels will be output.

nchans = 0
do while (channel_list(nchans+1) > 0 .and. nchans < MAXCHANS)
   nchans = nchans + 1
enddo

if (nchans > 1) then
   use_channels(:) = .false.

   do i=1,nchans
      use_channels(channel_list(i)) = .true.
   end do
else
   use_channels(:) = .true.
end if

! for each input file
do index=1, filecount

   ! read from HDF file into a derived type that holds all the information
   call gmi_load_l1c_granule(l1_files(index), granule)

   ! convert derived type information to DART sequence
   call make_obs_sequence(seq, granule, lon1, lon2, lat1, lat2, &
                          use_channels, along_track_thin,       &
                          cross_track_thin)

   ! write the sequence to a disk file
   call write_obs_seq(seq, outputfile)

   ! release the sequence memory
   call destroy_obs_sequence(seq)
enddo

call error_handler(E_MSG, source, 'Finished successfully.',source,revision,revdate)
call finalize_utilities()

end program convert_gmi_L1

