! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program convert_amsu_L1

! Program to convert the AMSU/A 'AIRABRAD' brightness temperatures 
! from netCDF to a DART observation sequence file.
!
! See the REAMDE.md in this directory for more informaiton, but the
! data citation information for the dataset is 
!
! Title: AIRS/Aqua L1B AMSU (A1/A2) geolocated and calibrated brightness temperatures V005
! Version: 005
! Creator: AIRS project
! Publisher: Goddard Earth Sciences Data and Information Services Center (GES DISC)
! Release Date: 2007-07-26T00:00:00.000Z
! Linkage: https://disc.gsfc.nasa.gov/datacollection/AIRABRAD_005.html
!
! The data are originally distributed in HDF-EOS format and can be converted to netCDF
! by the h4tonccf_nc4 program available from http://hdfeos.org/software/h4cflib.php

use         types_mod, only : r8, deg2rad, PI

use  obs_sequence_mod, only : obs_sequence_type, write_obs_seq, &
                              static_init_obs_sequence, destroy_obs_sequence

use    utilities_mod, only : initialize_utilities, register_module, &
                             error_handler, finalize_utilities, E_ERR, E_MSG, &
                             find_namelist_in_file, check_namelist_read, &
                             do_nml_file, do_nml_term, set_filename_list, &
                             nmlfileunit, get_next_filename

use     amsua_bt_mod, only : amsua_bt_granule

use amsua_netCDF_support_mod, only : initialize_amsua_netcdf, &
                                     channel_list_to_indices, &
                                     read_amsua_bt_netCDF_granule, &
                                     make_obs_sequence

implicit none

! ----------------------------------------------------------------------
! Declare namelist parameters
! ----------------------------------------------------------------------
        
integer, parameter :: MAXFILES = 512
integer, parameter :: AMSUA_BT_CHANNEL = 15

character(len=256) :: l1_files(MAXFILES) = ''
character(len=256) :: l1_file_list       = ''
character(len=256) :: outputfile         = ''

real(r8) :: lon1 =   0.0_r8,  &   !  lower longitude bound
            lon2 = 360.0_r8,  &   !  upper longitude bound 
            lat1 = -90.0_r8,  &   !  lower latitude bound
            lat2 =  90.0_r8       !  upper latitude bound

character(len=8) :: channel_list(AMSUA_BT_CHANNEL) = 'null'
integer  :: cross_track_thin = 0
integer  :: along_track_thin = 0

namelist /convert_amsu_L1_nml/ l1_files, l1_file_list, &
                               outputfile, &
                               lon1, lon2, lat1, lat2, &
                               channel_list, &
                               cross_track_thin, &
                               along_track_thin

! ----------------------------------------------------------------------
! Declare local parameters
! ----------------------------------------------------------------------

integer                  :: io, iunit, ifile
integer                  :: filecount
type(amsua_bt_granule)   :: granule
type(obs_sequence_type)  :: seq
logical                  :: use_channels(AMSUA_BT_CHANNEL) = .false.

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'convert_amsu_L1.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

! ----------------------------------------------------------------------
! start of executable program code
! ----------------------------------------------------------------------

call initialize_utilities('convert_amsu_L1')
call register_module(source,revision,revdate)
call static_init_obs_sequence()

!----------------------------------------------------------------------
! Read the namelist
!----------------------------------------------------------------------

call find_namelist_in_file('input.nml', 'convert_amsu_L1_nml', iunit)
read(iunit, nml = convert_amsu_L1_nml, iostat = io)
call check_namelist_read(iunit, io, 'convert_amsu_L1_nml')

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=convert_amsu_L1_nml)
if (do_nml_term()) write(    *      , nml=convert_amsu_L1_nml)

! when this routine returns, the l1_files variable will have
! all the filenames, regardless of which way they were specified.
filecount = set_filename_list(l1_files, l1_file_list, "convert_amsu_L1")

! use the first file to initialize the module
call initialize_amsua_netcdf(l1_files(1))

! FIXME ... maybe ... are all the files guaranteed to have the
! same channels - I think so.
call channel_list_to_indices(channel_list, use_channels)

! for each input file
do ifile = 1,filecount

   ! read from HDF file into a derived type that holds all the information
   call read_amsua_bt_netCDF_granule(l1_files(ifile), granule)   

   ! convert derived type information to DART sequence
   call make_obs_sequence(seq, granule, lon1, lon2, lat1, lat2, &
                          use_channels, along_track_thin, cross_track_thin)

   ! write the sequence to a disk file
   call write_obs_seq(seq, outputfile) 
 
   ! release the sequence memory
   call destroy_obs_sequence(seq)

enddo

call error_handler(E_MSG, source, 'Finished successfully.',source,revision,revdate)
call finalize_utilities()

end program convert_amsu_L1

