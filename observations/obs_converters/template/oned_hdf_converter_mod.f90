! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

program oned_hdf_converter_mod

!> title = "Generalized 1D Temperature HDF Observation Converter"
!> institution = " NCAR" ;
!> source = "NCAR/DAReS" ;
!> comment = "Generalized converter for 1D temperature data from NetCDF files" ;
!> references = "http://www.image.ucar.edu/DAReS/DART/DART_download" ;
!> dataset_title = "Generalized 1D Temperature HDF Data" ;

use types_mod, only : r8

use obs_sequence_mod, only : obs_sequence_type, write_obs_seq, &
                             static_init_obs_sequence, destroy_obs_sequence

use utilities_mod, only : initialize_utilities, register_module, &
                             error_handler, finalize_utilities, E_ERR, E_MSG, &
                             find_namelist_in_file, check_namelist_read, &
                             do_nml_file, do_nml_term, set_filename_list, &
                             nmlfileunit, get_next_filename

use hdf_data_mod, only : hdf_granule_type, hdf_ret_rdr

use obs_mod, only : make_obs_sequence, initialize_obs_sequence, &
                             compute_thin_factor

implicit none

! ----------------------------------------------------------------------
! Declare local parameters
! ----------------------------------------------------------------------

integer                 :: thin_factor, filecount
type(hdf_granule_type)  :: granule
type(obs_sequence_type) :: seq

integer :: io, iunit, index

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'convert_1d_hdf.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

! ----------------------------------------------------------------------
! Declare namelist parameters
! ----------------------------------------------------------------------

integer, parameter :: MAXFILES = 512

character(len=256) :: hdf_files(MAXFILES) = ''
character(len=256) :: hdf_file_list       = ''
character(len=256) :: outputfile          = ''

real(r8) :: min_position = -1.0e6_r8,  &   ! lower position bound
            max_position =  1.0e6_r8       ! upper position bound

real(r8) :: min_MMR_threshold = 1.0e-30
real(r8) :: top_pressure_level = 0.0001    ! no obs higher than this
integer  :: position_thin = 0
logical  :: use_NCEP_errs = .false.
integer  :: version = 1    ! HDF file format version

namelist /convert_1d_hdf_nml/ hdf_files, hdf_file_list, &
                               outputfile, &
                               min_position, max_position, &
                               min_MMR_threshold, top_pressure_level, &
                               position_thin, &
                               use_NCEP_errs, version

! ----------------------------------------------------------------------
! start of executable program code
! ----------------------------------------------------------------------

call initialize_utilities('convert_1d_hdf')
call register_module(source, revision, revdate)

! Initialize the obs_sequence module ...
call static_init_obs_sequence()

!----------------------------------------------------------------------
! Read the namelist
!----------------------------------------------------------------------

call find_namelist_in_file('input.nml', 'convert_1d_hdf_nml', iunit)
read(iunit, nml = convert_1d_hdf_nml, iostat = io)
call check_namelist_read(iunit, io, 'convert_1d_hdf_nml')

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=convert_1d_hdf_nml)
if (do_nml_term()) write(    *      , nml=convert_1d_hdf_nml)

! When this routine returns, the hdf_files variable will have
! all the filenames, regardless of which way they were specified.
filecount = set_filename_list(hdf_files, hdf_file_list, "convert_1d_hdf")

! Used to estimate the max size of the output sequence
thin_factor = compute_thin_factor(position_thin, 1)

! Initialize an empty obs_seq to start
seq = initialize_obs_sequence(filecount, thin_factor)

! For each input file
do index=1, filecount

   ! Read from HDF file into a derived type that holds all the information
   call hdf_ret_rdr(hdf_files(index), granule, version)   

   ! Convert derived type information to DART sequence
   call make_obs_sequence(seq, granule, min_position, max_position, &
                          min_MMR_threshold, top_pressure_level, &
                          position_thin, use_NCEP_errs, version)

enddo

! Write the sequence to a disk file
call write_obs_seq(seq, outputfile) 
 
! Release the sequence memory
call destroy_obs_sequence(seq)

call error_handler(E_MSG, source, 'Finished successfully.', source, revision, revdate)
call finalize_utilities()

end program oned_hdf_converter_mod