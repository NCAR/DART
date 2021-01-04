! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program update_mpas_states

!----------------------------------------------------------------------
! purpose: interface between DART and the model model
!
! method: Read DART analysis vector in netcdf and replace the corresponding
!         field in the mpas file to advance model after running this program.
!         Updated to process all ensemble members.
!
!         The update_mpas_states_nml namelist defines the input and output file
!         name lists for all ensemble members.
!         The input list should be matched with output_state_file_list in &filter_nml.
!         
! author: Soyoung Ha 23 Aug 16
!         Updated in 4 May 2017 for the Manhattan release
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, get_next_filename, E_ERR, error_handler
use time_manager_mod, only : time_type, print_time, print_date, operator(-), &
                             get_time, get_date, operator(/=)
use direct_netcdf_mod,only : read_transpose, read_variables
use        model_mod, only : static_init_model, statevector_to_analysis_file, &
                             get_model_size, &
                             get_num_vars, get_analysis_time, &
                             print_variable_ranges

use netcdf_utilities_mod, only : nc_open_file_readonly, &
                                 nc_open_file_readwrite, &
                                 nc_close_file

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'models/mpas_atm/update_mpas_states.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character(len=256)  ::  update_input_file_list = 'filter_out.txt'
character(len=256)  :: update_output_file_list = 'filter_in.txt'
logical             :: print_data_ranges        = .true.
integer             :: dom_id                   = 1

namelist /update_mpas_states_nml/ update_input_file_list, update_output_file_list, dom_id,   &
                                  print_data_ranges
!----------------------------------------------------------------------
character (len=256)   :: next_infile, next_outfile
character (len=512)   :: string1
integer               :: iunit, io, x_size, nvars
integer               :: ncAnlID, ncBckID, istatus
integer               :: filenum
real(r8), allocatable :: statevector(:)
type(time_type)       :: model_time
type(time_type)       :: state_time
!----------------------------------------------------------------------

call initialize_utilities(progname=source)

! Read the namelist to get the input filename. 

call find_namelist_in_file("input.nml", "update_mpas_states_nml", iunit)
read(iunit, nml = update_mpas_states_nml, iostat = io)
call check_namelist_read(iunit, io, "update_mpas_states_nml")

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the model namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------

call static_init_model()

nvars  = get_num_vars()
x_size = get_model_size()
allocate(statevector(x_size))

write(*,*)
write(*,*) 'update_mpas_states: Updating ',nvars,' variables'
write(*,*)

!----------------------------------------------------------------------
! Reads lists of input mpas (prior) and filter (analysis) files 
!----------------------------------------------------------------------
filenum = 1
fileloop: do        ! until out of files

  ! get a file name from the list, one at a time.
  next_infile  = get_next_filename( update_input_file_list, filenum)
  next_outfile = get_next_filename(update_output_file_list, filenum)
  if (next_infile == '' .or. next_outfile == '') exit fileloop

  ncAnlID = nc_open_file_readonly(next_infile, 'update_mpas_states - open readonly')
  ncBckID = nc_open_file_readwrite(next_outfile, 'update_mpas_states - open readwrite')

  !----------------------------------------------------------------------
  ! Read the model time
  !----------------------------------------------------------------------
  model_time = get_analysis_time(ncBckID, next_outfile)
  state_time = get_analysis_time(ncAnlID, next_infile)
  call print_time(state_time,'DART current time')
  call print_time(model_time,'mpas current time')

  if ( model_time /= state_time ) then
     call print_time(state_time,'DART current time',logfileunit)
     call print_time(model_time,'mpas current time',logfileunit)
     write(string1,*) trim(next_infile),' current time must equal model time'
     call error_handler(E_ERR,'update_mpas_states',string1,source,revision,revdate)
  endif

  !----------------------------------------------------------------------
  ! Read analysis state vector (assuming to be available at the model time)
  !----------------------------------------------------------------------
  call read_variables(ncAnlID, statevector, 1, nvars, dom_id)

  !----------------------------------------------------------------------
  ! if requested, print out the data ranges variable by variable
  ! (note if we are clamping data values, that happens in the conversion
  ! routine below and these values are before the clamping happens.)
  !----------------------------------------------------------------------
  if (print_data_ranges) then
      write(*,*) 
      write(*,*) ' Input: ', trim(next_infile) 
      write(*,*) 'Output: ', trim(next_outfile)
      call print_variable_ranges(statevector, 'Analysis states')
      write(*,*) 
  endif

  !----------------------------------------------------------------------
  ! update the current model state vector
  !----------------------------------------------------------------------
  write(*,*) 'Overwriting states in ',trim(next_outfile)
  call statevector_to_analysis_file(statevector, ncBckID, next_outfile)

  !----------------------------------------------------------------------
  ! Log what we think we're doing, and exit.
  !----------------------------------------------------------------------

  call print_date( model_time,'update_mpas_states:model date')
  call print_time( model_time,'update_mpas_states:model time')
  call print_date( model_time,'update_mpas_states:model date',logfileunit)
  call print_time( model_time,'update_mpas_states:model time',logfileunit)

  ! Because the files were open with the nc_open_file...() routines,
  ! it is not necessary to supply the filenames for error msg purposes.

  call nc_close_file(ncAnlID,'update_mpas_states')
  call nc_close_file(ncBckID,'update_mpas_states')
  
  filenum = filenum + 1

end do fileloop

call finalize_utilities()

end program update_mpas_states

