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
!  variables that are not wind, copied from one file to the another
!  variables that are wind, reconstructed
!  
! author: Soyoung Ha 23 Aug 16
!         Updated in 4 May 2017 for the Manhattan release
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, get_next_filename, E_ERR, E_MSG, error_handler
use time_manager_mod, only : time_type, print_time, print_date, operator(-), &
                             get_time, get_date, operator(/=)
use        model_mod, only : static_init_model, &
                             get_model_size, &
                             get_analysis_time

use state_structure_mod, only : get_num_variables, get_variable_name, &
                                get_variable_size

use netcdf_utilities_mod, only : nc_open_file_readonly, &
                                 nc_open_file_readwrite, &
                                 nc_get_variable, nc_put_variable, &
                                 nc_close_file
                                 
implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'models/mpas_atm/update_mpas_states.f90'

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------
character(len=256)  ::  update_input_file_list = 'filter_out.txt'
character(len=256)  :: update_output_file_list = 'filter_in.txt'

namelist /update_mpas_states_nml/ update_input_file_list, update_output_file_list

!----------------------------------------------------------------------
character (len=256)   :: next_infile, next_outfile
character (len=512)   :: string1
integer               :: iunit, io
integer               :: ncAnlID, ncBckID, istatus
integer               :: filenum, i
real(r8), allocatable :: variable(:)
type(time_type)       :: model_time
type(time_type)       :: state_time
!----------------------------------------------------------------------

call initialize_utilities(progname=source)

! Read the namelist to get the input filename. 
call find_namelist_in_file("input.nml", "update_mpas_states_nml", iunit)
read(iunit, nml = update_mpas_states_nml, iostat = io)
call check_namelist_read(iunit, io, "update_mpas_states_nml")

call static_init_model()

!----------------------------------------------------------------------
! Reads lists of input mpas (prior) and filter (analysis) files 
! HK @todo this is bad to have a serial loop around files.
!----------------------------------------------------------------------
filenum = 1
fileloop: do        ! until out of files

  ! get a file name from the list, one at a time.
  next_infile  = get_next_filename( update_input_file_list, filenum)
  next_outfile = get_next_filename(update_output_file_list, filenum)
  if (next_infile == '' .or. next_outfile == '') exit fileloop

  ncAnlID = nc_open_file_readonly(next_infile, 'update_mpas_states - open readonly')
  ncBckID = nc_open_file_readwrite(next_outfile, 'update_mpas_states - open readwrite')

  model_time = get_analysis_time(ncBckID, next_outfile)
  state_time = get_analysis_time(ncAnlID, next_infile)
  call print_time(state_time,'DART current time')
  call print_time(model_time,'mpas current time')

  if ( model_time /= state_time ) then
     call print_time(state_time,'DART current time',logfileunit)
     call print_time(model_time,'mpas current time',logfileunit)
     write(string1,*) trim(next_infile),' current time must equal model time'
     call error_handler(E_ERR,'update_mpas_states',string1,source)
  endif

  ! copy variables that are not wind from analysis to background
  varloop: do i = 1, get_num_variables(1)

      allocate(variable(get_variable_size(1, i)))

      call nc_get_variable(ncAnlID, get_variable_name(1,i), variable)
      call nc_put_variable(ncBckID, get_variable_name(1,i), variable)

      deallocate(variable)

  enddo varloop

  call error_handler(E_MSG, 'Overwriting states in ',trim(next_outfile), source)

  call print_date( model_time,'update_mpas_states:model date')
  call print_time( model_time,'update_mpas_states:model time')
  call print_date( model_time,'update_mpas_states:model date',logfileunit)
  call print_time( model_time,'update_mpas_states:model time',logfileunit)

  call nc_close_file(ncAnlID,'update_mpas_states')
  call nc_close_file(ncBckID,'update_mpas_states')

  filenum = filenum + 1

end do fileloop

call finalize_utilities()

end program update_mpas_states

