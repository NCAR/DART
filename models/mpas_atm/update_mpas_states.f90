! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

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
!         Updated in 4 May 2017 for the Manhatten release
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file, nc_check, &
                             get_next_filename
use time_manager_mod, only : time_type, print_time, print_date, operator(-), &
                             get_time, get_date
use direct_netcdf_mod,only : read_transpose, read_variables
use        model_mod, only : static_init_model, statevector_to_analysis_file, &
                             get_model_size, get_model_analysis_filename,     &
                             get_num_vars, get_analysis_time,                 &
                             print_variable_ranges

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character(len=256)  ::  update_input_file_list = 'filter_out.txt'
character(len=256)  :: update_output_file_list = 'filter_in.txt'
logical             :: print_data_ranges        = .true.
integer             :: dom_id                   = 1        ! Not needed now, but MPAS can be run in a regional mode later.

namelist /update_mpas_states_nml/ update_input_file_list, update_output_file_list, dom_id,   &
                                  print_data_ranges
!----------------------------------------------------------------------
character (len=256)   :: next_infile, next_outfile
character (len=256)   :: model_analysis_filename
integer               :: iunit, io, x_size, nvars
integer               :: ncAnlID, ncBckID, istatus
integer               :: filenum
real(r8), allocatable :: statevector(:)
type(time_type)       :: model_time
!----------------------------------------------------------------------

call initialize_utilities(progname='update_mpas_states')

! Read the namelist to get the input filename. 

call find_namelist_in_file("input.nml", "update_mpas_states_nml", iunit)
read(iunit, nml = update_mpas_states_nml, iostat = io)
call check_namelist_read(iunit, io, "update_mpas_states_nml")


!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the model namelists
! to set grid sizes, etc.
! Let us keep model_analysis_filename for now until model_mod.f90
! is updated not to use it any more (e.g. for consistency).
!----------------------------------------------------------------------

call static_init_model()
call get_model_analysis_filename(model_analysis_filename)

x_size = get_model_size()
allocate(statevector(x_size))

nvars = get_num_vars()
write(*,*)
write(*,*) 'update_mpas_states: Updating ',nvars,' variables in ',&
            trim(model_analysis_filename)

!----------------------------------------------------------------------
! Reads lists of input mpas (prior) and filter (analysis) files 
!----------------------------------------------------------------------
filenum = 1
fileloop: do        ! until out of files

  ! get a file name from the list, one at a time.
  next_infile  = get_next_filename( update_input_file_list, filenum)
  next_outfile = get_next_filename(update_output_file_list, filenum)
  if (next_infile == '' .or. next_outfile == '') exit fileloop

  call nc_check(nf90_open(trim(next_infile), NF90_NOWRITE, ncAnlID), &
             'update_mpas_states','open '//trim(next_infile))

  ! Overwrite this mpas file for state vector later
  call nc_check(nf90_open(trim(next_outfile), NF90_NOWRITE, ncBckID), &
             'update_mpas_states','open '//trim(next_outfile))

  !----------------------------------------------------------------------
  ! Read the model time
  !----------------------------------------------------------------------
  model_time = get_analysis_time(ncBckID, trim(next_outfile))

  !----------------------------------------------------------------------
  ! Read analysis state vector (assuming to be available at the model time)
  !----------------------------------------------------------------------
  !call read_transpose(state_ens_handle, name_handle, domain, dart_index, limit_mem)
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
  endif

  !----------------------------------------------------------------------
  ! update the current model state vector
  !----------------------------------------------------------------------
  call statevector_to_analysis_file(statevector, next_outfile, model_time)

  !----------------------------------------------------------------------
  ! Log what we think we're doing, and exit.
  !----------------------------------------------------------------------

  call print_date( model_time,'update_mpas_states:model date')
  call print_time( model_time,'update_mpas_states:model time')
  call print_date( model_time,'update_mpas_states:model date',logfileunit)
  call print_time( model_time,'update_mpas_states:model time',logfileunit)

  filenum = filenum + 1

end do fileloop

call finalize_utilities()

end program update_mpas_states

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
