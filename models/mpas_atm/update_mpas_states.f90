! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program update_mpas_states

!----------------------------------------------------------------------
! purpose: interface between DART and the model model
!
! method: Read DART state vector in netcdf and overwrite values in an mpas file
!         to advance model for each ensemble member after running this program.
!
!         The update_mpas_states_nml namelist setting for dart_analysis_filename
!         (e.g., the name of the output netcdf file from filter) should be matched
!         with restart_out_file_name in &filter_nml (for each member).
!         An input mpas (background) filename comes from model_analysis_filename 
!         in model_nml, thus no need to define it again.
!         
! author: Soyoung Ha 23 Aug 16
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file, nc_check
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

character(len=256)  :: dart_analysis_filename = 'filter_restart.nc'
logical             :: print_data_ranges        = .true.
integer             :: dom_id                   = 1           ! Not needed now, but MPAS can be run in a regional mode later.

namelist /update_mpas_states_nml/ dart_analysis_filename, dom_id,   &
                                  print_data_ranges

!----------------------------------------------------------------------

integer               :: iunit, io, x_size, nvars
integer               :: ncAnlID, ncBckID, istatus
real(r8), allocatable :: statevector(:)
character(len=256)    :: model_analysis_filename
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
! Reads input mpas (prior) and filter (analysis) files 
!----------------------------------------------------------------------

call nc_check(nf90_open(trim(dart_analysis_filename), NF90_NOWRITE, ncAnlID), &
             'update_mpas_states','open '//trim(dart_analysis_filename))

! Overwrite this mpas file for state vector later
call nc_check(nf90_open(trim(model_analysis_filename), NF90_NOWRITE, ncBckID), &
             'update_mpas_states','open '//trim(model_analysis_filename))

!----------------------------------------------------------------------
! Read the model time
!----------------------------------------------------------------------
model_time = get_analysis_time(ncBckID, trim(model_analysis_filename))

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
    call print_variable_ranges(statevector, 'Analysis states')
endif

!----------------------------------------------------------------------
! update the current model state vector
!----------------------------------------------------------------------
call statevector_to_analysis_file(statevector, model_analysis_filename, model_time)

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

call print_date( model_time,'update_mpas_states:model date')
call print_time( model_time,'update_mpas_states:model time')
call print_date( model_time,'update_mpas_states:model date',logfileunit)
call print_time( model_time,'update_mpas_states:model time',logfileunit)

call finalize_utilities()

end program update_mpas_states

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
