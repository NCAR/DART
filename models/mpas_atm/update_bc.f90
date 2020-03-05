! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: update_bc.f90 12937 2018-11-26 23:01:09Z nancy@ucar.edu $

program update_bc

!----------------------------------------------------------------------
! purpose: interface between DART and the model model
!
! method: Read DART analysis vector in netcdf and replace the corresponding
!         field in the mpas file to advance model after running this program.
!         Updated to process all ensemble members.
!
!         The update_bc_nml namelist defines the input and output file
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
                             get_next_filename, E_ERR, error_handler
use time_manager_mod, only : time_type, print_time, print_date, operator(-), &
                             get_time, get_date, operator(/=)
use direct_netcdf_mod,only : read_transpose, read_variables
use        model_mod, only : static_init_model, statevector_to_analysis_file, &
                             get_model_size, get_init_template_filename,      &
                             get_analysis_time, statevector_to_boundary_file, &
                             print_variable_ranges, &  ! , get_num_vars
                             set_lbc_variables, force_u_into_state
use state_structure_mod, only : get_num_variables, get_domain_size
use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://svn-dares-dart.cgd.ucar.edu/DART/releases/Manhattan/models/mpas_atm/update_bc.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 12937 $"
character(len=128), parameter :: revdate  = "$Date: 2018-11-26 16:01:09 -0700 (Mon, 26 Nov 2018) $"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character(len=256)  :: update_analysis_file_list = 'filter_in.txt'
character(len=256)  :: update_boundary_file_list = 'boundary_inout.txt'
integer             :: debug = 0


namelist /update_bc_nml/ update_analysis_file_list, update_boundary_file_list, debug

!----------------------------------------------------------------------
character (len=256)   :: next_infile, next_outfile
character (len=256)   :: bdy_template_filename
character (len=256)   :: static_filename
character (len=256)   :: string1
integer               :: iunit, io, x_size, nanlvars, nbdyvars
integer               :: d1size, d2size
integer               :: ncAnlID, ncBdyID, istatus
integer               :: filenum
real(r8), allocatable :: statevector(:)
type(time_type)       :: model_time
type(time_type)       :: state_time
!----------------------------------------------------------------------

call initialize_utilities(progname='update_bc')

! Read the namelist to get the input filename. 

call find_namelist_in_file("input.nml", "update_bc_nml", iunit)
read(iunit, nml = update_bc_nml, iostat = io)
call check_namelist_read(iunit, io, "update_bc_nml")

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the model namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------

! get the first member file to use as a template
bdy_template_filename = get_next_filename(update_boundary_file_list, 1)

call force_u_into_state()
call set_lbc_variables(bdy_template_filename)

call static_init_model()
call get_init_template_filename(static_filename)

x_size = get_model_size()
allocate(statevector(x_size))

! use get_num_variables() after setting the domains
! separate number for analysis file, boundary file
nanlvars = get_num_variables(1)
nbdyvars = get_num_variables(2)

write(*,*)
write(*,*) 'update_bc: Updating ',nbdyvars,' variables'

!----------------------------------------------------------------------
! Reads lists of input mpas (prior) and filter (analysis) files 
!----------------------------------------------------------------------
filenum = 1
fileloop: do        ! until out of files

  ! get a file name from the list, one at a time.
  next_infile  = get_next_filename(update_analysis_file_list, filenum)
  next_outfile = get_next_filename(update_boundary_file_list, filenum)
  if (next_infile == '' .or. next_outfile == '') exit fileloop

  !----------------------------------------------------------------------
  ! Reads input lbc (prior) and filter (analysis) files 
  !----------------------------------------------------------------------
  call nc_check(nf90_open(trim(next_infile), NF90_NOWRITE, ncAnlID), &
             'update_bc','open '//trim(next_infile))

  ! Overwrite this mpas file for state vector later
  call nc_check(nf90_open(trim(next_outfile), NF90_WRITE, ncBdyID), &
             'update_bc','open '//trim(next_outfile))

  !----------------------------------------------------------------------
  ! Read the model time
  !----------------------------------------------------------------------
  model_time = get_analysis_time(ncBdyID, trim(next_outfile))
  state_time = get_analysis_time(ncAnlID, trim(next_infile))
  call print_time(state_time,'DART current time')
  call print_time(model_time,'mpas current time')

  if ( model_time /= state_time ) then
   call print_time(state_time,'DART current time',logfileunit)
   call print_time(model_time,'mpas current time',logfileunit)
   write(string1,*) trim(next_infile),' current time must equal model time'
   call error_handler(E_ERR,'update_bc',string1,source,revision,revdate)
  endif

  !----------------------------------------------------------------------
  ! Read analysis state vector (assuming to be available at the model time)
  !----------------------------------------------------------------------
  !call read_transpose(state_ens_handle, name_handle, domain, dart_index, limit_mem)
  d1size = get_domain_size(1)
  d2size = get_domain_size(2)
  call read_variables(ncAnlID, statevector(1:d1size), 1, nanlvars, domain=1)
  call read_variables(ncBdyID, statevector(d1size+1:d1size+d2size), 1, nbdyvars, domain=2)

  !----------------------------------------------------------------------
  ! update the current model state vector
  !----------------------------------------------------------------------
  write(*,*) 'Updating boundary variables in ',trim(next_outfile)
  call statevector_to_boundary_file(statevector, ncAnlID, ncBdyID)

  !----------------------------------------------------------------------
  ! Log what we think we're doing, and exit.
  !----------------------------------------------------------------------

  call print_date( model_time,'update_bc:model date')
  call print_time( model_time,'update_bc:model time')
  call print_date( model_time,'update_bc:model date',logfileunit)
  call print_time( model_time,'update_bc:model time',logfileunit)

  call nc_check(nf90_close(ncAnlID), &
             'update_mpas_states','close '//trim(next_infile))
  call nc_check(nf90_close(ncBdyID), &
             'update_mpas_states','close '//trim(next_outfile))
  filenum = filenum + 1

end do fileloop

call finalize_utilities()

end program update_bc

! <next few lines under version control, do not edit>
! $URL: https://svn-dares-dart.cgd.ucar.edu/DART/releases/Manhattan/models/mpas_atm/update_bc.f90 $
! $Id: update_bc.f90 12937 2018-11-26 23:01:09Z nancy@ucar.edu $
! $Revision: 12937 $
! $Date: 2018-11-26 16:01:09 -0700 (Mon, 26 Nov 2018) $
