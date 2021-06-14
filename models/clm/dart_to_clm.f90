! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program dart_to_clm

!----------------------------------------------------------------------
! purpose: Update the CLM restart file with the DART posterior. If need
!          be, take the posterior SWE (a diagnostic variable) and
!          repartition it into prognostic snow layers.
!
! method: Read DART posterior and replace the valid values in the
!         CLM restart file. The 'clm_to_dart' program replaced the 
!         CLM indeterminate values with the DART MISSING_R8. Since 
!----------------------------------------------------------------------

use        types_mod, only : r8, MISSING_R8

use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             open_file, close_file, &
                             logfileunit, nmlfileunit, &
                             do_nml_file, do_nml_term, &
                             error_handler, E_MSG, E_ERR

use time_manager_mod, only : time_type, print_time, print_date, operator(/=)

use        model_mod, only : static_init_model, read_model_time

use   state_structure_mod, only : get_num_variables, &
                                  get_num_dims, get_dim_name,         &
                                  get_dim_length, get_variable_name,  &
                                  do_io_update, get_variable_size,    &
                                  get_model_variable_indices,         &
                                  get_has_missing_value, &
                                  get_missing_value, &
                                  get_FillValue, &
                                  get_domain_size

use netcdf_utilities_mod, only : nc_open_file_readonly, &
                                 nc_open_file_readwrite,&
                                 nc_get_variable_num_dimensions, &
                                 nc_get_variable_dimension_names, &
                                 nc_get_variable_size,  &
                                 nc_get_variable, &
                                 nc_put_variable, &
                                 nc_close_file, &
                                 NF90_MAX_NAME, NF90_MAX_VAR_DIMS

implicit none

character(len=*), parameter :: source = 'dart_to_clm.f90'

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character(len=256) :: dart_to_clm_input_file = 'dart_posterior.nc'
character(len=256) :: dart_to_clm_output_file = 'clm_restart.nc'
logical            :: repartition_swe = .false.
integer            :: verbose = 0

namelist /dart_to_clm_nml/ dart_to_clm_input_file, &
                           dart_to_clm_output_file, &
                           repartition_swe, verbose

!----------------------------------------------------------------------

integer            :: iunit, io, dom_restart, ivar, rank, irank
integer            :: ncid_dart, ncid_clm
type(time_type)    :: dart_time, clm_time

character(len=512) :: string1, string2, string3

!======================================================================

call initialize_utilities(progname='dart_to_clm')

! Call model_mod:static_init_model() to reads the clm namelist

call static_init_model()

! Read the namelist to get the input/output filenames. 

call find_namelist_in_file("input.nml", "dart_to_clm_nml", iunit)
read(iunit, nml = dart_to_clm_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_clm_nml")

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=dart_to_clm_nml)
if (do_nml_term()) write(     *     , nml=dart_to_clm_nml)

ncid_dart = nc_open_file_readonly( dart_to_clm_input_file,  'opening to read posterior')
ncid_clm  = nc_open_file_readwrite(dart_to_clm_output_file, 'opening to write restart')

write(string1,*)'converting DART file "'//trim(dart_to_clm_input_file)//'" to'
write(string2,*)'clm restart file "'//trim(dart_to_clm_output_file)//'"'
call error_handler(E_MSG,source,string1,text2=string2)

! Make sure we are updating the right restart file - 

dart_time = read_model_time(dart_to_clm_input_file)
clm_time  = read_model_time(dart_to_clm_output_file)

if (dart_time /= clm_time) then
   write(string1,*)'DART timestamp "'//trim(dart_to_clm_input_file)//'" /= '
   write(string2,*)'CLM  timestamp "'//trim(dart_to_clm_output_file)//'"'
   call error_handler(E_MSG,source,string1,text2=string2)
endif

call print_date( clm_time,'dart_to_clm:clm  model date')
call print_time( clm_time,'dart_to_clm:DART model time')
call print_date( clm_time,'dart_to_clm:clm  model date',logfileunit)
call print_time( clm_time,'dart_to_clm:DART model time',logfileunit)

! update the current clm restart file
!
! Get the DART posterior variable names from the restart domain (domain 1).
! The other domains do not need to be updated - they are not used for
! forecasting.

! There are some dependencies on the state structure, which come from the shapefile,
! and then the two files in question must also be compatible ...

dom_restart = 1

UPDATE : do ivar=1, get_num_variables(dom_restart)

   rank = get_num_dims(dom_restart,ivar)

   if (rank == 1) then

      call replace_values_1D(dom_restart, ivar, ncid_dart, ncid_clm)

   elseif (rank == 2) then

      call replace_values_2D(dom_restart, ivar, ncid_dart, ncid_clm)

   else
      write(string1, *) 'no support for data array of dimension ', rank
      call error_handler(E_ERR, source, string1, source)
   endif
enddo UPDATE

call nc_close_file(ncid_clm,  source)
call nc_close_file(ncid_dart, source)

! Log what we think we're doing, and exit.

call finalize_utilities('dart_to_clm')


!===============================================================================
contains
!===============================================================================

!> The DART posterior has had all the DART MISSING_R8 values replaced by
!> the CLM declared '_FillValue' and the clamping values have been honored. 
!> Get the 'original' variable from the netcdf file.

subroutine replace_values_1D(dom_id, ivar, ncid_dart, ncid_clm)

integer, intent(in) :: dom_id
integer, intent(in) :: ivar
integer, intent(in) :: ncid_dart
integer, intent(in) :: ncid_clm

character(len=*), parameter :: routine = 'replace_values_1D'

character(len=NF90_MAX_NAME) :: varname
real(r8), allocatable :: dart_array(:)
real(r8), allocatable :: clm_array(:)
real(r8) :: special
integer  :: varsize(1)

varname = get_variable_name(dom_id,ivar)

! I am not sure I like getting the special value from the state structure as
! opposed to getting it right from the value in the ncid_dart variable ...
! The state structure gets it from the 'shapefile' listed in model_nml ... 

! If it has both, the core DART routines have already made sure they are identical
! The problem here is that if only the missing_value -or- FillValue is used, the
! other may have a bogus value. We only have one logical query routine.  
if (get_has_missing_value(dom_id,ivar)) call get_missing_value(dom_id,ivar,special)
if (get_has_missing_value(dom_id,ivar)) call get_FillValue(    dom_id,ivar,special)

! Make sure variables in both files are identical in shape, etc.
call Compatible_Variables(varname, ncid_dart, ncid_clm, varsize)

allocate(dart_array(varsize(1)), clm_array(varsize(1)))

call nc_get_variable(ncid_clm,  varname,  clm_array)
call nc_get_variable(ncid_dart, varname, dart_array)

where(dart_array /= special) clm_array = dart_array

call nc_put_variable(ncid_clm, varname, clm_array, routine)

deallocate(dart_array, clm_array)

end subroutine replace_values_1D


!------------------------------------------------------------------
!>
!> The DART posterior has had all the DART MISSING_R8 values replaced by
!> the CLM declared '_FillValue' and the clamping values have been honored. 
!> Get the 'original' variable from the netcdf file.

subroutine replace_values_2D(dom_id, ivar, ncid_dart, ncid_clm)

integer, intent(in) :: dom_id
integer, intent(in) :: ivar
integer, intent(in) :: ncid_dart
integer, intent(in) :: ncid_clm

character(len=*), parameter :: routine = 'replace_values_2D'

character(len=NF90_MAX_NAME) :: varname
real(r8), allocatable :: dart_array(:,:)
real(r8), allocatable :: clm_array(:,:)
real(r8) :: special
integer  :: varsize(2)

varname = get_variable_name(dom_id,ivar)

! If it has both, the core DART routines have already made sure they are identical
! The problem here is that if only the missing_value -or- FillValue is used, the
! other may have a bogus value. We only have one logical query routine.  
if (get_has_missing_value(dom_id,ivar)) call get_missing_value(dom_id,ivar,special)
if (get_has_missing_value(dom_id,ivar)) call get_FillValue(    dom_id,ivar,special)

! Make sure variables in both files are identical in shape, etc.
call Compatible_Variables(varname, ncid_dart, ncid_clm, varsize)

allocate(dart_array(varsize(1),varsize(2)), clm_array(varsize(1),varsize(2)))

call nc_get_variable(ncid_clm,  varname,  clm_array)
call nc_get_variable(ncid_dart, varname, dart_array)

where(dart_array /= special) clm_array = dart_array

call nc_put_variable(ncid_clm, varname, clm_array, routine)

deallocate(dart_array, clm_array)

end subroutine replace_values_2D


!------------------------------------------------------------------
!>

subroutine Compatible_Variables(varname, ncid_dart, ncid_clm, varsize)

character(len=*), intent(in)  :: varname
integer,          intent(in)  :: ncid_dart
integer,          intent(in)  :: ncid_clm
integer,          intent(out) :: varsize(:)

character(len=*),parameter :: routine = 'Compatible_Variables'

integer ::  clm_numdims,  clm_dimlengths(NF90_MAX_VAR_DIMS)
integer :: dart_numdims, dart_dimlengths(NF90_MAX_VAR_DIMS)

character(len=NF90_MAX_NAME) :: clm_dimnames( NF90_MAX_VAR_DIMS)
character(len=NF90_MAX_NAME) :: dart_dimnames(NF90_MAX_VAR_DIMS)

call nc_get_variable_num_dimensions( ncid_clm, varname, clm_numdims,    routine)
call nc_get_variable_dimension_names(ncid_clm, varname, clm_dimnames,   routine)
call nc_get_variable_size(           ncid_clm, varname, clm_dimlengths, routine)

call nc_get_variable_num_dimensions( ncid_dart, varname, dart_numdims,    routine)
call nc_get_variable_dimension_names(ncid_dart, varname, dart_dimnames,   routine)
call nc_get_variable_size(           ncid_dart, varname, dart_dimlengths, routine)

if (clm_numdims /= dart_numdims) then
   write(string1,*)'Variable ranks not identical and must be.'
   write(string2,*)'CLM  "'//trim(varname)//'" has rank ',clm_numdims
   write(string3,*)'DART "'//trim(varname)//'" has rank ',dart_numdims
   call error_handler(E_ERR, routine, string1, source, text2=string2, text3=string3)
endif

! Ensure CLM netCDF variable is conformable with DART variable.

DimCheck : do irank = 1,clm_numdims

   if ( clm_dimlengths(irank) /= dart_dimlengths(irank) ) then
      write(string1,*)'Variable dimension lengths not identical and must be.'
      write(string2,*)'CLM  "'//trim(varname)//'" dim=',irank, ', length=',clm_dimlengths(irank)
      write(string2,*)'DART "'//trim(varname)//'" dim=',irank, ', length=',dart_dimlengths(irank)
      call error_handler(E_ERR, routine, string1, source, text2=string2, text3=string3)
   endif

   ! Do we want to check dimension name too?
   if ( clm_dimnames(irank) /= dart_dimnames(irank) ) then
      write(string1,*)'Variable dimension names not identical and must be.'
      write(string2,*)'CLM  "'//trim(varname)//'" dim=',irank, ', name=',clm_dimnames(irank)
      write(string2,*)'DART "'//trim(varname)//'" dim=',irank, ', name=',dart_dimnames(irank)
      call error_handler(E_ERR, routine, string1, source, text2=string2, text3=string3)
   endif

   varsize(irank) = clm_dimlengths(irank)

enddo DimCheck

if (verbose > 0) then
      write(string1,*)'Replacing indeterminate values ...'
      write(string2,*)'CLM  "'//trim(varname)//'" shape = ', clm_dimlengths(1:clm_numdims)
      write(string3,*)'DART "'//trim(varname)//'" shape = ',dart_dimlengths(1:clm_numdims)
      call error_handler(E_MSG, routine, string1, source, text2=string2, text3=string3)
endif

end subroutine Compatible_Variables



!===============================================================================

end program dart_to_clm
