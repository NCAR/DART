! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program test_state_structure

use           types_mod, only : r8, missing_r8
use        obs_kind_mod, only : get_raw_obs_kind_index
use       utilities_mod, only : register_module, error_handler, E_ERR, E_MSG,       &
                                open_file, close_file, nc_check, get_next_filename, &
                                find_namelist_in_file, check_namelist_read,         &
                                do_nml_file, do_nml_term, nmlfileunit,              &
                                initialize_utilities, finalize_utilities, do_output,&
                                to_upper
use   mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities,   &
                                task_sync, my_task_id
use state_structure_mod, only : static_init_state_type,     &
                                add_domain,                 &
                                get_domain_size,            &
                                get_num_domains,            &
                                get_variable_size,          &
                                get_variable_name,          &
                                get_kind_string,            &
                                get_kind_index,             &
                                get_varid_from_kind,        & 
                                get_varids_from_kind,       & 
                                get_num_variables,          &
                                get_num_dims,               &
                                get_dim_lengths,            &
                                get_dim_length,             &
                                get_dim_name,               &
                                get_io_num_dims,            &
                                get_io_dim_ids,             &
                                get_io_dim_lengths,         &
                                get_io_num_unique_dims,     &
                                get_io_unique_dim_name,     &
                                get_io_unique_dim_length,   &
                                add_time_unlimited,         &
                                get_unlimited_dimid,        &
                                set_var_id,                 &
                                get_io_clamping_maxval,     &
                                get_io_clamping_minval,     &
                                do_io_clamping,             &
                                do_io_update,               &
                                get_index_start,            &
                                get_index_end,              &
                                get_sum_variables,          &
                                get_sum_variables_below,    &
                                get_model_variable_indices, &
                                get_dart_vector_index,      &
                                get_num_varids_from_kind,   &
                                add_dimension_to_variable,  &
                                finished_adding_domain,     &
                                state_structure_info
use netcdf
implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=256) :: string1, string2, string3

logical, save :: module_initialized = .false.

integer :: iunit, io

integer :: nfields = -1
integer, parameter :: num_state_table_columns = 5
integer, parameter :: max_state_variables     = 10
character(len=NF90_MAX_NAME) :: variable_table( max_state_variables, num_state_table_columns )
integer  :: state_kind_list  ( max_state_variables )
logical  :: state_update_list( max_state_variables )
real(r8) :: state_clamp_list ( max_state_variables, 2 )

! Codes for interpreting the columns of the variable_table
integer, parameter :: VT_VARNAMEINDX  = 1 ! ... variable name
integer, parameter :: VT_KINDINDX     = 2 ! ... DART kind
integer, parameter :: VT_MINVALINDX   = 3 ! ... minimum value if any
integer, parameter :: VT_MAXVALINDX   = 4 ! ... maximum value if any
integer, parameter :: VT_STATEINDX    = 5 ! ... update (state) or not

! namelist variables
logical                      :: debug = .false.
character(len=256)           :: one_dom_restart_filename
character(len=NF90_MAX_NAME) :: one_dom_variables(max_state_variables*num_state_table_columns) = ' '

! namelist items we are going to create/overwrite
namelist /test_state_structure_nml/  &
   debug,                            &
   one_dom_restart_filename,         &
   one_dom_variables
                                
integer :: dom_id


! main code here
 
! initialize the dart libs
call initialize_module()

! Read back the namelist entry
call find_namelist_in_file("input.nml", "test_state_structure_nml", iunit)
read(iunit, nml = test_state_structure_nml, iostat = io)
call check_namelist_read(iunit, io, "test_state_structure_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=test_state_structure_nml)
if (do_nml_term()) write(     *     , nml=test_state_structure_nml)

call parse_variable_table( one_dom_variables, nfields, variable_table, &
                             kind_list=state_kind_list,                   &
                            clamp_list=state_clamp_list,                  &
                           update_list=state_update_list)

! add single domain from file. 
dom_id =  add_domain( one_dom_restart_filename, nfields,                          &
                      var_names   = variable_table   (1:nfields, VT_VARNAMEINDX), &
                      kind_list   = state_kind_list  (1:nfields),                 &
                      clamp_vals  = state_clamp_list (1:nfields,:),               &
                      update_list = state_update_list(1:nfields) )


! print out state structure information
call state_structure_info(dom_id)

call unit_test(one_dom_restart_filename, dom_id)

! finalize test_state_structure
call error_handler(E_MSG,'test_state_structure','Finished successfully.',source,revision,revdate)
call finalize_mpi_utilities()

! end of main code


contains

!----------------------------------------------------------------------

subroutine initialize_module()

  call initialize_mpi_utilities('test_state_structure')
  call register_module(source, revision, revdate)
  module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------
!>  This routine checks the user input against the variables available in the
!>  input netcdf file to see if it is possible to construct the DART state vector
!>  specified by the input.nml:model_nml:clm_variables  variable.
!>  Each variable must have 6 entries.
!>  1: variable name
!>  2: DART KIND
!>  3: minimum value - as a character string - if none, use 'NA'
!>  4: maximum value - as a character string - if none, use 'NA'
!>  5: does the variable get updated in the restart file or not ...
!>     only variables from restart files may be updated.
!>     'UPDATE'       => update the variable in the restart file
!>     'NO_COPY_BACK' => do not copy the variable back to the restart file
!>     all these variables will be updated INTERNALLY IN DART
!>     only variables marked '.r', 'UPDATE' will be modified for CLM.
!>
!>  The calling code should check to see if the variable exists.

subroutine parse_variable_table( state_variables, ngood, table, kind_list, clamp_list, update_list )

character(len=*), dimension(:),   intent(in)  :: state_variables
integer,                          intent(out) :: ngood
character(len=*), dimension(:,:), intent(out) :: table
integer,          dimension(:),   optional, intent(out) :: kind_list
real(r8),         dimension(:,:), optional, intent(out) :: clamp_list
logical,          dimension(:),   optional, intent(out) :: update_list

integer :: nrows, ncols, i
character(len=NF90_MAX_NAME) :: varname       ! column 1
character(len=NF90_MAX_NAME) :: dartstr       ! column 2
character(len=NF90_MAX_NAME) :: minvalstring  ! column 3
character(len=NF90_MAX_NAME) :: maxvalstring  ! column 4
character(len=NF90_MAX_NAME) :: update        ! column 5

real(r8) :: minvalue, maxvalue
integer  :: ios

nrows = size(table,1)
ncols = size(table,2)

! This loop just repackages the 1D array of values into a 2D array.
! We can do some miniminal checking along the way.
! Determining which file to check is going to be more complicated.

ngood = 0
MyLoop : do i = 1, nrows

   varname      = trim(state_variables(ncols*i - 4))
   dartstr      = trim(state_variables(ncols*i - 3))
   minvalstring = trim(state_variables(ncols*i - 2))
   maxvalstring = trim(state_variables(ncols*i - 1))
   update = trim(state_variables(ncols*i    ))

   call to_upper(update)

   table(i,VT_VARNAMEINDX) = trim(varname)
   table(i,VT_KINDINDX)    = trim(dartstr)
   table(i,VT_MINVALINDX)  = trim(minvalstring)
   table(i,VT_MAXVALINDX)  = trim(maxvalstring)
   table(i,VT_STATEINDX)   = trim(update)

   ! If the first element is empty, we have found the end of the list.
   if ( table(i,1) == ' ' ) exit MyLoop

   ! Any other condition is an error.
   if ( any(table(i,:) == ' ') ) then
      string1 = 'input.nml &model_nml:clm_variables not fully specified'
      string2 = 'must be 6 entries per variable. Last known variable name is'
      string3 = '['//trim(table(i,1))//'] ... (without the [], naturally)'
      call error_handler(E_ERR, 'parse_variable_table', string1, &
         source, revision, revdate, text2=string2, text3=string3)
   endif

   ! Make sure DART kind is valid

   if ( present(kind_list) )then
      kind_list(i) = get_raw_obs_kind_index(dartstr)
      if( kind_list(i)  < 0 ) then
         write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
         call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
      endif
   endif

   ! Set clamping ranges

   if ( present(clamp_list) )then
      clamp_list(:,:) = missing_r8
      
      ! If the character string can be interpreted as an r8, great.
      ! If not, there is no value to be used.
      
      read(table(i,VT_MINVALINDX),*,iostat=ios) minvalue
      if (ios == 0) clamp_list(i,1) = minvalue
      
      read(table(i,VT_MINVALINDX),*,iostat=ios) maxvalue
      if (ios == 0) clamp_list(i,2) = maxvalue
   endif

   ! Make sure the update variable has a valid string

   if ( present(update_list) )then
      select case (update)
         case ('UPDATE')
            update_list(i) = .true.
         case ('NO_COPY_BACK')
            update_list(i) = .false.
         case default
            write(string1,'(A)')  'only UPDATE or NO_COPY_BACK supported ', &
                                  'in model_state_variable namelist'
            write(string2,'(6A)') 'you provided : ', &
                                  trim(varname), ', ', trim(dartstr), ', ', trim(update)
            call error_handler(E_ERR,'verify_state_variables', &
                               string1,source,revision,revdate, text2=string2)
      end select
   endif

   ! Record the contents of the DART state vector

   if ( debug .and. do_output() ) then
      write(*,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2)),' ', &
                                     trim(table(i,3)), ' ', trim(table(i,4)),' ', &
                                     trim(table(i,5))
   endif

   ngood = ngood + 1
enddo MyLoop

end subroutine parse_variable_table

!----------------------------------------------------------------------
!>  This routine checks does unit testing
subroutine unit_test(filename, domid)
character(len=*), intent(in) :: filename
integer,          intent(in) :: domid

! netcdf variables
integer :: ncid
character(nf90_max_name) :: varname

! routine variables
integer :: i, varsize1, varsize2, model_size1, model_size2

if ( do_output() ) write(*,*) "Unit Testing"

call nc_check( nf90_open(trim(filename), NF90_NOWRITE, ncid), &
               'get_netcdf_variable_size', 'open '//trim(filename))

! testing get_variable_size()
model_size1 = 0
do i = 1, nfields
   varname = variable_table(i,VT_VARNAMEINDX)
   varsize1 = get_netcdf_variable_size(ncid, filename, varname)
   varsize2 = get_variable_size(domid, i)
   if ( varsize1 /= varsize2 ) then
      write(string1,*) 'variable ', trim(varname), ' : FAILED'
      write(string2,*) 'varsize from netcdf              = ', varsize1
      write(string3,*) 'varsize from state_structure_mod = ', varsize2
      call error_handler(E_MSG, 'unit_test', string1, &
         source, revision, revdate, text2=string2, text3=string3)
   else
      write(string1,*) 'variable ', trim(varname), ' : PASSED'
      call error_handler(E_MSG, 'unit_test', string1, &
         source, revision, revdate)
   endif    
   model_size1 = model_size1 + varsize1
enddo

! testing get_domain_size()
model_size2 = get_domain_size(domid)
if ( model_size1 /= model_size2 ) then
   write(string1,*) 'model_size : FAILED'
   write(string2,*) 'model_size from netcdf              = ', model_size1
   write(string3,*) 'model_size from state_structure_mod = ', model_size2
   call error_handler(E_MSG, 'unit_test', string1, &
      source, revision, revdate, text2=string2, text3=string3)
else
   write(string1,*) 'model_size : PASSED'
   call error_handler(E_MSG, 'unit_test', string1, &
      source, revision, revdate)
endif    

call nc_check( nf90_close(ncid), 'unit_test', 'close '//trim(filename) )

end subroutine unit_test

function get_netcdf_variable_size(ncid, filename, variable) result(var_size)
integer,          intent(in) :: ncid
character(len=*), intent(in) :: filename
character(len=*), intent(in) :: variable
integer :: var_size

integer :: numDims, dimLength, varId
integer, dimension(nf90_max_var_dims) :: varDimIds
character(len=nf90_max_name) :: dimName
integer :: i

call nc_check( nf90_inq_varid(ncid, trim(variable), varId), &
               'get_netcdf_variable_size', 'nf90_inq_varid'//trim(variable))

call nc_check( nf90_inquire_variable(ncid, varId, ndims=numDims, dimids=varDimIds), &
               'get_netcdf_variable_size', 'nf90_inquire_variable'//trim(variable))

var_size = 1
do i = 1, numDims
   call nc_check( nf90_inquire_dimension(ncid, varDimIds(i), len=dimLength, name=dimName), &
                  'get_netcdf_variable_size', 'nf90_inquire_dimension'//trim(dimName))
   var_size = var_size * dimLength      
enddo

end function get_netcdf_variable_size

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
