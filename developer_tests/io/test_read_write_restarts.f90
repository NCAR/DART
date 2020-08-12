! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program test_read_write_restarts

use            types_mod, only : r8, i8, vtablenamelength
use        utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                                 find_namelist_in_file, check_namelist_read,   &
                                 do_nml_file, do_nml_term, nmlfileunit, to_upper
use    mpi_utilities_mod, only : initialize_mpi_utilities,                     &
                                 finalize_mpi_utilities
use         obs_kind_mod, only : get_index_for_quantity
use  state_vector_io_mod, only : read_state, write_state
use  state_structure_mod, only : add_domain, get_sum_variables
use ensemble_manager_mod, only : init_ensemble_manager, ensemble_type
use     io_filenames_mod, only : io_filenames_init, file_info_type, READ_COPY, &
                                 set_file_metadata, set_io_copy_flag, WRITE_COPY
use     time_manager_mod, only : time_type
use           filter_mod, only : filter_set_initial_time

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer :: num_ens     = 1
integer :: num_domains = 1
integer :: imem, domid

type(file_info_type) :: file_input_handle
type(file_info_type) :: file_output_handle

character(len=256), allocatable :: file_array_input(:,:)
character(len=256), allocatable :: file_array_output(:,:)
character(len=512) :: my_base, my_desc, string1

integer :: iunit, io

type(ensemble_type)  :: ens_handle
type(time_type)      :: time1

logical :: single_file_in = .false.
logical :: read_time_from_file

integer(i8) :: model_size
integer     :: ens_size = 3

! Number of fields in the state vector
integer :: nfields

! DART state vector contents are specified in the input.nml:&model_nml namelist.
integer, parameter :: max_state_variables = 10 
integer, parameter :: num_state_table_columns = 3
character(len=vtablenamelength) :: variable_table( max_state_variables, num_state_table_columns )
integer :: state_kinds_list( max_state_variables )
character(len=vtablenamelength) :: model_variables(max_state_variables * num_state_table_columns ) = ' '

! identifiers for variable_table
integer, parameter :: VAR_NAME_INDEX = 1
integer, parameter :: VAR_QTY_INDEX = 2
integer, parameter :: VAR_UPDATE_INDEX = 3

! namelist variables
character(len=256) :: input_file(1)  = "cf_test.nc"
character(len=256) :: output_file(1) = "cf_test_out.nc"
logical            :: verbose = .false.

! namelist items we are going to create/overwrite
namelist /test_read_write_restarts_nml/ input_file, model_variables, verbose 

! main code here
 
! initialize the dart libs
call initialize_module()

! Read back the namelist entry
call find_namelist_in_file("input.nml", "test_read_write_restarts_nml", iunit)
read(iunit, nml = test_read_write_restarts_nml, iostat = io)
call check_namelist_read(iunit, io, "test_read_write_restarts_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=test_read_write_restarts_nml)
if (do_nml_term()) write(     *     , nml=test_read_write_restarts_nml)

! verify that the model_variables namelist was filled in correctly.  
! returns variable_table which has variable names, and kind strings.
call verify_state_variables(model_variables, nfields, variable_table, state_kinds_list)

domid = add_domain(input_file(1), &
                   num_vars  = nfields, &
                   var_names = variable_table(1:nfields, VAR_NAME_INDEX), &
                   kind_list = state_kinds_list)

! since we are calling add_domain directly instead of through
! static_assim_model_mod we need to get the total number of
! variables from the state_strucutre_mod instead of using
! get_model_size()
model_size = get_sum_variables(1, nfields, domid)

call init_ensemble_manager(ens_handle, ens_size, model_size)

call filter_set_initial_time(-1,-1,time1,read_time_from_file)

! set up the filename handle for reading
file_input_handle  = initialize_filenames(input_file,  READ_COPY)

call read_state(ens_handle, file_input_handle, read_time_from_file, time1)

! set up the filename handle for writing
file_output_handle = initialize_filenames(output_file, WRITE_COPY)

call write_state(ens_handle, file_output_handle)

! finalize test_read_write_restarts
call error_handler(E_MSG,'test_read_write_restarts','Finished successfully.',&
                   source,revision,revdate)

call finalize_mpi_utilities()

! end of main code


contains

!----------------------------------------------------------------------

subroutine initialize_module

  call initialize_mpi_utilities('test_read_write_restarts')
  call register_module(source, revision, revdate)

end subroutine initialize_module

!------------------------------------------------------------------
!> Verify that the namelist was filled in correctly, and check
!> that there are valid entries for the dart_kind. 
!> Returns a table with columns:  
!>
!>    netcdf_variable_name ; dart_kind_string

subroutine verify_state_variables( state_variables, ngood, table, kind_list)

character(len=*),  intent(inout) :: state_variables(:)
integer,           intent(out) :: ngood
character(len=*),  intent(out) :: table(:,:)
integer,           intent(out) :: kind_list(:) ! kind number

integer :: nrows, i
character(len=256) :: varname
character(len=32)  :: kindstr

nrows = size(table,1)

ngood = 0

if ( state_variables(1) == ' ' ) then ! no model_variables namelist provided
   string1 = 'model_nml:model_variables not specified'
   call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
endif

MyLoop : do i = 1, nrows

   varname = trim(state_variables(2*i -1))
   kindstr = trim(state_variables(2*i   ))
   
   table(i,1) = trim(varname)
   table(i,2) = trim(kindstr)

   if ( table(i,1) == ' ' .and. table(i,2) == ' ') exit MyLoop ! Found end of list.

   if ( table(i,1) == ' ' .and. table(i,2) == ' ') then
      string1 = 'model_nml:model_variables not fully specified'
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Make sure DART kind is valid

   kind_list(i) = get_index_for_quantity(kindstr)
   if( kind_list(i)  < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(kindstr)
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif
   
   ! Record the contents of the DART state vector

   if (verbose) then
      write(string1,'(A,I2,6A)') 'variable ',i,' is ',trim(varname), ', ', trim(kindstr)
      call error_handler(E_MSG,'verify_state_variables',string1,source,revision,revdate)
   endif

   ngood = ngood + 1
enddo MyLoop

end subroutine verify_state_variables

!----------------------------------------------------------------------

function initialize_filenames(filename, IO_STRING) result (file_handle)
character(len=*), intent(in) :: filename(:)
integer,          intent(in) :: IO_STRING
type(file_info_type) :: file_handle

character(len=256), allocatable :: file_array(:,:)

allocate(file_array(num_ens, num_domains))
file_array  = RESHAPE(filename,  (/num_ens,  num_domains/))

call io_filenames_init(file_handle,             &
                       ncopies       = 1,       &
                       cycling       = .false., &
                       single_file   = single_file_in, &
                       restart_files = file_array)

do imem = 1, num_ens
   write(my_base,'(A,I2)') 'outens_',    imem
   write(my_desc,'(A,I2)') 'output ens', imem
   call set_file_metadata(file_handle,                          &
                          cnum     = imem,                      &
                          fnames   = file_array(imem,:),  &
                          basename = my_base,                   &
                          desc     = my_desc)

   call set_io_copy_flag(file_handle,        &
                         cnum    = imem,     &
                         io_flag = IO_STRING)
enddo

end function

!----------------------------------------------------------------------

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
