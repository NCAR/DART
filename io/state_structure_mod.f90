! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!-------------------------------------------------------------------------------

module state_structure_mod
!> \defgroup state_structure_mod state_structure_mod
!> @{ @brief Contains the stucture of the state vector and routines to access this structure
!> 
!> Usage:
!>     call add_domain() in static_init_model()
!>         - must be called.
!>         - may be called more than once in static_init_model if a model has
!>           more than one domain (e.g. WRF).
!>     
!> The add_domain() call adds a 'domain' to the state. This may be a component in 
!> the case of XCESM or another coupled model.
!>  
!> There are three ways to add a domain (these are overloaded as add_domain):
!>    * add_domain_blank. This takes model size as an argument.
!>    * add_domain_from_file. This takes a netcdf file and a list of variables
!>    * add_domain_from_spec. This makes a skeleton structure for a domain. Dimensions
!>      for each variable must be added using add_dimension_to_variable(). This is intended
!>      to be used to create netcdf output for models like bgrid_solo that are spun up.
!>      Usage:  call add_domain_from_spec
!>              for each variable call add_dimension_to_variable() for each dimension in the variable
!>              call finished_adding_domain_from_spec
!>
!> There are optional arguements to add_domain_from_spec and add_domain_from_file. These
!> are:
!>     * kind_list     - dart kind of each variable
!>     * clamp_vals    - upper and lower bounds each variable (missing_r8 for no clamping)
!>     * update_list   - true/false. Write out the variable. Default is .true.
!>
!> The state_structure hierarchy is:
!>   state -> domain -> variable 
!> Each time add_domain is called another domain is added to the state.
!> Outside code can query with the state structure using domain_id and variable_id as
!> arguements to accessor funcitons.
!> Domain_id and variable_id are intergers currently.  It may be better to make domain_id
!> a private type so other information can be stored in the domain_id.
!>
!> Interaction with the state_stucture depends on which module is accessing the structure.
!> For example, model_mod querying the number of dimensions a variable is not concerned
!> with the unlimited dimension, but state_vector_io_mod is concerned with the unlimtied
!> dimension.
!> The variable type has an io_type inside it to hold io information
!> The io accessor functions have get_io_* in their name.
!> There is no support for diagnostic file structure at the moment.
!>
!> get_dart_vector_index() and its inverse get_model_variable_indices() link 
!> a model x,y,z to dart index. The order of the state vector is no longer under model_mod
!> control.  Beware when converting model_mods such as CAM that transform the order of
!> variables after reading from a netcdf file.  There can be many calculations in model_mod
!> that are assuming a transformed order which no longer exists.
use utilities_mod, only : E_ERR, error_handler, nc_check, do_output, &
                          to_upper
use  obs_kind_mod, only : paramname_length, get_raw_obs_kind_name, &
                          get_raw_obs_kind_index
use     types_mod, only : r8, i8, missing_r8
use      sort_mod, only : index_sort

use netcdf

implicit none

private
public :: static_init_state_type,     &
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

! diagnostic files
public :: create_diagnostic_structure, &
          end_diagnostic_structure

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=256) :: string1
character(len=256) :: string2

!-------------------------------------------------------------------------------
! global variables
!-------------------------------------------------------------------------------
integer, parameter :: max_num_domains = 10
integer, parameter :: diagnostic_domain = max_num_domains + 1 ! Need to separate this from state
logical, save      :: state_initialized = .false.

!-------------------------------------------------------------------------------
! Describes information pertaining to the IO portion of the variable information
!-------------------------------------------------------------------------------
type io_information
   private

   ! netcdf variable id 
   integer :: netcdf_id ! HK Do you want one for reading, one for writing?
   integer :: xtype     ! netCDF variable type (NF90_double, etc.)

   ! clamping variables
   logical  :: clamping = .false.    ! does variable need to be range-restricted before
   real(r8) :: minvalue = missing_r8 ! min value for clamping
   real(r8) :: maxvalue = missing_r8 ! max value for clamping
   ! logical  :: out_of_range_fail = .false. ! is out of range fatal if range-checking?
   
   ! dimension information, including unlimited dimensions
   integer :: io_numdims  = 0
   integer, dimension(NF90_MAX_VAR_DIMS) :: io_dimIds
   
   ! update information
   logical :: update = .true. ! default to update variables

   ! FUTURE VARIABLES
   ! character(len=NF90_MAX_NAME) :: long_name
   ! character(len=NF90_MAX_NAME) :: units
  
end type io_information


!-------------------------------------------------------------------------------
! state_type -> domain_type -> variable_type
! #domains      #variables     info
!-------------------------------------------------------------------------------
! Describes a variable, enough for netcdf io purposes to read and write
!-------------------------------------------------------------------------------
type variable_type
   private
    
   ! variable information
   character(len=NF90_MAX_NAME) :: varname
   integer     :: var_size = 0
   integer(i8) :: index_start = -1! first occurance of variable in state
   integer(i8) :: index_end = -1  ! last occurance of variable in state
   
   ! dimension information
   integer :: numdims  = 0 ! number of dims - excluding TIME?
   logical :: var_has_unlim = .false.
   character(len=NF90_MAX_NAME), dimension(NF90_MAX_VAR_DIMS) :: dimname
   integer,                      dimension(NF90_MAX_VAR_DIMS) :: dimlens

   ! dart information
   integer :: dart_kind
   character(len=paramname_length) :: kind_string = 'KIND_RAW_STATE_VARIABLE'

   type(io_information) :: io_info

end type variable_type

!-------------------------------------------------------------------------------
! Describes a domain
!-------------------------------------------------------------------------------
type domain_type
   private

   ! variables in domain
   integer :: num_variables
   type(variable_type), allocatable :: variable(:) 
   
   ! dimension informaion for domain
   integer :: num_unique_dims
   character(len=NF90_MAX_NAME), allocatable :: unique_dim_names(:)
   integer,                      allocatable :: unique_dim_length(:)
   integer :: unlimDimId = -1 ! initialize to no unlimited dimension
   logical :: has_unlimited = .false.

   ! number of elements in the domain
   integer(i8) :: dom_size

   ! netcdf file describing the shape of the variable
   character(len=256) :: info_file  ! do we need this?
   
end type domain_type

!-------------------------------------------------------------------------------
! Describes a state vector
!-------------------------------------------------------------------------------
type state_type
   private

   ! domains or separate files
   integer           :: num_domains = 0
   ! number of domains + 1 for diagnostic domain
   type(domain_type) :: domain(max_num_domains + 1)

   ! size of the state vector
   integer(i8) :: model_size

   !>@todo is this different to character(len=128) :: name?  Redundant?
   !character(len=1024)  :: domain_name(max_num_domains)
end type state_type

!-------------------------------------------------------------------------------
! module handle
!-------------------------------------------------------------------------------
type(state_type) :: state

! This is used when dart writes the state to diagnostic files
! It is a transformation of the state_type
! One diagnostic file contains all domains.
! It is not part of the state type
logical :: diagnostic_initialized = .false.

! debug flag for get_state_indices
logical :: debug = .false.

interface add_domain
   module procedure add_domain_blank
   module procedure add_domain_from_file
   module procedure add_domain_from_spec
end interface

interface get_index_start
   module procedure get_index_start_from_varname
   module procedure get_index_start_from_varid
end interface

interface get_index_end
   module procedure get_index_end_from_varname
   module procedure get_index_end_from_varid
end interface

contains

!-------------------------------------------------------------------------------
!> Initalizes state_type
!-------------------------------------------------------------------------------
subroutine static_init_state_type()

if(state_initialized) return !?

state%num_domains = 0
state%model_size  = 0
state_initialized = .true.

end subroutine static_init_state_type

!-------------------------------------------------------------------------------
!> Given an info_file, reads in a list var_names(num_vars) 
!> into the state_strucutre. 
!>
!> Returns a dom_id that can be used to harvest information of a particular domain
!>
!> Does this need to be a function or a subroutine?
!-------------------------------------------------------------------------------
function add_domain_from_file(info_file, num_vars, var_names, kind_list, clamp_vals, update_list) result(dom_id)

character(len=*), intent(in) :: info_file
integer,          intent(in) :: num_vars
character(len=*), intent(in) :: var_names(num_vars)
integer,          intent(in), optional :: kind_list(num_vars)
real(r8),         intent(in), optional :: clamp_vals(num_vars, 2)
logical,          intent(in), optional :: update_list(num_vars)
integer :: dom_id

integer :: ivar

if (.not. state_initialized) call static_init_state_type()

! add to domains
call assert_below_max_num_domains()
state%num_domains = state%num_domains + 1
dom_id = state%num_domains !>@todo this should be a handle.

! save information about the information file
state%domain(dom_id)%info_file = info_file

! set number of variables in this domain
state%domain(dom_id)%num_variables = num_vars

! load up the variable names
allocate(state%domain(dom_id)%variable(num_vars))

do ivar = 1, num_vars
   state%domain(dom_id)%variable(ivar)%varname = var_names(ivar)
enddo

! load up variable id's and sizes
call load_state_variable_info(state%domain(dom_id))

! load up the domain unique dimension info
call load_unique_dim_info(dom_id)

if ( present(kind_list)  )  call set_dart_kinds (dom_id, num_vars, kind_list)
if ( present(clamp_vals)  ) call set_clamping   (dom_id, num_vars, clamp_vals)
if ( present(update_list) ) call set_update_list(dom_id, num_vars, update_list)

end function add_domain_from_file

!-------------------------------------------------------------------------------
!> Defines a skeleton structure for the state structure.  Dimension can be
!> added to variables with add_dimension_to_variable.
!>
!> Returns a dom_id that can be used to harvest information of a particular domain
!-------------------------------------------------------------------------------
function add_domain_from_spec(num_vars, var_names, kind_list, clamp_vals, update_list) result(dom_id)

integer,          intent(in) :: num_vars
character(len=*), intent(in) :: var_names(num_vars)
integer,          intent(in), optional :: kind_list(num_vars)
real(r8),         intent(in), optional :: clamp_vals(num_vars, 2)
logical,          intent(in), optional :: update_list(num_vars)
integer :: dom_id

integer :: ivar

if (.not. state_initialized) call static_init_state_type()

! add to domains
call assert_below_max_num_domains()
state%num_domains = state%num_domains + 1
dom_id = state%num_domains

state%domain(dom_id)%info_file = 'NULL' ! no files

! set number of variables in this domain
state%domain(dom_id)%num_variables = num_vars
state%domain(dom_id)%num_unique_dims = 0

! load up the variable names
allocate(state%domain(dom_id)%variable(num_vars))

do ivar = 1, num_vars
   state%domain(dom_id)%variable(ivar)%varname = var_names(ivar)
enddo

if ( present(kind_list)   ) call set_dart_kinds (dom_id, num_vars, kind_list)
if ( present(clamp_vals)  ) call set_clamping   (dom_id, num_vars, clamp_vals)
if ( present(update_list) ) call set_update_list(dom_id, num_vars, update_list)

end function add_domain_from_spec

!-------------------------------------------------------------------------------
!> Add a blank domain - one variable called state, length = model_size
!-------------------------------------------------------------------------------
function add_domain_blank(model_size) result(dom_id)

integer(i8), intent(in) :: model_size
integer :: dom_id

if (.not. state_initialized) call static_init_state_type()

! add to domains
call assert_below_max_num_domains()
state%num_domains = state%num_domains + 1
dom_id = state%num_domains

! domain
state%domain(dom_id)%num_variables = 1
state%domain(dom_id)%dom_size      = model_size
state%model_size = state%model_size + model_size

! variable
allocate(state%domain(dom_id)%variable(1))
state%domain(dom_id)%variable(1)%varname   = 'state'
! state%domain(dom_id)%variable(1)%io_info%long_name = 'state'
! state%domain(dom_id)%variable(1)%io_info%units     = ''
state%domain(dom_id)%variable(1)%numdims   = 1
state%domain(dom_id)%variable(1)%io_info%io_numdims   = 1
state%domain(dom_id)%variable(1)%var_size  = model_size
state%domain(dom_id)%variable(1)%index_start = 1
state%domain(dom_id)%variable(1)%index_end = model_size

! dimension
state%domain(dom_id)%variable(1)%dimname(1) = 'model_size'
state%domain(dom_id)%variable(1)%dimlens(1) =  model_size

! load up the domain unique dimension info
state%domain(dom_id)%num_unique_dims = 1
allocate(state%domain(dom_id)%unique_dim_names(1))
allocate(state%domain(dom_id)%unique_dim_length(1))
state%domain(dom_id)%unique_dim_names(1)  = 'model_size'
state%domain(dom_id)%unique_dim_length(1) =  model_size
state%domain(dom_id)%variable(1)%io_info%io_dimids(1) = 1

end function add_domain_blank

!-------------------------------------------------------------------------------
!> Load metadata from netcdf file info state_strucutre
!-------------------------------------------------------------------------------
subroutine load_state_variable_info(domain)

type(domain_type), intent(inout) :: domain

! netcdf variables
integer :: ret, ncfile
character(len=512) :: nc_filename

nc_filename = domain%info_file

! open netcdf file - all restart files in a domain have the same info?
ret = nf90_open(nc_filename, NF90_NOWRITE, ncfile)
call nc_check(ret, 'load_state_variable_info nf90_open', trim(nc_filename))

! get the dimension id of the unlimited dimension if it exists
ret = nf90_inquire(ncfile, unlimitedDimId=domain%unlimDimId)
call nc_check(ret, 'load_state_variable_info, nf90_inquire')
if ( domain%unlimDimID /= -1 ) domain%has_unlimited = .true.

! get variable ids
call load_variable_ids(ncfile, domain)

! get all variable sizes, only readers store dimensions?
call load_variable_sizes(ncfile, domain)

! close netcdf file
ret = nf90_close(ncfile)
call nc_check(ret, 'load_state_variable_info nf90_close', trim(nc_filename))

end subroutine load_state_variable_info

!-------------------------------------------------------------------------------
!> Load netcdf variable ids
!-------------------------------------------------------------------------------
subroutine load_variable_ids(ncfile, domain)

integer,           intent(in)    :: ncfile ! netdcf file id - should this be part of the domain handle?
type(domain_type), intent(inout) :: domain

integer :: ret  ! netcdf retrun value
integer :: ivar, num_vars

num_vars = domain%num_variables

do ivar = 1, num_vars
   ! load netcdf id from variable name
   ret = nf90_inq_varid(ncfile, domain%variable(ivar)%varname,    &
                                domain%variable(ivar)%io_info%netcdf_id)

   call nc_check(ret, 'load_variable_ids, nf90_inq_var_id', &
                 trim(domain%variable(ivar)%varname))

enddo

end subroutine load_variable_ids

!-------------------------------------------------------------------------------
!> load dimension information and calculate variable and domain sizes
!-------------------------------------------------------------------------------
subroutine load_variable_sizes(ncfile, domain)

integer,           intent(in)    :: ncfile ! netdcf file id - should this be part of the domain handle?
type(domain_type), intent(inout) :: domain

integer :: ivar, jdim, num_vars, num_dims !< loop variables
integer :: variable_size, domain_size !< temporary sum variable
integer :: ret ! netcdf retrun value
integer(i8) :: index_start

domain_size = 0
num_vars    = domain%num_variables
index_start = state%model_size + 1

do ivar = 1, num_vars

   ! from netcdf id load variable dimension and ids
   ret = nf90_inquire_variable(ncfile,  domain%variable(ivar)%io_info%netcdf_id,  &
                                ndims = domain%variable(ivar)%io_info%io_numdims, &
                               dimids = domain%variable(ivar)%io_info%io_dimIds,  &
                                xtype = domain%variable(ivar)%io_info%xtype)

   call nc_check(ret, 'load_variable_sizes, inq_variable', &
                 trim(domain%variable(ivar)%varname))

   variable_size = 1
   num_dims      = domain%variable(ivar)%io_info%io_numdims

   do jdim = 1, num_dims

      ! load dimension names and lengths
      ret = nf90_inquire_dimension(ncfile, domain%variable(ivar)%io_info%io_dimIds(jdim),  &
                                    name = domain%variable(ivar)%dimname(jdim), &
                                     len = domain%variable(ivar)%dimlens(jdim))

      call nc_check(ret, 'load_variable_sizes, inq_dimension', &
                    trim(domain%variable(ivar)%dimname(jdim)))
      
      variable_size = variable_size * domain%variable(ivar)%dimlens(jdim)

   enddo

   ! subtract the unlimited dimension if it exists
   !> @todo : how to handle models with multiple unlimited dimensions?
   !>         nf90_inquire returns the first unlimited dimension id  which
   !>         is the slowest varying dimension.  For now am assuming that
   !>         there can only be one unlimited dimension. Just subtract 
   !>         to get 'spatial' dimensions.
   if ( any(domain%variable(ivar)%io_info%io_dimIds(:) == domain%unlimDimId) ) then
      domain%variable(ivar)%numdims = num_dims - 1 
      domain%variable(ivar)%var_has_unlim = .TRUE.
   else
      domain%variable(ivar)%numdims = num_dims
   endif

   domain%variable(ivar)%var_size = variable_size

   ! first and last location of variable in the state index
   domain%variable(ivar)%index_start = index_start
   domain%variable(ivar)%index_end   = index_start + variable_size - 1

   ! update counters
   domain_size = domain_size + variable_size
   index_start = domain%variable(ivar)%index_end + 1

enddo

domain%dom_size  = domain_size
state%model_size = state%model_size + domain_size

end subroutine load_variable_sizes

!-------------------------------------------------------------------------------
!> Identify the unique dimensions within a domain
!> I think you can sort dimension ids and remove duplicates
!> What if there are 0 dimensions?
!-------------------------------------------------------------------------------
subroutine load_unique_dim_info(dom_id)

integer, intent(in) :: dom_id ! domain identifier

integer, allocatable :: array_of_dimids(:)
integer, allocatable :: array_of_lengths(:)
integer, allocatable :: array_of_indices(:)
character(len=NF90_MAX_NAME), allocatable :: array_of_names(:)
logical, allocatable :: unique(:)

integer :: ivar, jdim, jdim_dom
integer :: count_dims, ndims, nvars, ndims_dom
integer :: num_unique

num_unique = 0

ndims_dom = get_domain_num_dims(dom_id)
nvars     = get_num_variables(dom_id)

allocate(array_of_dimids(ndims_dom))
allocate(array_of_names(ndims_dom))
allocate(array_of_lengths(ndims_dom))
allocate(array_of_indices(ndims_dom))
allocate(unique(ndims_dom))

count_dims = 1

do ivar = 1, nvars
   ndims = get_io_num_dims(dom_id, ivar)
   do jdim = 1, ndims
      array_of_dimids(count_dims)  = state%domain(dom_id)%variable(ivar)%io_info%io_dimIds(jdim)
      array_of_names(count_dims)   = state%domain(dom_id)%variable(ivar)%dimname(jdim)
      array_of_lengths(count_dims) = state%domain(dom_id)%variable(ivar)%dimlens(jdim)
      count_dims = count_dims + 1
   enddo
enddo

call index_sort(array_of_dimids, array_of_indices, ndims_dom)

count_dims = 1
unique(:) = .false.
unique(1) = .true.

do jdim_dom = 2, ndims_dom
   if( array_of_dimids(array_of_indices(jdim_dom)) /= &
       array_of_dimids(array_of_indices(jdim_dom-1)) ) then

      count_dims = count_dims + 1
      unique(jdim_dom) = .true.
   endif
enddo

state%domain(dom_id)%num_unique_dims = count_dims

allocate(state%domain(dom_id)%unique_dim_names(count_dims))
allocate(state%domain(dom_id)%unique_dim_length(count_dims))

count_dims = 1

do jdim_dom = 1, ndims_dom
   if(unique(jdim_dom)) then
      state%domain(dom_id)%unique_dim_names(count_dims)  = array_of_names(array_of_indices(jdim_dom))
      state%domain(dom_id)%unique_dim_length(count_dims) = array_of_lengths(array_of_indices(jdim_dom))
      count_dims = count_dims + 1
   endif
enddo

deallocate(array_of_dimids, array_of_names, array_of_lengths, array_of_indices, unique)

end subroutine load_unique_dim_info

!-------------------------------------------------------------------------------
!> Returns the number of domains being used in the state structure
!-------------------------------------------------------------------------------
function get_num_domains()

integer :: get_num_domains

get_num_domains = state%num_domains

end function get_num_domains

!-------------------------------------------------------------------------------
!> Returns the number of elements in the domain
!-------------------------------------------------------------------------------
function get_domain_size(dom_id)

integer, intent(in) :: dom_id
integer :: get_domain_size

get_domain_size = state%domain(dom_id)%dom_size

end function get_domain_size

!-------------------------------------------------------------------------------
!> Returns the number of variables in the domain
!-------------------------------------------------------------------------------
function get_num_variables(dom_id)

integer, intent(in) :: dom_id ! domain
integer :: get_num_variables

get_num_variables = state%domain(dom_id)%num_variables

end function get_num_variables

!-------------------------------------------------------------------------------
!> Returns the size of a variable in a specific domain
!-------------------------------------------------------------------------------
function get_variable_size(dom_id, ivar)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: ivar ! variable
integer :: get_variable_size

get_variable_size = state%domain(dom_id)%variable(ivar)%var_size

end function get_variable_size

!-------------------------------------------------------------------------------
!> Returns the number of dimensions for a variable
!-------------------------------------------------------------------------------
function get_num_dims(dom_id, ivar)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: ivar ! variable
integer :: get_num_dims

get_num_dims = state%domain(dom_id)%variable(ivar)%numdims

end function get_num_dims

!-------------------------------------------------------------------------------
!> Return and array containing the dimension lengths, excluding the UNLIMITED dim
!-------------------------------------------------------------------------------
function get_dim_lengths(dom_id, ivar)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: ivar ! variable
integer :: get_dim_lengths(state%domain(dom_id)%variable(ivar)%numdims)

integer :: num_dims

num_dims = get_num_dims(dom_id, ivar)

get_dim_lengths(1:num_dims) = state%domain(dom_id)%variable(ivar)%dimlens(1:num_dims)

end function get_dim_lengths

!-------------------------------------------------------------------------------
!> Returns the variable name
!-------------------------------------------------------------------------------
function get_variable_name(dom_id, ivar)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: ivar ! variable
character(len=NF90_MAX_NAME)  :: get_variable_name

get_variable_name = state%domain(dom_id)%variable(ivar)%varname

end function get_variable_name

!-------------------------------------------------------------------------------
! Should you even use this?
!-------------------------------------------------------------------------------
subroutine set_var_id(dom_id, ivar, new_varid)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: ivar ! variable
integer, intent(in) :: new_varid

state%domain(dom_id)%variable(ivar)%io_info%netcdf_id = new_varid

end subroutine set_var_id

!-------------------------------------------------------------------------------
!> Returns the dimension name
!-------------------------------------------------------------------------------
function get_dim_name(dom_id, ivar, jdim)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: ivar ! variable
integer, intent(in) :: jdim ! dimension
character(len=NF90_MAX_NAME) :: get_dim_name

get_dim_name = state%domain(dom_id)%variable(ivar)%dimname(jdim)

end function get_dim_name

!-------------------------------------------------------------------------------
!> Return dimension length
!-------------------------------------------------------------------------------
function get_dim_length(dom_id, ivar, jdim)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: ivar ! variable
integer, intent(in) :: jdim ! dimension
integer :: get_dim_length

get_dim_length = state%domain(dom_id)%variable(ivar)%dimlens(jdim)

end function get_dim_length

!-------------------------------------------------------------------------------
!> Return the number of dimensions in a domain
!> Repeat dimensions are allowed
!-------------------------------------------------------------------------------
function get_domain_num_dims(dom_id)

integer, intent(in) :: dom_id ! domain identifier
integer :: get_domain_num_dims

integer :: ivar, num_vars

num_vars = state%domain(dom_id)%num_variables

get_domain_num_dims = 0

do ivar = 1, num_vars
  get_domain_num_dims =  get_domain_num_dims + &
                            state%domain(dom_id)%variable(ivar)%io_info%io_numdims
enddo

end function get_domain_num_dims

!-------------------------------------------------------------------------------
!> Return the number of unique dimensions
!-------------------------------------------------------------------------------
function get_io_num_unique_dims(dom_id)

integer, intent(in) :: dom_id ! domain identifier
integer :: get_io_num_unique_dims

get_io_num_unique_dims = state%domain(dom_id)%num_unique_dims

end function get_io_num_unique_dims

!-------------------------------------------------------------------------------
!> Return the unique dimension names
!-------------------------------------------------------------------------------
function get_io_unique_dim_name(dom_id, ivar)

integer, intent(in) :: dom_id ! domain identifier
integer, intent(in) :: ivar ! index into array, not connected to dimId
character(len=NF90_MAX_NAME) :: get_io_unique_dim_name

get_io_unique_dim_name = state%domain(dom_id)%unique_dim_names(ivar)

end function get_io_unique_dim_name

!-------------------------------------------------------------------------------
!> Return the unique dimension lengths
!-------------------------------------------------------------------------------
function get_io_unique_dim_length(dom_id, ivar)

integer, intent(in) :: dom_id ! domain identifier
integer, intent(in) :: ivar ! index into array, not connected to dimId
integer :: get_io_unique_dim_length

get_io_unique_dim_length = state%domain(dom_id)%unique_dim_length(ivar)

end function get_io_unique_dim_length

!-------------------------------------------------------------------------------
!> Returns the starting dart index for variable
!-------------------------------------------------------------------------------
function get_index_start_from_varname(dom_id, varname)

integer,          intent(in) :: dom_id
character(len=*), intent(in) :: varname
integer(i8) :: get_index_start_from_varname

integer :: var_id

var_id = get_varid_from_varname(dom_id, varname)

get_index_start_from_varname = state%domain(dom_id)%variable(var_id)%index_start

end function get_index_start_from_varname

!-------------------------------------------------------------------------------
!> Returns the starting dart index for variable
!-------------------------------------------------------------------------------
function get_index_start_from_varid(dom_id, ivar)

integer, intent(in) :: dom_id
integer, intent(in) :: ivar
integer(i8) :: get_index_start_from_varid

get_index_start_from_varid = state%domain(dom_id)%variable(ivar)%index_start

end function get_index_start_from_varid

!-------------------------------------------------------------------------------
!> Returns the starting dart index for variable
!-------------------------------------------------------------------------------
function get_index_end_from_varname(dom_id, varname)

integer,          intent(in) :: dom_id
character(len=*), intent(in) :: varname
integer(i8) :: get_index_end_from_varname

integer :: var_id

var_id = get_varid_from_varname(dom_id, varname)

get_index_end_from_varname = state%domain(dom_id)%variable(var_id)%index_end

end function get_index_end_from_varname

!-------------------------------------------------------------------------------
!> Returns the ending dart index for variable
!-------------------------------------------------------------------------------
function get_index_end_from_varid(dom_id, ivar)

integer, intent(in) :: dom_id
integer, intent(in) :: ivar
integer(i8) :: get_index_end_from_varid

get_index_end_from_varid= state%domain(dom_id)%variable(ivar)%index_end

end function get_index_end_from_varid

!-------------------------------------------------------------------------------
!> Return unlimited dimension id
!-------------------------------------------------------------------------------
function get_unlimited_dimid(dom_id)

integer, intent(in) :: dom_id
integer :: get_unlimited_dimid

get_unlimited_dimid = state%domain(dom_id)%unlimDimId

end function get_unlimited_dimid

!-------------------------------------------------------------------------------
!> Adding space for an unlimited dimension in the dimesion arrays
!> The unlimited dimension needs to be last in the list for def_var
!-------------------------------------------------------------------------------
subroutine add_time_unlimited(unlimited_dimId)

integer, intent(in)  :: unlimited_dimId

call error_handler(E_ERR, 'add_time_unlimited', 'does not exist')

end subroutine add_time_unlimited

!-------------------------------------------------------------------------------
!> Returns the number of variables below start_var
!-------------------------------------------------------------------------------
function get_sum_variables_below(start_var, dom_id)

integer, intent(in) :: start_var
integer, intent(in) :: dom_id
integer :: get_sum_variables_below

integer :: ivar, idom

get_sum_variables_below = 0

do idom = 1, dom_id -1
   get_sum_variables_below = get_sum_variables_below + &
                                state%domain(idom)%dom_size
enddo

do ivar = 1, start_var -1
   get_sum_variables_below = get_sum_variables_below + &
                                state%domain(dom_id)%variable(ivar)%var_size
enddo

end function get_sum_variables_below

!-------------------------------------------------------------------------------
!> Returns the number of variables in range of start_var and end_var
!-------------------------------------------------------------------------------
function get_sum_variables(start_var, end_var, dom_id)

integer, intent(in) :: start_var, end_var
integer, intent(in) :: dom_id
integer             :: get_sum_variables

integer :: ivar

get_sum_variables = 0

do ivar = start_var, end_var
   get_sum_variables = get_sum_variables + &
                          state%domain(dom_id)%variable(ivar)%var_size
enddo

end function get_sum_variables

!-------------------------------------------------------------------------------
!> Given a dart state index, return the iloc, jloc, kloc location of the local variable
!-------------------------------------------------------------------------------
subroutine get_model_variable_indices(index_in, iloc, jloc, kloc, var_id,  dom_id, kind_index, kind_string)

integer(i8),      intent(in)  :: index_in
integer,          intent(out) :: iloc
integer,          intent(out) :: jloc
integer,          intent(out) :: kloc
integer,          intent(out), optional :: var_id
integer,          intent(out), optional :: dom_id
integer,          intent(out), optional :: kind_index
character(len=*), intent(out), optional :: kind_string

integer     :: local_ind
integer(i8) :: index_abs, index_start, index_end
integer     :: ndims, nvars, ndomains
integer     :: idom, ivar
integer     :: domid, varid
integer, allocatable :: dsize(:)

character(len=512) :: string1

! identity obs come in with a negative value - absolute
! value into the state vector.
index_abs = abs(index_in)

varid    = -1
ndomains = get_num_domains()

! get the local variable indicies
FindVariable : do idom = 1, ndomains
   nvars = get_num_variables(idom)
   do ivar = 1, nvars

      index_start = get_index_start(idom,ivar)
      index_end = get_index_end(idom,ivar)

      if ((index_abs .ge. index_start) .and. &
          (index_abs .le. index_end) ) then

         local_ind = index_abs - index_start
         varid     = ivar
         domid     = idom
         if(debug) print*, 'index_start, index_end', index_start, index_end 

         exit FindVariable

      endif
   enddo
enddo FindVariable

if( varid == -1 ) then
   write(string1,*) 'Problem, cannot find base offset, index_in is: ', index_in
   call error_handler(E_ERR, 'get_model_variable_indices',string1)
endif

ndims = get_num_dims(domid, varid)
allocate(dsize(ndims))
dsize = get_dim_lengths(domid,varid)

if(debug) print*, 'ndims', ndims, 'sizes', dsize(1:ndims)

! unfold from variable index
if (ndims == 1) then
   kloc = 1
   jloc = 1
   iloc =  local_ind + 1
elseif (ndims == 2) then
   kloc = 1
   jloc = local_ind/dsize(1) + 1
   iloc = local_ind - (jloc-1)*dsize(1) + 1
elseif (ndims == 3) then
   kloc =  local_ind / (dsize(1)*dsize(2)) + 1
   jloc = (local_ind - (kloc-1)*dsize(1)*dsize(2))/dsize(1) + 1
   iloc =  local_ind - (kloc-1)*dsize(1)*dsize(2) - (jloc-1)*dsize(1) + 1
else
   write(string1,*) 'can not calculate indices for variable ', &
        trim(get_variable_name(domid, varid)), ' ndims = ', ndims
   call error_handler(E_ERR, 'get_model_variable_indices',string1)
endif

if (debug) then
   write(*,'(A)')            '---------------------------------------------'
   write(*,'(A, I10)')       'index_in = ', index_in
   write(*,'(A, I2, A, I2)') '   domid = ', domid, '   varid = ', varid
   write(*,'(A, 3I10, A)')   '   (i, j, k) = (', iloc, jloc, kloc, ')'
   write(*,'(A)')            '---------------------------------------------'
endif

if (present(var_id)) &
   var_id = varid

if (present(dom_id)) &
   dom_id = domid

if (present(kind_string)) &
   kind_string = get_kind_string(domid,varid)

if (present(kind_index)) &
   kind_index = get_kind_index(domid,varid)

deallocate(dsize)

end subroutine get_model_variable_indices

!-------------------------------------------------------------------------------
!> Returns the dart index from local variable indices
!-------------------------------------------------------------------------------
function get_dart_vector_index(iloc, jloc, kloc, dom_id, var_id)

integer, intent(in) :: iloc, jloc, kloc
integer, intent(in) :: dom_id, var_id
integer(i8) :: get_dart_vector_index

integer :: ndims, offset
integer :: dsize(NF90_MAX_VAR_DIMS)
character(len=512) :: string1

ndims = get_num_dims(dom_id, var_id)

dsize(1:ndims) = get_dim_lengths(dom_id,var_id)

offset = get_index_start(dom_id, var_id)

if (ndims == 1) then
   get_dart_vector_index = offset + iloc - 1
else if(ndims == 2) then
   get_dart_vector_index = offset + iloc + (jloc-1)*dsize(1) - 1
else if(ndims == 3) then
   get_dart_vector_index = offset + iloc + (jloc-1)*dsize(1) + (kloc-1)*dsize(1)*dsize(2) - 1
else
   write(string1,*) 'Can not calculate indices for variable ', &
     trim(get_variable_name(dom_id, var_id)), ' ndims = ', ndims
   call error_handler(E_ERR, 'get_dart_vector_index',string1)
endif

end function get_dart_vector_index

!-------------------------------------------------------------------------------
!> Set clamping bounds for domain variables.
!>   missing_r8 values are used to set no-clamping
!-------------------------------------------------------------------------------
subroutine set_clamping(dom_id, num_vars, clamp_vals)

integer,  intent(in) :: dom_id
integer,  intent(in) :: num_vars
real(r8), intent(in) :: clamp_vals(num_vars, 2)

real(r8) :: min_value, max_value
integer  :: ivar

do ivar = 1, num_vars

   min_value = clamp_vals(ivar, 1)
   if (min_value /= missing_r8) then 
      state%domain(dom_id)%variable(ivar)%io_info%clamping = .true.   
      state%domain(dom_id)%variable(ivar)%io_info%minvalue = min_value
   endif

   max_value = clamp_vals(ivar, 2)
   if (max_value /= missing_r8) then 
      state%domain(dom_id)%variable(ivar)%io_info%clamping = .true.   
      state%domain(dom_id)%variable(ivar)%io_info%maxvalue = max_value
   endif

enddo

end subroutine set_clamping

!-------------------------------------------------------------------------------
!> Return clamping maximum for a given variable
!-------------------------------------------------------------------------------
function get_io_clamping_maxval(dom_id, var_id)

integer, intent(in) :: dom_id
integer, intent(in) :: var_id
real(r8) :: get_io_clamping_maxval

get_io_clamping_maxval= state%domain(dom_id)%variable(var_id)%io_info%maxvalue

end function get_io_clamping_maxval

!-------------------------------------------------------------------------------
!> Return clamping minimum for a given variable
!-------------------------------------------------------------------------------
function get_io_clamping_minval(dom_id, var_id)

integer, intent(in) :: dom_id
integer, intent(in) :: var_id
real(r8) :: get_io_clamping_minval

get_io_clamping_minval= state%domain(dom_id)%variable(var_id)%io_info%minvalue

end function get_io_clamping_minval

!-------------------------------------------------------------------------------
!> Returns whether a variable should be clamped or not
!-------------------------------------------------------------------------------
function do_io_clamping(dom_id, var_id)

integer, intent(in) :: dom_id ! domain identifier
integer, intent(in) :: var_id
logical :: do_io_clamping

do_io_clamping = state%domain(dom_id)%variable(var_id)%io_info%clamping

end function do_io_clamping

!-------------------------------------------------------------------------------
!> Used to add dimenions to a variable.
!> This allows the model to add meta data to state structure so a netcdf restart 
!> can be created which has the T,U,V etc. from a cold start (no existing netcdf
!> info file) e.g. the bgrid model. The number of variables has already been given 
!> in add_domain_from_spec.
!-------------------------------------------------------------------------------
subroutine add_dimension_to_variable(dom_id, var_id, dim_name, dim_size)

integer,          intent(in) :: dom_id
integer,          intent(in) :: var_id ! this is the order you gave in add_domain
character(len=*), intent(in) :: dim_name
integer,          intent(in) :: dim_size

integer :: d_new ! dimension you are adding

if (.not. state_initialized) then
   call error_handler(E_ERR, 'state_structure_mod', 'trying to add a dimension ',&
          'before initializing calling add_domain')
endif

state%domain(dom_id)%variable(var_id)%numdims = &
    state%domain(dom_id)%variable(var_id)%numdims + 1

state%domain(dom_id)%variable(var_id)%io_info%io_numdims = &
    state%domain(dom_id)%variable(var_id)%io_info%io_numdims + 1

d_new = state%domain(dom_id)%variable(var_id)%numdims

state%domain(dom_id)%variable(var_id)%dimname(d_new) = dim_name
state%domain(dom_id)%variable(var_id)%dimlens(d_new) = dim_size

end subroutine add_dimension_to_variable


!-------------------------------------------------------------------------------
! Finalize routine when using add_domain_from_spec
!-------------------------------------------------------------------------------
subroutine finished_adding_domain(dom_id)

integer, intent(in) :: dom_id ! domain identifier

integer :: ivar, jdim
integer :: num_vars, num_dims, variable_size
integer :: next_start, count_dims

state%domain(dom_id)%num_unique_dims = get_domain_num_dims(dom_id)

!> @todo sort by character to get unique dimensions. nancy has written this, just needs testing.

allocate(state%domain(dom_id)%unique_dim_names(get_domain_num_dims(dom_id)))
allocate(state%domain(dom_id)%unique_dim_length(get_domain_num_dims(dom_id)))

state%domain(dom_id)%dom_size = 0
count_dims = 1
next_start = 1

num_vars = get_num_variables(dom_id)
do ivar = 1, num_vars
   variable_size = 1
   num_dims = get_io_num_dims(dom_id, ivar)

   do jdim = 1, num_dims
      ! product of dimensions to get variable size
      variable_size = variable_size*get_dim_length(dom_id,ivar,jdim)

      ! unique dimensions - ignoring the fact that they are not unique
      state%domain(dom_id)%unique_dim_names(count_dims)  = get_dim_name(dom_id,ivar,jdim)
      state%domain(dom_id)%unique_dim_length(count_dims) = get_dim_length(dom_id,ivar,jdim)
      count_dims = count_dims + 1

   enddo

   state%domain(dom_id)%variable(ivar)%var_size = variable_size
   state%domain(dom_id)%dom_size = get_domain_size(dom_id) + variable_size

   state%domain(dom_id)%variable(ivar)%index_start = next_start 
   state%domain(dom_id)%variable(ivar)%index_end = next_start + variable_size -1
   next_start = next_start + variable_size

enddo

state%model_size = state%model_size + get_domain_size(dom_id)

end subroutine finished_adding_domain
 
!-------------------------------------------------------------------------------
!> Print information in the state structure
!-------------------------------------------------------------------------------
subroutine state_structure_info(dom_id)

integer, intent(in) :: dom_id ! domain identifier

integer :: ivar, jdim
integer :: num_vars
integer :: num_dims
integer :: array_ids(NF90_MAX_VAR_DIMS)
integer :: array_lengths(NF90_MAX_VAR_DIMS)
character(len=NF90_MAX_VAR_DIMS) :: dim_name

if ( .not. do_output() ) return

num_vars = get_num_variables(dom_id)

do ivar = 1, num_vars
   write(*,*) ' '
   write(*,*) 'varname     : ', trim(get_variable_name(dom_id,ivar))
   write(*,*) 'var_size    : ', get_variable_size(dom_id,ivar)
   write(*,*) 'index_start : ', get_index_start(dom_id,ivar)
   write(*,*) 'index_end   : ', get_index_end(dom_id,ivar)
   write(*,*) 'kind_string : ', get_kind_string(dom_id,ivar)
   write(*,'(A,I3)') ' dart_kind   : ', get_kind_index(dom_id,ivar)
   write(*,*) 'clamping    : ', do_io_clamping(dom_id,ivar)
   write(*,*) 'minvalue    : ', get_io_clamping_minval(dom_id,ivar)
   write(*,*) 'maxvalue    : ', get_io_clamping_maxval(dom_id,ivar)
   write(*,*) 'update      : ', do_io_update(dom_id,ivar)
   
   num_dims = get_num_dims(dom_id,ivar)
   write(*,*) 'numdims     : ', num_dims

   array_lengths(1:num_dims) = get_dim_lengths(dom_id,ivar)
   do jdim = 1, num_dims
       dim_name = get_dim_name(dom_id, ivar, jdim)
       write(*,'("  dim_id[",I2,"] ",A15," , length = ",I6)') array_ids(jdim), &
                                                              trim(dim_name), &
                                                              array_lengths(jdim)
   enddo

   num_dims = get_io_num_dims(dom_id,ivar)
   write(*,*) 'io_numdims     : ', num_dims

   array_ids(1:num_dims)     = get_io_dim_ids(dom_id,ivar)
   array_lengths(1:num_dims) = get_io_dim_lengths(dom_id,ivar)
   do jdim = 1, num_dims
       dim_name = get_dim_name(dom_id, ivar, jdim)
       write(*,'("  dim_id[",I2,"] ",A15," , length = ",I6)') array_ids(jdim), &
                                                              trim(dim_name), &
                                                              array_lengths(jdim)
   enddo

   write(*,*) ' '
   
enddo

end subroutine state_structure_info

!-------------------------------------------------------------------------------
!> Set DART kinds for domain variables
!-------------------------------------------------------------------------------
subroutine set_dart_kinds(dom_id, num_vars, kind_list)

integer, intent(in) :: dom_id
integer, intent(in) :: num_vars
integer, intent(in) :: kind_list(num_vars)

integer :: ivar

do ivar = 1,num_vars
   state%domain(dom_id)%variable(ivar)%dart_kind = kind_list(ivar)
   state%domain(dom_id)%variable(ivar)%kind_string = get_raw_obs_kind_name(kind_list(ivar))
enddo

end subroutine set_dart_kinds

!-------------------------------------------------------------------------------
!> Returns the variable dart kind index
!-------------------------------------------------------------------------------
function get_kind_index(dom_id, var_id)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: var_id ! variable
integer  :: get_kind_index

get_kind_index = state%domain(dom_id)%variable(var_id)%dart_kind

end function get_kind_index

!-------------------------------------------------------------------------------
!> Returns the variable dart kind string
!-------------------------------------------------------------------------------
function get_kind_string(dom_id, var_id)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: var_id ! variable
character(len=paramname_length)  :: get_kind_string

get_kind_string = state%domain(dom_id)%variable(var_id)%kind_string

end function get_kind_string

!-------------------------------------------------------------------------------
!> Return variable id given a domain number and variable name
!-------------------------------------------------------------------------------
function get_varid_from_varname(dom_id, varname) result(var_id)

integer,          intent(in) :: dom_id
character(len=*), intent(in) :: varname
integer :: var_id

integer ivar, num_vars

num_vars = get_num_variables(dom_id)

do ivar = 1, num_vars
   if ( trim(varname) == trim(get_variable_name(dom_id,ivar)) ) then
      var_id = ivar
      return
   endif
enddo

var_id = -1

end function get_varid_from_varname

!-------------------------------------------------------------------------------
!> Return variable id given a dart kind string
!-------------------------------------------------------------------------------
function get_varid_from_kind(dom_id, dart_kind_index) result(var_id)

integer, intent(in)  :: dom_id
integer, intent(in)  :: dart_kind_index
integer :: var_id

integer :: ivar, num_vars

var_id = -1

num_vars = get_num_variables(dom_id)
do ivar = 1, num_vars
   if ( dart_kind_index == get_kind_index(dom_id, ivar) ) then
      var_id = ivar
      return
   endif
enddo

end function get_varid_from_kind

!-------------------------------------------------------------------------------
!> Return variable id given a dart kind string
!-------------------------------------------------------------------------------
subroutine get_varids_from_kind(dom_id, dart_kind_index, varid_table)

integer, intent(in)  :: dom_id
integer, intent(in)  :: dart_kind_index
integer, intent(out) :: varid_table(:)

integer :: ivar, indx
integer :: num_vars

if ( size(varid_table) < get_num_varids_from_kind(dom_id, dart_kind_index) ) then
   write(string1,*) 'Found ', get_num_varids_from_kind(dom_id, dart_kind_index)
   write(string2,*) 'varid_table must be atleast this size, ', &
                    'you have size(varid_table) = ', size(varid_table)
   call error_handler(E_ERR,'get_varids_from_kind', string1, source, revision, revdate, text2=string2)
endif

! initalize to no variables found
varid_table(:) = -1

indx = 1
num_vars = get_num_variables(dom_id)
do ivar = 1, num_vars
   if ( dart_kind_index == get_kind_index(dom_id, ivar) ) then
      varid_table(indx) = ivar
      indx = indx + 1
   endif
enddo

end subroutine get_varids_from_kind

!-------------------------------------------------------------------------------
!> Set update information for domain variables
!>   UPDATE, NO_COPY_BACK
!-------------------------------------------------------------------------------
subroutine set_update_list(dom_id, num_vars, update_var_list)

integer, intent(in) :: dom_id
integer, intent(in) :: num_vars
logical, intent(in) :: update_var_list(num_vars)

state%domain(dom_id)%variable(:)%io_info%update = update_var_list

end subroutine set_update_list

!-------------------------------------------------------------------------------
!> Returns the number of dimensions for a variable
!-------------------------------------------------------------------------------
function get_io_num_dims(dom_id, ivar)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: ivar ! variable
integer :: get_io_num_dims

get_io_num_dims = state%domain(dom_id)%variable(ivar)%io_info%io_numdims

end function get_io_num_dims

!-------------------------------------------------------------------------------
!> Return an array of netcdf dimension ids
!-------------------------------------------------------------------------------
function get_io_dim_ids(dom_id, ivar)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: ivar ! variable
integer :: get_io_dim_ids(state%domain(dom_id)%variable(ivar)%io_info%io_numdims)

integer :: num_dims

num_dims = get_io_num_dims(dom_id,ivar) 
get_io_dim_ids = state%domain(dom_id)%variable(ivar)%io_info%io_dimIds(1:num_dims)

end function get_io_dim_ids

!-------------------------------------------------------------------------------
!> Return and array containing the dimension lengths
!-------------------------------------------------------------------------------
function get_io_dim_lengths(dom_id, ivar)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: ivar ! variable
integer :: get_io_dim_lengths(state%domain(dom_id)%variable(ivar)%io_info%io_numdims)

integer :: num_dims

num_dims = get_io_num_dims(dom_id,ivar)
get_io_dim_lengths = state%domain(dom_id)%variable(ivar)%dimlens(1:num_dims)

end function get_io_dim_lengths

!-------------------------------------------------------------------------------
!> Returns whether a variable should be updated or not
!-------------------------------------------------------------------------------
function do_io_update(dom_id, var_id)

integer, intent(in) :: dom_id ! domain identifier
integer, intent(in) :: var_id
logical :: do_io_update

do_io_update = state%domain(dom_id)%variable(var_id)%io_info%update

end function do_io_update

!-------------------------------------------------------------------------------
!> Returns number of variables that match a dart kind index within a domain
!-------------------------------------------------------------------------------
function get_num_varids_from_kind(dom_id, dart_kind_index) result(num_found)

integer, intent(in) :: dom_id
integer, intent(in) :: dart_kind_index
integer :: num_found

integer :: ivar

num_found = 0

do ivar = 1, get_num_variables(dom_id)
   if ( dart_kind_index == get_kind_index(dom_id,ivar) ) then
      num_found = num_found + 1
   endif
enddo

end function get_num_varids_from_kind

!-------------------------------------------------------------------------------
!> Assert that adding a domain will not cause the maximum number of domains
!> to be exceeded.
subroutine assert_below_max_num_domains()

if (state%num_domains + 1 > max_num_domains) then
   call error_handler(E_ERR, 'state_structure_mod max number of domains exceeded', &
                      'increase parameter max_num_domains')
endif

end subroutine

!-------------------------------------------------------------------------------
!> Transform the state type to a domain that can be used to write diagnostic 
!> files.
!> One domain for the whole state. If there is more than one domain:
!>   Need to change variable names
!>   Need to change dimension names
!-------------------------------------------------------------------------------
function create_diagnostic_structure() result (diag_id)

integer :: diag_id

integer :: i, j, k, var ! loop variables
character(len = 8)   :: dom_str = ''

diag_id = diagnostic_domain

if (diagnostic_initialized) return

diagnostic_initialized = .true.

! Find total number of variables in state
state%domain(diag_id)%num_variables = 0
do i = 1, state%num_domains
   state%domain(diag_id)%num_variables = state%domain(diag_id)%num_variables + state%domain(i)%num_variables
enddo

allocate(state%domain(diag_id)%variable(state%domain(diag_id)%num_variables))

! Loop around each domain in the state and add that domain's variables
! to the diagnostic domain
var = 0
do i = 1, state%num_domains
   do j = 1, state%domain(i)%num_variables
      var = var + 1

      ! Add variable to diagnostic domain
      state%domain(diag_id)%variable(var) = state%domain(i)%variable(j)
      ! Change variable name
      if (state%num_domains > 1) then
         write(dom_str, '(A, i2.2)') '_d', i
         state%domain(diag_id)%variable(var)%varname = &
             trim(state%domain(diag_id)%variable(var)%varname) // trim(dom_str)
         do k = 1, state%domain(i)%variable(j)%numdims
            state%domain(diag_id)%variable(var)%dimname(k) = &
              trim(state%domain(diag_id)%variable(var)%dimname(k)) //trim(dom_str)
         enddo
      endif
   enddo
enddo

state%domain(diag_id)%dom_size = sum(state%domain(1:state%num_domains)%dom_size)

end function create_diagnostic_structure
!-------------------------------------------------------------------------------
!> Clean up the diagnostic structure
subroutine end_diagnostic_structure()

if (diagnostic_initialized) then
   diagnostic_initialized = .false.
   deallocate(state%domain(diagnostic_domain)%variable)
   state%domain(diagnostic_domain)%num_variables = 0
   state%domain(diagnostic_domain)%dom_size = 0

endif

end subroutine end_diagnostic_structure

!-------------------------------------------------------------------------------

!> @}
end module state_structure_mod
