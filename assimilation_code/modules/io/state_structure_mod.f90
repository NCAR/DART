! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

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
!>
!>    * add_domain_from_file. This takes a netcdf file and a list of variables
!>
!>    * add_domain_from_spec. This makes a skeleton structure for a domain. Dimensions
!>      for each variable must be added using add_dimension_to_variable(). This is intended
!>      to be used to create netcdf output for models like bgrid_solo that are spun up.
!>      Usage:  call add_domain_from_spec
!>              for each variable call add_dimension_to_variable() for each dimension in the variable
!>              call finished_adding_domain_from_spec
!>
!> There are optional arguments to add_domain_from_spec and add_domain_from_file. These
!> are:
!>     * kind_list     - dart kind of each variable
!>     * clamp_vals    - upper and lower bounds each variable (missing_r8 for no clamping)
!>     * update_list   - true/false. Write out the variable. Default is .true.
!>
!> The state_structure hierarchy is:
!>   state -> domain -> variable 
!> Each time add_domain is called another domain is added to the state.
!> Outside code can query with the state structure using domain_id and variable_id as
!> arguments to accessor funcitons.
!> Domain_id and variable_id are integers currently.  It may be better to make domain_id
!> a private type so other information can be stored in the domain_id.
!>
!> Interaction with the state_stucture depends on which module is accessing the structure.
!> For example, model_mod querying the number of dimensions a variable is not concerned
!> with the unlimited dimension, but state_vector_io_mod is concerned with the unlimtied
!> dimension.
!> The variable type has an io_type inside it to hold io information
!> The io accessor functions have get_io_* in their name.
!>
!> get_dart_vector_index() and its inverse get_model_variable_indices() link 
!> a model x,y,z to dart index. The order of the state vector is no longer under model_mod
!> control.  Beware when converting model_mods such as CAM that transform the order of
!> variables after reading from a netcdf file.  There may be calculations in model_mod
!> that are assuming a transformed order which no longer exists.

use utilities_mod, only : E_ERR, error_handler, do_output

use  obs_kind_mod, only : get_name_for_quantity, get_index_for_quantity

use     types_mod, only : r8, r4, i8, digits12, MISSING_R8, MISSING_R4, MISSING_I, &
                          obstypelength, MAX_NUM_DOMS

use netcdf_utilities_mod, only : nc_check

use      sort_mod, only : index_sort

use netcdf

implicit none

private
public :: add_domain,                 &
          get_domain_size,            &
          get_num_domains,            &
          get_variable_size,          &
          get_variable_name,          &
          get_kind_string,            &
          get_kind_index,             &
          get_varid_from_varname,     & 
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
          get_io_dim_length,          &
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
          get_xtype,                  &
          get_units,                  &
          get_long_name,              &
          get_short_name,             &
          get_has_missing_value,      &
          get_FillValue,              &
          get_missing_value,          &
          get_has_FillValue,          &
          get_add_offset,             &
          get_scale_factor,           &
          set_dart_kinds,             &
          set_clamping,               &
          set_update_list,            &
          add_dimension_to_variable,  &
          finished_adding_domain,     &
          state_structure_info,       &
          hyperslice_domain,          &
          has_unlimited_dim

! diagnostic files
!>@todo FIXME these routines are deprecated because we are no supporting 'diagnostic' 
!> files, but they will likely be useful for the single file (multiple member) input.
public :: create_diagnostic_structure, &
          end_diagnostic_structure

character(len=*), parameter :: source = 'state_structure_mod.f90'

character(len=512) :: string1, string2, string3

!-------------------------------------------------------------------------------
! global variables
!-------------------------------------------------------------------------------
integer, parameter :: diagnostic_domain = MAX_NUM_DOMS + 1 ! Need to separate this from state

!-------------------------------------------------------------------------------
! Describes information pertaining to the IO portion of the variable information
!-------------------------------------------------------------------------------
type io_information
   private

   ! netcdf variable id 
   integer :: varid ! HK Do you want one for reading, one for writing?
   integer :: xtype     ! netCDF variable type (NF90_double, etc.)

   ! clamping variables
   logical  :: clamping = .false.    ! does variable need to be range-restricted before
   real(r8) :: minvalue = missing_r8 ! min value for clamping
   real(r8) :: maxvalue = missing_r8 ! max value for clamping
   ! logical  :: out_of_range_fail = .false. ! is out of range fatal if range-checking?
   
   ! dimension information, including unlimited dimensions
   integer :: io_numdims  = 0
   integer, dimension(NF90_MAX_VAR_DIMS) :: io_dimIds
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
   
   ! update information
   logical :: update = .true. ! default to update variables
   
   ! CF-Conventions
   character(len=NF90_MAX_NAME) :: units      = ' '
   character(len=NF90_MAX_NAME) :: short_name = ' '
   character(len=NF90_MAX_NAME) :: long_name  = ' '
   logical  :: has_missing_value = .false.
   logical  :: has_FillValue     = .false.
   integer  :: missingINT   = MISSING_I   ! missing values
   real(r4) :: missingR4    = MISSING_R4
   real(r8) :: missingR8    = MISSING_R8
   integer  :: spvalINT     = MISSING_I   ! fill values
   real(r4) :: spvalR4      = MISSING_R4
   real(r8) :: spvalR8      = MISSING_R8
   real(r8) :: scale_factor = MISSING_R8
   real(r8) :: add_offset   = MISSING_R8

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
   integer     :: var_size    = 0
   !>@todo FIXME : should we have a missingI8
   integer(i8) :: index_start = -1_i8 ! first occurance of variable in state
   integer(i8) :: index_end   = -1_i8 ! last occurance of variable in state
   
   ! dimension information
   integer :: numdims  = 0 ! number of dims - excluding TIME?
   logical :: var_has_unlim = .false.
   character(len=NF90_MAX_NAME), dimension(NF90_MAX_VAR_DIMS) :: dimname
   integer,                      dimension(NF90_MAX_VAR_DIMS) :: dimlens

   ! dart information
   integer :: dart_kind = -1
   character(len=obstypelength) :: kind_string = 'unknown'

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
   
   ! dimension information for domain
   integer :: num_unique_dims
   character(len=NF90_MAX_NAME), allocatable :: unique_dim_names(:)
   integer,                      allocatable :: unique_dim_length(:)
   integer,                      allocatable :: original_dim_IDs(:)
   integer :: unlimDimId = -1 ! initialize to no unlimited dimension
   logical :: has_unlimited = .false.

   ! number of elements in the domain
   integer(i8) :: dom_size = 0_i8

   ! netcdf file describing the shape of the variable
   character(len=256) :: info_file = 'NULL'

   ! string identifying the manner in which the domain was created
   ! 'blank', 'file', or 'spec'
   character(len=6) :: method = 'none'
   
end type domain_type

!-------------------------------------------------------------------------------
! Describes a state vector
!-------------------------------------------------------------------------------
type state_type
   private

   ! domains or separate files
   integer           :: num_domains = 0
   ! number of domains + 1 for diagnostic domain
   type(domain_type) :: domain(MAX_NUM_DOMS + 1)

   ! size of the state vector
   integer(i8) :: model_size = 0_i8

   !>@todo is this different to character(len=128) :: name?  Redundant?
   !character(len=1024)  :: domain_name(MAX_NUM_DOMS)
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

interface get_missing_value
   module procedure get_missing_value_r8
   module procedure get_missing_value_r4
   module procedure get_missing_value_int
end interface

interface get_FillValue
   module procedure get_spval_r8
   module procedure get_spval_r4
   module procedure get_spval_int
end interface

contains


!-------------------------------------------------------------------------------
!> Given an info_file, reads in a list var_names(num_vars) 
!> into the state_strucutre. 
!>
!> Returns a dom_id that can be used to harvest information of a particular domain
!>
!> Does this need to be a function or a subroutine?


function add_domain_from_file(info_file, num_vars, var_names, kind_list, clamp_vals, update_list) result(dom_id)

character(len=*), intent(in) :: info_file
integer,          intent(in) :: num_vars
character(len=*), intent(in) :: var_names(num_vars)
integer,          intent(in), optional :: kind_list(num_vars)
real(r8),         intent(in), optional :: clamp_vals(num_vars, 2)
logical,          intent(in), optional :: update_list(num_vars)
integer :: dom_id

integer :: ivar

! add to domains
call assert_below_max_num_domains('add_domain_from_file')
state%num_domains = state%num_domains + 1
!>@todo dom_id should be a handle.
dom_id = state%num_domains

! save information about the information file
state%domain(dom_id)%info_file = info_file
state%domain(dom_id)%method = 'file'

! set number of variables in this domain
state%domain(dom_id)%num_variables = num_vars

! load up the variable names
allocate(state%domain(dom_id)%variable(num_vars))

do ivar = 1, num_vars
   state%domain(dom_id)%variable(ivar)%varname = var_names(ivar)
enddo

! load up variable id's and sizes
call load_state_variable_info(state%domain(dom_id),dom_id)

! load up the domain unique dimension info
call load_unique_dim_info(dom_id)

! load up any cf-conventions if they exist
call load_common_cf_conventions(state%domain(dom_id))

if ( present(kind_list)  )  call set_dart_kinds (dom_id, num_vars, kind_list)
if ( present(clamp_vals)  ) call set_clamping   (dom_id, num_vars, clamp_vals)
if ( present(update_list) ) call set_update_list(dom_id, num_vars, update_list)

end function add_domain_from_file


!-------------------------------------------------------------------------------
!> Defines a skeleton structure for the state structure.  Dimension can be
!> added to variables with add_dimension_to_variable.
!>
!> Returns a dom_id that can be used to harvest information of a particular domain


function add_domain_from_spec(num_vars, var_names, kind_list, clamp_vals, update_list) result(dom_id)

integer,          intent(in) :: num_vars
character(len=*), intent(in) :: var_names(num_vars)
integer,          intent(in), optional :: kind_list(num_vars)
real(r8),         intent(in), optional :: clamp_vals(num_vars, 2)
logical,          intent(in), optional :: update_list(num_vars)
integer :: dom_id

integer :: ivar

! add to domains
call assert_below_max_num_domains('add_domain_from_spec')
state%num_domains = state%num_domains + 1
dom_id = state%num_domains

state%domain(dom_id)%method = 'spec'

! set number of variables in this domain
state%domain(dom_id)%num_variables = num_vars
state%domain(dom_id)%num_unique_dims = 0

! load up the variable names
allocate(state%domain(dom_id)%variable(num_vars))

do ivar = 1, num_vars
   state%domain(dom_id)%variable(ivar)%varname = var_names(ivar)
   state%domain(dom_id)%variable(ivar)%io_info%xtype = NF90_DOUBLE
enddo

if ( present(kind_list)   ) call set_dart_kinds (dom_id, num_vars, kind_list)
if ( present(clamp_vals)  ) call set_clamping   (dom_id, num_vars, clamp_vals)
if ( present(update_list) ) call set_update_list(dom_id, num_vars, update_list)

end function add_domain_from_spec


!-------------------------------------------------------------------------------
!> Add a blank domain - one variable called state, length = domain_size
! HK the above comment is not true, there are three dimensions created in this function.
! HK this should set has_unlimited = .true.
function add_domain_blank(domain_size) result(dom_id)

integer(i8), intent(in) :: domain_size
integer :: dom_id

integer(i8) :: domain_offset

! add to domains
call assert_below_max_num_domains('add_domain_blank')

state%num_domains = state%num_domains + 1
dom_id = state%num_domains

if (state%num_domains > 1 ) then
   domain_offset = get_index_end(dom_id-1,get_num_variables(dom_id-1))
else
   domain_offset = 0
endif

! domain
state%domain(dom_id)%method        = 'blank'
state%domain(dom_id)%num_variables = 1
state%domain(dom_id)%dom_size      = domain_size
state%model_size                   = state%model_size + domain_size

state%domain(dom_id)%num_unique_dims = 3
allocate(state%domain(dom_id)%original_dim_IDs(3))
allocate(state%domain(dom_id)%unique_dim_names(3))
allocate(state%domain(dom_id)%unique_dim_length(3))

state%domain(dom_id)%unique_dim_names(1)  = 'location'
state%domain(dom_id)%unique_dim_names(2)  = 'member'
state%domain(dom_id)%unique_dim_names(3)  = 'time'

state%domain(dom_id)%unique_dim_length(1) =  domain_size
state%domain(dom_id)%unique_dim_length(2) =  1
state%domain(dom_id)%unique_dim_length(3) =  1

state%domain(dom_id)%original_dim_IDs(1)    =  1
state%domain(dom_id)%original_dim_IDs(2)    =  2
state%domain(dom_id)%original_dim_IDs(3)    =  NF90_UNLIMITED

! variable
allocate(state%domain(dom_id)%variable(1))

state%domain(dom_id)%variable(1)%varname            = 'state'
state%domain(dom_id)%variable(1)%numdims            = 1
state%domain(dom_id)%variable(1)%var_size           = domain_size

state%domain(dom_id)%variable(1)%index_start = domain_offset + 1
state%domain(dom_id)%variable(1)%index_end   = domain_offset + domain_size

!>@todo FIXME : what is a good default for kind_string
state%domain(dom_id)%variable(1)%kind_string = 'QTY_STATE_VARIABLE'
state%domain(dom_id)%variable(1)%dart_kind   = &
       get_index_for_quantity(state%domain(dom_id)%variable(1)%kind_string)

state%domain(dom_id)%variable(1)%dimname(1) = 'location'
state%domain(dom_id)%variable(1)%dimname(2) = 'member'
state%domain(dom_id)%variable(1)%dimname(3) = 'time'

state%domain(dom_id)%variable(1)%dimlens(1) =  domain_size
state%domain(dom_id)%variable(1)%dimlens(2) =  1
state%domain(dom_id)%variable(1)%dimlens(3) =  1

state%domain(dom_id)%variable(1)%io_info%dimlens(1) =  domain_size
state%domain(dom_id)%variable(1)%io_info%dimlens(2) =  1
state%domain(dom_id)%variable(1)%io_info%dimlens(3) =  1

state%domain(dom_id)%variable(1)%io_info%xtype        = NF90_DOUBLE
state%domain(dom_id)%variable(1)%io_info%units        = 'none'
state%domain(dom_id)%variable(1)%io_info%io_numdims   = 3
state%domain(dom_id)%variable(1)%io_info%io_dimids(1) = 1
state%domain(dom_id)%variable(1)%io_info%io_dimids(2) = 2
state%domain(dom_id)%variable(1)%io_info%io_dimids(3) = NF90_UNLIMITED

end function add_domain_blank


!-------------------------------------------------------------------------------
!> Load metadata from netcdf file info state_strucutre


subroutine load_state_variable_info(domain, domain_index)

type(domain_type), intent(inout) :: domain
integer,           intent(in)    :: domain_index

! netcdf variables
integer :: ret, ncfile
character(len=256) :: nc_filename

nc_filename = domain%info_file

! open netcdf file - all restart files in a domain have the same info?
ret = nf90_open(nc_filename, NF90_NOWRITE, ncfile)
call nc_check(ret, 'load_state_variable_info nf90_open', trim(nc_filename))

! get the dimension id of the unlimited dimension if it exists
ret = nf90_inquire(ncfile, unlimitedDimId=domain%unlimDimId)
call nc_check(ret, 'load_state_variable_info, nf90_inquire')
if ( domain%unlimDimID /= -1 ) domain%has_unlimited = .true.

! get variable ids
call load_variable_ids(ncfile, domain, domain_index)

! get all variable sizes, only readers store dimensions?
call load_variable_sizes(ncfile, domain)

! close netcdf file
ret = nf90_close(ncfile)
call nc_check(ret, 'load_state_variable_info nf90_close', trim(nc_filename))

end subroutine load_state_variable_info


!-------------------------------------------------------------------------------
!> Load netcdf variable ids


subroutine load_variable_ids(ncfile, domain, domain_index)

integer,           intent(in)    :: ncfile ! netdcf file id - should this be part of the domain handle?
type(domain_type), intent(inout) :: domain
integer,           intent(in)    :: domain_index

integer :: ret  ! netcdf return value
integer :: ivar, num_vars

num_vars = domain%num_variables

do ivar = 1, num_vars
   ! load netcdf id from variable name
   ret = nf90_inq_varid(ncfile, domain%variable(ivar)%varname,    &
                                domain%variable(ivar)%io_info%varid)

   write(string1,*)'domain ',domain_index,', variable #',ivar,' "', &
       trim(domain%variable(ivar)%varname)//'" from file "'//trim(domain%info_file)//'"'
   call nc_check(ret, 'load_variable_ids, nf90_inq_var_id', string1) 

enddo

end subroutine load_variable_ids

!-------------------------------------------------------------------------------
!> load dimension information and calculate variable and domain sizes
!-------------------------------------------------------------------------------
subroutine load_variable_sizes(ncfile, domain)

integer,           intent(in)    :: ncfile ! netdcf file id - should this be part of the domain handle?
type(domain_type), intent(inout) :: domain

integer :: ivar, jdim, num_vars, num_dims !< loop variables
integer :: variable_size
integer :: ret ! netcdf return value
integer(i8) :: index_start, domain_size

domain_size = 0
num_vars    = domain%num_variables
index_start = state%model_size + 1

do ivar = 1, num_vars

   ! from netcdf id load variable dimension and ids
   ret = nf90_inquire_variable(ncfile,  domain%variable(ivar)%io_info%varid,  &
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
                                     len = domain%variable(ivar)%io_info%dimlens(jdim))

      domain%variable(ivar)%dimlens(jdim) = domain%variable(ivar)%io_info%dimlens(jdim)

      call nc_check(ret, 'load_variable_sizes, inq_dimension', &
                    trim(domain%variable(ivar)%dimname(jdim)))
      
      !>@todo FIXME we'll have to document that no user can have a 'member' dimension
      !>in their own netcdf files.  it's more reasonable to indicate 'time' is special
      !>but this code doesn't know if its a user file or one we wrote.  but without
      !>skipping these dimensions you get the wrong variable size because it includes
      !>all times and all members in the size.  this needs to be revisited later.

      if ((domain%variable(ivar)%dimname(jdim) == 'time') .or. &
          (domain%variable(ivar)%dimname(jdim) == 'member')) cycle

      variable_size = variable_size * domain%variable(ivar)%dimlens(jdim)

   enddo


   ! to be consistent this needs to ignore both 'time' and 'member' for
   ! files we write, and newer netcdf libs support multiple unlimited dims
   ! so this does need to change.  but it works, somehow, as-is so leave it
   ! for now.
   ! subtract the unlimited dimension if it exists
   !>@todo : how to handle models with multiple unlimited dimensions?
   !>         nf90_inquire returns the first unlimited dimension id  which
   !>         is the slowest varying dimension.  For now am assuming that
   !>         there can only be one unlimited dimension. Just subtract 
   !>         to get 'spatial' dimensions.
   if ( any(domain%variable(ivar)%io_info%io_dimIds(1:num_dims) == domain%unlimDimId) ) then
      domain%variable(ivar)%numdims = num_dims - 1 
      domain%variable(ivar)%var_has_unlim = .TRUE.
   else
      domain%variable(ivar)%numdims = num_dims
   endif

   ! member is not a spatial domain but could be included in a single file
   do jdim = 1, num_dims
      if ( domain%variable(ivar)%dimname(jdim) == 'member') then
         domain%variable(ivar)%numdims = domain%variable(ivar)%numdims - 1 
      endif
   enddo

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
allocate(state%domain(dom_id)%original_dim_IDs(count_dims))

count_dims = 1

do jdim_dom = 1, ndims_dom
   if(unique(jdim_dom)) then
      state%domain(dom_id)%unique_dim_names( count_dims) = array_of_names(  array_of_indices(jdim_dom))
      state%domain(dom_id)%unique_dim_length(count_dims) = array_of_lengths(array_of_indices(jdim_dom))
      state%domain(dom_id)%original_dim_IDs(   count_dims) = array_of_dimids( array_of_indices(jdim_dom))
      count_dims = count_dims + 1
   endif
enddo

deallocate(array_of_dimids, array_of_names, array_of_lengths, array_of_indices, unique)

end subroutine load_unique_dim_info

!-------------------------------------------------------------------------------
!> Check to see if the template file has any of the common cf-conventions
!>
!>  * units
!>  * short_name
!>  * long_name
!>  * short_name
!>  * _FillValue
!>  * missing_value
!>  * add_offset
!>  * scale_factor
!>
!> If they exist, load them up into the state structure.
!-------------------------------------------------------------------------------

subroutine load_common_cf_conventions(domain)

type(domain_type), intent(inout) :: domain

integer :: ivar
integer :: nvars

! netcdf variables
integer  :: ret, ncid, VarID
integer  :: var_xtype
integer  :: cf_spvalINT
real(r4) :: cf_spvalR4
real(digits12) :: cf_spvalR8
real(digits12) :: cf_scale_factor, cf_add_offset
character(len=512) :: ncFilename
character(len=NF90_MAX_NAME) :: var_name
character(len=NF90_MAX_NAME) :: cf_long_name, cf_short_name, cf_units

character(len=*), parameter :: routine = 'load_common_cf_conventions'

ncFilename = domain%info_file

ret = nf90_open(ncFilename, NF90_NOWRITE, ncid)
call nc_check(ret, routine,'nf90_open '//trim(ncFilename))

! determine attributes of each variable in turn

nvars = domain%num_variables

do ivar = 1, nvars
   var_name = domain%variable(ivar)%varname

   ret = nf90_inq_varid(ncid, trim(var_name), VarID)
   call nc_check(ret, routine, 'inq_varid '//trim(var_name))

   ! If the short_name, long_name and/or units attributes are set, get them.
   ! They are not REQUIRED by DART but are nice to keep around if they are present.

   if( nf90_inquire_attribute(    ncid, VarID, 'long_name') == NF90_NOERR ) then
      ret = nf90_get_att(ncid, VarID, 'long_name' , cf_long_name)
      call nc_check(ret, routine, 'get_att long_name '//trim(var_name))
      domain%variable(ivar)%io_info%long_name = cf_long_name
   endif

   if( nf90_inquire_attribute(    ncid, VarID, 'short_name') == NF90_NOERR ) then
      ret = nf90_get_att(ncid, VarID, 'short_name' , cf_short_name)
      call nc_check(ret, routine, 'get_att short_name '//trim(var_name))
      domain%variable(ivar)%io_info%short_name = cf_short_name
   endif

   if( nf90_inquire_attribute(    ncid, VarID, 'units') == NF90_NOERR )  then
      ret = nf90_get_att(ncid, VarID, 'units' , cf_units)
      call nc_check(ret, routine, 'get_att units '//trim(var_name))
      domain%variable(ivar)%io_info%units = cf_units
   endif

   ! Saving any FillValue, missing_value attributes ...
   ! Also stuff them into the 'r8' slots to facilitate simpler clamp_variable
   ! implementation. (Since we coerce the DART state to 'r8')

   var_xtype = domain%variable(ivar)%io_info%xtype
   select case (var_xtype)
      case ( NF90_INT )
          ! Sometimes the attributes are specified as GLOBAL attributes
          ret = nf90_get_att(ncid, NF90_GLOBAL, '_FillValue',    cf_spvalINT)
          if (ret == NF90_NOERR) then
             domain%variable(ivar)%io_info%spvalINT            = cf_spvalINT
             domain%variable(ivar)%io_info%spvalR8             = real(cf_spvalINT,r8)
             domain%variable(ivar)%io_info%has_FillValue       = .true.
          endif
          ret = nf90_get_att(ncid, NF90_GLOBAL, 'missing_value', cf_spvalINT)
          if (ret == NF90_NOERR) then
             domain%variable(ivar)%io_info%missingINT          = cf_spvalINT
             domain%variable(ivar)%io_info%spvalR8             = real(cf_spvalINT,r8)
             domain%variable(ivar)%io_info%has_missing_value   = .true.
          endif
          ! Usually the attributes are specified as variable attributes
          ret = nf90_get_att(ncid, VarID, '_FillValue',          cf_spvalINT)
          if (ret == NF90_NOERR) then
             domain%variable(ivar)%io_info%spvalINT            = cf_spvalINT
             domain%variable(ivar)%io_info%spvalR8             = real(cf_spvalINT,r8)
             domain%variable(ivar)%io_info%has_FillValue       = .true.
          endif
          ret = nf90_get_att(ncid, VarID, 'missing_value',       cf_spvalINT)
          if (ret == NF90_NOERR) then
             domain%variable(ivar)%io_info%missingINT          = cf_spvalINT
             domain%variable(ivar)%io_info%missingR8           = real(cf_spvalINT,r8)
             domain%variable(ivar)%io_info%has_missing_value   = .true.
          endif

      case ( NF90_FLOAT )
          ret = nf90_get_att(ncid, NF90_GLOBAL, '_FillValue',    cf_spvalR4)
          if (ret == NF90_NOERR) then
             domain%variable(ivar)%io_info%spvalR4             = cf_spvalR4
             domain%variable(ivar)%io_info%spvalR8             = real(cf_spvalR4,r8)
             domain%variable(ivar)%io_info%has_FillValue       = .true.
          endif
          ret = nf90_get_att(ncid, NF90_GLOBAL, 'missing_value', cf_spvalR4)
          if (ret == NF90_NOERR) then
             domain%variable(ivar)%io_info%missingR4           = cf_spvalR4
             domain%variable(ivar)%io_info%spvalR8             = real(cf_spvalR4,r8)
             domain%variable(ivar)%io_info%has_missing_value   = .true.
          endif
          ret = nf90_get_att(ncid, VarID, '_FillValue',          cf_spvalR4)
          if (ret == NF90_NOERR) then
             domain%variable(ivar)%io_info%spvalR4             = cf_spvalR4
             domain%variable(ivar)%io_info%spvalR8             = real(cf_spvalR4,r8)
             domain%variable(ivar)%io_info%has_FillValue       = .true.
          endif
          ret = nf90_get_att(ncid, VarID, 'missing_value',       cf_spvalR4)
          if (ret == NF90_NOERR) then
             domain%variable(ivar)%io_info%missingR4           = cf_spvalR4
             domain%variable(ivar)%io_info%missingR8           = real(cf_spvalR4,r8)
             domain%variable(ivar)%io_info%has_missing_value   = .true.
          endif

      case ( NF90_DOUBLE )

          ! If r8 = r4, 
          ! the missing_value must be present in both missingR4 and missingR8 
          ! ditto for _FillValue.
          ! This satisfies the overloaded operator 'get_missing_value, get_FillValue'

          ret = nf90_get_att(ncid, NF90_GLOBAL, '_FillValue',    cf_spvalR8)
          if (ret == NF90_NOERR) then
             domain%variable(ivar)%io_info%spvalR4             = cf_spvalR8
             domain%variable(ivar)%io_info%spvalR8             = cf_spvalR8
             domain%variable(ivar)%io_info%has_FillValue       = .true.
          endif
          ret = nf90_get_att(ncid, NF90_GLOBAL, 'missing_value', cf_spvalR8)
          if (ret == NF90_NOERR) then
             domain%variable(ivar)%io_info%missingR4           = cf_spvalR8
             domain%variable(ivar)%io_info%missingR8           = cf_spvalR8
             domain%variable(ivar)%io_info%has_missing_value   = .true.
          endif
          ret = nf90_get_att(ncid, VarID, '_FillValue',          cf_spvalR8)
          if (ret == NF90_NOERR) then
             domain%variable(ivar)%io_info%spvalR4             = cf_spvalR8
             domain%variable(ivar)%io_info%spvalR8             = cf_spvalR8
             domain%variable(ivar)%io_info%has_FillValue       = .true.
          endif
          ret = nf90_get_att(ncid, VarID, 'missing_value',       cf_spvalR8)
          if (ret == NF90_NOERR) then
             domain%variable(ivar)%io_info%missingR4           = cf_spvalR8
             domain%variable(ivar)%io_info%missingR8           = cf_spvalR8
             domain%variable(ivar)%io_info%has_missing_value   = .true.
          endif

      case DEFAULT
         write(string1,*) ' unsupported netcdf variable type : ', var_xtype
         call error_handler(E_ERR,routine,string1,source)
   end select

   ! If the variable has one or the other, no problem.
   ! If the variable has both _FillValue and missing_value attributes, the
   ! values must be the same or we are lost. DART only supports one missing
   ! value code. When we go to write, we have no way of knowing which value
   ! to use as a replacement for the DART missing code.

   if ( domain%variable(ivar)%io_info%has_missing_value .and. &
        domain%variable(ivar)%io_info%has_FillValue ) then

      if ( domain%variable(ivar)%io_info%missingR8 /= &
           domain%variable(ivar)%io_info%spvalR8 ) then

         write(string1, *) trim(var_name)//' missing_value /= _FillValue '
         write(string2,*) 'missing_value is ', domain%variable(ivar)%io_info%missingR8 
         write(string3,*) '_FillValue    is ', domain%variable(ivar)%io_info%spvalR8
         call error_handler(E_ERR,'set_dart_missing_value:',string1, &
                source, text2=string2, text3=string3)
      endif

   endif

   !>@todo FIXME : Not supporting scale factor or offset at the moment, so just error.
   !>              To fully support netCDF I/O we need to
   !>              pack and unpack the variable if these attributes exist.
   if (nf90_get_att(ncid, VarID, 'scale_factor',   cf_scale_factor) == NF90_NOERR) then
      domain%variable(ivar)%io_info%scale_factor = cf_scale_factor
      write(string1,*) 'scale_factor not supported at the moment'
      write(string2,*) 'contact DART if you would like to get this to work'
      call error_handler(E_ERR,routine,string1,source,text2=string2)
   endif

   if (nf90_get_att(ncid, VarID, 'add_offset', cf_add_offset) == NF90_NOERR) then
      domain%variable(ivar)%io_info%add_offset = cf_add_offset
      write(string1,*) 'add_offset not supported at the moment'
      write(string2,*) 'contact DART if you would like to get this to work'
      call error_handler(E_ERR,routine,string1,source,text2=string2)
   endif

enddo

! close netcdf file
ret = nf90_close(ncid)
call nc_check(ret, routine, 'nf90_close', trim(ncFilename))

end subroutine load_common_cf_conventions

!-------------------------------------------------------------------------------
!> TIEGCM, the top level of each variable is the boundary condition.
!> So the top level should not be part of the state.
!> How general does this need to be?  At the moment it is just a hack to
!> be able to proceed with TIEGCM.
!>
!> This slicing needs to be on the the DART state structure, not the io netcdf
!> structure.

subroutine hyperslice_domain(dom_id, dim_name, new_length)

integer,          intent(in) :: dom_id
character(len=*), intent(in) :: dim_name
integer,          intent(in) :: new_length

integer :: ivar, jdim
integer(i8) :: orig_domain_size, orig_model_size, index_start
integer(i8) :: domain_size, variable_size

orig_domain_size = get_domain_size(dom_id)
orig_model_size = state%model_size
! Start of the domain remains the same
index_start = state%domain(dom_id)%variable(1)%index_start

! Reduce the dimension size
do ivar = 1, get_num_variables(dom_id)
   do jdim = 1, get_num_dims(dom_id, ivar)
       if (get_dim_name(dom_id, ivar, jdim) == dim_name) then
          state%domain(dom_id)%variable(ivar)%dimlens(jdim) = new_length
       endif
   enddo
enddo

! need to modify state index accessing for each variable:
!   variable_size
!   index_start, index_end
!   model size
!   domain size

domain_size = 0
! variable size
do ivar = 1, get_num_variables(dom_id)

   variable_size = 1

   do jdim = 1, get_num_dims(dom_id, ivar)
      !HK there is a todo about this, what is going on with time and member?
      if ((state%domain(dom_id)%variable(ivar)%dimname(jdim) == 'time') .or. &
          (state%domain(dom_id)%variable(ivar)%dimname(jdim) == 'member')) cycle

      variable_size = variable_size *state%domain(dom_id)%variable(ivar)%dimlens(jdim)
   enddo

  state%domain(dom_id)%variable(ivar)%var_size = variable_size

  ! first and last location of variable in the state index
  state%domain(dom_id)%variable(ivar)%index_start = index_start
  state%domain(dom_id)%variable(ivar)%index_end   = index_start + variable_size - 1

  ! update counters
  domain_size = domain_size + variable_size
  index_start = state%domain(dom_id)%variable(ivar)%index_end + 1

enddo

! update state structure with new domain size and new model size
state%domain(dom_id)%dom_size  = domain_size
state%model_size = orig_model_size - orig_domain_size + domain_size

end subroutine hyperslice_domain
!-------------------------------------------------------------------------------
!> Returns the number of domains being used in the state structure


function get_num_domains()

integer :: get_num_domains

get_num_domains = state%num_domains

end function get_num_domains

!-------------------------------------------------------------------------------
!> Return whether the domain has an unlimited dimension

function has_unlimited_dim(dom_id)

integer, intent(in) :: dom_id
logical :: has_unlimited_dim

has_unlimited_dim = state%domain(dom_id)%has_unlimited

end function has_unlimited_dim

!-------------------------------------------------------------------------------
!> Returns the number of elements in the domain


function get_domain_size(dom_id)

integer, intent(in) :: dom_id
integer(i8) :: get_domain_size

call check_domain_id(dom_id,'get_domain_size')

get_domain_size = state%domain(dom_id)%dom_size

end function get_domain_size


!-------------------------------------------------------------------------------
!> Returns the number of variables in the domain


function get_num_variables(dom_id)

integer, intent(in) :: dom_id ! domain
integer :: get_num_variables

call check_domain_id(dom_id,'get_num_variables')

get_num_variables = state%domain(dom_id)%num_variables

end function get_num_variables


!-------------------------------------------------------------------------------
!> Returns the size of a variable in a specific domain


function get_variable_size(dom_id, ivar)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: ivar ! variable
integer :: get_variable_size

write(string1,*)'get_variable_size:requesting size of variable # ',ivar, ' for domain ',dom_id
call check_domain_id(dom_id, string1)

if (state%domain(dom_id)%num_variables < ivar) then
   write(string1,'(''domain '',i4,'' has '',i6,'' variable(s).'')') dom_id, &
                     state%domain(dom_id)%num_variables
   write(string2,*)'requested size of variable # ',ivar
   call error_handler(E_ERR,'get_variable_size',string1,source,text2=string2)
endif

get_variable_size = state%domain(dom_id)%variable(ivar)%var_size

end function get_variable_size


!-------------------------------------------------------------------------------
!> Returns the number of dimensions for a variable


function get_num_dims(dom_id, ivar)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: ivar ! variable
integer :: get_num_dims

get_num_dims = state%domain(dom_id)%variable(ivar)%numdims

end function get_num_dims


!-------------------------------------------------------------------------------
!> Return and array containing the dimension lengths, for the DART state
!> excluding the UNLIMITED dim

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


function get_variable_name(dom_id, ivar)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: ivar ! variable
character(len=NF90_MAX_NAME)  :: get_variable_name

get_variable_name = state%domain(dom_id)%variable(ivar)%varname

end function get_variable_name


!-------------------------------------------------------------------------------
! Should you even use this?
!>@todo no corresponding 'get', no way to check if it worked, is new_varid valid?


subroutine set_var_id(dom_id, ivar, new_varid)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: ivar ! variable
integer, intent(in) :: new_varid

state%domain(dom_id)%variable(ivar)%io_info%varid = new_varid

end subroutine set_var_id


!-------------------------------------------------------------------------------
!> Returns the dimension name


function get_dim_name(dom_id, ivar, jdim)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: ivar ! variable
integer, intent(in) :: jdim ! dimension
character(len=NF90_MAX_NAME) :: get_dim_name

get_dim_name = state%domain(dom_id)%variable(ivar)%dimname(jdim)

end function get_dim_name


!-------------------------------------------------------------------------------
!> Return dimension length


function get_dim_length(dom_id, ivar, jdim)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: ivar ! variable
integer, intent(in) :: jdim ! dimension
integer :: get_dim_length

get_dim_length = state%domain(dom_id)%variable(ivar)%dimlens(jdim)

end function get_dim_length

!-------------------------------------------------------------------------------
!> Return io dimension length

function get_io_dim_length(dom_id, ivar, jdim)
integer, intent(in) :: jdim ! dimension
integer :: get_io_dim_length

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: ivar ! variable

get_io_dim_length = state%domain(dom_id)%variable(ivar)%io_info%dimlens(jdim)

end function get_io_dim_length


!-------------------------------------------------------------------------------
!> Return the number of dimensions in a domain
!> Repeat dimensions are allowed


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


function get_io_num_unique_dims(dom_id)

integer, intent(in) :: dom_id ! domain identifier
integer :: get_io_num_unique_dims

get_io_num_unique_dims = state%domain(dom_id)%num_unique_dims

end function get_io_num_unique_dims


!-------------------------------------------------------------------------------
!> Return the unique dimension names


function get_io_unique_dim_name(dom_id, jdim)

integer, intent(in) :: dom_id ! domain identifier
integer, intent(in) :: jdim ! index into array, not connected to dimId
character(len=NF90_MAX_NAME) :: get_io_unique_dim_name

get_io_unique_dim_name = state%domain(dom_id)%unique_dim_names(jdim)

end function get_io_unique_dim_name


!-------------------------------------------------------------------------------
!> Return the unique dimension lengths


function get_io_unique_dim_length(dom_id, jdim)

integer, intent(in) :: dom_id ! domain identifier
integer, intent(in) :: jdim ! index into array, not connected to dimId
integer :: get_io_unique_dim_length

get_io_unique_dim_length = state%domain(dom_id)%unique_dim_length(jdim)

end function get_io_unique_dim_length


!-------------------------------------------------------------------------------
!> Return the original dimension ID from the source (blank,file, or spec)
!> This is intentionally not a public routine and is intended to be used
!> to summarize what is being used from the source (state_structure_info).

function get_original_dim_ID(dom_id, jdim)

integer, intent(in) :: dom_id
integer, intent(in) :: jdim
integer :: get_original_dim_ID

get_original_dim_ID = state%domain(dom_id)%original_dim_IDs(jdim)

end function get_original_dim_ID


!-------------------------------------------------------------------------------
!> Returns the starting dart index for variable


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


function get_index_start_from_varid(dom_id, ivar)

integer, intent(in) :: dom_id
integer, intent(in) :: ivar
integer(i8) :: get_index_start_from_varid

get_index_start_from_varid = state%domain(dom_id)%variable(ivar)%index_start

end function get_index_start_from_varid


!-------------------------------------------------------------------------------
!> Returns the starting dart index for variable


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


function get_index_end_from_varid(dom_id, ivar)

integer, intent(in) :: dom_id
integer, intent(in) :: ivar
integer(i8) :: get_index_end_from_varid

get_index_end_from_varid= state%domain(dom_id)%variable(ivar)%index_end

end function get_index_end_from_varid


!-------------------------------------------------------------------------------
!> Return unlimited dimension id


function get_unlimited_dimid(dom_id)

integer, intent(in) :: dom_id
integer :: get_unlimited_dimid

get_unlimited_dimid = state%domain(dom_id)%unlimDimId

end function get_unlimited_dimid


!-------------------------------------------------------------------------------
!> Adding space for an unlimited dimension in the dimesion arrays
!> The unlimited dimension needs to be last in the list for def_var
!>@todo this is a terrible name. The unlimited dimension can be for anything, not just time.


subroutine add_time_unlimited(unlimited_dimId)

integer, intent(in)  :: unlimited_dimId

call error_handler(E_ERR, 'add_time_unlimited', 'routine does not exist')

end subroutine add_time_unlimited


!-------------------------------------------------------------------------------
!> Returns the number of variables below start_var


function get_sum_variables_below(start_var, dom_id)

integer, intent(in) :: start_var
integer, intent(in) :: dom_id
integer(i8)         :: get_sum_variables_below

integer :: ivar, idom

get_sum_variables_below = 0_i8

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
!> Returns the number of elements in range of start_var and end_var


function get_sum_variables(start_var, end_var, dom_id)

integer, intent(in) :: start_var, end_var
integer, intent(in) :: dom_id
integer(i8)         :: get_sum_variables

integer :: ivar

get_sum_variables = 0_i8

do ivar = start_var, end_var
   get_sum_variables = get_sum_variables + &
                          state%domain(dom_id)%variable(ivar)%var_size
enddo

end function get_sum_variables


!-------------------------------------------------------------------------------
!> Given a dart state index, return the iloc, jloc, kloc location of the local variable


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

! identity obs come in with a negative value - absolute
! value into the state vector.
index_abs = abs(index_in)

varid    = -1
ndomains = get_num_domains()

! get the local variable indices
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
   write(string1,*) 'Can only calculate variable indices for 1<=ndims<=3 '
   write(string2,*) trim(get_variable_name(domid, varid)), 'has ndims = ', ndims
   call error_handler(E_ERR, 'get_model_variable_indices',string1, &
                      source, text2=string2)
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


function get_dart_vector_index(iloc, jloc, kloc, dom_id, var_id)

integer, intent(in) :: iloc, jloc, kloc
integer, intent(in) :: dom_id, var_id
integer(i8)         :: get_dart_vector_index

integer :: ndims
integer(i8) :: offset
integer :: dsize(NF90_MAX_VAR_DIMS)

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
!>   clamp_vals(ivar,1) must be the minimum value
!>   clamp_vals(ivar,2) must be the maximum value


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


function get_io_clamping_maxval(dom_id, var_id)

integer, intent(in) :: dom_id
integer, intent(in) :: var_id
real(r8) :: get_io_clamping_maxval

get_io_clamping_maxval= state%domain(dom_id)%variable(var_id)%io_info%maxvalue

end function get_io_clamping_maxval


!-------------------------------------------------------------------------------
!> Return clamping minimum for a given variable


function get_io_clamping_minval(dom_id, var_id)

integer, intent(in) :: dom_id
integer, intent(in) :: var_id
real(r8) :: get_io_clamping_minval

get_io_clamping_minval= state%domain(dom_id)%variable(var_id)%io_info%minvalue

end function get_io_clamping_minval


!-------------------------------------------------------------------------------
!> Returns whether a variable should be clamped or not


function do_io_clamping(dom_id, var_id)

integer, intent(in) :: dom_id ! domain identifier
integer, intent(in) :: var_id
logical :: do_io_clamping

do_io_clamping = state%domain(dom_id)%variable(var_id)%io_info%clamping

end function do_io_clamping


!-------------------------------------------------------------------------------
!> Used to add dimensions to a variable.
!> This allows the model to add meta data to state structure so a netcdf restart 
!> can be created which has the T,U,V etc. from a cold start (no existing netcdf
!> info file) e.g. the bgrid model. The number of variables has already been given 
!> in add_domain_from_spec.


subroutine add_dimension_to_variable(dom_id, var_id, dim_name, dim_size)

integer,          intent(in) :: dom_id
integer,          intent(in) :: var_id ! this is the order you gave in add_domain
character(len=*), intent(in) :: dim_name
integer,          intent(in) :: dim_size

integer :: d_new ! dimension you are adding

if ( state%domain(dom_id)%method /= 'spec') then
   write(string1,'(''domain '',i4,'' created with method "'',A,''"'')')dom_id, trim(state%domain(dom_id)%method)
   write(string2,*)'only domains created by add_domain_from_spec() may add dimensions to their variables.'
   call error_handler(E_ERR,'add_dimension_to_variable', string1, source, text2=string2)

endif

state%domain(dom_id)%variable(var_id)%numdims = &
    state%domain(dom_id)%variable(var_id)%numdims + 1

state%domain(dom_id)%variable(var_id)%io_info%io_numdims = &
    state%domain(dom_id)%variable(var_id)%io_info%io_numdims + 1

d_new = state%domain(dom_id)%variable(var_id)%numdims

state%domain(dom_id)%variable(var_id)%dimname(d_new) = dim_name
state%domain(dom_id)%variable(var_id)%dimlens(d_new) = dim_size
state%domain(dom_id)%variable(var_id)%io_info%dimlens(d_new) = dim_size
state%domain(dom_id)%variable(var_id)%io_info%io_dimids(d_new) = d_new

state%domain(dom_id)%num_unique_dims = d_new

end subroutine add_dimension_to_variable


!-------------------------------------------------------------------------------
!> Finalize routine when using add_domain_from_spec


subroutine finished_adding_domain(dom_id)

integer, intent(in) :: dom_id ! domain identifier

integer :: ivar, jdim
integer :: num_vars, num_dims, variable_size
integer(i8) :: next_start
integer :: count_dims
integer(i8) :: domain_offset

if ( state%domain(dom_id)%method /= 'spec') then
   write(string1,'(''domain '',i4,'' created with method "'',A,''"'')')dom_id, trim(state%domain(dom_id)%method)
   write(string2,*)'only domains created by add_domain_from_spec() may call finished_adding_domain.'
   call error_handler(E_ERR,'finished_adding_domain', string1, &
              source, text2=string2)

endif

state%domain(dom_id)%num_unique_dims = get_domain_num_dims(dom_id)

!>@todo sort by character to get unique dimensions. nancy has already written this, use and test

allocate(state%domain(dom_id)%unique_dim_names(get_domain_num_dims(dom_id)))
allocate(state%domain(dom_id)%unique_dim_length(get_domain_num_dims(dom_id)))

state%domain(dom_id)%dom_size = 0
count_dims = 1
next_start = 1

domain_offset = 0
if (state%num_domains > 1 ) domain_offset = get_index_end(dom_id-1,get_num_variables(dom_id-1))

num_vars = get_num_variables(dom_id)
do ivar = 1, num_vars
   variable_size = 1
   num_dims = get_num_dims(dom_id, ivar)

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

   state%domain(dom_id)%variable(ivar)%index_start = domain_offset + next_start 
   state%domain(dom_id)%variable(ivar)%index_end   = domain_offset + next_start + variable_size -1
   next_start = next_start + variable_size

enddo

state%model_size = state%model_size + get_domain_size(dom_id)

end subroutine finished_adding_domain

 
!-------------------------------------------------------------------------------
!> Print information in the state structure
!>@todo rename to print_domain_info or something

subroutine state_structure_info(dom_id)

integer, intent(in) :: dom_id ! domain identifier

integer :: ivar, jdim
integer :: num_vars
integer :: num_dims
integer :: array_ids(NF90_MAX_VAR_DIMS)
integer :: array_lengths(NF90_MAX_VAR_DIMS)
character(len=NF90_MAX_NAME) :: dim_name
integer  :: missingINT, spval_int
real(r4) :: missingR4,  spval_r4
real(r8) :: missingR8,  spval_r8
logical  :: has_missing, has_Fill

if ( .not. do_output() ) return

call check_domain_id(dom_id,'state_structure_info')

write(*,*) ' '
write(*,*) 'Reporting on domain # ',dom_id
write(*,*) 'Created by method "', trim(state%domain(dom_id)%method),'"'
if ( state%domain(dom_id)%info_file /= 'NULL' ) &
write(*,*) 'Origin file "', trim(state%domain(dom_id)%info_file), '"'

! summarize all the dimensions used by this domain

num_dims = get_io_num_unique_dims(dom_id)
write(*,'('' Number of dimensions  : '',I2)') num_dims
write(*,'('' unlimdimid            : '',I2)') get_unlimited_dimid(dom_id)
do jdim = 1, num_dims
   write(*,200) jdim, &
             get_io_unique_dim_length(dom_id,jdim), &
        trim(get_io_unique_dim_name(  dom_id,jdim))
enddo
write(*,*)

200 format(4x,i2,', length = ',I8,', name = "',A,'"')
201 format(4x,i2,':         ',2x,'  length = ',I8,', name = "',A,'"')
202 format(4x,i2,':         ',2x,'  ID = ',I8,', length = ',I8,', name = "',A,'"')

! report on each variable in this domain

num_vars = get_num_variables(dom_id)

do ivar = 1, num_vars
   write(*,*)         'VARNAME     : ', trim(get_variable_name(dom_id,ivar))
   write(*,*)         'var_size    : ', get_variable_size(dom_id,ivar)
   write(*,*)         'index_start : ', get_index_start(dom_id,ivar)
   write(*,*)         'index_end   : ', get_index_end(dom_id,ivar)
   write(*,*)         'quantity    : ', get_kind_string(dom_id,ivar)
   write(*,'(A,I3)') ' qty_index   : ', get_kind_index(dom_id,ivar)
   write(*,*)         'clamping    : ', do_io_clamping(dom_id,ivar)
   write(*,*)         'minvalue    : ', get_io_clamping_minval(dom_id,ivar)
   write(*,*)         'maxvalue    : ', get_io_clamping_maxval(dom_id,ivar)
   write(*,*)         'update      : ', do_io_update(dom_id,ivar)
   write(*,*)         'unlimdimid  : ', get_unlimited_dimid(dom_id)

   num_dims = get_num_dims(dom_id,ivar)
   write(*,'('' numdims     : '',I1)') num_dims

   array_lengths(1:num_dims) = get_dim_lengths(dom_id,ivar)
   do jdim = 1, num_dims
       dim_name = get_dim_name(dom_id, ivar, jdim)
       write(*,201) jdim, array_lengths(jdim), trim(dim_name)
   enddo

   ! Report on the native dimensions in the original netCDF file

   num_dims = get_io_num_dims(dom_id,ivar)
   write(*,'('' io_numdims  : '',I1)') num_dims

   array_ids(1:num_dims)     = get_io_dim_ids(dom_id,ivar)
   array_lengths(1:num_dims) = get_io_dim_lengths(dom_id,ivar)
   do jdim = 1, num_dims
       dim_name = get_dim_name(dom_id, ivar, jdim)
       write(*,202) jdim, array_ids(jdim), array_lengths(jdim), trim(dim_name)
   enddo

   if ( state%domain(dom_id)%info_file /= 'NULL' ) then
      write(*,*) 'units             : ', trim(get_units(dom_id,ivar))
      write(*,*) 'short_name        : ', trim(get_short_name(dom_id,ivar))
      write(*,*) 'long_name         : ', trim(get_long_name(dom_id,ivar))
      write(*,*) 'has_missing_value : ', get_has_missing_value(dom_id,ivar)
      write(*,*) 'has_FillValue     : ', get_has_FillValue(dom_id,ivar)

      has_missing = get_has_missing_value(dom_id,ivar)
      has_Fill    = get_has_FillValue(    dom_id,ivar)

      select case (get_xtype(dom_id,ivar))
         case (NF90_INT)
            write(*,*) 'xtype             : ', 'NF90_INT'
            if (has_missing) then
               call get_missing_value(dom_id,ivar,missingINT)
               write(*,*) 'missing_value     : ', missingINT
            endif
            if (has_Fill) then
               call get_FillValue(    dom_id,ivar,spval_int)
               write(*,*) '_FillValue        : ', spval_int
            endif
         case (NF90_FLOAT)
            write(*,*) 'xtype             : ', 'NF90_FLOAT'
            if (has_missing) then
               call get_missing_value(dom_id,ivar,missingR4)
               write(*,*) 'missing_value     : ', missingR4
            endif
            if (has_Fill) then
               call get_FillValue(    dom_id,ivar,spval_r4)
               write(*,*) '_FillValue        : ', spval_r4
            endif
         case (NF90_DOUBLE)
            write(*,*) 'xtype             : ', 'NF90_DOUBLE'
            if (has_missing) then
               call get_missing_value(dom_id,ivar,missingR8)
               write(*,*) 'missing_value     : ', missingR8
            endif
            if (has_Fill)    then
               call get_FillValue(    dom_id,ivar,spval_r8)
               write(*,*) '_FillValue        : ', spval_r8
            endif
      end select

   endif

   !>@todo FIXME : only storing r8 at the moment since DART is not using these values
   !>              to compress and uncompress file information
   if (get_add_offset(dom_id,ivar)   /= MISSING_R8 .and. &
       get_scale_factor(dom_id,ivar) /= MISSING_R8 ) then
      write(*,*) 'add_offset        : ', get_add_offset(dom_id,ivar)
      write(*,*) 'scale_factor      : ', get_scale_factor(dom_id,ivar)
   endif

   write(*,*)
   
enddo

end subroutine state_structure_info


!-------------------------------------------------------------------------------
!> Set DART kinds for domain variables


subroutine set_dart_kinds(dom_id, num_vars, kind_list)

integer, intent(in) :: dom_id
integer, intent(in) :: num_vars
integer, intent(in) :: kind_list(num_vars)

integer :: ivar

do ivar = 1,num_vars
   state%domain(dom_id)%variable(ivar)%dart_kind = kind_list(ivar)
   state%domain(dom_id)%variable(ivar)%kind_string = get_name_for_quantity(kind_list(ivar))
enddo

end subroutine set_dart_kinds


!-------------------------------------------------------------------------------
!> Returns the variable dart kind index

!>@todo need to switch all kinds to qty
function get_kind_index(dom_id, var_id)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: var_id ! variable
integer  :: get_kind_index

get_kind_index = state%domain(dom_id)%variable(var_id)%dart_kind

end function get_kind_index


!-------------------------------------------------------------------------------
!> Returns the variable dart kind string


function get_kind_string(dom_id, var_id)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: var_id ! variable
character(len=obstypelength)  :: get_kind_string

get_kind_string = state%domain(dom_id)%variable(var_id)%kind_string

end function get_kind_string


!-------------------------------------------------------------------------------
!> Return variable id given a domain number and variable name


function get_varid_from_varname(dom_id, varname) result(var_id)

integer,          intent(in) :: dom_id
character(len=*), intent(in) :: varname
integer :: var_id

integer ivar, num_vars

var_id   = -1
num_vars = get_num_variables(dom_id)

do ivar = 1, num_vars
   if ( varname == get_variable_name(dom_id,ivar) ) then
      var_id = ivar
      return
   endif
enddo

end function get_varid_from_varname


!-------------------------------------------------------------------------------
!> Return variable id given a dart kind index


function get_varid_from_kind(dom_id, dart_kind_index) result(var_id)

integer, intent(in)  :: dom_id
integer, intent(in)  :: dart_kind_index
integer :: var_id

integer :: ivar, num_vars

var_id   = -1

if ( get_num_varids_from_kind(dom_id, dart_kind_index) > 1 ) then
   write(string1,*) 'Found ', get_num_varids_from_kind(dom_id, dart_kind_index), &
                    ' > 1'
   write(string2,*) 'for dart kind : ', get_name_for_quantity(dart_kind_index)
   write(string3,*) 'Please use get_varids_from_kind to get a list of indices '
   call error_handler(E_ERR,'get_varid_from_kind', string1, &
              source, text2=string2, text3=string3)
endif

num_vars = get_num_variables(dom_id)

do ivar = 1, num_vars
   if ( dart_kind_index == get_kind_index(dom_id, ivar) ) then
      var_id = ivar
      return
   endif
enddo

end function get_varid_from_kind


!-------------------------------------------------------------------------------
!> Return a list of variable ids in a domain that match a given dart kind index


subroutine get_varids_from_kind(dom_id, dart_kind_index, varid_table)

integer, intent(in)  :: dom_id
integer, intent(in)  :: dart_kind_index
integer, intent(out) :: varid_table(:)

integer :: ivar, indx
integer :: num_vars

if ( size(varid_table) < get_num_varids_from_kind(dom_id, dart_kind_index) ) then
   write(string1,*) 'Found ', get_num_varids_from_kind(dom_id, dart_kind_index)
   write(string2,*) 'varid_table must be at least this size, ', &
                    'you have size(varid_table) = ', size(varid_table)
   call error_handler(E_ERR,'get_varids_from_kind', string1, source, text2=string2)
endif

! Initialize to no variables found
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


subroutine set_update_list(dom_id, num_vars, update_var_list)

integer, intent(in) :: dom_id
integer, intent(in) :: num_vars
logical, intent(in) :: update_var_list(num_vars)

state%domain(dom_id)%variable(:)%io_info%update = update_var_list

end subroutine set_update_list


!-------------------------------------------------------------------------------
!> Returns the number of dimensions for a variable


function get_io_num_dims(dom_id, ivar)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: ivar ! variable
integer :: get_io_num_dims

get_io_num_dims = state%domain(dom_id)%variable(ivar)%io_info%io_numdims

end function get_io_num_dims


!-------------------------------------------------------------------------------
!> Return an array of netcdf dimension ids


function get_io_dim_ids(dom_id, ivar)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: ivar ! variable
integer :: get_io_dim_ids(state%domain(dom_id)%variable(ivar)%io_info%io_numdims)

integer :: num_dims

num_dims = state%domain(dom_id)%variable(ivar)%io_info%io_numdims
get_io_dim_ids = state%domain(dom_id)%variable(ivar)%io_info%io_dimIds(1:num_dims)

end function get_io_dim_ids


!-------------------------------------------------------------------------------
!> Return and array containing the dimension lengths


function get_io_dim_lengths(dom_id, ivar)

integer, intent(in) :: dom_id ! domain
integer, intent(in) :: ivar ! variable
integer :: get_io_dim_lengths(state%domain(dom_id)%variable(ivar)%io_info%io_numdims)

integer :: num_dims

num_dims = get_io_num_dims(dom_id,ivar)
get_io_dim_lengths = state%domain(dom_id)%variable(ivar)%io_info%dimlens(1:num_dims)

end function get_io_dim_lengths


!-------------------------------------------------------------------------------
!> Returns whether a variable should be updated or not


function do_io_update(dom_id, var_id)

integer, intent(in) :: dom_id ! domain identifier
integer, intent(in) :: var_id
logical :: do_io_update

do_io_update = state%domain(dom_id)%variable(var_id)%io_info%update

end function do_io_update


!-------------------------------------------------------------------------------
!> Returns number of variables that match a dart kind index within a domain


function get_num_varids_from_kind(dom_id, dart_kind_index) result(num_found)

integer, intent(in) :: dom_id
integer, intent(in) :: dart_kind_index
integer             :: num_found

integer :: ivar

num_found = 0

do ivar = 1, get_num_variables(dom_id)
   if ( dart_kind_index == get_kind_index(dom_id,ivar) ) then
      num_found = num_found + 1
   endif
enddo

end function get_num_varids_from_kind


!-------------------------------------------------------------------------------
!> Return units of a variable if it exist


function get_units(dom_id, var_id)

integer, intent(in) :: dom_id
integer, intent(in) :: var_id
character(len=NF90_MAX_NAME) :: get_units

get_units = state%domain(dom_id)%variable(var_id)%io_info%units

end function get_units


!-------------------------------------------------------------------------------
!> Return long_name of a variable if it exist


function get_long_name(dom_id, var_id)

integer, intent(in) :: dom_id
integer, intent(in) :: var_id
character(len=NF90_MAX_NAME) :: get_long_name

get_long_name = state%domain(dom_id)%variable(var_id)%io_info%long_name

end function get_long_name


!-------------------------------------------------------------------------------
!> Return short_name of a variable if it exist


function get_short_name(dom_id, var_id)

integer, intent(in) :: dom_id
integer, intent(in) :: var_id
character(len=NF90_MAX_NAME) :: get_short_name

get_short_name = state%domain(dom_id)%variable(var_id)%io_info%short_name

end function get_short_name


!-------------------------------------------------------------------------------
!> Return .true. if the variable has a 'missing_value' 

function get_has_missing_value(dom_id, var_id)

integer, intent(in) :: dom_id
integer, intent(in) :: var_id
logical :: get_has_missing_value

get_has_missing_value = state%domain(dom_id)%variable(var_id)%io_info%has_missing_value

end function get_has_missing_value


!-------------------------------------------------------------------------------
!> Return .true. if the variable has a 'FillValue'

function get_has_FillValue(dom_id, var_id)

integer, intent(in) :: dom_id
integer, intent(in) :: var_id
logical :: get_has_FillValue

get_has_FillValue = state%domain(dom_id)%variable(var_id)%io_info%has_FillValue

end function get_has_FillValue


!-------------------------------------------------------------------------------
!> Return netCDF type of a variable (aka 'xtype')


function get_xtype(dom_id, var_id)

integer, intent(in) :: dom_id
integer, intent(in) :: var_id
integer :: get_xtype

get_xtype = state%domain(dom_id)%variable(var_id)%io_info%xtype

end function get_xtype


!-------------------------------------------------------------------------------
!> Return netCDF missing value of a 64bit real - if the variable has a missing value


subroutine get_missing_value_r8(dom_id, var_id, missing_value_r8)

integer,  intent(in)  :: dom_id
integer,  intent(in)  :: var_id
real(digits12), intent(out) :: missing_value_r8

missing_value_r8 = state%domain(dom_id)%variable(var_id)%io_info%missingR8

end subroutine get_missing_value_r8


!-------------------------------------------------------------------------------
!> Return netCDF missing value of a 32bit real - if the variable has a missing value


subroutine get_missing_value_r4(dom_id, var_id, missing_value_r4)

integer,  intent(in)  :: dom_id
integer,  intent(in)  :: var_id
real(r4), intent(out) :: missing_value_r4

missing_value_r4 = state%domain(dom_id)%variable(var_id)%io_info%missingR4

end subroutine get_missing_value_r4


!-------------------------------------------------------------------------------
!> Return netCDF missing value of a 32bit integer - if the variable has a missing value


subroutine get_missing_value_int(dom_id, var_id, missing_value_int)

integer, intent(in)  :: dom_id
integer, intent(in)  :: var_id
integer, intent(out) :: missing_value_int

missing_value_int = state%domain(dom_id)%variable(var_id)%io_info%missingINT

end subroutine get_missing_value_int


!-------------------------------------------------------------------------------
!> Return spvalR8 value of a variable if it exist


subroutine get_spval_r8(dom_id, var_id, spval_r8)

integer,  intent(in)  :: dom_id
integer,  intent(in)  :: var_id
real(digits12), intent(out) :: spval_r8

spval_r8 = state%domain(dom_id)%variable(var_id)%io_info%spvalR8

end subroutine get_spval_r8


!-------------------------------------------------------------------------------
!> Return spvalR4 value of a variable if it exist


subroutine get_spval_r4(dom_id, var_id, spval_r4)

integer,  intent(in)  :: dom_id
integer,  intent(in)  :: var_id
real(r4), intent(out) :: spval_r4

spval_r4 = state%domain(dom_id)%variable(var_id)%io_info%spvalR4

end subroutine get_spval_r4


!-------------------------------------------------------------------------------
!> Return spvalINT value of a variable if it exist


subroutine get_spval_int(dom_id, var_id, spval_int)

integer, intent(in)  :: dom_id
integer, intent(in)  :: var_id
integer, intent(out) :: spval_int

spval_int = state%domain(dom_id)%variable(var_id)%io_info%spvalINT

end subroutine get_spval_int


!-------------------------------------------------------------------------------
!> Return offset of a variable if it exists


function get_add_offset(dom_id, var_id)

integer, intent(in) :: dom_id
integer, intent(in) :: var_id
real(r8) :: get_add_offset

get_add_offset = state%domain(dom_id)%variable(var_id)%io_info%add_offset

end function get_add_offset


!-------------------------------------------------------------------------------
!> Return scale factor of a variable if it exists


function get_scale_factor(dom_id, var_id)

integer, intent(in) :: dom_id
integer, intent(in) :: var_id
real(r8) :: get_scale_factor

get_scale_factor = state%domain(dom_id)%variable(var_id)%io_info%scale_factor

end function get_scale_factor


!-------------------------------------------------------------------------------
!> Assert that adding a domain will not cause the maximum number of domains
!> to be exceeded.


subroutine assert_below_max_num_domains(context)
character(len=*), optional, intent(in) :: context

if (present(context)) then
   write(string1,*)trim(context), ':requesting to add domain #', &
                   state%num_domains + 1
else
   write(string1,*)'requesting to add domain #', state%num_domains + 1
endif

if (state%num_domains + 1 > MAX_NUM_DOMS) then
   write(string2,*)'maximum number of domains is ',MAX_NUM_DOMS
   write(string3,*)'increase "MAX_NUM_DOMS" in the common/types_mod.f90 and recompile'
   call error_handler(E_ERR, 'assert_below_max_num_domains', string1, &
              source, text2=string2, text3=string3)
endif

end subroutine assert_below_max_num_domains


!-------------------------------------------------------------------------------
!> Transform the state type to a domain that can be used to write diagnostic 
!> files.
!> One domain for the whole state. If there is more than one domain:
!>   Need to change variable names
!>   Need to change dimension names


function create_diagnostic_structure() result (diag_id)

integer :: diag_id

integer :: i, j, k, var ! loop variables
character(len=8)   :: dom_str = ''

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
!>


subroutine check_domain_id(dom_id, message)
integer,          intent(in) :: dom_id
character(len=*), intent(in) :: message

if (dom_id == diagnostic_domain) return

if (dom_id > state%num_domains .or. dom_id < 0 ) then
   write(string1,*)'number of known domains is ',state%num_domains
   write(string2,*)'requesting information for unknown domain ',dom_id
   call error_handler(E_ERR,'check_domain_id', message, &
              source, text2=string1, text3=string2)
endif

end subroutine check_domain_id

!> @}
end module state_structure_mod

