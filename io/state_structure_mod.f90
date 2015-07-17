! Aim: to have a progvar type that can be used for IO.
! How do you do domains?
! Different variables in different domains.  I think we could treat the components 
! of CESM as different domains.
! Note:
! (5) the format for 'wrf_state_variables' is an array of 5 strings:
!     wrf netcdf variable name, dart KIND_xxx string, type string (must be
!     unique, will soon be obsolete, we hope), 'UPDATE', and '999' if the
!     array is part of all domains.  otherwise, it is a string with the domain
!     numbers (e.g. '12' for domains 1 and 2, '13' for domains 1 and 3).

! Can we control the order of the state vector?

!-------------------------------------------------------------------------------

module state_structure_mod

use utilities_mod,        only : E_ERR, error_handler
use types_mod,            only : r8, i8
use utilities_mod,        only : nc_check
use sort_mod,             only : index_sort
use netcdf

implicit none

private
public :: static_init_state_type, &
          add_domain,             &
          get_domain_size,        &
          get_num_variables,      &
          get_variable_size,      &
          get_num_dims,           &
          get_dim_ids,            &
          get_variable_name,      &
          get_dim_lengths,        &
          get_num_domains,        &
          get_dim_name,           &
          get_dim_length,         &
          add_time_unlimited,     &
          set_var_id,             &
          get_unique_dim_name,    &
          get_num_unique_dims,    &
          get_unique_dim_length,  &
          get_sum_variables,      &
          get_sum_variables_below

!-------------------------------------------------------------------------------
! global variables
integer, parameter :: max_num_domains = 10
logical, save :: state_initialized = .false.


!-------------------------------------------------------------------------------
! state_type -> domain_type -> variable_type
! #domains      #variables     info
!-------------------------------------------------------------------------------
! Describes a variable, enough for netcdf io purposes to read and write
type netcdf_var_type
   private
   ! netcdf variable id 
   integer :: netcdf_id ! HK Do you want one for reading, one for writing?
   integer :: xtype     ! netCDF variable type (NF90_double, etc.)
    
   ! variable information
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   integer :: var_size
   integer :: index1 ! location in dart state vector of first occurrence
   integer :: indexN ! location in dart state vector of last  occurrence
   
   ! dimension information
   integer :: numdims ! number of dims - excluding TIME?
   character(len=NF90_MAX_NAME), dimension(NF90_MAX_VAR_DIMS) :: dimname
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimIds
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens

   ! FUTURE VARIABLES
   ! ! stride information
   ! integer, dimension(NF90_MAX_VAR_DIMS) :: start
   ! integer, dimension(NF90_MAX_VAR_DIMS) :: count
   ! integer, dimension(NF90_MAX_VAR_DIMS) :: stride
   ! integer, dimension(NF90_MAX_VAR_DIMS) :: scale_factor = 1
   ! ! can you check for attributes that aren't there?!
   ! integer, dimension(NF90_MAX_VAR_DIMS) :: add_offset
   !
   ! ! dart information
   ! integer :: dart_kind
   ! character(len=paramname_length) :: kind_string
   !
   ! ! restrict the range of information
   ! logical  :: clamping     ! does variable need to be range-restricted before
   ! real(r8) :: range(2)     ! being stuffed back into MPAS analysis file.
   ! logical  :: out_of_range_fail  ! is out of range fatal if range-checking?

end type netcdf_var_type

!-------------------------------------------------------------------------------
! Describes a domain
type domain_type
   private
   ! variables in domain
   integer :: num_variables
   type(netcdf_var_type), allocatable :: variable(:) 
   
   ! dimension informaion for domain
   integer :: num_unique_dims
   character(len=NF90_MAX_NAME), allocatable :: unique_dim_names(:)
   integer, allocatable :: unique_dim_length(:)
   integer :: unlimDimId

   ! number of elements in the domain
   integer(i8) :: size

   ! netcdf file describing the shape of the variable
   character(len=256) :: info_file  ! do we need this?
   
end type domain_type

!-------------------------------------------------------------------------------
! Describes a state vector
type state_type
   private
   ! domains or separate files
   integer           :: num_domains
   type(domain_type) :: domain(max_num_domains) ! number of domains

   ! size of the state vector
   integer(i8)       :: model_size

   !>@todo is this different to character(len=128) :: name?  Redundant?
   !character(len=1024)  :: domain_name(max_num_domains)
end type state_type

!-------------------------------------------------------------------------------
! module handle
type(state_type) :: state

interface add_domain
   module procedure add_domain_blank
   module procedure add_domain_from_file
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
! Accessor function for setting info
! The model is going to call these
!-------------------------------------------------------------------------------
!> Pass an indentifier to the model for this domain
!> Set up the domain type
!> Does this need to be a function or a subroutine
!-------------------------------------------------------------------------------
function add_domain_from_file(info_file, num_variables, var_names) result(dom_identifier)

character(len=*), intent(in) :: info_file
integer,          intent(in) :: num_variables
character(len=*), intent(in) :: var_names(num_variables)
integer           :: dom_identifier

integer :: ivar

if (.not. state_initialized) call static_init_state_type()

! add to domains
state%num_domains = state%num_domains + 1
dom_identifier    = state%num_domains !>@todo this should be a handle.

state%domain(dom_identifier)%info_file = info_file

! set number of variables in this domain
state%domain(dom_identifier)%num_variables = num_variables

! load up the variable names
allocate(state%domain(dom_identifier)%variable(num_variables))

do ivar = 1, num_variables
   state%domain(dom_identifier)%variable(ivar)%varname = var_names(ivar)
enddo

! load up variable id's and sizes
call load_state_variable_info(state%domain(dom_identifier))

! load up the domain unique dimension info
call load_unique_dim_info(dom_identifier)

end function add_domain_from_file

!-------------------------------------------------------------------------------
!> Add a blank domain - one variable called state, length = model_size
!-------------------------------------------------------------------------------
function add_domain_blank(model_size) result(dom_identifier)

integer(i8), intent(in) :: model_size
integer                 :: dom_identifier

if (.not. state_initialized) call static_init_state_type()

! add to domains
state%num_domains = state%num_domains + 1
dom_identifier = state%num_domains

! domain
state%domain(dom_identifier)%num_variables = 1
state%domain(dom_identifier)%size = model_size

! variable
allocate(state%domain(dom_identifier)%variable(1))
state%domain(dom_identifier)%variable(1)%varname   = 'state'
state%domain(dom_identifier)%variable(1)%long_name = 'state'
state%domain(dom_identifier)%variable(1)%units     = ''
state%domain(dom_identifier)%variable(1)%numdims   = 1
state%domain(dom_identifier)%variable(1)%var_size  = model_size

! dimension
state%domain(dom_identifier)%variable(1)%dimname(1) = 'model_size'
state%domain(dom_identifier)%variable(1)%dimlens(1) = model_size

! load up the domain unique dimension info
state%domain(dom_identifier)%num_unique_dims = 1
allocate(state%domain(dom_identifier)%unique_dim_names(1))
allocate(state%domain(dom_identifier)%unique_dim_length(1))
state%domain(dom_identifier)%unique_dim_names(1)  = 'model_size'
state%domain(dom_identifier)%unique_dim_length(1) =  model_size

end function add_domain_blank

!-------------------------------------------------------------------------------
!> Load metadata from netcdf file info state_strucutre
!-------------------------------------------------------------------------------
subroutine load_state_variable_info(domain)

type(domain_type), intent(in) :: domain

integer :: ret ! netcdf return value
character(len=512) :: netcdf_filename
integer :: ncfile ! netdcf file id - should this be part of the domain handle?

! open netcdf file - all restart files in a domain have the same info?
ret = nf90_open(domain%info_file, NF90_NOWRITE, ncfile)
call nc_check(ret, 'load_state_variable_info nf90_open', trim(domain%info_file))

! get variable ids
call load_variable_ids(ncfile, domain)

! get all variable sizes, only readers store dimensions?
call load_variable_sizes(ncfile, domain)

! close netcdf file
ret = nf90_close(ncfile)
call nc_check(ret, 'load_state_variable_info nf90_close', trim(domain%info_file))

end subroutine load_state_variable_info

!-------------------------------------------------------------------------------
!> Load netcdf variable ids
!-------------------------------------------------------------------------------
subroutine load_variable_ids(ncfile, domain)

integer,           intent(in)    :: ncfile ! netdcf file id - should this be part of the domain handle?
type(domain_type) :: domain

integer :: ret ! netcdf retrun value
integer :: ivar ! loop variable

do ivar = 1, domain%num_variables

   ! load netcdf id from variable name
   ret = nf90_inq_varid(ncfile, domain%variable(ivar)%varname, &
                                domain%variable(ivar)%netcdf_id)

   call nc_check(ret, 'load_variable_ids, nf90_inq_var_id', &
                 trim(domain%variable(ivar)%varname))

enddo

end subroutine load_variable_ids

!-------------------------------------------------------------------------------
!> load dimension information and calculate variable and domain sizes
!-------------------------------------------------------------------------------
subroutine load_variable_sizes(ncfile, domain)

integer :: ncfile ! netdcf file id - should this be part of the domain handle?
type(domain_type) :: domain

integer :: ivar, jdim !< loop variables
integer :: var_size, dom_size !< temporary sum variable
integer :: ret ! netcdf retrun value

dom_size = 0

do ivar = 1, domain%num_variables

   ! from netcdf id load variable dimension and ids
   ret = nf90_inquire_variable(ncfile,  domain%variable(ivar)%netcdf_id, &
                              ndims  = domain%variable(ivar)%numdims,    &
                              dimids = domain%variable(ivar)%dimIds,     &
                              xtype  = domain%variable(ivar)%xtype)

   call nc_check(ret, 'load_variable_sizes, inq_variable', &
                 trim(domain%variable(ivar)%varname))

   var_size = 1
   do jdim = 1, domain%variable(ivar)%numdims

      ! load dimension names and lengths
      ret = nf90_inquire_dimension(ncfile, domain%variable(ivar)%dimIds(jdim),  &
                                      name=domain%variable(ivar)%dimname(jdim), &
                                      len= domain%variable(ivar)%dimlens(jdim))

      call nc_check(ret, 'load_variable_sizes, inq_dimension', &
                    trim(domain%variable(ivar)%dimname(jdim)))
      
      var_size = var_size * domain%variable(ivar)%dimlens(jdim)

   enddo

   domain%variable(ivar)%var_size = var_size
   dom_size = dom_size + var_size

enddo

domain%size = dom_size

end subroutine load_variable_sizes

!-------------------------------------------------------------------------------
!> Identify the unique dimensions within a domain
!> I think you can sort dimension ids and remove duplicates
!> What if there are 0 dimensions?
!-------------------------------------------------------------------------------
subroutine load_unique_dim_info(idom)

integer :: idom ! domain index

integer, allocatable :: array_of_dimids(:)
integer, allocatable :: index(:)
character(len=NF90_MAX_NAME), allocatable :: array_of_names(:)
integer, allocatable :: array_of_lengths(:)
logical, allocatable :: u(:)

integer :: count, ivar, jdim
integer :: ndims
integer :: num_unique

num_unique = 0

ndims = get_domain_num_dims(idom)

allocate(array_of_dimids(ndims))
allocate(array_of_names(ndims))
allocate(array_of_lengths(ndims))
allocate(index(ndims))
allocate(u(ndims))

count = 1

do ivar = 1, state%domain(idom)%num_variables
   do jdim = 1, state%domain(idom)%variable(ivar)%numdims
      array_of_dimids(count)  = state%domain(idom)%variable(ivar)%dimIds(jdim)
      array_of_names(count)   = state%domain(idom)%variable(ivar)%dimname(jdim)
      array_of_lengths(count) = state%domain(idom)%variable(ivar)%dimlens(jdim)
      count = count + 1
   enddo
enddo

call index_sort(array_of_dimids, index, ndims)

count = 1
u(:) = .false.
u(1) = .true.

do ivar = 2, ndims
   if(array_of_dimids(index(ivar)) /= array_of_dimids(index(ivar-1))) then
      count = count + 1
      u(ivar) = .true.
   endif
enddo

state%domain(idom)%num_unique_dims = count

allocate(state%domain(idom)%unique_dim_names(count))
allocate(state%domain(idom)%unique_dim_length(count))

count = 1

do ivar = 1, ndims
   if(u(ivar)) then
      state%domain(idom)%unique_dim_names(count)  = array_of_names(index(ivar))
      state%domain(idom)%unique_dim_length(count) = array_of_lengths(index(ivar))
      count = count + 1
   endif
enddo

deallocate(array_of_dimids, array_of_names, array_of_lengths, index, u)

end subroutine load_unique_dim_info

!-------------------------------------------------------------------------------
! dump out the information for a variable
!-------------------------------------------------------------------------------
! Accessor functions for read/write functions
! This is state_vector_io_mod
!-------------------------------------------------------------------------------
function get_num_domains()

integer :: get_num_domains

get_num_domains = state%num_domains

end function get_num_domains

!-------------------------------------------------------------------------------
! Number of elements in the domain
!-------------------------------------------------------------------------------
function get_domain_size(idom)

integer, intent(in) :: idom
integer :: get_domain_size

get_domain_size = state%domain(idom)%size

end function get_domain_size

!-------------------------------------------------------------------------------
! Number of variables in the domain
!-------------------------------------------------------------------------------
function get_num_variables(idom)

integer, intent(in) :: idom
integer :: get_num_variables

get_num_variables = state%domain(idom)%num_variables

end function get_num_variables

!-------------------------------------------------------------------------------
! Size of specific variable
!-------------------------------------------------------------------------------
function get_variable_size(idom, ivar)

integer, intent(in) :: idom
integer, intent(in) :: ivar ! should you access this by index?
integer :: get_variable_size

get_variable_size = state%domain(idom)%variable(ivar)%var_size

end function get_variable_size

!-------------------------------------------------------------------------------
! Total number of dimension for variable
!-------------------------------------------------------------------------------
function get_num_dims(idom, ivar)

integer, intent(in) :: idom
integer, intent(in) :: ivar ! should you access this by index?
integer :: get_num_dims

get_num_dims = state%domain(idom)%variable(ivar)%numdims

end function get_num_dims

!-------------------------------------------------------------------------------
! Return an array of netcdf dimension ids
!-------------------------------------------------------------------------------
function get_dim_ids(idom, ivar)

integer, intent(in) :: idom
integer, intent(in) :: ivar ! should you access this by index?
integer :: get_dim_ids(state%domain(idom)%variable(ivar)%numdims)

integer :: num_dims

num_dims = state%domain(idom)%variable(ivar)%numdims

get_dim_ids = state%domain(idom)%variable(ivar)%dimIds(1:num_dims)

end function get_dim_ids

!-------------------------------------------------------------------------------
! Return and array containing the dimension lengths
!-------------------------------------------------------------------------------
function get_dim_lengths(idom, ivar)

integer, intent(in) :: idom
integer, intent(in) :: ivar
integer :: get_dim_lengths(state%domain(idom)%variable(ivar)%numdims)

integer :: num_dims

num_dims = state%domain(idom)%variable(ivar)%numdims

get_dim_lengths = state%domain(idom)%variable(ivar)%dimlens(1:num_dims)

end function get_dim_lengths

!-------------------------------------------------------------------------------
! Return variable name
!-------------------------------------------------------------------------------
function get_variable_name(idom, ivar)

integer, intent(in) :: idom
integer, intent(in) :: ivar ! should you access this by index?
character(len=NF90_MAX_NAME)  :: get_variable_name

get_variable_name = state%domain(idom)%variable(ivar)%varname

end function get_variable_name

!-------------------------------------------------------------------------------
! Should you even use this?
!-------------------------------------------------------------------------------
subroutine set_var_id(idom, ivar, new_varid)

integer, intent(in) :: idom ! domain
integer, intent(in) :: ivar ! variable
integer, intent(in) :: new_varid

state%domain(idom)%variable(ivar)%netcdf_id = new_varid

end subroutine set_var_id

!-------------------------------------------------------------------------------
! For dimensions
!-------------------------------------------------------------------------------
function get_dim_name(idom, ivar, jdim)

integer, intent(in) :: idom ! domain
integer, intent(in) :: ivar ! variable
integer, intent(in) :: jdim ! dimension
character(len=NF90_MAX_NAME) :: get_dim_name

get_dim_name = state%domain(idom)%variable(ivar)%dimname(jdim)

end function get_dim_name

!-------------------------------------------------------------------------------
! Return dimension length
!-------------------------------------------------------------------------------
function get_dim_length(idom, ivar, jdim)

integer, intent(in) :: idom ! domain
integer, intent(in) :: ivar ! variable
integer, intent(in) :: jdim ! dimension
integer :: get_dim_length

get_dim_length = state%domain(idom)%variable(ivar)%dimlens(jdim)

end function get_dim_length

!-------------------------------------------------------------------------------
!> Return the number of dimensions in a domain
!> Repeat dimensions are allowed
!-------------------------------------------------------------------------------
function get_domain_num_dims(idom)

integer, intent(in) :: idom ! domain index
integer :: get_domain_num_dims
integer :: ivar

get_domain_num_dims = 0

do ivar = 1, state%domain(idom)%num_variables
  get_domain_num_dims =  get_domain_num_dims + &
                         state%domain(idom)%variable(ivar)%numdims
enddo

end function get_domain_num_dims

!-------------------------------------------------------------------------------
! Return the number of unique dimensions
!-------------------------------------------------------------------------------
function get_num_unique_dims(idom)

integer, intent(in) :: idom
integer :: get_num_unique_dims

get_num_unique_dims = state%domain(idom)%num_unique_dims

end function get_num_unique_dims

!-------------------------------------------------------------------------------
! Return the unique dimension names
!-------------------------------------------------------------------------------
function get_unique_dim_name(idom, ivar)

integer, intent(in) :: idom ! domain identifier
integer, intent(in) :: ivar ! index into array, not connected to dimId
character(len=NF90_MAX_NAME) :: get_unique_dim_name

get_unique_dim_name = state%domain(idom)%unique_dim_names(ivar)

end function get_unique_dim_name

!-------------------------------------------------------------------------------
! Return the unique dimension lengths
!-------------------------------------------------------------------------------
function get_unique_dim_length(idom, ivar)

integer, intent(in) :: idom ! domain identifier
integer, intent(in) :: ivar ! index into array, not connected to dimId
integer :: get_unique_dim_length

get_unique_dim_length = state%domain(idom)%unique_dim_length(ivar)

end function get_unique_dim_length

!-------------------------------------------------------------------------------
!> Adding space for an unlimited dimension in the dimesion arrays
!> The unlimited dimension needs to be last in the list for def_var
!-------------------------------------------------------------------------------
subroutine add_time_unlimited(unlimited_dimId)

integer, intent(in)  :: unlimited_dimId

integer :: ivar !> loop variable

call error_handler(E_ERR, 'add_time_unlimited', 'does not exist')

end subroutine add_time_unlimited

!-------------------------------------------------------------------------------
! number of variables below start_var.
!-------------------------------------------------------------------------------
function get_sum_variables_below(start_var, domain)

integer, intent(in) :: start_var
integer, intent(in) :: domain
integer             :: get_sum_variables_below

integer :: ivar

get_sum_variables_below = 0

do ivar = 1, domain -1
   get_sum_variables_below = get_sum_variables_below + state%domain(ivar)%size
enddo

do ivar = 1, start_var -1
   get_sum_variables_below = get_sum_variables_below + &
                             state%domain(domain)%variable(ivar)%var_size
enddo

end function get_sum_variables_below

!-------------------------------------------------------------------------------
! number of variables in range of start_var and end_var
!-------------------------------------------------------------------------------
function get_sum_variables(start_var, end_var, domain)

integer, intent(in) :: start_var, end_var
integer, intent(in) :: domain
integer             :: get_sum_variables

integer :: ivar

get_sum_variables = 0

do ivar = start_var, end_var
   get_sum_variables = get_sum_variables + &
                       state%domain(domain)%variable(ivar)%var_size
enddo

end function get_sum_variables

!-------------------------------------------------------------------------------

end module state_structure_mod
