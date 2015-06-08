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

!------------------------------------

module state_structure_mod

use utilities_mod,        only : E_ERR, error_handler
use types_mod,            only : r8, i8
use utilities_mod,        only : nc_check
use sort_mod,             only : index_sort
use netcdf

implicit none

private
public :: static_init_state_type
public :: add_domain
public :: get_domain_size, get_num_variables, get_variable_size, get_num_dims, get_dim_ids, &
          get_variable_name, get_dim_lengths, get_num_domains, &
          get_dim_name, get_dim_length, add_time_unlimited, set_var_id, &
          get_unique_dim_name, get_num_unique_dims, get_unique_dim_length, &
          get_sum_variables, get_sum_variables_below


integer, parameter :: max_num_domains = 10
logical, save :: state_initialized = .false.


!------------------------------------
! state_type -> domain_type -> variable_type
! #domains      #variables     info

! Describes a variable, enough for netcdf io purposes to read and write
type netcdf_var_type
   private ! do you want these to be private? With accessor functions to set?
   integer :: netcdf_id ! Netcdf variable id HK Do you want one for reading, one for writing?
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   integer :: numdims       ! number of dims - excluding TIME?
   character(len=NF90_MAX_NAME), dimension(NF90_MAX_VAR_DIMS) :: dimname
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimIds
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
   integer, dimension(NF90_MAX_VAR_DIMS) :: start
   integer, dimension(NF90_MAX_VAR_DIMS) :: count
   integer, dimension(NF90_MAX_VAR_DIMS) :: stride
   integer, dimension(NF90_MAX_VAR_DIMS) :: scale_factor = 1
   integer, dimension(NF90_MAX_VAR_DIMS) :: add_offset ! can you check for attributes that aren't there?
   integer :: xtype         ! netCDF variable type (NF90_double, etc.)
   integer :: var_size
   ! do you need the start and end?
   integer :: index1        ! location in dart state vector of first occurrence
   integer :: indexN        ! location in dart state vector of last  occurrence
   integer :: dart_kind
   !character(len=paramname_length) :: kind_string
   logical  :: clamping     ! does variable need to be range-restricted before
   real(r8) :: range(2)     ! being stuffed back into MPAS analysis file.
   logical  :: out_of_range_fail  ! is out of range fatal if range-checking?
end type netcdf_var_type

! describes a domain
type domain_type
   private
   integer :: num_variables
   character(len=128) :: name
   character(len=256) :: info_file  ! do we need this?
   type(netcdf_var_type), allocatable :: variable(:) ! number of variables in each domain
   character(len=NF90_MAX_NAME), allocatable :: unique_dim_names(:) ! list of unique dimension ids in this domain
   integer, allocatable :: unique_dim_length(:) ! lengths of unique dimension ids in this domain - same order as unique_dim_names
   integer :: num_unique_dims ! number of unique dimensions in this domain
   integer(i8) :: size  ! number of elements in the domain
end type domain_type

! Describes a state vector
type state_type
   private
   integer :: num_domains ! domains or separate files.
   type(domain_type)    :: domain(max_num_domains) ! number of domains
   !character(len=1024)  :: domain_name(max_num_domains) !>@todo is this different to character(len=128) :: name?  Redundant?
end type state_type

type(state_type) :: state

interface add_domain
   module procedure add_domain_blank
   module procedure add_domain_from_info
end interface

contains

!------------------------------------
!> initalizes state_type
subroutine static_init_state_type()

if(state_initialized) return !?

state%num_domains = 0
state_initialized = .true.

end subroutine static_init_state_type

!------------------------------------
! Accessor function for setting info
! The model is going to call these
!------------------------------------

!> Pass an indentifier to the model for this domain
!> set up the domain type
!> Does this need to be a function or a subroutine
function add_domain_from_info(info_file, num_variables, var_names) result(dom_identifier)

character(len=*),             intent(in) :: info_file
integer,                      intent(in) :: num_variables
character(len=*),             intent(in) :: var_names(num_variables)
integer                                  :: dom_identifier

integer :: i

if (.not. state_initialized) call static_init_state_type()

! add to domains
state%num_domains = state%num_domains + 1
dom_identifier = state%num_domains

state%domain(dom_identifier)%info_file = info_file

! set number of variables in this domain
state%domain(dom_identifier)%num_variables = num_variables

! load up the variable names
allocate(state%domain(dom_identifier)%variable(num_variables))

do i = 1, num_variables
   state%domain(dom_identifier)%variable(i)%varname = var_names(i)
enddo

call get_state_variable_info(state%domain(dom_identifier))

! load up the domain unique dimension info
call get_unique_dimensions(dom_identifier)

end function add_domain_from_info

!------------------------------------
!> Add a blank domain - one variable called state, length = model_size
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
state%domain(dom_identifier)%variable(1)%varname = 'state'
state%domain(dom_identifier)%variable(1)%long_name = 'state'
state%domain(dom_identifier)%variable(1)%units = ''
state%domain(dom_identifier)%variable(1)%numdims = 1
state%domain(dom_identifier)%variable(1)%var_size = model_size

! dimension
state%domain(dom_identifier)%variable(1)%dimname(1) = 'model_size'
state%domain(dom_identifier)%variable(1)%dimlens(1) = model_size

! load up the domain unique dimension info
state%domain(dom_identifier)%num_unique_dims = 1
allocate(state%domain(dom_identifier)%unique_dim_names(1))
allocate(state%domain(dom_identifier)%unique_dim_length(1))
state%domain(dom_identifier)%unique_dim_names(1) = 'model_size'
state%domain(dom_identifier)%unique_dim_length(1) = model_size

end function add_domain_blank

!------------------------------------
!>
subroutine get_state_variable_info(domain)

type(domain_type) :: domain

integer :: ret ! netcdf return value
character(len=512) :: netcdf_filename
integer :: ncfile ! netdcf file id - should this be part of the domain handle?

! open netcdf file - all restart files in a domain have the same info?
ret = nf90_open(domain%info_file, NF90_NOWRITE, ncfile)
call nc_check(ret, 'get_state_variable_info opening', trim(domain%info_file))

! get variable ids
call get_variable_ids(ncfile, domain)

! get all variable sizes, only readers store dimensions?
call get_variable_sizes(ncfile, domain)

! close netcdf file
ret = nf90_close(ncfile)
call nc_check(ret, 'get_state_variable_info closing', trim(domain%info_file))


end subroutine get_state_variable_info

!------------------------------------
!> Get netcdf variable ids
subroutine get_variable_ids(ncfile, domain)

integer :: ncfile ! netdcf file id - should this be part of the domain handle?
type(domain_type) :: domain

integer :: ret ! netcdf retrun value
integer :: i ! loop variable

do i = 1, domain%num_variables

   ret = nf90_inq_varid(ncfile, domain%variable(i)%varname, domain%variable(i)%netcdf_id)
   call nc_check(ret, 'get_variable_ids, inq_var_id', trim(domain%variable(i)%varname))

enddo

end subroutine get_variable_ids

!------------------------------------
!> 
subroutine get_variable_sizes(ncfile, domain)

integer :: ncfile ! netdcf file id - should this be part of the domain handle?
type(domain_type) :: domain

integer :: i, j !< loop variable
integer :: sum !< temporary sum variable
integer :: ret ! netcdf retrun value

sum = 0

do i = 1, domain%num_variables

   ret = nf90_inquire_variable(ncfile, domain%variable(i)%netcdf_id, ndims=domain%variable(i)%numdims, dimids=domain%variable(i)%dimIds, xtype=domain%variable(i)%xtype)
   call nc_check(ret, 'totalsize, inq_var', trim(domain%variable(i)%varname))

   do j = 1, domain%variable(i)%numdims
      ret = nf90_inquire_dimension(ncfile, domain%variable(i)%dimIds(j), name=domain%variable(i)%dimname(j), len=domain%variable(i)%dimlens(j))
      call nc_check(ret, 'totalsize, inq_dimlen', trim(domain%variable(i)%dimname(j)))
   enddo

   domain%variable(i)%var_size = product(domain%variable(i)%dimlens(1:domain%variable(i)%numdims))

   sum = sum + domain%variable(i)%var_size

enddo

!domain_size = sum(variable_sizes(:, domain))
domain%size = sum

end subroutine get_variable_sizes

!------------------------------------
! dump out the information for a variable

!------------------------------------
! Accessor functions for read/write functions
! This is state_vector_io_mod
!
!------------------------------------
function get_num_domains()

integer :: get_num_domains

get_num_domains = state%num_domains

end function get_num_domains

!------------------------------------
! For domains
!------------------------------------
function get_domain_size(dom)

integer, intent(in) :: dom
integer :: get_domain_size

get_domain_size = state%domain(dom)%size

end function get_domain_size

!------------------------------------
function get_num_variables(dom)

integer, intent(in) :: dom
integer :: get_num_variables

get_num_variables = state%domain(dom)%num_variables

end function get_num_variables

!------------------------------------
! For variables
!------------------------------------
function get_variable_size(dom, i)

integer, intent(in) :: dom
integer, intent(in) :: i ! should you access this by index?
integer :: get_variable_size

get_variable_size = state%domain(dom)%variable(i)%var_size

end function get_variable_size

!------------------------------------
!> Returns the number of dimensions a variable has
function get_num_dims(dom, i)

integer, intent(in) :: dom
integer, intent(in) :: i ! should you access this by index?
integer :: get_num_dims

get_num_dims = state%domain(dom)%variable(i)%numdims

end function get_num_dims

!------------------------------------
function get_dim_ids(dom, i)

integer, intent(in) :: dom
integer, intent(in) :: i ! should you access this by index?
integer :: get_dim_ids(state%domain(dom)%variable(i)%numdims)

get_dim_ids = state%domain(dom)%variable(i)%dimIds(1:state%domain(dom)%variable(i)%numdims)

end function get_dim_ids

!------------------------------------
function get_dim_lengths(dom, i)

integer, intent(in) :: dom
integer, intent(in) :: i
integer :: get_dim_lengths(state%domain(dom)%variable(i)%numdims)

get_dim_lengths = state%domain(dom)%variable(i)%dimlens(1:state%domain(dom)%variable(i)%numdims)

end function get_dim_lengths

!------------------------------------
function get_variable_name(dom, i)

integer, intent(in) :: dom
integer, intent(in) :: i ! should you access this by index?
character(len=NF90_MAX_NAME)  :: get_variable_name

get_variable_name = state%domain(dom)%variable(i)%varname

end function get_variable_name

!------------------------------------
! Should you even use this?
subroutine set_var_id(dom, i, new_varid)

integer, intent(in) :: dom ! domain
integer, intent(in) :: i ! variable
integer, intent(in) :: new_varid

state%domain(dom)%variable(i)%netcdf_id = new_varid

end subroutine set_var_id

!------------------------------------
! For dimensions
!------------------------------------
function get_dim_name(dom, var, dim)

integer, intent(in) :: dom ! domain
integer, intent(in) :: var ! variable
integer, intent(in) :: dim ! dimension
character(len=NF90_MAX_NAME) :: get_dim_name

get_dim_name = state%domain(dom)%variable(var)%dimname(dim)

end function get_dim_name

!------------------------------------
function get_dim_length(dom, var, dim)

integer, intent(in) :: dom ! domain
integer, intent(in) :: var ! variable
integer, intent(in) :: dim ! dimension
integer :: get_dim_length

get_dim_length = state%domain(dom)%variable(var)%dimlens(dim)

end function get_dim_length

!------------------------------------
!> Return the number of dimensions in a domain
!> Repeat dimensions are allowed
function get_domain_num_dims(dom)

integer, intent(in) :: dom ! domain index
integer :: get_domain_num_dims
integer :: i

get_domain_num_dims = 0

do i = 1, state%domain(dom)%num_variables
  get_domain_num_dims =  get_domain_num_dims + state%domain(dom)%variable(i)%numdims
enddo

end function get_domain_num_dims

!------------------------------------
!> Identify the unique dimensions within a domain
!> I think you can sort dimension ids and remove duplicates
!> What if there are 0 dimensions?
subroutine get_unique_dimensions(dom)

integer :: dom ! domain index

integer, allocatable :: array_of_dimids(:)
integer, allocatable :: index(:)
character(len=NF90_MAX_NAME), allocatable :: array_of_names(:)
integer, allocatable :: array_of_lengths(:)
logical, allocatable :: u(:)

integer :: count, i, j
integer :: ndims
integer :: num_unique

num_unique = 0

ndims = get_domain_num_dims(dom)

allocate(array_of_dimids(ndims))
allocate(array_of_names(ndims))
allocate(array_of_lengths(ndims))
allocate(index(ndims))
allocate(u(ndims))

count = 1

do i = 1, state%domain(dom)%num_variables
   do j = 1, state%domain(dom)%variable(i)%numdims
      array_of_dimids(count) = state%domain(dom)%variable(i)%dimIds(j)
      array_of_names(count) = state%domain(dom)%variable(i)%dimname(j)
      array_of_lengths(count) = state%domain(dom)%variable(i)%dimlens(j)
      count = count + 1
   enddo
enddo

call index_sort(array_of_dimids, index, ndims)

count = 1
u(:) = .false.
u(1) = .true.

do i = 2, ndims
   if(array_of_dimids(index(i)) /= array_of_dimids(index(i-1))) then
      count = count + 1
      u(i) = .true.
   endif
enddo

state%domain(dom)%num_unique_dims = count

allocate(state%domain(dom)%unique_dim_names(state%domain(dom)%num_unique_dims))
allocate(state%domain(dom)%unique_dim_length(state%domain(dom)%num_unique_dims))

count = 1

do i = 1, ndims
   if(u(i)) then
      state%domain(dom)%unique_dim_names(count)  = array_of_names(index(i))
      state%domain(dom)%unique_dim_length(count) = array_of_lengths(index(i))
      count = count + 1
   endif
enddo

deallocate(array_of_dimids, array_of_names, array_of_lengths, index, u)

end subroutine get_unique_dimensions

!-------------------------------------------------------
function get_num_unique_dims(dom)

integer, intent(in) :: dom
integer :: get_num_unique_dims

get_num_unique_dims = state%domain(dom)%num_unique_dims

end function get_num_unique_dims

!-------------------------------------------------------
function get_unique_dim_name(dom, i)

integer, intent(in) :: dom ! domain identifier
integer, intent(in) :: i ! index into array, not connected to dimId
character(len=NF90_MAX_NAME) :: get_unique_dim_name

get_unique_dim_name = state%domain(dom)%unique_dim_names(i)

end function get_unique_dim_name

!-------------------------------------------------------
function get_unique_dim_length(dom, i)

integer, intent(in) :: dom ! domain identifier
integer, intent(in) :: i ! index into array, not connected to dimId
integer :: get_unique_dim_length

get_unique_dim_length = state%domain(dom)%unique_dim_length(i)

end function get_unique_dim_length

!-------------------------------------------------------
!> Adding space for an unlimited dimension in the dimesion arrays
!> The unlimited dimension needs to be last in the list for def_var
subroutine add_time_unlimited(unlimited_dimId)

integer, intent(in)  :: unlimited_dimId

integer :: i !> loop variable

call error_handler(E_ERR, 'add_time_unlimited', 'does not exist')

end subroutine add_time_unlimited

!------------------------------------
function get_sum_variables_below(start_var, domain)

integer, intent(in) :: start_var
integer, intent(in) :: domain
integer             :: get_sum_variables_below

integer :: i

get_sum_variables_below = 0

do i = 1, domain -1
   get_sum_variables_below = get_sum_variables_below + state%domain(i)%size
enddo

do i = 1, start_var -1
   get_sum_variables_below = get_sum_variables_below + state%domain(domain)%variable(i)%var_size
enddo

end function get_sum_variables_below

!------------------------------------
function get_sum_variables(start_var, end_var, domain)

integer, intent(in) :: start_var, end_var
integer, intent(in) :: domain
integer             :: get_sum_variables

integer :: i

get_sum_variables = 0

do i = start_var, end_var
   get_sum_variables = get_sum_variables + state%domain(domain)%variable(i)%var_size
enddo

end function get_sum_variables

!------------------------------------

end module state_structure_mod
