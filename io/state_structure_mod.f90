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
          get_sum_variables_below,&
          get_unlimited_dimid,    &
          get_model_variable_indices,      &
          get_dart_vector_index,         &
          add_dimension_to_variable, &
          finished_adding_domain

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
   integer :: var_size = 0
   integer(i8) :: index1 ! location in dart state vector of first occurrence
   integer(i8) :: indexN ! location in dart state vector of last  occurrence
   
   ! dimension information
   integer :: numdims  = 0! number of dims - excluding TIME?
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
   ! initialize all domains to no unlimited dimension
   integer :: unlimDimId = -1

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
   integer           :: num_domains = 0
   type(domain_type) :: domain(max_num_domains) ! number of domains

   ! size of the state vector
   integer(i8)       :: model_size

   !>@todo is this different to character(len=128) :: name?  Redundant?
   !character(len=1024)  :: domain_name(max_num_domains)
end type state_type

!-------------------------------------------------------------------------------
! module handle
type(state_type) :: state

! debug flag for get_state_indices
logical :: debug = .false.

interface add_domain
   module procedure add_domain_blank
   module procedure add_domain_from_file
   module procedure add_domain_from_spec
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
!> @todo dart kinds in state structure
function add_domain_from_spec(num_variables, var_names, dart_kinds) result(dom_identifier)

integer,          intent(in) :: num_variables
character(len=*), intent(in) :: var_names(num_variables)
integer,          intent(in) :: dart_kinds(num_variables)
integer           :: dom_identifier

integer :: ivar

if (.not. state_initialized) call static_init_state_type()

! add to domains
state%num_domains = state%num_domains + 1
dom_identifier    = state%num_domains !>@todo this should be a handle.

state%domain(dom_identifier)%info_file = 'NULL' ! no files

! set number of variables in this domain
state%domain(dom_identifier)%num_variables = num_variables

! load up the variable names
allocate(state%domain(dom_identifier)%variable(num_variables))

do ivar = 1, num_variables
   state%domain(dom_identifier)%variable(ivar)%varname = var_names(ivar)
enddo


end function add_domain_from_spec

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

type(domain_type), intent(inout) :: domain

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

! get the dimension id of the unlimited dimension if it exists
ret = nf90_inquire(ncfile, unlimitedDimId=domain%unlimDimId)
call nc_check(ret, 'set_unlimited_dim, nf90_inquire')

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

integer(i8) :: ind1 ! first index of variable
integer(i8) :: indN ! last index of variable

dom_size = 0
ind1 = state%model_size + 1

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

   ! first and last location of variable in the state index
   domain%variable(ivar)%index1 = ind1
   domain%variable(ivar)%indexN = ind1 + var_size - 1

   ! update counters
   dom_size = dom_size + var_size
   ind1     = domain%variable(ivar)%indexN + 1

enddo

domain%size = dom_size
state%model_size = state%model_size + dom_size

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
!> Starting dart index for variable
!-------------------------------------------------------------------------------
function get_ind1(idom,ivar)

integer, intent(in) :: idom
integer, intent(in) :: ivar
integer(i8) :: get_ind1

get_ind1 = state%domain(idom)%variable(ivar)%index1

end function get_ind1

!-------------------------------------------------------------------------------
!> Ending dart index for variable
!-------------------------------------------------------------------------------
function get_indN(idom,ivar)

integer, intent(in) :: idom
integer, intent(in) :: ivar
integer(i8) :: get_indN

get_indN = state%domain(idom)%variable(ivar)%indexN

end function get_indN

!-------------------------------------------------------------------------------
!> Return unlimited dimension id
!-------------------------------------------------------------------------------
function get_unlimited_dimid(idom)

integer, intent(in) :: idom
integer :: get_unlimited_dimid

get_unlimited_dimid = state%domain(idom)%unlimDimId

end function get_unlimited_dimid

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
! Given a dart state index, return the iloc, jloc, kloc location of the local variable
!-------------------------------------------------------------------------------
subroutine get_model_variable_indices(index_in, iloc, jloc, kloc, var_id,  dom_id)
integer(i8), intent(in)            :: index_in
integer,     intent(out)           :: iloc
integer,     intent(out)           :: jloc
integer,     intent(out)           :: kloc
integer,     intent(out), optional :: var_id
integer,     intent(out), optional :: dom_id

integer     :: local_ind
integer(i8) ::  index, ind1, indN
integer     :: ndims, ndomains
integer     :: idom, ivar
integer     :: domid, varid, unlimid
integer, allocatable :: dimids(:), dsize(:)

character(len=NF90_MAX_NAME) :: dname
character(len=NF90_MAX_NAME) :: vname

character(len=512) :: string1

! identity obs come in with a negative value - absolute
! value into the state vector.
index = abs(index_in)

varid    = -1
ndomains = get_num_domains()

! get the local variable indicies
FindVariable : do idom = 1, ndomains
   do ivar = 1, state%domain(idom)%num_variables

      ind1 = get_ind1(idom,ivar)
      indN = get_indN(idom,ivar)

      if ((index .ge. ind1) .and. &
          (index .le. indN) ) then

         local_ind = index - ind1
         varid     = ivar
         domid     = idom
         if(debug) print*, 'ind1, indN', ind1, indN 

         exit FindVariable

      endif
   enddo
enddo FindVariable

if( varid == -1 ) then
   write(string1,*) 'Problem, cannot find base offset, index_in is: ', index_in
   call error_handler(E_ERR, 'get_model_variable_indices',string1)
endif

ndims = get_num_dims(domid, varid)
allocate(dimids(ndims), dsize(ndims))

dimids  = get_dim_ids(domid,varid)
dsize   = get_dim_lengths(domid,varid)
unlimid = get_unlimited_dimid(domid)

if(debug) print*, 'ndims', ndims, 'sizes', dsize(1:ndims)

! substract the unlimited (TIME?) dimension if it exists
if(unlimid /= -1) ndims = ndims - 1

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

deallocate(dimids, dsize)

end subroutine get_model_variable_indices

!-------------------------------------------------------------------------------
! calculate dart index from local variable indices
!-------------------------------------------------------------------------------
function get_dart_vector_index(iloc, jloc, kloc, dom_id, var_id)

integer(i8) :: get_dart_vector_index
integer, intent(in) :: iloc, jloc, kloc
integer, intent(in) :: dom_id, var_id

integer :: ndims, offset
integer :: dsize(NF90_MAX_VAR_DIMS)
character(len=512) :: string1

ndims = get_num_dims(dom_id, var_id)

dsize(1:ndims) = get_dim_lengths(dom_id,var_id)

offset = get_ind1(dom_id, var_id)

if (ndims == 1) then
   get_dart_vector_index = offset + iloc - 1
else if(ndims == 2) then
   get_dart_vector_index = offset + iloc + (jloc-1)*dsize(1) - 1
else if(ndims == 3) then
   get_dart_vector_index = offset + iloc + (jloc-1)*dsize(1) + (kloc-1)*dsize(1)*dsize(2) - 1
else
   write(string1,*) 'can not calculate indices for variable ', &
     trim(get_variable_name(dom_id, var_id)), ' ndims = ', ndims
   call error_handler(E_ERR, 'get_dart_vector_index',string1)
endif

end function get_dart_vector_index

!-------------------------------------------------------------------------------
! Used to add dimenions to a variable.
! This allows the model to add meta data to state structure so a netcdf restart 
! can be created which has the T,U,V etc. from a cold start (no existing netcdf
! info file) e.g. the bgrid model. The number of variables has already been given in add_domain_from_spec.
!-------------------------------------------------------------------------------
subroutine add_dimension_to_variable(dom_id, var_id, dim_name, dim_size)

integer,          intent(in) :: dom_id
integer,          intent(in) :: var_id ! this is the order you gave in add_domain
character(len=*), intent(in) :: dim_name
integer,          intent(in) :: dim_size

integer :: d ! dimension you are adding

if (.not. state_initialized) then
   call error_handler(E_ERR, 'state_structure_mod', 'trying to add a dimension ',&
          'before initializing calling add_domain')
endif

state%domain(dom_id)%variable(var_id)%numdims = state%domain(dom_id)%variable(var_id)%numdims + 1

d = state%domain(dom_id)%variable(var_id)%numdims

state%domain(dom_id)%variable(var_id)%dimname(d) = dim_name
state%domain(dom_id)%variable(var_id)%dimlens(d) = dim_size

end subroutine add_dimension_to_variable

!-------------------------------------------------------------------------------
subroutine finished_adding_domain(dom_id)

integer, intent(in) :: dom_id ! domain identifier

integer :: i, j
integer :: count
integer :: next_start
state%domain(dom_id)%num_unique_dims = get_domain_num_dims(dom_id)

!> @todo sort by character to get unique dimensions. nancy has written this, just needs testing.

allocate(state%domain(dom_id)%unique_dim_names(get_domain_num_dims(dom_id)))
allocate(state%domain(dom_id)%unique_dim_length(get_domain_num_dims(dom_id)))

state%domain(dom_id)%size = 0
count = 1

next_start = 1

do i = 1, get_num_variables(dom_id)
   state%domain(dom_id)%variable(i)%var_size = 1
   do j = 1, state%domain(dom_id)%variable(i)%numdims

      ! product of dimensions to get variable sie
      state%domain(dom_id)%variable(i)%var_size = state%domain(dom_id)%variable(i)%dimlens(j)*state%domain(dom_id)%variable(i)%var_size

      ! unique dimensions - ignoring the fact that they are not unique
      state%domain(dom_id)%unique_dim_names(count) = state%domain(dom_id)%variable(i)%dimname(j)
      state%domain(dom_id)%unique_dim_length(count) = state%domain(dom_id)%variable(i)%dimlens(j)
      count = count + 1
   enddo

   state%domain(dom_id)%size = state%domain(dom_id)%size + state%domain(dom_id)%variable(i)%var_size

   state%domain(dom_id)%variable(i)%index1 = next_start 
   state%domain(dom_id)%variable(i)%indexN = next_start + state%domain(dom_id)%variable(i)%var_size -1
   next_start = next_start + state%domain(dom_id)%variable(i)%var_size 
enddo

state%model_size = state%model_size + state%domain(dom_id)%size

end subroutine finished_adding_domain

!-------------------------------------------------------------------------------

end module state_structure_mod
