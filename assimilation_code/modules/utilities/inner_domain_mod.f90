module inner_domain_mod

! Prototype for sharing inner domain information for PFF
! Each process has a certain number of observations, known at start
! Each observation has a certain number of inner domain variables, not easily known a priori

! This is a pragmatic implementation with limited error checking. 

use types_mod,     only : r8, i8, missing_r8

use utilities_mod, only : error_handler, E_ERR

implicit none

private

public :: inner_domain_type, init_inner_domain, clear_inner_domain,                    &
          add_obs_index_inner_domain, add_obs_ens_inner_domain, add_var_inner_domain,  &
          get_num_vars_inner_domain, get_obs_index_inner_domain,                       &
          get_obs_ens_inner_domain, get_var_index_inner_domain,                        &
          get_var_ens_inner_domain, output_inner_domain_info

! Have a fixed upper bound for the size of the inner domain for now
! Might need to make this smarter for big models with complex forward operators
integer, parameter :: max_num_vars = 10

type inner_domain_type
   integer                   :: max_num_obs, num_obs
   ! The observation index is dimensioned (max_num_obs)
   integer,      allocatable :: obs_index(:)

   ! Number of inner domain variables currently stored for each observation (max_num_obs)
   integer,      allocatable :: num_vars(:)

   ! For now, include the values of the forward operator ensemble; are these needed?
   ! Dimensioned (ens_size, max_num_obs)
   real(r8),     allocatable :: obs_ens(:, :)

   ! The inner domain indices are dimensioned (max_num_vars, max_num_obs)
   integer(i8),  allocatable :: var_index(:, :)

   ! ens is dimensioned (ens_size, max_num_vars, max_num_obs)
   real(r8),     allocatable :: ens(:, :, :)
end type inner_domain_type

type(inner_domain_type), save :: inner_domain

! Only allow inputs from get_state when computing prior forward operators
logical, save                 :: active

! Error message storage
character(len=512)          :: string1
character(len=*), parameter :: source = 'inner_domain_mod.f90'

contains

!--------------------------------------------------------------------

subroutine init_inner_domain(max_num_obs, ens_size)

! Input the maximum number of observations for this process and the ens_size
! Initialize the storage

integer, intent(in) :: max_num_obs
integer, intent(in) :: ens_size

! Now active
active = .true.

! Currently no variables stored
inner_domain%max_num_obs = max_num_obs
inner_domain%num_obs     = 0

allocate(inner_domain%obs_index(max_num_obs),         &
   inner_domain%num_vars(max_num_obs),                &
   inner_domain%obs_ens(ens_size, max_num_obs),       &
   inner_domain%var_index(max_num_vars, max_num_obs), &
   inner_domain%ens(ens_size, max_num_vars, max_num_obs))

inner_domain%num_vars  = 0
inner_domain%obs_index = -99
inner_domain%obs_ens   = missing_r8
inner_domain%var_index = -99
inner_domain%ens       = missing_r8

end subroutine init_inner_domain

!--------------------------------------------------------------------

subroutine clear_inner_domain

inner_domain%max_num_obs = 0
inner_domain%num_obs     = 0
deallocate(inner_domain%obs_index, inner_domain%num_vars, &
   inner_domain%obs_ens, inner_domain%var_index, inner_domain%ens)

! No longer active
active = .false.

end subroutine clear_inner_domain

!--------------------------------------------------------------------

subroutine add_obs_index_inner_domain(obs_index)

! Add an observation index to the inner_domain for this process

integer, intent(in) :: obs_index

! If not active don't allow input
if(.not. active) return

! Make sure that there is enough room
inner_domain%num_obs = inner_domain%num_obs + 1
if(inner_domain%num_obs > inner_domain%max_num_obs) then
   string1 = 'Exceeded max_num_obs in inner_domain mod'
   call error_handler(E_ERR, 'add_obs_index_inner_domain', string1, source)
endif

inner_domain%obs_index(inner_domain%num_obs) = obs_index

end subroutine add_obs_index_inner_domain

!--------------------------------------------------------------------

subroutine add_obs_ens_inner_domain(obs_index, obs_vals)

! Add an observation to the inner_domain
! Assumes that this observation index has been added most recently

integer,  intent(in) :: obs_index
real(r8), intent(in) :: obs_vals(:)

! If not active don't allow input
if(.not. active) return

! Make sure that the index agrees with the one to be loaded
if(inner_domain%obs_index(inner_domain%num_obs) .ne. obs_index) then
   string1 = 'Indexing error'
   call error_handler(E_ERR, 'add_obs_ens_inner_domain', string1, source)
endif

inner_domain%obs_ens(:, inner_domain%num_obs) = obs_vals

end subroutine add_obs_ens_inner_domain

!--------------------------------------------------------------------

subroutine add_var_inner_domain(var_index, vals)

! Add the index and the values of an inner domain variable for current ob
! Could add in an optional n_obs, too, but don't think this is known where
! vars will be added

integer(i8), intent(in)  :: var_index
real(r8),    intent(in)  :: vals(:)

integer :: cur_obs

! If not active don't allow input
if(.not. active) return

! Make sure there is room for another variable
cur_obs = inner_domain%num_obs
inner_domain%num_vars(cur_obs) = inner_domain%num_vars(cur_obs) + 1
if(inner_domain%num_vars(cur_obs) > max_num_vars) then
   string1 = 'Exceeded max_num_vars'
   call error_handler(E_ERR, 'add_var_inner_domain', string1, source)
endif

inner_domain%var_index(inner_domain%num_vars(cur_obs), cur_obs) = var_index
inner_domain%ens(:, inner_domain%num_vars(cur_obs), cur_obs)    = vals

end subroutine add_var_inner_domain

!--------------------------------------------------------------------

function get_num_vars_inner_domain(n_obs)

! Get the number of inner domain variables for the nth observation

integer             :: get_num_vars_inner_domain
integer, intent(in) :: n_obs

get_num_vars_inner_domain = inner_domain%num_vars(n_obs)

end function get_num_vars_inner_domain

!--------------------------------------------------------------------

subroutine get_obs_index_inner_domain(obs_index, num_obs)

integer, intent(out) :: obs_index(:)
integer, intent(out) :: num_obs

num_obs = inner_domain%num_obs

! For development, check to make sure there is enough room
if(size(obs_index) >= num_obs) then
   obs_index(1:num_obs) = inner_domain%obs_index(1:num_obs)
else
   string1 = 'obs_index too big'
   call error_handler(E_ERR, 'get_obs_index_inner_domain', string1, source)
endif

end subroutine get_obs_index_inner_domain

!--------------------------------------------------------------------

subroutine get_obs_ens_inner_domain(n_obs, obs_ens)

! Gets the ensemble of forward operator values for the nth observation in the list

integer,  intent(in)  :: n_obs
real(r8), intent(out) :: obs_ens(:)

! Check for out of range
if(n_obs > inner_domain%num_obs) then
   string1 = 'input n_obs greater than nunber of obs available'
   call error_handler(E_ERR, 'get_obs_ens_inner_domain', string1, source)
endif

obs_ens = inner_domain%obs_ens(:, n_obs)

end subroutine get_obs_ens_inner_domain

!--------------------------------------------------------------------

subroutine get_var_index_inner_domain(n_obs, var_index, num_vars)

! Returns the inner domain variable indices and the number of inner domain vars
! for the nth observation in the list

integer, intent(in)  :: n_obs
integer, intent(out) :: var_index(:)
integer, intent(out) :: num_vars

! Check for out of range
if(n_obs > inner_domain%num_obs) then
   string1 = 'input n_obs greater than nunber of obs available'
   call error_handler(E_ERR, 'get_var_index_inner_domain', string1, source)
endif

! Should put in more error checks for storage size for initial testing?
num_vars              = inner_domain%num_vars(n_obs)
var_index(1:num_vars) = inner_domain%var_index(1:num_vars, n_obs)

end subroutine get_var_index_inner_domain

!--------------------------------------------------------------------

subroutine get_var_ens_inner_domain(n_obs, n_var, ens)

! Returns the inner domain ensemble for the n_obs observation and its
! n_var inner domain variable

integer,  intent(in)  :: n_obs
integer,  intent(in)  :: n_var
real(r8), intent(out) :: ens(:)

! Lots more error checking for now
if(n_obs > inner_domain%num_obs) then
   string1 = 'n_obs exceeds num_obs'
   call error_handler(E_ERR, 'get_var_ens_inner_domain', string1, source)
endif

if(n_var > inner_domain%num_vars(n_obs)) then
   string1 = 'n_var exceeds num_vars'
   call error_handler(E_ERR, 'get_var_ens_inner_domain', string1, source)
endif

ens = inner_domain%ens(:, n_var, n_obs)

end subroutine get_var_ens_inner_domain

!--------------------------------------------------------------------

subroutine output_inner_domain_info(channel)

! Debugging routine to print out the information in an inner_domain

integer, intent(in) :: channel

integer :: i, j

! Output the sizes and max sizes
write(channel, *) 'Inner domain info'
write(channel, *) inner_domain%num_obs, ' observations with max_num_obs ', &
   inner_domain%max_num_obs

do i = 1, inner_domain%num_obs
   write(channel, *) 'observation ', i, ' has index ', inner_domain%obs_index(i)
   ! Loop through associated inner domain variables
   do j = 1, inner_domain%num_vars(i) 
      write(channel, *) 'inner domain var ', j, ' is index ', &
         inner_domain%var_index(j, i)
   end do
end do

end subroutine output_inner_domain_info

!--------------------------------------------------------------------

end module inner_domain_mod
