module obs_kind_mod

! This module is designed to provide general information about observation 
! types. It is not clear at present whether this is going to be viable or
! not given the requirement of underlying models and domains. For now, kind
! simply is an index that identifies what kind of observation this is 
! for instance temperature, pressure, satellite radiance, etc. This may be
! a good place to store other parameters associated with more complicated
! observation types.

use types_mod

private

public get_obs_kind, set_obs_kind, write_kind, read_kind, obs_kind_type, &
   interactive_kind, IDENTITY_OBSERVATION

type obs_kind_type
   private
   integer :: index
end type obs_kind_type

integer, parameter :: IDENTITY_OBSERVATION = -1

contains

!----------------------------------------------------------------------------

function get_obs_kind(obs_kind)

! Returns integer index (for now, may be string or something else later)
! from this obs_kind.

implicit none

integer :: get_obs_kind
type(obs_kind_type), intent(in) :: obs_kind

get_obs_kind = obs_kind%index

end function get_obs_kind

!----------------------------------------------------------------------------

function set_obs_kind(index)

! Constructor for obs_kind, takes integer index as input for now.

implicit none

type(obs_kind_type) :: set_obs_kind
integer, intent(in) :: index

set_obs_kind%index = index

end function set_obs_kind

!----------------------------------------------------------------------------

subroutine write_kind(file, kind)

! Writes out kind to file
implicit none

integer, intent(in) :: file
type(obs_kind_type), intent(in) :: kind

! For now output a character tag followed by the integer index.
write(file, 11)
11 format('kind ')
write(file, *) kind%index

end subroutine write_kind

!----------------------------------------------------------------------------

function read_kind(file)

! Reads a kind from file

type(obs_kind_type) :: read_kind
integer, intent(in) :: file

character*5 :: header

! Need additional error checks
read(file, 11) header
11 format(a5)
if(header /= 'kind ') then
   write(*, *) 'Error: Expected kind header "kind " in input file'
   stop
endif

! Now read the kind index
read(file, *) read_kind%index

end function read_kind


subroutine interactive_kind(kind)
!-------------------------------------------------------------------------
!
! Allows interactive input of kind. For now, this module only allows
! a single kind so there is no need to do anything except set to 1

implicit none

type(obs_kind_type), intent(out) :: kind

kind%index = 1

end subroutine interactive_kind



end module obs_kind_mod
