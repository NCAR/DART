module obs_kind_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!

! This module is designed to provide general information about observation 
! types. It is not clear at present whether this is going to be viable or
! not given the requirement of underlying models and domains. For now, kind
! simply is an index that identifies what kind of observation this is 
! for instance temperature, pressure, satellite radiance, etc. This may be
! a good place to store other parameters associated with more complicated
! observation types.

use types_mod, only : r8

implicit none
private

public get_obs_kind, set_obs_kind, write_kind, read_kind, obs_kind_type, &
   interactive_kind, IDENTITY_OBSERVATION, set_ncep_obs_kind

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

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
write(file, '(''kind'')' )
write(file, *) kind%index

end subroutine write_kind

!----------------------------------------------------------------------------

function read_kind(file)

! Reads a kind from file

type(obs_kind_type) :: read_kind
integer, intent(in) :: file

character(len=5) :: header

! Need additional error checks
read(file, '(a5)' ) header
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
! Allows interactive input of kind.

implicit none

type(obs_kind_type), intent(out) :: kind

! For bgrid, enter kind of observation
write(*, *) 'input obs kind: u = 1, v = 2, ps = 3, t = 4 q = 5'
read(*, *) kind%index



end subroutine interactive_kind


function set_ncep_obs_kind(obsindex)

implicit none

type(obs_kind_type) :: set_ncep_obs_kind
integer, intent(in) :: obsindex

! Now set the kind index
 set_ncep_obs_kind%index = obsindex

end function set_ncep_obs_kind


end module obs_kind_mod
