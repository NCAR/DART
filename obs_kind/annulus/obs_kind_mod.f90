! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module obs_kind_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! This module is designed to provide general information about observation 
! types for the MITgcm annulus model.

! Revised obs_kind treats identity observations (including the index) as 
! just another kind of observation. Identity observations of an extended 
! state vector are just represented by indices that go beyond the state
! vector size into the extended state vector.
! Initial revision 15 April, 2004

use     types_mod, only : r8
use utilities_mod, only : register_module, error_handler, E_ERR

implicit none
private

public :: get_obs_kind, set_obs_kind, write_kind, read_kind, obs_kind_type, &
       interactive_kind

public :: KIND_U, KIND_V, KIND_PS, KIND_T, KIND_QV, KIND_P, KIND_W, KIND_QR, KIND_TD, &
       KIND_VR, KIND_REF

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type obs_kind_type
   private
   integer :: index
end type obs_kind_type

! ADD A LONG TABLE OF DEFINED BUFR INDICES, ETC.

! have to include kind_vr (radar) for the moment, otherwise the generic
! obs_def_mod barfs
integer, parameter :: KIND_U = 1, KIND_V = 2, KIND_PS = 3, KIND_T = 4,   &
                      KIND_QV = 5, KIND_P = 6, KIND_W = 7, KIND_QR = 8, KIND_TD = 10, &
                      KIND_VR = 100, KIND_REF = 101

logical, save :: module_initialized = .false.

contains


!----------------------------------------------------------------------------

subroutine initialize_module

   call register_module(source, revision, revdate)
   module_initialized = .true.

end subroutine initialize_module



!----------------------------------------------------------------------------

function get_obs_kind(obs_kind)

! Returns integer index (for now, may be string or something else later)
! from this obs_kind.

implicit none

integer :: get_obs_kind
type(obs_kind_type), intent(in) :: obs_kind

if ( .not. module_initialized ) call initialize_module

get_obs_kind = obs_kind%index

end function get_obs_kind

!----------------------------------------------------------------------------

function set_obs_kind(index)

! Constructor for obs_kind, takes integer index as input for now.

implicit none

type(obs_kind_type) :: set_obs_kind
integer, intent(in) :: index

if ( .not. module_initialized ) call initialize_module

set_obs_kind%index = index

end function set_obs_kind

!----------------------------------------------------------------------------

subroutine write_kind(file, kind, fform)

! Writes out kind to file
implicit none

integer, intent(in) :: file
type(obs_kind_type), intent(in) :: kind
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    !supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! For now output a character tag followed by the integer index.
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      write(file) kind%index
   CASE DEFAULT
      write(file, '(''kind'')' )
      write(file, *) kind%index
END SELECT

end subroutine write_kind

!----------------------------------------------------------------------------

function read_kind(file, fform)

! Reads a kind from file

type(obs_kind_type) :: read_kind
integer, intent(in) :: file
character(len=*), intent(in), optional :: fform

character(len=5) :: header
character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    !supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Need additional error checks
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      read(file) read_kind%index
   CASE DEFAULT
      read(file, '(a5)' ) header
      if(header /= 'kind ') then
         call error_handler(E_ERR,'read_kind', 'Expected kind header "kind " in input file', &
                            source, revision, revdate)
      endif
! Now read the kind index
      read(file, *) read_kind%index
END SELECT


end function read_kind


subroutine interactive_kind(kind)
!-------------------------------------------------------------------------
!
! Allows interactive input of kind.

implicit none

type(obs_kind_type), intent(out) :: kind

if ( .not. module_initialized ) call initialize_module

write(*, *) 'input obs kind: u = ',kind_u,' v = ',kind_v,'  t = ',kind_t, &
            ' p = ',kind_p,' w = ',kind_w

write(*, *) 'input -1 times the state variable index for an identity observation'
read(*, *) kind%index

end subroutine interactive_kind

end module obs_kind_mod
