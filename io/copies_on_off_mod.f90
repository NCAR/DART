! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module copies_on_off_mod

implicit none

interface turn_read_copy_on
   module procedure turn_read_copy_on_single
   module procedure turn_read_copy_on_range
end interface

interface turn_write_copy_on
   module procedure turn_write_copy_on_single
   module procedure turn_write_copy_on_range
end interface

interface turn_write_copy_off
   module procedure turn_write_copy_off_single
   module procedure turn_write_copy_off_range
end interface

private

public :: setup_read_write
public :: query_read_copy, query_write_copy
public :: turn_read_copy_on, turn_read_copies_off
public :: turn_write_copy_on, turn_write_copy_off


! Stores which copies to read and write
logical, allocatable :: read_copies(:), write_copies(:)

contains

!-------------------------------------------------------
!> returns true/false depending on whether you should read this copy
function query_read_copy(c)

integer, intent(in) :: c !< copy number
logical             :: query_read_copy

if ( in_range(c) ) then
   query_read_copy = read_copies(c)
else
   query_read_copy = .false.
endif

end function query_read_copy

!-------------------------------------------------------
!> returns true/false depending on whether you should read this copy
function query_write_copy(c)

integer, intent(in) :: c !< copy number
logical             :: query_write_copy

if ( in_range(c) ) then
   query_write_copy = write_copies(c)
else
   query_write_copy = .false.
endif

end function query_write_copy

!-------------------------------------------------------
!> Make the arrays for which copies to read and write
!> Default to just the actual copies, no extras
subroutine setup_read_write(num_copies)

integer, intent(in) :: num_copies

if( .not. allocated(read_copies) ) allocate(read_copies(num_copies))
if( .not. allocated(write_copies) ) allocate(write_copies(num_copies))

read_copies(:) = .false.
write_copies(:) = .false.

end subroutine setup_read_write

!-------------------------------------------------------
!> Turn on copies to read
subroutine turn_read_copy_on_single(c)

integer, intent(in) :: c !< copy to read

if (in_range(c)) read_copies(c) = .true.

end subroutine turn_read_copy_on_single

!-------------------------------------------------------
!> Turn on copies to read
subroutine turn_write_copy_on_single(c)

integer, intent(in) :: c !< copy to write

if (in_range(c)) write_copies(c) = .true.

end subroutine turn_write_copy_on_single

!-------------------------------------------------------
!> Turn on copies to read
subroutine turn_read_copy_on_range(c1, c2)

integer, intent(in) :: c1 !< start copy to read
integer, intent(in) :: c2 !< end copy to read
integer             :: f1, f2

call force_range(c1, c2, f1, f2)
read_copies(f1:f2) = .true.

end subroutine turn_read_copy_on_range

!-------------------------------------------------------
!> Turn on copies to read
subroutine turn_write_copy_on_range(c1, c2)

integer, intent(in) :: c1 !< start copy to write
integer, intent(in) :: c2 !< end copy to write
integer             :: f1, f2

call force_range(c1, c2, f1, f2)
write_copies(f1:f2) = .true.

end subroutine turn_write_copy_on_range

!-------------------------------------------------------
!> Turn off copies to read
subroutine turn_read_copies_off(c1, c2)

integer, intent(in) :: c1 !< start copy to read
integer, intent(in) :: c2 !< end copy to read
integer             :: f1, f2

call force_range(c1, c2, f1, f2)
read_copies(f1:f2) = .false.

end subroutine turn_read_copies_off

!-------------------------------------------------------
!> Turn off copies to write
subroutine turn_write_copy_off_range(c1, c2)

integer, intent(in) :: c1 !< start copy to write
integer, intent(in) :: c2 !< end copy to write
integer             :: f1, f2

call force_range(c1, c2, f1, f2)
write_copies(f1:f2) = .false.

end subroutine turn_write_copy_off_range

!-------------------------------------------------------
!> Turn off copies to write
subroutine turn_write_copy_off_single(c)

integer, intent(in) :: c
if (in_range(c)) write_copies(c) = .false.

end subroutine turn_write_copy_off_single

!-------------------------------------------------------
!> Check the copy is in range
function in_range(c)

integer, intent(in) :: c
logical             :: in_range

in_range = .true.
if ( (c > size(read_copies)) .or. (c <= 0) ) in_range = .false.

end function in_range

!-------------------------------------------------------
!> Force range to be within size of read_copies
subroutine force_range(c1, c2, f1, f2)

integer, intent(in) :: c1, c2
integer, intent(out) :: f1, f2

f1 = c1
f2 = c2
if (c1 <= 0) f1 = 1
if (c2 > size(read_copies)) f2 = size(read_copies)

end subroutine force_range

!-------------------------------------------------------
end module copies_on_off_mod