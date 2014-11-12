! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module copies_on_off_mod

! Stores which copies to read and write
logical, allocatable :: read_copies(:), write_copies(:)

contains

!-------------------------------------------------------
!> returns true/false depending on whether you should read this copy
function query_read_copy(c)

integer, intent(in) :: c !< copy number
logical             :: query_read_copy

if (c > size(read_copies) ) then
   query_read_copy = .false.
else
   query_read_copy = read_copies(c)
endif

end function query_read_copy

!-------------------------------------------------------
!> returns true/false depending on whether you should read this copy
function query_write_copy(c)

integer, intent(in) :: c !< copy number
logical             :: query_write_copy

if ( c > size(write_copies) ) then
   query_write_copy = .false.
else
   query_write_copy = write_copies(c)
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

read_copies(c) = .true.

end subroutine turn_read_copy_on_single

!-------------------------------------------------------
!> Turn on copies to read
subroutine turn_write_copy_on_single(c)

integer, intent(in) :: c !< copy to write

write_copies(c) = .true.

end subroutine turn_write_copy_on_single

!-------------------------------------------------------
!> Turn on copies to read
subroutine turn_read_copy_on_range(c1, c2)

integer, intent(in) :: c1 !< start copy to read
integer, intent(in) :: c2 !< end copy to read

read_copies(c1:c2) = .true.

end subroutine turn_read_copy_on_range

!-------------------------------------------------------
!> Turn on copies to read
subroutine turn_write_copy_on_range(c1, c2)

integer, intent(in) :: c1 !< start copy to write
integer, intent(in) :: c2 !< end copy to write

write_copies(c1:c2) = .true.

end subroutine turn_write_copy_on_range

!-------------------------------------------------------
!> Turn off copies to read
subroutine turn_read_copies_off(c1, c2)

integer, intent(in) :: c1 !< start copy to read
integer, intent(in) :: c2 !< end copy to read

read_copies(c1:c2) = .false.

end subroutine turn_read_copies_off

!-------------------------------------------------------
!> Turn off copies to write
subroutine turn_write_copies_off(c1, c2)

integer, intent(in) :: c1 !< start copy to write
integer, intent(in) :: c2 !< end copy to write

write_copies(c1:c2) = .false.

end subroutine turn_write_copies_off

!-------------------------------------------------------

end module copies_on_off_mod