! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module copies_on_off_mod
!> \defgroup copies_on_off_mod copies_on_off_mod
!> @{
!> @brief IO storage module
!>
!> This module stores 2 things necessary for state IO:
!>    1. Which copies to read and write.
!>    2. Which copy number is associated with each named copy (e.g POST_INF_COPY = 5)
!>
!> Usage for read:
!>    call setup_read_write(num_copies)  
!>    call turn_read_copy_on(1,ens_size)
!>    call turn_read_copy_on(mean)
!>    --- IO is done ---
!>    call end_read_write
!>  
!>  The IO routines use 
!>     query_read_copy(copy) 
!> to find out whether a copy needs to be read.
!>
!> Usage for write:
!>    call setup_read_write(num_copies)  
!>    call turn_write_copy_on(1:ens_size)
!>    call turn_write_copy_on(mean)
!>    --- IO is done ---
!>    call end_read_write
!>  
!>  The IO routines use 
!>     query_write_copy(copy)
!> to find out whether a copy needs to be written.
!> turn_read/write_copy_on/off accepts single values or ranges (n,m) = (n:m)
!> Any numbers outside the range of 1:num_copies are ignored.
!> 
!> The named copy values are initialized to COPY_NOT_PRESENT
!> Filter is the only code that uses the named copies. The extra
!> copies are given values by set_state_copies (in filter_mod.f90).
!> Note there is no protection against code that uses this module 
!> changing the values of the named copy values. It may be beneficial
!> to make these routines part of io_filenames_mod. The named copies
!> could be part of the file_info_type.

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

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

public :: setup_read_write, end_read_write
public :: query_read_copy, query_write_copy
public :: turn_read_copy_on, turn_read_copies_off
public :: turn_write_copy_on, turn_write_copy_off
public :: is_inflation_copy, is_mean_copy, is_sd_copy, is_ensemble_copy

! These are public integers for now.  
! There is no protection against code using this module changing the
! value of these integers.
public :: ENS_MEAN_COPY, ENS_SD_COPY, &
          PRIOR_INF_COPY, PRIOR_INF_SD_COPY, &
          POST_INF_COPY, POST_INF_SD_COPY, &
          SPARE_PRIOR_MEAN, SPARE_PRIOR_SPREAD, &
          SPARE_PRIOR_INF_MEAN, SPARE_PRIOR_INF_SPREAD, &
          SPARE_POST_INF_MEAN, SPARE_POST_INF_SPREAD, &
          query_copy_present

! Used to test if a copy is not in use, e.g. the spare copies may not be in use.
integer, parameter :: COPY_NOT_PRESENT = -1

!>@todo FIXME : this should be a derived type with a copy number and a long name
!>@             associated with the file.

! Global indices into ensemble storage for state
! These are initialized to to COPY_NOT_PRESENT
integer :: ENS_MEAN_COPY          = COPY_NOT_PRESENT
integer :: ENS_SD_COPY            = COPY_NOT_PRESENT
integer :: PRIOR_INF_COPY         = COPY_NOT_PRESENT
integer :: PRIOR_INF_SD_COPY      = COPY_NOT_PRESENT
integer :: POST_INF_COPY          = COPY_NOT_PRESENT
integer :: POST_INF_SD_COPY       = COPY_NOT_PRESENT
! For large models with a single time step
integer :: SPARE_PRIOR_MEAN       = COPY_NOT_PRESENT
integer :: SPARE_PRIOR_SPREAD     = COPY_NOT_PRESENT
integer :: SPARE_PRIOR_INF_MEAN   = COPY_NOT_PRESENT
integer :: SPARE_PRIOR_INF_SPREAD = COPY_NOT_PRESENT
integer :: SPARE_POST_INF_MEAN    = COPY_NOT_PRESENT
integer :: SPARE_POST_INF_SPREAD  = COPY_NOT_PRESENT


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
!> returns true/false depending on whether you should write this copy
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
!> Destroys existing read_copies write_copies
subroutine setup_read_write(num_copies)

integer, intent(in) :: num_copies !< total size of the ensemble

if (allocated(read_copies))  deallocate(read_copies)
if (allocated(write_copies)) deallocate(write_copies)

allocate(read_copies(num_copies))
allocate(write_copies(num_copies))

read_copies(:) = .false.
write_copies(:) = .false.

end subroutine setup_read_write

!-------------------------------------------------------
!> Destroy the arrays for which copies to read and write
subroutine end_read_write()

if( allocated(read_copies) ) deallocate(read_copies)
if( allocated(write_copies) ) deallocate(write_copies)

end subroutine end_read_write

!-------------------------------------------------------
!> Turn on copies to read
subroutine turn_read_copy_on_single(c)

integer, intent(in) :: c !< copy to read

if (in_range(c)) read_copies(c) = .true.

end subroutine turn_read_copy_on_single

!-------------------------------------------------------
!> Turn on copies to write
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
!> Turn on copies to write
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
!> Assumes read_copies and write_copies are the same size
subroutine force_range(c1, c2, f1, f2)

integer, intent(in) :: c1, c2
integer, intent(out) :: f1, f2

f1 = c1
f2 = c2
if (c1 <= 0) f1 = 1
if (c2 > size(read_copies)) f2 = size(read_copies)

end subroutine force_range

!------------------------------------------------------------------
! Test whether a copy is part of the ensemble
function query_copy_present(copy)

integer, intent(in) :: copy
logical :: query_copy_present

if (copy == COPY_NOT_PRESENT) then
   query_copy_present = .false.
else
   query_copy_present = .true.
endif

end function

!------------------------------------------------------------------
! Test whether the copy is an ensemble member
!>@todo FIXME : There is probably a better way to do this. This is strictly 
!>              assuming the ensemble members are before ENS_MEAN_COPY
function is_ensemble_copy(copy)
integer, intent(in) :: copy
logical :: is_ensemble_copy

is_ensemble_copy= .false.

if( copy < ENS_MEAN_COPY ) is_ensemble_copy= .true.

end function is_ensemble_copy

!------------------------------------------------------------------
! Test whether the copy is ensemble sd
function is_sd_copy(copy)
integer, intent(in) :: copy
logical :: is_sd_copy

is_sd_copy = .false.

if( copy == ENS_SD_COPY ) is_sd_copy = .true.

end function is_sd_copy

!------------------------------------------------------------------
! Test whether the copy is ensemble mean
function is_mean_copy(copy)
integer, intent(in) :: copy
logical :: is_mean_copy

is_mean_copy = .false.

if( copy == ENS_MEAN_COPY ) is_mean_copy = .true.

end function is_mean_copy

!------------------------------------------------------------------
! Test whether the copy is inflation copy
function is_inflation_copy(copy)
integer, intent(in) :: copy
logical :: is_inflation_copy

is_inflation_copy = .false.

if( copy == PRIOR_INF_COPY         .or. &
    copy == PRIOR_INF_SD_COPY      .or. &
    copy == POST_INF_COPY          .or. &
    copy == POST_INF_SD_COPY       .or. &
    copy == SPARE_PRIOR_INF_MEAN   .or. &
    copy == SPARE_PRIOR_INF_SPREAD .or. &
    copy == SPARE_POST_INF_MEAN    .or. &
    copy == SPARE_POST_INF_SPREAD ) then
   is_inflation_copy = .true.
endif

end function is_inflation_copy

!-------------------------------------------------------
end module copies_on_off_mod
!> @}
