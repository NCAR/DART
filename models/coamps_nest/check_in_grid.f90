! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

! Program wrapper for ll2ij that takes its input from the file
! domain.dat - this allows us to quickly cal ll2ij many times for a
! single grid.  This is written to print out either a "T" or "F" if
! the lat/lon point described by the last two lines of domain.dat is
! in the grid or not - this is for use with actual observations to
! pare down the list to points that are actually within our domain.

program check_in_grid

  use coamps_intrinsic_mod, only : ll2ij

  use types_mod, only : r8

  implicit none

  character(len=10), parameter :: FILENAME = 'domain.dat'
  integer                      :: fileid
  real(kind=r8), dimension(16)  :: data

  integer :: ii
  real(kind=r8) :: grid_i, grid_j
  logical :: in_grid

  ! Pull the numbers out of the file - both domain information and
  ! also the lat/lon point to convert
  fileid = 13
  open(unit=fileid,file=FILENAME);
  read(fileid,*) data
  close(fileid)

  call ll2ij(int(data(1)), data(5), data(6), int(data(13)),       &
             int(data(14)), data(2), data(3), data(4), data(7),   &
             data(8), (/ grid_i /), (/ grid_j /), 1, (/data(15)/),&
             (/ data(16) /))

  ! Check using optimism
  in_grid = .true.
  if (grid_i > data(9))  in_grid = .false.
  if (grid_j > data(10)) in_grid = .false.
  if (grid_i < 1)        in_grid = .false.
  if (grid_j < 1)        in_grid = .false.
  print *, in_grid
end program check_in_grid

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
