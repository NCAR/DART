! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program location_test3

! Simple test program to exercise threed_sphere location module.

use location_mod
use types_mod,      only : r8
use utilities_mod,  only : get_unit, error_handler, E_ERR, &
                           initialize_utilities, finalize_utilities
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform


implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer, parameter :: nl = 20
type(location_type) :: loc1(nl), loc2, near1
integer  :: iunit, i, j, dummy(nl), num_close, close_ind(nl), rc, near1index
real(r8) :: loc2_val(3), dist(nl)
real(r8) :: minv(3), maxv(3), v(3), minr(3), maxr(3), maxdist
type(random_seq_type) :: ran_seq
logical :: ran_seq_init = .false.
type(get_close_type) :: cc_gc
character(len=100) :: buf


call initialize_utilities('location_test3')

! Open an output file
iunit = get_unit()
open(iunit, file = 'location_test_file')

if(.not. ran_seq_init) then
   call init_random_seq(ran_seq)
   ran_seq_init = .TRUE.
endif

minv(1) = -100.0_r8
minv(2) = -200.0_r8
minv(3) = -500.0_r8

maxv(1) = 10000.0_r8
maxv(2) = 20000.0_r8
maxv(3) = 50000.0_r8

minr(:) = 1e38_r8
maxr(:) = -1e38_r8

maxdist = 10000.0_r8

! set a known min/max at the box edges so we know any
! random points generated will be inside for now.
loc1(1)  = set_location(minv(1), minv(2), minv(3))
call write_location(0, loc1(1), charstring=buf)
write(*,*) 'location ', 1, trim(buf)

loc1(2)  = set_location(maxv(1), maxv(2), maxv(3))
call write_location(0, loc1(2), charstring=buf)
write(*,*) 'location ', 2, trim(buf)

! Set the location list
do i = 3, nl
   do j=1, 3
      v(j) = random_uniform(ran_seq) * (maxv(j)-minv(j)) + minv(j)
      if (v(j) < minr(j)) minr(j) = v(j)
      if (v(j) > maxr(j)) maxr(j) = v(j)
   enddo

   loc1(i)  = set_location(v(1), v(2), v(3))

   call write_location(0, loc1(i), charstring=buf)
   write(*,*) 'location ', i, trim(buf)
enddo

! Write this location to the file
do i = 1, nl
   call write_location(iunit, loc1(i))
enddo

call get_close_init(cc_gc, nl, maxdist, loc1)

call print_get_close_type(cc_gc)

dummy = 0

! generate another random and find nearest from list
do i = 1, nl
   do j=1, 3
      v(j) = random_uniform(ran_seq) * (maxr(j)-minr(j)) + minr(j)
   enddo

   loc2  = set_location(v(1), v(2), v(3))
   call write_location(iunit, loc2)
   call write_location(0, loc2, charstring=buf)
   print *, 'generated a random point at ', trim(buf)

   do j=1, nl
      dist(j) = get_dist(loc1(j), loc2)
      if (dist(j) <= maxdist) then
         print *, 'dist to point ', j, 'is less than maxdist', dist(j)
      endif
   enddo
 

   call get_close_obs(cc_gc, loc2, 0, loc1, dummy, dummy, num_close, close_ind, dist)
   if (num_close > 0) then
      print *, 'num close = ', num_close
      do j=1, min(num_close, nl)
         print *, j, close_ind(j)
         if (close_ind(j) >= 1 .and. close_ind(j) <= nl) then
            call write_location(0, loc1(close_ind(j)), charstring=buf)
            write(*,*) 'close box loc ', trim(buf), dist(j)
         endif
      enddo
   endif

   call find_nearest(cc_gc, loc2, loc1, near1index, rc)
   if (rc /= 0) then
      print *, 'bad return from find nearest, ', rc
      cycle
   endif
   if (near1index < 1) then
      print *, 'near1index < 1, ', near1index
      cycle
   endif
   print *, 'nearest location index = ', near1index
   call write_location(0, loc1(near1index), charstring=buf)
   write(*,*) 'near loc ', trim(buf)

enddo

close(iunit)

! Now read them back in and compute the distances from loc1
open(iunit, file = 'location_test_file')

do i = 1, nl
   loc1(i) = read_location(iunit)
enddo
do i = 1, nl
   loc2 = read_location(iunit)
   write(*, *) 'distance ', i, ' is ', get_dist(loc1(i), loc2)
enddo

close(iunit)

call finalize_utilities()

end program location_test3

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
