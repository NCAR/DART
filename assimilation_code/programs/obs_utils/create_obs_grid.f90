! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> program intended to create a set of observations located on
!> a regular grid.  obs have no data values - the output obs_seq
!> (suggested name obs_seq.in) must go through a program like
!> perfect_model_obs to give the observations values.

!> this code doesn't look like it was very polished; seems
!> quick and dirty.

program create_obs_grid

use    utilities_mod, only : initialize_utilities, finalize_utilities
use obs_sequence_mod, only : obs_sequence_type, interactive_obs, write_obs_seq, &
                             static_init_obs_sequence, obs_type, init_obs_sequence, &
                             set_copy_meta_data, init_obs, destroy_obs, insert_obs_in_seq, &
                             get_obs_def
use  assim_model_mod, only : static_init_assim_model
use time_manager_mod, only : time_type, operator(>)
use     location_mod, only : location_type, interactive_location

use obs_def_mod
use obs_kind_mod

implicit none

character(len=*), parameter :: source = 'create_obs_grid.f90'

type(obs_sequence_type) :: seq
character(len=256)      :: file_name

! Record the current time, date, etc. to the logfile
call initialize_utilities('create_obs_grid')

! Initialize the assim_model module, need this to get model
! state meta data for locations of identity observations
call static_init_assim_model()

! Initialize the obs_sequence module
call static_init_obs_sequence()

! Create grid of obs
seq = create_grid()

! Write the sequence to a file
write(*, *) 'Input filename for sequence (  obs_seq.in   usually works well)'
read(*, *) file_name
call write_obs_seq(seq, file_name)

! Clean up
call finalize_utilities('create_obs_grid')

contains

function create_grid()
 type(obs_sequence_type) :: create_grid

type(obs_type)     :: obs, prev_obs
type(obs_def_type) :: obs_def
type(time_type)    :: obs_time, prev_time
type(location_type) :: loc
integer            :: max_num_grids, num_copies, num_qc, end_it_all, max_num_obs
integer            :: num_dim, n(3), i, j, k, l

! these things aren't prompted for - they're fixed in the code
num_copies = 1
num_qc     = 0
max_num_obs = 1000000   ! FIXME: made up

write(*, *) 'Input upper bound on number of grids of observations in sequence'
read(*, *) max_num_grids

! Initialize an obs_sequence structure
call init_obs_sequence(create_grid, num_copies, num_qc, max_num_obs)

do i = 1, num_copies
   call set_copy_meta_data(create_grid, i, 'observations')
end do

! Initialize the obs variable
call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

! Loop to initialize each observation in turn; terminate by -1
do l = 1, max_num_grids
   write(*, *) 'input a -1 if there are no more grids'


   read(*, *) end_it_all
   if(end_it_all == -1) exit

   ! FIXME: this is the corner of a grid, need to prompt for
   ! extents in each dim (2d or 3d) and number of points in same.
   ! then loop below using same type and error, just bumping
   ! location each time.

   write(*, *) 'the location of the next observation defines the corner of a box'

   ! Need to have key available for specialized observation modules
   call interactive_obs(num_copies, num_qc, obs, i)

   write(*, *) 'enter the location of the opposite corner of the box'
   call interactive_location(loc)

   num_dim = -1
   do while (num_dim < 1 .and. num_dim > 3) 
      write(*, *) 'input 1, 2 or 3 for 1d, 2d or 3d grid'
      read(*, *) num_dim
   enddo

   n = 0
   do i=1, num_dim
      write(*,*) 'input nitems for dimension ', i
      read(*,*) n(i)
   enddo

   do k=1, n(3)
      do j=1, n(2)
         do i=1, n(1)
            ! set an obs based on the corner one
            ! compute new location and set it into an obs
! fixme here, too

            if(i == 1) then
               call insert_obs_in_seq(create_grid, obs)
            else
               ! if this is not the first obs, make sure the time is larger
               ! than the previous observation.  if so, we can start the
               ! linked list search at the location of the previous obs.
               ! otherwise, we have to start at the beginning of the entire
               ! sequence to be sure the obs are ordered correctly in
               ! monotonically increasing times. 
               call get_obs_def(obs, obs_def)
               obs_time = get_obs_def_time(obs_def)
               call get_obs_def(prev_obs, obs_def)
               prev_time = get_obs_def_time(obs_def)
               if(prev_time > obs_time) then
                  call insert_obs_in_seq(create_grid, obs)
               else
                  call insert_obs_in_seq(create_grid, obs, prev_obs)
               endif
            endif
            prev_obs = obs
         end do    ! i
      end do    ! j
   end do    ! k
end do    ! max_grids

call destroy_obs(obs)
call destroy_obs(prev_obs)

end function create_grid

end program create_obs_grid

