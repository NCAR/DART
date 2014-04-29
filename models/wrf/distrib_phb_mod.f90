!> Aim: to replace the array phb, with a distributed version
!> Current thinking: Only distribute the array as much as necessary.
!> Ideally, if memory was unliimited, every task would have the whole array.
!>
!> Splitting in x only, because this is simple, and I can't think
!> of a reason to split in more than one dimension yet.
module distrib_phb_mod

use mpi_utilities_mod,     only : datasize, my_task_id, task_count
use location_mod,          only : location_type
use types_mod,             only : r8
use pnetcdf_utilities_mod, only : pnet_check

use mpi

implicit none
private

public :: create_groups, read_phb, phb_filename, create_window, free_window
public :: phb ! for the model to see

! grid groups global to module

! grid window
real(r8), allocatable :: phb_array(:,:,:) !< local phb info
integer  :: phb_win
real(r8) :: duplicate_phb(*)
pointer(aa, duplicate_phb)

integer :: group_size = 4 !< should be namelist option
integer :: group_all !< mpi_comm_world group
integer :: subgroup !< subgroup for the grid
integer :: mpi_comm_grid
integer, allocatable :: group_members(:)
integer :: local_rank !< rank within group

! pntecdf variables
character*50                  :: phb_filename !< netcdf file containing grid info
integer                       :: ret !< return code for pnetcdf
integer                       :: ncfile !< ncfile
integer                       :: varId !< variable Id
integer(KIND=MPI_OFFSET_KIND) :: we_length, sn_length, bts_length !< lengths of dimensions
integer(KIND=MPI_OFFSET_KIND) :: start(4)
integer(KIND=MPI_OFFSET_KIND) :: count(4)
integer(KIND=MPI_OFFSET_KIND) :: stride(4)
integer(KIND=MPI_OFFSET_KIND) :: my_num_we !< splitting up by we

integer :: we_dimId, sn_dimId, bts_dimId !< dimension ids

contains

!-------------------------------------------------------------
!> create groups for grid
!> make a communicator for the distributed grid
subroutine create_groups

integer ierr ! all MPI errors are fatal anyway

allocate(group_members(group_size)) ! this is module global

call mpi_comm_group(mpi_comm_world, group_all, ierr)  ! get the word group from mpi_comm_world
call build_my_group(my_task_id(), group_size, group_members) ! create a list of processors in the grid group
call mpi_group_incl(group_all, group_size, group_members, subgroup, ierr)
call mpi_comm_create(mpi_comm_world, subgroup, mpi_comm_grid, ierr)
call mpi_comm_rank(mpi_comm_grid, local_rank, ierr) ! rank within group

end subroutine create_groups

!-----------------------------------------------------------
!> build the group to store the grid
subroutine build_my_group(myrank, group_size, group_members)

implicit none

integer, intent(in)     :: myrank ! why are you passing this in?
integer, intent(inout)  :: group_size ! need to modify this if your #tasks does not divide by group size
integer, intent(out)    :: group_members(group_size)

integer bottom, top !< start and end members of the group
integer i

bottom = (myrank / group_size ) * group_size
top = bottom + group_size - 1
if (top >= task_count()) then
   top = task_count() - 1
   group_size = top - bottom + 1
   print*, 'rank', myrank, 'bottom top', bottom, top, 'group_size', group_size
endif


! fill up group members
group_members(1) = bottom
do i = 2, group_size
   group_members(i) = group_members(i-1) + 1
enddo

end subroutine build_my_group

!-----------------------------------------------------------
!> create window for phb
!> Need to add a loop for domain.  Just single domain for now.
!> Window is filled with columns contiguous in memory, possibly
!> we will grab a whole column at once.
!> 
!> If you want to change the order in the window, be sure to 
!> change the subroutine who_has_grid_info
subroutine create_window

integer ii, jj, kk, ierr, sizedouble, count
integer(KIND=MPI_ADDRESS_KIND) :: window_size

call mpi_type_size(datasize, sizedouble, ierr) ! datasize comes from mpi_utilities_mod
window_size = my_num_we*sn_length*bts_length*sizedouble
aa = malloc(my_num_we*sn_length*bts_length)
call MPI_ALLOC_MEM(window_size, mpi_info_null, aa, ierr)

! can't do array assignment with a cray pointer, so you need to loop
count = 1
! which order should these be in? (z, x, y) so you can grab a column at a time?
do jj = 1, sn_length
   do kk = 1, my_num_we
      do ii = 1, bts_length
         duplicate_phb(count) = phb_array(kk, jj, ii)
         !print*, 'count ', count
         count = count + 1
      enddo
   enddo
enddo

call mpi_win_create(duplicate_phb,  window_size, sizedouble, MPI_INFO_NULL, mpi_comm_grid, phb_win, ierr)

end subroutine create_window

!-----------------------------------------------------------
!> read in phb in parallel
!> using pnetcdf. 
!> How does the synchronization cost change with increasing task number?
!> I think it should only vary with group size, and that total task number is
!> irrelevant, but we'll see.
subroutine read_phb
use pnetcdf

integer :: temp
integer(KIND=MPI_OFFSET_KIND) :: dummy_count
integer :: my_start !< start for we

ret = nfmpi_open(mpi_comm_grid, phb_filename, NF_NOWRITE, mpi_info_null, ncfile)
call pnet_check(ret, 'read_phb', 'opening file')

! get dimensions
ret = nfmpi_inq_dimid(ncfile, 'west_east', we_dimId)
call pnet_check(ret, 'read_phb', 'getting dimension id')

ret = nfmpi_inq_dimid(ncfile, 'south_north', sn_dimId)
call pnet_check(ret, 'read_phb', 'getting dimension id')

ret = nfmpi_inq_dimid(ncfile, 'bottom_top_stag', bts_dimId)
call pnet_check(ret, 'read_phb', 'getting dimension id')

! get dimension lengths
ret = nfmpi_inq_dimlen(ncfile, we_dimId, we_length)
call pnet_check(ret, 'read_phb', 'getting dimension length')

ret = nfmpi_inq_dimlen(ncfile, sn_dimId, sn_length)
call pnet_check(ret, 'read_phb', 'getting dimension length')

ret = nfmpi_inq_dimlen(ncfile, bts_dimId, bts_length)
call pnet_check(ret, 'read_phb', 'getting dimension length')

! allocate space to read the grid
! Q.how do you split a three d variable? Is this the place to use global arrays?
! split up in x dimension only? (or y or z)
temp = we_length / group_size
if (local_rank < group_size - 1) then
   my_num_we = temp
else
   my_num_we = we_length - temp*(group_size - 1)
endif

allocate(phb_array(my_num_we, sn_length, bts_length)) ! local phb array

!print*, 'lengths ', my_num_we, sn_length, bts_length
!print*, 'rank', my_task_id(), 'total num vars ', my_num_we*sn_length*bts_length

my_start = local_rank*temp + 1

dummy_count = 1

! Why do you need the fouth dimension for pnetcdf, it seems like netcdf will just read the first three
start = (/my_start, 1, 1, 1/)
count = (/my_num_we, sn_length, bts_length, dummy_count/)
stride = (/1, 1, 1, 1/)

! read in the grid
ret = nfmpi_inq_varid(ncfile, 'PHB', varId) ! get status of variable
call pnet_check(ret, 'read_phb', 'phb id')

ret = nfmpi_get_vars_real_all(ncfile, varId, start, count, stride, phb_array)
call pnet_check(ret, 'read_phb', 'reading phb')

! check against matlab read
!print*, phb_array(1,1, 1:5)

! close netcdf file
ret = nfmpi_close(ncfile) ! close netcdf file
call pnet_check(ret, 'read_phb', 'closing file')

end subroutine read_phb

!---------------------------------------------------------
!> Free the mpi windows you created
subroutine free_window
integer :: ierr

call mpi_win_free(phb_win, ierr)
call MPI_FREE_MEM(duplicate_phb, ierr)

deallocate(group_members, phb_array)

end subroutine free_window

!---------------------------------------------------------
!> Function to get phb
function phb(dom, i, j, k)

real(r8) :: phb !< geopotential height
integer, intent(in) :: i, j, k !< 3D location
integer, intent(in) :: dom !< domain number

integer                          :: owner !< which task has the part of phb we need
integer(KIND=MPI_ADDRESS_KIND)   :: target_disp !< displacement
integer                          :: ierr

! caluclate who has the info
call who_has_grid_info(dom, i, j, k, owner, target_disp)

! grab the info
call mpi_win_lock(MPI_LOCK_SHARED, owner, 0, phb_win, ierr)
call mpi_get(phb, 1, datasize, owner, target_disp, 1, datasize, phb_win, ierr)
call mpi_win_unlock(owner, phb_win, ierr)

end function

!---------------------------------------------------------
!> get the owner and location in memory of the grid point
!> @todo How do we split up domains across the window?
subroutine who_has_grid_info(dom, i, j, k, owner, target_disp)

integer,                        intent(in)  :: dom
integer,                        intent(in)  :: i
integer,                        intent(in)  :: j
integer,                        intent(in)  :: k
integer,                        intent(out) :: owner
integer(KIND=MPI_ADDRESS_KIND), intent(out) :: target_disp !< displacement

integer :: temp
integer :: local_we !< local index of my slice of we
integer :: owner_num_we  !< length of we slab on the owner
integer :: we

temp = we_length / group_size ! phb is split only in we (west-east)

we = i-1
owner  = we / temp ! we is i

! last task in the group has the most grid points
if (owner > group_size - 1) owner = group_size - 1

if (owner < group_size - 1) then
   owner_num_we = temp
else
   owner_num_we = we_length - temp*(group_size - 1)
endif

! phb is a 3d array. (z,x,y) window
local_we = i - temp*owner  
target_disp = (j - 1)*bts_length*owner_num_we + (local_we- 1)*bts_length + k -1 ! target disp starts at 0

end subroutine who_has_grid_info

end module distrib_phb_mod
