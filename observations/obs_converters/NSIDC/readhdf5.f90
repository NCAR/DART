
program readhdf5

use HDF5

implicit none

character(len=256) :: filename = 'data/SMAP_L2_SM_P_02522_D_20150723T002821_R12170_001.h5'

integer                       :: hdferr, rank
integer(HID_T)                :: file_id, dset_id, dspace_id
integer(HSIZE_T), allocatable :: dims(:), maxdims(:)
real, allocatable             :: longitude(:)

!-----------------------------------------------------------------------
! start of executable code

! initialize the Fortran interface
call h5open_f(hdferr)
write(*,*)'h5open_f error status is ',hdferr

! open the file
call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, hdferr)
write(*,*)'h5fopen_f error, file_id ',hdferr,file_id

! open the dataset
call h5dopen_f(file_id,'/Soil_Moisture_Retrieval_Data/longitude',dset_id, hdferr)
write(*,*)'h5dopen_f error, dset_id is ',hdferr,dset_id

! open the dataspace
call h5dget_space_f(dset_id, dspace_id, hdferr)
write(*,*)'h5dget_space_f error, dspace_id ',hdferr,dspace_id

! get the rank of the dataset
call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferr)
write(*,*)'h5sget_simple_extent_ndims_f error, rank is ',hdferr, rank

allocate(dims(rank),maxdims(rank))

! fill the dims array with the dimensions
call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, hdferr)
if (hdferr < 0) then
   write(*,*)'h5sget_simple_extent_dims_f error, dims(:) is ',hdferr, dims
else
   write(*,*)'h5sget_simple_extent_dims_f  rank, dims(:) is ',hdferr, dims
endif

!allocate input data to the dimensions
allocate(longitude(dims(1)))

! read the data using the dataset ID
call h5dread_f(dset_id, H5T_NATIVE_REAL, longitude, dims, hdferr)
write(*,*)'h5dread_f with dset_id error is ',hdferr

if (hdferr == 0) write(*,*)longitude(1:10)

deallocate(longitude, dims)

end program readhdf5
