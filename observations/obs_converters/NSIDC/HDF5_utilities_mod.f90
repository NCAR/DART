! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module HDF5_utilities_mod

use        types_mod, only : i2, i4, r4, r8, MISSING_R8, MISSING_I
use    utilities_mod, only : register_module, E_MSG, E_ERR, error_handler
use time_manager_mod, only : time_type, operator(>=), set_time, get_time

use HDF5
! or maybe .... include HDF5

implicit none
private

public :: h5_open, H5_CRTDAT, H5_RDWT, &
          h5_get_rank, &
          h5_get_dimensions, &
          h5_get_dset_dspace, &
          h5_check

!interface hf_get_var
!   module procedure hf_get_int_1d
!   module procedure hf_get_real_1d
!end interface

character(len=*), parameter :: source   = "HDF5_utilities_mod.f90"
character(len=512) :: string1, string2, string3

logical :: module_initialized = .false.


contains


!-----------------------------------------------------------------------
!> initialize the Fortran interface to HDF5

subroutine initialize_module(context)
character(len=*), optional :: context

integer :: hdferr

if (module_initialized) return

module_initialized = .true.

if (present(context)) then
   write(string1,*)'initializing Fortran interfaces: ',trim(context)
else
   write(string1,*)'initializing Fortran interfaces'
endif

! initialize the Fortran interface
call h5open_f(hdferr)
call h5_check(hdferr,'initialize_module','h5open_f',string1) 

end subroutine initialize_module 


!-----------------------------------------------------------------------
!> open the fortran interface to hdf5 libs 
!> open the file

function h5_open(filename, flag, context) result(file_id)

character(len=*),           intent(in) :: filename
integer(HID_T),             intent(in) :: flag
character(len=*), optional, intent(in) :: context
integer(HID_T)                         :: file_id

integer :: hdferr

if ( .not. module_initialized ) call initialize_module(context)

call h5fopen_f(filename, flag, file_id, hdferr)
call h5_check(hdferr,'h5_open','h5fopen_f', context, filename)

end function h5_open


!-----------------------------------------------------------------------
!>

subroutine h5_get_dset_dspace(file_id, dsetname, dsetid, dspaceid, context)

integer(HID_T),             intent(in)  :: file_id  !< hdf file ID
character(len=*),           intent(in)  :: dsetname !< dataset name
integer(HID_T),             intent(out) :: dsetid   !< dataset ID
integer(HID_T),             intent(out) :: dspaceid !< dataset dataspace ID
character(len=*), optional, intent(in)  :: context

integer :: hdferr
character(len=*), parameter :: routine = 'h5_get_dset_dspace'

if ( .not. module_initialized ) call initialize_module(context)

call h5dopen_f(file_id, dsetname, dsetid, hdferr)
call h5_check(hdferr, routine, 'h5dopen_f', context, dsetname)

call h5dget_space_f(dsetid, dspaceid, hdferr)
call h5_check(hdferr, routine, 'h5dget_space_f', context, dsetname)

end subroutine h5_get_dset_dspace


!-----------------------------------------------------------------------
!>

function h5_get_rank(dspaceid, context) result(rank)

integer(HID_T),             intent(in)  :: dspaceid !< dataset dataspace ID
character(len=*), optional, intent(in)  :: context
integer                                 :: rank     !< number of dimensions

integer :: hdferr
character(len=*), parameter :: routine = 'h5_get_rank'

call h5sget_simple_extent_ndims_f(dspaceid, rank, hdferr)

call h5_check(hdferr, routine, 'h5sget_simple_extent_ndims_f', context)

end function h5_get_rank


!-----------------------------------------------------------------------
!>

subroutine h5_get_dimensions(dspaceid, dims, maxdims, context)

integer(HID_T),             intent(in)  :: dspaceid
integer(HSIZE_T),           intent(out) :: dims(:)    !< actual   dimension lengths
integer(HSIZE_T), optional, intent(out) :: maxdims(:) !< declared dimension lengths 
character(len=*), optional, intent(in)  :: context

integer :: hdferr
character(len=*), parameter :: routine = 'h5_get_dimensions'

integer(HSIZE_T) :: declared_dimensions(size(dims))

! get the dimensions of the dataspace
call h5sget_simple_extent_dims_f(dspaceid, dims, declared_dimensions, hdferr)

! hdferr is 'Dataspace rank on success' and -1 on failure
! h5_check only thinks 0 is a success

if (hdferr < 0) &
   call h5_check(hdferr, routine, 'h5sget_simple_extent_dims_f', context)

! write(*,*)'TJH actual   dimensions ',dims
! write(*,*)'TJH declared dimensions ',declared_dimensions

if (present(maxdims)) maxdims = declared_dimensions

end subroutine h5_get_dimensions


!------------------------------------------------------------------
!> check return code from previous call. on error, print and stop.
!> if you want to continue after an error don't use this call. 
!> a negative hdferr is considered a failure

subroutine h5_check(hdferr, subr_name, h5routine, context, filename)

integer,          intent(in)           :: hdferr
character(len=*), intent(in)           :: subr_name
character(len=*), intent(in)           :: h5routine
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=512) :: string1, string2, string3

if (hdferr >= 0) return

! something wrong.  construct an error string, print and abort.

write(string1,*)'HDF5 ERROR from ', trim(h5routine)

if (present(context)) then
   write(string2,*)trim(context),', error code is ',hdferr
else
   write(string2,*)'error code is ',hdferr
endif

call error_handler(E_ERR, subr_name, string2, &
                   source, text2=string2, text3=filename)

end subroutine h5_check


subroutine H5_CRTDAT()

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!   Copyright by The HDF Group.                                               *
!   Copyright by the Board of Trustees of the University of Illinois.         *
!   All rights reserved.                                                      *
!                                                                             *
!   This routine is part of HDF5.  The full HDF5 copyright notice, including     *
!   terms governing use, modification, and redistribution, is contained in    *
!   the COPYING file, which can be found at the root of the source code       *
!   distribution tree, or in https://support.hdfgroup.org/ftp/HDF5/releases.  *
!   If you do not have access to either file, you may request a copy from     *
!   help@hdfgroup.org.                                                        *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! The following example shows how to create an empty dataset.
! It creates a file called 'dsetf.h5', defines the
! dataset dataspace, creates a dataset which is a 4x6 integer array,
! and then closes the dataspace, the dataset, and the file.
!
! This example is used in the HDF5 Tutorial.

character(len=8), parameter :: filename = "dsetf.h5" ! File name
character(len=4), parameter :: dsetname = "dset"     ! Dataset name

integer(HID_T) :: file_id       ! File identifier
integer(HID_T) :: dset_id       ! Dataset identifier
integer(HID_T) :: dspace_id     ! Dataspace identifier

integer(HSIZE_T), dimension(2) :: dims = (/4,6/) ! Dataset dimensions
integer     ::   rank = 2                        ! Dataset rank

integer     ::   error ! Error flag

! Initialize FORTRAN interface.
call h5open_f(error)

! Create a new file using default properties.
call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

! Create the dataspace.
call h5screate_simple_f(rank, dims, dspace_id, error)

! Create the dataset with default properties.
call h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)

! End access to the dataset and release resources used by it.
call h5dclose_f(dset_id, error)

! Terminate access to the data space.
call h5sclose_f(dspace_id, error)

! Close the file.
call h5fclose_f(file_id, error)

! Close FORTRAN interface.
call h5close_f(error)

end subroutine H5_CRTDAT


subroutine H5_RDWT

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!   Copyright by The HDF Group.                                               *
!   Copyright by the Board of Trustees of the University of Illinois.         *
!   All rights reserved.                                                      *
!                                                                             *
!   This routine is part of HDF5.  The full HDF5 copyright notice, including     *
!   terms governing use, modification, and redistribution, is contained in    *
!   the COPYING file, which can be found at the root of the source code       *
!   distribution tree, or in https://support.hdfgroup.org/ftp/HDF5/releases.  *
!   If you do not have access to either file, you may request a copy from     *
!   help@hdfgroup.org.                                                        *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! The following example shows how to write and read to/from an existing dataset.
! It opens the file created in the previous example, obtains the dataset
! identifier, writes the data to the dataset in the file,
! then reads the dataset  to memory.
!
! This example is used in the HDF5 Tutorial.

! Initialize the dset_data array.

character(len=8), parameter :: filename = "dsetf.h5" ! File name
character(len=4), parameter :: dsetname = "dset"     ! Dataset name

integer(HID_T) :: file_id       ! File identifier
integer(HID_T) :: dset_id       ! Dataset identifier

integer :: error ! Error flag
integer :: i, j

integer, dimension(4,6) :: dset_data, data_out ! Data buffers
integer(HSIZE_T), dimension(2) :: data_dims

DO i = 1, 4
   DO j = 1, 6
      dset_data(i,j) = (i-1)*6 + j
   END DO
END DO


! Initialize FORTRAN interface.
call h5open_f(error)

! Open an existing file.
call h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)

! Open an existing dataset.
call h5dopen_f(file_id, dsetname, dset_id, error)

! Write the dataset.

data_dims(1) = 4
data_dims(2) = 6
call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, dset_data, data_dims, error)

! Read the dataset.
call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data_out, data_dims, error)

! Close the dataset.
call h5dclose_f(dset_id, error)

! Close the file.
call h5fclose_f(file_id, error)

! Close FORTRAN interface.
call h5close_f(error)

end subroutine H5_RDWT


end module HDF5_utilities_mod

