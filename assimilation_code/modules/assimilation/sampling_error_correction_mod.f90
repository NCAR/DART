! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Correct covariances for fixed ensemble sizes.
!> Ref: Anderson, J., 2012: 
!> Localization and Sampling Error Correction in Ensemble Kalman Filter Data Assimilation.
!> Mon. Wea. Rev., 140, 2359-2371, doi: 10.1175/MWR-D-11-00013.1. 

!> query the needed table sizes, and read in the values for any
!> given ensemble size.  the two arrays of values returned are
!> the true_correl_mean and alpha. 

module sampling_error_correction_mod

use types_mod,           only : r8
use utilities_mod,       only : error_handler, E_ERR
use netcdf_utilities_mod, only : nc_check

use netcdf

implicit none
private

public :: get_sampling_error_table_size, &
          read_sampling_error_correction

character(len=*), parameter :: source = 'sampling_error_correction_mod.f90'

! Using hardcoded filename for ease of scripting.
! and for now, say where the default location in the dart distribution tree is
! since it's so obscure.
character(len=128) :: input_filename = 'sampling_error_correction_table.nc'
character(len=128) :: default_path = ' "assimilation_code/programs/gen_sampling_err_table/work"'

! module globals - nentries is the number of values per ensemble size,
! nens is how many different ensemble sizes this file contains.

logical :: module_initialized = .false.
integer :: nentries = -1
integer :: nens = -1

character(len=512) :: msgstring, msgstring1

contains

!----------------------------------------------------------------
!>

subroutine init_sampling_error_correction()

integer :: ncid

if (module_initialized) return

ncid = open_input_file(input_filename)

call read_input_info(ncid, nentries, nens)

call close_input_file(ncid, input_filename)

module_initialized = .true.

end subroutine init_sampling_error_correction

!----------------------------------------------------------------
!>

function get_sampling_error_table_size()

integer :: get_sampling_error_table_size

if (.not. module_initialized) call init_sampling_error_correction()

get_sampling_error_table_size = nentries

end function get_sampling_error_table_size

!----------------------------------------------------------------
!>

subroutine read_sampling_error_correction(requested_ens_size, true_correl_mean, alpha)

integer,  intent(in) :: requested_ens_size
real(r8), intent(out) :: true_correl_mean(:), alpha(:)

integer :: ncid, indx

if (.not. module_initialized) call init_sampling_error_correction()

ncid = open_input_file(input_filename)

indx = lookup_ens_index(ncid, nens, requested_ens_size)

if (indx < 0) then
   write(msgstring, *) 'file "'//trim(input_filename)//'" does not contain a entry for ensemble size ', &
                        requested_ens_size
   write(msgstring1, *) 'You can add one to the existing file with the "gen_sampling_err_table" program'
   call error_handler(E_ERR, 'read_sampling_error_correction:', 'unsupported ensemble size requested', &
                      source, text2=msgstring, text3=msgstring1)
endif

if (size(true_correl_mean(:)) /= nentries .or. &
    size(alpha(:)) /= nentries) then
   write(msgstring, *) 'one or both arrays "true_correl_mean" and "alpha" are not allocated correctly'
   write(msgstring1, *) 'they must be size ', nentries, ' but are ', size(true_correl_mean), ' and ', size(alpha)
   call error_handler(E_ERR, 'read_sampling_error_correction:', 'error in output array size', &
                      source, text2=msgstring, text3=msgstring1)
endif

call read_input_file(ncid, indx, true_correl_mean, alpha)

call close_input_file(ncid, input_filename)

end subroutine read_sampling_error_correction

!----------------------------------------------------------------
! support routines below
!----------------------------------------------------------------

function open_input_file(input_filename)

character(len=*), intent(in) :: input_filename
integer :: open_input_file

integer :: rc, ncid

rc = nf90_open(input_filename, NF90_NOWRITE, ncid)
if (rc /= nf90_noerr) then
   msgstring  = 'File "'//trim(input_filename)//'" not found in the current directory.'
   msgstring1 = 'This file can be copied from this location in the DART distribution: '
   call error_handler(E_ERR, 'read_sampling_error_correction:', msgstring, &
                      source, text2=msgstring1, text3=default_path)
endif

open_input_file = ncid

end function open_input_file

!----------------------------------------------------------------
!> get the 2 dims - number of entries for any given ensemble size,
!> and number of supported ensemble sizes.

subroutine read_input_info(ncid, nbins, nens)

integer, intent(in)  :: ncid
integer, intent(out) :: nbins
integer, intent(out) :: nens

call get_sec_dim(ncid, 'bins', nbins)
call get_sec_dim(ncid, 'ens_sizes',  nens)

end subroutine read_input_info

!----------------------------------------------------------------
!>

function lookup_ens_index(ncid, num_ens, requested_ens_size)

integer, intent(in) :: ncid
integer, intent(in) :: num_ens
integer, intent(in) :: requested_ens_size

integer :: lookup_ens_index
integer :: i, indx, id
integer, allocatable :: index_array(:)

allocate(index_array(num_ens))

call query_sec_data(ncid, 'ens_sizes', id)
call read_sec_data_int(ncid, 1, 'ens_sizes', id, index_array)

indx = -1
do i=1, num_ens
   if (index_array(i) == requested_ens_size) then
      indx = i
      exit
   endif
enddo

lookup_ens_index = indx
deallocate(index_array)

end function lookup_ens_index

!----------------------------------------------------------------
!> read the true correlation and correction arrays 

subroutine read_input_file(ncid, col, a1, a2)

integer,          intent(in)  :: ncid
integer,          intent(in)  :: col
real(r8),         intent(out) :: a1(:)
real(r8),         intent(out) :: a2(:)

integer :: id1, id2
character(len=64) :: c1, c2

c1 = 'true_corr_mean'
c2 = 'alpha'

call query_sec_data(ncid, c1, id1)
call query_sec_data(ncid, c2, id2)

call read_sec_data_real(ncid, col, c1, id1, a1)
call read_sec_data_real(ncid, col, c2, id2, a2)

end subroutine read_input_file

!----------------------------------------------------------------
!>

subroutine close_input_file(ncid, input_filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: input_filename

integer :: rc

rc = nf90_close(ncid)
call nc_check(rc, 'close_input_file', 'closing "'//trim(input_filename)//'"')

end subroutine close_input_file

!----------------------------------------------------------------
!> retrieve dimension for sampling error correction

subroutine get_sec_dim(ncid, c1, n1)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: c1
integer,          intent(out) :: n1

integer :: rc, id1

rc = nf90_inq_dimid(ncid, c1, id1)
call nc_check(rc, 'get_sec_dim', 'inq_dimid "'//trim(c1)//'"')

rc = nf90_inquire_dimension(ncid, id1, len=n1)
call nc_check(rc, 'get_sec_dim', 'inquire_dimension "'//trim(c1)//'"')

end subroutine get_sec_dim

!----------------------------------------------------------------
!> given a variable name, return variable id

subroutine query_sec_data(ncid, c1, id1)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: c1
integer,          intent(out) :: id1

integer :: rc

rc = nf90_inq_varid(ncid, name=c1, varid=id1)
call nc_check(rc, 'query_sec_data', 'querying variable "'//trim(c1)//'"')

end subroutine query_sec_data

!----------------------------------------------------------------
!>

subroutine read_sec_data_int(ncid, col, c1, id1, a1)

integer,          intent(in)  :: ncid
integer,          intent(in)  :: col
character(len=*), intent(in)  :: c1
integer,          intent(in)  :: id1
integer,          intent(out) :: a1(:)

integer :: rc

rc = nf90_get_var(ncid, id1, a1, start=(/ 1, col /), count=(/ size(a1), 1 /) )
call nc_check(rc, 'read_sec_data_int', 'reading variable "'//trim(c1)//'"')

end subroutine read_sec_data_int

!----------------------------------------------------------------
!>

subroutine read_sec_data_real(ncid, col, c1, id1, a1) 

integer,          intent(in)  :: ncid
integer,          intent(in)  :: col
character(len=*), intent(in)  :: c1
integer,          intent(in)  :: id1
real(r8),         intent(out) :: a1(:)

integer :: rc

rc = nf90_get_var(ncid, id1, a1, start=(/ 1, col /), count=(/ size(a1), 1 /) )
call nc_check(rc, 'read_sec_data_real', 'reading variable "'//trim(c1)//'"')

end subroutine read_sec_data_real

!----------------------------------------------------------------

end module sampling_error_correction_mod

