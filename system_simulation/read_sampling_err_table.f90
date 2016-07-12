! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> Correct covariances for fixed ensemble sizes.
!> See Anderson, J. L., 2011: Localization and Sampling Error Correction
!>   in Ensemble Kalman Filter Data Assimilation. 
!> Submitted for publication, Jan 2011.  Contact author.

!> read the entry for a single ensemble_size and print out the
!> two values, true_correl_mean and alpha.   mostly as a test for
!> accuracy, but also example code for assim_tools to use.

program read_sampling_err_table

use types_mod,      only : r8
use utilities_mod,  only : error_handler, E_ERR, nc_check,      &
                           initialize_utilities, finalize_utilities

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"


integer :: max_ens_size, nentries

character(len=128) :: input_filename = 'sampling_error_correction_table.nc'

real(r8), allocatable :: true_correl_mean(:), alpha(:)

integer  :: i, j, k
integer  :: iunit, io, ncid, ens_size, requested_ens_size

character(len=512) :: errstring, msgstring

!
! start of executable code
!

call initialize_utilities('read_sampling_err_table')

print *, 'Enter ensemble size: '
read *, requested_ens_size

ncid = open_input_file()

call setup_input_file(ncid, nentries, max_ens_size)

if (requested_ens_size < 2 .or. requested_ens_size > max_ens_size) then
   call error_handler(E_ERR, 'read_sampling_err_table:', 'bad ensemble size requested', &
                      source, revision, revdate)
endif

allocate(true_correl_mean(nentries), alpha(nentries))

call read_input_file(ncid, requested_ens_size, true_correl_mean, alpha)

print *, 'true correlation means, and alphas: '
do i=1, nentries
   print *, i, true_correl_mean(i), alpha(i)
enddo

call close_input_file(ncid)

call finalize_utilities()

! end of main program

contains

!----------------------------------------------------------------

!----------------------------------------------------------------
! main netcdf i/o routines
!----------------------------------------------------------------

!----------------------------------------------------------------

function open_input_file()

integer :: open_input_file

integer :: rc, ncid

rc = nf90_open(input_filename, NF90_NOWRITE, ncid)
call nc_check(rc, 'open_input_file', 'creating '//trim(input_filename))

rc = nf90_get_att(ncid, NF90_GLOBAL, 'ensemble_sizes', msgstring)
call nc_check(rc, 'open_input_file', 'getting global attributes')
print *, trim(msgstring)

rc = nf90_get_att(ncid, NF90_GLOBAL, 'bin_info', msgstring)
call nc_check(rc, 'open_input_file', 'getting global attributes')
print *, trim(msgstring)

open_input_file = ncid

end function open_input_file

!----------------------------------------------------------------

! get the 2 dims - number of entries for any given ensemble size,
! and valid ranges of ensemble sizes.

subroutine setup_input_file(ncid, nbins, nens)

integer, intent(in)  :: ncid
integer, intent(out) :: nbins, nens

integer :: rc
integer :: nbinsDimID, nensDimID

call get_sec_dim(ncid, 'bins', nbins)
call get_sec_dim(ncid, 'ens',  nens)

end subroutine setup_input_file

!----------------------------------------------------------------

! read 2 arrays

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

subroutine close_input_file(ncid)

integer, intent(in) :: ncid

integer :: rc

rc = nf90_close(ncid)
call nc_check(rc, 'close_input_file', 'closing '//trim(input_filename))

end subroutine close_input_file

!----------------------------------------------------------------

!----------------------------------------------------------------
! helper routines for above code
!----------------------------------------------------------------

!----------------------------------------------------------------

! retrieve a dimension

subroutine get_sec_dim(ncid, c1, n1)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: c1
integer,          intent(out) :: n1

integer :: rc, id1

rc = nf90_inq_dimid(ncid, c1, id1)
call nc_check(rc, 'get_sec_dim', 'querying dimension '//trim(c1))

rc = nf90_inquire_dimension(ncid, id1, len=n1)
call nc_check(rc, 'get_sec_dim', 'querying dimension '//trim(c1))

end subroutine get_sec_dim

!----------------------------------------------------------------

! given a name, return id for a 2d real array

subroutine query_sec_data(ncid, c1, id1)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: c1
integer,          intent(out) :: id1

integer :: rc

rc = nf90_inq_varid(ncid, name=c1, varid=id1)
call nc_check(rc, 'query_sec_data', 'querying variable '//trim(c1))

end subroutine query_sec_data

!----------------------------------------------------------------

subroutine read_sec_data_int(ncid, col, c1, id1, a1)

integer,          intent(in)  :: ncid
integer,          intent(in)  :: col
character(len=*), intent(in)  :: c1
integer,          intent(in)  :: id1
integer,          intent(out) :: a1(:)

integer :: rc

rc = nf90_get_var(ncid, id1, a1, start=(/ 1, col /), count=(/ size(a1), 1 /) )
call nc_check(rc, 'read_sec_data', 'reading variable "'//trim(c1)//'"')

end subroutine read_sec_data_int

!----------------------------------------------------------------

subroutine read_sec_data_real(ncid, col, c1, id1, a1) 
integer,          intent(in)  :: ncid
integer,          intent(in)  :: col
character(len=*), intent(in)  :: c1
integer,          intent(in)  :: id1
real(r8),         intent(out) :: a1(:)

integer :: rc

rc = nf90_get_var(ncid, id1, a1, start=(/ 1, col /), count=(/ size(a1), 1 /) )
call nc_check(rc, 'read_sec_data', 'reading variable "'//trim(c1)//'"')

end subroutine read_sec_data_real

!----------------------------------------------------------------


end program read_sampling_err_table

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
