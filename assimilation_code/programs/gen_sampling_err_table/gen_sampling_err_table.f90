! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Correct covariances for fixed ensemble sizes.
!>
!> Ref: Anderson, J., 2012: 
!> Localization and Sampling Error Correction in Ensemble Kalman Filter Data Assimilation.
!> Mon. Wea. Rev., 140, 2359-2371, doi: 10.1175/MWR-D-11-00013.1. 
!>
!> this version is a sparse array - the ens_sizes array lists the ensemble
!> sizes that have been computed, using the unlimited dimension.  to generate
!> a new ensemble size, it can be added at the end of the existing arrays.
!> the read code will search the unlimited dimension until it finds a matching
!> ensemble size.
!>
!> the initial table has the 40 sizes that the previous release of DART
!> came with.  to generate other ensemble sizes run this program with
!> the desired ensemble sizes listed in the namelist.
!>
!> the random number generator is seeded with the ensemble size
!> so repeated runs of the program will generate the same tables.
!>
!> the output file is named 'sampling_error_correction_table.nc'
!> and if one already exists it will be appended to instead of
!> created from scratch.


program gen_sampling_err_table

use types_mod,      only : r8, MISSING_I
use utilities_mod,  only : error_handler, E_ERR, file_exist,  &
                           initialize_utilities, finalize_utilities, &
                           find_namelist_in_file, check_namelist_read, &
                           do_nml_file, do_nml_term, nmlfileunit, E_MSG
use  netcdf_utilities_mod, only : nc_check
use random_seq_mod, only : random_seq_type, init_random_seq, twod_gaussians

use netcdf

implicit none

character(len=*), parameter :: source = 'gen_sampling_err_table.f90'

integer, parameter :: num_times   = 1
integer, parameter :: num_samples = 100000000 ! large number for statistical rigor

character(len=128) :: output_filename = 'sampling_error_correction_table.nc'
integer, parameter :: nentries = 200

integer, parameter :: MAX_LIST_LEN = 200
integer, parameter :: UNSET = -1


type (random_seq_type) :: ran_id

real(r8) :: alpha(nentries), true_correl_mean(nentries)
integer  :: bin_count(0:nentries+1)

integer, allocatable :: index_array(:)

integer  :: ncid, num_ens, add_ens, esize, current_count
integer  :: iunit, io

character(len=512) :: msgstring

!----------------------------------------------------------------
! namelist item(s)


integer :: ens_sizes(MAX_LIST_LEN) = UNSET
logical :: debug = .false.

namelist /gen_sampling_error_table_nml/ ens_sizes, debug


!----------------------------------------------------------------
! start of executable code

call initialize_utilities('gen_sampling_err_table')

! Read the namelist entry
call find_namelist_in_file("input.nml", "gen_sampling_error_table_nml", iunit)
read(iunit, nml = gen_sampling_error_table_nml, iostat = io)
call check_namelist_read(iunit, io, "gen_sampling_error_table_nml")

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=gen_sampling_error_table_nml)
if (do_nml_term()) write(     *     , nml=gen_sampling_error_table_nml)

! blank line to make it easier to see messages 
call error_handler(E_MSG, '', '')

! check that the ensemble sizes are all good - 3 is the minimum size accepted
! and there is no max limit.  returns the number of values found
num_ens = valid_entries(ens_sizes, 3, MISSING_I, 'namelist variable "ens_size"')
if (num_ens < 1) then
   call error_handler(E_ERR, 'gen_sampling_err_table:', 'no valid ensemble sizes specified', &
             source, text2='check values of namelist item "ens_sizes"')
endif

if (file_exist(output_filename)) then
   ncid = setup_to_append_output_file(nentries, num_samples, current_count)
else
   ncid = create_output_file(nentries, num_samples, current_count)
endif

! add the new values to an existing list, if present, remove values that 
! already exist or are repeated in the list.  add_ens returns with the actual
! number of new ensemble sizes which need to be computed and the ens_sizes
! array has been updated to remove duplicates.

! save original value and tell user if we're avoiding duplicates or
! existing values in file.

call merge_lists(ncid, current_count, index_array, num_ens, ens_sizes, add_ens)

if (add_ens <= 0) then

   call error_handler(E_MSG, 'gen_sampling_err_table:', &
             'no new ensemble sizes specified, nothing to do.', &
             source, text2='all sizes specified already exist in output file')

else
   if (num_ens /= add_ens) then
      write(msgstring, *) num_ens-add_ens, ' "ens_sizes" entries ignored; either duplicate values'
      call error_handler(E_MSG, 'gen_sampling_err_table:', msgstring, source, &
             text2='in namelist, or ensemble size already exists in output file')
      call error_handler(E_MSG, '', '')  ! blank line
   endif


   do esize=current_count+1, current_count+add_ens
      ! seed the generator with the ensemble size so the results are
      ! completely reproducible.
      call init_random_seq(ran_id, index_array(esize))
   
      call compute_table(index_array(esize), nentries, bin_count, true_correl_mean, alpha)
      call addto_output_file(ncid, esize, bin_count(1:nentries), true_correl_mean, &
                             alpha, index_array(esize))
   
      call sync_output_file(ncid)
   
   enddo
endif

call close_output_file(ncid)

call finalize_utilities()

! end of main program
!----------------------------------------------------------------

contains

!----------------------------------------------------------------
!> main computational routine

subroutine compute_table(this_size, nentries, bin_count, true_correl_mean, alpha)

integer,  intent(in)  :: this_size
integer,  intent(in)  :: nentries
integer,  intent(out) :: bin_count(0:nentries+1)
real(r8), intent(out) :: true_correl_mean(nentries)
real(r8), intent(out) :: alpha(nentries)

real(r8), allocatable  :: pairs(:,:), temp(:)

real(r8) :: zero_2(2) = 0.0_r8, cov(2, 2)
real(r8) :: t_correl, sample_correl, beta
real(r8) :: s_mean(2), s_var(2), reg_mean, reg_sd, t_sd1, t_sd2
real(r8) :: tcorrel_sum(nentries), reg_sum(nentries), reg_2_sum(nentries)

integer  :: i, j, k, bin_num

character(len=64) :: context = 'bin, count, mean, alpha'

write(msgstring, *) 'computing sampling error correction table for ensemble size of ', this_size
call error_handler(E_MSG, '', msgstring)

allocate(pairs(2, this_size), temp(this_size))

bin_count   = 0
tcorrel_sum = 0.0_r8
reg_sum     = 0.0_r8
reg_2_sum   = 0.0_r8

! Uniformly distributed sequence of true correlations  
do j = 1, num_samples
   ! Assume that true correlation is uniform [-1,1]
   t_correl = -1.0_r8 + 2.0_r8 * ((j - 1) / (num_samples - 1.0_r8))
   if(j > num_samples) t_correl = 0.0_r8

   do k = 1, num_times

      t_sd1 = 1.0_r8
      t_sd2 = 1.0_r8

      ! Generate the covariance matrix for this correlation
      cov(1, 1) = t_sd1**2
      cov(2, 2) = t_sd2**2
      cov(1, 2) = t_correl * t_sd1 * t_sd2
      cov(2, 1) = cov(1, 2)
   
      ! Loop to generate an ensemble size sample from this correl
      ! Generate a random sample of size ens_size from something with this correlation
      do i = 1, this_size
         call twod_gaussians(ran_id, zero_2, cov, pairs(:, i))
      enddo
   
      ! Compute the sample correlation
      call comp_correl(pairs, this_size, sample_correl)

      ! Bin statistics for each 0.01 sample correlation bin; 
      !  start with just mean true correlation
      bin_num = (sample_correl + 1.0_r8) / 0.01_r8   + 1
      if (bin_num < 1) then
         bin_count(0) = bin_count(0) + 1
         cycle
      endif
      if (bin_num > nentries) then
         bin_count(nentries+1) = bin_count(nentries+1) + 1
         cycle
      endif

      ! Keep a sum of the sample_correl and squared to compute standard deviation
      tcorrel_sum(bin_num) = tcorrel_sum(bin_num) + t_correl

      !-----------------
      ! Also interested in finding out what the spurious variance reduction factor is
      ! First need to compute sample s.d. for the obs and unobs variable
      do i = 1, 2
         temp = pairs(i, :)
         call sample_mean_var(temp, this_size, s_mean(i), s_var(i))
      enddo

      !-----------------
      reg_sum(bin_num) = reg_sum(bin_num) + t_correl * sqrt(s_var(2) / s_var(1))
      reg_2_sum(bin_num) = reg_2_sum(bin_num) + (t_correl * sqrt(s_var(2) / s_var(1)))**2

      bin_count(bin_num) = bin_count(bin_num) + 1

   enddo
enddo

! print out just to stdout how many counts, if any, fell below bin 0
if (debug) then
   write(msgstring,*) 0, bin_count(0), 0.0_r8, 0.0_r8
   call error_handler(E_MSG, context, msgstring) 
endif

! always generate exactly 'nentries' entries in the output file
do i = 1, nentries
   ! must have at least 2 counts to compute a std dev
   if(bin_count(i) <= 1) then
      write(msgstring, *) 'Bin ', i, ' has ', bin_count(i), ' counts'
      call error_handler(E_ERR, 'compute_table', msgstring, &
         source, text2='All bins must have at least 2 counts')
   endif
   
   ! Compute the standard deviation of the true correlations
   true_correl_mean(i) = tcorrel_sum(i) / bin_count(i)
   reg_mean = reg_sum(i) / bin_count(i)
   reg_sd = sqrt((reg_2_sum(i) - bin_count(i) * reg_mean**2) / (bin_count(i) - 1))

   if(reg_sd <= 0.0_r8) then
      alpha(i) = 1.0_r8
   else
      !!!beta = reg_mean**2 / reg_sd**2
      ! Correct for bias in the standard deviation for very small ensembles, too 
      beta = reg_mean**2 / (reg_sd**2 * (1.0_r8 + 1.0_r8 / this_size))
      alpha(i) = beta / (1.0_r8 + beta)
   endif

   if (debug) then
      write(msgstring,*) i, bin_count(i), true_correl_mean(i), alpha(i)
      call error_handler(E_MSG, context, msgstring) 
   endif

! original code had ascii output in a text file, using this format:
!      write(iunit, 10)   i, bin_count(i), true_correl_mean(i), alpha(i)
! 10 format (I4,I9,2G25.14)

enddo

if (debug) then
   write(msgstring,*) nentries+1, bin_count(nentries+1), 0.0_r8, 0.0_r8
   call error_handler(E_MSG, context, msgstring) 
endif

deallocate(pairs, temp)


end subroutine

!----------------------------------------------------------------
!> stat routines

subroutine comp_correl(ens, n, correl)

integer,  intent(in)  :: n
real(r8), intent(in)  :: ens(2, n)
real(r8), intent(out) :: correl

real(r8) :: sum_x, sum_y, sum_xy, sum_x2, sum_y2

sum_x  = sum(ens(2, :))
sum_y  = sum(ens(1, :))
sum_xy = sum(ens(2, :) * ens(1, :))
sum_x2 = sum(ens(2, :) * ens(2, :))

! Computation of correlation
sum_y2 = sum(ens(1, :) * ens(1, :))

correl = (n * sum_xy - sum_x * sum_y) / &
   sqrt((n * sum_x2 - sum_x**2) * (n * sum_y2 - sum_y**2))

end subroutine comp_correl

!----------------------------------------------------------------
!>

subroutine sample_mean_var(x, n, mean, var)

integer,  intent(in)  :: n
real(r8), intent(in)  :: x(n)
real(r8), intent(out) :: mean, var

real(r8) :: sx, s_x2

sx   = sum(x)
s_x2 = sum(x * x)
mean = sx / n
var  = (s_x2 - sx**2 / n) / (n - 1)

end subroutine sample_mean_var

!----------------------------------------------------------------
! netcdf i/o routines

!----------------------------------------------------------------
!> create a new netcdf file and add some global attributes

function create_output_file(nentries, num_samples, unlimlen)

integer, intent(in)  :: nentries
integer, intent(in)  :: num_samples
integer, intent(out) :: unlimlen

integer :: create_output_file
integer :: rc, fid

rc = nf90_create(output_filename, NF90_CLOBBER, fid)
call nc_check(rc, 'create_output_file', 'creating "'//trim(output_filename)//'"')

call set_global_int_att(fid, 'num_samples', num_samples)

call set_global_char_att(fid, 'title', 'Sampling Error Corrections for fixed ensemble sizes.' )

msgstring = 'Anderson, J., 2012: Localization and Sampling Error Correction &
     &in Ensemble Kalman Filter Data Assimilation. Mon. Wea. Rev., 140, 2359-2371, &
     &doi: 10.1175/MWR-D-11-00013.1. '

call set_global_char_att(fid, 'reference', msgstring)

msgstring = '$Id$'
call set_global_char_att(fid, 'version', msgstring)


call setup_output_file(fid)

! unlimited dim initially empty
unlimlen = 0

create_output_file = fid

end function create_output_file

!----------------------------------------------------------------
!> open an existing file and prepare to append to it

function setup_to_append_output_file(nentries, num_samples, unlimlen)

integer, intent(in)  :: nentries
integer, intent(in)  :: num_samples
integer, intent(out) :: unlimlen

integer :: setup_to_append_output_file

integer :: rc, fid, nbins, nsamp

! for netcdf, write means update 
rc = nf90_open(output_filename, NF90_WRITE, fid)
call nc_check(rc, 'setup_to_append_output_file', 'opening "'//trim(output_filename)//'"')

! before we start to update this file, make sure nbins and nentries
! match the existing values in the file.
  
call get_sec_dim_info(fid, 'bins', l1=nbins)
if (nbins /= nentries) then
   write(msgstring, *) 'existing file has ', nbins, ' bins, the program has ', nentries, ' bins.'
   call error_handler(E_ERR, 'setup_to_append_output_file', &
             'existing file used a different bin size', source, text2=msgstring)
endif

! also make sure num_samples matches
  
rc = nf90_get_att(fid, NF90_GLOBAL, 'num_samples', nsamp)
call nc_check(rc, 'setup_to_append_output_file', 'getting global attribute "num_samples"')
if (nsamp /= num_samples) then
   write(msgstring, *) 'existing file uses ', nsamp, ' samples, the program has ', &
                       num_samples, ' samples.'
   call error_handler(E_ERR, 'setup_to_append_output_file', &
             'existing file used a different number of samples', source, text2=msgstring)
endif

! get the current size of the unlimited dimension 
call get_sec_dim_info(fid, 'ens_sizes', l1=unlimlen)

setup_to_append_output_file = fid

end function setup_to_append_output_file

!----------------------------------------------------------------
!> define 2 dims and 4 arrays in output file

subroutine setup_output_file(ncid)

integer, intent(in)  :: ncid

integer :: rc
integer :: nbinsDimID, nensDimID

call setup_sec_dim(ncid, 'bins', nentries, nbinsDimID)
call setup_sec_unlimdim(ncid, 'ens_sizes', nensDimID)

call setup_sec_data_int (ncid, 'count', nbinsDimID, nensDimID)
call set_var_char_att(ncid, 'count', 'description','number of samples in each bin')

call setup_sec_data_real(ncid, 'true_corr_mean', nbinsDimID, nensDimID)
call setup_sec_data_real(ncid, 'alpha',          nbinsDimID, nensDimID)
call set_var_char_att(ncid, 'alpha', 'description','sampling error correction factors')

call setup_sec_data_int1d (ncid, 'ens_sizes', nensDimID)
call set_var_char_att(ncid, 'ens_sizes', 'description','ensemble size used for calculation')

rc = nf90_enddef(ncid)
call nc_check(rc, 'setup_output_file', 'enddef')

end subroutine setup_output_file

!----------------------------------------------------------------
!> write 3 arrays and an integer to output file

subroutine addto_output_file(ncid, slot, count_data, corrmean_data, alpha_data, index_num)

integer,  intent(in)  :: ncid
integer,  intent(in)  :: slot
integer,  intent(in)  :: count_data(:)
real(r8), intent(in)  :: corrmean_data(:)
real(r8), intent(in)  :: alpha_data(:)
integer,  intent(in)  :: index_num

call write_sec_data_int  (ncid, slot, 'count',          count_data)
call write_sec_data_real (ncid, slot, 'true_corr_mean', corrmean_data)
call write_sec_data_real (ncid, slot, 'alpha',          alpha_data)
call write_sec_data_int1d(ncid, slot, 'ens_sizes',      index_num)

end subroutine addto_output_file

!----------------------------------------------------------------
!>

subroutine sync_output_file(ncid)

integer, intent(in) :: ncid

integer :: rc

rc = nf90_sync(ncid)
call nc_check(rc, 'sync_output_file', 'syncing '//trim(output_filename))

end subroutine sync_output_file

!----------------------------------------------------------------
!>

subroutine close_output_file(ncid)

integer, intent(in) :: ncid

integer :: rc

rc = nf90_close(ncid)
call nc_check(rc, 'close_output_file', 'closing '//trim(output_filename))

end subroutine close_output_file

!----------------------------------------------------------------
!> define a dimension with name and a length and return the id

subroutine setup_sec_dim(ncid, c1, n1, id1)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: c1
integer,          intent(in)  :: n1
integer,          intent(out) :: id1

integer :: rc

rc = nf90_def_dim(ncid, name=c1, len=n1, dimid=id1)
call nc_check(rc, 'setup_sec_dim', 'adding dimension "'//trim(c1)//'"')

end subroutine setup_sec_dim

!----------------------------------------------------------------
!> define an unlimited dimension and return the id

subroutine setup_sec_unlimdim(ncid, c1, id1)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: c1
integer,          intent(out) :: id1

integer :: rc

rc = nf90_def_dim(ncid, name=c1, len=NF90_UNLIMITED, dimid=id1)
call nc_check(rc, 'setup_sec_unlimdim', 'adding dimension "'//trim(c1)//'"')

end subroutine setup_sec_unlimdim

!----------------------------------------------------------------
!> given a netCDF file ID, variable name, and 2 dimension IDs,
!> create a 2d integer variable 

subroutine setup_sec_data_int(ncid, c1, d1, d2)

integer,                    intent(in)  :: ncid
character(len=*),           intent(in)  :: c1
integer,                    intent(in)  :: d1, d2

integer :: rc, id1

rc = nf90_def_var(ncid, name=c1, xtype=nf90_int, dimids=(/ d1, d2 /), varid=id1)
call nc_check(rc, 'setup_sec_data_int', 'defining variable "'//trim(c1)//'"')

end subroutine setup_sec_data_int

!----------------------------------------------------------------
!> given a netCDF file ID, variable name, and a dimension ID,
!> create an integer variable

subroutine setup_sec_data_int1d(ncid, c1, d1)

integer,                    intent(in)  :: ncid
character(len=*),           intent(in)  :: c1
integer,                    intent(in)  :: d1

integer :: rc, id1

rc = nf90_def_var(ncid, name=c1, xtype=nf90_int, dimids=(/ d1 /), varid=id1)
call nc_check(rc, 'setup_sec_data_int1d', 'defining variable "'//trim(c1)//'"')

end subroutine setup_sec_data_int1d

!----------------------------------------------------------------
!> given a netCDF file ID, variable name, and 2 dimension IDs,
!> create a 2D real variable

subroutine setup_sec_data_real(ncid, c1, d1, d2)

integer,                    intent(in)  :: ncid
character(len=*),           intent(in)  :: c1
integer,                    intent(in)  :: d1, d2

integer :: rc, id1

rc = nf90_def_var(ncid, name=c1, xtype=nf90_double, dimids=(/ d1, d2 /), varid=id1)
call nc_check(rc, 'setup_sec_data_real', 'defining variable "'//trim(c1)//'"')

end subroutine setup_sec_data_real

!----------------------------------------------------------------
!>

subroutine write_sec_data_int(ncid, col, c1, a1)

integer,          intent(in) :: ncid
integer,          intent(in) :: col
character(len=*), intent(in) :: c1
integer,          intent(in) :: a1(:)

integer :: rc, id1

rc = nf90_inq_varid(ncid, c1, id1)
call nc_check(rc, 'write_sec_data_int', 'querying variable "'//trim(c1)//'"')

rc = nf90_put_var(ncid, id1, a1, start=(/ 1, col /), count=(/ size(a1), 1 /) )
call nc_check(rc, 'write_sec_data_int', 'writing variable "'//trim(c1)//'"')

end subroutine write_sec_data_int

!----------------------------------------------------------------

subroutine write_sec_data_int1d(ncid, col, c1, a1)

integer,          intent(in) :: ncid
integer,          intent(in) :: col
character(len=*), intent(in) :: c1
integer,          intent(in) :: a1

integer :: rc, id1

rc = nf90_inq_varid(ncid, c1, id1)
call nc_check(rc, 'write_sec_data_int1d', 'querying variable "'//trim(c1)//'"')

rc = nf90_put_var(ncid, id1, a1, start=(/ col /))
call nc_check(rc, 'write_sec_data_int1d', 'writing variable "'//trim(c1)//'"')

end subroutine write_sec_data_int1d

!----------------------------------------------------------------

subroutine write_sec_data_real(ncid, col, c1, a1) 

integer,          intent(in) :: ncid
integer,          intent(in) :: col
character(len=*), intent(in) :: c1
real(r8),         intent(in) :: a1(:)

integer :: rc, id1

rc = nf90_inq_varid(ncid, c1, id1)
call nc_check(rc, 'write_sec_data_real', 'querying variable "'//trim(c1)//'"')

rc = nf90_put_var(ncid, id1, a1, start=(/ 1, col /), count=(/ size(a1), 1 /) )
call nc_check(rc, 'write_sec_data_real', 'writing variable "'//trim(c1)//'"')

end subroutine write_sec_data_real

!----------------------------------------------------------------

subroutine read_sec_data_int(ncid, col, c1, a1)

integer,          intent(in)  :: ncid
integer,          intent(in)  :: col
character(len=*), intent(in)  :: c1
integer,          intent(out) :: a1(:)

integer :: rc, id1

rc = nf90_inq_varid(ncid, c1, id1)
call nc_check(rc, 'read_sec_data_int', 'querying variable "'//trim(c1)//'"')

rc = nf90_get_var(ncid, id1, a1, start=(/ 1, col /), count=(/ size(a1), 1 /) )
call nc_check(rc, 'read_sec_data_int', 'reading variable "'//trim(c1)//'"')

end subroutine read_sec_data_int

!----------------------------------------------------------------
!>

subroutine read_sec_data_real(ncid, col, c1, a1) 

integer,          intent(in)  :: ncid
integer,          intent(in)  :: col
character(len=*), intent(in)  :: c1
real(r8),         intent(out) :: a1(:)

integer :: rc, id1

rc = nf90_inq_varid(ncid, c1, id1)
call nc_check(rc, 'read_sec_data_real', 'querying variable "'//trim(c1)//'"')

rc = nf90_get_var(ncid, id1, a1, start=(/ 1, col /), count=(/ size(a1), 1 /) )
call nc_check(rc, 'read_sec_data_real', 'reading variable "'//trim(c1)//'"')

end subroutine read_sec_data_real

!----------------------------------------------------------------
!> retrieve either the length of a dimension, the dim id, or both

subroutine get_sec_dim_info(ncid, c1, l1, id1)

integer,           intent(in)  :: ncid
character(len=*),  intent(in)  :: c1
integer, optional, intent(out) :: l1
integer, optional, intent(out) :: id1

integer :: rc, ll1, lid1

rc = nf90_inq_dimid(ncid, c1, lid1)
call nc_check(rc, 'get_sec_dim_info', 'inq_dimid "'//trim(c1)//'"')

rc = nf90_inquire_dimension(ncid, lid1, len=ll1)
call nc_check(rc, 'get_sec_dim_info', 'inquire_dimension "'//trim(c1)//'"')

if (present(l1)) l1 = ll1
if (present(id1)) id1 = lid1

end subroutine get_sec_dim_info

!----------------------------------------------------------------
!> set a character attribute on a variable

subroutine set_var_char_att(fid, varname, attname, attvalue)

integer,          intent(in) :: fid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: attname
character(len=*), intent(in) :: attvalue

integer :: rc, id1

rc = nf90_inq_varid(fid, varname, id1)
call nc_check(rc, 'set_var_char_att', 'inquiring variable id "'//trim(varname)//'"')

rc = nf90_put_att(fid, id1, attname, attvalue)
call nc_check(rc, 'set_var_char_att', 'adding attribute "'//trim(attname)//'" to variable "'//&
                                      trim(varname)//'"')

end subroutine set_var_char_att

!----------------------------------------------------------------
!> set a global character attribute 

subroutine set_global_char_att(fid, attname, attvalue)

integer,          intent(in) :: fid
character(len=*), intent(in) :: attname
character(len=*), intent(in) :: attvalue

integer :: rc

rc = nf90_put_att(fid, NF90_GLOBAL, attname, attvalue)
call nc_check(rc, 'set_global_char_att', 'adding global attribute "'//trim(attname)//'"')

end subroutine set_global_char_att

!----------------------------------------------------------------
!> set a global integer attribute 

subroutine set_global_int_att(fid, attname, attvalue)

integer,          intent(in) :: fid
character(len=*), intent(in) :: attname
integer,          intent(in) :: attvalue

integer :: rc

rc = nf90_put_att(fid, NF90_GLOBAL, attname, attvalue)
call nc_check(rc, 'set_global_int_att', 'adding global attribute "'//trim(attname)//'"')

end subroutine set_global_int_att

!----------------------------------------------------------------
! misc utility routines

!----------------------------------------------------------------
!> figure out the list length.  -1 means no entries; a minimum
!> valid value can be specfied and is a fatal error if an entry
!> is below the threshold.
!>@todo FIXME this belongs in the utilities module

function valid_entries(list, min_valid, max_valid, where)

integer,          intent(in) :: list(:)
integer,          intent(in) :: min_valid
integer,          intent(in) :: max_valid
character(len=*), intent(in) :: where
integer :: valid_entries

integer :: i, val

val = -1
do i=1, MAX_LIST_LEN
   if (list(i) == UNSET) exit

   if (min_valid /= MISSING_I .and. list(i) < min_valid) then
      write(msgstring, *) 'minimum valid value is ', min_valid
      call error_handler(E_ERR, 'valid_entries', &
                 'illegal value found in '//trim(where), source, text2=msgstring)
   endif

   if (max_valid /= MISSING_I .and. list(i) > max_valid) then
      write(msgstring, *) 'maximum valid value is ', max_valid
      call error_handler(E_ERR, 'valid_entries', &
                 'illegal value found in '//trim(where), source, text2=msgstring)
   endif

   val = i
enddo

valid_entries = val

end function valid_entries

!----------------------------------------------------------------
!> if this routine returns, index_array is allocated and filled
!> with all ensemble sizes.  num_add_ens is the number of new
!> ensemble sizes that must be computed.  duplicates and existing
!> values are removed from the list.

subroutine merge_lists(fid, current_size, index_array, num_new_ens, ens_sizes, num_add_ens)

integer,              intent(in)    :: fid
integer,              intent(in)    :: current_size
integer, allocatable, intent(out)   :: index_array(:)
integer,              intent(in)    :: num_new_ens
integer,              intent(inout) :: ens_sizes(:)
integer,              intent(out)   :: num_add_ens

integer :: num_ens

! allocate list large enough to hold the max possible output size
! and read any existing entries.
allocate(index_array(current_size+num_new_ens))

if (current_size > 0) &
   call read_sec_data_int(fid, 1, 'ens_sizes', index_array(1:current_size))

! add new values to end of list and look for dups. 

if (num_new_ens > 0) &
   index_array(current_size+1:current_size+num_new_ens) = ens_sizes(1:num_new_ens)

! total list length
num_ens = current_size + num_new_ens

! this value will get updated if there are duplicates
! in the list or items already in the existing file.
call remove_dups(num_ens, index_array)

! return the count of new items to be computed
num_add_ens = num_ens - current_size

end subroutine merge_lists

!----------------------------------------------------------------
!> given a list and a count, remove any replicated values from the
!> list and update the count.   
!>@todo FIXME this belongs in a utilities module

subroutine remove_dups(num_vals, val_list)

integer, intent(inout) :: num_vals
integer, intent(inout) :: val_list(:)

integer :: i, j, k, newlen, next_slot
integer, allocatable :: new_list(:)

! empty list, nothing to do
if (num_vals == 0) return

! allocate enough temp space to fit a list with no duplicates
allocate(new_list(num_vals))

! iterate the input list and copy unique values to the temp list

next_slot = 0
outerloop: do i = 1, num_vals
   do j = 1, next_slot
      if (val_list(i) == new_list(j)) cycle outerloop
   enddo
   next_slot = next_slot + 1
   new_list(next_slot) = val_list(i)
enddo outerloop

! return values back into original list.
! guarenteed to be same size or shorter.
num_vals = next_slot
val_list(1:next_slot) = new_list(1:next_slot)

deallocate(new_list)

end subroutine remove_dups

!----------------------------------------------------------------
!> debugging routine

subroutine dump_list(count, list, label)

integer,          intent(in) :: count
integer,          intent(in) :: list(:)
character(len=*), intent(in) :: label

print *, 'dump_list: '//trim(label)
print *, ' len=', count
print *, ' values=', list(1:count)

end subroutine dump_list

!----------------------------------------------------------------

end program gen_sampling_err_table

