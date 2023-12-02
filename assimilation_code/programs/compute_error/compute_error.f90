! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Program to read true_state.nc and preassim.nc (or analysis.nc)
!> and print out the time-mean total error and spread.  The matlab macros 
!> do this with graphics, but if you want numerical output suitable for 
!> shell scripting/automation this should be easier to incorporate.
!> 
!> Unlike the Matlab routine, this program treats every value in the state
!> vector the same.  There is no distinction between variables of different
!> types (and possibly magnitudes), no weighting by area or volume based on
!> the underlying grid, and no computing values by model level.
!>
!> The output includes other text; to collect just the error and spread lines
!> grep standard output for 'Total'.  This program also uses the standard
!> DART utility routines that log both the output in dart_log.out and optionally
!> the namelist values in dart_log.nml.  If you are running this in a script
!> or repeatedly in a loop, you may want to delete or truncate those files
!> periodically.
!>
!> This program does compute any weighted errors.

program compute_error

use types_mod,     only : r8, metadatalength

use utilities_mod, only : initialize_utilities, &
                          error_handler, nmlfileunit, E_MSG, E_ERR,  &
                          find_namelist_in_file,                     &
                          check_namelist_read, finalize_utilities,   &
                          do_nml_file, do_nml_term
use    netcdf_utilities_mod, only : nc_check                                
use netcdf

implicit none

character(len=*), parameter :: source = 'compute_error.f90'

character(len = 512) :: message1, message2

integer :: iunit, model_size, io, ierr
integer :: ntimes, skip_start, i1, i2
integer :: true_copy, ens_mean_copy, ens_spread_copy
integer :: VarID

real(r8), allocatable :: truth(:,:), ens_mean(:,:), ens_spread(:,:), zeros(:,:)
real(r8), allocatable :: err_array(:), spread_array(:)

real(r8) :: total_error, total_spread

type fileinfo
  integer :: ncid
  integer :: numcopies
  integer :: model_size
  integer :: metadatacount
  integer :: ntimes
  integer :: sindex
  integer :: eindex
  real(r8), allocatable :: timevals(:)
  character(len=metadatalength), allocatable :: metadata(:)
  character(len=512) :: fname
end type fileinfo

type(fileinfo) :: t, d


! Namelist variables

character(len=256) :: diag_file_name  = 'preassim.nc'
character(len=256) :: truth_file_name = 'true_state.nc'
integer            :: skip_first_ntimes = 0

namelist /compute_error_nml/  &
   diag_file_name,            &
   truth_file_name,           &
   skip_first_ntimes

!----------------------------------------------------------------
! start of executable code
!----------------------------------------------------------------

call initialize_utilities('compute_error')


! Read the namelist entry and print it
call find_namelist_in_file("input.nml", "compute_error_nml", iunit)
read(iunit, nml = compute_error_nml, iostat = io)
call check_namelist_read(iunit, io, "compute_error_nml")

if (do_nml_file()) write(nmlfileunit, nml=compute_error_nml)
if (do_nml_term()) write(     *     , nml=compute_error_nml)

! record names for error messages, mostly.
t%fname = trim(truth_file_name)
d%fname = trim(diag_file_name)

! open truth, diag file
ierr = nf90_open(path=trim(truth_file_name), mode=nf90_nowrite, ncid=t%ncid)
call nc_check(ierr, 'compute_error', 'opening truth file')

ierr = nf90_open(path=trim(diag_file_name), mode=nf90_nowrite, ncid=d%ncid)
call nc_check(ierr, 'compute_error', 'opening diagnostic file')


! collect the dimension ids and lengths we are going to use from the truth file:  
!   copy, metadatalength, time and StateVariable.

call findlengths(t)
call findlengths(d)

if (t%ntimes <= 0) then
   write(message1, '(A)') 'file '//trim(t%fname)
   call error_handler(E_ERR, 'findlengths', 'time dimension must have at least 1 value', &
                      source, text2=message1)
endif
if (d%ntimes <= 0) then
   write(message1, '(A)') 'file '//trim(d%fname)
   call error_handler(E_ERR, 'findlengths', 'time dimension must have at least 1 value', &
                      source, text2=message1)
endif

call confirm_modelsizes_match(t, d)

! allocate space for data we are going to read, except for 'truth'.

allocate(t%metadata(t%numcopies), t%timevals(t%ntimes))
allocate(d%metadata(d%numcopies), d%timevals(d%ntimes))

call getchardata_nc(t%ncid, "CopyMetaData", t%metadata)
call getchardata_nc(d%ncid, "CopyMetaData", d%metadata)

call getrealdata1d_nc(t%ncid, "time", t%timevals)
call getrealdata1d_nc(d%ncid, "time", d%timevals)

call findtimesubset(t, d, ntimes)
model_size = t%model_size

! allocate space for data now that we know the time overlap.
allocate(truth(model_size, ntimes))
allocate(ens_mean(model_size, ntimes), ens_spread(model_size, ntimes))
allocate(zeros(model_size, ntimes), err_array(ntimes), spread_array(ntimes))

! findcopyindex() dies with an error if the copy name is not found.
call findcopyindex('CopyMetaData', t%metadata, t%numcopies, 'true state',      t%fname, true_copy)
call findcopyindex('CopyMetaData', d%metadata, d%numcopies, 'ensemble mean',   d%fname, ens_mean_copy)
call findcopyindex('CopyMetaData', d%metadata, d%numcopies, 'ensemble spread', d%fname, ens_spread_copy)


! get the data - truth, mean, spread - needed for the error calculation
ierr = nf90_inq_varid(t%ncid, "state", VarID)
call nc_check(ierr, 'compute_error', 'inq_varid state')

ierr = nf90_get_var(t%ncid, VarID, truth, &
                    start=(/ 1, true_copy, t%sindex /), &
                    count=(/ model_size, 1, t%ntimes /))
call nc_check(ierr, 'compute_error', 'get_var true state')

ierr = nf90_inq_varid(d%ncid, "state", VarID)
call nc_check(ierr, 'compute_error', 'inq_varid state')

ierr = nf90_get_var(d%ncid, VarID, ens_mean, &
                    start=(/ 1, ens_mean_copy, d%sindex /), &
                    count=(/ model_size, 1, d%ntimes /))
call nc_check(ierr, 'compute_error', 'get_var state mean')
ierr = nf90_get_var(d%ncid, VarID, ens_spread, &
                    start=(/ 1, ens_spread_copy, d%sindex /), &
                    count=(/ model_size, 1, d%ntimes /))
call nc_check(ierr, 'compute_error', 'get_var state spread')


! compute the error.  array is 'ntimes' long, each entry should be
!  sqrt of sum of (truth - diag) squared.
!  truth for spread is 0
! then avg error value is err array / ntimes
!
! weighted error would be:  sqrt of sum of ((truth - diag)^2 * weights)
! this codes does not do weights.

zeros(:,:) = 0.0_r8

call find_error(truth, ens_mean, err_array, total_error)
call find_error(zeros, ens_spread, spread_array, total_spread)

! output

call printout('Time-mean Ensemble Total  Error: ', total_error,  -1, -1.0_r8, &
              'Time-mean Ensemble Total Spread: ', total_spread, -1, -1.0_r8)
call error_handler(E_MSG,'', '')

call printout('First  Error Value, Time Index/Value: ', err_array(1),    1, t%timevals(t%sindex), &
              'First Spread Value, Time Index/Value: ', spread_array(1), 1, t%timevals(t%sindex))

i1 = minloc(err_array, 1)
i2 = minloc(spread_array, 1)
call printout('Min    Error Value, Time Index/Value: ', err_array(i1),    i1, t%timevals(t%sindex+i1), &
              'Min   Spread Value, Time Index/Value: ', spread_array(i2), i2, t%timevals(t%sindex+i2))

i1 = maxloc(err_array, 1)
i2 = maxloc(spread_array, 1)
call printout('Max    Error Value, Time Index/Value: ', err_array(i1),    i1, t%timevals(t%sindex+i1), &
              'Max   Spread Value, Time Index/Value: ', spread_array(i1), i2, t%timevals(t%sindex+i2))

call printout('Last   Error Value, Time Index/Value: ', err_array(ntimes),    ntimes, t%timevals(t%eindex), &
              'Last  Spread Value, Time Index/Value: ', spread_array(ntimes), ntimes, t%timevals(t%eindex))

! if the namelist 'skip_first_ntimes' is greater than 0, compute error again
! without those values.  intended use is to skip spinup spikes in error values.

if (skip_first_ntimes > 0 .and. ntimes > skip_first_ntimes+1) then
   skip_start = skip_first_ntimes+1
   call find_error(truth(:, skip_start:ntimes), &
                   ens_mean(:, skip_start:ntimes), err_array(1:ntimes-skip_start+1), total_error)
   call find_error(zeros(:, skip_start:ntimes), &
                   ens_spread(:, skip_start:ntimes), spread_array(1:ntimes-skip_start+1), total_spread)


   message1 = construct_msg_string("Skip First", skip_first_ntimes, "Times,")
   call error_handler(E_MSG,'', '')
   call printout(trim(message1)//' Mean  Error: ', total_error,  -1, -1.0_r8, &
                 trim(message1)//' Mean Spread: ', total_spread, -1, -1.0_r8)
endif

! FIXME: make an error handler call that just prints a blank line?  or error_handler(E_BLANK)
call error_handler(E_MSG,'', '')

! clean up

deallocate(truth)
deallocate(ens_mean, ens_spread, zeros, err_array, spread_array)

ierr = nf90_close(t%ncid)
call nc_check(ierr, 'compute_error', 'nc_close truth')
ierr = nf90_close(d%ncid) 
call nc_check(ierr, 'compute_error', 'nc_close diag')


call finalize_utilities('compute_error')

contains

!---------------------------------------------------

! the routine that computes the time-series error array,
! plus returning the mean of that array.

subroutine find_error(truth, diag, err_array, total_err)

real(r8), intent(in)  :: truth(:,:)
real(r8), intent(in)  :: diag(:,:)
real(r8), intent(out) :: err_array(:)
real(r8), intent(out) :: total_err

integer :: model_size, num_times, i, j
real(r8), allocatable :: sumsq(:)
 
model_size = size(truth, 1)
num_times = size(truth, 2)

!FIXME: are missing values in the truth or state possible?
! check pop or clm to test this.

allocate(sumsq(num_times))

! in matlab:
!   err = sqrt( sum( (pred - verif).^2, 2) );

sumsq(:) = 0.0_r8
do j=1, num_times
   do i=1, model_size
      sumsq(j) = sumsq(j) + (truth(i,j)-diag(i,j))**2
   enddo
enddo

err_array(:) = sqrt(sumsq)

total_err = sum(err_array) / num_times

deallocate(sumsq)

end subroutine find_error

!-----------------------------------------------------------------------

! given an array of string names and a string to find, return which index
! in that array it matches.  it is a fatal error to not find the given
! string in the array.

subroutine findcopyindex(copyarrayname, copyarray, copylen, tofind, infile, foundindex)

character(len=*), intent(in)  :: copyarrayname
character(len=*), intent(in)  :: copyarray(:)
integer,          intent(in)  :: copylen
character(len=*), intent(in)  :: tofind
character(len=*), intent(in)  :: infile
integer,          intent(out) :: foundindex

integer :: i

foundindex = -1
do i=1, copylen
   if (index(copyarray(i), tofind) /= 0) then
      foundindex = i
      exit
   endif
enddo

if (foundindex < 0) then
   write(message1, '(A)') 'does not contain a copy labeled "'//trim(tofind)//'"'
   write(message2, '(A)') 'Use "ncdump -v '//trim(copyarrayname)//' '//trim(infile)//'" to see existing copy names'
   call error_handler(E_ERR, 'findcopyindex', 'error in file "'//trim(infile)//'"', &
                      source, text2=message1, text3=message2)
endif

end subroutine findcopyindex

!-----------------------------------------------------------------------

! find the dimension lengths for a set of hardcoded dimension names 
! that i know should be in a dart diagnostic file.  it's a fatal error
! if any of these are missing.

subroutine findlengths(this)
 type(fileinfo), intent(inout) :: this

integer :: ierr, dimid

ierr = nf90_inq_dimid(this%ncid, "copy", dimid=dimid)
call nc_check(ierr, 'findlengths', 'inq_dimid copy')
ierr = nf90_inquire_dimension(this%ncid, dimid, len=this%numcopies)
call nc_check(ierr, 'findlengths', 'inquire_dimension copy')

ierr = nf90_inq_dimid(this%ncid, "metadatalength", dimid=dimid)
call nc_check(ierr, 'findlengths', 'inq_dimid metadatalength')
ierr = nf90_inquire_dimension(this%ncid, dimid, len=this%metadatacount)
call nc_check(ierr, 'findlengths', 'inquire_dimension metadatalength')

ierr = nf90_inq_dimid(this%ncid, "time", dimid=dimid)
call nc_check(ierr, 'findlengths', 'inq_dimid time')
ierr = nf90_inquire_dimension(this%ncid, dimid, len=this%ntimes)
call nc_check(ierr, 'findlengths', 'inquire_dimension time')

ierr = nf90_inq_dimid(this%ncid, "StateVariable", dimid=dimid)
call nc_check(ierr, 'findlengths', 'inq_dimid StateVariable')
ierr = nf90_inquire_dimension(this%ncid, dimid, len=this%model_size)
call nc_check(ierr, 'findlengths', 'inquire_dimension statelen')

end subroutine findlengths

!-----------------------------------------------------------------------

! make sure model sizes are the same between the two files.  we can deal
! with different length time series (assuming the times overlap at least
! some) but the model sizes have to match identically before we proceeed.

subroutine confirm_modelsizes_match(this, that)

type(fileinfo), intent(in) :: this
type(fileinfo), intent(in) :: that

if (this%model_size /= that%model_size) then
   write(message1, '(A)') 'model sizes must be identical in input files'
   write(message2, '(A,I10)') 'sizes are ', this%model_size, that%model_size
   call error_handler(E_ERR, 'confirm_modelsizes_match', message1, &
                      source, text2=message2)
endif

end subroutine confirm_modelsizes_match

!-----------------------------------------------------------------------

! small cover routine to inquire and then return character data.

subroutine getchardata_nc(ncid, varname, chardata)

integer,           intent(in)  :: ncid
character(len=*),  intent(in)  :: varname
character(len=*),  intent(out) :: chardata(:)

integer :: ierr, varid

ierr = nf90_inq_varid(ncid, trim(varname), varid=varid)
call nc_check(ierr, 'getchardata_nc', 'inq_varid chardata')
ierr = nf90_get_var(ncid, varid, chardata)
call nc_check(ierr, 'getchardata_nc', 'get_var chardata')

end subroutine getchardata_nc

!-----------------------------------------------------------------------

! small cover routine to inquire and then return a 1d integer array.

subroutine getintdata1d_nc(ncid, varname, intdata)

integer,           intent(in)  :: ncid
character(len=*),  intent(in)  :: varname
integer,           intent(out) :: intdata(:)

integer :: ierr, varid

ierr = nf90_inq_varid(ncid, varname, varid=varid)
call nc_check(ierr, 'getintdata_nc', 'inq_varid intdata')
ierr = nf90_get_var(ncid, varid, intdata)
call nc_check(ierr, 'getintdata_nc', 'get_var intdata')

end subroutine getintdata1d_nc

!-----------------------------------------------------------------------

! small cover routine to inquire and then return a 1d real array.

subroutine getrealdata1d_nc(ncid, varname, realdata)

integer,           intent(in)  :: ncid
character(len=*),  intent(in)  :: varname
real(r8),          intent(out) :: realdata(:)

integer :: ierr, varid

ierr = nf90_inq_varid(ncid, varname, varid=varid)
call nc_check(ierr, 'getrealdata1d_nc', 'inq_varid realdata')
ierr = nf90_get_var(ncid, varid, realdata)
call nc_check(ierr, 'getrealdata1d_nc', 'get_var realdata')

end subroutine getrealdata1d_nc

!-----------------------------------------------------------------------

! small cover routine to inquire and then return a 2d real array.

subroutine getrealdata2d_nc(ncid, varname, realdata)

integer,           intent(in)  :: ncid
character(len=*),  intent(in)  :: varname
real(r8),          intent(out) :: realdata(:,:)

integer :: ierr, varid

ierr = nf90_inq_varid(ncid, varname, varid=varid)
call nc_check(ierr, 'getrealdata2d_nc', 'inq_varid realdata')
ierr = nf90_get_var(ncid, varid, realdata)
call nc_check(ierr, 'getrealdata2d_nc', 'get_var realdata')

end subroutine getrealdata2d_nc

!-----------------------------------------------------------------------

! sets ntimes to the number of time vals in the common time array,
! and also sets sindex and eindex in both this and that to indicate
! the exact times which are overlapped.

subroutine findtimesubset(this, that, ntimes)

type(fileinfo), intent(inout) :: this
type(fileinfo), intent(inout) :: that
integer,        intent(out)   :: ntimes

integer :: i, j
real(r8) :: thistime
logical :: timematch

! set all the outputs here

this%sindex = -1
this%eindex = -1
that%sindex = -1
that%eindex = -1
ntimes = 0

!DEBUG print *, 'truth: ', this%ntimes, this%timevals(1), this%timevals(this%ntimes)
!DEBUG print *, ' diag: ', that%ntimes, that%timevals(1), that%timevals(that%ntimes)

! see if the common case is true: that both time arrays are the same
! length and equal. if so, return happily.

if (this%ntimes == that%ntimes) then
   timematch = .true.
   LOOP1: do i=1, this%ntimes
      if (abs(this%timevals(i) - that%timevals(i)) > epsilon(thistime)) then
         timematch = .false.
         exit LOOP1
      endif
   enddo LOOP1
   if (timematch) then
      ntimes = this%ntimes
      this%sindex = 1
      that%sindex = 1
      this%eindex = ntimes
      that%eindex = ntimes
      return
   endif 
endif

! see if another case is true - that these time lists are disjoint sets.

if (this%timevals(1) > that%timevals(that%ntimes)) then
   !DEBUG print *, this%timevals(1), that%timevals(that%ntimes), ' a > b '
   write(message1, '(A)') 'all times in the diag file are before times in the truth file'
   write(message2, '(A)') 'no common times in the two input files'
   call error_handler(E_ERR, 'compare_lengths', message1, &
                      source, text2=message2)
endif

if (this%timevals(this%ntimes) < that%timevals(1)) then
   !DEBUG print *, this%timevals(this%ntimes), that%timevals(1), ' a < b '
   write(message1, '(A)') 'all times in the diag file are after times in the truth file'
   write(message2, '(A)') 'no common times in the two input files'
   call error_handler(E_ERR, 'compare_lengths', message1, &
                      source, text2=message2)
endif

! find which one has the earliest time and loop over the other looking for
! a match of the start time.

if (this%timevals(1) < that%timevals(1)) then

   LOOP2: do i=1, this%ntimes
      if (abs(this%timevals(i) - that%timevals(1)) <= epsilon(thistime)) then
         this%sindex = i
         that%sindex = 1
         exit LOOP2 
      endif
   enddo LOOP2
   if (this%sindex < 0) goto 10  ! error

else

   LOOP3: do i=1, this%ntimes
      if (abs(this%timevals(1) - that%timevals(i)) <= epsilon(thistime)) then
         this%sindex = 1
         that%sindex = i
         exit LOOP3 
      endif
   enddo LOOP3
   if (this%sindex < 0) goto 10  ! error

endif

!DEBUG print *, 'set sindex ok, to vals', this%sindex, that%sindex

! sindex is now set.  start from sindex and stop when
! you run out of times or the values differ.

i = this%sindex
j = that%sindex
LOOP4: do 
   !DEBUG print *, i,j, this%timevals(i), that%timevals(j)

   ! check to make sure the time intervals are the same in both files
   if (abs(this%timevals(i) - that%timevals(j)) > epsilon(thistime)) goto 10 ! error

   if (i+1 > this%ntimes) exit LOOP4
   if (j+1 > that%ntimes) exit LOOP4

   i = i + 1
   j = j + 1

enddo LOOP4

this%eindex = i
that%eindex = j

!DEBUG print *, 'set eindex ok, to vals', this%eindex, that%eindex

this%ntimes = this%eindex - this%sindex + 1
that%ntimes = that%eindex - that%sindex + 1

ntimes = this%ntimes
!DEBUG print *, 'set ntimes to ', ntimes
return

! you only get here via an error goto

10 continue

write(message1, '(A)') 'there is an overlapping time period, but the exact times do not match'
write(message2, '(A)') 'mismatched times between the two input files'
call error_handler(E_ERR, 'compare_lengths', message1, &
                   source, text2=message2)

end subroutine findtimesubset

!-----------------------------------------------------------------------

subroutine printout(str1, rval1, ival1, rval1a, str2, rval2, ival2, rval2a)

character(len=*),  intent(in) :: str1
real(r8),          intent(in) :: rval1
integer,           intent(in) :: ival1
real(r8),          intent(in) :: rval1a
character(len=*),  intent(in) :: str2
real(r8),          intent(in) :: rval2
integer,           intent(in) :: ival2
real(r8),          intent(in) :: rval2a

call error_handler(E_MSG,'', '')
if (ival1 >= 0) then
   write(message1, *) str1, rval1, ival1, rval1a
else
   write(message1, *) str1, rval1
endif
call error_handler(E_MSG, '', message1)
if (ival2 >= 0) then
   write(message1, *) str2, rval2, ival2, rval2a
else
   write(message1, *) str2, rval2
endif
call error_handler(E_MSG, '', message1)

end subroutine printout

!-----------------------------------------------------------------------

! a bunch of messing around to construct a string from 2 strings
! plus an integer, where i want the integer to be left justified
! with no additional spaces.  you'd think there'd be a language
! construct to do this, but no.  and there's additional messing
! around because some compilers apparently have a problem with
! nesting adjustl() and trim() even though nesting these functions
! is perfectly legal in the language.

function construct_msg_string(s1, skipcount, s2)

integer,          intent(in) :: skipcount
character(len=*), intent(in) :: s1
character(len=*), intent(in) :: s2
character(len=80)            :: construct_msg_string

character(len=80) :: string1, string2

! some versions of some compilers (possibly absoft, pgi) have problems
! with some nested combinations of adjustl() and trim(). avoid the
! whole issue by bouncing between 2 strings and using them separately.

! the * format adds an initial blank space for the 'hollerith' character.
! the second write will fail because string1 plus 1 character is longer
! than string2.  it is ok with "(A)" because that keeps the original length.

write(string1, *) skip_first_ntimes
write(string2, '(A)') adjustl(string1)
write(string1, '(3(A,1X))') s1, trim(string2), s2

construct_msg_string = string1

end function construct_msg_string

!-----------------------------------------------------------------------

end program compute_error

