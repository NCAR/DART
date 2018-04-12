! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program ncinfo

! program to take 2 netCDF DART diagnostic files and compare
! the state variable, or whatever variables are listed in the
! namelist.  prints the min, max values and the min/max difference
! for each field listed.
!
! hopefully useful when comparing the results of two parallel
! experiments.  simple ncdiff balks when the input namelist variable
! is a different shape, and it cannot easily tell you when two
! variables are identical.   there are matlab functions that would
! do this with a short script, but not all platforms have matlab.

use     types_mod, only : r8
use utilities_mod, only : register_module, error_handler, E_ERR, E_MSG,       &
                          open_file, close_file, nc_check, get_next_filename, &
                          find_namelist_in_file, check_namelist_read,         &
                          do_nml_file, do_nml_term, nmlfileunit,              &
                          initialize_utilities, finalize_utilities
use parse_args_mod, only : get_args_from_string

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! variables used to read the netcdf info
integer, parameter :: maxd = 7
integer :: i, j, ndims, odims, ncrc, etype, nitems, nvars, xtype
integer :: ncinid1, invarid1
integer :: dimid(maxd), dimlen(maxd)
character(128) :: dimname(maxd)
integer :: nin1Dimensions, nin1Variables, nin1Attributes, in1unlimitedDimID
real(r8) ::  min1,  max1,  del1
integer  :: imin1, imax1, idel1

! data arrays
real(r8)              ::  zerod1
real(r8), allocatable ::   oned1(:)
integer               ::  izerod1
integer,  allocatable ::  ioned1(:)

logical, save :: module_initialized = .false.

! arg parsing code
character(len=256) :: argline
integer :: argcount = 1
character(len=NF90_MAX_NAME) :: argwords(2)

character(len=NF90_MAX_NAME) :: infile1
character(len=NF90_MAX_NAME) :: nextfield
logical :: from_file

character(len=512) :: msgstring, msgstring2, tmpstring
integer :: iunit, io
logical :: debug = .false.                   ! or .true.
logical :: fail_on_missing_field = .true.    ! or .false.
logical :: do_all_numeric_fields = .true.    ! or .false.
character(len=128) :: fieldnames(1000) = ''  ! something large
character(len=128) :: fieldlist_file = ''

! fieldnames here?
namelist /ncinfo_nml/  debug, fail_on_missing_field, &
                               do_all_numeric_fields,        &
                               fieldnames, fieldlist_file

! main code here
 
! flow:
!   initialization
call initialize_utilities('ncinfo')
call initialize_module()

! Read the namelist entry
call find_namelist_in_file("input.nml", "ncinfo_nml", iunit)
read(iunit, nml = ncinfo_nml, iostat = io)
call check_namelist_read(iunit, io, "ncinfo_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=ncinfo_nml)
if (do_nml_term()) write(     *     , nml=ncinfo_nml)

if (debug) then
   call error_handler(E_MSG, 'ncinfo', ' debug on')
endif

! whether to fail or just warn if a field is not found
if (fail_on_missing_field) then
   etype = E_ERR
else
   etype = E_MSG
endif

!   check inputs - get one filename string from stdin
!   eg:  echo infile.nc | ./ncinfo
read(*, '(A)') argline
call get_args_from_string(argline, argcount, argwords)

if (argcount /= 1) then
   msgstring = 'Usage: echo infile.nc | ./ncinfo'
   call error_handler(E_ERR, 'ncinfo', msgstring, &
                      source, revision, revdate) 
endif

infile1   = argwords(1)

call error_handler(E_MSG, 'ncinfo', ' reading file: '//trim(infile1))
if (do_all_numeric_fields) then
   call error_handler(E_MSG, 'ncinfo', ' doing all numeric fields')
else
   ! make sure the namelist specifies one or the other but not both
   ! only if we aren't trying to do all fields in the file.
   if (fieldnames(1) /= '' .and. fieldlist_file /= '') then
      call error_handler(E_ERR,'ncinfo', &
          'cannot specify both fieldnames and fieldlist_file', &
          source,revision,revdate)
   endif
   
   if (fieldlist_file /= '') then
      call error_handler(E_MSG, 'ncinfo', ' list of fields file: '//trim(fieldlist_file))
      from_file = .true.
   else
      call error_handler(E_MSG, 'ncinfo', ' field names specified in namelist.')
      from_file = .false.
   endif
endif

!   does it exist?  can it be opened?
!   fieldlist is ascii, one varname per line

! open the files
call nc_check(nf90_open(infile1, NF90_NOWRITE,   ncinid1), 'nf90_open', 'infile1')

call nc_check(nf90_inquire(ncinid1, nin1Dimensions, nin1Variables, &
              nin1Attributes, in1unlimitedDimID), 'nf90_inquire', 'infile1')

nvars = nin1Variables

if (debug) then
   write(msgstring, *) 'infile1 ndim, nvar, nattr:',nin1Dimensions, &
                      nin1Variables,nin1Attributes
   call error_handler(E_MSG, 'ncinfo', msgstring)
endif

!    input fields to get data from

fieldloop : do i=1, 100000

   if (do_all_numeric_fields) then
      if (i > nvars) exit fieldloop
      call nc_check(nf90_inquire_variable(ncinid1, i, nextfield, xtype), &
                    'nf90_inquire_variable', 'infile1')
      if (xtype /= NF90_INT .and. xtype /= NF90_FLOAT .and. xtype /= NF90_DOUBLE) then
         tmpstring = ' not integer, float, or double'
         msgstring = 'skipping variable '//trim(nextfield)//','//trim(tmpstring)
         call error_handler(E_MSG, 'ncinfo', msgstring, source, revision, revdate) 
         cycle fieldloop
      endif
   else
      if (from_file) then
         nextfield = get_next_filename(fieldlist_file, i)
      else
         nextfield = fieldnames(i)
      endif
      if (nextfield == '') exit fieldloop
   endif

   ! inquire in inputs for fieldname
   ncrc = nf90_inq_varid(ncinid1, trim(nextfield),  invarid1)
   if (ncrc /= NF90_NOERR) then
      tmpstring = ' not found in input file '//trim(infile1)
      if (etype == E_ERR) then
         msgstring = 'variable '//trim(nextfield)//trim(tmpstring)
      else
         msgstring = 'skipping variable '//trim(nextfield)//','//trim(tmpstring)
      endif
      call error_handler(etype, 'ncinfo', msgstring, source, revision, revdate) 
      cycle fieldloop
   endif
   if (debug) then
      write(msgstring, *) 'invarid1: ', trim(nextfield)//' ',  invarid1
      call error_handler(E_MSG, 'ncinfo', msgstring)
   endif

   ! get dimensions

   call nc_check(nf90_inquire_variable(ncinid1,  invarid1, ndims=ndims,  dimids=dimid), &
                 'nf90_inquire_variable', 'infile1/'//trim(nextfield))

   dimlen(:) = 1
   do j=1,ndims
      call nc_check( nf90_inquire_dimension(ncinid1,  dimid(j),  dimname(j),  dimlen(j)), &
                   'nf90_inquire_dimension', 'infile1/'//trim( dimname(j)) )
      if (debug) then
         write(msgstring, '(2A,I5,A,I8,2A)') trim(infile1), ' dim: ', j, ' len: ', dimlen(j), ' name: ', trim(dimname(j))
         call error_handler(E_MSG, 'ncinfo', msgstring)
      endif

      ! only possible if the unlimited dim is declared but hasn't been written to
      if (dimlen(j) == 0) then
         write(msgstring, *) trim(nextfield), 'will be skipped because it is empty in input file'
         call error_handler(E_MSG, 'ncinfo', msgstring)
         cycle fieldloop
      endif

   enddo

   if (ndims == 0) then
      write(tmpstring, '(2A)') trim(nextfield), ' [scalar value]'
   else
      write(tmpstring, '(2A,1I8)') trim(nextfield), '(', dimlen(1)
      do j=2, ndims
         write(tmpstring, '(2A,2I8)') trim(tmpstring), ',', dimlen(j)
      enddo
      write(tmpstring, '(2A)') trim(tmpstring), ')'
   endif

   ! allocate right dim array
   ! read/write and then deallocate

   select case(xtype)
     case(NF90_INT)
       select case(ndims)
         case (0)
           call nc_check(nf90_get_var(ncinid1, invarid1, izerod1), 'nf90_get_var', 'infile1')
           imin1 = izerod1
           imax1 = izerod1
         case (1:7)
           allocate(ioned1(product(dimlen)))
           call nc_check(nf90_get_var(ncinid1, invarid1, ioned1, count=dimlen), 'nf90_get_var', 'infile1')
           imin1 = minval(ioned1)
           imax1 = maxval(ioned1)
           deallocate(ioned1)
       end select

       ! common reporting code for integers
       write(msgstring, *) 'info for: ', trim(tmpstring)
       call error_handler(E_MSG, 'ncinfo', msgstring)
       write(msgstring, *) 'min/max/delta: ', imin1, imax1, imax1-imin1
       call error_handler(E_MSG, 'ncinfo', msgstring, source, revision, revdate)

     case default
       select case(ndims)
         case (0)
           call nc_check(nf90_get_var(ncinid1, invarid1, zerod1), 'nf90_get_var', 'infile1')
           min1 = zerod1
           max1 = zerod1
         case (1:7)
           allocate(oned1(product(dimlen)))
           call nc_check(nf90_get_var(ncinid1, invarid1, oned1, count=dimlen), 'nf90_get_var', 'infile1')
           min1 = minval(oned1)
           max1 = maxval(oned1)
           deallocate(oned1)
       end select

       ! common reporting code for reals
       write(msgstring, *) 'info for: ', trim(tmpstring)
       call error_handler(E_MSG, 'ncinfo', msgstring)
       write(msgstring, *) 'min/max/delta: ', min1, max1, max1-min1
       call error_handler(E_MSG, 'ncinfo', msgstring, source, revision, revdate)

   end select

enddo fieldloop

!  close up
call nc_check(nf90_close(ncinid1), 'nf90_close', 'infile1')

if (debug) then
   write(msgstring, *) 'closing file ',  trim(infile1)
   call error_handler(E_MSG, 'ncinfo', msgstring)
endif

call finalize_utilities('ncinfo')

! end of main code


contains

!----------------------------------------------------------------------

subroutine initialize_module

  call register_module(source, revision, revdate)
  module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
