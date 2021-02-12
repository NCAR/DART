! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> program to take 2 netCDF DART diagnostic files and compare
!> the state variable, or whatever variables are listed in the
!> namelist.  prints the min, max values and the min/max difference
!> for each field listed.
!>
!> hopefully useful when comparing the results of two parallel
!> experiments.  simple ncdiff balks when the input namelist variable
!> is a different shape, and it cannot easily tell you when two
!> variables are identical.   there are matlab functions that would
!> do this with a short script, but not all platforms have matlab.

program compare_states

use     types_mod, only : r8
use utilities_mod, only : error_handler, E_ERR, E_MSG,  &
                          open_file, close_file, get_next_filename,      &
                          find_namelist_in_file, check_namelist_read,    &
                          do_nml_file, do_nml_term, nmlfileunit,         &
                          initialize_utilities, finalize_utilities
use  netcdf_utilities_mod, only : nc_check
use parse_args_mod, only : get_args_from_string

use netcdf

implicit none

character(len=*), parameter :: source = 'compare_states.f90'

! variables used to read the netcdf info
integer, parameter :: maxd = 7
integer :: i, j, ndims, odims, ncrc, etype, nitems, nvars, xtype
integer :: ncinid1, ncinid2        ! netcdf id for file
integer :: invarid1, invarid2
integer :: dimid(maxd), dimlen(maxd), odimid(maxd), odimlen(maxd)
character(128) :: dimname(maxd), odimname(maxd)
integer :: nin1Dimensions, nin1Variables, nin1Attributes, in1unlimitedDimID
integer :: nin2Dimensions, nin2Variables, nin2Attributes, in2unlimitedDimID
real(r8) ::  min1,  min2,  max1,  max2,  delmin,  delmax
integer  :: imin1, imin2, imax1, imax2, idelmin, idelmax

! arrays and scalars for real and int
real(r8)              ::  zerod1,    zerod2
real(r8), allocatable ::  oned1(:),  oned2(:)       

integer               ::  izerod1,   izerod2
integer,  allocatable ::  ioned1(:), ioned2(:)       

logical, save :: module_initialized = .false.

! arg parsing code
character(len=256) :: argline
integer :: argcount = 2
character(len=NF90_MAX_NAME) :: argwords(3)

character(len=NF90_MAX_NAME) :: infile1, infile2
character(len=NF90_MAX_NAME) :: nextfield
logical :: from_file

character(len=512) :: msgstring, msgstring2, tmpstring
integer :: iunit, io
logical :: debug = .false.                   ! or .true.
logical :: fail_on_missing_field = .true.    ! or .false.
logical :: do_all_numeric_fields = .true.    ! or .false.
logical :: only_report_differences = .true.  ! or .false.
character(len=128) :: fieldnames(1000) = ''  ! something large
character(len=128) :: fieldlist_file = ''

! fieldnames here?
namelist /compare_states_nml/  debug, fail_on_missing_field, &
                               do_all_numeric_fields,        &
                               fieldnames, fieldlist_file,   &
                               only_report_differences

! main code here
 
! flow:
!   initialization
call initialize_utilities('compare_states')
call initialize_module()

! Read the namelist entry
call find_namelist_in_file("input.nml", "compare_states_nml", iunit)
read(iunit, nml = compare_states_nml, iostat = io)
call check_namelist_read(iunit, io, "compare_states_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=compare_states_nml)
if (do_nml_term()) write(     *     , nml=compare_states_nml)

if (debug) then
   call error_handler(E_MSG, 'compare_states', ' debug on')
endif

! whether to fail or just warn if a field is not found
if (fail_on_missing_field) then
   etype = E_ERR
else
   etype = E_MSG
endif

!   check inputs - get 2 white-space separated strings from stdin
!   eg:  echo infile1.nc infile2.nc | ./compare_states
read(*, '(A)') argline
call get_args_from_string(argline, argcount, argwords)

if (argcount /= 2) then
   msgstring = 'Usage: echo infile1.nc infile2.nc | ./compare_states'
   call error_handler(E_ERR, 'compare_states', msgstring, source) 
endif

infile1   = argwords(1)
infile2   = argwords(2)

call error_handler(E_MSG, 'compare_states', ' reading file: '//trim(infile1))
call error_handler(E_MSG, 'compare_states', '  and    file: '//trim(infile2))
if (do_all_numeric_fields) then
   call error_handler(E_MSG, 'compare_states', ' doing all numeric fields')
else
   ! make sure the namelist specifies one or the other but not both
   ! only if we aren't trying to do all fields in the file.
   if (fieldnames(1) /= '' .and. fieldlist_file /= '') then
      call error_handler(E_ERR,'compare_states', &
          'cannot specify both fieldnames and fieldlist_file', source)
   endif
   
   if (fieldlist_file /= '') then
      call error_handler(E_MSG, 'compare_states', ' list of fields file: '//trim(fieldlist_file))
      from_file = .true.
   else
      call error_handler(E_MSG, 'compare_states', ' field names specified in namelist.')
      from_file = .false.
   endif
endif

!   do they exist?  can they be opened?
!   infile1 & infile2 are netcdf
!   fieldlist is ascii, one varname per line

! open the files
call nc_check(nf90_open(infile1, NF90_NOWRITE,   ncinid1), 'nf90_open', 'infile1')
call nc_check(nf90_open(infile2, NF90_NOWRITE,   ncinid2), 'nf90_open', 'infile2')

call nc_check(nf90_inquire(ncinid1, nin1Dimensions, nin1Variables, &
              nin1Attributes, in1unlimitedDimID), 'nf90_inquire', 'infile1')
call nc_check(nf90_inquire(ncinid2, nin2Dimensions, nin2Variables, &
              nin2Attributes, in2unlimitedDimID), 'nf90_inquire', 'infile2')

! for now, loop over the number of vars in file 1.  at some point we
! should print out vars that are in 1 but not 2, and in 2 but not 1.
! but for that we need more logic to track which vars are processed.
! this code loops over the names in file 1 and finds them (or not)
! in file 2 but doesn't complain about unused vars in file 2.
nvars = nin1Variables

if (debug) then
   write(msgstring, *) 'infile1 ndim, nvar, nattr:',nin1Dimensions, &
                      nin1Variables,nin1Attributes
   call error_handler(E_MSG, 'compare_states', msgstring)
   write(msgstring, *) 'infile2 ndim, nvar, nattr:', nin2Dimensions, &
                       nin2Variables, nin2Attributes
   call error_handler(E_MSG, 'compare_states', msgstring)
endif

!    input files to get data from
!    list of netcdf fields to compare

fieldloop : do i=1, 100000

   if (do_all_numeric_fields) then
      if (i > nvars) exit fieldloop
      call nc_check(nf90_inquire_variable(ncinid1, i, nextfield, xtype), &
                    'nf90_inquire_variable', 'infile1')
      if (xtype /= NF90_INT .and. xtype /= NF90_FLOAT .and. xtype /= NF90_DOUBLE) then
         tmpstring = ' not integer, float, or double'
         msgstring = 'skipping variable '//trim(nextfield)//','//trim(tmpstring)
         call error_handler(E_MSG, 'compare_states', msgstring, source) 
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
      call error_handler(etype, 'compare_states', msgstring, source) 
      cycle fieldloop
   endif
   ncrc = nf90_inq_varid(ncinid2, trim(nextfield), invarid2)
   if (ncrc /= NF90_NOERR) then
      tmpstring = ' exists in input1 file '//trim(infile1)
      if (etype == E_ERR) then
         msgstring = 'variable '//trim(nextfield)//trim(tmpstring)
      else
         msgstring = 'skipping variable '//trim(nextfield)//','//trim(tmpstring)
      endif
      msgstring2 = 'but was not found in second file '//trim(infile2)
      call error_handler(etype, 'compare_states', msgstring, &
                         source, text2=msgstring2) 
      cycle fieldloop
   endif

   if (debug) then
      write(msgstring, *) 'invarid1: ', trim(nextfield)//' ',  invarid1
      call error_handler(E_MSG, 'compare_states', msgstring)
      write(msgstring, *) 'invarid2: ', trim(nextfield)//' ', invarid2
      call error_handler(E_MSG, 'compare_states', msgstring)
   endif

   ! get dimensions and make sure they match

   call nc_check(nf90_inquire_variable(ncinid1,  invarid1, ndims=ndims,  dimids=dimid), &
                 'nf90_inquire_variable', 'infile1/'//trim(nextfield))
   call nc_check(nf90_inquire_variable(ncinid2, invarid2, ndims=odims, dimids=odimid), &
                 'nf90_inquire_variable', 'infile2/'//trim(nextfield))

   if (ndims /= odims) then
      write(msgstring, *) 'variable ', trim(nextfield), &
         ' has different numbers of dimensions in the two files'
      call error_handler(E_MSG, 'compare_states', msgstring)
      write(msgstring, *) 'input dimension size ', ndims, ' does not match second file ', odims
      call error_handler(E_ERR, 'compare_states', msgstring, source)
   endif
   
   dimlen(:) = 1
   do j=1,ndims
      call nc_check( nf90_inquire_dimension(ncinid1,  dimid(j),  dimname(j),  dimlen(j)), &
                   'nf90_inquire_dimension', 'infile1/'//trim( dimname(j)) )
      if (debug) then
         write(msgstring, '(2A,I5,A,I8,2A)') trim(infile1), ' dim: ', j, ' len: ', dimlen(j), ' name: ', trim(dimname(j))
         call error_handler(E_MSG, 'compare_states', msgstring)
      endif
      call nc_check( nf90_inquire_dimension(ncinid2, odimid(j), odimname(j), odimlen(j)), &
                   'nf90_inquire_dimension', 'infile2/'//trim(odimname(j)) )
      if (debug) then
         write(msgstring, '(2A,I5,A,I8,2A)') trim(infile2), ' dim: ', j, ' len: ', odimlen(j), ' name: ', trim(odimname(j))
         call error_handler(E_MSG, 'compare_states', msgstring)
      endif
      
      if (dimlen(j) /= odimlen(j)) then
         write(msgstring, *) 'variable ', trim(nextfield), ' has different dimensions in the two files'
         call error_handler(E_MSG, 'compare_states', msgstring)
         write(msgstring, *) 'input dim length ', dimlen(j), ' does not match second file ', odimlen(j)
         call error_handler(E_ERR, 'compare_states', msgstring, source)
      endif

      ! only possible if the unlimited dim is declared but hasn't been written to
      if (dimlen(j) == 0) then
         write(msgstring, *) trim(nextfield), 'will be skipped because it is empty in input file'
         call error_handler(E_MSG, 'compare_states', msgstring)
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
   ! this section has lot of replicated code, but with the strongly typed
   ! arrays it is hard to do much else.  we could overload a subroutine
   ! with each dimension but we'd end up with the same amount of replicated code

   select case(xtype)
    case(NF90_INT)
     select case(ndims)
      case (0)
         call nc_check(nf90_get_var(ncinid1, invarid1, izerod1), 'nf90_get_var', 'infile1')
         call nc_check(nf90_get_var(ncinid2, invarid2, izerod2), 'nf90_get_var', 'infile2')
         imin1 = izerod1
         imax1 = izerod1
         imin2 = izerod2
         imax2 = izerod2
         idelmin = abs(izerod1-izerod2)
         idelmax = abs(izerod1-izerod2)
         if(izerod1 .ne. izerod2) then
            nitems = 1
         else
            nitems = 0
         endif
      case (1:7)
         allocate(ioned1(product(dimlen)), ioned2(product(dimlen)))
         call nc_check(nf90_get_var(ncinid1, invarid1, ioned1, count=dimlen), 'nf90_get_var', 'infile1')
         call nc_check(nf90_get_var(ncinid2, invarid2, ioned2, count=dimlen), 'nf90_get_var', 'infile2')
         imin1 = minval(ioned1)
         imax1 = maxval(ioned1)
         imin2 = minval(ioned2)
         imax2 = maxval(ioned2)
         idelmin = minval(abs(ioned1-ioned2))
         idelmax = maxval(abs(ioned1-ioned2))
         nitems = count(ioned1 .ne. ioned2)
         deallocate(ioned1, ioned2)
     end select
     ! common reporting code for integers
     if (nitems > 0 .or. .not. only_report_differences) then
        write(msgstring, *) 'checking equality of: ', trim(tmpstring)
        call error_handler(E_MSG, 'compare_states', msgstring)
     endif
     if (nitems > 0) then
        write(msgstring, *) 'arrays differ in ', nitems, ' places'
        call error_handler(E_MSG, 'compare_states', msgstring, source)
        write(msgstring, *) 'min/max file1: ', imin1, imax1
        call error_handler(E_MSG, 'compare_states', msgstring, source)
        write(msgstring, *) 'min/max file2: ', imin2, imax2
        call error_handler(E_MSG, 'compare_states', msgstring, source)
        write(msgstring, *) 'delta min/max: ', idelmin, idelmax
        call error_handler(E_MSG, 'compare_states', msgstring, source)
     else if (.not. only_report_differences) then
        write(msgstring, *) 'arrays same'
        call error_handler(E_MSG, 'compare_states', msgstring, source)
        write(msgstring, *) 'min/max value: ', imin1, imax1
        call error_handler(E_MSG, 'compare_states', msgstring, source)
     endif

    case default
     select case(ndims)
      case (0)
         call nc_check(nf90_get_var(ncinid1, invarid1, zerod1), 'nf90_get_var', 'infile1')
         call nc_check(nf90_get_var(ncinid2, invarid2, zerod2), 'nf90_get_var', 'infile2')
         min1 = zerod1
         max1 = zerod1
         min2 = zerod2
         max2 = zerod2
         delmin = abs(zerod1-zerod2)
         delmax = abs(zerod1-zerod2)
         if(zerod1 .ne. zerod2) then
            nitems = 1
         else
            nitems = 0
         endif
      case (1:7)
         allocate(oned1(product(dimlen)), oned2(product(dimlen)))
         call nc_check(nf90_get_var(ncinid1, invarid1, oned1, count=dimlen), 'nf90_get_var', 'infile1')
         call nc_check(nf90_get_var(ncinid2, invarid2, oned2, count=dimlen), 'nf90_get_var', 'infile2')
         min1 = minval(oned1)
         max1 = maxval(oned1)
         min2 = minval(oned2)
         max2 = maxval(oned2)
         delmin = minval(abs(oned1-oned2))
         delmax = maxval(abs(oned1-oned2))
         nitems = count(oned1 .ne. oned2)
         deallocate(oned1, oned2)
     end select

     ! common reporting code for reals
     if (nitems > 0 .or. .not. only_report_differences) then
        write(msgstring, *) 'checking equality of: ', trim(tmpstring)
        call error_handler(E_MSG, 'compare_states', msgstring)
     endif
     if (nitems > 0) then
        write(msgstring, *) 'arrays differ in ', nitems, ' places'
        call error_handler(E_MSG, 'compare_states', msgstring, source)
        write(msgstring, '(A,2(G25.16,1X))') 'Min/Max file1: ', min1, max1
        call error_handler(E_MSG, 'compare_states', msgstring, source)
        write(msgstring, '(A,2(G25.16,1X))') 'Min/Max file2: ', min2, max2
        call error_handler(E_MSG, 'compare_states', msgstring, source)
        write(msgstring, '(A,2(G25.16,1X))') 'delta Min/Max: ', delmin, delmax
        call error_handler(E_MSG, 'compare_states', msgstring, source)
     else if (.not. only_report_differences) then
        write(msgstring, *) 'arrays same'
        call error_handler(E_MSG, 'compare_states', msgstring, source)
        write(msgstring, *) 'min/max value: ', min1, max1
        call error_handler(E_MSG, 'compare_states', msgstring, source)
     endif
   end select

enddo fieldloop

!  close up
call nc_check(nf90_close(ncinid1), 'nf90_close', 'infile1')
call nc_check(nf90_close(ncinid2), 'nf90_close', 'infile2')

if (debug) then
   write(msgstring, *) 'closing files ',  trim(infile1)
   call error_handler(E_MSG, 'compare_states', msgstring)
   write(msgstring, *) 'and ', trim(infile2)
   call error_handler(E_MSG, 'compare_states', msgstring)
endif

call finalize_utilities('compare_states')

! end of main code


contains

!----------------------------------------------------------------------

subroutine initialize_module

  module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------

end program

