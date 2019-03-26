! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program replace_wrf_fields

use     types_mod, only : r8
use utilities_mod, only : register_module, error_handler, E_ERR, E_MSG,       &
                          open_file, close_file, get_next_filename,           &
                          find_namelist_in_file, check_namelist_read,         &
                          do_nml_file, do_nml_term, nmlfileunit,              &
                          initialize_utilities, finalize_utilities
use  netcdf_utilities_mod, only : nc_check
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
integer :: i, j, ndims, odims, ncrc, etype
integer :: ncinid, ncoutid        ! netcdf id for file
integer :: invarid, outvarid
integer :: dimid(maxd), dimlen(maxd), odimid(maxd), odimlen(maxd)
character(128) :: dimname(maxd), odimname(maxd)
integer ::  ninDimensions,  ninVariables,  ninAttributes,  inunlimitedDimID
integer :: noutDimensions, noutVariables, noutAttributes, outunlimitedDimID

! arrays for all possible dimensions
real(r8), pointer ::   oned(:)       
real(r8), pointer ::   twod(:,:)
real(r8), pointer :: threed(:,:,:)
real(r8), pointer ::  fourd(:,:,:,:)
real(r8), pointer ::  fived(:,:,:,:,:)
real(r8), pointer ::   sixd(:,:,:,:,:,:)
real(r8), pointer :: sevend(:,:,:,:,:,:,:)

logical, save :: module_initialized = .false.

! arg parsing code
character(len=256) :: argline
integer :: argcount = 2
character(len=NF90_MAX_NAME) :: argwords(3)

character(len=NF90_MAX_NAME) :: infile, outfile
character(len=NF90_MAX_NAME) :: nextfield
logical :: from_file

character(len=128) :: msgstring, msgstring2, tmpstring
integer :: iunit, io
logical :: debug = .false.                   ! or .true.
logical :: fail_on_missing_field = .true.    ! or .false.
character(len=128) :: fieldnames(1000) = ''  ! something large
character(len=128) :: fieldlist_file = ''

! fieldnames here?
namelist /replace_wrf_fields_nml/  debug, fail_on_missing_field, &
                                   fieldnames, fieldlist_file

! main code here
 
! flow:
!   initialization
call initialize_utilities('replace_wrf_fields')
call initialize_module()

! Read the namelist entry
call find_namelist_in_file("input.nml", "replace_wrf_fields_nml", iunit)
read(iunit, nml = replace_wrf_fields_nml, iostat = io)
call check_namelist_read(iunit, io, "replace_wrf_fields_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=replace_wrf_fields_nml)
if (do_nml_term()) write(     *     , nml=replace_wrf_fields_nml)

if (debug) then
   call error_handler(E_MSG, 'replace_wrf_fields', ' debug on')
endif

! whether to fail or just warn if a field is not found
if (fail_on_missing_field) then
   etype = E_ERR
else
   etype = E_MSG
endif

!   check inputs - get 2 white-space separated strings from stdin
!   eg:  echo infile.nc outfile.nc | ./replace_wrf_fields
read(*, '(A)') argline
call get_args_from_string(argline, argcount, argwords)

if (argcount /= 2) then
   msgstring = 'Usage: echo infile.nc outfile.nc | ./replace_wrf_fields'
   call error_handler(E_ERR, 'replace_wrf_fields', msgstring, &
                      source, revision, revdate) 
endif

infile    = argwords(1)
outfile   = argwords(2)

! make sure the namelist specifies one or the other but not both
if (fieldnames(1) /= '' .and. fieldlist_file /= '') then
   call error_handler(E_ERR,'replace_wrf_fields', &
       'cannot specify both fieldnames and fieldlist_file', &
       source,revision,revdate)
endif

call error_handler(E_MSG, 'replace_wrf_fields', ' reading file: '//trim(infile))
call error_handler(E_MSG, 'replace_wrf_fields', ' overwriting file: '//trim(outfile))
if (fieldlist_file /= '') then
   call error_handler(E_MSG, 'replace_wrf_fields', ' list of fields file: '//trim(fieldlist_file))
   from_file = .true.
else
   call error_handler(E_MSG, 'replace_wrf_fields', ' field names specified in namelist.')
   from_file = .false.
endif

!   do they exist?  can they be opened?
!   infile & outfile are netcdf
!   fieldlist is ascii, one wrf fieldname per line

! open the files
call nc_check(nf90_open( infile, NF90_NOWRITE,    ncinid), 'nf90_open',  'infile')
call nc_check(nf90_open(outfile, NF90_WRITE,     ncoutid), 'nf90_open', 'outfile')

if (debug) then
   call nc_check(nf90_inquire( ncinid,  ninDimensions,  ninVariables, &
                  ninAttributes,  inunlimitedDimID), 'nf90_inquire',  'infile')
   call nc_check(nf90_inquire(ncoutid, noutDimensions, noutVariables, &
                 noutAttributes, outunlimitedDimID), 'nf90_inquire', 'outfile')

   write(msgstring, *) ' infile ndim, nvar, nattr:', ninDimensions, &
                       ninVariables, ninAttributes
   call error_handler(E_MSG, 'replace_wrf_fields', msgstring)
   write(msgstring, *) 'outfile ndim, nvar, nattr:', noutDimensions, &
                       noutVariables, noutAttributes
   call error_handler(E_MSG, 'replace_wrf_fields', msgstring)
endif

!    input file to get data from
!    list of netcdf fields to copy over
!    output file to be updated in place

fieldloop : do i=1, 10000

   if (from_file) then
      nextfield = get_next_filename(fieldlist_file, i)
   else
      nextfield = fieldnames(i)
   endif
   if (nextfield == '') exit fieldloop

   ! inquire in input for fieldname
   ! inquire in output
   ncrc = nf90_inq_varid( ncinid, trim(nextfield),  invarid)
   if (ncrc /= NF90_NOERR) then
      tmpstring = ' not found in input file '//trim(infile)
      if (etype == E_ERR) then
         msgstring = 'variable '//trim(nextfield)//trim(tmpstring)
      else
         msgstring = 'skipping variable '//trim(nextfield)//','//trim(tmpstring)
      endif
      call error_handler(etype, 'replace_wrf_fields', msgstring, source, revision, revdate) 
      cycle fieldloop
   endif
   ncrc = nf90_inq_varid(ncoutid, trim(nextfield), outvarid)
   if (ncrc /= NF90_NOERR) then
      tmpstring = ' exists in input file '//trim(infile)
      if (etype == E_ERR) then
         msgstring = 'variable '//trim(nextfield)//trim(tmpstring)
      else
         msgstring = 'skipping variable '//trim(nextfield)//','//trim(tmpstring)
      endif
      msgstring2 = 'but was not found in output file '//trim(outfile)
      call error_handler(etype, 'replace_wrf_fields', msgstring, &
                         source, revision, revdate, text2=msgstring2) 
      cycle fieldloop
   endif

   if (debug) then
      write(msgstring, *) ' invarid: ', trim(nextfield)//' ',  invarid
      call error_handler(E_MSG, 'replace_wrf_fields', msgstring)
      write(msgstring, *) 'outvarid: ', trim(nextfield)//' ', outvarid
      call error_handler(E_MSG, 'replace_wrf_fields', msgstring)
   endif

   ! get dimensions and make sure they match

   call nc_check(nf90_inquire_variable( ncinid,  invarid, ndims=ndims,  dimids=dimid), &
                 'nf90_inquire_variable',  'infile/'//trim(nextfield))
   call nc_check(nf90_inquire_variable(ncoutid, outvarid, ndims=odims, dimids=odimid), &
                 'nf90_inquire_variable', 'outfile/'//trim(nextfield))

   if (ndims /= odims) then
      write(msgstring, *) 'variable ', trim(nextfield), &
         ' has different numbers of dimensions in the two files'
      call error_handler(E_MSG, 'replace_wrf_fields', msgstring)
      write(msgstring, *) 'input dimension size ', ndims, ' does not match output ', odims
      call error_handler(E_ERR, 'replace_wrf_fields', msgstring, source, revision, revdate)
   endif
   
   do j=1,ndims
      call nc_check( nf90_inquire_dimension( ncinid,  dimid(j),  dimname(j),  dimlen(j)), &
                   'nf90_inquire_dimension',  'infile/'//trim( dimname(j)) )
      if (debug) then
         write(msgstring, '(2A,I5,A,I8,2A)') trim(infile), ' dim: ', j, ' len: ', dimlen(j), ' name: ', trim(dimname(j))
         call error_handler(E_MSG, 'replace_wrf_fields', msgstring)
      endif
      call nc_check( nf90_inquire_dimension(ncoutid, odimid(j), odimname(j), odimlen(j)), &
                   'nf90_inquire_dimension', 'outfile/'//trim(odimname(j)) )
      if (debug) then
         write(msgstring, '(2A,I5,A,I8,2A)') trim(outfile), ' dim: ', j, ' len: ', odimlen(j), ' name: ', trim(odimname(j))
         call error_handler(E_MSG, 'replace_wrf_fields', msgstring)
      endif
      
      if (dimlen(j) /= odimlen(j)) then
         write(msgstring, *) 'variable ', trim(nextfield), ' has different dimensions in the two files'
         call error_handler(E_MSG, 'replace_wrf_fields', msgstring)
         write(msgstring, *) 'input dim length ', dimlen(j), ' does not match output ', odimlen(j)
         call error_handler(E_ERR, 'replace_wrf_fields', msgstring, source, revision, revdate)
      endif

      ! only possible if the unlimited dim is declared but hasn't been written to
      if (dimlen(j) == 0) then
         write(msgstring, *) trim(nextfield), 'will be skipped because it is empty in input file'
         call error_handler(E_MSG, 'replace_wrf_fields', msgstring)
         cycle fieldloop
      endif

   enddo


   select case(ndims)
      case (1)
         write(tmpstring, '(2A,1I5,A)') trim(nextfield), '(', dimlen(1),   ')'
      case (2)
         write(tmpstring, '(2A,2I5,A)') trim(nextfield), '(', dimlen(1:2), ')'
      case (3)
         write(tmpstring, '(2A,3I5,A)') trim(nextfield), '(', dimlen(1:3), ')'
      case (4)
         write(tmpstring, '(2A,4I5,A)') trim(nextfield), '(', dimlen(1:4), ')'
      case (5)
         write(tmpstring, '(2A,5I5,A)') trim(nextfield), '(', dimlen(1:5), ')'
      case (6)
         write(tmpstring, '(2A,6I5,A)') trim(nextfield), '(', dimlen(1:6), ')'
      case (7)
         write(tmpstring, '(2A,7I5,A)') trim(nextfield), '(', dimlen(1:7), ')'
      case default
         ! "can't happen"
         write(msgstring, *) 'array dimension is illegal value: ', ndims
         call error_handler(E_ERR, 'replace_wrf_fields', msgstring, source, revision, revdate)
   end select

   ! announce what we're about to do
   write(msgstring, *) 'copying ', trim(tmpstring)
   call error_handler(E_MSG, 'replace_wrf_fields', msgstring)

   ! allocate right dim array
   ! read/write and then deallocate

   select case(ndims)
      case (1)
         allocate(oned(dimlen(1)))
         call nc_check(nf90_get_var( ncinid,  invarid, oned), 'nf90_get_var',  'infile')
         call nc_check(nf90_put_var(ncoutid, outvarid, oned), 'nf90_put_var', 'outfile')
         deallocate(oned)
      case (2)
         allocate(twod(dimlen(1),dimlen(2)))
         call nc_check(nf90_get_var( ncinid,  invarid, twod), 'nf90_get_var',  'infile')
         call nc_check(nf90_put_var(ncoutid, outvarid, twod), 'nf90_put_var', 'outfile')
         deallocate(twod)
      case (3)
         allocate(threed(dimlen(1),dimlen(2),dimlen(3)))
         call nc_check(nf90_get_var( ncinid,  invarid, threed), 'nf90_get_var',  'infile')
         call nc_check(nf90_put_var(ncoutid, outvarid, threed), 'nf90_put_var', 'outfile')
         deallocate(threed)
      case (4)
         allocate(fourd(dimlen(1),dimlen(2),dimlen(3),dimlen(4)))
         call nc_check(nf90_get_var( ncinid,  invarid, fourd), 'nf90_get_var',  'infile')
         call nc_check(nf90_put_var(ncoutid, outvarid, fourd), 'nf90_put_var', 'outfile')
         deallocate(fourd)
      case (5)
         allocate(fived(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5)))
         call nc_check(nf90_get_var( ncinid,  invarid, fived), 'nf90_get_var',  'infile')
         call nc_check(nf90_put_var(ncoutid, outvarid, fived), 'nf90_put_var', 'outfile')
         deallocate(fived)
      case (6)
         allocate(sixd(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5),dimlen(6)))
         call nc_check(nf90_get_var( ncinid,  invarid, sixd), 'nf90_get_var',  'infile')
         call nc_check(nf90_put_var(ncoutid, outvarid, sixd), 'nf90_put_var', 'outfile')
         deallocate(sixd)
      case (7)
         allocate(sevend(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5),dimlen(6),dimlen(7)))
         call nc_check(nf90_get_var( ncinid,  invarid, sevend), 'nf90_get_var',  'infile')
         call nc_check(nf90_put_var(ncoutid, outvarid, sevend), 'nf90_put_var', 'outfile')
         deallocate(sevend)
      case default
         ! "really can't happen"
         write(msgstring, *) 'array dimension is illegal value: ', ndims
         call error_handler(E_ERR, 'replace_wrf_fields', msgstring, source, revision, revdate)
   end select


enddo fieldloop

!  close up
call nc_check(nf90_close( ncinid), 'nf90_close',  'infile')
call nc_check(nf90_close(ncoutid), 'nf90_close', 'outfile')

if (debug) then
   write(msgstring, *) 'closing files',  trim(infile), ' and ', trim(outfile)
   call error_handler(E_MSG, 'replace_wrf_fields', msgstring)
endif

call finalize_utilities('replace_wrf_fields')

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
