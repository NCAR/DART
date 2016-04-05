! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program netcdf_to_openggcm

!-----------------------------------------------------------------------
! purpose: interface between DART and the openggcm model
!
! method: Read DART state vector and overwrite values in a openggcm restart file.
!         If the DART state vector has an 'advance_to_time' present, a
!         file called openggcm_in.DART is created with a time_manager_nml namelist
!         appropriate to advance openggcm to the requested time.
!
!         The netcdf_to_openggcm_nml namelist setting for advance_time_present
!         determines whether or not the input file has an 'advance_to_time'.
!         Typically, only temporary files like 'assim_model_state_ic' have
!         an 'advance_to_time'.
!
! author: Tim Hoar 23 Feb 16
!-----------------------------------------------------------------------

use        types_mod, only : r4, r8, digits12, i8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, error_handler, E_ERR, E_MSG, nc_check, &
                             file_exist, open_file, close_file
use time_manager_mod, only : time_type, print_time, print_date, set_calendar_type, &
                             set_date, get_date, set_time, get_time, &
                             operator(+)

use netcdf_utilities, only : rd_netcdf

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character(len=256) :: netcdf_to_openggcm_input_file  = 'openggcm.nc'
character(len=256) :: netcdf_to_openggcm_output_file = 'openggcm.pot.dat'
logical            :: verbose  = .false.

namelist /netcdf_to_openggcm_nml/ netcdf_to_openggcm_input_file, &
                                  netcdf_to_openggcm_output_file, &
                                  verbose

!-----------------------------------------------------------------------
! global storage
!-----------------------------------------------------------------------

integer               :: iunit, io, dimi
type(time_type)       :: model_time
real(r4), allocatable :: statevector(:,:)

integer  :: nthe, nphi
integer  :: iyear, imonth, iday, ihour, iminute, isecond
integer  :: ncid, VarID
integer  :: ncNdims
integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimlens
character(len=NF90_MAX_NAME), dimension(NF90_MAX_VAR_DIMS) :: dimnames

character(len=512) :: string1, string2, string3

!=======================================================================

call initialize_utilities(progname='netcdf_to_openggcm')

!-----------------------------------------------------------------------
! Read the namelist to get the input and output filenames.

call find_namelist_in_file("input.nml", "netcdf_to_openggcm_nml", iunit)
read(iunit, nml = netcdf_to_openggcm_nml, iostat = io)
call check_namelist_read(iunit, io, "netcdf_to_openggcm_nml")

if (verbose) then
   write(string1,*)"..  converting  netcdf file '"//trim(netcdf_to_openggcm_input_file)//"'"
   write(string2,*)"to openggcm binary file '"//trim(netcdf_to_openggcm_output_file)//"'"
   call error_handler(E_MSG,'netcdf_to_openggcm',string1, text2=string2)
endif

if (file_exist(netcdf_to_openggcm_output_file)) then
   write(string1,*)"..  WARNING: openggcm binary file '"//trim(netcdf_to_openggcm_output_file)//"' EXISTS."
   write(string2,*)"WARNING:  '"//trim(netcdf_to_openggcm_output_file)//"' is being OVERWRITTEN."
   call error_handler(E_MSG,'netcdf_to_openggcm',string1, text2=string2)
endif

!-----------------------------------------------------------------------
! Reads the state, should check the shape if the output file already exists

call nc_check(nf90_open(netcdf_to_openggcm_input_file,NF90_NOWRITE, ncid), &
         'netcdf_to_openggcm','open '//trim(netcdf_to_openggcm_input_file))

call nc_check(nf90_inq_varid(ncid, 'pot', VarID), 'netcdf_to_openggcm', 'inq_varid "pot"')

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=ncNdims), &
            'netcdf_to_openggcm', 'inquire_variable "pot"')

if (ncNdims /= 2) then
   write(string1,*)"'pot' variable has unexpected rank."
   write(string2,*)'Expected a rank 2 variable, is rank ',ncNdims
   call error_handler(E_ERR,'netcdf_to_openggcm',string1, source, revision, revdate, text2=string2)
endif

nthe = 0
nphi = 0
do dimi = 1,ncNdims
   write(string1,'(''inquire "pot" dimension'',i2,A)') dimi
   call nc_check(nf90_inquire_dimension(ncid,dimIDs(dimi),name=dimnames(dimi),len=dimlens(dimi)), &
                    'netcdf_to_openggcm', string1)
   if (trim(dimnames(dimi)) == 'ig_lat') nthe = dimlens(dimi)
   if (trim(dimnames(dimi)) == 'ig_lon') nphi = dimlens(dimi)
enddo

if (nthe < 1) then
   write(string1,*)"dimension 'ig_lat' not found"
   call error_handler(E_ERR,'netcdf_to_openggcm',string1, source, revision, revdate)
endif
if (nphi < 1) then
   write(string1,*)"dimension 'ig_lon' not found"
   call error_handler(E_ERR,'netcdf_to_openggcm',string1, source, revision, revdate)
endif

allocate(statevector(nphi,nthe))

call rd_netcdf(ncid,'pot',nphi,nthe,statevector)

call set_calendar_type('Gregorian')
call get_model_time(ncid, model_time)
call get_date(model_time, iyear, imonth, iday, ihour, iminute, isecond)

call nc_check(nf90_close(ncid),'netcdf_to_openggcm', 'close')

! This is the inverse of the
iyear = iyear - 1900

iunit = open_file(netcdf_to_openggcm_output_file, form='unformatted',action='write')
write(iunit)nphi,nthe
write(iunit)iyear, imonth, iday, ihour, iminute, isecond
write(iunit)statevector
call close_file(iunit)

!-----------------------------------------------------------------------
! Log what we think we've done, and exit.
!-----------------------------------------------------------------------

if ( verbose ) then
   call print_time(model_time,'netcdf_to_openggcm:model time')
   call print_time(model_time,'netcdf_to_openggcm:model time',logfileunit)
   call print_date(model_time,'netcdf_to_openggcm:model date')
   call print_date(model_time,'netcdf_to_openggcm:model date',logfileunit)
endif

call finalize_utilities('netcdf_to_openggcm')

!=======================================================================
contains
!=======================================================================

subroutine get_model_time(ncid, model_time)
integer, intent(in)  :: ncid
type(time_type), intent(out) :: model_time

integer :: varid
real(digits12) :: length_of_run
character(len=NF90_MAX_NAME) :: unitstring
type(time_type) :: model_time_base, model_time_offset
integer :: iyear, imonth, iday, isecond
integer(i8) :: ibigday

call nc_check(nf90_inq_varid(ncid, 'time', varid),            'get_model_time','inq_varid time')
call nc_check(nf90_get_var(ncid, varid, values=length_of_run),'get_model_time','get_var time')
call nc_check(nf90_get_att(ncid, varid, 'units', unitstring), 'get_model_time','get_att time:units')

if (unitstring(1:14) == 'seconds since ') then
   read(unitstring,'(14x,i4,2(1x,i2))',iostat=io)iyear,imonth,iday
   if (io /= 0) then
      write(string1,*)'Unable to parse year/month/day from units. Error status was ',io
      write(string2,*)'expected something like "seconds since YYYY-MM-DD"'
      write(string3,*)'was                     "'//trim(unitstring)//'"'
      call error_handler(E_ERR, 'get_model_time', string1, &
          source, revision, revdate, text2=string2, text3=string3)
   endif
else
   write(string1,*)'Unable to read time units.'
   write(string2,*)'expected "seconds since YYYY-MM-DD"'
   write(string3,*)'was      "'//trim(unitstring)//'"'
   call error_handler(E_ERR, 'get_model_time', string1, &
          source, revision, revdate, text2=string2, text3=string3)
endif

model_time_base = set_date(iyear, imonth, iday)
ibigday         = int(length_of_run,i8)/86400_i8
iday            = int(ibigday)
isecond         = length_of_run - real(iday,digits12)*86400.0_digits12
model_time_offset = set_time(isecond, iday)
model_time = model_time_base + model_time_offset

call print_date(model_time,'get_model_time:netcdf model date')
call print_time(model_time,'get_model_time:DART   model time')

end subroutine get_model_time


end program netcdf_to_openggcm

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
