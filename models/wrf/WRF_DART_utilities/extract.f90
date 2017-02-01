! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$
 
PROGRAM extract

use    utilities_mod, only : file_exist, open_file, close_file, &
                             error_handler, E_ERR, initialize_utilities, &
                             register_module, logfileunit, finalize_utilities
use    netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

real, allocatable, dimension(:,:) :: psfc
real :: dt

integer :: iunit, n_file, sn, we

integer           :: ierr, ifile, timeindex
character(len=80) :: file_name, varname, ps_var

character(len=19) :: time_string, time_string_last

!----
!  misc stuff

integer  :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: var_id, outid, TimeDimID, TimeVarID, weDimID, snDimID, out_var_id, ncid, sn_id, we_id, DateDimID

call initialize_utilities('Extract')
call register_module(source, revision, revdate)

ps_var = 'MU'

read(5,*) n_file

iunit = open_file('wrfout.list', action = 'read')

do ifile=1,n_file

   read(iunit, *) file_name

   call check ( nf90_open(file_name, NF90_NOWRITE, ncid) )

   call check ( nf90_inq_dimid(ncid, "south_north", sn_id) )
   call check ( nf90_inquire_dimension(ncid, sn_id, varname, sn) )

   call check ( nf90_inq_dimid(ncid, "west_east", we_id) )
   call check ( nf90_inquire_dimension(ncid, we_id, varname, we) )
   allocate(psfc(we,sn))

   call check( nf90_get_att(ncid, nf90_global, 'DT', dt) )

   call check ( nf90_inq_varid(ncid, ps_var, var_id))
   call check ( nf90_get_var(ncid, var_id, psfc, start = (/ 1, 1, 1/)))
   call check ( nf90_inq_varid(ncid, "Times", var_id))
   call check ( nf90_get_var(ncid, var_id, time_string, start = (/ 1/)))
   ierr = NF90_close(ncid)

   if(file_exist('psfc.nc')) then
      call check( nf90_open('psfc.nc', nf90_write, outid) )
      call check(NF90_Inquire(outid, nDimensions, nVariables, nAttributes, unlimitedDimID))
      call check(NF90_Inq_Varid(outid, "time", TimeVarID))
      call check(NF90_Inquire_Dimension(outid, unlimitedDimID, varname, timeindex))
      call check ( nf90_get_var(outid, TimeVarID, time_string_last, start = (/ 1, timeindex/)))
      timeindex = timeindex + 1
      call check(NF90_Inq_Varid(outid, ps_var, out_var_id))
   else
      call check(nf90_create(path = 'psfc.nc', cmode = nf90_share, ncid = outid))
      call check(nf90_def_dim(ncid=outid, name="time", len = nf90_unlimited, dimid = TimeDimID))
      call check(nf90_def_dim(ncid=outid, name="DateStrLen", len = 19, dimid = DateDimID))
      call check(nf90_def_var(outid, name="time", xtype=nf90_char, &
           dimids = (/ DateDimID, TimeDimID /), varid = TimeVarID) )
      call check(nf90_def_dim(ncid=outid, name="west_east",   len = we, dimid = weDimID))
      call check(nf90_def_dim(ncid=outid, name="south_north", len = sn, dimid = snDimID))
      call check(nf90_def_var(ncid=outid, name=ps_var, xtype=nf90_real, &
           dimids = (/ weDimID, snDimID, TimeDimID /), varid  = out_var_id))
      call check(nf90_put_att(outid, NF90_GLOBAL, "DT", dt))
      call check(nf90_enddef(outid))
      call check(nf90_sync(outid))               ! sync to disk, but leave open
      timeindex = 1

   endif

!!$   if (psfc(1,1) /= 0.0 .and. time_string_last /= time_string) then
      call check(nf90_put_var( outid, TimeVarID, time_string, start=(/ 1, timeindex /) ))
      call check(nf90_put_var( outid, out_var_id, psfc, start=(/ 1, 1, timeindex /) ))
!!$   endif

   ierr = NF90_close(outid)

   deallocate(psfc)

enddo

call close_file(iunit)

write(logfileunit,*)'FINISHED extract.'
write(logfileunit,*)

call finalize_utilities('Extract')   ! closes log file.
 
contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent ( in) :: istatus 
    if(istatus /= nf90_noerr) call error_handler(E_ERR,'extract', &
         trim(nf90_strerror(istatus)), source, revision, revdate)
  end subroutine check

END PROGRAM extract

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
