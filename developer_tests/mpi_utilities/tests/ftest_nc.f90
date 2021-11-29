! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program ftest_nc

! fortran program which uses netcdf (and NOT mpi). netcdf 3.6.x and beyond
! comes with a fortran 90 module to define netcdf interfaces and
! parameters.  typically the module is found by giving the compiler
! a flag pointing to the netcdf include directory, but sometimes 
! it can be the lib dir.

! if you get an error about the netcdf or typesizes module not
! being found, you need to set NETCDF in your mkmf.template file
! to the right directory, *OR* you can comment NETCDF out of
! your template file and set it in your environment:
!  ksh & related: export NETCDF=xxx
!  csh & related: setenv NETCDF xxx
! if your system uses the 'module' command, be sure you have a
! compatible netCDF module loaded; typically the netCDF installation
! needs to be compiled with the same Fortran compiler as you are using.

use netCDF
use typeSizes

implicit none

character(len=32) :: filename = "ftestdata.nc"
integer :: ncfileid, istat, i
integer :: test1_length = 5
integer :: test1dimid, dataid

   print *, "program start"


! Typical sequence:
! NF90_OPEN             ! create netCDF dataset: enter define mode
!    NF90_def_dim       ! define dimensions: from name and length
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

   ! netcdf test to be sure this system can create all necessary kinds
   if(.not. byteSizesOK()) then
       print *, 'Compiler does not support required kinds of variables.'
       stop
   end if

!-------------------------------------------------------------------------------
! Open/Create file
!-------------------------------------------------------------------------------

   istat = nf90_create(path = trim(filename), cmode = nf90_share, &
                       ncid = ncfileid)
   if (istat /= nf90_noerr) call netcdf_error_exit(istat)

   print *, 'successfully opened ' // trim(filename)

!-------------------------------------------------------------------------------
! Define dimension(s)
!-------------------------------------------------------------------------------

   istat = nf90_def_dim(ncid=ncfileid, name="test1", &
                        len = test1_length, dimid = test1dimid)
   if (istat /= nf90_noerr) call netcdf_error_exit(istat)

!-------------------------------------------------------------------------------
! Write global attributes 
!-------------------------------------------------------------------------------

   istat = nf90_put_att(ncfileid, NF90_GLOBAL, "title", "netcdf test File")
   if (istat /= nf90_noerr) call netcdf_error_exit(istat)

!-------------------------------------------------------------------------------
! Create variables and attributes.
!-------------------------------------------------------------------------------

   istat = nf90_def_var(ncid=ncfileid, name="data", xtype=nf90_int, &
                 dimids=test1dimid, varid=dataid)
   if (istat /= nf90_noerr) call netcdf_error_exit(istat)

   istat = nf90_put_att(ncfileid, dataid, "long_name", "test data array")
   if (istat /= nf90_noerr) call netcdf_error_exit(istat)

!-------------------------------------------------------------------------------
! Leave define mode so we can fill
!-------------------------------------------------------------------------------
   istat = nf90_enddef(ncfileid)
   if (istat /= nf90_noerr) call netcdf_error_exit(istat)

!-------------------------------------------------------------------------------
! Fill the coordinate variables.
! The time variable is filled as time progresses.
!-------------------------------------------------------------------------------

   istat = nf90_put_var(ncfileid, dataid, (/ (i,i=1,5) /) )
   if (istat /= nf90_noerr) call netcdf_error_exit(istat)

!-------------------------------------------------------------------------------
! Close
!-------------------------------------------------------------------------------

   istat = nf90_close(ncfileid)
   if (istat /= nf90_noerr) call netcdf_error_exit(istat)

   print *, "netcdf file successfully closed"

   print *, "program end"

contains

subroutine netcdf_error_exit(istat)
   integer, intent(in) :: istat

   character(len=129) :: error_msg

   error_msg = 'netcdf error, string is: ' // nf90_strerror(istat)

   print *, error_msg
   stop

end subroutine netcdf_error_exit

end program ftest_nc

