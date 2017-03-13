! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> common routines used no matter which location module compiled with

!>@todo FIXME  define ALL the VERTISXXX her, add a call to has_vert_choice
!>to the xxx/location_mod, and if no, return without needing routines
!>in each of the specific location routines.  if yes, call query and 
!>match it with the type to return true/false.  remove vert_is_xxx from 
!>all the other location routines!!!  :)
!>will have to replicate the VERTISxxx in the location modules because
!>fortran doesn't allow circular 'use's between modules.  ugh.

!> the key for the netcdf routines is this sits above the location routines.  
!> but that makes a circular dep for VERTIS... if they are here it forces
!> a consistent set of values and removes all the vert_is_xxx routines from
!> modules that don't have any vert choices.  but you have to replicate the
!> parameters since you can't have this code use the publics and have the
!> specific location_mods use VERTISxxx from here.  and depending on which
!> specific module you use, only some or possibly none of the VERTs are defined.
!> best we can do, i guess.  fortran sucks sometimes.

!> can we do something about get_close here? simplify the code for
!> modules which don't have variable vertical options?

module location_utilitiess_mod

use      types_mod, only : r8, MISSING_R8, MISSING_I, PI, RAD2DEG, DEG2RAD, OBSTYPELENGTH
use  utilities_mod, only : error_handler, E_ERR, nc_check, E_MSG

!>@todo FIXME: there should be accessor functions for 5 LocationXXX variables below.
use location_mod, only: location_type, get_location, &
          LocationDims, LocationName, LocationLName, &
          LocationStorageOrder, LocationUnits, has_vert_choice, &
          vert_is_height, vert_is_pressure, vert_is_undef, vert_is_level, &
          vert_is_surface, vert_is_scale_height, query_location


use typeSizes
use netcdf

implicit none
private

public :: nc_write_location_atts, nc_get_location_varids, nc_write_location, &
          nc_write_location_vert, has_vert_choice

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! should import these but they don't exist in low order locations mods
!>@todo define them all here and replicate in the ones which use them.

integer, parameter :: VERTISUNDEF       = -2  ! has no specific vertical location (undefined)
integer, parameter :: VERTISSURFACE     = -1  ! surface value (value is surface elevation in m)
integer, parameter :: VERTISLEVEL       =  1  ! by level
integer, parameter :: VERTISPRESSURE    =  2  ! by pressure (in pascals)
integer, parameter :: VERTISHEIGHT      =  3  ! by height (in meters)
integer, parameter :: VERTISSCALEHEIGHT =  4  ! by scale height (unitless)

interface nc_write_location
   module procedure nc_write_single_location
   module procedure nc_write_multiple_locations
end interface

contains

!----------------------------------------------------------------------------
!> Create and add attributes to a 'location' dimension and variable.

subroutine nc_write_location_atts( ncFileID, fname, unlimDimID ) 
 
integer,           intent(in) :: ncFileID    ! handle to the netcdf file
character(len=*),  intent(in) :: fname       ! file name (for printing purposes)
integer, optional, intent(in) :: unlimDimID  ! handle to the dimension that grows

integer :: LocDimID
integer :: VarID

! define the rank/dimension of the location information
call nc_check(nf90_def_dim(ncid=ncFileID, name='location', len=LocationDims, &
       dimid = LocDimID), 'nc_write_location_atts', 'def_dim:location '//trim(fname))

! Define the location variable and attributes

if (present(unlimDimID)) then
   call nc_check(nf90_def_var(ncFileID, 'location', xtype=nf90_double, &
             dimids=(/ LocDimID, unlimDimID /), varid=VarID), &
            'nc_write_location_atts', 'location:def_var')
else
   call nc_check(nf90_def_var(ncFileID, 'location', xtype=nf90_double, &
             dimids=(/ LocDimID /), varid=VarID), &
            'nc_write_location_atts', 'location:def_var')
endif

call nc_check(nf90_put_att(ncFileID, VarID, 'description', 'location coordinates'), &
              'nc_write_location_basics', 'location:description')
call nc_check(nf90_put_att(ncFileID, VarID, 'location_type', trim(LocationName)), &
              'nc_write_location_atts', 'location:location_type')
call nc_check(nf90_put_att(ncFileID, VarID, 'long_name', trim(LocationLName)), &
              'nc_write_location_atts', 'location:long_name')
call nc_check(nf90_put_att(ncFileID, VarID, 'storage_order', trim(LocationStorageOrder)),  &
              'nc_write_location_atts', 'location:storage_order')
call nc_check(nf90_put_att(ncFileID, VarID, 'units', trim(LocationUnits)),   &
              'nc_write_location_atts', 'location:units')

end subroutine nc_write_location_atts

!----------------------------------------------------------------------------
!> Define the ancillary vertical array and attributes

subroutine nc_write_location_vert( ncFileID, fname )

integer,           intent(in) :: ncFileID    ! handle to the netcdf file
character(len=*),  intent(in) :: fname       ! file name (for printing purposes)

integer :: VarID

call nc_check(nf90_def_var(ncid=ncFileID, name='which_vert', xtype=nf90_int, &
          dimids=(/ nf90_unlimited /), varid=VarID), &
            'nc_write_location_atts', 'which_vert:def_var')

call nc_check(nf90_put_att(ncFileID, VarID, 'long_name', 'vertical coordinate system code'), &
           'nc_write_location_atts', 'which_vert:long_name')
call nc_check(nf90_put_att(ncFileID, VarID, 'VERTISUNDEF', VERTISUNDEF), &
           'nc_write_location_atts', 'which_vert:VERTISUNDEF')
call nc_check(nf90_put_att(ncFileID, VarID, 'VERTISSURFACE', VERTISSURFACE), &
           'nc_write_location_atts', 'which_vert:VERTISSURFACE')
call nc_check(nf90_put_att(ncFileID, VarID, 'VERTISLEVEL', VERTISLEVEL), &
           'nc_write_location_atts', 'which_vert:VERTISLEVEL')
call nc_check(nf90_put_att(ncFileID, VarID, 'VERTISPRESSURE', VERTISPRESSURE), &
           'nc_write_location_atts', 'which_vert:VERTISPRESSURE')
call nc_check(nf90_put_att(ncFileID, VarID, 'VERTISHEIGHT', VERTISHEIGHT), &
           'nc_write_location_atts', 'which_vert:VERTISHEIGHT')
call nc_check(nf90_put_att(ncFileID, VarID, 'VERTISSCALEHEIGHT', VERTISSCALEHEIGHT), &
           'nc_write_location_atts', 'which_vert:VERTISSCALEHEIGHT')

end subroutine nc_write_location_vert

!----------------------------------------------------------------------------
!> Return the LocationVarID and WhichVertVarID variables from a given netCDF file.
!>
!> ncFileId         the netcdf file descriptor
!> fname            the name of the netcdf file (for error messages only)
!> LocationVarID    the integer ID of the 'location' variable in the netCDF file
!> WhichVertVarID   the integer ID of the 'which_vert' variable in the netCDF file

subroutine nc_get_location_varids( ncFileID, fname, LocationVarID, WhichVertVarID )

integer,           intent(in)  :: ncFileID   ! handle to the netcdf file
character(len=*),  intent(in)  :: fname      ! file name (for printing purposes)
integer,           intent(out) :: LocationVarID
integer, optional, intent(out) :: WhichVertVarID

integer :: rc

call nc_check(nf90_inq_varid(ncFileID, 'location', varid=LocationVarID), &
          'nc_get_location_varids', 'inq_varid:location '//trim(fname))

if (present(WhichVertVarID)) then
  rc = nf90_inq_varid(ncFileID, 'which_vert', varid=WhichVertVarID)
  if (rc /= NF90_NOERR) then
      WhichVertVarID = MISSING_I
  endif
endif

end subroutine nc_get_location_varids

!----------------------------------------------------------------------------
!> Writes a SINGLE location to the specified netCDF variable and file.
!> The LocationVarID and WhichVertVarID must be the values returned from
!> the nc_get_location_varids call.

subroutine nc_write_single_location(ncFileID, LocationVarID, loc, locindex, WhichVertVarID)
 
integer,             intent(in) :: ncFileID, LocationVarID
type(location_type), intent(in) :: loc
integer,             intent(in) :: locindex
integer, optional,   intent(in) :: WhichVertVarID

real(r8), dimension(LocationDims) :: locations
integer,  dimension(1) :: intval

locations = get_location( loc ) 

call nc_check(nf90_put_var(ncFileID, LocationVarId, locations, &
              start=(/ 1, locindex /), count=(/ LocationDims, 1 /) ), &
              'nc_write_single_location', 'put_var:location')

if (present(WhichVertVarID)) then
   if (WhichVertVarID /= MISSING_I) then
      intval = query_location(loc, 'WHICH_VERT')
      call nc_check(nf90_put_var(ncFileID, WhichVertVarID, intval, &
                    start=(/ locindex /), count=(/ 1 /) ), &
                    'nc_write_single_location','put_var:vert' )
   endif
endif

end subroutine nc_write_single_location

!----------------------------------------------------------------------------
!> Writes an array of locations to the specified netCDF variable and file.
!> The LocationVarID and WhichVertVarID must be the values returned from
!> the nc_get_location_varids call.

subroutine nc_write_multiple_locations(ncFileID, LocationVarID, loc, loccount, startlocindex, WhichVertVarID)
 
integer,             intent(in) :: ncFileID, LocationVarID
type(location_type), intent(in) :: loc(:)
integer,             intent(in) :: loccount
integer,             intent(in) :: startlocindex
integer, optional,   intent(in) :: WhichVertVarID

real(r8), allocatable :: locations(:,:)
integer,  allocatable :: intvals(:)
logical :: dovert

dovert = .false.
if (present(WhichVertVarID)) then 
   if (WhichVertVarID /= MISSING_I) dovert = .true.
endif

allocate(locations(LocationDims,loccount))
if (dovert) allocate(intvals(loccount))

do i=1, loccount
   locations(:,i) = get_location( loc(i) ) 
   if (dovert) intvals(i) = query_location(loc, 'WHICH_VERT')
enddo

call nc_check(nf90_put_var(ncFileID, LocationVarId, locations, &
              start=(/ 1, startlocindex /), count=(/ LocationDims, loccount /) ), &
              'nc_write_multiple_locations', 'put_var:location')

if (present(WhichVertVarID)) then
   if (WhichVertVarID /= MISSING_I) then
      intval = query_location(loc, 'WHICH_VERT')
      call nc_check(nf90_put_var(ncFileID, WhichVertVarID, intval, &
                    start=(/ startlocindex /), count=(/ loccount /) ), &
                    'nc_write_multiple_locations','put_var:vert' )
   endif
endif

deallocate(locations)
if (dovert) deallocate(intvals)

end subroutine nc_write_multiple_locations

!----------------------------------------------------------------------------

end module location_utilitiess_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
