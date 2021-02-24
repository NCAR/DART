! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> common routines which can write netcdf arrays no matter which
!> location module is compiled in.

!>@todo FIXME  there is no documentation for this module.

!>@todo FIXME  define ALL the VERTISXXX here, add a call to has_vertical_choice
!>to the xxx/location_mod, and if no, return without needing routines
!>in each of the specific location routines.  if yes, call query and 
!>match it with the type to return true/false.  remove vert_is_xxx from 
!>all the other location routines!!!  :)
!>will have to replicate the VERTISxxx in the location modules because
!>fortran doesn't allow circular 'use's between modules.  ugh.

!>@todo There is no documentation for this module.

module location_io_mod

use            types_mod, only : r8, MISSING_I
use netcdf_utilities_mod, only : nc_check

!>@todo FIXME: should there be accessor functions for 5 LocationXXX variables below?
use        location_mod, only : location_type, get_location, &
                                LocationDims, LocationName, LocationLName, &
                                LocationStorageOrder, LocationUnits, &
                                has_vertical_choice, query_location


use typeSizes
use netcdf

implicit none
private

public :: nc_write_location_atts, nc_get_location_varids, &
          nc_write_location, nc_write_location_vert,      &
          nc_add_location_atts

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

!>@todo FIXME does the last arg need to be an optional actual dim?
!>or an additional dim?  check obs_seq_verify for usage

subroutine nc_write_location_atts(ncFileID, dimlen, use_dimID, fname) 
 
integer,                    intent(in) :: ncFileID    ! handle to the netcdf file
integer,                    intent(in) :: dimlen      ! number of locations to be created
integer,          optional, intent(in) :: use_dimID   ! if other than locations dim, use this
character(len=*), optional, intent(in) :: fname       ! file name (for error printing purposes)

integer :: LocDimID, LDimID, VarID
integer :: rc
character(len=*), parameter :: context = 'location_io_mod:nc_write_location_atts'

! get an id for the location dimension:  
!  if the user passes us in a dimension id, 
!     assume it is already created and use it.
!  if they pass in a valid length, create a 'location' dimension.
!  otherwise, use the unlimited dimension

if (present(use_dimID)) then
   LocDimID = use_dimID

else if (dimlen /= 0) then
   rc = nf90_def_dim(ncid=ncFileID, name='location', len=dimlen, dimid=LocDimID)
   call nc_check(rc, context, 'def_dim:location', fname)

else
   rc = nf90_inquire(ncFileID, UnlimitedDimID=LocDimID)
   call nc_check(rc, context, 'inquire:unlimited_dim', fname)
endif

! if there is more than a single number associated with each location
! define the dimension "across" the locations.  e.g. len 3 for 3d dimensions.
if (LocationDims > 1) then
   rc = nf90_def_dim(ncid=ncFileID, name='locdim', len=LocationDims, dimid=LDimID)
   call nc_check(rc, context, 'def_dim:locdim', fname)
endif

! Define the location variable and attributes

if (LocationDims > 1) then
   rc = nf90_def_var(ncFileID, 'location', xtype=nf90_double, &
                     dimids=(/ LDimID, LocDimID /), varid=VarID)
   call nc_check(rc, context, 'def_var2:location', fname)
else
   rc = nf90_def_var(ncFileID, 'location', xtype=nf90_double, &
                     dimids=(/ LocDimID /), varid=VarID)
   call nc_check(rc, context, 'def_var4:location', fname)
endif

call nc_add_location_atts(ncFileID, 'location', fname)

end subroutine nc_write_location_atts

!----------------------------------------------------------------------------
!> Add location-related attributes to the given variable.
!> (assumes that this variable will contain location info)

subroutine nc_add_location_atts(ncFileID, varname, fname) 
 
integer,                    intent(in) :: ncFileID    ! handle to the netcdf file
character(len=*),           intent(in) :: varname     ! name of location variable
character(len=*), optional, intent(in) :: fname       ! file name (for error printing purposes)

integer :: VarID
integer :: rc, ndims
character(len=*), parameter :: context = 'location_io_mod:nc_add_location_atts'

integer :: dimids(NF90_MAX_VAR_DIMS)

! find the id of the given variable name
rc = nf90_inq_varid(ncFileID, varname, varid=VarID)
call nc_check(rc, context, 'inq_varid:'//trim(varname), fname)

rc = nf90_put_att(ncFileID, VarID, 'description', 'location coordinates')
call nc_check(rc, context, 'put_att:description', fname)

rc = nf90_put_att(ncFileID, VarID, 'location_type', trim(LocationName))
call nc_check(rc, context, 'put_att:location_type', fname)

rc = nf90_put_att(ncFileID, VarID, 'long_name', trim(LocationLName))
call nc_check(rc, context, 'put_att:long_name', fname)

rc = nf90_put_att(ncFileID, VarID, 'storage_order', trim(LocationStorageOrder))
call nc_check(rc, context, 'put_att:storage_order', fname)

rc = nf90_put_att(ncFileID, VarID, 'units', trim(LocationUnits))
call nc_check(rc, context, 'put_att:units', fname)

! Some of the locations types need another variable to define
! the vertical coordinate system being used. If you are writing
! an array of locations, you also need to write an array of the vertical
! coordinate system. That array must be the same length as the location
! array being written.

if (has_vertical_choice()) then

   rc = nf90_inquire_variable(ncFileID, VarID, ndims=ndims, dimids=dimids)
   call nc_check(rc, context, 'inquire_variable:'//trim(varname), fname)
   call nc_write_location_vert(ncFileID, dimids(ndims), fname)

endif

end subroutine nc_add_location_atts

!----------------------------------------------------------------------------
!> Define the ancillary vertical array and attributes

subroutine nc_write_location_vert(ncFileID, use_dimID, fname)

integer,                     intent(in) :: ncFileID    ! handle to the netcdf file
integer,          optional,  intent(in) :: use_dimID   ! dimension to use 
character(len=*), optional,  intent(in) :: fname       ! file name (for printing purposes)

integer :: VarID, dimID, rc
character(len=*), parameter :: context = 'location_io_mod:nc_write_location_vert'

! if they give us the dimension the 'locations' array is using,
! use that.  otherwise default to the unlimited dimension for this variable.
if (present(use_dimID)) then
   dimID = use_dimID
else
   rc = nf90_inquire(ncFileID, UnlimitedDimID=dimID)
   call nc_check(rc, context, 'inquire:unlimited_dim', fname)
endif

rc = nf90_def_var(ncFileID, name='which_vert', xtype=nf90_int, &
                  dimids=(/ dimID /), varid=VarID)
call nc_check(rc, context, 'def_var:which_vert', fname)

rc = nf90_put_att(ncFileID, VarID, 'long_name', 'vertical coordinate system code')
call nc_check(rc, context, 'put_att:long_name', fname)

! these are the same across all location modules, for better
! or worse.  if we find location schemes which differ, this
! needs to come from a call to the location module proper.

rc = nf90_put_att(ncFileID, VarID, 'VERTISUNDEF', VERTISUNDEF)
call nc_check(rc, context, 'put_att:VERTISUNDEF', fname)

rc = nf90_put_att(ncFileID, VarID, 'VERTISSURFACE', VERTISSURFACE)
call nc_check(rc, context, 'put_att:VERTISSURFACE', fname)

rc = nf90_put_att(ncFileID, VarID, 'VERTISLEVEL', VERTISLEVEL)
call nc_check(rc, context, 'put_att:VERTISLEVEL', fname)

rc = nf90_put_att(ncFileID, VarID, 'VERTISPRESSURE', VERTISPRESSURE)
call nc_check(rc, context, 'put_att:VERTISPRESSURE', fname)

rc = nf90_put_att(ncFileID, VarID, 'VERTISHEIGHT', VERTISHEIGHT)
call nc_check(rc, context, 'put_att:VERTISHEIGHT', fname)

rc = nf90_put_att(ncFileID, VarID, 'VERTISSCALEHEIGHT', VERTISSCALEHEIGHT)
call nc_check(rc, context, 'put_att:VERTISSCALEHEIGHT', fname)

end subroutine nc_write_location_vert

!----------------------------------------------------------------------------
!> Return the LocationVarID and WhichVertVarID variables from a given netCDF file.
!>@todo FIXME do we still need this?
!>
!> ncFileId         the netcdf file descriptor
!> fname            the name of the netcdf file (for error messages only)
!> LocationVarID    the integer ID of the 'location' variable in the netCDF file
!> WhichVertVarID   the integer ID of the 'which_vert' variable in the netCDF file

subroutine nc_get_location_varids(ncFileID, LocationVarID, WhichVertVarID, fname)

integer,                    intent(in)  :: ncFileID   ! handle to the netcdf file
integer,                    intent(out) :: LocationVarID
integer,          optional, intent(out) :: WhichVertVarID
character(len=*), optional, intent(in)  :: fname      ! file name (for printing purposes)

integer :: rc
character(len=32) :: context = 'nc_write_location_varids'

rc = nf90_inq_varid(ncFileID, 'location', varid=LocationVarID)
call nc_check(rc, context, 'inq_varid:location ', fname)

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

subroutine nc_write_single_location(ncFileID, location, locindex, do_vert, fname)
 
integer,                    intent(in) :: ncFileID
type(location_type),        intent(in) :: location
integer,                    intent(in) :: locindex
logical, optional,          intent(in) :: do_vert
character(len=*), optional, intent(in) :: fname       ! file name (for error printing purposes)

integer :: LocationVarID
integer :: WhichVertVarID
real(r8), dimension(LocationDims) :: location_array
integer,  dimension(1) :: intval
logical :: write_vert
integer :: rc
character(len=32) :: context = 'nc_write_single_location'

! only let the user choose to write the vertical coordinate info
! if the location module has a vertical choice
write_vert = .false.
if (present(do_vert)) then
   if (has_vertical_choice()) write_vert = do_vert
endif

rc = nf90_inq_varid(ncFileID, 'location', varid=LocationVarID)
call nc_check(rc, context, 'location', fname)

location_array = get_location( location ) 

if (LocationDims > 1) then
   rc = nf90_put_var(ncFileID, LocationVarId, location_array, &
                     start=(/ 1, locindex /), count=(/ LocationDims, 1 /) )
else
   rc = nf90_put_var(ncFileID, LocationVarId, location_array, &
                     start=(/ locindex /), count=(/ 1 /) )
endif
call nc_check(rc, context, 'put_var:location', fname)

if (write_vert) then
   rc = nf90_inq_varid(ncFileID, 'which_vert', varid=WhichVertVarID)
   call nc_check(rc, context, 'inq_varid:which_vert', fname)

   intval = query_location(location, 'WHICH_VERT')
   rc = nf90_put_var(ncFileID, WhichVertVarID, intval, &
                     start=(/ locindex /), count=(/ 1 /) )
   call nc_check(rc, context, 'put_var:which_vert', fname)
endif

end subroutine nc_write_single_location

!----------------------------------------------------------------------------
!> Writes an array of locations to the specified netCDF variable and file.
!> The LocationVarID and WhichVertVarID must be the values returned from
!> the nc_get_location_varids call.

subroutine nc_write_multiple_locations(ncFileID, locations, loccount, startlocindex, do_vert, fname)

integer,             intent(in) :: ncFileID
type(location_type), intent(in) :: locations(:)
integer,             intent(in) :: loccount
integer, optional,   intent(in) :: startlocindex
logical, optional,   intent(in) :: do_vert
character(len=*), optional, intent(in) :: fname       ! file name (for error printing purposes)

integer :: LocationVarID
integer :: WhichVertVarID
real(r8), allocatable :: location_array(:,:)
integer,  allocatable :: intvals(:)
logical :: write_vert
integer :: rc, i, starthere
character(len=32) :: context = 'nc_write_single_location'

! only let the user choose to write the vertical coordinate info
! if the location module has a vertical choice
write_vert = .false.
if (present(do_vert)) then
   if (has_vertical_choice()) write_vert = do_vert
endif

starthere = 1
if (present(startlocindex)) starthere = startlocindex

rc = nf90_inq_varid(ncFileID, 'location', varid=LocationVarID)
call nc_check(rc, context, 'location', fname)

allocate(location_array(LocationDims,loccount))
if (write_vert) allocate(intvals(loccount))

do i=1, loccount
   location_array(:,i) = get_location( locations(i) ) 
   if (write_vert) intvals(i) = query_location(locations(i), 'WHICH_VERT')
enddo

if (LocationDims > 1) then
   rc = nf90_put_var(ncFileID, LocationVarId, location_array, &
                     start=(/ 1, starthere /), count=(/ LocationDims, loccount /) )
else
   rc = nf90_put_var(ncFileID, LocationVarId, location_array, &
                     start=(/ starthere /), count=(/ loccount /) )
endif
call nc_check(rc, context, 'put_var:location', fname)

if (write_vert) then
   rc = nf90_inq_varid(ncFileID, 'which_vert', varid=WhichVertVarID)
   call nc_check(rc, context, 'inq_varid:which_vert', fname)

   rc = nf90_put_var(ncFileID, WhichVertVarID, intvals, &
                     start=(/ starthere /), count=(/ loccount /) )
   call nc_check(rc, context, 'put_var:which_vert', fname)
endif

deallocate(location_array)
if (write_vert) deallocate(intvals)

end subroutine nc_write_multiple_locations

!----------------------------------------------------------------------------

end module location_io_mod

