! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! utility routines intended to make it easier to create the most
! common types of observations:  threed_sphere locations, which
! is what 99% of earth observations use.  

! added an alternative for obs locations specified in X,Y,Z cartesian
! coordinates.  used for subsurface water models, and some idealized
! atmosphere models.

module obs_utilities_mod


use        types_mod, only : i2, i4, r8, MISSING_R8, MISSING_I
use    utilities_mod, only : E_MSG, E_ERR, error_handler
use  netcdf_utilities_mod, only : nc_check
use      obs_def_mod, only : obs_def_type, set_obs_def_time, set_obs_def_type_of_obs, &
                             set_obs_def_error_variance, set_obs_def_location, &
                             get_obs_def_time, get_obs_def_location,           &
                             get_obs_def_type_of_obs, get_obs_def_error_variance,         &
                             set_obs_def_key
use obs_sequence_mod, only : obs_sequence_type, obs_type, insert_obs_in_seq, &
                             set_obs_values, set_qc, set_obs_def, get_qc,    &
                             get_obs_values, get_obs_def
use time_manager_mod, only : time_type, operator(>=), set_time, get_time
use     location_mod, only : set_location, location_type, get_location, &
                             query_location

use netcdf

implicit none
private

public :: create_3d_obs,    &
          query_3d_obs,     &
          add_obs_to_seq,   &
          getdimlen,        &
          getvarshape,      &
          getvar_real,      &
          getvar_int,       &
          getvar_char,      &
          get_2Dshort_as_r8,  &
          get_short_as_r8,  &
          get_int_as_r8,    &
          get_or_fill_real, &
          get_or_fill_int,  &
          get_or_fill_QC,   &
          getvar_real_2d,   &
          getvar_int_2d,    &
          getvar_real_1d_1val,     &
          getvar_int_1d_1val,      &
          getvar_real_2d_slice,    &
          get_or_fill_QC_2d_slice, &
          is_variable_integer,     &
          is_variable_real,        &
          query_varname,    &
          set_missing_name

!>@todo FIXME there is no documentation for this module

! module global storage - change with 'set_missing_name()' routine.
! if file is using a different name for this attribute. (nonstandard)
! note that '_FillValue' is also used for data that was never written,
! as opposed to data which is written but missing.  we don't care
! to distinguish between these cases - so what do we look for here?

character(len=NF90_MAX_NAME) :: missing_name = 'missing_value'

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'obs_utilities_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

! create an interface for both 3d sphere obs (most earth obs are this)
! and 3d cartesian obs (subsurface obs, idealized experiments)
interface create_3d_obs
   module procedure create_3d_sphere_obs
   module procedure create_3d_cartesian_obs
end interface

character(len=512) :: string1, string2
contains

!--------------------------------------------------------------------
!>  subroutine to create an observation type from observation data.  
!>
!>  NOTE: assumes the code is using the threed_sphere locations module, 
!>        that the observation has a single data value and a single
!>        qc value. there is a last optional argument now; if the obs
!>        has additional metadata associated with it, the calling code
!>        should get a unique 'key' from a set metadata routine, and
!>        pass that in as the last argument to this routine.  if present
!>        it will be set in the obs_def derived type.
!>
!>  lat   - latitude of observation
!>  lon   - longitude of observation
!>  vval  - vertical coordinate
!>  vkind - kind of vertical coordinate (pressure, level, etc)
!>  obsv  - observation value
!>  okind - observation kind
!>  oerr  - standard deviation of observation error
!>  day   - gregorian day
!>  sec   - gregorian second
!>  qc    - quality control value
!>  obs   - observation type
!>  key   - optional type-specific integer key from a set_metadata() routine
!>
!>  created Oct. 2007 Ryan Torn, NCAR/MMM
!>  adapted for more generic use 11 Mar 2010, nancy collins, ncar/image
!>  added optional metadata key   8 Nov 2013, nancy collins, ncar/image
!>  overloaded to allow 3d cartesian locations, 23 sept 2019, nsc

subroutine create_3d_sphere_obs(lat, lon, vval, vkind, obsv, okind, oerr, day, sec, qc, &
                         obs, key)

integer,           intent(in)    :: okind, vkind, day, sec
real(r8),          intent(in)    :: lat, lon, vval, obsv, oerr, qc
type(obs_type),    intent(inout) :: obs
integer, optional, intent(in)    :: key

real(r8) :: loc_info(4)
type(location_type) :: loc

! setting location info with an array will allow this code
! to compile with any location type.  this interface is only
! intended to work with threed_sphere locations.
loc_info(1) = lon
loc_info(2) = lat
loc_info(3) = vval
loc_info(4) = real(vkind)

loc = set_location(loc_info)
call create_obs(loc, obsv, okind, oerr, day, sec, qc, obs, key)

end subroutine create_3d_sphere_obs

!--------------------------------------------------------------------
!>  subroutine to create an observation type from observation data.  
!>
!>  NOTE: assumes the code is using the 3d cartesian locations module, 
!>        that the observation has a single data value and a single
!>        qc value. there is a last optional argument now; if the obs
!>        has additional metadata associated with it, the calling code
!>        should get a unique 'key' from a set metadata routine, and
!>        pass that in as the last argument to this routine.  if present
!>        it will be set in the obs_def derived type.
!>
!>  X, Y, Z - location of obs
!>  obsv  - observation value
!>  okind - observation kind
!>  oerr  - standard deviation of observation error
!>  day   - gregorian day
!>  sec   - gregorian second
!>  qc    - quality control value
!>  obs   - observation type
!>  key   - optional type-specific integer key from a set_metadata() routine
!>
!>  overloaded to allow 3d cartesian locations, 23 sept 2019, nsc

subroutine create_3d_cartesian_obs(X,Y,Z, obsv, okind, oerr, day, sec, qc, &
                         obs, key)

integer,           intent(in)    :: okind, day, sec
real(r8),          intent(in)    :: X,Y,Z, obsv, oerr, qc
type(obs_type),    intent(inout) :: obs
integer, optional, intent(in)    :: key

real(r8) :: loc_info(3)
type(location_type) :: loc

! setting location info with an array will allow this code
! to compile with any location type.  this interface is only
! intended to work with threed_cartesian locations.
loc_info(1) = X
loc_info(2) = Y
loc_info(3) = Z

loc = set_location(loc_info)
call create_obs(loc, obsv, okind, oerr, day, sec, qc, obs, key)

end subroutine create_3d_cartesian_obs

!--------------------------------------------------------------------
!>  private subroutine to create an observation type from observation data.  
!>
!>  loc   - observation location 
!>  obsv  - observation value
!>  okind - observation kind
!>  oerr  - standard deviation of observation error
!>  day   - gregorian day
!>  sec   - gregorian second
!>  qc    - quality control value
!>  obs   - observation type
!>  key   - optional type-specific integer key from a set_metadata() routine
!>

subroutine create_obs(loc, obsv, okind, oerr, day, sec, qc, obs, key)

type(location_type), intent(in)  :: loc
integer,           intent(in)    :: okind, day, sec
real(r8),          intent(in)    :: obsv, oerr, qc
type(obs_type),    intent(inout) :: obs
integer, optional, intent(in)    :: key

real(r8)              :: obs_val(1), qc_val(1)
type(obs_def_type)    :: obs_def

! doing it with an array will work for any location mod.
call set_obs_def_location(obs_def, loc)

call set_obs_def_type_of_obs(obs_def, okind)
call set_obs_def_time(obs_def, set_time(sec, day))
call set_obs_def_error_variance(obs_def, oerr * oerr)
if (present(key)) then
   call set_obs_def_key(obs_def, key)
endif
call set_obs_def(obs, obs_def)

obs_val(1) = obsv
call set_obs_values(obs, obs_val)
qc_val(1)  = qc
call set_qc(obs, qc_val)

end subroutine create_obs

!--------------------------------------------------------------------
!> takes a DART observation type and extracts the constituent pieces
!>
!> NOTE: assumes the code is using the threed_sphere locations module, 
!>       that the observation has a single data value and a single
!>       qc value, and that this obs type has no additional required
!>       data (e.g. gps and radar obs need additional data per obs)
!>
!> obs   - observation type (in)
!> lat   - latitude of observation (all the rest out)
!> lon   - longitude of observation
!> vval  - vertical coordinate
!> vkind - kind of vertical coordinate (pressure, level, etc)
!> obsv  - observation value
!> okind - observation kind
!> oerr  - observation error
!> day   - gregorian day
!> sec   - gregorian second
!> qc    - quality control value
!>
!> created Apr. 2010, nancy collins, ncar/image


subroutine query_3d_obs(obs, lat, lon, vval, vkind, obsv, okind, oerr, day, sec, qc)

 type(obs_type), intent(in)  :: obs
 integer,        intent(out) :: okind, vkind, day, sec
 real(r8),       intent(out) :: lat, lon, vval, obsv, oerr, qc

real(r8)              :: obs_val(1), qc_val(1), locvals(3)
type(obs_def_type)    :: obs_def
type(location_type)   :: location

! get the obs_def out of the obs.  it has the location, time, kind;
! basically everything except the actual obs value and qc.
call get_obs_def(obs, obs_def)

location = get_obs_def_location(obs_def)
locvals = get_location(location)
lon  = locvals(1)
lat  = locvals(2)
vval = locvals(3)
vkind = query_location(location)

okind = get_obs_def_type_of_obs(obs_def)
call get_time(get_obs_def_time(obs_def), sec, day)

oerr = sqrt(get_obs_def_error_variance(obs_def))

! now get the value and qc from the obs_type.
call get_obs_values(obs, obs_val)
obsv = obs_val(1) 
call get_qc(obs, qc_val)
qc = qc_val(1)

end subroutine query_3d_obs


!--------------------------------------------------------------------
!> adds an observation to a sequence.  inserts if first obs, inserts 
!> with a prev obs to save searching if that's possible.
!>
!> seq - observation sequence to add obs to
!> obs - observation, already filled in, ready to add
!> obs_time - time of this observation, in dart time_type format
!> prev_obs - the previous observation that was added to this sequence
!>            (will be updated by this routine)
!> prev_time - the time of the previously added observation (will also
!>            be updated by this routine)
!> first_obs - should be initialized to be .true., and then will be
!>            updated by this routine to be .false. after the first obs
!>            has been added to this sequence.
!>
!> created Mar 8, 2010   nancy collins, ncar/image

subroutine add_obs_to_seq(seq, obs, obs_time, prev_obs, prev_time, first_obs)

  type(obs_sequence_type), intent(inout) :: seq
  type(obs_type),          intent(inout) :: obs, prev_obs
  type(time_type),         intent(in)    :: obs_time
  type(time_type),         intent(inout) :: prev_time
  logical,                 intent(inout) :: first_obs

! insert(seq,obs) always works (i.e. it inserts the obs in
! proper time format) but it can be slow with a long file.
! supplying a previous observation that is older (or the same
! time) as the new one speeds up the searching a lot.

if(first_obs) then    ! for the first observation, no prev_obs
   call insert_obs_in_seq(seq, obs)
   first_obs = .false.
else               
   if(obs_time >= prev_time) then  ! same time or later than previous obs
      call insert_obs_in_seq(seq, obs, prev_obs)
   else                            ! earlier, search from start of seq
      call insert_obs_in_seq(seq, obs)
   endif
endif

! update for next time
prev_obs = obs
prev_time = obs_time

end subroutine add_obs_to_seq


!--------------------------------------------------------------------
!> given a netCDF file handle and a dimension name, return the dimension size
!>
!> ncid - open netcdf file handle
!> dimname - string name of netcdf dimension
!> dout - output value.  integer
!>
!> created 11 Mar 2010,  nancy collins,  ncar/image

subroutine getdimlen(ncid, dimname, dout)

 integer,            intent(in)   :: ncid
 character(len = *), intent(in)   :: dimname
 integer,            intent(out)  :: dout

integer :: dimid

! get the requested dimension size
call nc_check( nf90_inq_dimid(ncid, dimname, dimid), &
               'getdimlen', 'inq dimid '//trim(dimname))
call nc_check( nf90_inquire_dimension(ncid, dimid, len=dout), &
               'getdimlen', 'inquire dimension '//trim(dimname))

end subroutine getdimlen


!--------------------------------------------------------------------
!> subroutine that returns an array of the size of the variable
!>
!> ncid ....... open netcdf file handle
!> varname .... string name of netcdf variable
!> numdims .... the rank of the variable (untested on scalars)
!> dimlengths . an array specifying the length of each dimension
!>
!> created 4 Mar 2015,  tim hoar,  ncar/image

subroutine getvarshape(ncid, varname, numdims, dimlengths)

integer,          intent(in)   :: ncid
character(len=*), intent(in)   :: varname
integer,          intent(out)  :: numdims
integer,          intent(out)  :: dimlengths(NF90_MAX_VAR_DIMS)

integer :: varid, dimid
integer :: dimIDs(NF90_MAX_VAR_DIMS)

character(len=512) :: string1

call nc_check(nf90_inq_varid(ncid, varname, varid), &
              'getvarshape', 'inq_varid ['//trim(varname)//']')

call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimIDs,ndims=numdims), &
            'getvarshape', 'inquire_variable '//trim(varname))

do dimid = 1,numdims

   write(string1,'(''inquire_dimension'',i2,1x,A)') dimid,trim(varname)

   call nc_check(nf90_inquire_dimension(ncid, dimIDs(dimid), &
                  len=dimlengths(dimid)), 'getvarshape', string1)

enddo

end subroutine getvarshape


!--------------------------------------------------------------------
!> sets the name of the attribute that describes missing values.  
!>
!> the suggested rule is data that was never written to the netcdf
!> variable will have a value equal to "_FillValue", and this should
!> be added as an attribute on the variable when it is created.
!> the "missing_data" attribute indicates data that was written
!> but is not valid or present, e.g. land grid points in an ocean
!> data field.
!> 
!> but in practice there are problems.  one is that users may not
!> add any attributes and have a private convention for missing or
!> invalid data (e.g. 99999 or some such pattern).
!>
!> the netcdf files can be created with an option to NOT prefill
!> variables with the _FillValue.  this speeds execution when the
!> code expects to fully write all variables. but if there are
!> paths through the code where a variable isn't written by mistake,
!> the contents are undefined.
!>
!> and finally, some users mix 'missing_value' and "_FillValue" usage,
!> writing the value of "_FillValue" to indicate missing data. 
!> 
!> in this code all we care about is knowing what data value
!> indicates data is missing for any reason.  if the code uses 
!> a different attribute value, set it with this routine.
!>
!> name - string name of attribute that holds the missing value
!>
!> created 11 Mar 2010,  nancy collins,  ncar/image
!> (this usage comment updated 26 sep 2019, nsc)

subroutine set_missing_name(name)

 character(len = *), intent(in)   :: name

   if (len(name) > NF90_MAX_NAME) then
   write(string1,*) 'name must be less than ', NF90_MAX_NAME, ' chars long'
   write(string2,*) 'name is ', len(name), ' long'
   call error_handler(E_ERR,'set_missing_name:',string1,source,text2=string2)
   endif

   missing_name = name

end subroutine set_missing_name


!--------------------------------------------------------------------
!>inquires, gets the variable, and fills 
!>            in the missing value attribute if that arg is present.
!>            gets the entire array, no start or count specified.
!>
!>      ncid - open netcdf file handle
!>      varname - string name of netcdf variable
!>      darray - output array.  real(r8)
!>      dmiss - value that signals a missing value   real(r8), optional
!>
!>     created 11 Mar 2010,  nancy collins,  ncar/image

subroutine getvar_real(ncid, varname, darray, dmiss)

 integer,            intent(in)   :: ncid
 character(len = *), intent(in)   :: varname
 real(r8),           intent(out)  :: darray(:)
 real(r8), optional, intent(out)  :: dmiss

integer  :: varid, nfrc
real(r8) :: fill, miss
logical  :: found_miss

! read the data for the requested array, and optionally get the fill value
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'getvar_real', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, darray), &
               'getvar_real', 'getting var '// trim(varname))

! the logic here is: 
!  if the user hasn't asked about the fill value, don't look for any of
!  these attributes and just return.

!  if the user has told us another attribute name to look for, try that first.
!  it's currently NOT an error if it's not there.   then second, look for the
!  official '_FillValue' attr.  if it's found, return it as the missing value.
! 
!  if there are both, overwrite the missing value with the fill value and return
!  the fill value as the 'dmiss' argument.
!
!  if neither are there, set dmiss to the DART missing value.  this could also
!  be an error, but default to being permissive for now.

if (present(dmiss)) then
   dmiss = MISSING_R8
   found_miss = .false.

   ! if user defined another attribute name for missing vals
   ! look for it first, and make it an error if it's not there?
   if (missing_name /= '') then
      nfrc = nf90_get_att(ncid, varid, missing_name, miss)
      if (nfrc == NF90_NOERR) then 
         found_miss = .true.
         dmiss = miss
      endif
   endif

      ! this would make it a fatal error if not found
      !call nc_check( nf90_get_att(ncid, varid, missing_name', miss), &
      !         'getvar_real', 'getting attr "//trim(missing_name)//" for '//trim(varname))

   ! the predefined netcdf fill value.
   nfrc = nf90_get_att(ncid, varid, '_FillValue', fill)
   if (nfrc == NF90_NOERR) then
      if (.not. found_miss) then  
         found_miss = .true.
         dmiss = fill
      else
         ! found both, set all to fill value
         where(darray .eq. miss) darray = fill 
         dmiss = fill
      endif
   endif

   ! if you wanted an error if you specified dmiss and neither attr are
   ! there, you'd test found_miss here.  if it's still false, none were there.

endif
  
end subroutine getvar_real


!--------------------------------------------------------------------
!>   get_short_as_r8 - subroutine that inquires, gets the variable, and fills 
!>            in the missing value attribute if that arg is present.
!>            gets the entire array, no start or count specified.
!>
!> FIXME: this code handles scale and offset ok, but like most of the other
!> routines here it doesn't check for 'missing_value' which many obs files
!> have either in addition to the _FillValue or instead of it.
!> e.g. a dump from a MADIS mesonet file:
!>
!> double observationTime ( recNum )
!>  long_name : time of observation
!>  units : seconds since 1970-1-1 00:00:00.0
!>  _FillValue : 3.40282346e+38
!>  missing_value : -9999
!>  standard_name : time 
!> 
!> i've most commonly encountered the "missing_value" in the actual data 
!> for obs files, not the fill value.  we need a routine that looks for missing, 
!> then fill, and decides what to do if there are both.  or use the routine
!> that was already here 'set_missing_value' and the calling code tells us
!> which attribute name this particular netcdf file is using.
!>
!>  ncid    - open netcdf file handle
!>  varname - string name of netcdf variable
!>  darray  - output array.  real(r8)
!>  dmiss   - value that signals a missing value   real(r8), optional

subroutine get_short_as_r8(ncid, varname, darray, dmiss)

integer,            intent(in)  :: ncid
character(len=*),   intent(in)  :: varname
real(r8),           intent(out) :: darray(:)
real(r8), optional, intent(out) :: dmiss

integer     :: varid
integer(i2) :: shortarray(size(darray))
integer(i2) :: FillValue
real(r8)    :: scale_factor, add_offset
integer     :: offset_exists, scaling_exists, fill_exists
!>@todo FIXME need missing_exists or something

! read the data for the requested array, and get the fill value
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'get_short_as_r8', 'inquire var "'// trim(varname)//'"')

 offset_exists = nf90_get_att(ncid, varid, 'add_offset',   add_offset)
scaling_exists = nf90_get_att(ncid, varid, 'scale_factor', scale_factor)
   fill_exists = nf90_get_att(ncid, varid, '_FillValue',   FillValue)
!  miss_exists = nf90_get_att(ncid, varid, 'missing_value', miss_value)

call nc_check( nf90_get_var(ncid, varid, shortarray), &
               'get_short_as_r8', 'getting var '// trim(varname))

darray = real(shortarray,r8)


if (fill_exists == NF90_NOERR) then ! FillValue exists

   if ( (offset_exists == NF90_NOERR) .and. (scaling_exists == NF90_NOERR) ) then
      where (shortarray /= FillValue) darray = darray * scale_factor + add_offset
   elseif (offset_exists == NF90_NOERR) then
      where (darray /= FillValue) darray = darray * scale_factor
   elseif (scaling_exists == NF90_NOERR) then
      where (darray /= FillValue) darray = darray + add_offset
   endif

   if (present(dmiss)) dmiss = real(FillValue,r8)

else

   if ( (offset_exists == NF90_NOERR) .and. (scaling_exists == NF90_NOERR) ) then
      darray = darray * scale_factor + add_offset
   elseif (offset_exists == NF90_NOERR) then
      darray = darray * scale_factor
   elseif (scaling_exists == NF90_NOERR) then
      darray = darray + add_offset
   endif

   if (present(dmiss)) dmiss = MISSING_R8

endif

end subroutine get_short_as_r8


!--------------------------------------------------------------------
!>   get_2Dshort_as_r8 - subroutine that inquires, gets the variable, and fills 
!>            in the missing value attribute if that arg is present.
!>            gets the entire array, no start or count specified.
!>
!> FIXME: this code handles scale and offset ok, but like most of the other
!> routines here it doesn't check for 'missing_value' which many obs files
!> have either in addition to the _FillValue or instead of it.
!> e.g. a dump from a MADIS mesonet file:
!>
!> double observationTime ( recNum )
!>  long_name : time of observation
!>  units : seconds since 1970-1-1 00:00:00.0
!>  _FillValue : 3.40282346e+38
!>  missing_value : -9999
!>  standard_name : time 
!> 
!> i've most commonly encountered the "missing_value" in the actual data 
!> for obs files, not the fill value.  we need a routine that looks for missing, 
!> then fill, and decides what to do if there are both.  or use the routine
!> that was already here 'set_missing_value' and the calling code tells us
!> which attribute name this particular netcdf file is using.
!>
!>  ncid    - open netcdf file handle
!>  varname - string name of netcdf variable
!>  darray  - output array.  real(r8)
!>  dmiss   - value that signals a missing value   real(r8), optional

subroutine get_2Dshort_as_r8(ncid, varname, darray, dmiss)

integer,            intent(in)  :: ncid
character(len=*),   intent(in)  :: varname
real(r8),           intent(out) :: darray(:,:)
real(r8), optional, intent(out) :: dmiss

integer     :: varid
integer(i2), allocatable :: shortarray(:,:)
integer(i2) :: FillValue
real(r8)    :: scale_factor, add_offset
integer     :: offset_exists, scaling_exists, fill_exists
!>@todo FIXME need missing_exists or something

allocate(shortarray(size(darray,1),size(darray,2)))

! read the data for the requested array, and get the fill value
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'get_2Dshort_as_r8', 'inquire var "'// trim(varname)//'"')

 offset_exists = nf90_get_att(ncid, varid, 'add_offset',   add_offset)
scaling_exists = nf90_get_att(ncid, varid, 'scale_factor', scale_factor)
   fill_exists = nf90_get_att(ncid, varid, '_FillValue',   FillValue)
!  miss_exists = nf90_get_att(ncid, varid, 'missing_value', miss_value)

call nc_check( nf90_get_var(ncid, varid, shortarray), &
               'get_2Dshort_as_r8', 'getting var '// trim(varname))

darray = real(shortarray,r8)

if (fill_exists == NF90_NOERR) then ! FillValue exists

   if ( (offset_exists == NF90_NOERR) .and. (scaling_exists == NF90_NOERR) ) then
      where (shortarray /= FillValue) darray = darray * scale_factor + add_offset
   elseif (offset_exists == NF90_NOERR) then
      where (darray /= FillValue) darray = darray * scale_factor
   elseif (scaling_exists == NF90_NOERR) then
      where (darray /= FillValue) darray = darray + add_offset
   endif

   if (present(dmiss)) dmiss = real(FillValue,r8)

else

   if ( (offset_exists == NF90_NOERR) .and. (scaling_exists == NF90_NOERR) ) then
      darray = darray * scale_factor + add_offset
   elseif (offset_exists == NF90_NOERR) then
      darray = darray * scale_factor
   elseif (scaling_exists == NF90_NOERR) then
      darray = darray + add_offset
   endif

   if (present(dmiss)) dmiss = MISSING_R8

endif

deallocate(shortarray)

end subroutine get_2Dshort_as_r8


!--------------------------------------------------------------------
!>   get_int_as_r8 - subroutine that inquires, gets the variable, and fills 
!>            in the missing value attribute if that arg is present.
!>            gets the entire array, no start or count specified.
!>
!>  ncid    - open netcdf file handle
!>  varname - string name of netcdf variable
!>  darray  - output array.  real(r8)
!>  dmiss   - value that signals a missing value   real(r8), optional

subroutine get_int_as_r8(ncid, varname, darray, dmiss)

integer,            intent(in)  :: ncid
character(len=*),   intent(in)  :: varname
real(r8),           intent(out) :: darray(:)
real(r8), optional, intent(out) :: dmiss

integer     :: varid
integer(i4) :: shortarray(size(darray))
integer(i4) :: FillValue
real(r8)    :: scale_factor, add_offset
integer     :: offset_exists, scaling_exists, fill_exists

! read the data for the requested array, and get the fill value
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'get_int_as_r8', 'inquire var "'// trim(varname)//'"')

 offset_exists = nf90_get_att(ncid, varid, 'add_offset',   add_offset)
scaling_exists = nf90_get_att(ncid, varid, 'scale_factor', scale_factor)
   fill_exists = nf90_get_att(ncid, varid, '_FillValue',   FillValue)

call nc_check( nf90_get_var(ncid, varid, shortarray), &
               'get_int_as_r8', 'getting var '// trim(varname))

darray = real(shortarray,r8)

if (fill_exists == NF90_NOERR) then ! FillValue exists

   if ( (offset_exists == NF90_NOERR) .and. (scaling_exists == NF90_NOERR) ) then
      where (shortarray /= FillValue) darray = darray * scale_factor + add_offset
   elseif (offset_exists == NF90_NOERR) then
      where (darray /= FillValue) darray = darray * scale_factor
   elseif (scaling_exists == NF90_NOERR) then
      where (darray /= FillValue) darray = darray + add_offset
   endif

   if (present(dmiss)) dmiss = real(FillValue,r8)

else

   if ( (offset_exists == NF90_NOERR) .and. (scaling_exists == NF90_NOERR) ) then
      darray = darray * scale_factor + add_offset
   elseif (offset_exists == NF90_NOERR) then
      darray = darray * scale_factor
   elseif (scaling_exists == NF90_NOERR) then
      darray = darray + add_offset
   endif

   if (present(dmiss)) dmiss = MISSING_R8

endif

end subroutine get_int_as_r8


!--------------------------------------------------------------------
!>   getvar_int - subroutine that inquires, gets the variable, and fills 
!>            in the missing value attribute if that arg is present.
!>            gets the entire array, no start or count specified.
!>
!>      ncid - open netcdf file handle
!>      varname - string name of netcdf variable
!>      darray - output array.  integer
!>      dmiss - value that signals a missing value   integer, optional
!>
!>     created 11 Mar 2010,  nancy collins,  ncar/image

subroutine getvar_int(ncid, varname, darray, dmiss)

 integer,            intent(in)   :: ncid
 character(len = *), intent(in)   :: varname
 integer,            intent(out)  :: darray(:)
 integer,  optional, intent(out)  :: dmiss

integer  :: varid, nfrc
integer  :: fill, miss
logical  :: found_miss

! read the data for the requested array, and get the fill value
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'getvar_int', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, darray), &
               'getvar_int', 'getting var '// trim(varname))

! see long comment in getvar_real() about the logic here.
if (present(dmiss)) then
   dmiss = MISSING_I
   found_miss = .false.

   ! if user defined another attribute name for missing vals
   ! look for it first, and make it an error if it's not there?
   if (missing_name /= '') then
      nfrc = nf90_get_att(ncid, varid, missing_name, miss)
      if (nfrc == NF90_NOERR) then 
         found_miss = .true.
         dmiss = miss
      endif
   endif

      ! this would make it a fatal error if not found
      !call nc_check( nf90_get_att(ncid, varid, missing_name', miss), &
      !         'getvar_real', 'getting attr "//trim(missing_name)//" for '//trim(varname))

   ! the predefined netcdf fill value.
   nfrc = nf90_get_att(ncid, varid, '_FillValue', fill)
   if (nfrc == NF90_NOERR) then
      if (.not. found_miss) then  
         found_miss = .true.
         dmiss = fill
      else
         ! found both, set all to fill value
         where(darray .eq. miss) darray = fill 
         dmiss = fill
      endif
   endif

   ! if you wanted an error if you specified dmiss and neither attr are
   ! there, you'd test found_miss here.  if it's still false, none were there.

endif
  
end subroutine getvar_int


!--------------------------------------------------------------------
!>   getvar_char - subroutine that inquires, gets the variable, and fills 
!>            in the missing value attribute if that arg is present.
!>            gets the entire array, no start or count specified.
!>
!>      ncid - open netcdf file handle
!>      varname - string name of netcdf variable
!>      darray - output array.  character
!>      dmiss - value that signals a missing value  char, optional
!>
!>     created 11 Mar 2010,  nancy collins,  ncar/image

subroutine getvar_char(ncid, varname, darray, dmiss)

 integer,            intent(in)   :: ncid
 character(len = *), intent(in)   :: varname
 character(len = *), intent(out)  :: darray(:)
 character, optional, intent(out)  :: dmiss

integer   :: varid, nfrc
character :: fill, miss
logical   :: found_miss

! read the data for the requested array, and get the fill value
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'getvar_char', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, darray), &
               'getvar_char', 'getting var '// trim(varname))

! see long comment in getvar_real() about the logic here.
if (present(dmiss)) then
   dmiss = ''
   found_miss = .false.

   ! if user defined another attribute name for missing vals
   ! look for it first, and make it an error if it's not there?
   if (missing_name /= '') then
      nfrc = nf90_get_att(ncid, varid, missing_name, miss)
      if (nfrc == NF90_NOERR) then 
         found_miss = .true.
         dmiss = miss
      endif
   endif

   ! the predefined netcdf fill value.
   nfrc = nf90_get_att(ncid, varid, '_FillValue', fill)
   if (nfrc == NF90_NOERR) then
      if (.not. found_miss) then  
         found_miss = .true.
         dmiss = fill
      else
         ! found both, set all to fill value
         where(darray .eq. miss) darray = fill 
         dmiss = fill
      endif
   endif

   ! if you wanted an error if you specified dmiss and neither attr are
   ! there, you'd test found_miss here.  if it's still false, none were there.

endif
  
end subroutine getvar_char


!--------------------------------------------------------------------
!>   get_or_fill_real - subroutine which gets the requested netcdf variable
!>           but if it isn't there, it fills the array with 0s.  not an
!>           error if it's not present.  assumes real(r8) data array
!>
!>      ncid - open netcdf file handle
!>      varname - string name of netcdf variable
!>      darray - output array.  real(r8)
!>      fillval - real(r8) to fill with if not present.  optional.
!>      did_fill - reports back if fill was done.  logical, optional.
!>
!>     created Mar 8, 2010    nancy collins, ncar/image

subroutine get_or_fill_real(ncid, varname, darray, fillval, did_fill)

 integer,             intent(in)    :: ncid
 character(len = *),  intent(in)    :: varname
 real(r8),            intent(inout) :: darray(:)
 real(r8), optional,  intent(in)    :: fillval
 logical,  optional,  intent(out)   :: did_fill

integer :: varid, nfrc

! test to see if variable is present.  if yes, read it in.
! otherwise, set to fill value, or 0 if none given.

nfrc = nf90_inq_varid(ncid, varname, varid) 
if (nfrc == NF90_NOERR) then
   call nc_check( nf90_get_var(ncid, varid, darray), &
                  'get_or_fill_real', 'reading '//trim(varname) )
   if (present(did_fill)) did_fill = .false.
else
   if (present(fillval)) then
      darray = fillval
   else
      darray = 0.0_r8
   endif
   if (present(did_fill)) did_fill = .true.
endif
   
end subroutine get_or_fill_real


!--------------------------------------------------------------------
!>   get_or_fill_int - subroutine which gets the requested netcdf variable
!>           but if it isn't there, it fills the array with 0s.  not an
!>           error if it's not present.  assumes integer data array
!>
!>      ncid - open netcdf file handle
!>      varname - string name of netcdf variable
!>      darray - output array.  integer
!>      fillval - integer to fill with if not present.  optional.
!>      did_fill - reports back if fill was done.  logical, optional.
!>
!>     created Mar 8, 2010    nancy collins, ncar/image

subroutine get_or_fill_int(ncid, varname, darray, fillval, did_fill)

 integer,            intent(in)    :: ncid
 character(len = *), intent(in)    :: varname
 integer,            intent(inout) :: darray(:)
 integer, optional,  intent(in)    :: fillval
 logical, optional,  intent(out)   :: did_fill

integer :: varid, nfrc

! test to see if variable is present.  if yes, read it in.
! otherwise, set to fill value, or 0 if none given.

nfrc = nf90_inq_varid(ncid, varname, varid) 
if (nfrc == NF90_NOERR) then
   call nc_check( nf90_get_var(ncid, varid, darray), &
                  'get_or_fill_int', 'reading '//trim(varname) )
   if (present(did_fill)) did_fill = .false.
else
   if (present(fillval)) then
      darray = fillval
   else
      darray = 0
   endif
   if (present(did_fill)) did_fill = .true.
endif
   
end subroutine get_or_fill_int


!--------------------------------------------------------------------
!> get_or_fill_QC - get the expected QC netcdf variable.  if not found,
!>                  fill the array with 0s, and print out a message warning
!>                  that the input didn't have the expected array.
!>   
!>      minor tweak on the standard get_or_fill routine.  prints out a
!>      warning message, and don't return to the caller which one it did.
!>      would be more general if it didn't use 'QC' in the output string.
!>      if anyone else wants to use it, we could change 'QC' to 'warn' in
!>      the name, and remove the 'QC' from the output string, and make it
!>      slightly more general
!>
!>      ncid - open netcdf file handle
!>      varname - string name of netcdf variable
!>      darray - output array.  integer
!>
!>      created 17 mar 2010  nsc  ncar/image

subroutine get_or_fill_QC(ncid, varname, darray)

 integer,            intent(in)    :: ncid
 character(len = *), intent(in)    :: varname
 integer,            intent(inout) :: darray(:)

logical :: did_fill


call get_or_fill_int(ncid, varname, darray, 0, did_fill)

if (did_fill) &
  print *, 'QC field named ' // trim(varname) // ' was not found in input, 0 used instead'

end subroutine get_or_fill_QC


!--------------------------------------------------------------------
!>   getvar_real_2d - subroutine that inquires, gets the variable, and fills 
!>            in the missing value attribute if that arg is present.
!>            gets the entire array, no start or count specified.
!>            this version assumes you are reading an entire 2d array.
!>            see the slice versions for reading 1d from a 2d array.
!>
!>      ncid - open netcdf file handle
!>      varname - string name of netcdf variable
!>      darray - 2d output array.  real(r8)
!>      dmiss - value that signals a missing value   real(r8), optional
!>
!>     created 11 Mar 2010,  nancy collins,  ncar/image
!>     extended 1 May 2020 to support scale/offset

!> Repeated from the Unidata netCDF convention document:
!> scale_factor
!>
!> If present for a variable, the data are to be multiplied by this factor after the data
!> are read by the application that accesses the data.
!>
!> If valid values are specified using the valid_min, valid_max, valid_range, 
!> or _FillValue attributes, those values should be specified in the domain of the 
!> data in the file (the packed data), so that they can be interpreted before the 
!> scale_factor and add_offset are applied.
!>
!> add_offset
!>
!> If present for a variable, this number is to be added to the data after it is read 
!> by the application that accesses the data. If both scale_factor and add_offset 
!> attributes are present, the data are first scaled before the offset is added. 
!> The attributes scale_factor and add_offset can be used together to provide simple 
!> data compression to store low-resolution floating-point data as small integers in 
!> a netCDF dataset. When scaled data are written, the application should first subtract 
!> the offset and then divide by the scale factor, rounding the result to the nearest 
!> integer to avoid a bias caused by truncation towards zero.
!>
!> When scale_factor and add_offset are used for packing, the associated variable 
!> (containing the packed data) is typically of type byte or short, whereas the unpacked 
!> values are intended to be of type float or double. The attributes scale_factor and 
!> add_offset should both be of the type intended for the unpacked data, e.g. float or double.

subroutine getvar_real_2d(ncid, varname, darray, dmiss, &
   nc_start, nc_count, nc_stride, nc_map)

 integer,            intent(in)   :: ncid
 character(len = *), intent(in)   :: varname
 real(r8),           intent(out)  :: darray(:,:)
real(r8),         intent(out), optional :: dmiss
integer,          intent(in),  optional :: nc_start(:)
integer,          intent(in),  optional :: nc_count(:)
integer,          intent(in),  optional :: nc_stride(:)
integer,          intent(in),  optional :: nc_map(:)

integer  :: io, varid
integer  :: offset_exists, scaling_exists
integer  :: user_missing_exists, missing_exists, fill_exists
real(r8) :: user_missing_value, missing_value, FillValue
real(r8) :: scale_factor, add_offset

! read the data for the requested array

io = nf90_inq_varid(ncid, varname, varid)
call nc_check( io, 'getvar_real_2d', 'inquire var "'//trim(varname)//'"')

io = nf90_get_var(ncid, varid, darray, nc_start, nc_count, nc_stride, nc_map)
call nc_check( io, 'getvar_real_2d', 'getting var "'//trim(varname)//'"')

! apply the variable attributes

 offset_exists = nf90_get_att(ncid, varid, 'add_offset',    add_offset)
scaling_exists = nf90_get_att(ncid, varid, 'scale_factor',  scale_factor)
missing_exists = nf90_get_att(ncid, varid, 'missing_value', missing_value)
   fill_exists = nf90_get_att(ncid, varid, '_FillValue',    FillValue)

! the logic here is: 
!  if the user provides a missing value attribute name (via set_missing_name) , use it.
!  if the user does not provide  value, 'normal' stragegy is applied
!  This presents a problem when we go to write these variables, but since
!  this module does not write netCDF variables, ignore for now.

if (present(dmiss)) then
   if (missing_name /= '') then
      user_missing_exists = nf90_get_att(ncid, varid, missing_name, user_missing_value)
      endif
   dmiss = MISSING_R8
      else
   user_missing_exists = NF90_NOERR + 20 ! something so we know it does not exist
   endif

! We are going to replace FillValue, miss_value, user_missing_value with one value
! since they are all interpreted the same after being read.

if (user_missing_exists == NF90_NOERR) where(darray == user_missing_value) darray = MISSING_R8
if (     missing_exists == NF90_NOERR) where(darray ==      missing_value) darray = MISSING_R8
if (        fill_exists == NF90_NOERR) where(darray ==          FillValue) darray = MISSING_R8

if (scaling_exists == NF90_NOERR) then
   where (darray /= MISSING_R8) darray = darray * scale_factor
endif

if (offset_exists == NF90_NOERR) then
   where (darray /= MISSING_R8) darray = darray + add_offset
endif
  
end subroutine getvar_real_2d


!--------------------------------------------------------------------
!>   getvar_int_2d - subroutine that inquires, gets the variable, and fills 
!>            in the missing value attribute if that arg is present.
!>            gets the entire array, no start or count specified.
!>            this version assumes you are reading an entire 2d array.
!>            see the slice versions for reading 1d from a 2d array.
!>
!>      ncid - open netcdf file handle
!>      varname - string name of netcdf variable
!>      darray - 2d output array.  integer
!>      dmiss - value that signals a missing value   integer, optional
!>
!>     created 11 Mar 2010,  nancy collins,  ncar/image

subroutine getvar_int_2d(ncid, varname, darray, dmiss)

 integer,            intent(in)   :: ncid
 character(len = *), intent(in)   :: varname
 integer,            intent(out)  :: darray(:,:)
 integer,  optional, intent(out)  :: dmiss

integer  :: varid, nfrc
integer  :: fill, miss
logical  :: found_miss

! read the data for the requested array, and get the fill value
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'getvar_int_2d', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, darray), &
               'getvar_int_2d', 'getting var '// trim(varname))

! see long comment in getvar_real() about the logic here.
if (present(dmiss)) then
   dmiss = MISSING_I
   found_miss = .false.

   ! if user defined another attribute name for missing vals
   ! look for it first, and make it an error if it's not there?
   if (missing_name /= '') then
      nfrc = nf90_get_att(ncid, varid, missing_name, miss)
      if (nfrc == NF90_NOERR) then 
         found_miss = .true.
         dmiss = miss
      endif
   endif

      ! this would make it a fatal error if not found
      !call nc_check( nf90_get_att(ncid, varid, missing_name', miss), &
      !         'getvar_int_2d', 'getting attr "//trim(missing_name)//" for '//trim(varname))

   ! the predefined netcdf fill value.
   nfrc = nf90_get_att(ncid, varid, '_FillValue', fill)
   if (nfrc == NF90_NOERR) then
      if (.not. found_miss) then  
         found_miss = .true.
         dmiss = fill
      else
         ! found both, set all to fill value
         where(darray .eq. miss) darray = fill 
         dmiss = fill
      endif
   endif

   ! if you wanted an error if you specified dmiss and neither attr are
   ! there, you'd test found_miss here.  if it's still false, none were there.

endif
  
end subroutine getvar_int_2d


!--------------------------------------------------------------------
!>  getvar_real_1d_1val - subroutine that inquires, gets the variable, and fills 
!>            in the missing value attribute if that arg is present.
!>            takes a single start, uses count=1, returns a scalar
!>
!>      ncid - open netcdf file handle
!>      varname - string name of netcdf variable
!>      varstart - starting index in the 1d array
!>      dout - output value.  real(r8)
!>      dmiss - value that signals a missing value   real(r8), optional
!>
!>     created 11 Mar 2010,  nancy collins,  ncar/image

subroutine getvar_real_1d_1val(ncid, varname, varstart, dout, dmiss)

integer,            intent(in)   :: ncid
character(len = *), intent(in)   :: varname
integer,            intent(in)   :: varstart
real(r8),           intent(out)  :: dout
real(r8), optional, intent(out)  :: dmiss

integer :: varid, ret

! read the data for the requested array, and get the fill value
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'getvar_real_1d_val', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, dout, start = (/ varstart /) ), &
               'getvar_real_1d_val', 'getting var '// trim(varname))

if (present(dmiss)) then 
   ret = nf90_get_att(ncid, varid, '_FillValue', dmiss)
   if (ret /= NF90_NOERR) dmiss = MISSING_R8
endif

end subroutine getvar_real_1d_1val


!--------------------------------------------------------------------
!>   getvar_int_1d_1val - subroutine that inquires, gets the variable, and fills 
!>            in the missing value attribute if that arg is present.
!>            takes a single start, uses count=1, returns a scalar
!>
!>      ncid - open netcdf file handle
!>      varname - string name of netcdf variable
!>      varstart - starting index in the 1d array
!>      dout - output value.  int
!>      dmiss - value that signals a missing value   int, optional
!>
!>     created 11 Mar 2010,  nancy collins,  ncar/image

subroutine getvar_int_1d_1val(ncid, varname, varstart, dout, dmiss)

integer,            intent(in)   :: ncid
character(len = *), intent(in)   :: varname
integer,            intent(in)   :: varstart
integer,            intent(out)  :: dout
integer,  optional, intent(out)  :: dmiss

integer :: varid, ret

! read the data for the requested array, and get the fill value
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'getvar_int_1d_1val', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, dout, start = (/ varstart /) ), &
               'getvar_int_1d_1val', 'getting var '// trim(varname))

if (present(dmiss)) then
   ret = nf90_get_att(ncid, varid, '_FillValue', dmiss)
   if (ret /= NF90_NOERR) dmiss = MISSING_I
endif

end subroutine getvar_int_1d_1val

!--------------------------------------------------------------------
!>   getvar_real_2d_slice - subroutine that inquires, gets the variable, 
!>           and fills in the missing value attribute if that arg is present.
!>           assumes start = (/ 1, n /) and count = (/ m, 1 /)
!>           so takes a scalar start, count, returns a 1d_array
!>
!>      ncid - open netcdf file handle
!>      varname - string name of 2d netcdf variable
!>      varstart - starting index in the 2d array.  integer
!>      count - nitems to get. integer
!>      darray - 1d output array.  real(r8)
!>      dmiss - value that signals a missing value   real(r8), optional
!>
!>     created 11 Mar 2010,  nancy collins,  ncar/image
!>     updated 14 Jul 2011,  nancy collins,  ncar/image
!>         routine renamed and moved to the utilities module

subroutine getvar_real_2d_slice(ncid, varname, varstart, varcount, darray, dmiss)

integer,            intent(in)   :: ncid
character(len = *), intent(in)   :: varname
integer,            intent(in)   :: varstart
integer,            intent(in)   :: varcount
real(r8),           intent(out)  :: darray(:)
real(r8), optional, intent(out)  :: dmiss

integer :: varid, ret

! read the data for the requested array, and get the fill value
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'getvar_real_2d_slice', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, darray, &
                start=(/ 1, varstart /), count=(/ varcount, 1 /) ), &
               'getvar_real_2d_slice', 'getting var '// trim(varname))

if (present(dmiss)) then
   ret = nf90_get_att(ncid, varid, '_FillValue', dmiss)
   if (ret /= NF90_NOERR) dmiss = MISSING_R8
endif

end subroutine getvar_real_2d_slice


!--------------------------------------------------------------------
!>   get_or_fill_QC_2d_slice - subroutine which gets the requested netcdf variable
!>           but if it isn't there, it fills the array with 0s.  not an
!>           error if it's not present.  assumes integer data array
!>           assumes start = (/ 1, n /) and count = (/ m, 1 /)
!>           so takes a scalar start, count, returns a 1d_array
!>           also prints out a message if fill used.
!>
!>      ncid - open netcdf file handle
!>      varname - string name of 2d netcdf variable
!>      varstart - starting index in the 2d array.  integer
!>      count - nitems to get. integer
!>      darray - output array.  integer
!>
!>     created Mar 8, 2010    nancy collins, ncar/image

subroutine get_or_fill_QC_2d_slice(ncid, varname, varstart, varcount, darray)

integer,            intent(in)    :: ncid
character(len = *), intent(in)    :: varname
integer,            intent(in)    :: varstart
integer,            intent(in)    :: varcount
integer,            intent(inout) :: darray(:)

integer :: varid, nfrc

! test to see if variable is present.  if yes, read it in.
! otherwise, set to fill value, or 0 if none given.

nfrc = nf90_inq_varid(ncid, varname, varid)
if (nfrc == NF90_NOERR) then
   call nc_check( nf90_get_var(ncid, varid, darray, &
                  start=(/ 1, varstart /), count=(/ varcount, 1 /) ), &
                  'get_or_fill_int_2d_slice', 'reading '//trim(varname) )
else
   darray = 0
   if (varstart == 1) & 
     print *, 'QC field named ' // trim(varname) // ' was not found in input, 0 used instead'
endif

end subroutine get_or_fill_QC_2d_slice

!--------------------------------------------------------------------
!> given a list of variable names, check the netcdf file
!> for their existence and return the first one found.  if none 
!> of the given names are in the file, return -1 for index.
!> an optional arg can be used to force it to fail if a match is not found. 
!>
!> ncid - open netcdf file handle
!> nname - number of names in the namelist array
!> namelist - string array of netcdf variable names to test
!> index - index of first name which matched an array.  -1 if none.
!> force - if true, one of the names must match or it is a fatal error.
!>
!> created Mar 15, 2012    nancy collins, ncar/image

subroutine query_varname(ncid, nnames, namelist, index, force)

 integer,             intent(in)    :: ncid, nnames
 character(len = *),  intent(in)    :: namelist(:)
 integer,             intent(out)   :: index
 logical, optional,   intent(in)    :: force

integer :: varid, nfrc, i
character(512) :: msgstring

index = -1
do i=1, nnames
   nfrc = nf90_inq_varid(ncid, namelist(i), varid) 
   if (nfrc == NF90_NOERR) then
      index = i
      return
   endif
enddo
   
if (present(force)) then
   if (index == -1 .and. force) then
      msgstring = 'trying to find one of the following arrays in the input netcdf file'
      call error_handler(E_MSG, 'query_varname', msgstring)
      do i=1, nnames
         call error_handler(E_MSG, 'query_varname', namelist(i))
      enddo
      call error_handler(E_ERR, 'query_varname', 'fatal error, none are present', &
                         source, revision, revdate)
   
   endif
endif

end subroutine query_varname

!--------------------------------------------------------------------

function is_variable_integer(ncid, varname)
integer,            intent(in) :: ncid
character(len = *), intent(in) :: varname
logical :: is_variable_integer

integer :: varid, nfrc, xtype

nfrc = nf90_inq_varid(ncid, varname, varid)
if (nfrc == NF90_NOERR) then
   nfrc = nf90_inquire_variable(ncid, varid, xtype = xtype)
   is_variable_integer = (xtype == NF90_INT)
else
   is_variable_integer = .false.
endif

end function is_variable_integer

!--------------------------------------------------------------------

! this needs to accept either float or double because depending on
! how this code is compiled reals may be r4 or r8 at runtime.

function is_variable_real(ncid, varname)
integer,            intent(in) :: ncid
character(len = *), intent(in) :: varname
logical :: is_variable_real

integer :: varid, nfrc, xtype

nfrc = nf90_inq_varid(ncid, varname, varid)
if (nfrc == NF90_NOERR) then
   nfrc = nf90_inquire_variable(ncid, varid, xtype = xtype)
   is_variable_real = (xtype == NF90_FLOAT .or. xtype == NF90_DOUBLE)
else
   is_variable_real = .false.
endif

end function is_variable_real

!--------------------------------------------------------------------

end module obs_utilities_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
