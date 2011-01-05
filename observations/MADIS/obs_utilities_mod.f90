! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module obs_utilities_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$


use        types_mod, only : r8, MISSING_R8, MISSING_I
use    utilities_mod, only : nc_check
use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_kind, &
                             set_obs_def_error_variance, set_obs_def_location
use obs_sequence_mod, only : obs_sequence_type, obs_type, insert_obs_in_seq, &
                             set_obs_values, set_qc, set_obs_def
use time_manager_mod, only : time_type, operator(>=), set_time
use     location_mod, only : set_location


use netcdf

implicit none
private

public :: create_3d_obs,    &
          add_obs_to_seq,   &
          getdimlen,        &
          getvar_real,      &
          getvar_int,       &
          get_or_fill_real, &
          get_or_fill_int,  &
          get_or_fill_QC,   &
          set_missing_name


! module global storage
character(len=NF90_MAX_NAME) :: missing_name = ''

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   create_3d_obs - subroutine that is used to create an observation
!                   type from observation data.  
!
!       NOTE: assumes the code is using the threed_sphere locations module, 
!             that the observation has a single data value and a single
!             qc value, and that this obs type has no additional required
!             data (e.g. gps and radar obs need additional data per obs)
!
!    lat   - latitude of observation
!    lon   - longitude of observation
!    vval  - vertical coordinate
!    vkind - kind of vertical coordinate (pressure, level, etc)
!    obsv  - observation value
!    okind - observation kind
!    oerr  - observation error
!    day   - gregorian day
!    sec   - gregorian second
!    qc    - quality control value
!    obs   - observation type
!
!     created Oct. 2007 Ryan Torn, NCAR/MMM
!     adapted for more generic use 11 Mar 2010, nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine create_3d_obs(lat, lon, vval, vkind, obsv, okind, oerr, day, sec, qc, obs)
 integer, intent(in)         :: okind, vkind, day, sec
 real(r8), intent(in)        :: lat, lon, vval, obsv, oerr, qc
 type(obs_type), intent(inout) :: obs

real(r8)              :: obs_val(1), qc_val(1)
type(obs_def_type)    :: obs_def

call set_obs_def_location(obs_def, set_location(lon, lat, vval, vkind))
call set_obs_def_kind(obs_def, okind)
call set_obs_def_time(obs_def, set_time(sec, day))
call set_obs_def_error_variance(obs_def, oerr * oerr)
call set_obs_def(obs, obs_def)

obs_val(1) = obsv
call set_obs_values(obs, obs_val)
qc_val(1)  = qc
call set_qc(obs, qc_val)

end subroutine create_3d_obs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   add_obs_to_seq -- adds an observation to a sequence.  inserts if first
!           obs, inserts with a prev obs to save searching if that's possible.
!
!     seq - observation sequence to add obs to
!     obs - observation, already filled in, ready to add
!     obs_time - time of this observation, in dart time_type format
!     prev_obs - the previous observation that was added to this sequence
!                (will be updated by this routine)
!     prev_time - the time of the previously added observation (will also
!                be updated by this routine)
!     first_obs - should be initialized to be .true., and then will be
!                updated by this routine to be .false. after the first obs
!                has been added to this sequence.
!
!     created Mar 8, 2010   nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   getdimlen - subroutine that inquires, gets the dimension size
!
!      ncid - open netcdf file handle
!      dimname - string name of netcdf dimension
!      dout - output value.  integer
!
!     created 11 Mar 2010,  nancy collins,  ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   set_missing_name - subroutine that inquires, gets the variable, and fills 
!            in the missing value attribute if that arg is present.
!            gets the entire array, no start or count specified.
!
!      ncid - open netcdf file handle
!      varname - string name of netcdf variable
!      darray - output array.  real(r8)
!      dmiss - value that signals a missing value   real(r8), optional
!
!     created 11 Mar 2010,  nancy collins,  ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_missing_name(name)
 character(len = *), intent(in)   :: name

   if (len(name) > NF90_MAX_NAME) then
      print *, 'set_missing_name: name must be less than ', NF90_MAX_NAME, ' chars long'
      print *, 'set_missing_name: name is ', len(name), ' long'
      stop
   endif

   missing_name = name

end subroutine set_missing_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   getvar_real - subroutine that inquires, gets the variable, and fills 
!            in the missing value attribute if that arg is present.
!            gets the entire array, no start or count specified.
!
!      ncid - open netcdf file handle
!      varname - string name of netcdf variable
!      darray - output array.  real(r8)
!      dmiss - value that signals a missing value   real(r8), optional
!
!     created 11 Mar 2010,  nancy collins,  ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   getvar_int - subroutine that inquires, gets the variable, and fills 
!            in the missing value attribute if that arg is present.
!            gets the entire array, no start or count specified.
!
!      ncid - open netcdf file handle
!      varname - string name of netcdf variable
!      darray - output array.  integer
!      dmiss - value that signals a missing value   integer, optional
!
!     created 11 Mar 2010,  nancy collins,  ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getvar_int(ncid, varname, darray, dmiss)
 integer,            intent(in)   :: ncid
 character(len = *), intent(in)   :: varname
 integer,            intent(out)  :: darray(:)
 integer,  optional, intent(out)  :: dmiss

integer  :: varid, nfrc
real(r8) :: fill, miss
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_or_fill_real - subroutine which gets the requested netcdf variable
!           but if it isn't there, it fills the array with 0s.  not an
!           error if it's not present.  assumes real(r8) data array
!
!      ncid - open netcdf file handle
!      varname - string name of netcdf variable
!      darray - output array.  real(r8)
!      fillval - real(r8) to fill with if not present.  optional.
!      did_fill - reports back if fill was done.  logical, optional.
!
!     created Mar 8, 2010    nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_or_fill_int - subroutine which gets the requested netcdf variable
!           but if it isn't there, it fills the array with 0s.  not an
!           error if it's not present.  assumes integer data array
!
!      ncid - open netcdf file handle
!      varname - string name of netcdf variable
!      darray - output array.  integer
!      fillval - integer to fill with if not present.  optional.
!      did_fill - reports back if fill was done.  logical, optional.
!
!     created Mar 8, 2010    nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! get_or_fill_QC - get the expected QC netcdf variable.  if not found,
!                  fill the array with 0s, and print out a message warning
!                  that the input didn't have the expected array.
!   
!      minor tweak on the standard get_or_fill routine.  prints out a
!      warning message, and don't return to the caller which one it did.
!      would be more general if it didn't use 'QC' in the output string.
!      if anyone else wants to use it, we could change 'QC' to 'warn' in
!      the name, and remove the 'QC' from the output string, and make it
!      slightly more general
!
!      ncid - open netcdf file handle
!      varname - string name of netcdf variable
!      darray - output array.  integer
!
!      created 17 mar 2010  nsc  ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_or_fill_QC(ncid, varname, darray)
 integer,            intent(in)    :: ncid
 character(len = *), intent(in)    :: varname
 integer,            intent(inout) :: darray(:)

logical :: did_fill


call get_or_fill_int(ncid, varname, darray, 0, did_fill)

if (did_fill) &
  print *, 'QC field named ' // trim(varname) // ' was not found in input, 0 used instead'

end subroutine get_or_fill_QC



end module obs_utilities_mod
