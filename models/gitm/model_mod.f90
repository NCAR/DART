! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

! This is the interface between the gitm model and DART.

! Modules that are absolutely required for use are listed

use        types_mod, only : r4, r8, digits12, SECPERDAY, MISSING_R8,          &
                             rad2deg, deg2rad, PI, obstypelength, i8

use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time,&
                             print_time, print_date, set_calendar_type,        &
                             operator(*),  operator(+), operator(-),           &
                             operator(>),  operator(<), operator(/),           &
                             operator(/=), operator(<=)

use     location_mod, only : location_type, get_dist, query_location,          &
                             get_close_type, VERTISHEIGHT,                     &
                             set_location, get_location,                       &
                             loc_get_close_obs => get_close_obs, is_vertical,  &
                             vertical_localization_on

use    utilities_mod, only : register_module, error_handler,                   &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                             nc_check, do_output, to_upper,                    &
                             find_namelist_in_file, check_namelist_read,       &
                             open_file, file_exist, find_textfile_dims,        &
                             file_to_text, close_file

use netcdf_utilities_mod   ! later make an 'only' list once we know what we need

use  ensemble_manager_mod,  only : ensemble_type

use distributed_state_mod,  only : get_state

use     obs_kind_mod, only : get_index_for_quantity,  &
                             get_name_for_quantity,   &
                             QTY_GEOPOTENTIAL_HEIGHT

use     default_model_mod,  only : adv_1step, nc_write_model_vars, &
                                   pert_model_copies, get_close_state, &
                                   read_model_time, write_model_time, &
                                   convert_vertical_obs, convert_vertical_state, &
                                   init_time => fail_init_time, &
                                   init_conditions => fail_init_conditions

use     dart_gitm_mod, only: get_gitm_nLons, get_gitm_nLats, get_gitm_nAlts, &
                             get_nSpecies, get_nSpeciesTotal, get_nIons,     &
                             get_nSpeciesAll, decode_gitm_indices

use typesizes
use netcdf

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.
! this list has code in this module.
public :: get_model_size,         &
          get_state_meta_data,    &
          model_interpolate,      &
          shortest_time_between_assimilations, &
          static_init_model,      &
          end_model,              &
          get_close_obs,          &
          get_close_state,        &
          convert_vertical_obs,   &
          convert_vertical_state, &
          read_model_time,        &
          write_model_time 

! these routines also must be public.
! this list are names of routines where the code
! is passed through from other modules
public :: init_time,              &
          init_conditions,        &
          adv_1step,              &
          nc_write_model_atts,    &
          nc_write_model_vars,    &
          pert_model_copies


! generally useful routines for various support purposes.
! the interfaces here can be changed as appropriate.

public :: get_gitm_restart_dirname,    &
          get_state_time
          ! only used in a converter but this code
          ! has not been updated.  must go to/from
          ! gitm binary and netcdf - not 1d array.
          !restart_file_to_statevector, &
          !statevector_to_restart_file, &

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=256) :: string1, string2
logical, save :: module_initialized = .false.

!------------------------------------------------------------------
! things which can/should be in the model_nml
!
!  The DART state vector may consist of things like:
!
!  U    long_name = "X-WIND COMPONENT"      float   U(TIME, ALT, LAT, XE)
!  V    long_name = "Y-WIND COMPONENT"      float   V(TIME, ALT, YE, LON)
!  W    long_name = "Z-WIND COMPONENT"      float   W(TIME, ZE, LAT, LON)
!  TH   long_name = "POTENTIAL TEMPERATURE" float  TH(TIME, ALT, LAT, LON)
!  DBZ  long_name = "RADAR REFLECTIVITY"    float DBZ(TIME, ALT, LAT, LON)
!  WZ   long_name = "VERTICAL VORTICITY"    float  WZ(TIME, ALT, LAT, LON)
!  PI   long_name = "PERT. EXNER"	    float  PI(TIME, ALT, LAT, LON)
!  QV   long_name = "VAPOR MIXING RATIO"    float  QV(TIME, ALT, LAT, LON)
!  QC   long_name = "CLOUD MIXING RATIO"    float  QC(TIME, ALT, LAT, LON)
!  QR   long_name = "RAIN MIXING RATIO"     float  QR(TIME, ALT, LAT, LON)
!  QI   long_name = "ICE MIXING RATIO"      float  QI(TIME, ALT, LAT, LON)
!  QS   long_name = "SNOW MIXING RATIO"     float  QS(TIME, ALT, LAT, LON)
!  QH   long_name = "GRAUPEL MIXING RATIO"  float  QH(TIME, ALT, LAT, LON)
!
!  The variables in the gitm restart file that are used to create the
!  DART state vector are specified in input.nml:model_nml: gitm_state_variables
!
!------------------------------------------------------------------

integer, parameter :: max_state_variables = 80
integer, parameter :: num_state_table_columns = 2
character(len=NF90_MAX_NAME) :: gitm_state_variables(max_state_variables * num_state_table_columns ) = ' '
character(len=NF90_MAX_NAME) :: variable_table(max_state_variables, num_state_table_columns )

integer            :: assimilation_period_days = 0
integer            :: assimilation_period_seconds = 60
real(r8)           :: model_perturbation_amplitude = 0.2
logical            :: output_state_vector = .false.
integer            :: debug = 0   ! turn up for more and more debug messages
character(len=32)  :: calendar = 'Gregorian'
character(len=256) :: gitm_restart_dirname = 'gitm_restartdir'

namelist /model_nml/  &
   gitm_restart_dirname,        &
   output_state_vector,         &
   assimilation_period_days,    &  ! for now, this is the timestep
   assimilation_period_seconds, &
   model_perturbation_amplitude,&
   calendar,                    &
   debug,                       &
   gitm_state_variables

integer :: nfields

! Everything needed to describe a variable

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname       ! crazy species name
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=NF90_MAX_NAME) :: storder
   character(len=NF90_MAX_NAME) :: gitm_varname  ! NDensityS, IDensityS, ...
   integer :: gitm_dim                           ! dimension defining species
   integer :: gitm_index                         ! 'iSpecies' or u,v,w ...
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens ! nlons, nlats, nalts [, nspecies]
   integer :: posdef
   integer :: numdims
   integer :: varsize     ! prod(dimlens(1:numdims))
   integer :: index1      ! location in dart state vector of first occurrence
   integer :: indexN      ! location in dart state vector of last  occurrence
   integer :: dart_kind
   character(len=obstypelength) :: kind_string
end type progvartype

type(progvartype), dimension(max_state_variables) :: progvar

! nLons, nLats are the number of lons/lats PER block
!                  the number of blocks comes from UAM.in
! nAlts  is the one and only number of altitudes ... no block-dependence

integer :: nLons, nLats, nAlts

! "... keep in mind that if the model resolution is 5 deg latitude,
!  the model will actually go from -87.5 to 87.5 latitude
! (even though you specify -90 to 90 in the UAM.in file),
! since the latitudes/longitudes are at cell centers,
! while the edges are at the boundaries." -- Aaron Ridley

integer  :: NgridLon=-1, NgridLat=-1, NgridAlt=-1    ! scalar grid counts
integer  :: nBlocksLon=-1, nBlocksLat=-1             ! number of blocks along each dim
real(r8) :: LatStart=MISSING_R8, LatEnd=MISSING_R8, LonStart=MISSING_R8
integer  :: nSpeciesTotal=-1, nSpecies=-1, nIons=-1, nSpeciesAll=-1

! scalar grid positions

real(r8), allocatable :: LON(:)   ! longitude centers
real(r8), allocatable :: LAT(:)   ! latitude  centers
real(r8), allocatable :: ALT(:)   ! vertical level centers

integer               :: model_size      ! the state vector length
type(time_type)       :: model_time      ! valid time of the model state
type(time_type)       :: model_advance_time  ! smallest time to adv model
real(r8), allocatable :: ens_mean(:)     ! may be needed for forward ops

! set this to true if you want to print out the current time
! after each N observations are processed, for benchmarking.

logical :: print_timestamps = .false.
integer :: print_every_Nth  = 10000

integer, parameter :: nGhost = 2   ! number of ghost cells on all edges

!------------------------------------------------------------------
! The gitm restart manager namelist variables
!------------------------------------------------------------------

character(len= 64) :: ew_boundary_type, ns_boundary_type

INTERFACE get_index_range
      MODULE PROCEDURE get_index_range_string
      MODULE PROCEDURE get_index_range_int
END INTERFACE

contains


!==================================================================
! All the REQUIRED interfaces come first - just by convention.
!==================================================================


function get_model_size()
!------------------------------------------------------------------

! Returns the size of the model as an integer.
! Required for all applications.

integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size


!==================================================================


subroutine get_state_meta_data(index_in, location, var_type)
!------------------------------------------------------------------
! given an index into the state vector, return its location and
! if given, the var kind.   despite the name, var_type is a generic
! kind, like those in obs_kind/obs_kind_mod.f90, starting with QTY_

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_type

! Local variables

integer :: lat_index, lon_index, alt_index
integer :: n, nf, myindx, remainder, remainder2

if ( .not. module_initialized ) call static_init_model

! Find out which of the 3D fields index_in is part of
nf     = -1

FindIndex : do n = 1,nfields
   if( (progvar(n)%index1 <= index_in) .and. (index_in <= progvar(n)%indexN) ) then
      nf = n
      myindx = index_in - progvar(n)%index1 + 1
      exit FindIndex
   endif
enddo FindIndex

if( myindx == -1 ) then
   write(string1,*) 'Problem, cannot find base_offset, index_in is: ', index_in
   call error_handler(E_ERR,'get_state_meta_data',string1,source,revision,revdate)
endif

alt_index = 1 + (myindx - 1) / (NgridLon * NgridLat)
remainder = myindx - (alt_index-1) * NgridLon * NgridLat
lat_index = 1 + (remainder - 1) / NgridLon
remainder2 = remainder - (lat_index - 1) * NgridLon
lon_index = remainder2

location = set_location(LON(lon_index), LAT(lat_index), ALT(alt_index), VERTISHEIGHT)

if (present(var_type)) then
   var_type = progvar(nf)%dart_kind
endif

end subroutine get_state_meta_data


!==================================================================


subroutine model_interpolate(state_handle, ens_size, location, obs_type, interp_val, istatus)
!------------------------------------------------------------------
!     PURPOSE:
!
!     For a given lat, lon, and height, interpolate the correct state value
!     to that location for the filter from the gitm state vectors
!
!     Variables needed to be stored in the MODEL_MODULE data structure
!
!       LON   = 1D array storing the local grid center coords (degrees)
!       LAT   = 1D array storing the local grid center coords (degrees)
!       ALT   = 1D array storing the local grid center coords (meters)
!
!       ERROR codes:
!
!       ISTATUS = 99:  general error in case something terrible goes wrong...
!       ISTATUS = 15:  dont know what to do with vertical coord supplied
!       ISTATUS = 16:  longitude illegal
!       ISTATUS = 17:  latitude illegal
!       ISTATUS = 18:  altitude illegal
!       ISTATUS = 20:  asking to interpolate an unknown obs kind

! Passed variables

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_type
real(r8),            intent(out) :: interp_val(ens_size)
integer,             intent(out) :: istatus(ens_size)

! Local storage

real(r8) :: loc_array(3), llon, llat, lvert, lon_fract, lat_fract, alt_fract
integer  :: blon(2), blat(2), balt(2), i, j, k, ier, nhgt
integer(i8) :: base_offset, end_offset
real(r8) :: cube(2, 2, 2, ens_size), square(2, 2, ens_size), line(2, ens_size)

if ( .not. module_initialized ) call static_init_model

! Assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the interp_val will be set to the
! good value, and the last line here sets istatus to 0.
! make any error codes set here be in the 10s

interp_val = MISSING_R8     ! the DART bad value flag
istatus = 99                ! unknown error

! Get the individual locations values

loc_array = get_location(location)
llon      = loc_array(1)
llat      = loc_array(2)
lvert     = loc_array(3)

! IF (debug > 2) print *, 'requesting interpolation at ', llon, llat, lvert
print *, 'requesting interpolation at ', llon, llat, lvert

! Only height and level for vertical location type is supported at this point
if (.not. is_vertical(location, "HEIGHT") .and. .not. is_vertical(location, "LEVEL")) THEN
     istatus = 15
     return
endif

! Find the start and end offsets for this field in the state vector x(:)

! FIXME: this fails if you ask for a kind that doesn't exist in the state vector.
! in most cases you just want to set a bad istatus and return without stopping
! the entire assimilation.  if the error calls are commented out of the two
! index range subroutines, the offsets will come back as 0.  since the
! return values have already been set, just give it a more specific error
! code and return here.

if (obs_type == QTY_GEOPOTENTIAL_HEIGHT ) then
   ! ok to continue.  offsets unused in this case, but
   ! set them to something > 0 to indicate its ok.
   base_offset = 1
   end_offset = 1
else 
   call get_index_range(obs_type, base_offset, end_offset)
endif

if (debug > 2) print *, 'base offset now ', base_offset

! fail if this kind isn't in the state vector.
if (base_offset <= 0) then 
   istatus = 20
   return
endif

! Need to find bounding center indices and fractional offsets for lat, lon, alt
call find_lon_bounds(llon, blon(1), blon(2), lon_fract, ier)
if(ier /= 0) then
   istatus = 16
   return
endif

call find_lat_or_alt_bounds(llat, NgridLat, LAT, blat(1), blat(2), lat_fract, ier)
if(ier /= 0) then
   istatus = 17
   return
endif

if (is_vertical(location, "HEIGHT")) then ! call the height interpolation routine.

   call find_lat_or_alt_bounds(lvert, NgridAlt, ALT, balt(1), balt(2), alt_fract, ier)

else if (is_vertical(location, "LEVEL")) then ! set the levels and fraction.

   nhgt = int(lvert)
   if (nhgt < 1 .or. nhgt > NgridAlt) then
      istatus = 18
      return
   endif

   ! if we are below the top level, set the lower bound to the integer part of the
   ! level and set the fraction between it and the next level to the fractional part.
   ! if we are asking for the top level, set the upper bound to the requested level 
   ! and set the fraction to 1.
   if (nhgt < NGridAlt) then
      balt(1) = nhgt
      balt(2) = nhgt + 1
      alt_fract = lvert - real(nhgt,r8)
   else
      balt(1) = nhgt - 1
      balt(2) = nhgt
      alt_fract = 1.0_r8
   endif
   ier = 0
else
   ! shouldn't happen
   istatus = 99
   return
endif
if(ier /= 0) then
   istatus = 18
   return
endif

! if we're asking about height, we have the alt arrays directly.
if (obs_type == QTY_GEOPOTENTIAL_HEIGHT) then

   ! Interpolate to the given altitude - lat/lon doesn't matter here.
   interp_val = (1 - alt_fract) * ALT(balt(1)) + alt_fract * ALT(balt(2))
   
   istatus = 0
   return
endif

! for the rest of the state vector contents - find the offset to the start
! of the requested kind, and do a tri-linear interpolation inside the enclosing
! cube from the grid.

! Get the grid values for the first
do i = 1, 2
   do j = 1, 2
      do k = 1, 2
         cube(i, j, k, :) = get_grid_value(base_offset, blon(i), blat(j), balt(k), state_handle, ens_size)
      end do
   end do
end do

! Interpolate to the given altitude
do i = 1, 2
   do j = 1, 2
      square(i, j, :) = (1 - alt_fract) * cube(i, j, 1, :) + alt_fract * cube(i, j, 2, :)
   end do
end do

! Interpolate to the given latitude
do i = 1, 2
   line(i, :) = (1 - lat_fract) * square(i, 1, :) + lat_fract * square(i, 2, :)
end do

! Interpolate to the given longitude
interp_val(:) = (1 - lon_fract) * line(1, :) + lon_fract * line(2, :)

! All good.
istatus(:) = 0

end subroutine model_interpolate


!==================================================================


function shortest_time_between_assimilations()
!------------------------------------------------------------------
! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.

type(time_type) :: shortest_time_between_assimilations

if ( .not. module_initialized ) call static_init_model

shortest_time_between_assimilations = model_advance_time

end function shortest_time_between_assimilations


!==================================================================


subroutine static_init_model()
!------------------------------------------------------------------
! Called to do one time initialization of the model.
!
! All the grid information comes from the initialization of
! the dart_gitm_mod module.

! Local variables - all the important ones have module scope

character(len=NF90_MAX_NAME)    :: varname
character(len=obstypelength) :: kind_string
integer :: varsize
integer :: iunit, io, ivar, index1, indexN
integer :: ss, dd

if ( module_initialized ) return ! only need to do this once.

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.
module_initialized = .true.

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
call error_handler(E_MSG,'static_init_model','model_nml values are',' ',' ',' ')
if (do_output()) write(logfileunit, nml=model_nml)
if (do_output()) write(     *     , nml=model_nml)

! Get the GITM variables in a restricted scope setting.

nLons         = get_gitm_nLons()
nLats         = get_gitm_nLats()
nAlts         = get_gitm_nAlts()
nSpecies      = get_nSpecies()
nSpeciesTotal = get_nSpeciesTotal()
nIons         = get_nIons()
nSpeciesAll   = get_nSpeciesAll()

if ((debug > 0) .and.  do_output() ) then
   write(*,*)
   write(*,*)'nLons         is ',nLons
   write(*,*)'nLats         is ',nLats
   write(*,*)'nAlts         is ',nAlts
   write(*,*)'nSpecies      is ',nSpecies
   write(*,*)'nSpeciesTotal is ',nSpeciesTotal
   write(*,*)'nIons         is ',nIons
   write(*,*)'nSpeciesAll   is ',nSpeciesAll
endif

!---------------------------------------------------------------
! Set the time step ... causes gitm namelists to be read.
! Ensures model_advance_time is multiple of 'dynamics_timestep'

call set_calendar_type( calendar )   ! comes from model_mod_nml

model_advance_time = set_model_time_step()

call get_time(model_advance_time,ss,dd) ! set_time() assures the seconds [0,86400)

write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model',string1,source,revision,revdate)

!---------------------------------------------------------------
! 1) get grid dimensions
! 2) allocate space for the grids
! 3) read them from the block restart files, could be stretched ...

call get_grid_info(NgridLon, NgridLat, NgridAlt, nBlocksLon, nBlocksLat, &
                   LatStart, LatEnd, LonStart)

if( debug  > 0 ) then
   write(*,*)'grid dims are ',NgridLon,NgridLat,NgridAlt
endif

allocate( LON( NgridLon ))
allocate( LAT( NgridLat ))
allocate( ALT( NgridAlt ))

call get_grid(gitm_restart_dirname, nBlocksLon, nBlocksLat, &
              nLons, nLats, nAlts, LON, LAT, ALT )

!---------------------------------------------------------------
! Compile the list of gitm variables to use in the creation
! of the DART state vector. Required to determine model_size.
!
! Verify all variables are in the gitm restart file
!
! Compute the offsets into the state vector for the start of each
! different variable type. Requires reading shapes from the gitm
! restart file. As long as TIME is the LAST dimension, we're OK.
!
! Record the extent of the data type in the state vector.

call verify_state_variables( gitm_state_variables, nfields, variable_table)

index1  = 1;
indexN  = 0;

do ivar = 1, nfields

   varname                   = trim(variable_table(ivar,1))
   kind_string               = trim(variable_table(ivar,2))
   progvar(ivar)%varname     = varname
   progvar(ivar)%kind_string = kind_string
   progvar(ivar)%dart_kind   = get_index_for_quantity( progvar(ivar)%kind_string )
   progvar(ivar)%dimlens     = 0

   ! I would really like decode_gitm_indices to set the following (on a per-variable basis)
   ! progvar(ivar)%storder
   ! progvar(ivar)%numdims
   ! progvar(ivar)%dimlens

   ! This routine also checks to make sure user specified accurate GITM variables
   call decode_gitm_indices( varname, &
                             progvar(ivar)%gitm_varname, &
                             progvar(ivar)%gitm_dim,     &
                             progvar(ivar)%gitm_index,   &
                             progvar(ivar)%long_name,    &
                             progvar(ivar)%units)

   if (progvar(ivar)%varname == 'f107') then ! if we are dealing with f107
      varsize = 1
      progvar(ivar)%storder     = '0d'
      progvar(ivar)%numdims     = 1
      progvar(ivar)%dimlens     = 1

   else !anything but f107
      varsize = NgridLon * NgridLat * NgridAlt
      progvar(ivar)%storder     = 'xyz3d'
      progvar(ivar)%numdims     = 3
      progvar(ivar)%dimlens(1:progvar(ivar)%numdims) = (/ NgridLon, NgridLat, NgridAlt /)

   endif

   progvar(ivar)%varsize     = varsize
   progvar(ivar)%index1      = index1
   progvar(ivar)%indexN      = index1 + varsize - 1
   index1                    = index1 + varsize      ! sets up for next variable

   if ( debug > 0 ) then
      write(logfileunit,*)
      write(logfileunit,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(logfileunit,*) ' storage      ',trim(progvar(ivar)%storder)
      write(logfileunit,*) ' long_name    ',trim(progvar(ivar)%long_name)
      write(logfileunit,*) ' units        ',trim(progvar(ivar)%units)
      write(logfileunit,*) ' numdims      ',progvar(ivar)%numdims
      write(logfileunit,*) ' dimlens      ',progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
      write(logfileunit,*) ' varsize      ',progvar(ivar)%varsize
      write(logfileunit,*) ' index1       ',progvar(ivar)%index1
      write(logfileunit,*) ' indexN       ',progvar(ivar)%indexN
      write(logfileunit,*) ' dart_kind    ',progvar(ivar)%dart_kind
      write(logfileunit,*) ' kind_string  ',trim(progvar(ivar)%kind_string)
      write(logfileunit,*) ' gitm_varname ',trim(progvar(ivar)%gitm_varname)
      write(logfileunit,*) ' gitm_dim     ',progvar(ivar)%gitm_dim
      write(logfileunit,*) ' gitm_index   ',progvar(ivar)%gitm_index

      write(     *     ,*)
      write(     *     ,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(     *     ,*) ' storage      ',trim(progvar(ivar)%storder)
      write(     *     ,*) ' long_name    ',trim(progvar(ivar)%long_name)
      write(     *     ,*) ' units        ',trim(progvar(ivar)%units)
      write(     *     ,*) ' numdims      ',progvar(ivar)%numdims
      write(     *     ,*) ' dimlens      ',progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
      write(     *     ,*) ' varsize      ',progvar(ivar)%varsize
      write(     *     ,*) ' index1       ',progvar(ivar)%index1
      write(     *     ,*) ' indexN       ',progvar(ivar)%indexN
      write(     *     ,*) ' dart_kind    ',progvar(ivar)%dart_kind
      write(     *     ,*) ' kind_string  ',trim(progvar(ivar)%kind_string)
      write(     *     ,*) ' gitm_varname ',trim(progvar(ivar)%gitm_varname)
      write(     *     ,*) ' gitm_dim     ',progvar(ivar)%gitm_dim
      write(     *     ,*) ' gitm_index   ',progvar(ivar)%gitm_index
   endif

enddo

model_size = progvar(nfields)%indexN

if ( debug > 0 ) then
  write(logfileunit,'("grid: NgridLon, NgridLat, NgridAlt =",3(1x,i5))') NgridLon, NgridLat, NgridAlt
  write(     *     ,'("grid: NgridLon, NgridLat, NgridAlt =",3(1x,i5))') NgridLon, NgridLat, NgridAlt
  write(logfileunit, *)'model_size = ', model_size
  write(     *     , *)'model_size = ', model_size
endif

allocate( ens_mean(model_size) )

end subroutine static_init_model


!==================================================================


subroutine end_model()
!------------------------------------------------------------------
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

if (allocated(LON)) deallocate(LON, LAT, ALT)

end subroutine end_model


!==================================================================


subroutine nc_write_model_atts(ncid, domain_id )
  integer, intent(in) :: ncid 
  integer, intent(in) :: domain_id

integer :: i, myndims
character(len=129), allocatable, dimension(:) :: textblock
integer :: nlines, linelen
logical :: has_gitm_namelist



call nc_begin_define_mode(ncid)
call nc_add_global_attribute(ncid, 'model', 'gitm')

!-------------------------------------------------------------------------------
! Determine shape of most important namelist
!-------------------------------------------------------------------------------

call find_textfile_dims('gitm_vars.nml', nlines, linelen)
if (nlines > 0) then
   has_gitm_namelist = .true.

   allocate(textblock(nlines))
   textblock = ''

   call nc_define_dimension(ncid, 'nlines',  nlines)
   call nc_define_dimension(ncid, 'linelen', linelen)
   call nc_define_character_variable(ncid, 'gitm_in', (/ 'nlines ', 'linelen' /))
   call nc_add_attribute_to_variable(ncid, 'gitm_in', 'long_name', 'contents of gitm_in namelist')

else
  has_gitm_namelist = .false.
endif

!----------------------------------------------------------------------------
! output only grid info - state vars will be written by other non-model_mod code
!----------------------------------------------------------------------------

call nc_define_dimension(ncid, 'LON', NgridLon)
call nc_define_dimension(ncid, 'LAT', NgridLat)
call nc_define_dimension(ncid, 'ALT', NgridAlt)
call nc_define_dimension(ncid, 'WL',  1)  ! wavelengths - currently only 1?

!----------------------------------------------------------------------------
! Create the (empty) Coordinate Variables and the Attributes
!----------------------------------------------------------------------------

! Grid Longitudes
call nc_define_real_variable(ncid, 'LON', (/ 'LON' /) )
call nc_add_attribute_to_variable(ncid, 'LON', 'type',           'x1d')
call nc_add_attribute_to_variable(ncid, 'LON', 'long_name',      'grid longitudes')
call nc_add_attribute_to_variable(ncid, 'LON', 'cartesian_axis', 'X')
call nc_add_attribute_to_variable(ncid, 'LON', 'units',          'degrees_east')
call nc_add_attribute_to_variable(ncid, 'LON', 'valid_range',     (/ 0.0_r8, 360.0_r8 /) )

! Grid Latitudes
call nc_define_real_variable(ncid, 'LAT', (/ 'LAT' /) )
call nc_add_attribute_to_variable(ncid, 'LAT', 'type',           'y1d')
call nc_add_attribute_to_variable(ncid, 'LAT', 'long_name',      'grid latitudes')
call nc_add_attribute_to_variable(ncid, 'LAT', 'cartesian_axis', 'Y')
call nc_add_attribute_to_variable(ncid, 'LAT', 'units',          'degrees_north')
call nc_add_attribute_to_variable(ncid, 'LAT', 'valid_range',     (/ -90.0_r8, 90.0_r8 /) )

! Grid Altitudes
call nc_define_real_variable(ncid, 'ALT', (/ 'ALT' /) )
call nc_add_attribute_to_variable(ncid, 'ALT', 'type',           'z1d')
call nc_add_attribute_to_variable(ncid, 'ALT', 'long_name',      'grid altitudes')
call nc_add_attribute_to_variable(ncid, 'ALT', 'cartesian_axis', 'Z')
call nc_add_attribute_to_variable(ncid, 'ALT', 'units',          'meters')
call nc_add_attribute_to_variable(ncid, 'ALT', 'positive',       'up')

! Grid wavelengths
call nc_define_real_variable(ncid, 'WL', (/ 'WL' /) )
call nc_add_attribute_to_variable(ncid, 'WL', 'type',           'x1d')
call nc_add_attribute_to_variable(ncid, 'WL', 'long_name',      'grid wavelengths')
call nc_add_attribute_to_variable(ncid, 'WL', 'cartesian_axis', 'X')
call nc_add_attribute_to_variable(ncid, 'WL', 'units',          'wavelength_index')
call nc_add_attribute_to_variable(ncid, 'WL', 'valid_range',     (/ 0.9_r8, 38.1_r8 /) )

!----------------------------------------------------------------------------
! Finished with dimension/variable definitions, must end 'define' mode to fill.
!----------------------------------------------------------------------------

call nc_end_define_mode(ncid)

!----------------------------------------------------------------------------
! Fill the coordinate variables
!----------------------------------------------------------------------------

call nc_put_variable(ncid, 'LON', LON)
call nc_put_variable(ncid, 'LAT', LAT)
call nc_put_variable(ncid, 'ALT', ALT)
! what about WL?


if (has_gitm_namelist) then
   call file_to_text('gitm_vars.nml', textblock)
   call nc_put_variable(ncid, 'gitm_in', textblock)
   deallocate(textblock)
endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts


!==================================================================

subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, &
                         obs_loc, obs_kind, obs_type, &
                         num_close, close_ind, dist, ens_handle)
!------------------------------------------------------------------
! Given a DART location (referred to as "base") and a set of candidate
! locations & kinds (obs, obs_kind), returns the subset close to the
! "base", their indices, and their distances to the "base" ...
!
! For vertical distance computations, general philosophy is to convert all
! vertical coordinates to a common coordinate. This coordinate type is defined
! in the namelist with the variable "vert_localization_coord".

type(get_close_type),              intent(in)    :: gc
type(location_type),               intent(inout) :: base_obs_loc
integer,                           intent(in)    :: base_obs_kind
type(location_type), dimension(:), intent(inout) :: obs_loc
integer,             dimension(:), intent(in)    :: obs_kind
integer,             dimension(:), intent(in)    :: obs_type
integer,                           intent(out)   :: num_close
integer,             dimension(:), intent(out)   :: close_ind
real(r8),            dimension(:), intent(out)   :: dist
type(ensemble_type), optional,     intent(in)    :: ens_handle

integer                :: t_ind, istatus1, istatus2, k, i, is_in_close_ind, is_in_obs_kind, f107_ind
integer                :: base_which, local_obs_which
real(r8), dimension(3) :: base_array, local_obs_array
type(location_type)    :: local_obs_loc

! Initialize variables to missing status

num_close       = 0
close_ind       = -99
dist            = 1.0e9_r8   !something big and positive (far away)
istatus1        = 0
istatus2        = 0
is_in_obs_kind  = 0
is_in_close_ind = 0
f107_ind        = -37 !a bad index, hopefully out of bounds of obs_kind


! Convert base_obs vertical coordinate to requested vertical coordinate if necessary

base_array = get_location(base_obs_loc)
base_which = nint(query_location(base_obs_loc))

! fixme ...
if (vertical_localization_on()) then
!  if (base_which /= wrf%dom(1)%vert_coord) then
!     call vert_interpolate(ens_mean, base_obs_loc, base_obs_kind, istatus1)
!  elseif (base_array(3) == MISSING_R8) then
!     istatus1 = 1
!  endif
endif

if (istatus1 == 0) then

   ! Loop over potentially close subset of obs priors or state variables
   ! This way, we are decreasing the number of distance computations that will follow.
   ! This is a horizontal-distance operation and we don't need to have the relevant
   ! vertical coordinate information yet (for obs_loc).

   ! FIXME:::::
   !call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
   !                       num_close, close_ind, dist)

!!!! THE following 20-ish+ lines are implementing the search (if f107's dist to obs is to be calculated)
!!!! Alex 03/07/2012
   do i = 1, size(obs_kind) !have to go over the whole size because these are all the candidates
      if (obs_kind(i) .eq. get_index_for_quantity('QTY_1D_PARAMETER')) then !so right now any QTY_1D_PARAMETER will match.
!+ right now the only parameter is f107, but if you add more parameters, you might want to change their localizations, as
!+ right now they will be either all at the meas. location or all far (depending on est_f107 setting in pbs_file.sh)
         is_in_obs_kind = 1 !true
         f107_ind = i !its index
      endif
   enddo
   if (is_in_obs_kind == 1) then !only check the close_ind if f107 needs to be added
      do k = 1, num_close !go only as far as the data is reasonable (not -99 = data missing)
         if (close_ind(k) .eq. f107_ind) then !if is already in close_ind, take note of it
            is_in_close_ind = 1
         endif
      enddo
   endif
   if ((is_in_obs_kind == 1) .and. (is_in_close_ind == 0)) then !if it needs to be added (is in obs_kind), but is not added yet
      num_close = num_close + 1
      close_ind(num_close) = f107_ind
!      write(*,*) "F107 ADDED, n_c, f107_i ", num_close, f107_ind
   endif


   do k = 1, num_close

      t_ind = close_ind(k)
      local_obs_loc   = obs_loc(t_ind)
      local_obs_which = nint(query_location(local_obs_loc))

      ! Convert vertical coordinate to requested vertical coordinate if necessary.
      ! Should only be necessary for obs priors, as state location information already
      ! contains the correct vertical coordinate
      ! (filter_assim's call to get_state_meta_data).
      if (vertical_localization_on()) then
 !fixme       if (local_obs_which /= wrf%dom(1)%vert_coord) then
 !fixme           call vert_interpolate(ens_mean, local_obs_loc, obs_kind(t_ind), istatus2)
            ! Store the "new" location into the original full local array
            obs_loc(t_ind) = local_obs_loc
 !fixme        endif
      endif

      ! Compute distance - set distance to a very large value if vert coordinate is
      ! missing or vert_interpolate returned error (istatus2=1)
      local_obs_array = get_location(local_obs_loc)


      ! TJH FIXME ... the pbs_file script actually modifies the value of dist in
      ! TJH FIXME ... this file and recompiles. NOT APPROPRIATE.
      if (((vertical_localization_on())        .and. &
           (local_obs_array(3) == MISSING_R8)) .or.  &
           (istatus2 == 1)                   ) then
         dist(k) = 1.0e9_r8
      else
         if (close_ind(k) .eq. f107_ind) then !check if we came across the parameter
            dist(k) = 0 !changed by pbs_file script
         else
            dist(k) = get_dist(base_obs_loc, local_obs_loc, base_obs_kind, obs_kind(t_ind))
         endif
      endif
   enddo
endif

end subroutine get_close_obs


!==================================================================
! The remaining PUBLIC interfaces come next
!==================================================================


subroutine restart_file_to_statevector(dirname, state_vector, model_time)
!------------------------------------------------------------------
! Reads the current time and state variables from a gitm restart
! file and packs them into a dart state vector.
!
! FIXME:
! this routine needs:
! 1.  a base dirname for the restart files.
! they will have the format 'dirname/bNNNN.rst'  where NNNN has
! leading 0s and is the block number.   blocks start in the
! southwest corner of the lat/lon grid and go west first, then
! north and end in the northeast corner.   (assuming var 'dirname')
! the other info is in 'dirname/header.rst'
!
! 2. the overall grid size, lon/lat/alt when you've read in all
!    the blocks.  (nGridLon, nGridLat, nGridAlt, will compute totalVarSize)
!
! 3. the number of blocks in Lon and Lat (nBlocksLon, nBlocksLat)
!
! 4. the number of lon/lats in a single grid block  (nLons, nLats, nAlts)
!
! 5. the number of neutral species (and probably a mapping between
!    the species number and the variable name)  (nSpeciesTotal, nSpecies)
!
! 6. the number of ion species (ditto - numbers <-> names) (nIons)
!
! we assume that the 'UseTopography' flag is false - that all columns
! have the same altitude arrays.  this is true on earth but not on
! other planets.
!
! in addition to reading in the state data, it fills Longitude,
! Latitude, and Altitude arrays with the grid spacing.  this grid
! is orthogonal and rectangular but can have irregular spacing along
! any or all of the three dimensions.

character(len=*), intent(in)  :: dirname
real(r8),         intent(out) :: state_vector(:)
type(time_type),  intent(out) :: model_time

if ( .not. module_initialized ) call static_init_model

state_vector = MISSING_R8

! this is going to have to loop over all the blocks, both to get
! the data values and to get the full grid spacings.

model_time = get_state_time(dirname)

if (do_output()) &
    call print_time(model_time,'time in restart file '//trim(dirname)//'/header.rst')
if (do_output()) &
    call print_date(model_time,'date in restart file '//trim(dirname)//'/header.rst')

! sort the required fields into the order they exist in the
! binary restart files and fill in the state vector as you
! read each field.  when this routine returns all the data has
! been read.

call get_data(trim(dirname), state_vector)

end subroutine restart_file_to_statevector


!==================================================================


subroutine statevector_to_restart_file(state_vector, dirname, statedate)
!------------------------------------------------------------------
! Writes the current time and state variables from a dart state
! vector (1d array) into a gitm netcdf restart file.
!
real(r8),         intent(in) :: state_vector(:)
character(len=*), intent(in) :: dirname
type(time_type),  intent(in) :: statedate

character(len=256) :: dirnameout

if ( .not. module_initialized ) call static_init_model

print *, 'in statevector_to_restart_file, debug, nfields = ', debug, nfields

! sort the required fields into the order they exist in the
! binary restart files and write out the state vector data
! field by field.  when this routine returns all the data has
! been written.

dirnameout = trim(dirname) // '.out'
call put_data(dirname, dirnameout, state_vector)

if (do_output()) &
    call print_time(statedate,'time in restart file '//trim(dirname)//'/header.rst')
if (do_output()) &
    call print_date(statedate,'date in restart file '//trim(dirname)//'/header.rst')

end subroutine statevector_to_restart_file


!==================================================================


subroutine get_gitm_restart_dirname( dirname )
!------------------------------------------------------------------
character(len=*), intent(OUT) :: dirname

if ( .not. module_initialized ) call static_init_model

dirname = trim(gitm_restart_dirname)

end subroutine get_gitm_restart_dirname


!==================================================================


function get_state_time( dirname )
!------------------------------------------------------------------
! the static_init_model ensures that the gitm namelists are read.
!
type(time_type)              :: get_state_time
character(len=*), intent(in) :: dirname

type(time_type) :: model_offset, base_time

integer  :: iunit, i, ios
integer  :: istep
real(r8) :: tsimulation
integer  :: iyear, imonth, iday, ihour, imin, isec
integer  :: ndays,nsec

character(len=256) :: filename
character(len=100) :: cLine

if ( .not. module_initialized ) call static_init_model

tsimulation = MISSING_R8
istep       = -1
iyear       = -1
imonth      = -1
iday        = -1
ihour       = -1
imin        = -1
isec        = -1

write(filename,'(a,''/header.rst'')') trim(dirname)

iunit = open_file(trim(filename), action='read')

FILEREAD : do i = 1, 100

   read(iunit,'(a)',iostat=ios) cLine

   if (ios < 0) exit FILEREAD  ! end of file

   if (ios /= 0) then
      write(string1,*) 'cannot read ',trim(filename)
      call error_handler(E_ERR,'get_grid_info',string1,source,revision,revdate)
   endif

   select case( cLine(1:6) )
      case('#ISTEP')
         read(iunit,*)istep
      case('#TSIMU')
         read(iunit,*)tsimulation
      case('#TIMES')
         read(iunit,*)iyear
         read(iunit,*)imonth
         read(iunit,*)iday
         read(iunit,*)ihour
         read(iunit,*)imin
         read(iunit,*)isec
      case default
   end select

enddo FILEREAD

call close_file(iunit)

base_time      = set_date(iyear, imonth, iday, ihour, imin, isec)
ndays          = tsimulation/86400
nsec           = tsimulation - ndays*86400
model_offset   = set_time(nsec,ndays)
get_state_time = base_time + model_offset

if (debug > 8) then
   write(*,*)'get_state_time : iyear       ',iyear
   write(*,*)'get_state_time : imonth      ',imonth
   write(*,*)'get_state_time : iday        ',iday
   write(*,*)'get_state_time : ihour       ',ihour
   write(*,*)'get_state_time : imin        ',imin
   write(*,*)'get_state_time : isec        ',isec
   write(*,*)'get_state_time : tsimulation ',tsimulation
   write(*,*)'get_state_time : ndays       ',ndays
   write(*,*)'get_state_time : nsec        ',nsec

   call print_date(     base_time, 'get_state_time:model base date')
   call print_time(     base_time, 'get_state_time:model base time')
   call print_time(  model_offset, 'get_state_time:model offset')
   call print_date(get_state_time, 'get_state_time:model date')
   call print_time(get_state_time, 'get_state_time:model time')
endif

end function get_state_time


!==================================================================
! The remaining interfaces come last
!==================================================================


function get_grid_value(base_offset, ilon, ilat, ialt, state_handle, ens_size)
!------------------------------------------------------------------

integer(i8),         intent(in) :: base_offset
integer,             intent(in) :: ilon, ilat, ialt
type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
real(r8)                        :: get_grid_value(ens_size)

! Returns the value for the given lon,lat,alt point in the field that
! starts at offset base_offset

integer(i8) :: offset

offset = (ilon - 1) + (ilat - 1) * NgridLon + (ialt - 1) * (NgridLon * NgridLat)
get_grid_value = get_state(base_offset + offset, state_handle)

end function get_grid_value


!==================================================================


subroutine find_lon_bounds(llon, lower, upper, fract, ier)
!------------------------------------------------------------------

! Finds position of a given longitude in an array of longitude grid points and returns
! the index of the lower and upper bounds and the fractional offset. Assumes that the
! first longitude in the list is the smallest and that the largest is less than
! 360 degrees.  ier returns 0 unless there is an error.

real(r8), intent(in)  :: llon
integer,  intent(out) :: lower, upper
real(r8), intent(out) :: fract
integer,  intent(out) :: ier

! For now, assume that the spacing on longitudes is arbitrary.
! Do a silly linear search. Probably not worth any fancier searching unless
! models get to be huge.

integer  :: i
real(r8) :: width

if(llon < 0.0_r8 .or. llon > 360.0_r8) then
   ier = 2
   return
endif

! Look for case where longitude is less than smallest in list
if(llon <= LON(1)) then
   lower = NgridLon
   upper = 1
   width = 360.0_r8 - LON(NgridLon) + LON(1)
   fract = (llon + 360.0_r8 - LON(NgridLon)) / width
   ier = 0
   return
endif

! Look for case where longitude is greater than largest in list
if(llon >= LON(NgridLon)) then
  lower = NgridLon
  upper = 1
  width = 360.0 - LON(NgridLon) + LON(1)
  fract = (llon - LON(NgridLon)) / width
  ier = 0
  return
endif

! Otherwise in the interior
do i = 2, NgridLon
   if(llon < LON(i)) then
      lower = i - 1
      upper = i
      fract = (llon - LON(i-1)) / (LON(i) - LON(i - 1))
      ier = 0
      return
   endif
end do

! Shouldn't ever fall off end of loop
ier = 2

end subroutine find_lon_bounds


!==================================================================


subroutine find_lat_or_alt_bounds(llat, nbounds, bounds, lower, upper, fract, ier)
!------------------------------------------------------------------
! Finds position of a given latitude in an array of latitude grid points and returns
! the index of the lower and upper bounds and the fractional offset. Used for both
! latitude and altitude which have similar linear arrays. ier returns 0 unless there
! is an error.

real(r8), intent(in)  :: llat
integer,  intent(in)  :: nbounds
real(r8), intent(in)  :: bounds(nbounds)
integer,  intent(out) :: lower, upper
real(r8), intent(out) :: fract
integer,  intent(out) :: ier

! For now, assume that the spacing on latitudes or altitudes is arbitrary
! Do a silly linear search. Probably not worth any fancier searching unless
! models get to be huge.

integer :: i

if(llat < bounds(1) .or. llat > bounds(nbounds)) then
   ier = 2
   return
endif

do i = 2, nbounds
   if(llat <= bounds(i)) then
      lower = i - 1
      upper = i
      fract = (llat - bounds(i-1)) / (bounds(i) - bounds(i - 1))
      ier = 0
      return
   endif
end do

! Shouldn't ever fall off end of loop
ier = 2

end subroutine find_lat_or_alt_bounds


!==================================================================


subroutine get_grid_info(NgridLon, NgridLat, NgridAlt, &
                nBlocksLon, nBlocksLat, LatStart, LatEnd, LonStart)
!------------------------------------------------------------------
!
! Read the grid dimensions from the restart netcdf file.
!
! The file name comes from module storage ... namelist.

integer,  intent(out) :: NgridLon   ! Number of Longitude centers
integer,  intent(out) :: NgridLat   ! Number of Latitude  centers
integer,  intent(out) :: NgridAlt   ! Number of Vertical grid centers
integer,  intent(out) :: nBlocksLon, nBlocksLat
real(r8), intent(out) :: LatStart, LatEnd, LonStart

character(len=10) :: filename = 'UAM.in'

character(len=100) :: cLine  ! iCharLen_ == 100
character(len=256) :: fileloc

integer :: i, iunit, ios

if ( .not. module_initialized ) call static_init_model

! get the ball rolling ...

nBlocksLon = 0
nBlocksLat = 0
LatStart   = 0.0_r8
LatEnd     = 0.0_r8
LonStart   = 0.0_r8

write(fileloc,'(a,''/'',a)') trim(gitm_restart_dirname),trim(filename)

iunit = open_file(trim(fileloc), action='read')

UAMREAD : do i = 1, 1000000

   read(iunit,'(a)',iostat=ios) cLine

   if (ios /= 0) then
      ! If we get to the end of the file or hit a read error without
      ! finding what we need, die.
      write(string1,*) 'cannot find #GRID in ',trim(fileloc)
      call error_handler(E_ERR,'get_grid_info',string1,source,revision,revdate)
   endif

   if (cLine(1:5) .ne. "#GRID") cycle UAMREAD

   nBlocksLon = read_in_int( iunit,'NBlocksLon',trim(fileloc))
   nBlocksLat = read_in_int( iunit,'NBlocksLat',trim(fileloc))
   LatStart   = read_in_real(iunit,'LatStart',  trim(fileloc))
   LatEnd     = read_in_real(iunit,'LatEnd',    trim(fileloc))
   LonStart   = read_in_real(iunit,'LonStart',  trim(fileloc))

   exit UAMREAD

enddo UAMREAD

call close_file(iunit)

NgridLon = nBlocksLon * nLons
NgridLat = nBlocksLat * nLats
NgridAlt = nAlts

end subroutine get_grid_info


!==================================================================


subroutine get_grid(dirname, nBlocksLon, nBlocksLat, &
                  nLons, nLats, nAlts, LON, LAT, ALT )
!------------------------------------------------------------------
! open enough of the restart files to read in the lon, lat, alt arrays
!
character(len=*), intent(in) :: dirname
integer, intent(in) :: nBlocksLon ! Number of Longitude blocks
integer, intent(in) :: nBlocksLat ! Number of Latitude  blocks
integer, intent(in) :: NLons      ! Number of Longitude centers per block
integer, intent(in) :: NLats      ! Number of Latitude  centers per block
integer, intent(in) :: NAlts      ! Number of Vertical grid centers

real(r8), dimension( : ), intent(inout) :: LON, LAT, ALT

integer :: ios, nb, offset, iunit, nboff
character(len=256) :: filename
real(r8), allocatable :: temp(:)

! a temp array large enough to hold any of the
! Lon,Lat or Alt array from a block plus ghost cells
allocate(temp(1-nGhost:max(nLons,nLats,nAlts)+nGhost))

! get the dirname, construct the filenames inside

if (debug > 9) then
   write(*,*)'TJH DEBUG: get_grid : size(temp)   is ',size(temp)
   write(*,*)'TJH DEBUG: get_grid : nBlocksLon   is ',nBlocksLon
   write(*,*)'TJH DEBUG: get_grid : nLons,nGhost is ',nLons,nGhost
endif

! go across the south-most block row picking up all longitudes
do nb = 1, nBlocksLon

   iunit = open_block_file(dirname, nb, 'read', filename)

   read(iunit,iostat=ios) temp(1-nGhost:nLons+nGhost)
   if ( ios /= 0 ) then
      write(string1,*)'ERROR reading file ', trim(filename)
      write(string2,*)'longitude block ',nb,' of ',nBlocksLon
      call error_handler(E_ERR,'get_grid',string1, &
                 source,revision,revdate,text2=string2)
   endif

   offset = (nLons * (nb - 1))
   LON(offset+1:offset+nLons) = temp(1:nLons)

   call close_file(iunit)
enddo

! go up west-most block row picking up all latitudes
do nb = 1, nBlocksLat

   nboff = ((nb - 1) * nBlocksLon) + 1
   iunit = open_block_file(dirname, nboff, 'read', filename)

   ! get past lon array and read in lats
   read(iunit) temp(1-nGhost:nLons+nGhost)

   read(iunit,iostat=ios) temp(1-nGhost:nLats+nGhost)
   if ( ios /= 0 ) then
      write(string1,*)'ERROR reading file ', trim(filename)
      write(string2,*)'latitude block ',nb,' of ',nBlocksLat
      call error_handler(E_ERR,'get_grid',string1, &
                 source,revision,revdate,text2=string2)
   endif

   offset = (nLats * (nb - 1))
   LAT(offset+1:offset+nLats) = temp(1:nLats)

   call close_file(iunit)
enddo

! this code assumes UseTopography is false - that all columns share
! the same altitude array, so we can read it from the first block.
! if this is not the case, this code has to change.

iunit = open_block_file(dirname, 1, 'read', filename)

! get past lon and lat arrays and read in alt array
read(iunit) temp(1-nGhost:nLons+nGhost)
read(iunit) temp(1-nGhost:nLats+nGhost)
read(iunit) temp(1-nGhost:nAlts+nGhost)

ALT(1:nAlts) = temp(1:nAlts)

call close_file(iunit)

deallocate(temp)

! convert from radians into degrees
LON = LON * rad2deg
LAT = LAT * rad2deg

if (debug > 4) then
   print *, 'All LONs ', LON
   print *, 'All LATs ', LAT
   print *, 'All ALTs ', ALT
endif

if ( debug > 1 ) then ! A little sanity check
   write(*,*)'LON range ',minval(LON),maxval(LON)
   write(*,*)'LAT range ',minval(LAT),maxval(LAT)
   write(*,*)'ALT range ',minval(ALT),maxval(ALT)
endif

end subroutine get_grid


!==================================================================


function open_block_file(dirname, blocknum, rw, filename)
!------------------------------------------------------------------
! open the requested block number restart file and return the
! open file unit

integer                       :: open_block_file
character(len=*), intent(in)  :: dirname
integer,          intent(in)  :: blocknum
character(len=*), intent(in)  :: rw   ! 'read' or 'readwrite'
character(len=*), intent(out) :: filename

write(filename, '(A,i4.4,A)') trim(dirname)//'/b', blocknum, '.rst'

if ( rw == 'read' .and. .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'open_block_file',string1,source,revision,revdate)
endif

!print *, 'opening file ', trim(filename), ' for ', trim(rw)

open_block_file = open_file(filename, 'unformatted', rw)

!print *, 'returned file descriptor is ', open_block_file

end function open_block_file


!==================================================================


subroutine get_data(dirname, statevector)
!------------------------------------------------------------------
! open all restart files and read in the requested data item
!
character(len=*), intent(in)  :: dirname
real(r8),         intent(out) :: statevector(:)

integer :: ib, jb, nb, iunit, blockoffset, i

character(len=256) :: filename

! get the dirname, construct the filenames inside open_block_file

do jb = 1, nBlocksLat
do ib = 1, nBlocksLon

   nb = (jb-1) * nBlocksLon + ib

   blockoffset = nLats * ngridLon * (jb-1) + nLons * (ib-1)

!print *, 'ib,jb = ', ib, jb
!print *, 'blockoffset, nb = ', blockoffset, nb

   iunit = open_block_file(dirname, nb, 'read', filename)

   call read_data(iunit, blockoffset, statevector)

   call close_file(iunit)
enddo
enddo

if ( debug > 4 ) then ! A little sanity check
   write (*,*) 'variable data after read: '
   do i = 1, nfields
      write(*,*) trim(progvar(i)%varname), ' range ', &
                 minval(statevector(progvar(i)%index1:progvar(i)%indexN)), &
                 maxval(statevector(progvar(i)%index1:progvar(i)%indexN))
   enddo
endif

end subroutine get_data


!==================================================================


subroutine put_data(dirname, dirnameout, statevector)
!------------------------------------------------------------------
! open all restart files and write out the requested data item
!
 character(len=*), intent(in) :: dirname, dirnameout
 real(r8),         intent(in) :: statevector(:)

integer :: ib, jb, nb, iunit, ounit
integer :: i, blockoffset
character(len=256) :: readfilename, writefilename

! get the dirname, construct the filenames inside open_block_file

if ( debug > 4 ) then ! A little sanity check
   write (*,*) 'variable data to be written: '
   do i = 1, nfields
      write(*,*) trim(progvar(i)%varname), 'range ', &
                 minval(statevector(progvar(i)%index1:progvar(i)%indexN)), &
                 maxval(statevector(progvar(i)%index1:progvar(i)%indexN))
   enddo
endif

do jb = 1, nBlocksLat
do ib = 1, nBlocksLon

   nb = (jb-1) * nBlocksLon + ib

   blockoffset = nLats * ngridLon * (jb-1) + &
                 nLons * (ib-1)

!print *, 'ib,jb = ', ib, jb
!print *, 'blockoffset, nb = ', blockoffset, nb

   iunit = open_block_file(dirname,    nb, 'read',  readfilename)
   ounit = open_block_file(dirnameout, nb, 'write', writefilename)

   call write_data(iunit, ounit, blockoffset, statevector, readfilename, writefilename)

   call close_file(iunit)
   call close_file(ounit)
enddo
enddo

end subroutine put_data


!==================================================================


subroutine unpack_data(data3d, ivar, blockoffset, statevector)
!------------------------------------------------------------------
! put the requested data into the state vector
!
real(r8), intent(in)    :: data3d(:,:,:)
integer,  intent(in)    :: ivar         ! index into progvar struct
integer,  intent(in)    :: blockoffset
real(r8), intent(inout) :: statevector(:)

integer :: i, j, k, offset, base

!print *, 'ivar = ', ivar
base = progvar(ivar)%index1 - 1
!print *, 'blockoffset, base = ', blockoffset, base

do k=1,nAlts
 do j=1,nLats
  do i=1,nLons

      offset = ((k-1) * ngridLat * ngridLon) +  &
               ((j-1) * ngridLon) +             &
               i
    if (base+blockoffset+offset < 1 .or. &
        base+blockoffset+offset > model_size) then
      print *, 'i,j,k, index: ', i, j, k, base+blockoffset+offset
    else
      statevector(base + blockoffset + offset) = data3d(nGhost+i, nGhost+j, nGhost+k)
      !print *, 'i,j,k,varoffset = ', i,j,k,blockoffset + offset
    endif

  enddo
 enddo
enddo

end subroutine unpack_data


!==================================================================


subroutine unpack_data0d(data0d, ivar, blockoffset, statevector)
!------------------------------------------------------------------
! put the f107 estimate (a scalar, hence 0d) into the state vector.
! Written specifically
! for f107 since f107 is the same for all blocks. So what it does
! is take f107 from the first block (blockoffset = 0) and disregard
! f107 values from all other blocks (hopefully they are the same).
! written by alex

real(r8), intent(in)    :: data0d
integer,  intent(in)    :: ivar         ! index into progvar struct
integer,  intent(in)    :: blockoffset
real(r8), intent(inout) :: statevector(:)

integer :: offset, base

base = progvar(ivar)%index1 - 1
offset = 1

if (blockoffset > 0) then
   ! if not the first block, don't put this value into the state vector
   ! print *, 'u: BO > 0, NOT updating SV, throwing this f107 value AWAY! ', &
   ! base, blockoffset, offset, data0d
else
   ! if the first block (blockoffset = 0), then put this value into the state vector
   ! print *, 'u: BASE+BO+O is fine', base, blockoffset, offset, data0d
   ! blockoffset is 0 (f107 does not depend on blocks), offset is 1
   statevector(base + blockoffset + offset) = data0d
endif

end subroutine unpack_data0d


!==================================================================


subroutine pack_data(statevector, ivar, blockoffset, data3d)
!------------------------------------------------------------------
! put the state vector data into a 3d array
!
real(r8), intent(in)    :: statevector(:)
integer,  intent(in)    :: ivar         ! index into progvar struct
integer,  intent(in)    :: blockoffset
real(r8), intent(inout) :: data3d(:,:,:)

integer :: i, j, k, offset, base

base = progvar(ivar)%index1 - 1

do k=1,nAlts
 do j=1,nLats
  do i=1,nLons

      offset = ((k-1) * ngridLat * ngridLon) +  &
               ((j-1) * ngridLon) +             &
               i
      data3d(nGhost+i, nGhost+j, nGhost+k) = statevector(base + blockoffset + offset)
      !print *, 'i,j,k,varoffset = ', i,j,k,blockoffset + offset

  enddo
 enddo
enddo

end subroutine pack_data


!==================================================================


subroutine pack_data0d(statevector, ivar, blockoffset, data0d)
!------------------------------------------------------------------
! put the f107 estimate (scalar) from the statevector into a 0d container
! the only trick this routine does is give all blocks the same f107 (the
! f107 value from block 1 state vector goes to block 1,2,3,4 restart files)
! so no matter what, always grab the f107 from block 1 (manipulate
! the blockoffset variable).
! written by alex

real(r8), intent(in)    :: statevector(:)
integer,  intent(in)    :: ivar         ! index into progvar struct
integer,  intent(in)    :: blockoffset
real(r8), intent(inout) :: data0d

integer :: offset, base

base   = progvar(ivar)%index1 - 1
offset = 1

if (blockoffset > 0) then
   ! block > 1
   ! print *, 'p: BO>0, updating data0d w f107 from BLOCK 1 (ONE) !!!', &
   !           base, blockoffset, offset, statevector(base + 0 + offset)
   ! 0 blockoffset corresponds to block 1
   data0d = statevector(base + 0 + offset)
else
   ! block = 1
   ! print *, 'p: BASE+BO+O is fine', &
   !           base, blockoffset, offset, statevector(base + blockoffset + offset)
   data0d = statevector(base + blockoffset + offset)
endif

end subroutine pack_data0d


!==================================================================


subroutine read_data(iunit, blockoffset, statevector)
!------------------------------------------------------------------
! open all restart files and read in the requested data items
!
integer,  intent(in)    :: iunit
integer,  intent(in)    :: blockoffset
real(r8), intent(inout) :: statevector(:)

real(r8), allocatable :: temp1d(:), temp3d(:,:,:), temp4d(:,:,:,:)
real(r8) :: temp0d !Alex: single parameter has "zero dimensions"
integer :: i, j, inum, maxsize, ivals(NSpeciesTotal)

! a temp array large enough to hold any of the
! Lon,Lat or Alt array from a block plus ghost cells
allocate(temp1d(1-nGhost:max(nLons,nLats,nAlts)+nGhost))

! temp array large enough to hold 1 species, temperature, etc
allocate(temp3d(1-nGhost:nLons+nGhost, &
                1-nGhost:nLats+nGhost, &
                1-nGhost:nAlts+nGhost))

! temp array large enough to hold velocity vect, etc
maxsize = max(3, nSpecies)
allocate(temp4d(1-nGhost:nLons+nGhost, &
                1-nGhost:nLats+nGhost, &
                1-nGhost:nAlts+nGhost, maxsize))

! get past lon and lat arrays and read in alt array
read(iunit) temp1d(1-nGhost:nLons+nGhost)
read(iunit) temp1d(1-nGhost:nLats+nGhost)
read(iunit) temp1d(1-nGhost:nAlts+nGhost)

call get_index_from_gitm_varname('NDensityS', inum, ivals)
if (inum > 0) then
   ! if i equals ival, use the data from the state vect
   ! otherwise read/write what's in the input file
   j = 1
   do i = 1, nSpeciesTotal
      read(iunit)  temp3d
      if (j <= inum) then
         if (i == progvar(ivals(j))%gitm_index) then
            call unpack_data(temp3d, ivals(j), blockoffset, statevector)
            j = j + 1
         endif
      endif
   enddo
else
   ! nothing at all from this variable in the state vector.
   ! copy all data over from the input file to output file
   do i = 1, nSpeciesTotal
      read(iunit)  temp3d
   enddo
endif

call get_index_from_gitm_varname('IDensityS', inum, ivals)
if (inum > 0) then
   ! one or more items in the state vector need to replace the
   ! data in the output file.  loop over the index list in order.
   j = 1
   do i = 1, nIons
      read(iunit)  temp3d
      if (j <= inum) then
         if (i == progvar(ivals(j))%gitm_index) then
            ! read from input but write from state vector
            call unpack_data(temp3d, ivals(j), blockoffset, statevector)
            j = j + 1
         endif
      endif
   enddo
else
   ! nothing at all from this variable in the state vector.
   ! read past this variable
   do i = 1, nIons
      read(iunit)  temp3d
   enddo
endif

read(iunit)  temp3d
call get_index_from_gitm_varname('Temperature', inum, ivals)
if (inum > 0) then
   call unpack_data(temp3d, ivals(1), blockoffset, statevector)
endif

read(iunit) temp3d
call get_index_from_gitm_varname('ITemperature', inum, ivals)
if (inum > 0) then
   call unpack_data(temp3d, ivals(1), blockoffset, statevector)
endif

read(iunit) temp3d
call get_index_from_gitm_varname('eTemperature', inum, ivals)
if (inum > 0) then
   call unpack_data(temp3d, ivals(1), blockoffset, statevector)
endif

!print *, 'reading in temp4d for vel'
read(iunit) temp4d(:,:,:,1:3)
call get_index_from_gitm_varname('Velocity', inum, ivals)
if (inum > 0) then
   ! copy out any requested bits into state vector
   j = 1
   do i = 1, 3
      if (j <= inum) then
         if (i == progvar(ivals(j))%gitm_index) then
            temp3d = temp4d(:,:,:,i)
            call unpack_data(temp3d, ivals(j), blockoffset, statevector)
            j = j + 1
         endif
      endif
   enddo
endif

!print *, 'reading in temp4d for ivel'
read(iunit) temp4d(:,:,:,1:3)
call get_index_from_gitm_varname('IVelocity', inum, ivals)
if (inum > 0) then
   ! copy out any requested bits into state vector
   j = 1
   do i = 1, 3
      if (j <= inum) then
         if (i == progvar(ivals(j))%gitm_index) then
            ! read from input but write from state vector
            temp3d = temp4d(:,:,:,i)
            call unpack_data(temp3d, ivals(j), blockoffset, statevector)
            j = j + 1
         endif
      endif
   enddo
endif

!print *, 'reading in temp4d for vvel'
read(iunit) temp4d(:,:,:,1:nSpecies)
call get_index_from_gitm_varname('VerticalVelocity', inum, ivals)
if (inum > 0) then
   ! copy out any requested bits into state vector
   j = 1
   do i = 1, nSpecies
      if (j <= inum) then
         if (i == progvar(ivals(j))%gitm_index) then
            temp3d = temp4d(:,:,:,i)
            call unpack_data(temp3d, ivals(j), blockoffset, statevector)
            j = j + 1
         endif
      endif
   enddo
endif

!alex begin
read(iunit)  temp0d
call get_index_from_gitm_varname('f107', inum, ivals)
if (inum > 0) then
   call unpack_data0d(temp0d, ivals(1), blockoffset, statevector) !see comments in the body of the subroutine
endif

read(iunit)  temp3d
call get_index_from_gitm_varname('Rho', inum, ivals)
if (inum > 0) then
   call unpack_data(temp3d, ivals(1), blockoffset, statevector)
endif
!alex end

!print *, 'calling dealloc'
deallocate(temp1d, temp3d, temp4d)

end subroutine read_data


!==================================================================


subroutine write_data(iunit, ounit, blockoffset, statevector, infile, outfile)
!------------------------------------------------------------------
! open all restart files and write out the requested data item
!
integer,          intent(in) :: iunit, ounit
integer,          intent(in) :: blockoffset
real(r8),         intent(in) :: statevector(:)
character(len=*), intent(in) :: infile, outfile

real(r8), allocatable :: temp1d(:), temp3d(:,:,:), temp4d(:,:,:,:), data3d(:,:,:)
real(r8) :: data0d, temp0d !Alex !parameter is technically zero-dimensional
integer :: ios
integer :: i, j, inum, maxsize, ivals(NSpeciesTotal)

! a temp array large enough to hold any of the
! Lon,Lat or Alt array from a block plus ghost cells
allocate(temp1d(1-nGhost:max(nLons,nLats,nAlts)+nGhost))

! temp array large enough to hold 1 species, temperature, etc
allocate(temp3d(1-nGhost:nLons+nGhost, 1-nGhost:nLats+nGhost, 1-nGhost:nAlts+nGhost))
allocate(data3d(1-nGhost:nLons+nGhost, 1-nGhost:nLats+nGhost, 1-nGhost:nAlts+nGhost))

! temp array large enough to hold velocity vect, etc
maxsize = max(3, nSpecies)
allocate(temp4d(1-nGhost:nLons+nGhost, 1-nGhost:nLats+nGhost, 1-nGhost:nAlts+nGhost, maxsize))

! copy over lat, lon, alt arrays verbatim
read(iunit,iostat=ios) temp1d(1-nGhost:nLons+nGhost)
if (ios /= 0) then
   write(string1,*)'unable to read lons from ',trim(infile)
   call error_handler(E_ERR,'write_data',string1,source,revision,revdate)
endif
write(ounit,iostat=ios) temp1d(1-nGhost:nLons+nGhost)
if (ios /= 0) then
   write(string1,*)'unable to write lons to ',trim(outfile)
   call error_handler(E_ERR,'write_data',string1,source,revision,revdate)
endif

read(iunit,iostat=ios) temp1d(1-nGhost:nLats+nGhost)
if (ios /= 0) then
   write(string1,*)'unable to read lats from ',trim(infile)
   call error_handler(E_ERR,'write_data',string1,source,revision,revdate)
endif
write(ounit,iostat=ios) temp1d(1-nGhost:nLats+nGhost)
if (ios /= 0) then
   write(string1,*)'unable to write lats to ',trim(outfile)
   call error_handler(E_ERR,'write_data',string1,source,revision,revdate)
endif

read(iunit,iostat=ios) temp1d(1-nGhost:nAlts+nGhost)
if (ios /= 0) then
   write(string1,*)'unable to read alts from ',trim(infile)
   call error_handler(E_ERR,'write_data',string1,source,revision,revdate)
endif
write(ounit,iostat=ios) temp1d(1-nGhost:nAlts+nGhost)
if (ios /= 0) then
   write(string1,*)'unable to write alts to ',trim(infile)
   call error_handler(E_ERR,'write_data',string1,source,revision,revdate)
endif

call get_index_from_gitm_varname('NDensityS', inum, ivals)
if (inum > 0) then
   ! if i equals ival, use the data from the state vect
   ! otherwise read/write what's in the input file
   j = 1
   do i = 1, nSpeciesTotal
      read(iunit)  temp3d
      if (j <= inum) then
         if (i == progvar(ivals(j))%gitm_index) then

            ! FIXME: if the program restart is really resetting the ghost zones
            ! correctly, then we shouldn't need to initialize the array with the
            ! temp3d data (which has pre-assimilation values in it).  but alexey
            ! says this causes problems, which is suspicious and should be looked
            ! at more.  this line might make it run, but the ghost zones were not
            ! updated by the assimilation.
            !alex: the horizontal ghost cells at middle altitudes (just not in the
            !altitude top and bottom ghost cells) should be fine as they get overwritten
            !in GITM. The top and bottom altitude ghost cells are the culprits.
            !Ideally, they should be extrapolated to once DART provides the posterior estimates,
            !but this is not implemented yet. This lack should not affect the
            !assimilation too much if the observations come from middle altitudes
            ! (as is the case with CHAMP and GRACE).
            data3d = temp3d

            call pack_data(statevector, ivals(j), blockoffset, data3d)

            ! FIXME: also needs fixing.  if we have made some value negative
            ! where the model doesn't support it, and it can't be 0 either, then
            ! this should be the smallest positive value that the model will accept.
            ! the original data divided by 2 is going to change the distribution of
            ! values and is certainly not right.  leave it here for now to get the
            ! assimilation running, but this needs looking at and changing soon.
            !alex: fixed on 5/20/13. How? Well, the limits of the variables are as
            !follows (taken from a GITM initialization on 12/1/2002 via gitm/matlab/rst2mat.m):
! MINIMA AND MAXIMA
! _
! LonT 0.17453 6.1087
! LatT -1.4835 1.4835
! AltT 100000 630038.9261
! TempT                  163.0163 1223.5239
! ITempT                 163.0154 1967.9977
! eTempT                 184.665 2710.9351
! _
! NDST_1,iO_3P_ 607871671694.2802 624953710511309568
! NDST_2,iO2_        1554285.5124 2977090934472271872
! NDST_3,iN2_      261275675.713  12920995180058857472
! NDST_4,iN_4S_   2725259865.9408 51174404943040.94
! NDST_5,iNO_             91.5983 137340019620842.8
! NDST_6,iN_2D_       490627.9878 656673564758.9565
! NDST_7,iN_2P_       135692.3204 2582963359.5952
! NDST_8,iH_    297189394877.7289 160285753329765.5
! NDST_9,iHe_    39396601335.7323 31530483811658.93
! NDST_10,iCO2_           51.3449 5237123737981628
! NDST_11,iO_1D_          32.2011 26604279344.3065
! _
! IDST_1,iO_4SP_         100 2345587357569.55
! IDST_2,iO2P_             4.0622 121043204145.427
! IDST_3,iN2P_             2.3259e-05 6408254083.7036
! IDST_4,iNP_              1.6073e-05 725487968.9667
! IDST_5,iNOP_            15.9515 182204005544.7968
! IDST_6,iO_2DP_           2.6996e-11 798313237.9133
! IDST_7,iO_2PP_           9.5018e-11 365561613.5574
! IDST_8,iHP_              1 250583438981.8537
! IDST_9,iHeP_             1 13445167727.3174
! IDST_10,ie_      543075134.2391 2346712140512.865
! _
! VT -253.1168 236.8601
! IVT -809.0201 1382.4808
! VVT -82.1731 633.1406
! RhoT 1.6352e-14 7.7801e-07
            !
            !So it makes sense to saturate NDST, IDST and Rho by 1.0e-16 from below,
            ! Temp, ITemp and eTemp by 100.0 from below, and F107 by 60.0 from below

!            where (data3d < 0.0_r8) data3d = temp3d/2 !alex, old - bad because might change distr
            where (data3d < 0.0_r8) data3d = 1.0e-16_r8 !alex, new

            write(ounit) data3d
            j = j + 1
         else
            write(ounit) temp3d
         endif
      else
         ! this one not in state vector, copy over from input
         write(ounit) temp3d
      endif
   enddo
else
   ! nothing at all from this variable in the state vector.
   ! copy all data over from the input file to output file
   do i = 1, nSpeciesTotal
      read(iunit)  temp3d
      write(ounit) temp3d
   enddo
endif

call get_index_from_gitm_varname('IDensityS', inum, ivals)
if (inum > 0) then
   ! one or more items in the state vector need to replace the
   ! data in the output file.  loop over the index list in order.
   j = 1
   do i = 1, nIons
      read(iunit)  temp3d
      if (j <= inum) then
         if (i == progvar(ivals(j))%gitm_index) then
            ! read from input but write from state vector
            data3d = temp3d
            call pack_data(statevector, ivals(j), blockoffset, data3d)
            where (data3d < 0.0_r8) data3d = 1.0e-16_r8 !alex
            write(ounit) data3d
            j = j + 1
         else
            write(ounit) temp3d
         endif
      else
         ! this one not in state vector, copy over from input
         write(ounit) temp3d
      endif
   enddo
else
   ! nothing at all from this variable in the state vector.
   ! copy all data over from the input file to output file
   do i = 1, nIons
      read(iunit)  temp3d
      write(ounit) temp3d
   enddo
endif

read(iunit)  temp3d
data3d = temp3d
call get_index_from_gitm_varname('Temperature', inum, ivals)
if (inum > 0) then
   call pack_data(statevector, ivals(1), blockoffset, data3d)
   where (data3d < 0.0_r8) data3d = 100.0_r8 !alex
   write(ounit) data3d
else
   write(ounit) temp3d
endif


read(iunit) temp3d
data3d = temp3d
call get_index_from_gitm_varname('ITemperature', inum, ivals)
if (inum > 0) then
   call pack_data(statevector, ivals(1), blockoffset, data3d)
   where (data3d < 0.0_r8) data3d = 100.0_r8 !alex
   write(ounit) data3d
else
   write(ounit) temp3d
endif

read(iunit) temp3d
data3d = temp3d
call get_index_from_gitm_varname('eTemperature', inum, ivals)
if (inum > 0) then
   call pack_data(statevector, ivals(1), blockoffset, data3d)
   where (data3d < 0.0_r8) data3d = 100.0_r8 !alex
   write(ounit) data3d
else
   write(ounit) temp3d
endif

!print *, 'reading in temp4d for vel'
read(iunit) temp4d(:,:,:,1:3)
call get_index_from_gitm_varname('Velocity', inum, ivals)
if (inum > 0) then
   ! one or more items in the state vector need to replace the
   ! data in the output file.  loop over the index list in order.
   j = 1
   do i = 1, 3
      if (j <= inum) then
         if (i == progvar(ivals(j))%gitm_index) then
            ! read from input but write from state vector
            data3d = temp4d(:,:,:,i)
            call pack_data(statevector, ivals(j), blockoffset, data3d)
            temp4d(:,:,:,i) = data3d
            j = j + 1
         endif
      endif
   enddo
endif
write(ounit) temp4d(:,:,:,1:3)

!print *, 'reading in temp4d for ivel'
read(iunit) temp4d(:,:,:,1:3)
call get_index_from_gitm_varname('IVelocity', inum, ivals)
if (inum > 0) then
   ! one or more items in the state vector need to replace the
   ! data in the output file.  loop over the index list in order.
   j = 1
   do i = 1, 3
      if (j <= inum) then
         if (i == progvar(ivals(j))%gitm_index) then
            ! read from input but write from state vector
            data3d = temp4d(:,:,:,i)
            call pack_data(statevector, ivals(j), blockoffset, data3d)
            temp4d(:,:,:,i) = data3d
            j = j + 1
         endif
      endif
   enddo
endif
write(ounit) temp4d(:,:,:,1:3)

!print *, 'reading in temp4d for vvel'
read(iunit) temp4d(:,:,:,1:nSpecies)
call get_index_from_gitm_varname('VerticalVelocity', inum, ivals)
if (inum > 0) then
   ! one or more items in the state vector need to replace the
   ! data in the output file.  loop over the index list in order.
   j = 1
   do i = 1, nSpecies
      if (j <= inum) then
         if (i == progvar(ivals(j))%gitm_index) then
            ! read from input but write from state vector
            data3d = temp4d(:,:,:,i)
            call pack_data(statevector, ivals(j), blockoffset, data3d)
            temp4d(:,:,:,i) = data3d
            j = j + 1
         endif
      endif
   enddo
endif
write(ounit) temp4d(:,:,:,1:nSpecies)


!alex begin: added f107 and Rho to the restart files:
read(iunit) temp0d
data0d = temp0d
call get_index_from_gitm_varname('f107', inum, ivals)
if (inum > 0) then
   call pack_data0d(statevector, ivals(1), blockoffset, data0d)
   if (data0d < 0.0_r8) data0d = 60.0_r8 !alex
   write(ounit) data0d
else
   write(ounit) temp0d
endif

read(iunit)  temp3d
data3d = temp3d
call get_index_from_gitm_varname('Rho', inum, ivals)
if (inum > 0) then
   call pack_data(statevector, ivals(1), blockoffset, data3d)
   where (data3d < 0.0_r8) data3d = 1.0e-16_r8 !alex
   write(ounit) data3d
else
   write(ounit) temp3d
endif
!alex end

deallocate(temp1d, temp3d, temp4d, data3d)

end subroutine write_data


!==================================================================


subroutine get_index_range_string(string,index1,indexN)
!------------------------------------------------------------------
! Determine where a particular DART kind (string) exists in the
! DART state vector.

character(len=*), intent(in)  :: string
integer,          intent(out) :: index1,indexN

integer :: i

index1 = 0
indexN = 0

FieldLoop : do i=1,nfields
   if (progvar(i)%kind_string /= trim(string)) cycle FieldLoop
   index1 = progvar(i)%index1
   indexN = progvar(i)%indexN
   exit FieldLoop
enddo FieldLoop

! FIXME: generally we don't halt the program when
! asking for an obs kind that isn't in the state vector.
! we just return an error code from the interpolate() call.
! to make that happen, comment this block out and just
! let it return with the indices set to 0.
if ((index1 == 0) .or. (indexN == 0)) then
   write(string1,*) 'Problem, cannot find indices for '//trim(string)
   call error_handler(E_ERR,'get_index_range_string',string1,source,revision,revdate)
endif

end subroutine get_index_range_string


!==================================================================


subroutine get_index_range_int(dartkind,index1,indexN)
!------------------------------------------------------------------
! Determine where a particular DART kind (integer) exists in the
! DART state vector.

integer, intent(in) :: dartkind
integer(i8), intent(out) :: index1,indexN

integer :: i
character(len=obstypelength) :: string

index1 = 0
indexN = 0

FieldLoop : do i=1,nfields
   if (progvar(i)%dart_kind /= dartkind) cycle FieldLoop
   index1 = progvar(i)%index1
   indexN = progvar(i)%indexN
   exit FieldLoop
enddo FieldLoop

string = get_name_for_quantity(dartkind)

if ((index1 == 0) .or. (indexN == 0)) then
   write(string1,*) 'Problem, cannot find indices for kind ',dartkind,trim(string)
   call error_handler(E_ERR,'get_index_range_int',string1,source,revision,revdate)
endif

end subroutine get_index_range_int


!==================================================================


subroutine get_index_from_gitm_varname(gitm_varname, inum, ivals)
!------------------------------------------------------------------
! Determine where any data from a given gitm_varname lies in the
! DART state vector.

character(len=*), intent(in) :: gitm_varname
integer, intent(out) :: inum, ivals(:)

integer :: gindex(nfields)
integer :: i, limit

inum = 0
limit = size(ivals)

FieldLoop : do i=1,nfields
   if (progvar(i)%gitm_varname /= gitm_varname) cycle FieldLoop
   inum = inum + 1
   if (inum > limit) then
      write(string1,*) 'found too many matches, ivals needs to be larger than ', limit
      call error_handler(E_ERR,'get_index_from_gitm_varname',string1,source,revision,revdate)
   endif
   ! i is index into progvar array - the order of the fields in the sv
   ! gitm_index is index into the specific variable in the gitm restarts
   ivals(inum) = i
   gindex(inum) = progvar(i)%gitm_index
enddo FieldLoop

!if (inum > 0) then
!   print *, 'before sort, inum: ', inum
!   print *, 'before sort, gindex: ', gindex(1:inum)
!   print *, 'before sort, ivals: ', ivals(1:inum)
!endif

! return the vals sorted by gitm_index order if more than 1
if (inum > 1) call sortindexlist(gindex, ivals, inum)

!if (inum > 0) then
!   print *, 'after  sort, inum: ', inum
!   print *, 'after  sort, gindex: ', gindex(1:inum)
!   print *, 'after  sort, ivals: ', ivals(1:inum)
!endif


end subroutine get_index_from_gitm_varname


!==================================================================


function set_model_time_step()
!------------------------------------------------------------------
! the static_init_model ensures that the gitm namelists are read.
!
type(time_type) :: set_model_time_step

if ( .not. module_initialized ) call static_init_model

   set_model_time_step = set_time(assimilation_period_seconds, &
                           assimilation_period_days) ! (seconds, days)

end function set_model_time_step


!==================================================================


subroutine verify_state_variables( state_variables, ngood, table )
!------------------------------------------------------------------
character(len=*), dimension(:),   intent(in)  :: state_variables
integer,                          intent(out) :: ngood
character(len=*), dimension(:,:), intent(out) :: table

integer :: nrows, i
character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME) :: dartstr

if ( .not. module_initialized ) call static_init_model

nrows = size(table,1)

ngood = 0
MyLoop : do i = 1, nrows

   varname    = trim(state_variables(2*i -1))
   dartstr    = trim(state_variables(2*i   ))
   table(i,1) = trim(varname)
   table(i,2) = trim(dartstr)

   if ( table(i,1) == ' ' .and. table(i,2) == ' ' ) exit MyLoop ! Found end of list.

   if ( table(i,1) == ' ' .or. table(i,2) == ' ' ) then
      string1 = 'model_nml:gitm_state_variables not fully specified'
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Normally sure variable exists in GITM restart variable list
   ! This is done in dart_gitm_mod:decode_gitm_indices() 


   ! Make sure DART kind is valid

   if( get_index_for_quantity(dartstr) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Record the contents of the DART state vector

   if ( debug > 0 ) then
      write(logfileunit,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2))
      write(     *     ,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2))
   endif

   ngood = ngood + 1
enddo MyLoop

if (ngood == nrows) then
   string1 = 'WARNING: There is a possibility you need to increase ''max_state_variables'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'verify_state_variables',string1,source,revision,revdate,text2=string2)
endif

end subroutine verify_state_variables


!==================================================================


function find_index(x, xa, n)
!------------------------------------------------------------------
! This function returns the array index (here, the value returned by
! find_index is designated as i) such that x is between xa(i) and xa(i+1).
! If x is less than xa(1), then i=-1 is returned.  If x is greater than
! xa(n), then i=-1 is returned.  It is assumed that the values of
! xa increase monotonically with increasing i.

    integer              :: find_index
    integer,  intent(in) :: n             ! array size
    real(r8), intent(in) :: xa(n)         ! array of locations
    real(r8), intent(in) :: x             ! location of interest

    integer :: lower, upper, mid          ! lower and upper limits, and midpoint
    integer :: order

    lower = 0
    upper = n + 1

    IF( xa(1) .lt. xa(n) ) THEN
      order = 1
    ELSE
      order = -1
    ENDIF

    IF ( x .gt. maxval(xa) ) THEN
      mid = -1
    ELSEIF ( x .lt. minval(xa) ) THEN
      mid = -1
    ELSE

10    IF ((upper-lower).gt.1) THEN
        mid=(lower+upper)/2
        IF( order .eq. 1 ) THEN
          IF (x .ge. xa(mid)) THEN
            lower = mid
          ELSE
            upper = mid
          ENDIF
          go to 10
        ELSE
          IF (x .lt. xa(mid)) THEN
            lower = mid
          ELSE
            upper = mid
          ENDIF
          go to 10
        ENDIF
      ENDIF

    ENDIF

    find_index = lower

return
end function find_index


!==================================================================


subroutine define_var_dims(myprogvar, ndims, DimIDs, memberDimID, unlimitedDimID, &
               LONDimID, LATDimID, ALTDimID, WLDimID)
!------------------------------------------------------------------

type(progvartype),     intent(in)  :: myprogvar
integer,               intent(out) :: ndims
integer, dimension(:), intent(out) :: DimIDs
integer,               intent(in)  :: memberDimID, unlimitedDimID
integer,               intent(in)  :: LONDimID, LATDimID, ALTDimID
integer,               intent(in)  :: WLDimID

select case( myprogvar%storder )
case('xyz3d')

      ndims = 5

      DimIDs(1) = LONDimID
      DimIDs(2) = LATDimID
      DimIDs(3) = ALTDimID
      DimIDs(4) = memberDimID
      DimIDs(5) = unlimitedDimID

case('xy2d')

      ndims = 4

      DimIDs(1) = LONDimID
      DimIDs(2) = LATDimID
      DimIDs(3) = memberDimID
      DimIDs(4) = unlimitedDimID

case('0d')

      ndims = 3

      DimIDs(1) = WLDimID
      DimIDs(2) = memberDimID
      DimIDs(3) = unlimitedDimID

case('x1d')

      ndims = 3

      DimIDs(1) = LONDimID
      DimIDs(2) = memberDimID
      DimIDs(3) = unlimitedDimID

case('y1d')

      ndims = 3

      DimIDs(1) = LATDimID
      DimIDs(2) = memberDimID
      DimIDs(3) = unlimitedDimID

case('z1d')

      ndims = 3

      DimIDs(1) = ALTDimID
      DimIDs(2) = memberDimID
      DimIDs(3) = unlimitedDimID

case default

      write(string1,*)'unknown storage order '//trim(myprogvar%storder)//&
                              ' for variable '//trim(myprogvar%varname)
      call error_handler(E_ERR,'define_var_dims',string1,source,revision,revdate)

end select

return
end subroutine define_var_dims


!==================================================================


function read_in_real(iunit,varname,filename)
!------------------------------------------------------------------
integer,          intent(in) :: iunit
character(len=*), intent(in) :: varname,filename
real(r8)                     :: read_in_real

character(len=100) :: cLine
integer :: i, ios

! Read a line and remove anything after a space or TAB
read(iunit,'(a)',iostat=ios) cLine
if (ios /= 0) then
   write(string1,*) 'cannot find '//trim(varname)//' in '//trim(filename)
   call error_handler(E_ERR,'get_grid_dims',string1,source,revision,revdate)
endif

i=index(cLine,' ');     if( i > 0 ) cLine(i:len(cLine))=' '
i=index(cLine,char(9)); if( i > 0 ) cLine(i:len(cLine))=' '

! Now that we have a line with nothing else ... parse it
read(cLine,*,iostat=ios)read_in_real

if(ios /= 0) then
   write(string1,*)'unable to read '//trim(varname)//' in '//trim(filename)
   call error_handler(E_ERR,'read_in_real',string1,source,revision,revdate)
endif

end function read_in_real


!==================================================================


function read_in_int(iunit,varname,filename)
!------------------------------------------------------------------
integer,          intent(in) :: iunit
character(len=*), intent(in) :: varname,filename
integer                      :: read_in_int

character(len=100) :: cLine
integer :: i, ios

! Read a line and remove anything after a space or TAB
read(iunit,'(a)',iostat=ios) cLine
if (ios /= 0) then
   write(string1,*) 'cannot find '//trim(varname)//' in '//trim(filename)
   call error_handler(E_ERR,'get_grid_dims',string1,source,revision,revdate)
endif

i=index(cLine,' ');     if( i > 0 ) cLine(i:len(cLine))=' '
i=index(cLine,char(9)); if( i > 0 ) cLine(i:len(cLine))=' '

read(cLine,*,iostat=ios)read_in_int

if(ios /= 0) then
   write(string1,*)'unable to read '//trim(varname)//' in '//trim(filename)
   call error_handler(E_ERR,'read_in_int',string1,source,revision,revdate,&
             text2=cLine)
endif

end function read_in_int


!==================================================================


subroutine sortindexlist(list, x, inum)
!------------------------------------------------------------------
! sort list x into order based on values in list.
! should only be called on short ( < hundreds) of values or will be slow

integer, intent(inout) :: list(:)
integer, intent(inout) :: x(:)
integer, intent(in)    :: inum

integer :: tmp
integer :: j, k

!  DO A N^2 SORT - only use for short lists
do j = 1, inum - 1
   do k = j + 1, inum
      ! if list() is in wrong order, exchange both list items and
      ! items in x array.
      if(list(j) .gt. list(k)) then
         tmp = list(k)
         list(k) = list(j)
         list(j) = tmp
         tmp = x(k)
         x(k) = x(j)
         x(j) = tmp
      end if
   end do
end do
end subroutine sortindexlist


!===================================================================
end module model_mod
!===================================================================

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
