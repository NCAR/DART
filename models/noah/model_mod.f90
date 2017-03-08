! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

use        types_mod, only : r8, MISSING_R8, obstypelength

use time_manager_mod, only : time_type, set_time, set_date, get_time,          &
                             print_time, print_date, set_calendar_type,        &
                             operator(*),  operator(+), operator(-),           &
                             operator(>),  operator(<), operator(/),           &
                             operator(/=), operator(<=)

use     location_mod, only : location_type, get_dist, query_location,          &
                             get_close_maxdist_init, get_close_type,           &
                             set_location, get_location, horiz_dist_only,      &
                             vert_is_surface,  VERTISSURFACE,                  &
                             vert_is_height,   VERTISHEIGHT,                   &
                             vert_is_level,    VERTISLEVEL,                    &
                             get_close_obs_init, get_close_obs,                &
                             set_location_missing, write_location

use    utilities_mod, only : register_module, error_handler, nc_check,         &
                             E_ERR, E_MSG, logfileunit, get_unit,              &
                             nmlfileunit, do_output, do_nml_file, do_nml_term, &
                             find_namelist_in_file, check_namelist_read,       &
                             file_exist, find_textfile_dims, file_to_text

use     obs_kind_mod, only : QTY_SOIL_TEMPERATURE,    &
                             QTY_LIQUID_WATER,        &
                             QTY_ICE,                 &
                             QTY_SNOWCOVER_FRAC,      &
                             QTY_SNOW_THICKNESS,      &
                             QTY_LEAF_CARBON,         &
                             QTY_WATER_TABLE_DEPTH,   &
                             QTY_GEOPOTENTIAL_HEIGHT, &
                             get_index_for_quantity

use mpi_utilities_mod, only: my_task_id
use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

use typesizes
use netcdf

implicit none
private

! required by DART code - will be called from filter and other
! DART executables.  interfaces to these routines are fixed and
! cannot be changed in any way.
public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          get_model_time_step,    &
          end_model,              &
          static_init_model,      &
          init_time,              &
          init_conditions,        &
          nc_write_model_atts,    &
          nc_write_model_vars,    &
          pert_model_state,       &
          get_close_maxdist_init, &
          get_close_obs_init,     &
          get_close_obs,          &
          ens_mean_for_model

! not required by DART but for larger models can be useful for
! utility programs that are tightly tied to the other parts of
! the model_mod code.
public :: noah_to_dart_vector, &
          dart_vector_to_model_file, &
          get_noah_restart_filename, &
          get_noah_timestepping,     &
          get_debug_level

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------
! The NSOLDX parameter comes from the NOAH source code. We need it
! because we have to read the NOAH namelist for timestep information.
!------------------------------------------------------------------

integer :: nfields
integer, parameter :: NSOLDX = 100
integer, parameter :: MAX_STATE_VARIABLES = 40
integer, parameter :: NUM_STATE_TABLE_COLUMNS = 2

!------------------------------------------------------------------
! things which can/should be in the DART model_nml
! The variables in the noah restart file that are used to create the
! DART state vector are specified in the input.nml:model_nml namelist.
! For example:
!
! noah_state_variables  = 'STC',    'QTY_SOIL_TEMPERATURE',
!                         'SMC',    'QTY_SOIL_MOISTURE',
!                         'SH2O',   'QTY_LIQUID_SOIL_MOISTURE',
!                         'T1',     'QTY_SKIN_TEMPERATURE',
!                         'SNOWH',  'QTY_SNOW_DEPTH',
!                         'SNEQV',  'QTY_LIQUID_EQUIVALENT',
!                         'CMC',    'QTY_CANOPY_WATER'
!
!------------------------------------------------------------------

character(len=128)    :: noah_netcdf_filename   = 'restart.nc'
integer               :: assimilation_period_days     = 0
integer               :: assimilation_period_seconds  = 60
real(r8)              :: model_perturbation_amplitude = 0.2
logical               :: output_state_vector          = .true.
integer               :: debug    = 0  ! turn up for more and more debug messages
character(len=obstypelength) :: noah_state_variables(NUM_STATE_TABLE_COLUMNS,MAX_STATE_VARIABLES) = ' '

namelist /model_nml/ noah_netcdf_filename, &
          assimilation_period_days, assimilation_period_seconds,   &
          model_perturbation_amplitude, output_state_vector,       &
          debug, noah_state_variables

!------------------------------------------------------------------
! Everything needed to recreate the NOAH METADTA_NAMELIST
! We are going to require that the namelist be in a file whose name
! is (I believe) standard .... namelist.hrldas
!
! To restart the file, we write a new namelist.
! DART needs to write a NOAH-compatible namelist.
!------------------------------------------------------------------

! ZSOIL, set through the namelist, is the BOTTOM of each soil layer (m)
! Values are negative, implying depth below the surface.

character(len=128)    :: noah_namelist_filename = 'namelist.hrldas' ! mandate

real(r8), dimension(NSOLDX) :: zsoil
integer                     :: nsoil

CHARACTER(len=256) :: indir = "."
character(len=256) :: outdir = "."
character(len=256) :: hrldas_constants_file = " "
character(len=256) :: external_fpar_filename_template = " "
character(len=256) :: external_lai_filename_template = " "
character(len=256) :: restart_filename_requested = " "
integer            :: split_output_count = 1
integer            :: restart_frequency_hours = -999
integer            :: output_timestep  = -999
integer            :: subwindow_xstart = 1
integer            :: subwindow_ystart = 1
integer            :: subwindow_xend = 0
integer            :: subwindow_yend = 0
integer            :: sfcdif_option = 0
integer            :: iz0tlnd = 0
logical            :: update_snow_from_forcing = .TRUE.

integer  :: start_year, start_month, start_day, start_hour, start_min
integer  :: noah_timestep = -999
integer  :: forcing_timestep = -999
integer  :: khour = 0
integer  :: kday  = 0
real(r8) :: zlvl, zlvl_wind

namelist / NOAHLSM_OFFLINE/ indir, nsoil, zsoil, forcing_timestep, noah_timestep, &
       start_year, start_month, start_day, start_hour, start_min, &
       restart_frequency_hours, output_timestep, &
       split_output_count, sfcdif_option, iz0tlnd, update_snow_from_forcing, &
       khour, kday, zlvl, zlvl_wind, hrldas_constants_file, outdir, restart_filename_requested, &
       external_fpar_filename_template, external_lai_filename_template, &
       subwindow_xstart, subwindow_xend, subwindow_ystart, subwindow_yend

!------------------------------------------------------------------

! define model parameters here
type(time_type)     :: time_step
type(location_type),allocatable, dimension(:) :: state_loc

! Everything needed to describe a variable

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=NF90_MAX_NAME) :: MemoryOrder
   character(len=NF90_MAX_NAME) :: description
   character(len=NF90_MAX_NAME) :: stagger
   character(len=obstypelength), dimension(NF90_MAX_VAR_DIMS) :: dimnames
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
   integer  :: numdims
   integer  :: maxlevels
   integer  :: xtype
   integer  :: varsize     ! prod(dimlens(1:numdims))
   integer  :: index1      ! location in dart state vector of first occurrence
   integer  :: indexN      ! location in dart state vector of last  occurrence
   integer  :: dart_kind
   integer  :: rangeRestricted
   real(r8) :: maxvalue
   real(r8) :: minvalue
   character(len=obstypelength) :: kind_string
end type progvartype

type(progvartype), dimension(MAX_STATE_VARIABLES) :: progvar

!------------------------------------------------------------------------------
! These are the metadata arrays that are the same size as the state vector.

real(r8), allocatable, dimension(:) :: ens_mean ! may be needed for forward ops
real(r8), allocatable, dimension(:) :: levels   ! depth

!------------------------------------------------------------------
! module storage
!------------------------------------------------------------------

integer            :: model_size       ! the state vector length
type(time_type)    :: model_time       ! valid time of the model state
type(time_type)    :: model_time_step  ! smallest time to adv model
character(len=256) :: string1, string2, string3
logical, save      :: module_initialized = .false.
character(len=32)  :: calendar = 'Gregorian'

real(r8), allocatable, dimension(:,:) :: xlong, xlat
integer :: south_north, west_east

INTERFACE vector_to_prog_var
      MODULE PROCEDURE vector_to_1d_prog_var
      MODULE PROCEDURE vector_to_2d_prog_var
      MODULE PROCEDURE vector_to_3d_prog_var
END INTERFACE

!==================================================================
contains
!==================================================================



subroutine static_init_model()
!------------------------------------------------------------------
! one time initialization of the model

! Local variables - all the important ones have module scope

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname
character(len=NF90_MAX_NAME)          :: dimname
character(len=obstypelength)          :: kind_string

integer  :: VarID, dimlen, varsize
integer  :: iunit, io, ivar
integer  :: i, index1, nLayers

if ( module_initialized ) return ! only need to do this once.

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.
module_initialized = .true.

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the DART namelist
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the DART namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Check to make sure the required NOAH input files exist
if ( .not. file_exist(noah_netcdf_filename) ) then
   write(string1,*) 'cannot open NOAH restart file ', trim(noah_netcdf_filename),' for reading.'
   call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
endif
if ( .not. file_exist(noah_namelist_filename) ) then
   write(string1,*) 'cannot open NOAH namelist file ', trim(noah_namelist_filename),' for reading.'
   call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
endif

! Read the NOAH namelist
call find_namelist_in_file(trim(noah_namelist_filename), 'NOAHLSM_OFFLINE', iunit)
read(iunit, nml = NOAHLSM_OFFLINE, iostat = io)
call check_namelist_read(iunit, io, 'NOAHLSM_OFFLINE')

! Record the NOAH namelist
if (do_nml_file()) write(nmlfileunit, nml=NOAHLSM_OFFLINE)
if (do_nml_term()) write(     *     , nml=NOAHLSM_OFFLINE)

! Check to make sure the NOAH constants file exists
if ( .not. file_exist(hrldas_constants_file) ) then
   write(string1,*) 'cannot open file ', trim(hrldas_constants_file),' for reading.'
   call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
endif

! Check to make sure the required NOAH namelist items are set:
if ( (kday             < 0    ) .or. &
     (khour            < 0    ) .or. &
     (forcing_timestep /= 3600) .or. &
     (noah_timestep    /= 3600) .or. &
     (output_timestep  /= 3600) .or. &
     (restart_frequency_hours /= 1) ) then
   write(string3,*)'the only configuration supported is for hourly timesteps'
   write(string2,*)'restart_frequency_hours must be equal to the noah_timstep'
   write(string1,*)'unsupported noah namelist settings'
   call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate,&
                            text2=string2,text3=string3)
endif

call get_hrldas_constants(hrldas_constants_file)

! The time_step in terms of a time type must also be initialized.

call set_calendar_type( calendar )

call nc_check(nf90_open(adjustl(noah_netcdf_filename), NF90_NOWRITE, iunit), &
                   'static_init_model', 'open '//trim(noah_netcdf_filename))

model_time      = get_state_time(iunit, trim(noah_netcdf_filename))

! FIXME ... make sure model_time_step is attainable given OUTPUT_TIMESTEP
model_time_step = set_time(assimilation_period_seconds, assimilation_period_days)

if (do_output() .and. (debug > 0)) then
   call print_date(model_time     ,' static_init_model:model date')
   call print_time(model_time     ,' static_init_model:model time')
   call print_time(model_time_step,' static_init_model:model timestep')
endif

! Make sure the number of soil layers is as we expect

call nc_check(nf90_inq_dimid(iunit, 'soil_layers_stag', dimIDs(1)), &
                  'static_init_model','inq_dimid soil_layers_stag '//trim(noah_netcdf_filename))
call nc_check(nf90_inquire_dimension(iunit, dimIDs(1), len=nLayers), &
                  'static_init_model','inquire_dimension soil_layers_stag '//trim(noah_netcdf_filename))

if (nsoil /= nLayers) then
   write(string1,*) 'Expected ',nsoil,' soil layers ', &
                       trim(noah_netcdf_filename),' has ',nLayers
   call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
endif

! FIXME ... extend to 2D case ... should be in get_state_meta_data()
! Convert soil thicknesses (from namelist.hrldas) to "heights".
! Closer to the center of the earth is an increasingly large negative number
allocate(state_loc(0:nsoil))
state_loc(0)   = set_location(xlong(1,1), xlat(1,1), 0.0_r8, VERTISHEIGHT)
do i = 1,nsoil
   state_loc(  i) = set_location(xlong(1,1), xlat(1,1), zsoil(i), VERTISHEIGHT)
enddo

!---------------------------------------------------------------
! Compile the list of NOAH variables to use in the creation
! of the DART state vector. Required to determine model_size.
!
! Verify all variables are in the NOAH netcdf file.
! Compute the offsets into the state vector for each variable type.
! Record the extent of the variable type in the state vector.

call verify_state_variables(iunit, noah_netcdf_filename, nfields)

index1  = 1
FILL_PROGVAR : do ivar = 1, nfields

   varname                   = trim(noah_state_variables(1,ivar))
   kind_string               = trim(noah_state_variables(2,ivar))
   progvar(ivar)%varname     = varname
   progvar(ivar)%kind_string = kind_string
   progvar(ivar)%dart_kind   = get_index_for_quantity( progvar(ivar)%kind_string )
   progvar(ivar)%dimlens     = 0
   progvar(ivar)%dimnames    = ' '
   progvar(ivar)%maxlevels   = 0

   string2 = trim(noah_netcdf_filename)//' '//trim(varname)

   call nc_check(nf90_inq_varid(iunit, trim(varname), VarID), &
            'static_init_model', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(iunit, VarID, dimids=dimIDs, &
                 ndims=progvar(ivar)%numdims, xtype=progvar(ivar)%xtype), &
            'static_init_model', 'inquire '//trim(string2))

   ! If the long_name and/or units attributes are set, get them.
   ! They are not REQUIRED to exist but are nice to use if they are present.

   if( nf90_inquire_attribute(    iunit, VarID, 'long_name') == NF90_NOERR ) then
      call nc_check( nf90_get_att(iunit, VarID, 'long_name' , progvar(ivar)%long_name), &
                  'static_init_model', 'get_att long_name '//trim(string2))
   else
      progvar(ivar)%long_name = varname
   endif

   if( nf90_inquire_attribute(    iunit, VarID, 'units') == NF90_NOERR )  then
      call nc_check( nf90_get_att(iunit, VarID, 'units' , progvar(ivar)%units), &
                  'static_init_model', 'get_att units '//trim(string2))
   else
      progvar(ivar)%units = '-'
   endif

   if( nf90_inquire_attribute(    iunit, VarID, 'MemoryOrder') == NF90_NOERR )  then
      call nc_check( nf90_get_att(iunit, VarID, 'MemoryOrder' , progvar(ivar)%MemoryOrder), &
                  'static_init_model', 'get_att MemoryOrder '//trim(string2))
   else
      progvar(ivar)%MemoryOrder = '-'
   endif

   if( nf90_inquire_attribute(    iunit, VarID, 'stagger') == NF90_NOERR )  then
      call nc_check( nf90_get_att(iunit, VarID, 'stagger' , progvar(ivar)%stagger), &
                  'static_init_model', 'get_att stagger '//trim(string2))
   else
      progvar(ivar)%stagger = '-'
   endif

   ! This is fundamentally hardcoded clamping. See WRF/model_mod.f90 for a namelist-driven
   ! example.  It would be nice to use the netCDF file valid_range attribute ...
   !
   ! if the variable is bounded, then we need to know how to restrict it.
   ! rangeRestricted == 0 is unbounded
   ! rangeRestricted == 1 is bounded below
   ! rangeRestricted == 2 is bounded above           ( TJH unsupported )
   ! rangeRestricted == 3 is bounded above and below ( TJH unsupported )

   if ( varname == 'QFX' ) then
      progvar(ivar)%rangeRestricted = 0
      progvar(ivar)%minvalue        = -1.0_r8*huge(0.0_r8)
      progvar(ivar)%maxvalue        = huge(0.0_r8)
   elseif ( varname == 'HFX' ) then
      progvar(ivar)%rangeRestricted = 0
      progvar(ivar)%minvalue        = -1.0_r8*huge(0.0_r8)
      progvar(ivar)%maxvalue        = huge(0.0_r8)
   else
      progvar(ivar)%rangeRestricted = 1
      progvar(ivar)%minvalue        = 0.0_r8
      progvar(ivar)%maxvalue        = huge(0.0_r8)
   endif

   ! These variables have a Time dimension. We only want the most recent time.

   varsize = 1
   dimlen  = 1
   DimensionLoop : do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(iunit, dimIDs(i), name=dimname, len=dimlen), &
                                          'static_init_model', string1)

      if ((trim(dimname) == 'Time') .or. (trim(dimname) == 'time')) dimlen = 1
      progvar(ivar)%dimlens( i) = dimlen
      progvar(ivar)%dimnames(i) = trim(dimname)
      varsize = varsize * dimlen

   enddo DimensionLoop

   progvar(ivar)%varsize     = varsize
   progvar(ivar)%index1      = index1
   progvar(ivar)%indexN      = index1 + varsize - 1
   index1                    = index1 + varsize      ! sets up for next variable

   if ( (do_output()) .and. debug > 10 ) then
      write(logfileunit,*)
      write(logfileunit,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(logfileunit,*) '  long_name   ',trim(progvar(ivar)%long_name)
      write(logfileunit,*) '  MemoryOrder ',trim(progvar(ivar)%MemoryOrder)
      write(logfileunit,*) '  stagger     ',trim(progvar(ivar)%stagger)
      write(logfileunit,*) '  units       ',trim(progvar(ivar)%units)
      write(logfileunit,*) '  xtype       ',progvar(ivar)%xtype
      write(logfileunit,*) '  dimnames    ',progvar(ivar)%dimnames(1:progvar(ivar)%numdims)
      write(logfileunit,*) '  dimlens     ',progvar(ivar)%dimlens( 1:progvar(ivar)%numdims)
      write(logfileunit,*) '  numdims     ',progvar(ivar)%numdims
      write(logfileunit,*) '  varsize     ',progvar(ivar)%varsize
      write(logfileunit,*) '  index1      ',progvar(ivar)%index1
      write(logfileunit,*) '  indexN      ',progvar(ivar)%indexN
      write(logfileunit,*) '  dart_kind   ',progvar(ivar)%dart_kind
      write(logfileunit,*) '  kind_string ',progvar(ivar)%kind_string

      write(     *     ,*)
      write(     *     ,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(     *     ,*) '  long_name   ',trim(progvar(ivar)%long_name)
      write(     *     ,*) '  MemoryOrder ',trim(progvar(ivar)%MemoryOrder)
      write(     *     ,*) '  stagger     ',trim(progvar(ivar)%stagger)
      write(     *     ,*) '  units       ',trim(progvar(ivar)%units)
      write(     *     ,*) '  xtype       ',progvar(ivar)%xtype
      write(     *     ,*) '  dimnames    ',progvar(ivar)%dimnames(1:progvar(ivar)%numdims)
      write(     *     ,*) '  dimlens     ',progvar(ivar)%dimlens( 1:progvar(ivar)%numdims)
      write(     *     ,*) '  numdims     ',progvar(ivar)%numdims
      write(     *     ,*) '  varsize     ',progvar(ivar)%varsize
      write(     *     ,*) '  index1      ',progvar(ivar)%index1
      write(     *     ,*) '  indexN      ',progvar(ivar)%indexN
      write(     *     ,*) '  dart_kind   ',progvar(ivar)%dart_kind
      write(     *     ,*) '  kind_string ',progvar(ivar)%kind_string
   endif

enddo FILL_PROGVAR

call nc_check(nf90_close(iunit), 'static_init_model', 'close '//trim(noah_netcdf_filename))

model_size = progvar(nfields)%indexN

if ( (do_output()) .and. debug > 99 ) then

   write(*,*)
   do i=0,nsoil
      call write_location(iunit,state_loc(i),charstring=string1)
      write(*,*)'location ',i,' is ',trim(string1)
   enddo

endif

end subroutine static_init_model



subroutine init_conditions(x)
!------------------------------------------------------------------
!
! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no
! synthetic data experiments using perfect_model_obs are planned,
! this can be a NULL INTERFACE.

real(r8), intent(out) :: x(:)

if ( .not. module_initialized ) call static_init_model

write(string1,*) 'PROBLEM: no known way to set arbitrary initial conditions.'
write(string2,*) 'start_from_restart must be .true. in this model.'
call error_handler(E_ERR,'init_conditions',string1,source,revision,revdate, &
                               text2=string2)

x = MISSING_R8

end subroutine init_conditions



subroutine adv_1step(x, time)
!------------------------------------------------------------------
!
! Does a single timestep advance of the model. The input value of
! the vector x is the starting condition and x is updated to reflect
! the changed state after a timestep. The time argument is intent
! in and is used for models that need to know the date/time to
! compute a timestep, for instance for radiation computations.
! This interface is only called if the namelist parameter
! async is set to 0 in perfect_model_obs of filter or if the
! program integrate_model is to be used to advance the model
! state as a separate executable. If one of these options
! is not going to be used (the model will only be advanced as
! a separate model-specific executable), this can be a
! NULL INTERFACE.

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

if ( .not. module_initialized ) call static_init_model

write(string1,*) 'PROBLEM: cannot advance model with async == 0.'
write(string2,*) 'async == 2 is a good choice.'
call error_handler(E_ERR,'adv_1step',string1,source,revision,revdate, &
                               text2=string2)

end subroutine adv_1step



function get_model_size()
!------------------------------------------------------------------
!
! Returns the size of the model as an integer. Required for all
! applications.

integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size



subroutine init_time(time)
!------------------------------------------------------------------
!
! Companion interface to init_conditions. Returns a time that is somehow
! appropriate for starting up a long integration of the model.
! At present, this is only used if the namelist parameter
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no
! synthetic data experiments using perfect_model_obs are planned,
! this can be a NULL INTERFACE.

type(time_type), intent(out) :: time

if ( .not. module_initialized ) call static_init_model

time = model_time

end subroutine init_time



subroutine model_interpolate(x, location, itype, obs_val, istatus)
!------------------------------------------------------------------
!
! Given a state vector, a location, and a model state variable type,
! interpolates the state variable field to that location and returns
! the value in obs_val. The istatus variable should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case an alternate value should be returned. The itype variable
! is a model specific integer that specifies the type of field (for
! instance temperature, zonal wind component, etc.).

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: itype
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

real(r8)               :: loc_lon, loc_lat, loc_depth
real(r8), dimension(3) :: loc
integer,  dimension(1) :: loninds, latinds
integer                :: gridloni, gridlatj, zlev, n, ivar, indx

if ( .not. module_initialized ) call static_init_model

! FIXME - for the single column case - there is no obvious way to
! determine the extent of the domain ... EVERYTHING matches.

if( west_east*south_north /= 1 ) then
     write(string1,*) 'PROBLEM: not set up for a case with multiple locations yet.'
     call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate)
endif

! The return code for successful return should be 0.
! Any positive number is an error.
! Negative values are reserved for use by the DART framework.
! Using distinct positive values for different types of errors can be
! useful in diagnosing problems.

istatus = 1
obs_val = MISSING_R8

! if observation is outside region encompassed in the history file - fail

loc       = get_location(location) ! loc is in DEGREES
loc_lon   = loc(1)
loc_lat   = loc(2)
loc_depth = loc(3)
! latinds   = minloc(abs(XLAT  - loc_lat))   ! these return 'arrays' ...
! loninds   = minloc(abs(XLONG - loc_lon))   ! these return 'arrays' ...
! gridlatj  = latinds(1)
! gridloni  = loninds(1)

! one use of model_interpolate is to allow other modules/routines
! the ability to 'see' the model levels. To do this, we can create
! locations with model levels and 'interpolate' them to
! QTY_GEOPOTENTIAL_HEIGHT

if ( (itype == QTY_GEOPOTENTIAL_HEIGHT) .and. vert_is_level(location) ) then
   if (nint(loc_depth) > nsoil) then
      obs_val = MISSING_R8
      istatus = 1
   else
      obs_val = zsoil(nint(loc_depth))
      istatus = 0
   endif
   return
endif

! need to know what variable we are interpolating.

ivar = 0
FindVariable : do n = 1,nfields
   if( progvar(n)%dart_kind == itype ) then
      ivar = n
      exit FindVariable
   endif
enddo FindVariable

if (ivar == 0) then
   write(string1,*) 'unable to find state vector component matching type ',itype
   call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate)
endif

if ( progvar(ivar)%varsize == 1 ) then
   indx = progvar(ivar)%index1
else
   ! This assumes the soil layer has a constant value.
   DEPTH : do n = 1,nsoil
      if (loc_depth >= zsoil(n)) then
         zlev = n
         exit DEPTH
      endif
   enddo DEPTH
   ! FIXME what if zlev is never set
   indx = progvar(ivar)%index1 + zlev - 1
endif

obs_val = x( indx )
if (obs_val /= MISSING_R8) istatus = 0

if ( (do_output()) .and. debug > 20 ) then
   write(*,*)'model_interpolate : progvar%kind_string is ',trim(progvar(ivar)%kind_string)
   write(*,*)'model_interpolate : state index         is ',indx
   write(*,*)'model_interpolate : value               is ',obs_val
   call write_location(n,location,charstring=string1)
   write(*,*)'observation location ',trim(string1)
endif

end subroutine model_interpolate



function get_model_time_step()
!------------------------------------------------------------------
!
! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.

type(time_type) :: get_model_time_step

if ( .not. module_initialized ) call static_init_model

! The NOAH model can only be advanced in multiples of the restart frequency.

get_model_time_step = set_time(khour*3600,kday)

end function get_model_time_step



subroutine get_state_meta_data(index_in, location, var_type)
!------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated location. A second intent(out) optional argument kind
! can be returned if the model has more than one type of field (for
! instance temperature and zonal wind component). This interface is
! required for all filter applications as it is required for computing
! the distance between observations and state variables.

integer,             intent(in)            :: index_in
type(location_type), intent(out)           :: location
integer,             intent(out), optional :: var_type

integer :: n, layer, ivar

if ( .not. module_initialized ) call static_init_model

if( west_east*south_north /= 1 ) then
     write(string1,*) 'PROBLEM: not set up for a case with multiple locations yet.'
     call error_handler(E_ERR,'get_state_meta_data',string1,source,revision,revdate)
endif

layer = -1

FindIndex : do n = 1,nfields
   if( (progvar(n)%index1 <= index_in) .and. (index_in <= progvar(n)%indexN) ) then
      layer    = index_in - progvar(n)%index1 + 1
      var_type = progvar(n)%dart_kind
      ivar     = n
      exit FindIndex
   endif
enddo FindIndex

if ( (do_output()) .and. debug > 30 ) then
   write(*,*)'get_state_meta_data: index_in is ',index_in
   write(*,*)'get_state_meta_data: ivar     is ',ivar
   write(*,*)'get_state_meta_data: layer    is ',layer
   write(*,*)
endif

if( layer == -1 ) then
     write(string1,*) 'Problem, cannot find base_offset, index_in is: ', index_in
     call error_handler(E_ERR,'get_state_meta_data',string1,source,revision,revdate)
endif

if (progvar(ivar)%varsize == 1) layer = 0

location = state_loc(layer)

end subroutine get_state_meta_data



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

! good style ... perhaps you could deallocate stuff (from static_init_model?).
! deallocate(state_loc)
if ( .not. module_initialized ) call static_init_model

end subroutine end_model



function nc_write_model_atts( ncFileID ) result (ierr)
!------------------------------------------------------------------
! This routine writes all the netCDF 'infrastructure' and sets up the
! global attributes, dimensions, coordinate variables, and output variables.
! The actuall filling of the output variables is done by
! nc_write_model_vars() which can be called repeatedly for each
! assimilation cycle.
!
! All errors are fatal, so the return code is always '0 == normal'
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

!----------------------------------------------------------------------
! variables if we just blast out one long state vector
!----------------------------------------------------------------------

integer :: StateVarDimID    ! netCDF pointer to state variable dimension (model size)
integer :: MemberDimID      ! netCDF pointer to dimension of ensemble    (ens_size)
integer :: TimeDimID        ! netCDF pointer to time dimension           (unlimited)

integer :: StateVarVarID   ! netCDF pointer to state variable coordinate array
integer :: StateVarID      ! netCDF pointer to 3D [state,copy,time] array

!----------------------------------------------------------------------
! variables if we parse the state vector into prognostic variables
!----------------------------------------------------------------------

integer :: weDimID
integer :: snDimID
integer :: nsoilDimID
integer :: myndims
integer :: ivar, varID
integer, dimension(NF90_MAX_VAR_DIMS) :: mydimids
character(len=NF90_MAX_NAME) :: varname

!----------------------------------------------------------------------
! variables for the namelist output
!----------------------------------------------------------------------

character(len=129), allocatable, dimension(:) :: textblock
integer :: LineLenDimID, nlinesDimID, nmlVarID
integer :: nlines, linelen
logical :: has_noah_namelist

!----------------------------------------------------------------------
! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.
!----------------------------------------------------------------------

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

integer :: i

character(len=128) :: filename

if ( .not. module_initialized ) call static_init_model

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file,
! and then put into define mode.
!-------------------------------------------------------------------------------

ierr = -1 ! assume things go poorly

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename,*) 'ncFileID', ncFileID

call nc_check(nf90_inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID), &
                     'nc_write_model_atts', 'inquire '//trim(filename))
call nc_check(nf90_redef(ncFileID), 'nc_write_model_atts', 'redef '//trim(filename))

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension.
! Our job is create the 'model size' dimension.
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncid=ncFileID, name='NMLlinelen', dimid = linelenDimID), &
                  'nc_write_model_atts', 'inq_dimid NMLlinelen '//trim(filename))
call nc_check(nf90_inq_dimid(ncid=ncFileID, name='copy', dimid=MemberDimID), &
                  'nc_write_model_atts', 'inq_dimid copy '//trim(filename))
call nc_check(nf90_inq_dimid(ncid=ncFileID, name='time', dimid= TimeDimID), &
                  'nc_write_model_atts', 'inq_dimid time '//trim(filename))

if ( TimeDimID /= unlimitedDimId ) then
   write(string1,*)'Time Dimension ID ',TimeDimID, &
                     ' should equal Unlimited Dimension ID',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts', string1, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size / state variable dimension / whatever ...
!-------------------------------------------------------------------------------

call nc_check(nf90_def_dim(ncid=ncFileID, name='StateVariable',  &
                           len=model_size, dimid=StateVarDimID), &
                           'nc_write_model_atts', 'def_dim state '//trim(filename))

!-------------------------------------------------------------------------------
! Write Global Attributes
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'creation_date' ,str1), &
                          'nc_write_model_atts', 'put_att creation_date '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_source'  ,source), &
                          'nc_write_model_atts', 'put_att model_source '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revision',revision), &
                          'nc_write_model_atts', 'put_att model_revision '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revdate' ,revdate), &
                          'nc_write_model_atts', 'put_att model_revdate '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model','NOAH'), &
                          'nc_write_model_atts', 'put_att model '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'calendar',trim(calendar)), &
                          'nc_write_model_atts', 'put_att calendar '//trim(filename))

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'HRLDAS_CONSTANTS_FILE',trim(hrldas_constants_file)), &
                          'nc_write_model_atts', 'put_att HRLDAS_CONSTANTS_FILE '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'HRLDAS_INDIR',trim(INDIR)), &
                          'nc_write_model_atts', 'put_att HRLDAS_INDIR '//trim(filename))

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'RESTART_FILENAME_REQUESTED ',trim(restart_filename_requested)), &
                          'nc_write_model_atts', 'put_att RESTART_FILENAME_REQUESTED '//trim(filename))

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'KDAY',KDAY), &
                          'nc_write_model_atts', 'put_att KDAY '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'KHOUR',KHOUR), &
                          'nc_write_model_atts', 'put_att KHOUR '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'FORCING_TIMESTEP',forcing_timestep), &
                          'nc_write_model_atts', 'put_att FORCING_TIMESTEP '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'NOAH_TIMESTEP',noah_timestep), &
                          'nc_write_model_atts', 'put_att NOAH_TIMESTEP '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'OUTPUT_TIMESTEP',output_timestep), &
                          'nc_write_model_atts', 'put_att OUTPUT_TIMESTEP '//trim(filename))

! call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'Soil_type_index',Soil_type_index), &
!                           'nc_write_model_atts', 'put_att Soil_type_index'//trim(filename))
! call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'Vegetation_type_index',Vegetation_type_index), &
!                           'nc_write_model_atts', 'put_att Vegetation_type_index')
! call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'Urban_veg_category',Urban_veg_category), &
!                           'nc_write_model_atts', 'put_att Urban_veg_category'//trim(filename))

!-------------------------------------------------------------------------------
! Determine shape of namelist.
! long lines are truncated when read into textblock
!-------------------------------------------------------------------------------

call find_textfile_dims(noah_namelist_filename, nlines, linelen)
if (nlines > 0) then
  has_noah_namelist = .true.
else
  has_noah_namelist = .false.
endif

if (has_noah_namelist) then
   allocate(textblock(nlines))
   textblock = ''

   call nc_check(nf90_def_dim(ncid=ncFileID, name='noahNMLnlines', &
          len = nlines, dimid = nlinesDimID), &
          'nc_write_model_atts', 'def_dim noahNMLnlines '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name=trim(noah_namelist_filename), xtype=nf90_char,    &
          dimids = (/ linelenDimID, nlinesDimID /),  varid=nmlVarID), &
          'nc_write_model_atts', 'def_var noah_namelist '//trim(filename))

   call nc_check(nf90_put_att(ncFileID, nmlVarID, 'long_name', 'contents of '//trim(noah_namelist_filename)), &
          'nc_write_model_atts', 'put_att noah_namelist '//trim(filename))

endif

!-------------------------------------------------------------------------------
! Here is the extensible part. The simplest scenario is to output the state vector,
! parsing the state vector into model-specific parts is complicated, and you need
! to know the geometry, the output variables (PS,U,V,T,Q,...) etc. We're skipping
! complicated part.
!-------------------------------------------------------------------------------

if ( output_state_vector ) then

   !----------------------------------------------------------------------------
   ! Create a variable for the state vector
   !----------------------------------------------------------------------------

  ! Define the state vector coordinate variable and some attributes.
   call nc_check(nf90_def_var(ncid=ncFileID,name='StateVariable', xtype=NF90_INT, &
                              dimids=StateVarDimID, varid=StateVarVarID), &
                             'nc_write_model_atts', 'def_var StateVariable')
   call nc_check(nf90_put_att(ncFileID, StateVarVarID,'long_name','State Variable ID'), &
                             'nc_write_model_atts', 'put_att StateVariable long_name')
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, 'units',     'indexical'), &
                             'nc_write_model_atts', 'put_att StateVariable units')
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, 'valid_range', (/ 1, model_size /)), &
                             'nc_write_model_atts', 'put_att StateVariable valid_range')

   ! Define the actual (3D) state vector, which gets filled as time goes on ...
   call nc_check(nf90_def_var(ncid=ncFileID, name='state', xtype=NF90_REAL, &
                 dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), &
                 varid=StateVarID), 'nc_write_model_atts', 'def_var state')
   call nc_check(nf90_put_att(ncFileID, StateVarID, 'long_name', 'model state or fcopy'), &
                             'nc_write_model_atts', 'put_att state long_name')

   ! Leave define mode so we can fill the coordinate variable.
   call nc_check(nf90_enddef(ncfileID),'nc_write_model_atts', 'state_vector enddef')

   ! Fill the state variable coordinate variable
   call nc_check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /)), &
                                    'nc_write_model_atts', 'put_var state')

else

   !----------------------------------------------------------------------------
   ! We need to output the prognostic variables.
   !----------------------------------------------------------------------------
   ! Define the additional dimensions IDs
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_dim(ncid=ncFileID, name='west_east', &
          len=west_east, dimid=weDimID),'nc_write_model_atts','west_east def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='south_north', &
          len=south_north, dimid=snDimID),'nc_write_model_atts', 'south_north def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='soil_layers_stag',  &
          len=nsoil, dimid=nsoilDimID), 'nc_write_model_atts', 'def_dim soil_layers_stag'//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Coordinate Variables and the Attributes
   !----------------------------------------------------------------------------

   ! Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name='XLONG', xtype=nf90_real, &
                 dimids=(/ weDimID, snDimID /), varid=VarID),&
                 'nc_write_model_atts', 'XLONG def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'coordinate longitude'), &
                 'nc_write_model_atts', 'XLONG long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'coordinates', 'XLONG XLAT'),  &
                 'nc_write_model_atts', 'XLONG coordinates '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'FieldType', 104),  &
                 'nc_write_model_atts', 'XLONG FieldType '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'MemoryOrder', 'XY'),  &
                 'nc_write_model_atts', 'XLONG MemoryOrder '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, '_FillValue', -999.0 ), &
                 'nc_write_model_atts', 'XLONG _FillValue '//trim(filename))

   ! Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name='XLAT', xtype=nf90_real, &
                 dimids=(/ weDimID, snDimID /), varid=VarID),&
                 'nc_write_model_atts', 'XLAT def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'coordinate latitude'), &
                 'nc_write_model_atts', 'XLAT long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'coordinates', 'XLONG XLAT'),  &
                 'nc_write_model_atts', 'XLAT coordinates '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'FieldType', 104),  &
                 'nc_write_model_atts', 'XLAT FieldType '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'MemoryOrder', 'XY'),  &
                 'nc_write_model_atts', 'XLAT MemoryOrder '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, '_FillValue', -999.0 ), &
                 'nc_write_model_atts', 'XLAT _FillValue '//trim(filename))

   ! subsurface levels
   call nc_check(nf90_def_var(ncFileID,name='soil_layers_stag', xtype=nf90_real, &
                 dimids=(/ nsoilDimID /), varid=VarID),&
                 'nc_write_model_atts', 'soil_layers_stag def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'coordinate soil levels'), &
                 'nc_write_model_atts', 'soil_layers_stag long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'm'),  &
                 'nc_write_model_atts', 'soil_layers_stag units '//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Prognostic Variables and their Attributes
   !----------------------------------------------------------------------------

   do ivar=1, nfields

      varname = trim(progvar(ivar)%varname)
      string1 = trim(filename)//' '//trim(varname)

      ! match shape of the variable to the dimension IDs

      call define_var_dims(ivar, ncFileID, MemberDimID, unlimitedDimID, myndims, mydimids)

      ! define the variable and set the attributes

      call nc_check(nf90_def_var(ncid=ncFileID, name=trim(varname), xtype=progvar(ivar)%xtype, &
                    dimids = mydimids(1:myndims), varid=VarID),&
                    'nc_write_model_atts', trim(string1)//' def_var' )

      call nc_check(nf90_put_att(ncFileID, VarID, 'long_name', trim(progvar(ivar)%long_name)), &
           'nc_write_model_atts', trim(string1)//' put_att long_name' )

      call nc_check(nf90_put_att(ncFileID, VarID, 'DART_kind', trim(progvar(ivar)%kind_string)), &
           'nc_write_model_atts', trim(string1)//' put_att dart_kind' )
      call nc_check(nf90_put_att(ncFileID, VarID, 'units', trim(progvar(ivar)%units)), &
           'nc_write_model_atts', trim(string1)//' put_att units' )

      ! Preserve the original missing_value/_FillValue code.

!     if (  progvar(ivar)%xtype == NF90_INT ) then
!        call nc_check(nf90_put_att(ncFileID, VarID, 'missing_value', progvar(ivar)%spvalINT), &
!             'nc_write_model_atts', trim(string1)//' put_att missing_value' )
!        call nc_check(nf90_put_att(ncFileID, VarID, '_FillValue',  progvar(ivar)%spvalINT), &
!             'nc_write_model_atts', trim(string1)//' put_att _FillValue' )

!     elseif (  progvar(ivar)%xtype == NF90_FLOAT ) then
!        call nc_check(nf90_put_att(ncFileID, VarID, 'missing_value', progvar(ivar)%spvalR4), &
!             'nc_write_model_atts', trim(string1)//' put_att missing_value' )
!        call nc_check(nf90_put_att(ncFileID, VarID, '_FillValue',  progvar(ivar)%spvalR4), &
!             'nc_write_model_atts', trim(string1)//' put_att _FillValue' )

!     elseif (  progvar(ivar)%xtype == NF90_DOUBLE ) then
!        call nc_check(nf90_put_att(ncFileID, VarID, 'missing_value', progvar(ivar)%spvalR8), &
!             'nc_write_model_atts', trim(string1)//' put_att missing_value' )
!        call nc_check(nf90_put_att(ncFileID, VarID, '_FillValue',  progvar(ivar)%spvalR8), &
!             'nc_write_model_atts', trim(string1)//' put_att _FillValue' )
!     endif

   enddo

   !----------------------------------------------------------------------------
   ! Finished with dimension/variable definitions, must end 'define' mode to fill.
   !----------------------------------------------------------------------------

   call nc_check(nf90_enddef(ncfileID), 'nc_write_model_atts', 'prognostic enddef')

   !----------------------------------------------------------------------------
   ! Fill the coordinate variables
   !----------------------------------------------------------------------------

   call nc_check(nf90_inq_varid(ncFileID, 'XLONG', VarID), &
                'nc_write_model_atts', 'inq_varid XLONG '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, xlong ), &
                'nc_write_model_atts', 'put_var XLONG '//trim(filename))

   call nc_check(nf90_inq_varid(ncFileID, 'XLAT', VarID), &
                'nc_write_model_atts', 'inq_varid XLAT '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, xlat ), &
                'nc_write_model_atts', 'put_var XLAT '//trim(filename))

   call nc_check(nf90_inq_varid(ncFileID, 'soil_layers_stag', VarID), &
                'nc_write_model_atts', 'inq_varid soil_layers_stag '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, zsoil(1:nsoil)), &
                'nc_write_model_atts', 'put_var soil_layers_stag '//trim(filename))
endif

!-------------------------------------------------------------------------------
! Fill the variables we can
!-------------------------------------------------------------------------------

if (has_noah_namelist) then
   call file_to_text(noah_namelist_filename, textblock)
   call nc_check(nf90_put_var(ncFileID, nmlVarID, textblock ), &
                 'nc_write_model_atts', 'put_var nmlVarID')
   deallocate(textblock)
endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID),'nc_write_model_atts', 'sync')

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts




function nc_write_model_vars( ncFileID, state_vec, copyindex, timeindex ) result (ierr)
!------------------------------------------------------------------
! All errors are fatal, so the return code is always '0 == normal'
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: state_vec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME)          :: varname
character(len=NF90_MAX_NAME),dimension(NF90_MAX_VAR_DIMS) :: dimnames
integer :: i, ivar, VarID, ncNdims, dimlen, numdims, timedimcounter
integer :: TimeDimID, CopyDimID

real(r8), allocatable, dimension(:)       :: data_1d_array
real(r8), allocatable, dimension(:,:)     :: data_2d_array
real(r8), allocatable, dimension(:,:,:)   :: data_3d_array

character(len=128) :: filename

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename,*) 'ncFileID', ncFileID

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file,
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncFileID, 'copy', dimid=CopyDimID), &
            'nc_write_model_vars', 'inq_dimid copy '//trim(filename))

call nc_check(nf90_inq_dimid(ncFileID, 'time', dimid=TimeDimID), &
            'nc_write_model_vars', 'inq_dimid time '//trim(filename))

if ( output_state_vector ) then

   !----------------------------------------------------------------------------
   ! simply blast out the vector
   !----------------------------------------------------------------------------

   call nc_check(nf90_inq_varid(ncFileID, 'state', VarID), &
                               'nc_write_model_vars', 'inq_varid state' )
   call nc_check(nf90_put_var(ncFileID, VarID, state_vec,  &
                              start=(/ 1, copyindex, timeindex /)), &
                             'nc_write_model_vars', 'put_var state')

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   !----------------------------------------------------------------------------

   do ivar = 1,nfields

      varname = trim(progvar(ivar)%varname)
      string2 = trim(filename)//' '//trim(varname)

      ! Ensure netCDF variable is conformable with progvar quantity.
      ! The TIME and Copy dimensions are intentionally not queried.
      ! This requires that Time is the unlimited dimension (the last one in Fortran),
      ! and that 'copy' is the second-to-last. The variables declared in the DART
      ! diagnostic files are required to have the same shape as in the source
      ! restart file. If Time is present there, it must also be the 'last' one.

      ! FIXME ... somewhere I should ensure that IF time is present in the original
      ! prognostic variable from the model, it is the last/unlimited dimension.

      call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
            'nc_write_model_vars', 'inq_varid '//trim(string2))

      call nc_check(nf90_inquire_variable(ncFileID,VarID,dimids=dimIDs,ndims=ncNdims), &
            'nc_write_model_vars', 'inquire '//trim(string2))

      timedimcounter = 0
      mystart(:) = 1
      mycount(:) = 1
      DimCheck : do i = 1,ncNdims

         write(string1,'(A,i2,A)') 'inquire dimension ',i,trim(string2)
         call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), name=dimnames(i), len=dimlen), &
               'nc_write_model_vars', trim(string1))

         if (dimIDs(i) == CopyDimID) cycle DimCheck
         if (dimIDs(i) == TimeDimID) then
            timedimcounter = 1
            cycle DimCheck
         endif

         if ( dimlen /= progvar(ivar)%dimlens(i) ) then
            write(string1,*) trim(string2),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
            write(string2,*)' but it should be.'
            call error_handler(E_ERR, 'nc_write_model_vars', trim(string1), &
                            source, revision, revdate, text2=trim(string2))
         endif

         mycount(i) = dimlen

      enddo DimCheck

      where(dimIDs == CopyDimID) mystart = copyindex
      where(dimIDs == CopyDimID) mycount = 1
      where(dimIDs == TimeDimID) mystart = timeindex
      where(dimIDs == TimeDimID) mycount = 1

      if ( (do_output()) .and. debug > 10 ) then
         write(*,*)'nc_write_model_vars '//trim(varname)//' start is ',mystart(1:ncNdims)
         write(*,*)'nc_write_model_vars '//trim(varname)//' count is ',mycount(1:ncNdims)
         write(*,'(A20)')'nc_write_model_vars ',dimnames(1:progvar(ivar)%numdims)
      endif

      ! If the original variable is shaped:
      !     XXXXXX(time, somedimension) -or- XXXXXX(somedimension)
      ! then it is a 1D variable in our context.
      ! If it is shaped
      !     XXXXXX(time, south_north, west_east) -or- XXXXXX(south_north, west_east)
      ! it really is 2D ...
      !
      ! this adjustment to numdims below is to remove the Time dimension

      numdims = progvar(ivar)%numdims - timedimcounter

      if ( numdims == 1 ) then

         if ( ncNdims /= 3 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 3 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate, text2=string2)
         endif

         allocate(data_1d_array( progvar(ivar)%dimlens(1) ))
         call vector_to_prog_var(state_vec, ivar, data_1d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_1d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_1d_array)

      elseif ( numdims == 2 ) then

         if ( ncNdims /= 4 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 4 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate, text2=string2)
         endif

         allocate(data_2d_array( progvar(ivar)%dimlens(1),  &
                                 progvar(ivar)%dimlens(2) ))
         call vector_to_prog_var(state_vec, ivar, data_2d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_2d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_2d_array)

      elseif ( numdims == 3 ) then

         if ( ncNdims /= 5 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 5 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate, text2=string2)
         endif

         allocate(data_3d_array( progvar(ivar)%dimlens(1),  &
                                 progvar(ivar)%dimlens(2),  &
                                 progvar(ivar)%dimlens(3) ))
         call vector_to_prog_var(state_vec, ivar, data_3d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_3d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_3d_array)

      else

         write(string1,*)'do not know how to handle NOAH variables with more than 3 non-time dimensions'
         write(string2,*)trim(progvar(ivar)%varname),' has dimensions ', progvar(ivar)%dimnames(1:progvar(ivar)%numdims)
         call error_handler(E_ERR,'nc_write_model_vars',string1, &
                    source,revision,revdate,text2=string2)

      endif

   enddo

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID), 'nc_write_model_vars', 'sync')

ierr = 0 ! If we got here, things went well.

end function nc_write_model_vars



subroutine pert_model_state(state, pert_state, interf_provided)
!------------------------------------------------------------------
!
! Perturbs a model state for generating initial ensembles.
! The perturbed state is returned in pert_state.
! A model may choose to provide a NULL INTERFACE by returning
! .false. for the interf_provided argument. This indicates to
! the filter that if it needs to generate perturbed states, it
! may do so by adding an O(0.1) magnitude perturbation to each
! model state variable independently. The interf_provided argument
! should be returned as .true. if the model wants to do its own
! perturbing of states.  The returned pert_state should in any
! case be valid, since it will be read by filter even if
! interf_provided is .false.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

if ( .not. module_initialized ) call static_init_model

call error_handler(E_ERR,'pert_model_state', &
                  'NOAH cannot be started from a single vector', &
                  source, revision, revdate, &
                  text2='see comments in noah/model_mod.f90::pert_model_state()',&
                  text3='or noah/model_mod.html#pert_model_state')

interf_provided = .false.

end subroutine pert_model_state



subroutine ens_mean_for_model(ens_mean)
!------------------------------------------------------------------

real(r8), intent(in) :: ens_mean(:)

if ( .not. module_initialized ) call static_init_model

end subroutine ens_mean_for_model


!==================================================================
! PUBLIC interfaces that aren't required by the DART code but are
! generally useful for other related utility programs.
! (less necessary for small models; generally used for larger models
! with predefined file formats and control structures.)
!==================================================================

function get_debug_level()
   integer :: get_debug_level

   get_debug_level = debug

end function



subroutine verify_state_variables( ncid, filename, ngood )
!------------------------------------------------------------------

integer,                          intent(in)  :: ncid
character(len=*),                 intent(in)  :: filename
integer,                          intent(out) :: ngood

integer :: i, VarID
character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME) :: dartstr

if ( .not. module_initialized ) call static_init_model

ngood = 0
MyLoop : do i = 1, MAX_STATE_VARIABLES

   varname    = trim(noah_state_variables(1,i))
   dartstr    = trim(noah_state_variables(2,i))

   if ( varname == ' ' .and. dartstr == ' ' ) exit MyLoop ! Found end of list.

   if ( varname == ' ' .or. dartstr == ' ' ) then
      string1 = 'model_nml:noah_state_variables not fully specified'
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Make sure variable exists in netCDF file

   write(string1,'(''there is no variable '',a,'' in '',a)') trim(varname), trim(filename)
   call nc_check(NF90_inq_varid(ncid, trim(varname), VarID), &
                 'verify_state_variables', trim(string1))

   ! Make sure DART kind is valid

   if( get_index_for_quantity(dartstr) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Record the contents of the DART state vector

   if (do_output() .and. (debug > 0)) then
      write(logfileunit,*)'variable ',i,' is ',trim(varname), '   ', trim(dartstr)
      write(     *     ,*)'variable ',i,' is ',trim(varname), '   ', trim(dartstr)
   endif

   ngood = ngood + 1
enddo MyLoop

if (ngood == MAX_STATE_VARIABLES) then
   string1 = 'WARNING: There is a possibility you need to increase ''MAX_STATE_VARIABLES'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'verify_state_variables',string1,source,revision,revdate,text2=string2)
endif

end subroutine verify_state_variables



subroutine get_noah_restart_filename( noah_restart_filename )
!------------------------------------------------------------------
character(len=*), intent(out) :: noah_restart_filename

if ( .not. module_initialized ) call static_init_model

noah_restart_filename = noah_netcdf_filename

end subroutine get_noah_restart_filename



subroutine noah_to_dart_vector(filename, state_vector, restart_time)
!------------------------------------------------------------------
! Reads the current time and state variables from a model data
! file and packs them into a dart state vector.

character(len=*), intent(in)  :: filename
real(r8),         intent(out) :: state_vector(:)
type(time_type),  intent(out) :: restart_time

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, ncstart, nccount
character(len=NF90_MAX_NAME) :: dimname, varname

integer  :: ncid, ncNdims, dimlen, VarID
integer  :: i, indx1, indx2, indx3, indx4, indx, ivar, ntimes

real(r8), allocatable, dimension(:)       :: data_1d_array
real(r8), allocatable, dimension(:,:)     :: data_2d_array
real(r8), allocatable, dimension(:,:,:)   :: data_3d_array
real(r8), allocatable, dimension(:,:,:,:) :: data_4d_array

if ( .not. module_initialized ) call static_init_model

state_vector(:) = MISSING_R8

! Check that the input file exists ...

if ( .not. file_exist(filename) ) then
   write(string1,*) 'file <', trim(filename),'> does not exist.'
   call error_handler(E_ERR,'noah_to_dart_vector',string1,source,revision,revdate)
endif

call nc_check(nf90_open(adjustl(filename), NF90_NOWRITE, ncid), &
                   'noah_to_dart_vector', 'open '//trim(filename))

restart_time = get_state_time(ncid, trim(filename))

if ( (do_output()) .and. debug > 99 ) then
   call print_date(restart_time,'noah_to_dart_vector:date of restart file '//trim(filename))
   call print_time(restart_time,'noah_to_dart_vector:time of restart file '//trim(filename))
endif

! Start counting and filling the state vector one item at a time,
! repacking the Nd arrays into a single 1d list of numbers.

do ivar=1, nfields

   ntimes     = -1
   ncstart(:) = -1
   nccount(:) = -1
   varname    = trim(progvar(ivar)%varname)
   string3    = trim(filename)//' '//trim(varname)

   call nc_check(nf90_inq_varid(ncid, varname, VarID), &
            'noah_to_dart_vector', 'inq_varid '//trim(string3))
   call nc_check(nf90_inquire_variable(ncid,VarID,dimids=dimIDs,ndims=ncNdims), &
            'noah_to_dart_vector', 'inquire '//trim(string3))

   ! Check the shape of the variable

   if ( ncNdims /= progvar(ivar)%numdims ) then
      write(string1, *) 'netCDF rank of '//trim(varname)//' does not agree with internal rank.'
      write(string2, *) 'netCDF rank is ',ncNdims,' expected ',progvar(ivar)%numdims
      write(string3, *) 'should not happen'
      call error_handler(E_ERR,'noah_to_dart_vector', string1, &
                        source,revision,revdate,text2=string2,text3=string3)
   endif

   ! Check the memory order of the variable
   ! making sure we only compare the last timestep ...

   do i = 1,ncNdims
      write(string1,'(''inquire dimension'',i2,A)') i,trim(string3)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), name=dimname, len=dimlen), &
                        'noah_to_dart_vector',string1)

      if (trim(dimname) /= trim(progvar(ivar)%dimnames(i))) then
         write(string1, *) 'netCDF dimnames of '//trim(varname)//' does not match expected dimname'
         write(string2, *) 'netCDF dimname is ',trim(dimname),' expected ',trim(progvar(ivar)%dimnames(i))
         write(string3, *) 'should not happen'
         call error_handler(E_ERR,'noah_to_dart_vector', string1, &
                           source,revision,revdate,text2=string2,text3=string3)
      endif

      ncstart(i) = 1
      nccount(i) = dimlen

      if ( trim(dimname) == 'Time' ) then
         ntimes     = dimlen
         ncstart(i) = dimlen
         nccount(i) = 1
      elseif ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(string3),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
         call error_handler(E_ERR,'noah_to_dart_vector',string1,source,revision,revdate)
      endif
   enddo

   if ( (do_output()) .and. debug > 99 ) then
      write(*,*)
      write(*,*)'noah_to_dart_vector: variable ',trim(varname)
      write(*,*)'noah_to_dart_vector: ncstart ',ncstart(1:ncNdims)
      write(*,*)'noah_to_dart_vector: nccount ',nccount(1:ncNdims)
      write(*,*)
   endif

   ! FIXME - this is probably the place to ensure that the Time dimension is the last
   ! dimension if it is present. unlimited dimension

   ! Pack the variable into the DART state vector

   indx = progvar(ivar)%index1

   if (ncNdims == 1) then

      allocate(data_1d_array(nccount(1)))

      call nc_check(nf90_get_var(ncid, VarID, data_1d_array,  &
                     start=ncstart(1:ncNdims), count=nccount(1:ncNdims)), &
                  'noah_to_dart_vector', 'get_var '//trim(string3))

      do indx1 = 1, nccount(1)
         state_vector(indx) = data_1d_array(indx1)
         indx = indx + 1
      enddo
      deallocate(data_1d_array)

   elseif (ncNdims == 2) then

      allocate(data_2d_array(nccount(1), nccount(2)))

      call nc_check(nf90_get_var(ncid, VarID, data_2d_array,  &
                     start=ncstart(1:ncNdims), count=nccount(1:ncNdims)), &
                  'noah_to_dart_vector', 'get_var '//trim(string3))

      do indx2 = 1, nccount(2)
      do indx1 = 1, nccount(1)
         state_vector(indx) = data_2d_array(indx1,indx2)
         indx = indx + 1
      enddo
      enddo
      deallocate(data_2d_array)

   elseif (ncNdims == 3) then

      allocate(data_3d_array(nccount(1), nccount(2), nccount(3)))

      call nc_check(nf90_get_var(ncid, VarID, data_3d_array,  &
                     start=ncstart(1:ncNdims), count=nccount(1:ncNdims)), &
                  'noah_to_dart_vector', 'get_var '//trim(string3))

      do indx3 = 1, nccount(3)
      do indx2 = 1, nccount(2)
      do indx1 = 1, nccount(1)
         state_vector(indx) = data_3d_array(indx1,indx2,indx3)
         indx = indx + 1
      enddo
      enddo
      enddo
      deallocate(data_3d_array)

   elseif (ncNdims == 4) then

      allocate(data_4d_array(nccount(1), &
                             nccount(2), &
                             nccount(3), &
                             nccount(4)))

      call nc_check(nf90_get_var(ncid, VarID, data_4d_array,  &
                     start=ncstart(1:ncNdims), count=nccount(1:ncNdims)), &
                  'noah_to_dart_vector', 'get_var '//trim(string3))

      ! TJH RAFAEL transform goes here.

      do indx4 = 1, nccount(4)
      do indx3 = 1, nccount(3)
      do indx2 = 1, nccount(2)
      do indx1 = 1, nccount(1)
         state_vector(indx) = data_4d_array(indx1,indx2,indx3,indx4)
         indx = indx + 1
      enddo
      enddo
      enddo
      enddo
      deallocate(data_4d_array)

   else

      write(string1, *)'Variable '//trim(varname)//' has ',ncNdims,' dimensions.'
      write(string2, *)'cannot handle that.'
      call error_handler(E_ERR,'noah_to_dart_vector', string1, &
                        source,revision,revdate,text2=string2)

   endif

   indx = indx - 1
   if ( indx /= progvar(ivar)%indexN ) then
      write(string1, *)'Variable '//trim(varname)//' filled wrong.'
      write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',indx
      call error_handler(E_ERR,'noah_to_dart_vector', string1, &
                        source,revision,revdate,text2=string2)
   endif

enddo

if ( (do_output()) .and. debug > 99 ) then
   write(*,*)'newest time is ',ntimes
   do i = 1,size(state_vector)
      write(*,*)'state vector(',i,') is',state_vector(i)
   enddo
endif

end subroutine noah_to_dart_vector



subroutine dart_vector_to_model_file(state_vector, filename, dart_time, skip_variables)
!------------------------------------------------------------------
! Writes the current time and state variables from a dart state
! vector (1d array) into a noah netcdf restart file.
!
! This is VERY similar to nc_write_model_vars() for this model.
! If it were not for the 'copy' dimension, it would be identical, I think.

real(r8),         intent(in) :: state_vector(:)
character(len=*), intent(in) :: filename
type(time_type),  intent(in) :: dart_time
character(len=*), intent(in) :: skip_variables(:)

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME)          :: varname
character(len=NF90_MAX_NAME),dimension(NF90_MAX_VAR_DIMS) :: dimnames

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: ncFileID, VarID, ncNdims, TimeDimID
integer :: timeindex, dimlen, numdims, timedimcounter

type(time_type) :: file_time

! temp space to hold data while we are writing it
integer :: i, ivar
real(r8), allocatable, dimension(:)         :: data_1d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array
real(r8), allocatable, dimension(:,:,:)     :: data_3d_array

if ( .not. module_initialized ) call static_init_model

! Check that the output file exists ...

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for writing.'
   call error_handler(E_ERR,'dart_vector_to_model_file',string1,source,revision,revdate)
endif

call nc_check(nf90_open(trim(filename), NF90_WRITE, ncFileID), &
                  'dart_vector_to_model_file','open '//trim(filename))

call nc_check(nf90_inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID), &
                  'dart_vector_to_model_file', 'inquire '//trim(filename))

call nc_check(nf90_inq_dimid(ncFileID, 'Time', TimeDimID), &
                  'dart_vector_to_model_file','inq_dimid Time '//trim(filename))

! make sure the time in the file is the same as the time on the data
! we are trying to insert.  we are only updating part of the contents
! of the NOAH restart file, and state vector contents from a different
! time won't be consistent with the rest of the file.

file_time = get_state_time(ncFileID, trim(filename), timeindex)

if ( file_time /= dart_time ) then
   call print_time(dart_time,'DART current time',logfileunit)
   call print_time(file_time,'NOAH current time',logfileunit)
   call print_time(dart_time,'DART current time')
   call print_time(file_time,'NOAH current time')
   write(string1,*)trim(filename),' current time /= model time. FATAL error.'
   call error_handler(E_ERR,'dart_vector_to_model_file',string1,source,revision,revdate)
endif

if (do_output() .and. (debug > 0)) then
   call print_date(file_time,'dart_vector_to_model_file: date of restart file '//trim(filename))
   call print_time(file_time,'dart_vector_to_model_file: time of restart file '//trim(filename))
endif

! The DART prognostic variables are only defined for a single time.
! IF the netCDF variable has a TIME dimension, it must be the last dimension.

UPDATE : do ivar=1, nfields

   varname = trim(progvar(ivar)%varname)
   string2 = trim(filename)//' '//trim(varname)

   ! If this variable is on the skip list ... skip it.

   SKIPME : do i = 1,size(skip_variables)
      if (len_trim(skip_variables(i)) < 1) cycle SKIPME
      if (skip_variables(i) == varname) cycle UPDATE
   enddo SKIPME

   ! Ensure netCDF variable is conformable with DART progvar quantity.

   call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
            'dart_vector_to_model_file', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncFileID,VarID,dimids=dimIDs,ndims=ncNdims), &
            'dart_vector_to_model_file', 'inquire '//trim(string2))

   timedimcounter = 0
   mystart(:) = 1
   mycount(:) = 1
   DimCheck : do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), name=dimnames(i), len=dimlen), &
            'dart_vector_to_model_file', string1)

      if (dimIDs(i) == TimeDimID) timedimcounter = 1

      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(string2),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
         write(string2,*)' but it should be.'
         call error_handler(E_ERR, 'dart_vector_to_model_file', string1, &
                         source, revision, revdate, text2=string2)
      endif

      mycount(i) = dimlen

   enddo DimCheck

   if (dimIDs(ncNdims) /= TimeDimID) then
      write(string1,*) trim(string2),' required to have "Time" as the last/unlimited dimension'
      write(string2,*)' last dimension is ',trim(dimnames(ncNdims))
      call error_handler(E_ERR, 'dart_vector_to_model_file', string1, &
                      source, revision, revdate, text2=string2)
   endif

   where(dimIDs == TimeDimID) mystart = timeindex
   where(dimIDs == TimeDimID) mycount = 1

   if ( (do_output()) .and. debug > 99 ) then
      write(*,*)'dart_vector_to_model_file '//trim(varname)//' start is ',mystart(1:ncNdims)
      write(*,*)'dart_vector_to_model_file '//trim(varname)//' count is ',mycount(1:ncNdims)
      write(*,*)'dart_vector_to_model_file ',dimnames(1:progvar(ivar)%numdims)
   endif

   numdims = progvar(ivar)%numdims - timedimcounter

   if ( numdims == 1 ) then

      allocate(data_1d_array(progvar(ivar)%dimlens(1)))
      call vector_to_prog_var(state_vector, ivar, data_1d_array,limit=.true.)
      call nc_check(nf90_put_var(ncFileID, VarID, data_1d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'dart_vector_to_model_file', 'put_var '//trim(string2))
      deallocate(data_1d_array)

   elseif ( numdims == 2 ) then

      allocate(data_2d_array(progvar(ivar)%dimlens(1), &
                             progvar(ivar)%dimlens(2)))
      call vector_to_prog_var(state_vector, ivar, data_2d_array,limit=.true.)
      call nc_check(nf90_put_var(ncFileID, VarID, data_2d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'dart_vector_to_model_file', 'put_var '//trim(string2))
      deallocate(data_2d_array)

   elseif ( numdims == 3) then

      allocate(data_3d_array(progvar(ivar)%dimlens(1), &
                             progvar(ivar)%dimlens(2), &
                             progvar(ivar)%dimlens(3)))
      call vector_to_prog_var(state_vector, ivar, data_3d_array,limit=.true.)
      ! TJH RAFAEL undo transform goes here.
      call nc_check(nf90_put_var(ncFileID, VarID, data_3d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'dart_vector_to_model_file', 'put_var '//trim(string2))
      deallocate(data_3d_array)

   else

      write(string1, *) 'no support for data array of dimension ', ncNdims
      call error_handler(E_ERR,'dart_vector_to_model_file', string1, &
                        source,revision,revdate)
   endif

   ! Make note that the variable has been updated by DART
   call nc_check(nf90_Redef(ncFileID),'dart_vector_to_model_file', 'redef '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VarID,'DART','variable modified by DART'),&
                 'dart_vector_to_model_file', 'modified '//trim(varname))
   call nc_check(nf90_enddef(ncfileID),'dart_vector_to_model_file','state enddef '//trim(filename))

enddo UPDATE

call nc_check(nf90_close(ncFileID),'dart_vector_to_model_file','close '//trim(filename))
ncFileID = 0

end subroutine dart_vector_to_model_file



function get_state_time(ncid, filename, timeindex)
!------------------------------------------------------------------
! The restart netcdf files have the time of the state.
! We are always using the 'most recent' which is, by defn, the last one.
!
! The way the HRLDAS driver works is a bit wonky.
! The time in the restart file is NOT the time at which the state is valid.
! It is one noah_timestep AHEAD of the valid time.
!
! for instance, if the noah_timestep is 3600 seconds, the restart_frequency_hours is 1,
! and the filename is RESTART.2004010102_DOMAIN1 the
!
!        Time = UNLIMITED ; // (blah_blah_blah currently)
!        DateStrLen = 19 ;
!variables:
!        char Times(Time, DateStrLen) ;
!
! Times =
!  '2004-01-01_02:00:00' ;
!
! BUT - the data is for the previous noah_timestep ... i.e. 2004-01-01_01:00:00
! No kidding.

type(time_type) :: get_state_time
integer,           intent(in)  :: ncid
character(len=*),  intent(in)  :: filename
integer, optional, intent(out) :: timeindex

character(len=19), allocatable, dimension(:) :: datestring
integer               :: year, month, day, hour, minute, second
integer               :: DimID, VarID, strlen, ntimes
type(time_type)       :: filetime, timestep

if ( .not. module_initialized ) call static_init_model

! Get the dimensions for the strings of times

call nc_check(nf90_inq_dimid(ncid, 'Time', DimID), &
                  'get_state_time','inq_dimid Time '//trim(filename))
call nc_check(nf90_inquire_dimension(ncid, DimID, len=ntimes), &
                  'get_state_time','inquire_dimension Time '//trim(filename))

call nc_check(nf90_inq_dimid(ncid, 'DateStrLen', DimID), &
                  'get_state_time','inq_dimid DateStrLen '//trim(filename))
call nc_check(nf90_inquire_dimension(ncid, DimID, len=strlen), &
                  'get_state_time','inquire_dimension DateStrLen '//trim(filename))

if (strlen /= len(datestring)) then
   write(string1,*)'DatStrLen string length ',strlen,' /= ',len(datestring)
   call error_handler(E_ERR,'get_state_time', string1, source, revision, revdate)
endif

! Get all the Time strings, use the last one.

call nc_check(nf90_inq_varid(ncid, 'Times', VarID), &
                   'get_state_time', 'inq_varid Times '//trim(filename))

allocate(datestring(ntimes))

call nc_check(nf90_get_var(ncid, VarID, datestring), &
                   'get_state_time', 'get_var Times '//trim(filename))

read(datestring(ntimes),'(i4,5(1x,i2))')year, month, day, hour, minute, second

timestep       = set_time(noah_timestep, 0)
filetime       = set_date(year, month, day, hours=hour, minutes=minute, seconds=second)
get_state_time = filetime - timestep

if (present(timeindex)) timeindex = ntimes

if ( (do_output()) .and. debug > 99 ) write(*,*)'get_state_time: Last time string is '//trim(datestring(ntimes))
if ( (do_output()) .and. debug > 99 ) call print_date(get_state_time,' get_state_time: means valid time is ')

deallocate(datestring)

end function get_state_time



subroutine get_hrldas_constants(filename)
!------------------------------------------------------------------
! Read the 'wrfinput' netCDF file for grid information, etc.
! This is all time-invariant, so we can mostly ignore the Time coordinate.
!
! MODULE variables set by this routine:
!    south_north
!    west_east
!    xlong
!    xlat

character(len=*), intent(in) :: filename

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, ncstart, nccount
character(len=NF90_MAX_NAME)          :: dimname

integer :: i, iunit, DimID, VarID, numdims, dimlen, xtype

if ( .not. module_initialized ) call static_init_model

call nc_check(nf90_open(adjustl(filename), NF90_NOWRITE, iunit), &
                   'get_hrldas_constants', 'open '//trim(filename))

call nc_check(nf90_inq_dimid(iunit, 'south_north', DimID), &
                  'get_hrldas_constants','inq_dimid south_north '//trim(filename))
call nc_check(nf90_inquire_dimension(iunit, DimID, len=south_north), &
                  'get_hrldas_constants','inquire_dimension south_north '//trim(filename))

call nc_check(nf90_inq_dimid(iunit, 'west_east', DimID), &
                  'get_hrldas_constants','inq_dimid west_east '//trim(filename))
call nc_check(nf90_inquire_dimension(iunit, DimID, len=west_east), &
                  'get_hrldas_constants','inquire_dimension west_east '//trim(filename))

! Require that the xlong and xlat are the same shape.

allocate(xlong(west_east,south_north), xlat(west_east,south_north))

call nc_check(nf90_inq_varid(iunit, 'XLONG', VarID), &
                  'get_hrldas_constants','inq_varid XLONG '//trim(filename))

call nc_check(nf90_inquire_variable(iunit, VarID, dimids=dimIDs, &
                  ndims=numdims, xtype=xtype), &
                  'get_hrldas_constants', 'inquire_variable XLONG '//trim(filename))

! Form the start/count such that we always get the 'latest' time.

ncstart(:) = 0
nccount(:) = 0

do i = 1,numdims

   write(string1,'(''inquire dimension'',i2,A)') i,trim(filename)
   call nc_check(nf90_inquire_dimension(iunit, dimIDs(i), name=dimname, len=dimlen), &
                                          'get_hrldas_constants', string1)
   ncstart(i) = 1
   nccount(i) = dimlen

   if ((trim(dimname) == 'Time') .or. (trim(dimname) == 'time')) then
      ncstart(i) = dimlen
      nccount(i) = 1
   endif

enddo

if ( (do_output()) .and. debug > 99 ) write(*,*)'DEBUG get_hrldas_constants ncstart is',ncstart(1:numdims)
if ( (do_output()) .and. debug > 99 ) write(*,*)'DEBUG get_hrldas_constants nccount is',nccount(1:numdims)

! finally get the longitudes

call nc_check(nf90_get_var(iunit, VarID, xlong, &
                  start=ncstart(1:numdims), count=nccount(1:numdims)), &
                  'get_hrldas_constants', 'get_var XLONG '//trim(filename))

where(xlong < 0.0_r8) xlong = xlong + 360.0_r8
where(xlong == 360.0_r8) xlong = 0.0_r8

! finally get the latitudes

call nc_check(nf90_inq_varid(iunit, 'XLAT', VarID), &
                  'get_hrldas_constants','inq_varid XLAT '//trim(filename))

call nc_check(nf90_get_var(iunit, VarID, xlat, &
                  start=ncstart(1:numdims), count=nccount(1:numdims)), &
                  'get_hrldas_constants', 'get_var XLAT '//trim(filename))

write(string1,*) 'get_hrldas_constants() not verified for multiple gridcells.'
call error_handler(E_MSG,'get_hrldas_constants',string1,source,revision,revdate)

end subroutine get_hrldas_constants



subroutine define_var_dims(ivar, ncid, memberdimid, unlimiteddimid, ndims, dimids)
!------------------------------------------------------------------
! I am trying to preserve the shape of the NOAH variable as much as possible.
!
! noah_variable(Time,       soil_layers_stag, south_north, west_east) becomes
! noah_variable(time, copy, soil_layers_stag, south_north, west_east)
!
! Since 'time' or 'Time' is the unlimited dimension in both ... I can skip it
! in the DEFDIM loop.

integer,               intent(in)  :: ivar, ncid, memberdimid, unlimiteddimid
integer,               intent(out) :: ndims
integer, dimension(:), intent(out) :: dimids

integer :: i, mydimid

ndims = 0

DEFDIM : do i = 1,progvar(ivar)%numdims

   if ((trim(progvar(ivar)%dimnames(i)) == 'Time') .or. &
       (trim(progvar(ivar)%dimnames(i)) == 'time')) cycle DEFDIM

   call nc_check(nf90_inq_dimid(ncid=ncid, name=progvar(ivar)%dimnames(i), dimid=mydimid), &
                           'define_var_dims','inq_dimid '//trim(progvar(ivar)%dimnames(i)))

   ndims         = ndims + 1
   dimids(ndims) = mydimid

enddo DEFDIM

ndims = ndims + 1
dimids(ndims) = memberdimid
ndims = ndims + 1
dimids(ndims) = unlimitedDimid

if ( (do_output()) .and. debug > 99 ) then
   write(logfileunit,*)
   write(logfileunit,*)'define_var_dims knowledge'
   write(logfileunit,*)trim(progvar(ivar)%varname),' has dimnames ', &
                       progvar(ivar)%dimnames(1:progvar(ivar)%numdims)
   write(logfileunit,*)' thus dimids ',dimids(1:ndims)
   write(     *     ,*)
   write(     *     ,*)'define_var_dims knowledge'
   write(     *     ,*)trim(progvar(ivar)%varname),' has dimnames ', &
                       progvar(ivar)%dimnames(1:progvar(ivar)%numdims)
   write(     *     ,*)' thus dimids ',dimids(1:ndims)

endif

return
end subroutine define_var_dims



subroutine vector_to_1d_prog_var(x, ivar, data_1d_array, limit)
!------------------------------------------------------------------
! convert the values from a 1d array, starting at an offset, into a 1d array.
!
! If the optional argument (ncid) is specified, some additional
! processing takes place. The variable in the netcdf is read.
! This must be the same shape as the intended output array.
! Anywhere the DART MISSING code is encountered in the input array,
! the corresponding (i.e. original) value from the netCDF file is
! used.

real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: ivar
real(r8), dimension(:),   intent(out) :: data_1d_array
logical,  OPTIONAL,       intent(in)  :: limit

integer :: i,ii

if ( .not. module_initialized ) call static_init_model

! unpack the right part of the DART state vector into a 1D array.

ii = progvar(ivar)%index1

do i = 1, progvar(ivar)%dimlens(1)
   data_1d_array(i) = x(ii)
   ii = ii + 1
enddo

if (present(limit)) then
if ( limit ) then
   if (    progvar(ivar)%rangeRestricted == 1) then
      where(data_1d_array < progvar(ivar)%minvalue) data_1d_array = progvar(ivar)%minvalue
   elseif (progvar(ivar)%rangeRestricted == 2) then
      where(data_1d_array > progvar(ivar)%maxvalue) data_1d_array = progvar(ivar)%maxvalue
   elseif (progvar(ivar)%rangeRestricted == 3) then
      where(data_1d_array < progvar(ivar)%minvalue) data_1d_array = progvar(ivar)%minvalue
      where(data_1d_array > progvar(ivar)%maxvalue) data_1d_array = progvar(ivar)%maxvalue
   endif
endif
endif

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' packed wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_1d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine vector_to_1d_prog_var



subroutine vector_to_2d_prog_var(x, ivar, data_2d_array, limit)
!------------------------------------------------------------------
! convert the values from a 1d array, starting at an offset,
! into a 2d array.
!
real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: ivar
real(r8), dimension(:,:), intent(out) :: data_2d_array
logical,  OPTIONAL,       intent(in)  :: limit

integer :: i,j,ii

if ( .not. module_initialized ) call static_init_model

! unpack the right part of the DART state vector into a 2D array.

ii = progvar(ivar)%index1

do j = 1,progvar(ivar)%dimlens(2)
do i = 1,progvar(ivar)%dimlens(1)
   data_2d_array(i,j) = x(ii)
   ii = ii + 1
enddo
enddo

if (present(limit)) then
if ( limit ) then
   if (    progvar(ivar)%rangeRestricted == 1) then
      where(data_2d_array < progvar(ivar)%minvalue) data_2d_array = progvar(ivar)%minvalue
   elseif (progvar(ivar)%rangeRestricted == 2) then
      where(data_2d_array > progvar(ivar)%maxvalue) data_2d_array = progvar(ivar)%maxvalue
   elseif (progvar(ivar)%rangeRestricted == 3) then
      where(data_2d_array < progvar(ivar)%minvalue) data_2d_array = progvar(ivar)%minvalue
      where(data_2d_array > progvar(ivar)%maxvalue) data_2d_array = progvar(ivar)%maxvalue
   endif
endif
endif

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_2d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine vector_to_2d_prog_var



subroutine vector_to_3d_prog_var(x, ivar, data_3d_array, limit)
!------------------------------------------------------------------
! convert the values from a 1d array, starting at an offset,
! into a 3d array.
!
real(r8), dimension(:),     intent(in)  :: x
integer,                    intent(in)  :: ivar
real(r8), dimension(:,:,:), intent(out) :: data_3d_array
logical,  OPTIONAL,         intent(in)  :: limit

integer :: i,j,k,ii

if ( .not. module_initialized ) call static_init_model

! unpack the right part of the DART state vector into a 3D array.

ii = progvar(ivar)%index1

do k = 1,progvar(ivar)%dimlens(3)
do j = 1,progvar(ivar)%dimlens(2)
do i = 1,progvar(ivar)%dimlens(1)
   data_3d_array(i,j,k) = x(ii)
   ii = ii + 1
enddo
enddo
enddo

if (present(limit)) then
if ( limit ) then
   if (    progvar(ivar)%rangeRestricted == 1) then
      where(data_3d_array < progvar(ivar)%minvalue) data_3d_array = progvar(ivar)%minvalue
   elseif (progvar(ivar)%rangeRestricted == 2) then
      where(data_3d_array > progvar(ivar)%maxvalue) data_3d_array = progvar(ivar)%maxvalue
   elseif (progvar(ivar)%rangeRestricted == 3) then
      where(data_3d_array < progvar(ivar)%minvalue) data_3d_array = progvar(ivar)%minvalue
      where(data_3d_array > progvar(ivar)%maxvalue) data_3d_array = progvar(ivar)%maxvalue
   endif
endif
endif

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_3d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine vector_to_3d_prog_var



subroutine get_noah_timestepping(day,hour,dynamical,output,forcing,restart)
integer,          intent(out) :: day,hour,dynamical,output,forcing,restart

day       = kday
hour      = khour
dynamical = noah_timestep
output    = output_timestep
forcing   = forcing_timestep
restart   = restart_frequency_hours*3600

end subroutine


!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
