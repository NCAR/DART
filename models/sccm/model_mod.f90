! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module model_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
! $Name$ 

! Model mod for the single column climate model for Yuqiong Liu. The model
! is a single column of CCM 3.6 combined with a simple land surface model
! column under it. The land surface column has a number of tiles and
! a number of levels. An arbitrary number of atmospheric tracers, the first
! one being q, are carried.

! Modules that are absolutely required for use are listed
use        types_mod, only : r8
use time_manager_mod, only : time_type, set_time, get_time, get_date, &
                             set_date, operator(+), set_calendar_type, &
                             GREGORIAN
use     location_mod, only : location_type, set_location, query_location
use    utilities_mod, only : file_exist, open_file, close_file, &
                             register_module, error_handler, E_ERR, E_MSG, logfileunit

implicit none
private

public :: get_model_size, &
          adv_1step, &
          get_state_meta_data, &
          model_interpolate, &
          get_model_time_step, &
          end_model, &
          static_init_model, &
          init_time, &
          init_conditions, &
          model_get_close_states, &
          nc_write_model_atts, &
          nc_write_model_vars, &
          pert_model_state, &
          write_sccm_state, &
          read_sccm_state

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"


! Basic model parameters controlled by nameslist; have defaults

!---------------------------------------------------------------
! Namelist with default values

integer  :: kpt        = 2          ! Number of tiles per land surface box
integer  :: msl        = 6          ! Number of soil levels
integer  :: plev       = 18         ! Number of atmospheric levels
integer  :: pcnst      = 1          ! Number of tracers
integer  :: time_step_days = 0
integer  :: time_step_seconds = 1200
logical  :: output_state_vector = .true.

namelist /model_nml/ kpt, msl, plev, pcnst, &
     time_step_days, time_step_seconds, output_state_vector

!----------------------------------------------------------------
! Define the location of the state variables in module storage

type(location_type), allocatable :: state_loc(:)
type(time_type) :: time_step
integer :: model_size          ! Global storage for size of model

! Storage for keeping updated pressure levels and soil level locations
real(r8), allocatable :: pressure(:), depth(:), temp_state(:)


contains

!==================================================================



subroutine static_init_model()
!------------------------------------------------------------------
!
! Called to do one time initialization of the model. 
! For sccm, this entails reading the namelist to get the size and 
! configuration details of the model.

real(r8) :: x_loc
integer  :: i, j, iunit, io, indx, g_offset, s_offset, tile
type(time_type) :: temp_time
character(len=129) :: err_string, nml_string

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Begin by reading the namelist input
if(file_exist('input.nml')) then
   iunit = open_file('input.nml', action = 'read')
   read(iunit, nml = model_nml, iostat = io)
   if(io /= 0) then
      ! A non-zero return means a bad entry was found for this namelist
      ! Reread the line into a string and print out a fatal error message.
      BACKSPACE iunit
      read(iunit, '(A)') nml_string
      write(err_string, *) 'INVALID NAMELIST ENTRY: ', trim(adjustl(nml_string))
      call error_handler(E_ERR, 'static_init_model:&model_nml problem', &
                         err_string, source, revision, revdate)
   endif
   call close_file(iunit)
endif

! Record the namelist values used for the run ...
call error_handler(E_MSG,'static_init_model','model_nml values are',' ',' ',' ')
write(logfileunit, nml=model_nml)
write(     *     , nml=model_nml)

! Compute the model state vector size
model_size = (3 + pcnst) * plev + 1 + 4*kpt + 2*kpt*msl
! Create storage for locations
! Also storage to keep track of pressure of the atmospheric levels and 
! depth of soil levels
allocate(state_loc(model_size), pressure(plev), depth(msl), temp_state(model_size))

! Yuqiong: NOTICE: NEED TO GET THE ABILITY TO COMPUTE THE PRESSURE LEVELS SOMEHOW
! SAME FOR SOIL IN THE LONG RUN
! COULD READ A's AND B's FROM NETCDF FILE FOR MODEL AND COMPUTE EACH TIME???


! Define the locations of the model state variables
! For small model like this assume no localization will be used
! Many variables are at surface, so begin by setting all to surface
do i = 1, model_size
   state_loc(i) =  set_location(0.0_r8, 0.0_r8, 1.0_r8, -1)
end do

! Put in vertical level for u, v, t, tracers
do i = 1, 3 + pcnst
   do j = 1, plev
      indx = (i - 1) * plev + j
      state_loc(indx) = set_location(0.0_r8, 0.0_r8, 1.0_r8 * j, 1)
   end do
end do


! The soil moisture and temperature have soil levels for vertical
! Offsets for where the ground and soil state variables start in dart
g_offset = (3 + pcnst) * plev + 2
s_offset = g_offset + 4 * kpt
do i = 1, 2
   do tile = 1, kpt
      do j = 1, msl
         indx = s_offset + (i - 1) * kpt * msl + (tile - 1) * msl + (j - 1)
         state_loc(indx) = set_location(0.0_r8, 0.0_r8, 1.0_r8 * j, 1)
      end do
   end do
end do

! The time_step in terms of a time type must also be initialized. 
time_step = set_time(time_step_seconds, time_step_days)

end subroutine static_init_model



subroutine init_conditions(x)
!------------------------------------------------------------------
! subroutine init_conditions(x)
!
! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no 
! synthetic data experiments using perfect_model_obs are planned, 
! this can be a NULL INTERFACE.

real(r8), intent(out) :: x(:)


end subroutine init_conditions



subroutine adv_1step(x, time)
!------------------------------------------------------------------
! subroutine adv_1step(x, time)
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

real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: time

end subroutine adv_1step



function get_model_size()
!------------------------------------------------------------------
!
! Returns the size of the model as an integer.  Model size is computed
! in static_init_model from namelist parameters.

integer :: get_model_size

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
! instance temperature, zonal wind component, etc.). For this model,
! the kinds are specified as in CAM
!		U is itype 1
!		V is itype 2
! 		PS is itype 3
!		T is itype 4
!		Q is itype 5

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: itype
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

integer :: lo, hi, i, indx, base_indx, which_vert
real(r8) :: fract, obs_press

! Default for successful return
istatus = 0

! Only types 1 to 5 currently supported
if(itype < 1 .or. itype > 5) then
   istatus = 1
   obs_val = -1_r8
   return
endif

! If itype is pressure, just return ps
if(itype == 3) then
   indx = (3 + pcnst) * plev + 1
   obs_val = x(indx)
   return
endif

! Only relevant part of location is vertical, need pressure for kind
which_vert = query_location(location, 'which_vert')
! YUQIONG: NEED TO PUT CHECKS FOR VERTICAL LOCATION TYPE HERE
obs_press = query_location(location, 'vloc')

write(*, *) 'obs_press is ', obs_press
write(*, *) 'pressure(1), plev ', pressure(1), pressure(plev)
if(1 == 1) stop

! Search through the pressure level array and get upper, lower and fraction
! If above the top or below the bottom, can't use for now
if(obs_press < pressure(1) .or. obs_press > pressure(plev)) then
   istatus = 1
   obs_val = -1_r8
   return
endif

do i = 2, plev
   if(obs_press < pressure(i)) then
      lo = i - 1
      hi = i
      fract = (obs_press - pressure(lo)) / (pressure(hi) - pressure(lo))
      goto 111
   endif
end do

! Pressure location found, now get correct type of observation
111 continue
if(itype == 1) then
   base_indx = 0
else if(itype == 2) then
   base_indx = plev
else if(itype == 4) then
   base_indx = 2*plev
else if(itype == 5) then
   base_indx = 3*plev
endif

lo = lo + base_indx
hi = lo + 1

! Do the pressure interpolation (linear) and return    
obs_val = x(lo) + fract * (x(hi) - x(lo))

write(*, *) 'in model_interpolate ', lo, hi, fract, obs_val
if(1 == 1) stop

end subroutine model_interpolate



function get_model_time_step()
!------------------------------------------------------------------
!
! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.

type(time_type) :: get_model_time_step

get_model_time_step = time_step

end function get_model_time_step



subroutine get_state_meta_data(index_in, location, var_type)
!------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated location.  For the SCCM, we are assuming for now that no
! localization will be used and all the locations are the same. 
! The type is returned to distinguish the different variable types
! with the following values (these should later be changed to
! use the global kind module: JLA):
!
!   1 = u,  2 = v,  3 = t,   4 = tracer,   5 = ps,  
!   6 = h2osno,   7 = h2ocan,   8 = tv,   9 = tg,   
!   10 = h2osoi,   11 = tsoi
!   

integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

integer :: g_index, s_index

location = state_loc(index_in)
if(present(var_type)) then
   ! Set indices for the ground part of the model
   g_index = index_in - (3 + pcnst)*plev + 1
   s_index = g_index - 4*kpt 
   if(index_in <= plev) then
      var_type = 1    ! u column comes first
   else if(index_in <= 2*plev) then
      var_type = 2    ! Then v column
   else if(index_in <= 3*plev) then
      var_type = 3    ! Then t column
   else if(index_in <= (3 + pcnst) *plev) then
      var_type = 4    ! Then pcnst tracer columns
   else if(index_in == 4*plev + 1) then
      var_type = 5    ! ps
   else if(g_index <= kpt) then
      var_type = 6    ! Snow water per tile
   else if(g_index <= 2*kpt) then
      var_type = 7    ! Water in canopy
   else if(g_index <= 3*kpt) then
      var_type = 8    ! Temperature of vegetation
   else if(g_index <= 4*kpt) then
      var_type = 9    ! Surface skin temperature
   else if(s_index <= kpt*msl) then
      var_type = 10   ! Soil moisture, by level and tile
   else if(s_index <= 2*kpt*msl) then
      var_type = 11   ! Temperature of soil by level and tile
   else
      call error_handler(E_ERR,'get_state_meta_data',&   
         'index_in is not in range of model size', source, revision, revdate)
   endif

endif


end subroutine get_state_meta_data



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

end subroutine end_model



subroutine model_get_close_states(o_loc, radius, inum, indices, dist, x)
!------------------------------------------------------------------
! 
! Computes a list of model state variable indices that are within 
! distance radius of a given location, o_loc. The units of the radius
! and the metric for computing distances is defined by the location module
! that is in  use. The number of state variables that are within radius
! of o_loc is returned in inum. The indices of each of these is 
! returned in indices and the corresponding distance in dist. The model
! state is available in x because it is sometimes required to determine
! the distance (for instance, the current model surface pressure field
! is required to compute the location of state variables in a sigma
! vertical coordinate model). A model can choose to do no computation
! here and return a value of -1 in inum. If this happens, the filter
! will do a naive search through ALL state variables for close states.
! This can work fine in low-order models, but can be far too expensive
! in large models.

type(location_type), intent(in) :: o_loc
real(r8), intent(in) :: radius
integer, intent(out) :: inum, indices(:)
real(r8), intent(out) :: dist(:)
real(r8), intent(in) :: x(:)

! Simplest interface just sets inum to -1 and returns
inum = -1

end subroutine model_get_close_states



function nc_write_model_atts( ncFileID ) result (ierr)
!------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! TJH Jan 24 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! For the lorenz_96 model, each state variable is at a separate location.
! that's all the model-specific attributes I can think of ...
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

use typeSizes
use netcdf

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID


integer :: StateVarDimID, StateVarVarID, StateVarID
integer :: MemberDimID,     MemberVarID   ! one per ensemble member
integer :: TimeDimID,         TimeVarID
integer ::   TileDimID,   TileVarID   ! tiles per land surface box
integer ::   SoilDimID,   SoilVarID   ! soil levels
integer ::  AtmosDimID,  AtmosVarID   ! atmospheric levels
integer :: TracerDimID, TracerVarID   ! tracers
integer :: uVarID         ! var 1     ! U          [plev]
integer :: vVarID         ! var 2     ! V          [plev]
integer :: tVarID         ! var 3     ! T          [plev]
integer :: trcrVarID      ! var 4     !            [plev,pcnst]
integer :: psVarID        ! var 5     ! surface pressure
integer ::  h2osnoVarID   ! var 6     ! snow water per tile [kpt]
integer ::  h2ocanVarID   ! var 7     ! water in canopy     [kpt]
integer :: TvegVarID      ! var 8     ! temperature of vegetation   [kpt]
integer :: TskinVarID     ! var 9     ! surface skin temperature    [kpt]
integer :: h2osoilVarID   ! var 10    ! soil moisture      [plev, kpt]
integer :: TsoilVarID     ! var 11    ! soil temperature   [plev, kpt]

character(len=129) :: errstring
character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1
integer :: i

!------------------------------------------------------------------

ierr = 0  ! assume normal termination 

!--------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file
!--------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))
call check(nf90_Redef(ncFileID))

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies
!-------------------------------------------------------------------------------
                                                                                                           
call check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID))
call check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid=  TimeDimID))
                                                                                                           
if ( TimeDimID /= unlimitedDimId ) then
   write(errstring,*)'Time Dimension ID ',TimeDimID,' should equal Unlimited Dimension ID',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts', errstring, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
!-------------------------------------------------------------------------------
call check(nf90_def_dim(ncid=ncFileID, name="StateVariable", &
                        len=model_size, dimid = StateVarDimID))

!-------------------------------------------------------------------------------
! Write Global Attributes
!-------------------------------------------------------------------------------
                                                                                                           
call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date",str1))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source",source))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate",revdate))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model","SCCM"))

!-------------------------------------------------------------------------------
! Define the dimensions IDs
!-------------------------------------------------------------------------------

call check(nf90_def_dim(ncid=ncFileID, name="tiles",    len = kpt,   dimid = TileDimID))
call check(nf90_def_dim(ncid=ncFileID, name="soillev",  len = msl,   dimid = SoilDimID))
call check(nf90_def_dim(ncid=ncFileID, name="atmoslev", len = plev,  dimid = AtmosDimID))
call check(nf90_def_dim(ncid=ncFileID, name="tracers",  len = pcnst, dimid = TracerDimID))

!-------------------------------------------------------------------------------
! Create the (empty) Variables and the Attributes
!-------------------------------------------------------------------------------

call check(nf90_def_var(ncFileID, name="tiles", xtype=nf90_int, &
                                                dimids=TileDimID, varid=TileVarID) ) 
call check(nf90_put_att(ncFileID, TileVarID, "long_name", "tiles per land surface box"))
!call check(nf90_put_att(ncFileID, TileVarID, "units", "n/a"))
!call check(nf90_put_att(ncFileID, TileVarID, "valid_range", (/ 0.0_r8, 360.0_r8 /)))


call check(nf90_def_var(ncFileID, name="soillev", xtype=nf90_int, &
                                                dimids=SoilDimID, varid=SoilVarID) ) 
call check(nf90_put_att(ncFileID, SoilVarID, "long_name", "soil levels"))
!call check(nf90_put_att(ncFileID, SoilVarID, "units", "degrees_north"))
!call check(nf90_put_att(ncFileID, SoilVarID, "valid_range", (/ -90.0_r8, 90.0_r8 /)))

call check(nf90_def_var(ncFileID, name="atmoslev", xtype=nf90_double, &
                                                dimids=AtmosDimID, varid=AtmosVarID) ) 
call check(nf90_put_att(ncFileID, AtmosVarID, "long_name", "atmospheric levels"))
!call check(nf90_put_att(ncFileID, AtmosVarID, "units", "hPa"))
!call check(nf90_put_att(ncFileID, AtmosVarID, "positive", "down"))

call check(nf90_def_var(ncFileID, name="tracers", xtype=nf90_int, &
                                                dimids=TracerDimID, varid=TracerVarID) ) 
call check(nf90_put_att(ncFileID, TracerVarID, "long_name", "tracers"))


if ( output_state_vector ) then

   !----------------------------------------------------------------------------
   ! Create attributes for the state vector
   !----------------------------------------------------------------------------

   ! Define the state vector coordinate variable
   call check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=nf90_int, &
              dimids=StateVarDimID, varid=StateVarVarID))
   call check(nf90_put_att(ncFileID, StateVarVarID, "long_name", "State Variable ID"))
   call check(nf90_put_att(ncFileID, StateVarVarID, "units",     "indexical") )
   call check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", (/ 1, model_size /)))


   ! Define the actual state vector
   call check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_real, &
              dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), &
              varid = StateVarID))
   call check(nf90_put_att(ncFileID, StateVarID, "long_name", "model state or fcopy"))


   ! Leave define mode so we can fill
   call check(nf90_enddef(ncfileID))

   ! Fill the state variable coordinate variable
   call check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ))


else

   ! U wind component
   call check(nf90_def_var(ncid=ncFileID, name="u", xtype=nf90_real, &
          dimids = (/ AtmosDimID, MemberDimID, unlimitedDimID /), varid  = uVarID))
   call check(nf90_put_att(ncFileID, uVarID, "long_name", "zonal wind component"))
   call check(nf90_put_att(ncFileID, uVarID, "units", "m/s"))


   ! V wind component
   call check(nf90_def_var(ncid=ncFileID, name="v", xtype=nf90_real, &
          dimids = (/ AtmosDimID, MemberDimID, unlimitedDimID /), varid  = vVarID))
   call check(nf90_put_att(ncFileID, vVarID, "long_name", "meridional wind component"))
   call check(nf90_put_att(ncFileID, vVarID, "units", "m/s"))


   ! Temperature
   call check(nf90_def_var(ncid=ncFileID, name="t", xtype=nf90_real, &
          dimids = (/ AtmosDimID, MemberDimID, unlimitedDimID /), varid  = tVarID))
   call check(nf90_put_att(ncFileID, tVarID, "long_name", "temperature"))
   !call check(nf90_put_att(ncFileID, tVarID, "units", "degrees Kelvin"))


   ! Tracers
   call check(nf90_def_var(ncid=ncFileID, name="tracer", xtype=nf90_real, &
          dimids = (/ AtmosDimID, TracerDimID, MemberDimID, unlimitedDimID /), &
          varid  = trcrVarID))
   !call check(nf90_put_att(ncFileID, trcrVarID, "long_name", "zonal wind component"))
   !call check(nf90_put_att(ncFileID, trcrVarID, "units", "m/s"))

   ! Surface Pressure
   call check(nf90_def_var(ncid=ncFileID, name="ps", xtype=nf90_real, &
          dimids = (/ MemberDimID, unlimitedDimID /), varid  = psVarID))
   call check(nf90_put_att(ncFileID, psVarID, "long_name", "surface pressure"))
   !call check(nf90_put_att(ncFileID, psVarID, "units", "m/s"))


   ! Snow water (per tile)
   call check(nf90_def_var(ncid=ncFileID, name="h2osno", xtype=nf90_real, &
          dimids = (/ TileDimID, MemberDimID, unlimitedDimID /), varid  = h2osnoVarID))
   call check(nf90_put_att(ncFileID, h2osnoVarID, "long_name", "snow water"))
   !call check(nf90_put_att(ncFileID, h2osnoVarID, "units", "degrees Kelvin"))


   ! Water in Canopy (per tile)
   call check(nf90_def_var(ncid=ncFileID, name="h2ocan", xtype=nf90_real, &
          dimids = (/ TileDimID, MemberDimID, unlimitedDimID /), varid  = h2ocanVarID))
   call check(nf90_put_att(ncFileID, h2ocanVarID, "long_name", "water in canopy"))
   !call check(nf90_put_att(ncFileID, h2ocanVarID, "units", "?"))


   ! Temperature of the Vegetation (per tile)
   call check(nf90_def_var(ncid=ncFileID, name="tv", xtype=nf90_real, &
          dimids = (/ TileDimID, MemberDimID, unlimitedDimID /), varid  = TvegVarID))
   call check(nf90_put_att(ncFileID, TvegVarID, "long_name", "temperature of vegetation"))
   !call check(nf90_put_att(ncFileID, TvegVarID, "units", "?"))


   ! Surface Skin Temperature
   call check(nf90_def_var(ncid=ncFileID, name="tg", xtype=nf90_real, &
          dimids = (/ TileDimID, MemberDimID, unlimitedDimID /), varid  = TskinVarID))
   call check(nf90_put_att(ncFileID, TskinVarID, "long_name", "surface skin temperature"))
   call check(nf90_put_att(ncFileID, TskinVarID, "units", "?"))


   ! Soil Moisture, by level and tile  (one tile-all levels, ... next tile)
   call check(nf90_def_var(ncid=ncFileID, name="h2osoi", xtype=nf90_real, &
          dimids = (/ SoilDimID, TileDimID, MemberDimID, unlimitedDimID /), &
          varid  = h2osoilVarID))
   call check(nf90_put_att(ncFileID, h2osoilVarID, "long_name", "Soil Moisture"))
   !call check(nf90_put_att(ncFileID, h2osoilVarID, "units", "?"))


   ! Soil Temperature, by level and tile  (one tile-all levels, ... next tile)
   call check(nf90_def_var(ncid=ncFileID, name="tsoi", xtype=nf90_real, &
          dimids = (/ SoilDimID, TileDimID, MemberDimID, unlimitedDimID /), &
          varid  = TsoilVarID))
   call check(nf90_put_att(ncFileID, TsoilVarID, "long_name", "Soil Temperature"))
   !call check(nf90_put_att(ncFileID, TsoilVarID, "units", "?"))


   call check(nf90_enddef(ncfileID))

endif

!-------------------------------------------------------------------------------
! Fill the coordinate variables
!-------------------------------------------------------------------------------

call check(nf90_put_var(ncFileID, TileVarID,   (/ (i, i=1,  kpt) /) ))
call check(nf90_put_var(ncFileID, SoilVarID,   (/ (i, i=1,  msl) /) ))
call check(nf90_put_var(ncFileID, AtmosVarID,  (/ (i, i=1, plev) /) ))
call check(nf90_put_var(ncFileID, TracerVarID, (/ (i, i=1,pcnst) /) ))

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call check(nf90_sync(ncFileID))
write (*,*)'nc_write_model_atts: netCDF file ',ncFileID,' is synched ...'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains
                                                                                                           
! Internal subroutine - checks error status after each netcdf, prints
!                       text message each time an error code is returned.

subroutine check(istatus)
    integer, intent (in) :: istatus
    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'nc_write_model_atts', &
          trim(nf90_strerror(istatus)), source, revision, revdate)
end subroutine check


end function nc_write_model_atts


function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         
!------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! TJH 23 May 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! For the lorenz_96 model, each state variable is at a separate location.
! that's all the model-specific attributes I can think of ...
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

use typeSizes
use netcdf

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID
integer :: u1VarID, v1VarID, t1VarID, uVarID, vVarID, tVarID
integer :: qnHVarID, qnOHVarID, qnOVarID
!type(model_type)                   :: Var

ierr = 0  ! assume normal termination

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))


write(*,*)'nc_write_model_vars: time index is ',timeindex 

if ( output_state_vector ) then

   call check(NF90_inq_varid(ncFileID, "state", StateVarID) )
   call check(NF90_put_var(ncFileID, StateVarID, statevec,  &
                start=(/ 1, copyindex, timeindex /)))

else

   ! This block is STRONGLY dependent on the information I gleaned from 
   ! get_state_meta_data. Should anything in get_state_meta_data change,
   ! those changes must occur here also. TJH 31 Aug 2005


!  call vector_to_prog_var(statevec, Var)
!  
!  call check(NF90_inq_varid(ncFileID,  "u1",  u1VarID))
!  call check(nf90_put_var( ncFileID,  u1VarId, var%vars_3d(:,:,:, 1), &
!             start=(/ 1, 1, 1, copyindex, timeindex /) ))
!                                                                                                                                       
!  call check(NF90_inq_varid(ncFileID,  "v1",  v1VarID))
!  call check(nf90_put_var( ncFileID,  v1VarId, var%vars_3d(:,:,:, 2), &
!             start=(/ 1, 1, 1, copyindex, timeindex /) ))
!                                                                                                                                       
!  call check(NF90_inq_varid(ncFileID,  "t1",  t1VarID))
!  call check(nf90_put_var( ncFileID,  t1VarId, var%vars_3d(:,:,:, 3), &
!             start=(/ 1, 1, 1, copyindex, timeindex /) ))
!  
!  call check(NF90_inq_varid(ncFileID,  "u",  uVarID))
!  call check(nf90_put_var( ncFileID,  uVarId, var%vars_3d(:,:,:, 4), &
!             start=(/ 1, 1, 1, copyindex, timeindex /) ))
!  
!  call check(NF90_inq_varid(ncFileID,  "v",  vVarID))
!  call check(nf90_put_var( ncFileID,  vVarId, var%vars_3d(:,:,:, 5), &
!             start=(/ 1, 1, 1, copyindex, timeindex /) ))
!  
!  call check(NF90_inq_varid(ncFileID,  "t",  tVarID))
!  call check(nf90_put_var( ncFileID,  tVarId, var%vars_3d(:,:,:, 6), &
!             start=(/ 1, 1, 1, copyindex, timeindex /) ))
!  
!  call check(NF90_inq_varid(ncFileID,  "qnH",  qnHVarID))
!  call check(nf90_put_var( ncFileID,  qnHVarId, var%vars_3d(:,:,:, 7), &
!             start=(/ 1, 1, 1, copyindex, timeindex /) ))
!  
!  call check(NF90_inq_varid(ncFileID,  "qnOH",  qnOHVarID))
!  call check(nf90_put_var( ncFileID,  qnOHVarId, var%vars_3d(:,:,:, 8), &
!             start=(/ 1, 1, 1, copyindex, timeindex /) ))
!  
!  call check(NF90_inq_varid(ncFileID,  "qnH",  qnOVarID))
!  call check(nf90_put_var( ncFileID,  qnOVarId, var%vars_3d(:,:,:, 9), &
!             start=(/ 1, 1, 1, copyindex, timeindex /) ))

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

write (*,*)'Finished filling variables ...'
call check(nf90_sync(ncFileID))
write (*,*)'netCDF file is synched ...'

!call end_model_instance(Var)   ! should avoid any memory leaking

contains
                                                                                                         
 ! Internal subroutine - checks error status after each netcdf, prints
 !                       text message each time an error code is returned.
 subroutine check(istatus)
  integer, intent ( in) :: istatus
    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'nc_write_model_vars', &
          trim(nf90_strerror(istatus)), source, revision, revdate)
                                                                                                         
 end subroutine check

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
! perturbing of states.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

interf_provided = .false.

end subroutine pert_model_state


subroutine write_sccm_state(state, time, delta_seconds, file_id)
!------------------------------------------------------------------
!

real(r8), intent(in) :: state(:)
integer, intent(in) :: delta_seconds, file_id
type(time_type), intent(in) :: time

real(r8) :: atmos(plev), sfc(kpt), soil(msl * kpt)
integer :: g_offset, s_offset, i, date_stamp
integer :: year, month, day, hour, minute, second, days, seconds

! Write the delta seconds
write(file_id, *) "advance-time ", delta_seconds

! First write day as date integer and seconds
call get_date(time, year, month, day, hour, minute, second)
write(*, *) 'get_date', year, month, day, hour, minute, second
! Create and write the date stamp for the sccm
year = year - year / 100 * 100
date_stamp = year * 10000 + month * 100 + day
write(file_id, *) 'date ', date_stamp
call get_time(time, seconds, days)
write(file_id, *) 'sec ', seconds

! Write out the model size descriptors
write(file_id, *) 'KPT ', kpt
write(file_id, *) 'MSL ', msl
write(file_id, *) 'PLEV ', plev
write(file_id, *) 'PCNST ', pcnst

! Write the pressure level values
write(file_id, *) 'press_lev', pressure

! Write the soil level values
write(file_id, *) 'soi_lev', depth


! Write out the atmospheric t
atmos(:) = state(2*plev + 1 : 3*plev)
write(file_id, *) 't3(plev) ', atmos
! Write out atmospheric u
atmos(:) = state(1 : plev) 
write(file_id, *) 'u3(plev) ', atmos
! Write out atmospheric v
atmos(:) = state(plev + 1 : 2*plev)
write(file_id, *) 'v3(plev) ', atmos
! Write out the tracer fields
do i = 1, pcnst
   atmos(:) = state((2 + i)*plev + 1: (3 + i) * plev)
   write(file_id, *) 'q3(pcnst,plev) ', atmos
end do
! Write out the ps
write(file_id, *) 'ps ', state((3 + pcnst)* plev + 1)


! Offsets for where the ground and soil state variables start in dart
g_offset = (3 + pcnst) * plev + 2
s_offset = g_offset + 4 * kpt
! Snow moisture content is first
sfc = state(g_offset : g_offset + kpt - 1)
write(file_id, *) 'h2osno(kpt) ', sfc
! Canopy moisture next
sfc = state(g_offset + kpt : g_offset + 2*kpt - 1)
write(file_id, *) 'h20can(kpt) ', sfc
! Now soil moisture
soil = state(s_offset : s_offset + msl*kpt - 1)
write(file_id, *) 'h2osoi(kpt,msl) ', soil
! Now temperature of vegetation
sfc = state(g_offset + 2*kpt : g_offset + 3*kpt - 1)
write(file_id, *) 'tv(kpt) ', sfc
! Now skin temperature
sfc = state(g_offset + 3*kpt : g_offset + 4*kpt - 1)
write(file_id, *) 'tg(kpt) ', sfc
! Finally, soil temperature
soil = state(s_offset + msl*kpt : s_offset + 2*msl*kpt - 1)
write(file_id, *) 'tsoi(kpt,msl) ', soil

end subroutine write_sccm_state



subroutine read_sccm_state(state, time, file_id)
!------------------------------------------------------------------
!

real(r8), intent(out) :: state(:)
type(time_type), intent(out) :: time
integer, intent(in) :: file_id

real(r8) :: atmos(plev), sfc(kpt), soil(msl * kpt)
integer :: g_offset, s_offset, i, date_stamp
integer :: year, month, day, hour, minute, second, days, seconds
integer :: kpt_in, msl_in, plev_in, pcnst_in, delta_seconds
character * 40 :: temp_string, temp_string2

! Read the delta time entry which is not used
read(file_id, *) temp_string, delta_seconds

! First read day as date integer and seconds
read(file_id, *) temp_string, date_stamp
read(file_id, *) temp_string, seconds
write(*, *) 'date stamp is ', date_stamp
year = date_stamp / 10000
month = (date_stamp - year * 10000) / 100
day = (date_stamp - year * 10000 - month * 100)
! Year needs a 19 or 20 in front
if(year > 50) then
   year = year + 1900
else
   year = year + 2000
endif
write(*, *) 'year month day ', year, month, day

! First, set a date for just the day part with 0 seconds
time = set_date(year, month, day, 0, 0, 0)
! Then add in the seconds
time = time + set_time(seconds, 0)

! Read the model size descriptors, check for errors
read(file_id, *) temp_string, kpt_in
read(file_id, *) temp_string, msl_in
read(file_id, *) temp_string, plev_in
read(file_id, *) temp_string, pcnst_in
if(kpt_in /= kpt .or. msl_in /= msl .or. plev_in /= plev .or. pcnst_in /= pcnst)then
   write(*, *) 'ERROR ON MODELS SIZE FOR READ'
   stop
endif

! Read in the pressure level values
read(file_id, *) temp_string, pressure
! DART USES PRESSURE IN PA, NEED TO MULTIPLY BY 100
pressure = pressure * 100.0_r8

! Read in the soil level values
read(file_id, *) temp_string, depth

! Read the atmospheric t
read(file_id, *) temp_string, atmos
state(2*plev + 1 : 3*plev) = atmos
! Read atmospheric u
read(file_id, *) temp_string, atmos
state(1 : plev) = atmos
! Read atmospheric v
read(file_id, *) temp_string, atmos
state(plev + 1 : 2*plev) = atmos
! Read the tracer fields
do i = 1, pcnst
   read(file_id, *) temp_string, temp_string2, atmos
   state((2 + i)*plev + 1: (3 + i) * plev) = atmos
end do
! Read the ps
read (file_id, *) temp_string, state((3 + pcnst)* plev + 1)


! Offsets for where the ground and soil state variables start in dart
g_offset = (3 + pcnst) * plev + 2
s_offset = g_offset + 4 * kpt
! Snow moisture content is first
read(file_id, *) temp_string, sfc
state(g_offset : g_offset + kpt - 1) = sfc
! Canopy moisture next
read(file_id, *) temp_string, sfc
state(g_offset + kpt : g_offset + 2*kpt - 1) = sfc
! Now soil moisture
read(file_id, *) temp_string, temp_string2, soil
state(s_offset : s_offset + msl*kpt - 1) = soil
! Now temperature of vegetation
read(file_id, *) temp_string, sfc
state(g_offset + 2*kpt : g_offset + 3*kpt - 1) = sfc
! Now skin temperature
read(file_id, *) temp_string, sfc
state(g_offset + 3*kpt : g_offset + 4*kpt - 1) = sfc
! Finally, soil temperature
read(file_id, *) temp_string, temp_string2, soil
state(s_offset + msl*kpt : s_offset + 2*msl*kpt - 1) = soil

end subroutine read_sccm_state


!===================================================================
! End of model_mod
!===================================================================
end module model_mod
