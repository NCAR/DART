! vars_3d still has order of variables; lon,lev,lat,field in some(?) routines
!  (use vars_3d as a search tag if I decide to change over to all lon,lat,lev)
! ARE THEY SEPARATE?
!  no; Tim's routines call vector_to_prog_var to fill Var, as do mine (write_cam_init)
!      these are used for P{oste}rior_Diag.nc; don't change until I verify that
!      assim_model/assim_model_mod.f90 doesn't rely on coordinate order.
!    do we want to change that to a more standard?  nc_write_model_atts, trans_xxx...?

! NOVERT marks modifications for fields with no vertical location,
! i.e. GWD parameters.

!!!!!!!!!!!!!!!!!!!!!!
! WARNING; dimension of 1d state component is set to the largest likely value; num_lons.
!          This may cause problems if the component is size num_lats or num_levs
! FIX  or HOOK (see more below)
!!!!!!!!!!!!!!!!!!!!!!


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
!
! changes for reading hybrid coefficients from initial file are marked with 'coef'
! changes for calculating pressures on cam vertical levels are marked with  'plevs'
! 
! netCDF filename; where will this come from in DART?
!                  it's created by CAM, and written out to a file somewhere; read it in?
!
! Do we need other functionality form bgrid model_mod?

!----------------------------------------------------------------------
! purpose: interface between CAM and DART
!              translate to/from state_vector and caminput.nc,
!              initialize model,
!              write out CAM fields to Prior and Posterior_Diag.nc,
!              generate expected obs from model state

!
! method: Read CAM 'initial' file (netCDF format).
!         Reform fields into a state vector.
!         (Let DART modify values state vector.)
!         Reform state vector back into CAM fields.
!         Replace those fields on the CAM initial file with the new values,
!         preserving all other information on the file.
!         Also read hybrid coordinate coefficients from CAM input file (for plevs_dart)
!
! author: Kevin Raeder 2/14/03  and 8/1/03
!         based on prog_var_to_vector and vector_to_prog_var by Jeff Anderson
!
! modified: Tim Hoar 02 Sep 03 
!         nc_write_model_atts, nc_write_model_vars now write out "prognostic"
!         files instead of a nondescript state variable vector glom
!
! augmented; Kevin Raeder 7/1/04 for CAM3.0 and to use namelist input for
!         lists of fields to include in state vector.
!         'CAM3' marks changes
!         Later;?
!         Also; read field attributes from netcdf file and write them out 
!         from nc_write_model_atts, instead of hard-coded there.
!----------------------------------------------------------------------

use netcdf
use        types_mod, only : r8
use time_manager_mod, only : time_type, set_time, print_time, set_calendar_type, &
                             THIRTY_DAY_MONTHS, JULIAN, GREGORIAN, NOLEAP, NO_CALENDAR
use    utilities_mod, only : file_exist, open_file, check_nml_error, close_file, &
                             register_module, error_handler, E_ERR, E_MSG, logfileunit
use     location_mod, only : location_type, get_location, set_location, &
                             get_dist, vert_is_level, query_location, &
                             LocationDims, LocationName, LocationLName
use     obs_kind_mod, only : KIND_U, KIND_V, KIND_PS, KIND_T, KIND_QV, KIND_P
!    and maybe KIND_W, KIND_QR, KIND_TD, KIND_VR, KIND_REF, KIND_U10, KIND_V10, 
!              KIND_T2, KIND_Q2, KIND_TD2
use    random_nr_mod, only : init_ran1
use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian


implicit none
private

public model_type, prog_var_to_vector, vector_to_prog_var, read_cam_init, &
   read_cam_init_size, init_model_instance, end_model_instance, &
   write_cam_init, get_model_size, static_init_model, &
   get_state_meta_data, get_model_time_step, model_interpolate, &
   init_conditions, init_time, adv_1step, end_model, &
   model_get_close_states, nc_write_model_atts, nc_write_model_vars, &
   TYPE_PS, TYPE_T, TYPE_U, TYPE_V, TYPE_Q, pert_model_state

!-----------------------------------------------------------------------
! CVS Generated file description for error handling, do not edit
character(len=128) :: version = "$Id$"
character(len=128) :: tag = "$Name$"
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Public definition of variable types
! Values will be defined in order_state_fields
! All fields will have entries in the TYPE_xD corresponding to their orders
!   in state_names_Xd.  The explicitly named TYPE_s are for convenience
integer :: TYPE_PS, TYPE_T, TYPE_U, TYPE_V, TYPE_Q
integer, allocatable :: TYPE_1D(:), TYPE_2D(:), TYPE_3D(:)

! CAM3; additional types, or can TRACER handle the new ones?

!----------------------------------------------------------------------

! A type for cam model, very simple for now for conversion only
type model_type
    private
   real(r8), pointer :: vars_0d(:)
   real(r8), pointer :: vars_1d(:, :)
   real(r8), pointer :: vars_2d(:, :, :)
   real(r8), pointer :: vars_3d(:, :, :, :)
end type model_type

!----------------------------------------------------------------------
! Global storage for describing cam model class
integer :: model_size, num_lons, num_lats, num_levs

type(time_type) :: Time_step_atmos

! Random sequence and init for pert_model_state
logical :: first_pert_call = .true.
type(random_seq_type)   :: random_seq

!----------------------------------------------------------------------
! Namelist variables with default values follow

! output_state_vector = .true.     results in a "state-vector" netCDF file
! output_state_vector = .false.    results in a "prognostic-var" netCDF file
logical :: output_state_vector = .false.

! File where basic info about model configuration can be found
character(len = 128) :: model_config_file = 'caminput.nc', &
                        model_version = '3.0'

! Define highest pressure for which obserations are used in mb
real(r8) :: highest_obs_pressure_mb = 30.0

! Namelist variables for defining state vector, and default values
! read in sizes from first namelist, allocate, set default values, then get values from
! second namelist.  
! Or, allocate with defaults values, read in namelist, deallocate and reallocate.

! Better names than CAM code
integer :: state_num_0d = 0              ! # of scalars fields to read from file
integer :: state_num_1d = 0              ! # of 1d fields to read from file
integer :: state_num_2d = 1              ! # of 2d fields to read from file
integer :: state_num_3d = 4              ! # of 3d fields to read from file

! make these allocatable?  Where would I define the default values?
integer :: iii
character (len=8),dimension(20) :: state_names_0d = (/('        ',iii=1,20)/)
character (len=8),dimension(20) :: state_names_1d = (/('        ',iii=1,20)/)
character (len=8),dimension(20) :: state_names_2d = (/'PS      ',('        ',iii=1,19)/)
character (len=8),dimension(20) :: state_names_3d =  &
          (/'T       ','U       ','V       ','Q       ',('        ',iii=1,16)/)

! NOVERT  need (a) namelist parameter(s) to define which_vert for each 2D (xD?) field
!         There's a danger of having a mismatch with the state_names_Xd; should this
!         definition be part of state_names_Xd, which is parsed into a name and which_vert
!         after being read?  Not for now.

! Are these defaults good? or should they all be 0?

integer , dimension(20) :: which_vert_1d = (/(-2,iii=1,20)/)
integer , dimension(20) :: which_vert_2d = (/(-1,iii=1,20)/)
integer , dimension(20) :: which_vert_3d = (/( 1,iii=1,20)/)


! Is there a way to exclude stat_nums from namelist and have those filled in 
! the  subroutine which sorts state_names?
! Yes, use two namelists model_nml_1 and model_nml_2 at future date

! list of fields which trans_pv_sv_pert0 needs to perturb because they're
! constant valued model parameters and show no spread when start_from_restart = .true.
character (len=8),dimension(20) :: state_names_pert = (/('        ',iii=1,20)/)
real(r8), dimension(20) :: state_names_sd = (/(0.,iii=1,20)/)


! Specify shortest time step that the model will support
! This is limited below by CAMs fixed time step but is also impacted
! by numerical stability concerns for repeated restarting in leapfrog.
integer :: Time_step_seconds = 21600, Time_step_days = 0

namelist /model_nml/ output_state_vector , model_version , model_config_file & 
                     ,state_num_0d   ,state_num_1d   ,state_num_2d   ,state_num_3d &
                     ,state_names_0d ,state_names_1d ,state_names_2d ,state_names_3d &
                     ,                which_vert_1d  ,which_vert_2d  ,which_vert_3d &
                     ,state_names_pert ,state_names_sd &
                     ,highest_obs_pressure_mb ,Time_step_seconds ,Time_step_days

!----------------------------------------------------------------------
! derived parameters
integer :: nflds         ! # fields to read

! Arrays to store lats, lons, and gaussian weights
real(r8), allocatable :: lons(:), lats(:), gw(:)

!coef hybrid coeffs at interfaces and mid-levels
real(r8), allocatable :: hyai(:), hybi(:), hyam(:), hybm(:)

!reference pressure for pressure calculations using hybrid coeff 
! should be read from same file as hybrid coeffs?
real(r8):: P0               ! reference pressure

! from ncdump of CAM standard configuration initial file:
!   P0 ; lat = 64 ; lon = 128 ; lev = 26 ; time
!   hyai, hybi, hyam, hybm, gw
!   U, V, T, Q, PS, 
!   (names of advected and nonadv constituents)
!   PHIS, SGH, PBLH, LANDM, 
!   TPERT, QPERT, CLOUD, QCWAT, TCWAT, LCWAT,
!   TSICERAD, TS, TSICE, TS1, TS2, TS3, TS4,
!   SNOWHICE, LANDFRAC,
!   TBOT, ICEFRAC, SICTHK, TSOCN, CLDLIQ, CLDICE,
! list of variables to be included in the state vector,
! according to the category to which they belong, 
! in the order of 
! FIX;
! WARNING; Is this the right order if we include 1d and/or 0d fields in state vector?
! SEE prog_var_to_vec ...
!   (state_num_0d, state_num_1d, state_num_2d, state_num_3d, 
!    state_num_adv_tracers, state_num_notadv_tracers).

! CAM3 array 'cflds' is filled with simple loops over state_names_xxx, 
! which requires that the field names be provided in the right order.
! Later I should replace that with code which orders namelist input field names
! into cflds, regardless of their original order, and tallies how many of each.
! Is there a way to exclude state_nums from namelist and have those filled in 
! the same subroutine?

character (len=8), allocatable :: cflds(:)

! Attribute values for the fields which comprise the state vector
! These are filled by nc_read_model_atts
character (len=128), allocatable :: state_long_names(:)
character (len=128), allocatable :: state_units(:)
! character (len=128), allocatable :: state_units_long_names(:)

! array for the linking of obs_kinds (KIND_) to model field TYPE_s
! It's filled in obs_field_location
! The max size of KIND_ should come from obs_kind_mod
integer, dimension(1000) :: obs_loc_in_sv = (/(0,iii=1,1000)/)
!
!---- namelist (saved in file input.nml) ----
!-----------------------------------------------------------------------

!#######################################################################

contains

!#######################################################################



  subroutine read_cam_init_size(file_name, num_lons, num_lats, num_levs)
!=======================================================================
! subroutine read_cam_init_size(file_name, num_lons, num_lats, num_levs)
!
! Gets the number of lons, lats and levels from a CAM init netcdf file
! in file_name

character(len = *), intent(in) :: file_name
integer, intent(out) :: num_lons, num_lats, num_levs

character (len=NF90_MAX_NAME) :: clon,clat,clev
integer :: londimid, levdimid, latdimid, ncfileid, ncfldid

write(*, *) 'file_name in read_cam_init_size is ', trim(file_name)

!------------------------------------
! read CAM 'initial' file domain info
call check(nf90_open(path = trim(file_name), mode = nf90_write, ncid = ncfileid))

! get dimension 'id's
call check(nf90_inq_dimid(ncfileid, 'lon', londimid))
call check(nf90_inq_dimid(ncfileid, 'lat', latdimid))
call check(nf90_inq_dimid(ncfileid, 'lev', levdimid))

! get dimension sizes
call check(nf90_inquire_dimension(ncfileid, londimid, clon , num_lons ))
call check(nf90_inquire_dimension(ncfileid, latdimid, clat , num_lats ))
call check(nf90_inquire_dimension(ncfileid, levdimid, clev , num_levs ))

contains 

   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus) 
   integer, intent ( in) :: istatus 
   if (istatus /= nf90_noerr) call error_handler(E_ERR, 'read_cam_init_size', &
          trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check

end subroutine read_cam_init_size


  subroutine read_cam_init(file_name, var)
!=======================================================================
! subroutine read_cam_init(file_name, var)
!

character(len = *), intent(in) :: file_name
type(model_type), intent(out) :: var

! Local workspace
integer :: i, j, k, n, m, ifld  ! grid and constituent indices
integer :: plon, plat, plev, num_coord, coord_order

character (len=NF90_MAX_NAME) :: clon,clat,clev
integer :: londimid, levdimid, latdimid, ncfileid, ncfldid
integer :: coord_3d(3), coord_2d(2)
real(r8), allocatable :: temp_3d(:,:,:)

!----------------------------------------------------------------------
! read CAM 'initial' file domain info
call check(nf90_open(path = trim(file_name), mode = nf90_write, ncid = ncfileid))

call size_coord(coord_3d,coord_2d, coord_order)
allocate (temp_3d(coord_3d(1),coord_3d(2),coord_3d(3)))

! read CAM 'initial' file fields desired

ifld = 0
!0d fields; scalars are recognized and handled differently than vectors
!           by NetCDF
! kdr; Tim says that nf90_put_var was probably bombing because
!      netcdf recognized that it was writing a scalar, called an f77
!      scalar put_var, and then choked when it saw the count = ...
!      So, this padding is probably not necessary.
do i= 1, state_num_0d
   ifld = ifld + 1
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
!  fields on file are 1D; TIME(=1)
   call check(nf90_get_var(ncfileid, ncfldid, var%vars_0d(i) ))
end do
             

!1d fields
do i= 1, state_num_1d
   ifld = ifld + 1
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
!  fields on file are 2D; (height, lat, or lon!?) TIME(=1)
!  assume longest for now; lon
   num_coord = coord_3d(1)
   call check(nf90_get_var(ncfileid, ncfldid, var%vars_1d(:, i) &
             ,start=(/1,1/) ,count=(/num_coord, 1/) ))
end do

!2d fields
do i= 1, state_num_2d
   ifld = ifld + 1
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
!  fields on file are 3D; lon, lat, TIME(=1)
   call check(nf90_get_var(ncfileid, ncfldid, var%vars_2d(:, :, i) &
             ,start=(/1,1,1/) ,count=(/coord_2d(1), coord_2d(2), 1/) ))
end do

! 3d fields
do i=1, state_num_3d
   ifld = ifld + 1
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
!  fields on file are 4D; lon, lev, lat, TIME(=1) 
!                     or; lon, lat, lev, TIME
   call check(nf90_get_var(ncfileid, ncfldid, temp_3d  &
             ,start=(/1,1,1,1/) ,count=(/coord_3d(1), coord_3d(2), coord_3d(3),1/) ))
! read into dummy variable, then repackage depending on coord_order
! put it into lon,lev,lat order for continuity with Tim's/historical.
   if (coord_order == 1) then
!     lon,lev,lat as in original CAM
      var%vars_3d(:, :, :, i) = temp_3d
   elseif (coord_order == 2) then
!     lon,lat,lev as in new CAM
      do k=1,coord_3d(3)
      do n=1,coord_3d(2)
      do m=1,coord_3d(1)
         var%vars_3d(m,k,n,i) = temp_3d(m,n,k)
      end do
      end do
      end do
   end if
end do

deallocate (temp_3d)

contains 

   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus) 
   integer, intent ( in) :: istatus 
   if (istatus /= nf90_noerr) call error_handler(E_ERR, 'read_cam_init', &
          trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check

end subroutine read_cam_init


  subroutine size_coord(coord_3d, coord_2d, coord_order)
!=======================================================================
! subroutine size_coord(coord_3d, coord_2d)
!

! Figure out which coordinates are lon, lat, lev, based on CAM version
! from the namelist, which has form #.#[.#[.#]]
integer,           intent(out) :: coord_3d(3), coord_2d(2)
integer, optional, intent(out) :: coord_order

! local workspace
character (len=4)      :: form_version = '(I0)'
character (len=4)      :: char_version
integer                :: part, nchars, tot_chars, i
integer, dimension(4)  :: int_version = (/(0,i=1,4)/)


! Choose order of coordinates based on CAM version
part = 1
nchars=0
char_version = '    '
tot_chars = len_trim(model_version)
do i=1,tot_chars+1
   if ( i == tot_chars+1 .or.  model_version(i:i) == '.' ) then
      write(form_version(3:3),'(I1)') nchars
      read(char_version,form_version) int_version(part)
      part = part + 1
      nchars = 0
      char_version = '    '
   else
      nchars = nchars + 1
      char_version(nchars:nchars) = model_version(i:i)
   endif
enddo
WRITE(*,'(A,A10,4(I3,2X))') 'model_version, version(1:4) = ' &
                            ,model_version,(int_version(i),i=1,4)
   
! assume cam3.0.7 format to start
! test on version cam3.0.3
coord_order = 2
if (int_version(1) < 3) then
   coord_order = 1
else if (int_version(1) == 3 .and. int_version(2) == 0 .and. int_version(3) < 3) then
   coord_order = 1
endif

! pick the order
if (coord_order == 1) then
   coord_3d(1) = num_lons
   coord_3d(2) = num_levs
   coord_3d(3) = num_lats
else if (coord_order == 2) then
   coord_3d(1) = num_lons
   coord_3d(2) = num_lats
   coord_3d(3) = num_levs
endif

coord_2d(1) = num_lons
coord_2d(2) = num_lats

return

end subroutine size_coord


  subroutine nc_read_model_atts(att, att_vals, nflds)
!=======================================================================
! subroutine nc_read_model_atts(att, att_vals, nflds)
!
! reads the value of an attribute for each of the fields in cflds.
!
! should be called with att = one of the attributes from the program variable
! input file, which will be written to the Posterior and Prior.nc files

!----------------------------------------------------------------------
! Local workspace
integer :: i, nchars, ierr
integer :: ncfileid, ncfldid, ncattid

!----------------------------------------------------------------------
integer  :: att_type
integer, intent(in)  :: nflds
character (len=*), intent(in)  :: att 
character (len=128), dimension(nflds), intent(out)  :: att_vals 

! open CAM 'initial' file 
! DONE ALREADY?
call check(nf90_open(path = trim(model_config_file), mode = nf90_nowrite, &
          ncid = ncfileid))

! read CAM 'initial' file attribute desired
PRINT*,'reading ',att
do i = 1,nflds
   call check(nf90_inq_varid(ncfileid, trim(cflds(i)), ncfldid))
! could be inquire_attribute
! 
   ierr = nf90_inquire_attribute(ncfileid, ncfldid, trim(att), att_type, nchars, ncattid) 

   if (ierr == nf90_noerr) then
      att_vals(i)(1:64) = '                                                                '
      att_vals(i)(65:128) = '                                                                '
      call check(nf90_get_att(ncfileid, ncfldid, trim(att) ,att_vals(i) ))
      WRITE(*,'(I6,2A)') ncfldid, cflds(i), trim(att_vals(i))
   else
      WRITE(*,*) ncfldid, cflds(i), 'NOT AVAILABLE'
   end if
enddo

contains 

   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus) 
   integer, intent ( in) :: istatus 
   if (istatus /= nf90_noerr) call error_handler(E_ERR, 'nc_read_model_atts', &
          trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check

end subroutine nc_read_model_atts



  subroutine read_cam_coord(var, idim, cfield)
!=======================================================================
! subroutine read_cam_coord(var, idim, cfield)
!
! should be called with cfield = one of :
!          (/'lat     ','lon     ','gw      '
!           ,'hyai    ','hybi    ','hyam    ','hybm    '/)

!----------------------------------------------------------------------
! Local workspace
integer :: i,ifld             ! grid indices
integer :: ncfileid, ncfldid, idim

!----------------------------------------------------------------------
real(r8), dimension(idim), intent(out) :: var
character (len=8), intent(in)  :: cfield 

! read CAM 'initial' file domain info
call check(nf90_open(path = trim(model_config_file), mode = nf90_nowrite, &
           ncid = ncfileid))

! read CAM 'initial' file field desired
call check(nf90_inq_varid(ncfileid, trim(cfield), ncfldid))
call check(nf90_get_var(ncfileid, ncfldid, var ,start=(/1/) ,count=(/idim/) ))
PRINT*,'reading ',cfield,' using id ',ncfldid
WRITE(*,*) (var(i),i=1,idim)

contains 

   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus) 
   integer, intent ( in) :: istatus 
   if (istatus /= nf90_noerr) call error_handler(E_ERR, 'read_cam_coord', &
          trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check

end subroutine read_cam_coord



  subroutine read_cam_scalar (var, cfield)
!=======================================================================
! subroutine read_cam_scalar (var, cfield)
!
! should be called with cfield = one of :
!          (/'P0      '/)

real(r8), intent(out) :: var
character (len=8), intent(in)  :: cfield 

! Local workspace
integer :: ncfileid, ncfldid

! read CAM 'initial' file domain info
call check(nf90_open(path = trim(model_config_file), mode = nf90_nowrite, &
           ncid = ncfileid))

! read CAM 'initial' file field desired
call check(nf90_inq_varid(ncfileid, trim(cfield), ncfldid))
call check(nf90_get_var(ncfileid, ncfldid, var  ))
PRINT*,'reading ',cfield,' using id ',ncfldid
WRITE(*,*) var

contains 

   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus) 
   integer, intent ( in) :: istatus 
   if (istatus /= nf90_noerr) call error_handler(E_ERR, 'read_cam_scalar', &
          trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check

end subroutine read_cam_scalar



subroutine plevs_cam (ncol    , ncold   ,ps      ,pmid    )
!=======================================================================
!
! Purpose:
! Define the pressures of the interfaces and midpoints from the
! coordinate definitions and the surface pressure.
!
! Method:
!
! Author: B. Boville (plevs0), 
!         Kevin Raeder modified  8/1/03 to use hy[ab][im] from within module
!         rather than in a common block, for use in DART,
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

! coef; commented these out
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use pmgrid

! coef; commented these out
! #include <comhyb.h>

!-----------------------------------------------------------------------
integer , intent(in)  :: ncol               ! Longitude dimension
integer , intent(in)  :: ncold              ! Declared longitude dimension
real(r8), intent(in)  :: ps(ncold)          ! Surface pressure (pascals)
real(r8), intent(out) :: pmid(ncold,num_levs)   ! Pressure at model levels
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
  integer i,k             ! Longitude, level indices
!-----------------------------------------------------------------------
!
! Set midpoint pressures and layer thicknesses
!
! coef
do k=1,num_levs
   do i=1,ncol
      pmid(i,k) = hyam(k)*P0 + hybm(k)*ps(i)
!      pdel(i,k) = pint(i,k+1) - pint(i,k)
   end do
end do

return
end subroutine plevs_cam




  subroutine prog_var_to_vector(var, x)
!=======================================================================
! subroutine prog_var_to_vector(var, x)
!

type(model_type), intent(in) :: var
real(r8), intent(out) :: x(:)

integer :: i, j, k, nf, nt, indx, coord_3d(3), coord_2d(2), coord_order
character(len=129) :: errstring

! Do order as ps, t, u, v, q, tracers to be consistent with b-grid

! HOOK FOR FUTURE 0d and 1d components of state vector
if (state_num_0d .ne. 0 .or. state_num_1d .ne. 0) then
   write(errstring, *) 'scalar and vector components of state vector are not coded ',&
                       'into prog_var_to_vector '
   call error_handler(E_ERR, 'prog_var_to_vector', errstring, source, revision, revdate)
endif

! vars_3d; need this if rearranging coordinates here (?)
! Get order of coordinates to be consistent with read/write_cam_init
! call size_coord(coord_3d, coord_2d, coord_order)

! Start copying fields to straight vector
indx = 0

!  0d variables
do nf = 1, state_num_0d
   indx = indx + 1
   x(indx) = var%vars_0d(nf)
end do

do i = 1, num_lons
!  1d variables
!  FIX for flexible dimension; move out of this loop
   do nf = 1, state_num_1d
      indx = indx + 1
      x(indx) = var%vars_1d(i, nf)
   end do
   do j = 1, num_lats
!  Surface pressure and other 2d 
      do nf = 1, state_num_2d
         indx = indx + 1
         x(indx) = var%vars_2d(i, j, nf)
      end do
!     u,v,t,q, and tracers at successively lower levels
      do k = 1, num_levs
         do nf= 1, state_num_3d
            indx = indx + 1
            x(indx) = var%vars_3d(i, k, j, nf)
         end do
      end do
   end do
end do

! Temporary check
if(indx /= model_size) then
   write(errstring, *) 'indx ',indx,' model_size ',model_size,' must be equal '
   call error_handler(E_ERR, 'prog_var_to_vector', errstring, source, revision, revdate)
endif

end subroutine prog_var_to_vector




  subroutine vector_to_prog_var(x, var) 
!=======================================================================
! subroutine vector_to_prog_var(x, var) 
!

real(r8), intent(in) :: x(:)
type(model_type), intent(out) :: var

integer :: i, j, k, nf, nt, n0, indx
character(len=129) :: errstring

! Start copying fields from straight vector
indx = 0

! HOOK FOR FUTURE 0d and 1d components of state vector
if (state_num_0d .ne. 0 .or. state_num_1d .ne. 0) then
   write(errstring, *) 'scalar and vector components of state vector are not coded ',&
                       'into vector _to_prog_var'
   call error_handler(E_ERR, 'vector_to_prog_var', errstring, source, revision, revdate)
endif


! 0d arrays
do nf = 1, state_num_0d
   indx = indx + 1
   var%vars_0d(nf) = x(indx)
end do

do i = 1, num_lons
!  1d arrays
!  FIX for flexible dimension; move out of this loop?
   do nf = 1, state_num_1d
      indx = indx + 1
      var%vars_1d(i, nf) = x(indx)
   end do
   do j = 1, num_lats
!  Surface pressure and other 2d fields 
      do nf = 1, state_num_2d
         indx = indx + 1
         var%vars_2d(i, j, nf) = x(indx)
      end do
!     u,v,t,q  and tracers at successive levels
      do k = 1, num_levs
         do nf = 1, state_num_3d
            indx = indx + 1
            var%vars_3d(i, k, j, nf) = x(indx)
         end do 
      end do
   end do
end do

! Temporary check
if(indx /= model_size) then
   write(errstring, *) 'indx ',indx,' model_size ',model_size,' must be equal '
   call error_handler(E_ERR, 'vector_to_prog_var', errstring, source, revision, revdate)
endif

end subroutine vector_to_prog_var




  subroutine write_cam_init(file_name, var)
!=======================================================================
! subroutine write_cam_init(file_name, var)

! write CAM 'initial' file fields that have been updated

character (len = *), intent(in) :: file_name
type(model_type), intent(in) :: var

integer i, k, n, m, ifld, ncfileid, ncfldid
integer :: coord_3d(3), coord_2d(2), coord_order
real(r8), allocatable :: temp_3d(:,:,:)

! Read CAM 'initial' file domain info
call check(nf90_open(path = trim(file_name), mode = nf90_write, ncid = ncfileid))
 
! Figure out coordinate sizes based on CAM version
call size_coord(coord_3d, coord_2d, coord_order)
allocate (temp_3d(coord_3d(1),coord_3d(2),coord_3d(3)))

ifld = 0
! 0d fields are first
! kdr; Tim says that nf90_put_var was probably bombing because
!      netcdf recognized that it was writing a scalar, called an f77
!      scalar put_var, and then choked when it saw the count = ...
!      So, this padding is probably not necessary.
do i = 1, state_num_0d
   ifld = ifld + 1
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   call check(nf90_put_var(ncfileid, ncfldid, var%vars_0d(i) ))
end do 

! 1d fields 
do i = 1, state_num_1d
   ifld = ifld + 1
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
! FIX; for flexible dimension
   call check(nf90_put_var(ncfileid, ncfldid, var%vars_1d(:, i), &
      start=(/1, 1/), count = (/coord_3d(1), 1/)))
end do 

! 2d fields 
do i = 1, state_num_2d
   ifld = ifld + 1
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   call check(nf90_put_var(ncfileid, ncfldid, var%vars_2d(:, :, i), &
      start=(/1, 1, 1/), count = (/coord_2d(1), coord_2d(2), 1/)))
end do 

! 3d fields
do i = 1, state_num_3d
   ifld = ifld + 1
!  repackage depending on coord_order then write the dummy variable
!  put it into lon,lev,lat order for continuity with Tim's/historical.
   if (coord_order == 1) then
!     lon,lev,lat as in original CAM
      temp_3d = var%vars_3d(:, :, :, i)
   else if (coord_order == 2) then
!     lon,lat,lev as in new CAM
      do k=1,coord_3d(3)
      do n=1,coord_3d(2)
      do m=1,coord_3d(1)
         temp_3d(m,n,k) = var%vars_3d(m,k,n,i)
      end do
      end do
      end do
   end if
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   call check(nf90_put_var(ncfileid, ncfldid, temp_3d &
             ,start=(/1,1,1,1/) ,count=(/coord_3d(1), coord_3d(2), coord_3d(3), 1/) ))
end do

deallocate (temp_3d)
call check(nf90_close(ncfileid))

contains 
   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus)

   integer, intent ( in) :: istatus

   if (istatus /= nf90_noerr) call error_handler(E_ERR, 'write_cam_init', &
          trim(nf90_strerror(istatus)), source, revision, revdate)

   end subroutine check

end subroutine write_cam_init




  function get_model_size()
!=======================================================================
! function get_model_size()
!

integer :: get_model_size

get_model_size = model_size

end function get_model_size



  subroutine static_init_model()
!=======================================================================
! subroutine static_init_model()
!
! Initializes class data for CAM model (all the stuff that needs to
! be done once. For now, does this by reading info from a fixed
! name netcdf file. 

integer :: i, j, iunit, ierr, io
! calendar types listed in time_manager_mod.f90
integer :: calendar_type = GREGORIAN

! Register the module
call register_module(source, revision, revdate)

! setting calendar type
! this information is NOT passed to CAM; it must be set in the CAM namelist
call set_calendar_type(calendar_type)

! Reading the namelist input
if(file_exist('input.nml')) then
   iunit = open_file('input.nml', action = 'read')
   ierr = 1
   do while(ierr /= 0)
      read(iunit, nml = model_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'model_nml')
   end do
 11 continue
   call close_file(iunit)
! kdr
else
   write(logfileunit, '(A)') 'WARNING; input.nml not available for read of model_nml'
endif

! Record the namelist values 
write(logfileunit, nml=model_nml)

! Get num lons, lats and levs from netcdf and put in global storage
call read_cam_init_size(model_config_file, num_lons, num_lats, num_levs)

! Set the model minimum time step from the namelist seconds and days input
Time_step_atmos = set_time(Time_step_seconds, Time_step_days)
! kdr debug
call print_time(Time_step_atmos)

! Compute overall model size and put in global storage
! FIX for flexible 1D size; switch num_lons and num_lats depending on case
model_size = state_num_0d + num_lons * (state_num_1d + num_lats *  &
            (state_num_2d + num_levs * state_num_3d) )

! Allocate space for longitude and latitude global arrays
! and Allocate space for hybrid vertical coord coef arrays
allocate(lons(num_lons), lats(num_lats), gw(num_lats), hyai(num_levs+1), &
   hybi(num_levs+1), hyam(num_levs), hybm(num_levs))

! values for num_lons and num_lats should come from netcdf file in read_cam_init_size
call read_cam_coord(lons, num_lons, 'lon     ')
call read_cam_coord(lats, num_lats, 'lat     ')
call read_cam_coord(gw  , num_lats, 'gw      ')

! read hybrid vert coord coefs
call read_cam_coord(hyai, num_levs+1, 'hyai    ')
call read_cam_coord(hybi, num_levs+1, 'hybi    ')
call read_cam_coord(hyam, num_levs  , 'hyam    ')
call read_cam_coord(hybm, num_levs  , 'hybm    ')
call read_cam_scalar(P0, 'P0      ')    ! thats a p-zero

write(*, *) 'CAM size initialized as ', model_size

! CAM3 subroutine to order the state vector parts into cflds 
nflds = state_num_0d + state_num_1d + state_num_2d + state_num_3d      
! # fields to read
allocate (cflds(nflds))
call order_state_fields (cflds, nflds)

! GWD; array for the linking of obs_kinds (KIND_) to model field TYPE_s
call obs_field_location(obs_loc_in_sv)

! CAM3 get field attributes needed by nc_write_model_atts from caminput.nc
allocate (state_long_names(nflds), state_units(nflds))    ! , state_units_long_names(nflds))
call nc_read_model_atts('long_name', state_long_names, nflds)
call nc_read_model_atts('units', state_units, nflds)
! call nc_read_model_atts('units_long_name', state_units_long_names, nflds)

end subroutine static_init_model



  subroutine init_model_instance(var)
!=======================================================================
! subroutine init_model_instance(var)
!
! Initializes an instance of a cam model state variable

type(model_type), intent(out) :: var

! Initialize the storage space and return

! FIX for flexible 1d dimension
allocate( &
   var%vars_0d(                              state_num_0d), &
   var%vars_1d(num_lons,                     state_num_1d), &
   var%vars_2d(num_lons,           num_lats, state_num_2d), &
   var%vars_3d(num_lons, num_levs, num_lats, state_num_3d))

end subroutine init_model_instance



  subroutine end_model_instance(var)
!=======================================================================
! subroutine end_model_instance(var)
!
! Ends an instance of a cam model state variable


type(model_type), intent(inout) :: var

deallocate(var%vars_0d, var%vars_1d, var%vars_2d, var%vars_3d)

end subroutine end_model_instance



  subroutine adv_1step(x, Time)
!=======================================================================
! subroutine adv_1step(x, Time)
!
! Does single time-step advance for B-grid model with vector state as
! input and output. This is a modified version of subroutine atmosphere
! in original bgrid_solo_atmosphere driver.

real(r8), intent(inout) :: x(:)

! Time is needed for more general models like this; need to add in to
! low-order models
type(time_type), intent(in) :: Time

! This is a no-op for CAM; only asynch integration
! Can be used to test the assim capabilities with a null advance

end subroutine adv_1step



  subroutine get_state_meta_data(index_in, location, var_type)
!=======================================================================
! subroutine get_state_meta_data(index_in, location, var_type)
!
! Given an integer index into the state vector structure, returns the
! associated location. This is not a function because the more general
! form of the call has a second intent(out) optional argument kind.
! Maybe a functional form should be added?
! Types for this CAM model are, TYPE_PS, TYPE_T, TYPE_U, TYPE_V, TYPE_Q
! see order_state_fields

integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_type

integer  :: indx, num_per_col, col_num, col_elem, col_elem_3dm1, lon_index, lat_index
real(r8) :: lon, lat, lev
integer  :: var_type_temp, which_vert
character(len=129) :: errstring

! HOOK FOR FUTURE 0d and 1d components of state vector
if (state_num_0d .ne. 0 .or. state_num_1d .ne. 0) then
   write(errstring, *) 'scalar and vector components of state vector are not coded ',&
                       'into get_state_meta_data '
   call error_handler(E_ERR, 'get_state_meta_data', errstring, source, revision, revdate)
endif

! Easier to compute with a 0 to size - 1 index
! No it's not (harumph).  Change def of col_num to compensate.
! indx = index_in - 1
indx = index_in

! FIX; where would 1d and 0d variables fit in here?

! Compute number of items per column
num_per_col = num_levs * state_num_3d + state_num_2d

! How many complete columns are before this one
col_num  = (indx - 1) / num_per_col 
col_elem = indx - col_num * num_per_col
col_elem_3dm1 = col_elem - state_num_2d - 1

! What lon and lat index for this column
lon_index = col_num / num_lats
lat_index = col_num - lon_index * num_lats

! Get actual lon lat values from static_init arrays ???
lon = lons(lon_index + 1)
lat = lats(lat_index + 1)

! Since CAM has pure pressure at the top, pure sigma at the
! bottom and hybrid in-between, we are just referring to
! the LEVEL index for the vertical. The coefficients to reconstruct
! pressure, height, etc. are available in the netCDF files.

! same test as 'if (col_elem <= state_num_2d) then', but more efficient
! NOVERT
if (col_elem_3dm1 < 0) then
   lev = -1
   which_vert = which_vert_2d(TYPE_2D(col_elem))
else
   lev = col_elem_3dm1 / state_num_3d  + 1
   var_type_temp = mod(col_elem_3dm1, state_num_3d ) + 1
   which_vert = which_vert_3d(TYPE_3D(var_type_temp))
endif

location = set_location(lon, lat, lev, which_vert)  

! Now figure out which beast in column this is
! If the type is wanted, return it
if(present(var_type)) then
   if (col_elem_3dm1 < 0) then
      ! 2D fields Surface pressure is the first element, 
      var_type = TYPE_2D(col_elem)
      var_type_temp = -99
   else
      ! 3D variables
      var_type_temp = mod(col_elem_3dm1, state_num_3d ) + 1
      var_type = TYPE_3D(var_type_temp)
   endif
end if

end subroutine get_state_meta_data




  subroutine model_interpolate(x, location, type, interp_val, istatus)
!=======================================================================
!

real(r8),           intent(out) :: interp_val
real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: type
integer,            intent(out) :: istatus

integer :: lon_below, lon_above, lat_below, lat_above, i, vstatus, which_vert
real(r8) :: bot_lon, top_lon, delta_lon, bot_lat, top_lat
real(r8) :: lon_fract, lat_fract, val(2, 2), temp_lon, a(2)
real(r8) :: lon, lat, level, lon_lat_lev(3), pressure

! Would it be better to pass state as prog_var_type (model state type) to here?
! As opposed to the stripped state vector. YES. This would give time interp.
! capability here; do we really want this or should it be pushed up?

! Start with no errors in 
istatus = 0
vstatus = 0


! Get the position, determine if it is model level or pressure in vertical
lon_lat_lev = get_location(location)
lon = lon_lat_lev(1); lat = lon_lat_lev(2);
if(vert_is_level(location)) then
   level = lon_lat_lev(3)
else
   pressure = lon_lat_lev(3)
endif

! Get lon and lat grid specs, num_lons, num_lats are globally defined for cam
   bot_lon = lons(1)
   top_lon = lons(num_lons)
   delta_lon = lons(2) - lons(1)
   bot_lat = lats(1)
   top_lat = lats(num_lats)

! Compute bracketing lon indices
if(lon >= bot_lon .and. lon <= top_lon) then
   lon_below = int((lon - bot_lon) / delta_lon) + 1
   lon_above = lon_below + 1
   lon_fract = (lon - ((lon_below - 1) * delta_lon + bot_lon)) / delta_lon
else
! At wraparound point
   lon_below = num_lons
   lon_above = 1
   if(lon < bot_lon) then
      temp_lon = lon + 360.0
   else
      temp_lon = lon
   endif
   lon_fract = (temp_lon - top_lon) / delta_lon
endif


! Next, compute neighboring lat rows
! NEED TO BE VERY CAREFUL ABOUT POLES; WHAT'S BEING DONE MAY BE WRONG
! Inefficient search used for latitudes in Gaussian grid. Might want to speed up.
if(lat >= bot_lat .and. lat <= top_lat) then

   do i = 2, num_lats
      if(lat <= lats(i)) then
         lat_above = i
         lat_below = i - 1
         lat_fract = (lat - lats(lat_below)) / (lats(lat_above) - lats(lat_below))
         goto 20
      end if 
   end do

else if(lat <= bot_lat) then
! South of bottom lat NEED TO DO BETTER: NOT REALLY REGULAR
   lat_below = 1
   lat_above = 2
   lat_fract = 0.0
else
! North of top lat NEED TO DO BETTER: NOT REALLY REGULAR
   lat_below = num_lats - 1
   lat_above = num_lats
   lat_fract = 1.0
endif

20 continue

! Case 1: model level specified in vertical
if(vert_is_level(location)) then
! Now, need to find the values for the four corners
   call get_val(val(1, 1), x, lon_below, lat_below, nint(level), type, vstatus)
   if (vstatus /= 1) call get_val(val(1, 2), x, lon_below, lat_above, nint(level), type, vstatus)
   if (vstatus /= 1) call get_val(val(2, 1), x, lon_above, lat_below, nint(level), type, vstatus)
   if (vstatus /= 1) call get_val(val(2, 2), x, lon_above, lat_above, nint(level), type, vstatus)

else
! Case of pressure specified in vertical
   call get_val_pressure(val(1, 1), x, lon_below, lat_below, pressure, type, vstatus)
   if (vstatus /= 1) call get_val_pressure(val(1, 2), x, lon_below, lat_above, pressure, type, vstatus)
   if (vstatus /= 1) call get_val_pressure(val(2, 1), x, lon_above, lat_below, pressure, type, vstatus)
   if (vstatus /= 1) call get_val_pressure(val(2, 2), x, lon_above, lat_above, pressure, type, vstatus)
endif

! if (pflag > 0) write(53,'(A,2F7.2/)') '  subsurface obs lat, lon = ',lat,lon

! Do the weighted average for interpolation
!write(*, *) 'fracts ', lon_fract, lat_fract
! kdr Guam;
! istatus   meaning                  return expected obs?   assimilate?
! 0         obs and model are fine;  yes                    yes
! 1         fatal problem;           no                     no
! 2         exclude valid obs        yes                    no
!
istatus = vstatus
if(istatus /= 1) then
   do i = 1, 2
      a(i) = lon_fract * val(2, i) + (1.0 - lon_fract) * val(1, i)
   end do
   interp_val = lat_fract * a(2) + (1.0 - lat_fract) * a(1)
else
   interp_val = 0.
endif

end subroutine model_interpolate



subroutine get_val_pressure(val, x, lon_index, lat_index, pressure, obs_kind, istatus)
!=======================================================================
!
! Gets the vertically interpolated value on pressure for variable obs_kind
! at lon_index, lat_index horizontal grid point
!
! This version excludes observations below lowest level pressure and above
! highest level pressure.


real(r8), intent(out) :: val
real(r8), intent(in) :: x(:), pressure
integer, intent(in) :: lon_index, lat_index, obs_kind 
integer, intent(out) :: istatus

real(r8) :: ps(1), pfull(1, num_levs), fraction
type(location_type) :: ps_location
integer :: top_lev, bot_lev, i, k, vstatus
real(r8) :: bot_val, top_val, ps_lon

! No errors to start with
istatus = 0
vstatus = 0

! Need to get the surface pressure at this point. Easy for A-grid.
! kdr debug
!if (lat_index == 6 .and. lon_index == 12 ) &
!   write(*,*) 'get_val_press; calling get_val for surf press'
call get_val(ps(1), x, lon_index, lat_index, -1, KIND_PS, vstatus)
! call get_val(ps(1), x, lon_index, lat_index, -1, 3, vstatus)
if (vstatus > 0) then
   val = 0.
   istatus = 1
   return
endif

! Next, get the values on the levels for this ps
call plevs_cam (1, 1, ps, pfull)
!write(*, *) 'surface pressure is ', ps(1)
!write(*, *) 'pressure levs in model are ', pfull

! Interpolate in vertical to get two bounding levels
!if (obs_kind == 4 .and. lat_index == 6 .and. lon_index == 12 .and. pressure > 90000.) &
!   write(*, *)'    pressure top bottom ', pressure, pfull(1, 1), pfull(1, num_levs)

if(pressure <= pfull(1, 1) .or. pressure >= pfull(1, num_levs)) then
   istatus = 1
   val = 0.

! Exclude obs if obs_kind is found in list from namelist
! (obs_kinds found in obs_kind/obs_kind_mod.f90; KIND_x)
! This should be done in obs_???, but not for now.
! do i=1,exclude_obs_kinds_num
!    if (include .and. obs_kind == exclude_obs_kinds(i)) include = .false.
! enddo
elseif (obs_kind == 3 .or. obs_kind == 5) then
   istatus = 1
   val = 0.
! else if (.not.include) then
!   istatus = 1
!   val = 0.

else 
!   if(pressure < 20000.) then
   if(pressure < highest_obs_pressure_mb * 100.0) then
      istatus = 2
   else
      istatus = 0
   endif
! Search down through pressures
   do i = 2, num_levs 
      if(pressure < pfull(1, i)) then
         top_lev = i -1
         bot_lev = i
         fraction = (pfull(1, i) - pressure) / &
            (pfull(1, i) - pfull(1, i - 1))
         goto 21
      endif
   end do

21 continue

!if (obs_kind == 4 .and. lat_index == 6 .and. lon_index == 12 ) &
!   write(*,*) '              calling get_val for bot val'
   call get_val(bot_val, x, lon_index, lat_index, bot_lev, obs_kind, vstatus)
!if (obs_kind == 4 .and. lat_index == 6 .and. lon_index == 12 ) &
!   write(*,*) '              bot_val,vstatus =',bot_val,vstatus
   if (vstatus == 0) call get_val(top_val, x, lon_index, lat_index, top_lev, obs_kind, vstatus)
   if (vstatus == 0) then
      val = (1.0 - fraction) * bot_val + fraction * top_val
   else
     istatus = 1
      val = 0.
   endif
end if

end subroutine get_val_pressure




  subroutine get_val(val, x, lon_index, lat_index, level, obs_kind, istatus)
!=======================================================================
! function get_val(x, lon_index, lat_index, level, obs_kind, istatus)
!

real(r8), intent(out) :: val
real(r8), intent(in) :: x(:)
integer, intent(in) :: lon_index, lat_index, level, obs_kind
integer, intent(out) :: istatus

integer :: per_col, indx, field_type
character(len=129) :: errstring

! No errors to start with
istatus = 0

! Compute size of grid storage in a column; includes tracers
! Single 2D state vector is pressure

! HOOK FOR FUTURE 0d and 1d components of state vector
if (state_num_0d .ne. 0 .or. state_num_1d .ne. 0) then
   write(errstring, *) 'scalar and vector components of state vector are not coded ',&
                       'into get_val ' 
   istatus = 1 
   call error_handler(E_ERR, 'get_val', errstring, source, revision, revdate)
endif


per_col = state_num_2d + num_levs * state_num_3d

! Find the starting index for this column
! FIX add in 1d and 0d fields here
indx = per_col * (lat_index - 1 + (lon_index - 1) * num_lats)

! Increment index to variable's position within the column

field_type = obs_loc_in_sv(obs_kind)

if (field_type <= 0) then
   istatus = 1
   val = 0.
   return
else if (field_type <= state_num_2d) then
   indx = indx + field_type
else if (field_type <= state_num_2d + state_num_3d) then
   indx = indx + state_num_2d &
               + (level - 1) * state_num_3d &
               + (field_type - state_num_2d)
else
   istatus = 1
   val = 0.
   return
end if

val = x(indx)


end subroutine get_val



  function get_model_time_step()
!=======================================================================
! function get_model_time_step()
!
! Returns the the time step of the model. In the long run should be repalced
! by a more general routine that returns details of a general time-stepping
! capability.

type(time_type) :: get_model_time_step

! Time_step_atmos is global static storage
get_model_time_step =  Time_step_atmos

end function get_model_time_step



  subroutine end_model()
!=======================================================================
! subroutine end_model()
!
! At some point, this stub should coordinate with atmosphere_end but
! that requires an instance variable.

end subroutine end_model



  subroutine init_time(time)
!=======================================================================
! subroutine init_time(time)
!
! For now returns value of Time_init which is set in initialization routines.

type(time_type), intent(out) :: time

! Where should initial time come from here?
! WARNING: CURRENTLY SET TO 0
time = set_time(0, 0)

end subroutine init_time



  subroutine init_conditions(x)
!=======================================================================
! subroutine init_conditions(x)
!
! Reads in restart initial conditions  -- noop for CAM

real(r8), intent(inout) :: x(:)

end subroutine init_conditions


  subroutine model_get_close_states(o_loc, radius, nfound, indices, dist, x)
!=======================================================================
! subroutine model_get_close_states(o_loc, radius, nfound, indices, dist)
!

type(location_type), intent(in)  :: o_loc
real(r8),            intent(in)  :: radius
integer,             intent(out) :: nfound, indices(:)
real(r8),            intent(out) :: dist(:)
real(r8),            intent(in)  :: x(:)

type(location_type) :: s_loc
real(r8) :: loc_array(3), o_lon, o_lat, sloc_array(3), ps(1)
real(r8) :: pfull(1, num_levs), m_press, t_dist
integer  :: num, max_size, i, j, num1
integer  :: hsize, num_per_col, col_base_index, istatus, which_vert
integer,  allocatable :: lon_ind(:), lat_ind(:)
real(r8), allocatable :: close_dist(:)
character(len=129)    :: errstring

! HOOK FOR FUTURE 0d and 1d components of state vector
if (state_num_0d .ne. 0 .or. state_num_1d .ne. 0) then
   write(errstring, *) 'scalar and vector components of state vector are not coded ',&
                       'into model_get_close_states '
   call error_handler(E_ERR, 'model_get_close_states', errstring, source, revision, revdate)
endif

!write(*, *) 'in model_get_close_states', radius
loc_array = get_location(o_loc)
!write(*, *) 'oloc is ', loc_array(:)

! Number found starts at 0
nfound = 0

! Assume that grid size is known from static initialized storage

! Num of close horizontal grid points starts at 0, too
num = 0
! For now, just allocate enough space for all grid points, may want
! to make this smaller at some point for big models.
max_size = num_lons * num_lats
allocate(lon_ind(max_size), lat_ind(max_size), close_dist(max_size))

! Look for close grid points in horizontal only
call grid_close_states2(o_loc, lons, lats, num_lons, num_lats, radius, &
   num, lon_ind, lat_ind, close_dist)
! kdr write(*, *) 'back from grid_close_states num = ', num

! Compute size of grid storage for full levels
hsize = num_lons * num_lats
! FIX; would 0d and 1d fields need to go in here?
num_per_col = num_levs * state_num_3d + state_num_2d

! For vertical localization need the vertical pressure structure for this column
do i = 1, num
   col_base_index = ((lon_ind(i) - 1) * num_lats + lat_ind(i) - 1) * num_per_col

! Compute the ps from the state for this column
! DOES NOT WORK FOR > 1 REGION
!   call get_val(ps(1), x, lon_ind(i), lat_ind(i), -1, 3, istatus)
! TEMP FIX
    ps(1) = 100000.
!   write(*, *) 'ps in model_get_close_states is ', i, ps(1)
! Next get the values on the levels for this ps
   call plevs_cam(1, 1, ps, pfull)

   do j = 1, num_per_col

! Added for vertical localization, 17 May, 2004
      call get_state_meta_data(col_base_index + j, s_loc) 
      sloc_array = get_location(s_loc)
      which_vert = nint(query_location(s_loc))
! Surface pressure has ps as vertical, others have their level's pressure
! Put the appropriate pressure into a location type for computing distance
      if (which_vert == -2) then
         ! NOVERT; field with no vertical location; get_dist will calculate 
         ! horiz dist only based on which_vert of s_loc
      else if(which_vert == -1 ) then       
         ! surface field; change which_vert for the distance calculation
         s_loc = set_location(sloc_array(1), sloc_array(2), ps(1), 2)
      else if(which_vert == 1 ) then
         m_press = pfull(1, int(sloc_array(3)))
         s_loc = set_location(sloc_array(1), sloc_array(2), m_press, 2)
      else
         write(errstring, *) 'which_vert = ',which_vert,' not handled in model_get_close_states '
         call error_handler(E_ERR, 'model_get_close_states', errstring, source, revision, revdate)
      endif

      t_dist = get_dist(s_loc, o_loc)

! kdr
! reduce influence of obs on points obove 150 hPa, which is cap of obs too
! linear with pressure
! ERROR; PS model points have no m_press, so they may be included in this;
!        add which_vert condition
      if (which_vert == 1 .and. m_press < (highest_obs_pressure_mb*100._r8)) then
!         WRITE(*,*) 'model_get_close_states; increasing t_dist from ',t_dist
         t_dist = t_dist / (m_press/(highest_obs_pressure_mb*100._r8))
!         WRITE(*,'(2(A,F10.4),A,3F10.4)') '   to ',t_dist,' for m_press = ',m_press, &
!                    ' and o_loc = ',(loc_array(i),i=1,3)
!         if (t_dist > radius) WRITE(*,*) '   Point ',col_base_index+j,' should be excluded'
      endif
! kdr end

! radius is 2x cutoff 
! These distances are passed through, so recomputation without the tapering is not a problem.
! They're used in cov_cutoff/cov_cutoff_mod.f90.

      if(t_dist < radius) then
         nfound = nfound + 1
         if(nfound <= size(indices)) indices(nfound) = col_base_index + j
         if(nfound <= size(dist)) dist(nfound) = t_dist
      endif
   end do
end do

deallocate(lon_ind, lat_ind, close_dist)

end subroutine model_get_close_states




  subroutine grid_close_states2(o_loc, lons, lats, nlon, nlat, radius, &
                   num, close_lon_ind, close_lat_ind, close_dist)
!=======================================================================
! subroutine grid_close_states2(o_loc, lons, lats, nlon, nlat, radius, &
!                  num, close_lon_ind, close_lat_ind, close_dist)
!
!
! Finds close state points from a particular grid; Just uses horizontal
! distance by setting pressure of state location to same as observation.


type(location_type), intent(in)    :: o_loc
integer,             intent(in)    :: nlon, nlat
real(r8),            intent(in)    :: lons(nlon), lats(nlat), radius
integer,             intent(inout) :: num
integer,             intent(inout) :: close_lon_ind(:), close_lat_ind(:)
real(r8),            intent(out)   :: close_dist(:)

real(r8) :: glat, glon, loc_array(3), o_lon, o_lat, o_lev
real(r8) :: gdist, diff, row_dist(nlon)
integer  :: blat_ind, blon_ind, i, j, lat_ind, lon_ind
integer  :: row_lon_ind(nlon), row_num
real(r8), parameter :: glev = 1.0_r8
type(location_type) :: loc

! Get the lat and lon from the loc
loc_array = get_location(o_loc)
o_lon = loc_array(1)
o_lat = loc_array(2)

! Get index to closest lat and lon for this observation
blat_ind = get_closest_lat_index(o_lat, lats, nlat)
!write(*, *) 'closest latitude in grid is ', blat_ind, lats(blat_ind)
blon_ind = get_closest_lon_index(o_lon, lons, nlon)
!write(*, *) 'closest longitude in grid is ', blon_ind, lons(blon_ind)

! Begin a search along the latitude axis in the positive direction
do lat_ind = blat_ind, nlat
   glat = lats(lat_ind)
   ! Take care of storage round-off
   if(glat < -90.0_r8) glat =  0.0_r8
   if(glat >  90.0_r8) glat = 90.0_r8

   ! Search all the contiguous close longitudes around the base longitude
   call lon_search(glat, glev, blon_ind, o_loc, radius, lons, &
      row_lon_ind, row_dist, row_num)

   ! If none are found, it's time to search in the negative latitude direction
   if(row_num == 0) goto 11

   ! Copy the points found in the row into summary storage
   close_lon_ind(num+1 : num+row_num) = row_lon_ind(1:row_num)
   close_lat_ind(num+1 : num+row_num) = lat_ind
   close_dist(num+1 : num+row_num) = row_dist(1:row_num)
   num = num + row_num
end do

! Search in the negative lat direction
11 continue
do lat_ind = blat_ind - 1, 1, -1
   glat = lats(lat_ind)
   ! Take care of storage round-off
   if(glat < -90.0_r8) glat =  0.0_r8
   if(glat >  90.0_r8) glat = 90.0_r8

   ! Search all the contiguous close longitudes around the base longitude
   call lon_search(glat, glev, blon_ind, o_loc, radius, lons, &
      row_lon_ind, row_dist, row_num)

   ! If none are found, it's time to give up
   if(row_num == 0) return

   ! Copy the points found in the row into summary storage
   close_lon_ind(num+1 : num+row_num) = row_lon_ind(1:row_num)
   close_lat_ind(num+1 : num+row_num) = lat_ind
   close_dist(num+1 : num+row_num) = row_dist(1:row_num)
   num = num + row_num
end do

end subroutine grid_close_states2



  subroutine lon_search(glat, glev, blon_ind, o_loc, radius, lons, &
                      close_lon_ind, close_dist, num)
!=======================================================================
!  subroutine lon_search(glat, glev, blon_ind, o_loc, radius, lons, &
!                     close_lon_ind, close_dist, num)
!
! Given an observation location and radius and a latitude row from a grid,
! searches to find all longitude points in this row that are within radius
! of the observation location and returns their latitude index, longitude
! index, and the distance between them and the observation.

real(r8),            intent(in)  :: glat, glev, radius, lons(:)
integer,             intent(in)  :: blon_ind
type(location_type), intent(in)  :: o_loc
integer,             intent(out) :: close_lon_ind(:), num
real(r8),            intent(out) :: close_dist(:)

type(location_type) :: loc
integer  :: nlon, j, max_pos, lon_ind, which_vert
real(r8) :: glon, gdist, olev

! Total number found is 0 at start
num = 0
nlon = size(lons)

! Search as far as possible in the positive direction
do j = 0, nlon - 1
   max_pos = j
   lon_ind = blon_ind + j
   if(lon_ind > nlon) lon_ind = lon_ind - nlon
   glon = lons(lon_ind)

   ! Correct for longitude storage round-off
   if(glon > 360.0_r8) glon = 360.0_r8
   if(glon <   0.0_r8) glon =   0.0_r8

   ! Use same vertical "philosophy" as the existing location object.
   ! As of April, 2004 -- the vertical is (still) ignored in get_dist.
   which_vert = nint(query_location(o_loc))
   ! Only looking for horizontal distance here, use same level as obs.
   olev = query_location(o_loc, 'VLOC')  
   loc        = set_location(glon, glat, olev, which_vert)
   gdist      = get_dist(loc, o_loc)
   if(gdist <= radius) then
      num = num + 1
      close_lon_ind(num) = lon_ind
      close_dist(num) = gdist
      ! If radius is too far for closest longitude, 
      ! no need to search further or to search other side
   else if (j == 0) then
      return
   else
      ! Look in negative longitude offset direction next
      goto 21
   endif
end do
! Falling off end means the whole longitude circle has been searched; move along
return

! Search around the other way
21 continue
do j = 1, nlon - 1 - max_pos

   lon_ind = blon_ind - j
   if(lon_ind < 1) lon_ind = nlon + lon_ind
   glon = lons(lon_ind)

   ! Correct for longitude storage round-off
   if(glon > 360.0_r8) glon = 360.0_r8
   if(glon <   0.0_r8) glon =   0.0_r8

   ! Use same vertical "philosophy" as the existing location object.
   ! As of April, 2004 -- the vertical is (still) ignored in get_dist.
   which_vert = nint(query_location(o_loc))
   ! Only looking for horizontal distance here, use same level as obs.
   olev = query_location(o_loc, 'VLOC')  
   loc        = set_location(glon, glat, olev, which_vert)
   gdist      = get_dist(loc, o_loc)

   if(gdist <= radius) then
      num = num + 1
      close_lon_ind(num) = lon_ind
      close_dist(num) = gdist
   else
      ! No more longitudes in negative direction
      return
   endif
end do

end subroutine lon_search
                                              


  function nc_write_model_atts( ncFileID ) result (ierr)
!=======================================================================
! function nc_write_model_atts( ncFileID ) result (ierr)
!
! Writes the model-specific attributes to a netCDF file.
! TJH Fri Aug 29 MDT 2003
!
use typeSizes
use netcdf

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

!-----------------------------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: lonDimID, latDimID, ilevDimID, ScalarDimID
integer :: MemberDimID, StateVarDimID, TimeDimID
integer :: lonVarID, latVarID, ilevVarID, hyaiVarID, hybiVarID, P0VarID, gwVarID
integer :: xVarID,StateVarID, StateVarVarID 
integer :: i, ifld
character(len=129)    :: errstring 
character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1,str2

ierr = 0     ! assume normal termination

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file,
! and then put into define mode.
!-------------------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))
call check(nf90_Redef(ncFileID))

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies
!-------------------------------------------------------------------------------

call check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID))
call check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid=  TimeDimID))

if ( TimeDimID /= unlimitedDimId ) then
  write(errstring,*)'Time dimension ID ',TimeDimID,'must match Unlimited Dimension ID ',unlimitedDimId
  call error_handler(E_ERR,'nc_write_model_atts', errstring, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
!-------------------------------------------------------------------------------
call check(nf90_def_dim(ncid=ncFileID, name="StateVariable",  &
                        len=model_size, dimid = StateVarDimID))

!-------------------------------------------------------------------------------
! Write Global Attributes
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date",str1))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate",revdate))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model","CAM"))

! how about namelist input? might be nice to save ...

!-------------------------------------------------------------------------------
! Define the new dimensions IDs
!-------------------------------------------------------------------------------

call check(nf90_def_dim(ncid=ncFileID, name="scalar",   len = 1,   dimid = ScalarDimID))
call check(nf90_def_dim(ncid=ncFileID, name="lat",  len = num_lats,   dimid = latDimID))
call check(nf90_def_dim(ncid=ncFileID, name="lon",  len = num_lons,   dimid = lonDimID))
call check(nf90_def_dim(ncid=ncFileID, name="ilev", len = num_levs+1, dimid = ilevDimID))

!-------------------------------------------------------------------------------
! Create the (empty) Coordinate Variables and their attributes
!-------------------------------------------------------------------------------

! Grid Longitudes
call check(nf90_def_var(ncFileID, name="lon", xtype=nf90_double, &
                                               dimids=lonDimID, varid=lonVarID) )
call check(nf90_put_att(ncFileID, lonVarID, "long_name", "longitude"))
call check(nf90_put_att(ncFileID, lonVarID, "units", "degrees_east"))
call check(nf90_put_att(ncFileID, lonVarID, "valid_range", (/ 0.0_r8, 360.0_r8 /)))

! Grid Latitudes
call check(nf90_def_var(ncFileID, name="lat", xtype=nf90_double, &
                                               dimids=latDimID, varid=latVarID) )
call check(nf90_put_att(ncFileID, latVarID, "long_name", "latitude"))
call check(nf90_put_att(ncFileID, latVarID, "units", "degrees_north"))
call check(nf90_put_att(ncFileID, latVarID, "valid_range", (/ -90.0_r8, 90.0_r8 /)))

! Hybrid grid levels
call check(nf90_def_var(ncFileID, name="ilev", xtype=nf90_double, &
                                                dimids=ilevDimID, varid=ilevVarID) )
call check(nf90_put_att(ncFileID, ilevVarID, "long_name", "hybrid level at interfaces (1000*(A+B))"))
call check(nf90_put_att(ncFileID, ilevVarID, "units", "level"))
call check(nf90_put_att(ncFileID, ilevVarID, "positive", "down"))
call check(nf90_put_att(ncFileID, ilevVarID, "standard_name", "atmosphere_hybrid_sigma_pressure_coordinate"))
call check(nf90_put_att(ncFileID, ilevVarID, "formula_terms", "a: hyai b: hybi P0: P0 ps: PS"))

! Hybrid grid level coefficients, parameters
call check(nf90_def_var(ncFileID, name="hyai", xtype=nf90_double, &
                                                dimids=ilevDimID, varid=hyaiVarID) )
call check(nf90_put_att(ncFileID, hyaiVarID, "long_name", "hybrid A coefficient at layer interfaces" ))

call check(nf90_def_var(ncFileID, name="hybi", xtype=nf90_double, &
                                                dimids=ilevDimID, varid=hybiVarID) )
call check(nf90_put_att(ncFileID, hybiVarID, "long_name", "hybrid B coefficient at layer interfaces"))
call check(nf90_def_var(ncFileID, name="P0", xtype=nf90_double, &
                                                dimids=ScalarDimID, varid=P0VarID) )
call check(nf90_put_att(ncFileID, P0VarID, "long_name", "reference pressure"))
call check(nf90_put_att(ncFileID, P0VarID, "units", "Pa"))

! Gaussian weights -- because they're there.
call check(nf90_def_var(ncFileID, name="gw", xtype=nf90_double, &
                                                dimids=latDimID, varid=gwVarID) )
call check(nf90_put_att(ncFileID, gwVarID, "long_name", "gauss weights"))


if ( output_state_vector ) then

! CAM3; need to adapt this to state_long_names, state_units, etc?


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
              dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), varid=StateVarID))
   call check(nf90_put_att(ncFileID, StateVarID, "long_name", "model state or fcopy"))
   call check(nf90_put_att(ncFileID, StateVarId, "vector_to_prog_var","CAM"))

   ! Leave define mode so we can fill 
   call check(nf90_enddef(ncfileID))

   ! Fill the state variable coordinate variable
   call check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ))

else

!  write(*,*)'ERROR:CAM:model_mod:    trying to output the prognostic variables.'
!  write(*,*)'      That is not implemented yet.'
!  write(*,*)'      TJH 27 June 2003'
!  stop

   !----------------------------------------------------------------------------
   ! TJH;
   ! We need to process the prognostic variables.
   ! I like the CAM philosophy of using nflds to declare the number of prognostic
   ! variables and an array of characters to specify them. This is clearly
   ! the way it must be for models with "lots"/variable numbers of params. 
   ! 
   ! I'd like to see the metadata handled the same way.
   !
   !----------------------------------------------------------------------------
   ! repeated for reference
   !
   !character (len=8),dimension(nflds) :: cflds = &
   !       (/'PS      ','T       ','U       ','V       ','Q       ' /)
   !
   ! cflds(1) ... PS ... (lon,     lat,time)
   ! cflds(2) ... T  ... (lon,ilev,lat,time)
   ! cflds(3) ... U  ... (lon,ilev,lat,time)
   ! cflds(4) ... V  ... (lon,ilev,lat,time)
   ! cflds(5) ... Q  ... (lon,ilev,lat,time)

   !----------------------------------------------------------------------------
   ! Create the (empty) Variables and the Attributes
   !----------------------------------------------------------------------------

   ! 0-d fields
   ifld = 0
   do i = 1,state_num_0d
      ifld = ifld + 1
      call check(nf90_def_var(ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real, &
                 dimids = (/ MemberDimID, unlimitedDimID /), &
                 varid  = xVarID))
      call check(nf90_put_att(ncFileID, xVarID, "long_name", state_long_names(ifld)))
      call check(nf90_put_att(ncFileID, xVarID, "units", state_units(ifld)))
!       call check(nf90_put_att(ncFileID, xVarID, "units_long_name", state_units_long_names(ifld)))
   enddo

   ! 1-d fields
! FIX; handle flexible 1d variables
! kdr; how to handle 1D dimension?
!                 dimids = (/ lonDimID, MemberDimID, unlimitedDimID /), &
!                             ^^^
   do i = 1,state_num_1d
      ifld = ifld + 1
      call check(nf90_def_var(ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real, &
                 dimids = (/ lonDimID, MemberDimID, unlimitedDimID /), &
                 varid  = xVarID))
      call check(nf90_put_att(ncFileID, xVarID, "long_name", state_long_names(ifld)))
      call check(nf90_put_att(ncFileID, xVarID, "units", state_units(ifld)))
!       call check(nf90_put_att(ncFileID, xVarID, "units_long_name", state_units_long_names(ifld)))
   enddo

   ! 2-d fields
   do i = 1,state_num_2d
      ifld = ifld + 1
      call check(nf90_def_var(ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real, &
                 dimids = (/ lonDimID, latDimID, MemberDimID, unlimitedDimID /), &
                 varid  = xVarID))
      call check(nf90_put_att(ncFileID, xVarID, "long_name", state_long_names(ifld)))
      call check(nf90_put_att(ncFileID, xVarID, "units", state_units(ifld)))
!       call check(nf90_put_att(ncFileID, xVarID, "units_long_name", state_units_long_names(ifld)))
   enddo

   ! 3-d fields
   do i = 1,state_num_3d
      ifld = ifld + 1
      call check(nf90_def_var(ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real, &
            dimids = (/ lonDimID, ilevDimID, latDimID, MemberDimID, unlimitedDimID /), &
            varid  = xVarID))
      call check(nf90_put_att(ncFileID, xVarID, "long_name", state_long_names(ifld)))
      call check(nf90_put_att(ncFileID, xVarID, "units", state_units(ifld)))
!       call check(nf90_put_att(ncFileID, xVarID, "units_long_name", state_units_long_names(ifld)))
   enddo

   ! Leave define mode so we can fill 
   call check(nf90_enddef(ncfileID))

endif

!-------------------------------------------------------------------------------
! Fill the coordinate variables
!-------------------------------------------------------------------------------

call check(nf90_put_var(ncFileID,  ilevVarID, (/ (i,i=1, num_levs+1) /) ))
call check(nf90_put_var(ncFileID,   latVarID, lats ))
call check(nf90_put_var(ncFileID,   lonVarID, lons ))
call check(nf90_put_var(ncFileID,  hyaiVarID, hyai ))
call check(nf90_put_var(ncFileID,  hybiVarID, hybi ))
call check(nf90_put_var(ncFileID,    gwVarID,   gw ))
call check(nf90_put_var(ncFileID,    P0VarID,   P0 ))

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call check(nf90_sync(ncFileID))

write (*,*)'nc_write_model_atts: netCDF file ',ncFileID,' is synched ...' 

contains
   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus)
   integer, intent ( in) :: istatus
   if (istatus /= nf90_noerr) call error_handler(E_ERR, 'nc_write_model_atts', &
          trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check

end function nc_write_model_atts



  function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)
!=======================================================================
! function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)
!
! Writes the model-specific variables to a netCDF file
! TJH 25 June 2003
!
! There are two different (staggered) 3D grids being used simultaneously here. 
! The routines "prog_var_to_vector" and "vector_to_prog_var", 
! packs the prognostic variables into
! the requisite array for the data assimilation routines. That routine
! is the basis for the information stored in the netCDF files.
!
! TemperatureGrid : surface pressure  vars%ps(tis:tie, tjs:tje) 
!                 : temperature       vars%t (tis:tie, tjs:tje, klb:kup)
!                 : tracers           vars%r (tis:tie, tjs:tje, klb:kub, 1:vars%ntrace)
! VelocityGrid    : u                 vars%u (vis:vie, vjs:vje, klb:kub) 
!                 : v                 vars%v (vis:vie, vjs:tje, klb:kup)
!
! So there are six different dimensions and five different variables as long as
! simply lump "tracers" into one. 

use typeSizes
use netcdf

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

!-----------------------------------------------------------------------------------------
! real, dimension(SIZE(statevec)) :: x     CAM is r8, no need to make a r4 copy ...
type(model_type) :: Var

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID, ncfldid
integer :: ifld, i

ierr = 0     ! assume normal termination

!-------------------------------------------------------------------------------
! CAM storage bounds are 'tight' -- no indices needed
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! then get all the Variable ID's we need.
!-------------------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

if ( output_state_vector ) then

   call check(NF90_inq_varid(ncFileID, "state", StateVarID) )
   call check(NF90_put_var(ncFileID, StateVarID, statevec,  &
                start=(/ 1, copyindex, timeindex /)))                               

else

   !----------------------------------------------------------------------------
   ! Fill the variables
   !----------------------------------------------------------------------------

   call init_model_instance(Var)     ! Explicity released at end of routine. 

   call vector_to_prog_var(statevec,  Var)
   
   ifld = 0
   ZeroDVars : do i = 1, state_num_0d    
      ifld = ifld + 1
      call check(NF90_inq_varid(ncFileID, trim(cflds(ifld)), ncfldid))
      call check(nf90_put_var( ncFileID, ncfldid, Var%vars_0d(i), &
                 start=(/ copyindex, timeindex /) ))
!, &
!                 count=(/1, 1/) ))
   end do ZeroDVars

! FIX flexible 1d variable size
   OneDVars : do i = 1, state_num_1d    
      ifld = ifld + 1
      call check(NF90_inq_varid(ncFileID, trim(cflds(ifld)), ncfldid))
      call check(nf90_put_var( ncFileID, ncfldid, Var%vars_1d(:,i), &
                 start=(/ 1, copyindex, timeindex /), &
                 count=(/num_lons, 1, 1/) ))
   end do OneDVars

   TwoDVars : do i = 1, state_num_2d    
      ifld = ifld + 1
      call check(NF90_inq_varid(ncFileID, trim(cflds(ifld)), ncfldid))
      call check(nf90_put_var( ncFileID, ncfldid, Var%vars_2d(:,:,i), &
                 start=(/ 1, 1, copyindex, timeindex /), &
                 count=(/num_lons, num_lats, 1, 1/) ))
   end do TwoDVars

   ThreeDVars : do i = 1,state_num_3d
      ifld = ifld + 1
      call check(NF90_inq_varid(ncFileID, trim(cflds(ifld)), ncfldid))
! changing to lon,lat,lev might cause P{oste}rior_Diag.nc files to be different
      call check(nf90_put_var( ncFileID, ncfldid, Var%vars_3d(:,:,:,i), &
                 start=(/ 1,1,1,copyindex,timeindex /), &
                 count=(/num_lons,num_levs,num_lats,1,1/) ))
   end do ThreeDVars

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

write (*,*)'Finished filling variables ...'
call check(nf90_sync(ncFileID))
write (*,*)'netCDF file is synched ...'

call end_model_instance(Var)   ! should avoid any memory leaking

contains
   ! Internal subroutine - checks error status after each netcdf, prints 
   !                       text message each time an error code is returned. 
   subroutine check(istatus)
   integer, intent ( in) :: istatus
   if (istatus /= nf90_noerr) call error_handler(E_ERR, 'nc_write_model_vars', &
          trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check

end function nc_write_model_vars





  function get_closest_lat_index(o_lat, lats, nlat)
!=======================================================================
! function get_closest_lat_index(o_lat, lats, nlat)
!

integer, intent(in) :: nlat
real(r8), intent(in) :: o_lat, lats(nlat)
integer :: get_closest_lat_index

real(r8) :: lat_bot, lat_top, lat_int, diff
integer :: lower_ind

! Find closest lat
lat_bot = lats(1)
lat_top = lats(nlat)
lat_int = lats(2) - lats(1)
if(o_lat <= lat_bot) then
   get_closest_lat_index = 1
else if(o_lat >= lat_top) then
   get_closest_lat_index = nlat
else
   diff = (o_lat - lat_bot) / lat_int
   lower_ind = int(diff) + 1
   if(diff - int(diff) < 0.5) then
      get_closest_lat_index = lower_ind
   else
      get_closest_lat_index = lower_ind + 1
   endif
endif

end function get_closest_lat_index



  function get_closest_lon_index(o_lon, lons, nlon)
!=======================================================================
! function get_closest_lon_index(o_lon, lons, nlon)

integer, intent(in) :: nlon
real(r8), intent(in) :: o_lon, lons(nlon)
integer :: get_closest_lon_index

real(r8) :: diff, lon_bot, lon_top, lon_int
integer :: lower_ind, blon_ind

! Find closest longitude on grid to given longitude
lon_bot = lons(1)
lon_top = lons(nlon)
lon_int = lons(2) - lons(1)
if(o_lon <= lon_bot) then
   diff = (lon_bot - o_lon) / lon_int
   if(diff > 0.5) then
      get_closest_lon_index = nlon
   else
      get_closest_lon_index = 1
   end if
else if(o_lon >= lon_top) then
   diff = (o_lon - lon_top) / lon_int
   if(diff > 0.5) then
      get_closest_lon_index = 1
   else
      get_closest_lon_index = nlon
   end if
else
   diff = (o_lon - lon_bot) / lon_int
   lower_ind = int(diff) + 1
   if(diff - int(diff) < 0.5) then
      get_closest_lon_index = lower_ind
   else
      get_closest_lon_index = lower_ind + 1
   end if
end if

end function get_closest_lon_index


  subroutine pert_model_state(state, pert_state, interf_provided)
!       call pert_model_state(ens_mean, temp_ens, interf_provided)
!=======================================================================
! subroutine pert_model_state(state, pert_state, interf_provided)
!
! Perturbs a model state for generating initial ensembles
! Returning interf_provided means go ahead and do this with uniform
! small independent perturbations.

! added to give each ens member a different sequence when perturbing model parameter fields

real(r8), intent(in)    :: state(:)
real(r8), intent(out)   :: pert_state(:)
logical,  intent(out)   :: interf_provided

type(random_seq_type)   :: random_seq
type(model_type)        :: var_temp
integer                 :: i, ipert, j, k, m, field_num, ens_member, iunit
real(r8)                :: pert_val

! FIX for 1D 0D  fields?

! trans_pv_sv_pert0.f90 needs to perturb model parameters for the filter_ics.
! Use the (single) state value as the "ens_mean" here.

interf_provided = .true.

call init_model_instance(var_temp)
call vector_to_prog_var(state,var_temp)

! If first call initialize random sequence for perturbations.
if(first_pert_call) then 
   call init_random_seq(random_seq)
   first_pert_call = .false.
endif

! init_random_seq calls init_ran1, but I need to call init_ran1 with a different seed/temp 
! for each ens_member.  init_random_seq only needed for documentation; initializing 
! the module
if(file_exist('ens_member')) then
   iunit = open_file('ens_member', action = 'read')
   read(iunit, *) ens_member
   call init_ran1(random_seq,-1*ens_member)
   call close_file(iunit)
else
   WRITE(*,*) 'no ens_member file available; perturbing next ensemble member (in filter) '
endif

ipert = 1
do while (state_names_pert(ipert) /= '        ')
   WRITE(*,*) 'Perturbing ',state_names_pert(ipert)
   do m=1,nflds
      if (state_names_pert(ipert) == cflds(m)) then
         WRITE(*,*) '   Found match  ',cflds(m)
           
         if (m <= state_num_2d + state_num_1d + state_num_0d) then
            field_num = m - state_num_1d - state_num_0d
            WRITE(*,'(A,1p2E12.4)') '    first and last state = ', &
                   var_temp%vars_2d(1,1,field_num),var_temp%vars_2d(num_lons,num_lats,field_num)
            do j = 1, num_lats
            do i = 1, num_lons
               pert_val = random_gaussian(random_seq, var_temp%vars_2d(i,j,field_num), &
                                          state_names_sd(ipert)) 
               var_temp%vars_2d(i,j,field_num) = pert_val
            end do
            end do
            WRITE(*,'(A,1p2E12.4)') ' new first and last state = ', &
                   var_temp%vars_2d(1,1,field_num),var_temp%vars_2d(num_lons,num_lats,field_num)
         else  
            field_num = m - state_num_2d - state_num_1d - state_num_0d
            do j = 1, num_lats
            do k = 1, num_levs
            do i = 1, num_lons
!              pert_val = rand#(O(0-1)) * standard dev  + mean
               pert_val = random_gaussian(random_seq, var_temp%vars_3d(i,k,j,field_num), &
                                          state_names_sd(ipert)) 
               var_temp%vars_3d(i,k,j,field_num) = pert_val
            end do
            end do
            end do
         endif
      end if
   end do
   ipert = ipert + 1
end do
call prog_var_to_vector(var_temp,pert_state)
call end_model_instance(var_temp)

end subroutine pert_model_state



  subroutine order_state_fields(cflds,nflds)
!=======================================================================
! subroutine order_state_fields(cflds,nflds)
! 
! fills cflds with state_names for use in I/O of caminput.nc
! Could eventually tally the number of each kind of field; 2D,3D
! and compare each entry against a master list.
! Sort by class of variable too? So user could provide one unordered list?
! Also assigns TYPE_s for use by get_state_meta_data, and other routines

integer :: i, nfld
integer, intent(in) :: nflds
character (len = *), dimension(nflds), intent(out) :: cflds 
character (len = 129) :: errstring

allocate (TYPE_1D(state_num_1d),TYPE_2D(state_num_2d),TYPE_3D(state_num_3d))
! kdr where should these be deallocated?


nfld = 0

! 0D fields
do i=1,state_num_0d
   nfld = nfld + 1
   cflds(nfld)(:) = state_names_0d(i)
end do

! 1D fields
do i=1,state_num_1d
   nfld = nfld + 1
   cflds(nfld)(:) = state_names_1d(i)
   TYPE_1D(i) = nfld
end do

! 2D fields
do i=1,state_num_2d
   nfld = nfld + 1
   cflds(nfld)(:) = state_names_2d(i)
   TYPE_2D(i) = nfld
   if (state_names_2d(i) == 'PS      ') TYPE_PS = nfld
end do

! 3D fields (including q)
do i=1,state_num_3d
   nfld = nfld + 1
   cflds(nfld)(:) = state_names_3d(i)
   TYPE_3D(i) = nfld
   if (state_names_3d(i) == 'T       ') TYPE_T = nfld
   if (state_names_3d(i) == 'U       ') TYPE_U = nfld
   if (state_names_3d(i) == 'V       ') TYPE_V = nfld
   if (state_names_3d(i) == 'Q       ') TYPE_Q = nfld
end do

if (nfld .ne. nflds) then
   write(errstring, *) 'nfld = ',nfld,', nflds = ',nflds,' must be equal '
   call error_handler(E_ERR, 'order_state_fields', errstring, source, revision, revdate)
else
   write(logfileunit,'(/A/)') 'State vector is composed of '
!   write(logfileunit,'((8(A8,1X)))') (cflds(i),i=1,nflds)
   do i=1,state_num_2d
      write(logfileunit,'(/A,I4)') cflds(i), TYPE_2D(i)
   enddo
   do i=1,state_num_3d
      write(logfileunit,'(/A,I4)') cflds(state_num_2d+i), TYPE_3D(i)
   enddo
   write(logfileunit,'(/A)') 'TYPE_PS, TYPE_T, TYPE_U, TYPE_V, TYPE_Q = ' 
   write(logfileunit,'((8(I8,1X)))') TYPE_PS, TYPE_T, TYPE_U, TYPE_V, TYPE_Q
endif

return

end subroutine order_state_fields


  subroutine obs_field_location(obs_loc_in_sv)
!=======================================================================
! subroutine obs_field_location(obs_loc_in_sv)

! ? Should this be a function instead; removes need to dimension obs_loc_in arbitrarily
!   and wastefully.  But then it's called millions of times, instead of accessing an
!   array that's defined once.

! Makes an array of 'locations w/i the state vector'
! of  all the available obs kinds that come from obs_kind_mod. 
! The obs kind that's needed will be the index into this array,
! the corresponding value will be the position of that field (not individual variable) 
! within the state vector according to state_name_Xd.  
! There will be lots of empty array elements, since KIND_x has a lot of "missing" values.
! This subroutine will be called from static_init_model, so it will not have to be 
! recomputed for every obs.

! use     obs_kind_mod, only : KIND_U, KIND_V, KIND_PS, KIND_T, KIND_QV, KIND_P

integer :: i, nfld
integer, intent(out) :: obs_loc_in_sv(:)
character (len = 129) :: errstring

! 2D fields
obs_loc_in_sv(KIND_PS) = TYPE_PS

! 3D fields
obs_loc_in_sv(KIND_T) = TYPE_T
obs_loc_in_sv(KIND_U) = TYPE_U
obs_loc_in_sv(KIND_V) = TYPE_V
obs_loc_in_sv(KIND_QV) = TYPE_Q

write(*,*) 'OBS_KIND   FIELD_TYPE'
do i=1,1000
   if (obs_loc_in_sv(i) /= 0) write(*,'(2I8)') i, obs_loc_in_sv(i)
enddo

! In the future, if fields are not ordered nicely, or if users are specifying
! correspondence of obs fields with state fields, I may want code like:
! The max size of KIND_ should come from obs_kind_mod
! do i=1,state_num_3d
!    if (state_name_3d(i)(1:1) == 'T' .and. &
!        KIND_T <= 1000) ) obs_loc_in_sv(KIND_T) = TYPE_3D(i)
! enddo 

return

end subroutine obs_field_location

!#######################################################################
! end of cam model_mod
!#######################################################################

end module model_mod
