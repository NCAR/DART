! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module model_mod

! <next three lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
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

implicit none
private

public model_type, prog_var_to_vector, vector_to_prog_var, read_cam_init, &
   read_cam_init_size, init_model_instance, end_model_instance, &
   write_cam_init, get_model_size, static_init_model, &
   get_state_meta_data, get_model_time_step, model_interpolate, &
   init_conditions, init_time, adv_1step, end_model, &
   model_get_close_states, nc_write_model_atts, nc_write_model_vars, &
   TYPE_PS, TYPE_T, TYPE_U, TYPE_V, TYPE_Q, TYPE_TRACER, pert_model_state

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
integer, parameter :: TYPE_PS = 0, TYPE_T = 1, TYPE_U = 2, TYPE_V = 3, TYPE_Q = 4, TYPE_TRACER = 5

!----------------------------------------------------------------------

! A type for cam model, very simple for now for conversion only
type model_type
    private
   real(r8), pointer :: vars_2d(:, :, :)
   real(r8), pointer :: vars_3d(:, :, :, :)
   real(r8), pointer :: tracers(:, :, :, :)
end type model_type

!-----------------------------------------------------------------------                       
! need a model namelist
! output_state_vector = .true.     results in a "state-vector" netCDF file
! output_state_vector = .false.    results in a "prognostic-var" netCDF file
                                                                                               
logical :: output_state_vector = .false.
                                                                                               
namelist /model_nml/ output_state_vector

!----------------------------------------------------------------------
! File where basic info about model configuration can be found; should be namelist

character(len = 128) :: model_config_file = 'caminput.nc'
!character(len = 128) :: model_config_file = 'T5H0-12icl.cam2.i.0001-09-01-43200.nc'
!----------------------------------------------------------------------

!
! Global storage for describing cam model class
integer :: model_size, num_lons, num_lats, num_levs

type(time_type) :: Time_step_atmos

integer, parameter :: n3dflds=4     ! # of 3d fields to read from file
!                                   including Q, but not tracers
integer, parameter :: pcnst  =0     ! advected tracers (don't include Q in this)
integer, parameter :: pnats  =0     ! nonadvected tracers
integer, parameter :: num_tracers = pcnst + pnats
integer, parameter :: n2dflds=1     ! # of 2d fields to read from file
! derived parameters
integer, parameter :: n3tflds  = n3dflds+pcnst+pnats   ! # fields to read
integer, parameter :: nflds  = n3tflds+n2dflds         ! # fields to read

! Arrays to store lats, lons, and gaussian weights
real(r8), allocatable :: lons(:), lats(:), gw(:)

!coef hybrid coeffs at interfaces and mid-levels
real(r8), allocatable :: hyai(:), hybi(:), hyam(:), hybm(:)

!reference pressure for pressure calculations using hybrid coeff 
! should be read from same file as hybrid coeffs?
real(r8):: P0               ! reference pressure

! list variables according to the category to which they belong, 
! in the order the categories appear above (n3dflds,pcnst,pnats,n2dflds).
! from ncdump of CAM standard configuration initial file:
!   lat = 64 ; lon = 128 ; lev = 26 ;
!   hyai, hybi, hyam, hybm, gw
!   U, V, T, Q, PS, 
!   (names of advected and nonadv constituents)
!   PHIS, SGH, LANDM, 
!   TS, TSICE, TS1, TS2, TS3, TS4
!   SNOWHICE, LANDFRAC, CWAT
character (len=8),dimension(nflds) :: cflds = &
          (/'PS      ','T       ','U       ','V       ','Q       ' /)
!          (/'U       ','V       ','T       ','Q       ','PS      ' /)

!  ??? Need to understand what Kevin was doing with nlevs; I may be missing something essential
integer, dimension(nflds) :: nlevs 
data nlevs/ n3tflds*0, n2dflds*1/

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

write(*, *) 'file_name in read_cam_init is ', trim(file_name)

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
integer :: i,j,ifld  ! grid and constituent indices
integer :: plon, plat, plev

character (len=NF90_MAX_NAME) :: clon,clat,clev
integer :: londimid, levdimid, latdimid, ncfileid, ncfldid

!----------------------------------------------------------------------
! read CAM 'initial' file domain info
call check(nf90_open(path = trim(file_name), mode = nf90_write, ncid = ncfileid))

! Could do this for storage size error check later
!call check(nf90_inq_dimid(ncfileid, 'lon', londimid))
!call check(nf90_inq_dimid(ncfileid, 'lat', latdimid))
!call check(nf90_inq_dimid(ncfileid, 'lev', levdimid))

! get dimension sizes
!call check(nf90_inquire_dimension(ncfileid, londimid, clon , plon ))
!call check(nf90_inquire_dimension(ncfileid, latdimid, clat , plat ))
!call check(nf90_inquire_dimension(ncfileid, levdimid, clev , plev ))

! Get the Gaussian Weights
! call check(nf90_inq_varid(ncfileid, 'gw', ncfldid))
! call check(nf90_get_var(ncfileid, ncfldid, gw))

!!! ??? I'm not sure what this next part is really doing, verify
! specify # vertical levels for 3D fields (2d have already been filled)
!i=1
!do while (nlevs(i)==0 .and. i<=nflds)
!   nlevs(i) = plev
!   i=i+1
!end do
!WRITE(*,*) 'nlevs = ',(nlevs(j),j=1,nflds)

! read CAM 'initial' file fields desired

!2d fields; hard coded for only 1 2d field here (generalize???)
do ifld= 1, 1
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
!  fields on file are 3D; lon, lat, TIME(=1)
   call check(nf90_get_var(ncfileid, ncfldid, var%vars_2d(:, :, 1) &
             ,start=(/1,1,1/) ,count=(/num_lons, num_lats, 1/) ))
end do

! 3d fields; hard coded for only 4 here now (generalize???)
do ifld=2, 5
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
!  fields on file are 4D; lon, lev, lat, TIME(=1) 
   call check(nf90_get_var(ncfileid, ncfldid, var%vars_3d(:, :, :, ifld - 1) &
             ,start=(/1,1,1,1/) ,count=(/num_lons, num_levs, num_lats,1/) ))

!!! ?? WARNING: DOES THE NUMBER OF VERTICAL LEVELS PER 3D FIELD VARY
! AS KEVIN"S CODE SUGGESTED. IF SO NEED TO FIX
end do

contains 

   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus) 
   integer, intent ( in) :: istatus 
   if (istatus /= nf90_noerr) call error_handler(E_ERR, 'read_cam_init', &
          trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check

end subroutine read_cam_init



  subroutine read_cam_coord(var, idim, cfield)
!=======================================================================
! subroutine read_cam_coord(var, idim, cfield)
!
! should be called with cfield = one of :
!          (/'lat     ','lon     ','gw      ','P0      '
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
!integer , intent(in)  :: nver               ! vertical dimension
! coef
!real(r8), intent(in) :: hyai(nver+1)        ! hybrid As at interface levels
!real(r8), intent(in) :: hybi(nver+1)        ! hybrid Bs at interface levels
!real(r8), intent(in) :: hyam(nver)          ! hybrid As at layer mid points
!real(r8), intent(in) :: hybm(nver)          ! hybrid Bs at layer mid points
!real(r8), parameter  :: P0 = 1.0E+05       ! reference pressure
! end coef
real(r8), intent(in)  :: ps(ncold)          ! Surface pressure (pascals)
!real(r8), intent(out) :: pint(ncold,num_levs+1) ! Pressure at model interfaces
real(r8), intent(out) :: pmid(ncold,num_levs)   ! Pressure at model levels
!real(r8), intent(out) :: pdel(ncold,num_levs)   ! Layer thickness (pint(k+1) - pint(k))
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
  integer i,k             ! Longitude, level indices
!-----------------------------------------------------------------------
!
! Set interface pressures
!
!do k=1,num_levs+1
!   do i=1,ncol
!      pint(i,k) = hyai(k)*P0 + hybi(k)*ps(i)
!   end do
!end do
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

integer :: i, j, k, nf, nt, indx
character(len=129) :: errstring

! Do order as ps, t, u, v, q, tracers to be consistent with b-grid

! Start copying fields to straight vector
indx = 0
do i = 1, num_lons
   do j = 1, num_lats
!  Surface pressure and other 2d flds are first
      do nf = 1, n2dflds
         indx = indx + 1
         x(indx) = var%vars_2d(i, j, nf)
      end do
!     u,v,t,q, and tracers at successively lower levels
      do k = 1, num_levs
         do nf= 1, n3dflds
            indx = indx + 1
            x(indx) = var%vars_3d(i, k, j, nf)
         end do
         do nt = 1, num_tracers
            IF (i==1 .and. j==1) PRINT*,'filling tracers'
            indx = indx + 1
            x(indx) = var%tracers(i, k, j, nt)
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
do i = 1, num_lons
   do j = 1, num_lats
! Surface pressure and other 2d fields are first
      do nf = 1, n2dflds
         indx = indx + 1
         var%vars_2d(i, j, nf) = x(indx)
      end do
!     u,v,t,q  and tracers at successive levels
      do k = 1, num_levs
         do nf = 1, n3dflds
            indx = indx + 1
            var%vars_3d(i, k, j, nf) = x(indx)
         end do 
         do nt = 1, num_tracers
            indx = indx + 1
            var%tracers(i, k, j, nt) = x(indx)
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


character (len = *), intent(in) :: file_name
type(model_type), intent(in) :: var

integer ifld, ncfileid, ncfldid

! Read CAM 'initial' file domain info
call check(nf90_open(path = trim(file_name), mode = nf90_write, ncid = ncfileid))

! Try doing this in revised order
! 2d fields are first
do ifld = 1, 1
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   call check(nf90_put_var(ncfileid, ncfldid, var%vars_2d(:, :, ifld), &
      start=(/1, 1, 1/), count = (/num_lons, num_lats, 1/)))
end do 

! write CAM 'initial' file fields that have been updated
! ????? WARNING: PROBLEMS WITH TRACERS HERE
! ??? THIS WHOLE MODULE NEEDS FURTHER REVISION WITH KEVIN"S HELP
! Now do 3d fields, ignoring tracers for now
do ifld= 2, 5
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   call check(nf90_put_var(ncfileid, ncfldid, var%vars_3d(:,:,:,ifld - 1) &
             ,start=(/1,1,1,1/) ,count=(/num_lons, num_levs, num_lats,1/) ))
end do

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
! name netcdf file. Need to make this file a namelist parameter
! at some point.

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
   enddo
 11 continue
   call close_file(iunit)
endif

! Record the namelist values 
write(logfileunit, nml=model_nml)

! Get num lons, lats and levs from netcdf and put in global storage
call read_cam_init_size(model_config_file, num_lons, num_lats, num_levs)

! Compute the model size (Need to be sure where all these come from
! and if they themselves need to be initialized here)

! Need to set Time_step_atmos here to 1 hour for now 
! kdr automate this? (setting model time"step")
! Time_step_atmos = set_time(43200, 0)
Time_step_atmos = set_time(3600, 0)
! kdr debug
call print_time(Time_step_atmos)


! Compute overall model size and put in global storage
model_size = num_lons * num_lats * (n2dflds + num_levs * &
   (n3dflds + pcnst + pnats))

! Allocate space for longitude and latitude global arrays
allocate(lons(num_lons), lats(num_lats), gw(num_lats))
! Allocate space for hybrid vertical coord coef arrays
allocate(hyai(num_levs+1), hybi(num_levs+1), hyam(num_levs), hybm(num_levs))

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

end subroutine static_init_model



  subroutine init_model_instance(var)
!=======================================================================
! subroutine init_model_instance(var)
!
! Initializes an instance of a cam model state variable

type(model_type), intent(out) :: var

! Initialize the storage space and return

allocate(var%vars_2d(num_lons, num_lats, n2dflds), &
   var%vars_3d(num_lons, num_levs, num_lats, n3dflds), &
   var%tracers(num_lons, num_levs, num_lats, num_tracers))

end subroutine init_model_instance



  subroutine end_model_instance(var)
!=======================================================================
! subroutine end_model_instance(var)
!
! Ends an instance of a cam model state variable


type(model_type), intent(inout) :: var

deallocate(var%vars_2d, var%vars_3d, var%tracers)

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

integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_type

integer  :: indx, num_per_col, col_num, col_elem, lon_index, lat_index
real(r8) :: lon, lat, lev
integer  :: local_var_type, var_type_temp

! Easier to compute with a 0 to size - 1 index
indx = index_in - 1

! Compute number of items per column
num_per_col = num_levs * (n3dflds + pcnst + pnats) + n2dflds

! What column is this index in
col_num  = indx / num_per_col 
col_elem = indx - col_num * num_per_col

! What lon and lat index for this column
lon_index = col_num / num_lats
lat_index = col_num - lon_index * num_lats

! Get actual lon lat values from static_init arrays ???
lon = lons(lon_index + 1)
lat = lats(lat_index + 1)

! Now figure out which beast in column this is
! Surface pressure is the first element
lev = (col_elem - 1) / (n3dflds + pcnst + pnats) + 1
if(col_elem == 0) then
   local_var_type = TYPE_PS
   lev = -1
else
   var_type_temp = mod(col_elem - 1, n3dflds + pcnst + pnats)
   if(var_type_temp == 0) then
      local_var_type = TYPE_T
   else if(var_type_temp == 1) then
      local_var_type = TYPE_U
   else if(var_type_temp == 2) then
      local_var_type = TYPE_V
   else
      local_var_type = TYPE_Q
   endif
endif

!write(*, '(1x,3(f6.2,1x),i3)') lon, lat, lev, local_var_type

! Since CAM has pure pressure at the top, pure sigma at the
! bottom and hybrid in-between, we are just referring to
! the LEVEL index for the vertical. The coefficients to reconstruct
! pressure, height, etc. are available in the netCDF files.

location = set_location(lon, lat, lev, 1)  ! 1 == level (indexical)

! If the type is wanted, return it
if(present(var_type)) var_type = local_var_type


end subroutine get_state_meta_data



!function model_interpolate(x, location, type)
!=======================================================================
!
!real :: model_interpolate
!real, intent(in) :: x(:)
!type(location_type), intent(in) :: location
!integer, intent(in) :: type
!
!integer :: lon_below, lon_above, lat_below, lat_above, i
!real :: bot_lon, top_lon, delta_lon, bot_lat, top_lat
!real :: lon_fract, lat_fract, val(2, 2), temp_lon, a(2)
!real :: lon, lat, level, lon_lat_lev(3)
!
!! First version only allows (level specifies vertical???)
!lon_lat_lev = get_location(location)
!lon = lon_lat_lev(1); lat = lon_lat_lev(2); level = lon_lat_lev(3)
!
!! Get lon and lat grid specs, num_lons, num_lats are globally defined for cam
!   bot_lon = lons(1)
!   top_lon = lons(num_lons)
!   delta_lon = lons(2) - lons(1)
!   bot_lat = lats(1)
!   top_lat = lats(num_lats)
!
!! Compute bracketing lon indices
!if(lon >= bot_lon .and. lon <= top_lon) then
!   lon_below = int((lon - bot_lon) / delta_lon) + 1
!   lon_above = lon_below + 1
!   lon_fract = (lon - ((lon_below - 1) * delta_lon + bot_lon)) / delta_lon
!else
!! At wraparound point
!   lon_below = num_lons
!   lon_above = 1
!   if(lon < bot_lon) then
!      temp_lon = lon + 360.0
!   else
!      temp_lon = lon
!   endif
!   lon_fract = (temp_lon - top_lon) / delta_lon
!endif
!
!
!! Next, compute neighboring lat rows
!! NEED TO BE VERY CAREFUL ABOUT POLES; WHAT'S BEING DONE MAY BE WRONG
!! Inefficient search used for latitudes in Gaussian grid. Might want to speed up.
!if(lat >= bot_lat .and. lat <= top_lat) then
!
!   do i = 2, num_lats
!      if(lat <= lats(i)) then
!         lat_above = i
!         lat_below = i - 1
!         lat_fract = (lat - lats(lat_below)) / (lats(lat_above) - lats(lat_below))
!         goto 20
!      end if 
!   end do
!
!else if(lat <= bot_lat) then
!! South of bottom lat NEED TO DO BETTER: NOT REALLY REGULAR
!   lat_below = 1
!   lat_above = 2
!   lat_fract = 0.0
!else
!! North of top lat NEED TO DO BETTER: NOT REALLY REGULAR
!   lat_below = num_lats - 1
!   lat_above = num_lats
!   lat_fract = 1.0
!endif
!
!! Level is obvious for now
!
!! Now, need to find the values for the four corners
!20 continue
!val(1, 1) =  get_val(x, lon_below, lat_below, int(level), type)
!val(1, 2) =  get_val(x, lon_below, lat_above, int(level), type)
!val(2, 1) =  get_val(x, lon_above, lat_below, int(level), type)
!val(2, 2) =  get_val(x, lon_above, lat_above, int(level), type)
!
!! Do the weighted average for interpolation
!!write(*, *) 'fracts ', lon_fract, lat_fract
!do i = 1, 2
!   a(i) = lon_fract * val(2, i) + (1.0 - lon_fract) * val(1, i)
!end do
!
!model_interpolate = lat_fract * a(2) + (1.0 - lat_fract) * a(1)
!
!end function model_interpolate



  subroutine model_interpolate(x, location, type, interp_val, istatus, rstatus)
!=======================================================================
! subroutine model_interpolate(x, location, type)
!

real(r8) :: interp_val
real(r8), intent(in) :: x(:)
type(location_type), intent(in) :: location
integer, intent(in) :: type
integer, optional, intent(out) :: istatus
real(r8), optional, intent(out) :: rstatus

integer :: lon_below, lon_above, lat_below, lat_above, i
real(r8) :: bot_lon, top_lon, delta_lon, bot_lat, top_lat
real(r8) :: lon_fract, lat_fract, val(2, 2), temp_lon, a(2)
real(r8) :: lon, lat, level, lon_lat_lev(3), pressure

! Would it be better to pass state as prog_var_type (model state type) to here?
! As opposed to the stripped state vector. YES. This would give time interp.
! capability here; do we really want this or should it be pushed up?

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
   if(present(rstatus)) then
      call get_val(val(1, 1), x, lon_below, lat_below, nint(level), type, istatus, rstatus)
      if (istatus /= 1) call get_val(val(1, 2), x, lon_below, lat_above, nint(level), type, &
          istatus, rstatus)
      if (istatus /= 1) call get_val(val(2, 1), x, lon_above, lat_below, nint(level), type, &
          istatus, rstatus)
      if (istatus /= 1) call get_val(val(2, 2), x, lon_above, lat_above, nint(level), type, &
          istatus, rstatus)
   else if(present(istatus)) then
      call get_val(val(1, 1), x, lon_below, lat_below, nint(level), type, istatus)
      if (istatus /= 1) call get_val(val(1, 2), x, lon_below, lat_above, nint(level), type, istatus)
      if (istatus /= 1) call get_val(val(2, 1), x, lon_above, lat_below, nint(level), type, istatus)
      if (istatus /= 1) call get_val(val(2, 2), x, lon_above, lat_above, nint(level), type, istatus)
   else
      call get_val(val(1, 1), x, lon_below, lat_below, nint(level), type)
      call get_val(val(1, 2), x, lon_below, lat_above, nint(level), type)
      call get_val(val(2, 1), x, lon_above, lat_below, nint(level), type)
      call get_val(val(2, 2), x, lon_above, lat_above, nint(level), type)
   endif

else
! Case of pressure specified in vertical
   if(present(rstatus)) then
      call get_val_pressure(val(1, 1), x, lon_below, lat_below, pressure, type, istatus, rstatus)
      if (istatus /= 1) call get_val_pressure(val(1, 2), x, lon_below, lat_above, pressure, type, &
          istatus, rstatus)
      if (istatus /= 1) call get_val_pressure(val(2, 1), x, lon_above, lat_below, pressure, type, &
          istatus, rstatus)
      if (istatus /= 1) call get_val_pressure(val(2, 2), x, lon_above, lat_above, pressure, type, &
          istatus, rstatus)
   else if(present(istatus)) then
      call get_val_pressure(val(1, 1), x, lon_below, lat_below, pressure, type, istatus)
      if (istatus /= 1) call get_val_pressure(val(1, 2), x, lon_below, lat_above, pressure, type, &
          istatus)
      if (istatus /= 1) call get_val_pressure(val(2, 1), x, lon_above, lat_below, pressure, type, &
          istatus)
      if (istatus /= 1) call get_val_pressure(val(2, 2), x, lon_above, lat_above, pressure, type, &
          istatus)
   else
      call get_val_pressure(val(1, 1), x, lon_below, lat_below, pressure, type)
      call get_val_pressure(val(1, 2), x, lon_below, lat_above, pressure, type)
      call get_val_pressure(val(2, 1), x, lon_above, lat_below, pressure, type)
      call get_val_pressure(val(2, 2), x, lon_above, lat_above, pressure, type)
   endif
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
if( (present(istatus) .and. istatus /= 1) .or. .not.present(istatus) ) then
   do i = 1, 2
      a(i) = lon_fract * val(2, i) + (1.0 - lon_fract) * val(1, i)
   end do
   interp_val = lat_fract * a(2) + (1.0 - lat_fract) * a(1)
else
   interp_val = 0.
endif

end subroutine model_interpolate



  subroutine get_val_pressure(val, x, lon_index, lat_index, pressure, type, istatus, rstatus)
!=======================================================================
! function get_val_pressure(x, lon_index, lat_index, pressure, type)
!
! Gets the vertically interpolated value on pressure for variable type
! at lon_index, lat_index horizontal grid point
!
! This version excludes observations below lowest level pressure and above
! highest level pressure.


real(r8) :: val
real(r8), intent(in) :: x(:), pressure
integer, intent(in) :: lon_index, lat_index, type 
integer, optional, intent(out) :: istatus
real(r8), optional, intent(out) :: rstatus

real(r8) :: ps(1), pfull(1, num_levs), fraction
type(location_type) :: ps_location
integer :: top_lev, bot_lev, i, k
real(r8) :: bot_val, top_val, ps_lon

! Need to get the surface pressure at this point. Easy for A-grid.
if (present(rstatus)) then
   !get_val does not return a value of istatus or rstatus yet.  Should it? then test here.
   call get_val(ps(1), x, lon_index, lat_index, -1, 3, istatus, rstatus)
else if (present(istatus)) then
   call get_val(ps(1), x, lon_index, lat_index, -1, 3, istatus)
else
   call get_val(ps(1), x, lon_index, lat_index, -1, 3)
endif

! Next, get the values on the levels for this ps
call plevs_cam (1, 1, ps, pfull)

! Interpolate in vertical to get two bounding levels

if(pressure < pfull(1, 1) ) then
   if(present(istatus)) then
      istatus = 1
      if (present(rstatus)) rstatus = pfull(1, 1)
   endif
   val = 0.
else if(pressure > pfull(1, num_levs)) then
   if(present(istatus)) then
      istatus = 1
      if (present(rstatus)) rstatus = pfull(1, num_levs)
   endif
   val = 0.
else
   if(pressure < 20000.) then
      if(present(istatus)) then
         istatus = 2
         if (present(rstatus)) rstatus = pressure
      endif
   else
      if(present(istatus)) then
         istatus = 0
         if (present(rstatus)) rstatus = 0.
      endif
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

21 if (present(rstatus)) then
      !get_val does not return a value of istatus or rstatus yet.  Should it? then test here.
      call get_val(bot_val, x, lon_index, lat_index, bot_lev, type, istatus, rstatus)
      call get_val(top_val, x, lon_index, lat_index, top_lev, type, istatus, rstatus)
   else if (present(istatus)) then
      call get_val(bot_val, x, lon_index, lat_index, bot_lev, type, istatus)
      call get_val(top_val, x, lon_index, lat_index, top_lev, type, istatus)
   else
      call get_val(bot_val, x, lon_index, lat_index, bot_lev, type)
      call get_val(top_val, x, lon_index, lat_index, top_lev, type)
   endif 
   val = (1.0 - fraction) * bot_val + fraction * top_val
end if

end subroutine get_val_pressure



! function get_val_pressure(x, lon_index, lat_index, pressure, type, pflag)
!=======================================================================
! function get_val_pressure(x, lon_index, lat_index, pressure, type, pflag)
!
!! This version does extrapolations in log pressure below surface and above
!!  top model layer.
!
!
!real :: get_val_pressure
!real, intent(in) :: x(:), pressure
!integer, intent(in) :: lon_index, lat_index, type 
!integer, intent(inout) :: pflag
!
!real :: ps(1), pfull(1, num_levs), fraction
!type(location_type) :: ps_location
!integer :: top_lev, bot_lev, i, k
!real :: bot_val, top_val, ps_lon
!
!! Gets the vertically interpolated value on pressure for variable type
!! at lon_index, lat_index horizontal grid point
!
!! Need to get the surface pressure at this point. Easy for A-grid.
!ps = get_val(x, lon_index, lat_index, -1, 3)
!
! Next, get the values on the levels for this ps
!!! Kevin, you need to insert the appropriate call for you routine here
!! ps and pfull are currently 2 and 3D arrays to support Bgrid interface
!! I assume that you might want to change to a scalar and a 1D array?
!! call compute_pres_full(Dynam%Vgrid, ps, pfull)
!call plevs_cam (1, 1, ps, pfull)
!
!! Interpolate in vertical to get two bounding levels
!! What to do about pressures above top??? Just use the top for now.
!! Could extrapolate, but that would be very tricky. Might need to
!! reject somehow.

!!! kdr 10/3/03
!! We'll extrapolate for now, but may need to reject some data instead.
!
!if(pressure > ps(1)*1.05 .and. pflag == 0) then
!   write(53,'(/A,1p3E14.5)')  &
!        'WARNING; obs press > CAM surf press (&pfull)' &
!        ,pressure,ps(1),pfull(1,num_levs)
!   pflag = pflag + 1
!   write(53,'(A,1p3E14.5)')  &
!!!        'WARNING; obs press > CAM surf press (&pfull)' &
!!        ,pressure,ps(1),pfull(1,num_levs)
!endif
!if(pressure < pfull(1,1)*0.9 .and. pflag == 0) then
!   write(53,'(/A,1p3E14.5)')  &
!        'WARNING; obs press < .9xCAM top layer pressure ' &
!        ,pressure,pfull(1,1)
!   pflag = pflag + 1
!endif
!
!!if(pressure < pfull(1, 1) .and. pfull(1,top_lev) .ne. 0.) then
!! Extrapolate linearly in log(pressure)
!   top_lev = 1
!   bot_lev = 2
!   bot_val = get_val(x, lon_index, lat_index, bot_lev, type)
!   top_val = get_val(x, lon_index, lat_index, top_lev, type)
!   fraction = (bot_val - top_val)/(alog(pfull(1,bot_lev)/pfull(1,top_lev)))
!   get_val_pressure = top_val + fraction * (alog(pressure/pfull(1,top_lev)))
!   if (pressure > ps(1)*1.05) then
!   write(53,'(A,1p4E12.5)') '     pressure, pfull(1:2) =         ' &
!        ,pressure,pfull(1,top_lev),pfull(1,bot_lev)
!   write(53, '(A,2I3,1p4E12.5)') 'bot_lev, top_lev, fraction', bot_lev, top_lev, fraction
!!   write(53, '(A,1p3E12.5/)') '     bot_val, top_val, get_val_press', bot_val, top_val, get_val_pressure
!   endif

!!else if(pressure > pfull(1, num_levs)) then
!! Extrapolate linearly in log(pressure)
!   bot_lev = num_levs 
!   top_lev = bot_lev - 1
!   if(pfull(1,top_lev)*pfull(1,bot_lev) .ne. 0.  ) then
!!      bot_val = get_val(x, lon_index, lat_index, bot_lev, type)
!      top_val = get_val(x, lon_index, lat_index, top_lev, type)
!      fraction = (bot_val - top_val)/(alog(pfull(1,bot_lev)/pfull(1,top_lev)))
!      get_val_pressure = bot_val + fraction * (alog(pressure/pfull(1,bot_lev)))
!   if (pressure > ps(1)*1.05) then
!      write(53,'(A,1p4E12.5)') '     pressure, pfull(bot:top)  =         ' &
!           ,pressure, pfull(1,bot_lev),pfull(1,top_lev)
!!   write(53, '(A,2I3,1p4E12.5)') 'bot_lev, top_lev, fraction', bot_lev, top_lev, fraction
!      write(53, '(A,1p3E12.5)') '     get_val_press, bot_val, top_val ' &
!                                      ,get_val_pressure,bot_val,top_val
!!   endif
!    else
!      write(53, '(A,1p4E12.5)') 'pfull(1,top_lev)*pfull(1,bot_lev) = 0.' &
!           ,pfull(1,top_lev),pfull(1,bot_lev)
!    endif
!
!else
!! Search down through pressures
!   do i = 2, num_levs 
!      if(pressure < pfull(1, i)) then
!!         top_lev = i -1
!         bot_lev = i
!         fraction = (pfull(1, i) - pressure) / &
!            (pfull(1, i) - pfull(1, i - 1))
!         goto 21
!      endif
!   end do
!21 bot_val = get_val(x, lon_index, lat_index, bot_lev, type)
!   top_val = get_val(x, lon_index, lat_index, top_lev, type)
!   get_val_pressure = (1.0 - fraction) * bot_val + fraction * top_val
!end if
!
!!! Get the value at these two points
!!21 bot_val = get_val(x, lon_index, lat_index, bot_lev, type)
!!top_val = get_val(x, lon_index, lat_index, top_lev, type)
!!write(53, *) 'bot_lev, top_lev, fraction', bot_lev, top_lev, fraction
!!
!!get_val_pressure = (1.0 - fraction) * bot_val + fraction * top_val
!!write(53, *) 'bot_val, top_val, val', bot_val, top_val, get_val_pressure
!
!end function get_val_pressure




  subroutine get_val(val, x, lon_index, lat_index, level, type, istatus, rstatus)
!=======================================================================
! function get_val(x, lon_index, lat_index, level, type)
!

real(r8) :: val
real(r8), intent(in) :: x(:)
integer, intent(in) :: lon_index, lat_index, level, type
integer, optional, intent(out) :: istatus
real(r8), optional, intent(out) :: rstatus

integer :: per_col, indx

! Compute size of grid storage in a column; includes tracers
! Single 2D state vector is pressure
per_col = 1 + num_levs * n3tflds

! Find the starting index for this column
indx = per_col * (lat_index - 1 + (lon_index - 1) * num_lats)

! Pressure is first 
if(type == 3) then
   indx = indx + 1
else
! For interior fields compute the base for their level and add offset
   indx = indx + 1 + (level - 1) * n3tflds
! Temperature
   if(type == 4) then
      indx = indx + 1
! U wind component
   else if(type == 1) then
      indx = indx + 2
! V wind component
   else if(type == 2) then
      indx = indx + 3
! Tracers
   else if(type > 4) then
      indx = indx + type - 1
   end if
endif
   
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



  subroutine init_time(i_time)
!=======================================================================
! subroutine init_time(i_time)
!
! For now returns value of Time_init which is set in initialization routines.

type(time_type), intent(out) :: i_time

! Where should initial time come from here?
! WARNING: CURRENTLY SET TO 0
i_time = set_time(0, 0)

end subroutine init_time



  subroutine init_conditions(x)
!=======================================================================
! subroutine init_conditions(x)
!
! Reads in restart initial conditions  -- noop for CAM

real(r8), intent(inout) :: x(:)

end subroutine init_conditions



  subroutine model_get_close_states(o_loc, radius, nfound, indices, dist)
!=======================================================================
! subroutine model_get_close_states(o_loc, radius, nfound, indices, dist)
!

type(location_type), intent(in)  :: o_loc
real(r8),            intent(in)  :: radius
integer,             intent(out) :: nfound, indices(:)
real(r8),            intent(out) :: dist(:)

real(r8) :: loc_array(3), o_lon, o_lat
integer  :: num, max_size, i, j, num1
integer  :: hsize, num_per_col, col_base_index
integer,  allocatable :: lon_ind(:), lat_ind(:)
real(r8), allocatable :: close_dist(:)

write(*, *) 'in model_get_close_states', radius
loc_array = get_location(o_loc)
write(*, *) 'oloc is ', loc_array(:)

! Number found starts at 0
nfound = 0

! Assume that grid size is known from static initialized storage

! Num of close horizontal grid points starts at 0, too
num = 0
! For now, just allocate enough space for all grid points, may want
! to make this smaller at some point for big models.
max_size = num_lons * num_lats
allocate(lon_ind(max_size), lat_ind(max_size), close_dist(max_size))

! Look for close grid points on the 
call grid_close_states2(o_loc, lons, lats, num_lons, num_lats, radius, &
   num, lon_ind, lat_ind, close_dist)
! kdr write(*, *) 'back from grid_close_states num = ', num

! Compute size of grid storage for full levels
hsize = num_lons * num_lats
num_per_col = num_levs * (n3dflds + pcnst + pnats) + n2dflds

! Add all variables in this column to the close list with this distance
! kdr write(*, *) 'available space is size(indices) ', size(indices)
do i = 1, num
   col_base_index = ((lon_ind(i) - 1) * num_lats + lat_ind(i) - 1) * num_per_col
   do j = 1, num_per_col
      nfound = nfound + 1
      if(nfound <= size(indices)) indices(nfound) = col_base_index + j
      if(nfound <= size(dist)) dist(nfound) = close_dist(i)
   end do
end do
! kdr write(*, *) 'number at end is ', nfound

deallocate(lon_ind, lat_ind, close_dist)

end subroutine model_get_close_states




  subroutine grid_close_states2(o_loc, lons, lats, nlon, nlat, radius, &
                   num, close_lon_ind, close_lat_ind, close_dist)
!=======================================================================
! subroutine grid_close_states2(o_loc, lons, lats, nlon, nlat, radius, &
!                  num, close_lon_ind, close_lat_ind, close_dist)
!
!
! Finds close state points from a particular grid;


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
real(r8) :: glon, gdist

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
   loc        = set_location(glon, glat, glev, which_vert)
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
   loc        = set_location(glon, glat, glev, which_vert)
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
integer :: lonDimID, latDimID, ilevDimID, ScalarDimID, TracerDimID
integer :: MemberDimID, StateVarDimID, TimeDimID
integer :: lonVarID, latVarID, ilevVarID, hyaiVarID, hybiVarID, P0VarID, gwVarID
integer :: psVarID, TVarID, UVarID, VVarID, QVarID, ifld
integer :: TracerVarID, StateVarID, StateVarVarID
integer :: i
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
if ( num_tracers > 0 ) then
call check(nf90_def_dim(ncid=ncFileID, name="tracers",len = num_tracers, dimid = TracerDimID))
endif

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

! Number of Tracers
if ( num_tracers > 0 ) then
   call check(nf90_def_var(ncFileID, name="tracer", xtype=nf90_int, &
                                                  dimids=TracerDimID, varid=tracerVarID) )
   call check(nf90_put_att(ncFileID, tracerVarID, "long_name", "tracer identifier"))
endif

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
   ! We need to process the prognostic variables.
   ! I like the CAM philosophy of using nflds to declare the number of prognostic
   ! variables and an array of characters to specify them. This is clearly
   ! the way it must be for models with "lots"/variable numbers of params. 
   ! 
   ! I'd like to see the metadata handled the same way.
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

   ! PS ... ifld == 1 of nflds
   ifld = 1
   call check(nf90_def_var(ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real, &
              dimids = (/ lonDimID, latDimID, MemberDimID, unlimitedDimID /), &
              varid  = psVarID))
   call check(nf90_put_att(ncFileID, psVarID, "long_name", "surface pressure"))
   call check(nf90_put_att(ncFileID, psVarID, "units", "Pa"))
   call check(nf90_put_att(ncFileID, psVarID, "units_long_name", "pascals"))

   ! T ... ifld == 2 of nflds
   ifld = 2
   call check(nf90_def_var(ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real, &
         dimids = (/ lonDimID, ilevDimID, latDimID, MemberDimID, unlimitedDimID /), &
         varid  = TVarID))
   call check(nf90_put_att(ncFileID, TVarID, "long_name", "Temperature"))
   call check(nf90_put_att(ncFileID, TVarID, "units", "K"))

! kdr these ilevDimIDs were levDimIDs in TJH code
   ! U ... ifld == 3 of nflds
   ifld = 3
   call check(nf90_def_var(ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real, &
         dimids = (/ lonDimID, ilevDimID, latDimID, MemberDimID, unlimitedDimID /), &
         varid  = UVarID))
   call check(nf90_put_att(ncFileID, UVarID, "long_name", "Zonal Wind"))
   call check(nf90_put_att(ncFileID, UVarID, "units", "m/s"))

   ! V ... ifld == 4 of nflds
   ifld = 4
   call check(nf90_def_var(ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real, &
         dimids = (/ lonDimID, ilevDimID, latDimID, MemberDimID, unlimitedDimID /), &
         varid  = VVarID))
   call check(nf90_put_att(ncFileID, VVarID, "long_name", "Meridional Wind"))
   call check(nf90_put_att(ncFileID, VVarID, "units", "m/s"))

   ! Q ... ifld == 5 of nflds
   ifld = 5
   call check(nf90_def_var(ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real, &
         dimids = (/ lonDimID, ilevDimID, latDimID, MemberDimID, unlimitedDimID /), &
         varid  = QVarID))
   call check(nf90_put_att(ncFileID, QVarID, "long_name", "Specific Humidity"))
   call check(nf90_put_att(ncFileID, QVarID, "units", "kg/kg"))

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
if ( num_tracers > 0 ) then
   call check(nf90_put_var(ncFileID, tracerVarID, (/ (i,i=1,num_tracers) /) ))
endif

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
integer :: ifld

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
   
   TwoDVars : do ifld = 1, 1    ! PS always first? of one?
   
!   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   call check(NF90_inq_varid(ncFileID, trim(cflds(ifld)), ncfldid))
   call check(nf90_put_var( ncFileID, ncfldid, Var%vars_2d(:,:,1), &
             start=(/ 1, 1, copyindex, timeindex /), count=(/num_lons, num_lats, 1, 1/) ))
   enddo TwoDVars

   ThreeDVars : do ifld = 2,5
   call check(NF90_inq_varid(ncFileID, trim(cflds(ifld)), ncfldid))
   call check(nf90_put_var( ncFileID, ncfldid, Var%vars_3d(:,:,:,ifld-1), &
        start=(/ 1,1,1,copyindex,timeindex /), count=(/num_lons,num_levs,num_lats,1,1/) ))
   enddo ThreeDVars

   if ( num_tracers > 0 ) then
      call check(NF90_inq_varid(ncFileID,  "tracers",  ncfldid))
      call check(nf90_put_var( ncFileID,  ncfldid, Var%tracers(:,:,:,:), & 
           start=(/ 1,1,1,1,copyindex,timeindex /), &
           count=(/ num_lons, num_levs, num_lats, num_tracers,1,1/) ))
   endif
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
!=======================================================================
! subroutine pert_model_state(state, pert_state, interf_provided)
!
! Perturbs a model state for generating initial ensembles
! Returning interf_provided means go ahead and do this with uniform
! small independent perturbations.

real(r8),  intent(in) :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

interf_provided = .false.

end subroutine pert_model_state


!#######################################################################
! end of cam model_mod
!#######################################################################

end module model_mod
