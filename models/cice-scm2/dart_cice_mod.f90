! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module dart_cice_mod

use        types_mod, only : r8, rad2deg, PI, SECPERDAY, digits12
use time_manager_mod, only : time_type, get_date, set_date, get_time, set_time, &
                             set_calendar_type, get_calendar_string, &
                             print_date, print_time, operator(==), operator(-)
use    utilities_mod, only : get_unit, open_file, close_file, file_exist, &
                             register_module, error_handler, &
                             find_namelist_in_file, check_namelist_read, &
                             E_ERR, E_MSG, find_textfile_dims

use  netcdf_utilities_mod, only : nc_check


use typesizes
use netcdf

implicit none
private

public :: set_model_time_step,get_horiz_grid_dims, &
          get_ncat_dim, read_horiz_grid

character(len=*), parameter :: source   = "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

character(len=512) :: msgstring
logical, save :: module_initialized = .false.

character(len=256) :: ic_filename      = 'cice.r.nc'

contains

subroutine initialize_module

integer :: iunit, io

! Read calendar information
! In 'restart' mode, this is primarily the calendar type and 'stop'
! information. The time attributes of the restart file override
! the namelist time information.

! FIXME : Real observations are always GREGORIAN dates ...
! but stomping on that here gets in the way of running
! a perfect_model experiment for pre-1601 AD cases.
call set_calendar_type('gregorian')

! Make sure we have a cice restart file (for grid dims)
if ( .not. file_exist(ic_filename) ) then
   msgstring = 'dart_cice_mod: '//trim(ic_filename)//' not found'
   call error_handler(E_ERR,'initialize_module', &
          msgstring, source, revision, revdate)
endif

module_initialized = .true.

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

end subroutine initialize_module
!!!!!!!!!!!!!!!!
function set_model_time_step()

! the initialize_module ensures that the cice namelists are read.
! The restart times in the cice_in&restart_nml are used to define
! appropriate assimilation timesteps.
!
type(time_type) :: set_model_time_step

if ( .not. module_initialized ) call initialize_module

! Check the 'restart_option' and 'restart_n' to determine
! when we can stop the model
! CMB not sure if nday is actually different than ndays, no matter here though
!if ( (trim(restart_option) == 'ndays') .or. (trim(restart_option) == 'nday' ) ) then
!   set_model_time_step = set_time(0, restart_n) ! (seconds, days)
!else if ( trim(restart_option) == 'nyears' ) then
   ! FIXME ... CMB I guess we ignore it and make the freq 1 day anyway?
   set_model_time_step = set_time(0, 1) ! (seconds, days)
!else
!   call error_handler(E_ERR,'set_model_time_step', &
!              'restart_option must be ndays or nday', source, revision, revdate)
!endif

end function set_model_time_step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_horiz_grid_dims(Nx)

!
! Read the lon, lat grid size from the restart netcdf file.
! The actual grid file is a binary file with no header information.
!
! The file name comes from module storage ... namelist.

integer, intent(out) :: Nx   ! Number of Longitudes

integer :: grid_id, dimid, nc_rc

if ( .not. module_initialized ) call initialize_module

call nc_check(nf90_open(trim(ic_filename), nf90_nowrite, grid_id), &
         'get_horiz_grid_dims','open '//trim(ic_filename))

! Longitudes : get dimid for 'ni' or 'nlon', and then get value
nc_rc = nf90_inq_dimid(grid_id, 'ni', dimid)
if (nc_rc /= nf90_noerr) then
  msgstring = "unable to find either 'ni' or 'nlon' in file "//trim(ic_filename)
  call error_handler(E_ERR, 'get_horiz_grid_dims', msgstring, &
                         source,revision,revdate)
endif

call nc_check(nf90_inquire_dimension(grid_id, dimid, len=Nx), &
         'get_horiz_grid_dims','inquire_dimension ni '//trim(ic_filename))

call nc_check(nf90_close(grid_id), &
         'get_horiz_grid_dims','close '//trim(ic_filename) )

end subroutine get_horiz_grid_dims
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_ncat_dim(Ncat)

!
! Read the ncat size from the restart netcdf file.

integer, intent(out) :: Ncat   ! Number of categories in ice-thick dist

integer :: grid_id, dimid, nc_rc

if ( .not. module_initialized ) call initialize_module

call nc_check(nf90_open(trim(ic_filename), nf90_nowrite, grid_id), &
         'get_ncat_dim','open '//trim(ic_filename))

! ncat : get dimid for 'ncat' and then get value
nc_rc = nf90_inq_dimid(grid_id, 'ncat', dimid)
if (nc_rc /= nf90_noerr) then
   nc_rc = nf90_inq_dimid(grid_id, 'Ncat', dimid)
   if (nc_rc /= nf90_noerr) then
      msgstring = "unable to find either 'ncat' or 'Ncat' in file "//trim(ic_filename)
      call error_handler(E_ERR, 'get_horiz_grid_dims', msgstring, &
                         source,revision,revdate)
   endif
endif

call nc_check(nf90_inquire_dimension(grid_id, dimid, len=Ncat), &
         'get_ncat_dim','inquire_dimension ni '//trim(ic_filename))

! tidy up

call nc_check(nf90_close(grid_id), &
         'get_ncat_dim','close '//trim(ic_filename) )

end subroutine get_ncat_dim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_horiz_grid(nx, TLAT, TLON)

integer,                    intent(in)  :: nx
real(r8), dimension(nx), intent(out) :: TLAT, TLON

integer :: grid_id, reclength,VarId,status

if ( .not. module_initialized ) call initialize_module

! Check to see that the file exists.

if ( .not. file_exist(ic_filename) ) then
   msgstring = 'cice grid '//trim(ic_filename)//' not found'
   call error_handler(E_ERR,'read_horiz_grid', &
          msgstring, source, revision, revdate)
endif

! Open it and read them in the EXPECTED order.
! Actually, we only need the first two, so I'm skipping the rest.

call nc_check(nf90_open(trim(ic_filename), nf90_nowrite, grid_id), &
         'read_horiz_grid','open '//trim(ic_filename))
! Latitude
call nc_check(nf90_inq_varid(grid_id, 'tlat', VarId), &
         'read_horiz_grid','inquiring tlat from '//trim(ic_filename))
call nc_check(nf90_get_var(grid_id, VarId, TLAT, &
               start=(/1/), &
               count=(/nx/)), &
'read_horiz_grid','getting tlat from '//trim(ic_filename))
!Longitude
call nc_check(nf90_inq_varid(grid_id, 'tlon', VarId), &
'read_horiz_grid','inquiring tlon from '//trim(ic_filename))
call nc_check(nf90_get_var(grid_id, VarId, TLON, &
               start=(/1/), &
               count=(/nx/)), &
     'read_horiz_grid','getting tlon from '//trim(ic_filename))

call nc_check(nf90_close(grid_id), &
         'read_horiz_grid','close '//trim(ic_filename) )

TLAT = TLAT * rad2deg
TLON = TLON * rad2deg

! ensure [0,360) [-90,90]

where (TLON <   0.0_r8) TLON = TLON + 360.0_r8
where (TLON > 360.0_r8) TLON = TLON - 360.0_r8

where (TLAT < -90.0_r8) TLAT = -90.0_r8
where (TLAT >  90.0_r8) TLAT =  90.0_r8

end subroutine read_horiz_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module dart_cice_mod
