! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!>@todo FIXME we need to redefine GPSRO_REFRACTIVITY to be the non-local
!>operator and define a new obs type that is the local operator.  in that
!>case we don't need any additional metadata and can save time and space.
!> (or redefine them both, and make GPSRO_REFRACTIVITY an alias for the
!> non-local version)

program convert_cosmic_gps_cdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_cosmic_gps_cdf - program that reads a COSMIC GPS observation 
!                            profile and writes the data to a DART 
!                            observation sequence file. 
!
!     created June 2008 Ryan Torn, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use          types_mod, only : r8
use   time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, set_time,&
                               increment_time, get_time, set_date, operator(-),  &
                               print_date
use      utilities_mod, only : initialize_utilities, find_namelist_in_file,    &
                               check_namelist_read, nmlfileunit, do_nml_file,  &
                               get_next_filename, error_handler, E_ERR, E_MSG, &
                               find_textfile_dims, do_nml_term, finalize_utilities
use  netcdf_utilities_mod, only : nc_check
use       location_mod, only : VERTISHEIGHT, set_location
use   obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq,         &
                               static_init_obs_sequence, init_obs, destroy_obs,   &
                               write_obs_seq, init_obs_sequence,                  &
                               destroy_obs_sequence, set_obs_values, set_obs_def, &
                               set_copy_meta_data, set_qc, set_qc_meta_data
use   obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_type_of_obs, &
                               set_obs_def_error_variance, set_obs_def_location,        &
                               set_obs_def_key
use    obs_def_gps_mod, only : set_gpsro_ref
use       obs_kind_mod, only : GPSRO_REFRACTIVITY     ! also GPSRO_BENDING_ANGLE
use  obs_utilities_mod, only : add_obs_to_seq

use           netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"


integer, parameter ::   num_copies = 1,   &   ! number of copies in sequence
                        num_qc     = 1        ! number of QC entries

character (len=512) :: msgstring, msgstring2, next_infile
character (len=NF90_MAX_NAME)  :: name
character (len=6)   :: subset
integer :: ncid, varid, nlevels, k, nfiles, num_new_obs, oday, osec, &
           iyear, imonth, iday, ihour, imin, isec, glat, glon, zloc, obs_num, &
           io, iunit, nobs, filenum, dummy, numrejected
character (len=1) :: badqc
logical :: file_exist, first_obs, from_list = .false.
real(r8) :: hght_miss, refr_miss, azim_miss, benda_miss, oerr,   & 
            qc, lato, lono, hghto, refro, azimo, wght, nx, ny,   & 
            nz, rfict, obsval, phs, obs_val(1), qc_val(1)

real(r8), allocatable :: lat(:), lon(:), hght(:), refr(:), azim(:), & 
                         hghtp(:), refrp(:), benda(:)

type(obs_def_type)      :: obs_def
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time

!------------------------------------------------------------------------
!  Declare namelist parameters
!------------------------------------------------------------------------

integer, parameter :: NMAXLEVELS = 200   !  max number of observation levels

logical  :: local_operator         = .true.   ! see html file for more on non/local
logical  :: use_original_kuo_error = .false.  ! alternative is use lidia c's version
real(r8) :: obs_levels(NMAXLEVELS) = -1.0_r8
real(r8) :: ray_ds                 = 5000.0_r8    ! delta stepsize (m) along ray, nonlocal op
real(r8) :: ray_htop               = 15000.0_r8 ! max height (m) for nonlocal op
character(len=256) :: gpsro_netcdf_file     = 'cosmic_gps_input.nc'
character(len=256) :: gpsro_netcdf_filelist = ''
character(len=256) :: gpsro_out_file        = 'obs_seq.gpsro'

namelist /convert_cosmic_gps_nml/ obs_levels, local_operator, ray_ds,   &
                                  ray_htop, gpsro_netcdf_file,          &
                                  gpsro_netcdf_filelist, gpsro_out_file, &
                                  use_original_kuo_error

! initialize some values
obs_num = 1
qc = 0.0_r8
first_obs = .true.
call set_calendar_type(GREGORIAN)

!  read the necessary parameters from input.nml
call initialize_utilities()

call find_namelist_in_file("input.nml", "convert_cosmic_gps_nml", iunit)
read(iunit, nml = convert_cosmic_gps_nml, iostat = io)
call check_namelist_read(iunit, io, "convert_cosmic_gps_nml")

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=convert_cosmic_gps_nml)
if (do_nml_term()) write(     *     , nml=convert_cosmic_gps_nml)

! namelist checks for sanity
if (.not. use_original_kuo_error .and. .not. local_operator) then
  call error_handler(E_ERR, 'convert_cosmic_gps_cdf',                     &
                     'New error values only implemented for local operator', &
                     source, revision, revdate)
endif

!  count observation levels, make sure observation levels increase from 0
nlevels = 0
LEVELLOOP: do k = 1, NMAXLEVELS
  if ( obs_levels(k) == -1.0_r8 )  exit LEVELLOOP
  nlevels = k
end do LEVELLOOP
do k = 2, nlevels
  if ( obs_levels(k-1) >= obs_levels(k) ) then
    call error_handler(E_ERR, 'convert_cosmic_gps_cdf',       &
                       'Observation levels should increase',  &
                       source, revision, revdate)
  end if
end do

! cannot have both a single filename and a list; the namelist must
! shut one off.
if (gpsro_netcdf_file /= '' .and. gpsro_netcdf_filelist /= '') then
  call error_handler(E_ERR, 'convert_cosmic_gps_cdf',                     &
                     'One of gpsro_netcdf_file or filelist must be NULL', &
                     source, revision, revdate)
endif
if (gpsro_netcdf_filelist /= '') from_list = .true.

! need to know a reasonable max number of obs that could be added here.
if (from_list) then
   call find_textfile_dims(gpsro_netcdf_filelist, nfiles, dummy)
   num_new_obs = nlevels * nfiles
else
   num_new_obs = nlevels
endif

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
inquire(file=gpsro_out_file, exist=file_exist)
if ( file_exist ) then

   write(msgstring, *) "found existing obs_seq file, appending to ", trim(gpsro_out_file)
   write(msgstring2, *) "adding up to a maximum of ", num_new_obs, " new observations"
   call error_handler(E_MSG, 'convert_cosmic_gps_cdf', msgstring, &
                      source, revision, revdate, text2=msgstring2)
   call read_obs_seq(gpsro_out_file, 0, 0, num_new_obs, obs_seq)

else

   write(msgstring,  *) "no previously existing obs_seq file, creating ", trim(gpsro_out_file)
   write(msgstring2, *) "with up to a maximum of ", num_new_obs, " observations"
   call error_handler(E_MSG, 'convert_cosmic_gps_cdf', msgstring, &
                      source, revision, revdate, text2=msgstring2)

   call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)
   call set_copy_meta_data(obs_seq, 1, 'COSMIC GPS Observation')
   call set_qc_meta_data(obs_seq, 1, 'COSMIC QC')

end if


! main loop that does either a single file or a list of files.
! the data is currently distributed as a single profile per file.

filenum = 1
numrejected = 0

fileloop: do      ! until out of files

   ! get the single name, or the next name from a list
   if (from_list) then 
      next_infile = get_next_filename(gpsro_netcdf_filelist, filenum)
   else
      next_infile = gpsro_netcdf_file
      if (filenum > 1) next_infile = ''
   endif
   if (next_infile == '') exit fileloop
  
   !  open the next occultation profile
   call nc_check( nf90_open(next_infile, nf90_nowrite, ncid), 'file open', next_infile)

   ! this is new with the cosmic 2013 reprocessed obs.  they have included profiles
   ! which have failed their QC checking but added a global attribute 'bad' which
   ! is char "1" for bad profiles.   if there is no 'bad' attribute, or if there is one
   ! and the value is "0", then it is ok to proceed.  otherwise close the file and
   ! cycle the loop to pick up the next file.

   io = nf90_get_att(ncid,nf90_global,'bad', badqc)
   if ((io == nf90_noerr) .and. (badqc /= "0")) then
      ! profile with a bad QC
      numrejected = numrejected + 1
      call nc_check( nf90_close(ncid) , 'close file', next_infile)
      filenum = filenum + 1
      cycle fileloop 
   endif

   ! process the profile
   call nc_check( nf90_get_att(ncid,nf90_global,'year',  iyear) ,'get_att year',   next_infile)
   call nc_check( nf90_get_att(ncid,nf90_global,'month', imonth),'get_att month',  next_infile)
   call nc_check( nf90_get_att(ncid,nf90_global,'day',   iday)  ,'get_att day',    next_infile)
   call nc_check( nf90_get_att(ncid,nf90_global,'hour',  ihour) ,'get_att hour',   next_infile)
   call nc_check( nf90_get_att(ncid,nf90_global,'minute',imin)  ,'get_att minute', next_infile)
   call nc_check( nf90_get_att(ncid,nf90_global,'second',isec)  ,'get_att second', next_infile)
   
   time_obs = set_date(iyear, imonth, iday, ihour, imin, isec)
   call get_time(time_obs,  osec, oday)
   
   call nc_check( nf90_inq_dimid(ncid, "MSL_alt", varid),          'inq dimid MSL_alt', next_infile)
   call nc_check( nf90_inquire_dimension(ncid, varid, name, nobs), 'inq dim MSL_alt',   next_infile)
   
   call nc_check( nf90_get_att(ncid,nf90_global,'lon',  glon) ,'get_att lon',   next_infile)
   call nc_check( nf90_get_att(ncid,nf90_global,'lat',  glat) ,'get_att lat',   next_infile)
   call nc_check( nf90_get_att(ncid,nf90_global,'rfict',rfict),'get_att rfict', next_infile)
   rfict = rfict * 1000.0_r8
   
   allocate(hghtp(nobs)) ;  allocate(refrp(nobs))
   allocate( lat(nobs))  ;  allocate( lon(nobs))
   allocate(hght(nobs))  ;  allocate(refr(nobs))
   allocate(azim(nobs))  ;  allocate(benda(nobs))
   
   ! read the latitude array
   call nc_check( nf90_inq_varid(ncid, "Lat", varid) ,'inq varid Lat', next_infile)
   call nc_check( nf90_get_var(ncid, varid, lat)     ,'get var   Lat', next_infile)
   
   ! read the latitude array
   call nc_check( nf90_inq_varid(ncid, "Lon", varid) ,'inq varid Lon', next_infile)
   call nc_check( nf90_get_var(ncid, varid, lon)     ,'get var   Lon', next_infile)
   
   ! read the altitude array
   call nc_check( nf90_inq_varid(ncid, "MSL_alt", varid) ,'inq varid MSL_alt', next_infile)
   call nc_check( nf90_get_var(ncid, varid, hght)        ,'get_var   MSL_alt', next_infile)
   call nc_check( nf90_get_att(ncid, varid, '_FillValue', hght_miss) ,'get_att _FillValue MSL_alt', next_infile)
   
   ! read the refractivity
   call nc_check( nf90_inq_varid(ncid, "Ref", varid) ,'inq varid Ref', next_infile)
   call nc_check( nf90_get_var(ncid, varid, refr)    ,'get var   Ref', next_infile)
   call nc_check( nf90_get_att(ncid, varid, '_FillValue', refr_miss) ,'get_att _FillValue Ref', next_infile)
   
   ! note about bending angle:
   !  it is currently unused but people are interested in trying to assimilate this 
   !  instead of refractivity.  i don't think we have a forward operator for it yet.

   ! read the bending angle
   call nc_check( nf90_inq_varid(ncid, "Bend_ang", varid) ,'inq varid Bend_ang', next_infile)
   call nc_check( nf90_get_var(ncid, varid, benda)    ,'get var   Bend_ang', next_infile)
   call nc_check( nf90_get_att(ncid, varid, '_FillValue', benda_miss) ,'get_att _FillValue benda', next_infile)
   
   ! read the azimuth of occultation plane - only used for non-local operator
   call nc_check( nf90_inq_varid(ncid, "Azim", varid) ,'inq varid Azim', next_infile)
   call nc_check( nf90_get_var(ncid, varid, azim)     ,'get var   Azim', next_infile)
   call nc_check( nf90_get_att(ncid, varid, '_FillValue', azim_miss) ,'get_att _FillValue Azim', next_infile)
   
   call nc_check( nf90_close(ncid) , 'close file', next_infile)
   
   !where(hght(1: < -998)) then
   !   write(*,*) 'height array:'
   !   write(*,*) hght
   !endwhere
   !where (any(refr < -998)) then
   !  write(*,*) 'refractivity array:'
   !  write(*,*) refr
   !endwhere
   !where (any(benda < -998)) then
   !  write(*,*) 'bending angle array:'
   !  write(*,*) benda
   !endwhere

   ! convert units here.
   hghtp(:) = hght(:) * 1000.0_r8
   refrp(:) = refr(:) * 1.0e-6_r8
    
   obsloop2: do k = 1, nlevels
   
     call interp_height_wght(hght, obs_levels(k), nobs, zloc, wght)
     if ( zloc < 1 )  cycle obsloop2
   
     ! lon(zloc) and lon(zloc+1) range from -180 to +180
     ! call a subroutine to handle the wrap point.
     lono = compute_lon_wrap(lon(zloc), lon(zloc+1), wght)
     lato  = wght * lat(zloc)  + (1.0_r8 - wght) * lat(zloc+1)
     hghto = wght * hght(zloc) + (1.0_r8 - wght) * hght(zloc+1)
     hghto = hghto * 1000.0_r8
     refro = wght * refr(zloc) + (1.0_r8 - wght) * refr(zloc+1)
     azimo = wght * azim(zloc) + (1.0_r8 - wght) * azim(zloc+1)
   
     if ( local_operator ) then
   
        nx        = 0.0_r8
        ny        = 0.0_r8
        nz        = 0.0_r8
        ray_ds    = 0.0_r8
        ray_htop  = 0.0_r8
        rfict     = 0.0_r8
   
        obsval = refro
        if (use_original_kuo_error) then
           oerr   = 0.01_r8 * ref_obserr_kuo_percent(hghto * 0.001_r8) * obsval
        else
           oerr   = gsi_refractivity_error(hghto, lato, is_it_global=.true., factor=1.0_r8)
        endif
        subset = 'GPSREF'
   
     else
   
       !  compute tangent unit vector
       call tanvec01(lono, lato, azimo, nx, ny, nz)
   
       !  compute the excess phase
       call excess(refrp, hghtp, lono, lato, hghto, nx, & 
                   ny, nz, rfict, ray_ds, ray_htop, phs, nobs)
   
       !  if too high, phs will return as 0.  cycle loop here.
       if (phs <= 0) cycle obsloop2
   
       obsval = phs
       oerr   = 0.01_r8 * excess_obserr_percent(hghto * 0.001_r8) * obsval
       !DEBUG print *, 'hghto,obsval,perc,err=', hghto, obsval, &
       !DEBUG          excess_obserr_percent(hghto * 0.001_r8), oerr
       subset = 'GPSEXC'
   
     end if
   
     call set_gpsro_ref(obs_num, nx, ny, nz, rfict, ray_ds, ray_htop, subset)
     call set_obs_def_location(obs_def,set_location(lono,lato,hghto,VERTISHEIGHT))
     call set_obs_def_type_of_obs(obs_def, GPSRO_REFRACTIVITY)
     call set_obs_def_time(obs_def, set_time(osec, oday))
     call set_obs_def_error_variance(obs_def, oerr * oerr)
     call set_obs_def_key(obs_def, obs_num)
     call set_obs_def(obs, obs_def)
   
     obs_val(1) = obsval
     call set_obs_values(obs, obs_val)
     qc_val(1)  = qc
     call set_qc(obs, qc_val)

     call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

     obs_num = obs_num+1

   end do obsloop2

  ! clean up and loop if there is another input file
  deallocate( lat, lon, hght, refr, azim, benda, hghtp, refrp)

  filenum = filenum + 1

end do fileloop

! done with main loop.  if we added any new obs to the sequence, write it out.
if (obs_num > 1) call write_obs_seq(obs_seq, gpsro_out_file)

! cleanup memory
call destroy_obs_sequence(obs_seq)
call destroy_obs(obs)   ! do not destroy prev_obs, which is same as obs

write(msgstring, *) 'processed ', filenum-1, ' total profiles'
call error_handler(E_MSG, 'convert_cosmic_gps_obs', msgstring, source, revision, revdate)

if (numrejected > 0) then
   write(msgstring,  *) numrejected, ' profiles rejected for bad incoming quality control'
   call error_handler(E_MSG, 'convert_cosmic_gps_obs', msgstring, source, revision, revdate)
endif

call finalize_utilities()

! END OF MAIN ROUTINE

contains

! local subroutines/functions follow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   excess - subroutine that computes the excess phase based on the 
!            method of Sergey et al. (2004), eq. 1 and 2
!
!    refp   - refractivity profile
!    hghtp  - height profile of refractivity observations
!    lon    - longitude of GPS perigee point
!    lat    - latitude of GPS perigee point
!    height - height of GPS perigee point
!    nx     - x component of tangent point of ray
!    ny     - y component of tangent point of ray
!    nz     - z component of tangent point of ray
!    rfict  - local curvature radius
!    ds     - increment of excess computation - intent(in) now
!    htop   - heighest level of computation - intent(in) now
!    phs    - excess phase
!    nobs   - number of points in profile
!
!     created June 2004, Hui Liu NCAR/MMM
!     modified Ryan Torn, NCAR/MMM  
!     updated Nancy Collins, NCAR/IMAGe
!       quit before computing segment when one endpoint is above top, not after
!       pass delta step and htop in as input args instead of out.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine excess(refp, hghtp, lon, lat, height, nx, ny, nz, &
                                    rfict, ds, htop, phs, nobs)

use  types_mod, only : r8

implicit none

integer, intent(in)   :: nobs
real(r8), intent(in)  :: refp(nobs), hghtp(nobs), lon, lat, height, nx, & 
                         ny, nz, rfict, ds, htop
real(r8), intent(out) :: phs

integer :: iter
real(r8) :: xo, yo, zo, ref00, ref1, ref2, ss, xx, yy, zz, height1, &
            lat1, lon1, dphs1, dphs2

! make sure height is not already above requested top of integration
if( height >= htop ) then
   phs = 0.0_r8           ! excess phase
   return
endif

!  convert location of perigee from geodetic to Cartesian coordinate
call geo2car(height, lat, lon, xo, yo, zo, rfict)

!  get refractivity at perigee
call ref_int(height, refp, hghtp, ref00, nobs)

!  integrate refractivity along a straight line path in cartesian coordinate
!   (x-xo)/a = (y-yo)/b = (z-zo)/c,  (a,b,c) is the line direction
ref1 = ref00  ;  ref2 = ref00

!    initialization
phs = 0.0_r8           ! excess phase
ss  = 0.0_r8           ! distance to perigee from a ray point
iter = 0

do 

  iter = iter + 1
  ss   = ss + ds

  !  integrate to one direction of the line for one step
  xx = xo + ss * nx
  yy = yo + ss * ny
  zz = zo + ss * nz

  !  convert the location of the point to geodetic coordinates
  !   height(m), lat, lon(deg)
  call car2geo(xx, yy, zz, height1, lat1, lon1, rfict)
  if( height1 >= htop ) exit  ! break out of loop if above level

  call ref_int(height1, refp, hghtp, ref00, nobs)
  dphs1 = (ref1 + ref00) * ds / 2.0_r8

  ref1 = ref00

  !  integrate to the other direction of the line for one step
  xx = xo - ss * nx
  yy = yo - ss * ny
  zz = zo - ss * nz

  call car2geo(xx, yy, zz, height1, lat1, lon1, rfict)
  call ref_int(height1, refp, hghtp, ref00, nobs)
  dphs2 = (ref2 + ref00) * ds / 2.0_r8

  ref2 = ref00

  phs = phs + dphs1 + dphs2

end do

return
end subroutine excess

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   car2geo - subroutine that converts cartesian coordinates to
!             geodetical.
!
!    x1    - x coordinate
!    x2    - y coordinate
!    x3    - z coordinate
!    s1    - geodetical height
!    s2    - geodetical latitude
!    s3    - geodetical longitude
!    rfict - local curvature radius
!
!     created Hui Liu NCAR/MMM
!     Modified June 2008, Ryan Torn NCAR/MMM 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine car2geo (x1, x2, x3, s1, s2, s3, rfict)

use types_mod, only : r8, rad2deg, pi

implicit none

real(r8), intent(in)  :: x1, x2, x3, rfict
real(r8), intent(out) :: s1, s2, s3

real(r8) :: rho, azmth

rho   = sqrt (x1 * x1 + x2 * x2 + x3 * x3)
s1    = rho - rfict
s2    = asin(x3 / rho)
azmth = atan2(x2, x1)
s3    = mod((azmth + 4.0_r8 * pi), 2.0_r8 * pi)

s2 = s2 * rad2deg
s3 = s3 * rad2deg

return
end subroutine car2geo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   geo2car - subroutine that converts geodetical coordinates to 
!             cartesian.
!
!    s1    - geodetical height
!    s2    - geodetical latitude
!    s3    - geodetical longitude
!    x1    - x coordinate
!    x2    - y coordinate
!    x3    - z coordinate
!    rfict - local curvature radius
!
!     created Hui Liu NCAR/MMM
!     modified June 2008, Ryan Torn NCAR/MMM 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine geo2car(s1, s2, s3, x1, x2, x3, rfict)

use types_mod, only : r8, deg2rad

implicit none

real(r8), intent(in)  :: s1, s2, s3, rfict
real(r8), intent(out) :: x1, x2, x3

real(r8) :: g3, g4

g3  = s1 + rfict
g4  = g3 * cos(s2 * deg2rad)

x1  = g4 * cos(s3 * deg2rad)
x2  = g4 * sin(s3 * deg2rad)
x3  = g3 * sin(s2 * deg2rad)

return
end subroutine geo2car

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   interp_height_wght - subroutine that finds the vertical levels
!                        closest to the desired height level and the
!                        weights to give to these levels to perform
!                        vertical interpolation.
!
!    hght  - height levels in column
!    level - height level to interpolate to
!    zgrid - index of lowest level for interpolation
!    wght  - weight to give to the lower level in interpolation
!    iz    - number of vertical levels
!
!     created June 2008, Ryan Torn NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interp_height_wght(hght, level, iz, zgrid, wght)

use        types_mod, only : r8

implicit none

integer,  intent(in)  :: iz
real(r8), intent(in)  :: hght(iz), level
integer,  intent(out) :: zgrid
real(r8), intent(out) :: wght

integer :: k, klev, kbot, ktop, kinc, kleva

if ( hght(1) > hght(iz) ) then
  kbot = iz  ;   ktop = 1   ;  kinc = -1   ;  kleva = 0
else
  kbot = 1   ;   ktop = iz  ;  kinc = 1    ;  kleva = 1
end if

if ( (hght(kbot) <= level) .AND. (hght(ktop) >= level) ) then

  do k = kbot, ktop, kinc  !  search for appropriate level
    if ( hght(k) > level ) then
      klev = k - kleva
      exit
    endif
  enddo

  ! compute the weights
  zgrid = klev
  wght  = (level-hght(klev+1)) / (hght(klev) - hght(klev+1))

else

  zgrid = -1
  wght  = 0.0_r8

endif

return
end subroutine interp_height_wght

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   ref_int - subroutine that calculates the refractivity at any 
!             particular point.
!
!    height - GPS observation height (m)
!    refp   - refractivity profile
!    hghtp  - height profile
!    lref   - local refractivity value
!    nobs   - number of observations in the profile
!
!     created Hui Liu NCAR/MMM
!     updated June 2008, Ryan TORN NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ref_int(height, refp, hghtp, lref, nobs)

use   types_mod, only : r8

implicit none

integer, intent(in)   :: nobs
real(r8), intent(in)  :: height, refp(nobs), hghtp(nobs)
real(r8), intent(out) :: lref

integer  :: bot_lev, k
real(r8) :: fract

bot_lev = -1
fract = 0.0_r8

! make sure it's not higher than the highest available level.
if (height > hghtp(1)) then
   write(msgstring, *) 'requested height is ', height, '; highest available is ', hghtp(1)
   call error_handler(E_ERR, 'bad level, above highest in file', &
                      source, revision, revdate, text2=msgstring)
endif

!  Search down through height levels
heights: do k = 2, nobs
  if ( height >= hghtp(k) ) then
    bot_lev = k
    fract = (hghtp(k) - height) / (hghtp(k) - hghtp(k-1))
    exit heights
  endif
end do heights

! the hghtp() array is currently an interpolated list of levels
! and on at least 1 PGI compiler version computing the lowest value 
! rounds enough that a height exactly equal to the lowest level 
! compares as less than instead of equal.  so test and if it's very very 
! close to the lowest level then return it as equal; otherwise it's 
! an internally inconsistent input file.
if (bot_lev < 0) then
  if (abs(height - hghtp(nobs)) < 0.00001_r8) then
     bot_lev = nobs
     fract = 0.0_r8
  else
     write(msgstring, *) 'requested height is ', height, '; lowest available is ', hghtp(nobs)
     call error_handler(E_ERR, 'bad level, below lowest in file', &
                       source, revision, revdate, text2=msgstring)
  endif
endif

lref = (1.0_r8 - fract) * refp(bot_lev) + fract * refp(bot_lev-1)

return
end subroutine ref_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   tanvec01 - subroutine that computes the unit vector tangent of the 
!              ray at the perigee.
!
!    lon0  - longitude of the tangent point
!    lat0  - latitude of the tangent point
!    azim0 - angle between occultation plane from north
!    uz    - x component of tangent vector at tangent point of array
!    uy    - y component of tangent vector at tangent point of array
!    uz    - z component of tangent vector at tangent point of array
!
!     created Hui Liu, NCAR/MMM
!     modified June 2008, Ryan Torn NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tanvec01(lon0, lat0, azim0, ux, uy, uz)

use   types_mod, only : r8, deg2rad

implicit none

real(r8), intent(in)  :: lon0, lat0, azim0
real(r8), intent(out) :: ux, uy, uz

real(r8) :: zz0(3), lon, lat, azim, rtp(3), rnm(3), rno(3), uon(3), vlen0

zz0(1) = 0.0_r8  ;  zz0(2) = 0.0_r8  ;  zz0(3) = 1.0_r8
lon = lon0 * deg2rad  ;  lat = lat0 * deg2rad  ;  azim = azim0 * deg2rad

rtp(1) = cos(lat) * cos(lon)
rtp(2) = cos(lat) * sin(lon)
rtp(3) = sin(lat)

!  compute unit vector normal to merdion plane through tangent point
call vprod(rtp, zz0, rnm)
vlen0 = sqrt(rnm(1)*rnm(1) + rnm(2)*rnm(2) + rnm(3)*rnm(3))
rnm(:) = rnm(:) / vlen0

!  compute unit vector toward north from perigee point
call vprod(rnm, rtp, rno)
vlen0 = sqrt(rno(1)*rno(1) + rno(2)*rno(2) + rno(3)*rno(3))
rno(:) = rno(:) / vlen0

!  rotate the vector rno around rtp for a single azim to get tangent vector
call spin(rno, rtp, azim, uon)
vlen0 = sqrt(uon(1)*uon(1) + uon(2)*uon(2) + uon(3)*uon(3))
ux = uon(1) / vlen0  ;  uy = uon(2) / vlen0  ;  uz = uon(3) / vlen0

return
end subroutine tanvec01

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   spin - subroutine that rotates vector v1 around vs clockwise by 
!          by a specified angle.
!
!    v1 - vector to rotate
!    vs - vector to rotate about
!     a - angle to rotate v1 around
!    v2 - output vector after rotation
!
!     created Hui Liu NCAR/MMM
!     modified June 2008, Ryan Torn NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spin(v1, vs, a, v2)

use   types_mod, only : r8, deg2rad

implicit none

real(r8), intent(in)  :: v1(3), vs(3), a
real(r8), intent(out) :: v2(3)

real(r8) :: vsabs, vsn(3), a1, a2, a3, s(3,3) 

! Calculation of the unit vector for the rotation
vsabs  = sqrt(vs(1)*vs(1) + vs(2)*vs(2) + vs(3)*vs(3))
vsn(:) = vs(:) / vsabs

! Calculation the rotation matrix
a1 = cos(a)  ;  a2 = 1.0_r8 - a1  ;  a3 = sin(a)
s(1,1) = a2 * vsn(1) * vsn(1) + a1
s(1,2) = a2 * vsn(1) * vsn(2) - a3 * vsn(3)
s(1,3) = a2 * vsn(1) * vsn(3) + a3 * vsn(2)
s(2,1) = a2 * vsn(2) * vsn(1) + a3 * vsn(3)
s(2,2) = a2 * vsn(2) * vsn(2) + a1
s(2,3) = a2 * vsn(2) * vsn(3) - a3 * vsn(1)
s(3,1) = a2 * vsn(3) * vsn(1) - a3 * vsn(2)
s(3,2) = a2 * vsn(3) * vsn(2) + a3 * vsn(1)
s(3,3) = a2 * vsn(3) * vsn(3) + a1

!  Compute the rotated vector
v2(:) = s(:,1) * v1(1) + s(:,2) * v1(2) + s(:,3) * v1(3)

return
end subroutine spin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   vprod - subroutine that computes the vector product of two vectors
!
!    x - first vector to take product of
!    y - second vector to take product of
!    z - vector product
!
!     created Hui Liu NCAR/MMM
!     modified June 2008, Ryan Torn NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vprod(x, y, z)
 
use   types_mod, only : r8

implicit none

real(r8), intent(in)  :: x(3), y(3)
real(r8), intent(out) :: z(3)

z(1) = x(2)*y(3) - x(3)*y(2)
z(2) = x(3)*y(1) - x(1)*y(3)
z(3) = x(1)*y(2) - x(2)*y(1)

return
end subroutine vprod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   compute_lon_wrap - interpolate between 2 longitude values taking
!                      into account the wrap at -180 degrees
!
!    lon1, lon2 - longitude in degrees between -180 and +180
!    weight - interpolation weight between lon1 and lon2 (0 to 1)
! 
!    returns interpolated longitude between 0 and 360 degrees.
!
! if the longitudes are the same sign (both negative or both positive)
! then do the interpolation with the original values.  if the signs
! are different then we need to decide if they are crossing 0 (where we
! still use the original values) or if they are crossing the -180/180 line
! and we have to wrap the negative value.

! to decide between the 0 and 180 cases, take the positive value and subtract 
! the negative value (which adds it on) and see if the sum is > 180.  if not, 
! we're at the 0 crossing and we do nothing.  if yes, then we add 360 to the 
! negative value and interpolate between two positive values.   in either case
! once we have the result, if it's < 0 add 360 so the longitude returned is
! between 0 and 360 in longitude.
!
! this does not try to do anything special if the profile is tracking directly
! over one of the poles.  this is because at the exact poles all longitudes are
! identical, so being off in the longitude in any direction won't be a large
! difference in real distance.
!
!     created nancy collins NCAR/IMAGe
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function compute_lon_wrap(lon1, lon2, weight)

real(r8), intent(in) :: lon1, lon2, weight
real(r8) :: compute_lon_wrap

real(r8) :: lon1a, lon2a, lono

! r/w temporaries in case we have to change the value.  
lon1a = lon1
lon2a = lon2


! if different signs and crossing the -180/180 boundary, add 360
! to the negative value.
if (lon1 <= 0.0_r8 .and. lon2 >= 0.0_r8) then
   if (lon2 - lon1 > 180.0_r8) lon1a = lon1a + 360.0_r8
else if (lon1 >= 0.0_r8 .and. lon2 <= 0.0_r8) then
   if (lon1 - lon2 > 180.0_r8) lon2a = lon2a + 360.0_r8
endif

! linear interpolation, and make return value between 0 and 360.
lono  = weight * lon1a  + (1.0_r8 - weight) * lon2a
if (lono < 0.0_r8) lono = lono + 360.0_r8

compute_lon_wrap = lono

end function compute_lon_wrap

! three different error options:

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   excess_obserr_percent - function that computes the observation 
!                           error percentage for a given height.
!
!    hght - height of excess observation (km)
!
!     created Hui Liu NCAR/MMM
!     updated June 2008, Ryan Torn NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function excess_obserr_percent(hght)

use types_mod, only : r8

implicit none

integer, parameter :: nobs_level = 18

real(r8), intent(in) :: hght

integer  :: k, k0
real(r8) :: hght0, exc_err(nobs_level), obs_ht(nobs_level), excess_obserr_percent

!-------------------------------------------------------------------------------
! obs height in km
data obs_ht/17.0_r8, 16.0_r8, 15.0_r8, 14.0_r8, 13.0_r8, 12.0_r8, & 
            11.0_r8, 10.0_r8,  9.0_r8,  8.0_r8,  7.0_r8,  6.0_r8, &
             5.0_r8,  4.0_r8,  3.0_r8,  2.0_r8,  1.0_r8,  0.0_r8/

data exc_err/0.2_r8,  0.2_r8,  0.2_r8,  0.2_r8,  0.2_r8,  0.2_r8, & 
             0.2_r8,  0.3_r8,  0.4_r8,  0.5_r8,  0.6_r8,  0.7_r8, &
             0.8_r8,  0.9_r8,  1.0_r8,  1.1_r8,  1.2_r8,  1.3_r8/

!-------------------------------------------------------------------------------

hght0 = max(hght,0.0_r8)

do k = 1, nobs_level
  if ( hght >= obs_ht(k) ) then
    k0 = k
    exit
  end if
end do

excess_obserr_percent = exc_err(k0) + (exc_err(k0) - exc_err(k0-1)) / & 
                        (obs_ht(k0) - obs_ht(k0-1)) * (hght - obs_ht(k0))

return
end function excess_obserr_percent

! routines below here were lifted from soyoung ha's version of
! the ncep bufr format converter for gps obs, including lidia's error spec

! soyoung says the error spec from Bill Kuo is only 1 of 5 possible profiles
! and most appropriate for the CONUS domain.  if used for global obs, it
! should be updated to include the other functions.  nsc jan 2016

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   ref_obserr_kuo_percent - function that computes the observation 
!                            error for a refractivity observation.
!                            These numbers are taken from a Kuo paper.
!
!    hght - height of refractivity observation
!
!     created Hui Liu NCAR/MMM
!     modified June 2008, Ryan Torn NCAR/MMM
!     modified August 2015, Soyoung Ha NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function ref_obserr_kuo_percent(hght)

use   types_mod, only : r8

implicit none

integer, parameter :: nobs_level = 22  !  maximum number of obs levels

real(r8), intent(in)  :: hght

integer  :: k, k0
real(r8) :: hght0, ref_err(nobs_level), obs_ht(nobs_level), ref_obserr_kuo_percent

!   observation error heights (km) and errors
data obs_ht/17.98_r8, 16.39_r8, 14.97_r8, 13.65_r8, 12.39_r8, 11.15_r8, &
             9.95_r8,  8.82_r8,  7.82_r8,  6.94_r8,  6.15_r8,  5.45_r8, & 
             4.83_r8,  4.27_r8,  3.78_r8,  3.34_r8,  2.94_r8,  2.59_r8, & 
             2.28_r8,  1.99_r8,  1.00_r8,  0.00_r8/

data ref_err/0.48_r8,  0.56_r8,  0.36_r8,  0.28_r8,  0.33_r8,  0.41_r8, & 
             0.57_r8,  0.73_r8,  0.90_r8,  1.11_r8,  1.18_r8,  1.26_r8, & 
             1.53_r8,  1.85_r8,  1.81_r8,  2.08_r8,  2.34_r8,  2.43_r8, & 
             2.43_r8,  2.43_r8,  2.43_r8,  2.43_r8/

hght0 = max(hght, 0.0_r8)

do k = 1, nobs_level
  
  k0 = k
  if ( hght0 >= obs_ht(k) )  exit

end do 

if(k0.eq.1) then    ! HA
   ref_obserr_kuo_percent = ref_err(k0)
else
   ref_obserr_kuo_percent = ref_err(k0) + (ref_err(k0) - ref_err(k0-1)) / &
                           (obs_ht(k0)-obs_ht(k0-1)) * (hght0-obs_ht(k0))
endif

return
end function ref_obserr_kuo_percent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   gsi_refractivity_error - function that computes the observation 
!                            error as in GSI
!
!    Input:
!    H   -- input real value geometric height [m]
!    lat -- latitude in degree 
!    is_it_global -- is your forecast model regional or global? (T or F)
!    factor -- quick way to scale all the errors up or down for testing
!
!    Output:
!    gsi_refractivity_error -- output refractivity observation error [N]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function gsi_refractivity_error(H, lat, is_it_global, factor)
 real(r8), intent(in)  :: H
 real(r8), intent(in)  :: lat
 logical,  intent(in)  :: is_it_global
 real(r8), intent(in)  :: factor
 real(r8)              :: gsi_refractivity_error

 real(r8) :: zkm, rerr
 
 zkm = H * 0.001       ! height in km
 rerr = 1.0_r8

 if(is_it_global) then    ! for global

 if((lat >= 20.0_r8) .or. (lat <= -20.0_r8)) then
    rerr = -1.321_r8 + 0.341_r8 * zkm - 0.005_r8 * zkm ** 2
 else
    if(zkm > 10.0_r8) then
       rerr = 2.013_r8 - 0.060_r8 * zkm + 0.0045_r8 * zkm ** 2
    else
       rerr = -1.18_r8 + 0.058_r8 * zkm + 0.025_r8 * zkm ** 2
    endif
 endif

 else     ! for regional 

 if((lat >= 20.0_r8) .or. (lat <= -20.0_r8)) then
    if(zkm > 10.0_r8) then
       rerr = -1.321_r8 + 0.341_r8 * zkm - 0.005_r8 * zkm ** 2
    else
       rerr = -1.2_r8 + 0.065_r8 * zkm + 0.021_r8 * zkm ** 2
    endif
 endif

 endif    ! if(is_it_global) then

 if (factor /= 1.0_r8) then
    gsi_refractivity_error = 1./(abs(exp(rerr))*factor)
 else
    gsi_refractivity_error = 1./abs(exp(rerr))
 endif

 return
end function gsi_refractivity_error

! end error options

! next function is currently unused in this converter:

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    compute_geopotential_height 
!    subroutine converts geometric height to geopotential height
!
!    Input:
!    H   -- input real value geometric height [m]
!    lat -- latitude in degree 
!
!    Output:
!    Z -- output real value geopotential height [m]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function compute_geopotential_height(H, lat)
 real(r8), intent(in)  :: H 
 real(r8), intent(in)  :: lat
 real(r8)              :: compute_geopotential_height

! -----------------------------------------------------------------------*/
   real(r8) :: pi2, latr
   real(r8) :: semi_major_axis, semi_minor_axis, grav_polar, grav_equator
   real(r8) :: earth_omega, grav_constant, flattening, somigliana
   real(r8) :: grav_ratio, sin2, termg, termr, grav, eccentricity

!  Parameters below from WGS-84 model software inside GPS receivers.
   parameter(semi_major_axis = 6378.1370d3)    ! (m)
   parameter(semi_minor_axis = 6356.7523142d3) ! (m)
   parameter(grav_polar = 9.8321849378)        ! (m/s2)
   parameter(grav_equator = 9.7803253359)      ! (m/s2)
   parameter(earth_omega = 7.292115d-5)        ! (rad/s)
   parameter(grav_constant = 3.986004418d14)   ! (m3/s2)
   parameter(grav = 9.80665d0)                 ! (m/s2) WMO std g at 45 deg lat
   parameter(eccentricity = 0.081819d0)        ! unitless
   parameter(pi2 = 3.14159265358979d0/180.d0)

!  Derived geophysical constants
   parameter(flattening = (semi_major_axis-semi_minor_axis) / semi_major_axis)

   parameter(somigliana = (semi_minor_axis/semi_major_axis)*(grav_polar/grav_equator)-1.d0)

   parameter(grav_ratio = (earth_omega*earth_omega * &
                semi_major_axis*semi_major_axis * semi_minor_axis)/grav_constant)

!  Sanity Check
   if(lat.gt.90 .or. lat.lt.-90) then
      print*, 'compute_geopotential_height: Latitude is not between -90 and 90 degrees: ',lat
      return
   endif

   latr = lat * (pi2)        ! in radians
   sin2  = sin(latr) * sin(latr)
   termg = grav_equator * ( (1.d0+somigliana*sin2) / &
           sqrt(1.d0-eccentricity*eccentricity*sin2) )
   termr = semi_major_axis / (1.d0 + flattening + grav_ratio - 2.d0*flattening*sin2)

   compute_geopotential_height = (termg/grav)*((termr*H)/(termr+H))
   !print '(A,3f15.5)','compute_geopotential_height:',&
   !                   H,compute_geopotential_height,H-compute_geopotential_height

end function compute_geopotential_height

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
