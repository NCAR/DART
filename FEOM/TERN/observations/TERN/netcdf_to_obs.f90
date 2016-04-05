! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program netcdf_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  netcdf_to_obs - input is a netCDF file with the following ncdump:
!
!  float Latitude(Latitude) ;
!          Latitude:units = "degrees_north" ;
!          Latitude:scale_factor = 1.f ;
!          Latitude:add_offset = 0.f ;
!          Latitude:long_name = "Latitude" ;
!          Latitude:_CoordinateAxisType = "Lat" ;
!
!  float Longitude(Longitude) ;
!          Longitude:units = "degrees_east" ;
!          Longitude:scale_factor = 1.f ;
!          Longitude:add_offset = 0.f ;
!          Longitude:long_name = "Longitude" ;
!          Longitude:_CoordinateAxisType = "Lon" ;
!
!  short soil_moisture_c(Longitude, Latitude) ;
!          soil_moisture_c:long_name = "Volumetric Soil Moisture from C-band" ;
!          soil_moisture_c:units = "percent" ;
!          soil_moisture_c:coordinates = "Latitude Longitude" ;
!          soil_moisture_c:_FillValue = -32767s ;
!          soil_moisture_c:scale_factor = 1.f ;
!          soil_moisture_c:add_offset = 0.f ;
!
!  short sm_c_error(Longitude, Latitude) ;
!          sm_c_error:long_name = "Uncertainty of Soil moisture in C-band";
!          sm_c_error:units = "percent" ;
!          sm_c_error:coordinates = "Latitude Longitude" ;
!          sm_c_error:_FillValue = -32767s ;
!          sm_c_error:scale_factor = 0.01f ;
!          sm_c_error:add_offset = 0.f ;
!
!  created  15 Feb 2014   Tim Hoar NCAR/IMAGe
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, PI, DEG2RAD, MISSING_R8

use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              find_namelist_in_file, check_namelist_read,  &
                              do_nml_file, do_nml_term, nmlfileunit, logfileunit, &
                              nc_check, error_handler, E_ERR

use  time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                              set_date, get_time, &
                              print_date, print_time

use      location_mod, only : location_type, is_location_in_region, &
                              set_location, VERTISHEIGHT, VERTISUNDEF

use  obs_sequence_mod, only : obs_sequence_type, obs_type, &
                              static_init_obs_sequence, init_obs,            &
                              write_obs_seq, init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data

use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq, getdimlen

use      obs_kind_mod, only : SOIL_MOISTURE

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! Namelist with default values
!-----------------------------------------------------------------------

character(len=256) :: input_file  = 'LPRM-AMSR_E_L3_D_SOILM3_V002_20020619.nc'
character(len=256) :: output_file = 'obs_seq.out'
real(r8)           :: min_lat     = -90.0_r8
real(r8)           :: max_lat     =  90.0_r8
real(r8)           :: min_lon     =   0.0_r8
real(r8)           :: max_lon     = 360.0_r8
logical :: verbose = .false.  ! set to .true. for more print info

namelist /netcdf_to_obs_nml/ input_file, output_file, &
            min_lat, max_lat, min_lon, max_lon, verbose

!-----------------------------------------------------------------------
! globally-scoped variables
!-----------------------------------------------------------------------

integer :: nlon, nlat, ilon, ilat
real(r8), allocatable, dimension(:)   :: lat, lon
real(r8), allocatable, dimension(:,:) :: datmat
real(r8), allocatable, dimension(:,:) :: obserr

logical :: first_obs
integer :: i, oday, osec, iunit, io
integer :: num_copies, num_qc, max_obs, ncid, VarID

real(r8) :: qc, vert, obs_value, obs_error

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time
type(location_type)     :: obs_loc, minl, maxl

character(len=256) :: string1, string2, string3

!-----------------------------------------------------------------------
! start of executable code
!-----------------------------------------------------------------------

call initialize_utilities('netcdf_to_obs')

call find_namelist_in_file('input.nml', 'netcdf_to_obs_nml', iunit)
read(iunit, nml = netcdf_to_obs_nml, iostat = io)
call check_namelist_read(iunit, io, 'netcdf_to_obs_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=netcdf_to_obs_nml)
if (do_nml_term()) write(     *     , nml=netcdf_to_obs_nml)

! time setup
call set_calendar_type(GREGORIAN)

! Decode time from the file name ... they have the form:
! LPRM-AMSR_E_L3_D_SOILM3_V002_20020619.nc

time_obs = decode_filename(adjustl(input_file))

! Figure out how many variables are going to be read ...

! open netcdf file
! get dimensions, read metadata arrays
call nc_check( nf90_open(input_file, nf90_nowrite, ncid), &
               'netcdf_to_obs', 'opening file '//trim(input_file))

call getdimlen(ncid, 'Latitude',  nlat)
call getdimlen(ncid, 'Longitude', nlon)

allocate(lat(nlat),lon(nlon))
allocate(datmat(nlat, nlon))
allocate(obserr(nlat, nlon))

call nc_check(nf90_inq_VarID(ncid,'Latitude',VarID),'netcdf_to_obs','inq_varid Latitude')
call nc_check(nf90_get_var(ncid,VarID,lat),         'netcdf_to_obs','get_var   Latitude')

call nc_check(nf90_inq_VarID(ncid,'Longitude',VarID),'netcdf_to_obs','inq_varid Longitude')
call nc_check(nf90_get_var(ncid,VarID,lon),          'netcdf_to_obs','get_var   Longitude')

! convert from -180/180 to 0/360 if necessary
where (lon < 0.0_r8) lon = lon + 360.0_r8

! FIXME ... bulletproof the input ... min < max, for example
! Set boundaries for geographic selection.
! Supports crossing the dateline.
minl = set_location(min_lon, min_lat, 0.0_r8, VERTISUNDEF)
maxl = set_location(max_lon, max_lat, 0.0_r8, VERTISUNDEF)

if (verbose) then
   write(     *     ,*)
   write(     *     ,*)'Applying geographic subsetting of '
   write(     *     ,*)'minimum location is (lat,lon) ',min_lat, min_lon
   write(     *     ,*)'maximum location is (lat,lon) ',max_lat, max_lon
   write(logfileunit,*)
   write(logfileunit,*)'Applying geographic subsetting of '
   write(logfileunit,*)'minimum location is (lat,lon) ',min_lat, min_lon
   write(logfileunit,*)'maximum location is (lat,lon) ',max_lat, max_lon
endif

!----------------------------------------------------------------------
! read the data for the requested arrays.
!       float Soil_Moisture_from_C_band(Xtrack, Track) ;
!                Soil_Moisture_from_C_band:long_name = "AMSRv05 Soil Moisture (6.9GHz)" ;
!                Soil_Moisture_from_C_band:units = "[m3 m-3]" ;
!                Soil_Moisture_from_C_band:_FillValue = -9999.f ;
!                Soil_Moisture_from_C_band:actual_range = -99.f, 99.f ;
!                Soil_Moisture_from_C_band:coordinates = "lat lon" ;
!                Soil_Moisture_from_C_band:grid_mapping = "wgs84" ;
!                Soil_Moisture_from_C_band:gain = "0.01" ;
!        float Soil_Moisture_Error_from_C_band(Xtrack, Track) ;
!                Soil_Moisture_Error_from_C_band:long_name = "AMSRv05 Soil Moisture Error (6.9GHz)" ;
!                Soil_Moisture_Error_from_C_band:units = "[m3 m-3]" ;
!                Soil_Moisture_Error_from_C_band:_FillValue = -9999.f ;
!                Soil_Moisture_Error_from_C_band:actual_range = -99.f, 50.f ;
!                Soil_Moisture_Error_from_C_band:coordinates = "lat lon" ;
!                Soil_Moisture_Error_from_C_band:grid_mapping = "wgs84" ;
!                Soil_Moisture_Error_from_C_band:gain = "0.01" ;
!
!----------------------------------------------------------------------

! call get_variable(input_file, ncid, 'soil_moisture_c', datmat)
! call get_variable(input_file, ncid,      'sm_c_error', obserr)
call get_variable(input_file, ncid, 'Soil_Moisture_from_C_band', datmat)
call get_variable(input_file, ncid, 'Soil_Moisture_Error_from_C_band', obserr)

! extract time of observation into gregorian day, sec.
call get_time(time_obs, osec, oday)

! Depth of the observation is a couple centimeters
vert = 0.05_r8

! each observation in this series will have a single observation value
! and a quality control flag.  the max possible number of obs needs to
! be specified but it will only write out the actual number created.
max_obs    = nlat * nlon
num_copies = 1
num_qc     = 1

! call the initialization code, and initialize two empty observation types
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
first_obs = .true.

! create a new, empty obs_seq file.  you must give a max limit
! on number of obs.  increase the size if too small.
call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)

! the first one needs to contain the string 'observation' and the
! second needs the string 'QC'.
call set_copy_meta_data(obs_seq, 1, 'observation')
call   set_qc_meta_data(obs_seq, 1, 'Data QC')

! we have to assign a quality control flag
qc = 0.0_r8

Longitude: do ilon = 1, nlon
Latitude:  do ilat = 1, nlat

   if (datmat(ilat,ilon) == MISSING_R8) cycle Latitude

   obs_loc = set_location(lon(ilon), lat(ilat), vert, VERTISUNDEF)

   if (.not. is_location_in_region(obs_loc, minl, maxl)) cycle Latitude

   ! FIXME ... this comes in some percent ... NOAH is m^3/m^3, CLM is kg/m^2
   obs_value = datmat(ilat,ilon) ! FIXME ... convert to
   obs_error = obserr(ilat,ilon) ! FIXME ... units

   ! make an obs derived type, and then add it to the sequence
   call create_3d_obs(lat(ilat), lon(ilon), vert, VERTISHEIGHT, obs_value, &
                               SOIL_MOISTURE, obs_error, oday, osec, qc, obs)
   call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

enddo Latitude
enddo Longitude

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (verbose) then
      write(     *     ,*)'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
      write(logfileunit,*)'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   endif
   call write_obs_seq(obs_seq, output_file)
endif

! end of main program
call finalize_utilities()


!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------


subroutine get_variable(filename, ncid, varname, slab)
character(len=*),         intent(in)  :: filename
integer,                  intent(in)  :: ncid
character(len=*),         intent(in)  :: varname
real(r8), dimension(:,:), intent(out) :: slab

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimlens
character(len=NF90_MAX_NAME)          :: dimname
character(len=NF90_MAX_NAME)          :: dimnames(NF90_MAX_VAR_DIMS)
integer                               :: numdims, xtype, dimlen
character(len=256)                    :: long_name, units
real(r8)                              :: scale_factor, add_offset, FillValue
logical                               :: filled

string3 = trim(filename)//' '//trim(varname)

call nc_check( nf90_inq_VarID(ncid, trim(varname), VarID), &
               'netcdf_to_obs', 'inquire var '// trim(string3))

call nc_check(nf90_inquire_variable(ncid,VarID,dimids=dimIDs,ndims=numdims,xtype=xtype), &
               'netcdf_to_obs', 'inquire '//trim(string3))

if (numdims /= 2) then
   write(string1,*)'variable ['//trim(varname)//'] is unsupported shape.'
   write(string2,*)'Can only support 2D variables, this is ',numdims
   call error_handler(E_ERR, 'netcdf_to_obs', string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif

dimlen  = 1
DimensionLoop : do i = 1,numdims
   write(string1,'(''inquire dimension'',i2,A)') i,trim(string3)
   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), name=dimname, len=dimlen), &
                                       'netcdf_to_obs', string1)
   dimlens( i) = dimlen
   dimnames(i) = dimname
enddo DimensionLoop

if ((dimlens(1) /= nlat) .or. (dimlens(2) /= nlon)) then
   write(string1,*)'variable ['//trim(varname)//'] is unsupported size.'
   write(string2,*)'Expected ',nlat,' by ',nlon,' ... got ...  ',dimlens(1),' by ',dimlens(2)
   call error_handler(E_ERR, 'netcdf_to_obs', string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif

! If the long_name and/or units attributes are set, get them.
! They are not REQUIRED to exist but are nice to use if they are present.
long_name = varname
units     = 'unknown'
scale_factor = 1.0_r8
add_offset   = 0.0_r8
FillValue    = 0.0_r8
filled       = .false.

if( nf90_inquire_attribute(    ncid, VarID, 'long_name') == NF90_NOERR ) then
   call nc_check( nf90_get_att(ncid, VarID, 'long_name', long_name), &
               'static_init_model', 'get_att long_name '//trim(string3))
endif

if( nf90_inquire_attribute(    ncid, VarID, 'units') == NF90_NOERR )  then
   call nc_check( nf90_get_att(ncid, VarID, 'units', units), &
               'static_init_model', 'get_att units '//trim(string3))
endif

if( nf90_inquire_attribute(    ncid, VarID, 'scale_factor') == NF90_NOERR )  then
   call nc_check( nf90_get_att(ncid, VarID, 'scale_factor', scale_factor), &
               'static_init_model', 'get_att scale_factor '//trim(string3))
endif

if( nf90_inquire_attribute(    ncid, VarID, 'add_offset') == NF90_NOERR )  then
   call nc_check( nf90_get_att(ncid, VarID, 'add_offset', add_offset), &
               'static_init_model', 'get_att add_offset '//trim(string3))
endif

if( nf90_inquire_attribute(    ncid, VarID, '_FillValue') == NF90_NOERR )  then
   call nc_check( nf90_get_att(ncid, VarID, '_FillValue', FillValue), &
               'static_init_model', 'get_att _FillValue '//trim(string3))
   filled = .true.
endif

call nc_check( nf90_get_var(ncid, VarID, slab), &
               'netcdf_to_obs', 'getting var '// trim(varname))

! Apply the scale/offset/fill settings
! Check what happens with the fill value

if (filled) then
   where (slab /= FillValue) slab = slab * scale_factor + add_offset
   where (slab == FillValue) slab = MISSING_R8
else
   slab = slab * scale_factor + add_offset
endif

end subroutine get_variable



function decode_filename(filename)
! Extract the date from the filename when the filenames are like:
! LPRM-AMSR_E_L3_D_SOILM3_V002_20020619.nc

character(len=*), intent(in) :: filename
type(time_type) :: decode_filename

integer :: strlen, strstart
integer :: iyear, imonth, iday, rcio

! Need to find the YYYYMMDD part, right before the '.nc' end.
strlen   = len_trim(filename) - 3
strstart = strlen - 7

read(filename(strstart:strlen),'(i4,i2,i2)',iostat=rcio)iyear,imonth,iday

if (rcio /= 0 ) then
   write(string1,*)'Unable to parse date from input filename.'
   write(string2,*)'Trying to read YYYYMMDD from <'//filename(strstart:strlen)//'>'
   call error_handler(E_ERR, 'netcdf_to_obs', string1, &
              source, revision, revdate, text2=string2)
endif

decode_filename = set_date(iyear,imonth,iday,0,0,0)

call print_date(decode_filename,'Date from filename is ')
call print_time(decode_filename,'Time from filename is ')
call print_date(decode_filename,'Date from filename is ',logfileunit)
call print_time(decode_filename,'Time from filename is ',logfileunit)

end function decode_filename


end program netcdf_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
