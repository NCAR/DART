! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program MIDAS_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! MIDAS_to_obs - reads the MIDAS data as created by Alex in a netCDF file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use          types_mod, only : r8, MISSING_R8, metadatalength

use      utilities_mod, only : initialize_utilities, finalize_utilities, &
                               register_module, error_handler, E_ERR, E_MSG, &
                               do_nml_file, do_nml_term, &
                               check_namelist_read, find_namelist_in_file, &
                               nmlfileunit, file_exist

use  netcdf_utilities_mod, only : nc_check

use   time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                               set_time, get_time, print_time, print_date

use       location_mod, only : location_type, VERTISUNDEF

use   obs_sequence_mod, only : obs_sequence_type, obs_type, &
                               static_init_obs_sequence, init_obs, write_obs_seq, &
                               init_obs_sequence, get_num_obs, &
                               set_copy_meta_data, set_qc_meta_data

use        obs_def_mod, only : obs_def_type

use  obs_utilities_mod, only : add_obs_to_seq, create_3d_obs

use       obs_kind_mod, only : MIDAS_TEC

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

character(len=256) :: input_file    = 'infile.nc'
character(len=256) :: obs_out_file  = 'obs_seq.out'
logical            :: verbose       = .false.

namelist /MIDAS_to_obs_nml/ input_file, obs_out_file, verbose

!-----------------------------------------------------------------------
! globally-scoped variables
!-----------------------------------------------------------------------

character(len=256)      :: string1
integer                 :: itime
logical                 :: first_obs
integer                 :: oday, osec, rcio, iunit
integer                 :: num_copies, num_qc, max_obs
real(r8)                :: qc
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: prev_time, time_obs

integer :: nlon, nlat, ntimes
real(r8), allocatable, dimension(:)   :: latitude
real(r8), allocatable, dimension(:)   :: longitude
real(r8), allocatable, dimension(:)   :: time
real(r8), allocatable, dimension(:,:) :: TEC
real(r8), allocatable, dimension(:,:) :: ObsErrVar

integer  :: ncid, ilon, ilat

!-----------------------------------------------------------------------
! start of executable code
!-----------------------------------------------------------------------

call initialize_utilities('MIDAS_to_obs')

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "MIDAS_to_obs_nml", iunit)
read(iunit, nml = MIDAS_to_obs_nml, iostat = rcio)
call check_namelist_read(iunit, rcio, "MIDAS_to_obs_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=MIDAS_to_obs_nml)
if (do_nml_term()) write(     *     , nml=MIDAS_to_obs_nml)

! time setup
call set_calendar_type(GREGORIAN)
prev_time = set_time(0, 0)

! Read the basic MIDAS netCDF information
ntimes   = read_midas_metadata(input_file)

num_copies = 1
num_qc     = 1
first_obs  = .true.

max_obs = ntimes * nlon * nlat
allocate(TEC(nlon,nlat),ObsErrVar(nlon,nlat))

call static_init_obs_sequence()
call init_obs(        obs,      num_copies, num_qc)
call init_obs(        prev_obs, num_copies, num_qc)
call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)

! the first one needs to contain the string 'observation' and the
! second needs the string 'QC'.
call set_copy_meta_data(obs_seq, 1, 'observation')
call set_qc_meta_data(  obs_seq, 1, 'MIDAS QC')

call nc_check(nf90_open(input_file, NF90_NOWRITE, ncid), &
                   'main', 'open '//trim(input_file))

! Actually get the observations for each timestep
obsloop: do itime = 1,ntimes

   TEC       = MISSING_R8 ! just in case the get_slab fails
   ObsErrVar = MISSING_R8 ! just in case the get_slab fails

   call get_slab(ncid, itime,      'TEC', TEC, time_obs)
   call get_slab(ncid, itime, 'Variance', ObsErrVar)

   if (verbose) call print_date(time_obs, 'obs time is')

   call get_time(time_obs, osec, oday)

   ! make an obs derived type, and then add it to the sequence
   ! If the QC value is good, use the observation.
   ! Increasingly larger QC values are more questionable quality data.

   qc   = 0.0_r8   ! all observations are assumed to be GREAT, i.e. 0.0

   do ilon = 1,nlon
   do ilat = 1,nlat

      if (verbose) then
         write(*,*)ilon, ilat, latitude(ilat), longitude(ilon), &
                   TEC(ilon,ilat), ObsErrVar(ilon,ilat)
      endif

      call create_3d_obs(latitude(ilat), longitude(ilon), 0.0_r8, VERTISUNDEF, &
           TEC(ilon,ilat), MIDAS_TEC, ObsErrVar(ilon,ilat), oday, osec, qc, obs)

      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
   enddo
   enddo

enddo obsloop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (verbose) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
else
   write(string1,*)'There are no observations to write.'
   call error_handler(E_MSG,'MIDAS_to_obs',string1,source,revision,revdate)
endif

! end of main program
call finalize_utilities()

contains



function read_midas_metadata(input_file)
!----------------------------------------------------------------------------
! Read the list of parameters for every site we know about and
! return the number of sites we know about.

integer                      :: read_midas_metadata
character(len=*), intent(in) :: input_file

integer :: ncid, VarID
integer :: time_dimid
integer :: lon_dimid
integer :: lat_dimid

! Check to make sure the required parameter file exists

if ( .not. file_exist(input_file) ) then
   write(string1,*) 'MIDAS data file [', trim(input_file),'] does not exist.'
   call error_handler(E_ERR,'read_midas_metadata',string1,source,revision,revdate)
endif

call nc_check(nf90_open(input_file, NF90_NOWRITE, ncid), &
                   'read_midas_metadata', 'open '//trim(input_file))

! Read the dimensions for everything of interest

call nc_check(nf90_inq_dimid(ncid, 'time', time_dimid), &
                  'read_midas_metadata','inq_dimid time '//trim(input_file))
call nc_check(nf90_inquire_dimension(ncid, time_dimid, len=read_midas_metadata), &
                  'read_midas_metadata','inquire_dimension time '//trim(input_file))

call nc_check(nf90_inq_dimid(ncid, 'latitude', lat_dimid), &
                  'read_midas_metadata','inq_dimid latitude '//trim(input_file))
call nc_check(nf90_inquire_dimension(ncid, lat_dimid, len=nlat), &
                  'read_midas_metadata','inquire_dimension latitude '//trim(input_file))

call nc_check(nf90_inq_dimid(ncid, 'longitude', lon_dimid), &
                  'read_midas_metadata','inq_dimid longitude '//trim(input_file))
call nc_check(nf90_inquire_dimension(ncid, lon_dimid, len=nlon), &
                  'read_midas_metadata','inquire_dimension longitude '//trim(input_file))

allocate(longitude(nlon), latitude(nlat), time(read_midas_metadata))

! Read the coordinate variables

call nc_check(nf90_inq_varid(ncid, 'longitude', VarID), 'read_midas_metadata', 'inq_varid longitude')
call nc_check(nf90_get_var(ncid, VarID, longitude), 'get_var longitude')

where(longitude  <   0.0_r8) longitude = longitude + 360.0_r8
where(longitude == 360.0_r8) longitude = 0.0_r8

call nc_check(nf90_inq_varid(ncid, 'latitude', VarID), 'read_midas_metadata', 'inq_varid latitude')
call nc_check(nf90_get_var(ncid, VarID, latitude), 'get_var latitude')

call nc_check(nf90_inq_varid(ncid, 'time', VarID), 'read_midas_metadata', 'inq_varid time')
call nc_check(nf90_get_var(ncid, VarID, time), 'get_var time')

call nc_check(nf90_close(ncid), 'read_midas_metadata', 'close '//trim(input_file))

if (verbose) then
   write(*,*)
   write(*,*)'There are ',read_midas_metadata,' timesteps in file '//trim(input_file)
   write(*,*)'longitude ',longitude 
   write(*,*)'latitude  ',latitude
endif

end function read_midas_metadata



subroutine get_slab(ncid, itime, varname, datmat, time_obs)
!----------------------------------------------------------------------------
! Read an entire 2D slab for a single timestep from an open netCDF file.
! We know the variables in the netCDF file do not have a scale/offset.

integer,                   intent(in)  :: ncid
integer,                   intent(in)  :: itime
character(len=*),          intent(in)  :: varname
real(r8), dimension(:,:),  intent(out) :: datmat
type(time_type), optional, intent(out) :: time_obs

integer :: VarID, numdims
integer :: i, days, seconds

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
integer, dimension(NF90_MAX_VAR_DIMS) :: mystart, mycount
character(len   =  NF90_MAX_NAME)     :: dimname

mystart(:) = 1
mycount(:) = 1

! Check to make sure the required variable exists
call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), &
         'get_slab', 'inq_varid '//trim(varname))

! Construct the hyperslabbing indices independent of the storage order

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
         'get_slab', 'inquire '//trim(varname))

DimensionLoop : do i = 1,numdims

   write(string1,'(''inquire dimension'',i2,1x,A)') i,trim(varname)
   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), name=dimname), &
          'get_slab', string1)

   if (trim(dimname) == 'Time')      mystart(i) = itime
   if (trim(dimname) == 'time')      mystart(i) = itime
   if (trim(dimname) == 'latitude')  mycount(i) = nlat
   if (trim(dimname) == 'longitude') mycount(i) = nlon

enddo DimensionLoop

! Get hyperslab from file

call nc_check(nf90_get_var(ncid, VarID, datmat, start=mystart(1:numdims), &
         count=mycount(1:numdims)), 'get_slab', 'get_var '//trim(varname))

! Get the time relating to this timestep if needed

if (present(time_obs)) then

   days     = time(itime) ! check for integer truncation
   seconds  = nint((time(itime) - real(days,r8)) * 86400.0_r8)
   time_obs = set_time(seconds,days)

   if (verbose) then
      call print_date(time_obs, 'observation date is')
      call print_time(time_obs, 'observation time is')
   endif

endif

if (verbose) write(*,'(''Read slab index '',i4,'' for variable '',A)')itime,trim(varname)

end subroutine get_slab

!----------------------------------------------------------------------------

end program MIDAS_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
