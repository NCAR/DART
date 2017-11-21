! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_aviso

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! convert_aviso - program that reads a netCDF file from CMEMS/Aviso+
!                 containing L3 (along-track) sea level anomalies from
!                 any of three supported platforms and writes a DART
!                 observation sequence file.
!
! http://marine.copernicus.eu/services-portfolio/access-to-products/
! product: SEALEVEL_GLO_SLA_L3_REP_OBSERVATIONS_008_018
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8

use     utilities_mod, only : nc_check, initialize_utilities, finalize_utilities, &
                              error_handler, E_ERR, E_MSG, find_namelist_in_file, &
                              check_namelist_read, do_output, logfileunit

use  time_manager_mod, only : time_type, set_calendar_type, set_date, set_time, &
                              increment_time, get_time, operator(+), GREGORIAN

use      location_mod, only : VERTISSURFACE

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data

use      obs_kind_mod, only : AltiKa_SEA_SURFACE_ANOMALY, &     ! AltiKa
                              Cryosat2_SEA_SURFACE_ANOMALY, &   ! Cryosat-2
                              Envisat_SEA_SURFACE_ANOMALY, &    ! Envisat
                              Envisatn_SEA_SURFACE_ANOMALY, &   ! Envisat Extension
                              ERS1_SEA_SURFACE_ANOMALY, &       ! ERS-1
                              ERS2_SEA_SURFACE_ANOMALY, &       ! ERS-2
                              GFO_SEA_SURFACE_ANOMALY, &        ! Geosat Follow On
                              Haiyang2A_SEA_SURFACE_ANOMALY, &  ! Haiyang-2A
                              J1_SEA_SURFACE_ANOMALY, &         ! Jason-1
                              J1g_SEA_SURFACE_ANOMALY, &        ! Jason-1 Geodetic Phase
                              J1n_SEA_SURFACE_ANOMALY, &        ! Jason-1 New Orbit
                              J2_SEA_SURFACE_ANOMALY, &         ! OSTM/Jason-2
                              TP_SEA_SURFACE_ANOMALY, &         ! Topex/Poseidon
                              TPn_SEA_SURFACE_ANOMALY           ! Topex/Poseidon New Orbit

use          sort_mod, only : index_sort

use obs_utilities_mod, only : getvar_real, add_obs_to_seq, &
                              create_3d_obs, getvar_int, getdimlen, &
                              getvar_int_2d, query_varname, &
                              get_short_as_r8, get_int_as_r8

use           netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, parameter :: use_input_qc = .false.

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

integer  :: iunit, io
integer  :: ncid, varid, ios
integer  :: ntime, n, i, oday, osec, nused
integer  :: iyear, imonth, iday, ihour, imin, isec
logical  :: first_obs

real(r8) :: qc, missing

real(r8), allocatable :: lat(:), lon(:), sla(:), RE(:), time(:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time

character(len=100) :: platform, attvalue

character(len=256) :: RE_file
character(len=256) :: profiler_out_file = 'obs_seq.aviso'
character(len=256) :: input_file = '/glade/p/work/fredc/OBS/SLA_alongtrack/dt_global_j1_sla_vfec_20080101_20140106.nc'

character(len=512) :: string1, string2, string3

! namelist with default settings

real(r8) :: observation_error_std = 0.03_r8

namelist / convert_aviso_nml / observation_error_std 

!-----------------------------------------------------------------------
! start of executable code
!-----------------------------------------------------------------------

call initialize_utilities('convert_aviso')

! read command line argument(s) - the second is optional.
call getarg(1, input_file)
call getarg(2, RE_file)

! Read the namelist.
call find_namelist_in_file('input.nml', 'convert_aviso_nml', iunit)
read(iunit, nml = convert_aviso_nml, iostat = io)
call check_namelist_read(iunit, io, 'convert_aviso_nml')

! Record the namelist values used for the run.
call error_handler(E_MSG,'convert_aviso:','convert_aviso_nml values are',' ',' ',' ')
if (do_output()) write(logfileunit, nml=convert_aviso_nml)
if (do_output()) write(     *     , nml=convert_aviso_nml)

if (input_file == '') then
   write(string1,*)'..  Require a command-line argument specifying the input file.'
   write(string2,*)'Must be a fully-qualified file name.'
   write(string3,*)'USAGE: convert_aviso  <input_file_name>'
   call error_handler(E_ERR, 'convert_aviso', string1, &
          source, revision, revdate, text2=string2, text3=string3)
endif

write(string1,*)'Converting input file "'//trim(input_file)//'"'
call error_handler(E_MSG, 'convert_aviso', string1)

first_obs = .true.

call nc_check( nf90_open(input_file, nf90_nowrite, ncid), &
               'convert_aviso', 'opening file "'//trim(input_file)//'"')

call nc_check( nf90_get_att(ncid, NF90_GLOBAL, "platform", platform), &
               'convert_aviso', 'reading global attribute platform')

call getdimlen(ncid, "time", ntime)

allocate( lat(ntime))
allocate( lon(ntime))
allocate( sla(ntime))
allocate(time(ntime))

call get_int_as_r8(  ncid, "latitude",  lat)
call get_int_as_r8(  ncid, "longitude", lon)
call get_short_as_r8(ncid, "SLA",       sla, missing)

!>@todo should check to see that the RE in the file is the same time
!>  as the data, the same length, etc.

allocate(RE(ntime))

if (RE_file == '') then
   write(*,*)'..  Use a constant representativeness error from namelist.'
   RE(:) = observation_error_std
else
   write(*,*)'..  Use a specified representativeness error.'
   call nc_check( nf90_open(RE_file, nf90_nowrite, ncid), &
                  'convert_aviso', 'opening file '//trim(RE_file)) 
   call get_short_as_r8(ncid, "RE", RE, missing)
endif 

! Convert the time to DART type

ios = nf90_inq_varid(ncid, 'time', varid)
call nc_check(ios, 'convert_aviso', 'inq_varid "time"')

ios = nf90_get_var(ncid, varid, time)
call nc_check(ios, 'convert_aviso', 'get_var "time"')

ios = nf90_get_att(ncid, varid, 'units', attvalue)
call nc_check(ios, 'convert_aviso', 'get_att time "units"')

! time:units = "days since 1950-01-01 00:00:00" ;
!               1234567890

if (attvalue(1:10) /= 'days since') then
   write(string1,*)'expecting time units of [days since ... ]'
   write(string2,*)'read time units of ['//trim(attvalue)//']'
   call error_handler(E_ERR, 'convert_aviso', string1, &
          source, revision, revdate, text2=string2)
endif

read(attvalue,'(11x,i4,5(1x,i2))',iostat=ios)iyear,imonth,iday,ihour,imin,isec
if (ios /= 0) then
   write(string1,*)'Unable to read time units. Error status was ',ios
   write(string2,*)'expected "days since YYYY-MM-DD HH:MM:SS"'
   write(string3,*)'was      "'//trim(attvalue)//'"'
   call error_handler(E_ERR, 'convert_aviso', string1, &
          source, revision, revdate, text2=string2, text3=string3)
endif

call set_calendar_type(GREGORIAN)
comp_day0 = set_date(iyear, imonth, iday, ihour, imin, isec)

call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

! create a new one
call init_obs_sequence(obs_seq, num_copies, num_qc, ntime)
do i = 1, num_copies
  call set_copy_meta_data(obs_seq, i, 'observation')
end do
do i = 1, num_qc
  call set_qc_meta_data(obs_seq, i, 'Data QC')
end do

! Set the DART data quality control.  Be consistent with NCEP codes;
! 0 is 'must use', 1 is good, no reason not to use it.
qc = 1.0_r8

! add sla component data to obs_seq
nused = 0
staloop: do n = 1, ntime

  if ( sla(n) == missing)                          cycle staloop
  if ( lat(n) >  90.0_r8 .or. lat(n) <  -90.0_r8 ) cycle staloop
  if ( lon(n) > 360.0_r8 .or. lon(n) <    0.0_r8 ) cycle staloop

  ! compute time of observation
  oday     = int(time(n))
  osec     = nint((time(n) - real(oday,r8))*86400.0_r8)
  time_obs = set_time(osec, oday) + comp_day0
  call get_time(time_obs, osec, oday)

  select case (platform)

    case("AltiKa")
      call create_3d_obs(lat(n), lon(n), 0.0_r8, VERTISSURFACE, sla(n), &
                       AltiKa_SEA_SURFACE_ANOMALY, RE(n), oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
    case("Cryosat-2")
      call create_3d_obs(lat(n), lon(n), 0.0_r8, VERTISSURFACE, sla(n), &
                       Cryosat2_SEA_SURFACE_ANOMALY, RE(n), oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
    case("ERS-1")
      call create_3d_obs(lat(n), lon(n), 0.0_r8, VERTISSURFACE, sla(n), &
                       ERS1_SEA_SURFACE_ANOMALY, RE(n), oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
    case("ERS-2")
      call create_3d_obs(lat(n), lon(n), 0.0_r8, VERTISSURFACE, sla(n), &
                       ERS2_SEA_SURFACE_ANOMALY, RE(n), oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
    case("Envisat")
      call create_3d_obs(lat(n), lon(n), 0.0_r8, VERTISSURFACE, sla(n), &
                       Envisat_SEA_SURFACE_ANOMALY, RE(n), oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
    case("Envisat Extension")
      call create_3d_obs(lat(n), lon(n), 0.0_r8, VERTISSURFACE, sla(n), &
                       Envisatn_SEA_SURFACE_ANOMALY, RE(n), oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
    case("Geosat Follow On")
      call create_3d_obs(lat(n), lon(n), 0.0_r8, VERTISSURFACE, sla(n), &
                       GFO_SEA_SURFACE_ANOMALY, RE(n), oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
    case("Haiyang-2A")
      call create_3d_obs(lat(n), lon(n), 0.0_r8, VERTISSURFACE, sla(n), &
                       Haiyang2A_SEA_SURFACE_ANOMALY, RE(n), oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
    case("Jason-1")
      call create_3d_obs(lat(n), lon(n), 0.0_r8, VERTISSURFACE, sla(n), &
                       J1_SEA_SURFACE_ANOMALY, RE(n), oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
    case("Jason-1 Geodetic Phase")
      call create_3d_obs(lat(n), lon(n), 0.0_r8, VERTISSURFACE, sla(n), &
                       J1g_SEA_SURFACE_ANOMALY, RE(n), oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
    case("Jason-1 New Orbit")
      call create_3d_obs(lat(n), lon(n), 0.0_r8, VERTISSURFACE, sla(n), &
                       J1n_SEA_SURFACE_ANOMALY, RE(n), oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
    case("OSTM/Jason-2")
      call create_3d_obs(lat(n), lon(n), 0.0_r8, VERTISSURFACE, sla(n), &
                       J2_SEA_SURFACE_ANOMALY, RE(n), oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
    case("Topex/Poseidon")
      call create_3d_obs(lat(n), lon(n), 0.0_r8, VERTISSURFACE, sla(n), &
                       TP_SEA_SURFACE_ANOMALY, RE(n), oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
    case("Topex/Poseidon New Orbit")
      call create_3d_obs(lat(n), lon(n), 0.0_r8, VERTISSURFACE, sla(n), &
                       TPn_SEA_SURFACE_ANOMALY, RE(n), oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

    case default
      string1 = 'Unknown "platform" from "'//trim(input_file)//'"'
      string2 = 'read as "'//trim(platform)//'"'
      string3 = 'must be "AltiKa", "Cryosat-2", "ERS-1", "ERS-2", "Envisat", "Envisat Extension", &
                 "Geosat Follow On", "Haiyang-2A", "Jason-1", "Jason-1 New Orbit", "Jason-1 Geodetic Phase", &
                 "OSTM/Jason-2", "Topex/Poseidon", or "Topex/Poseidon New Orbit"'
      call error_handler(E_ERR, 'query_varname', string1, &
                  source, revision, revdate, text2=string2, text3=string3)
  end select

  nused = nused + 1

end do staloop

! need to wait to close file because in the loop it queries the
! report types.
call nc_check(nf90_close(ncid),'convert_aviso','closing "'//trim(input_file)//'"')

! if we added any obs to the sequence, write it now.
if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, profiler_out_file)

call finalize_utilities()

end program convert_aviso

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
