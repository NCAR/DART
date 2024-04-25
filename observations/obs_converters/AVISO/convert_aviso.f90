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

use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              error_handler, E_ERR, E_MSG

use  netcdf_utilities_mod, only : nc_check

use  time_manager_mod, only : time_type, set_calendar_type, set_date, set_time, &
                              increment_time, get_time, operator(+), GREGORIAN

use      location_mod, only : VERTISSURFACE

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, & 
                              init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data

use      obs_kind_mod, only : J1_SEA_SURFACE_ANOMALY, &  ! Jason-1
                              EN_SEA_SURFACE_ANOMALY, &  ! Envisat
                              GFO_SEA_SURFACE_ANOMALY    ! Geosat Follow On

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

integer  :: ncid, varid, ios
integer  :: ntime, n, i, oday, osec, nused, idx
integer  :: iyear, imonth, iday, ihour, imin, isec
logical  :: first_obs

real(r8) :: oerr, qc, missing

real(r8), allocatable :: lat(:), lon(:), sla(:), RE(:), time(:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time

character(len=100) :: platform, attvalue

character(len=256) :: RE_file
character(len=256) :: profiler_out_file = 'obs_seq.aviso'
character(len=256) :: input_file = '/glade/p/work/fredc/OBS/SLA_alongtrack/dt_global_j1_sla_vfec_20080101_20140106.nc'

character(len=512) :: string1, string2, string3

!-----------------------------------------------------------------------
! start of executable code
!-----------------------------------------------------------------------

call initialize_utilities('convert_aviso')

! command line argument
call GET_COMMAND_ARGUMENT(1, input_file)

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

idx = index(input_file,'sla')
RE_file = input_file(1:idx-1) // 'RE' // input_file(idx+3:len(input_file))

call nc_check( nf90_open(RE_file, nf90_nowrite, ncid), &
               'convert_aviso', 'opening file '//trim(RE_file))

allocate(RE(ntime))

call get_short_as_r8(ncid, "RE", RE, missing)

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

  oerr = 0.03_r8  ! this works mostly for surface to 100mb

  select case (platform)

    case("Jason-1")
      call create_3d_obs(lat(n), lon(n), 0.0_r8, VERTISSURFACE, sla(n), &
                       J1_SEA_SURFACE_ANOMALY, RE(n), oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
    case("Envisat")
      call create_3d_obs(lat(n), lon(n), 0.0_r8, VERTISSURFACE, sla(n), &
                       EN_SEA_SURFACE_ANOMALY, RE(n), oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
    case("Geosat Follow On")
      call create_3d_obs(lat(n), lon(n), 0.0_r8, VERTISSURFACE, sla(n), &
                       GFO_SEA_SURFACE_ANOMALY, RE(n), oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
    case default
      string1 = 'Unknown "platform" from "'//trim(input_file)//'"'
      string2 = 'read as "'//trim(platform)//'"'
      string3 = 'must be "Jason-1", "Envisat", or "Geosat Follow On"'
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
