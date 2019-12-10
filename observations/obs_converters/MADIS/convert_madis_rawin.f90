! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_madis_rawin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_madis_rawin - program that reads a netCDF file from the 
!                         MADIS database that contains rawinsonde data 
!                         and writes a DART obs_seq file using the DART 
!                         library routines.
!
!     created Dec. 2007 Ryan Torn, NCAR/MMM
!     modified Dec. 2008 Soyoung Ha and David Dowell, NCAR/MMM
!     - added dewpoint as an output variable
!     - added relative humidity as an output variable
!
!     modified to use a common set of utilities, better netcdf error checks,
!     able to insert obs with any time correctly (not only monotonically
!     increasing times)    nancy collins,  ncar/image   11 march 2010
!     
!     keep original obs times, make source for all converters as similar
!     as possbile.   nancy collins,  ncar/image   26 march 2010
! 
!     add code to use the (new?) pressure levels for vertical coord in the
!     significant wind obs.  also make an optional namelist section; if
!     enabled, the selection of mandatory/significant/both and height/pressure
!     can be from a namelist and avoid the prompts and reads from the console.
!     nancy collins, ncar/image, 23 mar 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, missing_r8
use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              find_namelist_in_file, check_namelist_read,         &
                              do_nml_file, do_nml_term, nmlfileunit
use  netcdf_utilities_mod, only : nc_open_file_readonly, nc_close_file, nc_check
use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                                  increment_time, get_time, operator(-), GREGORIAN
use      location_mod, only : VERTISSURFACE, VERTISPRESSURE, VERTISHEIGHT
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, & 
                              init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data
use        meteor_mod, only : sat_vapor_pressure, specific_humidity, & 
                              wind_dirspd_to_uv, pres_alt_to_pres, &
                              temp_and_dewpoint_to_rh
use       obs_err_mod, only : rawin_temp_error, rawin_wind_error, &
                              rawin_pres_error, rawin_rel_hum_error
use dewpoint_obs_err_mod, only : dewpt_error_from_rh_and_temp, &
                                 rh_error_from_dewpt_and_temp
use obs_def_altimeter_mod, only : compute_altimeter
use          obs_kind_mod, only : RADIOSONDE_U_WIND_COMPONENT,      & 
                                  RADIOSONDE_V_WIND_COMPONENT,      & 
                                  RADIOSONDE_TEMPERATURE,           & 
                                  RADIOSONDE_SPECIFIC_HUMIDITY,     &
                                  RADIOSONDE_RELATIVE_HUMIDITY,     &
                                  RADIOSONDE_DEWPOINT,              &
                                  RADIOSONDE_SURFACE_ALTIMETER 
use             sort_mod,  only : index_sort
use     obs_utilities_mod, only : add_obs_to_seq, create_3d_obs, &
                                  getdimlen, getvar_int, set_missing_name, &
                                  get_or_fill_int

use           netcdf

implicit none

character(len=19),  parameter :: rawin_in_file  = 'rawin_input.nc'
character(len=129), parameter :: rawin_out_file = 'obs_seq.rawin'

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

integer :: oday, osec, nsig, nsound, nused, io, iunit, &
           nmaxml, nmaxsw, nmaxswp, nmaxswh, nmaxst, maxobs, nvars_man, nvars_sigt, k, n, i, ncid

integer,  allocatable :: obscnt(:), used(:), sorted_used(:), tused(:), nman(:)
real(r8), allocatable :: lon(:), lat(:), elev(:), otime(:)

logical :: fexist, first_obs

real(r8) :: uwnd, vwnd, qobs, qsat, dptk, oerr, &
            pres_miss, wdir_miss, wspd_miss, tair_miss, tdew_miss, prespa, & 
            qc, altim, rh, qerr, hght_miss
real(r8), parameter   :: HGHT_THRESHOLD = 40000.0_r8   

real(r8), allocatable :: pres(:), wdir(:), wspd(:), tair(:), tdew(:), hght(:)
integer,  allocatable :: qc_pres(:), qc_wdir(:), qc_wspd(:), qc_tair(:), qc_tdew(:), qc_hght(:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time


! the following logical parameters control which water-vapor variables appear in the output file,
! whether to use the NCEP error or Lin and Hubbard (2004) moisture error model, and if the
! input file has data quality control fields, whether to use or ignore them.  more recent
! MADIS files contain the significant winds in both height and pressure vertical coordinates
! and the last variable controls which one is used.
! if you are using the namelist (see 'use_namelist' below) you can change these at runtime.
! if not, they are compile-time settings.
logical :: LH_err                    = .false.
logical :: include_specific_humidity = .true.
logical :: include_relative_humidity = .false.
logical :: include_dewpoint          = .false.
logical :: use_input_qc              = .true. 
logical :: wind_use_vert_pressure    = .true.  ! use what you find.  if both, use this one.
logical :: verbose                   = .false. ! get more output

! mandatory levels are always converted.  significant levels are controlled
! by these vars.  if you are using the namelist they are set by the
! namelist in input.nml.  if not, you are prompted at runtime for T/F responses.
logical :: do_significant_level_temps = .true.
logical :: do_significant_level_winds = .true.
character(len=32) :: sig_wind_type = 'unknown'

logical :: use_namelist = .false.

! THIS WILL NOT BE READ IN UNLESS use_namelist IS SET TO .true.
! if true, the code will NOT prompt for input from the console.

namelist /convert_madis_rawin_nml/ &
   do_significant_level_temps, do_significant_level_winds,  &
   wind_use_vert_pressure, LH_err, include_specific_humidity, &
   include_relative_humidity, include_dewpoint, use_input_qc, &
   verbose

!-------- start of executable code -----------

call initialize_utilities('convert_madis_rawin')

if (use_namelist) then
   ! Read the namelist entry
   call find_namelist_in_file("input.nml", "convert_madis_rawin_nml", iunit)
   read(iunit, nml = convert_madis_rawin_nml, iostat = io)
   call check_namelist_read(iunit, io, "convert_madis_rawin_nml")
   
   ! Record the namelist values used for the run ...
   if (do_nml_file()) write(nmlfileunit, nml=convert_madis_rawin_nml)
   if (do_nml_term()) write(     *     , nml=convert_madis_rawin_nml)
else
   ! prompt for optional significant level values
   print*,'Include significant level winds, temperature?: '
   read*, do_significant_level_winds, do_significant_level_temps
endif

! put the reference date into DART format
call set_calendar_type(GREGORIAN)
comp_day0 = set_date(1970, 1, 1, 0, 0, 0)

first_obs = .true.


ncid = nc_open_file_readonly(rawin_in_file, 'convert_madis_rawin')

call getdimlen(ncid, "recNum", nsound)
call set_missing_name("missing_value")

if (verbose) print *, 'total number of soundings = ', nsound

allocate(obscnt(nsound))  
allocate(lon(nsound), lat(nsound), elev(nsound), otime(nsound))
allocate(nman(nsound), used(nsound), sorted_used(nsound), tused(nsound))

! apparently the code uses 99999 as missing data.
! don't include that in the counts.  get all the values
! so you can print them if asked, and then set the optional
! ones to 0 if not requested.

! mandatory levels
call getvar_int(ncid, "numMand", obscnt)
nmaxml = maxval(obscnt, mask=obscnt<1000)
if (verbose) print *, 'max mandatory levels per sounding, over all soundings = ', nmaxml

! significant level temps
call getvar_int(ncid, "numSigT", obscnt)
nmaxst = maxval(obscnt, mask=obscnt<1000)
if (verbose) print *, 'max significant level temperatures per sounding, over all soundings = ', nmaxst

if ( .not. do_significant_level_temps ) nmaxst = 0

! significant level winds
! slightly more complex in that they could be on pressure levels,
! heights, or both.  find the max of either.

! this call fills the array with 0s if not present
call get_or_fill_int(ncid, "numSigPresW", obscnt)
nmaxswp = maxval(obscnt, mask=obscnt<1000)
if (verbose) print *, 'max significant level winds per sounding, pressure vertical, over all soundings = ', nmaxswp

call get_or_fill_int(ncid, "numSigW", obscnt)
nmaxswh = maxval(obscnt, mask=obscnt<1000)
if (verbose) print *, 'max significant level winds per sounding, height vertical, over all soundings = ', nmaxswh

if ( .not. do_significant_level_winds ) then
  nmaxsw = 0
else
  nmaxsw = max(nmaxswp, nmaxswh)
endif



! compute the largest possible number of obs that can be created from this input

nvars_man = 4
nvars_sigt = 1
if (include_specific_humidity) then
  nvars_man  = nvars_man  + 1
  nvars_sigt = nvars_sigt + 1
endif
if (include_relative_humidity) then
  nvars_man  = nvars_man  + 1
  nvars_sigt = nvars_sigt + 1
endif
if (include_dewpoint) then
  nvars_man  = nvars_man  + 1
  nvars_sigt = nvars_sigt + 1
endif

maxobs = nsound * (nvars_man * nmaxml + 2 * nmaxsw + nvars_sigt * nmaxst + 1)
deallocate(obscnt)

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

inquire(file=rawin_out_file, exist=fexist)
if ( fexist ) then

  ! existing file found, append to it
  call read_obs_seq(rawin_out_file, 0, 0, maxobs, obs_seq)
  
  if (verbose) print *, 'existing file found, appending to it'

else

  ! create a new one
  call init_obs_sequence(obs_seq, num_copies, num_qc, maxobs)
  do n = 1, num_copies
    call set_copy_meta_data(obs_seq, n, 'MADIS observation')
  enddo
  do n = 1, num_qc
    call set_qc_meta_data(obs_seq, n, 'Data QC')
  enddo

  if (verbose) print *, 'creating new obs_seq file'
endif

! Set the DART data quality control.  Be consistent with NCEP codes;
! 0 is 'must use', 1 is good, no reason not to use it.
qc = 1.0_r8

nused = 0
sondeloop1 : do n = 1, nsound !  loop over all soundings in the file 

  call getvar_real_1d_1val(ncid, "staLat",  n, lat(n)  )
  call getvar_real_1d_1val(ncid, "staLon",  n, lon(n)  )
  call getvar_real_1d_1val(ncid, "staElev", n, elev(n) )
  call getvar_int_1d_1val (ncid, "numMand", n, nman(n) )
  call getvar_real_1d_1val(ncid, "synTime", n, otime(n))  !! , time_miss)
  ! the original code had a line to get the fill value here but it
  ! was commented out.  is there one?  do we need it?

  if (nman(n) < 0 .or. nman(n) > nmaxml) cycle sondeloop1

  if ( otime(n) < 0.0_r8 ) cycle sondeloop1

  ! check the lat/lon values to see if they are ok
  if ( lat(n) >  90.0_r8 .or. lat(n) <  -90.0_r8 ) cycle sondeloop1
  if ( lon(n) > 180.0_r8 .or. lon(n) < -180.0_r8 ) cycle sondeloop1

  ! Check for duplicate soundings
  do i = 1, nused
    if ( lon(n) == lon(used(i)) .and. &
         lat(n) == lat(used(i)) ) cycle sondeloop1
  enddo

  nused = nused + 1
  used(nused) = n
  tused(nused) = nint(otime(n))

enddo sondeloop1

! sort by obs time.  this will make it much faster to create
! the output file because it won't have to search the linked list
! for the right insertion point.

call index_sort(tused, sorted_used, nused)


sondeloop2: do i = 1, nused

  ! get the next unique observation in sorted time order
  n = used(sorted_used(i))

  ! change lon from -180 to 180 into 0-360
  if (lon(n) < 0.0_r8) lon(n) = lon(n) + 360.0_r8

  ! compute time of observation
  time_obs = increment_time(comp_day0, nint(otime(n)))

  ! extract actual time of observation in file into oday, osec.
  call get_time(time_obs, osec, oday)

  if (nman(n) > 0) then
    allocate(pres(nman(n)))  ;  allocate(tair(nman(n)))  ;  allocate(tdew(nman(n)))
    allocate(wdir(nman(n)))  ;  allocate(wspd(nman(n)))
  
    allocate(qc_pres(nman(n)))  ;  allocate(qc_tair(nman(n)))  ;  allocate(qc_tdew(nman(n)))
    allocate(qc_wdir(nman(n)))  ;  allocate(qc_wspd(nman(n)))
  
    call getvar_real_2d(ncid, "prMan", n, nman(n), pres, pres_miss)
    call getvar_real_2d(ncid, "tpMan", n, nman(n), tair, tair_miss)
    call getvar_real_2d(ncid, "tdMan", n, nman(n), tdew, tdew_miss)
    call getvar_real_2d(ncid, "wdMan", n, nman(n), wdir, wdir_miss)
    call getvar_real_2d(ncid, "wsMan", n, nman(n), wspd, wspd_miss)
  
    ! if user says to use QC, read them in or fill if not there
    if (use_input_qc) then
       call get_or_fill_QC_2d(ncid, "prManQCR", n, nman(n), qc_pres)
       call get_or_fill_QC_2d(ncid, "tpManQCR", n, nman(n), qc_tair)
       call get_or_fill_QC_2d(ncid, "tdManQCR", n, nman(n), qc_tdew)
       call get_or_fill_QC_2d(ncid, "wdManQCR", n, nman(n), qc_wdir)
       call get_or_fill_QC_2d(ncid, "wsManQCR", n, nman(n), qc_wspd)
    else
       qc_pres = 0
       qc_tair = 0 ;  qc_tdew = 0
       qc_wdir = 0 ;  qc_wspd = 0
    endif
  
    if ( pres(1) /= pres_miss .and. qc_pres(1) == 0 .and. elev(n) < 9999.0) then
  
      altim = compute_altimeter(pres(1), elev(n))
      oerr  = rawin_pres_error(pres_alt_to_pres(elev(n)) * 0.01_r8)
      if ( altim >= 890.0_r8 .and. altim <= 1100.0_r8 .and. oerr /= missing_r8 ) then
  
        call create_3d_obs(lat(n), lon(n), elev(n), VERTISSURFACE, altim, &
                           RADIOSONDE_SURFACE_ALTIMETER, oerr, oday, osec, qc, obs)
        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
  
      endif
  
    endif
  
    do k = 2, nman(n)   ! obtain the mandatory level data
  
      prespa = pres(k) * 100.0_r8
  
      if ( wdir(k) /= wdir_miss .and. qc_wdir(k) == 0 .and. &
           wspd(k) /= wspd_miss .and. qc_wspd(k) == 0  ) then
  
        call wind_dirspd_to_uv(wdir(k), wspd(k), uwnd, vwnd)
        oerr = rawin_wind_error(pres(k))
        if ( abs(uwnd) <= 150.0_r8 .and. & 
             abs(vwnd) <= 150.0_r8 .and. oerr /= missing_r8 ) then
  
          call create_3d_obs(lat(n), lon(n), prespa, VERTISPRESSURE, uwnd, &
                             RADIOSONDE_U_WIND_COMPONENT, oerr, oday, osec, qc, obs)
          call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
   
          call create_3d_obs(lat(n), lon(n), prespa, VERTISPRESSURE, vwnd, &
                             RADIOSONDE_V_WIND_COMPONENT, oerr, oday, osec, qc, obs)
          call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
  
        endif
  
      endif
  
      if ( tair(k) /= tair_miss .and. qc_tair(k) == 0 ) then
  
        oerr = rawin_temp_error(pres(k))
        if ( tair(k) >= 180.0_r8 .and. &
             tair(k) <= 330.0_r8 .and. oerr /= missing_r8 ) then
  
          call create_3d_obs(lat(n), lon(n), prespa, VERTISPRESSURE, tair(k), &
                             RADIOSONDE_TEMPERATURE, oerr, oday, osec, qc, obs)
          call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
  
        endif
  
      endif
    
      ! if the air and dewpoint obs are both ok, then see which of the possible
      ! three types of moisture obs to generate.
      if ( tair(k) /= tair_miss .and. qc_tair(k) == 0 .and. &
           tdew(k) /= tdew_miss .and. qc_tdew(k) == 0  ) then
  
        ! before we start computing things based on the dewpoint,
        ! make sure it isn't larger than the air temp.  if it is
        ! more than a degree larger, skip it completely.  if it is
        ! less, set them equal and continue.
        if (tdew(k) > tair(k)) then
           if (tdew(k) > tair(k) + 1.0_r8) goto 100
           tdew(k) = tair(k)
        endif

        ! tdew is the dewpoint depression
        dptk = tair(k) - tdew(k)
  
        if ( include_specific_humidity ) then
  
          qobs = specific_humidity(sat_vapor_pressure(dptk),    prespa)
          qsat = specific_humidity(sat_vapor_pressure(tair(k)), prespa)
          if ( LH_err ) then
            qerr = rh_error_from_dewpt_and_temp(tair(k), dptk)
          else
            qerr = rawin_rel_hum_error(pres(k), tair(k), qobs / qsat)
          endif
          oerr = max(qerr * qsat, 0.0001_r8)
    
          if ( qobs >  0.0_r8  .and. &
               qobs <= 0.07_r8 .and. qerr /= missing_r8 ) then
    
            call create_3d_obs(lat(n), lon(n), prespa, VERTISPRESSURE, qobs, &
                               RADIOSONDE_SPECIFIC_HUMIDITY, oerr, oday, osec, qc, obs)
            call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
    
          endif
    
        endif
    
        if ( include_relative_humidity ) then
    
          rh = temp_and_dewpoint_to_rh(tair(k), dptk)
          if ( LH_err ) then
            oerr = rh_error_from_dewpt_and_temp(tair(k), dptk)
          else
            oerr = rawin_rel_hum_error(pres(k), tair(k), rh)
          endif
    
          if ( rh >  0.0_r8 .and. &
               rh <= 1.5_r8 .and. oerr /= missing_r8 ) then
  
            call create_3d_obs(lat(n), lon(n), prespa, VERTISPRESSURE, rh, &
                               RADIOSONDE_RELATIVE_HUMIDITY, oerr, oday, osec, qc, obs)
            call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
          endif
  
        endif
    
        if ( include_dewpoint ) then
  
          rh = temp_and_dewpoint_to_rh(tair(k), dptk)
          oerr = dewpt_error_from_rh_and_temp(tair(k), rh)
    
          if ( rh >  0.0_r8 .and. &
               rh <= 1.5_r8 .and. oerr /= missing_r8 ) then
  
            call create_3d_obs(lat(n), lon(n), prespa, VERTISPRESSURE, dptk, &
                               RADIOSONDE_DEWPOINT, oerr, oday, osec, qc, obs)
            call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
  
          endif
  
        endif
  
      endif  ! quality control/missing check on tair, tdew
  
100 continue

    enddo

    deallocate(pres, wdir, wspd, tair, tdew, qc_pres, qc_wdir, qc_wspd, qc_tair, qc_tdew)
  endif


  !  If desired, read the significant-level temperature data, write to obs_seq
  !  skip soundings with nsig of 99999
  call getvar_int_1d_1val(ncid, "numSigT", n, nsig )

  if ( do_significant_level_temps .and. nsig > 0 .and. nsig <= nmaxst) then

    allocate(pres(nsig))     ;  allocate(tair(nsig))     ;  allocate(tdew(nsig))
    allocate(qc_pres(nsig))  ;  allocate(qc_tair(nsig))  ;  allocate(qc_tdew(nsig))

    !  read significant level data
    call getvar_real_2d(ncid, "prSigT", n, nsig, pres, pres_miss)
    call getvar_real_2d(ncid, "tpSigT", n, nsig, tair, tair_miss)
    call getvar_real_2d(ncid, "tdSigT", n, nsig, tdew, tdew_miss)

    if (use_input_qc) then
       call get_or_fill_QC_2d(ncid, "prSigTQCR", n, nsig, qc_pres)
       call get_or_fill_QC_2d(ncid, "tpSigTQCR", n, nsig, qc_tair)
       call get_or_fill_QC_2d(ncid, "tdSigTQCR", n, nsig, qc_tdew)
    else
       qc_pres = 0
       qc_tair = 0
       qc_tdew = 0
    endif
  
    do k = 1, nsig

      prespa = pres(k) * 100.0_r8

      if ( tair(k) /= tair_miss .and. qc_tair(k) == 0 ) then
  
        oerr = rawin_temp_error(pres(k))
        if ( tair(k) >= 180.0_r8 .and. &
             tair(k) <= 330.0_r8 .and. oerr /= missing_r8 ) then

          call create_3d_obs(lat(n), lon(n), prespa, VERTISPRESSURE, tair(k), &
                             RADIOSONDE_TEMPERATURE, oerr, oday, osec, qc, obs)
          call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
  
        endif
    
      endif
  
      ! if the air and dewpoint obs are both ok, then see which of the possible
      ! three types of moisture obs to generate.
      if ( tair(k) /= tair_miss .and. qc_tair(k) == 0 .and. &
           tdew(k) /= tdew_miss .and. qc_tdew(k) == 0  ) then

        ! before we start computing things based on the dewpoint,
        ! make sure it isn't larger than the air temp.  if it is
        ! more than a degree larger, skip it completely.  if it is
        ! less, set them equal and continue.
        if (tdew(k) > tair(k)) then
           if (tdew(k) > tair(k) + 1.0_r8) goto 200
           tdew(k) = tair(k)
        endif

        ! tdew is the dewpoint depression
        dptk = tair(k) - tdew(k)

        if ( include_specific_humidity ) then

          qobs = specific_humidity(sat_vapor_pressure(dptk),    prespa)
          qsat = specific_humidity(sat_vapor_pressure(tair(k)), prespa)
          if ( LH_err ) then
            qerr = rh_error_from_dewpt_and_temp(tair(k), dptk)
          else
            qerr = rawin_rel_hum_error(pres(k), tair(k), qobs / qsat)
          endif
          oerr = max(qerr * qsat, 0.0001_r8)
          if ( qobs >  0.0_r8  .and. &
               qobs <= 0.07_r8 .and. qerr /= missing_r8 ) then
  
            call create_3d_obs(lat(n), lon(n), prespa, VERTISPRESSURE, qobs, &
                               RADIOSONDE_SPECIFIC_HUMIDITY, oerr, oday, osec, qc, obs)
            call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
  
          endif
  
        endif
  
        if ( include_relative_humidity ) then

          rh = temp_and_dewpoint_to_rh(tair(k), dptk)
          if ( LH_err ) then
            oerr = rh_error_from_dewpt_and_temp(tair(k), dptk)
          else
            oerr = rawin_rel_hum_error(pres(k), tair(k), rh)
          endif
          
          if ( rh >  0.0_r8 .and. &
               rh <= 1.5_r8 .and. oerr /= missing_r8 ) then

            call create_3d_obs(lat(n), lon(n), prespa, VERTISPRESSURE, rh, &
                               RADIOSONDE_RELATIVE_HUMIDITY, oerr, oday, osec, qc, obs)
            call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
  
          endif

        endif
  
        if ( include_dewpoint ) then
  
          rh = temp_and_dewpoint_to_rh(tair(k), dptk)
          oerr = dewpt_error_from_rh_and_temp(tair(k), rh)
  
          if ( rh >  0.0_r8 .and. &
               rh <= 1.5_r8 .and. oerr /= missing_r8 ) then

            call create_3d_obs(lat(n), lon(n), prespa, VERTISPRESSURE, dptk, &
                               RADIOSONDE_DEWPOINT, oerr, oday, osec, qc, obs)
            call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
  
          endif

        endif

      endif  ! quality control/missing check on tair and tdew

200 continue

    enddo
    deallocate(pres, tair, tdew, qc_pres, qc_tair, qc_tdew)

  endif

  !  If desired, read the pressure significant-level wind data, write to obs_seq
  !  Convert whichever is present - with vert units of pressure or height.

  ! select pressure or height, depending on which one is there.
  if ( do_significant_level_winds ) &
    call single_vert_type_sigwinds(ncid, nsound, n, sig_wind_type)


  !  note that pressure vs height are two separate sections because
  !  the netcdf arrays with data have slightly different field names,
  !  and the vertical coordinate is different in the create obs call.

  nsig = 0
  if ( do_significant_level_winds .and. sig_wind_type == "pressure") then
    call getvar_int_1d_1val(ncid, "numSigPresW", n, nsig )
  endif

  !  nsig is 0 if not doing sig winds, and skip soundings with nsig of 99999
  if ( nsig > 0 .and. nsig <= nmaxsw ) then

    allocate(pres(nsig))     ;  allocate(wdir(nsig))     ;  allocate(wspd(nsig))
    allocate(qc_pres(nsig))  ;  allocate(qc_wdir(nsig))  ;  allocate(qc_wspd(nsig))

    !  read significant level data
    call getvar_real_2d(ncid,   "prSigW", n, nsig, pres, pres_miss)
    call getvar_real_2d(ncid, "wdSigPrW", n, nsig, wdir, wdir_miss)
    call getvar_real_2d(ncid, "wsSigPrW", n, nsig, wspd, wspd_miss)

    if (use_input_qc) then
       call get_or_fill_QC_2d(ncid,   "prSigWQCR", n, nsig, qc_pres)
       call get_or_fill_QC_2d(ncid, "wdSigPrWQCR", n, nsig, qc_wdir)
       call get_or_fill_QC_2d(ncid, "wsSigPrWQCR", n, nsig, qc_wspd)
    else
       qc_pres = 0
       qc_wdir = 0
       qc_wspd = 0
    endif
    
    do k = 1, nsig

      prespa = pres(k) * 100.0_r8
 
      !  add data to the observation sequence here.
      if ( wdir(k) /= wdir_miss .and. qc_wdir(k) == 0 .and. &
           wspd(k) /= wspd_miss .and. qc_wspd(k) == 0 .and. &
           pres(k) /= pres_miss .and. qc_pres(k) == 0 ) then


        call wind_dirspd_to_uv(wdir(k), wspd(k), uwnd, vwnd)

        ! expects press in hPa, which is what we already have
        oerr = rawin_wind_error(pres(k))
        if ( abs(uwnd) <= 150.0_r8 .and. &
             abs(vwnd) <= 150.0_r8 .and. oerr /= missing_r8 ) then

          call create_3d_obs(lat(n), lon(n), prespa, VERTISPRESSURE, uwnd, &
                             RADIOSONDE_U_WIND_COMPONENT, oerr, oday, osec, qc, obs)
          call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

          call create_3d_obs(lat(n), lon(n), prespa, VERTISPRESSURE, vwnd, &
                             RADIOSONDE_V_WIND_COMPONENT, oerr, oday, osec, qc, obs)
          call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

        endif

      endif

    enddo

    deallocate(pres, wdir, wspd, qc_pres, qc_wdir, qc_wspd)

  endif

  !  If desired, read the height significant-level wind data, write to obs_seq
  nsig = 0
  if ( do_significant_level_winds .and. sig_wind_type == "height") then
    call getvar_int_1d_1val(ncid, "numSigW", n, nsig )
  endif

  !  nsig is 0 if not doing sig winds, and skip soundings with nsig of 99999
  if ( nsig > 0 .and. nsig <= nmaxsw ) then

    allocate(hght(nsig))     ;  allocate(wdir(nsig))     ;  allocate(wspd(nsig))
    allocate(qc_hght(nsig))  ;  allocate(qc_wdir(nsig))  ;  allocate(qc_wspd(nsig))

    !  read significant level data
    call getvar_real_2d(ncid, "htSigW", n, nsig, hght, hght_miss)
    call getvar_real_2d(ncid, "wdSigW", n, nsig, wdir, wdir_miss)
    call getvar_real_2d(ncid, "wsSigW", n, nsig, wspd, wspd_miss)

    if (use_input_qc) then
       call get_or_fill_QC_2d(ncid, "htSigWQCR", n, nsig, qc_hght)
       call get_or_fill_QC_2d(ncid, "wdSigWQCR", n, nsig, qc_wdir)
       call get_or_fill_QC_2d(ncid, "wsSigWQCR", n, nsig, qc_wspd)
    else
       qc_hght = 0
       qc_wdir = 0
       qc_wspd = 0
    endif

    levelloop: do k = 1, nsig

      !  add data to the observation sequence here.
      if ( wdir(k) /= wdir_miss .and. qc_wdir(k) == 0 .and. &
           wspd(k) /= wspd_miss .and. qc_wspd(k) == 0 .and. &
           hght(k) /= hght_miss .and. qc_hght(k) == 0) then

         ! make sure elevation is a reasonable value if height comes
         ! in as zero and use it instead of 0.0
         if (hght(k) == 0.0) then
            if (elev(n) < 9999.0) then
               hght(k) = elev(n)
            else
               cycle levelloop
            endif
         endif

        call wind_dirspd_to_uv(wdir(k), wspd(k), uwnd, vwnd)
        ! the pres_alt_to_pres() fails above 44km, so give any obs above
        ! this height the same wind error value as 40km.
        oerr = rawin_wind_error(pres_alt_to_pres(min(hght(k),HGHT_THRESHOLD)) * 0.01_r8)
        if ( abs(uwnd) <= 150.0_r8 .and. & 
             abs(vwnd) <= 150.0_r8 .and. oerr /= missing_r8 ) then

          call create_3d_obs(lat(n), lon(n), hght(k), VERTISHEIGHT, uwnd, &
                             RADIOSONDE_U_WIND_COMPONENT, oerr, oday, osec, qc, obs)
          call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

          call create_3d_obs(lat(n), lon(n), hght(k), VERTISHEIGHT, vwnd, &
                             RADIOSONDE_V_WIND_COMPONENT, oerr, oday, osec, qc, obs)
          call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
        endif
      endif

    enddo levelloop

    deallocate(hght, wdir, wspd, qc_hght, qc_wdir, qc_wspd)

  endif

enddo sondeloop2

! have to close at end of loop, unlike other versions of the madis converters
call nc_close_file(ncid, 'convert_madis_rawin')

! if we added any obs to the sequence, write it now.
if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, rawin_out_file)

!end of main program
call finalize_utilities()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! make sure this sounding doesn't have significant level winds on both
! pressure and height.  glen has files where some non-conus radiosondes
! do have both in the same sounding.  we aren't sure if these are disjoint
! obs or replicated obs, so to err on the safe side we're only converting
! one or the other.

subroutine single_vert_type_sigwinds(ncid, max_sound, this_sound, sig_wind_type)
 integer, intent(in) :: ncid
 integer, intent(in) :: max_sound
 integer, intent(in) :: this_sound
 character(len=*), intent(out) :: sig_wind_type

integer, allocatable :: p(:), h(:)
integer :: pcount, hcount
integer :: MAX_CNT = 1000

allocate(p(max_sound), h(max_sound))

! avoid 9999 or 99999 which is used as a fill
call get_or_fill_int(ncid, "numSigPresW", p)
pcount = min(p(this_sound), MAX_CNT)
if (pcount == MAX_CNT) pcount = 0

call get_or_fill_int(ncid, "numSigW", h)
hcount = min(h(this_sound), MAX_CNT)
if (hcount == MAX_CNT) hcount = 0

if (pcount > 0 .and. hcount > 0) then
  if ( wind_use_vert_pressure ) then
    sig_wind_type = 'pressure'
    if (verbose) print *, 'sounding ', this_sound, ' has sigwinds on both; will use pressure levels based on setting of "wind_use_vert_pressure"'
  else
    sig_wind_type = 'height'
    if (verbose) print *, 'sounding ', this_sound, ' has sigwinds on both; will use height levels based on setting of "wind_use_vert_pressure"'
  endif
else if (pcount > 0) then
  sig_wind_type = 'pressure'
  if (verbose) print *, 'sounding ', this_sound, ' has significant level winds on pressure levels'
else if (hcount > 0) then
  sig_wind_type = 'height'
  if (verbose) print *, 'sounding ', this_sound, ' has significant level winds on height levels'
else
  sig_wind_type = 'none'
endif

deallocate(p, h)

end subroutine single_vert_type_sigwinds


! specialized versions of the netcdf get routines that seem to be
! pretty specific to this version of the code, so i didn't put them
! in the general observations utilities file.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   getvar_real_1d_1val - subroutine that inquires, gets the variable, and fills 
!            in the missing value attribute if that arg is present.
!            takes a single start, uses count=1, returns a scalar
!
!      ncid - open netcdf file handle
!      varname - string name of netcdf variable
!      start - starting index in the 1d array
!      dout - output value.  real(r8)
!      dmiss - value that signals a missing value   real(r8), optional
!
!     created 11 Mar 2010,  nancy collins,  ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getvar_real_1d_1val(ncid, varname, start, dout, dmiss)
 integer,            intent(in)   :: ncid
 character(len = *), intent(in)   :: varname
 integer,            intent(in)   :: start
 real(r8),           intent(out)  :: dout
 real(r8), optional, intent(out)  :: dmiss

integer :: varid

! read the data for the requested array, and get the fill value
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'getvar_real', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, dout, start = (/ start /) ), &
               'getvar_real', 'getting var '// trim(varname))

if (present(dmiss)) &
   call nc_check( nf90_get_att(ncid, varid, '_FillValue', dmiss), &
               'getvar_real', 'getting attr "_FillValue" for '//trim(varname))

end subroutine getvar_real_1d_1val

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   getvar_int_1d_1val - subroutine that inquires, gets the variable, and fills 
!            in the missing value attribute if that arg is present.
!            takes a single start, uses count=1, returns a scalar
!
!      ncid - open netcdf file handle
!      varname - string name of netcdf variable
!      start - starting index in the 1d array
!      dout - output value.  int
!      dmiss - value that signals a missing value   int, optional
!
!     created 11 Mar 2010,  nancy collins,  ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getvar_int_1d_1val(ncid, varname, start, dout, dmiss)
 integer,            intent(in)   :: ncid
 character(len = *), intent(in)   :: varname
 integer,            intent(in)   :: start
 integer,            intent(out)  :: dout
 integer,  optional, intent(out)  :: dmiss

integer :: varid

! read the data for the requested array, and get the fill value
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'getvar_int_1d_1val', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, dout, start = (/ start /) ), &
               'getvar_int_1d_1val', 'getting var '// trim(varname))

if (present(dmiss)) &
   call nc_check( nf90_get_att(ncid, varid, '_FillValue', dmiss), &
               'getvar_int_1d_1val', 'getting attr "_FillValue" for '//trim(varname))

end subroutine getvar_int_1d_1val

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   getvar_real_2d - subroutine that inquires, gets the variable, and fills 
!            in the missing value attribute if that arg is present.
!     SPECIALIZED for this use - assumes start = (/ 1, n /) and count = (/ m, 1 /)
!           so takes a scalar start, count, returns a 1d_array
!
!      ncid - open netcdf file handle
!      varname - string name of netcdf variable
!      start - starting index in the 2d array.  integer
!      count - nitems to get. integer
!      darray - output array.  real(r8)
!      dmiss - value that signals a missing value   real(r8), optional
!
!     created 11 Mar 2010,  nancy collins,  ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getvar_real_2d(ncid, varname, start, count, darray, dmiss)
 integer,            intent(in)   :: ncid
 character(len = *), intent(in)   :: varname
 integer,            intent(in)   :: start
 integer,            intent(in)   :: count
 real(r8),           intent(out)  :: darray(:)
 real(r8), optional, intent(out)  :: dmiss

integer :: varid

! read the data for the requested array, and get the fill value
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'getvar_real_2d', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, darray, &
                start=(/ 1, start /), count=(/ count, 1 /) ), &
               'getvar_real_2d', 'getting var '// trim(varname))

if (present(dmiss)) &
   call nc_check( nf90_get_att(ncid, varid, '_FillValue', dmiss), &
               'getvar_real_2d', 'getting attr "_FillValue" for '//trim(varname))

end subroutine getvar_real_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_or_fill_QC_2d - subroutine which gets the requested netcdf variable
!           but if it isn't there, it fills the array with 0s.  not an
!           error if it's not present.  assumes integer data array
!     SPECIALIZED for this use - assumes start = (/ 1, n /) and count = (/ m, 1 /)
!           so takes a scalar start, count, returns a 1d_array
!           also prints out a message if fill used.
!
!      ncid - open netcdf file handle
!      varname - string name of netcdf variable
!      start - starting index in the 2d array.  integer
!      count - nitems to get. integer
!      darray - output array.  integer
!
!     created Mar 8, 2010    nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_or_fill_QC_2d(ncid, varname, start, count, darray)
 integer,            intent(in)    :: ncid
 character(len = *), intent(in)    :: varname
 integer,            intent(in)    :: start
 integer,            intent(in)    :: count
 integer,            intent(inout) :: darray(:)

integer :: varid, nfrc

! test to see if variable is present.  if yes, read it in.
! otherwise, set to fill value, or 0 if none given.

nfrc = nf90_inq_varid(ncid, varname, varid)
if (nfrc == NF90_NOERR) then
   call nc_check( nf90_get_var(ncid, varid, darray, &
                  start=(/ 1, start /), count=(/ count, 1 /) ), &
                  'get_or_fill_int_2d', 'reading '//trim(varname) )
else
   darray = 0
   if (start == 1 .and. verbose) & 
     print *, 'QC field named ' // trim(varname) // ' was not found in input, 0 used instead'
endif

end subroutine get_or_fill_QC_2d

end program convert_madis_rawin

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
