module dart_obs_seq_mod

use kinds, only: r_kind, r_single, i_kind
use params
use mpisetup

use obs_utilities_mod, only : add_obs_to_seq
use  obs_sequence_mod, only : obs_type, obs_sequence_type, init_obs_sequence, insert_obs_in_seq, &
                              set_copy_meta_data, set_qc_meta_data, write_obs_seq, assignment(=), &
                              init_obs, static_init_obs_sequence, set_obs_def, set_obs_values, set_qc, &
                              destroy_obs, destroy_obs_sequence
use       obs_def_mod, only : set_obs_def_location, set_obs_def_error_variance, &
                              set_obs_def_type_of_obs, set_obs_def_time, set_obs_def_key, &
                              obs_def_type, set_obs_def_external_FO, destroy_obs_def
use   obs_def_gps_mod, only : set_gpsro_ref
use         types_mod, only : obstypelength
use      obs_kind_mod
use      location_mod, only : location_type, set_location, VERTISSURFACE, VERTISPRESSURE, VERTISHEIGHT
use  time_manager_mod, only : time_type, set_date, set_time, set_calendar_type, GREGORIAN, &
                              increment_time, decrement_time, get_time
use           radinfo, only : nuchan
use     utilities_mod, only : to_upper, error_handler, E_ERR, E_MSG

implicit none
private 

! Public subroutine
public :: dart_obs_seq, set_debug

character(len=512) :: string1
logical :: debug = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

subroutine dart_obs_seq (datestring,                              &
                         nobs_conv, nobs_oz, nobs_sat, nobstot,   &
                         ens_size, obs_seq_out_filename, nobs_start, nobs_end)
   implicit none

   character(len=10),     intent(in) :: datestring
   integer,               intent(in) :: nobs_conv, nobs_oz, nobs_sat, nobstot
   integer,               intent(in) :: ens_size
   character(len = 129),  intent(in) :: obs_seq_out_filename
   integer,               intent(in) :: nobs_start, nobs_end

   logical                 :: is_conv, is_sat

   integer                 :: num_obs, num_obs_last, num_copies, num_qc, max_num_obs=8000000
   type(obs_sequence_type) :: seq
   type(obs_type)          :: obs, prev_obs
   type(obs_def_type)      :: obs_def
   type(location_type)     :: location
   type(time_type)         :: time_obs, prev_time
   integer                 :: i, k, n, which_vert, ie, obskind, obs_quantity
   character(len = 129)    :: copy_meta_data, qc_meta_data
   logical                 :: first_obs
   integer                 :: year, month, day, hour, days, seconds
   integer                 :: gps_key = 0
   real(r_kind)            :: gnx, gny, gnz, ds, htop, rfict
   real(r_kind)            :: lat, lon, vloc, obsv, oerr
   real(r_kind)            :: this_mean, this_stdv, this_variance, test_val
   real(r_kind)            :: qc_val(2)
   real(r_kind)            :: missing_val  = -888888.000
   real(r_kind),allocatable  :: ens_copy(:)

   logical                 :: write_external
   logical                 :: has_external = .true.

   ! Number of obs to process
   num_obs = nobs_end - nobs_start + 1

! add back the unbiascorrected mean that was extracted in enkf
   do i = 1, ens_size
      anal_ob(i,:) = anal_ob(i,:) + ensmean_obnobc(:)
!     anal_ob_chunk(i,:) = anal_ob_chunk(i,:) + ensmean_obnobc(nobs_start:nobs_end)  ! CSS for scatterv
   end do

   ! Figure out the number of copies/metadata to output
   if ( write_prior_copies ) then
      num_copies = 5 + 2*ens_size! Obs,prior/post mean, prior/post spread, plus ensemble of priors and posteriors
      allocate(ens_copy(num_copies))
      num_qc     = 2   ! Number of QC--1 for input, the 2nd for DART QC flag, which is bogus, but needed so other DART programs can work
   else
      num_copies = 1                 ! Just 1 true copy for the observation. Priors are metadata, NOT copies
      allocate(ens_copy(ens_size+2)) ! Ensemble of priors, plus obs, plus a mean value
      num_qc     = 1   ! Number of QC--just 1, coming from GSI
   endif

   ! Initialize some values
   ens_copy(:) = missing_val

! Initialize an obs_sequence structure
   call init_obs_sequence(seq, num_copies, num_qc, num_obs)  ! Last entry can be > num_obs (e.g., max_num_obs) but that wastes memory

!  Set QC meta data
   qc_meta_data = 'GSI Quality Control'
   call set_qc_meta_data(seq, 1, qc_meta_data)

!  Set copy meta data
   copy_meta_data = 'observation'
   call set_copy_meta_data(seq, 1, copy_meta_data)

   ! Set more copy and QC meta data if write_prior_copies
   if ( write_prior_copies ) then
      copy_meta_data = 'prior ensemble mean'       ; call set_copy_meta_data(seq, 2, copy_meta_data)
      copy_meta_data = 'posterior ensemble mean'   ; call set_copy_meta_data(seq, 3, copy_meta_data)
      copy_meta_data = 'prior ensemble spread'     ; call set_copy_meta_data(seq, 4, copy_meta_data)
      copy_meta_data = 'posterior ensemble spread' ; call set_copy_meta_data(seq, 5, copy_meta_data)
      do ie=1,ens_size
         write(copy_meta_data,'(a21,1x,i6)') 'prior ensemble member', ie
         call set_copy_meta_data(seq,5+ie*2-1,copy_meta_data)
         write(copy_meta_data,'(a25,1x,i6)') 'posterior ensemble member', ie
         call set_copy_meta_data(seq,5+ie*2,copy_meta_data)
      enddo
      ! Additional QC metadata
      qc_meta_data = 'DART quality control'
      call set_qc_meta_data(seq, 2, qc_meta_data)
   endif

! initialize the obs variable
   call init_obs(obs, num_copies, num_qc)
   call init_obs(prev_obs, num_copies, num_qc)

! set observation time type
   call set_calendar_type(GREGORIAN)

   first_obs = .true.
   obsloop: do i = nobs_start, nobs_end

      ! Initialize QC values inside the loop, because we might modify these values later if modify_dart_qc_flag_for_big_ob_error = true
      qc_val(1)  = 1.0  ! For Input GSI QC
      qc_val(2)  = 0.0  ! For DART QC, not used if write_prior_copies == .false.

      ! Find out if the ob is a conventional ob or radiance
      if ( i > nobs_conv + nobs_oz ) then
         is_sat  = .true.  ; is_conv = .false.
         if ( .not. convert_sat ) cycle obsloop
      else
         is_sat  = .false. ; is_conv = .true.
         if ( .not. convert_conv ) cycle obsloop
      endif

      ! here, obskind is really obs_type, obs_quantity is the generic quantity
      if ( is_conv ) then
         call prepbufr_to_dart_obs_kind(obtype(i),stattype(i),obskind,which_vert,obs_quantity)
      else if ( is_sat .and. convert_sat ) then
         call radiance_to_dart_obs_kind(obtype(i),nuchan(indxsat(i-nobs_conv-nobs_oz)),obskind,which_vert,obs_quantity) ! indxsat has size "nobs_sat"
      end if

      ! skip if unsupported by prepbufr_to_dart_obs_kind or radiance_to_dart_obs_kind
      if ( obskind < 0 ) cycle obsloop

      ! time info
      read(datestring, '(i4,3i2)') year, month, day, hour
      if ( obtime(i) >= 0.0 ) then
         call get_time(increment_time(set_date(year, month, day, hour), &
                                      nint(obtime(i))*3600, 0),         &
                                      seconds, days)
      else if ( obtime(i) < 0.0 ) then
         call get_time(decrement_time(set_date(year, month, day, hour), &
                                      nint(abs(obtime(i)))*3600, 0),    &
                                      seconds, days)
      end if
      time_obs = set_time(seconds, days)

      if ( lie_about_ob_times ) then
        !To make DART assimilate observations in order obs are provided, rather
        !than in time-order, need to lie about the observation times.
        ! When writing out more than 1 obs_seq file, they need to be combined in proper order, so
         ! make obs time dependent on which processor we're on ("nproc")
         call get_time(increment_time(set_date(year, month, day, hour), &
                                      nproc, 0),         &
                                      seconds, days)  ! 
         time_obs = set_time(seconds, days)
      endif 

      !transfer to r_kind (r8) variables
      vloc = obpress(i)*100.0 !Convert to Pa
      lon  = obloclon(i)
      lat  = obloclat(i)
      oerr = oberrvar(i)  ! Already ob error variances

      ens_copy(1) = ob(i) ! Observation value

      !Recenter radiance or surface pressure prior perturbations (which don't
      !include bias-corrections) about the ensemble mean
      !prior H(x_bar), which does include BC.
      ! Really only need to do for obs with bias-correction, like surface
      ! pressure or radiance, but doesn't hurt to do for all obs
      if ( recenter_about_mean_prior ) then ! Default is true
         anal_ob(1:ens_size,i) = anal_ob(1:ens_size,i) + (ensmean_ob(i)-ensmean_obnobc(i))
        !anal_ob_chunk(1:ens_size,i) = anal_ob_chunk(1:ens_size,i) + (ensmean_ob(i)-ensmean_obnobc(i))  ! CSS scatterv
      endif

      ! modify observations and observation errors if necessary
      if( obs_quantity == QTY_SURFACE_PRESSURE) then
         ens_copy(1) = ens_copy(1)*100.0 ! Get into Pa
         oerr = oerr * 10000.0             ! Get into Pa^2--this is variance
         anal_ob(1:ens_size,i) = anal_ob(1:ens_size,i) * 100.        ! Get into Pa
         ensmean_ob(i) = ensmean_ob(i) * 100. ! Get into Pa
        !anal_ob_chunk(1:ens_size,i) = anal_ob_chunk(1:ens_size,i) * 100.        ! CSS scatterv
      endif

      ! Fill array with priors 
      if ( write_prior_copies ) then
         ens_copy(6:num_copies:2) = anal_ob(1:ens_size,i) ! Fill the prior elements
!        ens_copy(6:num_copies:2) = anal_ob_chunk(1:ens_size,i) ! CSS scatterv
         call compute_mean_std(ens_size,ens_copy(6:num_copies:2),ens_copy(2),ens_copy(4)) ! Calculate prior mean/stddev
      else
         ens_copy(2:ens_size+1) = anal_ob(1:ens_size,i)    
!        ens_copy(2:ens_size+1) = anal_ob_chunk(1:ens_size,i)   ! CSS scatterv
         ens_copy(ens_size+2) = ensmean_ob(i)
      endif

      ! modify QC flag if really big observation error
      ! note that obs_sprd_prior(i) also has spread information but it is safer to compute the spread 
      ! from exactly what is going into the obs_seq file
      modify_dart_qc_flag_for_big_ob_error = .false. ! for now this doesn't work properly. so force to false.
      if ( modify_dart_qc_flag_for_big_ob_error ) then
         if ( write_prior_copies ) then
            call compute_mean_std(ens_size,ens_copy(6:num_copies:2),this_mean,this_stdv) ! Calculate prior mean/stddev
         else
            call compute_mean_std(ens_size,ens_copy(2:ens_size+1),this_mean,this_stdv)
         endif
         this_variance = this_stdv * this_stdv ! get into variance
         test_val = variance_coef * this_variance
         if ( oerr .gt. test_val ) then
       !    print *, sqrt(oerr), this_stdv, vloc*0.01
            qc_val(1) = 10.0
            if ( write_prior_copies ) qc_val(2) = 7.0 ! tell DART/obs_diag this observation failed the outlier check so obs_diag won't process
         endif
      endif

      location = set_location(lon, lat, vloc, which_vert)

      call set_obs_def_location(obs_def, location)
      call set_obs_def_type_of_obs(obs_def, obskind)
      call set_obs_def_time(obs_def, time_obs)
      if ( obskind == GPSRO_REFRACTIVITY ) then
         gnx   = 0.0
         gny   = 0.0
         gnz   = 0.0
         ds    = 0.0
         htop  = 0.0
         rfict = 0.0
         call set_gpsro_ref(gps_key, gnx, gny, gnz, rfict, ds, htop, 'GPSREF')
         call set_obs_def_key(obs_def, gps_key)
      end if

      ! Metadata stuff
      if ( .not. write_prior_copies ) then ! No need for metadata if we're writing prior copies
         write_external = write_this_ob_type_external_FO(obskind)
         call set_obs_def_external_FO(obs_def, has_external, write_external,i,ens_size,real(ens_copy(2:ens_size+2),r_kind)) ! CSS sets obs_def%has_exteral_FO
        !call set_obs_def_external_FO(obs_def, has_external, write_external,i,ens_size,real(anal_ob_chunk(1:ens_size,i),r_kind)) ! CSS scatterv
      endif

      ! If not using precomputed obs, GPSRO_REFRACTIVITY must have height
      ! coordinate, but that is unavailable from GSI diag files.  So don't
      ! include GPSRO obs in the obs_sequence file--make user get it from
      ! elsewhere, I guess.
      !if ( .not. write_external ) then
      !  if ( obskind == GPSRO_REFRACTIVITY ) cycle obsloop
      !endif

      call set_obs_def_error_variance(obs_def, oerr)  ! oerr is the ob error variance
      !call set_obs_def_key(obs_def, key)
      call set_obs_def(obs, obs_def)
      call set_obs_values(obs, ens_copy(1:num_copies))
      call set_qc(obs, qc_val(1:num_qc))
      call add_obs_to_seq(seq, obs, time_obs, prev_obs, prev_time, first_obs)

   end do obsloop

   call write_obs_seq(seq, trim(adjustl(obs_seq_out_filename)) )

   call destroy_obs_def(obs_def)  ! CSS added, is this needed?
   call destroy_obs(obs)  
   call destroy_obs_sequence(seq)

   if (allocated(ens_copy)) deallocate(ens_copy)

end subroutine dart_obs_seq

subroutine set_debug(debug_in)
    logical, intent(in) :: debug_in
    debug = debug_in
end subroutine set_debug


!-----------------------------------------------------------------------
! Everything below here is private
!-----------------------------------------------------------------------

subroutine prepbufr_to_dart_obs_kind (obtype, obstype, obs_kind, which_vert, obs_quantity)

! based on DART/observations/NCEP/ascii_to_obs/real_obs_mod.f90

   implicit none
   character(len=20), intent(in)  :: obtype  !'ps','u','v','oz','t','amsua_n15', etc.
   integer,           intent(in)  :: obstype !prepbufr report type
   integer,           intent(out) :: obs_kind
   integer,           intent(out) :: which_vert
   integer,           intent(out) :: obs_quantity

!   assign each observation the correct observation type
!------------------------------------------------------------------------------

   ! make sure we do not fall through the code below without setting
   ! a valid obs kind (e.g. the obstype is one not listed)
   obs_kind = -1
   obs_quantity = -1

   if(obtype(1:3) == 'gps') then
     obs_quantity = QTY_GPSRO
     obs_kind     = GPSRO_REFRACTIVITY
   endif

   if(obtype(1:3) == ' pw') then
      obs_quantity =  QTY_PRECIPITABLE_WATER
      if(obstype == 153                    ) obs_kind = GPS_PRECIPITABLE_WATER
   endif

   if(obtype(1:3) == '  t') then
     obs_quantity = QTY_TEMPERATURE
     if(obstype == 120 .or. obstype == 132) obs_kind = RADIOSONDE_TEMPERATURE
     if(obstype == 130 .or. obstype == 131) obs_kind = AIRCRAFT_TEMPERATURE
     if(obstype == 133 .or. obstype == 134) obs_kind = ACARS_TEMPERATURE ! CSS HRRRE calls this AMDAR_TEMPERATURE
     if(obstype == 135                    ) obs_kind = AMDAR_TEMPERATURE
     if(obstype == 161 .or. obstype == 163) obs_kind = ATOV_TEMPERATURE
     if(obstype == 171 .or. obstype == 173) obs_kind = ATOV_TEMPERATURE
     if(obstype == 180 .or. obstype == 182) obs_kind = MARINE_SFC_TEMPERATURE
     if(obstype == 181 .or. obstype == 183 .or. obstype == 187 ) obs_kind = LAND_SFC_TEMPERATURE
     if(obstype == 188 )                    obs_kind = MESONET_TEMPERATURE
   endif

   if(obtype(1:3) == '  q') then
     obs_quantity = QTY_SPECIFIC_HUMIDITY
     if(obstype == 120 .or. obstype == 132) obs_kind = RADIOSONDE_RELATIVE_HUMIDITY  ! RADIOSONDE_SPECIFIC_HUMIDITY
     if(obstype == 130 .or. obstype == 131) obs_kind =   AIRCRAFT_RELATIVE_HUMIDITY  ! AIRCRAFT_SPECIFIC_HUMIDITY
     if(obstype == 133 .or. obstype == 134) obs_kind =      ACARS_RELATIVE_HUMIDITY  ! ACARS_SPECIFIC_HUMIDITY  ! Really should be AMDAR per HRRRE changes, but that requires some obs_def/preprocess changes.
     if(obstype == 180 .or. obstype == 182) obs_kind = MARINE_SFC_RELATIVE_HUMIDITY  ! MARINE_SFC_SPECIFIC_HUMIDITY
     if(obstype == 181 .or. obstype == 183  .or. obstype == 187 ) obs_kind = LAND_SFC_RELATIVE_HUMIDITY    ! LAND_SFC_SPECIFIC_HUMIDITY
     if(obstype == 188 )                    obs_kind =    MESONET_RELATIVE_HUMIDITY
   endif

   if(obtype(1:3) == ' ps') then
     obs_quantity = QTY_SURFACE_PRESSURE
     ! what to use: PRESSURE or ALTIMETER ?
     if(obstype == 120                    ) obs_kind = RADIOSONDE_SURFACE_PRESSURE ! RADIOSONDE_SURFACE_ALTIMETER
     if(obstype == 180 .or. obstype == 182) obs_kind =         MARINE_SFC_PRESSURE ! MARINE_SFC_ALTIMETER
     if(obstype == 181 .or. obstype == 187) obs_kind =           LAND_SFC_PRESSURE ! LAND_SFC_ALTIMETER
     if(obstype == 188 )                    obs_kind =    MESONET_SURFACE_PRESSURE
   endif

   if(obtype(1:3) == '  u') then
     obs_quantity = QTY_U_WIND_COMPONENT
     if(obstype == 220 .or. obstype == 232)  obs_kind = RADIOSONDE_U_WIND_COMPONENT
     if(obstype == 221                    )  obs_kind = RADIOSONDE_U_WIND_COMPONENT
     if(obstype == 223 .or. obstype == 227)  obs_kind = PROFILER_U_WIND_COMPONENT
     if(obstype == 224                    )  obs_kind = VADWND_U_WIND_COMPONENT ! CSS added new type ! PILOT_U_WIND_COMPONENT !VADWND actually
     if(obstype == 229                    )  obs_kind = PILOT_U_WIND_COMPONENT
     if(obstype == 230 .or. obstype == 231)  obs_kind = AIRCRAFT_U_WIND_COMPONENT
     if(obstype == 233 .or. obstype == 234)  obs_kind = ACARS_U_WIND_COMPONENT
     if(obstype == 235                    )  obs_kind = AMDAR_U_WIND_COMPONENT
     if(obstype == 280 .or. obstype == 282)  obs_kind = MARINE_SFC_U_WIND_COMPONENT
     if(obstype == 281 .or. obstype == 284  .or. obstype == 287 )  obs_kind = LAND_SFC_U_WIND_COMPONENT
     if(obstype == 288                    )  obs_kind = MESONET_U_WIND_COMPONENT
     if(obstype >= 240 .and. obstype <= 259) obs_kind = SAT_U_WIND_COMPONENT
     ! 285 is QSCAT
     if(obstype == 285                    )  obs_kind = QKSWND_U_WIND_COMPONENT
     ! 289 is WINDSAT, 290 is ASCAT, 291 is OSCAT
     if(obstype >= 289 .and. obstype <= 291) obs_kind = QKSWND_U_WIND_COMPONENT
   endif

   if(obtype(1:3) == '  v') then
     obs_quantity = QTY_V_WIND_COMPONENT
     if(obstype == 220 .or. obstype == 232)  obs_kind = RADIOSONDE_V_WIND_COMPONENT
     if(obstype == 221                    )  obs_kind = RADIOSONDE_V_WIND_COMPONENT
     if(obstype == 223 .or. obstype == 227)  obs_kind = PROFILER_V_WIND_COMPONENT
     if(obstype == 224                    )  obs_kind = VADWND_V_WIND_COMPONENT ! CSS added new type !PILOT_V_WIND_COMPONENT !VADWND actually
     if(obstype == 229                    )  obs_kind = PILOT_V_WIND_COMPONENT
     if(obstype == 230 .or. obstype == 231)  obs_kind = AIRCRAFT_V_WIND_COMPONENT
     if(obstype == 233 .or. obstype == 234)  obs_kind = ACARS_V_WIND_COMPONENT
     if(obstype == 235                    )  obs_kind = AMDAR_V_WIND_COMPONENT
     if(obstype == 280 .or. obstype == 282)  obs_kind = MARINE_SFC_V_WIND_COMPONENT
     if(obstype == 281 .or. obstype == 284 .or. obstype == 287 )  obs_kind = LAND_SFC_V_WIND_COMPONENT
     if(obstype == 288                    )  obs_kind = MESONET_V_WIND_COMPONENT ! CSS added
     if(obstype >= 240 .and. obstype <= 259) obs_kind = SAT_V_WIND_COMPONENT
     ! 285 is QSCAT
     if(obstype == 285                    )  obs_kind = QKSWND_V_WIND_COMPONENT
     ! 289 is WINDSAT, 290 is ASCAT, 291 is OSCAT
     if(obstype >= 289 .and. obstype <= 291) obs_kind = QKSWND_V_WIND_COMPONENT
   endif

   if (debug .and. obs_kind < 0) then
      ! the "real" fix if the record type is not found might actually be to
      ! accept all record types within valid ranges, and depend on the first
      ! preprocessing steps (in the prepbufr converter) to remove obs record
      ! types which are not desired.  for now, avoid giving them the wrong type
      ! and quietly loop. (the looping happens just after this function returns)
      write(string1,*) 'unsupported obtype,report type combination: ', obtype, obstype
      call error_handler(E_MSG,'prepbufr_to_dart_obs_kind',string1,'dart_obs_seq_mod.f90')
   endif

   which_vert = VERTISPRESSURE
!  if ( obstype >= 180 .and. obstype <= 191 ) which_vert = VERTISSURFACE ! CSS maybe keep as pressure
!  if ( obstype >= 280 .and. obstype <= 291 ) which_vert = VERTISSURFACE
end subroutine prepbufr_to_dart_obs_kind

subroutine radiance_to_dart_obs_kind(obtype, channel, obs_kind, which_vert, obs_quantity)
   character(len=*),  intent(in)  :: obtype  !'ps','u','v','oz','t','amsua_n15', etc.
   integer,           intent(in)  :: channel ! radiance channel
   integer,           intent(out) :: obs_kind, which_vert
   integer,           intent(out) :: obs_quantity

   character(len=256)             :: str_channel, this_string

   ! Make a string containing satellite,sensor,and channel to match obs_def_radiance_mod.f90 entries
   write(str_channel,fmt='(i8)') channel
   this_string = trim(adjustl(obtype))//'_ch'//trim(adjustl(str_channel))

   ! Make uppercase to match the case in obs_def_radiance_mod.f90
   call to_upper(this_string)

   ! Bad things will happen in DART preprocess "obs_def" files if there are hyphens
   !  in names. Replace with underscores
   call replace_hyphen(this_string)
   
   ! Be careful about upper/lower case.  Make sure this matches the obs_def
   obs_kind     = get_index_for_type_of_obs(this_string)    ! from obs_kind_mod
   obs_quantity = QTY_TEMPERATURE ! QTY_RADIANCE
   which_vert   = VERTISPRESSURE

end subroutine radiance_to_dart_obs_kind

subroutine replace_hyphen(string)
   character(*), intent(inout) :: string
   integer :: j
   do j = 1,len_trim(string)
      if ( string(j:j) == '-') string(j:j) = '_'
   enddo
end subroutine replace_hyphen

function write_this_ob_type_external_FO(ob_type)
   integer, intent(in)           :: ob_type
   logical                       :: write_this_ob_type_external_FO

   integer                       :: i
   logical                       :: is_all
   character(len=obstypelength)  :: ob_type_string

   ! Initialize to false
   write_this_ob_type_external_FO = .false.
   is_all = .false.

   ! Get the ob_type string corresponding to the ob_type integer
   ob_type_string = get_name_for_type_of_obs(ob_type)   ! from obs_kind_mod

   ! Determine whether we want to write out this ob type
   ! ntypes_to_compute_FO_max,write_FO_for_these_obs_types from params
   do i=1,ntypes_to_compute_FO_max
     if( trim(adjustl(write_FO_for_these_obs_types(i))) == trim(adjustl(ob_type_string)) .or. &
         trim(adjustl(write_FO_for_these_obs_types(i))) == 'all' ) then   ! If any element is 'all', it gets output
        write_this_ob_type_external_FO = .true.
        is_all = trim(adjustl(write_FO_for_these_obs_types(i))) == 'all'
     endif
   enddo

   if ( is_all ) then
      do i=1,ntypes_to_compute_FO_max
        if( trim(adjustl(exclude_these_obs_types(i))) == trim(adjustl(ob_type_string)) ) then
           write_this_ob_type_external_FO = .false.
        endif
      enddo
   endif

end function write_this_ob_type_external_FO


subroutine compute_mean_std(esize,evalue,emean,estd)

    integer,      intent(in)    :: esize
    real(r_kind), intent(in)    :: evalue(esize)
    real(r_kind), intent(inout) :: emean, estd

    emean = sum(evalue)/esize
    estd  = sqrt( sum((evalue-emean)**2)/max(1,(esize-1)) )

end subroutine compute_mean_std

end module dart_obs_seq_mod
