! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program Fluxnetfull_to_obs

!------------------------------------------------------------------------
!
!   Fluxnetfull_to_obs - a program that converts Ameriflux/Fluxnet FULLSET eddy
!                        covariance tower data of carbon, water and energy into
!                        DART obs_seq formatted files
!        ## Fixme  Add more description here once completed

use         types_mod, only : r8, MISSING_R8

use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              register_module, error_handler, E_MSG, E_ERR, &
                              open_file, close_file, do_nml_file, do_nml_term, &
                              check_namelist_read, find_namelist_in_file, &
                              nmlfileunit, logfileunit

use  time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                              set_date, set_time, get_time, print_time, &
                              print_date, operator(-), operator(+), operator(>), &
                              operator(<), operator(==), operator(<=), operator(>=)

use      location_mod, only : VERTISHEIGHT

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data

use      obs_kind_mod, only : TOWER_SENSIBLE_HEAT_FLUX, &
                              TOWER_NETC_ECO_EXCHANGE,  &
                              TOWER_LATENT_HEAT_FLUX,   &
                              TOWER_ER_FLUX,            &
                              TOWER_GPP_FLUX

implicit none

character(len=*), parameter :: source   = 'Fluxnetfull_to_obs.f90'
character(len=512), parameter :: string1,string2,string3

!-----------------------------------------------------------------------
! Namelist with default values
!-----------------------------------------------------------------------

character(len=256) :: text_input_file = 'textdata.input'
character(len=256) :: obs_out_file    = 'obs_seq.out'
real(r8)           :: timezoneoffset  = -1.0_r8
real(r8)           :: latitude        = -1.0_r8
real(r8)           :: longitude       = -1.0_r8
real(r8)           :: elevation       = -1.0_r8
real(r8)           :: flux_height     = -1.0_r8
! A maxgooqc=3 allows for good,medium and poor quality gap-filled data
real(r8)           :: maxgoodqc       = 3.0_r8
! Always 'true' except for latent,sensible heat and NEE for hourly time periods
! Fixme This option must be worked into the code later on
logical            :: gap_filled      = .true.
logical            :: verbose         = .false.

namelist /Fluxnetfull_to_obs_nml/ text_input_file, obs_out_file, &
             timezoneoffset, latitude, longitude, elevation, &
             flux_height, maxgoodqc, gap_filled, verbose

!-----------------------------------------------------------------------
! globally-scoped variables
!-----------------------------------------------------------------------

character(len=300)      :: input_line, bigline
character(len=512)      :: string1, string2, string3
integer                 :: iline, nlines, nwords
logical                 :: first_obs
integer                 :: oday, osec, rcio, iunit
integer                 :: num_copies, num_qc, max_obs
real(r8)                :: oerr, qc
real(r8)                :: sig_gppdt, sig_gppnt, sig_recodt, sig_recont
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: prev_time, offset
real(r8), parameter     :: umol_to_gC = (1.0_r8/1000000.0_r8) * 12.0_r8

! Fixme: These are for high resolution format HH or HR
! Fixme: If aggregrated (DD,WW,MM) need to edit

type towerdata
  type(time_type)   :: time_obs
  character(len=20) :: startstring    = 'TIMESTAMP_START'
  character(len=20) :: endstring      = 'TIMESTAMP_END'
  character(len=20) :: neestring      = 'NEE_VUT_REF'
  character(len=20) :: neeUNCstring   = 'NEE_VUT_REF_JOINTUNC'
  character(len=20) :: neeQCstring    = 'NEE_VUT_REF_QC'
  character(len=20) :: lestring       = 'LE_F_MDS'
  character(len=20) :: leUNCstring    = 'LE_RANDUNC'
  character(len=20) :: leQCstring     = 'LE_F_MDS_QC'
  character(len=20) :: hstring        = 'H_F_MDS'
  character(len=20) :: hUNCstring     = 'H_RANDUNC'
  character(len=20) :: hQCstring      = 'H_F_MDS_QC'
  character(len=20) :: gppDTstring    = 'GPP_DT_VUT_REF'
  character(len=20) :: gppNTstring    = 'GPP_NT_VUT_REF'
  character(len=20) :: recoDTstring    = 'RECO_DT_VUT_REF'
  character(len=20) :: recoNTstring    = 'RECO_NT_VUT_REF'
  character(len=20) :: gppUNCNT16string   = 'GPP_NT_VUT_16'
  character(len=20) :: gppUNCNT84string   = 'GPP_NT_VUT_84'
  character(len=20) :: gppUNCDT16string   = 'GPP_DT_VUT_16'
  character(len=20) :: gppUNCDT84string   = 'GPP_DT_VUT_84'
  character(len=20) :: recoUNCNT16string   = 'RECO_NT_VUT_16'
  character(len=20) :: recoUNCNT84string   = 'RECO_NT_VUT_84'
  character(len=20) :: recoUNCDT16string   = 'RECO_DT_VUT_16'
  character(len=20) :: recoUNCDT84string   = 'RECO_DT_VUT_84'
  integer  :: startindex
  integer  :: endindex
  integer  :: neeindex
  integer  :: neeUNCindex
  integer  :: neeQCindex
  integer  :: leindex
  integer  :: leUNCindex
  integer  :: leQCindex
  integer  :: hindex
  integer  :: hUNCindex
  integer  :: hQCindex
  integer  :: gppDTindex
  integer  :: gppNTindex
  integer  :: recoDTindex
  integer  :: recoNTindex
  integer  :: gppUNCDT16index
  integer  :: gppUNCDT84index
  integer  :: gppUNCNT16index
  integer  :: gppUNCNT84index
  integer  :: recoUNCDT16index
  integer  :: recoUNCDT84index
  integer  :: recoUNCNT16index
  integer  :: recoUNCNT84index  

  character(len=12) :: start_time    
  character(len=12) :: end_time      
  real(r8) :: nee
  real(r8) :: neeUNC
  integer  :: neeQC
  real(r8) :: le
  real(r8) :: leUNC
  integer  :: leQC
  real(r8) :: h
  real(r8) :: hUNC
  integer  :: hQC
  real(r8) :: gppNT
  real(r8) :: gppDT
  real(r8) :: gpp
  real(r8) :: gppNTQC
  real(r8) :: gppDTQC
  real(r8) :: gppUNCNT84
  real(r8) :: gppUNCNT16
  real(r8) :: gppUNCDT84
  real(r8) :: gppUNCDT16
  real(r8) :: recoNT
  real(r8) :: recoDT
  real(r8) :: reco
  real(r8) :: recoNTQC
  real(r8) :: recoDTQC
  real(r8) :: recoUNCNT84
  real(r8) :: recoUNCNT16
  real(r8) :: recoUNCDT84
  real(r8) :: recoUNCDT16
end type towerdata

type(towerdata) :: tower

!-----------------------------------------------------------------------
! start of executable code
!-----------------------------------------------------------------------

call initialize_utilities('Fluxnetfull_to_obs')

! Read the namelist entry
call find_namelist_in_file("input.nml", "Fluxnetfull_to_obs_nml", iunit)
read(iunit, nml = Fluxnetfull_to_obs_nml, iostat = rcio)
call check_namelist_read(iunit, rcio, "Fluxnetfull_to_obs_nml")

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=Fluxnetfull_to_obs_nml)
if (do_nml_term()) write(     *     , nml=Fluxnetfull_to_obs_nml)

! time setup
call set_calendar_type(GREGORIAN)
offset    = set_time(nint(abs(timezoneoffset)*3600.0_r8),0)
prev_time = set_time(0, 0)

write(string1, *) 'tower located at lat, lon, elev  =', latitude, longitude, elevation
write(string2, *) 'flux observations taken at       =', flux_height,'m'
if (verbose) call error_handler(E_MSG,'Fluxnetfull_to_obs',string1,text2=string2)

! check the lat/lon values to see if they are ok
if (longitude < 0.0_r8) longitude = longitude + 360.0_r8

if (( latitude > 90.0_r8 .or. latitude  <  -90.0_r8 ) .or. &
    (longitude <  0.0_r8 .or. longitude >  360.0_r8 )) then

   write (string2,*)'latitude  should be [-90, 90] but is ',latitude
   write (string3,*)'longitude should be [  0,360] but is ',longitude

   string1 ='tower location error in input.nml&Fluxnetfull_to_obs_nml'
   call error_handler(E_ERR,'Fluxnetfull_to_obs', source, string1, &
                      text2=string2,text3=string3)

endif

! Specify the maximum number of observations in the input file,
! but only the actual number created will be written out.
! Each line has 5 flux observations available.
! Each observation in this series will have a single
! observation value and a quality control flag.  
! Initialize two empty observations - one to track location
! in observation sequence - the other is for the new observation.

iunit = open_file(text_input_file, 'formatted', 'read')
if (verbose) call error_handler(E_MSG,'Fluxnetfull_to_obs','opened input file '//trim(text_input_file))

nlines     = count_file_lines(iunit)
max_obs    = 5*nlines
num_copies = 1
num_qc     = 1
first_obs  = .true.

call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)

! the first one needs to contain the string 'observation' and the
! second needs the string 'QC'.
call set_copy_meta_data(obs_seq, 1, 'observation')
call set_qc_meta_data(  obs_seq, 1, 'Fluxnet QC')

! Subroutine decode_header reads obs file header and identifies the columns 
! where flux variables of interest are located
rewind(iunit)
call decode_header(iunit, nwords)

obsloop: do iline = 2,nlines

   ! read in entire text line into a buffer
   read(iunit,'(A)',iostat=rcio) bigline
   if (rcio < 0) exit obsloop
   if (rcio > 0) then
      write (string1,'(''Cannot read (error '',i3,'') line '',i8,'' in '',A)') &
                    rcio, iline, trim(text_input_file)
      call error_handler(E_ERR,'main', string1, source)
   endif

   input_line = adjustl(bigline)

   ! Parse the header line into the tower structure (including the observation time)
   call stringparse(input_line, nwords, iline)

   if (iline <= 2) then
      write(*,*)''
      call print_date(tower%time_obs, ' first observation date (local time) is')
      call print_time(tower%time_obs, ' first observation time (local time) is')
      write(*,*)'first observation raw values: (column,string,value) timezone not applied'
      write(*,*)tower%hindex          , tower%hstring          , tower%h
      write(*,*)tower%hUNCindex       , tower%hUNCstring       , tower%hUNC
      write(*,*)tower%hQCindex        , tower%hQCstring        , tower%hQC
      write(*,*)tower%leindex         , tower%lestring         , tower%le
      write(*,*)tower%leUNCindex      , tower%leUNCstring      , tower%leUNC
      write(*,*)tower%leQCindex       , tower%leQCstring       , tower%leQC
      write(*,*)tower%neeindex        , tower%neestring        , tower%nee
      write(*,*)tower%neeUNCindex     , tower%neeUNCstring     , tower%UNCnee
      write(*,*)tower%neeQCindex      , tower%neeQCstring      , tower%neeQC
      write(*,*)tower%gppDTindex      , tower%gppDTstring      , tower%gppDT
      write(*,*)tower%gppDTUNC16index , tower%gppDTUNC16string , tower%gppDTUNC16
      write(*,*)tower%gppDTUNC84index , tower%gppDTUNC84string , tower%gppDTUNC84
      write(*,*)tower%gppDTQC 
      write(*,*)tower%gppNTindex      , tower%gppNTstring      , tower%gppNT  
      write(*,*)tower%gppNTUNC16index , tower%gppNTUNC16string , tower%gppNTUNC16
      write(*,*)tower%gppNTUNC84index , tower%gppNTUNC84string , tower%gppNTUNC84
      write(*,*)tower%gppNTQC 
      write(*,*)tower%recoDTindex     , tower%recoDTstring     , tower%recoDT  
      write(*,*)tower%recoDTUNC16index, tower%recoDTUNC16string, tower%recoDTUNC16
      write(*,*)tower%recoDTUNC84index, tower%recoDTUNC84string, tower%recoDTUNC84
      write(*,*)tower%recoDTQC 
      write(*,*)tower%recoNTindex     , tower%recoNTstring     , tower%recoNT      
      write(*,*)tower%recoNTUNC16index, tower%recoNTUNC16string, tower%recoNTUNC16
      write(*,*)tower%recoNTUNC84index, tower%recoNTUNC84string, tower%recoNTUNC84
      write(*,*)tower%recoNTQC
      write(*,*)''

      write(logfileunit,*)''
      call print_date(tower%time_obs, ' first observation date (local time) is',logfileunit)
      call print_time(tower%time_obs, ' first observation time (local time) is',logfileunit)
      write(logfileunit,*)'first observation raw values: (column,string,value) timezone not applied'
      write(logfileunit,*)tower%hindex          , tower%hstring          , tower%h
      write(logfileunit,*)tower%hUNCindex       , tower%hUNCstring       , tower%hUNC
      write(logfileunit,*)tower%hQCindex        , tower%hQCstring        , tower%hQC
      write(logfileunit,*)tower%leindex         , tower%lestring         , tower%le
      write(logfileunit,*)tower%leUNCindex      , tower%leUNCstring      , tower%leUNC
      write(logfileunit,*)tower%leQCindex       , tower%leQCstring       , tower%leQC
      write(logfileunit,*)tower%neeindex        , tower%neestring        , tower%nee
      write(logfileunit,*)tower%neeQCindex      , tower%neeQCstring      , tower%neeQC
      write(logfileunit,*)tower%gppDTUNC16index , tower%gppDTUNC16string , tower%gppDTUNC16
      write(logfileunit,*)tower%gppDTUNC84index , tower%gppDTUNC84string , tower%gppDTUNC84
      write(logfileunit,*)tower%gppDTQC 
      write(logfileunit,*)tower%gppNTindex      , tower%gppNTstring      , tower%gppNT
      write(logfileunit,*)tower%gppNTUNC16index , tower%gppNTUNC16string , tower%gppNTUNC16
      write(logfileunit,*)tower%gppNTUNC84index , tower%gppNTUNC84string , tower%gppNTUNC84
      write(logifleunit,*)tower%gppNTQC
      write(logfileunit,*)tower%recoDTindex     , tower%recoDTstring     , tower%recoDT
      write(logfileunit,*)tower%recoDTUNC16index, tower%recoDTUNC16string, tower%recoDTUNC16
      write(logfileunit,*)tower%recoDTUNC84index, tower%recoDTUNC84string, tower%recoDTUNC84
      write(logfileunit,*)tower%recoDTQC
      write(logfileunit,*)tower%recoNTindex     , tower%recoNTstring     , tower%recoNT
      write(logfileunit,*)tower%recoNTUNC16index, tower%recoNTUNC16string, tower%recoNTUNC16
      write(logfileunit,*)tower%recoNTUNC84index, tower%recoNTUNC84string, tower%recoNTUNC84
      write(logfileunit,*)tower%recoNTQC
      write(logfileunit,*)''
   end if

   call get_time(tower%time_obs, osec, oday)

   if (verbose) then
      write(string1,*)'obs time is (seconds,days) ',osec, oday,' obs date is '
      call print_date(tower%time_obs, trim(string1))
      call print_date(tower%time_obs, trim(string1),logfileunit)
   endif

   ! Create and add observation and uncertainty (1 SD) to obs_seq file
   ! Assign the observation the appropriate obs type
   if (tower%hQC <= maxgoodqc) then   ! Sensible Heat Flux [W m-2]
      oerr = tower%hUNC
      qc   = real(tower%hQC,r8)
      call create_3d_obs(latitude, longitude, flux_height, VERTISHEIGHT, tower%h, &
                         TOWER_SENSIBLE_HEAT_FLUX, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, tower%time_obs, prev_obs, prev_time, first_obs)
   endif

   if (tower%leQC <= maxgoodqc) then   ! Latent Heat Flux [W m-2]
      oerr = tower%leUNC
      qc   = real(tower%leQC,r8)
      call create_3d_obs(latitude, longitude, flux_height, VERTISHEIGHT, tower%le, &
                         TOWER_LATENT_HEAT_FLUX, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, tower%time_obs, prev_obs, prev_time, first_obs)
   endif

   if (tower%neeQC <= maxgoodqc) then       ! Net Ecosystem Exchange  [umol m-2 s-1]
      oerr      = tower%neeUNC * umol_to_gC
      tower%nee = -tower%nee * umol_to_gC   ! Matches units in CLM [gC m-2 s-1]
      qc        = real(tower%neeQC,r8)
      call create_3d_obs(latitude, longitude, flux_height, VERTISHEIGHT, tower%nee, &
                         TOWER_NETC_ECO_EXCHANGE, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, tower%time_obs, prev_obs, prev_time, first_obs)
   endif

   if (tower%gppDTQC .and. tower%gppNTQC <= maxgoodqc) then    ! Gross Primary Production  [umol m-2 s-1]
      sig_gppdt = (((tower%gppDTUNC84-tower%gppDTUNC16 / 2))^2)^0.5  ! Ustar unc contribution
      sig_gppnt = (((tower%gppNTUNC84-tower%gppNTUNC16 / 2))^2)^0.5  
      
      oerr      = (0.25 * (sig_gppdt)^2 + 0.25 * (sig_gppnt)^2)^0.5  ! Combine Ustar and partitioning unc
      oerr      = oerr * umol_to_gC
      ! Take average of night and day partitioning methods
      tower%gpp = ((tower%gppDT + tower%gppNT) / 2)  * umol_to_gC    ! Matches units in CLM [gC m-2 s-1]
      qc        = maxval(real(tower%gppDTQC,r8),real(tower%gppNTQC,r8))
      call create_3d_obs(latitude, longitude, flux_height, VERTISHEIGHT, tower%gpp, &
                         TOWER_GPP_FLUX, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, tower%time_obs, prev_obs, prev_time, first_obs)
   endif

   if (tower%recoDTQC .and. tower%recoNTQC <= maxgoodqc) then   ! Gross Primary Production  [umol m-2 s-1]
      sig_recodt = (((tower%recoDTUNC84-tower%recoDTUNC16 / 2))^2)^0.5  ! Ustar unc contribution
      sig_recont = (((tower%recoNTUNC84-tower%recoNTUNC16 / 2))^2)^0.5  

      oerr      = (0.25 * (sig_recodt)^2 + 0.25 * (sig_recont)^2)^0.5  ! Combine Ustar and partitioning unc
      oerr      = oerr * umol_to_gC
      ! Take average of night and day partitioning methods
      tower%reco = ((tower%recoDT + tower%recoNT) / 2)  * umol_to_gC    ! Matches units in CLM [gC m-2 s-1]
      qc        = maxval(real(tower%recoDTQC,r8),real(tower%recoNTQC,r8))
      call create_3d_obs(latitude, longitude, flux_height, VERTISHEIGHT, tower%reco, &
                         TOWER_ER_FLUX, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, tower%time_obs, prev_obs, prev_time, first_obs)
   endif

end do obsloop

! If obs added to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   write(string1,*)'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   if (verbose) call error_handler(E_MSG,'Fluxnetfull_to_obs',string1)
   call write_obs_seq(obs_seq, obs_out_file)
endif

! end of main program
call finalize_utilities()

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   create_3d_obs - subroutine that is used to create an observation
!                   type from observation data.
!
!       NOTE: assumes the code is using the threed_sphere locations module,
!             that the observation has a single data value and a single
!             qc value, and that this obs type has no additional required
!             data (e.g. gps and radar obs need additional data per obs)
!
!    lat   - latitude of observation
!    lon   - longitude of observation
!    vval  - vertical coordinate
!    vkind - kind of vertical coordinate (pressure, level, etc)
!    obsv  - observation value
!    okind - observation kind
!    oerr  - observation error (in units of standard deviation)
!    day   - gregorian day
!    sec   - gregorian second
!    qc    - quality control value
!    obs   - observation type
!
!     created Oct. 2007 Ryan Torn, NCAR/MMM
!     adapted for more generic use 11 Mar 2010, nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine create_3d_obs(lat, lon, vval, vkind, obsv, okind, oerr, day, sec, qc, obs)
use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_type_of_obs, &
                             set_obs_def_error_variance, set_obs_def_location
use obs_sequence_mod, only : obs_type, set_obs_values, set_qc, set_obs_def
use time_manager_mod, only : time_type, set_time
use     location_mod, only : set_location

 integer,        intent(in)    :: okind, vkind, day, sec
 real(r8),       intent(in)    :: lat, lon, vval, obsv, oerr, qc
 type(obs_type), intent(inout) :: obs

real(r8)           :: obs_val(1), qc_val(1)
type(obs_def_type) :: obs_def

call set_obs_def_location(obs_def, set_location(lon, lat, vval, vkind))
call set_obs_def_type_of_obs(obs_def, okind)
call set_obs_def_time(obs_def, set_time(sec, day))
call set_obs_def_error_variance(obs_def, oerr * oerr)
call set_obs_def(obs, obs_def)

obs_val(1) = obsv
call set_obs_values(obs, obs_val)
qc_val(1)  = qc
call set_qc(obs, qc_val)

end subroutine create_3d_obs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   add_obs_to_seq -- adds an observation to a sequence.  inserts if first
!           obs, inserts with a prev obs to save searching if that's possible.
!
!     seq - observation sequence to add obs to
!     obs - observation, already filled in, ready to add
!     obs_time - time of this observation, in dart time_type format
!     prev_obs - the previous observation that was added to this sequence
!                (will be updated by this routine)
!     prev_time - the time of the previously added observation (will also
!                be updated by this routine)
!     first_obs - should be initialized to be .true., and then will be
!                updated by this routine to be .false. after the first obs
!                has been added to this sequence.
!
!     created Mar 8, 2010   nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine add_obs_to_seq(seq, obs, obs_time, prev_obs, prev_time, first_obs)

use obs_sequence_mod, only : obs_sequence_type, obs_type, insert_obs_in_seq
use time_manager_mod, only : time_type, operator(>=)

type(obs_sequence_type), intent(inout) :: seq
type(obs_type),          intent(inout) :: obs, prev_obs
type(time_type),         intent(in)    :: obs_time
type(time_type),         intent(inout) :: prev_time
logical,                 intent(inout) :: first_obs

! insert(seq,obs) always works (i.e. it inserts the obs in
! proper time format) but it can be slow with a long file.
! supplying a previous observation that is older (or the same
! time) as the new one speeds up the searching a lot.

if(first_obs) then    ! for the first observation, no prev_obs
   call insert_obs_in_seq(seq, obs)
   first_obs = .false.
else
   if(obs_time >= prev_time) then  ! same time or later than previous obs
      call insert_obs_in_seq(seq, obs, prev_obs)
   else                            ! earlier, search from start of seq
      call insert_obs_in_seq(seq, obs)
   endif
endif

! update for next time
prev_obs  = obs
prev_time = obs_time

end subroutine add_obs_to_seq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   count_file_lines --
!           count the lines in a text file.
!           rewinds the unit after counting.
!
!     iunit - handle to the already-open text file
!
!     created May 2, 2012   Tim Hoar, NCAR/IMAGe
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function count_file_lines(iunit)

integer, intent(in) :: iunit
integer :: count_file_lines

integer :: i
character(len=128) :: oneline

integer, parameter :: tenmillion = 10000000
rewind(iunit)

count_file_lines = 0
countloop : do i = 1,tenmillion

   read(iunit,'(A)',iostat=rcio) oneline

   if (rcio < 0) exit countloop ! end of file
   if (rcio > 0) then
      write (string1,'('' read around line '',i8)')i
      call error_handler(E_ERR,'count_file_lines', string1, source)
   endif
   count_file_lines = count_file_lines + 1

enddo countloop
rewind(iunit)

if (count_file_lines >= tenmillion) then
   write (string1,'('' suspiciously large number of lines '',i8)')count_file_lines
   call error_handler(E_MSG,'count_file_lines', string1, source)
endif

end function count_file_lines




subroutine decode_header(iunit,ncolumns)
! Reads the first line of the obs header and identifies which columns
! the flux variable of interest is located

integer, intent(in) :: iunit
integer, intent(out) :: ncolumns

integer, parameter :: maxwordlength = 30
integer :: i,charcount,columncount,wordlength
character(len=maxwordlength), dimension(:), allocatable :: columns
integer, dimension(10) :: qc = 0

! Read the line and strip off any leading whitespace.

read(iunit,'(A)',iostat=rcio) bigline
if (rcio /= 0) then
  write(string1,*)'Cannot parse header. Begins <',trim(bigline(1:40)),'>'
  call error_handler(E_ERR,'decode_header',string1, source)
endif

input_line = adjustl(bigline)

! Comma separated file, thus count commas to determine the number of columns

charcount = CountChar(input_line,',')
ncolumns  = charcount + 1
allocate(columns(ncolumns))

columncount  = 0  ! track the number of columns
wordlength   = 0  ! number of characters in the column descriptor
charcount    = 0  ! the position of the (last) comma
do i = 1,len_trim(input_line)
   if (input_line(i:i) == ',') then
      columncount = columncount + 1
      if (wordlength > maxwordlength) then
         write(string1,*)'unexpected long word ... starts <',&
                           input_line((i-wordlength):(i-1)),'>'
         call error_handler(E_ERR,'decode_header',string1, source)
      endif
      columns(columncount) = input_line((i-wordlength):(i-1)) 
      write(string1,*) 'word(',columncount,') is ',columns(columncount)
      if (verbose) call error_handler(E_MSG,'decode_header',string1)
      wordlength = 0
      charcount = i
   else
      wordlength = wordlength + 1
   endif
enddo

! There is one more column after the last comma

if ((columncount+1) /= ncolumns) then
    write(string1,*)'parsed wrong number of words ...'
    write(string2,*)'expected ',ncolumns,' got ',columncount+1
    call error_handler(E_ERR,'decode_header',string1,source, &
                       text2=trim(string2), text3=trim(input_line))
endif

columns(ncolumns) = input_line((charcount+1):len_trim(input_line))

write(string1,*)'word(',ncolumns,') is ',columns(ncolumns)
if (verbose) call error_handler(E_MSG,'decode_header',string1)

! Finally, identify column index based on string name
tower%startindex       = Match(columns, tower%startstring)
tower%endindex         = Match(columns, tower%endstring)
tower%neeindex         = Match(columns, tower%neestring)
tower%neeUNCindex      = Match(columns, tower%neeUNCstring)
tower%neeQCindex       = Match(columns, tower%neeQCstring)
tower%neeUNCQCindex    = Match(columns, tower%neeUNCstring)
tower%leindex          = Match(columns, tower%lestring)
tower%leQCindex        = Match(columns, tower%leQCstring)
tower%leUNCindex       = Match(columns, tower%leUNCstring)
tower%hindex           = Match(columns, tower%hstring)
tower%hQCindex         = Match(columns, tower%hQCstring)
tower%hUNCindex        = Match(columns, tower%hUNCstring)
tower%gppDTindex       = Match(columns, tower%gppDTstring)
tower%gppNTindex       = Match(columns, tower%gppNTstring)
tower%recoDTindex      = Match(columns, tower%recoDTstring)
tower%recoNTindex      = Match(columns, tower%recoNTstring)
tower%gppUNCDT16index  = Match(columns, tower%gppUNCDT16string)
tower%gppUNCDT84index  = Match(columns, tower%gppUNCDT84string)
tower%gppUNCNT16index  = Match(columns, tower%gppUNCNT16string)
tower%gppUNCNT84index  = Match(columns, tower%gppUNCNT84string)
tower%recoUNCDT16index = Match(columns, tower%recoUNCDT16string)
tower%recoUNCDT84index = Match(columns, tower%recoUNCDT84string)
tower%recoUNCNT16index = Match(columns, tower%recoUNCNT16string)
tower%recoUNCNT84index = Match(columns, tower%recoUNCNT84string)

! Confirm indices were found successfully
qc( 1) = CheckIndex( tower%startindex      , tower%startstring)
qc( 2) = CheckIndex( tower%endindex        , tower%endstring)
qc( 3) = CheckIndex( tower%neeindex        , tower%neestring)
qc( 4) = CheckIndex( tower%neeUNCindex     , tower%neeUNCstring)
qc( 5) = CheckIndex( tower%neeQCindex      , tower%neeQCstring)
qc( 6) = CheckIndex( tower%neeUNCQCindex   , tower%neeUNCstring)
qc( 7) = CheckIndex( tower%leindex         , tower%lestring)
qc( 8) = CheckIndex( tower%leQCindex       , tower%leQCstring)
qc( 9) = CheckIndex( tower%leUNCindex      , tower%leUNCstring)
qc(10) = CheckIndex( tower%hindex          , tower%hstring)
qc(11) = CheckIndex( tower%hQCindex        , tower%hQCstring)
qc(12) = CheckIndex( tower%hUNCindex       , tower%hUNCstring)
qc(13) = CheckIndex( tower%gppDTindex      , tower%gppDTstring)
qc(14) = CheckIndex( tower%gppNTindex      , tower%gppNTstring)
qc(15) = CheckIndex( tower%recoDTindex     , tower%recoDTstring)
qc(16) = CheckIndex( tower%recoNTindex     , tower%recoNTstring)
qc(17) = CheckIndex( tower%gppUNCDT16index , tower%gppUNCDT16string)
qc(18) = CheckIndex( tower%gppUNCDT84index , tower%gppUNCDT84string)
qc(19) = CheckIndex( tower%gppUNCNT16index , tower%gppUNCNT16string)
qc(20) = CheckIndex( tower%gppUNCNT84index , tower%gppUNCNT84string)
qc(21) = CheckIndex( tower%recoUNCDT16index, tower%recoUNCDT16string)
qc(22) = CheckIndex( tower%recoUNCDT84index, tower%recoUNCDT84string)
qc(23) = CheckIndex( tower%recoUNCNT16index, tower%recoUNCNT16string)
qc(24) = CheckIndex( tower%recoUNCNT84index, tower%recoUNCNT84string)

if (any(qc /= 0) ) then
  write(string1,*)'Did not find all the required column indices.'
  call error_handler(E_ERR,'decode_header',string1, source)
endif

if (verbose) then
110 format('index for ',A20,' is ',i3)
   write(*,110)tower%startstring      ,tower%startindex
   write(*,110)tower%endstring        , tower%endindex
   write(*,110)tower%neestring        , tower%neeindex
   write(*,110)tower%neeUNCstring     , tower%neeUNCindex
   write(*,110)tower%neeQCstring      , tower%neeQCindex
   write(*,110)tower%neeUNCstring     , tower%neeUNCindex
   write(*,110)tower%lestring         , tower%leindex
   write(*,110)tower%leQCstring       , tower%leQCindex
   write(*,110)tower%leUNCstring      , tower%leUNCindex
   write(*,110)tower%hstring          , tower%hindex
   write(*,110)tower%hQCstring        , tower%hQCindex
   write(*,110)tower%hUNCstring       , tower%hUNCindex
   write(*,110)tower%gppDTstring      , tower%gppDTindex
   write(*,110)tower%gppNTstring      , tower%gppNTindex
   write(*,110)tower%recoDTstring     , tower%recoDTindex
   write(*,110)tower%recoNTstring     , tower%recoNTindex
   write(*,110)tower%gppUNCDT16string , tower%gppUNCDT16index
   write(*,110)tower%gppUNCDT84string , tower%gppUNCDT84index
   write(*,110)tower%gppUNCNT16string , tower%gppUNCNT16index
   write(*,110)tower%gppUNCNT84string , tower%gppUNCNT84index
   write(*,110)tower%recoUNCDT16string, tower%recoUNCDT16index
   write(*,110)tower%recoUNCDT84string, tower%recoUNCDT84index
   write(*,110)tower%recoUNCNT16string, tower%recoUNCNT16index
   write(*,110)tower%recoUNCNT84string, tower%recoUNCNT84index
endif

deallocate(columns)

end subroutine decode_header



function CountChar(str1,solo)
! Count the number of instances of the single character in a character string.
! useful when parsing a comma-separated list, for example.
! Count the commas and add 1 to get the number of items in the list.

integer                      :: CountChar
character(len=*), intent(in) :: str1
character,        intent(in) :: solo

integer :: i

CountChar = 0
do i = 1,len_trim(str1)
   if (str1(i:i) == solo) CountChar = CountChar + 1
enddo

end function CountChar



function Match(sentence,word)
! Determine the first occurrence of the 'word' in a sentence.
! In this context, a sentence is a character array, the dimension
! of the array is the number of words in the sentence.
! This is a case-sensitive match. Trailing blanks are removed.

integer :: Match
character(len=*), dimension(:), intent(in) :: sentence
character(len=*),               intent(in) :: word

integer :: i

Match = 0
WordLoop : do i = 1,size(sentence)
   if (trim(sentence(i)) == trim(word)) then
      Match = i
      return
   endif
enddo WordLoop

end function Match



function CheckIndex( myindex, context )
! Routine to issue a warning if the index was not found.
! Returns an error code ... 0 means the index WAS found
! a negative number means the index was NOT found - an error condition.
! ALL indices checked before fatally ending.

integer                       :: CheckIndex
integer,          intent(in)  :: myindex
character(len=*), intent(in)  :: context

if (myindex == 0) then
   write(string1,*)'Did not find column header matching ',trim(context)
   call error_handler(E_MSG,'decode_header:CheckIndex',string1, source)
   CheckIndex = -1 ! not a good thing
else
   CheckIndex = 0  ! Good to go
endif

end function CheckIndex



subroutine stringparse(str1, nwords, linenum)
! just declare everything as reals and chunk it

character(len=*), intent(in) :: str1
integer         , intent(in) :: nwords
integer         , intent(in) :: linenum

real(r8), allocatable, dimension(:) :: values
integer :: iyear0, imonth0, iday0, ihour0, imin0
integer :: iyear1, imonth1, iday1, ihour1, imin1
type(time_type) :: date_start, date_end

allocate(values(nwords))

values = MISSING_R8

read(str1,*,iostat=rcio) values
if (rcio /= 0) then
   write(string1,*)'Cannot parse line',linenum,'. Begins <',trim(str1(1:40)),'>'
   call error_handler(E_ERR,'stringparse',string1, source)
endif

! Assign flux observations, uncertainty and QC to tower structure
!
! Fixme: These are for high resolution format HH or HR only
! Fixme: If aggregrated (DD,WW,MM) need to edit, expand this

! start_time,end_time     format   YYYYMMDDHHMM
! nee           units    [umolCO2 m-2 s-1], VUT_REF
! neeUNC        units    [umolCO2 m-2 s-1], joint
! neeQC         no dim   0=measured,gap_filled: 1=good,2=medium,3=poor
! le            units    [W m-2]
! leUNC         units    [W m-2], random
! leQC          no dim   0=measured,gap_filled: 1=good gf,2=medium,3=poor
! h             units    [W m-2]
! hUNC          units    [W m-2], random
! hQC           no dim   0=measured,gap_filled: 1=good gf,2=medium,3=poor
! gppDT         units    [umolCO2 m-2 s-1] VUT_REF daytime partition
! gppNT         units    [umolCO2 m-2 s-1] VUT_REF nighttime partition
! gppUNCDT[xx]  units    [umolCO2 m-2 s-1] 16,84 percentile
! gppUNCNT[xx]  units    [umolCO2 m-2 s-1] 16,84 percentile
! recoDT        units    [umolCO2 m-2 s-1] VUT_REF daytime partition
! recoNT        units    [umolCO2 m-2 s-1] VUT_REF nighttime partition
! recoUNCDT[xx] units    [umolCO2 m-2 s-1] 16,84 percentile
! recoUNCNT[xx] units    [umolCO2 m-2 s-1] 16,84 percentile

! Convert to 'CLM-friendly' units AFTER we determine observation error variance.
! That happens in the main routine.

! Fixme: Double check these defs and units
! (CLM) NEE,GPP,ER  units     [gC m-2 s-1]
! (CLM) LE,SH       units     [W m-2]

tower%start_time  =      values(tower%startindex      )
tower%end_end     =      values(tower%endindex        )
tower%nee         =      values(tower%neeindex        )
tower%nee         =      values(tower%neeUNCindex     )
tower%neeQC       = nint(values(tower%neeQCindex     ))
tower%le          =      values(tower%leindex         )
tower%leUNC       =      values(tower%leUNCindex      )
tower%leQC        = nint(values(tower%leQCindex      ))
tower%h           =      values(tower%hindex          )
tower%hUNC        =      values(tower%hUNCindex       )
tower%hQC         = nint(values(tower%hQCindex       ))
tower%gppDT       =      values(tower%gppDTindex      )
tower%gppNT       =      values(tower%gppNTindex      )
tower%recoDT      =      values(tower%recoDTindex     )
tower%recoNT      =      values(tower%recoNTindex     )
tower%gppUNCDT16  =      values(tower%gppUNCDT16index )
tower%gppUNCDT84  =      values(tower%gppUNCDT84index )
tower%gppUNCNT16  =      values(tower%gppUNCNT16index )
tower%gppUNCNT84  =      values(tower%gppUNCNT84index )
tower%recoUNCDT16 =      values(tower%recoUNCDT16index)
tower%recoUNCDT84 =      values(tower%recoUNCDT84index)
tower%recoUNCNT16 =      values(tower%recoUNCNT16index)
tower%recoUNCNT84 =      values(tower%recoUNCNT84index)
deallocate(values)

! The observation time must be assigned through the flux timestamp string

write(tower%start_time,'(i4, 4i2)') iyear0,imonth0,iday0,ihour0,imin0
write(tower%end_time,  '(i4, 4i2)') iyear1,imonth1,iday1,ihour1,imin1

date_start= set_date(iyear0,imonth0,iday0,ihour0,imin0,0)
date_end=   set_date(iyear1,imonth1,iday1,ihour1,imin0,0)

! Assign average of flux time window to obs_seq (DART time)
tower%time_obs = (date_start+date_end) / 2

! Covert from Fluxnet provided LTC to UTC, UTC is standard for DART and CLM 
! For example, EST = UTC-5 
if (timezoneoffset < 0.0_r8) then
   tower%time_obs = tower%time_obs + offset
else
   tower%time_obs = tower%time_obs - offset
endif

! If missing value (-9999) manually assign poor QC value
! such that value is excluded in obs_seq
if (tower%neeQC < 0) tower%neeQC = maxgoodqc + 1000 
if (tower%leQC  < 0) tower%leQC  = maxgoodqc + 1000
if (tower%hQC   < 0) tower%hQC   = maxgoodqc + 1000

tower%gppNTQC = 1
tower%gppDTQC = 1
tower%recoNTQC = 1
tower%recoDTQC = 1

if (tower%gppNT < 0) tower%gppNTQC = maxgoodqc + 1000
if (tower%gppDT < 0) tower%gppDTQC = maxgoodqc + 1000
if (tower%recoNT < 0) tower%recoNTQC = maxgoodqc + 1000
if (tower%recoDT < 0) tower%recoDTQC = maxgoodqc + 1000
 
end subroutine stringparse



end program Fluxnetfull_to_obs


