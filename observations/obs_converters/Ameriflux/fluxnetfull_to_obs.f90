! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program fluxnetfull_to_obs

!------------------------------------------------------------------------
!
!   fluxnetfull_to_obs - a program that converts Ameriflux/Fluxnet FULLSET eddy
!                        covariance tower data of NEE, GPP, RE, latent heat and sensible
!                        heat fluxes into DART obs_seq formatted files. Works on native time 
!                        resolution (HH, HR) that ED tower data or collected or at coarser
!                        aggregated time resolution files (DD,WW,MM) based on ONEFLUX
!                        methodology (Pastorello et al., 2020)
    

use         types_mod, only : r8, MISSING_R8

use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              error_handler, E_MSG, E_ERR, &
                              open_file, close_file, do_nml_file, do_nml_term, &
                              check_namelist_read, find_namelist_in_file, &
                              nmlfileunit, logfileunit

use  time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                              set_date, set_time, get_time, print_time, &
                              print_date, operator(-), operator(+), operator(>), &
                              operator(<), operator(==), operator(<=), operator(>=), &
                              operator(/)

use      location_mod, only : VERTISHEIGHT

use  obs_sequence_mod, only : obs_sequence_type, obs_type, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data

use obs_utilities_mod, only : add_obs_to_seq, create_3d_obs

use      obs_kind_mod, only : TOWER_SENSIBLE_HEAT_FLUX, &
                              TOWER_NETC_ECO_EXCHANGE,  &
                              TOWER_LATENT_HEAT_FLUX,   &
                              TOWER_ER_FLUX,            &
                              TOWER_GPP_FLUX

implicit none

character(len=*), parameter :: source   = 'fluxnetfull_to_obs.f90'

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
! A maxgooqc=3 allows for good=1, medium=2, and poor=3 quality gap-filled data
real(r8)           :: maxgoodqc       = 3.0_r8
! Always true except for latent,sensible heat and NEE for hourly time periods
! Can only be false of HH or HR time resolution
logical            :: gap_filled      = .true.
! Option for energy balance correction for latent and sensible heat
! Recommend to keep false as these values are typically missing
logical            :: energy_balance  = .false.
character(len=2)   :: time_resolution = 'HH'
logical            :: verbose         = .false.

namelist /fluxnetfull_to_obs_nml/ text_input_file, obs_out_file, &
             timezoneoffset, latitude, longitude, elevation, &
             flux_height, maxgoodqc, gap_filled, energy_balance, &
             time_resolution, verbose

!-----------------------------------------------------------------------
! globally-scoped variables
!-----------------------------------------------------------------------

character(len=5000)     :: input_line, bigline
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
logical                 :: res 

! Initialize with default tower strings, modify later

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
  character(len=20) :: gppNTUNC16string   = 'GPP_NT_VUT_16'
  character(len=20) :: gppNTUNC84string   = 'GPP_NT_VUT_84'
  character(len=20) :: gppDTUNC16string   = 'GPP_DT_VUT_16'
  character(len=20) :: gppDTUNC84string   = 'GPP_DT_VUT_84'
  character(len=20) :: recoNTUNC16string   = 'RECO_NT_VUT_16'
  character(len=20) :: recoNTUNC84string   = 'RECO_NT_VUT_84'
  character(len=20) :: recoDTUNC16string   = 'RECO_DT_VUT_16'
  character(len=20) :: recoDTUNC84string   = 'RECO_DT_VUT_84'
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
  integer  :: gppDTUNC16index
  integer  :: gppDTUNC84index
  integer  :: gppNTUNC16index
  integer  :: gppNTUNC84index
  integer  :: recoDTUNC16index
  integer  :: recoDTUNC84index
  integer  :: recoNTUNC16index
  integer  :: recoNTUNC84index  
  
  character(len=12) :: start_time    
  character(len=12) :: end_time
  real(r8) :: nee
  real(r8) :: neeUNC
  integer  :: neeQC
  real(r8) :: neeQCfrac
  real(r8) :: le
  real(r8) :: leUNC
  integer  :: leQC
  real(r8) :: h
  real(r8) :: hUNC
  integer  :: hQC
  real(r8) :: gppNT
  real(r8) :: gppDT
  real(r8) :: gpp
  integer  :: gppNTQC
  integer  :: gppDTQC
  real(r8) :: gppNTUNC84
  real(r8) :: gppNTUNC16
  real(r8) :: gppDTUNC84
  real(r8) :: gppDTUNC16
  real(r8) :: recoNT
  real(r8) :: recoDT
  real(r8) :: reco
  integer  :: recoNTQC
  integer  :: recoDTQC
  real(r8) :: recoNTUNC84
  real(r8) :: recoNTUNC16
  real(r8) :: recoDTUNC84
  real(r8) :: recoDTUNC16
end type towerdata

type(towerdata) :: tower


!-----------------------------------------------------------------------
! start of executable code
!-----------------------------------------------------------------------

call initialize_utilities('fluxnetfull_to_obs')

! Read the namelist entry
call find_namelist_in_file("input.nml", "fluxnetfull_to_obs_nml", iunit)
read(iunit, nml = fluxnetfull_to_obs_nml, iostat = rcio)
call check_namelist_read(iunit, rcio, "fluxnetfull_to_obs_nml")

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=fluxnetfull_to_obs_nml)
if (do_nml_term()) write(     *     , nml=fluxnetfull_to_obs_nml)

! time setup
call set_calendar_type(GREGORIAN)
offset    = set_time(nint(abs(timezoneoffset)*3600.0_r8),0)
prev_time = set_time(0, 0)

write(string1, *) 'tower located at lat, lon, elev  =', latitude, longitude, elevation
write(string2, *) 'flux observations taken at       =', flux_height,'m'
if (verbose) call error_handler(E_MSG,'fluxnetfull_to_obs',string1,text2=string2)

! check the lat/lon values to see if they are ok
if (longitude < 0.0_r8) longitude = longitude + 360.0_r8

if (( latitude > 90.0_r8 .or. latitude  <  -90.0_r8 ) .or. &
    (longitude <  0.0_r8 .or. longitude >  360.0_r8 )) then

   write (string2,*)'latitude  should be [-90, 90] but is ',latitude
   write (string3,*)'longitude should be [  0,360] but is ',longitude

   string1 ='tower location error in input.nml &fluxnetfull_to_obs_nml'
   call error_handler(E_ERR, source, string1, &
                      text2=string2,text3=string3)
endif

! Provide gap-filling warning, and provide error if gap-filling
! turned on with DD or coarser time resolution

if (gap_filled .eqv. .false.) then

   write(string1,*) 'WARNING!: gap filling is turned off. This is only recommended for'
   write(string2,*) 'NEE,LE,SH  observations 0=measured,gap_filled: 1=good gf,2=medium,3=poor. This approach' 
   write(string3,*) 'assigns QC values > maxgoodQC thus gap-filled data will not be written to obs_seq'
   call error_handler(E_MSG, source, string1, &
                      text2=string2,text3=string3)

   if (time_resolution == 'DD' .or. time_resolution == 'MM' .or. &
       time_resolution == 'WW') then
      write(string1,*) 'ERROR!: gap filling can only be turned off for native time'
      write(string2,*) 'resolution data (HH, HR)'
      write(string3,*) 'Coarser time resolution (DD, WW, MM)  must have gap_filled = .true.'
      call error_handler(E_ERR, source, string1, &
                      text2=string2,text3=string3)
   endif

endif


! Modify values to towerdata strings based on user input.nml
! 1) time_resolution: DD(daily),MM(monthly),YY(yearly) uses 'TIMESTAMP' header

if (time_resolution == 'DD' .or. time_resolution == 'MM' .or. &
    time_resolution == 'YY') then

  tower%startstring    = 'TIMESTAMP'
  tower%endstring      = 'TIMESTAMP'
  res                  = .false.  

  write(string1, *) 'Time resolution is set to =', time_resolution
  write(string2, *) 'Using TIMESTAMP to set DART observation time'
  if (verbose) call error_handler(E_MSG,'fluxnetfull_to_obs',string1,text2=string2)

elseif (time_resolution == 'HH' .or. time_resolution == 'HR' .or. &
        time_resolution == 'WW') then

  res = .true.
  write(string1, *) 'Time resolution is set to =', time_resolution
  write(string2, *) 'Using TIMESTAMP_START and TIMESTAMP_END to set: DART observation time'
  if (verbose) call error_handler(E_MSG,'fluxnetfull_to_obs',string1,text2=string2)

else
  write(string1,*) 'time_resolution set incorrectly within input.nml'
  write(string2,*) 'time_resolution must be HR,HH,DD,WW,MM, or YY'
  call error_handler(E_ERR, source, string1,text2=string2)

endif

! 2) energy_balance: .true. changes sensible and latent heat strings
!    Note: There are no qc values for energy_balance correction
!    The qc values are manually set later in code 
if (energy_balance .eqv. .true.) then

   tower%lestring       = 'LE_CORR'
   tower%leUNCstring    = 'LE_CORR_JOINTUNC'
   tower%hstring        = 'H_CORR'
   tower%hUNCstring     = 'H_CORR_JOINTUNC'
   write(string1,*) 'WARNING! Energy balance correction data turned on for LE and H'
   write(string2,*) 'Check to make sure LE_CORR and H_CORR data is not missing'
   call error_handler(E_MSG, source, string1,text2=string2)
else

   write(string1,*) 'Standard LE (LE_F_MDS) and H (H_F_MDS) data being used'
   call error_handler(E_MSG, source, string1)
endif

! Specify the maximum number of observations in the input file,
! but only the actual number created will be written out.
! Each line has 5 flux observations available.
! Each observation in this series will have a single
! observation value and a quality control flag.  
! Initialize two empty observations - one to track location
! in observation sequence - the other is for the new observation.

iunit = open_file(text_input_file, 'formatted', 'read')
if (verbose) call error_handler(E_MSG,'fluxnetfull_to_obs','opened input file '//trim(text_input_file))

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
      write(*,*)tower%neeUNCindex     , tower%neeUNCstring     , tower%neeUNC
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
      write(logfileunit,*)tower%neeQCindex      , tower%neeQCstring      , tower%neeQCfrac
      write(logfileunit,*)tower%gppDTUNC16index , tower%gppDTUNC16string , tower%gppDTUNC16
      write(logfileunit,*)tower%gppDTUNC84index , tower%gppDTUNC84string , tower%gppDTUNC84
      write(logfileunit,*)tower%gppDTQC 
      write(logfileunit,*)tower%gppNTindex      , tower%gppNTstring      , tower%gppNT
      write(logfileunit,*)tower%gppNTUNC16index , tower%gppNTUNC16string , tower%gppNTUNC16
      write(logfileunit,*)tower%gppNTUNC84index , tower%gppNTUNC84string , tower%gppNTUNC84
      write(logfileunit,*)tower%gppNTQC
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
      write(string1,*)'obs time (UTC) is (seconds,days) ',osec, oday,' obs date (UTC) is '
      call print_date(tower%time_obs, trim(string1))
      call print_date(tower%time_obs, trim(string1),logfileunit)
   

      write(*,*)''
      write(string1, *) 'Display tower%start_time and tower%end_time (LTC)  =', tower%start_time,' ', tower%end_time
      call error_handler(E_MSG,'fluxnetfull_to_obs',string1)
      write(*,*)''
      write(string1, *) 'Display tower%nee tower%neeQC  tower%neeQCfrac =', tower%nee, tower%neeQC, tower%neeQCfrac
      call error_handler(E_MSG,'fluxnetfull_to_obs',string1)
      write(*,*)''
      write(string1, *) 'Display tower%le tower%leQC  =', tower%le, tower%leQC
      call error_handler(E_MSG,'fluxnetfull_to_obs',string1)
      write(*,*)''
      write(string1, *) 'Display tower%h tower%hQC  =', tower%h, tower%hQC
      call error_handler(E_MSG,'fluxnetfull_to_obs',string1)
      write(*,*)''
      write(string1, *) 'Display tower%gppDT tower%gppDTQC  =', tower%gppDT, tower%gppDTQC
      call error_handler(E_MSG,'fluxnetfull_to_obs',string1)
      write(*,*)''
      write(string1, *) 'Display tower%gppNT tower%gppNTQC  =', tower%gppNT, tower%gppNTQC
      call error_handler(E_MSG,'fluxnetfull_to_obs',string1)
      write(*,*)''
      write(string1, *) 'Display tower%recoDT tower%recoDTQC  =', tower%recoDT, tower%recoDTQC
      call error_handler(E_MSG,'fluxnetfull_to_obs',string1)
      write(*,*)''
      write(string1, *) 'Display tower%recoNT tower%recoNTQC  =', tower%recoNT, tower%recoNTQC
      call error_handler(E_MSG,'fluxnetfull_to_obs',string1)

   endif

   ! Create and add observation and uncertainty (1 SD) to obs_seq file
   ! Assign the observation the appropriate obs type
   if (tower%hQC <= maxgoodqc) then   ! Sensible Heat Flux [W m-2]
      oerr = tower%hUNC
      qc   = real(tower%hQC,r8)
      ! Check for missing uncertainty value, if needed assign % error
      ! based on empirical examination of data      
      if (oerr <=0) then
         select case(time_resolution)
           case ('HR', 'HH')
              oerr= tower%h*0.2
           case ('DD', 'WW')
              oerr= tower%h*0.1
           case ('MM')
              oerr= tower%h*0.05 
           case default
              write(string1, *) 'ERROR, time_resolution must be HH,HR,DD,WW,MM, value is:', time_resolution
              call error_handler(E_ERR,'fluxnetfull_to_obs',string1)
         end select      
      endif
      call create_3d_obs(latitude, longitude, flux_height, VERTISHEIGHT, tower%h, &
                         TOWER_SENSIBLE_HEAT_FLUX, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, tower%time_obs, prev_obs, prev_time, first_obs)

      if (verbose) then
         write(*,*)''
         write(string1, *) 'Display tower%h, tower%hUNC (1 SD)  =', tower%h,' ', tower%hUNC
         call error_handler(E_MSG,'fluxnetfull_to_obs',string1)
      endif

   endif

   if (tower%leQC <= maxgoodqc) then   ! Latent Heat Flux [W m-2]
      oerr = tower%leUNC
      qc   = real(tower%leQC,r8)
      if (oerr <=0) then
         select case( time_resolution )
           case ('HR', 'HH')
              oerr= tower%le*0.2
           case ('DD', 'WW')
              oerr= tower%le*0.1
           case ('MM')
              oerr= tower%le*0.05
           case default
              write(string1, *) 'ERROR, time_resolution must be HH,HR,DD,WW,MM, value is:', time_resolution 
              call error_handler(E_ERR,'fluxnetfull_to_obs',string1)
         end select      
      endif

      call create_3d_obs(latitude, longitude, flux_height, VERTISHEIGHT, tower%le, &
                         TOWER_LATENT_HEAT_FLUX, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, tower%time_obs, prev_obs, prev_time, first_obs)

      if (verbose) then
         write(*,*)''
         write(string1, *) 'Display tower%le, tower%leUNC  (1 SD) =', tower%le,' ', tower%leUNC
         call error_handler(E_MSG,'fluxnetfull_to_obs',string1)
      endif

   endif

   if (tower%neeQC <= maxgoodqc) then       ! Net Ecosystem Exchange  [umol m-2 s-1]
      oerr      = tower%neeUNC * umol_to_gC
      tower%nee = -tower%nee * umol_to_gC   ! Matches units in CLM [gC m-2 s-1]
      qc        = real(tower%neeQC,r8)
      if (oerr <=0) then
         select case( time_resolution )
           case ('HR', 'HH')
              oerr= abs(tower%nee)*0.2
           case ('DD', 'WW')
              oerr= abs(tower%nee)*0.1
           case ('MM')
              oerr= abs(tower%nee)*0.05
           case default
              write(string1, *) 'ERROR, time_resolution must be HH,HR,DD,WW,MM, value is:', time_resolution 
              call error_handler(E_ERR,'fluxnetfull_to_obs',string1)
         end select      
      endif
      call create_3d_obs(latitude, longitude, flux_height, VERTISHEIGHT, tower%nee, &
                         TOWER_NETC_ECO_EXCHANGE, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, tower%time_obs, prev_obs, prev_time, first_obs)

      if (verbose) then
         write(*,*)''
         write(string1, *) 'Display tower%nee, tower%neeUNC  (1 SD) =', tower%nee,' ', tower%neeUNC
         call error_handler(E_MSG,'fluxnetfull_to_obs',string1)
      endif

   endif



   if (tower%gppDTQC <=maxgoodqc .and. tower%gppNTQC <= maxgoodqc) then    ! Gross Primary Production  [umol m-2 s-1]
      sig_gppdt = (((tower%gppDTUNC84-tower%gppDTUNC16) / 2)**2)**0.5  ! Ustar unc contribution
      sig_gppnt = (((tower%gppNTUNC84-tower%gppNTUNC16) / 2)**2)**0.5  
      
      oerr      = (0.25 * (sig_gppdt)**2 + 0.25 * (sig_gppnt)**2)**0.5  ! Combine Ustar and partitioning unc
      oerr      = oerr * umol_to_gC
      ! Take average of night and day partitioning methods
      tower%gpp = ((tower%gppDT + tower%gppNT) / 2)  * umol_to_gC    ! Matches units in CLM [gC m-2 s-1]
      qc        = maxval((/real(tower%gppDTQC,r8),real(tower%gppNTQC,r8)/))
      if (oerr <=0) then
         select case( time_resolution )
           case ('HR', 'HH')
              oerr= tower%gpp*0.2
           case ('DD', 'WW')
              oerr= tower%gpp*0.1
           case ('MM')
              oerr= tower%gpp*0.05
           case default
              write(string1, *) 'ERROR, time_resolution must be HH,HR,DD,WW,MM, value is:', time_resolution 
              call error_handler(E_ERR,'fluxnetfull_to_obs',string1)
         end select      
      endif
      call create_3d_obs(latitude, longitude, flux_height, VERTISHEIGHT, tower%gpp, &
                         TOWER_GPP_FLUX, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, tower%time_obs, prev_obs, prev_time, first_obs)

      if (verbose) then
         write(*,*)''
         write(string1, *) 'Display tower%gpp, gpp uncertainty (1 SD)  =', tower%gpp,' ', oerr
         call error_handler(E_MSG,'fluxnetfull_to_obs',string1)
      endif

   endif

   if (tower%recoDTQC <= maxgoodqc .and. tower%recoNTQC <= maxgoodqc) then   ! Gross Primary Production  [umol m-2 s-1]
      sig_recodt = (((tower%recoDTUNC84-tower%recoDTUNC16) / 2)**2)**0.5  ! Ustar unc contribution
      sig_recont = (((tower%recoNTUNC84-tower%recoNTUNC16) / 2)**2)**0.5  

      oerr      = (0.25 * (sig_recodt)**2 + 0.25 * (sig_recont)**2)**0.5  ! Combine Ustar and partitioning unc
      oerr      = oerr * umol_to_gC
      ! Take average of night and day partitioning methods
      tower%reco = ((tower%recoDT + tower%recoNT) / 2)  * umol_to_gC    ! Matches units in CLM [gC m-2 s-1]
      qc        = maxval((/real(tower%recoDTQC,r8),real(tower%recoNTQC,r8)/))
      if (oerr <=0) then
         select case( time_resolution )
           case ('HR', 'HH')
              oerr= tower%reco*0.2
           case ('DD', 'WW')
              oerr= tower%reco*0.1
           case ('MM')
              oerr= tower%reco*0.05
           case default
              write(string1, *) 'ERROR, time_resolution must be HH,HR,DD,WW,MM, value is:', time_resolution 
              call error_handler(E_ERR,'fluxnetfull_to_obs',string1)
         end select      
      endif
      call create_3d_obs(latitude, longitude, flux_height, VERTISHEIGHT, tower%reco, &
                         TOWER_ER_FLUX, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, tower%time_obs, prev_obs, prev_time, first_obs)

      if (verbose) then
         write(*,*)''
         write(string1, *) 'Display tower%reco, reco uncertainty (1 SD)  =', tower%reco,' ', oerr
         call error_handler(E_MSG,'fluxnetfull_to_obs',string1)
      endif

   endif




end do obsloop

call close_file(iunit)

! If obs added to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   write(string1,*)'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   if (verbose) call error_handler(E_MSG,'fluxnetfull_to_obs',string1)
   call write_obs_seq(obs_seq, obs_out_file)
endif

! end of main program
call finalize_utilities()

contains


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
integer, dimension(23) :: qc = 0

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
tower%neeUNCindex      = Match(columns, tower%neeUNCstring)
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
tower%gppDTUNC16index  = Match(columns, tower%gppDTUNC16string)
tower%gppDTUNC84index  = Match(columns, tower%gppDTUNC84string)
tower%gppNTUNC16index  = Match(columns, tower%gppNTUNC16string)
tower%gppNTUNC84index  = Match(columns, tower%gppNTUNC84string)
tower%recoDTUNC16index = Match(columns, tower%recoDTUNC16string)
tower%recoDTUNC84index = Match(columns, tower%recoDTUNC84string)
tower%recoNTUNC16index = Match(columns, tower%recoNTUNC16string)
tower%recoNTUNC84index = Match(columns, tower%recoNTUNC84string)

! Confirm indices were found successfully
qc( 1) = CheckIndex( tower%startindex      , tower%startstring)
qc( 2) = CheckIndex( tower%endindex        , tower%endstring)
qc( 3) = CheckIndex( tower%neeindex        , tower%neestring)
qc( 4) = CheckIndex( tower%neeUNCindex     , tower%neeUNCstring)
qc( 5) = CheckIndex( tower%neeQCindex      , tower%neeQCstring)
qc( 6) = CheckIndex( tower%leindex         , tower%lestring)
qc( 7) = CheckIndex( tower%leQCindex       , tower%leQCstring)
qc( 8) = CheckIndex( tower%leUNCindex      , tower%leUNCstring)
qc( 9) = CheckIndex( tower%hindex          , tower%hstring)
qc(10) = CheckIndex( tower%hQCindex        , tower%hQCstring)
qc(11) = CheckIndex( tower%hUNCindex       , tower%hUNCstring)
qc(12) = CheckIndex( tower%gppDTindex      , tower%gppDTstring)
qc(13) = CheckIndex( tower%gppNTindex      , tower%gppNTstring)
qc(14) = CheckIndex( tower%recoDTindex     , tower%recoDTstring)
qc(15) = CheckIndex( tower%recoNTindex     , tower%recoNTstring)
qc(16) = CheckIndex( tower%gppDTUNC16index , tower%gppDTUNC16string)
qc(17) = CheckIndex( tower%gppDTUNC84index , tower%gppDTUNC84string)
qc(18) = CheckIndex( tower%gppNTUNC16index , tower%gppNTUNC16string)
qc(19) = CheckIndex( tower%gppNTUNC84index , tower%gppNTUNC84string)
qc(20) = CheckIndex( tower%recoDTUNC16index, tower%recoDTUNC16string)
qc(21) = CheckIndex( tower%recoDTUNC84index, tower%recoDTUNC84string)
qc(22) = CheckIndex( tower%recoNTUNC16index, tower%recoNTUNC16string)
qc(23) = CheckIndex( tower%recoNTUNC84index, tower%recoNTUNC84string)

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
   write(*,110)tower%gppDTUNC16string , tower%gppDTUNC16index
   write(*,110)tower%gppDTUNC84string , tower%gppDTUNC84index
   write(*,110)tower%gppNTUNC16string , tower%gppNTUNC16index
   write(*,110)tower%gppNTUNC84string , tower%gppNTUNC84index
   write(*,110)tower%recoDTUNC16string, tower%recoDTUNC16index
   write(*,110)tower%recoDTUNC84string, tower%recoDTUNC84index
   write(*,110)tower%recoNTUNC16string, tower%recoNTUNC16index
   write(*,110)tower%recoNTUNC84string, tower%recoNTUNC84index
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

character(len=*), intent(in) :: str1
integer         , intent(in) :: nwords
integer         , intent(in) :: linenum

real(r8), allocatable, dimension(:) :: values


integer :: yeara, yearb, montha, monthb, daya, dayb
integer :: houra, hourb, mina, minb

integer :: time_adjust
type(time_type) :: date_start, date_end

if (res .eqv. .true.) then
   time_adjust = 2 
else
   time_adjust = 1  
endif



allocate(values(nwords-time_adjust))
values = MISSING_R8

! First two/one element(s) of string read in as character (time stamp)
! remainder of line read in as reals.

if (res .eqv. .true.) then

   read(str1,*,iostat=rcio) tower%start_time,tower%end_time,values
   if (rcio /= 0) then
      write(string1,*)'Cannot parse line',linenum,'. Begins <',trim(str1(1:40)),'>'
      call error_handler(E_ERR,'stringparse',string1, source)
   endif
else
   read(str1,*,iostat=rcio) tower%start_time,values
   if (rcio /= 0) then
      write(string1,*)'Cannot parse line',linenum,'. Begins <',trim(str1(1:40)),'>'
      call error_handler(E_ERR,'stringparse',string1, source)
   endif
   tower%end_time=tower%start_time
endif


! Assign flux observations, uncertainty and QC to tower structure
! Description of FLUXNET FULLSET variables

! start_time,end_time    format   YYYYMMDDHHMM
! nee           units    [umolCO2 m-2 s-1], VUT_REF
! neeUNC        units    [umolCO2 m-2 s-1], joint
! neeQC (HH,HR) no dim   0=measured, gap_filled:1=good,2=medium,3=poor,-9999=missing
! neeQC (DD-MM) no dim   fraction between 0-1, % of good quality filled data
! le            units    [W m-2]
! leUNC         units    [W m-2], random
! leQC (HH,HR) no dim   0=measured, gap_filled:1=good,2=medium,3=poor,-9999=missing
! leQC (DD-MM) no dim   fraction between 0-1, % of good quality filled data
! h             units    [W m-2]
! hUNC          units    [W m-2], random
! hQC           no dim   0=measured, gap_filled:1=good,2=medium,3=poor,-9999=missing
! hQC (DD-MM) no dim   fraction between 0-1, % of good quality filled data
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

! CLM history file names and units
! (CLM) NEP,GPP,ER  units     [gC m-2 s-1]
! (CLM) EFLX_LH_TOT_R,FSH     units     [W m-2]

tower%nee         =           values(tower%neeindex        -time_adjust )
tower%neeUNC      =           values(tower%neeUNCindex     -time_adjust)
tower%neeQC       =      nint(values(tower%neeQCindex      -time_adjust))
tower%neeQCfrac   =           values(tower%neeQCindex      -time_adjust)
tower%le          =           values(tower%leindex         -time_adjust )
tower%leUNC       =           values(tower%leUNCindex      -time_adjust )
tower%leQC        =      nint(values(tower%leQCindex       -time_adjust))
tower%h           =           values(tower%hindex          -time_adjust )
tower%hUNC        =           values(tower%hUNCindex       -time_adjust )
tower%hQC         =      nint(values(tower%hQCindex        -time_adjust))
tower%gppDT       =           values(tower%gppDTindex      -time_adjust)
tower%gppNT       =           values(tower%gppNTindex      -time_adjust)
tower%recoDT      =           values(tower%recoDTindex     -time_adjust)
tower%recoNT      =           values(tower%recoNTindex     -time_adjust)
tower%gppDTUNC16  =           values(tower%gppDTUNC16index -time_adjust)
tower%gppDTUNC84  =           values(tower%gppDTUNC84index -time_adjust)
tower%gppNTUNC16  =           values(tower%gppNTUNC16index -time_adjust)
tower%gppNTUNC84  =           values(tower%gppNTUNC84index -time_adjust)
tower%recoDTUNC16 =           values(tower%recoDTUNC16index-time_adjust)
tower%recoDTUNC84 =           values(tower%recoDTUNC84index-time_adjust)
tower%recoNTUNC16 =           values(tower%recoNTUNC16index-time_adjust)
tower%recoNTUNC84 =           values(tower%recoNTUNC84index-time_adjust)
deallocate(values)



read(tower%start_time(1:12), fmt='(i4, 4i2)') yeara,montha,daya,houra,mina 
read(tower%end_time(1:12),   fmt='(i4, 4i2)') yearb,monthb,dayb,hourb,minb


! Certain time resolutions (MM) does not define days
if (time_resolution == 'MM') then
   daya = 1
   dayb = 1
endif



write(*,*)''
write(string1, *) 'Display tower%start_time,yeara,montha,daya,houra,mina (LTC)  =', tower%start_time, yeara, montha, daya, houra, mina
write(string2, *) 'Display tower%end_time,yearb,monthb,dayb,hourb,minb   (LTC)  =', tower%end_time, yearb, monthb, dayb, hourb, minb
if (verbose) call error_handler(E_MSG,'fluxnetfull_to_obs',string1,text2=string2)



date_start= set_date(yeara,montha,daya,houra,mina,0)
date_end=   set_date(yearb,monthb,dayb,hourb,minb,0)



! Assign average of flux time window to obs_seq (DART time)
tower%time_obs = (date_start+date_end) / 2

! Covert from Fluxnet provided LTC to UTC, UTC is standard for DART and CLM 
! For example, EST = UTC-5 
if (timezoneoffset < 0.0_r8) then
   tower%time_obs = tower%time_obs + offset
else
   tower%time_obs = tower%time_obs - offset
endif


! Reject NEE data where neeQC is missing
if (tower%neeQC < 0) tower%neeQC = maxgoodqc + 100

! If NEE (DD,WW,MM) convert from  % filled QC to integer QC  good/fair/poor
if (time_resolution == 'DD' .or. time_resolution == 'MM' .or. &
    time_resolution == 'WW') then
   if (tower%neeQCfrac >= 0.90) tower%neeQC = 1 ! >90% fill is good
   if (tower%neeQCfrac < 0.90 .and. tower%neeQC >= 0.60) tower%neeQC = 2 ! 60-90% fill is fair
   if (tower%neeQCfrac < 0.60 .and. tower%neeQC >= 0.0_r8) tower%neeQC = 3 ! <60% fill is poor
endif


! The QC values are typically missing for le and h (-9999)
! Thus  manually assign fair  QC values in these cases
! When energy balance is turned off
if (energy_balance .eqv. .false.) then
   if (tower%leQC  < 0.0_r8) tower%leQC  = 2
   if (tower%hQC   < 0.0_r8) tower%hQC   = 2
else  ! No QC values for energy balance corrected le and h.  Assign fair QC.
   tower%leQC  = 2
   tower%hQC  = 2
endif


! No qc values for gpp/reco, thus assign fair qc unless
! the gpp/reco values are negative values, then reject.
tower%gppNTQC = 2
tower%gppDTQC = 2
tower%recoNTQC = 2
tower%recoDTQC = 2

if (tower%gppNT < 0.0_r8) tower%gppNTQC = maxgoodqc + 100
if (tower%gppDT < 0.0_r8) tower%gppDTQC = maxgoodqc + 100
if (tower%recoNT < 0.0_r8) tower%recoNTQC = maxgoodqc + 100
if (tower%recoDT < 0.0_r8) tower%recoDTQC = maxgoodqc + 100

! Assign very bad qc to gap_filled data if user requests it
! such that maxgoodqc threshold does not add gap_filled data to obs_seq file
! If leQC and hQC are missing values already forced to 2, thus rejected
if (gap_filled .eqv. .false.) then
   if (tower%neeQC >0) tower%neeQC = maxgoodqc + 100
   if (tower%leQC >0)  tower%leQC =  maxgoodqc + 100
   if (tower%hQC >0)   tower%hQC =   maxgoodqc + 100
   ! GPP and RECO are modeled data and considered gap_filled
   tower%gppNTQC = maxgoodqc + 100
   tower%gppDTQC = maxgoodqc + 100
   tower%recoNTQC = maxgoodqc + 100
   tower%recoDTQC = maxgoodqc + 100
endif

end subroutine stringparse



end program fluxnetfull_to_obs


