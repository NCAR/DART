
!xlf90 -g -qrealsize=8 ens_obs_var.f90 -o ens_obs_var.exe

program ens_obs_to_dart

use         types_mod, only : r8, missing_r8, missing_data, DEG2RAD, earth_radius
use     utilities_mod, only : open_file, close_file, initialize_utilities, &
                              register_module, logfileunit, E_MSG, timestamp, &
                              error_handler, find_namelist_in_file, check_namelist_read
use  obs_sequence_mod, only : obs_type, obs_sequence_type, init_obs_sequence, insert_obs_in_seq, &
                              set_copy_meta_data, set_qc_meta_data, write_obs_seq, assignment(=), &
                              init_obs, static_init_obs_sequence, set_obs_def, set_obs_values, set_qc
use       obs_def_mod, only : set_obs_def_location, set_obs_def_error_variance, &
                              set_obs_def_kind, set_obs_def_time, set_obs_def_key, &
                              obs_def_type
use   obs_def_gps_mod, only : set_gpsro_ref
use      obs_kind_mod, only : SAT_U_WIND_COMPONENT, SAT_V_WIND_COMPONENT, &
                              QKSWND_U_WIND_COMPONENT, QKSWND_V_WIND_COMPONENT, &
                              RADIOSONDE_U_WIND_COMPONENT, RADIOSONDE_V_WIND_COMPONENT, &
                              RADIOSONDE_TEMPERATURE, RADIOSONDE_VAPOR_MIXING_RATIO, &
                              METAR_U_10_METER_WIND, METAR_V_10_METER_WIND, METAR_TEMPERATURE_2_METER, &
                              METAR_VAPOR_MIXING_RATIO_2_METER, &
                              BUOY_U_WIND_COMPONENT, BUOY_V_WIND_COMPONENT, &
                              BUOY_TEMPERATURE, BUOY_VAPOR_MIXING_RATIO, &
                              SHIP_U_WIND_COMPONENT, SHIP_V_WIND_COMPONENT, &
                              SHIP_TEMPERATURE, SHIP_VAPOR_MIXING_RATIO, &
                              SYNOP_U_WIND_COMPONENT, SYNOP_V_WIND_COMPONENT, &
                              SYNOP_TEMPERATURE, SYNOP_VAPOR_MIXING_RATIO, &
                              AIREP_U_WIND_COMPONENT, AIREP_V_WIND_COMPONENT, &
                              AIREP_TEMPERATURE, AIREP_VAPOR_MIXING_RATIO, &
                              AMDAR_U_WIND_COMPONENT, AMDAR_V_WIND_COMPONENT, &
                              AMDAR_TEMPERATURE, AMDAR_VAPOR_MIXING_RATIO, &
                              PILOT_U_WIND_COMPONENT, PILOT_V_WIND_COMPONENT, &
                              PILOT_TEMPERATURE, PILOT_VAPOR_MIXING_RATIO, &
                              PROFILER_U_WIND_COMPONENT, PROFILER_V_WIND_COMPONENT, &
                              BOGUS_U_WIND_COMPONENT, BOGUS_V_WIND_COMPONENT, &
                              BOGUS_TEMPERATURE, BOGUS_VAPOR_MIXING_RATIO, GPSRO_REFRACTIVITY
use      location_mod, only : location_type, set_location, VERTISSURFACE, VERTISPRESSURE, VERTISHEIGHT
use  time_manager_mod, only : time_type, set_date, set_calendar_type, GREGORIAN

implicit none

type(obs_sequence_type) :: seq
type(obs_type)          :: obs
type(obs_def_type)      :: obs_def
type(location_type)     :: location
type(time_type)         :: time
integer                 :: i, k, n, which_vert, q
integer                 :: num_obs, num_obs_last, num_copies, num_qc, max_num_obs=8000000
character(len = 129)    :: copy_meta_data, qc_meta_data

!    parameter (ne=16)
integer                 :: ens_size, ne, ie
integer                 :: ista, istamax=250000   !max. num. of stations
CHARACTER (LEN = 40)    :: name          ! Station name
CHARACTER (LEN = 12)    :: platform      ! Instrument platform
CHARACTER (LEN =  5)    :: id            ! 5 digit station identifer
CHARACTER (LEN = 19)    :: date          ! CCYY-MM-DD_HH:MM:SS date

CHARACTER (LEN = 80)    :: filename1
CHARACTER (LEN = 129)   :: obs_seq_out_filename = 'obs_seq.out'
CHARACTER (LEN = 120)   :: fmt_info, &
                           fmt_srfc, &
                           fmt_each, fmt_each_iv
character (len = 160)   :: header_string, eheader_string

real(r8)                :: missing_v, missing_qc, missing_e
integer                 :: iqc, iost, fm, ilev
real(r8)                :: xlon,xlat,elevation
real(r8)                :: data1,error1, data2,error2, &
                           data3,error3, data4,error4, data5,error5, &
                           data6,error6, data7,error7, data8,error8, &
                           data9,error9
real(r8)                :: uerr, verr, terr, qerr, gpserr
real(r8)                :: uqc(1), vqc(1), tqc(1), qqc(1), gpsqc(1)
integer                 :: iqc1, iqc2, iqc3, iqc4, iqc5, iqc6, iqc7, iqc8, iqc9

integer,dimension(:),allocatable  :: iunit
real(r8), dimension(:), allocatable   :: edata1,eerror1, edata2,eerror2, &
                                     edata3,eerror3, edata4,eerror4, edata5,eerror5, &
                                     edata6,eerror6, edata7,eerror7, edata8,eerror8, &
                                     edata9,eerror9, euinv,evinv,etinv,eqinv
integer,dimension(:), allocatable :: eiqc1, eiqc2, eiqc3, eiqc4, eiqc5, eiqc6, eiqc7, eiqc8, eiqc9
real(r8), dimension(:), allocatable   :: uens, vens, tens, qens, gpsens
logical                           :: valid_u, valid_v, valid_t, valid_q, valid_gps

type obs_info_type
     CHARACTER (LEN = 30)   :: name          ! Station name
     CHARACTER (LEN = 12)   :: platform      ! Instrument platform
     CHARACTER (LEN =  5)   :: id            ! 5 digit station identifer
     CHARACTER (LEN = 19)   :: date_char     ! CCYY-MM-DD_HH:MM:SS date
     INTEGER                :: levels        ! number of levels
     REAL(r8)               :: lat           ! Latitude in degree
     REAL(r8)               :: lon           ! Longitude in degree
     REAL(r8)               :: elv           ! Elevation in m
     REAL(r8)               :: p_station     ! station pressure in Pa
     REAL(r8)               :: p_lev         ! Pressure in Pa at the observation level 
     REAL(r8)               :: h_lev         ! Height in m at the observation level
end type obs_info_type

type(obs_info_type) :: obs_info

integer :: unitnml, io
namelist /ens_obs_to_dart_nml/ ens_size, obs_seq_out_filename

!!!!!!!!!!!!!!!!!!!!!!!!!!! Start !!!!!!!!!!!!!!!!!!!
!    print*,'ens_size = '
!    read(*,*) ens_size

     call find_namelist_in_file("input.nml", "ens_obs_to_dart_nml", unitnml)
     read(unitnml, nml = ens_obs_to_dart_nml, iostat = io)
     call check_namelist_read(unitnml, io, "ens_obs_to_dart_nml")

     ne=ens_size
     print*,'ens_size = ', ne
     num_copies = 5+ens_size*2
     print*,'num_copies = ', num_copies

     allocate(iunit(ne))
     allocate(edata1(ne), eerror1(ne), edata2(ne),eerror2(ne), &
              edata3(ne),eerror3(ne), edata4(ne),eerror4(ne), edata5(ne),eerror5(ne), &
              edata6(ne),eerror6(ne), edata7(ne),eerror7(ne), edata8(ne),eerror8(ne), &
              edata9(ne),eerror9(ne), euinv(ne),evinv(ne),etinv(ne),eqinv(ne))
     allocate(eiqc1(ne),eiqc2(ne),eiqc3(ne),eiqc4(ne),eiqc5(ne),eiqc6(ne),eiqc7(ne),eiqc8(ne),eiqc9(ne))

     allocate(uens(num_copies),vens(num_copies),tens(num_copies),qens(num_copies),gpsens(num_copies))

     fmt_info =    '(A12,1X,A19,1X,A40,1X,I6,3(F12.3,11X),6X,A5)'
     fmt_srfc =    '(F12.3,I4,F7.2,F12.3,I4,F7.3)'
     fmt_each =    '(3(F12.3,I4,F7.2),11X,3(F12.3,I4,F7.2),11X,1(F12.3,I4,F7.2)))'
     fmt_each_iv = '(3(F12.3,I4,F7.2),11X,3(F12.3,I4,F7.2),11X,1(F12.3,I4,F7.2),4f17.7)'

     missing_v  = -888888.000
     missing_qc = -88.0
     missing_e  = -888.0

! create obs_sequence
     call static_init_obs_sequence()

     call set_calendar_type(GREGORIAN)

! Initialize an obs_sequence structure, set meta data

     num_qc = 1
     call init_obs_sequence(seq, num_copies, num_qc, max_num_obs)

     copy_meta_data = 'observation'
     call set_copy_meta_data(seq, 1, copy_meta_data)

     qc_meta_data = 'WRFVAR Quality Control'
     call set_qc_meta_data(seq, 1, qc_meta_data)

     copy_meta_data = 'prior ensemble mean'
     call set_copy_meta_data(seq, 2, copy_meta_data)
     copy_meta_data = 'posterior ensemble mean'
     call set_copy_meta_data(seq, 3, copy_meta_data)
     copy_meta_data = 'prior ensemble spread'
     call set_copy_meta_data(seq, 4, copy_meta_data)
     copy_meta_data = 'posterior ensemble spread'
     call set_copy_meta_data(seq, 5, copy_meta_data)
     do ie=1,ens_size
        write(copy_meta_data,'(a21,1x,i6)') 'prior ensemble member', ie
        call set_copy_meta_data(seq,5+ie*2-1,copy_meta_data)
        write(copy_meta_data,'(a25,1x,i6)') 'posterior ensemble member', ie
        call set_copy_meta_data(seq,5+ie*2,copy_meta_data)
     enddo

     call init_obs(obs, num_copies, num_qc)

! open WRFVAR modified filtered obs file for one member, same format as all others
     open(12,file='filtered_obs',status='old')

     ! open ensemble of modified filtered obs files
     do n=1,ne
        iunit(n)=100+n
        write(filename1,'(a,i3.3)') 'filtered_obs.e',n
        open(iunit(n),file=filename1,status='old')
     enddo

     ! Loop over all files until end of header information is reached
     do
       read(unit=12,fmt='(a)',iostat=iost) header_string
       do n=1,ne
          read(unit=iunit(n),fmt='(a)',iostat=iost) eheader_string
       enddo
       if (iost /= 0) then
          print*, 'Problem reading obs header, exit'
          stop
       endif
       if (header_string(1:6) == '#-----' .or. header_string(1:21) == 'MODIFIED FILTERED OBS') exit
     enddo

 obs_loop: do ista=1,istamax  !CSS
        read(12,fmt=trim(fmt_info),end=999) platform,date,name,ilev,xlat,xlon,elevation,id
        read(platform(4:6),'(I3)') fm

        if(mod(ista,1000)==0)  print*, ista, platform, date, name
        read(12,fmt=trim(fmt_srfc),end=999) data1,iqc1,error1,data2,iqc2,error2

        ! assign obs_info structure
        obs_info%name=name
        obs_info%platform=platform
        obs_info%id=id
        obs_info%date_char=date
        obs_info%levels=ilev
        obs_info%elv=elevation
        obs_info%p_station=data1

        do n=1,ne
           read(iunit(n),fmt=trim(fmt_info),end=999)platform,date,name,ilev,xlat,xlon,elevation,id
           read(iunit(n),fmt=trim(fmt_srfc),end=999)edata1(n),eiqc1(n),eerror1(n),edata2(n),eiqc2(n),eerror2(n)
        enddo     !CSS...read 2 lines for each member

        do n=1,ne !CSS...loop over all the members to check for the error
           if(edata1(n).ne.data1 .or. eiqc1(n).ne.iqc1 .or. eerror1(n).ne.error1) then
             print*,'Observation mismatch!, ensemble member',n
             write(*,fmt=trim(fmt_info))platform,date,name,ilev,xlat,xlon,elevation,id
             write(*,fmt=trim(fmt_srfc))data1,iqc1,error1,data2,iqc2,error2
             write(*,fmt=trim(fmt_srfc))edata1(n),eiqc1(n),eerror1(n),edata2(n),eiqc2(n),eerror2(n)
! CSS         stop
	     print*,' '
	     print*,'skip this station...cycling through lines to get to start of next station' ! CSS
	     print*,' '
	     do k=1,ilev  ! CSS cycle through the appropriate number of lines (levels) before going to next station
                 read(12,fmt=trim(fmt_each),end=999) data3,iqc3,error3, data4,iqc4,error4, data5,iqc5,error5, &  
                                               data6,iqc6,error6, data7,iqc7,error7, data8,iqc8,error8, &
                                               data9,iqc9,error9

                 do q=1,ne   ! CSS just cycling through the lines of the ensemble members
                    read(iunit(q),fmt=trim(fmt_each_iv),end=999)   &
                       edata3(q),eiqc3(q),eerror3(q), edata4(q),eiqc4(q),eerror4(q), edata5(q),eiqc5(q),eerror5(q), &
                       edata6(q),eiqc6(q),eerror6(q), edata7(q),eiqc7(q),eerror7(q), edata8(q),eiqc8(q),eerror8(q), &
                       edata9(q),eiqc9(q),eerror9(q), euinv(q),evinv(q),etinv(q),eqinv(q)
                 enddo
             enddo

	     cycle obs_loop ! CSS break out of current cycle and go to next station
           endif
        enddo

        do k=1,ilev
           read(12,fmt=trim(fmt_each),end=999) data3,iqc3,error3, data4,iqc4,error4, data5,iqc5,error5, &
                                               data6,iqc6,error6, data7,iqc7,error7, data8,iqc8,error8, &
                                               data9,iqc9,error9
           obs_info%p_lev=data3
           obs_info%h_lev=data6

           ! 1st copy is obs value
           uens(1)=data4
           vens(1)=data5
           tens(1)=data7
           qens(1)=data9

           uerr=error4
           verr=error5
           terr=error7
           qerr=error9

           uqc(1)=real(iqc4)
           vqc(1)=real(iqc5)
           tqc(1)=real(iqc7)
           qqc(1)=real(iqc9)


           do n=1,ne
              read(iunit(n),fmt=trim(fmt_each_iv),end=999)   &
                  edata3(n),eiqc3(n),eerror3(n), edata4(n),eiqc4(n),eerror4(n), edata5(n),eiqc5(n),eerror5(n), &
                  edata6(n),eiqc6(n),eerror6(n), edata7(n),eiqc7(n),eerror7(n), edata8(n),eiqc8(n),eerror8(n), &
                  edata9(n),eiqc9(n),eerror9(n), euinv(n),evinv(n),etinv(n),eqinv(n)
           enddo

           if ( fm.eq.116 ) then ! GPS-RO refractivity
              gpsens(1) = data8
              gpserr = error8
              gpsqc(1) = real(iqc8)
              call get_ens_stat(edata8,eiqc8,eerror8,etinv,ne,gpsens,valid_gps) ! GPSRO, innovation is in "T" column
           else
              call get_ens_stat(edata4,eiqc4,eerror4,euinv,ne,uens,valid_u) ! U
              call get_ens_stat(edata5,eiqc5,eerror5,evinv,ne,vens,valid_v) ! V 
              call get_ens_stat(edata7,eiqc7,eerror7,etinv,ne,tens,valid_t) ! T 
              call get_ens_stat(edata9,eiqc9,eerror9,eqinv,ne,qens,valid_q) ! Q (g/kg)
           end if

           do n=1,num_copies
              if(qens(n) .gt. missing_r8) qens(n)=qens(n)/1e3  !convert back to kg/kg
           enddo
           if(qerr .gt. missing_r8) qerr=qerr/1e3

           if (valid_u .or. valid_v .or. valid_t .or. valid_q .or. valid_gps) then

              select case (fm)

              case (12) ;     !SYNOP
                   if (valid_u) call insert_ens_obs(SYNOP_U_WIND_COMPONENT, uens, uerr, uqc, &
                                                    obs_info, VERTISSURFACE, ne, obs, seq)
                   if (valid_v) call insert_ens_obs(SYNOP_V_WIND_COMPONENT, vens, verr, vqc, &
                                                    obs_info, VERTISSURFACE, ne, obs, seq)
                   if (valid_t) call insert_ens_obs(SYNOP_TEMPERATURE, tens, terr, tqc, &
                                                    obs_info, VERTISSURFACE, ne, obs, seq)
                   if (valid_q) call insert_ens_obs(SYNOP_VAPOR_MIXING_RATIO, qens, qerr, qqc, &
                                                    obs_info, VERTISSURFACE, ne, obs, seq)
              case (13) ;     !SHIP
                   if (valid_u) call insert_ens_obs(SHIP_U_WIND_COMPONENT, uens, uerr, uqc, &
                                                    obs_info, VERTISSURFACE, ne, obs, seq)
                   if (valid_v) call insert_ens_obs(SHIP_V_WIND_COMPONENT, vens, verr, vqc, &
                                                    obs_info, VERTISSURFACE, ne, obs, seq)
                   if (valid_t) call insert_ens_obs(SHIP_TEMPERATURE, tens, terr, tqc, &
                                                    obs_info, VERTISSURFACE, ne, obs, seq)
                   if (valid_q) call insert_ens_obs(SHIP_VAPOR_MIXING_RATIO, qens, qerr, qqc, &
                                                    obs_info, VERTISSURFACE, ne, obs, seq)
              case (15,16) ;  !METAR
                   if (valid_u) call insert_ens_obs(METAR_U_10_METER_WIND, uens, uerr, uqc, &
                                                    obs_info, VERTISSURFACE, ne, obs, seq)
                   if (valid_v) call insert_ens_obs(METAR_V_10_METER_WIND, vens, verr, vqc, &
                                                    obs_info, VERTISSURFACE, ne, obs, seq)
                   if (valid_t) call insert_ens_obs(METAR_TEMPERATURE_2_METER, tens, terr, tqc, &
                                                    obs_info, VERTISSURFACE, ne, obs, seq)
                   if (valid_q) call insert_ens_obs(METAR_VAPOR_MIXING_RATIO_2_METER, qens, qerr, qqc, &
                                                    obs_info, VERTISSURFACE, ne, obs, seq)
              case (18,19) ;  !BUOY
                   if (valid_u) call insert_ens_obs(BUOY_U_WIND_COMPONENT, uens, uerr, uqc, &
                                                    obs_info, VERTISSURFACE, ne, obs, seq)
                   if (valid_v) call insert_ens_obs(BUOY_V_WIND_COMPONENT, vens, verr, vqc, &
                                                    obs_info, VERTISSURFACE, ne, obs, seq)
                   if (valid_t) call insert_ens_obs(BUOY_TEMPERATURE, tens, terr, tqc, &
                                                    obs_info, VERTISSURFACE, ne, obs, seq)
                   if (valid_q) call insert_ens_obs(BUOY_VAPOR_MIXING_RATIO, qens, qerr, qqc, &
                                                    obs_info, VERTISSURFACE, ne, obs, seq)
              case (32:34) ;  !PILOT
                   if (valid_u) call insert_ens_obs(PILOT_U_WIND_COMPONENT, uens, uerr, uqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
                   if (valid_v) call insert_ens_obs(PILOT_V_WIND_COMPONENT, vens, verr, vqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
                   if (valid_t) call insert_ens_obs(PILOT_TEMPERATURE, tens, terr, tqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
                   if (valid_q) call insert_ens_obs(PILOT_VAPOR_MIXING_RATIO, qens, qerr, qqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
              case (35:38) ;  !TEMP (sounding)
                   if (valid_u) call insert_ens_obs(RADIOSONDE_U_WIND_COMPONENT, uens, uerr, uqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
                   if (valid_v) call insert_ens_obs(RADIOSONDE_V_WIND_COMPONENT, vens, verr, vqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
                   if (valid_t) call insert_ens_obs(RADIOSONDE_TEMPERATURE, tens, terr, tqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
                   if (valid_q) call insert_ens_obs(RADIOSONDE_VAPOR_MIXING_RATIO, qens, qerr, qqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
              case (42) ;     !AMDAR
                   if (valid_u) call insert_ens_obs(AMDAR_U_WIND_COMPONENT, uens, uerr, uqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
                   if (valid_v) call insert_ens_obs(AMDAR_V_WIND_COMPONENT, vens, verr, vqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
                   if (valid_t) call insert_ens_obs(AMDAR_TEMPERATURE, tens, terr, tqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
                   if (valid_q) call insert_ens_obs(AMDAR_VAPOR_MIXING_RATIO, qens, qerr, qqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
              case (88) ;     !SATOB
                   if (valid_u) call insert_ens_obs(SAT_U_WIND_COMPONENT, uens, uerr, uqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
                   if (valid_v) call insert_ens_obs(SAT_V_WIND_COMPONENT, vens, verr, vqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
              case (96:97) ;  !AIREP
                   if (valid_u) call insert_ens_obs(AIREP_U_WIND_COMPONENT, uens, uerr, uqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
                   if (valid_v) call insert_ens_obs(AIREP_V_WIND_COMPONENT, vens, verr, vqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
                   if (valid_t) call insert_ens_obs(AIREP_TEMPERATURE, tens, terr, tqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
                   if (valid_q) call insert_ens_obs(AIREP_VAPOR_MIXING_RATIO, qens, qerr, qqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
              case (116) ;    !GPS REFRACTIVITY
                   if (valid_gps) call insert_ens_obs(GPSRO_REFRACTIVITY, gpsens, gpserr, gpsqc, &
                                                    obs_info, VERTISHEIGHT, ne, obs, seq)
              case (132) ;    !PROFI
                   if (valid_u) call insert_ens_obs(PROFILER_U_WIND_COMPONENT, uens, uerr, uqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
                   if (valid_v) call insert_ens_obs(PROFILER_V_WIND_COMPONENT, vens, verr, vqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
              case (135) ;    !BOGUS
                   if (valid_u) call insert_ens_obs(BOGUS_U_WIND_COMPONENT, uens, uerr, uqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
                   if (valid_v) call insert_ens_obs(BOGUS_V_WIND_COMPONENT, vens, verr, vqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
                   if (valid_t) call insert_ens_obs(BOGUS_TEMPERATURE, tens, terr, tqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
                   if (valid_q) call insert_ens_obs(BOGUS_VAPOR_MIXING_RATIO, qens, qerr, qqc, &
                                                    obs_info, VERTISPRESSURE, ne, obs, seq)
              case (281) ;    !Quiks
                   if (valid_u) call insert_ens_obs(QKSWND_U_WIND_COMPONENT, uens, uerr, uqc, &
                                                    obs_info, VERTISSURFACE, ne, obs, seq)
                   if (valid_v) call insert_ens_obs(QKSWND_V_WIND_COMPONENT, vens, verr, vqc, &
                                                    obs_info, VERTISSURFACE, ne, obs, seq)
              case default ;
                   !not supported types
              end select

           endif  ! end if over "valid"

        enddo !end loop k level

     enddo obs_loop ! end loop over obs

999  print*, ista-1

     if(ista-1.eq.istamax)then
        print*,'you should give larger istamax'
     endif

!    write out the sequence
     call write_obs_seq(seq, obs_seq_out_filename)

     ! clean-up
     close(12)
     do n = 1,ne
        close(iunit(n))
     enddo

     deallocate(iunit)
     deallocate(edata1,eerror1,edata2,eerror2, &
              edata3,eerror3,edata4,eerror4,edata5,eerror5, &
              edata6,eerror6,edata7,eerror7,edata8,eerror8, &
              edata9,eerror9,euinv,evinv,etinv,eqinv)
     deallocate(eiqc1,eiqc2,eiqc3,eiqc4,eiqc5,eiqc6,eiqc7,eiqc8,eiqc9)
     deallocate(uens,vens,tens,qens,gpsens)


contains

subroutine get_ens_stat(edata,eqc,eerr,einv,esize,ecopy,valid)
    integer :: esize, i
    real(r8):: edata(esize), eerr(esize), einv(esize)
    integer :: eqc(esize)
    real(r8)    :: ecopy(esize*2+5)
    real(r8):: emean, evar, esp
    logical :: valid
    integer :: num_copies

    valid=.true.
    num_copies=esize*2+5

    ! initialize all copies except the 1st copy (obs value)
    do i=2,num_copies
       ecopy(i)=missing_r8
    enddo

    if (ecopy(1)==missing_v) then
       valid=.false.
       return
    endif

    do i=1,esize
       if (edata(i)==missing_v .or. eqc(i)==missing_qc .or. eerr(i)==missing_e) then
          valid=.false.
          return
       endif
       if (edata(i).ne.ecopy(1)) then  ! Actual ob value
          print*,'Ensemble obs value mismatch: ',ecopy(1), edata(i)
          valid=.false.
          return
       !  stop
       endif
       if (eerr(i).ne.eerr(1)) then ! Ob error
          print*,'Ensemble obs error mismatch: ',eerr(1), eerr(i)
          valid=.false.
          return
       !  stop
       endif
       ecopy(5+2*i-1)=edata(i)-einv(i)  ! Model-simulated observation
    enddo
    call compute_mean_var(ecopy(6:5+2*esize-1:2),esize,emean,evar)
    ecopy(2)=emean
    ecopy(4)=sqrt(evar)
    return
end subroutine get_ens_stat

subroutine compute_mean_var(evalue,esize,emean,evar)
    integer     :: esize
    real(r8)    :: evalue(esize), emean, evar
    emean=sum(evalue)/esize 
    evar=sum((evalue-emean)**2)/(esize-1)
end subroutine compute_mean_var

subroutine insert_ens_obs(obskind, ens_copy, obserr, obsqc, &
                          obs_info, which_vert, esize, obs, seq)
    type(obs_info_type) :: obs_info
    integer             :: obskind, which_vert, esize
    type(obs_type)      :: obs
    type(obs_sequence_type) :: seq
    real(r8)                :: obserr
    real(r8)                :: obsqc(1)
!   real(r8), dimension(esize+3) :: ens_copy
    real(r8), dimension(esize*2+5) :: ens_copy

    integer                 :: year, month, day, hours, minutes, seconds
    CHARACTER (len = 80)    :: dummy
    real(r8)                :: lat,lon
    real(r8)                :: gnx, gny, gnz, ds, htop, rfict
    character(len=6)        :: subset
    integer                 :: gpsref_key = 0 ! Initialize it, but it doesn't do anything

    ! time info
    read(obs_info%date_char, '(I4,5(A1,I2))')  &
        year, dummy, month, dummy, day, dummy, hours, dummy, minutes, dummy, seconds
    time = set_date(year, month, day, hours, minutes, seconds)
    call set_obs_def_time(obs_def, time)

    ! location info, location obs_def to be set later
    lat = xlat
    lon = xlon
    ! Dart longitude from 0 to 360
    if(lon < 0.0_r8) lon = lon + 360.0_r8
    if (which_vert.eq.VERTISSURFACE) then
       if (obs_info%elv == missing_v) return
       location = set_location(lon, lat, obs_info%elv, which_vert)
    else if (which_vert.eq.VERTISHEIGHT) then
       if (obs_info%h_lev == missing_v) return
       location = set_location(lon, lat, obs_info%h_lev, which_vert)
    else if (which_vert.eq.VERTISPRESSURE) then
       if (obs_info%p_lev == missing_v) return
       location = set_location(lon, lat, obs_info%p_lev, which_vert)
    end if

    call set_obs_def_location(obs_def, location)
    call set_obs_def_kind(obs_def, obskind)

    if (obs_info%platform(4:6).eq.'116' ) then ! GPSRO refractivity
         gnx = 0.0
         gny = 0.0
         gnz = 0.0
         ds = 0.0
         htop = 0.0
         rfict = 0.0
         subset = 'GPSREF'
         call set_gpsro_ref(gpsref_key, gnx, gny, gnz, rfict, ds, htop, subset)
         call set_obs_def_key(obs_def, gpsref_key)
    endif

    call set_obs_def_error_variance(obs_def, obserr*obserr)
    call set_obs_def(obs, obs_def)
    call set_obs_values(obs, ens_copy)
    call set_qc(obs,obsqc)
    call insert_obs_in_seq(seq,obs)

end subroutine insert_ens_obs

end program ens_obs_to_dart


