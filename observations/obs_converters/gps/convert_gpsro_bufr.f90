! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_gpsro_bufr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_gpsro_bufr - program that reads a COSMIC GPS observation 
!                        in BUFR format and writes the data to a DART 
!                        observation sequence file. 
!
!   For platforms with little-endian binary file format (e.g. Intel, AMD®, and
!   non-MIPS SGI processors), you may need to first run the byte-swapping program
!   DART/observations/NCEP/prep_bufr/exe/grabbufr.x over the bufr file 
!   before running this program. For the info, check
!   DART/observations/NCEP/prep_bufr/prep_bufr.html. 
!
!   created August 2015 Soyoung Ha, NCAR/MMM
!
!   Tested with gdas.gpsro.2012060100 for local operator only.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use          types_mod, only : r8
use   time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, set_time,&
                               increment_time, get_time, set_date, operator(-),  &
                               operator(>=), operator(<), operator(>), operator(<=), &
                               print_date, decrement_time
use      utilities_mod, only : initialize_utilities, find_namelist_in_file,    &
                               check_namelist_read, nmlfileunit, do_nml_file,   &
                               get_next_filename, error_handler, E_ERR, E_MSG, &
                               nc_check, find_textfile_dims, do_nml_term
use       location_mod, only : VERTISHEIGHT, set_location
use   obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq,       &
                               static_init_obs_sequence, init_obs, destroy_obs, &
                               write_obs_seq, init_obs_sequence, get_num_obs,   &
                               insert_obs_in_seq, destroy_obs_sequence,         &
                               set_copy_meta_data, set_qc_meta_data, set_qc,    & 
                               set_obs_values, set_obs_def
use   obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_type_of_obs, &
                               set_obs_def_error_variance, set_obs_def_location, &
                               set_obs_def_key
use    obs_def_gps_mod, only : set_gpsro_ref
use       obs_kind_mod, only : GPSRO_REFRACTIVITY
use  obs_utilities_mod, only : add_obs_to_seq


implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"


integer, parameter ::   num_copies = 1,   &   ! number of copies in sequence
                        num_qc     = 1        ! number of QC entries

character (len=129) :: msgstring, next_infile
character (len=80)  :: name
character (len=6)   :: subset
integer :: nlevels, nfiles, num_new_obs, oday, osec, &
           iyear, imonth, iday, ihour, imin, isec, obs_num, &
           io, iunit, filenum, dummy
logical :: file_exist, first_obs, did_obs, from_list = .false.
real(r8) :: oerr, qc, nx, ny, nz, &
            rfict, obsval, obs_val(1), qc_val(1)

!real(r8), allocatable :: lat(:), lon(:), hght(:), refr(:), azim(:), & 
!                         hghtp(:), refrp(:)

type(obs_def_type)      :: obs_def
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: obs_time, prev_time
type(time_type)         :: anal_time, window_min, window_max

integer :: gday, gsec, dsec, bsec, bday, esec, eday, num_excluded_bytime
logical :: use_bending_angle = .false.      ! Reserve for later.
!------------------------------------------------------------------------
! Array for bufr data
!------------------------------------------------------------------------
 character(len=80)  hdr1a
 character(len=10)  nemo
 character(len=8)   subsetb
 character(len=1)   cflg
 character(len=7)   obstype
 character(len=120) crecord

 integer, parameter:: unit_in=91
 integer, parameter:: unit_aux=92
 integer, parameter:: n1ahdr=10
 integer, parameter:: maxlevs=500    ! max number of observation levels (no longer namelist parameter)
 integer, parameter:: maxinfo=16
 integer, parameter:: mxib=31
 integer, dimension(50):: isats

 integer :: num_satid, nk
 integer :: ireadmg,ireadsb
 integer :: idate,iret,num_message,num_subset
 integer :: ikx, nlevs, nreps_ROSEQ1
 integer :: i,k,m,said,ptid
 integer :: nread,ndata,nprof_gps
 integer :: ibit(mxib),nib
 integer :: nsatid,nlines,istat
 integer :: isatid,isub,iuse

 real(r8) :: qfro, pcc, roc, impact, geopot
 real(r8) :: height,rlat,rlon,ref,referr,ref_pccf
 real(r8) :: freq_chk,freq,bend,bend_error,bend_pccf
 real(r8), dimension(n1ahdr):: hdr
 real(r8), dimension(50,maxlevs):: data1b
 real(r8), dimension(50,maxlevs):: data2a
 real(r8), dimension(maxlevs):: nreps_this_ROSEQ2

 integer, allocatable, dimension(:):: gps_satid
 !real(r8), allocatable, dimension(:,:):: cdata_all

 integer, parameter:: max_num_obs=200000    ! max number of observations 

 logical  lone, good

!------------------------------------------------------------------------
!  Declare namelist parameters
!------------------------------------------------------------------------
real(r8) :: ray_hbot =  3000.0_r8     ! lowest level (meter)
real(r8) :: ray_htop = 30000.0_r8     ! highest level (m) for observation
real(r8) :: ray_ds                    ! No longer used. Just keep for set_gpsro_ref now.
real(r8) :: obs_window_hr = -1.0      ! observation time window (+/-hours); A negative value means taking all the observations regardless of their report time.
real(r8) :: gsi_error_inv_factor = 1.0     ! scale factor for the inverse of gsi error
logical  :: if_global = .true.               ! Is your forecast model a regional or a global model?
logical  :: obs_error_in_gsi = .true.        ! obs error is overwritten as one percent of the refractivity obs value
logical  :: convert_to_geopotential_height = .false.  ! convert geometric height in bufr data to geopotential height?
logical  :: debug = .false.
character(len=128) :: gpsro_bufr_file     = 'gdas.gpsro.bufr'
character(len=128) :: gpsro_bufr_filelist = 'list_gpsro_bufr'
character(len=128) :: gpsro_out_file      = 'obs_seq.gpsro'
character(len=128) :: gpsro_aux_file      = 'convinfo.txt'         ! gps satellite info (from GSI)

namelist /convert_gpsro_bufr_nml/ ray_htop, ray_hbot, obs_window_hr, if_global, &
                                  obs_error_in_gsi, convert_to_geopotential_height, &
                                  gsi_error_inv_factor, &
                                  gpsro_bufr_file, gpsro_bufr_filelist, &
                                  gpsro_out_file,  gpsro_aux_file, debug
!------------------------------------------------------------------------
 nx        = 0.0_r8
 ny        = 0.0_r8
 nz        = 0.0_r8
 rfict     = 0.0_r8
 nprof_gps = 0

 data hdr1a / 'YEAR MNTH DAYS HOUR MINU PCCF ELRC SAID PTID GEODU' /
 data nemo /'QFRO'/

!------------------------------------------------------------------------
! Open and read the auxiliary info (mainly for satellite id to be used)
!------------------------------------------------------------------------
 open(unit_aux,file=gpsro_aux_file,action='read',form='formatted')

 rewind(unit_aux)
 nsatid=0
 nlines=0
 read1: do
        read(unit_aux,1030,iostat=istat,end=1130) cflg,obstype,crecord
        if (istat /= 0) exit
        nlines=nlines+1
        if(cflg == '!') cycle
        if(trim(obstype) /= 'gps') cycle    ! We process gps only.
        read(crecord,*) isatid,isub,iuse
        if (iuse < 0) cycle       ! Soyoung: -1 for monitoring. We skip it.
        nsatid=nsatid+1
        isats(nsatid)=isatid
 enddo read1
 1030 format(a1,a7,2x,a120)
 1130 continue
 if (istat>0) then
    write(*,*) '***ERROR*** error reading convinfo, istat=',istat
    stop
 endif

 num_satid=nsatid
 allocate(gps_satid(num_satid))
 gps_satid(1:num_satid)=isats(1:num_satid)
 write(*,'(A,20I5)') 'convert_gpsro_bufr: satellite id:',gps_satid

!------------------------------------------------------------------------
! initialize some values
!------------------------------------------------------------------------
obs_num = 1
qc = 0.0_r8
first_obs = .true.
call set_calendar_type(GREGORIAN)

!  read the necessary parameters from input.nml
call initialize_utilities()

call find_namelist_in_file("input.nml", "convert_gpsro_bufr_nml", iunit)
read(iunit, nml = convert_gpsro_bufr_nml, iostat = io)
call check_namelist_read(iunit, io, "convert_gpsro_bufr_nml")

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=convert_gpsro_bufr_nml)
if (do_nml_term()) write(     *     , nml=convert_gpsro_bufr_nml)

!  FIXME: We take all the native observation levels (instead of interpolating to
!  the levels at regular intervals. 
!  Then, how to count observation levels from each profile? (e.g., nlevels)
!  For now, we take maxlevs as num_new_obs, and are not sure if gpsro_bufr_filelist 
!  would work with it.
!nlevels = maxlevs

! namelist checks for sanity

! cannot have both a single filename and a list; the namelist must
! shut one off.
if (gpsro_bufr_file /= '' .and. gpsro_bufr_filelist /= '') then
  call error_handler(E_ERR, 'convert_gpsro_bufr',                     &
                    'One of gpsro_bufr_file or filelist must be NULL', &
                     source, revision, revdate)
endif
if (gpsro_bufr_filelist /= '') from_list = .true.

if ( obs_window_hr > 0.0 ) then
     print '(A,F5.2)', "Limit gpsro data within the time window (+/- hours): ", obs_window_hr
     print*, "Enter target assimilation time (gregorian day, second): "
     read*, gday,gsec

     call set_calendar_type(GREGORIAN)
     anal_time = set_time(gsec, gday)
     dsec = nint(obs_window_hr * 3600.)
     window_min = decrement_time(anal_time, dsec)
     window_max = increment_time(anal_time, dsec)
     call get_time(window_min,bsec,bday)
     call get_time(window_max,esec,eday)

     print '(2(A,2I8))', 'Observations are taken between ',bday,bsec,' and ',eday,esec
     num_excluded_bytime = 0   ! total number of obs beyond the time window
endif

! need to know a reasonable max number of obs that could be added here.
if (from_list) then
   call find_textfile_dims(gpsro_bufr_filelist, nfiles, dummy)
   num_new_obs = max_num_obs * nfiles
else
   num_new_obs = max_num_obs           !maxlevs
endif

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
inquire(file=gpsro_out_file, exist=file_exist)
if ( file_exist ) then

   print *, "found existing obs_seq file, appending to ", trim(gpsro_out_file)
   call read_obs_seq(gpsro_out_file, 0, 0, num_new_obs, obs_seq)

else

  print *, "no existing obs_seq file, creating ", trim(gpsro_out_file)
  print *, "max entries = ", num_new_obs
  call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)

  do k = 1, num_copies
    call set_copy_meta_data(obs_seq, k, 'COSMIC GPS observation')
  end do
  do k = 1, num_qc
    call set_qc_meta_data(obs_seq, k, 'COSMIC QC')
  end do

end if


did_obs = .false.

!------------------------------------------------------------------------
! main loop that does either a single file or a list of files
!------------------------------------------------------------------------

filenum = 1
fileloop: do      ! until out of files

   ! get the single name, or the next name from a list
   if (from_list) then 
      next_infile = get_next_filename(gpsro_bufr_filelist, filenum)
   else
      next_infile = gpsro_bufr_file
      if (filenum > 1) next_infile = ''
   endif
   if (next_infile == '') exit fileloop
  
   !------------------------------------------------------------------------
   ! Open file for input
   !------------------------------------------------------------------------
   open(unit_in,file=next_infile,action='read',form='unformatted')

   !------------------------------------------------------------------------
   ! Reading bufr data
   !------------------------------------------------------------------------
   call openbf(unit_in,'IN',unit_in)
   call datelen(10)
   call readmg(unit_in,subsetb,idate,iret)
   if (iret /= 0) goto 1010

   ! Big loop over the bufr file
   num_message=0
   msg_report: do while (ireadmg(unit_in,subsetb,idate) == 0)
     num_message = num_message + 1
     num_subset = 0
     write(*,'(I10,I5,a10)') idate,num_message,subsetb

     read_loop: do while (ireadsb(unit_in) == 0)
       num_subset = num_subset+1

       ! Extract header information
       call ufbint(unit_in,hdr,n1ahdr,1,iret,hdr1a)
       call ufbint(unit_in,qfro,1,1,iret,nemo)

       ! observation time in minutes
       iyear  = hdr(1) ! year
       imonth = hdr(2) ! month
       iday   = hdr(3) ! day
       ihour  = hdr(4) ! hour
       imin   = hdr(5) ! minute
       pcc=hdr(6)         ! profile per cent confidence
       roc=hdr(7)         ! Earth local radius of curvature
       said=hdr(8)        ! Satellite identifier (receiver)
       ptid=hdr(9)        ! Platform transmitter ID number

       write(*,*) iyear, imonth, iday, ihour, imin
       obs_time = set_date(iyear, imonth, iday, ihour, imin, 0)
       call get_time(obs_time,  osec, oday)

       ! Check the sateliite list in convinfo
       ikx = 0
       find_loop: do i = 1, num_satid
          if ( said == gps_satid(i) ) then
               ikx = 1
               exit find_loop
          endif
       enddo find_loop
       if (ikx==0) then
           cycle read_loop
       endif

       ! Check profile quality flags
       if ( ((said > 739).and.(said < 746)).or.(said == 820).or.(said == 786)) then  !CDAAC processing
           if(pcc==0.0) then
              write(*,*)'convert_gpsro_bufr: 1 bad profile said=',said,'ptid=',ptid,&
                  ' SKIP this report'
              cycle read_loop
           endif
        endif

        if ((said == 4).or.(said == 3).or.(said == 421).or.(said == 440).or.&
            (said == 821)) then ! GRAS SAF processing
           call upftbv(unit_in,nemo,qfro,mxib,ibit,nib)
           lone = .false.
           if(nib > 0) then
              do i=1,nib
                 if(ibit(i)== 6) then
                    lone = .true.
                    exit
                 endif
              enddo
           endif
           if(lone) then
              write(*,*) 'convert_gpsro_bufr: 2 bad profile said=',said,'ptid=',ptid,&
                         ' SKIP this report'
              cycle read_loop
           endif
        endif

       ! Read observations
       call ufbint(unit_in,nreps_this_ROSEQ2,1,maxlevs,nreps_ROSEQ1,'{ROSEQ2}')
       call ufbseq(unit_in,data1b,50,maxlevs,  nlevs,'ROSEQ1')    ! bending angle
       call ufbseq(unit_in,data2a,50,maxlevs,nlevels,'ROSEQ3')    ! refractivity
       if(nlevs/=nlevels) then
          write(*,'(A,2I5,2(A,I5),A)') 'convert_gpsro_bufr:  *WARNING* said,ptid=',said,ptid,&
                     ' with gps_bnd nlevs=', nlevs,&
                     ' and gps_ref nlevels=',nlevels,&
                     ' SKIP this report'
          cycle read_loop
       endif

       ! Exclude observations outside the time window.
       if ( obs_window_hr > 0.0 ) then
          if (obs_time <= window_min .or. obs_time > window_max ) then
              num_excluded_bytime = num_excluded_bytime + nlevels
              cycle read_loop
          end if
       end if
   
       ! Increment report counters
       nprof_gps = nprof_gps + 1      ! count reports in bufr file
       write(*,'(5(A8,I5))') ' nprof =',nprof_gps,' said =',said,' ptid=',ptid,&
                             ' nlevs =',nlevs,' nlevels =',nlevels

       !first_obs = .true.

       ! Loop over nlevels in profile
       obsloop2: do k = 1, nlevels

          rlat=data1b(1,k)      ! earth relative latitude (degrees)
          rlon=data1b(2,k)      ! earth relative longitude (degrees)
          height=data2a(1,k)    ! geometric height above geoid (m)
          ref=data2a(2,k)       ! refractivity obs (units of N)
          referr=data2a(4,k)    ! gps ref obs error (units of N)
          ref_pccf=data2a(6,k)  ! percent confidence (%)

          if(use_bending_angle) then        ! for later

          ! Loop over number of replications of ROSEQ2 nested inside ROSEQ1
          nk = nreps_this_ROSEQ2(k)
          do i=1,nk
              m=(6*i)-2
              freq_chk=data1b(m,k)          ! frequency (hertz)
              if(nint(freq_chk).ne.0) cycle ! do not want non-zero freq., go on
                                            ! to next replication of ROSEQ2
              freq=data1b(m,k)
              impact=data1b(m+1,k)      ! impact parameter (m)
              bend=data1b(m+2,k)        ! bending angle (rad)
              bend_error=data1b(m+4,k)  ! RMSE in bending angle (rad)
          enddo
          bend_pccf=data1b((6*nk)+4,k)  ! percent confidence for this ROSEQ1 replication

          endif     !if(use_bending_angle) then

          ! Preliminary (sanity) checks for bad and missing data
          good=.true.
          if ((abs(rlat)>90.0_r8).or.(abs(rlon)>360.0_r8)) good=.false.
          if ((ref>=1.e+9_r8).or.(ref<=0.0_r8))            good=.false.
          if ((height<=ray_hbot).or.(height>=ray_htop))    good=.false.
          ! Remove MetOP/GRAS ata below 8 km (following GSI)
          if ((height<=8000.0_r8).and.((said==4).or.(said==3))) good=.false.

          if(use_bending_angle) then
          if ((bend>=1.e+9_r8).or.(bend<=0.0).or.(impact>=1.e+9_r8).or.(impact<roc)) then
               good=.false.
          endif
          endif

          ! If observation is "good" print out sample data
          if(good) then

            obsval = ref    ! unit "N"

            !if(convert_to_geopotential_height) &
            geopot = compute_geopotential_height(height,rlat)
          
            if(obs_error_in_gsi) then
               !oerr   = 0.01_r8 * obsval
               oerr   = gsi_refractivity_error(height, rlat, if_global, gsi_error_inv_factor)
            else
               oerr   = 0.01_r8 * ref_obserr_kuo_percent(height * 0.001_r8) * obsval
            endif
            subset = 'GPSREF'
   
            if(use_bending_angle) then
               obsval = bend
               if((impact-roc) <= 10000.0_r8) then
                   oerr=(-bend*0.09/10000.0_r8)*(impact-roc)+bend*1.e-1_r8
               else
                   oerr=max(7.e-6_r8,bend*1.e-2_r8)
               endif
            endif

            if(rlon==360.0_r8) rlon=0.0_r8
            if(rlon < 0.0_r8)  rlon=rlon+360.0_r8

            if(debug.and.(k.lt.10)) then
               if(use_bending_angle) then
               write(*,'(A2,I5,4f9.2,f20.3)') 'k=',k,rlat,rlon,height,bend,oerr
               else
               write(*,'(A2,I5,6f9.2)') 'k=',k,rlat,rlon,height,ref,oerr,referr
               endif
            endif

            if(convert_to_geopotential_height) height = geopot
            call set_gpsro_ref(obs_num, nx, ny, nz, rfict, ray_ds, ray_htop, subset)
            call set_obs_def_location(obs_def,set_location(rlon,rlat,height,VERTISHEIGHT))
            call set_obs_def_type_of_obs(obs_def, GPSRO_REFRACTIVITY)
            call set_obs_def_time(obs_def, set_time(osec, oday))
            call set_obs_def_error_variance(obs_def, oerr * oerr)
            call set_obs_def_key(obs_def, obs_num)
            call set_obs_def(obs, obs_def)
   
            obs_val(1) = obsval
            call set_obs_values(obs, obs_val)
            qc_val(1)  = qc
            call set_qc(obs, qc_val)
       
            call add_obs_to_seq(obs_seq, obs, obs_time, prev_obs, prev_time, first_obs)
       
            obs_num = obs_num+1

            if (.not. did_obs) did_obs = .true.
   
          endif    !(good) then
       end do obsloop2                   ! End of k loop over nlevels

       write(*,*)
     enddo read_loop            ! subsets
   enddo msg_report             ! messages

1010 continue

   !------------------------------------------------------------------------
   ! Close unit to input file
   !------------------------------------------------------------------------
   call closbf(unit_in)
   write(6,1020)'convert_gpsro_bufr:  nprof_gps = ',nprof_gps,' obs_num=',obs_num
   if ( obs_window_hr > 0.0 ) write(6,*) 'num_excluded_bytime=', num_excluded_bytime
1020 format(2(A,i8))

  filenum = filenum + 1

end do fileloop

! clean up and loop if there is another input file
deallocate(gps_satid)

! done with main loop.  if we added any obs to the sequence, write it out.
if (did_obs) then
!print *, 'ready to write, nobs = ', get_num_obs(obs_seq)
   if (get_num_obs(obs_seq) > 0) &
      call write_obs_seq(obs_seq, gpsro_out_file)

   ! minor stab at cleanup, in the off chance this will someday get turned
   ! into a subroutine in a module.  probably not all that needs to be done,
   ! but a start.
   call destroy_obs(obs)
! if obs == prev_obs then you can't delete the same obs twice.
! buf if they differ, then it's a leak.  for now, don't delete prev
! since the program is exiting here anyway.
   !call destroy_obs(prev_obs)
   if (get_num_obs(obs_seq) > 0) call destroy_obs_sequence(obs_seq)
endif

! END OF MAIN ROUTINE

contains

! local subroutines/functions follow

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
!
!   gsi_refractivity_error - function that computes the observation 
!                            error as in GSI
!
!    Input:
!    H   -- input real value geometric height [m]
!    lat -- latitude in degree 
!    is_it_global -- is your forecast model regional or global? (T or F)
!
!    Output:
!    gsi_refractivity_error -- output refractivity observation error [N]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function gsi_refractivity_error(H, lat, is_it_global, factor)
 real(r8), intent(in)  :: H
 real(r8), intent(in)  :: lat
 real(r8), intent(in)  :: factor
 logical,  intent(in)  :: is_it_global
 real(r8)              :: gsi_refractivity_error

 real(r8) :: zkm, rerr
 integer  :: kk
 
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

 gsi_refractivity_error = 1./(abs(exp(rerr))*factor)
 !gsi_refractivity_error = 1./abs(exp(rerr))

 return
end function gsi_refractivity_error

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
