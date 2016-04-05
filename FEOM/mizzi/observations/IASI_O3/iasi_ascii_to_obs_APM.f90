! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
program iasi_ascii_to_obs
!
! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   iasi_ascii_to_obs - a program that only needs minor customization to read
!      in a iasi_ascii-based dataset - either white-space separated values or
!      fixed-width column data.
!      
!     this is work in progress. IASI dataset are in HDF format. I do not
!     have HDF libraries for now, so Gabi Pfister reads the hdf file in IDL and
!     did some processing before she dumped the data in ascii. what you are
!     reading here is a 'processed' dataset of IASI O3
!
!     created 29 Mar 2010   nancy collins NCAR/IMAGe
!     modified 16 Mar 2012  ave arellano UArizona
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   use         types_mod, only : r8, PI, DEG2RAD
   use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                                 open_file, close_file
   use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                                 operator(>=), increment_time, get_time, &
                                 operator(-), GREGORIAN, operator(+), print_date
   use      location_mod, only : VERTISUNDEF, VERTISHEIGHT,VERTISPRESSURE
   use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                                 static_init_obs_sequence, init_obs, write_obs_seq, & 
                                 init_obs_sequence, get_num_obs, & 
                                 set_copy_meta_data, set_qc_meta_data
   use      obs_kind_mod, only : IASI_O3_RETRIEVAL
!
   implicit none
!
   character(len=64), parameter :: iasi_ascii_input_file = 'iasi_asciidata.input'
   character(len=64), parameter :: obs_out_file    = 'iasi_obs_seq.out'
!
   logical, parameter :: debug = .true.  ! set to .true. to print info
!
   character (len=84) :: input_line
   character (len=10) :: otype_char
!
   integer :: n, i, j, oday, osec, rcio, iunit, otype, ilev
   integer :: year, month, day, hour, minute, second
   integer :: num_copies, num_qc, max_obs
!           
   logical  :: file_exist, first_obs
!
   integer, parameter :: nlevels = 40
   integer  :: qc_count
   real(r8) :: iasi_col, iasi_vmr, sza, cloud, dfs100, dfs300
   real(r8) :: lat, lon, vert, retlev2, qc, rnlev_use
   real(r8) :: altretlev(nlevels) = 0.0_r8
   real(r8) :: akcol(nlevels,nlevels) = 0.0_r8
!
   real(r8) :: iasi_err(nlevels) = 0.0_r8
   real(r8) :: apcol(nlevels) = 0.0_r8
   real(r8) :: iasi_col_prof(nlevels) = 0.0_r8
   real(r8) :: iasi_vmr_prof(nlevels) = 0.0_r8
   real(r8) :: aircol_val(nlevels) = 0.0_r8
   real(r8) :: apcol_val(nlevels) = 0.0_r8
   real(r8) :: QOpt_prof(nlevels) = 0.0_r8
   real(r8) :: seconds, trm
   integer  :: sellev, nlev_use
!
   type(obs_sequence_type) :: obs_seq
   type(obs_type)          :: obs, prev_obs
   type(time_type)         :: comp_day0, time_obs, prev_time
!
! Start of executable code
   call initialize_utilities('iasi_ascii_to_obs')
!
! Time setup
   call set_calendar_type(GREGORIAN)
! 
! Open IASI ascii input file
   iunit = open_file(iasi_ascii_input_file, 'formatted', 'read')
   if (debug) print *, 'opened input file ' // trim(iasi_ascii_input_file)
!
! Each observation in this series will have a single observation value 
! and a quality control flag.  the max possible number of obs needs to
! be specified but it will only write out the actual number created.
   max_obs    = 10000000
   num_copies = 1
   num_qc     = 1
!
! Call the initialization code, and initialize two empty observation types
   call static_init_obs_sequence()
   call init_obs(obs,      num_copies, num_qc)
   call init_obs(prev_obs, num_copies, num_qc)
   first_obs = .true.
!
! Create a new, empty obs_seq file.  you must give a max limit
! on number of obs.  increase the size if too small.
   call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)
!
! The first one needs to contain the string 'observation' and the
! second needs the string 'QC'.
   call set_copy_meta_data(obs_seq, 1, 'IASI O3 observation')
   call set_qc_meta_data(obs_seq, 1, 'IASI O3 QC index')
!
! If you want to append to existing files (e.g. you have a lot of
! small iasi_ascii files you want to combine), you can do it this way,
! or you can use the obs_sequence_tool to merge a list of files 
! once they are in DART obs_seq format.
!
! Set the DART data quality control.   0 is good data. 
! increasingly larger QC values are more questionable quality data.
   qc = 0.0_r8
   qc_count = 0
!
! Loop for reading IASI ascii data file
   obsloop: do    ! No end limit - have the loop break when input ends
!
! Read in a line from the iasi_ascii file.   What you need to create an obs:
!  location: lat, lon, and height in pressure or meters
!  time: when the observation was taken
!  type: from the DART list of obs types
!  error: very important - the instrument error plus representativeness error
!        (see html file for more info)
!  averaging kernel and a priori profile
!  for now, we chose 2 'retrieval levels' corresponding to highest sensitivity
!  assume to be independent from each other
!  assume here a line is a type (1/2), location, time, value, obs error
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! READ IASI ASCII CODE BLOCK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! READ RECORD ONE
      read(iunit, "(A)", iostat=rcio) input_line
!      print *, input_line
      if (rcio /= 0) then 
         if (debug) print *, 'got bad read code from input file, rcio = ', rcio
         exit obsloop
      endif
!
! Convert to date/time data
      read(input_line(1:10), *, iostat=rcio) otype_char
      if (rcio /= 0) then 
         if (debug) print *, 'got bad read code trying to get obs type, rcio = ', rcio
         exit obsloop
      endif
!   
! otype is fixed to 1 for IASI O3
      if (debug) print *, 'next observation type = ', otype_char
      read(input_line(12:15), *, iostat=rcio) year
      read(input_line(16:17), *, iostat=rcio) month
      read(input_line(18:19), *, iostat=rcio) day
!
! READ RECORD TWO
      read(iunit,"(9F14.4)", iostat=rcio) seconds, lat, lon, sza, cloud, retlev2, &
      rnlev_use, dfs100, dfs300 
!
! Check lon and convert to 0-360
      if ( lon > 180.0_r8 .or. lon < -180.0_r8 ) cycle obsloop
      if ( lon < 0.0_r8 )  lon = lon + 360.0_r8
!
! Assign date/time data
      nlev_use = int(rnlev_use)
      retlev2 = retlev2*1000.0_r8
      hour =  seconds/3600
      minute = (seconds - hour*3600)/60
      second = (seconds - hour*3600 - minute*60)
!      print *, "APM(2): hour, minute, second = ",year,month,day,hour,minute,second
!
! Put date/time data into DART format
      time_obs = set_date(year, month, day, hour, minute, second)
      if (debug) call print_date(time_obs, 'next obs time is')
!
! Put day and time into Gregorian format
      call get_time(time_obs, osec, oday)
!
! READ RECORD THREE
      read(iunit, *, iostat=rcio) (altretlev(i),i=1,nlev_use)
      print *, "APM(3): altretlev = ",(altretlev(i),i=1,nlev_use)
      altretlev(:) = altretlev(:)*1000.0_r8
!
! READ RECORD FOUR
      read(iunit,*, iostat=rcio) iasi_col, iasi_vmr
!      print *, "APM(4): iasi_col,iasi_vmr = ", iasi_col, iasi_vmr
!
! READ RECORD FIVE
      read(iunit,*, iostat=rcio) (iasi_err(i),i=1,nlev_use)
      print *, "APM(5):iasi_err = ",(iasi_err(i),i=1,nlev_use) 
!
! READ RECORD SIX
! APM: This read matches the IDL (column,row) convention
!      The averaging kernel dot product sum over the rows (as opposed to over
!      the columns as in MOPITT).
      read(iunit,*, iostat=rcio) ((akcol(i,j),j=1,nlev_use),i=1,nlev_use)
      do i=1,nlev_use
         print *, "APM(6):akcol =  ", (akcol(j,i),j=1,nlev_use)
      enddo
!
! READ RECORD SEVEN
      read(iunit,*, iostat=rcio) (apcol(i),i=1,nlev_use)
      print *, "APM(7): apcol = ",(apcol(i),i=1,nlev_use)
!
! Calculate (I-A)xa
      apcol_val(:) = 0.0_r8
      do j = 1, nlev_use
         do i = 1, nlev_use
            if(i.eq.j) then
               trm = (1. - akcol(i,j)) * apcol(i)
            else
               trm = -1. * akcol(i,j) * apcol(i)
            endif
            apcol_val(j) = apcol_val(j) + trm
         enddo
      enddo
!
! READ RECORD EIGHT
      read(iunit,*, iostat=rcio) (iasi_col_prof(i), i=1,nlev_use)
      print *, "APM(8): iasi_col_prof  = ",(iasi_col_prof(i),i=1,nlev_use)
!
! READ RECORD NINE
      read(iunit,*, iostat=rcio) (iasi_vmr_prof(i),i=1,nlev_use)
      if (rcio /= 0) then 
         if (debug) print *, 'got bad read code getting rest of iasi o3 obs, rcio = ', rcio
         exit obsloop
      endif
      print *, "APM(9): iasi_vmr_prof = ",(iasi_vmr_prof(i),i=1,nlev_use)
! 
! Convert error from % to retrievel error
      iasi_err(:) = iasi_err(:)*iasi_col
!
! Calculate air column
      do ilev = 1, nlev_use
         if (iasi_vmr_prof(ilev)>0.0_r8) then
            aircol_val(ilev) = iasi_col_prof(ilev)/iasi_vmr_prof(ilev)
         else
            aircol_val(ilev) = 0.0_r8
         endif
      enddo
!
! Calculate QOpt Retriaval
      QOpt_prof(:) = iasi_col_prof(:) - apcol_val(:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! WRITE DATA TO OBS_SEQ FILE  CODE BLOCK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Loop over all IASI vertical levels
      do ilev = 1, nlev_use
         qc_count = qc_count + 1
!
! This example assumes there is an obs type, where otype=1 is
! a temperature measured in height, and if otype=2, there's a wind
! speed and direction and height is pressure.  any kind of observation
! can use any of the vertical types; this is just an example.
!
! Fixed otype=1 for IASI O3
! No height since it is a column integrated quantity
!
! Make an obs derived type, and then add it to the sequence
         call create_3d_obs(lat, lon, altretlev(ilev), VERTISHEIGHT, QOpt_prof(ilev), &
                            IASI_O3_RETRIEVAL, iasi_err(ilev), oday, osec, qc, obs, &
                            akcol(1:nlev_use,ilev), apcol_val(ilev), altretlev(1:nlev_use), &
                            aircol_val(1:nlev_use), qc_count, nlev_use)
         call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
!
         if (debug) print *, 'added iasi obs to output seq'
      enddo 
   end do obsloop
!
! If we added obs to the sequence, write it out to a file
!   print *,"APM: obsloop exit num_obs ",get_num_obs(obs_seq)
   if ( get_num_obs(obs_seq) > 0 ) then
      if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
      call write_obs_seq(obs_seq, obs_out_file)
   endif
!
! End of main program
   call finalize_utilities()
!
   contains
!
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
!    oerr  - observation error
!    day   - gregorian day
!    sec   - gregorian second
!    qc    - quality control value
!    obs   - observation type
!
!    metadata (AFAJ)
!    akcol - averaging kernel
!    apcol_val - a priori sub-column (I-A)xa
!    nlevels - number of levels
!    aircol_val - air sub column
!    qc_count - obs count (key)
!
!     created Oct. 2007 Ryan Torn, NCAR/MMM
!     adapted for more generic use 11 Mar 2010, nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine create_3d_obs(lat, lon, vval, vkind, obsv, okind, oerr, day, sec, qc, &
                            obs, akcol, apcol_val, altretlev, aircol_val, qc_count, nlev_use)
      use        types_mod, only : r8
      use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_kind, &
                             set_obs_def_error_variance, set_obs_def_location, set_obs_def_key
      use obs_def_iasi_O3_mod, only : set_obs_def_iasi_o3
      use obs_sequence_mod, only : obs_type, set_obs_values, set_qc, set_obs_def
      use time_manager_mod, only : time_type, set_time
      use     location_mod, only : set_location
!
      integer, parameter            :: nlevels = 40
      integer,        intent(in)    :: nlev_use
      integer,        intent(in)    :: okind, vkind, day, sec
      integer,        intent(in)    :: qc_count
      real(r8),       intent(in)    :: lat, lon, vval, obsv, oerr, qc
      real(r8),       intent(in)    :: akcol(nlevels)
      real(r8),       intent(in)    :: altretlev(nlevels)
      real(r8),       intent(in)    :: aircol_val(nlevels)
      real(r8),       intent(in)    :: apcol_val
      type(obs_type), intent(inout) :: obs
!
      real(r8)           :: obs_val(1), qc_val(1)
      type(obs_def_type) :: obs_def
!
      call set_obs_def_location(obs_def, set_location(lon, lat, vval, vkind))
      call set_obs_def_kind(obs_def, okind)
      call set_obs_def_time(obs_def, set_time(sec, day))
      call set_obs_def_error_variance(obs_def, oerr * oerr)
      call set_obs_def_iasi_o3(qc_count, apcol_val, altretlev(1:nlev_use), akcol(1:nlev_use), &
      aircol_val(1:nlev_use), nlev_use)
      call set_obs_def_key(obs_def, qc_count)
!
      obs_val(1) = obsv
      call set_obs_values(obs, obs_val)
      qc_val(1)  = qc
      call set_qc(obs, qc_val)
!
      call set_obs_def(obs, obs_def)
   end subroutine create_3d_obs
!
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
      use types_mod, only        : r8
      use obs_sequence_mod, only : obs_sequence_type, obs_type, insert_obs_in_seq
      use time_manager_mod, only : time_type, operator(>=)
!
      type(obs_sequence_type), intent(inout) :: seq
      type(obs_type),          intent(inout) :: obs, prev_obs
      type(time_type),         intent(in)    :: obs_time
      type(time_type),         intent(inout) :: prev_time
      logical,                 intent(inout) :: first_obs
!
! insert(seq,obs) always works (i.e. it inserts the obs in
! proper time format) but it can be slow with a long file.
! supplying a previous observation that is older (or the same
! time) as the new one speeds up the searching a lot.
!
      if(first_obs) then
         call insert_obs_in_seq(seq, obs)
         first_obs = .false.
      else               
         if(obs_time >= prev_time) then
            call insert_obs_in_seq(seq, obs, prev_obs)
         else
            call insert_obs_in_seq(seq, obs)
         endif
      endif
!
! update for next time
      prev_obs = obs
      prev_time = obs_time
   end subroutine add_obs_to_seq
end program iasi_ascii_to_obs
