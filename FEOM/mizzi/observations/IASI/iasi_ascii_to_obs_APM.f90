! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
program iasi_ascii_to_obs
!
! <next few lines under version control, do not edit>
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/trunk/observations/iasi_ascii/iasi_ascii_to_obs.f90 $
! $Id: iasi_ascii_to_obs.f90 4865 2011-04-18 20:56:25Z nancy $
! $Revision: 4865 $
! $Date: 2011-04-18 13:56:25 -0700 (Mon, 18 Apr 2011) $
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
integer :: n, i, oday, osec, rcio, iunit, otype, ilev
integer :: year, month, day, hour, minute, second
integer :: num_copies, num_qc, max_obs
!           
logical  :: file_exist, first_obs
!
integer, parameter :: nlevels = 40
integer  :: qc_count
real(r8) :: iasi_err(nlevels) = 0.0_r8
real(r8) :: iasi_col, iasi_vmr, sza, cloud, dfs100, dfs300
real(r8) :: lat, lon, vert, retlev2, qc, rnlev_use
real(r8) :: altretlev2(nlevels) = 0.0_r8
real(r8) :: akcol(nlevels,nlevels) = 0.0_r8
!
real(r8) :: apcol(nlevels) = 0.0_r8
real(r8) :: iasi_col_prof(nlevels) = 0.0_r8
real(r8) :: iasi_vmr_prof(nlevels) = 0.0_r8
real(r8) :: aircol_val(nlevels) = 0.0_r8
real(r8) :: apcol_val(nlevels)
real(r8) :: seconds
integer  :: sellev, nlev_use
double precision, allocatable, dimension(:) :: egnval
real(r8), allocatable, dimension(:) :: egnval_A, egnval_R
real(r8), allocatable, dimension(:) :: iasi_col_prof_rot, apcol_rot, iasi_err_rot
real(r8), allocatable, dimension(:) :: aircol_val_rot, apcol_val_rot
double precision, allocatable, dimension(:,:) :: egnvec
real(r8), allocatable, dimension(:,:) :: R, A_rot, R_rot, egnvec_A, egnvec_R
double precision, allocatable, dimension (:) :: zwrk
integer :: lwrk,info,idfs_A,idfs_R   
integer  :: iretlev2, j
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
call set_copy_meta_data(obs_seq, 1, 'observation')
call set_qc_meta_data(obs_seq, 1, 'Data QC')
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
   print *, input_line
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
   hour =  seconds/3600
   minute = (seconds - hour*3600)/60
   second = (seconds - hour*3600 - minute*60)
!   print *, "APM: hour, minute, second = ",hour,minute,second
!
! Put date/time data into DART format
   time_obs = set_date(year, month, day, hour, minute, second)
   if (debug) call print_date(time_obs, 'next obs time is')
!
! Put day and time into Gregorian format
   call get_time(time_obs, osec, oday)
!
! READ RECORD THREE
   read(iunit, *, iostat=rcio) (altretlev2(i),i=1,nlev_use)
!   print *, "APM: altretlev2(1),altretlev2(nlev_use) = ",altretlev2(1),altretlev2(nlev_use)
   altretlev2(1:nlev_use) = altretlev2(1:nlev_use)*1000.0_r8
! 
! Asign sellev
   sellev = int(retlev2+.01)     
!   print *, "APM: sellev = ", sellev
!
! READ RECORD FOUR
   read(iunit,*, iostat=rcio) iasi_col, iasi_vmr
!   print *, "APM: iasi_col,iasi_vmr = ", iasi_col, iasi_vmr
!
! READ RECORD FIVE
   read(iunit,*, iostat=rcio) (iasi_err(i),i=1,nlev_use)
!   print *, "APM: iasi_err(1),iasi_err(nlev_use) = ",iasi_err(1),iasi_err(nlev_use)
!
! READ RECORD SIX
   read(iunit,*, iostat=rcio) ((akcol(i,j),i=1,nlev_use),j=1,nlev_use)
!   print *, "APM: akcol(1,1),akcol(2,1),akcol(3,1) = ",akcol(1,1),akcol(2,1),akcol(3,1)
!
! READ RECORD SEVEN
   read(iunit,*, iostat=rcio) (apcol(i),i=1,nlev_use)
!   print *, "APM: apcol(1),apcol(nlev_use) = ",apcol(1),apcol(nlev_use)
!
! Transform to (I-A)xa
   apcol_val(:) = 0.0_r8
   do i = 1, nlev_use
      do j = 1, nlev_use
        apcol_val(j) = apcol_val(j) + akcol(i,j)*apcol(j)
      enddo
   enddo
   apcol_val(:) = apcol(:) - apcol_val(:)
!
! READ RECORD EIGHT
   read(iunit,*, iostat=rcio) (iasi_col_prof(i), i=1,nlev_use)
!   print *, "APM: iasi_col_prof(1),iasi_col_prof(nlev_use) = ",iasi_col_prof(1),iasi_col_prof(nlev_use)
!
! READ RECORD NINE
!   read(iunit,*, iostat=rcio) (iasi_vmr_prof(i), i=1,nlev_use)
!   print *, "APM: iasi_vmr_prof(1),iasi_vmr_prof(nlev_use) = ",iasi_vmr_prof(1),iasi_vmr_prof(nlev_use)
! 
! Convert the error term from % to retrievel error
   iasi_err(:) = iasi_err(:)*iasi_col_prof(:)
!
! Scale retrieval by VMR
   do ilev = 1, nlevels
      !print *, iasi_vmr_prof(i), iasi_col_prof(i)
      if (iasi_vmr_prof(ilev)>0.0_r8) then
      aircol_val(ilev) = iasi_col_prof(ilev)/iasi_vmr_prof(ilev)
      else
      aircol_val(ilev) = 0.0_r8
      endif
   enddo
!
! Check for read exist code error
   if (rcio /= 0) then 
      if (debug) print *, 'got bad read code getting rest of iasi o3 obs, rcio = ', rcio
      exit obsloop
   endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! CALCULATE EIGENSTRUCTURES CODE BLOCK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Allocate arrays 
   allocate (zwrk(1))
   allocate (R(nlev_use,nlev_use))
   allocate (egnval(nlev_use),egnval_A(nlev_use),egnval_R(nlev_use))
   allocate (egnvec(nlev_use,nlev_use),egnvec_A(nlev_use,nlev_use),egnvec_R(nlev_use,nlev_use))
!
! Calculate eigenvalues and eigenvector for A
   egnvec(:,:)= akcol(:,:) 
   call dsyev('V', 'L', nlev_use, egnvec, nlev_use, egnval, zwrk, -1, info)
   lwrk = int(zwrk(1))
   deallocate (zwrk)
   allocate (zwrk(lwrk))
   call dsyev('V', 'L', nlev_use, egnvec, nlev_use, egnval, zwrk, lwrk, info)
!
! Reorder eigenvalues
   do i=1,nlev_use
      egnval_A(i) = egnval(nlev_use-i+1)
      print *, "APM: i,egnval_A(i) = ",i,egnval_A(i)
   enddo
!
! Reorder eigenvectors
   do i=1,nlev_use
      egnvec_A(:,i) = egnvec(:,nlev_use-i+1)
   enddo
!
! Calculate R matrix
   do i = 1, nlev_use
      do j = 1, nlev_use
         R(i,j) = iasi_err(i)*iasi_err(j)
!         print *, i,j,R(i,j)
     enddo
   enddo
!
! Calculate eigenvalues and eigenvector for R
   egnvec(:,:)= R(:,:)
   deallocate (zwrk)
   allocate (zwrk(1))
   call dsyev('V', 'L', nlev_use, egnvec, nlev_use, egnval, zwrk, -1, info)
   lwrk = int(zwrk(1))
   deallocate (zwrk)
   allocate (zwrk(lwrk))
   call dsyev('V', 'L', nlev_use, egnvec, nlev_use, egnval, zwrk, lwrk, info)
   deallocate(R,zwrk)
!
! Reorder eigenvalues
   do i=1,nlev_use
      egnval_R(i) = egnval(nlev_use-i+1)
      print *, "APM: i,egnval_R(i) = ",i,egnval_R(i)
   enddo
!
! Reorder eigenvectors
   do i=1,nlev_use
      egnvec_R(:,i) = egnvec(:,nlev_use-i+1)
   enddo
   deallocate(egnval,egnvec)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! CALCULATE ROTATION CODE BLOCK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Truncate the eigenspectrum
   idfs_A=3
   idfs_R=3
   do i=idfs_A+1,nlev_use
      egnval_A(i)=0.
   enddo
   do j=idfs_A+1,nlev_use
      do i=1,nlev_use
         egnvec_A(i,j)=0.
      enddo
   enddo   
   do i=idfs_R+1,nlev_use
      egnval_R(i)=0.
   enddo
   do j=idfs_R+1,nlev_use
      do i=1,nlev_use
         egnvec_R(i,j)=0.
      enddo
   enddo   
!
! Calculate A_rot
   allocate(A_rot(nlev_use,nlev_use),R_rot(nlev_use,nlev_use))
   allocate(iasi_col_prof_rot(nlev_use),apcol_rot(nlev_use),iasi_err_rot(nlev_use))
   allocate(aircol_val_rot(nlev_use),apcol_val_rot(nlev_use))

   A_rot(:,:)=0.
   do i=1,nlev_use
      do j=1,nlev_use
         A_rot(i,j)=A_rot(i,j)+egnvec_A(i,j)*akcol(j,i)
      enddo
   enddo
!
! Calculate iasi_col_prof_rot, apcol_rot, and iasi_err_rot
   iasi_col_prof_rot(:)=0.
   apcol_rot(:)=0.
   iasi_err_rot(:)=0.
   do i=1,nlev_use
      do j=1,nlev_use
         iasi_col_prof_rot(i)=iasi_col_prof_rot(i)+egnvec_A(i,j)*iasi_col_prof(j)
         apcol_rot(i)=apcol_rot(i)+egnvec_A(i,j)*apcol(j)
         iasi_err_rot(i)=iasi_err_rot(i)+egnvec_A(i,j)*iasi_err(j)
      enddo
   enddo
!
! Scale rotated retrieval by VMR
   do ilev = 1, nlevels
      !print *, iasi_vmr_prof(i), iasi_col_prof(i)
      if (iasi_vmr_prof(ilev)>0.0_r8) then
      aircol_val_rot(ilev) = iasi_col_prof_rot(ilev)/iasi_vmr_prof(ilev)
      else
      aircol_val_rot(ilev) = 0.0_r8
      endif
   enddo
!
! Transform to rotated (I-A)xa
   apcol_val_rot(:) = 0.0_r8
   do i = 1, nlev_use
      do j = 1, nlev_use
        apcol_val_rot(j) = apcol_val_rot(j) + A_rot(i,j)*apcol(j)
      enddo
   enddo
   apcol_val_rot(:) = apcol_rot(:) - apcol_val_rot(:)
   deallocate(egnval_A,egnvec_A,egnval_R,egnvec_R)
stop
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! WRITE DATA TO OBS_SEQ FILE  CODE BLOCK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
   qc_count = qc_count + 1
   call create_3d_obs(lat, lon, retlev2, VERTISHEIGHT, iasi_col_prof_rot(sellev), &
                      IASI_O3_RETRIEVAL, iasi_err_rot(sellev), oday, osec, qc, obs, &
                      A_rot(sellev,:), apcol_val_rot(sellev), altretlev2, &
                      aircol_val_rot, qc_count, nlev_use, apcol_rot)

   call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
   deallocate(iasi_col_prof_rot,apcol_rot,iasi_err_rot,aircol_val_rot,apcol_val_rot)
   if (debug) print *, 'added iasi obs to output seq'
end do obsloop
!
! If we added obs to the sequence, write it out to a file
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
subroutine create_3d_obs(lat, lon, vval, vkind, obsv, okind, oerr, day, sec, qc, obs, akcol, apcol_val, altretlev, aircol_val, qc_count, nlev_use, apcol)
use        types_mod, only : r8
use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_kind, &
                             set_obs_def_error_variance, set_obs_def_location, set_obs_def_key
use obs_def_iasi_mod, only : set_obs_def_iasi_o3
use obs_sequence_mod, only : obs_type, set_obs_values, set_qc, set_obs_def
use time_manager_mod, only : time_type, set_time
use     location_mod, only : set_location

 integer,        intent(in)    :: okind, vkind, day, sec
 real(r8),       intent(in)    :: lat, lon, vval, obsv, oerr, qc
 type(obs_type), intent(inout) :: obs
 integer, parameter            :: nlevels = 40
 real(r8),       intent(in)    :: akcol(nlevels)
 real(r8),       intent(in)    :: altretlev(nlevels)
 real(r8),       intent(in)    :: apcol(nlevels)
 real(r8),       intent(in)    :: aircol_val(nlevels)
 real(r8),       intent(in)    :: apcol_val ! aircol_val
 integer,        intent(in)    :: qc_count
 integer,        intent(in)    :: nlev_use

 real(r8)           :: obs_val(1), qc_val(1)
 type(obs_def_type) :: obs_def

call set_obs_def_location(obs_def, set_location(lon, lat, vval, vkind))
call set_obs_def_kind(obs_def, okind)
call set_obs_def_time(obs_def, set_time(sec, day))
call set_obs_def_error_variance(obs_def, oerr * oerr)
call set_obs_def_iasi_o3(qc_count, akcol(1:nlev_use), apcol_val, altretlev(1:nlev_use), aircol_val(1:nlev_use), nlev_use, apcol(1:nlev_use))
call set_obs_def_key(obs_def, qc_count)
!call set_obs_def(obs, obs_def)

obs_val(1) = obsv
call set_obs_values(obs, obs_val)
qc_val(1)  = qc
call set_qc(obs, qc_val)

call set_obs_def(obs, obs_def)
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
 use        types_mod, only : r8
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
prev_obs = obs
prev_time = obs_time

end subroutine add_obs_to_seq

end program iasi_ascii_to_obs
