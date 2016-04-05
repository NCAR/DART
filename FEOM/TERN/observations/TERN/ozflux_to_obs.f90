! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program ozflux_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   ozflux_to_obs - a program that only needs minor customization to read
!      in a text-based dataset - either white-space separated values or
!      fixed-width column data.
!
!     created 3 May 2012   Tim Hoar NCAR/IMAGe
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, MISSING_R8

use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              register_module, error_handler, E_MSG, E_ERR, &
                              open_file, close_file, do_nml_file, do_nml_term, &
                              check_namelist_read, find_namelist_in_file, &
                              nmlfileunit

use  time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                              set_date, set_time, get_time, print_time, &
                              print_date, operator(-), operator(+), operator(>), &
                              operator(<), operator(==), operator(<=), operator(>=)

use      location_mod, only : VERTISHEIGHT

use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data

use      obs_kind_mod, only : TOWER_SENSIBLE_HEAT_FLUX, &
                              TOWER_NETC_ECO_EXCHANGE,  &
                              TOWER_LATENT_HEAT_FLUX, &
                              LEAF_AREA_INDEX

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! Namelist with default values
!-----------------------------------------------------------------------

character(len=128) :: text_input_file = 'textdata.input'
character(len=128) :: obs_out_file    = 'obs_seq.out'
real(r8)           :: timezoneoffset  = -1.0_r8
real(r8)           :: latitude        = -1.0_r8
real(r8)           :: longitude       = -1.0_r8
real(r8)           :: elevation       = -1.0_r8
real(r8)           :: flux_height     = -1.0_r8
real(r8)           :: maxgoodqc       = 3.0_r8
logical            :: verbose         = .TRUE.

namelist /ozflux_to_obs_nml/ text_input_file, obs_out_file, &
             timezoneoffset, latitude, longitude, elevation, &
             flux_height, maxgoodqc, verbose

!-----------------------------------------------------------------------
! globally-scoped variables
!-----------------------------------------------------------------------

character(len=4000)     :: input_line, bigline
character(len=256)      :: string1, string2, string3
integer                 :: iline, nlines, nwords
logical                 :: first_obs
integer                 :: oday, osec, rcio, iunit
integer                 :: num_copies, num_qc, max_obs
real(r8)                :: oerr, qc
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: prev_time, offset
real(r8), parameter     :: umol_to_gC = (1.0_r8/1000000.0_r8) * 12.0_r8

type towerdata
  type(time_type)   :: time_obs
  character(len=20) :: datestring  = 'DT'
  character(len=20) :: neestring   = 'Fc'             ! CarbonFlux
  character(len=20) :: neeQCstring = 'Fc_QCFlag'
  character(len=20) :: lestring    = 'Fe'             ! LatentHeat
  character(len=20) :: leQCstring  = 'Fe_QCFlag'
  character(len=20) :: hstring     = 'Fh'             ! SensibleHeat
  character(len=20) :: hQCstring   = 'Fh_QCFlag'
  character(len=20) :: laistring   = 'Lai_1km_new'
  character(len=20) :: gppstring   = 'GPP_CABLE'
  integer  :: dateindex
  integer  :: neeindex
  integer  :: neeQCindex
  integer  :: leindex
  integer  :: leQCindex
  integer  :: hindex
  integer  :: hQCindex
  integer  :: laiindex
  integer  :: gppindex
  real(r8) :: nee
  integer  :: neeQC
  real(r8) :: le
  integer  :: leQC
  real(r8) :: h
  integer  :: hQC
  real(r8) :: lai
  real(r8) :: gpp
end type towerdata

type(towerdata) :: tower

!-----------------------------------------------------------------------
! start of executable code
!-----------------------------------------------------------------------

call initialize_utilities('ozflux_to_obs')

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "ozflux_to_obs_nml", iunit)
read(iunit, nml = ozflux_to_obs_nml, iostat = rcio)
call check_namelist_read(iunit, rcio, "ozflux_to_obs_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=ozflux_to_obs_nml)
if (do_nml_term()) write(     *     , nml=ozflux_to_obs_nml)

! time setup
call set_calendar_type(GREGORIAN)
offset    = set_time(nint(abs(timezoneoffset)*3600.0_r8),0)
prev_time = set_time(0, 0)

if (verbose) print *, 'tower located at lat, lon, elev  =', latitude, longitude, elevation
if (verbose) print *, 'flux observations taken at       =', flux_height,'m'

! check the lat/lon values to see if they are ok
if (longitude < 0.0_r8) longitude = longitude + 360.0_r8

if (( latitude > 90.0_r8 .or. latitude  <  -90.0_r8 ) .or. &
    (longitude <  0.0_r8 .or. longitude >  360.0_r8 )) then

   write (string2,*)'latitude  should be [-90, 90] but is ',latitude
   write (string3,*)'longitude should be [  0,360] but is ',longitude

   string1 ='tower location error in input.nml&ozflux_to_obs_nml'
   call error_handler(E_ERR,'ozflux_to_obs', string1, source, revision, &
                      revdate, text2=string2,text3=string3)

endif

! We need to know the maximum number of observations in the input file.
! Each line has info for the 3 observations we want.
! The max possible number of obs needs to be specified but it
! will only write out the actual number created.
! Each observation in this series will have a single
! observation value and a quality control flag.  
! Initialize two empty observations - one to track location
! in observation sequence - the other is for the new observation.

iunit = open_file(text_input_file, 'formatted', 'read')
if (verbose) print *, 'opened input file ' // trim(text_input_file)

nlines     = count_file_lines(iunit)
max_obs    = 3*nlines
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
call set_qc_meta_data(  obs_seq, 1, 'original QC')

! The first line describes all the fields ... column headers, if you will

rewind(iunit)
call decode_header(iunit, nwords)

obsloop: do iline = 2,nlines

   ! read in entire text line into a buffer
   read(iunit,'(A)',iostat=rcio) bigline
   if (rcio < 0) exit obsloop
   if (rcio > 0) then
      write (string1,'(''Cannot read (error '',i3,'') line '',i8,'' in '',A)') &
                    rcio, iline, trim(text_input_file)
      call error_handler(E_ERR,'main', string1, source, revision, revdate)
   endif

   input_line = adjustl(bigline)

   ! parse the line into the tower structure (including the observation time)
   call stringparse(input_line, nwords, iline, ',')

   if (iline <= 83) then
      write(*,*)''
      write(*,*)'Check of the first observation: (column,string,value)'
      write(*,*)tower%hindex    , tower%hstring     , tower%h
      write(*,*)tower%hQCindex  , tower%hQCstring   , tower%hQC
      write(*,*)tower%leindex   , tower%lestring    , tower%le
      write(*,*)tower%leQCindex , tower%leQCstring  , tower%leQC
      write(*,*)tower%neeindex  , tower%neestring   , tower%nee
      write(*,*)tower%neeQCindex, tower%neeQCstring , tower%neeQC
      write(*,*)tower%laiindex  ,                 0 , tower%lai
      write(*,*)tower%gppindex  ,                 0 , tower%gpp
      call print_date(tower%time_obs, 'observation date is')
      call print_time(tower%time_obs, 'observation time is')
   end if

   if (verbose) call print_date(tower%time_obs, 'obs time is')

   call get_time(tower%time_obs, osec, oday)

   ! make an obs derived type, and then add it to the sequence
   ! If the QC value is good, use the observation.
   ! Increasingly larger QC values are more questionable quality data.
   ! The observation errors are from page 183, Table 7.1(A) in 
   ! Chapter 7 of a book by A.D. Richardson et al. via Andy Fox.  

   ! Sensible Heat Flux [W m-2]
   if ((tower%hQC <= maxgoodqc) .and. (tower%h /= MISSING_R8)) then
      oerr = 10.0_r8 + abs(tower%h)*0.22_r8
      qc   = real(tower%hQC,r8)
      call create_3d_obs(latitude, longitude, flux_height, VERTISHEIGHT, tower%h, &
                         TOWER_LATENT_HEAT_FLUX, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, tower%time_obs, prev_obs, prev_time, first_obs)
   endif

   ! Latent Heat Flux [W m-2]
   if ((tower%leQC <= maxgoodqc) .and. (tower%le /= MISSING_R8)) then
      oerr = 10.0_r8 + abs(tower%le)*0.32_r8
      qc   = real(tower%leQC,r8)
      call create_3d_obs(latitude, longitude, flux_height, VERTISHEIGHT, tower%le, &
                         TOWER_SENSIBLE_HEAT_FLUX, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, tower%time_obs, prev_obs, prev_time, first_obs)
   endif

   ! Net Ecosystem Exchange  [umol m-2 s-1]
   if ((tower%neeQC <= maxgoodqc) .and. (tower%nee /= MISSING_R8)) then
      if (tower%nee <= 0) then
         oerr = (2.0_r8 + abs(tower%nee)*0.1_r8) * umol_to_gC
      else
         oerr = (2.0_r8 + abs(tower%nee)*0.4_r8) * umol_to_gC
      endif
      tower%nee = -tower%nee * umol_to_gC ! to match convention in CLM [gC m-2 s-1]
      qc        = real(tower%neeQC,r8)
      call create_3d_obs(latitude, longitude, flux_height, VERTISHEIGHT, tower%nee, &
                         TOWER_NETC_ECO_EXCHANGE, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, tower%time_obs, prev_obs, prev_time, first_obs)
   endif

   ! Leaf Area Index   FIXME ... errors, qcs etc.
   if (tower%lai /= MISSING_R8) then
      oerr = 0.3_r8   ! FIXME
      qc   = 0                 ! FIXME
      call create_3d_obs(latitude, longitude, flux_height, VERTISHEIGHT, tower%lai, &
                         LEAF_AREA_INDEX, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, tower%time_obs, prev_obs, prev_time, first_obs)
   endif

!  FIXME ... add GPP at at later date.
!  ! Gross Primary Productivity  g Carbon per M^2/s
!  if (tower%gpp <= maxgoodqc) then
!     oerr = 10.0_r8 + abs(tower%gpp)*0.32_r8
!     qc   = real(tower%lQC,r8)
!     call create_3d_obs(latitude, longitude, flux_height, VERTISHEIGHT, tower%gpp, &
!                        GROSS_PRIMARY_PRODUCTIVITY, oerr, oday, osec, qc, obs)
!     call add_obs_to_seq(obs_seq, obs, tower%time_obs, prev_obs, prev_time, first_obs)
!  endif

end do obsloop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (verbose) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
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
      call error_handler(E_ERR,'count_file_lines', string1, &
                         source, revision, revdate)
   endif
   count_file_lines = count_file_lines + 1

enddo countloop
rewind(iunit)

if (count_file_lines >= tenmillion) then
   write (string1,'('' suspiciously large number of lines '',i8)')count_file_lines
   call error_handler(E_MSG,'count_file_lines', string1, &
                         source, revision, revdate)
endif

end function count_file_lines




subroutine decode_header(iunit,ncolumns)
! Reads the first line of the header and parses the information.
! And by parse, I mean determine which columns are the columns
! of interest.

integer, intent(in) :: iunit
integer, intent(out) :: ncolumns

integer, parameter :: maxwordlength = 40
integer :: i,charcount,columncount,wordlength,maxlength
character(len=maxwordlength), dimension(:), allocatable :: columns
integer, dimension(7) :: qc = 0

! Read the line and strip off any leading whitespace.

read(iunit,'(A)',iostat=rcio) bigline
if (rcio /= 0) then
  write(string1,*)'Cannot parse header. Begins <',trim(bigline(1:40)),'>'
  call error_handler(E_ERR,'decode_header',string1, source, revision, revdate)
endif

input_line = adjustl(bigline)

! Count how many commas are in the line - use this to determine how many columns

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
         call error_handler(E_ERR,'decode_header',string1, source, revision, revdate)
      endif
      columns(columncount) = input_line((i-wordlength):(i-1)) 
      if (verbose) write(*,*)'word(',columncount,') is ',columns(columncount)
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
    call error_handler(E_ERR,'decode_header',string1,source,revision,revdate, &
                       text2=trim(string2), text3=trim(input_line))
endif

columns(ncolumns) = input_line((charcount+1):len_trim(input_line))

if (verbose) write(*,*)'(last) word(',ncolumns,') is ',columns(ncolumns)

! get to the task at hand

tower%dateindex  = Match(columns, tower%datestring)
tower%hindex     = Match(columns, tower%hstring)     ! Sensible Heat
tower%hQCindex   = Match(columns, tower%hQCstring)
tower%leindex    = Match(columns, tower%lestring)    ! Latent Heat
tower%leQCindex  = Match(columns, tower%leQCstring)
tower%neeindex   = Match(columns, tower%neestring)   ! Carbon Flux
tower%neeQCindex = Match(columns, tower%neeQCstring)
tower%laiindex   = Match(columns, tower%laistring)
tower%gppindex   = Match(columns, tower%gppstring)

write(*,110)tower%datestring , tower%dateindex
write(*,110)tower%neestring  , tower%neeindex
write(*,110)tower%neeQCstring, tower%neeQCindex
write(*,110)tower%lestring   , tower%leindex
write(*,110)tower%leQCstring , tower%leQCindex
write(*,110)tower%hstring    , tower%hindex
write(*,110)tower%hQCstring  , tower%HQCindex
write(*,110)tower%laistring  , tower%laiindex
write(*,110)tower%gppstring  , tower%gppindex
110 format('index for ',A15,' is ',i4)

! Check to make sure we got all the indices we need

qc( 1) = CheckIndex( tower%dateindex  , tower%datestring )
qc( 2) = CheckIndex( tower%hindex     , tower%hstring    )
qc( 3) = CheckIndex( tower%hQCindex   , tower%hQCstring  )
qc( 4) = CheckIndex( tower%leindex    , tower%lestring   )
qc( 5) = CheckIndex( tower%leQCindex  , tower%leQCstring )
qc( 6) = CheckIndex( tower%neeindex   , tower%neestring  )
qc( 7) = CheckIndex( tower%neeQCindex , tower%neeQCstring)

if (any(qc /= 0)) then
  write(string1,*)'Did not find all the required column indices.'
  call error_handler(E_ERR,'decode_header',string1, source, revision, revdate)
endif

deallocate(columns)

end subroutine decode_header



function CountChar(str1,separator)
! Count the number of instances of the single character in a character string.
! useful when parsing a comma-separated list, for example.
! Count the commas and add 1 to get the number of items in the list.

integer                      :: CountChar
character(len=*), intent(in) :: str1
character,        intent(in) :: separator

integer :: i

CountChar = 0
do i = 1,len_trim(str1)
   if (str1(i:i) == separator) then
      CountChar = CountChar + 1
   endif
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
! I want to check ALL the indexes before fatally ending.

integer                       :: CheckIndex
integer,          intent(in)  :: myindex
character(len=*), intent(in)  :: context

if (myindex == 0) then
   write(string1,*)'Did not find column header matching ',trim(context)
   call error_handler(E_MSG,'decode_header:CheckIndex',string1, source, revision, revdate)
   CheckIndex = -1 ! not a good thing
else
   CheckIndex = 0  ! Good to go
endif

end function CheckIndex



subroutine stringparse(str1, nwords, linenum, separator)
! The first 'word' is a complicated string:
! just declare everything as reals and chunk it

character(len=*), intent(in) :: str1
integer         , intent(in) :: nwords
integer         , intent(in) :: linenum
character,        intent(in) :: separator

real(r8), allocatable, dimension(:) :: values
integer :: i,CountChar, iword, istart,ifinish
integer :: iyear, imonth, iday, ihour, imin, isec, seconds
type(time_type) :: time0, time1, time2

character(len=40) :: numberstring
character(len=20) :: datestring
allocate(values(nwords-1))

! initialize all values
values(:) = MISSING_R8

! The first word is a date ... YYYY-MM-DD hh:mm:ss
CountChar = 0
FirstWord : do i = 1,len_trim(str1)
   if (str1(i:i) == separator)  then
      CountChar = i
      exit FirstWord
   endif
enddo FirstWord

read(str1(1:CountChar),900,iostat=rcio) iyear,imonth,iday,ihour,imin,isec
900 format(i4,5(1x,i2),1x)
if (rcio /= 0) then
   write(string1,*)'Cannot parse date at line ',linenum
   write(string2,*)'String is <',trim(str1(1:CountChar)),'>'
   call error_handler(E_ERR,'stringparse',string1, source, revision, revdate, &
   text2=trim(string2))
endif

! TJH DEBUG write(*,*)'date read as', iyear,imonth,iday,ihour,imin,isec

! Parse the remainder of the line 
iword = 1
istart = CountChar+1 ! where the next word starts

AllWords: do i = istart,len_trim(str1)

   if (str1(i:i) == separator) then  ! we have a word to process
      iword        = iword + 1

      if (iword == 177) then ! embedded is another 'YYYY-MM-DD HH:mm:ss'
         istart = i
         cycle AllWords
      endif

      if ((i-istart) > 1) then
!        write(*,*)'word ',iword,' is ', istart+1, i-1
         numberstring = str1((istart+1):(i-1))
!        write(*,*)'word ',iword,' is ', istart, i, '<',trim(numberstring),'>'
         read(numberstring,*,iostat=rcio) values(iword)
         if (rcio /= 0) then
            write(string1,*)'Cannot parse values at line',linenum, 'word', iword
            write(string2,*)'word is <',trim(numberstring),'>'
            write(string3,*)'word boundaries are columns',istart+1,i-1
            call error_handler(E_ERR, 'stringparse', string1, &
                 source, revision, revdate, text2=trim(string2), text3=string3)
         endif
      endif
      istart = i
   endif

enddo AllWords

! Stuff what we want into the tower structure
!
! Convert to 'CLM-friendly' units AFTER we determine observation error variance.
! That happens in the main routine.
!
! NEE_or_fMDS    has units     [umolCO2 m-2 s-1] 
! H_f            has units     [W m-2]
! LE_f           has units     [W m-2]
!
! (CLM) NEE      has units     [gC m-2 s-1]

tower%nee   = values(tower%neeindex)
tower%le    = values(tower%leindex)
tower%h     = values(tower%hindex)
tower%lai   = values(tower%laiindex)
tower%gpp   = values(tower%gppindex)

tower%neeQC = nint(values(tower%neeQCindex))
tower%leQC  = nint(values(tower%leQCindex))
tower%hQC   = nint(values(tower%hQCindex))

deallocate(values)

tower%time_obs = set_date(iyear, imonth, iday, ihour, imin, isec)

! 8AM East Coast is 1PM Greenwich 
if (timezoneoffset < 0.0_r8) then
   tower%time_obs = tower%time_obs + offset
else
   tower%time_obs = tower%time_obs - offset
endif

! The QC values can be 'missing' ... in which case the values are too

if (tower%neeQC < 0) tower%neeQC = maxgoodqc + 1000 
if (tower%leQC  < 0) tower%leQC  = maxgoodqc + 1000
if (tower%hQC   < 0) tower%hQC   = maxgoodqc + 1000

! if (tower%neeQC < maxgoodqc) then
!    write(*,*)'nee umol m-2 s-1 ',tower%nee
!    write(*,*)'nee   gC m-2 s-1 ',tower%nee*umol_to_gC
! endif
 
end subroutine stringparse



end program ozflux_to_obs


! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
