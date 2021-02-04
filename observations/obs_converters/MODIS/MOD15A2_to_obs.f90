! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program MOD15A2_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   MOD15A2_to_obs - 
!
!     created 21 April 2014   Tim Hoar NCAR/IMAGe
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, MISSING_R8

use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              register_module, error_handler, E_MSG, E_ERR, &
                              open_file, close_file, do_nml_file, do_nml_term, &
                              check_namelist_read, find_namelist_in_file, &
                              nmlfileunit, logfileunit, find_textfile_dims

use  time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                              set_date, set_time, get_time, print_time, &
                              print_date, operator(-), operator(+), &
                              operator(>), operator(<), &
                              operator(==), operator(/=), operator(<=), operator(>=)

use      location_mod, only : VERTISSURFACE

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data

use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq

use      obs_kind_mod, only : MODIS_LEAF_AREA_INDEX, &
                              MODIS_FPAR

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! Namelist with default values
!-----------------------------------------------------------------------

character(len=128) :: text_input_file = 'MOD15A2.fn_usbouldr.txt'
character(len=128) :: metadata_file   = 'MOD15A2_site_metadata.txt'
character(len=128) :: obs_out_file    = 'obs_seq.out'
real(r8)           :: maxgoodqc       = 10
logical            :: verbose         = .false.

namelist /MOD15A2_to_obs_nml/ text_input_file, metadata_file, obs_out_file, maxgoodqc, verbose

!-----------------------------------------------------------------------
! globally-scoped variables
!-----------------------------------------------------------------------

character(len=600)      :: input_line, bigline
character(len=256)      :: string1, string2, string3
integer                 :: iline
logical                 :: first_obs
integer                 :: rcio, iunit
integer                 :: num_copies, num_qc, max_obs
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: prev_time

type modis_type
  character(len=20) :: HDFnamestring       = 'HDFname'
  character(len=20) :: Productstring       = 'Product'
  character(len=20) :: Datestring          = 'Date'
  character(len=20) :: Sitestring          = 'Site'
  character(len=20) :: ProcessDatestring   = 'ProcessDate'
  character(len=20) :: Bandstring          = 'Band'
  character(len=20) :: FparExtra_QCstring  = 'FparExtra_QC'
  character(len=20) :: FparLai_QCstring    = 'FparLai_QC'
  character(len=20) :: FparStdDevstring    = 'FparStdDev_1km'
  character(len=20) :: Fparstring          = 'Fpar_1km'
  character(len=20) :: LaiStdDevstring     = 'LaiStdDev_1km'
  character(len=20) :: Laistring           = 'Lai_1km'
  integer  :: HDFnameindex
  integer  :: Productindex
  integer  :: Dateindex
  integer  :: Siteindex
  integer  :: ProcessDateindex
  integer  :: Bandindex
  integer  :: FparExtra_QCindex
  integer  :: FparLai_QCindex
  integer  :: FparStdDevindex
  integer  :: Fparindex
  integer  :: LaiStdDev_index
  integer  :: Lai_index
  character(len=80) :: HDFname
  character(len=80) :: Product
  character(len=80) :: Date
  character(len=80) :: Site
  character(len=80) :: ProcessDate
  integer  :: FparExtra_QC
  integer  :: FparLai_QC
  real(r8) :: FparStdDev
  real(r8) :: Fpar ! the fraction of incident photosynthetically active 
                   ! radiation (400-700nm) absorbed by the green elements 
                   ! of a vegetation canopy.
  real(r8) :: LaiStdDev
  real(r8) :: Lai  ! one-sided green leaf area per unit ground area in broadland canopies
                   ! and one-half the total needle surface area per ground area 
                   ! in coniferous  canopies.
  type(time_type) :: time_obs
  type(time_type) :: FparExtra_QC_time
  type(time_type) :: FparLai_QC_time
  type(time_type) :: FparStdDev_time
  type(time_type) :: Fpar_time
  type(time_type) :: LaiStdDev_time
  type(time_type) :: Lai_time
end type modis_type
    type(modis_type) :: modis

! The preselected "Field Site and Flux Tower" data files have a
! preliminary record describing the file format.
! the first line/record  contains:
! HDFname,Product,Date,Site,ProcessDate,Band,1,2,3, ...,47,48,49

type record_type
  character(len=80) :: HDFname
  character(len=80) :: Product
  character(len=80) :: Date
  character(len=80) :: Site
  character(len=80) :: ProcessDate
  character(len=20) :: Band
  integer, dimension(49) :: ivalues
  type(time_type) :: time_obs
end type record_type
    type(record_type) :: modisrecord

! The preselected "Field Site and Flux Tower" data files have a
! companion file with the metadata for the site.
type metadata_type
   integer :: site_name_indx      ! Match(columns,'SITE NAME')
   integer :: country_indx        ! Match(columns,'COUNTRY')
   integer :: latitude_indx       ! Match(columns,'LATITUDE')
   integer :: longitude_indx      ! Match(columns,'LONGITUDE')
   integer :: modis_site_id_indx  ! Match(columns,'MODIS SITE ID')
   integer :: fluxnet_id_indx     ! Match(columns,'FLUXNET ID')
   integer :: fluxnet_key_id_indx ! Match(columns,'FLUXNET KEY ID')
   integer :: nstations = 0
   character(len=80), allocatable, dimension(:) :: site_name
   character(len=80), allocatable, dimension(:) :: country
   real(r8),          allocatable, dimension(:) :: latitude
   real(r8),          allocatable, dimension(:) :: longitude
   character(len=80), allocatable, dimension(:) :: modis_site_id
   character(len=80), allocatable, dimension(:) :: fluxnet_id
   character(len=80), allocatable, dimension(:) :: fluxnet_key_id
end type metadata_type
    type(metadata_type) :: metadata

integer  :: siteindex
real(r8) :: latitude
real(r8) :: longitude
integer  :: landcover

!-----------------------------------------------------------------------
! start of executable code
!-----------------------------------------------------------------------

call initialize_utilities('MOD15A2_to_obs')

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "MOD15A2_to_obs_nml", iunit)
read(iunit, nml = MOD15A2_to_obs_nml, iostat = rcio)
call check_namelist_read(iunit, rcio, "MOD15A2_to_obs_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=MOD15A2_to_obs_nml)
if (do_nml_term()) write(     *     , nml=MOD15A2_to_obs_nml)

! time setup
call set_calendar_type(GREGORIAN)
prev_time = set_time(0, 0)

! Check to make sure the metadata filename exists.

call Get_MODIS_subset_metadata()

! We need to know the maximum number of observations in the input file.
! There are multiple lines for each observation, so file line count
! is more than enough. This program will only write out the actual number 
! created.

call find_textfile_dims(text_input_file, max_obs)

iunit = open_file(text_input_file, 'formatted', 'read')
if (verbose) print *, 'opened input file ' // trim(text_input_file)

call decode_data_header(iunit) ! which data csv columns are which

num_copies = 1
num_qc     = 1
first_obs  = .true.

call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)
call init_record(modisrecord) ! one line of data from the csv file

! the first one needs to contain the string 'observation' and the
! second needs the string 'QC'.
call set_copy_meta_data(obs_seq, 1, 'observation')
call set_qc_meta_data(  obs_seq, 1, 'MODIS QC')

obsloop: do iline = 2,max_obs

   ! read in entire text line into a buffer
   read(iunit,'(A)',iostat=rcio) bigline
   if (rcio < 0) exit obsloop
   if (rcio > 0) then
      write (string1,'(''Cannot read (error '',i3,'') line '',i8,'' in '',A)') &
                    rcio, iline, trim(text_input_file)
      call error_handler(E_ERR,'main', string1, source, revision, revdate)
   endif

   input_line = adjustl(bigline)

   ! parse the line into temporary modis structure (including the observation time)
   ! If the time is an existing one, just add contents to modis structure and
   ! continue. If time is a new time, flush the buffer, reset and continue.
   call stringparse( input_line, iline )

   if (iline <= 2) then

      ! Only needs to be done once since file is for one site.
      ! Extract the site from the filename and determine metadata.

      siteindex = get_metadata(modisrecord%Site)
      latitude  = metadata%latitude( siteindex)
      longitude = metadata%longitude(siteindex)
      landcover = -23 ! FIXME ... no landcover metadata included in this version.

      ! check the lat/lon values to see if they are ok
      if (longitude < 0.0_r8) longitude = longitude + 360.0_r8

      if (( latitude > 90.0_r8 .or. latitude  <  -90.0_r8 ) .or. &
          (longitude <  0.0_r8 .or. longitude >  360.0_r8 )) then

         write (string2,*)'latitude  should be [-90, 90] but is ',latitude
         write (string3,*)'longitude should be [  0,360] but is ',longitude

         string1 ='modis location error in input.nml&MOD15A2_to_obs_nml'
         call error_handler(E_ERR,'MOD15A2_to_obs', string1, source, revision, &
                      revdate, text2=string2,text3=string3)
      endif

      ! First observation can define the following
      modis%HDFname     = modisrecord%HDFname
      modis%Product     = modisrecord%Product
      modis%Date        = modisrecord%Date
      modis%Site        = modisrecord%Site
      modis%ProcessDate = modisrecord%ProcessDate
      modis%time_obs    = modisrecord%time_obs

      if (verbose) then
         write(*,*)
         write(*,*)'Check of the first observation: (column,string,value)'
         write(*,*)modis%HDFnameindex,     modis%HDFnamestring,     trim(modis%HDFname)
         write(*,*)modis%Productindex,     modis%Productstring,     trim(modis%Product)
         write(*,*)modis%Dateindex,        modis%Datestring,        trim(modis%Date)
         write(*,*)modis%Siteindex,        modis%Sitestring,        trim(modis%Site)
         write(*,*)modis%ProcessDateindex, modis%ProcessDatestring, trim(modis%ProcessDate)
         write(*,*)modis%Bandindex,        modis%Bandstring,        trim(modisrecord%Band)
         write(*,*)'   site # ',siteindex,' is located at lat, lon =', latitude, longitude
         write(*,*)'   data(25) is ',modisrecord%ivalues(25)
      endif
   endif

                call print_date(modisrecord%time_obs, modisrecord%Band//' obs date is')
   if (verbose) call print_time(modisrecord%time_obs, modisrecord%Band//' obs time is')
   if (verbose) write(*,*)modisrecord%Band,' data(25) is ',modisrecord%ivalues(25)
   if (verbose) write(*,*)

   if (modisrecord%time_obs /= modis%time_obs) then
      ! time to write the current modis structure to the obs sequence file
      call stow_observation()
      modis%time_obs = modisrecord%time_obs
   endif

   ! determine which values were read and update 

   select case ( modisrecord%Band )
   case ( 'FparExtra_QC' )
     ! This is of no value to us at this time.
     modis%FparExtra_QC      = modisrecord%ivalues(25) ! 25 is the center pixel of the 49 in the record
     modis%FparExtra_QC_time = modisrecord%time_obs
   case ( 'FparLai_QC' )
     modis%FparLai_QC        = qc_and_scale(modisrecord)
     modis%FparLai_QC_time   = modisrecord%time_obs
   case ( 'FparStdDev_1km' )
     modis%FparStdDev        = qc_and_scale(modisrecord)
     modis%FparStdDev_time   = modisrecord%time_obs
   case ( 'Fpar_1km' )
     modis%Fpar              = qc_and_scale(modisrecord)
     modis%Fpar_time         = modisrecord%time_obs
   case ( 'LaiStdDev_1km' )
     modis%LaiStdDev         = qc_and_scale(modisrecord)
     modis%LaiStdDev_time    = modisrecord%time_obs
   case ( 'Lai_1km' )
     modis%Lai               = qc_and_scale(modisrecord)
     modis%Lai_time          = modisrecord%time_obs
   case default
   end select

end do obsloop

call stow_observation() ! Flush that last observation 

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (verbose) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
endif

! end of main program
call finalize_utilities()

contains


subroutine Get_MODIS_subset_metadata()
! Reads the list of station names and locations

character(len=80), dimension(7) :: columns

integer :: i, nlines, iline, nquotes, nwords, s1, s2
integer, allocatable, dimension(:,:) :: positions

! The station metadata file has two header records. The first is important,
! the second is fundamentally a 'blank' line.
!SITE_NAME,COUNTRY,LATITUDE,LONGITUDE,MODIS_SITE_ID,FLUXNET_ID,FLUXNET_KEY_ID
! ,,,,,,
!Tomakomai_National_Forest,Japan,42.73697,141.51864,fn_jptomako,JP-Tom,jp.tomakomai.01

call find_textfile_dims(metadata_file, nlines)
metadata%nstations = nlines - 2

if (verbose) print *, 'there are ',metadata%nstations,' lines in the file.'

allocate(metadata%site_name(     metadata%nstations))
allocate(metadata%country(       metadata%nstations))
allocate(metadata%latitude(      metadata%nstations))
allocate(metadata%longitude(     metadata%nstations))
allocate(metadata%modis_site_id( metadata%nstations))
allocate(metadata%fluxnet_id(    metadata%nstations))
allocate(metadata%fluxnet_key_id(metadata%nstations))

iunit = open_file(metadata_file, 'formatted', 'read')
if (verbose) print *, 'opened input file ' // trim(metadata_file)
call decode_metadata(iunit)

! Skip the first two lines 
rewind(iunit)
read(iunit,'(A)',iostat=rcio) bigline
read(iunit,'(A)',iostat=rcio) bigline

readloop: do iline = 1,metadata%nstations

   ! read in entire text line into a buffer
   read(iunit,'(A)',iostat=rcio) bigline
   if (rcio < 0) exit readloop
   if (rcio > 0) then
      write (string1,'(''Cannot read (error '',i3,'') line '',i8,'' in '',A)') &
                    rcio, iline+2, trim(text_input_file)
      call error_handler(E_ERR,'Get_MODIS_subset_metadata', string1, source, revision, revdate)
   endif

   input_line = adjustl(bigline)

   ! Some of the site names are quoted strings with embedded commas that screw up the
   ! rest of the parsing logic. I am going to search for matching (double) quotes and
   ! replace any embedded commas with whitespace. bother ...

   nquotes = CountChar(input_line,'"')

   if ( nquotes == 2) then
      allocate(positions(1,2))
      positions(1,1) = scan(input_line,'"',back=.FALSE.)
      positions(1,2) = scan(input_line,'"',back=.TRUE.)
      do i = positions(1,1),positions(1,2)
         if (input_line(i:i) == ',') input_line(i:i) = ' '
      enddo
      deallocate(positions)
   elseif ( nquotes == 0 ) then
      ! normal situation, do nothing
   else
      write(string1,*)'unsupported number of double quotes on line ',iline+2
      write(string2,*)'in file ',trim(text_input_file)
      write(string3,*)'found ',nquotes,' double quotes. expected a pair, if any.'
      call error_handler(E_ERR,'Get_MODIS_subset_metadata', string1, &
                 source, revision, revdate, text2=string2, text3=string3)
   endif

   ! Find the commas and the corresponding word lengths.

   nwords = CountChar(input_line,',')
   nwords = nwords + 1           ! count the last word, too

   allocate(positions(nwords,2)) ! each word has a start and stop index in string
   call SetWordPositions(input_line, ',', positions)

   s1 = positions(metadata%site_name_indx,1);
   s2 = positions(metadata%site_name_indx,2);
   metadata%site_name(iline)      = trim(input_line(s1:s2))

   s1 = positions(metadata%country_indx,1);
   s2 = positions(metadata%country_indx,2);
   metadata%country(iline)        = trim(input_line(s1:s2))

   s1 = positions(metadata%modis_site_id_indx,1);
   s2 = positions(metadata%modis_site_id_indx,2);
   metadata%modis_site_id(iline)  = trim(input_line(s1:s2))

   s1 = positions(metadata%fluxnet_id_indx,1);
   s2 = positions(metadata%fluxnet_id_indx,2);
   metadata%fluxnet_id(iline)     = trim(input_line(s1:s2))

   s1 = positions(metadata%fluxnet_key_id_indx,1);
   s2 = positions(metadata%fluxnet_key_id_indx,2);
   metadata%fluxnet_key_id(iline) = trim(input_line(s1:s2))
      
   s1 = positions(metadata%latitude_indx,1);
   s2 = positions(metadata%latitude_indx,2);
   read(input_line(s1:s2),*,iostat=rcio) metadata%latitude(iline)
   if (rcio/=0) metadata%latitude(iline) = MISSING_R8

   s1 = positions(metadata%longitude_indx,1);
   s2 = positions(metadata%longitude_indx,2);
   read(input_line(s1:s2),*,iostat=rcio) metadata%longitude(iline)
   if (rcio/=0) metadata%longitude(iline) = MISSING_R8

   if (verbose) then
      write(*,*)'line number ',iline
      write(*,*)'   site_name      is ', trim(metadata%site_name(iline))
      write(*,*)'   country        is ', trim(metadata%country(iline))
      write(*,*)'   latitude       is ', metadata%latitude( iline)
      write(*,*)'   longitude      is ', metadata%longitude(iline)
      write(*,*)'   modis_site_id  is ', trim(metadata%modis_site_id(iline))
      write(*,*)'   fluxnet_id     is ', trim(metadata%fluxnet_id(iline))
      write(*,*)'   fluxnet_key_id is ', trim(metadata%fluxnet_key_id(iline))
      write(*,*)
   endif
   deallocate(positions)

enddo readloop

end subroutine Get_MODIS_subset_metadata


subroutine decode_metadata(iunit)
! Reads the first line of the header and parses the information.
! And by parse, I mean determine which columns are the columns
! of interest.

integer,              intent(in)  :: iunit

integer, parameter :: maxwordlength = 80
integer :: i,charcount,columncount,wordlength
character(len=maxwordlength), dimension(:), allocatable :: columns
integer, dimension(7) :: qc = 0

! Read the line and strip off any leading whitespace.

read(iunit,'(A)',iostat=rcio) bigline
if (rcio /= 0) then
  write(string1,*)'Cannot parse header. Begins <',trim(bigline(1:40)),'>'
  call error_handler(E_ERR,'decode_metadata',string1, source, revision, revdate)
endif

input_line = adjustl(bigline)

! Count how many commas are in the line - use this to determine how many columns

charcount = CountChar(input_line,',')
columncount = charcount + 1
allocate(columns(columncount))

columncount  = 0  ! track the number of columns
wordlength   = 0  ! number of characters in the column descriptor
charcount    = 0  ! the position of the (last) comma
do i = 1,len_trim(input_line)
   if (input_line(i:i) == ',') then
      columncount = columncount + 1
      if (wordlength > maxwordlength) then
         write(string1,*)'unexpected long word ... starts <',&
                           input_line((i-wordlength):(i-1)),'>'
         call error_handler(E_ERR,'decode_metadata',string1, source, revision, revdate)
      endif
      columns(columncount) = input_line((i-wordlength):(i-1)) 
      if (verbose) write(*,*)'word(',columncount,') is ',trim(columns(columncount))
      wordlength = 0
      charcount = i
   else
      wordlength = wordlength + 1
   endif
enddo

! There is one more column after the last comma

columns(columncount+1) = input_line((charcount+1):len_trim(input_line))

metadata%site_name_indx      = Match(columns, 'SITE NAME' )
metadata%country_indx        = Match(columns, 'COUNTRY' )
metadata%latitude_indx       = Match(columns, 'LATITUDE' )
metadata%longitude_indx      = Match(columns, 'LONGITUDE' )
metadata%modis_site_id_indx  = Match(columns, 'MODIS SITE ID' )
metadata%fluxnet_id_indx     = Match(columns, 'FLUXNET ID' )
metadata%fluxnet_key_id_indx = Match(columns, 'FLUXNET KEY ID' )

! Check to make sure we got all the indices we need

qc( 1) = CheckIndex( metadata%site_name_indx,     'SITE NAME' )
qc( 2) = CheckIndex( metadata%country_indx,       'COUNTRY' )
qc( 3) = CheckIndex( metadata%latitude_indx,      'LATITUDE' )
qc( 4) = CheckIndex( metadata%longitude_indx,     'LONGITUDE' )
qc( 5) = CheckIndex( metadata%modis_site_id_indx, 'MODIS SITE ID' )
qc( 6) = CheckIndex( metadata%fluxnet_id_indx,    'FLUXNET ID' )
qc( 7) = CheckIndex( metadata%fluxnet_key_id_indx,'FLUXNET KEY ID' )

if (any(qc /= 0) ) then
  write(string1,*)'Did not find all the required column indices.'
  call error_handler(E_ERR,'decode_metadata',string1, source, revision, revdate)
endif

deallocate(columns)

end subroutine decode_metadata



subroutine decode_data_header(iunit)
! Reads the first line of the header and parses the information.
! And by parse, I mean determine which columns are the columns
! of interest.
!
! The first line/record should be:
! HDFname,Product,Date,Site,ProcessDate,Band,1,2,3 ... 47,48,49
! We search for the 'Band' column and use that to decode whether
! the observation is FparExtra_QC or FparLai_QC or Lai_1km or LaiStdDev_1km or ...

integer, intent(in) :: iunit

integer, parameter :: maxwordlength = 80
integer :: i,charcount,columncount,wordlength
character(len=maxwordlength), dimension(:), allocatable :: columns
integer, dimension(6) :: qc = 0

! Read the line and strip off any leading whitespace.

read(iunit,'(A)',iostat=rcio) bigline
if (rcio /= 0) then
  write(string1,*)'Cannot parse header. Begins <',trim(bigline(1:40)),'>'
  call error_handler(E_ERR,'decode_data_header',string1, source, revision, revdate)
endif

input_line = adjustl(bigline)

! Count how many commas are in the line - use this to determine how many columns

charcount = CountChar(input_line,',')
columncount = charcount + 1
allocate(columns(columncount))

columncount  = 0  ! track the number of columns
wordlength   = 0  ! number of characters in the column descriptor
charcount    = 0  ! the position of the (last) comma
do i = 1,len_trim(input_line)
   if (input_line(i:i) == ',') then
      columncount = columncount + 1
      if (wordlength > maxwordlength) then
         write(string1,*)'unexpected long word ... starts <',&
                           input_line((i-wordlength):(i-1)),'>'
         call error_handler(E_ERR,'decode_data_header',string1, source, revision, revdate)
      endif
      columns(columncount) = input_line((i-wordlength):(i-1)) 
      if (verbose) write(*,*)'word(',columncount,') is ',trim(columns(columncount))
      wordlength = 0
      charcount = i
   else
      wordlength = wordlength + 1
   endif
enddo

! check to make sure there are 6 metadata items and 49 data items 
! in the data file metadata record.
if ( columncount /= 54 ) then
  write(string1,*)'Did not find all the required column separators.'
  write(string2,*)'Found ',columncount,'expected 54.'
  write(string3,*)'Require 6 to separate metadata and 48 to identify the 49 data records.'
  call error_handler(E_ERR,'decode_data_header',string1, &
             source, revision, revdate, text2=string2, text3=string3)
endif

! There is one more column after the last comma

columns(columncount+1) = input_line((charcount+1):len_trim(input_line))

! Finally get to the task at hand

modis%HDFnameindex     = Match( columns, modis%HDFnamestring )
modis%Productindex     = Match( columns, modis%Productstring )
modis%Dateindex        = Match( columns, modis%Datestring )
modis%Siteindex        = Match( columns, modis%Sitestring )
modis%ProcessDateindex = Match( columns, modis%ProcessDatestring )
modis%Bandindex        = Match( columns, modis%Bandstring )

! Check to make sure we got all the indices we need

qc( 1) = CheckIndex( modis%HDFnameindex,     trim(modis%HDFnamestring) )
qc( 2) = CheckIndex( modis%Productindex,     trim(modis%Productstring) )
qc( 3) = CheckIndex( modis%Dateindex,        trim(modis%Datestring) )
qc( 4) = CheckIndex( modis%Siteindex,        trim(modis%Sitestring) )
qc( 5) = CheckIndex( modis%ProcessDateindex, trim(modis%ProcessDatestring) )
qc( 6) = CheckIndex( modis%Bandindex,        trim(modis%Bandstring) )

if (any(qc /= 0) ) then
  write(string1,*)'Did not find all the required column indices.'
  call error_handler(E_ERR,'decode_data_header',string1, source, revision, revdate)
endif

deallocate(columns)

end subroutine decode_data_header



function   get_metadata(observation_site)
character(len=*), intent(in) :: observation_site
integer :: get_metadata

integer :: istation
character(len=40) :: x,y

get_metadata = -1

StationLoop: do istation = 1,metadata%nstations

   x = trim(metadata%modis_site_id(istation))
   y = trim(observation_site)

   if (x == y) then
      get_metadata = istation
      exit StationLoop
   endif

enddo StationLoop

end function get_metadata



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
! I want to check ALL the indexes before fatally ending.

integer                       :: CheckIndex
integer,          intent(in)  :: myindex
character(len=*), intent(in)  :: context

if (myindex == 0) then
   write(string1,*)'Did not find column header matching ',trim(context)
   call error_handler(E_MSG,'decode_data_header:CheckIndex',string1 )
   CheckIndex = -1 ! not a good thing
else
   CheckIndex = 0  ! Good to go
endif

end function CheckIndex



subroutine stringparse( myline, linenum )
! just chunk it

character(len=*), intent(in) :: myline
integer,          intent(in) :: linenum

character(len=80), dimension(6) :: column_names
integer :: iyear, idoy
type(time_type) :: time0, time1

modisrecord%ivalues(:) = MISSING_R8

read(myline,*,iostat=rcio) column_names, modisrecord%ivalues

if (rcio /= 0) then
   write(string1,*)'Cannot parse line',linenum,'. Begins <',trim(myline(1:40)),'>'
   call error_handler(E_ERR,'stringparse',string1, source, revision, revdate)
endif

modisrecord%HDFname     = trim(column_names(modis%HDFnameindex))
modisrecord%Product     = trim(column_names(modis%Productindex))
modisrecord%Date        = trim(column_names(modis%Dateindex))
modisrecord%Site        = trim(column_names(modis%Siteindex))
modisrecord%ProcessDate = trim(column_names(modis%ProcessDateindex))
modisrecord%Band        = trim(column_names(modis%Bandindex))

! These are 8 day composite values. I cannot determine if the Date is
! the center of the 8 day composite, the beginning, the end ... 
! is [DDD] modulo [0,365] or [1,366] ... FIXME

read(modisrecord%Date,'(1x,i4,i3)')iyear,idoy
time0   = set_date(iyear, 1, 1, 0, 0, 0)
time1   = set_time(0, idoy)
modisrecord%time_obs = time0 + time1

end subroutine stringparse





function qc_and_scale( myrecord )
! Checks the appropriate QC values for the value in question
! and converts 'good' values to proper units as referenced in:
! https://lpdaac.usgs.gov/sites/default/files/public/modis/docs/MODIS-LAI-FPAR-User-Guide.pdf

type(record_type), intent(in) :: myrecord
real(r8) :: qc_and_scale
real(r8) :: value

qc_and_scale = MISSING_R8

value = real(myrecord%ivalues(25),r8)

select case ( myrecord%Band )
case ('FparLai_QC')     ! Class flag
   qc_and_scale = value
case ('FparStdDev_1km') ! dimensionless
   if (value < 248) qc_and_scale = value * 0.01_r8
case ('Fpar_1km')       ! dimensionless
   if (value < 249) qc_and_scale = value * 0.01_r8
case ('LaiStdDev_1km')  ! dimensionless
   if (value < 248) qc_and_scale = value * 0.1_r8
case ('Lai_1km')        ! dimensionless
   if (value < 249) qc_and_scale = value * 0.1_r8
case default
   write(string1,*)'Did not match <',trim(myrecord%Band),'>'
   call error_handler(E_ERR,'qc_and_scale',string1, source, revision, revdate)
end select

end function qc_and_scale



subroutine init_record( myrecord )
type(record_type), intent(out) :: myrecord
myrecord%HDFname     = 'null'
myrecord%Product     = 'null'
myrecord%Date        = 'null'
myrecord%Site        = 'null'
myrecord%ProcessDate = 'null'
myrecord%Band        = 'null'
myrecord%ivalues(:)  = MISSING_R8
myrecord%time_obs    = set_time(0,0)
end subroutine init_record



subroutine stow_observation()
! so by now, everything should be in the modis structure

real(r8) :: oerr, qc, obval
integer  :: oday, osec
integer  :: bob, bits567

if (verbose) then
   call print_date(modis%time_obs)
   write(*,*) 'FparExtra_QC is ', modis%FparExtra_QC
   write(*,*) 'FparLai_QC   is ', modis%FparLai_QC
   write(*,*) 'FparStdDev   is ', modis%FparStdDev
   write(*,*) 'Fpar         is ', modis%Fpar
   write(*,*) 'LaiStdDev    is ', modis%LaiStdDev
   write(*,*) 'Lai          is ', modis%Lai
   write(*,*)
endif

! QC flags are binary-coded ascii strings e.g., 10011101
! bits 5,6,7 (the last three) are decoded as follows:
! 000 ... Main(RT) method used, best result possible (no saturation)
! 001 ... Main(RT) method used with saturation, Good, very usable
! 010 ... Main(RT) method failed due to bad geometry, empirical algorithm used
! 011 ... Main(RT) method failed due to other problems
! 100 ... pixel not produced at all

! relying on integer arithmetic
bob = 1000 * (modis%FparLai_QC / 1000)
bits567   =   modis%FparLai_QC - bob

! TJH-OK write(*,*)'modis%FparLai_QC and last 3',modis%FparLai_QC, bits567

if (bits567 > maxgoodqc) then
   modis%Fpar = MISSING_R8
   modis%Lai  = MISSING_R8
endif

call get_time(modis%time_obs, osec, oday)

if ( (modis%Fpar /= MISSING_R8) .and. (modis%FparStdDev /= MISSING_R8) ) then
   ! Check to make sure all Fpar observations are for the same time
   if ( (modis%Fpar_time /= modis%FparLai_QC_time)   .or. &
        (modis%Fpar_time /= modis%FparStdDev_time) ) then
      call print_time(modis%Fpar_time,        str='Fpar_time')
      call print_time(modis%FparLai_QC_time,  str='FparLai_QC_time')
      call print_time(modis%FparStdDev_time,  str='FparStdDev_time')
      write(string1,*)'All Fpar times did not match.'
      call error_handler(E_ERR,'stow_observation',string1, source, revision, revdate)
   endif
   qc    = real(bits567,r8)
   obval = modis%Fpar
   oerr  = min(0.1_r8, (modis%FparStdDev*modis%Fpar))

   call create_3d_obs(latitude, longitude, 0.0_r8, VERTISSURFACE, obval, &
               MODIS_FPAR, oerr, oday, osec, qc, obs)
   call add_obs_to_seq(obs_seq, obs, modis%time_obs, prev_obs, prev_time, first_obs)
endif

if ( (modis%Lai /= MISSING_R8) .and. (modis%LaiStdDev /= MISSING_R8) ) then
   ! Check to make sure all Lai observations are for the same time
   if ( (modis%Lai_time /= modis%FparLai_QC_time) .or. &
        (modis%Lai_time /= modis%LaiStdDev_time) ) then
      call print_time(modis%Lai_time,        str='Lai_time')
      call print_time(modis%FparLai_QC_time, str='FparLai_QC_time')
      call print_time(modis%LaiStdDev_time,  str='LaiStdDev_time')
      write(string1,*)'All Lai times did not match.'
      call error_handler(E_ERR,'stow_observation',string1, source, revision, revdate)
   endif
   qc    = real(bits567,r8)
   obval = modis%Lai
   oerr  = min(0.1_r8, modis%LaiStdDev)

   call create_3d_obs(latitude, longitude, 0.0_r8, VERTISSURFACE, obval, &
               MODIS_LEAF_AREA_INDEX, oerr, oday, osec, qc, obs)
   call add_obs_to_seq(obs_seq, obs, modis%time_obs, prev_obs, prev_time, first_obs)
endif

! Time to reinitialize the modis structure data values

modis%FparExtra_QC      = MISSING_R8
modis%FparExtra_QC_time = set_time(0,0)
modis%FparLai_QC        = MISSING_R8
modis%FparLai_QC_time   = set_time(0,0)
modis%FparStdDev        = MISSING_R8
modis%FparStdDev_time   = set_time(0,0)
modis%Fpar              = MISSING_R8
modis%Fpar_time         = set_time(0,0)
modis%LaiStdDev         = MISSING_R8
modis%LaiStdDev_time    = set_time(0,0)
modis%Lai               = MISSING_R8
modis%Lai_time          = set_time(0,0)

end subroutine stow_observation



subroutine SetWordPositions( bigstring, solo, positions )
character(len=*),        intent(in)  :: bigstring
character,               intent(in)  :: solo
integer, dimension(:,:), intent(out) :: positions

integer :: nwords
integer :: i, iword

nwords = size(positions,1)

iword = 1 
positions(iword,1) = 1

! First, simply count up the number of delimiters 
do i = 1,len_trim(bigstring)
   if (bigstring(i:i) == solo) then
      positions(iword,2) = i-1
      iword = iword + 1
      positions(iword,1) = i+1
   endif
enddo
positions(iword,2) = len_trim(bigstring)

if ( 0 == 99 ) then  ! Test. Passes with gfortran compiler
   write(*,*)'There are ',nwords,' in this line: '
   write(*,*)'[',trim(bigstring),']'
   do i = 1,nwords
      write(*,*)'word ',i,' is [',bigstring(positions(i,1):positions(i,2)),']'
   enddo
endif

end subroutine SetWordPositions



end program MOD15A2_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
