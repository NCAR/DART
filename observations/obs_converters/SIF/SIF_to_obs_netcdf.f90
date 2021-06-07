! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program SIF_to_obs_netcdf

!-------------------------------------------------------------------------------
!
! SIF_to_obs_netcdf - reads monthly 'harmonized' SIF  
!                     observations and writes a DART obs_seq file.
!
! SIF observation description: Harmonized long term SIF data product(SIF005)
! which combines GOME-2 and SCIAMACHY SIF retrievals, combined with MODIS
! data to produce single, continuous product (global,monthly, 0.05 degree,
! 2002-2018; Wen et al., 2020, RSE).   Downloaded from:
! https://climatesciences.jpl.nasa.gov/sif/download-data/level-3
!  
! SIF variables provided: Provides daily average SIF value (mW/m^2/nm/sr)
! centered on 740 nm wavelength, and SIF standard deviation with same units.


use         types_mod, only : r8, missing_r8

use      location_mod, only : VERTISUNDEF, set_location

use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              nmlfileunit, do_nml_file, do_nml_term, &
                              find_namelist_in_file, check_namelist_read, &
                              error_handler, E_ERR, E_MSG, &
                              open_file, close_file, file_exist

use  time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, set_time, &
                              set_date, print_date, get_time, print_time

use  obs_sequence_mod, only : obs_sequence_type, obs_type, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data, &
                              destroy_obs_sequence, set_obs_def, &
                              set_obs_values, set_qc

use      obs_kind_mod, only : HARMONIZED_SIF

use  obs_def_land_mod, only : set_SIF_wavelength

use       obs_def_mod, only : obs_def_type, set_obs_def_type_of_obs, &
                              set_obs_def_location, set_obs_def_time, &
                              set_obs_def_error_variance, set_obs_def_key

use obs_utilities_mod, only : add_obs_to_seq

use netcdf_utilities_mod, only : nc_open_file_readonly, nc_close_file, &
                                 nc_get_dimension_size, nc_get_variable, &
                                 nc_variable_attribute_exists, &
                                 nc_get_attribute_from_variable

implicit none

character(len=*), parameter :: source  = 'SIF_to_obs_netcdf.f90'
character(len=512) :: string1,string2,string3

!-------------------------------------------------------------------------------
! Namelist input with default values

character(len=256) :: input_file_list = 'file_list.txt'
integer            :: verbose = 0
integer            :: wavelength = 740

namelist /SIF_to_obs_netcdf_nml/ input_file_list, verbose, wavelength

!-------------------------------------------------------------------------------

integer, parameter :: MAX_NUM_INPUT_FILES = 500
integer            :: num_input_files = 0
integer            :: ifile

character(len=256), dimension(MAX_NUM_INPUT_FILES) :: filename_seq_list

character(len=256) :: output_name

integer :: iocode, iunit
integer :: osec, oday
integer :: year, month, day
integer :: num_copies, num_qc
integer :: key
           
logical  :: first_obs

real(r8), allocatable, dimension(:,:) :: SIF
real(r8), allocatable, dimension(:,:) :: SIF_ERR
integer,  allocatable, dimension(:,:) :: EVI_Quality
real(r8), allocatable, dimension(:) :: lon
real(r8), allocatable, dimension(:) :: lat

integer :: ncid, nlat, nlon, i, j, num_new_obs

real(r8) :: qc
real(r8) :: rlat, rlon
real(r8) :: FillValue 
logical  :: has_FillValue = .false.

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time

!-------------------------------------------------------------------------------
! start of executable code

call initialize_utilities(source)

! time setup
call set_calendar_type(GREGORIAN)

! Read the namelist entry
call find_namelist_in_file("input.nml", "SIF_to_obs_netcdf_nml", iunit)
read(iunit, nml = SIF_to_obs_netcdf_nml, iostat = iocode)
call check_namelist_read(iunit, iocode, "SIF_to_obs_netcdf_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=SIF_to_obs_netcdf_nml)
if (do_nml_term()) write(     *     , nml=SIF_to_obs_netcdf_nml)

num_input_files = Check_Input_Files(input_file_list, filename_seq_list) 
write(*,*)' There are ',num_input_files,' input files.'

! call the initialization code, and initialize two empty observation types
! each observation in this series will have a single observation value 
! and a quality control flag.  Set the DART data quality control, 0 is good. 
! increasingly larger QC values are more questionable quality data

key        = 0
num_copies = 1
num_qc     = 2
qc         = 0.0_r8

call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

!-------------------------------------------------------------------------------
! Loop over all the input data files.

FileLoop: do ifile = 1,num_input_files
   
   first_obs = .true.

   ! Fills up matrix with global SIF Harmonized data
   ncid = nc_open_file_readonly(filename_seq_list(ifile), source)
   
   nlat  = nc_get_dimension_size(ncid,'lat')
   nlon  = nc_get_dimension_size(ncid,'lon')

   allocate( SIF(nlon,nlat) )
   allocate( SIF_ERR(nlon,nlat) )
   allocate( EVI_Quality(nlon,nlat) )
   allocate( lon(nlon))
   allocate( lat(nlat))

   call nc_get_variable(ncid, 'lat', lat, source)
   call nc_get_variable(ncid, 'lon', lon, source)
   call nc_get_variable(ncid, 'SIF_740_daily_corr', SIF, source)
   call nc_get_variable(ncid, 'SIF_740_daily_corr_SD', SIF_ERR, source)
   call nc_get_variable(ncid, 'EVI_Quality', EVI_Quality, source)

   has_FillValue = nc_variable_attribute_exists(ncid,'SIF_740_daily_corr','_FillValue')

   if ( has_FillValue ) then
      call nc_get_attribute_from_variable(ncid,'SIF_740_daily_corr','_FillValue',FillValue)
   endif

   ! Convert data units etc.
   where(lon < 0.0_r8) lon = lon + 360.0_r8
   SIF_ERR=SIF_ERR**2

   call nc_close_file(ncid, source)
   
   ! create a new, empty obs_seq file.  you must give a max limit
   ! on number of obs.  increase the size if too small.
   num_new_obs=nlon*nlat
   
   call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)

   ! the first one needs to contain the string 'observation' and the
   ! second needs the string 'QC'.
   call set_copy_meta_data(obs_seq, 1, 'observation')
   call   set_qc_meta_data(obs_seq, 1,     'Data QC')
   call   set_qc_meta_data(obs_seq, 2, 'EVI_Quality')
   
   call get_yyyymmdd(filename_seq_list(ifile), year, month, day)

   write(output_name,'(''obs_seq.'',i4.4,''-'',i2.2,''-'',i2.2,''.out'')') &
                     year, month, day
 
   time_obs = set_date(year,month,day,0,0,0)
   call get_time(time_obs, osec, oday)

   ! A little helpful logging
   if (verbose > 0) then
      write(*,*)'converting "',trim(filename_seq_list(ifile)),'"'
      write(string1,*)'"'//trim(output_name)//'" date is'
      write(string2,*)'"'//trim(output_name)//'" time is'
      call print_date(time_obs, str=string1)
      call print_time(time_obs, str=string2)
   endif

   LATLOOP: do j=1,nlat
   
      ! ensure the lat values are in range
      rlat = lat(j)
      if ( rlat >  90.0_r8 .or. rlat <  -90.0_r8 ) cycle LATLOOP
   
      LONLOOP: do i=1,nlon
     
         ! ensure the lat values are in range 
         rlon = lon(i)   
         if ( rlon > 360.0_r8 .or. rlon <    0.0_r8 ) cycle LONLOOP

         ! Skip anything marked as such
         if ( has_FillValue .and. SIF(i,j) == FillValue ) cycle LONLOOP

         ! reject obviously bad observations
         if ( SIF(i,j) <  0.0_r8 ) cycle LONLOOP

         ! reject observations with questionable QC values
         qc  = DecipherEVI(EVI_Quality(i,j),i,j)
         if ( qc >= 50 ) cycle LONLOOP

         ! Update the metadata array holding the wavelengths.
         ! Return the position in that array.
         key = set_SIF_wavelength(wavelength)
 
         ! make an obs derived type, and then add it to the sequence
         ! VERTISUNDEF declares that the entire column will have the same
         ! localization distance - only dependent on the horizontal distance.
         ! VERTISSURFACE is also horizontal distance only.

         call create_3d_obs(lat(j), lon(i), 0.0_r8, VERTISUNDEF, SIF(i,j), &
                 HARMONIZED_SIF, SIF_ERR(i,j), oday, osec, qc, EVI_Quality(i,j), key, obs)

         call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
      
      enddo LONLOOP 
   enddo LATLOOP

   ! if we added any obs to the sequence, write it out to a file now.
   if ( get_num_obs(obs_seq) > 0 ) then
      if (verbose > 0) then
         write(string1,*)'"'//trim(output_name)//'" has ', &
                        get_num_obs(obs_seq), ' observations.'
         call error_handler(E_MSG, source, string1)
      endif
      call write_obs_seq(obs_seq, output_name)
   endif

   deallocate(SIF)
   deallocate(SIF_ERR)
   deallocate(EVI_Quality)
   deallocate(lon)
   deallocate(lat)

   call destroy_obs_sequence(obs_seq)

enddo FileLoop

! end of main program
call finalize_utilities()

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------


! Read a file containing a list of filenames and checks that they exist 

function Check_Input_Files(input_list, output_list)

character(len=*),               intent(in)  :: input_list
character(len=*), dimension(:), intent(out) :: output_list
integer                                     :: Check_Input_Files

character(len=256) :: input_line
character(len=256) :: ladjusted
integer, parameter :: MAXLINES = 500
integer :: iline

iunit = open_file(trim(input_list), 'formatted', 'read')

Check_Input_Files = 0
FileNameLoop: do iline = 1,MAXLINES ! a lot of lines 

   ! read in entire text line into a buffer
   read(iunit, "(A)", iostat=iocode) input_line
   if (iocode > 0) then 
      write(string1,*) 'While reading ', trim(input_list)
      write(string2,*) 'got read code (iostat) = ', iocode,' around line ',iline
      call error_handler(E_ERR, 'Check_Input_Files', string1, &
                    source,text2=string2)
   elseif (iocode < 0) then 
      ! Normal end of file
      exit FileNameLoop
   else
      Check_Input_Files = Check_Input_Files + 1

      ladjusted = adjustl(input_line)
      if ( file_exist(trim(ladjusted)) ) then
         output_list(Check_Input_Files) = trim(ladjusted)
      else
         write(string1,*)'file does not exist.'
         call error_handler(E_ERR,'Check_Input_Files',&
          string1,source,text2=trim(ladjusted))
      endif
   endif

enddo FileNameLoop

if (Check_Input_Files >= MAXLINES-1 ) then
   write(string1,*)'Too many files to process. Increase MAXLINES and try again.'
   call error_handler(E_ERR,'Check_Input_Files',string1,source)
endif

call close_file(iunit)

end function Check_Input_Files


!-------------------------------------------------------------------------------
! get the year,month,day from the name of input files

subroutine get_yyyymmdd(filename, year, month, day)

character(len=*), intent(in) :: filename
integer,          intent(out) :: year, month, day

! The filenames END with "SIF005_201808.nc"

integer :: strlen, index1, io

strlen = len_trim(filename)
index1 = strlen - 9 + 1  !    YYYYMM.nc is 9 chars

! Force day to be mid-month
day=15

read(filename(index1:strlen),"(I4,I2,3X)",iostat=io) year, month   
if (io /= 0) then
   write(string1,*)'Unable to read year, month, day from filename.'
   write(string2,*)'filename was "'//trim(filename)//'"'
   write(string3,*)'Expected last 10 chars to contain the date, they are "' &
                   & //filename(index1:strlen)//'"'
   call error_handler(E_ERR, 'get_yyymmdd', string1, &
              source, text2=string2, text3=string3)
endif

end subroutine get_yyyymmdd


!-------------------------------------------------------------------------------
!   create an observation type from all the pieces
!
!   NOTE: assumes the code is using the threed_sphere locations module, 
!         that the observation has a single data value and TWO qc values.
!         Extended from the observations/utilities/obs_utilities_mod.f90
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
!    qc    - (converted, NCEP-like) quality control value
!    evi   - original (binary-coded-decimal) quality control value 
!    key   - index into obs_def_land_mod arrays holding the metadata
!    obs   - observation type

subroutine create_3d_obs(lat, lon, vval, vkind, obsv, okind, oerr, day, sec, &
                         qc, evi, key, obs)

integer,        intent(in)    :: okind, vkind, day, sec
real(r8),       intent(in)    :: lat, lon, vval, obsv, oerr, qc
integer,        intent(in)    :: evi
integer,        intent(in)    :: key
type(obs_type), intent(inout) :: obs

real(r8)           :: obs_val(1), qc_val(2)
type(obs_def_type) :: obs_def

call set_obs_def_location(obs_def, set_location(lon, lat, vval, vkind))
call set_obs_def_type_of_obs(obs_def, okind)
call set_obs_def_time(obs_def, set_time(sec, day))
call set_obs_def_error_variance(obs_def, oerr * oerr)
call set_obs_def_key(obs_def, key)
call set_obs_def(obs, obs_def)

obs_val(1) = obsv
call set_obs_values(obs, obs_val)
qc_val(1)  = qc
qc_val(2)  = evi
call set_qc(obs, qc_val)

end subroutine create_3d_obs

!-------------------------------------------------------------------------------
! Decipher the binary-coded-integer into a DART-compatible error value.
! The DART-compatible error scheme is based on NCEP-like error codes.
! 0.0 = best quality
! 1.0 = lesser
! ... etc.
! By not rejecting the observations at this step, the 'input_qc_threshold'
! namelist value can be used to test whether or not 'suspect' observations
! improve the result or not.

function DecipherEVI(evi,ilon,jlat) result (qc)

integer, intent(in) :: evi
integer, intent(in) :: ilon,jlat
real(r8)            :: qc

character(len=64) :: binary_string
character(len=6) :: bits543210
integer :: ndigits

! QC Interpretation
! https://cornell.app.box.com/s/gkp4moy4grvqsus1q5oz7u5lc30i7o41/file/687233625170
! https://lpdaac.usgs.gov/documents/103/MOD13_User_Guide_V6.pdf
! From Table 5:
! Bits 0-1
! 00 produced with good quality
! 01 produced, but check other QA
! 10 produced, but most likely cloudy
! 11 not produced, due to other reasons than clouds
! Bits 2-5
! 0000 Highest quality
! 0001 Lower quality
! 0010 Decreasing quality
! 0100 Decreasing quality
! 1000 Decreasing quality
! 1001 Decreasing quality
! 1010 Decreasing quality
! 1100 Lowest     quality
! 1101 Quality so low that it is not useful
! 1110 L1B data faulty
! 1111 Not useful for any other reason/not processed

! Determine if least significant bits are right or left
! write(*, '(''8 is '',B4)')8
! write(*, '(''7 is '',B4)')7
! write(*, '(''6 is '',B4)')6
! write(*, '(''5 is '',B4)')5
! write(*, '(''4 is '',B4)')4
! write(*, '(''3 is '',B4)')3
! write(*, '(''2 is '',B4)')2
! write(*, '(''1 is '',B4)')1
! write(*, '(''0 is '',B4)')0

! write(*,*)
! write(*,'(''1234567890123456'')')
! write(*, '(I0,1x,B0)') evi,evi

write(binary_string, '(B0)') evi

ndigits = len_trim(binary_string)

! Grab the least significant six digits (as a character string)
bits543210 = binary_string((ndigits-5):ndigits)

! This interpretation is ENTIRELY SUBJECTIVE ... 

select case (bits543210)
         !543210   ! bits 2-5                 | bits 0-1
   case ('000000') ! 0000 ... Highest Quality | 00 V1 produced with qood quality
      qc = 0.0_r8
   case ('000100') ! 0001 ... Lower Quality   | 00 V1 produced with good quality
      qc = 1.0_r8
   case ('001000') ! 0010 ... Lower Quality   | 00 ditto
      qc = 2.0_r8
   case ('010000') ! 0100 ... Lower Quality   | 00 ditto
      qc = 3.0_r8
   case ('100000') ! 1000 ... Lower Quality   | 00 ditto
      qc = 4.0_r8
   case ('100100') ! 1001 ... Lower Quality   | 00 ditto
      qc = 5.0_r8  
   case ('101000') ! 1010 ... Lower Quality   | 00 ditto
      qc = 6.0_r8
   case ('110000') ! 1100 ... Lowest Quality  | 00 ditto
      qc = 7.0_r8



   case ('000001') ! 0000 ... Highest Quality | 01 V1 produced, but check with other QA
      qc =10.0_r8
   case ('000101') ! 0001 ... Lower Quality   | 01 V1 produced, but check with other QA
      qc =11.0_r8
   case ('001001') ! 0010 ... Lower Quality   | 01 ditto
      qc =12.0_r8
   case ('010001') ! 0100 ... Lower Quality   | 01 ditto
      qc =13.0_r8
   case ('100001') ! 1000 ... Lower Quality   | 01 ditto
      qc =14.0_r8
   case ('100101') ! 1001 ... Lower Quality   | 01 ditto
      qc =15.0_r8
   case ('101001') ! 1010 ... Lower Quality   | 01 ditto
      qc =16.0_r8
   case ('110001') ! 1100 ... Lowest Quality  | 01 ditto
      qc =17.0_r8

   case default
      qc = 50.0_r8 ! do not use
end select

! Large negative values will be thrown out

if (evi < 0) qc=100 

if (verbose > 1) &
   write(*,'(''ilon,jlat,EVI,... '',3(I7,1x),A,F10.2)')ilon, jlat, evi, bits543210, qc

end function DecipherEVI

 
end program SIF_to_obs_netcdf 

