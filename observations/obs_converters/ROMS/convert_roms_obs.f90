! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> program that converts ROMS observations PLUS FORWARD OPERATOR VALUES
!> (estimated values computed during the ROMS run).  this program
!> needs to read 1 obs input file to get the original obs value, time,
!> location, error, etc; plus N (ensemble size) ROMS "mod" files to
!> extract the forward operators.  possibly also a third type of file
!> to get the grid information - we are currently hoping that we can
!> get lat/lon/depth directly from the ROMS input obs file, but if not
!> we have to add more code. the input obs files currently has X,Y,Z
!> values relative to the ROMS grid - we'd have to open a grid file and
!> compute the lat/lon values from that.

program convert_roms_obs

use         types_mod, only : r8, missing_r8, obstypelength

use     utilities_mod, only : initialize_utilities, finalize_utilities,            &
                              error_handler, do_nml_term, do_nml_file,             &
                              E_ERR, E_WARN, E_MSG, logfileunit, nmlfileunit,      &
                              find_namelist_in_file, check_namelist_read,          &
                              open_file, close_file, find_textfile_dims,           &
                              file_to_text, do_output, set_filename_list

use netcdf_utilities_mod, only : nc_check

use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                              increment_time, get_time, operator(-), GREGORIAN, &
                              print_time, print_date

use      location_mod, only : VERTISHEIGHT, location_type, get_location

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, get_obs_def, &
                              set_copy_meta_data, set_qc_meta_data, set_obs_def

use      obs_kind_mod, only : get_quantity_for_type_of_obs, &
                              get_name_for_type_of_obs, get_index_for_type_of_obs

use       obs_def_mod, only : set_obs_def_external_FO, obs_def_type

use          sort_mod, only : index_sort

use obs_utilities_mod, only : getvar_real, add_obs_to_seq, &
                              create_3d_obs, getvar_int, getdimlen, &
                              query_varname, set_missing_name

use    random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

use         model_mod, only : static_init_model, end_model, &
                              get_time_information, get_location_from_ijk

use           netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer, parameter :: MAX_ENS = 100

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

integer  :: iunit, io, darttype, num_input_files, rc
integer  :: ncid, nobs, n, i, oday, osec, nobs2
logical  :: file_exist, first_obs

character(len=512) :: string1, string2, string3

integer, parameter :: MAX_TYPES = 100
character(len=128) :: roms_prov_strings(MAX_TYPES)
integer            :: roms_prov_integers(MAX_TYPES) = -1
integer            :: dart_type_integers(MAX_TYPES) = -1
integer            :: roms_prov_count = 0

! mapping between roms type numbers and dart types.
! use the roms number as the index; the value is the dart kind

type(obs_def_type) :: obs_def

real(r8), allocatable :: lat(:), lon(:), depth(:), ovar(:),  &
                         oval(:), raw_qc(:), combined_qc(:), &
                         temp_fo(:), forw_ops(:,:)
integer, allocatable  :: oprov(:)
real(r8) :: missing

! roms grid in i,j,k space - needs model_mod routine to convert
! to lat/lon/depth

real(r8), allocatable :: Xgrid(:), Ygrid(:), Zgrid(:)

type(random_seq_type) :: my_rand

type(time_type), allocatable :: otim(:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time

! Namelist information and defaults

integer  :: ens_size                                  = 1
character(len=256) :: roms_mod_obs_files(MAX_ENS)     = ''
character(len=256) :: roms_mod_obs_filelist           = 'filelist.txt'
character(len=256) :: dart_output_obs_file            = 'obs_seq.out'
logical  :: append_to_existing                        = .false.
logical  :: use_precomputed_values                    = .true.
logical  :: locations_in_IJK                          = .false.
logical  :: add_random_noise                          = .false.
real(r8) :: pert_amplitude                            = 0.01
integer  :: verbose                                   = 0
character(len=256) :: type_translations(2, MAX_TYPES) = 'NULL'

!>@todo FIXME -
!> we could get the ens size from the number of input files
!> but we can also ask users to tell us how many should be there
!> so we can error out if the counts don't match.

!>@todo FIXME -
!> TJH: I vote to remove the roms_obs_files() mechanism since roms_mod_obs_filelist
!> is so much more flexible. The roms_obs_files variable just mucks up the 
!> namelist print and log files.

!>@todo FIXME -
!> TJH: The type_translations table is currently implemented as an explicit
!> set of type translations ... it _could_ be made to exist as a superset with
!> multiple ROMS variants that relate to the same DART TYPE. hernan and tim
!> prefer the superset, Chris and Nancy prefer the subset ... Andy?

!>@todo FIXME -
!> TJH: could remove the locations_in_IJK variable by simply checking 
!> if the obs_lon,obs_lat variables exist. 
!> If they do exist, use them; if not - try the IJK method.

namelist /convert_roms_obs_nml/ &
   ens_size, &
   roms_mod_obs_files, &
   roms_mod_obs_filelist, &
   dart_output_obs_file, &
   append_to_existing, &
   use_precomputed_values, &
   add_random_noise, &
   pert_amplitude, &
   verbose, &
   type_translations

!------------
! start of executable code
!------------

call initialize_utilities('convert_roms_obs')

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'convert_roms_obs_nml', iunit)
read(iunit, nml = convert_roms_obs_nml, iostat = io)
call check_namelist_read(iunit, io, 'convert_roms_obs_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=convert_roms_obs_nml)
if (do_nml_term()) write(     *     , nml=convert_roms_obs_nml)

if (locations_in_IJK) then
   write(string1,*)'Using locations in IJK is currently unsupported.'
   write(string2,*)'The routine to convert from IJK to lat-lon-depth is not finished.'
   call error_handler(E_ERR, 'convert_roms_obs', string1, &
              source, revision, revdate, text2=string2)
endif

! when set_filename_list() returns, roms_mod_obs_files contains the file list
! whether they were specified directly in the namelist or in a separate file.
num_input_files = set_filename_list(roms_mod_obs_files, roms_mod_obs_filelist, 'convert_roms_obs')

! Initialize repeatable random sequence
call init_random_seq(my_rand)

! use the machinery in the ROMS model_mod to read in the grid
! and handle any needed conversions

call static_init_model()

call nc_check( nf90_open(roms_mod_obs_files(1), nf90_nowrite, ncid), &
               'convert_roms_obs', 'opening first file "'//trim(roms_mod_obs_files(1))//'"')

call getdimlen(ncid, "datum", nobs)

allocate(oprov(nobs))
allocate(oval(nobs))
allocate(ovar(nobs))
allocate(otim(nobs))
allocate(raw_qc(nobs))
allocate(combined_qc(nobs))
allocate(temp_fo(nobs))
allocate(forw_ops(ens_size, nobs))
allocate(lat(nobs))
allocate(lon(nobs))
allocate(depth(nobs))

! read in the data arrays
!>@todo FIXME make sure the missing values in the depth arrays are handled
! correctly further downstream

call set_missing_name('missing_value')

if (locations_in_IJK) then
   allocate(Xgrid(nobs))
   allocate(Ygrid(nobs))
   allocate(Zgrid(nobs))
   call getvar_real(ncid, "obs_Xgrid", Xgrid)
   call getvar_real(ncid, "obs_Ygrid", Ygrid)
   call getvar_real(ncid, "obs_Zgrid", Zgrid, dmiss=missing)
   where(Zgrid == missing) Zgrid = MISSING_R8
else
   call getvar_real(ncid, "obs_lat", lat)
   call getvar_real(ncid, "obs_lon", lon)
   call getvar_real(ncid, "obs_depth", depth, dmiss=missing)
   where(depth == missing)
       depth = MISSING_R8
   end where
endif

call getvar_real(ncid, "obs_value", oval)
call getvar_real(ncid, "obs_error", ovar)  ! already squared, so variance not stddev
call getvar_int( ncid, "obs_provenance", oprov)

! these come back as dart time types
call get_time_information(roms_mod_obs_files(1), ncid, "obs_time", "datum", all_times=otim)

! right now i'm getting the mapping from the string attribute on the
! obs_provenance netcdf variable, NOT the global attribute "obs_provenance"
! because that one is harder to parse, and will be harder for users to give
! us matching strings in the namelist.

call get_provenance_maps(ncid, roms_prov_count, roms_prov_strings, roms_prov_integers)
call get_translation_table()

! we can close the input obs file here, and then loop over the
! ensemble of mod files, reading in the expected values and QCs only

call nc_check( nf90_close(ncid) , &
               'convert_roms_obs', 'closing file '//trim(roms_mod_obs_files(1)))

! Set the incoming data quality control.  this combines
! the qcs (scale) from all ensembles to make a combined value.

combined_qc(:) = 0.0_r8

! loop over the ensemble of mod files, opening each file and
! reading in both the expected values per obs, and also the qc

call set_missing_name('_FillValue')

do i = 1, ens_size

   write(string1,'(''opening file '',i3,'' "'',A,''"'')')i,trim(roms_mod_obs_files(i))

   call nc_check( nf90_open(roms_mod_obs_files(i), nf90_nowrite, ncid), &
                  'convert_roms_obs', string1)

   call getdimlen(ncid, "datum", nobs2)

   if (nobs2 /= nobs) then
      write(string1,*)'require "datum" dimensions to be equal.'
      write(string2,*)'"datum" dimension is ',nobs,  ' in ',trim(roms_mod_obs_files(1))
      write(string3,*)'"datum" dimension in ',nobs2, ' in ',trim(roms_mod_obs_files(i))
      call error_handler(E_ERR, 'convert_roms_obs', string1, &
                 source, revision, revdate, text2=string2, text3=string3)
   endif

   call getvar_real(ncid, "NLmodel_value", temp_fo, dmiss=missing)
   where (temp_fo == missing) temp_fo = MISSING_R8
   call getvar_real(ncid, "obs_scale",     raw_qc,  dmiss=missing)
   where (raw_qc == missing) raw_qc = huge(missing)

   forw_ops(i, :) = temp_fo(:)

   ! 'max' works on arrays ... who knew ...
   combined_qc = max(combined_qc, raw_qc)

   call nc_check( nf90_close(ncid) , &
                  'convert_roms_obs', 'closing file '//trim(roms_mod_obs_files(i)))

enddo

first_obs = .true.

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

inquire(file=dart_output_obs_file, exist=file_exist)

if ( file_exist .and. append_to_existing ) then

  ! existing file found, append to it
  call read_obs_seq(dart_output_obs_file, 0, 0, nobs, obs_seq)

else

  ! create a new one
  call init_obs_sequence(obs_seq, num_copies, num_qc, nobs)
  do i = 1, num_copies
    call set_copy_meta_data(obs_seq, i, 'observation')
  end do
  do i = 1, num_qc
    call set_qc_meta_data(obs_seq, i, 'Data QC')
  end do

endif

obsloop: do n = 1, nobs

   if (combined_qc(n) == huge(missing)) cycle obsloop

   ! time of observation
   time_obs = otim(n)

   if (first_obs) then
      call print_time(time_obs, str='obs time is ')
      call print_date(time_obs, str='obs date is ')
   endif

   call convert_romstype_to_darttype(oprov(n), darttype)
   if (darttype < 0) then
      write(string1, *) 'Observation ', n, ' failed to convert ROMS type to DART type, skipping'
      write(string2, *) 'ROMS provenance value was: ', oprov(n)
      call error_handler(E_MSG, 'convert_roms_obs', string1, text2=string2)
      cycle obsloop
   endif

   if (locations_in_IJK) then
      rc = convert_ijk_to_latlondepth(Xgrid(n), Ygrid(n), Zgrid(n), darttype, &
                                      lon(n), lat(n), depth(n))
      if (rc /= 0) then
         write(string1, *) 'Observation ', n, ' failed to convert I,J,K to lat/lon, skipping'
         write(string2, *) 'I,J,K values and return code were: ', Xgrid(n), Ygrid(n), Zgrid(n), rc
         call error_handler(E_MSG, 'convert_roms_obs', string1, text2=string2)
         cycle obsloop
      endif
   endif

   if (depth(n) == MISSING_R8) cycle obsloop

   ! check the lat/lon values to see if they are ok
   if (( lat(n) >  90.0_r8 .or. lat(n) <  -90.0_r8 ) .or. &
       ( lon(n) > 360.0_r8 .or. lon(n) < -180.0_r8 )) then
      write(string1, *) 'Observation ', n, ' has out-of range latitude/longitudes, skipping'
      write(string2, *) 'latitude/longitude values were: ', lat(n), lon(n)
      call error_handler(E_MSG, 'convert_roms_obs', string1, text2=string2)
      cycle obsloop
   endif

   if ( lon(n) < 0.0_r8 )  lon(n) = lon(n) + 360.0_r8

   !>@todo FIXME ... depth may be MISSING_R8 ... in that the parent variable
   !> "obs_depth" has a 'missing_value' attribute

   ! extract actual time of observation in file into oday, osec.
   call get_time(time_obs, osec, oday)

   call create_3d_obs(lat(n), lon(n), depth(n), VERTISHEIGHT, oval(n), &
                      darttype, sqrt(ovar(n)), oday, osec, combined_qc(n), obs)

   if (add_random_noise) then
      do i=1, ens_size
         forw_ops(i, n) = random_gaussian(my_rand, forw_ops(i, n), pert_amplitude)
      enddo
   endif

   call get_obs_def(obs, obs_def)
   call set_obs_def_external_FO(obs_def, .true., use_precomputed_values, &
                                n, ens_size, forw_ops(:, n))
   call set_obs_def(obs, obs_def)

   call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

end do obsloop

deallocate(oprov, oval, ovar, otim, raw_qc, combined_qc, temp_fo, forw_ops, lat, lon, depth)

if (allocated(Xgrid)) deallocate(Xgrid)
if (allocated(Ygrid)) deallocate(Ygrid)
if (allocated(Zgrid)) deallocate(Zgrid)

! if we added any obs to the sequence, write it now.
if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, dart_output_obs_file)

call end_model()

!>@todo ... write out a summary of how many observations were successfully converted.

! end of main program
call finalize_utilities()

contains

!-----------------------------------------------------------------------
!>
!> roms observation types are in the 'obs_provenance' integer array
!> and the table of contents in the global attribute 'obs_provenance'
!> (formatted as a single string with embedded newlines)
!>
!> dart observation specific types are numeric values - see the comment
!> below for the current list.  we can easily add more specific types
!> by editing DART/obs_def/obs_def_ocean_mod.f90 file.  follow the
!> examples already there.  the first column is where you define a new
!> obs name.  the second column has to be an existing type.  the third
!> column should be 'COMMON_CODE' like most of the others.  it also must
!> be a fortran comment, so lines start with !

subroutine convert_romstype_to_darttype(roms_prov, dart_type)

integer, intent(in)  :: roms_prov
integer, intent(out) :: dart_type

integer :: i

if (roms_prov < 1) then
   write(string1,*)'unknown roms provenance integer ',roms_prov
   call error_handler(E_ERR, 'convert_romstype_to_darttype',string1, &
                      source, revision, revdate)
endif

DUMBLOOP: do i = 1,roms_prov_count

   if (roms_prov == roms_prov_integers(i)) then
       dart_type  = dart_type_integers(i)
       return
   endif

enddo DUMBLOOP

if (roms_prov > roms_prov_count) then
   call error_handler(E_ERR, 'convert_romstype_to_darttype', 'bad roms type, greater than limit', &
                      source, revision, revdate)
endif

dart_type = dart_type_integers(roms_prov)

end subroutine convert_romstype_to_darttype

!-----------------------------------------------------------------------
!>
!> ISTATUS : 0 - all went well
!> ISTATUS : 3 - vertical index is a missing_value
!> all other error codes inherited from model_mod:get_location_from_ijk()

function convert_ijk_to_latlondepth(r_iloc, r_jloc, r_kloc, dart_type, lon, lat, depth)

real(r8), intent(in)  :: r_iloc
real(r8), intent(in)  :: r_jloc
real(r8), intent(in)  :: r_kloc
integer,  intent(in)  :: dart_type
real(r8), intent(out) :: lat
real(r8), intent(out) :: lon
real(r8), intent(out) :: depth
integer :: convert_ijk_to_latlondepth

integer :: dart_kind, status
type(location_type) :: dart_location
real(r8) :: dummy(3)

!>@todo FIXME ... r_kloc may be MISSING_R8 ... in that the parent variable
!> "obs_depth" in the IJK world has a 'missing_value' attribute
if (r_kloc == MISSING_R8) then
   convert_ijk_to_latlondepth = 3
   return 
endif

! convert from specific obs type to a generic obs kind
dart_kind = get_quantity_for_type_of_obs(dart_type)

! if we didn't get lats/lons directly, call the ROMS model_mod code
! to convert the grid indices into lat/lon/depth

status = get_location_from_ijk(r_iloc, r_jloc, r_kloc, dart_kind, dart_location)

if (status /= 0) then
   convert_ijk_to_latlondepth = status
   return
endif

dummy = get_location(dart_location)
lon   = dummy(1)
lat   = dummy(2)
depth = dummy(3)

! things went ok
convert_ijk_to_latlondepth = 0

end function convert_ijk_to_latlondepth

!-----------------------------------------------------------------------
!> the global attribute 'obs_provenance' relates the integer codes with
!> the character strings that humans prefer.
!> The character string MUST IDENTICALLY match the 
!> input.nml:&convert_roms_obs_nml:type_translations that relate those
!> character strings to DART observation TYPES. 

subroutine get_provenance_maps(ncid, flag_count, flag_names, flag_vals)

integer,          intent(in)  :: ncid
integer,          intent(out) :: flag_count
character(len=*), intent(out) :: flag_names(:)
integer,          intent(out) :: flag_vals(:)

integer :: io

integer, parameter :: MAX_PROV_STRLEN = 1000
character(len=MAX_PROV_STRLEN) :: provenance_str
character(len=MAX_PROV_STRLEN) :: provenance_str2

integer :: stringlength, num_newlines, pos1, pos2
integer :: newline_pos(MAX_PROV_STRLEN)

io = nf90_get_att(ncid, NF90_GLOBAL, 'obs_provenance', provenance_str)
call nc_check(io, 'get_provenance_maps', 'read global attribute "obs_provenance"')

io = nf90_inquire_attribute(ncid, NF90_GLOBAL, 'obs_provenance', len=stringlength)
call nc_check(io, 'get_provenance_maps', 'inquire_attribute "obs_provenance:flag_values"')

if (stringlength > len(provenance_str)) then
   write(string1,*)'GLOBAL attribute obs_provenance has length ',stringlength
   write(string2,*)'MAX_PROV_STRLEN is only ',MAX_PROV_STRLEN
   write(string3,*)'Must increase MAX_PROV_STRLEN and recompile.' 
   call error_handler(E_ERR, 'get_provenance_maps', string1, &
                         source, revision, revdate, text2=string2, text3=string3)
endif

io = nf90_get_att(ncid, NF90_GLOBAL, 'obs_provenance', provenance_str)
call nc_check(io, 'get_provenance_maps', 'get_att obs_provenance '//trim(roms_mod_obs_files(1)))

! The challenge is now to parse this one long array into lines with arbitrary lengths
! 1) count up the newline characters
! 2) record the location of the newlines in an array
! 3) parse

num_newlines = 0
do i=1, stringlength
   if (provenance_str(i:i) == new_line('A')) then
      num_newlines = num_newlines + 1
      newline_pos(num_newlines) = i
   endif
enddo

if (maxval(newline_pos) == stringlength) then
   ! There is a newline at the end - so there is nothing to do.
else
   ! The last line is unterminated.
   num_newlines = num_newlines + 1
   newline_pos(num_newlines) = stringlength + 1
endif

pos1 = 1
flag_count = 0
PARSELOOP : do i=1, num_newlines

   pos2 = newline_pos(i)

   if (pos2 == pos1) cycle PARSELOOP ! first line may be nothing

   flag_count = flag_count + 1

   provenance_str2 = trim(provenance_str(pos1+1:pos2-1))
   call get_flag_value(provenance_str2,flag_vals(flag_count))  
   call get_flag_name( provenance_str2,flag_names(flag_count))  

   pos1 = pos2

enddo PARSELOOP

if (verbose > 0) then

   string1 = 'Dumping out the obs_provenance map as read from the global attribute:'
   call error_handler(E_MSG, '', string1)

   do i=1, flag_count
      write(string1,'(''provenance integer '',i10,'' means "'',A,''"'')')flag_vals(i),trim(flag_names(i))
      call error_handler(E_MSG, '', string1)
   enddo
endif

end subroutine get_provenance_maps

!-----------------------------------------------------------------------
!> Given a string of the form
!>
!>  :obs_provenance = "\n",
!>          "    1: gridded AVISO sea level anomaly \n",
!>          "    2: gridded Aquarius SSS \n",
!>          "    4: XBT from Met Office \n",
!>          "    8: CTD from Met Office \n",
!>          "   16: ARGO floats \n",
!>          "   32: glider UCSD \n",
!>          "   64: gridded MODIS AQUA SCHL \n",
!>          "  128: blended satellite SST \n",
!>          "  256: CTD from CalCOFI \n",
!>          "  512: HF radar current UCSD \n",
!>          " 1024: glider La Push\n",
!>          " 2048: CTD from GLOBEC \n",
!>          " 4096: buoy, thermistor from Met Office \n",
!>          " 8192: ARGO floats from Met Office \n",
!>          "16384: glider MBARI \n",
!>          "32768: glider temperature NSF OOI \n",
!>          "65536: NOAA tide station sea level \n",
!>          "all others: mixed; based on sum(unique(obs_provenance))" ;
!>
!> the flag value 'ARGO floats' is the integer '16'
!> Make no assumptions about how many digits there are
!>
!> For the 'all others ...' assign a negative number.
!> That number should not appear as a provenance value and so those
!> observations will not match anything. 

subroutine get_flag_value(mystring,ivalue)

character(len=*), intent(in)  :: mystring
integer,          intent(out) :: ivalue

integer :: io, iend

iend    =   index(mystring,':') - 1
string1 = adjustl(mystring(1:iend))

read(string1,*,iostat=io) ivalue

! In at least one instance, there was no unique integer code referenced

if (io /= 0) ivalue = MISSING_R8 ! this will coerce to an integer

if (verbose > 1) then
   write(string1,*)' given ['//trim(mystring)//'] the flag value is ',ivalue
   call error_handler(E_MSG, '', string1)
endif

end subroutine get_flag_value

!-----------------------------------------------------------------------
!> Given a string of the form
!>  6: ARGO floats temperature from Met Office 
!> the flag name is everything beyond the colon 

subroutine get_flag_name(mystring,myflag)

character(len=*), intent(in)  :: mystring
character(len=*), intent(out) :: myflag

integer :: istart

istart  =   index(mystring,':') + 1
string1 = adjustl(mystring(istart:))
myflag  = trim(string1)

if (verbose > 1) then
   write(string1,*)' given ['//trim(mystring)//'] the flag name is "'//trim(myflag)//'"'
   call error_handler(E_MSG, '', string1)
endif

end subroutine get_flag_name

!-----------------------------------------------------------------------

function match_roms_flag(str_to_match)

character(len=*), intent(in) :: str_to_match
integer :: match_roms_flag

do i=1, roms_prov_count
   if ( str_to_match == roms_prov_strings(i)) then
      match_roms_flag = roms_prov_integers(i)
      return
   endif
enddo

match_roms_flag = -1

end function match_roms_flag

!-----------------------------------------------------------------------
!> user gives us a translation table from the namelist
!> each pair of strings should be:
!>    a roms string that matches the flags in the netcdf file
!>    a dart type string like "ARGO_TEMPERATURE"
!>
!> integer, intent(out) :: dart_type_integers(:)

subroutine get_translation_table()

integer :: i, j, entries, romstype, darttype

entries = 0

parseloop: do i=1, MAX_TYPES
   if (type_translations(1, i) == 'NULL') exit parseloop

   ! type_translations(1, i) is roms name - string
   ! type_translations(2, i) is dart type - string

   ! convert strings to corresponding IDs.
   ! there are no assumptions about the IDs other than they must be positive
   romstype = match_roms_flag(   type_translations(1, i))
   darttype = get_index_for_type_of_obs(type_translations(2, i))

   if (romstype < 0) then
      call error_handler(E_ERR, 'get_translation_table', &
                 'unrecognized roms description: '//trim(type_translations(1, i)), &
                         source, revision, revdate)
   endif

   if (darttype < 0) then
      call error_handler(E_ERR, 'get_translation_table', &
                 'unknown DART Observation TYPE: '//trim(type_translations(2, i)), &
                         source, revision, revdate)
   endif

   ! Now we have to fill the dart_type_integers with the DART type in the same order
   ! find the index of the romstype so we can stuff the darttype in the right slot

   DUMBLOOP : do j = 1,roms_prov_count

      if (romstype == roms_prov_integers(j)) then
          dart_type_integers(j) = darttype
          exit DUMBLOOP
      endif

   enddo DUMBLOOP

   entries = entries + 1
enddo parseloop

if (type_translations(2, entries+1) /= 'NULL') then
   call error_handler(E_ERR, 'get_translation_table', 'roms/dart kinds table must be specified in pairs', &
                      source, revision, revdate)
endif

! Print a report of the obs_provenance strings that match a DART type
! i.e. the supported observations

if (verbose > 0) then
   call error_handler(E_MSG, '', '')
   call error_handler(E_MSG, '', 'DART types that map to ROMS provenance values:')

   do i = 1,roms_prov_count
      if (dart_type_integers(i) > 0) then
   
         write(string1,'(''DART '',i4,1x,A34,'' equates to ROMS '',i8,1x,A)') &
            dart_type_integers(i), trim(get_name_for_type_of_obs(dart_type_integers(i))), &
            roms_prov_integers(i), trim(roms_prov_strings(i))
    
         call error_handler(E_MSG, '', string1)
      endif
   enddo
endif

! Print a report that indicates that there were obs_provenance strings not matched by
! DART namelist 'type_translations' entries.

call error_handler(E_MSG, '', '')
call error_handler(E_MSG, '', 'ROMS provenance values that do not map to a DART type:')

if (all(dart_type_integers(1:roms_prov_count) > 0)) then
   call error_handler(E_MSG, '','none - all provenance values have a DART counterpart - good.')
   call error_handler(E_MSG, '', '')
else
   do i = 1,roms_prov_count
      if ( dart_type_integers(i) < 1 )  then
         write(string1,*)'..  WARNING: provenance integer,string ',roms_prov_integers(i), ',"'//trim(roms_prov_strings(i))//'"'
         write(string2,*)'WARNING: has no matching entry in the input.nml:&convert_roms_obs_nml:type_translations table.'
         write(string3,*)'WARNING: These observations WILL BE IGNORED.'
         call error_handler(E_MSG, '', string1, text2=string2, text3=string3)
         call error_handler(E_MSG, '', '')
      endif
   enddo

   write(string1,*)'At least one observation provenance does not correspond to a DART type.'
   write(string2,*)'Specify the provenance in the input.nml:&convert_roms_obs_nml:type_translations table -or- '
   write(string3,*)'if you want to continue, set input.nml:&utilities_nml:termlevel = 2'
   call error_handler(E_WARN, 'convert_roms_obs:get_translation_table', string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif

end subroutine get_translation_table


!-----------------------------------------------------------------------

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
