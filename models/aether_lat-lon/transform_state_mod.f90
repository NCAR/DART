! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module transform_state_mod

!-----------------------------------------------------------------------
!
! Routines used by aether_to_dart and dart_to_aether
!
!-----------------------------------------------------------------------

use types_mod, only : &
    r4, r8, MISSING_R4, MISSING_R8, vtablenamelength, MISSING_I, RAD2DEG

use time_manager_mod, only : &
    time_type, set_calendar_type, set_time, get_time, set_date, &
    print_date, print_time


use utilities_mod, only : &
    open_file, close_file, file_exist, &
    error_handler, E_ERR, E_MSG, E_WARN, &
    nmlfileunit, do_output, do_nml_file, do_nml_term,  &
    find_namelist_in_file, check_namelist_read

use netcdf_utilities_mod, only : &
    nc_open_file_readonly, nc_open_file_readwrite, nc_create_file, &
    nc_get_dimension_size, nc_get_variable, &
    nc_close_file

implicit none
private

public :: static_init_blocks, &
          state_time, &
          block_file_name, open_block_file, aether_name_to_dart, &
          nblocks_lon,  nblocks_lat,  nblocks_lev, &
          lons, lats, levs, &
          nlon, nlat, nlev, &
          nx_per_block, ny_per_block, nz_per_block, nghost, &
          variables, VT_ORIGININDX, VT_VARNAMEINDX, &
          nvar, nvar_neutral, nvar_ion, &
          aether_restart_dirname, &
          purge_chars, debug

character(len=256), parameter :: source   = 'aether_lat-lon/transform_state_mod.f90'

logical :: module_initialized = .false.

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character(len=256) :: aether_restart_dirname    = '.'
! An ensemble of file names is created using this root and $member in it,

integer, parameter              :: MAX_STATE_VARIABLES     = 100
integer, parameter              :: NUM_STATE_TABLE_COLUMNS = 2
character(len=vtablenamelength) :: variables(NUM_STATE_TABLE_COLUMNS,MAX_STATE_VARIABLES) = ' ' 

! number of blocks along each dim
integer  :: nblocks_lon=MISSING_I, nblocks_lat=MISSING_I, nblocks_lev=MISSING_I
! These are not used in DA, and lon_start is used only for 1D modeling
! real(r8) :: lat_start  =MISSING_I, lat_end    =MISSING_I, lon_start=MISSING_I

integer  :: debug = 0

namelist /transform_state_nml/ aether_restart_dirname, variables, debug, &
                               nblocks_lon, nblocks_lat, nblocks_lev

!-----------------------------------------------------------------------
! Dimensions

! To be assigned get_grid_from_blocks (aether_to_dart, dart_to_aether).
integer                :: nlev, nlat, nlon
real(r8), allocatable  :: levs(:), lats(:), lons(:)

! Aether block parameters (nblocks_{lon,lat,lev} are read from a namelist)
integer  :: nx_per_block, ny_per_block, nz_per_block

integer, parameter :: nghost = 2   ! number of ghost cells on all edges

!-----------------------------------------------------------------------
! Codes for interpreting the NUM_STATE_TABLE_COLUMNS of the variables table
! VT_ORIGININDX is used differently from the usual domains context.
integer, parameter :: VT_VARNAMEINDX  = 1 ! ... variable name
integer, parameter :: VT_ORIGININDX   = 2 ! file of origin

!-----------------------------------------------------------------------
! Day 0 in Aether's calendar is (+/1 a day) -4710/11/24 0 UTC
! integer               :: aether_ref_day = 2451545  ! cJULIAN2000 in Aether = day of date 2000/01/01.
character(len=32)     :: calendar = 'GREGORIAN'

! But what we care about is the ref time for the times in the files, which is 1965-1-1 00:00
integer, dimension(:) :: aether_ref_date(5) = (/1965,1,1,0,0/)  ! y,mo,d,h,m (secs assumed 0)
type(time_type)       :: aether_ref_time, state_time 
integer               :: aether_ref_ndays, aether_ref_nsecs

!-----------------------------------------------------------------------
! to be assigned in the verify_variables subroutine
integer  :: nvar, nvar_neutral, nvar_ion

!-----------------------------------------------------------------------
character(len=512) :: error_string_1, error_string_2

contains

!-----------------------------------------------------------------------
! Like static_init_model, but for aether_to_dart and dart_to_aether.
!   Read the namelist, 
!   parse the 'variables' table,
!   get the Aether grid information
!   convert the Aether time into a DART time.

subroutine static_init_blocks(member)

integer, intent(in) :: member

character(len=128)  :: aether_filename
integer             :: iunit, io

character(len=*), parameter :: routine = 'static_init_blocks'

if (module_initialized) return ! only need to do this once

! This prevents subroutines called from here from calling static_init_mod.
module_initialized = .true.

!------------------
! Read the namelist

call find_namelist_in_file("input.nml", 'transform_state_nml', iunit)
read(iunit, nml = transform_state_nml, iostat = io)
! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=transform_state_nml)
if (do_nml_term()) write(     *     , nml=transform_state_nml)
call check_namelist_read(iunit, io, 'transform_state_nml') ! closes, too.


! error-check, convert namelist input 'variables' to global variables.
call verify_variables(variables)

! Aether uses Julian time internally, andor a Julian calendar
! (days from the start of the calendar), depending on the context)
call set_calendar_type( calendar )   

!--------------------------------
! 1) get grid dimensions
! 2) allocate space for the grids
! 3) read them from the block restart files, could be stretched ...
! Opens and closes the grid block file, but not the filter netcdf file.
call get_grid_from_blocks()

if( debug  > 0 ) then
    write(error_string_1,'(A,3I5)') 'grid dims are ', nlon, nlat, nlev
    call error_handler(E_MSG, routine, error_string_1, source)
endif

! Convert the Aether reference date (not calendar day = 0 date)
! to the days and seconds of the calendar set in model_mod_nml.
aether_ref_time = set_date(aether_ref_date(1), aether_ref_date(2), aether_ref_date(3), &
                     aether_ref_date(4), aether_ref_date(5))
call get_time(aether_ref_time, aether_ref_nsecs, aether_ref_ndays)

! Get the model time from a restart file.
aether_filename = block_file_name(variables(VT_ORIGININDX,1), member, 0)
state_time = read_aether_time(trim(aether_restart_dirname)//'/'//trim(aether_filename))

if ( debug > 0 ) then
  write(error_string_1,'("grid: nlon, nlat, nlev =",3(1x,i5))') nlon, nlat, nlev
  call error_handler(E_MSG, routine, error_string_1, source)
endif
    
end subroutine static_init_blocks

!-----------------------------------------------------------------------
! Parse the table of variables' characteristics.

subroutine verify_variables(variables)

character(len=*), intent(in) :: variables(:,:)

character(len=vtablenamelength) :: varname, rootstr
integer  :: i

character(len=*), parameter  :: routine = 'verify_variables'

nvar = 0
MY_LOOP : do i = 1, size(variables,2)

   varname = variables(VT_VARNAMEINDX,i)
   rootstr = variables(VT_ORIGININDX,i)

   if ( varname == ' ' .and. rootstr == ' ' ) exit MY_LOOP ! Found end of list.

   if ( varname == ' ' .or.  rootstr == ' ' ) then
      error_string_1 = 'variable list not fully specified'
      call error_handler(E_ERR, routine, error_string_1, source)
   endif
  
   if (i > 1) then
   if (variables(VT_ORIGININDX,i-1) == 'ions' .and. rootstr /= 'ions'  ) then
      write(error_string_1,'(A,I1,A)') ' File type (',i, &
           ') in transform_state_nml:variables is out of order or invalid.'
      call error_handler(E_ERR, routine, error_string_1, source)
   endif
   endif

   ! The internal DART routines check if the variable name is valid.

   ! All good to here - fill the output variables

   nvar = nvar + 1
   if (variables(VT_ORIGININDX,i) == 'neutrals') nvar_neutral = nvar_neutral + 1
   if (variables(VT_ORIGININDX,i) == 'ions') nvar_ion = nvar_ion + 1


enddo MY_LOOP

if (nvar == MAX_STATE_VARIABLES) then
   error_string_1 = 'WARNING: you may need to increase "MAX_STATE_VARIABLES"'
   write(error_string_2,'(''you have specified at least '',i4,'' perhaps more.'')') nvar
   call error_handler(E_MSG, routine, error_string_1, source, text2=error_string_2)
endif

end subroutine verify_variables

!-----------------------------------------------------------------------
! ? Will this need to open the grid_{below,corners,down,left} filetypes?
!   This code can handle it; a longer filetype passed in, and no member.
! ? Aether output files?

function block_file_name(filetype, memnum, blocknum)

character(len=*), intent(in)  :: filetype  ! one of {grid,ions,neutrals}
integer,          intent(in)  :: blocknum
integer,          intent(in)  :: memnum
character(len=128)            :: block_file_name

character(len=*), parameter :: routine = 'block_file_name'

block_file_name = trim(filetype)
if (memnum   >= 0) write(block_file_name, '(A,A2,I0.4)') trim(block_file_name), '_m', memnum
if (blocknum >= 0) write(block_file_name, '(A,A2,I0.4)') trim(block_file_name), '_g', blocknum
block_file_name = trim(block_file_name)//'.nc'
if ( debug > 0 ) then
   write(error_string_1,'("filename, memnum, blocknum = ",A,2(1x,i5))') &
        trim(block_file_name), memnum, blocknum
   call error_handler(E_MSG, routine, error_string_1, source)
endif
    
end function block_file_name

!-----------------------------------------------------------------------
! Read block grid values (2D arrays) from a grid NetCDF file.
! Allocate and fill the full-domain 1-D dimension arrays (lon, lat, levs)

! This routine needs:
!
! 1.  A base dirname for the restart files (aether_restart_dirname).
! The filenames have the format 'dirname/{neutrals,ions}_mMMMM_gBBBB.rst'  
! where BBBB is the block number, MMMM is the member number,
! and they have leading 0s.   Blocks start in the
! southwest corner of the lat/lon grid and go east first, 
! then to the west end of the next row north and end in the northeast corner. 
! 
! In the process, the routine will find:
!
! 1. The number of blocks in Lon and Lat (nblocks_lon, nblocks_lat)
!
! 2. The number of lons and lats in a single grid block  (nx_per_block, ny_per_block, nz_per_block)
!
! 3. The overall grid size, {nlon,nlat,nalt} when you've read in all the blocks. 
!
! 4. The number of neutral species (and probably a mapping between
!    the species number and the variable name)  (nvar_neutral)
!
! 5. The number of ion species (ditto - numbers <-> names) (nvar_ion)
!
! In addition to reading in the state data, it fills Longitude, Latitude, and Altitude arrays.  
! This grid is orthogonal and rectangular but can have irregular spacing along
! any of the three dimensions.

subroutine get_grid_from_blocks()

integer               :: nb, offset, ncid, nboff
integer               :: starts(3), ends(3), xcount, ycount, zcount
character(len=256)    :: filename
real(r4), allocatable :: temp(:,:,:)

character(len=*), parameter :: routine = 'get_grid_from_blocks'

! Read the x,y,z  from a NetCDF block file(s),
! in order to calculate the n[xyz]_per_block dimensions. 
! grid_g0000.nc looks like a worthy candidate, but a restart could be used.
write (filename,'(2A)')  trim(aether_restart_dirname),'/grid_g0000.nc'
ncid = nc_open_file_readonly(filename, routine)

! The grid (and restart) file variables have halos, so strip them off
! to get the number of actual data values in each dimension of the block.
nx_per_block = nc_get_dimension_size(ncid, 'x', routine) - (2 * nghost)
ny_per_block = nc_get_dimension_size(ncid, 'y', routine) - (2 * nghost)
nz_per_block = nc_get_dimension_size(ncid, 'z', routine)

nlon = nblocks_lon * nx_per_block
nlat = nblocks_lat * ny_per_block
nlev = nblocks_lev * nz_per_block     

write(error_string_1,'(3(A,I5))') 'nlon = ', nlon, ', nlat = ', nlat, ', nlev = ', nlev
call error_handler(E_MSG, routine, error_string_1, source)

allocate( lons( nlon ))
allocate( lats( nlat ))
allocate( levs( nlev ))

if (debug > 4) then
   write(error_string_1,'(2A)') 'Successfully read Aether grid file:', trim(filename)
   call error_handler(E_MSG, routine, error_string_1, source)
   write(error_string_1,'(A,I5)') '   nx_per_block:', nx_per_block, &
        ' ny_per_block:', ny_per_block, ' nz_per_block:', nz_per_block
   call error_handler(E_MSG, routine, error_string_1, source)
endif

! A temp array large enough to hold any of the 3D
! Lon, Lat or Alt arrays from a block plus ghost cells.
! The restart files have C-indexing (fastest changing dim is the last).
allocate(temp( 1:nz_per_block, &
               1-nghost:ny_per_block+nghost, &
               1-nghost:nx_per_block+nghost))
temp = MISSING_R4 

starts(1) = 1 - nghost
starts(2) = 1 - nghost
starts(3) = 1
ends(1)   = nx_per_block + nghost
ends(2)   = ny_per_block + nghost
ends(3)   = nz_per_block
xcount = nx_per_block + (2 * nghost)
ycount = ny_per_block + (2 * nghost)
zcount = nz_per_block
if ( debug > 0 ) then
   write(error_string_1,'(2(A,3i5),A,3(1X,i5))') &
        'starts = ',starts, 'ends = ',ends, '[xyz]counts = ',xcount,ycount,zcount
   call error_handler(E_MSG, routine, error_string_1, source)
endif

! go across the south-most block row picking up all longitudes
do nb = 1, nblocks_lon

   ! filename is trimmed by passage to open_block_file + "len=*" there.
   filename = trim(aether_restart_dirname)//'/'//block_file_name('grid', -1, nb-1)
   ncid = open_block_file(filename, 'read')

   ! Read 3D array and extract the longitudes of the non-halo data of this block.
   ! The restart files have C-indexing (fastest changing dim is the last),
   ! So invert the dimension bounds.
   call nc_get_variable(ncid, 'Longitude', &
        temp(starts(3):ends(3), starts(2):ends(2), starts(1):ends(1)), &
        context=routine, &
        nc_count=(/ zcount,ycount,xcount /))

   offset = (nx_per_block * (nb - 1))
   lons(offset+1:offset+nx_per_block) = temp(1,1,1:nx_per_block)

   call nc_close_file(ncid)
enddo

! go up west-most block row picking up all latitudes
do nb = 1, nblocks_lat

   ! Aether's block name counter start with 0, but the lat values can come from 
   ! any lon=const column of blocks. 
   nboff = ((nb - 1) * nblocks_lon)
   filename = trim(aether_restart_dirname)//'/'//block_file_name('grid', -1, nboff)
   ncid = open_block_file(filename, 'read')

   call nc_get_variable(ncid, 'Latitude', &
        temp(starts(3):ends(3), starts(2):ends(2), starts(1):ends(1)), &
        context=routine, nc_count=(/zcount,ycount,xcount/))


   offset = (ny_per_block * (nb - 1))
   lats(offset+1:offset+ny_per_block) = temp(1,1:ny_per_block,1)

   call nc_close_file(ncid)
enddo


! this code assumes all columns share the same altitude array,
! so we can read it from the first block.
! if this is not the case, this code has to change.

filename = trim(aether_restart_dirname)//'/'//block_file_name('grid', -1, 0)
ncid = open_block_file(filename, 'read')

temp = MISSING_R8
call nc_get_variable(ncid, 'Altitude', &
     temp(starts(3):ends(3), starts(2):ends(2), starts(1):ends(1)), &
     context=routine, nc_count=(/zcount,ycount,xcount/))

levs(1:nz_per_block) = temp(1:nz_per_block,1,1)

call nc_close_file(ncid)

deallocate(temp)

! convert from radians into degrees
lons = lons * RAD2DEG
lats = lats * RAD2DEG

if (debug > 4) then
   print *, routine, 'All lons ', lons
   print *, routine, 'All lats ', lats
   print *, routine, 'All levs ', levs
endif

if ( debug > 1 ) then ! Check dimension limits
   write(error_string_1,'(A,2F15.4)') 'LON range ', minval(lons), maxval(lons)
   call error_handler(E_MSG, routine, error_string_1, source)
   write(error_string_1,'(A,2F15.4)') 'LAT range ', minval(lats), maxval(lats)
   call error_handler(E_MSG, routine, error_string_1, source)
   write(error_string_1,'(A,2F15.4)') 'ALT range ', minval(levs), maxval(levs)
   call error_handler(E_MSG, routine, error_string_1, source)
endif

end subroutine get_grid_from_blocks

!-----------------------------------------------------------------------
! Read the Aether restart file time and convert to a DART time.

function read_aether_time(filename)
type(time_type)              :: read_aether_time
character(len=*), intent(in) :: filename

integer  :: ncid
integer  :: tsimulation   ! the time read from a restart file; seconds from aether_ref_date.
integer  :: ndays, nsecs

character(len=*), parameter :: routine = 'read_aether_time'

tsimulation = MISSING_I

ncid = open_block_file(filename, 'read')
call nc_get_variable(ncid, 'time', tsimulation, context=routine)
call nc_close_file(ncid, routine, filename)

! Calculate the DART time of the file time.
ndays     = tsimulation / 86400
nsecs     = tsimulation - (ndays * 86400)
! The ref day is not finished, but don't need to subtract 1 because 
! that was accounted for in the integer calculation of ndays.
ndays     = aether_ref_ndays + ndays
read_aether_time = set_time(nsecs, ndays)

if (do_output()) &
    call print_time(read_aether_time, routine//': time in restart file '//filename)
if (do_output()) &
    call print_date(read_aether_time, routine//': date in restart file '//filename)

if (debug > 8) then
   write(error_string_1,'(A,I5)')'tsimulation ', tsimulation
   call error_handler(E_MSG, routine, error_string_1, source)
   write(error_string_1,'(A,I5)')'ndays       ', ndays
   call error_handler(E_MSG, routine, error_string_1, source)
   write(error_string_1,'(A,I5)')'nsecs       ', nsecs
   call error_handler(E_MSG, routine, error_string_1, source)

   call print_date(aether_ref_time, routine//':model base date')
   call print_time(aether_ref_time, routine//':model base time')
endif

end function read_aether_time

!-----------------------------------------------------------------------
! Convert Aether's non-CF-compliant names into CF-compliant names for filter.
! For the ions, it moves the name of the ion from the end of the variable names
! to the beginning.

function aether_name_to_dart(varname)

character(len=vtablenamelength), intent(in) :: varname

character(len=vtablenamelength) :: aether_name_to_dart, aether
character(len=64)               :: parts(8), var_root
integer                         :: char_num, first, i_parts, aether_len, end_str

aether     = trim(varname)
aether_len = len_trim(varname)
parts = ''

! Look for the last ' '.  The characters after that are the species.
! If there's no ' ', the whole string is the species.
char_num = 0
char_num = scan(trim(aether),' ', back=.true.)
var_root = aether(char_num+1:aether_len)
! purge_chars removes unwanted [()\] 
parts(1) = purge_chars( trim(var_root),')(\', plus_minus=.true.)
end_str  = char_num

! Tranform remaining pieces of varname into DART versions.
char_num = MISSING_I
first = 1
i_parts = 2
P_LOOP: do
   ! This returns the position of the first blank *within the substring* passed in.
   char_num = scan(aether(first:end_str),' ', back=.false.)
   if (char_num > 0 .and. first < aether_len) then
      parts(i_parts) = purge_chars(aether(first:first+char_num-1), '.)(\', plus_minus=.true.)

      first   = first + char_num 
      i_parts = i_parts + 1
   else
      exit P_LOOP
   endif
enddo P_LOOP

! Construct the DART field name from the parts
aether_name_to_dart = trim(parts(1))
i_parts = 2
Build : do
   if (trim(parts(i_parts)) /= '') then
      aether_name_to_dart = trim(aether_name_to_dart)//'_'//trim(parts(i_parts))
      i_parts = i_parts + 1
   else
      exit Build
   endif
enddo Build
   
end function aether_name_to_dart
   
!-----------------------------------------------------------------------
! Replace undesirable characters with better.
   
function purge_chars(ugly_string, chars, plus_minus)
   
character (len=*), intent(in) :: ugly_string, chars
logical,           intent(in) :: plus_minus
character (len=64)            :: purge_chars

character (len=256) :: temp_str

integer :: char_num, end_str, pm_num

! Trim is not needed here
temp_str = ugly_string
end_str  = len_trim(temp_str)
char_num = MISSING_I
Squeeze : do 
   ! Returns 0 if chars are not found
   char_num = scan(temp_str, chars)
   ! Need to change it to a char that won't be found by scan in the next iteration,
   ! and can be easily removed.
   if (char_num > 0) then
      ! Squeeze out the character
      temp_str(char_num:end_str-1) = temp_str(char_num+1:end_str)
      temp_str(end_str:end_str) = ''
!       temp_str(char_num:char_num) = ' '
   else
      exit Squeeze
   endif
enddo Squeeze

! Replace + and - with pos and neg.  Assume there's only 1.
temp_str = trim(adjustl(temp_str))
end_str  = len_trim(temp_str)
pm_num   = scan(trim(temp_str),'+-', back=.false.)
if (pm_num == 0 .or. .not. plus_minus) then
   purge_chars = trim(temp_str)
else
   if (temp_str(pm_num:pm_num) == '+') then
      purge_chars = temp_str(1:pm_num-1)//'pos'
   else if (temp_str(pm_num:pm_num) == '-') then
      purge_chars = temp_str(1:pm_num-1)//'neg'
   endif
   if (pm_num + 1 <= end_str) &
      purge_chars = trim(purge_chars)//temp_str(pm_num+1:end_str)
endif
   
end function purge_chars

!-----------------------------------------------------------------------
! Open an Aether restart block file (neutral, ion, ...?)

function open_block_file(filename, rw)

! filename is trimmed by this definition
character(len=*), intent(in) :: filename
character(len=*), intent(in) :: rw   ! 'read' or 'readwrite'
integer                      :: open_block_file

character(len=*), parameter :: routine = 'open_block_file'

if ( .not. file_exist(filename) ) then
   write(error_string_1,'(4A)') 'cannot open file ', filename,' for ', rw
   call error_handler(E_ERR, routine, error_string_1, source)
endif

if (debug > 0) then
   write(error_string_1,'(4A)') 'Opening file ', trim(filename), ' for ', rw
   call error_handler(E_MSG, routine, error_string_1, source)
end if


if (rw == 'read') then
   open_block_file = nc_open_file_readonly(filename, routine)
else if (rw == 'readwrite') then
   open_block_file = nc_open_file_readwrite(filename, routine)
else
   error_string_1 = ': must be called with rw={read,readwrite}, not '//rw
   call error_handler(E_ERR, routine, error_string_1, source)
endif


if (debug > 80) then
   write(error_string_1,'(4A)') 'Returned file descriptor is ', open_block_file
   call error_handler(E_MSG, routine, error_string_1, source)
end if
    
end function open_block_file

end module transform_state_mod
