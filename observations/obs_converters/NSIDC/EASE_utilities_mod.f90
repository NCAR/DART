
module EASE_utilities_mod

implicit none
private

public :: ezlh_convert, ezlh_inverse
public :: get_grid_dims
public :: deconstruct_filename
public :: read_ease_Tb
public :: read_ease_TIM
public :: EASE_MISSING

! This routine came from http://nsidc.org/data/ease/tools.html#geo_data_files
! as ezlconv.f ... and was converted to F90 and put in a module by TJH.

! RE_km    ... radius of the earth (km), authalic sphere based on International datum 
! CELL_km  ... nominal cell size in kilometers
! COS_PHI1 ... scale factor for standard paralles at +/-30.00 degrees

real, parameter :: RE_km    = 6371.228
real, parameter :: CELL_km  = 25.067525
real, parameter :: COS_PHI1 = 0.866025403
real, parameter :: PI       = 3.141592653589793

integer, parameter :: l_nrows = 721
integer, parameter :: l_ncols = 721
integer, parameter :: h_nrows = 586
integer, parameter :: h_ncols = 1383
integer, parameter :: EASE_MISSING = -32768

contains

!==========================================================================
! ezlhconv.for - FORTRAN routines for conversion of azimuthal 
!		equal area and equal area cylindrical grid coordinates
!
!	30-Jan.-1992 H.Maybee
!	20-Mar-1992 Ken Knowles  303-492-0644  knowles@kryos.colorado.edu
!       16-Dec-1993 MJ Brodzik   303-492-8263  brodzik@jokull.colorado.edu
!                   Copied from nsmconv.f, changed resolutions from 
!                   40-20-10 km to 25-12.5 km
!       21-Dec-1993 MJ Brodzik   303-492-8263  brodzik@jokull.colorado.edu
!                   Fixed sign of Southern latitudes in ease_inverse.
!	12-Sep-1994 David Hoogstrate 303-492-4116 hoogstra@jokull.colorado.edu
!		    Changed grid cell size. Changed "c","f" to "l","h"
!	25-Oct-1994 David Hoogstrate 303-492-4116 hoogstra@jokull.colorado.edu
!		    Changed row size from 587 to 586 for Mercator projection
!		    Changed function names to "ezlh-.."
!$Log: ezlhconv.f,v $
!Revision 1.3  1994/11/01 23:40:43  brodzik
!Replaced all references to 'ease' with 'ezlh'
!
!==========================================================================


!--------------------------------------------------------------------------
function ezlh_convert (grid, lat, lon, r, s)
character(len=*), intent(in)  :: grid
real,             intent(in)  :: lat, lon
real,             intent(out) :: r, s
integer                       :: ezlh_convert
!
!	convert geographic coordinates (spherical earth) to 
!	azimuthal equal area or equal area cylindrical grid coordinates
!
!	status = ezlh_convert (grid, lat, lon, r, s)
!
!	input : grid - projection name '[NSM][lh]'
!               where l = "low"  = 25km resolution
!                     h = "high" = 12.5km resolution
!		lat, lon - geo. coords. (decimal degrees)
!
!	output: r, s - column, row coordinates
!
!	result: status = 0 indicates normal successful completion
!			-1 indicates error status (point not on grid)
!
!--------------------------------------------------------------------------
integer :: cols, rows, scale
real    :: Rg, phi, lam, rho, r0, s0

ezlh_convert = -1

if ((grid(1:1) == 'N') .or. (grid(1:1) == 'S')) then
   cols = l_ncols
   rows = l_nrows
else if (grid(1:1) == 'M') then
   cols = h_ncols
   rows = h_nrows
else
   print *, 'ezlh_convert: unknown projection: ', grid
   return
endif

if      ((grid(2:2) == 'l') .or. (grid(2:2) == 'L')) then
   scale = 1
else if ((grid(2:2) == 'h') .or. (grid(2:2) == 'H')) then
   scale = 2
else
   print *, 'ezlh_convert: unknown resolution: ', grid
   return
endif

Rg = scale * RE_km/CELL_km

! r0,s0 are defined such that cells at all scales
! have coincident center points

r0  = (cols-1)/2.0 * scale
s0  = (rows-1)/2.0 * scale
phi = lat * PI/180.0
lam = lon * PI/180.0

if (grid(1:1) == 'N') then ! Northern
   rho = 2.0 * Rg * sin(PI/4.0 - phi/2.0)
   r   = r0 + rho * sin(lam)
   s   = s0 + rho * cos(lam)
else if (grid(1:1) == 'S') then ! Southern
   rho = 2.0 * Rg * cos(PI/4.0 - phi/2.0)
   r   = r0 + rho * sin(lam)
   s   = s0 - rho * cos(lam)
else if (grid(1:1) == 'M') then ! Global
   r = r0 + Rg * lam * COS_PHI1
   s = s0 - Rg * sin(phi) / COS_PHI1
endif

ezlh_convert = 0
return
end function ezlh_convert



!--------------------------------------------------------------------------
function ezlh_inverse (grid, r, s, lat, lon)
character(len=*), intent(in)  :: grid
real,             intent(in)  :: r, s
real,             intent(out) :: lat, lon
integer                       :: ezlh_inverse
!	convert azimuthal equal area or equal area cylindrical 
!	grid coordinates to geographic coordinates (spherical earth)
!
!	status = ezlh_inverse (grid, r, s, lat, lon)
!
!	input : grid - projection name '[NSM][lh]'
!               where l = "low"  = 25km resolution
!                     h = "high" = 12.5km resolution
!		r, s - grid column and row coordinates
!
!	output: lat, lon - geo. coords. (decimal degrees)
!
!	result: status = 0 indicates normal successful completion
!			-1 indicates error status (point not on grid)
!
!--------------------------------------------------------------------------
integer :: cols, rows, scale
real    :: Rg, phi, lam, rho, r0, s0, c
real    :: gamma, beta, epsilon, x, y
real    :: sinphi1, cosphi1

ezlh_inverse = -1

if ((grid(1:1) == 'N') .or. (grid(1:1) == 'S')) then
   cols = l_ncols
   rows = l_nrows
else if (grid(1:1) == 'M') then
   cols = h_ncols
   rows = h_nrows
else
   print *, 'ezlh_inverse: unknown projection: ', grid
   return
endif

if      ((grid(2:2) == 'l') .or. (grid(2:2) == 'L')) then
   scale = 1
else if ((grid(2:2) == 'h') .or. (grid(2:2) == 'H')) then
   scale = 2
else
   print *, 'ezlh_inverse: unknown resolution: ', grid
   return
endif

Rg = scale * RE_km/CELL_km

r0 = (cols-1)/2.0 * scale
s0 = (rows-1)/2.0 * scale
x  =   r - r0
y  = -(s - s0)

if ((grid(1:1) == 'N').or.(grid(1:1) == 'S')) then 
   rho = sqrt(x*x + y*y)
   if (rho == 0.0) then
      if (grid(1:1) == 'N') lat =  90.0 
      if (grid(1:1) == 'S') lat = -90.0 
      lon = 0.0
   else
      if (grid(1:1) == 'N') then
         sinphi1 = sin(PI/2.0)
         cosphi1 = cos(PI/2.0)
         if (y == 0.0) then
            if (r <= r0) lam = -PI/2.0
            if (r  > r0) lam =  PI/2.0
         else
            lam = atan2(x,-y)
         endif
      else if (grid(1:1) == 'S') then
         sinphi1 = sin(-PI/2.0)
         cosphi1 = cos(-PI/2.0)
         if (y == 0.0) then
            if (r <= r0) lam = -PI/2.0
            if (r  > r0) lam =  PI/2.0
         else
            lam = atan2(x,y)
         endif
      endif
      gamma = rho/(2.0 * Rg)
      if (abs(gamma) > 1.0) return
      c    = 2.0 * asin(gamma)
      beta = cos(c) * sinphi1 + y * sin(c) * (cosphi1/rho)
      if (abs(beta) > 1.0) return
      phi = asin(beta)
      lat = phi * 180.0 / PI
      lon = lam * 180.0 / PI
   endif

else if (grid(1:1) == 'M') then

   ! allow .5 cell tolerance in arcsin function
   ! so that grid coordinates which are less than .5 cells
   ! above 90.00N or below 90.00S are given a lat of 90.00

   epsilon = 1.0 + 0.5/Rg
   beta = y*COS_PHI1/Rg
   if (abs(beta) > epsilon) return
   if (beta <= -1.0) then
      phi = -PI/2.0
   else if (beta >= 1.0) then
      phi =  PI/2.0
   else
      phi = asin(beta)
   endif
   lam = x/COS_PHI1/Rg
   lat = phi * 180.0 / PI
   lon = lam * 180.0 / PI

endif

ezlh_inverse = 0
return

end function ezlh_inverse




function deconstruct_filename(filename, gridarea, iyear, idoy, &
                passdir, freq, polarization, is_time_file)
character(len=*), intent(in)    :: filename
character(len=2), intent(out)   :: gridarea
integer,          intent(out)   :: iyear, idoy
character(len=1), intent(out)   :: passdir
real,             intent(out)   :: freq
character(len=1), intent(out)   :: polarization
logical,          intent(out)   :: is_time_file
integer                         :: deconstruct_filename

! http://nsidc.org/data/docs/daac/nsidc0301_amsre_gridded_tb.gd.html#namingconvention
! EASE_utilities_modSE-Grid brightness temperature data files are named 
! according to the following convention and as described in Table 1:
!
! ID2rx-AMSRE-aayyyydddp.vnn.ccc

! Where:

! Table 1. EASE-Grid Brightness Temperature Data Files' Naming Convention
! Variable	Description
! ID2           Inverse Distance Squared
! rx            Resolution number of swath input data (r1, r3)
! AMSRE         Identifies this a file containing AMSR-E data
! aa            Area of coverage (NL = north, SL = south, ML = global)
! yyyy          Four-digit year
! ddd           Three-digit day of year
! p             Pass direction (A = ascending, D = descending)
! vnn           Data version number (for example, v01, v02)
! ccc           AMSR-E channel indicator: numeric frequency 
!               (06, 10, 18, 23, 36, or 89) followed by polarization (H or V)
!
! For example, ID2r3-AMSRE-SL2005135D.v01.89V
!
! The EASE-Grid time files follow the same naming convention, except they 
! end in the extension .TIM instead of the channel indicator.
!
! For example, ID2r3-AMSRE-SL2005135D.v01.TIM

! character(len=2), intent(out) :: gridarea
! integer,          intent(out) :: iyear, idoy
! character(len=1), intent(out) :: passdir
! character(len=1), intent(out) :: polarization
! integer                       :: deconstruct_filename

integer :: i, ipos, nchars, iocode, channel
character(len=512) :: myfile
integer, parameter :: nfreq = 6
real, parameter, dimension(nfreq) :: actualF = (/ 6.9, 10.7, 18.7, 23.8, 36.5, 89.0 /)

deconstruct_filename = 0

! set failed values
gridarea     = 'xx'
iyear        = 0
idoy         = 0
passdir      = 'x'
freq         = 0.0
polarization = 'x'
is_time_file = .false.

! Find the last occurrence of AMSRE, can parse from there.
myfile = adjustl(filename)
nchars = len_trim(myfile)
ipos   = -1
LocateAMSRE : do i = (nchars-4),1,-1
   if (myfile(i:i+4) == 'AMSRE') then
      ipos = i
      exit LocateAMSRE
   endif
enddo LocateAMSRE

if (ipos < 0) then
   write(*,*)'Unable to parse filename [',trim(myfile),']'
   deconstruct_filename = -99
   return
endif

if (myfile(nchars-2:nchars) == 'TIM') then
   is_time_file = .true.
else
   read(myfile(ipos:nchars),200,iostat=iocode)gridarea,iyear,idoy,passdir,channel,polarization
endif

deconstruct_filename = iocode

! Check for failure and exit early 
if (deconstruct_filename /= 0 ) return

! Relate the integer frequencies to the actual frequencies.
! The actual frequencies are:
! (6.9 GHz, 10.7 GHz, 18.7 GHz, 23.8 GHz, 36.5 GHz, and 89.0 GHz)
! The 'channel' values: 06, 10, 18, 23, 36, 89 

freq = 0.0
FREQUENCY: do i = 1,nfreq
   if ( abs(real(channel) - actualF(i)) < 1.0 ) then
      freq = actualF(i)
      exit FREQUENCY
   endif
enddo FREQUENCY

if (freq == 0.0) then
   deconstruct_filename = -99
   return
endif

! Sanity check
if ( 1 == 0 ) then
   write(*,*)'The salient part of the filename is ',myfile(ipos:nchars)
   write(*,*)'gridarea     is ',gridarea
   write(*,*)'iyear        is ',iyear
   write(*,*)'idoy         is ',idoy
   write(*,*)'passdir      is ',passdir
   write(*,*)'channel      is ',channel
   write(*,*)'frequency    is ',freq
   write(*,*)'polarization is ',polarization
   write(*,*)'time file    is ',is_time_file
   write(*,*)
endif

return

200 format (6x,a2,i4,i3,a1,5x,i2,a1)

! ID2r3-AMSRE-NL2011001D.v03.89H
! ID2r3-AMSRE-NL2011001D.v03.89V
! ID2r3-AMSRE-NL2011001D.v03.TIM

end function deconstruct_filename



function read_ease_Tb(filename, iunit, Tb)
! The EASE Tb files are (binary) flat files.
! The brightness temperatures are in unsigned 16 bit integers
! representing tenths of a degree Kelvin
! 0 is a MISSING value.

character(len=*),        intent(in)  :: filename
integer,                 intent(in)  :: iunit
integer, dimension(:,:), intent(out) :: Tb
integer                              :: read_ease_Tb

integer, parameter :: i2 = SELECTED_INT_KIND(4) ! +/- 32767

integer(i2), allocatable, dimension(:,:) :: datmat
integer :: irow, icol, nrows, ncols, iocode

nrows = size(Tb,1)
ncols = size(Tb,2)
allocate(datmat(nrows,ncols))

open(unit=iunit, file=trim(filename), access='direct', &
          form='unformatted', recl=2*nrows*ncols, iostat=iocode)

if (iocode /= 0) then
   write(*,*)'read_ease_Tb failed to open [',trim(filename)//']'
   stop
endif

read(iunit,rec=1,iostat=iocode) datmat
if (iocode /= 0) then
   write(*,*)'read_ease_Tb failed to read [',trim(filename)//']'
   stop
endif

close(iunit)

! Convert unsigned 16 bit integers to signed 32 bit integers
! unsigned 16 bit [0,65535]  2^16 = 65536
do irow = 1,nrows
do icol = 1,ncols
   Tb(irow,icol) = IAND(INT(datmat(irow,icol)),65535)
enddo
enddo

!--------------------------------------------------------------------------------
! Sanity check block ... 
! write it out as an unformatted write so I can read it in matlab
! PASSED 8 Nov 2013 ... see matlab/Check_EASE_data.m
!
! open(unit=iunit, file="sanity_check.ieee", form='unformatted')
! write(iunit) nrows,ncols
! write(iunit) Tb
! close(iunit)
! write(*,*)'i2 max possible is ',HUGE(datmat)
! if ( HUGE(datmat) /= 32767 ) then
!    write(*,*)'ERROR ... run matlab/CheckEASE_data.m to see if Tb read correctly.'
!    write(*,*)'ERROR ... run matlab/CheckEASE_data.m to see if Tb read correctly.'
!    stop
! endif
!--------------------------------------------------------------------------------

read_ease_Tb = 0 ! successful finish

end function read_ease_Tb




function read_ease_TIM(filename, iunit, time_matrix)
! Data are 2-byte signed integers, little-endian byte-order, 
! indicating time of data acquisition as minutes since midnight 
! (0:00 UTC) of the date of the enclosing file. The values 
! in the time files range from 0 to 1440, with the 
! value -32768 indicating missing data.
! See more at
! http://nsidc.org/data/docs/daac/nsidc0301_amsre_gridded_tb.gd.html

character(len=*),        intent(in)  :: filename
integer,                 intent(in)  :: iunit
integer, dimension(:,:), intent(out) :: time_matrix
integer                              :: read_ease_TIM

integer, parameter :: i2 = SELECTED_INT_KIND(4)

integer(i2), allocatable, dimension(:,:) :: datmat
integer :: nrows, ncols, iocode

nrows = size(time_matrix,1)
ncols = size(time_matrix,2)
allocate(datmat(nrows,ncols))

open(unit=iunit, file=trim(filename), access='direct', &
          form='unformatted', recl=2*nrows*ncols, iostat=iocode)

if (iocode /= 0) then
   write(*,*)'read_ease_TIM failed to open [',trim(filename)//']'
   stop
endif

read(iunit,rec=1,iostat=iocode) datmat
if (iocode /= 0) then
   write(*,*)'read_ease_TIM failed to read [',trim(filename)//']'
   stop
endif

close(iunit)

! preserve to output so I can mess around with datmat to check
! for endian issues. If checks pass, 'time_matrix' will be returned.

time_matrix = datmat

! Put in rudimentary check to see if the endian is correct
! by checking values/bounds, etc. Replace any missing values
! with a positive value and then check for ANY negative values
! or any values greater than 1440. Not bulletproof, but pretty good.

where (datmat == EASE_MISSING) datmat = 0

if (any(datmat < 0)) then
   write(*,*)'ERROR: read_ease_TIM read [',trim(filename)//']'
   write(*,*)'ERROR:   encountered out-of-range values (negative)'
   write(*,*)'ERROR:   most likely an endian problem.' 
   write(*,*)'ERROR:   The files are little-endian.'
   stop
endif

if (any(datmat > 1440)) then
   write(*,*)'ERROR: read_ease_TIM read [',trim(filename)//']'
   write(*,*)'ERROR:   encountered out-of-range values ( > 1440 )'
   write(*,*)'ERROR:   most likely an endian problem.' 
   write(*,*)'ERROR:   The files are little-endian.'
   stop
endif

read_ease_TIM = 0 ! successful finish

end function read_ease_TIM



function get_grid_dims( grid )
character(len=*), intent(in) :: grid
integer, dimension(2)        :: get_grid_dims

integer :: cols, rows

get_grid_dims(:) = -1  ! a bad error code

if ((grid(1:1) == 'N') .or. (grid(1:1) == 'S')) then
   cols = l_ncols
   rows = l_nrows
else if (grid(1:1) == 'M') then
   cols = h_ncols
   rows = h_nrows
else
   print *, 'get_grid_dims: unknown projection: ', grid
   return
endif

get_grid_dims(1) = rows
get_grid_dims(2) = cols

end function get_grid_dims


end module



