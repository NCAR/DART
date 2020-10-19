! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!     A short main program to test reading AMSU brightness temperature.
!     March, 2012
!     Feng Ding

program L1_AMSUA_to_netcdf

use amsua_support_mod

! INCLUDE 'amsua_bt_typ.inc'
! INCLUDE 'amsua_bt_struct.inc'

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'L1_AMSUA_to_netcdf.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

! ----------------------------------------------------------------------
! Declare namelist parameters
! ----------------------------------------------------------------------

character(len=256) :: file_name  = ''
character(len=256) :: outputfile = ''
integer            :: track  = 1 ! 0-based index along track
integer            :: xtrack = 0 ! 0-based index across-track

namelist /L1_AMSUA_to_netcdf_nml/ file_name, outputfile, &
                               track, xtrack

TYPE(amsua_bt_gran_t) amsua_bt_gran

integer :: chan    !0-based channel index.

!----------------------------------------------------------------------
! Read the namelist
!----------------------------------------------------------------------

call find_namelist_in_file('input.nml', 'convert_airs_L2_nml', iunit)
read(iunit, nml = convert_airs_L2_nml, iostat = io)
call check_namelist_read(iunit, io, 'convert_airs_L2_nml')

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=convert_airs_L2_nml)
if (do_nml_term()) write(    *      , nml=convert_airs_L2_nml)


if (iargc().ne.3) then
    print *, "This code extracts a single profile from a specified"
    print *, " input file to stdout. It requires exactly three "
    print *, "arguments."
    print *, "  1) scan line number [1, 45]"
    print *, "  2) field-of-view number [1, 30]"
    print *, "  3) file name"
    STOP
 end if 

if (track.lt.1.OR.track.gt.45) then
  print *, "Error: along-track scan line number [1, 45]"
  print *, "got ", track
  STOP
endif

if (xtrack.lt.1.OR.xtrack.gt.30) then
  print *, "Error: second argument must be scan line number [1, 30]"
  print *, "got ", xtrack
  STOP
endif

CALL amsua_bt_rdr(file_name, amsua_bt_gran)

! Each AMSU-A scan has 2 "state"s, indicating whether the AMSU-A1 and
! AMSU-A2 instruments were in science mode when the data
! was taken and whether the data was successfully transmitted.

if (amsua_bt_gran%state1(track).ne.0) then
   if (amsua_bt_gran%state1(track).EQ.1) then 
      print *, "Warning, AMSU-A1 state for this profile is SPECIAL"
   else if (amsua_bt_gran%state1(track).EQ.2) then 
      print *, "Warning, AMSU-A1 state for this profile is ERRONEOUS"
   else if (amsua_bt_gran%state1(track).EQ.3) then 
      print *, "Warning, AMSU-A1 state for this profile is MISSING"
   else
      print *, "Warning, AMSU-A1 state for this profile is UNKNOWN"
   endif 
              
   print *, "NOT PROCESS"

endif

if (amsua_bt_gran%state2(track).ne.0) then
   if (amsua_bt_gran%state2(track).EQ.1) then 
      print *, "Warning, AMSU-A2 state for this profile is SPECIAL"
   else if (amsua_bt_gran%state2(track).EQ.2) then 
      print *, "Warning, AMSU-A2 state for this profile is ERRONEOUS"
   else if (amsua_bt_gran%state2(track).EQ.3) then 
      print *, "Warning, AMSU-A2 state for this profile is MISSING"
   else
      print *, "Warning, AMSU-A2 state for this profile is UNKNOWN"
   endif 
              
   print *, "NOT PROCESS"

endif

print *, "# AMSU Brightness Temperatures (Kelvins)"
print *, "# Channels 1-15"
print *, "# -9999 flags bad value"

DO chan = 1, AMSUA_BT_CHANNEL 
   WRITE(*, "(f8.2)") amsua_bt_gran%brightness_temp(chan,xtrack,track)
ENDDO

STOP
end program L1_AMSUA_to_netcdf
