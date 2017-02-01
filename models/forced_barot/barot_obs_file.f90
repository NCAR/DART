! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$
 
program barot_obs_file

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! Generates a regular set, currently hardcoded to 40 by 30,
! of observation file locations for the forced barot model.

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

write(*, *) 1200
do i = 1, 40
   do j = 1, 30
      lon = ((i - 0.5) / 40.0) * 360.0 
! Don't want to have to interpolate lats
      lat = -86.0 + (j / 30.0) * 170.0
      write(*, *) lon, lat, (1e6)**2
   end do
end do
end program barot_obs_file

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
