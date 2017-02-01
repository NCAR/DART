! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program read_geo

implicit none

! Given a station id, find the lat, lon and elevation from the geo file
! provided by the Oklahoma Mesonent website
!
! the geoinfo.csv file has data in the following order:
!
! stnm - station number
! stid - CHAR len=4
! name - CHAR of variable size
! city - CHAR of variable size
! rang - distance of station from city, float
! cdir - CHAR of variable size
! cnty - CHAR of variable size
! nlat - latitude, float
! nlon - longitude, float
! elev - elevation, float
! + other stuff we don't need
!

real :: nlat, nlon
integer :: elev
character(len=4) :: stid
character(len=4) :: search_stid

! temp vars
integer :: stnm
character(len=15) :: name, city, cdir, cnty
real :: rang
character(len=330) :: line
integer :: n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, nend
logical :: stn_found
! passed variable
   search_stid = 'BBOW'
  stn_found = .false.
  open(1, file='geoinfo.csv', status='old')

 read(1,*) ! header

  do while (1 .eq. 1)
! read through each line and look for a match 
 read(1,'(A330)',ERR=200) line
! write (*,*) line

! Determine locations of the first 10 delimiters
n1  = index(line, ',')
nend = len_trim(line)
n2  = n1 + index(line(n1+1:nend), ',')
n3  = n2 + index(line(n2+1:nend), ',')
n4  = n3 + index(line(n3+1:nend), ',')
n5  = n4 + index(line(n4+1:nend), ',')
n6  = n5 + index(line(n5+1:nend), ',')
n7  = n6 + index(line(n6+1:nend), ',')
n8  = n7 + index(line(n7+1:nend), ',')
n9  = n8 + index(line(n8+1:nend), ',')
n10 = n9 + index(line(n9+1:nend), ',')

!write (*,*) 'index values ',n1, n2, n3
!
 read (line(1:n1-1),'(I3)') stnm
 read (line(n1+1:n2-1),'(A4)') stid
 read (line(n2+1:n3-1),'(A15)') name
 read (line(n3+1:n4-1),'(A15)') city
 read (line(n4+1:n5-1),'(f6.3)') rang
 read (line(n5+1:n6-1),'(A3)') cdir
 read (line(n6+1:n7-1),'(A15)') cnty
 read (line(n7+1:n8-1),'(f10.7)') nlat
 read (line(n8+1:n9-1),'(f10.7)') nlon
 read (line(n9+1:n10-1),'(I4)') elev
! read (line,'(I3, A4, A15, A15, F6.3, A3, A15, &
! read (1,'(I3, A4, A15, A15, F6.3, A3, A15, &
!       F10.7, F10.7, F8.3)') stnm,stid,name,city &
!       ,rang,cdir,cnty,nlat,nlon,elev
! write(*,*) 'found ',stnm, nlat, nlon, elev
  if (search_stid .eq. stid) then
    write(*,*) 'station found ',stid, nlat, nlon, elev
    stn_found = .true.
    goto 200
  end if

  end do
200  continue
  if (stn_found) then
    write(*,*) ' done' 
  else
    write(*,*) ' station ',stid,' not found'
  end if


end program read_geo

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
