! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program ftest_nml

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! very simple fortran program which reads in an external namelist file.
! if successful, will print a message and exit.

integer :: iunit, errcode

integer :: array_size, array_data(10)
namelist / ftest / array_size, array_data

   print *, "program start"
  
   array_size = -1
   array_data = -99

   iunit = 11
   open (iunit, name="ftest_nml.nml", iostat=errcode)
   if (errcode /= 0) then
       print *, "cannot open namelist file, error = ", errcode
       stop
   endif

   read (iunit, nml=ftest, iostat=errcode)
   if (errcode /= 0) then
       print *, "cannot read namelist file, error = ", errcode
       stop
   endif

   close (iunit, iostat=errcode)

   print *, "array size should be 10, value is ", array_size
   print *, "array contents should be 1-10, values are: ", array_data

   print *, "program end"

end program ftest_nml

