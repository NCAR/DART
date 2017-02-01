! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$ 

! THIS PROGRAM HAS BEEN OBSOLETED by converting our restart files to
! NetCDF format.  you can assign a value with 'ncap2', one of the NCO
! utility program.  for now i'm going to leave the mkmf and path_names
! files around and the program will just print out a helpful message.
! in a few months/years we should remove this program entirely.
!
! here's an example of using ncap2 to set the T,U and V inf values:
!  ncap2 -s 'T=1.0;U=1.0;V=1.0' wrfinput_d01 prior_inf.nc
!  ncap2 -s 'T=0.6;U=0.6;V=0.6' wrfinput_d01 prior_sd.nc 
!
! this works as long as you have at least version 4.4.2 of the NCO utils.
! some earlier versions change the full 3d arrays into a single scalar.
! if you see this, get a more recent version of the nco tools.
!
! nsc 27jul2016

program fill_inflation_restart

write(*,*) ''
write(*,*) 'This program is OBSOLETE since DART inflation files are now in NetCDF format.' 
write(*,*) 'To fill an initial inflation file use one of the standard NCO utilities like "ncap2" with'
write(*,*) 'a copy of a model restart file to set the initial inflation mean, and a second file'
write(*,*) 'for the initial inflation standard deviation.  Inflation mean and sd values now'
write(*,*) 'are formatted in the files exactly like restart values, e.g. arranged by variable'
write(*,*) 'type like T, U, V, etc.'
write(*,*) ''
write(*,*) 'Here is an example using version 4.4.2 or later of the NCO tools:'
write(*,*) '  ncap2 -s "T=1.0;U=1.0;V=1.0" wrfinput_d01 prior_inf.nc'
write(*,*) '  ncap2 -s "T=0.6;U=0.6;V=0.6" wrfinput_d01 prior_sd.nc'
write(*,*) ''

end program fill_inflation_restart

! <next few lines under version control, do not edit>
! $URL$ 
! $Id$ 
! $Revision$ 
! $Date$ 
