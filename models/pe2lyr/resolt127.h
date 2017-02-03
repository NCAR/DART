! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

! physical constants for earth

 real, parameter :: cp = 1004
 real, parameter :: g = 9.80
 real, parameter :: r = 287. 
 real, parameter :: p0 = 1.e5          
 real, parameter :: rot = 7.292e-5     
 real, parameter :: rsphere = 6.3712e6
 
!==> ztop is the pseudo-height at upper boundary (free surface).
!    tdrag and tdiab are timescales for lower layer drag and
!    thermal relaxation.
!    delth is static stability (difference in theta between two
!    levels).
!    ndiss is order of hyper-diffusion, acts to damp largest
!    total wavenumber on a timescale of efold seconds.
 
 integer, parameter :: ndiss = 8
 real, parameter :: efold = 6*3600.
 real, parameter :: tdrag = 4.0*86400.
 real, parameter :: tdiab = 20.0*86400.
 real, parameter :: theta1 = 280.
 real, parameter :: delth = 30.
 real, parameter :: ztop = 1.5e4
 real, parameter :: delta_pi = 0.9
 real, parameter :: hmax = 2000.

! time step info
! istpday=72 works for t31

 integer, parameter :: istpday = 72*4
 real, parameter :: dt = 86400./istpday

! resolution info

 integer, parameter :: ntrunc = 127
 integer, parameter :: nlons = 96*4
 integer, parameter :: nlats = nlons/2
 integer, parameter :: nlevs = 2
 integer, parameter :: model_size = 2*3*nlats*nlons

! <next few lines under version control, do not edit>
! $URL$
! $Revision$
! $Date$
