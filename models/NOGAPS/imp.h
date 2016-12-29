!.............................START PROLOGUE............................
!
! This code may (or may not) be part of the NOGAPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$
!
! SCCS IDENTIFICATION:   $orgHeadURL: https://svn2.assembla.com/svn/ngdart/imp.h $
! SCCS IDENTIFICATION:   @(#)$orgId: imp.h 28 2010-03-03 20:09:43Z thoar $
!                        15:46:32 /h/cm/library/nognav/src/inc/fcst/imp.h_v
!
! CONFIGURATION IDENTIFICATION:  Version 4.0
!
! MODULE NAME:            imp.h
! DESCRIPTION:
!             This include file parameterizes real variables, e.g. three = 3.0
!
! CONTRACT NUMBER
!    AND TITLE:           None
! REFERENCES:             None
! CLASSIFICATION:         Unclassified
! RESTRICTIONS:           None
! COMPUTER/OPERATING
!    SYSTEM DEPENDENCIES: None
!
! LIBRARIES OF RESIDENCE: None
!
! USAGE:
! include 'imp.h'
!
! PARAMETERS:             None
!
! COMMON BLOCKS:          None
!
! FILES:                  None
! DATA BASES:             None
! NON-FILE INPUT/OUTPUT:  None
! ERROR CONDITIONS:       None
! ADDITIONAL COMMENTS:    None
!
!...................... MAINTENANCE SECTION ..............................
!
! MODULES CALLED:         None
!
! LOCAL VARIABLES AND
!            STRUCTURES:
! Name			Type	Description
! -------------------	------	-----------
! four			Real    constant, 4.0
! haf			Real    constant, 0.5
! one			Real    constant, 1.0
! three			Real    constant, 3.0
! two			Real    constant, 2.0
! zero			Real    constant, 0.0
!
! METHOD:
! Local variables:
! END OF DEFINITIONS
!
! INCLUDE FILES:          None
!
! COMPILER DEPENDENCIES:  FORTRAN 77
! COMPILE OPTIONS:        None
! MAKEFILE:               makefile
!
!  RECORD OF CHANGES:
! <<CHANGE NOTICE>> version 2.0 (1 Nov 1995) -- Rosmond, T.
!   FORTRAN 77 version
! <<CHANGE NOTICE>> version 2.1 (1 Dec 1995) -- Caudle, K.
!   Standardization modifications.
!
!........................ END PROLOGUE ...................................
 
! Local variables:
      real four
      real haf
      real one
      real three
      real two
      real zero
 
! END OF DEFINITIONS
 
!
      parameter(zero=0.0)
      parameter(one=1.0)
      parameter(two=2.0)
      parameter(three=3.0)
      parameter(four=4.0)
      parameter(haf=0.5)
!
