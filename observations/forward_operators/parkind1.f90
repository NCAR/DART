! Description:
!> @file
!!   Defines the KINDs used by RTTOV
!
!> @brief
!!   Defines the KINDs used by RTTOV
!!
!! @details
!!   The user-level interface only requires jprb, jpim and jplm.
!!   This module can be edited, for example, to run RTTOV with
!!   single precision, but such changes are not recommended.
!
MODULE PARKIND1
!
!     *** Define usual kinds for strong typing ***
!     Option of Kinds modified for NEC interface Roger Saunders
!
IMPLICIT NONE
SAVE
!
!     Integer Kinds
!     -------------
!
INTEGER, PARAMETER :: JPIT = SELECTED_INT_KIND(2)       !< Tiny integer (only used internally by RTTOV)
INTEGER, PARAMETER :: JPIS = SELECTED_INT_KIND(4)       !< Small integer (only used internally by RTTOV)
INTEGER, PARAMETER :: JPIB = SELECTED_INT_KIND(12)      !< Big integer (only used internally by RTTOV)
INTEGER, PARAMETER :: JPIM = SELECTED_INT_KIND(9)       !< Standard integer type
!INTEGER, PARAMETER :: JPIM = SELECTED_INT_KIND(15)     !Required for METO IBM OPS code (not standalone)

! Special integer type to be used for sensitive adress calculations
! should be *8 for a machine with 8byte adressing for optimum performance
!ifdef ADDRESS64
INTEGER, PARAMETER :: JPIA = JPIB                       !< Not used by RTTOV
!#else
!INTEGER, PARAMETER :: JPIA = JPIM
!#endif
!
!     Real Kinds
!     ----------
!
INTEGER, PARAMETER :: JPRT = SELECTED_REAL_KIND(2,1)    !< Tiny real (not used by RTTOV)
INTEGER, PARAMETER :: JPRS = SELECTED_REAL_KIND(4,2)    !< Small real (not used by RTTOV)
INTEGER, PARAMETER :: JPRM = SELECTED_REAL_KIND(6,37)   !< Medium real (only used internally by RTTOV)
INTEGER, PARAMETER :: JPRB = SELECTED_REAL_KIND(13,300) !< Standard real type
!
!     Logical Kinds
!     -------------
INTEGER, PARAMETER :: JPLM = KIND(.TRUE.)               !< Standard logical type
!
END MODULE PARKIND1
