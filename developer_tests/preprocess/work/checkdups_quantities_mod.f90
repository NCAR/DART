! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! ! in this section, define the quantity of interest.  it must
! ! start with QTY_xxx and be less than 32 characters total.
! ! if there are units, min/max values, pdf, etc you can add one or
! ! more name=value pairs.  at this point there *cannot* be spaces
! ! around the equals sign.  the value can have spaces if it is
! ! surrounded by double quotes.  e.g. desc="assumes dry air density"
! ! comment lines (like these) can be added if they start with an additional !
 
! BEGIN DART PREPROCESS QUANTITY DEFINITIONS
!
!  ! QTY_STATE_VARIABLE is predefined in the preprocess program.
!  ! please do not change this first entry.
!
!  QTY_STATE_VARIABLE               desc="basic item in a state"
!  QTY_1D_INTEGRAL                  desc="compute value with an integral"
!  QTY_STATE_VAR_POWER     pdf=log  desc="raising a state value to a power"
!  QTY_LARGE_SCALE_STATE            desc="state varies with large time/space scale"
!  QTY_SMALL_SCALE_STATE            desc="state varies with small time/space scale"
! 
! ! tests for duplicate difference detection:
!    QTY_STATE_VARIABLE
!    QTY_STATE_VARIABLE               desc="bob"
!    QTY_STATE_VARIABLE               start="bob"
!    QTY_STATE_VARIABLE   start=bob   desc="basic item in a state"
!    QTY_STATE_VARIABLE               desc="basic item in a state"  end=bob
!    QTY_STATE_VARIABLE   start=bob   desc="basic item in a state"  end=bob
!    QTY_STATE_VARIABLE               desc="basic item in a state"
!
! END DART PREPROCESS QUANTITY DEFINITIONS


