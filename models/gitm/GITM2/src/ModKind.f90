! This code may (or may not) be part of the GITM distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

!^CFG COPYRIGHT UM
!BOP
!MODULE: ModKind - define various precisions in a machine independent way

!DESCRIPTION:
! The Fortran 77 style real*4 and real*8 declarations are obsolete,
! and compilers often issue warnings. The real and double precision
! types are machine and compiler flag dependent.
! The Fortran 90 way is to define the {\bf kind} parameter.
! Typical usage:
!\begin{verbatim}
! real(Real8_) :: CpuTime  ! variable declaration
! CpuTime = 0.0_Real8_     ! 8 byte real constant
!\end{verbatim}

!INTERFACE:
module ModKind

  !PUBLIC DATA MEMBERS:
  integer, parameter :: Real4_ = selected_real_kind(6,30)
  integer, parameter :: Real8_ = selected_real_kind(12,100)
  integer, parameter :: Int8_  = selected_int_kind(16)

  ! Number of bytes in the default real number (precision)
  ! This is standard F90 initialization expression but may give warnings:
  integer, parameter :: nByteReal = 4 + (1.00000000041 - 1.0)*10000000000.0

  !EOP
end module ModKind

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
