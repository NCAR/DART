! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program Precision_Check

      implicit none

      integer, parameter :: i4 = SELECTED_INT_KIND(8)
      integer, parameter :: i8 = SELECTED_INT_KIND(13)
      integer, parameter :: r4 = SELECTED_REAL_KIND(6,30)
      integer, parameter :: c4 = SELECTED_REAL_KIND(6,30)
      integer, parameter :: r8 = SELECTED_REAL_KIND(12)
      integer, parameter :: c8 = SELECTED_REAL_KIND(12)

      integer     :: inative
      integer(i4) :: int4
      integer(i8) :: int8

      real        :: rnative
      real(r4)    :: realr4
      real(r8)    :: realr8

      double precision :: dp ! deliberate use of 'double precision'

      write(*,*)
      write(*,*)'This explores the use of the intrinisc &
                &SELECTED_[REAL,INT]_KIND() functions'
      write(*,*)'and the interplay with the compiler options. &
                &You are encouraged to use the'
      write(*,*)'"autopromotion" flags on your compiler and compare &
                &the results.'
      write(*,*)

      write(*,'(''----------------------------------------------'')')
      write(*,*)'"integer"'
      write(*,*)'DIGITS    = ',   digits(inative)
      write(*,*)'HUGE      = ',     huge(inative)
      write(*,*)'KIND      = ',     kind(inative)

      write(*,'(''----------------------------------------------'')')
      write(*,*)'"integer(i4)" i4 = SELECTED_INT_KIND(8)'
      write(*,*)'DIGITS    = ',   digits(int4)
      write(*,*)'HUGE      = ',     huge(int4)
      write(*,*)'KIND      = ',     kind(int4)

      write(*,'(''----------------------------------------------'')')
      write(*,*)'"integer(i8)" i8 = SELECTED_INT_KIND(13)'
      write(*,*)'DIGITS    = ',   digits(int8)
      write(*,*)'HUGE      = ',     huge(int8)
      write(*,*)'KIND      = ',     kind(int8)

      write(*,'(''----------------------------------------------'')')
      write(*,*)'"real"'
      write(*,*)'DIGITS    = ',   digits(rnative)
      write(*,*)'EPSILON   = ',  epsilon(rnative)
      write(*,*)'HUGE      = ',     huge(rnative)
      write(*,*)'KIND      = ',     kind(rnative)
      write(*,*)'PRECISION = ',precision(rnative)

      write(*,'(''----------------------------------------------'')')
      write(*,*)'"real(r4)" r4 = SELECTED_REAL_KIND(6,30)'
      write(*,*)'DIGITS    = ',   digits(realr4)
      write(*,*)'EPSILON   = ',  epsilon(realr4)
      write(*,*)'HUGE      = ',     huge(realr4)
      write(*,*)'KIND      = ',     kind(realr4)
      write(*,*)'PRECISION = ',precision(realr4)

      write(*,'(''----------------------------------------------'')')
      write(*,*)'"real(r8)" r8 = SELECTED_REAL_KIND(13)'
      write(*,*)'DIGITS    = ',   digits(realr8)
      write(*,*)'EPSILON   = ',  epsilon(realr8)
      write(*,*)'HUGE      = ',     huge(realr8)
      write(*,*)'KIND      = ',     kind(realr8)
      write(*,*)'PRECISION = ',precision(realr8)

      write(*,'(''----------------------------------------------'')')
      write(*,*)'"double precision"'
      write(*,*)'DIGITS    = ',   digits(dp)
      write(*,*)'EPSILON   = ',  epsilon(dp)
      write(*,*)'HUGE      = ',     huge(dp)
      write(*,*)'KIND      = ',     kind(dp)
      write(*,*)'PRECISION = ',precision(dp)

      end program Precision_Check

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
