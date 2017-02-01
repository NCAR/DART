! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$
 
MODULE module_couple_uv

  use         types_mod, only : r8

CONTAINS

!------------------------------------------------------------------------

  SUBROUTINE couple_uvw ( u, v, w, mu, mub,  msfu, msfv, msfm, &
       ids, ide, jds, jde, kds, kde )

    IMPLICIT NONE

!--Input data.

    INTEGER, INTENT(IN) :: ids, ide, jds, jde, kds, kde

    REAL(r8), DIMENSION(ids:ide+1, jds:jde  , kds:kde),   INTENT(INOUT) :: u
    REAL(r8), DIMENSION(ids:ide  , jds:jde+1, kds:kde),   INTENT(INOUT) :: v
    REAL(r8), DIMENSION(ids:ide  , jds:jde,   kds:kde+1), INTENT(INOUT) :: w

    REAL(r8), DIMENSION(ids:ide+1, jds:jde  ),          INTENT(IN) :: msfu
    REAL(r8), DIMENSION(ids:ide  , jds:jde+1),          INTENT(IN) :: msfv
    REAL(r8), DIMENSION(ids:ide  , jds:jde  ),          INTENT(IN) :: msfm

    REAL(r8), DIMENSION(ids:ide  , jds:jde  ),          INTENT(IN) :: mu, mub

    REAL(r8), ALLOCATABLE, DIMENSION(:, :) :: muu, muv, muw

    allocate(muu(ids:ide+1, jds:jde  ))
    allocate(muv(ids:ide  , jds:jde+1))
    allocate(muw(ids:ide  , jds:jde  ))

!--couple variables u, v

    CALL calc_mu_uvw ( mu, mub, muu, muv, muw, &
         ids, ide, jds, jde )

    CALL couple ( muu, u, msfu, &
         ids, ide+1, jds, jde, kds, kde )

    CALL couple ( muv, v, msfv, &
         ids, ide, jds, jde+1, kds, kde )

    CALL couple ( muw, w, msfm, &
         ids, ide, jds, jde, kds, kde+1 )

    deallocate(muu)
    deallocate(muv)
    deallocate(muw)

  END SUBROUTINE couple_uvw

!-------------------------------------------------------------------------------

  SUBROUTINE calc_mu_uvw ( mu, mub, muu, muv, muw, &
       ids, ide, jds, jde )

    IMPLICIT NONE

!--Input data

    INTEGER, INTENT(IN) :: ids, ide, jds, jde

    REAL(r8), DIMENSION(ids:ide, jds:jde),     INTENT(IN)  :: mu,  mub

    REAL(r8), DIMENSION(ids:ide+1, jds:jde  ), INTENT(OUT) :: muu
    REAL(r8), DIMENSION(ids:ide  , jds:jde+1), INTENT(OUT) :: muv
    REAL(r8), DIMENSION(ids:ide  , jds:jde  ), INTENT(OUT) :: muw

    REAL(r8), DIMENSION(ids-1:ide+1, jds-1:jde+1) :: mut

    INTEGER :: i, j

    DO j=jds,jde
       DO i=ids,ide
          mut(i,j) = mu(i,j)+mub(i,j)
       ENDDO

       mut(ids-1,j) = mut(ids,j)
       mut(ide+1,j) = mut(ide,j)
    ENDDO

    DO i=ids-1,ide+1
       mut(i,jds-1)=mut(i,jds)
       mut(i,jde+1)=mut(i,jde)
    ENDDO

    DO j=jds,jde
       DO i=ids,ide+1
          muu(i,j) = 0.5_r8*(mut(i,j)+mut(i-1,j))
       ENDDO
    ENDDO

    DO j=jds,jde+1
       DO i=ids,ide
          muv(i,j) = 0.5_r8*(mut(i,j)+mut(i,j-1))
       ENDDO
    ENDDO

    DO j=jds,jde
       DO i=ids,ide
          muw(i,j) = mut(i,j)
       ENDDO
    ENDDO

  END SUBROUTINE calc_mu_uvw

!-------------------------------------------------------------------------------

  SUBROUTINE couple ( mut, field, msf, &
       ids, ide, jds, jde, kds, kde )

    IMPLICIT NONE

!--Input data

    INTEGER, INTENT(IN) :: ids, ide, jds, jde, kds, kde

    REAL(r8), DIMENSION(ids:ide, jds:jde),          INTENT(IN   ) :: mut, msf
  
    REAL(r8), DIMENSION(ids:ide, jds:jde, kds:kde), INTENT(INOUT) :: field
  
!--Local data
  
    INTEGER :: i, j, k
  
    DO j=jds,jde
       DO k=kds,kde
          DO i=ids,ide
             field(i,j,k)=field(i,j,k)*mut(i,j)/msf(i,j)
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE couple

END MODULE module_couple_uv

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
