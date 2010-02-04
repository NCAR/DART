!
MODULE module_stochastic_cloud
!
! DART $Id$
!
CONTAINS
!---------------------------------------------------------------------
! 
! By   : Robert Tardif, NCAR/RAL
! Date : June 2007
!
!---------------------------------------------------------------------
!      SUBROUTINE read_stochastic_cloud_sim()
!
!      END SUBROUTINE read_stochastic_cloud
!---------------------------------------------------------------------
      SUBROUTINE cloud_stochastic(moist,p_qc,icloud,z,PBLH,rho,t_phy, &
      & p_phy,stoch_yn,stoch_cbh,stoch_lwp,itimestep,ntime,         &
      & ims,ime,jms,jme,kms,kme,n_moist,indexrealization )
!---------------------------------------------------------------------
      IMPLICIT NONE
!---------------------------------------------------------------------
      INTEGER,      INTENT(IN   )   ::   ims,ime, jms,jme, kms,kme,  &
                                         ntime, itimestep, p_qc,     &
                                         n_moist, indexrealization

      INTEGER,      INTENT(INOUT)   ::   icloud
      
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme, n_moist ), &
            INTENT(INOUT) :: moist

      INTEGER, DIMENSION(ntime) :: stoch_yn
      REAL, DIMENSION(ntime) :: stoch_cbh, stoch_lwp

      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),          &
         INTENT(IN ) :: z, rho, t_phy, p_phy

      REAL, DIMENSION( ims:ime, jms:jme ),          &
         INTENT(IN ) :: PBLH

      INTEGER :: icbh
      INTEGER :: i, j, k

      INTEGER :: indexrun !...RT
      REAL :: randraw !...RT

      REAL :: threshold_lwp, threshold_cbh
      REAL :: bwgt, twgt

!---------------------------------------------------------------------

      threshold_lwp = 700.0
      threshold_cbh = 6.0
! code stochastic process here... (later)

! assume prior runs of stochastic model composed of 48 hours 
! w/ data every 20 secs (same as timestep of model)

      if (itimestep.eq.1) then
         write(99,*) 'in stochastic: indexrealization= ', indexrun
      endif

! did the stochastic model produce a cloud located above the model-diagnosed PBL?
         IF ((stoch_yn(itimestep).EQ.1).AND.&
         & (stoch_cbh(itimestep)*1000.0.GT.PBLH(1,1))) then
            icloud = 1

! find level corresponding to cloud base
            DO i=ims,ime
               DO j=jms,jme
                  DO k=kms,kme
                     IF ((z(i,k,j)/1000.0).GT.stoch_cbh(itimestep)) THEN
                        icbh = k-1
                        twgt = (stoch_cbh(itimestep)*1000.-z(i,icbh,j))/(z(i,k,j)-z(i,icbh,j))
                        bwgt = 1.0-twgt
                        GOTO 111
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
 111        CONTINUE

! clipping of "unrealistic" values

            IF (stoch_lwp(itimestep).gt.threshold_lwp) THEN
               stoch_lwp(itimestep) = threshold_lwp
            ENDIF

! building the cloud in the vertical

            CALL build_stochastic_cloud(moist,p_qc,icbh,stoch_lwp, &
            & itimestep, bwgt, twgt, z,rho,t_phy, p_phy, &
            & ims,ime,jms,jme,kms,kme,n_moist,ntime)

         ELSE
            icloud = 0
            DO i=ims,ime
               DO j=jms,jme
                  DO k=kms,kme
                     moist(i,k,j,p_qc) = 0.0
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         
      END SUBROUTINE cloud_stochastic

!---------------------------------------------------------------------
!---------------------------------------------------------------------

      SUBROUTINE build_stochastic_cloud(moist,p_qc,icbh,stoch_lwp,itimestep, &
      & bwgt,twgt, z,rho,t_phy,p_phy, ims,ime,jms,jme,kms,kme,n_moist,ntime)
!---------------------------------------------------------------------
      IMPLICIT NONE
!---------------------------------------------------------------------
      INTEGER,      INTENT(IN   )    ::   ims,ime, jms,jme, kms,kme, &
                                          p_qc, n_moist, itimestep, ntime

      REAL, INTENT(IN)               ::   bwgt, twgt
      
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme, n_moist ), &
            INTENT(INOUT) :: moist

      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),          &
         INTENT(IN ) :: z, rho, t_phy, p_phy

      REAL, DIMENSION(ntime) :: stoch_lwp

      INTEGER, INTENT(IN   ) :: icbh

      INTEGER :: i, j, k, icth

      REAL :: qbase, qs, PseudoLapse, tt, ews, lwp, excess, excessql, pp, &
              zz

! building the cloud in the vertical assuming adiabatic profile
! input: cloud base height and LWP

      tt =  t_phy(1,icbh,1) ! temp at cloud base
      tt = bwgt*t_phy(1,icbh,1) + twgt*t_phy(1,icbh+1,1)
      ews=satvappres(tt-273.15)*100.0
      pp = bwgt*p_phy(1,icbh,1) + twgt*p_phy(1,icbh+1,1)
      qbase = 0.622*(ews/(pp-ews))
      PseudoLapse = pseudoadiabat(tt,pp) ! pseudo-adiabat lapse rate
      lwp = 0.0
      DO i=ims,ime
         DO j=jms,jme
            DO k=icbh+1,kme                        
               zz = z(i,k-1,j)
               if ( k == icbh+1 ) zz = bwgt*z(1,icbh,1) + twgt*z(1,icbh+1,1)
               tt = tt+PseudoLapse*(z(i,k,j)-zz)
               ews=satvappres(tt-273.15)*100.0
               qs = 0.622*(ews/(p_phy(1,k,1)-ews))
               moist(i,k,j,P_QC) = (qbase - qs)
               lwp = lwp + (moist(i,k,j,P_QC)*rho(i,k,j)*(z(i,k,j)-zz))
               PseudoLapse = pseudoadiabat(tt,p_phy(i,k,j))
               if ((lwp*1000.0).gt.stoch_lwp(itimestep)) then
                  icth = k ! cloud top index
                  goto 112
               endif
            ENDDO
         ENDDO
      ENDDO
 112  continue

! determining the excess water at the upper level w.r.t. the LWP
! and taking it out at cloud top to match observed LWP
      excess = lwp - stoch_lwp(itimestep)/1000.0
      excessql = excess/(rho(1,icth,1)*(z(1,icth,1)-z(1,icth-1,1)))
      moist(1,icth,1,P_QC) = moist(1,icth,1,P_QC) - excessql

      lwp = 0.0
      DO i=ims,ime
         DO j=jms,jme
            DO k=icbh+1,icth
               zz = z(i,k-1,j)
               if ( k == icbh+1 ) zz = bwgt*z(1,icbh,1) + twgt*z(1,icbh+1,1)
               lwp = lwp + (moist(i,k,j,P_QC)*rho(i,k,j)*(z(i,k,j)-zz))
            ENDDO
         ENDDO
      ENDDO

! make sure cloud water is set to zero outside of cloud layer between [icbh, icth]
               DO i=ims,ime
                  DO j=jms,jme
                     DO k=kms,kme
                        IF ((k.lt.icbh).or.(k.gt.icth)) then
                           moist(i,k,j,P_QC) = 0.0
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO

      END SUBROUTINE build_stochastic_cloud

!---------------------------------------------------------------------
!---------------------------------------------------------------------
      FUNCTION satvappres(t)
      IMPLICIT NONE

      real :: satvappres
      real :: t

! t in deg. C and ews in hPa

      satvappres = 6.112*EXP((17.67*t)/(t+243.5)) ! from Buck (JAM, 1981)

      END FUNCTION satvappres
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      FUNCTION pseudoadiabat(t,p)
      IMPLICIT NONE

      real, INTENT(IN   ) :: t, p

      real :: pseudoadiabat

      real :: ra,rv,rva,cpv,gmd,rs,alv,alr,all,rw,a1,a2,factor,pp


! units of p from Pa to hPa

! taken from the COBEL code
      pp = p/100.0

      ra=287.0
      cpv=1.85*1004.5
      rv=461.5
      rva=rv/ra

      gmd=9.8/1004.5
      rs=satvappres(t-273.15)

      alv=1000.*(2500.-2.34*(t-273.15))
      alr=alv/ra
      all=alv/rv
      rw=0.622*rs/(pp-rs)
      a1=(1.+rw)*(1.+alr*rw/t)
      a2=cpv+all*alv*(1.+rw*rva)/(t*t)
      factor=a1/(1.+rw*a2/1004.5)

      pseudoadiabat=-gmd*(1.0 + rw)*factor


      END FUNCTION pseudoadiabat
!---------------------------------------------------------------------
!---------------------------------------------------------------------
END MODULE module_stochastic_cloud
