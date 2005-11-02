MODULE module_initialize

  USE module_model_constants
  USE module_ideal
  USE module_interpolations
  
CONTAINS
  
  SUBROUTINE initgen(bamex,&
       &u_init,v_init,t_init,q_init,&
       &u_g_init,v_g_init,&
       &glw_init,gsw_init,&
       &p_init,th2_init,q2_init,u10_init,v10_init,&
       &tsk_init,qsfc_init,z_gen,nz_gen,z_g,ntimes,times,&
       &nsplinetimes,splinetimes,&
       &pblh_ref,Stepbl,Lowlyr,Ht,Pblh,&
       &Th_phy,T_phy,Moist,U_phy,V_phy,Tke_myj,&
       &U_g_gen,V_g_gen,Glw_gen,Gsw_gen,P_gen,P8w_gen,&
       &Z,Z8w,Dz8w,&
       &n_moist,P_QV,ims,ime,jms,jme,kms,kme)

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: bamex
    INTEGER, INTENT(IN) :: nz_gen,ntimes,nsplinetimes
    INTEGER, INTENT(IN) :: n_moist,P_QV,ims,ime,jms,jme,kms,kme
    INTEGER, INTENT(out) :: stepbl

    REAL, DIMENSION(1:nz_gen), INTENT(IN) :: &
         &u_init,v_init,t_init,q_init
    REAL, DIMENSION(1:nz_gen,ntimes), INTENT(IN) :: &
         &u_g_init,v_g_init
    REAL, DIMENSION(1:nz_gen,ntimes), INTENT(IN) :: &
         z_gen,p_init
    REAL, DIMENSION(ntimes) :: glw_init,gsw_init
    REAL, DIMENSION(ntimes), INTENT(IN) :: times
    REAL, DIMENSION(nsplinetimes), INTENT(IN) :: splinetimes    
    REAL, DIMENSION(kms:kme,nsplinetimes), INTENT(OUT) :: &
         &U_g_gen,V_g_gen,P_gen,P8w_gen
    REAL, DIMENSION(nsplinetimes), INTENT(OUT) :: &
         &Glw_gen,Gsw_gen

    REAL, INTENT(IN) :: tsk_init,qsfc_init,&
         &th2_init,q2_init,u10_init,v10_init

    REAL, INTENT(IN) :: pblh_ref,z_g

    INTEGER, DIMENSION(IMS:IME,JMS:JME), INTENT(OUT) :: lowlyr

    REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) ::   &
         &ht,pblh

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(OUT) ::&
         &t_phy,th_phy,u_phy,v_phy,tke_myj,&
         &z,dz8w
    
    REAL, DIMENSION(ims:ime, kms:kme+1, jms:jme) :: z8w


    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, n_moist), INTENT(OUT) :: &
         &moist

    INTEGER :: nz 
    REAL :: p2,psfc,t2,dtdz
    REAL, DIMENSION(nz_gen+2) :: z_var, var,atmp,btmp,ctmp
    REAL, DIMENSION(ntimes) :: aatmp,bbtmp,cctmp
    REAL, DIMENSION(kms:kme,ntimes) :: ptmp,p8wtmp,utmp,vtmp
    INTEGER :: i,j,k,kk

! do temporal spline on gsw, glw

    CALL spline(ntimes,times,glw_init,aatmp,bbtmp,cctmp)
    DO i=1,nsplinetimes
       glw_gen(i)=seval(ntimes,splinetimes(i),times,glw_init,&
            &aatmp,bbtmp,cctmp)
    ENDDO

    CALL spline(ntimes,times,gsw_init,aatmp,bbtmp,cctmp)
    DO i=1,nsplinetimes
       gsw_gen(i)=MAX(0.,seval(ntimes,splinetimes(i),times,gsw_init,&
            &aatmp,bbtmp,cctmp))
    ENDDO

    nz = kme

    OPEN(55,file='grid_wrf1d.ascii') 

    z8w(:,1,:)=0.
    DO k=1,nz
       READ(55,*)i,z8w(:,k+1,:)
    ENDDO

    CLOSE(55)

    DO k=1,nz
       z(:,k,:)=.5*(z8w(:,k,:)+z8w(:,k+1,:))
    ENDDO


! estimate p2 and psfc using splines
! and do spline interpolation on p - don't bother with the 
! hydrostatic equation 


    DO i=1,ntimes
       DO k=1,nz_gen
          z_var(k)=z_gen(k,i)
          var(k)=p_init(k,i)
       ENDDO

       CALL spline(nz_gen,z_var,var,atmp,btmp,ctmp)
       
       DO k=1,nz
          ptmp(k,i)=&
!            &linear(z_var,var,nz_gen,z(i,k,j))
               &seval(nz_gen,z(1,k,1),z_var,var,atmp,btmp,ctmp)
          p8wtmp(k,i)=&
!            &linear(z_var,var,nz_gen,z(i,k,j))
               &seval(nz_gen,z8w(1,k,1),z_var,var,atmp,btmp,ctmp)
          
          IF (i==1) THEN
             p2=seval(nz_gen,2.,z_var,var,atmp,btmp,ctmp)
             psfc=seval(nz_gen,0.,z_var,var,atmp,btmp,ctmp)
          ENDIF
       ENDDO


    ENDDO

    DO k=1,nz
!       CALL spline(ntimes,times,p_init(k,:),aatmp,bbtmp,cctmp)
       DO i=1,nsplinetimes
          p_gen(k,i)=&
!               &seval(ntimes,splinetimes(i),times,ptmp(k,:),&
!               &aatmp,bbtmp,cctmp)
               &linear(times,ptmp(k,:),ntimes,splinetimes(i))
          p8w_gen(k,i)=&
!               &seval(ntimes,splinetimes(i),times,p8wtmp(k,:),&
!               &aatmp,bbtmp,cctmp)
               &linear(times,p8wtmp(k,:),ntimes,splinetimes(i))
       ENDDO
    ENDDO

! do linear interpolation t

    IF (bamex) THEN
       t2=th2_init
    ELSE
       t2=th2_init*(p2/p1000mb)**rcp
    ENDIF

    z_var(1)=0.
    var(1)=tsk_init
    z_var(2)=2.
    var(2)=t2

    DO k=1,nz_gen
       z_var(k+2)=z_gen(k,1)
       var(k+2)=t_init(k)
    ENDDO

!    CALL spline(nz_gen+2,z_var,var,atmp,btmp,ctmp)

    DO i=ims,ime
       DO j=jms,jme
          DO k=kms,kme
             t_phy(i,k,j)=&
                  &linear(z_var,var,nz_gen+2,z(i,k,j))
!                  &seval(nz_gen+2,z(i,k,j),z_var,var,atmp,btmp,ctmp)
             
             th_phy(i,k,j)=t_phy(i,k,j)*(p1000mb/p_gen(k,1))**rcp
          ENDDO
       ENDDO
    ENDDO

!    DO k=1,nz_gen+2
!       WRITE(150,*)var(k),z_var(k)
!    ENDDO
!    DO k=1,nz
!       WRITE(160,*)t(:,k,:),z(:,k,:)
!    ENDDO
!    RETURN
   

! do linear interpolation q

    IF (bamex) THEN

       z_var(1)=2.
       var(1)=q2_init
       
       DO k=1,nz_gen
          z_var(k+1)=z_gen(k,1)
          var(k+1)=q_init(k)
       ENDDO

    ELSE

       z_var(1)=0.
       var(1)=qsfc_init
       z_var(2)=2.
       var(2)=q2_init
       
       DO k=1,nz_gen
          z_var(k+2)=z_gen(k,1)
          var(k+2)=q_init(k)
       ENDDO
       
    ENDIF

!    CALL spline(nz_gen+2,z_var,var,atmp,btmp,ctmp)
       
    moist=0.
    DO i=ims,ime
       DO j=jms,jme
          DO k=kms,kme
             moist(i,k,j,P_QV)=&
                  &linear(z_var,var,nz_gen+2,z(i,k,j))
!                  &seval(nz_gen+2,z(i,k,j),z_var,var,atmp,btmp,ctmp)
          ENDDO
       ENDDO
    ENDDO

! do linear interpolation u

    z_var(1)=0.
    var(1)=0.
    z_var(2)=10.
    var(2)=u10_init

    DO k=1,nz_gen
       z_var(k+2)=z_gen(k,1)
       var(k+2)=u_init(k)
    ENDDO


!    CALL spline(nz_gen+2,z_var,var,atmp,btmp,ctmp)

    DO i=ims,ime
       DO j=jms,jme
          DO k=kms,kme
             u_phy(i,k,j)=&
                  &linear(z_var,var,nz_gen+2,z(i,k,j))
!                  &seval(nz_gen+2,z(i,k,j),z_var,var,atmp,btmp,ctmp)
          ENDDO
       ENDDO
    ENDDO

! do linear interpolation v

    z_var(1)=0.
    var(1)=0.
    z_var(2)=10.
    var(2)=v10_init

    DO k=1,nz_gen
       z_var(k+2)=z_gen(k,1)
       var(k+2)=v_init(k)
    ENDDO

!    CALL spline(nz_gen+2,z_var,var,atmp,btmp,ctmp)

    DO i=ims,ime
       DO j=jms,jme
          DO k=kms,kme
             v_phy(i,k,j)=&
                  &linear(z_var,var,nz_gen+2,z(i,k,j))
!                  &seval(nz_gen+2,z(i,k,j),z_var,var,atmp,btmp,ctmp)
          ENDDO
       ENDDO
    ENDDO

!    DO k=1,nz_gen+2
!       WRITE(150,*)var(k),z_var(k)
!    ENDDO
!    DO k=1,nz
!       WRITE(160,*)v_phy(:,k,:),z(:,k,:)
!    ENDDO
!    RETURN

    kk=1
    DO WHILE(z(1,kk,1) < z_g)
       kk=kk+1
    ENDDO

!u_g 

    DO i=1,ntimes
       DO k=1,nz_gen
          z_var(k)=z_gen(k,i)
          var(k)=u_g_init(k,i)
       ENDDO

       DO k=kk,nz
          utmp(k,i)=&
               &linear(z_var,var,nz_gen,z(1,k,1))
       ENDDO

       utmp(1:kk-1,i)=utmp(kk,i)

    ENDDO

    DO k=1,nz
       DO i=1,nsplinetimes
          u_g_gen(k,i)=&
               &linear(times,utmp(k,:),ntimes,splinetimes(i))
       ENDDO
    ENDDO

! v_g

    DO i=1,ntimes
       DO k=1,nz_gen
          z_var(k)=z_gen(k,i)
          var(k)=v_g_init(k,i)
       ENDDO

       DO k=kk,nz
          vtmp(k,i)=&
               &linear(z_var,var,nz_gen,z(1,k,1))
       ENDDO

       vtmp(1:kk-1,i)=vtmp(kk,i)

    ENDDO

    DO k=1,nz
       DO i=1,nsplinetimes
          v_g_gen(k,i)=&
               &linear(times,vtmp(k,:),ntimes,splinetimes(i))
       ENDDO
    ENDDO

    stepbl=1
    ht=0.
    pblh=pblh_ref
    lowlyr=1
    tke_myj=1.e-3

    DO k=kms,kme-1
       DO i=ims,ime
          DO j=jms,jme
             DZ8W(i,k,j)=z8w(i,k+1,j)-z8w(i,k,j)
             IF (k==kme-1) THEN
                DZ8W(i,k+1,j)=2.*DZ8W(i,k,j)-DZ8W(i,k-1,j)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    
  END SUBROUTINE initgen


  SUBROUTINE initgrid(U,U_g,V_g,V,T,TH,Exn,Q,P,P8w,Rho,Z,Z8w,Zo,&
       &tsk)
    
    USE module_namelist
    
    IMPLICIT NONE
    
    REAL, DIMENSION(1:nz) :: u,u_g,v,v_g,t,th,p,exn,q,rho,z

    REAL, DIMENSION(nz+1):: p8w,z8w
    REAL :: tsk

!    INTEGER :: nz

    REAL :: zo,dtdz
    
    INTEGER :: i,k

    REAL :: temp

    OPEN(55,file='grid_wrf1d.ascii') 

    z8w(1)=0.
    DO k=1,nz
       READ(55,*)i,z8w(k+1)
    ENDDO
    
    CLOSE(55)

    DO k=1,nz
       z(k)=.5*(z8w(k)+z8w(k+1))
    ENDDO

    dtdz=dtdz_ref
    
    p8w(1)=ps_ref
    DO k=1,nz
       t(k)=tsk-dtdz*z(k)
       p(k)=ps_ref*(t(k)/ts_ref)**(g/(dtdz*r_d))
       p8w(k+1)=ps_ref*&
            &((t(k)-dtdz*(z8w(k+1)-z(k)))/ts_ref)**(g/(dtdz*r_d))
       exn(k)=(p(k)/ps_ref)**rcp
       q(k)=qsrat_ref*(1.-z(k)/z8w(nz+1))*qsat(p(k),t(k))
       u_g(k)=u_g_ref
       v_g(k)=v_g_ref
       u(k)=u_g_ref
       v(k)=v_g_ref
       th(k)=t(k)/exn(k)
       rho(k)=p(k)/(r_d*t(k))
    ENDDO

    zo=zo_ref
    
  END SUBROUTINE initgrid

  SUBROUTINE initvar(z_init,z8w_init,u_init,v_init,t_init,th_init,&
       &exn_init,q_init,p_init,p8w_init,rho_init,tsoil,qsoil,&
       &zo_init,itimestep,stepbl,&
       &lowlyr,ht,znt,pblh,tsk,qsfc,mavail,&
       &p_phy,p8w,th_phy,t_phy,moist,u_phy,v_phy,pi_phy,rho,tke_myj,&
       &z,z8w,dz8w,num_soil_layers,&
       &ims,ime,jms,jme,kms,kme)

    USE module_namelist

    IMPLICIT NONE

    INTEGER :: ims,ime,jms,jme,kms,kme,num_soil_layers

    INTEGER :: itimestep,stepbl

    INTEGER,DIMENSION(IMS:IME,JMS:JME) :: LOWLYR

    REAL :: tsoil,thsoil,qsoil,zo_init

    REAL, DIMENSION(kms:kme) :: z_init,z8w_init,u_init,v_init,t_init,&
         &th_init,exn_init,q_init,p_init,p8w_init,rho_init

    REAL, DIMENSION( ims:ime, jms:jme ) ::   &
         &QSFC,TSK,xland,ht,znt,pblh,mavail

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme) ::P_PHY,p8w,PI_PHY,&
         &T_PHY,TH_PHY,U_PHY,V_PHY,RHO,tke_myj,&
         &Z,z8w,DZ8W

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, n_moist) :: moist

    INTEGER :: i,j,k

    itimestep=1
    stepbl=1

    ht=0.
    znt=zo_init
    qsfc=qsoil
    tsk=tsoil
    pblh=pblh_ref
    mavail=mavail_ref
    lowlyr=1
    moist=0.

    DO k=kms,kme
       DO i=ims,ime
          DO j=jms,jme
             P_PHY(i,k,j)=p_init(k)
             p8w(i,k,j)=p8w_init(k)
             PI_PHY(i,k,j)=exn_init(k)
             T_PHY(i,k,j)=t_init(k)
             TH_PHY(i,k,j)=th_init(k)
             moist(i,k,j,P_QV)=q_init(k)
             U_PHY(i,k,j)=u_init(k)
             V_PHY(i,k,j)=v_init(k)
             RHO(i,k,j)=rho_init(k)
             tke_myj(i,k,j)=1.e-2
             z(i,k,j)=z_init(k)
             z8w(i,k,j)=z8w_init(k)
          ENDDO
       ENDDO
    ENDDO
 
    DO k=kms,kme-1
       DO i=ims,ime
          DO j=jms,jme
             DZ8W(i,k,j)=z8w(i,k+1,j)-z8w(i,k,j)
             IF (k==kme-1) THEN
                DZ8W(i,k+1,j)=2.*DZ8W(i,k,j)-DZ8W(i,k-1,j)
             ENDIF
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE initvar

END MODULE module_initialize


