MODULE module_initialize

  USE module_model_constants
  USE module_ideal
  USE module_interpolations
  USE module_namelist, only : force_uvg, force_flux, t_advection, n_moist, &
                              P_QV, init_f_type, qv_advection
  
CONTAINS
  
  SUBROUTINE initf(u_init,v_init, &
       t_init,q_init,&
       u_g_init,v_g_init, &
       t_init_upstream_x, t_init_upstream_y, &
       qv_init_upstream_x, qv_init_upstream_y, &
       tau_init_u, tau_init_v, &
       glw_init,gsw_init,&
       precip_init,&
       ts_init, &
       qvs_init, &
       ustar_init, &
       hflux_init, &
       qvflux_init,&
       p_init, &
       th2_init,q2_init, &
       u10_init,v10_init,&
       tsk_init,qsfc_init, &
       z_f,nz_f,z_g,ntimes,ntimes_flux,ntimes_smos,ntimes_sfc,&
       times,times_flux,&
       times_smos,&
       times_sfc,&
       nsplinetimes,splinetimes,&
       nsplinetimes_advection,splinetimes_advection,&
       nsplinetimes_flux,splinetimes_flux,&
       nsplinetimes_smos,splinetimes_smos,&
       nsplinetimes_sfc,splinetimes_sfc,&
       pblh_ref,stepbl,lowlyr,ht,pblh,&
       th_phy,t_phy,moist,u_phy,v_phy,p_phy, p8w, tke_myj,&
       u_g_f,v_g_f,glw_f,gsw_f,precip_f,p_f,p8w_f,&
       th_upstream_x, th_upstream_y, &
       qv_upstream_x, qv_upstream_y, &
       tau_u, tau_v, &
       ts_f,qvs_f,ustar_f,hflux_f,qvflux_f,&
       z,z8w,dz8w,&
       ims,ime,jms,jme,kms,kme)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nz_f,ntimes,ntimes_flux, ntimes_smos,&
         &ntimes_sfc,&
         nsplinetimes,nsplinetimes_flux, nsplinetimes_smos,&
         &nsplinetimes_sfc, nsplinetimes_advection


    INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
    INTEGER, INTENT(out) :: stepbl

    REAL, DIMENSION(1:nz_f), INTENT(IN) :: u_init,v_init,t_init,q_init

    REAL, DIMENSION(1:nz_f,ntimes), INTENT(IN) :: u_g_init,v_g_init, &
         z_f,p_init, t_init_upstream_x, t_init_upstream_y, &
         qv_init_upstream_x, qv_init_upstream_y, &
         tau_init_u, tau_init_v

    REAL, DIMENSION(ntimes_flux) :: glw_init,gsw_init
    REAL, DIMENSION(ntimes_smos) :: precip_init 
    REAL, DIMENSION(ntimes_sfc) :: ts_init,qvs_init,&
         ustar_init,hflux_init,qvflux_init

    REAL, DIMENSION(ntimes), INTENT(IN) :: times
    REAL, DIMENSION(ntimes_flux), INTENT(IN) :: times_flux
    REAL, DIMENSION(ntimes_smos), INTENT(IN) :: times_smos
    REAL, DIMENSION(ntimes_sfc), INTENT(IN) :: times_sfc

    REAL, DIMENSION(nsplinetimes), INTENT(IN) :: splinetimes    
    REAL, DIMENSION(nsplinetimes_advection), INTENT(IN) :: splinetimes_advection    
    REAL, DIMENSION(nsplinetimes_flux), INTENT(IN) :: splinetimes_flux    
    REAL, DIMENSION(nsplinetimes_smos), INTENT(IN) :: splinetimes_smos    
    REAL, DIMENSION(nsplinetimes_sfc), INTENT(IN) :: splinetimes_sfc    

    REAL, DIMENSION(kms:kme,nsplinetimes), INTENT(OUT) :: &
         u_g_f,v_g_f,p_f,p8w_f, th_upstream_x, th_upstream_y, &
         qv_upstream_x, qv_upstream_y, tau_u, tau_v

    REAL, DIMENSION(nsplinetimes_flux), INTENT(OUT) :: &
         glw_f,gsw_f

    REAL, DIMENSION(nsplinetimes_smos), INTENT(OUT) :: precip_f


    REAL, DIMENSION(nsplinetimes_sfc), INTENT(OUT) :: &
         ts_f,qvs_f,ustar_f,hflux_f,qvflux_f


    REAL, INTENT(IN) :: tsk_init,qsfc_init,&
         th2_init,q2_init,u10_init,v10_init

    REAL, INTENT(IN) :: pblh_ref,z_g

    INTEGER, DIMENSION(IMS:IME,JMS:JME), INTENT(OUT) :: lowlyr

    REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) ::   &
         ht,pblh

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(OUT) ::&
         t_phy,th_phy,u_phy,v_phy,p_phy, p8w, tke_myj,&
         z,dz8w
    
    REAL, DIMENSION(ims:ime, kms:kme+1, jms:jme) :: z8w

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, n_moist), INTENT(OUT) :: &
         moist

! local
    INTEGER :: nz, near_sfc_num, skip_sfc_f, ntimes_snd, nsplinetimes_snd
    REAL :: p2,psfc,t2
    REAL, DIMENSION(nz_f+2) :: z_var, var,atmp,btmp,ctmp
    REAL, DIMENSION(ntimes_flux) :: aatmp_flux,bbtmp_flux,cctmp_flux
    REAL, DIMENSION(kms:kme,ntimes) :: ptmp,p8wtmp,utmp,vtmp,ttmp
    INTEGER :: i,j,k,kk

! only need to interpolate first vertical profiles if we have UVG
    IF ( init_f_type == "OBS" ) THEN
       ntimes_snd = 1
       nsplinetimes_snd = 1
    ELSE
       ntimes_snd = ntimes
       nsplinetimes_snd = nsplinetimes
    ENDIF

    IF ( ntimes_flux > 1 ) THEN
       CALL spline(ntimes_flux,times_flux,glw_init,aatmp_flux,bbtmp_flux,cctmp_flux)
       DO i=1,nsplinetimes_flux
          glw_f(i)=seval(ntimes_flux,splinetimes_flux(i),times_flux,glw_init,&
               &aatmp_flux,bbtmp_flux,cctmp_flux)
!print*,i,splinetimes_flux(i),glw_f(i)
       ENDDO

       CALL spline(ntimes_flux,times_flux,gsw_init,aatmp_flux,bbtmp_flux,cctmp_flux)
       DO i=1,nsplinetimes_flux
          gsw_f(i)=MAX(0.,seval(ntimes_flux,splinetimes_flux(i),times_flux,gsw_init,&
               &aatmp_flux,bbtmp_flux,cctmp_flux))
       ENDDO

    ELSE

       gsw_f(1) = gsw_init(1)
       glw_f(1) = glw_init(1)

    ENDIF

! do temporal linear on precip
    IF ( ntimes_smos > 1 ) THEN
       DO i=1,nsplinetimes_smos
          precip_f(i)=linear(times_smos,precip_init,ntimes_smos,&
               &splinetimes_smos(i))
       ENDDO

    ELSE

       precip_f(1) = precip_init(1)

    ENDIF

    IF (force_flux) THEN

       IF ( ntimes_sfc > 1 ) THEN
          DO i=1,nsplinetimes_sfc
             ts_f(i)=linear(times_sfc,ts_init,ntimes_sfc,&
                  &splinetimes_sfc(i))
             qvs_f(i)=linear(times_sfc,qvs_init,ntimes_sfc,&
                  &splinetimes_sfc(i))
             ustar_f(i)=linear(times_sfc,ustar_init,ntimes_sfc,&
                  &splinetimes_sfc(i))
             hflux_f(i)=linear(times_sfc,hflux_init,ntimes_sfc,&
                  &splinetimes_sfc(i))
             qvflux_f(i)=linear(times_sfc,qvflux_init,ntimes_sfc,&
                  &splinetimes_sfc(i))
          ENDDO
          
       ELSE
          
          ts_f(1) = ts_init(1)
          qvs_f(1) = qvs_init(1)
          ustar_f(1) = ustar_init(1)
          hflux_f(1) = hflux_init(1)
          qvflux_f(1) = qvflux_init(1)
          
       ENDIF

    ENDIF

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


    p8wtmp = -9999.0
    ptmp = -9999.0
    DO i=1,ntimes_snd
       DO k=1,nz_f
          z_var(k)=z_f(k,i)
          var(k)=p_init(k,i)
       ENDDO

       CALL spline(nz_f,z_var,var,atmp,btmp,ctmp)
       
       DO k=1,nz
          ptmp(k,i)=&
!            &linear(z_var,var,nz_f,z(i,k,j))
               &seval(nz_f,z(1,k,1),z_var,var,atmp,btmp,ctmp)
          p8wtmp(k,i)=&
!            &linear(z_var,var,nz_f,z(i,k,j))
               &seval(nz_f,z8w(1,k,1),z_var,var,atmp,btmp,ctmp)
          
          IF (i==1) THEN
             p2=seval(nz_f,2.,z_var,var,atmp,btmp,ctmp)
             psfc=seval(nz_f,0.,z_var,var,atmp,btmp,ctmp)
          ENDIF
       ENDDO

    ENDDO

    DO k=1,nz
       DO i=1,nsplinetimes_snd
          p_f(k,i)=&
               &linear(times,ptmp(k,:),ntimes,splinetimes(i))
          p8w_f(k,i)=&
               &linear(times,p8wtmp(k,:),ntimes,splinetimes(i))
       ENDDO
    ENDDO

    DO i=ims,ime
       DO j=jms,jme
          DO k=kms,kme
             p_phy(i,k,j)= p_f(k,1) 
             p8w(i,k,j)  = p8w_f(k,1) 
          ENDDO
       ENDDO
    ENDDO

! do linear interpolation t
    near_sfc_num = 0
! first fills in lowest level with tsk, then if available subs in th2 

    z_var(1:2) = z_f(1:2,1)
    IF ( tsk_init > -273 ) THEN
       near_sfc_num = near_sfc_num+1
       z_var(near_sfc_num)=0.
       var(near_sfc_num)=tsk_init
    ENDIF

    IF ( th2_init > -273 ) THEN
       IF ( z_var(near_sfc_num+1) > 2.0001 ) THEN
          near_sfc_num = near_sfc_num+1
          t2=th2_init
          var(near_sfc_num)=t2
          z_var(near_sfc_num)=2.
       ELSE
          t2=th2_init
          var(near_sfc_num+1)=t2
          z_var(near_sfc_num+1)=2.
       ENDIF
    ENDIF

    DO k=1,nz_f
       z_var(k+near_sfc_num)=z_f(k,1)
       var(k+near_sfc_num)=t_init(k)
    ENDDO

    DO i=ims,ime
       DO j=jms,jme
          DO k=kms,kme
             t_phy(i,k,j)=&
                  &linear(z_var,var,nz_f+near_sfc_num,z(i,k,j))
             th_phy(i,k,j)=t_phy(i,k,j)*(p1000mb/p_f(k,1))**rcp
          ENDDO
       ENDDO
    ENDDO

! replace missing value in th
    WHERE ( th_phy < -273 ) th_phy = -9999.0 


!    DO k=1,nz_f+2
!       WRITE(150,*)var(k),z_var(k)
!    ENDDO
!    DO k=1,nz
!       WRITE(160,*)t(:,k,:),z(:,k,:)
!    ENDDO
!    RETURN
   

! do linear interpolation q
! use qsfc if we have it

    z_var(1:2) = z_f(1:2,1)
    near_sfc_num = 0
    IF ( qsfc_init > 0.0 ) THEN
       near_sfc_num = near_sfc_num+1
       z_var(near_sfc_num)=0.
       var(near_sfc_num)=qsfc_init
    ENDIF

    IF ( q2_init > 0.0 ) THEN 
       IF ( z_var(near_sfc_num+1) > 2.0001 ) THEN
          near_sfc_num = near_sfc_num+1
          var(near_sfc_num)=q2_init
          z_var(near_sfc_num)=2.
       ELSE
          var(near_sfc_num+1)=q2_init
          z_var(near_sfc_num+1)=2.
       ENDIF
    ENDIF

    DO k=1,nz_f
       z_var(k+near_sfc_num)=z_f(k,1)
       var(k+near_sfc_num)=q_init(k)
    ENDDO

    moist=0.
    DO i=ims,ime
       DO j=jms,jme
          DO k=kms,kme
             moist(i,k,j,P_QV)=&
                  &linear(z_var,var,nz_f+near_sfc_num,z(i,k,j))
          ENDDO
       ENDDO
    ENDDO

! do linear interpolation u

    z_var(1:2) = z_f(1:2,1)
    near_sfc_num = 1
    z_var(1)=0.
    var(1)=0.

    IF ( u10_init > -500 ) THEN
       IF ( z_var(near_sfc_num+1) > 10.0001 ) THEN
          near_sfc_num = near_sfc_num+1
          var(near_sfc_num)=u10_init
          z_var(near_sfc_num)=10.0
       ELSE
          var(near_sfc_num+1)=u10_init
          z_var(near_sfc_num+1)=10.0
       ENDIF
    ENDIF
 
    skip_sfc_f = 0
    IF ( z_f(1,1) == 2.0 ) skip_sfc_f = 1
 
    DO k=1,nz_f-skip_sfc_f
       z_var(k+near_sfc_num)=z_f(k+skip_sfc_f,1)
       var(k+near_sfc_num)=u_init(k+skip_sfc_f)
    ENDDO

    DO i=ims,ime
       DO j=jms,jme
          DO k=kms,kme
             u_phy(i,k,j)=&
                  &linear(z_var,var,nz_f+near_sfc_num,z(i,k,j))
          ENDDO
       ENDDO
    ENDDO


! do linear interpolation v

    z_var(1:2) = z_f(1:2,1)
    near_sfc_num = 1
    z_var(1)=0.
    var(1)=0.

    IF ( v10_init > -500 ) THEN
       IF ( z_var(near_sfc_num+1) > 10.0001 ) THEN
          near_sfc_num = near_sfc_num+1
          var(near_sfc_num)=v10_init
          z_var(near_sfc_num)=10.0
       ELSE
          var(near_sfc_num+1)=v10_init
          z_var(near_sfc_num+1)=10.0
       ENDIF
    ENDIF

    skip_sfc_f = 0
    IF ( z_f(1,1) == 2.0 ) skip_sfc_f = 1
 
    DO k=1,nz_f-skip_sfc_f
       z_var(k+near_sfc_num)=z_f(k+skip_sfc_f,1)
       var(k+near_sfc_num)=v_init(k+skip_sfc_f)
    ENDDO

    DO i=ims,ime
       DO j=jms,jme
          DO k=kms,kme
             v_phy(i,k,j)=&
                  &linear(z_var,var,nz_f+near_sfc_num,z(i,k,j))
          ENDDO
       ENDDO
    ENDDO

!    DO k=1,nz_f+2
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

    DO i=1,ntimes_snd
       DO k=1,nz_f
          z_var(k)=z_f(k,i)
          var(k)=u_g_init(k,i)
       ENDDO

       IF ( minval(var) > -500 ) THEN

          DO k=kk,nz
             utmp(k,i)=&
                  &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO

       ELSE
          utmp(:,i) = -9999.
       ENDIF ! missing?
       utmp(1:kk-1,i)=utmp(kk,i)

    ENDDO

    DO k=1,nz
       DO i=1,nsplinetimes_snd
          u_g_f(k,i)=&
               &linear(times,utmp(k,:),ntimes,splinetimes(i))
       ENDDO
    ENDDO
! v_g

    DO i=1,ntimes_snd
       DO k=1,nz_f
          z_var(k)=z_f(k,i)
          var(k)=v_g_init(k,i)
       ENDDO

       IF ( minval(var) > -500 ) THEN

          DO k=kk,nz
             vtmp(k,i)=&
                  &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO

       ELSE
          vtmp(:,i) = -9999.
       ENDIF ! missing?

       vtmp(1:kk-1,i)=vtmp(kk,i)

    ENDDO

    DO k=1,nz
       DO i=1,nsplinetimes_snd
          v_g_f(k,i)=&
               &linear(times,vtmp(k,:),ntimes,splinetimes(i))
       ENDDO
    ENDDO

! temperature tendencies
    IF ( t_advection ) then

! t upstream x
       DO i=1,ntimes_snd
          DO k=1,nz_f
             z_var(k)=z_f(k,i)
             var(k)=t_init_upstream_x(k,i)
          ENDDO

          DO k=1,nz
             ttmp(k,i)=&
                     &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO
       ENDDO
       DO k=1,nz
          DO i=1,nsplinetimes_advection
             th_upstream_x(k,i)=&
                  &linear(times,ttmp(k,:),ntimes,splinetimes_advection(i))
          ENDDO
       ENDDO
       th_upstream_x=th_upstream_x*(p1000mb/p_f)**rcp
       WHERE(th_upstream_x<-100) th_upstream_x = -9999.0 

! t upstream y
       DO i=1,ntimes_snd
          DO k=1,nz_f
             z_var(k)=z_f(k,i)
             var(k)=t_init_upstream_y(k,i)
          ENDDO

          DO k=1,nz
             ttmp(k,i)=&
                     &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO
       ENDDO
       DO k=1,nz
          DO i=1,nsplinetimes_advection
             th_upstream_y(k,i)=&
                  &linear(times,ttmp(k,:),ntimes,splinetimes_advection(i))
          ENDDO
       ENDDO
       th_upstream_y=th_upstream_y*(p1000mb/p_f)**rcp
       WHERE(th_upstream_y<-100) th_upstream_y = -9999.0 

    ENDIF
! moisture tendencies
    IF ( qv_advection ) then

! qv upstream x
       DO i=1,ntimes_snd
          DO k=1,nz_f
             z_var(k)=z_f(k,i)
             var(k)=qv_init_upstream_x(k,i)
          ENDDO

          DO k=1,nz
             ttmp(k,i)=&
                     &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO
       ENDDO
       DO k=1,nz
          DO i=1,nsplinetimes_advection
             qv_upstream_x(k,i)=&
                  &linear(times,ttmp(k,:),ntimes,splinetimes_advection(i))
          ENDDO
       ENDDO
       WHERE(qv_upstream_x<-100) qv_upstream_x = -9999.0 

! qv upstream y
       DO i=1,ntimes_snd
          DO k=1,nz_f
             z_var(k)=z_f(k,i)
             var(k)=qv_init_upstream_y(k,i)
          ENDDO

          DO k=1,nz
             ttmp(k,i)=&
                     &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO
       ENDDO
       DO k=1,nz
          DO i=1,nsplinetimes_advection
             qv_upstream_y(k,i)=&
                  &linear(times,ttmp(k,:),ntimes,splinetimes_advection(i))
          ENDDO
       ENDDO
       WHERE(qv_upstream_y<-100) qv_upstream_y = -9999.0 

    ENDIF

    IF ( t_advection .or. qv_advection ) THEN
! tau u 
       DO i=1,ntimes_snd
          DO k=1,nz_f
             z_var(k)=z_f(k,i)
             var(k)=tau_init_u(k,i)
          ENDDO

          DO k=1,nz
             ttmp(k,i)=&
                     &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO
       ENDDO
       DO k=1,nz
          DO i=1,nsplinetimes_advection
             tau_u(k,i)=&
                  &linear(times,ttmp(k,:),ntimes,splinetimes_advection(i))
          ENDDO
       ENDDO

! tau v 
       DO i=1,ntimes_snd
          DO k=1,nz_f
             z_var(k)=z_f(k,i)
             var(k)=tau_init_v(k,i)
          ENDDO

          DO k=1,nz
             ttmp(k,i)=&
                     &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO
       ENDDO
       DO k=1,nz
          DO i=1,nsplinetimes_advection
             tau_v(k,i)=&
                  &linear(times,ttmp(k,:),ntimes,splinetimes_advection(i))
          ENDDO
       ENDDO

    ENDIF ! t_advection

 
    stepbl=1
    ht=0.
    pblh=pblh_ref
    lowlyr=1
    tke_myj=.5*epsq2

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
    
  END SUBROUTINE initf


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

    REAL :: tsoil,qsoil,zo_init

    REAL, DIMENSION(kms:kme) :: z_init,z8w_init,u_init,v_init,t_init,&
         &th_init,exn_init,q_init,p_init,p8w_init,rho_init

    REAL, DIMENSION( ims:ime, jms:jme ) ::   &
         &QSFC,TSK,ht,znt,pblh,mavail

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
             tke_myj(i,k,j)=.5*epsq2
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


