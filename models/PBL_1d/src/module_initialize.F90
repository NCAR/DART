MODULE module_initialize
!
! DART $Id$
!

  USE module_model_constants
  USE module_ideal
  USE module_interpolations, only : seval, spline, linear
  USE module_namelist, only : force_uvg, force_flux, t_advection, n_moist, &
                              init_f_type, qv_advection, u_advection,      &
                              P_QV,P_QC,P_QR,P_QI,P_QG, force_soil
CONTAINS
  
  SUBROUTINE initf(itime_f,u_init_f,v_init_f, &
       t_init_f,qv_init_f,qc_init_f,qr_init_f,   &
       qi_init_f,qg_init_f,                     &
       p_init_f,&
       t_init_upstream_x, t_init_upstream_y, &
       qv_init_upstream_x, qv_init_upstream_y, &
       qc_init_upstream_x, qc_init_upstream_y, &
       qr_init_upstream_x, qr_init_upstream_y, &
       qi_init_upstream_x, qi_init_upstream_y, &
       qg_init_upstream_x, qg_init_upstream_y, &
       u_init_upstream_x, u_init_upstream_y, &
       v_init_upstream_x, v_init_upstream_y, &
       t2_init_upstream_x, t2_init_upstream_y, &
       q2_init_upstream_x, q2_init_upstream_y, &
       u10_init_upstream_x, u10_init_upstream_y, &
       v10_init_upstream_x, v10_init_upstream_y, &
       tau_init_u, tau_init_v, &
       tau_init_u10, tau_init_v10, &
       glw_init,gsw_init,&
       precip_init,&
       tslb_init, smois_init, &
       ts_init, &
       qvs_init, &
       ustar_init, &
       hflux_init, &
       qvflux_init,&
       th2_init,q2_init, &
       u10_init,v10_init,&
       tsk_init,qsfc_init, &
       z_f,nz_f,z_g,ntimes,ntimes_flux,ntimes_smos,ntimes_sfc,&
       zs_f, ns_f, ntimes_soil, &
       times,times_flux,&
       times_smos,&
       times_sfc,times_soil, &
       nsplinetimes,splinetimes,&
       nsplinetimes_advection,splinetimes_advection,&
       nsplinetimes_flux,splinetimes_flux,&
       nsplinetimes_smos,splinetimes_smos,&
       nsplinetimes_sfc,splinetimes_sfc,&
       nsplinetimes_soil,splinetimes_soil, &
       pblh_ref,stepbl,lowlyr,ht,pblh,&
       th_phy,t_phy,moist,u_phy,v_phy,p_phy, p8w, tke_myj,&
       u_g_f,v_g_f,glw_f,gsw_f,precip_f,u_f,v_f,t_f,p_f,p8w_f,q_f,&
       tslb_f, smois_f, &
       th_upstream_x, th_upstream_y, &
       qv_upstream_x, qv_upstream_y, &
       qc_upstream_x, qc_upstream_y, &
       qr_upstream_x, qr_upstream_y, &
       qi_upstream_x, qi_upstream_y, &
       qg_upstream_x, qg_upstream_y, &
       u_upstream_x, u_upstream_y, &
       v_upstream_x, v_upstream_y, &
       tau_u, tau_v, &
       ts_f,qvs_f,ustar_f,hflux_f,qvflux_f,&
       z,z8w,dz8w, zs, num_soil, &
       ims,ime,jms,jme,kms,kme)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nz_f,ntimes,ntimes_flux, ntimes_smos,&
         ntimes_sfc, ns_f, ntimes_soil, &
         nsplinetimes,nsplinetimes_flux, nsplinetimes_smos,&
         nsplinetimes_sfc, nsplinetimes_advection, &
         nsplinetimes_soil, num_soil


    INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme, itime_f
    INTEGER, INTENT(out) :: stepbl

    REAL, DIMENSION(1:nz_f,ntimes), INTENT(IN) :: u_init_f,&
                                  v_init_f,t_init_f,qv_init_f, &
                                  qc_init_f,qr_init_f,qi_init_f, &
                                  qg_init_f

    REAL, DIMENSION(1:nz_f,ntimes), INTENT(IN) :: z_f,p_init_f, &
         t_init_upstream_x, t_init_upstream_y, &
         qv_init_upstream_x, qv_init_upstream_y, &
         u_init_upstream_x, u_init_upstream_y, &
         v_init_upstream_x, v_init_upstream_y, &
         tau_init_u, tau_init_v

    REAL, DIMENSION(1:ns_f), INTENT(IN)     :: zs_f
    REAL, DIMENSION(1:num_soil), INTENT(IN) :: zs

    REAL, DIMENSION(:,:), INTENT(IN)           :: &
         qc_init_upstream_x, qc_init_upstream_y, &
         qr_init_upstream_x, qr_init_upstream_y, &
         qi_init_upstream_x, qi_init_upstream_y, &
         qg_init_upstream_x, qg_init_upstream_y

    REAL, DIMENSION(ntimes), INTENT(IN) ::          &
         t2_init_upstream_x, t2_init_upstream_y,    &
         q2_init_upstream_x, q2_init_upstream_y,    &
         u10_init_upstream_x, u10_init_upstream_y,  &
         v10_init_upstream_x, v10_init_upstream_y,  &
         tau_init_u10, tau_init_v10

    REAL, DIMENSION(ntimes_flux) :: glw_init,gsw_init
    REAL, DIMENSION(ntimes_smos) :: precip_init 
    REAL, DIMENSION(:,:), INTENT(IN) :: tslb_init, smois_init
    REAL, DIMENSION(ntimes_sfc) :: ts_init,qvs_init,&
         ustar_init,hflux_init,qvflux_init

    REAL, DIMENSION(ntimes), INTENT(IN) :: times
    REAL, DIMENSION(ntimes_flux), INTENT(IN) :: times_flux
    REAL, DIMENSION(ntimes_smos), INTENT(IN) :: times_smos
    REAL, DIMENSION(ntimes_sfc), INTENT(IN) :: times_sfc
    REAL, DIMENSION(ntimes_soil), INTENT(IN) :: times_soil

    REAL, DIMENSION(nsplinetimes), INTENT(IN) :: splinetimes    
    REAL, DIMENSION(nsplinetimes_advection), INTENT(IN) :: splinetimes_advection    
    REAL, DIMENSION(nsplinetimes_flux), INTENT(IN) :: splinetimes_flux    
    REAL, DIMENSION(nsplinetimes_smos), INTENT(IN) :: splinetimes_smos    
    REAL, DIMENSION(nsplinetimes_sfc), INTENT(IN) :: splinetimes_sfc    
    REAL, DIMENSION(nsplinetimes_soil), INTENT(IN) :: splinetimes_soil    

    REAL, DIMENSION(kms:kme,nsplinetimes), INTENT(OUT) :: &
         u_g_f,v_g_f,u_f,v_f,t_f,p_f,p8w_f,q_f,           &
         th_upstream_x, th_upstream_y,                    &
         qv_upstream_x, qv_upstream_y,                    &
         u_upstream_x, u_upstream_y,                      &
         v_upstream_x, v_upstream_y,                      &
         tau_u, tau_v
    REAL, DIMENSION(:,:),                  INTENT(OUT) :: &
         qc_upstream_x, qc_upstream_y,                    &
         qr_upstream_x, qr_upstream_y,                    &
         qi_upstream_x, qi_upstream_y,                    &
         qg_upstream_x, qg_upstream_y

    REAL, DIMENSION(nsplinetimes_flux), INTENT(OUT) :: &
         glw_f,gsw_f

    REAL, DIMENSION(nsplinetimes_smos), INTENT(OUT) :: precip_f
    REAL, DIMENSION(num_soil, nsplinetimes_soil), INTENT(OUT) :: tslb_f,smois_f


    REAL, DIMENSION(nsplinetimes_sfc), INTENT(OUT) :: &
         ts_f,qvs_f,ustar_f,hflux_f,qvflux_f


    REAL, INTENT(IN), dimension(ntimes_smos) :: tsk_init,qsfc_init,&
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
    REAL, DIMENSION(ns_f+1) :: depth_var, sl_var
    REAL, DIMENSION(ntimes_flux) :: aatmp_flux,bbtmp_flux,cctmp_flux
    REAL, DIMENSION(kms:kme,ntimes) :: ptmp,p8wtmp,utmp,vtmp,ttmp,qtmp
    REAL, DIMENSION(num_soil,ntimes_soil) :: soil_tmp
    INTEGER :: i,j,k,kk,ii

! only need to interpolate first vertical profiles if we have UVG
    IF ( init_f_type == "OBS" .or. init_f_type == "SFC" ) THEN
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

!-------------------------------------------------------
! soil linear
! missing will be filled in so that no gradients exist
!-------------------------------------------------------

    IF ( force_soil ) THEN

      IF ( zs(1) .lt. zs_f(1)/100.0)  THEN
        PRINT*,'NEED SOIL INPUT AT SHALLOWER LAYER THAN TOP MODEL SOIL LAYER'
        STOP 'module_initialize: S/R initf'
      ENDIF
      DO i = 1, ntimes_soil
        DO k=1,ns_f
          depth_var(k)=zs_f(k)
          sl_var(k)=tslb_init(k,i)
          IF ( sl_var(k) <= 0.0 ) sl_var(k) = sl_var(k-1)
        ENDDO
        if ( maxval(depth_var) .gt. 10.0 ) depth_var(k) = depth_var(k)/100.0
        depth_var(ns_f+1) = 3.0
        sl_var(ns_f+1)    = sl_var(ns_f)

        ! vertical interp
        DO k=1,num_soil
          soil_tmp(k,i) = linear(depth_var,sl_var,ns_f+1,zs(k))
        ENDDO

      ENDDO ! ntimes_soil

      ! time interp
      IF ( ntimes_soil == nsplinetimes_soil ) THEN
        tslb_f = soil_tmp
      ELSE
        DO k = 1,num_soil
          DO i = 1,nsplinetimes_soil
            tslb_f(k,i) = &
            linear(times_soil,soil_tmp(k,:),ntimes_soil,splinetimes_soil(i))
          ENDDO
        ENDDO
      ENDIF

      DO i = 1, ntimes_soil
        DO k=1,ns_f
          depth_var(k)=zs_f(k)
          sl_var(k)=smois_init(k,i)
          IF ( sl_var(k) <= 0.0 ) sl_var(k) = sl_var(k-1)
        ENDDO
        if ( maxval(depth_var) .gt. 10.0 ) depth_var(k) = depth_var(k)/100.0
        depth_var(ns_f+1) = 3.0
        sl_var(ns_f+1)    = sl_var(ns_f)

        ! vertical interp
        DO k=1,num_soil
          soil_tmp(k,i) = linear(depth_var,sl_var,ns_f+1,zs(k))
        ENDDO

      ENDDO ! ntimes_soil

      ! time interp
      IF ( ntimes_soil == nsplinetimes_soil ) THEN
        smois_f = soil_tmp
      ELSE
        DO k = 1,num_soil
          DO i = 1,nsplinetimes_soil
            smois_f(k,i) = &
            linear(times_soil,soil_tmp(k,:),ntimes_soil,splinetimes_soil(i))
          ENDDO
        ENDDO
      ENDIF

    ENDIF !soil

!--------------------------------------------------------------------
! just skip everything else if we only want sfc obs sequences
    IF ( trim(init_f_type) == 'SFC' ) return
!--------------------------------------------------------------------
    
! estimate p2 and psfc using splines
! and do spline interpolation on p - don't bother with the 
! hydrostatic equation 

    p8wtmp = -9999.0
    ptmp = -9999.0
    DO i=1,ntimes_snd
       DO k=1,nz_f
          z_var(k)=z_f(k,i)
          var(k)=p_init_f(k,i)
       ENDDO

       CALL spline(nz_f,z_var,var,atmp,btmp,ctmp)
       
       DO k=1,nz
          ptmp(k,i)=&
               &seval(nz_f,z(1,k,1),z_var,var,atmp,btmp,ctmp)
          p8wtmp(k,i)=&
               &seval(nz_f,z8w(1,k,1),z_var,var,atmp,btmp,ctmp)
          IF (i==1) THEN
             p2=seval(nz_f,2.,z_var,var,atmp,btmp,ctmp)
             psfc=seval(nz_f,0.,z_var,var,atmp,btmp,ctmp)
          ENDIF
       ENDDO

    ENDDO

    DO k=1,nz
       DO i=1,nsplinetimes_snd
          p_f(k,i)=linear(times,ptmp(k,:),ntimes,splinetimes(i))
          p8w_f(k,i)=linear(times,p8wtmp(k,:),ntimes,splinetimes(i))
       ENDDO
    ENDDO

    DO i=ims,ime
       DO j=jms,jme
          DO k=kms,kme
             p_phy(i,k,j)= p_f(k,itime_f) 
             p8w(i,k,j)= p8w_f(k,itime_f) 
          ENDDO
       ENDDO
    ENDDO

! do linear interpolation t
    print*,'Initializing T'
    DO ii = 1, ntimes_snd

    near_sfc_num = 0
! first fills in lowest level with tsk, then if available subs in th2 

    z_var(1:2) = z_f(1:2,ii)
    IF ( tsk_init(ii) > -273 ) THEN
       near_sfc_num = near_sfc_num+1
       z_var(near_sfc_num)=0.
       var(near_sfc_num)=tsk_init(ii)
    ENDIF

    IF ( th2_init(ii) > -273 ) THEN
       IF ( z_var(near_sfc_num+1) > 2.0001 ) THEN
          near_sfc_num = near_sfc_num+1
          t2=th2_init(ii)
          var(near_sfc_num)=t2
          z_var(near_sfc_num)=2.
       ELSE
          t2=th2_init(ii)
          var(near_sfc_num+1)=t2
          z_var(near_sfc_num+1)=2.
       ENDIF
    ENDIF

    DO k=1,nz_f
       z_var(k+near_sfc_num)=z_f(k,ii)
       var(k+near_sfc_num)=t_init_f(k,ii)
    ENDDO

    DO k=kms,kme
             t_f(k,ii) = linear(z_var,var,nz_f+near_sfc_num,z(1,k,1))
    ENDDO

    ENDDO ! time loop

! replace missing value in th
    WHERE ( th_phy < -273 ) th_phy = -9999.0 

    DO i=ims,ime
       DO j=jms,jme
          DO k=kms,kme
             t_phy(i,k,j) = t_f(k,itime_f)
             th_phy(i,k,j)= t_f(k,itime_f) * (p1000mb/p_f(k,itime_f))**rcp
          ENDDO
       ENDDO
    ENDDO

    print*,'Initializing QV'
! do linear interpolation q
! use qsfc if we have it

    moist=0.
    DO ii = 1, ntimes_snd

    z_var(1:2) = z_f(1:2,ii)
    near_sfc_num = 0
    IF ( qsfc_init(ii) >= 0.0 ) THEN
       near_sfc_num = near_sfc_num+1
       z_var(near_sfc_num)=0.
       var(near_sfc_num)=qsfc_init(ii)
    ENDIF

    IF ( q2_init(ii) > 0.0 ) THEN 
       IF ( z_var(near_sfc_num+1) > 2.0001 ) THEN
          near_sfc_num = near_sfc_num+1
          var(near_sfc_num)=q2_init(ii)
          z_var(near_sfc_num)=2.
       ELSE
          var(near_sfc_num+1)=q2_init(ii)
          z_var(near_sfc_num+1)=2.
       ENDIF
    ENDIF

! for QV an extra check to ensure initialization.  if we have not 
! put in a value, use the first model level for q2
    if ( near_sfc_num == 0 ) then
      near_sfc_num = near_sfc_num + 1
      z_var(near_sfc_num)=2.
      var(near_sfc_num)=qv_init_f(near_sfc_num+1,ii)
    endif

    DO k=1,nz_f
       z_var(k+near_sfc_num)=z_f(k,ii)
       var(k+near_sfc_num)=qv_init_f(k,ii)
    ENDDO

    DO k=kms,kme
             q_f(k,ii) = linear(z_var,var,nz_f+near_sfc_num,z(1,k,1))
    ENDDO

    ENDDO ! time loop

    DO i=ims,ime
       DO j=jms,jme
          DO k=kms,kme
             moist(i,k,j,P_QV) = q_f(k,itime_f)
          ENDDO
       ENDDO
    ENDDO

    IF ( P_QC > 1 ) THEN
    print*,'Initializing QC'
! do linear interpolation q
    DO ii = 1, ntimes_snd

    DO k=1,nz_f
       z_var(k)=z_f(k,ii)
       var(k)=qc_init_f(k,ii)
    ENDDO

    q_f(1,ii) = 0.0 ! first level is 0 for all hydrometeors
    DO k=kms+1,kme
             q_f(k,ii) = linear(z_var,var,nz_f,z(1,k,1))
    ENDDO

    ENDDO ! time loop

    DO i=ims,ime
       DO j=jms,jme
          DO k=kms,kme
             moist(i,k,j,P_QC) = q_f(k,itime_f)
          ENDDO
       ENDDO
    ENDDO

    ENDIF ! init QC?

    IF ( P_QR > 1 ) THEN
    print*,'Initializing QR'
! do linear interpolation q
    DO ii = 1, ntimes_snd

    DO k=1,nz_f
       z_var(k)=z_f(k,ii)
       var(k)=qr_init_f(k,ii)
    ENDDO

    q_f(1,ii) = 0.0 ! first level is 0 for all hydrometeors
    DO k=kms+1,kme
             q_f(k,ii) = linear(z_var,var,nz_f,z(1,k,1))
    ENDDO

    ENDDO ! time loop

    DO i=ims,ime
       DO j=jms,jme
          DO k=kms,kme
             moist(i,k,j,P_QR) = q_f(k,itime_f)
          ENDDO
       ENDDO
    ENDDO

    ENDIF ! init QR?

    IF ( P_QI > 1 ) THEN
    print*,'Initializing QI'
! do linear interpolation q
    DO ii = 1, ntimes_snd

    DO k=1,nz_f
       z_var(k)=z_f(k,ii)
       var(k)=qi_init_f(k,ii)
    ENDDO

    q_f(1,ii) = 0.0 ! first level is 0 for all hydrometeors
    DO k=kms+1,kme
             q_f(k,ii) = linear(z_var,var,nz_f,z(1,k,1))
    ENDDO

    ENDDO ! time loop

    DO i=ims,ime
       DO j=jms,jme
          DO k=kms,kme
             moist(i,k,j,P_QI) = q_f(k,itime_f)
          ENDDO
       ENDDO
    ENDDO

    ENDIF ! init QI?

    IF ( P_QG > 1 ) THEN
    print*,'Initializing QG'
! do linear interpolation q
    DO ii = 1, ntimes_snd

    DO k=1,nz_f
       z_var(k)=z_f(k,ii)
       var(k)=qg_init_f(k,ii)
    ENDDO

    q_f(1,ii) = 0.0 ! first level is 0 for all hydrometeors
    DO k=kms+1,kme
             q_f(k,ii) = linear(z_var,var,nz_f,z(1,k,1))
    ENDDO

    ENDDO ! time loop

    DO i=ims,ime
       DO j=jms,jme
          DO k=kms,kme
             moist(i,k,j,P_QG) = q_f(k,itime_f)
          ENDDO
       ENDDO
    ENDDO

    ENDIF ! init QG?

    print*,'Initializing U'

    DO ii = 1, ntimes_snd

    z_var(1:2) = z_f(1:2,ii)
    near_sfc_num = 1
    z_var(1)=0.
    var(1)=0.

    IF ( u10_init(ii) > -500 ) THEN
       IF ( z_var(near_sfc_num+1) > 10.0001 ) THEN
          near_sfc_num = near_sfc_num+1
          var(near_sfc_num)=u10_init(ii)
          z_var(near_sfc_num)=10.0
       ELSE
          var(near_sfc_num+1)=u10_init(ii)
          z_var(near_sfc_num+1)=10.0
       ENDIF
    ENDIF
 
    skip_sfc_f = 0
    IF ( z_f(1,ii) == 2.0 ) skip_sfc_f = 1
 
    DO k=1,nz_f-skip_sfc_f
       z_var(k+near_sfc_num)=z_f(k+skip_sfc_f,ii)
       var(k+near_sfc_num)=u_init_f(k+skip_sfc_f,ii)
    ENDDO

    DO k=kms,kme
         u_f(k,ii) = linear(z_var,var,nz_f+near_sfc_num,z(1,k,1))
    ENDDO

    ENDDO ! time loop u

    DO i=ims,ime
       DO j=jms,jme
          DO k=kms,kme
             u_phy(i,k,j)=u_f(k,itime_f)
          ENDDO
       ENDDO
    ENDDO


! do linear interpolation v

    DO ii = 1, ntimes_snd

    z_var(1:2) = z_f(1:2,ii)
    near_sfc_num = 1
    z_var(1)=0.
    var(1)=0.

    IF ( v10_init(ii) > -500 ) THEN
       IF ( z_var(near_sfc_num+1) > 10.0001 ) THEN
          near_sfc_num = near_sfc_num+1
          var(near_sfc_num)=v10_init(ii)
          z_var(near_sfc_num)=10.0
       ELSE
          var(near_sfc_num+1)=v10_init(ii)
          z_var(near_sfc_num+1)=10.0
       ENDIF
    ENDIF

    skip_sfc_f = 0
    IF ( z_f(1,ii) == 2.0 ) skip_sfc_f = 1
 
    DO k=1,nz_f-skip_sfc_f
       z_var(k+near_sfc_num)=z_f(k+skip_sfc_f,ii)
       var(k+near_sfc_num)=v_init_f(k+skip_sfc_f,ii)
    ENDDO

    DO k=kms,kme
             v_f(k,ii) = linear(z_var,var,nz_f+near_sfc_num,z(1,k,1))
    ENDDO

    ENDDO ! time loop v

    DO i=ims,ime
       DO j=jms,jme
          DO k=kms,kme
             v_phy(i,k,j) = v_f(k,itime_f)
          ENDDO
       ENDDO
    ENDDO

! u_g
    kk=1
    DO WHILE(z(1,kk,1) < z_g)
       kk=kk+1
    ENDDO

    DO i=1,ntimes_snd
       DO k=1,nz_f
          z_var(k)=z_f(k,i)
          var(k)=u_init_f(k,i)
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
          var(k)=v_init_f(k,i)
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

       near_sfc_num = 0
       IF ( t2_init_upstream_x(i) > -273 ) THEN
         near_sfc_num = near_sfc_num+1
         t2=t2_init_upstream_x(i)
         var(near_sfc_num)=t2
         z_var(near_sfc_num)=2.
       ENDIF

       DO k=1,nz_f
             z_var(k+near_sfc_num)=z_f(k,i)
             var(k+near_sfc_num)=t_init_upstream_x(k,i)
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

       i = 1
       IF ( near_sfc_num == 1 ) THEN
         th_upstream_x(1,:)=th_upstream_x(1,:)*(p1000mb/psfc)**rcp
         i = 2
       ENDIF
       DO k=i,nz
         th_upstream_x(k,:)=th_upstream_x(k,:)*(p1000mb/p_f(k,:))**rcp
       ENDDO
       WHERE(th_upstream_x<-100) th_upstream_x = -9999.0 

! t upstream y
       DO i=1,ntimes_snd

          near_sfc_num = 0
          IF ( t2_init_upstream_y(i) > -273 ) THEN
            near_sfc_num = near_sfc_num+1
            t2=t2_init_upstream_y(i)
            var(near_sfc_num)=t2
            z_var(near_sfc_num)=2.
          ENDIF
          DO k=1,nz_f
             z_var(k+near_sfc_num)=z_f(k,i)
             var(k+near_sfc_num)=t_init_upstream_y(k,i)
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
       i = 1
       IF ( near_sfc_num == 1 ) THEN
         th_upstream_y(1,:)=th_upstream_y(1,:)*(p1000mb/psfc)**rcp
         i = 2
       ENDIF
       DO k=i,nz
         th_upstream_y(k,:)=th_upstream_y(k,:)*(p1000mb/p_f(k,:))**rcp
       ENDDO
       WHERE(th_upstream_y<-100) th_upstream_y = -9999.0 

    ENDIF
! moisture tendencies
    IF ( qv_advection ) then

! qv upstream x
       DO i=1,ntimes_snd
          near_sfc_num = 0
          IF ( q2_init_upstream_x(i) > 0.0 ) THEN
            near_sfc_num = near_sfc_num+1
            var(near_sfc_num)=q2_init_upstream_x(i)
            z_var(near_sfc_num)=2.
          ENDIF

          DO k=1,nz_f
             z_var(k+near_sfc_num)=z_f(k,i)
             var(k+near_sfc_num)=qv_init_upstream_x(k,i)
          ENDDO

          DO k=1,nz
             qtmp(k,i)=&
                     &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO
       ENDDO
       DO k=1,nz
          DO i=1,nsplinetimes_advection
             qv_upstream_x(k,i)=&
                  &linear(times,qtmp(k,:),ntimes,splinetimes_advection(i))
          ENDDO
       ENDDO
       WHERE(qv_upstream_x<-100) qv_upstream_x = -9999.0 

! qv upstream y
       DO i=1,ntimes_snd
          near_sfc_num = 0
          IF ( q2_init_upstream_y(i) > 0.0 ) THEN
            near_sfc_num = near_sfc_num+1
            var(near_sfc_num)=q2_init_upstream_y(i)
            z_var(near_sfc_num)=2.
          ENDIF

          DO k=1,nz_f
             z_var(k+near_sfc_num)=z_f(k,i)
             var(k+near_sfc_num)=qv_init_upstream_y(k,i)
          ENDDO

          DO k=1,nz
             qtmp(k,i)=&
                     &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO
       ENDDO
       DO k=1,nz
          DO i=1,nsplinetimes_advection
             qv_upstream_y(k,i)=&
                  &linear(times,qtmp(k,:),ntimes,splinetimes_advection(i))
          ENDDO
       ENDDO
       WHERE(qv_upstream_y<-100) qv_upstream_y = -9999.0 

       IF ( P_QC > 1 ) THEN

! qc upstream x
       DO i=1,ntimes_snd

          DO k=1,nz_f
             z_var(k)=z_f(k,i)
             var(k)=qc_init_upstream_x(k,i)
          ENDDO

          qtmp(1,i) = 0.0 ! all hydrometeors are 0 on first level
          DO k=2,nz
             qtmp(k,i)=&
                     &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO
       ENDDO
       DO k=1,nz
          DO i=1,nsplinetimes_advection
             qc_upstream_x(k,i)=&
                  &linear(times,qtmp(k,:),ntimes,splinetimes_advection(i))
          ENDDO
       ENDDO
       WHERE(qc_upstream_x<-100) qc_upstream_x = -9999.0 

! qc upstream y
       DO i=1,ntimes_snd

          DO k=1,nz_f
             z_var(k)=z_f(k,i)
             var(k)=qc_init_upstream_y(k,i)
          ENDDO

          qtmp(1,i) = 0.0 ! all hydrometeors are 0 on first level
          DO k=2,nz
             qtmp(k,i)=&
                     &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO
       ENDDO

       DO k=1,nz
          DO i=1,nsplinetimes_advection
             qc_upstream_y(k,i)=&
                  &linear(times,qtmp(k,:),ntimes,splinetimes_advection(i))
          ENDDO
       ENDDO
       WHERE(qc_upstream_y<-100) qc_upstream_y = -9999.0 
 
       ENDIF ! QC ?

       IF ( P_QR > 1 ) THEN

! qr upstream x
       DO i=1,ntimes_snd

          DO k=1,nz_f
             z_var(k)=z_f(k,i)
             var(k)=qr_init_upstream_x(k,i)
          ENDDO

          qtmp(1,i) = 0.0 ! all hydrometeors are 0 on first level
          DO k=2,nz
             qtmp(k,i)=&
                     &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO
       ENDDO
       DO k=1,nz
          DO i=1,nsplinetimes_advection
             qr_upstream_x(k,i)=&
                  &linear(times,qtmp(k,:),ntimes,splinetimes_advection(i))
          ENDDO
       ENDDO
       WHERE(qr_upstream_x<-100) qr_upstream_x = -9999.0 

! qr upstream y
       DO i=1,ntimes_snd

          DO k=1,nz_f
             z_var(k)=z_f(k,i)
             var(k)=qr_init_upstream_y(k,i)
          ENDDO

          qtmp(1,i) = 0.0 ! all hydrometeors are 0 on first level
          DO k=2,nz
             qtmp(k,i)=&
                     &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO
       ENDDO
       DO k=1,nz
          DO i=1,nsplinetimes_advection
             qr_upstream_y(k,i)=&
                  &linear(times,qtmp(k,:),ntimes,splinetimes_advection(i))
          ENDDO
       ENDDO
       WHERE(qr_upstream_y<-100) qr_upstream_y = -9999.0 
 
       ENDIF ! QR ?

       IF ( P_QI > 1 ) THEN

! qi upstream x
       DO i=1,ntimes_snd

          DO k=1,nz_f
             z_var(k)=z_f(k,i)
             var(k)=qi_init_upstream_x(k,i)
          ENDDO

          qtmp(1,i) = 0.0 ! all hydrometeors are 0 on first level
          DO k=2,nz
             qtmp(k,i)=&
                     &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO
       ENDDO
       DO k=1,nz
          DO i=1,nsplinetimes_advection
             qi_upstream_x(k,i)=&
                  &linear(times,qtmp(k,:),ntimes,splinetimes_advection(i))
          ENDDO
       ENDDO
       WHERE(qi_upstream_x<-100) qi_upstream_x = -9999.0 

! qi upstream y
       DO i=1,ntimes_snd

          DO k=1,nz_f
             z_var(k)=z_f(k,i)
             var(k)=qi_init_upstream_y(k,i)
          ENDDO

          qtmp(1,i) = 0.0 ! all hydrometeors are 0 on first level
          DO k=2,nz
             qtmp(k,i)=&
                     &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO
       ENDDO
       DO k=1,nz
          DO i=1,nsplinetimes_advection
             qi_upstream_y(k,i)=&
                  &linear(times,qtmp(k,:),ntimes,splinetimes_advection(i))
          ENDDO
       ENDDO
       WHERE(qi_upstream_y<-100) qi_upstream_y = -9999.0 
 
       ENDIF ! QI ?

       IF ( P_QG > 1 ) THEN

! qg upstream x
       DO i=1,ntimes_snd

          DO k=1,nz_f
             z_var(k)=z_f(k,i)
             var(k)=qg_init_upstream_x(k,i)
          ENDDO

          qtmp(1,i) = 0.0 ! all hydrometeors are 0 on first level
          DO k=2,nz
             qtmp(k,i)=&
                     &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO
       ENDDO
       DO k=1,nz
          DO i=1,nsplinetimes_advection
             qg_upstream_x(k,i)=&
                  &linear(times,qtmp(k,:),ntimes,splinetimes_advection(i))
          ENDDO
       ENDDO
       WHERE(qg_upstream_x<-100) qg_upstream_x = -9999.0 

! qg upstream y
       DO i=1,ntimes_snd

          DO k=1,nz_f
             z_var(k)=z_f(k,i)
             var(k)=qg_init_upstream_y(k,i)
          ENDDO

          qtmp(1,i) = 0.0 ! all hydrometeors are 0 on first level
          DO k=2,nz
             qtmp(k,i)=&
                     &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO
       ENDDO
       DO k=1,nz
          DO i=1,nsplinetimes_advection
             qg_upstream_y(k,i)=&
                  &linear(times,qtmp(k,:),ntimes,splinetimes_advection(i))
          ENDDO
       ENDDO
       WHERE(qg_upstream_y<-100) qg_upstream_y = -9999.0 
 
       ENDIF ! QG ?

    ENDIF ! q_advection

! wind tendencies ( need for taus )
    IF ( t_advection .or. u_advection .or. qv_advection ) then

! u upstream x
       DO i=1,ntimes_snd
          near_sfc_num = 0
          IF ( u10_init_upstream_x(i) > -999.0 ) THEN
            near_sfc_num = near_sfc_num+1
            var(near_sfc_num)=u10_init_upstream_x(i)
            z_var(near_sfc_num)=10.
          ENDIF

          DO k=1,nz_f
             z_var(k+near_sfc_num)=z_f(k,i)
             var(k+near_sfc_num)=u_init_upstream_x(k,i)
          ENDDO

          DO k=1,nz
             utmp(k,i)=&
                     &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO
       ENDDO

       DO k=1,nz
          DO i=1,nsplinetimes_advection
             u_upstream_x(k,i)=&
                  &linear(times,utmp(k,:),ntimes,splinetimes_advection(i))
          ENDDO
       ENDDO

! u upstream y
       DO i=1,ntimes_snd
          near_sfc_num = 0
          IF ( u10_init_upstream_y(i) > -999.0 ) THEN
            near_sfc_num = near_sfc_num+1
            var(near_sfc_num)=u10_init_upstream_y(i)
            z_var(near_sfc_num)=10.
          ENDIF

          DO k=1,nz_f
             z_var(k+near_sfc_num)=z_f(k,i)
             var(k+near_sfc_num)=u_init_upstream_y(k,i)
          ENDDO

          DO k=1,nz
             utmp(k,i)=&
                     &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO
       ENDDO
       DO k=1,nz
          DO i=1,nsplinetimes_advection
             u_upstream_y(k,i)=&
                  &linear(times,utmp(k,:),ntimes,splinetimes_advection(i))
          ENDDO
       ENDDO

! v upstream x
       DO i=1,ntimes_snd
          near_sfc_num = 0
          IF ( v10_init_upstream_x(i) > -999.0 ) THEN
            near_sfc_num = near_sfc_num+1
            var(near_sfc_num)=v10_init_upstream_x(i)
            z_var(near_sfc_num)=10.
          ENDIF

          DO k=1,nz_f
             z_var(k+near_sfc_num)=z_f(k,i)
             var(k+near_sfc_num)=v_init_upstream_x(k,i)
          ENDDO

          DO k=1,nz
             vtmp(k,i)=&
                     &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO
       ENDDO
       DO k=1,nz
          DO i=1,nsplinetimes_advection
             v_upstream_x(k,i)=&
                  linear(times,vtmp(k,:),ntimes,splinetimes_advection(i))
          ENDDO
       ENDDO

! v upstream y
       DO i=1,ntimes_snd
          near_sfc_num = 0
          IF ( v10_init_upstream_y(i) > -999.0 ) THEN
            near_sfc_num = near_sfc_num+1
            var(near_sfc_num)=v10_init_upstream_y(i)
            z_var(near_sfc_num)=10.
          ENDIF

          DO k=1,nz_f
             z_var(k+near_sfc_num)=z_f(k,i)
             var(k+near_sfc_num)=v_init_upstream_y(k,i)
          ENDDO

          DO k=1,nz
             vtmp(k,i)=&
                     &linear(z_var,var,nz_f,z(1,k,1))
          ENDDO
       ENDDO
       DO k=1,nz
          DO i=1,nsplinetimes_advection
             v_upstream_y(k,i)=&
                  linear(times,vtmp(k,:),ntimes,splinetimes_advection(i))
          ENDDO
       ENDDO

    ENDIF

! advecting winds
    IF ( t_advection .or. qv_advection .or. u_advection ) THEN
! tau u 
       DO i=1,ntimes_snd
          near_sfc_num = 0
          IF ( tau_init_u10(i) /= -9999.0 ) THEN
            near_sfc_num = near_sfc_num+1
            var(near_sfc_num)=tau_init_u10(i)
            z_var(near_sfc_num)=10.
          ENDIF

          DO k=1,nz_f
             z_var(k+near_sfc_num)=z_f(k,i)
             var(k+near_sfc_num)=tau_init_u(k,i)
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
          near_sfc_num = 0
          IF ( tau_init_v10(i) /= -9999.0 ) THEN
            near_sfc_num = near_sfc_num+1
            var(near_sfc_num)=tau_init_v10(i)
            z_var(near_sfc_num)=10.
          ENDIF

          DO k=1,nz_f
             z_var(k+near_sfc_num)=z_f(k,i)
             var(k+near_sfc_num)=tau_init_v(k,i)
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

    ENDIF ! advecting winds

 
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


