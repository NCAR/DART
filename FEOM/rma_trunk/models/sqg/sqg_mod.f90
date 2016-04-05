! This code is not protected by the DART copyright agreement.
! DART $Id$

!========================================================================
! Uniform-PV two-surface QG+1 model in spectral form (Hakim 2000)
!========================================================================

MODULE sqg_mod

    USE netcdf
    USE spectral_mod

CONTAINS

!========================================================================
SUBROUTINE sqg_main()

    implicit none

    ! spectral values
    complex, dimension(2*kmax,2*lmax) :: thbB,thspB,thspB1,thspB2
    complex, dimension(2*kmax,2*lmax) :: thbT,thspT,thspT1,thspT2
    complex, dimension(2*kmax,2*lmax) :: sB,sBold
    ! spectral work arrays
    real,    dimension(2*kmax,2*lmax) :: thxyB,thxyT
    ! grid point values
    real,    dimension(mmax,nmax)     :: thxB,thyB,uB,vB
    real,    dimension(mmax,nmax)     :: thxT,thyT,uT,vT
    ! tendencies
    real,    dimension(mmax,nmax)     :: lap
    real,    dimension(mmax,nmax)     :: hx,hy,hu,hv 
    real,    dimension(mmax,nmax)     :: tthB,tthT
    complex, dimension(mmax,nmax)     :: tthspB,tthspT
    ! basic state
    real,    dimension(mmax,nmax)     :: thbyB,thbyT
    real,    dimension(mmax,nmax)     :: ulinB,ulinT
    ! input and output
    character(len=64), parameter      :: inpfile = 'sqgInput.nc'
    character(len=64), parameter      :: outfile = 'sqgOutput.nc'
    character(len=64), parameter      :: rstfile = 'sqgRestart.nc'
    ! misc
    real                              :: dco,lam
    real                              :: cxB,cyB,cxT,cyT
    integer                           :: itime
    logical                           :: first,bot,top
    real,    dimension(mmax,nmax), parameter :: Rblank = 0.0

    if (verbose .gt. 0) print *,'Running with Rossby number: ', Ross

    ! initialize diffusion: 
    call diffusion(dco)

    ! initialize derivative operators:
    call d_setup()

    ! initialize theta fields:
    call init(inpfile,thxyB,thxyT)

    ! initialize base-state jet:
    if ( hw ) then 
        call init_jet(thbB,thbT,thbyB,thbyT,ulinB,ulinT,lam)
        ulinB = ulinB + 0; ulinT = ulinT - 0*H*lam     ! barotropic wind (Ross = 0!)
    else
        thbB = 0;thbT = 0;thbyB=0.;thbyT=0.;ulinB=0.;ulinT=0.;lam=0.
    endif

    if (verbose .gt. 1) print*,'lam = ',lam
    if (verbose .gt. 1) print*,'extrema ulinT = ',maxval(ulinT),minval(ulinT)

    ! write basic state to disk (linear shear:  ON : 1, OFF : 0)
    !call write_diag('base.nc',0,thbB,thbT)
    !call dump(thbB,thbT,.TRUE.,lam,1,'base.nc')

    ! create output file
    call write_diag(outfile,0,Rblank,Rblank)

    ! advection flags
    if     (model .eq. 0) then     ! 2D
    top = .FALSE.; bot = .TRUE.
    elseif (model .eq. 1) then     ! 2sQG
    top = .TRUE.; bot = .TRUE.
    elseif (model .eq. 2) then     ! tropo sQG
    top = .TRUE.; bot = .FALSE.
    elseif (model .eq. 3) then     ! surface sQG
    top = .FALSE.; bot = .TRUE.
    elseif (model .eq. 4) then     ! tropo HsQG
    top = .TRUE.; bot = .FALSE.
    endif

    if (verbose .gt. 0) print*,'max BOTTOM initial value=',maxval(abs(thxyB)), bot
    if (verbose .gt. 0) print*,'max TOPPOM initial value=',maxval(abs(thxyT)), top

    ! map into spectral space at the same resolution:
    call xy_to_sp(cmplx(thxyB,0.),thspB,2*kmax,2*lmax,kmax,lmax)
    call xy_to_sp(cmplx(thxyT,0.),thspT,2*kmax,2*lmax,kmax,lmax)

    ! initialize terrain:
    hu = 0.; hv = 0.; hx = 0.; hy = 0.
    if (iterr .and. Ross .eq. 0) then 
        if (verbose .gt. 1)   print*,'initializing terrain'
        call terrain(hx,hy,hu,hv)
        call invert(thspB,0.0*thspT,thxB,thxT,thyB,thyT,vB,vT,uB,uT, &
                    thbB,thbT,thbyB,thbyT,ulinB,ulinT,               &
                    .TRUE.,.TRUE.,.TRUE.,0.0,sB,sBold,lap)
        hu = uT; hv = vT
        if (verbose .gt. 1) print*,'extrema hx = ',maxval(hx),minval(hx)
        if (verbose .gt. 1) print*,'extrema hy = ',maxval(hy),minval(hy)
        if (verbose .gt. 1) print*,'extrema hu = ',maxval(hu),minval(hu)
        if (verbose .gt. 1) print*,'extrema hv = ',maxval(hv),minval(hv)
    endif

    first = .TRUE.

    ! option to zero k = 1 modes
    !thspB(2,:) = 0.; thspB(2*kmax,:) = 0.
    !thspT(2,:) = 0.; thspT(2*kmax,:) = 0.

    !!!!!!!!!!!!!!!!!!!!
    ! BEGIN: Time loop !
    !!!!!!!!!!!!!!!!!!!!
    do itime = 1, ntims

        ! save old streamFUNCTION for Ekman calculation.
        sBold = sB; if (first) sBold = 0.0

        ! option to zero l = 0 modes, which are nonlinear solutions
        !thspB(:,1) = 0; thspT(:,1) = 0;

        ! this is for growing the most unstable mode
        if ( (grow) .and. (cyT .gt. 0.001) ) then
            if (verbose .gt. 1) print*,'RESCALING'
            thspB = thspB/2.0 ; thspT = thspT/2.0 ; cyT = 0.0
        endif

        ! Invert theta for streamFUNCTION; compute gradients for advection:
        call invert(thspB,thspT,thxB,thxT,thyB,thyT,vB,vT,uB,uT, &
                    thbB,thbT,thbyB,thbyT,ulinB,ulinT,&
                    first,bot,top,lam,sB,sBold,lap)

        ! option to compute potential enstrophy norm and growth rate:
        if ( inorm ) call norm(thspB,thspT,itime)

        ! write data to file for plots
        if ( mod((itime-1),iplot) .eq. 0 )  &
            call dump(thspB+thbB,thspT+thbT,.TRUE.,lam,itime,outfile)

        ! spectral advection:
        if (bot) call advect(uB,   vB,   thxB,thyB,thbyB,hx,    hy,    ulinB,tthB,lam,lap  )
        if (top) call advect(uT+hu,vT+hv,thxT,thyT,thbyT,Rblank,Rblank,ulinT,tthT,lam,Rblank)

        ! Write Courant numbers to stdout
        if ( verbose .gt. 0 ) then
            if (mod((itime-1),10) .eq. 0) then 
                cxB = maxval(abs(uB+ulinB))*dt/(XL/real(2*kmax))
                cyB = maxval(abs(vB))*dt/(YL/real(2*lmax))
                cxT = maxval(abs(uT+ulinT))*dt/(XL/real(2*kmax))
                cyT = maxval(abs(vT))*dt/(YL/real(2*lmax))
                write(*,'(A23,F10.3,4F8.3)') 'time,cxB,cyB,cxT,cyT = ', &
                                  real(itime-1)*dt,cxB,cyB,cxT,cyT
            endif
        endif

        ! FFT back to spectral space:
        if (bot) then 
            tthspB = cmplx(tthB,0.); call ft_2d(tthspB,mmax,nmax,-1)
        endif
        if (top) then 
            tthspT = cmplx(tthT,0.); call ft_2d(tthspT,mmax,nmax,-1)
        endif

        ! advance one time step with explicit (hyper-)diffusion:
        if ( first ) then
            thspB1 = 0.0 ; thspB2 = 0.0
            thspT1 = 0.0 ; thspT2 = 0.0
        endif
        if (bot) call tadv(thspB,tthspB,thspB1,thspB2,dco,first)
        if (top) call tadv(thspT,tthspT,thspT1,thspT2,dco,first)

        thspB(kmax,:) = 0.; thspB(lmax,:) = 0.

        first = .FALSE.

    end do ! itime
    !!!!!!!!!!!!!!!!!!
    ! END: Time loop !
    !!!!!!!!!!!!!!!!!!

    ! write (final + 1) time to disk for future restart:
    call write_diag(rstfile,0,thxyB,thxyT)
    call dump(thspB,thspT,.FALSE.,lam,1,rstfile)

    return
END SUBROUTINE sqg_main
!========================================================================

!========================================================================
SUBROUTINE invert(ithspB,ithspT,thxBr,thxTr,thyBr,thyTr,vBr,vTr,uBr,uTr, &
                  thbB,thbT,thbyB,thbyT,ulinB,ulinT, &
                  first,bot,top,lam,sB,sBold,sblre)
! Invert PV and transform to a larger grid for de-aliasing.

    implicit none

    complex, dimension(2*kmax,2*lmax), intent(in)  :: ithspB,ithspT
    complex, dimension(2*kmax,2*lmax), intent(in)  :: thbB,thbT
    complex, dimension(2*kmax,2*lmax), intent(in)  :: sBold
    real,    dimension(mmax,nmax),     intent(in)  :: thbyB,thbyT
    real,    dimension(mmax,nmax),     intent(in)  :: ulinB,ulinT
    real,                              intent(in)  :: lam
    logical,                           intent(in)  :: first,bot,top
    real,    dimension(mmax,nmax),     intent(out) :: thxBr,thxTr,thyBr,thyTr
    real,    dimension(mmax,nmax),     intent(out) :: uBr,uTr,vBr,vTr
    real,    dimension(mmax,nmax),     intent(out) :: sblre
    complex, dimension(2*kmax,2*lmax), intent(out) :: sB

    complex, dimension(2*kmax,2*lmax) :: temps
    complex, dimension(2*kmax,2*lmax) :: thspB,thspT
    complex, dimension(2*kmax,2*lmax) :: tempspB,tempspT
    complex, dimension(2*kmax,2*lmax) :: tempxyB,tempxyT
    complex, dimension(2*kmax,2*lmax) :: htempB,htempT
    complex, dimension(2*kmax,2*lmax) :: u1spB,u1spT,v1spB,v1spT
    complex, dimension(mmax,nmax)     :: thxB,thxT,thyB,thyT
    complex, dimension(mmax,nmax)     :: uB,uT,vB,vT
    complex, dimension(mmax,nmax)     :: szspB,szspT,szzspB,szzspT
    complex, dimension(mmax,nmax)     :: u1B,v1B,u1T,v1T
    complex, dimension(mmax,nmax)     :: temp

    logical :: correct

    ! flag for next-order corrections
    correct = .FALSE.
    if (Ross .gt. 1.e-5) correct = .TRUE. 

    thspB = ithspB; thspT = ithspT ! local copies of spectral boundary theta
!!!!!!!!!!!!!!!!!!!!!!!!!!! Grad Theta !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    thxT = 0.; thyT = 0.; vT = 0.; uT = 0.
    thxB = 0.; thyB = 0.; vB = 0.; uB = 0.

    call d_s2b(thspB,thxB,1,d_oper%dx) ! x derivative: small to big domain.
    call d_s2b(thspB,thyB,1,d_oper%dy) ! y derivative: small to big domain.
    call d_s2b(thspT,thxT,1,d_oper%dx) ! x derivative: small to big domain.
    call d_s2b(thspT,thyT,1,d_oper%dy) ! y derivative: small to big domain.

!!!!!!!!!!!!!!!!!!!!!!!!!!! sQG Bottom !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (bot .or. correct) then ! velocity contribution from lower boundary
        call d_b2b( thxB,vB,1,-d_oper%iz) 
        call d_b2b(-thyB,uB,1,-d_oper%iz)
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!! sQG Toppom !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (top .or. correct) then ! velocity contribution from upper boundary
        call d_b2b( thxT,vT,1,d_oper%iz) ! z integral: big to big domain.
        call d_b2b(-thyT,uT,1,d_oper%iz) ! z integral: big to big domain.
    endif
!!!!!!!!!!!!!!!!!!!!!!!! sQG Cross Boundary !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (bot .and. top .or. correct) then ! lower velocity from upper boundary
        call d_b2b( thxT,temp,1,d_oper%izo); vB = vB + temp;
        call d_b2b(-thyT,temp,1,d_oper%izo); uB = uB + temp;
    endif
    if (bot .and. top .or. correct) then ! upper velocity from lower boundary
        call d_b2b( thxB,temp,1,-d_oper%izo); vT = vT + temp;
        call d_b2b(-thyB,temp,1,-d_oper%izo); uT = uT + temp;
    endif
    ! FFT back to xy grid.
    if (bot .or. correct) then 
        call ft_2d(thxB,mmax,nmax,1) ! t_x
        call ft_2d(vB,  mmax,nmax,1) ! v0
        call ft_2d(uB,  mmax,nmax,1) ! u0
        call ft_2d(thyB,mmax,nmax,1) ! t_y
    endif
    if (top .or. correct) then 
        call ft_2d(thxT,mmax,nmax,1) ! t_x
        call ft_2d(thyT,mmax,nmax,1) ! t_y
        call ft_2d(vT,  mmax,nmax,1) ! v0
        call ft_2d(uT,  mmax,nmax,1) ! u0         
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!! sQG +1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (correct) then
        szzspB=0.;szspB=0.;u1spB=0.;v1spB=0. ! these _MUST_ be zeroed.
        szzspT=0.;szspT=0.;u1spT=0.;v1spT=0. 
        u1spB=0.;u1spT=0.;v1spB=0.;v1spT=0.
        u1B=0.;u1T=0.;v1B=0.;v1T=0.

        ! basic-state terms: add periodic parts here; linear parts are handled last.
        if (hw) then 
            thspB = thspB + thbB; thspT = thspT + thbT
            uB = uB + ulinB; uT = uT + ulinT - (lam*H) ! don't double count!
            thyB = thyB + thbyB; thyT = thyT + thbyT
        endif

        ! big-grid theta (szsp) and theta_z (szzsp):
        call d_s2b(thspB,szspB,1,d_oper%Id);  call ft_2d(szspB,mmax,nmax,1)
        call d_s2b(thspT,szspT,1,d_oper%Id);  call ft_2d(szspT,mmax,nmax,1)

        call d_s2b(thspB,szzspB,1,-d_oper%dz); call d_s2b(thspT,temp,1, d_oper%dzo)
        szzspB = szzspB + temp
        call ft_2d(szzspB,mmax,nmax,1)
        call d_s2b(thspT,szzspT,1, d_oper%dz); call d_s2b(thspB,temp,1,-d_oper%dzo)
        szzspT = szzspT + temp
        call ft_2d(szzspT,mmax,nmax,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Homogeneous contributions (Note: xy_to_sp dealiases; no d_b2s)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Compute u1 = -F1_z contributions (nonlinearities)
        temp = -uB*szspB
        call xy_to_sp(temp,tempspB,mmax,nmax,kmax,lmax)
        temp = -uT*szspT
        call xy_to_sp(temp,tempspT,mmax,nmax,kmax,lmax)
        ! u1 bottom
        call d_s2s(tempspB,htempB,1, d_oper%dz)
        call d_s2s(tempspT,htempT,1,-d_oper%dzo)
        u1spB = - (htempB + htempT) ! -F1_z
        ! u1 toppom
        call d_s2s(tempspB,htempB,1, d_oper%dzo)
        call d_s2s(tempspT,htempT,1,-d_oper%dz)
        u1spT = - (htempB + htempT) ! -F1_z

        ! Compute v1 = -G1_z contributions (nonlinearities)
        temp = vB*szspB
        call xy_to_sp(temp,tempspB,mmax,nmax,kmax,lmax)
        temp = vT*szspT
        call xy_to_sp(temp,tempspT,mmax,nmax,kmax,lmax)
        ! v1 bottom
        call d_s2s(tempspB,htempB,1,-d_oper%dz)
        call d_s2s(tempspT,htempT,1, d_oper%dzo)
        v1spB = - (htempB + htempT) ! -G1_z
        ! v1 toppom
        call d_s2s(tempspB,htempB,1,-d_oper%dzo)
        call d_s2s(tempspT,htempT,1, d_oper%dz)
        v1spT = - (htempB + htempT) ! -G1_z

        ! Compute Phi1 contributions (nonlinearities)
        temp = szzspB*szspB
        call xy_to_sp(temp,tempspB,mmax,nmax,kmax,lmax)
        temp = szzspT*szspT
        call xy_to_sp(temp,tempspT,mmax,nmax,kmax,lmax)
        ! u1 bottom
        call d_s2s(tempspB,htempB,1, d_oper%iz  * d_oper%dy)
        call d_s2s(tempspT,htempT,1,-d_oper%izo * d_oper%dy)
        u1spB = u1spB - (htempB + htempT) ! add -P1_y
        ! u1 toppom
        call d_s2s(tempspB,htempB,1, d_oper%izo * d_oper%dy)
        call d_s2s(tempspT,htempT,1,-d_oper%iz  * d_oper%dy)
        u1spT = u1spT - (htempB + htempT) ! add -P1_y
        ! v1 bottom
        call d_s2s(tempspB,htempB,1, d_oper%iz  * d_oper%dx)
        call d_s2s(tempspT,htempT,1,-d_oper%izo * d_oper%dx)
        v1spB = v1spB + (htempB + htempT) ! add P1_x
        ! v1 toppom
        call d_s2s(tempspB,htempB,1, d_oper%izo * d_oper%dx)
        call d_s2s(tempspT,htempT,1,-d_oper%iz  * d_oper%dx)
        v1spT = v1spT + (htempB + htempT) ! add P1_x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Particular contributions 
     ! Notes: \Phi contrib gives 2 factor!!!
     !        Dealiasing is crucial here (no shortcuts, DJM!)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! u1 bottom
        temp = (-2.*thyB*szspB) + (uB*szzspB)
        call xy_to_sp(temp,tempspB,mmax,nmax,kmax,lmax)
        u1spB = u1spB + tempspB !add (-P1_y - F1_z)
        ! u1 toppom
        temp = (-2.*thyT*szspT) + (uT*szzspT)
        call xy_to_sp(temp,tempspT,mmax,nmax,kmax,lmax)
        u1spT = u1spT + tempspT !add (-P1_y - F1_z)
        ! v1 bottom
        temp = (2.*thxB*szspB) + (vB*szzspB)
        call xy_to_sp(temp,tempspB,mmax,nmax,kmax,lmax)
        v1spB = v1spB + tempspB !add (P1_x - G1_z)
        ! v1 toppom
        temp = (2.*thxT*szspT) + (vT*szzspT)
        call xy_to_sp(temp,tempspT,mmax,nmax,kmax,lmax)
        v1spT = v1spT + tempspT !add (P1_x - G1_z)

        ! linear-shear contributions (note: recycling var names)
        if (hw) then 
            call d_s2s(thspB,tempspB,1,-d_oper%izo*d_oper%dx*d_oper%dx) !phi_xx^B
            call d_s2s(thspT,tempspT,1, d_oper%izo*d_oper%dx*d_oper%dx) !phi_xx^T
            call d_s2s(thspB,htempB, 1,-d_oper%izo*d_oper%dy*d_oper%dy) !phi_yy^B
            call d_s2s(thspT,htempT, 1, d_oper%izo*d_oper%dy*d_oper%dy) !phi_yy^T
            call d_s2s(thspB,tempxyB,1,-d_oper%izo*d_oper%dx*d_oper%dy) !phi_xy^B
            call d_s2s(thspT,tempxyT,1, d_oper%izo*d_oper%dx*d_oper%dy) !phi_xy^T

            !u1spB = u1spB - lam*(H*(htempT - tempspT) + (thspT - thspB) )
            !u1spT = u1spT + lam*(H*(htempB - tempspB) + (thspT - thspB) )
            ! as per DM email 12/24/02:
            u1spB = u1spB - lam*(H*(htempT - tempspT) - thspB )
            u1spT = u1spT + lam*(H*(htempB - tempspB) + thspT )
            v1spB = v1spB + lam*2*H*tempxyT
            v1spT = v1spT - lam*2*H*tempxyB
        endif

        ! map from small to large spectral array
        call d_s2b(u1spB,u1B,1,d_oper%Id); call d_s2b(v1spB,v1B,1,d_oper%Id)
        call d_s2b(u1spT,u1T,1,d_oper%Id); call d_s2b(v1spT,v1T,1,d_oper%Id)

        ! FFT u1 and v1 to grid point space
        call ft_2d(u1B,mmax,nmax,1); call ft_2d(v1B,mmax,nmax,1)
        call ft_2d(u1T,mmax,nmax,1); call ft_2d(v1T,mmax,nmax,1)

        ! remove periodic base state (added above):
        uB = uB - ulinB; uT = uT - ulinT + (lam*H)
        thyB = thyB - thbyB; thyT = thyT - thbyT
    endif ! correct

    ! return u = u_0 + Ro * u_1; v = v_0 + Ro * v_1
    ! (look into using the f90 "transfer" FUNCTION to return real parts...)

    thxBr = real(thxB);          thyBr = real(thyB)
    thxTr = real(thxT);          thyTr = real(thyT)
    uBr   = real(uB + Ross*u1B); uTr   = real(uT + Ross*u1T)
    vBr   = real(vB + Ross*v1B); vTr   = real(vT + Ross*v1T)

    ! Ekman layer calculations (need streamFUNCTION and Laplacian.
    if (gamma .gt. 0) then
        ! reset these since they were changed for next-order calculations
        thspB = ithspB; thspT = ithspT ! local copies of spectral boundary theta
        sb = 0. ! surface O(1) streamFUNCTION
        call d_s2s(thspB,sB,   1,-d_oper%iz)  ! bottom contribution
        call d_s2s(thspT,temps,1, d_oper%izo) ! toppom contribution
        sB = sB + temps;
        ! now compute Laplacian from previous time-step's sB
        temp = 0.
        call d_s2b(sBold,temp,1,d_oper%dx*d_oper%dx + d_oper%dy*d_oper%dy)
        call ft_2d(temp,mmax,nmax,1)
        sblre = real(temp)
    endif

    return
END SUBROUTINE invert
!========================================================================

!========================================================================
SUBROUTINE xy_to_sp(xy,sp,mx,ny,km,lm)
! Map an (x,y) array onto a _smaller_ spectral array.
! Input: xy(mx,ny) --- a grid point array.
! Output: sp(2*km,2*lm) --- a spectral array.

    implicit none

    integer,                       intent(in)  :: km,lm,mx,ny
    complex, dimension(:,:),       intent(in)  :: xy
    complex, dimension(2*km,2*lm), intent(out) :: sp

    complex, dimension(mx,ny) :: copy
    integer                   :: kmp1,lmp1,k,l,k2,l2,kk,ll

    ! initialize arrays:
    sp = 0.0 ; copy = xy

    kmp1 = km + 1; lmp1 = lm + 1; k2 = mx - km; l2 = ny - lm

    call ft_2d(copy,mx,ny,-1)

    !do k=1,kmp1; do l=1,lmp1
    do k=1,km; do l=1,lm

        kk = km + k; ll = lm + l

        ! waves: 0 <= k <= +/-kmax; 0 <= l <= +/-lmax:
        sp(k,l) = copy(k,l)/real(mx*ny)

        ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
        if (k .gt. 1) then
            sp(kk,l) = copy(k2+k,l)/real(mx*ny)
        endif

        ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
        if (l .gt. 1) then
            sp(k,ll) = copy(k,l2+l)/real(mx*ny)
        endif

        ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
        if ((k .gt. 1) .and. (l .gt. 1)) then
            sp(kk,ll) = copy(k2+k,l2+l)/real(mx*ny)
        endif

    enddo; enddo

    return
END SUBROUTINE xy_to_sp
!========================================================================

!========================================================================
SUBROUTINE sp_to_xy(sp,xy,km,lm,mx,ny)
! Map an (km,lm) spectral array onto a _bigger_ grid point array.
! Input: sp(2*km,2*lm) --- a spectral array.
! Output: xy(mx,ny) --- a grid point array.

    implicit none

    integer,                       intent(in)  :: km,lm,mx,ny
    complex, dimension(:,:),       intent(in)  :: sp
    real,    dimension(mx,ny),     intent(out) :: xy

    complex :: copy(mx,ny)
    integer :: kmp1,lmp1,k,l,k2,l2,kk,ll

    copy = 0.0

    kmp1 = km + 1; lmp1 = lm + 1; k2 = mx - km; l2 = ny - lm

    !do 10 k = 1,kmp1; do 10 l = 1,lmp1
    do k = 1,km ; do l = 1,lm

        ! waves: 0 <= k <= +/-kmax; 0 <= l <= +/-lmax:
        copy(k,l) = sp(k,l)

        kk = km + k; ll = lm + l

        ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
        if (k .ge. 1) then
           copy(k2+k,l) = sp(kk,l)
        endif

        ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
        if (l .ge. 1) then
           copy(k,l2+l) = sp(k,ll)
        endif

        ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
        if ((k .ge. 1) .and. (l .ge. 1)) then
           copy(k2+k,l2+l) = sp(kk,ll)
        endif

     enddo ; enddo

    call ft_2d(copy,mx,ny,1)
    xy = real(copy)

    return
END SUBROUTINE sp_to_xy
!========================================================================

!========================================================================
SUBROUTINE diffusion(dco)
! Compute the diffusion coefficient for del^n diffusion (n even>0).
! Tau gives the e-folding time scale for damping modes at the 
! Nyquist frequency.

    implicit none

    real, intent(out) :: dco
    real              :: temp

    ! 10/21/99 gjh mods...
    !temp = (cmplx(0.,1.)**n)
    temp = 1

    !dco = (((real(facx*kmax)**n) + (real(facy*lmax)**n))*tau)**(-1)
    dco = 1.0/(((real(facx*kmax)**2) &
              + (real(facy*lmax)**2))**(n/2)*tau)

    dco = temp*dco
    if (dco .lt. 0.0) dco = -1.0*dco

    if (verbose .gt. 1) print*,'diffusion coefficient:',dco

    return
END SUBROUTINE diffusion
!========================================================================

!========================================================================
SUBROUTINE init(infile,thB,thT)
! Initialize perturbation fields.

    implicit none

    character(len=*),               intent(in)  :: infile
    real, dimension(2*kmax,2*lmax), intent(out) :: thB,thT
    real, dimension(2*kmax,2*lmax,2)            :: theta

    real    :: fac
    integer :: twokmax,twolmax,levmax
    integer :: ncid, dimid, varid

    if (verbose .gt. 0)  print*,'initializing theta fields...'
    thB = 0.0 ; thT = 0.0

    call nc_check( nf90_open(trim(infile), NF90_NOWRITE, ncid), 'init', 'open, ' // trim(infile) )
    call nc_check( nf90_inq_dimid(ncid, 'nx', dimid), 'init', 'inq_dimid nx, ' // trim(infile) )
    call nc_check( nf90_inquire_dimension(ncid, dimid, len = twokmax), 'init ', 'inquire_dimension nx,' // trim(infile) )
    call nc_check( nf90_inq_dimid(ncid, 'ny', dimid), 'init', 'inq_dimid ny, ' // trim(infile) )
    call nc_check( nf90_inquire_dimension(ncid, dimid, len = twolmax), 'init ', 'inquire_dimension ny,' // trim(infile) )
    call nc_check( nf90_inq_dimid(ncid, 'nz', dimid), 'init', 'inq_dimid nz, ' // trim(infile) )
    call nc_check( nf90_inquire_dimension(ncid, dimid, len = levmax),  'init ', 'inquire_dimension nz,' // trim(infile) )

    if (twokmax .ne. 2*kmax .and. twolmax .ne. 2*lmax) then  
        print*,'x and y resolution from sqgInputt.nc do not match spectral_mod.f90!!!'
        print*,'x and y from spectral_mod.f90 : ', 2*kmax, 2*lmax
        print*,'x and y from sqgInput.nc      : ', twokmax,twolmax
        call nc_check( nf90_close(ncid), 'init', 'close, ' // trim(infile) )
        stop
    endif

    if ( levmax .ne. 2 ) then
        print*,'nz in sqgInput.nc is not 2!!!'
        print*,'nz must be          : ', 2
        print*,'nz from sqgInput.nc : ', levmax
        call nc_check( nf90_close(ncid), 'init', 'close, ' // trim(infile) )
        stop
    endif

    call nc_check( nf90_inq_varid(ncid, 'theta', varid), 'init', 'inq_varid theta, ' // trim(infile) )
    call nc_check( nf90_get_var(ncid, varid, theta), 'init', 'get_var theta, ' // trim(infile) )
    call nc_check( nf90_close(ncid), 'init', 'close, ' // trim(infile) )

    thB = theta(:,:,1)
    thT = theta(:,:,2)

    return
END SUBROUTINE init
!========================================================================

!========================================================================
SUBROUTINE init_jet(thbB,thbT,thbyB,thbyT,ulinB,ulinT,lam)
! Initialize basic state jet.
! thbB  = periodic basic state theta; bottom boundary (spectral grid).
! thbT  = periodic basic state theta; top boundary (spectral grid).
! thbyB = y derivative of periodic basic state theta; bottom boundary (grid point).
! thbyT = y derivative of periodic basic state theta; top boundary (grid point).
! ulinB = basic state wind; bottom boundary (grid point).
! ulinT = basic state wind; top boundary (grid point).

    implicit none

    complex, dimension(2*kmax,2*lmax), intent(out) :: thbB,thbT
    real, dimension(mmax,nmax),        intent(out) :: thbyB,thbyT
    real, dimension(mmax,nmax),        intent(out) :: ulinB,ulinT
    real,                              intent(out) :: lam

    real,    dimension(2*kmax,2*lmax) :: thbxyB,thbxyT
    complex, dimension(mmax,nmax)     :: thyB,thyT
    complex, dimension(mmax,nmax)     :: uB,uT
    complex, dimension(mmax,nmax)     :: temp
    real                              :: y,yp,dyy
    integer                           :: j

    if (verbose .gt. 0) print*,'initializing basic state...'

    ! first determine lamda, the linear shear parameter
    lam = (HW_theta(0.,0.) - HW_theta(YL,0.)) / YL
    if (verbose .gt. 1) print*,'lam = ',lam

    ! spectral variables
    dyy = YL/real(2*lmax)
    do j=1, 2*lmax
        y = real(j-1)*dyy; yp = y - (0.5*(YL - hwp))
        thbxyB(:,j) = HW_theta(y,0.) + lam*yp
        thbxyT(:,j) = HW_theta(y,H)  + lam*yp
        !print*,'periodic theta grid point:',y,yp,thbxyT(1,j)
    enddo
    ! map into spectral space at the same resolution:
    call xy_to_sp(cmplx(thbxyB,0.),thbB,2*kmax,2*lmax,kmax,lmax)
    call xy_to_sp(cmplx(thbxyT,0.),thbT,2*kmax,2*lmax,kmax,lmax)

    ! grid point variables
    dyy = YL/real(nmax)
    do j=1, nmax
        y = real(j-1)*dyy; yp = y - (0.5*(YL - hwp))
        thbyB(:,j) = HW_thetay(y,0.) + lam
        thbyT(:,j) = HW_thetay(y,H)  + lam

        ! old U
        !ulinB(:,j) = HW_ubar(y,0.)
        !ulinT(:,j) = HW_ubar(y,H)
    enddo

    ! new U: solve numerically given theta
    uB = 0.; uT = 0.
    call d_s2b( thbB,thyB,1, d_oper%dy); call d_s2b( thbT,thyT,1,d_oper%dy)
    call d_b2b(-thyB,uB,  1,-d_oper%iz); call d_b2b(-thyT,uT,  1,d_oper%iz)
    call d_b2b(-thyT,temp,1, d_oper%izo); uB = uB + temp;
    call d_b2b(-thyB,temp,1,-d_oper%izo); uT = uT + temp;
    call ft_2d(uB,mmax,nmax,1); call ft_2d(uT,mmax,nmax,1)
    ulinB = real(uB); ulinT = real(uT)
    ulinT = ulinT + lam*H

    return
END SUBROUTINE init_jet
!========================================================================

!========================================================================
FUNCTION ran1(idum)
! Numerical Recipes random number generator:

    implicit none

    REAL    ran1
    INTEGER idum

    INTEGER IA,IM,IQ,IR,NTAB,NDIV
    REAL    AM,EPS,RNMX
    PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
               NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
    INTEGER j,k,iv(NTAB),iy
    SAVE iv,iy
    DATA iv /NTAB*0/, iy /0/
    if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do j=NTAB+8,1,-1
        k=idum/IQ
        idum=IA*(idum-k*IQ)-IR*k
        if (idum.lt.0) idum=idum+IM
        if (j.le.NTAB) iv(j)=idum
        enddo
        iy=iv(1)
    end if
    k=idum/IQ
    idum=IA*(idum-k*IQ)-IR*k
    if (idum.lt.0) idum=idum+IM
    j=1+iy/NDIV
    iy=iv(j)
    iv(j)=idum
    ran1=min(AM*iy,RNMX)

    return
END FUNCTION ran1
!========================================================================

!========================================================================
SUBROUTINE dump(thspB,thspT,ilam,lam,it,outfile)
! Write theta or streamFUNCTION to disk.

    implicit none

    complex, dimension(:,:), intent(in) :: thspB,thspT
    logical,                 intent(in) :: ilam
    real,                    intent(in) :: lam
    integer,                 intent(in) :: it
    character(len=*),        intent(in) :: outfile

    integer                           :: j
    real,    dimension(2*kmax,2*lmax) :: thxyB,thxyT
    complex, dimension(2*kmax,2*lmax) :: copy

    copy = thspB
    call sp_to_xy(copy,thxyB,kmax,lmax,2*kmax,2*lmax)
    copy = thspT
    call sp_to_xy(copy,thxyT,kmax,lmax,2*kmax,2*lmax)

    ! add in linear shear
    if (ilam) then 
        do j=1, 2*lmax
            ! fixed 08/10/2005 GJH & RBM
            thxyB(:,j) = thxyB(:,j) - lam*real(j-1)*YL/real(2*lmax)
            thxyT(:,j) = thxyT(:,j) - lam*real(j-1)*YL/real(2*lmax)
        enddo
    endif

    if (verbose .gt. 0) print*,'Writing to disk...'
    call write_diag(outfile,it,thxyB,thxyT)
  
    return
END SUBROUTINE dump
!========================================================================

!========================================================================
SUBROUTINE advect(u,v,f_x,f_y,fb_y,h_x,h_y,ub,tf,lam,lap)
! Spectral advection on the (x,y) grid.

    implicit none

    real, dimension(:,:),       intent(in)  :: u,v,f_x,f_y,fb_y,ub,lap,h_x,h_y
    real,                       intent(in)  :: lam
    real, dimension(:,:),       intent(out) :: tf

    real    :: x,y,dx,dy
    real    :: terr
    real    :: famp,rw
    real    :: rphasex,rphasey
    integer :: i,j
    integer, save :: iseed

    dx = XL/real((mmax)); dy = YL/real((nmax))

    ! random wave-one forcing; random walk in phase (30 August 2006)
    rw   = 4.0   ! random wavenumber (32.0, 4.0)
    famp = 0.0   ! forcing amplitude (0.0, 0.1*H)
    rphasex = 2.0*pi*(ran1(iseed)-0.5)
    rphasey = 2.0*pi*(ran1(iseed)-0.5)

    do i=1,mmax; do j=1,nmax

        x = real(i-1)*dx; y = real(j-1)*dy

        if (linear) then 
            tf(i,j) = -((v(i,j) * (fb_y(i,j) - lam)) + &
                        (ub(i,j) * f_x(i,j)))
        else
            tf(i,j) = -((v(i,j) * (f_y(i,j) + fb_y(i,j) - lam)) + &
                       ((u(i,j) + ub(i,j)) * f_x(i,j)))
        endif

        tf(i,j) = tf(i,j) - (ekman(x,y)*lap(i,j)) ! Ekman layer

        ! terrain:
        if (iterr) then
            terr = -( (v(i,j) * h_y(i,j)) + ((u(i,j) + ub(i,j)) * h_x(i,j)) )
            tf(i,j) = tf(i,j) + terr 
        endif

        ! random wave-one forcing; random walk in phase (30 August 2006)
        tf(i,j) = tf(i,j) - famp*sin((rw*x*2*pi/XL)-rphasex)*sin((rw*y*2*pi/YL)-rphasey)

    enddo; enddo

    return
END SUBROUTINE advect
!========================================================================

!========================================================================
SUBROUTINE tadv(dat,tend,told,told2,dco,first)
! Time-advance SUBROUTINE.

    implicit none

    complex, dimension(2*kmax,2*lmax), intent(inout) :: dat, told, told2
    complex, dimension(mmax,nmax),     intent(in)    :: tend
    real,                              intent(in)    :: dco
    logical,                           intent(in)    :: first
    integer :: k,l,kk,ll
    real    :: ak,bl,dts
    complex :: ttsp

    if (first) then
        dts = dt*(12.0/23.0)
    else
        dts = dt
    endif

    do k=1,kmax; do l=1,lmax

        ! 0 <= k <= +/-kmax; 0 <= l <= +/-lmax:
        ak = facx*real(k - 1); bl = facy*real(l - 1)
        ttsp = tend(k,l)/real(mmax*nmax)
        call tstep_ab(ttsp,dat(k,l),told(k,l),told2(k,l), &
                      ak,bl,dco,dts)

        kk = kmax + k; ll = lmax + l

        ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
        !         if (k .gt. 1 .and. k .lt. kmax) then
        if (k .gt. 1) then
            ak = -1.*facx*real(kmaxp1 - k); bl = facy*real(l - 1)
            ttsp = tend(k2+k,l)/real(mmax*nmax)
            call tstep_ab(ttsp,dat(kk,l),told(kk,l), &
                          told2(kk,l),ak,bl,dco,dts)
        endif

        ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
        !         if (l .le. lmax) then
        if (l .gt. 1) then
            ak = facx*real(k-1); bl = -1.*facy*real(lmaxp1 - l)
            ttsp = tend(k,l2+l)/real(mmax*nmax)
            call tstep_ab(ttsp,dat(k,ll),told(k,ll), &
                          told2(k,ll),ak,bl,dco,dts)
        endif

        ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
        !         if (k .le. kmax .and. l .le. lmax) then
        if ((k .gt. 1) .and. (l .gt. 1)) then
            ak=-1.*facx*real(kmaxp1 - k); bl=-1.*facy*real(lmaxp1 - l)
            ttsp = tend(k2+k,l2+l)/real(mmax*nmax)
            call tstep_ab(ttsp,dat(kk,ll),told(kk,ll), &
                          told2(kk,ll),ak,bl,dco,dts)
        endif

    enddo; enddo

    return
END SUBROUTINE tadv
!========================================================================

!========================================================================
SUBROUTINE tstep_ab(tend,dat,told,told2,ak,bl,dco,dts)
! 3rd-order Adams-Bashforth time step (Durran 1991).

    implicit none

    complex, intent(inout) :: dat, told, told2
    complex, intent(in)    :: tend
    real, intent(in)       :: ak,bl,dco,dts

    complex :: new, temp
    real    :: tfac,relax

    ! changed 02/07/03 to include relaxation to jet
    relax = 0.0

    if (trl .lt. 1.e3 .and. trl .gt. 0) relax = 1.0/trl 
    ! 12/17/2008: damps l=0 faster for stability (used to zero in main block)
    if (bl .eq. 0) relax = relax*4

    tfac = dts*((dco*(((ak**2)+(bl**2))**(n/2))) + relax)

    ! new t-step:
    !tfac = dco*dts*((ak**n)+(bl**n))
    !tfac = dco*dts*((ak**2)+(bl**2))**(n/2)

    told = told*exp(-1.0*tfac)
    told2 = told2*exp(-tfac)

    temp = (23.0*tend) - (16.0*told) + (5.0*told2)
    new = dat + ( (dts/12.)*(temp) )

    !dat=new*exp(-1.0*tfac); told2=told*exp(tfac); told=tend
    dat=new*exp(-1.0*tfac); told2=told; told=tend

    return
END SUBROUTINE tstep_ab
!========================================================================

!========================================================================
SUBROUTINE ft_2d(f,ni,nj,isign)
! FFT-calling SUBROUTINE.

    implicit none

    integer,                   intent(in)    :: ni, nj, isign
    complex, dimension(ni,nj), intent(inout) :: f

    real, dimension(ni,nj) :: re, im

    re = real(f); im = aimag(f)
    call fft(re,im,ni*nj,ni,ni,   isign)
    call fft(re,im,ni*nj,nj,nj*ni,isign)
    f = cmplx(re,im)

    return
END SUBROUTINE ft_2d
!========================================================================

!========================================================================
SUBROUTINE d_setup()
! Set up matrices for derivatives and integrals.

    implicit none

    complex, dimension(2*kmax,2*lmax) :: dx,dy,dz,dzo,iz,izo,Id
    real    :: m
    integer :: k,l

    ! dx
    dx(1,:) = 0.; dx(kmax+1,:) = 0.
    do k=2,kmax
        dx(k,:) = facx*real(k-1)*cmplx(0.,1.)
        dx(2*kmax-k+2,:) = -1.*dx(k,:)
    enddo

    ! dy
    dy(:,1) = 0.; dy(:,lmax+1) = 0.
    do l=2,lmax
        dy(:,l) = facy*real(l-1)*cmplx(0.,1.)
        dy(:,2*lmax-l+2) = -1.*dy(:,l)
    enddo

    ! dz,dzo
    do k=1,2*kmax; do l=1,2*lmax
        if (k .eq. 1 .and. l .eq. 1) then 
            dz(k,l) = 0.; dzo(k,l) = 0.; iz(k,l) = 0.; izo(k,l) = 0.
        elseif (k .eq. 1 .and. l .eq. lmax+1) then
            dz(k,l) = 0.; dzo(k,l) = 0.; iz(k,l) = 0.; izo(k,l) = 0.
        elseif (k .eq. kmax+1 .and. l .eq. lmax+1) then
            dz(k,l) = 0.; dzo(k,l) = 0.; iz(k,l) = 0.; izo(k,l) = 0.
        elseif (l .eq. 1 .and. k .eq. kmax+1) then
            dz(k,l) = 0.; dzo(k,l) = 0.; iz(k,l) = 0.; izo(k,l) = 0.
        else
            ! must force real, since small complex residual may remain:
            m = real(((-(dx(k,1)**2)) - (dy(1,l)**2))**0.5)
            dz(k,l) = m / tanh(m*H) ! h->\infty: tanh->1, dz->m
            iz(k,l) = 1./ (m * tanh(m*H)) ! h->\infty: iz-> 1/m
            ! operators for opposing boundary
            if (m*H .lt. 50 .and. model .ne. 0) then
                dzo(k,l) = m / sinh(m*H) ! h->\infty: dzo->0
                izo(k,l) = 1./ (m * sinh(m*H)) ! h->\infty: izo->0
            elseif (model .eq. 0) then 
                print *,'SELECTED 2D EULER DYNAMICS....'
                dz(k,l) = m**2
                iz(k,l) = 1./dz(k,l)
            else
                dzo(k,l) = 0.
                izo(k,l) = 0.
            endif
        endif
    enddo; enddo

    ! Id
    Id = 1.0

    d_oper%dx  = dx ; d_oper%dy  = dy
    d_oper%dz  = dz ; d_oper%dzo = dzo
    d_oper%iz  = iz ; d_oper%izo = izo
    d_oper%Id  = Id

    return
END SUBROUTINE d_setup
!========================================================================

!========================================================================
SUBROUTINE d_s2b(temp_in,temp_out,dflag,dn)
! small array to big array.
! dflag =  n: n derivatives. dflag = -n: n integrations.

    implicit none

    complex, dimension(:,:), intent(in)  :: temp_in,dn
    integer,                 intent(in)  :: dflag
    complex, dimension(:,:), intent(out) :: temp_out

    integer :: k,l,kk,ll

    do k = 1,kmax; do l = 1,lmax

        ! waves: 0 <= k <= +/-kmax; 0 <= l <= +/-lmax:
        temp_out(k,l) = (dn(k,l)**dflag)*temp_in(k,l)

        kk = kmax + k; ll = lmax + l

        ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
        if (k .gt. 1) temp_out(k2+k,l) = (dn(kk,l)**dflag)*temp_in(kk,l)

        ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
        if (l .gt. 1) temp_out(k,l2+l) = (dn(k,ll)**dflag)*temp_in(k,ll)

        ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
        if ((k .gt. 1) .and. (l .gt. 1)) temp_out(k2+k,l2+l) = (dn(kk,ll)**dflag)*temp_in(kk,ll)

    enddo; enddo

    return
END SUBROUTINE d_s2b
!========================================================================

!========================================================================
SUBROUTINE d_s2s(temp_in,temp_out,dflag,dn)
! small array to small array.
! dflag =  n: n derivatives. dflag = -n: n integrations.

    implicit none

    complex, dimension(:,:), intent(in)  :: temp_in,dn
    integer,                 intent(in)  :: dflag
    complex, dimension(:,:), intent(out) :: temp_out

    integer :: k,l

    temp_out = 0.0

    do k = 1,2*kmax; do l = 1,2*lmax
        if (dn(k,l) .ne. 0) temp_out(k,l) = (dn(k,l)**dflag)*temp_in(k,l)
    enddo; enddo

    return
END SUBROUTINE d_s2s
!========================================================================

!========================================================================
SUBROUTINE d_b2b(temp_in,temp_out,dflag,dn)
! big array to big array.
! dflag =  n: n derivatives. dflag = -n: n integrations.

    implicit none

    complex, dimension(:,:), intent(in)  :: temp_in
    complex, dimension(:,:), intent(in)  :: dn
    integer,                 intent(in)  :: dflag
    complex, dimension(:,:), intent(out) :: temp_out

    integer :: k,l,kk,ll

    temp_out = 0.0

    do k = 1,kmax; do l = 1,lmax

        ! waves: 0 <= k <= +/-kmax; 0 <= l <= +/-lmax:
        if (dn(k,l) .ne. 0) temp_out(k,l) = (dn(k,l)**dflag)*temp_in(k,l)

        kk = kmax + k; ll = lmax + l

        ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
        if (k .gt. 1) then
            if (dn(kk,l) .ne. 0) temp_out(k2+k,l) = (dn(kk,l)**dflag)*temp_in(k2+k,l)
        endif

        ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
        if (l .gt. 1) then
            if (dn(k,ll) .ne. 0) temp_out(k,l2+l) = (dn(k,ll)**dflag)*temp_in(k,l2+l)
        endif

        ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
        if ((k .gt. 1) .and. (l .gt. 1)) then
            if (dn(kk,ll) .ne. 0) temp_out(k2+k,l2+l) = (dn(kk,ll)**dflag)*temp_in(k2+k,l2+l)
        endif

    enddo; enddo

    return
END SUBROUTINE d_b2b
!========================================================================

!========================================================================
SUBROUTINE d_b2s(temp_in,temp_out,dflag,dn)
! big array to small array.
! dflag =  n: n derivatives. dflag = -n: n integrations.

    implicit none

    complex, dimension(:,:), intent(in)  :: temp_in
    complex, dimension(:,:), intent(in)  :: dn
    integer,                 intent(in)  :: dflag
    complex, dimension(:,:), intent(out) :: temp_out

    integer :: k,l,kk,ll

    do k = 1,kmax; do l = 1,lmax

        ! waves: 0 <= k <= +/-kmax; 0 <= l <= +/-lmax:
        temp_out(k,l) = (dn(k,l)**dflag)*temp_in(k,l)

        kk = kmax + k; ll = lmax + l

        ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
        if (k .ge. 1) temp_out(kk,l) = (dn(kk,l)**dflag)*temp_in(k2+k,l)

        ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
        if (l .ge. 1) temp_out(k,ll) = (dn(k,ll)**dflag)*temp_in(k,l2+l)

        ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
        if ((k .ge. 1) .and. (l .ge. 1)) temp_out(kk,ll) = (dn(kk,ll)**dflag)*temp_in(k2+k,l2+l)

    enddo; enddo

    return
END SUBROUTINE d_b2s
!========================================================================

!========================================================================
FUNCTION HW_ubar(y,z)
! Hoskins-West base state u wind

    implicit none

    real, intent(in) :: y,z
    real :: HW_ubar
    real :: yp

    yp = y - (0.5*(YL - hwp))
    if ( (yp .lt. 0.0) .and. (amu .ne. 0.0) ) yp = 0.0
    if ( (yp .gt. hwp) .and. (amu .ne. 0.0) ) yp = hwp
    HW_UBAR = shear*z - (amu/2.0)*(z + ((sinh(ryl*z)/sinh(ryl))*cos(yp*ryl)))

    return
END FUNCTION HW_ubar

!========================================================================
FUNCTION HW_theta(y,z)
! Hoskins-West base state potential temperature
! 9/2006: added amp and correction factor for HW 'hiccup'
! NOW INCOMPATIBLE WITH HW_ubar!!! MUST solve for U numerically!!!

    implicit none

    real, intent(in)  :: y,z

    real :: HW_theta
    real :: yp

    real, parameter :: amp    = 1.0 !control jet strength (GJH's value 1.50)
    real, parameter :: hiccup = 1.0 !                     (GJH's value 0.75)

    yp = y - (0.5*(YL - hwp))
    if ( (yp .lt. 0.0) .and. (amu .ne. 0.0) ) yp = 0.0
    if ( (yp .gt. hwp) .and. (amu .ne. 0.0) ) yp = hwp
    HW_theta = (-shear*yp) + (amu/2.) * (yp+(hiccup*(cosh(ryl*z)/sinh(ryl))*sin(yp*ryl)))
    HW_theta = amp*HW_theta

    return
END FUNCTION HW_theta
!========================================================================

!========================================================================
FUNCTION HW_thetay(y,z)
! Hoskins-West base state potential temperature gradient

    implicit none

    real, intent(in) :: y,z

    real :: HW_thetay
    real :: yp

    yp = y - (0.5*(YL - hwp))
    ! fixed this block 02/20/03

    if ( amu .eq. 0.0 ) then
        HW_thetay = -shear
     else
        if ( (yp .lt. 0.0) .or. (yp .gt. hwp) ) then 
            HW_thetay = 0.0
        else
            HW_thetay = -shear + (amu/2.0) + ((amu*ryl/2.0)*(cosh(ryl*z)*cos(yp*ryl))/sinh(ryl))
        endif
    endif

    return
END FUNCTION HW_thetay
!========================================================================

!========================================================================
SUBROUTINE norm(thB,thT,itime)
! hacked from qg code to check mode growth rates

    implicit none

    complex, dimension(2*kmax,2*lmax), intent(in) :: thB,thT
    integer,                           intent(in) :: itime

    real :: V,Vgr

    real, save :: V_old,V_0

    ! enstrophy norm (for normal modes, any norm will do)
    V = (sum(abs(thB)**2) + sum(abs(thT)**2))
    V = 0.5*(V**0.5)

    ! wait 5 time steps for norms to stabilize
    if (itime .eq. 5) V_0 = V

    if (itime .gt. 2) then
        Vgr = (log(V) - log(V_old))/dt
        !Vgr = ((V - V_old)/(dt*0.5*(V+V_old)))
    endif

    if ( verbose .gt. 0 ) then
        print*,'V inst growth rates:',itime*dt,Vgr
        !if (itime .gt. 5) print*,'V mean growth rates:',(alog(V/V_0))/(real(itime-5)*dt)
    endif

    V_old = V

    return
END SUBROUTINE norm
!========================================================================

!========================================================================
FUNCTION ekman(x,y)
! Ekman pumping at the bottom surface

    implicit none

    real, intent(in) :: x,y
    real :: ekman

    if (x .ge. 2.*XL/3.) then
        ekman = gamma
    else
        ekman = 0.
        ekman = gamma
        !ekman = gamma/4.
    endif

    ekman = gamma
    !if (abs(y-(YL/2)) .ge. 3) ekman = 0.5

    return
END FUNCTION ekman
!========================================================================

!========================================================================
SUBROUTINE write_diag(output_file,it,thB,thT)
! write to disk

    implicit none

    character(len=*),               intent(in) :: output_file
    integer,                        intent(in) :: it
    real, dimension(2*kmax,2*lmax), intent(in) :: thB, thT
    real, dimension(2*kmax,2*lmax,2)           :: theta

    integer :: ncid, vardim(4), varid
    integer :: k, count(4), start(4), ierr
    real    :: time

    if ( it .eq. 0 ) then

        ! Create a new NetCDF file
        call nc_check( nf90_create(output_file, ior(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid), 'write_diag', 'create ' // trim(output_file) )

        ! Define dimensions
        call nc_check( nf90_def_dim(ncid, "nx",   2*kmax,         vardim(1)), 'write_diag', 'def_dim, nx '   // trim(output_file))
        call nc_check( nf90_def_dim(ncid, "ny",   2*lmax,         vardim(2)), 'write_diag', 'def_dim, ny '   // trim(output_file) )
        call nc_check( nf90_def_dim(ncid, "nz",   2,              vardim(3)), 'write_diag', 'def_dim, nz '   // trim(output_file) )
        call nc_check( nf90_def_dim(ncid, "time", NF90_UNLIMITED, vardim(4)), 'write_diag', 'def_dim, time ' // trim(output_file) )

        ! Define variables
        call nc_check( nf90_def_var(ncid, "time",   NF90_FLOAT, vardim(4), varid), 'write_diag', 'def_var, time '  // trim(output_file) )
        call nc_check( nf90_def_var(ncid, "theta",  NF90_FLOAT, vardim,    varid), 'write_diag', 'def_var, theta ' // trim(output_file) )
        call nc_check( nf90_put_att(ncid, varid, "description", "potential temperature"), 'write_diag', 'put_att, description ' // trim(output_file) )

        ! Put global attributes
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "model", model), 'write_diag', 'put_att, model ' // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "ntims", ntims), 'write_diag', 'put_att, ntims ' // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "dt",    dt),    'write_diag', 'put_att, dt '    // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "iplot", iplot), 'write_diag', 'put_att, iplot ' // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "XL",    XL),    'write_diag', 'put_att, XL '    // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "YL",    YL),    'write_diag', 'put_att, YL '    // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "H",     H),     'write_diag', 'put_att, H '     // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "Ross",  Ross),  'write_diag', 'put_att, Ross '  // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "gamma", gamma), 'write_diag', 'put_att, gamma ' // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "n",     n),     'write_diag', 'put_att, n '     // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "tau",   tau),   'write_diag', 'put_att, tau '   // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "trl",   trl),   'write_diag', 'put_att, trl '   // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "amu",   amu),   'write_diag', 'put_att, amu '   // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "shear", shear), 'write_diag', 'put_att, shear ' // trim(output_file) )

        call nc_check( nf90_enddef(ncid), 'write_diag', 'enddef, ' // trim(output_file) )
        call nc_check( nf90_close(ncid),  'write_diag', 'close, '  // trim(output_file) )

    else 

        time = (it-1)*dt
        k = 1 + (it-1)/iplot

        count(1) = 2*kmax;  start(1) = 1
        count(2) = 2*lmax;  start(2) = 1
        count(3) = 2;       start(3) = 1
        count(4) = 1;       start(4) = k

        theta(:,:,1) = thB
        theta(:,:,2) = thT

        ! Open the netCDF file, write variables and close
        call nc_check( nf90_open(output_file, NF90_WRITE, ncid),       'write_diag', 'open, '            // trim(output_file) )
        call nc_check( nf90_inq_varid(ncid, "time",   varid) ,         'write_diag', 'inq_varid, time '  // trim(output_file) )
        call nc_check( nf90_put_var(ncid, varid, time, (/k/)),         'write_diag', 'put_var, time '    // trim(output_file) )
        call nc_check( nf90_inq_varid(ncid, "theta", varid),           'write_diag', 'inq_varid, theta ' // trim(output_file) )
        call nc_check( nf90_put_var(ncid, varid, theta, start, count), 'write_diag', 'put_var, theta '   // trim(output_file) )
        call nc_check( nf90_close(ncid),                               'write_diag', 'close, '           // trim(output_file) )

    endif

    return
END SUBROUTINE write_diag
!========================================================================

!========================================================================
SUBROUTINE terrain(hx,hy,hu,hv)
! Terrain and barotropic wind on the _advection_ grid:
! feb 2005: added tropo winds from z = 0 for Hsqg

    implicit none

    real, dimension(mmax,nmax), intent(out) :: hx,hy,hu,hv

    complex, dimension(mmax,nmax)     :: hxsp,hysp
    complex, dimension(mmax,nmax)     :: hh
    complex, dimension(mmax,nmax)     :: hxs,hys
    real,    dimension(2*kmax,2*lmax) :: hspR
    complex, dimension(2*kmax,2*lmax) :: hspC

    real :: ubaro(nmax),x,y,ddx,ddy,xcen,ycen,amdx,amdy
    real :: rr,c1,c2,ak,bl,lam
    integer :: i,j,k,l,LL,MM,kk,mmd2,nmd2,mmaxp1,nmaxp1
    real,    dimension(mmax,nmax)     :: Rblank
    complex, dimension(2*kmax,2*lmax) :: Cblank

    ! barotropic wind
    hu = 0.; hv = 0.

    ddx = XL/real((mmax)); ddy = YL/real((nmax))
    xcen = mmax*ddx/2; ycen = nmax*ddy/2

    if ( verbose .gt. 1) print*,'terrain center : ',xcen,ycen

    ! gaussian topography:
    do i=1,mmax; do j=1,nmax
        x = (i-1)*ddx; y = (j-1)*ddy
        amdx=min(abs(x-xcen),abs(XL+x-xcen),abs(XL-x+xcen))
        amdy=min(abs(y-ycen),abs(YL+y-ycen),abs(YL-y+ycen))
        rr = (((amdx/asx)**2) + ((amdy/asy)**2))**0.5
        hh(i,j) = hamp*exp(-1.0*((rr/asig)**2))
    enddo; enddo

    if ( verbose .gt. 1 ) print*,'terrain height : ',real(hh(:,nmax/2))

    call ft_2d(hh,mmax,nmax,-1)
    hh = hh / (mmax*nmax)

    ! form spectral h_x, h_y:
    call d_b2b(hh,hxs,1,d_oper%dx) ! x derivative: small to big domain.
    call d_b2b(hh,hys,1,d_oper%dy) ! y derivative: small to big domain.

    ! back to grid point space
    call ft_2d(hxs,mmax,nmax,1)
    call ft_2d(hys,mmax,nmax,1)
    hx = real(hxs)
    hy = real(hys)

    print*,hx(:,nmax/2)
    if ( verbose .gt. 1 ) print*,'terrain height : ',hx(:,nmax/2)

    if (model .eq. 4) then ! HsQG needs topo winds on the tropo

        Rblank = 0.; Cblank = 0. ! for HsQG inversion call

        !first make spectral topo height on small grid
        ddx = XL/real((2*kmax)); ddy = YL/real((2*lmax))
        do i=1,2*kmax; do j=1,2*lmax
            x = (i-1)*ddx; y = (j-1)*ddy
            amdx=min(abs(x-xcen),abs(XL+x-xcen),abs(XL-x+xcen))
            amdy=min(abs(y-ycen),abs(YL+y-ycen),abs(YL-y+ycen))
            rr = (((amdx/asx)**2) + ((amdy/asy)**2))**0.5
            hspR(i,j) = hamp*exp(-1.0*((rr/asig)**2))
        enddo; enddo

        if ( verbose .gt. 1 ) print*,'max topo = ',maxval(abs(hspR))
        call xy_to_sp(cmplx(hspR,0.),hspC,2*kmax,2*lmax,kmax,lmax)
        if ( verbose .gt. 1 ) print*,'max topo spectral = ',maxval(abs(hspC))
        call invert(-hspC,Cblank,Rblank,Rblank,Rblank,Rblank,Rblank, &
                    hv,Rblank,hu,Cblank,Cblank,Rblank,Rblank,Rblank, &
                    Rblank,.TRUE.,.TRUE.,.TRUE.,lam,Cblank,Cblank,Rblank)
        ! hu and hv have the tropo winds due to topography
        if ( verbose .gt. 1 ) print*,'max tropo winds due to topography: ',maxval(abs(hu)),maxval(abs(hv))

    endif

    return
END SUBROUTINE terrain
!========================================================================

!========================================================================
SUBROUTINE scaling(grav,tnot,Ns,km,Cor,Ls,Hs,Us,Ws,Ps,Ts)
! scaling parameters:

    implicit none

    real, intent(out) :: grav,tnot,Ns,km,Cor,Ls,Hs,Us,Ws,Ps,Ts

    grav = 9.81        ! gravitational constant (m s^-2)
    tnot = 300.        ! theta constant at Z = 0 (K)
    Ns   = 1.e-2       ! bouyancy frequency
    km   = 1000.       ! 1000 m in a kilometer
    Cor  = 1.e-4       ! Coriolis parameter (s^-1)
    Hs   = 10.*km      ! vertical length scale (m)
    Ls   = Ns*Hs/Cor   ! horizontal length scale (m)

    if (Ross .ne. 0) then 
        Us = Ross*Cor*Ls   ! Horizontal wind scale (m s^-1)
        if (verbose .gt. 1)  print*,'Us = ',Us
    else
        Us = 1.0
    endif

    Ws = Ross*Hs*Us/Ls     ! Vertical wind scale (m s^-1)
    Ps = Us*Cor*Ls         ! geopotential scale (m^2 s^-2)
    Ts = Ps*tnot/(grav*Hs) ! potential temperature (K)

    return
END SUBROUTINE scaling
!========================================================================

!========================================================================
SUBROUTINE dx_echo(dxs,dys,dzs,dts)
! send back dimensional values for dx, dy, dz, and dt.

    implicit none

    real, intent(out) :: dxs,dys,dzs,dts
    real :: grav,tnot,Ns,km,Cor,Ls,Hs,Us,Ws,Ps,Ts

    call scaling(grav,tnot,Ns,km,Cor,Ls,Hs,Us,Ws,Ps,Ts)

    dxs = Ls*XL/real(2*kmax) ! meters
    dys = Ls*YL/real(2*kmax) ! meters
    dzs = Hs*ZH/real(pmax)   ! meters
    dts = Ls*dt/Us           ! seconds

    return
END SUBROUTINE dx_echo
!========================================================================

!========================================================================
SUBROUTINE nc_check(ierr,subr_name,context)

  ! check for netcdf errors

  implicit none

  integer,          intent(in) :: ierr
  character(len=*), intent(in) :: subr_name, context

  character(len=129) :: error_msg

  if (ierr /= nf90_noerr) then
    error_msg = trim(subr_name) // ': ' // trim(context) // ': ' // trim(nf90_strerror(ierr))
    print*,trim(adjustl(error_msg))
    stop
  end if

  return
END SUBROUTINE nc_check
!========================================================================

END MODULE sqg_mod
!========================================================================

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
