! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module pe2lyr_mod

! 2-layer isentropic primitive equation model on a sphere.
!
! See 
! Zou, X., Barcilon, A., Navon, I.M., Whitaker, J., Cacuci, D.G.. 1993:
! An Adjoint Sensitivity Study of Blocking in a Two-Layer Isentropic Model.
! Monthly Weather Review: Vol. 121, No. 10, pp. 2833-2857.
! for details. 
!
! Contact: <Jeffrey.S.Whitaker@noaa.gov>

use spharmt_mod

use utilities_mod, only : register_module, error_handler, E_ERR

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"


type modelvars

   real cp,g,r,p0,efold,dt,delpi,rot,pitop,tdrag,tdiab,theta1,delth,hmax
   integer ndiss
   complex, dimension(:,:), pointer :: vrtnm,dvrtdtnm,divnm,ddivdtnm, &
   delpinm,ddelpidtnm
   real, dimension(:), pointer :: dissnm,pibg
   real, dimension(:,:), pointer :: hg

   logical :: isinitialized

end type modelvars

logical, save :: module_initialized = .false.



 contains


subroutine initialize_module

   call register_module(source, revision, revdate)
   module_initialized = .true.

end subroutine initialize_module



 subroutine modelvars_init(model_dat,sphere_dat,&
 cp,g,r,p0,efold,ndiss,dt,delta_pi,ztop,rot,tdrag,tdiab,theta1,delth,hmax)

! initialize a modelvars object.

   type (modelvars) :: model_dat
   type (sphere) :: sphere_dat
   real del((sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2)
   real datag(sphere_dat%nlons,sphere_dat%nlats)

   if ( .not. module_initialized ) call initialize_module

   if (.not. sphere_dat%isinitialized) call error_handler(E_ERR, &
      'modelvars_init','sphere_dat object must be initialized before call', &
      source, revision, revdate)

   ntrunc = sphere_dat%ntrunc
   nmdim = (ntrunc+1)*(ntrunc+2)/2
   nlats = sphere_dat%nlats
   nlons = sphere_dat%nlons

   model_dat%dt = dt
   model_dat%cp = cp
   model_dat%g = g
   model_dat%r = r 
   model_dat%p0 = p0
   model_dat%efold = efold
   model_dat%ndiss = ndiss
   model_dat%delth = delth
   model_dat%theta1 = theta1
   model_dat%rot = rot
!  ztop = (r*theta1/g)*(1. + ptop/p0)
   pitop = cp - g*ztop/theta1
   model_dat%pitop = pitop
   model_dat%tdrag = tdrag
   model_dat%tdiab = tdiab
   delpi = delta_pi*(cp-pitop)
   model_dat%delpi = delpi
   model_dat%hmax = hmax

   allocate(model_dat%vrtnm(nmdim,2))
   allocate(model_dat%dvrtdtnm(nmdim,2))
   allocate(model_dat%divnm(nmdim,2))
   allocate(model_dat%ddivdtnm(nmdim,2))
   allocate(model_dat%delpinm(nmdim,2))
   allocate(model_dat%ddelpidtnm(nmdim,2))
   allocate(model_dat%dissnm(nmdim))
   allocate(model_dat%pibg(nlats))
   allocate(model_dat%hg(nlons,nlats))

 
!==> define damping coefficients.

   mwaves = ntrunc+1
   delmax = (mwaves-1)*mwaves
   del = real(sphere_dat%indxn)*real(sphere_dat%indxn+1.)
   model_dat%dissnm = -(1./efold)*(del/delmax)**(ndiss/2)
 
!==> define radiative state interface exnr function.
 
   pi = 4.*atan(1.0)
   pibgav = 0.
   do j=1,nlats
      rlat = asin(sphere_dat%gaulats(j))
      f1 = cos(2.*rlat)
      f3 = (1./2.)*f1*( (sin(2.*rlat))**2 + 2.)
      df3 = -sin(2.*rlat)* &
            ((sin(2.*rlat))**2+2.) + (1./2.)* &
            cos(2.*rlat)*4.*sin(2.*rlat)*cos(2.*rlat)
      model_dat%pibg(j) = 0.5*(cp+pitop) + 0.5*f3*delpi
      pibgav = pibgav + 0.5*sphere_dat%weights(j)*model_dat%pibg(j)
      s = sin(rlat)
      c = cos(rlat)
      ct = c/s
      dpidphi = 0.5*delpi*df3
      dmdphi = delth*dpidphi
      rad = (sphere_dat%rsphere*model_dat%rot*c)**2 - ct*dmdphi
      ubg = c*(-sphere_dat%rsphere*model_dat%rot*c + sqrt( rad ) )
      zint = (theta1/g)*(cp - model_dat%pibg(j))
!     print *,(180./pi)*rlat,zint/1000.,ubg/c
   enddo

!==> define orography (zonal wave #2).

   do j=1,nlats
      rmu = sphere_dat%gaulats(j)
      rlatdep=4.*(rmu**2-rmu**4)
      do i=1,nlons
         rlon = 2.*pi*float(i-1)/float(nlons)
         model_dat%hg(i,j) = hmax*rlatdep*sin(2.*rlon)
      enddo
   enddo
 
!==> default initial state is state of rest plus barotropic vort pert.

   model_dat%delpinm = (0.,0.)
   model_dat%divnm = (0.,0.)
   model_dat%delpinm(1,1) = (2./sqrt(2.))*cmplx(cp-pibgav,0.)
   model_dat%delpinm(1,2) = (2./sqrt(2.))*cmplx(pibgav-pitop,0.)
!==> scale of perturbation is 20 +/- 15 degrees.
! use random number so runs are not identical.
   call random_seed
   call random_number(delon)
   delon = 20. + 15.*(2.*delon-1.)
   delon = delon*(pi/180.)
   amppsi = -1.e7
   do j=1,nlats
      rlat = asin(sphere_dat%gaulats(j))
      rlatdep = (sin(2.*rlat))**12
      do i=1,nlons
         rlon = float(i-1)*2.*pi/float(nlons)
         rlondep = exp(-((rlon-pi)/delon)**2) 
         datag(i,j) = amppsi*rlatdep*rlondep
      enddo
   enddo
   call spharm(sphere_dat,datag,model_dat%vrtnm(:,1),1)
   model_dat%vrtnm(:,1) = sphere_dat%lap(:)*model_dat%vrtnm(:,1)
   model_dat%vrtnm(:,2) = model_dat%vrtnm(:,1)

   model_dat%isinitialized = .TRUE.

 end subroutine modelvars_init



 subroutine modelvars_destroy(model_dat)

! deallocate pointers in modelvars object.

   type (modelvars) :: model_dat

   if ( .not. module_initialized ) call initialize_module

   if (.not. model_dat%isinitialized) return

   deallocate(model_dat%vrtnm)
   deallocate(model_dat%dvrtdtnm)
   deallocate(model_dat%divnm)
   deallocate(model_dat%ddivdtnm)
   deallocate(model_dat%delpinm)
   deallocate(model_dat%ddelpidtnm)
   deallocate(model_dat%dissnm)
   deallocate(model_dat%pibg)
   deallocate(model_dat%hg)

 end subroutine modelvars_destroy



subroutine gettend(model_dat,sphere_dat)

 type (sphere) :: sphere_dat
 type (modelvars) :: model_dat

 COMPLEX, dimension(sphere_dat%ntrunc+1,sphere_dat%nlats) :: coeff1,coeff2
 COMPLEX, dimension((sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2) :: scrnm
 REAL, dimension(sphere_dat%nlons,sphere_dat%nlats) :: &
 vec1,vec2,sinlats,coslats_sq,thtadot,totdelpig
 REAL, dimension(sphere_dat%nlons,sphere_dat%nlats,2) :: &
 ug,vg,mg,delpig,vrtg

 if ( .not. module_initialized ) call initialize_module

 if (.not. model_dat%isinitialized) call error_handler(E_ERR, &
    'gettend','model_dat object must be initialized before call', &
    source, revision, revdate) 

!==> compute u,v,vorticity,div,delpi on grid.
 
 do k=1,2
    call getuv(sphere_dat,model_dat%vrtnm(:,k),model_dat%divnm(:,k),&
               ug(:,:,k),vg(:,:,k))
    call spharm(sphere_dat,vrtg(:,:,k),model_dat%vrtnm(:,k),-1)
    call spharm(sphere_dat,delpig(:,:,k),model_dat%delpinm(:,k),-1)
 enddo

!==> compute montgomery potential

 mg(:,:,1) = model_dat%theta1*(delpig(:,:,1)+delpig(:,:,2)+model_dat%pitop) + &
             model_dat%g*model_dat%hg
 mg(:,:,2) = mg(:,:,1) + model_dat%delth*(delpig(:,:,2)+model_dat%pitop)

!==> compute diabatic heating and planetary vort.

 do j=1,sphere_dat%nlats
    totdelpig(:,j) = delpig(:,j,1)+delpig(:,j,2)
    thtadot(:,j) = &
    2.*(model_dat%pibg(j)-(delpig(:,j,2)+model_dat%pitop))*model_dat%delth/ &
    (model_dat%tdiab*totdelpig(:,j))
    sinlats(:,j) = sphere_dat%gaulats(j)
    coslats_sq(:,j) = 1.-sphere_dat%gaulats(j)**2
 enddo
 
 do k=1,2

!==> vort. flux

    vec1 = ug(:,:,k)*(vrtg(:,:,k)+2.*model_dat%rot*sinlats) + & 
    (float(2-k)/model_dat%tdrag)*vg(:,:,k) + &
    (thtadot*totdelpig/(4.*model_dat%delth*delpig(:,:,k)))* &
    (ug(:,:,2)-ug(:,:,1))
    vec2 = vg(:,:,k)*(vrtg(:,:,k)+2.*model_dat%rot*sinlats) - &
    (float(2-k)/model_dat%tdrag)*ug(:,:,k) - &
    (thtadot*totdelpig/(4.*model_dat%delth*delpig(:,:,k)))* &
    (vg(:,:,2)-vg(:,:,1))
    call rfft(sphere_dat, vec1, coeff1, 1)
    call rfft(sphere_dat, vec2, coeff2, 1)
    call sumnm(sphere_dat,coeff1,coeff2,model_dat%dvrtdtnm(:,k),-1,1)
    call sumnm(sphere_dat,coeff2,coeff1,model_dat%ddivdtnm(:,k),1,-1)

!==> add dissipation to vort. eqn.

    model_dat%dvrtdtnm(:,k) =  &
    model_dat%dvrtdtnm(:,k) + model_dat%dissnm*model_dat%vrtnm(:,k)

!==> mass flux.

    vec1 = ug(:,:,k)*delpig(:,:,k)
    vec2 = vg(:,:,k)*delpig(:,:,k)
    call rfft(sphere_dat, vec1, coeff1, 1)
    call rfft(sphere_dat, vec2, coeff2, 1)
    call sumnm(sphere_dat,coeff1,coeff2,model_dat%ddelpidtnm(:,k),-1,1)

!==> laplacian and dissipation terms to div eqn.

    vec1 = mg(:,:,k) + (0.5*(ug(:,:,k)**2+vg(:,:,k)**2)/coslats_sq)
    call spharm(sphere_dat,vec1,scrnm,1)
    model_dat%ddivdtnm(:,k) = model_dat%ddivdtnm(:,k) + &
    model_dat%dissnm*model_dat%divnm(:,k) - &
    sphere_dat%lap*scrnm

!==> diabatic mass flux.

    vec2 = &
    ((-1)**k*0.5*thtadot/model_dat%delth)*totdelpig

!==> add diabatic term to continuity eqn.

     call spharm(sphere_dat,vec2,scrnm,1)
     model_dat%ddelpidtnm(:,k) = model_dat%ddelpidtnm(:,k) + scrnm
 enddo

end subroutine gettend



subroutine rk4(model_dat,sphere_dat)
! step model_dat structure forward with fourth order runge-kutta scheme.

 type (sphere) :: sphere_dat
 type (modelvars) :: model_dat

 complex, dimension((sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2,2) :: &
  vrtnm,dvrtdtnm1,dvrtdtnm2,&
  divnm,ddivdtnm1,ddivdtnm2,& 
  delpinm,ddelpidtnm1,ddelpidtnm2

  if ( .not. module_initialized ) call initialize_module

 h = model_dat%dt
 hh = 0.5*h
 h6 = h/6.


 vrtnm = model_dat%vrtnm
 divnm = model_dat%divnm
 delpinm = model_dat%delpinm

 call gettend(model_dat,sphere_dat)

 dvrtdtnm1 = model_dat%dvrtdtnm
 ddivdtnm1 = model_dat%ddivdtnm
 ddelpidtnm1 = model_dat%ddelpidtnm


 model_dat%vrtnm=vrtnm+hh*model_dat%dvrtdtnm
 model_dat%divnm=divnm+hh*model_dat%ddivdtnm
 model_dat%delpinm=delpinm+hh*model_dat%ddelpidtnm

 call gettend(model_dat,sphere_dat)


 model_dat%vrtnm=vrtnm+hh*model_dat%dvrtdtnm
 model_dat%divnm=divnm+hh*model_dat%ddivdtnm
 model_dat%delpinm=delpinm+hh*model_dat%ddelpidtnm

 dvrtdtnm2 = model_dat%dvrtdtnm
 ddivdtnm2 = model_dat%ddivdtnm
 ddelpidtnm2 = model_dat%ddelpidtnm

 call gettend(model_dat,sphere_dat)

 model_dat%vrtnm=vrtnm+h*model_dat%dvrtdtnm
 model_dat%divnm=divnm+h*model_dat%ddivdtnm
 model_dat%delpinm=delpinm+h*model_dat%ddelpidtnm

 dvrtdtnm2 = 2.*(dvrtdtnm2 + model_dat%dvrtdtnm)
 ddivdtnm2 = 2.*(ddivdtnm2 + model_dat%ddivdtnm)
 ddelpidtnm2 = 2.*(ddelpidtnm2 + model_dat%ddelpidtnm)


 call gettend(model_dat,sphere_dat)


 model_dat%vrtnm = vrtnm + h6*(dvrtdtnm1+model_dat%dvrtdtnm+dvrtdtnm2)
 model_dat%divnm = divnm + h6*(ddivdtnm1+model_dat%ddivdtnm+ddivdtnm2)
 model_dat%delpinm = delpinm + h6*(ddelpidtnm1+model_dat%ddelpidtnm+ddelpidtnm2)


end subroutine rk4

end module pe2lyr_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
