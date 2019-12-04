
! FESOM 2 (Finite-volumE Sea ice-Ocean Model)
! multi-resolution ocean general circulation model
! FESOM/fesom2 is licensed under the GNU General Public License v2.0
! Copyright (C) 2018  FESOM team
!
! This source file was taken from  the FESOM V1.4 modules


subroutine fcn_density(t,s,z,rho)
  !
  ! - calculates insitu density as a function of potential temperature
  !   (t is relative to the surface)
  !   using the Jackett and McDougall equation of state (1992?)
  ! Qiang 02,07,2010:  Should this be updated (1995 or 2003)? The current 
  !   version is also different to the international equation of state 
  !   (Unesco 1983). What is the exact reference for this version then? 

  use o_PARAM
  implicit none

  real(kind=8), intent(IN)       :: t, s, z
  real(kind=8), intent(OUT)      :: rho                 
  real(kind=8)                   :: rhopot, bulk

  if(density_linear) then
     rho = -rho0 * 2.0e-4 * t
  else
     bulk = 19092.56 + t*(209.8925 				&
          - t*(3.041638 - t*(-1.852732e-3			&
          - t*(1.361629e-5))))				&
          + s*(104.4077 - t*(6.500517			&
          -  t*(.1553190 - t*(-2.326469e-4))))		&
          + sqrt(s**3)*(-5.587545				&
          + t*(0.7390729 - t*(1.909078e-2)))		&
          - z *(4.721788e-1 + t*(1.028859e-2		&
          + t*(-2.512549e-4 - t*(5.939910e-7))))		&
          - z*s*(-1.571896e-2				&
          - t*(2.598241e-4 + t*(-7.267926e-6)))		&
          - z*sqrt(s**3)					&
          *2.042967e-3 + z*z*(1.045941e-5			&
          - t*(5.782165e-10 - t*(1.296821e-7)))		&
          + z*z*s						&
          *(-2.595994e-7					&
          + t*(-1.248266e-9 + t*(-3.508914e-9)))

     rhopot = ( 999.842594					&
          + t*( 6.793952e-2			&
          + t*(-9.095290e-3			&
          + t*( 1.001685e-4			&
          + t*(-1.120083e-6			&
          + t*( 6.536332e-9)))))			&
          + s*( 0.824493				&
          + t *(-4.08990e-3			&
          + t *( 7.64380e-5			&
          + t *(-8.24670e-7			&
          + t *( 5.38750e-9)))))			&
          + sqrt(s**3)*(-5.72466e-3		&
          + t*( 1.02270e-4			&
          + t*(-1.65460e-6)))			&
          + 4.8314e-4*s**2)
     rho = rhopot / (1.0 + 0.1*z/bulk)
  end if
end subroutine fcn_density
!
!----------------------------------------------------------------------------
!
subroutine fcn_dens0(t,s,rho)
  ! - calculate density at ocean surface.
  ! - if input t is potential temperature, rho is then
  !   potential density (both referenced to surface).
  ! - This routine is used when diagnosing mixed layer thickness.
  !   Using this instead of fcn_density(t,s,0,rho) reduces computation load.

  implicit none
  real(kind=8), intent(IN)       :: t, s
  real(kind=8), intent(OUT)      :: rho                 
  real(kind=8)                    :: tn

  tn=t*1.00024

  rho = ( 999.842594			        & 
       + tn*( 6.793952e-2			&
       + tn*(-9.095290e-3			&
       + tn*( 1.001685e-4			&
       + tn*(-1.120083e-6			&
       + tn*( 6.536332e-9)))))			&
       + s*( 0.824493				&
       + tn *(-4.08990e-3			&
       + tn *( 7.64380e-5			&
       + tn *(-8.24670e-7			&
       + tn *( 5.38750e-9)))))			&
       + sqrt(s**3)*(-5.72466e-3		&
       + tn*( 1.02270e-4			&
       + tn*(-1.65460e-6)))			&
       + 4.8314e-4*s**2)

end subroutine fcn_dens0
!
!----------------------------------------------------------------------------
!
subroutine compute_ref_density
  use o_MESH
  use o_PARAM
  use o_array
  use g_PARFE
  implicit none
  !
  integer         :: n2, n, k
  real(kind=8)    :: T, S, z
  
  ! this should be modified later to get level mean t and s
  
  density_ref=0.0
!  S=34.  
!  T=4.0   
  S=38.5  ! marmara sea deep water salinity
  T=14.5  ! marmara sea deep water temperature

  do n2=1,myDim_nod2d+eDim_nod2d
     do k=1,num_layers_below_nod2d(n2)+1
        n=nod3d_below_nod2d(k,n2)
        z=min(coord_nod3D(3,n),0.)
        call fcn_density(T, S, z, density_ref(n))
     end do
  end do
  if(mype==0) write(*,*) 'Reference density computed.'
end subroutine compute_ref_density
!
!----------------------------------------------------------------------------
!
subroutine compute_density
  use o_MESH
  use o_PARAM
  use o_array
  use g_PARFE
  implicit none
  !
  integer         :: n2, n, k
  real(kind=8)    :: z
  !
  do n2=1,myDim_nod2d+eDim_nod2d
     do k=1,num_layers_below_nod2d(n2)+1
        n=nod3d_below_nod2d(k,n2)
        z=min(coord_nod3D(3,n), 0.0) 
        call fcn_density(tracer(n,1), tracer(n,2), z, density_insitu(n))
     end do
  end do
end subroutine compute_density
!
!----------------------------------------------------------------------------
!
!subroutine compute_pressure
!  ! Qiang, 16.12.2010: add the nonlinear free surface option
!  use o_MESH
!  use o_ELEMENTS
!  use o_PARAM
!  use o_array
!  use o_MATRICES
!  use g_PARFE
!#ifdef use_fullfreesurf
!  use i_therm_parms
!  use i_array
!#endif
!  implicit none
!  !
!  integer         :: i, n, node_lo, node_hi
!  real(kind=8)    :: dens, denscor, z_up, z_lo, wd_ice_eff
!  
!  do n=1,myDim_nod2d+eDim_nod2d    
!     node_hi = nod3D_below_nod2D(1,n)
!     z_up=0.0 ! already re-checked, qiang, 15.03.2011
!#ifdef use_fullfreesurf
!     wd_ice_eff=(rhoice*m_ice(n)+rhosno*m_snow(n))*rho0r
!     hpressure(node_hi)=g*min(wd_ice_eff,max_ice_loading)
!#else
!     hpressure(node_hi)=0.0
!#endif
!     do i=2, num_layers_below_nod2D(n)+1
!        node_lo = nod3D_below_nod2D(i,n) 
!        dens=0.5_8*(density_insitu(node_hi)-density_ref(node_hi) &
!             +density_insitu(node_lo)-density_ref(node_lo))
!        z_lo=coord_nod3d(3,node_lo)
!        hpressure(node_lo)=hpressure(node_hi)+g*(z_up-z_lo)*dens*rho0r
!        node_hi=node_lo
!        z_up=z_lo
!     end do
!  end do
!
!end subroutine compute_pressure
!!
!!----------------------------------------------------------------------------
!!
!subroutine init_pressure_force
!  ! Initialization for PGF calculation on sigma grids
!  ! -- assume that the difference in the number of layers between neighbour nodes
!  ! is not more than one on sigma grids
!  ! -- in case under cavities, assume that the cavity draft does not go deeper than
!  ! the highest bottom depth under each 2d element 
!  ! breaking these assumptions the model can not work!
!  ! Qiang, 9,4,2010
!  use o_mesh
!  use o_elements
!  use o_array
!  use g_parfe
!  implicit none
!
!  integer       :: k, j, kk, elem, n2, n3, mloc(1), lay_el(3), lay_m
!  integer       :: elnodes2(3), bottnodes(3), surfnodes(3), elnodes(6)
!  real(kind=8)  :: zmean, bottdep(3), surfdep(3)
!
!  allocate(dens_interp_nodes(4, max_num_layers-1, myDim_elem2D))
!
!  !find nodes for the (4 points) cubic interpolation 
!  do elem=1,myDim_elem2D                   
!     if(grid_type_elem2d(elem)==0) cycle !not sigma grid part
!     elnodes2=elem2d_nodes(:,elem)
!     lay_el=num_layers_below_nod2d(elnodes2)+1
!     lay_m=minval(lay_el)
!     if(any(lay_el>lay_m)) lay_m=lay_m+1
!     lay_el(1)=min(lay_el(1),lay_m)
!     lay_el(2)=min(lay_el(2),lay_m)
!     lay_el(3)=min(lay_el(3),lay_m)
!     bottnodes(1)=nod3d_below_nod2d(lay_el(1),elnodes2(1))
!     bottnodes(2)=nod3d_below_nod2d(lay_el(2),elnodes2(2))
!     bottnodes(3)=nod3d_below_nod2d(lay_el(3),elnodes2(3))
!     bottdep=coord_nod3d(3,bottnodes)
!     surfnodes=nod3d_below_nod2d(1,elnodes2)
!     surfdep=coord_nod3d(3,surfnodes)
!     do k=1,lay_m-1
!        elnodes(1:3)=nod3d_below_nod2d(k,elnodes2)
!        if(k==lay_m-1) then
!           elnodes(4:6)=bottnodes
!        else
!           elnodes(4:6)=nod3d_below_nod2d(k+1,elnodes2)
!        end if
!
!        zmean=sum(coord_nod3d(3,elnodes))/6.0
!        if(any(surfdep<zmean)) then
!           mloc=minloc(surfdep)
!           dens_interp_nodes(4, k, elem)=mloc(1)+3
!           zmean=minval(surfdep)
!        elseif(any(bottdep>zmean)) then
!           mloc=maxloc(bottdep)
!           dens_interp_nodes(4, k, elem)=mloc(1)    
!           zmean=maxval(bottdep)
!        else
!           dens_interp_nodes(4, k, elem)=0        
!        end if
!
!        do j=1,3
!           n2=elnodes2(j)
!           do kk=1,lay_el(j)
!              n3=nod3d_below_nod2d(kk,n2)
!              if(coord_nod3d(3,n3)<=zmean) then
!                 dens_interp_nodes(j, k, elem)=kk-1 
!                 if(kk==1) dens_interp_nodes(j,k,elem)=1
!                 exit
!              end if
!           end do ! kk   
!        end do ! j
!     end do ! k
!  end do ! elem 
!
!end subroutine init_pressure_force
!!
!!----------------------------------------------------------------------------
!!
!subroutine compute_pressure_force
!  ! Calculating PGF on sigma grids
!  ! -- assume that the difference in the number of layers between neighbour nodes
!  ! is not more than one on sigma grids
!  ! -- in case under cavities, assume that the cavity draft does not go deeper than
!  ! the highest bottom depth under each 2d element 
!  ! breaking these assumptions the model can not work!
!  ! Qiang, 9,4,2010
!  use o_param
!  use o_mesh
!  use o_elements
!  use o_array
!  use g_parfe
!  implicit none
!
!  integer            :: j, k, elem, bn, lay_el(3), lay_m
!  integer            :: elnodes(6), elnodes2(3), bottnodes(3), surfnodes(3)
!  integer            :: n2, ind, flag, s_nodes(4), s_ind(4)
!  real(kind=8)       :: dx2d(3), dy2d(3), z, dz, bottdep(3), surfdep(3)
!  real(kind=8)       :: p_grad(2), intz(2), intratio, intdens(2)
!  real(kind=8)       :: rho(3), rho_grad(2)
!  real(kind=8)       :: s_z(4), s_dens(4), s_H, aux1, aux2, aux(2)
!  real(kind=8)       :: s_dup, s_dlo, a, b, c, d
!
!  do elem=1, myDim_elem2D                  
!     if(grid_type_elem2d(elem)==0) cycle !not sigma grid part
!     elnodes2=elem2d_nodes(:,elem)
!     dx2d=bafux_2d(:,elem)
!     dy2d=bafuy_2d(:,elem)
!     lay_el=num_layers_below_nod2d(elnodes2)+1
!     lay_m=minval(lay_el)
!     if(any(lay_el>lay_m)) lay_m=lay_m+1
!     lay_el(1)=min(lay_el(1),lay_m)
!     lay_el(2)=min(lay_el(2),lay_m)
!     lay_el(3)=min(lay_el(3),lay_m)
!     bottnodes(1)=nod3d_below_nod2d(lay_el(1),elnodes2(1))
!     bottnodes(2)=nod3d_below_nod2d(lay_el(2),elnodes2(2))
!     bottnodes(3)=nod3d_below_nod2d(lay_el(3),elnodes2(3))
!     bottdep=coord_nod3d(3,bottnodes)
!     surfnodes=nod3d_below_nod2d(1,elnodes2)
!     surfdep=coord_nod3d(3,surfnodes)
!     p_grad=0.0
!
!     do k=1,lay_m-1
!        elnodes(1:3)=nod3d_below_nod2d(k,elnodes2)
!        if(k==lay_m-1) then
!           elnodes(4:6)=bottnodes
!        else
!           elnodes(4:6)=nod3d_below_nod2d(k+1,elnodes2)
!        end if
!        
!        bn=dens_interp_nodes(4,k,elem)         
!        if(bn==0) then
!           z=sum(coord_nod3d(3,elnodes))/6.0
!        elseif(bn<4) then
!           z=bottdep(bn)
!        else
!           z=surfdep(bn-3)
!        end if
!
!        do j=1,3
!           n2=elnodes2(j)
!           ind=dens_interp_nodes(j,k,elem)    
!           !cubic-spline interpolation
!           s_ind=(/ind-1, ind, ind+1, ind+2/) 
!           flag=0
!           if(ind==1) then !assume the mesh has at lease 3 layers
!              flag=-1
!              s_ind(1)=1
!           elseif(ind==lay_el(j)-1) then
!              flag=1
!              s_ind(4)=ind+1
!           end if
!           s_nodes=nod3d_below_nod2d(s_ind,n2)
!           s_z=coord_nod3d(3,s_nodes)
!           s_dens=density_insitu(s_nodes)-density_ref(s_nodes)
!           s_H=s_z(3)-s_z(2)
!           aux1=(s_dens(3)-s_dens(2))/s_H
!           ! derivatives computed in a way to get monotonic profile
!           if(flag==0) then
!              aux2=(s_dens(2)-s_dens(1))/(s_z(2)-s_z(1))
!              s_dup=0.0
!              if(aux1*aux2>0.)  s_dup=2.0*aux1*aux2/(aux1+aux2)
!              aux2=(s_dens(4)-s_dens(3))/(s_z(4)-s_z(3))
!              s_dlo=0.0
!              if(aux1*aux2>0.) s_dlo=2.0*aux1*aux2/(aux1+aux2)
!           elseif(flag==1) then
!              aux2=(s_dens(2)-s_dens(1))/(s_z(2)-s_z(1))
!              s_dup=0.0
!              if(aux1*aux2>0.)  s_dup=2.0*aux1*aux2/(aux1+aux2)
!              s_dlo=1.5*aux1-0.5*s_dup
!           else
!              aux2=(s_dens(4)-s_dens(3))/(s_z(4)-s_z(3))
!              s_dlo=0.0
!              if(aux1*aux2>0.) s_dlo=2.0*aux1*aux2/(aux1+aux2)
!              s_dup=1.5*aux1-0.5*s_dlo
!           end if
!           ! cubic polynomial coefficients
!           a=s_dens(2)
!           b=s_dup
!           c=-(2.0*s_dup+s_dlo)/s_H + 3.0*(s_dens(3)-s_dens(2))/s_H**2
!           d=(s_dup+s_dlo)/s_H**2 - 2.0*(s_dens(3)-s_dens(2))/s_H**3
!           ! interploate
!           dz=z-s_z(2)
!           rho(j)=a+b*dz+c*dz**2+d*dz**3
!        end do
!        rho_grad(1)=sum(rho*dx2d) !elementwise density gradient
!        rho_grad(2)=sum(rho*dy2d)
!
!        !next, compute pressure gradient
!        dz=(sum(coord_nod3d(3,elnodes(1:3)))-sum(coord_nod3d(3,elnodes(4:6))))/3.0
!        aux=g*dz*rho_grad/rho0
!        PGF(:,k,elem) = p_grad + aux*0.5    
!        p_grad=p_grad + aux
!
!     end do !k
!  end do ! elem           
!
!end subroutine compute_pressure_force
!
!----------------------------------------------------------------------------
!
function theta(s,t,p,pr)
  ! to compute local potential temperature at pr
  ! using bryden 1973 polynomial for adiabatic lapse rate
  ! and runge-kutta 4-th order integration algorithm.
  ! ref: bryden,h.,1973,deep-sea res.,20,401-408
  ! fofonoff,n.,1977,deep-sea res.,24,489-491
  ! units:
  !       pressure        p        decibars
  !       temperature     t        deg celsius (ipts-68)
  !       salinity        s        (ipss-78)
  !       reference prs   pr       decibars
  !       potential tmp.  theta    deg celsius
  ! checkvalue: theta= 36.89073 c,s=40 (ipss-78),t=40 deg c,
  ! p=10000 decibars,pr=0 decibars
  implicit none
  real*8 			:: theta, s, t, p, pr
  real*8 			:: h, xk, q
  real*8, external	        :: atg

  h = pr - p
  xk = h*atg(s,t,p)
  t = t + 0.5*xk
  q = xk
  p = p + 0.5*h
  xk = h*atg(s,t,p)
  t = t + 0.29289322*(xk-q)
  q = 0.58578644*xk + 0.121320344*q
  xk = h*atg(s,t,p)
  t = t + 1.707106781*(xk-q)
  q = 3.414213562*xk - 4.121320344*q
  p = p + 0.5*h
  xk = h*atg(s,t,p)
  theta = t + (xk-2.0*q)/6.0
  return
end function theta
!
!-----------------------------------------------------------------
!
function atg(s,t,p)
  ! adiabatic temperature gradient deg c per decibar
  ! ref: bryden,h.,1973,deep-sea res.,20,401-408
  ! units:
  !       pressure        p        decibars
  !       temperature     t        deg celsius (ipts-68)
  !       salinity        s        (ipss-78)
  !       adiabatic       atg      deg. c/decibar
  ! checkvalue: atg=3.255976e-4 c/dbar for s=40 (ipss-78),
  ! t=40 deg c,p0=10000 decibars
  implicit none
  real*8  atg, s, t, p, ds

  ds = s - 35.0
  atg = (((-2.1687e-16*t+1.8676e-14)*t-4.6206e-13)*p   &
       +((2.7759e-12*t-1.1351e-10)*ds+((-5.4481e-14*t        &
       +8.733e-12)*t-6.7795e-10)*t+1.8741e-8))*p             &
       +(-4.2393e-8*t+1.8932e-6)*ds                          &
       +((6.6228e-10*t-6.836e-8)*t+8.5258e-6)*t+3.5803e-5

  return
end function atg

