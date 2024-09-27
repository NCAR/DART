module apm_model_fields_vertloc_mod

   use apm_mapping_mod, only    : w3fb13
  
   implicit none
   private

   public :: vertical_locate, &
             get_model_profile, &
             get_DART_diag_data, &
             handle_err, &
             interp_hori_vert, &
             interp_to_obs

   contains
  
subroutine vertical_locate(prs_loc,prs_obs,nlev_obs,locl_prf,nlay_obs,kend)
!
! This subroutine identifies a vertical location for 
! vertical positioning/localization 
! 
   implicit none
   integer                         :: nlay_obs,nlev_obs
   integer                         :: k,kstr,kmax,kend
   real                            :: prs_loc
   real                            :: wt_ctr,wt_end
   real                            :: zmax
   real,dimension(nlev_obs)        :: prs_obs
   real,dimension(nlay_obs)        :: locl_prf,locl_prf_sm
!
! apply vertical smoother
!   wt_ctr=2.
!   wt_end=1.
!   avgk_sm(:)=0.
!   do k=1,nlay
!      if(k.eq.1) then
!         avgk_sm(k)=(wt_ctr*avgk(k)+wt_end*avgk(k+1))/(wt_ctr+wt_end)
!         cycle
!      elseif(k.eq.nlay) then
!         avgk_sm(k)=(wt_end*avgk(k-1)+wt_ctr*avgk(k))/(wt_ctr+wt_end)
!         cycle
!      else
!         avgk_sm(k)=(wt_end*avgk(k-1)+wt_ctr*avgk(k)+wt_end*avgk(k+1))/(wt_ctr+2.*wt_end)
!      endif
!   enddo
!
   locl_prf_sm(:)=locl_prf(:)   
! locate maximum
   zmax=-1.e10
   kmax=0
   do k=1,kend
      if(abs(locl_prf_sm(k)).gt.zmax) then
         zmax=abs(locl_prf_sm(k))
         kmax=k
      endif
   enddo
   if(kmax.eq.1) kmax=kmax+1
   prs_loc=(prs_obs(kmax)+prs_obs(kmax+1))/2.
end subroutine vertical_locate

!-------------------------------------------------------------------------------

subroutine get_model_profile(prf_locl,prf_full,nz_mdl,prs_obs,prs_mdl, &
   tmp_mdl,qmr_mdl,fld_mdl,nlev_obs,v_wgts,prior,kend)
   implicit none
   integer                                :: nz_mdl
   integer                                :: nlev_obs
   integer                                :: i,j,k,kend
   real                                   :: Ru,Rd,cp,eps,AvogN,msq2cmsq,grav
   real,dimension(nz_mdl)                 :: prs_mdl,tmp_mdl,qmr_mdl,fld_mdl
   real,dimension(nz_mdl)                 :: tmp_prf,vtmp_prf,fld_prf
   real,dimension(nlev_obs-1)             :: thick,v_wgts,prior,prf_locl,prf_full
   real,dimension(nlev_obs)               :: fld_prf_mdl,vtmp_prf_mdl,prs_obs
!
! Constants (mks units)
   Ru=8.316
   Rd=286.9
   cp=1004.
   eps=0.61
   AvogN=6.02214e23
   msq2cmsq=1.e4
   grav=9.8
!
! calculate temperature from potential temperature
   do k=1,nz_mdl
      tmp_prf(k)=tmp_mdl(k)*((prs_mdl(k)/ &
      100000.)**(Rd/cp))
   enddo         
! calculate virtual temperature
   do k=1,nz_mdl
      vtmp_prf(k)=tmp_prf(k)*(1.+eps*qmr_mdl(k))
   enddo         
! convert to molar density         
   do k=1,nz_mdl
      fld_prf(k)=fld_mdl(k)*prs_mdl(k)/Ru/tmp_prf(k)
   enddo
! Vertical interpolation
   fld_prf_mdl(:)=-9999.  
   vtmp_prf_mdl(:)=-9999.   
   call interp_to_obs(fld_prf_mdl,fld_prf,prs_mdl,prs_obs,nz_mdl,nlev_obs,kend)
   call interp_to_obs(vtmp_prf_mdl,vtmp_prf,prs_mdl,prs_obs,nz_mdl,nlev_obs,kend)
!   
! calculate number density times vertical weighting
   prf_locl(:)=-9999.
   prf_full(:)=-9999.
   do k=1,nlev_obs-1
      thick(k)=Rd*(vtmp_prf_mdl(k)+vtmp_prf_mdl(k+1))/2./grav* &
      log(prs_obs(k)/prs_obs(k+1))     
   enddo
!
! apply scattering weights
   do k=1,nlev_obs-1
! full term
      prf_full(k)=thick(k) * (fld_prf_mdl(k)+fld_prf_mdl(k+1))/2.* &
      v_wgts(k) + (1.-v_wgts(k))*prior(k)
!
! no thicknesses      
      prf_locl(k)=(fld_prf_mdl(k)+fld_prf_mdl(k+1))/2.* &
      v_wgts(k) + (1.-v_wgts(k))*prior(k)
   enddo
!   print *, 'prf_full  ',prf_full(:)
!   print *, 'fld fld   ',fld_prf_mdl(:)
!   print *, 'avgk_obs ',v_wgts(:)
end subroutine get_model_profile

!-------------------------------------------------------------------------------

subroutine get_DART_diag_data(file_in,name,data,nx,ny,nz,nc)
   use netcdf
   implicit none
   integer, parameter                    :: maxdim=6
!
   integer                               :: i,icycle,rc,fid,typ,natts
   integer                               :: v_id
   integer                               :: v_ndim
   integer                               :: nx,ny,nz,nc
   integer,dimension(maxdim)             :: one
   integer,dimension(maxdim)             :: v_dimid
   integer,dimension(maxdim)             :: v_dim
!
   character(len=180)                    :: vnam
   character*(*)                         :: name
   character*(*)                         :: file_in
!
   real,dimension(nx,ny,nz,nc)           :: data
!
   ! open netcdf file
   rc = nf90_open(trim(file_in),NF90_NOWRITE,fid)
   if(rc.ne.nf90_noerr) call handle_err(rc,'nf90_open')
!
! get variables identifiers
   rc = nf90_inq_varid(fid,trim(name),v_id)
   if(rc.ne.nf90_noerr) call handle_err(rc,'nf90_inq_varid')
!
! get dimension identifiers
   v_dimid(:)=0
   rc = nf90_inquire_variable(fid,v_id,vnam,typ,v_ndim,v_dimid,natts)
   if(rc.ne.nf90_noerr) call handle_err(rc,'nf90_inquire_variable')
   if(maxdim.lt.v_ndim) then
      print *, 'ERROR: maxdim is too small ',maxdim,v_ndim
      call abort
   endif            
!
! get dimensions
   v_dim(:)=1
   do i=1,v_ndim
      rc = nf90_inquire_dimension(fid,v_dimid(i),len = v_dim(i))
   if(rc.ne.nf90_noerr) call handle_err(rc,'nf90_inquire_dimension')
   enddo
!
! check dimensions
   if(nx.ne.v_dim(1)) then
      print *, 'ERROR: nx dimension conflict ',nx,v_dim(1)
      call abort
   else if(ny.ne.v_dim(2)) then             
      print *, 'ERROR: ny dimension conflict ',ny,v_dim(2)
      call abort
   else if(nz.ne.v_dim(3)) then             
      print *, 'ERROR: nz dimension conflict ',nz,v_dim(3)
      call abort
   else if(nc.ne.v_dim(4)) then             
      print *, 'ERROR: nc dimension conflict ',nc,v_dim(4)
      call abort
   endif
!
! get data
   one(:)=1
   rc = nf90_get_var(fid,v_id,data,one,v_dim)
   if(rc.ne.nf90_noerr) call handle_err(rc,'nf90_get_var')
   rc = nf90_close(fid)
   if(rc.ne.nf90_noerr) call handle_err(rc,'nf90_close')
   return
end subroutine get_DART_diag_data

!-------------------------------------------------------------------------------

subroutine handle_err(rc,text)
   implicit none
   integer         :: rc
   character*(*)   :: text
   print *, 'APM: NETCDF ERROR ',trim(text),' ',rc
   call abort
end subroutine handle_err

!-------------------------------------------------------------------------------

subroutine interp_hori_vert(fld1_prf,fld2_prf,fld1_mdl,fld2_mdl,x_mdl,y_mdl, &
   x_obs,y_obs,prs_mdl,prs_obs,nx_mdl,ny_mdl,nz_mdl,nlev_obs,reject,kend)
   implicit none
   integer                                :: nx_mdl,ny_mdl,nz_mdl,nlev_obs
   integer                                :: i,j,k,i_min,j_min,kend
   integer                                :: im,ip,jm,jp,quad,reject
   real                                   :: re,pi,rad2deg
   real                                   :: rad,rad_crit,rad_min,mdl_x,mdl_y,obs_x,obs_y
   real                                   :: dx_dis,dy_dis
   real                                   :: x_obs,y_obs,x_obs_temp
   real                                   :: x_obser,y_obser
   real                                   :: w_q1,w_q2,w_q3,w_q4,wt
   real,dimension(nlev_obs)               :: fld1_prf,fld2_prf,prs_obs
   real,dimension(nz_mdl)                 :: fld1_prf_mdl,fld2_prf_mdl,prs_prf_mdl
   real,dimension(nx_mdl,ny_mdl)          :: x_mdl,y_mdl
   real,dimension(nx_mdl,ny_mdl,nz_mdl)   :: fld1_mdl,fld2_mdl,prs_mdl
   real                                   :: cone_fac,cen_lat,cen_lon,truelat1, &
                                             truelat2,moad_cen_lat,stand_lon, &
                                             pole_lat,pole_lon,xi,xj,zi,zj,zlon,zlat
   integer                                :: ierr
!
! Set constants
   pi=4.*atan(1.)
   rad2deg=360./(2.*pi)
   re=6371000.
   rad_crit=10000.
   reject=0.
!
   cone_fac=.715567   
   cen_lat=40.
   cen_lon=-97.
   truelat1=33.
   truelat2=45.
   moad_cen_lat=40.0000
   stand_lon=-97.
   pole_lat=90.
   pole_lon=0.
!
! Code to test projection and inverse projection
!   zi=nx_mdl
!   zj=ny_mdl
!   call w3fb13(real(y_mdl(zi,zj)),real(x_mdl(zi,zj)+360.), &
!   real(y_mdl(1,1)),real(x_mdl(1,1)+360.),12000.,cen_lon+360.,truelat1,truelat2,xi,xj)
!   print *,'i,j ',zi,zj
!   print *,'lon lat ',x_mdl(zi,zj)+360.,y_mdl(zi,zj)
!   print *, 'xi,xj ',xi,xj
!
!   call w3fb14(xi,xj,real(y_mdl(1,1)),real(x_mdl(1,1)+360.),12000.,cen_lon+360.,truelat1, &
!   truelat2,zlat,zlon,ierr)
!   print *, 'zlon,zlat ',zlon,zlat
!   print *, 'ierr ',ierr
!
! The input grids need to be in degrees
   x_obser=x_obs
   y_obser=y_obs
   if(x_obser.lt.0.) x_obser=(360.+x_obser)
   call w3fb13(y_obser,x_obser,real(y_mdl(1,1)),real(x_mdl(1,1)+360.), &
   12000.,cen_lon+360.,truelat1,truelat2,xi,xj)
   i_min = nint(xi)
   j_min = nint(xj)
!   print *,' '
!   print *,'i_min,j_min ',i_min,j_min
!   print *,'i_max,j_max ',nx_mdl,ny_mdl
!   x_obs_temp=x_obs
!   if(x_obs.gt.180.) x_obs_temp=x_obs-360.
!   print *,'x_obs,y_obs ',x_obs_temp,y_obs
!   print *,'x_mdl(1,1),x_mdl(1,ny_mdl),x_mdl(nx_mdl,1),x_mdl(nx_mdl,ny_mdl) ', &
!   x_mdl(1,1),x_mdl(1,ny_mdl),x_mdl(nx_mdl,1),x_mdl(nx_mdl,ny_mdl)
!   print *,'y_mdl(1,1),y_mdl(1,ny_mdl),y_mdl(nx_mdl,1),y_mdl(nx_mdl,ny_mdl) ', &
!   y_mdl(1,1),y_mdl(1,ny_mdl),y_mdl(nx_mdl,1),y_mdl(nx_mdl,ny_mdl)
!
! Check lower bounds
   if(i_min.lt.1 .and. int(xi).eq.0) then
      i_min=1
   elseif (i_min.lt.1 .and. int(xi).lt.0) then
      i_min=-9999
      j_min=-9999
      reject=1
   endif
   if(j_min.lt.1 .and. int(xj).eq.0) then
      j_min=1
   elseif (j_min.lt.1 .and. int(xj).lt.0) then
      i_min=-9999
      j_min=-9999
      reject=1
   endif
!
! Check upper bounds
   if(i_min.gt.nx_mdl .and. int(xi).eq.nx_mdl) then
      i_min=nx_mdl
   elseif (i_min.gt.nx_mdl .and. int(xi).gt.nx_mdl) then
      i_min=-9999
      j_min=-9999
      reject=1
   endif
   if(j_min.gt.ny_mdl .and. int(xj).eq.ny_mdl) then
      j_min=ny_mdl
   elseif (j_min.gt.ny_mdl .and. int(xj).gt.ny_mdl) then
      i_min=-9999
      j_min=-9999
      reject=1
   endif
!   print *,'i_min,j_min,reject ',i_min,j_min,reject
   if(reject.eq.1) return  
!
! Use model point closest to the observatkon   
   fld1_prf(:)=-9999.
   fld1_prf_mdl(:)=fld1_mdl(i_min,j_min,:)
   fld2_prf_mdl(:)=fld2_mdl(i_min,j_min,:)
   prs_prf_mdl(:)=prs_mdl(i_min,j_min,:)
!
! Vertical interpolation
   call interp_to_obs(fld1_prf,fld1_prf_mdl,prs_prf_mdl,prs_obs,nz_mdl,nlev_obs,kend)
   call interp_to_obs(fld2_prf,fld2_prf_mdl,prs_prf_mdl,prs_obs,nz_mdl,nlev_obs,kend)
   return
!
! This part of subroutine is not used  
! Do horizontal interpolation
   im=i_min-1
   if(im.eq.0) im=1
   ip=i_min+1
   if(ip.eq.nx_mdl+1) ip=nx_mdl
   jm=j_min-1
   if(jm.eq.0) jm=1
   jp=j_min+1
   if(jp.eq.ny_mdl+1) jp=ny_mdl
!
! Find quadrant
   quad=0
   mdl_x=x_mdl(i_min,j_min)
   if(x_mdl(i_min,j_min).lt.0.) mdl_x=360.+x_mdl(i_min,j_min)
   mdl_y=y_mdl(i_min,j_min)
   if(mdl_x.ge.x_obser.and.mdl_y.ge.y_obser) quad=1 
   if(mdl_x.le.x_obser.and.mdl_y.ge.y_obser) quad=2 
   if(mdl_x.le.x_obser.and.mdl_y.le.y_obser) quad=3 
   if(mdl_x.ge.x_obser.and.mdl_y.le.y_obser) quad=4
   if(quad.eq.0) then
      print *, 'APM:ERROR IN PROCEDURE INTERP_HORIONTAL quad = 0 '
      stop
   endif
!
! Quad 1
   if (quad.eq.1) then
      mdl_x=x_mdl(i_min,j_min)
      if(x_mdl(i_min,j_min).lt.0.) mdl_x=360.+x_mdl(i_min,j_min) 
      dx_dis=abs(x_obser-mdl_x)/rad2deg*cos((y_obser+y_mdl(i_min,j_min))/rad2deg/2.)*re
      dy_dis=abs(y_obser-y_mdl(i_min,j_min))/rad2deg*re
      w_q1=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mdl_x=x_mdl(im,j_min)
      if(x_mdl(im,j_min).lt.0.) mdl_x=360.+x_mdl(im,j_min) 
      dx_dis=abs(x_obser-mdl_x)/rad2deg*cos((y_obser+y_mdl(im,j_min))/rad2deg/2.)*re
      dy_dis=abs(y_obser-y_mdl(im,j_min))/rad2deg*re
      w_q2=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mdl_x=x_mdl(im,jm)
      if(x_mdl(im,jm).lt.0.) mdl_x=360.+x_mdl(im,jm) 
      dx_dis=abs(x_obser-mdl_x)/rad2deg*cos((y_obser+y_mdl(im,jm))/rad2deg/2.)*re
      dy_dis=abs(y_obser-y_mdl(im,jm))/rad2deg*re
      w_q3=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mdl_x=x_mdl(i_min,jm)
      if(x_mdl(i_min,jm).lt.0.) mdl_x=360.+x_mdl(i_min,jm) 
      dx_dis=abs(x_obser-mdl_x)/rad2deg*cos((y_obser+y_mdl(i_min,jm))/rad2deg/2.)*re
      dy_dis=abs(y_obser-y_mdl(i_min,jm))/rad2deg*re
      w_q4=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
! Quad 2
   else if (quad.eq.2) then
      mdl_x=x_mdl(ip,j_min)
      if(x_mdl(ip,j_min).lt.0.) mdl_x=360.+x_mdl(ip,j_min) 
      dx_dis=abs(x_obser-mdl_x)/rad2deg*cos((y_obser+y_mdl(ip,j_min))/rad2deg/2.)*re
      dy_dis=abs(y_obser-y_mdl(ip,j_min))/rad2deg*re
      w_q1=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mdl_x=x_mdl(i_min,j_min)
      if(x_mdl(i_min,j_min).lt.0.) mdl_x=360.+x_mdl(i_min,j_min) 
      dx_dis=abs(x_obser-mdl_x)/rad2deg*cos((y_obser+y_mdl(i_min,j_min))/rad2deg/2.)*re
      dy_dis=abs(y_obser-y_mdl(i_min,j_min))/rad2deg*re
      w_q2=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mdl_x=x_mdl(i_min,jm)
      if(x_mdl(i_min,jm).lt.0.) mdl_x=360.+x_mdl(i_min,jm) 
      dx_dis=abs(x_obser-mdl_x)/rad2deg*cos((y_obser+y_mdl(i_min,jm))/rad2deg/2.)*re
      dy_dis=abs(y_obser-y_mdl(i_min,jm))/rad2deg*re
      w_q3=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mdl_x=x_mdl(ip,jm)
      if(x_mdl(ip,jm).lt.0.) mdl_x=360.+x_mdl(ip,jm) 
      dx_dis=abs(x_obser-mdl_x)/rad2deg*cos((y_obser+y_mdl(ip,jm))/rad2deg/2.)*re
      dy_dis=abs(y_obser-y_mdl(ip,jm))/rad2deg*re
      w_q4=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
! Quad 3
   else if (quad.eq.3) then
      mdl_x=x_mdl(ip,jp)
      if(x_mdl(ip,jp).lt.0.) mdl_x=360.+x_mdl(ip,jp) 
      dx_dis=abs(x_obser-mdl_x)/rad2deg*cos((y_obser+y_mdl(ip,jp))/rad2deg/2.)*re
      dy_dis=abs(y_obser-y_mdl(ip,jp))/rad2deg*re
      w_q1=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mdl_x=x_mdl(i_min,jp)
      if(x_mdl(i_min,jp).lt.0.) mdl_x=360.+x_mdl(i_min,jp) 
      dx_dis=abs(x_obser-mdl_x)/rad2deg*cos((y_obser+y_mdl(i_min,jp))/rad2deg/2.)*re
      dy_dis=abs(y_obser-y_mdl(i_min,jp))/rad2deg*re
      w_q2=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mdl_x=x_mdl(i_min,j_min)
      if(x_mdl(i_min,j_min).lt.0.) mdl_x=360.+x_mdl(i_min,j_min) 
      dx_dis=abs(x_obser-mdl_x)/rad2deg*cos((y_obser+y_mdl(i_min,j_min))/rad2deg/2.)*re
      dy_dis=abs(y_obser-y_mdl(i_min,j_min))/rad2deg*re
      w_q3=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mdl_x=x_mdl(ip,j_min)
      if(x_mdl(ip,j_min).lt.0.) mdl_x=360.+x_mdl(ip,j_min) 
      dx_dis=abs(x_obser-mdl_x)/rad2deg*cos((y_obser+y_mdl(ip,j_min))/rad2deg/2.)*re
      dy_dis=abs(y_obser-y_mdl(ip,j_min))/rad2deg*re
      w_q4=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
! Quad 4
   else if (quad.eq.4) then
      mdl_x=x_mdl(i_min,jp)
      if(x_mdl(i_min,jp).lt.0.) mdl_x=360.+x_mdl(i_min,jp) 
      dx_dis=abs(x_obser-mdl_x)/rad2deg*cos((y_obser+y_mdl(i_min,jp))/rad2deg/2.)*re
      dy_dis=abs(y_obser-y_mdl(i_min,jp))/rad2deg*re
      w_q1=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mdl_x=x_mdl(im,jp)
      if(x_mdl(im,jp).lt.0.) mdl_x=360.+x_mdl(im,jp) 
      dx_dis=abs(x_obser-mdl_x)/rad2deg*cos((y_obser+y_mdl(im,jp))/rad2deg/2.)*re
      dy_dis=abs(y_obser-y_mdl(im,jp))/rad2deg*re
      w_q2=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mdl_x=x_mdl(im,jm)
      if(x_mdl(im,jm).lt.0.) mdl_x=360.+x_mdl(im,jm) 
      dx_dis=abs(x_obser-mdl_x)/rad2deg*cos((y_obser+y_mdl(im,jm))/rad2deg/2.)*re
      dy_dis=abs(y_obser-y_mdl(im,jm))/rad2deg*re
      w_q3=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mdl_x=x_mdl(i_min,j_min)
      if(x_mdl(i_min,j_min).lt.0.) mdl_x=360.+x_mdl(i_min,j_min) 
      dx_dis=abs(x_obser-mdl_x)/rad2deg*cos((y_obser+y_mdl(i_min,j_min))/rad2deg/2.)*re
      dy_dis=abs(y_obser-y_mdl(i_min,j_min))/rad2deg*re
      w_q4=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
   endif
   if(x_obser.ne.x_mdl(i_min,j_min).or.y_obser.ne.y_mdl(i_min,j_min)) then
      wt=1./w_q1+1./w_q2+1./w_q3+1./w_q4
   endif
!
   fld1_prf_mdl(:)=0.
   fld2_prf_mdl(:)=0.
   prs_prf_mdl(:)=0.
   do k=1,nz_mdl
      if(x_obser.eq.x_mdl(i_min,j_min).and.y_obser.eq.y_mdl(i_min,j_min)) then
         fld1_prf_mdl(k)=fld1_mdl(i_min,j_min,k)
         fld2_prf_mdl(k)=fld2_mdl(i_min,j_min,k)
         prs_prf_mdl(k)=prs_mdl(i_min,j_min,k)
      else if(quad.eq.1) then
         fld1_prf_mdl(k)=(1./w_q1*fld1_mdl(i_min,j_min,k)+1./w_q2*fld1_mdl(im,j_min,k)+ &
         1./w_q3*fld1_mdl(im,jm,k)+1./w_q4*fld1_mdl(i_min,jm,k))/wt
         fld2_prf_mdl(k)=(1./w_q1*fld2_mdl(i_min,j_min,k)+1./w_q2*fld2_mdl(im,j_min,k)+ &
         1./w_q3*fld2_mdl(im,jm,k)+1./w_q4*fld2_mdl(i_min,jm,k))/wt
         prs_prf_mdl(k)=(1./w_q1*prs_mdl(i_min,j_min,k)+1./w_q2*prs_mdl(im,j_min,k)+ &
         1./w_q3*prs_mdl(im,jm,k)+1./w_q4*prs_mdl(i_min,jm,k))/wt
      else if(quad.eq.2) then
         fld1_prf_mdl(k)=(1./w_q1*fld1_mdl(ip,j_min,k)+1./w_q2*fld1_mdl(i_min,j_min,k)+ &
         1./w_q3*fld1_mdl(i_min,jm,k)+1./w_q4*fld1_mdl(ip,jm,k))/wt
         fld2_prf_mdl(k)=(1./w_q1*fld2_mdl(ip,j_min,k)+1./w_q2*fld2_mdl(i_min,j_min,k)+ &
         1./w_q3*fld2_mdl(i_min,jm,k)+1./w_q4*fld2_mdl(ip,jm,k))/wt
         prs_prf_mdl(k)=(1./w_q1*prs_mdl(ip,j_min,k)+1./w_q2*prs_mdl(i_min,j_min,k)+ &
         1./w_q3*prs_mdl(i_min,jm,k)+1./w_q4*prs_mdl(ip,jm,k))/wt
      else if(quad.eq.3) then
         fld1_prf_mdl(k)=(1./w_q1*fld1_mdl(ip,jp,k)+1./w_q2*fld1_mdl(i_min,jp,k)+ &
         1./w_q3*fld1_mdl(i_min,j_min,k)+1./w_q4*fld1_mdl(ip,j_min,k))/wt
         fld2_prf_mdl(k)=(1./w_q1*fld2_mdl(ip,jp,k)+1./w_q2*fld2_mdl(i_min,jp,k)+ &
         1./w_q3*fld2_mdl(i_min,j_min,k)+1./w_q4*fld2_mdl(ip,j_min,k))/wt
         prs_prf_mdl(k)=(1./w_q1*prs_mdl(ip,jp,k)+1./w_q2*prs_mdl(i_min,jp,k)+ &
         1./w_q3*prs_mdl(i_min,j_min,k)+1./w_q4*prs_mdl(ip,j_min,k))/wt
      else if(quad.eq.4) then
         fld1_prf_mdl(k)=(1./w_q1*fld1_mdl(i_min,jp,k)+1./w_q2*fld1_mdl(im,jp,k)+ &
         1./w_q3*fld1_mdl(im,j_min,k)+1./w_q4*fld1_mdl(i_min,j_min,k))/wt
         fld2_prf_mdl(k)=(1./w_q1*fld2_mdl(i_min,jp,k)+1./w_q2*fld2_mdl(im,jp,k)+ &
         1./w_q3*fld2_mdl(im,j_min,k)+1./w_q4*fld2_mdl(i_min,j_min,k))/wt
         prs_prf_mdl(k)=(1./w_q1*prs_mdl(i_min,jp,k)+1./w_q2*prs_mdl(im,jp,k)+ &
         1./w_q3*prs_mdl(im,j_min,k)+1./w_q4*prs_mdl(i_min,j_min,k))/wt
      endif 
   enddo
!
! Vertical interpolation
   call interp_to_obs(fld1_prf,fld1_prf_mdl,prs_prf_mdl,prs_obs,nz_mdl,nlev_obs,kend)
   call interp_to_obs(fld2_prf,fld2_prf_mdl,prs_prf_mdl,prs_obs,nz_mdl,nlev_obs,kend)
end subroutine interp_hori_vert

!-------------------------------------------------------------------------------

subroutine interp_to_obs(prf_mdl,fld_mdl,prs_mdl,prs_obs,nz_mdl,nlev_obs,kend)
! Assumes prs_obs and prs_mdl are bottom to top
   implicit none
   integer                          :: nz_mdl,nlev_obs
   integer                          :: k,kk,kend
   real                             :: wt_dw,wt_up 
   real,dimension(nz_mdl)           :: fld_mdl,prs_mdl
   real,dimension(nlev_obs)         :: prs_obs
   real,dimension(nlev_obs)         :: prf_mdl
!
   prf_mdl(:)=-9999.
   kend=-9999
   do k=1,nlev_obs-1
      if((prs_obs(k)+prs_obs(k+1))/2..lt.prs_mdl(nz_mdl) .and. &
      kend.eq.-9999) then
         kend=k
         exit
      endif
   enddo
   do k=1,nlev_obs
      if(prs_obs(k) .gt. prs_mdl(1)) then
         prf_mdl(k)=fld_mdl(1)
         cycle
      endif
      if(prs_obs(k) .lt. prs_mdl(nz_mdl)) then
         prf_mdl(k)=fld_mdl(nz_mdl)
         if(kend.eq.-9999) kend=k-1
         cycle
      endif
      do kk=1,nz_mdl-1
         if(prs_mdl(kk).ge.prs_obs(k) .and. prs_mdl(kk+1).lt.prs_obs(k)) then
            wt_up=log(prs_mdl(kk))-log(prs_obs(k))
            wt_dw=log(prs_obs(k))-log(prs_mdl(kk+1))
            prf_mdl(k)=(wt_up*fld_mdl(kk)+wt_dw*fld_mdl(kk+1))/(wt_dw+wt_up)
            exit
         endif
      enddo               
   enddo
end subroutine interp_to_obs

end module apm_model_fields_vertloc_mod
