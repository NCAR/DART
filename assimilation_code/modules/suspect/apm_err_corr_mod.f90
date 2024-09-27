module apm_err_corr_mod
   implicit none
   private
   public :: vertical_transform, &
             get_WRFINPUT_land_mask, &
             get_WRFINPUT_lat_lon, &
             get_WRFINPUT_geo_ht, &
             get_WRFCHEM_emiss_data, &
             put_WRFCHEM_emiss_data, &
             get_WRFCHEM_icbc_data, &
             put_WRFCHEM_icbc_data, &
             init_random_seed, &
             init_const_random_seed, &
             apm_pack, &
             apm_unpack, &
             apm_pack_2d, &
             apm_unpack_2d, &
             recenter_factors
   contains

subroutine vertical_transform(A_chem,geo_ht,nx,ny,nz,nz_chem,corr_lngth_vt)
   implicit none
   integer,                               intent(in)   :: nx,ny,nz,nz_chem
   real,                                  intent(in)   :: corr_lngth_vt
   real,dimension(nx,ny,nz),              intent(in)   :: geo_ht
   real,dimension(nx,ny,nz_chem,nz_chem), intent(out)  :: A_chem
!
   integer             :: i,j,k,l,ll
   real                :: vcov
!
   A_chem(:,:,:,:)=0. 
   do k=1,nz_chem
      do l=1,nz_chem
         do i=1,nx
            do j=1,ny
               vcov=1.-abs(geo_ht(i,j,k)-geo_ht(i,j,l))/corr_lngth_vt
               if(vcov.lt.0.) vcov=0.
! row 1         
               if(k.eq.1 .and. l.eq.1) then
                  A_chem(i,j,k,l)=1.
               elseif(k.eq.1 .and. l.gt.1) then
                  A_chem(i,j,k,l)=0.
               endif
! row 2         
               if(k.eq.2 .and. l.eq.1) then
                  A_chem(i,j,k,l)=vcov
               elseif(k.eq.2 .and. l.eq.2) then
                  A_chem(i,j,k,l)=sqrt(1.-A_chem(i,j,k,l-1)*A_chem(i,j,k,l-1))
               elseif (k.eq.2 .and. l.gt.2) then
                  A_chem(i,j,k,l)=0.
               endif
! row 3 and greater         
               if(k.ge.3) then
                  if(l.eq.1) then
                     A_chem(i,j,k,l)=vcov
                  elseif(l.lt.k .and. l.ne.1) then
                     do ll=1,l-1
                        A_chem(i,j,k,l)=A_chem(i,j,k,l)+A_chem(i,j,l,ll)*A_chem(i,j,k,ll)
                     enddo
                     if(A_chem(i,j,l,l).ne.0) A_chem(i,j,k,l)=(vcov-A_chem(i,j,k,l))/A_chem(i,j,l,l)
                  elseif(l.eq.k) then
                     do ll=1,l-1
                        A_chem(i,j,k,l)=A_chem(i,j,k,l)+A_chem(i,j,k,ll)*A_chem(i,j,k,ll)
                     enddo
                     A_chem(i,j,k,l)=sqrt(1.-A_chem(i,j,k,l))
                  endif
               endif
            enddo
         enddo
      enddo
   enddo
end subroutine vertical_transform

!-------------------------------------------------------------------------------

subroutine perturb_fields(chem_fac_mem_old,chem_fac_mem_new, &
lat,lon,A_chem,nx,ny,nz_chem,nchem_spcs,ngrid_corr,sw_corr_tm, &
corr_lngth_hz)

   use apm_utilities_mod,  only :get_dist
  
   implicit none
   integer,                               intent(in)   :: nx,ny,nz_chem
   integer,                               intent(in)   :: nchem_spcs,sw_corr_tm
   integer,                               intent(in)   :: ngrid_corr
   real,                                  intent(in)   :: corr_lngth_hz
   real,dimension(nx,ny),                 intent(in)   :: lat,lon
   real,dimension(nx,ny,nz_chem,nz_chem), intent(in)   :: A_chem
   real,dimension(nx,ny,nz_chem,nchem_spcs),   intent(out)  :: chem_fac_mem_old
   real,dimension(nx,ny,nz_chem,nchem_spcs),   intent(out)  :: chem_fac_mem_new
!
   integer                             :: i,j,k,isp,ii,jj,kk
   integer                             :: ii_str,ii_end,jj_str,jj_end
   real                                :: pi,wgt
   real                                :: u_ran_1,u_ran_2,zdist
   real,allocatable,dimension(:)       :: pert_chem_sum_old,pert_chem_sum_new
   real,allocatable,dimension(:,:)     :: wgt_sum
   real,allocatable,dimension(:,:,:,:) :: pert_chem_old, pert_chem_new
!
! Constants
   pi=4.*atan(1.)
   chem_fac_mem_old(:,:,:,:)=0.
   chem_fac_mem_new(:,:,:,:)=0.
!
! Define horizontal perturbations
   if(sw_corr_tm) then
      allocate(pert_chem_old(nx,ny,nz_chem,nchem_spcs))
      do i=1,nx
         do j=1,ny
            do k=1,nz_chem
               do isp=1,nchem_spcs
                  call random_number(u_ran_1)
                  if(u_ran_1.eq.0.) call random_number(u_ran_1)
                  call random_number(u_ran_2)
                  if(u_ran_2.eq.0.) call random_number(u_ran_2)
                  pert_chem_old(i,j,k,isp)=sqrt(-2.*log(u_ran_1))*cos(2.*pi*u_ran_2)
               enddo
            enddo
         enddo
      enddo
   endif
   allocate(pert_chem_new(nx,ny,nz_chem,nchem_spcs))
   do i=1,nx
      do j=1,ny
         do k=1,nz_chem
            do isp=1,nchem_spcs
               call random_number(u_ran_1)
               if(u_ran_1.eq.0.) call random_number(u_ran_1)
               call random_number(u_ran_2)
               if(u_ran_2.eq.0.) call random_number(u_ran_2)
               pert_chem_new(i,j,k,isp)=sqrt(-2.*log(u_ran_1))*cos(2.*pi*u_ran_2)
            enddo
         enddo
      enddo
   enddo
!
! Impose horizontal correlations
   allocate(wgt_sum(nx,ny))
   if(sw_corr_tm) then
      chem_fac_mem_old(:,:,:,:)=0.
      wgt_sum(:,:)=0.
      do i=1,nx
         do j=1,ny
            ii_str=max(1,i-ngrid_corr)
            ii_end=min(nx,i+ngrid_corr)
            jj_str=max(1,j-ngrid_corr)
            jj_end=min(ny,j+ngrid_corr)
            do ii=ii_str,ii_end
               do jj=jj_str,jj_end
                  zdist=get_dist(lat(ii,jj),lat(i,j),lon(ii,jj),lon(i,j))
                  if(zdist.le.2.0*corr_lngth_hz) then
                     wgt=1./exp(zdist*zdist/corr_lngth_hz/corr_lngth_hz)
                     wgt_sum(i,j)=wgt_sum(i,j)+wgt
                     do k=1,nz_chem
                        do isp=1,nchem_spcs
                           chem_fac_mem_old(i,j,k,isp)=chem_fac_mem_old(i,j,k,isp)+wgt*pert_chem_old(ii,jj,k,isp)
                        enddo
                     enddo
                  endif
               enddo
            enddo
            if(wgt_sum(i,j).gt.0) then
               do k=1,nz_chem
                  do isp=1,nchem_spcs
                     chem_fac_mem_old(i,j,k,isp)=chem_fac_mem_old(i,j,k,isp)/wgt_sum(i,j)
                  enddo
               enddo
            else
               do k=1,nz_chem
                  do isp=1,nchem_spcs
                     chem_fac_mem_old(i,j,k,isp)=pert_chem_old(i,j,k,isp)
                  enddo
               enddo
            endif
         enddo
      enddo
   endif
   chem_fac_mem_new(:,:,:,:)=0.
   wgt_sum(:,:)=0.
   do i=1,nx
      do j=1,ny
         ii_str=max(1,i-ngrid_corr)
         ii_end=min(nx,i+ngrid_corr)
         jj_str=max(1,j-ngrid_corr)
         jj_end=min(ny,j+ngrid_corr)
         do ii=ii_str,ii_end
            do jj=jj_str,jj_end
               zdist=get_dist(lat(ii,jj),lat(i,j),lon(ii,jj),lon(i,j))
               if(zdist.le.2.0*corr_lngth_hz) then
                  wgt=1./exp(zdist*zdist/corr_lngth_hz/corr_lngth_hz)
                  wgt_sum(i,j)=wgt_sum(i,j)+wgt
                  do k=1,nz_chem
                     do isp=1,nchem_spcs
                        chem_fac_mem_new(i,j,k,isp)=chem_fac_mem_new(i,j,k,isp)+wgt*pert_chem_new(ii,jj,k,isp)
                     enddo
                  enddo
               endif
            enddo
         enddo
         if(wgt_sum(i,j).gt.0) then
            do k=1,nz_chem
               do isp=1,nchem_spcs
                  chem_fac_mem_new(i,j,k,isp)=chem_fac_mem_new(i,j,k,isp)/wgt_sum(i,j)
               enddo
            enddo
         else   
            do k=1,nz_chem
               do isp=1,nchem_spcs
                  chem_fac_mem_new(i,j,k,isp)=pert_chem_new(i,j,k,isp)
               enddo
            enddo
         endif
      enddo
   enddo
   deallocate(wgt_sum)
   if(sw_corr_tm) then
      deallocate(pert_chem_old)
   endif
   deallocate(pert_chem_new)
!
! Impose vertical correlations
   if(sw_corr_tm) then
      allocate(pert_chem_sum_old(nz_chem))
      do i=1,nx
         do j=1,ny
            do isp=1,nchem_spcs
               pert_chem_sum_old(:)=0.
               do k=1,nz_chem
                  do kk=1,nz_chem 
                     pert_chem_sum_old(k)=pert_chem_sum_old(k)+A_chem(i,j,k,kk)*chem_fac_mem_old(i,j,kk,isp)
                  enddo
               enddo
               do k=1,nz_chem
                  chem_fac_mem_old(i,j,k,isp)=pert_chem_sum_old(k)
               enddo 
            enddo
         enddo
      enddo
      deallocate(pert_chem_sum_old)
   endif
   allocate(pert_chem_sum_new(nz_chem))
   do i=1,nx
      do j=1,ny
         do isp=1,nz_chem
            pert_chem_sum_new(:)=0.
            do k=1,nz_chem
               do kk=1,nz_chem 
                  pert_chem_sum_new(k)=pert_chem_sum_new(k)+A_chem(i,j,k,kk)*chem_fac_mem_new(i,j,kk,isp)
               enddo
            enddo
         enddo
         do k=1,nz_chem
            chem_fac_mem_new(i,j,k,isp)=pert_chem_sum_new(k)
         enddo
      enddo
   enddo
   deallocate(pert_chem_sum_new)
end subroutine perturb_fields

!-------------------------------------------------------------------------------

subroutine get_WRFINPUT_land_mask(xland,nx,ny)
   implicit none
   include 'netcdf.inc'
   integer, parameter                    :: maxdim=6
   integer                               :: nx,ny
   integer                               :: i,rc
   integer                               :: f_id
   integer                               :: v_id,v_ndim,typ,natts
   integer,dimension(maxdim)             :: one
   integer,dimension(maxdim)             :: v_dimid
   integer,dimension(maxdim)             :: v_dim
   real,dimension(nx,ny)                 :: xland
   character(len=150)                    :: v_nam
   character*(80)                         :: name
   character*(80)                         :: file
!
! open netcdf file
   file='wrfinput_d01.e001'
   name='XLAND'
   rc = nf_open(trim(file),NF_NOWRITE,f_id)
!   print *, trim(file)
   if(rc.ne.0) then
      print *, 'nf_open error ',trim(file)
      stop
   endif
!
! get variables identifiers
   rc = nf_inq_varid(f_id,trim(name),v_id)
!  print *, v_id
   if(rc.ne.0) then
      print *, 'nf_inq_varid error ', v_id
      stop
   endif
!
! get dimension identifiers
   v_dimid=0
   rc = nf_inq_var(f_id,v_id,v_nam,typ,v_ndim,v_dimid,natts)
!   print *, v_dimid
   if(rc.ne.0) then
      print *, 'nf_inq_var error ', v_dimid
      stop
   endif
!
! get dimensions
   v_dim(:)=1
   do i=1,v_ndim
      rc = nf_inq_dimlen(f_id,v_dimid(i),v_dim(i))
   enddo
!   print *, v_dim
   if(rc.ne.0) then
      print *, 'nf_inq_dimlen error ', v_dim
      stop
   endif
!
! check dimensions
   if(nx.ne.v_dim(1)) then
      print *, 'ERROR: nx dimension conflict ',nx,v_dim(1)
      stop
   else if(ny.ne.v_dim(2)) then
      print *, 'ERROR: ny dimension conflict ',ny,v_dim(2)
      stop
   else if(1.ne.v_dim(3)) then             
      print *, 'ERROR: nz dimension conflict ','1',v_dim(3)
      stop
!   else if(1.ne.v_dim(4)) then             
!      print *, 'ERROR: time dimension conflict ',1,v_dim(4)
!      stop
   endif
!
! get data
   one(:)=1
   rc = nf_get_vara_real(f_id,v_id,one,v_dim,xland)
   if(rc.ne.0) then
      print *, 'nf_get_vara_real ', xland(1,1)
      stop
   endif
   rc = nf_close(f_id)
   return
end subroutine get_WRFINPUT_land_mask   

!-------------------------------------------------------------------------------

subroutine get_WRFINPUT_lat_lon(lat,lon,nx,ny)
   implicit none
   include 'netcdf.inc'
   integer, parameter                    :: maxdim=6
   integer                               :: nx,ny
   integer                               :: i,rc
   integer                               :: f_id
   integer                               :: v_id,v_ndim,typ,natts
   integer,dimension(maxdim)             :: one
   integer,dimension(maxdim)             :: v_dimid
  integer,dimension(maxdim)             :: v_dim
   real,dimension(nx,ny)                 :: lat,lon
   character(len=150)                    :: v_nam
   character*(80)                         :: name
   character*(80)                         :: file
!
! open netcdf file
   file='wrfinput_d01.e001'
   name='XLAT'
   rc = nf_open(trim(file),NF_NOWRITE,f_id)
!   print *, trim(file)
   if(rc.ne.0) then
      print *, 'nf_open error ',trim(file)
      stop
   endif
!
! get variables identifiers
   rc = nf_inq_varid(f_id,trim(name),v_id)
!   print *, v_id
   if(rc.ne.0) then
      print *, 'nf_inq_varid error ', v_id
      stop
   endif
!
! get dimension identifiers
   v_dimid=0
   rc = nf_inq_var(f_id,v_id,v_nam,typ,v_ndim,v_dimid,natts)
!   print *, v_dimid
   if(rc.ne.0) then
      print *, 'nf_inq_var error ', v_dimid
      stop
   endif
!
! get dimensions
   v_dim(:)=1
   do i=1,v_ndim
      rc = nf_inq_dimlen(f_id,v_dimid(i),v_dim(i))
   enddo
!   print *, v_dim
   if(rc.ne.0) then
      print *, 'nf_inq_dimlen error ', v_dim
      stop
   endif
!
! check dimensions
   if(nx.ne.v_dim(1)) then
      print *, 'ERROR: nx dimension conflict ',nx,v_dim(1)
      stop
   else if(ny.ne.v_dim(2)) then
      print *, 'ERROR: ny dimension conflict ',ny,v_dim(2)
      stop
   else if(1.ne.v_dim(3)) then             
      print *, 'ERROR: nz dimension conflict ','1',v_dim(3)
      stop
!   else if(1.ne.v_dim(4)) then             
!      print *, 'ERROR: time dimension conflict ',1,v_dim(4)
!      stop
   endif
!
! get data
   one(:)=1
   rc = nf_get_vara_real(f_id,v_id,one,v_dim,lat)
   if(rc.ne.0) then
      print *, 'nf_get_vara_real ', lat(1,1)
      stop
   endif
   
   name='XLONG'
   rc = nf_inq_varid(f_id,trim(name),v_id)
   if(rc.ne.0) then
      print *, 'nf_inq_varid error ', v_id
      stop
   endif
   rc = nf_get_vara_real(f_id,v_id,one,v_dim,lon)
   if(rc.ne.0) then
      print *, 'nf_get_vara_real ', lon(1,1)
      stop
   endif
   rc = nf_close(f_id)
   return
end subroutine get_WRFINPUT_lat_lon

!-------------------------------------------------------------------------------

subroutine get_WRFINPUT_geo_ht(geo_ht,nx,ny,nz,nzp,nmem)
   implicit none
   include 'netcdf.inc'
   integer, parameter                    :: maxdim=6
   integer                               :: k,nx,ny,nz,nzp,nmem
   integer                               :: i,imem,rc
   integer                               :: f_id
   integer                               :: v_id_ph,v_id_phb,v_ndim,typ,natts
   integer,dimension(maxdim)             :: one
   integer,dimension(maxdim)             :: v_dimid
   integer,dimension(maxdim)             :: v_dim
   real,dimension(nx,ny,nzp)             :: ph,phb
   real,dimension(nx,ny,nz)              :: geo_ht
   character(len=150)                    :: v_nam
   character*(80)                        :: name,cmem
   character*(80)                        :: file
!
! Loop through members to find ensemble mean geo_ht
   geo_ht(:,:,:)=0.
   do imem=1,nmem
      if(imem.ge.0.and.imem.lt.10) write(cmem,"('.e00',i1)"),imem
      if(imem.ge.10.and.imem.lt.100) write(cmem,"('.e0',i2)"),imem
      if(imem.ge.100.and.imem.lt.1000) write(cmem,"('.e',i3)"),imem
!
! open netcdf file
      file='wrfinput_d01'//trim(cmem)
      rc = nf_open(trim(file),NF_NOWRITE,f_id)
      if(rc.ne.0) then
         print *, 'nf_open error ',trim(file)
         stop
      endif
!
! get variables identifiers
      name='PH'
      rc = nf_inq_varid(f_id,trim(name),v_id_ph)
      if(rc.ne.0) then
         print *, 'nf_inq_varid error ', v_id_ph
         stop
      endif
      name='PHB'
      rc = nf_inq_varid(f_id,trim(name),v_id_phb)
      if(rc.ne.0) then
         print *, 'nf_inq_varid error ', v_id_phb
         stop
      endif
!
! get dimension identifiers
      v_dimid=0
      rc = nf_inq_var(f_id,v_id_ph,v_nam,typ,v_ndim,v_dimid,natts)
      if(rc.ne.0) then
         print *, 'nf_inq_var error ', v_dimid
         stop
      endif
!
! get dimensions
      v_dim(:)=1
      do i=1,v_ndim
         rc = nf_inq_dimlen(f_id,v_dimid(i),v_dim(i))
      enddo
      if(rc.ne.0) then
         print *, 'nf_inq_dimlen error ', v_dim
         stop
      endif
!
! check dimensions
      if(nx.ne.v_dim(1)) then
         print *, 'ERROR: nx dimension conflict ',nx,v_dim(1)
         stop
      else if(ny.ne.v_dim(2)) then
         print *, 'ERROR: ny dimension conflict ',ny,v_dim(2)
         stop
      else if(nzp.ne.v_dim(3)) then             
         print *, 'ERROR: nzp dimension conflict ','nzp',v_dim(3)
         stop
      endif
!
! get data
      one(:)=1
      rc = nf_get_vara_real(f_id,v_id_ph,one,v_dim,ph)
      if(rc.ne.0) then
         print *, 'nf_get_vara_real ', ph(1,1,1)
         stop
      endif
      rc = nf_get_vara_real(f_id,v_id_phb,one,v_dim,phb)
      if(rc.ne.0) then
         print *, 'nf_get_vara_real ', phb(1,1,1)
         stop
      endif
!
! get mean geo_ht
      do k=1,nz
         geo_ht(:,:,k)=geo_ht(:,:,k)+(ph(:,:,k)+phb(:,:,k)+ph(:,:,k+1)+ &
         phb(:,:,k+1))/2./float(nmem)
      enddo
      rc = nf_close(f_id)
   enddo
end subroutine get_WRFINPUT_geo_ht

!-------------------------------------------------------------------------------

subroutine get_WRFCHEM_emiss_data(file,name,data,nx,ny,nz_chem)
   implicit none
   include 'netcdf.inc'
   integer, parameter                    :: maxdim=6
   integer                               :: nx,ny,nz_chem
   integer                               :: i,rc
   integer                               :: f_id
   integer                               :: v_id,v_ndim,typ,natts
   integer,dimension(maxdim)             :: one
   integer,dimension(maxdim)             :: v_dimid
   integer,dimension(maxdim)             :: v_dim
   real,dimension(nx,ny,nz_chem)         :: data
   character(len=150)                    :: v_nam
   character*(*)                         :: name
   character*(*)                         :: file
!
! open netcdf file
   rc = nf_open(trim(file),NF_SHARE,f_id)
!   print *, trim(file)
   if(rc.ne.0) then
      print *, 'nf_open error in get ',rc, trim(file)
      stop
   endif
!
! get variables identifiers
   rc = nf_inq_varid(f_id,trim(name),v_id)
!   print *, v_id
   if(rc.ne.0) then
      print *, 'nf_inq_varid error ', v_id
      stop
   endif
!
! get dimension identifiers
   v_dimid=0
   rc = nf_inq_var(f_id,v_id,v_nam,typ,v_ndim,v_dimid,natts)
!   print *, v_dimid
   if(rc.ne.0) then
      print *, 'nf_inq_var error ', v_dimid
      stop
   endif
!
! get dimensions
   v_dim(:)=1
   do i=1,v_ndim
      rc = nf_inq_dimlen(f_id,v_dimid(i),v_dim(i))
   enddo
!   print *, v_dim
   if(rc.ne.0) then
      print *, 'nf_inq_dimlen error ', v_dim
      stop
   endif
!
! check dimensions
   if(nx.ne.v_dim(1)) then
      print *, 'ERROR: nx dimension conflict ',nx,v_dim(1)
      stop
   else if(ny.ne.v_dim(2)) then
      print *, 'ERROR: ny dimension conflict ',ny,v_dim(2)
      stop
   else if(nz_chem.ne.v_dim(3)) then             
      print *, 'ERROR: nz_chem dimension conflict ',nz_chem,v_dim(3)
      stop
   else if(1.ne.v_dim(4)) then             
      print *, 'ERROR: time dimension conflict ',1,v_dim(4)
      stop
   endif
!
! get data
   one(:)=1
   rc = nf_get_vara_real(f_id,v_id,one,v_dim,data)
   if(rc.ne.0) then
      print *, 'nf_get_vara_real ', data(1,1,1)
      stop
   endif
   rc = nf_close(f_id)
   return
end subroutine get_WRFCHEM_emiss_data

!-------------------------------------------------------------------------------

subroutine put_WRFCHEM_emiss_data(file,name,data,nx,ny,nz_chem)
   implicit none
   include 'netcdf.inc'
   integer, parameter                    :: maxdim=6
   integer                               :: nx,ny,nz_chem
   integer                               :: i,rc
   integer                               :: f_id
   integer                               :: v_id,v_ndim,typ,natts
   integer,dimension(maxdim)             :: one
   integer,dimension(maxdim)             :: v_dimid
   integer,dimension(maxdim)             :: v_dim
   real,dimension(nx,ny,nz_chem)         :: data
   character(len=150)                    :: v_nam
   character*(*)                         :: name
   character*(*)                         :: file
!
! open netcdf file
   rc = nf_open(trim(file),NF_WRITE,f_id)
   if(rc.ne.0) then
      print *, 'nf_open error in put ',rc, trim(file)
      stop
   endif
!   print *, 'f_id ',f_id
!
! get variables identifiers
    rc = nf_inq_varid(f_id,trim(name),v_id)
!    print *, v_id
    if(rc.ne.0) then
       print *, 'nf_inq_varid error ', v_id
       stop
    endif
!    print *, 'v_id ',v_id
!
! get dimension identifiers
    v_dimid=0
    rc = nf_inq_var(f_id,v_id,v_nam,typ,v_ndim,v_dimid,natts)
!    print *, v_dimid
    if(rc.ne.0) then
       print *, 'nf_inq_var error ', v_dimid
       stop
    endif
!    print *, 'v_ndim, v_dimid ',v_ndim,v_dimid      
!
! get dimensions
    v_dim(:)=1
    do i=1,v_ndim
       rc = nf_inq_dimlen(f_id,v_dimid(i),v_dim(i))
    enddo
!    print *, v_dim
    if(rc.ne.0) then
       print *, 'nf_inq_dimlen error ', v_dim
       stop
    endif
!
! check dimensions
    if(nx.ne.v_dim(1)) then
       print *, 'ERROR: nx dimension conflict ',nx,v_dim(1)
       stop
    else if(ny.ne.v_dim(2)) then
       print *, 'ERROR: ny dimension conflict ',ny,v_dim(2)
       stop
    else if(nz_chem.ne.v_dim(3)) then             
       print *, 'ERROR: nz_chem dimension conflict ',nz_chem,v_dim(3)
       stop
    else if(1.ne.v_dim(4)) then             
       print *, 'ERROR: time dimension conflict ',1,v_dim(4)
       stop
    endif
!
! put data
    one(:)=1
!   rc = nf_close(f_id)
!   rc = nf_open(trim(file),NF_WRITE,f_id)
   rc = nf_put_vara_real(f_id,v_id,one(1:v_ndim),v_dim(1:v_ndim),data)
   if(rc.ne.0) then
      print *, 'nf_put_vara_real return code ',rc
      print *, 'f_id,v_id ',f_id,v_id
      print *, 'one ',one(1:v_ndim)
      print *, 'v_dim ',v_dim(1:v_ndim)
      stop
   endif
   rc = nf_close(f_id)
   return
end subroutine put_WRFCHEM_emiss_data

!-------------------------------------------------------------------------------

subroutine init_random_seed()
   implicit none
   integer, allocatable :: aseed(:)
   integer :: i, n, un, istat, dt(8), pid, t(2), s
   integer(8) :: count, tms, ierr

   call random_seed(size = n)
   allocate(aseed(n))
!
! Fallback to XOR:ing the current time and pid. The PID is
! useful in case one launches multiple instances of the same
! program in parallel.                                                  
   call system_clock(count)
   if (count /= 0) then
      t = transfer(count, t)
   else
      call date_and_time(values=dt)
      tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
           + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
           + dt(3) * 24 * 60 * 60 * 60 * 1000 &
           + dt(5) * 60 * 60 * 1000 &
           + dt(6) * 60 * 1000 + dt(7) * 1000 &
           + dt(8)
      t = transfer(tms, t)
   end if
   s = ieor(t(1), t(2))
!   pid = getpid() + 1099279 ! Add a prime
   call pxfgetpid(pid,ierr)
   s = ieor(s, pid)
   if (n >= 3) then
      aseed(1) = t(1) + 36269
      aseed(2) = t(2) + 72551
      aseed(3) = pid
      if (n > 3) then
         aseed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
      end if
   else
      aseed = s + 37 * (/ (i, i = 0, n - 1 ) /)
   end if
   call random_seed(put=aseed)
end subroutine init_random_seed

!-------------------------------------------------------------------------------

subroutine init_const_random_seed(rank,date)
   implicit none
   integer                          :: rank,date,primes_dim
   integer                          :: n,at,found,i,str
   integer,allocatable,dimension(:) :: primes,aseed
   logical                          :: is_prime
    
   call random_seed(size=n)
   primes_dim=rank*n
   allocate (aseed(n))
   allocate (primes(primes_dim))
   primes(1)=2
   at=2
   found=1
   do
      is_prime=.true.
      do i=1,found
         if(mod(at,primes(i)).eq.0) then
            is_prime=.false.
            exit
         endif
      enddo
      if(is_prime) then
         found=found+1
         primes(found)=at
      endif
      at=at+1
      if(found.eq.primes_dim) then
         exit
      endif
   enddo
   str=(rank-1)*n+1
   do i=str,primes_dim
      aseed(i-str+1)=date*primes(i)
   enddo
   call random_seed(put=aseed)
   deallocate(aseed,primes)
end subroutine init_const_random_seed

!-------------------------------------------------------------------------------

subroutine apm_pack(A_pck,A_unpck,nx,ny,nz,nl)
   implicit none
   integer                      :: nx,ny,nz,nl
   integer                      :: i,j,k,l,idx
   real,dimension(nx,ny,nz,nl)  :: A_unpck
   real,dimension(nx*ny*nz*nl)  :: A_pck
   idx=0
   do l=1,nl
      do k=1,nz
         do j=1,ny
            do i=1,nx
               idx=idx+1
               A_pck(idx)=A_unpck(i,j,k,l)
            enddo
         enddo
      enddo
   enddo
end subroutine apm_pack

!-------------------------------------------------------------------------------

subroutine apm_unpack(A_pck,A_unpck,nx,ny,nz,nl)
   implicit none
   integer                      :: nx,ny,nz,nl
   integer                      :: i,j,k,l,idx
   real,dimension(nx,ny,nz,nl)  :: A_unpck
   real,dimension(nx*ny*nz*nl)  :: A_pck
   idx=0
   do l=1,nl
      do k=1,nz
         do j=1,ny
            do i=1,nx
               idx=idx+1
               A_unpck(i,j,k,l)=A_pck(idx)
            enddo
         enddo
      enddo
   enddo
end subroutine apm_unpack

!-------------------------------------------------------------------------------

subroutine apm_pack_2d(A_pck,A_unpck,nx,ny,nz,nl)
   implicit none
   integer                      :: nx,ny,nz,nl
   integer                      :: i,j,k,l,idx
   real,dimension(nx,ny)        :: A_unpck
   real,dimension(nx*ny)        :: A_pck
   idx=0
   do j=1,ny
      do i=1,nx
         idx=idx+1
         A_pck(idx)=A_unpck(i,j)
      enddo
   enddo
end subroutine apm_pack_2d

!-------------------------------------------------------------------------------

subroutine apm_unpack_2d(A_pck,A_unpck,nx,ny,nz,nl)
   implicit none
   integer                      :: nx,ny,nz,nl
   integer                      :: i,j,k,l,idx
   real,dimension(nx,ny)        :: A_unpck
   real,dimension(nx*ny)        :: A_pck
   idx=0
   do j=1,ny
      do i=1,nx
         idx=idx+1
         A_unpck(i,j)=A_pck(idx)
      enddo
   enddo
end subroutine apm_unpack_2d

!-------------------------------------------------------------------------------

subroutine recenter_factors(chem_fac,nx,ny,nz_chem,nchem_spcs, &
num_mem,sprd_chem)
   implicit none
   integer,           intent(in)       :: nx,ny,nz_chem,nchem_spcs,num_mem
   real,              intent(in)       :: sprd_chem
   real,dimension(nx,ny,nz_chem,nchem_spcs,num_mem),intent(inout) :: chem_fac
   integer                                                        :: i,j,k,isp,imem
   real                                                           :: mean,std
   real,dimension(num_mem)                                        :: mems,pers
!
! Recenter about ensemble mean
   do i=1,nx
      do j=1,ny
         do k=1,nz_chem
            do isp=1,nchem_spcs
               mems(:)=chem_fac(i,j,k,isp,:)
               mean=sum(mems)/real(num_mem)
               pers=(mems-mean)*(mems-mean)
               std=sqrt(sum(pers)/real(num_mem-1))
               do imem=1,num_mem
                  chem_fac(i,j,k,isp,imem)=(chem_fac(i,j,k,isp,imem)-mean)*sprd_chem/std
               enddo
            enddo
         enddo
      enddo
   enddo
end subroutine recenter_factors

!-------------------------------------------------------------------------------

subroutine get_WRFCHEM_icbc_data(file,name,data,nx,ny,nz,nt)
   implicit none
   include 'netcdf.inc'
   integer, parameter                    :: maxdim=6
   integer                               :: nx,ny,nz,nt
   integer                               :: i,rc
   integer                               :: f_id
   integer                               :: v_id,v_ndim,typ,natts
   integer,dimension(maxdim)             :: one
   integer,dimension(maxdim)             :: v_dimid
   integer,dimension(maxdim)             :: v_dim
   real,dimension(nx,ny,nz,nt)           :: data
   character(len=200)                    :: v_nam
   character*(*)                         :: name
   character*(*)                         :: file
!
! open netcdf file
   rc = nf_open(trim(file),NF_SHARE,f_id)
   if(rc.ne.0) then
      print *, 'nf_open error in get ',rc, trim(file)
      stop
   endif
!
! get variables identifiers
   rc = nf_inq_varid(f_id,trim(name),v_id)
!   print *, v_id
   if(rc.ne.0) then
      print *, 'nf_inq_varid error ', v_id
      stop
   endif
!
! get dimension identifiers
   v_dimid=0
   rc = nf_inq_var(f_id,v_id,v_nam,typ,v_ndim,v_dimid,natts)
!   print *, v_dimid
   if(rc.ne.0) then
      print *, 'nf_inq_var error ', v_dimid
      stop
   endif
!
! get dimensions
   v_dim(:)=1
   do i=1,v_ndim
      rc = nf_inq_dimlen(f_id,v_dimid(i),v_dim(i))
   enddo
!   print *, v_dim
   if(rc.ne.0) then
      print *, 'nf_inq_dimlen error ', v_dim
      stop
   endif
!
! check dimensions
   if(nx.ne.v_dim(1)) then
      print *, 'ERROR: nx dimension conflict ',nx,v_dim(1)
      stop
   else if(ny.ne.v_dim(2)) then
      print *, 'ERROR: ny dimension conflict ',ny,v_dim(2)
      stop
   else if(nz.ne.v_dim(3)) then             
      print *, 'ERROR: nz dimension conflict ',nz,v_dim(3)
      stop
   else if(nt.ne.v_dim(4)) then             
      print *, 'ERROR: time dimension conflict ',1,v_dim(4)
      stop
   endif
!
! get data
   one(:)=1
   rc = nf_get_vara_real(f_id,v_id,one,v_dim,data)
   rc = nf_close(f_id)
   return
end subroutine get_WRFCHEM_icbc_data

!-------------------------------------------------------------------------------

subroutine put_WRFCHEM_icbc_data(file,name,data,nx,ny,nz,nt)
   implicit none
   include 'netcdf.inc'
   integer, parameter                    :: maxdim=6
   integer                               :: nx,ny,nz,nt
   integer                               :: i,rc
   integer                               :: f_id
   integer                               :: v_id,v_ndim,typ,natts
   integer,dimension(maxdim)             :: one
   integer,dimension(maxdim)             :: v_dimid
   integer,dimension(maxdim)             :: v_dim
   real,dimension(nx,ny,nz,nt)           :: data
   character(len=200)                    :: v_nam
   character*(*)                         :: name
   character*(*)                         :: file
!
! open netcdf file
   rc = nf_open(trim(file),NF_WRITE,f_id)
   if(rc.ne.0) then
      print *, 'nf_open error in put ',rc, trim(file)
      stop
   endif
!   print *, 'f_id ',f_id
!
! get variables identifiers
   rc = nf_inq_varid(f_id,trim(name),v_id)
!   print *, v_id
   if(rc.ne.0) then
      print *, 'nf_inq_varid error ', v_id
      stop
   endif
!   print *, 'v_id ',v_id
!
! get dimension identifiers
   v_dimid=0
   rc = nf_inq_var(f_id,v_id,v_nam,typ,v_ndim,v_dimid,natts)
!   print *, v_dimid
   if(rc.ne.0) then
      print *, 'nf_inq_var error ', v_dimid
      stop
   endif
!   print *, 'v_ndim, v_dimid ',v_ndim,v_dimid      
!
! get dimensions
   v_dim(:)=1
   do i=1,v_ndim
      rc = nf_inq_dimlen(f_id,v_dimid(i),v_dim(i))
   enddo
!   print *, v_dim
   if(rc.ne.0) then
      print *, 'nf_inq_dimlen error ', v_dim
      stop
   endif
!
! check dimensions
   if(nx.ne.v_dim(1)) then
      print *, 'ERROR: nx dimension conflict ',nx,v_dim(1)
      stop
   else if(ny.ne.v_dim(2)) then
      print *, 'ERROR: ny dimension conflict ',ny,v_dim(2)
      stop
   else if(nz.ne.v_dim(3)) then             
      print *, 'ERROR: nz dimension conflict ',nz,v_dim(3)
      stop
   else if(nt.ne.v_dim(4)) then             
      print *, 'ERROR: time dimension conflict ',1,v_dim(4)
      stop
   endif
!
! put data
   one(:)=1
   rc = nf_put_vara_real(f_id,v_id,one(1:v_ndim),v_dim(1:v_ndim),data)
   rc = nf_close(f_id)
   return
end subroutine put_WRFCHEM_icbc_data

end module apm_err_corr_mod
