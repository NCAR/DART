module apm_upper_bdy_mod
  
   use        types_mod, only     : r8  
   use time_manager_mod, only     : set_date,                   &
                                    set_calendar_type,          &
                                    time_type,                  &
                                    get_time
   use    utilities_mod, only     : timestamp,                  &
                                    register_module,            &
                                    open_file,                  &
                                    close_file,                 &
                                    initialize_utilities,       &
                                    find_namelist_in_file,      &
                                    check_namelist_read,        &
                                    error_handler,              &
                                    E_ERR,                      & 
                                    E_WARN,                     & 
                                    E_MSG,                      &
                                    E_DBG
   use     location_mod, only     : location_type,              &
                                    set_location

   implicit none
   private

   public :: get_upper_bdy_height, &
             get_upper_bdy_fld, &
             get_upper_bdy_single_fld, &
             get_MOZART_INT_DATA, &
             get_MOZART_REAL_DATA, &
             wrf_dart_ubval_interp, &
             apm_get_exo_coldens, &
             apm_get_upvals, &
             apm_interpolate

   contains

!-------------------------------------------------------------------------------

subroutine get_upper_bdy_single_fld(fld,model,data_file,nx,ny,nz,ntim,lon_obs, &
lat_obs,prs_obs,nprs_obs,fld_prf_mdl,date_obs,datesec_obs)  

   integer,                           intent(in)    :: nx,ny,nz,ntim
   integer,                           intent(in)    :: nprs_obs
   real*8,                            intent(in)    :: lon_obs,lat_obs
   real*8,dimension(nprs_obs),        intent(in)    :: prs_obs
   real*8,dimension(nprs_obs),        intent(out)   :: fld_prf_mdl
   integer                                          :: i,j,k,kk,itim
   integer                                          :: indx,jndx,kndx
   integer                                          :: date_obs,datesec_obs
   integer                                          :: itim_sav,year,month,day,hour,minute,second
   type(time_type)                                  :: time_var
   integer                                          :: jdate_obs,jdate_bck,jdate_fwd,yrleft,jday
   integer,dimension(ntim)                          :: date,datesec
   real                                             :: pi,rad2deg
   real                                             :: bck_xwt,fwd_xwt
   real                                             :: bck_ywt,fwd_ywt
   real                                             :: zwt_up,zwt_dw
   real                                             :: twtx,twty,twt
   real                                             :: ztrp_jbck,ztrp_jfwd
   real                                             :: wt_bck,wt_fwd   
   real,dimension(nx)                               :: lon_glb
   real,dimension(ny)                               :: lat_glb
   real,dimension(nz)                               :: prs_glb,ztrp_fld
   real,dimension(nz)                               :: fld_glb_xmym,fld_glb_xpym,fld_glb_xmyp,fld_glb_xpyp
   real,dimension(nx,ny,nz,ntim)                    :: fld_glb
   character(len=120)                               :: string1
   character(len=*)                                 :: fld,model
   character(len=*)                                 :: data_file
   character(len=*), parameter                      :: routine = 'get_upper_bdy_single_fld'
   character(len=*), parameter                      :: source = 'get_upper_bdy_single_fld.f90'
!
!______________________________________________________________________________________________   
!
! Read the upper boundary large scale data (do this once)
!______________________________________________________________________________________________   
!
   pi=4.*atan(1.)
   rad2deg=360./(2.*pi)
   fld_prf_mdl(:)=0.
!
   call get_MOZART_INT_DATA(data_file,'date',ntim,1,1,1,date)
   call get_MOZART_INT_DATA(data_file,'datesec',ntim,1,1,1,datesec)
   call get_MOZART_REAL_DATA(data_file,'lev',nz,1,1,1,prs_glb)
   call get_MOZART_REAL_DATA(data_file,'lat',ny,1,1,1,lat_glb)
   call get_MOZART_REAL_DATA(data_file,'lon',nx,1,1,1,lon_glb)
! mozart
   if(trim(model).eq.'mozart' .or. trim(model).eq.'MOZART') then
      call get_MOZART_REAL_DATA(data_file,trim(fld),nx,ny,nz,ntim,fld_glb)
! waccm
   elseif (trim(model).eq.'waccm' .or. trim(model).eq.'WACCM') then
      call get_MOZART_REAL_DATA(data_file,trim(fld),nx,ny,nz,ntim,fld_glb)
   else
      print *, 'APM: Large scale model type does not exist '
      call exit_all(-77)
   endif
   lon_glb(:)=lon_glb(:)/rad2deg
   lat_glb(:)=lat_glb(:)/rad2deg
!
!______________________________________________________________________________________________   
!
! Find large scale data correspondeing to the observation time
!______________________________________________________________________________________________   
!
   jdate_obs=date_obs*24*60*60+datesec_obs   
   year=date(1)/10000
   yrleft=mod(date(1),10000)
   month=yrleft/100
   day=mod(yrleft,100)
   time_var=set_date(year,month,day,0,0,0)
   call get_time(time_var,second,jday)
   jdate_bck=jday*24*60*60+datesec(1)
!
   year=date(2)/10000
   yrleft=mod(date(2),10000)
   month=yrleft/100
   day=mod(yrleft,100)
   time_var=set_date(year,month,day,0,0,0)
   call get_time(time_var,second,jday)
   jdate_fwd=jday*24*60*60+datesec(2)
!   
   wt_bck=0
   wt_fwd=0
   itim_sav=0
   do itim=1,ntim-1
      if(jdate_obs.gt.jdate_bck .and. jdate_obs.le.jdate_fwd) then
         wt_bck=real(jdate_fwd-jdate_obs)
         wt_fwd=real(jdate_obs-jdate_bck)
         itim_sav=itim
         exit
      endif
      jdate_bck=jdate_fwd
      year=date(itim+1)/10000
      yrleft=mod(date(itim+1),10000)
      month=yrleft/100
      day=mod(yrleft,100)
      time_var=set_date(year,month,day,0,0,0)
      call get_time(time_var,second,jday)
      jdate_fwd=jday*24*60*60+datesec(itim+1)
   enddo
   if(itim_sav.eq.0) then
      write(string1, *) 'APM: upper bdy data not found for this time '
      call error_handler(E_MSG, routine, string1, source)
      call exit_all(-77)
   endif
!______________________________________________________________________________________________   
!
! Find large scale grid box containing the observation location
!______________________________________________________________________________________________   
!
   indx=-9999
   do i=1,nx-1
      if(lon_obs .le. lon_glb(1)) then
         indx=1
         bck_xwt=1.
         fwd_xwt=0.
         twtx=bck_xwt+fwd_xwt
         exit
      elseif(lon_obs .ge. lon_glb(nx)) then
         indx=nx-1
         bck_xwt=0.
         fwd_xwt=1.
         twtx=bck_xwt+fwd_xwt
         exit
      elseif(lon_obs.gt.lon_glb(i) .and. &
         lon_obs.le.lon_glb(i+1)) then
         indx=i
         bck_xwt=lon_glb(i+1)-lon_obs
         fwd_xwt=lon_obs-lon_glb(i)
         twtx=bck_xwt+fwd_xwt
         exit
      endif
   enddo
   if(indx.lt.0) then
      write(string1, *) 'APM: Obs E/W location outside large scale domain'
      call error_handler(E_MSG, routine, string1, source)
      call exit_all(-77)
   endif
!
   jndx=-9999   
   do j=1,ny-1
      if(lat_obs .le. lat_glb(1)) then
         jndx=1
         bck_ywt=1.
         fwd_ywt=0.
         twty=bck_ywt+fwd_ywt
         exit
      elseif(lat_obs .ge. lat_glb(ny)) then
         jndx=ny-1
         bck_ywt=0.
         fwd_ywt=1.
         twty=bck_ywt+fwd_ywt
         exit
      elseif(lat_obs.gt.lat_glb(j) .and. &
         lat_obs.le.lat_glb(j+1)) then
         jndx=j
         bck_ywt=lat_glb(j+1)-lat_obs
         fwd_ywt=lat_obs-lat_glb(j)
         twty=bck_ywt+fwd_ywt
         exit
      endif
   enddo
   if(jndx.lt.0) then
      write(string1, *) 'APM: Obs N/S location outside large scale domain'
      call error_handler(E_MSG, routine, string1, source)
      call exit_all(-77)
   endif
!
!______________________________________________________________________________________________   
!
! Interpolate large scale field to observation location
!______________________________________________________________________________________________   
!
! Temporal
   do k=1,nz
      fld_glb_xmym(k)=(wt_bck*fld_glb(indx,jndx,k,itim_sav) + &
      wt_fwd*fld_glb(indx,jndx,k,itim_sav+1))/(wt_bck+wt_fwd)
      fld_glb_xpym(k)=(wt_bck*fld_glb(indx+1,jndx,k,itim_sav) + &
      wt_fwd*fld_glb(indx+1,jndx,k,itim_sav+1))/(wt_bck+wt_fwd)
      fld_glb_xmyp(k)=(wt_bck*fld_glb(indx,jndx+1,k,itim_sav) + &
      wt_fwd*fld_glb(indx,jndx+1,k,itim_sav+1))/(wt_bck+wt_fwd)
      fld_glb_xpyp(k)=(wt_bck*fld_glb(indx+1,jndx+1,k,itim_sav) + &
      wt_fwd*fld_glb(indx+1,jndx+1,k,itim_sav+1))/(wt_bck+wt_fwd)
   enddo
!
! Horizontal   
   do k=1,nz
      ztrp_jbck=(bck_xwt*fld_glb_xmym(k) + fwd_xwt*fld_glb_xpym(k))/twtx
      ztrp_jfwd=(bck_xwt*fld_glb_xmyp(k) + fwd_xwt*fld_glb_xpyp(k))/twtx
      ztrp_fld(k)=(bck_ywt*ztrp_jbck + fwd_ywt*ztrp_jfwd)/twty      
   enddo
!
! Vertical   
   do k=1,nprs_obs
      kndx=-9999
      do kk=1,nz-1
         if(prs_obs(k).le.prs_glb(kk)) then
            kndx=1
            zwt_up=1.            
            zwt_dw=0.            
            twt=zwt_up+zwt_dw
            exit
         elseif(prs_obs(k).ge.prs_glb(nz)) then
            kndx=nz-1
            zwt_up=0.            
            zwt_dw=1.            
            twt=zwt_up+zwt_dw
            exit
         elseif(prs_obs(k).gt.prs_glb(kk) .and. &
         prs_obs(k).le.prs_glb(kk+1)) then
            kndx=kk
            zwt_up=prs_glb(kk+1)-prs_obs(k)            
            zwt_dw=prs_obs(k)-prs_glb(kk)
            twt=zwt_up+zwt_dw
            exit
         endif
      enddo
      if(kndx.le.0) then
         write(string1, *) 'APM: Obs vertical location outside large scale domain' 
         call error_handler(E_MSG, routine, string1, source)
         call exit_all(-77)
      endif
      fld_prf_mdl(k)=(zwt_up*ztrp_fld(kndx) + zwt_dw*ztrp_fld(kndx+1))/twt
   enddo
end subroutine get_upper_bdy_single_fld

!-------------------------------------------------------------------------------

subroutine get_upper_bdy_height(prs_model,prs_ref,hgt_ref,model,data_file,nx,ny,nz,ntim, &
    lon_obs,lat_obs,alt_obs,date_obs,datesec_obs)  

   integer,                           intent(in)    :: nx,ny,nz,ntim
   integer,                           intent(in)    :: date_obs,datesec_obs
   real*8,                            intent(in)    :: prs_ref,hgt_ref
   real,                              intent(in)    :: alt_obs
   real,                              intent(in)    :: lon_obs,lat_obs
   real*8,                            intent(out)   :: prs_model
   integer                                          :: i,j,k,kk,itim
   integer                                          :: indx,jndx,kndx
   integer                                          :: itim_sav,year,month,day,hour,minute,second
   type(time_type)                                  :: time_var
   integer                                          :: jdate_obs,jdate_bck,jdate_fwd,yrleft,jday
   integer,dimension(ntim)                          :: date,datesec
   real                                             :: pi,rad2deg
   real                                             :: bck_xwt,fwd_xwt
   real                                             :: bck_ywt,fwd_ywt
   real                                             :: zwt_up,zwt_dw
   real                                             :: twtx,twty,twt
   real                                             :: ztrp_jbck,ztrp_jfwd
   real                                             :: wt_bck,wt_fwd
   real                                             :: geoht_model
   real,dimension(nx)                               :: lon_glb
   real,dimension(ny)                               :: lat_glb
   real,dimension(nz)                               :: prs_glb,ztrp_geoht
   real,dimension(nz)                               :: geoht_glb_xmym,geoht_glb_xpym,geoht_glb_xmyp,geoht_glb_xpyp
   real,dimension(nx,ny,nz,ntim)                    :: geoht_glb
   character(len=120)                               :: string1
   character(len=*)                                 :: model
   character(len=*)                                 :: data_file
   character(len=*), parameter                      :: routine = 'get_upper_bdy_height'
   character(len=*), parameter                      :: source = 'get_upper_bdy_height.f90'
!
!______________________________________________________________________________________________   
!
! Read the upper boundary large scale data (do this once)
!______________________________________________________________________________________________   
!
   pi=4.*atan(1.)
   rad2deg=360./(2.*pi)
!
   call get_MOZART_INT_DATA(data_file,'date',ntim,1,1,1,date)
   call get_MOZART_INT_DATA(data_file,'datesec',ntim,1,1,1,datesec)
   call get_MOZART_REAL_DATA(data_file,'lev',nz,1,1,1,prs_glb)
   call get_MOZART_REAL_DATA(data_file,'lat',ny,1,1,1,lat_glb)
   call get_MOZART_REAL_DATA(data_file,'lon',nx,1,1,1,lon_glb)
! mozart
   if(trim(model).eq.'mozart' .or. trim(model).eq.'MOZART') then
      call get_MOZART_REAL_DATA(data_file,'Z3',nx,ny,nz,ntim,geoht_glb)
! waccm
   elseif (trim(model).eq.'waccm' .or. trim(model).eq.'WACCM') then
      call get_MOZART_REAL_DATA(data_file,'Z3',nx,ny,nz,ntim,geoht_glb)
   else
      print *, 'APM: Large scale model type does not exist '
      call exit_all(-77)
   endif

   lon_glb(:)=lon_glb(:)/rad2deg
   lat_glb(:)=lat_glb(:)/rad2deg
!
!______________________________________________________________________________________________   
!
! Find large scale data correspondeing to the observation time
!______________________________________________________________________________________________   
!
   jdate_obs=date_obs*24*60*60+datesec_obs   
   year=date(1)/10000
   yrleft=mod(date(1),10000)
   month=yrleft/100
   day=mod(yrleft,100)
   time_var=set_date(year,month,day,0,0,0)
   print *, 'year.month,day 1 ',year,month,day
   call get_time(time_var,second,jday)
   jdate_bck=jday*24*60*60+datesec(1)
!
   year=date(2)/10000
   yrleft=mod(date(2),10000)
   month=yrleft/100
   day=mod(yrleft,100)
   time_var=set_date(year,month,day,0,0,0)
   print *, 'year.month,day 2 ',year,month,day
   call get_time(time_var,second,jday)
   jdate_fwd=jday*24*60*60+datesec(2)

   print *, 'back forward ',jdate_bck,jdate_fwd
!   
   wt_bck=0
   wt_fwd=0
   itim_sav=0
   do itim=1,ntim-1
      if(jdate_obs.gt.jdate_bck .and. jdate_obs.le.jdate_fwd) then
         wt_bck=real(jdate_fwd-jdate_obs)
         wt_fwd=real(jdate_obs-jdate_bck)
         itim_sav=itim
         exit
      endif
      jdate_bck=jdate_fwd
      year=date(itim+1)/10000
      yrleft=mod(date(itim+1),10000)
      month=yrleft/100
      day=mod(yrleft,100)
      time_var=set_date(year,month,day,0,0,0)
      call get_time(time_var,second,jday)
      jdate_fwd=jday*24*60*60+datesec(itim+1)
   enddo
   if(itim_sav.eq.0) then
      write(string1, *) 'APM: upper bdy data not found for this time '
      call error_handler(E_MSG, routine, string1, source)
      call exit_all(-77)
   endif
!______________________________________________________________________________________________   
!
! Find large scale grid box containing the observation location
!______________________________________________________________________________________________   
!
   indx=-9999
   do i=1,nx-1
      if(lon_obs .le. lon_glb(1)) then
         indx=1
         bck_xwt=1.
         fwd_xwt=0.
         twtx=bck_xwt+fwd_xwt
         exit
      elseif(lon_obs .ge. lon_glb(nx)) then
         indx=nx-1
         bck_xwt=0.
         fwd_xwt=1.
         twtx=bck_xwt+fwd_xwt
         exit
      elseif(lon_obs.gt.lon_glb(i) .and. &
         lon_obs.le.lon_glb(i+1)) then
         indx=i
         bck_xwt=lon_glb(i+1)-lon_obs
         fwd_xwt=lon_obs-lon_glb(i)
         twtx=bck_xwt+fwd_xwt
         exit
      endif
   enddo
   if(indx.lt.0) then
      write(string1, *) 'APM: Obs E/W location outside large scale domain'
      call error_handler(E_MSG, routine, string1, source)
      call exit_all(-77)
   endif
!
   jndx=-9999   
   do j=1,ny-1
      if(lat_obs .le. lat_glb(1)) then
         jndx=1
         bck_ywt=1.
         fwd_ywt=0.
         twty=bck_ywt+fwd_ywt
         exit
      elseif(lat_obs .ge. lat_glb(ny)) then
         jndx=ny-1
         bck_ywt=0.
         fwd_ywt=1.
         twty=bck_ywt+fwd_ywt
         exit
      elseif(lat_obs.gt.lat_glb(j) .and. &
         lat_obs.le.lat_glb(j+1)) then
         jndx=j
         bck_ywt=lat_glb(j+1)-lat_obs
         fwd_ywt=lat_obs-lat_glb(j)
         twty=bck_ywt+fwd_ywt
         exit
      endif
   enddo
   if(jndx.lt.0) then
      write(string1, *) 'APM: Obs N/S location outside large scale domain'
      call error_handler(E_MSG, routine, string1, source)
      call exit_all(-77)
   endif
!
!______________________________________________________________________________________________   
!
! Interpolate large scale field to observation location
!______________________________________________________________________________________________   
!
! Temporal
   do k=1,nz
      geoht_glb_xmym(k)=(wt_bck*geoht_glb(indx,jndx,k,itim_sav) + &
      wt_fwd*geoht_glb(indx,jndx,k,itim_sav+1))/(wt_bck+wt_fwd)
      geoht_glb_xpym(k)=(wt_bck*geoht_glb(indx+1,jndx,k,itim_sav) + &
      wt_fwd*geoht_glb(indx+1,jndx,k,itim_sav+1))/(wt_bck+wt_fwd)
      geoht_glb_xmyp(k)=(wt_bck*geoht_glb(indx,jndx+1,k,itim_sav) + &
      wt_fwd*geoht_glb(indx,jndx+1,k,itim_sav+1))/(wt_bck+wt_fwd)
      geoht_glb_xpyp(k)=(wt_bck*geoht_glb(indx+1,jndx+1,k,itim_sav) + &
      wt_fwd*geoht_glb(indx+1,jndx+1,k,itim_sav+1))/(wt_bck+wt_fwd)
   enddo
!
! Horizontal   
   do k=1,nz
      ztrp_jbck=(bck_xwt*geoht_glb_xmym(k) + fwd_xwt*geoht_glb_xpym(k))/twtx
      ztrp_jfwd=(bck_xwt*geoht_glb_xmyp(k) + fwd_xwt*geoht_glb_xpyp(k))/twtx
      ztrp_geoht(k)=(bck_ywt*ztrp_jbck + fwd_ywt*ztrp_jfwd)/twty
   enddo
!
! Vertical   
   kndx=-9999
   do kk=1,nz-1
      if(alt_obs.le.ztrp_geoht(kk)) then
         kndx=1
         zwt_up=1.            
         zwt_dw=0.            
         twt=zwt_up+zwt_dw
         exit
      elseif(alt_obs.ge.ztrp_geoht(nz)) then
         kndx=nz-1
         zwt_up=0.            
         zwt_dw=1.            
         twt=zwt_up+zwt_dw
         exit
      elseif(alt_obs.gt.ztrp_geoht(kk) .and. &
      alt_obs.le.ztrp_geoht(kk+1)) then
         kndx=kk
         zwt_up=ztrp_geoht(kk+1)-alt_obs            
         zwt_dw=alt_obs-ztrp_geoht(kk)
         twt=zwt_up+zwt_dw
         exit
      endif
   enddo
   if(kndx.le.0) then
      write(string1, *) 'APM: Obs vertical location outside large scale domain' 
      call error_handler(E_MSG, routine, string1, source)
      call exit_all(-77)
   endif
! geoht_model should equal alt_obs
   geoht_model=(zwt_up*ztrp_geoht(kndx) + zwt_dw*ztrp_geoht(kndx+1))/twt
   prs_model=(zwt_up*prs_glb(kndx) + zwt_dw*prs_glb(kndx+1))/twt
 end subroutine get_upper_bdy_height

!-------------------------------------------------------------------------------
      
subroutine get_upper_bdy_fld(fld,model,data_file,nx,ny,nz,ntim,lon_obs,lat_obs,prs_obs,nprs_obs, &
fld_prf_mdl,tmp_prf_mdl,qmr_prf_mdl,date_obs,datesec_obs)  

   integer,                           intent(in)    :: nx,ny,nz,ntim
   integer,                           intent(in)    :: nprs_obs
   real*8,                            intent(in)    :: lon_obs,lat_obs
   real*8,dimension(nprs_obs),        intent(in)    :: prs_obs
   real*8,dimension(nprs_obs),        intent(out)   :: fld_prf_mdl,tmp_prf_mdl,qmr_prf_mdl
   integer                                          :: i,j,k,kk,itim
   integer                                          :: indx,jndx,kndx
   integer                                          :: date_obs,datesec_obs
   integer                                          :: itim_sav,year,month,day,hour,minute,second
   type(time_type)                                  :: time_var
   integer                                          :: jdate_obs,jdate_bck,jdate_fwd,yrleft,jday
   integer,dimension(ntim)                          :: date,datesec
   real                                             :: pi,rad2deg
   real                                             :: bck_xwt,fwd_xwt
   real                                             :: bck_ywt,fwd_ywt
   real                                             :: zwt_up,zwt_dw
   real                                             :: twtx,twty,twt
   real                                             :: ztrp_jbck,ztrp_jfwd
   real                                             :: wt_bck,wt_fwd   
   real,dimension(nx)                               :: lon_glb
   real,dimension(ny)                               :: lat_glb
   real,dimension(nz)                               :: prs_glb,ztrp_fld,ztrp_tmp,ztrp_qmr
   real,dimension(nz)                               :: fld_glb_xmym,fld_glb_xpym,fld_glb_xmyp,fld_glb_xpyp
   real,dimension(nz)                               :: tmp_glb_xmym,tmp_glb_xpym,tmp_glb_xmyp,tmp_glb_xpyp
   real,dimension(nz)                               :: qmr_glb_xmym,qmr_glb_xpym,qmr_glb_xmyp,qmr_glb_xpyp
   real,dimension(nx,ny,nz,ntim)                    :: fld_glb,tmp_glb,qmr_glb
   character(len=120)                               :: string1
   character(len=*)                                 :: fld,model
   character(len=*)                                 :: data_file
   character(len=*), parameter                      :: routine = 'get_upper_bdy_fld'
   character(len=*), parameter                      :: source = 'get_upper_bdy_fld.f90'
!
!______________________________________________________________________________________________   
!
! Read the upper boundary large scale data (do this once)
!______________________________________________________________________________________________   
!
   pi=4.*atan(1.)
   rad2deg=360./(2.*pi)
   fld_prf_mdl(:)=0.
   tmp_prf_mdl(:)=0.
   qmr_prf_mdl(:)=0.
!
   call get_MOZART_INT_DATA(data_file,'date',ntim,1,1,1,date)
   call get_MOZART_INT_DATA(data_file,'datesec',ntim,1,1,1,datesec)
   call get_MOZART_REAL_DATA(data_file,'lev',nz,1,1,1,prs_glb)
   call get_MOZART_REAL_DATA(data_file,'lat',ny,1,1,1,lat_glb)
   call get_MOZART_REAL_DATA(data_file,'lon',nx,1,1,1,lon_glb)
! mozart
   if(trim(model).eq.'mozart' .or. trim(model).eq.'MOZART') then
      call get_MOZART_REAL_DATA(data_file,trim(fld),nx,ny,nz,ntim,fld_glb)
! waccm
   elseif (trim(model).eq.'waccm' .or. trim(model).eq.'WACCM') then
      call get_MOZART_REAL_DATA(data_file,trim(fld),nx,ny,nz,ntim,fld_glb)
   else
      print *, 'APM: Large scale model type does not exist '
      call exit_all(-77)
   endif
   call get_MOZART_REAL_DATA(data_file,'T',nx,ny,nz,ntim,tmp_glb)
   call get_MOZART_REAL_DATA(data_file,'Q',nx,ny,nz,ntim,qmr_glb)
   lon_glb(:)=lon_glb(:)/rad2deg
   lat_glb(:)=lat_glb(:)/rad2deg
!
!______________________________________________________________________________________________   
!
! Find large scale data correspondeing to the observation time
!______________________________________________________________________________________________   
!
   jdate_obs=date_obs*24*60*60+datesec_obs   
   year=date(1)/10000
   yrleft=mod(date(1),10000)
   month=yrleft/100
   day=mod(yrleft,100)
   time_var=set_date(year,month,day,0,0,0)
   call get_time(time_var,second,jday)
   jdate_bck=jday*24*60*60+datesec(1)
!
   year=date(2)/10000
   yrleft=mod(date(2),10000)
   month=yrleft/100
   day=mod(yrleft,100)
   time_var=set_date(year,month,day,0,0,0)
   call get_time(time_var,second,jday)
   jdate_fwd=jday*24*60*60+datesec(2)
!   
   wt_bck=0
   wt_fwd=0
   itim_sav=0
   do itim=1,ntim-1
      if(jdate_obs.gt.jdate_bck .and. jdate_obs.le.jdate_fwd) then
         wt_bck=real(jdate_fwd-jdate_obs)
         wt_fwd=real(jdate_obs-jdate_bck)
         itim_sav=itim
         exit
      endif
      jdate_bck=jdate_fwd
      year=date(itim+1)/10000
      yrleft=mod(date(itim+1),10000)
      month=yrleft/100
      day=mod(yrleft,100)
      time_var=set_date(year,month,day,0,0,0)
      call get_time(time_var,second,jday)
      jdate_fwd=jday*24*60*60+datesec(itim+1)
   enddo
   if(itim_sav.eq.0) then
      write(string1, *) 'APM: upper bdy data not found for this time '
      call error_handler(E_MSG, routine, string1, source)
      call exit_all(-77)
   endif
!______________________________________________________________________________________________   
!
! Find large scale grid box containing the observation location
!______________________________________________________________________________________________   
!
   indx=-9999
   do i=1,nx-1
      if(lon_obs .le. lon_glb(1)) then
         indx=1
         bck_xwt=1.
         fwd_xwt=0.
         twtx=bck_xwt+fwd_xwt
         exit
      elseif(lon_obs .ge. lon_glb(nx)) then
         indx=nx-1
         bck_xwt=0.
         fwd_xwt=1.
         twtx=bck_xwt+fwd_xwt
         exit
      elseif(lon_obs.gt.lon_glb(i) .and. &
         lon_obs.le.lon_glb(i+1)) then
         indx=i
         bck_xwt=lon_glb(i+1)-lon_obs
         fwd_xwt=lon_obs-lon_glb(i)
         twtx=bck_xwt+fwd_xwt
         exit
      endif
   enddo
   if(indx.lt.0) then
      write(string1, *) 'APM: Obs E/W location outside large scale domain'
      call error_handler(E_MSG, routine, string1, source)
      call exit_all(-77)
   endif
!
   jndx=-9999   
   do j=1,ny-1
      if(lat_obs .le. lat_glb(1)) then
         jndx=1
         bck_ywt=1.
         fwd_ywt=0.
         twty=bck_ywt+fwd_ywt
         exit
      elseif(lat_obs .ge. lat_glb(ny)) then
         jndx=ny-1
         bck_ywt=0.
         fwd_ywt=1.
         twty=bck_ywt+fwd_ywt
         exit
      elseif(lat_obs.gt.lat_glb(j) .and. &
         lat_obs.le.lat_glb(j+1)) then
         jndx=j
         bck_ywt=lat_glb(j+1)-lat_obs
         fwd_ywt=lat_obs-lat_glb(j)
         twty=bck_ywt+fwd_ywt
         exit
      endif
   enddo
   if(jndx.lt.0) then
      write(string1, *) 'APM: Obs N/S location outside large scale domain'
      call error_handler(E_MSG, routine, string1, source)
      call exit_all(-77)
   endif
!
!______________________________________________________________________________________________   
!
! Interpolate large scale field to observation location
!______________________________________________________________________________________________   
!
! Temporal
   do k=1,nz
      fld_glb_xmym(k)=(wt_bck*fld_glb(indx,jndx,k,itim_sav) + &
      wt_fwd*fld_glb(indx,jndx,k,itim_sav+1))/(wt_bck+wt_fwd)
      fld_glb_xpym(k)=(wt_bck*fld_glb(indx+1,jndx,k,itim_sav) + &
      wt_fwd*fld_glb(indx+1,jndx,k,itim_sav+1))/(wt_bck+wt_fwd)
      fld_glb_xmyp(k)=(wt_bck*fld_glb(indx,jndx+1,k,itim_sav) + &
      wt_fwd*fld_glb(indx,jndx+1,k,itim_sav+1))/(wt_bck+wt_fwd)
      fld_glb_xpyp(k)=(wt_bck*fld_glb(indx+1,jndx+1,k,itim_sav) + &
      wt_fwd*fld_glb(indx+1,jndx+1,k,itim_sav+1))/(wt_bck+wt_fwd)
!
      tmp_glb_xmym(k)=(wt_bck*tmp_glb(indx,jndx,k,itim_sav) + &
      wt_fwd*tmp_glb(indx,jndx,k,itim_sav+1))/(wt_bck+wt_fwd)
      tmp_glb_xpym(k)=(wt_bck*tmp_glb(indx+1,jndx,k,itim_sav) + &
      wt_fwd*tmp_glb(indx+1,jndx,k,itim_sav+1))/(wt_bck+wt_fwd)
      tmp_glb_xmyp(k)=(wt_bck*tmp_glb(indx,jndx+1,k,itim_sav) + &
      wt_fwd*tmp_glb(indx,jndx+1,k,itim_sav+1))/(wt_bck+wt_fwd)
      tmp_glb_xpyp(k)=(wt_bck*tmp_glb(indx+1,jndx+1,k,itim_sav) + &
      wt_fwd*tmp_glb(indx+1,jndx+1,k,itim_sav+1))/(wt_bck+wt_fwd)
!
      qmr_glb_xmym(k)=(wt_bck*qmr_glb(indx,jndx,k,itim_sav) + &
      wt_fwd*qmr_glb(indx,jndx,k,itim_sav+1))/(wt_bck+wt_fwd)
      qmr_glb_xpym(k)=(wt_bck*qmr_glb(indx+1,jndx,k,itim_sav) + &
      wt_fwd*qmr_glb(indx+1,jndx,k,itim_sav+1))/(wt_bck+wt_fwd)
      qmr_glb_xmyp(k)=(wt_bck*qmr_glb(indx,jndx+1,k,itim_sav) + &
      wt_fwd*qmr_glb(indx,jndx+1,k,itim_sav+1))/(wt_bck+wt_fwd)
      qmr_glb_xpyp(k)=(wt_bck*qmr_glb(indx+1,jndx+1,k,itim_sav) + &
      wt_fwd*qmr_glb(indx+1,jndx+1,k,itim_sav+1))/(wt_bck+wt_fwd)
   enddo
!
! Horizontal   
   do k=1,nz
      ztrp_jbck=(bck_xwt*fld_glb_xmym(k) + fwd_xwt*fld_glb_xpym(k))/twtx
      ztrp_jfwd=(bck_xwt*fld_glb_xmyp(k) + fwd_xwt*fld_glb_xpyp(k))/twtx
      ztrp_fld(k)=(bck_ywt*ztrp_jbck + fwd_ywt*ztrp_jfwd)/twty
!      
      ztrp_jbck=(bck_xwt*tmp_glb_xmym(k) + fwd_xwt*tmp_glb_xpym(k))/twtx
      ztrp_jfwd=(bck_xwt*tmp_glb_xmyp(k) + fwd_xwt*tmp_glb_xpyp(k))/twtx
      ztrp_tmp(k)=(bck_ywt*ztrp_jbck + fwd_ywt*ztrp_jfwd)/twty
!      
      ztrp_jbck=(bck_xwt*qmr_glb_xmym(k) + fwd_xwt*qmr_glb_xmyp(k))/twtx
      ztrp_jfwd=(bck_xwt*qmr_glb_xmyp(k) + fwd_xwt*qmr_glb_xpyp(k))/twtx
      ztrp_qmr(k)=(bck_ywt*ztrp_jbck + fwd_ywt*ztrp_jfwd)/twty      
   enddo
!
! Vertical   
   do k=1,nprs_obs
      kndx=-9999
      do kk=1,nz-1
         if(prs_obs(k).le.prs_glb(kk)) then
            kndx=1
            zwt_up=1.            
            zwt_dw=0.            
            twt=zwt_up+zwt_dw
            exit
         elseif(prs_obs(k).ge.prs_glb(nz)) then
            kndx=nz-1
            zwt_up=0.            
            zwt_dw=1.            
            twt=zwt_up+zwt_dw
            exit
         elseif(prs_obs(k).gt.prs_glb(kk) .and. &
         prs_obs(k).le.prs_glb(kk+1)) then
            kndx=kk
            zwt_up=prs_glb(kk+1)-prs_obs(k)            
            zwt_dw=prs_obs(k)-prs_glb(kk)
            twt=zwt_up+zwt_dw
            exit
         endif
      enddo
      if(kndx.le.0) then
         write(string1, *) 'APM: Obs vertical location outside large scale domain' 
         call error_handler(E_MSG, routine, string1, source)
         call exit_all(-77)
      endif
      fld_prf_mdl(k)=(zwt_up*ztrp_fld(kndx) + zwt_dw*ztrp_fld(kndx+1))/twt
      tmp_prf_mdl(k)=(zwt_up*ztrp_tmp(kndx) + zwt_dw*ztrp_tmp(kndx+1))/twt
      qmr_prf_mdl(k)=(zwt_up*ztrp_qmr(kndx) + zwt_dw*ztrp_qmr(kndx+1))/twt
   enddo
 end subroutine get_upper_bdy_fld

!-------------------------------------------------------------------------------

subroutine get_MOZART_INT_DATA(file,name,nx,ny,nz,ntim,fld)
   implicit none
   include 'netcdf.inc'
   integer,parameter                                :: maxdim=7000
   integer                                          :: nx,ny,nz,ntim
   integer                                          :: i,rc
   integer                                          :: f_id
   integer                                          :: v_id,v_ndim,typ,natts
   integer,dimension(maxdim)                        :: one
   integer,dimension(maxdim)                        :: v_dimid
   integer,dimension(maxdim)                        :: v_dim
   integer,dimension(ntim)                          :: fld
   character(len=*)                                 :: file
   character(len=*)                                 :: name
   character(len=120)                               :: v_nam
!
! open netcdf data file
   rc = nf_open(trim(file),NF_NOWRITE,f_id)
!
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
   else if(nz.ne.v_dim(3)) then             
      print *, 'ERROR: nz dimension conflict ','1',v_dim(3)
      stop
   else if(ntim.ne.v_dim(4)) then             
      print *, 'ERROR: time dimension conflict ',1,v_dim(4)
      stop
   endif
!
! get data
   one(:)=1
   rc = nf_get_vara_int(f_id,v_id,one,v_dim,fld)
   if(rc.ne.0) then
      print *, 'nf_get_vara_real ', fld(1)
      stop
   endif
   rc = nf_close(f_id)
   return
     
end subroutine get_MOZART_INT_DATA

!-------------------------------------------------------------------------------

subroutine get_MOZART_REAL_DATA(file,name,nx,ny,nz,ntim,fld)
   implicit none
   include 'netcdf.inc'   
   integer,parameter                                :: maxdim=7000
   integer                                          :: nx,ny,nz,ntim
   integer                                          :: i,rc
   integer                                          :: f_id
   integer                                          :: v_id,v_ndim,typ,natts
   integer,dimension(maxdim)                        :: one
   integer,dimension(maxdim)                        :: v_dimid
   integer,dimension(maxdim)                        :: v_dim
   real,dimension(nx,ny,nz,ntim)                    :: fld
   character(len=*)                                 :: file
   character(len=*)                                 :: name
   character(len=120)                               :: v_nam
!
! open netcdf data file
   rc = nf_open(trim(file),NF_NOWRITE,f_id)
!   print *, 'f_id ',f_id
!
   if(rc.ne.0) then
      print *, 'nf_open error ',trim(file)
      stop
   endif
!
! get variables identifiers
   rc = nf_inq_varid(f_id,trim(name),v_id)
!   print *, 'v_id ',v_id
!
   if(rc.ne.0) then
      print *, 'nf_inq_varid error ', v_id
      stop
   endif
   !
! get dimension identifiers
   v_dimid=0
   rc = nf_inq_var(f_id,v_id,v_nam,typ,v_ndim,v_dimid,natts)
!   print *, 'v_dimid ',v_dimid
!
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
!   print *, 'v_dim ',v_dim
!
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
      print *, 'ERROR: nz dimension conflict ','1',v_dim(3)
      stop
   else if(ntim.ne.v_dim(4)) then             
      print *, 'ERROR: time dimension conflict ',1,v_dim(4)
      stop
   endif
!
! get data
   one(:)=1
   rc = nf_get_vara_real(f_id,v_id,one,v_dim,fld)
!   print *, 'fld ', fld(1,1,1,1),fld(nx/2,ny/2,nz/2,ntim/2),fld(nx,ny,nz,ntim)
!
   if(rc.ne.0) then
      print *, 'nf_get_vara_real ', fld(1,1,1,1)
      stop
   endif
   rc = nf_close(f_id)
   return
     
end subroutine get_MOZART_REAL_DATA
!
! APM FROM IASI
subroutine wrf_dart_ubval_interp(obs_val,del_prs,domain,species,lon,lat,lev,im2,istatus)
   use netcdf
   use        types_mod, only     : r8  
   implicit none
   integer,parameter                                 :: nx1=2,ny1=96,nz1=38,nm1=12,nmspc1=8,nchr1=20
   integer,parameter                                 :: nx2=100,ny2=40,nz2=66,nm2=12
   integer                                           :: fid1,fid2,domain,rc,im2
   integer                                           :: i,j,k,imn1,istatus
   character(len=20)                                 :: species
   character(len=180)                                :: file_nam,file_in1,file_in2,path
   real(r8)                                          :: lon,lat,lev,del_lon1
   real(r8)                                          :: obs_val,del_prs    
   real,dimension(nx1,ny1)                           :: xlon1,xlat1
   real,dimension(ny1)                               :: xlat_tmp1
   real,dimension(nz1)                               :: xlev1,prs_tmp1
   real,dimension(ny1,nz1)                           :: fld1,fld_tmp1
   real,dimension(nx1,ny1,nz1)                       :: fldd1
   real,dimension(nx2,ny2)                           :: xlat2,xlon2
   real,dimension(nz2)                               :: xlev2,prs_tmp2
   real,dimension(nm2)                               :: ddoyr
   real,dimension(nx2,ny2,nz2,nm2)                   :: o3_col_dens
   real,dimension(nx2,ny2,nz2)                       :: fld_tmp2
   logical                                           :: use_interp_1
!
! decide with upper boundary data to use
   use_interp_1=.true.
   del_lon1=2.
!
! assign upper boundary profile files
   path='./'
!
! this file has OX NOX HNO3 CH4 CO N2X N2O5 H20
   file_in1='ubvals_b40.20th.track1_1996-2005.nc'
!
! this file has o3 only 
   file_in2='exo_coldens_d01'
   if(domain.eq.2) then
      file_in2='exo_coldens_d02'
   endif
!
! open upper boundary profile files
   if (use_interp_1) then
      file_nam=trim(path)//trim(file_in1)
      rc = nf90_open(trim(file_nam),NF90_NOWRITE,fid1)
      if(rc.ne.0) then
         print *, 'APM: nc_open error file=',trim(file_nam)
         call abort
      endif
!      print *, 'opened ',trim(file_nam)
   else
      file_nam=trim(path)//trim(file_in2)
      rc = nf90_open(trim(file_nam),NF90_NOWRITE,fid2)
      if(rc.ne.0) then
         print *, 'APM: nc_open error file=',trim(file_nam)
         call abort
      endif
!      print *, 'opened ',trim(file_nam)
   endif
!
! select upper boundary data from ubvals_b40.20th.track1_1996-2005.nc
   if (use_interp_1) then
      imn1=6
      call apm_get_upvals(fid1,species,imn1,fld1,xlat_tmp1,xlev1)
      rc=nf90_close(fid1)
   else
!
! select upper boundary data from exo_coldens_dxx
      call apm_get_exo_coldens(fid2,'XLAT',xlat2,nx2,ny2,1,1)
!      print *, 'XLAT',xlat2(1,1),xlat2(nx2,ny2)
      call apm_get_exo_coldens(fid2,'XLONG',xlon2,nx2,ny2,1,1)
!      print *, 'XLON',xlon2(1,1),xlon2(nx2,ny2)
      call apm_get_exo_coldens(fid2,'coldens_levs',xlev2,nz2,1,1,1)
!      print *, 'coldens_levs',xlev2(:)
      call apm_get_exo_coldens(fid2,'days_of_year',ddoyr,nm2,1,1,1)
!      print *, 'ddoyr',ddoyr(1),ddoyr(nm2)
      call apm_get_exo_coldens(fid2,'o3_column_density',o3_col_dens,nx2,ny2,nz2,nm2)
!      print *, 'o3_coldens',o3_col_dens(1,1,1,1),o3_col_dens(nx2,ny2,nz2,nm2)
      rc=nf90_close(fid2)
   endif
!   print *, 'ny1,nz1 ',ny1,nz1
!   print *, 'fld1 ',fld1
!   print *, 'xlat1 ',xlat1
!   print *, 'xlev1 ',xlev1
!
! convert longitude to 0 - 360
   if (.not.  use_interp_1) then
      do i=1,nx2
         do j=1,ny2
            if(xlon2(i,j).lt.0.) then
               xlon2(i,j)=xlon2(i,j)+360.
            endif
         enddo
      enddo
   endif
!
! invert the pressure grid and data
   if (use_interp_1) then
      do k=1,nz1
         prs_tmp1(nz1-k+1)=xlev1(k)
         do j=1,ny1
            fld_tmp1(j,nz1-k+1)=fld1(j,k)
         enddo
      enddo
      xlev1(1:nz1)=prs_tmp1(1:nz1)*100.
      fldd1(1,1:ny1,1:nz1)=fld_tmp1(1:ny1,1:nz1)
      fldd1(2,1:ny1,1:nz1)=fld_tmp1(1:ny1,1:nz1)
!
! interpolate data1 to (lat,lev) point
      do j=1,ny1
         xlon1(1,j)=lon-del_lon1
         xlon1(2,j)=lon+del_lon1
         if(lon.lt.0.) then
            xlon1(1,j)=lon+360.-del_lon1
            xlon1(2,j)=lon+360.+del_lon1
         endif
         do i=1,nx1
            xlat1(i,j)=xlat_tmp1(j)
         enddo
      enddo
!      print *, 'IN UPVAL SUB: lon,lat,lev ',lon,lat,lev
!      print *, 'IN UPVAL SUB: xlon,xlat,xlev ',xlon1(1,48),xlat1(1,48)
!      do j=1,nz1
!        print *, 'IN UPVAL SUB: fldd1 ',j,xlev1(j),fldd1(1,48,j)
!      enddo
      call apm_interpolate(obs_val,del_prs,lon,lat,lev,xlon1,xlat1,xlev1, &
      fldd1,nx1,ny1,nz1,istatus)
!      print *, 'IN UPVAL SUB: obs_val,del_prs ',obs_val,del_prs
   else
      do k=1,nz2
         prs_tmp2(nz2-k+1)=xlev2(k)
         do i=1,nx2
            do j=1,ny2
               fld_tmp2(i,j,nz2-k+1)=o3_col_dens(i,j,k,im2)
            enddo
         enddo
      enddo
      xlev2(1:nz2)=prs_tmp2(1:nz2)
      o3_col_dens(1:nx2,1:ny2,1:nz2,im2)=fld_tmp2(1:nx2,1:ny2,1:nz2)
!
! interpolate data2 to (lat,lon,lev) point
      call apm_interpolate(obs_val,del_prs,lon,lat,lev,xlon2,xlat2,xlev2, &
      o3_col_dens(1,1,1,im2),nx2,ny2,nz2,istatus)
   endif
end subroutine wrf_dart_ubval_interp
!
subroutine apm_get_exo_coldens(fid,fldname,dataf,nx,ny,nz,nm)
   use netcdf
   implicit none
   integer,parameter                      :: maxdim=4
   integer                                :: nx,ny,nz,nm
   integer                                :: i,rc,v_ndim,natts,fid
   integer                                :: v_id,typ
   integer,dimension(maxdim)              :: v_dimid,v_dim,one
   character(len=*)                       :: fldname
   character(len=180)                     :: vnam
   real,dimension(nx,ny,nz,nm)            :: dataf
!
! get variables identifiers
   rc = nf90_inq_varid(fid,trim(fldname),v_id)
   if(rc.ne.0) then
      print *, 'APM: nf_inq_varid error'
      call abort
   endif
!
! get dimension identifiers
   v_dimid=0
   rc = nf90_inquire_variable(fid,v_id,vnam,typ,v_ndim,v_dimid,natts)
   if(rc.ne.0) then
      print *, 'APM: nc_inq_var error'
      call abort
   endif
   if(maxdim.lt.v_ndim) then
      print *, 'ERROR: maxdim is too small ',maxdim,v_ndim
      call abort
   endif 
!
! get dimensions
   v_dim(:)=1
   do i=1,v_ndim
      rc = nf90_inquire_dimension(fid,v_dimid(i),len=v_dim(i))
      if(rc.ne.0) then
         print *, 'APM: nf_inq_dimlen error'
         call abort
      endif
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
   else if(nm.ne.v_dim(4)) then             
      print *, 'ERROR: nm dimension conflict ',nm,v_dim(4)
      call abort
   endif
!
! get data
   one(:)=1
   rc = nf90_get_var(fid,v_id,dataf,one,v_dim)
end subroutine apm_get_exo_coldens
!
subroutine apm_get_upvals(fid,species,imn,dataf,lats,levs)
   use netcdf
   implicit none
   integer,parameter                                 :: maxdim=4
   integer,parameter                                 :: ny1=96,nz1=38,nm1=12,nmspc1=8,nmchr1=20
   integer                                           :: i,j,idx,fid,rc,typ,natts,imn
   integer                                           :: vid1,vid2,vid3,vid4,vid5
   integer                                           :: vndim1,vndim2,vndim3,vndim4,vndim5
   integer,dimension(maxdim)                         :: vdimid1,vdimid2,vdimid3,vdimid4,vdimid5
   integer,dimension(maxdim)                         :: one,vdim1,vdim2,vdim3,vdim4,vdim5
   integer,dimension(nm1)                            :: mths
   character(len=20)                                 :: species
   character(len=nmchr1),dimension(nmchr1,nmspc1)    :: spcs    
   character(len=180)                                :: vnam1,vnam2,vnam3,vnam4,vnam5
   real,dimension(ny1)                               :: lats
   real,dimension(nz1)                               :: levs
   real,dimension(ny1,nmspc1,nm1,nz1)                :: vmrs
   real,dimension(ny1,nz1)                           :: dataf
!
! get variables identifiers
   rc = nf90_inq_varid(fid,'lat',vid1)
   if(rc.ne.0) then
      print *, 'APM: nf_inq_varid error_1'
      call abort
   endif
   rc = nf90_inq_varid(fid,'lev',vid2)
   if(rc.ne.0) then
      print *, 'APM: nf_inq_varid error_2'
      call abort
   endif
   rc = nf90_inq_varid(fid,'month',vid3)
   if(rc.ne.0) then
      print *, 'APM: nf_inq_varid error_3'
      call abort
   endif
   rc = nf90_inq_varid(fid,'specname',vid4)
   if(rc.ne.0) then
      print *, 'APM: nf_inq_varid error_4'
      call abort
   endif
   rc = nf90_inq_varid(fid,'vmr',vid5)
   if(rc.ne.0) then
      print *, 'APM: nf_inq_varid error_5'
      call abort
   endif
!
! get dimension identifiers
   vdimid1=0
   rc = nf90_inquire_variable(fid,vid1,vnam1,typ,vndim1,vdimid1,natts)
   if(rc.ne.0) then
      print *, 'APM: nc_inq_var error_1'
      call abort
   endif
   vdimid2=0
   rc = nf90_inquire_variable(fid,vid2,vnam2,typ,vndim2,vdimid2,natts)
   if(rc.ne.0) then
      print *, 'APM: nc_inq_var error_2'
      call abort
   endif
   vdimid3=0
   rc = nf90_inquire_variable(fid,vid3,vnam3,typ,vndim3,vdimid3,natts)
   if(rc.ne.0) then
      print *, 'APM: nc_inq_var error_3'
      call abort
   endif
   vdimid4=0
   rc = nf90_inquire_variable(fid,vid4,vnam4,typ,vndim4,vdimid4,natts)
   if(rc.ne.0) then
      print *, 'APM: nc_inq_var error_4'
      call abort
   endif
   vdimid5=0
   rc = nf90_inquire_variable(fid,vid5,vnam5,typ,vndim5,vdimid5,natts)
   if(rc.ne.0) then
      print *, 'APM: nc_inq_var error_5'
      call abort
   endif
!
! test the number of dimensions
   if(1.lt.vndim1) then
      print *, 'ERROR: maxdim is too small 1 ',1,vndim1
      call abort
   endif 
   if(1.lt.vndim2) then
      print *, 'ERROR: maxdim is too small 2 ',1,vndim2
      call abort
   endif 
   if(1.lt.vndim3) then
      print *, 'ERROR: maxdim is too small 3 ',1,vndim3
      call abort
   endif 
   if(2.lt.vndim4) then
      print *, 'ERROR: maxdim is too small 4',1,vndim4
      call abort
   endif 
   if(4.lt.vndim5) then
      print *, 'ERROR: maxdim is too small 5',1,vndim5
      call abort
   endif 
!
! get dimensions
   vdim1(:)=1
   do i=1,vndim1
      rc = nf90_inquire_dimension(fid,vdimid1(i),len=vdim1(i))
      if(rc.ne.0) then
         print *, 'APM: nf_inq_dimlen error_1'
         call abort
      endif
   enddo
   vdim2(:)=1
   do i=1,vndim2
      rc = nf90_inquire_dimension(fid,vdimid2(i),len=vdim2(i))
      if(rc.ne.0) then
         print *, 'APM: nf_inq_dimlen error_2'
         call abort
      endif
   enddo
   vdim3(:)=1
   do i=1,vndim3
      rc = nf90_inquire_dimension(fid,vdimid3(i),len=vdim3(i))
      if(rc.ne.0) then
         print *, 'APM: nf_inq_dimlen error_3'
         call abort
      endif
   enddo
   vdim4(:)=1
   do i=1,vndim4
      rc = nf90_inquire_dimension(fid,vdimid4(i),len=vdim4(i))
      if(rc.ne.0) then
         print *, 'APM: nf_inq_dimlen error_4'
         call abort
      endif
   enddo
   vdim5(:)=1
   do i=1,vndim5
      rc = nf90_inquire_dimension(fid,vdimid5(i),len=vdim5(i))
      if(rc.ne.0) then
         print *, 'APM: nf_inq_dimlen error_5'
         call abort
      endif
   enddo
!
! check dimensions
   if(ny1.ne.vdim1(1)) then
      print *, 'ERROR: ny1 dimension conflict 1 ',ny1,vdim1(1)
      call abort
   else if(nz1.ne.vdim2(1)) then             
      print *, 'ERROR: nz1 dimension conflict 2 ',nz1,vdim2(1)
      call abort
   else if(nm1.ne.vdim3(1)) then             
      print *, 'ERROR: nm1 dimension conflict 3 ',nm1,vdim3(1)
      call abort
   endif
   if(nmchr1.ne.vdim4(1)) then             
      print *, 'ERROR: nmchr1 dimension conflict 4 ',nmchr1,vdim4(1)
      call abort
   else if(nmspc1.ne.vdim4(2)) then             
      print *, 'ERROR: nmspc1 dimension conflict 4 ',nmspc1,vdim4(2)
      call abort
   endif
   if(ny1.ne.vdim5(1)) then
      print *, 'ERROR: ny1 dimension conflict 5 ',ny1,vdim5(1)
      call abort
   else if(nmspc1.ne.vdim5(2)) then             
      print *, 'ERROR: nmspc1 dimension conflict 5 ',nmspc1,vdim5(2)
      call abort
   else if(nm1.ne.vdim5(3)) then             
      print *, 'ERROR: nm1 dimension conflict 5 ',nm1,vdim5(3)
      call abort
   else if(nz1.ne.vdim5(4)) then             
      print *, 'ERROR: nz1 dimension conflict 5 ',nz1,vdim5(4)
      call abort
   endif
!
! get data
   one(:)=1
   rc = nf90_get_var(fid,vid1,lats,one,vdim1)
   if(rc.ne.0) then
      print *, 'APM: get_var error_1'
      call abort
   endif
!   print *, 'lats ',lats
   one(:)=1
   rc = nf90_get_var(fid,vid2,levs,one,vdim2)
   if(rc.ne.0) then
      print *, 'APM: get_var error_2'
      call abort
   endif
!   print *, 'levs ',levs
   one(:)=1
   rc = nf90_get_var(fid,vid3,mths,one,vdim3)
   if(rc.ne.0) then
      print *, 'APM: get_var error_3'
      call abort
   endif
!   print *, 'mths ',mths
   one(:)=1
   rc = nf90_get_var(fid,vid4,spcs,one,vdim4)
   if(rc.ne.0) then
      print *, 'APM: get_var error_4'
      call abort
   endif
!   print *, 'spcs ',spcs
   one(:)=1
   rc = nf90_get_var(fid,vid5,vmrs,one,vdim5)
   if(rc.ne.0) then
      print *, 'APM: get_var error_5'
      call abort
   endif
!   print *, 'vmrs ',vmrs
!
! locate requested field
  do i=1,nmspc1
     if(trim(species).eq.trim(spcs(i,1))) then
        idx=i
        exit
     endif
  enddo
   do i=1,ny1
      do j=1,nz1
         dataf(i,j)=vmrs(i,idx,imn,j)
      enddo
   enddo
end subroutine apm_get_upvals
!
subroutine apm_interpolate(obs_val,del_prs,lon,lat,lev,xlon,xlat,xlev,dataf,nx,ny,nz,istatus)
!
! longitude and latitude must be in degrees
! pressure grid must be in hPa and go from bottom to top
!
   use        types_mod, only     : r8
   implicit none
   integer                                :: nx,ny,nz,nzm,istatus
   integer                                :: i,j,k,im,ip,jm,jp,quad
   integer                                :: k_lw,k_up,i_min,j_min 
   real(r8)                               :: obs_val,del_prs
   real(r8)                               :: lon,lat,lev
   real                                   :: l_lon,l_lat,l_lev
   real                                   :: fld_lw,fld_up
   real                                   :: xlnp_lw,xlnp_up,xlnp_pt
   real                                   :: dz_lw,dz_up
   real                                   :: mop_x,mop_y
   real                                   :: re,pi,rad2deg
   real                                   :: rad,rad_crit,rad_min,mod_x,mod_y
   real                                   :: dx_dis,dy_dis
   real                                   :: w_q1,w_q2,w_q3,w_q4,wt
   real,dimension(nz)                     :: xlev
   real,dimension(nx,ny)                  :: xlon,xlat
   real,dimension(nx,ny,nz)               :: dataf
!
! set constants
   pi=4.*atan(1.)
   rad2deg=360./(2.*pi)
   re=6371000.
   rad_crit=200000.
   quad=0
!
! find the closest point            
   rad_min=1.e10
   l_lon=lon
   l_lat=lat
   l_lev=lev
   if(l_lon.lt.0.) l_lon=l_lon+360.
!   print *, 'lon,lat,lev ',l_lon,l_lat,l_lev
!
   do i=1,nx
      do j=1,ny
         mod_x=(xlon(i,j))/rad2deg
         if(xlon(i,j).lt.0.) mod_x=(360.+xlon(i,j))/rad2deg
         mod_y=xlat(i,j)/rad2deg
         mop_x=l_lon/rad2deg
         mop_y=l_lat/rad2deg
         dx_dis=abs(mop_x-mod_x)*cos((mop_y+mod_y)/2.)*re
         dy_dis=abs(mop_y-mod_y)*re
         rad=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)  
         rad_min=min(rad_min,rad)
         if(rad.eq.rad_min) then
            i_min=i
            j_min=j
         endif 
      enddo
   enddo
   if(rad_min.gt.rad_crit) then            
      print *, 'APM: ERROR in intrp - min dist exceeds threshold ',rad_min, rad_crit
      print *, 'grid ',i_min,j_min,xlon(i_min,j_min),xlat(i_min,j_min)
      print *, 'point ',l_lon,l_lat
      istatus=2
      return
!      call abort
   endif
!
! do interpolation
   im=i_min-1
   if(im.eq.0) im=1
   ip=i_min+1
   if(ip.eq.nx+1) ip=nx
   jm=j_min-1
   if(jm.eq.0) jm=1
   jp=j_min+1
   if(jp.eq.ny+1) jp=ny
!
! find quadrant and interpolation weights
   quad=0
   mod_x=xlon(i_min,j_min)
   if(xlon(i_min,j_min).lt.0.) mod_x=xlon(i_min,j_min)+360.
   mod_y=xlat(i_min,j_min)
   if(mod_x.ge.l_lon.and.mod_y.ge.l_lat) quad=1 
   if(mod_x.le.l_lon.and.mod_y.ge.l_lat) quad=2 
   if(mod_x.le.l_lon.and.mod_y.le.l_lat) quad=3 
   if(mod_x.ge.l_lon.and.mod_y.le.l_lat) quad=4
   if(quad.eq.0) then
      print *, 'APM: ERROR IN INTERPOLATE quad = 0 '
      call abort
   endif
!
! Quad 1
   if (quad.eq.1) then
      mod_x=xlon(i_min,j_min)
      if(xlon(i_min,j_min).lt.0.) mod_x=360.+xlon(i_min,j_min) 
      dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(i_min,j_min))/rad2deg/2.)*re
      dy_dis=abs(l_lat-xlat(i_min,j_min))/rad2deg*re
      w_q1=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mod_x=xlon(im,j_min)
      if(xlon(im,j_min).lt.0.) mod_x=360.+xlon(im,j_min) 
      dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(im,j_min))/rad2deg/2.)*re
      dy_dis=abs(l_lat-xlat(im,j_min))/rad2deg*re
      w_q2=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mod_x=xlon(im,jm)
      if(xlon(im,jm).lt.0.) mod_x=360.+xlon(im,jm) 
      dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(im,jm))/rad2deg/2.)*re
      dy_dis=abs(l_lat-xlat(im,jm))/rad2deg*re
      w_q3=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mod_x=xlon(i_min,jm)
      if(xlon(i_min,jm).lt.0.) mod_x=360.+xlon(i_min,jm) 
      dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(i_min,jm))/rad2deg/2.)*re
      dy_dis=abs(l_lat-xlat(i_min,jm))/rad2deg*re
      w_q4=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
! Quad 2
   else if (quad.eq.2) then
      mod_x=xlon(ip,j_min)
      if(xlon(ip,j_min).lt.0.) mod_x=360.+xlon(ip,j_min) 
      dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(ip,j_min))/rad2deg/2.)*re
      dy_dis=abs(l_lat-xlat(ip,j_min))/rad2deg*re
      w_q1=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mod_x=xlon(i_min,j_min)
      if(xlon(i_min,j_min).lt.0.) mod_x=360.+xlon(i_min,j_min) 
      dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(i_min,j_min))/rad2deg/2.)*re
      dy_dis=abs(l_lat-xlat(i_min,j_min))/rad2deg*re
      w_q2=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mod_x=xlon(i_min,jm)
      if(xlon(i_min,jm).lt.0.) mod_x=360.+xlon(i_min,jm) 
      dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(i_min,jm))/rad2deg/2.)*re
      dy_dis=abs(l_lat-xlat(i_min,jm))/rad2deg*re
      w_q3=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mod_x=xlon(ip,jm)
      if(xlon(ip,jm).lt.0.) mod_x=360.+xlon(ip,jm) 
      dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(ip,jm))/rad2deg/2.)*re
      dy_dis=abs(l_lat-xlat(ip,jm))/rad2deg*re
      w_q4=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
! Quad 3
   else if (quad.eq.3) then
      mod_x=xlon(ip,jp)
      if(xlon(ip,jp).lt.0.) mod_x=360.+xlon(ip,jp) 
      dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(ip,jp))/rad2deg/2.)*re
      dy_dis=abs(l_lat-xlat(ip,jp))/rad2deg*re
      w_q1=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mod_x=xlon(i_min,jp)
      if(xlon(i_min,jp).lt.0.) mod_x=360.+xlon(i_min,jp) 
      dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(i_min,jp))/rad2deg/2.)*re
      dy_dis=abs(l_lat-xlat(i_min,jp))/rad2deg*re
      w_q2=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mod_x=xlon(i_min,j_min)
      if(xlon(i_min,j_min).lt.0.) mod_x=360.+xlon(i_min,j_min) 
      dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(i_min,j_min))/rad2deg/2.)*re
      dy_dis=abs(l_lat-xlat(i_min,j_min))/rad2deg*re
      w_q3=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mod_x=xlon(ip,j_min)
      if(xlon(ip,j_min).lt.0.) mod_x=360.+xlon(ip,j_min) 
      dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(ip,j_min))/rad2deg/2.)*re
      dy_dis=abs(l_lat-xlat(ip,j_min))/rad2deg*re
      w_q4=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
! Quad 4
   else if (quad.eq.4) then
      mod_x=xlon(i_min,jp)
      if(xlon(i_min,jp).lt.0.) mod_x=360.+xlon(i_min,jp) 
      dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(i_min,jp))/rad2deg/2.)*re
      dy_dis=abs(l_lat-xlat(i_min,jp))/rad2deg*re
      w_q1=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mod_x=xlon(im,jp)
      if(xlon(im,jp).lt.0.) mod_x=360.+xlon(im,jp) 
      dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(im,jp))/rad2deg/2.)*re
      dy_dis=abs(l_lat-xlat(im,jp))/rad2deg*re
      w_q2=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mod_x=xlon(im,jm)
      if(xlon(im,jm).lt.0.) mod_x=360.+xlon(im,jm) 
      dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(im,jm))/rad2deg/2.)*re
      dy_dis=abs(l_lat-xlat(im,jm))/rad2deg*re
      w_q3=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      mod_x=xlon(i_min,j_min)
      if(xlon(i_min,j_min).lt.0.) mod_x=360.+xlon(i_min,j_min) 
      dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(i_min,j_min))/rad2deg/2.)*re
      dy_dis=abs(l_lat-xlat(i_min,j_min))/rad2deg*re
      w_q4=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
   endif
   if(l_lon.ne.xlon(i_min,j_min).or.l_lat.ne.xlat(i_min,j_min)) then
      wt=1./w_q1+1./w_q2+1./w_q3+1./w_q4
   endif
!
! find vertical indexes
   nzm=nz-1
   k_lw=-1
   k_up=-1
   do k=1,nzm
      if(k.eq.1 .and. l_lev.gt.xlev(k)) then
         k_lw=k
         k_up=k
         exit
      endif
      if(l_lev.le.xlev(k) .and. l_lev.gt.xlev(k+1)) then
         k_lw=k
         k_up=k+1
         exit
      endif
      if(k.eq.nzm .and. l_lev.ge.xlev(k+1)) then
         k_lw=k+1
         k_up=k+1
         exit
      endif
   enddo
   if(k_lw.le.0 .or. k_up.le.0) then
      print *, 'APM: ERROR IN K_LW OR K_UP ',k_lw,k_up
      call abort
   endif
!
! horizontal interpolation             
   fld_lw=0.
   fld_up=0.   
   if(l_lon.eq.xlon(i_min,j_min).and.l_lat.eq.xlat(i_min,j_min)) then
      fld_lw=dataf(i_min,j_min,k_lw)
      fld_up=dataf(i_min,j_min,k_up)
   else if(quad.eq.1) then
      fld_lw=(1./w_q1*dataf(i_min,j_min,k_lw)+1./w_q2*dataf(im,j_min,k_lw)+ &
      1./w_q3*dataf(im,jm,k_lw)+1./w_q4*dataf(i_min,jm,k_lw))/wt
      fld_up=(1./w_q1*dataf(i_min,j_min,k_up)+1./w_q2*dataf(im,j_min,k_up)+ &
      1./w_q3*dataf(im,jm,k_up)+1./w_q4*dataf(i_min,jm,k_up))/wt
   else if(quad.eq.2) then
      fld_lw=(1./w_q1*dataf(ip,j_min,k_lw)+1./w_q2*dataf(i_min,j_min,k_lw)+ &
      1./w_q3*dataf(i_min,jm,k_lw)+1./w_q4*dataf(ip,jm,k_lw))/wt
      fld_up=(1./w_q1*dataf(ip,j_min,k_up)+1./w_q2*dataf(i_min,j_min,k_up)+ &
      1./w_q3*dataf(i_min,jm,k_up)+1./w_q4*dataf(ip,jm,k_up))/wt
   else if(quad.eq.3) then
      fld_lw=(1./w_q1*dataf(ip,jp,k_lw)+1./w_q2*dataf(i_min,jp,k_lw)+ &
      1./w_q3*dataf(i_min,j_min,k_lw)+1./w_q4*dataf(ip,j_min,k_lw))/wt
      fld_up=(1./w_q1*dataf(ip,jp,k_up)+1./w_q2*dataf(i_min,jp,k_up)+ &
      1./w_q3*dataf(i_min,j_min,k_up)+1./w_q4*dataf(ip,j_min,k_up))/wt
   else if(quad.eq.4) then
      fld_lw=(1./w_q1*dataf(i_min,jp,k_lw)+1./w_q2*dataf(im,jp,k_lw)+ &
      1./w_q3*dataf(im,j_min,k_lw)+1./w_q4*dataf(i_min,j_min,k_lw))/wt
      fld_up=(1./w_q1*dataf(i_min,jp,k_up)+1./w_q2*dataf(im,jp,k_up)+ &
      1./w_q3*dataf(im,j_min,k_up)+1./w_q4*dataf(i_min,j_min,k_up))/wt
   endif 
!   print *,'fld_lw ',fld_lw
!   print *,'fld_up ',fld_up
!
! vertical interpolation
!   print *,'p_lw,p_up,p ',xlev(k_lw),xlev(k_up),l_lev

   xlnp_lw=log(xlev(k_lw))
   xlnp_up=log(xlev(k_up))
   xlnp_pt=log(l_lev)
   dz_lw=xlnp_lw-xlnp_pt
   dz_up=xlnp_pt-xlnp_up
   if(dz_lw.eq.0.) then
      obs_val=fld_lw
   else if(dz_up.eq.0.) then
      obs_val=fld_up
   else if(dz_lw.ne.0. .and. dz_up.ne.0.) then
      obs_val=(1./dz_lw*fld_lw+1./dz_up*fld_up)/(1./dz_lw+1./dz_up)
   endif
   del_prs=xlev(k_lw)-xlev(k_up)
end subroutine apm_interpolate

end module apm_upper_bdy_mod
