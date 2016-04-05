! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program create_mopitt_obs_sequence

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!=============================================
! Created by AFA for MOPITT CO retrieval obs
! Based from create_obs_sequence.f90
!=============================================
!
use    utilities_mod, only : timestamp, 		&
                             register_module, 		&
                             open_file, 		&
                             close_file, 		&
                             initialize_utilities, 	&
                             open_file, 		&
                             close_file, 		&
                             find_namelist_in_file,  	&
                             check_namelist_read,    	&
                             error_handler, 		&
                             E_ERR,			& 
                             E_WARN,			& 
                             E_MSG, 			&
                             E_DBG

use obs_sequence_mod, only : obs_sequence_type, 	&
                             interactive_obs, 		&
                             write_obs_seq, 		&
                             interactive_obs_sequence,  &
                             static_init_obs_sequence,  &
                             init_obs_sequence,         &
                             init_obs,                  &
                             set_obs_values,            &
                             set_obs_def,               &
                             set_qc,                    &
                             set_qc_meta_data,          &
                             set_copy_meta_data,        &
                             insert_obs_in_seq,         &
                             obs_type
                    
use obs_def_mod, only      : set_obs_def_kind,          &
                             set_obs_def_location,      &
                             set_obs_def_time,          &
                             set_obs_def_key,           &
                             set_obs_def_error_variance,&
                             obs_def_type,              &
                             init_obs_def,              &
                             get_obs_kind

use obs_def_mopitt_mod, only :  set_obs_def_mopitt_co

use  assim_model_mod, only : static_init_assim_model

use location_mod, only  : location_type, 		&
                          set_location

use time_manager_mod, only : set_date, 			&
                             set_calendar_type, 	&
                             time_type, 		&
                             get_time

use obs_kind_mod, only  : KIND_CO, 		&
                          MOPITT_CO_RETRIEVAL, 		&
                          get_kind_from_menu

use random_seq_mod, only : random_seq_type, 		&
                           init_random_seq, 		&
                           random_uniform

use sort_mod, only       : index_sort


implicit none

! version controlled file description for error handling, do not edit                          
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"


! add variables AFA
type(obs_sequence_type) :: seq
type(obs_type)          :: obs
type(obs_type)          :: obs_old
type(obs_def_type)      :: obs_def
type(location_type)     :: obs_location
type(time_type)         :: obs_time
integer                 :: obs_kind
integer                 :: obs_key

character*4             :: chr_year
character*2             :: chr_month, chr_day,chr_hour
character*129           :: filen, filedir, fileapm
integer,parameter       :: fileid=88
integer,parameter       :: fileidapm=99
integer                 :: calendar_type, year, month, day, hour, day1,month1,year1
integer                 :: seconds, days
integer                 :: ios, i, k1, k2 , icopy ! kmr, kmc
!real                   :: pi=3.1415926535898

real                    :: bin_beg, bin_end
integer                 :: days_in_month(12) =(/ &
                            31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31  /)
!real                    :: period_1, period_2
!integer                 :: period

!============================================================
!mopitt variables
!============================================================

integer, parameter      :: mop_dim=10
integer, parameter      :: nlev_cam=10

integer         :: k,kk
real	        :: nlevels
real            :: sec
real            :: lon
real            :: lat
real            :: psurf
real            :: err
real            :: nprs(mop_dim) =(/ &
                                1000.,900.,800.,700.,600.,500.,400.,300.,200.,100. /)


!============================================================
!obs sequence extra variables
!============================================================
integer,parameter   :: max_num_obs=1000000
integer,parameter   :: num_copies=1
integer,parameter   :: num_qc=1
character*129   :: copy_meta_data
character*129   :: qc_meta_data='BASED ON MOPITT QUALITY STATEMENT'
character*129   :: file_name='mopitt_obs_seq'
character*129   :: filename

!============================================================
!dummy variables
!============================================================
real*8          :: longitude
real*8          :: latitude
real*8          :: level
integer         :: which_vert
integer         :: hour1
integer         :: minute
integer         :: second
integer         :: qc_count
real*8, dimension(num_qc)       :: co_qc
real*8, dimension(nlev_cam)     :: co_avgker
real*8, dimension(num_copies)   :: co_vmr
real*8                          :: co_error
real*8                          :: co_prior
real*8                          :: co_psurf

integer,parameter	:: nlon_qc=144
integer,parameter	:: nlat_qc=96
integer,parameter	:: nqc_obs=40

real*8 :: dlon_qc=2.50
real*8 :: dlat_qc=1.9149
real*8 :: eps_tol=1.e-5

integer :: index_qc, xg_count(nlon_qc,nlat_qc), xg(nlon_qc,nlat_qc,500)
integer :: qc_mopitt(1000000), qc_thinning(1000000)
integer :: index_20(1000), index , j,ik, ikk
integer :: lon_qc, lat_qc, io, iunit

integer :: mlev

integer :: line_count

type (random_seq_type) :: inc_ran_seq
real*8                 :: unif(1000)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! APM: variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
character*129 :: transform_typ
!
integer :: nlvls,lwrk,info,nlev      
integer :: nrows,nrows_c,nrows_d
!
real*8 :: dofs,qstatus,sum,cum
!
real*8, allocatable, dimension(:) :: mop_prs,x_r,x_p,var
real*8, allocatable, dimension(:) :: mopittlev,xcomp,xcomperr,xapr
real*8, allocatable, dimension(:,:) :: avg_k,cov_p,cov_r,cov_m,cov_a
real*8, allocatable, dimension(:,:) :: avgker
!
double precision, allocatable, dimension(:) :: wrk,SV_cvr,SV_cva,rs_avg_xp
double precision, allocatable, dimension(:) :: s_x_r,s_x_p,SV_s_avg_k
double precision, allocatable, dimension(:) :: rs_x_r,rs_x_p,w,err2_rs_r,err2_rs_a
!
double precision, allocatable, dimension(:,:) :: Z,ZT,U_cvr,V_cvr,U_cva,V_cva
double precision, allocatable, dimension(:,:) :: UT_cvr,VT_cvr,UT_cva,VT_cva
double precision, allocatable, dimension(:,:) :: SV,s_cvr,s_cva
double precision, allocatable, dimension(:,:) :: s_cov_r,s_cov_a,s_avg_k
double precision, allocatable, dimension(:,:) :: U_s_avg_k,UT_s_avg_k,VT_s_avg_k,rs_avg_k
double precision, allocatable, dimension(:,:) :: rs_cov_a,rs_cov_r

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
namelist /create_mopitt_obs_nml/filedir,filename,year,month,day,hour,bin_beg, bin_end

call find_namelist_in_file("input.nml", "create_mopitt_obs_nml", iunit)
read(iunit, nml = create_mopitt_obs_nml, iostat = io)
call check_namelist_read(iunit, io, "create_mopitt_obs_nml")

! Record the namelist values used for the run ...
call error_handler(E_MSG,'init_create_mopitt_obs','create_mopitt_obs_nml values are',' ',' ',' ')
write(     *     , nml=create_mopitt_obs_nml)

! Record the current time, date, etc. to the logfile
call initialize_utilities('create_obs_sequence')
call register_module(source,revision,revdate)


! Initialize the assim_model module, need this to get model
! state meta data for locations of identity observations
!call static_init_assim_model()

! Initialize the obs_sequence module
call static_init_obs_sequence()

! Initialize an obs_sequence structure
call init_obs_sequence(seq, num_copies, num_qc, max_num_obs)

! Initialize the obs variable
call init_obs(obs, num_copies, num_qc)

do icopy =1, num_copies
   if (icopy == 1) then
       copy_meta_data='MOPITT CO observation'
   else
       copy_meta_data='Truth'
   endif
   call set_copy_meta_data(seq, icopy, copy_meta_data)
enddo

call set_qc_meta_data(seq, 1, qc_meta_data)

qc_mopitt(:)=100
qc_thinning(:)=100


!-------------------------------------------------------
! Read MOPITT obs
!-------------------------------------------------------

! Set dates and initialize qc_count
  calendar_type=3                          !Gregorian
  call set_calendar_type(calendar_type)
  qc_count=0

! Perhaps make this a time loop for later runs
  write(chr_year,'(i4.4)') year
  write(chr_month,'(i2.2)') month
  write(chr_day,'(i2.2)') day
  write(chr_hour,'(i2.2)') hour

  if ( mod(year,4) == 0 ) then
       days_in_month(2) = days_in_month(2) + 1
  endif
  if ( mod(year,100) == 0 ) then
       days_in_month(2) = days_in_month(2) - 1
  endif
  if ( mod(year,400) == 0 ) then
       days_in_month(2) = days_in_month(2) + 1
  endif

! Open MOPITT binary file
  filen=chr_year//chr_month//chr_day//chr_hour//'.dat'
  write(6,*)'opening ',TRIM(filedir)//TRIM(filen)

! Open file to save raw MOPITT data corresponding to the obs_seq file
!  fileapm=chr_year//chr_month//chr_day//chr_hour//'.obsseq.dat'
!  open(fileidapm,file=TRIM(filedir)//TRIM(fileapm),                     &
!       form='formatted', status='unknown',  &
!       iostat=ios)

! Read MOPITT file 1
  index_qc=0
  line_count = 0
  open(fileid,file=TRIM(filedir)//TRIM(filen),                     &
       form='formatted', status='old',  &
       iostat=ios)

! Error Check
  if (ios /=0) then
      write(6,*) 'no mopitt file for the day ', day
      go to 999
  endif

! Read MOPITT
  read(fileid,*,iostat=ios) transform_typ, sec, lat, lon, nlevels, dofs 
!  print *, 'trans_typ, sec, lat, lon, nlevels, dofs ',trim(transform_typ),sec,lat,lon,nlevels,dofs
  nlvls=int(nlevels)

! Error Check
  if (ios /=0) then
      write(6,*) 'no data on file ', TRIM(filen)
      go to 999
  endif

!-------------------------------------------------------
! MAIN LOOP FOR MOPITT OBS
!-------------------------------------------------------
  do while(ios == 0)

       ! Read MOPITT variables
       allocate (mop_prs(nlvls))
       allocate (x_r(nlvls))
       allocate (x_p(nlvls))
       allocate (avg_k(nlvls,nlvls))
       allocate (cov_a(nlvls,nlvls))
       allocate (cov_r(nlvls,nlvls))
       allocate (cov_m(nlvls,nlvls))
       read(fileid,*) mop_prs
       read(fileid,*) x_r
       read(fileid,*) x_p
       read(fileid,*) avg_k
       read(fileid,*) cov_a
       read(fileid,*) cov_r
       read(fileid,*) cov_m
       deallocate (mop_prs,x_r,x_p,avg_k,cov_a,cov_r,cov_m)
!       
       index_qc = index_qc + 1
       qc_mopitt(index_qc)=0

       !-------------------------------------------------------
       ! Bin to nlat_qcxnlon_qc
       !-------------------------------------------------------

       ! find lon_qc, lat_qc
       lon_qc=nint((lon+180)/dlon_qc) + 1
       lat_qc=nint((lat+90)/dlat_qc) + 1
       if (lat>89.5) then
           lat_qc=nlat_qc
       elseif (lat<-89.5) then
           lat_qc=1
       endif

       !print *, 'lon lat ',lon_qc, lat_qc
       xg_count(lon_qc,lat_qc)=xg_count(lon_qc,lat_qc)+1
       xg(lon_qc,lat_qc,xg_count(lon_qc,lat_qc))=index_qc

       !read next data point
       read(fileid,*,iostat=ios) transform_typ, sec, lat, lon, nlevels, dofs 
       nlvls=int(nlevels)
       ! print *, nlevels, sec, nlevels, mopittlev, lat, lon
       !if (sec < 0.0) then 
       !    print *, line_count, nlev, sec 
       !    stop
       !endif
       !stop

  enddo !ios

9999   continue

close(fileid)

! Now do the thinning
  call init_random_seq(inc_ran_seq)
  do i=1,nlon_qc
     do j=1,nlat_qc
        if (xg_count(i,j)>nqc_obs) then

            ! draw nqc_obs
              do ik=1,xg_count(i,j)
                  unif(ik)=random_uniform(inc_ran_seq)
              enddo

              call index_sort(unif,index_20,xg_count(i,j))

              do ik=1,nqc_obs
                    index=xg(i,j,index_20(ik))
                    qc_thinning(index)=0
              enddo

        else
              do k=1,xg_count(i,j)
                   index=xg(i,j,k)
                   qc_thinning(index)=0
              enddo

        endif !xg_count
     enddo !j
  enddo !i

!===================================================================================

! Read MOPITT file AGAIN

  index_qc=0
! NOTE NOTE NOTE Check if it should be BIG_ENDIAN
  open(fileid,file=TRIM(filedir)//TRIM(filen),                     &
       form='formatted', status='old',   &
       iostat=ios)

! Error Check
  if (ios /=0) then
      write(6,*) 'no mopitt file for the day ', day
      go to 999
  endif

! Read MOPITT
    read(fileid,*,iostat=ios) transform_typ, sec, lat, lon, nlevels, dofs 
    nlvls=int(nlevels)
! Error Check
  if (ios /=0) then
      write(6,*) 'no data on file ', TRIM(filen)
      go to 999
  endif

!-------------------------------------------------------
! MAIN LOOP FOR MOPITT OBS
!-------------------------------------------------------
  do while(ios == 0)

     ! Read MOPITT variables
       allocate (mop_prs(nlvls))
       allocate (x_r(nlvls))
       allocate (x_p(nlvls))
       allocate (avg_k(nlvls,nlvls))
       allocate (cov_a(nlvls,nlvls))
       allocate (cov_r(nlvls,nlvls))
       allocate (cov_m(nlvls,nlvls))
       read(fileid,*) mop_prs
       read(fileid,*) x_r
       read(fileid,*) x_p
       read(fileid,*) avg_k
       read(fileid,*) cov_a
       read(fileid,*) cov_r
       read(fileid,*) cov_m
!
!       print *, 'mop_prs ',mop_prs        
!       print *, 'x_r ',x_r        
!       print *, 'x_p ',x_p
!       print *, 'avg_k ',avg_k
!       print *, 'cov_a ',cov_a
!       print *, 'cov_r ',cov_r
!       print *, 'cov_m ',cov_m
       index_qc=index_qc+1
       if ( (qc_mopitt(index_qc)==0).and.(qc_thinning(index_qc)==0) ) then
             co_qc(1)=0
       else
             co_qc(1)=100
       endif
       !print *, qc_mopitt(index_qc),qc_thinning(index_qc)

       if ( co_qc(1) == 0 )  then
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Scale the foreward operator with the inverse of cov_r  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate SVD of cov_r (Z=U_xxx * SV_xxx * VT_xxx)
          lwrk=max(4*nlvls,5*nlvls)
          allocate (Z(nlvls,nlvls),SV(nlvls,nlvls),U_cvr(nlvls,nlvls),V_cvr(nlvls,nlvls))
          allocate (UT_cvr(nlvls,nlvls),VT_cvr(nlvls,nlvls),wrk(lwrk))
          allocate (SV_cvr(nlvls),s_cvr(nlvls,nlvls))
          Z(:,:)=cov_r(:,:)
          call dgesvd('A','A',nlvls,nlvls,Z,nlvls,SV_cvr,U_cvr,nlvls,VT_cvr,nlvls,wrk,lwrk,info)
          deallocate(wrk)
!
! Calculate inverse square root of cov_r
          call mat_transpose(U_cvr,UT_cvr,nlvls,nlvls)
          call mat_transpose(VT_cvr,V_cvr,nlvls,nlvls)
          call vec_to_mat(SV_cvr,SV,nlvls)
          call diag_inv_sqrt(SV,nlvls)
          call mat_tri_prd(V_cvr,SV,UT_cvr,s_cvr,nlvls,nlvls,nlvls,nlvls,nlvls,nlvls)
!
! Calculate SVD of cov_a (Z=U_xxx * S_xxx * VT_xxx)
          lwrk=max(4*nlvls,5*nlvls)
          allocate (SV_cva(nlvls),U_cva(nlvls,nlvls),V_cva(nlvls,nlvls))
          allocate (UT_cva(nlvls,nlvls),VT_cva(nlvls,nlvls),wrk(lwrk))
          allocate (s_cva(nlvls,nlvls))
          Z=cov_a(:,:)
          call dgesvd('A','A',nlvls,nlvls,Z,nlvls,SV_cva,U_cva,nlvls,VT_cva,nlvls,wrk,lwrk,info)
          deallocate(wrk)
!
! Calculate square root of cov_a
          call mat_transpose(U_cva,UT_cva,nlvls,nlvls)
          call mat_transpose(VT_cva,V_cva,nlvls,nlvls)
          call vec_to_mat(SV_cva,SV,nlvls)
          call diag_sqrt(SV,nlvls)
          call mat_tri_prd(V_cva,SV,UT_cva,s_cva,nlvls,nlvls,nlvls,nlvls,nlvls,nlvls)
!
          deallocate(SV_cvr,SV_cva,U_cvr,U_cva,V_cvr,V_cva)
          deallocate(Z,SV,UT_cvr,UT_cva,VT_cvr,VT_cva)
!
! Precondition terms in the forward operator
          allocate(s_avg_k(nlvls,nlvls),s_cov_r(nlvls,nlvls),s_cov_a(nlvls,nlvls))
          allocate(s_x_r(nlvls),s_x_p(nlvls))
!       
          call mat_prd(s_cvr,avg_k,s_avg_k,nlvls,nlvls,nlvls,nlvls)
! s_cov_r should be identity, it is.
          call mat_tri_prd(s_cvr,cov_r,s_cvr,s_cov_r,nlvls,nlvls,nlvls,nlvls,nlvls,nlvls)
          call mat_tri_prd(s_cvr,cov_a,s_cvr,s_cov_a,nlvls,nlvls,nlvls,nlvls,nlvls,nlvls)
          call lh_mat_vec_prd(s_cvr,x_r,s_x_r,nlvls)
          call lh_mat_vec_prd(s_cvr,x_p,s_x_p,nlvls)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Rotate the forward operator based on the SVD of 
! the averaging kernal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate SVD of the scaled averaging kernal
          lwrk=max(4*nlvls,5*nlvls)
          allocate (Z(nlvls,nlvls),ZT(nlvls,nlvls),SV_s_avg_k(nlvls))
          allocate (U_s_avg_k(nlvls,nlvls),VT_s_avg_k(nlvls,nlvls),wrk(lwrk))
          allocate (UT_s_avg_k(nlvls,nlvls),rs_avg_k(nlvls,nlvls))
          Z(:,:)=s_avg_k(:,:)
          call dgesvd('A','A',nlvls,nlvls,Z,nlvls,SV_s_avg_k,U_s_avg_k,nlvls,VT_s_avg_k,nlvls,wrk,lwrk,info)
          deallocate(wrk)
!
! Rotate s_avg_k with transpose U_s_avg_k  
          call mat_transpose(U_s_avg_k,UT_s_avg_k,nlvls,nlvls)
          call mat_prd(UT_s_avg_k,s_avg_k,rs_avg_k,nlvls,nlvls,nlvls,nlvls)
!
! Rotate remaining terms in the forward operator
          allocate(rs_cov_a(nlvls,nlvls),rs_cov_r(nlvls,nlvls),rs_x_r(nlvls),rs_x_p(nlvls))    
          call mat_tri_prd(UT_s_avg_k,s_cov_a,U_s_avg_k,rs_cov_a,nlvls,nlvls,nlvls,nlvls,nlvls,nlvls)
          call mat_tri_prd(UT_s_avg_k,s_cov_r,U_s_avg_k,rs_cov_r,nlvls,nlvls,nlvls,nlvls,nlvls,nlvls)
          call lh_mat_vec_prd(UT_s_avg_k,s_x_r,rs_x_r,nlvls)
          call lh_mat_vec_prd(UT_s_avg_k,s_x_p,rs_x_p,nlvls)
!
! Calculate rotated/scaled prior term (rs_avg_ker*x_p)
          allocate (rs_avg_xp(nlvls))
          call lh_mat_vec_prd(rs_avg_k,x_p,rs_avg_xp,nlvls)
!
! Get new errors (check if err2_rs_r < 0 the qstatus=1)
          allocate(err2_rs_r(nlvls),err2_rs_a(nlvls))
          call diag_vec(rs_cov_r,err2_rs_r,nlvls)
          call diag_vec(rs_cov_a,err2_rs_a,nlvls)
          qstatus=0.0
          do i=1,nlvls
             if (rs_cov_r(i,i) .lt. 0.) then
!
! APM need to add skip for when qstatus = 1.0
                qstatus = 1.0
             else
                err2_rs_r(i)=sqrt(rs_cov_r(i,i))
                err2_rs_a(i)=sqrt(rs_cov_a(i,i))
             endif
          enddo
!
! APM count number of non-zero singular values
!          allocate(var(nlvls)) 
!          nrows_d=int(dofs)
!          if (dofs-nrows.gt.0.) nrows=nrows+1
!          sum=0.
!          do i=1,nlvls
!             sum=sum+rs_avg_k(i,i)
!          enddo
!          var(:)=0.
!          do i=1,nlvls
!             var(i)=rs_avg_k(i,i)/sum
!          enddo
!          cum=0
!          do i=1,nlvls
!             cum=cum+var(i)
!             if (cum .ge. 0.95) then
!                nrows_c=i
!                exit
!             endif
!          enddo
!          nrows=nrows_c
!          if (nrows_d .gt. nrows_c) then
!             print *, 'WARNING: nrows_c < nrows_d ',nrows_c,nrows_d
!          endif
!          deallocate(var)
!
           nrows=0
           do i=1,nlvls
              if(abs(SV_s_avg_k(i)).gt.eps_tol) then
                 nrows=nrows+1
              endif
           enddo 
!
! APM results match Ave's results to this point for the leading components  
! APM make assignments to Ave's scaled variables
          if (qstatus .eq. 0) then
             nlev=nrows
             allocate(xcomp(nlev),xcomperr(nlev),xapr(nlev),avgker(nlev,nlvls))
             allocate(mopittlev(nlvls))
             mopittlev(:)=mop_prs(:)
             do i=1,nlev
                xcomp(i)=rs_x_r(i)
                xcomperr(i)=err2_rs_r(i)
                xapr(i)=rs_x_p(i)-rs_avg_xp(i)
                do j=1,nlvls
                   avgker(i,j)=rs_avg_k(i,j)
                enddo
             enddo
!
! APM clean up allocatables from my code
             deallocate(s_cvr,s_cva,s_cov_r,s_cov_a,s_avg_k,s_x_r,s_x_p)
             deallocate(rs_cov_r,rs_cov_a,rs_avg_k,rs_x_r,rs_x_p)
             deallocate(rs_avg_xp,err2_rs_r,err2_rs_a)
             deallocate(Z,ZT,SV_s_avg_k,U_s_avg_k)
             deallocate(UT_s_avg_k,VT_s_avg_k)
       !
       !--------------------------------------------------------
       ! assign obs variables for obs_sequence
       !--------------------------------------------------------
       !
       ! location
             latitude=lat 
             if (lon<0) then
                longitude=lon+360
             else
                longitude=lon
             endif

        ! time (get time from sec MOPITT variable)
             hour1 = int(sec/3600d0)
             minute = int( (sec-hour1*3600d0)/60d0)
             second = int(sec - hour1*3600d0 - minute*60d0)
             if ( hour == 24 ) then
                if (sec <= 21.01*3600d0) then
                   day1 = day+1
                   if (day1 > days_in_month(month)) then
                      day1 = 1
                      if (month < 12) then
                         month1 = month + 1
                         year1 = year
                      else
                         month1 = 1
                         year1  = year+1
                      endif
                   else
                      month1 = month
                      year1 = year
                   endif
                else
                   day1 = day
                   month1 = month
                   year1 = year
               endif
             else
                day1 = day
                month1 = month
                year1 = year
             endif
             obs_time=set_date(year1,month1,day1,hour1,minute,second)

             call get_time(obs_time, seconds, days)
!
! APM: Write raw MOPITT data to file
!             write(fileidapm,'(a15,2i10,16e14.6)') transform_typ,days,seconds,latitude,longitude,nlevels,dofs
!             write(fileidapm,'(10e14.6)') mop_prs 
!             write(fileidapm,'(10e14.6)') x_r
!             write(fileidapm,'(10e14.6)') x_p
!             write(fileidapm,'(100e14.6)') avg_k
!             write(fileidapm,'(100e14.6)') cov_a
!             write(fileidapm,'(100e14.6)') cov_r
!             write(fileidapm,'(100e14.6)') cov_m
!
        !--------------------------------------------------------
        ! Loop through the mop_dim levels for now
        ! Use each mixing ratio as a separate obs
        !--------------------------------------------------------
        
             do kk = 1, nlev
 
                qc_count=qc_count+1
                level = mopittlev(kk)*100

            ! since this is already in spectral space, there is no vertical location
            ! set which_vert to VERTISUNDEF=-2
            ! level doesnt have value here, replace level passed to set_location to kk
                which_vert=-2
                obs_location=set_location(longitude, latitude, level, which_vert)
                co_psurf=mopittlev(1)*100               !hPa to Pa
                co_avgker(1:nlvls)=avgker(kk,1:nlvls)               !unitless
                co_prior=xapr(kk)  
                co_vmr(1)=xcomp(kk)
                err = xcomperr(kk)
                co_error=err*err

                obs_kind             = MOPITT_CO_RETRIEVAL

                call set_obs_def_kind(obs_def, obs_kind)
                call set_obs_def_location(obs_def, obs_location)
                call set_obs_def_time(obs_def, obs_time)
                call set_obs_def_error_variance(obs_def, co_error)
                call set_obs_def_mopitt_co(qc_count, co_avgker, co_prior, co_psurf, &
                                     nlvls)
                call set_obs_def_key(obs_def, qc_count)

                call set_obs_values(obs, co_vmr, 1)
                call set_qc(obs, co_qc, num_qc)
                call set_obs_def(obs, obs_def)

                if ( qc_count == 1 ) then
                   call insert_obs_in_seq(seq, obs)
                else
                   call insert_obs_in_seq(seq, obs, obs_old )
                endif
                obs_old=obs

             enddo    !k loop
             deallocate(xcomp,xcomperr,xapr,avgker,mopittlev)
          else   
             deallocate(s_cvr,s_cva,s_cov_r,s_cov_a,s_avg_k,s_x_r,s_x_p)
             deallocate(rs_cov_r,rs_cov_a,rs_avg_k,rs_x_r,rs_x_p)
             deallocate(rs_avg_xp,err2_rs_r,err2_rs_a)
             deallocate(Z,ZT,SV_s_avg_k,U_s_avg_k)
             deallocate(UT_s_avg_k,VT_s_avg_k)
          endif    ! qstatus
       endif    ! co_qc(1)

    ! read next data point
      deallocate(mop_prs,x_r,x_p,avg_k,cov_a,cov_r,cov_m)
      read(fileid,*,iostat=ios) transform_typ, sec, lat, lon, nlevels, dofs 
      nlvls=int(nlevels)
   enddo    !ios

  !----------------------------------------------------------------------
  ! Write the sequence to a file
  !----------------------------------------------------------------------
   if  (bin_beg == 3.01) then
      file_name=trim(file_name)//chr_year//chr_month//chr_day//'06'
   elseif (bin_beg == 9.01) then
      file_name=trim(file_name)//chr_year//chr_month//chr_day//'12'
   elseif (bin_beg == 15.01) then
      file_name=trim(file_name)//chr_year//chr_month//chr_day//'18'
   elseif (bin_beg == 21.01) then
      file_name=trim(file_name)//chr_year//chr_month//chr_day//'24'
   endif !bin

   call write_obs_seq(seq, file_name)

999 continue
   close(fileid)
!   close(fileidapm)

!-----------------------------------------------------------------------------
! Clean up
!-----------------------------------------------------------------------------
   call timestamp(string1=source,string2=revision,string3=revdate,pos='end')

end program create_mopitt_obs_sequence
!
    subroutine mat_prd(A_mat,B_mat,C_mat,na,ma,nb,mb)
!
! compute dot product of two matrics
    integer :: ma,na,mb,nb,i,j,k
    double precision :: A_mat(na,ma),B_mat(nb,mb),C_mat(na,mb)
!
! check that na=mb
    if(ma .ne. nb) then
       print *, 'Error in matrix dimension ma (cols) must equal nb (rows) ',ma,' ',nb
       stop
    endif
!
! initialze the product array
    C_mat(:,:)=0.
!
! calculate inner product
    do i=1,na
       do j=1,mb
          do k=1,mb
             C_mat(i,j)=C_mat(i,j)+A_mat(i,k)*B_mat(k,j) 
          enddo
       enddo
    enddo
    return
    end subroutine mat_prd
!
    subroutine mat_tri_prd(A_mat,B_mat,C_mat,D_mat,na,ma,nb,mb,nc,mc)
!
! compute dot product of three matrics D=A*B*C
    integer :: na,ma,nb,mb,nc,mc,i,j,k
    double precision :: A_mat(na,ma),B_mat(nb,mb),C_mat(nc,mc),D_mat(na,mc)
    double precision :: Z_mat(nb,mc)
!
! check that na=mb
    if(ma .ne. nb) then
       print *, 'Error in matrix dimension ma (cols) must equal nb (rows) ',ma,' ',nb
       stop
    endif
    if(mb .ne. nc) then
       print *, 'Error in matrix dimension mb (cols) must equal nc (rows) ',mb,' ',nc
       stop
    endif
!
! initialze the product array
    Z_mat(:,:)=0.
    D_mat(:,:)=0.
!
! calculate first inner product Z=B*C
    do i=1,nb
       do j=1,mc
          do k=1,mb
             Z_mat(i,j)=Z_mat(i,j)+B_mat(i,k)*C_mat(k,j) 
          enddo
       enddo
    enddo
!
! calculate second inner product D=A*Z
    do i=1,na
       do j=1,mc
          do k=1,ma
             D_mat(i,j)=D_mat(i,j)+A_mat(i,k)*Z_mat(k,j) 
          enddo
       enddo
    enddo
    return
    end subroutine mat_tri_prd
!
    subroutine vec_to_mat(a_vec,A_mat,n)
!
! compute dot product of two matrics
    integer :: n,i
    double precision :: a_vec(n),A_mat(n,n)
!
! initialze the product array
    A_mat(:,:)=0.
!
! calculate inner product
    do i=1,n
       A_mat(i,i)=a_vec(i) 
    enddo
    return
    end subroutine vec_to_mat
!
    subroutine diag_inv_sqrt(A_mat,n)
!
! calculate inverse square root of diagonal elements
    integer :: n,i
    double precision :: A_mat(n,n)
    do i=1,n
       if(A_mat(i,i).le.0.) then
          print *, 'Error in Subroutine vec_to_mat arg<=0 ',i,' ',A_mat(i,i)
          call abort
       endif
       A_mat(i,i)=1./sqrt(A_mat(i,i)) 
    enddo
    return
    end subroutine diag_inv_sqrt
!
    subroutine diag_sqrt(A_mat,n)
!
! calculate square root of diagonal elements
    integer :: n,i
    double precision :: A_mat(n,n)
    do i=1,n
       if(A_mat(i,i).lt.0.) then
          print *, 'Error in Subroutine vec_to_mat arg<0 ',i,' ',A_mat(i,i)
          call abort
       endif
       A_mat(i,i)=sqrt(A_mat(i,i)) 
    enddo
    return
    end subroutine diag_sqrt
!
    subroutine lh_mat_vec_prd(SCL_mat,a_vec,s_a_vec,n)
!
! calculate left hand side scaling of column vector
    integer :: n,i,j
    double precision :: SCL_mat(n,n),a_vec(n),s_a_vec(n)
!
! initialize s_a_vec
    s_a_vec(:)=0.
!
! conduct scaling
    do i=1,n
       do j=1,n
          s_a_vec(i)=s_a_vec(i)+SCL_mat(i,j)*a_vec(j)
       enddo 
    enddo
    return
    end subroutine lh_mat_vec_prd
!
    subroutine rh_vec_mat_prd(SCL_mat,a_vec,s_a_vec,n)
!
! calculate right hand side scaling of a row vector
    integer :: n,i,j
    double precision :: SCL_mat(n,n),a_vec(n),s_a_vec(n)
!
! initialize s_a_vec
    s_a_vec(:)=0.
!
! conduct scaling
    do i=1,n
       do j=1,n
          s_a_vec(i)=s_a_vec(i)+a_vec(j)*SCL_mat(j,i) 
       enddo
    enddo
    return
    end subroutine rh_vec_mat_prd
!
    subroutine mat_transpose(A_mat,AT_mat,n,m)
!
! calculate matrix transpose
    integer :: n,m,i,j
    double precision :: A_mat(n,m),AT_mat(m,n)
    do i=1,n
       do j=1,m
          AT_mat(j,i)=A_mat(i,j) 
       enddo
    enddo
    return
    end subroutine mat_transpose
!
    subroutine diag_vec(A_mat,a_vec,n)
!
! calculate square root of diagonal elements
    integer :: n,i
    double precision :: A_mat(n,n),a_vec(n)
    do i=1,n
       a_vec(i)=A_mat(i,i) 
    enddo
    return
    end subroutine diag_vec
!
