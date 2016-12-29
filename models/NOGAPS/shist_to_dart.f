c This code may (or may not) be part of the NOGAPS distribution,
c So it is not protected by the DART copyright agreement.
c
c DART $Id$
c
      subroutine shist_to_dart(lgtrdy,lnmode,lfcst,loutp,ptime,
     &             dart_out_size,dart_out,pt_mean_size,pt_mean)
c
cc
      use gcvars
      use dyncons
      use cubic
      use fields
      use navdas
      use times
      use types_mod, only : r8   

c
      implicit none
c
      include 'imp.h'
      include 'dbms.h'
      include 'const.h'
      include 'phctrl.h'
      include 'dgdir.h'
 
c Parameters:
      logical lgtrdy
      logical lnmode
      logical lfcst
      logical loutp
      logical ptime
      integer dart_out_size
      real(r8) :: dart_out(dart_out_size)
      integer pt_mean_size
      real(r8) :: pt_mean(pt_mean_size)
 
c Local variables:
      integer i
      integer it
      integer istat
      integer itaux
      integer itdum
      integer ith
      integer itx
      integer k
      integer locf
      integer lpath
      integer ml
      integer nspec
      real cpuoel
      real ograv
      real osquv
      real tpcnt
      real rts2pt(2)
c     real rtime1(2,25)
c     real rtime2(2,25)
      real rtime1(2,30)
      real rtime2(2,30)
      real rtimmpi1(2,10)
      real sr(3)
      logical onlyo
c      logical restart
      character(len=10) ipdto
c     character(len=8) tnm(25)
      character(len=8) tnm(30)
 
c----------------------------------------------------------------
c----------------------------------------------------------------

      INTEGER ij,np,iijj, itime, jpt, ipt, k1,k2,k3
      INTEGER j, mem, gridpt, io, index, ierr, nmem 
      INTEGER hlev, index_loc

      INTEGER lwr_bnd, upr_bnd,lwr_flip,upr_flip,nobs

      INTEGER ijpts, numproc, max_mem, nZfcols,Nmem_loc
      PARAMETER( ijpts = im*jm )
      INTEGER ngridpts, Ndim, Ndimloc
      REAL dsm, dsv, lat, lon, temp1, angle, ampl
      PARAMETER( ngridpts = ijpts*lm )
      PARAMETER( Ndim = ngridpts )

      INTEGER lendir
      CHARACTER(len=10) yyyymmddhh
      CHARACTER(len=90) infile 
      CHARACTER(len=72) shistm
      CHARACTER(len=20) file_out
      CHARACTER(len=2) mem_id

      REAL, DIMENSION(:,:), ALLOCATABLE ::th
      REAL, DIMENSION(:), ALLOCATABLE ::rh

      REAL, DIMENSION(:,:), ALLOCATABLE :: u_mean
      REAL, DIMENSION(:,:), ALLOCATABLE :: v_mean
      REAL, DIMENSION(:,:), ALLOCATABLE :: te_mean
      REAL, DIMENSION(:,:), ALLOCATABLE :: th_mean
      REAL, DIMENSION(:,:), ALLOCATABLE :: sht_mean
      REAL, DIMENSION(:,:), ALLOCATABLE :: phi_mean
c      REAL, DIMENSION(:),   ALLOCATABLE :: pt_mean

      REAL, DIMENSION(:), ALLOCATABLE :: glob_u
      REAL, DIMENSION(:), ALLOCATABLE :: glob_v
      REAL, DIMENSION(:), ALLOCATABLE :: glob_te
      REAL, DIMENSION(:), ALLOCATABLE :: glob_th
      REAL, DIMENSION(:), ALLOCATABLE :: glob_sht
      REAL, DIMENSION(:), ALLOCATABLE :: glob_pt
c      REAL, DIMENSION(:), ALLOCATABLE :: dart_out

      REAL, DIMENSION(:),   ALLOCATABLE :: rnorm

      REAL, DIMENSION(:,:,:), ALLOCATABLE :: temc

      REAL, DIMENSION(:,:), ALLOCATABLE :: dsqgeot

c----------------------------------------------------------------
c----------------------------------------------------------------

c END OF DEFINITIONS
 
c
c     This subroutine controls the four basic functions of the spectral
c     model: (1) interpolation of analysis fields from p surfaces to the
c     model's coodinates, (2) normal mode initialization, (3) forecast,
c     and (4) transformation and interpolation of model fields back to
c     standard fnoc mandatory pressure surface field format.
c
c     data tnm/'job time', 'fcst mod', 'nmode   ', 'sig2p   ',
c    &     'p2sig   ', 'getrdy  ', 'rstran  ', 'srtran  ', 'miscell ',
c    &     'bcubiy  ', 'bcubij  ', 'diabatic', 'dnostics', 'implicit',
c    &     'bcubver ', 'lsprec  ', 'cumulus ', 'longwrad', 'shrtwrad',
c    &     'pbl     ', 'shalocum', 'histrdwt', 'allocate', 'dbasread',
c    &     'dbaswrit'/
c
      data tnm/'job time', 'fcst mod', 'nmode   ', 'sig2p   ',
     &     'p2sig   ', 'getrdy  ', 'rstran  ', 'srtran  ', 'miscell ',
     &     'bcubiy  ', 'bcubij  ', 'diabatic', 'dnostics', 'implicit',
     &     'bcubver ', 'lsprec  ', 'cumulus ', 'longwrad', 'shrtwrad',
     &     'pbl     ', 'shalocum', 'histrdwt', 'allocate', 'dbasread',
     &     'dbaswrit', 'histsort', 'histgath', 'histscat',
     &     'gwdrag  ', 'cld&scmx'/
c
c
c
c **********************************************************************
c **********************************************************************
c
c start the model !!!!!
c

      call clokon(rjobtm)
c
c nproc is the number of message-passing nodes
c jsplit is the number of multi-tasking 'chunks' per node
c jsplit can also be used to reduce size of 'chunks' for better cache use
c e.g., T3E : nproc > 1, jsplit = 1
c       C90 : nproc = 1, jsplit > 1
c total no. of chunks = nproc*jsplit
c
c
c  mpi constants and pointer arrays
c
      call clokon(rtmisc)
c
      call consio(ptime,loutp,lgtrdy,lnmode)
c
      i = (locf(cpad(10))-locf(tofday))
      call mpi_bcast(tofday,i,mpi_byte,0,mpicomm,istat)
c
      call consallpe
c
c loutp and itau might have been changed in cons
c
      call mpi_bcast(loutp,1,mpilogical,0,mpicomm,istat)
      call mpi_bcast(itau,1,mpiint,0,mpicomm,istat)
c
c write model grid geometry info to eye readable file
c
      if( ir == 1 )  call writasc(hstdir,istat)
c
      nspec = 4*lm + 1
c
      call clokof(rtmisc)
c
c **********************************************************************
      call clokon(rtalloc)
c
c  allocate memory
c
      allocate(gaunow(ijsize,6*lm+3))
      allocate(c3vars(ijsize,nc3grd))
      allocate(specvar(mlsize,2,nshist))
      allocate(specten(mlsize,2,nspec))
      allocate(tendsnl(ijsize,8*lm+1))
crp      allocate(prcpvars(ijsize,4))
c
c **********************************************************************
c
c assign pointers for two time levels of spectral variables
c
      vornow => specvar(1:mlsize,1:2,1:lm)
      divnow => specvar(1:mlsize,1:2,1+lm:2*lm)
      temnow => specvar(1:mlsize,1:2,1+lm*2:3*lm)
      shnow => specvar(1:mlsize,1:2,1+lm*3:4*lm)
      plnow => specvar(:,:,1+4*lm)
c
      vorold = >specvar(1:mlsize,1:2,2+4*lm:1+5*lm)
      divold = >specvar(1:mlsize,1:2,2+5*lm:1+6*lm)
      temold = >specvar(1:mlsize,1:2,2+6*lm:1+7*lm)
      shold = >specvar(1:mlsize,1:2,2+7*lm:1+8*lm)
      plold = >specvar(:,:,2+8*lm)
      dsqgeo = >specvar(:,:,3+8*lm)
c
      vorten = >specten(1:mlsize,1:2,1:lm)
      divten = >specten(1:mlsize,1:2,1+lm:2*lm)
      temten = >specten(1:mlsize,1:2,1+lm*2:3*lm)
      shten = >specten(1:mlsize,1:2,1+lm*3:4*lm)
      plten = >specten(:,:,1+4*lm)
c
c
c assign pointers for current time gaussian grid prognostic variables
c
      ut = >gaunow(1:ijsize,1:lm)
      vt = >gaunow(1:ijsize,1+lm:2*lm)
      rvor = >gaunow(1:ijsize,1+2*lm:3*lm)
      rdiv = >gaunow(1:ijsize,1+3*lm:4*lm)
      tt = >gaunow(1:ijsize,1+4*lm:5*lm)
      sht = >gaunow(1:ijsize,1+5*lm:6*lm)
      pt = >gaunow(:,1+6*lm)
      dlpl = >gaunow(:,2+6*lm)
      dtpl = >gaunow(:,3+6*lm)
c
c *********************************************************
c
c assign pointers to gaussian grid  non-linear tendency terms
c
      utenvd = >tendsnl(1:ijsize,1:lm)
      uttem = >tendsnl(1:ijsize,lm+1:2*lm)
      utsht = >tendsnl(1:ijsize,2*lm+1:3*lm)
      vtenvd = >tendsnl(1:ijsize,3*lm+1:4*lm)
      vttem = >tendsnl(1:ijsize,4*lm+1:5*lm)
      vtsht = >tendsnl(1:ijsize,5*lm+1:6*lm)
      stten = >tendsnl(1:ijsize,6*lm+1:7*lm)
      sshten = >tendsnl(1:ijsize,7*lm+1:8*lm)
      pten = >tendsnl(1:ijsize,8*lm+1)
c
c  assign pointers for precip variables
crp      cuprec3 = >prcpvars(:,1)
crp      precls3 = >prcpvars(:,2)
crp      cuprec6 = >prcpvars(:,3)
crp      precls6 = >prcpvars(:,4)
c
c  assign pointers for surface variables
c
      sgeo = >c3vars(:,1)
      sgphis = >c3vars(:,2)
      gt = >c3vars(:,3)
      gwclim = >c3vars(:,4)
      gtclim = >c3vars(:,5)
      gw = >c3vars(:,6)
      alb = >c3vars(:,7)
      z0 = >c3vars(:,8)
      sn = >c3vars(:,9)
      ss = >c3vars(:,10)
      rs = >c3vars(:,11)
      fws = >c3vars(:,12)
      fhs = >c3vars(:,13)
      rmb = >c3vars(:,14)
      drag = >c3vars(:,15)
      cuprec = >c3vars(:,16)
      precls = >c3vars(:,17)
      prrate = >c3vars(:,18)
      dswtop = >c3vars(:,19)
      ulwtop = >c3vars(:,20)
      plcl = >c3vars(:,21)
      cumtop = >c3vars(:,22)
      cufrac = >c3vars(:,23)
      conice = >c3vars(:,24)
      uzmlxy = >c3vars(:,25)
      dtt = >c3vars(1:ijsize,26:25+lm)
      pstend = >c3vars(:,26+lm)
      sr_terrain = >c3vars(:,27+lm)
      rgc = >c3vars(:,28+lm)
      ulm = >c3vars(:,29+lm)
      vlm = >c3vars(:,30+lm)
      qlm = >c3vars(:,31+lm)
      zlm = >c3vars(:,32+lm)
c  assign pointers for precip variables
      cuprec3 = >c3vars(:,33+lm)
      precls3 = >c3vars(:,34+lm)
      cuprec6 = >c3vars(:,35+lm)
      precls6 = >c3vars(:,36+lm)
c
      ts = >c3vars(:,37+lm)
c
      tsmin = >c3vars(:,38+lm)
      tsmax = >c3vars(:,39+lm)
      freezing_rain = >c3vars(:,40+lm)
      frozen_rain = >c3vars(:,41+lm)
      snowtot = >c3vars(:,42+lm)
      surdpd = >c3vars(:,43+lm)
      slpmean = >c3vars(:,44+lm)
      snowage = >c3vars(:,45+lm)
c
c *****************************************
c
      allocate(phi(ijsize,lm))
      allocate(plt(ijsize,lm))
      allocate(pk(ijsize,lm))
      allocate(pk2(ijsize,lm))
      allocate(sd(ijsize,lm))
      allocate(igc(ijsize))
c
      call clokof(rtalloc)

c------------------------------------------------------------------------

c------------------------------------------------------------------------
  
      call chlen(pstdir1,lendir)
      if(pstdir1(lendir:lendir) == '/')lendir = lendir - 1

c     if(ir == 1) then
c        infile = pstdir0(1:lendir)//'/yyyymmddhh.out'
c        OPEN(UNIT = 120, FILE = infile, STATUS = 'OLD')
c        READ(120,'(1A10)') yyyymmddhh
c        CLOSE(UNIT = 120)
c     end if 

c     call mpi_bcast(yyyymmddhh,10,mpi_byte,0,mpicomm,istat)

c -- allocate variables for ensemble

      ALLOCATE(glob_u(ijpts))
      ALLOCATE(glob_v(ijpts))
      ALLOCATE(glob_te(ijpts))
      ALLOCATE(glob_th(ijpts))
      ALLOCATE(glob_sht(ijpts))
      ALLOCATE(glob_pt(ijpts))
             
      ALLOCATE(u_mean(ijpts,lm))
      ALLOCATE(v_mean(ijpts,lm))
      ALLOCATE(te_mean(ijpts,lm))
      ALLOCATE(th_mean(ijpts,lm))
      ALLOCATE(sht_mean(ijpts,lm))
c      ALLOCATE(pt_mean(ijpts))

      ALLOCATE(th(ijsize,lm))
      ALLOCATE(rh(ijsize))

      ALLOCATE(temc(mlsize,2,lm))

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  read in spectral history file and convert to Gaussian grid

      shistm=pstdir1(1:lendir)//'/shist'

      itau = 0

      call nfreads(shistm,mnf,lpakr,1,nshist,mlsiz,mlstrt,specvar,
     &                itau,cdtg,istat)

      call tranuv(im,jm,lm,jtrun,trigs,coslt,reps4,cim,poly,dpoly,
     &               vornow,divnow,ut,vt)  
      call transr(im,jm,lm,jtrun,trigs,poly,temnow,tt)  
      call transr(im,jm,lm,jtrun,trigs,poly,shnow,sht) 
      call transr_1lev(im,jm,jtrun,trigs,poly,plnow,pt)
      call prexp(ijsize,lm,asig,bsig,oddfac,pt,pk,pk2,plt) 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c convert angular velocity to winds and theta-v to temperature

       do k = 1,lm

            do i = 1,ijsize

               ut(i,k) = rad*sqrt(onocos(i))*ut(i,k)
               vt(i,k) = rad*sqrt(onocos(i))*vt(i,k)

               tt(i,k) = tt(i,k)*pk(i,k)/( 1.0 + p608
     &                 *( sht(i,k)/( 1. - sht(i,k) )  ) )

            end do

         end do ! end loop over k = 1,lm

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c gather data into global array

          do k = 1,lm

            call mpi_allgatherv(ut(1,k),ijsize,mpireal,glob_u,
     &                          ijsiz,ijstrt,mpireal,mpicomm,istat)
            call mpi_allgatherv(vt(1,k),ijsize,mpireal,glob_v,
     &                          ijsiz,ijstrt,mpireal,mpicomm,istat)
            call mpi_allgatherv(tt(1,k),ijsize,mpireal,glob_te,
     &                          ijsiz,ijstrt,mpireal,mpicomm,istat)
            call mpi_allgatherv(sht(1,k),ijsize,mpireal,glob_sht,
     &                          ijsiz,ijstrt,mpireal,mpicomm,istat)

            do i = 1, ijpts

               u_mean(i,k)  = glob_u(i)
               v_mean(i,k)  = glob_v(i)
               te_mean(i,k) = glob_te(i)
               sht_mean(i,k) = glob_sht(i)
   
            end do

         end do ! end loop over k = 1,lm

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        call mpi_allgatherv(pt,ijsize,mpireal,glob_pt,
     &                       ijsiz,ijstrt,mpireal,mpicomm,istat)

         pt_mean = glob_pt

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c unscramble each field and write to file

      if ( ir == 1 ) then

c       ! U, V, T, Q first ...


c        ALLOCATE(dart_out(4*ngridpts+ijpts))
        dart_out = 0.0
     
        do k = 1,lm   
      
         gridpt = 0

         do np = 1, nproc
          do iijj = 1, ijsiz(np)

            ij = (k-1)*ijpts + ijnx(iijj,np)

            gridpt = gridpt + 1

            dart_out(           ij) = u_mean(gridpt,k)
            dart_out(  ngridpts+ij) = v_mean(gridpt,k)
            dart_out(2*ngridpts+ij) = te_mean(gridpt,k)
            dart_out(3*ngridpts+ij) = sht_mean(gridpt,k) 

          end do
         end do
        end do
     
c       ! then the surface pressure

        gridpt = 0       
        do np = 1, nproc
          do iijj = 1, ijsiz(np)

            ij = ijnx(iijj,np)
            gridpt = gridpt + 1
            dart_out(4*ngridpts+ij) = pt_mean(gridpt)

          end do
         end do

       OPEN(unit=390, File='dart_state_vector.dat',
     &        access="direct",
     &        recl=4*( 4*ngridpts+ijpts ) )

       write(390,rec=1,iostat=io) dart_out
       close(390)

c      DEALLOCATE(dart_out)

       end if

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Read surface geopotential height from c3grid file and place into
c  the variable glob_pt ... then unscramble the vector and place
c  into pt_mean
c
         itau = 0
         shistm=pstdir1(1:lendir)//'/c3grid'

       call nfread(shistm,mnf,lpakr,1,nc3grd,ijsiz,ijstrt,c3vars,itau,
     &              idx,istat)

       call mpi_allgatherv(sgeo(1),ijsize,mpireal,glob_pt,  
     &                       ijsiz,ijstrt,mpireal,mpicomm,istat)

        gridpt = 0
        do np = 1, nproc
          do iijj = 1, ijsiz(np)

            ij = ijnx(iijj,np)
            gridpt = gridpt + 1
            pt_mean( ij ) = glob_pt( gridpt )

          end do
         end do

c ----------------------------------------------------------------
c ----------------------------------------------------------------
c ----------------------------------------------------------------

      return
      end

c <next few lines under version control, do not edit>
c $URL$
c $Id$
c $Revision$
c $Date$

