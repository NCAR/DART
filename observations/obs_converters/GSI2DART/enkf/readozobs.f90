module readozobs
!$$$  module documentation block
!
! module: readozobs                    read ozone data from diag_sbuv2* files.
!
! prgmmr: whitaker         org: esrl/psd               date: 2009-02-23
!
! abstract: read ozone data from diag_sbuv2* files written out
!  by GSI forward operator code.
!
! Public Subroutines:
!  get_num_ozobs: determine the number of observations to read.
!  get_ozobs_data: read the data.
!   
! Public Variables: None
!
! program history log:
!   2009-02-23  Initial version.
!
! attributes:
!   language: f95
!
!$$$

use kinds, only: r_single,i_kind,r_kind
use params, only: nsats_oz,sattypes_oz
  
implicit none

private
public :: get_num_ozobs, get_ozobs_data

contains

subroutine get_num_ozobs(obspath,datestring,num_obs_tot,id)
    character (len=500), intent(in) :: obspath
    character (len=10), intent(in) :: datestring
    character(len=500) obsfile
    character(len=8), intent(in) :: id
    integer(i_kind) :: nlevs  ! number of levels (layer amounts + total column) per obs   
    character(20) :: isis     ! sensor/instrument/satellite id
    character(10) :: obstype  !  type of ozone obs
    character(10) :: dplat    ! sat sensor
    real(r_single), allocatable, dimension(:) :: err,grs,pob
    real(r_single),allocatable,dimension(:,:)::diagbuf
    real(r_single),allocatable,dimension(:,:,:)::rdiagbuf
    real(r_kind) :: errorlimit,errorlimit2
    integer(i_kind),allocatable,dimension(:,:)::idiagbuf
    integer(i_kind) iunit,jiter,ii,ireal,iint,iextra,idate,ios,nsat,n,k
    integer(i_kind), intent(out) :: num_obs_tot
    integer(i_kind), allocatable, dimension(:) :: iouse
    integer(i_kind):: nread,nkeep
    logical :: fexist
    iunit = 7
    num_obs_tot = 0
!   make consistent with screenobs
    errorlimit=1._r_kind/sqrt(1.e9_r_kind)
    errorlimit2=1._r_kind/sqrt(1.e-6_r_kind)
    do nsat=1,nsats_oz
        nread = 0
        nkeep = 0
        obsfile = trim(adjustl(obspath))//"diag_"//trim(sattypes_oz(nsat))//"_ges."//datestring//'_'//trim(adjustl(id))
        inquire(file=obsfile,exist=fexist)
        if (.not. fexist .or. datestring .eq. '0000000000') then
        obsfile = trim(adjustl(obspath))//"diag_"//trim(sattypes_oz(nsat))//"_ges."//trim(adjustl(id))
        endif
        !print *,obsfile
        open(iunit,form="unformatted",file=obsfile,iostat=ios)
        read(iunit,err=20,end=30) isis,dplat,obstype,jiter,nlevs,idate,iint,ireal,iextra
        allocate(pob(nlevs),grs(nlevs),err(nlevs),iouse(nlevs))
        read(iunit,err=20,end=30) pob,grs,err,iouse
10      continue
        read(iunit,err=20,end=30) ii
        allocate(idiagbuf(iint,ii))
        allocate(diagbuf(ireal,ii))
        allocate(rdiagbuf(6,nlevs,ii))
        read(iunit,err=20,end=30) idiagbuf,diagbuf,rdiagbuf
        do k=1,nlevs
          nread=nread+ii
          if (iouse(k) < 0 .or. pob(k) <= 0.001 .or. &
              pob(k) > 1200._r_kind) cycle
          do n=1,ii
            if (rdiagbuf(3,k,n) <= errorlimit .or.  &
                rdiagbuf(3,k,n) >= errorlimit2 .or.  &
                abs(rdiagbuf(1,k,n)) > 1.e9_r_kind) cycle
            nkeep = nkeep + 1
            num_obs_tot = num_obs_tot + 1
          end do
        end do
        deallocate(idiagbuf,diagbuf,rdiagbuf)
        go to 10
20      continue
        print *,'error reading diag_sbuv file'
30      continue
        if(allocated(pob))deallocate(pob,grs,err,iouse)
        close(iunit)
        write(6,100) nsat,trim(sattypes_oz(nsat)),nread,nkeep,num_obs_tot
100     format(2x,i3,2x,a20,2x,'nread= ',i9,2x,'nkeep= ',i9,2x,'num_obs_tot= ',i9)
    enddo
    print *,num_obs_tot,' ozone obs'
    close(iunit)
end subroutine get_num_ozobs

subroutine get_ozobs_data(obspath, datestring, nobs_max, h_x, h_xnobc, x_obs, x_err, &
           x_lon, x_lat, x_press, x_time, x_code, x_errorig, x_type, id,id2)

  character*500, intent(in) :: obspath
  character*500 obsfile,obsfile2
  character*10, intent(in) :: datestring
  character(len=8), intent(in) :: id,id2

  integer(i_kind) :: nlevs  ! number of levels (layer amounts + total column) per obs   
  character(20) :: isis,isis2     ! sensor/instrument/satellite id
  character(10) :: obstype,obstype2  !  type of ozone obs
  character(10) :: dplat,dplat2    ! sat sensor
  integer(i_kind) iunit,jiter,ii,ireal,iint,iextra,idate,nob,n,ios,nobs_max,nsat,k
  integer(i_kind) iunit2,jiter2,nlevs2,idate2,iint2,ireal2,iextra2,ii2

  real(r_single), dimension(nobs_max) :: h_x,h_xnobc,x_obs,x_err,x_lon,&
                               x_lat,x_press,x_time,x_errorig
  integer(i_kind), dimension(nobs_max) :: x_code
  character(len=20), dimension(nobs_max) ::  x_type

  real(r_single),allocatable,dimension(:,:)::diagbuf,diagbuf2
  real(r_single),allocatable,dimension(:,:,:)::rdiagbuf,rdiagbuf2
  integer(i_kind),allocatable,dimension(:,:)::idiagbuf,idiagbuf2
  real(r_single), allocatable, dimension(:) :: err,grs,pob
  real(r_single), allocatable, dimension(:) :: err2,grs2,pob2
  integer(i_kind), allocatable, dimension(:) :: iouse,iouse2
  logical twofiles, fexist
  real(r_kind) :: errorlimit,errorlimit2
 
! make consistent with screenobs
  errorlimit=1._r_kind/sqrt(1.e9_r_kind)
  errorlimit2=1._r_kind/sqrt(1.e-6_r_kind)

  twofiles = id /= id2
  iunit = 7
  iunit2 = 17
  nob = 0
  do nsat=1,nsats_oz
      obsfile = trim(adjustl(obspath))//"diag_"//trim(sattypes_oz(nsat))//"_ges."//datestring//'_'//trim(adjustl(id))
      inquire(file=obsfile,exist=fexist)
      if (.not. fexist .or. datestring .eq. '0000000000') then
      obsfile = trim(adjustl(obspath))//"diag_"//trim(sattypes_oz(nsat))//"_ges."//trim(adjustl(id))
      endif
      !print *,obsfile
      open(iunit,form="unformatted",file=obsfile,iostat=ios)
      read(iunit,err=20,end=30) isis,dplat,obstype,jiter,nlevs,idate,iint,ireal,iextra
      allocate(pob(nlevs),grs(nlevs),err(nlevs),iouse(nlevs))
      read(iunit,err=20,end=30) pob,grs,err,iouse
      if(twofiles)then
        obsfile2 = trim(adjustl(obspath))//"diag_"//trim(sattypes_oz(nsat))//"_ges."//datestring//'_'//trim(adjustl(id2))
        inquire(file=obsfile,exist=fexist)
        if (.not. fexist .or. datestring .eq. '0000000000') then
        obsfile2 = trim(adjustl(obspath))//"diag_"//trim(sattypes_oz(nsat))//"_ges."//trim(adjustl(id2))
        endif
        !print *,obsfile2
        open(iunit2,form="unformatted",file=obsfile2,iostat=ios)
        read(iunit2,err=20,end=30) isis2,dplat2,obstype2,jiter2,nlevs2,idate2,iint2,ireal2,iextra2
        if(isis /= isis2 .or. dplat /= dplat2 .or. obstype /= obstype2 .or. jiter /= jiter2 .or. &
           nlevs /= nlevs2 .or. idate /= idate2 .or. iint /= iint2 .or. ireal /= ireal2)then
           write(6,*) 'inconsistency in ozone files'
           write(6,*) 'isis',isis,isis2
           write(6,*) 'dplat',dplat,dplat2
           write(6,*) 'obstype',obstype,obstype2
           write(6,*) 'jiter',jiter,jiter2
           write(6,*) 'nlevs',nlevs,nlevs2
           write(6,*) 'idate',idate,idate2
           write(6,*) 'iint',iint,iint2
           write(6,*) 'ireal',ireal,ireal2
           call stop2(66)
        end if
        allocate(pob2(nlevs),grs2(nlevs),err2(nlevs),iouse2(nlevs))
        read(iunit2,err=20,end=30) pob2,grs2,err2,iouse2
        do k=1,nlevs
          if(pob(k) /= pob2(k) .or. grs(k) /= grs2(k) .or. err(k) /= err2(k) .or. &
             iouse(k) /= iouse2(k))then
             write(6,*) ' ozone file vertical inconsistency level = ',k
             write(6,*) 'pob',pob(k),pob2(k)
             write(6,*) 'grs',grs(k),grs2(k)
             write(6,*) 'err',err(k),err2(k)
             write(6,*) 'iouse',iouse(k),iouse2(k)
             call stop2(67)
           end if
        end do
      end if
10    continue
      read(iunit,err=20,end=30) ii
      allocate(idiagbuf(iint,ii))
      allocate(diagbuf(ireal,ii))
      allocate(rdiagbuf(6,nlevs,ii))
      allocate(rdiagbuf2(6,nlevs,ii))
      read(iunit,err=20,end=30) idiagbuf,diagbuf,rdiagbuf
      if(twofiles)then
        read(iunit2,err=20,end=30) ii2
        if(ii /= ii2)then
          write(6,*) 'ii inconsistency in ozone ',ii,ii2
          call stop2(68)
        end if
        allocate(idiagbuf2(iint,ii2))
        allocate(diagbuf2(ireal,ii2))
        read(iunit2,err=20,end=30) idiagbuf2,diagbuf2,rdiagbuf2
      else
         do n=1,ii
           do k=1,nlevs
             rdiagbuf2(2,k,n)=rdiagbuf(2,k,n)
           end do
         end do
      end if

      do n=1,ii
         if(twofiles)then
           if(diagbuf(1,n) /= diagbuf2(1,n) .or. diagbuf(2,n) /= diagbuf2(2,n))then
             write(6,*) 'lat lon inconsistency in ozone '
             write(6,*) 'lat',diagbuf(1,n),diagbuf2(1,n)
             write(6,*) 'lon',diagbuf(2,n),diagbuf2(2,n)
           end if
         end if
         do k=1,nlevs
           ! if pressure is zero, it's a total column ob so don't use it.
           if(iouse(k) < 0 .or. pob(k) <= 0.001_r_kind .or. &
              pob(k) > 1200._r_kind)cycle
           if (rdiagbuf(3,k,n) <= errorlimit .or. &
               rdiagbuf(3,k,n) >= errorlimit2 .or. &
               abs(rdiagbuf(1,k,n)) > 1.e9_r_kind) cycle
           nob = nob + 1
           x_code(nob) = 700 + k ! made up code for ozone level k
           x_lat(nob) = diagbuf(1,n)
           x_lon(nob) = diagbuf(2,n)
           !print *,n,k,pob(k)
           x_press(nob) = pob(k)
           x_time(nob) = diagbuf(3,n)
           x_err(nob) = (1./rdiagbuf(3,k,n))**2
           x_errorig(nob) = x_err(nob)
           x_obs(nob) = rdiagbuf(1,k,n)
           h_xnobc(nob) = rdiagbuf(1,k,n)-rdiagbuf2(2,k,n)
           h_x(nob) = rdiagbuf(1,k,n)-rdiagbuf(2,k,n)
           x_type(nob) = ' oz                 '
         enddo
      enddo
      deallocate(diagbuf,rdiagbuf,rdiagbuf2,idiagbuf)
      if(twofiles)deallocate(diagbuf2,idiagbuf2)
      go to 10
20    continue
      print *,'error reading diag_sbuv file'
30    continue
      close(iunit)
      if(twofiles) close(iunit2)
      if (allocated(pob)) deallocate(pob,err,grs,iouse)
      if (allocated(pob2)) deallocate(pob2,err2,grs2,iouse2)
  enddo
  if (nob /= nobs_max) then
      print *,'number of obs not what expected in get_ozobs_data',nob,nobs_max
      call stop2(93)
  end if

 end subroutine get_ozobs_data

end module readozobs
