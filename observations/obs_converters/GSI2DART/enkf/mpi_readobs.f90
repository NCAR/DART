module mpi_readobs
!$$$  module documentation block
!
! module: mpi_readobs                  read obs, ob priors and associated
!                                      metadata if called from root task, 
!                                      otherwise receive data from root task.
!
! prgmmr: whitaker         org: esrl/psd               date: 2009-02-23
!
! abstract:
!
! Public Subroutines:
!  mpi_readobs: called by subroutine readobs in module enkf_obsmod. 
!   Read obs, ob priors and metadata from diag* files
!   created by GSI forward operator code and broadcast to all tasks.
!   
! Public Variables: None
!
! Modules Used:
!  readsatobs: to read satellite radiance diag* files.
!  readconvobs: to read diag_conv* files (obs from prepbufr file).
!  readozobs: to read diag_sbuv* ozone files.
!  mpisetup
!
! program history log:
!   2009-02-23  Initial version.
!
! attributes:
!   language: f95
!
!$$$
  
use kinds, only: r_kind, r_single, i_kind
use radinfo, only: npred
use readconvobs
use readsatobs
use readozobs
use mpisetup

implicit none

private
public :: mpi_getobs

contains

subroutine mpi_getobs(obspath, datestring, nobs_conv, nobs_oz, nobs_sat, nobs_tot, &
                      sprd_ob, ensmean_ob, ensmean_obbc, ob, &
                      oberr, oblon, oblat, obpress, &
                      obtime, oberrorig, obcode, obtype, &
                      biaspreds, anal_ob, indxsat, nanals)
    character*500, intent(in) :: obspath
    character*10, intent(in) :: datestring
    character(len=10) :: id,id2
    real(r_single), allocatable, dimension(:) :: ensmean_ob,ob,oberr,oblon,oblat,obpress,obtime,oberrorig,ensmean_obbc,sprd_ob
    integer(i_kind), allocatable, dimension(:) :: obcode,indxsat
    real(r_single), allocatable, dimension(:,:) :: biaspreds
    real(r_single), allocatable, dimension(:,:) :: anal_ob
    real(r_single), allocatable, dimension(:) :: h_xnobc
    real(r_single) :: analsi,analsim1
    character(len=20), allocatable,  dimension(:) ::  obtype
    integer(i_kind) nob, ierr, iozproc, isatproc, &
            nobs_conv, nobs_oz, nobs_sat, nobs_tot, nanal
    integer(i_kind), intent(in) :: nanals
    !include 'mpif.h'
    !integer mpi_status(mpi_status_size)
    iozproc=max(0,min(1,numproc-1))
    isatproc=max(0,min(2,numproc-2))
! get total number of conventional and sat obs for ensmean.
    id = 'ensmean'
    if(nproc == 0)call get_num_convobs(obspath,datestring,nobs_conv,id)
    if(nproc == iozproc)call get_num_ozobs(obspath,datestring,nobs_oz,id)
    if(nproc == isatproc)call get_num_satobs(obspath,datestring,nobs_sat,id)
    call mpi_bcast(nobs_conv,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(nobs_oz,1,mpi_integer,iozproc,mpi_comm_world,ierr)
    call mpi_bcast(nobs_sat,1,mpi_integer,isatproc,mpi_comm_world,ierr)
    if(nproc == 0)print *,'nobs_conv, nobs_oz, nobs_sat = ',nobs_conv,nobs_oz,nobs_sat
    nobs_tot = nobs_conv + nobs_oz + nobs_sat
! if nobs_tot != 0 (there were some obs to read)
    if (nobs_tot > 0) then
       if (nproc == 0) then
          ! this array only needed on root.
          allocate(anal_ob(nanals,nobs_tot))
       end if
       ! these arrays needed on all processors.
       allocate(h_xnobc(nobs_tot))
       allocate(sprd_ob(nobs_tot),ob(nobs_tot),oberr(nobs_tot),oblon(nobs_tot),&
       oblat(nobs_tot),obpress(nobs_tot),obtime(nobs_tot),oberrorig(nobs_tot),obcode(nobs_tot),&
       obtype(nobs_tot),ensmean_ob(nobs_tot),ensmean_obbc(nobs_tot),&
       biaspreds(npred+1, nobs_sat),indxsat(nobs_sat))
    else
! stop if no obs found (must be an error somewhere).
       print *,'no obs found!'
       call stop2(11)
    end if

! loop over ensmean and all ens members.
    nanal = nproc
    id = 'ensmean'
    id2=id
    if (nanal /= 0 .and. nanal <= nanals) then
        write(id2,'(a3,(i3.3))') 'mem',nanal
    end if
! read obs.
    if (nobs_conv > 0) then
! first nobs_conv are conventional obs.
      call get_convobs_data(obspath, datestring, nobs_conv,               &
        ensmean_obbc(1:nobs_conv), h_xnobc(1:nobs_conv), ob(1:nobs_conv), &
        oberr(1:nobs_conv), oblon(1:nobs_conv), oblat(1:nobs_conv),       &
        obpress(1:nobs_conv), obtime(1:nobs_conv), obcode(1:nobs_conv),   &
        oberrorig(1:nobs_conv), obtype(1:nobs_conv), id,id2)
!   if (nanal <= nanals) print *,id,id2, 'read conv obs'
    end if
    if (nobs_oz > 0) then
! second nobs_oz are conventional obs.
      call get_ozobs_data(obspath, datestring, nobs_oz,  &
        ensmean_obbc(nobs_conv+1:nobs_conv+nobs_oz),     &
        h_xnobc(nobs_conv+1:nobs_conv+nobs_oz),          &
        ob(nobs_conv+1:nobs_conv+nobs_oz),               &
        oberr(nobs_conv+1:nobs_conv+nobs_oz),            &
        oblon(nobs_conv+1:nobs_conv+nobs_oz),            &
        oblat(nobs_conv+1:nobs_conv+nobs_oz),            &
        obpress(nobs_conv+1:nobs_conv+nobs_oz),          &
        obtime(nobs_conv+1:nobs_conv+nobs_oz),           &
        obcode(nobs_conv+1:nobs_conv+nobs_oz),           &
        oberrorig(nobs_conv+1:nobs_conv+nobs_oz),        &
        obtype(nobs_conv+1:nobs_conv+nobs_oz), id,id2)
!   if (nanal <= nanals) print *,id,id2, 'read oz obs'
    end if
    if (nobs_sat > 0) then
      biaspreds = 0. ! initialize bias predictor array to zero.
! last nobs_sat are satellite radiance obs.
      call get_satobs_data(obspath, datestring, nobs_sat, &
        ensmean_obbc(nobs_conv+nobs_oz+1:nobs_tot),       &
        h_xnobc(nobs_conv+nobs_oz+1:nobs_tot),            &
        ob(nobs_conv+nobs_oz+1:nobs_tot),                 &
        oberr(nobs_conv+nobs_oz+1:nobs_tot),              &
        oblon(nobs_conv+nobs_oz+1:nobs_tot),              &
        oblat(nobs_conv+nobs_oz+1:nobs_tot),              &
        obpress(nobs_conv+nobs_oz+1:nobs_tot),            &
        obtime(nobs_conv+nobs_oz+1:nobs_tot),             &
        obcode(nobs_conv+nobs_oz+1:nobs_tot),             &
        oberrorig(nobs_conv+nobs_oz+1:nobs_tot),          &
        obtype(nobs_conv+nobs_oz+1:nobs_tot), biaspreds,indxsat,id,id2)
!     if (nanal <= nanals) print *,id,id2,'read sat obs'
    end if

    if (nanal <= nanals) then
     if (nanal == 0) then
        do nanal=1,nanals
           call mpi_recv(h_xnobc,nobs_tot,mpi_real4,nanal, &
                         1,mpi_comm_world,mpi_status,ierr)
           anal_ob(nanal,:) = h_xnobc(:)
        enddo
        analsi=1._r_single/float(nanals)
        analsim1=1._r_single/float(nanals-1)
!$omp parallel do private(nob,nanal)
        do nob=1,nobs_tot
! remove ensemble mean from each member.
           ensmean_ob(nob)  = sum(anal_ob(:,nob))*analsi
! ensmean_ob is unbiascorrected ensemble mean (anal_ob
           anal_ob(:,nob) = anal_ob(:,nob)-ensmean_ob(nob)
! compute sprd
           sprd_ob(nob) = sum(anal_ob(:,nob)**2)*analsim1
        enddo
!$omp end parallel do
       else ! nanal/nproc != 0
        ! send to root.
        call mpi_send(h_xnobc,nobs_tot,mpi_real4,0,1,mpi_comm_world,ierr)
       end if ! if nanal == 0
    end if ! nproc <= nanals
    call mpi_bcast(ensmean_ob,nobs_tot,mpi_real4,0,mpi_comm_world,ierr)
    call mpi_bcast(sprd_ob,nobs_tot,mpi_real4,0,mpi_comm_world,ierr)
    deallocate(h_xnobc)

 end subroutine mpi_getobs

end module mpi_readobs
