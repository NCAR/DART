program gsi_to_dart

use kinds
use mpisetup
use params, only : datestring, datapath, convert_conv, lie_about_ob_times, recenter_about_mean_prior, &
                   convert_sat, sattypes_rad, dsis, sattypes_oz,dsis, ens_size, obs_seq_out_filename, &
                   nsats_rad,nsats_oz,nsatmax_rad,nsatmax_oz, modify_dart_qc_flag_for_big_ob_error, &
                   obsprd_prior, ensmean_obnobc, ensmean_ob, ob, oberrvar, oberrvar_orig, &
                   obloclon, obloclat, obpress, obtime, biaspreds, anal_ob, stattype, indxsat, &
                   obtype, obsmod_cleanup, write_FO_for_these_obs_types, write_prior_copies, &
                   exclude_these_obs_types, output_option, anal_ob_chunk, variance_coef

use           radinfo, only : radinfo_read, radinfo_clean
use       mpi_readobs, only : mpi_getobs
use  dart_obs_seq_mod, only : dart_obs_seq, set_debug
use     utilities_mod, only : find_namelist_in_file, check_namelist_read
use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities

implicit none

logical         :: debug = .false.
integer(i_kind) :: nobs_sat, nobs_oz, nobs_conv, nobstot
integer(i_kind) :: nobs_start, nobs_end
integer(i_kind) :: i, irank
integer(i_kind) :: ierr
integer(i_kind) :: unitnml, io
integer(i_kind) :: pe_write_conv,pe_write_rad
integer(i_kind),allocatable,dimension(:) :: ista, iend, displs, scount
character(len=4) :: char_proc
character(len=256) :: my_output_filename
character(len=256) :: obs_seq_out_filename_conv, obs_seq_out_filename_sat
real(r_single),    allocatable, dimension(:) :: workgrid_in, workgrid_out

! namelist variables are declared and initialized in params and radinfo
namelist /gsi_to_dart_nml/ ens_size, &
   convert_conv, convert_sat, datestring, datapath, sattypes_rad, dsis, sattypes_oz, &
   write_FO_for_these_obs_types, write_prior_copies, &
   exclude_these_obs_types, output_option, obs_seq_out_filename, &
   lie_about_ob_times, recenter_about_mean_prior, debug, &
   modify_dart_qc_flag_for_big_ob_error, variance_coef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize MPI (from mpisetup module)
! provides nproc           process number
! provides numproc         total number of processes
! provides mpi_comm_world
! provides mpi_status
! provides mpi_real4   NOTE: many variables are declared 'r_single', which requires
!                      the use of mpi_real4, even when the CPP directive is _REAL8_.

call mpi_initialize

! Print out some info
call initialize_mpi_utilities('gsi_to_dart',communicator=mpi_comm_world)

! The barrier makes sure the 'starting' banner from initialize_mpi_utilities is
! printed before any subsequent print statements. Some tasks were printing the
! error messages before the starting banner ... 
call mpi_barrier(mpi_comm_world,ierr)

! Read namelist on all PEs
call find_namelist_in_file("input.nml", "gsi_to_dart_nml", unitnml)
read(unitnml, nml = gsi_to_dart_nml, iostat = io)
call check_namelist_read(unitnml, io, "gsi_to_dart_nml")

! Pass the debug value to the dart_obs_seq_mod
call set_debug(debug)

! Do some error checking
if (ens_size < 2) then
   print *,'ERROR: ens_size ',ens_size,' must be > 1'
   call stop2(18)
endif

! Do some error checking
if (numproc < ens_size+1) then
   print *,'total number of mpi tasks must be >= ens_size+1'
   print *,'tasks, ens_size+1 = ',numproc,ens_size+1
   call stop2(19)
endif

if ( output_option .gt. 3 .or. output_option .le. 0) then
   write(*,*)' output_option must be 1, 2, or 3, not ',output_option
   call stop2(20)
endif

if ( .not. convert_sat .and. .not. convert_conv ) then
   write(*,*)' Both convert_sat and convert_conv are false. This program will not do anything, so exiting here.'
   call stop2(21)
endif

! assign processes for writing files
pe_write_conv=0
pe_write_rad=max(0,min(1,numproc-1))

! find number of satellite files
nsats_rad=0
if ( convert_sat ) then
   do i=1,nsatmax_rad
     if(sattypes_rad(i) == ' ') cycle
     nsats_rad=nsats_rad+1
   end do
   ! initialize satellite stuff
   call radinfo_read
endif
if(nproc == 0) write(*,*) 'number of satellite radiance files used ',nsats_rad

! find number of ozone files
nsats_oz=0
do i=1,nsatmax_oz
  if(sattypes_oz(i) == ' ') cycle
  nsats_oz=nsats_oz+1
end do
if(nproc == 0) write(*,*) 'number of satellite ozone files used ',nsats_oz

! read the diag files to get obs, ob priors, and observation errors/location/time
call mpi_getobs(datapath, datestring, nobs_conv, nobs_oz, nobs_sat, nobstot, &
                obsprd_prior, ensmean_obnobc, ensmean_ob, ob, &
                oberrvar, obloclon, obloclat, obpress, &
                obtime, oberrvar_orig, stattype, obtype, biaspreds,&
                anal_ob,indxsat,ens_size)  ! anal_ob only allocated on nproc==0
                                           ! everything else on all PEs

if(nproc == 0) write(*,*) 'total number of obs ',nobstot

! output obs_sequence file
if ( output_option .eq. 3 ) then
   if ( .not. allocated(anal_ob)) allocate(anal_ob(ens_size,nobstot)) ! may have been allocated in mpi_getobs
   call mpi_bcast(anal_ob,nobstot*ens_size,mpi_real4,0,mpi_comm_world,ierr) !  should really do mpi_scatterv...but doing so doesn't seem to save memory overall...

   if ( convert_conv ) then
      nobs_start = 1
      nobs_end = nobs_conv
      if ( convert_sat ) then
         nobs_end = nobstot
      end if
   else
      if ( convert_sat ) then
         nobs_start = nobs_conv+nobs_oz+1
         nobs_end = nobstot
      end if
   end if

   allocate(ista(0:numproc-1), iend(0:numproc-1))
!  allocate(displs(0:numproc-1), scount(0:numproc-1))
   do irank = 0,numproc-1
!hcl      call para_range(1, nobstot, numproc, irank, ista(irank), iend(irank))
     call para_range(nobs_start, nobs_end, numproc, irank, ista(irank), iend(irank))
!     scount(irank) = iend(irank)-ista(irank) + 1  ! number of obs on this PE
   enddo
   write(*,*)'The processor ', nproc, ' will write obs from ', ista(nproc), ' to ', iend(nproc)
!  displs(0) = 0
!  do i = 1, numproc-1
!     displs(i) = scount(i-1) + displs(i-1)
!  end do

!  allocate( anal_ob_chunk(ens_size,ista(nproc):iend(nproc)) )
!  if ( nproc.eq.0 ) allocate( workgrid_in(nobstot))
!  allocate( workgrid_out(scount(nproc)) )

!  do i = 1,ens_size
!     if ( nproc .eq. 0 ) workgrid_in(:) = anal_ob(i,:)
!     call mpi_scatterv(workgrid_in(:),scount,displs, &
!             & mpi_real4,workgrid_out(:),scount(nproc),mpi_real4,0,&
!             & mpi_comm_world,ierr)
!     anal_ob_chunk(i,:) = workgrid_out(:)
!  enddo

   write(char_proc,'(i4.4)') nproc
   my_output_filename = trim(adjustl(obs_seq_out_filename))//'.'//char_proc
   write(*,*)'About to write ',trim(my_output_filename),' from processor',nproc
   if ( ista(nproc).le.iend(nproc)) then  ! Make sure the starting ob to process <= ending ob to process (probably should never happen)
      call dart_obs_seq (datestring,            &
         nobs_conv, nobs_oz, nobs_sat, nobstot, &
         ens_size, my_output_filename,          &
         ista(nproc),iend(nproc) ) 
   endif
   deallocate(ista,iend)
!  deallocate(displs,scount)
!  deallocate(anal_ob_chunk,workgrid_out)
!  if ( nproc.eq.0 ) deallocate(workgrid_in)
else if ( output_option .eq. 2 ) then
   ! if writing radiance file, send anal_ob to pe_write_rad
   if ( convert_sat .and. nproc == pe_write_rad) then
      if ( .not. allocated(anal_ob)) allocate(anal_ob(ens_size,nobstot)) ! may have been allocated in mpi_getobs
   endif
   if ( convert_sat .and. nproc==0)            call mpi_send(anal_ob,nobstot*ens_size,mpi_real4,pe_write_rad,1,mpi_comm_world,ierr) ! send from nproc=0
   if ( convert_sat .and. nproc==pe_write_rad) call mpi_recv(anal_ob,nobstot*ens_size,mpi_real4,0,1,mpi_comm_world,mpi_status,ierr) ! recieve on nproc=pe_write_rad
   ! what if pe_write_rad == 0?  will the mpi_send/recv fail (sending and receiving from same processor)?

   obs_seq_out_filename_conv = trim(adjustl(obs_seq_out_filename))//'.conv'
   obs_seq_out_filename_sat  = trim(adjustl(obs_seq_out_filename))//'.rad'

   ! write out DART obs_sequence files
   if ( convert_conv .and. nobs_conv > 0 ) then
      if ( nproc == pe_write_conv ) then
         write(*,*)'About to write "'//trim(obs_seq_out_filename_conv)//'"'
         call dart_obs_seq (datestring,                            &
                            nobs_conv, nobs_oz, nobs_sat, nobstot, &
                            ens_size, obs_seq_out_filename_conv,     &
                            1, nobs_conv) 
      end if
   endif

   if ( convert_sat .and. nobs_sat > 0 ) then
      if ( nproc == pe_write_rad ) then
         write(*,*)'About to write "'//trim(obs_seq_out_filename_sat)//'"'
         call dart_obs_seq (datestring,                            &
                            nobs_conv, nobs_oz, nobs_sat, nobstot, &
                            ens_size, obs_seq_out_filename_sat,      &
                            nobs_conv+nobs_oz+1, nobstot) 
      end if
   end if
else if ( output_option .eq. 1 ) then
   ! write out DART obs_sequence files
   if ( nproc == pe_write_conv ) then
      write(*,*)'About to write "'//trim(obs_seq_out_filename)//'"'
      call dart_obs_seq (datestring,                            &
                         nobs_conv, nobs_oz, nobs_sat, nobstot, &
                         ens_size, obs_seq_out_filename,     &
                         1, nobstot) 
   end if
endif

! deallocate arrays...these were allocated in mpi_getobs 
call obsmod_cleanup ! from params
if ( convert_sat ) call radinfo_clean  ! from radinfo

! print ending info
call finalize_mpi_utilities(callfinalize=.false.)

! finalize MPI
call mpi_cleanup

end program gsi_to_dart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine stop2(ierror_code)
! adapted from GSI/src/main/stop1.f90

  use kinds, only: i_kind
  use mpisetup, only : mpi_comm_world, MPI_SUCCESS
  implicit none

  integer(i_kind), intent(in) :: ierror_code
  integer(i_kind)             :: ierr

  write(*,*)'****STOP2****  ABORTING EXECUTION w/code=',ierror_code
  write(0,*)'****STOP2****  ABORTING EXECUTION w/code=',ierror_code
  call mpi_abort(mpi_comm_world,ierror_code,ierr)
  stop
  return
end subroutine stop2

subroutine para_range(n1, n2, nprocs, irank, ista, iend) 

  use kinds, only: i_kind
  implicit none

  integer(i_kind), intent(in)  :: n1, n2, nprocs, irank
  integer(i_kind), intent(out) :: ista, iend

  integer(i_kind) :: iwork1, iwork2

  iwork1 = (n2 - n1 + 1) / nprocs
  iwork2 = mod(n2 - n1 + 1, nprocs)
  ista = irank * iwork1 + n1 + min(irank, iwork2)
  iend = ista + iwork1 - 1
  if (iwork2 > irank) iend = iend + 1
  return
end subroutine para_range 

