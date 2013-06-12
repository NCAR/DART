! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

program create_increment

! This program pulls pieces out of the large COAMPS restart file,
! then assembles them into a state vector that can be used by DART.
! This includes two pieces of information - the current time and
! the actual state

  use coamps_util_mod,      only : check_alloc_status,           &
                                   check_dealloc_status,         &
                                   check_io_status,              &
                                   read_flat_file,               &
                                   write_flat_file

  use utilities_mod,        only : E_ERR,                        &
                                   error_handler,                &
                                   file_exist,                   &
                                   get_unit

  use mpi_utilities_mod,    only : initialize_mpi_utilities,     &
                                   finalize_mpi_utilities,       &
                                   my_task_id,                   &
                                   task_count

  use types_mod,            only : r8, r4

  implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

  character(len=*), parameter :: routine = 'create_increment'
  character(len=64)           :: coamps_file_name
  character(len=64)           :: coamps_bg_file_name
  character(len=64)           :: coamps_out_file_name

  integer, parameter :: fld_max=3000

  character(len=16), parameter :: nml_file = 'compute_incr.nml'
  integer                      :: kka      = 30
  integer, dimension(7)        :: m        = 101
  integer, dimension(7)        :: n        = 101
  integer                      :: nfiles   = 0
  integer                      :: ens_size = 96
  integer                      :: icycle   = 6
  character(len=10)            :: cdtg
  character(len=10)            :: cdtgm1
  character(len=64)            :: flist(fld_max)
  character(len=80)            :: dsnrff='./'
  namelist /compute_incr/ kka, m, n
  namelist /field_proc/ ens_size, nfiles, cdtg, cdtgm1, icycle, dsnrff, flist

  real(kind=r8), allocatable  :: xdata(:)
  real(kind=r8), allocatable  :: xmean_post(:)
  real(kind=r8), allocatable  :: xmean_prior(:)

  integer            :: nz
  integer            :: cur_file_index
  integer            :: fcnt

  integer            :: file_unit
  integer            :: nml_unit

  integer            :: io_status
  character(len=4)   :: icycle_str
  character(len=5)   :: fileaction
  character(len=7)   :: filestatus
  character(len=180) :: dsnrff_anal
  character(len=180) :: dsnrff_fcst
  character(len=180) :: dsnrff_incr

  integer :: fieldsize, nsize
  integer :: filesize
  integer :: r4_length

  character(len=6) :: fldname
  character(len=3) :: lvlname
  integer          :: nestnum, nprocs
  integer :: nfile_per_proc, nfile_lo, ishift, eshift, ifile, efile
  integer, parameter :: proc = 3
  integer :: ii

  integer :: alloc_status, dealloc_status

  call initialize_mpi_utilities()
  nprocs = task_count()

  nml_unit = get_unit()
  if (file_exist(nml_file)) then
    open(nml_unit, file=nml_file, status='old', iostat=io_status)
    call check_io_status(io_status, routine, source, revision, revdate,   &
                         'Opening ' // nml_file)
      read(nml_unit, nml=compute_incr, iostat=io_status)
      call check_io_status(io_status, routine, source, revision, revdate, &
                          'Reading ' // nml_file // ' compute_incr')
      read(nml_unit, nml=field_proc,   iostat=io_status)
      call check_io_status(io_status, routine, source, revision, revdate, &
                          'Reading ' // nml_file // ' field_proc')
    close(nml_unit)
  else
    call error_handler(E_ERR, 'create_increment', 'compute_increment name' // &
                      'list read failed - target file not found',       &
                       source, revision, revdate)
  end if


  nfile_per_proc = int(nfiles/nprocs)
  nfile_lo       = mod(nfiles,nprocs)

  if(my_task_id().lt.nfile_lo) then
    ishift = my_task_id()
    eshift = 0
  else
    ishift = nfile_lo
    eshift = -1
  end if

  ifile = my_task_id()*(nfile_per_proc) + 1 + ishift
  efile = ifile + nfile_per_proc + eshift

  FILE_LIST: do fcnt=ifile,efile
    file_unit = 700+my_task_id()
    
    coamps_file_name = flist(fcnt)


    coamps_bg_file_name=coamps_file_name
    write(icycle_str,'(I4.4)') icycle
    coamps_bg_file_name(38:52) = cdtgm1//'_'//icycle_str
    coamps_bg_file_name(58:64) = 'fcstfld'

    coamps_out_file_name = coamps_file_name
    coamps_out_file_name(58:64) = 'fcstfld'

    fldname=coamps_file_name(1:6)
    if(fldname == 'terrht') cycle
    lvlname=coamps_file_name(8:10)
    read(coamps_file_name(26:26),'(I1)') nestnum
    fieldsize=m(nestnum)*n(nestnum)


    select case(lvlname)
      case ('sig')
        select case(fldname)
          case ('wwwind')
            nz = kka+1
          case default
            nz = kka
        end select
      case default
        nz = 1
    end select

    fileaction = 'read'
    filestatus = 'old'

    inquire(IOLENGTH=r4_length) 0_r4
    filesize=r4_length*nz*fieldsize

    nsize=fieldsize*nz
    allocate(xdata(nsize),stat=alloc_status)
      call check_alloc_status(alloc_status, routine, source, revision, &
                              revdate, 'allocate xdata')
    allocate(xmean_post(nsize),stat=alloc_status)
      call check_alloc_status(alloc_status, routine, source, revision, &
                              revdate, 'allocate xdata_post')
    allocate(xmean_prior(nsize),stat=alloc_status)
      call check_alloc_status(alloc_status, routine, source, revision, &
                              revdate, 'allocate xdata_prior')

    xmean_post(:)  = 0.0
    xmean_prior(:) = 0.0
    do cur_file_index=1,ens_size
      !write(dsnrff_anal,'(A,I5.5,A)') trim(dsnrff),cur_file_index,'/data/analyses/'
      write(dsnrff_anal,'(A,I5.5,A)') trim(dsnrff),cur_file_index,'/data/'
      write(dsnrff_fcst,'(A,I5.5,A)') trim(dsnrff),cur_file_index,'/data/'
      
      if(.not. file_exist(trim(dsnrff_anal)//coamps_file_name)) print '(A)',trim(dsnrff_anal)//coamps_file_name
      if(.not. file_exist(trim(dsnrff_fcst)//coamps_bg_file_name)) print '(A)',trim(dsnrff_fcst)//coamps_bg_file_name

      if(.not. file_exist(trim(dsnrff_fcst)//coamps_bg_file_name) .or. &
         .not. file_exist(trim(dsnrff_anal)//coamps_file_name)) then
           deallocate(xdata,stat=dealloc_status)
             call check_dealloc_status(dealloc_status, routine, source, revision, &
                                       revdate, 'deallocate xdata_p')
           deallocate(xmean_post,stat=dealloc_status)
             call check_dealloc_status(dealloc_status, routine, source, revision, &
                                       revdate, 'deallocate xdata_post')
           deallocate(xmean_prior,stat=dealloc_status)
             call check_dealloc_status(dealloc_status, routine, source, revision, &
                                       revdate, 'deallocate xdata_prior')

            cycle FILE_LIST
      end if

      open( unit=file_unit,                                         &
            file=trim(dsnrff_anal)//coamps_file_name,               &
            status=filestatus, access='direct', action=fileaction,  &
            form='unformatted', recl=filesize, iostat=io_status)
      call check_io_status(io_status, routine, source, revision,    &
                           revdate, 'Opening ' //                   &
                           coamps_file_name // ' for read')
      call read_flat_file(file_unit,xdata(:))
      do ii=1,nsize
        xmean_post(ii) = xmean_post(ii) + (xdata(ii))/real(ens_size,kind=r8)
      end do
      close(file_unit)

      open( unit=file_unit,                                         &
            file=trim(dsnrff_fcst)//coamps_bg_file_name,            &
            status=filestatus, access='direct', action=fileaction,  &
            form='unformatted', recl=filesize, iostat=io_status)
      call check_io_status(io_status, routine, source, revision,    &
                           revdate, 'Opening ' //                   &
                           coamps_bg_file_name // ' for read')
      call read_flat_file(file_unit,xdata(:))
      do ii=1,nsize
        xmean_prior(ii) = xmean_prior(ii) + (xdata(ii))/real(ens_size,kind=r8)
      end do
      close(file_unit)
    end do

    do ii=1,nsize
      xdata(ii) = xmean_post(ii) - xmean_prior(ii)
    end do

    fileaction = 'write'
    filestatus = 'unknown'

    write(dsnrff_incr,'(A,A)') trim(dsnrff),'incr/data/'
    print '(A)','Writing '//trim(dsnrff_incr)//coamps_out_file_name
    open( unit=file_unit,                                         &
          file=trim(dsnrff_incr)//coamps_out_file_name,               &
          status=filestatus, access='direct', action=fileaction,  &
          form='unformatted', recl=filesize, iostat=io_status)
    call check_io_status(io_status, routine, source, revision,    &
                         revdate, 'Opening ' //                   &
                         coamps_out_file_name // ' for write')
    call write_flat_file(file_unit,xdata)
    close(file_unit)

    deallocate(xdata,stat=dealloc_status)
      call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'deallocate xdata')
    deallocate(xmean_post,stat=dealloc_status)
      call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'deallocate xdata_post')
    deallocate(xmean_prior,stat=dealloc_status)
      call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'deallocate xdata_prior')

  end do FILE_LIST

100 continue
  call finalize_mpi_utilities()

end program create_increment

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
