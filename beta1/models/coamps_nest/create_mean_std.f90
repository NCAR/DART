! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

program create_mean_std

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
                                   do_output,                    &
                                   get_unit

  use mpi_utilities_mod,    only : initialize_mpi_utilities,     &
                                   finalize_mpi_utilities,       &
                                   my_task_id,                   &
                                   task_count

  use types_mod,            only : r4, r8

  implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

  character(len=*), parameter :: routine = 'create_mean_std'
  character(len=64)           :: coamps_file_name

  integer, parameter :: fld_max=3000

  character(len=16), parameter :: nml_file = 'compute_mean.nml'
  integer                      :: kka      = 30
  integer, dimension(7)        :: m        = 101
  integer, dimension(7)        :: n        = 101
  integer                      :: nfiles   = 0
  integer                      :: ens_size = 96
  character(len=10)            :: cdtg
  character(len=64)            :: flist(fld_max)
  character(len=80)            :: dsnrff='./'
  namelist /compute_mean/ kka, m, n
  namelist /field_proc/ ens_size, nfiles, cdtg, dsnrff, flist

  real(kind=r8), allocatable  :: xdata(:,:)
  real(kind=r8), allocatable  :: xmean(:)
  real(kind=r8), allocatable  :: xsdev(:)

  integer            :: nz
  integer            :: cur_file_index
  integer            :: fcnt

  integer            :: file_unit
  integer            :: nml_unit

  integer            :: io_status
  character(len=5)   :: fileaction
  character(len=7)   :: filestatus
  character(len=180) :: dsnrff1

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
      read(nml_unit, nml=compute_mean, iostat=io_status)
      call check_io_status(io_status, routine, source, revision, revdate, &
                          'Reading ' // nml_file // ' compute_mean')
      read(nml_unit, nml=field_proc,   iostat=io_status)
      call check_io_status(io_status, routine, source, revision, revdate, &
                          'Reading ' // nml_file // ' field_proc')
    close(nml_unit)
  else
    call error_handler(E_ERR, 'create_mean_std', 'compute_mean name' // &
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
    !if(my_task_id==proc) print *,file_unit,coamps_file_name

    fldname=coamps_file_name(1:6)
    if(fldname == 'terrht') cycle FILE_LIST
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
    allocate(xdata(nsize,ens_size),stat=alloc_status)
      call check_alloc_status(alloc_status, routine, source, revision, &
                              revdate, 'allocate xdata')
    allocate(xmean(nsize),stat=alloc_status)
      call check_alloc_status(alloc_status, routine, source, revision, &
                              revdate, 'allocate xmean')
    allocate(xsdev(nsize),stat=alloc_status)
      call check_alloc_status(alloc_status, routine, source, revision, &
                              revdate, 'allocate xsdev')

    xmean(:) = 0.0
    xsdev(:) = 0.0
    do cur_file_index=1,ens_size
      write(dsnrff1,'(A,I5.5,A)') trim(dsnrff),cur_file_index,'/data/'
      if(.not. file_exist(trim(dsnrff1)//coamps_file_name)) then 
         print '(A)',trim(dsnrff1)//coamps_file_name
         deallocate(xdata,stat=dealloc_status)
           call check_dealloc_status(dealloc_status, routine, source, revision, &
                                     revdate, 'deallocate xdata')
         deallocate(xmean,stat=dealloc_status)
           call check_dealloc_status(dealloc_status, routine, source, revision, &
                                     revdate, 'deallocate xmean')
         deallocate(xsdev,stat=dealloc_status)
           call check_dealloc_status(dealloc_status, routine, source, revision, &
                                     revdate, 'deallocate xsdev')
         
         cycle FILE_LIST
      end if

      if(do_output())  print '(A)',trim(dsnrff1)//coamps_file_name

      open( unit=file_unit,                                         &
            file=trim(dsnrff1)//coamps_file_name,                   &
            status=filestatus, access='direct', action=fileaction,  &
            form='unformatted', recl=filesize, iostat=io_status)
      call check_io_status(io_status, routine, source, revision,    &
                           revdate, 'Opening ' //                   &
                           coamps_file_name // ' for read')
      call read_flat_file(file_unit,xdata(:,cur_file_index))
      do ii=1,nsize
        xmean(ii) = xmean(ii) + (xdata(ii,cur_file_index)   )/real(ens_size)
      end do
      close(file_unit)
    end do

    do ii=1,nsize
      do cur_file_index=1,ens_size
        xsdev(ii) = xsdev(ii) + ((xdata(ii,cur_file_index) - xmean(ii))**2)/real(ens_size-1)
      end do
    end do

    xsdev(:) = sqrt(xsdev(:))

    fileaction = 'write'
    filestatus = 'unknown'

    write(dsnrff1,'(A,A)') trim(dsnrff),'mean/data/'
    open( unit=file_unit,                                         &
          file=trim(dsnrff1)//coamps_file_name,                   &
          status=filestatus, access='direct', action=fileaction,  &
          form='unformatted', recl=filesize, iostat=io_status)
    call check_io_status(io_status, routine, source, revision,    &
                         revdate, 'Opening ' //                   &
                         coamps_file_name // ' for write')
    call write_flat_file(file_unit,xmean)
    close(file_unit)

    write(dsnrff1,'(A,A)') trim(dsnrff),'sdev/data/'
    open( unit=file_unit,                                         &
          file=trim(dsnrff1)//coamps_file_name,                   &
          status=filestatus, access='direct', action=fileaction,  &
          form='unformatted', recl=filesize, iostat=io_status)
    call write_flat_file(file_unit,xsdev)
    close(file_unit)

    deallocate(xdata,stat=dealloc_status)
      call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'deallocate xdata')
    deallocate(xmean,stat=dealloc_status)
      call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'deallocate xmean')
    deallocate(xsdev,stat=dealloc_status)
      call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'deallocate xsdev')

  end do FILE_LIST

100 continue
  call finalize_mpi_utilities()

end program create_mean_std

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
