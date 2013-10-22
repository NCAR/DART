! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

program recntr_bndyperts

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

  character(len=*), parameter :: routine = 'recntr_bndyperts'
  character(len=64)           :: coamps_file_name
  character(len=64)           :: coamps_new_file_name

  integer, parameter :: fld_max=6000

  character(len=*), parameter  :: nml_file = 'recntr_bndyperts.nml'
  integer                      :: nfiles   = 1
  integer                      :: ens_size = 96
  character(len=64)            :: flist(fld_max)
  character(len=10)            :: dtg_new = '1970010100'
  character(len=80)            :: dsnrff='./'
  character(len=80)            :: dsnrff_out='./'
  character(len=80)            :: dsnrff_new='./detrm_bdy'
  namelist /field_proc/ ens_size, nfiles, dsnrff, flist, dsnrff_out, dsnrff_new, dtg_new

  real(kind=r8), allocatable  :: xdata(:,:)
  real(kind=r8), allocatable  :: xmean(:)
  real(kind=r8), allocatable  :: xmean_new(:)

  integer            :: nmems
  integer            :: nestnum
  integer            :: nz
  integer            :: cur_file_index
  integer            :: fcnt

  integer            :: file_unit
  integer            :: nml_unit

  integer            :: io_status
  character(len=5)   :: fileaction
  character(len=7)   :: filestatus
  character(len=10)  :: dtg
  character(len=180) :: dsnrff1

  integer :: fieldsize, nsize
  integer :: filesize
  integer :: r4_length

  character(len=7) :: fldtype
  character(len=6) :: fldname
  character(len=3) :: lvlname
  integer          :: nprocs
  integer :: nfile_per_proc, nfile_lo, ishift, eshift, ifile, efile
  integer, parameter :: proc = 3
  integer :: ii

  call initialize_mpi_utilities()
  nprocs = task_count()

  nml_unit = get_unit()
  if (file_exist(nml_file)) then
    open(nml_unit, file=nml_file, status='old', iostat=io_status)
    call check_io_status(io_status, routine, source, revision, revdate,   &
                        'Opening ' // nml_file)

    read(nml_unit, nml=field_proc,   iostat=io_status)
    call check_io_status(io_status, routine, source, revision, revdate, &
                        'Reading ' // nml_file // ' field_proc')
    close(nml_unit)
  else
    call error_handler(E_ERR, 'recntr_bndyperts', 'compute_mean name' // &
                      'list read failed - target file not found',       &
                       source, revision, revdate)
  end if

  nfile_per_proc = int(nfiles/nprocs)
  nfile_lo       = mod(nfiles,nprocs)
  print *,nfile_per_proc,nfile_lo

  if(my_task_id().lt.nfile_lo) then
    ishift = my_task_id()
    eshift = 0
  else
    ishift = nfile_lo
    eshift = -1
  end if

  ifile = my_task_id()*(nfile_per_proc) + 1 + ishift
  efile = ifile + nfile_per_proc + eshift

  FILE_LOOP: do fcnt=ifile,efile
    file_unit = 700+my_task_id()
    
    coamps_file_name = flist(fcnt)
    coamps_new_file_name = coamps_file_name
    coamps_new_file_name(38:47) = dtg_new

    fldname=coamps_file_name(1:6)
    if(fldname == 'terrht') cycle FILE_LOOP

    fldtype = coamps_file_name(58:64)
    lvlname = coamps_file_name(8:10)
    dtg     = coamps_file_name(38:47)

    read(coamps_file_name(26:26),'(I1)') nestnum
    read(coamps_file_name(28:31),'(I4)') fieldsize
    read(coamps_file_name(33:36),'(I4)') nz

    print '(2(I3,1x),2(A,1x))',fcnt,file_unit,coamps_file_name, coamps_new_file_name
    !if(my_task_id()==proc) print '(1(I3,1x),A)',file_unit,coamps_file_name

    fileaction = 'read'
    filestatus = 'old'

    inquire(IOLENGTH=r4_length) 0_r4
    filesize=r4_length*nz*fieldsize

    nsize=fieldsize*nz
    allocate(xdata(nsize,ens_size))
    allocate(xmean(nsize))
    allocate(xmean_new(nsize))

    xmean(:) = 0.0
    xmean_new(:) = 0.0
    if(my_task_id()==proc) print *,fcnt,coamps_file_name
    nmems = ens_size
    READ_LOOP: do cur_file_index=1,ens_size
      write(dsnrff1,'(A,I5.5,A)') trim(dsnrff),cur_file_index,'/data/'

      open( unit=file_unit,                                         &
            file=trim(dsnrff1)//coamps_file_name,                   &
            status=filestatus, access='direct', action=fileaction,  &
            form='unformatted', recl=filesize, iostat=io_status)
      call check_io_status(io_status, routine, source, revision,    &
                           revdate, 'Opening ' //                   &
                           coamps_file_name // ' for read')
      call read_flat_file(file_unit,xdata(:,cur_file_index))
      do ii=1,nsize
        xmean(ii) = xmean(ii) + (xdata(ii,cur_file_index))
      end do

      close(file_unit)
    end do READ_LOOP

    xmean(:) = xmean(:)/real(nmems, kind=r8)

    write(dsnrff1,'(2A)') trim(dsnrff_new),'/data/'
    open( unit=file_unit,                                         &
          file=trim(dsnrff1)//coamps_new_file_name,               &
          status=filestatus, access='direct', action=fileaction,  &
          form='unformatted', recl=filesize, iostat=io_status)
    call read_flat_file(file_unit,xmean_new(:))
    close(file_unit)

    fileaction = 'write'
    filestatus = 'unknown'

    WRITE_LOOP: do cur_file_index=1,ens_size

      xdata(:, cur_file_index) =  (xdata(:, cur_file_index) - xmean(:) ) + xmean_new(:)

      write(dsnrff1,'(A,I5.5,A)') trim(dsnrff_out),cur_file_index,'/data/'
      open( unit=file_unit,                                         &
            file=trim(dsnrff1)//coamps_new_file_name,               &
            status=filestatus, access='direct', action=fileaction,  &
            form='unformatted', recl=filesize, iostat=io_status)
      call check_io_status(io_status, routine, source, revision,    &
                           revdate, 'Opening ' //                   &
                           coamps_new_file_name // ' for write')
      call write_flat_file(file_unit,xdata(:, cur_file_index))
      close(file_unit)
    end do WRITE_LOOP

    deallocate(xdata)
    deallocate(xmean)
    deallocate(xmean_new)

  end do FILE_LOOP

100 continue
  call finalize_mpi_utilities()

end program recntr_bndyperts

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
