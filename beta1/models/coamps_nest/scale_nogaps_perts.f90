! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

program scale_nogaps_perts

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

  character(len=*), parameter :: routine = 'scale_nogaps_perts'
  character(len=64)           :: nogaps_file_name

  integer, parameter :: fld_max=6000

  character(len=*), parameter  :: nml_file = 'scale_nogaps_perts.nml'
  integer                      :: nfiles   = 1
  integer                      :: ens_size = 2
  character(len=64)            :: flist(fld_max)
  character(len=80)            :: dsngff='./'
  character(len=80)            :: dsngff_out='./scaled50'
  real(kind=r8)                :: pert_scale=0.5
  namelist /field_proc/ ens_size, nfiles, dsngff, flist, pert_scale, dsngff_out

  real(kind=r8), allocatable  :: xdata(:,:)
  real(kind=r8), allocatable  :: xmean(:)

  integer            :: nx
  integer            :: ny
  integer            :: nz
  integer            :: cur_file_index
  integer            :: fcnt

  integer            :: file_unit
  integer            :: nml_unit

  integer            :: io_status
  character(len=5)   :: fileaction
  character(len=7)   :: filestatus
  character(len=180) :: dsngff1

  integer :: fieldsize, nsize
  integer :: filesize
  integer :: r4_length

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
    call error_handler(E_ERR, 'scale_nogaps_perts', 'compute_mean name' // &
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

  FILE_LOOP: do fcnt=ifile,efile
    file_unit = 700+my_task_id()
    
    nogaps_file_name = flist(fcnt)

    fldname=nogaps_file_name(1:6)
    lvlname=nogaps_file_name(8:10)
    read(nogaps_file_name(30:32),'(I3)') nx
    read(nogaps_file_name(34:36),'(I3)') ny
    fieldsize=nx*ny
    nz = 1

    !if(my_task_id==proc) print '(3(I3,1x),A)',file_unit,nx,ny,nogaps_file_name

    fileaction = 'read'
    filestatus = 'old'

    inquire(IOLENGTH=r4_length) 0_r4
    filesize=r4_length*nz*fieldsize

    nsize=fieldsize*nz
    allocate(xdata(nsize,ens_size))
    allocate(xmean(nsize))

    xmean(:) = 0.0
    if(my_task_id()==proc) print *,fcnt,nogaps_file_name
    READ_LOOP: do cur_file_index=1,ens_size
      write(dsngff1,'(A,A,I3.3,A)') trim(dsngff),'mbr',cur_file_index,'/outp/'

      if(.not. file_exist(trim(dsngff1)//nogaps_file_name) ) then
        print '(2A)', 'missing: ',trim(dsngff1)//nogaps_file_name
        deallocate(xdata)
        deallocate(xmean)
        cycle FILE_LOOP
      end if

      open( unit=file_unit,                                         &
            file=trim(dsngff1)//nogaps_file_name,                   &
            status=filestatus, access='direct', action=fileaction,  &
            form='unformatted', recl=filesize, iostat=io_status)
      call check_io_status(io_status, routine, source, revision,    &
                           revdate, 'Opening ' //                   &
                           nogaps_file_name // ' for read')
      call read_flat_file(file_unit,xdata(:,cur_file_index))
      do ii=1,nsize
        xmean(ii) = xmean(ii) + (xdata(ii,cur_file_index)   )/real(ens_size)
      end do
      close(file_unit)
    end do READ_LOOP

    fileaction = 'write'
    filestatus = 'unknown'

    WRITE_LOOP: do cur_file_index=1,ens_size

      xdata(:, cur_file_index) = pert_scale * (xdata(:, cur_file_index) - xmean(:) ) + xmean(:)

      write(dsngff1,'(A,A,I3.3,A)') trim(dsngff_out),'mbr',cur_file_index,'/outp/'
      !print '(A)',trim(dsngff1)//nogaps_file_name
      open( unit=file_unit,                                         &
            file=trim(dsngff1)//nogaps_file_name,                   &
            status=filestatus, access='direct', action=fileaction,  &
            form='unformatted', recl=filesize, iostat=io_status)
      call check_io_status(io_status, routine, source, revision,    &
                           revdate, 'Opening ' //                   &
                           nogaps_file_name // ' for write')
      call write_flat_file(file_unit,xdata(:, cur_file_index))
      close(file_unit)
    end do WRITE_LOOP

    deallocate(xdata)
    deallocate(xmean)

  end do FILE_LOOP

100 continue
  call finalize_mpi_utilities()

end program scale_nogaps_perts

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
