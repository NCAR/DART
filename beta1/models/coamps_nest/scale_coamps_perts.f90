! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

program scale_coamps_perts

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

  character(len=*), parameter :: routine = 'scale_coamps_perts'
  character(len=64)           :: coamps_file_name

  integer, parameter :: fld_max=6000

  character(len=*), parameter  :: nml_file = 'scale_coamps_perts.nml'
  integer                      :: kka      = 30
  integer, dimension(7)        :: m        = 101
  integer, dimension(7)        :: n        = 101
  integer                      :: nfiles   = 1
  integer, dimension(10)       :: bad_mems = -1
  integer                      :: ens_size = 96
  character(len=64)            :: flist(fld_max)
  character(len=80)            :: dsnrff='./'
  character(len=80)            :: dsnrff_out='./scaled50'
  real(kind=r8)                :: pert_scale=0.5
  namelist /domain_dims/ kka, m, n
  namelist /field_proc/ ens_size, nfiles, dsnrff, flist, pert_scale, dsnrff_out, &
                        bad_mems

  real(kind=r8), allocatable  :: xdata(:,:)
  real(kind=r8), allocatable  :: xmean(:)

  logical, allocatable  :: mem_available_mask(:)

  integer            :: nmems
  integer            :: nestnum
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
  character(len=180) :: dsnrff1

  integer :: fieldsize, nsize
  integer :: filesize
  integer :: r4_length

  real(kind=r8) :: xmax, xmin

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
    call error_handler(E_ERR, 'scale_coamps_perts', 'compute_mean name' // &
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
    
    coamps_file_name = flist(fcnt)

    fldname=coamps_file_name(1:6)
    if(fldname == 'terrht') cycle FILE_LOOP
    fldtype=coamps_file_name(58:64)
    lvlname=coamps_file_name(8:10)
    read(coamps_file_name(26:26),'(I1)') nestnum

    fieldsize=m(nestnum)*n(nestnum)

    select case(fldtype)
     case ('bndyfld','bndyprt')

     case default
       select case(lvlname)
        case ('sig')
          select case(fldname)
            case ('wwwind','bdwwnd')
              nz = kka+1
            case default
              nz = kka
          end select
        case default
          nz = 1
       end select
    end select

    !if(my_task_id==proc) print '(3(I3,1x),A)',file_unit,nx,ny,coamps_file_name

    fileaction = 'read'
    filestatus = 'old'

    inquire(IOLENGTH=r4_length) 0_r4
    filesize=r4_length*nz*fieldsize

    nsize=fieldsize*nz
    allocate(xdata(nsize,ens_size))
    allocate(xmean(nsize))
!    allocate(xvar(nsize))
    allocate(mem_available_mask(ens_size))
    mem_available_mask(:) = .true.

    xmean(:) = 0.0
    if(my_task_id()==proc) print *,fcnt,coamps_file_name
    nmems = ens_size
    READ_LOOP: do cur_file_index=1,ens_size
      write(dsnrff1,'(A,I5.5,A)') trim(dsnrff),cur_file_index,'/data/analyses/'

      if(.not. file_exist(trim(dsnrff1)//coamps_file_name) ) then
        !print '(2A)', 'missing: ',trim(dsnrff1)//coamps_file_name
        nmems = nmems - 1
        mem_available_mask(cur_file_index) = .false.
        cycle READ_LOOP
      end if

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

      xmax=maxval(xdata(:, cur_file_index))
      xmin=minval(xdata(:, cur_file_index))
      if(cur_file_index==1) then
      print '(2(A,1x),I03,1x,A,2(1x,F10.3))' ,coamps_file_name,'Member',cur_file_index,'max/min:',xmax,xmin
      end if

      close(file_unit)
    end do READ_LOOP

    xmean(:) = xmean(:)/real(nmems, kind=r8)

!    xvar(:) = 0_r8
!    VAR_LOOP1: do cur_file_index=1,ens_size
!      if(.not. mem_available_mask(cur_file_index)) cycle VAR_LOOP1
!      xvar(:) = xvar(:) + (xdata(:, cur_file_index) - xmean(:) )**2
!    end do VAR_LOOP1
!    xvar(:) = sqrt(xvar(:)/real(nmems-1, kind=r8))
!    xmax=maxval(xvar(:))
!    xmin=minval(xvar(:))
!    print '(72x,A,2(1x,F10.3))' ,'var max/min:',xmax,xmin


    fileaction = 'write'
    filestatus = 'unknown'

    WRITE_LOOP: do cur_file_index=1,ens_size
      if(.not. mem_available_mask(cur_file_index)) cycle WRITE_LOOP

      xdata(:, cur_file_index) = pert_scale * (xdata(:, cur_file_index) - xmean(:) ) + xmean(:)

      xmax=maxval(xdata(:, cur_file_index))
      xmin=minval(xdata(:, cur_file_index))
      if(cur_file_index==1) then
      print '(2(A,1x),I03,1x,A,2(1x,F10.3))' ,coamps_file_name,'Member',cur_file_index,'max/min:',xmax,xmin
      end if

      write(dsnrff1,'(A,I5.5,A)') trim(dsnrff_out),cur_file_index,'/data/analyses/'
      open( unit=file_unit,                                         &
            file=trim(dsnrff1)//coamps_file_name,                   &
            status=filestatus, access='direct', action=fileaction,  &
            form='unformatted', recl=filesize, iostat=io_status)
      call check_io_status(io_status, routine, source, revision,    &
                           revdate, 'Opening ' //                   &
                           coamps_file_name // ' for write')
      call write_flat_file(file_unit,xdata(:, cur_file_index))
      close(file_unit)
    end do WRITE_LOOP

!    xvar(:) = 0_r8
!    VAR_LOOP2: do cur_file_index=1,ens_size
!      if(.not. mem_available_mask(cur_file_index)) cycle VAR_LOOP2
!      xvar(:) = xvar(:) + (xdata(:, cur_file_index) - xmean(:) )**2
!    end do VAR_LOOP2
!    xvar(:) = sqrt(xvar(:)/real(nmems-1, kind=r8))
!    xmax=maxval(xvar(:))
!    xmin=minval(xvar(:))
!    print '(72x,A,2(1x,F10.3))' ,'var max/min:',xmax,xmin

!    deallocate(xvar)
    deallocate(xdata)
    deallocate(xmean)
    deallocate(mem_available_mask)

  end do FILE_LOOP

100 continue
  call finalize_mpi_utilities()

end program scale_coamps_perts

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
