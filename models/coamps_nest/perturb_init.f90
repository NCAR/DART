! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

program perturb_init

! This program pulls pieces out of the large COAMPS restart file,
! then assembles them into a state vector that can be used by DART.
! This includes two pieces of information - the current time and
! the actual state

  use coamps_nest_mod,      only : coamps_nest,                  &
                                   initialize_nest,              &
                                   set_nest_id,                  &
                                   get_nest_i_width,             &
                                   get_nest_j_width

  use coamps_vertical_mod,  only : coamps_vertical,              &
                                   initialize_vertical,          &
                                   get_num_levels,               &
                                   get_msigma,                   &
                                   get_wsigma

  use coamps_util_mod,      only : check_alloc_status,           &
                                   check_dealloc_status,         &
                                   check_io_status,              &
                                   generate_flat_file_name,      &
                                   read_flat_file,               &
                                   write_flat_file,              &
                                   read_datahd_file,             &
                                   DATAHD_LEN

  use utilities_mod,        only : E_ERR,                        &
                                   error_handler,                &
                                   file_exist,                   &
                                   get_unit

  use random_seq_mod, only: random_seq_type, &
                            init_random_seq, &
                            random_gaussian, &
                            random_uniform

  use types_mod,            only : r8, r4

  implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

  character(len=*), parameter :: routine = 'perturb_init'
  character(len=64)           :: coamps_file_name
  character(len=64)           :: coamps_new_file_name

  character(len=64)           :: coamps_init_pert
  character(len=64)           :: coamps_init_file

  integer, parameter :: nflds = 5

  character(len=6), dimension(nflds), parameter :: &
  bndy_fldnames=(/'bduwnd', 'bdvwnd', 'bdpott', 'bdperp', 'bdmixr'/)

  character(len=6), dimension(nflds), parameter :: &
  full_fldnames=(/'uuwind', 'vvwind', 'pottmp', 'perprs', 'wvapor'/)

  character(len=*),  parameter :: nml_file = 'pert_init.nml'
  integer,           parameter :: maxgrds  = 7
  integer,           parameter :: NEST_ID  = 1
  integer,           parameter :: SECPERHR = 3600
  character(len=10), parameter :: dtg_pert = '0000000000'

  integer                      :: current_mem 
  integer                      :: xpts 
  integer                      :: ypts 
  integer                      :: gridnum 
  character(len=6)             :: fldname 

  real(kind=r8)                :: alpha       = 0.75_r8  
  real(kind=r8)                :: fcp_scale   = 1.75_r8  
  integer                      :: ens_size    = 10 
  integer                      :: npert_tot   = 400
  integer                      :: kgetbc      = 6 
  integer, dimension(3,maxgrds):: ktauf = &
  reshape( (/6,0,0,6,0,0,6,0,0,6,0,0,6,0,0,6,0,0,6,0,0/),(/3,maxgrds/))
  integer                      :: nbdypt      = 7
  character(len=10)            :: cdtg        = '1970010100'
  character(len=80)            :: dsnrff      ='./'
  character(len=80)            :: dsnrff_pert ='./'

  namelist /pert_init/ alpha, fcp_scale, ens_size, npert_tot, nbdypt, &
                       cdtg, dsnrff, dsnrff_pert, kgetbc, ktauf

  real(kind=r8), allocatable :: xpert(:)
  real(kind=r8), allocatable :: xpert_m1(:)
  real(kind=r8), allocatable :: xinit(:)
  integer,       allocatable :: ipert_nums(:, :)

  integer            :: file_unit
  integer            :: nml_unit

  integer            :: io_status
  character(len=5)   :: fileaction
  character(len=7)   :: filestatus

  integer :: fieldsize
  integer :: filesize
  integer :: r4_length


  real(kind=r8), dimension(DATAHD_LEN) :: coamps_datahd

  type(random_seq_type) :: random_sequence
  type(coamps_nest)     :: static_nest
  type(coamps_vertical) :: static_vgrid

  integer               :: rand_seq_seed
  integer               :: nx, ny, nz, nt
  integer               :: ii, nn, tt, ipert, tau_hr
  real(kind=r8)         :: ztop, zbot

!------------------------------------------------------------------------------
! initialize variables
!------------------------------------------------------------------------------

  inquire(IOLENGTH=r4_length) 0_r4
  file_unit = get_unit()
  nml_unit  = get_unit()

  print '(A)', 'Enter the ensemble member: '
  read(*, *) current_mem
!------------------------------------------------------------------------------
! Read namelist data
!------------------------------------------------------------------------------
  if (file_exist(nml_file)) then
    open(nml_unit, file=nml_file, status='old', iostat=io_status)
    call check_io_status(io_status, routine, source, revision, revdate,   &
                        'Opening ' // nml_file)

    read(nml_unit, nml=pert_init,   iostat=io_status)
    call check_io_status(io_status, routine, source, revision, revdate, &
                        'Reading ' // nml_file // ' pert_init')
    close(nml_unit)
  else
    call error_handler(E_ERR, 'recntr_initperts', 'compute_mean name' // &
                      'list read failed - target file not found',       &
                       source, revision, revdate)
  end if

!------------------------------------------------------------------------------
! Read data header file and set sizes
!------------------------------------------------------------------------------
  call read_datahd_file(cdtg, coamps_datahd)
  call initialize_vertical(coamps_datahd, static_vgrid)

  call set_nest_id(static_nest, NEST_ID)
  call initialize_nest(static_nest, cdtg, coamps_datahd)

  gridnum = 1
  !nt = ktauf(1, NEST_ID)/kgetbc + 1
  nt = 1
  nx = get_nest_i_width(static_nest)
  ny = get_nest_j_width(static_nest)
  !fieldsize = nx*nbdypt*2 + (ny - nbdypt*2) * nbdypt*2
  fieldsize = nx*ny
  
!------------------------------------------------------------------------------
! Initialize and set random sequence. Integer seed based on dtg.
!------------------------------------------------------------------------------
  read(cdtg(3:10),'(I8)') rand_seq_seed
  rand_seq_seed=-rand_seq_seed
  call init_random_seq(r=random_sequence, seed=rand_seq_seed)

  allocate(ipert_nums(ens_size, nt))
  do tt=1,nt
    do ii=1,ens_size
      !ipert_nums(ii, tt) = ceiling(random_uniform(random_sequence) * real(npert_tot, kind=r8))
      ipert_nums(ii, tt) = ii
    end do
  end do

!------------------------------------------------------------------------------
! Loop over all fields
!------------------------------------------------------------------------------
  VAR_LOOP : do nn=1,nflds
    fldname = full_fldnames(nn)

    print '(3(A,1x),I3.3)','Perturbing',fldname,'for member',current_mem

    select case(fldname)
    case('bdwwnd', 'wwwind')
      nz = get_num_levels(static_vgrid) + 1 
      ztop = get_wsigma(static_vgrid,  1)
      zbot = get_wsigma(static_vgrid, nz)
    case default
      nz = get_num_levels(static_vgrid)
      ztop = get_msigma(static_vgrid,  1)
      zbot = get_msigma(static_vgrid, nz)
    end select

    allocate(xinit(nz*fieldsize))
    allocate(xpert(nz*fieldsize))
    allocate(xpert_m1(nz*fieldsize))

    filesize=r4_length*nz*fieldsize

   TIME_LOOP : do tt=1,nt

    fileaction = 'read'
    filestatus = 'old'

    tau_hr = (tt-1)*kgetbc
!------------------------------------------------------------------------------
! Read all perts and compute mean.
!------------------------------------------------------------------------------
    xinit(:) = 0.0_r8
    xpert(:) = 0.0_r8

    xpts = nx ; ypts = ny

    PERT_LOOP : do ii=1,ens_size
      ipert = ipert_nums(ii, tt)

      call generate_flat_file_name( var_name   = fldname,      &
                                    level_type = 'sig',         &
                                    level1     = int(ztop),     &
                                    level2     = int(zbot),     &
                                    gridnum    = gridnum,       &
                                    aoflag     = 'a',           &
                                    xpts       = xpts,          &
                                    ypts       = ypts,          &
                                    dtg        = dtg_pert,      &
                                    tau_hh     = ipert,         &
                                    tau_mm     = 0,             &
                                    tau_ss     = 0,             &
                                    field_type = 'pertfld',     &
                                    file_name  = coamps_init_pert )

      open( unit=file_unit,                                         &
            file=trim(dsnrff_pert)//coamps_init_pert,               &
            status=filestatus, access='direct', action=fileaction,  &
            form='unformatted', recl=filesize, iostat=io_status)
      call check_io_status(io_status, routine, source, revision,    &
                           revdate, 'Opening ' //                   &
                           coamps_init_pert // ' for read')

      call read_flat_file(file_unit, xpert)
      close(file_unit)


      xinit(:) = xinit(:) + xpert(:)/real(ens_size, kind=r8)

    end do PERT_LOOP

!------------------------------------------------------------------------------
! Read perturbation for this member and recenter around zero mean.
!------------------------------------------------------------------------------
    ipert = ipert_nums(current_mem, tt)

    xpts = nx ; ypts = ny
    !xpts = fieldsize ; ypts = nz

    call generate_flat_file_name( var_name   = fldname,  &
                                  level_type = 'sig',         &
                                  level1     = int(ztop),     &
                                  level2     = int(zbot),     &
                                  gridnum    = gridnum,       &
                                  aoflag     = 'a',           &
                                  xpts       = xpts,          &
                                  ypts       = ypts,          &
                                  dtg        = dtg_pert,      &
                                  tau_hh     = ipert,         &
                                  tau_mm     = 0,             &
                                  tau_ss     = 0,             &
                                  field_type = 'pertfld',     &
                                  file_name  = coamps_init_pert )

    open( unit=file_unit,                                         &
          file=trim(dsnrff_pert)//coamps_init_pert,               &
          status=filestatus, access='direct', action=fileaction,  &
          form='unformatted', recl=filesize, iostat=io_status)
    call check_io_status(io_status, routine, source, revision,    &
                         revdate, 'Opening ' //                   &
                         coamps_init_pert // ' for read')

    call read_flat_file(file_unit, xpert)
    close(file_unit)

    xpert(:) = xpert(:) - xinit(:)

!------------------------------------------------------------------------------
! Read deterministic fields.
!------------------------------------------------------------------------------
    xpts = nx ; ypts = ny
    !xpts = fieldsize ; ypts = nz

    call generate_flat_file_name( var_name   = fldname,  &
                                  level_type = 'sig',         &
                                  level1     = int(ztop),     &
                                  level2     = int(zbot),     &
                                  gridnum    = gridnum,       &
                                  aoflag     = 'a',           &
                                  xpts       = xpts,          &
                                  ypts       = ypts,          &
                                  dtg        = cdtg,          &
                                  tau_hh     = tau_hr,        &
                                  tau_mm     = 0,             &
                                  tau_ss     = 0,             &
                                  field_type = 'analfld',     &
                                  file_name  = coamps_init_file )

    open( unit=file_unit,                                         &
          file=trim(dsnrff)//'analyses/'//coamps_init_file,       &
          status=filestatus, access='direct', action=fileaction,  &
          form='unformatted', recl=filesize, iostat=io_status)
    call check_io_status(io_status, routine, source, revision,    &
                         revdate, 'Opening ' //                   &
                         coamps_init_file // ' for read')

    call read_flat_file(file_unit, xinit)
    close(file_unit)

!------------------------------------------------------------------------------
! Perturb boundaries condition based on Torn et al MWR 2006.
! For tau_hr=0 perturb the mean field, for all other perturb the tendencies.
!------------------------------------------------------------------------------
    !if(tau_hr == 0) then
       xinit(:) = xinit(:) + xpert(:)
    !else
    !   xpert(:) = alpha * xpert_m1(:) + sqrt(1.0_r8 - alpha**2) * xpert(:)
    !   xinit(:) = xinit(:) + (xpert(:) - xpert_m1(:))/real(kgetbc*SECPERHR, kind=r8) 
    !end if
    !xpert_m1(:) = xpert(:)

!------------------------------------------------------------------------------
! Write perturbed boundary to file
!------------------------------------------------------------------------------
    fileaction = 'write'
    filestatus = 'unknown'
    open( unit=file_unit,                                         &
          file=trim(dsnrff)//'analyses/'//coamps_init_file,       &
          status=filestatus, access='direct', action=fileaction,  &
          form='unformatted', recl=filesize, iostat=io_status)
    call check_io_status(io_status, routine, source, revision,    &
                         revdate, 'Opening ' //                   &
                         coamps_init_file // ' for write')
    call write_flat_file(file_unit, xinit)
    close(file_unit)

  end do TIME_LOOP

  deallocate(xinit)
  deallocate(xpert)
  deallocate(xpert_m1)

  end do VAR_LOOP

  deallocate(ipert_nums)

end program perturb_init

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
