! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

!------------------------------
! MODULE:       coamps_flatfile_mod
! AUTHOR:       T. R. Whitcomb and P. A. Reinecke
!               Naval Research Laboratory
! DART VERSION: Jamaica
!
! Module containing the data structure and routines for dealing with
! COAMPS flat files.
!------------------------------ 

module coamps_flat_file_mod

  use coamps_util_mod,      only : C_REAL,               & 
                                   C_REAL4,              & 
                                   check_io_status,      &
                                   check_alloc_status,   &
                                   check_dealloc_status, &
                                   fix_for_platform4,    &
                                   uppercase,            &
                                   lowercase

  implicit none
  private

  !------------------------------
  ! BEGIN PUBLIC INTERFACE
  !------------------------------

  ! Initialization
  public :: read_flat_file
  public :: write_flat_file
  public :: generate_one_flat_file_name
  !------------------------------
  ! END PUBLIC INTERFACE
  !------------------------------

  !------------------------------
  ! BEGIN EXTERNAL INTERFACES
  !------------------------------
  ! [none]
  !------------------------------
  ! END EXTERNAL INTERFACES
  !------------------------------


  !------------------------------
  ! BEGIN TYPES AND CONSTANTS 
  !------------------------------
  ! [none]
  !------------------------------
  ! END TYPES AND CONSTANTS 
  !------------------------------

  !------------------------------
  ! BEGIN MODULE VARIABLES
  !------------------------------

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source = &
   "$URL$"
character(len=32 ), parameter :: revision = "Revision: 4371"
character(len=128), parameter :: revdate = "Date: 2010-05-21 16:23:38 -0600 (Fri, 21 May 2010) "

  logical :: module_initialized = .false.
  !------------------------------
  ! END MODULE VARIABLES
  !------------------------------
contains

  !------------------------------
  ! BEGIN PUBLIC ROUTINES
  !------------------------------

  ! write_flat_file
  ! ----------------
  ! Given the unit number of an *open* COAMPS flat 
  ! file, read the file into an array. 
  !  PARAMETERS
  !   IN  flat_unit      unit number of an open flat file
  !   OUT flat_array     coamps_grid structure to be filled
  subroutine write_flat_file(flat_unit, flat_array)
    integer,                         intent(in) :: flat_unit
    real(kind=C_REAL), dimension(:), intent(in) :: flat_array

    real(kind=C_REAL4), dimension(:), allocatable :: flat_array_tmp
    character(len=*), parameter :: routine = 'write_flat_file'
    integer :: io_status, alloc_status, dealloc_status
    integer :: field_size

    field_size=size(flat_array)
    allocate(flat_array_tmp(field_size), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'flat_array_tmp')

    ! COAMPS flat files are real(kind=4)

    flat_array_tmp(:) = real(flat_array(:), kind=C_REAL4)
    call fix_for_platform4(flat_array_tmp, field_size, C_REAL4)

    write(unit=flat_unit, rec=1, iostat=io_status) flat_array_tmp
    call check_io_status(io_status, routine, source, revision, &
                         revdate, 'writing flat file')

    deallocate(flat_array_tmp, stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'flat_array_tmp')
  end subroutine write_flat_file

  ! read_flat_file
  ! ----------------
  ! Given the unit number of an *open* COAMPS flat 
  ! file, read the file into an array. 
  !  PARAMETERS
  !   IN  flat_unit      unit number of an open flat file
  !   OUT flat_array     coamps_grid structure to be filled
  subroutine read_flat_file(flat_unit, flat_array)
    integer, intent(in)                          :: flat_unit
    real(kind=C_REAL), dimension(:), intent(out) :: flat_array

    real(kind=C_REAL4), dimension(:), allocatable :: flat_array_tmp
    character(len=*), parameter :: routine = 'read_flat_file'
    integer :: io_status, alloc_status, dealloc_status
    integer :: field_size

    field_size=size(flat_array)
    allocate(flat_array_tmp(field_size), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'flat_array_tmp')

    ! Read in the data - COAMPS writes flat files as C_REAL4's
    read(unit=flat_unit, rec=1, iostat=io_status) flat_array_tmp
    call check_io_status(io_status, routine, source, revision, &
                         revdate, 'Reading flat file')
    call fix_for_platform4(flat_array_tmp, field_size, C_REAL4)
    flat_array(:)=real(flat_array_tmp(:) , kind=C_REAL)

    deallocate(flat_array_tmp, stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'flat_array_tmp')
  end subroutine read_flat_file

  ! generate_one_flat_file_name
  ! -----------------------
  ! Given field, level, and grid information, generate the properly
  ! formatted 64-character COAMPS flat file name.  Note that this
  ! does *not* generate any path information - it only returns the
  ! file name.
  !  PARAMETERS
  !   IN  var_name          the field the file contains
  !   IN  level_type        vertical level type (height/pressure/etc)
  !   IN  level1            lowest vertical level in the file
  !   IN  level2            highest vertical level in the file
  !                         (for files for a single level, level1 is 
  !                          that level and level2 is left to 0)
  !   IN  gridnum           nest number (only 1 supported for now)
  !   IN  aoflag            field type: (a)tmosphere or (o)cean
  !   IN  xpts              number of points in the x direction
  !   IN  ypts              number of points in the y direction
  !   IN  dtg               base date-time group
  !   IN  tau_hh            forecast lead time - hour component
  !   IN  tau_mm            forecast lead time - minute component
  !   IN  tau_ss            forecast lead time - second component
  !   IN  field_type        type of field (e.g. fcstfld, infofld)
  !   OUT file_name         COAMPS flat file name
  subroutine generate_one_flat_file_name(var_name, level_type, level1, &
                                         level2, gridnum, aoflag, xpts,&
                                         ypts, dtg, tau_hh, tau_mm,    &
                                         tau_ss, field_type,           &
                                         file_name)
    character(len=6),  intent(in)  :: var_name
    character(len=3),  intent(in)  :: level_type
    integer,           intent(in)  :: level1
    integer,           intent(in)  :: level2
    integer,           intent(in)  :: gridnum
    character(len=1),  intent(in)  :: aoflag
    integer,           intent(in)  :: xpts
    integer,           intent(in)  :: ypts
    character(len=10), intent(in)  :: dtg
    integer,           intent(in)  :: tau_hh
    integer,           intent(in)  :: tau_mm
    integer,           intent(in)  :: tau_ss
    character(len=7),  intent(in)  :: field_type
    character(len=64), intent(out) :: file_name

    write(file_name, 100) var_name, level_type,               &
         &  level1, level2, gridnum, aoflag, xpts, ypts, dtg, &
         &  tau_hh, tau_mm, tau_ss, field_type

    ! make sure the file name is lower case
    file_name=lowercase(file_name)

100 format(A6,'_',A3,'_',I6.6,'_',I6.6,'_',I1,A1,I4.4,'x',I4.4,'_', &
           A10,'_',I4.4,I2.2,I2.2,'_',A7)
  end subroutine generate_one_flat_file_name

  !------------------------------
  ! END PUBLIC ROUTINES
  !------------------------------

  !------------------------------
  ! BEGIN PRIVATE ROUTINES
  !------------------------------

  !------------------------------
  ! END PRIVATE ROUTINES
  !------------------------------
end module

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
