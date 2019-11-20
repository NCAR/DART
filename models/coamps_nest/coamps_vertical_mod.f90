! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

!------------------------------
! MODULE:       coamps_vertical_mod
! AUTHOR:       T. R. Whitcomb
!               Naval Research Laboratory
! DART VERSION: Jamaica
!               Manhattan (updated jun 2017)
!
! Module containing data structures and routines for dealing
! with the vertical component of a coamps domain
!------------------------------ 

module coamps_vertical_mod

    use coamps_util_mod, only : check_alloc_status

    use types_mod, only : r8

    use utilities_mod, only : error_handler, E_ERR

    implicit none

    private

    !------------------------------
    ! BEGIN PUBLIC INTERFACE
    !------------------------------
    public :: coamps_vertical
    public :: initialize_vertical

    public :: get_num_levels

    public :: get_msigma
    public :: get_wsigma
    public :: get_dsigmaw

    public :: dump_vertical_info

    !------------------------------
    ! END PUBLIC INTERFACE
    !------------------------------

    !------------------------------
    ! BEGIN EXTERNAL INTERFACES
    !------------------------------

    interface get_msigma
        module procedure get_msigma_array, get_msigma_value
    end interface get_msigma

    interface get_wsigma
        module procedure get_wsigma_array, get_wsigma_value
    end interface get_wsigma

    !------------------------------
    ! END EXTERNAL INTERFACES
    !------------------------------

    !------------------------------
    ! BEGIN TYPES AND CONSTANTS 
    !------------------------------

    type :: coamps_vertical
        private

        integer                              :: sigma_levels
        real(kind=r8), dimension(:), pointer :: dsigma
        real(kind=r8), dimension(:), pointer :: dsigmaw
        real(kind=r8), dimension(:), pointer :: msigma
        real(kind=r8), dimension(:), pointer :: wsigma
    end type coamps_vertical

    integer, parameter :: DATAHD_NUM_S_LEVELS   = 2
    integer, parameter :: DATAHD_DSIGMA_OFFSET  = 500
    integer, parameter :: DATAHD_MSIGMA_OFFSET  = 800

    ! Model vertical coordinate indexes top-down
    integer, parameter :: TOP_LEVEL    = 1

    !------------------------------
    ! END TYPES AND CONSTANTS 
    !------------------------------

    !------------------------------
    ! BEGIN MODULE VARIABLES
    !------------------------------

    ! version controlled file description for error handling, do not edit
    character(len=*), parameter :: source   = &
       "$URL$"
    character(len=*), parameter :: revision = "$Revision$"
    character(len=*), parameter :: revdate  = "$Date$"

    !------------------------------
    ! END MODULE VARIABLES
    !------------------------------

contains

    !------------------------------
    ! BEGIN PUBLIC ROUTINES
    !------------------------------

    ! initialize_vertical
    ! -------------------
    ! Initializes the COAMPS vertical coordinate data
    !  PARAMETERS
    !   IN  datahd            contents of COAMPS datahd file
    ! INOUT vertical          COAMPS vertical coordinate to populate
    subroutine initialize_vertical(datahd, vertical)
        real(kind=r8), dimension(:), intent(in)    :: datahd
        type(coamps_vertical),       intent(inout) :: vertical

        character(len=*), parameter :: routine = 'initialize_vertical'
        integer :: alloc_status

        ! Number of mass and w sigma levels
        integer :: m_levels, bottom_m_level
        integer :: w_levels, bottom_w_level
        integer :: cur_level 

        nullify(vertical%dsigma)
        nullify(vertical%dsigmaw)
        nullify(vertical%msigma)
        nullify(vertical%wsigma)

        vertical%sigma_levels = int(datahd(DATAHD_NUM_S_LEVELS))

        m_levels = vertical%sigma_levels
        w_levels = vertical%sigma_levels + 1

        bottom_m_level = m_levels
        bottom_w_level = w_levels

        allocate( vertical%dsigma(m_levels),  stat = alloc_status )
        call check_alloc_status(alloc_status, routine, source, revision, &
                                revdate, 'vertical%dsigma')
        allocate( vertical%dsigmaw(w_levels), stat=alloc_status)
        call check_alloc_status(alloc_status, routine, source, revision, &
                                revdate, 'vertical%dsigmaw')
        allocate( vertical%msigma(m_levels),  stat=alloc_status)
        call check_alloc_status(alloc_status, routine, source, revision, &
                                revdate, 'vertical%msigma')
        allocate( vertical%wsigma(w_levels),  stat=alloc_status)
        call check_alloc_status(alloc_status, routine, source, revision, &
                                revdate, 'vertical%wsigma')

        ! Data about the mass sigma levels can be directly read in...
        do cur_level = TOP_LEVEL, bottom_m_level
            vertical%dsigma(cur_level) = datahd(DATAHD_DSIGMA_OFFSET + &
                                                cur_level)
            vertical%msigma(cur_level) = datahd(DATAHD_MSIGMA_OFFSET + &
                                                cur_level)
        end do 

        ! ... but the w sigma levels are calculated.  This computation is
        ! based on the COAMPS utility package source code.
        vertical%dsigmaw(TOP_LEVEL)       = vertical%dsigma(TOP_LEVEL)*0.5
        vertical%dsigmaw(bottom_w_level)  = vertical%dsigma(bottom_m_level)*0.5
        do cur_level = 1, (m_levels-1)
            vertical%dsigmaw(cur_level + 1) = vertical%msigma(cur_level    )-&
                                              vertical%msigma(cur_level + 1)
        end do

        ! The w sigma levels are built from the lowest level up - remember that
        ! the "top" is actually level 1
        vertical%wsigma(bottom_w_level) = 0.0
        do cur_level = bottom_m_level, TOP_LEVEL, -1
            vertical%wsigma(cur_level) = vertical%wsigma(cur_level + 1) + &
                                         vertical%dsigma(cur_level    )
        end do
    end subroutine initialize_vertical

    ! get_num_levels
    ! --------------
    ! Returns the level of sigma levels
    !  PARAMETERS
    !   IN  vertical        COAMPS vertical structure to pull from
    function get_num_levels(vertical)
        type(coamps_vertical), intent(in)  :: vertical
        integer                            :: get_num_levels

        get_num_levels = vertical%sigma_levels
    end function get_num_levels

    ! get_msigma_array
    ! ----------------
    ! Returns the mass sigma levels for a given COAMPS domain
    !  PARAMETERS
    !   IN  vertical          COAMPS vertical structure  to pull from
    function get_msigma_array(vertical)
        type(coamps_vertical),       intent(in)  :: vertical
        real(kind=r8), dimension(:), pointer     :: get_msigma_array

        ! Don't attempt a copy if the data in the vertical won't fit
        !if ( size(vertical%msigma) > size(get_msigma_array) ) then
        !    call error_handler(E_ERR, 'get_msigma_array', 'msigma ' // &
        !                       'passed into function is too small',    &
        !                       source, revision, revdate)
        !end if

        get_msigma_array => vertical%msigma
    end function get_msigma_array

    ! get_msigma_value
    ! ----------------
    ! Returns a single mass sigma level at the given index
    !  PARAMETERS
    !   IN  vertical        COAMPS vertical structure to pull from
    !   IN  sigma_index     Index of the sigma level to extract
    function get_msigma_value(vertical, sigma_index)
        type(coamps_vertical), intent(in)  :: vertical
        integer,               intent(in)  :: sigma_index
        real(kind=r8)                      :: get_msigma_value

        get_msigma_value = vertical%msigma(sigma_index)
    end function get_msigma_value

    ! get_wsigma_array
    ! ----------
    ! Returns the w sigma levels for a given COAMPS domain
    !  PARAMETERS
    !   IN  vertical          COAMPS vertical structure to pull from
    function get_wsigma_array(vertical)
        type(coamps_vertical),       intent(in)  :: vertical
        real(kind=r8), dimension(:), pointer     :: get_wsigma_array

        ! Bail if the data in the vertical won't fit
        !if ( size(vertical%wsigma) > size(get_wsigma_array) ) then
        !    call error_handler(E_ERR, 'get_wsigma_array', 'wsigma ' // &
        !                       'passed into function is too small',    &
        !                       source, revision, revdate)
        !end if

        get_wsigma_array => vertical%wsigma
    end function get_wsigma_array

    ! get_wsigma_value
    ! ----------------
    ! Returns a single w sigma level at the given index
    !  PARAMETERS
    !   IN  vertical        COAMPS vertical structure to pull from
    !   IN  sigma_index     Index of the sigma level to extract
    function get_wsigma_value(vertical, sigma_index)
        type(coamps_vertical), intent(in)  :: vertical
        integer,               intent(in)  :: sigma_index
        real(kind=r8)                      :: get_wsigma_value

        get_wsigma_value = vertical%wsigma(sigma_index)
    end function get_wsigma_value

    ! get_dsigmaw
    ! -----------
    ! Returns the distance between w sigma levels for a given COAMPS
    ! vertical
    !  PARAMETERS
    !   IN  vertical          COAMPS vertical structure to pull from
    function get_dsigmaw(vertical)
        type(coamps_vertical),       intent(in)  :: vertical
        real(kind=r8), dimension(:), pointer     :: get_dsigmaw

        ! Bail if the data in the vertical won't fit
        !if ( size(vertical%dsigmaw) > size(get_dsigmaw) ) then
        !    call error_handler(E_ERR, 'get_dsigmaw', 'dsigmaw ' // &
        !    'passed into function is too small',     &
        !    source, revision, revdate)
        !end if

        get_dsigmaw => vertical%dsigmaw
    end function get_dsigmaw

    ! dump_vertical_info
    ! ------------------
    ! Dumps the COAMPS vertical structure in human-readable format
    !  PARAMETERS
    !   IN  vertical        COAMPS vertical structure to dump
    subroutine dump_vertical_info(vertical)
        type(coamps_vertical), intent(in) :: vertical

        integer :: kk

        write (*,*) "**** COAMPS VERTICAL ****"
        write (*,*) "Number of sigma levels:   ", vertical%sigma_levels

        write (*,*) "D(sigma levels)           "
        do kk = 1, vertical%sigma_levels
            write (*,*) vertical%dsigma(kk) 
        end do

        write (*,*) "Mass sigma levels:        "
        do kk = 1, vertical%sigma_levels
            write (*,*) vertical%msigma(kk) 
        end do

        write (*,*) "W sigma levels"
        do kk = 1, vertical%sigma_levels + 1
            write (*,*) vertical%wsigma(kk)
        end do

        write (*,*) "D(w sigma levels)"
        do kk = 1, vertical%sigma_levels + 1
            write (*,*) vertical%dsigmaw(kk)
        end do
    end subroutine dump_vertical_info



!TJH !***********************************************************************
!TJH ! interpolate from sigma levels to p-levels
!TJH 
!TJH subroutine s2pint(din,dout,zin,zout,kin,kout,nn)
!TJH 
!TJH ! #include <s2pint.prol>
!TJH ! CONFIGURATION IDENTIFICATION $HeadURL: https://coamps3.nrlmry.navy.mil/svn/atmos/branches/beta/libsrc/amlib/s2pint.F $ 
!TJH ! CONFIGURATION IDENTIFICATION @(#)$Id: s2pint.F 689 2010-12-07 21:29:08Z chen $
!TJH 
!TJH !  use domdec
!TJH implicit none
!TJH 
!TJH !***********************************************************************
!TJH ! arguments:
!TJH !***********************************************************************
!TJH 
!TJH integer, intent(in)  :: kin
!TJH integer, intent(in)  :: kout
!TJH integer, intent(in)  :: nn
!TJH real,    intent(in)  ::  din(iminf(nn):imaxf(nn),jminf(nn):jmaxf(nn),kin)
!TJH real,    intent(out) :: dout(iminf(nn):imaxf(nn),jminf(nn):jmaxf(nn),kout)
!TJH real,    intent(in)  ::  zin(iminf(nn):imaxf(nn),jminf(nn):jmaxf(nn),kin)
!TJH real,    intent(in)  :: zout(kout)
!TJH 
!TJH !***********************************************************************
!TJH ! local variables and dynamic storage:
!TJH !***********************************************************************
!TJH 
!TJH integer :: i
!TJH integer :: j
!TJH integer :: ki
!TJH integer :: ko
!TJH real    :: temp
!TJH 
!TJH !***********************************************************************
!TJH ! end of definitions
!TJH !***********************************************************************
!TJH 
!TJH do ko=1,kout
!TJH   do j = jminf(nn), jmaxf(nn)
!TJH     do i= iminf(nn), imaxf(nn)
!TJH       if (zout(ko) .le. zin(i,j,  1)) dout(i,j,ko)=din(i,j,  1)
!TJH       if (zout(ko) .ge. zin(i,j,kin)) dout(i,j,ko)=din(i,j,kin)
!TJH     enddo
!TJH   enddo
!TJH enddo
!TJH do ki=2,kin
!TJH   do ko=1,kout
!TJH     do j = jminf(nn), jmaxf(nn)
!TJH       do i= iminf(nn), imaxf(nn)
!TJH         if ((zout(ko) .gt. zin(i,j,ki-1)) .and.  &
!TJH             (zout(ko) .le. zin(i,j,ki  ))) then
!TJH           temp=(alog(zout(ko))-alog(zin(i,j,ki-1))) / (alog(zin(i,j,ki))-alog(zin(i,j,ki-1)))
!TJH           dout(i,j,ko)=din(i,j,ki-1)+temp*(din(i,j,ki) - din(i,j,ki-1))
!TJH         endif
!TJH       enddo
!TJH     enddo
!TJH   enddo
!TJH enddo
!TJH 
!TJH return
!TJH end
!TJH 
!TJH 
!TJH 
!TJH !***********************************************************************
!TJH ! interpolate from sigma-z levels to constant z levels
!TJH 
!TJH subroutine z2zint(din,dout,zin,zout,kin,kout,nn)
!TJH  
!TJH !  use domdec
!TJH 
!TJH ! #include <z2zint.prol>
!TJH ! CONFIGURATION IDENTIFICATION $HeadURL:
!TJH ! https://coamps3.nrlmry.navy.mil/svn/atmos/branches/beta/libsrc/amlib/z2zint.F $ 
!TJH ! CONFIGURATION IDENTIFICATION @(#)$Id: z2zint.F 689 2010-12-07 21:29:08Z chen $
!TJH 
!TJH implicit none
!TJH 
!TJH !***********************************************************************
!TJH !        parameters:
!TJH !***********************************************************************
!TJH 
!TJH integer, intent(in) :: kin
!TJH integer, intent(in) :: kout
!TJH integer, intent(in) :: nn
!TJH real,    intent(in) ::  din(iminf(nn):imaxf(nn),jminf(nn):jmaxf(nn),kin)
!TJH real,    intent(in) :: dout(iminf(nn):imaxf(nn),jminf(nn):jmaxf(nn),kout)
!TJH real,    intent(in) ::  zin(iminf(nn):imaxf(nn),jminf(nn):jmaxf(nn),kin)
!TJH real,    intent(in) :: zout(kout)
!TJH 
!TJH !***********************************************************************
!TJH !       local variables and dynamic storage:
!TJH !***********************************************************************
!TJH 
!TJH integer :: i
!TJH integer :: j
!TJH integer :: ki
!TJH integer :: ko
!TJH 
!TJH real :: temp
!TJH 
!TJH do ko=1,kout
!TJH   do j = jminf(nn), jmaxf(nn)
!TJH     do i= iminf(nn), imaxf(nn)
!TJH       if (zout(ko) .ge. zin(i,j,  1)) dout(i,j,ko)=din(i,j,  1)
!TJH       if (zout(ko) .le. zin(i,j,kin)) dout(i,j,ko)=din(i,j,kin)
!TJH     enddo
!TJH   enddo
!TJH enddo
!TJH 
!TJH do ki=2,kin
!TJH   do ko=1,kout
!TJH     do j = jminf(nn), jmaxf(nn)
!TJH       do i= iminf(nn), imaxf(nn)
!TJH         if ((zout(ko) .lt. zin(i,j,ki-1)) .and. &
!TJH             (zout(ko) .ge. zin(i,j,ki))) then
!TJH         temp = (zout(ko)-zin(i,j,ki-1))/(zin(i,j,ki)   - zin(i,j,ki-1))
!TJH         dout(i,j,ko) = din(i,j,ki-1)+temp*(din(i,j,ki) - din(i,j,ki-1))
!TJH         endif
!TJH       enddo
!TJH     enddo
!TJH   enddo
!TJH enddo
!TJH 
!TJH return
!TJH end


!------------------------------
! END PUBLIC ROUTINES
!------------------------------

!------------------------------
! BEGIN PRIVATE ROUTINES
!------------------------------
!------------------------------
! END PRIVATE ROUTINES
!------------------------------

end module coamps_vertical_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
