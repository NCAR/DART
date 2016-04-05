! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

!------------------------------
! MODULE:       coamps_nest_mod
! AUTHOR:       T. R. Whitcomb
!               Naval Research Laboratory
! DART VERSION: Jamaica
!
! Module containing data structures and routines for dealing
! with the nest component of a coamps domain
!------------------------------ 

module coamps_nest_mod

    use coamps_util_mod, only : generate_flat_file_name,    &
                                read_flat_file,             &
                                check_alloc_status,         &
                                check_io_status,            &
                                check_dealloc_status

    use coamps_map_mod,  only : coamps_grid, grid_point_to_latlon

    use types_mod, only : r4, r8

    use utilities_mod, only: E_ERR, error_handler, get_unit

    implicit none

    private

    !------------------------------
    ! BEGIN PUBLIC INTERFACE
    !------------------------------
    public :: nest_point
    public :: coamps_nest

    public :: initialize_nest
    public :: initialize_nest_latlon

    public :: register_child_nest
    public :: register_parent_nest
    public :: get_parent_nest_id
    public :: set_nest_id
    public :: get_nest_id
    public :: get_nest
    public :: get_nest_latlon

    public :: normalize_x_distance
    public :: normalize_y_distance

    public :: nest_point_to_coarse_point
    public :: coarse_point_to_nest_point

    public :: get_nest_size
    public :: get_nest_level
    public :: get_nest_i_width
    public :: get_nest_j_width
    public :: get_nest_delta_x
    public :: get_nest_delta_y
    public :: nest_index_2d_to_1d
    public :: nest_index_1d_to_3d

    public :: get_i_coord
    public :: get_j_coord

    public :: in_this_nest

    public :: make_nest_point

    public :: get_terrain
    public :: get_terrain_height_at_points

    public :: decompose_nest
    public :: get_num_subnests
    public :: get_subnest_i_width
    public :: get_subnest_j_width
    public :: get_subnest_iminf
    public :: get_subnest_imaxf
    public :: get_subnest_jminf
    public :: get_subnest_jmaxf

    public :: child_to_parent
    public :: parent_to_child

    public :: dump_nest_info

    public :: assignment(=)

    !------------------------------
    ! END PUBLIC INTERFACE
    !------------------------------

    !------------------------------
    ! BEGIN EXTERNAL INTERFACES
    !------------------------------

    interface make_nest_point
        module procedure make_nest_point_r8, make_nest_point_int
    end interface make_nest_point

    interface assignment(=)
        module procedure assign_nest_point
    end interface

    interface initialize_list
        module procedure initialize_coamps_nest_list
    end interface initialize_list

    interface add_to_list
        module procedure  add_coamps_nest_to_list
    end interface add_to_list

    interface print_list
        module procedure  print_coamps_nest_list
    end interface print_list

    interface is_empty_list
        module procedure  is_coamps_nest_list_empty
    end interface is_empty_list

    interface get_iterator
        module procedure  get_coamps_nest_list_iterator
    end interface get_iterator

    interface has_next
        module procedure  coamps_nest_iterator_has_next
    end interface has_next

    interface get_next
        module procedure  coamps_nest_iterator_get_next
    end interface get_next

    !------------------------------
    ! END EXTERNAL INTERFACES
    !------------------------------

    !------------------------------
    ! BEGIN TYPES AND CONSTANTS 
    !------------------------------
    type :: coamps_nest_list
        private
        type(coamps_nest_list_node), pointer :: head
        type(coamps_nest_list_node), pointer :: tail
    end type coamps_nest_list

    type :: coamps_nest_list_node
        private
        type(coamps_nest),           pointer :: value
        type(coamps_nest_list_node), pointer :: next
    end type coamps_nest_list_node

    type :: coamps_nest_iterator
        private
        type(coamps_nest_list_node), pointer :: current
    end type coamps_nest_iterator

    type :: coamps_nest
        private

        integer :: id

        integer            :: pts_x
        integer            :: pts_y
        real(kind=r8)      :: delta_x
        real(kind=r8)      :: delta_y
        integer            :: anchor_i
        integer            :: anchor_j
        integer            :: nest_level

        integer                    :: parent_nest_id
        type(coamps_nest), pointer :: parent_nest
        type(coamps_nest_list)     :: child_nests

        ! Domain decomposition
        integer, dimension(:), pointer :: iminf
        integer, dimension(:), pointer :: imaxf
        integer, dimension(:), pointer :: jminf
        integer, dimension(:), pointer :: jmaxf

        real(kind=r8), dimension(:,:), pointer :: terrain
        real(kind=r8), dimension(:,:), pointer :: lat
        real(kind=r8), dimension(:,:), pointer :: lon
    end type coamps_nest

    type :: nest_point
        !TESTONLY private
        type(coamps_nest), pointer :: nest
        real(kind=r8)              :: ii
        real(kind=r8)              :: jj
    end type nest_point

    ! Nest data - the main offset is given by
    ! (nest_num * nest_offset) and the data offsets are as follows
    integer, parameter :: DATAHD_NEST_OFFSET    = 30
    integer, parameter :: DATAHD_NUM_X_POINTS   = 0
    integer, parameter :: DATAHD_NUM_Y_POINTS   = 1
    integer, parameter :: DATAHD_ANCHOR_POINT_I = 2
    integer, parameter :: DATAHD_ANCHOR_POINT_J = 3
    integer, parameter :: DATAHD_PARENT_NEST    = 6
    integer, parameter :: DATAHD_DELTA_X        = 7
    integer, parameter :: DATAHD_DELTA_Y        = 8
    integer, parameter :: DATAHD_NEST_LEVEL     = 17

    integer, parameter :: COARSE_NEST = 1
    !------------------------------
    ! END TYPES AND CONSTANTS 
    !------------------------------

    !------------------------------
    ! BEGIN MODULE VARIABLES
    !------------------------------

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

    !------------------------------
    ! END MODULE VARIABLES
    !------------------------------

contains

    !------------------------------
    ! BEGIN PUBLIC ROUTINES
    !------------------------------

    ! initialize_nest
    ! ---------------
    ! Populates a COAMPS nest based on the given nest number and datahd
    ! record
    ! INOUT nest              COAMPS nest to initialize
    !   IN  dtg               COAMPS date-time-group
    !   IN  datahd            datahd record to source
    subroutine initialize_nest(nest, dtg, datahd)
        type(coamps_nest),           intent(inout) :: nest
        character(len=10),           intent(in)    :: dtg
        real(kind=r8), dimension(:), intent(in)    :: datahd

        integer :: nest_offset

        nest_offset = nest%id * DATAHD_NEST_OFFSET

        nest%pts_x = int(datahd(DATAHD_NUM_X_POINTS + nest_offset))
        nest%pts_y = int(datahd(DATAHD_NUM_Y_POINTS + nest_offset))

        nest%delta_x = datahd(DATAHD_DELTA_X + nest_offset)
        nest%delta_y = datahd(DATAHD_DELTA_Y + nest_offset)

        nest%anchor_i = int(datahd(DATAHD_ANCHOR_POINT_I + nest_offset))
        nest%anchor_j = int(datahd(DATAHD_ANCHOR_POINT_J + nest_offset))

        nest%parent_nest_id = int(datahd(DATAHD_PARENT_NEST + nest_offset))

        nest%nest_level = int(datahd(DATAHD_NEST_LEVEL + nest_offset))

        call initialize_list(nest%child_nests)
        call read_terrain_height(nest, dtg)

        nullify(nest%iminf)
        nullify(nest%imaxf)
        nullify(nest%jminf)
        nullify(nest%jmaxf)
    end subroutine initialize_nest

    ! register_parent_nest
    ! --------------------
    ! Given a child nest and a parent nest, registers the parent nest
    ! in the child nest as its parent
    !  PARAMETERS
    !   IN  parent_nest     The parent nest to register
    ! INOUT child_nest      The nest to register in
    subroutine register_parent_nest(parent_nest, child_nest)
        type(coamps_nest), target, intent(in)    :: parent_nest
        type(coamps_nest),         intent(inout) :: child_nest

        child_nest%parent_nest => parent_nest
    end subroutine register_parent_nest

    ! register_child_nest
    ! -------------------
    ! Given a child nest and a parent nest, registers the child nest in
    ! the parent nest's list of child nests
    ! as a child nest
    !  PARAMETERS
    !   IN  child_nest      The child nest to register
    ! INOUT parent_nest     The nest to register in
    subroutine register_child_nest(child_nest, parent_nest)
        type(coamps_nest), intent(in)    :: child_nest
        type(coamps_nest), intent(inout) :: parent_nest

        ! Skip writing the coarse nest as its own child since this messes with
        ! the coarse -> fine grid conversion
        if (child_nest%id .eq. COARSE_NEST) return

        call add_to_list(child_nest, parent_nest%child_nests) 
    end subroutine register_child_nest

    ! get_parent_nest_id
    ! ------------------
    ! Returns a nest's parent nest ID
    !  PARAMETERS
    !   IN  nest            The nest to query
    function get_parent_nest_id(nest)
        type(coamps_nest), intent(in) :: nest
        integer                       :: get_parent_nest_id

        get_parent_nest_id = nest%parent_nest_id
    end function get_parent_nest_id

    ! set_parent_nest
    ! ---------------
    ! Sets a given nest as the current nest's parent 
    !  PARAMETERS
    ! INOUT this_nest       The nest to set the parent of
    !   IN  parent_nest     The parent nest
    subroutine set_parent_nest(this_nest, parent_nest)
        type(coamps_nest),         intent(inout) :: this_nest
        type(coamps_nest), target, intent(in)    :: parent_nest

        this_nest%parent_nest => parent_nest 
    end subroutine set_parent_nest
    
    ! set_nest_id
    ! -----------
    ! Sets the nest's ID - the nest number in the domain
    !  PARAMETERS
    ! INOUT nest            The nest to set the ID of
    !   IN  id              The ID to set
    subroutine set_nest_id(nest, id)
        type(coamps_nest), intent(inout) :: nest
        integer,           intent(in)    :: id

        nest%id = id
    end subroutine set_nest_id

    ! get_nest_id
    ! -----------
    ! Returns the nests's ID - the nest number in the domain
    !  PARAMETERS
    !   IN nest             The nest to query the ID of
    function get_nest_id(nest)
        type(coamps_nest), intent(in)  :: nest
        integer                        :: get_nest_id 

        get_nest_id = nest%id
    end function get_nest_id

    ! normalize_x_distance
    ! --------------------
    ! Rescales an x distance using the grid spacing
    function normalize_x_distance(nest, x_distance)
        type(coamps_nest), intent(in)  :: nest
        real(kind=r8),     intent(in)  :: x_distance
        real(kind=r8)                  :: normalize_x_distance

        normalize_x_distance = x_distance
    end function normalize_x_distance

    ! normalize_y_distance
    ! --------------------
    ! Rescales an y distance using the grid spacing
    function normalize_y_distance(nest, y_distance)
        type(coamps_nest), intent(in)  :: nest
        real(kind=r8),     intent(in)  :: y_distance
        real(kind=r8)                  :: normalize_y_distance

        normalize_y_distance = y_distance
    end function normalize_y_distance

    ! get_terrain
    ! ----------------------------
    ! Given a COAMPS nest returns a pointer to the terrain on the nest
    !   IN  nest              coamps_nest to pull terrht from
    function get_terrain(nest) result(terrain)
        type(coamps_nest), intent(in)           :: nest
        real(kind=r8), dimension(:,:), pointer  :: terrain 

        terrain => nest%terrain
    end function get_terrain

    ! get_terrain_height_at_points
    ! ----------------------------
    ! Given a COAMPS nest and a vector representation of a set of 
    ! (i,j) points, returns the terrain height for those points
    !  PARAMETERS
    !   IN  grid              coamps_grid to pull terrht from
    !   IN  points_i          points's x-direction index
    !   IN  points_j          points's y-direction index
    !   OUT terrain_height    points's terrain height
    subroutine get_terrain_height_at_points(nest, points_i, points_j,&
                                            terrain_height)
        type(coamps_nest),           intent(in)  :: nest
        integer, dimension(:),       intent(in)  :: points_i
        integer, dimension(:),       intent(in)  :: points_j
        real(kind=r8), dimension(:), intent(out) :: terrain_height

        integer :: nn

        ! Sanity check
        if (size(points_i) .ne. size(points_j)) then
            call error_handler(E_ERR, 'get_terrain_height_at_points',   &
                               'Error in getting terrain height - ' //  &
                               'i/j point arrays are different sizes!', &
                               source, revision, revdate)
        end if

        do nn=1,size(points_i)
            terrain_height(nn) = nest%terrain(points_i(nn), points_j(nn))
        end do
    end subroutine get_terrain_height_at_points

    ! decompose_nest
    ! --------------
    ! Given the nest, the number of I/O processors, the number of
    ! processors in x and y, and the number of boundary points (nbnam
    ! in the COAMPS namelist), mimics the domain decomposition
    ! performed by COAMPS routine domdec.F (the math is taken from that
    ! module and variables renamed for use here).  The maximum and
    ! minimum [ij] points for each processor are given in arrays
    ! [ij](min|max)f within the coamps_nest structure.  This routine
    ! populates those arrays.  
    !  PARAMETERS
    ! INOUT nest              coamps_nest to decompose
    !   IN  num_io_procs      number of I/O processors (0 or 1)
    !   IN  dom_proc_x        number of processors in the x direction
    !   IN  dom_proc_y        number of processors in the y direction
    !   IN  nbound            number of halo boundary points
    subroutine decompose_nest(nest, num_io_procs, dom_proc_x, &
                              dom_proc_y, nbound) 
        type(coamps_nest), intent(inout) :: nest
        integer, intent(in)              :: num_io_procs
        integer, intent(in)              :: dom_proc_x
        integer, intent(in)              :: dom_proc_y
        integer, intent(in)              :: nbound

        ! Total number of processors we are iterating over, along with
        ! the current processor and some counters
        integer :: nproc_x, nproc_y
        integer :: totalprocs, curproc, x_count, y_count

        ! Number of points per processor (not including the boundary)
        integer :: xperproc, yperproc
        integer :: xextra, yextra

        ! The initial boundaries for the domain (not counting the
        ! halo points)
        integer, dimension(:), pointer :: begin_i
        integer, dimension(:), pointer :: begin_j
        integer, dimension(:), pointer :: end_i
        integer, dimension(:), pointer :: end_j

        ! Error handling
        character(len=*), parameter :: routine = 'decompose_nest'
        integer :: alloc_status, dealloc_status

        ! If we are only doing single processor I/O, it doesn't matter
        ! how many processors there are - we treat it as a single
        ! decomposition
        if (num_io_procs .gt. 0) then
            nproc_x = 1
            nproc_y = 1
        else
            nproc_x = dom_proc_x
            nproc_y = dom_proc_y
        endif
        totalprocs = nproc_x * nproc_y

        ! Temporary storage
        nullify(begin_i)
        nullify(begin_j)
        nullify(end_i)
        nullify(end_j)
        allocate( begin_i(totalprocs), stat=alloc_status )
        call check_alloc_status(alloc_status, routine, source, revision, &
                                revdate, 'begin_i')
        allocate( begin_j(totalprocs), stat=alloc_status )
        call check_alloc_status(alloc_status, routine, source, revision, &
                                revdate, 'begin_j')
        allocate( end_i(totalprocs), stat=alloc_status )
        call check_alloc_status(alloc_status, routine, source, revision, &
                                revdate, 'end_i')
        allocate( end_j(totalprocs), stat=alloc_status )
        call check_alloc_status(alloc_status, routine, source, revision, &
                                revdate, 'end_j')

        ! Permanent storage
        call initialize_decomposition(nest, totalprocs)

        ! Get the processor-point counts
        xperproc = (nest%pts_x - 2) / nproc_x
        yperproc = (nest%pts_y - 2) / nproc_y
        xextra   = (nest%pts_x - 2) - (xperproc * nproc_x)
        yextra   = (nest%pts_y - 2) - (yperproc * nproc_y)

        do y_count = 1,nproc_y
        do x_count = 1,nproc_x

        ! 2D grid -> 1D number
        curproc = (y_count - 1) * nproc_x + x_count

        ! X points
        if (x_count .eq. 1) then
            begin_i(curproc) = 2
        else
            begin_i(curproc) = end_i(curproc-1) + 1
        endif

        end_i(curproc) = begin_i(curproc) + xperproc - 1
        if (x_count .gt. nproc_x-xextra) then
            end_i(curproc) = end_i(curproc) + 1
        endif

        ! Y Points
        if (y_count .eq. 1) then
            begin_j(curproc) = 2
        else
            begin_j(curproc) = end_j(curproc-nproc_x) + 1
        endif
        end_j(curproc) = begin_j(curproc) + yperproc - 1
        if (y_count .gt. nproc_y-yextra) then
            end_j(curproc) = end_j(curproc) + 1
        endif

        ! Take the boundaries into consideration
        nest%iminf(curproc) = begin_i(curproc) - nbound
        nest%imaxf(curproc) = end_i(curproc)   + nbound
        nest%jminf(curproc) = begin_j(curproc) - nbound
        nest%jmaxf(curproc) = end_j(curproc)   + nbound
        enddo
        enddo

        deallocate( begin_i, begin_j, end_i, end_j, stat=dealloc_status )
        call check_dealloc_status(dealloc_status, routine, source, &
                                  revision, revdate, '[ij](min|max)f')
    end subroutine decompose_nest

    
    ! get_subnest_iminf
    ! -----------------
    ! Accessor method for the index of a subnest's left extent
    !  PARAMETERS
    !   IN  nest              coamps_nest to get data from
    !   IN  subnest_num       decomposed subnest number
    function get_subnest_iminf(nest, subnest_num)
        type(coamps_nest), intent(in)  :: nest
        integer,           intent(in)  :: subnest_num
        integer                        :: get_subnest_iminf

        get_subnest_iminf = nest%iminf(subnest_num)
    end function get_subnest_iminf

    ! get_subnest_imaxf
    ! -----------------
    ! Accessor method for the index of a subnest's right extent
    !  PARAMETERS
    !   IN  nest              coamps_nest to get data from
    !   IN  subnest_num       decomposed subnest number
    function get_subnest_imaxf(nest, subnest_num)
        type(coamps_nest), intent(in)  :: nest
        integer,           intent(in)  :: subnest_num
        integer                        :: get_subnest_imaxf

        get_subnest_imaxf = nest%imaxf(subnest_num)
    end function get_subnest_imaxf

    ! get_subnest_jminf
    ! -----------------
    ! Accessor method for the index of a subnest's right extent
    !  PARAMETERS
    !   IN  nest              coamps_nest to get data from
    !   IN  subnest_num       decomposed subnest number
    function get_subnest_jminf(nest, subnest_num)
        type(coamps_nest), intent(in)  :: nest
        integer,           intent(in)  :: subnest_num
        integer                        :: get_subnest_jminf

        get_subnest_jminf = nest%jminf(subnest_num)
    end function get_subnest_jminf

    ! get_subnest_jmaxf
    ! -----------------
    ! Accessor method for the index of a subnest's right extent
    !  PARAMETERS
    !   IN  nest              coamps_nest to get data from
    !   IN  subnest_num       decomposed subnest number
    function get_subnest_jmaxf(nest, subnest_num)
        type(coamps_nest), intent(in)  :: nest
        integer,           intent(in)  :: subnest_num
        integer                        :: get_subnest_jmaxf

        get_subnest_jmaxf = nest%jmaxf(subnest_num)
    end function get_subnest_jmaxf

    ! get_num_subnests
    ! ----------------
    ! Returns the number of subnests computed for a particular nest
    ! Note that this only makes sense after the nest has 
    !  PARAMETERS
    !   IN  nest            coamps_nest to query 
    function get_num_subnests(nest)
        type(coamps_nest), intent(in)  :: nest
        integer                        :: get_num_subnests

        if (associated(nest%iminf)) then
            get_num_subnests = size(nest%iminf)
        else
            call error_handler(E_ERR, 'get_num_subnests', 'Error in ' // &
                               'computing number of subnests - nest ' // &
                               'has not been decomposed yet!', source,   &
                               revision, revdate)
        end if
    end function get_num_subnests

    ! in_this_nest
    ! ------------
    ! Returns true if a given point is valid (the nest contains the point)
    !  PARAMETERS
    !   IN  test_pt           Point to test the validity of
    function in_this_nest(test_pt)
        type(nest_point),    intent(in)  :: test_pt
        logical                          :: in_this_nest

        !FIXME fixed for now, should vary with the coamps parameter.
        integer, parameter :: nbdypt = 5 

        real(kind=r8) :: i_min
        real(kind=r8) :: j_min
        real(kind=r8) :: i_max
        real(kind=r8) :: j_max

        if(.not.associated(test_pt%nest)) then
         in_this_nest = .false.
         return
        end if

        i_min = real(1 + nbdypt, kind=r8)
        j_min = real(1 + nbdypt, kind=r8)
        i_max = real(test_pt%nest%pts_x - nbdypt, kind=r8)
        j_max = real(test_pt%nest%pts_y - nbdypt, kind=r8)

        if ((test_pt%ii .ge. i_min) .and. &
            (test_pt%jj .ge. j_min) .and. &
            (test_pt%ii .le. i_max) .and. &
            (test_pt%jj .le. j_max)) then
            in_this_nest = .true.
        else
            in_this_nest = .false.
        end if
    end function in_this_nest

    ! make_nest_point_r8
    ! ---------------
    ! Construct a nest point structure
    !  PARAMETERS
    !   IN  nest            The COAMPS nest the point is in
    !   IN  i_coord         The i-coordinate
    !   IN  j_coord         The j-coordinate
    function make_nest_point_r8(nest, i_coord, j_coord) result(make_nest_point)
        type(coamps_nest), target, intent(in)  :: nest
        real(kind=r8),             intent(in)  :: i_coord
        real(kind=r8),             intent(in)  :: j_coord
        type(nest_point)                       :: make_nest_point

        make_nest_point%nest => nest
        make_nest_point%ii   =  i_coord
        make_nest_point%jj   =  j_coord
    end function make_nest_point_r8

    ! make_nest_point_int
    ! ---------------
    ! Construct a nest point structure
    !  PARAMETERS
    !   IN  nest            The COAMPS nest the point is in
    !   IN  i_coord         The i-coordinate
    !   IN  j_coord         The j-coordinate
    function make_nest_point_int(nest, i_coord, j_coord) result(make_nest_point)
        type(coamps_nest), target, intent(in)  :: nest
        integer,                   intent(in)  :: i_coord
        integer,                   intent(in)  :: j_coord
        type(nest_point)                       :: make_nest_point

        make_nest_point%nest => nest
        make_nest_point%ii   =  real(i_coord, kind=r8)
        make_nest_point%jj   =  real(j_coord, kind=r8)
    end function make_nest_point_int

    ! get_i_coord
    ! -----------
    ! Accessor for nest point structure
    function get_i_coord(nest_pt)
        type(nest_point), intent(in)  :: nest_pt
        real(kind=r8)                 :: get_i_coord

        get_i_coord = nest_pt%ii
    end function get_i_coord

    ! get_j_coord
    ! -----------
    ! Accessor for nest point structure
    function get_j_coord(nest_pt)
        type(nest_point), intent(in)  :: nest_pt
        real(kind=r8)                 :: get_j_coord

        get_j_coord = nest_pt%jj
    end function get_j_coord

    ! get_nest
    ! --------
    ! Accessor for nest point structure
    function get_nest(nest_pt)
        type(nest_point),           intent(in)  :: nest_pt
        type(coamps_nest), pointer              :: get_nest

        get_nest => nest_pt%nest
    end function get_nest

    ! child_to_parent
    ! ---------------
    ! Convert an i or j coordinate from a child nest to its equivalent 
    ! coordinate on its parent nest
    !  PARAMETERS
    !   IN  point             i or j point on child nest       
    !   IN  anchor            parent point corresponding to child's "1" 
    function child_to_parent(point, anchor)
        real(kind=r8), intent(in)  :: point
        integer,       intent(in)  :: anchor    
        real(kind=r8)              :: child_to_parent

        child_to_parent = (point - 1) / 3_r8 + anchor
    end function child_to_parent

    ! parent_to_child
    ! ---------------
    ! Convert an i or j coordinate from a parent nest to its equivalent
    ! coordinate within the child nest
    ! N.B. Since the parent nest generally contains areas not covered by
    !      child nests, it is possible that this function will return a 
    !      value outside the child nest boundaries
    !  PARAMETERS
    !   IN  point             i or j point on parent nest
    !   IN  anchor            parent point corresponding to child's "1"
    function parent_to_child(point, anchor)
        real(kind=r8), intent(in)  :: point
        integer,       intent(in)  :: anchor    
        real(kind=r8)              :: parent_to_child

        parent_to_child = (point - anchor) * 3_r8 + 1
    end function parent_to_child

    ! assign_nest_point
    ! -----------------
    ! Allow overloading of assignment operator for nest_point type
    ! and just copy everything over
    !  PARAMETERS
    !   OUT this_pt           point to assign to
    !   IN  new_pt            point to assign from
    subroutine assign_nest_point(this_pt, new_pt)
        type(nest_point), intent(out) :: this_pt
        type(nest_point), intent(in)  :: new_pt

        this_pt%ii   =  new_pt%ii
        this_pt%jj   =  new_pt%jj

        ! There is only one nest
        this_pt%nest => new_pt%nest
    end subroutine assign_nest_point

    ! dump_nest_info
    ! --------------
    ! Dumps the nest information in human-readable form
    !  PARAMETERS
    !   IN  nest            COAMPS nest to dump
    !   IN  dump_terrain    (Optional) Write the terrain to a file
    subroutine dump_nest_info(nest, dump_terrain)
        type(coamps_nest), intent(in) :: nest
        logical, optional, intent(in) :: dump_terrain

        character(len=8) :: terrain_name
        integer          :: terrain_unit
        integer          :: reclen
        integer          :: io_status

        write (*,*) "  ** NEST", nest%id, "**"
        write (*,*) "  X-direction spacing (m):  ", &
                    nest%delta_x
        write (*,*) "  Y-direction spacing (m):  ", &
                    nest%delta_y
        write (*,*) "  Number of X grid points:  ", &
                    nest%pts_x
        write (*,*) "  Number of Y grid points:  ", &
                    nest%pts_y
        write (*,*) "  i point on parent nest:   ", &
                    nest%anchor_i
        write (*,*) "  j point on parent nest:   ", &
                    nest%anchor_j
        write (*,*) "  Nest level:",                &
                    nest%nest_level
        write (*,*) "  Parent nest:              ", &
                    nest%parent_nest_id     
        write (*,*) "  Child nest(s):"
                    call print_list(nest%child_nests)

        if (present(dump_terrain)) then
            if (dump_terrain) then
                inquire(iolength=reclen) nest%terrain 
                write(terrain_name, '(A7,I1)') 'terrain', nest%id
                terrain_unit = get_unit()
                open(unit = terrain_unit, file = terrain_name, &
                     status ='replace', access='direct',       &
                     action = 'write', form = 'unformatted',   &
                     recl = reclen, iostat = io_status)
                if (io_status .ne. 0) return 

                write (unit=terrain_unit, iostat = io_status) nest%terrain
                if (io_status .ne. 0) return 

                close(terrain_unit)
            end if
        end if
        
    end subroutine dump_nest_info

    ! nest_point_to_coarse_point
    ! --------------------------
    ! Converts an i/j point on a COAMPS nest to the corresponding i/j
    ! point on the coarsest nest.
    !  PARAMETERS
    !   IN  domain            COAMPS domain that contains the nests
    !   IN  nest_pt           point to map onto the coarse grid
    !   OUT coarse_pt         point on the coarsest grid 
    recursive subroutine nest_point_to_coarse_point(nest_pt, coarse_pt)
        type(nest_point),    intent(in)  :: nest_pt
        type(nest_point),    intent(out) :: coarse_pt

        type(nest_point) :: parent_pt

        if (nest_pt%nest%id .eq. COARSE_NEST) then
            coarse_pt = nest_pt
        else
            parent_pt%nest => nest_pt%nest%parent_nest

            parent_pt%ii   = child_to_parent(nest_pt%ii, nest_pt%nest%anchor_i)
            parent_pt%jj   = child_to_parent(nest_pt%jj, nest_pt%nest%anchor_j)

            call nest_point_to_coarse_point(parent_pt, coarse_pt)
        end if
    end subroutine nest_point_to_coarse_point

    ! coarse_point_to_nest_point
    ! --------------------------
    ! Converts an i/j point on a coarser nest to an i/j point on the finest nest
    ! available at that particular point
    !  PARAMETERS
    !   IN  coarse_pt         Point on a coarse nest to try to improve         
    !   OUT finest_pt         The highest-res point possible at this level
    !   OUT in_nest           (optional) true if coarse_pt is a valid point 
    recursive subroutine coarse_point_to_nest_point(coarse_pt, finest_pt, &
                                                    in_nest)
        type(nest_point),    intent(in)  :: coarse_pt
        type(nest_point),    intent(out) :: finest_pt
        logical, optional,   intent(out) :: in_nest

        type(nest_point) :: child_pt
        type(nest_point) :: trial_pt
        logical          :: found_point

        integer :: finest_nest_level

        type(coamps_nest_iterator) :: child_iterator

        ! Assume that all child nests are within the boundaries of the parent
        ! nest, so we can stop searching this branch if the point is outside
        ! this nest
        if (.not. in_this_nest(coarse_pt)) then
            if (present(in_nest)) in_nest = .false.
            return
        else
            if (present(in_nest)) in_nest = .true.
        end if

        ! Base case: this nest has the finest point
        finest_pt = coarse_pt

        ! No child nests means that the current nest is the finest possible
        ! nest for this search path
        if (is_empty_list(coarse_pt%nest%child_nests)) return

        ! See if any of the child nests can do a better job at this point
        ! Use the built-in "nest level" as a measure of how deep we've reached
        finest_nest_level   = coarse_pt%nest%nest_level
        child_iterator      = get_iterator(coarse_pt%nest%child_nests)
        do while (has_next(child_iterator))

            child_pt%nest => get_next(child_iterator)

            child_pt%ii   = parent_to_child(coarse_pt%ii,         &
                                            child_pt%nest%anchor_i)
            child_pt%jj   = parent_to_child(coarse_pt%jj,         &
                                            child_pt%nest%anchor_j)

            call coarse_point_to_nest_point(child_pt, trial_pt, found_point)

            ! No need for further checking since the point's not in that nest
            if (.not. found_point) cycle

            if (trial_pt%nest%nest_level > finest_nest_level) then
                finest_nest_level = trial_pt%nest%nest_level
                finest_pt         = trial_pt
            end if
        end do
    end subroutine coarse_point_to_nest_point

    ! get_subnest_i_width
    ! -------------------
    ! Returns the i-direction width of the specified subnest
    !  PARAMETERS
    !   IN  nest            nest to query
    !   IN  subnest_num     The particular subnest to query
    function get_subnest_i_width(nest, subnest_num)
        type(coamps_nest), intent(in)  :: nest
        integer,           intent(in)  :: subnest_num
        integer                        :: get_subnest_i_width

        get_subnest_i_width = get_subnest_imaxf(nest, subnest_num) - &
                              get_subnest_iminf(nest, subnest_num) + 1       
    end function get_subnest_i_width

    ! get_subnest_j_width
    ! -------------------
    ! Returns the j-direction width of the specified subnest
    !  PARAMETERS
    !   IN  nest            nest to query
    !   IN  subnest_num     The particular subnest to query
    function get_subnest_j_width(nest, subnest_num)
        type(coamps_nest), intent(in)  :: nest
        integer,           intent(in)  :: subnest_num
        integer                        :: get_subnest_j_width

        get_subnest_j_width = get_subnest_jmaxf(nest, subnest_num) - &
                              get_subnest_jminf(nest, subnest_num) + 1       
    end function get_subnest_j_width

    ! get_nest_size
    ! -------------
    ! Returns the number of points contained in the given nest
    !  PARAMETERS
    !   IN  nest            nest to query
    function get_nest_size(nest)
        type(coamps_nest), intent(in)  :: nest
        integer                        :: get_nest_size

        get_nest_size = nest%pts_x * nest%pts_y
    end function get_nest_size
    
    ! get_nest_level
    ! -------------
    ! Returns the nest level
    !  PARAMETERS
    !   IN  nest            nest to query
    function get_nest_level(nest_pt)
        type(nest_point), intent(in)  :: nest_pt
        integer                        :: get_nest_level

        get_nest_level = nest_pt%nest%nest_level
    end function get_nest_level
    
    ! get_nest_i_width
    ! ----------------
    ! Returns the maximum i-coordinate for this nest
    !  PARAMETERS
    !   IN  nest            nest to query
    function get_nest_i_width(nest)
        type(coamps_nest), intent(in)  :: nest
        integer                        :: get_nest_i_width

        get_nest_i_width = nest%pts_x
    end function get_nest_i_width

    ! get_nest_j_width
    ! ------------
    ! Returns the maximum j-coordinate for this nest
    !  PARAMETERS
    !   IN  nest            nest to query
    function get_nest_j_width(nest)
        type(coamps_nest), intent(in)  :: nest
        integer                        :: get_nest_j_width

        get_nest_j_width = nest%pts_y
    end function get_nest_j_width
    
    ! get_nest_delta_x
    ! ----------------
    ! Returns the x grid spacing of the nest
    !  PARAMETERS
    !   IN  nest            nest to query
    function get_nest_delta_x(nest)
        type(coamps_nest), intent(in)  :: nest
        integer                        :: get_nest_delta_x

        get_nest_delta_x = nest%delta_x
    end function get_nest_delta_x
    
    ! get_nest_dely
    ! ----------------
    ! Returns the y grid spacing of the nest
    !  PARAMETERS
    !   IN  nest            nest to query
    function get_nest_delta_y(nest)
        type(coamps_nest), intent(in)  :: nest
        integer                        :: get_nest_delta_y

        get_nest_delta_y = nest%delta_y
    end function get_nest_delta_y

    ! nest_index_2d_to_1d
    ! -------------------
    ! Converts a 2d (i,j) coordinate representation to a 1-d version
    !  PARAMETERS
    !   IN  nest            Nest to base the conversion on
    !   IN  ii              2D i-coordinate
    !   IN  jj              2D j-coordinate
    function nest_index_2d_to_1d(nest, ii, jj)
        type(coamps_nest), intent(in)  :: nest
        integer,           intent(in)  :: ii
        integer,           intent(in)  :: jj
        integer                        :: nest_index_2d_to_1d

        nest_index_2d_to_1d = (jj - 1) * nest%pts_x + ii

    end function nest_index_2d_to_1d

    ! nest_index_1d_to_3d
    ! -------------------
    ! Converts a 1d index to a 3d (i,j,k) coordinate
    !  PARAMETERS
    !   IN  nest            Nest to base the conversion on
    !   IN  ii              2D i-coordinate
    !   IN  jj              2D j-coordinate
    function nest_index_1d_to_3d(nest, indx) result(ijk)
        type(coamps_nest), intent(in)  :: nest
        integer,           intent(in)  :: indx
        integer, dimension(3)          :: ijk

        integer :: nsize2D, indx2D

        nsize2D = get_nest_size(nest)

        ijk(3) = (indx-1)/nsize2D + 1
        indx2D = indx - (ijk(3)-1)*nsize2D 

        ijk(2) = (indx2D-1)/nest%pts_x + 1

        ijk(1) = indx2D - (ijk(2)-1)*nest%pts_x

    end function nest_index_1d_to_3d

    !------------------------------
    ! END PUBLIC ROUTINES
    !------------------------------

    !------------------------------
    ! BEGIN PRIVATE ROUTINES
    !------------------------------

    ! initialize_decomposition
    ! ------------------------
    ! Given a COAMPS nest the total number of non-I/O processors for
    ! that nest, allocate space for each processor's [ij](min|max)f 
    ! extent
    ! N.B. The nest is not defined as intent(out) here since that 
    !      could mean that the other components of the coamps_nest 
    !      are undefined after this routine call.  This is avoided 
    !      by using intent(inout).
    !  PARAMETERS
    ! INOUT nest              coamps_nest to allocate space in
    !   IN  totalprocs        total number of non-I/O processors
    subroutine initialize_decomposition(nest, totalprocs)
        type(coamps_nest), intent(inout) :: nest
        integer, intent(in)              :: totalprocs

        character(len=*), parameter :: routine = 'initialize_decomposition'
        integer :: alloc_status

        nullify(nest%iminf)
        nullify(nest%imaxf)
        nullify(nest%jminf)
        nullify(nest%jmaxf)

        allocate( nest%iminf(totalprocs), stat=alloc_status )
        call check_alloc_status(alloc_status, routine, source, revision, &
                                revdate, 'nest%iminf')
        allocate( nest%imaxf(totalprocs), stat=alloc_status )
        call check_alloc_status(alloc_status, routine, source, revision, &
                                revdate, 'nest%imaxf')
        allocate( nest%jminf(totalprocs), stat=alloc_status )
        call check_alloc_status(alloc_status, routine, source, revision, &
                                revdate, 'nest%jminf')
        allocate( nest%jmaxf(totalprocs), stat=alloc_status )
        call check_alloc_status(alloc_status, routine, source, revision, &
                                revdate, 'nest%jmaxf')
    end subroutine initialize_decomposition

    ! generate_terrht_filename
    ! ------------------------
    ! Generates the COAMPS terrain height flat file name for a particular
    ! nest of a grid at a given date-time group.  This does *not* 
    ! generate any path information - it only returns the file name.
    !
    ! Assumes that the coamps_nest structure has already been populated
    ! with the domain information, since we need the number of x and y
    ! grid points to generate the file name
    !  PARAMETERS
    !   IN  nest              nest that we need terrain info for
    !   IN  dtg               base date-time group for model run 
    !   OUT terrht_filename   name of terrain height flat file
    subroutine generate_terrht_filename(nest, dtg, terrht_filename)
        type(coamps_nest),   intent(in)  :: nest
        character(len=10),   intent(in)  :: dtg
        character(len=64),   intent(out) :: terrht_filename

        ! The format of the terrain height file is fixed except for
        ! the horizontal grid size and date-time group
        call generate_flat_file_name( var_name   = 'terrht',      &
                                      level_type = 'sfc',         &
                                      level1     = 0,             &
                                      level2     = 0,             &
                                      gridnum    = nest%id,       &
                                      aoflag     = 'a',           &
                                      xpts       = nest%pts_x,    &
                                      ypts       = nest%pts_y,    &
                                      dtg        = dtg,           &
                                      tau_hh     = 0,             &
                                      tau_mm     = 0,             &
                                      tau_ss     = 0,             &
                                      field_type = 'fcstfld',     &
                                      file_name  = terrht_filename )
    end subroutine generate_terrht_filename

    ! read_terrain_height
    ! -------------------
    ! Reads the terrain height for a given nest and date-time group
    !  PARAMETERS
    ! INOUT nest              The nest to add terrain data to
    !   IN  dtg               COAMPS date-time group
    subroutine read_terrain_height(nest, dtg)
        type(coamps_nest), intent(inout) :: nest
        character(len=10), intent(in)    :: dtg

        character(len=64) :: terrht_name
        integer           :: terrht_unit
        integer           :: terrht_record_len
        integer           :: r4_len

        character(len=*), parameter :: routine = 'read_terrain_height'
        integer :: io_status
        integer :: alloc_status
        integer :: dealloc_status

        ! Terrain data is stored in a flat file, so it uses the COAMPS real
        ! number size.  Storing it as a single dimension instead of 2-D makes
        ! (possible) byteswapping easier 
        real(kind=r8), dimension(:), pointer :: temp_terrain 

        ! Temporary storage
        allocate(temp_terrain(nest%pts_x * nest%pts_y), stat=alloc_status)
        call check_alloc_status(alloc_status, routine, source, revision, &
                                revdate, 'temporary terrain')

        ! Permanent storage
        nullify(nest%terrain)
        allocate(nest%terrain(nest%pts_x, nest%pts_y), stat = alloc_status)
        call check_alloc_status(alloc_status, routine, source, revision,  &
                                revdate, 'nest terrain')

        inquire(IOLENGTH=r4_len) 0_r4

        terrht_record_len = nest%pts_x * nest%pts_y * r4_len

        terrht_unit = get_unit()
        call generate_terrht_filename(nest, dtg, terrht_name)
        open(unit = terrht_unit, file = terrht_name, status = 'old',    & 
             access = 'direct', action = 'read', form = 'unformatted',  &
             recl = terrht_record_len, iostat = io_status)  
        call check_io_status(io_status, routine, source, revision, &
                             revdate, 'Opening terrain height file')

        call read_flat_file(terrht_unit, temp_terrain)
        close(terrht_unit)

        nest%terrain = reshape(real(temp_terrain, kind=r8), (/ nest%pts_x, nest%pts_y /))

        deallocate(temp_terrain, stat=dealloc_status)
        call check_dealloc_status(dealloc_status, routine, source, revision, &
                                  revdate, 'temporary terrain')
    end subroutine read_terrain_height

    ! initialize_nest_latlon
    ! -------------------
    ! Initialize the lat/lon arrays for the nest
    ! INOUT nest              The nest to add terrain data to
    !   IN  dtg               COAMPS date-time group
    subroutine initialize_nest_latlon(nest, grid)
        type(coamps_nest), intent(inout) :: nest
        type(coamps_grid), intent(in)    :: grid

        type(nest_point) :: coarse_point
        real(kind=r8)    :: lat, lon
        integer          :: ii, jj

        character(len=*), parameter :: routine = 'initialize_nest_latlon'
        integer :: alloc_status

        nullify(nest%lat)
        allocate(nest%lat(nest%pts_x, nest%pts_y), stat = alloc_status)
        call check_alloc_status(alloc_status, routine, source, revision, revdate, 'lat')

        nullify(nest%lon)
        allocate(nest%lon(nest%pts_x, nest%pts_y), stat = alloc_status)
        call check_alloc_status(alloc_status, routine, source, revision, revdate, 'lon')

        do ii=1,nest%pts_x
          do jj=1,nest%pts_y
            call nest_point_to_coarse_point(make_nest_point(nest, ii, jj), coarse_point)
            call grid_point_to_latlon(grid, coarse_point%ii, coarse_point%jj,          &
                                      nest%lat(ii,jj), nest%lon(ii,jj))
          end do
        end do

    end subroutine initialize_nest_latlon

    ! get_nest_latlon
    ! -------------------
    ! Initialize the lat/lon arrays for the nest
    ! IN nest              The nest
    ! IN ii, jj            ii, jj points
    subroutine get_nest_latlon(nest, ii, jj, lat, lon)
        type(coamps_nest), intent(in)  :: nest
        integer,           intent(in)  :: ii
        integer,           intent(in)  :: jj
        real(kind=r8),     intent(out) :: lat
        real(kind=r8),     intent(out) :: lon

        lat = nest%lat(ii,jj)
        lon = nest%lon(ii,jj)

    end subroutine get_nest_latlon

    ! LINKED LIST STUFF

    ! initialize_coamps_nest_list
    ! ---------------------------
    ! Constructor for coamps_nest linked list
    !  PARAMETERS
    ! INOUT list    List to initialize
    subroutine initialize_coamps_nest_list(list)
        type(coamps_nest_list), intent(inout) :: list

        nullify(list%head)
        nullify(list%tail)
    end subroutine initialize_coamps_nest_list

    ! add_coamps_nest_to_list
    ! -------------------
    ! Adds the given coamps_nest to a given linked list
    !  PARAMETERS
    !   IN  new_entry       Integer to add to the list
    ! INOUT list            List to add the entry to
    subroutine add_coamps_nest_to_list(new_entry, list)
        type(coamps_nest), target, intent(in)    :: new_entry
        type(coamps_nest_list),    intent(inout) :: list

        ! Add a new node
        if (.not. associated(list%head)) then
            allocate(list%head)
            list%tail => list%head
        else
            allocate(list%tail%next)
            list%tail => list%tail%next
        end if

        ! Assign the contents
        nullify(list%tail%next)
        list%tail%value => new_entry
    end subroutine add_coamps_nest_to_list

    ! print_coamps_nest_list
    ! ------------------
    ! Traverses the linked list and prints the nest ID of each element
    !  PARAMETERS
    !   IN  list            linked list to iterate over
    subroutine print_coamps_nest_list(list)
        type(coamps_nest_list), intent(in) :: list

        type(coamps_nest_iterator) :: iterator

        iterator = get_iterator(list)

        do while (has_next(iterator))
            write (*,*) get_nest_id(get_next(iterator))
        end do
    end subroutine print_coamps_nest_list

    ! is_coamps_nest_list_empty
    ! -------------------------
    ! Returns true if the given linked list contains no elements
    !  PARAMETERS
    !   IN  list            linked list to query
    function is_coamps_nest_list_empty(list)
        type(coamps_nest_list), intent(in)  :: list
        logical                             :: is_coamps_nest_list_empty

        is_coamps_nest_list_empty = .not. associated(list%head)
    end function is_coamps_nest_list_empty

    ! get_coamps_nest_list_iterator
    ! -----------------------------
    ! Returns an iterator for the given linked list
    !  PARAMETERS
    !   IN  list            list to iterate over
    function get_coamps_nest_list_iterator(list)
        type(coamps_nest_list),      intent(in) :: list
        type(coamps_nest_iterator)              :: get_coamps_nest_list_iterator

        get_coamps_nest_list_iterator%current => list%head
    end function get_coamps_nest_list_iterator

    ! coamps_nest_iterator_has_next
    ! -----------------------------
    ! Returns true if the iterator can return another value
    !  PARAMETERS
    !   IN  iterator        the iterator to check the state of
    function coamps_nest_iterator_has_next(iterator)
        type(coamps_nest_iterator), intent(in) :: iterator
        logical                            :: coamps_nest_iterator_has_next

        if (associated(iterator%current)) then
            coamps_nest_iterator_has_next = .true.
        else
            coamps_nest_iterator_has_next = .false.
        end if
    end function coamps_nest_iterator_has_next

    ! coamps_nest_iterator_get_next
    ! ------------------------------
    ! Returns the current value at the iterator location, ***and advances the
    ! iterator to the nest entry***
    !  PARAMETERS
    !   IN  iterator        List iterator to get value of and advance
    function coamps_nest_iterator_get_next(iterator)
        type(coamps_nest_iterator), intent(inout) :: iterator
        type(coamps_nest), pointer                :: &
                coamps_nest_iterator_get_next

        coamps_nest_iterator_get_next => iterator%current%value
        iterator%current               => iterator%current%next  
    end function coamps_nest_iterator_get_next


    !------------------------------
    ! END PRIVATE ROUTINES
    !------------------------------

end module coamps_nest_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
