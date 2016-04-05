! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

!------------------------------
! MODULE:       coamps_netcdf_mod
! AUTHOR:       T. R. Whitcomb
!               Naval Research Laboratory
! DART VERSION: Jamaica
!
! Write COAMPS state vector to a NetCDF file
!------------------------------ 

module coamps_netcdf_mod

    use netcdf

    use location_mod,      only : location_type, get_location

    use location_mod,      only : VERTISUNDEF, VERTISSURFACE, VERTISLEVEL,   &
                                  VERTISPRESSURE, VERTISHEIGHT

    use types_mod,         only : r8

    use utilities_mod,     only : nc_check,                     &
                                  do_output,                    &
                                  error_handler,                &
                                  E_ERR

    use coamps_util_mod,   only : check_alloc_status,           &
                                  check_dealloc_status

    use coamps_statevec_mod, only : state_vector,               &
                                  get_domain,                   &
                                  state_iterator,               &
                                  get_iterator,                 &
                                  has_next,                     &
                                  get_next

    use coamps_statevar_mod, only : state_variable,             &
                                  get_var_substate,             &
                                  get_var_name,                 &
                                  get_nc_varid,                 &
                                  set_nc_varid,                 &
                                  get_vert_value,               &
                                  get_vert_type,                &
                                  get_var_stagger,              &
                                  get_var_nlevs,                &
                                  get_nest_number,              &
                                  dump_state_variable

    use coamps_domain_mod, only : coamps_domain,                &
                                  nest_point_to_latlon,         &
                                  get_domain_num_levels,        &
                                  get_domain_msigma,            &
                                  get_domain_wsigma,            &
                                  get_nest_count,               &
                                  get_domain_nest

    use coamps_nest_mod,   only : coamps_nest,                  &
                                  make_nest_point,              &
                                  get_nest_i_width,             &
                                  get_nest_j_width

    implicit none

    private

    !------------------------------
    ! BEGIN PUBLIC INTERFACE
    !------------------------------

    public :: nc_write_statearray_atts
    public :: nc_write_prognostic_atts
    public :: nc_write_statearray_data
    public :: nc_write_prognostic_data

    !------------------------------
    ! END PUBLIC INTERFACE
    !------------------------------

    !------------------------------
    ! BEGIN EXTERNAL INTERFACES
    !------------------------------

    interface populate_variable
        module procedure populate_integer_variable, populate_real_variable
    end interface populate_variable

    !------------------------------
    ! END EXTERNAL INTERFACES
    !------------------------------

    !------------------------------
    ! BEGIN TYPES AND CONSTANTS 
    !------------------------------

    type :: coordinate_var
        private

        integer          :: dim_id
        integer          :: var_id
        integer          :: nest_number
        integer          :: dim_length
        integer          :: coordinate_direction 
        real(kind=r8)    :: stagger
    end type coordinate_var

    integer, parameter :: XDIR_COORD = 1
    integer, parameter :: YDIR_COORD = 2
    integer, parameter :: ZDIR_COORD = 3
    integer, parameter :: NULL_COORD = 0
    integer, parameter :: LON_DIM    = 1
    integer, parameter :: LAT_DIM    = 2

    real(kind=r8), dimension(2) :: LAT_LIMITS   = (/ -90.0_r8, 90.0_r8 /)
    real(kind=r8), dimension(2) :: LON_LIMITS   = (/ 0.0_r8, 360.0_r8  /)

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

    character(len=128) :: msgstring   ! general purpose string for printing

    logical, save :: module_initialized = .false.

    !------------------------------
    ! END MODULE VARIABLES
    !------------------------------

contains
    !------------------------------
    ! BEGIN PUBLIC ROUTINES
    !------------------------------

    ! nc_write_statearray_data
    ! ----------------------------
    ! Writes the model state vector to the supplied location in the supplied
    ! NetCDF file
    !  PARAMETERS
    !   IN  ncFileID      NetCDF file to write to
    !   IN  state           model state vector
    !   IN  copy_index      Which copy (i.e. member) this is 
    !   IN  time_index      Which time this is
    subroutine nc_write_statearray_data(ncFileID, state, copy_index, time_index)
        integer, intent(in)                     :: ncFileID
        real(kind=r8), dimension(:), intent(in) :: state
        integer, intent(in)                     :: copy_index
        integer, intent(in)                     :: time_index

        integer :: state_var_id

        ! NetCDF variables                                               
        integer :: n_dims, n_vars, n_atts  ! NetCDF counts               
        integer :: ulim_dimid              ! unlimited dimension         

        ! ------------------------------------------------------------------
        ! Error handling
        ! ------------------------------------------------------------------
        character(len=*), parameter :: routine = 'nc_write_statearray_data'

        ! ------------------------------------------------------------------
        ! Make sure we have the most recent additions to the file before adding more
        ! ------------------------------------------------------------------
        call nc_check(nf90_inquire(ncFileID, n_dims, n_vars, n_atts, ulim_dimid), routine)
        call nc_check(nf90_sync(ncFileID), routine) 

        ! ------------------------------------------------------------------
        ! Get the numerical index of the state variable, then dump the
        ! state vector to the corresponding ensemble member and time
        ! ------------------------------------------------------------------
        call nc_check(nf90_inq_varid(ncFileID, 'state', state_var_id), routine) 
        call nc_check(nf90_put_var(ncFileID, state_var_id, state,       &
                      start=(/ 1, copy_index, time_index /)), routine)
     
        ! ------------------------------------------------------------------
        ! Flush the buffer to disk                                       
        ! ------------------------------------------------------------------
        call nc_check( nf90_sync(ncFileID), routine) 
        if (do_output()) write (*,*) 'Model data written and file synced...' 

    end subroutine nc_write_statearray_data

    ! nc_write_prognostic_data
    ! ----------------------------
    ! Writes the models prognostic variables to the supplied NetCDF file
    !  PARAMETERS
    !   IN  ncFileID      NetCDF file to write to
    !   IN  state           model state vector
    !   IN  copy_index      Which copy (i.e. member) this is 
    !   IN  time_index      Which time this is
    subroutine nc_write_prognostic_data(ncFileID, state_list, statevec, copy_index, time_index)
       integer,                     intent(in)      :: ncFileID
       type(state_vector),          intent(in)      :: state_list
       real(kind=r8), dimension(:), intent(in)      :: statevec
       integer,                     intent(in)      :: copy_index
       integer,                     intent(in)      :: time_index
       integer                                      :: ierr   

       type(state_iterator)                         :: iterator
       type(state_variable)                         :: cur_var
       type(coamps_nest)                            :: nest
       type(coamps_domain)                          :: domain
       integer, dimension(3)                        :: var_dims
       character(len=NF90_MAX_NAME)                 :: var_name
       integer                                      :: alloc_status
       integer                                      :: varid 
       integer                                      :: nx, ny

       real(kind=r8), allocatable, dimension(:,:,:) :: var3d
       real(kind=r8), allocatable, dimension(:,:)   :: var2d
       real(kind=r8), pointer,     dimension(:)     :: var_substate

       character(len=*), parameter :: routine = 'nc_write_prognostic_data'

       domain   = get_domain(state_list)
       iterator = get_iterator(state_list)
       output_vars:  do while(has_next(iterator))

          cur_var  = get_next(iterator)
          var_name = get_var_name(cur_var)
          var_dims = get_var_dims(cur_var, domain)
          varid    = get_nc_varid(cur_var)

          nest = get_domain_nest(domain, get_nest_number(cur_var))
          nx   = get_nest_i_width(nest)
          ny   = get_nest_j_width(nest)

          var_substate => get_var_substate(cur_var, statevec)

          select case(get_var_rank(cur_var, domain))
             case(2) ! Two dimensional variables         

               allocate(var2d(nx,ny), stat=alloc_status)
               call check_alloc_status(alloc_status, routine, source, revision,   &
                                       revdate, '2D allocate '//trim(var_name))

               var2d = reshape(var_substate, (/nx, ny/))
               call nc_check(nf90_put_var(ncFileID, varid, var2d,                 &
                             start=(/ 1, 1, copy_index, time_index /)),           &
                             routine, 'put_var '//trim(var_name))
             
               deallocate(var2d, stat=alloc_status)
               call check_dealloc_status(alloc_status, routine, source, revision, &
                                         revdate, '2D deallocate '//trim(var_name))

             case(3)  ! Three dimensional variables

               allocate(var3d(nx, ny, var_dims(3)), stat=alloc_status)
               call check_alloc_status(alloc_status, routine, source, revision,   &
                                       revdate, '3D allocate '//trim(var_name))

               var3d = reshape(var_substate, (/nx, ny, var_dims(3)/))
               call nc_check(nf90_put_var(ncFileID, varid,                 &
                             var3d(1:var_dims(1), 1:var_dims(2), :),       &
                             start=(/ 1, 1, 1, copy_index, time_index /)), &
                             routine, 'put_var '//trim(var_name))

             deallocate(var3d, stat=alloc_status)
             call check_dealloc_status(alloc_status, routine, source, revision, &
                              revdate, '3D deallocate '//trim(var_name))

           case default

             write(msgstring,*) 'not built to handle ndims == ', get_var_rank(cur_var, domain)
             call error_handler(E_ERR, routine, msgstring, source, revision, revdate) 
          end select

      end do output_vars

    end subroutine nc_write_prognostic_data

    ! nc_write_statearray_atts
    ! --------------------------
    ! Writes the model attributes (metadata, etc.) to a NetCDF file
    !  PARAMETERS
    !   IN  ncFileID      NetCDF file to write to
    !   IN  model_size      Length of the state vector
    subroutine nc_write_statearray_atts(ncFileID, locations, types)
        integer,                           intent(in) :: ncFileID
        type(location_type), dimension(:), intent(in) :: locations
        integer,             dimension(:), intent(in) :: types
               
        type(coordinate_var) :: state_coord
        type(coordinate_var) :: lat_coord
        type(coordinate_var) :: lon_coord
        type(coordinate_var) :: vert_coord
        integer              :: state_var_id
        integer              :: type_var_id

        integer :: model_size

        ! NetCDF variables                                               
        integer :: n_dims, n_vars, n_atts  ! NetCDF counts               
        integer :: ulim_dimid              ! unlimited dimension         

        character(len=*), parameter :: routine = 'nc_write_state_array_atts'

        model_size = size(types)

        ! ------------------------------------------------------------------
        ! Make sure we have the most recent additions to the file before adding more
        ! ------------------------------------------------------------------
        call nc_check(nf90_inquire(ncFileID, n_dims, n_vars, n_atts, ulim_dimid), routine)
        call nc_check(nf90_sync(ncFileID), routine) 

        ! ------------------------------------------------------------------
        ! Enter into define mode
        ! ------------------------------------------------------------------
        call nc_check(nf90_redef(ncFileID), routine) 
     
        call set_global_attributes(ncFileID)

        lat_coord   = new_coordinate_var(ncFileID, 'lat',            model_size, nf90_double, &
                      'Latitude',     'degrees_north', LAT_LIMITS, cv_dir = XDIR_COORD)
        lon_coord   = new_coordinate_var(ncFileID, 'lon',            model_size, nf90_double, &
                      'Longitude',    'degrees_east', LON_LIMITS, cv_dir = YDIR_COORD)
        vert_coord  = new_coordinate_var(ncFileID, 'height',         model_size, nf90_double, &
                      'Sigma Height', 'meters', cv_dir = ZDIR_COORD)
        state_coord = new_coordinate_var(ncFileID, 'StateVariable',  model_size, nf90_int,    &
                      'State variable Coordinate', 'index', (/1.0_r8,                         &
					  real(model_size,kind=r8)/), cv_dir = NULL_COORD)

        ! Data variables
        state_var_id = new_variable(ncFileID, 'state', NF90_DOUBLE,   &
                                    (/ get_dim_id(state_coord),       &
                                       get_copy_dim_id(ncFileID),     &
                                       get_time_dim_id(ncFileID) /),  &
                                    var_long_name = 'Model State') 

        type_var_id = new_variable(ncFileID, 'type', NF90_INT,       &
                                   (/ get_dim_id(state_coord) /),      &
                                   var_long_name = 'Variable Type')

        !-------------------------------------------------------------------
        ! Need to get out of define mode to fill the variables 
        !-------------------------------------------------------------------
        call nc_check(nf90_enddef(ncFileID), routine) 

        ! ------------------------------------------------------------------
        ! Fill the coordinate variables
        ! ------------------------------------------------------------------
        call write_state_coordinate_data(ncFileID, state_coord, lat_coord, &
                                         lon_coord, vert_coord, locations)

        ! ------------------------------------------------------------------
        ! Fill the observation types (kinds?)
        ! ------------------------------------------------------------------
        call nc_check(nf90_put_var(ncFileID, type_var_id, types), routine)
     
        ! ------------------------------------------------------------------
        ! Sync the NetCDF file
        ! ------------------------------------------------------------------
        call nc_check(nf90_sync(ncFileID), routine) 
        if (do_output()) write (*,*) 'Model attributes written and file synced...' 

    end subroutine nc_write_statearray_atts

    ! nc_write_prognostic_atts
    ! -------------------
    ! Write model-specific global attributes to a NetCDF file.
    !  PARAMETERS
    !   IN  ncFileID          numeric ID of an *open* NetCDF file
    !   OUT ierr              0 if writing was successful
    !
    ! If the dimension of the domain is specified to by nx, ny, nx ...
    ! then var_stagger is related to var_dims in the following way: 
    !
    ! T-stagger has dimensions (nx,ny,nz); 
    ! W-stagger has dimensions (nx,ny,nz+1); 
    ! U-stagger has dimensions (nx-1,ny,nz); and 
    ! V-stagger has dimensions (nx,ny-1,nz) 
    ! (COAMPS is on a C-grid with p-u-p ordering in the horizontal).
 
    subroutine nc_write_prognostic_atts( ncFileID, state_list)
      integer,            intent(in)      :: ncFileID      ! netCDF file 
      type(state_vector), intent(in)      :: state_list

      integer                        :: nx, ny, nz
      integer                        :: ii, jj, n
      integer                        :: nest_count
      integer                        :: ncoord_total
      integer                        :: ncoord
      integer                        :: ndims
      integer, dimension(3)          :: dimids
      type(state_iterator)           :: iterator
      type(state_variable), pointer  :: cur_var
      type(coamps_domain)            :: domain

      type(coordinate_var), allocatable, dimension(:) :: coords

      integer,              dimension(2) :: varidT
      integer,              dimension(2) :: varidU
      integer,              dimension(2) :: varidV
      type(coordinate_var), dimension(2) :: latlon_coord

      ! Error handling
      character(len=*), parameter  :: routine = 'nc_write_prognostic_atts'
      integer                      :: alloc_status
      integer                      :: dealloc_status
 
      character(len=NF90_MAX_NAME) :: coord_str, var_name 

      integer, parameter :: X_COORD_INDEX   = 0
      integer, parameter :: XM1_COORD_INDEX = 1
      integer, parameter :: Y_COORD_INDEX   = 2
      integer, parameter :: YM1_COORD_INDEX = 3


      ! Make sure we have the most recent additions to the file before adding more
      call nc_check(nf90_sync(ncFileID), routine) 
      call nc_check(nf90_redef(ncFileID), routine) 

      domain = get_domain(state_list)

      ! Total number of coordinate variables 
      nest_count = get_nest_count(domain)
      ncoord_total = 4 * nest_count + 2

      allocate(coords(ncoord_total), stat=alloc_status)
      call check_alloc_status(alloc_status, routine, source, revision, revdate, 'coords')

      nz = get_domain_num_levels(domain)

400 format(A,'_g',I2.2)
402 format(A,'_',A,I5.5,'_g',I2.2)

      define_coord_vars: do n=1,nest_count
        ncoord = 1 + 4*(n-1)
        nx   = get_nest_i_width(get_domain_nest(domain, n))
        ny   = get_nest_j_width(get_domain_nest(domain, n))

        write(coord_str, 400) 'nx', n
        coords(ncoord + X_COORD_INDEX)   =                                         &
                   new_coordinate_var(ncFileID, trim(coord_str), nx, nf90_double,  &
                   'X Grid Points', '', cv_stagger = 0.0_r8, cv_dir = XDIR_COORD,  &
                   cv_nest = n)

        write(coord_str, 400) 'nxm1', n
        coords(ncoord + XM1_COORD_INDEX)   =                                                &
                   new_coordinate_var(ncFileID,   trim(coord_str),   nx-1, nf90_double,     &
                   'Staggered X Grid Points', '', cv_stagger = 0.5_r8, cv_dir = XDIR_COORD, &
                   cv_nest = n)

        write(coord_str, 400) 'ny', n
        coords(ncoord + Y_COORD_INDEX)   =                                          &
                   new_coordinate_var(ncFileID, trim(coord_str), ny, nf90_double,     &
                   'Y Grid Points', '', cv_stagger = 0.0_r8, cv_dir = YDIR_COORD,   &
                   cv_nest = n)

        write(coord_str, 400) 'nym1', n
        coords(ncoord + YM1_COORD_INDEX)   =                                                &
                   new_coordinate_var(ncFileID, trim(coord_str), ny-1, nf90_double,         &
                   'Staggered Y Grid Points', '', cv_stagger = 0.5_r8, cv_dir = YDIR_COORD, &
                   cv_nest = n)

      end do define_coord_vars

      coords(ncoord_total-1) = &
                   new_coordinate_var(ncFileID, 'sigw', nz+1, nf90_double, &
                   'W Level Sigma', 'meters', cv_stagger = 0.0_r8, cv_dir = ZDIR_COORD)

      coords(ncoord_total)   = &
                   new_coordinate_var(ncFileID, 'sigm',   nz, nf90_double, & 
                   'Mass Level Sigma', 'meters',  cv_stagger = 0.5_r8, cv_dir = ZDIR_COORD)
      
      ! ------------------------------------------------------------------
      ! Loop over all the variables in the state vector
      ! ------------------------------------------------------------------
      iterator = get_iterator(state_list)
      define_atts:  do while(has_next(iterator))

        cur_var => get_next(iterator)

        select case (get_vert_type(cur_var))
          case(VERTISHEIGHT)
            write(var_name, 402) trim(get_var_name(cur_var)), 'Z', &
                                 int(get_vert_value(cur_var)), get_nest_number(cur_var)
          case(VERTISPRESSURE)
            write(var_name, 402) trim(get_var_name(cur_var)), 'P', &
                                 int(get_vert_value(cur_var)), get_nest_number(cur_var)
          case(VERTISSURFACE, VERTISUNDEF)
            write(var_name, 402) trim(get_var_name(cur_var)), 'S', &
                                 get_vert_value(cur_var), get_nest_number(cur_var)
          case default
            write(var_name, 400) trim(get_var_name(cur_var)), get_nest_number(cur_var)
        end select

        ndims  = get_var_rank(cur_var, domain)

        dimids = get_dimids(ndims        = ndims,                         &
                            var_dims     = get_var_dims(cur_var, domain), & 
                            nest_number  = get_nest_number(cur_var),      & 
                            local_coords = coords) 

        call set_nc_varid(cur_var,                             &
                  new_variable(                                &
                  ncFileID    = ncFileID,                      &
                  var_name    = var_name,                      &
                  var_type    = nf90_double,                   &
                  dim_ids     = (/ dimids(1:ndims),            &
                                 get_copy_dim_id(ncFileID),    &
                                 get_time_dim_id(ncFileID) /), &
                  var_stagger = get_var_stagger(cur_var)))

      end do define_atts

      ! ------------------------------------------------------------------
      ! Need to get out of define mode to fill the coordinate variables 
      ! ------------------------------------------------------------------
      call nc_check(nf90_enddef(ncFileID), routine//': nf_enddef') 

      fill_coords: do ii=1,ncoord_total-2
        call nc_check(nf90_put_var(ncFileID, get_var_id(coords(ii)),       &
                     (/ (real(jj-1,kind=r8)+get_coord_stagger(coords(ii)), &
                        jj=1,get_dim_length(coords(ii))) /) ), routine) 
      end do fill_coords

      call nc_check(nf90_put_var(ncFileID, get_var_id(coords(ncoord_total)),  & 
                    get_domain_msigma(domain)), routine//': put_var msigma') 

      call nc_check(nf90_put_var(ncFileID, get_var_id(coords(ncoord_total-1)),  & 
                    get_domain_wsigma(domain)), routine//': put_var wsigma') 


      ! ------------------------------------------------------------------
      ! Define and fill the lat/lon grid at mass and momentum points
      ! ------------------------------------------------------------------
      fill_latlon: do n=1,nest_count
        ncoord = 1 + 4*(n-1)

        nx   = get_nest_i_width(get_domain_nest(domain, n))
        ny   = get_nest_j_width(get_domain_nest(domain, n))

        call nc_check(nf90_redef(ncFileID), routine//': nf_redef(fill_latlon)') 
        varidT = define_latlon_vars(ncFileID, 'T', (/nx,   ny/),   n, coords)
        varidU = define_latlon_vars(ncFileID, 'U', (/nx-1, ny/),   n, coords)
        varidV = define_latlon_vars(ncFileID, 'V', (/nx,   ny-1/), n, coords)

        call nc_check(nf90_enddef(ncFileID), routine//': nf_enddef(fill_latlon)') 
        latlon_coord(LAT_DIM) = coords(ncoord + Y_COORD_INDEX)
        latlon_coord(LON_DIM) = coords(ncoord + X_COORD_INDEX)
        call fill_nc_latlon(ncFileID, domain, n, varidT, latlon_coord)

        latlon_coord(LAT_DIM) = coords(ncoord + Y_COORD_INDEX)
        latlon_coord(LON_DIM) = coords(ncoord + XM1_COORD_INDEX)
        call fill_nc_latlon(ncFileID, domain, n, varidU, latlon_coord)

        latlon_coord(LAT_DIM) = coords(ncoord + YM1_COORD_INDEX)
        latlon_coord(LON_DIM) = coords(ncoord + X_COORD_INDEX)
        call fill_nc_latlon(ncFileID, domain, n, varidV, latlon_coord)

      end do fill_latlon

      deallocate(coords, stat=dealloc_status)
      call check_dealloc_status(dealloc_status, routine, source, revision,   &
                              revdate, 'coords')

    end subroutine nc_write_prognostic_atts

    !------------------------------
    ! END PUBLIC ROUTINES
    !------------------------------

    !------------------------------
    ! BEGIN PRIVATE ROUTINES
    !------------------------------

    ! set_global_attributes
    ! ---------------------
    ! Set the single-use attributes common to the entire model in the NetCDF
    ! file
    subroutine set_global_attributes(ncFileID)
        integer, intent(in) :: ncFileID

        integer, parameter :: DT_YEAR   = 1
        integer, parameter :: DT_MONTH  = 2
        integer, parameter :: DT_DAY    = 3
        integer, parameter :: DT_HOUR   = 5
        integer, parameter :: DT_MINUTE = 6
        integer, parameter :: DT_SECOND = 7

        integer, dimension(8)        :: dt_values 
        character(len=NF90_MAX_NAME) :: dt_string 

        character(len=*), parameter :: routine = 'set_global_attributes'

        call DATE_AND_TIME(values=dt_values) 
        write (dt_string, '(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))')   & 
             dt_values(DT_YEAR), dt_values(DT_MONTH), dt_values(DT_DAY),  &
             dt_values(DT_HOUR), dt_values(DT_MINUTE), dt_values(DT_SECOND) 
     
        call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'creation_date', dt_string), routine) 
        call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_source', source),     routine) 
        call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revision', revision), routine)  
        call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revdate', revdate),   routine) 
        call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model', 'COAMPS'),          routine)

    end subroutine set_global_attributes

    ! new_dimension
    ! -------------
    ! Define a new dimension in the NetCDF file
    function new_dimension(ncFileID, dim_name, length)
        integer,          intent(in)  :: ncFileID
        character(len=*), intent(in)  :: dim_name
        integer,          intent(in)  :: length
        integer                       :: new_dimension
        
        character(len=*), parameter :: routine = 'new_dimension'

        call nc_check(nf90_def_dim(ncid=ncFileID, name=dim_name,   &
                                   len=length, dimid=new_dimension), routine)  
    end function new_dimension

    ! new_variable
    ! ------------
    ! Defines a variable in the NetCDF file
    function new_variable(ncFileID, var_name, var_type, dim_ids, &
                          var_stagger, var_units, var_long_name, var_range)
        integer,                               intent(in)  :: ncFileID
        character(len=*),                      intent(in)  :: var_name
        integer,                               intent(in)  :: var_type
        integer, dimension(:),                 intent(in)  :: dim_ids
        character(len=*), optional,            intent(in)  :: var_long_name
        character(len=*), optional,            intent(in)  :: var_units
        character(len=*), optional,            intent(in)  :: var_stagger
        real(kind=r8), dimension(:), optional, intent(in)  :: var_range
        integer                                            :: new_variable

        character(len=*), parameter :: routine = 'new_variable'

        call nc_check(nf90_def_var(ncFileID, trim(var_name), &
                      var_type, dim_ids, new_variable), routine//' '//trim(var_name)) 

        if (present(var_long_name)) then
            call nc_check(nf90_put_att(ncFileID, new_variable, 'long_name', &
                          trim(var_long_name)), routine//' '//trim(var_long_name)) 
        end if

        if (present(var_units)) then
            call nc_check(nf90_put_att(ncFileID, new_variable, 'units', &
                          trim(var_units)), routine//': units')
        end if

        if (present(var_range)) then
            call nc_check(nf90_put_att(ncFileID, new_variable, &
                          'valid_range',var_range), routine//': range')
        end if

        if (present(var_stagger)) then
            call nc_check(nf90_put_att(ncFileID, new_variable, &
                          'staggering',var_stagger), routine//': stagger')
        end if

    end function new_variable

    ! new_coordinate_var
    ! ------------------
    ! Creates a new "coordinate variable" - in the NetCDF parlance, this is 
    ! a variable and a dimension with the same name.  Note that this function
    ! Requires that a long name and units be supplied
    function new_coordinate_var(ncFileID, cv_name, cv_length, cv_type, &
                                cv_long_name, cv_units, cv_range, cv_dir, cv_stagger, &
                                cv_nest)
        integer,                               intent(in)  :: ncFileID
        character(len=*),                      intent(in)  :: cv_name
        integer,                               intent(in)  :: cv_length
        integer,                               intent(in)  :: cv_type
        character(len=*),                      intent(in)  :: cv_long_name
        character(len=*),                      intent(in)  :: cv_units
        real(kind=r8), dimension(:), optional, intent(in)  :: cv_range
        integer,                     optional, intent(in)  :: cv_dir
        real(kind=r8),               optional, intent(in)  :: cv_stagger
        integer,                     optional, intent(in)  :: cv_nest
        type(coordinate_var)                               :: new_coordinate_var


        call set_dim_id(new_coordinate_var, new_dimension(ncFileID, cv_name, cv_length))

        call set_dim_length(new_coordinate_var, cv_length)

        if (present(cv_range)) then
            call set_var_id(new_coordinate_var,                           &
                            new_variable(ncFileID, cv_name, cv_type,    &
                                         (/get_dim_id(new_coordinate_var)/), &
                                         var_range     = cv_range,           &
                                         var_long_name = cv_long_name,       &
                                         var_units     = cv_units))
        else
            call set_var_id(new_coordinate_var,                          &
                            new_variable(ncFileID, cv_name, cv_type,   &
                                         (/get_dim_id(new_coordinate_var)/), &
                                         var_long_name = cv_long_name,       &
                                         var_units     = cv_units))
        end if

        if(present(cv_stagger)) call set_coord_stagger(new_coordinate_var, cv_stagger)
        if(present(cv_dir)    ) call set_coord_dir    (new_coordinate_var, cv_dir    )
        if(present(cv_nest)   ) call set_coord_nest   (new_coordinate_var, cv_nest   )

    end function new_coordinate_var

    ! define_latlon_vars
    ! -------------------------
    ! Defines that lat/lon variable on a specified mesh for a given nest
    function define_latlon_vars(ncFileID, staggering, var_dims, nest_number, coords) result(varids)
      integer,                            intent(in)  :: ncFileID
      character(len=*),                   intent(in)  :: staggering     
      integer,              dimension(:), intent(in)  :: var_dims
      integer,                            intent(in)  :: nest_number
      type(coordinate_var), dimension(:), intent(in)  :: coords

      integer, dimension(2)                           :: varids
      integer, dimension(2)                           :: dimids

      character(len=NF90_MAX_NAME)                    :: coord_str

      dimids = get_dimids(2, var_dims, nest_number, coords)

401 format(A,'_g',I2.2)

      write(coord_str, 401) 'lat'//trim(staggering), nest_number
      varids(LAT_DIM) = new_variable(ncFileID, coord_str, nf90_double, dimids,         &
                  trim(staggering), 'degrees_north','Latitude', LAT_LIMITS)

      write(coord_str, 401) 'lon'//trim(staggering), nest_number
      varids(LON_DIM) = new_variable(ncFileID, coord_str, nf90_double, dimids,           &
                  trim(staggering), 'degrees_east', 'Longitude',LON_LIMITS)

    end function define_latlon_vars

    ! fill_grid_latlon
    ! -------------------
    ! Fills the latitude and longitude netcdf fields on the staggered and unstaggered grids
    subroutine fill_nc_latlon(ncFileID, domain, nest_number, latlon_varid, latlon_coord)
      integer,                            intent(in) :: ncFileID      ! netCDF file 
      type(coamps_domain),                intent(in) :: domain
      integer,                            intent(in) :: nest_number 
      integer,              dimension(2), intent(in) :: latlon_varid
      type(coordinate_var), dimension(2), intent(in) :: latlon_coord

      integer            :: nlat, nlon
      integer            :: ii, jj
      real(kind=r8)      :: lat_stagger
      real(kind=r8)      :: lon_stagger
      type(coamps_nest)  :: nest

      real(kind=r8), allocatable, dimension(:,:)  :: lat
      real(kind=r8), allocatable, dimension(:,:)  :: lon

      ! Error handling
      integer                     :: alloc_status, dealloc_status
      character(len=*), parameter :: routine = 'fill_nc_latlon'

      nest = get_domain_nest(domain, nest_number)

      lat_stagger  = get_coord_stagger(latlon_coord(LAT_DIM))
      lon_stagger  = get_coord_stagger(latlon_coord(LON_DIM))

      nlat = get_dim_length(latlon_coord(LAT_DIM))
      nlon = get_dim_length(latlon_coord(LON_DIM))

      allocate(lon(nlon, nlat), stat=alloc_status)
      call check_alloc_status(alloc_status, routine, source, revision,   &
           revdate, 'lon(nlon, nlat)')

      allocate(lat(nlon, nlat), stat=alloc_status)
      call check_alloc_status(alloc_status, routine, source, revision,   &
           revdate, 'lat(nlon, nlat)')

      do ii = 1, nlon 
        do jj = 1, nlat 
          call nest_point_to_latlon(                     &
               domain  = domain,                         &
               nest_pt = make_nest_point(nest,           &
                  (real(ii,kind=r8) + lon_stagger),  &
                  (real(jj,kind=r8) + lat_stagger)), &
               lat = lat(ii,jj), lon = lon(ii,jj) )
        end do
      end do

      call nc_check(nf90_put_var(ncFileID, latlon_varid(LON_DIM), lon), routine)
      call nc_check(nf90_put_var(ncFileID, latlon_varid(LAT_DIM), lat), routine)

      deallocate(lon, stat=dealloc_status)
      call check_dealloc_status(alloc_status, routine, source, revision, &
                                revdate, 'lon')

      deallocate(lat, stat=dealloc_status)
      call check_dealloc_status(alloc_status, routine, source, revision, &
                                revdate, 'lat')

    end subroutine fill_nc_latlon

    ! populate_integer_variable
    ! -------------------------
    ! Writes a given *scalar* integer value to a variable
    subroutine populate_integer_variable(ncFileID, var_id, value, location)
        integer, intent(in)                     :: ncFileID
        integer, intent(in)                     :: var_id
        integer, dimension(:), intent(in)       :: value
        integer,       optional,     intent(in) :: location

        character(len=*), parameter :: routine = 'populate_integer_variable'

        call nc_check(nf90_put_var(ncFileID, var_id, value), routine) 
    end subroutine populate_integer_variable

    ! populate_real_variable
    ! ----------------------
    ! Writes a given *scalar* integer value to a variable
    subroutine populate_real_variable(ncFileID, var_id, values, location)
        integer,       intent(in) :: ncFileID
        integer,       intent(in) :: var_id
        real(kind=r8), dimension(:), intent(in) :: values
        integer,       optional,     intent(in) :: location

        character(len=*), parameter :: routine = 'populate_real_variable'

        call nc_check(nf90_put_var(ncFileID, var_id, values), routine) 
    end subroutine populate_real_variable

    ! write_state_coordinate_data
    ! ---------------------
    ! Write coordinate variable data to the file
    subroutine write_state_coordinate_data(ncFileID, state_coord, lon_coord, &
                                     lat_coord, vert_coord, locations)
        integer,                            intent(in) :: ncFileID
        type(coordinate_var),               intent(in) :: state_coord
        type(coordinate_var),               intent(in) :: lat_coord
        type(coordinate_var),               intent(in) :: lon_coord
        type(coordinate_var),               intent(in) :: vert_coord
        type(location_type),  dimension(:), intent(in) :: locations

        type(location_type)         :: loc        
        real(kind=r8), dimension(3) :: loc3d      ! lon/lat/hgt
        real(kind=r8), allocatable  :: lon_location(:)
        real(kind=r8), allocatable  :: lat_location(:)
        real(kind=r8), allocatable  :: vert_location(:)

        integer                     :: ii         

        allocate( lon_location(size(locations)))
        allocate( lat_location(size(locations)))
        allocate(vert_location(size(locations)))

        do ii=1, size(locations)
          loc3d             = get_location(locations(ii))
          lon_location(ii)  = loc3d(1)
          lat_location(ii)  = loc3d(2)
          vert_location(ii) = loc3d(3)
        end do

        call populate_variable(ncFileID, get_var_id(state_coord), (/(ii,ii=1,size(locations))/) )
        call populate_variable(ncFileID, get_var_id(lon_coord),   lon_location)
        call populate_variable(ncFileID, get_var_id(lat_coord),   lat_location)
        call populate_variable(ncFileID, get_var_id(vert_coord), vert_location)

        deallocate(lon_location)
        deallocate(lat_location)
        deallocate(vert_location)

    end subroutine write_state_coordinate_data

    ! set_dim_id
    ! -------------
    ! Sets the dimension ID of a coordinate variable
    subroutine set_dim_id(coord_var, new_dim_id)
        type(coordinate_var), intent(inout) :: coord_var
        integer,              intent(in)    :: new_dim_id

        coord_var%dim_id = new_dim_id
    end subroutine set_dim_id

    ! set_coord_nest
    ! -------------
    ! Sets the nest number of the coordinate variable
    subroutine set_coord_nest(coord_var, coordinate_nest)
        type(coordinate_var), intent(inout) :: coord_var
        integer,              intent(in)    :: coordinate_nest

        coord_var%nest_number = coordinate_nest
    end subroutine set_coord_nest

    ! set_coord_dir
    ! -------------
    ! Sets the dirction of the coordinate variable
    subroutine set_coord_dir(coord_var, coordinate_direction)
        type(coordinate_var), intent(inout) :: coord_var
        integer,              intent(in)    :: coordinate_direction

        coord_var%coordinate_direction = coordinate_direction
    end subroutine set_coord_dir

    ! set_dim_length
    ! -------------
    ! Sets the dimension length of a coordinate variable
    subroutine set_dim_length(coord_var, new_dim_length)
        type(coordinate_var), intent(inout) :: coord_var
        integer,              intent(in)    :: new_dim_length

        coord_var%dim_length = new_dim_length
    end subroutine set_dim_length

    ! set_var_id
    ! -------------
    ! Sets the variable ID of a coordinate variable
    subroutine set_var_id(coord_var, new_var_id)
        type(coordinate_var), intent(inout) :: coord_var
        integer,              intent(in)    :: new_var_id

        coord_var%var_id = new_var_id
    end subroutine set_var_id

    ! set_coord_stagger
    ! -------------
    ! Sets the variable staggering of a coordinate variable
    subroutine set_coord_stagger(coord_var, new_var_stagger)
        type(coordinate_var), intent(inout) :: coord_var
        real(kind=r8),        intent(in)    :: new_var_stagger

        coord_var%stagger = new_var_stagger
    end subroutine set_coord_stagger

    ! get_dim_id
    ! -------------
    ! Returns the dimension ID of a coordinate variable
    function get_dim_id(coord_var)
        type(coordinate_var), intent(in)  :: coord_var
        integer                           :: get_dim_id

        get_dim_id = coord_var%dim_id
    end function get_dim_id

    ! get_dim_length
    ! -------------
    ! Returns the dimension length of a coordinate variable
    function get_dim_length(coord_var)
        type(coordinate_var), intent(in)  :: coord_var
        integer                           :: get_dim_length

        get_dim_length = coord_var%dim_length
    end function get_dim_length

    ! get_coord_nest
    ! -------------
    ! Gets the nest number of the coordinate variable
    function get_coord_nest(coord_var)
        type(coordinate_var), intent(in) :: coord_var
        integer                          :: get_coord_nest
        get_coord_nest=coord_var%nest_number
    end function get_coord_nest

    ! get_coord_dir
    ! -------------
    ! Gets the dirction of the coordinate variable
    function get_coord_dir(coord_var)
        type(coordinate_var), intent(in) :: coord_var
        integer                          :: get_coord_dir
        get_coord_dir=coord_var%coordinate_direction
    end function get_coord_dir

    ! get_var_id
    ! -------------
    ! Returns the variable ID of a coordinate variable
    function get_var_id(coord_var)
        type(coordinate_var), intent(in)  :: coord_var
        integer                           :: get_var_id

        get_var_id = coord_var%var_id
    end function get_var_id

    ! get_coord_stagger
    ! -------------
    ! Returns the variable staggering of a coordinate variable
    function get_coord_stagger(coord_var)
        type(coordinate_var), intent(in)  :: coord_var
        real(kind=r8)                     :: get_coord_stagger

        get_coord_stagger = coord_var%stagger
    end function get_coord_stagger

    ! get_copy_dim_id
    ! ---------------
    ! Get the "copy" dimension from a NetCDF file
    function get_copy_dim_id(ncFileID)
        integer, intent(in) :: ncFileID
        integer             :: get_copy_dim_id

        call nc_check(nf90_inq_dimid(ncFileID, 'copy', dimid=get_copy_dim_id), 'get_copy_dim_id') 
    end function get_copy_dim_id

    ! get_time_dim_id
    ! ---------------
    ! Get the "time" dimension from a NetCDF file
    function get_time_dim_id(ncFileID)
        integer, intent(in) :: ncFileID
        integer             :: get_time_dim_id

        call nc_check(nf90_inq_dimid(ncFileID, 'time', dimid=get_time_dim_id), 'get_time_dim_id') 
    end function get_time_dim_id

  ! get_dimids
  ! -------------------
  ! Given a list of variable dimensions, try to figure out
  ! the matching netCDF dimension IDs (the table of netCDF
  ! dimension IDs and sizes is stored in a local variable
  ! to avoid incessant netCDF queries).
  !  PARAMETERS
  !   IN  ndims           the rank of the variable
  !   IN  var_dims        the array of the dimension of the variable
  !   IN  localdimds      the table of netCDF dim IDs and sizes
  !   OUT mydims          the array of netCDF dimension IDs 
  function get_dimids( ndims, var_dims, nest_number, local_coords ) result(mydims)
    integer,                            intent(in) :: ndims
    integer,              dimension(:), intent(in) :: var_dims
    integer,                            intent(in) :: nest_number
    type(coordinate_var), dimension(:), intent(in) :: local_coords
    integer,              dimension(SIZE(var_dims))  :: mydims
    integer,              dimension(3) :: COORD_DIRS = (/XDIR_COORD, YDIR_COORD, ZDIR_COORD/)
    integer :: i, j

    mydims = -999  ! set array to a bad value

    do i=1,ndims

    dimloop: do j=1,size(local_coords)
        if(var_dims(i)   == get_dim_length(local_coords(j)) .and. &
           COORD_DIRS(i) == get_coord_dir(local_coords(j))  .and. &
           (nest_number  == get_coord_nest(local_coords(j)) .or.  COORD_DIRS(i) == ZDIR_COORD) &                           &
          ) then
             mydims(i) = get_dim_id(local_coords(j))
             exit dimloop
        endif
    enddo dimloop

    enddo

    if ( any(mydims(1:ndims) < 0 )) then
       write(msgstring,*) 'ERROR : cannot find dimids to match ',mydims(1:ndims)
       call error_handler(E_ERR, 'get_dimids', msgstring, source, revision, revdate) 
    endif
  end function get_dimids

  ! get_var_dims
  ! -------------------
  ! Gets flag indicating if var is mean field
  !  PARAMETERS
  !   IN     var               state variable
  function get_var_dims(var, domain)
    type(state_variable), intent(in) :: var
    type(coamps_domain),  intent(in) :: domain
    integer, dimension(3)            :: get_var_dims

    integer           :: nx, ny, nz

    nx = get_nest_i_width(get_domain_nest(domain, get_nest_number(var)))
    ny = get_nest_j_width(get_domain_nest(domain, get_nest_number(var)))
    nz = get_var_nlevs(var, domain)
     
    select case (get_var_stagger(var))
      case ('U')
        get_var_dims = (/nx-1, ny, nz/)
      case ('V')
        get_var_dims = (/nx, ny-1, nz/)
      case default
        get_var_dims = (/nx, ny, nz/)
    end select
  end function get_var_dims

  ! get_var_rank
  ! -------------------
  ! Returns the number of non-singleton dimensions
  !  PARAMETERS
  !   OUT ndims             number of non-singleton dimensions
  !   IN  iterator          iterator for variable list
  function get_var_rank(var, domain) result(ndims)
    type(state_variable), intent(in) :: var
    type(coamps_domain),  intent(in) :: domain
    integer, dimension(3)            :: var_dims
    integer                          :: ndims, i

    var_dims = get_var_dims(var, domain)
    ndims=3
    dimloop: do i=1,3 
      if(var_dims(i) == 1) ndims=ndims-1
    enddo dimloop
  end function get_var_rank
    
    !------------------------------
    ! END PRIVATE ROUTINES
    !------------------------------

end module coamps_netcdf_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
