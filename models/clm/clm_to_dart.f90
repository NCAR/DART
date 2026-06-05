! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program clm_to_dart

!-------------------------------------------------------------------------------
! purpose: Update CLM output files to be suitable as input to filter.
!          Current functionality:
!          1. Replaces random values in empty snow layers with a consistent _FillValue
!             in the CLM restart file
!          2. When estimate_params = .true. in model_nml, expands the compact CLM
!             parameter file (clm_param_file) from scalar/1D-PFT form onto the full
!             2D lat/lon grid, writing clm_param_expanded_file.  This spatial
!             expansion (augmentation) is required so DART can assign a 
!             unique lat/lon location to each parameter element and apply spatial
!             localization correctly.
!             After assimilation, dart_to_clm averages the posterior 2D fields back
!             to a single representative value per PFT and writes it back to the
!             original parameter file.
!
! USAGE:  The clm filename is read from the clm_to_dart_nml namelist in input.nml.
!         Parameter expansion settings are read from model_nml in input.nml.
!         clm_to_dart
!
! Tim Hoar added empty snow layer functionality (2 June 2021)
! Brett Raczka  added parameter augmentation functionality (5 June 2026)
!-------------------------------------------------------------------------------

use        types_mod, only : r8, obstypelength

use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             error_handler, E_MSG, E_ERR, to_upper

use netcdf_utilities_mod, only : nc_check, &
                                 nc_open_file_readwrite, &
                                 nc_open_file_readonly,  &
                                 nc_close_file, &
                                 nc_synchronize_file, &
                                 nc_begin_define_mode, &
                                 nc_end_define_mode, &
                                 nc_variable_exists, &
                                 nc_get_variable, &
                                 nc_put_variable, &
                                 nc_get_variable_size, &
                                 nc_get_variable_dimension_names, &
                                 nc_get_attribute_from_variable, &
                                 nc_get_dimension_size

use netcdf

implicit none

character(len=*), parameter :: source = 'clm_to_dart.f90'

!-------------------------------------------------------------------------------
! namelist variables
!-------------------------------------------------------------------------------

character(len=256) :: clm_restart_file = 'clm_restart.nc'

! Parameter file options (only used when estimate_params = .true. in model_nml).
! clm_param_file         : input compact parameter file (scalar or 1D PFT arrays).
!                          Set to (empty) to skip parameter expansion entirely.
! clm_param_expanded_file: output file with parameter values mapped onto the full
!                          2D lat/lon state space so DART can assign locations and
!                          apply spatial localization. This is a global parameter
!                          value mapped in 2D state space -- all grid points initially
!                          hold the same value, which then evolves through the DA.
character(len=256) :: clm_param_file          = ''
character(len=256) :: clm_param_expanded_file = 'clm_params_expanded.nc'

integer            :: verbose = 0

namelist /clm_to_dart_nml/ clm_restart_file,          &
                            clm_param_file,             &
                            clm_param_expanded_file,    &
                            verbose

!-------------------------------------------------------------------------------
! Variables read from model_nml to support parameter expansion.
! Only the variables declared here are updated; other model_nml entries
! present in input.nml are safely ignored by Fortran namelist I/O.
!-------------------------------------------------------------------------------

integer, parameter :: max_state_variables    = 40
integer, parameter :: num_state_table_columns = 6
integer, parameter :: max_param_pfts          = 8

logical            :: estimate_params     = .false.
character(len=256) :: clm_history_filename = 'clm_history.nc'
integer            :: assimilate_pfts(max_param_pfts) = -1
character(len=obstypelength) :: &
   clm_variables(max_state_variables * num_state_table_columns) = ' '

! Partial read of model_nml -- only the four variables above are populated.
namelist /model_nml/ clm_history_filename, estimate_params, &
                     assimilate_pfts, clm_variables

!-------------------------------------------------------------------------------
! global storage
!-------------------------------------------------------------------------------

integer :: iunit, ncid, ncolumn, xtype
integer :: io, nvariables, ivar, ndims
integer :: i, j, nlevsno, numsnowlevels
integer :: dimids(NF90_MAX_VAR_DIMS)
integer :: dimlen(NF90_MAX_VAR_DIMS)
character(len=NF90_MAX_NAME) :: dimnames(NF90_MAX_VAR_DIMS)
character(len=NF90_MAX_NAME) :: varname

real(r8), allocatable :: SNOW_DEPTH(:) ! "snow depth"
real(r8), allocatable :: H2OSNO(:)     ! "snow water" (in column - includes traces)
real(r8), allocatable :: frac_sno(:)   ! "fraction of ground covered by snow (0 to 1)"
integer,  allocatable :: SNLSNO(:)     ! "negative number of snow layers"

real(r8), allocatable :: variable(:,:)
real(r8)              :: FillValue
real(r8)              :: missingValue

character(len=512) :: string1, string2

!===============================================================================

call initialize_utilities(progname='clm_to_dart')

! Read the clm_to_dart namelist to get filenames and options.

call find_namelist_in_file("input.nml", "clm_to_dart_nml", iunit)
read(iunit, nml = clm_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "clm_to_dart_nml") ! closes, too.

ncid = nc_open_file_readwrite(clm_restart_file, 'open for unused snow layer value replacement')

call get_snow_metadata()

nlevsno = nc_get_dimension_size(ncid,'levsno') ! The number of snow layers

! Get the number of variables in the netCDF file.
! We will query each of them to see if they use the 'levsno' or 'levtot'
! dimension. If so, we have to replace the bogus values - being careful not
! to clobber any 'trace' snow amounts.

io = nf90_inquire(ncid, nvariables=nvariables)
call nc_check(io, source, 'determining number of variables in file')

if (verbose > 0) write(*,*)'There are ',nvariables,' variables in the file.'

VARIABLES : do ivar = 1,nvariables

   write(string1,*)'inquire variable number ',ivar
   io = nf90_inquire_variable(ncid, ivar, varname, xtype, ndims, dimids)
   call nc_check(io, source, string1)

   ! immediately skip the variables that cannot contain snow layers 
   if (xtype /= NF90_DOUBLE) cycle VARIABLES
   if (ndims /= 2)           cycle VARIABLES

   ! Now that we are guaranteed to have a 2D variable - skip anything
   ! that is not dimensioned (*levels*,column)
   call nc_get_variable_dimension_names(ncid, varname, dimnames)
   if (trim(dimnames(2)) /= 'column') cycle VARIABLES

   if ( verbose > 1 ) &
      write(*,*)trim(string1),' varname ',trim(varname), &
             ' dimensions are ', trim(dimnames(1)), ' ', trim(dimnames(2))

   ! For 2D variables, the levels are always the first dimension

   select case (dimnames(1))
      case ( 'levtot', 'levsno')
         if (verbose > 0) then
            write(string1,*)'variable # ',ivar,' is "',trim(varname), &
                            '" - replacing indeterminate values'
            call error_handler(E_MSG,'clm_to_dart',string1)
         endif

         call nc_get_variable_size(ncid, varname, dimlen)
         allocate( variable( dimlen(1),dimlen(2) ) )
         call nc_get_variable(ncid, varname, variable)
         call nc_get_attribute_from_variable(ncid,varname,'_FillValue',FillValue)
         call nc_get_attribute_from_variable(ncid,varname,'missing_value',missingValue)

         ! Replace the bogus values for layers we KNOW to be unused.
         ! The SNLSNO has the negative number of snow layers, so ...

         do j = 1, dimlen(2)  ! loop over columns
            numsnowlevels = abs(SNLSNO(j))

            if (verbose > 2) then
               ! Debug block to check what happens in related variables
               ! when there is no snow in the layer closest to the ground.
               if (numsnowlevels == 0 .and. &
                   variable(nlevsno,j) /= FillValue .and. &
                   variable(nlevsno,j) > 0.0_r8) then
                  write(*,*)trim(varname), nlevsno, j, &
                      variable(nlevsno,j), H2OSNO(j), frac_sno(j), SNOW_DEPTH(j)
               endif
            endif

            
                ! Prevent unused snow layers from being updated by the assimilation.
                ! Unused snow layers have indeterminate values. The indeterminate values are replaced with 
                ! FillValue which prevents filter from updating the unused snow layers during the assimilation.
                do i = 1, nlevsno - numsnowlevels  ! loop over unused layers

                    
                     ! trace amounts of snow are in the level closest to the ground
                     ! frac_sno(j) seems to be a reliable indicator of a trace of snow
                     if (frac_sno(j) > 0.0_r8 .and. i == nlevsno) cycle
                       
                     variable(i,j) = FillValue
            
                enddo

         enddo

         ! update the netCDF file - in place

         call nc_put_variable(ncid,varname,variable,'replacing bogus values')
         deallocate(variable)

      case default
         continue
   end select

enddo VARIABLES

call nc_close_file(ncid)

deallocate(SNLSNO, H2OSNO, frac_sno, SNOW_DEPTH)

!-------------------------------------------------------------------------------
! Block 2: Parameter file spatial expansion.
! Only runs if clm_param_file is non-empty in clm_to_dart_nml.
!-------------------------------------------------------------------------------

if (len_trim(clm_param_file) > 0) then

   ! Read the relevant portions of model_nml.
   call find_namelist_in_file("input.nml", "model_nml", iunit)
   read(iunit, nml = model_nml, iostat = io)
   call check_namelist_read(iunit, io, "model_nml")

   if (estimate_params) then
      write(string1,*)'Expanding parameter file "'//trim(clm_param_file)//'"'
      write(string2,*)'onto 2D lat/lon grid -> "'//trim(clm_param_expanded_file)//'"'
      call error_handler(E_MSG, 'clm_to_dart', string1, text2=string2)
      call expand_params_to_grid()
   else
      write(string1,*)'clm_param_file is set but estimate_params = .false. in model_nml.'
      write(string2,*)'Parameter expansion skipped. Set estimate_params = .true. to enable.'
      call error_handler(E_MSG, 'clm_to_dart', string1, text2=string2)
   endif

endif

call finalize_utilities('clm_to_dart')

!===============================================================================
contains
!===============================================================================


!-------------------------------------------------------------------------------
! Expands a compact CLM parameter file (scalar or 1D PFT arrays) onto the full
! 2D lat/lon grid, writing the result to clm_param_expanded_file.
!
! For PFT-indexed variables (dimension 'pft'): one 2D (lat,lon) field is created
! per active PFT index in assimilate_pfts, named varname_pftNN (e.g. jmaxb0_pft07).
! For scalar variables (0 dimensions): a single 2D field is created retaining
! the original variable name.
!
! The 2D grid is read from clm_history_filename.
! Metadata (units, long_name, _FillValue) is copied from the source variable.
! Global attributes record the source file, grid source, and expanded PFT indices.

subroutine expand_params_to_grid()

character(len=*), parameter :: routine = 'expand_params_to_grid'

! netCDF IDs
integer :: ncid_src, ncid_hist, ncid_out
integer :: dimid_lon, dimid_lat
integer :: varid_lon_out, varid_lat_out
integer :: varid_src, varid_out
integer :: io_nc

! Grid dimensions and coordinate arrays
integer :: nlon, nlat, npft_src
real(r8), allocatable :: lon_vals(:), lat_vals(:)
real(r8), allocatable :: expanded_2d(:,:)   ! (nlon, nlat) -- Fortran order

! Source variable inquiry
integer :: ndims_src, src_dimids(NF90_MAX_VAR_DIMS), src_xtype
character(len=NF90_MAX_NAME) :: src_dim1_name

! Parsed param variable names from clm_variables namelist
character(len=NF90_MAX_NAME) :: param_varnames(max_state_variables)
integer :: num_param_vars

! Active PFT indices (0-based) extracted from assimilate_pfts
integer :: active_pfts(max_param_pfts), num_active_pfts

! Output variable registry -- built during define phase, used during write phase.
! Stores the output variable name, source variable name, and PFT index for each
! output field. PFT index = -1 indicates a scalar variable.
integer, parameter :: max_out_vars = max_state_variables * max_param_pfts
character(len=NF90_MAX_NAME) :: out_varname(max_out_vars)
character(len=NF90_MAX_NAME) :: out_srcname(max_out_vars)
integer                      :: out_pft(max_out_vars)
integer                      :: nout_vars

! Loop indices and temp variables
integer :: ivar, p, k, ovar
integer :: pft0, pft1                       ! 0-based and 1-based PFT index
real(r8), allocatable :: pft_array(:)       ! 1D source PFT array
real(r8)              :: scalar_val, fill_val
character(len=obstypelength) :: origin_str
character(len=NF90_MAX_NAME) :: varname_out
character(len=512)    :: attr_units, attr_longname
character(len=512)    :: pft_index_list

! -----------------------------------------------------------------------
! Step 1: Parse clm_variables to collect names of 'param'-origin variables
! -----------------------------------------------------------------------

num_param_vars = 0

do ivar = 1, max_state_variables
   ! Each variable occupies num_state_table_columns consecutive slots.
   ! col1=varname (offset -5), col5=origin (offset -1)
   if (trim(clm_variables(num_state_table_columns*ivar - 5)) == ' ') exit
   origin_str = trim(clm_variables(num_state_table_columns*ivar - 1))
   call to_upper(origin_str)
   if (trim(origin_str) == 'PARAM') then
      num_param_vars = num_param_vars + 1
      param_varnames(num_param_vars) = &
         trim(clm_variables(num_state_table_columns*ivar - 5))
   endif
enddo

if (num_param_vars == 0) then
   write(string1,*)'estimate_params = .true. but no ''param'' entries in clm_variables.'
   write(string2,*)'No parameter variables to expand. Check model_nml:clm_variables.'
   call error_handler(E_MSG, routine, string1, text2=string2)
   return
endif

write(string1,'(a,i0,a)')'Found ',num_param_vars,' param variable(s) to expand to 2D grid.'
call error_handler(E_MSG, routine, string1)

! -----------------------------------------------------------------------
! Step 2: Collect active PFT indices from assimilate_pfts (0-based)
! -----------------------------------------------------------------------

num_active_pfts = 0
do k = 1, max_param_pfts
   if (assimilate_pfts(k) >= 0) then
      num_active_pfts = num_active_pfts + 1
      active_pfts(num_active_pfts) = assimilate_pfts(k)
   endif
enddo

if (verbose > 0 .and. num_active_pfts > 0) then
   write(string1,'(a,i0,a)')'Active PFT indices (0-based): ', num_active_pfts,' PFT(s):'
   call error_handler(E_MSG, routine, string1)
   do k = 1, num_active_pfts
      write(string1,'(a,i0)')'  PFT index: ', active_pfts(k)
      call error_handler(E_MSG, routine, string1)
   enddo
endif

! -----------------------------------------------------------------------
! Step 3: Read grid dimensions and coordinate values from history file
! -----------------------------------------------------------------------

ncid_hist = nc_open_file_readonly(clm_history_filename, routine)

nlon = nc_get_dimension_size(ncid_hist, 'lon')
nlat = nc_get_dimension_size(ncid_hist, 'lat')

allocate(lon_vals(nlon), lat_vals(nlat), expanded_2d(nlon, nlat))

call nc_get_variable(ncid_hist, 'lon', lon_vals, routine)
call nc_get_variable(ncid_hist, 'lat', lat_vals, routine)

call nc_close_file(ncid_hist, routine)

write(string1,'(a,i0,a,i0)')'Grid from history file: nlon=',nlon,' nlat=',nlat
call error_handler(E_MSG, routine, string1)

! -----------------------------------------------------------------------
! Step 4: Open source param file and create output expanded file.
!         DEFINE PHASE: define all dimensions, coordinate vars, and
!         parameter output variables in a single pass before nf90_enddef.
! -----------------------------------------------------------------------

ncid_src = nc_open_file_readonly(clm_param_file, routine)

io_nc = nf90_create(trim(clm_param_expanded_file), NF90_CLOBBER, ncid_out)
call nc_check(io_nc, routine, 'creating '//trim(clm_param_expanded_file))

! Define dimensions (lon, lat) matching CLM history file convention
io_nc = nf90_def_dim(ncid_out, 'lon', nlon, dimid_lon)
call nc_check(io_nc, routine, 'defining lon dimension')
io_nc = nf90_def_dim(ncid_out, 'lat', nlat, dimid_lat)
call nc_check(io_nc, routine, 'defining lat dimension')

! Define coordinate variables
io_nc = nf90_def_var(ncid_out, 'lon', NF90_DOUBLE, [dimid_lon], varid_lon_out)
call nc_check(io_nc, routine, 'defining lon coordinate variable')
io_nc = nf90_put_att(ncid_out, varid_lon_out, 'long_name',    'coordinate longitude')
call nc_check(io_nc, routine, 'adding lon long_name')
io_nc = nf90_put_att(ncid_out, varid_lon_out, 'units',         'degrees_east')
call nc_check(io_nc, routine, 'adding lon units')
io_nc = nf90_put_att(ncid_out, varid_lon_out, 'cartesian_axis','X')
call nc_check(io_nc, routine, 'adding lon cartesian_axis')

io_nc = nf90_def_var(ncid_out, 'lat', NF90_DOUBLE, [dimid_lat], varid_lat_out)
call nc_check(io_nc, routine, 'defining lat coordinate variable')
io_nc = nf90_put_att(ncid_out, varid_lat_out, 'long_name',    'coordinate latitude')
call nc_check(io_nc, routine, 'adding lat long_name')
io_nc = nf90_put_att(ncid_out, varid_lat_out, 'units',         'degrees_north')
call nc_check(io_nc, routine, 'adding lat units')
io_nc = nf90_put_att(ncid_out, varid_lat_out, 'cartesian_axis','Y')
call nc_check(io_nc, routine, 'adding lat cartesian_axis')

! Define 2D parameter variables.
! Variables are stored as (lon, lat) in Fortran memory = (lat, lon) in netCDF convention,
! matching the CLM history file layout so that the existing fill_rank2_metadata CASE("lat")
! in model_mod.f90 correctly assigns lon/lat locations to each state vector element.

nout_vars = 0

do ivar = 1, num_param_vars

   varname_out = trim(param_varnames(ivar))

   if (.not. nc_variable_exists(ncid_src, varname_out)) then
      write(string1,*)'Parameter variable "'//trim(varname_out)//'" not found in "'// &
                       trim(clm_param_file)//'"'
      call error_handler(E_ERR, routine, string1)
   endif

   io_nc = nf90_inq_varid(ncid_src, trim(varname_out), varid_src)
   call nc_check(io_nc, routine, 'inquiring varid for '//trim(varname_out))
   io_nc = nf90_inquire_variable(ncid_src, varid_src, xtype=src_xtype, &
                                  ndims=ndims_src, dimids=src_dimids)
   call nc_check(io_nc, routine, 'inquiring variable '//trim(varname_out))

   ! Get fill value (non-fatal if absent; use a large negative default)
   fill_val = -9999.0_r8
   io_nc = nf90_get_att(ncid_src, varid_src, '_FillValue', fill_val)

   ! Get optional metadata attributes
   attr_units    = ''
   attr_longname = ''
   io_nc = nf90_get_att(ncid_src, varid_src, 'units',     attr_units)
   io_nc = nf90_get_att(ncid_src, varid_src, 'long_name', attr_longname)

   if (ndims_src == 1) then

      ! PFT-indexed variable: check that the dimension is 'pft'
      io_nc = nf90_inquire_dimension(ncid_src, src_dimids(1), name=src_dim1_name)
      call nc_check(io_nc, routine, 'inquiring dimension for '//trim(varname_out))

      if (trim(src_dim1_name) /= 'pft') then
         write(string1,*)'Variable "'//trim(varname_out)//'" has 1D dimension "'// &
                          trim(src_dim1_name)//'" (expected ''pft'').'
         write(string2,*)'Only ''pft'' and scalar (0D) dimensions are supported for param variables.'
         call error_handler(E_ERR, routine, string1, text2=string2)
      endif

      if (num_active_pfts == 0) then
         write(string1,*)'PFT-indexed variable "'//trim(varname_out)//'" found but'
         write(string2,*)'no valid PFT indices in assimilate_pfts. Skipping this variable.'
         call error_handler(E_MSG, routine, string1, text2=string2)
         cycle
      endif

      ! Define one 2D output variable per active PFT
      do p = 1, num_active_pfts
         pft0 = active_pfts(p)

         ! Output variable name: varname_pftNN (0-based, zero-padded to 2 digits)
         write(varname_out, '(a,a,i2.2)') trim(param_varnames(ivar)), '_pft', pft0

         io_nc = nf90_def_var(ncid_out, trim(varname_out), NF90_DOUBLE, &
                              [dimid_lon, dimid_lat], varid_out)
         call nc_check(io_nc, routine, 'defining variable '//trim(varname_out))

         if (len_trim(attr_longname) > 0) then
            io_nc = nf90_put_att(ncid_out, varid_out, 'long_name', trim(attr_longname))
            call nc_check(io_nc, routine, 'adding long_name for '//trim(varname_out))
         endif
         if (len_trim(attr_units) > 0) then
            io_nc = nf90_put_att(ncid_out, varid_out, 'units', trim(attr_units))
            call nc_check(io_nc, routine, 'adding units for '//trim(varname_out))
         endif
         io_nc = nf90_put_att(ncid_out, varid_out, '_FillValue', fill_val)
         call nc_check(io_nc, routine, 'adding _FillValue for '//trim(varname_out))

         write(string2, '(i0)') pft0
         io_nc = nf90_put_att(ncid_out, varid_out, 'source_pft_index_0based', pft0)
         call nc_check(io_nc, routine, 'adding source_pft_index for '//trim(varname_out))
         io_nc = nf90_put_att(ncid_out, varid_out, 'source_variable', trim(param_varnames(ivar)))
         call nc_check(io_nc, routine, 'adding source_variable for '//trim(varname_out))
         io_nc = nf90_put_att(ncid_out, varid_out, 'expansion_type', 'PFT-indexed')
         call nc_check(io_nc, routine, 'adding expansion_type for '//trim(varname_out))

         ! Register in output variable table for write phase
         nout_vars = nout_vars + 1
         out_varname(nout_vars) = trim(varname_out)
         out_srcname(nout_vars) = trim(param_varnames(ivar))
         out_pft(    nout_vars) = pft0

         if (verbose > 0) then
            write(string1,*)'Defined: ',trim(param_varnames(ivar)), &
               ' PFT ',pft0,' -> ',trim(varname_out)
            call error_handler(E_MSG, routine, string1)
         endif

      enddo  ! active PFTs

   else if (ndims_src == 0) then

      ! Scalar variable: replicate across the full 2D grid under original name
      io_nc = nf90_def_var(ncid_out, trim(param_varnames(ivar)), NF90_DOUBLE, &
                           [dimid_lon, dimid_lat], varid_out)
      call nc_check(io_nc, routine, 'defining scalar variable '//trim(param_varnames(ivar)))

      if (len_trim(attr_longname) > 0) then
         io_nc = nf90_put_att(ncid_out, varid_out, 'long_name', trim(attr_longname))
         call nc_check(io_nc, routine, 'adding long_name for '//trim(param_varnames(ivar)))
      endif
      if (len_trim(attr_units) > 0) then
         io_nc = nf90_put_att(ncid_out, varid_out, 'units', trim(attr_units))
         call nc_check(io_nc, routine, 'adding units for '//trim(param_varnames(ivar)))
      endif
      io_nc = nf90_put_att(ncid_out, varid_out, '_FillValue', fill_val)
      call nc_check(io_nc, routine, 'adding _FillValue for '//trim(param_varnames(ivar)))
      io_nc = nf90_put_att(ncid_out, varid_out, 'source_variable', trim(param_varnames(ivar)))
      call nc_check(io_nc, routine, 'adding source_variable for '//trim(param_varnames(ivar)))
      io_nc = nf90_put_att(ncid_out, varid_out, 'expansion_type', 'scalar')
      call nc_check(io_nc, routine, 'adding expansion_type for '//trim(param_varnames(ivar)))

      ! Register in output variable table for write phase (-1 flags as scalar)
      nout_vars = nout_vars + 1
      out_varname(nout_vars) = trim(param_varnames(ivar))
      out_srcname(nout_vars) = trim(param_varnames(ivar))
      out_pft(    nout_vars) = -1

      if (verbose > 0) then
         write(string1,*)'Defined: scalar ',trim(param_varnames(ivar)),' -> 2D field'
         call error_handler(E_MSG, routine, string1)
      endif

   else
      write(string1,*)'Parameter variable "'//trim(param_varnames(ivar))//'" has rank ', &
                       ndims_src,'. Only scalar (0D) and 1D PFT-indexed variables are supported.'
      call error_handler(E_ERR, routine, string1)
   endif

enddo  ! param variables (define phase)

! Add global attributes before ending define mode
io_nc = nf90_put_att(ncid_out, NF90_GLOBAL, 'source_param_file', trim(clm_param_file))
call nc_check(io_nc, routine, 'adding global attr source_param_file')
io_nc = nf90_put_att(ncid_out, NF90_GLOBAL, 'grid_source', trim(clm_history_filename))
call nc_check(io_nc, routine, 'adding global attr grid_source')
io_nc = nf90_put_att(ncid_out, NF90_GLOBAL, 'description', &
   'Global CLM parameter values mapped onto 2D lat/lon state space for DART assimilation. '// &
   'Each grid point initially holds the same value (replicated from the source PFT slot). '// &
   'After assimilation, dart_to_clm averages the posterior 2D field back to a single '// &
   'representative value per PFT and writes it back to the original parameter file.')
call nc_check(io_nc, routine, 'adding global attr description')

! Build a string listing the active PFT indices for the global attribute
pft_index_list = ''
do k = 1, num_active_pfts
   write(string1, '(i0)') active_pfts(k)
   if (k == 1) then
      pft_index_list = trim(string1)
   else
      pft_index_list = trim(pft_index_list)//' '//trim(string1)
   endif
enddo
io_nc = nf90_put_att(ncid_out, NF90_GLOBAL, 'assimilated_pft_indices_0based', trim(pft_index_list))
call nc_check(io_nc, routine, 'adding global attr assimilated_pft_indices_0based')

! End define mode -- all dimensions and variables are now defined
io_nc = nf90_enddef(ncid_out)
call nc_check(io_nc, routine, 'ending define mode for '//trim(clm_param_expanded_file))

! -----------------------------------------------------------------------
! Step 5: Write coordinate data
! -----------------------------------------------------------------------

io_nc = nf90_inq_varid(ncid_out, 'lon', varid_lon_out)
call nc_check(io_nc, routine, 'inquiring lon varid for writing')
io_nc = nf90_put_var(ncid_out, varid_lon_out, lon_vals)
call nc_check(io_nc, routine, 'writing lon coordinate values')

io_nc = nf90_inq_varid(ncid_out, 'lat', varid_lat_out)
call nc_check(io_nc, routine, 'inquiring lat varid for writing')
io_nc = nf90_put_var(ncid_out, varid_lat_out, lat_vals)
call nc_check(io_nc, routine, 'writing lat coordinate values')

! -----------------------------------------------------------------------
! Step 6: Write 2D parameter fields (data phase, after nf90_enddef)
! -----------------------------------------------------------------------

do ovar = 1, nout_vars

   io_nc = nf90_inq_varid(ncid_out, trim(out_varname(ovar)), varid_out)
   call nc_check(io_nc, routine, 'inquiring varid for writing '//trim(out_varname(ovar)))

   io_nc = nf90_inq_varid(ncid_src, trim(out_srcname(ovar)), varid_src)
   call nc_check(io_nc, routine, 'inquiring source varid for '//trim(out_srcname(ovar)))

   if (out_pft(ovar) >= 0) then

      ! PFT-indexed: read full PFT array, extract the one value at Fortran index pft0+1
      npft_src = nc_get_dimension_size(ncid_src, 'pft')
      allocate(pft_array(npft_src))
      call nc_get_variable(ncid_src, trim(out_srcname(ovar)), pft_array, routine)

      pft0 = out_pft(ovar)
      pft1 = pft0 + 1   ! convert 0-based to Fortran 1-based

      if (pft1 > npft_src) then
         write(string1,*)'PFT index ',pft0,' (0-based, Fortran index ',pft1, &
            ') exceeds pft dimension size ',npft_src,' in "'//trim(clm_param_file)//'"'
         call error_handler(E_ERR, routine, string1)
      endif

      scalar_val = pft_array(pft1)   ! save before dealloc for diagnostic logging
      expanded_2d(:,:) = scalar_val
      deallocate(pft_array)

      if (verbose > 0) then
         write(string1,'(a,a,a,i0,a,es12.5)') &
            'Writing ',trim(out_varname(ovar)),': PFT ',pft0,' value = ',scalar_val
         call error_handler(E_MSG, routine, string1)
      endif

   else

      ! Scalar: read the single value directly
      io_nc = nf90_get_var(ncid_src, varid_src, scalar_val)
      call nc_check(io_nc, routine, 'reading scalar variable '//trim(out_srcname(ovar)))
      expanded_2d(:,:) = scalar_val

      if (verbose > 0) then
         write(string1,'(a,a,a,es12.5)') &
            'Writing ',trim(out_varname(ovar)),': scalar value = ',scalar_val
         call error_handler(E_MSG, routine, string1)
      endif

   endif

   io_nc = nf90_put_var(ncid_out, varid_out, expanded_2d)
   call nc_check(io_nc, routine, 'writing '//trim(out_varname(ovar)))

enddo  ! output variables (write phase)

! -----------------------------------------------------------------------
! Step 7: Finalize
! -----------------------------------------------------------------------

call nc_close_file(ncid_src, routine)

io_nc = nf90_close(ncid_out)
call nc_check(io_nc, routine, 'closing '//trim(clm_param_expanded_file))

deallocate(lon_vals, lat_vals, expanded_2d)

write(string1,*)'Parameter expansion complete: ',nout_vars,' 2D field(s) written to "'// &
                 trim(clm_param_expanded_file)//'"'
call error_handler(E_MSG, routine, string1)

end subroutine expand_params_to_grid


!-------------------------------------------------------------------------------
!> The SNLSNO variable contains the information about how many snow layers are
!> being used in each column.
!> H2OSNO is used to detect the columns that have trace amounts of snow.

subroutine get_snow_metadata()

! If SNLSNO does not exist, there is nothing in 
! the file that has bogus values in the snow layers.

if ( .not. nc_variable_exists(ncid,'SNLSNO') ) then
   write(string1,*)'"'//trim(clm_restart_file)//'" has no "SNLSNO" variable.'
   string2 = 'CLM restart files have "SNLSNO", clm_to_dart only works with restart files.'
   call error_handler(E_ERR,'get_snow_metadata',string1,source,text2=string2)
endif

! All variables are dimensioned the same size
call nc_get_variable_size(ncid, 'SNLSNO', ncolumn)
allocate(SNLSNO(ncolumn),H2OSNO(ncolumn),frac_sno(ncolumn),SNOW_DEPTH(ncolumn))

call nc_get_variable(ncid,'SNLSNO',SNLSNO)
call nc_get_variable(ncid,'frac_sno',frac_sno)
call nc_get_variable(ncid,'SNOW_DEPTH',SNOW_DEPTH)

! ctsm5.1.dev043 replaced the H2OSNO variable with H2OSNO_NO_LAYERS
if ( nc_variable_exists(ncid,'H2OSNO')) then
   call nc_get_variable(ncid,'H2OSNO',H2OSNO)

else if (nc_variable_exists(ncid,'H2OSNO_NO_LAYERS')) then
   call nc_get_variable(ncid,'H2OSNO_NO_LAYERS',H2OSNO)

else
   string1 = 'neither "H2OSNO" nor "H2OSNO_NO_LAYERS" exists in file'
   string2 = 'one or the other is needed'
   call error_handler(E_ERR,'get_snow_metadata',string1,source,text2=string2)
endif

if (verbose > 1) write(*,*)'minval of SNLSNO is ',minval(SNLSNO)

!SNLSNO has a _FillValue of -9999 but there should never be any of these 

call nc_get_attribute_from_variable(ncid,'SNLSNO','_FillValue',FillValue)
if (any(SNLSNO == FillValue)) then
   write(string1,*)'SNLSNO has at least one _FillValue ... unexpected, unable to proceed.'
   call error_handler(E_ERR,'get_snow_metadata',string1,source)
endif

end subroutine get_snow_metadata


end program clm_to_dart
