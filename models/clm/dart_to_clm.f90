! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program dart_to_clm

!----------------------------------------------------------------------
! purpose: Update the CLM restart file with the DART posterior. If need
!          be, take the posterior SWE (a diagnostic variable) and
!          repartition it into prognostic snow layers.
!          When estimate_params = .true. in model_nml, also averages the
!          posterior 2D parameter fields (from clm_params_expanded.nc)
!          back to a single representative value per PFT and writes them
!          back to the compact parameter file (clm_params.nc).
!
! method: Read DART posterior and replace the valid values in the
!         CLM restart file. Anything with a DART posterior _FillValue
!         is replaced with the original CLM value.
!         For parameters: read the posterior 2D spatial field, compute
!         an area-weighted (or PFT-weighted) spatial mean per ensemble
!         member, and write back to the original compact parameter file.
!         This averaging is done independently for each ensemble member,
!         preserving ensemble spread in parameter space.
!
! discussion: Some variables in the CLM restart file have variables that
!         have neither the declared _FillValue nor a predictable value.
!         The 'clm_to_dart' program replaced those CLM indeterminate
!         values with the _FillValue so the DART netCDF readers correctly
!         replace those values with the DART MISSING value. The DART
!         netCDF write routines replace the DART MISSING values with the
!         variables declared _FillValue. This routine replaces the
!         _FillValue with whatever was originally in the CLM restart file.
!
!         This routine also gives the user the option to manually repartition
!         the prognostic snow-related variable updates from the DART update
!         of the diagnostic total column SWE variable H2OSNO. The default
!         behavior is for the prognostic snow related layers to be updated
!         directly from DART.
!----------------------------------------------------------------------

use        types_mod, only : r8, MISSING_R8, obstypelength

use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             open_file, close_file, &
                             logfileunit, nmlfileunit, &
                             do_nml_file, do_nml_term, &
                             error_handler, E_MSG, E_ERR, to_upper

use time_manager_mod, only : time_type, print_time, print_date, operator(/=)

use        model_mod, only : static_init_model, read_model_time

use   state_structure_mod, only : get_num_variables, &
                                  get_num_dims, get_dim_name,         &
                                  get_dim_length, get_variable_name,  &
                                  do_io_update, get_variable_size,    &
                                  get_model_variable_indices,         &
                                  get_has_missing_value, &
                                  get_missing_value, &
                                  get_FillValue, &
                                  get_domain_size, &
                                  get_varid_from_varname

use netcdf_utilities_mod, only : nc_check,                        &
                                 nc_open_file_readonly,            &
                                 nc_open_file_readwrite,           &
                                 nc_get_variable_num_dimensions,   &
                                 nc_get_variable_dimension_names,  &
                                 nc_get_variable_size,             &
                                 nc_get_variable,                  &
                                 nc_put_variable,                  &
                                 nc_close_file,                    &
                                 nc_variable_exists,               &
                                 nc_get_dimension_size,            &
                                 NF90_MAX_NAME, NF90_MAX_VAR_DIMS

use netcdf

implicit none

character(len=*), parameter :: source = 'dart_to_clm.f90'

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character(len=256) :: dart_to_clm_input_file  = 'dart_posterior.nc'
character(len=256) :: dart_to_clm_output_file = 'clm_restart.nc'
integer            :: repartition_swe          = 0
character(len=256) :: repartition_vhist_file   = 'clm_vector_history.nc'
character(len=256) :: repartition_analysis_file = 'dart_posterior_vector.nc'
integer            :: verbose                  = 0

! Parameter estimation options (only used when estimate_params = .true. in model_nml).
!
! dart_to_clm_param_expanded_file: DART posterior for the parameter domain.
!   This 2D lat/lon file was created by clm_to_dart and updated by filter.
!   Each grid point holds the posterior parameter value at that location.
!
! dart_to_clm_param_output_file: compact parameter file to update.
!   The spatially-averaged posterior values are written back here
!   (at the correct PFT slot) for use in the next CLM forecast.
!
! param_spatial_avg_method: collapses the 2D posterior field to one value.
!   'AREA_WEIGHTED' (default): weighted mean using gridcell area; grid cells
!      where the target PFT is absent (weight = 0) are excluded automatically.
!      Correct geometric weighting; recommended for general use.
!   'PFT_WEIGHTED': weighted mean using fractional spatial coverage of the
!      specific target PFT at each grid cell. Grid cells with higher PFT
!      coverage contribute more. Requires pfts1d_itype_veg in restart file;
!      falls back to AREA_WEIGHTED if unavailable.
!   Note: each ensemble member is averaged independently -- the spatially-
!   averaged parameter value for each member is unique, preserving spread.
character(len=256) :: dart_to_clm_param_expanded_file = 'clm_params_expanded.nc'
character(len=256) :: dart_to_clm_param_output_file   = 'clm_params.nc'
character(len=32)  :: param_spatial_avg_method         = 'AREA_WEIGHTED'

namelist /dart_to_clm_nml/ dart_to_clm_input_file,          &
                            dart_to_clm_output_file,         &
                            repartition_swe,                 &
                            repartition_vhist_file,          &
                            repartition_analysis_file,       &
                            dart_to_clm_param_expanded_file, &
                            dart_to_clm_param_output_file,   &
                            param_spatial_avg_method,        &
                            verbose

!------------------------------------------------------------------
! Variables for model_nml read.
! ALL variables present in input.nml:model_nml must be declared here;
! Fortran namelist I/O sets IOSTAT non-zero for any undeclared entry,
! which DART's check_namelist_read treats as a fatal error.
! Variables not needed by this program are declared with defaults and unused.
!------------------------------------------------------------------

integer, parameter :: max_state_variables_d2c    = 40
integer, parameter :: num_state_table_columns_d2c = 6
integer, parameter :: max_param_pfts_d2c          = 8

! Variables actively used by dart_to_clm
logical            :: estimate_params      = .false.
character(len=256) :: clm_history_filename = 'clm_history.nc'
character(len=256) :: clm_restart_filename = 'clm_restart.nc'
integer            :: assimilate_pfts(max_param_pfts_d2c) = -1
character(len=obstypelength) :: &
   clm_variables(max_state_variables_d2c * num_state_table_columns_d2c) = ' '

! Variables declared to satisfy namelist read; not used by this program
character(len=256) :: clm_vector_history_filename = 'clm_vector_history.nc'
character(len=256) :: clm_parameter_filename      = 'clm_params_expanded.nc'
integer            :: assimilation_period_days     = 0
integer            :: assimilation_period_seconds  = 21600
integer            :: debug                        = 0
character(len=32)  :: calendar                     = 'Gregorian'

namelist /model_nml/ clm_restart_filename, clm_history_filename,         &
                     clm_vector_history_filename, clm_parameter_filename, &
                     assimilation_period_days, assimilation_period_seconds, &
                     calendar, debug,                                       &
                     estimate_params, assimilate_pfts, clm_variables

!----------------------------------------------------------------------

integer            :: iunit, io, dom_restart, ivar, rank, irank
integer            :: nlevsno, ICE_varsize(2)
integer            :: ncid_dart, ncid_clm
type(time_type)    :: dart_time, clm_time

character(len=512) :: string1, string2, string3
character(len=NF90_MAX_NAME) :: varname
!======================================================================

call initialize_utilities(progname='dart_to_clm')

! Call model_mod:static_init_model() to reads the clm namelist

call static_init_model()

! Read the namelist to get the input/output filenames. 

call find_namelist_in_file("input.nml", "dart_to_clm_nml", iunit)
read(iunit, nml = dart_to_clm_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_clm_nml")

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=dart_to_clm_nml)
if (do_nml_term()) write(     *     , nml=dart_to_clm_nml)

ncid_dart = nc_open_file_readonly( dart_to_clm_input_file,  'opening to read posterior')
ncid_clm  = nc_open_file_readwrite(dart_to_clm_output_file, 'opening to write restart')

write(string1,*)'reading updated values from DART file "'//trim(dart_to_clm_input_file)//'"'
write(string2,*)'and writing to CLM restart file "'//trim(dart_to_clm_output_file)//'"'
call error_handler(E_MSG,source,string1,text2=string2)

! Make sure we are updating the right restart file - 

dart_time = read_model_time(dart_to_clm_input_file)
clm_time  = read_model_time(dart_to_clm_output_file)

if (dart_time /= clm_time) then
   write(string1,*)'DART timestamp "'//trim(dart_to_clm_input_file)//'" /= '
   write(string2,*)'CLM  timestamp "'//trim(dart_to_clm_output_file)//'"'
   call error_handler(E_MSG,source,string1,text2=string2)
endif

call print_date( clm_time,'dart_to_clm:clm  model date')
call print_time( clm_time,'dart_to_clm:DART model time')
call print_date( clm_time,'dart_to_clm:clm  model date',logfileunit)
call print_time( clm_time,'dart_to_clm:DART model time',logfileunit)

! update the current clm restart file
!
! Get the DART posterior variable names from the restart domain (domain 1).
! The other domains do not need to be updated - they are not used for
! forecasting.

! There are some dependencies on the state structure, which come from the shapefile,
! and then the two files in question must also be compatible ...


dom_restart = 1

UPDATE : do ivar=1, get_num_variables(dom_restart)
  
  ! If repartitioning SWE is enabled, then skip applying snow variable updates
  ! within UPDATE loop which applies un-repartitioned DART posterior
  ! The repartitioned snow variable updates will be performed in 'update_snow'
  if (repartition_swe > 0) then


   varname = get_variable_name(dom_restart,ivar)
   select case (varname)
      case ('SNOWDP', 'SNOW_DEPTH', 'DZSNO', 'ZSNO', 'ZISNO', 'H2OSOI_LIQ', 'H2OSOI_ICE')
         write(string1,*)'re-partitioning of SWE is enabled for '//trim(varname)
         write(string2,*)'posterior is coming from repartitioning of H2OSNO.'
         call error_handler(E_MSG, 'dart_to_clm', string1, text2=string2)
         cycle UPDATE
      case default
      ! do nothing -- proceed to default DART update subroutines:
      ! replace_values_1D and replace_values_2D  
   end select
   endif
  
   rank = get_num_dims(dom_restart,ivar)

   if (rank == 1) then
      call replace_values_1D(dom_restart, ivar, ncid_dart, ncid_clm)

   elseif (rank == 2) then
      call replace_values_2D(dom_restart, ivar, ncid_dart, ncid_clm)

   else
      write(string1, *) 'no support for data array of dimension ', rank
      call error_handler(E_ERR, source, string1, source)
   endif

enddo UPDATE

 ! Manually repartition snow layer variables with H2OSNO variable

  if (repartition_swe > 0) then
 
 ! Pass in the soil/snow and column dimensions to update_snow subroutine
     if (nc_variable_exists(ncid_clm, 'H2OSOI_ICE')) then
        call nc_get_variable_size(ncid_clm, 'H2OSOI_ICE', ICE_varsize, source)
        nlevsno = nc_get_dimension_size(ncid_clm, 'levsno')
     else
        write(string1,*)'Snow repartitioning requires "H2OSOI_ICE" variable'
        call error_handler(E_ERR,source,string1,source)
     endif
     
     call update_snow(dom_restart, ncid_dart, ncid_clm, ICE_varsize(2), &
                      ICE_varsize(1), nlevsno, repartition_swe)

  endif



call nc_close_file(ncid_clm,  source)
call nc_close_file(ncid_dart, source)

!------------------------------------------------------------------
! Parameter domain update (domain 4).
! Read model_nml to determine whether parameter estimation is active.
! If so, average the posterior 2D parameter fields back to a single
! representative value per PFT and write to the compact param file.
! This is done per ensemble member (dart_to_clm is called once per
! member by assimilate.csh), so each member's result is independent.
!------------------------------------------------------------------

call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

if (estimate_params) then
   write(string1,*)'Parameter estimation active: averaging posterior 2D fields'
   write(string2,*)'from "'//trim(dart_to_clm_param_expanded_file)// &
                   '" -> "'//trim(dart_to_clm_param_output_file)//'"'
   call error_handler(E_MSG, source, string1, text2=string2)
   call average_params_to_scalar()
endif

! Log what we think we're doing, and exit.

call finalize_utilities('dart_to_clm')


!===============================================================================
contains
!===============================================================================

!------------------------------------------------------------------
!  Reads the posterior 2D parameter fields from dart_to_clm_param_expanded_file,
!  computes a spatial mean for each field (per ensemble member, independently),
!  and writes the result back to the correct PFT slot in dart_to_clm_param_output_file.
!
! Two averaging methods are supported (param_spatial_avg_method namelist option):
!   AREA_WEIGHTED: weighted mean using gridcell area from the CLM history file.
!      Grid cells where the target PFT is absent (pfts1d_wtxy = 0) are excluded.
!   PFT_WEIGHTED:  weighted mean using the fractional spatial coverage of the
!      target PFT at each grid cell. Requires pfts1d_itype_veg in the CLM restart
!      file; falls back to AREA_WEIGHTED if that variable is not available.
!
!  For PFT-indexed parameters: the source_pft_index_0based attribute on each
!    expanded variable identifies which PFT slot in the compact file to update.
!  For scalar parameters: the mean is written back as a single scalar value.
! 
!  Diagnostics: the prior value (currently in the compact file) and the posterior
!   mean are both logged for each parameter variable.

subroutine average_params_to_scalar()

character(len=*), parameter :: routine = 'average_params_to_scalar'

! netCDF file handles
integer :: ncid_expanded, ncid_compact, ncid_hist, ncid_rst
integer :: varid_expanded, varid_compact, varid_area
integer :: io_nc

! Grid dimensions
integer :: nlon, nlat, npft_compact, npft_src

! Param variable info parsed from clm_variables namelist
character(len=NF90_MAX_NAME) :: param_varnames(max_state_variables_d2c)
integer :: num_param_vars

! Active PFTs (0-based)
integer :: active_pfts(max_param_pfts_d2c), num_active_pfts
integer :: pft0, pft1

! Spatial data arrays
real(r8), allocatable :: posterior_2d(:,:)   ! (nlon, nlat) posterior field
real(r8), allocatable :: area_2d(:,:)        ! (nlon, nlat) gridcell area (km^2)
real(r8), allocatable :: pft_cover_2d(:,:)   ! (nlon, nlat) PFT fractional coverage
real(r8), allocatable :: pft_array_compact(:) ! 1D PFT array from compact file
real(r8), allocatable :: pfts1d_wtxy(:)      ! PFT weights from restart file (sparse)
integer,  allocatable :: pfts1d_ixy(:)       ! PFT lon indices from restart file
integer,  allocatable :: pfts1d_jxy(:)       ! PFT lat indices from restart file
integer,  allocatable :: pfts1d_itype_veg(:) ! PFT vegetation type (0-based) if available

! Averaging variables
real(r8) :: weighted_sum, weight_sum, posterior_mean, prior_val
real(r8) :: fill_val
integer  :: i, j, k, ivar, npft_rst
character(len=obstypelength) :: origin_str
character(len=NF90_MAX_NAME) :: expanded_varname
character(len=32)  :: avg_method
logical :: pft_type_available

! Variable type detection -- mirrors expand_params_to_grid in clm_to_dart.f90.
! Each variable is checked against the compact param file to determine whether
! it is PFT-indexed (ndims=1, dim='pft') or scalar (ndims=0), independently of
! whether num_active_pfts > 0.  This correctly handles mixed lists of PFT and
! scalar params in clm_variables.
integer :: ndims_src, src_dimids(NF90_MAX_VAR_DIMS)
integer :: varid_compact_chk
character(len=NF90_MAX_NAME) :: src_dim1_name
logical :: is_pft_indexed

! -----------------------------------------------------------------------
! Step 1: Parse clm_variables to collect 'param'-origin variable names
! -----------------------------------------------------------------------

num_param_vars = 0
do ivar = 1, max_state_variables_d2c
   if (trim(clm_variables(num_state_table_columns_d2c*ivar - 5)) == ' ') exit
   origin_str = trim(clm_variables(num_state_table_columns_d2c*ivar - 1))
   call to_upper(origin_str)
   if (trim(origin_str) == 'PARAM') then
      num_param_vars = num_param_vars + 1
      param_varnames(num_param_vars) = &
         trim(clm_variables(num_state_table_columns_d2c*ivar - 5))
   endif
enddo

if (num_param_vars == 0) then
   write(string1,*)'estimate_params = .true. but no ''param'' entries in clm_variables.'
   write(string2,*)'No parameter averaging performed.'
   call error_handler(E_MSG, routine, string1, text2=string2)
   return
endif

! -----------------------------------------------------------------------
! Step 2: Collect active PFT indices (0-based) from assimilate_pfts
! -----------------------------------------------------------------------

num_active_pfts = 0
do k = 1, max_param_pfts_d2c
   if (assimilate_pfts(k) >= 0) then
      num_active_pfts = num_active_pfts + 1
      active_pfts(num_active_pfts) = assimilate_pfts(k)
   endif
enddo

! -----------------------------------------------------------------------
! Step 3: Open files
! -----------------------------------------------------------------------

ncid_expanded = nc_open_file_readonly( dart_to_clm_param_expanded_file, routine)
ncid_compact  = nc_open_file_readwrite(dart_to_clm_param_output_file,    routine)
ncid_hist     = nc_open_file_readonly( clm_history_filename,             routine)
ncid_rst      = nc_open_file_readonly( clm_restart_filename,             routine)

! -----------------------------------------------------------------------
! Step 4: Read grid dimensions and gridcell area from history file
! -----------------------------------------------------------------------

nlon = nc_get_dimension_size(ncid_hist, 'lon')
nlat = nc_get_dimension_size(ncid_hist, 'lat')

allocate(posterior_2d(nlon, nlat), area_2d(nlon, nlat), pft_cover_2d(nlon, nlat))

call nc_get_variable(ncid_hist, 'area', area_2d, routine)

call nc_close_file(ncid_hist, routine)

! -----------------------------------------------------------------------
! Step 5: Read PFT spatial data from restart file for weighting.
!         pfts1d_wtxy: weight of each PFT instance relative to its gridcell.
!         pfts1d_ixy/jxy: lon/lat grid index of each PFT instance.
!         pfts1d_itype_veg: PFT vegetation type (0-based). Optional -- used
!           only for PFT_WEIGHTED mode to identify which instances belong to
!           the target PFT. Falls back to AREA_WEIGHTED if unavailable.
! -----------------------------------------------------------------------

npft_rst = nc_get_dimension_size(ncid_rst, 'pft')
allocate(pfts1d_wtxy(npft_rst), pfts1d_ixy(npft_rst), pfts1d_jxy(npft_rst))

call nc_get_variable(ncid_rst, 'pfts1d_wtxy', pfts1d_wtxy, routine)
call nc_get_variable(ncid_rst, 'pfts1d_ixy',  pfts1d_ixy,  routine)
call nc_get_variable(ncid_rst, 'pfts1d_jxy',  pfts1d_jxy,  routine)

! Check whether PFT vegetation type is available (needed for PFT_WEIGHTED)
pft_type_available = nc_variable_exists(ncid_rst, 'pfts1d_itype_veg')
if (pft_type_available) then
   allocate(pfts1d_itype_veg(npft_rst))
   call nc_get_variable(ncid_rst, 'pfts1d_itype_veg', pfts1d_itype_veg, routine)
else
   write(string1,*)'pfts1d_itype_veg not found in "'//trim(clm_restart_filename)//'".'
   write(string2,*)'PFT_WEIGHTED averaging requires this variable. Falling back to AREA_WEIGHTED.'
   call error_handler(E_MSG, routine, string1, text2=string2)
endif

call nc_close_file(ncid_rst, routine)

! Validate/normalise the averaging method string
avg_method = trim(param_spatial_avg_method)
call to_upper(avg_method)
if (trim(avg_method) /= 'AREA_WEIGHTED' .and. trim(avg_method) /= 'PFT_WEIGHTED') then
   write(string1,*)'Unrecognised param_spatial_avg_method: "'//trim(param_spatial_avg_method)//'"'
   write(string2,*)'Valid options: AREA_WEIGHTED, PFT_WEIGHTED. Defaulting to AREA_WEIGHTED.'
   call error_handler(E_MSG, routine, string1, text2=string2)
   avg_method = 'AREA_WEIGHTED'
endif
if (trim(avg_method) == 'PFT_WEIGHTED' .and. .not. pft_type_available) then
   avg_method = 'AREA_WEIGHTED'
endif

write(string1,*)'Spatial averaging method: '//trim(avg_method)
call error_handler(E_MSG, routine, string1)

! -----------------------------------------------------------------------
! Step 6: Loop over param variables and compute spatial averages
! -----------------------------------------------------------------------

do ivar = 1, num_param_vars

   ! Determine the actual variable type from the compact param file.
   ! This mirrors the logic in expand_params_to_grid (clm_to_dart.f90) and ensures
   ! scalar params are handled correctly even when num_active_pfts > 0.
   io_nc = nf90_inq_varid(ncid_compact, trim(param_varnames(ivar)), varid_compact_chk)
   call nc_check(io_nc, routine, 'inquiring varid for '//trim(param_varnames(ivar)))
   io_nc = nf90_inquire_variable(ncid_compact, varid_compact_chk, &
                                  ndims=ndims_src, dimids=src_dimids)
   call nc_check(io_nc, routine, 'inquiring variable '//trim(param_varnames(ivar)))

   is_pft_indexed = .false.
   if (ndims_src == 1) then
      io_nc = nf90_inquire_dimension(ncid_compact, src_dimids(1), name=src_dim1_name)
      call nc_check(io_nc, routine, 'inquiring dim for '//trim(param_varnames(ivar)))
      if (trim(src_dim1_name) == 'pft') is_pft_indexed = .true.
   endif

   if (is_pft_indexed) then

      ! --- PFT-indexed variables: one 2D field per active PFT ---
      if (num_active_pfts == 0) then
         write(string1,*)'PFT-indexed variable "'//trim(param_varnames(ivar))// &
                         '" found but no valid PFT indices in assimilate_pfts. Skipping.'
         call error_handler(E_MSG, routine, string1)
         cycle
      endif
      do k = 1, num_active_pfts
         pft0 = active_pfts(k)
         pft1 = pft0 + 1   ! Fortran 1-based index into compact file

         ! Construct the expanded field name (e.g. jmaxb0_pft07)
         write(expanded_varname, '(a,a,i2.2)') &
            trim(param_varnames(ivar)), '_pft', pft0

         if (.not. nc_variable_exists(ncid_expanded, expanded_varname)) then
            write(string1,*)'Expanded variable "'//trim(expanded_varname)//'" not found.'
            write(string2,*)'Skipping. Check that clm_to_dart ran successfully.'
            call error_handler(E_MSG, routine, string1, text2=string2)
            cycle
         endif

         ! Read 2D posterior field
         call nc_get_variable(ncid_expanded, trim(expanded_varname), posterior_2d, routine)

         ! Build 2D weight field based on selected method
         pft_cover_2d(:,:) = 0.0_r8
         if (trim(avg_method) == 'PFT_WEIGHTED') then
            ! Weight by PFT fractional coverage of target PFT type
            do j = 1, npft_rst
               if (pfts1d_itype_veg(j) == pft0) then
                  pft_cover_2d(pfts1d_ixy(j), pfts1d_jxy(j)) = &
                     pft_cover_2d(pfts1d_ixy(j), pfts1d_jxy(j)) + pfts1d_wtxy(j)
               endif
            enddo
         else
            ! AREA_WEIGHTED: use PFT presence mask (any instance of target PFT)
            do j = 1, npft_rst
               if (pfts1d_itype_veg(j) == pft0) then
                  pft_cover_2d(pfts1d_ixy(j), pfts1d_jxy(j)) = 1.0_r8  ! mark as present
               endif
            enddo
         endif

         ! For AREA_WEIGHTED: multiply presence mask by gridcell area
         ! For PFT_WEIGHTED:  multiply PFT coverage by gridcell area (for area correction)
         ! In both cases: weight(i,j) = pft_cover_2d(i,j) * area_2d(i,j)
         !   AREA_WEIGHTED: pft_cover_2d is 0/1 presence mask  -> area of cells with PFT
         !   PFT_WEIGHTED:  pft_cover_2d is fractional coverage -> area * frac

         weighted_sum = 0.0_r8
         weight_sum   = 0.0_r8
         do j = 1, nlat
            do i = 1, nlon
               if (pft_cover_2d(i,j) > 0.0_r8) then
                  weighted_sum = weighted_sum + posterior_2d(i,j) * pft_cover_2d(i,j) * area_2d(i,j)
                  weight_sum   = weight_sum   +                      pft_cover_2d(i,j) * area_2d(i,j)
               endif
            enddo
         enddo

         if (weight_sum > 0.0_r8) then
            posterior_mean = weighted_sum / weight_sum
         else
            ! No grid cells with target PFT found -- use unweighted mean as fallback
            write(string1,*)'No grid cells with PFT ',pft0,' found for weighting.'
            write(string2,*)'Using unweighted mean of all grid points as fallback.'
            call error_handler(E_MSG, routine, string1, text2=string2)
            posterior_mean = sum(posterior_2d) / real(nlon*nlat, r8)
         endif

         ! Read prior value from compact file for diagnostic logging
         npft_compact = nc_get_dimension_size(ncid_compact, 'pft')
         allocate(pft_array_compact(npft_compact))
         call nc_get_variable(ncid_compact, trim(param_varnames(ivar)), &
                              pft_array_compact, routine)
         prior_val = pft_array_compact(pft1)

         write(string1,'(a,a,a,i0,a,es12.5,a,es12.5)') &
            'Param update: ', trim(param_varnames(ivar)), &
            ' PFT ', pft0, &
            '  prior=', prior_val, '  posterior=', posterior_mean
         call error_handler(E_MSG, routine, string1)

         ! Write posterior mean back to compact file at the correct PFT slot
         pft_array_compact(pft1) = posterior_mean
         call nc_put_variable(ncid_compact, trim(param_varnames(ivar)), &
                              pft_array_compact, routine)

         deallocate(pft_array_compact)

      enddo  ! active PFTs

   else  ! scalar (ndims=0) or non-pft 1D -- treated as single global value

      ! --- Scalar variable: single 2D field, no PFT index ---
      if (.not. nc_variable_exists(ncid_expanded, trim(param_varnames(ivar)))) then
         write(string1,*)'Scalar expanded variable "'//trim(param_varnames(ivar))//'" not found.'
         write(string2,*)'Skipping. Check that clm_to_dart ran successfully.'
         call error_handler(E_MSG, routine, string1, text2=string2)
         cycle
      endif

      call nc_get_variable(ncid_expanded, trim(param_varnames(ivar)), posterior_2d, routine)

      ! Simple area-weighted mean for scalars (no PFT mask needed)
      weighted_sum = 0.0_r8
      weight_sum   = 0.0_r8
      do j = 1, nlat
         do i = 1, nlon
            weighted_sum = weighted_sum + posterior_2d(i,j) * area_2d(i,j)
            weight_sum   = weight_sum   +                      area_2d(i,j)
         enddo
      enddo
      posterior_mean = weighted_sum / weight_sum

      ! Read prior value (scalar)
      io_nc = nf90_inq_varid(ncid_compact, trim(param_varnames(ivar)), varid_compact)
      call nc_check(io_nc, routine, 'inquiring varid for '//trim(param_varnames(ivar)))
      io_nc = nf90_get_var(ncid_compact, varid_compact, prior_val)
      call nc_check(io_nc, routine, 'reading prior scalar '//trim(param_varnames(ivar)))

      write(string1,'(a,a,a,es12.5,a,es12.5)') &
         'Param update (scalar): ', trim(param_varnames(ivar)), &
         '  prior=', prior_val, '  posterior=', posterior_mean

      call error_handler(E_MSG, routine, string1)

      ! Write back as scalar
      io_nc = nf90_put_var(ncid_compact, varid_compact, posterior_mean)
      call nc_check(io_nc, routine, 'writing posterior scalar '//trim(param_varnames(ivar)))

   endif

enddo  ! param variables

! -----------------------------------------------------------------------
! Step 7: Clean up
! -----------------------------------------------------------------------

call nc_close_file(ncid_expanded, routine)
call nc_close_file(ncid_compact,  routine)

deallocate(posterior_2d, area_2d, pft_cover_2d)
deallocate(pfts1d_wtxy, pfts1d_ixy, pfts1d_jxy)
if (allocated(pfts1d_itype_veg)) deallocate(pfts1d_itype_veg)

write(string1,*)'Parameter averaging complete: ',num_param_vars,' variable(s) updated in "'// &
                 trim(dart_to_clm_param_output_file)//'"'
call error_handler(E_MSG, routine, string1)

end subroutine average_params_to_scalar


!------------------------------------------------------------------
!> The DART posterior has had all the DART MISSING_R8 values replaced by
!> the CLM declared '_FillValue' and the clamping values have been honored.
!> Get the 'original' variable from the netcdf file.

subroutine replace_values_1D(dom_id, ivar, ncid_dart, ncid_clm)

integer, intent(in) :: dom_id
integer, intent(in) :: ivar
integer, intent(in) :: ncid_dart
integer, intent(in) :: ncid_clm

character(len=*), parameter :: routine = 'replace_values_1D'

character(len=NF90_MAX_NAME) :: varname
real(r8), allocatable :: dart_array(:)
real(r8), allocatable :: clm_array(:)
real(r8) :: special
integer  :: varsize(1)

varname = get_variable_name(dom_id,ivar)

! I am not sure I like getting the special value from the state structure as
! opposed to getting it right from the value in the ncid_dart variable ...
! The state structure gets it from the 'shapefile' listed in model_nml ... 

! If it has both, the core DART routines have already made sure they are identical
! The problem here is that if only the missing_value -or- FillValue is used, the
! other may have a bogus value. We only have one logical query routine.  
if (get_has_missing_value(dom_id,ivar)) call get_missing_value(dom_id,ivar,special)
if (get_has_missing_value(dom_id,ivar)) call get_FillValue(    dom_id,ivar,special)

! Make sure variables in both files are identical in shape, etc.
call Compatible_Variables(varname, ncid_dart, ncid_clm, varsize)

allocate(dart_array(varsize(1)), clm_array(varsize(1)))

call nc_get_variable(ncid_clm,  varname,  clm_array)
call nc_get_variable(ncid_dart, varname, dart_array)

where(dart_array /= special) clm_array = dart_array

call nc_put_variable(ncid_clm, varname, clm_array, routine)

deallocate(dart_array, clm_array)

end subroutine replace_values_1D


!------------------------------------------------------------------
!>
!> The DART posterior has had all the DART MISSING_R8 values replaced by
!> the CLM declared '_FillValue' and the clamping values have been honored. 
!> Get the 'original' variable from the netcdf file.

subroutine replace_values_2D(dom_id, ivar, ncid_dart, ncid_clm)

integer, intent(in) :: dom_id
integer, intent(in) :: ivar
integer, intent(in) :: ncid_dart
integer, intent(in) :: ncid_clm

character(len=*), parameter :: routine = 'replace_values_2D'

character(len=NF90_MAX_NAME) :: varname
real(r8), allocatable :: dart_array(:,:)
real(r8), allocatable :: clm_array(:,:)
real(r8) :: special
integer  :: varsize(2)

varname = get_variable_name(dom_id,ivar)

! If it has both, the core DART routines have already made sure they are identical
! The problem here is that if only the missing_value -or- FillValue is used, the
! other may have a bogus value. We only have one logical query routine.  
if (get_has_missing_value(dom_id,ivar)) call get_missing_value(dom_id,ivar,special)
if (get_has_missing_value(dom_id,ivar)) call get_FillValue(    dom_id,ivar,special)

! Make sure variables in both files are identical in shape, etc.
call Compatible_Variables(varname, ncid_dart, ncid_clm, varsize)

allocate(dart_array(varsize(1),varsize(2)), clm_array(varsize(1),varsize(2)))

call nc_get_variable(ncid_clm,  varname,  clm_array)
call nc_get_variable(ncid_dart, varname, dart_array)

where(dart_array /= special) clm_array = dart_array

call nc_put_variable(ncid_clm, varname, clm_array, routine)

deallocate(dart_array, clm_array)

end subroutine replace_values_2D

!------------------------------------------------------------------

subroutine update_snow(dom_id, ncid_dart, ncid_clm, ncolumn, nlevel, nlevsno, repartition_swe)

!> This repartitions snow layer variables based upon DART update of diagnostic H2OSNO 
!> variable.  This subroutine expects H2OSNO variable is available, as well as 
!> prognostic clm variables: SNOWDP/SNOW_DEPTH, DZSNO, H2OSOI_LIQ, H2OSOI_ICE
!>                           ZSNO and ZISNO

! The posterior snow water equivalent (H2OSNO) cannot be zero or negative
! because H2OSNO is used to calculate the snow density -- which is used
! to calculate the snow layer depth.
! However the posterior H2OSNO values received from DART can be zero/negative
! thus  negative values of H2OSNO produced by DART, are set to prior value
! when calculating the snow density

! @fixme  Alternatively could automatically set all snow related variables to
!  zero, when H2OSNO posterior is zero/negative.

! Pseudo-code below:  Note this reparitioning is only applied to snow layers
! and subsurface layers are left unchanged
!
! SnowDensity(layer,column)  = (H2OSOI_LIQ(layer,column) + H2OSOI_ICE(layer,column))
!                              / DZSNO(layer,column)
! 
! wt_swe(layer,column)       = H2OSOI_LIQ(layer,column)+H2OSOI_ICE(layer,column))
!                              / H2OSNO(column)
! wt_liq(layer,column)       = H2OSOI_LIQ(layer,column)/(H2OSOI_LIQ(layer,column)+
!                               H2OSOI_ICE(layer,column))
! 
! wt_ice(layer,column)       = 1-wt_liq(layer,column)
! 
! GainH2OSNO_l(layer,column) = GainH2OSNO(column)*wt_swe(layer,column)
! 
! GainH2O_LIQ(layer,column)  = GainH2OSNO_l(layer,column)*wt_liq(layer,column)
! 
! GainH2O_ICE(layer,column)  = GainH2OSNO_l(layer,column)*(1-wt_liq(layer,column))
! 
! GainDZSNO(layer,column)    = GainH2OSNO_l(layer,column)/SnowDensity(layer,column)
! 
! H2OSOI_LIQ_update(layer,column) = H2OSOI_LIQ(layer,column)+GainH2O_LIQ(layer,column)
! 
! H2OSOI_ICE_update(layer,column) = H2OSOI_ICE(layer,column)+GainH2O_ICE(layer,column)
! 
! GainSNOWDP(column)         = sum(GainDZSNO(layer,column))
! 
! SNOWDP_update(column)      = SNOWDP(column)+GainSNOWDP(column)
!
!=========================================================================
!    double frac_sno(column) ;
!            frac_sno:long_name = "fraction of ground covered by snow (0 to 1)" ;
!            frac_sno:units = "unitless" ;
!    int    SNLSNO(column) ;
!            SNLSNO:long_name = "number of snow layers" ;
!            SNLSNO:units = "unitless" ;
!    double SNOWDP(column) ;
!            SNOWDP:long_name = "snow depth" ;
!            SNOWDP:units = "m" ;
!    double H2OSNO(column) ;
!            H2OSNO:long_name = "snow water" ;
!            H2OSNO:units = "kg/m2" ;
!    double H2OSOI_LIQ(column, levtot) ;
!            H2OSOI_LIQ:long_name = "liquid water" ;
!            H2OSOI_LIQ:units = "kg/m2" ;
!    double H2OSOI_ICE(column, levtot) ;
!            H2OSOI_ICE:long_name = "ice lens" ;
!            H2OSOI_ICE:units = "kg/m2" ;
!    double DZSNO(column, levsno) ;
!            DZSNO:long_name = "snow layer thickness" ;
!            DZSNO:units = "m" ;
!    double ZSNO(column, levsno) ;
!	     ZSNO:long_name = "snow layer depth" ;
!	     ZSNO:units = "m" ;
!    double ZISNO(column, levsno) ;
!            ZISNO:long_name = "snow interface depth at the top of the given layer" ;
!            ZISNO:units = "m" ;

integer, intent(in) :: dom_id
integer, intent(in) :: ncid_dart
integer, intent(in) :: ncid_clm
integer             :: ncid_clm_vector, ncid_dart_vector
character(len=*), parameter :: routine = 'update_snow'

! Total columns (ncolumn), soil/snow levels (nlevel) and snow levels (nlevsno)
! Need repartition_swe to distinguish between all-layer(1) and bottom layer(2)
! SWE partitioning
integer, intent(in) :: ncolumn, nlevel, nlevsno, repartition_swe
integer  :: icolumn, ilevel, c
integer  :: VarID, varsize(2)
real(r8) :: snowden, wt_swe, wt_liq, wt_ice

integer  :: ivar
character(len=NF90_MAX_NAME) :: varname
character(len=512) :: string1, string2
real(r8) :: special

real(r8) :: dart_H2OSNO(ncolumn)
real(r8) :: dart_SNOWDP(ncolumn)
real(r8) :: dart_DZSNO(nlevsno,ncolumn)
real(r8) :: dart_ZSNO(nlevsno,ncolumn)
real(r8) :: dart_ZISNO(nlevsno,ncolumn)
real(r8) :: dart_H2OLIQ(nlevel,ncolumn)
real(r8) :: dart_H2OICE(nlevel,ncolumn)

real(r8) :: clm_H2OSNO(ncolumn)  !(column,time) for vector history
integer  :: clm_SNLSNO(ncolumn)
real(r8) :: clm_SNOWDP(ncolumn)
real(r8) :: clm_DZSNO(nlevsno,ncolumn)
real(r8) :: clm_ZSNO(nlevsno,ncolumn)
real(r8) :: clm_ZISNO(nlevsno,ncolumn)
real(r8) :: clm_H2OLIQ(nlevel,ncolumn)
real(r8) :: clm_H2OICE(nlevel,ncolumn)

integer  :: snlsno(ncolumn)
real(r8) :: h2osno_pr(ncolumn)
real(r8) :: h2osno_po(ncolumn)
real(r8) :: snowdp_pr(ncolumn)
real(r8) :: snowdp_po(ncolumn)
real(r8) :: dzsno_pr(nlevsno,ncolumn)
real(r8) :: dzsno_po(nlevsno,ncolumn)
real(r8) :: zsno_pr(nlevsno,ncolumn)
real(r8) :: zsno_po(nlevsno,ncolumn)
real(r8) :: zisno_pr(nlevsno,ncolumn)
real(r8) :: zisno_po(nlevsno,ncolumn)
real(r8) :: h2oliq_pr(nlevel,ncolumn)
real(r8) :: h2oliq_po(nlevel,ncolumn)
real(r8) :: h2oice_pr(nlevel,ncolumn)
real(r8) :: h2oice_po(nlevel,ncolumn)
real(r8) :: gain_dzsno(nlevsno,ncolumn)
real(r8) :: gain_h2oice, gain_h2oliq, gain_h2osno

! Check existence for snow related variables required to be
! in DART state and CLM domain (restart,history,vector history)
! If they do not exist, throw error immediately, provide guidance.
! Also use this opportunity to locate proper SWE (H2OSNO) and snow
! depth clm variables. The H2OSNO variable (vectory history) is 
! required for repartitioning.  There are mutliple 'snow depth' 
! variable names depending on clm version.


if (verbose > 0) then

   ! restart domain    
   domain1 : do ivar=1, get_num_variables(1)
   varname = get_variable_name(1,ivar)
   write(string1,*)'CLM domain 1 variable:  "'//trim(varname)//'"' 
   call error_handler(E_MSG, routine, string1, source)
   enddo domain1
  
   ! history domain   
   domain2 : do ivar=1, get_num_variables(2)
   varname = get_variable_name(2,ivar)
   write(string1,*)'CLM domain 2 variable: "'//trim(varname)//'"' 
   call error_handler(E_MSG, routine, string1, source)
   enddo domain2

   ! history vector domain   
   domain3 : do ivar=1, get_num_variables(3)
   varname = get_variable_name(3,ivar)
   write(string1,*)'CLM domain 3 variable: "'//trim(varname)//'"'
   call error_handler(E_MSG, routine, string1, source)
   enddo domain3

endif


! BEGIN clm variable and DART state checks for snow repartitioning

if (nc_variable_exists(ncid_clm, 'SNLSNO') .and. nc_variable_exists(ncid_clm, 'DZSNO') .and. &
    nc_variable_exists(ncid_clm, 'H2OSOI_LIQ') .and. nc_variable_exists(ncid_clm, 'H2OSOI_ICE') &
    .and. nc_variable_exists(ncid_clm, 'ZSNO') .and. nc_variable_exists(ncid_clm, 'ZISNO')) then                             
    call nc_get_variable(ncid_clm, 'SNLSNO',  clm_SNLSNO)
    call nc_get_variable(ncid_clm, 'DZSNO',  clm_DZSNO)
    call nc_get_variable(ncid_clm, 'ZSNO',  clm_ZSNO)
    call nc_get_variable(ncid_clm, 'ZISNO',  clm_ZISNO)
    call nc_get_variable(ncid_clm, 'H2OSOI_LIQ',  clm_H2OLIQ)
    call nc_get_variable(ncid_clm, 'H2OSOI_ICE',  clm_H2OICE)
else
    write(string1,*)'Snow repartitioning requires clm variables: '
    write(string2,*)'SNLSNO, DZSNO, ZSNO, ZISNO, H2OSOI_LIQ, H2OSOI_ICE,'
    call error_handler(E_ERR,routine,string1,source,text2=string2)
endif

ncid_clm_vector= nc_open_file_readonly(repartition_vhist_file, &
                 'confirm H2OSNO is in clm_vector_history.nc')

if (nc_variable_exists(ncid_clm_vector, 'H2OSNO')) then
   call nc_get_variable(ncid_clm_vector, 'H2OSNO',  clm_H2OSNO)
else
   write(string1,*)'Snow repartitioning requires clm SWE variable'
   write(string2,*)'Check for H2OSNO within clm vector history'
   call error_handler(E_ERR,routine,string1,source,text2=string2)
endif

! Multiple options for clm snow depth variable
if (nc_variable_exists(ncid_clm, 'SNOW_DEPTH')) then
   call nc_get_variable(ncid_clm, 'SNOW_DEPTH',  clm_SNOWDP)
elseif (nc_variable_exists(ncid_clm, 'SNOWDP')) then
   call nc_get_variable(ncid_clm, 'SNOWDP',  clm_SNOWDP)
else
   write(string1,*)'Snow repartitioning requires clm snow depth variable'
   write(string2,*)'Check restart/history files for "SNOW_DEPTH" or "SNOWDP"'
   call error_handler(E_ERR,routine,string1,source,text2=string2)
endif


! Check for DART state posterior variables.
! H2OSNO posterior comes from dart_posterior_vector.nc (filter analysis stage)
! whereas all other posterior variables come from dart_posterior.nc

if (nc_variable_exists(ncid_dart, 'H2OSOI_LIQ') .and. nc_variable_exists(ncid_dart, 'DZSNO') &
    .and. nc_variable_exists(ncid_dart, 'H2OSOI_ICE') .and. nc_variable_exists(ncid_dart, 'ZSNO') &
    .and. nc_variable_exists(ncid_dart, 'ZISNO')) then
    call nc_get_variable(ncid_dart, 'DZSNO',  dart_DZSNO)
    call nc_get_variable(ncid_dart, 'ZSNO',  dart_ZSNO)
    call nc_get_variable(ncid_dart, 'ZISNO',  dart_ZISNO)
    call nc_get_variable(ncid_dart, 'H2OSOI_LIQ',  dart_H2OLIQ)
    call nc_get_variable(ncid_dart, 'H2OSOI_ICE',  dart_H2OICE)
else
    write(string1,*)'Snow repartitioning requires dart state variables: '
    write(string2,*)'DZSNO, ZSNO, ZISNO, H2OSOI_LIQ, H2OSOI_ICE'
    call error_handler(E_ERR,routine,string1,source,text2=string2)
endif

ncid_dart_vector= nc_open_file_readonly(repartition_analysis_file, &
                 'confirm H2OSNO is in dart_posterior_vector.nc')

if (nc_variable_exists(ncid_dart_vector, 'H2OSNO')) then
   call nc_get_variable(ncid_dart_vector, 'H2OSNO',  dart_H2OSNO, routine)
else
   write(string1,*)'Snow repartitioning requires H2OSNO in DART state'
   write(string2,*)'Check that H2OSNO is in DART model_nml as vector history'
   call error_handler(E_ERR,routine,string1,source,text2=string2)
endif

! Multiple options for snow depth variable
if (nc_variable_exists(ncid_dart, 'SNOW_DEPTH')) then
   call nc_get_variable(ncid_dart, 'SNOW_DEPTH',  dart_SNOWDP)
elseif (nc_variable_exists(ncid_dart, 'SNOWDP')) then
   call nc_get_variable(ncid_dart, 'SNOWDP',  dart_SNOWDP)
else
   write(string1,*)'Snow repartitioning requires clm snow depth variable'
   write(string2,*)'Check model_nml within input.nml for "SNOW_DEPTH" or "SNOWDP"'
   call error_handler(E_ERR,routine,string1,source,text2=string2)
endif

! END clm variable and DART state checks

h2osno_pr=clm_H2OSNO
snlsno=clm_SNLSNO
snowdp_pr=clm_SNOWDP
dzsno_pr=clm_DZSNO
zsno_pr=clm_ZSNO
zisno_pr=clm_ZISNO
h2oliq_pr=clm_H2OLIQ
h2oice_pr=clm_H2OICE
gain_dzsno(:,:)= 0.0_r8

! Likely unnecessary check coming from clm restart file
where (h2osno_pr < 0.0_r8) h2osno_pr = 0.0_r8
where (snowdp_pr < 0.0_r8) snowdp_pr = 0.0_r8
where ( dzsno_pr < 0.0_r8)  dzsno_pr = 0.0_r8
where (h2oliq_pr < 0.0_r8) h2oliq_pr = 0.0_r8
where (h2oice_pr < 0.0_r8) h2oice_pr = 0.0_r8

! The  H2OSNO posterior is taken from DART analysis stage (dart_posterior_vector.nc)
! The repartitioning of the updated H2OSNO is based on prior snow relative distribution

h2osno_po = dart_H2OSNO  
snowdp_po = snowdp_pr
dzsno_po =  dzsno_pr
zsno_po = zsno_pr
zisno_po = zisno_pr
h2oliq_po = h2oliq_pr
h2oice_po = h2oice_pr


! Adjust the variables to be consistent with the updated H2OSNO
! Remember: snlsno has the NEGATIVE of the number of snow layers.

! The posterior snow water equivalent (h2osno_po) cannot be zero or negative
! because H2OSNO is used to calculate the snow density -- which is used
! to calculate the snow layer depth.
! However the posterior H2OSNO values received from DART can be zero/negative
! *therefore*  negative values of H2OSNO produced by DART, are set to prior value
! when calculating the snow density

! Alternatively could automatically set all snow related variables to
! zero, when H2OSNO posterior is zero/negative. Currently not implmented.


PARTITION: do icolumn = 1,ncolumn
      
   if (h2osno_po(icolumn) < 0.0_r8) then  ! If there IS NOT snow in column...

      write(*,*)'negative value found: column ',icolumn
      
      ! Force any existing layers to be near zero but not exactly zero
      ! Leave layer aggregation/initialization to CLM
      do ilevel=1,-snlsno(icolumn)

         h2oliq_po(ilevel,icolumn) = 0.00000001_r8
         h2oice_po(ilevel,icolumn) = 0.00000001_r8
          dzsno_po(ilevel,icolumn) = 0.00000001_r8   

      enddo

      snowdp_po(icolumn) = sum( dzsno_po(:,icolumn))
      h2osno_po(icolumn) = sum(h2oice_po(:,icolumn))

   else
      
      if (snlsno(icolumn) < 0) then  ! If there IS snow in the column ...


         do c = 1,-snlsno(icolumn)   ! Identify total snow layers to repartition
            
            ! In CLM4 nlevsno=5, but in CLM5 nlevsno=12
            ilevel = nlevsno-c+1     ! Loop through each snow layer 

            ! Calculate the snow density for each layer
            if (dzsno_pr(ilevel,icolumn) > 0.0_r8) then
                snowden = (h2oliq_pr(ilevel,icolumn) + h2oice_pr(ilevel,icolumn)) / &
                dzsno_pr(ilevel,icolumn)
            else
                snowden = 0.0_r8
            endif

            ! Calculate the fraction of the SWE in each active layer
            if (h2osno_pr(icolumn) > 0.0_r8 .and. repartition_swe==1) then
                wt_swe = (h2oliq_pr(ilevel,icolumn) + h2oice_pr(ilevel,icolumn)) / &
                          h2osno_pr(icolumn)
            ! Force all SWE to bottom layer if repartition_swe = 2
            ! *AND* if level is bottom snow layer
            elseif (h2osno_pr(icolumn) > 0.0_r8 .and. repartition_swe==2 &
                    .and. ilevel==nlevsno ) then
                wt_swe = 1.0_r8
            ! For all other cases do not allow snow layer to be adjusted
            else
                wt_swe = 0.0_r8
            endif
 
            ! Calculate the fraction of liquid (and ice) water in this layer  
            if ((h2oliq_pr(ilevel,icolumn) + h2oice_pr(ilevel,icolumn)) > 0.0_r8) then
                 wt_liq = h2oliq_pr(ilevel,icolumn) / &
                (h2oliq_pr(ilevel,icolumn) + h2oice_pr(ilevel,icolumn))
                 wt_ice = 1.0_r8 - wt_liq
            else
                 wt_liq = 0.0_r8
                 wt_ice = 0.0_r8
            endif
            
            ! Calculate the increment for each layer for SWE, assuming identical
            ! layer distribution (liq,ice) as the prior 
            ! If there is no increment of column SWE (h2osno_pr-h2osno_po=0) then force gain =0
            ! and the column will not be re-partitioned for any layers
            if (abs(h2osno_po(icolumn) - h2osno_pr(icolumn)) > 0.0_r8) then 
               gain_h2osno = (h2osno_po(icolumn) - h2osno_pr(icolumn)) * wt_swe
               gain_h2oliq = gain_h2osno * wt_liq
               gain_h2oice = gain_h2osno * wt_ice
            else
               gain_h2osno = 0.0_r8
               gain_h2oliq = 0.0_r8
               gain_h2oice = 0.0_r8
            endif
            
            ! Calculate increment for each layer depth, using prior snow density
            if (snowden > 0.0_r8 ) then
                gain_dzsno(ilevel,icolumn) = gain_h2osno / snowden
            else
                gain_dzsno(ilevel,icolumn) = 0.0_r8
            endif
            
            ! Apply the increment for liquid, ice and depth for each layer.
            h2oliq_po(ilevel,icolumn) = h2oliq_pr(ilevel,icolumn) + gain_h2oliq
            h2oice_po(ilevel,icolumn) = h2oice_pr(ilevel,icolumn) + gain_h2oice

            if (h2oliq_po(ilevel,icolumn) < 0.0_r8) h2oliq_po(ilevel,icolumn) = 0.00000001_r8
            if (h2oice_po(ilevel,icolumn) < 0.0_r8) h2oice_po(ilevel,icolumn) = 0.00000001_r8

            ! Important to update snow layer dimensions because CLM code relies
            ! on snow layer thickness for compaction/aggregation snow algorithm
            ! to function properly            
            
            ! Only update snow layer dimensions if column H2OSNO has been adjusted 
            if (abs(h2osno_po(icolumn) - h2osno_pr(icolumn)) > 0.0_r8) then
                
                dzsno_po(ilevel,icolumn) =  dzsno_pr(ilevel,icolumn) + gain_dzsno(ilevel,icolumn)
                if (dzsno_po(ilevel,icolumn) < 0.0_r8) dzsno_po(ilevel,icolumn) = 0.00000001_r8

                ! For consistency  with updated dzsno_po (thickness)
                ! also update zsno_po (middle depth) and zisno (top interface depth)
             
                zisno_po(ilevel,icolumn) = sum(dzsno_po(ilevel:nlevsno,icolumn))*-1.0_r8
            
                if (ilevel ==  nlevsno) then
                    zsno_po(ilevel,icolumn) = zisno_po(ilevel,icolumn)/2.0_r8 
                else
                    zsno_po(ilevel,icolumn) = sum(zisno_po(ilevel:ilevel+1,icolumn))/2.0_r8
                endif
            endif
            
            if (abs(h2osno_po(icolumn) - h2osno_pr(icolumn)) &
                > 0.0_r8 .and. verbose > 2) then
               ! Diagnostic output for active snow columns in which SWE is updated
               ! by DART.  These columns undergo re-partitioning.
               ! column,level,active snow layers,SWE,ice mass,liq mass,
               ! thickness,interface,middle 
               ! PRIOR: icolumn,ilevel,snlsno,h2osno_pr,h2oice_pr,h2oliq_pr,
               !        dzsno_pr, zisno_pr, zsno_pr
               ! POSTERIOR: icolumn,ilevel,snlsno,h2osno_po,h2oice_po,h2oliq_po,
               !             dzsno_po, zisno_po, zsno_po

               call error_handler(E_MSG,routine,'  ',source)
               if (repartition_swe==1) then               
                  write (string1,*)'Repartitioning SWE adjustment to all snow layers'
                  call error_handler(E_MSG,routine,string1)
               elseif (repartition_swe==2) then
                  write (string1,*)'Repartitioning SWE adjustment to bottom snow layer only'
                  call error_handler(E_MSG,routine,string1)
               endif
                             
               write (string1,*)'PRIOR: ',icolumn,ilevel,snlsno(icolumn),h2osno_pr(icolumn),&
                      h2oice_pr(ilevel,icolumn),h2oliq_pr(ilevel,icolumn),&
                      dzsno_pr(ilevel,icolumn),zisno_pr(ilevel,icolumn),&
                      zsno_pr(ilevel,icolumn)
               write (string2,*)'POST : ',icolumn,ilevel,snlsno(icolumn),h2osno_po(icolumn),&        
                      h2oice_po(ilevel,icolumn),h2oliq_po(ilevel,icolumn),&
                      dzsno_po(ilevel,icolumn),zisno_po(ilevel,icolumn),&
                      zsno_po(ilevel,icolumn)
               call error_handler(E_MSG,routine,string1,source,text2=string2)
            endif
             
         enddo
      
      endif
      
      ! Update the total snow depth variable to match updates to layers 
      ! Only sum the gain of layers that are known to be active
      ! and recieved an H2OSNO adjustment.  This eliminates operating
      ! on columns that do not require re-partitioning
      if (abs(h2osno_po(icolumn) - h2osno_pr(icolumn)) & 
          > 0.0_r8 .and. snlsno(icolumn) < 0) then
         snowdp_po(icolumn) = snowdp_pr(icolumn) + sum(gain_dzsno(nlevsno+1+snlsno(icolumn):nlevsno,icolumn))
      else
         snowdp_po(icolumn) = snowdp_pr(icolumn)
      endif

      if (snowdp_po(icolumn) < 0.0_r8) snowdp_po(icolumn) = 0.00000001_r8

   endif

enddo PARTITION


! Update the 'dart_array' (dart_posterior.nc) with the repartitioned values
 dart_SNOWDP = snowdp_po
 dart_DZSNO  = dzsno_po
 dart_ZSNO   = zsno_po
 dart_ZISNO  = zisno_po
! Important to update only the manually repartitioned above surface (snow) layers
! Keep the subsurface layer values that come from (unpartitioned) DART posterior
 dart_H2OLIQ(1:nlevsno,:) =h2oliq_po(1:nlevsno,:)  
 dart_H2OICE(1:nlevsno,:) =h2oice_po(1:nlevsno,:)

! Update the 'clm_array' with repartitioned 'dart_array' similar to 
! subroutine replace_values_2D, but with snow variables

! UPDATE SNOW DEPTH 
if (nc_variable_exists(ncid_dart, 'SNOW_DEPTH')) then
   VarID=get_varid_from_varname(dom_id, 'SNOW_DEPTH')
else
   VarID=get_varid_from_varname(dom_id, 'SNOWDP')
endif

! Identify location of missing_values and FillValues to prevent updates (masking)  
! Is this right using 'get_has_missing_value' instead of 'get_has_FillValue' ??
! I am following the same usage as in subroutine 'replace_values_2D'
if (get_has_missing_value(dom_id,VarID)) call get_missing_value(dom_id,VarID,special)
if (get_has_missing_value(dom_id,VarID)) call get_FillValue(    dom_id,VarID,special)

if (nc_variable_exists(ncid_dart, 'SNOW_DEPTH')) then
   call Compatible_Variables('SNOW_DEPTH', ncid_dart, ncid_clm, varsize)
   where(dart_SNOWDP /= special) clm_SNOWDP = dart_SNOWDP
   call nc_put_variable(ncid_clm, 'SNOW_DEPTH', clm_SNOWDP, routine)
else
   call Compatible_Variables('SNOWDP', ncid_dart, ncid_clm, varsize)
   where(dart_SNOWDP /= special) clm_SNOWDP = dart_SNOWDP
   call nc_put_variable(ncid_clm, 'SNOWDP', clm_SNOWDP, routine)
endif

! UPDATE SNOW LAYER THICKNESS
VarID=get_varid_from_varname(dom_id, 'DZSNO')
! Identify location of missing_values and FillValues to prevent updates (masking)  
if (get_has_missing_value(dom_id,VarID)) call get_missing_value(dom_id,VarID,special)
if (get_has_missing_value(dom_id,VarID)) call get_FillValue(    dom_id,VarID,special)
! Make sure variables in both files are identical in shape, etc.
call Compatible_Variables('DZSNO', ncid_dart, ncid_clm, varsize)
where(dart_DZSNO /= special) clm_DZSNO = dart_DZSNO
call nc_put_variable(ncid_clm, 'DZSNO', clm_DZSNO, routine)

VarID=get_varid_from_varname(dom_id, 'ZSNO')
! Identify location of missing_values and FillValues to prevent updates (masking)  
if (get_has_missing_value(dom_id,VarID)) call get_missing_value(dom_id,VarID,special)
if (get_has_missing_value(dom_id,VarID)) call get_FillValue(    dom_id,VarID,special)
! Make sure variables in both files are identical in shape, etc.
call Compatible_Variables('ZSNO', ncid_dart, ncid_clm, varsize)
where(dart_ZSNO /= special) clm_ZSNO = dart_ZSNO
call nc_put_variable(ncid_clm, 'ZSNO', clm_ZSNO, routine)

VarID=get_varid_from_varname(dom_id, 'ZISNO')
! Identify location of missing_values and FillValues to prevent updates (masking)  
if (get_has_missing_value(dom_id,VarID)) call get_missing_value(dom_id,VarID,special)
if (get_has_missing_value(dom_id,VarID)) call get_FillValue(    dom_id,VarID,special)
! Make sure variables in both files are identical in shape, etc.
call Compatible_Variables('ZISNO', ncid_dart, ncid_clm, varsize)
where(dart_ZISNO /= special) clm_ZISNO = dart_ZISNO
call nc_put_variable(ncid_clm, 'ZISNO', clm_ZISNO, routine)

! UPDATE LIQUID LAYER
VarID=get_varid_from_varname(dom_id, 'H2OSOI_LIQ')
! Identify location of missing_values and FillValues to prevent updates (masking)  
if (get_has_missing_value(dom_id,VarID)) call get_missing_value(dom_id,VarID,special)
if (get_has_missing_value(dom_id,VarID)) call get_FillValue(    dom_id,VarID,special)
! Make sure variables in both files are identical in shape, etc.
call Compatible_Variables('H2OSOI_LIQ', ncid_dart, ncid_clm, varsize)
where(dart_H2OLIQ /= special) clm_H2OLIQ = dart_H2OLIQ
call nc_put_variable(ncid_clm, 'H2OSOI_LIQ', clm_H2OLIQ, routine)

! UPDATE ICE LAYER
VarID=get_varid_from_varname(dom_id, 'H2OSOI_ICE')
! Identify location of missing_values and FillValues to prevent updates (masking)  
if (get_has_missing_value(dom_id,VarID)) call get_missing_value(dom_id,VarID,special)
if (get_has_missing_value(dom_id,VarID)) call get_FillValue(    dom_id,VarID,special)
! Make sure variables in both files are identical in shape, etc.
call Compatible_Variables('H2OSOI_ICE', ncid_dart, ncid_clm, varsize)
where(dart_H2OICE /= special) clm_H2OICE = dart_H2OICE
call nc_put_variable(ncid_clm, 'H2OSOI_ICE', clm_H2OICE, routine)


end subroutine update_snow

!------------------------------------------------------------------
!>

subroutine Compatible_Variables(varname, ncid_dart, ncid_clm, varsize)

character(len=*), intent(in)  :: varname
integer,          intent(in)  :: ncid_dart
integer,          intent(in)  :: ncid_clm
integer,          intent(out) :: varsize(:)

character(len=*),parameter :: routine = 'Compatible_Variables'

integer ::  clm_numdims,  clm_dimlengths(NF90_MAX_VAR_DIMS)
integer :: dart_numdims, dart_dimlengths(NF90_MAX_VAR_DIMS)

character(len=NF90_MAX_NAME) :: clm_dimnames( NF90_MAX_VAR_DIMS)
character(len=NF90_MAX_NAME) :: dart_dimnames(NF90_MAX_VAR_DIMS)

call nc_get_variable_num_dimensions( ncid_clm, varname, clm_numdims,    routine)
call nc_get_variable_dimension_names(ncid_clm, varname, clm_dimnames,   routine)
call nc_get_variable_size(           ncid_clm, varname, clm_dimlengths, routine)

call nc_get_variable_num_dimensions( ncid_dart, varname, dart_numdims,    routine)
call nc_get_variable_dimension_names(ncid_dart, varname, dart_dimnames,   routine)
call nc_get_variable_size(           ncid_dart, varname, dart_dimlengths, routine)

if (clm_numdims /= dart_numdims) then
   write(string1,*)'Variable ranks not identical and must be.'
   write(string2,*)'CLM  "'//trim(varname)//'" has rank ',clm_numdims
   write(string3,*)'DART "'//trim(varname)//'" has rank ',dart_numdims
   call error_handler(E_ERR, routine, string1, source, text2=string2, text3=string3)
endif

! Ensure CLM netCDF variable is conformable with DART variable.

DimCheck : do irank = 1,clm_numdims

   if ( clm_dimlengths(irank) /= dart_dimlengths(irank) ) then
      write(string1,*)'Variable dimension lengths not identical and must be.'
      write(string2,*)'CLM  "'//trim(varname)//'" dim=',irank, ', length=',clm_dimlengths(irank)
      write(string3,*)'DART "'//trim(varname)//'" dim=',irank, ', length=',dart_dimlengths(irank)
      call error_handler(E_ERR, routine, string1, source, text2=string2, text3=string3)
   endif

   ! Do we want to check dimension name too?
   if ( clm_dimnames(irank) /= dart_dimnames(irank) ) then
      write(string1,*)'Variable dimension names not identical and must be.'
      write(string2,*)'CLM  "'//trim(varname)//'" dim=',irank, ', name=',clm_dimnames(irank)
      write(string3,*)'DART "'//trim(varname)//'" dim=',irank, ', name=',dart_dimnames(irank)
      call error_handler(E_ERR, routine, string1, source, text2=string2, text3=string3)
   endif

   varsize(irank) = clm_dimlengths(irank)

enddo DimCheck

if (verbose > 0) then
      write(string1,*)'Replacing indeterminate values ...'
      write(string2,*)'CLM  "'//trim(varname)//'" shape = ', clm_dimlengths(1:clm_numdims)
      write(string3,*)'DART "'//trim(varname)//'" shape = ',dart_dimlengths(1:clm_numdims)
      call error_handler(E_MSG, routine, string1, source, text2=string2, text3=string3)
endif

end subroutine Compatible_Variables



!===============================================================================

end program dart_to_clm
