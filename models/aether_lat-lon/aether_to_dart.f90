! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

program aether_to_dart

!----------------------------------------------------------------------
! purpose: Transform the Aether model restarts into a DART filter_input.nc.
!
! method: Read aether "restart" files of model state (multiple files, 
!         one block per aether mpi task)
!         Reform fields into a DART netcdf file
!
! USAGE:  The aether restart dirname and output filename are read from 
!         the aether_to_dart_nml namelist.
!         <edit input.nml:aether_to_dart_nml>
!----------------------------------------------------------------------
! Converts Aether restart files to a netCDF file

use types_mod, only : r4, MISSING_I, vtablenamelength

use time_manager_mod, only: time_type

use utilities_mod, only :                                    &
    finalize_utilities, error_handler, E_ERR, E_MSG, E_WARN, &
    initialize_utilities, do_output

use default_model_mod, only : write_model_time

use transform_state_mod, only :                              &
    static_init_blocks, aether_name_to_dart,                 &
    nghost, open_block_file, aether_restart_dirname,         &
    VT_ORIGININDX, VT_VARNAMEINDX, nvar_neutral,  nvar_ion,  &
    nx_per_block, ny_per_block, nz_per_block,                &
    nblocks_lon, nblocks_lat, variables,                     &
    lats, levs, lons, debug, state_time,                     &
    block_file_name, nlat, nlon, nlev, purge_chars

use netcdf_utilities_mod, only :                                  &
    nc_create_file, nc_close_file,                                &
    nc_begin_define_mode, nc_end_define_mode,                     &
    nc_define_dimension,                                          &
    nc_add_global_attribute, nc_add_global_creation_time,         &
    nc_get_attribute_from_variable, nc_add_attribute_to_variable, &
    nc_define_real_variable, nc_define_real_scalar,               &
    nc_get_variable, nc_put_variable,                             &
    nc_synchronize_file

implicit none

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer :: member = MISSING_I, &
           num_args, ncid
character(len=3)   :: char_mem
character(len=31)  :: filter_io_root = 'filter_input'
character(len=64)  :: filter_io_file = ''
character(len=512) :: error_string_1, error_string_2
character(len=31),  parameter :: progname = 'aether_to_dart'
character(len=256), parameter :: source   = 'aether_lat-lon/aether_to_dart.f90'

character(len=4), parameter :: LEV_DIM_NAME = 'alt'
character(len=4), parameter :: LAT_DIM_NAME = 'lat'
character(len=4), parameter :: LON_DIM_NAME = 'lon'
character(len=4), parameter :: TIME_DIM_NAME = 'time'

character(len=4), parameter :: LEV_VAR_NAME = 'alt'
character(len=4), parameter :: LAT_VAR_NAME = 'lat'
character(len=4), parameter :: LON_VAR_NAME = 'lon'
character(len=4), parameter :: TIME_VAR_NAME = 'time'

!======================================================================

call initialize_utilities(progname)

!----------------------------------------------------------------------
! Get the ensemble member
!----------------------------------------------------------------------
num_args = command_argument_count()
if (num_args == 0) then
   write(error_string_1,*) 'Usage: ./aether_to_dart member_number (0-based)' 
   call error_handler(E_ERR, progname, error_string_1)
endif

call get_command_argument(1,char_mem)
read(char_mem,'(I3)') member

!----------------------------------------------------------------------
! Convert the files
!----------------------------------------------------------------------

call static_init_blocks(member)

! Must be after static_init_blocks, which provides filter_io_root from the namelist.
write(filter_io_file,'(2A, I0.4, A3)') trim(filter_io_root),'_', member + 1,'.nc'
call error_handler(E_MSG, '', '')
write(error_string_1,'(A,I3,2A)') 'Converting Aether member ',member, &
     ' restart files to the NetCDF file ', trim(filter_io_file)
write(error_string_2,'(3A)') ' in directory ', trim(aether_restart_dirname)
call error_handler(E_MSG, progname, error_string_1, text2=error_string_2)
call error_handler(E_MSG, '', '')

! nc_create_file does not leave define mode.
ncid = nc_create_file(trim(aether_restart_dirname)//'/'//trim(filter_io_file))
! def_fill_dimvars does leave define mode.
call def_fill_dimvars(ncid)

! Write_model_time will make a time variable, if needed, which it is not.
! state_time is read in transform_state_mod and is available by the use statement.
call write_model_time(ncid, state_time)

! Define (non-time) variables
call restarts_to_filter(ncid, member, define=.true.)

! Read and convert (non-time) variables
call restarts_to_filter(ncid, member, define=.false.)
! subr. called by this routine closes the file only if define = .true.
call nc_close_file(ncid)

call error_handler(E_MSG, '', '')
write(error_string_1,'(3A)') 'Successfully converted the Aether restart files to ', &
     "'"//trim(filter_io_file)//"'"
call error_handler(E_MSG, progname, error_string_1)
call error_handler(E_MSG, '', '')
    
! end - close the log, etc
call finalize_utilities()

!-----------------------------------------------------------------------
contains

!-----------------------------------------------------------------------
! Open all restart files (blocks x {neutrals,ions}) for 1 member 
! and transfer the requested variable contents to the filter input file.
! This is called with 'define' = 
! .true.  define variables in the file or 
! .false. transfer the data from restart files to a filter_inpu.nc file.

subroutine restarts_to_filter(ncid_output, member, define)

integer,          intent(in)  :: ncid_output, member
logical,          intent(in)  :: define

integer :: ib, jb, ib_loop, jb_loop

if (define) then
   ! if define, run one block.
   ! the block_to_filter_io call defines the variables in the whole domain netCDF file.
   ib_loop = 1
   jb_loop = 1
   call nc_begin_define_mode(ncid_output)
else
   ! if not define, run all blocks.
   ! the block_to_filter_io call adds the (ib,jb) block to the netCDF variables 
   ! in order to make a file containing the data for all the blocks.
   ib_loop = nblocks_lon
   jb_loop = nblocks_lat
end if

do jb = 1, jb_loop
   do ib = 1, ib_loop
      call block_to_filter_io(ncid_output, ib, jb, member, define)
   enddo
enddo

if (define) then
   call nc_end_define_mode(ncid_output)
endif
    
end subroutine restarts_to_filter

!-----------------------------------------------------------------------
! Transfer variable data from a block restart file to the filter_input.nc file. 
! It's called with 2 modes:
! define = .true.  define the NC variables in the filter_input.nc 
! define = .false. write the data from a block to the NC file using write_filter_io.

subroutine block_to_filter_io(ncid_output, ib, jb, member, define)

integer,          intent(in) :: ncid_output
integer,          intent(in) :: ib, jb
integer,          intent(in) :: member
logical,          intent(in) :: define

real(r4), allocatable :: temp1d(:), temp2d(:,:), temp3d(:,:,:)
! real(r4), allocatable :: alt1d(:), density_ion_e(:,:,:)
integer               :: ivar, nb, ncid_input
! TEC? integer               :: maxsize
!      logical               :: no_idensity
!      real(r4)              :: temp0d 
character(len=32)     :: att_val 
character(len=128)    :: file_root 
character(len=256)    :: filename
character(len=vtablenamelength) :: varname, dart_varname

character(len=*), parameter :: routine = 'block_to_filter_io'

! The block number, as counted in Aether.
! Lower left is 0, increase to the East, then 1 row farther north, West to East.
nb = (jb - 1) * nblocks_lon + ib - 1

! a temp array large enough to hold any of the
! Lon,Lat or Alt array from a block plus ghost cells
allocate(temp1d(1-nghost:max(nx_per_block, ny_per_block, nz_per_block) + nghost))

! treat alt specially since we want to derive TEC here
! TODO: See density_ion_e too.
! allocate( alt1d(1-nghost:max(nx_per_block, ny_per_block, nz_per_block) + nghost))

! temp array large enough to hold any 2D field 
allocate(temp2d(1-nghost:ny_per_block+nghost, &
                1-nghost:nx_per_block+nghost))

! TODO: We need all altitudes, but there might be vertical blocks in the future.
!       But there would be no vertical halos.
!       Make transform_state_mod: zcount adapt to whether there are blocks.
!       Temp needs to have C-ordering, which is what the restart files have.
! temp array large enough to hold 1 species, temperature, etc
allocate(temp3d(1:nz_per_block, &
                1-nghost:ny_per_block+nghost, &
                1-nghost:nx_per_block+nghost))

! TODO: Waiting for e- guidance from Aaron.
! save density_ion_e to compute TEC
! allocate(density_ion_e(1:nz_per_block, &
!                        1-nghost:ny_per_block+nghost, &
!                        1-nghost:nx_per_block+nghost))

! TODO: Aether gives a unique name to each (of 6) velocity components.
!       Do we want to use a temp4d array to handle them?  
!       They are independent variables in the block files (and state).
! ! temp array large enough to hold velocity vect, etc
! maxsize = max(3, nvar_ion)
! allocate(temp4d(1-nghost:nx_per_block+nghost, &
!                 1-nghost:ny_per_block+nghost, &
!                 1-nghost:nz_per_block+nghost, maxsize))


! TODO; Does Aether need a replacement for these Density fields?  Yes.
!       But they are probably read by the loops below.
!       Don't need to fetch index because Aether has NetCDF restarts,
!       so just loop over the field names to read.
! 
! ! assume we could not find the electron density for VTEC calculations
! no_idensity = .true.
! 
! if (inum > 0) then
!    ! one or more items in the state vector need to replace the
!    ! data in the output file.  loop over the index list in order.
!    j = 1
! ! TODO:   electron density is not in the restart files, but it's needed for TEC
!           In Aether they will be from an ions file, but now only from an output file (2023-10-30).
!           Can that be handled like the neutrals and ions files, using variables(VT_ORIGININDX,:)
!           to build an output file name?  Are outputs in block form?
!                ! save the electron density for TEC computation
!                density_ion_e(:,:,:) = temp3d(:,:,:)

! Handle the 2 restart file types (ions and neutrals).
! Each field has a file type associated with it: variables(VT_ORIGININDX,f_index)

file_root = variables(VT_ORIGININDX,1)
write(filename,'(A,"/",A)') trim(aether_restart_dirname), &
     trim(block_file_name(trim(file_root), member, nb))
ncid_input = open_block_file(filename, 'read')

do ivar = 1, nvar_neutral
   ! The nf90 functions cannot read the variable names with the '\'s in them.
   varname = purge_chars(trim(variables(VT_VARNAMEINDX,ivar)), '\', plus_minus=.false.)
   if (debug >= 100 .and. do_output()) print*, routine,'varname = ', varname
   ! Translate the Aether field name into a CF-compliant DART field name.
   dart_varname = aether_name_to_dart(varname)

   ! TODO: Given the subroutine name, perhaps these definition sections should be 
   !       one call higher up, with the same loop around it.
   if (define) then
      ! Define the variable in the filter_input.nc file (the output from this program).
      ! The calling routine entered define mode.

      if (debug > 10 .and. do_output()) then
         write(error_string_1,'(A,I0,2A)') 'Defining ivar = ', ivar,':', dart_varname
         call error_handler(E_MSG, routine, error_string_1, source)
      end if
   
      call nc_define_real_variable(ncid_output, dart_varname, &
           (/ LEV_DIM_NAME, LAT_DIM_NAME, LON_DIM_NAME/) )
      call nc_get_attribute_from_variable(ncid_input, varname, 'units', att_val, routine)
      call nc_add_attribute_to_variable(ncid_output, dart_varname, 'units', att_val, routine)

   else if (file_root == 'neutrals') then
   ! Read 3D array and extract the non-halo data of this block.
      call nc_get_variable(ncid_input, varname, temp3d, context=routine)
      call write_filter_io(temp3d, dart_varname, ib, jb, ncid_output)
   else
      write(error_string_1,'(A,I3,A)') 'Trying to read neutrals, but variables(', &
           VT_ORIGININDX,ivar , ') /= "neutrals"'
      call error_handler(E_ERR, routine, error_string_1, source)
   endif

enddo
call nc_close_file(ncid_input)

file_root = variables(VT_ORIGININDX,nvar_neutral+1)
write(filename,'(A,"/",A)') trim(aether_restart_dirname), &
     trim(block_file_name(trim(file_root), member, nb))
ncid_input = open_block_file(filename, 'read')

do ivar = nvar_neutral +1, nvar_neutral + nvar_ion
   ! Purging \ from aether name.
   varname = purge_chars(trim(variables(VT_VARNAMEINDX,ivar)), '\', plus_minus=.false.)
   dart_varname = aether_name_to_dart(varname)

   if (define) then

      if (debug > 10 .and. do_output()) then
         write(error_string_1,'(A,I0,2A)') 'Defining ivar = ', ivar,':', dart_varname
         call error_handler(E_MSG, routine, error_string_1, source)
      end if
   
      call nc_define_real_variable(ncid_output, dart_varname, &
                                   (/ LEV_DIM_NAME, LAT_DIM_NAME, LON_DIM_NAME/) )
      call nc_get_attribute_from_variable(ncid_input, varname, 'units', att_val, routine)
      call nc_add_attribute_to_variable(ncid_output, dart_varname, 'units', att_val, routine)
      print*, routine,': defined ivar, dart_varname, att = ', &
              ivar, trim(dart_varname), trim(att_val)

   else if (file_root == 'ions') then
      call nc_get_variable(ncid_input, varname, temp3d, context=routine)
      call write_filter_io(temp3d, dart_varname, ib, jb, ncid_output)
   else
      write(error_string_1,'(A,I3,A)') 'Trying to read ions, but variables(', &
           VT_ORIGININDX,ivar , ') /= "ions"'
      call error_handler(E_ERR, routine, error_string_1, source)
   endif

enddo

! Leave file open if fields were just added (define = .false.),
! so that time can be added.
if (define) call nc_close_file(ncid_input)

! TODO: Does Aether need TEC to be calculated? Yes
! ! add the VTEC as an extended-state variable
! ! NOTE: This variable will *not* be written out to the Aether restart files 
! 
! if (no_idensity) then
!    write(error_string_1,*) 'Cannot compute the VTEC without the electron density'
!    call error_handler(E_ERR, routine, error_string_1, source)
! end if
! 
!       temp2d = 0._r8
!       ! compute the TEC integral
!       do i =1,nz_per_block-1 ! approximate the integral over the altitude as a sum of trapezoids
!          ! area of a trapezoid: A = (h2-h1) * (f2+f1)/2
!          temp2d(:,:) = temp2d(:,:) + ( alt1d(i+1)-alt1d(i) )  * &
!                        ( density_ion_e(:,:,i+1)+density_ion_e(:,:,i) ) /2.0_r8
!       end do  
!       ! convert temp2d to TEC units
!       temp2d = temp2d/1e16_r8
!    call write_block_to_filter2d(temp2d, ivals(1), block, ncid, define) 

! TODO: Does Aether need f10_7 to be calculated or processed? Yes
! !gitm_index = get_index_start(domain_id, 'VerticalVelocity')
! call get_index_from_gitm_varname('f107', inum, ivals)
! if (inum > 0) then
!   call write_block_to_filter0d(temp0d, ivals(1), ncid, define) !see comments in the body of the subroutine
! endif
! 

deallocate(temp1d, temp2d, temp3d)
! deallocate(alt1d, density_ion_e)
   
end subroutine block_to_filter_io

!-----------------------------------------------------------------------
! Open all restart files (neutrals,ions) for a block and read in the requested data items.
! The write_filter_io calls will write the data to the filter_input.nc.

subroutine write_filter_io(data3d, varname, ib, jb, ncid)

real(r4), intent(in) :: data3d(1:nz_per_block, &
                               1-nghost:ny_per_block+nghost, &
                               1-nghost:nx_per_block+nghost)

character(len=vtablenamelength), intent(in) :: varname
integer,  intent(in) :: ib, jb
integer,  intent(in) :: ncid

integer :: starts(3)

character(len=*), parameter :: routine = 'write_filter_io'

! write(varname,'(A)') trim(variables(VT_VARNAMEINDX,ivar))

! to compute the start, consider (ib-1)*nx_per_block+1
starts(1) = 1
starts(2) = (jb-1) * ny_per_block + 1
starts(3) = (ib-1) * nx_per_block + 1

call nc_put_variable(ncid, varname, &
                     data3d(1:nz_per_block,1:ny_per_block,1:nx_per_block), &
                     context=routine, nc_start=starts, &
                     nc_count=(/nz_per_block,ny_per_block,nx_per_block/))
   
end subroutine write_filter_io

!-----------------------------------------------------------------------
! Add dimension variable contents to the file.

subroutine def_fill_dimvars(ncid)

integer, intent(in) :: ncid

character(len=*), parameter :: routine = 'def_fill_dimvars'

! File is still in define mode from nc_create_file
! call nc_begin_define_mode(ncid)

! Global atts for aether_to_dart and dart_to_aether.
call nc_add_global_creation_time(ncid, routine)
call nc_add_global_attribute(ncid, "model_source", source, routine)
call nc_add_global_attribute(ncid, "model", "aether", routine)

! define grid dimensions
call nc_define_dimension(ncid, trim(LEV_DIM_NAME),  nlev,  routine)
call nc_define_dimension(ncid, trim(LAT_DIM_NAME),  nlat,  routine)
call nc_define_dimension(ncid, trim(LON_DIM_NAME),  nlon,  routine)

! define grid variables
! z
call nc_define_real_variable(     ncid, trim(LEV_VAR_NAME), (/ trim(LEV_DIM_NAME) /), routine)
call nc_add_attribute_to_variable(ncid, trim(LEV_VAR_NAME), 'units',     'm', routine)
call nc_add_attribute_to_variable &
     (ncid, trim(LEV_VAR_NAME), 'long_name', 'height above mean sea level', routine)

! latitude
call nc_define_real_variable(     ncid, trim(LAT_VAR_NAME), (/ trim(LAT_DIM_NAME) /),  routine)
call nc_add_attribute_to_variable(ncid, trim(LAT_VAR_NAME), 'units',     'degrees_north', routine)
call nc_add_attribute_to_variable(ncid, trim(LAT_VAR_NAME), 'long_name', 'latitude', routine)

! longitude
call nc_define_real_variable(     ncid, trim(LON_VAR_NAME), (/ trim(LON_VAR_NAME) /), routine)
call nc_add_attribute_to_variable(ncid, trim(LON_VAR_NAME), 'units', 'degrees_east', routine)
call nc_add_attribute_to_variable(ncid, trim(LON_VAR_NAME), 'long_name', 'longitude',  routine)

! Dimension 'time' will no longer be created by write_model_time,
! or by nc_define_unlimited_dimension.  It will be a scalar variable.
! time
call nc_define_real_scalar(       ncid, trim(TIME_VAR_NAME), routine)
call nc_add_attribute_to_variable(ncid, trim(TIME_VAR_NAME), 'calendar', 'gregorian', routine)
call nc_add_attribute_to_variable &
     (ncid, trim(TIME_VAR_NAME), 'units', 'days since 1601-01-01 00:00:00', routine)
call nc_add_attribute_to_variable &
     (ncid, trim(TIME_VAR_NAME), 'long_name', 'gregorian_days',  routine)

call nc_end_define_mode(ncid)

call nc_put_variable(ncid, trim(LEV_VAR_NAME),  levs,  routine)
call nc_put_variable(ncid, trim(LAT_VAR_NAME),  lats,  routine)
call nc_put_variable(ncid, trim(LON_VAR_NAME),  lons,  routine)
! time will be written elsewhere.

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine def_fill_dimvars

!-----------------------------------------------------------------------
end program aether_to_dart

