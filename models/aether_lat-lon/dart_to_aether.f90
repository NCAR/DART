! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

program dart_to_aether

!----------------------------------------------------------------------
! purpose: Transform a DART filter_output.nc into the Aether model restarts.
!
! method: Read DART state netcdf file and overwrite values in Aether restart files.
!
! this version assumes that the DART grid is global and the data needs to be
! blocked into one block per Aether mpi task.  there is a different converter
! for when Aether only needs a single input/output file.
!
!----------------------------------------------------------------------

use types_mod, only : r4, MISSING_I,  MISSING_R4, vtablenamelength

use utilities_mod, only :                                    &
    finalize_utilities, error_handler, E_ERR, E_MSG, E_WARN, &
    initialize_utilities, do_output

use default_model_mod, only : write_model_time

use transform_state_mod, only :                                   &
    debug, aether_restart_dirname, nblocks_lat,                   &
    nblocks_lon, nghost, nlat, nlon, nlev,                        &
    nx_per_block, ny_per_block, nz_per_block,                     &
    nvar_ion, nvar_neutral, VT_ORIGININDX, VT_VARNAMEINDX,        &
    block_file_name, open_block_file, aether_name_to_dart,        &
    variables, purge_chars, static_init_blocks

use netcdf_utilities_mod, only :                                  &
    nc_open_file_readonly, nc_close_file,                         &
    nc_begin_define_mode, nc_end_define_mode,                     &
    nc_define_dimension,                                          &
    nc_add_global_attribute, nc_add_global_creation_time,         &
    nc_get_attribute_from_variable, nc_add_attribute_to_variable, &
    nc_define_real_variable, nc_define_real_scalar,               &
    nc_get_variable, nc_put_variable, nc_variable_exists,         &
    nc_synchronize_file, NF90_FILL_REAL

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
character(len=31),  parameter :: progname = 'dart_to_aether'
character(len=256), parameter :: source   = 'aether_lat-lon/dart_to_aether.f90'

!======================================================================

call initialize_utilities(progname)

!----------------------------------------------------------------------
! Get the ensemble member
!----------------------------------------------------------------------
num_args = command_argument_count()
if (num_args == 0) then
   write(error_string_1,*) 'Usage: ./dart_to_aether member_number (0-based)' 
   call error_handler(E_ERR, progname, error_string_1)
endif

call get_command_argument(1,char_mem)
read(char_mem,'(I3)') member

!----------------------------------------------------------------------
! Convert the files
!----------------------------------------------------------------------

call static_init_blocks(member)

write(filter_io_file,'(2A,I0.4,A3)') trim(filter_io_root),'_',member + 1,'.nc'

call error_handler(E_MSG, source, '', '')
write(error_string_1,'(3A)') 'Extracting fields from DART file ',trim(filter_io_file)
write(error_string_2,'(A,I3,2A)') 'into Aether restart member ',member, &
     ' in directory ', trim(aether_restart_dirname)
call error_handler(E_MSG, progname, error_string_1, text2=error_string_2)
call error_handler(E_MSG, '', '')

ncid = nc_open_file_readonly(trim(aether_restart_dirname)//'/'//trim(filter_io_file), source)

call filter_to_restarts(ncid, member)

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------
call error_handler(E_MSG, source,'','')
write(error_string_1,'(3A)') 'Successfully converted to the Aether restart files in directory'
write(error_string_2,'(3A)') "'"//trim(aether_restart_dirname)//"'"
call error_handler(E_MSG, source, error_string_1, source, text2=error_string_2)

call nc_close_file(ncid)
    
! end - close the log, etc
call finalize_utilities()

!-----------------------------------------------------------------------
contains 
!-----------------------------------------------------------------------
! Extract (updated) variables from a filter_output.nc file
! and write to existing block restart files.

subroutine filter_to_restarts(ncid, member)

integer, intent(in) :: member, ncid

real(r4), allocatable           :: fulldom3d(:,:,:)
character(len=64)               :: file_root
integer                         :: ivar
character(len=vtablenamelength) :: varname, dart_varname

character(len=*), parameter :: routine = 'filter_to_restarts'

! Space for full domain field (read from filter_output.nc)
! and halo around the full domain
allocate(fulldom3d(1:nlev, &
                   1-nghost:nlat+nghost, &
                   1-nghost:nlon+nghost))

! get the dirname, construct the filenames inside open_block_file

! Not all fields have halos suitable for calculating gradients.  
! These do (2023-11-8): neutrals; temperature, O, O2, N2, and the horizontal winds. 
!                       ions; none.
! The current transform_state will fill all neutral halos anyway, 
! since that's simpler and won't break the model.
! TODO: add an attribute to the variables (?) to denote whether a field 
!       should have its halo filled?
do ivar = 1, nvar_neutral
   varname = purge_chars(trim(variables(VT_VARNAMEINDX,ivar)), '\', plus_minus=.false.)
   if (debug >= 0 .and. do_output()) then
      write(error_string_1,'("varname = ",A)') trim(varname)
      call error_handler(E_MSG, routine, error_string_1, source)
   endif
   dart_varname = aether_name_to_dart(varname)

   file_root = trim(variables(VT_ORIGININDX,ivar))
   if (trim(file_root) == 'neutrals') then
      ! This parameter is available through the `use netcdf` command.
      fulldom3d = NF90_FILL_REAL
      
      call nc_get_variable(ncid, dart_varname, fulldom3d(1:nlev,1:nlat,1:nlon), &
                           context=routine)
      ! Copy updated field values to full domain halo.
      ! Block domains+halos will be easily read from this.
      call add_halo_fulldom3d(fulldom3d)

      call filter_io_to_blocks(fulldom3d, varname, file_root, member)
   else
      write(error_string_1,'(3A)') "file_root of varname = ",trim(varname), &
           ' expected to be "neutrals"'
      call error_handler(E_ERR, routine, error_string_1, source)
   endif

enddo

do ivar = nvar_neutral + 1, nvar_neutral + nvar_ion
   varname = purge_chars(trim(variables(VT_VARNAMEINDX,ivar)), '\', plus_minus=.false.)
   dart_varname = aether_name_to_dart(varname)

   file_root = trim(variables(VT_ORIGININDX,ivar))
   if (debug >= 0 .and. do_output()) then
      write(error_string_1,'("varname, dart_varname, file_root = ",3(2x,A))') &
             trim(varname), trim(dart_varname), trim(file_root)
      call error_handler(E_MSG, routine, error_string_1, source)
   endif

   if (trim(file_root) == 'ions') then
      fulldom3d = NF90_FILL_REAL
      call nc_get_variable(ncid, dart_varname, fulldom3d(1:nlev,1:nlat,1:nlon), &
                           context=routine)
      ! 2023-11: ions do not have real or used data in their halos.
      !          Make this clear by leaving the halos filled with MISSING_R4
      !          TODO: Will this be translated into NetCDF missing_value?
      ! call add_halo_fulldom3d(fulldom3d)

      call filter_io_to_blocks(fulldom3d, varname, file_root, member)

   else
      write(error_string_1,'(3A)') "file_root of varname = ",trim(varname), &
           ' expected to be "ions"'
      call error_handler(E_ERR, routine, error_string_1, source)
   endif
enddo

deallocate(fulldom3d)
   
end subroutine filter_to_restarts


!-----------------------------------------------------------------------
! Copy updated data from the full domain into the halo regions,
! in preparation for extracting haloed blocks into the block restart files.
! First, the halos past the East and West edges are taken from the wrap-around points.
! Then, the halos beyond the edge latitudes in the North and South 
! are taken by reaching over the pole to a longitude that's half way around the globe.
! This is independent of the number of blocks.

subroutine add_halo_fulldom3d(fulldom3d)

! Space for full domain field (read from filter_output.nc)
! and halo around the full domain
real(r4), intent(inout) :: fulldom3d(1:nz_per_block,       &
                                     1-nghost:nlat+nghost, &
                                     1-nghost:nlon+nghost)  

integer               :: g, i, j, haflat, haflon
real(r4), allocatable :: normed(:,:)
character(len=16)     :: debug_format

character(len=*), parameter :: routine = 'add_halo_fulldom3d'

! An array for debugging by renormalizing an altitude of fulldom3d.
allocate(normed(1-nghost:nlat+nghost, &
                1-nghost:nlon+nghost))

haflat = nlat / 2
haflon = nlon / 2

do g = 1,nghost
   ! left; reach around the date line.
   !         There's no data at the ends of the halos for this copy.
   fulldom3d  (:,1:nlat,     1-g) &
   = fulldom3d(:,1:nlat,nlon+1-g)

   ! right
   fulldom3d  (:,1:nlat,nlon+g) &
   = fulldom3d(:,1:nlat,g)

   ! bottom; reach over the S Pole for halo values.
   !         There is data at the ends of the halos for these.)

   fulldom3d  (:, 1-g ,1-nghost       :haflon) &
   = fulldom3d(:,   g ,1-nghost+haflon:nlon)
   fulldom3d  (:, 1-g ,haflon+1:nlon) &
   = fulldom3d(:,   g ,       1:haflon)
   ! Last 2 (halo) points on the right edge (at the bottom)
   fulldom3d  (:, 1-g ,  nlon+1:  nlon+nghost) &
   = fulldom3d(:,   g ,haflon+1:haflon+nghost)

   ! top
   fulldom3d  (:, nlat  +g, 1-nghost       :haflon) &
   = fulldom3d(:, nlat+1-g, 1-nghost+haflon:nlon)
   fulldom3d  (:, nlat  +g, haflon+1:nlon) &
   = fulldom3d(:, nlat+1-g,        1:haflon)
   ! Last 2 (halo) points on the right edge (at the top)
   fulldom3d  (:, nlat  +g,   nlon+1:  nlon+nghost) &
   = fulldom3d(:, nlat+1-g, haflon+1:haflon+nghost)
enddo

if (any(fulldom3d == MISSING_R4)) then
   error_string_1 = 'ERROR: some fulldom3d contain MISSING_R4 after halos'
   call error_handler(E_ERR, routine, error_string_1, source)
endif

! TODO: Keep halo corners check for future use?
!       Add more robust rescaling.
! Print the 4x4 arrays (corners & middle) to see whether values are copied correctly.
! Level 44 values range from 800-eps to 805.  I don't want to see the 80.
! For O+ range from 0 to 7e+11, but are close to 1.1082e+10 near the corners.
! 2023-12-20; Aaron sent new files with 54 levels.
if (debug >= 100 .and. do_output()) then
   if (fulldom3d(54,10,10) > 1.e+10_r4) then
      normed = fulldom3d(54,:,:) - 1.1092e+10_r4
      debug_format = '(3(4E10.4,2X))'
   else if (fulldom3d(54,10,10) < 1000._r4) then
      normed = fulldom3d(54,:,:) - 800._r4
      debug_format = '(3(4F10.5,2X))'
   endif
   
   ! Debug HDF5 
   write(error_string_1,'("normed_field(10,nlat+1,nlon+2) = ",3(1x,i5))') normed(nlat+1,nlon+2)
   call error_handler(E_MSG, routine, error_string_1, source)
   
   ! 17 format debug_format
   print*,'top'
   do j = nlat+2, nlat-1, -1
      write(*,debug_format) (normed(j,i), i=      -1,       2), &
                            (normed(j,i), i=haflon-1,haflon+2), &
                            (normed(j,i), i=  nlon-1,  nlon+2)
   enddo
   print*,'middle'
   do j = haflat+2, haflat-1 , -1
      write(*,debug_format) (normed(j,i), i=      -1,       2), &
                            (normed(j,i), i=haflon-1,haflon+2), &
                            (normed(j,i), i=  nlon-1,  nlon+2)
   enddo
   print*,'bottom'
   do j = 2,-1, -1
      write(*,debug_format) (normed(j,i), i=      -1,       2), &
                            (normed(j,i), i=haflon-1,haflon+2), &
                            (normed(j,i), i=  nlon-1,  nlon+2)
   enddo
endif

deallocate(normed)

end subroutine add_halo_fulldom3d

!-----------------------------------------------------------------------
! Transfer part of the full field into a block restart file.

subroutine filter_io_to_blocks(fulldom3d, varname, file_root, member)

real(r4), intent(in) :: fulldom3d(1:nz_per_block,       &
                                  1-nghost:nlat+nghost, &
                                  1-nghost:nlon+nghost  )
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: file_root
integer,          intent(in) :: member

! Don't collect velocity components (6 of them)
!   real(r4) :: temp0d 
! , temp1d(:)   ?
integer            :: ncid_output
integer            :: ib, jb, nb
integer            :: starts(3), ends(3), xcount, ycount, zcount
character(len=256) :: block_file

character(len=*), parameter :: routine = 'filter_io_to_blocks'

! a temp array large enough to hold any of the
! Lon,Lat or Alt array from a block plus ghost cells
! allocate(temp1d(1-nghost:max(nx_per_block,ny_per_block,nz_per_block)+nghost))

zcount = nz_per_block
ycount = ny_per_block + (2 * nghost)
xcount = nx_per_block + (2 * nghost)

if (debug > 0 .and. do_output()) then
   write(error_string_1,'(A,I0,A,I0,A)') 'Now putting the data for ', nblocks_lon, &
        ' blocks lon by ',nblocks_lat,' blocks lat'
   call error_handler(E_MSG, routine, error_string_1, source)
end if

starts(1) = 1
ends(1)   = nz_per_block

do jb = 1, nblocks_lat
   starts(2) = (jb - 1) * ny_per_block - nghost + 1
   ends(2)   =  jb      * ny_per_block + nghost

   do ib = 1, nblocks_lon
      starts(3) = (ib - 1) * nx_per_block - nghost + 1
      ends(3)   =  ib      * nx_per_block + nghost

      nb = (jb - 1) * nblocks_lon + ib - 1

      write(block_file,'(A,"/",A)') trim(aether_restart_dirname), &
           trim(block_file_name(trim(file_root), member, nb))
      ncid_output = open_block_file(block_file, 'readwrite')
      if (.not.nc_variable_exists(ncid_output,varname)) then
         write(error_string_1,'(4A)') 'variable ', varname, ' does not exist in ',block_file
         call error_handler(E_ERR, routine, error_string_1, source)
      endif

      if ( debug > 0 .and. do_output()) then
        write(error_string_1,'(A,3(2X,i5))') "block, ib, jb = ", nb, ib, jb
        call error_handler(E_MSG, routine, error_string_1, source)
        write(error_string_1,'(3(A,3i5))') &
             'starts = ',starts, 'ends = ',ends, '[xyz]counts = ',xcount,ycount,zcount
        call error_handler(E_MSG, routine, error_string_1, source)
      endif      

      call nc_put_variable(ncid_output, trim(varname), &
           fulldom3d(starts(1):ends(1), starts(2):ends(2), starts(3):ends(3)), &
           context=routine, nc_count=(/ zcount,ycount,xcount /) )

      call nc_close_file(ncid_output)

   enddo
enddo

! 
! TODO: ? Add f107 and Rho to the restart files
!    call read_filter_io_block0d(ncid, ivals(1), data0d)
!    if (data0d < 0.0_r8) data0d = 60.0_r8 !alex
!    write(ounit) data0d

end subroutine filter_io_to_blocks

!-----------------------------------------------------------------------
end program dart_to_aether

