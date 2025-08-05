! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

program dart_to_cice

!----------------------------------------------------------------------
! purpose: implement a 'partition function' to modify the cice state 
!          to be consistent with the states from assimilation
!
! method: Read in restart (restart with prior) and out restart (restart 
!         with posterior) written by DART after filter. 
!
! author: M Wieringa (2023) based on C Bitz (2016) and C Riedel (2022)
!----------------------------------------------------------------------

use  types_mod, only : r8
use  utilities_mod, only : initialize_utilities, finalize_utilities, &
                           find_namelist_in_file, check_namelist_read, &
                           file_exist, error_handler, E_ERR, E_MSG
use ice_postprocessing_mod, only : cice_rebalancing, area_simple_squeeze, &
                                   volume_simple_squeeze, read_cice_state_variable
use  netcdf_utilities_mod,  only : nc_open_file_readwrite, nc_put_variable, &
                                   nc_close_file, nc_open_file_readonly, &
                                   nc_get_variable_size
use  netcdf

implicit none

character(len=*), parameter :: source   = 'dart_to_cice'

!------------------------------------------------------------------
! SET NAMELIST AND ALLOCATE VARIABLES    
!------------------------------------------------------------------
character(len=256) :: dart_to_cice_input_file = 'post_filter_restart.nc'
character(len=256) :: original_cice_restart_file = 'pre_filter_restart.nc'
character(len=256) :: postprocessed_output_file = 'postprocessed_restart.nc'
character(len=128) :: postprocess = 'cice'

namelist /dart_to_cice_nml/ dart_to_cice_input_file,    &
                            original_cice_restart_file, &
                            postprocessed_output_file,  &
                            postprocess

! general variable initialization
character(len=512) :: string1

integer :: iunit, io, ncid, Ncat, Nx, Ny, nsize(3)
real(r8), allocatable :: aicen_original(:,:,:), vicen_original(:,:,:), vsnon_original(:,:,:)
real(r8), allocatable :: aicen(:,:,:), vicen(:,:,:), vsnon(:,:,:), Tsfcn(:,:,:)
real(r8), allocatable :: qice001(:,:,:), qice002(:,:,:), qice003(:,:,:), qice004(:,:,:), qice005(:,:,:), qice006(:,:,:), qice007(:,:,:), qice008(:,:,:)
real(r8), allocatable :: sice001(:,:,:), sice002(:,:,:), sice003(:,:,:), sice004(:,:,:), sice005(:,:,:), sice006(:,:,:), sice007(:,:,:), sice008(:,:,:)
real(r8), allocatable :: qsno001(:,:,:), qsno002(:,:,:), qsno003(:,:,:)

!------------------------------------------------------------------
! INITIALIZE AND PERFORM CHECKS ON FILES                
!------------------------------------------------------------------
call initialize_utilities(progname=source)

call find_namelist_in_file("input.nml", "dart_to_cice_nml", iunit)
read(iunit, nml = dart_to_cice_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_cice_nml")

if ( .not. file_exist(dart_to_cice_input_file) ) then
   write(string1,*) 'cannot open "', trim(dart_to_cice_input_file),'" for updating.'
   call error_handler(E_ERR,source,'filename not found ',trim(dart_to_cice_input_file))
endif

if ( .not. file_exist(original_cice_restart_file) ) then
   write(string1,*) 'cannot open "', trim(original_cice_restart_file),'" for reading.'
   call error_handler(E_ERR,source,'filename not found ',trim(original_cice_restart_file))
endif

!------------------------------------------------------------------
! READ VARIABLES FROM RESTART FILES               
!------------------------------------------------------------------
ncid = nc_open_file_readonly(original_cice_restart_file, source)
call nc_get_variable_size(ncid, 'aicen', nsize, source)
Nx = nsize(1)
Ny = nsize(2)
Ncat = nsize(3)
call nc_close_file(ncid, source)

! Read the pre-assim variables
call read_cice_state_variable('aicen', aicen_original, original_cice_restart_file)
call read_cice_state_variable('vicen', vicen_original, original_cice_restart_file)
call read_cice_state_variable('vsnon', vsnon_original, original_cice_restart_file)

! Read the post-assim variables
call read_cice_state_variable('aicen', aicen, dart_to_cice_input_file)
call read_cice_state_variable('vicen', vicen, dart_to_cice_input_file)
call read_cice_state_variable('vsnon', vsnon, dart_to_cice_input_file)
call read_cice_state_variable('Tsfcn', Tsfcn, dart_to_cice_input_file)
call read_cice_state_variable('qice001', qice001, dart_to_cice_input_file)
call read_cice_state_variable('qice002', qice002, dart_to_cice_input_file)
call read_cice_state_variable('qice003', qice003, dart_to_cice_input_file)
call read_cice_state_variable('qice004', qice004, dart_to_cice_input_file)
call read_cice_state_variable('qice005', qice005, dart_to_cice_input_file)
call read_cice_state_variable('qice006', qice006, dart_to_cice_input_file)
call read_cice_state_variable('qice007', qice007, dart_to_cice_input_file)
call read_cice_state_variable('qice008', qice008, dart_to_cice_input_file)
call read_cice_state_variable('sice001', sice001, dart_to_cice_input_file)
call read_cice_state_variable('sice002', sice002, dart_to_cice_input_file)
call read_cice_state_variable('sice003', sice003, dart_to_cice_input_file)
call read_cice_state_variable('sice004', sice004, dart_to_cice_input_file)
call read_cice_state_variable('sice005', sice005, dart_to_cice_input_file)
call read_cice_state_variable('sice006', sice006, dart_to_cice_input_file)
call read_cice_state_variable('sice007', sice007, dart_to_cice_input_file)
call read_cice_state_variable('sice008', sice008, dart_to_cice_input_file)
call read_cice_state_variable('qsno001', qsno001, dart_to_cice_input_file)
call read_cice_state_variable('qsno002', qsno002, dart_to_cice_input_file)
call read_cice_state_variable('qsno003', qsno003, dart_to_cice_input_file)

!------------------------------------------------------------------
! PERFORM POSTPROCESSING 
!------------------------------------------------------------------
if (postprocess == 'cice') then
   call error_handler(E_MSG, source, 'calling cice postprocessing...')
   call cice_rebalancing(qice001, qice002,     &
                         qice003, qice004,     &
                         qice005, qice006,     &
                         qice007, qice008,     &
                         sice001, sice002,     &
                         sice003, sice004,     &
                         sice005, sice006,     &
                         sice007, sice008,     &
                         qsno001, qsno002,     &
                         qsno003, aicen,       &
                         vicen, vsnon,         &
                         aicen_original,       &
                         vicen_original,       &
                         vsnon_original,       &
                         Tsfcn,                &
                         Ncat, Nx, Ny)
   call error_handler(E_MSG, source, 'cice postprocessing function completed...')

else if (postprocess == 'aice') then
   call error_handler(E_MSG, source, 'calling aice postprocessing...')
   call area_simple_squeeze(qice001, qice002,     &
                            qice003, qice004,     &
                            qice005, qice006,     &
                            qice007, qice008,     &
                            sice001, sice002,     &
                            sice003, sice004,     &
                            sice005, sice006,     &
                            sice007, sice008,     &
                            qsno001, qsno002,     &
                            qsno003, aicen,       &
                            vicen, vsnon,         &
                            aicen_original,       &
                            vicen_original,       &
                            vsnon_original,       &
                            Tsfcn,                &
                            Ncat, Nx, Ny)
   call error_handler(E_MSG, source, 'aice postprocessing function completed...')
else if (postprocess == 'vice') then
   call error_handler(E_MSG, source, 'calling vice postprocessing...')
   call volume_simple_squeeze(qice001, qice002,     &
                              qice003, qice004,     &
                              qice005, qice006,     &
                              qice007, qice008,     &
                              sice001, sice002,     &
                              sice003, sice004,     &
                              sice005, sice006,     &
                              sice007, sice008,     &
                              qsno001, qsno002,     &
                              qsno003, aicen,       &
                              vicen, vsnon,         &
                              aicen_original,       &
                              vicen_original,       &
                              vsnon_original,       &
                              Tsfcn,                &
                              Ncat, Nx, Ny)

   call error_handler(E_MSG, source, 'vice postprocessing function completed...')
else
   call error_handler(E_MSG, source, 'No valid postprocessing method called. No adjustments will be made.')
end if

!------------------------------------------------------------------
! WRITE VARIABLES TO RESTART FILE
!------------------------------------------------------------------
ncid = nc_open_file_readwrite(postprocessed_output_file)
call nc_put_variable(ncid, 'aicen', aicen)
call nc_put_variable(ncid, 'vicen', vicen)
call nc_put_variable(ncid, 'vsnon', vsnon)
call nc_put_variable(ncid, 'Tsfcn', Tsfcn)
call nc_put_variable(ncid, 'qice001', qice001)
call nc_put_variable(ncid, 'qice002', qice002)
call nc_put_variable(ncid, 'qice003', qice003)
call nc_put_variable(ncid, 'qice004', qice004)
call nc_put_variable(ncid, 'qice005', qice005)
call nc_put_variable(ncid, 'qice006', qice006)
call nc_put_variable(ncid, 'qice007', qice007)
call nc_put_variable(ncid, 'qice008', qice008)
call nc_put_variable(ncid, 'sice001', sice001)
call nc_put_variable(ncid, 'sice002', sice002)
call nc_put_variable(ncid, 'sice003', sice003)
call nc_put_variable(ncid, 'sice004', sice004)
call nc_put_variable(ncid, 'sice005', sice005)
call nc_put_variable(ncid, 'sice006', sice006)
call nc_put_variable(ncid, 'sice007', sice007)
call nc_put_variable(ncid, 'sice008', sice008)
call nc_put_variable(ncid, 'qsno001', qsno001)
call nc_put_variable(ncid, 'qsno002', qsno002)
call nc_put_variable(ncid, 'qsno003', qsno003)
call nc_close_file(ncid)

!------------------------------------------------------------------
! DEALLOCATE AND FINALIZE                
!------------------------------------------------------------------
deallocate(aicen, vicen, vsnon, Tsfcn, aicen_original, vicen_original, vsnon_original)
deallocate(qice001, qice002, qice003, qice004, qice005, qice006, qice007, qice008)
deallocate(sice001, sice002, sice003, sice004, sice005, sice006, sice007, sice008)
deallocate(qsno001, qsno002, qsno003)

call finalize_utilities('dart_to_cice')
!------------------------------------------------------------------
! END             
!------------------------------------------------------------------
end program dart_to_cice

