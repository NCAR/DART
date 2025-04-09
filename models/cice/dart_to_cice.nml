! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

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
                           file_exist, error_handler, E_ERR, E_MSG, to_upper
use ice_postprocessing_mod, only : cice_rebalancing, area_simple_squeeze, &
                                   volume_simple_squeeze, get_3d_variable, &
                                   write_3d_variable
use  netcdf_utilities_mod,  only : nc_check
use  netcdf 

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------
! SET NAMELIST AND ALLOCATE VARIABLES    
!------------------------------------------------------------------
character(len=256) :: dart_to_cice_input_file = 'post_filter_restart.nc'
character(len=256) :: original_cice_restart_file = 'pre_filter_restart.nc'
character(len=256) :: postprocessed_output_file = 'postprocessed_restart.nc'
character(len=128) :: balance_method = 'simple_squeeze'
character(len=128) :: postprocess = 'cice'

namelist /dart_to_cice_nml/ dart_to_cice_input_file,    &
                            original_cice_restart_file, &
                            postprocessed_output_file,  &
                            balance_method,             &
                            postprocess

! general variable iniatlization
character(len=512) :: string1, string2, msgstring
character(len=128) :: method
character(len=3)   :: nchar

integer :: iunit, io, ncid, dimid, Ncat, Nx, Ny
real(r8), allocatable :: aicen_original(:,:,:), vicen_original(:,:,:), vsnon_original(:,:,:)
real(r8), allocatable :: aicen(:,:,:), vicen(:,:,:), vsnon(:,:,:), Tsfcn(:,:,:)
real(r8), allocatable :: qice001(:,:,:), qice002(:,:,:), qice003(:,:,:), qice004(:,:,:), qice005(:,:,:), qice006(:,:,:), qice007(:,:,:), qice008(:,:,:)
real(r8), allocatable :: sice001(:,:,:), sice002(:,:,:), sice003(:,:,:), sice004(:,:,:), sice005(:,:,:), sice006(:,:,:), sice007(:,:,:), sice008(:,:,:)
real(r8), allocatable :: qsno001(:,:,:), qsno002(:,:,:), qsno003(:,:,:)

!------------------------------------------------------------------
! INIALIZE AND PERFORM CHECKS ON FILES                
!------------------------------------------------------------------
call initialize_utilities(progname='dart_to_cice')

call find_namelist_in_file("input.nml", "dart_to_cice_nml", iunit)
read(iunit, nml = dart_to_cice_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_cice_nml")

method = balance_method
call to_upper(method)

write(string1,*) 'converting DART output file "'// &
                 &trim(dart_to_cice_input_file)//'" to one CICE will like'
write(string2,*) 'using the "'//trim(balance_method)//'" method.'
call error_handler(E_MSG,'dart_to_cice',string1,text2=string2)

if ( .not. file_exist(dart_to_cice_input_file) ) then
   write(string1,*) 'cannot open "', trim(dart_to_cice_input_file),'" for updating.'
   call error_handler(E_ERR,'dart_to_cice:filename not found ',trim(dart_to_cice_input_file))
endif

if ( .not. file_exist(original_cice_restart_file) ) then
   write(string1,*) 'cannot open "', trim(original_cice_restart_file),'" for reading.'
   call error_handler(E_ERR,'dart_to_cice:filename not found ',trim(original_cice_restart_file))
endif

!------------------------------------------------------------------
! READ VARIABLES FROM RESTART FILES               
!------------------------------------------------------------------
! Read the pre-assim variables
call nc_check( nf90_open(trim(original_cice_restart_file), NF90_NOWRITE, ncid), &
                  'dart_to_cice', 'open "'//trim(original_cice_restart_file)//'"')

! get dimension information
call nc_check(nf90_inq_dimid(ncid, "ncat", dimid), &
              'dart_to_cice', 'inquire ncat dimid from "'//trim(original_cice_restart_file)//'"')
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Ncat), &
              'dart_to_cice', 'inquire ncat from "'//trim(original_cice_restart_file)//'"')
call nc_check(nf90_inq_dimid(ncid,"ni",dimid), &
               'dart_to_cice', 'inquire ni dimid from "'//trim(original_cice_restart_file)//'"')
call nc_check(nf90_inquire_dimension(ncid,dimid,len=Nx),&
               'dart_to_cice', 'inquire ni from "'//trim(original_cice_restart_file)//'"')
call nc_check(nf90_inq_dimid(ncid,"nj",dimid), &
               'dart_to_cice', 'inquire nj dimid from "'//trim(original_cice_restart_file)//'"')
call nc_check(nf90_inquire_dimension(ncid,dimid,len=Ny),&
               'dart_to_cice', 'inquire nj from "'//trim(original_cice_restart_file)//'"')

call get_3d_variable(ncid, 'aicen', aicen_original, original_cice_restart_file)
call get_3d_variable(ncid, 'vicen', vicen_original, original_cice_restart_file)
call get_3d_variable(ncid, 'vsnon', vsnon_original, original_cice_restart_file)

call nc_check(nf90_close(ncid),'dart_to_cice', 'close '//trim(original_cice_restart_file))

! Read the post-assim variables 
call nc_check( nf90_open(trim(dart_to_cice_input_file), NF90_WRITE, ncid), &
                  'dart_to_cice', 'open "'//trim(dart_to_cice_input_file)//'"')

! get the key restart variables post-assimilation (allocated in routine)
call get_3d_variable(ncid, 'aicen', aicen, dart_to_cice_input_file)
call get_3d_variable(ncid, 'vicen', vicen, dart_to_cice_input_file)
call get_3d_variable(ncid, 'vsnon', vsnon, dart_to_cice_input_file)
call get_3d_variable(ncid, 'Tsfcn', Tsfcn, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice001', sice001, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice002', sice002, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice003', sice003, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice004', sice004, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice005', sice005, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice006', sice006, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice007', sice007, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice008', sice008, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice001', qice001, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice002', qice002, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice003', qice003, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice004', qice004, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice005', qice005, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice006', qice006, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice007', qice007, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice008', qice008, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qsno001', qsno001, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qsno002', qsno002, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qsno003', qsno003, dart_to_cice_input_file)

call nc_check(nf90_close(ncid),'dart_to_cice', 'close '//trim(dart_to_cice_input_file))

! write(*,*) 'aicen dimensions are ', shape(aicen)
!------------------------------------------------------------------
! PERFORM POSTPROCESSING 
!------------------------------------------------------------------
if (postprocess == 'cice') then
   write(*,*) 'calling cice postprocessing...'
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
   write(*,*) 'cice postprocessing function completed...'

else if (postprocess == 'aice') then 
   write(*,*) 'calling aice postprocessing...'
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
   write(*,*) 'aice postprocessing function completed...'
else if (postprocess == 'vice') then
   write(*,*) 'calling vice postprocessing...'
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
   write(*,*) 'vice postprocessing function completed...'
else
   write(*,*) 'No valid postprocessing method called. No adjustments will be made.'
end if

!------------------------------------------------------------------
! WRITE VARIABLES TO RESTART FILE
!------------------------------------------------------------------
call nc_check( nf90_open(trim(postprocessed_output_file), NF90_WRITE, ncid), &
                  'dart_to_cice', 'open "'//trim(postprocessed_output_file)//'"')

call write_3d_variable(ncid, 'aicen', aicen, postprocessed_output_file)
call write_3d_variable(ncid, 'vicen', vicen, postprocessed_output_file)
call write_3d_variable(ncid, 'vsnon', vsnon, postprocessed_output_file)
call write_3d_variable(ncid, 'Tsfcn', Tsfcn, postprocessed_output_file)
call write_3d_variable(ncid, 'qice001', qice001, postprocessed_output_file)
call write_3d_variable(ncid, 'qice002', qice002, postprocessed_output_file)
call write_3d_variable(ncid, 'qice003', qice003, postprocessed_output_file)
call write_3d_variable(ncid, 'qice004', qice004, postprocessed_output_file)
call write_3d_variable(ncid, 'qice005', qice005, postprocessed_output_file)     
call write_3d_variable(ncid, 'qice006', qice006, postprocessed_output_file)
call write_3d_variable(ncid, 'qice007', qice007, postprocessed_output_file)
call write_3d_variable(ncid, 'qice008', qice008, postprocessed_output_file)
call write_3d_variable(ncid, 'sice001', sice001, postprocessed_output_file)
call write_3d_variable(ncid, 'sice002', sice002, postprocessed_output_file)
call write_3d_variable(ncid, 'sice003', sice003, postprocessed_output_file)
call write_3d_variable(ncid, 'sice004', sice004, postprocessed_output_file)
call write_3d_variable(ncid, 'sice005', sice005, postprocessed_output_file)
call write_3d_variable(ncid, 'sice006', sice006, postprocessed_output_file)     
call write_3d_variable(ncid, 'sice007', sice007, postprocessed_output_file)
call write_3d_variable(ncid, 'sice008', sice008, postprocessed_output_file)
call write_3d_variable(ncid, 'qsno001', qsno001, postprocessed_output_file)
call write_3d_variable(ncid, 'qsno002', qsno002, postprocessed_output_file)
call write_3d_variable(ncid, 'qsno003', qsno003, postprocessed_output_file)

call nc_check(nf90_close(ncid),'dart_to_cice', 'close '//trim(postprocessed_output_file))


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

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
