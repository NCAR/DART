! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: radiance_bias_mod.f90 11427 2017-03-31 21:24:30Z nancy@ucar.edu $

!>  Module to perform variational bias correction
module radiance_bias_mod

use      types_mod,       only : r8, i8, digits12, PI, missing_r8, obstypelength
use  utilities_mod,       only : file_exist, get_unit, check_namelist_read, do_output,    &
                                 find_namelist_in_file, register_module, error_handler,   &
                                 E_ERR, E_MSG, nmlfileunit, do_nml_file, do_nml_term,     &
                                 open_file, close_file, timestamp, to_upper

use obs_sequence_mod,     only : obs_sequence_type, obs_type, get_obs_from_key, get_obs_def, get_obs_values
   
use          obs_def_mod, only : obs_def_type, get_obs_def_location, get_obs_def_time, &
                                 get_obs_def_error_variance, get_obs_def_type_of_obs,  &
                                 get_obs_def_biaspreds

use         obs_kind_mod, only : get_num_types_of_obs, get_index_for_type_of_obs,           &
                                 get_quantity_for_type_of_obs, assimilate_this_type_of_obs, &
                                 QTY_BRIGHTNESS_TEMPERATURE, get_name_for_type_of_obs

use ensemble_manager_mod, only : ensemble_type, get_my_num_vars, get_my_vars,       & 
                                 get_var_owner_index,  map_pe_to_task, all_vars_to_all_copies, &
                                 compute_copy_mean, all_copies_to_all_vars, get_copy, map_task_to_pe

use mpi_utilities_mod,    only : my_task_id, sum_across_tasks, task_count, start_mpi_timer, &
                                 read_mpi_timer, array_broadcast, task_sync, send_sum_to_all

use quality_control_mod, only : good_dart_qc

use radinfo, only: adp_anglebc, emiss_bc, use_edges, newpc4pred, angord, & 
                    jpch_rad, nusis, nuchan, npred, predx, inew_rad, varA, ostats, &
                    gsi_variable_init, init_rad, init_rad_vars, radinfo_write, radinfo_read
!use mpi ! CSS ... the GSI EnSRF has an explicit mpi_allreduce, DART didn't have a proper wrapper
         !             initially, but send_sum_to_all was added in mpi_utilities_mod to fix that.


implicit none


private

! public subroutines...call order from assim_tools_mod:
! 1) radiance_bias_init
! 2) update_biascorr
! 3) if numiter > 1, then iterate between calling
!      apply_biascorr and update_biascorr
! 4) radiance_bias_finalize
public :: radiance_bias_init, radiance_bias_finalize
public :: apply_biascorr, update_biascorr

! public variables
public :: numiter, lupd_satbiasc, nobs_sat

! Below variables visible to subroutines in this module only
logical :: module_initialized = .false. ! Indicates if module initialization subroutine has been called yet
integer :: my_num_radiance_obs = 0  ! Number of radiance obs on this processor (unique to each PE)
integer :: nobs_sat = 0 ! Total number of radiance obs across all processors (all PEs have this)
integer, allocatable   :: my_radiance_ob_keys(:), my_radiance_ob_indices(:), radiance_ob_keys(:)
integer, allocatable   :: numobspersat(:), indxsat(:)
real(r8), allocatable  :: oberrvarmean(:), radiance_ob_errors(:), obfit_post(:)
real(r8), allocatable  :: deltapredx(:,:), biaspreds(:,:)
character(len=obstypelength), allocatable :: base_obs_names(:)
character(len = 255) :: msgstring, msgstring2, msgstring3
type(obs_type)       :: observation
type(obs_def_type)   :: obs_def

! Namelist variables defined in this module
integer(i8)        :: numiter = 1
logical            :: lupd_satbiasc = .false.
character(len=100) :: bias_correction_source = ''
real(r8)           :: biasvar = 0.1_r8

!npred, adp_anglebc, emiss_bc, use_edges, newpc4pred, angord defined and given defaults in "init_rad" from radinfo
namelist / radiance_bias_nml /  numiter, lupd_satbiasc, bias_correction_source, biasvar, &
                                 npred,adp_anglebc,angord,use_edges,emiss_bc,newpc4pred

!============================================================================
! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://svn-dares-dart.cgd.ucar.edu/DART/releases/Manhattan/assimilation_code/modules/assimilation/radiance_bias_mod.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 11427 $"
character(len=128), parameter :: revdate  = "$Date: 2017-03-31 15:24:30 -0600 (Fri, 31 Mar 2017) $"

contains

!============================================================================

! initialialization routines 
subroutine radiance_bias_init(obs_seq,obs_ens_handle,keys,OBS_GLOBAL_QC_COPY)
   type(obs_sequence_type), intent(in)  :: obs_seq
   type(ensemble_type),     intent(in)  :: obs_ens_handle
   integer,                 intent(in)  :: keys(:)
   integer,                 intent(in)  :: OBS_GLOBAL_QC_COPY

   integer :: iunit, io

   call register_module(source, revision, revdate)

   ! do this up front
   module_initialized = .true.

   ! Call these before reading namelist in case we want to overwrite defaults.
   ! The two routines basically set a bunch of constants needed to make
   !   routines in radinfo work and some default namelist variables
   call gsi_variable_init
   call init_rad

   ! Read the namelist. This will overwrite some defaults set in radinfo
   call find_namelist_in_file("input.nml", "radiance_bias_nml", iunit)
   read(iunit, nml = radiance_bias_nml, iostat = io)
   call check_namelist_read(iunit, io, "radiance_bias_nml")

   ! Write the namelist values to the log file
   if (do_nml_file()) write(nmlfileunit, nml=radiance_bias_nml)
   if (do_nml_term()) write(     *     , nml=radiance_bias_nml)

   ! Update some variables declared in radinfo (angord, npred) based on namelist we just read
   call init_rad_vars

   if ( lupd_satbiasc ) then

      numiter = max(1,numiter) ! force numiter >= 1 

      ! Get the number, keys, and indices of the radiance obs ON THIS PROCESSOR
      !   Sets my_num_radiance_obs, my_radiance_ob_keys, and my_radiance_ob_indices, which are accessible to all routines in this module
      call get_my_num_sat_obs(obs_seq,obs_ens_handle,keys,OBS_GLOBAL_QC_COPY)

      ! Now get the total number of radiance obs over all processors; set nobs_sat, which is a public variable
      call sum_across_tasks(my_num_radiance_obs,nobs_sat)

      if (my_task_id() == 0 ) then
         write(msgstring, '(A,1x,I8,1x,A)') 'There are ', nobs_sat, ' radiance obs with QC indicating they will be assimilated'
         call error_handler(E_MSG,'radiance_bias_init:',msgstring)
      endif

      ! Read bias correction files (allocate and fill predx) and allocate deltapredx.
      !  predx is background bias correction coefficients. deltapredx is the increment of bias correction coefficient we ultimately want to find.
      !  We still want to output a bias correction file if lupd_satbiasc == true and nobs_sat == 0; then deltapredx = 0 and output file == input file
      call read_bias_correction

      ! Once nobs_sat has been determined, allocate some variables that are visible by all routines in this module
      if ( nobs_sat > 0 ) then
         allocate(radiance_ob_keys(nobs_sat)) ! currently not used
         allocate(radiance_ob_errors(nobs_sat))
         allocate(base_obs_names(nobs_sat)) 
         allocate(obfit_post(nobs_sat))
         allocate(biaspreds(npred+1, nobs_sat)) ! predictors for bias correction
         biaspreds = 0.0_r8 ! initialize bias predictor array to zero
         allocate(indxsat(nobs_sat)) ! really a GSI-specific variable
         indxsat = 0 ! default of 0. if it's 0 later, that means a particular radiance observation did not have a match in the satinfo file (shouldn't happen)
         
         ! print output
         write(msgstring, '(A,1X,I2,1X,A,1X,I2,1X,A)') 'Updating bias correction with ',numiter, 'iterations and ',npred,' predictors'
         call error_handler(E_MSG,'radiance_bias_init:', msgstring)
      
         ! Check source of bias correction. For now, only GSI okay, but can add new options
         if ( trim(adjustl(bias_correction_source)) == 'gsi') then
            write(msgstring,  '(A)') 'Bias correction source is '//trim(adjustl(bias_correction_source))
            call error_handler(E_MSG,'radiance_bias_init:', msgstring)
         else
            call error_handler(E_ERR, 'radiance_bias_init:', 'Illegal bias_correction_source', &
               source, revision, revdate, text2='bias correction source is: '//trim(bias_correction_source))
         endif
      else
         write(msgstring, '(A)') 'No radiance obs so force numiter = 1 and do not update bias correction coefficients but output a coefficient file.'
         call error_handler(E_MSG,'radiance_bias_init:', msgstring)
         numiter = 1
      endif

   else
      write(msgstring, '(A)') 'No radiance bias correction update so force numiter = 1'
      call error_handler(E_MSG,'radiance_bias_init:', msgstring)
      numiter = 1
   endif
end subroutine radiance_bias_init

subroutine get_my_num_sat_obs(obs_seq,obs_ens_handle,keys,OBS_GLOBAL_QC_COPY)
   type(obs_sequence_type),  intent(in) :: obs_seq
   type(ensemble_type),      intent(in) :: obs_ens_handle
   integer,                  intent(in) :: keys(:)
   integer,                  intent(in) :: OBS_GLOBAL_QC_COPY

   integer  :: i, nn, base_obs_kind, base_obs_type, this_obs_key
   integer  :: radiance_keys_temp(obs_ens_handle%my_num_vars)
   integer  :: radiance_indices_temp(obs_ens_handle%my_num_vars)
   real(r8) :: obs_qc

   ! This routine returns the number, keys, and indices of radiance obs on this processor
   nn = 0
   do i = 1,obs_ens_handle%my_num_vars

      ! Only value of 0 for DART QC field should be assimilated
      obs_qc = obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, i)
      if(nint(obs_qc) /=0) cycle

      this_obs_key = keys(obs_ens_handle%my_vars(i))
      call get_obs_from_key(obs_seq, this_obs_key, observation) 
      call get_obs_def(observation, obs_def)
      base_obs_type = get_obs_def_type_of_obs(obs_def)
      if (base_obs_type > 0) then
         base_obs_kind = get_quantity_for_type_of_obs(base_obs_type)
         if (base_obs_kind == QTY_BRIGHTNESS_TEMPERATURE) then
            nn = nn + 1
            radiance_keys_temp(nn) = this_obs_key
            radiance_indices_temp(nn) = i
         endif
      endif
   enddo

   ! set the relevant variables visible to all routines in this module
   allocate( my_radiance_ob_keys(nn))
   allocate( my_radiance_ob_indices(nn))
   my_radiance_ob_keys(:) = radiance_keys_temp(1:nn) ! keys of radiance obs on this PE; keys is in relation to the global obs_sequence
   my_radiance_ob_indices(:) = radiance_indices_temp(1:nn) ! indices of radiance obs on this PE
   my_num_radiance_obs = nn ! number of radiance obs on this PE

end subroutine get_my_num_sat_obs

!!! generic subroutines below. in theory, these will call different routines based on $bias_correction_source !!

subroutine read_bias_correction
   if ( trim(adjustl(bias_correction_source)) == 'gsi') then
      call radinfo_read ! from radinfo. reads bias correction coefficients, sets many variables, including predx, the background bias correction coefficients
      allocate(deltapredx(npred,jpch_rad)) ! can't be allocated until after calling radinfo_read, which sets jpch_rad to something appropriate
   endif
   deltapredx = 0.0_r8 ! initialize bias correction increment to 0
   call error_handler(E_MSG,'read_bias_correction: ','Done reading bias correction file')
end subroutine read_bias_correction

subroutine write_bias_correction
   if ( trim(adjustl(bias_correction_source)) == 'gsi') then
      call radinfo_write ! from radinfo, also deallocates variables set in radinfo_read
   endif
   call error_handler(E_MSG,'write_bias_correction: ','Done writing bias correction file')
end subroutine write_bias_correction

subroutine apply_biascorr(obs_seq,obs_ens_handle,ens_size)
   type(obs_sequence_type),  intent(in)    :: obs_seq
   type(ensemble_type),      intent(inout) :: obs_ens_handle
   integer,                  intent(in)    :: ens_size
   ! Update the radiance priors based on the latest bias correction coefficients that have 
   !  been updated in update_biascorr
   if ( trim(adjustl(bias_correction_source)) == 'gsi') then
      call apply_biascorr_gsi(obs_seq,obs_ens_handle,ens_size)
   endif
   call error_handler(E_MSG, 'apply_biascorr:', 'Done with apply_biascorr') 
end subroutine apply_biascorr

subroutine update_biascorr(niter,obs_seq, obs_ens_handle, keys, obs_val_index, OBS_MEAN_COPY, ens_size, OBS_GLOBAL_QC_COPY)
   integer(i8), intent(in)               :: niter
   type(obs_sequence_type),  intent(in)  :: obs_seq
   type(ensemble_type), intent(inout)    :: obs_ens_handle
   integer, intent(in)                   :: keys(:)
   integer, intent(in)                   :: obs_val_index, OBS_MEAN_COPY, ens_size, OBS_GLOBAL_QC_COPY

   ! Get information about ALL radiance obs onto ALL processors for bias correction update, which is a global operation
   ! In other words, we need to gather all radiance obs on each processor, which get_global_radiance_info does
   !      get_global_radiance_info also fills some important variables local to this module
   call get_global_radiance_info(obs_seq, obs_ens_handle, keys, obs_val_index, OBS_MEAN_COPY, ens_size, OBS_GLOBAL_QC_COPY)

   ! Now that we have the global information, update the bias correction coefficients
   if ( trim(adjustl(bias_correction_source)) == 'gsi') then
      call update_biascorr_gsi(niter)
   endif
   call error_handler(E_MSG, 'update_biascorr:',  'Done with update_biascorr') 
end subroutine update_biascorr

subroutine get_global_radiance_info(obs_seq, obs_ens_handle, keys, obs_val_index, OBS_MEAN_COPY, ens_size, OBS_GLOBAL_QC_COPY)
   ! This subroutine puts information about radiance obs onto ALL processors for bias correction update, which is a global operation
   type(obs_sequence_type),  intent(in)  :: obs_seq
   type(ensemble_type), intent(inout)    :: obs_ens_handle
   integer, intent(in)                   :: keys(:)
   integer, intent(in)                   :: obs_val_index, OBS_MEAN_COPY, ens_size, OBS_GLOBAL_QC_COPY

   integer :: i, nn, base_obs_kind, base_obs_type
   real(digits12) :: base_rad, elapsed
   real(r8) :: obs_value(1)
   real(r8), allocatable :: mean_temp(:), qc_temp(:)

   ! Populate mean posterior observation for ALL obs (including non-radiances) on this processor
   call compute_copy_mean(obs_ens_handle, 1, ens_size, OBS_MEAN_COPY) ! Mean of the members sum(H(x_i))/N

   call all_copies_to_all_vars(obs_ens_handle) ! get all observations onto each processor for a subset of copies

   ! Get the mean copy on process 0 and then broadcast to all processors
   allocate(mean_temp(obs_ens_handle%num_vars))
   call get_copy(map_task_to_pe(obs_ens_handle, 0), obs_ens_handle, OBS_MEAN_COPY, mean_temp)
   call array_broadcast(mean_temp, 0)  ! now all processors have mean posterior for all obs (including non-radiances)

!  do i = 1, obs_ens_handle%num_vars
     !call replace_obs_values(seq, keys(i), mean_temp(i), ens_mean_index)
!   end do

   ! Get the QC on process 0 and then broadcast to all processors
   allocate(qc_temp(obs_ens_handle%num_vars))
   call get_copy(map_task_to_pe(obs_ens_handle, 0), obs_ens_handle, OBS_GLOBAL_QC_COPY, qc_temp)
   call array_broadcast(qc_temp, 0) ! now all processors have QC for all obs (including non-radiances)

   call start_mpi_timer(base_rad)

   nn = 0
   do i = 1, obs_ens_handle%num_vars ! All processors loop over all variables
      ! we only want radiance obs that were assimilated to be counted
      if(nint(qc_temp(i)) /=0) cycle 

      ! Every pe has information about the global obs sequence
      call get_obs_from_key(obs_seq, keys(i), observation)
      call get_obs_def(observation, obs_def)
      base_obs_type = get_obs_def_type_of_obs(obs_def)
      if (base_obs_type > 0) then
         base_obs_kind = get_quantity_for_type_of_obs(base_obs_type)
         if (base_obs_kind == QTY_BRIGHTNESS_TEMPERATURE) then
            nn = nn + 1
            radiance_ob_keys(nn) = keys(i) ! currently not used
            radiance_ob_errors(nn) = get_obs_def_error_variance(obs_def)
            base_obs_names(nn) = get_name_for_type_of_obs(base_obs_type)
            biaspreds(1:npred,nn) = get_obs_def_biaspreds(obs_def,npred,int(1,8),npred)
            call get_obs_values(observation, obs_value(1:1), obs_val_index)
            obfit_post(nn) = obs_value(1) - mean_temp(i) ! y - H(x_bar)
           !obfit_post(nn) = obs_value(1) - mean_temp(keys(i)) ! y - H(x_bar)
            if (my_task_id() == 0 ) write(*,*)'CSS ob_val,Hx_a = ',obs_value(1),mean_temp(keys(i)),mean_temp(i)
         endif
      endif
   enddo

   elapsed = read_mpi_timer(base_rad)
   if (my_task_id() == 0 ) then
      write(msgstring, '(A29,1x,F12.4,1x,A,1x,I5)') 'time to query radiance obs = ', elapsed, 'rank ', my_task_id()
      call error_handler(E_MSG,'get_radiance_info_global:',msgstring)
   endif
   
   ! Sanity check
   if ( nn /= nobs_sat ) then
      write(msgstring,  '(A)') 'Mismatching numbers'
      write(msgstring2, '(A,I5)') 'nobs_sat is', nobs_sat
      write(msgstring3, '(A,I5)') 'nn is ',nn
      call error_handler(E_ERR, 'get_radiance_info_global:', msgstring, &
         source, revision, revdate, text2=msgstring2, text3=msgstring3)
   endif

   ! Reset to the original input not mess-up other things later
   call all_vars_to_all_copies(obs_ens_handle) ! get all COPIES on each processor for a subset of observations

   ! Clean-up
   deallocate(mean_temp,qc_temp)

end subroutine get_global_radiance_info

!! Routines for GSI-based bias correction !! 
subroutine map_to_gsi_satinfo
   character(len=256)   :: str_channel, full_ob_name
   integer :: i, j, nproc, numproc, i1, i2, ierr
   real(digits12) :: base_rad, elapsed
   integer :: indxsat1(nobs_sat)

   ! figure out which line in GSI "satinfo" each radiance ob belongs to. Output is in indxsat
   call start_mpi_timer(base_rad)
  
   nproc = my_task_id()
   numproc = task_count()
   i1 = nproc*real(nobs_sat/numproc) + 1
   i2 = (nproc+1)*real(nobs_sat/numproc)
   if (nproc == numproc-1) i2 = nobs_sat

   indxsat1 = 0
!  do i = 1, nobs_sat ! loop over ALL radiance obs
   do i = i1, i2 ! loop over a subset of observations on this PE ... indxsat, indxsat1 both initialized to 0 and on all PEs, so we can do this 
                 !    parallelization hack. base_obs_names also available on all processors
      do j = 1,jpch_rad
         ! Make a string containing satellite,sensor,and channel to match obs_def_radiance_mod.f90 (the "obs_types") entries
         write(str_channel,fmt='(i5)') nuchan(j)
         full_ob_name = trim(adjustl(nusis(j)))//'_ch'//trim(adjustl(str_channel)) ! e.g., amsua_n19_ch7, amsua_metop-a_ch9
         call to_upper(full_ob_name) ! Make uppercase to match the case in obs_def_radiance_mod.f90
         call replace_hyphen(full_ob_name)  ! Bad things will happen in DART preprocess "obs_def" files if there are hyphens in names. Replace with underscores
         if ( trim(adjustl(full_ob_name)) == trim(adjustl(base_obs_names(i))) ) then ! we found a match for this satellite/sensor/channel.
            indxsat1(i) = j 
            exit ! all done. break out of loop and go to next ob
         endif
      enddo
      ! make sure indxsat /= 0 ... this really shouldn't happen
      if ( indxsat1(i) == 0 ) then
         write(msgstring, '(A)') 'indxsat1 == 0 for '//trim(adjustl(base_obs_names(i)))//' which should not occur'
         call error_handler(E_ERR, 'map_to_gsi_satinfo:', msgstring, source, revision, revdate)
      endif
   enddo

   ! now get indxsat1 across all processors and store in indxsat
   call send_sum_to_all(indxsat1,indxsat)
!  call mpi_allreduce(indxsat1,indxsat,int(nobs_sat),mpi_integer,mpi_sum,mpi_comm_world,ierr)

   elapsed = read_mpi_timer(base_rad)
   if (my_task_id() == 0 ) then
      write(msgstring, '(A,1x,F12.4,1x,A,1x,I5)') 'time to map to GSI satinfo : ', elapsed, 'rank ', my_task_id()
      call error_handler(E_MSG,'map_to_gsi_satinfo:',msgstring)
   endif
end subroutine map_to_gsi_satinfo

subroutine apply_biascorr_gsi(obs_seq,obs_ens_handle,ens_size)
   type(obs_sequence_type),  intent(in)    :: obs_seq
   type(ensemble_type),      intent(inout) :: obs_ens_handle
   integer,                  intent(in)    :: ens_size

   integer              :: i, j, np, base_obs_kind, base_obs_type, ii
   character(len=256)   :: str_channel, full_ob_name
   character(len=obstypelength) :: base_obs_name
   real(r8) :: these_biaspreds(npred+2), bias_amount, before(ens_size), after(ens_size), intermed(ens_size)

   ! Module must be initialized before calling this routine
   if ( .not. module_initialized ) call error_handler(E_ERR, 'apply_bias_corr_gsi:', 'radiance_bias_init not initialized. Very bad.', source, revision, revdate)

   ! my_num_radiance_obs is the number of radiance obs with acceptable QC on this processor defined in get_my_num_sat_obs
   do i = 1, my_num_radiance_obs ! could also be size(my_radiance_ob_keys)
      ii = my_radiance_ob_indices(i) ! my_radiance_ob_indices allocated/filled in get_my_num_sat_obs
      call get_obs_from_key(obs_seq, my_radiance_ob_keys(i), observation) ! my_radiance_ob_keys allocated/filled in get_my_num_sat_obs
      call get_obs_def(observation, obs_def)
      base_obs_type = get_obs_def_type_of_obs(obs_def)
      before(1:ens_size) = obs_ens_handle%copies(1:ens_size,ii)
      if (base_obs_type > 0) then
         base_obs_kind = get_quantity_for_type_of_obs(base_obs_type)
         if (base_obs_kind == QTY_BRIGHTNESS_TEMPERATURE) then
            base_obs_name = get_name_for_type_of_obs(base_obs_type)
            do j = 1,jpch_rad
               ! Make a string containing satellite,sensor,and channel to match obs_def_radiance_mod.f90 entries
               write(str_channel,fmt='(i5)') nuchan(j)
               full_ob_name = trim(adjustl(nusis(j)))//'_ch'//trim(adjustl(str_channel)) ! e.g., amsua_n19_ch7, amsua_metop-a_ch9
               call to_upper(full_ob_name) ! Make uppercase to match the case in obs_def_radiance_mod.f90
               call replace_hyphen(full_ob_name)  ! Bad things will happen in DART preprocess "obs_def" files if there are hyphens in names. Replace with underscores
               if ( trim(adjustl(full_ob_name)) == trim(adjustl(base_obs_name)) ) then ! we found a match for this satellite/sensor/channel.
                  these_biaspreds(1:npred+2) = get_obs_def_biaspreds(obs_def,npred+2,int(1,8),npred+2)
                  ! subtract bias_amount to start with unbiascorrected values before adding the updated correction
                  ! Note that obs_ens_handle%copies has all obs on this PE, not just radiances. Fortunately, we saved the indices we need to update (my_radiance_ob_indices)
                  bias_amount = these_biaspreds(npred+2)
                  if (.not. adp_anglebc) then
                     obs_ens_handle%copies(1:ens_size,ii) = obs_ens_handle%copies(1:ens_size,ii) - bias_amount + these_biaspreds(npred+1)  ! incorporate total angle-dependent BC
                  else
                     obs_ens_handle%copies(1:ens_size,ii) = obs_ens_handle%copies(1:ens_size,ii) - bias_amount
                  endif
                  intermed(1:ens_size) = obs_ens_handle%copies(1:ens_size,ii)
                  do np=1,npred
                     obs_ens_handle%copies(1:ens_size,ii) = obs_ens_handle%copies(1:ens_size,ii) + &
                                                             these_biaspreds(np)*(predx(np,j)+deltapredx(np,j))
                  enddo
                  after(1:ens_size) = obs_ens_handle%copies(1:ens_size,ii)
                  write(*,fmt='(a,5f8.2)')' CSS BC1/2 = ',bias_amount,sum( these_biaspreds(1:npred)*(predx(1:npred,j)+deltapredx(1:npred,j)) ),maxval(before),maxval(intermed),maxval(after)
                 !write(*,fmt='(a)')' CSS name = ',trim(adjustl(base_obs_name))
                  exit ! all done for this ob. break out of loop and go to next ob
               endif
            enddo
         else
            write(msgstring,  '(A)') 'Not a radiance ob QTY. Very bad.'
            call error_handler(E_ERR, 'apply_bias_corr_gsi:', msgstring, &
               source, revision, revdate)
         endif
      else
         write(msgstring,  '(A)') 'Not a radiance ob TYPE. Very bad.'
         call error_handler(E_ERR, 'apply_bias_corr_gsi:', msgstring, &
            source, revision, revdate)
      endif
   enddo

   call task_sync ! Call MPI barrier. Don't move on until all processors are done updating their respective obs (mostly here for debugging)

end subroutine apply_biascorr_gsi

subroutine update_biascorr_gsi(niter)

   ! Code taken directly from GSI/enkf/radbias.f90 with minor modifications
   use types_mod, only : i_kind =>i8, r_kind =>r8, r_double =>digits12  ! CSS added

   integer :: nproc, numproc, mpi_realkind ! CSS added this line and next 4
   real(r_kind) :: two = 2.0_r_kind
   real(r_kind) :: r10 = 10.0_r_kind
   real(r_kind) :: zero = 0.0_r_kind
   real(r_double) :: elapsed ! CSS replaced mpi_wtime with DART mpi calls below

   integer(i_kind) i,m,i1,i2,nn,n
   real(r_kind) increment(npred),biaserrvar,a(npred,npred),atmp(npred,npred)
   real(r_kind) inctmp(npred)
   real(r_kind) buffertmp(npred,jpch_rad)
   real(r_kind), allocatable, dimension(:,:) :: biaspredtmp
   real(r_kind), allocatable, dimension(:) :: obinc
   real(r_kind) deltapredx1(npred,jpch_rad)
   real(r_double) t1
   integer(i_kind), intent(in) :: niter
   !integer(i_kind) ierr
   integer :: ierr ! CSS removed "i_kind" designation
   character(len=72) fmt
   write(fmt, '("(i2,1x,i4,1x,a20,1x,i4,",I0,"(1x,e10.3))")') npred

   ! Begin CSS added to make below code work
   nproc = my_task_id()
   numproc = task_count()
  !if (r_kind == r_double) then
!  if (r8 == digits12) then ! No need to define mpi_realkind if using send_sum_to_all rather than explicit mpi_allreduce
!     mpi_realkind = mpi_real8
!  else ! if (r_kind == r_single) then
!     mpi_realkind = mpi_real4
!  endif

   ! these routines only need to be called once
   if ( niter .eq. 1 ) then
      call map_to_gsi_satinfo ! fill indxsat, which, for the ith radiance ob, is the corresponding line number in satinfo
      call channelstats       ! using indxsat, fill oberrvarmean, numobspersat, (mean oberr, num obs for each channel)
   endif
   ! End CSS

   ! below code taken directly from GSI/enkf/radbias.f90 with very minor modifications
! if (nproc == 0) t1 = mpi_wtime() ! do this loop in parallel (a chunk of channels/sensors on each processor).
  if (nproc == 0) call start_mpi_timer(t1) ! CSS do this loop in parallel (a chunk of channels/sensors on each processor).
  i1 = nproc*real(jpch_rad/numproc) + 1
  i2 = (nproc+1)*real(jpch_rad/numproc)
  if (nproc == numproc-1) i2 = jpch_rad
  deltapredx1=0._r_kind
  do i=i1,i2
   allocate( biaspredtmp(npred,numobspersat(i)) )
   allocate( obinc(numobspersat(i)) )
   if (newpc4pred) then
       ! set variances for bias predictor coeff. based on diagonal info
       ! of previous analysis error variance. This code copied from berror.f90.
       do n=1,npred
          if (inew_rad(i)) then
             biaserrvar = 10000.0_r_kind
          else
             if (numobspersat(i) == 0) then
                ! channel missing, increase to twice previous analysis cycle
                ! add 10**-6 to prevent vanishly small error variance
                varA(n,i)=two*varA(n,i)+1.0e-6_r_kind
                biaserrvar=varA(n,i)
             else
                biaserrvar=1.1_r_kind*varA(n,i)+1.0e-6_r_kind
             end if
             if (biaserrvar > r10) biaserrvar = r10
             if (varA(n,i)>10000.0_r_kind) varA(n,i)=10000.0_r_kind
          endif
       end do
   else
       if (biasvar < 0.) then
          ! if biasvar set < 0 in namelist, background err variance
          ! is inversely proportional to number of obs (with -biasvar
          ! as proportionality constant). Maximum allowed value 0.1.
          if (numobspersat(i) > int(-biasvar/0.1_r_kind)) then
             biaserrvar = -biasvar/numobspersat(i)
          else
             biaserrvar = 0.1_r_kind
          endif
          !if (niter .eq. 2) print *,'biaserrvar:',numobspersat(i),trim(adjustl(nusis(i))),nuchan(i),biaserrvar
       else
          biaserrvar = biasvar ! single constant value
       endif
   endif
   if (oberrvarmean(i) > 1.e10 .or. numobspersat(i) == 0) then
      deallocate(biaspredtmp,obinc)
      cycle
   endif
   ! compute B**-1 + p * R**-1 * pT
   a = 0._r_kind
   do n=1,npred
      nn = 0
      do m=1,nobs_sat
          ! only use the numobspersat(i) obs associated with this channel/instrument
          if (indxsat(m) == i) then
             nn = nn + 1
             biaspredtmp(n,nn) = biaspreds(n,m)/sqrt(radiance_ob_errors(m)) ! CSS modified sqrt(oberrvar(nobs_conv+nobs_oz+m))
          end if
      enddo
      a(n,n) = 1._r_kind/biaserrvar
   enddo
   if (r_kind == kind(1.d0)) then
      call dgemm('N', 'T', npred, npred, nn, 1.d0, &
                 biaspredtmp, npred, biaspredtmp, npred, 0.d0, atmp, npred)
   else
      call sgemm('N', 'T', npred, npred, nn, 1.e0, &
                 biaspredtmp, npred, biaspredtmp, npred, 0.e0, atmp, npred)
   endif
   a = a + atmp
   ! compute inverse of symmetric matrix a = B**-1 + p * R**-1 * pT.
   call symminv(a,npred)
   ! a now contains analysis error covariance matrix (inverse of Hessian).
   ! update bias predictor variance info (varA) assuming a is diagonal.
   if (niter == numiter .and. newpc4pred) then
       do n=1,npred
          varA(n,i) = a(n,n)
       enddo
   endif
   ! p * R**-1
   do n=1,npred
      nn = 0
      do m=1,nobs_sat
          if (indxsat(m) == i) then
             nn = nn + 1
             biaspredtmp(n,nn) = biaspreds(n,m)/radiance_ob_errors(m) ! CSS modified oberrvar(nobs_conv+nobs_oz+m)
          end if
      enddo
   enddo
   nn = 0
   do m=1,nobs_sat
      if (indxsat(m) == i) then
         nn = nn + 1
         obinc(nn) = obfit_post(m) ! CSS modified obfit_post(nobs_conv+nobs_oz+m)
      end if
   enddo
   ! update bias correction coefficients for this channel/sensor.
   ! b = b +  (B**-1 + p * R**-1 * pT)**-1 * (p * R**-1) * (y - Hx)
   if (r_kind == kind(1.d0)) then
      call dgemv('N',npred,numobspersat(i),1.d0,biaspredtmp,npred,obinc,1,0.d0,inctmp,1)
      call dgemv('N',npred,npred,1.d0,a,npred,inctmp,1,0.d0,increment,1)
   else
      call sgemv('N',npred,numobspersat(i),1.e0,biaspredtmp,npred,obinc,1,0.e0,inctmp,1)
      call sgemv('N',npred,npred,1.e0,a,npred,inctmp,1,0.e0,increment,1)
   endif
   ! bias correction increment.
   !deltapredx1(:,i) = increment 
   ! blend of previous iteration and new estimate to improve convergence.
   if (niter .eq. 1) then
      deltapredx1(:,i) = increment
   else
      deltapredx1(:,i) = 0.5*(deltapredx(:,i) + increment)
   end if
   if (index(nusis(i),'amsua') .gt. 0) then
       write(6,fmt) niter,i,trim(adjustl(nusis(i))),nuchan(i),deltapredx1(:,i)
   end if
   deallocate(biaspredtmp)
   deallocate(obinc)
  enddo
! if (nproc == 0)  print *,'time to update bias correction on root',mpi_wtime()-t1,'secs'
  elapsed = read_mpi_timer(t1) ! CSS
  if (nproc == 0)  print *,'time to update bias correction on root',elapsed,'secs' ! CSS
! t1 = mpi_wtime()
  call start_mpi_timer(t1) ! CSS
  call send_sum_to_all(deltapredx1,deltapredx)
! call mpi_allreduce(deltapredx1,deltapredx,int(npred*jpch_rad),mpi_realkind,mpi_sum,mpi_comm_world,ierr) ! CSS added int
  if (niter == numiter .and. newpc4pred) then
     ! distribute updated varA to all processors.
     buffertmp=zero
     do i=i1,i2
     do n=1,npred
       buffertmp(n,i) = varA(n,i)
     enddo
     enddo
     call send_sum_to_all(buffertmp,varA)
!    call mpi_allreduce(buffertmp,varA,int(jpch_rad*npred),mpi_realkind,mpi_sum,mpi_comm_world,ierr) ! CSS added int
     ! update ostats
     do i=1,jpch_rad
        ostats(i) = numobspersat(i)
     enddo
  endif
! if (nproc == 0) print *,'time in update_biascorr mpi_allreduce on root = ',mpi_wtime()-t1
  elapsed = read_mpi_timer(t1) ! CSS
  if (nproc == 0) print *,'time in update_biascorr mpi_allreduce on root = ',elapsed ! CSS
  !CSS added below for diagnostics
! if (nproc == 0) then
!    do nn=1,nobs_sat
!       write(*,*)'CSS nn/y-Hx,R = ',nn,obfit_post(nn),radiance_ob_errors(nn)
!       write(*,*)'CSS nn/P1,P4,P9 = ',nn,biaspreds(1,nn),biaspreds(4,nn),biaspreds(9,nn)
!     enddo
!  endif
  ! End CSS 
end subroutine update_biascorr_gsi

subroutine channelstats ! slightly modified from GSI/enkf/enkf_obsmod.f90
   use types_mod, only : i_kind =>i8, r_kind =>r8, r_double =>digits12  ! CSS
   integer(i_kind) :: nob, i
   real(r_kind) :: zero = 0.0_r_kind ! CSS added
   ! count number of obs per channel/sensor.
   allocate(numobspersat(jpch_rad))
   allocate(oberrvarmean(jpch_rad))
   numobspersat = 0
   oberrvarmean = zero
   do nob=1,nobs_sat
      i=indxsat(nob) ! line in satinfo for this observation
      numobspersat(i) = numobspersat(i) + 1
      oberrvarmean(i) = oberrvarmean(i) + radiance_ob_errors(nob)
   enddo
   ! average ob error for each channel.
   do i=1,jpch_rad
      if (numobspersat(i) > 0) then
         oberrvarmean(i) = oberrvarmean(i)/real(numobspersat(i),r_kind)
      else
         oberrvarmean(i) = 9.9e31_r8
      end if
   enddo
end subroutine channelstats

! ----------------
! utility routines
! ----------------
subroutine replace_hyphen(string)
   character(len=*), intent(inout) :: string
   integer :: j
   do j = 1,len_trim(string)
      if ( string(j:j) == '-') string(j:j) = '_'
   enddo
end subroutine replace_hyphen

subroutine symminv(a,n)  ! from GSI/src/enkf/radbias.f90
  use types_mod,          only : i_kind =>i8, r_kind =>r8   ! CSS
  ! cholesky decomp inverse of a symm. matrix.
  integer(i_kind), intent(in) :: n  ! CSS made integer "i_kind"
  real(r_kind), intent(inout) :: a(n,n)
  integer ierr,i,j
  if (r_kind == kind(1.d0)) then
     call dpotrf('L',n,a,n,ierr)
     if (ierr /= 0) then
        print *,'ierr=',ierr,'in dpotrf!'
     end if
     call dpotri('L',n,a,n,ierr)
     if (ierr /= 0) then
        print *,'ierr=',ierr,'in dpotri!'
     end if
  else if (r_kind == kind(1.e0)) then
     call spotrf('L',n,a,n,ierr)
     if (ierr /= 0) then
        print *,'ierr=',ierr,'in spotrf!'
     end if
     call spotri('L',n,a,n,ierr)
     if (ierr /= 0) then
        print *,'ierr=',ierr,'in spotri!'
     end if
  else
     print *,'kind specification must be default real or double precision'
     stop
  endif
  do j=2,n
   do i=1,j-1
     a(i,j) = a(j,i)
   enddo
  enddo
end subroutine symminv

! ------------------
! final routine
! ------------------
subroutine radiance_bias_finalize

   ! Bias correction coeffs increment just updated in last call to update_biascorr, but not yet added to predx, so do so
   ! Note that we don't update predx until here; predx = prior + updated increment
   predx = predx + deltapredx
   deltapredx = 0.0 ! probably unneccesary, but reset for safety
   if (my_task_id() == 0 ) call write_bias_correction ! output from root even if nobs_sat == 0

   ! clean-up
   if( allocated(my_radiance_ob_keys))    deallocate(my_radiance_ob_keys)
   if( allocated(my_radiance_ob_indices)) deallocate(my_radiance_ob_indices)
   if( allocated(radiance_ob_keys))       deallocate(radiance_ob_keys)
   if( allocated(numobspersat))           deallocate(numobspersat)
   if( allocated(oberrvarmean))           deallocate(oberrvarmean)
   if( allocated(indxsat))                deallocate(indxsat)
   if( allocated(radiance_ob_errors))     deallocate(radiance_ob_errors)
   if( allocated(deltapredx))             deallocate(deltapredx)
   if (allocated(base_obs_names))         deallocate(base_obs_names)
   if (allocated(obfit_post))             deallocate(obfit_post)
   if (allocated(biaspreds))              deallocate(biaspreds)
end subroutine radiance_bias_finalize

!========================================================================
! end module radiance_bias_mod
!========================================================================

end module radiance_bias_mod

! <next few lines under version control, do not edit>
! $URL: https://svn-dares-dart.cgd.ucar.edu/DART/releases/Manhattan/assimilation_code/modules/assimilation/radiance_bias_mod.f90 $
! $Id: radiance_bias_mod.f90 11427 2017-03-31 21:24:30Z nancy@ucar.edu $
! $Revision: 11427 $
! $Date: 2017-03-31 15:24:30 -0600 (Fri, 31 Mar 2017) $
