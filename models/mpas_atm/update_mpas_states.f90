! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program update_mpas_states

!----------------------------------------------------------------------
! purpose: interface between DART and the model model
!
! method: Read DART analysis vector in netcdf and replace the corresponding
!         field in the mpas file to advance model after running this program.
!         Updated to process all ensemble members.
!
!         The update_mpas_states_nml namelist defines the input and output file
!         name lists for all ensemble members.
!         The input list should be matched with output_state_file_list in &filter_nml.
!
!  variables that are not wind are copied from one file to the another
!  variables that are wind are reconstructed
!  
! author: Soyoung Ha 23 Aug 16
!         Updated in 4 May 2017 for the Manhattan release
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, get_next_filename, E_ERR, E_MSG, error_handler
use time_manager_mod, only : time_type, print_time, print_date, operator(-), &
                             get_time, get_date, operator(/=)
use        model_mod, only : static_init_model, &
                             get_model_size, &
                             get_analysis_time, update_u_from_reconstruct, &
                             use_increments_for_u_update, uv_increments_cell_to_edges, &
                             uv_field_cell_to_edges, dom_id => anl_domid

use state_structure_mod, only : get_num_variables, get_variable_name, &
                                get_variable_size, get_varid_from_varname, &
                                get_dim_lengths

use netcdf_utilities_mod, only : nc_open_file_readonly, &
                                 nc_open_file_readwrite, &
                                 nc_get_variable, nc_put_variable, &
                                 nc_close_file, nc_get_variable_size, &
                                 nc_get_variable_info
                                 
implicit none

character(len=*), parameter :: source   = 'models/mpas_atm/update_mpas_states.f90'

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------
character(len=256)  ::  update_input_file_list = 'filter_out.txt'
character(len=256)  :: update_output_file_list = 'filter_in.txt'

namelist /update_mpas_states_nml/ update_input_file_list, update_output_file_list

!----------------------------------------------------------------------
character (len=256)   :: next_infile, next_outfile
character (len=512)   :: string1
integer               :: iunit, io
integer               :: ncAnlID, ncBckID, istatus
integer               :: filenum, i
real(r8), allocatable :: variable(:)
type(time_type)       :: model_time
type(time_type)       :: state_time
integer               :: file_dims(3) ! File (nVertLevels, nEdges | nCells, Time)
integer               :: file_dimlens(3) 
integer               :: dims(2) ! State (nVertLevels, nEdges | nCells)
integer               :: dimlens(2) 
real(r8), allocatable :: u(:,:), ucell(:,:), vcell(:,:) 
real(r8), allocatable :: ucell_dart(:,:), vcell_dart(:,:)
!----------------------------------------------------------------------

call initialize_utilities(progname=source)

! Read the namelist to get the input filename. 
call find_namelist_in_file("input.nml", "update_mpas_states_nml", iunit)
read(iunit, nml = update_mpas_states_nml, iostat = io)
call check_namelist_read(iunit, io, "update_mpas_states_nml")

call static_init_model()

!----------------------------------------------------------------------
! Reads lists of input mpas (prior) and filter (analysis) files 
! HK @todo this is bad to have a serial loop around files.
!----------------------------------------------------------------------
filenum = 1
fileloop: do        ! until out of files

  ! get a file name from the list, one at a time.
  next_infile  = get_next_filename( update_input_file_list, filenum)
  next_outfile = get_next_filename(update_output_file_list, filenum)
  if (next_infile == '' .or. next_outfile == '') exit fileloop

  ncAnlID = nc_open_file_readonly(next_infile, 'update_mpas_states - open readonly') ! analysis from DART
  ncBckID = nc_open_file_readwrite(next_outfile, 'update_mpas_states - open readwrite') ! background, original mpas file

  model_time = get_analysis_time(ncBckID, next_outfile)
  state_time = get_analysis_time(ncAnlID, next_infile)
  call print_time(state_time,'DART current time')
  call print_time(model_time,'mpas current time')

  if ( model_time /= state_time ) then
     call print_time(state_time,'DART current time',logfileunit)
     call print_time(model_time,'mpas current time',logfileunit)
     write(string1,*) trim(next_infile),' current time must equal model time'
     call error_handler(E_ERR,'update_mpas_states',string1,source)
  endif

  ! copy variables that are not wind from analysis to background
  varloop: do i = 1, get_num_variables(dom_id)

      if (get_variable_name(dom_id,i) == 'uReconstructZonal' .or. &
          get_variable_name(dom_id,i) == 'uReconstructMeridional'.or. &
          get_variable_name(dom_id,i) == 'u') cycle varloop

      allocate(variable(get_variable_size(dom_id, i)))

      call nc_get_variable_size(ncAnlId, get_variable_name(dom_id,i), file_dimlens)
      call nc_get_variable(ncAnlID, get_variable_name(dom_id,i), variable, nc_count=file_dimlens)
      call nc_put_variable(ncBckID, get_variable_name(dom_id,i), variable, nc_count=file_dimlens)

      deallocate(variable)

  enddo varloop

  ! deal with wind
  if (update_u_from_reconstruct) then

     if (use_increments_for_u_update) then
        !  Read in the previous reconstructed winds from the original mpas netcdf file
        !  and compute what increments (changes in values) were added by the assimilation.
        !  Read in the original edge normal wind 'u' field from that same mpas netcdf
        !  file and add the interpolated increments to compute the updated 'u' values.

        ! read in u, uReconstrtuctZonal, uReconstructMeridional from background
        ! read in uReconstrtuctZonal, uReconstructMeridional from analysis

        call nc_get_variable_info(ncBckID, 'u', dimlens=file_dimlens) ! not in state structure
        allocate(u(file_dimlens(1), file_dimlens(2)))
        call nc_get_variable(ncBckID, 'u', u)

        dims = get_dim_lengths(dom_id, get_varid_from_varname(dom_id, 'uReconstructZonal'))
        allocate(ucell(dims(1), dims(2)), ucell_dart(dims(1), dims(2)))
        call nc_get_variable(ncBckID, 'uReconstructZonal', ucell)
        call nc_get_variable(ncAnlID, 'uReconstructZonal', ucell_dart)
        ucell = ucell_dart - ucell ! u increments

        dims = get_dim_lengths(dom_id, get_varid_from_varname(dom_id, 'uReconstructMeridional'))
        allocate(vcell(dims(1), dims(2)), vcell_dart(dims(1), dims(2)))
        call nc_get_variable(ncBckID, 'uReconstructMeridional', vcell)
        call nc_get_variable(ncAnlID, 'uReconstructMeridional', vcell_dart)
        vcell = vcell_dart - vcell ! v increments
   
        u = u + uv_increments_cell_to_edges(ucell, vcell)

        call nc_put_variable(ncBckID, 'u', u)
        call nc_put_variable(ncBckID, 'uReconstructZonal', ucell_dart)
        call nc_put_variable(ncBckID, 'uReconstructMeridional', vcell_dart)

        deallocate(u, ucell, vcell, ucell_dart, vcell_dart)
         
     else
        ! The state vector has updated zonal and meridional wind components.
        ! put them directly into the arrays.  These are the full values, not
        ! just increments.
        dims = get_dim_lengths(dom_id, get_varid_from_varname(dom_id, 'u'))
        allocate(u(dims(1), dims(2)))
        dims = get_dim_lengths(dom_id, get_varid_from_varname(dom_id, 'uReconstructZonal'))
        allocate(ucell_dart(dims(1), dims(2)))
        call nc_get_variable(ncAnlID, 'uReconstructZonal', ucell_dart)
        dims = get_dim_lengths(dom_id, get_varid_from_varname(dom_id, 'uReconstructMeridional'))
        allocate(vcell_dart(dims(1), dims(2)))
        call nc_get_variable(ncAnlID, 'uReconstructMeridional', vcell_dart)
 
        call uv_field_cell_to_edges(ucell_dart, vcell_dart, u)

        call nc_put_variable(ncBckID, 'u', variable)
        call nc_put_variable(ncBckID, 'uReconstructZonal', variable)
        call nc_put_variable(ncBckID, 'uReconstructMeridional', variable)

        deallocate(u, ucell_dart, vcell_dart)

     endif

  else
     ! copy u from analysis to background
    allocate(variable(get_variable_size(dom_id, get_varid_from_varname(dom_id, 'u'))))
    call nc_get_variable_size(ncAnlId, get_variable_name(dom_id,i), dimlens)
    call nc_get_variable(ncAnlID, 'u', variable, nc_count=dimlens)
    call nc_put_variable(ncBckID, 'u', variable, nc_count=dimlens)  

    deallocate(variable)

  endif

  call error_handler(E_MSG, 'Overwritten states in ',trim(next_outfile), source)

  call print_date( model_time,'update_mpas_states:model date')
  call print_time( model_time,'update_mpas_states:model time')
  call print_date( model_time,'update_mpas_states:model date',logfileunit)
  call print_time( model_time,'update_mpas_states:model time',logfileunit)

  call nc_close_file(ncAnlID,'update_mpas_states')
  call nc_close_file(ncBckID,'update_mpas_states')

  filenum = filenum + 1

end do fileloop

call finalize_utilities()

end program update_mpas_states

