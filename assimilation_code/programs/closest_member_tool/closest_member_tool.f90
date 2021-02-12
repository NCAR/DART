! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Select the member closest to the ensemble mean.
!>
!> This program has options to compute <em> distance </em> in several different ways
!> and returns the ensemble member which has the smallest total distance from
!> the ensemble mean.

program closest_member_tool

use types_mod,            only : r8, i8, obstypelength, MAX_NUM_DOMS, MAX_FILES

use time_manager_mod,     only : time_type, set_time_missing, operator(/=), &
                                 print_time
 
use utilities_mod,        only : find_namelist_in_file,        &
                                 error_handler, nmlfileunit, E_MSG, E_ERR,      &
                                 check_namelist_read, do_nml_file, do_nml_term, &
                                 open_file, close_file, set_multiple_filename_lists, &
                                 get_next_filename

use  location_mod,        only : location_type

use  obs_kind_mod,        only : get_num_quantities, get_index_for_quantity, &
                                 get_name_for_quantity

use  sort_mod,            only : index_sort

use assim_model_mod,      only : static_init_assim_model, get_model_size, &
                                 get_state_meta_data

use state_vector_io_mod,  only : read_state

use io_filenames_mod,     only : file_info_type, io_filenames_init,        &
                                 set_io_copy_flag, set_file_metadata,      &
                                 set_member_file_metadata, file_info_dump, &
                                 stage_metadata_type, get_stage_metadata,  &
                                 get_restart_filename, READ_COPY

use state_structure_mod,  only : get_num_domains

use mpi_utilities_mod,    only : initialize_mpi_utilities, task_count, &
                                 finalize_mpi_utilities, my_task_id,   &
                                 send_sum_to, sum_across_tasks

use ensemble_manager_mod, only : ensemble_type, init_ensemble_manager, compute_copy_mean, &
                                 get_my_vars, get_my_num_vars, end_ensemble_manager

implicit none

character(len=*), parameter :: source = 'closest_member_tool.f90'

integer               :: iunit, io, ens, i, j, total_j, qtyindex
integer               :: num_qtys, stype

integer(i8)           :: ii, model_size
integer, allocatable  :: index_list(:)
integer, parameter    :: max_list_len = 500
character(len=512)    :: msgstring, msgstring1
logical               :: allqtys, done
logical, allocatable  :: useqty(:), useindex(:)
type(location_type)   :: loc
type(time_type)       :: member_time
type(file_info_type)  :: ens_file_info


character(len=64)     :: method_name(4) = (/     &
   "Simple Difference    ", &
   "Normalized Difference", &
   "Simple RMS Diff      ", &
   "Normalized RMS Diff  "  /)

!----------------------------------------------------------------
! These variables are namelist-controllable.
!
integer                        :: ens_size               = 1
integer                        :: difference_method      = 4
logical                        :: single_restart_file_in = .false.
character(len=256)             :: input_restart_file_list(MAX_NUM_DOMS) = ''
character(len=256)             :: input_restart_files(MAX_FILES)        = '' 
character(len=256)             :: output_file_name                      = 'closest_restart'
character(len=obstypelength)   :: use_only_qtys(max_list_len)           = ''

!----------------------------------------------------------------
! different methods to compute 'distance' from mean:
!  1 = simple absolute difference
!  2 = normalized absolute difference
!  3 = simple rms difference
!  4 = normalized rms difference
!
!  (suggest more...)
!----------------------------------------------------------------

namelist /closest_member_tool_nml/  &
   input_restart_files,     &
   input_restart_file_list, &
   output_file_name,        &
   ens_size,                &
   single_restart_file_in,  &
   difference_method,       &
   use_only_qtys         

type(ensemble_type)             :: ens_handle
character(len=256), allocatable :: file_array_input(:,:)
character(len=256)              :: my_base, my_desc
integer(i8), allocatable        :: vars_array(:)
integer(i8)                     :: owners_index
integer                         :: num_domains, imem
integer                         :: ENS_MEAN_COPY 
integer                         :: copies, my_num_vars, num_copies
real(r8), allocatable           :: total_diff(:)

!----------------------------------------------------------------
! program start
!----------------------------------------------------------------

call initialize_mpi_utilities('closest_member_tool')


! Read the namelist entry and print it
call find_namelist_in_file("input.nml", "closest_member_tool_nml", iunit)
read(iunit, nml = closest_member_tool_nml, iostat = io)
call check_namelist_read(iunit, io, "closest_member_tool_nml")

if (do_nml_file()) write(nmlfileunit, nml=closest_member_tool_nml)
if (do_nml_term()) write(     *     , nml=closest_member_tool_nml)

if (single_restart_file_in) then
   write(msgstring,  *) 'single_restart_file_in is not supported.'
   write(msgstring1, *) 'Please contact DART if you would like to use this capability.'
   call error_handler(E_ERR,msgstring,msgstring1)
endif

! Initialize the model so we can get the size.

call static_init_assim_model()
model_size = get_model_size()

write(msgstring, *) 'Model size/restart data length =', model_size
call error_handler(E_MSG,'',msgstring)
write(msgstring, *) 'Ensemble member count = ', ens_size
call error_handler(E_MSG,'',msgstring)
write(msgstring, *) 'Computing difference using method: '//trim(method_name(difference_method))
call error_handler(E_MSG,'',msgstring)

! Make space that is ensemble size and an extra copy for the mean
call init_ensemble_manager(ens_handle, ens_size+1, model_size)

num_domains = get_num_domains()

! Given either a vector of input_state_files or a text file containing
! a list of files, return a vector of files containing the filenames.
call set_multiple_filename_lists(input_restart_files(:), &
                                 input_restart_file_list(:), &
                                 num_domains, &
                                 ens_size, &
                                 'closest_member_tool', &
                                 'input_restart_files', &
                                 'input_state_file_list')

! be ens_size but rather a single file (or multiple files if more than one domain)
allocate(file_array_input(ens_size, num_domains))

file_array_input  = RESHAPE(input_restart_files,  (/ens_size,  num_domains/))

! read in the ensemble and the mean - always in a separate file
call io_filenames_init(ens_file_info, &
                       ncopies       = ens_size, &
                       cycling       = single_restart_file_in, &
                       single_file   = single_restart_file_in, &
                       restart_files = file_array_input)

do imem = 1, ens_size
   write(my_base,'(A,I0.2)') 'inens_',                 imem
   write(my_desc,'(A,I0.2)') 'input ensemble member ', imem
   call set_file_metadata(ens_file_info,                       &
                          cnum     = imem,                     &
                          fnames   = file_array_input(imem,:), &
                          basename = my_base,                  &
                          desc     = my_desc)

   call set_io_copy_flag(ens_file_info,      &
                         cnum    = imem,     &
                         io_flag = READ_COPY)
enddo

! Read the ensemble from files
member_time = set_time_missing()
call read_state(ens_handle, ens_file_info, read_time_from_file=.false., model_time=member_time)

! Compute mean
ENS_MEAN_COPY = ens_size + 1
call compute_copy_mean(ens_handle, 1, ens_size, ENS_MEAN_COPY)

allocate(index_list(ens_size))

! are we adding up only differences from particular quantities, or the entire 
! vector?
if (use_only_qtys(1) /= '') then
   allqtys = .false.

   num_qtys = get_num_quantities()
   allocate(useqty(0:num_qtys))
   useqty = .false.

   done = .false.
   QtyList:do i=1, max_list_len
      if (use_only_qtys(i) == '') then
         done = .true.
         exit QtyList
      endif
      qtyindex = get_index_for_quantity(use_only_qtys(i))   
      if (qtyindex < 0) then
         write(msgstring, *) 'unrecognized QTY string: '//trim(use_only_qtys(i))
         call error_handler(E_ERR,'closest_member_tool', msgstring, source)
      endif
      useqty(qtyindex) = .true.

   enddo QtyList

   if (.not. done) then
      write(msgstring, *) 'cannot have more than ', max_list_len, ' qtys'
      call error_handler(E_ERR,'closest_member_tool', msgstring, source)
   endif

   write(msgstring, *) 'Computing difference based only on items in state vector items of quantity:'
   call error_handler(E_MSG,'',msgstring)
   do i=0, num_qtys
      if (useqty(i)) then
         write(msgstring, *) '   ', trim(get_name_for_quantity(i))
         call error_handler(E_MSG,'',msgstring)
      endif
   enddo

else
   allqtys = .true.
endif

! if we are not processing all qtys of state vector items, set up a mask
! for the qtys we are.  do this once at the start so we don't replicate
! work.  the useqty(number_of_different_qtys) array says whether we are
! going to add differences for this type.   the useindex(state_vector_len)
! array is being precomputed here so it's fast to loop over the ensemble
! members and only add up differences for the qtys of interest.

my_num_vars = get_my_num_vars(ens_handle)
num_copies  = ens_handle%num_copies
allocate(useindex(my_num_vars), vars_array(my_num_vars))
call get_my_vars(ens_handle, vars_array)

if (.not. allqtys) then
   useindex(:) = .false.

   j = 0
   do ii=1, my_num_vars
      owners_index = vars_array(ii)
      call get_state_meta_data(owners_index, loc, stype)
      if (stype < 0 .or. stype > num_qtys) then
         write(msgstring, *) 'bad QTY from get_state_meta_data, ', stype, ' for index ', owners_index
         write(msgstring1, *) 'must be between 0 and ', num_qtys
         call error_handler(E_ERR,'closest_member_tool', msgstring, &
                            source, text2=msgstring1)

      endif
      
      if (useqty(stype)) then 
         useindex(ii) = .true.
         j = j + 1
      endif
   enddo

   ! compute the total across all members
   call sum_across_tasks(j, total_j)
   write(msgstring, *) 'using ', total_j, ' of ', model_size, ' items in the state vector'
   call error_handler(E_MSG,'closest_member_tool', msgstring)
else
   ! use everything.
   useindex(:) = .true.
endif

allocate(total_diff(ens_size))

total_diff = compute_diff(ens_handle%copies(:,:), ens_handle%copies(ENS_MEAN_COPY,:))

!------------------- Print out results     -----------------------

if (my_task_id() == 0) then
   call index_sort(total_diff, index_list, ens_size)
   call error_handler(E_MSG, '', ' ')
   write(msgstring, "(A,I5)") 'Member with the minimum difference from the mean is ', index_list(1)
   call error_handler(E_MSG, '', msgstring)
   call error_handler(E_MSG, '', ' ')

   do ens=1, ens_size
      write(msgstring, "(A,I5,A,G18.6)") "Member ", index_list(ens), " difference ", total_diff(index_list(ens))
      call error_handler(E_MSG, '', msgstring)
   enddo

   !------------------- Write results to file -----------------------
   
   ! if the input is a single file, write the ensemble member number to a file.
   ! if the input is separate files, write the full filename to a file.
   
   iunit = open_file(output_file_name, 'formatted', 'write')
   
   if (single_restart_file_in) then
      write(iunit, "(I6)") index_list(1)
   else
      !> @todo FIXME is this domain by domain?  if so, need to loop over domains?
      msgstring = get_next_filename(input_restart_file_list(1), index_list(1))
      write(iunit, "(A)") trim(msgstring)
   endif
   
   call close_file(iunit)
     
   call error_handler(E_MSG, '', ' ')
   write(msgstring, *) 'Writing closest member information to file: ', trim(output_file_name)
   call error_handler(E_MSG, '', msgstring)

endif

!------------------- Write results to file -----------------------

deallocate(index_list, useindex, vars_array)
deallocate(total_diff)
if (.not. allqtys) deallocate(useqty)

call end_ensemble_manager(ens_handle)
call finalize_mpi_utilities()

!----------------------------------------------------------------
!----------------------------------------------------------------

contains

function compute_diff(ens_mems, ens_mean)
real(r8), intent(in) :: ens_mems(num_copies,my_num_vars)
real(r8), intent(in) :: ens_mean(my_num_vars)
real(r8) :: compute_diff(ens_size)

real(r8), allocatable :: local_diffs(:), adiff(:), global_diff(:)

allocate(adiff(my_num_vars))
allocate(local_diffs(ens_size), global_diff(ens_size))

do copies = 1, ens_size
   select case (difference_method)
   
      ! simple absolute difference
      case (1)
         where(useindex) 
            adiff(:) = abs(ens_mean(:) - ens_mems(copies,:))
         endwhere

         local_diffs(copies) = sum(adiff(:))

      ! normalized absolute difference
      case (2)
   
         where (useindex) 
            where (ens_mean(:) /= 0.0_r8) 
               adiff(:) = abs((ens_mean(:) - ens_mems(copies,:))/ens_mean(:))
            elsewhere
               adiff(:) = abs(ens_mems(copies,:))
            endwhere
         endwhere
   
         local_diffs(copies) = sum(adiff(:))

       ! simple rms difference
       case (3)
   
          where (useindex) adiff(:) = ens_mean(:) - ens_mems(copies, :)
   
          local_diffs(copies) = sum(adiff * adiff)

       ! normalized rms difference
       case (4)
   
          where (useindex) 
             where (ens_mean(:) /= 0.0_r8) 
                adiff = (ens_mean(:) - ens_mems(copies,:)) / ens_mean(:)
             elsewhere
                adiff = ens_mems(copies,:)
             endwhere
          endwhere
   
          local_diffs(copies) = sum(adiff * adiff)
   
      case default
         write(msgstring, *) 'Valid values for difference_method are 1-4, value is', difference_method
         call error_handler(E_ERR,'closest_member_tool','Bad value for difference_method', &
                            source, text2=msgstring)
   end select
enddo

call send_sum_to(local_diffs, 0, global_diff)

if(my_task_id() == 0) then
   ! need to square total difference after local sums are computed
   if (difference_method == 3 .or. difference_method  == 4) then 
      where(global_diff > 0) global_diff = sqrt(global_diff)
   endif
endif

compute_diff = global_diff

deallocate(adiff, global_diff, local_diffs)

end function compute_diff

end program closest_member_tool

