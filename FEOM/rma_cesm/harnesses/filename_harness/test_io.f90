!> Test harness for a limited transpose code
program test_io_read_transpose

use filter_mod,           only : filter_read_restart_direct, filter_set_initial_time, &
                                 filter_write_restart_direct
use types_mod,            only : r8, i8
use time_manager_mod,     only : time_type
use mpi_utilities_mod,    only : get_dart_mpi_comm, initialize_mpi_utilities, &
                                 finalize_mpi_utilities, datasize, my_task_id, task_sync
use ensemble_manager_mod, only : ensemble_type, init_ensemble_manager, end_ensemble_manager
use model_mod,            only : get_model_size, construct_file_name_in
use assim_model_mod,      only : static_init_assim_model
use utilities_mod,        only : find_namelist_in_file, check_namelist_read
use io_filenames_mod,     only : io_filenames_init,restart_files_in,set_filenames

use state_vector_io_mod
use state_structure_mod, only: get_num_domains, get_domain_size,&
                               get_num_variables, &
                               ! get_state_indices, get_num_variables, &
                               get_variable_name, static_init_state_type

use mpi

implicit none

type(ensemble_type) :: state_ens_handle
type(time_type)     :: time1

integer(i8)         :: model_size ! size of the whole state vector
integer             :: ens_size
logical             :: pio_read = .false.
integer              :: ind_in =  2993

! timing variables
real(r8)            ::read_time, write_time, start 

character(len = 129) :: inflation_in(2), inflation_out(2)

character(len = 129) :: filewrite
character(len = 129) :: info_file

integer             :: num_domains
integer             :: icopy
integer             :: domain_num
integer             :: ivar, idom

integer :: io, iunit

! get_state_indices variables
integer(i8) :: index_in
integer     :: lat_index, lon_index, depth_index
integer     :: var_type

namelist /test_io_nml/ ens_size, pio_read, ind_in

! initialize utilities
call initialize_mpi_utilities('test_io')

! read test_io namelist for ensemble size
call find_namelist_in_file('input.nml', 'test_io_nml', iunit)
read(iunit, nml = test_io_nml, iostat = io)
call check_namelist_read(iunit, io, 'test_io_nml')

write(*,*) "1"
! intialize model mod to get number of domains and model size
call static_init_assim_model()

write(*,*) "2"
call state_vector_io_init()

write(*,*) "3"
! initalize state type for netcdf variables
call static_init_state_type()

write(*,*) "4"
call init_ensemble_manager(state_ens_handle, ens_size , model_size)


! info_file = 'pop.r.nc'

! add_domain_from_info called in init_model
! domain_num = add_domain_from_info(info_file, num_variables, var_names)

! write(*,*) "---------------------------------------------------"
! 
! write(*,'(A)') "STATE VARIABLE INFO"
! write(*,'(A,i2.2)') "num_domains : ", get_num_domains()
! do idom=1,get_num_domains()
!     write(*,'(A,i2.2,A)') "domain[", idom, "]"
!     write(*,'(A,i10)')    "size = ", get_domain_size(idom)
!     write(*,'(A,i2.2)')   "num_variables : ", get_num_variables(idom)
!     do ivar=1,get_num_variables(idom)
!        write(*,'(A,i2.2,2A)') " variable[",ivar,"] = ", trim(get_variable_name(idom,ivar))
!     enddo
! enddo
! 
! write(*,*) "---------------------------------------------------"
! 
! write(*,*) "ENS SIZE =  " , ens_size
! write(*,*) "INFO FILE : ", info_file

! write(*,*) 'here'
! call get_state_indices(ind_in, lat_index, lon_index, depth_index, var_type)

! get filenames
write(*,*) "5"
call state_vector_io_init()

write(*,*) "6"
call io_filenames_init()

inflation_in(1) = 'inf_in1'
inflation_in(2) = 'inf_in2'
inflation_out(1) = 'inf_out1'
inflation_out(2) = 'inf_out2'

write(*,*) "7"
call set_filenames(state_ens_handle, ens_size, inflation_in, inflation_out)

!do icopy = 1,ens_size
!   write(*,*) trim(restart_files_in(1,1))
!enddo

call end_ensemble_manager(state_ens_handle)

write(*,*) "8"
call finalize_mpi_utilities(async=0)

end program test_io_read_transpose
