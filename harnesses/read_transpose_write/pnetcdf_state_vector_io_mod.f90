module pnetcdf_state_vector_io_mod
! Aim to read a netcdf restart file. 
! the local directory. This is because you don't want to allocate
! state_ens_handle%vars because you will run out of memory

use types_mod,            only : r8, missing_r8 
use mpi_utilities_mod,    only : initialize_mpi_utilities,     &
                                 finalize_mpi_utilities,       &
                                 task_count, my_task_id, task_sync
use utilities_mod,        only : nc_check
use model_mod,            only : fill_variable_list,           &
                                 info_file_name,               &
                                 construct_file_name_in,       &
                                 get_model_size,               &
                                 pert_model_state
use ensemble_manager_mod, only : ensemble_type, print_ens_handle
use io_filenames_mod,     only : restart_files_in
use state_vector_io_mod
use assim_tools_mod,      only : test_state_copies
use data_structure_mod,   only : copies_in_window

use mpi
use pnetcdf

implicit none

public :: single_read_transpose

contains

subroutine single_read_transpose(state_ens_handle, domain, dart_index)

integer varid0
type(ensemble_type), intent(inout) :: state_ens_handle
integer,             intent(in)    :: domain
integer,             intent(inout) :: dart_index !< This is for mulitple domains

! netcdf variables
integer            :: ncfile, nvars, ngatts, unlimited 
integer            :: err, arr


! variables for vard_read
integer, allocatable :: offsets(:), blocklens(:)
integer :: remainder, buftype
! variables for navigating dimensions and variables
integer            :: numdims ! number of dimensions read
integer            :: numvars ! number of variables read
integer            :: num_ens ! number of ensembles

integer            :: var_ndims ! variable dimensions
integer            :: var_size  ! variable size
integer            :: var_size_total(5) ! total size of each variable
integer            :: var_dimids(3) ! dimensions ids (i,j,k)->(1,2,3)
integer            :: var_ids(5)    ! netcdf ids for variables   
character(len=256) :: var_names(5)  ! variables names

integer(KIND=MPI_OFFSET_KIND),allocatable :: dim_sizes(:)
character(len=256)           ,allocatable :: dim_names(:)

! blocking variables for netcdf
integer(KIND=MPI_OFFSET_KIND) :: bufcount
integer(KIND=MPI_OFFSET_KIND),allocatable :: start(:)
integer(KIND=MPI_OFFSET_KIND),allocatable :: count(:)
integer(KIND=MPI_OFFSET_KIND), allocatable :: stride(:)

! files we want to read in
!character(len=256)  :: dir
character(len=1024) :: filename

! for timing
real(r8) :: begin_time
real(r8), allocatable:: buf(:)

! processor specific variables
integer :: num_copies
integer :: my_pe

! indexing variables
integer :: i, j, k
integer :: icopy
integer :: istart
integer :: iend

integer :: start_task, start_index

! timing variables
real :: send_time, read_time, read_time_start, total_time
real :: bandwidth
real :: size_file

logical :: interf_provided
integer :: domainW, dart_indexW
character(len=129) :: file_out

logical :: debug_print = .false.
 
   read_time  = 0.0

   num_ens      = copies_in_window(state_ens_handle)
   num_copies   = state_ens_handle%num_copies ! have the extras, incase you
   my_pe        = state_ens_handle%my_pe
   
   ! specific for POP
   ! SATL_CUR, TEMP,CUR, UVEL_VUR, VVEL_CUR, PSURC_CUR
   numvars   = 5  ! need to change to be more generic
   var_names = fill_variable_list(numvars)

   do icopy = 1,num_ens
      ! single restart file in
      filename = restart_files_in(icopy, 1)
      if (my_task_id()==0) write(*,*) 'OPENING: ' , trim(filename)

      ! open netcdf file 
      call nc_check(nf90mpi_open(MPI_COMM_WORLD,trim(filename),       &
                    NF90_NOWRITE+NF90_64BIT_OFFSET, MPI_INFO_NULL,  &
                    ncfile),'open error')

      ! inquire general information
      call nc_check(nfmpi_inq(ncfile, ndims=numdims, nvars=nvars,       &
                    ngatts=ngatts,unlimdimid=unlimited), 'nf90mpi_inq')
      
      if(my_task_id()==0 .and. debug_print) then
         121 format('numdims [',I2,'], nvars [',I3,'], ngatts[',I4,']')
         print 121, numdims, nvars, ngatts
         print *, " "
      endif
      
      ! inquire dimension information
      allocate(dim_sizes(numdims))
      allocate(dim_names(numdims))
      
      do i = 1,numdims
        call nc_check(nf90mpi_inquire_dimension(ncfile, i,            &
                      name=dim_names(i), len=dim_sizes(i)),'dim_len')
        if(my_task_id()==0 .and. debug_print) then
          120 format('dim ',A,' size = ',I4) 
          print 120, trim(dim_names(i)), dim_sizes(i) 
        endif
      enddo
      
      if(my_task_id() == 0 .and. debug_print) print*, " "
      
      istart = 1

      ! get the var_names ids
      do i = 1,numvars
        call nc_check(nf90mpi_inq_varid(ncfile, var_names(i), &
                      varid=var_ids(i)), 'inq_varid')
        if(my_task_id() == 0 .and. debug_print) then
          123 format('var_id[',I2,'] variable_name = ',A9)
          print 123, var_ids(i), trim(var_names(i))
        endif
      enddo
      
      if(my_task_id() == 0 .and. debug_print) print*, " "
  
      ! start round robbin at PE0  
      start_task = 0 

      ! read the variables in one at a time 
      do i = 1,numvars
        call nc_check(nf90mpi_inquire_variable(ncfile, var_ids(i), &
                      ndims=var_ndims, dimids=var_dimids),'inq_var')
       
        allocate(count(var_ndims))

        var_size = 1   
        do j = 1,var_ndims
          count(j)  = dim_sizes(var_dimids(j))
          var_size  = var_size*count(j)
        enddo
        
        if(my_task_id() ==0 .and. debug_print) print *, 'var_size : ', var_size, ' var_name : ', trim(var_names(i))

        bufcount  = var_size/task_count()
        remainder = mod(var_size,task_count())

        start_index = mod(task_count()+my_task_id()-start_task,task_count()) 
        
        if (start_index < remainder) bufcount = bufcount + 1

        allocate(blocklens(bufcount))
        allocate(offsets(bufcount))
        
        do j = 1,bufcount
           blocklens(j) = 1
           offsets(j) = start_index+task_count()*(j-1) 
        enddo

        buftype = MPI_REAL8 

        call mpi_type_indexed(bufcount,blocklens,offsets,buftype,arr,err)
        call mpi_type_commit(arr,err)

        iend = istart+bufcount-1

        call nc_check(nf90mpi_get_vard_all(ncfile,var_ids(i), arr, &
              state_ens_handle%copies(icopy,istart:iend), &
              bufcount, buftype),'nf90mpi_get_vard_all')
        
        call mpi_type_free(arr,err)

        istart = iend + 1
        
        start_task = mod(start_task+remainder,task_count())

        deallocate(blocklens)
        deallocate(offsets)
        deallocate(count)
      
      enddo
      
      deallocate(dim_sizes)
      deallocate(dim_names)

      ! close the netcdf file
      call nc_check(nf90mpi_close(ncfile), 'closing file')
   enddo 

   !! copy state to other copies
   !do i = 2,num_ens
   !   state_ens_handle%copies(i,:) = state_ens_handle%copies(1,:)
   !      
   !   !call pert_model_state(state_ens_handle%copies(i,:),&
   !   !                      state_ens_handle%copies(i,:), interf_provided)
   !enddo
   
   !domainW     = 1
   !dart_indexW = 1
   !file_out    = 'transpose_write.nc'
   
   !write(*,*) 'PNETCDF ------------------'
   
   !call transpose_write(state_ens_handle, file_out,domainW,dart_indexW,.false.)

end subroutine single_read_transpose

end module pnetcdf_state_vector_io_mod
