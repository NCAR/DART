! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module pio_transpose_mod
! Aim to read a netcdf restart file. 
! the local directory. This is because you don't want to allocate
! state_ens_handle%vars because you will run out of memory

use types_mod,            only : r8, missing_r8 
use mpi_utilities_mod,    only : task_count, my_task_id, get_dart_mpi_comm  
use ensemble_manager_mod, only : ensemble_type
use io_filenames_mod,     only : io_filenames_init, get_input_file, get_output_file
use ensemble_manager_mod,   only : copies_in_window

use state_structure_mod, only :  get_variable_name
use state_vector_io_mod

use mpi
use pio

implicit none

public :: pio_read_transpose

contains

subroutine pio_read_transpose(state_ens_handle)

type(ensemble_type), intent(inout) :: state_ens_handle

! netcdf variables
integer :: nvars, ngatts, unlimited !ncfile1

! variables for navigating dimensions and variables
integer :: numdims ! number of dimensions read
integer :: numvars ! number of variables read
integer :: numens ! number of ensembles

integer :: var_ndims ! variable dimensions
integer :: var_size  ! variable size
integer :: var_dimids(3) ! dimensions ids (i,j,k)->(1,2,3)
integer :: var_ids(5)    ! netcdf ids for variables   

character(len=256) :: var_name  ! variables names
character(len=256) :: filename  

integer,            allocatable :: dim_sizes(:)
character(len=256), allocatable :: dim_names(:)

! blocking variables for netcdf
integer :: bufcount
real(r8), allocatable:: var_block(:)

! indexing variables
integer :: i, j, k
integer :: icopy, istart, iend
integer :: start_task, start_index
integer :: remainder 
integer, allocatable :: start(:), count(:)


!PIO variables
integer :: num_iotasks
integer :: num_aggregator
integer :: stride
integer :: ret
integer, allocatable :: compdof(:)

type(var_desc_t) :: varid  
type(file_desc_t),     allocatable :: ncfile(:)
type(io_desc_t),       allocatable :: iodesc(:)
type(iosystem_desc_t), allocatable :: iosystem(:)

! MPI variables
integer :: dart_comm
integer :: num_tasks
integer :: rank

dart_comm  = get_dart_mpi_comm()
num_tasks  = task_count()
rank       = my_task_id()

! specific for POP
! SATL_CUR, TEMP,CUR, UVEL_VUR, VVEL_CUR, PSURC_CUR

numvars   = 5 
numens   = copies_in_window(state_ens_handle)

allocate(ncfile(numens))
allocate(iodesc(numens))
allocate(iosystem(numens))

! FIXME : need to find optimal processer layout for reads 
num_iotasks    = 1
num_aggregator = 1
stride         = 1

do icopy = 1,numens 
   call pio_init(rank, dart_comm, num_iotasks, num_aggregator, &
                 stride, PIO_rearr_box, iosystem(icopy),icopy-1)

   filename = get_input_file(icopy,1)
   if(rank == 0) print*, 'OPENING: ' , trim(filename)

   ret = pio_openfile(iosystem(icopy), ncfile(icopy), &
                      PIO_iotype_pnetcdf, filename)
enddo 

! find the number of dimensions in the file
ret = pio_inquire(ncfile(1), numdims)
      
allocate(dim_sizes(numdims))
allocate(dim_names(numdims))

! find the length of each dimension 
do i = 1,numdims
  ret = pio_inquire_dimension(ncfile(1), i, dim_names(i), dim_sizes(i))
  if(rank == 0) then
    120 format('dim ',A,' size = ',I4) 
    print 120, trim(dim_names(i)), dim_sizes(i) 
  endif
enddo

if(rank == 0) print*, " "
      
! index into copy array
istart = 1

! starting task for round robin
start_task = 0 

! read the variables in one at a time 
VAR_LOOP: do i = 1,numvars
     var_name = get_variable_name(1,i)
     ! get the var_name ids
     ret = pio_inq_varid(ncfile(1), var_name, var_ids(i))

     ! find number of dimensions and dimension id for variable
     ret  = pio_inquire_variable(ncfile(1), var_ids(i), &
                   ndims=var_ndims, dimids=var_dimids)
 
     allocate(count(var_ndims))
     
     ! calculate variable size
     var_size = 1   
     do j = 1,var_ndims
       count(j)  = dim_sizes(var_dimids(j))
       var_size  = var_size*count(j)
     enddo
     
     if(rank == 0) write(*,*) trim(var_name), var_size

     ! calculate amt information on each processor
     bufcount  = var_size/task_count()
     remainder = mod(var_size,task_count())
     
     ! index where to start in netcdf variable
     start_index = mod(num_tasks+rank-start_task, num_tasks) 
     
     ! add space for leftovers
     if (start_index < remainder) bufcount = bufcount + 1

     allocate(compdof(bufcount))
     allocate(var_block(bufcount))

     ! set up decomposition
     do j = 1,bufcount
        compdof(j) = start_index + task_count()*(j-1) + 1
     enddo

     ! grab varid for read 
     ret = pio_inq_varid(ncfile(1), var_name, varid)
        
     iend = istart+bufcount-1

     do icopy = 1,numens !COPY_LOOP

        ! round robin decomposition
        call pio_initdecomp(iosystem(icopy), PIO_double, count, &
                            compdof, iodesc(icopy))   
        
        ! index into copy array 
        call pio_read_darray(ncfile(icopy), varid, iodesc(icopy), &
               state_ens_handle%copies(icopy,istart:iend), ret)
        
        call pio_freedecomp(ncfile(icopy), iodesc(icopy))
        
     enddo

     ! update copy index and start task for next iteration
     istart     = iend + 1
     start_task = mod(start_task+remainder,task_count())

     deallocate(count)
     deallocate(compdof)
     deallocate(var_block)

enddo VAR_LOOP

deallocate(dim_sizes)
deallocate(dim_names)

! close the netcdf files
do icopy = 1,numens !COPY_LOOP
   call pio_closefile(ncfile(icopy))
enddo

end subroutine pio_read_transpose

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
subroutine pio_transpose_write(state_ens_handle)

type(ensemble_type), intent(inout) :: state_ens_handle

! netcdf variables
integer :: nvars, ngatts, unlimited

character(len=256) :: filename

! variables for navigating dimensions and variables
integer            :: numdims ! number of dimensions read
integer            :: numvars ! number of variables read
integer            :: numens ! number of ensembles

integer            :: var_ndims ! variable dimensions
integer            :: var_size  ! variable size
integer            :: var_dimids(3) ! dimensions ids (i,j,k)->(1,2,3)
integer            :: var_ids(5)    ! netcdf ids for variables   
character(len=256) :: var_name  ! variables names

integer,            allocatable :: dim_sizes(:)
character(len=256), allocatable :: dim_names(:)

! blocking variables for netcdf
integer :: bufcount
real(r8), allocatable:: var_block(:)

! indexing variables
integer :: remainder 
integer :: i, j, k
integer :: icopy, istart, iend
integer :: start_task, start_index

integer, allocatable :: start(:), count(:)

!PIO variables
integer :: num_iotasks
integer :: num_aggregator
integer :: ret
integer :: stride
integer, allocatable :: compdof(:)
real(r8) :: fillval

type(var_desc_t) :: varid  
type(file_desc_t),     allocatable :: ncfile(:)
type(io_desc_t),       allocatable :: iodesc(:)
type(iosystem_desc_t), allocatable :: iosystem(:)

! MPI variables
integer :: dart_comm
integer :: num_tasks
integer :: rank

dart_comm  = get_dart_mpi_comm()
num_tasks  = task_count()
rank       = my_task_id()
 
! specific for POP
! SATL_CUR, TEMP,CUR, UVEL_VUR, VVEL_CUR, PSURC_CUR
numvars   = 5 
numens   = copies_in_window(state_ens_handle)

allocate(ncfile(numens))
allocate(iodesc(numens))
allocate(iosystem(numens))

! FIXME : need to find optimal processer layout for reads 
num_iotasks    = 1
num_aggregator = 1
stride         = 1

do icopy = 1, numens
   call pio_init(rank, dart_comm, num_iotasks, num_aggregator, &
           stride, PIO_rearr_box, iosystem(icopy),icopy-1)
   
   filename = get_output_file(icopy,1,.true.)
   if(rank == 0) print*, 'WRITTING: ' , trim(filename)

   ret = pio_openfile(iosystem(icopy), ncfile(icopy), PIO_iotype_pnetcdf, filename, &
            PIO_WRITE)
enddo

! find the number of dimensions in the file
ret = pio_inquire(ncfile(1), numdims)
      
allocate(dim_sizes(numdims))
allocate(dim_names(numdims))

! find the length of each dimension 
do i = 1,numdims
  ret = pio_inquire_dimension(ncfile(1), i, dim_names(i), dim_sizes(i))
  if(rank == 0) then
    120 format('dim ',A,' size = ',I4) 
    print 120, trim(dim_names(i)), dim_sizes(i) 
  endif
enddo

if(rank == 0) print*, " "
      
! index into copy array
istart = 1

! starting task for round robin
start_task = 0 

! read the variables in one at a time 
VAR_LOOP: do i = 1,numvars
     
     var_name = get_variable_name(1,i)
     ! get the var_name ids
     ret = pio_inq_varid(ncfile(1), var_name, var_ids(i))

     ! find number of dimensions and dimension id for variable
     ret  = pio_inquire_variable(ncfile(1), var_ids(i), &
                   ndims=var_ndims, dimids=var_dimids)
 
     allocate(count(var_ndims))
     
     ! calculate variable size
     var_size = 1   
     do j = 1,var_ndims
       count(j)  = dim_sizes(var_dimids(j))
       var_size  = var_size*count(j)
     enddo
     
     if(rank == 0) write(*,*) trim(var_name), var_size

     ! calculate amt information on each processor
     bufcount  = var_size/task_count()
     remainder = mod(var_size,task_count())
     
     ! index where to start in netcdf variable
     start_index = mod(num_tasks+rank-start_task, num_tasks) 
     
     ! add space for leftovers
     if (start_index < remainder) bufcount = bufcount + 1

     allocate(compdof(bufcount))
     allocate(var_block(bufcount))

     do j = 1,bufcount
        compdof(j) = start_index + task_count()*(j-1) + 1
     enddo

     ! index into copy array 
     iend = istart+bufcount-1

     ! grab varid for write 
     ret = pio_inq_varid(ncfile(1), var_name, varid)
     
     do icopy = 1, numens

         ! round robin decomposition
         call pio_initdecomp(iosystem(icopy), PIO_double, count, &
                             compdof, iodesc(icopy))   
         
         fillval = -99999.0
         call pio_write_darray(ncfile(icopy), varid, iodesc(icopy), &
                 state_ens_handle%copies(icopy,istart:iend), ret)
     
         ! close the netcdf file
         call pio_freedecomp(ncfile(icopy), iodesc(icopy))
     enddo

     ! update copy index and start task for next iteration
     istart     = iend + 1
     start_task = mod(start_task+remainder,task_count())
     
     deallocate(count)
     deallocate(compdof)
     deallocate(var_block)

enddo VAR_LOOP

deallocate(dim_sizes)
deallocate(dim_names)

do icopy = 1, numens
   call pio_closefile(ncfile(icopy))
enddo

end subroutine pio_transpose_write

end module pio_transpose_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
