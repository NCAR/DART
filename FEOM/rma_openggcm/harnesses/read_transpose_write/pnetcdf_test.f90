module pnetcdf_test
! Aim to read a netcdf restart file. 
! the local directory. This is because you don't want to allocate
! state_ens_handle%vars because you will run out of memory

use types_mod,            only : r8, missing_r8 
use mpi_utilities_mod,    only : initialize_mpi_utilities,     &
                                 finalize_mpi_utilities,       &
                                 task_count, my_task_id
use utilities_mod,        only : nc_check
use model_mod,            only : fill_variable_list,           &
                                 info_file_name,               &
                                 construct_file_name_in
use ensemble_manager_mod, only : ensemble_type, print_ens_handle
use io_filenames_mod,     only : restart_files_in
use assim_tools_mod,      only : test_state_copies

use mpi
use pnetcdf

implicit none

public :: parallel_read_transpose

contains

subroutine parallel_read_transpose(state_ens_handle, restart_in_file_name, domain, dart_index)

type(ensemble_type), intent(inout) :: state_ens_handle
character(len=129),  intent(in)    :: restart_in_file_name
integer,             intent(in)    :: domain
integer,             intent(inout) :: dart_index !< This is for mulitple domains


integer(KIND=MPI_OFFSET_KIND), allocatable :: stride(:)


real(r8), allocatable ::  data_block(:)

! ncfile
integer            :: ncfile, ndims, nvars, ngatts, unlimited 
integer            :: var_ndims
integer            :: var_size
integer            :: var_dimids(3) !dimensions ids (i,j,k)
integer            :: var_ids(5)
character(len=256) :: var_names(5)

integer :: j, i, k

integer(KIND=MPI_OFFSET_KIND),allocatable :: dim_sizes(:)
character(len=256)           ,allocatable :: dim_names(:)

integer :: total_var_size(5)
integer :: numdims
integer :: numvars

integer(KIND=MPI_OFFSET_KIND),allocatable :: start(:)
integer(KIND=MPI_OFFSET_KIND),allocatable :: count(:)

!character(len=256)  :: dir
character(len=1024) :: filename
character(len=256)  :: variable_list(5)

real(r8)                      :: begin

integer :: num_copies
integer :: my_pe
integer :: icopy
integer :: istart
integer :: iend

!call initialize_mpi_utilities('test io')
if(my_task_id()==0) then
call print_ens_handle(state_ens_handle)
endif

num_copies = state_ens_handle%num_copies
my_pe      = state_ens_handle%my_pe

!dir      = '/glade/p/work/hendric/test_cases/pop_large_test_case/'
filename = trim(info_file_name(1))

! specific for POP
! SATL_CUR, TEMP,CUR, UVEL_VUR, VVEL_CUR, PSURC_CUR
numvars   = 5 
var_names = fill_variable_list(numvars)

istart = 1

COPIES: do icopy = 1,num_copies
   begin = MPI_WTIME()

   filename = restart_files_in(icopy, domain)

   if(my_task_id()==0) write(*,*) "opening : " , trim(filename)
   
   ! Open netcdf file
   call nc_check(nf90mpi_open(MPI_COMM_WORLD,trim(filename),       &
                 NF90_NOWRITE + NF90_64BIT_OFFSET, MPI_INFO_NULL,  &
                 ncfile),'open error')
   
   ! inquire general information
   call nc_check(nfmpi_inq(ncfile, ndims=ndims, nvars=nvars,       &
                 ngatts=ngatts,unlimdimid=unlimited), 'nf90mpi_inq')
   
   if(my_task_id()==0) then
      121 format('ndims [',I2,'], nvars [',I3,'], ngatts[',I4,']')
      print 121, ndims, nvars, ngatts
      print *, " "
   endif
   
   ! inquire dimension information
   allocate(dim_sizes(ndims))
   allocate(dim_names(ndims))
   
   do i = 1,ndims
     call nc_check(nf90mpi_inquire_dimension(ncfile, i,            &
                   name=dim_names(i), len=dim_sizes(i)),'dim_len')
     if(my_task_id()==0) then
       120 format('dim ',A,' size = ',I4) 
       print 120, trim(dim_names(i)), dim_sizes(i) 
     endif
   enddo
   
   if(my_task_id() == 0) print*, " "
   
   ! get the var_names ids
   do i = 1,numvars
     call nc_check(nf90mpi_inq_varid(ncfile, var_names(i), &
                   varid=var_ids(i)), 'inq_varid')
     if(my_task_id() == 0) then
       123 format('var_id[',I2,'] variable_name = ',A9)
       print 123, var_ids(i), trim(var_names(i))
     endif
   enddo
   
   if(my_task_id() == 0) print*, " "
   
   do i = 1,numvars
     call nc_check(nf90mpi_inquire_variable(ncfile, var_ids(i), &
                   ndims=var_ndims, dimids=var_dimids),'inq_var')
   
     allocate(start(var_ndims))
     allocate(count(var_ndims))
     allocate(stride(var_ndims))
   
     !! hyperslab read
     !k = 2
     !start(k) = (dim_sizes(var_dimids(k))/task_count())*(my_task_id())+1
     !count(k) = (dim_sizes(var_dimids(k))/task_count())
     !var_size = count(k)
   
     !do j = 1,var_ndims
     !  if (j /= k) then
     !    stride(j) = task_count()
     !    start(j) = 1
     !    count(j) = dim_sizes(var_dimids(j))
     !    var_size = var_size*count(j)
     !  endif
     !enddo
   
     ! strided read
     k = 1
     start(k)  = my_task_id()+1
     count(k)  = (dim_sizes(var_dimids(k))/task_count())
     stride(k) = task_count()
     var_size  = count(k)
   
     do j = 1,var_ndims
       if (j /= k) then
         stride(j) = 1
         start(j)  = 1
         count(j)  = dim_sizes(var_dimids(j))
         var_size  = var_size*count(j)
       endif
     enddo

     iend = istart+var_size
     
     allocate(data_block(var_size))
     call nc_check(nfmpi_get_vars_double_all(ncfile,var_ids(i), &
                   start=start,count=count,stride=stride,       &
                   dvals=data_block),'get')
   
     state_ens_handle%copies(icopy,istart:iend) = data_block 
     !if (my_task_id()==0) then
     !  open(unit=12,FILE="x2.txt",FORM="FORMATTED",STATUS="OLD", &
     !       ACTION="WRITE")
     !  do j = 1, var_size
     !     write(12,*) data_block(j) 
     !  end do
     !endif
     istart = iend
   
     deallocate(start)
     deallocate(count)
     deallocate(stride)
     deallocate(data_block)
   
   enddo
   
   if(my_task_id() == 0) then
     888 format('MPI Wall Time = ',F5.2, 'sec')
     print 888, MPI_WTIME()-begin
     print *, " " 
   endif
   
   
   ! Get the variable dimensions
   do i = 1,numvars
     call nc_check(nf90mpi_inquire_variable(ncfile, var_ids(i), &
                   ndims=numdims), 'inq_dimid')
   
     total_var_size(i) = product(dim_sizes(1:numdims))
   
     if(my_task_id() == 0) then
       124 format('var_id[',I2,'] num_dims = ',I1,' total_var_size = ',I9)
       print 124, var_ids(i), numdims, total_var_size(i)
     endif
   enddo
   
   deallocate(dim_sizes)
   deallocate(dim_names)

enddo COPIES
     
  i = my_task_id()+1
  write(filename,'(A,I1.1)') 'p.',i
  open(unit=i,FILE=filename,FORM="FORMATTED",STATUS="OLD", &
       ACTION="WRITE")
  do j = 1, state_ens_handle%my_num_vars
     write(i,*) state_ens_handle%copies(1,j)
  end do
  close(unit=i)

end subroutine parallel_read_transpose

end module pnetcdf_test
