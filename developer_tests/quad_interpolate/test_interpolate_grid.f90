! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download


program test_interpolate_grid

use regular_grid_mod, only : create_grid0, create_field
use types_mod, only : r8, MISSING_R8, metadatalength
use functions_mod, only : f_sine, f_sum
use grid_mod, only : grid_type, dump_grid
use target_grid_mod, only : create_gridT
use utilities_mod, only : initialize_utilities, finalize_utilities, &
                          array_dump, register_module, do_nml_file, do_nml_term, &
                          nmlfileunit, find_namelist_in_file, &
                          check_namelist_read

use quad_interp_mod


implicit none

! irregular grids have lon(:,:) lat(:,:)
! so the split between reg and irreg has to
! be somewhere.  do we need separate progs
! for interp to reg and interp to irreg ?


! maybe grid0 should also be a grid type?

real(r8), allocatable :: lon(:)
real(r8), allocatable :: lat(:)
real(r8), allocatable :: field0(:,:), fieldT(:,:), field1(:,:)

real(r8) :: resolution = 30
integer  :: n
integer  :: i,j ! loop variables

integer :: case = 1
logical :: is_regular = .true.
integer :: debug = 0

character(len=256) :: target_filename
character(len=metadatalength) :: lon_name, lat_name
type(grid_type) :: grid0, gridT, grid1

namelist /test_interpolate_grid_nml/ is_regular, target_filename, &
         case, debug, resolution, lon_name, lat_name


call initialize_utilities("test_interpolate_grid")
call read_namelist()


! make source grid and give it data

call create_grid0(resolution, grid0, n)

print *, 'allocating field ', n, ' by ', n
allocate(field0(n,n))
call create_field(grid0, field0, f_sine)
call array_dump(field0, 4)

! (export FORT_FMT_RECL=10240 for easy print with ifort)


! make target grid, no data

call create_gridT(target_filename, lon_name, lat_name, is_regular, gridT)
call dump_grid(gridT)
allocate(fieldT(gridT%nlon, gridT%nlat))


! call quad utils to move data from src to dst
call do_reg_interp(grid0, field0, gridT, fieldT)
call array_dump(fieldT, 4)


! create a copy of the src grid with no data

allocate(field1(n,n))
field1(:,:) = MISSING_R8

! call quad utils to move data back from dst to src2
!call do_reg_interp(gridT, fieldT, grid1, field1)


! compare field0 and field1 and evaluate success or failure
! call evaluate_results(grid0, field0, grid1, field1)


call finalize_utilities()


contains

!------------------------------------------------------------------
subroutine read_namelist()          
                                 
integer :: iunit, io

! Read the namelist 
call find_namelist_in_file("input.nml", "test_interpolate_grid_nml", iunit)
read(iunit, nml = test_interpolate_grid_nml, iostat = io)
call check_namelist_read(iunit, io, "test_interpolate_grid_nml")

! Output the namelist values if requested
if (do_nml_file()) write(nmlfileunit, nml=test_interpolate_grid_nml)
if (do_nml_term()) write(     *     , nml=test_interpolate_grid_nml)

end subroutine read_namelist




end program test_interpolate_grid
