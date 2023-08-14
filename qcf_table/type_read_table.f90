program read_table

implicit none
type obs_error_info_type
   logical               :: bounded_below, bounded_above
   real                  :: lower_bound,   upper_bound
end type

type probit_inflation_type
   integer               :: dist_type
   logical               :: bounded_below, bounded_above
   real                  :: lower_bound,   upper_bound
end type

type probit_state_type
   integer               :: dist_type
   logical               :: bounded_below, bounded_above
   real                  :: lower_bound,   upper_bound
end type

type probit_extended_state_type
   integer               :: dist_type
   logical               :: bounded_below, bounded_above
   real                  :: lower_bound,   upper_bound
end type

type obs_inc_info_type
   integer :: filter_kind
   logical :: rectangular_quadrature, gaussian_likelihood_tails
   logical :: sort_obs_inc, spread_restoration
   logical :: bounded_below, bounded_above
   real :: lower_bound,   upper_bound
end type

type qcf_table_data_type
   type(obs_error_info_type) :: obs_error_info
   type(probit_inflation_type) :: probit_inflation
   type(probit_state_type) :: probit_state
   type(probit_extended_state_type) :: probit_extended_state
   type(obs_inc_info_type) :: obs_inc_info
end type

! Reads in the QCEFF input options from tabular data file
!character(len=50), intent(in) :: qcf_table_filename
!real(r8), intent(out) :: qcf_table_data
!real, dimension(:, :), allocatable :: qcf_table_data_rows
type(qcf_table_data_type), allocatable :: qcf_table_data(:)
character(len=129), dimension(:), allocatable :: rowheaders !!!!! might need to change len=30

integer, parameter :: fileid = 10 !file identifier
character(len=30), parameter :: tester_QTY = 'QTY_GPSRO'
integer :: QTY_loc(1)

!integer, parameter :: num_columns = 28
integer :: nlines
integer :: io
integer :: numrows
integer :: row

!real, dimension(1:num_columns, 1:num_rows) :: table_data
!integer :: table_data_1, table_data_2
character(len=30), dimension(4) :: header1
character(len=30), dimension(29) :: header2
!variables for table values ^^^ 

open(unit=fileid, file='cam_qcf_table.txt')
nlines = 0

do !do loop to get number of rows (or QTY's) in the table
  read(fileid,*,iostat=io)
  if(io/=0) exit
  nlines = nlines + 1
end do
close(fileid)

print*, nlines

numrows = nlines - 2
print *, 'numrows: ', numrows

allocate(qcf_table_data(numrows))
allocate(rowheaders(numrows))
write(*,*) shape(qcf_table_data)

open(unit=fileid, file='cam_qcf_table.txt')

read(fileid, *) header1
read(fileid, *) header2 !! skip the headers
Write(*, *) "header1: ", header1
Write(*, *) "header2: ", header2

do row = 1, numrows
   read(fileid, *) rowheaders(row), qcf_table_data(row)%obs_error_info%bounded_below, qcf_table_data(row)%obs_error_info%bounded_above, & 
                   qcf_table_data(row)%obs_error_info%lower_bound, qcf_table_data(row)%obs_error_info%upper_bound, qcf_table_data(row)%probit_inflation%dist_type, &
                   qcf_table_data(row)%probit_inflation%bounded_below, qcf_table_data(row)%probit_inflation%bounded_above, &
                   qcf_table_data(row)%probit_inflation%lower_bound, qcf_table_data(row)%probit_inflation%upper_bound, qcf_table_data(row)%probit_state%dist_type, &
                   qcf_table_data(row)%probit_state%bounded_below, qcf_table_data(row)%probit_state%bounded_above, &
                   qcf_table_data(row)%probit_state%lower_bound, qcf_table_data(row)%probit_state%upper_bound, qcf_table_data(row)%probit_extended_state%dist_type, &
                   qcf_table_data(row)%probit_extended_state%bounded_below, qcf_table_data(row)%probit_extended_state%bounded_above, &
                   qcf_table_data(row)%probit_extended_state%lower_bound, qcf_table_data(row)%probit_extended_state%upper_bound, &
                   qcf_table_data(row)%obs_inc_info%filter_kind, qcf_table_data(row)%obs_inc_info%rectangular_quadrature, &
                   qcf_table_data(row)%obs_inc_info%gaussian_likelihood_tails, qcf_table_data(row)%obs_inc_info%sort_obs_inc, &
                   qcf_table_data(row)%obs_inc_info%spread_restoration, qcf_table_data(row)%obs_inc_info%bounded_below, qcf_table_data(row)%obs_inc_info%bounded_above, &
                   qcf_table_data(row)%obs_inc_info%lower_bound, qcf_table_data(row)%obs_inc_info%upper_bound
                   
   write(*, *) "rowheader(", row, "): ", rowheaders(row)
   write(*, *) "qcf_table_data(", row, "): "
   write(*, *) qcf_table_data(row)%obs_error_info%bounded_below, qcf_table_data(row)%obs_error_info%bounded_above, &
                   qcf_table_data(row)%obs_error_info%lower_bound, qcf_table_data(row)%obs_error_info%upper_bound, qcf_table_data(row)%probit_inflation%dist_type, &
                   qcf_table_data(row)%probit_inflation%bounded_below, qcf_table_data(row)%probit_inflation%bounded_above, &
                   qcf_table_data(row)%probit_inflation%lower_bound, qcf_table_data(row)%probit_inflation%upper_bound, qcf_table_data(row)%probit_state%dist_type, &
                   qcf_table_data(row)%probit_state%bounded_below, qcf_table_data(row)%probit_state%bounded_above, &
                   qcf_table_data(row)%probit_state%lower_bound, qcf_table_data(row)%probit_state%upper_bound, qcf_table_data(row)%probit_extended_state%dist_type, &
                   qcf_table_data(row)%probit_extended_state%bounded_below, qcf_table_data(row)%probit_extended_state%bounded_above, &
                   qcf_table_data(row)%probit_extended_state%lower_bound, qcf_table_data(row)%probit_extended_state%upper_bound, &
                   qcf_table_data(row)%obs_inc_info%filter_kind, qcf_table_data(row)%obs_inc_info%rectangular_quadrature, &
                   qcf_table_data(row)%obs_inc_info%gaussian_likelihood_tails, qcf_table_data(row)%obs_inc_info%sort_obs_inc, &
                   qcf_table_data(row)%obs_inc_info%spread_restoration, qcf_table_data(row)%obs_inc_info%bounded_below, qcf_table_data(row)%obs_inc_info%bounded_above, &
                   qcf_table_data(row)%obs_inc_info%lower_bound, qcf_table_data(row)%obs_inc_info%upper_bound 
end do

close(fileid)

QTY_loc = findloc(rowheaders, tester_QTY)
write(*, *) 'findloc of GPSRO: ', QTY_loc(1)

deallocate(qcf_table_data, rowheaders)

end program read_table
