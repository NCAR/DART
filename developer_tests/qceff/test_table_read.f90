program test_table_read

use algorithm_info_mod, only : init_algorithm_info_mod, end_algorithm_info_mod
use utilities_mod,      only : initialize_utilities

implicit none

character(len=129) :: qcf_table_filename

call initialize_utilities('test_table_read')

!n = command_argument_count()
call get_command_argument(1,qcf_table_filename)


!qcf_table_filename = 'qcf_table_v2.txt'

call init_algorithm_info_mod(qcf_table_filename)

call end_algorithm_info_mod()


end program test_table_read
