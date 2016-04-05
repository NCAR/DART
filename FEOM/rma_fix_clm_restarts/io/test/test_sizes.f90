! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> start of test of setting and getting variable sizes from the state struct.

program test_sizes

use        types_mod, only : r8, i8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             error_handler, E_ERR

use     state_structure_mod, only : add_domain,                 &
                                    get_domain_size,            &
                                    get_num_domains,            &
                                    get_variable_size,          &
                                    get_variable_name,          &
                                    get_num_variables,          &
                                    get_num_dims,               &
                                    get_dim_lengths,            &
                                    get_dim_length,             &
                                    add_dimension_to_variable,  &
                                    finished_adding_domain,     &
                                    state_structure_info



! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer :: ds, did1, did2


! main code here
 
! initialize the dart libs
call initialize_utilities('test_sizes')

did1 = add_domain(100_i8)

ds = get_domain_size(did1)
call failed((ds == 100), 'domain 1 size not 100')

did2 = add_domain(10_i8)

ds = get_domain_size(did2)
call failed((ds == 10), 'domain 2 size not 10')

call finalize_utilities()

! end of main code


contains

!----------------------------------------------------------------------

!> print a message and halt if the first argument isn't .true.

subroutine failed(is_ok, str)

logical,          intent(in) :: is_ok 
character(len=*), intent(in) :: str

if (is_ok) return

call error_handler(E_ERR, 'test code', str, source, revision, revdate)

end subroutine failed

!----------------------------------------------------------------------

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
