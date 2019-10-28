! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> Aim: Unit test for create_diagnostic_structure and end_diagnostic_structure
!> Tests:
!>   * diag_id should be greater than the number of domains
!>   * The number of variables in the diangostic domain should be equal to the 
!>     sum of the number of variables in every domain
!>   * The variable and dimension names should stay the same if there is 1 domain 
!>     in the state
!>   * The variable and dimension names should be appended with _d0# where # is the 
!>     domain number if there are multiple domains.
!>   * If you end the diagnostic structure the number of variables should be 0.
!>   * If you end the diagnostic structure the size should be 0.
!>   * Test that you can call create and destroy multiple times without error.
!>   * Check that you can add_domains until you read max_num_domains then the code
!>     should error out (so you can't overwrite the diagnostic domain).

program test_diag_structure

use types_mod,           only : i8
use utilities_mod,       only : initialize_utilities, finalize_utilities
use state_structure_mod, only : create_diagnostic_structure, end_diagnostic_structure, &
                                add_domain, get_num_variables, get_num_dims, &
                                get_num_domains, get_dim_name, get_variable_name, &
                                get_domain_size
use assert_mod,          only : assert_equal, assert_greater, assert_not_equal

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

integer :: diag_id, domain_id
integer(i8) :: model_size, m
integer :: i, j, n

model_size = 44

! init library
call initialize_utilities('test_diag_structure')

! Add domain
domain_id = add_domain(model_size)

! Create diagnostic domain
diag_id = create_diagnostic_structure()

! Test diag_id is greater than num_domains
call assert_greater(diag_id, get_num_domains(), 'diag id')

! Test number of variables in diag domain is equal to the total 
! number of variables in the state
n = 0
do i = 1, get_num_domains()
   do j = 1, get_num_variables(i)
      n = n + 1
   enddo
enddo
call assert_equal(n, get_num_variables(diag_id), 'num vars in use')

print*, '---- single domain -----'
! Print out variable names
do i = 1, get_num_variables(diag_id)
   print*, 'var name: ', trim(get_variable_name(diag_id, i))
   print*, '   dim names: ', (trim(get_dim_name(diag_id, i, j)), j = 1, get_num_dims(diag_id, i))
enddo

! Test domain_size of diag domain is equal to SUM(domain_sizes)
m = 0
do i = 1, get_num_domains()
      m = m + get_domain_size(i)
enddo
call assert_equal(get_domain_size(diag_id), m, 'size of diag')

! End diagnostic domain
call end_diagnostic_structure()

! Try and access diagnostic domain
! Number of variables should be 0
call assert_equal(get_num_variables(diag_id), 0, 'num vars after ended')

! Size should be equal to zero
call assert_equal(get_domain_size(diag_id), int(0,i8), 'domain size after ended')

! Add more domains to the state
print*, '---- multiple domains -----'
domain_id = add_domain(model_size)
domain_id = add_domain(model_size)

! Create diagnostic domain
diag_id = create_diagnostic_structure()

! Test number of variables in diag domain is equal to the total 
! number of variables in the state
n = 0
do i = 1, get_num_domains()
   do j = 1, get_num_variables(i)
      n = n + 1
   enddo
enddo
call assert_equal(n, get_num_variables(diag_id), 'num vars in use')

! Test domain_size of diag domain is equal to SUM(domain_sizes)
m = 0
do i = 1, get_num_domains()
      m = m + get_domain_size(i)
enddo
call assert_equal(get_domain_size(diag_id), m, 'size of diag')


! Print out variable names
do i = 1, get_num_variables(diag_id)
   print*, 'var name: ', trim(get_variable_name(diag_id, i))
   print*, '   dim names: ', (trim(get_dim_name(diag_id, i, j)), j = 1, get_num_dims(diag_id, i))
enddo

!-----------------------------------------------------------------------------
! Other tests
!-----------------------------------------------------------------------------
! Create diagnostic domain more than once 
diag_id = create_diagnostic_structure()
diag_id = create_diagnostic_structure()

! Destroy diagnostic domain more than once
call end_diagnostic_structure()
call end_diagnostic_structure()

! Note this errors out so you should have this at the end
! Add domains until you reach max_num_domains
print *, 'this last test is expected to cause a fatal error:'
do i = 1, 20 ! so you don't end up in an infinite loop
   domain_id = add_domain(model_size)
enddo
call assert_not_equal(i, 21, 'not reached max_num_domains')

! This should die
!print*, get_num_dims(diag_id, 1)

call finalize_utilities()

!-----------------------------------------------------------------------------

end program test_diag_structure

