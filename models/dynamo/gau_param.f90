! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module gau_param

! Update vector with a Gaussian noise.
! x(i) *= (1. + G * gaussian)
 
use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

implicit none

integer, parameter   :: r8 = selected_real_kind(14,30)

integer,  parameter  :: model_size = 1
real(r8), parameter  :: G          = 0.005_r8
real(r8), parameter  :: deltat     = 0.01_r8
real(r8)             :: x(model_size)
character (len = 128):: file_in      = 'model.in'
character (len = 128):: file_out      = 'model.out'
integer              :: i, ts, td, errval
type(random_seq_type), save :: sr


integer, parameter   :: nr = 100



contains
!===================================================================



subroutine read_state( )

integer :: size, i, errval

open(11,file=trim(file_in),status='old',iostat=errval)
if ( errval .eq. 0 ) then
   read(11,*)
   read(11,*)ts, td
   do i = 1, model_size
     read(11,*) x(i)
   enddo
   close(11)
else
   ts = 0
   td = 0
   x = -77.6_r8
endif

end subroutine read_state

subroutine march_ahead()

integer :: i, getpid, tempr

tempr = mod(getpid(),54000)
call init_random_seq(sr, tempr)
do i = 1, model_size
   ! Advance the parameter with a Gaussian noise
   x(i) = x(i)*(1._r8 + G * random_gaussian(sr, 0.0_r8, 1.0_r8)) 
end do

end  subroutine march_ahead

subroutine model_output()

open(11,file=trim(file_out),status='unknown',iostat=errval)
write(11,*)ts, td
do i = 1, model_size
   write(11,*) x(i)
enddo
close(11)
end subroutine model_output

end module gau_param


program gau_param_main
use gau_param
implicit none

call read_state( )
call march_ahead( )
call model_output( )

end program gau_param_main

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
