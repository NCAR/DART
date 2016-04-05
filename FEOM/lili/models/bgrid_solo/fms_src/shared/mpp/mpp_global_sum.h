!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function MPP_GLOBAL_SUM_( domain, field, flags )
    MPP_TYPE_ :: MPP_GLOBAL_SUM_
    type(domain2D), intent(in) :: domain
    MPP_TYPE_, intent(in) :: field(:,: MPP_EXTRA_INDICES_ )
    integer, intent(in), optional :: flags
    MPP_TYPE_, allocatable, dimension(:,:) :: field2D, global2D
    integer :: i,j, ioff,joff

    if( size(field,1).EQ.domain%x%compute%size .AND. size(field,2).EQ.domain%y%compute%size )then
!field is on compute domain
        ioff = -domain%x%compute%begin + 1
        joff = -domain%y%compute%begin + 1
    else if( size(field,1).EQ.domain%x%data%size .AND. size(field,2).EQ.domain%y%data%size )then
!field is on data domain
        ioff = -domain%x%data%begin + 1
        joff = -domain%y%data%begin + 1
    else
        call mpp_error( FATAL, 'MPP_GLOBAL_SUM_: incoming field array must match either compute domain or data domain.' )
    end if
    if( PRESENT(flags) )then
        if( flags.NE.BITWISE_EXACT_SUM )call mpp_error( FATAL, 'MPP_GLOBAL_SUM_: only valid flag is BITWISE_EXACT_SUM.' )
!this is bitwise exact across different PE counts.
        allocate( field2D(domain%x%compute%begin:domain%x%compute%end,domain%y%compute%begin:domain%y%compute%end) )
        allocate( global2D(domain%x%global%size,domain%y%global%size) )
        do j = domain%y%compute%begin, domain%y%compute%end
           do i = domain%x%compute%begin, domain%x%compute%end
              field2D(i,j) = sum( field(i+ioff:i+ioff,j+joff:j+joff MPP_EXTRA_INDICES_) )
           end do
        end do

        call mpp_global_field( domain, field2D, global2D )
        MPP_GLOBAL_SUM_ = sum(global2D)
        deallocate( field2D)
        deallocate(global2D)
    else
!this is not bitwise-exact across different PE counts
        MPP_GLOBAL_SUM_ = sum( field(domain%x%compute%begin+ioff:domain%x%compute%end+ioff, &
                                     domain%y%compute%begin+joff:domain%y%compute%end+joff MPP_EXTRA_INDICES_) )
        call mpp_sum( MPP_GLOBAL_SUM_, domain%list(:)%pe )
    end if

    return
  end function MPP_GLOBAL_SUM_
