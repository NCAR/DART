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
  function MPP_GLOBAL_REDUCE_2D_( domain, field, locus )
    MPP_TYPE_ :: MPP_GLOBAL_REDUCE_2D_
    type(domain2D), intent(in) :: domain
    MPP_TYPE_, intent(in) :: field(:,:)
    integer, intent(out), optional :: locus(2)
    MPP_TYPE_ :: field3D(size(field,1),size(field,2),1)
    integer :: locus3D(3)
#ifdef use_CRI_pointers
    pointer( ptr, field3D )
    ptr = LOC(field)
    if( PRESENT(locus) )then
        MPP_GLOBAL_REDUCE_2D_ = MPP_GLOBAL_REDUCE_3D_( domain, field3D, locus3D )
        locus = locus3D(1:2)
    else
        MPP_GLOBAL_REDUCE_2D_ = MPP_GLOBAL_REDUCE_3D_( domain, field3D )
    end if
#else
    field3D = RESHAPE( field, SHAPE(field3D) )
    if( PRESENT(locus) )then
        MPP_GLOBAL_REDUCE_2D_ = MPP_GLOBAL_REDUCE_3D_( domain, field3D, locus3D )
        locus = locus3D(1:2)
    else
        MPP_GLOBAL_REDUCE_2D_ = MPP_GLOBAL_REDUCE_3D_( domain, field3D )
    end if
#endif
    return
  end function MPP_GLOBAL_REDUCE_2D_

  function MPP_GLOBAL_REDUCE_3D_( domain, field, locus )
    MPP_TYPE_ :: MPP_GLOBAL_REDUCE_3D_
    type(domain2D), intent(in) :: domain
    MPP_TYPE_, intent(in) :: field(0:,0:,:)
    integer, intent(out), optional :: locus(3)
    MPP_TYPE_ :: local
    integer :: here, ioff, joff

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GLOBAL_REDUCE: You must first call mpp_domains_init.' )
    if( size(field,1).EQ.domain%x%compute%size .AND. size(field,2).EQ.domain%y%compute%size )then
!field is on compute domain
        ioff = domain%x%compute%begin
        joff = domain%y%compute%begin
    else if( size(field,1).EQ.domain%x%data%size .AND. size(field,2).EQ.domain%y%data%size )then
!field is on data domain
        ioff = domain%x%data%begin
        joff = domain%y%data%begin
    else
        call mpp_error( FATAL, 'MPP_GLOBAL_REDUCE_: incoming field array must match either compute domain or data domain.' )
    end if

!get your local max/min
    local = REDUCE_VAL_(field)
!find the global
    MPP_GLOBAL_REDUCE_3D_ = local
    call MPP_REDUCE_( MPP_GLOBAL_REDUCE_3D_, domain%list(:)%pe )
!find locus of the global max/min
    if( PRESENT(locus) )then
!which PE is it on? min of all the PEs that have it
        here = mpp_npes()+1
        if( MPP_GLOBAL_REDUCE_3D_.EQ.local )here = pe
        call mpp_min( here, domain%list(:)%pe )
!find the locus here
        if( pe.EQ.here )locus = REDUCE_LOC_(field)
        locus(1) = locus(1) + ioff
        locus(2) = locus(2) + joff
        call mpp_broadcast( locus, 3, here, domain%list(:)%pe )
    end if
    return
  end function MPP_GLOBAL_REDUCE_3D_

  function MPP_GLOBAL_REDUCE_4D_( domain, field, locus )
    MPP_TYPE_ :: MPP_GLOBAL_REDUCE_4D_
    type(domain2D), intent(in) :: domain
    MPP_TYPE_, intent(in) :: field(:,:,:,:)
    integer, intent(out), optional :: locus(4)
    MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4))
    integer :: locus3D(3)
#ifdef use_CRI_pointers
    pointer( ptr, field3D )
    ptr = LOC(field)
    if( PRESENT(locus) )then
        MPP_GLOBAL_REDUCE_4D_ = MPP_GLOBAL_REDUCE_3D_( domain, field3D, locus3D )
        locus(1:2) = locus3D(1:2)
        locus(3) = modulo(locus3D(3),size(field,3))
        locus(4) = (locus3D(3)-locus(3))/size(field,3) + 1
        if( locus(3).EQ.0 )then
            locus(3) = size(field,3)
            locus(4) = locus(4) - 1
        end if
    else
        MPP_GLOBAL_REDUCE_4D_ = MPP_GLOBAL_REDUCE_3D_( domain, field3D )
    end if
#else
    field3D = RESHAPE( field, SHAPE(field3D) )
    if( PRESENT(locus) )then
        MPP_GLOBAL_REDUCE_4D_ = MPP_GLOBAL_REDUCE_3D_( domain, field3D, locus3D )
        locus(1:2) = locus3D(1:2)
        locus(3) = modulo(locus3D(3),size(field,3))
        locus(4) = (locus3D(3)-locus(3))/size(field,3) + 1
        if( locus(3).EQ.0 )then
            locus(3) = size(field,3)
            locus(4) = locus(4) - 1
        end if
    else
        MPP_GLOBAL_REDUCE_4D_ = MPP_GLOBAL_REDUCE_3D_( domain, field3D )
    end if
#endif
    return
  end function MPP_GLOBAL_REDUCE_4D_

  function MPP_GLOBAL_REDUCE_5D_( domain, field, locus )
    MPP_TYPE_ :: MPP_GLOBAL_REDUCE_5D_
    type(domain2D), intent(in) :: domain
    MPP_TYPE_, intent(in) :: field(:,:,:,:,:)
    integer, intent(out), optional :: locus(5)
    MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4)*size(field,5))
    integer :: locus3D(3)
#ifdef use_CRI_pointers
    pointer( ptr, field3D )
    ptr = LOC(field)
    if( PRESENT(locus) )then
        MPP_GLOBAL_REDUCE_5D_ = MPP_GLOBAL_REDUCE_3D_( domain, field3D, locus3D )
        locus(1:2) = locus3D(1:2)
        locus(3) = modulo(locus3D(3),size(field,3))
        locus(4) = modulo(locus3D(3),size(field,3)*size(field,4))
        locus(5) = (locus3D(3)-locus(4))/size(field,3)/size(field,4) + 1
        if( locus(3).EQ.0 )then
            locus(3) = size(field,3)
            locus(4) = locus(4) - 1
        end if
    else
        MPP_GLOBAL_REDUCE_5D_ = MPP_GLOBAL_REDUCE_3D_( domain, field3D )
    end if
#else
    field3D = RESHAPE( field, SHAPE(field3D) )
    if( PRESENT(locus) )then
        MPP_GLOBAL_REDUCE_5D_ = MPP_GLOBAL_REDUCE_3D_( domain, field3D, locus3D )
        locus(1:2) = locus3D(1:2)
        locus(3) = modulo(locus3D(3),size(field,3))
        locus(4) = modulo(locus3D(3),size(field,3)*size(field,4))
        locus(5) = (locus3D(3)-locus(4))/size(field,3)/size(field,4) + 1
        if( locus(3).EQ.0 )then
            locus(3) = size(field,3)
            locus(4) = locus(4) - 1
        end if
    else
        MPP_GLOBAL_REDUCE_5D_ = MPP_GLOBAL_REDUCE_3D_( domain, field3D )
    end if
#endif
    return
  end function MPP_GLOBAL_REDUCE_5D_
