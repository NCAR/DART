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
    subroutine MPP_READ_2DDECOMP_2D_( unit, field, domain, data, tindex )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(inout) :: data(:,:)
      integer, intent(in), optional :: tindex
      MPP_TYPE_ :: data3D(size(data,1),size(data,2),1)
#ifdef use_CRI_pointers
      pointer( ptr, data3D )
      ptr = LOC(data)
      call mpp_read( unit, field, domain, data3D, tindex )
#else
      data3D = RESHAPE( data, SHAPE(data3D) )
      call mpp_read( unit, field, domain, data3D, tindex )
      data   = RESHAPE( data3D, SHAPE(data) )
#endif
      return
    end subroutine MPP_READ_2DDECOMP_2D_

    subroutine MPP_READ_2DDECOMP_3D_( unit, field, domain, data, tindex )
!mpp_read reads <data> which has the domain decomposition <domain>
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(inout) :: data(:,:,:)
      integer, intent(in), optional :: tindex
      MPP_TYPE_, allocatable :: cdata(:,:,:)
      MPP_TYPE_, allocatable :: gdata(:)
      integer :: lngth, lenx,leny,lenz,i,j,k,n
!NEW: data may be on compute OR data domain
      logical :: data_has_halos, halos_are_global, x_is_global, y_is_global
      integer :: is, ie, js, je, isd, ied, jsd, jed, isg, ieg, jsg, jeg, ioff, joff

      if (.NOT. present(tindex) .AND. mpp_file(unit)%time_level .ne. -1) &
      call mpp_error(FATAL, 'MPP_READ: need to specify a time level for data with time axis')

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_READ: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_READ: invalid unit number.' )

      call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
      call mpp_get_data_domain   ( domain, isd, ied, jsd, jed, x_is_global=x_is_global, y_is_global=y_is_global )
      call mpp_get_global_domain ( domain, isg, ieg, jsg, jeg )
      if( size(data,1).EQ.ie-is+1 .AND. size(data,2).EQ.je-js+1 )then
          data_has_halos = .FALSE.
      else if( size(data,1).EQ.ied-isd+1 .AND. size(data,2).EQ.jed-jsd+1 )then
          data_has_halos = .TRUE.
      else
          call mpp_error( FATAL, 'MPP_READ: data must be either on compute domain or data domain.' )
      end if
      halos_are_global = x_is_global .AND. y_is_global
      if( npes.GT.1 .AND. mpp_file(unit)%threading.EQ.MPP_SINGLE )then
          if( halos_are_global )then !you can read directly into data array
              if( pe.EQ.0 )call read_record( unit, field, size(data), data, tindex )
          else
              lenz=size(data,3)
              lngth=lenx*leny*lenz
              allocate(gdata(lngth))          
! read field on pe 0 and pass to all pes
              if( pe.EQ.0 ) call read_record( unit, field, lngth, gdata, tindex )
! broadcasting global array, this can be expensive!          
              call mpp_transmit( gdata, lngth, ALL_PES, gdata, lngth, 0 )
              ioff = is; joff = js
              if( data_has_halos )then
                  ioff = isd; joff = jsd
              end if
              do k=1,size(data,3)
                 do j=js,je
                    do i=is,ie
                       n=(i-isg+1) + (j-jsg)*lenx + (k-1)*lenx*leny
                       data(i-ioff+1,j-joff+1,k)=gdata(n)
                    enddo
                 enddo
              enddo
              deallocate(gdata)
          end if
      else if( data_has_halos )then
! for uniprocessor or multithreaded read
! read compute domain as contiguous data

          allocate( cdata(is:ie,js:je,size(data,3)) )
          call read_record(unit,field,size(cdata),cdata,tindex,domain)

          data(is-isd+1:ie-isd+1,js-jsd+1:je-jsd+1,:) = cdata(:,:,:)
          deallocate(cdata)
      else
          call read_record(unit,field,size(data),data,tindex,domain)
      end if

      return
    end subroutine MPP_READ_2DDECOMP_3D_
