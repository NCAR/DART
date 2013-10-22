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
    subroutine MPP_WRITE_2DDECOMP_2D_( unit, field, domain, data, tstamp )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain2D), intent(inout) :: domain
      MPP_TYPE_, intent(inout) :: data(:,:)
      real(DOUBLE_KIND), intent(in), optional :: tstamp
      MPP_TYPE_ :: data3D(size(data,1),size(data,2),1)
#ifdef use_CRI_pointers
      pointer( ptr, data3D )
      ptr = LOC(data)
#else
      data3D = RESHAPE( data, SHAPE(data3D) )
#endif
      call mpp_write( unit, field, domain, data3D, tstamp )
      return
    end subroutine MPP_WRITE_2DDECOMP_2D_

    subroutine MPP_WRITE_2DDECOMP_3D_( unit, field, domain, data, tstamp )
!mpp_write writes <data> which has the domain decomposition <domain>
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain2D), intent(inout) :: domain !must have intent(out) as well because active domain might be reset
      MPP_TYPE_, intent(inout) :: data(:,:,:)
      real(DOUBLE_KIND), intent(in), optional :: tstamp
!cdata is used to store compute domain as contiguous data
!gdata is used to globalize data for multi-PE single-threaded I/O
      MPP_TYPE_, allocatable, dimension(:,:,:) :: cdata, gdata
!NEW: data may be on compute OR data domain
      logical :: data_has_halos, halos_are_global, x_is_global, y_is_global
      integer :: is, ie, js, je, isd, ied, jsd, jed, isg, ieg, jsg, jeg

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_WRITE: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE: invalid unit number.' )

      call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
      call mpp_get_data_domain   ( domain, isd, ied, jsd, jed, x_is_global=x_is_global, y_is_global=y_is_global )
      call mpp_get_global_domain ( domain, isg, ieg, jsg, jeg )
      if( size(data,1).EQ.ie-is+1 .AND. size(data,2).EQ.je-js+1 )then
          data_has_halos = .FALSE.
      else if( size(data,1).EQ.ied-isd+1 .AND. size(data,2).EQ.jed-jsd+1 )then
          data_has_halos = .TRUE.
      else
          call mpp_error( FATAL, 'MPP_WRITE: data must be either on compute domain or data domain.' )
      end if
      halos_are_global = x_is_global .AND. y_is_global
      if( npes.GT.1 .AND. mpp_file(unit)%threading.EQ.MPP_SINGLE )then
          if( halos_are_global )then
              call mpp_update_domains( data, domain )
!all non-0 PEs have passed their data to PE 0 and may now exit
              if( pe.NE.0 )return
              call write_record( unit, field, size(data), data, tstamp )
          else
!put field onto global domain
              allocate( gdata(isg:ieg,jsg:jeg,size(data,3)) )
              call mpp_global_field( domain, data, gdata )
!all non-0 PEs have passed their data to PE 0 and may now exit
              if( pe.NE.0 )return
              call write_record( unit, field, size(gdata), gdata, tstamp )
          end if
      else if( data_has_halos )then
!store compute domain as contiguous data and pass to write_record
          allocate( cdata(is:ie,js:je,size(data,3)) )
          cdata(:,:,:) = data(is-isd+1:ie-isd+1,js-jsd+1:je-jsd+1,:)
          call write_record( unit, field, size(cdata), cdata, tstamp, domain )
      else
!data is already contiguous
          call write_record( unit, field, size(data), data, tstamp, domain )
      end if

      return
    end subroutine MPP_WRITE_2DDECOMP_3D_
