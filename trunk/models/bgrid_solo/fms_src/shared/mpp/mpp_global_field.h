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
    subroutine MPP_GLOBAL_FIELD_2D_( domain, local, global, flags )
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(:,:)
      MPP_TYPE_, intent(out) :: global(:,:)
      integer, intent(in), optional :: flags
      MPP_TYPE_ :: local3D (size( local,1),size( local,2),1)
      MPP_TYPE_ :: global3D(size(global,1),size(global,2),1)
#ifdef use_CRI_pointers
      pointer( lptr,  local3D )
      pointer( gptr, global3D )
      lptr = LOC( local)
      gptr = LOC(global)
      call mpp_global_field( domain, local3D, global3D, flags )
#else
      local3D = RESHAPE( local, SHAPE(local3D) )
      call mpp_global_field( domain, local3D, global3D, flags )
      global = RESHAPE( global3D, SHAPE(global) )
#endif
    end subroutine MPP_GLOBAL_FIELD_2D_

    subroutine MPP_GLOBAL_FIELD_3D_( domain, local, global, flags )
!get a global field from a local field
!local field may be on compute OR data domain
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(:,:,:)
      MPP_TYPE_, intent(out) :: global(domain%x%global%begin:,domain%y%global%begin:,:)
      integer, intent(in), optional :: flags
      integer :: i, j, k, m, n, nd, nwords, lpos, rpos, ioff, joff
      logical :: xonly, yonly
      MPP_TYPE_ :: clocal (domain%x%compute%size*    domain%y%compute%size*    size(local,3))
      MPP_TYPE_ :: cremote(domain%x%compute%max_size*domain%y%compute%max_size*size(local,3))
      integer :: words_per_long, stackuse
      character(len=8) :: text
#ifdef use_CRI_pointers
      pointer( ptr_local,  clocal  )
      pointer( ptr_remote, cremote )

      ptr_local  = LOC(mpp_domains_stack)
      ptr_remote = LOC(mpp_domains_stack(size(clocal)+1))
#endif

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GLOBAL_FIELD: must first call mpp_domains_init.' )
      xonly = .FALSE.
      yonly = .FALSE.
      if( PRESENT(flags) )then
          xonly = flags.EQ.XUPDATE
          yonly = flags.EQ.YUPDATE
          if( .NOT.xonly .AND. .NOT.yonly )call mpp_error( WARNING, 'MPP_GLOBAL_FIELD: you must have flags=XUPDATE or YUPDATE.' )
      end if

#ifdef use_CRI_pointers
      words_per_long = size(transfer(local(1,1,1),mpp_domains_stack))
      stackuse = (size(clocal)+size(cremote))*words_per_long
      if( stackuse.GT.mpp_domains_stack_size )then
          write( text, '(i8)' )stackuse
          call mpp_error( FATAL, &
               'MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
      end if
      ptr_local  = LOC(mpp_domains_stack)
      ptr_remote = LOC( mpp_domains_stack((size(clocal)+1)*words_per_long) )
      mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, stackuse )
#endif
      if( size(global,1).NE.domain%x%global%size .OR. size(global,2).NE.domain%y%global%size .OR. &
          size(local,3).NE.size(global,3) ) &
           call mpp_error( FATAL, 'MPP_GLOBAL_FIELD: incoming arrays do not match domain.' )
    if( size(local,1).EQ.domain%x%compute%size .AND. size(local,2).EQ.domain%y%compute%size )then
!local is on compute domain
        ioff = -domain%x%compute%begin + 1
        joff = -domain%y%compute%begin + 1
    else if( size(local,1).EQ.domain%x%data%size .AND. size(local,2).EQ.domain%y%data%size )then
!local is on data domain
        ioff = -domain%x%data%begin + 1
        joff = -domain%y%data%begin + 1
    else
        call mpp_error( FATAL, 'MPP_GLOBAL_FIELD_: incoming field array must match either compute domain or data domain.' )
    end if

! make contiguous array from compute domain
      m = 0
      do k = 1,size(local,3)
         do j = domain%y%compute%begin,domain%y%compute%end
            do i = domain%x%compute%begin,domain%x%compute%end
               m = m + 1
               clocal(m) = local(i+ioff,j+joff,k)
               global(i,j,k) = clocal(m) !always fill local domain directly
            end do
         end do
      end do

!fill off-domains (note loops begin at an offset of 1)
      if( xonly )then
          nd = size(domain%x%list)
          do n = 1,nd-1
             lpos = mod(domain%x%pos+nd-n,nd)
             rpos = mod(domain%x%pos   +n,nd)
             nwords = domain%x%list(rpos)%compute%size * domain%x%list(rpos)%compute%size * size(local,3)
             call mpp_transmit( clocal, size(clocal), domain%x%list(lpos)%pe, cremote, nwords, domain%x%list(rpos)%pe )
             m = 0
             do k = 1,size(global,3)
                do j = domain%y%compute%begin,domain%y%compute%end
                   do i = domain%x%list(rpos)%compute%begin,domain%x%list(rpos)%compute%end
                      m = m + 1
                      global(i,j,k) = cremote(m)
                   end do
                end do
             end do
          end do
          call mpp_sync_self(domain%x%list(:)%pe)
      else if( yonly )then
          nd = size(domain%y%list)
          do n = 1,nd-1
             lpos = mod(domain%y%pos+nd-n,nd)
             rpos = mod(domain%y%pos   +n,nd)
             nwords = domain%y%list(rpos)%compute%size * domain%y%list(rpos)%compute%size * size(local,3)
             call mpp_transmit( clocal, size(clocal), domain%y%list(lpos)%pe, cremote, nwords, domain%y%list(rpos)%pe )
             m = 0
             do k = 1,size(global,3)
                do j = domain%y%list(rpos)%compute%begin,domain%y%list(rpos)%compute%end
                   do i = domain%x%compute%begin,domain%x%compute%end
                      m = m + 1
                      global(i,j,k) = cremote(m)
                   end do
                end do
             end do
          end do
          call mpp_sync_self(domain%y%list(:)%pe)
      else
          nd = size(domain%list)
          do n = 1,nd-1
             lpos = mod(domain%pos+nd-n,nd)
             rpos = mod(domain%pos   +n,nd)
             nwords = domain%list(rpos)%x%compute%size * domain%list(rpos)%y%compute%size * size(local,3)
             call mpp_transmit( clocal, size(clocal), domain%list(lpos)%pe, cremote, nwords, domain%list(rpos)%pe )
             m = 0
             do k = 1,size(global,3)
                do j = domain%list(rpos)%y%compute%begin,domain%list(rpos)%y%compute%end
                   do i = domain%list(rpos)%x%compute%begin,domain%list(rpos)%x%compute%end
                      m = m + 1
                      global(i,j,k) = cremote(m)
                   end do
                end do
             end do
          end do
!          call mpp_sync_self(domain%list(:)%pe)
          call mpp_sync_self()
      end if
          
      return
    end subroutine MPP_GLOBAL_FIELD_3D_

    subroutine MPP_GLOBAL_FIELD_4D_( domain, local, global, flags )
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(:,:,:,:)
      MPP_TYPE_, intent(out) :: global(:,:,:,:)
      integer, intent(in), optional :: flags
      MPP_TYPE_ :: local3D (size( local,1),size( local,2),size( local,3)*size(local,4))
      MPP_TYPE_ :: global3D(size(global,1),size(global,2),size(global,3)*size(local,4))
#ifdef use_CRI_pointers
      pointer( lptr, local3D  )
      pointer( gptr, global3D )
      lptr = LOC(local)
      gptr = LOC(global)
      call mpp_global_field( domain, local3D, global3D, flags )
#else
      local3D = RESHAPE( local, SHAPE(local3D) )
      call mpp_global_field( domain, local3D, global3D, flags )
      global = RESHAPE( global3D, SHAPE(global) )
#endif
    end subroutine MPP_GLOBAL_FIELD_4D_

    subroutine MPP_GLOBAL_FIELD_5D_( domain, local, global, flags )
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(:,:,:,:,:)
      MPP_TYPE_, intent(out) :: global(:,:,:,:,:)
      integer, intent(in), optional :: flags
      MPP_TYPE_ :: local3D (size( local,1),size( local,2),size( local,3)*size( local,4)*size(local,5))
      MPP_TYPE_ :: global3D(size(global,1),size(global,2),size(global,3)*size(global,4)*size(local,5))
#ifdef use_CRI_pointers
      pointer( lptr, local3D  )
      pointer( gptr, global3D )
      lptr = LOC(local)
      gptr = LOC(global)
      call mpp_global_field( domain, local3D, global3D, flags )
#else
      local3D = RESHAPE( local, SHAPE(local3D) )
      call mpp_global_field( domain, local3D, global3D, flags )
      global = RESHAPE( global3D, SHAPE(global) )
#endif
    end subroutine MPP_GLOBAL_FIELD_5D_

    subroutine MPP_GLOBAL1D_FIELD_2D_( domain, local, global )
      type(domain1D), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(domain%data%begin:  ,:)
      MPP_TYPE_, intent(out) :: global(domain%global%begin:,:)
      integer :: i, n, nwords, lpos, rpos
      MPP_TYPE_, allocatable, dimension(:) :: clocal, cremote

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GLOBAL_FIELD: must first call mpp_domains_init.' )

      if( size(local,1).NE.domain%data%size   .OR. size(global,1).NE.domain%global%size .OR. &
          size(local,2).NE.size(global,2) ) &
           call mpp_error( FATAL, 'MPP_GLOBAL_FIELD: incoming arrays do not match domain.' )

      allocate( clocal (domain%compute%size    *size(local,2)) )
      allocate( cremote(domain%compute%max_size*size(local,2)) )
      clocal = TRANSFER( local(domain%compute%begin:domain%compute%end,:), clocal )
!always fill local directly
      global(domain%compute%begin:domain%compute%end,:) = local(domain%compute%begin:domain%compute%end,:)
!fill off-domains (note loops begin at an offset of 1)
      n = size(domain%list)
      do i = 1,n-1
         lpos = mod(domain%pos+n-i,n)
         rpos = mod(domain%pos  +i,n)
         nwords = domain%list(rpos)%compute%size * size(local,2)
         call mpp_transmit( clocal, size(clocal), domain%list(lpos)%pe, cremote, nwords, domain%list(rpos)%pe )
         global(domain%list(rpos)%compute%begin:domain%list(rpos)%compute%end,:) = &
              RESHAPE( cremote(1:nwords), (/domain%list(rpos)%compute%size,size(local,2)/) )
      end do
      call mpp_sync_self(domain%list(:)%pe)
!PV786667: the deallocate stmts can be removed when fixed (7.3.1.3m)
      deallocate( clocal, cremote )
          
      return
    end subroutine MPP_GLOBAL1D_FIELD_2D_
