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
    subroutine MPP_UPDATE_DOMAINS_2D_( field, domain, flags )
!updates data domain of 2D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:)
      type(domain2D), intent(in) :: domain
      integer, intent(in), optional :: flags
      MPP_TYPE_ :: field3D(size(field,1),size(field,2),1)
#ifdef use_CRI_pointers
      pointer( ptr, field3D )
      ptr = LOC(field)
      call mpp_update_domains( field3D, domain, flags )
#else
      field3D = RESHAPE( field, SHAPE(field3D) )
      call mpp_update_domains( field3D, domain, flags )
      field = RESHAPE( field3D, SHAPE(field) )
#endif
      return
    end subroutine MPP_UPDATE_DOMAINS_2D_

    subroutine MPP_UPDATE_DOMAINS_4D_( field, domain, flags )
!updates data domain of 4D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:,:,:)
      type(domain2D), intent(in) :: domain
      integer, intent(in), optional :: flags
      MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4))
#ifdef use_CRI_pointers
      pointer( ptr, field3D )
      ptr = LOC(field)
      call mpp_update_domains( field3D, domain, flags )
#else
      field3D = RESHAPE( field, SHAPE(field3D) )
      call mpp_update_domains( field3D, domain, flags )
      field = RESHAPE( field3D, SHAPE(field) )
#endif
      return
    end subroutine MPP_UPDATE_DOMAINS_4D_

    subroutine MPP_UPDATE_DOMAINS_5D_( field, domain, flags )
!updates data domain of 5D field whose computational domains have been computed
      MPP_TYPE_, intent(inout) :: field(:,:,:,:,:)
      type(domain2D), intent(in) :: domain
      integer, intent(in), optional :: flags
      MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4)*size(field,5))
#ifdef use_CRI_pointers
      pointer( ptr, field3D )
      ptr = LOC(field)
      call mpp_update_domains( field3D, domain, flags )
#else
      field3D = RESHAPE( field, SHAPE(field3D) )
      call mpp_update_domains( field3D, domain, flags )
      field = RESHAPE( field3D, SHAPE(field) )
#endif
      return
    end subroutine MPP_UPDATE_DOMAINS_5D_

    subroutine MPP_UPDATE_DOMAINS_3D_( field, domain, flags )
!updates data domain of 3D field whose computational domains have been computed
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(inout) :: field(domain%x%data%begin:,domain%y%data%begin:,:)
      integer, intent(in), optional :: flags
      integer :: update_flags
!equate to mpp_domains_stack
      integer :: wordlen        !#words of MPP_TYPE_ fit in 1 word of mpp_domains_stack
      MPP_TYPE_ :: buffer(size(mpp_domains_stack))
#ifdef use_CRI_pointers
      pointer( ptr, buffer )
#endif
      integer :: buffer_pos
      integer :: i, j, k, m, n
      integer :: is, ie, js, je, ke
!receive domains saved here for unpacking
!for non-blocking version, could be recomputed
      integer, dimension(8) :: isr, ier, jsr, jer
      integer :: to_pe, from_pe, list, pos, msgsize
      logical :: recv_e, recv_se, recv_s, recv_sw, recv_w, recv_nw, recv_n, recv_ne
      logical :: send_e, send_se, send_s, send_sw, send_w, send_nw, send_n, send_ne
      logical :: folded
      character(len=8) :: text

      update_flags = XUPDATE+YUPDATE   !default
      if( PRESENT(flags) )update_flags = flags
      recv_w = BTEST(update_flags,WEST)
      recv_e = BTEST(update_flags,EAST)
      recv_s = BTEST(update_flags,SOUTH)
      recv_n = BTEST(update_flags,NORTH)
      recv_ne = recv_e .AND. recv_n
      recv_se = recv_e .AND. recv_s
      recv_sw = recv_w .AND. recv_s
      recv_nw = recv_w .AND. recv_n
      send_w = recv_e
      send_e = recv_w
      send_s = recv_n
      send_n = recv_s
      send_ne = send_e .AND. send_n
      send_se = send_e .AND. send_s
      send_sw = send_w .AND. send_s
      send_nw = send_w .AND. send_n
      if( recv_w .AND. BTEST(domain%fold,WEST)  .AND. BTEST(grid_offset_type,EAST)  ) &
           call mpp_error( FATAL, 'Incompatible grid offset and fold.' )
      if( recv_s .AND. BTEST(domain%fold,SOUTH) .AND. BTEST(grid_offset_type,NORTH) ) &
           call mpp_error( FATAL, 'Incompatible grid offset and fold.' )
      if( recv_e .AND. BTEST(domain%fold,EAST)  .AND. BTEST(grid_offset_type,WEST)  ) &
           call mpp_error( FATAL, 'Incompatible grid offset and fold.' )
      if( recv_n .AND. BTEST(domain%fold,NORTH) .AND. BTEST(grid_offset_type,SOUTH) ) &
           call mpp_error( FATAL, 'Incompatible grid offset and fold.' )


      n = size(domain%list)
      ke = size(field,3)
      buffer_pos = 0        !this initialization goes away if update_domains becomes non-blocking
#ifdef use_CRI_pointers
      ptr = LOC(mpp_domains_stack)
      wordlen = size(transfer(buffer(1),mpp_domains_stack))
#endif
!send
      do list = 0,n-1
         m = mod( domain%pos+list, n )
         if( .NOT.domain%list(m)%overlap )cycle
         call mpp_clock_begin(pack_clock)
         to_pe = domain%list(m)%pe
         pos = buffer_pos
         if( send_w .AND. domain%list(m)%send_w%overlap )then
             is = domain%list(m)%send_w%is; ie = domain%list(m)%send_w%ie
             js = domain%list(m)%send_w%js; je = domain%list(m)%send_w%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_w_off%is; ie = domain%list(m)%send_w_off%ie
                 js = domain%list(m)%send_w_off%js; je = domain%list(m)%send_w_off%je
             end if
             call mpp_clock_begin(pack_loop_clock)
             do k = 1,ke
                do j = js,je
                   do i = is,ie
                      pos = pos + 1
                      buffer(pos) = field(i,j,k)
                   end do
                end do
             end do
             call mpp_clock_end(pack_loop_clock)
         end if
         if( send_nw .AND. domain%list(m)%send_nw%overlap )then
             is = domain%list(m)%send_nw%is; ie = domain%list(m)%send_nw%ie
             js = domain%list(m)%send_nw%js; je = domain%list(m)%send_nw%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_nw_off%is; ie = domain%list(m)%send_nw_off%ie
                 js = domain%list(m)%send_nw_off%js; je = domain%list(m)%send_nw_off%je
             end if
             call mpp_clock_begin(pack_loop_clock)
             do k = 1,ke
                do j = js,je
                   do i = is,ie
                      pos = pos + 1
                      buffer(pos) = field(i,j,k)
                   end do
                end do
             end do
             call mpp_clock_end(pack_loop_clock)
         end if
         if( send_n .AND. domain%list(m)%send_n%overlap )then
             is = domain%list(m)%send_n%is; ie = domain%list(m)%send_n%ie
             js = domain%list(m)%send_n%js; je = domain%list(m)%send_n%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_n_off%is; ie = domain%list(m)%send_n_off%ie
                 js = domain%list(m)%send_n_off%js; je = domain%list(m)%send_n_off%je
             end if
             call mpp_clock_begin(pack_loop_clock)
             do k = 1,ke
                do j = js,je
                   do i = is,ie
                      pos = pos + 1
                      buffer(pos) = field(i,j,k)
                   end do
                end do
             end do
             call mpp_clock_end(pack_loop_clock)
         end if
         if( send_ne .AND. domain%list(m)%send_ne%overlap )then
             is = domain%list(m)%send_ne%is; ie = domain%list(m)%send_ne%ie
             js = domain%list(m)%send_ne%js; je = domain%list(m)%send_ne%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_ne_off%is; ie = domain%list(m)%send_ne_off%ie
                 js = domain%list(m)%send_ne_off%js; je = domain%list(m)%send_ne_off%je
             end if
             call mpp_clock_begin(pack_loop_clock)
             do k = 1,ke
                do j = js,je
                   do i = is,ie
                      pos = pos + 1
                      buffer(pos) = field(i,j,k)
                   end do
                end do
             end do
             call mpp_clock_end(pack_loop_clock)
         end if
         if( send_e .AND. domain%list(m)%send_e%overlap )then
             is = domain%list(m)%send_e%is; ie = domain%list(m)%send_e%ie
             js = domain%list(m)%send_e%js; je = domain%list(m)%send_e%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_e_off%is; ie = domain%list(m)%send_e_off%ie
                 js = domain%list(m)%send_e_off%js; je = domain%list(m)%send_e_off%je
             end if
             call mpp_clock_begin(pack_loop_clock)
             do k = 1,ke
                do j = js,je
                   do i = is,ie
                      pos = pos + 1
                      buffer(pos) = field(i,j,k)
                   end do
                end do
             end do
             call mpp_clock_end(pack_loop_clock)
         end if
         if( send_se .AND. domain%list(m)%send_se%overlap )then
             is = domain%list(m)%send_se%is; ie = domain%list(m)%send_se%ie
             js = domain%list(m)%send_se%js; je = domain%list(m)%send_se%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_se_off%is; ie = domain%list(m)%send_se_off%ie
                 js = domain%list(m)%send_se_off%js; je = domain%list(m)%send_se_off%je
             end if
             call mpp_clock_begin(pack_loop_clock)
             do k = 1,ke
                do j = js,je
                   do i = is,ie
                      pos = pos + 1
                      buffer(pos) = field(i,j,k)
                   end do
                end do
             end do
             call mpp_clock_end(pack_loop_clock)
         end if
         if( send_s .AND. domain%list(m)%send_s%overlap )then
             is = domain%list(m)%send_s%is; ie = domain%list(m)%send_s%ie
             js = domain%list(m)%send_s%js; je = domain%list(m)%send_s%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_s_off%is; ie = domain%list(m)%send_s_off%ie
                 js = domain%list(m)%send_s_off%js; je = domain%list(m)%send_s_off%je
             end if
             call mpp_clock_begin(pack_loop_clock)
             do k = 1,ke
                do j = js,je
                   do i = is,ie
                      pos = pos + 1
                      buffer(pos) = field(i,j,k)
                   end do
                end do
             end do
             call mpp_clock_end(pack_loop_clock)
         end if
         if( send_sw .AND. domain%list(m)%send_sw%overlap )then
             is = domain%list(m)%send_sw%is; ie = domain%list(m)%send_sw%ie
             js = domain%list(m)%send_sw%js; je = domain%list(m)%send_sw%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%send_sw_off%is; ie = domain%list(m)%send_sw_off%ie
                 js = domain%list(m)%send_sw_off%js; je = domain%list(m)%send_sw_off%je
             end if
             call mpp_clock_begin(pack_loop_clock)
             do k = 1,ke
                do j = js,je
                   do i = is,ie
                      pos = pos + 1
                      buffer(pos) = field(i,j,k)
                   end do
                end do
             end do
             call mpp_clock_end(pack_loop_clock)
         end if
         call mpp_clock_end(pack_clock)
         call mpp_clock_begin(send_clock)
         msgsize = pos - buffer_pos
         if( msgsize.GT.0 )then
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos*wordlen )
             if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                 write( text,'(i8)' )mpp_domains_stack_hwm
                 call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                      //trim(text)//') from all PEs.' )
             end if
             call mpp_send( buffer(buffer_pos+1:buffer_pos+msgsize), msgsize, to_pe )
             buffer_pos = pos
         end if
         call mpp_clock_end(send_clock)
      end do
             
!recv
      do list = 0,n-1
         m = mod( domain%pos+n-list, n )
         if( .NOT.domain%list(m)%overlap )cycle
         call mpp_clock_begin(recv_clock)
         from_pe = domain%list(m)%pe
         msgsize = 0
         if( recv_e .AND. domain%list(m)%recv_e%overlap )then
             is = domain%list(m)%recv_e%is; ie = domain%list(m)%recv_e%ie
             js = domain%list(m)%recv_e%js; je = domain%list(m)%recv_e%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_e_off%is; ie = domain%list(m)%recv_e_off%ie
                 js = domain%list(m)%recv_e_off%js; je = domain%list(m)%recv_e_off%je
             end if
             msgsize = msgsize + (ie-is+1)*(je-js+1)*ke
         end if
         if( recv_se .AND. domain%list(m)%recv_se%overlap )then
             is = domain%list(m)%recv_se%is; ie = domain%list(m)%recv_se%ie
             js = domain%list(m)%recv_se%js; je = domain%list(m)%recv_se%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_se_off%is; ie = domain%list(m)%recv_se_off%ie
                 js = domain%list(m)%recv_se_off%js; je = domain%list(m)%recv_se_off%je
             end if
             msgsize = msgsize + (ie-is+1)*(je-js+1)*ke
         end if
         if( recv_s .AND. domain%list(m)%recv_s%overlap )then
             is = domain%list(m)%recv_s%is; ie = domain%list(m)%recv_s%ie
             js = domain%list(m)%recv_s%js; je = domain%list(m)%recv_s%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_s_off%is; ie = domain%list(m)%recv_s_off%ie
                 js = domain%list(m)%recv_s_off%js; je = domain%list(m)%recv_s_off%je
             end if
             msgsize = msgsize + (ie-is+1)*(je-js+1)*ke
         end if
         if( recv_sw .AND. domain%list(m)%recv_sw%overlap )then
             is = domain%list(m)%recv_sw%is; ie = domain%list(m)%recv_sw%ie
             js = domain%list(m)%recv_sw%js; je = domain%list(m)%recv_sw%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_sw_off%is; ie = domain%list(m)%recv_sw_off%ie
                 js = domain%list(m)%recv_sw_off%js; je = domain%list(m)%recv_sw_off%je
             end if
             msgsize = msgsize + (ie-is+1)*(je-js+1)*ke
         end if
         if( recv_w .AND. domain%list(m)%recv_w%overlap )then
             is = domain%list(m)%recv_w%is; ie = domain%list(m)%recv_w%ie
             js = domain%list(m)%recv_w%js; je = domain%list(m)%recv_w%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_w_off%is; ie = domain%list(m)%recv_w_off%ie
                 js = domain%list(m)%recv_w_off%js; je = domain%list(m)%recv_w_off%je
             end if
             msgsize = msgsize + (ie-is+1)*(je-js+1)*ke
         end if
         if( recv_nw .AND. domain%list(m)%recv_nw%overlap )then
             is = domain%list(m)%recv_nw%is; ie = domain%list(m)%recv_nw%ie
             js = domain%list(m)%recv_nw%js; je = domain%list(m)%recv_nw%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_nw_off%is; ie = domain%list(m)%recv_nw_off%ie
                 js = domain%list(m)%recv_nw_off%js; je = domain%list(m)%recv_nw_off%je
             end if
             msgsize = msgsize + (ie-is+1)*(je-js+1)*ke
         end if
         if( recv_n .AND. domain%list(m)%recv_n%overlap )then
             is = domain%list(m)%recv_n%is; ie = domain%list(m)%recv_n%ie
             js = domain%list(m)%recv_n%js; je = domain%list(m)%recv_n%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_n_off%is; ie = domain%list(m)%recv_n_off%ie
                 js = domain%list(m)%recv_n_off%js; je = domain%list(m)%recv_n_off%je
             end if
             msgsize = msgsize + (ie-is+1)*(je-js+1)*ke
         end if
         if( recv_ne .AND. domain%list(m)%recv_ne%overlap )then
             is = domain%list(m)%recv_ne%is; ie = domain%list(m)%recv_ne%ie
             js = domain%list(m)%recv_ne%js; je = domain%list(m)%recv_ne%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_ne_off%is; ie = domain%list(m)%recv_ne_off%ie
                 js = domain%list(m)%recv_ne_off%js; je = domain%list(m)%recv_ne_off%je
             end if
             msgsize = msgsize + (ie-is+1)*(je-js+1)*ke
         end if
         if( msgsize.GT.0 )then
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize)*wordlen )
             if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                 write( text,'(i8)' )mpp_domains_stack_hwm
                 call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                      //trim(text)//') from all PEs.' )
             end if
             call mpp_recv( buffer(buffer_pos+1:buffer_pos+msgsize), msgsize, from_pe )
             buffer_pos = buffer_pos + msgsize
         end if
         call mpp_clock_end(recv_clock)
      end do
             
!unpack recv
!unpack halos in reverse order
      do list = n-1,0,-1
         m = mod( domain%pos+n-list, n )
         if( .NOT.domain%list(m)%overlap )cycle
         call mpp_clock_begin(unpk_clock)
         from_pe = domain%list(m)%pe
         pos = buffer_pos
         if( recv_ne .AND. domain%list(m)%recv_ne%overlap )then
             is = domain%list(m)%recv_ne%is; ie = domain%list(m)%recv_ne%ie
             js = domain%list(m)%recv_ne%js; je = domain%list(m)%recv_ne%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_ne_off%is; ie = domain%list(m)%recv_ne_off%ie
                 js = domain%list(m)%recv_ne_off%js; je = domain%list(m)%recv_ne_off%je
             end if
             msgsize = (ie-is+1)*(je-js+1)*ke
             pos = buffer_pos - msgsize
             buffer_pos = pos
             if( domain%list(m)%recv_ne%folded )then
                 do k = 1,ke
                    do j = je,js,-1
                       do i = ie,is,-1
                          pos = pos + 1
                          field(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
             else
                 do k = 1,ke
                    do j = js,je
                       do i = is,ie
                          pos = pos + 1
                          field(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
             end if
         end if
         if( recv_n .AND. domain%list(m)%recv_n%overlap )then
             is = domain%list(m)%recv_n%is; ie = domain%list(m)%recv_n%ie
             js = domain%list(m)%recv_n%js; je = domain%list(m)%recv_n%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_n_off%is; ie = domain%list(m)%recv_n_off%ie
                 js = domain%list(m)%recv_n_off%js; je = domain%list(m)%recv_n_off%je
             end if
             msgsize = (ie-is+1)*(je-js+1)*ke
             pos = buffer_pos - msgsize
             buffer_pos = pos
             if( domain%list(m)%recv_n%folded )then
                 do k = 1,ke
                    do j = je,js,-1
                       do i = ie,is,-1
                          pos = pos + 1
                          field(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
             else
                 do k = 1,ke
                    do j = js,je
                       do i = is,ie
                          pos = pos + 1
                          field(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
             end if
         end if
         if( recv_nw .AND. domain%list(m)%recv_nw%overlap )then
             is = domain%list(m)%recv_nw%is; ie = domain%list(m)%recv_nw%ie
             js = domain%list(m)%recv_nw%js; je = domain%list(m)%recv_nw%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_nw_off%is; ie = domain%list(m)%recv_nw_off%ie
                 js = domain%list(m)%recv_nw_off%js; je = domain%list(m)%recv_nw_off%je
             end if
             msgsize = (ie-is+1)*(je-js+1)*ke
             pos = buffer_pos - msgsize
             buffer_pos = pos
             if( domain%list(m)%recv_nw%folded )then
                 do k = 1,ke
                    do j = je,js,-1
                       do i = ie,is,-1
                          pos = pos + 1
                          field(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
             else
                 do k = 1,ke
                    do j = js,je
                       do i = is,ie
                          pos = pos + 1
                          field(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
             end if
         end if
         if( recv_w .AND. domain%list(m)%recv_w%overlap )then
             is = domain%list(m)%recv_w%is; ie = domain%list(m)%recv_w%ie
             js = domain%list(m)%recv_w%js; je = domain%list(m)%recv_w%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_w_off%is; ie = domain%list(m)%recv_w_off%ie
                 js = domain%list(m)%recv_w_off%js; je = domain%list(m)%recv_w_off%je
             end if
             msgsize = (ie-is+1)*(je-js+1)*ke
             pos = buffer_pos - msgsize
             buffer_pos = pos
             if( domain%list(m)%recv_w%folded )then
                 do k = 1,ke
                    do j = je,js,-1
                       do i = ie,is,-1
                          pos = pos + 1
                          field(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
             else
                 do k = 1,ke
                    do j = js,je
                       do i = is,ie
                          pos = pos + 1
                          field(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
             end if
         end if
         if( recv_sw .AND. domain%list(m)%recv_sw%overlap )then
             is = domain%list(m)%recv_sw%is; ie = domain%list(m)%recv_sw%ie
             js = domain%list(m)%recv_sw%js; je = domain%list(m)%recv_sw%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_sw_off%is; ie = domain%list(m)%recv_sw_off%ie
                 js = domain%list(m)%recv_sw_off%js; je = domain%list(m)%recv_sw_off%je
             end if
             msgsize = (ie-is+1)*(je-js+1)*ke
             pos = buffer_pos - msgsize
             buffer_pos = pos
             if( domain%list(m)%recv_sw%folded )then
                 do k = 1,ke
                    do j = je,js,-1
                       do i = ie,is,-1
                          pos = pos + 1
                          field(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
             else
                 do k = 1,ke
                    do j = js,je
                       do i = is,ie
                          pos = pos + 1
                          field(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
             end if
         end if
         if( recv_s .AND. domain%list(m)%recv_s%overlap )then
             is = domain%list(m)%recv_s%is; ie = domain%list(m)%recv_s%ie
             js = domain%list(m)%recv_s%js; je = domain%list(m)%recv_s%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_s_off%is; ie = domain%list(m)%recv_s_off%ie
                 js = domain%list(m)%recv_s_off%js; je = domain%list(m)%recv_s_off%je
             end if
             msgsize = (ie-is+1)*(je-js+1)*ke
             pos = buffer_pos - msgsize
             buffer_pos = pos
             if( domain%list(m)%recv_s%folded )then
                 do k = 1,ke
                    do j = je,js,-1
                       do i = ie,is,-1
                          pos = pos + 1
                          field(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
             else
                 do k = 1,ke
                    do j = js,je
                       do i = is,ie
                          pos = pos + 1
                          field(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
             end if
         end if
         if( recv_se .AND. domain%list(m)%recv_se%overlap )then
             is = domain%list(m)%recv_se%is; ie = domain%list(m)%recv_se%ie
             js = domain%list(m)%recv_se%js; je = domain%list(m)%recv_se%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_se_off%is; ie = domain%list(m)%recv_se_off%ie
                 js = domain%list(m)%recv_se_off%js; je = domain%list(m)%recv_se_off%je
             end if
             msgsize = (ie-is+1)*(je-js+1)*ke
             pos = buffer_pos - msgsize
             buffer_pos = pos
             if( domain%list(m)%recv_se%folded )then
                 do k = 1,ke
                    do j = je,js,-1
                       do i = ie,is,-1
                          pos = pos + 1
                          field(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
             else
                 do k = 1,ke
                    do j = js,je
                       do i = is,ie
                          pos = pos + 1
                          field(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
             end if
         end if
         if( recv_e .AND. domain%list(m)%recv_e%overlap )then
             is = domain%list(m)%recv_e%is; ie = domain%list(m)%recv_e%ie
             js = domain%list(m)%recv_e%js; je = domain%list(m)%recv_e%je
             if( grid_offset_type.NE.AGRID )then
                 is = domain%list(m)%recv_e_off%is; ie = domain%list(m)%recv_e_off%ie
                 js = domain%list(m)%recv_e_off%js; je = domain%list(m)%recv_e_off%je
             end if
             msgsize = (ie-is+1)*(je-js+1)*ke
             pos = buffer_pos - msgsize
             buffer_pos = pos
             if( domain%list(m)%recv_e%folded )then
                 do k = 1,ke
                    do j = je,js,-1
                       do i = ie,is,-1
                          pos = pos + 1
                          field(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
             else
                 do k = 1,ke
                    do j = js,je
                       do i = is,ie
                          pos = pos + 1
                          field(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
             end if
         end if
         call mpp_clock_end(unpk_clock)
      end do

      call mpp_clock_begin(wait_clock)
      call mpp_sync_self( domain%list(:)%pe )
      call mpp_clock_end(wait_clock)
      return
    end subroutine MPP_UPDATE_DOMAINS_3D_

    subroutine MPP_REDISTRIBUTE_2D_( domain_in, field_in, domain_out, field_out )
      type(domain2D), intent(in) :: domain_in, domain_out
      MPP_TYPE_, intent(in)  :: field_in (:,:)
      MPP_TYPE_, intent(out) :: field_out(:,:)
      MPP_TYPE_ :: field3D_in (size(field_in, 1),size(field_in, 2),1)
      MPP_TYPE_ :: field3D_out(size(field_out,1),size(field_out,2),1)
#ifdef use_CRI_pointers
      pointer( ptr_in,  field3D_in  )
      pointer( ptr_out, field3D_out )
      ptr_in  = LOC(field_in )
      ptr_out = LOC(field_out)
      call mpp_redistribute( domain_in, field3D_in, domain_out, field3D_out )
#else
      field3D_in = RESHAPE( field_in, SHAPE(field3D_in) )
      call mpp_redistribute( domain_in, field3D_in, domain_out, field3D_out )
      field_out = RESHAPE( field3D_out, SHAPE(field_out) )
#endif
      return
    end subroutine MPP_REDISTRIBUTE_2D_

    subroutine MPP_REDISTRIBUTE_3D_( domain_in, field_in, domain_out, field_out )
      type(domain2D), intent(in) :: domain_in, domain_out
!      MPP_TYPE_, intent(in)  :: field_in ( domain_in%x%data%begin:, domain_in%y%data%begin:,:)
!      MPP_TYPE_, intent(out) :: field_out(domain_out%x%data%begin:,domain_out%y%data%begin:,:)
      MPP_TYPE_, intent(in)  :: field_in (:,:,:)
      MPP_TYPE_, intent(out) :: field_out(:,:,:)
      integer :: is, ie, js, je, ke, isc, iec, jsc, jec
      integer :: i, j, k
      integer :: list, m, n, pos, msgsize
      integer :: to_pe, from_pe
!      MPP_TYPE_, dimension(domain_in%x%compute%size*domain_in%y%compute%size*size(field_in,3)) :: send_buf, recv_buf
      MPP_TYPE_ :: buffer(size(mpp_domains_stack))
#ifdef use_CRI_pointers
      pointer( ptr, buffer )
#endif
      integer :: buffer_pos, wordlen
      character(len=8) :: text

!      ke = size(field_in,3)
!      if( ke.NE.size(field_out,3) )call mpp_error( FATAL, 'MPP_REDISTRIBUTE: mismatch between field_in and field_out.' )
!      if( UBOUND(field_in,1).NE.domain_in%x%data%end .OR. &
!          UBOUND(field_in,2).NE.domain_in%y%data%end ) &
!          call mpp_error( FATAL, 'MPP_REDISTRIBUTE: field_in must be on data domain of domain_in.' )
!      if( UBOUND(field_out,1).NE.domain_out%x%data%end .OR. &
!          UBOUND(field_out,2).NE.domain_out%y%data%end ) &
!          call mpp_error( FATAL, 'MPP_REDISTRIBUTE: field_out must be on data domain of domain_out.' )

!fix ke
      ke = 0
      if( domain_in%pe.NE.NULL_PE )ke = size(field_in,3)
      if( domain_out%pe.NE.NULL_PE )then
          if( ke.NE.0 .AND. ke.NE.size(field_out,3) ) &
               call mpp_error( FATAL, 'MPP_REDISTRIBUTE: mismatch between field_in and field_out.' )
          ke = size(field_out,3)
      end if
      if( ke.EQ.0 )call mpp_error( FATAL, 'MPP_REDISTRIBUTE: either domain_in or domain_out must be native.' )
!check sizes
      if( domain_in%pe.NE.NULL_PE )then
          if( size(field_in,1).NE.domain_in%x%data%size .OR. size(field_in,2).NE.domain_in%y%data%size ) &
               call mpp_error( FATAL, 'MPP_REDISTRIBUTE: field_in must be on data domain of domain_in.' )
      end if
      if( domain_out%pe.NE.NULL_PE )then
          if( size(field_out,1).NE.domain_out%x%data%size .OR. size(field_out,2).NE.domain_out%y%data%size ) &
               call mpp_error( FATAL, 'MPP_REDISTRIBUTE: field_out must be on data domain of domain_out.' )
      end if

      buffer_pos = 0
#ifdef use_CRI_pointers
      ptr = LOC(mpp_domains_stack)
      wordlen = size(TRANSFER(buffer(1),mpp_domains_stack))
#endif
!send
      call mpp_get_compute_domain( domain_in, isc, iec, jsc, jec )
      n = size(domain_out%list)
      do list = 0,n-1
         m = mod( domain_out%pos+list+n, n )
         to_pe = domain_out%list(m)%pe
         call mpp_get_compute_domain( domain_out%list(m), is, ie, js, je )
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             msgsize = (ie-is+1)*(je-js+1)*ke
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize)*wordlen )
             if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                 write( text,'(i8)' )mpp_domains_stack_hwm
                 call mpp_error( FATAL, 'MPP_REDISTRIBUTE: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                      //trim(text)//') from all PEs.' )
             end if
             pos = buffer_pos
             do k = 1,ke
                do j = js-domain_in%y%data%begin+1,je-domain_in%y%data%begin+1
                   do i = is-domain_in%x%data%begin+1,ie-domain_in%x%data%begin+1
                      pos = pos+1
                      buffer(pos) = field_in(i,j,k)
                   end do
                end do
             end do
             if( debug )write( stderr(),* )'PE', pe, ' to PE ', to_pe, 'is,ie,js,je=', is, ie, js, je
             call mpp_send( buffer(buffer_pos+1:buffer_pos+msgsize), msgsize, to_pe )
             buffer_pos = pos
         end if
      end do
!recv
      call mpp_get_compute_domain( domain_out, isc, iec, jsc, jec )
      n = size(domain_in%list)
      do list = 0,n-1
         m = mod( domain_in%pos+n-list, n )
         from_pe = domain_in%list(m)%pe
         call mpp_get_compute_domain( domain_in%list(m), is, ie, js, je )
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             msgsize = (ie-is+1)*(je-js+1)*ke
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize)*wordlen )
             if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                 write( text,'(i8)' )mpp_domains_stack_hwm
                 call mpp_error( FATAL, 'MPP_REDISTRIBUTE: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                      //trim(text)//') from all PEs.' )
             end if
             if( debug )write( stderr(),* )'PE', pe, ' from PE ', from_pe, 'is,ie,js,je=', is, ie, js, je
             call mpp_recv( buffer(buffer_pos+1:buffer_pos+msgsize), msgsize, from_pe )
             pos = buffer_pos
             do k = 1,ke
                do j = js-domain_out%y%data%begin+1,je-domain_out%y%data%begin+1
                   do i = is-domain_out%x%data%begin+1,ie-domain_out%x%data%begin+1
                      pos = pos+1
                      field_out(i,j,k) = buffer(pos)
                   end do
                end do
             end do
             buffer_pos = pos
         end if
      end do

!      call mpp_sync_self( domain_in%list(:)%pe )
      call mpp_sync_self()
      return
    end subroutine MPP_REDISTRIBUTE_3D_

    subroutine MPP_REDISTRIBUTE_4D_( domain_in, field_in, domain_out, field_out )
      type(domain2D), intent(in) :: domain_in, domain_out
      MPP_TYPE_, intent(in)  :: field_in (:,:,:,:)
      MPP_TYPE_, intent(out) :: field_out(:,:,:,:)
      MPP_TYPE_ :: field3D_in (size(field_in, 1),size(field_in, 2),size(field_in ,3)*size(field_in ,4))
      MPP_TYPE_ :: field3D_out(size(field_out,1),size(field_out,2),size(field_out,3)*size(field_out,4))
#ifdef use_CRI_pointers
      pointer( ptr_in,  field3D_in  )
      pointer( ptr_out, field3D_out )
      ptr_in  = LOC(field_in )
      ptr_out = LOC(field_out)
      call mpp_redistribute( domain_in, field3D_in, domain_out, field3D_out )
#else
      field3D_in = RESHAPE( field_in, SHAPE(field3D_in) )
      call mpp_redistribute( domain_in, field3D_in, domain_out, field3D_out )
      field_out = RESHAPE( field3D_out, SHAPE(field_out) )
#endif
      return
    end subroutine MPP_REDISTRIBUTE_4D_

    subroutine MPP_REDISTRIBUTE_5D_( domain_in, field_in, domain_out, field_out )
      type(domain2D), intent(in) :: domain_in, domain_out
      MPP_TYPE_, intent(in)  :: field_in (:,:,:,:,:)
      MPP_TYPE_, intent(out) :: field_out(:,:,:,:,:)
      MPP_TYPE_ :: field3D_in (size(field_in, 1),size(field_in, 2),size(field_in ,3)*size(field_in ,4)*size(field_in ,5))
      MPP_TYPE_ :: field3D_out(size(field_out,1),size(field_out,2),size(field_out,3)*size(field_out,4)*size(field_out,5))
#ifdef use_CRI_pointers
      pointer( ptr_in,  field3D_in  )
      pointer( ptr_out, field3D_out )
      ptr_in  = LOC(field_in )
      ptr_out = LOC(field_out)
      call mpp_redistribute( domain_in, field3D_in, domain_out, field3D_out )
#else
      field3D_in = RESHAPE( field_in, SHAPE(field3D_in) )
      call mpp_redistribute( domain_in, field3D_in, domain_out, field3D_out )
      field_out = RESHAPE( field3D_out, SHAPE(field_out) )
#endif
      return
    end subroutine MPP_REDISTRIBUTE_5D_
#ifdef VECTOR_FIELD_
!VECTOR_FIELD_ is set to false for MPP_TYPE_ integer or logical.
!vector fields
    subroutine MPP_UPDATE_DOMAINS_2D_V_( fieldx, fieldy, domain, flags, gridtype )
!updates data domain of 2D field whose computational domains have been computed
      MPP_TYPE_, intent(inout), dimension(:,:) :: fieldx, fieldy
      type(domain2D), intent(inout) :: domain
      integer, intent(in), optional :: flags, gridtype
      MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),1)
      MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),1)
#ifdef use_CRI_pointers
      pointer( ptrx, field3Dx )
      pointer( ptry, field3Dy )
      ptrx = LOC(fieldx)
      ptry = LOC(fieldy)
      call mpp_update_domains( field3Dx, field3Dy, domain, flags, gridtype )
#else
      field3Dx = RESHAPE( fieldx, SHAPE(field3Dx) )
      field3Dy = RESHAPE( fieldy, SHAPE(field3Dy) )
      call mpp_update_domains( field3Dx, field3Dy, domain, flags, gridtype )
      fieldx = RESHAPE( field3Dx, SHAPE(fieldx) )
      fieldy = RESHAPE( field3Dy, SHAPE(fieldy) )
#endif
      return
    end subroutine MPP_UPDATE_DOMAINS_2D_V_

    subroutine MPP_UPDATE_DOMAINS_3D_V_( fieldx, fieldy, domain, flags, gridtype )
!updates data domain of 3D field whose computational domains have been computed
      type(domain2D), intent(inout) :: domain
      MPP_TYPE_, intent(inout), dimension(domain%x%data%begin:,domain%y%data%begin:,:) :: fieldx, fieldy
      integer, intent(in), optional :: flags, gridtype
      integer :: update_flags
      integer :: i,j,k,n, is, ie, js, je, ke, pos
      integer :: ioff, joff
      MPP_TYPE_ :: buffer(size(mpp_domains_stack))
#ifdef use_CRI_pointers
      pointer( ptr, buffer )
      MPP_TYPE_ :: field(size(fieldx,1),size(fieldx,2),2*size(fieldx,3))
      pointer( ptrf, field )
#endif
      integer :: buffer_pos, msgsize, wordlen
      character(len=8) :: text

!gridtype
      grid_offset_type = AGRID
      if( PRESENT(gridtype) )then
          if( gridtype.NE.AGRID .AND. &
              gridtype.NE.BGRID_NE .AND. gridtype.NE.BGRID_SW .AND. &
              gridtype.NE.CGRID_NE .AND. gridtype.NE.CGRID_SW ) &
               call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: gridtype must be one of AGRID|BGRID_NE|BGRID_SW|CGRID_NE|CGRID_SW.' )
!grid_offset_type used by update domains to determine shifts.
          grid_offset_type = gridtype
          call compute_overlaps(domain)
          if( grid_offset_type.NE.domain%gridtype ) &
               call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: gridtype cannot be changed during run.' )
      end if
!need to add code for EWS boundaries
      if( BTEST(domain%fold,WEST) .AND. BTEST(update_flags,WEST) ) &
           call mpp_error( FATAL, 'velocity stencil not yet active for WEST fold, contact author.' )
      if( BTEST(domain%fold,EAST) .AND. BTEST(update_flags,EAST) ) &
           call mpp_error( FATAL, 'velocity stencil not yet active for EAST fold, contact author.' )
      if( BTEST(domain%fold,SOUTH) .AND. BTEST(update_flags,SOUTH) ) &
           call mpp_error( FATAL, 'velocity stencil not yet active for SOUTH fold, contact author.' )

#ifdef use_CRI_pointers
!optimization for BGRID when fieldx and fieldy are contiguous
      ptrf = LOC(fieldx)
      if( LOC(field(1,1,size(fieldx,3)+1)).EQ.LOC(fieldy) .AND. &
           ( domain%gridtype.EQ.BGRID_NE .OR. domain%gridtype.EQ.BGRID_SW ) )then
          call mpp_update_domains( field, domain, flags )
      else
          call mpp_update_domains( fieldx, domain, flags )
          call mpp_update_domains( fieldy, domain, flags )
      end if
#else
      call mpp_update_domains( fieldx, domain, flags )
      call mpp_update_domains( fieldy, domain, flags )
#endif

#ifdef use_CRI_pointers
      ptr = LOC(mpp_domains_stack)
#endif
      wordlen = size(TRANSFER(buffer(1),mpp_domains_stack))
      buffer_pos = 0
!for all gridtypes
      update_flags = XUPDATE+YUPDATE   !default
      if( PRESENT(flags) )update_flags = flags
      ke = size(fieldx,3)
      call mpp_get_global_domain( domain, xsize=ioff, ysize=joff )
!northern boundary fold
      if( BTEST(domain%fold,NORTH) .AND. BTEST(update_flags,NORTH) )then
          js = domain%y%global%end + 1
          je = domain%y%data%end
          if( je.GE.js )then
!on offset grids, we need to move data leftward by one point
              pos = domain%x%pos - 1 !the one on your left
              if( pos.GE.0 )then
                  is = domain%x%list(pos)%data%end+1; ie=is
              else if( domain%x%cyclic )then
                  pos = pos + size(domain%x%list)
                  is = domain%x%list(pos)%data%end+1 - ioff; ie=is
              else
                  is=1; ie=0
              end if
              n = buffer_pos
              if( ie.EQ.is )then
                  msgsize = (je-js+1)*ke*2 !only half this on CGRID actually
                  mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize)*wordlen )
                  if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                      write( text,'(i8)' )mpp_domains_stack_hwm
                      call mpp_error( FATAL, 'MPP_UPDATE: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                           //trim(text)//') from all PEs.' )
                  end if
                  select case(grid_offset_type)
                  case(BGRID_NE)
                      do k = 1,ke
                         do j = js,je
                            n = n + 2
                            buffer(n-1) = fieldx(is,j,k)
                            buffer(n  ) = fieldy(is,j,k)
                         end do
                      end do
                      call mpp_send( buffer(buffer_pos+1:buffer_pos+n), n, domain%x%list(pos)%pe )
                      buffer_pos = buffer_pos + n
                  case(CGRID_NE)
                      do k = 1,ke
                         do j = js,je
                            n = n + 1
                            buffer(n) = fieldx(is,j,k)
                         end do
                      end do
                      call mpp_send( buffer(buffer_pos+1:buffer_pos+n), n, domain%x%list(pos)%pe )
                      buffer_pos = buffer_pos + n
                  end select
!receive data at x%data%end
                  pos = domain%x%pos + 1 !the one on your right
                  if( pos.LT.size(domain%x%list) )then
                      n = (je-js+1)*ke
                  else if( domain%x%cyclic )then
                      pos = pos - size(domain%x%list)
                      n = (je-js+1)*ke
                  else
                      n = 0
                  end if
                  if( n.GT.0 )then
                      mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+n)*wordlen )
                      if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                          write( text,'(i8)' )mpp_domains_stack_hwm
                          call mpp_error( FATAL, 'MPP_UPDATE: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                               //trim(text)//') from all PEs.' )
                      end if
                      select case(grid_offset_type)
                      case(BGRID_NE)
                          do k = 1,ke
                             do j = js,je
                                do i = domain%x%data%begin,domain%x%data%end-1
                                   fieldx(i,j,k) = fieldx(i+1,j,k)
                                   fieldy(i,j,k) = fieldy(i+1,j,k)
                                end do
                             end do
                          end do
                          n = n*2
                          call mpp_recv( buffer(buffer_pos+1:buffer_pos+n), n, domain%x%list(pos)%pe )
                          i = domain%x%data%end
                          n = buffer_pos
                          do k = 1,ke
                             do j = js,je
                                n = n + 2
                                fieldx(i,j,k) = buffer(n-1)
                                fieldy(i,j,k) = buffer(n  )
                             end do
                          end do
                      case(CGRID_NE)
                          do k = 1,ke
                             do j = js,je
                                do i = domain%x%data%begin,domain%x%data%end-1
                                   fieldx(i,j,k) = fieldx(i+1,j,k)
                                   fieldy(i,j,k) = fieldy(i+1,j,k)
                                end do
                             end do
                          end do
                          call mpp_recv( buffer(buffer_pos+1:buffer_pos+n), n, domain%x%list(pos)%pe )
                          i = domain%x%data%end
                          n = buffer_pos
                          do k = 1,ke
                             do j = js,je
                                n = n + 1
                                fieldx(i,j,k) = buffer(n)
                             end do
                          end do
                      end select
                  end if
              end if
!flip the sign
              is = domain%x%data%begin
              ie = domain%x%data%end
              do k = 1,ke
                 do j = js,je
                    do i = is,ie
                       fieldx(i,j,k) = -fieldx(i,j,k)
                       fieldy(i,j,k) = -fieldy(i,j,k)
                    end do
                 end do
              end do
          end if
!eliminate redundant vector data at fold
          j = domain%y%global%end
          if( domain%y%data%begin.LE.j .AND. j.LE.domain%y%data%end )then !fold is within domain
!ship left-half data to right half: on BGRID_NE the x%data%end point is not in mirror domain and must be done separately.
              if( domain%x%pos.LT.(size(domain%x%list)+1)/2 )then
                  is = domain%x%data%begin
                  ie = min(domain%x%data%end,(domain%x%global%begin+domain%x%global%end)/2)
                  n = buffer_pos
                  select case(grid_offset_type)
                  case(BGRID_NE)
                      do k = 1,ke
                         do i = is,ie-1
                            n = n + 2
                            buffer(n-1) = fieldx(i,j,k)
                            buffer(n)   = fieldy(i,j,k)
                         end do
                      end do
                      call mpp_send( buffer(buffer_pos+1:n), n-buffer_pos, domain%x%list(size(domain%x%list)-domain%x%pos-1)%pe )
                      buffer_pos = n
                  case(CGRID_NE)
                      do k = 1,ke
                         do i = is,ie
                            n = n + 1
                            buffer(n) = fieldy(i,j,k)
                         end do
                      end do
                      call mpp_send( buffer(buffer_pos+1:n), n-buffer_pos, domain%x%list(size(domain%x%list)-domain%x%pos-1)%pe )
                      buffer_pos = n
                  end select
              end if
              if( domain%x%pos.GE.size(domain%x%list)/2 )then
                  is = max(domain%x%data%begin,(domain%x%global%begin+domain%x%global%end)/2+1)
                  ie = domain%x%data%end
                  select case(grid_offset_type)
                  case(BGRID_NE)
                      n = (ie-is+1)*ke*2
                      call mpp_recv( buffer(buffer_pos+1:buffer_pos+n), n, domain%x%list(size(domain%x%list)-domain%x%pos-1)%pe )
                      n = buffer_pos
!get all values except at x%data%end
                      do k = 1,ke
                         do i = ie-1,is,-1
                            n = n + 2
                            fieldx(i,j,k) = -buffer(n-1)
                            fieldy(i,j,k) = -buffer(n)
                         end do
                      end do
!now get the value at domain%x%data%end
                      pos = domain%x%pos - 1
                      if( pos.GE.size(domain%x%list)/2 )then
                          i = domain%x%list(pos)%data%end
                          buffer_pos = n
                          do k = 1,ke
                             n = n + 2
                             buffer(n-1) = fieldx(i,j,k)
                             buffer(n  ) = fieldy(i,j,k)
                          end do
                          call mpp_send( buffer(buffer_pos+1:n), n-buffer_pos, domain%x%list(pos)%pe )
                          buffer_pos = n
                      end if
                      pos = domain%x%pos + 1
                      if( pos.LT.size(domain%x%list) )then
                          n = ke*2
                          call mpp_recv( buffer(buffer_pos+1:buffer_pos+n), n, domain%x%list(pos)%pe )
                          n = buffer_pos
                          i = domain%x%data%end
                          do k = 1,ke
                             n = n + 2
                             fieldx(i,j,k) = buffer(n-1)
                             fieldy(i,j,k) = buffer(n  )
                          end do
                      end if
                  case(CGRID_NE)
                      n = (ie-is+1)*ke
                      call mpp_recv( buffer(buffer_pos+1:buffer_pos+n), n, domain%x%list(size(domain%x%list)-domain%x%pos-1)%pe )
                      n = buffer_pos
                      do k = 1,ke
                         do i = ie,is,-1
                            n = n + 1
                            fieldy(i,j,k) = -buffer(n)
                         end do
                      end do
                  end select
              end if
!poles set to 0: BGRID only
              if( grid_offset_type.EQ.BGRID_NE )then
                  do i = domain%x%global%begin-1,domain%x%global%end,(domain%x%global%begin+domain%x%global%end)/2
                     if( domain%x%data%begin.LE.i .AND. i.LE.domain%x%data%end )then
                         do k = 1,ke
                            fieldx(i,j,k) = 0.
                            fieldy(i,j,k) = 0.
                         end do
                     end if
                  end do
              end if
!these last three code blocks correct an error where the data in your halo coming from other half may have the wrong sign
!off west edge
              select case(grid_offset_type)
              case(BGRID_NE)
                  is = domain%x%global%begin - 1
                  if( is.GT.domain%x%data%begin )then
                      if( 2*is-domain%x%data%begin.GT.domain%x%data%end ) &
                           call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_V: BGRID_NE west edge ubound error.' )
                      do k = 1,ke
                         do i = domain%x%data%begin,is-1
                            fieldx(i,j,k) = fieldx(2*is-i,j,k)
                            fieldy(i,j,k) = fieldy(2*is-i,j,k)
                         end do
                      end do
                  end if
              case(CGRID_NE)
                  is = domain%x%global%begin
                  if( is.GT.domain%x%data%begin )then
                      if( 2*is-domain%x%data%begin-1.GT.domain%x%data%end ) &
                           call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_V: CGRID_NE west edge ubound error.' )
                      do k = 1,ke
                         do i = domain%x%data%begin,is-1
                            fieldy(i,j,k) = fieldy(2*is-i-1,j,k)
                         end do
                      end do
                  end if
              end select
!right of midpoint
              is = (domain%x%global%begin+domain%x%global%end)/2
              if( domain%x%compute%begin.LE.is .AND. is.LT.domain%x%data%end )then
                  select case(grid_offset_type)
                  case(BGRID_NE)
                      ie = domain%x%data%end
                      if( 2*is-ie.LT.domain%x%data%begin )ie = ie - 1
                      if( 2*is-ie.LT.domain%x%data%begin ) &
                           call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_V: BGRID_NE midpoint lbound error.' )
                      do k = 1,ke
                         do i = is+1,ie
                            fieldx(i,j,k) = -fieldx(2*is-i,j,k)
                            fieldy(i,j,k) = -fieldy(2*is-i,j,k)
                         end do
                      end do
                  case(CGRID_NE)
                      if( 2*is-domain%x%data%end+1.LT.domain%x%data%begin ) &
                           call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_V: CGRID_NE midpoint lbound error.' )
                     do k = 1,ke
                         do i = is+1,domain%x%data%end
                            fieldy(i,j,k) = -fieldy(2*is-i+1,j,k)
                         end do
                      end do
                  end select
              end if
!off east edge
              is = domain%x%global%end
              if( is.LT.domain%x%data%end )then
                  select case(grid_offset_type)
                  case(BGRID_NE)
                      if( 2*is-domain%x%data%end.LT.domain%x%data%begin ) &
                           call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_V: BGRID_NE east edge lbound error.' )
                      do k = 1,ke
                         do i = is+1,domain%x%data%end
                            fieldx(i,j,k) = fieldx(2*is-i,j,k)
                            fieldy(i,j,k) = fieldy(2*is-i,j,k)
                         end do
                      end do
                  case(CGRID_NE)
                      if( 2*is-domain%x%data%end+1.LT.domain%x%data%begin ) &
                           call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS_V: CGRID_NE east edge lbound error.' )
                      do k = 1,ke
                         do i = is+1,domain%x%data%end
                            fieldy(i,j,k) = fieldy(2*is-i+1,j,k)
                         end do
                      end do
                  end select
              end if
          end if
      end if

      grid_offset_type = AGRID  !reset
      call mpp_sync_self()
      return
    end subroutine MPP_UPDATE_DOMAINS_3D_V_

    subroutine MPP_UPDATE_DOMAINS_4D_V_( fieldx, fieldy, domain, flags, gridtype )
!updates data domain of 4D field whose computational domains have been computed
      MPP_TYPE_, intent(inout), dimension(:,:,:,:) :: fieldx, fieldy
      type(domain2D), intent(inout) :: domain
      integer, intent(in), optional :: flags, gridtype
      MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),size(fieldx,3)*size(fieldx,4))
      MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),size(fieldy,3)*size(fieldy,4))
#ifdef use_CRI_pointers
      pointer( ptrx, field3Dx )
      pointer( ptry, field3Dy )
      ptrx = LOC(fieldx)
      ptry = LOC(fieldy)
      call mpp_update_domains( field3Dx, field3Dy, domain, flags, gridtype )
#else
      field3Dx = RESHAPE( fieldx, SHAPE(field3Dx) )
      field3Dy = RESHAPE( fieldy, SHAPE(field3Dy) )
      call mpp_update_domains( field3Dx, field3Dy, domain, flags, gridtype )
      fieldx = RESHAPE( field3Dx, SHAPE(fieldx) )
      fieldy = RESHAPE( field3Dy, SHAPE(fieldy) )
#endif
      return
    end subroutine MPP_UPDATE_DOMAINS_4D_V_

    subroutine MPP_UPDATE_DOMAINS_5D_V_( fieldx, fieldy, domain, flags, gridtype )
!updates data domain of 5D field whose computational domains have been computed
      MPP_TYPE_, intent(inout), dimension(:,:,:,:,:) :: fieldx, fieldy
      type(domain2D), intent(inout) :: domain
      integer, intent(in), optional :: flags, gridtype
      MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),size(fieldx,3)*size(fieldx,4)*size(fieldx,5))
      MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),size(fieldy,3)*size(fieldy,4)*size(fieldy,5))
#ifdef use_CRI_pointers
      pointer( ptrx, field3Dx )
      pointer( ptry, field3Dy )
      ptrx = LOC(fieldx)
      ptry = LOC(fieldy)
      call mpp_update_domains( field3Dx, field3Dy, domain, flags, gridtype )
#else
      field3Dx = RESHAPE( fieldx, SHAPE(field3Dx) )
      field3Dy = RESHAPE( fieldy, SHAPE(field3Dy) )
      call mpp_update_domains( field3Dx, field3Dy, domain, flags, gridtype )
      fieldx = RESHAPE( field3Dx, SHAPE(fieldx) )
      fieldy = RESHAPE( field3Dy, SHAPE(fieldy) )
#endif
      return
    end subroutine MPP_UPDATE_DOMAINS_5D_V_
#endif VECTOR_FIELD_
