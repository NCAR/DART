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

   if (.not.module_is_initialized) call fms_io_init ( )

   if (.not.associated(Current_domain)) call mpp_error &
        (FATAL,'set_domain not called')

       if (present(end)) end = .false. 

       if (do_read()) read (unit,end=10) gdata
       if (.not.read_all_pe) then
         ! did not read on all PEs need to broadcast data
           lngth = size(gdata,1)*size(gdata,2)*size(gdata,3)*size(gdata,4)
!           call mpp_transmit ( gdata, lngth, ALL_PES, &
!                               gdata, lngth, 0        )       
           call mpp_broadcast ( gdata, lngth, mpp_root_pe() )
       endif   
     ! return data for compute domain
       data(is:ie,js:je,:,:) = gdata(is:ie,js:je,:,:)
       return  

   10  call read_eof (end)

