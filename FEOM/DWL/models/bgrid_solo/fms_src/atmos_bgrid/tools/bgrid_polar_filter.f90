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

module bgrid_polar_filter_mod
!-----------------------------------------------------------------------

use types_mod, only : r8
use bgrid_horiz_mod, only: horiz_grid_type, bgrid_type
use         fft_mod, only: fft_init, fft_grid_to_fourier,  &
                                     fft_fourier_to_grid
use         fms_mod, only: error_mesg, FATAL, write_version_number


implicit none
private

!-----------------------------------------------------------------------

   public   polar_filter, polar_filter_wind, polar_filter_init

   public pfilt_control_type
   public pfilt_index_type

   type pfilt_index_type
        private
        integer    :: len, len1, leng, lenc
        integer    :: ilb, iub, jlb, jub
        integer    :: is, ie, js, je
        integer    :: jss, jes, jsn, jen
        integer    :: nlpf_s, nlpf_n, nlpf
        integer    :: maxlen, weight
        integer    :: iloc, jloc
        logical    :: sigma
        integer,pointer :: pelist(:), sizelist(:)
        integer,pointer :: nlpf2d(:,:), pelist2d(:,:)
        real(r8)            :: dlm, cph0
        real(r8),   pointer :: cph(:), sklm(:)
   end type pfilt_index_type

   type pfilt_control_type
        private
        integer    :: nlpf, npes
        real(r8), pointer, dimension(:,:) :: slm, clm
        type (pfilt_index_type) :: Tmp, Vel
   end type pfilt_control_type

!    pelist   = PE numbers along x direction
!    sizelist = number of grid points on each PE along x direction
!    pelist2d = PE numbers arranged by the grid's decomposition
!    nlpf2d   = number of latitude rows of polar filter for all PEs
!               arranged by the grid's decomposition

!---------------------- private data -----------------------------------

   real(r8)   , allocatable, dimension(:) :: avg
   integer, allocatable, dimension(:) :: icnt

   logical :: do_load_balance = .true.

!-----------------------------------------------------------------------
!----------------- interfaces ------------------------------------------

   interface polar_filter
     module procedure polar_filter_2d, polar_filter_3d, polar_filter_4d
     module procedure polar_filter_two_3d
   end interface

   interface polar_filter_wind
     module procedure polar_filter_wind_2d, polar_filter_wind_3d
   end interface

!-----------------------------------------------------------------------
!  timing data

 integer, parameter :: time_level = 7
 integer :: id_total
 logical :: do_clock_init = .true.

! version id info

 character(len=128) :: version='$Revision$'
 character(len=128) :: tagname='$Id$'
 logical :: do_log = .true.

!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine filter_field (Index, max_fft, nlpf_zonal,       &
                          put_pe, put_len, get_pe, get_len, &
                          data, mask                        )

!-----------------------------------------------------------------------

 type(pfilt_index_type), intent(in)     :: Index
integer, intent(in)                     :: max_fft
integer, intent(in),    dimension(0:)   :: nlpf_zonal
integer, intent(in),    dimension(:)    :: put_pe, put_len, &
                                           get_pe, get_len

   real(r8), intent(inout), dimension(Index%ilb:,Index%jlb:,:) :: data
   real(r8), intent(in),    dimension(Index%ilb:,Index%jlb:,:), &
                                                     optional :: mask

!-----------------------------------------------------------------------
!      ----------- local storage -------------

    real(r8), dimension(Index%len1,Index%nlpf*size(data,3)) :: g, gm
    real(r8), dimension(Index%leng+1,max_fft) :: z, zm
    real(r8), dimension(Index%lenc  ,max_fft) :: ss
 complex(r8), dimension(Index%lenc  ,max_fft) :: c

 integer :: is, ie, jss, jes, jsn, jen, i, k, n, kd, lngth, lenc,  &
            nlpf, nlpfs, nlpfn, nfft, n0, n1, n2, len1, leng
 real(r8)    :: zm_min

!-----------------------------------------------------------------------
!     ---- return if no rows are filtered ----

      !call mpp_clock_begin (id_total)
      if (max_fft == 0) then
      !    call mpp_clock_end (id_total)
          return
      endif

!       --- check the dimensions of the input array ---

  if (ubound(data,1) /= Index%iub .or.  &
      ubound(data,2) /= Index%jub ) then
!!!      print *, 'pe, size data = ',mpp_pe(),size(data,1),size(data,2)
!!!      call error_mesg ('polar_filter',  &
!!!                       'input array has the wrong dimensions.', FATAL)
  endif

!-----------------------------------------------------------------------
!----------------- setup grid constants --------------------------------

       lngth  = Index%len;  len1 = Index%len1
       leng = Index%leng; lenc = Index%lenc
       is   = Index%is;   ie   = Index%ie

       jss = Index%jss
       jes = Index%jes
       jsn = Index%jsn
       jen = Index%jen
       nlpfs = Index%nlpf_s
       nlpfn = Index%nlpf_n
       nlpf  = Index%nlpf
       kd = size(data,3)

!print *, 'PE,max_fft,nlpfs,nlpfn,nlpf,kd=',mpp_pe(), &
!             max_fft,nlpfs,nlpfn,nlpf,kd
!=======================================================================
!----------- gather contiguous data needed for filtering ---------------

  do k = 1, size(data,3)
     n0 = (k-1)*nlpf
     n1 = n0 + nlpfs
     n2 = n1 + nlpfn

     if (nlpfs > 0) g (1:lngth,n0+1:n1) = data(is:ie,jss:jes:+1,k)
     if (nlpfn > 0) g (1:lngth,n1+1:n2) = data(is:ie,jen:jsn:-1,k)

     if (.not.Index%sigma .and. present(mask)) then
        if (nlpfs > 0) gm(1:lngth,n0+1:n1) = mask(is:ie,jss:jes:+1,k)
        if (nlpfn > 0) gm(1:lngth,n1+1:n2) = mask(is:ie,jen:jsn:-1,k)
     endif

  enddo

!=======================================================================
!----- add cos lat to end of row (needed for filter) -------------------

  if ( nlpf > 0 .and. len1 == lngth+1) then
      do k = 1, size(data,3)
         n0 = (k-1)*nlpf
         n2 = n0 + nlpf
         g (len1,n0+1:n2) = Index%cph (1:nlpf)
         gm(len1,n0+1:n2) = Index%cph (1:nlpf)
      enddo
  endif

!=======================================================================
!-------------------- gather full zonal rows ---------------------------
  nfft = nlpf_zonal(0)

  if ( nlpf > 0 ) then

      call gather_zonal_rows (Index, nlpf_zonal, g, z(:,1:nfft) )

      if (.not.Index%sigma .and. present(mask)) then
        call gather_zonal_rows (Index, nlpf_zonal, gm, zm(:,1:nfft) )
      endif

  endif

!=======================================================================
!---------------- fill in missing values -------------------------------
!---------------- find  minimum value of mask --------------------------

    zm_min = 1.0

    if (.not.Index%sigma .and. present(mask)) then
        zm_min = minval(zm(1:leng,1:nfft))
        if (zm_min < 0.5) call fill_missing (zm(1:leng,1:nfft),  &
                                              z(1:leng,1:nfft))
    endif

!=======================================================================
!------- load balancing of fft rows ---------
    call fft_transmit_to ( z, nfft, put_pe, put_len, get_pe, get_len )

!=======================================================================
!-------- compute filter response ---------

    call set_filter_response ( Index, z(leng+1,1:nfft), ss(:,1:nfft) )

!=======================================================================
!---------------- transform to fourier coefficients --------------------
       c (:,1:nfft) = fft_grid_to_fourier (z(:,1:nfft))

!------------------------- filter --------------------------------------

       c (:,1:nfft) = c (:,1:nfft) * ss (:,1:nfft)

!---------------- transform back to grid spce --------------------------

       z (:,1:nfft) = fft_fourier_to_grid (c(:,1:nfft))

!=======================================================================
!----- return data to original pe ----

    nfft = nlpf_zonal(0)
    call fft_transmit_back ( z, nfft, put_pe, put_len, get_pe, get_len)

!=======================================================================
!--------------- restore zonal mean if mask used -----------------------


    if (zm_min < 0.5) call fix_missing (zm(1:leng,1:nfft),  &
                                         z(1:leng,1:nfft))

!=======================================================================
!------------------- scatter full zonal rows ---------------------------

  if ( nlpf > 0 ) then
      z(leng+1,:) = 0.0
      call scatter_zonal_rows ( Index, nlpf_zonal, z(:,1:nfft), g )
  endif

!=======================================================================
!-------------- return filtered data to original array -----------------
  do k = 1, size(data,3)
     n0 = (k-1)*nlpf
     n1 = n0 + nlpfs
     n2 = n1 + nlpfn

     if (nlpfs > 0) data(is:ie,jss:jes:+1,k) = g(1:lngth,n0+1:n1)
     if (nlpfn > 0) data(is:ie,jen:jsn:-1,k) = g(1:lngth,n1+1:n2)

  enddo

  !call mpp_clock_end (id_total)
!-----------------------------------------------------------------------

 end subroutine filter_field

!#######################################################################

 subroutine filter_two_fields (Index, max_fft, nlpf_zonal,       &
                               put_pe, put_len, get_pe, get_len, &
                               u, v, slm, clm, mask)

!-----------------------------------------------------------------------
 type(pfilt_index_type), intent(in)     :: Index
integer, intent(in)                     :: max_fft
integer, intent(in),    dimension(0:)   :: nlpf_zonal
integer, intent(in),    dimension(:)    :: put_pe, put_len, &
                                           get_pe, get_len

 real(r8), intent(inout), dimension(Index%ilb:,Index%jlb:,:) :: u,v
 real(r8), intent(in),    dimension(:,:), optional   :: slm, clm
 real(r8), intent(in),    dimension(Index%ilb:,Index%jlb:,:), &
                                                    optional :: mask
!-----------------------------------------------------------------------
!      ----------- local storage -------------

     real(r8), dimension(Index%len1,Index%nlpf*size(u,3)*2) :: g, gm
     real(r8), dimension(Index%leng+1,max_fft) :: z, zm
     real(r8), dimension(Index%lenc  ,max_fft) :: ss
  complex(r8), dimension(Index%lenc  ,max_fft) :: c

  integer :: is, ie, jss, jes, jsn, jen, k, kd, lngth, leng, lenc, len1, &
             ns, nn, n0, n1, n2, n3, n4
  integer :: nlpfs, nlpfn, nlpf, nfft
  real(r8)    :: gm_min

!-----------------------------------------------------------------------
!     ---- return if no rows are filtered ----

      !call mpp_clock_begin (id_total)
      if (max_fft == 0) then
          !call mpp_clock_end (id_total)
          return
      endif

!          --- check the dimensions of the input array ---

  if (ubound(u,1) /= Index%iub .or. ubound(u,2) /= Index%jub)  &
           call error_mesg ('polar_filter_wind',  &
                         'input array has the wrong dimensions.', FATAL)

!-----------------------------------------------------------------------
!----------------- setup grid constants --------------------------------

      lngth  = Index%len;   len1 = Index%len1
      leng = Index%leng;  lenc = Index%lenc
      is   = Index%is;    ie   = Index%ie

      jss = Index%jss
      jes = Index%jes
      jsn = Index%jsn
      jen = Index%jen
      nlpfs = Index%nlpf_s
      nlpfn = Index%nlpf_n
      nlpf  = Index%nlpf
      kd    = size(u,3)

      ns = Index%nlpf_s
      nn = Index%nlpf_n

!=======================================================================
!----------- gather contiguous data needed for filtering ---------------

!----------------- gather rows needed for filtering --------------------
!-----------------  southern & northern hemisphere  --------------------

  do k = 1, size(u,3)
     n0 = (k-1)*nlpf*2
     n1 = n0 + nlpfs
     n2 = n1 + nlpfn
     n3 = n2 + nlpfs
     n4 = n3 + nlpfn

     if (present(slm) .and. present(clm)) then
        g (1:lngth,n0+1:n1) = -u(is:ie,jss:jes:+1,k)*slm(:,1:ns) &
                            +v(is:ie,jss:jes:+1,k)*clm(:,1:ns)
        g (1:lngth,n1+1:n2) = -u(is:ie,jen:jsn:-1,k)*slm(:,1:nn) &
                            -v(is:ie,jen:jsn:-1,k)*clm(:,1:nn)
        g (1:lngth,n2+1:n3) = +u(is:ie,jss:jes:+1,k)*clm(:,1:ns) &
                            +v(is:ie,jss:jes:+1,k)*slm(:,1:ns)
        g (1:lngth,n3+1:n4) = +u(is:ie,jen:jsn:-1,k)*clm(:,1:nn) &
                            -v(is:ie,jen:jsn:-1,k)*slm(:,1:nn)
     else
        g (1:lngth,n0+1:n1) = u(is:ie,jss:jes:+1,k)
        g (1:lngth,n1+1:n2) = u(is:ie,jen:jsn:-1,k)
        g (1:lngth,n2+1:n3) = v(is:ie,jss:jes:+1,k)
        g (1:lngth,n3+1:n4) = v(is:ie,jen:jsn:-1,k)
     endif

     if (.not.Index%sigma .and. present(mask)) then
         gm(1:lngth,n0+1:n1) = mask(is:ie,jss:jes:+1,k)
         gm(1:lngth,n1+1:n2) = mask(is:ie,jen:jsn:-1,k)
         gm(1:lngth,n2+1:n4) = gm(1:lngth,n0+1:n2)
     endif

  enddo

!=======================================================================
!----- add cos lat to end of rows (needed for filter) ------------------

  if ( nlpf > 0 .and. len1 == lngth+1 ) then
     do k = 1, size(u,3)
        n0 = (k-1)*nlpf*2
        n2 = n0 + nlpf
        n4 = n2 + nlpf
        g (len1,n0+1:n2) = Index%cph (1:nlpf)
        g (len1,n2+1:n4) = Index%cph (1:nlpf)
     enddo
  endif

!=======================================================================
!-------------------- gather full zonal rows ---------------------------

  nfft = nlpf_zonal(0)

  if ( nlpf > 0 ) then

      call gather_zonal_rows (Index, nlpf_zonal, g, z(:,1:nfft) )

      if (.not.Index%sigma .and. present(mask)) then
        call gather_zonal_rows (Index, nlpf_zonal, gm, zm(:,1:nfft) )
      endif

  endif

!=======================================================================
!---------------- fill in missing values -------------------------------
!---------------- find  minimum value of mask --------------------------

    gm_min = 1.0
    if (.not.Index%sigma .and. present(mask)) then
        gm_min = minval(zm(1:leng,1:nfft))
        if (gm_min < 0.5) call fill_missing (zm(1:leng,1:nfft),  &
                                              z(1:leng,1:nfft))
    endif

!=======================================================================
!------- load balancing of fft rows ---------

  call fft_transmit_to ( z, nfft, put_pe, put_len, get_pe, get_len )

  if (nfft == 0) then
     !call mpp_clock_end (id_total)
     return      ! ?????
  endif

!=======================================================================
!-------- compute filter response ---------

    call set_filter_response ( Index, z(leng+1,1:nfft), ss(:,1:nfft) )

!=======================================================================
!---------------- transform to fourier coefficients --------------------

       c(:,1:nfft) = fft_grid_to_fourier (z(:,1:nfft))

!------------------------- filter --------------------------------------

       c (:,1:nfft) = c (:,1:nfft)* ss (:,1:nfft)

!---------------- transform back to grid spce --------------------------

       z(:,1:nfft) = fft_fourier_to_grid (c(:,1:nfft))

!=======================================================================
!----- return data to original pe ----

    nfft = nlpf_zonal(0)

    call fft_transmit_back ( z, nfft, put_pe, put_len, get_pe, get_len)

!=======================================================================
!--------------- restore zonal mean if mask used -----------------------

    if (gm_min < 0.5) call fix_missing (zm(1:leng,1:nfft),  &
                                         z(1:leng,1:nfft))

!=======================================================================
!------------------- scatter full zonal rows ---------------------------

    z(leng+1,:) = 0.0
    call scatter_zonal_rows ( Index, nlpf_zonal, z(:,1:nfft), g )

!=======================================================================
!-------------- return filtered data to original array -----------------

 do k=1,kd
     n0 = (k-1)*nlpf*2
     n1 = n0 + nlpfs
     n2 = n1 + nlpfn
     n3 = n2 + nlpfs
     n4 = n3 + nlpfn

     if (present(slm) .and. present(clm)) then
        u(is:ie,jss:jes:+1,k) = -g(1:lngth,n0+1:n1)*slm(:,1:ns) &
                                +g(1:lngth,n2+1:n3)*clm(:,1:ns)
        u(is:ie,jen:jsn:-1,k) = -g(1:lngth,n1+1:n2)*slm(:,1:nn) &
                                +g(1:lngth,n3+1:n4)*clm(:,1:nn)
        v(is:ie,jss:jes:+1,k) = +g(1:lngth,n0+1:n1)*clm(:,1:ns) &
                                +g(1:lngth,n2+1:n3)*slm(:,1:ns)
        v(is:ie,jen:jsn:-1,k) = -g(1:lngth,n1+1:n2)*clm(:,1:nn) &
                                -g(1:lngth,n3+1:n4)*slm(:,1:nn)
     else
        u(is:ie,jss:jes:+1,k) = g(1:lngth,n0+1:n1)
        u(is:ie,jen:jsn:-1,k) = g(1:lngth,n1+1:n2)
        v(is:ie,jss:jes:+1,k) = g(1:lngth,n2+1:n3)
        v(is:ie,jen:jsn:-1,k) = g(1:lngth,n3+1:n4)
     endif

 enddo

 !call mpp_clock_end (id_total)
!-----------------------------------------------------------------------

 end subroutine filter_two_fields

!#######################################################################

 subroutine polar_filter_4d (Control, data, grid, mask)

!-----------------------------------------------------------------------
   type(pfilt_control_type), intent(in)               :: Control
   real(r8),    intent(inout), dimension(:,:,:,:)         :: data
   integer, intent(in)                                :: grid
   real(r8),    intent(in),    dimension(:,:,:), optional :: mask
!-----------------------------------------------------------------------
   integer :: k, mxfft
   integer, dimension(Control%npes) :: put_pe, put_len, get_pe, get_len
   integer, dimension(0:Control%npes-1) :: nlpf_zonal

    select case (grid)

!   ------------------------------------------------------------------

      case (1)

         do k = 1, size(data,4), 2

            if ( k+1 <= size(data,4) ) then
                if ( k == 1 ) then
                   call set_zonal_balance (2*size(data,3),       &
                                           Control%Tmp%pelist2d, &
                                           Control%Tmp%nlpf2d,   &
                                           nlpf_zonal            )
                   call load_balance_filter (nlpf_zonal,         &
                          put_pe, put_len, get_pe, get_len, mxfft)
                 endif
                 mxfft = max (mxfft, Control%Tmp%nlpf*size(data,3))

                 call filter_two_fields (Control%Tmp, mxfft, nlpf_zonal,&
                                   put_pe, put_len, get_pe, get_len, &
                                   data(:,:,:,k), data(:,:,:,k+1),   &
                                   mask=mask)
            else
                 call set_zonal_balance (size(data,3),         &
                                         Control%Tmp%pelist2d, &
                                         Control%Tmp%nlpf2d,   &
                                         nlpf_zonal            )
                 call load_balance_filter (nlpf_zonal,   &
                              put_pe, put_len, get_pe, get_len, mxfft)
                 mxfft = max (mxfft, Control%Tmp%nlpf*size(data,3))
 
                 call filter_field (Control%Tmp, mxfft, nlpf_zonal,  &
                                    put_pe, put_len, get_pe, get_len, &
                                    data(:,:,:,k), mask)
            endif

         enddo

!   ------------------------------------------------------------------

      case (2:3)
         call set_zonal_balance (size(data,3), Control%Vel%pelist2d, &
                                 Control%Vel%nlpf2d, nlpf_zonal      )
         call load_balance_filter (nlpf_zonal, put_pe, put_len, &
                                   get_pe, get_len, mxfft)
         mxfft = max (mxfft, Control%Vel%nlpf*size(data,3))
         do k = 1, size(data,4)
            call filter_field (Control%Vel ,  mxfft, nlpf_zonal, &
                               put_pe, put_len, get_pe, get_len, &
                               data(:,:,:,k), mask)
         enddo

      case default
         call error_mesg ('polar_filter_4d', 'invalid grid arg', FATAL)

    end select

!-----------------------------------------------------------------------

 end subroutine polar_filter_4d

!#######################################################################

 subroutine polar_filter_3d (Control, data, grid, mask)

!-----------------------------------------------------------------------
   type(pfilt_control_type), intent(in)               :: Control
   real(r8),    intent(inout), dimension(:,:,:)           :: data
   integer, intent(in)                                :: grid
   real(r8),    intent(in),    dimension(:,:,:), optional :: mask
!-----------------------------------------------------------------------
   integer :: mxfft
   integer, dimension(Control%npes) :: put_pe, put_len, get_pe, get_len
   integer, dimension(0:Control%npes-1) :: nlpf_zonal
    select case (grid)

      case (1)
         call set_zonal_balance (size(data,3), Control%Tmp%pelist2d, &
                                 Control%Tmp%nlpf2d, nlpf_zonal      )
         call load_balance_filter (nlpf_zonal, put_pe, put_len, &
                                   get_pe, get_len, mxfft)
         mxfft = max (mxfft, Control%Tmp%nlpf*size(data,3))
         call filter_field (Control%Tmp,  mxfft, nlpf_zonal,  &
                            put_pe, put_len, get_pe, get_len, &
                            data, mask)

      case (2:3)
         call set_zonal_balance (size(data,3), Control%Vel%pelist2d, &
                                 Control%Vel%nlpf2d, nlpf_zonal      )
         call load_balance_filter (nlpf_zonal, put_pe, put_len, &
                                   get_pe, get_len, mxfft)
         mxfft = max (mxfft, Control%Vel%nlpf*size(data,3))
         call filter_field (Control%Vel ,  mxfft, nlpf_zonal, &
                            put_pe, put_len, get_pe, get_len, &
                            data, mask)

      case default
         call error_mesg ('polar_filter_3d', 'invalid grid arg', FATAL)

    end select

!-----------------------------------------------------------------------

 end subroutine polar_filter_3d

!#######################################################################

 subroutine polar_filter_2d (Control, data, grid, mask)

!-----------------------------------------------------------------------
   type(pfilt_control_type), intent(in)             :: Control
   real(r8),    intent(inout), dimension(:,:)           :: data
   integer, intent(in)                              :: grid
   real(r8),    intent(in),    dimension(:,:), optional :: mask
!-----------------------------------------------------------------------
   real(r8), dimension(size(data,1),size(data,2),1) :: data3, mask3
   integer :: mxfft
   integer, dimension(Control%npes) :: put_pe, put_len, get_pe, get_len
   integer, dimension(0:Control%npes-1) :: nlpf_zonal
!-----------------------------------------------------------------------

    data3(:,:,1) = data(:,:)

    select case (grid)

      case (1)
         call set_zonal_balance ( 1, Control%Tmp%pelist2d,      &
                                 Control%Tmp%nlpf2d, nlpf_zonal )
         call load_balance_filter (nlpf_zonal, put_pe, put_len, &
                                   get_pe, get_len, mxfft)
         mxfft = max (mxfft, Control%Tmp%nlpf)
         if (present(mask)) then
             mask3(:,:,1) = mask(:,:)
             call filter_field (Control%Tmp,  mxfft, nlpf_zonal,  &
                                put_pe, put_len, get_pe, get_len, &
                                data3, mask3)
         else
             call filter_field (Control%Tmp,  mxfft, nlpf_zonal,  &
                                put_pe, put_len, get_pe, get_len, &
                                data3)
         endif

      case (2:3)
         call set_zonal_balance ( 1, Control%Vel%pelist2d,      &
                                 Control%Vel%nlpf2d, nlpf_zonal )
         call load_balance_filter (nlpf_zonal, put_pe, put_len, &
                                   get_pe, get_len, mxfft)
         mxfft = max (mxfft, Control%Vel%nlpf)
         if (present(mask)) then
             mask3(:,:,1) = mask(:,:)
             call filter_field (Control%Vel ,  mxfft, nlpf_zonal, &
                                put_pe, put_len, get_pe, get_len, &
                                data3, mask3)
         else
             call filter_field (Control%Vel ,  mxfft, nlpf_zonal, &
                                put_pe, put_len, get_pe, get_len, &
                                data3)
         endif

      case default
         call error_mesg ('polar_filter_2d', 'invalid grid arg', FATAL)

    end select

    data(:,:) = data3(:,:,1)

!-----------------------------------------------------------------------

 end subroutine polar_filter_2d

!#######################################################################

 subroutine polar_filter_two_3d (Control, u, v, grid, mask)

!-----------------------------------------------------------------------
   type(pfilt_control_type), intent(in)               :: Control
   real(r8),    intent(inout), dimension(:,:,:)           :: u, v
   integer, intent(in)                                :: grid
   real(r8),    intent(in),    dimension(:,:,:), optional :: mask
!-----------------------------------------------------------------------
   integer :: mxfft
   integer, dimension(Control%npes) :: put_pe, put_len, get_pe, get_len
   integer, dimension(0:Control%npes-1) :: nlpf_zonal

    select case (grid)

      case (1)
         call set_zonal_balance (2*size(u,3), Control%Tmp%pelist2d, &
                                 Control%Tmp%nlpf2d, nlpf_zonal     )
         call load_balance_filter (nlpf_zonal, put_pe, put_len, &
                                   get_pe, get_len, mxfft)
         mxfft = max (mxfft, Control%Tmp%nlpf*size(u,3))
         call filter_two_fields (Control%Tmp, mxfft, nlpf_zonal, &
                   put_pe, put_len, get_pe, get_len, u, v, mask=mask)

      case (2:3)
         call set_zonal_balance (2*size(u,3), Control%Vel%pelist2d, &
                                 Control%Vel%nlpf2d, nlpf_zonal     )
         call load_balance_filter (nlpf_zonal, put_pe, put_len, &
                                   get_pe, get_len, mxfft)
         mxfft = max (mxfft, Control%Vel%nlpf*size(u,3))
         call filter_two_fields (Control%Vel, mxfft, nlpf_zonal, &
                            put_pe, put_len, get_pe, get_len, &
                            u, v, mask=mask)

      case default
         call error_mesg ('polar_filter_two_3d', 'invalid grid arg', &
                          FATAL)

    end select

!-----------------------------------------------------------------------

 end subroutine polar_filter_two_3d

!#######################################################################

 subroutine polar_filter_wind_3d (Control, u, v, mask)

!-----------------------------------------------------------------------
   type(pfilt_control_type), intent(in)             :: Control
   real(r8),  intent(inout), dimension(:,:,:)           :: u, v
   real(r8),  intent(in),    dimension(:,:,:), optional :: mask
!-----------------------------------------------------------------------
   integer :: mxfft
   integer, dimension(Control%npes) :: put_pe, put_len, get_pe, get_len
   integer, dimension(0:Control%npes-1) :: nlpf_zonal


     call set_zonal_balance (2*size(u,3), Control%Vel%pelist2d, &
                             Control%Vel%nlpf2d, nlpf_zonal     )
     call load_balance_filter (nlpf_zonal, put_pe, put_len, &
                               get_pe, get_len, mxfft)
     mxfft = max (mxfft, Control%Vel%nlpf*size(u,3))

     call filter_two_fields (Control%Vel, mxfft, nlpf_zonal,     &
                             put_pe, put_len, get_pe, get_len,   &
                             u, v, Control%slm, Control%clm, mask)

!-----------------------------------------------------------------------

 end subroutine polar_filter_wind_3d

!#######################################################################

 subroutine polar_filter_wind_2d (Control, u, v, mask)

!-----------------------------------------------------------------------
   type(pfilt_control_type), intent(in)           :: Control
   real(r8),  intent(inout), dimension(:,:)           :: u, v
   real(r8),  intent(in),    dimension(:,:), optional :: mask
!-----------------------------------------------------------------------
   real(r8), dimension(size(u,1),size(u,2),1) :: u3, v3, mask3
   integer :: mxfft
   integer, dimension(Control%npes) :: put_pe, put_len, get_pe, get_len
   integer, dimension(0:Control%npes-1) :: nlpf_zonal
!-----------------------------------------------------------------------

     call set_zonal_balance ( 2, Control%Vel%pelist2d,      &
                             Control%Vel%nlpf2d, nlpf_zonal )
     call load_balance_filter (nlpf_zonal, put_pe, put_len, &
                               get_pe, get_len, mxfft)
     mxfft = max (mxfft, Control%Vel%nlpf)

      u3(:,:,1) = u(:,:)
      v3(:,:,1) = v(:,:)

      if (present(mask)) then
         mask3(:,:,1) = mask(:,:)
         call filter_two_fields (Control%Vel, mxfft, nlpf_zonal,   &
                                 put_pe, put_len, get_pe, get_len, &
                                 u3, v3, Control%slm, Control%clm, &
                                 mask=mask3)
      else
         call filter_two_fields (Control%Vel, mxfft, nlpf_zonal,   &
                                 put_pe, put_len, get_pe, get_len, &
                                 u3, v3, Control%slm, Control%clm  )
      endif

      u(:,:) = u3(:,:,1)
      v(:,:) = v3(:,:,1)

!-----------------------------------------------------------------------

 end subroutine polar_filter_wind_2d

!#######################################################################

 subroutine fill_missing (mask,data)

!-----------------------------------------------------------------------
      real(r8), intent(in),    dimension(:,:) :: mask
      real(r8), intent(inout), dimension(:,:) :: data
!-----------------------------------------------------------------------
      integer :: k, id, kd
      real(r8)    :: rcnt

      id = size(data,1); kd = size(data,2)

      allocate (icnt(kd), avg(kd))

!---------------- compute mean of each transform row -------------------
!---------------- then fill in missing values --------------------------

      do k=1,kd
         icnt(k)=count(mask(1:id,k) > 0.50)
         if (icnt(k) == id .or. icnt(k) == 0) cycle
         rcnt=1.0/(1.0 * icnt(k))
         avg(k)=sum(data(1:id,k)*mask(1:id,k))*rcnt
 
         call intrp (mask(1:id,k),data(1:id,k))
      enddo

!-----------------------------------------------------------------------

 end subroutine fill_missing

!#######################################################################

 subroutine fix_missing (mask,data)

!-----------------------------------------------------------------------
   real(r8), intent(in),    dimension(:,:) :: mask
   real(r8), intent(inout), dimension(:,:) :: data
!-----------------------------------------------------------------------
   integer :: k, id, kd
   real(r8)    :: rcnt, dif

      id = size(data,1); kd = size(data,2)

!---------------- restore zonal mean values (avg) ----------------------

      do k=1,kd
         if (icnt(k) == id .or. icnt(k) == 0) cycle

         rcnt=1.0/icnt(k)
         dif=avg(k)-sum(data(1:id,k)*mask(1:id,k))*rcnt
         data(1:id,k)=(data(1:id,k)+dif)*mask(1:id,k)
      enddo

      deallocate (icnt, avg)

!-----------------------------------------------------------------------

 end subroutine fix_missing

!#######################################################################

 subroutine intrp (mask,a)

!-----------------------------------------------------------------------
   real(r8), intent(in),    dimension(:) :: mask
   real(r8), intent(inout), dimension(:) :: a
!-----------------------------------------------------------------------
   integer, dimension(size(a,1)) :: m1,m2,mbas
   real(r8),    dimension(size(a,1)) :: base,slop
   integer :: lngth,last,nseg,i,m,n
!-----------------------------------------------------------------------
!  fill in missing values by linear interpolating

      m1(:)=99999
      m2(:)=99999
      mbas(:)=99999
      base(:)=99999.
      slop(:)=99999.

      lngth=size(a,1)
      last=-1
      nseg=1

      do i=1,lngth
          if ( mask(i) < 0.50 ) then
              if (last == 1) then
                 m1(nseg)=i
                 mbas(nseg)=i-1
                 base(nseg)=a(i-1)
              endif
              last=0
          else if ( mask(i) > 0.50 ) then
              if (last == 0) then
                 m2(nseg)=i-1
                 slop(nseg)=(a(i)-base(nseg))/float(i-m1(nseg))
                 nseg=nseg+1
              endif
              last=1
          endif
      enddo

      if ( m1(nseg) == 99999 .and. m2(nseg) == 99999) nseg=nseg-1

      if ( m1(1) == 99999 .and. m2(nseg) == 99999) then
          m1(1)=1
          m2(nseg)=lngth
          mbas(1)=mbas(nseg)-lngth
          base(1)=base(nseg)
          slop(1)=(a(m2(1)+1)-base(1))/float(m2(1)+1-mbas(1))
          slop(nseg)=slop(1)
      endif
      if ( m1(1) == 99999 ) then
          m1(1)=1
          mbas(1)=0
          base(1)=a(lngth)
          slop(1)=(a(m2(1)+1)-base(1))/float(m2(1)+1-mbas(1))
      endif
      if ( m2(nseg) == 99999 ) then
          m2(nseg)=lngth
          slop(nseg)=(a(1)-base(nseg))/float(lngth+1-mbas(nseg))
      endif

      do n=1,nseg
          do m=m1(n),m2(n)
              a(m)=base(n)+float(m-mbas(n))*slop(n)
          enddo
      enddo

!-----------------------------------------------------------------------

 end subroutine intrp

!#######################################################################

 function polar_filter_init (Hgrid, reflat, weight, sigma, verbose)  &
                     result (Control)

   type(horiz_grid_type), intent(in)           :: Hgrid
   real(r8),                  intent(in), optional :: reflat
   integer,               intent(in), optional :: weight, verbose
   logical,               intent(in), optional :: sigma

   type(pfilt_control_type) :: Control

!-----------------------------------------------------------------------
!  reflat = reference latitiude in degrees (default = 60.)
!  weight = weight applied to filter response function (ss)
!           e.g., ss ** weight (default: weight = 1)
!  sigma  = sigma or eta/step-mtn coordinate (default: sigma = .false.)
!-----------------------------------------------------------------------
   integer :: k, n, km, lngth, len1, lenc, iverbose, nlpf
   real(r8)    :: rlat,  dtr, hpi, filter_lats

!  integer :: hsg, heg, hs, he, vsg, veg, vs, ve
!  integer :: hss, hes, hsn, hen, vss, ves, vsn, ven
!  integer :: nlpf_hs, nlpf_hn, nlpf_h
!  integer :: nlpf_vs, nlpf_vn, nlpf_v
!  integer :: npes, pe, mxfft, leng, iloc, jloc
   integer :: mxfft, leng
!-----------------------------------------------------------------------

 if (do_log) then
   call write_version_number (version,tagname)
   do_log = .false.
 endif

 if (do_clock_init) then
   !id_total = mpp_clock_init ('BGRID: polar_filter (TOTAL)', time_level, flags=MPP_CLOCK_SYNC)
   do_clock_init = .false.
 endif

      iverbose = 0;  if (present(verbose)) iverbose = verbose

!   --- define variables used by all routines in this module ---

      filter_lats = 30.
      if (present(reflat)) then
          if ( reflat-epsilon(reflat) <= 90. .and. reflat >= 0. ) then
             filter_lats = 90. - reflat
          else
             call error_mesg  ('polar_filter_init',      &
                'reflat must lie between 0 and 90.', FATAL)
          endif
      endif

!   --- number of latitude rows of filtering per hemisphere ---

      nlpf = int(float(Hgrid%Tmp%jeg-Hgrid%Tmp%jsg+1)*filter_lats/180.)
      Control%nlpf = nlpf

!-----------------------------------------------------------------------
!----- optional (additional) weight for filter reponse function ------

      Control%Tmp%weight = 1
      Control%Vel%weight = 1
  if (present(weight)) then
      Control%Tmp%weight = weight
      Control%Vel%weight = weight
  endif

!----- optional sigma coordinate flag -----

      Control%Tmp%sigma = .false.
      Control%Tmp%sigma = .false.
  if (present(sigma)) then
      Control%Tmp%sigma = sigma
      Control%Vel%sigma = sigma
  endif

!---------------------------------------------------------------------
  Control%npes = 1
  call setup_index_type ( nlpf, Hgrid%Tmp, Control%Tmp )
  call setup_index_type ( nlpf, Hgrid%Vel, Control%Vel )

!-----------------------------------------------------------------------
!     ----- more indices (the computational domain) -----

      Control % Tmp % ilb = Hgrid%ilb 
      Control % Tmp % jlb = Hgrid%jlb 
      Control % Tmp % iub = Hgrid%iub 
      Control % Tmp % jub = Hgrid%jub 

      Control % Tmp % is = Hgrid%Tmp%is
      Control % Tmp % ie = Hgrid%Tmp%ie
      Control % Tmp % js = Hgrid%Tmp%js
      Control % Tmp % je = Hgrid%Tmp%je

      Control % Vel % ilb = Hgrid%ilb 
      Control % Vel % jlb = Hgrid%jlb 
      Control % Vel % iub = Hgrid%iub 
      Control % Vel % jub = Hgrid%jub 

      Control % Vel % is = Hgrid%Vel%is
      Control % Vel % ie = Hgrid%Vel%ie
      Control % Vel % js = Hgrid%Vel%js
      Control % Vel % je = Hgrid%Vel%je

!     ----- num lon points and max wave number ----
!     (assume same number of x-points for both Tmp and Vel grids)

      lngth  = Hgrid%Tmp%ie - Hgrid%Tmp%is + 1
      len1 = lngth
!     east-most box in zonal row needs one more point
!     this is used to store the latitude value
      if (Hgrid%Tmp%ie == Hgrid%Tmp%ieg) len1 = len1 + 1

      leng = Hgrid%nlon
      km   = leng/2
      lenc = km+1

      Control % Tmp % len  = lngth
      Control % Tmp % len1 = len1
      Control % Tmp % leng = leng
      Control % Tmp % lenc = lenc
      Control % Vel % len  = lngth
      Control % Vel % len1 = len1
      Control % Vel % leng = leng
      Control % Vel % lenc = lenc


!  ---- trigometric constants ----

      hpi = acos(0.0); dtr = hpi/90.

!   --- fourier transform initialization ----

      call fft_init (leng)

!-----------------------------------------------------------------------
! **** trig constants for converting u,v to polar stereographic ****

   mxfft = max (Control%Vel%nlpf_s, Control%Vel%nlpf_n)

   allocate (Control%slm(lngth,mxfft), Control%clm(lngth,mxfft))

   if (mxfft > 0) then
         Control%slm(:,1) = sin(Hgrid%Vel%tlm(Hgrid%Vel%is:Hgrid%Vel%ie,Hgrid%Vel%js))
         Control%clm(:,1) = cos(Hgrid%Vel%tlm(Hgrid%Vel%is:Hgrid%Vel%ie,Hgrid%Vel%js))
      do k = 2, mxfft
         Control%slm(:,k) = Control%slm(:,1)
         Control%clm(:,k) = Control%clm(:,1)
      enddo
   endif

!-----------------------------------------------------------------------
!----- initialize geography constants for polar filter function --------

!     *** ref lat (first lat w/o filtering) ****

     rlat = hpi - (float(Control%nlpf)+0.50)*Hgrid%dph

     Control%Tmp%cph0 = cos(rlat)
     Control%Vel%cph0 = cos(rlat)

     Control%Tmp%dlm = Hgrid%dlm
     Control%Vel%dlm = Hgrid%dlm

    call set_filter_latitudes ( Control%Tmp, Hgrid%Tmp%tph(Hgrid%Tmp%is,:) )
    call set_filter_latitudes ( Control%Vel, Hgrid%Vel%tph(Hgrid%Vel%is,:) )

!-----------------------------------------------------------------------

 end function polar_filter_init

!#######################################################################

 subroutine set_indexing (nlpf, jsg, jeg, js, je, jss, jes, jsn, jen)

   integer, intent(in)  :: nlpf, jsg, jeg, js,  je
   integer, intent(out) ::       jss, jes, jsn, jen

   integer :: nlpf_s, nlpf_n

!  note: nlpf_s = jes-jss+1
!        nlpf_n = jen-jsn+1

!     ----- starting/ending local indices -----
!     ----- adjust end index for zero length ----

!     --- southern hemisphere ---

      jss    = max (jsg, js)
      jes    = min (jsg+nlpf-1, je)
      nlpf_s = max (0, jes - jss + 1)
      jes    = jss + nlpf_s - 1

!     --- northern hemisphere ---

      jsn    = max (jeg-nlpf+1, js)
      jen    = min (jeg, je)
      nlpf_n = max (0, jen - jsn + 1)
      jen    = jsn + nlpf_n - 1

 end subroutine set_indexing

!#######################################################################

 subroutine set_filter_latitudes (Index, ph)

  type(pfilt_index_type), intent(inout)  :: Index
  real(r8)                  , intent(in)     :: ph (Index%jlb:)

  integer ::  j, k

   allocate ( Index%cph (Index%nlpf), &
              Index%sklm(Index%lenc) )


!--- southern hemisphere (s -> n) ---

    do j = 1, Index%nlpf_s
        Index % cph(j) = cos( ph(Index%jss+j-1) )
    enddo

!--- northern hemisphere (n -> s) ---

    do j = 1, Index%nlpf_n
        Index % cph(Index%nlpf_s+j) = cos( ph(Index%jen-j+1) )
    enddo

!--- sin of wave number ---

    do k = 2, Index%lenc
        Index % sklm(k) = sin( 0.5*float(k-1)*Index%dlm )
    enddo

 end subroutine set_filter_latitudes

!#######################################################################

 subroutine set_filter_response (Index, cph, ss)

  type(pfilt_index_type), intent(in)  :: Index
  real(r8)                  , intent(in)  :: cph (:)
  real(r8)                  , intent(out) :: ss (:,:)

  real(r8), dimension(Index%lenc) :: cph0_sklm
  integer :: j, k

! **** standard filter response (after arakawa & lamb) ****

if (size(cph) == 0) return

    do k = 2, Index%lenc
      cph0_sklm (k) = Index%cph0 * Index%sklm(k)
    enddo

!--- mean always one ---

    ss (1,:) = 1.0

!--- compute response function (range 0. to 1.) -----

    if ( Index%weight == 1 ) then
       do j = 1, size(ss,2)
       do k = 2, Index%lenc
           ss (k,j) = min( 1.0_r8, cph(j)/cph0_sklm(k) )
       enddo
       enddo
     else
       do j = 1, size(ss,2)
       do k = 2, Index%lenc
           ss (k,j) = min( 1.0_r8, (cph(j)/cph0_sklm(k))**Index%weight )
       enddo
       enddo
     endif


 end subroutine set_filter_response

!#######################################################################

 subroutine set_zonal_balance (nlev, pelist2d, nlpf2d, nlpf)
  integer, intent(in)  :: nlev, pelist2d(:,:), nlpf2d(:,:)
  integer, intent(out) :: nlpf(0:)

  integer, dimension(size(nlpf2d,1)) :: nlpf2d_zonal
  integer :: i, j, remain

! the returned value is the number of zonal rows on each processor
! initialize this value to zero
    nlpf = 0

    do j=1,size(nlpf2d,2)
!      skip pe rows with no filter rows
       if (nlpf2d(1,j) == 0) cycle

       nlpf2d_zonal(:) = (nlpf2d(1,j)*nlev)/size(nlpf2d,1)
       remain = nlpf2d(1,j)*nlev - sum(nlpf2d_zonal(:))
       do i=1,size(nlpf2d,1)
         if (remain == 0) exit
         nlpf2d_zonal(i) = nlpf2d_zonal(i) + 1
         remain = remain - 1
       enddo

!      return 1-dimensinal structure
       do i=1,size(nlpf2d,1)
         nlpf(pelist2d(i,j)) = nlpf2d_zonal(i)
       enddo
    enddo

 end subroutine set_zonal_balance

!#######################################################################

subroutine load_balance_filter (nlpf, put_pe, put_len, &
                                      get_pe, get_len, mxfft)

    integer, intent(in)  :: nlpf(0:)
    integer, intent(out) :: put_pe(:), put_len(:), &
                            get_pe(:), get_len(:), mxfft

    integer, dimension(0:size(nlpf)-1) :: nlpf_done, nlpf_need
    integer :: i, j, pe, nlpf0, mm, nn, need, get, ig, ip
    integer :: js, je, ji

!--- initialize ---

      put_pe = -1 ;  put_len = 0
      get_pe = -1 ;  get_len = 0

      pe = 0

!--- skip load balance ??? ---

      if (.not.do_load_balance) then
         mxfft = nlpf(pe)
         return
      endif

!--------------------- load balancing of FFTs --------------------------
!      each PE should have approximately the same number of FFTs

      nlpf0 = sum(nlpf)
      nlpf_need = nlpf0/size(nlpf)

      nn = nlpf0 - sum(nlpf_need)
      mm = count( nlpf > nlpf0/size(nlpf) )

      if ( nn >= mm ) then
          where ( nlpf > nlpf0/size(nlpf) ) nlpf_need = nlpf_need+1
          nn = nn - mm
      endif

      do i = 0, size(nlpf)-1
         if (nn == 0) exit
         if ( nlpf_need(i) == nlpf0/size(nlpf) ) then
              nlpf_need(i) = nlpf_need(i) + 1
              nn = nn - 1
         endif
      enddo

       ig=0; ip=0
       put_pe = -1;  put_len = 0
       get_pe = -1;  get_len = 0

       nlpf_done  = nlpf

      do i = 0, size(nlpf)-1
          if ( nlpf_done(i) >= nlpf_need(i) ) cycle

!         --- pe i needs rows ---
          need = nlpf_need(i) - nlpf_done(i)

!         --- search pe's for rows ---
          if ( 2*i < size(nlpf) ) then
             js=0; je=size(nlpf)-1; ji=1
          else
!            --- reverse search order for high number pe's ----
             js=size(nlpf)-1; je=0; ji=-1
          endif

          do j = js, je, ji
              if ( nlpf_done(j) <= nlpf_need(j) ) cycle
!             --- pe j has rows ---

              get = min( need, nlpf_done(j) - nlpf_need(j) )
              nlpf_done(i) = nlpf_done(i) + get
              nlpf_done(j) = nlpf_done(j) - get
              if (i == pe) then
                  ig = ig+1
                  get_pe (ig) = j
                  get_len(ig) = get
              endif
              if (j == pe) then
                  ip = ip+1
                  put_pe (ip) = i
                  put_len(ip) = get
              endif
              need = nlpf_need(i) - nlpf_done(i)
              if (need == 0) exit
          enddo
      enddo

      mxfft = max( nlpf_done(pe), nlpf(pe) )

      if ( sum(put_len) /=0 .and. sum(get_len) /=0 ) then
          call error_mesg  ('load_balance_filter in bgrid_polar_filter_mod',&
                'cannot get and put fft rows on the same pe', FATAL)
      endif

 end subroutine load_balance_filter

!#######################################################################

 subroutine fft_transmit_to ( g, nfft, put_pe, put_len, &
                                       get_pe, get_len  )

 real(r8),    intent(inout), dimension(:,:) :: g
 integer, intent(inout)                 :: nfft
 integer, intent(in),    dimension(:)   :: put_pe, put_len, &
                                           get_pe, get_len

 integer :: i, n, len1, pfft, plen1, gfft, glen1
!-----------------------------------------------------------------------

   len1 = size(g ,1)


!----- put fft rows ------

   nfft = nfft - sum(put_len)
   n = nfft

!  call mpp_sync ()

   do i = 1, size(put_pe)
      if (put_len(i) == 0) exit
      pfft = put_len(i);  plen1 = pfft*len1
!print *, 'pe (put) = ',mpp_pe(),plen1

!!!      call mpp_send ( g(:,n+1), plen1, put_pe(i) )
!!!      call mpp_sync_self ()

      n = n + pfft
   enddo

!----- get fft rows ------

   do i = 1, size(get_pe)
      if (get_len(i) == 0) exit
      gfft = get_len(i);  glen1 = gfft*len1
!print *, 'pe (get) = ',mpp_pe(),glen1

!!!      call mpp_recv ( g(:,nfft+1), glen1, get_pe(i) )
!!!      call mpp_sync_self ()

      nfft = nfft + gfft
   enddo

!  call mpp_sync ()
!-----------------------------------------------------------------------

 end subroutine fft_transmit_to

!#######################################################################

 subroutine fft_transmit_back ( g, nfft, put_pe, put_len, &
                                         get_pe, get_len  )

 real(r8),    intent(inout), dimension(:,:) :: g
 integer, intent(in)                    :: nfft
 integer, intent(in),    dimension(:)   :: put_pe, put_len, &
                                           get_pe, get_len

 integer :: i, n, pfft, plen1, gfft, glen1
 integer :: len1
!-----------------------------------------------------------------------

 len1 = size(g ,1)

!----- return data to original pe ----

   n = nfft - sum(put_len)

!----- put fft rows (that were got) ------

!  call mpp_sync ()

   do i = 1, size(get_pe)
      if (get_len(i) == 0) exit
      gfft = get_len(i);  glen1 = gfft*len1

!!!      call mpp_send ( g(:,n+1), glen1, get_pe(i) )
!!!      call mpp_sync_self ()

      n = n + gfft
   enddo

!----- get fft rows (that were put) ------
   do i = 1, size(put_pe)
      if (put_len(i) == 0) exit
      pfft = put_len(i);  plen1 = pfft*len1

!!!      call mpp_recv ( g(:,n+1), plen1, put_pe(i) )
!!!      call mpp_sync_self ()

      n = n + pfft
   enddo

!  call mpp_sync ()
!-----------------------------------------------------------------------

 end subroutine fft_transmit_back

!#######################################################################

 subroutine gather_zonal_rows (Index, nlpf, local, zonal )

    type(pfilt_index_type), intent(in)  :: Index
    integer, intent(in),  dimension(0:) :: nlpf
    real(r8),    intent(in),  dimension(:,:) :: local
    real(r8),    intent(out), dimension(:,:) :: zonal

    real(r8), dimension(1:Index%maxlen*size(zonal,2)) :: temp

    integer :: i, j, is, ie, js, je, ks, ke, glen, plen, get_pe

!-----------------------------------------------------------------------
!-------------------- gather full zonal rows ---------------------------

  if ( size(local,2) == 0 ) return

  is = 1; js = 1

  do i = 1, size(Index%pelist)

    plen = size(local,1) * nlpf(Index%pelist(i))
    glen = Index%sizelist(i) * nlpf(0)
    get_pe = Index%pelist(i)
    if (glen == 0) get_pe = -3
    ie = is + Index%sizelist(i) - 1
    je = js + nlpf(Index%pelist(i)) - 1

    if ( Index%pelist(i) == 0 ) then
       zonal (is:ie, 1:nlpf(0)) = local(:,js:je)
    else
! This block of code is never reached with single pe
       !!!call mpp_transmit ( local(:,js), plen, Index%pelist(i), &
       !!!                    temp    , glen,       get_pe     )
!       call mpp_sync_self ()
!      reshape
       do j = 1, nlpf(0)
          ks = (j-1)*Index%sizelist(i)+1
          ke =  j   *Index%sizelist(i)
          zonal (is:ie,j) = temp(ks:ke)
       enddo
    endif

    is = ie+1
    js = je+1
  enddo

 end subroutine gather_zonal_rows

!#######################################################################

 subroutine scatter_zonal_rows (Index, nlpf, zonal, local )

    type(pfilt_index_type), intent(in)  :: Index
    integer, intent(in),  dimension(0:) :: nlpf
    real(r8),    intent(in),  dimension(:,:) :: zonal
    real(r8),    intent(out), dimension(:,:) :: local

    real(r8), dimension(1:Index%maxlen*size(zonal,2)) :: temp

    integer :: i, j, is, ie, js, je, ks, ke, glen, plen, put_pe

!-----------------------------------------------------------------------
!-------------------- scatter full zonal rows --------------------------

  if ( size(local,2) == 0 ) return

  is = 1; js = 1

  do i = 1, size(Index%pelist)

    glen = size(local,1) * nlpf(Index%pelist(i))
    plen = Index%sizelist(i) * nlpf(0)
    put_pe = Index%pelist(i)
    if (plen == 0) put_pe = -3
    ie = is + Index%sizelist(i) - 1
    je = js + nlpf(Index%pelist(i)) - 1

    if ( Index%pelist(i) == 0 ) then
       local(:,js:je) = zonal (is:ie, 1:nlpf(0))
    else
!      reshape
       do j = 1, nlpf(0)
          ks = (j-1)*Index%sizelist(i)+1
          ke =  j   *Index%sizelist(i)
          temp(ks:ke) = zonal (is:ie,j)
       enddo
! For single pe this block is never reached'
       !!!call mpp_transmit ( temp    , plen,       put_pe  , &
       !!!                    local(:,js), glen, Index%pelist(i) )
!       call mpp_sync_self ()
    endif

    is = ie+1
    js = je+1
  enddo

 end subroutine scatter_zonal_rows

!#######################################################################

 subroutine setup_index_type ( nlpfg, Bgrid, Index )
   integer,                intent(in)    :: nlpfg
   type(bgrid_type),         intent(in)    :: Bgrid
   type(pfilt_index_type), intent(inout) :: Index

!   type(domain1d) :: X, Y
   integer :: npes, layout(2), pe, iloc, jloc
   integer :: i, j, is, js, je, isg, ieg, jsg, jeg
   integer :: jss, jes, jsn, jen, nlpf_s, nlpf_n, nlpf

   integer, allocatable, dimension(:) :: ypelist, ysizelist,  &
                           xbeg, xend, xsize, ybeg, yend, ysize


     npes = 1
     !call mpp_get_layout ( Domain, layout )
     ! Layout is 1 domain each direction for 1 pe
     layout(1) = 1; layout(2) = 1

     allocate ( Index%pelist(layout(1)), Index%sizelist(layout(1)) )
     allocate ( Index%nlpf2d  (layout(1),layout(2)), &
                Index%pelist2d(layout(1),layout(2))  )

     allocate ( ypelist(layout(2)), ysizelist(layout(2)) )
     allocate ( xbeg(0:npes-1), xend(0:npes-1), xsize(0:npes-1), &
                ybeg(0:npes-1), yend(0:npes-1), ysize(0:npes-1)  )

     !call mpp_get_domain_components ( Domain, X, Y )
     !call mpp_get_pelist ( X, Index%pelist )
     !call mpp_get_pelist ( Y,      ypelist )
     ! For one pe, the pelists are just 0
     Index%pelist = 0
     ypelist = 0

     !!!call mpp_get_compute_domains ( Domain, xbeg, xend, xsize, &
     !!!                                       ybeg, yend, ysize  )
     ! For one pe can get sizes for compute domain directly
     xbeg = Bgrid%is; xend = Bgrid%ie; xsize = Bgrid%ie - Bgrid%is + 1
     ybeg = Bgrid%js; yend = Bgrid%je; ysize = Bgrid%je - Bgrid%js + 1
   
! compute size lists along each axis

     do i = 1, layout(1)
       Index%sizelist(i) = xsize(Index%pelist(i))
     enddo
       Index%sizelist(layout(1)) = Index%sizelist(layout(1)) + 1
       Index%maxlen = maxval( Index%sizelist )

     do j = 1, layout(2)
       ysizelist(j) = ysize(ypelist(j))
     enddo

     !call mpp_get_global_domain ( Domain, isg, ieg, jsg, jeg )
     ! For one pe, global and compute domains are the same
     isg = xbeg(0); ieg = xend(0); jsg = ybeg(0); jeg = yend(0)

! determine filtering on all PEs

     do pe = 0, npes-1

        is  = xbeg(pe)
        js  = ybeg(pe)
        je  = yend(pe)

        iloc = getloc ( is, isg, Index%sizelist )
        jloc = getloc ( js, jsg,      ysizelist )
        Index%pelist2d(iloc,jloc) = pe

        call set_indexing ( nlpfg, jsg, jeg, js, je, jss, jes, jsn, jen )
        nlpf_s = jes-jss+1         ! s hemis rows
        nlpf_n = jen-jsn+1         ! n hemis rows
        nlpf   = nlpf_s + nlpf_n   ! total rows of polar filtering

        Index%nlpf2d(iloc,jloc) = nlpf

        !if ( pe == mpp_pe() ) then
            Index % jss    = jss
            Index % jes    = jes
            Index % jsn    = jsn
            Index % jen    = jen
            Index % nlpf_s = nlpf_s
            Index % nlpf_n = nlpf_n
            Index % nlpf   = nlpf
            Index % iloc   = iloc
            Index % jloc   = jloc
        !endif

     enddo

     deallocate ( ypelist, ysizelist, xbeg, xend, xsize, ybeg, yend, ysize )

 end subroutine setup_index_type

!#######################################################################
!  this function return the array index i, where value = array(i)

 function getloc ( is, isg, sizelist )
 integer, intent(in) :: is, isg, sizelist(:)
 integer :: getloc, iloc, n
   getloc = 0
   iloc = isg
   do n = 1, size(sizelist)
      if ( is == iloc ) then
          getloc = n
          exit
      else
          iloc = iloc + sizelist(n)
      endif
   enddo
 end function getloc

!#######################################################################

end module bgrid_polar_filter_mod

