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

module bgrid_horiz_mod

!-----------------------------------------------------------------------
!
!        allocates storage and initializes grid constants
!              for the FMS B-grid dynamical core
!
!-----------------------------------------------------------------------

use types_mod, only : r8
use  constants_mod, only: RADIUS
use        fms_mod, only: error_mesg, FATAL, NOTE

implicit none
private

!-----------------------------------------------------------------------
!------- public interfaces -------

public  horiz_grid_init,       &
        get_horiz_grid_bound,  &
        get_horiz_grid_size,   &
        update_np, update_sp,  &
        TGRID, VGRID

!-----------------------------------------------------------------------
!             public derived data types
!-----------------------------------------------------------------------

public bgrid_type
public horiz_grid_type

!-----------------------------------------------------------------------
!
!    NOTE: all horizontal indexing uses global indices
!                    ( 1:nlon, 1:nlat )

type bgrid_type
   integer :: is,  ie,  js,  je         ! compute domain indices
   integer :: isd, ied, jsd, jed        ! data    domain indices
   integer :: isg, ieg, jsg, jeg        ! global  domain indices
   real(r8), pointer, dimension(:)   :: blong, blatg   ! global grid edges
   real(r8), pointer, dimension(:,:) :: dx, rdx, area, rarea
   real(r8), pointer, dimension(:,:) :: tph, tlm, aph, alm
   real(r8)                          :: dy, rdy
!!!   type(domain2D) :: Domain, Domain_nohalo
end type

!    is   = first x-axis index in the compute domain
!    ie   = last  x-axis index in the compute domain
!    js   = first y-axis index in the compute domain
!    je   = last  y-axis index in the compute domain
!    isd  = first x-axis index in the data    domain
!    ied  = last  x-axis index in the data    domain
!    jsd  = first y-axis index in the data    domain
!    jed  = last  y-axis index in the data    domain
!    isg  = first x-axis index in the global  domain
!    ieg  = last  x-axis index in the global  domain
!    jsg  = first y-axis index in the global  domain
!    jeg  = last  y-axis index in the global  domain
!
!    dx    = grid spacing for x-axis (in meters)
!    rdx   = reciprocal of dx (1/m)
!    dy    = grid spacing for y-axis (in meters)
!    rdy   = reciprocal of dy (1/m)
!    area  = area of a grid box (in m2)
!    rarea = reciprocal of area (1/m2)
!
!    tph  = latitude at the center of grid box (in radians)
!    tlm  = longitude at the center of grid box (in radians)
!    aph  = actual latitude at the center of grid box (in radians)
!    alm  = actual longitude at the center of grid box (in radians)
!
!    blong = longitude grid box boundaries along the global x/longitude axis (in radians)
!    blatg = latitude  grid box boundaries along the global y/latitude  axis (in radians)
!
!    Domain        = domain2D variable with halos of size ihalo,jhalo
!    Domain_nohalo = domain2D variable without halos
!                    used for outputing diagnostic fields
!
!-----------------------------------------------------------------------

type horiz_grid_type
   type(bgrid_type) :: Tmp, Vel

   integer :: nlon, nlat, isize, jsize
   integer :: ilb, iub, jlb, jub, ihalo, jhalo
   logical :: channel, decompx
   real(r8)    :: dlmd, dphd, dlm, dph

   real(r8), pointer, dimension(:,:) :: sinphv, tanphv
end type horiz_grid_type

!    Tmp = grid constants for the temperature/tracer/mass grid
!    Vel = grid constants for the u/v wind component grid
!
!    sinphv = sine of Vel%aph
!    tanphv = tangent of Vel%aph
!
!    nlon = number of grid points along the global x-axis (excluding halo points)
!    nlon = number of grid points along the global y-axis (excluding halo points)
!
!    isize = total number of grid points along the x-axis for the local domain
!             (including halo points)
!    jsize = total number of grid points along the y-axis for the local domain
!             (including halo points)
!
!    ilb  = lower bound x-axis
!    iub  = upper bound x-axis
!    jlb  = lower bound y-axis
!    jub  = upper bound y-axis
!
!    ihalo = number of halo points along the east and west boundaries
!    jhalo = number of halo points along the south and north boundaries
!
!    dlm  = grid spacing for x-axis (in radians)
!    dph  = grid spacing for y-axis (in radians)
!    dlmd = grid spacing for x-axis (in degrees of longitude)
!    dphd = grid spacing for y-axis (in degrees of latitude)
!
!-----------------------------------------------------------------------
!-------- public parameters -------------

   integer, parameter :: TGRID = 51, VGRID = 52

!-----------------------------------------------------------------------
!-------- private parameters ------------

   real(r8),    parameter :: eps=0.0001

   real(r8)    :: lon_beg = 0.0

!-----------------------------------------------------------------------

contains

!#######################################################################

function horiz_grid_init (nlon, nlat, ihalo, jhalo, decomp,  &
                          channel,  tph0d, tlm0d)  result (Hgrid)

!-----------------------------------------------------------------------
!
!    all arguments are input only
!    ----------------------------
!
!    nlon, nlat  = number of grid points in a global 2-d field along
!                  the x (longitude) and y (latitude) axes, respectively
!    ihalo       = number of halo points along the west and east 
!                  boundaries (default = 1)
!    jhalo       = number of halo points along the south and north 
!                  boundaries (default = 1)
!    decomp      = domain decomposition (num X pes by num Y pes),
!                  default decomposition along y-axis then x-axis
!    tph0d,tlm0d = grid/globe transformation
!                  (set both to zero for no transformation, the default)
!

     integer, intent (in)           :: nlon, nlat
     integer, intent (in), optional :: ihalo, jhalo, decomp(2)
     logical, intent (in), optional :: channel
     real(r8)   , intent (in), optional :: tph0d, tlm0d

     type(horiz_grid_type), target :: Hgrid

!-----------------------------------------------------------------------
!------------------- local/private declarations ------------------------

real(r8),    allocatable :: tlmi(:), tphj(:), slat(:), dxj(:)
integer, allocatable :: xrows(:), yrows(:)

real(r8)    :: tph0d_local, tlm0d_local
real(r8)    :: hpi, dtr
integer :: i, j, ilb, iub, jlb, jub, isize, jsize, npes, pe
integer :: is, ie, hs, he, vs, ve
integer :: isd, ied, hsd, hed, vsd, ved
integer :: isg, ieg, hsg, heg, vsg, veg
integer :: domain_layout(2)

!-----------------------------------------------------------------------

                            Hgrid % channel = .false.
      if (present(channel)) Hgrid % channel = channel

!------------- halo region ---------------------------------------------

      Hgrid % ihalo = 1; if (present(ihalo)) Hgrid % ihalo = ihalo
      Hgrid % jhalo = 1; if (present(jhalo)) Hgrid % jhalo = jhalo

      if (Hgrid % ihalo < 1 .or. Hgrid % jhalo < 1) then
        call error_mesg ('horiz_grid_init', 'negative halo size', FATAL)
      endif

!------------- parallel interface --------------------------------------
!        ---- domain decomposition -----

!      call mpp_domains_init

      npes = 1

 
!---- set-up x- & y-axis decomposition -----

 domain_layout = (/ 0, 0 /)
 if (present(decomp)) domain_layout = decomp

 call set_domain_layout ( npes, nlon, nlat, domain_layout )

! flag to indicate x-axis (2d) decomposition
      if ( domain_layout(1) > 1 ) then
         Hgrid % decompx = .true.
      else
         Hgrid % decompx = .false.
      endif


!    ---- mass/temperature grid domain with and without halos ----

!     call mpp_define_domains ( (/1,nlon,1,nlat/), domain_layout,   &
!                               Hgrid % Tmp % Domain,               &
!                               xflags = CYCLIC_GLOBAL_DOMAIN,      &
!                               xhalo = Hgrid % ihalo,              &
!                               yhalo = Hgrid % jhalo               )
!     call mpp_define_domains ( (/1,nlon,1,nlat/), domain_layout,   &
!                               Hgrid % Tmp % Domain_nohalo,        &
!                               xflags = CYCLIC_GLOBAL_DOMAIN,      &
!                               xhalo = 0, yhalo = 0                )

!   ---- compute exact decomposition ----
!   ---- compute 2d layout of PEs ----

     allocate ( xrows(domain_layout(1)), yrows(domain_layout(2)) )
!     call compute_xy_extent ( npes, Hgrid%Tmp%Domain, xrows, yrows )
! For single pe
xrows = nlon
yrows = nlat

!    ---- velocity grid has one less latitude row ----

     yrows(domain_layout(2)) = yrows(domain_layout(2)) - 1

!    ---- velocity grid domain with and without halos ----

!     call mpp_define_domains ( (/1,nlon,1,nlat-1/), domain_layout, &
!                               Hgrid % Vel % Domain,               &
!                               xflags = CYCLIC_GLOBAL_DOMAIN,      &
!                               xhalo = Hgrid % ihalo,              &
!                               yhalo = Hgrid % jhalo,              &
!                               xextent = xrows, yextent = yrows    )
!     call mpp_define_domains ( (/1,nlon,1,nlat-1/), domain_layout, &
!                               Hgrid % Vel % Domain_nohalo,        &
!                               xflags = CYCLIC_GLOBAL_DOMAIN,      &
!                               xhalo = 0, yhalo = 0,               &
!                               xextent = xrows, yextent = yrows    )

     deallocate ( xrows, yrows )


!------------- indices for global compute domain -----------------------

!!!     call mpp_get_global_domain ( Hgrid%Tmp%Domain, isg, ieg, hsg, heg )
!!!     call mpp_get_global_domain ( Hgrid%Vel%Domain, isg, ieg, vsg, veg )

! Compute global domain for 1 process
     isg = 1; ieg = nlon
     hsg = 1; heg = nlat
     vsg = 1; veg = nlat - 1

     Hgrid % Tmp % isg = isg;   Hgrid % Tmp % ieg = ieg
     Hgrid % Tmp % jsg = hsg;   Hgrid % Tmp % jeg = heg
     Hgrid % Vel % isg = isg;   Hgrid % Vel % ieg = ieg
     Hgrid % Vel % jsg = vsg;   Hgrid % Vel % jeg = veg

!------------- indices for data domain -----------------------

!!!     call mpp_get_data_domain ( Hgrid%Tmp%Domain, isd, ied, hsd, hed )
!!!     call mpp_get_data_domain ( Hgrid%Vel%Domain, isd, ied, vsd, ved )

! Compute data domain for 1 process
isd = isg - Hgrid%ihalo; ied = ieg + Hgrid%ihalo
hsd = hsg - Hgrid%jhalo; hed = heg + Hgrid%jhalo
vsd = vsg - Hgrid%jhalo; ved = veg + Hgrid%jhalo

     Hgrid % Tmp % isd = isd;   Hgrid % Tmp % ied = ied
     Hgrid % Tmp % jsd = hsd;   Hgrid % Tmp % jed = hed
     Hgrid % Vel % isd = isd;   Hgrid % Vel % ied = ied
     Hgrid % Vel % jsd = vsd;   Hgrid % Vel % jed = ved

!------------- indices for computational domain ------------------------

!!!     call mpp_get_compute_domain ( Hgrid%Tmp%Domain, is, ie, hs, he )
!!!     call mpp_get_compute_domain ( Hgrid%Vel%Domain, is, ie, vs, ve )

is = isg; ie = ieg; hs = hsg; he = heg
vs = vsg; ve = veg

! Compute domain is same as global domain for 1 process
     Hgrid % Tmp % is = isg;   Hgrid % Tmp % ie = ieg
     Hgrid % Tmp % js = hsg;   Hgrid % Tmp % je = heg
     Hgrid % Vel % is = isg;   Hgrid % Vel % ie = ieg
     Hgrid % Vel % js = vsg;   Hgrid % Vel % je = veg

!------------- indices including halo regions --------------------------

!!!     call mpp_get_data_domain ( Hgrid%Tmp%Domain, ilb, iub, jlb, jub )

! For one pe, just use TMP data domain bounds
     ilb = isd;   iub = ied
     jlb = hsd;   jub = hed

     Hgrid % ilb = ilb;  Hgrid % iub = iub
     Hgrid % jlb = jlb;  Hgrid % jub = jub

!------------ other old values -----------------------------------------

     Hgrid % nlon = nlon
     Hgrid % nlat = nlat

     isize = Hgrid % iub - Hgrid % ilb + 1
     jsize = Hgrid % jub - Hgrid % jlb + 1

     Hgrid % isize = isize
     Hgrid % jsize = jsize

!-----------------------------------------------------------------------
!   check to make sure halos are not too big for domain decomposition
!   this is only a problem for polar halos
!   we do not want the polar halos to go off processor
!   only need to check the mass grid

  if (update_np(Hgrid,TGRID) .or. update_sp(Hgrid,TGRID)) then
     if (Hgrid%jhalo > Hgrid%Tmp%je-Hgrid%Tmp%js+1) call error_mesg &
     ('horiz_grid_init', 'j halo size too big for decomposition', FATAL)
  endif

!-----------------------------------------------------------------------
!    -------- allocate space for arrays ---------

    call alloc_array_space ( ilb, iub, jlb, jub, Hgrid%Tmp )
    call alloc_array_space ( ilb, iub, jlb, jub, Hgrid%Vel )
                 
    allocate ( Hgrid % sinphv(ilb:iub,jlb:jub), &
               Hgrid % tanphv(ilb:iub,jlb:jub)  )

    allocate ( tlmi (is-1:ie+2), tphj (vs-1:ve+1), &
               slat (vs-1:ve+2), dxj  (vs-1:ve+1)  ) 

!-----------------------------------------------------------------------
!--------------derived geometrical constants----------------------------

      Hgrid % dlmd = 360./float(nlon)
      Hgrid % dphd = 180./float(nlat)

      hpi = acos(0.0)
      dtr = hpi/90.
      Hgrid % dlm = Hgrid % dlmd*dtr
      Hgrid % dph = Hgrid % dphd*dtr

!     --- dy is the same on both grids ---

      Hgrid % Tmp %  dy = RADIUS * Hgrid % dph
      Hgrid % Vel %  dy = RADIUS * Hgrid % dph
      Hgrid % Tmp % rdy = 1.0 / Hgrid % Tmp % dy
      Hgrid % Vel % rdy = 1.0 / Hgrid % Vel % dy

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!---- get grid box boundaries for computational global domain ----

     Hgrid % Tmp % blong(isg-1:ieg+2) = get_grid (nlon+3, lon_beg-Hgrid%dlm, &
                                                                  Hgrid%dlm  )
     Hgrid % Tmp % blatg(hsg)   = -hpi
     Hgrid % Tmp % blatg(heg+1) =  hpi
     Hgrid % Tmp % blatg(hsg+1:heg) = get_grid (nlat-1, -hpi+Hgrid%dph, &
                                                             Hgrid%dph  )

     tlmi(is-1:ie+1) = 0.5*(Hgrid%Tmp%blong(is-1:ie+1)+Hgrid%Tmp%blong(is  :ie+2))
     tphj(hs  :he  ) = 0.5*(Hgrid%Tmp%blatg(hs  :he  )+Hgrid%Tmp%blatg(hs+1:he+1))

     if (Hgrid % channel) then
!    --- channel model ---
         dxj = RADIUS * Hgrid % dlm
     else
!    --- sphere ---
         slat(hs:he+1) = sin(Hgrid%Tmp%blatg(hs:he+1))
         dxj (hs:he)   = RADIUS * Hgrid % dlm / Hgrid % dph *  &
                                  (slat(hs+1:he+1)-slat(hs:he))
     endif

     Hgrid % Tmp % tlm = 0.0;  Hgrid % Tmp % tph = 0.0;  Hgrid % Tmp % dx = 0.0

     Hgrid % Tmp % tlm (is-1:ie+1,:) = spread (tlmi(is-1:ie+1), 2, jsize)
     Hgrid % Tmp % tph (:,hs:he)     = spread (tphj(hs  :he  ), 1, isize)
     Hgrid % Tmp % dx  (:,hs:he)     = spread (dxj (hs  :he  ), 1, isize)

!-------------initialize lat/lon at velocity points---------------------

     Hgrid % Vel % blong(isg-1:ieg+2) = get_grid (nlon+3,            &
                                              lon_beg-0.5*Hgrid%dlm, &
                                                          Hgrid%dlm  )
     Hgrid % Vel % blatg(vsg-1) = -hpi
     Hgrid % Vel % blatg(veg+2) =  hpi
     Hgrid % Vel % blatg(vsg:veg+1) = get_grid (nlat, -hpi+0.5*Hgrid%dph, &
                                                               Hgrid%dph  )

     tlmi(is-1:ie+1) = 0.5*(Hgrid%Vel%blong(is-1:ie+1)+Hgrid%Vel%blong(is:ie+2))
     tphj(vs-1:ve+1) = 0.5*(Hgrid%Vel%blatg(vs-1:ve+1)+Hgrid%Vel%blatg(vs:ve+2))

     if (Hgrid % channel) then
!    --- channel model ---
         dxj = RADIUS * Hgrid % dlm
     else
!    --- sphere ---
         slat(vs-1:ve+2) = sin(Hgrid%Vel%blatg(vs-1:ve+2))
         dxj (vs-1:ve+1) = RADIUS * Hgrid % dlm / Hgrid % dph *  &
                                   (slat(vs:ve+2)-slat(vs-1:ve+1))
     endif

     Hgrid % Vel % tlm = 0.0;  Hgrid % Vel % tph = 0.0
     Hgrid % Vel % dx  = 0.0;  Hgrid % Vel % rdx = 0.0

     Hgrid % Vel % tlm (is-1:ie+1,:) = spread (tlmi(is-1:ie+1), 2, jsize)
     Hgrid % Vel % tph (:,vs-1:ve+1) = spread (tphj(vs-1:ve+1), 1, isize)
     Hgrid % Vel % dx  (:,vs-1:ve+1) = spread (dxj (vs-1:ve+1), 1, isize)

     Hgrid % Vel % rdx(:,vs-1:ve+1) = 1.0 / Hgrid % Vel % dx(:,vs-1:ve+1)

     deallocate ( tlmi , tphj , slat , dxj )

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!------------- grid box areas ------------------------------------------

      Hgrid % Tmp % area = Hgrid % Tmp % dx * Hgrid % Tmp % dy
      Hgrid % Vel % area = Hgrid % Vel % dx * Hgrid % Vel % dy

!--- area of velocity pole (time 2 for computational purposes) ----

      if (vs-1 == vsg-1) &
      Hgrid % Vel % area (:,vs-1) = Hgrid % Vel % area (:,vs-1) * 2.0
      if (ve+1 == veg+1) &
      Hgrid % Vel % area (:,ve+1) = Hgrid % Vel % area (:,ve+1) * 2.0

!--- reciprocal of area ----

      where (Hgrid % Tmp % area > 0.0)
         Hgrid % Tmp % rarea = 1.0 / Hgrid % Tmp % area
      elsewhere
         Hgrid % Tmp % rarea = 0.0
      endwhere

      where (Hgrid % Vel % area > 0.0)
         Hgrid % Vel % rarea = 1.0 / Hgrid % Vel % area
      elsewhere
         Hgrid % Vel % rarea = 0.0
      endwhere

!-------------------compute "actual" lat/lon----------------------------

      tph0d_local = 0.0; if (present(tph0d)) tph0d_local = tph0d
      tlm0d_local = 0.0; if (present(tlm0d)) tlm0d_local = tlm0d

      if (tph0d_local > eps .or. tlm0d_local > eps) then
           call trans_latlon (tph0d_local,  tlm0d_local,            &
                              Hgrid % Tmp % tlm, Hgrid % Tmp % tph, &
                              Hgrid % Tmp % alm, Hgrid % Tmp % aph)
           call trans_latlon (tph0d_local,  tlm0d_local,            &
                              Hgrid % Vel % tlm, Hgrid % Vel % tph, &
                              Hgrid % Vel % alm, Hgrid % Vel % aph)
      else
           Hgrid % Tmp % aph = Hgrid % Tmp % tph
           Hgrid % Tmp % alm = Hgrid % Tmp % tlm
           Hgrid % Vel % aph = Hgrid % Vel % tph
           Hgrid % Vel % alm = Hgrid % Vel % tlm
      endif

!------------- trigometric constants -----------------------------------

      if (Hgrid % channel) then
!     --- channel model ---
          Hgrid % sinphv = 0.7071
          Hgrid % tanphv = 0.0
      else
!     --- sphere ---
          Hgrid % sinphv = sin(Hgrid % Vel % aph)
          Hgrid % tanphv = tan(Hgrid % Vel % aph)
      endif

!-----------------------------------------------------------------------
!----- write (to standard output?) domain decomposition -------

!      call print_decomp ( Hgrid%Tmp%Domain, npes )

!-----------------------------------------------------------------------

end function horiz_grid_init

!#######################################################################

function get_grid (npts, start, space) result (grid)

integer, intent(in) :: npts
real(r8),    intent(in) :: start, space
real(r8)                :: grid(npts)
integer :: j

!---- compute equally spaced grid ----

      do j = 1, npts
         grid(j) = start + (j-1.0)*space
      enddo

end function get_grid

!#######################################################################

   subroutine trans_latlon (tph0d, tlm0d, tlm, tph, alm, aph)

   real(r8), intent(in)  :: tph0d, tlm0d, tlm(:,:), tph(:,:)
   real(r8), intent(out) :: alm(:,:), aph(:,:)

   real(r8) :: dtr, pie, tph0, tlm0, stph0, ctph0, ttph0
   real(r8), dimension(size(tlm,1),size(tlm,2)) :: stph, cc, ee

!-----------------------------------------------------------------------

      pie  = 2.*acos(0.0)
      dtr  = pie/180.
      tph0 = tph0d*dtr;  tlm0 = tlm0d*dtr
      stph0 = sin(tph0); ctph0 = cos(tph0)
      ttph0 = stph0/ctph0

!     ------- compute actual lat and lon -------
      stph = sin(tph)
      cc = cos(tph)*cos(tlm)
      aph = asin(ctph0*stph+stph0*cc)

      ee = cc/(ctph0*cos(aph)) - tan(aph)*ttph0
      where (ee > 1.0)  ee = 1.0
      where (tlm > pie)
         alm = tlm0-acos(ee)
      elsewhere
         alm = tlm0+acos(ee)
      endwhere

!-----------------------------------------------------------------------

end subroutine trans_latlon

!#######################################################################

subroutine get_horiz_grid_bound ( Hgrid, grid, blon, blat, global )

   type (horiz_grid_type), intent(in)  :: Hgrid
   integer,                intent(in)  ::  grid
   real(r8),                   intent(out) :: blon(:), blat(:)
   logical, optional,      intent(in)  :: global

!       returns the grid box boundaries for either
!           the compute or global grid

    select case (grid)
       case (TGRID)
          call horiz_grid_bound ( Hgrid%Tmp, blon, blat, global )
       case (VGRID)
          call horiz_grid_bound ( Hgrid%Vel, blon, blat, global )
       case default
          call error_mesg ('get_horiz_grid_bound', 'invalid grid', FATAL)
    end select

!-----------------------------------------------------------------------

end subroutine get_horiz_grid_bound

!#######################################################################

subroutine get_horiz_grid_size ( Hgrid, grid, nlon, nlat, global )

   type (horiz_grid_type), intent(in)  :: Hgrid
   integer,                intent(in)  ::  grid
   integer,                intent(out) :: nlon, nlat
   logical, optional,      intent(in)  :: global

!       returns the number of longitude and latitude grid boxes
!              for either the compute or global grid

    select case (grid)
       case (TGRID)
          call horiz_grid_size ( Hgrid%Tmp, nlon, nlat, global )
       case (VGRID)
          call horiz_grid_size ( Hgrid%Vel, nlon, nlat, global )
       case default
          call error_mesg ('get_horiz_grid_size', 'invalid grid', FATAL)
    end select

end subroutine get_horiz_grid_size

!#######################################################################

subroutine horiz_grid_bound ( Grid, blon, blat, global )

   type (bgrid_type), intent(in)  :: Grid
   real(r8),              intent(out) :: blon(:), blat(:) 
   logical, optional, intent(in)  :: global

!      private routine that returns the grid box boundaries
!          for either the compute or global grid

    integer :: is, ie, ip, js, je, jp
    logical :: lglobal 

    lglobal = .false.;  if (present(global)) lglobal = global

!---------- define longitudinal grid box edges ------------

    if (lglobal) then
      is = Grid % isg; ie = Grid % ieg
    else    
      is = Grid % is ; ie = Grid % ie
    endif   
      ip = ie-is+1 

      if (size(blon) /= ip+1) call error_mesg  &
                          ('get_horiz_grid_bound', &
                           'invalid argument dimension for blon', FATAL)

      blon = Grid % blong (is:ie+1)

!---------- define latitudinal grid box edges ------------

    if (lglobal) then
      js = Grid % jsg; je = Grid % jeg
    else    
      js = Grid % js ; je = Grid % je
    endif   
      jp = je-js+1 

      if (size(blat) /= jp+1) call error_mesg  &
                          ('get_horiz_grid_bound', &
                           'invalid argument dimension for blat', FATAL)

      blat = Grid %  blatg (js:je+1)

!-----------------------------------------------------------------------

end subroutine horiz_grid_bound

!#######################################################################

subroutine horiz_grid_size ( Grid, nlon, nlat, global )

   type (bgrid_type), intent(in)  :: Grid
   integer,                intent(out) :: nlon, nlat
   logical, optional,      intent(in)  :: global

   logical :: lglobal

   lglobal = .false.;  if (present(global)) lglobal = global

!---- return the size of the grid used for physics computations ----

   if (lglobal) then
      nlon = Grid % ieg - Grid % isg + 1
      nlat = Grid % jeg - Grid % jsg + 1
   else
      nlon = Grid % ie - Grid % is + 1
      nlat = Grid % je - Grid % js + 1
   endif

end subroutine horiz_grid_size

!#######################################################################

function update_np (Hgrid, grid) result (answer)

  type (horiz_grid_type), intent(in)  :: Hgrid
  integer,                intent(in)  ::  grid
  logical :: answer
  integer :: je, jeg

     select case (grid)
       case (TGRID)
         je = Hgrid%Tmp%je; jeg=Hgrid%Tmp%jeg
       case(VGRID)
         je = Hgrid%Vel%je; jeg=Hgrid%Vel%jeg
       case default
          call error_mesg ('update_np', 'invalid grid', FATAL)
     end select

     answer = .false.
     if ( je + Hgrid%jhalo > jeg ) answer = .true.

end function update_np

!#######################################################################

function update_sp (Hgrid, grid) result (answer)

  type (horiz_grid_type), intent(in)  :: Hgrid
  integer,                intent(in)  ::  grid
  logical :: answer
  integer :: js, jsg

     select case (grid)
       case (TGRID)
         js = Hgrid%Tmp%js; jsg=Hgrid%Tmp%jsg
       case(VGRID)
         js = Hgrid%Vel%js; jsg=Hgrid%Vel%jsg
       case default
          call error_mesg ('update_sp', 'invalid grid', FATAL)
     end select

     answer = .false.
     if ( js - Hgrid%jhalo < jsg ) answer = .true.

end function update_sp

!#######################################################################

 subroutine alloc_array_space ( ilb, iub, jlb, jub, Grid )
 integer, intent(in) :: ilb, iub, jlb, jub
 type(bgrid_type), intent(inout) :: Grid
         
      allocate ( Grid % dx   (ilb:iub,jlb:jub), &
                 Grid % rdx  (ilb:iub,jlb:jub), &
                 Grid % tph  (ilb:iub,jlb:jub), &
                 Grid % tlm  (ilb:iub,jlb:jub), &
                 Grid % aph  (ilb:iub,jlb:jub), &
                 Grid % alm  (ilb:iub,jlb:jub), &
                 Grid % area (ilb:iub,jlb:jub), &
                 Grid % rarea(ilb:iub,jlb:jub)  )

      allocate ( Grid % blong (Grid%isg-1:Grid%ieg+2), &
                 Grid % blatg (Grid%jsg-1:Grid%jeg+2)  )

 end subroutine alloc_array_space

!#######################################################################

! subroutine compute_xy_extent ( npes, Domain, xrows, yrows )
! integer,        intent(in) :: npes
! type(domain2D), intent(in) :: Domain
! integer       , intent(out) :: xrows(:), yrows(:)
! integer, dimension(0:npes-1) :: xbegin, xend, xsize, &
!                                 ybegin, yend, ysize
! integer :: i, j, pe, m, n, isg, ieg, jsg, jeg
!
!   call mpp_get_global_domain   ( Domain, isg, ieg, jsg, jeg )
!   call mpp_get_compute_domains ( Domain, xbegin, xend, xsize, &
!                                          ybegin, yend, ysize  )
!
!!  compute the exact x-axis decomposition
!   i = isg
!   m = 0
!   xrows = 0
!   do pe = 0, npes-1
!     if ( xbegin(pe) == i ) then
!       m = m+1
!       xrows (m) = xsize (pe)
!       i = i + xsize(pe)
!       if ( m == size(xrows) .or. i > ieg ) exit
!     endif
!   enddo
!
!!  compute the exact y-axis decomposition
!   j = jsg 
!   n = 0   
!   yrows = 0
!   do pe = 0, npes-1
!     if ( ybegin(pe) == j ) then
!       n = n+1 
!       yrows (n) = ysize (pe)
!       j = j + ysize(pe)
!       if ( n == size(yrows) .or. j > jeg ) exit
!     endif   
!   enddo   
!
!
! end subroutine compute_xy_extent

!#######################################################################

 subroutine set_domain_layout ( npes, nlon, nlat, domain_layout )
  integer, intent(in)    :: npes, nlon, nlat
  integer, intent(inout) :: domain_layout(2)

!---- set up number of processors along each axis ----

!----- set-up x- & y-axis decomposition -----
!----- (along y-axis only, then x-axis) -----

!      ****  this scheme will fail if ....  ****
!          the number of PEs (npes) is prime 
!          and npes > 2*nlat and npes > nlon

   if ( domain_layout(1)+domain_layout(2) == 0 ) then
      domain_layout = (/ 1, npes /)
      do while ( 2*domain_layout(2) > nlat )
          do
            domain_layout(1) = domain_layout(1) + 1
            domain_layout(2) = npes/domain_layout(1)
            if ( domain_layout(1)*domain_layout(2) == npes ) exit
          enddo
          if ( domain_layout(1) == nlon ) exit
          if ( domain_layout(1) == npes ) exit
      enddo
   else
     if (domain_layout(1) == 0) domain_layout(1) = npes/domain_layout(2)
     if (domain_layout(2) == 0) domain_layout(2) = npes/domain_layout(1)
   endif

   if ( domain_layout(1)*domain_layout(2) /= npes ) call error_mesg &
              ('horiz_grid_init', 'number of processors requested not &
                                  &compatible with grid', FATAL )

 end subroutine set_domain_layout

!#######################################################################

! subroutine print_decomp ( Domain, npes ) 
! type(domain2D), intent(in) :: Domain
! integer,        intent(in) :: npes
! integer, dimension(0:npes-1) :: xsize, ysize
! integer :: i, j, xlist(npes), ylist(npes), layout(2)
! type (domain1D) :: Xdom, Ydom
!
!   call mpp_get_layout ( Domain, layout )
!   call mpp_get_compute_domains   ( Domain, xsize=xsize, ysize=ysize )
!   call mpp_get_domain_components ( Domain, Xdom, Ydom )
!   call mpp_get_pelist ( Xdom, xlist(1:layout(1)) ) 
!   call mpp_get_pelist ( Ydom, ylist(1:layout(2)) ) 
!
!   write (*,100) 
!   write (*,110) (xsize(xlist(i)),i=1,layout(1))
!   write (*,120) (ysize(ylist(j)),j=1,layout(2))
!
!   100 format ('ATMOSPHERIC (B-GRID) MODEL DOMAIN DECOMPOSITION')
!   110 format ('  X-AXIS = ',24i4,/,(11x,24i4))
!   120 format ('  Y-AXIS = ',24i4,/,(11x,24i4))
!
! end subroutine print_decomp

!#######################################################################

end module bgrid_horiz_mod

