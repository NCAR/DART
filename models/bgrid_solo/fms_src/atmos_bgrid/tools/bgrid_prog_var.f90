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

module bgrid_prog_var_mod

!-----------------------------------------------------------------------
!
!       allocates storage for the basic dynamical variables
!
!-----------------------------------------------------------------------
!--------- public defined data type prog_var_type ----------------------
!
!     nlon = number of longitude points (first dimension)
!            includes 2 halo points (1 west, 1 east)
!     nlat = number of latitude points (second dimension)
!            includes 3 halo points
!     nlev = number of vertical levels
!
!     ntrace = number of tracers
!
!     u    = zonal wind component
!     v    = meridional wind component
!     t    = temperature
!     r    = arbitrary number of tracers (includes specific humidity)
!
!     ps   = surface pressure
!     pssl = surface pressure adjust to eta=1. (for eta coordinate)
!
!-----------------------------------------------------------------------

use types_mod, only : r8
use      bgrid_horiz_mod, only: horiz_grid_type
use       bgrid_vert_mod, only: vert_grid_type
use       bgrid_halo_mod, only: update_halo, TEMP, UWND, VWND
use bgrid_cold_start_mod, only: cold_start_resol, cold_start
use              fms_mod, only: file_exist, open_restart_file, error_mesg, &
                                FATAL, close_file,  &
                                write_version_number
use    field_manager_mod, only: MODEL_ATMOS
use   tracer_manager_mod, only: get_tracer_names, set_tracer_profile

implicit none
private

public :: prog_var_type, prog_var_init, var_init,  &
          prog_var_time_diff, prog_var_times_scalar, &
          prog_var_equals_scalar
public :: open_prog_var_file, read_prog_var, write_prog_var

! data structure that contains all prognostic fields and tracers
type prog_var_type
     integer       :: nlon, nlat, nlev, ntrace
     integer       :: ilb, iub, jlb, jub, klb, kub
     real(r8), pointer :: ps(:,:), pssl(:,:)
     real(r8), pointer :: u(:,:,:), v(:,:,:), t(:,:,:), r(:,:,:,:)
end type prog_var_type

! overloaded interface for initializing real(r8) model arrays
interface var_init
    module procedure var_init_type_4d, var_init_bound_4d, &
                     var_init_type_3d, var_init_bound_3d, &
                     var_init_type_2d, var_init_bound_2d
end interface

! private data

logical :: do_log = .true.
character(len=128) :: version='$Revision$'
character(len=128) :: tagname='$Id$'

integer :: unit_in
logical :: read_pssl

character(len=64) :: res_file_name = 'bgrid_prog_var.res'

character(len=80) :: restart_format = &
              'bgrid grid atmospheric dynamical core: restart format 05'

contains

!#######################################################################
! creates a prog_var_type variable

 subroutine prog_var_init (Hgrid, nlev, ntrs, Vars)

  type(horiz_grid_type), intent(in)  :: Hgrid
  integer,               intent(in)  :: nlev, ntrs
  type(prog_var_type)  , intent(out) :: Vars
!-----------------------------------------------------------------------
! write version info to logfile
  if (do_log) then
    call write_version_number (version,tagname)
    do_log = .false.
  endif

! all arrays have the same horizontal dimensions regardless of
! whether the field is on the temperature or velocity grid

    Vars % ilb = Hgrid % ilb
    Vars % iub = Hgrid % iub
    Vars % jlb = Hgrid % jlb
    Vars % jub = Hgrid % jub
    Vars % klb = 1
    Vars % kub = nlev

    Vars % nlon = Hgrid % nlon
    Vars % nlat = Hgrid % nlat
    Vars % nlev = nlev
    Vars % ntrace = ntrs

    Vars % ps   => var_init_bound_2d (Vars % ilb, Vars % iub, &
                                      Vars % jlb, Vars % jub)

    Vars % pssl => var_init_bound_2d (Vars % ilb, Vars % iub, &
                                      Vars % jlb, Vars % jub)

    Vars % u => var_init_bound_3d (Vars % ilb, Vars % iub, &
                                   Vars % jlb, Vars % jub, &
                                   nlev)
    Vars % v => var_init_bound_3d (Vars % ilb, Vars % iub, &
                                   Vars % jlb, Vars % jub, &
                                   nlev)
    Vars % t => var_init_bound_3d (Vars % ilb, Vars % iub, &
                                   Vars % jlb, Vars % jub, &
                                   nlev)
    Vars % r => var_init_bound_4d (Vars % ilb, Vars % iub, &
                                   Vars % jlb, Vars % jub, &
                                   nlev, ntrs)

 end subroutine prog_var_init

!#######################################################################
!##### overloaded functions that allocate a single real(r8) variable #######
!#######################################################################
!
!      variables must be declard as pointers
!      real(r8), pointer :: field(:,:,:)
!      field => var_init (Hgrid,nlev)
!
!#######################################################################

 function var_init_bound_2d (ilb, iub, jlb, jub) result (var)

  integer, intent(in)           :: ilb, iub, jlb, jub
  real(r8), dimension(:,:), pointer :: var

    allocate ( var (ilb:iub, jlb:jub) )
    var = 0.0

 end function var_init_bound_2d

!#######################################################################

 function var_init_type_2d (Hgrid) result (var)

  type(horiz_grid_type), intent(in) :: Hgrid
  real(r8), dimension(:,:), pointer     :: var

    var => var_init_bound_2d (Hgrid % ilb, Hgrid % iub, &
                              Hgrid % jlb, Hgrid % jub)

 end function var_init_type_2d

!#######################################################################

 function var_init_bound_3d (ilb, iub, jlb, jub, kdim) result (var)

  integer, intent(in)             :: ilb, iub, jlb, jub, kdim
  real(r8), dimension(:,:,:), pointer :: var

    allocate ( var (ilb:iub, jlb:jub, 1:kdim) )
    var = 0.0

 end function var_init_bound_3d

!#######################################################################

 function var_init_type_3d (Hgrid, kdim) result (var)

  type(horiz_grid_type), intent(in) :: Hgrid
  integer, intent(in)               :: kdim
  real(r8), dimension(:,:,:), pointer   :: var

    var => var_init_bound_3d (Hgrid % ilb, Hgrid % iub, &
                              Hgrid % jlb, Hgrid % jub, kdim)

 end function var_init_type_3d

!#######################################################################

 function var_init_bound_4d (ilb, iub, jlb, jub, kdim, ntrace) result (var)

  integer, intent(in)               :: ilb, iub, jlb, jub, kdim, ntrace
  real(r8), dimension(:,:,:,:), pointer :: var

    allocate ( var (ilb:iub, jlb:jub, 1:kdim, 1:ntrace) )
    var = 0.0

 end function var_init_bound_4d

!#######################################################################

 function var_init_type_4d (Hgrid, kdim, ntrace) result (var)

  type(horiz_grid_type), intent(in)   :: Hgrid
  integer, intent(in)                 :: kdim, ntrace
  real(r8), dimension(:,:,:,:), pointer   :: var

    var => var_init_bound_4d (Hgrid % ilb, Hgrid % iub, &
                              Hgrid % jlb, Hgrid % jub, kdim, ntrace)

 end function var_init_type_4d

!#######################################################################
!#######################################################################
! sets all prognostic variables to a scalar

 subroutine prog_var_equals_scalar (Var, scalar)

  type(prog_var_type), intent(inout) :: Var
  real(r8)               , intent(in)    :: scalar

     Var % u    = scalar
     Var % v    = scalar
     Var % t    = scalar
     Var % r    = scalar
     Var % ps   = scalar
     Var % pssl = scalar

 end subroutine prog_var_equals_scalar

!#######################################################################
! multiplies all prognostic variables by a scalar

 subroutine prog_var_times_scalar (Var, scalar)

  type(prog_var_type), intent(inout) :: Var
  real(r8)               , intent(in)    :: scalar

     Var % u    = Var % u    * scalar
     Var % v    = Var % v    * scalar
     Var % t    = Var % t    * scalar
     Var % r    = Var % r    * scalar
     Var % ps   = Var % ps   * scalar
     Var % pssl = Var % pssl * scalar

 end subroutine prog_var_times_scalar

!#######################################################################
! performs time differencing on all prognostic variables
!             Var = Var + dt * Var_dt
! all tracers are used unless argument nt is supplied

 subroutine prog_var_time_diff (dt, Var_dt, Var, nt)

  real(r8),                intent(in)    :: dt
  type(prog_var_type), intent(inout) :: Var_dt, Var
  integer, optional,   intent(in)    :: nt

  integer :: ntp

!----- explicit differencing with two time levels -----

   ntp = Var_dt % ntrace
   if (present(nt)) ntp = min(Var_dt%ntrace, nt)

   Var % ps   = Var % ps   + dt * Var_dt % ps
   Var % pssl = Var % pssl + dt * Var_dt % pssl

   Var % u = Var % u + dt * Var_dt % u
   Var % v = Var % v + dt * Var_dt % v
   Var % t = Var % t + dt * Var_dt % t

   Var % r(:,:,:,1:ntp) = Var % r(:,:,:,1:ntp) + &
                                dt * Var_dt % r(:,:,:,1:ntp)

!----- zero out tendencies -----

   Var_dt % ps   = 0.0
   Var_dt % pssl = 0.0

   Var_dt % u = 0.0
   Var_dt % v = 0.0
   Var_dt % t = 0.0

   Var_dt % r(:,:,:,1:ntp) = 0.0


 end subroutine prog_var_time_diff

!#######################################################################
!########## routines for reading and writing restart files #############
!#######################################################################
! this routine returns the model resolution and number of tracers
! reads the control record and resolution of restart file

subroutine open_prog_var_file (ix, jx, kx)

 integer, intent(out) :: ix, jx, kx

 integer :: ic, vers, day, sec, ntsd, nt, ntp
 character(len=80) :: control
 character(len=2) :: avers

! write version info to logfile
  if (do_log) then
    call write_version_number (version,tagname)
    do_log = .false.
  endif

! when restart file does not exist
! set up simple initial conditions

  if ( .not.file_exist('INPUT/'//trim(res_file_name)) ) then
       call cold_start_resol ( ix, jx, kx )
       read_pssl = .false.
       return
  endif

! open restart file and get restart version number
! if control record cannot be read then file uses older format (1 or 2)

  unit_in = open_restart_file ( 'INPUT/'//trim(res_file_name), 'read' )
  read  (unit_in,err=2)  control

! extract version number
  ic = index(control,'restart format ')
  !!!if (ic == 0) call mpp_error ('bgrid_prog_var_mod', &
    !!!           'problem extracting restart version number', FATAL)
  if (ic == 0) call error_mesg ('bgrid_prog_var_mod', &
               'problem extracting restart version number', FATAL)
  avers = control(ic+15:ic+16)
  read (avers,'(i2.2)') vers
  go to 3

! read version number from old format (first rewind file)
2 rewind (unit_in)
  read   (unit_in) vers
  write  (avers,'(i2.2)') vers

! read first (non-control) record of restart file
! note: ntsd,day,sec are no longer used and
!       number of time levels (nvar) is not read or used

3 continue
  select case (vers)
    case (1:2)
         read  (unit_in) ntsd, day, sec, ix, jx, kx, nt, ntp
         read_pssl = .false.
    case (3)
         read  (unit_in) ntsd, day, sec, ix, jx, kx, nt, ntp
         read_pssl = .true.
    case (4)
         read  (unit_in) ix, jx, kx, nt, ntp
         read_pssl = .true.
    case (5)
         read  (unit_in) ix, jx, kx
         read_pssl = .true.
    case default
       !call mpp_error ('bgrid_prog_var_mod', &
       !                'cannot not read restart version '//avers, FATAL)
       call error_mesg ('bgrid_prog_var_mod', &
                       'cannot not read restart version '//avers, FATAL)
  end select

end subroutine open_prog_var_file

!#######################################################################

subroutine read_prog_var (Hgrid, Var, eta, peta, fis, res)

 type(horiz_grid_type), intent(inout) :: Hgrid
 type  (prog_var_type), intent(inout) :: Var
 real(r8), intent(out), dimension(:)      :: eta, peta
 real(r8), intent(out), dimension(Hgrid%ilb:Hgrid%iub, &
                              Hgrid%jlb:Hgrid%jub) :: fis, res
 integer :: n, unit
 integer :: isd,ied,vsd,ved
 character(len=64) :: tr_name
 real(r8) :: tr_surf, tr_mult

  if ( .not.file_exist('INPUT/'//trim(res_file_name)) ) then

     ! set up simple initial conditions
     ! when restart file does not exist

       call cold_start ( Hgrid, eta, peta, fis, res, Var%ps, Var%pssl, &
                                           Var%u, Var%v, Var%t )
  else

     ! must pass fields to read_data on data domain
     ! mass fields are on data domain
     ! set up indexing for velocity fields on data domain

       isd = Hgrid%Vel%isd;  ied = Hgrid%Vel%ied
       vsd = Hgrid%Vel%jsd;  ved = Hgrid%Vel%jed

     ! read non-distributed data from root pe (vertical coordinate info)
       read (unit_in) eta, peta

     !  --- read variables ---

     ! initialize domain for temperature grid 
     ! read surf pres, topog, and more
       !!! call set_domain ( Hgrid%Tmp%Domain )

write(*, *) 'in read prog_var'
write(*, *) 'Read prog_var calls to read_data are not supported in'
write(*, *) 'this non-mpi version of the bgrid model'
write(*, *) 'MOdify as appropriate'
if (1 == 1) stop
       !!!call read_data ( unit_in, Var%ps  )
       if (read_pssl) &
       !!!call read_data ( unit_in, Var%pssl)
       !!!call read_data ( unit_in,     res )
       !!!call read_data ( unit_in,     fis )
     
     ! initialize domain for velocity grid 
     ! read u and v wind components
       !!!call set_domain ( Hgrid%Vel%Domain )

     ! pass velocity fields on data domain
       !!!call read_data ( unit_in, Var%u(isd:ied,vsd:ved,:) )
       !!!call read_data ( unit_in, Var%v(isd:ied,vsd:ved,:) )

     ! re-initialize domain for temperature grid 
     ! read temperature and tracers
       !!! call set_domain ( Hgrid%Tmp%Domain )
       !!!call read_data ( unit_in, Var%t )

     ! done reading B-grid dynamics restart
       call close_file (unit_in)

  endif

! read tracer restart file(s)
  do n = 1, Var%ntrace
     call get_tracer_names ( MODEL_ATMOS, n, tr_name )
     if (file_exist('INPUT/tracer_'//trim(tr_name)//'.res')) then
         unit = open_restart_file( 'INPUT/tracer_'//trim(tr_name)//'.res', 'read' )
         !!!call read_data ( unit, Var%r(:,:,:,n) )
         call close_file (unit)
     else
       ! initialize new tracers (apply surface value only)
         call set_tracer_profile ( MODEL_ATMOS, n, tr_surf, tr_mult )
         Var%r(:,:,:,n) = tr_surf
     endif
  enddo

! update all boundaries for restart variables

  call update_halo (Hgrid, TEMP, res)
  call update_halo (Hgrid, TEMP, fis)
  call update_halo (Hgrid, TEMP, Var%ps)
  if (read_pssl) &
  call update_halo (Hgrid, TEMP, Var%pssl)
  call update_halo (Hgrid, TEMP, Var%t)
  call update_halo (Hgrid, TEMP, Var%r)
  call update_halo (Hgrid, UWND, Var%u)
  call update_halo (Hgrid, VWND, Var%v)

! for old restart formats, initialize pssl
  if (.not.read_pssl) Var%pssl = Var%ps * res


end subroutine read_prog_var

!#######################################################################

 subroutine write_prog_var (Var, Hgrid, Vgrid, fis, res)

 type  (prog_var_type), intent(in) :: Var
 type(horiz_grid_type), intent(in) :: Hgrid
 type (vert_grid_type), intent(in) :: Vgrid
   real(r8), intent(in), dimension(Hgrid%ilb:Hgrid%iub, &
                               Hgrid%jlb:Hgrid%jub) :: fis, res

 integer :: n, unit
 integer :: isd,ied, vsd,ved
 character(len=64) :: rname

!-----------------------------------------------------------------------
! open output restart file


  unit = open_restart_file ( 'RESTART/'//trim(res_file_name), 'write' )

! write non-distributed data from root pe
       write (unit)  restart_format
       write (unit)  Hgrid%nlon, Hgrid%nlat, Vgrid%nlev
       write (unit)  Vgrid%eta, Vgrid%peta ! vertical coordinate info

! must pass fields to write_data on data domain
! mass fields are on data domain
! set up indexing for velocity fields on data domain

  isd = Hgrid%Vel%isd;  ied = Hgrid%Vel%ied
  vsd = Hgrid%Vel%jsd;  ved = Hgrid%Vel%jed

! initialize domain for temperature grid (save surf pres, topog)
  !!! call set_domain ( Hgrid%Tmp%Domain )

write(*, *) 'write_prog_var in bgrid_prog_var is not supported '
write(*, *) 'in this non-mpi bgrid version. Need to modify the'
write(*, *) 'calls to write_data in fms_io'
if (1 == 1) stop


!  call write_data ( unit, Var%ps  )
!  call write_data ( unit, Var%pssl)
!  call write_data ( unit,     res )
!  call write_data ( unit,     fis )

! initialize domain for velocity grid (save u and v components)
  !!! call set_domain ( Hgrid%Vel%Domain )

! pass velocity fields on data domain
!  call write_data ( unit, Var%u(isd:ied,vsd:ved,:) )
!  call write_data ( unit, Var%v(isd:ied,vsd:ved,:) )

! re-initialize domain for temperature grid (save temp and tracers)
  !!! call set_domain ( Hgrid%Tmp%Domain )
!  call write_data ( unit, Var%t   )

! done writing B-grid dynamics restart
  call close_file (unit)

! write tracer restart file(s)
  do n = 1, Var%ntrace
     call get_tracer_names ( MODEL_ATMOS, n, rname )
     unit = open_restart_file( 'RESTART/tracer_'//trim(rname)//'.res', 'write' )
!     call write_data ( unit, Var%r(:,:,:,n) )
     call close_file (unit)
  enddo

 end subroutine write_prog_var

!#######################################################################

end module bgrid_prog_var_mod

