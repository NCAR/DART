! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

PROGRAM add_pert_where_high_refl

! Add 3D random but smooth perturbations to WRF fields in/near where the observations
! indicate high reflectivity values.  The technique is somewhat like that of
! Caya et al. 2005, Monthly Weather Review, 3081-3094.
!
! David Dowell 27 June 2007
!
! input parameters from command line:
! (1) refl_ob_file  -- name of text file containing WRF grid indices where observed reflectivity is high
! (2) wrf_file      -- path name of WRF netcdf file
! (3) lh            -- horizontal length scale (m) for perturbations
! (4) lv            -- vertical length scale (m) for perturbations
! (5) u_sd          -- std. dev. of u noise (m/s), before smoothing
! (6) v_sd          -- std. dev. of v noise (m/s), before smoothing
! (7) w_sd          -- std. dev. of w noise (m/s), before smoothing
! (8) t_sd          -- std. dev. of potential temperature noise (K), before smoothing
! (9) td_sd         -- std. dev. of dewpoint noise (K), before smoothing
! (10) qv_sd        -- std. dev. of water vapor mixing ratio noise, before smoothing
!                         (input value is in g/kg, value after conversion is in kg/kg)
! (11) ens_num      -- ensemble number for consistently seeding the random number generator
! (12) gdays        -- gregorian day number of wrf analysis time
! (13) gsecs        -- gregorian seconds of day of wrf analysis time
!
! output:

use        types_mod, only : r8, gravity, t_kelvin, ps0, gas_constant, gas_constant_v
use    utilities_mod, only : error_handler, E_ERR, initialize_utilities, finalize_utilities
use   random_seq_mod, only : random_gaussian, random_seq_type, init_random_seq
use    netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! command-line parameters
character(len=129) :: refl_ob_file
character(len=129) :: wrf_file
real(r8)           :: lh
real(r8)           :: lv
real(r8)           :: u_sd
real(r8)           :: v_sd
real(r8)           :: w_sd
real(r8)           :: t_sd
real(r8)           :: td_sd
real(r8)           :: qv_sd
integer            :: ens_num
integer            :: gdays
integer            :: gsecs

! local variables
integer               :: n_obs                  ! number of gridded reflectivity observations
integer, allocatable  :: i_ob(:)                ! grid indices of observations
integer, allocatable  :: j_ob(:)                ! "                          "
integer, allocatable  :: k_ob(:)                ! "                          "
real(r8), allocatable :: refl_ob(:)             ! observed reflectivity (dBZ)

real(r8), allocatable :: phb(:,:,:)             ! base-state geopotential (m^2 s^-2)
real(r8), allocatable :: ht(:,:,:)              ! height MSL of mass grid points (m)
real(r8), allocatable :: ht_u(:,:,:)            ! height MSL of u grid points (m)
real(r8), allocatable :: ht_v(:,:,:)            ! height MSL of v grid points (m)
real(r8), allocatable :: ht_w(:,:,:)            ! height MSL of w grid points (m)
real(r8), allocatable :: f(:,:,:)               ! WRF field
real(r8), allocatable :: f2(:,:,:)              ! Extra WRF field
real(r8), allocatable :: sd(:,:,:)              ! standard deviations of grid-point noise

real(r8), allocatable :: dnw(:)                 ! d(eta) values between full (w) levels
real(r8), allocatable :: ph(:,:,:)              ! perturbation geopotential (m^2 s^-2)
real(r8), allocatable :: qv(:,:,:)              ! water vapor (kg/kg)
real(r8), allocatable :: t(:,:,:)               ! perturbation potential temperature (K)
real(r8)              :: rho                    ! density (kg m^-3)
real(r8), allocatable :: mu(:,:)                ! perturbation dry air mass in column (Pa)
real(r8), allocatable :: mub(:,:)               ! base state dry air mass in column (Pa)
real(r8)              :: ph_e
real(r8)              :: qvf1
real(r8), PARAMETER   :: ts0 = 300.0_r8        ! Base potential temperature for all levels.
real(r8), PARAMETER   :: kappa = 2.0_r8/7.0_r8 ! gas_constant / cp
real(r8), PARAMETER   :: rd_over_rv = gas_constant / gas_constant_v
real(r8), PARAMETER   :: cpovcv = 1.4_r8        ! cp / (cp - gas_constant)
real(r8), allocatable :: p(:,:,:)               ! pressure (mb)


real(r8)              :: dx, dy                 ! horizontal grid spacings (m)
integer               :: bt, sn, we             ! WRF grid dimensions
integer               :: i, j, k, o

! netcdf stuff
integer :: var_id, ncid, ierr
character(len=80) :: varname

! f2kcli stuff
integer :: status, length
character(len=120) :: string

! random number generator stuff
type (random_seq_type) :: rs


call initialize_utilities('add_pert_where_high_refl')

 ! Get command-line parameters, using the F2KCLI interface.  See f2kcli.f90 for details.

if( COMMAND_ARGUMENT_COUNT() .ne. 13 ) then
  print*, 'INCORRECT # OF ARGUMENTS ON COMMAND LINE:  ', COMMAND_ARGUMENT_COUNT()
  call exit(1)
else

  call GET_COMMAND_ARGUMENT(1,refl_ob_file,length,status)
  if( status .ne. 0 ) then
    print*, 'refl_ob_file NOT RETRIEVED FROM COMMAND LINE:  ', status
    call exit(1)
  endif

  call GET_COMMAND_ARGUMENT(2,wrf_file,length,status)
  if( status .ne. 0 ) then
    print*, 'wrf_file NOT RETRIEVED FROM COMMAND LINE:  ', status
    call exit(1)
  endif

  call GET_COMMAND_ARGUMENT(3,string,length,status)
  if( status .ne. 0 ) then
    print*,  'lh NOT RETRIEVED FROM COMMAND LINE:  ', status
    call exit(1)
  else
    read(string,*) lh
  endif
  
  call GET_COMMAND_ARGUMENT(4,string,length,status)
  if( status .ne. 0 ) then
    print*,  'lv NOT RETRIEVED FROM COMMAND LINE:  ', status
    call exit(1)
  else
    read(string,*) lv
  endif

  call GET_COMMAND_ARGUMENT(5,string,length,status)
  if( status .ne. 0 ) then
    print*,  'u_sd NOT RETRIEVED FROM COMMAND LINE:  ', status
    call exit(1)
  else
    read(string,*) u_sd
  endif
  
  call GET_COMMAND_ARGUMENT(6,string,length,status)
  if( status .ne. 0 ) then
    print*,  'v_sd NOT RETRIEVED FROM COMMAND LINE:  ', status
    call exit(1)
  else
    read(string,*) v_sd
  endif

  call GET_COMMAND_ARGUMENT(7,string,length,status)
  if( status .ne. 0 ) then
    print*,  'w_sd NOT RETRIEVED FROM COMMAND LINE:  ', status
    call exit(1)
  else
    read(string,*) w_sd
  endif

  call GET_COMMAND_ARGUMENT(8,string,length,status)
  if( status .ne. 0 ) then
    print*,  't_sd NOT RETRIEVED FROM COMMAND LINE:  ', status
    call exit(1)
  else
    read(string,*) t_sd
  endif

  call GET_COMMAND_ARGUMENT(9,string,length,status)
  if( status .ne. 0 ) then
    print*,  'td_sd NOT RETRIEVED FROM COMMAND LINE:  ', status
    call exit(1)
  else
    read(string,*) td_sd
  endif

  call GET_COMMAND_ARGUMENT(10,string,length,status)
  if( status .ne. 0 ) then
    print*,  'qv_sd NOT RETRIEVED FROM COMMAND LINE:  ', status
    call exit(1)
  else
    read(string,*) qv_sd
  endif

  call GET_COMMAND_ARGUMENT(11,string,length,status)
  if( status .ne. 0 ) then
    print*,  'ens_num NOT RETRIEVED FROM COMMAND LINE:  ', status
    call exit(1)
  else
    read(string,*) ens_num
  endif

  call GET_COMMAND_ARGUMENT(12,string,length,status)
  if( status .ne. 0 ) then
    print*,  'gdays NOT RETRIEVED FROM COMMAND LINE:  ', status
    call exit(1)
  else
    read(string,*) gdays
  endif

  call GET_COMMAND_ARGUMENT(13,string,length,status)
  if( status .ne. 0 ) then
    print*,  'gsecs NOT RETRIEVED FROM COMMAND LINE:  ', status
    call exit(1)
  else
    read(string,*) gsecs
  endif

endif

qv_sd = 0.001_r8 * qv_sd    ! convert g/kg to kg/kg


! Read locations where high reflectivity was observed.

! first, count observations

open(unit=11, file=refl_ob_file, status='old')
n_obs = 0
998 read(11,*,end=999)
  n_obs = n_obs + 1
go to 998
999 close(11)

! now allocate storage and read the observations

allocate(i_ob(n_obs))
allocate(j_ob(n_obs))
allocate(k_ob(n_obs))
allocate(refl_ob(n_obs))

open(unit=11, file=refl_ob_file, status='old')
do o=1, n_obs
  read(11,*) i_ob(o), j_ob(o), k_ob(o), refl_ob(o) 
enddo
close(11)


! Open WRF file and obtain miscellaneous values.

call check ( nf90_open(wrf_file, NF90_WRITE, ncid) )

call check ( nf90_inq_dimid(ncid, 'bottom_top', var_id) )
call check ( nf90_inquire_dimension(ncid, var_id, varname, bt) )

call check ( nf90_inq_dimid(ncid, 'south_north', var_id) )
call check ( nf90_inquire_dimension(ncid, var_id, varname, sn) )

call check ( nf90_inq_dimid(ncid, 'west_east', var_id) )
call check ( nf90_inquire_dimension(ncid, var_id, varname, we) )

call check( nf90_get_att(ncid, nf90_global, 'DX', dx) )
call check( nf90_get_att(ncid, nf90_global, 'DY', dy) )


! Read WRF base-state geopotential height field and compute height (m MSL)
! of each grid point.

allocate(phb(we,sn,bt+1))
allocate(ht(we,sn,bt))
allocate(ht_u(we+1,sn,bt))
allocate(ht_v(we,sn+1,bt))
allocate(ht_w(we,sn,bt+1))

call check ( nf90_inq_varid(ncid, 'PHB', var_id))
call check ( nf90_get_var(ncid, var_id, phb, start = (/ 1, 1, 1, 1/)))

do k=1, bt
  do j=1, sn
    do i=1, we
      ht(i,j,k) = ( phb(i,j,k) + phb(i,j,k+1) ) / (2.0_r8*gravity)
    enddo
  enddo
enddo
do k=1, bt
  do j=1, sn
    do i=2, we
      ht_u(i,j,k) = ( phb(i-1,j,k) + phb(i-1,j,k+1) + phb(i,j,k) + phb(i,j,k+1) ) / (4.0_r8*gravity)
    enddo
    ht_u(1,j,k) = ht_u(2,j,k)
    ht_u(we+1,j,k) = ht_u(we,j,k)
  enddo
enddo
do k=1, bt
  do i=1, we
    do j=2, sn
      ht_v(i,j,k) = ( phb(i,j-1,k) + phb(i,j-1,k+1) + phb(i,j,k) + phb(i,j,k+1) ) / (4.0_r8*gravity)
    enddo
    ht_v(i,1,k) = ht_v(i,2,k)
    ht_v(i,sn+1,k) = ht_v(i,sn,k)
  enddo
enddo
do k=1, bt+1
  do j=1, sn
    do i=1, we
      ht_w(i,j,k) = phb(i,j,k) / gravity
    enddo
  enddo
enddo


! Initialize random number generator based on the analysis time and
! the ensemble number so repeated runs have reproducible values.
call init_random_seq(rs, (gdays + gsecs)*1000 + ens_num)


! Add perturbations.

if (u_sd .gt. 0.0_r8) then

  allocate(f(we+1,sn,bt))
  call check ( nf90_inq_varid(ncid, 'U', var_id))
  call check ( nf90_get_var(ncid, var_id, f, start = (/ 1, 1, 1, 1/)))
  allocate(sd(we+1,sn,bt))  

  sd(:,:,:) = 0.0_r8
  do o=1, n_obs
    sd(i_ob(o), j_ob(o), k_ob(o)) = u_sd
    sd(i_ob(o)+1, j_ob(o), k_ob(o)) = u_sd
  enddo

  call add_smooth_perturbations(f, sd, we+1, sn, bt, lh, lv, dx, dy, ht_u)
  
  call check ( nf90_put_var(ncid, var_id, f, start = (/ 1, 1, 1, 1/)))
  deallocate(f)
  deallocate(sd)

end if


if (v_sd .gt. 0.0_r8) then

  allocate(f(we,sn+1,bt))
  call check ( nf90_inq_varid(ncid, 'V', var_id))
  call check ( nf90_get_var(ncid, var_id, f, start = (/ 1, 1, 1, 1/)))
  allocate(sd(we,sn+1,bt))  

  sd(:,:,:) = 0.0_r8
  do o=1, n_obs
    sd(i_ob(o), j_ob(o), k_ob(o)) = v_sd
    sd(i_ob(o), j_ob(o)+1, k_ob(o)) = v_sd
  enddo

  call add_smooth_perturbations(f, sd, we, sn+1, bt, lh, lv, dx, dy, ht_v)

  call check ( nf90_put_var(ncid, var_id, f, start = (/ 1, 1, 1, 1/)))
  deallocate(f)
  deallocate(sd)

end if


! note:  Perturbing w is not advised currently because there is no enforcement
!        that w perturbations should be small near the lower and upper boundaries.
if (w_sd .gt. 0.0_r8) then

  allocate(f(we,sn,bt+1))
  call check ( nf90_inq_varid(ncid, 'W', var_id))
  call check ( nf90_get_var(ncid, var_id, f, start = (/ 1, 1, 1, 1/)))
  allocate(sd(we,sn,bt+1))

  sd(:,:,:) = 0.0_r8
  do o=1, n_obs
    sd(i_ob(o), j_ob(o), k_ob(o)) = w_sd
    sd(i_ob(o), j_ob(o), k_ob(o)+1) = w_sd
  enddo

  call add_smooth_perturbations(f, sd, we, sn, bt+1, lh, lv, dx, dy, ht_w)

  call check ( nf90_put_var(ncid, var_id, f, start = (/ 1, 1, 1, 1/)))
  deallocate(f)
  deallocate(sd)

end if


if (t_sd .gt. 0.0_r8) then

  allocate(f(we,sn,bt))
  call check ( nf90_inq_varid(ncid, 'T', var_id))
  call check ( nf90_get_var(ncid, var_id, f, start = (/ 1, 1, 1, 1/)))
  allocate(sd(we,sn,bt))

  sd(:,:,:) = 0.0_r8
  do o=1, n_obs
    sd(i_ob(o), j_ob(o), k_ob(o)) = t_sd
  enddo

  call add_smooth_perturbations(f, sd, we, sn, bt, lh, lv, dx, dy, ht)

  call check ( nf90_put_var(ncid, var_id, f, start = (/ 1, 1, 1, 1/)))
  deallocate(f)
  deallocate(sd)

end if


! note:  Perturbing dewpoint can produce supersaturation.
if (td_sd .gt. 0.0_r8) then

  allocate(dnw(bt))
  allocate(ph(we,sn,bt+1))
  allocate(qv(we,sn,bt))
  allocate(t(we,sn,bt))
  allocate(mu(we,sn))
  allocate(mub(we,sn))
  allocate(p(we,sn,bt))
  call check ( nf90_inq_varid(ncid, 'DNW', var_id))
  call check ( nf90_get_var(ncid, var_id, dnw, start = (/ 1, 1/)))
  call check ( nf90_inq_varid(ncid, 'PH', var_id))
  call check ( nf90_get_var(ncid, var_id, ph, start = (/ 1, 1, 1, 1/)))
  call check ( nf90_inq_varid(ncid, 'QVAPOR', var_id))
  call check ( nf90_get_var(ncid, var_id, qv, start = (/ 1, 1, 1, 1/)))
  call check ( nf90_inq_varid(ncid, 'T', var_id))
  call check ( nf90_get_var(ncid, var_id, t, start = (/ 1, 1, 1, 1/)))
  call check ( nf90_inq_varid(ncid, 'MU', var_id))
  call check ( nf90_get_var(ncid, var_id, mu, start = (/ 1, 1, 1/)))
  call check ( nf90_inq_varid(ncid, 'MUB', var_id))
  call check ( nf90_get_var(ncid, var_id, mub, start = (/ 1, 1, 1/)))
  allocate(f(we,sn,bt))
  allocate(f2(we,sn,bt))
! f2 is the sensible temperature array
  allocate(sd(we,sn,bt))

  ! compute pressure
  ! see model_rho_t and model_pressure_t functions in model_mod.f90
  do k=1, bt
    do j=1, sn
      do i=1, we
        ph_e = ( (ph(i,j,k+1) + phb(i,j,k+1)) - (ph(i,j,k) + phb(i,j,k)) ) / dnw(k)
        rho = - ( mub(i,j)+mu(i,j) ) / ph_e
        qvf1 = 1.0_r8 + qv(i,j,k) / rd_over_rv
        p(i,j,k) = 0.01_r8 * ps0 * ( (gas_constant*(ts0+t(i,j,k))*qvf1) / (ps0/rho) )**cpovcv
        f2(i,j,k) = (ts0 + t(i,j,k))*(100.0*p(i,j,k)/ps0)**kappa
      enddo
    enddo
  enddo

  ! compute dewpoint
  do k=1, bt
    do j=1, sn
      do i=1, we
        call compute_td(f(i,j,k), p(i,j,k), qv(i,j,k))
      enddo
    enddo
  enddo

  ! perturb dewpoint
  sd(:,:,:) = 0.0_r8
  do o=1, n_obs
    sd(i_ob(o), j_ob(o), k_ob(o)) = td_sd
  enddo


  call add_smooth_perturbations(f, sd, we, sn, bt, lh, lv, dx, dy, ht)

! check to make sure the perturbed dewpoint is <= temperature + 4 K
  do k=1, bt
    do j=1, sn
      do i=1, we
        if ( f(i,j,k) .gt. f2(i,j,k)+4.0 ) then
!          write(*,*) 'supersaturation violation i,j,k ', i, j, k, f(i,j,k), f2(i,j,k)
          f(i,j,k) = f2(i,j,k)+4.0
        end if
      enddo
    enddo
  enddo

  ! compute qv
  do k=1, bt
    do j=1, sn
      do i=1, we
        call compute_qv(qv(i,j,k), p(i,j,k), f(i,j,k))
      enddo
    enddo
  enddo

  call check ( nf90_inq_varid(ncid, 'QVAPOR', var_id))
  call check ( nf90_put_var(ncid, var_id, qv, start = (/ 1, 1, 1, 1/)))
  deallocate(dnw)
  deallocate(ph)
  deallocate(qv)
  deallocate(t)
  deallocate(mu)
  deallocate(mub)
  deallocate(p)
  deallocate(f)
  deallocate(f2)
  deallocate(sd)

end if


! note:  Since negative qv values are set to 0 after the perturbations are added,
!        the following procedure produces a net increase in qv in the domain.
if (qv_sd .gt. 0.0_r8) then

  allocate(f(we,sn,bt))
  call check ( nf90_inq_varid(ncid, 'QVAPOR', var_id))
  call check ( nf90_get_var(ncid, var_id, f, start = (/ 1, 1, 1, 1/)))
  allocate(sd(we,sn,bt))

  sd(:,:,:) = 0.0_r8
  do o=1, n_obs
    sd(i_ob(o), j_ob(o), k_ob(o)) = qv_sd
  enddo

  call add_smooth_perturbations(f, sd, we, sn, bt, lh, lv, dx, dy, ht)

  f(:,:,:) = max(0.0_r8, f(:,:,:))          ! require qv to be nonnegative

  call check ( nf90_put_var(ncid, var_id, f, start = (/ 1, 1, 1, 1/)))
  deallocate(f)
  deallocate(sd)

end if


! Close file and deallocate arrays.

ierr = NF90_close(ncid)

deallocate(i_ob)
deallocate(j_ob)
deallocate(k_ob)
deallocate(refl_ob)
deallocate(phb)
deallocate(ht)
deallocate(ht_u)
deallocate(ht_v)
deallocate(ht_w)


call finalize_utilities('add_pert_where_high_refl')

contains


  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent ( in) :: istatus 
    if(istatus /= nf90_noerr) call error_handler(E_ERR,'add_pert_where_high_refl', &
         trim(nf90_strerror(istatus)), source, revision, revdate)
  end subroutine check

!------------------------------------------------------------------------------------
  
  ! Compute dewpoint (in Kelvin) from the specified values of pressure and water vapor mixing ratio.
  ! Author:  David Dowell
  ! Date:  July 9, 2007

  subroutine compute_td(td, p, qv)
    implicit none

!-- returned parameter
    real(r8), intent(out) :: td              ! dewpoint (K)

!-- passed parameters
    real(r8), intent(in) :: p                ! pressure (mb)
    real(r8), intent(in) :: qv               ! water vapor mixing ratio (kg/kg)

!-- local variables
    real(r8) :: e                            ! water vapor pressure (mb)
    real(r8) :: qv_ckd                       ! checked mixing ratio
    real(r8), PARAMETER :: e_min = 0.001_r8  ! threshold for minimum vapor pressure (mb),
                                             !   to avoid problems near zero in Bolton's equation
    real(r8), PARAMETER :: qv_min = 1.0e-12  ! threshold for minimum water vapor mixing ratio
    real(r8), PARAMETER :: qv_max = 0.050    ! threshold for maximum water vapor mixing ratio

! ensure qv is a reasonable number
!    qv_ckd = max (qv, qv_min)
!    qv_ckd = min (qv_ckd, qv_max)
    qv_ckd = qv
    e = qv_ckd * p / (0.622_r8 + qv_ckd)                                         ! vapor pressure
    e = max(e, e_min)                                                            ! avoid problems near zero
    td = t_kelvin + (243.5_r8 / ((17.67_r8 / log(e/6.112_r8)) - 1.0_r8) )        ! Bolton's approximation

  end subroutine compute_td

!------------------------------------------------------------------------------------
  
  ! Compute water vapor mixing ratio (in kg/kg) from the specified values of pressure and dewpoint.
  ! Author:  David Dowell
  ! Date:  July 9, 2007

  subroutine compute_qv(qv, p, td)
    implicit none

!-- returned parameter
    real(r8), intent(out) :: qv           ! water vapor mixing ratio (kg/kg)

!-- passed parameters
    real(r8), intent(in) :: p             ! pressure (mb)
    real(r8), intent(in) :: td            ! dewpoint (K)

!-- local variables
    real(r8) :: tdc                       ! dewpoint (Celsius)
    real(r8) :: e                         ! water vapor pressure (mb)

    tdc = td - t_kelvin
    e = 6.112_r8 * exp(17.67_r8 * tdc / (tdc+243.5_r8) )       ! Bolton's approximation
    qv = 0.622_r8 * e / (p-e)

    return

  end subroutine compute_qv

!------------------------------------------------------------------------------------

  ! Add smooth perturbations to an array.  The technique is based on
  ! Caya et al. 2005, Monthly Weather Review, 3081-3094.
  ! Author:  David Dowell
  ! Date:  July 9, 2007

  subroutine add_smooth_perturbations(f, sd, nx, ny, nz, lh, lv, dx, dy, ht)
    implicit none

!-- passed parameters
    integer, intent(in) :: nx, ny, nz        ! grid dimensions
    real(r8), intent(in) :: lh, lv           ! horizontal and vertical length scales (m)
    real(r8), intent(in) :: dx, dy           ! horizontal grid spacings (m)
    real(r8), intent(in) :: ht(nx,ny,nz)     ! heights MSL of f and sd grid points (m)
    real(r8), intent(in) :: sd(nx,ny,nz)     ! standard deviation of grid-point noise

!-- passed and returned variable
    real(r8), intent(inout) :: f(nx,ny,nz)   ! field to be perturbed

!-- local variables

    real(r8) :: r(nx,ny,nz)                  ! realization of random, normally distributed noise
    integer :: i, i0, j, j0, k, k0           ! grid indices
    integer :: i1, i2, j1, j2, k1, k2        ! more grid indices
    real(r8) :: rlh, rlv                     ! reciprocals of lh and lv
    integer, parameter :: nl = 5             ! number of length scales for computing exponential weight


    rlh = 1.0_r8 / lh
    rlv = 1.0_r8 / lv

!   generate random, normally-distributed gridpoint noise

    r(:,:,:) = 0.0_r8
    do k0=1, nz
      do j0=1, ny
        do i0=1, nx
          if (sd(i0,j0,k0) .ne. 0.0_r8) then
            r(i0,j0,k0) = random_gaussian(rs, 0.0_r8, sd(i0,j0,k0))
          endif
        enddo
      enddo
    enddo

!   smooth the perturbations with an inverse exponential function

    do k0=1, nz
      do j0=1, ny
        do i0=1, nx

          if (r(i0,j0,k0).ne.0.0_r8) then

            i1 = max(1, nint(i0-nl*lh/dx))
            i2 = min(nx, nint(i0+nl*lh/dx))
            j1 = max(1, nint(j0-nl*lh/dy))
            j2 = min(ny, nint(j0+nl*lh/dy))
            k1 = k0
            do while ( (k1.gt.1) .and. ( ht(i0,j0,k1) .gt. (ht(i0,j0,k0)-nl*lv) ) )
              k1 = k1 - 1
            enddo
            k2 = k0
            do while ( (k2.lt.nz) .and. ( ht(i0,j0,k2) .lt. (ht(i0,j0,k0)+nl*lv) ) )
              k2 = k2 + 1
            enddo

            do k=k1, k2
              do j=j1, j2
                do i=i1, i2
                  f(i,j,k) = f(i,j,k)                                               &
                           + r(i0,j0,k0)*exp( -dx*abs(i0-i)*rlh                     &
                                              -dy*abs(j0-j)*rlh                     &
                                              -abs(ht(i0,j0,k0)-ht(i,j,k))*rlv )
                enddo
              enddo
            enddo

          endif

        enddo
      enddo
    enddo

  end subroutine add_smooth_perturbations


END PROGRAM add_pert_where_high_refl

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
