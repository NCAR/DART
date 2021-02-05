! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! BEGIN DART PREPROCESS TYPE DEFINITIONS
! TEMPERATURE,          QTY_TEMPERATURE,         COMMON_CODE
! SURFACE_PRESSURE,     QTY_SURFACE_PRESSURE,    COMMON_CODE
! SKIN_TEMPERATURE,     QTY_SKIN_TEMPERATURE,    COMMON_CODE
! TES_NADIR_OBS,        QTY_NADIR_RADIANCE
! END DART PREPROCESS TYPE DEFINITIONS


! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_TES_nadir_mod, only : write_TES_nadir_obs, read_TES_nadir_obs, &
!                                    interactive_TES_nadir_obs, get_expected_TES_nadir_obs
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE


! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(TES_NADIR_OBS)                                                         
!            call get_expected_TES_nadir_obs(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF


! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(TES_NADIR_OBS)
!         call read_TES_nadir_obs(obs_def%key, ifile, fform)
! END DART PREPROCESS READ_OBS_DEF


! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(TES_NADIR_OBS)
!         call write_TES_nadir_obs(obs_def%key, ifile, fform)
! END DART PREPROCESS WRITE_OBS_DEF


! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(TES_NADIR_OBS)
!         call interactive_TES_nadir_obs(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF


! BEGIN DART PREPROCESS MODULE CODE
module obs_def_TES_nadir_mod

use        types_mod, only : r8, missing_r8, PI, DEG2RAD
use    utilities_mod, only : register_module, error_handler, E_ERR, E_WARN, E_MSG, &
                             check_namelist_read, find_namelist_in_file, &
                             logfileunit, nmlfileunit, do_output, file_exist, &
                             open_file, close_file, get_unit, do_nml_file, do_nml_term, &
                             ascii_file_format
use     location_mod, only : location_type, set_location, get_location, &
                             is_vertical, VERTISSURFACE, VERTISPRESSURE
use  assim_model_mod, only : interpolate

use     obs_kind_mod, only : QTY_SURFACE_PRESSURE, &
                             QTY_TEMPERATURE, &
                             QTY_SKIN_TEMPERATURE, &
                             QTY_NADIR_RADIANCE

use ensemble_manager_mod, only : ensemble_type

use obs_def_utilities_mod, only : track_status

implicit none
private

public :: set_TES_nadir, get_TES_nadir, write_TES_nadir_obs, read_TES_nadir_obs, &
          interactive_TES_nadir_obs, get_expected_TES_nadir_obs

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.


! All variables declared here are in MODULE STORAGE
character(len=129) :: errstring

! parameters for dealing with the Planck function
real(r8), parameter :: c_light = 2.99792458e10_r8     ! cm s^-1
real(r8), parameter :: h_planck = 6.626069e-34_r8     ! J s
real(r8), parameter :: k_boltz = 1.3806662e-23_r8     ! J K^-1
real(r8) :: alfa, beta

! other parameters required for code
real(r8), parameter :: co2frac = 0.9532_r8            ! percentage
real(r8), parameter :: gascon = 1.9114e6_r8           ! cm^2 K^-1 s^-2
real(r8), parameter :: mb_to_cgs = 1.0e3_r8           ! g cm^-1 s^-1
real(r8), parameter :: stpconv = 3709.5_r8            ! g cm^-1 s^-1 K^-1
real(r8), parameter :: tes_g = 374.0_r8               ! g cm^-2
real(r8), parameter :: amu_to_kgram = 1.672e-27_r8    ! kg amu^-1
real(r8), parameter :: cm_to_m = 0.01_r8              ! m cm^-1
integer,  parameter :: molwt = 44                     ! amu (mean Mol Wgt of CO2)
real(r8) :: fac_coeff1

! module-wide bits read in by init_corrk
integer  :: nv, np, nt, ng
real(r8), allocatable :: ck_v(:), ck_p(:), ck_t(:), gw(:)
real(r8), allocatable :: values(:,:,:,:)

! Create a private module derived type to store extra observation metadata in
!   This is modeled after Nancy's construction for obs_def_gps_mod
type TES_nadir_type
   private
   integer         :: scan_length
   real(r8)        :: wavenumber
   real(r8)        :: emission_angle
   real(r8)        :: l_sub_s
end type TES_nadir_type
integer  :: keycount

! The idea here is that TES_nadir_data will be allocated to the namelist 
!   specified value of max_TES_nadir_obs
type(TES_nadir_type), allocatable :: TES_data(:)

!----------------------------------------------------------------------------
! Namelist with default values
!
!  1.          N_layers :: number of vertical layers used in code
!  2.           n_gauss :: number of Gaussian quadrature abscissas
!  3.         corrk_dir :: path to directory containing the k-table files
!  4.          fn_corrk :: name of file within corrk_dir containing overview data
!  5.              ptop :: RT assumed pressure at top of domain (mbar)
!  6.       isothermalP :: RT assumed level above which atmosphere is isothermal
!  7.        fixed_dust :: logical for whether dust is parameterized (MCD-MGS)
!  8. max_TES_nadir_obs :: maximum number of TES nadir observations to process
!
! Default namelist values
integer              :: N_layers              = 30
integer              :: n_gauss               = 20
character(len = 129) :: corrk_dir             = '/home/wglawson/lee/TES_ktables/'
character(len = 129) :: fn_corrk              = 'lattice_coordinates.dat'
real(r8)             :: ptop                  = 1.0e-4_r8
real(r8)             :: isothermalP           = 4.0e-4_r8
logical              :: fixed_dust            = .true.
integer              :: max_TES_nadir_obs     = 100000

namelist /obs_def_TES_nadir_mod_nml/ N_layers, n_gauss, corrk_dir, fn_corrk, &
                                     ptop, isothermalP, fixed_dust, &
                                     max_TES_nadir_obs

contains

!----------------------------------------------------------------------------

  subroutine initialize_module
!----------------------------------------------------------------------------
! subroutine initialize_module
!
! Here we want to do the one-time overhead associated with using this observation
!   type, including:  registering module,  reading and writing out namelist,  
!   define used constants,  allocate storage for arrays,  read in necessary 
!   external data (possibly via other defined subroutines & functions)

! local variables; most "useful" variables here should be declared above in 
!   module storage
integer  :: iunit, io, k, status, rc


! a DART tradition
call register_module(source, revision, revdate)

! Global count of all TES nadir observations from any input file
keycount = 0

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_def_TES_nadir_mod_nml", iunit)
read(iunit, nml = obs_def_TES_nadir_mod_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_def_TES_nadir_mod_nml")

if (do_nml_file()) write(nmlfileunit, nml=obs_def_TES_nadir_mod_nml)
if (do_nml_term()) write(    *      , nml=obs_def_TES_nadir_mod_nml)

! Allocate arrays whose dimensions depend on namelist values
allocate( gw(n_gauss) )
allocate( TES_data(max_TES_nadir_obs), stat = rc )
if (rc /= 0) then
   write(errstring, *) 'Initial allocation failed for TES nadir observation data,', &
                       'itemcount = ', max_TES_nadir_obs
   call error_handler(E_ERR,'initialize_module', errstring, &
                      source, revision, revdate)
endif

! Calculate and store some derived "conglomerate" parameters:
!   alfa = 2.h.c^2 with units of J cm^2 s^-1 -- for Planck Law
alfa = 2.0_r8 * h_planck * c_light * c_light
!   beta = h.c/k with units of K cm -- for Planck Law
beta = h_planck * c_light / k_boltz
!   fac_coeff1 = co2frac * gascon * mb_to_cgs / (tes_g * stpconv) -- for 
!     RT of optical depth once divided by the cosine of emission angle
fac_coeff1 = co2frac * gascon * mb_to_cgs / ( tes_g * stpconv )

! Get the Gaussian weights for the Correlated K approach (size n_gauss)
call get_gaussian_weights( n_gauss,gw )

! Read in the K-tables -- should I read in both the scan_length = 1 tables
!   AND the scan_length = 2 tables (i.e., lo-res and hi-res?)
! init_corrk will allocate arrays ck_v, ck_p, ck_t, and values
call init_corrk

module_initialized = .true.

end subroutine initialize_module


  subroutine get_gaussian_weights( n,w )
!----------------------------------------------------------------------------
! subroutine get_gaussian_weights( n,w )
! 
! This my Fortran 90 version of "gauleg_LN_CO2.f90", which is in the original
!   TES forward operator code bundle I got from NASA GSFC (M. Smith et al.).
!   gauleg_LN_CO2.f90's lineage appears to have its roots from Numerical Recipes
!   and was originally coded by Mike Smith.  Here I have stripped out some
!   apparently unneeded and unused stuff.  I have tested this directly against
!   gauleg_LN_CO2.f90, and it performs identically.  

integer,  intent( in) :: n
real(r8), intent(out) :: w(n)

! local variables
real(r8), parameter :: small = 5.0e-10_r8
integer  :: m, j, k
real(r8) :: xm, xl, z, z1, p1, p2, p3, pp, rn

! initialize weights and set other initial values
w = 0.0_r8
m = int( (n + 1)/2 )
xm = 0.5_r8
xl = 0.5_r8
rn = real( n,r8 )

do j = 1, m
   z = cos( PI * ( real( j,r8 ) - 0.25_r8 ) / ( rn + 0.5_r8 ) ) 
   z1 = 1.0e6_r8
   small_test: do
      p1 = 1.0_r8
      p2 = 0.0_r8
      do k = 1, n
         p3 = p2
         p2 = p1
         p1 = ( real( 2*k-1,r8 ) * z * p2 - real( k-1,r8 ) * p3 ) / real( k,r8 )
      end do
      pp = rn * ( z * p1 - p2 ) / ( z * z - 1.0_r8 )
      z1 = z
      z = z1 - p1 / pp
      if ( abs( z - z1 ) < small ) exit small_test
   end do small_test
   w(j) = 2.0_r8 * xl / ( ( 1.0_r8 - z * z ) * pp * pp )
   w(n+1-j) = w(j)
end do

end subroutine get_gaussian_weights


  subroutine init_corrk
!----------------------------------------------------------------------------
! subroutine init_corrk
! 
! This my Fortran 90 version of "init_corrk_LN_CO2.f90", which is in the original
!   TES forward operator code bundle I got from NASA GSFC (M. Smith et al.).
!   This code reads in the Temperatures, Pressures, and Wavenumbers valid for
!   the k-tables stored in a specified directory (namelist items : corrk_dir).
!   It then goes from k-table file to file and reads in the full k-tables.
! I have tested this subroutine in a stand-alone mode and it works.

integer  :: k_unit, strlen1, strlen2, strlen, istat
integer  :: iflag, k, tmp, iv, ig, it, ip
character(len = 129) :: fn, msgstr
character(len = 5)   :: x

! In the future we will likely have to point this subroutine to both lo-res and
!   hi-res k-tables
! Use DART utilities for opening file and getting a unit number
! Read executive lattice file -- in namelist as fn_corrk, stored in corrk_dir

! 1. Construct file name for lattice file
strlen1 = len_trim( corrk_dir )
strlen2 = len_trim( fn_corrk )
strlen = strlen1 + strlen2
fn(:strlen1) = trim( corrk_dir )
fn((strlen1+1):strlen) = trim( fn_corrk )

! 2. Get a free I/O file unit from DART
k_unit = get_unit()

! 3. Open the file and confirm success
open( UNIT=k_unit, FILE=trim(fn), STATUS='old', IOSTAT=istat )
if ( istat > 0 ) then
   write(msgstr,*)'istat > 0 for opening fn_corrk :: istat = ',istat
   call error_handler(E_ERR, 'obs_def_TES_nadir_mod:init_corrk', msgstr, &
                      source, revision, revdate)
end if

! 4. Begin reading process

! 4. a. iflag -- ought to be equal to 1
read( UNIT=k_unit, FMT=*, IOSTAT=istat ) iflag
if ( istat > 0 ) then
   write(msgstr,*)'istat > 0 for reading iflag :: istat = ',istat
   call error_handler(E_ERR, 'obs_def_TES_nadir_mod:init_corrk', msgstr, &
                      source, revision, revdate)
end if
if ( iflag /= 1 ) then
   write(msgstr,*)'iflag .ne. 1 in fn_corrk :: iflag = ',iflag
   call error_handler(E_ERR, 'obs_def_TES_nadir_mod:init_corrk', msgstr, &
                      source, revision, revdate)
end if

! 4. b. nv -- number of wavenumbers
read( UNIT=k_unit, FMT=*, IOSTAT=istat ) nv
if ( istat > 0 ) then
   write(msgstr,*)'istat > 0 for reading nv :: istat = ',istat
   call error_handler(E_ERR, 'obs_def_TES_nadir_mod:init_corrk', msgstr, &
                      source, revision, revdate)
end if

! 4. c. np -- number of pressures
read( UNIT=k_unit, FMT=*, IOSTAT=istat ) np
if ( istat > 0 ) then
   write(msgstr,*)'istat > 0 for reading np :: istat = ',istat
   call error_handler(E_ERR, 'obs_def_TES_nadir_mod:init_corrk', msgstr, &
                      source, revision, revdate)
end if

! 4. d. nt -- number of temperatures  
read( UNIT=k_unit, FMT=*, IOSTAT=istat ) nt
if ( istat > 0 ) then
   write(msgstr,*)'istat > 0 for reading nt :: istat = ',istat
   call error_handler(E_ERR, 'obs_def_TES_nadir_mod:init_corrk', msgstr, &
                      source, revision, revdate)
end if

! 4. e. ng -- number of gaussian abscissae
read( UNIT=k_unit, FMT=*, IOSTAT=istat ) ng
if ( istat > 0 ) then
   write(msgstr,*)'istat > 0 for reading ng :: istat = ',istat
   call error_handler(E_ERR, 'obs_def_TES_nadir_mod:init_corrk', msgstr, &
                      source, revision, revdate)
end if

! 5. Allocate array sizes based on what we just read in
allocate( ck_v( nv ) )
allocate( ck_p( np ) )
allocate( ck_t( nt ) )
allocate( values( ng,nt,np,nv ) )

! 6. Read more data into the arrays we just allocated

! 6. a. ck_v -- array of wavenumber files in corrk_dir
read( UNIT=k_unit, FMT=*, IOSTAT=istat ) ( ck_v(k), k=1,nv )
if ( istat > 0 ) then
   write(msgstr,*)'istat > 0 for reading ck_v :: istat = ',istat
   call error_handler(E_ERR, 'obs_def_TES_nadir_mod:init_corrk', msgstr, &
                      source, revision, revdate)
end if

! 6. b. ck_p -- array of pressures included in k-tables
read( UNIT=k_unit, FMT=*, IOSTAT=istat ) ( ck_p(k), k=1,np )
if ( istat > 0 ) then
   write(msgstr,*)'istat > 0 for reading ck_p :: istat = ',istat
   call error_handler(E_ERR, 'obs_def_TES_nadir_mod:init_corrk', msgstr, &
                      source, revision, revdate)
end if

! 6. c. ck_t -- array of temperatures included in k-tables
read( UNIT=k_unit, FMT=*, IOSTAT=istat ) ( ck_t(k), k=1,nt )
if ( istat > 0 ) then
   write(msgstr,*)'istat > 0 for reading ck_t :: istat = ',istat
   call error_handler(E_ERR, 'obs_def_TES_nadir_mod:init_corrk', msgstr, &
                      source, revision, revdate)
end if

! 7. Close the lattice file
close( UNIT=k_unit )

! 8. Now loop through the k-tables files -- just reuse the same unit number
do iv = 1, nv
   ! Build the filename for the k-table file
   fn = ' '
   write(x,'(i5)') nint( ck_v(iv) * 100.0_r8 )
   fn = trim(  corrk_dir ) // 'ln_distribution_' // x // '.dat'
   ! Open the file & check for success
   open( UNIT=k_unit, FILE=trim(fn), STATUS='old', IOSTAT=istat )
   if ( istat > 0 ) then
      write(msgstr,*)'istat > 0 for opening k-table file ',iv,' :: istat = ',istat
      call error_handler(E_ERR, 'obs_def_TES_nadir_mod:init_corrk', msgstr, &
                      source, revision, revdate)
   end if
   ! Read in the k-table data into the array "values"
   read( UNIT=k_unit, FMT=*, IOSTAT=istat ) (((values(ig,it,ip,iv), ig=1,ng), it=1,nt), ip=1,np)
   if ( istat > 0 ) then
      write(msgstr,*)'istat > 0 for reading k-table file ',iv,' :: istat = ',istat
      call error_handler(E_ERR, 'obs_def_TES_nadir_mod:init_corrk', msgstr, &
                      source, revision, revdate)
   end if
   ! Close the file
   close( UNIT=k_unit )
end do

end subroutine init_corrk


  subroutine set_TES_nadir(teskey, scan_length, wavenumber, emission_angle, l_sub_s)
!------------------------------------------------------------------------------
!
! subroutine set_TES_nadir(teskey, scan_length, wavenumber, emission_angle, l_sub_s)
!
! Increment key and set all private data for this observation

integer,          intent(out) :: teskey
integer,          intent(in)  :: scan_length
real(r8),         intent(in)  :: wavenumber, emission_angle, l_sub_s

if ( .not. module_initialized ) call initialize_module

keycount = keycount + 1
teskey = keycount

if(teskey > max_TES_nadir_obs) then
   write(errstring, *) 'key (',teskey,') exceeds max_TES_nadir_obs (',max_TES_nadir_obs,')'
   call error_handler(E_ERR,'set_TES_nadir', errstring, &
                      source, revision, revdate)
endif

TES_data(teskey)%scan_length     = scan_length
TES_data(teskey)%wavenumber      = wavenumber
TES_data(teskey)%emission_angle  = emission_angle
TES_data(teskey)%l_sub_s         = l_sub_s

end subroutine set_TES_nadir


  subroutine get_TES_nadir(teskey, scan_length, wavenumber, emission_angle, l_sub_s)
!------------------------------------------------------------------------------
!
! subroutine get_TES_nadir(teskey, scan_length, wavenumber, emission_angle, l_sub_s)
!
! Increment key and set all private data for this observation

integer,          intent(in)  :: teskey
integer,          intent(out) :: scan_length
real(r8),         intent(out) :: wavenumber, emission_angle, l_sub_s

if ( .not. module_initialized ) call initialize_module

if (teskey < 1 .or. teskey > keycount) then
   write(errstring, *) 'key (',teskey,') out of valid range (1<=key<=',keycount,')'
   call error_handler(E_ERR,'get_TES_nadir', errstring, &
                      source, revision, revdate)
endif

scan_length    = TES_data(teskey)%scan_length
wavenumber     = TES_data(teskey)%wavenumber
emission_angle = TES_data(teskey)%emission_angle
l_sub_s        = TES_data(teskey)%l_sub_s

end subroutine get_TES_nadir


  subroutine write_TES_nadir_obs(teskey, ifile, fform)
!----------------------------------------------------------------------------
! subroutine write_TES_nadir_obs(teskey, ifile, fform)
!
! The following is largely lifted from obs_def_gps_mod.f90 and obs_def_radar_mod.f90

integer,          intent(in)           :: teskey, ifile
character(len=*), intent(in), optional :: fform

integer :: i

if ( .not. module_initialized ) call initialize_module

! write out obs_seq info 
if (ascii_file_format(fform)) then
   write(ifile,98) teskey
   write(ifile, *) TES_data(teskey)%scan_length, &
                   TES_data(teskey)%wavenumber, &
                   TES_data(teskey)%emission_angle, &
                   TES_data(teskey)%l_sub_s
else
   write(ifile) teskey
   write(ifile) TES_data(teskey)%scan_length, &
                TES_data(teskey)%wavenumber, &
                TES_data(teskey)%emission_angle, &
                TES_data(teskey)%l_sub_s
endif

98 format('TES_nadir_obs', i8)

end subroutine write_TES_nadir_obs


  subroutine read_TES_nadir_obs(teskey, ifile, fform)
!----------------------------------------------------------------------------
! subroutine read_TES_nadir_obs(teskey, ifile, fform)
!
! The following is largely lifted from obs_def_gps_mod.f90 and obs_def_radar_mod.f90

integer,          intent(out)          :: teskey
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

integer             :: keyin    ! the metadata key in the current obs sequence

integer             :: scan_length
real(r8)            :: wavenumber, emission_angle, l_sub_s
character(len=13)   :: header

if ( .not. module_initialized ) call initialize_module

if (ascii_file_format(fform)) then
   read(ifile, fmt='(a13, i8)') header, keyin    ! throw away keyin
   if(header /= 'TES_nadir_obs') then
      write(errstring,*)'Expected header "TES_nadir_obs" in input file'
      call error_handler(E_ERR,'read_TES_nadir_obs',errstring, &
                               source, revision, revdate)
   end if
   read(ifile,*) scan_length, wavenumber, emission_angle, l_sub_s
else
   read(ifile) keyin          ! read and throw away
   read(ifile) scan_length, wavenumber, emission_angle, l_sub_s
endif

! increment key and set all private data for this observation
call set_TES_nadir(teskey, scan_length, wavenumber, emission_angle, l_sub_s)

end subroutine read_TES_nadir_obs


  subroutine interactive_TES_nadir_obs(teskey)
!----------------------------------------------------------------------------
! subroutine interactive_TES_nadir_obs(teskey)
!
! The following is largely lifted from obs_def_1d_state_mod.f90 and obs_def_radar_mod.f90

integer, intent(out) :: teskey
integer              :: scan_length
real(r8)             :: wavenumber, emission_angle, l_sub_s

if ( .not. module_initialized ) call initialize_module

if( keycount >= max_TES_nadir_obs ) then
   write(errstring, *)'Not enough room left for another TES_nadir_obs.'
   call error_handler(E_MSG,'interactive_TES_nadir_obs',errstring, source, revision, revdate)
   write(errstring, *)'Can only have max_TES_nadir_obs (currently ',max_TES_nadir_obs,')' 
   call error_handler(E_ERR,'interactive_TES_nadir_obs',errstring, source, revision, revdate)
endif

! get obs metadata
write(*, *) 'Creating an interactive TES nadir observation'
write(*, *)
write(*, *) 'Input the TES scan length:'
write(*, *) '   Enter 1 for low-res (nominally 10/cm)'
write(*, *) '   Enter 2 for high-res (nominally 5/cm)'
read(*,*) scan_length

write(*, *)
write(*, *) 'Input the TES channel wavenumber in cm^-1'
read(*,*) wavenumber

write(*, *)
write(*, *) 'Input the TES emission angle in degrees'
read(*,*) emission_angle

write(*, *)
write(*, *) 'Input the L_sub_s in degrees for this observation'
read(*,*) l_sub_s

! increment key and set all private data for this observation
call set_TES_nadir(teskey, scan_length, wavenumber, emission_angle, l_sub_s)

write(*, *)
write(*, *) '   End of specialized section for TES nadir observations.'
write(*, *) '   You will now have to enter the regular obs information.'
write(*, *)

end subroutine interactive_TES_nadir_obs


  subroutine get_expected_TES_nadir_obs(state_handle, ens_size, location, teskey, expected_obs, istatus)
!----------------------------------------------------------------------------
! subroutine get_expected_TES_nadir_obs(state, location, teskey, val, istatus)
!
! We are going to try to capture the code here from the original
!   TES forward operator code bundle I got from NASA GSFC (M. Smith et al.).
!   The code here is split over several routines and functions including:
!     guess_LN_CO2.f90
!     trans_LN_CO2.f90
!     inv_LN_CO2.f90
!     atmosrad_LN_CO2.f90
!     planck.f90
!
! In order to calculate the expected TES radiance, we will need:
!   state vector
!   latitude on planet
!   longitude on planet
!   local time
!   UTC -- determined by longitude + local time
!   sol/date -- determined by L_sub_s (?)
!   L_sub_s
!   emission angle
!   TES channel -- specified as wavenumber 
!   scan length (1 = lo-res, 2 = hi-res)
!   uncertainty

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer, intent(in)             :: teskey
real(r8), intent(out)           :: expected_obs(ens_size)
integer, intent(out)            :: istatus(ens_size)


! The first item of business is to get the vertical columns of temperature 
!   and pressure, as well as the surface values, all interpolated to the 
!   correct location for the nadir observation.  This is not quite as simple
!   as it may sound because the nominal time and location may not synch 
!   correctly with the cited local time value.  I think the issue is the 
!   martian equation of time -- the GCM approximates time as mean solar time,
!   whereas the satellite is viewing w.r.t. true local solar time.  Worse,
!   the GCM isn't exactly on mean solar time because it assumes 669 sols per
!   year instead of 668.599, or whatever the appropriate value should be.
! What is most important, I imagine, is that the location and the local time
!   be as close to matching as possible -- since local time is a proxy for
!   longitude, the location really specifies the local time FOR A GIVEN 
!   value of UTC.  To match location and local time, we may need to shift
!   the nominal time tag of the observation either forward or backward.  
!   This is the essence of the equation of time, but I don't think even that
!   is accurate enough for our faked GCM clocks and calendars **check this**
! Okay, TES doesn't include a Date or UTC field, just L_s and Local_time 
!   (and also Mars J2000 positions and Ephemeris_time) -- we will need to
!   use L_s and Local_time in conjunction with longitude to back-calculate
!   UTC fields and Date fields so that time_manager will know at what time
!   step each observation is valid for -- this is something to be handled 
!   during the construction of the obs_seq files, NOT here.

! istatus = 0  :: nominally okay
!         = 1  :: did not find k-table for specified wavenumber
!         = 2  :: trouble with interpolating to find PSFC
!         = 3  :: trouble with interpolating to find TSK
!         = 4  :: could not find acceptable value for isothermalP
!         = 5  :: trouble with interpolating to find T below isothermalP
!         = 11 :: corrk trouble -- pressure out of range
!         = 12 :: corrk trouble -- temperature too low
!         = 13 :: corrk trouble -- temperature too high & outside of 15 micron center
!         = 14 :: corrk trouble -- should never return 4
!         = 15 :: corrk trouble -- corner entries in array values are = -1

integer  :: k, j, istat, iv
real(r8) :: obsloc(3), lon, lat
real(r8), dimension(ens_size) :: psfc, tsk, delz, isoP
real(r8) :: vtmp, ptmp(ens_size), ttmp(ens_size), isothermalT(ens_size), coeff1
real(r8) :: atmos_now, surface_now, rad_now(ens_size)
real(r8), parameter :: tolerance = 0.05_r8
logical,  parameter :: debug = .false., debug_ng = .false.

type(location_type) :: location2
integer  :: which_vert

!> @todo See how much how stack this is using
real(r8) :: ret_p(0:N_layers, ens_size), ret_t(0:N_layers, ens_size)
real(r8) :: t_mid(N_layers, ens_size)
real(r8) :: corrk(0:N_layers, n_gauss, ens_size)
real(r8) :: corrk_tmp(n_gauss, ens_size)
real(r8) :: dust_opt_dep(0:N_layers, ens_size)
real(r8) :: opt_dep(0:N_layers,n_gauss, ens_size)
real(r8) :: trans10k(0:N_layers,n_gauss, ens_size)
real(r8) :: trans(0:N_layers, ens_size)

real(r8) :: isothermalT_dummy(ens_size)
integer  :: original_tmp_istat(ens_size)
real(r8) :: ret_t_dummy(ens_size)

integer  :: imem ! loop index for ensemble members
! Istatus values for the various calls to interpolation
integer, dimension(ens_size) :: psfc_istat, ret_t_istat, tmp_istat, tsk_istat
! Return codes for k-table look-up
integer, dimension(ens_size) :: rc1_istat, rc2_istat, rc3_istat, rc4_istat, rc5_istat, rc
logical  :: return_now ! whether or not you can return after interpolation (e.g. all ens 
                       ! members failed)

! Make sure the module is initialized -- we need the k-tables in memory!
if ( .not. module_initialized ) call initialize_module

istatus(:) = 0

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Find the wavenumber in k-table!

! It is important to know not only the wavenumber (TES channel) of the 
!   observation in question but also the index of that wavenumber in the
!   array of possible k-table wavenumbers
vtmp = TES_data(teskey)%wavenumber
! Find index of this wavenumber in array of k-table wavenumber values
iv = 0
do k = 1, nv
   if ( abs( ck_v(k) - vtmp ) <= tolerance ) then
      iv = k
      exit
   end if
end do
if (debug) print*, 'iv = ', iv, ';  vtmp = ', vtmp
! If we haven't found a wavenumber in ck_v that is within tolerance of this
!   observation's cited wavenumber value, then we have a problem -- set
!   val to missing and the return code to a non-zero value
if ( iv == 0 ) then
   write(errstring, *)'Did not find wavenumber for ob!; TES obs key ', teskey
   call error_handler(E_WARN,'get_expected_TES_nadir_obs', errstring, &
                      source, revision, revdate)
   expected_obs(:) = missing_r8
   istatus(:) = 1
   return
end if

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Vertical profile interpolation

! Since nadir radiances are related to vertical integrals, a vertical
!   location does not really make sense -- for threed_sphere/location, we 
!   will use and expect a vertical location of "undefined".  Don't let this
!   necessarily break the code though.
if ( .not. is_vertical(location, "UNDEFINED") ) then
   write(errstring, *)'vertical location should be undefined; TES obs key ', teskey
   call error_handler(E_WARN,'get_expected_TES_nadir_obs', errstring, &
                      source, revision, revdate)
endif

! Convert from location_type to real numbers
obsloc = get_location( location )

! Strip out lat and lon for interpolation purposes
lon = obsloc( 1 )    ! EAST longitude in degrees: 0 to 360
lat = obsloc( 2 )    ! latitude in degrees: -90 to 90
if (debug) print*, 'lon = ', lon, ';  lat = ', lat

! Get surface pressure at location
which_vert = VERTISSURFACE
location2 = set_location( lon, lat, missing_r8, which_vert )
call interpolate(state_handle, ens_size, location2, QTY_SURFACE_PRESSURE, psfc, psfc_istat )
call track_status(ens_size, psfc_istat, expected_obs, istatus, return_now)
if ( any(psfc_istat > 0) ) then
   !> @todo Are warnings written by tasks other than task 0?
   write(errstring, *)'trouble with PSFC interpolation, key = ', teskey, '; istat = ', psfc_istat
   call error_handler(E_WARN,'get_expected_TES_nadir_obs', errstring, &
                      source, revision, revdate)
endif
where (psfc_istat /= 0) istatus = 2
if (return_now) return


! Make sure all pressures are actually in mbars, as the forward operator
!   code (i.e., k-tables) is expecting.  The model uses Pascals.
! ptop is already in mbars
where (istatus == 0) psfc = psfc * 1.0e-2_r8

! Now surface temperature at locations
call interpolate(state_handle, ens_size, location2, QTY_SKIN_TEMPERATURE, tsk, tsk_istat)
call track_status(ens_size, tsk_istat, expected_obs, istatus, return_now)

if ( any(tsk_istat > 0) ) then
   write(errstring, *)'trouble with TSK interpolation, key = ', teskey, '; istat = ', tsk_istat
   call error_handler(E_WARN,'get_expected_TES_nadir_obs', errstring, &
                      source, revision, revdate)
endif
where (tsk_istat /= 0 .and. psfc_istat == 0) istatus = 3
if(return_now) return


! We want to construct the vertical temperature and pressure profiles -- 
!   Philosophically, we do not want to know what the vertical resolution of 
!   our model is -- rather we want to specify actual vertical levels on
!   which to interpolate the values we are looking for.  If we directly
!   specify the pressure levels we want, then we will only have to interpolate
!   temperature.
! We will use the namelist item N_layers to accomplish this, though in a sigma
!   coordinate fashion where we have N_layers evenly spaced in log pressure
!   over the pressure difference between p_surf and p_top.
! p([0:N]) = psfc * exp( -( log(psfc/ptop)/N ) * [0:N] )
do imem=1, ens_size
   if (istatus(imem) == 0) then
      ret_p(1:N_layers, imem) = 0.0_r8
      ret_p(0, imem) = psfc(imem)
      delz(imem) = log( psfc(imem) / ptop ) / real(N_layers,r8)
   else
      ret_p(1:N_layers, imem) = missing_r8
      delz(imem) = missing_r8
   endif
enddo

where (istatus == 0)
      ret_p(k, :) = psfc * exp( -delz * real(k,r8) )
end where

if (debug) print*, 'ret_p = ', ret_p
if (debug) print*, ' '

! Use model's interpolate routine to find isothermalT, the temperature when the
!   model attains height isothermalP -- keep in mind that interpolate assumes 
!   pressure is in Pa, not mbars.
! Store namelist item isothermalP in isoP so that we can change it if necessary
isoP(:) = isothermalP  ! At this point, isoP is the same across the ensemble
which_vert = VERTISPRESSURE
ptmp = isoP * 1.0e2_r8

location2 = set_location( lon, lat, ptmp(1), which_vert ) ! At this point, same across the ensemble

call interpolate(state_handle, ens_size, location2, QTY_TEMPERATURE, isothermalT, tmp_istat)
call track_status(ens_size, tmp_istat, expected_obs, istatus, return_now)
if (debug) print*, 'isothermalT = ', isothermalT

! If we had trouble it is MOST LIKELY because isothermalP is above the model's
!   value of ptop -- try lowering (in altitude) isothermalP.  For now we will
!   (arbitrarily) try lowering it 6 times, and if we haven't found an acceptable
!   value by then, then return with missing value because we won't be able to
!   construct a reasonable T(p) profile for vertical intergration.
!> @todo This needs refactoring at some point
if ( any(tmp_istat > 0) ) then
   original_tmp_istat = tmp_istat ! original pass/fails
   LOWERING : do k = 1, 6

      where (tmp_istat /= 0)  ! keep going if still not had an interpolation pass
         isoP = 2.0_r8 * isoP
         ptmp = 2.0_r8 * ptmp
      end where
      location2 = set_location( lon, lat, ptmp(1), which_vert )

      call interpolate( state_handle, ens_size, location2, QTY_TEMPERATURE, isothermalT_dummy, tmp_istat )

      if (debug) print*, 'k = ', k, '  isothermalT = ', isothermalT
      write(errstring, *)'lowered isothermalP, key = ',teskey,'; new isoP = ', isoP
      call error_handler(E_MSG,'get_expected_TES_nadir_obs', errstring, &
                            source, revision, revdate)

      where (original_tmp_istat /= 0 .and. tmp_istat == 0)
         ! update the original because we now have an interpolation that passed.
         original_tmp_istat = 0
         isothermalT = isothermalT_dummy
      end where

      if ( all(original_tmp_istat == 0 )) exit LOWERING! we have interpolations that passed for all ensemble members

   end do LOWERING

   if ( k == 6 ) then
      write(errstring, *)'could not find acceptable isothermalP level for a least one ensemble member, key = ',teskey, &
                           '; last isoP = ', isoP
      call error_handler(E_WARN,'get_expected_TES_nadir_obs', errstring, &
                           source, revision, revdate)
   end if

   where(original_tmp_istat /= 0) ! did not find an acceptable isothermalP level
      istatus = 4
      expected_obs = missing_r8
   end where
   if (all(original_tmp_istat /=0)) return

end if

! Now set the corresponding temperature profile -- TSK is surface; interpolate
!   for others -- which_vert has already been set to VERTISPRESSURE
ret_t(:, :) = 0.0_r8
ret_t(0, :) = tsk

!>@todo FIXME Possibly we could overload interpolate and set_location to take
!> an array of locations
!> then have interpolate sort the list into unique values, then call model_interpolate
!> for only the unique locations. Or sort into unique values of isoP here.
!> For now this is a loop around ensemble_size calling interpolate - overcommunicating
do imem = 1, ens_size
   do k = 1, N_layers
      ! Test whether we are surface-ward of the namelist specified isothermalP,
      !   but compare to isoP in case we had to bring isothermalP closer to the
      !   surface to abide by the model's own value for p_top.
      if ( (ret_p(k, imem) > isoP(imem)) )  then
         ! For model interpolation, make sure we are back in Pa
         ptmp(imem) = ret_p(k, imem) * 1.0e2_r8

         location2 = set_location( lon, lat, ptmp(imem), which_vert )

         call interpolate( state_handle, ens_size, location2, QTY_TEMPERATURE, ret_t_dummy( k ), ret_t_istat )

         if ( ret_t_istat(imem) > 0 ) then
            write(errstring, *)'trouble with T interpolation, key = ', teskey, '; istat = ', &
                                 istat,'; level k = ', k
            call error_handler(E_WARN,'get_expected_TES_nadir_obs', errstring, &
                                 source, revision, revdate)
            istatus(imem) = 5
            expected_obs(imem) = missing_r8
         end if

         if (all(istatus /= 0)) return

      else
         ret_t(k, imem) = isothermalT(imem)
      end if
   end do

enddo


if (debug) print*, 'ret_t = ', ret_t
if (debug) print*, ' '

! We will also want T and P on the mid-points between the layers above.  The
!   original TES code just interpolates linearly in log(P), but we could use
!   our model_mod to interpolate more closely with what the model says.  The
!   TES approach is certainly faster -- I wonder how much of a difference it 
!   makes?  
! Here we will use the TES approach...
!   If we are not explicitly interpolating from the model state vector, then
!   don't actually need p_mid because t_mid is simply the mid-points of ret_t
!p_mid = 0.0_r8
t_mid = 0.0_r8
do k = 1, N_layers
   !p_mid(k) = exp( 0.5_r8 * ( log(ret_p(k)) + log(ret_p(k-1)) ) )
   where (istatus == 0)
      t_mid(k, :) = 0.5_r8 * ( ret_t(k, :) + ret_t(k-1, :) )
   end where
end do
if (debug) print*, 't_mid = ', t_mid
if (debug) print*, ' '

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Dust stuff

! We will need access to the vertical dust distribution for the (lon,lat)
!   location as well.  We can either make TAU_OD part of the state vector
!   (including TRC01 and TRC02 if "active dust") or we can try to give 
!   the forward operator access to the fixed dust distribution subroutine.
!   This second option is definitely against the philosophy of DART.

! For now we will stick with the logical "fixed_dust" construction and call
!   the function mcd_mgs below -- keep in mind that mcd_mgs assumes inputs 
!   are in Pa, NOT mbars!
if ( fixed_dust ) then
   do imem = 1, ens_size
      if (istatus(imem) /= 0) cycle
      do k = 0, N_layers
         dust_opt_dep(k, imem) = mcd_mgs( TES_data(teskey)%l_sub_s, lat, ret_p(k, imem) )
      end do
   enddo
else
   ! INSERT DUST CODE HERE FOR ACTIVE DUST RUNS
end if
if (debug) print*, 'VIS dust_opt_dep = ', dust_opt_dep
if (debug) print*, ' '

! dust_opt_dep above is the optical depth due to dust in the visible!  We now
!   need a way to convert the visible tau to IR tau relevant for whatever
!   wavenumber we are currently treating
! Push this task into a subroutine below, and go and get the multiplier data
!   from the Forget 1998 dust paper later.
do imem = 1, ens_size
   if(istatus(imem) /= 0) cycle
   call dust_vis_to_ir( dust_opt_dep(:, imem), N_layers+1, vtmp )
enddo
if (debug) print*, ' IR dust_opt_dep = ', dust_opt_dep
if (debug) print*, ' '

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Do k-table look-up / interpolation

! Now grab the appropriate entries in the correlated-k tables for the 
!   pressures and temperatures of the vertical column we've created
corrk = 0.0_r8


! The original (non-ensemble code) got the first failure, then returned
! Tracking the status here.
LAYERS: do k = 0, N_layers
   ptmp = ret_p(k, :)
   ttmp = ret_t(k, :)
   do imem = 1, ens_size
     if (istatus(imem) /= 0 ) cycle  ! to avoid missing_r8 going into this function
     call get_corrk(iv, vtmp, ptmp(imem), ttmp(imem), n_gauss, corrk_tmp(:, imem), rc(imem))
   enddo
   ! check return code:
   ! 1 = pressure out of range
   ! 2 = temperature too low
   ! 3 = temperature too high & outside of 15 micron center
   ! 4 = should never return 4
   ! 5 = corner entries in array values are = -1

! Return code 1
   if ( any(rc == 1) ) then  ! does every task print warnings?
      write(errstring, *)'Pressure is outside of k-table bounds; key = ',teskey, &
                           '; p(mbar) = ',ptmp,'; level = ',k
      call error_handler(E_WARN,'get_expected_TES_nadir_obs', errstring, &
                           source, revision, revdate)
   endif
   where (rc == 1 .and. istatus == 0) ! don't want to overwrite earlier failures
      expected_obs = missing_r8
      rc1_istat  = 10+rc
   endwhere
   call track_status(ens_size, rc1_istat, expected_obs, istatus, return_now)

! Return code 2
   if ( any(rc == 2)) then
      write(errstring, *)'Temperature is lower than k-table bounds; key = ',teskey, &
                           '; T(K) = ',ttmp,'; level = ',k
      call error_handler(E_WARN,'get_expected_TES_nadir_obs', errstring, &
                           source, revision, revdate)
   endif
   where (rc == 2 .and. istatus == 0)
      expected_obs = missing_r8
      rc2_istat = 10+rc
   endwhere
   call track_status(ens_size, rc2_istat, expected_obs, istatus, return_now)

! Return code 3
   if (any(rc == 3)) then
      write(errstring, *)'Temperature is higher than k-table bounds; key = ',teskey, &
                           '; T(K) = ',ttmp,'; level = ',k
      call error_handler(E_WARN,'get_expected_TES_nadir_obs', errstring, &
                              source, revision, revdate)
   endif
   where (rc == 3 .and. istatus == 0)
      expected_obs = missing_r8
      rc3_istat = 10+rc
   endwhere
   call track_status(ens_size, rc3_istat, expected_obs, istatus, return_now)

! Return code 4
   if (any(rc == 4)) then
      write(errstring, *)'This should not return; key = ',teskey
      call error_handler(E_WARN,'get_expected_TES_nadir_obs', errstring, &
                              source, revision, revdate)
   endif
   where (rc == 4 .and. istatus == 0)
      expected_obs = missing_r8
      rc4_istat = 10+rc
   endwhere
   call track_status(ens_size, rc4_istat, expected_obs, istatus, return_now)

   if (any(rc == 5)) then
      write(errstring, *)'Trying to interpolate where k-tables have no data; key = ',teskey
      call error_handler(E_WARN,'get_expected_TES_nadir_obs', errstring, &
                           source, revision, revdate)
   endif
   where (rc == 5 .and. istatus == 0)
      expected_obs = missing_r8
      rc5_istat = 10+rc
   endwhere
   call track_status(ens_size, rc4_istat, expected_obs, istatus, return_now)

   do imem = 1, ens_size
      if (istatus(imem) == 0) then
         corrk(k, :, imem) = corrk_tmp(:,imem)
      else
         corrk(k, :, imem) = missing_r8
      endif
   enddo

enddo LAYERS

if (debug .and. debug_ng) print*, 'corrk = ', corrk
if (debug .and. debug_ng) print*, ' '

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Start the radiative transfer part

ENSEMBLE_MEMBERS: do imem = 1, ens_size

   if (istatus(imem) /= 0) cycle ENSEMBLE_MEMBERS

   !> @todo This should go into a subroutine
   ! Take observation's emission angle into consideration here -- not very
   !   sensitive because all obs are nadir viewing, which means emission
   !   angle is within 2 degrees of nadir.  But still, it's a good idea.
   coeff1 = fac_coeff1 * COS( DEG2RAD * TES_data(teskey)%emission_angle )

   ! Evaluate opt_dep and trans10k
   opt_dep = 0.0_r8
   trans10k = 0.0_r8
   ! Top layer first
   opt_dep(N_layers,:, imem) = coeff1 * corrk(N_layers,:, imem) * ret_p(N_layers, imem) + dust_opt_dep(N_layers, imem)
   trans10k(N_layers,:, imem) = exp( -opt_dep(N_layers,:, imem) )
   ! Then fill in layers below
   do k = N_layers, 1, -1
      opt_dep(k-1,:, imem) = opt_dep(k,:, imem) + coeff1 * 0.5_r8 * ( corrk(k-1,:, imem) + corrk(k,:, imem) ) * &
                        ( ret_p(k-1, imem) - ret_p(k, imem) ) + dust_opt_dep(k-1, imem) - dust_opt_dep(k, imem)
      trans10k(k-1,:, imem) = exp( -opt_dep(k-1,:, imem) )
   end do
   if (debug .and. debug_ng) print*, 'opt_dep = ', opt_dep(:, :, imem)
   if (debug .and. debug_ng) print*, ' '
   if (debug .and. debug_ng) print*, 'trans10k = ', trans10k(:, :, imem)
   if (debug .and. debug_ng) print*, ' '

   ! Calculate transmittance at each level
   trans = 0.0_r8
   do k = 0, N_layers
      trans(k, imem) = sum( trans10k(k,:, imem) * gw )
   end do
   if (debug) print*, 'trans = ', trans(:, imem)
   if (debug) print*, ' '

   ! Get atmospheric contribution -- depends on t_mid
   atmos_now = 0.0_r8
   do k = 1, N_layers
      atmos_now = atmos_now + ( alfa * vtmp**3 )/( exp( beta * vtmp / t_mid(k, imem) ) - &
                                 1.0_r8 ) * ( trans(k, imem) - trans(k-1, imem) )
   end do
   atmos_now = atmos_now + ( alfa * vtmp**3 )/( exp( beta * vtmp / isothermalT(imem) ) - &
                                 1.0_r8 ) * ( 1 - trans(N_layers, imem) )
   if (debug) print*, 'atmos_now = ', atmos_now
   if (debug) print*, ' '

   ! Get surface contribution
   surface_now = 0.0_r8
   surface_now = ( alfa * vtmp**3 )/( exp( beta * vtmp / ret_t(0, imem) ) - 1.0_r8 ) * trans(0, imem)
   if (debug) print*, 'surface_now = ', surface_now
   if (debug) print*, ' '

   ! Total radiance measured at spacecraft
   rad_now(imem) = atmos_now + surface_now
   if (debug) print*, 'rad_now = ', rad_now
   if (debug) print*, ' '

enddo ENSEMBLE_MEMBERS

where (istatus == 0)
   expected_obs = rad_now
elsewhere
   expected_obs = missing_r8
endwhere

end subroutine get_expected_TES_nadir_obs


  function mcd_mgs(l_sub_s, lat, pressure) result(optical_depth)
!------------------------------------------------------------------------
! function mcd_mgs(l_sub_s, lat, pressure) result(optical_depth)
!
! This is a function for calculating the fixed dust amount according to
!   a specified function that supposedly matches the MGS era well and is
!   used in the Mars Climate Database (MCD) -- this is taken from Mars WRF.
! This has been tested successfully against the Matlab version of this code.

real(r8)              :: optical_depth
real(r8), intent( in) :: l_sub_s, lat, pressure

real(r8) :: glat, gls, zls, taueq, tauS, tauN, tauref, topdust, zp

glat = DEG2RAD * lat
gls = DEG2RAD * l_sub_s
zls = sin( gls - 2.76_r8 )
taueq = 0.2_r8 + (0.5_r8 - 0.2_r8) * (cos( 0.5_r8 * (gls - 4.363_r8) ))**14
tauS = 0.1_r8 + (0.5_r8 - 0.1_r8) * (cos( 0.5_r8 * (gls - 4.363_r8) ))**14
tauN = 0.1_r8

if ( lat >= 0.0_r8 ) then
   tauref = tauN + 0.5_r8 * (taueq - tauN) * (1.0_r8 + tanh( 0.1_r8 * (45.0_r8 - lat) ))
else
   tauref = tauS + 0.5_r8 * (taueq - tauS) * (1.0_r8 + tanh( 0.1_r8 * (45.0_r8 + lat) ))
end if

topdust = 60.0_r8 + 18.0_r8 * zls - (32.0_r8 + 18.0_r8 * zls) * (sin( glat ))**4 - &
              8.0_r8 * zls * (sin(glat))**5 
zp = ( 700.0_r8 / pressure )**( 70.0_r8 / topdust )

optical_depth = pressure * (tauref / 700.0_r8) * &
                   max( 1.0e-3_r8, exp( 7.0e-3_r8 * (1.0_r8 - max( zp,1.0_r8 )) ) )

end function mcd_mgs


  subroutine dust_vis_to_ir( tau, ltau, nu )
!------------------------------------------------------------------------
! subroutine dust_vis_to_ir( tau, ltau, nu )
!
! This subroutine takes in a vector of optical depths in the visible due to 
!   dust and outputs that same optical depth in the specified band in the IR
!   The conversion is through the top panel of figure 3 in Forget's 1998 GRL
!   dust paper.  The basic approach is to compare the relative magnitudes of
!   Qext throughout the 15 um band to Qext at 9 um, then to use Forget's 
!   canonical conversion of tau_vis / tau_9um = 2.

integer,  intent(  in ) :: ltau
real(r8), intent(inout) :: tau(ltau)
real(r8), intent(  in ) :: nu

integer,  parameter :: nv_lores = 24
real(r8), dimension(nv_lores) :: v_lores, vis_to_ir, fac

real(r8), parameter :: tolerance = 0.05_r8
integer  :: k, iv

! Set up the v_lores array -- this is just data entry from the 
!   lattice_coordinates.dat file:
v_lores(1:5)   = (/ 572.46_r8,  583.05_r8,  593.62_r8,  604.20_r8,  614.78_r8 /)
v_lores(6:10)  = (/ 625.37_r8,  635.94_r8,  646.52_r8,  657.11_r8,  667.68_r8 /)
v_lores(11:15) = (/ 678.26_r8,  688.85_r8,  699.43_r8,  710.00_r8,  720.59_r8 /)
v_lores(16:20) = (/ 731.17_r8,  741.75_r8,  752.32_r8,  762.91_r8,  773.49_r8 /)
v_lores(21:24) = (/ 784.06_r8,  794.65_r8,  805.23_r8,  815.80_r8 /)

! Find the index of nu
do k = 1, nv_lores
   if ( abs( v_lores(k) - nu ) <= tolerance ) then
      iv = k
      exit
   end if
end do

! Set up the vis_to_ir array -- this is just data entry from the estimates
!   from Forget 1998
! fac stores the conversion factors
fac(1:5)   = (/ 5.0_r8, 5.0_r8, 5.0_r8, 5.0_r8, 5.0_r8 /)
fac(6:10)  = (/ 6.2_r8, 7.6_r8, 7.1_r8, 6.2_r8, 6.6_r8 /)
fac(11:15) = (/ 6.0_r8, 5.4_r8, 5.5_r8, 5.2_r8, 5.0_r8 /)
fac(16:20) = (/ 5.0_r8, 5.0_r8, 5.0_r8, 5.0_r8, 5.0_r8 /)
fac(21:24) = (/ 5.0_r8, 5.0_r8, 5.0_r8, 5.0_r8 /)

vis_to_ir = 1.0_r8 / fac

tau = vis_to_ir(iv) * tau

end subroutine dust_vis_to_ir


  subroutine get_corrk(iv, v, p, t, ng, c, rc)
!------------------------------------------------------------------------
! subroutine get_corrk(iv, v, p, t, ng, c, rc)
!
! This subroutine is the one that actually interfaces with the correlated-k
!   tables.  It will take in a specified (wavenumber, pressure, temperature)
!   pairing and return the n_gauss k-coefficients for that bin -- it also
!   does some amount of interpolation within the tables
! This is largely modeled on the original TES forward operator code we 
!   received from Mike Smith et al. (get_corrk_LN_CO2.f90)

integer,  intent( in) :: iv, ng
real(r8), intent( in) :: v, p, t
real(r8), intent(out) :: c(ng)
integer,  intent(out) :: rc

integer  :: k, j, stat, p_ind, t_ind
logical  :: t_flag
real(r8) :: del_p, del_t, del_p_mb
real(r8), dimension(ng) :: vtp, vtp1, vt1p, vt1p1


! Set up logical flags to let us know if we're having problems with k-table
!   bound issues -- initially these are false
t_flag = .false.

! Check to make sure that pressure and temperature values are within the 
!   table bounds (we have already checked for wavenumber value back in main
!   get_expected_ routine).
! Pressure first
call tablesearch(np, ck_p, p, p_ind, del_p_mb, stat)
if ( stat == 1 ) then
   write(errstring, *)'Pressure is outside of lattice bounds -- setting rc to 1'
   call error_handler(E_MSG,'get_corrk', errstring, source, revision, revdate)
   rc = 1
   return
end if

! Assuming we've gotten beyond the above code, we can now find del_p in terms of
!   of scale height (log p) instead of mbars
del_p = ( log( p ) - log( ck_p(p_ind) ) )/( log( ck_p(p_ind+1) ) - log( ck_p(p_ind) ) )

! Now check temperature
call tablesearch(nt, ck_t, t, t_ind, del_t, stat)
if ( stat == 1 ) then
   if ( t < ck_t(1) ) then
      write(errstring, *)'Temperature is lower than lattice bound -- setting rc to 2'
      call error_handler(E_MSG,'get_corrk', errstring, source, revision, revdate)
      rc = 2
      return
   elseif ( t > ck_t(nt) ) then
      ! If we are in the dead center of the 15 um band, then we can abide by
      !   surface temperatures that are warmer than 280 K -- the surface
      !   contribution will not be measured from space so interpolation doesn't
      !   really matter.  Otherwise, return and complain.
      if ( v >= 657.0_r8 .and. v <= 678.5_r8 ) then
         t_flag = .true.
         t_ind = nt
      else
         write(errstring, *)'Temperature is higher than lattice bound -- setting rc to 3'
         call error_handler(E_MSG,'get_corrk', errstring, source, revision, revdate)
         rc = 3
         return
      end if
   else
      write(errstring, *)'We should not be able to get here -- setting rc to 4'
      call error_handler(E_MSG,'get_corrk', errstring, source, revision, revdate)
      rc = 4
      return
   end if
end if

! Time to find our location within k-table space and interpolate for the 
!   k-coefficients
! The TES team's k-tables have "-1" entries whereever they failed to include 
!   data within the tables.  The tables generated by Chris Lee do not have such
!   a placeholder -- all data have been filled in.
! Pre-store the potentially relevant entries in the array values
vtp   = values(:,t_ind,p_ind,iv)
vtp1  = values(:,t_ind,p_ind+1,iv)

! Now check whether we set t_flag to .true. above
if ( t_flag ) then
   c = exp( (1.0_r8 - del_p) * vtp + del_p * vtp1 )
else
   vt1p  = values(:,t_ind+1,p_ind,iv)
   vt1p1 = values(:,t_ind+1,p_ind+1,iv)

   ! Within get_corrk_LN_CO2, there is a long, scary bit where one evaluates lots and
   !   lots of if statements to one's place within the k-table space.  ALL good cases
   !   satisfy the condition that vtp and vt1p1 are not equal to -1.
   if ( vtp(1) /= -1.0_r8 .and. vt1p1(1) /= -1.0_r8 ) then

      ! Need to check further which acceptable case we are in to know how to
      !   properly interpolate
      if ( vt1p(1) == -1.0_r8 .and. vtp1(1) /= -1.0_r8 ) then         
         c = exp( (1.0_r8 - del_p) * vtp + (del_p - del_t) * vtp1 + del_t * vt1p1 )
      elseif ( vt1p(1) /= -1.0_r8 .and. vtp1(1) == -1.0_r8 ) then
         c = exp( (1.0_r8 - del_t) * vtp + (del_t - del_p) * vt1p + del_p * vt1p1 )
      else
         c = exp( (1.0_r8 - del_p) * (1.0_r8 - del_t) * vtp + &
                  (1.0_r8 - del_p) * del_t * vt1p + &
                  del_p * (1.0_r8 - del_t) * vtp1 + &
                  del_p * del_t * vt1p1 )
      end if
   else
      write(errstring, *)'Either vtp(1) or vt1p1(1) is equal to -1 -- setting rc to 5'
      call error_handler(E_MSG,'get_corrk', errstring, source, revision, revdate)
      rc = 5
      return
   end if

end if

! If we have returned already, then we were successful
rc = 0

end subroutine get_corrk


  subroutine tablesearch(n, table, x, ind, d_ind, status)
!------------------------------------------------------------------------------
! subroutine tablesearch(n, table, x, ind, d_ind, status)
! 
! This subroutine is basically a copy of Mike Smith's routine tablesearch_LN_CO2.f90
!   which is a part of the original TES retrieval code.
!
! Routine to determine if a value is within the bounds of the input table
! Tables must be ordered, either assending or descending
!
! Input:
!   integer n  - number of values in input table
!   real table - values to be search through
!   real x     - value to be found in table
!
! Output:
!   integer ind    - x is between table(ind) and table(ind+1)
!   real d_ind     - normalized distance from table(ind)
!   integer status - status of operation, 0 if found, 1 if out of bounds.

integer,  intent( in) :: n
real(r8), intent( in) :: table(n), x
integer,  intent(out) :: ind, status
real(r8), intent(out) :: d_ind

integer :: k

! Is value within bounds of table?
!   assume the code fails until we find that it doesn't
status = 1

! If x is greater than whole table return with status = 1
if ( (x > table(1)) .AND. (x > table(n)) ) then
  if (table(1) > table(n)) then
    ind = 1
  else
    ind = n
  end if
  d_ind = 0.0_r8
  return
end if
! If x is less than whole table return with status = 1
if ( (x < table(1)) .AND. (x < table(n)) ) then
  if (table(1) < table(n)) then
    ind = 1
  else
    ind = n
  end if
  d_ind = 0.0_r8
  return
end if

!Look for ind such that table(ind) <= x < table(ind+1)
!           *OR*        table(ind) >= x > table(ind+1)
do k = 1, n
  if ( ((x >= table(k)) .AND. (x < table(k+1)) ) .OR. &
       ((x <= table(k)) .AND. (x > table(k+1)) ) ) then
    ind = k
    d_ind = ( x - table(k) ) / ( table(k+1) - table(k) )
    status = 0 !success
    return
  end if
end do

end subroutine tablesearch

end module obs_def_TES_nadir_mod

! END DART PREPROCESS MODULE CODE

