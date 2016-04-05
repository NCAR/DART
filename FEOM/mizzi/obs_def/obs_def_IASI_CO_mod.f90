! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
!
! BEGIN DART PREPROCESS KIND LIST
! IASI_CO_RETRIEVAL, KIND_CO
! END DART PREPROCESS KIND LIST
!
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_iasi_co_mod, only : write_iasi_co, read_iasi_co, &
!                                  interactive_iasi_co, get_expected_iasi_co, &
!                                  set_obs_def_iasi_co
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(IASI_CO_RETRIEVAL)                                                           
!            call get_expected_iasi_co(state, location, obs_def%key, obs_val, istatus)  
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!
! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(IASI_CO_RETRIEVAL)
!         call read_iasi_co(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF
!
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(IASI_CO_RETRIEVAL)
!         call write_iasi_co(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF
!
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(IASI_CO_RETRIEVAL)
!         call interactive_iasi_co(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!
! BEGIN DART PREPROCESS SET_OBS_DEF_IASI_CO
!      case(IASI_CO_RETRIEVAL)
!         call set_obs_def_iasi_co(obs_def%key)
! END DART PREPROCESS SET_OBS_DEF_IASI_CO
!
! BEGIN DART PREPROCESS MODULE CODE
module obs_def_iasi_CO_mod
use        types_mod, only : r8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type, set_location, get_location, VERTISPRESSURE, VERTISSURFACE
use  assim_model_mod, only : interpolate
use    obs_kind_mod, only  : KIND_CO, KIND_SURFACE_PRESSURE
implicit none
public :: write_iasi_co, read_iasi_co, interactive_iasi_co, &
          get_expected_iasi_co, set_obs_def_iasi_co

! Storage for the special information required for observations of this type
integer, parameter               :: max_iasi_co_obs = 10000000
integer, parameter               :: iasi_dim = 19
integer, parameter               :: iasi_dimp = 20
integer                          :: num_iasi_co_obs = 0
real(r8), dimension(max_iasi_co_obs,19) :: avg_kernel
real(r8), dimension(max_iasi_co_obs,20) :: pressure
real(r8), dimension(max_iasi_co_obs)	:: iasi_prior
real(r8), dimension(max_iasi_co_obs)	:: iasi_psurf	
integer,  dimension(max_iasi_co_obs)    :: iasi_nlevels
integer,  dimension(max_iasi_co_obs)    :: iasi_nlevelsp
!
! For now, read in all info on first read call, write all info on first write call
logical :: already_read = .false., already_written = .false.
!
! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source: /home/thoar/CVS.REPOS/DART/obs_def/obs_def_iasi_CO_mod.f90,v $", &
revision = "$Revision$", &
revdate  = "$Date$"
logical, save :: module_initialized = .false.
integer  :: counts1 = 0
contains
!
!----------------------------------------------------------------------
subroutine initialize_module
!----------------------------------------------------------------------------
! subroutine initialize_module
call register_module(source, revision, revdate)
module_initialized = .true.
end subroutine initialize_module
!
subroutine read_iasi_co(key, ifile, fform)
!----------------------------------------------------------------------
!subroutine read_iasi_co(key, ifile, fform)
integer, intent(out)            :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional    :: fform
character(len=32) 		:: fileformat
integer			        :: iasi_nlevels_1
integer			        :: iasi_nlevelsp_1
real(r8)			:: iasi_prior_1
real(r8)			:: iasi_psurf_1
real(r8), dimension(iasi_dim)	:: avg_kernels_1
real(r8), dimension(iasi_dimp)	:: pressure_1
integer 			:: keyin
if ( .not. module_initialized ) call initialize_module
fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading
avg_kernels_1(:) = 0.0_r8
pressure_1(:) = 0.0_r8
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
   iasi_nlevels_1 = read_iasi_nlevels(ifile, fileformat)
   iasi_nlevelsp_1 = iasi_nlevels_1+1
   iasi_prior_1 = read_iasi_prior(ifile, fileformat)
   iasi_psurf_1 = read_iasi_psurf(ifile, fileformat)
   avg_kernels_1(1:iasi_nlevels_1) = read_iasi_avg_kernels(ifile, iasi_nlevels_1, fileformat)
   pressure_1(1:iasi_nlevelsp_1) = read_iasi_pressure(ifile, iasi_nlevelsp_1, fileformat)
   read(ifile) keyin
   CASE DEFAULT
   iasi_nlevels_1 = read_iasi_nlevels(ifile, fileformat)
   iasi_nlevelsp_1 = iasi_nlevels_1+1
   iasi_prior_1 = read_iasi_prior(ifile, fileformat)
   iasi_psurf_1 = read_iasi_psurf(ifile, fileformat)
   avg_kernels_1(1:iasi_nlevels_1) = read_iasi_avg_kernels(ifile, iasi_nlevels_1, fileformat)
   pressure_1(1:iasi_nlevelsp_1) = read_iasi_pressure(ifile, iasi_nlevelsp_1, fileformat)
   read(ifile, *) keyin
END SELECT
counts1 = counts1 + 1
key = counts1
call set_obs_def_iasi_co(key, avg_kernels_1, pressure_1, iasi_prior_1, iasi_psurf_1, &
                           iasi_nlevels_1, iasi_nlevelsp_1)
end subroutine read_iasi_co
!
subroutine write_iasi_co(key, ifile, fform)
!----------------------------------------------------------------------
!subroutine write_iasi_co(key, ifile, fform)
integer, intent(in)             :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional 	:: fform
character(len=32) 		:: fileformat
real(r8), dimension(iasi_dim)   :: avg_kernels_temp
real(r8), dimension(iasi_dimp)  :: pressure_temp
if ( .not. module_initialized ) call initialize_module
fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))
!
! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading
avg_kernels_temp=avg_kernel(key,:)
pressure_temp=pressure(key,:)
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
   call write_iasi_nlevels(ifile, iasi_nlevels(key), fileformat)
   call write_iasi_prior(ifile, iasi_prior(key), fileformat)
   call write_iasi_psurf(ifile, iasi_psurf(key), fileformat)
   call write_iasi_avg_kernels(ifile, avg_kernels_temp, iasi_nlevels(key), fileformat)
   call write_iasi_pressure(ifile, pressure_temp, iasi_nlevelsp(key), fileformat)
   write(ifile) key
   CASE DEFAULT
   call write_iasi_nlevels(ifile, iasi_nlevels(key), fileformat)
   call write_iasi_prior(ifile, iasi_prior(key), fileformat)
   call write_iasi_psurf(ifile, iasi_psurf(key), fileformat)
   call write_iasi_avg_kernels(ifile, avg_kernels_temp, iasi_nlevels(key), fileformat)
   call write_iasi_pressure(ifile, pressure_temp, iasi_nlevelsp(key), fileformat)
   write(ifile, *) key
END SELECT 
end subroutine write_iasi_co
!
subroutine interactive_iasi_co(key)
!----------------------------------------------------------------------
!subroutine interactive_iasi_co(key)
!
! Initializes the specialized part of a IASI observation
! Passes back up the key for this one
integer, intent(out) :: key
character(len=129) :: msgstring
if ( .not. module_initialized ) call initialize_module
!
! Make sure there's enough space, if not die for now (clean later)
if(num_iasi_co_obs >= max_iasi_co_obs) then
   ! PUT IN ERROR HANDLER CALL
   write(msgstring, *)'Not enough space for a iasi CO obs.'
   call error_handler(E_MSG,'interactive_iasi_co',msgstring,source,revision,revdate)
   write(msgstring, *)'Can only have max_iasi_co obs (currently ',max_iasi_co_obs,')'
   call error_handler(E_ERR,'interactive_iasi_co',msgstring,source,revision,revdate)
endif
!
! Increment the index
num_iasi_co_obs = num_iasi_co_obs + 1
key = num_iasi_co_obs
!
! Otherwise, prompt for input for the three required beasts
write(*, *) 'Creating an interactive_iasi_co observation'
write(*, *) 'Input the IASI Prior '
read(*, *) iasi_prior
write(*, *) 'Input IASI Surface Pressure '
read(*, *) iasi_psurf(num_iasi_co_obs)
write(*, *) 'Input the 19 Averaging Kernel Weights '
read(*, *) avg_kernel(num_iasi_co_obs,:)
write(*, *) 'Input the 20 Averaging Pressure Levels '
read(*, *) pressure(num_iasi_co_obs,:)
end subroutine interactive_iasi_co
!
subroutine get_expected_iasi_co(state, location, key, val, istatus)
!----------------------------------------------------------------------
!subroutine get_expected_iasi_co(state, location, key, val, istatus)
real(r8), intent(in)            :: state(:)
type(location_type), intent(in) :: location
integer, intent(in)             :: key
real(r8), intent(out)           :: val
integer, intent(out)            :: istatus
integer :: i,kstr
type(location_type) :: loc2
real(r8)            :: mloc(3)
real(r8)	    :: obs_val,wrf_psf,level,missing
real(r8)            :: co_min,iasi_prs_mid,iasi_psf,iasi_psf_save
integer             :: nlevels,nlevelsp,nnlevels
integer             :: iflg
character(len=129)  :: msgstring
!
! Initialize DART
if ( .not. module_initialized ) call initialize_module
!
! Initialize variables
co_min=1.e-4
!co_min=-9.
missing=-888888.0_r8
!
! Get iasi data
nlevels = iasi_nlevels(key)
nlevelsp = iasi_nlevelsp(key)
iasi_psf = iasi_psurf(key)
iasi_psf_save = iasi_psurf(key)
!
! Get location infomation
mloc = get_location(location)
if (mloc(2)>90.0_r8) then
    mloc(2)=90.0_r8
elseif (mloc(2)<-90.0_r8) then
    mloc(2)=-90.0_r8
endif
!
! Get wrf surface pressure
wrf_psf = 0.0_r8
istatus = 0
loc2 = set_location(mloc(1), mloc(2), 0.0_r8, VERTISSURFACE)
call interpolate(state, loc2, KIND_SURFACE_PRESSURE, wrf_psf, istatus)  
!write(msgstring, *)'APM ERROR: wrf_psf, iasi_psf, status ',wrf_psf,iasi_psf,istatus 
!call error_handler(E_MSG,'set_obs_def_iasi_co',msgstring,source,revision,revdate)
!
! Correct iasi surface pressure
if(iasi_psf.gt.wrf_psf) then
   iasi_psf=wrf_psf
endif
!
! Find kstr - the surface level index
kstr=0
do i=1,iasi_dim
   if (i.eq.1 .and. iasi_psf.gt.pressure(key,2)) then
      kstr=1
      exit
   endif
   if (i.ne.1 .and. i.ne.iasi_dim .and. iasi_psf.le.pressure(key,i) .and. &
      iasi_psf.gt.pressure(key,i+1)) then
      kstr=i
      exit   
   endif
enddo
if (kstr.eq.0) then
   write(msgstring, *)'APM: ERROR in IASI obs def kstr=0: iasi_psf=',iasi_psf
   call error_handler(E_MSG,'set_obs_def_iasi_co',msgstring,source,revision,revdate)
   stop
elseif (kstr.gt.6) then
   write(msgstring, *)'APM: ERROR surface pressure is unrealistic: iasi_psf=',iasi_psf
   call error_handler(E_MSG,'set_obs_def_iasi_co',msgstring,source,revision,revdate)
   stop
endif
!
! Reject ob when number of IASI levels from WRF cannot equal actual number of IASI levels
nnlevels=iasi_dim-kstr+1
if(nnlevels.ne.nlevels) then
   obs_val=missing
   istatus=2
   write(msgstring, *)'APM: NOTICE reject ob - # of WRF IASI levels .ne. # of IASI levels  '
   call error_handler(E_MSG,'set_obs_def_iasi_co',msgstring,source,revision,revdate)
   return
endif   
!
! Find the lowest pressure level midpoint
iasi_prs_mid=(iasi_psf+pressure(key,kstr+1))/2.
!
! Apply IASI Averaging kernel A and IASI Prior (I-A)xa
! x = Axm + (I-A)xa , where x is a 10 element vector 
val = 0.0_r8
do i=1,nlevels
!
! APM: remove the if test to use layer average data
   if (i .eq.1) then
      loc2 = set_location(mloc(1),mloc(2),iasi_prs_mid, VERTISPRESSURE)
   else
      iasi_prs_mid=(pressure(key,kstr+i-1)+pressure(key,kstr+i))/2.
      loc2 = set_location(mloc(1),mloc(2),iasi_prs_mid, VERTISPRESSURE)
   endif
!
! Interpolate WRF CO data to IASI pressure level midpoint
   obs_val = 0.0_r8
   istatus = 0
   call interpolate(state, loc2, KIND_CO, obs_val, istatus)  
   if (istatus /= 0) then
!      write(msgstring, *)'APM ERROR: istatus,kstr,obs_val ',istatus,kstr,obs_val 
!      call error_handler(E_MSG,'set_obs_def_iasi_co',msgstring,source,revision,revdate)
!      write(msgstring, *)'APM ERROR: iasi_prs_mid ',iasi_prs_mid
!      call error_handler(E_MSG,'set_obs_def_iasi_co',msgstring,source,revision,revdate)
!      write(msgstring, *)'APM ERROR: wrf_psf,iasi_psurf,iasi_psf ', &
!      wrf_psf,iasi_psf_save,iasi_psf
!      call error_handler(E_MSG,'set_obs_def_iasi_co',msgstring,source,revision,revdate)
!      write(msgstring, *)'APM ERROR: i, nlevels ',i,nlevels
!      call error_handler(E_MSG,'set_obs_def_iasi_co',msgstring,source,revision,revdate)
      write(msgstring, *)'APM NOTICE: WRF extrapolation needed reject ob '
      call error_handler(E_MSG,'set_obs_def_iasi_co',msgstring,source,revision,revdate)
      return
!      stop
   endif
!
! Check for WRF CO lower bound
   if (obs_val.lt.co_min) then
      obs_val=co_min
      write(msgstring, *)'APM NOTICE: in obs_def_iasi resetting minimum co value '
      call error_handler(E_MSG,'set_obs_def_iasi_co',msgstring,source,revision,revdate)
   endif
!
! apply averaging kernel
!   val = val + avg_kernel(key,i) * log10(obs_val*1.e-6)  
   val = val + avg_kernel(key,i) * obs_val*1.e3  
!   val = val + avg_kernel(key,i) * obs_val  
!   write(msgstring, *)'APM DATA: i,iasi_prs_mid,val ',i,iasi_prs_mid,val
!   call error_handler(E_MSG,'set_obs_def_iasi_co',msgstring,source,revision,revdate)
!   write(msgstring, *)'APM DATA: i,obs_val,avg_kernel ',i,obs_val,avg_kernel(key,i)
!   call error_handler(E_MSG,'set_obs_def_iasi_co',msgstring,source,revision,revdate)
enddo
!write(msgstring, *)'APM DATA: iasi_prior ',i,obs_val,iasi_prior(key)
!call error_handler(E_MSG,'set_obs_def_iasi_co',msgstring,source,revision,revdate)
!
val = val + iasi_prior(key)
end subroutine get_expected_iasi_co
!
!----------------------------------------------------------------------
!
subroutine set_obs_def_iasi_co(key, co_avgker, co_press, co_prior, co_psurf, co_nlevels, co_nlevelsp)
!----------------------------------------------------------------------
! Allows passing of obs_def special information 
integer,	 	intent(in)	:: key, co_nlevels, co_nlevelsp
real*8,dimension(19),	intent(in)	:: co_avgker	
real*8,dimension(20),	intent(in)	:: co_press	
real*8,			intent(in)	:: co_prior
real*8,			intent(in)	:: co_psurf
character(len=129) 			:: msgstring
if ( .not. module_initialized ) call initialize_module
if(num_iasi_co_obs >= max_iasi_co_obs) then
   ! PUT IN ERROR HANDLER CALL
   write(msgstring, *)'Not enough space for a iasi CO obs.'
   call error_handler(E_MSG,'set_obs_def_iasi_co',msgstring,source,revision,revdate)
   write(msgstring, *)'Can only have max_iasi_co_obs (currently ',max_iasi_co_obs,')'
   call error_handler(E_ERR,'set_obs_def_iasi_co',msgstring,source,revision,revdate)
endif
avg_kernel(key,:) 	= co_avgker(:)
pressure(key,:) 	= co_press(:)
iasi_prior(key)   	= co_prior
iasi_psurf(key)	        = co_psurf
iasi_nlevels(key)       = co_nlevels
iasi_nlevelsp(key)      = co_nlevelsp
end subroutine set_obs_def_iasi_co
!
function read_iasi_prior(ifile, fform)
integer,                    intent(in) :: ifile
real(r8)                               :: read_iasi_prior
character(len=*), intent(in), optional :: fform
character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat
if ( .not. module_initialized ) call initialize_module
fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_prior
   CASE DEFAULT
      read(ifile, *) read_iasi_prior
END SELECT
end function read_iasi_prior
!
function read_iasi_nlevels(ifile, fform)
integer,                    intent(in) :: ifile
integer                               :: read_iasi_nlevels
character(len=*), intent(in), optional :: fform
character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat
if ( .not. module_initialized ) call initialize_module
fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_nlevels
   CASE DEFAULT
      read(ifile, *) read_iasi_nlevels
END SELECT
end function read_iasi_nlevels
!
function read_iasi_nlevelsp(ifile, fform)
integer,                    intent(in) :: ifile
integer                               :: read_iasi_nlevelsp
character(len=*), intent(in), optional :: fform
character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat
if ( .not. module_initialized ) call initialize_module
fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_nlevelsp
   CASE DEFAULT
      read(ifile, *) read_iasi_nlevelsp
END SELECT
end function read_iasi_nlevelsp
!
subroutine write_iasi_prior(ifile, iasi_prior_temp, fform)
integer,                    intent(in) :: ifile
real(r8), 		    intent(in) :: iasi_prior_temp
character(len=32),          intent(in) :: fform
character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat
if ( .not. module_initialized ) call initialize_module
fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) iasi_prior_temp
   CASE DEFAULT
      write(ifile, *) iasi_prior_temp
END SELECT
end subroutine write_iasi_prior
!
subroutine write_iasi_nlevels(ifile, iasi_nlevels_temp, fform)
integer,                    intent(in) :: ifile
integer,                    intent(in) :: iasi_nlevels_temp
character(len=32),          intent(in) :: fform
character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat
if ( .not. module_initialized ) call initialize_module
fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) iasi_nlevels_temp
   CASE DEFAULT
      write(ifile, *) iasi_nlevels_temp
END SELECT
end subroutine write_iasi_nlevels
!
subroutine write_iasi_nlevelsp(ifile, iasi_nlevelsp_temp, fform)
integer,                    intent(in) :: ifile
integer,                    intent(in) :: iasi_nlevelsp_temp
character(len=32),          intent(in) :: fform
character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat
if ( .not. module_initialized ) call initialize_module
fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) iasi_nlevelsp_temp
   CASE DEFAULT
      write(ifile, *) iasi_nlevelsp_temp
END SELECT
end subroutine write_iasi_nlevelsp
!
function read_iasi_psurf(ifile, fform)
integer,                    intent(in) :: ifile
real(r8)                               :: read_iasi_psurf
character(len=*), intent(in), optional :: fform
character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat
if ( .not. module_initialized ) call initialize_module
fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_psurf
   CASE DEFAULT
      read(ifile, *) read_iasi_psurf
END SELECT
end function read_iasi_psurf
!
subroutine write_iasi_psurf(ifile, iasi_psurf_temp, fform)
integer,                    intent(in) :: ifile
real(r8),		    intent(in) :: iasi_psurf_temp
character(len=32),          intent(in) :: fform
character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat
if ( .not. module_initialized ) call initialize_module
fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) iasi_psurf_temp
   CASE DEFAULT
      write(ifile, *) iasi_psurf_temp
END SELECT
end subroutine write_iasi_psurf
!
function read_iasi_avg_kernels(ifile, nlevels, fform)
integer,                    intent(in) :: ifile, nlevels
real(r8), dimension(19)        :: read_iasi_avg_kernels
character(len=*), intent(in), optional :: fform
character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat
read_iasi_avg_kernels(:) = 0.0_r8
if ( .not. module_initialized ) call initialize_module
fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_avg_kernels(1:nlevels)
   CASE DEFAULT
      read(ifile, *) read_iasi_avg_kernels(1:nlevels)
END SELECT
end function read_iasi_avg_kernels
!
subroutine write_iasi_avg_kernels(ifile, avg_kernels_temp, nlevels_temp, fform)
integer,                    intent(in) :: ifile, nlevels_temp
real(r8), dimension(19), intent(in)  :: avg_kernels_temp
character(len=32),          intent(in) :: fform
character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat
if ( .not. module_initialized ) call initialize_module
fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) avg_kernels_temp(1:nlevels_temp)
   CASE DEFAULT
      write(ifile, *) avg_kernels_temp(1:nlevels_temp)
END SELECT
end subroutine write_iasi_avg_kernels
!
function read_iasi_pressure(ifile, nlevelsp, fform)
integer,                    intent(in) :: ifile, nlevelsp
real(r8), dimension(20)        :: read_iasi_pressure
character(len=*), intent(in), optional :: fform
character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat
read_iasi_pressure(:) = 0.0_r8
if ( .not. module_initialized ) call initialize_module
fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_pressure(1:nlevelsp)
   CASE DEFAULT
      read(ifile, *) read_iasi_pressure(1:nlevelsp)
END SELECT
end function read_iasi_pressure
!
subroutine write_iasi_pressure(ifile, pressure_temp, nlevelsp_temp, fform)
integer,                    intent(in) :: ifile, nlevelsp_temp
real(r8), dimension(20), intent(in)  :: pressure_temp
character(len=32),          intent(in) :: fform
character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat
if ( .not. module_initialized ) call initialize_module
fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) pressure_temp(1:nlevelsp_temp)
   CASE DEFAULT
      write(ifile, *) pressure_temp(1:nlevelsp_temp)
END SELECT
end subroutine write_iasi_pressure
!
end module obs_def_iasi_CO_mod
! END DART PREPROCESS MODULE CODE
