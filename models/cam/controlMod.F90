! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

#include <misc.h>
#include <preproc.h>

module controlMod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar        !parameter statements 
  use clm_varctl        !run control variables 
  use histFileMod       !history file variables
  use spmdMod           !spmd routines and variables

  implicit none


! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

  save

! Namelist variables only used locally

  integer :: hist_ndens                         !output density of netcdf history files
  logical :: mkfsurdat                          !true => make surface data from raw data
  character(len=256) :: rpntpath                !full UNIX pathname of restart pointer file
  character(len=  7) :: runtyp(4)               !run type
  character(len  =8) :: hist_fldaux1(maxalflds) !fields for first  auxillary history file
  character(len  =8) :: hist_fldaux2(maxalflds) !fields for second auxillary history file
#if (defined COUP_CSM)
  integer :: csm_dtime  !value passed to coupler on initialization set to namelist input 
                        !consistency check done later that restart file gives same value
#endif

!=======================================================================
CONTAINS
!=======================================================================

  subroutine control_init (cam_caseid , cam_ctitle, cam_irad , cam_nsrest, &
                           cam_crtinic, cam_nhtfrq, cam_mfilt, cam_irt )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! initialize run control variables 
! 
! Method: 
! When running in cam mode, the base calendar info, nstep, nestep, 
! nsrest, and time step are input to the land model from CAM. 
! The values in the clmexp namelist are not used. The minimum 
! namelist parameters are:
!    o fsurdat
!    o finidat
!    o fpftcon 
! When running in offline or csm mode, the minimum namelist parameters are:
!    o nsrest
!    o nestep or nelapse
!    o fsurdat
!    o finidat
!    o dtime
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------

#if (defined OFFLINE) || (defined COUP_CSM)
  use time_manager, only : calendar, dtime, nestep, nelapse, start_ymd, &
                           start_tod, stop_ymd, stop_tod, ref_ymd, ref_tod
#else
  use time_manager, only : get_step_size
#endif

! ------------------------ includes ------------------------------------
    include 'netcdf.inc'
! ----------------------------------------------------------------------

! ------------------------ arguments -----------------------------------
    character(len=*), optional, intent(in) :: cam_caseid    !cam caseid
    character(len=*), optional, intent(in) :: cam_ctitle    !cam title
    integer         , optional, intent(in) :: cam_irad      !cam radiation frequency
    integer         , optional, intent(in) :: cam_nsrest    !cam run type
    character(len=*), optional, intent(in) :: cam_crtinic   !cam initial dataset frequency
    integer         , optional, intent(in) :: cam_nhtfrq    !cam history write freq for tape 1
    integer         , optional, intent(in) :: cam_mfilt     !cam number of files per tape for tape 1
    integer         , optional, intent(in) :: cam_irt       !cam mss retention time
! ----------------------------------------------------------------------

! ------------------------ local variables -----------------------------
    character(len=256) :: homedir   !full UNIX filepath name of home directory
    character(len=256) :: logid     !logid part of file path name
    character(len=256) :: cap       !upper case logid
    character(len=  1) :: ctmp      !character temporary
    integer :: i,j,n                !loop indices
    integer :: iundef               !integer undefined value
    real(r8):: rundef               !real undefined value
    integer :: ierr                 !error code
! -----------------------------------------------------------------


#if (defined COUP_CAM)
! The following can be read in but are overwritten with values from the
! cam time_manager module - consequently they are only declared as local 
! variables here

    integer :: dtime      ! timestep in seconds
    integer :: nestep     ! final timestep (or day if negative) number
    integer :: nelapse    ! number of timesteps (or days if negative) to extend a run
    integer :: start_ymd  ! starting date for run in yearmmdd format
    integer :: start_tod  ! starting time of day for run in seconds
    integer :: stop_ymd   ! stopping date for run in yearmmdd format
    integer :: stop_tod   ! stopping time of day for run in seconds
    integer :: ref_ymd    ! reference date for time coordinate in yearmmdd format
    integer :: ref_tod    ! reference time of day for time coordinate in seconds
    character(len=32) :: calendar ! Calendar in date calculations ('NO_LEAP' or 'GREGORIAN')
#endif    
! ----------------------------------------------------------------------

! ------------------------ namelist variables --------------------------
    namelist /clmexp/  &
         ctitle, caseid, nsrest,  &
         calendar, dtime, nelapse, nestep, start_ymd, start_tod,  &
         stop_ymd, stop_tod, ref_ymd, ref_tod, &
         nrevsn, rpntpath, hist_ndens, hist_dov2xy, &
         hist_nhtfrq, hist_mfilt, hist_fldadd, &
         hist_chntyp, hist_fldaux1, hist_fldaux2, hist_crtinic, &
         archive_dir, mss_wpass, mss_irt, &
         finidat, fsurdat, fpftcon, frivinp_rtm, offline_atmdir, &
         mksrf_fvegtyp, mksrf_fsoitex, mksrf_fsoicol, mksrf_flanwat, &
         mksrf_fglacier, mksrf_furban, mksrf_flai, mksrf_offline_fgrid, &
         mksrf_offline_edgen, mksrf_offline_edgee, mksrf_offline_edges, &
         mksrf_offline_edgew, mksrf_offline_fnavyoro, &
         conchk, irad, wrtdia, csm_doflxave, rtm_nsteps 
         

! === define run =======================
!
!    o caseid     = 32 character case name
!    o ctitle     = 80 character case title
!    o nsrest     = integer flag. 0: initial run. 1: restart: 3: branch
!
! === model time =======================
!
!    o dtime      = real model time step (s)
!    o calendar   = Calendar to use in date calculations. 
!                  'no_leap' (default) or 'gregorian'
!    o start_ymd  = Starting date for run encoded in yearmmdd format.  
!                   Default value is read from initial conditions file.
!    o start_tod  = Starting time of day for run in seconds since 0Z.  
!                   Default value is read from initial conditions file.
!    o stop_ymd   = Stopping date for run encoded in yearmmdd format.  
!                   No default.
!    o stop_tod   = Stopping time of day for run in seconds since 0Z.  
!                   Default: 0.
!    o nelapse    = nnn, Specify the ending time for the run as an interval
!                   starting at the current time in either timesteps
!                   (if positive) or days (if negative).
!                   Either nestep or (stop_ymd,stop_tod) take precedence.
!    o nestep     = nnnn, Specify the ending time for the run as an interval
!                   starting at (start_ymd,start_tod) in either timesteps
!                   (if positive) or days (if negative).
!                   (stop_ymd,stop_tod) takes precedence if set.
!    o ref_ymd    = Reference date for time coordinate encoded in yearmmdd format.
!                   Default value is start_ymd.
!    o ref_tod    = Reference time of day for time coordinate in seconds since 0Z.
!                   Default value is start_tod.
!
! === input data ===
!
!    o finidat         = 256 character initial conditions file name
!    o fsurdat         = 256 character surface data file name
!    o fpftcon         = 256 character data file with PFT physiological constants
!    o frivinp_rtm     = 256 character input data file for rtm
!    o nrevsn          = 256 character restart file name for use with branch run
!
! === offline forcing data ===
!
!    o offline_atmdir  = 256 character directory for input atm data files (can be Mass Store)
!
! === input data when making surface data [fsurdat] ===
!
!    o mksrf_offline_fgrid   = offline - land grid dataset to use instead of generating grid 
!    o mksrf_offline_fnavyoro= offline - 20 min navy orography dataset
!    o mksrf_offline_edgen   = offline - northern edge of grid (degrees): >  -90 and <= 90
!    o mksrf_offline_edgee   = offline - eastern edge of grid (degrees) : see following notes
!    o mksrf_offline_edges   = offline - southern edge of grid (degrees): >= -90 and <  90
!    o mksrf_offline_edgew   = offline - western edge of grid (degrees) : see following notes
!    o mksrf_fvegtyp         = 256 character vegetation type data file name
!    o mksrf_fsoitex         = 256 character soil texture data file name
!    o mksrf_fsoicol         = 256 character soil color data file name
!    o mksrf_flanwat         = 256 character inland water data file name
!    o mksrf_furban          = 256 character urban data file name
!    o mksrf_fglacier        = 256 character glacier data file name
!    o mksrf_flai            = 256 character lai data file file name
!
! === history and restart files ===
!
!    o hist_ndens    = integer, can have value of 1 (nc_double) or 2 (nf_float) 
!    o hist_dov2xy   = true if want grid-average history field (false = vector)
!    o hist_nhtfrq   = integer history interval (+ = iterations,  - = hours, 0=monthly ave)
!    o hist_mfilt    = integer number of time samples per history file
!    o hist_fldadd   = 8 character name of fields to change to active 
!    o hist_chntyp   = paired 8 character field name and field type. 
!                      OVERRIDES default settings in routine histlst(e.g., 'TV','maximum')
!    o hist_fldaux1  = 8 character name of fields for first  auxillary history file
!    o hist_fldaux2  = 8 character name of fields for second auxillary history file
!    o hist_crtinic  = 8 character frequency to generate initial dataset
!                      can be set to 'MONTHLY', 'YEARLY' or 'NONE'.
! kdr                  added 'ENDOFRUN' option for data assim, as in cam source
!    o rpntpath      = 256 character full UNIX pathname of the local restart pointer 
!                      file. This file must exist when the model is restarted. This 
!                      file is overwritten and updated every time new restart data 
!                      files are output. 
!
! === long term archiving =====
!
!    o archive_dir = 256 character long term archive directory (can be MSS directory)
!    o mss_irt     = integer mass store retention period (days)
!    o mss_wpass   = 8 character mass store write password for output data sets
!
! === model physics ===
!
!    o conchk     = true if want error energy and water conservation checks
!    o irad       = integer solar radiation frequency (+ = iteration. - = hour)
!    o wrtdia     = true if want output written
!    o csm_doflxave = true => flux averaging is to be performed (only used for csm mode)
!
! === rtm control variables ===
!
!    o rtm_nsteps  = if > 1, average rtm over rtm_nsteps time steps  
!
! ----------------------------------------------------------------------

    if (masterproc) then          
       write(6,*) 'Attempting to initialize run control settings .....'
    endif

    runtyp(0 + 1) = 'initial'
    runtyp(1 + 1) = 'restart'
    runtyp(3 + 1) = 'branch '

    iundef = -9999999
    rundef = -9999999.

! ----------------------------------------------------------------------
! Default values
! ----------------------------------------------------------------------

! control variables

    caseid  = ' '
    ctitle  = ' '
    nsrest  = iundef

! initial data

    fsurdat = ' '
    finidat = ' '
    fpftcon = ' '
    frivinp_rtm = ' '
    nrevsn  = ' '

! offline mode

    offline_atmdir   = ' '

! surface generation 

    mksrf_offline_fgrid    = ' '
    mksrf_offline_fnavyoro = ' '
    mksrf_offline_edgen    =   90.
    mksrf_offline_edgee    =  180.
    mksrf_offline_edges    =  -90.
    mksrf_offline_edgew    = -180.

    mksrf_fvegtyp  = ' '
    mksrf_fsoitex  = ' '
    mksrf_fsoicol  = ' '
    mksrf_flanwat  = ' '
    mksrf_furban   = ' '
    mksrf_fglacier = ' '
    mksrf_flai     = ' '

! long term archive settings

    archive_dir = ' '
    mss_irt = 0
    mss_wpass = ' '
    
! history file variables

    hist_ndens  = 1	
    hist_dov2xy(1) = .true.
    hist_dov2xy(2:maxhist) = .true.
    hist_nhtfrq(1) = -24
    hist_nhtfrq(2:maxhist) = iundef
    hist_mfilt(1) = 1
    hist_mfilt(2:maxhist)  = iundef
    hist_fldadd(:) = ' '
    hist_chntyp(:,:) = ' '
    hist_fldaux1(:) = ' '
    hist_fldaux2(:) = ' '
    hist_crtinic = 'YEARLY'
    rpntpath = 'not_specified' 

! other namelist variables

    irad = -1
    conchk = .true.
    wrtdia = .false.
    csm_doflxave = .true.

! RTM control variables 

#if (defined RTM)
    rtm_nsteps = (3600*3)/dtime ! 3 hours
#endif

! ----------------------------------------------------------------------
! Read namelist from standard input. Override if coupled to CAM
! ----------------------------------------------------------------------

    if (masterproc) then
       read(5, clmexp, iostat=ierr)
       if (ierr /= 0) then
          if (masterproc) then
             write(6,*)'error: namelist input resulted in error code ',ierr
          endif
          call endrun
       endif
    endif

#if (defined COUP_CAM)
    caseid  = cam_caseid
    ctitle  = cam_ctitle 
    irad    = cam_irad
    nsrest  = cam_nsrest
    hist_crtinic  = cam_crtinic
#if (defined PERGRO)    
    hist_nhtfrq(1) = -48
#else
    hist_nhtfrq(1) = cam_nhtfrq
#endif
    hist_mfilt(1) = cam_mfilt
    mss_irt = cam_irt    
#endif

#if (defined OFFLINE) 
    if (fsurdat == ' ') then
       if (mksrf_offline_fnavyoro /= ' ' .and. mksrf_offline_fgrid /= ' ') then
          if (masterproc) then
             write(6,*) 'cannot set both MKSRF_OFFLINE_FNAVYORO and MKSRF_OFFLINE_FGRID'
          endif
          call endrun
       endif
       if (mksrf_offline_fgrid /= ' ') then  ! must have global grid (overwrite namelist values)
          mksrf_offline_edgen =   90.
          mksrf_offline_edgee =  180.
          mksrf_offline_edges =  -90.
          mksrf_offline_edgew = -180.
       endif
    endif
#endif

!if archive directory not input in namelist - set default from caseid

    if (archive_dir == ' ') then
       logid  = ' '
       call getenv ('LOGNAME', logid)
       if (logid(1:1) == ' ') then
          write (6,*) 'error: logname not defined'
          call endrun
       end if
       cap = ' '
       do i = 1, len_trim(logid)
          cap(i:i) = logid(i:i)
          ctmp = cap(i:i)
          if (ichar(logid(i:i))>=97 .and. ichar(logid(i:i))<=122) then
             cap(i:i) = char(ichar(ctmp) - 32)
          endif
       end do
       archive_dir = '/' // trim(cap) // '/csm/' // trim(caseid) // '/lnd' 
    end if
       
#if (defined SPMD)
    call control_spmd()
#endif

! ----------------------------------------------------------------------
! Define run 
! ----------------------------------------------------------------------

! Determine run type

    if (nsrest == iundef) then
       if (masterproc) write(6,*) 'error: must set nsrest'
       call endrun
    end if

! ----------------------------------------------------------------------
! Surface data
! ----------------------------------------------------------------------

    if (fsurdat == ' ') then
       mkfsurdat = .true.
    else
       mkfsurdat = .false.
    endif

    if (mkfsurdat) then

       if (mksrf_fvegtyp  == ' ' .or. &
           mksrf_fsoitex  == ' ' .or. &
           mksrf_fsoicol  == ' ' .or. &
           mksrf_flanwat  == ' ' .or. &
           mksrf_furban   == ' ' .or. &
           mksrf_fglacier == ' ' .or. &
           mksrf_flai     == ' ') then
          if (masterproc) then
             write(6,*) 'error: need to set data file name'
             write(6,*)'mksrf_fvegtyp = ',mksrf_fvegtyp
             write(6,*)'mksrf_fsoitex = ',mksrf_fsoitex
             write(6,*)'mksrf_fsoicol = ',mksrf_fsoicol
             write(6,*)'mksrf_flanwat = ',mksrf_flanwat
             write(6,*)'mksrf_furban  = ',mksrf_furban
             write(6,*)'mksrf_fglacier= ',mksrf_fglacier
             write(6,*)'mksrf_flai    = ',mksrf_flai
          endif
          call endrun
       end if

       if (nsrest > 0) then
          if (masterproc) then
             write(6,*) 'error: can not make surface data ', &
                  'during a continuation run'
          endif
          call endrun
       end if

       if (finidat /= ' ') then
          if (masterproc) then
             write(6,*) 'error: can not make surface data ', &
                  'when finidat is already specified'
             write(6,*) 'set finidat to empty string in namelist' 
          endif
          call endrun
       end if
       
    end if

#if (defined OFFLINE)
! ----------------------------------------------------------------------
! Offline data
! ----------------------------------------------------------------------

    if (masterproc .and. offline_atmdir == ' ') then 
       if (masterproc) then
          write(6,*)'error: atmos  input data file must be specified'
       endif
       call endrun
    endif

#endif

! ----------------------------------------------------------------------
! Model physics
! ----------------------------------------------------------------------

#if (defined COUP_CAM)
!time manager initialization done in cam code
    dtime = get_step_size() 
#elif (defined COUP_CSM)
!upon restart dtime not read in from restart file until after
!initial information sent to coupler - need to set it to namelist
!value here and do consistency check in routine restrd.
    csm_dtime = dtime
#endif

    if (irad < 0) irad = nint(-irad*3600./dtime)

#if (defined COUP_CSM)
    if (csm_doflxave .and. irad ==1 ) then
       if (masterproc) then
          write(6,*)'error: irad must be greater that one if', &
            ' flux averaging option is enabled'
       endif
       call endrun
    endif
#endif

! ----------------------------------------------------------------------
! History and restart files 
! ----------------------------------------------------------------------

    mss_irt = min(mss_irt,1825)

    fldaux(:,:) = ' '
    if (maxhist-1 > 3) then
       if (masterproc) write(6,*) 'CLM_CTLI error: must create additional fldaux'
       call endrun
    end if
    do j = 1, maxhist-1
       do i = 1, maxalflds
          if (j == 1) fldaux(i,j) = hist_fldaux1(i)
          if (j == 2) fldaux(i,j) = hist_fldaux2(i)
       end do
    end do

    nhist = 0
    do j = 1, maxhist-1
       do i = 1, maxalflds
          if (fldaux(i,j) /= ' ') nhist = j
       end do
    end do
    nhist = nhist + 1

    do i = 1, nhist
       if (hist_mfilt(i) == iundef ) then
          if (masterproc) then
             write(6,*)'error: must set hist_mfilt for file ',i
          endif
          call endrun
       end if
       if (hist_nhtfrq(i) == iundef) then
          if (masterproc) then
             write(6,*)'error: must set hist_nhtfrq for file ',i
          endif
          call endrun
       else if (hist_nhtfrq(i) < 0) then
          hist_nhtfrq(i) = nint(-hist_nhtfrq(i)*SHR_CONST_CDAY/(24.*dtime))
       endif
    end do

    if (hist_ndens == 1) then
       ncprec = nf_double
    else if (hist_ndens == 2) then
       ncprec = nf_float
    else
       if (masterproc) then
          write(6,*)'error: history tape hist_ndens must be 1 or 2'
       endif
       call endrun
    end if

    if (rpntpath == 'not_specified') then
       call getenv ('HOME', homedir)
       rpntpath = trim(homedir)//'/lnd.'//trim(caseid)//'.rpointer'
    endif

    do i = 1, nhist
       if (hist_nhtfrq(i)==0) then
          hist_mfilt(i) = 1
       endif
    end do

    if (nsrest == 0) nrevsn = ' '
    if (nsrest == 1) nrevsn = 'set by restart pointer file file'
    if (nsrest == 3 .and. nrevsn == ' ') then
       if (masterproc) write(6,*) 'error: need to set restart data file name' 
       call endrun
    end if

!kdr orig;    if (trim(hist_crtinic) /= 'MONTHLY' .and. trim(hist_crtinic) /= 'YEARLY') then 
    if (trim(hist_crtinic) /= 'MONTHLY' .and. trim(hist_crtinic) /= 'YEARLY' &
        .and. trim(hist_crtinic) /= 'ENDOFRUN') then 
       hist_crtinic = 'NONE'
    endif

! ----------------------------------------------------------------------
! Restart pointer file
! ----------------------------------------------------------------------

! split the full pathname of the restart pointer file into a 
! directory name and a file name
! check if the directory exists and if not, make it 

    rpntdir = ' '
    rpntfil = ' '
    do n = len_trim(rpntpath),0,-1     
       if (rpntpath(n:n) ==  '/') then
          rpntdir = rpntpath(1:n-1)
          rpntfil = rpntpath(n+1:len_trim(rpntpath))
          go to 100
       endif
    enddo
    rpntdir = '.'        ! no "/" found, set path = "."
    rpntfil = rpntpath   ! no "/" found, use whole input string.
100 continue

    if (masterproc) then
       write(6,*) 'Successfully initialized run control settings'
       write(6,*)
    endif

    return
  end subroutine control_init

!=======================================================================

#if (defined SPMD)

  subroutine control_spmd()

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Distribute namelist data all processors. The cpp SPMD definition 
! provides for the funnelling of all program i/o through the master 
! processor. Processor 0 either reads restart/history data from the 
! disk and distributes it to all processors, or collects data from 
! all processors and writes it to disk.
! 
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

#if (defined OFFLINE) || (defined COUP_CSM)
    use time_manager, only : calendar, dtime, nestep, nelapse, start_ymd, &
                             start_tod, stop_ymd, stop_tod, ref_ymd, ref_tod
#endif
    use mpishorthand

! ------------------------ local variables -----------------------------
    integer ier       !error code
!-----------------------------------------------------------------------

! run control variables

    call mpi_bcast (caseid, len(caseid), mpichar, 0, mpicom, ier)
    call mpi_bcast (ctitle, len(ctitle), mpichar, 0, mpicom, ier)
    call mpi_bcast (nsrest, 1, mpiint, 0, mpicom, ier)

#if (defined OFFLINE) || (defined COUP_CSM)
    call mpi_bcast (nestep, 1, mpiint, 0, mpicom, ier)
    call mpi_bcast (nelapse, 1, mpiint, 0, mpicom, ier)
    call mpi_bcast (dtime, 1, mpiint, 0, mpicom, ier)
    call mpi_bcast (calendar, 32, mpichar, 0, mpicom, ier)
    call mpi_bcast (start_ymd, 1, mpiint, 0, mpicom, ier)
    call mpi_bcast (start_tod, 1, mpiint, 0, mpicom, ier)
    call mpi_bcast (stop_ymd, 1, mpiint, 0, mpicom, ier)
    call mpi_bcast (stop_tod, 1, mpiint, 0, mpicom, ier)
    call mpi_bcast (ref_ymd, 1, mpiint, 0, mpicom, ier)
    call mpi_bcast (ref_tod, 1, mpiint, 0, mpicom, ier)
#endif

! initial file variables

    call mpi_bcast (nrevsn, len(nrevsn), mpichar, 0, mpicom, ier)
    call mpi_bcast (finidat, len(finidat), mpichar, 0, mpicom, ier)
    call mpi_bcast (fsurdat, len(fsurdat), mpichar, 0, mpicom, ier)
#if (defined RTM)
    call mpi_bcast (frivinp_rtm, len(frivinp_rtm), mpichar, 0, mpicom, ier)
#endif

! surface dataset generation variables

    if (fsurdat == ' ') then
       call mpi_bcast (mksrf_fvegtyp, len(mksrf_fvegtyp), mpichar, 0, mpicom, ier)
       call mpi_bcast (mksrf_fsoitex, len(mksrf_fsoitex), mpichar, 0, mpicom, ier)
       call mpi_bcast (mksrf_fsoicol, len(mksrf_fsoicol), mpichar, 0, mpicom, ier)
       call mpi_bcast (mksrf_flanwat, len(mksrf_flanwat), mpichar, 0, mpicom, ier)
       call mpi_bcast (mksrf_furban, len(mksrf_furban), mpichar, 0, mpicom, ier)
       call mpi_bcast (mksrf_fglacier, len(mksrf_fglacier), mpichar, 0, mpicom, ier)
       call mpi_bcast (mksrf_flai, len(mksrf_flai), mpichar, 0, mpicom, ier)
    endif

! physics variables

    call mpi_bcast (conchk, 1, mpilog, 0, mpicom, ier)
    call mpi_bcast (irad, 1, mpiint, 0, mpicom, ier)
    call mpi_bcast (csm_doflxave, 1, mpilog, 0, mpicom, ier)
    call mpi_bcast (rtm_nsteps, 1, mpiint, 0, mpicom, ier)
    call mpi_bcast (wrtdia, 1, mpilog, 0, mpicom, ier)

! history and restart file variables
    
    call mpi_bcast (hist_dov2xy, size(hist_dov2xy), mpilog, 0, mpicom, ier)
    call mpi_bcast (hist_nhtfrq, size(hist_nhtfrq), mpiint, 0, mpicom, ier)
    call mpi_bcast (hist_mfilt, size(hist_mfilt), mpiint, 0, mpicom, ier)
    call mpi_bcast (hist_ndens, 1, mpiint, 0, mpicom, ier)
    call mpi_bcast (hist_chntyp, len(hist_chntyp(1,1))*size(hist_chntyp), mpichar, 0, mpicom, ier)
    call mpi_bcast (hist_fldadd, len(hist_fldadd(1))*size(hist_fldadd), mpichar, 0, mpicom, ier)
    call mpi_bcast (hist_fldaux1, len(hist_fldaux1(1))*size(hist_fldaux1), mpichar, 0, mpicom, ier)
    call mpi_bcast (hist_fldaux2, len(hist_fldaux2(1))*size(hist_fldaux2), mpichar, 0, mpicom, ier)
    call mpi_bcast (hist_crtinic, len(hist_crtinic), mpichar, 0, mpicom, ier)
    call mpi_bcast (rpntpath, len(rpntpath), mpichar, 0, mpicom, ier)

! long term archiving variables

    call mpi_bcast (mss_irt, 1, mpiint, 0, mpicom, ier)
    call mpi_bcast (mss_wpass, len(mss_wpass), mpichar, 0, mpicom, ier)
    call mpi_bcast (archive_dir, len(archive_dir), mpichar, 0, mpicom, ier)

    return
  end subroutine control_spmd
#endif

!=======================================================================

  subroutine control_print

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Write out run control variables
!
! Method: 
!
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

! ------------------------ local variables -----------------------------
    integer i  !loop index
!-----------------------------------------------------------------------

    write(6,*) 'define run:'
    write(6,*) '   run type              = ',runtyp(nsrest+1)
    write(6,*) '   case title            = ',trim(ctitle)
    write(6,*) 'input data files:'
    write(6,*) '   PFT physiology = ',trim(fpftcon)
    if (mkfsurdat) then
       write(6,*) '   generated surface dataset using raw data'
       write(6,*) '     plant types  = ',trim(mksrf_fvegtyp)
       write(6,*) '     inland water = ',trim(mksrf_flanwat)
       write(6,*) '     glacier      = ',trim(mksrf_fglacier)
       write(6,*) '     urban        = ',trim(mksrf_furban)
       write(6,*) '     soil texture = ',trim(mksrf_fsoitex)
       write(6,*) '     soil color   = ',trim(mksrf_fsoicol)
       write(6,*) '     lai and sai  = ',trim(mksrf_flai)
#if (defined OFFLINE)
       if (mksrf_offline_fgrid /= ' ') then
          write (6,*)'   land grid and mask obtained from = ',trim(mksrf_offline_fgrid)
       endif
       if (mksrf_offline_fnavyoro  /= ' ') then
          write (6,*)'   land mask obtained obtained from = ',trim(mksrf_offline_fnavyoro)
          write (6,*)'   regular grid is generated by model with'
          write (6,*)'      northern edge (degrees)  = ',mksrf_offline_edgen
          write (6,*)'      southern edge (degrees)  = ',mksrf_offline_edges
          write (6,*)'      western  edge (degrees)  = ',mksrf_offline_edgew
          write (6,*)'      eastern  edge (degrees)  = ',mksrf_offline_edgee
       endif
#endif
    else
       write(6,*) '   surface data   = ',trim(fsurdat)
    end if
    if (nsrest == 0 .and. finidat == ' ') write(6,*) '   initial data created by model'
    if (nsrest == 0 .and. finidat /= ' ') write(6,*) '   initial data   = ',trim(finidat)
    if (nsrest /= 0) write(6,*) '   restart data   = ',trim(nrevsn)
#if (defined OFFLINE)
    if (offline_atmdir /= ' ') then
       write(6,*) '   atmosperic forcing data    = ',trim(offline_atmdir) 
    end if
#elif (defined COUP_CAM)
    write(6,*) '   atmosperhic forcing data is from cam model'
#elif (defined COUP_CSM)
    write(6,*) '   atmospheric forcint data is from csm flux coupler'
#endif
#if (defined RTM)
    if (frivinp_rtm /= ' ') write(6,*) '   RTM river data       = ',trim(frivinp_rtm)
#endif
    write(6,*) 'history and restart parameters:'
    if (hist_ndens == 1) then
       write(6,*)'   history tape data will be double precision'
    else if (hist_ndens == 2) then
       write(6,*)'   history tape data will be single precision'
    end if
    write(6,101) (i,hist_dov2xy(i),i=1,nhist)
    write(6,*) '   there will be ',nhist,' history files'
    do i = 1, nhist
       if (hist_nhtfrq(i)==0) then
          write(6,*) '   history file ',i,' is monthly averaged'
       else
          write(6,*) '   history file ',i,' time interval (iterations)= ', hist_nhtfrq(i)
       endif
    end do
    write(6,104) (i,hist_mfilt(i),i=1,nhist)
    if (mss_irt /= 0) then 
       write(6,*)'   mass store path                    = ',trim(archive_dir) 
       write(6,*)'   mass store retention (days)        = ',mss_irt
       write(6,*)'   mass store write password          = ',mss_wpass
    endif
    write(6,*)'   restart pointer file directory     = ',trim(rpntdir)
    write(6,*)'   restart pointer file name          = ',trim(rpntfil)
    if (hist_crtinic == 'MONTHLY') then
       write(6,*)'initial datasets will be written monthly'
    else if (hist_crtinic == 'YEARLY') then
       write(6,*)'initial datasets will be written yearly'
! kdr added for data assimilation
    else if (hist_crtinic == 'ENDOFRUN') then
       write(6,*)'initial datasets will be written at the end of the run'
    else
       write(6,*)'initial datasets will not be produced'
    endif
    write(6,*) 'model physics parameters:'
#if (defined PERGRO)
    write(6,*) '   flag for random perturbation test is set'
#else
    write(6,*) '   flag for random perturbation test is not set'
#endif
    write(6,*) '   energy and water conservation checks   = ',conchk
    write(6,*) '   solar radiation frequency (iterations) = ',irad
#if (defined COUP_CSM)
    write(6,*) 'communication with the flux coupler'
    if (csm_doflxave) then
       write(6,*)'    data will be sent to the flux coupler ', &
            'only when an albedo calculation is performed '
       write(6,*)'     fluxes will be averaged on steps where ', &
            'communication with the flux coupler does not occur'
    else
       write(6,*)'    data will be sent and received to/from ', &
            'the flux coupler at every time step except for nstep=1'
    endif
#endif
#if (defined RTM)
    if (rtm_nsteps > 1) then
       write(6,*)'river runoff calculation performed only every ',rtm_nsteps,' nsteps'
    else
       write(6,*)'river runoff calculation performed every time step'
    endif
#endif       
    if (nsrest == 1) then
       write(6,*) 'restart warning:'
       write(6,*) '   Namelist not checked for agreement with initial run.'
       write(6,*) '   Namelist should not differ except for ending time step and run type'
    end if
    if (nsrest == 3) then
       write(6,*) 'branch warning:'
       write(6,*) '   Namelist not checked for agreement with initial run.'
       write(6,*) '   Surface data set and reference date should not differ from initial run'
    end if
#if (defined COUP_CSM)       
    write(6,*) '   last time step determined by flux coupler'
#endif
       
101 format (1x,'   history fields are grid-average    = ',4(i1,':',l1,' '))
104 format (1x,'   time samples per history file      = ',4(i1,':',i2,' '))

  end subroutine control_print

!=======================================================================

end module controlMod


