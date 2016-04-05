! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research

#include <misc.h>
#include <preproc.h>
#if ( defined SCAM )
#include <max.h>
#endif

module controlMod

! <next few lines under version control, do not edit>
! $URL: http://subversion.ucar.edu/DAReS/DART/trunk/models/cam/model_mod.f90 $
! $Id: model_mod.f90 2721 2007-03-27 00:08:01Z thoar $
! $Revision: 2721 $
! $Date: 2007-03-26 18:08:01 -0600 (Mon, 26 Mar 2007) $

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: controlMod
!
! !DESCRIPTION:
! Module which initializes run control variables. The following possible
! namelist variables are set default values and possibly read in on startup
!
! === define run =======================
!
!    o caseid     = 256 character case name
!    o ctitle     = 256 character case title
!    o nsrest     = integer flag. 0: initial run. 1: restart: 3: branch
!
! === model time =======================
!
!    o dtime      = integer model time step (s)
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
!    o hist_fincl1   = 10 character name of fields for first  auxillary history file
!    o hist_fincl2   = 10 character name of fields for second auxillary history file
!    o hist_fincl3   = 10 character name of fields for first  auxillary history file
!    o hist_fincl4   = 10 character name of fields for second auxillary history file
!    o hist_fincl5   = 10 character name of fields for first  auxillary history file
!    o hist_fincl6   = 10 character name of fields for second auxillary history file
!    o hist_fexcl1   = 8  character name of fields for first  auxillary history file
!    o hist_fexcl2   = 8  character name of fields for second auxillary history file
!    o hist_fexcl3   = 8  character name of fields for first  auxillary history file
!    o hist_fexcl4   = 8  character name of fields for second auxillary history file
!    o hist_fexcl5   = 8  character name of fields for first  auxillary history file
!    o hist_fexcl6   = 8  character name of fields for second auxillary history file
!    o hist_crtinic  = 8  character frequency to generate initial dataset
!                         ['6-HOURLY','DAILY','MONTHLY','YEARLY','NONE']
! kdr                     added 'ENDOFRUN' option for DART as in CAM
!    o rpntpath      = 256 character full UNIX pathname of the local restart pointer file.
!                      This file must exist when the model is restarted.
!                      This file is overwritten every time new restart data files are output.
!
! === long term archiving =====
!
!    o archive_dir = 256 character long term archive directory (can be MSS directory)
!    o mss_irt     = integer mass store retention period (days)
!    o mss_wpass   = 8 character mass store write password for output data sets
!
! === model physics ===
!
!    o irad         = integer solar radiation frequency (+ = iteration. - = hour)
!    o wrtdia       = true if want output written
!    o csm_doflxave = true => flux averaging is to be performed (only used for csm mode)
!
! === rtm control variables ===
!
!    o rtm_nsteps  = if > 1, average rtm over rtm_nsteps time steps
!
! When coupled to CAM: base calendar info, nstep, nestep, nsrest, and time
! step are input to the land model from CAM. The values in the clmexp namelist
! are not used. The minimum namelist parameters are:
! o fsurdat
! o finidat
! o fpftcon
! When running in offline mode, the minimum namelist parameters are:
! o nsrest
! o nestep or nelapse
! o fsurdat
! o finidat
! o dtime
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varctl
  use spmdMod
  use shr_sys_mod, only : shr_sys_getenv
  use decompMod  , only : clump_pproc
  use histFileMod, only : max_tapes, max_namlen, &
                          hist_empty_htapes, hist_dov2xy, &
                          hist_avgflag_pertape, hist_type1d_pertape, &
                          hist_nhtfrq, hist_ndens, hist_mfilt, &
                          hist_fincl1, hist_fincl2, hist_fincl3, &
                          hist_fincl4, hist_fincl5, hist_fincl6, &
                          hist_fexcl1, hist_fexcl2, hist_fexcl3, &
                          hist_fexcl4, hist_fexcl5, hist_fexcl6
  use shr_const_mod, only : SHR_CONST_CDAY
  use abortutils, only : endrun
#if (defined COUP_CSM)
  use clm_csmMod   , only : csm_dtime
#endif
#if ( defined SCAM )
  use scamMod, only :lsmsurffile,lsminifile
#endif
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: control_init  ! initial run control information
  public :: control_print ! print run control information
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! PRIVATE TYPES:
! Namelist variables only used locally
  logical            :: mkfsurdat            ! true => make surface data from raw data
  character(len=256) :: rpntpath             ! full UNIX pathname of restart pointer file
  character(len=  7) :: runtyp(4)            ! run type
#if (defined _OPENMP)
   integer, external :: omp_get_max_threads  ! max number of threads that can execute
                                             ! concurrently in a single parallel region
#endif
!-----------------------------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: control_init
!
! !INTERFACE:
  subroutine control_init (cam_caseid , cam_ctitle, cam_irad , cam_nsrest, &
                           cam_crtinic, cam_nhtfrq, cam_mfilt, cam_irt )
!
! !DESCRIPTION:
! Initialize CLM run control information
!
! !USES:
#if (defined OFFLINE) || (defined COUP_CSM)
  use time_manager, only : calendar, dtime, nestep, nelapse, start_ymd, &
                           start_tod, stop_ymd, stop_tod, ref_ymd, ref_tod
#else
  use time_manager, only : get_step_size, is_perpetual
#endif
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'

    character(len=*), optional, intent(in) :: cam_caseid    ! cam caseid
    character(len=*), optional, intent(in) :: cam_ctitle    ! cam title
    integer         , optional, intent(in) :: cam_irad      ! cam radiation frequency
    integer         , optional, intent(in) :: cam_nsrest    ! cam run type
    character(len=*), optional, intent(in) :: cam_crtinic   ! cam initial dataset frequency
    integer         , optional, intent(in) :: cam_nhtfrq    ! cam history write freq for tape 1
    integer         , optional, intent(in) :: cam_mfilt     ! cam number of files per tape for tape 1
    integer         , optional, intent(in) :: cam_irt       ! cam mss retention time
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: homedir   ! full UNIX filepath name of home directory
    character(len=256) :: logid     ! logid part of file path name
    character(len=256) :: cap       ! upper case logid
    character(len=  1) :: ctmp      ! character temporary
    integer :: i,j,n                ! loop indices
    integer :: iundef               ! integer undefined value
    real(r8):: rundef               ! real undefined value
    integer :: ierr                 ! error code
!------------------------------------------------------------------------

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

    ! ----------------------------------------------------------------------
    ! Namelist Variables
    ! ----------------------------------------------------------------------

    namelist /clmexp/  &
         ctitle, caseid, nsrest,  &
         calendar, dtime, nelapse, nestep, start_ymd, start_tod,  &
         stop_ymd, stop_tod, ref_ymd, ref_tod, &
         hist_empty_htapes, hist_dov2xy, &
         hist_avgflag_pertape, hist_type1d_pertape, &
         hist_nhtfrq, hist_ndens, hist_mfilt, &
         hist_fincl1, hist_fincl2, hist_fincl3, &
         hist_fincl4, hist_fincl5, hist_fincl6, &
         hist_fexcl1, hist_fexcl2, hist_fexcl3, &
         hist_fexcl4, hist_fexcl5, hist_fexcl6, &
         hist_crtinic, archive_dir, mss_wpass, mss_irt, &
         nrevsn, rpntpath, offline_atmdir, &
         finidat, fsurdat, fpftcon, frivinp_rtm, &
         mksrf_all_pfts, mksrf_fvegtyp, mksrf_fsoitex, mksrf_fsoicol, mksrf_flanwat, &
         mksrf_fglacier, mksrf_furban, mksrf_flai, mksrf_offline_fgrid, &
         mksrf_offline_edgen, mksrf_offline_edgee, mksrf_offline_edges, &
         mksrf_offline_edgew, mksrf_offline_fnavyoro, &
         irad, wrtdia, csm_doflxave, rtm_nsteps, pertlim, &
         clump_pproc, brnch_retain_casename

    ! ----------------------------------------------------------------------
    ! Default values
    ! ----------------------------------------------------------------------

    if (masterproc) then
       write(6,*) 'Attempting to initialize run control settings .....'
    endif

    runtyp(0 + 1) = 'initial'
    runtyp(1 + 1) = 'restart'
    runtyp(3 + 1) = 'branch '

    iundef = -9999999
    rundef = -9999999.

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

    mksrf_all_pfts         = .false.
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

    hist_crtinic = 'YEARLY'
    rpntpath = 'not_specified'

    ! other namelist variables

    irad = -1
    wrtdia = .false.
    csm_doflxave = .true.
    pertlim = 0.

#if (defined RTM)
    ! If rtm_nsteps is not set in the namelist then
    ! will be given default value below

    rtm_nsteps = -999
#endif

    ! Set clumps per procoessor

    clump_pproc = 1
#if (defined _OPENMP)
    clump_pproc = omp_get_max_threads()
#else
#if (defined UNICOSMP)
#if (defined SSP)
    clump_pproc = 1
#else
    clump_pproc = 4
#endif
#endif
#endif

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input. Override if coupled to CAM
    ! ----------------------------------------------------------------------

    if (masterproc) then
#if ( defined SCAM )
       fsurdat=lsmsurffile
       finidat=lsminifile
#else
       read(5, clmexp, iostat=ierr)
       if (ierr /= 0) then
          if (masterproc) then
             write(6,*)'error: namelist input resulted in error code ',ierr
          endif
          call endrun
       endif
#endif

#if (defined COUP_CAM)
       ! Override select set of namelist values with CAM input

       caseid = cam_caseid
       ctitle = cam_ctitle
       irad   = cam_irad
       nsrest = cam_nsrest
       hist_crtinic   = cam_crtinic
       hist_mfilt(1)  = cam_mfilt
       hist_nhtfrq(1) = cam_nhtfrq
#if (defined PERGRO)
       hist_empty_htapes = .true.
       hist_fincl1 = 'TSA'
#endif
       mss_irt = cam_irt
       if (is_perpetual()) then
#if (defined RTM) || (defined DGVM)
          write(6,*)'RTM or DGVM cannot be defined in perpetual mode'
          call endrun()
#endif
          if (finidat == ' ') then
             write(6,*)'must specify initial dataset for perpetual mode'
             call endrun()
          end if
    end if
#endif

#if (defined OFFLINE)
       ! Consistency checks

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

       ! If archive directory not input in namelist - set default from caseid

       if (archive_dir == ' ') then
          logid  = ' '
          call shr_sys_getenv('LOGNAME', logid, ierr)
          if (ierr /= 0) then
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

#if (defined RTM)
       ! If rtm_nsteps was not entered in the namelist, give it the
       ! following default value

       if (rtm_nsteps == -999) then
#if (defined COUP_CAM)
          rtm_nsteps = (3600*3)/get_step_size() ! 3 hours
#else
          rtm_nsteps = (3600*3)/dtime           ! 3 hours
#endif
       endif
#endif

       ! Check that hist_type_1d is not set for primary tape

       if (hist_type1d_pertape(1) /= ' ') then
          write(6,*)'CONTROL_INIT error: hist_type1d_pertape can only be set for tapes 2-6'
          call endrun()
       end if

    end if   ! end of if-masterproc block

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

    do i = 1, max_tapes
       if (hist_nhtfrq(i) == 0) then
          hist_mfilt(i) = 1
       else if (hist_nhtfrq(i) < 0) then
          hist_nhtfrq(i) = nint(-hist_nhtfrq(i)*SHR_CONST_CDAY/(24.*dtime))
       endif
    end do

    if (rpntpath == 'not_specified') then
       call shr_sys_getenv('HOME', homedir, ierr)
       rpntpath = trim(homedir)//'/lnd.'//trim(caseid)//'.rpointer'
    endif

    if (nsrest == 0) nrevsn = ' '
    if (nsrest == 1) nrevsn = 'set by restart pointer file file'
    if (nsrest == 3 .and. nrevsn == ' ') then
       if (masterproc) write(6,*) 'error: need to set restart data file name'
       call endrun
    end if

! kdr added ENDOFRUN
    if (trim(hist_crtinic) /= 'MONTHLY'  .and. trim(hist_crtinic) /= 'YEARLY' .and. &
        trim(hist_crtinic) /= '6-HOURLY' .and. trim(hist_crtinic) /= 'DAILY'  .and. &
        trim(hist_crtinic) /= 'ENDOFRUN') then
       hist_crtinic = 'NONE'
    endif
#if (defined DGVM)
    !only permit yearly initial datasets for DGVM
    if (trim(hist_crtinic) == 'MONTHLY') hist_crtinic = 'YEARLY'
#endif

    ! ----------------------------------------------------------------------
    ! Restart pointer file
    ! ----------------------------------------------------------------------

    ! split the full pathname of the restart pointer file into a
    ! directory name and a file name
    ! check if the directory exists and if not, make it

    rpntdir = ' '
    rpntfil = ' '
    do n = len_trim(rpntpath),1,-1
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

  end subroutine control_init

#if (defined SPMD)

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: control_spmd
!
! !INTERFACE:
  subroutine control_spmd()
!
! !DESCRIPTION:
! Distribute namelist data all processors. The cpp SPMD definition
! provides for the funnelling of all program i/o through the master
! processor. Processor 0 either reads restart/history data from the
! disk and distributes it to all processors, or collects data from
! all processors and writes it to disk.
!
! !USES:
!
#if (defined OFFLINE) || (defined COUP_CSM)
    use time_manager, only : calendar, dtime, nestep, nelapse, start_ymd, &
         start_tod, stop_ymd, stop_tod, ref_ymd, ref_tod
#endif
    use spmdMod, only : mpicom
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer ier       !error code
!-----------------------------------------------------------------------

    ! run control variables

    call mpi_bcast (caseid, len(caseid), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (ctitle, len(ctitle), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (nsrest,           1, MPI_INTEGER  , 0, mpicom, ier)

#if (defined OFFLINE) || (defined COUP_CSM)
    call mpi_bcast (nestep   , 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (nelapse  , 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (dtime    , 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (start_ymd, 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (start_tod, 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (stop_ymd , 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (stop_tod , 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (ref_ymd  , 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (ref_tod  , 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (calendar ,len(calendar), MPI_CHARACTER, 0, mpicom, ier)
#endif

    ! initial file variables

    call mpi_bcast (nrevsn, len(nrevsn), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (finidat, len(finidat), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fsurdat, len(fsurdat), MPI_CHARACTER, 0, mpicom, ier)
#if (defined RTM)
    call mpi_bcast (frivinp_rtm, len(frivinp_rtm), MPI_CHARACTER, 0, mpicom, ier)
#endif

    ! surface dataset generation variables

    if (fsurdat == ' ') then
       call mpi_bcast (mksrf_all_pfts, 1                  , MPI_LOGICAL  , 0, mpicom, ier)
       call mpi_bcast (mksrf_fvegtyp , len(mksrf_fvegtyp) , MPI_CHARACTER, 0, mpicom, ier)
       call mpi_bcast (mksrf_fsoitex , len(mksrf_fsoitex) , MPI_CHARACTER, 0, mpicom, ier)
       call mpi_bcast (mksrf_fsoicol , len(mksrf_fsoicol) , MPI_CHARACTER, 0, mpicom, ier)
       call mpi_bcast (mksrf_flanwat , len(mksrf_flanwat) , MPI_CHARACTER, 0, mpicom, ier)
       call mpi_bcast (mksrf_furban  , len(mksrf_furban)  , MPI_CHARACTER, 0, mpicom, ier)
       call mpi_bcast (mksrf_fglacier, len(mksrf_fglacier), MPI_CHARACTER, 0, mpicom, ier)
       call mpi_bcast (mksrf_flai    , len(mksrf_flai)    , MPI_CHARACTER, 0, mpicom, ier)
    endif

    ! physics variables

    call mpi_bcast (irad, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (csm_doflxave, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (rtm_nsteps, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (wrtdia, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! history file variables

    call mpi_bcast (hist_empty_htapes, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (hist_dov2xy, size(hist_dov2xy), MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (hist_nhtfrq, size(hist_nhtfrq), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (hist_mfilt, size(hist_mfilt), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (hist_ndens, size(hist_ndens), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (hist_crtinic, len(hist_crtinic), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_avgflag_pertape, size(hist_avgflag_pertape), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_type1d_pertape, max_namlen*size(hist_type1d_pertape), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl1, max_namlen*size(hist_fexcl1), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl2, max_namlen*size(hist_fexcl2), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl3, max_namlen*size(hist_fexcl3), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl4, max_namlen*size(hist_fexcl4), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl5, max_namlen*size(hist_fexcl5), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl6, max_namlen*size(hist_fexcl6), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl1, (max_namlen+2)*size(hist_fincl1), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl2, (max_namlen+2)*size(hist_fincl2), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl3, (max_namlen+2)*size(hist_fincl3), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl4, (max_namlen+2)*size(hist_fincl4), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl5, (max_namlen+2)*size(hist_fincl5), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl6, (max_namlen+2)*size(hist_fincl6), MPI_CHARACTER, 0, mpicom, ier)

    ! restart file variables

    call mpi_bcast (rpntpath, len(rpntpath), MPI_CHARACTER, 0, mpicom, ier)

    ! clump decomposition variables

    call mpi_bcast (clump_pproc, 1, MPI_INTEGER, 0, mpicom, ier)

    ! long term archiving variables

    call mpi_bcast (mss_irt, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (mss_wpass, len(mss_wpass), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (archive_dir, len(archive_dir), MPI_CHARACTER, 0, mpicom, ier)

    ! error growth perturbation limit
    call mpi_bcast (pertlim, 1, MPI_REAL8, 0, mpicom, ier)

  end subroutine control_spmd
#endif

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: control_print
!
! !INTERFACE:
  subroutine control_print ()
!
! !DESCRIPTION:
! Write out run control variables
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer i  !loop index
!------------------------------------------------------------------------

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
    if (mss_irt /= 0) then
       write(6,*) 'Mass store control values'
       write(6,*)'   mass store path                    = ',trim(archive_dir)
       write(6,*)'   mass store retention (days)        = ',mss_irt
       write(6,*)'   mass store write password          = ',mss_wpass
    endif
    write(6,*) 'Restart parameters:'
    write(6,*)'   restart pointer file directory     = ',trim(rpntdir)
    write(6,*)'   restart pointer file name          = ',trim(rpntfil)
    if (hist_crtinic == 'MONTHLY') then
       write(6,*)'initial datasets will be written monthly'
    else if (hist_crtinic == 'YEARLY') then
       write(6,*)'initial datasets will be written yearly'
    else if (hist_crtinic == 'DAILY') then
       write(6,*)'initial datasets will be written daily'
    else if (hist_crtinic == '6-HOURLY') then
       write(6,*)'initial datasets will be written 6-hourly'
    else
       write(6,*)'initial datasets will not be produced'
    endif
    write(6,*) 'model physics parameters:'
#if (defined PERGRO)
    write(6,*) '   flag for random perturbation test is set'
#else
    write(6,*) '   flag for random perturbation test is not set'
#endif
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
#if (defined PERGRO)
    write(6,*) '   perturbation limit = ',pertlim
#endif

  end subroutine control_print

end module controlMod
