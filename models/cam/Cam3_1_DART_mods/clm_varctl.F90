! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research

#include <misc.h>
#include <preproc.h>

module clm_varctl

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varctl
!
! !DESCRIPTION:
! Module containing run control variables
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! Run control variables
!
  character(len=256) :: caseid                  ! case id
  character(len=256) :: ctitle                  ! case title
  integer :: nsrest                             ! 0: initial run. 1: restart: 3: branch
  logical, public :: brnch_retain_casename = .false. ! true => allow case name to remain the same for branch run
                                                     ! by default this is not allowed
!
! Initial file variables
!
  character(len= 8) :: hist_crtinic             ! if set to '6-HOURLY', 'MONTHLY' or 'YEARLY', write initial cond. file
! kdr                                             add 'ENDOFRUN' for DART, as in CAM
!
! Long term archive variables
!
  character(len=256) :: archive_dir             ! long term archive directory (can be mass store)
  character(len=  8) :: mss_wpass               ! mass store write password for output files
  integer            :: mss_irt                 ! mass store retention period
!
! Run input files
!
  character(len=256) :: finidat                 ! initial conditions file name
  character(len=256) :: fsurdat                 ! surface data file name
  character(len=256) :: fpftcon                 ! ASCII data file with PFT physiological constants
  character(len=256) :: nrevsn                  ! restart data file name for branch run
  character(len=256) :: frivinp_rtm             ! RTM input data file name
  character(len=256) :: offline_atmdir          ! directory for input offline model atm data forcing files (Mass Store ok)
!
! Files and logical variables for generating surface dataset
!
  logical            :: mksrf_all_pfts           ! true => surface dataset with all pft types will be generated
  real(r8)           :: mksrf_offline_edgen      ! northern edge of grid (degrees): >  -90 and < 90
  real(r8)           :: mksrf_offline_edgee      ! eastern edge of grid (degrees) : see following notes
  real(r8)           :: mksrf_offline_edges      ! southern edge of grid (degrees): >= -90 and <  90
  real(r8)           :: mksrf_offline_edgew      ! western edge of grid (degrees) : see following notes
  character(len=256) :: mksrf_offline_fgrid      ! land grid file name to use instead of generating grid
  character(len=256) :: mksrf_offline_fnavyoro   ! directory for 20 min navy orography dataset
  character(len=256) :: mksrf_fvegtyp            ! when making [fsurdat]: vegetation data file name
  character(len=256) :: mksrf_fsoitex            ! when making [fsurdat]: soil texture data file name
  character(len=256) :: mksrf_fsoicol            ! when making [fsurdat]: soil color data file name
  character(len=256) :: mksrf_flanwat            ! when making [fsurdat]: inland water data file name
  character(len=256) :: mksrf_furban             ! when making [fsurdat]: urban data file name
  character(len=256) :: mksrf_fglacier           ! when making [fsurdat]: glacier data file name
  character(len=256) :: mksrf_flai               ! when making [fsurdat]: lai data filename
!
! Physics
!
  integer :: irad         ! solar radiation frequency (iterations)
  logical :: wrtdia       ! true => write global average diagnostics to std out
  logical :: csm_doflxave ! true => only communicate with flux coupler on albedo calc time steps
!
! Rtm control variables
!
  integer :: rtm_nsteps   ! if > 1, average rtm over rtm_nsteps time steps
!
! Derived variables (run, history and restart file)
!
  character(len=256) :: rpntdir          ! directory name for local restart pointer file
  character(len=256) :: rpntfil          ! file name for local restart pointer file
  character(len=256) :: version          ! model version number
!
! Error growth perturbation limit
!
  real(r8) :: pertlim                    ! perturbation limit when doing error growth test
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein and Gordon Bonan
!
!EOP
!-----------------------------------------------------------------------

end module clm_varctl
