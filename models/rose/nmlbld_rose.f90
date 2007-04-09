! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program nmlbld_rose

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

   use     types_mod, only: r8, pi
   use utilities_mod, only: open_file, close_file,  &
                            error_handler, E_ERR, E_MSG, E_WARN, logfileunit, &
                            initialize_utilities, register_module, &
                            find_namelist_in_file, check_namelist_read

   implicit none

   ! version controlled file description for error handling, do not edit
   character(len=128), parameter :: &
      source   = "$URL$", &
      revision = "$Revision$", &
      revdate  = "$Date$"

   logical            :: output_prog_diag = .false.
   character (len=50) :: input_dir = '../input_current/'
   character (len=50) :: out_dir   = '/ptmp/tmatsuo/rose/'
   character (len=30) :: restart_file = 'dyn_restart_001-1999.nc'  
   real(kind=r8)      :: amp_tune = 1.
   real(kind=r8)      :: pha_tune = 0.
   real(kind=r8)      :: target_time = 0.125 ![hr] 
   integer            :: ntime = 8
   integer            :: ens_element = 1

   namelist /rose_nml/ target_time, &
                       input_dir, out_dir, restart_file,&
                       output_prog_diag, &
                       amp_tune, pha_tune, &
                       ntime, ens_element

   integer :: iunit, io
   character(len=129) :: err_string, nml_string

   call initialize_utilities('nmlbld_rose')
   call register_module(source,revision,revdate)

   ! Read the namelist entry
   call find_namelist_in_file("rose.nml_default", "rose_nml", iunit)
   read(iunit, nml = rose_nml, iostat = io)
   call check_namelist_read(iunit, io, "rose_nml")

   ! Read another piece of information and add to namelist.

   read(*,*)  target_time
   read(*,*)  ens_element

   iunit = open_file('rose.nml',action = 'write')

   write(iunit, '("&rose_nml")') 
   write(iunit, '("ntime = ", i2)')             ntime
   write(iunit, '("output_prog_diag = ", l)')   output_prog_diag
   write(iunit, '("input_dir = ", 3a)')         "'",trim(input_dir),"'"
   write(iunit, '("out_dir = ", 3a)')           "'",trim(out_dir),"'"
   write(iunit, '("restart_file = ", 3a)')      "'",trim(restart_file),"'"
   write(iunit, '("amp_tune = ", f15.10)')      amp_tune
   write(iunit, '("pha_tune = ", f15.10)')      pha_tune
   write(iunit, '("ens_element = ", i3)')       ens_element
   write(iunit, '("target_time = ", f15.10,"/")')   target_time
   

end program nmlbld_rose
