! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program nmlbld_rose

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$
! $Name$

   use     types_mod, only: r8, pi
   use utilities_mod, only: open_file, close_file,  &
                            error_handler, E_ERR, E_MSG, E_WARN, logfileunit, &
                            initialize_utilities, register_module, &
                            find_namelist_in_file, check_namelist_read

   implicit none

   ! CVS Generated file description for error handling, do not edit
   character(len=128) :: &
   source   = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

   logical            :: output_prog_diag = .false.
   character (len=50) :: input_dir = '../DAinput/'
   character (len=50) :: out_dir   = '../DAoutput/'
   character (len=30) :: ncep_file = 'nmc_lbc.02.nc'
   character (len=30) :: restart_file = 'NMC_SOC.day151_2002.dat'  
   real(kind=r8)      :: h_tune = pi
   real(kind=r8)      :: z_tune = 1.0
   real(kind=r8)      :: target_time =  0.125 ! e.g, 168 = 7 days * 24 [hr]  
   integer            :: ntime = 8

   namelist /rose_nml/ target_time, &
                       input_dir, out_dir, ncep_file, restart_file,&
                       output_prog_diag, &
                       h_tune, z_tune, &
                       ntime

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

   iunit = open_file('rose.nml',action = 'write')

   write(iunit, '("&rose_nml")') 
   write(iunit, '("ntime = ", i2)')             ntime
   write(iunit, '("output_prog_diag = ", l)')   output_prog_diag
   write(iunit, '("input_dir = ", 3a)')         "'",trim(input_dir),"'"
   write(iunit, '("out_dir = ", 3a)')           "'",trim(out_dir),"'"
   write(iunit, '("ncep_file = ", 3a)')         "'",trim(ncep_file),"'"
   write(iunit, '("restart_file = ", 3a)')      "'",trim(restart_file),"'"
   write(iunit, '("h_tune = ", f15.10)')        h_tune
   write(iunit, '("z_tune = ", f15.10)')        z_tune
   write(iunit, '("target_time = ", f15.10,"/")')   target_time

end program nmlbld_rose
