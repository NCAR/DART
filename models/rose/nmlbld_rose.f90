! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program nmlbld_rose

! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$

   use types_mod, only: r8, pi
   use utilities_mod, only: open_file, close_file, &
                            check_nml_error, file_exist
   implicit none
   logical            :: old_restart = .true.
   logical            :: output_prog_diag = .false.
   character (len=50) :: input_dir = '../DAinput/'
   character (len=50) :: out_dir   = '../DAoutput/'
   character (len=30) :: ncep_file = 'nmc_lbc.02.nc'
   character (len=30) :: restart_file = 'NMC_SOC.day151_2002.dat'  
   real(kind=r8)      :: h_tune = pi
   real(kind=r8)      :: z_tune = 1.0
   real(kind=r8)      :: target_time =  0.125 ! e.g, 168 = 7 days * 24 [hr]  
   integer            :: nstart = 0
   integer            :: ntime = 8

   namelist /rose_nml/ old_restart, nstart, target_time, &
                       input_dir, out_dir, ncep_file, restart_file,&
                       output_prog_diag, &
                       h_tune, z_tune, &
                       ntime

   integer :: iunit, io, ierr

   if(file_exist('rose.nml_default')) then
      iunit = open_file('rose.nml_default', action = 'read')
      ierr = 1
      do while(ierr /= 0)
         read(iunit, nml = rose_nml, iostat = io)
         ierr = check_nml_error(io, 'rose_nml')
      end do  
         call close_file(iunit)
   endif

      read(*,*)  old_restart, nstart, target_time

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
      write(iunit, '("old_restart = ", l)')        old_restart
      write(iunit, '("nstart = ", i2)')            nstart
      write(iunit, '("target_time = ", f15.10,"/")')   target_time

end program nmlbld_rose
