! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program pert_sounding

use pert_sounding_mod

implicit none

real, allocatable, dimension(:,:,:) :: z1d        ! base sounding profile to be read
real, allocatable, dimension(:,:)   :: sfc1d      ! base sounding sfc data to be read
real, allocatable, dimension(:,:)   :: p1d        ! perturbations to be added to any variable
integer                             :: nz         ! number of vertical levels

integer                             :: pert_kind  ! type of perturbation, 1:sinusoidal, 2:uncorrelated gaussian, 3:read from file
integer                             :: n_modes    ! for pert_kind=1, number of modes
integer, allocatable, dimension(:)  :: k_modes    ! for pert_kind=1, wavenumber of each mode
integer                             :: nens       ! number of ensemble members
integer                             :: npertvar   ! number of variables to be perturbed
integer, allocatable, dimension(:)  :: pertvar    ! indices of perturbation variables
real, allocatable, dimension(:)     :: magvar     ! perturbation magnitudes for perturbation variables
character(len=10)                   :: dummy
real, allocatable, dimension(:)     :: std_bot    ! lowest-level std deviation for perturbation variables
real, allocatable, dimension(:)     :: std_top    ! highest-level std deviation for perturbation variables

real, allocatable, dimension(:)     :: hgt_mid    ! transitioning from std_bot to std_top begins at this height (m)
real, allocatable, dimension(:)     :: hgt_top    ! transitioning from std_bot to std_top end at this height (m)

integer, allocatable, dimension(:)  :: lvl_mid    ! level index corresponding to hgt_mid
integer, allocatable, dimension(:)  :: lvl_top    ! level index corresponding to hgt_top

integer                             :: i, ivar, curvar, ipert, iens, ik
character(len=3)                    :: sndgvar(5)
logical                             :: this_pertvar

! Descriptions of sounding variables in the order they appear in the sounding file
sndgvar(1) = 'HGT'
sndgvar(2) = 'TH '
sndgvar(3) = 'QV '
sndgvar(4) = ' U '
sndgvar(5) = ' V '

! Read user definitions from the file "pert_sounding.input"
! --------------------------------------------------------
open(unit=10,file='pert_sounding.input',status='old')

read(10,*) pert_kind
if ((pert_kind.le.0).or.(pert_kind.gt.3)) then
   print *, "Error: pert_kind should be 1-3."
   print *, "Read:  ", pert_kind
   stop
endif
read(10,*) n_modes
if ((pert_kind.eq.1).and.(n_modes.le.0)) then
   print *, "Error: n_modes should be positive."
   print *, "Read:  ", n_modes
   stop
endif
allocate(k_modes(n_modes))
read(10,*) (k_modes(i),i=1,n_modes)
k_modes(:) = abs(k_modes)
read(10,*) nens
if (nens.le.0) then
   print *, "Error: nens should be positive."
   print *, "Read:  ", nens
   stop
endif
read(10,*) npertvar
if ((npertvar.le.0) .or. (npertvar.ge.5)) then
   print *, "Error: npertvar should be 1-4."
   print *, "Read:  ", npertvar
   stop
endif
allocate(pertvar(npertvar))
allocate(magvar(npertvar))
allocate(std_bot(npertvar))
allocate(std_top(npertvar))
allocate(hgt_mid(npertvar))
allocate(hgt_top(npertvar))
allocate(lvl_mid(npertvar))
allocate(lvl_top(npertvar))
std_bot = 0.0
std_top = 0.0
hgt_mid = 0.0
hgt_top = 0.0
lvl_mid = 0
lvl_top = 0
read(10,*) (pertvar(i),i=1,npertvar)
if (any(pertvar.le.0) .or. any(pertvar.ge.5)) then
   print *, "Error: Some or all pertvar values outside 1-4 range."
   print *, ""
   stop
endif
read(10,*) (magvar(i),i=1,npertvar)
magvar(:) = abs(magvar(:))
read(10,*) dummy
read(10,*) (std_bot(i),i=1,npertvar)
std_bot(:) = abs(std_bot(:))
read(10,*) (std_top(i),i=1,npertvar)
std_top(:) = abs(std_top(:))
read(10,*) (hgt_mid(i),i=1,npertvar)
hgt_mid(:) = abs(hgt_mid(:))
read(10,*) (hgt_top(i),i=1,npertvar)
hgt_top(:) = abs(hgt_top(:))
do i = 1, npertvar
   if (hgt_top(i).lt.hgt_mid(i)) then
      print *, "Variable: ", sndgvar(pertvar(i)+1) 
      print *, "Error: hgt_top must be equal to or greater than hgt_mid"
      print *, ""
      stop
   endif
enddo

close(unit=10)

! --------------------------------------------------------
! Determine number of height levels in base sounding file
call read_numlevels(nz)
if (nz.eq.0) then
   print *, "Error: No levels found in base sounding file."
   print *, ""
   stop
endif

! --------------------------------------------------------
! Based on number of vert. levels and ensemble size, allocate variables
allocate(sfc1d(nens+1,4))
allocate(z1d(nens+1,5,nz))
allocate(p1d(nens,nz))
sfc1d = 0.0
z1d   = 0.0
p1d   = 0.0

! --------------------------------------------------------
! Read the base sounding
print *, 'Reading base sounding'
call read_sounding(nz,sfc1d(1,:),z1d(1,:,:))

! --------------------------------------------------------
! Determine levels where top standard deviations are defined
! This is only necessary if standard deviation bounds (at any level)
! are requested
lvl_mid = 0
lvl_top = 0
do ivar = 1, npertvar
   if ((std_bot(ivar).ne.0).or.(std_top(ivar).ne.0)) then
      call level1d(nz,z1d(1,1,:),hgt_mid(ivar),lvl_mid(ivar))
      call level1d(nz,z1d(1,1,:),hgt_top(ivar),lvl_top(ivar))
   endif
enddo

! --------------------------------------------------------
! Initialize random number seed
call rand_init()

! --------------------------------------------------------
! Loop through all variables and add perturbations
! (zero perturbations added to unperturbed variables)
ipert = 0
do ivar = 1, 5
   p1d = 0
   curvar = ivar-1  ! index of requested variable (skipping height)
   this_pertvar = .false.
   do i = 1, npertvar
      if (pertvar(i).eq.curvar) this_pertvar = .true.
   enddo
   if (this_pertvar) then
      ipert = ipert + 1
      print *, 'Adding perturbations to variable ',sndgvar(ivar)
      ! Obtain perturbations for each member
      do iens = 1, nens
         call pert1d(pert_kind,p1d(iens,:),1.0,z1d(1,1,:),nz,k_modes,n_modes)
      enddo
      ! Normalize spread to specified standard deviation - if either bound is non-zero
      if ((std_bot(ipert).ne.0).or.(std_top(ipert).ne.0)) then
         call norm1d(p1d,std_bot(ipert),std_top(ipert),lvl_mid(ipert),lvl_top(ipert),z1d(1,1,:),nz,nens)
      endif
   else
      print *, 'No perturbations added to variable ', sndgvar(ivar)
   endif
   ! Add perturbations to base state for each ensemble member
   do iens = 1, nens
      ! Check if this is applicable to surface variables
      if ((ivar.eq.1).or.(ivar.eq.2).or.(ivar.eq.3)) then
         sfc1d(1+iens,ivar) = sfc1d(1,ivar) + p1d(iens,1)
      endif
      ! Now add perturbations at each level
      do ik = 1, nz
         z1d(1+iens,ivar,ik) = z1d(1,ivar,ik) + p1d(iens,ik)
      enddo
   enddo
enddo
sfc1d(2:1+iens,4) = sfc1d(1,4)

open(unit=9,file='test_sounding.output',status='replace')
do ivar = 2,5
   do ik = 1, nz
      do iens = 1, nens+1
         write(9,*) z1d(iens,ivar,ik)
      enddo
   enddo
enddo
close(unit=9)

! Write to ensemble sounding files
do iens = 1, nens
   call write_sounding_ens(nz,sfc1d(1+iens,:),z1d(1+iens,:,:),iens)
enddo


end program pert_sounding

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
