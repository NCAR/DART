! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: dart_to_cice.f90 8565 2015-09-11 17:16:08Z hkershaw $

program dart_to_cice

!----------------------------------------------------------------------
! purpose: muck with cice state vector after filter
!
! method: Read in restart (restart with prior) and out restart (restart 
!         with posterior) written by DART after filter. 
!
! author: C Bitz June 2016
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, nc_check, E_ERR, file_exist, &
                             error_handler
use netcdf 

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/rma_cice/models/CICE/dart_to_cice.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 8565 $"
character(len=128), parameter :: revdate  = "$Date: 2015-09-11 10:16:08 -0700 (Fri, 11 Sep 2015) $"

!------------------------------------------------------------------

character (len = 128) :: dart_to_cice_input_file   = 'dart_restart.nc'
character (len = 128) :: balance_method = 'simple_squeeze'
namelist /dart_to_cice_nml/ dart_to_cice_input_file, balance_method

character(len=512) :: msgstring
character (len = 15) :: varname
integer, parameter :: Nx=320, Ny=384  ! hardwire for now
integer, parameter :: Ncat=5          ! number of categories in ice-thickness dist, hardwire
real(r8) :: aicen(Nx,Ny,Ncat), vicen(Nx,Ny,Ncat), vsnon(Nx,Ny,Ncat)
real(r8) :: aice(Nx,Ny)
integer  :: i, j, k
integer :: VarID, ncid, iunit, io
real(r8) :: squeeze
logical  :: update_restart = .false.

!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_cice')

call find_namelist_in_file("input.nml", "dart_to_cice_nml", iunit)
read(iunit, nml = dart_to_cice_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_cice_nml")

write(*,*)
write(*,'(''dart_to_cice:converting DART output restart file '',A, &
      &'' to one CICE will like'')') &
     trim(dart_to_cice_input_file)

if ( .not. file_exist(dart_to_cice_input_file) ) then
   write(msgstring,*) 'cannot open file ', trim(dart_to_cice_input_file),' for updating.'
   call error_handler(E_ERR,'dart_to_cice:filename not found ',trim(dart_to_cice_input_file))
endif

! open file with read and write 
call nc_check( nf90_open(trim(dart_to_cice_input_file), NF90_WRITE, ncid), &
                  'dart_to_cice', 'open '//trim(dart_to_cice_input_file))

! get the key restart variables

varname='aicen'
call nc_check(NF90_inq_varid(ncid, trim(varname), VarID), &
         'dart_to_cice', trim(varname)//' inq_varid '//trim(dart_to_cice_input_file))
call nc_check(nf90_get_var(ncid, VarID, aicen), 'dart_to_cice', &
         'get_var '//trim(varname))

varname='vicen'
call nc_check(NF90_inq_varid(ncid, trim(varname), VarID), &
         'dart_to_cice', trim(varname)//' inq_varid '//trim(dart_to_cice_input_file))
call nc_check(nf90_get_var(ncid, VarID, vicen), 'dart_to_cice', &
         'get_var '//trim(varname))

varname='vsnon'
call nc_check(NF90_inq_varid(ncid, trim(varname), VarID), &
         'dart_to_cice', trim(varname)//' inq_varid '//trim(dart_to_cice_input_file))
call nc_check(nf90_get_var(ncid, VarID, vsnon), 'dart_to_cice', &
         'get_var '//trim(varname))

aice = aicen(:,:,1)
do k = 2, Ncat  
  aice = aice+aicen(:,:,k)
enddo

! fix me if (.not. present(balance_method)) balance_method='simple_squeeze'

SELECT CASE (balance_method)
  CASE ('simple_squeeze')
     do j = 1, Ny   ! size(data_3d_array,2)
        do i = 1, Nx   ! size(data_3d_array,1)
           if (aice(i,j).lt.0.) then
              aicen(i,j,:)=0.
              vicen(i,j,:)=0.
              vsnon(i,j,:)=0.
              update_restart = .true.
           elseif (aice(i,j).gt.1.) then
              squeeze=1./aice(i,j)
              aicen(i,j,:)=aicen(i,j,:)*squeeze
              vicen(i,j,:)=vicen(i,j,:)*squeeze
              vsnon(i,j,:)=vsnon(i,j,:)*squeeze
              update_restart = .true.
           endif
        enddo
     enddo
  CASE ('tendency_weight')  ! this is identical to the above for now
     do j = 1, Ny   ! size(data_3d_array,2)
        do i = 1, Nx   ! size(data_3d_array,1)
           if (aice(i,j).lt.0.) then
              aicen(i,j,:)=0.
              vicen(i,j,:)=0.
              vsnon(i,j,:)=0.
              update_restart = .true.
           elseif (aice(i,j).gt.1.) then
              squeeze=1./aice(i,j)
              aicen(i,j,:)=aicen(i,j,:)*squeeze
              vicen(i,j,:)=vicen(i,j,:)*squeeze
              vsnon(i,j,:)=vsnon(i,j,:)*squeeze
              update_restart = .true.
           endif
        enddo
     enddo
END SELECT

!for testing make something to fix
!  update_restart = .true.
!  aicen(10,10,1)=1.1
!  write(*,*) (aicen(10,10,k), k=1,5)

if (update_restart) then
! now update the variables

varname='aicen'
call nc_check(NF90_inq_varid(ncid, trim(varname), VarID), &
         'dart_to_cice', trim(varname)//' inq_varid '//trim(dart_to_cice_input_file))
call nc_check(nf90_put_var(ncid, VarID, aicen), 'dart_to_cice', &
         'get_var '//trim(varname))

varname='vicen'
call nc_check(NF90_inq_varid(ncid, trim(varname), VarID), &
         'dart_to_cice', trim(varname)//' inq_varid '//trim(dart_to_cice_input_file))
call nc_check(nf90_put_var(ncid, VarID, vicen), 'dart_to_cice', &
         'get_var '//trim(varname))

varname='vsnon'
call nc_check(NF90_inq_varid(ncid, trim(varname), VarID), &
         'dart_to_cice', trim(varname)//' inq_varid '//trim(dart_to_cice_input_file))
call nc_check(nf90_put_var(ncid, VarID, vsnon), 'dart_to_cice', &
         'get_var '//trim(varname))

endif

call finalize_utilities('dart_to_cice')

end program dart_to_cice

! <next few lines under version control, do not edit>
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/rma_cice/models/CICE/dart_to_cice.f90 $
! $Id: dart_to_cice.f90 8565 2015-09-11 17:16:08Z hkershaw $
! $Revision: 8565 $
! $Date: 2015-09-11 10:16:08 -0700 (Fri, 11 Sep 2015) $
