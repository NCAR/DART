! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

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
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------

character(len=256) :: dart_to_cice_input_file   = 'dart_restart.nc'
character(len=256) :: balance_method = 'simple_squeeze'
namelist /dart_to_cice_nml/ dart_to_cice_input_file, balance_method

character(len=512) :: string1, string2, string3, msgstring
character(len=15) :: varname

integer :: Nx, Ny
integer :: Ncat   ! number of categories in ice-thickness dist

real(r8), allocatable :: aicen(:,:,:)
real(r8), allocatable :: vicen(:,:,:)
real(r8), allocatable :: vsnon(:,:,:)
real(r8), allocatable ::  aice(:,:)

! real(r8) :: aicen(Nx,Ny,Ncat), vicen(Nx,Ny,Ncat), vsnon(Nx,Ny,Ncat)
! real(r8) :: aice(Nx,Ny)

integer  :: i, j, k
integer  :: DimID, VarID, ncid, iunit, io, ndims
real(r8) :: squeeze

! TODO FIXME ... CC ... don't understand the purpose of this.
logical  :: update_restart = .false.

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimLengths
character(len=NF90_MAX_NAME)          :: dimName

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
   write(string1,*) 'cannot open file ', trim(dart_to_cice_input_file),' for updating.'
   call error_handler(E_ERR,'dart_to_cice:filename not found ',trim(dart_to_cice_input_file))
endif

! open file with read and write 
call nc_check( nf90_open(trim(dart_to_cice_input_file), NF90_WRITE, ncid), &
                  'dart_to_cice', 'open '//trim(dart_to_cice_input_file))

! get the key restart variables
! the variables are allocated in the routine

call get_3d_variable(ncid, 'aicen', aicen, dart_to_cice_input_file)
call get_3d_variable(ncid, 'vicen', vicen, dart_to_cice_input_file)
call get_3d_variable(ncid, 'vsnon', vsnon, dart_to_cice_input_file)

Nx   = size(aicen,1)
Ny   = size(aicen,2)
Ncat = size(aicen,3)

allocate( aice(Nx,Ny) )

aice = aicen(:,:,1)
do k = 2, Ncat  
  aice = aice+aicen(:,:,k)
enddo

! fix me if (.not. present(balance_method)) balance_method='simple_squeeze'

SELECT CASE (balance_method)
  CASE ('simple_squeeze')
     do j = 1, Ny
        do i = 1, Nx
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
     do j = 1, Ny
        do i = 1, Nx
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
   aicen(10,10,1)=1.1
   write(*,*) (aicen(10,10,k), k=1,5)

if (update_restart) then
   ! now update the variables
   
   varname='aicen'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, aicen)
   call nc_check(io, 'dart_to_cice', 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))
   
   varname='vicen'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, vicen)
   call nc_check(io, 'dart_to_cice', 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))
   
   varname='vsnon'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, vsnon)
   call nc_check(io, 'dart_to_cice', 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))
   
endif

call nc_check(nf90_close(ncid),'dart_to_cice', 'close '//trim(dart_to_cice_input_file))

call finalize_utilities('dart_to_cice')

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

subroutine get_3d_variable(ncid, varname, var, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
real(r8),         intent(in) :: var(:,:,:)

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimLengths
character(len=NF90_MAX_NAME)          :: dimName

write(msgstring,*) trim(varname)//' '//trim(filename)

io = nf90_inq_varid(ncid, trim(varname), VarID)
call nc_check(io, 'dart_to_cice', 'inq_varid '//trim(msgstring))

io = nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=ndims)
call nc_check(io, 'dart_to_cice', 'inquire_variable '//trim(msgstring))

if (ndims /= 3) then
   write(string2,*) 'expected 3 dimension, got ', ndims
   call error_handler(E_ERR,'dart_to_cice',msgstring,text2=string2)
endif

dimLengths = 1
DimensionLoop : do i = 1,ndims

   write(string1,'(''inquire dimension'',i2,A)') i,trim(msgstring)
   io = nf90_inquire_dimension(ncid, dimIDs(i), name=dimname, len=dimLengths(i))
   call nc_check(io, 'dart_to_cice', string1)

enddo DimensionLoop

allocate( var(dimLengths(1), dimLengths(2), dimLengths(3)) )

call nc_check(nf90_get_var(ncid, VarID, var), 'dart_to_cice', &
         'get_var '//trim(msgstring))

end subroutine get_3d_variable

end program dart_to_cice

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
