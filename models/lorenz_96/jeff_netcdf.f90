! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program jeff_netcdf

!----------------------------------------------------------------------
! purpose: humor jeff
!----------------------------------------------------------------------

use        types_mod, only : r8

use    utilities_mod, only : initialize_utilities, finalize_utilities, nc_check, &
                             open_file, close_file, find_namelist_in_file, &
                             check_namelist_read, nmlfileunit, do_nml_file, do_nml_term, &
                             E_ERR, error_handler, get_unit

use time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                             read_time, get_time, set_time,  &
                             print_date, get_date, &
                             print_time, write_time, &
                             operator(-), operator(+)

use netcdf
use typesizes

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=256) :: string1, string2

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character(len=129)           :: truth_file  = 'true_state.nc'
character(len=129)           :: prior_file  = 'preassim.nc'
character(len=129)           :: poste_file  = 'analysis.nc'
character(len=NF90_MAX_NAME) :: varstring   = 'state'
logical                      :: verbose     = .TRUE.

namelist /jeff_netcdf_nml/ truth_file, prior_file, poste_file, &
                           varstring, verbose

!----------------------------------------------------------------------

type netcdfmetadata
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=NF90_MAX_NAME), dimension(NF90_MAX_VAR_DIMS) :: dimnames
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimids
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
   integer :: ncid
   integer :: varid
   integer :: numdims
   integer :: xtype
   integer :: varsize     ! prod(dimlens(1:numdims))
end type netcdfmetadata

type(netcdfmetadata) :: truthmeta, priormeta, postemeta

real(r8), allocatable, dimension(:,:,:)     :: truth3Dvar  ! time, copy, loc
real(r8), allocatable, dimension(:,:,:)     :: prior3Dvar
real(r8), allocatable, dimension(:,:,:)     :: poste3DVar

real(r8), allocatable, dimension(:,:,:,:)   :: truth4Dvar  ! time, copy, lon, lat
real(r8), allocatable, dimension(:,:,:,:)   :: prior4Dvar
real(r8), allocatable, dimension(:,:,:,:)   :: poste4DVar

real(r8), allocatable, dimension(:,:,:,:,:) :: truth5Dvar  ! time, copy, lon, lat, lev
real(r8), allocatable, dimension(:,:,:,:,:) :: prior5Dvar
real(r8), allocatable, dimension(:,:,:,:,:) :: poste5DVar

integer :: iunit, io

!----------------------------------------------------------------------
!----------------------------------------------------------------------

call initialize_utilities(progname='jeff_netcdf')

write(*,*)
write(*,*)'Reading the namelist to get the input filename.'

call find_namelist_in_file("input.nml", "jeff_netcdf_nml", iunit)
read(iunit, nml = jeff_netcdf_nml, iostat = io)
call check_namelist_read(iunit, io, "jeff_netcdf_nml")

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=jeff_netcdf_nml)
if (do_nml_term()) write(     *     , nml=jeff_netcdf_nml)

! This harvests all kinds of initialization information

truthmeta = netcdfquery(truth_file, varstring)
priormeta = netcdfquery(prior_file, varstring)
postemeta = netcdfquery(poste_file, varstring)

! Get the variable 'truth' 

if (truthmeta%numdims == 3) then
   allocate(truth3Dvar(truthmeta%dimlens(1), truthmeta%dimlens(2), truthmeta%dimlens(3)))
   call nc_check(nf90_get_var(truthmeta%ncid, truthmeta%varid, truth3Dvar), &
                              'humor jeff', 'get_var truth_file')
elseif (truthmeta%numdims == 4) then
   allocate(truth4Dvar(truthmeta%dimlens(1), truthmeta%dimlens(2), &
                       truthmeta%dimlens(3), truthmeta%dimlens(4)))
   call nc_check(nf90_get_var(truthmeta%ncid, truthmeta%varid, truth4Dvar), &
                              'humor jeff', 'get_var truth_file')
elseif (truthmeta%numdims == 5) then
   allocate(truth5Dvar(truthmeta%dimlens(1), truthmeta%dimlens(2), &
                       truthmeta%dimlens(3), truthmeta%dimlens(4), truthmeta%dimlens(5)))
   call nc_check(nf90_get_var(truthmeta%ncid, truthmeta%varid, truth5Dvar), &
                              'humor jeff', 'get_var truth_file')
else
   write(string1,'(a,1x,a,'' is of rank '',i3)')trim(truth_file),trim(varstring),truthmeta%numdims
   write(string2,*)'Only supported ranks are 3,4,5'
   call error_handler(E_ERR,'humor jeff',string1,source,revision,revdate,text2=string2)
endif
call nc_check( nf90_close(truthmeta%ncid), 'humor jeff', 'close '//trim(truth_file))

! Get the variable prior estimate
  
if (priormeta%numdims == 3) then
   allocate(prior3Dvar(priormeta%dimlens(1), priormeta%dimlens(2), priormeta%dimlens(3)))
   call nc_check(nf90_get_var(priormeta%ncid, priormeta%varid, prior3Dvar), &
                              'humor jeff', 'get_var prior_file')
elseif (priormeta%numdims == 4) then
   allocate(prior4Dvar(priormeta%dimlens(1), priormeta%dimlens(2), &
                       priormeta%dimlens(3), priormeta%dimlens(4)))
   call nc_check(nf90_get_var(priormeta%ncid, priormeta%varid, prior4Dvar), &
                              'humor jeff', 'get_var prior_file')
elseif (priormeta%numdims == 5) then
   allocate(prior5Dvar(priormeta%dimlens(1), priormeta%dimlens(2), &
                       priormeta%dimlens(3), priormeta%dimlens(4), priormeta%dimlens(5)))
   call nc_check(nf90_get_var(priormeta%ncid, priormeta%varid, prior5Dvar), &
                              'humor jeff', 'get_var prior_file')
else
   write(string1,'(a,1x,a,'' is of rank '',i3)')trim(prior_file),trim(varstring),priormeta%numdims
   write(string2,*)'Only supported ranks are 3,4,5'
   call error_handler(E_ERR,'humor jeff',string1,source,revision,revdate,text2=string2)
endif
call nc_check( nf90_close(priormeta%ncid), 'humor jeff', 'close '//trim(prior_file))

! Get the variable posterior estimate

if (postemeta%numdims == 3) then
   allocate(poste3Dvar(postemeta%dimlens(1), postemeta%dimlens(2), postemeta%dimlens(3)))
   call nc_check(nf90_get_var(postemeta%ncid, postemeta%varid, poste3Dvar), &
                              'humor jeff', 'get_var poste_file')
elseif (postemeta%numdims == 4) then
   allocate(poste4Dvar(postemeta%dimlens(1), postemeta%dimlens(2), &
                       postemeta%dimlens(3), postemeta%dimlens(4)))
   call nc_check(nf90_get_var(postemeta%ncid, postemeta%varid, poste4Dvar), &
                              'humor jeff', 'get_var poste_file')
elseif (postemeta%numdims == 5) then
   allocate(poste5Dvar(postemeta%dimlens(1), postemeta%dimlens(2), &
                       postemeta%dimlens(3), postemeta%dimlens(4), postemeta%dimlens(5)))
   call nc_check(nf90_get_var(postemeta%ncid, postemeta%varid, poste5Dvar), &
                              'humor jeff', 'get_var poste_file')
else
   write(string1,'(a,1x,a,'' is of rank '',i3)')trim(poste_file),trim(varstring),postemeta%numdims
   write(string2,*)'Only supported ranks are 3,4,5'
   call error_handler(E_ERR,'humor jeff',string1,source,revision,revdate,text2=string2)
endif
call nc_check( nf90_close(postemeta%ncid), 'humor jeff', 'close '//trim(poste_file))

!----------------------------------------------------------------------
! jeff ... here you go ... the useful variables are:
! truth3Dvar(:,:,:) -or- truth4Dvar(:,:,:,:) -or- truth5Dvar(:,:,:,:,:)
! prior3Dvar(:,:,:) -or- prior4Dvar(:,:,:,:) -or- prior5Dvar(:,:,:,:,:)
! poste3Dvar(:,:,:) -or- poste4Dvar(:,:,:,:) -or- poste5Dvar(:,:,:,:,:)
! depending on the location geometry, predominantly.
!----------------------------------------------------------------------

write(*,*)'mean true state ',sum(truth3Dvar)/truthmeta%varsize
write(*,*)'mean prior      ',sum(prior3Dvar(:,3:22,:))/priormeta%varsize
write(*,*)'mean posterior  ',sum(poste3Dvar(:,3:22,:))/postemeta%varsize

call finalize_utilities()

!----------------------------------------------------------------------
!----------------------------------------------------------------------

contains

function netcdfquery( fname, varname ) result(myinfo)
character(len=*), intent(in) :: fname
character(len=*), intent(in) :: varname
type(netcdfmetadata)         :: myinfo

integer :: ncid, i

write(string1,'(a,'' in '',a)') trim(varname),trim(fname)

call nc_check(nf90_open(trim(fname), nf90_nowrite, ncid), 'myinfo','open '//trim(fname))
call nc_check(nf90_inq_varid(ncid, trim(varname), myinfo%varid), 'myinfo','inq_varid '//trim(string1))
call nc_check(nf90_inquire_variable(ncid, myinfo%varid, &
                                   dimids=myinfo%dimids, &
                                    ndims=myinfo%numdims, &
                                    xtype=myinfo%xtype), &
              'myinfo', 'inquire '//trim(string1))

myinfo%varsize = 1
DimensionLoop : do i = 1,myinfo%numdims
   write(string2,'(''inquire dimension'',i2,A)') i,trim(string1)
   call nc_check(nf90_inquire_dimension(ncid, myinfo%dimids(i), &
                                         name=myinfo%dimnames(i), &
                                          len=myinfo%dimlens(i)), &
                                       'myinfo', string2)
   myinfo%varsize = myinfo%varsize * myinfo%dimlens(i)
enddo DimensionLoop

myinfo%ncid    = ncid
myinfo%varname = trim(varname)

! Get supporting attributes if they exist

if( nf90_inquire_attribute(    ncid, myinfo%varid, 'units') == NF90_NOERR ) then
   call nc_check( nf90_get_att(ncid, myinfo%varid, 'units', myinfo%units), &
               'myinfo', 'get_att units '//trim(string1))
else
   myinfo%units = 'none'
endif

if( nf90_inquire_attribute(    ncid, myinfo%varid, 'long_name') == NF90_NOERR ) then
   call nc_check( nf90_get_att(ncid, myinfo%varid, 'long_name', myinfo%long_name), &
               'myinfo', 'get_att long_name '//trim(string1))
else
   myinfo%long_name = varname
endif

if ( verbose ) then

   write(*,*)
   write(*,*)trim(fname)
   write(*,*)'varname              is "',trim(myinfo%varname),'"'
   write(*,*)'long_name        may be "',trim(myinfo%long_name),'"'
   write(*,*)'units            may be "',trim(myinfo%units),'"'
   write(*,*)'variable          ID is ',myinfo%varid
   write(*,*)'xtype             is is ',myinfo%xtype
   write(*,*)'number of dimensions is ',myinfo%numdims
   write(*,*)'total size           is ',myinfo%varsize
   do i = 1,myinfo%numdims
      write(*,'('' dim '',i3,'' is dimid '',i3,'' dimlen '',i10,1x,a)')i,myinfo%dimids(i),myinfo%dimlens(i),trim(myinfo%dimnames(i))
   enddo
   write(*,*)'variable shape is ',myinfo%dimlens(1:myinfo%numdims)

endif

end function netcdfquery

end program jeff_netcdf

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
