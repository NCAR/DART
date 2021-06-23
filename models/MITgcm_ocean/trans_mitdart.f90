!Program to Convert DART-netcdf files to MITgcm-binary and vice-versa.
!Author: S Siva Reddy, sivamtech07@gmail.com/sivareddy.sanikommu@kaust.edu.sa
!Date: 30-Jul-2020
!*******************************************************************************
!When converting from DART-netcdf to MITgcm: 
!  Inputs:
!       MODE.txt		! echo "D2M" > MODE.txt 
!              data             !GRID INFORMATION IS READ FROM HERE
!              INPUT.nc	!PSAL,PTMP,UVEL,VVEL, and ETA data is read from here        
!  Outputs:
!              PSAL.data
!              PTMP.data
!              UVEL.data
!              VVEL.data
!              ETA.data
!
!
!When converting from MITgcm to DART-netcdf : 
!  Inputs:
!       MODE.txt		! echo "M2D" > MODE.txt
!       data		!GRID INFORMATION IS READ FROM HERE
!              PSAL.data
!              PTMP.data
!              UVEL.data
!              VVEL.data
!              ETA.data
!  Outputs:
!       OUTPUT.nc
!*********************************************************************
!Main Program
!---------------------------------------------------------------------

program trans_mitdart

use        types_mod, only : r4, r8, i8, SECPERDAY

use    utilities_mod, only : initialize_utilities, register_module, &
                             get_unit, &
                             error_handler, E_ERR

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'trans_mitdart.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

logical, save :: module_initialized = .false.

character(len=1024) :: msgstring
integer :: iunit

!------------------------------------------------------------------
!
! MITgcm namelist section:  we want to share the 'data' namelist file
! with the model, so we must declare all possible namelist entries to
! avoid getting an error from a valid namelist file.  Most of these
! values are unused in this model_mod code; only a few are needed and
! those are indicated in comments below.
!------------------------------------------------------------------
integer, parameter :: MAX_LEN_FNAM = 512

integer, parameter :: max_nx = 1024
integer, parameter :: max_ny = 1024
integer, parameter :: max_nz = 512
integer, parameter :: max_nr = 512

!--   Gridding parameters variable declarations 
logical :: usingCartesianGrid, usingCylindricalGrid, &
           usingSphericalPolarGrid, usingCurvilinearGrid, &
           deepAtmosphere

real(r8) :: dxSpacing, dySpacing, delX(max_nx), delY(max_ny), &
            ygOrigin, xgOrigin, rSphere, &
            Ro_SeaLevel, delZ(max_nz), delP, delR(max_nr), delRc(max_nr+1), &
            rkFac, groundAtK1

character(len=MAX_LEN_FNAM) :: delXFile, delYFile, &
                               delRFile, delRcFile, &
                               horizGridFile

!--   Gridding parameters namelist
NAMELIST /PARM04/ &
      usingCartesianGrid, usingCylindricalGrid, &
      dxSpacing, dySpacing, delX, delY, delXFile, delYFile, &
      usingSphericalPolarGrid, ygOrigin, xgOrigin, rSphere, &
      usingCurvilinearGrid, horizGridFile, deepAtmosphere, &
      Ro_SeaLevel, delZ, delP, delR, delRc, delRFile, delRcFile, &
      rkFac, groundAtK1

! Grid parameters - the values will be read from a
! standard MITgcm namelist and filled in here.

integer :: Nx=-1, Ny=-1, Nz=-1    ! grid counts for each field

! locations of cell centers (C) and edges (G) for each axis.
real(r8), allocatable :: XC(:), XG(:), YC(:), YG(:), ZC(:), ZG(:)
character(3) :: RWFLAG

!=======================================================================
! Get the party started
!=======================================================================

call initialize_utilities(source)
call register_module(source,revision,revdate)

iunit = get_unit()
open(iunit,file="MODE.txt",status="old")
read(iunit,*)RWFLAG

if (RWFLAG .eq. "D2M") call DART2MIT()
if (RWFLAG .eq. "M2D") call MIT2DART()

close(iunit)

contains

!==================================================================


subroutine static_init_trans()
!------------------------------------------------------------------
!
! Called to do one time initialization of the trans_mitdart. In this case,
! it reads in the grid information and then the model data.

integer :: i, io

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.
module_initialized = .true.

! Grid-related variables are in PARM04
delX(:) = 0.0_r4
delY(:) = 0.0_r4
delZ(:) = 0.0_r4
delR(:) = 0.0_r4

open(unit = 11, file = "data", iostat = io)
if (io /= 0) then
   print *, 'Error while opening ', trim("data")
   STOP 5
endif

read(11, nml = PARM04, iostat = io)
if (io /= 0) then
   print *, 'Error while reading ', trim("data")
   STOP 5
endif

! we use either delR or delZ in mitgcm
if (delR(1) /= 0.0_r4) then
   delZ = delR
endif

! The only way I know to compute the number of
! levels/lats/lons is to set the default value of delZ to 0.0
! before reading the namelist.  now loop until you get back
! to zero and that is the end of the list.
! Not a very satisfying/robust solution ...

Nx = -1
do i=1, size(delX)
 if (delX(i) == 0.0_r4) then
    Nx = i-1
    exit
 endif
enddo
if (Nx == -1) then
   write(msgstring,*)'could not figure out number of longitudes from delX in namelist'
endif

Ny = -1
do i=1, size(delY)
 if (delY(i) == 0.0_r4) then
    Ny = i-1
    exit
 endif
enddo
if (Ny == -1) then
   write(msgstring,*)'could not figure out number of latitudes from delY in namelist'
endif

Nz = -1
do i=1, size(delZ)
 if (delZ(i) == 0.0_r4) then
    Nz = i-1
    exit
 endif
enddo
if (Nz == -1) then
   write(msgstring,*)'could not figure out number of depth levels from delZ in namelist'
endif

! We know enough to allocate grid variables. 

if (.not. allocated(XC)) allocate(XC(Nx))
if (.not. allocated(YC)) allocate(YC(Ny))
if (.not. allocated(ZC)) allocate(ZC(Nz))
if (.not. allocated(XG)) allocate(XG(Nx))
if (.not. allocated(YG)) allocate(YG(Ny))
if (.not. allocated(ZG)) allocate(ZG(Nz))


! XG (the grid edges) and XC (the grid centroids) must be computed.

XG(1) = xgOrigin
XC(1) = xgOrigin + 0.5_r8 * delX(1)
do i=2, Nx
 XG(i) = XG(i-1) + delX(i-1)
 XC(i) = XC(i-1) + 0.5_r8 * delX(i-1) + 0.5_r8 * delX(i) 
enddo

! YG (the grid edges) and YC (the grid centroids) must be computed.

YG(1) = ygOrigin
YC(1) = ygOrigin + 0.5_r8 * delY(1)
do i=2, Ny
 YG(i) = YG(i-1) + delY(i-1)
 YC(i) = YC(i-1) + 0.5_r8 * delY(i-1) + 0.5_r8 * delY(i) 
enddo

! the namelist contains a list of thicknesses of each depth level (delZ)
! ZG (the grid edges) and ZC (the grid centroids) must be computed.

ZG(1) = 0.0_r8
ZC(1) = -0.5_r8 * delZ(1)
do i=2, Nz
 ZG(i) = ZG(i-1) - delZ(i-1)
 ZC(i) = ZC(i-1) - 0.5_r8 * delZ(i-1) - 0.5_r8 * delZ(i) 
enddo

end subroutine static_init_trans

!------------------------------------------------------------------
!> converts the binary input files to a netCDF file

subroutine MIT2DART()

integer  :: ncid, iunit

! for the dimensions and coordinate variables
integer :: XGDimID, XCDimID, YGDimID, YCDimID, ZGDimID, ZCDimID
integer :: XGVarID, XCVarID, YGVarID, YCVarID, ZGVarID, ZCVarID

! for the prognostic variables
integer :: SVarID, TVarID, UVarID, VVarID, EtaVarID 

real(r4), allocatable :: data_3d(:,:,:),data_2d(:,:)
real :: FVAL=-999.0

if ( .not. module_initialized ) call static_init_trans

ALLOCATE(data_3d(Nx,Ny,Nz))
ALLOCATE(data_2d(Nx,Ny))

call check(nf90_create(path="OUTPUT.nc",cmode=NF90_CLOBBER,ncid=ncid))

! Define the new dimensions IDs
   
call check(nf90_def_dim(ncid=ncid, name="XG", len = Nx, dimid = XGDimID))
call check(nf90_def_dim(ncid=ncid, name="XC", len = Nx, dimid = XCDimID))
call check(nf90_def_dim(ncid=ncid, name="YG", len = Ny, dimid = YGDimID))
call check(nf90_def_dim(ncid=ncid, name="YC", len = Ny, dimid = YCDimID))
call check(nf90_def_dim(ncid=ncid, name="ZC", len = Nz, dimid = ZCDimID))
   
! Create the (empty) Coordinate Variables and the Attributes

! U Grid Longitudes

call check(nf90_def_var(ncid,name="XG",xtype=nf90_real,dimids=XGDimID,varid=XGVarID))
call check(nf90_put_att(ncid,  XGVarID, "units", "degrees_east"))
call check(nf90_put_att(ncid,  XGVarID, "modulo", (/ 360.0_r8 /)))
call check(nf90_put_att(ncid,  XGVarID, "point_spacing", "even"))
call check(nf90_put_att(ncid,  XGVarID, "axis", "X"))
call check(nf90_put_att(ncid,  XGVarID, "standard_name", "longitude"))

! S,T,V,Eta Grid Longitudes

call check(nf90_def_var(ncid,name="XC",xtype=nf90_real,dimids=XCDimID,varid=XCVarID))
call check(nf90_put_att(ncid,  XCVarID, "units", "degrees_east"))
call check(nf90_put_att(ncid,  XCVarID, "modulo", (/ 360.0_r8 /)))
call check(nf90_put_att(ncid,  XCVarID, "point_spacing", "even"))
call check(nf90_put_att(ncid,  XCVarID, "axis", "X"))
call check(nf90_put_att(ncid,  XCVarID, "standard_name", "longitude"))

! V Grid Latitudes

call check(nf90_def_var(ncid,name="YG",xtype=nf90_real,dimids=YGDimID,varid=YGVarID))
call check(nf90_put_att(ncid, YGVarID, "units", "degrees_north"))
call check(nf90_put_att(ncid, YGVarID, "point_spacing", "even"))
call check(nf90_put_att(ncid, YGVarID, "axis", "Y"))
call check(nf90_put_att(ncid,YGVarID,"standard_name","latitude"))

! S,T,U,Eta Grid Latitudes

call check(nf90_def_var(ncid,name="YC",xtype=nf90_real,dimids=YCDimID,varid=YCVarID))
call check(nf90_put_att(ncid, YCVarID, "units", "degrees_north"))
call check(nf90_put_att(ncid, YCVarID, "point_spacing", "even"))
call check(nf90_put_att(ncid, YCVarID, "axis", "Y"))
call check(nf90_put_att(ncid,YCVarID,"standard_name","latitude"))

! Depths

call check(nf90_def_var(ncid,name="ZC",xtype=nf90_double,dimids=ZCDimID,varid=ZCVarID))
call check(nf90_put_att(ncid, ZCVarID, "units", "meters"))
call check(nf90_put_att(ncid, ZCVarID, "positive", "up"))
call check(nf90_put_att(ncid, ZCVarID, "point_spacing", "uneven"))
call check(nf90_put_att(ncid, ZCVarID, "axis", "Z"))
call check(nf90_put_att(ncid, ZCVarID, "standard_name", "depth"))

! Create the (empty) Prognostic Variables and the Attributes

call check(nf90_def_var(ncid=ncid, name="PSAL", xtype=nf90_real, &
     dimids = (/XCDimID,YCDimID,ZCDimID/),varid=SVarID))
call check(nf90_put_att(ncid, SVarID, "long_name", "potential salinity"))
call check(nf90_put_att(ncid, SVarID, "missing_value", FVAL))
call check(nf90_put_att(ncid, SVarID, "_FillValue", FVAL))
call check(nf90_put_att(ncid, SVarID, "units", "psu"))
call check(nf90_put_att(ncid, SVarID, "units_long_name", "practical salinity units"))

call check(nf90_def_var(ncid=ncid, name="PTMP", xtype=nf90_real, &
     dimids=(/XCDimID,YCDimID,ZCDimID/),varid=TVarID))
call check(nf90_put_att(ncid, TVarID, "long_name", "Potential Temperature"))
call check(nf90_put_att(ncid, TVarID, "missing_value", FVAL))
call check(nf90_put_att(ncid, TVarID, "_FillValue", FVAL))
call check(nf90_put_att(ncid, TVarID, "units", "C"))
call check(nf90_put_att(ncid, TVarID, "units_long_name", "degrees celsius"))

call check(nf90_def_var(ncid=ncid, name="UVEL", xtype=nf90_real, &
     dimids=(/XGDimID,YCDimID,ZCDimID/),varid=UVarID))
call check(nf90_put_att(ncid, UVarID, "long_name", "Zonal Velocity"))
call check(nf90_put_att(ncid, UVarID, "mssing_value", FVAL))
call check(nf90_put_att(ncid, UVarID, "_FillValue", FVAL))
call check(nf90_put_att(ncid, UVarID, "units", "m/s"))
call check(nf90_put_att(ncid, UVarID, "units_long_name", "meters per second"))

call check(nf90_def_var(ncid=ncid, name="VVEL", xtype=nf90_real, &
     dimids=(/XCDimID,YGDimID,ZCDimID/),varid=VVarID))
call check(nf90_put_att(ncid, VVarID, "long_name", "Meridional Velocity"))
call check(nf90_put_att(ncid, VVarID, "missing_value", FVAL))
call check(nf90_put_att(ncid, VVarID, "_FillValue", FVAL))
call check(nf90_put_att(ncid, VVarID, "units", "m/s"))
call check(nf90_put_att(ncid, VVarID, "units_long_name", "meters per second"))

call check(nf90_def_var(ncid=ncid, name="ETA", xtype=nf90_real, &
     dimids=(/XCDimID,YCDimID/),varid=EtaVarID))
call check(nf90_put_att(ncid, EtaVarID, "long_name", "sea surface height"))
call check(nf90_put_att(ncid, EtaVarID, "missing_value", FVAL))
call check(nf90_put_att(ncid, EtaVarID, "_FillValue", FVAL))
call check(nf90_put_att(ncid, EtaVarID, "units", "meters"))

! Finished with dimension/variable definitions, must end 'define' mode to fill.

call check(nf90_enddef(ncid))

! Fill the coordinate variables

call check(nf90_put_var(ncid, XGVarID, XG ))
call check(nf90_put_var(ncid, XCVarID, XC ))
call check(nf90_put_var(ncid, YGVarID, YG ))
call check(nf90_put_var(ncid, YCVarID, YC ))
call check(nf90_put_var(ncid, ZCVarID, ZC ))

! Fill the data

iunit = get_unit()
open(iunit, file='PSAL.data', form='UNFORMATTED', status='OLD', &
            access='DIRECT', recl=4*Nx*Ny*Nz, convert='BIG_ENDIAN')
read(iunit,rec=1)data_3d
close(iunit)
where (data_3d == 0.0_r4) data_3d = FVAL
call check(nf90_put_var(ncid,SVarID,data_3d,start=(/1,1,1/)))

open(iunit, file='PTMP.data', form='UNFORMATTED', status='OLD', &
            access='DIRECT', recl=4*Nx*Ny*Nz,  convert='BIG_ENDIAN')
read(iunit,rec=1)data_3d
close(iunit)
where (data_3d == 0.0_r4) data_3d = FVAL
call check(nf90_put_var(ncid,TVarID,data_3d,start=(/1,1,1/)))

open(iunit, file='UVEL.data', form='UNFORMATTED', status='OLD', &
            access='DIRECT', recl=4*Nx*Ny*Nz,  convert='BIG_ENDIAN')
read(iunit,rec=1)data_3d
close(iunit)
where (data_3d == 0.0_r4) data_3d = FVAL
call check(nf90_put_var(ncid,UVarID,data_3d,start=(/1,1,1/)))

open(iunit, file='VVEL.data', form='UNFORMATTED', status='OLD', &
            access='DIRECT', recl=4*Nx*Ny*Nz,  convert='BIG_ENDIAN')
read(iunit,rec=1)data_3d
close(iunit)
where (data_3d == 0.0_r4) data_3d = FVAL
call check(nf90_put_var(ncid,VVarID,data_3d,start=(/1,1,1/)))
 
open(iunit, file='ETA.data', form='UNFORMATTED', status='OLD', &
            access='DIRECT', recl=4*Nx*Ny,  convert='BIG_ENDIAN')
read(iunit,rec=1)data_2d
close(iunit)
where (data_2d == 0.0_r4) data_2d = FVAL
call check(nf90_put_var(ncid,EtaVarID,data_2d,start=(/1,1/)))

call check(nf90_close(ncid))

DEALLOCATE(data_3d)
DEALLOCATE(data_2d)

end subroutine MIT2DART

!------------------------------------------------------------------
!> Subroutine for Reading netCDF and writing in binary

subroutine DART2MIT()

integer :: ncid, varid, iunit
real(r4), allocatable :: data_3d(:,:,:),data_2d(:,:)
real :: FVAL

if ( .not. module_initialized ) call static_init_trans

ALLOCATE(data_3d(Nx,Ny,Nz))
ALLOCATE(data_2d(Nx,Ny))

iunit = get_unit()
call check(nf90_open("INPUT.nc",NF90_NOWRITE,ncid))

!Fill the data
call check( NF90_INQ_VARID(ncid,'PSAL',varid) )
call check( NF90_GET_VAR(ncid,varid,data_3d))
call check( nf90_get_att(ncid,varid,"_FillValue",FVAL))        
where (data_3d == FVAL) data_3d = 0.0_r4

open(iunit, file='PSAL.data', form="UNFORMATTED", status='UNKNOWN', &
            access='DIRECT', recl=4*Nx*Ny*Nz, convert='BIG_ENDIAN')
write(iunit,rec=1)data_3d
close(iunit)

call check( NF90_INQ_VARID(ncid,'PTMP',varid) )
call check( NF90_GET_VAR(ncid,varid,data_3d))
call check( nf90_get_att(ncid,varid,"_FillValue",FVAL))        
where (data_3d == FVAL) data_3d = 0.0_r4
open(iunit, file='PTMP.data', form="UNFORMATTED", status='UNKNOWN', &
            access='DIRECT', recl=4*Nx*Ny*Nz, convert='BIG_ENDIAN')
write(iunit,rec=1)data_3d
close(iunit)

call check( NF90_INQ_VARID(ncid,'UVEL',varid) )
call check( NF90_GET_VAR(ncid,varid,data_3d))
call check( nf90_get_att(ncid,varid,"_FillValue",FVAL))        
where (data_3d == FVAL) data_3d = 0.0_r4
open(iunit, file='UVEL.data', form="UNFORMATTED", status='UNKNOWN', &
            access='DIRECT', recl=4*Nx*Ny*Nz, convert='BIG_ENDIAN')
write(iunit,rec=1)data_3d
close(iunit)

call check( NF90_INQ_VARID(ncid,'VVEL',varid) )
call check( NF90_GET_VAR(ncid,varid,data_3d))
call check( nf90_get_att(ncid,varid,"_FillValue",FVAL))        
where (data_3d == FVAL) data_3d = 0.0_r4
open(iunit, file='VVEL.data', form="UNFORMATTED", status='UNKNOWN', &
            access='DIRECT', recl=4*Nx*Ny*Nz, convert='BIG_ENDIAN')
write(iunit,rec=1)data_3d
close(iunit)

call check( NF90_INQ_VARID(ncid,'ETA',varid) )
call check( NF90_GET_VAR(ncid,varid,data_2d))
call check( nf90_get_att(ncid,varid,"_FillValue",FVAL))        
where (data_2d == FVAL) data_2d = 0.0_r4
open(iunit, file='ETA.data', form="UNFORMATTED", status='UNKNOWN', &
            access='DIRECT', recl=4*Nx*Ny, convert='BIG_ENDIAN')
write(iunit,rec=1)data_2d
close(iunit)

call check( NF90_CLOSE(ncid) )

DEALLOCATE(data_3d)
DEALLOCATE(data_2d)

end subroutine DART2MIT

!===============================================================================
!> Subroutine that checks error status on NC file 
!>  Check the error status of the netcdf command

SUBROUTINE check(status)

integer, intent (in) :: status

if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop "Stopped"
end if

END SUBROUTINE check

!===================================================================
! End of trans_mitdart
!===================================================================
end program trans_mitdart

