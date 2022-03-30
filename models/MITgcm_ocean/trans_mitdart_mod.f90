! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download


module trans_mitdart_mod

use types_mod,     only: r4, r8
use utilities_mod, only: initialize_utilities, register_module, &
                         get_unit, find_namelist_in_file, file_exist, &
                         check_namelist_read
use netcdf

implicit none

character(len=*), parameter :: source   = 'trans_mitdart_mod.f90'

logical             :: module_initialized = .false.
character(len=1024) :: msgstring
integer             :: io, iunit

logical             :: do_bgc        = .false.
logical             :: log_transform = .false. 

namelist /trans_mitdart_nml/ do_bgc, log_transform

!------------------------------------------------------------------
!
! MITgcm namelist section:  we want to share the 'data' namelist file
! with the model, so we must declare all possible namelist entries to
! avoid getting an error from a valid namelist file.  Most of these
! values are unused in this model_mod code; only a few are needed and
! those are indicated in comments below.
!------------------------------------------------------------------
integer, parameter :: MAX_LEN_FNAM = 512

integer, parameter :: max_nx = 2048
integer, parameter :: max_ny = 2048
integer, parameter :: max_nz = 512
integer, parameter :: max_nr = 512

!-- record lengths for reading/writing binary files
integer :: recl3d
integer :: recl2d

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

private

public :: static_init_trans, mit2dart, dart2mit

contains

!==================================================================


subroutine static_init_trans()
!------------------------------------------------------------------
!
! Called to do one time initialization of the trans_mitdart. In this case,
! it reads in the grid information and then the model data.

integer :: i, io

if (module_initialized) return

module_initialized = .true.

call find_namelist_in_file('input.nml', 'trans_mitdart_nml', iunit)
read(iunit, nml = trans_mitdart_nml, iostat = io)
call check_namelist_read(iunit, io, 'trans_mitdart_nml')


! Grid-related variables are in PARM04
delX(:) = 0.0_r4
delY(:) = 0.0_r4
delZ(:) = 0.0_r4
delR(:) = 0.0_r4

call find_namelist_in_file('data', 'PARM04', iunit)
read(iunit, nml = PARM04, iostat = io)
call check_namelist_read(iunit, io, 'data')

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

! set record lengths
recl3d = Nx*Ny*Nz*4
recl2d = Nx*Ny*4

! MEG Better have that as inout namelist parameter
! Are we also doing bgc on top of physics?
! If we found nitrate then the rest of the binaries (for the 
! remaining 9 variables) should be also there.
! TODO may also enhance this functionality
! if (file_exist('NO3.data')) do_bgc = .true.

end subroutine static_init_trans

!------------------------------------------------------------------
!> converts the binary input files to a netCDF file

subroutine mit2dart()

integer  :: ncid, iunit

! for the dimensions and coordinate variables
integer :: XGDimID, XCDimID, YGDimID, YCDimID, ZGDimID, ZCDimID
integer :: XGVarID, XCVarID, YGVarID, YCVarID, ZGVarID, ZCVarID

! for the prognostic variables
integer :: SVarID, TVarID, UVarID, VVarID, EtaVarID
integer :: no3_varid, po4_varid, o2_varid, phy_varid, alk_varid 
integer :: dic_varid, dop_varid, don_varid, fet_varid

! diagnostic variable
integer :: chl_varid  

real(r4), allocatable :: data_3d(:,:,:), data_2d(:,:)

real(r4) :: FVAL

if (.not. module_initialized) call static_init_trans

FVAL=-999.0_r4

allocate(data_3d(Nx,Ny,Nz))
allocate(data_2d(Nx,Ny))

call check(nf90_create(path="OUTPUT.nc",cmode=or(nf90_clobber,nf90_64bit_offset),ncid=ncid))

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
call check(nf90_put_att(ncid, EtaVarID, "units", "m"))
call check(nf90_put_att(ncid, EtaVarID, "units_long_name", "meters"))

!> Add BLING data:

if (do_bgc) then 
   ! 1. BLING tracer: nitrate NO3
   call check(nf90_def_var(ncid=ncid, name="NO3", xtype=nf90_real, &
        dimids=(/XCDimID,YCDimID,ZCDimID/),varid=no3_varid))
   call check(nf90_put_att(ncid, no3_varid, "long_name"      , "Nitrate"))
   call check(nf90_put_att(ncid, no3_varid, "missing_value"  , FVAL))
   call check(nf90_put_att(ncid, no3_varid, "_FillValue"     , FVAL))
   call check(nf90_put_att(ncid, no3_varid, "units"          , "mol N/m3"))
   call check(nf90_put_att(ncid, no3_varid, "units_long_name", "moles Nitrogen per cubic meters"))
   
   ! 2. BLING tracer: phosphate PO4
   call check(nf90_def_var(ncid=ncid, name="PO4", xtype=nf90_real, &
        dimids=(/XCDimID,YCDimID,ZCDimID/),varid=po4_varid))
   call check(nf90_put_att(ncid, po4_varid, "long_name"      , "Phosphate"))
   call check(nf90_put_att(ncid, po4_varid, "missing_value"  , FVAL))
   call check(nf90_put_att(ncid, po4_varid, "_FillValue"     , FVAL))
   call check(nf90_put_att(ncid, po4_varid, "units"          , "mol P/m3"))
   call check(nf90_put_att(ncid, po4_varid, "units_long_name", "moles Phosphorus per cubic meters"))
   
   ! 3. BLING tracer: oxygen O2
   call check(nf90_def_var(ncid=ncid, name="O2", xtype=nf90_real, &
        dimids=(/XCDimID,YCDimID,ZCDimID/),varid=o2_varid))
   call check(nf90_put_att(ncid, o2_varid, "long_name"      , "Dissolved Oxygen"))
   call check(nf90_put_att(ncid, o2_varid, "missing_value"  , FVAL))
   call check(nf90_put_att(ncid, o2_varid, "_FillValue"     , FVAL))
   call check(nf90_put_att(ncid, o2_varid, "units"          , "mol O/m3"))
   call check(nf90_put_att(ncid, o2_varid, "units_long_name", "moles Oxygen per cubic meters"))
   
   ! 4. BLING tracer: phytoplankton PHY
   call check(nf90_def_var(ncid=ncid, name="PHY", xtype=nf90_real, &
        dimids=(/XCDimID,YCDimID,ZCDimID/),varid=phy_varid))
   call check(nf90_put_att(ncid, phy_varid, "long_name"      , "Phytoplankton Biomass"))
   call check(nf90_put_att(ncid, phy_varid, "missing_value"  , FVAL))
   call check(nf90_put_att(ncid, phy_varid, "_FillValue"     , FVAL))
   call check(nf90_put_att(ncid, phy_varid, "units"          , "mol C/m3"))
   call check(nf90_put_att(ncid, phy_varid, "units_long_name", "moles Carbon per cubic meters"))
   
   ! 5. BLING tracer: alkalinity ALK
   call check(nf90_def_var(ncid=ncid, name="ALK", xtype=nf90_real, &
        dimids=(/XCDimID,YCDimID,ZCDimID/),varid=alk_varid))
   call check(nf90_put_att(ncid, alk_varid, "long_name"      , "Alkalinity"))
   call check(nf90_put_att(ncid, alk_varid, "missing_value"  , FVAL))
   call check(nf90_put_att(ncid, alk_varid, "_FillValue"     , FVAL))
   call check(nf90_put_att(ncid, alk_varid, "units"          , "mol eq/m3"))
   call check(nf90_put_att(ncid, alk_varid, "units_long_name", "moles equivalent per cubic meters"))
   
   ! 6. BLING tracer: dissolved inorganic carbon DIC
   call check(nf90_def_var(ncid=ncid, name="DIC", xtype=nf90_real, &
        dimids=(/XCDimID,YCDimID,ZCDimID/),varid=dic_varid))
   call check(nf90_put_att(ncid, dic_varid, "long_name"      , "Dissolved Inorganic Carbon"))
   call check(nf90_put_att(ncid, dic_varid, "missing_value"  , FVAL))
   call check(nf90_put_att(ncid, dic_varid, "_FillValue"     , FVAL))
   call check(nf90_put_att(ncid, dic_varid, "units"          , "mol C/m3"))
   call check(nf90_put_att(ncid, dic_varid, "units_long_name", "moles Carbon per cubic meters"))
   
   ! 7. BLING tracer: dissolved organic phosphorus DOP
   call check(nf90_def_var(ncid=ncid, name="DOP", xtype=nf90_real, &
        dimids=(/XCDimID,YCDimID,ZCDimID/),varid=dop_varid))
   call check(nf90_put_att(ncid, dop_varid, "long_name"      , "Dissolved Organic Phosphorus"))
   call check(nf90_put_att(ncid, dop_varid, "missing_value"  , FVAL))
   call check(nf90_put_att(ncid, dop_varid, "_FillValue"     , FVAL))
   call check(nf90_put_att(ncid, dop_varid, "units"          , "mol P/m3"))
   call check(nf90_put_att(ncid, dop_varid, "units_long_name", "moles Phosphorus per cubic meters"))
   
   ! 8. BLING tracer: dissolved organic nitrogen DON
   call check(nf90_def_var(ncid=ncid, name="DON", xtype=nf90_real, &
        dimids=(/XCDimID,YCDimID,ZCDimID/),varid=don_varid))
   call check(nf90_put_att(ncid, don_varid, "long_name"      , "Dissolved Organic Nitrogen"))
   call check(nf90_put_att(ncid, don_varid, "missing_value"  , FVAL))
   call check(nf90_put_att(ncid, don_varid, "_FillValue"     , FVAL))
   call check(nf90_put_att(ncid, don_varid, "units"          , "mol N/m3"))
   call check(nf90_put_att(ncid, don_varid, "units_long_name", "moles Nitrogen per cubic meters"))
   
   ! 9. BLING tracer: dissolved inorganic iron FET
   call check(nf90_def_var(ncid=ncid, name="FET", xtype=nf90_real, &
        dimids=(/XCDimID,YCDimID,ZCDimID/),varid=fet_varid))
   call check(nf90_put_att(ncid, fet_varid, "long_name"      , "Dissolved Inorganic Iron"))
   call check(nf90_put_att(ncid, fet_varid, "missing_value"  , FVAL))
   call check(nf90_put_att(ncid, fet_varid, "_FillValue"     , FVAL))
   call check(nf90_put_att(ncid, fet_varid, "units"          , "mol Fe/m3"))
   call check(nf90_put_att(ncid, fet_varid, "units_long_name", "moles Iron per cubic meters"))
   
   ! 10. BLING tracer: Surface Chlorophyl CHL
   call check(nf90_def_var(ncid=ncid, name="CHL", xtype=nf90_real, &
        dimids=(/XCDimID,YCDimID/),varid=chl_varid))
   call check(nf90_put_att(ncid, chl_varid, "long_name"      , "Surface Chlorophyll"))
   call check(nf90_put_att(ncid, chl_varid, "missing_value"  , FVAL))
   call check(nf90_put_att(ncid, chl_varid, "_FillValue"     , FVAL))
   call check(nf90_put_att(ncid, chl_varid, "units"          , "mg/m3"))
   call check(nf90_put_att(ncid, chl_varid, "units_long_name", "milligram per cubic meters"))
endif   

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
            access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
read(iunit,rec=1)data_3d
close(iunit)
where (data_3d == 0.0_r4) data_3d = FVAL
call check(nf90_put_var(ncid,SVarID,data_3d,start=(/1,1,1/)))

open(iunit, file='PTMP.data', form='UNFORMATTED', status='OLD', &
            access='DIRECT', recl=recl3d,  convert='BIG_ENDIAN')
read(iunit,rec=1)data_3d
close(iunit)
where (data_3d == 0.0_r4) data_3d = FVAL
call check(nf90_put_var(ncid,TVarID,data_3d,start=(/1,1,1/)))

open(iunit, file='UVEL.data', form='UNFORMATTED', status='OLD', &
            access='DIRECT', recl=recl3d,  convert='BIG_ENDIAN')
read(iunit,rec=1)data_3d
close(iunit)
where (data_3d == 0.0_r4) data_3d = FVAL
call check(nf90_put_var(ncid,UVarID,data_3d,start=(/1,1,1/)))

open(iunit, file='VVEL.data', form='UNFORMATTED', status='OLD', &
            access='DIRECT', recl=recl3d,  convert='BIG_ENDIAN')
read(iunit,rec=1)data_3d
close(iunit)
where (data_3d == 0.0_r4) data_3d = FVAL
call check(nf90_put_var(ncid,VVarID,data_3d,start=(/1,1,1/)))
 
open(iunit, file='ETA.data', form='UNFORMATTED', status='OLD', &
            access='DIRECT', recl=recl2d,  convert='BIG_ENDIAN')
read(iunit,rec=1)data_2d
close(iunit)
where (data_2d == 0.0_r4) data_2d = FVAL
call check(nf90_put_var(ncid,EtaVarID,data_2d,start=(/1,1/)))

if (do_bgc) then 
   open(iunit, file='NO3.data', form='UNFORMATTED', status='OLD', &
               access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
   read(iunit,rec=1)data_3d
   close(iunit)
   call fill_var_md(data_3d, FVAL)
   call check(nf90_put_var(ncid,no3_varid,data_3d,start=(/1,1,1/)))
   
   open(iunit, file='PO4.data', form='UNFORMATTED', status='OLD', &
               access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
   read(iunit,rec=1)data_3d
   close(iunit)
   call fill_var_md(data_3d, FVAL)
   call check(nf90_put_var(ncid,po4_varid,data_3d,start=(/1,1,1/)))
   
   open(iunit, file='O2.data', form='UNFORMATTED', status='OLD', &
               access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
   read(iunit,rec=1)data_3d
   close(iunit)
   call fill_var_md(data_3d, FVAL)
   call check(nf90_put_var(ncid,o2_varid,data_3d,start=(/1,1,1/)))
   
   open(iunit, file='PHY.data', form='UNFORMATTED', status='OLD', &
               access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
   read(iunit,rec=1)data_3d
   close(iunit)
   call fill_var_md(data_3d, FVAL)
   call check(nf90_put_var(ncid,phy_varid,data_3d,start=(/1,1,1/)))
   
   open(iunit, file='ALK.data', form='UNFORMATTED', status='OLD', &
               access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
   read(iunit,rec=1)data_3d
   close(iunit)
   call fill_var_md(data_3d, FVAL)
   call check(nf90_put_var(ncid,alk_varid,data_3d,start=(/1,1,1/)))
   
   open(iunit, file='DIC.data', form='UNFORMATTED', status='OLD', &
               access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
   read(iunit,rec=1)data_3d
   close(iunit)
   call fill_var_md(data_3d, FVAL)
   call check(nf90_put_var(ncid,dic_varid,data_3d,start=(/1,1,1/)))
   
   open(iunit, file='DOP.data', form='UNFORMATTED', status='OLD', &
               access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
   read(iunit,rec=1)data_3d
   close(iunit)
   call fill_var_md(data_3d, FVAL)
   call check(nf90_put_var(ncid,dop_varid,data_3d,start=(/1,1,1/)))
   
   open(iunit, file='DON.data', form='UNFORMATTED', status='OLD', &
               access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
   read(iunit,rec=1)data_3d
   close(iunit)
   call fill_var_md(data_3d, FVAL)
   call check(nf90_put_var(ncid,don_varid,data_3d,start=(/1,1,1/)))
   
   open(iunit, file='FET.data', form='UNFORMATTED', status='OLD', &
               access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
   read(iunit,rec=1)data_3d
   close(iunit)
   call fill_var_md(data_3d, FVAL)
   call check(nf90_put_var(ncid,fet_varid,data_3d,start=(/1,1,1/)))
   
   open(iunit, file='CHL.data', form='UNFORMATTED', status='OLD', &
               access='DIRECT', recl=recl2d,  convert='BIG_ENDIAN')
   read(iunit,rec=1)data_2d
   close(iunit)
   where (data_2d == 0.0_r4) 
       data_2d = FVAL
   elsewhere
       data_2d = log10(data_2d)
   endwhere 
   call check(nf90_put_var(ncid,chl_varid,data_2d,start=(/1,1/)))
endif

call check(nf90_close(ncid))

deallocate(data_3d)
deallocate(data_2d)

end subroutine mit2dart

!------------------------------------------------------------------
!> Subroutine for Reading netCDF and writing in binary

subroutine dart2mit()

integer :: ncid, varid, iunit
real(r4), allocatable :: data_3d(:,:,:),data_2d(:,:)
real(r4) :: FVAL

allocate(data_3d(Nx,Ny,Nz))
allocate(data_2d(Nx,Ny))

if (.not. module_initialized) call static_init_trans

iunit = get_unit()
call check(nf90_open("INPUT.nc",NF90_NOWRITE,ncid))

!Fill the data
call check( NF90_INQ_VARID(ncid,'PSAL',varid) )
call check( NF90_GET_VAR(ncid,varid,data_3d))
call check( nf90_get_att(ncid,varid,"_FillValue",FVAL))        
where (data_3d == FVAL) data_3d = 0.0_r4

open(iunit, file='PSAL.data', form="UNFORMATTED", status='UNKNOWN', &
            access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
write(iunit,rec=1)data_3d
close(iunit)

call check( NF90_INQ_VARID(ncid,'PTMP',varid) )
call check( NF90_GET_VAR(ncid,varid,data_3d))
call check( nf90_get_att(ncid,varid,"_FillValue",FVAL))        
where (data_3d == FVAL) data_3d = 0.0_r4

open(iunit, file='PTMP.data', form="UNFORMATTED", status='UNKNOWN', &
            access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
write(iunit,rec=1)data_3d
close(iunit)

call check( NF90_INQ_VARID(ncid,'UVEL',varid) )
call check( NF90_GET_VAR(ncid,varid,data_3d))
call check( nf90_get_att(ncid,varid,"_FillValue",FVAL))        
where (data_3d == FVAL) data_3d = 0.0_r4

open(iunit, file='UVEL.data', form="UNFORMATTED", status='UNKNOWN', &
            access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
write(iunit,rec=1)data_3d
close(iunit)

call check( NF90_INQ_VARID(ncid,'VVEL',varid) )
call check( NF90_GET_VAR(ncid,varid,data_3d))
call check( nf90_get_att(ncid,varid,"_FillValue",FVAL))        
where (data_3d == FVAL) data_3d = 0.0_r4

open(iunit, file='VVEL.data', form="UNFORMATTED", status='UNKNOWN', &
            access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
write(iunit,rec=1)data_3d
close(iunit)

call check( NF90_INQ_VARID(ncid,'ETA',varid) )
call check( NF90_GET_VAR(ncid,varid,data_2d))
call check( nf90_get_att(ncid,varid,"_FillValue",FVAL))        
where (data_2d == FVAL) data_2d = 0.0_r4

open(iunit, file='ETA.data', form="UNFORMATTED", status='UNKNOWN', &
            access='DIRECT', recl=recl2d, convert='BIG_ENDIAN')
write(iunit,rec=1)data_2d
close(iunit)

if (do_bgc) then 
   call check( NF90_INQ_VARID(ncid,'NO3',varid) )
   call check( NF90_GET_VAR(ncid,varid,data_3d))
   call check( nf90_get_att(ncid,varid,"_FillValue",FVAL))
   call fill_var_dm(data_3d, FVAL)
 
   open(iunit, file='NO3.data', form="UNFORMATTED", status='UNKNOWN', &
               access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
   write(iunit,rec=1)data_3d
   close(iunit)
   
   call check( NF90_INQ_VARID(ncid,'PO4',varid) )
   call check( NF90_GET_VAR(ncid,varid,data_3d))
   call check( nf90_get_att(ncid,varid,"_FillValue",FVAL))
   call fill_var_dm(data_3d, FVAL)
 
   open(iunit, file='PO4.data', form="UNFORMATTED", status='UNKNOWN', &
               access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
   write(iunit,rec=1)data_3d
   close(iunit)
   
   call check( NF90_INQ_VARID(ncid,'O2',varid) )
   call check( NF90_GET_VAR(ncid,varid,data_3d))
   call check( nf90_get_att(ncid,varid,"_FillValue",FVAL))
   call fill_var_dm(data_3d, FVAL)
 
   open(iunit, file='O2.data', form="UNFORMATTED", status='UNKNOWN', &
               access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
   write(iunit,rec=1)data_3d
   close(iunit)
   
   call check( NF90_INQ_VARID(ncid,'PHY',varid) )
   call check( NF90_GET_VAR(ncid,varid,data_3d))
   call check( nf90_get_att(ncid,varid,"_FillValue",FVAL))
   call fill_var_dm(data_3d, FVAL)
 
   open(iunit, file='PHY.data', form="UNFORMATTED", status='UNKNOWN', &
               access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
   write(iunit,rec=1)data_3d
   close(iunit)
   
   call check( NF90_INQ_VARID(ncid,'ALK',varid) )
   call check( NF90_GET_VAR(ncid,varid,data_3d))
   call check( nf90_get_att(ncid,varid,"_FillValue",FVAL))
   call fill_var_dm(data_3d, FVAL)
 
   open(iunit, file='ALK.data', form="UNFORMATTED", status='UNKNOWN', &
               access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
   write(iunit,rec=1)data_3d
   close(iunit)
   
   call check( NF90_INQ_VARID(ncid,'DIC',varid) )
   call check( NF90_GET_VAR(ncid,varid,data_3d))
   call check( nf90_get_att(ncid,varid,"_FillValue",FVAL))
   call fill_var_dm(data_3d, FVAL)
 
   open(iunit, file='DIC.data', form="UNFORMATTED", status='UNKNOWN', &
               access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
   write(iunit,rec=1)data_3d
   close(iunit)
   
   call check( NF90_INQ_VARID(ncid,'DOP',varid) )
   call check( NF90_GET_VAR(ncid,varid,data_3d))
   call check( nf90_get_att(ncid,varid,"_FillValue",FVAL))
   call fill_var_dm(data_3d, FVAL)
 
   open(iunit, file='DOP.data', form="UNFORMATTED", status='UNKNOWN', &
               access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
   write(iunit,rec=1)data_3d
   close(iunit)
   
   call check( NF90_INQ_VARID(ncid,'DON',varid) )
   call check( NF90_GET_VAR(ncid,varid,data_3d))
   call check( nf90_get_att(ncid,varid,"_FillValue",FVAL))
   call fill_var_dm(data_3d, FVAL)
 
   open(iunit, file='DON.data', form="UNFORMATTED", status='UNKNOWN', &
               access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
   write(iunit,rec=1)data_3d
   close(iunit)
   
   call check( NF90_INQ_VARID(ncid,'FET',varid) )
   call check( NF90_GET_VAR(ncid,varid,data_3d))
   call check( nf90_get_att(ncid,varid,"_FillValue",FVAL))
   call fill_var_dm(data_3d, FVAL)
 
   open(iunit, file='FET.data', form="UNFORMATTED", status='UNKNOWN', &
               access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
   write(iunit,rec=1)data_3d
   close(iunit)
endif

call check( NF90_CLOSE(ncid) )

deallocate(data_3d)
deallocate(data_2d)

end subroutine dart2mit

!===============================================================================
!> Subroutine that checks error status on NC file 
!>  Check the error status of the netcdf command

subroutine check(status)

integer, intent (in) :: status

if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop "Stopped"
end if

end subroutine check


!===============================================================================
!> Check the tracer variables after reading from the binaries
!> Make sure they are non-negative
!> Do the transform if requested
!> md: mit2dart; dm: dart2mit

subroutine fill_var_md(var, fillval)

real(r4), intent(inout) :: var(:, :, :)
real(r4), intent(in)    :: fillval

real(r4) :: low_conc

if (.not. module_initialized) call static_init_trans

low_conc = 1.0e-12

! Make sure the tracer concentration is positive 
where(var < 0.0_r4) var = low_conc

if (log_transform) then
   where (var == 0.0_r4)
       var = fillval
   elsewhere
       var = log(var)
   endwhere
else
   where (var == 0.0_r4) var = fillval
endif

end subroutine

!------------------------------------------------------------------

subroutine fill_var_dm(var, fillval)

real(r4), intent(inout) :: var(:, :, :)
real(r4), intent(in)    :: fillval

if (.not. module_initialized) call static_init_trans

if (log_transform) then
   where (var == fillval)
       var = 0.0_r4
   elsewhere
       var = exp(var)
   endwhere
else
   where (var == fillval) var = 0.0_r4
endif

end subroutine

!------------------------------------------------------------------

end module trans_mitdart_mod

