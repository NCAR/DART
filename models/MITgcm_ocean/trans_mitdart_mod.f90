! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download


module trans_mitdart_mod

use types_mod,     only: r4, r8
use utilities_mod, only: initialize_utilities, register_module, &
                         get_unit, find_namelist_in_file, file_exist, &
                         check_namelist_read
use netcdf_utilities_mod, only : nc_get_variable, nc_get_dimension_size
use netcdf

implicit none

character(len=*), parameter :: source   = 'trans_mitdart_mod.f90'

logical             :: module_initialized = .false.
character(len=1024) :: msgstring
integer             :: io, iunit

logical             :: do_bgc        = .false.
logical             :: log_transform = .false.
logical             :: compress      = .false.
! set compress = .true. remove missing values from state
logical             :: output_chl_data = .false.
! CHL.data is not written to mit .data files by default

namelist /trans_mitdart_nml/ do_bgc, log_transform, compress

real(r4), parameter :: FVAL=-999.0_r4 ! may put this as a namelist option
real(r4), parameter :: binary_fill=0.0_r4

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
integer :: ncomp2 = -1  ! length of 2D compressed dim
integer :: ncomp3 = -1, ncomp3U = -1, ncomp3V = -1   ! length of 3D compressed dim

integer, parameter :: MITgcm_3D_FIELD   = 1
integer, parameter :: MITgcm_3D_FIELD_U = 2
integer, parameter :: MITgcm_3D_FIELD_V = 3

! locations of cell centers (C) and edges (G) for each axis.
real(r8), allocatable :: XC(:), XG(:), YC(:), YG(:), ZC(:), ZG(:)
real(r8), allocatable :: XCcomp(:), XGcomp(:), YCcomp(:), YGcomp(:), ZCcomp(:), ZGcomp(:)

integer, allocatable  :: Xcomp_ind(:), Ycomp_ind(:), Zcomp_ind(:) !HK are the staggered grids compressed the same?
!MEG: For staggered grids
integer, allocatable  :: Xcomp_indU(:), Ycomp_indU(:), Zcomp_indU(:)    
integer, allocatable  :: Xcomp_indV(:), Ycomp_indV(:), Zcomp_indV(:) 

! 3D variables, 3 grids:
!
! XC, YC, ZC  1  PSAL, PTMP, NO3, PO4, O2, PHY, ALK, DIC, DOP, DON, FET
! XC, YC, ZG  2  UVEL
! XC, YG, ZC  3  VVEL

! MEG: For compression, especially if we're doing Arakawa C-grid, 
! we will need 3 different compressions for the above variables

! 2D variables, 1 grid:
!
! YC, XC  ETA, CHL

private

public :: static_init_trans, mit2dart, dart2mit

interface write_compressed
   module procedure write_compressed_2d
   module procedure write_compressed_3d
end interface write_compressed

interface read_compressed
   module procedure read_compressed_2d
   module procedure read_compressed_3d
end interface read_compressed

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

end subroutine static_init_trans

!------------------------------------------------------------------
!> converts the binary input files to a netCDF file

subroutine mit2dart()

integer  :: ncid

! for the dimensions and coordinate variables
integer :: XGDimID, XCDimID, YGDimID, YCDimID, ZGDimID, ZCDimID
integer :: XGVarID, XCVarID, YGVarID, YCVarID, ZGVarID, ZCVarID
integer :: comp2ID, comp3ID, comp3UD, comp3VD ! compressed dim
integer :: XGcompVarID, XCcompVarID, YGcompVarID, YCcompVarID, ZGcompVarID, ZCcompVarID
integer :: XindID, YindID, ZindID
integer :: XindUD, YindUD, ZindUD
integer :: XindVD, YindVD, ZindVD
integer :: all_dimids(9) ! store the 9 dimension ids that are used

! for the prognostic variables
integer :: SVarID, TVarID, UVarID, VVarID, EtaVarID
integer :: no3_varid, po4_varid, o2_varid, phy_varid, alk_varid 
integer :: dic_varid, dop_varid, don_varid, fet_varid

! diagnostic variable
integer :: chl_varid  

if (.not. module_initialized) call static_init_trans

call check(nf90_create(path="OUTPUT.nc",cmode=or(nf90_clobber,nf90_64bit_offset),ncid=ncid))

! Define the new dimensions IDs

call check(nf90_def_dim(ncid=ncid, name="XC", len = Nx, dimid = XCDimID))
call check(nf90_def_dim(ncid=ncid, name="YC", len = Ny, dimid = YCDimID))
call check(nf90_def_dim(ncid=ncid, name="ZC", len = Nz, dimid = ZCDimID))

call check(nf90_def_dim(ncid=ncid, name="XG", len = Nx, dimid = XGDimID))
call check(nf90_def_dim(ncid=ncid, name="YG", len = Ny, dimid = YGDimID))

print *, ''

if (compress) then
   ncomp2 = get_compressed_size_2d()

   write(*, '(A, I12, A, I8)') '2D: ', Nx*Ny, ', COMP2D: ', ncomp2

   ncomp3  = get_compressed_size_3d(MITgcm_3D_FIELD)
   ncomp3U = get_compressed_size_3d(MITgcm_3D_FIELD_U)
   ncomp3V = get_compressed_size_3d(MITgcm_3D_FIELD_V) 

   write(*, '(A, I12, A, 3I8)') '3D: ', Nx*Ny*Nz, ', COMP3D [T-S, U, V]: ', ncomp3, ncomp3U, ncomp3V

   ! Put the compressed dimensions in the restart file
   call check(nf90_def_dim(ncid=ncid, name="comp2d",  len = ncomp2,  dimid = comp2ID))
   call check(nf90_def_dim(ncid=ncid, name="comp3d",  len = ncomp3,  dimid = comp3ID))
   call check(nf90_def_dim(ncid=ncid, name="comp3dU", len = ncomp3U, dimid = comp3UD))
   call check(nf90_def_dim(ncid=ncid, name="comp3dV", len = ncomp3V, dimid = comp3VD)) 
else
  comp2ID = -1
  comp3ID = -1
endif

all_dimids = (/XCDimID, YCDimID, ZCDimID, XGDimID, YGDimID, &
               comp2ID, comp3ID, comp3UD, comp3VD/)

print *, ''

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

! Compressed grid variables
if (compress) then
   call check(nf90_def_var(ncid,name="XGcomp",xtype=nf90_real,dimids=comp3ID,varid=XGcompVarID))
   call check(nf90_def_var(ncid,name="XCcomp",xtype=nf90_real,dimids=comp3ID,varid=XCcompVarID))
   call check(nf90_def_var(ncid,name="YGcomp",xtype=nf90_real,dimids=comp3ID,varid=YGcompVarID))
   call check(nf90_def_var(ncid,name="YCcomp",xtype=nf90_real,dimids=comp3ID,varid=YCcompVarID))
   call check(nf90_def_var(ncid,name="ZCcomp",xtype=nf90_double,dimids=comp3ID,varid=ZCcompVarID))

   call check(nf90_def_var(ncid,name="Xcomp_ind",xtype=nf90_int,dimids=comp3ID,varid=XindID))
   call check(nf90_def_var(ncid,name="Ycomp_ind",xtype=nf90_int,dimids=comp3ID,varid=YindID))
   call check(nf90_def_var(ncid,name="Zcomp_ind",xtype=nf90_int,dimids=comp3ID,varid=ZindID))

   call check(nf90_def_var(ncid,name="Xcomp_indU",xtype=nf90_int,dimids=comp3UD,varid=XindUD))
   call check(nf90_def_var(ncid,name="Ycomp_indU",xtype=nf90_int,dimids=comp3UD,varid=YindUD))
   call check(nf90_def_var(ncid,name="Zcomp_indU",xtype=nf90_int,dimids=comp3UD,varid=ZindUD))

   call check(nf90_def_var(ncid,name="Xcomp_indV",xtype=nf90_int,dimids=comp3VD,varid=XindVD))
   call check(nf90_def_var(ncid,name="Ycomp_indV",xtype=nf90_int,dimids=comp3VD,varid=YindVD))
   call check(nf90_def_var(ncid,name="Zcomp_indV",xtype=nf90_int,dimids=comp3VD,varid=ZindVD))
endif

! The size of these variables will depend on the compression

! Create the (empty) Prognostic Variables and the Attributes

SVarID = define_variable(ncid,"PSAL", nf90_real, all_dimids, MITgcm_3D_FIELD)
call add_attributes_to_variable(ncid, SVarID, "potential salinity", "psu", "practical salinity units")

TVarID = define_variable(ncid,"PTMP", nf90_real, all_dimids, MITgcm_3D_FIELD)
call add_attributes_to_variable(ncid, TVarID, "Potential Temperature", "C", "degrees celsius")

UVarID = define_variable(ncid,"UVEL", nf90_real, all_dimids, MITgcm_3D_FIELD_U)
call add_attributes_to_variable(ncid, UVarID, "Zonal Velocity", "m/s", "meters per second")

VVarID = define_variable(ncid,"VVEL", nf90_real, all_dimids, MITgcm_3D_FIELD_V)
call add_attributes_to_variable(ncid, VVarID, "Meridional Velocity", "m/s", "meters per second")

EtaVarID = define_variable_2d(ncid,"ETA", nf90_real, all_dimids)
call add_attributes_to_variable(ncid, EtaVarID, "sea surface height", "m", "meters")

! Create the BLING netcdf variables:

if (do_bgc) then 
   ! 1. BLING tracer: nitrate NO3
   no3_varid = define_variable(ncid,"NO3", nf90_real, all_dimids, MITgcm_3D_FIELD)
   call add_attributes_to_variable(ncid, no3_varid, "Nitrate", "mol N/m3", "moles Nitrogen per cubic meters")
     
   ! 2. BLING tracer: phosphate PO4
   po4_varid = define_variable(ncid,"PO4", nf90_real, all_dimids, MITgcm_3D_FIELD)
   call add_attributes_to_variable(ncid, po4_varid, "Phosphate", "mol P/m3", "moles Phosphorus per cubic meters")
  
   ! 3. BLING tracer: oxygen O2
   o2_varid = define_variable(ncid,"O2", nf90_real, all_dimids, MITgcm_3D_FIELD)
   call add_attributes_to_variable(ncid, o2_varid, "Dissolved Oxygen", "mol O/m3", "moles Oxygen per cubic meters")
  
   ! 4. BLING tracer: phytoplankton PHY
   phy_varid = define_variable(ncid,"PHY", nf90_real, all_dimids, MITgcm_3D_FIELD)
   call add_attributes_to_variable(ncid, phy_varid, "Phytoplankton Biomass", "mol C/m3", "moles Carbon per cubic meters")
   
   ! 5. BLING tracer: alkalinity ALK
   alk_varid = define_variable(ncid,"ALK", nf90_real, all_dimids, MITgcm_3D_FIELD)
   call add_attributes_to_variable(ncid, alk_varid, "Alkalinity", "mol eq/m3", "moles equivalent per cubic meters")

   ! 6. BLING tracer: dissolved inorganic carbon DIC
   dic_varid = define_variable(ncid,"DIC", nf90_real, all_dimids, MITgcm_3D_FIELD)
   call add_attributes_to_variable(ncid, dic_varid, "Dissolved Inorganic Carbon", "mol C/m3", "moles Carbon per cubic meters")

      ! 7. BLING tracer: dissolved organic phosphorus DOP
   dop_varid = define_variable(ncid,"DOP", nf90_real, all_dimids, MITgcm_3D_FIELD)
   call add_attributes_to_variable(ncid, dop_varid, "Dissolved Organic Phosphorus", "mol P/m3", "moles Phosphorus per cubic meters")
   
   ! 8. BLING tracer: dissolved organic nitrogen DON
   don_varid = define_variable(ncid,"DON", nf90_real, all_dimids, MITgcm_3D_FIELD)
   call add_attributes_to_variable(ncid, don_varid, "Dissolved Organic Nitrogen", "mol N/m3", "moles Nitrogen per cubic meters")
   
   ! 9. BLING tracer: dissolved inorganic iron FET
   fet_varid = define_variable(ncid,"FET", nf90_real, all_dimids, MITgcm_3D_FIELD)
   call add_attributes_to_variable(ncid, fet_varid, "Dissolved Inorganic Iron", "mol Fe/m3", "moles Iron per cubic meters")
    
   ! 10. BLING tracer: Surface Chlorophyl CHL
   chl_varid = define_variable_2d(ncid,"CHL", nf90_real, all_dimids)
   call add_attributes_to_variable(ncid, chl_varid, "Surface Chlorophyll", "mg/m3", "milligram per cubic meters" )
endif

! Finished with dimension/variable definitions, must end 'define' mode to fill.

call check(nf90_enddef(ncid))

! Fill the coordinate variables

call check(nf90_put_var(ncid, XGVarID, XG ))
call check(nf90_put_var(ncid, XCVarID, XC ))
call check(nf90_put_var(ncid, YGVarID, YG ))
call check(nf90_put_var(ncid, YCVarID, YC ))
call check(nf90_put_var(ncid, ZCVarID, ZC ))

if (compress) then
   allocate(XCcomp(ncomp3))
   allocate(XGcomp(ncomp3))
   allocate(YCcomp(ncomp3))
   allocate(YGcomp(ncomp3))
   allocate(ZCcomp(ncomp3))
   allocate(ZGcomp(ncomp3))
   allocate(Xcomp_ind(ncomp3))
   allocate(Ycomp_ind(ncomp3))
   allocate(Zcomp_ind(ncomp3))

   allocate(Xcomp_indU(ncomp3U))
   allocate(Ycomp_indU(ncomp3U))
   allocate(Zcomp_indU(ncomp3U))

   allocate(Xcomp_indV(ncomp3V))
   allocate(Ycomp_indV(ncomp3V))
   allocate(Zcomp_indV(ncomp3V))

   call fill_compressed_coords()

   call check(nf90_put_var(ncid, XGcompVarID, XGcomp ))
   call check(nf90_put_var(ncid, XCcompVarID, XCcomp ))
   call check(nf90_put_var(ncid, YGcompVarID, YGcomp ))
   call check(nf90_put_var(ncid, YCcompVarID, YCcomp ))
   call check(nf90_put_var(ncid, ZCcompVarID, ZCcomp ))

   call check(nf90_put_var(ncid, XindID, Xcomp_ind ))
   call check(nf90_put_var(ncid, YindID, Ycomp_ind ))
   call check(nf90_put_var(ncid, ZindID, Zcomp_ind ))

   call check(nf90_put_var(ncid, XindUD, Xcomp_indU ))
   call check(nf90_put_var(ncid, YindUD, Ycomp_indU ))
   call check(nf90_put_var(ncid, ZindUD, Zcomp_indU ))

   call check(nf90_put_var(ncid, XindVD, Xcomp_indV ))
   call check(nf90_put_var(ncid, YindVD, Ycomp_indV ))
   call check(nf90_put_var(ncid, ZindVD, Zcomp_indV ))
endif

! Fill the netcdf variables
call from_mit_to_netcdf_3d('PSAL.data', ncid, SVarID, MITgcm_3D_FIELD)
call from_mit_to_netcdf_3d('PTMP.data', ncid, TVarID, MITgcm_3D_FIELD)
call from_mit_to_netcdf_3d('UVEL.data', ncid, UVarID, MITgcm_3D_FIELD_U)
call from_mit_to_netcdf_3d('VVEL.data', ncid, VVarID, MITgcm_3D_FIELD_V)
call from_mit_to_netcdf_2d('ETA.data' , ncid, EtaVarID)

print *, 'Done writing physical variables'

if (do_bgc) then
   call from_mit_to_netcdf_tracer_3d('NO3.data', ncid, no3_varid)
   call from_mit_to_netcdf_tracer_3d('PO4.data', ncid, po4_varid)
   call from_mit_to_netcdf_tracer_3d('O2.data' , ncid, o2_varid)
   call from_mit_to_netcdf_tracer_3d('PHY.data', ncid, phy_varid)
   call from_mit_to_netcdf_tracer_3d('ALK.data', ncid, alk_varid)
   call from_mit_to_netcdf_tracer_3d('DIC.data', ncid, dic_varid)
   call from_mit_to_netcdf_tracer_3d('DOP.data', ncid, dop_varid)
   call from_mit_to_netcdf_tracer_3d('DON.data', ncid, don_varid)
   call from_mit_to_netcdf_tracer_3d('FET.data', ncid, fet_varid)
   call from_mit_to_netcdf_tracer_2d('CHL.data', ncid, chl_varid)

   print *, 'Done writing biogeochemical variables'
endif

call check(nf90_close(ncid))

end subroutine mit2dart

!------------------------------------------------------------------
!> Subroutine for Reading netCDF and writing in binary

subroutine dart2mit()

integer :: ncid
recl2d = Nx*Ny*8

if (.not. module_initialized) call static_init_trans

call check(nf90_open("INPUT.nc",NF90_NOWRITE,ncid))

if (compress) then
   ncomp2 = nc_get_dimension_size(ncid,'comp2d')

   ncomp3  = nc_get_dimension_size(ncid,'comp3d')
   ncomp3U = nc_get_dimension_size(ncid,'comp3dU')
   ncomp3V = nc_get_dimension_size(ncid,'comp3dV')

   allocate(Xcomp_ind(ncomp3))
   allocate(Ycomp_ind(ncomp3))
   allocate(Zcomp_ind(ncomp3))
   
   allocate(Xcomp_indU(ncomp3U))
   allocate(Ycomp_indU(ncomp3U))
   allocate(Zcomp_indU(ncomp3U))

   allocate(Xcomp_indV(ncomp3V))
   allocate(Ycomp_indV(ncomp3V))
   allocate(Zcomp_indV(ncomp3V))
   
   call nc_get_variable(ncid, 'Xcomp_ind', Xcomp_ind)
   call nc_get_variable(ncid, 'Ycomp_ind', Ycomp_ind)
   call nc_get_variable(ncid, 'Zcomp_ind', Zcomp_ind)

   call nc_get_variable(ncid, 'Xcomp_indU', Xcomp_indU)
   call nc_get_variable(ncid, 'Ycomp_indU', Ycomp_indU)
   call nc_get_variable(ncid, 'Zcomp_indU', Zcomp_indU)

   call nc_get_variable(ncid, 'Xcomp_indV', Xcomp_indV)
   call nc_get_variable(ncid, 'Ycomp_indV', Ycomp_indV)
   call nc_get_variable(ncid, 'Zcomp_indV', Zcomp_indV)
endif

!Fill the data
iunit = get_unit()
open(iunit, file='PICKUP.OUTPUT', form="UNFORMATTED", status='UNKNOWN', &
            access='DIRECT', recl=recl2d, convert='BIG_ENDIAN')

call from_netcdf_to_mit_3d_pickup(ncid, 'UVEL', 1, MITgcm_3D_FIELD_U)
call from_netcdf_to_mit_3d_pickup(ncid, 'VVEL', 2, MITgcm_3D_FIELD_V)
call from_netcdf_to_mit_3d_pickup(ncid, 'PTMP', 3, MITgcm_3D_FIELD)
call from_netcdf_to_mit_3d_pickup(ncid, 'PSAL', 4, MITgcm_3D_FIELD)
call from_netcdf_to_mit_2d_pickup(ncid, 'ETA')

print *, 'Done writing physical variables into model binary files'

if (do_bgc) then
   call from_netcdf_to_mit_tracer_pickup(ncid, 'DIC', 1)
   call from_netcdf_to_mit_tracer_pickup(ncid, 'ALK', 2)
   call from_netcdf_to_mit_tracer_pickup(ncid, 'O2' , 3)
   call from_netcdf_to_mit_tracer_pickup(ncid, 'NO3', 4)
   call from_netcdf_to_mit_tracer_pickup(ncid, 'PO4', 5)
   call from_netcdf_to_mit_tracer_pickup(ncid, 'FET', 6)
   call from_netcdf_to_mit_tracer_pickup(ncid, 'DON', 7)
   call from_netcdf_to_mit_tracer_pickup(ncid, 'DOP', 8)
   call from_netcdf_to_mit_tracer_pickup(ncid, 'PHY', 9)
   print *, 'Done writing biogeochemical variables into model binary files'
endif

call check( NF90_CLOSE(ncid) )

if (compress) deallocate(Xcomp_ind, Ycomp_ind, Zcomp_ind)

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
! 3D variable
function define_variable(ncid, VARname, nc_type, all_dimids, field) result(varid)

integer, intent(in)          :: ncid
character(len=*), intent(in) :: VARname  ! variable name
integer, intent(in)          :: nc_type
integer, intent(in)          :: all_dimids(9) ! possible dimension ids
integer, intent(in)          :: field
integer                      :: varid ! netcdf variable id

integer ::  dimids(3)

if (compress) then
   if (field == MITgcm_3D_FIELD) then 
      call check(nf90_def_var(ncid=ncid, name=VARname, xtype=nc_type, &
           dimids=all_dimids(7),varid=varid))
   elseif (field == MITgcm_3D_FIELD_U) then
      call check(nf90_def_var(ncid=ncid, name=VARname, xtype=nc_type, &
           dimids=all_dimids(8),varid=varid))
   elseif (field == MITgcm_3D_FIELD_V) then
      call check(nf90_def_var(ncid=ncid, name=VARname, xtype=nc_type, &
           dimids=all_dimids(9),varid=varid))
   endif
else
   dimids = which_dims(VARname, all_dimids)
   call check(nf90_def_var(ncid=ncid, name=VARname, xtype=nc_type, &
        dimids=dimids, varid=varid))
endif

end function define_variable

!------------------------------------------------------------------
! For the non-compressed variables, X,Y,Z dimesnions vary
! depending on the variable
function which_dims(VARname, all_dimids) result(dimids)

character(len=*), intent(in) :: VARname  ! variable name
integer,          intent(in) :: all_dimids(9)
integer                      :: dimids(3)
! 3D variables, 3 grids:
! XC, YC, ZC  1  PSAL, PTMP, NO3, PO4, O2, PHY, ALK, DIC, DOP, DON, FET
! XG, YC, ZC  2  UVEL
! XC, YG, ZC  3  VVEL

if (VARname == 'UVEL') then
   dimids = (/all_dimids(4),all_dimids(2),all_dimids(3)/)
   return
endif
if (VARname == 'VVEL') then
   dimids = (/all_dimids(1),all_dimids(5),all_dimids(3)/)
   return
endif

dimids = (/all_dimids(1),all_dimids(2),all_dimids(3)/)

end function

!------------------------------------------------------------------
! 2D variable
function define_variable_2d(ncid, name, nc_type, all_dimids) result(varid)

integer, intent(in)          :: ncid
character(len=*), intent(in) :: name  ! variable name
integer, intent(in)          :: nc_type
integer, intent(in)          :: all_dimids(9)
integer                      :: varid ! netcdf variable id

! 2D variables, 1 grid:
! YC, XC 1 ETA, CHL

if (compress) then
   call check(nf90_def_var(ncid=ncid, name=name, xtype=nc_type, &
        dimids = (/all_dimids(6)/),varid=varid))
else
   call check(nf90_def_var(ncid=ncid, name=name, xtype=nc_type, &
        dimids = (/all_dimids(1),all_dimids(2)/),varid=varid))
endif

end function define_variable_2d

!------------------------------------------------------------------
subroutine add_attributes_to_variable(ncid, varid, long_name, units, units_long_name)

integer,          intent(in) :: ncid, varid ! which file, which variable
character(len=*), intent(in) :: long_name, units, units_long_name

call check(nf90_put_att(ncid, varid, "long_name"      , long_name))
call check(nf90_put_att(ncid, varid, "missing_value"  , FVAL))
call check(nf90_put_att(ncid, varid, "_FillValue"     , FVAL))
call check(nf90_put_att(ncid, varid, "units"          , units))
call check(nf90_put_att(ncid, varid, "units_long_name", units_long_name))

end subroutine

!------------------------------------------------------------------
subroutine from_mit_to_netcdf_3d(mitfile, ncid, varid, field)

character(len=*), intent(in) :: mitfile
integer,          intent(in) :: ncid, varid, field ! which file, which variable, grid type

integer  :: iunit
real(r4) :: var_data(Nx,Ny,Nz)

iunit = get_unit()
! HK are the mit files big endian by default?
open(iunit, file=mitfile, form='UNFORMATTED', status='OLD', &
            access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
read(iunit,rec=1) var_data
close(iunit)

where (var_data == binary_fill) var_data = FVAL !HK do we also need a check for nans here?

if (compress) then
  call write_compressed(ncid, varid, var_data, field)
else
  call check(nf90_put_var(ncid,varid,var_data))
endif

end subroutine from_mit_to_netcdf_3d

!------------------------------------------------------------------
subroutine from_mit_to_netcdf_2d(mitfile, ncid, varid)

character(len=*), intent(in) :: mitfile
integer,          intent(in) :: ncid, varid ! which file, which variable

integer  :: iunit
real(r4) :: var_data(Nx,Ny), var_T_data(Nx,Ny,Nz)

iunit = get_unit()
! HK are the mit files big endian by default?
open(iunit, file=mitfile, form='UNFORMATTED', status='OLD', &
            access='DIRECT', recl=recl2d, convert='BIG_ENDIAN')
read(iunit,rec=1) var_data
close(iunit)

! Manually get PTMP surface layer
open(iunit, file='PTMP.data', form='UNFORMATTED', status='OLD', &
            access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
read(iunit,rec=1) var_T_data
close(iunit)

where (var_T_data(:,:,1) == binary_fill) var_data = FVAL !HK do we also need a check for nans here?

if (compress) then
  call write_compressed(ncid, varid, var_data)
else
  call check(nf90_put_var(ncid,varid,var_data))
endif

end subroutine from_mit_to_netcdf_2d


!------------------------------------------------------------------
subroutine from_mit_to_netcdf_tracer_3d(mitfile, ncid, varid)

character(len=*), intent(in) :: mitfile
integer,          intent(in) :: ncid, varid ! which file, which variable

integer  :: iunit
real(r4) :: var_data(Nx,Ny,Nz)
real(r4) :: low_conc

low_conc = 1.0e-12

iunit = get_unit()
! HK are the mit files big endian by default?
open(iunit, file=mitfile, form='UNFORMATTED', status='OLD', &
            access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
read(iunit,rec=1) var_data
close(iunit)

! CHL is treated differently - HK CHL is 2d so you will not enter this
if (mitfile=='CHL.data') then
   where (var_data == binary_fill)
       var_data = FVAL
   elsewhere
       var_data = log10(var_data)
   endwhere
else
   ! Make sure the tracer concentration is positive
   where(var_data < binary_fill) var_data = low_conc
   
   if (log_transform) then
      where (var_data == binary_fill)
          var_data = FVAL
      elsewhere
          var_data = log(var_data)
      endwhere
   else
      where (var_data == binary_fill) var_data = FVAL
   endif
endif

if (compress) then
   call write_compressed(ncid, varid, var_data, MITgcm_3D_FIELD)
else
   call check(nf90_put_var(ncid,varid,var_data))
endif

end subroutine from_mit_to_netcdf_tracer_3d

!------------------------------------------------------------------
subroutine from_mit_to_netcdf_tracer_2d(mitfile, ncid, varid)

character(len=*), intent(in) :: mitfile
integer,          intent(in) :: ncid, varid ! which file, which variable

integer  :: iunit
real(r4) :: var_data(Nx,Ny)
real(r4) :: low_conc

low_conc = 1.0e-12

iunit = get_unit()
! HK are the mit files big endian by default?
open(iunit, file=mitfile, form='UNFORMATTED', status='OLD', &
            access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
read(iunit,rec=1) var_data
close(iunit)

! CHL is treated differently
if (mitfile=='CHL.data') then
   where (var_data == binary_fill)
       var_data = FVAL
   elsewhere
       var_data = log10(var_data)
   endwhere
else
   ! Make sure the tracer concentration is positive
   where(var_data < binary_fill) var_data = low_conc
   
   if (log_transform) then
      where (var_data == binary_fill)
          var_data = FVAL
      elsewhere
          var_data = log(var_data)
      endwhere
   else
      where (var_data == binary_fill) var_data = FVAL
   endif
endif

if (compress) then
   call write_compressed(ncid, varid, var_data)
else
   call check(nf90_put_var(ncid,varid,var_data))
endif

end subroutine from_mit_to_netcdf_tracer_2d

!------------------------------------------------------------------
subroutine from_netcdf_to_mit_2d(ncid, name)

integer,          intent(in) :: ncid ! which file,
character(len=*), intent(in) :: name ! which variable

integer  :: iunit
real(r4) :: var(Nx,Ny)
integer  :: varid
real(r4) :: local_fval

call check( NF90_INQ_VARID(ncid,name,varid) )
call check( nf90_get_att(ncid,varid,"_FillValue",local_fval))
! initialize var to netcdf fill value
var(:,:) = local_fval

if (compress) then
   call read_compressed(ncid, varid, var)
else
   call check(nf90_get_var(ncid,varid,var))
endif

where (var == local_fval) var = binary_fill

iunit = get_unit()
open(iunit, file=trim(name)//'.data', form="UNFORMATTED", status='UNKNOWN', &
            access='DIRECT', recl=recl2d, convert='BIG_ENDIAN')
write(iunit,rec=1)var
close(iunit)

end subroutine from_netcdf_to_mit_2d

!------------------------------------------------------------------
subroutine from_netcdf_to_mit_3d(ncid, name, field)

integer,          intent(in) :: ncid ! which file,
character(len=*), intent(in) :: name ! which variable

integer  :: iunit, field
real(r4) :: var(Nx,Ny,Nz)
integer  :: varid
real(r4) :: local_fval

call check( NF90_INQ_VARID(ncid,name,varid) )
call check( nf90_get_att(ncid,varid,"_FillValue",local_fval))
! initialize var to netcdf fill value
var(:,:,:) = local_fval

if (compress) then
   call read_compressed(ncid, varid, var, field)
else
   call check(nf90_get_var(ncid,varid,var))
endif

where (var == local_fval) var = binary_fill

iunit = get_unit()
open(iunit, file=trim(name)//'.data', form="UNFORMATTED", status='UNKNOWN', &
            access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
write(iunit,rec=1)var
close(iunit)

end subroutine from_netcdf_to_mit_3d

!------------------------------------------------------------------
subroutine from_netcdf_to_mit_2d_pickup(ncid, name)

integer,          intent(in) :: ncid ! which file,
character(len=*), intent(in) :: name ! which variable

integer  :: iunit
real(r4) :: var(Nx,Ny)
real(r8) :: var8(Nx,Ny)
integer  :: varid
real(r4) :: local_fval



call check( NF90_INQ_VARID(ncid,name,varid) )
call check( nf90_get_att(ncid,varid,"_FillValue",local_fval))

! initialize var to netcdf fill value
var(:,:) = local_fval

if (compress) then
   call read_compressed(ncid, varid, var)
else
   call check(nf90_get_var(ncid,varid,var))
endif

where (var == local_fval) var = binary_fill
var8 = var


if (do_bgc) then
  write(iunit,rec=401) var8
else
  write(iunit,rec=481) var8
endif
close(iunit)

end subroutine from_netcdf_to_mit_2d_pickup

!------------------------------------------------------------------
subroutine from_netcdf_to_mit_3d_pickup(ncid, name, lev, field)

integer,          intent(in) :: ncid ! which file,
character(len=*), intent(in) :: name ! which variable

integer  :: iunit, lev, field
real(r4) :: var(Nx,Ny,Nz)
real(r8) :: var8(Nx,Ny,Nz)
integer  :: varid, i
real(r4) :: local_fval
integer  :: LB, RB, RF
 
  
call check( NF90_INQ_VARID(ncid,name,varid) )
call check( nf90_get_att(ncid,varid,"_FillValue",local_fval))

! initialize var to netcdf fill value
var(:,:,:) = local_fval
   
if (compress) then
   call read_compressed(ncid, varid, var, field)
else
   call check(nf90_get_var(ncid,varid,var))
endif

where (var == local_fval) var = binary_fill
var8 = var



LB = Nz * (lev-1) + 1
RB = Nz * lev
RF = Nz * (lev-1)
do i = LB, RB
   write(iunit,rec=i) var8(:, :, i - RF)
enddo
close(iunit)

end subroutine from_netcdf_to_mit_3d_pickup

!------------------------------------------------------------------
subroutine from_netcdf_to_mit_tracer(ncid, name)

integer,          intent(in) :: ncid ! which file
character(len=*), intent(in) :: name ! which variable

integer  :: iunit
real(r4) :: var(Nx,Ny,Nz)
integer  :: varid
real(r4) :: local_fval

call check( NF90_INQ_VARID(ncid,name,varid) )
call check( nf90_get_att(ncid,varid,"_FillValue",local_fval))
! initialize var to netcdf fill value
var(:,:,:) = local_fval

if (compress) then
   call read_compressed(ncid, varid, var, MITgcm_3D_FIELD)
else
  call check(nf90_get_var(ncid,varid,var))
endif

if (log_transform) then
   where (var == local_fval)
       var = binary_fill
   elsewhere
       var = exp(var)
   endwhere
else
   where (var == local_fval) var = binary_fill
endif

iunit = get_unit()
open(iunit, file=trim(name)//'.data', form="UNFORMATTED", status='UNKNOWN', &
            access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
write(iunit,rec=1)var
close(iunit)

end subroutine from_netcdf_to_mit_tracer

!------------------------------------------------------------------
subroutine from_netcdf_to_mit_tracer_pickup(ncid, name, lev)

integer,          intent(in) :: ncid ! which file
character(len=*), intent(in) :: name ! which variable

integer  :: iunit, lev
real(r4) :: var(Nx,Ny,Nz)
real(r8) :: var8(Nx,Ny,Nz)
integer  :: varid
real(r4) :: local_fval
real(r4) :: low_conc, large_conc = 5.0 ! From Siva's old code

low_conc = 1.0e-12

call check( NF90_INQ_VARID(ncid,name,varid) )
call check( nf90_get_att(ncid,varid,"_FillValue",local_fval))

! initialize var to netcdf fill value
var(:,:,:) = local_fval

if (compress) then 
   call read_compressed(ncid, varid, var, MITgcm_3D_FIELD) 
else
  call check(nf90_get_var(ncid,varid,var))
endif

if (log_transform) then 
   where (var == local_fval)
       var = binary_fill
   elsewhere
       var = exp(var)
   endwhere
else
   where (var == local_fval) var = binary_fill
   where (var > large_conc) var = low_conc
endif

var8 = var

iunit = get_unit()
open(iunit, file='PICKUP_PTRACERS.OUTPUT', form="UNFORMATTED", status='UNKNOWN', &
            access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
write(iunit,rec=lev) var8
close(iunit)

end subroutine from_netcdf_to_mit_tracer_pickup

!------------------------------------------------------------------
subroutine from_netcdf_to_mit_tracer_chl(ncid, name)

integer,          intent(in) :: ncid ! which file
character(len=*), intent(in) :: name ! which variable

integer  :: iunit
real(r4) :: var(Nx,Ny)
integer  :: varid
real(r4) :: local_fval

call check( NF90_INQ_VARID(ncid,name,varid) )
call check( nf90_get_att(ncid,varid,"_FillValue",local_fval))
! initialize var to netcdf fill value
var(:,:) = local_fval

if (compress) then
   call read_compressed(ncid, varid, var)
else
  call check(nf90_get_var(ncid,varid,var))
endif

where (var == local_fval)
    var = binary_fill
elsewhere
    var = 10**(var)
endwhere


iunit = get_unit()
open(iunit, file=trim(name)//'.data', form="UNFORMATTED", status='UNKNOWN', &
            access='DIRECT', recl=recl2d, convert='BIG_ENDIAN')
write(iunit,rec=1)var
close(iunit)

end subroutine from_netcdf_to_mit_tracer_chl


!------------------------------------------------------------------
! Assumes all 3D variables are masked in the 
! same location
function get_compressed_size_3d(field) result(n3)

integer   :: n3, field
integer   :: iunit
real(r4)  :: var3d(NX,NY,NZ)
integer   :: i, j, k
character(len=MAX_LEN_FNAM) :: source

if (field == MITgcm_3D_FIELD)   source = 'PSAL.data'
if (field == MITgcm_3D_FIELD_U) source = 'UVEL.data'
if (field == MITgcm_3D_FIELD_V) source = 'VVEL.data' 

iunit = get_unit()
open(iunit, file=trim(source), form='UNFORMATTED', status='OLD', &
            access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
read(iunit,rec=1) var3d
close(iunit)

n3 = 0

! Get compressed size
do i=1,NX
   do j=1,NY
      do k=1,NZ
         if (var3d(i,j,k) /= binary_fill) then !HK also NaN?
           n3 = n3 + 1
         endif
      enddo
   enddo
enddo

end function get_compressed_size_3d

!------------------------------------------------------------------
! Assumes all 2D variables are masked in the
! same location
function get_compressed_size_2d() result(n2)

integer  :: n2
integer  :: iunit
real(r4) :: var3d(NX,NY,NZ)
integer  :: i,j

iunit = get_unit()
open(iunit, file='PTMP.data', form='UNFORMATTED', status='OLD', &
            access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
read(iunit,rec=1) var3d
close(iunit)

n2 = 0

! Get compressed size
do i=1,NX
   do j=1,NY
      if (var3d(i,j,1) /= binary_fill) then !HK also NaN?
        n2 = n2 + 1
      endif
   enddo
enddo

end function get_compressed_size_2d

!------------------------------------------------------------------
subroutine fill_compressed_coords()

!XG,etc read from PARAM04 in static_init_trans
real(r4) :: var3d(NX,NY,NZ)
integer :: n, i, j, k

iunit = get_unit()
open(iunit, file='PSAL.data', form='UNFORMATTED', status='OLD', &
            access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
read(iunit,rec=1) var3d
close(iunit)

n = 1

do k=1,NZ  ! k first so 2d is first
   do i=1,NX
       do j=1,NY
         if (var3d(i,j,k) /= binary_fill) then !HK also NaN?
            XCcomp(n) = XC(i)
            YCcomp(n) = YC(j)
            ZCcomp(n) = ZC(k)
            XGcomp(n) = XG(i)
            YGcomp(n) = YG(j)
            ZGcomp(n) = ZG(k)

            Xcomp_ind(n) = i  ! Assuming grids are compressed the same
            Ycomp_ind(n) = j
            Zcomp_ind(n) = k

            n = n + 1
         endif
      enddo
   enddo
enddo

! UVEL: 
iunit = get_unit()
open(iunit, file='UVEL.data', form='UNFORMATTED', status='OLD', &
            access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
read(iunit,rec=1) var3d
close(iunit)

n = 1
   
do k=1,NZ  ! k first so 2d is first
   do i=1,NX
       do j=1,NY
         if (var3d(i,j,k) /= binary_fill) then !HK also NaN?
            Xcomp_indU(n) = i  
            Ycomp_indU(n) = j
            Zcomp_indU(n) = k
   
            n = n + 1
         endif
      enddo
   enddo
enddo

! VVEL: 
iunit = get_unit()
open(iunit, file='VVEL.data', form='UNFORMATTED', status='OLD', &
            access='DIRECT', recl=recl3d, convert='BIG_ENDIAN')
read(iunit,rec=1) var3d
close(iunit)

n = 1

do k=1,NZ  ! k first so 2d is first
   do i=1,NX
       do j=1,NY
         if (var3d(i,j,k) /= binary_fill) then !HK also NaN?
            Xcomp_indV(n) = i
            Ycomp_indV(n) = j
            Zcomp_indV(n) = k
  
            n = n + 1
         endif
      enddo
   enddo
enddo

end subroutine fill_compressed_coords

!------------------------------------------------------------------
subroutine write_compressed_2d(ncid, varid, var_data)

integer,  intent(in) :: ncid, varid
real(r4), intent(in) :: var_data(Nx,Ny)

real(r4) :: comp_var(ncomp2)
integer  :: n
integer  :: i,j ! loop variables

n = 1
do i = 1, NX
   do j = 1, NY
      if (var_data(i,j) /= FVAL) then
         comp_var(n) = var_data(i,j)
         n = n + 1
      endif
   enddo
enddo

call check(nf90_put_var(ncid,varid,comp_var))

end subroutine write_compressed_2d

!------------------------------------------------------------------
subroutine write_compressed_3d(ncid, varid, var_data, field)

integer,  intent(in) :: ncid, varid, field
real(r4), intent(in) :: var_data(Nx,Ny,Nz)

real(r4), allocatable :: comp_var(:)
integer  :: n
integer  :: i,j,k ! loop variables

if (field == MITgcm_3D_FIELD_U) then 
  allocate(comp_var(ncomp3U))
  do i = 1,ncomp3U
     comp_var(i) = var_data(Xcomp_indU(i), Ycomp_indU(i), Zcomp_indU(i))
  enddo

elseif (field == MITgcm_3D_FIELD_V) then
  allocate(comp_var(ncomp3V))
  do i = 1,ncomp3V
     comp_var(i) = var_data(Xcomp_indV(i), Ycomp_indV(i), Zcomp_indV(i))
  enddo

else
  allocate(comp_var(ncomp3))
  do i = 1,ncomp3
     comp_var(i) = var_data(Xcomp_ind(i), Ycomp_ind(i), Zcomp_ind(i))
  enddo 
endif

!n = 1
!do k = 1, NZ !k first so 2d is first
!   do i = 1, NX
!      do j = 1, NY
!         if (var_data(i,j,k) /= FVAL) then
!            print *, 'n: ', n, ', var_data(i,j,k): ', var_data(i,j,k)
!            comp_var(n) = var_data(i,j,k)
!            n = n + 1
!         endif
!      enddo
!   enddo
!enddo

call check(nf90_put_var(ncid,varid,comp_var))

deallocate(comp_var)

end subroutine write_compressed_3d

!------------------------------------------------------------------
subroutine read_compressed_2d(ncid, varid, var)

integer,  intent(in)    :: ncid, varid
real(r4), intent(inout) :: var(NX,NY)

real(r4) :: comp_var(ncomp2)
integer  :: n ! loop variable
integer  :: i,j,k ! x,y,z
integer  :: c

c = 1

call check(nf90_get_var(ncid,varid,comp_var))

do n = 1, ncomp3 
   i = Xcomp_ind(n)
   j = Ycomp_ind(n)
   k = Zcomp_ind(n)
   if (k == 1) then
     var(i,j) = comp_var(c)
     c = c + 1
   endif
enddo

end subroutine read_compressed_2d

!------------------------------------------------------------------
subroutine read_compressed_3d(ncid, varid, var, field)

integer,  intent(in)    :: ncid, varid, field
real(r4), intent(inout) :: var(NX,NY,NZ)

real(r4), allocatable :: comp_var(:)
integer  :: n ! loop variable
integer  :: i,j,k ! x,y,k

if (field == MITgcm_3D_FIELD_U) then
  allocate(comp_var(ncomp3U))
  call check(nf90_get_var(ncid,varid,comp_var))
  do n = 1, ncomp3U
     i = Xcomp_indU(n)
     j = Ycomp_indU(n)
     k = Zcomp_indU(n)
     var(i,j,k) = comp_var(n)
  enddo

elseif (field == MITgcm_3D_FIELD_V) then
  allocate(comp_var(ncomp3V))
  call check(nf90_get_var(ncid,varid,comp_var))
  do n = 1, ncomp3V
     i = Xcomp_indV(n)
     j = Ycomp_indV(n)
     k = Zcomp_indV(n)
     var(i,j,k) = comp_var(n)
  enddo

else
  allocate(comp_var(ncomp3))
  call check(nf90_get_var(ncid,varid,comp_var))
  do n = 1, ncomp3
     i = Xcomp_ind(n)
     j = Ycomp_ind(n)
     k = Zcomp_ind(n)
     var(i,j,k) = comp_var(n)
  enddo
endif

deallocate(comp_var)

end subroutine read_compressed_3d

end module trans_mitdart_mod

