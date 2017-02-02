function geotiff2netcdf(fname,ofname)
% geotiff2netcdf converts a geotiff file to a netcdf file
%
% Just to be safe, the function aborts if the output file exists.
% It tells you as nicely as possible
%
% EXAMPLE:
%
% fname  = '/glade/proj3/image/romine/geo/n0r_201105250000.tif';
% ofname = 'bob.nc';
% geotiff2netcdf(fname,ofname)
%

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (exist(fname,'file') ~= 2)
   fprintf('\n')
   error('%s does not exist.',fname)
end

if (exist(ofname,'file') == 2)
   fprintf('\n')
   error('%s exists ... please remove and rerun or pick another name.',ofname)
end

%             Filename: '/glade/proj3/image/romine/geo/n0r_201105250000.tif'
%          FileModDate: '27-Jan-2012 15:01:40'
%             FileSize: 15622714
%               Format: 'tif'
%        FormatVersion: []
%               Height: 2600
%                Width: 6000
%             BitDepth: 8
%            ColorType: 'indexed'
%            ModelType: 'ModelTypeGeographic'
%                  PCS: ''
%           Projection: ''
%               MapSys: ''
%                 Zone: []
%         CTProjection: ''
%             ProjParm: []
%           ProjParmId: {}
%                  GCS: 'WGS 84'
%                Datum: 'World Geodetic System 1984'
%            Ellipsoid: 'WGS 84'
%            SemiMajor: 6378137
%            SemiMinor: 6.3568e+06
%                   PM: 'Greenwich'
%    PMLongToGreenwich: 0
%            UOMLength: ''
%    UOMLengthInMeters: 1
%             UOMAngle: 'degree'
%    UOMAngleInDegrees: 1
%            TiePoints: [1x1 struct]
%           PixelScale: [3x1 double]
%            RefMatrix: [3x2 double]
%          BoundingBox: [2x2 double]
%         CornerCoords: [1x1 struct]
%         GeoTIFFCodes: [1x1 struct]
%
% bob.TiePoints
%
%    ImagePoints: [1x1 struct]
%    WorldPoints: [1x1 struct]
%
% bob.TiePoints.ImagePoints
%
%    Row: 0.5000
%    Col: 0.5000
%
% bob.TiePoints.WorldPoints
%
%    X: -126.0050
%    Y: 50.0050
%
% bob.PixelScale           
%
%    0.0100
%    0.0100
%         0
%
% bob.RefMatrix 
%
%         0   -0.0100
%    0.0100         0
% -126.0100   50.0100
%
% bob.BoundingBox
%
% -126.0050   24.0050
%  -66.0050   50.0050
%
% bob.CornerCoords
%
%      X: [-126.0050 -66.0050 -66.0050 -126.0050]
%      Y: [50.0050 50.0050 24.0050 24.0050]
%    Row: [0.5000 0.5000 2.6005e+03 2.6005e+03]
%    Col: [0.5000 6.0005e+03 6.0005e+03 0.5000]
%    Lat: [50.0050 50.0050 24.0050 24.0050]
%    Lon: [-126.0050 -66.0050 -66.0050 -126.0050]
%
% bob.GeoTIFFCodes
%
%           Model: 2
%             PCS: 32767
%             GCS: 4326
%       UOMLength: 32767
%        UOMAngle: 9122
%           Datum: 6326
%              PM: 8901
%       Ellipsoid: 7030
%        ProjCode: 32767
%      Projection: 32767
%    CTProjection: 32767
%          MapSys: 32767
%      ProjParmId: []

% help netcdf
% NETCDF Summary of MATLAB NETCDF capabilities.
%    MATLAB provides low-level access to netCDF files via direct access to 
%    more than 40 functions in the netCDF library.  To use these MATLAB 
%    functions, you must be familiar with the netCDF C interface.  The 
%    "NetCDF C Interface Guide" for version 4.0.1 may be consulted at 
%    <http://www.unidata.ucar.edu/software/netcdf/old_docs/docs_4_0_1/>.
% 
%    In most cases, the syntax of the MATLAB function is similar to the 
%    syntax of the netCDF library function.  The functions are implemented 
%    as a package called "netcdf".  To use these functions, one needs to 
%    prefix the function name with package name "netcdf", i.e. 
% 
%       ncid = netcdf.open ( ncfile, mode );
% 
%    The following table lists all the netCDF library functions supported by 
%    the netCDF package.
% 
%       File Functions
%       --------------
%       abort            - Revert recent netCDF file definitions.
%       close            - Close netCDF file.
%       create           - Create new netCDF file.
%       endDef           - End netCDF file define mode.
%       inq              - Return information about netCDF file.
%       inqFormat        - Return netCDF file format.
%       inqLibVers       - Return netCDF library version information.
%       open             - Open netCDF file.
%       reDef            - Set netCDF file into define mode.
%       setDefaultFormat - Change default netCDF file format.
%       setFill          - Set netCDF fill mode.
%       sync             - Synchronize netCDF dataset to disk.  
%       
%       Dimension Functions
%       -------------------
%       defDim           - Create netCDF dimension.
%       inqDim           - Return netCDF dimension name and length.
%       inqDimID         - Return dimension ID.
%       inqUnlimDims     - Return unlimited dimensions visible in group.
%       renameDim        - Change name of netCDF dimension.
%       
%       Group Functions
%       ---------------
%       defGrp           - Create group.
%       inqNcid          - Return ID of named group.
%       inqGrps          - Return IDs of child groups.
%       inqVarIDs        - Return all variable IDs for group.
%       inqDimIDs        - Return all dimension IDs visible from group.
%       inqGrpName       - Return relative name of group.
%       inqGrpNameFull   - Return complete name of group.
%       inqGrpParent     - Find ID of parent group.
% 
%       Variable Functions
%       ------------------
%       defVar           - Create netCDF variable.
%       defVarChunking   - Set chunking layout.
%       defVarDeflate    - Set variable compression.
%       defVarFill       - Set fill parameters for variable.
%       defVarFletcher32 - Set checksum mode.
%       getVar           - Return data from netCDF variable.
%       inqVar           - Return information about variable.
%       inqVarChunking   - Return chunking layout for variable.
%       inqVarDeflate    - Return variable compression information.
%       inqVarFill       - Return fill value setting for variable.
%       inqVarFletcher32 - Return checksum settings.
%       inqVarID         - Return ID associated with variable name.
%       putVar           - Write data to netCDF variable.
%       renameVar        - Change name of netCDF variable.
%       
%       Attribute Functions
%       -------------------
%       copyAtt          - Copy attribute to new location.
%       delAtt           - Delete netCDF attribute.
%       getAtt           - Return netCDF attribute.
%       inqAtt           - Return information about netCDF attribute.
%       inqAttID         - Return ID of netCDF attribute.
%       inqAttName       - Return name of netCDF attribute.
%       putAtt           - Write netCDF attribute.
%       renameAtt        - Change name of attribute.
% 
%  
%    The following functions have no equivalents in the netCDF library.
% 
%       getConstantNames - Return list of constants known to netCDF library.
%       getConstant      - Return numeric value of named constant
%  
%    Please read the files netcdfcopyright.txt and mexnccopyright.txt for 
%    more information.

x     = geotiffread(fname);
bob   = geotiffinfo(fname);
ncid  = netcdf.create(ofname, 'NC_NOCLOBBER');

     HeightDimID = netcdf.defDim(ncid,'Height',        bob.Height);
      WidthDimID = netcdf.defDim(ncid,'Width',         bob.Width);
   HeightP1DimID = netcdf.defDim(ncid,'HeightP1',      bob.Height+1);
    WidthP1DimID = netcdf.defDim(ncid,'WidthP1',       bob.Width+1);
 NPixelRowsDimID = netcdf.defDim(ncid,'PixelRows',     size(bob.PixelScale,1));
NRefMatRowsDimID = netcdf.defDim(ncid,'RefMatrixRows', size(bob.RefMatrix,1));
NRefMatColsDimID = netcdf.defDim(ncid,'RefMatrixCols', size(bob.RefMatrix,2));
  NBboxRowsDimID = netcdf.defDim(ncid,'BBoxRows',      size(bob.BoundingBox,1));
  NBboxColsDimID = netcdf.defDim(ncid,'BBoxCols',      size(bob.BoundingBox,2));
   NCCorLenDimID = netcdf.defDim(ncid,'CCorLen',       length(bob.CornerCoords.X));

VarID = netcdf.getConstant('NC_GLOBAL');
netcdf.putAtt(ncid,VarID,'Filename',         bob.Filename);
netcdf.putAtt(ncid,VarID,'FileModDate',      bob.FileModDate);
netcdf.putAtt(ncid,VarID,'FileSize',         bob.FileSize);
netcdf.putAtt(ncid,VarID,'OrgFormat',        bob.Format);
netcdf.putAtt(ncid,VarID,'BitDepth',         bob.BitDepth);
netcdf.putAtt(ncid,VarID,'ColorType',        bob.ColorType);
netcdf.putAtt(ncid,VarID,'ModelType',        bob.ModelType);
netcdf.putAtt(ncid,VarID,'GCS',              bob.GCS);
netcdf.putAtt(ncid,VarID,'Datum',            bob.Datum);
netcdf.putAtt(ncid,VarID,'Ellipsoid',        bob.Ellipsoid);
netcdf.putAtt(ncid,VarID,'SemiMajor',        bob.SemiMajor);
netcdf.putAtt(ncid,VarID,'SemiMinor',        bob.SemiMinor);
netcdf.putAtt(ncid,VarID,'PM',               bob.PM);
netcdf.putAtt(ncid,VarID,'PMLongToGreenwich',bob.PMLongToGreenwich);
netcdf.putAtt(ncid,VarID,'UOMLengthInMeters',bob.UOMLengthInMeters);
netcdf.putAtt(ncid,VarID,'UOMAngle',         bob.UOMAngle);
netcdf.putAtt(ncid,VarID,'UOMAngleInDegrees',bob.UOMAngleInDegrees);
netcdf.putAtt(ncid,VarID,'TiePoints.ImagePoints.Row',bob.TiePoints.ImagePoints.Row);
netcdf.putAtt(ncid,VarID,'TiePoints.ImagePoints.Col',bob.TiePoints.ImagePoints.Col);
netcdf.putAtt(ncid,VarID,'TiePoints.WorldPoints.X',bob.TiePoints.WorldPoints.X);
netcdf.putAtt(ncid,VarID,'TiePoints.WorldPoints.Y',bob.TiePoints.WorldPoints.Y);

VarID    = netcdf.defVar(ncid,'datmat'     ,'NC_BYTE', [WidthDimID HeightDimID]);
VarIDps  = netcdf.defVar(ncid,'PixelScale' ,'NC_FLOAT',NPixelRowsDimID);
VarIDref = netcdf.defVar(ncid,'RefMatrix'  ,'NC_FLOAT',[NRefMatColsDimID NRefMatRowsDimID]);
VarIDbb  = netcdf.defVar(ncid,'BoundingBox','NC_FLOAT',[NBboxColsDimID NBboxRowsDimID]);

VarIDccX   = netcdf.defVar(ncid,'CornerCoords.X'  ,'NC_FLOAT',NCCorLenDimID);
VarIDccY   = netcdf.defVar(ncid,'CornerCoords.Y'  ,'NC_FLOAT',NCCorLenDimID);
VarIDccRow = netcdf.defVar(ncid,'CornerCoords.Row','NC_FLOAT',NCCorLenDimID);
VarIDccCol = netcdf.defVar(ncid,'CornerCoords.Col','NC_FLOAT',NCCorLenDimID);
VarIDccLat = netcdf.defVar(ncid,'CornerCoords.Lat','NC_FLOAT',NCCorLenDimID);
VarIDccLon = netcdf.defVar(ncid,'CornerCoords.Lon','NC_FLOAT',NCCorLenDimID);

VarIdLat   = netcdf.defVar(ncid,'Lat'      ,'NC_FLOAT',HeightDimID);
VarIdLon   = netcdf.defVar(ncid,'Lon'      ,'NC_FLOAT', WidthDimID);
VarIdLatE  = netcdf.defVar(ncid,'Lat_edges','NC_FLOAT',HeightP1DimID);
VarIdLonE  = netcdf.defVar(ncid,'Lon_edges','NC_FLOAT', WidthP1DimID);

netcdf.endDef(ncid)

% Create some useful coordinate arrays.
% Must decide if N-S or S-N, no convention.

if ( bob.TiePoints.WorldPoints.Y > min(bob.CornerCoords.Lat) ) 
   % Tie Point is the max lat ... we are going North-to-South
   LatEdges = bob.TiePoints.WorldPoints.Y:-bob.PixelScale(1):min(bob.CornerCoords.Lat);
else
   % Tie Point is the min lat ... we are going South-to-North
   LatEdges = bob.TiePoints.WorldPoints.Y: bob.PixelScale(1):max(bob.CornerCoords.Lat);
end 
LonEdges = bob.TiePoints.WorldPoints.X: bob.PixelScale(1):bob.CornerCoords.Lon(2);
latdiff  = diff(LatEdges);
londiff  = diff(LonEdges);
Lats     = LatEdges(1:bob.Height) + latdiff/2;
Lons     = LonEdges(1:bob.Width ) + londiff/2;

% fill the variables

netcdf.putVar(ncid,VarID     ,x')
netcdf.putVar(ncid,VarIDps   ,bob.PixelScale)
netcdf.putVar(ncid,VarIDref  ,bob.RefMatrix')
netcdf.putVar(ncid,VarIDbb   ,bob.BoundingBox')
netcdf.putVar(ncid,VarIDccX  ,bob.CornerCoords.X)
netcdf.putVar(ncid,VarIDccY  ,bob.CornerCoords.Y)
netcdf.putVar(ncid,VarIDccRow,bob.CornerCoords.Row)
netcdf.putVar(ncid,VarIDccCol,bob.CornerCoords.Col)
netcdf.putVar(ncid,VarIDccLat,bob.CornerCoords.Lat)
netcdf.putVar(ncid,VarIDccLon,bob.CornerCoords.Lon)
netcdf.putVar(ncid,VarIdLat  ,Lats     );
netcdf.putVar(ncid,VarIdLon  ,Lons     );
netcdf.putVar(ncid,VarIdLatE ,LatEdges );
netcdf.putVar(ncid,VarIdLonE ,LonEdges );

netcdf.sync(ncid) 
netcdf.close(ncid)

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
