function obsstruct = plot_obs_netcdf_diffs(fname, ObsTypeString, region,  ...
                      CopyString1, CopyString2, QCString, maxQC, verbose)
%
% fname         = 'obs_sequence_001.nc';
% ObsTypeString = 'RADIOSONDE_U_WIND_COMPONENT';
% region        = [0 360 -90 90 -Inf Inf];
% CopyString1   = 'NCEP BUFR observation';
% CopyString2   = 'prior ensemble mean';
% QCString      = 'DART quality control';
% maxQC         = 1;
% verbose       = 1;   % anything > 0 == 'true'
%
% bob = plot_obs_netcdf_diffs(fname, ObsTypeString, region, CopyString1, CopyString2, ...
%                       QCString, maxQC, verbose);

% record the user input

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2007, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

obsstruct.fname         = fname;
obsstruct.ObsTypeString = ObsTypeString;
obsstruct.region        = region;
obsstruct.CopyString1   = CopyString1;
obsstruct.CopyString2   = CopyString2;
obsstruct.QCString      = QCString;
obsstruct.maxQC         = maxQC;
obsstruct.verbose       = verbose;

% get going

ObsTypes       = nc_varget(fname,'ObsTypes');
ObsTypeStrings = nc_varget(fname,'ObsTypesMetaData');
CopyStrings    = nc_varget(fname,'CopyMetaData');
QCStrings      = nc_varget(fname,'QCMetaData');

t              = nc_varget(fname,'time');
obs_type       = nc_varget(fname,'obs_type');
z_type         = nc_varget(fname,'which_vert');

loc            = nc_varget(fname,'location');
obs            = nc_varget(fname,'observations');
qc             = nc_varget(fname,'qc');

my_types   = unique(obs_type);  % only ones in the file, actually.
timeunits  = nc_attget(fname,'time','units');
timerange  = nc_attget(fname,'time','valid_range');
calendar   = nc_attget(fname,'time','calendar');
timebase   = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin = datenum(timebase(1),timebase(2),timebase(3));
timestring = datestr(timerange + timeorigin);

% Echo summary if requested

if ( verbose > 0 ) 
   for i = 1:length(my_types)
      obtype = my_types(i);
      inds   = find(obs_type == obtype);
      myz    = loc(inds,3);
     
      disp(sprintf('N = %6d %s obs (type %3d) between levels %.2f and %.2f', ...
               length(inds), ObsTypeStrings(obtype,:), obtype, ...
               unique(min(myz)), unique(max(myz))))
   end

%  uniquelevels = unique(loc(:,3));

%  for i = 1:length(uniquelevels)
%     mylevel = uniquelevels(i);
%     inds    = find(loc(:,3) == mylevel);
%     disp(sprintf('level %2d %f has %d observations',i,mylevel,length(inds)))
%  end

end

% Find observations of the correct types.

myind     = strmatch(ObsTypeString,ObsTypeStrings);

if ( isempty(myind) )           
   error('no %s observations ... stopping',obsstruct.ObsTypeString)
end

mytype1   = get_copy_index(fname, CopyString1);
mytype2   = get_copy_index(fname, CopyString2);
inds      = find(obs_type == myind);
mylocs    = loc(inds,:);
myobs1    = obs(inds,mytype1);
myobs2    = obs(inds,mytype2);
myobs     = myobs2 - myobs1;

if ~ isempty(QCString)
   myQCind = get_qc_index(fname,  QCString);
   myqc    = qc(inds,myQCind);
else
   myqc    = [];
end

clear myobs1 myobs2 obs loc qc

% geographic subset if needed

inds = locations_in_region(mylocs,region);

obsstruct.lons = mylocs(inds,1);
obsstruct.lats = mylocs(inds,2);
obsstruct.z    = mylocs(inds,3);
obsstruct.obs  =  myobs(inds);
obsstruct.Ztyp = z_type(inds);
obsstruct.numbadqc = 0;

if (isempty(myqc))
   obsstruct.qc = [];
else
   obsstruct.qc = myqc(inds);
end

% subset based on qc value

if ( (~ isempty(myqc)) & (~ isempty(maxQC)) )

   inds = find(obsstruct.qc > maxQC);

   obsstruct.numbadqc = length(inds);
   
   if (~isempty(inds))
       badobs.lons = obsstruct.lons(inds);
       badobs.lats = obsstruct.lats(inds);
       badobs.Ztyp = obsstruct.Ztyp(inds);
       badobs.z    = obsstruct.z(   inds);
       badobs.obs  = obsstruct.obs(inds);
       badobs.qc   = obsstruct.qc(inds);       
   end
   
   disp(sprintf('Removing %d obs with a %s value greater than %f', ...
                length(inds),QCString,maxQC))

   inds = find(obsstruct.qc <= maxQC);

   bob = obsstruct.lons(inds); obsstruct.lons = bob;
   bob = obsstruct.lats(inds); obsstruct.lats = bob;
   bob = obsstruct.Ztyp(inds); obsstruct.Ztyp = bob;
   bob = obsstruct.z(   inds); obsstruct.z    = bob;
   bob = obsstruct.obs( inds); obsstruct.obs  = bob;
   bob = obsstruct.qc(  inds); obsstruct.qc   = bob;

end

%-------------------------------------------------------------------------------
% Create graphic with area-weighted symbols for the good observations.
%-------------------------------------------------------------------------------

figure(1); clf

xmin = min(region(1:2));
xmax = max(region(1:2));
ymin = min(region(3:4));
ymax = max(region(3:4));
zmin = min(obsstruct.z);
zmax = max(obsstruct.z);

scalearray = scaleme(obsstruct.obs,36);
scalearray = 128 * ones(size(obsstruct.obs));

scatter3(obsstruct.lons, obsstruct.lats, obsstruct.z, ...
              scalearray, obsstruct.obs,'d','filled');

axis([xmin xmax ymin ymax zmin zmax])

str1 = sprintf('%s level (%.2f - %.2f)',ObsTypeString,zmin,zmax);
str2 = sprintf('%s - %s (%d locations)',CopyString2,CopyString1,length(obsstruct.obs));
str3 = sprintf('%s - %s',timestring(1,:),timestring(2,:));

title( {str1, str3, str2}, 'Interpreter','none','FontSize',16);
xlabel('longitude')
ylabel('latitude')

if     (obsstruct.Ztyp(1) == -2) % VERTISUNDEF     = -2
   zlabel('curious ... undefined')
elseif (obsstruct.Ztyp(1) == -1) % VERTISSURFACE   = -1
   zlabel('surface')
elseif (obsstruct.Ztyp(1) ==  1) % VERTISLEVEL     =  1
   zlabel('level')
elseif (obsstruct.Ztyp(1) ==  2) % VERTISPRESSURE  =  2
   set(gca,'ZDir','reverse')
   zlabel('pressure')
elseif (obsstruct.Ztyp(1) ==  3) % VERTISHEIGHT    =  3
   zlabel('height')
end

myworldmap;
set(gca,'CLim',[min(obsstruct.obs) max(obsstruct.obs)])
h = colorbar;
set(get(h,'YLabel'),'String',ObsTypeString,'Interpreter','none')

%-------------------------------------------------------------------------------
% Create graphic of spatial distribution of 'bad' observations & their QC value.
%-------------------------------------------------------------------------------

if (obsstruct.numbadqc > 0 )

   figure(2); clf
   
   subplot('position',[0.1 0.20 0.8 0.65])
   scalearray = 128 * ones(size(badobs.obs));
   
   zmin = min(badobs.z);
   zmax = max(badobs.z);
   
   scatter3(badobs.lons, badobs.lats, badobs.z, scalearray, badobs.qc,'filled')
   
   title( {str1, str3, 'Bad Observations'}, 'Interpreter','none','FontSize',16);
   xlabel('longitude')
   ylabel('latitude')
   
   if     (badobs.Ztyp(1) == -2) % VERTISUNDEF     = -2
      zlabel('curious ... undefined')
   elseif (badobs.Ztyp(1) == -1) % VERTISSURFACE   = -1
      zlabel('surface')
   elseif (badobs.Ztyp(1) ==  1) % VERTISLEVEL     =  1
      zlabel('level')
   elseif (badobs.Ztyp(1) ==  2) % VERTISPRESSURE  =  2
      set(gca,'ZDir','reverse')
      zlabel('pressure')
   elseif (badobs.Ztyp(1) ==  3) % VERTISHEIGHT    =  3
      zlabel('height')
   end
   
   axis([region(1) region(2) ymin ymax zmin zmax])
   
   myworldmap;
   set(gca,'CLim',[min(badobs.qc) max(badobs.qc)])
   h = colorbar;
   set(get(h,'YLabel'),'String',QCString,'Interpreter','none')
   
   subplot('position',[0.1 0.05 0.8 0.10])
   axis off
   
   qcvals  = unique(badobs.qc);
   qccount = zeros(size(qcvals));
   for i = 1:length(qcvals)
      qccount(i) = sum(badobs.qc == qcvals(i));
      s{i} = sprintf('%d obs with qc == %d',qccount(i),qcvals(i));
   end
   
   dy = 1.0/length(s);
   for i = 1:length(s)
      text(0.0, (i-1)*dy ,s{i})
   end

end

function h = myworldmap

%---------------------------------------------------------------------------
% GET THE ELEVATION DATA AND SET UP THE ASSOCIATED COORDINATE DATA
%---------------------------------------------------------------------------

load topo;                 % GET Matlab-native [180x360] ELEVATION DATASET
lats = [-89.5:89.5];       % CREATE LAT ARRAY FOR TOPO MATRIX
lons = [0.5:359.5];        % CREATE LON ARRAY FOR TOPO MATRIX
nlon = length(lons);
nlat = length(lats);

%---------------------------------------------------------------------------
% IF WE NEED TO SWAP HEMISPHERES, DO SO NOW.
% If we didn't explicitly tell it, make a guess.
%---------------------------------------------------------------------------

ax   = axis;

if (ax(1) < -2)
   lons = lons - 180.0;
   topo = [ topo(:,nlon/2+1:nlon) topo(:,1:nlon/2) ];
end

%---------------------------------------------------------------------------
% We need to determine the geographic subset of the elevation matrix.
%---------------------------------------------------------------------------

lon_ind1 = min(find(ax(1) <= lons));
lon_ind2 = min(find(ax(2) <= lons));
lat_ind1 = min(find(ax(3) <= lats));
lat_ind2 = min(find(ax(4) <= lats));

if (isempty(lon_ind1)) lon_ind1 = 1;    end;
if (isempty(lon_ind2)) lon_ind2 = nlon; end;
if (isempty(lat_ind1)) lat_ind1 = 1;    end;
if (isempty(lat_ind2)) lat_ind2 = nlat; end;

elev = topo(lat_ind1:lat_ind2,lon_ind1:lon_ind2);
x    = lons(lon_ind1:lon_ind2);
y    = lats(lat_ind1:lat_ind2);

%---------------------------------------------------------------------------
% Contour the "subset"
% There are differences between 6.5 and 7.0 that make changing the colors
% of the filled contours a real pain. Providing both solutions.
%---------------------------------------------------------------------------

orgholdstate = ishold;
hold on;

switch  get(gca,'ZDir')
   case 'reverse'
      zlevel = max(ax(5:6));
   otherwise
      zlevel = min(ax(5:6));
end

fcolor = [0.7 0.7 0.7];    % light grey

[c,h] = contourf(x,y,elev,[0.0 0.0],'k-');

new_level = 1000;

h_patch = get(h, 'Children');

for i = 1:numel(h_patch)
    y = get(h_patch(i), 'YData');
    s = size(y);
    set(h_patch(i), 'ZData', zlevel*ones(s),'FaceColor',fcolor);
end

if (orgholdstate == 0) hold off; end;


function s = scaleme(x,minsize)
% scaleme returns a uniformly scaled array the same size as the input
% array where the maximum is 10 times the minimum 
maxsize = 10*minsize;
minx    = min(x);
maxx    = max(x);
slope   = (maxsize-minsize)/(maxx-minx);
b       = minsize - slope*minx;

s = x*slope + b;

