function PlotCubedSpherePatches (pentagons_file, cs_file, fields, levels, out_form, obs_locs)
%% PlotCubedSpherePatches reads in arrays containing a 1D fields on single levels of CAM-SE output,
% and uses grid data from the pentagons_file to plot the field as a patch object.
% This is a collection of polygons (quadrilaterals in this case, with shared vertices),
% each having a color corresponding to the field value inside it.
%
% This version identifies any patch that straddles the dateline, 
% instead of checking whether the Center_lon of the patch is within
% some small (guessed) distance of the dateline.
%
% Example:
%
% pentagons_file = '/glade/u/home/raeder/Homme/ne30np4_091226_pentagons.nc'
% cs_file        = 'caminput.nc'  (if processing P{oste}rior files, extract 1 copy first)
% fields         = {'T','V'};    % Cell array of strings.  Put PS last?
% levels         = [20 30];
% out_form         = 'pdf';        FORMAT RES(dpi) SIZE  ppt_compat.  VIEWERS
%                                  none   Don't write out files.
%                                  png    300      '100' yes          gimp or gthumb
%                                  png    600      '220' yes          gimp or gthumb
%                                  pdf    300      '350' aliasing     gv
%                                  pdf    600      '550' aliasing     gv
%                                  eps    300      '750' aliasing     gv
% Optional obs locations structure:
% obs_locs(1) = struct('type',{'T'},'lon',280,'lat',85)
% (not used/needed yet: obs_locs.count(1) = 1)
% So far, only 1 'type' of obs can be displayed at a time.
% ...
%                                    
% PlotCubedSpherePatches (pentagons_file, cs_file, fields, levels, out_form, [obs_locs])

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Upgrade?; can 'print -append' only be used with [e]ps files?  It would be nice to put
%           all levels of a field in a file.

% Set symbol for observation location scatter plot markers.
if (nargin == 6)
   if     (obs_locs(1).type == 'T') 
      marker = '+'
   elseif (obs_locs(1).type == 'Q')
      marker = 'h'
   elseif (obs_locs(1).type == 'V')
      marker = '^'
   elseif (obs_locs(1).type == 'U')
      marker = '>'
   else
      marker = '*'
   end
end

if (exist(pentagons_file,'file') ~= 2), error('%s does not exist.',pentagons_file); end
if (exist(cs_file,'file')        ~= 2), error('%s does not exist.',cs_file); end

%----------------------------------
% Set up a log file for run-time messages
% strfind appears to find the location(s?) of '.' within the cs_file string.
cs_dots = strfind(cs_file,'.');
dots_size = length(cs_dots);
cs_root = cs_file(1:cs_dots(dots_size)-1);

lgfname = sprintf('%s.log',cs_root);
logfid = fopen(lgfname,'wt');
% sprintf(logfid,'%s\n',lgfname)

% Overestimate of the node longitudinal spacing, in degrees, not accounting for latitude.
coarse_grid = 1.5

%--------------------------------------------------------------------------------------
% Read in and package the patch vertex information from the pentagons file.

% File with nodes-to-patches mapping information
info_pent = ncinfo(pentagons_file);
ncid_pent = netcdf.open(pentagons_file,'NC_NOWRITE')

% Get dimensions of corner lon and lat arrays
if (info_pent.Dimensions(1).Name == 'grid_size')
   ncols    = info_pent.Dimensions(1).Length
   ncorners = info_pent.Dimensions(2).Length
elseif (info_pent.Dimensions(1).Name == 'grid_corners')
   ncols    = info_pent.Dimensions(2).Length
   ncorners = info_pent.Dimensions(1).Length
else
   error('Dimension 1 of %s is neither grid_size nor grid_corners.',pentagons_file);
end

% # of 'corners' in the pentagons file.  
% 5  for unrefined, 
% 8  for refined to 1/4 degree (newgulf_30_x4.g_scrip), 
% 12 for refined to 1/8 degree (old refined grid).
if (ncorners ~= 5 & ncorners ~= 8 & ncorners ~= 12) 
   error('%d ~= {5,8,12}, the correct number of "corners".',ncorners); end

%- - - - - - - - -
% Get info from the pent file about the lons.
lon_ind_pent = get_info_index(info_pent, 'fld','grid_corner_lon')
% lon_ind_pent = find(info_pent.Variables == 'grid_corner_lon') - 1 ;
[lon_name_pent, xtype_pent, lon_dimids_pent, numatts_pent] = netcdf.inqVar(ncid_pent,lon_ind_pent);

% Read the field
FillValue = has_att(ncid_pent, lon_ind_pent, numatts_pent, '_FillValue');
data = my_getVar(ncid_pent,lon_ind_pent,FillValue);
% Repack the 1D data array into the ncorners corners of each ncol.
corner_lons = reshape(data,ncorners,ncols) ;
clear data

%- - - - - - - - -
% Repeat for lats.
lat_ind_pent = get_info_index(info_pent, 'fld', 'grid_corner_lat') 
[lat_name_pent, xtype_pent, lat_dimids_pent, numatts_pent] = netcdf.inqVar(ncid_pent,lat_ind_pent);
FillValue = has_att(ncid_pent, lat_ind_pent, numatts_pent, '_FillValue');
data = my_getVar(ncid_pent,lat_ind_pent,FillValue);
corner_lats = reshape(data,ncorners,ncols) ;
clear data

% Some corners are redundant.
% A NaN at a (redundant) corner tells patch to leave the polygon open, and unfilled.
% Having first and last corners with equal locations allows patch to fill the polygon.

%- - - - - - - - -
% Also get centers(nodes) lons 
clon_ind_pent = get_info_index(info_pent, 'fld', 'grid_center_lon') 
[clon_name_pent, xtype_pent, clon_dimids_pent, numatts_pent] = netcdf.inqVar(ncid_pent,clon_ind_pent);
FillValue = has_att(ncid_pent, clon_ind_pent, numatts_pent, '_FillValue');
center_lons = my_getVar(ncid_pent,clon_ind_pent,FillValue);

%- - - - - - - - -
% Also get centers(nodes) lats
clat_ind_pent = get_info_index(info_pent, 'fld', 'grid_center_lat') 
[clat_name_pent, xtype_pent, clat_dimids_pent, numatts_pent] = netcdf.inqVar(ncid_pent,clat_ind_pent);
FillValue = has_att(ncid_pent, clat_ind_pent, numatts_pent, '_FillValue');
center_lats = my_getVar(ncid_pent,clat_ind_pent,FillValue);

%- - - - - - - - - -

pause on
for n = 1:ncols
    % Patches that span the Greenwich meridian need special treatment
    % so that those patches don't span the picture.
    if (abs(center_lats(n)) ~= 90. )
       % Calculate a boundary of the region in which patches must be fixed.
       % The longitudinal extent of the region, in degrees, depends on the latitude.
       lon_lim = 360. - abs(coarse_grid ./ cos(degtorad(center_lats(n)))) ;
    else
       lon_lim = 0. ;
    end
    lon_diff = max(corner_lons(:,n)) - min(corner_lons(:,n)) ;
    if (lon_diff > lon_lim) 
       % Put all the corners on the same side of the Greenwich meridian.
       corner_lons(find(corner_lons(:,n)>lon_lim),n) = 0.;
       fprintf(logfid,'\n   Corner_lons(:,%d) reduced to 0. ',n);
       % fprintf(logfid,'%f', find(corner_lons(:,n)==0.,n)
       find(corner_lons(:,n)==0.,n)
    end
end

%--------------------------------------------------------------------------------------
% Extract the desired fields from the cubed sphere file
% transform them to the patch objects, and write them to files.

% File with fields on it.
info_cs   = ncinfo(cs_file);
ncid_cs   = netcdf.open(cs_file,       'NC_NOWRITE');
dum = fprintf('info_cs variables %s %s %s %s %s %s %s \n',info_cs.Variables.Name);

% Statistics header in log file
dum = fprintf(logfid,'\n  Statistics of each level of each field \n');

% Get dimension of level 
nlevs_ind = get_info_index(info_cs,'dim','lev')
nlevs = info_cs.Dimensions(nlevs_ind).Length

right = 360. + 2*coarse_grid ;
axlims = [-coarse_grid right -90. 90.];

% Search for each desired field in list of fields on cs_file.
nfld_plot   = length(fields)
% kluge for comparison
% frame = 1;
frame = 0;
for ifld = 1:nfld_plot
    % Find the index in the list of available fields of the variable to be plotted.
    % There really should be only 1.
    field_char = char(fields(ifld));
    dum = fprintf('Field %2s \n',field_char);
    plot_fld_ind = get_info_index(info_cs,'fld',fields(ifld))

    % Get info from the cs file about the variable.
    [fldname_cs, xtype_cs, fld_dimids_cs, numatts_cs] = netcdf.inqVar(ncid_cs,plot_fld_ind);


    % Read the field
    FillValue = has_att(ncid_cs, plot_fld_ind, numatts_cs, '_FillValue');
    data = my_getVar(ncid_cs,plot_fld_ind,FillValue) ;

    % The order of the dimensions; 
    fld_dimids_cs
    % fld_dimids_cs lists the dimids of fld FROM RIGHT TO LEFT as printed by ncdump.
    % Dimensions in fld_dimids_cs are in fortran-order; most rapidly varying first.
    % The values still use 0-based indexing, which applies to dimids, 
    % but not the referencing of Dimensions (and Variables).
    % I want to know what dimension is in the 2nd position as listed by ncdump,
    % which is the 2nd from last number in fld_dimids_cs.
    info_cs.Dimensions.Name
    len_fld_dimids = length(fld_dimids_cs)
    % 2nd from the last.
    name_ind = fld_dimids_cs(len_fld_dimids - 1)
    % But Dimensions has fortran indexing; 1,...  so add 1 to the fld_dimids_cs val
    name = info_cs.Dimensions(name_ind+1).Name

    if (strcmp(field_char,'PS')) 
       data2D = ones(ncols,1);
       data2D(:,1) = data ;
    elseif (strcmp(name,'copy'))
       % ncdump -h caminput.nc :  float T(time, copy, ncol, lev) ;
       % Note the ' at the end of the next expression.
       data2D = reshape(data,nlevs,ncols)';
    elseif (strcmp(name,'lev'))
       % ncdump -h *.i.*;         double T(time[,lev], ncol) ;
       data2D = reshape(data,ncols,nlevs);
    else
       error('data2D not filled for %s, 2nd2last %s, last %s',field_char, ...
             info_cs.Dimensions(name_ind).Name,info_cs.Dimensions(len_fld_dimids-1).Name);
    end

% See diagnostics/matlab/plot_rmse_xxx_evolution.m for the template
% for this multi-level file creation.

% Add variable label(s) here
    if (strcmp(field_char,'PS'))
       lvl_max = 1
    else
       lvl_max = length(levels)
    end
    dum = fprintf(logfid,'Field %s \n', field_char);
    for lvl = 1:lvl_max
        if (strcmp(field_char,'PS'))
           field = data2D(:,1);
        else
           field = data2D(:,levels(lvl));
        end

        % Output file name; a file for each field.
        outfname = sprintf('%s_%s_l%i.%s',cs_root,field_char,levels(lvl),out_form)
 
        frame = frame + 1;
        figure(frame)
        % Personal kluge to make windows big and well-positioned for Mac (and Matlab window) 
        % on the left, extra screen on the right.
        % The units are pixels, and the origin is the lower left corner of the computer screen.
        set (frame,'Position',[1000., 800., 1500., 1000.]) 

        clear cdata 
        p = patch(corner_lons,corner_lats,field', ...
                  'FaceColor', 'flat', ...
                  'EdgeColor', 'none',   ...
                  'CData', field,  ...
                  'CDataMapping','scaled');
        hold on
        if (nargin == 6)
           s = scatter([obs_locs.lon],[obs_locs.lat],100,'white',marker)
        end
        % Modify the (default=jet) color scheme
        % my_colors = jet;
        % my_colors (1,:) = 1.;
        % colormap (my_colors);
        colorbar('SouthOutside')
        % Custom colorbar can be used for plotting subregions using:
        % set (gca,'Clim',[.00 .14]);
        %           Color [low high]
        %            limit
        % pause

% Add map overlay here.
        continents;
% Add state outlines
%         USoverlay;

% Add level label(s; pressure as well as model level?) here.
% Add data date? 
       maintitle = sprintf('%s  level %d',field_char, levels(lvl));
       subtitle = sprintf('ne30np4') ;
       % title(maintitle, ...
       title({maintitle, subtitle}, ...
             'Interpreter', 'none', 'Fontsize', 20, 'FontWeight', 'bold')
       StatAnnotation(field,logfid,levels(lvl))
        axis(axlims)
       BottomAnnotation(cs_file)


       
       % If it's the first, check for existence of this file.  Abort, or continue (without 'append'?).
       % -painters is the renderer needed to produce vector output.
       % -depsc is a MATLAB format for producing color encapsulated postscript files.
       % -r sets the resolution, dots/inch.  The default for rasters is 150 (too low) 
       %    and 'otherwise' is 864 (too high?)
       % -append adds pictures to an existing outfname.   BUT can't be used with eps.
       switch ( out_form )
       case ( 'pdf' )
          print(gcf,'-painters','-dpdf', '-r600',outfname);
       case ( 'eps' )
          print(gcf,'-painters','-depsc','-r300',outfname);
       case ( 'png' )
          % Enough resolution for ppt
          % print(gcf,'-painters','-dpng', '-r600',outfname);
          print(gcf,'-painters','-dpng', '-r300',outfname);
       otherwise
       end   
    end

% Add footer label(s) here.
end

fclose(logfid);

%----------------------------------------------------
function data = has_att(ncid,varid,numatts,attstring)

data = [];

for iatt = 1:numatts
   attid = iatt - 1;
   attname = netcdf.inqAttName(ncid, varid, attid);

   switch( attname )
      case (attstring)
         data = netcdf.getAtt(ncid, varid, attstring);
      otherwise
   end
end

%----------------------------------------------
function data = my_getVar(ncid,varid,FillValue)

data = netcdf.getVar(ncid, varid);
if ( ~ isempty(FillValue) )
   data(data == FillValue) = NaN;
end

%----------------------------------------------
function get_info_index = get_info_index(info,kind,field)

% field is a cell array (size 1) of character(s).

get_info_index = -999;

if (kind == 'fld') 
   dum = fprintf('get_info_index variables %s %s %s %s %s %s %s \n',info.Variables.Name);
   num = length(info.Variables) 
elseif (kind == 'dim')
   dum = fprintf('get_info_index dimensions %s %s %s %s %s %s %s \n',info.Dimensions.Name);
   num = length(info.Dimensions)
end

for i = 1:num
    if (kind == 'fld') 
       % msg = ['Does' field ' match ' info.Variables(i).Name]
       % info_Var_Name = info.Variables(i).Name
       % field
       same = strcmp(info.Variables(i).Name,char(field)) ;
    elseif (kind == 'dim')
       same = strcmp(info.Dimensions(i).Name,char(field)) ;
    end
    if (same == 1)
       get_info_index = i ;
       % NetCDF functions uses C indexing
       if (kind == 'fld') 
          get_info_index = get_info_index - 1 ;
       end
       break
    end
end

%------------------------------------------------------
function StatAnnotation(field,logfid,level)
% annotates the directory containing the data being plotted
%                  [left bottom width height]
% Scoot down to make room for colorbar 'SouthOutside'
% subplot('position',[0.48 0.05 0.04 0.04])
subplot('position',[0.48 0.05 0.04 0.04])
axis off
string1 = sprintf('level=%3i  max=%-g  min=%-g  mean=%-g  median=%-g  stdev=%-g ',      ...
                   level,    max(field),min(field),mean(field),median(field),std(field));
h = text(0.05, 0.5, string1);
set(h,'HorizontalAlignment','center', ...
      'VerticalAlignment','middle',...
      'Interpreter','none',...
      'FontSize',12)
fprintf (logfid, '%s \n',string1);

%------------------------------------------------------
function BottomAnnotation(main)
% annotates the directory containing the data being plotted
% Scoot down to make room for colorbar 'SouthOutside'
% subplot('position',[0.48 0.01 0.04 0.04])
subplot('position',[0.48 0.01 0.04 0.04])
axis off
fullname = which(main);   % Could be in MatlabPath
if( isempty(fullname) )
   if ( main(1) == '/' )  % must be a absolute pathname
      string1 = sprintf('data file: %s',main);
   else                   % must be a relative pathname
      mydir = pwd;
      string1 = sprintf('data file: %s/%s',mydir,main);
   end
else
   string1 = sprintf('data file: %s',fullname);
end

h = text(0.0, 0.5, string1);
set(h,'HorizontalAlignment','center', ...
      'VerticalAlignment','middle',...
      'Interpreter','none',...
      'FontSize',8)

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
