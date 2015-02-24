function plotdat = plot_bias_xxx_profile(fname, copy, varargin)
%% plot_bias_xxx_profile plots the vertical profile of the observation-space quantities for all possible levels, all possible variables.
% Part of the observation-space diagnostics routines.
%
% 'obs_diag' produces a netcdf file containing the diagnostics.
% obs_diag condenses the obs_seq.final information into summaries for a few specified
% regions - on a level-by-level basis.
%
% The number of observations possible reflects only those observations
% that have incoming QC values of interest. Any observation with a DART
% QC of 5 or 6 is not considered 'possible' for the purpose of this graphic.
%
% NOTE: if the observation was designated as a TRUSTED observation in the
%       obs_diag program, the observations that were rejected by the outlier
%       threshhold STILL PARTICIPATE in the calculation of the rmse, spread, etc.
%       The _values_ plotted by plot_profile reflect that. The number of observations
%       "used" becomes unclear. The number of observations used (designated by the
%       asterisk) is ALWAYS the number of observations successfully assimilated.
%       For TRUSTED observations, this is different than the number used to calculate
%       bias, rmse, spread, etc.
%
% USAGE: plotdat = plot_bias_xxx_profile(fname, copy);
%
% fname    :  netcdf file produced by 'obs_diag'
%
% copy     : string defining the metric of interest. 'rmse', 'spread', etc.
%            Possible values are available in the netcdf 'CopyMetaData' variable.
%            (ncdump -v CopyMetaData obs_diag_output.nc)
%
% obsname  : Optional. If present, The strings of each observation type to plot.
%            Each observation type will be plotted in a separate graphic.
%            Default is to plot all available observation types.
%
% OUTPUT: 'plotdat' is a structure containing what was plotted.
%         A .pdf of each graphic is created. Each .pdf has a name that 
%         reflects the variable, quantity, and region being plotted.
%
% EXAMPLE 1: All the observation types possible are plotted in separate figures.
%
% fname   = 'obs_diag_output.nc';
% copy    = 'totalspread';
% plotdat = plot_bias_xxx_profile(fname, copy);
%
% EXAMPLE 2: Just a single observation type.
%
% fname   = 'obs_diag_output.nc';
% copy    = 'totalspread';
% obsname = 'RADIOSONDE_U_WIND_COMPONENT';
% plotdat = plot_bias_xxx_profile(fname, copy, 'obsname', obsname);

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%%--------------------------------------------------------------------
% Decode,Parse,Check the input
%---------------------------------------------------------------------

default_obsname = 'none';
p = inputParser;

addRequired(p,'fname',@ischar);
addRequired(p,'copy',@ischar);
addParamValue(p,'obsname',default_obsname,@ischar);
parse(p, fname, copy, varargin{:});

% if you want to echo the input
% disp(['fname   : ', p.Results.fname])
% disp(['copy    : ', p.Results.copy])
% disp(['obsname : ', p.Results.obsname])

if ~isempty(fieldnames(p.Unmatched))
   disp('Extra inputs:')
   disp(p.Unmatched)
end

if (exist(fname,'file') ~= 2)
   error('file/fname <%s> does not exist',fname)
end

%%--------------------------------------------------------------------
% Harvest plotting info/metadata from netcdf file.
%---------------------------------------------------------------------

plotdat.fname         = fname;
plotdat.copystring    = copy;

plotdat.binseparation = nc_attget(fname, nc_global, 'bin_separation');
plotdat.binwidth      = nc_attget(fname, nc_global, 'bin_width');
time_to_skip          = nc_attget(fname, nc_global, 'time_to_skip');
plotdat.lonlim1       = nc_attget(fname, nc_global, 'lonlim1');
plotdat.lonlim2       = nc_attget(fname, nc_global, 'lonlim2');
plotdat.latlim1       = nc_attget(fname, nc_global, 'latlim1');
plotdat.latlim2       = nc_attget(fname, nc_global, 'latlim2');
plotdat.biasconv      = nc_attget(fname, nc_global, 'bias_convention');

plotdat.bincenters    = nc_varget(fname, 'time');
plotdat.binedges      = nc_varget(fname, 'time_bounds');
plotdat.mlevel        = nc_varget(fname, 'mlevel');
plotdat.plevel        = nc_varget(fname, 'plevel');
plotdat.plevel_edges  = nc_varget(fname, 'plevel_edges');
plotdat.hlevel        = nc_varget(fname, 'hlevel');
plotdat.hlevel_edges  = nc_varget(fname, 'hlevel_edges');

diminfo               = nc_getdiminfo(fname,'region');
plotdat.nregions      = diminfo.Length;
plotdat.region_names  = nc_varget(fname,'region_names');

% Matlab wants character matrices to be Nx1 instead of 1xN.

if (plotdat.nregions == 1 && (size(plotdat.region_names,2) == 1) )
   plotdat.region_names = deblank(plotdat.region_names');
end

% Coordinate between time types and dates

calendar              = nc_attget(fname,'time','calendar');
timeunits             = nc_attget(fname,'time','units');
timebase              = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin            = datenum(timebase(1),timebase(2),timebase(3));
timefloats            = zeros(size(time_to_skip));  % stupid int32 type conversion
timefloats(:)         = time_to_skip(:);
skip_seconds          = timefloats(4)*3600 + timefloats(5)*60 + timefloats(6);
iskip                 = timefloats(3) + skip_seconds/86400.0;

% Set up a structure to use for plotting

plotdat.bincenters    = plotdat.bincenters + timeorigin;
plotdat.binedges      = plotdat.binedges   + timeorigin;
plotdat.Nbins         = length(plotdat.bincenters);
plotdat.toff          = plotdat.binedges(1) + iskip;
plotdat.timespan      = sprintf('%s through %s', datestr(plotdat.toff), ...
                        datestr(max(plotdat.binedges(:))));
plotdat.xlabel        = sprintf('bias (%s) and %s',plotdat.biasconv,copy);

[plotdat.allvarnames, plotdat.allvardims] = get_varsNdims(fname);
[plotdat.varnames,    plotdat.vardims]    = FindVerticalVars(plotdat);

plotdat.nvars         = length(plotdat.varnames);
plotdat.copyindex     = get_copy_index(fname,copy);
plotdat.biasindex     = get_copy_index(fname,'bias');
plotdat.Npossindex    = get_copy_index(fname,'Nposs');
plotdat.Nusedindex    = get_copy_index(fname,'Nused');
plotdat.NQC5index     = get_copy_index(fname,'N_DARTqc_5');
plotdat.NQC6index     = get_copy_index(fname,'N_DARTqc_6');

figuredata = setfigure();

%%---------------------------------------------------------------------
% Loop around (copy-level-region) observation types
%----------------------------------------------------------------------

% Either use all the variables or just the one optionally specified.

if strcmp(p.Results.obsname,'none')
   varlist = 1:plotdat.nvars;
else
   varlist = find (strcmpi(p.Results.obsname,plotdat.varnames));
   if isempty(varlist)
      error('%s is not in the list of observations',p.Results.obsname)
   end
end

for ivar = varlist

   % create the variable names of interest.

   plotdat.myvarname = plotdat.varnames{ivar};
   plotdat.guessvar  = sprintf('%s_VPguess',plotdat.varnames{ivar});
   plotdat.analyvar  = sprintf('%s_VPanaly',plotdat.varnames{ivar});

   plotdat.trusted   = nc_read_att(fname, plotdat.guessvar, 'TRUSTED');
   if (isempty(plotdat.trusted)), plotdat.trusted = 'NO'; end

   % get appropriate vertical coordinate variable

   guessdims = nc_var_dims(  fname, plotdat.guessvar);
   analydims = nc_var_dims(  fname, plotdat.analyvar);
   varinfo   = nc_getvarinfo(fname, plotdat.analyvar);

   % this is a superfluous check ... FindVerticalVars already weeds out
   % variables only present on surface or undef because obs_diag
   % does not time-average statistics for these.

   if (~ isempty(strfind(guessdims{2},'surface')))
      fprintf('%s is a surface field.\n',plotdat.guessvar)
      fprintf('Cannot display a surface field this way.\n')
      continue
   elseif (~ isempty(strfind(guessdims{2},'undef')))
      fprintf('%s has no vertical definition.\n',plotdat.guessvar)
      fprintf('Cannot display this field this way.\n')
      continue
   end

   [level_org, level_units, nlevels, level_edges, Yrange] = FindVerticalInfo(fname, plotdat.guessvar);
   plotdat.level_org   = level_org;
   plotdat.level_units = level_units;
   plotdat.nlevels     = nlevels;
   plotdat.level_edges = level_edges;
   plotdat.Yrange      = Yrange;

   % Matlab likes strictly ASCENDING order for things to be plotted,
   % then you can impose the direction. The data is stored in the original
   % order, so the sort indices are saved to reorder the data.

   if (plotdat.level_org(1) > plotdat.level_org(plotdat.nlevels))
      plotdat.YDir = 'reverse';
   else
      plotdat.YDir = 'normal';
   end

   % Add error-checking for output from older versions of obs_diag.

   [levels, indices]   = sort(plotdat.level_org);
   plotdat.level       = unique(levels);
   if (length(plotdat.level) ~= length(levels))
      error('There is a duplicated value in the array specifying the levels - must change your input.nml and rerun obs_diag')
   end

   plotdat.indices     = indices;
   level_edges         = sort(plotdat.level_edges);
   plotdat.level_edges = level_edges;

   guess = nc_varget(fname, plotdat.guessvar);
   analy = nc_varget(fname, plotdat.analyvar);
   n = size(analy);

   % singleton dimensions are auto-squeezed - which is unfortunate.
   % We want these things to be 3D. [copy-level-region]
   % Sometimes there is one region, sometimes one level, ...
   % To complicate matters, the stupid 'ones' function does not allow
   % the last dimension to be unity ... so you have double the size
   % of the array ...

   if ( plotdat.nregions == 1 )
      bob = NaN*ones(varinfo.Size(1),varinfo.Size(2),1);
      ted = NaN*ones(varinfo.Size(1),varinfo.Size(2),1);
      bob(:,:,1) = guess;
      ted(:,:,1) = analy;
      guess = bob; clear bob
      analy = ted; clear ted
   elseif ( plotdat.nlevels == 1 )
      bob = NaN*ones(varinfo.Size);
      ted = NaN*ones(varinfo.Size);
      bob(:,1,:) = guess;
      ted(:,1,:) = analy;
      guess = bob; clear bob
      analy = ted; clear ted
   end

   % check to see if there is anything to plot
   % The number possible is decreased by the number of observations
   % rejected by namelist control.

   fprintf('%d %s observations had DART QC of 5 (all regions).\n', ...
           sum(sum(guess(plotdat.NQC5index, :,:))),plotdat.myvarname)
   fprintf('%d %s observations had DART QC of 6 (all regions).\n', ...
           sum(sum(guess(plotdat.NQC6index, :,:))),plotdat.myvarname)

   nposs = sum(guess(plotdat.Npossindex,:,:)) - ...
           sum(guess(plotdat.NQC5index ,:,:)) - ...
           sum(guess(plotdat.NQC6index ,:,:));

   if ( sum(nposs(:)) < 1 )
      fprintf('No obs for %s...  skipping\n', plotdat.varnames{ivar})
      continue
   end

   plotdat.ges_copy   = guess(plotdat.copyindex,  :, :);
   plotdat.anl_copy   = analy(plotdat.copyindex,  :, :);
   plotdat.ges_bias   = guess(plotdat.biasindex,  :, :);
   plotdat.anl_bias   = analy(plotdat.biasindex,  :, :);
   plotdat.ges_Nqc5   = guess(plotdat.NQC5index,  :, :);
   plotdat.anl_Nqc5   = analy(plotdat.NQC5index,  :, :);
   plotdat.ges_Nqc6   = guess(plotdat.NQC6index,  :, :);
   plotdat.anl_Nqc6   = analy(plotdat.NQC6index,  :, :);
   plotdat.ges_Nused  = guess(plotdat.Nusedindex, :, :);
   plotdat.anl_Nused  = guess(plotdat.Nusedindex, :, :);
   plotdat.ges_Nposs  = guess(plotdat.Npossindex, :, :) - ...
                        plotdat.ges_Nqc5 - plotdat.ges_Nqc6;
   plotdat.anl_Nposs  = analy(plotdat.Npossindex, :, :) - ...
                        plotdat.anl_Nqc5 - plotdat.anl_Nqc6;
   plotdat.Xrange     = FindRange(plotdat);

   % plot by region - each in its own figure.

   for iregion = 1:plotdat.nregions
      figure(iregion); clf(iregion); orient(figuredata.orientation); wysiwyg
      plotdat.region   = iregion;
      plotdat.myregion = deblank(plotdat.region_names(iregion,:));
      myplot(plotdat, figuredata);
      BottomAnnotation(fname)

      psfname = sprintf('%s_bias_%s_profile_region%d', ...
                plotdat.varnames{ivar}, plotdat.copystring, iregion);
      print(gcf,'-dpdf',psfname);
   end

end


%=====================================================================
% 'Helper' functions
%=====================================================================


function myplot(plotdat,figdata)

%% Interlace the [ges,anl] to make a sawtooth plot.
% By this point, the middle two dimensions are singletons.
% The data must be sorted to match the order of the levels.
cg = plotdat.ges_copy(:,:,plotdat.region); CG = cg(plotdat.indices);
ca = plotdat.anl_copy(:,:,plotdat.region); CA = ca(plotdat.indices);

mg = plotdat.ges_bias(:,:,plotdat.region); MG = mg(plotdat.indices);
ma = plotdat.anl_bias(:,:,plotdat.region); MA = ma(plotdat.indices);

g = plotdat.ges_Nposs(:,:,plotdat.region); G = g(plotdat.indices);
a = plotdat.anl_Nposs(:,:,plotdat.region); A = a(plotdat.indices);

nobs_poss   = G;
nposs_delta = G - A;

g = plotdat.ges_Nused(:,:,plotdat.region); G = g(plotdat.indices);
a = plotdat.anl_Nused(:,:,plotdat.region); A = a(plotdat.indices);
nobs_used   = G;
nused_delta = G - A;

% Determine some quantities for the legend
nobs = sum(nobs_used);
if ( nobs > 1 )
   bias_guess  = mean(MG(isfinite(MG)));
   bias_analy  = mean(MA(isfinite(MA)));
   other_guess = mean(CG(isfinite(CG)));
   other_analy = mean(CA(isfinite(CA)));
else
   bias_guess  = NaN;
   bias_analy  = NaN;
   other_guess = NaN;
   other_analy = NaN;
end

str_bias_pr  = sprintf('%s pr=%.5g','bias',bias_guess);
str_bias_po  = sprintf('%s po=%.5g','bias',bias_analy);
str_other_pr = sprintf('%s pr=%.5g',plotdat.copystring,other_guess);
str_other_po = sprintf('%s po=%.5g',plotdat.copystring,other_analy);

% Plot the bias and 'xxx' on the same (bottom) axis.
% The observation count will use the axis on the top.
% Ultimately, we want to suppress the 'auto' feature of the
% axis labelling, so we manually set some values that normally
% don't need to be set.

ax1 = subplot('position',figdata.position);

% add type of vertical coordinate info for adjusting axes to accomodate legend

Stripes(plotdat.Xrange, plotdat.level_edges, plotdat.level_units);
set(ax1, 'YDir', plotdat.YDir, 'YTick', plotdat.level, 'Layer', 'top')
set(ax1,'YAxisLocation','left','FontSize',figdata.fontsize)

% draw the result of the experiment

hold on;
h1 = plot(MG,plotdat.level,'k+-',MA,plotdat.level,'k+--', ...
          CG,plotdat.level,'ro-',CA,plotdat.level,'ro--');
set(h1,'LineWidth',figdata.linewidth);
hold off;

zeroline = line([0 0],plotdat.Yrange,'Color',[0 100 0]/255,'Parent',ax1);
set(zeroline,'LineWidth',2.5,'LineStyle','-')

h = legend(h1, str_bias_pr, str_bias_po, str_other_pr, str_other_po, 'Location', 'NorthWest');
set(h,'Interpreter','none','Box','off')

% If the observation is trusted, reference that somehow
switch lower(plotdat.trusted)
   case 'true'
      axlims = axis;
      tx = axlims(2) + (axlims(2) - axlims(1))/20;
      if  strcmpi('normal',plotdat.YDir)
         ty = plotdat.Yrange(1);
      else 
         ty = plotdat.Yrange(2);
      end
      h = text(tx,ty,'TRUSTED. Values include outlying observations.');
      set(h,'FontSize',20,'Rotation',90,'VerticalAlignment','middle')
   otherwise
end

% Create another axes to use for plotting the observation counts

ax2 = axes('position',get(ax1,'Position'), ...
   'XAxisLocation','top', ...
   'YAxisLocation','right', ...
   'Color','none', ...
   'XColor','b', ...
   'YColor',get(ax1,'YColor'), ...
   'YLim',get(ax1,'YLim'), ...
   'YDir',get(ax1,'YDir'), ...
   'FontSize',get(ax1,'FontSize'));

h2 = line(nobs_poss,plotdat.level,'Color','b','Parent',ax2);
h3 = line(nobs_used,plotdat.level,'Color','b','Parent',ax2);
set(h2,'LineStyle','none','Marker','o');
set(h3,'LineStyle','none','Marker','*');

% use same Y ticks - but no labels.
set(ax2,'YTick',get(ax1,'YTick'), 'YTicklabel',[]);

% use the same X ticks, but find the right label values
xscale = matchingXticks(ax1,ax2);

set(get(ax1,'Ylabel'),'String',plotdat.level_units, ...
                      'Interpreter','none','FontSize',figdata.fontsize)
set(get(ax1,'Xlabel'),'String',{plotdat.xlabel, plotdat.timespan}, ...
                      'Interpreter','none','FontSize',figdata.fontsize)
set(get(ax2,'Xlabel'),'String', ...
    ['# of obs (o=possible, \ast=assimilated) x' int2str(uint32(xscale))],'FontSize',figdata.fontsize)

title({plotdat.myregion, plotdat.myvarname},  ...
      'Interpreter', 'none', 'FontSize', figdata.fontsize, 'FontWeight', 'bold')


%=====================================================================


function BottomAnnotation(main)
%% annotates the filename containing the data being plotted
subplot('position',[0.10 0.01 0.8 0.04])
axis off

if ( main(1) == '/' )   % must be an absolute pathname
   string1 = sprintf('data file: %s',main);
else
   mydir = pwd;
   string1 = sprintf('data file: %s/%s',mydir,main);
end

h = text(0.5, 0.33, string1);
set(h, 'HorizontalAlignment', 'center', ...
   'VerticalAlignment','middle', ...
   'Interpreter', 'none', ...
   'FontSize', 8);


%=====================================================================


function [y,ydims] = FindVerticalVars(x)
%% Returns UNIQUE (i.e. base) vertical variable names
% In this context, if the variable has a 'time' dimension
% it cannot be a variable of interest.

if ( ~(isfield(x,'allvarnames') && isfield(x,'allvardims')))
   error('Doh! no ''allvarnames'' and ''allvardims'' components')
end

j = 0;
basenames = struct([]);
basedims  = struct([]);

for i = 1:length(x.allvarnames)
   dimnames = lower(x.allvardims{i});
   if (isempty(strfind(dimnames,'time')))
      j = j + 1;

      basenames{j} = ReturnBase(x.allvarnames{i});
      basedims{j}  = x.allvardims{i};
   end
end

[~,i,j] = unique(basenames);
y       = struct([]);
ydims   = struct([]);

for k = 1:length(i)
   fprintf('%2d is %s\n',k,basenames{i(k)})
   y{k} = basenames{i(k)};
   ydims{k} = basedims{i(k)};
end


%=====================================================================


function [level_org, level_units, nlevels, level_edges, Yrange] = FindVerticalInfo(fname,varname)
%% Find the vertical dimension and harvest some info

varinfo  = nc_getvarinfo(fname,varname);
leveldim = [];

for i = 1:length(varinfo.Dimension)
   inds = strfind(varinfo.Dimension{i},'level');
   if ( ~ isempty(inds)), leveldim = i; end
end

if ( isempty(leveldim) )
   error('There is no level information for %s in %s',varname,fname)
end

level_org   = nc_varget(fname,varinfo.Dimension{leveldim});
level_units = nc_attget(fname,varinfo.Dimension{leveldim},'units');
nlevels     = varinfo.Size(leveldim);
edgename    = sprintf('%s_edges',varinfo.Dimension{leveldim});
level_edges = nc_varget(fname, edgename);
Yrange      = [min(level_edges) max(level_edges)];


%=====================================================================


function s = ReturnBase(string1)
%% Pick off the variable name.
s = [];
inds = strfind(string1,'_guess');
if (inds > 0 )
   s = string1(1:inds-1);
   return
end

inds = strfind(string1,'_analy');
if (inds > 0 )
   s = string1(1:inds-1);
   return
end

inds = strfind(string1,'_VPguess');
if (inds > 0 )
   s = string1(1:inds-1);
   return
end

inds = strfind(string1,'_VPanaly');
if (inds > 0 )
   s = string1(1:inds-1);
   return
end


%=====================================================================


function x = FindRange(y)
%% Trying to pick 'nice' limits for plotting.
% Completely ad hoc ... and not well posed.
%
% In this scope, y is bounded from below by 0.0
%
% If the numbers are very small ...

bob  = [y.ges_copy(:) ; y.ges_bias(:); y.anl_copy(:) ; y.anl_bias(:)];
inds = find(isfinite(bob));

if ( isempty(inds) )
   x = [0 1];
else
   glommed = bob(inds);
   ymin    = min(glommed);
   ymax    = max(glommed);

   if ( ymax > 1.0 )
      ymin = floor(min(glommed));
      ymax =  ceil(max(glommed));
   end

   if (ymin == 0 && ymax == 0)
      ymax = 1;
   end

   if (ymin == ymax)
      ymin = ymin - 0.1*ymin;
      ymax = ymax + 0.1*ymax;
   end

   Yrange = [ymin ymax];

   % Make sure a zero bias is visible on plot
   if  ymax < 0
      Yrange = [ ymin 0.0 ];
   elseif  ymin > 0
      Yrange = [ 0.0 ymax ];
   end

   x = sort([min([Yrange(1) 0.0]) Yrange(2)] ,'ascend');
end


%=====================================================================


function h = Stripes(x,edges,units)
%% plot the subtle background stripes that indicate the vertical
%  extent of what was averaged.
%
% FIXME:
% This really should be modified to add a percentage of the data
% range to provide space for the legend. Right now it is hardwired
% to assume that we are plotting hPa, on a 'reverse' axis.
% kdr axlims(3) should be conditional on the observation vertical coordinate:
%     values for pressure coordinates are inappropriate for height coord.
%     It also assumes 4 plots/page, but 2 works better for plotting all levels of CAM5.
%     That requires a smaller % of vertical range for the legend.

% plot two little dots at the corners to make Matlab
% choose some plot limits. Given those nice limits and
% tick labels ... KEEP THEM. Later, make the dots invisible.

h = plot([min(x) max(x)],[min(edges) max(edges)]);
axlims          = axis;
legend_fraction = 0.22;

% partial fix to legend space; add in option for vert coord = height.

switch lower(units)
   case 'hpa'
      axlims(4) = max(edges);
      axlims(3) = min(edges) - legend_fraction*(axlims(4)-min(edges));
   case 'm'
      axlims(3) = min(edges);
      axlims(4) = max(edges) + legend_fraction*(max(edges)-axlims(3));
   otherwise
end
axis(axlims)

% set up list of x,y values defining corner of every other stripe.

xc = [ axlims(1) axlims(2) axlims(2) axlims(1) axlims(1) ];

hold on;
for i = 1:2:(length(edges)-1)
   yc = [ edges(i) edges(i) edges(i+1) edges(i+1) edges(i) ];
   hf = fill(xc,yc,[0.8 0.8 0.8],'EdgeColor','none');
end
hold off;

set(gca,'XGrid','on')
set(h,'Visible','off') % make the dots invisible


%=====================================================================


function figdata = setfigure()
%%
%  figure out a page layout
%  extra space at the bottom for the date/file annotation
%  extra space at the top because the titles have multiple lines

orientation = 'tall';
fontsize    = 16;
position    = [0.15 0.12 0.7 0.75];
linewidth   = 2.0;

figdata = struct('expcolors',  {{'k','r','b','m','g','c','y'}}, ...
   'expsymbols', {{'o','s','d','p','h','s','*'}}, ...
   'prpolines',  {{'-','--'}}, 'position', position, ...
   'fontsize',fontsize, 'orientation',orientation, ...
   'linewidth',linewidth);


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

