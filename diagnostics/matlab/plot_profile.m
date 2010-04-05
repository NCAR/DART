function plotdat = plot_profile(fname,copystring)
%% plot_profile plots the vertical profile of the observation-space quantities for all possible levels, all possible variables.
% Part of the observation-space diagnostics routines.
%
% 'obs_diag' produces a netcdf file containing the diagnostics.
%
% USAGE: plotdat = plot_profile(fname,copystring);
%
% fname  :  netcdf file produced by 'obs_diag'
% copystring :  'copy' string == quantity of interest. These
%            can be any of the ones available in the netcdf 
%            file 'CopyMetaData' variable.
%            (ncdump -v CopyMetaData obs_diag_output.nc)
%
% EXAMPLE:
%
% fname = 'obs_diag_output.nc';   % netcdf file produced by 'obs_diag'
% copystring = 'totalspread';   % 'copy' string == quantity of interest
% plotdat = plot_profile(fname,copystring);


%% DART software - Copyright � 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL:
% https://proxy.subversion.ucar.edu/DAReS/DART/trunk/diagnostics/matlab/plot_profile.m $
% $Id$
% $Revision$
% $Date$

if (exist(fname,'file') ~= 2)
   error('file/fname <%s> does not exist',fname)
end

% Harvest plotting info/metadata from netcdf file.

plotdat.fname         = fname;
plotdat.copystring    = copystring;

plotdat.binseparation = nc_attget(fname, nc_global, 'bin_separation');
plotdat.binwidth      = nc_attget(fname, nc_global, 'bin_width');
time_to_skip          = nc_attget(fname, nc_global, 'time_to_skip');
plotdat.rat_cri       = nc_attget(fname, nc_global, 'rat_cri');
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
region_names          = nc_varget(fname,'region_names');
plotdat.region_names  = deblank(region_names);

% Coordinate between time types and dates

calendar              = nc_attget(fname,'time','calendar');
timeunits             = nc_attget(fname,'time','units');
timebase              = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin            = datenum(timebase(1),timebase(2),timebase(3));
timefloats            = zeros(size(time_to_skip));  % stupid int32 type conversion
timefloats(:)         = time_to_skip(:);
skip_seconds          = timefloats(4)*3600 + timefloats(5)*60 + timefloats(6);
iskip                 = timefloats(3) + skip_seconds/86400.0;

plotdat.bincenters    = plotdat.bincenters + timeorigin;
plotdat.binedges      = plotdat.binedges   + timeorigin;
plotdat.Nbins         = length(plotdat.bincenters);
plotdat.toff          = plotdat.binedges(1) + iskip;

plotdat.timespan      = sprintf('%s through %s', datestr(plotdat.toff), ...
                        datestr(max(plotdat.binedges(:))));

% set up a structure with all static plotting components

plotdat.xlabel    = {sprintf('%s',copystring), plotdat.timespan};
plotdat.linewidth = 2.0;

[plotdat.allvarnames, plotdat.allvardims] = get_varsNdims(fname);
[plotdat.varnames,    plotdat.vardims]    = FindVerticalVars(plotdat);

plotdat.nvars       = length(plotdat.varnames);

plotdat.copyindex   = get_copy_index(fname,copystring); 
plotdat.Npossindex  = get_copy_index(fname,'Nposs');
plotdat.Nusedindex  = get_copy_index(fname,'Nused');

%----------------------------------------------------------------------
% Loop around (copy-level-region) observation types
%----------------------------------------------------------------------

for ivar = 1:plotdat.nvars
    
   % create the variable names of interest.
    
   plotdat.myvarname = plotdat.varnames{ivar};
   plotdat.guessvar  = sprintf('%s_VPguess',plotdat.varnames{ivar});
   plotdat.analyvar  = sprintf('%s_VPanaly',plotdat.varnames{ivar});

   % remove any existing postscript file - will simply append each
   % level as another 'page' in the .ps file.
   
   psfname = sprintf('%s_%s_profile.ps',plotdat.varnames{ivar},plotdat.copystring);
   fprintf('Removing %s from the current directory.\n',psfname)
   system(sprintf('rm %s',psfname));

   % get appropriate vertical coordinate variable

   guessdims = nc_var_dims(  fname, plotdat.guessvar);
   analydims = nc_var_dims(  fname, plotdat.analyvar);
   varinfo   = nc_getvarinfo(fname, plotdat.analyvar);

   if ( findstr('surface',guessdims{2}) > 0 )
      fprintf('%s is a surface field.\n',plotdat.guessvar)
      fprintf('Cannot display a surface field this way.\n')
   elseif ( findstr('undef',guessdims{2}) > 0 )
      fprintf('%s has no vertical definition.\n',plotdat.guessvar)
      fprintf('Cannot display this field this way.\n')
   end

   [level_org level_units nlevels level_edges Yrange] = FindVerticalInfo(fname, plotdat.guessvar);
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
   [levels, indices]   = sort(plotdat.level_org);
   plotdat.level       = levels;
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
       bob = NaN*ones(varinfo.Size(1),varinfo.Size(2),2);
       ted = NaN*ones(varinfo.Size(1),varinfo.Size(2),2);
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
   nposs = sum(guess(plotdat.Npossindex,:,:));

   if ( sum(nposs(:)) < 1 )
      fprintf('No obs for %s...  skipping\n', plotdat.varnames{ivar})
      continue
   end

   plotdat.ges_copy   = guess(plotdat.copyindex,  :, :);
   plotdat.anl_copy   = analy(plotdat.copyindex,  :, :);
   plotdat.ges_Nposs  = guess(plotdat.Npossindex, :, :);
   plotdat.anl_Nposs  = analy(plotdat.Npossindex, :, :);
   plotdat.ges_Nused  = guess(plotdat.Nusedindex, :, :);
   plotdat.anl_Nused  = guess(plotdat.Nusedindex, :, :);
   plotdat.Xrange     = FindRange(plotdat);

   % plot by region

   clf; orient tall

   for iregion = 1:plotdat.nregions
      plotdat.region = iregion;  
      plotdat.myregion = deblank(plotdat.region_names(iregion,:));

      myplot(plotdat);
   end

   if (plotdat.nregions > 2)
      CenterAnnotation(plotdat.myvarname)
   end
   BottomAnnotation(fname)

   % create a postscript file
   print(gcf,'-dpsc','-append',psfname);

end

%----------------------------------------------------------------------
% 'Helper' functions
%----------------------------------------------------------------------

function myplot(plotdat)

   % Interlace the [ges,anl] to make a sawtooth plot.
   % By this point, the middle two dimensions are singletons.
   % The data must be sorted to match the order of the levels.
   cg = plotdat.ges_copy(:,:,plotdat.region); CG = cg(plotdat.indices);
   ca = plotdat.anl_copy(:,:,plotdat.region); CA = ca(plotdat.indices);

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
      other_guess = mean(CG(isfinite(CG))); 
      other_analy = mean(CA(isfinite(CA))); 
   else
      other_guess = NaN;
      other_analy = NaN;
   end

   str_other_pr = sprintf('%s pr=%.5g',plotdat.copystring,other_guess);
   str_other_po = sprintf('%s po=%.5g',plotdat.copystring,other_analy);

   % Plot 'xxx' on the bottom axis.
   % The observation count will use the axis on the top.
   % Ultimately, we want to suppress the 'auto' feature of the
   % axis labelling, so we manually set some values that normally
   % don't need to be set.
   
   % if more then 4 regions, this will not work (well) ... 
   if ( plotdat.nregions > 2 )
       ax1 = subplot(2,2,plotdat.region);
   else
       ax1 = subplot(1,plotdat.nregions,plotdat.region);
       axpos = get(ax1,'Position');
       axpos(4) = 0.925*axpos(4);
       set(ax1,'Position',axpos);
   end

   Stripes(plotdat.Xrange, plotdat.level_edges);
   set(ax1,'YDir', plotdat.YDir,'YTick',plotdat.level)
   ylabel(plotdat.level_units)

   %% draw the result of the experiment

   hold on;
   h1 = plot(CG,plotdat.level,'k+-',CA,plotdat.level,'k+:');

   set(h1,'LineWidth',plotdat.linewidth);
   h = legend(h1, str_other_pr, str_other_po);
   legend(h,'boxoff')
   set(h,'Interpreter','none')

   switch plotdat.copystring
      case 'bias'
         biasline = line([0 0],plotdat.Yrange,'Color','k','Parent',ax1);
         set(biasline,'LineWidth',2.0,'LineStyle','-.')
         plotdat.xlabel = {sprintf('bias (%s)',plotdat.biasconv), plotdat.timespan};
      otherwise
   end

   hold off;

   %% Create another axes to use for plotting the observation counts

   ax2 = axes('position',get(ax1,'Position'), ...
           'XAxisLocation','top', ...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','b','YColor','b',...
           'YLim',get(ax1,'YLim'), ...
           'YDir',get(ax1,'YDir'));

   h2 = line(nobs_poss,plotdat.level,'Color','b','Parent',ax2);
   h3 = line(nobs_used,plotdat.level,'Color','b','Parent',ax2);
   set(h2,'LineStyle','none','Marker','o');
   set(h3,'LineStyle','none','Marker','+');   

   % use same Y ticks
   set(ax2,'YTick',     get(ax1,'YTick'), ...
           'YTicklabel',get(ax1,'YTicklabel'));

   % use the same X ticks, but find the right label values
   [xticks, newticklabels] = matchingXticks(ax1,ax2);
   set(ax2,'XTick', xticks, 'XTicklabel', newticklabels)

   set(get(ax2,'Xlabel'),'String','# of obs (o=poss, +=used)')
   set(get(ax1,'Xlabel'),'String',plotdat.xlabel,'Interpreter','none')

   if (plotdat.nregions <=2 )
      title({plotdat.myvarname, plotdat.myregion},  ...
        'Interpreter', 'none', 'Fontsize', 12, 'FontWeight', 'bold')
   else
      title(plotdat.myregion, ...
        'Interpreter', 'none', 'Fontsize', 12, 'FontWeight', 'bold')
   end



function CenterAnnotation(main)
subplot('position',[0.48 0.48 0.04 0.04])
axis off
h = text(0.5,0.5,main);
set(h,'HorizontalAlignment','center', ...
      'VerticalAlignment','bottom', ...
      'Interpreter','none', ...
      'FontSize',12, ...
      'FontWeight','bold')



function BottomAnnotation(main)
% annotates the filename containing the data being plotted
subplot('position',[0.10 0.01 0.8 0.04])
axis off
if ( main(1) == '/' )   % must be a absolute pathname
   string1 = sprintf('data file: %s',main);
else
   mydir = pwd;
   string1 = sprintf('data file: %s/%s',mydir,main);
end
h = text(0.5, 0.5, string1);
set(h, 'Interpreter', 'none', 'FontSize', 8);
set(h, 'HorizontalAlignment','center');



function [y,ydims] = FindVerticalVars(x)
% Returns UNIQUE (i.e. base) vertical variable names
if ( ~(isfield(x,'allvarnames') && isfield(x,'allvardims')))
   error('Doh! no ''allvarnames'' and ''allvardims'' components')
end

j = 0;
basenames = struct([]);
basedims  = struct([]);

for i = 1:length(x.allvarnames)
   indx = findstr('time',x.allvardims{i});
   if (isempty(indx)) 
      j = j + 1;

      basenames{j} = ReturnBase(x.allvarnames{i});
      basedims{j}  = x.allvardims{i};
   end
end

[b,i,j] = unique(basenames);
y       = struct([]);
ydims   = struct([]);

for j = 1:length(i)
   fprintf('%2d is %s\n',j,basenames{j})
    y{j} = basenames{j};
ydims{j} = basedims{j};
end



function [level_org level_units nlevels level_edges Yrange] = FindVerticalInfo(fname,varname)
% Find the vertical dimension and harvest some info

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


function s = ReturnBase(string1)
inds = findstr('_guess',string1);
if (inds > 0 )
   s = string1(1:inds-1);
end

inds = findstr('_analy',string1);
if (inds > 0 )
   s = string1(1:inds-1);
end

inds = findstr('_VPguess',string1);
if (inds > 0 )
   s = string1(1:inds-1);
end

inds = findstr('_VPanaly',string1);
if (inds > 0 )
   s = string1(1:inds-1);
end



function x = FindRange(y)
% Trying to pick 'nice' limits for plotting.
% Completely ad hoc ... and not well posed.
%
% In this scope, y is bounded from below by 0.0
%
% If the numbers are very small ... 

bob  = [y.ges_copy(:); ...
        y.anl_copy(:)];
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
   switch y.copystring
   case {'bias'}
      if  ymax < 0
         Yrange = [ ymin 0.0 ];
      elseif  ymin > 0
         Yrange = [ 0.0 ymax ];
      end
   end

   x = sort([min([Yrange(1) 0.0]) Yrange(2)] ,'ascend');
end



function h = Stripes(x,edges)
%% plot the subtle background stripes that indicate the vertical
%  extent of what was averaged.
%
% FIXME:
% This really should be modified to add a percentage of the data
% range to provide space for the legend. Right now it is hardwired
% to assume that we are plotting hPa, on a 'reverse' axis. 

hold on;

% plot two little dots at the corners to make Matlab
% choose some plot limits. Given those nice limits and
% tick labels ... KEEP THEM. Later, make the dots invisible.

h = plot([min(x) max(x)],[min(edges) max(edges)]);
axlims    = axis;
axlims(4) = max(edges);
axlims(3) = -100;   % This gives extra space for the legend.
axis(axlims)

xc = [ axlims(1) axlims(2) axlims(2) axlims(1) axlims(1) ];

for i = 1:2:(length(edges)-1)
  yc = [ edges(i) edges(i) edges(i+1) edges(i+1) edges(i) ];
  hf = fill(xc,yc,[0.8 0.8 0.8],'EdgeColor','none');
  set(hf,'FaceAlpha',0.3,'EdgeAlpha',0.3)
  set(hf,'AlphaDataMapping','none','FaceVertexAlphaData',0.3)
end
set(gca,'XGrid','on')
hold off;
set(h,'Visible','off')



function [xticks newticklabels] = matchingXticks(ax1, ax2)
%% This takes the existing X ticks from ax1 (presumed nice)
% and determines the matching labels for ax2 so we can keep
% at least one of the axes looking nice.
   
Dlimits = get(ax1,'XLim');
DXticks = get(ax1,'XTick');
nXticks = length(DXticks);
xlimits = get(ax2,'XLim');

slope   = (xlimits(2) - xlimits(1))/(Dlimits(2) - Dlimits(1));
xtrcpt  = xlimits(2) -slope*Dlimits(2);

xticks        = slope*DXticks + xtrcpt;
newticklabels = num2str(round(10*xticks')/10);

