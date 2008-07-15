function plotdat = plot_rmse_xxx_profile(fname,copystring)
% plot_rmse_xxx_profile plots the vertical profile of the observation-space RMSE and any other quantity for all possible levels, all possible variables.
% Part of the observation-space diagnostics routines.
%
% 'obs_diag' produces a netcdf file containing the diagnostics.
%
% USAGE: plotdat = plot_rmse_xxx_profile(fname,copystring);
%
% fname  :  netcdf file produced by 'obs_diag'
% copystring :  'copy' string == quantity of interest. These
%            can be any of the ones available in the netcdf 
%            file 'CopyMetaData' variable.
%            (ncdump -v CopyMetaData obs_diag_output.nc)
%
% EXAMPLE: plot the rmse and totalspread on the same axis:
%
% fname = 'obs_diag_output.nc';   % netcdf file produced by 'obs_diag'
% copystring = 'totalspread';   % 'copy' string == quantity of interest
% plotdat = plot_rmse_xxx_profile(fname,copystring);


% Data Assimilation Research Testbed -- DART
% Copyright 2004-2007, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL:
% https://proxy.subversion.ucar.edu/DAReS/DART/trunk/diagnostics/matlab/plot_rmse_xxx_profile.m $
% $Id$
% $Revision$
% $Date$

%----------------------------------------------------------------------
% Get plotting metadata from obs_diag run.
%----------------------------------------------------------------------

% Harvest info from netcdf file.

plotdat.fname         = fname;
plotdat.copystring    = copystring;
plotdat.bincenters    = getnc(fname,'time');
plotdat.binedges      = getnc(fname,'time_bounds');
plotdat.mlevel        = getnc(fname,'mlevel');
plotdat.plevel        = getnc(fname,'plevel');
plotdat.plevel_edges  = getnc(fname,'plevel_edges');
plotdat.hlevel        = getnc(fname,'hlevel');
plotdat.hlevel_edges  = getnc(fname,'hlevel_edges');
plotdat.nregions      = length(getnc(fname,'region'));
plotdat.region_names  = getnc(fname,'region_names');

if (plotdat.nregions == 1)
   plotdat.region_names = deblank(plotdat.region_names');
end

f = netcdf(fname,'nowrite');
plotdat.binseparation      = f.bin_separation(:);
plotdat.binwidth           = f.bin_width(:);
time_to_skip               = f.time_to_skip(:);
plotdat.rat_cri            = f.rat_cri(:);
plotdat.input_qc_threshold = f.input_qc_threshold(:);
plotdat.lonlim1            = f.lonlim1(:);
plotdat.lonlim2            = f.lonlim2(:);
plotdat.latlim1            = f.latlim1(:);
plotdat.latlim2            = f.latlim2(:);
plotdat.rmseconv           = f.rmse_convention(:);

% Coordinate between time types and dates

timeunits    = f{'time'}.units(:);
calendar     = f{'time'}.calendar(:);
timebase     = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin   = datenum(timebase(1),timebase(2),timebase(3));
skip_seconds = time_to_skip(4)*3600 + time_to_skip(5)*60 + time_to_skip(6);
iskip        = time_to_skip(3) + skip_seconds/86400;

plotdat.bincenters = plotdat.bincenters + timeorigin;
plotdat.binedges   = plotdat.binedges   + timeorigin;
plotdat.Nbins      = length(plotdat.bincenters);
plotdat.toff       = plotdat.bincenters(1) + iskip;

% set up a structure with all static plotting components

plotdat.xlabel    = sprintf('rmse and %s',copystring);
plotdat.linewidth = 2.0;

[plotdat.allvarnames, plotdat.allvardims] = get_varsNdims(f);
[plotdat.varnames,    plotdat.vardims]    = FindVerticalVars(plotdat);

plotdat.nvars = length(plotdat.varnames);

plotdat.copyindex   = get_copy_index(fname,copystring); 
plotdat.rmseindex   = get_copy_index(fname,'rmse');
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
   
   psfname = sprintf('%s_rmse_%s_profile.ps',plotdat.varnames{ivar},plotdat.copystring);
   disp(sprintf('Removing %s from the current directory.',psfname))
   system(sprintf('rm %s',psfname));

   % get appropriate vertical coordinate variable

   guessdims = nc_var_dims(f, plotdat.guessvar);
   analydims = nc_var_dims(f, plotdat.analyvar);

   if ( findstr('surface',guessdims{2}) > 0 )
      disp(sprintf('%s is a surface field.',plotdat.guessvar))
      error('Cannot display a surface field this way.')
   elseif ( findstr('undef',guessdims{2}) > 0 )
      disp(sprintf('%s has no vertical definition.',plotdat.guessvar))
      error('Cannot display this field this way.')
   else
      plotdat.level_org     = getnc(fname, guessdims{2});
      plotdat.level_units   = f{guessdims{2}}.units(:);
      plotdat.nlevels       = length(f{guessdims{2}});
      edgename              = sprintf('%s_edges',guessdims{2});
      plotdat.level_edges   = getnc(fname,edgename);
      plotdat.Yrange        = [min(plotdat.level_edges) max(plotdat.level_edges)];
   end

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
   
   guess = getnc(fname, plotdat.guessvar,-1,-1,-1,-1,-1,-1,0);  
   analy = getnc(fname, plotdat.analyvar,-1,-1,-1,-1,-1,-1,0); 
   
   % Check for one region ... if the last dimension is a singleton 
   % dimension, it is auto-squeezed  - which is bad.
   % We want these things to be 3D.

   n = size(guess);
   if ( length(n) < 3 )
       bob = NaN*ones(n(1),n(2),2);
       ted = NaN*ones(n(1),n(2),2);
       bob(:,:,1) = guess;
       ted(:,:,1) = analy;
       guess = bob; clear bob
       analy = ted; clear ted
   end
   
   % check to see if there is anything to plot
   nposs = sum(guess(plotdat.Npossindex,:,:));

   if ( sum(nposs(:)) < 1 )
      disp(sprintf('No obs for %s...  skipping', plotdat.varnames{ivar}))
      continue
   end

   plotdat.ges_copy   = guess(plotdat.copyindex,  :, :);
   plotdat.anl_copy   = analy(plotdat.copyindex,  :, :);
   plotdat.ges_rmse   = guess(plotdat.rmseindex,  :, :);
   plotdat.anl_rmse   = analy(plotdat.rmseindex,  :, :);

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
   % BottomAnnotation(ges)

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

   mg = plotdat.ges_rmse(:,:,plotdat.region); MG = mg(plotdat.indices);
   ma = plotdat.anl_rmse(:,:,plotdat.region); MA = ma(plotdat.indices);

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
      rmse_guess  = mean(MG(isfinite(MG)));   
      rmse_analy  = mean(MA(isfinite(MA)));   
      other_guess = mean(CG(isfinite(CG))); 
      other_analy = mean(CA(isfinite(CA))); 
   else
      rmse_guess  = NaN;
      rmse_analy  = NaN;
      other_guess = NaN;
      other_analy = NaN;
   end

   str_rmse_pr  = sprintf('%s pr=%.3f','rmse',rmse_guess);
   str_rmse_po  = sprintf('%s po=%.3f','rmse',rmse_analy);
   str_other_pr = sprintf('%s pr=%.3f',plotdat.copystring,other_guess);
   str_other_po = sprintf('%s po=%.3f',plotdat.copystring,other_analy);

   % Plot the rmse and 'xxx' on the same (bottom) axis.
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
   hold on;
   h1 = plot(MG,plotdat.level,'k+-',MA,plotdat.level,'k+:', ...
             CG,plotdat.level,'ro-',CA,plotdat.level,'ro:');
   hold off;
   set(h1,'LineWidth',plotdat.linewidth);
   h = legend(h1, str_rmse_pr, str_rmse_po, str_other_pr, str_other_po, ...
          'Location','East');
   legend(h,'boxoff')

   axlims = [plotdat.Xrange plotdat.Yrange];
   axis(axlims)
   set(gca,'YDir', plotdat.YDir)
   hold on; plot([0 0],plotdat.Yrange,'k-')

%   disp(plotdat.myregion)
%   [plotdat.level CG' CA' MG' MA' ma' mg' ca' cg' plotdat.level_org]

   set(gca,'YTick',plotdat.level,'Ylim',plotdat.Yrange)
   ylabel(plotdat.level_units)
   
   % use same X,Y limits for all plots in this region
   nXticks = length(get(ax1,'XTick'));
   xlimits = plotdat.Xrange;
   xinc    = (xlimits(2)-xlimits(1))/(nXticks-1);
   xticks  = xlimits(1):xinc:xlimits(2);
   set(ax1,'XTick',xticks,'Xlim',xlimits)
   
   % create a separate scale for the number of observations
   ax2 = axes('position',get(ax1,'Position'), ...
           'XAxisLocation','top', ...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','b','YColor','b',...
           'YLim',plotdat.Yrange, ...
           'YDir',plotdat.YDir);
   h2 = line(nobs_poss,plotdat.level,'Color','b','Parent',ax2);
   h3 = line(nobs_used,plotdat.level,'Color','b','Parent',ax2);
   set(h2,'LineStyle','none','Marker','o');
   set(h3,'LineStyle','none','Marker','+');   

   % use same number of X ticks and the same Y ticks
  
   xlimits = get(ax2,'XLim');
   xinc   = (xlimits(2)-xlimits(1))/(nXticks-1);
   xticks = xlimits(1):xinc:xlimits(2);
   nicexticks = round(10*xticks')/10;
   set(ax2,'YTick',get(ax1,'YTick'),'YTicklabel',[], ...
           'XTick',          xticks,'XTicklabel',num2str(nicexticks))
       
   set(get(ax2,'Xlabel'),'String','# of obs (dashed) o=poss, +=used')
   set(get(ax1,'Xlabel'),'String',plotdat.xlabel)
   set(ax1,'Position',get(ax2,'Position'))
   grid

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
% annotates the directory containing the data being plotted
subplot('position',[0.48 0.01 0.04 0.04])
axis off
bob = which(main);
[pathstr,name,ext,versn] = fileparts(bob);
h = text(0.0,0.5,pathstr);
set(h,'HorizontalAlignment','center', ...
      'VerticalAlignment','middle',...
      'Interpreter','none',...
      'FontSize',8)



function [y,ydims] = FindVerticalVars(x)
% Returns UNIQUE (i.e. base) vertical variable names
if ( ~(isfield(x,'allvarnames') && isfield(x,'allvardims')))
   error('Doh! no ''allvarnames'' and ''allvardims'' components')
end

j = 0;

for i = 1:length(x.allvarnames)
   indx = findstr('time',x.allvardims{i});
   if (isempty(indx)) 
      j = j + 1;

      basenames{j} = ReturnBase(x.allvarnames{i});
      basedims{j}  = x.allvardims{i};
   end
end

[b,i,j] = unique(basenames);

for j = 1:length(i)
   disp(sprintf('%2d is %s',j,basenames{j}))
    y{j} = basenames{j};
ydims{j} = basedims{j};
end



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
% In this scope, y is bounded from below by 0.0

glommed = [y.ges_copy(:); y.ges_rmse(:); ...
           y.anl_copy(:); y.anl_rmse(:)];
      
Yrange = [floor(min(glommed)) ceil(max(glommed))];

if ( isfinite(Yrange) )
   x = [min([Yrange(1) 0.0]) Yrange(2)];
else
   x = [0 1];
end



function Stripes(x,edges)
% EraseMode: [ {normal} | background | xor | none ]

xc = [ x(1) x(2) x(2) x(1) x(1) ];

hold on;
for i = 1:2:length(edges)
  yc = [ edges(i) edges(i) edges(i+1) edges(i+1) edges(i) ];
  h = fill(xc,yc,[0.8 0.8 0.8], ...
  'EraseMode','background','EdgeColor','none');
end
hold off;
