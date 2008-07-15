function plotdat = plot_evolution(fname,copystring)
% plot_evolution plots the temporal evolution of the observation-space quantities for all possible levels, all possible variables.
% Part of the observation-space diagnostics routines.
%
% 'obs_diag' produces a netcdf file containing the diagnostics.
%
% USAGE: plotdat = plot_evolution(fname,copystring);
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
% copystring = 'bias';   % 'copy' string == quantity of interest
% plotdat = plot_evolution(fname,copystring);

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
plotdat.biasconv           = f.bias_convention(:);

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

plotdat.linewidth = 2.0;

[plotdat.allvarnames, plotdat.allvardims] = get_varsNdims(f);
[plotdat.varnames,    plotdat.vardims]    = FindTemporalVars(plotdat);

plotdat.nvars = length(plotdat.varnames);

plotdat.copyindex   = get_copy_index(fname,copystring); 
plotdat.Npossindex  = get_copy_index(fname,'Nposs');
plotdat.Nusedindex  = get_copy_index(fname,'Nused');

%----------------------------------------------------------------------
% Loop around (time-copy-level-region) observation types
%----------------------------------------------------------------------

for ivar = 1:plotdat.nvars
    
   % create the variable names of interest.
    
   plotdat.myvarname = plotdat.varnames{ivar};  
   plotdat.guessvar  = sprintf('%s_guess',plotdat.varnames{ivar});
   plotdat.analyvar  = sprintf('%s_analy',plotdat.varnames{ivar});

   % remove any existing postscript file - will simply append each
   % level as another 'page' in the .ps file.
   
   psfname = sprintf('%s_%s_evolution.ps',plotdat.varnames{ivar},plotdat.copystring);
   disp(sprintf('Removing %s from the current directory.',psfname))
   system(sprintf('rm %s',psfname));

   % get appropriate vertical coordinate variable

   guessdims = nc_var_dims(f, plotdat.guessvar);
   analydims = nc_var_dims(f, plotdat.analyvar);

   if ( findstr('surface',guessdims{3}) > 0 )
      plotdat.level = 1;
      plotdat.level_units = 'surface';
   elseif ( findstr('undef',guessdims{3}) > 0 )
      plotdat.level = 1;
      plotdat.level_units = 'undefined';
   else
      plotdat.level = getnc(fname, guessdims{3});
      plotdat.level_units = f{guessdims{3}}.units(:);
   end

   guess = getnc(fname, plotdat.guessvar,-1,-1,-1,-1,-1,-1,0);  
   analy = getnc(fname, plotdat.analyvar,-1,-1,-1,-1,-1,-1,0); 
   
   % check to see if there is anything to plot
   nposs = sum(guess(:,plotdat.Npossindex,:,:));

   if ( sum(nposs(:)) < 1 )
      disp(sprintf('%s no obs for %s...  skipping', plotdat.varnames{ivar}))
      continue
   end
   
   % here is where you have to check for one region ... the last
   % dimension (if it is a singleton dimension) is auto-squeezed ... 

   for ilevel = 1:length(plotdat.level)

      plotdat.ges_copy   = guess(:,plotdat.copyindex,  ilevel,:);
      plotdat.anl_copy   = analy(:,plotdat.copyindex,  ilevel,:);

      plotdat.ges_Nposs  = guess(:,plotdat.Npossindex, ilevel,:);
      plotdat.anl_Nposs  = analy(:,plotdat.Npossindex, ilevel,:);
      plotdat.ges_Nused  = guess(:,plotdat.Nusedindex, ilevel,:);
      plotdat.anl_Nused  = guess(:,plotdat.Nusedindex, ilevel,:);
      plotdat.Yrange     = FindRange(plotdat);
      
      % plot by region

      if (plotdat.nregions > 2)
         clf; orient tall
      else 
         clf; orient landscape
      end

      for iregion = 1:plotdat.nregions

         plotdat.region   = iregion;  
         plotdat.myregion = deblank(plotdat.region_names(iregion,:));
         plotdat.title    = sprintf('%s @ %d %s',    ...
                              plotdat.myvarname,     ...
                              plotdat.level(ilevel), ...
                              plotdat.level_units);

         myplot(plotdat);
      end

      % create a postscript file
      print(gcf,'-dpsc','-append',psfname);

   end
end

%----------------------------------------------------------------------
% 'Helper' functions
%----------------------------------------------------------------------

function myplot(plotdat)

   % Interlace the [ges,anl] to make a sawtooth plot.
   % By this point, the middle two dimensions are singletons.
   cg = plotdat.ges_copy(:,:,:,plotdat.region);
   ca = plotdat.anl_copy(:,:,:,plotdat.region);

   g = plotdat.ges_Nposs(:,:,:,plotdat.region);
   a = plotdat.anl_Nposs(:,:,:,plotdat.region);
   nobs_poss = reshape([g a]',2*plotdat.Nbins,1);

   g = plotdat.ges_Nused(:,:,:,plotdat.region);
   a = plotdat.anl_Nused(:,:,:,plotdat.region);
   nobs_used = reshape([g a]',2*plotdat.Nbins,1);

   tg = plotdat.bincenters;
   ta = plotdat.bincenters;
   t = reshape([tg ta]',2*plotdat.Nbins,1);

   % Determine some quantities for the legend
   nobs = sum(nobs_used);
   if ( nobs > 1 )
      mean_prior = mean(cg(isfinite(cg))); 
      mean_post  = mean(ca(isfinite(ca))); 
   else
      mean_prior = NaN;
      mean_post  = NaN;
   end

   string_guess = sprintf('guess:    mean=%.3f', mean_prior);
   string_analy = sprintf('analysis: mean=%.3f', mean_post);

   % Plot the requested quantity on the left axis.
   % The observation count will use the axis on the right.
   % We want to suppress the 'auto' feature of the axis labelling, 
   % so we manually set some values that normally
   % don't need to be set.
   
   ax1 = subplot(plotdat.nregions,1,plotdat.region);

   h1 = plot(ta,ca,'ro-',tg,cg,'k+-','LineWidth',plotdat.linewidth);
   h = legend(string_analy, string_guess);
   legend(h,'boxoff')

   axlims = axis;
   axlims = [axlims(1:2) plotdat.Yrange];
   axis(axlims)

   switch lower(plotdat.copystring)
      case 'bias'
         % plot a zero-bias line
         h4 = line(t,0*t, 'Color','r','Parent',ax1);
         set(h4,'LineWidth',1.5,'LineSTyle',':')
         plotdat.ylabel = sprintf('%s (%s)',plotdat.copystring,plotdat.biasconv);
      otherwise
         plotdat.ylabel = sprintf('%s',plotdat.copystring);
   end
   
   % hokey effort to decide to plot months/days vs. daynum vs.
   ttot = plotdat.bincenters(plotdat.Nbins) - plotdat.bincenters(1) + 1;
   
   if ((plotdat.bincenters(1) > 1000) && (ttot > 5))
      datetick('x',6,'keeplimits','keepticks');
      monstr = datestr(plotdat.bincenters(1),21);
      xlabel(sprintf('month/day - %s start',monstr))
   elseif (plotdat.bincenters(1) > 1000)
      datetick('x',15,'keeplimits','keepticks')
      monstr = datestr(plotdat.bincenters(1),21);
      xlabel(sprintf('%s start',monstr))
   else
      xlabel('days')
   end
   
   % use same X,Y limits for all plots in this region
   nYticks = length(get(ax1,'YTick'));
   ylimits = plotdat.Yrange;
   yinc    = (ylimits(2)-ylimits(1))/(nYticks-1);
   yticks  = ylimits(1):yinc:ylimits(2);
   set(ax1,'YTick',yticks,'Ylim',ylimits)
   
   % create a separate scale for the number of observations
   ax2 = axes('position',get(ax1,'Position'), ...
           'XAxisLocation','top', ...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','b','YColor','b');
   h2 = line(t,nobs_poss,'Color','b','Parent',ax2);
   h3 = line(t,nobs_used,'Color','b','Parent',ax2);
   set(h2,'LineStyle','none','Marker','o');
   set(h3,'LineStyle','none','Marker','+');   

   % use same number of Y ticks and the same X ticks
   
   ylimits = get(ax2,'YLim');
   yinc   = (ylimits(2)-ylimits(1))/(nYticks-1);
   yticks = ylimits(1):yinc:ylimits(2);
   niceyticks = round(10*yticks')/10;
   set(ax2,'XTick',get(ax1,'XTick'),'XTicklabel',get(ax1,'XTicklabel'), ...
           'YTick',          yticks,'YTicklabel',num2str(niceyticks))
       
   set(get(ax2,'Ylabel'),'String','# of obs : o=poss, +=used')
   set(get(ax1,'Ylabel'),'String',plotdat.ylabel)
   set(ax1,'Position',get(ax2,'Position'))
   grid

   if (plotdat.region == 1)
   title({plotdat.title, plotdat.myregion}, ...
         'Interpreter', 'none', 'Fontsize', 12, 'FontWeight', 'bold')
   else
   title(plotdat.myregion, 'Interpreter', 'none', ...
         'Fontsize', 12, 'FontWeight', 'bold')
   end



function BottomAnnotation(main)
% annotates the directory containing the data being plotted
subplot('position',[0.48 0.01 0.04 0.04])
axis off
h = text(0.0,0.5,main);
set(h,'HorizontalAlignment','center', ...
      'VerticalAlignment','middle',...
      'Interpreter','none',...
      'FontSize',8)



function [y,ydims] = FindTemporalVars(x)
% Returns UNIQUE (i.e. base) temporal variable names
if ( ~(isfield(x,'allvarnames') && isfield(x,'allvardims')))
   error('Doh! no ''allvarnames'' and ''allvardims'' components')
end

j = 0;

for i = 1:length(x.allvarnames)
   indx = findstr('time',x.allvardims{i});
   if (indx > 0) 
      j = j + 1;

      basenames{j} = ReturnBase(x.allvarnames{i});
      basedims{j}  = x.allvardims{i};
   end
end

[b,i,j] = unique(basenames);
y     = cell(length(i),1);
ydims = cell(length(i),1);
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

glommed = [y.ges_copy(:);  ...
           y.anl_copy(:) ];
ymin = floor(min(glommed));
ymax =  ceil(max(glommed));
Yrange = [ymin ymax];

if ( isfinite(Yrange) )
   x = [min([Yrange(1) 0.0]) Yrange(2)];
else
   x = [0 1];
end

