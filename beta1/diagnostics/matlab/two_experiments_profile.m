function two_experiments_profile(files, titles, varnames, qtty, prpo)
%
% Each variable gets its own figure.
% Each region gets its own axis.
%
% I'm having some problems preserving the ticks on the left,
% so I used a transparent object. Because of THAT - OpenGL is the
% default rendered, which is normally a pretty low resolution.
% For good reproducibility, manually specify the 'painters' option when
% printing.
%
% files = {'/ptmp/nancy/CSL/Base5/032-061s0_def_reg/obs_diag_output.nc',
%          '/ptmp/thoar/GPS+AIRS/Sep_032-061/obs_diag_output.nc'};
% titles = {'Base5', 'GPS+AIRS'};
% files = {'/fs/image/home/hliu/DART/models/wrf/work/ernesto/ctl/obs_diag_output.nc',
%          '/fs/image/home/hliu/DART/models/wrf/work/ernesto/airsQ/obs_diag_output.nc'};
% titles   = {'Control', 'airsQ'};
% varnames = {'RADIOSONDE_U_WIND_COMPONENT', 'RADIOSONDE_TEMPERATURE'};
% qtty     = 'rmse';     % rmse, spread, totalspread, bias, etc.
% prpo     = 'analysis'; % [analy, analysis, posterior ] == posterior
% prpo     = 'forecast'; % [guess, forecast, prior     ] == prior
%
% two_experiments_profile(files, titles, varnames, qtty, prpo)
% print -dpdf myplot.pdf

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%%--------------------------------------------------------------------
% Decode,Parse,Check the input
%---------------------------------------------------------------------

if (length(files) ~= length(titles))
   error('each file must have a title')
end

NumExp = length(files);

for i = 1:NumExp
   if (exist(files{i},'file') ~= 2)
      error('File %s does not exist',files{i})
   end
end

% set up all the stuff that is common.

commondata = check_compatibility(files, varnames, qtty);
figuredata = setfigure(commondata);

%%--------------------------------------------------------------------
% Set some static data
%---------------------------------------------------------------------

nvars = length(varnames);

for ivar = 1:nvars
   fprintf('Working on %s ...\n',varnames{ivar})
   clf;
   
   for ireg = 1:commondata.nregions
      
      %---------------------------------------------------------------------
      % Getting the data for each experiment
      %---------------------------------------------------------------------
      
      Nlimits = zeros(NumExp,2);  % range of observation count - min, then max
      Dlimits = zeros(NumExp,2);  % range of the data
      Ylimits = zeros(NumExp,2);  % range of the vertical coords
      plotobj = cell(1,NumExp);
      
      for iexp = 1:NumExp
         
         plotobj{iexp} = getvals(files{iexp}, varnames{ivar}, qtty, prpo, ireg);
         plotobj{iexp}.title  = titles{iexp};
         
         Nlimits(iexp,:) = plotobj{iexp}.Nrange;
         Dlimits(iexp,:) = plotobj{iexp}.Drange;
         Ylimits(iexp,:) = plotobj{iexp}.Yrange;
         
      end
      
      %---------------------------------------------------------------------
      % Find nice limits that encompass all experiments
      %---------------------------------------------------------------------
      
      Nrange = [min(Nlimits(:,1)) max(Nlimits(:,2))];
      Drange = [min(Dlimits(:,1)) max(Dlimits(:,2))];
      Yrange = [min(Ylimits(:,1)) max(Ylimits(:,2))];
      
      %---------------------------------------------------------------------
      % Plot all regions - one month to a page
      %---------------------------------------------------------------------
      
      myplot( plotobj, Drange, Yrange, figuredata);
      
   end % of loop around regions
  
   if (commondata.nregions > 3)
     CenterAnnotation(varnames{ivar}) 
   end

   BottomAnnotation(plotobj{1}.timespan,files)
   
   if ( ivar ~= nvars )
      disp('Pausing, hit any key to continue ...')
      pause
   end
   
end  % of loop around variable



%=====================================================================
% End of main function body. Helper functions below.
%=====================================================================



function common = check_compatibility(filenames, varnames, copystring)
%% Trying to prevent the comparison of apples and oranges.
% make sure the diagnostics were generated the same way.

% need to check that the timeframe is the same for all files

% need to check that the region definitions are the same for all files

mystat     = 0;
nexp       = length(filenames);
commondata = cell(1,nexp);
priornames = struct([]);
postenames = struct([]);

for i = 1:length(varnames)
   priornames{i} = sprintf('%s_VPguess',varnames{i});
   postenames{i} = sprintf('%s_VPanaly',varnames{i});
end

for i = 1:nexp
   
   varexist(filenames{i}, {priornames{:}, postenames{:}, 'time', 'time_bounds'})
   
   diminfo = nc_getdiminfo(filenames{i},    'copy'); ncopies   = diminfo.Length;
   diminfo = nc_getdiminfo(filenames{i},'obstypes'); nobstypes = diminfo.Length;
   diminfo = nc_getdiminfo(filenames{i},  'region'); nregions  = diminfo.Length;
   
   commondata{i}.copyindex    = get_copy_index(filenames{i},copystring);
   commondata{i}.ncopies      = ncopies;
   commondata{i}.nobstypes    = nobstypes;
   commondata{i}.nregions     = nregions;
   commondata{i}.times        = nc_varget(filenames{i},'time');
   commondata{i}.time_bnds    = nc_varget(filenames{i},'time_bounds');
   commondata{i}.time_to_skip = nc_attget(filenames{i},nc_global,'time_to_skip');
   commondata{i}.lonlim1      = nc_attget(filenames{i},nc_global,'lonlim1');
   commondata{i}.lonlim2      = nc_attget(filenames{i},nc_global,'lonlim2');
   commondata{i}.latlim1      = nc_attget(filenames{i},nc_global,'latlim1');
   commondata{i}.latlim2      = nc_attget(filenames{i},nc_global,'latlim2');
   
end

% error checking - compare everything to the first experiment
for i = 2:nexp
   
   if (any(commondata{i}.lonlim1 ~= commondata{1}.lonlim1))
      fprintf('The left longitudes of the regions (i.e. lonlim1) are not compatible.\n')
      mystat = 1;
   end
   
   if (any(commondata{i}.lonlim2 ~= commondata{1}.lonlim2))
      fprintf('The right longitudes of the regions (i.e. lonlim2) are not compatible.\n')
      mystat = 1;
   end
   
   if (any(commondata{i}.latlim1 ~= commondata{1}.latlim1))
      fprintf('The bottom latitudes of the regions (i.e. latlim1) are not compatible.\n')
      mystat = 1;
   end
   
   if (any(commondata{i}.latlim2 ~= commondata{1}.latlim2))
      fprintf('The top latitudes of the regions (i.e. latlim2) are not compatible.\n')
      mystat = 1;
   end
   
   if (any(commondata{i}.time_bnds ~= commondata{1}.time_bnds))
      fprintf('The time boundaries of the experiments (i.e. time_bnds) are not compatible.\n')
      mystat = 1;
   end

   if (any(commondata{i}.time_to_skip ~= commondata{1}.time_to_skip))
      fprintf('The time skipped in the experiments (i.e. time_to_skip) are not compatible.\n')
      mystat = 1;
   end
   
end

if mystat > 0
   error('The experiments are not compatible ... stopping.')
end

common = commondata{1};

% Coordinate between time types and dates

timeunits         = nc_attget(filenames{1},'time','units');
calendar          = nc_attget(filenames{1},'time','calendar');
timebase          = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin        = datenum(timebase(1),timebase(2),timebase(3));
timefloats        = zeros(size(commondata{1}.time_to_skip));  % stupid int32 type conversion
timefloats(:)     = commondata{1}.time_to_skip(:);
skip_seconds      = timefloats(4)*3600 + timefloats(5)*60 + timefloats(6);
iskip             = timefloats(3) + skip_seconds/86400;

common.bincenters = commondata{1}.times     + timeorigin;
common.binedges   = commondata{1}.time_bnds + timeorigin;
common.Nbins      = length(common.bincenters);
common.toff       = common.binedges(1) + iskip;

common.timespan   = sprintf('%s through %s', datestr(common.toff), ...
   datestr(max(common.binedges(:))));


%=====================================================================


function plotdat = getvals(fname, varname, copystring, prpo, regionindex )
%% basic function to retrieve plotting data 

if (exist(fname,'file') ~= 2)
   error('%s does not exist',fname)
end

plotdat.fname         = fname;
plotdat.varname       = varname;
plotdat.copystring    = copystring;
plotdat.region        = regionindex;

plotdat.binseparation = nc_attget(fname,nc_global,'bin_separation');
plotdat.binwidth      = nc_attget(fname,nc_global,'bin_width');
time_to_skip          = nc_attget(fname,nc_global,'time_to_skip');
plotdat.lonlim1       = nc_attget(fname,nc_global,'lonlim1');
plotdat.lonlim2       = nc_attget(fname,nc_global,'lonlim2');
plotdat.latlim1       = nc_attget(fname,nc_global,'latlim1');
plotdat.latlim2       = nc_attget(fname,nc_global,'latlim2');
plotdat.biasconv      = nc_attget(fname,nc_global,'bias_convention');

diminfo               = nc_getdiminfo(fname,'region');
plotdat.nregions      = diminfo.Length;
plotdat.region_names  = nc_varget(fname,'region_names');

% Matlab wants character matrices to be Nx1 instead of 1xN.

if (plotdat.nregions == 1 && (size(plotdat.region_names,2) == 1) )
   plotdat.region_names = deblank(plotdat.region_names');
end

% Coordinate between time types and dates

timeunits             = nc_attget(fname,'time','units');
calendar              = nc_attget(fname,'time','calendar');
timebase              = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin            = datenum(timebase(1),timebase(2),timebase(3));
timefloats            = zeros(size(time_to_skip));  % stupid int32 type conversion
timefloats(:)         = time_to_skip(:);
skip_seconds          = timefloats(4)*3600 + timefloats(5)*60 + timefloats(6);
iskip                 = timefloats(3) + skip_seconds/86400;

plotdat.bincenters    = nc_varget(fname,'time');
plotdat.binedges      = nc_varget(fname,'time_bounds');
plotdat.bincenters    = plotdat.bincenters + timeorigin;
plotdat.binedges      = plotdat.binedges   + timeorigin;
plotdat.Nbins         = length(plotdat.bincenters);
plotdat.toff          = plotdat.binedges(1) + iskip;

plotdat.timespan      = sprintf('%s through %s', datestr(plotdat.toff), ...
   datestr(max(plotdat.binedges(:))));

% Get the right indices for the intended variable, regardless of the storage order

plotdat.copyindex     = get_copy_index(fname, copystring);
plotdat.priorvar      = sprintf('%s_VPguess',plotdat.varname);
plotdat.postevar      = sprintf('%s_VPanaly',plotdat.varname);

myinfo.diagn_file     = fname;
myinfo.copyindex      = plotdat.copyindex;
myinfo.regionindex    = plotdat.region;
[start, count]        = GetNCindices(myinfo,'diagn',plotdat.priorvar);
plotdat.prior         = nc_varget(fname, plotdat.priorvar, start, count);
[start, count]        = GetNCindices(myinfo,'diagn',plotdat.postevar);
plotdat.poste         = nc_varget(fname, plotdat.postevar, start, count);

% Now that we know the variable ... get the appropriate vertical information

priordims             = nc_getvarinfo(fname,plotdat.priorvar);
plotdat.levels        = nc_varget(fname,priordims.Dimension{2});
plotdat.level_units   = nc_attget(fname,priordims.Dimension{2},'units');
plotdat.nlevels       = length(plotdat.levels);
plotdat.level_edges   = nc_varget(fname,sprintf('%s_edges',priordims.Dimension{2}));

if (plotdat.levels(1) > plotdat.levels(plotdat.nlevels))
   plotdat.YDir = 'reverse';
else
   plotdat.YDir = 'normal';
end

%% Determine data limits - Do we use prior and/or posterior
%  always make sure we have a zero bias line ...

plotdat.useposterior = 0;
plotdat.useprior     = 0;
switch lower(prpo)
   case {'analy','analysis','posterior'}
      plotdat.useposterior = 1;
      bob = plotdat.poste(:);
   case {'guess','forecast','prior'}
      plotdat.useprior = 1;
      bob = plotdat.prior(:);
   otherwise
      plotdat.useposterior = 1;
      plotdat.useprior = 1;
      bob = [plotdat.prior(:) ; plotdat.poste(:)];   % one long array
end

switch copystring
   case {'bias'}
      dmin = min( [ min(bob) 0.0 ] );
      dmax = max( [ max(bob) 0.0 ] );
      plotdat.Drange = [ dmin dmax ];
      plotdat.xlabel = sprintf('%s %s',copystring, plotdat.biasconv);
   otherwise
      plotdat.Drange = [min(bob) max(bob)];
      plotdat.xlabel = copystring;
end

%% Get the right indices for the number of observations possible
%  Get the right indices for the number of observations used
%  FIXME - should the number rejected because of incoming QC be disqualified

myinfo.diagn_file = fname;
myinfo.copyindex  = get_copy_index(fname, 'Nposs');
[start, count]    = GetNCindices(myinfo,'diagn',plotdat.priorvar);
plotdat.nposs     = nc_varget(fname, plotdat.priorvar, start, count);

if ( plotdat.useprior )
   myinfo.copyindex = get_copy_index(fname, 'Nused');
   [start, count]   = GetNCindices(myinfo,'diagn',plotdat.priorvar);
   plotdat.nused    = nc_varget(fname, plotdat.priorvar, start, count);
else
   myinfo.copyindex = get_copy_index(fname, 'Nused');
   [start, count]   = GetNCindices(myinfo,'diagn',plotdat.postevar);
   plotdat.nused    = nc_varget(fname, plotdat.postevar, start, count);
end

%% Set the last of the ranges

plotdat.Yrange = [min(plotdat.level_edges) max(plotdat.level_edges)];
plotdat.Nrange = [min(plotdat.nused(:))    max(plotdat.nposs(:))];


%=====================================================================


function myplot( plotdat, Drange, Yrange, figdata)
%% Create graphic for one region - for all experiments.

Nexp    = length(plotdat);
iregion = plotdat{1}.region;
ax1     = subplot('position',figdata.plotlims(iregion,:));

Stripes(Drange, plotdat{1}.level_edges, plotdat{1}.level_units, plotdat{1}.nregions);
set(ax1,'YDir',plotdat{1}.YDir,'YTick',sort(plotdat{1}.levels),'Layer','top')
set(ax1,'YAxisLocation','left')

% draw the results of the experiments, priors and posteriors -
% each with their own line type.
iexp   = 0;
hd     = [];   % handle to an unknown number of data lines
legstr = {[]}; % strings for the legend

hold on
for i = 1:Nexp
   
   if ( plotdat{i}.useprior )
      iexp         = iexp + 1;
      lty          = sprintf('%s%s%s',figdata.expcolors{i},figdata.prpolines{1}, ...
         figdata.expsymbols{i});
      hd(iexp)     = plot(plotdat{i}.prior, plotdat{i}.levels, lty,'LineWidth', 2.0);
      legstr{iexp} = sprintf('%s forecast',plotdat{i}.title);
   end
   
   if ( plotdat{i}.useposterior )
      iexp         = iexp + 1;
      lty          = sprintf('%s%s%s',figdata.expcolors{i},figdata.prpolines{2}, ...
         figdata.expsymbols{i});
      hd(iexp)     = plot(plotdat{i}.poste, plotdat{i}.levels, lty,'LineWidth', 2.0);
      legstr{iexp} = sprintf('%s analysis',plotdat{i}.title);
   end
end
hold off;

switch plotdat{1}.copystring
   case {'bias'}
      zeroline = line([0 0],Yrange,'Color','k','Parent',ax1);
      set(zeroline,'LineWidth',1.0,'LineStyle','-.')
   otherwise
end

% Create another axes to use for plotting the observation counts

ax2 = axes('position',get(ax1,'Position'), ...
   'XAxisLocation','top', ...
   'YAxisLocation','right',...
   'Color','none',...
   'XColor','b','YColor','b',...
   'YLim',get(ax1,'YLim'), ...
   'YDir',get(ax1,'YDir'));

% Plot the data, which sets the range of the axis
for i = 1:Nexp
   h2 = line(plotdat{i}.nposs, plotdat{i}.levels,'Color',figdata.expcolors{i},'Parent',ax2);
   h3 = line(plotdat{i}.nused, plotdat{i}.levels,'Color',figdata.expcolors{i},'Parent',ax2);
   set(h2,'LineStyle','none','Marker',figdata.expsymbols{i});
   set(h3,'LineStyle','none','Marker','+');
end

% use same Y ticks
set(ax2,'YTick',get(ax1,'YTick'), ...
   'YTicklabel',get(ax1,'YTicklabel'));

% use the same X ticks, but find the right label values
xscale = matchingXticks(ax1,ax2);

% Annotate the whole thing - gets pretty complicated for multiple
% regions on one page. Trying to maximize content, minimize clutter.
% Any plot object will do for annotating region,levels,etc

annotate( ax1, ax2, plotdat{1}, figdata, xscale)

lh = legend(hd,legstr,'Location','NorthWest');
legend(lh,'boxoff','Interpreter','none');

% The legend linesizes should match - 2 is hardwired - suprises me.

set(lh,'FontSize',figdata.fontsize);
kids = get(lh,'Children');
set(kids,'LineWidth',2.0);


%=====================================================================


function h = Stripes(x, edges, units, nregions)
%% plot the subtle background stripes that indicate the vertical
%  extent of what was averaged.
%
%  axlims(3) should be conditional on the observation vertical coordinate:
%  values for pressure coordinates are inappropriate for height coord.
%  It also assumes 4 plots/page, but 2 works better for plotting all levels of CAM5.
%  That requires a smaller % of vertical range for the legend.

% plot two little dots at the corners to make Matlab
% choose some plot limits. Given those nice limits and
% tick labels ... KEEP THEM. Later, make the dots invisible.

h = plot([min(x) max(x)],[min(edges) max(edges)]);
axlims          = axis;
region_factor   = floor((nregions+1)/2.0);
legend_fraction = 0.1375;

% partial fix to legend space; add in option for vert coord = height.

switch lower(units)
   case 'hpa'
      axlims(4) = max(edges);
      axlims(3) = min(edges) - region_factor*legend_fraction*(axlims(4)-min(edges));
   case 'm'
      axlims(3) = min(edges) ;
      axlims(4) = max(edges) + region_factor*legend_fraction*(max(edges)-axlims(3));
   otherwise
end
axis(axlims)

% set up list of x,y values defining corner of every other stripe

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


function annotate(ax1, ax2, plotobj, figdata, xscale)
%% Each configuration of subplots is exploited.

if ( plotobj.nregions == 1 )
   %% One figure ... everything gets annotated.
   set(get(ax1,'Ylabel'),'String',plotobj.level_units,'Interpreter','none')
   set(get(ax2,'Ylabel'),'String',plotobj.level_units,'Interpreter','none')
   set(get(ax2,'Xlabel'),'String',['# of obs (o=poss, +=used) x' int2str(uint32(xscale))])
   set(get(ax1,'Xlabel'),'String',plotobj.xlabel,'Interpreter','none')
   
   th = title({deblank(plotobj.region_names(plotobj.region,:)), plotobj.varname});
   set(th,'Interpreter','none','FontSize',figdata.fontsize);
   
elseif ( plotobj.nregions == 2 )
   %% Two figures ... side-by-side
   if (plotobj.region == 1) % figure on left ... no anno on the right
      set(ax2,'YTickLabel',[])  
      set(get(ax1,'Ylabel'),'String',plotobj.level_units,'Interpreter','none')
   else                     % figure on right .. no anno on the left
      set(ax1,'YTickLabel',[])
      set(get(ax2,'Ylabel'),'String',plotobj.level_units,'Interpreter','none')
   end
   set(get(ax1,'Xlabel'),'String',plotobj.xlabel,'Interpreter','none')
   set(get(ax2,'Xlabel'),'String',['# of obs (o=poss, +=used) x' int2str(uint32(xscale))])
   th = title({deblank(plotobj.region_names(plotobj.region,:)), plotobj.varname});
   set(th,'Interpreter','none','FontSize',figdata.fontsize);
   
elseif ( plotobj.nregions == 3 )
   
   if (plotobj.region == 1)
      set(ax2,'YTickLabel',[])
      set(get(ax1,'Ylabel'),'String',plotobj.level_units,'Interpreter','none')
      set(get(ax1,'Xlabel'),'String',plotobj.copystring,'Interpreter','none')
      titlestring = {deblank(plotobj.region_names(plotobj.region,:))};
   elseif (plotobj.region == 2)
      set(ax1,'YTickLabel',[])
      set(ax2,'YTickLabel',[])
      set(get(ax1,'Xlabel'),'String',plotobj.xlabel,'Interpreter','none')
      titlestring = {deblank(plotobj.region_names(plotobj.region,:)), plotobj.varname};
   else
      set(ax1,'YTickLabel',[])
      set(get(ax2,'Ylabel'),'String',plotobj.level_units,'Interpreter','none')
      set(get(ax1,'Xlabel'),'String',plotobj.copystring,'Interpreter','none')
      titlestring = {deblank(plotobj.region_names(plotobj.region,:))};
   end
   set(get(ax2,'Xlabel'),'String',['# of obs (o=poss, +=used) x' int2str(uint32(xscale))])
   th = title(titlestring);
   set(th,'Interpreter','none','FontSize',figdata.fontsize);
   
elseif ( plotobj.nregions == 4 )

   titlestring = deblank(plotobj.region_names(plotobj.region,:));
   th = title(titlestring);
   set(th,'Interpreter','none','FontSize',figdata.fontsize);

   set(get(ax1,'Xlabel'),'String',plotobj.copystring,'Interpreter','none')

   xh = get(ax2,'Xlabel');  % handle to number of obs possible string
   set(xh,'String',['# of obs (o=poss, +=used) x' int2str(uint32(xscale))])
   set(xh,'FontSize',8)

   % tile them 2-by-2  - adjust Y label 

   if (    plotobj.region == 1) % figure on left ... no anno on the right
      set(ax2,'YTickLabel',[])  
      set(get(ax1,'Ylabel'),'String',plotobj.level_units,'Interpreter','none')

   elseif (plotobj.region == 3) % figure on left ... no anno on the right
      set(ax2,'YTickLabel',[])  
      set(get(ax1,'Ylabel'),'String',plotobj.level_units,'Interpreter','none')

   elseif (plotobj.region == 2) % figure on  right ... no anno on the left
      set(ax1,'YTickLabel',[])  
      set(get(ax2,'Ylabel'),'String',plotobj.level_units,'Interpreter','none')

   else                         % figure on right .. no anno on the left
      set(ax1,'YTickLabel',[])
      set(get(ax2,'Ylabel'),'String',plotobj.level_units,'Interpreter','none')

   end

end


%=====================================================================


function CenterAnnotation(main)
%% annotates the variable being plotted
%  the bottom of the top row of figures is at 0.580
%  the top of the bottom row of figures is at 0.480 (0.125+0.355)
%  so the center is 0.53:

subplot('position',[0.48 0.51 0.04 0.04])
axis off
h = text(0.5,0.5,main);
set(h,'HorizontalAlignment','center', ...
   'VerticalAlignment','bottom', ...
   'Interpreter','none', ...
   'FontSize',12, ...
   'FontWeight','bold')


%=====================================================================


function BottomAnnotation(timespan,filenames)
%% annotates the timespan of the data being plotted as well as
%  the filenames containing the data being plotted

subplot('position',[0.10 0.01 0.8 0.05])
axis off

nfiles = length(filenames);

% vertical spacing for text

if isempty(timespan)
   dy = 1.0/(nfiles+1);
else
   dy = 1.0/(nfiles+2);
end

% The first thing is the timeframe, then all the files.

ifile = 1;
ty = 1.0 - (ifile-1)*dy;
h  = text(0.5, ty, timespan);
set(h, 'Interpreter', 'none', 'FontSize', 8);
set(h, 'HorizontalAlignment','center');

for ifile = 1:nfiles
   main = filenames{ifile};

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
   
   ty = 1.0 - ifile*dy;
   h = text(0.5, ty, string1);
   set(h, 'Interpreter', 'none', 'FontSize', 8);
   set(h, 'HorizontalAlignment','center');
end


%=====================================================================


function varexist(filename, varnames)
%% We already know the file exists by this point.
% Lets check to make sure that file contains all needed variables.

nvars  = length(varnames);
gotone = ones(1,nvars);

for i = 1:nvars
   gotone(i) = nc_isvar(filename,varnames{i});
   if ( ~ gotone(i) )
      fprintf('\n%s is not a variable in %s\n',varnames{i},filename)
   end
end

if ~ all(gotone)
   error('missing required variable ... exiting')
end


%=====================================================================


function figdata = setfigure(commondata)
%%
%  the hardest part is figuring out the gaps, etc.
%  gapwidth = (1 - (nregions*axiswidth + 0.14))/(nregions - 1)
%  0.14 is actually left margin + right margin ... 0.07 + 0.07

plotlims = zeros(commondata.nregions,4);

if (commondata.nregions > 4)
   error('Cannot handle %d regions - 4 is the max.',commondata.nregions)
   
elseif (commondata.nregions == 4)
   orientation   = 'tall';
   fontsize      = 10;
   plotlims(1,:) = [0.100 0.580 0.375 0.355];
   plotlims(2,:) = [0.525 0.580 0.375 0.355];
   plotlims(3,:) = [0.100 0.125 0.375 0.355];
   plotlims(4,:) = [0.525 0.125 0.375 0.355];
   
elseif (commondata.nregions == 3)
   orientation   = 'landscape';
   fontsize      = 12;
   plotlims(1,:) = [0.0700 0.14 0.265 0.725];
   plotlims(2,:) = [0.3675 0.14 0.265 0.725];
   plotlims(3,:) = [0.6650 0.14 0.265 0.725];
   
elseif (commondata.nregions == 2)
   orientation   = 'landscape';
   fontsize      = 12;
   plotlims(1,:) = [0.07 0.14 0.41 0.725];
   plotlims(2,:) = [0.52 0.14 0.41 0.725];
   
elseif (commondata.nregions == 1)
   orientation   = 'tall';
   fontsize      = 16;
   plotlims(1,:) = [0.10 0.12 0.8 0.75];
end

figdata = struct('expcolors',  {{'k','r','g','m','b','c','y'}}, ...
   'expsymbols', {{'o','s','d','p','h','s','*'}}, ...
   'prpolines',  {{'-',':'}}, 'plotlims', plotlims, ...
   'fontsize',fontsize, 'orientation',orientation);

clf; orient(gcf, orientation); wysiwyg


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

