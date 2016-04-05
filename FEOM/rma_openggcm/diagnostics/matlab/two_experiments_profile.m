function two_experiments_profile(files, titles, obsnames, copy, prpo, varargin)
% Plot two or more experiments on the same axis.
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
% USAGE: two_experiments_profile(files, titles, obsnames, copy, prpo, 'level', 1)
%
% files    : Cell array containing the locations of the obs_diag_output.nc 
%            files to compare. Each file is presumed to be the results from
%            a single experiment.
%
% titles   : Cell array containing the titles used to annotate each of the experiments. 
%
% obsnames : Cell array containing the strings of each observation type to plot.
%            Each observation type will be plotted in a separate graphic.
%
% copy     : string defining the metric of interest. 'rmse', 'spread', etc.
%            Possible values are available in the netcdf 'CopyMetaData' variable.
%            (ncdump -v CopyMetaData obs_diag_output.nc)
%
% prpo     : string defining whether to plot the prior or posterior metrics.
%            Due to the amount of information already plotted, we made a
%            conscious decision not to support plotting both prior and posterior
%            on the same plot.
%
% level    : The index of the level to plot. Defaults to level 1.
%
% 
% OUTPUT: A .pdf of each graphic is created. Each .pdf has a name that 
%         reflects the variable, quantity, and region being plotted.
%
% EXAMPLE
%
% files = {'/ptmp/nancy/CSL/Base5/032-061s0_def_reg/obs_diag_output.nc',
%          '/ptmp/thoar/GPS+AIRS/Sep_032-061/obs_diag_output.nc'};
% titles = {'Base5', 'GPS+AIRS'};
% files = {'/glade/scratch/nancy/fvdiags/obs_diag_2005_08.nc', ...
%  '/glade/scratch/raeder/SE_NCEP_assim1/Diag_hemi_poles_2005.8.1-30/obs_diag_output.nc'};
% titles   = {'FV', 'SE'};
% obsnames = {'RADIOSONDE_TEMPERATURE', ...
%             'RADIOSONDE_U_WIND_COMPONENT', 'RADIOSONDE_V_WIND_COMPONENT'};
% copy     = 'rmse';     % rmse, spread, totalspread, bias, etc.
% prpo     = 'analysis'; % [analy, analysis, posterior ] == posterior
% prpo     = 'forecast'; % [guess, forecast, prior     ] == prior
% prpo     = 'both';
%
% two_experiments_profile(files, titles, obsnames, copy, prpo)
%
% Example 2: restrict the data limits to the data in a certain vertical area. 
%            In this case, the observations using a pressure vertical coordinate
%            between +Inf (the surface) and 100hPa (inclusive) are used to determine
%            the scale. All values will be plotted, but the highest levels may be
%            clipped. The optional argument pairs at the end consist of a string
%            and a length 2 array specifying the [bottom top] levels to consider.
%
% two_experiments_profile(files, titles, obsnames, copy, prpo,'plevel',[Inf 100])
%
% two_experiments_profile(files, titles, obsnames, copy, prpo, ...
%            'plevel',[Inf 100],'mlevel',[1 10],'hlevel',[0 20000])

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%%--------------------------------------------------------------------
% Decode,Parse,Check the input
%---------------------------------------------------------------------
defaultPlevels = [ Inf  0 ];
defaultHlevels = [-Inf Inf];
defaultMlevels = [  1  Inf];
p = inputParser;
addRequired(p,'files',@iscell);
addRequired(p,'titles',@iscell);
addRequired(p,'obsnames',@iscell);
addRequired(p,'copy',@ischar);
addRequired(p,'prpo',@ischar);
addParamValue(p,'plevel',defaultPlevels,@isnumeric);
addParamValue(p,'hlevel',defaultHlevels,@isnumeric);
addParamValue(p,'mlevel',defaultMlevels,@isnumeric);
parse(p, files, titles, obsnames, copy, prpo, varargin{:});

% if you want to echo the input
% disp(['files   : ', p.Results.files])
% disp(['titles  : ', p.Results.titles])
% disp(['obsnames: ', p.Results.obsnames])
% disp(['copy    : ', p.Results.copy])
% disp(['prpo    : ', p.Results.prpo])
% fprintf( 'plevel : %f %f \n', p.Results.plevel)
% fprintf( 'hlevel : %f %f \n', p.Results.hlevel)
% fprintf( 'mlevel : %f %f \n', p.Results.mlevel)

if ~isempty(fieldnames(p.Unmatched))
   disp('Extra inputs:')
   disp(p.Unmatched)
end

% if ~isempty(p.UsingDefaults)
%    disp('Using defaults: ')
%    disp(p.UsingDefaults)
% end

if (numel(p.Results.plevel) ~= 2)
   error('plevel must be an array of length two ... [bottom top]')
end

if (numel(p.Results.hlevel) ~= 2)
   error('hlevel must be an array of length two ... [bottom top]')
end

if (numel(p.Results.mlevel) ~= 2)
   error('mlevel must be an array of length two ... [bottom top]')
end

% Now that the input passes sanity checks ...

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

commondata = check_compatibility(files, obsnames, copy);
figuredata = setfigure(NumExp);

%%--------------------------------------------------------------------
% Set some static data
%---------------------------------------------------------------------

nvars = length(obsnames);

for ivar = 1:nvars
   fprintf('Working on %s ...\n',obsnames{ivar})

   for iregion = 1:commondata.nregions

      %---------------------------------------------------------------------
      % Getting the data for each experiment
      %---------------------------------------------------------------------

      Nlimits = zeros(NumExp,2);  % range of observation count - min, then max
      Dlimits = zeros(NumExp,2);  % range of the data
      Ylimits = zeros(NumExp,2);  % range of the vertical coords
      plotobj = cell(1,NumExp);

      for iexp = 1:NumExp

         plotobj{iexp} = getvals(files{iexp}, obsnames{ivar}, copy, prpo, iregion, p);
         plotobj{iexp}.title  = titles{iexp};

         Nlimits(iexp,:) = plotobj{iexp}.Nrange;
         Dlimits(iexp,:) = plotobj{iexp}.Drange;
         Ylimits(iexp,:) = plotobj{iexp}.Yrange;

      end

      %---------------------------------------------------------------------
      % Find nice limits that encompass all experiments
      % Note that Dlimits has been constructed by ignoring the top levels.
      %---------------------------------------------------------------------

      Nrange = [min(Nlimits(:,1)) max(Nlimits(:,2))];
      Drange = [min(Dlimits(:,1)) max(Dlimits(:,2))];
      Yrange = [min(Ylimits(:,1)) max(Ylimits(:,2))];
      span = abs(Drange(2) - Drange(1))* 0.05;
      Drange(1) = Drange(1) - span;
      Drange(2) = Drange(2) + span;

      %---------------------------------------------------------------------
      % Plot all regions - one region to a page
      %---------------------------------------------------------------------

      myplot(plotobj, Drange, Yrange, figuredata);

      BottomAnnotation(files)

      psfname = sprintf('%s_%s_region%d_profile_%dexp', ...
                obsnames{ivar}, plotobj{1}.copystring, iregion, NumExp);
      print(iregion,'-dpdf',psfname)

   end % of loop around regions

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


function plotdat = getvals(fname, varname, copystring, prpo, regionindex, opt )
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
plotdat.Npossindex    = get_copy_index(fname, 'Nposs');
plotdat.Nusedindex    = get_copy_index(fname, 'Nused');
plotdat.NQC4index     = get_copy_index(fname, 'N_DARTqc_4');
plotdat.NQC5index     = get_copy_index(fname, 'N_DARTqc_5');
plotdat.NQC6index     = get_copy_index(fname, 'N_DARTqc_6');
plotdat.NQC7index     = get_copy_index(fname, 'N_DARTqc_7');

plotdat.priorvar      = sprintf('%s_VPguess',plotdat.varname);
plotdat.postevar      = sprintf('%s_VPanaly',plotdat.varname);

myinfo.diagn_file     = fname;
myinfo.copyindex      = plotdat.copyindex;
myinfo.regionindex    = plotdat.region;
[start, count]        = GetNCindices(myinfo,'diagn',plotdat.priorvar);
plotdat.prior         = nc_varget(fname, plotdat.priorvar, start, count);
[start, count]        = GetNCindices(myinfo,'diagn',plotdat.postevar);
plotdat.poste         = nc_varget(fname, plotdat.postevar, start, count);
plotdat.trusted       = nc_read_att(fname, plotdat.priorvar, 'TRUSTED');
if (isempty(plotdat.trusted)), plotdat.trusted = 'NO'; end

% Now that we know the variable ... get the appropriate vertical information

priordims             = nc_getvarinfo(fname,plotdat.priorvar);
plotdat.levels        = nc_varget(fname,priordims.Dimension{2});
plotdat.level_units   = nc_attget(fname,priordims.Dimension{2},'units');
plotdat.nlevels       = length(plotdat.levels);
plotdat.level_edges   = nc_varget(fname,sprintf('%s_edges',priordims.Dimension{2}));

plotdat.YDir = 'normal';
inds = 1:plotdat.nlevels;

% find the levels of interest for setting the data limits
switch lower(priordims.Dimension{2})
    case {'plevel'}
        plotdat.YDir = 'reverse';
        inds = find((plotdat.levels <= opt.Results.plevel(1)) & ...
                    (plotdat.levels >= opt.Results.plevel(2)));
    case {'hlevel'}
        inds = find((plotdat.levels >= opt.Results.hlevel(1)) & ...
                    (plotdat.levels <= opt.Results.hlevel(2)));
    case {'mlevel'}
        inds = find((plotdat.levels >= opt.Results.mlevel(1)) & ...
                    (plotdat.levels <= opt.Results.mlevel(2)));
    otherwise
end

%% Determine data limits - Do we use prior and/or posterior
%  always make sure we have a zero bias line ...

plotdat.useposterior = 0;
plotdat.useprior     = 0;

switch lower(prpo)
   case {'analy','analysis','posterior'}
      plotdat.useposterior = 1;
      plotdat.prpo = 'analysis';
      bob = plotdat.poste(inds);
   case {'guess','forecast','prior'}
      plotdat.useprior = 1;
      plotdat.prpo = 'forecast';
      bob = plotdat.prior(inds);
   otherwise
      plotdat.useposterior = 1;
      plotdat.useprior = 1;
      plotdat.prpo = 'forecast and analysis';
      bob = [plotdat.prior(inds) ; plotdat.poste(inds)];   % one long array
end

switch copystring
   case {'bias'}
      dmin = min( [ min(bob) 0.0 ] );
      dmax = max( [ max(bob) 0.0 ] );
      plotdat.Drange = [ dmin dmax ];
      plotdat.xlabel = sprintf('%s (%s)',copystring, plotdat.biasconv);
   case {'rmse'}
      plotdat.Drange = [0.0 max(bob)];
      plotdat.xlabel = copystring;
   otherwise
      plotdat.Drange = [min(bob) max(bob)];
      plotdat.xlabel = copystring;
end

%% Get the indices for the number of observations possible
%  Get the indices for the number of observations used
%  The number of obs possible is affected by namelist selection of
%  which observations to assimilate, and what incoming QC is 'good'.

myinfo.diagn_file = fname;
myinfo.copyindex  = plotdat.Npossindex;
[start, count]    = GetNCindices(myinfo,'diagn',plotdat.priorvar);
plotdat.nposs     = nc_varget(fname, plotdat.priorvar, start, count);

myinfo.copyindex  = plotdat.NQC5index;
[start, count]    = GetNCindices(myinfo,'diagn',plotdat.priorvar);
plotdat.Nqc5      = nc_varget(fname, plotdat.priorvar, start, count);
plotdat.nposs     = plotdat.nposs - plotdat.Nqc5;

myinfo.copyindex  = plotdat.NQC6index;
[start, count]    = GetNCindices(myinfo,'diagn',plotdat.priorvar);
plotdat.Nqc6      = nc_varget(fname, plotdat.priorvar, start, count);
plotdat.nposs     = plotdat.nposs - plotdat.Nqc6;

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

figure(iregion);
clf(iregion); orient(figdata.orientation); wysiwyg
ax1 = subplot('position',figdata.position);

Stripes(Drange, plotdat{1}.level_edges, plotdat{1}.level_units, Nexp);
set(ax1,'YDir',plotdat{1}.YDir,'YTick',sort(plotdat{1}.levels),'Layer','top')
set(ax1,'YAxisLocation','left','FontSize',figdata.fontsize)

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
      hd(iexp)     = plot(plotdat{i}.prior, plotdat{i}.levels, lty,'LineWidth', ...
                          figdata.linewidth);
      legstr{iexp} = sprintf('%s Prior',plotdat{i}.title);
   end

   if ( plotdat{i}.useposterior )
      iexp         = iexp + 1;
      lty          = sprintf('%s%s%s',figdata.expcolors{i},figdata.prpolines{2}, ...
                             figdata.expsymbols{i});
      hd(iexp)     = plot(plotdat{i}.poste, plotdat{i}.levels, lty,'LineWidth', ...
                          figdata.linewidth);
      legstr{iexp} = sprintf('%s Posterior',plotdat{i}.title);
   end
end
hold off;

switch plotdat{1}.copystring
   case {'bias','rmse'}
      zeroline = line([0 0],get(ax1,'YLim'),'Color',[0 100 0]/255,'Parent',ax1);
      set(zeroline,'LineWidth',2.5,'LineStyle','-')
   otherwise
end

set(ax1,'XLim',Drange)

axlims = axis;
dx = (axlims(2) - axlims(1))/20;
if strcmpi('normal',get(ax1,'YDir'))
   ty = axlims(3);
else
   ty = axlims(4);
end

for i = 1:Nexp
   % If the observation is trusted, reference that somehow
   switch lower(plotdat{i}.trusted)
      case 'true'
         tx = axlims(2) + (i*dx);
         h = text(tx,ty,sprintf('TRUSTED OBSERVATION in %s',plotdat{i}.title));
         set(h, 'FontSize', 20, 'Rotation', 90, ...
            'VerticalAlignment', 'middle', 'Interpreter', 'none')
      otherwise
   end
end

% Create another axes to use for plotting the observation counts

ax2 = axes('position',get(ax1,'Position'), ...
   'XAxisLocation','top', ...
   'YAxisLocation','right',...
   'Color','none', ...
   'XColor','b', ...
   'YColor',get(ax1,'YColor'), ...
   'YLim',get(ax1,'YLim'), ...
   'YDir',get(ax1,'YDir'), ...
   'FontSize',get(ax1,'FontSize'));

% Plot the data, which sets the range of the axis
for i = 1:Nexp
   h2 = line(plotdat{i}.nposs, plotdat{i}.levels,'Color',figdata.expcolors{i},'Parent',ax2);
   h3 = line(plotdat{i}.nused, plotdat{i}.levels,'Color',figdata.expcolors{i},'Parent',ax2);
   set(h2,'LineStyle','none','Marker','o','MarkerSize',10);
   set(h3,'LineStyle','none','Marker','*','MarkerSize',10);
end

% use same Y ticks but no labels
set(ax2,'YTick',get(ax1,'YTick'), 'YTicklabel',[]);

% use the same X ticks, but find the right label values
xscale = matchingXticks(ax1,ax2);

% Annotate the whole thing - gets pretty complicated for multiple
% regions on one page. Trying to maximize content, minimize clutter.
% Any plot object will do for annotating region,levels,etc

annotate( ax1, ax2, plotdat{1}, figdata, xscale)

lh = legend(hd,legstr,'Location','Best');
set(lh,'Interpreter','none','Box','off');

% The legend linesizes should match - 2 is hardwired - suprises me.

set(lh,'FontSize',figdata.fontsize);
kids = get(lh,'Children');
set(kids,'LineWidth',figdata.linewidth);


%=====================================================================


function h = Stripes(x, edges, units, nexp)
%% plot the subtle background stripes that indicate the vertical
%  extent of what was averaged.
%
%  axlims(3) should be conditional on the observation vertical coordinate:
%  values for pressure coordinates are inappropriate for height coord.

% plot two little dots at the corners to make Matlab
% choose some plot limits. Given those nice limits and
% tick labels ... KEEP THEM. Later, make the dots invisible.

h = plot([min(x) max(x)],[min(edges) max(edges)]);
axlims          = axis;
legend_fraction = 0.05 * nexp + 0.005;

% partial fix to legend space; add in option for vert coord = height.

switch lower(units)
   case 'hpa'
      axlims(4) = max(edges);
      axlims(3) = min(edges) - legend_fraction*(axlims(4)-min(edges));
   case 'm'
      axlims(3) = min(edges) ;
      axlims(4) = max(edges) + legend_fraction*(max(edges)-axlims(3));
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

%% One figure ... everything gets annotated.
set(get(ax1,'Ylabel'),'String',plotobj.level_units, ...
                      'Interpreter','none','FontSize',figdata.fontsize)
set(get(ax1,'Xlabel'),'String',{plotobj.xlabel,plotobj.timespan}, ...
                      'Interpreter','none','FontSize',figdata.fontsize)
set(get(ax2,'Xlabel'),'String', ...
   ['# of obs (o=possible, \ast=assimilated) x' int2str(uint32(xscale))],'FontSize',figdata.fontsize)

th = title({deblank(plotobj.region_names(plotobj.region,:)), plotobj.varname});
set(th,'Interpreter','none','FontSize',figdata.fontsize,'FontWeight','bold');


%=====================================================================


function BottomAnnotation(filenames)
%% annotates the filenames containing the data being plotted

nfiles  = length(filenames);
dy      = 1.0/(nfiles+3);
yheight = 0.0225*(nfiles+1);

subplot('position',[0.10 0.01 0.8 yheight])
axis off

% list all the files.

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

   ty = 1.0 - (ifile+2)*dy;
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


function figdata = setfigure(nexp)
%%
%  figure out a page layout
%  extra space at the bottom for the date/file annotation
%  extra space at the top because the titles have multiple lines

ybot = 0.06 + nexp*0.035;  % room for dates/files
ytop = 0.125;              % room for title (always 2 lines)
dy   = 1.0 - ytop - ybot;
orientation = 'tall';
fontsize    = 16;
position    = [0.15 ybot 0.7 dy];
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

