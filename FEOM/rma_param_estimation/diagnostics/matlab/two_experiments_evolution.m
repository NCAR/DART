function two_experiments_evolution(files, titles, obsnames, copy, prpo, varargin)
%% two_experiments_evolution  compares multiple netcdf files created by obs_diag.
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
% USAGE: two_experiments_evolution(files, titles, obsnames, copy, prpo, 'level', 1)
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
% OUTPUT: A .pdf of each graphic is created. Each .pdf has a name that 
%         reflects the variable, quantity, and region being plotted.
%
% EXAMPLE
%
% files = {'/ptmp/nancy/CSL/Base5/032-061s0_def_reg/obs_diag_output.nc',
%          '/ptmp/thoar/GPS+AIRS/Sep_032-061/obs_diag_output.nc'};
% titles = {'Base5', 'GPS+AIRS'};
% obsnames = {'RADIOSONDE_U_WIND_COMPONENT', 'RADIOSONDE_TEMPERATURE'};
%
% copy  = 'spread';   % rmse, spread, totalspread, bias, etc.
% prpo  = 'analysis'; % [analy, analysis, posterior ] == posterior
% prpo  = 'forecast'; % [guess, forecast, prior     ] == prior
%
% two_experiments_evolution(files, titles, obsnames, copy, prpo, 'level', 1)

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%%--------------------------------------------------------------------
% Decode,Parse,Check the input
%---------------------------------------------------------------------

default_level = 1;
p = inputParser;

addRequired(p,'files',@iscell);
addRequired(p,'titles',@iscell);
addRequired(p,'obsnames',@iscell);
addRequired(p,'copy',@ischar);
addRequired(p,'prpo',@ischar);
addParamValue(p,'level',default_level,@isnumeric);
parse(p, files, titles, obsnames, copy, prpo, varargin{:});

% if you want to echo the input
% disp(['files   : ', p.Results.files])
% disp(['titles  : ', p.Results.titles])
% disp(['obsnames: ', p.Results.obsnames])
% fprintf('level : %d \n', p.Results.level)

if ~isempty(fieldnames(p.Unmatched))
   disp('Extra inputs:')
   disp(p.Unmatched)
end

NumExp = length(files);

for i = 1:NumExp
   if (exist(files{i},'file') ~= 2)
      error('File %s does not exist',files{i})
   end
end

if (NumExp ~= length(titles))
   error('each file must have an experiment title')
end

%% set up all the stuff that is common.

commondata = check_compatibility(files, obsnames, copy);
figuredata = setfigure(NumExp);

%%--------------------------------------------------------------------
% Set some static data
%---------------------------------------------------------------------

nvars = length(obsnames);

for ivar = 1:nvars
   fprintf('Working on %s ...\n',obsnames{ivar})

   %------------------------------------------------------------------------
   % Plot each region in a separate figure window.
   %------------------------------------------------------------------------

   for iregion = 1:commondata.nregions

      figure(iregion);
      clf(iregion);
      orient(figuredata.orientation);
      wysiwyg;

      %---------------------------------------------------------------------
      % 1) Get the data for each experiment
      % 2) plot the data
      % 3) annotate
      %---------------------------------------------------------------------

      for iexp = 1:NumExp

         plotobj{iexp} = getvals(files{iexp}, obsnames{ivar}, copy, prpo, iregion, p.Results.level);
         plotobj{iexp}.title        = titles{iexp};
         plotobj{iexp}.nregions     = commondata.nregions;
         plotobj{iexp}.region_names = commondata.region_names;

      end

      myplot(plotobj, figuredata);

      BottomAnnotation(plotobj)

      psfname = sprintf('%s_%s_region%d_ilev%d_evolution_%dexp', ...
                obsnames{ivar}, plotobj{1}.copystring, iregion, p.Results.level, NumExp);
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
   priornames{i} = sprintf('%s_guess',varnames{i});
   postenames{i} = sprintf('%s_analy',varnames{i});
end

for i = 1:nexp

   varexist(filenames{i}, {priornames{:}, postenames{:}, 'time', 'time_bounds'})

   diminfo = nc_getdiminfo(filenames{i},    'copy'); ncopies   = diminfo.Length;
   diminfo = nc_getdiminfo(filenames{i},'obstypes'); nobstypes = diminfo.Length;
   diminfo = nc_getdiminfo(filenames{i},  'region'); nregions  = diminfo.Length;

   commondata{i}.ncopies   = ncopies;
   commondata{i}.nobstypes = nobstypes;
   commondata{i}.nregions  = nregions;
   commondata{i}.times     = nc_varget(filenames{i},'time');
   commondata{i}.time_bnds = nc_varget(filenames{i},'time_bounds');
   commondata{i}.copyindex = get_copy_index(filenames{i},copystring);
   commondata{i}.lonlim1   = nc_attget(filenames{i},nc_global,'lonlim1');
   commondata{i}.lonlim2   = nc_attget(filenames{i},nc_global,'lonlim2');
   commondata{i}.latlim1   = nc_read_att(filenames{i}, nc_global,'latlim1');
   commondata{i}.latlim2   = nc_read_att(filenames{i}, nc_global,'latlim2');

   commondata{i}.region_names = nc_varget(filenames{i},'region_names');

   if (commondata{i}.nregions == 1 && (size(commondata{i}.region_names,2) == 1) )
      commondata{i}.region_names = deblank(commondata{i}.region_names');
   end

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

end

if mystat > 0
   error('The experiments are not compatible ... stopping.')
end

common = commondata{1};


%=====================================================================


function plotdat = getvals(fname, varname, copystring, prpo, regionindex, levelindex )
%% Get the data for each experiment
if (exist(fname,'file') ~= 2)
   error('%s does not exist',fname)
end

plotdat.fname         = fname;
plotdat.varname       = varname;
plotdat.copystring    = copystring;
plotdat.region        = regionindex;
plotdat.levelindex    = levelindex;
plotdat.bincenters    = nc_varget(fname,'time');
plotdat.binedges      = nc_varget(fname,'time_bounds');
plotdat.mlevel        = local_nc_varget(fname,'mlevel');
plotdat.plevel        = local_nc_varget(fname,'plevel');
plotdat.plevel_edges  = local_nc_varget(fname,'plevel_edges');
plotdat.hlevel        = local_nc_varget(fname,'hlevel');
plotdat.hlevel_edges  = local_nc_varget(fname,'hlevel_edges');
plotdat.ncopies       = length(nc_varget(fname,'copy'));

dimensionality        = nc_read_att(fname, nc_global, 'LocationRank');
plotdat.biasconv      = nc_attget(fname, nc_global, 'bias_convention');
plotdat.binseparation = nc_attget(fname, nc_global, 'bin_separation');
plotdat.binwidth      = nc_attget(fname, nc_global, 'bin_width');
plotdat.lonlim1       = nc_attget(fname, nc_global, 'lonlim1');
plotdat.lonlim2       = nc_attget(fname, nc_global, 'lonlim2');
plotdat.latlim1       = nc_read_att(fname, nc_global, 'latlim1');
plotdat.latlim2       = nc_read_att(fname, nc_global, 'latlim2');

% Coordinate between time types and dates

timeunits             = nc_attget(fname,'time','units');
calendar              = nc_attget(fname,'time','calendar');
timebase              = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin            = datenum(timebase(1),timebase(2),timebase(3));

plotdat.bincenters    = plotdat.bincenters + timeorigin;
plotdat.binedges      = plotdat.binedges   + timeorigin;
plotdat.Nbins         = length(plotdat.bincenters);

plotdat.timespan      = sprintf('%s through %s', ...
   datestr(min(plotdat.binedges(:))), ...
   datestr(max(plotdat.binedges(:))));

% Get the right indices for the intended variable, regardless of the storage order
% as well as some indices of other quantities of interest for future use.

plotdat.copyindex     = get_copy_index(fname, copystring);
plotdat.Npossindex    = get_copy_index(fname, 'Nposs');
plotdat.Nusedindex    = get_copy_index(fname, 'Nused');
plotdat.NQC4index     = get_copy_index(fname, 'N_DARTqc_4');
plotdat.NQC5index     = get_copy_index(fname, 'N_DARTqc_5');
plotdat.NQC6index     = get_copy_index(fname, 'N_DARTqc_6');
plotdat.NQC7index     = get_copy_index(fname, 'N_DARTqc_7');

plotdat.priorvar      = sprintf('%s_guess',plotdat.varname);
plotdat.postevar      = sprintf('%s_analy',plotdat.varname);

plotdat.trusted       = nc_read_att(fname, plotdat.priorvar, 'TRUSTED');
if (isempty(plotdat.trusted)), plotdat.trusted = 'NO'; end

myinfo.diagn_file     = fname;
myinfo.copyindex      = plotdat.copyindex;
myinfo.regionindex    = plotdat.region;
myinfo.levelindex     = plotdat.levelindex;

% get appropriate vertical coordinate variable

guessdims = nc_var_dims(fname, plotdat.priorvar);
analydims = nc_var_dims(fname, plotdat.postevar);

if ( dimensionality == 1 ) % observations on a unit circle, no level
   plotdat.level = 1;
   plotdat.level_units = [];
elseif ( strfind(guessdims{3},'surface') > 0 )
   plotdat.level       = 1;
   plotdat.level_units = 'surface';
   plotdat.level_edges = [];
elseif ( strfind(guessdims{3},'undef') > 0 )
   plotdat.level       = 1;
   plotdat.level_units = 'undefined';
   plotdat.level_edges = [];
else
   plotdat.level       = nc_varget(fname, guessdims{3});
   plotdat.level_units = nc_attget(fname, guessdims{3}, 'units');
   plotdat.level_edges = nc_varget(fname,sprintf('%s_edges',guessdims{3}));
end

[start, count]        = GetNCindices(myinfo,'diagn',plotdat.priorvar);
plotdat.prior         = nc_varget(fname, plotdat.priorvar, start, count);

[start, count]        = GetNCindices(myinfo,'diagn',plotdat.postevar);
plotdat.poste         = nc_varget(fname, plotdat.postevar, start, count);

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
      plotdat.ylabel = sprintf('%s (%s)',copystring, plotdat.biasconv);
   otherwise
      plotdat.Drange = [min(bob) max(bob)];
      plotdat.ylabel = copystring;
end

%% Get the number of observations possible and the number used.
%  N_DARTqc_5 is the number ignored because of namelist control.
%  N_DARTqc_6 is the number ignored because of incoming QC values.
%  It doesn't matter which prior/poste variable you get this information
%  from - they are both the same.

myinfo.diagn_file = fname;
myinfo.copyindex  = plotdat.Npossindex;
myinfo.levelindex = plotdat.levelindex;
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
   myinfo.copyindex = plotdat.Nusedindex;
   [start, count]   = GetNCindices(myinfo,'diagn',plotdat.priorvar);
   plotdat.nused    = nc_varget(fname, plotdat.priorvar, start, count);
else
   myinfo.copyindex = plotdat.Nusedindex;
   [start, count]   = GetNCindices(myinfo,'diagn',plotdat.postevar);
   plotdat.nused    = nc_varget(fname, plotdat.postevar, start, count);
end

%% Set the last of the ranges

plotdat.Nrange = [min(plotdat.nused(:))    max(plotdat.nposs(:))];


%=====================================================================


function myplot(plotobj, figdata)
%% myplot Creates a graphic for one region

Nexp    = length(plotobj);

%% Create the background

ax1   = subplot('position',figdata.position);
set(ax1,'YAxisLocation','left','FontSize',figdata.fontsize)

%% draw the results of the experiments, priors and posteriors
%  each with their own line type.
iexp   = 0;
hd     = [];   % handle to an unknown number of data lines
legstr = {[]}; % strings for the legend

for i = 1:Nexp
   if ( plotobj{i}.useprior )
      iexp         = iexp + 1;
      hd(iexp)     = line(plotobj{i}.bincenters, plotobj{i}.prior, ...
         'Color',    figdata.expcolors{i}, ...
         'Marker',   figdata.expsymbols{i}, ...
         'LineStyle',figdata.prpolines{1}, ...
         'LineWidth', figdata.linewidth,'Parent',ax1);
      legstr{iexp} = sprintf('%s Prior',plotobj{i}.title);
   end

   if ( plotobj{i}.useposterior )
      iexp         = iexp + 1;
      hd(iexp)     = line(plotobj{i}.bincenters, plotobj{i}.poste, ...
         'Color',    figdata.expcolors{i}, ...
         'Marker',   figdata.expsymbols{i}, ...
         'LineStyle',figdata.prpolines{2}, ...
         'LineWidth',figdata.linewidth,'Parent',ax1);
      legstr{iexp} = sprintf('%s Posterior',plotobj{i}.title);
   end
end

% Plot a bias line.

switch plotobj{1}.copystring
   case {'bias'}
      zeroline = line(get(ax1,'XLim'),[0 0],'Color',[0 100 0]/255,'Parent',ax1);
      set(zeroline,'LineWidth',2.5,'LineStyle','-')
   otherwise
end

% hokey effort to decide to plot months/days vs. daynum vs.
ttot = plotobj{1}.bincenters(plotobj{1}.Nbins) - plotobj{1}.bincenters(1) + 1;

if ((plotobj{1}.bincenters(1) > 1000) && (ttot > 5))
   datetick('x',6,'keeplimits','keepticks');
   monstr = datestr(plotobj{1}.bincenters(1),21);
   xlabelstring = sprintf('month/day - %s start',monstr);
elseif (plotobj{1}.bincenters(1) > 1000)
   datetick('x',15,'keeplimits','keepticks')
   monstr = datestr(plotobj{1}.bincenters(1),21);
   xlabelstring = sprintf('%s start',monstr);
else
   xlabelstring = 'days';
end

% Create another axes to use for plotting the observation counts

ax2 = axes('position',get(ax1,'Position'), ...
   'XAxisLocation','top', ...
   'YAxisLocation','right', ...
   'Color','none', ...
   'XColor',get(ax1,'XColor'), ...
   'YColor','b', ...
   'XLim',get(ax1,'XLim'), ...
   'YDir',get(ax1,'YDir'), ...
   'FontSize',get(ax1,'FontSize'));

% Plot the data, which sets the range of the axis
for i = 1:Nexp
   h2 = line(plotobj{i}.bincenters, plotobj{i}.nposs, ...
      'Color',figdata.expcolors{i},'Parent',ax2);
   h3 = line(plotobj{i}.bincenters, plotobj{i}.nused, ...
      'Color',figdata.expcolors{i},'Parent',ax2);
   set(h2,'LineStyle','none','Marker','o','MarkerSize',10);
   set(h3,'LineStyle','none','Marker','*','MarkerSize',10);
end

% use same X ticks
% use the same Y ticks, but find the right label values
set(ax2,'XTick', get(ax1,'XTick'), 'XTicklabel', []);
matchingYticks(ax1,ax2);

% Annotate. Trying to maximize content, minimize clutter.
annotate( ax1, ax2, plotobj{1}, figdata)

lh = legend(hd,legstr);
set(lh,'Interpreter','none','Box','off');

% The legend linesizes should match - 2 is hardwired - suprises me.

set(lh,'FontSize',figdata.fontsize);
kids = get(lh,'Children');
set(kids,'LineWidth',figdata.linewidth);


%=====================================================================


function annotate(ax1, ax2, plotobj, figdata)
%% One figure ... everything gets annotated.

set(get(ax1,'Xlabel'),'String',plotobj.timespan, ...
   'Interpreter','none','FontSize',figdata.fontsize)

if ( plotobj.useprior )
   ylabel = sprintf('forecast %s',plotobj.ylabel);
else
   ylabel = sprintf('analysis %s',plotobj.ylabel);
end

set(get(ax1,'Ylabel'),'String',ylabel, ...
    'Interpreter','none','FontSize',figdata.fontsize)
set(get(ax2,'Ylabel'),'String','# of obs (o=possible, \ast=assimilated)', ...
   'FontSize',figdata.fontsize)

if ( isempty(plotobj.level_units) )
   th = title({deblank(plotobj.region_names(plotobj.region,:)), ...
               sprintf('%s @ %d', plotobj.varname, plotobj.level(plotobj.levelindex))});
else
   th = title({deblank(plotobj.region_names(plotobj.region,:)), ...
               sprintf('%s @ %d %s', plotobj.varname, plotobj.level(plotobj.levelindex), ...
               plotobj.level_units)});
end
set(th,'Interpreter','none','FontSize',figdata.fontsize,'FontWeight','bold');


%=====================================================================


function BottomAnnotation(plotstruct)
%% annotates the filename containing the data being plotted

Nexp    = length(plotstruct);
dy      = 1.0/(Nexp+2);
yheight = 0.0225*(Nexp+1);

subplot('position',[0.10 0.01 0.8 yheight])
axis off

for ifile = 1:Nexp
   main = plotstruct{ifile}.fname;

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

   ty = 1.0 - (ifile+1)*dy;
   h = text(0.0, ty, string1);
   set(h, 'Interpreter', 'none', ...
          'HorizontalAlignment','left', ...
          'FontSize', 8);
      
   % If the observation is trusted for this experiment, annotate as such.

   switch lower(plotstruct{ifile}.trusted)
      case 'true'
         ty = 1.0 - Nexp*dy + ifile*dy*2;
         h = text(0.0,ty,sprintf('TRUSTED OBSERVATION in %s',plotstruct{ifile}.title));
         set(h,'FontSize',20,'Interpreter','none','HorizontalAlignment','left')
      otherwise
   end
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
%% try to set the axes into nicer-looking sizes with room for annotation.
%  figure out a page layout
%  extra space at the bottom for the date/file annotation
%  extra space at the top because the titles have multiple lines

ybot = 0.06 + nexp*0.075;  % room for dates/files
ytop = 0.125;              % room for title (always 2 lines)
dy   = 1.0 - ytop - ybot;
orientation   = 'landscape';
fontsize      = 16;
position      = [0.10 ybot 0.8 dy];
linewidth     = 2.0;

figdata = struct('expcolors',  {{'k','r','b','m','g','c','y'}}, ...
   'expsymbols', {{'o','s','d','p','h','s','*'}}, ...
   'prpolines',  {{'-','--'}}, 'position', position, ...
   'fontsize',fontsize, 'orientation',orientation, ...
   'linewidth',linewidth);


%=====================================================================


function value = local_nc_varget(fname,varname)
%% If the variable exists in the file, return the contents of the variable.
% if the variable does not exist, return empty value instead of error-ing
% out.

[variable_present, varid] = nc_var_exists(fname,varname);
if (variable_present)
   ncid  = netcdf.open(fname,'NOWRITE');
   value = netcdf.getVar(ncid, varid);
   netcdf.close(ncid)
else
   value = [];
end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

