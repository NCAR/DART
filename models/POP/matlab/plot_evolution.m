function plotdat = plot_evolution(fname, copystring, varargin)
%% plot_evolution plots the temporal evolution of the observation-space quantities for all possible levels, all possible variables.
% Part of the observation-space diagnostics routines.
%
% 'obs_diag' produces a netcdf file containing the diagnostics.
% obs_diag condenses the obs_seq.final information into summaries for a few specified
% regions - on a level-by-level basis.
%
% USAGE: plotdat = plot_evolution(fname, copystring [,varargin]);
%
% fname         : netcdf file produced by 'obs_diag'
% copystring    : 'copy' string == quantity of interest. These
%                 can be any of the ones available in the netcdf
%                 file 'CopyMetaData' variable.
%                 (ncdump -v CopyMetaData obs_diag_output.nc)
% obstypestring : 'observation type' string. Optional.
%                 Must match something in the netcdf
%                 file 'ObservationTypes' variable.
%                 (ncdump -v ObservationTypes obs_diag_output.nc)
%                 If specified, only this observation type will be plotted.
%                 If not specified, all observation types incluced in the netCDF file
%                 will be plotted.
%
% OUTPUT: two files will result for each observation type plotted. One is a
%         postscript file containing a page for each level - all regions.
%         The other file is a simple text file containing summary information
%         about how many observations were assimilated, how many were available, etc.
%         Both of these filenames contain the observation type as part of the name.
%
%
% EXAMPLE 1 - plot the evolution of the bias for all observation types, all levels
%
% fname      = 'obs_diag_output.nc';   % netcdf file produced by 'obs_diag'
% copystring = 'bias';                 % 'copy' string == quantity of interest
% plotdat    = plot_evolution(fname, copystring);
%
%
% EXAMPLE 2 - plot the evolution of the rmse for just the mooring temperature observations
%             This requires that the 'MOORING_TEMPERATURE' is one of the known observation
%             types in the netCDF file.
%
% fname      = 'obs_diag_output.nc';   % netcdf file produced by 'obs_diag'
% copystring = 'rmse';                 % 'copy' string == quantity of interest
% plotdat    = plot_evolution(fname, copystring, 'MOORING_TEMPERATURE');

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if nargin == 2
   nvars = 0;
elseif nargin == 3
   varname = varargin{1};
   nvars = 1;
else
   error('wrong number of arguments ... ')
end

if (exist(fname,'file') ~= 2)
   error('file/fname <%s> does not exist',fname)
end

% Harvest plotting info/metadata from netcdf file.

plotdat.fname         = fname;
plotdat.copystring    = copystring;
plotdat.bincenters    = nc_varget(fname,'time');
plotdat.binedges      = nc_varget(fname,'time_bounds');
plotdat.mlevel        = local_nc_varget(fname,'mlevel');
plotdat.plevel        = local_nc_varget(fname,'plevel');
plotdat.plevel_edges  = local_nc_varget(fname,'plevel_edges');
plotdat.hlevel        = local_nc_varget(fname,'hlevel');
plotdat.hlevel_edges  = local_nc_varget(fname,'hlevel_edges');
plotdat.ncopies       = length(nc_varget(fname,'copy'));
plotdat.nregions      = length(nc_varget(fname,'region'));
plotdat.region_names  = nc_varget(fname,'region_names');

% Matlab wants character matrices to be Nx1 instead of 1xN.

if (plotdat.nregions == 1 && (size(plotdat.region_names,2) == 1) )
   plotdat.region_names = deblank(plotdat.region_names');
end

dimensionality        = local_nc_attget(fname, nc_global, 'LocationRank');
plotdat.binseparation = local_nc_attget(fname, nc_global, 'bin_separation');
plotdat.binwidth      = local_nc_attget(fname, nc_global, 'bin_width');
time_to_skip          = local_nc_attget(fname, nc_global, 'time_to_skip');
plotdat.lonlim1       = local_nc_attget(fname, nc_global, 'lonlim1');
plotdat.lonlim2       = local_nc_attget(fname, nc_global, 'lonlim2');
plotdat.latlim1       = local_nc_attget(fname, nc_global, 'latlim1');
plotdat.latlim2       = local_nc_attget(fname, nc_global, 'latlim2');
plotdat.biasconv      = local_nc_attget(fname, nc_global, 'bias_convention');

% Coordinate between time types and dates

calendar     = nc_attget(fname,'time','calendar');
timeunits    = nc_attget(fname,'time','units');
timebase     = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin   = datenum(timebase(1),timebase(2),timebase(3));
if ( isempty(time_to_skip) == 1)
   iskip = 0;
elseif ( numel(time_to_skip) == 6)
   skip_seconds = time_to_skip(4)*3600 + time_to_skip(5)*60 + time_to_skip(6);
   iskip        = time_to_skip(3) + skip_seconds/86400;
else
   error('time_to_skip variable has unusual length. Should be either 0 or 6.')
end

plotdat.bincenters = plotdat.bincenters + timeorigin;
plotdat.binedges   = plotdat.binedges   + timeorigin;
plotdat.Nbins      = length(plotdat.bincenters);
plotdat.toff       = plotdat.bincenters(1) + iskip;

% set up a structure with all static plotting components

plotdat.linewidth = 2.0;

if (nvars == 0)
   [plotdat.allvarnames, plotdat.allvardims] = get_varsNdims(fname);
   [plotdat.varnames,    plotdat.vardims]    = FindTemporalVars(plotdat);
   plotdat.nvars       = length(plotdat.varnames);
else
   plotdat.varnames{1} = varname;
   plotdat.nvars       = nvars;
end

plotdat.copyindex   = get_copy_index(fname,copystring);
plotdat.Npossindex  = get_copy_index(fname,'Nposs');
plotdat.Nusedindex  = get_copy_index(fname,'Nused');
plotdat.NQC4index   = get_copy_index(fname,'N_DARTqc_4');
plotdat.NQC5index   = get_copy_index(fname,'N_DARTqc_5');
plotdat.NQC6index   = get_copy_index(fname,'N_DARTqc_6');
plotdat.NQC7index   = get_copy_index(fname,'N_DARTqc_7');

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
   fprintf('Removing %s from the current directory.\n',psfname)
   system(sprintf('rm %s',psfname));

   % remove any existing log file -

   lgfname = sprintf('%s_%s_obscount.txt',plotdat.varnames{ivar},plotdat.copystring);
   fprintf('Removing %s from the current directory.\n',lgfname)
   system(sprintf('rm %s',lgfname));
   logfid = fopen(lgfname,'wt');
   fprintf(logfid,'%s\n',lgfname);

   % get appropriate vertical coordinate variable

   guessdims = nc_var_dims(fname, plotdat.guessvar);
   analydims = nc_var_dims(fname, plotdat.analyvar);

   if ( dimensionality == 1 ) % observations on a unit circle, no level
      plotdat.level = 1;
      plotdat.level_units = [];
   elseif ( strfind(guessdims{3},'surface') > 0 )
      plotdat.level       = 1;
      plotdat.level_units = 'surface';
   elseif ( strfind(guessdims{3},'undef') > 0 )
      plotdat.level       = 1;
      plotdat.level_units = 'undefined';
   else
      plotdat.level       = nc_varget(fname, guessdims{3});
      plotdat.level_units = nc_attget(fname, guessdims{3}, 'units');
   end
   plotdat.nlevels = length(plotdat.level);

   % Here is the tricky part. Singleton dimensions are auto-squeezed ...
   % single levels, single regions ...

   guess_raw = nc_varget(fname, plotdat.guessvar);
   guess = reshape(guess_raw, plotdat.Nbins,   plotdat.ncopies, ...
      plotdat.nlevels, plotdat.nregions);

   analy_raw = nc_varget(fname, plotdat.analyvar);
   analy = reshape(analy_raw, plotdat.Nbins,   plotdat.ncopies, ...
      plotdat.nlevels, plotdat.nregions);

   % check to see if there is anything to plot
   nextrap = sum(guess(:,plotdat.NQC4index  ,:,:)) + ...
             sum(guess(:,plotdat.NQC6index  ,:,:));
   nposs   = sum(guess(:,plotdat.Npossindex, :,:)) - nextrap;

   if ( sum(nposs(:)) < 1 )
      fprintf('%s no obs for %s...  skipping\n', plotdat.varnames{ivar})
      continue
   end

   for ilevel = 1:plotdat.nlevels

      fprintf(logfid,'\nlevel %d %f %s\n',ilevel,plotdat.level(ilevel),plotdat.level_units);
      plotdat.ges_Nqc4  = guess(:,plotdat.NQC4index  ,ilevel,:);
      plotdat.anl_Nqc4  = analy(:,plotdat.NQC4index  ,ilevel,:);
      fprintf(logfid,'DART QC == 4, prior/post %d %d\n',sum(plotdat.ges_Nqc4(:)), ...
         sum(plotdat.anl_Nqc4(:)));

      plotdat.ges_Nqc5  = guess(:,plotdat.NQC5index  ,ilevel,:);
      plotdat.anl_Nqc5  = analy(:,plotdat.NQC5index  ,ilevel,:);
      fprintf(logfid,'DART QC == 5, prior/post %d %d\n',sum(plotdat.ges_Nqc5(:)), ...
         sum(plotdat.anl_Nqc5(:)));

      plotdat.ges_Nqc6  = guess(:,plotdat.NQC6index  ,ilevel,:);
      plotdat.anl_Nqc6  = analy(:,plotdat.NQC6index  ,ilevel,:);
      fprintf(logfid,'DART QC == 6, prior/post %d %d\n',sum(plotdat.ges_Nqc6(:)), ...
         sum(plotdat.anl_Nqc6(:)));

      plotdat.ges_Nqc7  = guess(:,plotdat.NQC7index  ,ilevel,:);
      plotdat.anl_Nqc7  = analy(:,plotdat.NQC7index  ,ilevel,:);
      fprintf(logfid,'DART QC == 7, prior/post %d %d\n',sum(plotdat.ges_Nqc7(:)), ...
         sum(plotdat.anl_Nqc7(:)));

      plotdat.ges_Nposs = guess(:,plotdat.Npossindex, ilevel,:);
      plotdat.anl_Nposs = analy(:,plotdat.Npossindex, ilevel,:);
      fprintf(logfid,'# obs     ,   prior/post %d %d\n',sum(plotdat.ges_Nposs(:)), ...
                                                        sum(plotdat.anl_Nposs(:)));

      plotdat.ges_Nposs = plotdat.ges_Nposs - plotdat.ges_Nqc4 - plotdat.ges_Nqc6;
      plotdat.anl_Nposs = plotdat.anl_Nposs - plotdat.anl_Nqc4 - plotdat.anl_Nqc6;
      fprintf(logfid,'# obs poss,   prior/post %d %d\n',sum(plotdat.ges_Nposs(:)), ...
         sum(plotdat.anl_Nposs(:)));

      plotdat.ges_Nused = guess(:,plotdat.Nusedindex, ilevel,:);
      plotdat.anl_Nused = analy(:,plotdat.Nusedindex, ilevel,:);
      fprintf(logfid,'# obs used,   prior/post %d %d\n',sum(plotdat.ges_Nused(:)), ...
         sum(plotdat.anl_Nused(:)));

      plotdat.ges_copy  = guess(:,plotdat.copyindex,  ilevel,:);
      plotdat.anl_copy  = analy(:,plotdat.copyindex,  ilevel,:);

      plotdat.Yrange    = FindRange(plotdat);

      % plot by region

      if (plotdat.nregions > 2)
         clf; orient tall; wysiwyg
      else
         clf; orient landscape; wysiwyg
      end

      for iregion = 1:plotdat.nregions

         plotdat.region   = iregion;
         plotdat.myregion = deblank(plotdat.region_names(iregion,:));
         if ( isempty(plotdat.level_units) )
            plotdat.title    = plotdat.myvarname;
         else
            plotdat.title    = sprintf('%s @ %d %s',    ...
               plotdat.myvarname,     ...
               plotdat.level(ilevel), ...
               plotdat.level_units);
         end

         myplot(plotdat);
      end

      % create a postscript file
      print(gcf,'-dpsc','-append',psfname);

      % block to go slow and look at each one ...
      % disp('Pausing, hit any key to continue ...')
      % pause

   end
end

%=====================================================================
% 'Helper' functions
%=====================================================================


function myplot(plotdat)

%% The prior and posterior are plotted as separate items.
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

string_guess = sprintf('forecast: mean=%.5g', mean_prior);
string_analy = sprintf('analysis: mean=%.5g', mean_post);
plotdat.subtitle = sprintf('%s     %s     %s',plotdat.myregion, string_guess, string_analy);

% Plot the requested quantity on the left axis.
% The observation count will use the axis on the right.
% We want to suppress the 'auto' feature of the axis labelling,
% so we manually set some values that normally
% don't need to be set.

ax1 = subplot(plotdat.nregions,1,plotdat.region);

h1 = plot(tg,cg,'k+-',ta,ca,'ro-','LineWidth',plotdat.linewidth);
h  = legend('forecast', 'analysis');
legend(h,'boxoff','Interpreter','none')

% get the range of the existing axis and replace with
% replace y axis values
% reset the axes limits

axlims = axis;
axlims = [axlims(1:2) plotdat.Yrange];
axis(axlims)

plotdat.ylabel{1} = plotdat.myregion;
switch lower(plotdat.copystring)
   case 'bias'
      % plot a zero-bias line
      zeroline = line(axlims(1:2),[0 0], 'Color','r','Parent',ax1);
      set(zeroline,'LineWidth',1.0,'LineSTyle','-.')
      plotdat.ylabel{2} = sprintf('%s (%s)',plotdat.copystring,plotdat.biasconv);
   otherwise
      plotdat.ylabel{2} = sprintf('%s',plotdat.copystring);
end

% hokey effort to decide to plot months/days vs. daynum vs.
ttot = plotdat.bincenters(plotdat.Nbins) - plotdat.bincenters(1) + 1;

if ((plotdat.bincenters(1) > 1000) && (ttot > 5))
   datetick('x',6,'keeplimits','keepticks');
   monstr = datestr(plotdat.bincenters(1),21);
   xlabelstring = sprintf('month/day - %s start',monstr);
elseif (plotdat.bincenters(1) > 1000)
   datetick('x',15,'keeplimits','keepticks')
   monstr = datestr(plotdat.bincenters(1),21);
   xlabelstring = sprintf('%s start',monstr);
else
   xlabelstring = 'days';
end

% only put x axis on last/bottom plot
if (plotdat.region == plotdat.nregions)
   xlabel(xlabelstring)
end

% more annotation ...
if (plotdat.region == 1)
   title({plotdat.title, plotdat.subtitle}, ...
      'Interpreter', 'none', 'Fontsize', 12, 'FontWeight', 'bold')
   BottomAnnotation(plotdat.fname)
else
   title(plotdat.subtitle, ...
      'Interpreter', 'none', 'Fontsize', 12, 'FontWeight', 'bold')
end

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

% use same X ticks - with no labels
set(ax2,'XTick',get(ax1,'XTick'),'XTickLabel',[])

% use the same Y ticks, but find the right label values
% Make sure that values at ticks are whole numbers.
matchingYticks(ax1,ax2);

set(get(ax2,'Ylabel'),'String','# of obs : o=poss, +=used')
set(get(ax1,'Ylabel'),'String',plotdat.ylabel,'Interpreter','none')


%=====================================================================


function BottomAnnotation(main)
%% annotates the full path of the file being plotted
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


%=====================================================================


function [y,ydims] = FindTemporalVars(x)
%% Returns UNIQUE (i.e. base) temporal variable names
if ( ~(isfield(x,'allvarnames') && isfield(x,'allvardims')))
   error('Doh! no ''allvarnames'' and ''allvardims'' components')
end

j = 0;

for i = 1:length(x.allvarnames)
   indx = strfind(x.allvardims{i},'time');
   if (indx > 0)
      j = j + 1;

      basenames{j} = ReturnBase(x.allvarnames{i});
      basedims{ j} = x.allvardims{i};
   end
end

[~,i,j] = unique(basenames);
y     = cell(length(i),1);
ydims = cell(length(i),1);
for k = 1:length(i)
   fprintf('%2d is %s\n',k,basenames{i(k)})
   y{k}     = basenames{i(k)};
   ydims{k} = basedims{ i(k)};
end


%=====================================================================


function s = ReturnBase(string1)
%% Pick off the variable name.
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

bob  = [y.ges_copy(:) ; ...
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
   elseif ( ymax < 0.0 && strcmp(y.copystring,'bias') )
      ymax = 0.0;
   end

   Yrange = [ymin ymax];

   x = [min([Yrange(1) 0.0]) Yrange(2)];
end


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


%=====================================================================


function value = local_nc_attget(fname,varid,varname)
%% If the (global) attribute exists, return the value.
% If it does not, do not throw a hissy-fit.

value = [];
if (varid == nc_global)
   finfo = ncinfo(fname);
   for iatt = 1:length(finfo.Attributes)
      if (strcmp(finfo.Attributes(iatt).Name, deblank(varname)))
         value = finfo.Attributes(iatt).Value;
         return
      end
   end
else
   fprintf('function not supported for local variables, only global atts.\n')
end



% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
