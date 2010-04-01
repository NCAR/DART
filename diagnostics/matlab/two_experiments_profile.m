function two_experiments_profile(files, titles, varnames, qtty, prpo)
%
% Each variable gets its own figure.
% Each region gets its own axis.
% Multiple quantities (rmse, bias) can be plotted on same axis.
% 
% files = {'/fs/image/home/hliu/DART/models/wrf/work/ernesto/gpsonly/obs_diag_output.nc',
%          '/fs/image/home/hliu/DART/models/wrf/work/ernesto/ctl/obs_diag_output.nc'};
% titles   = {'GPS only', 'Control'};
% varnames = {'RADIOSONDE_U_WIND_COMPONENT', 'RADIOSONDE_TEMPERATURE'};
% qtty     = {'rmse','bias'};     % rmse, spread, totalspread, bias, etc.
% prpo     = 'analysis'; % [analy, analysis, posterior ] == posterior
% prpo     = 'forecast'; % [guess, forecast, prior     ] == prior
% prpo     = 'both'; % prior & posterior
%
% two_experiments_profile(files, titles, varnames, qtty, prpo)
%

% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

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

commondata = check_compatibility(files, varnames);

%%--------------------------------------------------------------------
% Set some static data
%---------------------------------------------------------------------

nvars = length(varnames);
nqtty = length(qtty);


for ivar = 1:nvars
   fprintf('Working on %s ...\n',varnames{ivar})

   % FIXME might want to set up a loop over regions here ...

   % set up all the stuff that is common. 

   %---------------------------------------------------------------------
   % Getting the data for each experiment
   %---------------------------------------------------------------------

   Nlimits = zeros(NumExp,2);  % range of observation count - min, then max
   Dlimits = zeros(NumExp,2);  % range of the data
   Ylimits = zeros(NumExp,2);  % range of the vertical coords

   for iexp = 1:NumExp

      plotobj{iexp} = getvals(files{iexp}, varnames{ivar}, qtty{1}, prpo);
      plotobj{iexp}.title = titles{iexp};

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

   myplot( plotobj, Nrange, Drange, Yrange);
   
   if ( ivar ~= nvars )
      disp('Pausing, hit any key to continue ...')
      pause
   end

end  % of loop around variable

%=====================================================================
% End of main function body. Helper functions below.
%=====================================================================



function commondata = check_compatibility(filenames, varnames)
%----------------------------------------------------------------------
% Trying to prevent the comparison of apples and oranges.
% make sure the diagnostics were generated the same way.

% need to check that the variables exist in all files

% need to check that the copies exist in all files

% need to check that the timeframe is the same for all files

% need to check that the region definitions are the same for all files 

% nice to check that the number of possible observations is the same 

nexp = length(filenames);
commondata = struct('ncopies',   zeros(1,nexp), ...
                    'nobstypes', zeros(1,nexp), ...
                    'nregions',  zeros(1,nexp), ...
                    'times',     zeros(1,nexp));

for i = 1,length(filenames)
   diminfo   = nc_getdiminfo(filenames{i},    'copy'); ncopies   = diminfo.Length;  
   diminfo   = nc_getdiminfo(filenames{i},'obstypes'); nobstypes = diminfo.Length;  
   diminfo   = nc_getdiminfo(filenames{i},  'region'); nregions  = diminfo.Length;  
   times     = nc_varget(filenames{i},'time');
   time_bnds = nc_varget(filenames{i},'time_bounds');

   commondata.ncopies(i)   = ncopies;
   commondata.nobstypes(i) = nobstypes;
   commondata.nregions(i)  = nregions;

   % FIXME check to make sure variables exist in the files

end

mystat = 0;

% FIXME ... more realistic error checking, please
if ( 1 == 2 )
   fprintf('There are different numbers of copies in the experiments.\n')
%   fprintf('one experiment had %d, the other had %d\n',copyA,copyB)
   mystat = 1;
end



function plotdat = getvals(fname, varname, copystring, prpo )
%----------------------------------------------------------------------
if (exist(fname,'file') ~= 2)
   error('%s does not exist',fname)
end

plotdat.fname         = fname;
plotdat.varname       = varname;  
plotdat.copystring    = copystring;
plotdat.bincenters    = nc_varget(fname,'time');
plotdat.binedges      = nc_varget(fname,'time_bounds');
plotdat.region_names  = nc_varget(fname,'region_names');
plotdat.nregions      = size(plotdat.region_names,1);
plotdat.binseparation = nc_attget(fname,nc_global,'bin_separation');
plotdat.binwidth      = nc_attget(fname,nc_global,'bin_width');
time_to_skip          = nc_attget(fname,nc_global,'time_to_skip');
plotdat.lonlim1       = nc_attget(fname,nc_global,'lonlim1');
plotdat.lonlim2       = nc_attget(fname,nc_global,'lonlim2');
plotdat.latlim1       = nc_attget(fname,nc_global,'latlim1');
plotdat.latlim2       = nc_attget(fname,nc_global,'latlim2');
plotdat.biasconv      = nc_attget(fname,nc_global,'bias_convention');

% Coordinate between time types and dates

timeunits             = nc_attget(fname,'time','units');
calendar              = nc_attget(fname,'time','calendar');
timebase              = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin            = datenum(timebase(1),timebase(2),timebase(3));
skip_seconds          = time_to_skip(4)*3600 + time_to_skip(5)*60 + time_to_skip(6);
iskip                 = time_to_skip(3) + skip_seconds/86400;

plotdat.bincenters    = plotdat.bincenters + timeorigin;
plotdat.binedges      = plotdat.binedges   + timeorigin;
plotdat.Nbins         = length(plotdat.bincenters);
plotdat.toff          = plotdat.bincenters(1) + iskip;

plotdat.timespan      = sprintf('%s through %s',  ...
                        datestr(min(plotdat.binedges(:))), ...
                        datestr(max(plotdat.binedges(:))));

% Get the right indices for the intended variable, regardless of the storage order

plotdat.copyindex = get_copy_index(fname, copystring);
plotdat.priorvar  = sprintf('%s_VPguess',plotdat.varname);
plotdat.postevar  = sprintf('%s_VPanaly',plotdat.varname);

myinfo.diagn_file = fname;
myinfo.copyindex  = plotdat.copyindex;
[start, count]    = GetNCindices(myinfo,'diagn',plotdat.priorvar);
count(count < 1)  = -1;
plotdat.prior     = nc_varget(fname, plotdat.priorvar, start, count);
plotdat.poste     = nc_varget(fname, plotdat.postevar, start, count);

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
   plotdat.Drange = [    0    max(bob)];
   plotdat.xlabel = {sprintf('%s %s',copystring, plotdat.biasconv), plotdat.timespan};
otherwise
   plotdat.Drange = [min(bob) max(bob)];
   plotdat.xlabel = {copystring, plotdat.timespan};
end

%% Get the right indices for the number of observations possible
%  Get the right indices for the number of observations used
%  FIXME - should the number rejected because of incoming QC be disqualified

myinfo.diagn_file = fname;
myinfo.copyindex  = get_copy_index(fname, 'Nposs');
[start, count]    = GetNCindices(myinfo,'diagn',plotdat.priorvar);
count(count < 1)  = -1;
plotdat.nposs     = nc_varget(fname, plotdat.priorvar, start, count);

if ( plotdat.useprior ) 
   myinfo.copyindex  = get_copy_index(fname, 'Nused');
   [start, count]    = GetNCindices(myinfo,'diagn',plotdat.priorvar);
   count(count < 1)  = -1;
   plotdat.nused     = nc_varget(fname, plotdat.priorvar, start, count);
else
   myinfo.copyindex  = get_copy_index(fname, 'Nused');
   [start, count]    = GetNCindices(myinfo,'diagn',plotdat.postevar);
   count(count < 1)  = -1;
   plotdat.nused     = nc_varget(fname, plotdat.postevar, start, count);
end

%% Set the last of the ranges

plotdat.Yrange = [min(plotdat.level_edges) max(plotdat.level_edges)];
plotdat.Nrange = [min(plotdat.nused(:))    max(plotdat.nposs(:))];



function h = Stripes(x,edges)
% EraseMode: [ {normal} | background | xor | none ]

hold on;

% plot two little dots at the corners to make Matlab
% choose some plot limits. Given those nice limits and
% tick labels ... KEEP THEM. Later, make the dots invisible.

h = plot([min(x) max(x)],[min(edges) max(edges)]);
axlims = axis;
axlims(4) = max(edges);
axlims(3) = -100;
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



function kids = myplot( plotobj, Nrange, Drange, Yrange)

orientation = 'tall';
fontsize    = 16;
expcolors   = {'k','r'};
prpolines   = {'-',':'};
expsymbols  = {'+','o'};

clf; orient(gcf, orientation) 

dx = 0.8;
plotlims = [0.10 0.15 dx 0.7;
            0.38 0.15 dx 0.7;
            0.66 0.15 dx 0.7];

Nexp = size(Nrange,2);

iplot = 0;
%for iregion = plotobj.regions
for iregion = 1
    
   %% Create the background stripes, etc.
   
   iplot = iplot + 1;
   ax1   = subplot('position',plotlims(iplot,:));
   hd    = zeros(1,2*Nexp);   % handle to the data lines
   
   Stripes(Drange, plotobj{1}.level_edges);
   set(ax1,'YDir',plotobj{1}.YDir,'YTick',sort(plotobj{1}.levels))
   set(ax1,'YAxisLocation','left')
   hold on
     
   %% draw the results of the experiments, priors and posteriors
   %  each with their own line type.
   iexp = 0;
   legstr = {[]};
   
   for i = 1:Nexp

      if ( plotobj{i}.useprior )
         iexp     = iexp + 1;
         lty = sprintf('%s%s%s',expcolors{i},prpolines{1},expsymbols{i});
         hd(iexp) = plot(plotobj{i}.prior, plotobj{i}.levels, lty,'LineWidth', 2.0);
         legstr{iexp} = sprintf('%s forecast',plotobj{i}.title);
      end
     
      if ( plotobj{i}.useposterior ) 
         iexp     = iexp + 1;
         lty = sprintf('%s%s%s',expcolors{i},prpolines{2},expsymbols{i});
         hd(iexp) = plot(plotobj{i}.poste, plotobj{i}.levels, lty,'LineWidth', 2.0);
         legstr{iexp} = sprintf('%s analysis',plotobj{i}.title);
      end
   end

   switch plotobj{1}.copystring
   case {'bias'}
         biasline = line([0 0],Yrange,'Color','k','Parent',ax1);
         set(biasline,'LineWidth',2.0,'LineStyle','-.')
   otherwise   % draw the results of the experiments, priors and posteriors
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

   % Plot the data, which sets the range of the axis
   for i = 1:Nexp
      h2 = line(plotobj{i}.nposs, plotobj{i}.levels,'Color','b','Parent',ax2);
      h3 = line(plotobj{i}.nused, plotobj{i}.levels,'Color',expcolors{i},'Parent',ax2);
      set(h2,'LineStyle','none','Marker','o');
      set(h3,'LineStyle','none','Marker','+');
      hn(i) = h3;
   end
   
   % use same Y ticks
   set(ax2,'YTick',     get(ax1,'YTick'), ...
           'YTicklabel',get(ax1,'YTicklabel'));
   set(get(ax1,'Ylabel'),'String',plotobj{i}.level_units,'Interpreter','none')
   set(get(ax2,'Ylabel'),'String',plotobj{i}.level_units,'Interpreter','none')

   % use the same X ticks, but find the right label values
   [xticks, newticklabels] = matchingXticks(ax1,ax2);
   set(ax2,'XTick', xticks, 'XTicklabel', newticklabels)

   set(get(ax2,'Xlabel'),'String','# of obs (o=poss, +=used)')
   set(get(ax1,'Xlabel'),'String',plotobj{i}.xlabel,'Interpreter','none')

   set(ax1,'Position',get(ax2,'Position'))

   % Annotate the whole thing

   th = title({deblank(plotobj{1}.region_names(iregion,:)), plotobj{1}.varname});
   set(th,'Interpreter','none','FontSize',fontsize);

   lh = legend(hd,legstr);     
   legend(lh,'boxoff');
   
   set(lh,'FontSize',16);
   kids = get(lh,'Children');
   if (length(kids) < 8)
      set(kids([2 5]),'LineWidth',2.0);
   else
      set(kids([2 5 8 11]),'LineWidth',2.0);
   end

end


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

