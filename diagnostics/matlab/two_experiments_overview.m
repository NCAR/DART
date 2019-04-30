function two_experiments_overview(file1,label1,file2,label2,varargin)
%% two_experiments_overview create a trellis-plot-like view of the comparison of two experiments.
%
% EXAMPLE :
%
% file1 = 'Diags_2010.08.15-31_0-500m/obs_diag_output.nc';
% file2 = 'Diags_2010.08.15-31_Fixed_0-500m/obs_diag_output.nc';
% label1 = 'Experiment 1';
% label2 = 'Experiment 2';
% two_experiments_overview(file1, label1, file2, label2)
%
% EXAMPLE : specifying a 'flag level' to indicate 'weak' areas:
% In this case any region/level with only 10% of the max obs is 'weak'.
%
% two_experiments_overview(file1,label1,file2,label2,'FlagLevel',0.10)
%
% verticalobs = { ...
%     'RADIOSONDE_U_WIND_COMPONENT', 'RADIOSONDE_V_WIND_COMPONENT', ...
%     'RADIOSONDE_TEMPERATURE',      'RADIOSONDE_SPECIFIC_HUMIDITY', ...
%     'ACARS_U_WIND_COMPONENT',      'ACARS_V_WIND_COMPONENT', ...
%     'ACARS_TEMPERATURE',           'GPSRO_REFRACTIVITY', ...
%     'AIRCRAFT_U_WIND_COMPONENT',   'AIRCRAFT_V_WIND_COMPONENT', ...
%     'AIRCRAFT_TEMPERATURE',        'SAT_HORIZONTAL_WIND'};
%
% two_experiments_overview(file1,label1,file2,label2,'obsnames',verticalobs)

% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%%--------------------------------------------------------------------
% Decode,Parse,Check the input
%---------------------------------------------------------------------

defaultFlagLevel   = 0.00;   % ten percent would be 0.10
defaultVarCheck    = -1;
default_verbosity  = true;
default_markersize = 12;
default_pause      = false;
default_range      = [NaN NaN];
default_obsnames   = {'all'};
default_papertype  = 'uslegal';

p = inputParser;

addRequired(p,'file1' , @ischar);
addRequired(p,'label1', @ischar);
addRequired(p,'file2' , @ischar);
addRequired(p,'label2', @ischar);

if (exist('inputParser/addParameter','file') == 2)
    addParameter(p, 'FlagLevel',  defaultFlagLevel,   @isnumeric);
    addParameter(p, 'VarCheck',   defaultVarCheck,    @isnumeric);
    addParameter(p, 'verbose',    default_verbosity,  @islogical);
    addParameter(p, 'MarkerSize', default_markersize, @isnumeric);
    addParameter(p, 'pause',      default_pause,      @islogical);
    addParameter(p, 'range',      default_range,      @isnumeric);
    addParameter(p, 'obsnames',   default_obsnames,   @iscell);
    addParameter(p, 'papertype',  default_papertype,  @ischar);
else
    addParamValue(p, 'FlagLevel', defaultFlagLevel,   @isnumeric); %#ok<NVREPL>
    addParamValue(p, 'VarCheck',  defaultVarCheck,    @isnumeric); %#ok<NVREPL>
    addParamValue(p, 'verbose',   default_verbosity,  @islogical); %#ok<NVREPL>
    addParamValue(p, 'MarkerSize',default_markersize, @isnumeric); %#ok<NVREPL>
    addParamValue(p, 'pause',     default_pause,      @islogical); %#ok<NVREPL>
    addParamValue(p, 'range',     default_range,      @isnumeric); %#ok<NVREPL>
    addParamValue(p, 'obsnames',  default_obsnames,   @iscell);    %#ok<NVREPL>
    addParamValue(p, 'papertype', default_papertype,  @ischar);    %#ok<NVREPL>
end

p.parse(file1,label1,file2,label2,varargin{:}) % parse inputs

global figuredata
figuredata.MarkerSize = p.Results.MarkerSize;
figuredata.range      = p.Results.range;
figuredata.verbose    = p.Results.verbose;

%% collect the results of parsing (makes code easier to read)

FileA     = p.Results.file1;
FileB     = p.Results.file2;
LabelA    = p.Results.label1;
LabelB    = p.Results.label2;
VarCheck  = p.Results.VarCheck;
FlagLevel = p.Results.FlagLevel;

if (exist(FileA,'file') ~= 2), error('File %s does not exist.',FileA); end
if (exist(FileB,'file') ~= 2), error('File %s does not exist.',FileB); end

%% prepare the figures
% These positions fit 12 variables on a single page.

f1 = gcf;    set_figure(f1, p.Results.papertype);
f2 = figure; set_figure(f2, p.Results.papertype);
f3 = figure; set_figure(f3, p.Results.papertype);
f4 = figure; set_figure(f4, p.Results.papertype);

dx = 0.2000;
dy = 0.2400;

x1 = 0.0300;
x2 = x1 + 0.2500;
x3 = x2 + 0.2500;
x4 = x3 + 0.2500;

positions = [ x1 0.7000 dx dy; ...
    x2 0.7000 dx dy; ...
    x3 0.7000 dx dy; ...
    x4 0.7000 dx dy; ...
    x1 0.3700 dx dy; ...
    x2 0.3700 dx dy; ...
    x3 0.3700 dx dy; ...
    x4 0.3700 dx dy; ...
    x1 0.0400 dx dy; ...
    x2 0.0400 dx dy; ...
    x3 0.0400 dx dy; ...
    x4 0.0400 dx dy];
ntiles = size(positions,1);

bluered = flipud(redblue(15));
figure(f1); colormap(bluered)
figure(f2); colormap(bluered)

%% Create the list of vertical profile prior observation types in both files.

verticalobs = parse_DART_vars(FileA, FileB, p.Results.obsnames);

nvariables = length(verticalobs);

if nvariables <= ntiles
    printformat = '-dpdf';
else
    printformat = '-dpsc';
end

%% plot some reference plot just to make sure we're not upside down or ...

if ( VarCheck > nvariables )
    fprintf('\nThere are only %d possible variables in %s\n',nvariables,FileA)
    for ivar=1:nvariables
        fprintf('%40s is VarCheck %d\n',verticalobs{ivar},ivar)
    end
    error('VarCheck must be less than %d.',nvariables)
    
elseif ( VarCheck > 0 )
    close all
    files  = { FileA,  FileB};
    titles = {LabelA, LabelB};
    obsnames{1} = verticalobs{VarCheck};
    copy = 'bias';
    prpo = 'forecast';
    % fires up N figure windows ... for each region
    two_experiments_profile(files, titles, obsnames, copy, prpo)
    figure
end

%% plot the new stuff

ncidA = netcdf.open(FileA,'NOWRITE');
ncidB = netcdf.open(FileB,'NOWRITE');

pressure_dimid = netcdf.inqDimID(ncidA,'plevel');

plevels = ncread(FileA,'plevel');
hlevels = ncread(FileA,'hlevel');
regions = ncread(FileA,'region_names')';
nregions = size(regions,1);

region_name = cell(1,nregions);
for iregion=1:nregions
    region_name{iregion} = deblank(regions(iregion,:));
end

biasindex = get_copy_index(FileA,'bias');
rmseindex = get_copy_index(FileA,'rmse');
nusedindx = get_copy_index(FileA,'Nused');

plotnum = 0;
firstpage = true;
for ivar = 1:nvariables
    
    varname = sprintf('%s_VPguess',verticalobs{ivar});
    datasetA = ncread(FileA, varname);
    datasetB = ncread(FileB, varname);
    
    % datasetA,B are shaped  (region-level-copy)
    
    fprintf('Working on %s\n',varname)
    plotnum = plotnum + 1;
    
    diffs   = datasetA - datasetB;
    rmsemat = diffs(:,:,rmseindex)';
    nummat  = diffs(:,:,nusedindx)';
    biasmat = (abs(datasetA(:,:,biasindex)) - abs(datasetB(:,:,biasindex)))';
    
    % ignore the levels where there are no observations in the first dataset
    
    oldmask = (datasetA(:,:,nusedindx) < 1)';
    rmsemat(oldmask) = 0.0;
    biasmat(oldmask) = 0.0;
    
    varid = netcdf.inqVarID(ncidA,varname);
    [~, ~, dimids, ~] = netcdf.inqVar(ncidA,varid);
    
    if (any(dimids == pressure_dimid))
        % fprintf('Pressure vertical coords for %s\n',varname)
        levels = plevels;
    else
        % fprintf('Height vertical coords for %s\n',varname)
        levels = hlevels;
    end
    nlevels = length(levels);
    
    %----------------------------------------------------------------------
    % Plot the number of observations used and create the mask of which
    % cells are weakest. FIXME ... remove color from these cells.
    
    figure(f3);
    ha = axes('position',positions(plotnum,:));
    
    [hall, legendstr, weak] = obsplot(datasetA(:,:,nusedindx)', datasetB(:,:,nusedindx)', ...
        verticalobs{ivar}, levels, region_name, '# observations used', FlagLevel);
    
    if(ivar == nvariables)
        L = legend(hall,legendstr);
        xlabel('max observations per region')
        set(L,'Location','Best')
    end
    
    %----------------------------------------------------------------------
    % Plot the rmse
    
    figure(f1);
    axes('position',positions(plotnum,:))
    string1 = sprintf('rmse(%s)-rmse(%s)',LabelA,LabelB);
    myplot(rmsemat, verticalobs{ivar}, levels, region_name, string1, weak, LabelB)
    
    %----------------------------------------------------------------------
    % Plot the bias
    
    figure(f2);
    axes('position',positions(plotnum,:))
    string1 = sprintf('abs(bias(%s))-abs(bias(%s))',LabelA,LabelB);
    myplot(biasmat, verticalobs{ivar}, levels, region_name, string1, weak, LabelB)
    
    %----------------------------------------------------------------------
    % Plot the difference of the number of observations used
    % also spew these to the command line for context ...
    
    figure(f4);
    axes('position',positions(plotnum,:));
    
    %  offset each region, would need to scale for this to make sense.
    %     nummat   = nummat + ones(nlevels,1) * [0:nregions-1];
    
    line(nummat, 1:nlevels, 'LineWidth', 3.0);
    axis([-Inf Inf 0 nlevels])
    
    % Just convert the nuisance scientific notation to regular numbers
    xtick = get(gca,'XTick');
    xticklabel = cell(1,numel(xtick));
    for ix = 1:numel(xtick)
        xticklabel{ix} = sprintf('%d',xtick(ix));
    end
    
    yticklabel = cell(1,nlevels);
    for i=1:length(levels)
        yticklabel{i} = sprintf('%d',round(levels(i)));
    end
    
    set(gca,'YDir','normal', ...
        'XTick',xtick,'XTickLabel',xticklabel, ...
        'YTick',1:length(levels),'YTickLabel',yticklabel);
    
    string1 = sprintf('obs used; %s - %s',LabelA,LabelB);
    h = title({verticalobs{ivar},string1});
    set(h,'interpreter','none');
    grid on
    set(gca,'GridColor',[0 0 0],'GridAlpha',0.30);  % Darken the grid lines
    
    if(ivar == nvariables)
        L = legend(region_name);
        set(L,'Location','Best');
        %         set(ha,'Visible','off')
        %         set(hp,'Visible','off')
    end
    
    if figuredata.verbose
        fprintf('%s observation counts\n',LabelA)
        datasetA(:,:,nusedindx)'
        fprintf('%s observation counts\n',LabelB)
        datasetB(:,:,nusedindx)'
        fprintf('difference of observation counts\n')
        diffs(:,:,nusedindx)'
    end
    
    %----------------------------------------------------------------------
    % determine if we need to print another page
        
    if (ivar == ntiles) || (ivar == nvariables)
        if firstpage
            print(f1, printformat, 'rmse_diff')
            print(f2, printformat, 'abs_bias_diff')
            print(f3, printformat, 'num_obs_used')
            print(f4, printformat, 'num_obs_used_diff')
        else
            % can only get here with printformat = dpsc
            print(f1, printformat, '-append', 'rmse_diff')
            print(f2, printformat, '-append', 'abs_bias_diff')
            print(f3, printformat, '-append', 'num_obs_used')
            print(f4, printformat, '-append', 'num_obs_used_diff')
        end

        plotnum = 0;
        firstpage = false;
        if ivar ~= nvariables
            if p.Results.pause
                disp('pausing, hit any key to continue ...')
                pause
            end
            clf(f1)
            clf(f2)
            clf(f3)
            clf(f4)
        end
    end
    
end

netcdf.close(ncidA)
netcdf.close(ncidB)

end

%%----------------------------------------------------------------------

function vars = parse_DART_vars(file1, file2, verticalobs)

% rips off the _guess _analy _VPguess _VPanaly, and finds the unique ones,

if strncmpi(verticalobs{1},'all',3)
    obstype    = sort(get_DARTvars(file1));
    replace    = '';
    expression = '_VPguess';
    obstype    = regexprep(obstype,expression,replace);
    allvars    = sort(unique(obstype));
else
    allvars    = verticalobs;
end

% see if the _VPguess version of the variable exists in both files.

info1 = ncinfo(file1);
info2 = ncinfo(file2);

mask = false(1,length(allvars));

for ivar = 1:length(allvars)
    
    varname = sprintf('%s_VPguess',allvars{ivar});
    
    in1 = 0;
    for icandidate = 1:length(info1.Variables)
        in1 = in1 + strcmp(varname,info1.Variables(icandidate).Name);
    end
    
    in2 = 0;
    for icandidate = 1:length(info2.Variables)
        in2 = in2 + strcmp(varname,info2.Variables(icandidate).Name);
    end
    
    if (in1 > 0) && (in2 > 0)
        mask(ivar) = true;
    end
    
end

vars = allvars(mask);

end

%%----------------------------------------------------------------------

function myplot(datmat, varname, levels, regions, titlestring, weak, better)

% remove the color from the 'weak' items
% by using a symmetric colormap, we are forcing zero to be in the middle,
% and the colormap we are using has no color at zero.

global figuredata

[i,j] = find(weak > 0);
datmat(i(:),j(:)) = 0.0;

[ny, nx] = size(datmat);
xticks   = 1:nx;
yticks   = 1:ny;
imagesc(xticks,yticks,datmat)

% make a symmetric colormap ...
datmax = max(abs(datmat(:)));
if (datmax == 0.0)
    datmax = 1.0;
end

% xticklabel = {'North','Tropics','South'};
xticklabel = cell(1,length(regions));
for iregion=1:length(regions)
    bob = mysplit(regions{iregion});
    xticklabel{iregion} = bob{1};
end

yticklabel = cell(1,length(levels));
for ilevel=1:length(levels)
    yticklabel{ilevel} = sprintf('%d',round(levels(ilevel)));
end

set(gca,'YDir','normal', ...
    'XTick',xticks,'XTickLabel',xticklabel, ...
    'YTick',yticks,'YTickLabel',yticklabel, ...
    'CLim',[-datmax datmax]);

h = title({varname, titlestring});
set(h,'interpreter','none');
h = colorbar;
xh = get(h,'YLabel');
set(xh,'String',sprintf('positive means %s is better',better));

% plot a symbol at each of the weak region/levels

line(j,i,'LineStyle','none', ...
    'Marker','x', ...
    'MarkerSize',figuredata.MarkerSize, ...
    'MarkerFaceColor','k', ...
    'MarkerEdgeColor','k')

end

%%----------------------------------------------------------------------

function [hall, legendstr, weak] = obsplot(old, new, varname, levels, ...
    regions, titlestring, FlagLevel)

% first step is to identify the maximum value of each column/region
% since the known minimum is 0 (no observations).

nregions = length(regions);
nlevels  = length(levels);
maxs = max([old; new]);

% scale each column to have equal dynamic range
xticklabel = cell(1,nregions);
xticklabel{1} = ' ';
for iregion = 1:nregions
    old(:,iregion) = old(:,iregion) / maxs(iregion);
    new(:,iregion) = new(:,iregion) / maxs(iregion);
    xticklabel{iregion+1} = sprintf('%d',maxs(iregion));
end

yticklabel = cell(1,nlevels);
for i=1:nlevels
    yticklabel{i} = sprintf('%d',round(levels(i)));
end

% create a mask for the 'weak, sparsely-observed' cells
weak = (old < FlagLevel) | (new < FlagLevel);

% then offset from one another
shift = ones(nlevels,1) * [0:nregions-1];
old   = old + shift;
new   = new + shift;

h1 = line(old, 1:nlevels, 'LineWidth',1.0, 'Marker','*');
h2 = line(new, 1:nlevels, 'LineWidth',1.0, 'Marker','o');

axis([-Inf Inf 0 nlevels])

set(gca,'YDir','normal', ...
    'XTick',0:nregions,'XTickLabel',xticklabel, ...
    'YTick',1:nlevels,'YTickLabel',yticklabel);

h = title({varname,titlestring});
set(h,'interpreter','none');
grid on
set(gca,'GridColor',[0 0 0],'GridAlpha',0.30);  % Darken the grid lines
hall = [h1; h2];

legendstr = cell(1,nregions*2);
for iregion = 1:nregions
    legendstr{         iregion} = sprintf('Exp1 %s',regions{iregion});
    legendstr{nregions+iregion} = sprintf('Exp2 %s',regions{iregion});
end

end


%%----------------------------------------------------------------------


function colors = rwb()
% Creates a red-to-white-to-blue colormap without yellows.
%
% Example:
% bob = rwb;
% rgbplot(bob);
%
% -or, more simply-
%
% rgbplot(rwb)

num_points = 96;
mymean= 1.0;
mystd = 1.0;
x_min = mymean - 3*mystd;
x_max = mymean + 3*mystd;
x     = linspace(x_min, x_max, num_points);
e     = exp(1);
basen = (1.0 / (mystd * sqrt(2*pi)));
expon = -0.5 * (((x-mymean) / mystd).^2 );
y     = basen * (e .^ expon);

green = y./max(y);
red   = [1.0-green(1:2:end) ones(1,48)];
blue  = fliplr(red);

bob    = [red' green' blue'];
colors = bob(17:80,:);

end


%%----------------------------------------------------------------------


function c = redblue(m)

n = fix(m/2);
x = n~=(m/2);
r = [(0:1:n-1)/n,ones(1,n+x)];
g = [(0:1:n-1)/n,ones(1,x),(n-1:-1:0)/n];
b = [ones(1,n+x),(n-1:-1:0)/n];
c = [r(:),g(:),b(:)];

end


%%----------------------------------------------------------------------


function firstword = mysplit(inputstring)
% Because Mathworks refuses to leave well enough alone ...

existisstupid = which('strsplit');
if (isempty(existisstupid))
    firstword = {strtok(inputstring)};
else
    firstword = strsplit(inputstring);
end

end


%%----------------------------------------------------------------------


function set_figure(fignum,paper)

clf(fignum);
orient(fignum,'landscape');

set(fignum,'Units','inches');
set(fignum,'PaperPositionMode','auto');

switch lower(paper)
    case 'uslegal'
        disp('Setting paper to US Legal')
        set(fignum,'PaperType','uslegal','PaperUnits','inches');
        set(fignum,'PaperPosition',[0.25 0.25 13.5 8.0]);
        set(fignum,'Position',[0.25 0.25 13.5 8.0]);
    otherwise
        disp('Setting paper to US Letter')
        set(fignum,'PaperType','usletter','PaperUnits','inches');
        set(fignum,'PaperPosition',[0.25 0.25 10.5 8.0]);
        set(fignum,'Position',[0.25 0.25 10.5 8.0]);
end

end


%%----------------------------------------------------------------------
% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
