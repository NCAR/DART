function two_experiments_overview(varargin)
% two_experiments_overview create a trellis-plot-like view of all variables.
% With no arguments, two_experiments_overview will create 4 figures that attempt
% to provide an overview of two experiments. Each experiment must have been
% summarized by obs_diag() and have its own 'obs_diag_output.nc' file.
% The 'old' experiment has a default name of 'old_obs_diag_output.nc'
% The 'new' experiment has a default name of 'new_obs_diag_output.nc'
%
%
% EXAMPLE using defaults:
%
% two_experiments_overview
%
%
% EXAMPLE specifying filenames:
%
% file1 = '/glade/scratch/raeder/POP_force/POP15/Diags_2010.08.15-31_0-500m/obs_diag_output.nc';
% file2 = '/glade/scratch/raeder/ATM_spinup2/Diags_2010.08.15-31_Fixed_0-500m/obs_diag_output.nc';
% two_experiments_overview('OldFile',file1,'NewFile',file2)
%
%
% EXAMPLE specifying filenames and a 'flag level' to indicate 'weak' areas:
%
% two_experiments_overview('OldFile',file1,'NewFile',file2,'FlagLevel',0.10)

%% DART software - Copyright 2004 - 2016 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

matver = ver('Matlab');

if str2double(matver.Version) < 8.3
    error('Need at least Matlab 2014a for this function.')
end

%% handle the user input - somehow.

p = inputParser;
p.FunctionName = 'input parser :: no required arguments, optional input is FlagLevel (percentage)';

% set defaults for optional parameters
defaultOldFile = 'old_obs_diag_output.nc';
defaultNewFile = 'new_obs_diag_output.nc';
defaultFlagLevel = 0.00;   % ten percent would be 0.10
defaultVarCheck = -1;

% Thanks Mathworks, for making another unwarranted change.
existisstupid = which('addParameter');
if (isempty(existisstupid))
   addParamValue(p, 'FlagLevel', defaultFlagLevel, @isnumeric);
   addParamValue(p, 'OldFile',   defaultOldFile,   @ischar);
   addParamValue(p, 'NewFile',   defaultNewFile,   @ischar);
   addParamValue(p, 'VarCheck',  defaultVarCheck,  @isnumeric);
else
   addParameter(p, 'FlagLevel', defaultFlagLevel, @isnumeric);
   addParameter(p, 'OldFile',   defaultOldFile,   @ischar);
   addParameter(p, 'NewFile',   defaultNewFile,   @ischar);
   addParameter(p, 'VarCheck',  defaultVarCheck,  @isnumeric);
end

p.parse(varargin{:}) % parse inputs

% collect the results of parsing (makes code easier to read)

FlagLevel = p.Results.FlagLevel;
OldFile   = p.Results.OldFile;
NewFile   = p.Results.NewFile;
VarCheck  = p.Results.VarCheck;

if (exist(OldFile,'file') ~= 2), error('File %s does not exist.',OldFile); end
if (exist(NewFile,'file') ~= 2), error('File %s does not exist.',NewFile); end

%% Create the list of vertical profile prior observation types in both files.

verticalobs = parse_DART_vars(OldFile, NewFile); 
nvariables = length(verticalobs);

%% plot some reference plot just to make sure we're not upside down or ...

if ( VarCheck > nvariables )
   fprintf('\nThere are only %d possible variables in %s\n',nvariables,OldFile)
   for ivar=1:nvariables
      fprintf('%40s is VarCheck %d\n',verticalobs{ivar},ivar)
   end
   error('VarCheck must be less than %d.',nvariables)

elseif ( VarCheck > 0 )
   close all
   files = {OldFile, NewFile};
   titles = {'old','new'};
   obsnames{1} = verticalobs{VarCheck};
   copy = 'bias';
   prpo = 'forecast';
   % fires up N figure windows ... for each region
   two_experiments_profile(files, titles, obsnames, copy, prpo)
   figure
end

%% plot the new stuff 

oldncid = netcdf.open(OldFile,'NOWRITE');
newncid = netcdf.open(NewFile,'NOWRITE');

dimid = netcdf.inqDimID(oldncid,'region');
[~, nregions] = netcdf.inqDim(oldncid,dimid);
varid = netcdf.inqVarID(oldncid,'region_names');
full_region_names = netcdf.getVar(oldncid,varid)';
region_name = cell(1,nregions);
for iregion=1:nregions
   region_name{iregion} = deblank(full_region_names(iregion,:));
end

pressure_dimid = netcdf.inqDimID(oldncid,'plevel'); 

plevels     = getvar(oldncid,'plevel');
hlevels     = getvar(oldncid,'hlevel');

biasindex = get_copy_index(OldFile,'bias');
rmseindex = get_copy_index(OldFile,'rmse');
nusedindx = get_copy_index(OldFile,'Nused');

f1 = gcf;    clf(f1); orient landscape
f2 = figure; clf(f2); orient landscape
f3 = figure; clf(f3); orient landscape
f4 = figure; clf(f4); orient landscape

for ivar = 1:nvariables
    
    varname = sprintf('%s_VPguess',verticalobs{ivar});
    olddata = getvar(oldncid, varname);
    newdata = getvar(newncid, varname);
    
    fprintf('Working on %s\n',varname)
    
    diffs   = olddata - newdata;
    rmsemat = squeeze(diffs(:,:,rmseindex))';
    nummat  = squeeze(diffs(:,:,nusedindx))';
    biasmat = squeeze(abs(olddata(:,:,biasindex)) - abs(newdata(:,:,biasindex)))';
    
    % ignore the levels where there are no 'old' observations.
    
    oldmask = (olddata(:,:,nusedindx) < 1)';
    rmsemat(oldmask) = 0.0;
    biasmat(oldmask) = 0.0;
    
    varid = netcdf.inqVarID(oldncid,varname);
    [~, ~, dimids, ~] = netcdf.inqVar(oldncid,varid);
    
    if (any(dimids == pressure_dimid))
        % fprintf('Pressure vertical coords for %s\n',varname)
        levels = plevels;
    else
        % fprintf('Height vertical coords for %s\n',varname)
        levels = hlevels;
    end
    
    % Plot the number of observations used and create the mask of which
    % cells are weakest.
    
    figure(f3);
    subplot(5,4,ivar)
    
    [hall, legendstr, weak] = obsplot(olddata(:,:,nusedindx)', newdata(:,:,nusedindx)', ...
        verticalobs{ivar}, levels, region_name, '% observations used', FlagLevel);
    
    if(ivar == nvariables)
        legend(hall,legendstr)
        xlabel('max observations per region')
    end
    
    % Plot the rmse
    
    figure(f1);
    subplot(5,4,ivar)
    myplot(rmsemat, verticalobs{ivar}, levels, region_name, 'rmse(old) - rmse(new)', weak)
    
    % Plot the bias
    
    figure(f2);
    subplot(5,4,ivar)
    myplot(biasmat, verticalobs{ivar}, levels, region_name, 'abs(bias(old)) - abs(bias(new))', weak)
    
    % Plot the difference of the number of observations used
    % also spew these to the command line for context ...
    
    figure(f4);
    subplot(5,4,ivar)
    
    plot(nummat, 1:length(levels), 'LineWidth',3.0);
    axis([-Inf Inf 0 length(levels)])
    set(gca,'YTick',1:length(levels))
    h = title({verticalobs{ivar},'obs used; old - new'});
    set(h,'interpreter','none');
    grid on
    
    if(ivar == nvariables)
        legend(region_name)
    end
    
    disp('old observation counts')
    olddata(:,:,nusedindx)'
    disp('new observation counts, then (old-new)')
    [newdata(:,:,nusedindx)' diffs(:,:,nusedindx)']
    
end

netcdf.close(oldncid)
netcdf.close(newncid)

bob = flipud(redblue(64));
figure(f1); colormap(bob)
figure(f2); colormap(bob)

end

%%----------------------------------------------------------------------

function vars = parse_DART_vars(file1, file2)

% rips off the _guess _analy _VPguess _VPanaly, finds the unique ones,
% and sorts them more-or-less alphabetically.

obstype    = get_DARTvars(file1);
replace    = '';
expression =   '_guess'; obstype = regexprep(obstype,expression,replace);
expression =   '_analy'; obstype = regexprep(obstype,expression,replace);
expression = '_VPguess'; obstype = regexprep(obstype,expression,replace);
expression = '_VPanaly'; obstype = regexprep(obstype,expression,replace);

% see if the _VPguess version of the variable exists in both files.

info1 = ncinfo(file1);
info2 = ncinfo(file2);

mask = false(1,length(obstype));

for ivar = 1:length(obstype)
   
   varname    = sprintf('%s_VPguess',obstype{ivar});
   
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

vars = sort(unique(obstype(mask)));

end

%%----------------------------------------------------------------------

function myplot(datmat, varname, levels, regions, titlestring, weak)

[ny, nx] = size(datmat);
xticks   = 1:nx;
yticks   = 1:ny;
imagesc(xticks,yticks,datmat)

% make a symmetric colormap ...
datmax = max(abs(datmat(:)));

% xticklabel = {'North','Tropics','South'};
xticklabel = cell(1,length(regions));
for i=1:length(regions)
    bob = mysplit(regions{i});
    xticklabel{i} = bob{1};
end

yticklabel = cell(1,length(levels));
for i=1:length(levels)
    yticklabel{i} = sprintf('%d',round(levels(i)));
end

set(gca,'YDir','normal', ...
    'XTick',xticks,'XTickLabel',xticklabel, ...
    'YTick',yticks,'YTickLabel',yticklabel, ...
    'CLim',[-datmax datmax]);
h = title({varname, titlestring});
set(h,'interpreter','none');
h = colorbar;
xh = get(h,'YLabel');
set(xh,'String','positive means new is better');

% plot a symbol at each of the weak region/levels

[i,j] = find(weak > 0);

line(j,i,'LineStyle','none','Marker','x','MarkerSize',10)

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

% create a mask for the 'weak, sparsely-observed' cells
weak = (old < FlagLevel) | (new < FlagLevel);

% then offset from one another
shift = ones(nlevels,1) * [0:nregions-1];
old   = old + shift;
new   = new + shift;

h1 = plot(old, 1:nlevels, '*-','LineWidth',1.0);
h2 = line(new, 1:nlevels, 'LineWidth',1.0, 'Marker','o');

axis([-Inf Inf 0 nlevels])
set(gca,'YTick',1:nlevels,'XTick',0:nregions,'XTickLabel',xticklabel)
h = title({varname,titlestring});
set(h,'interpreter','none');
grid on
hall = [h1; h2];

legendstr = cell(1,nregions*2);
for iregion = 1:nregions
    legendstr{         iregion} = sprintf('old %s',regions{iregion});
    legendstr{nregions+iregion} = sprintf('new %s',regions{iregion});
end

end

%%----------------------------------------------------------------------

function var = getvar(ncid, varname)
varid = netcdf.inqVarID(ncid,varname);
var   = netcdf.getVar(ncid,varid);
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

% subplot(2,1,1)
% rgbplot(bob)

% subplot(2,1,2)
% rgbplot(colors)

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

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$



