function figdata = set_obsdiag_figure(orientation,varargin)
%%
%  figure out a page layout
%  extra space at the bottom for the date/file annotation
%  extra space at the top because the titles have multiple lines
%
%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

default_nexp = 1; % The number of experiments
p = inputParser;

addRequired(p,'orientation',@ischar);

if (exist('inputParser/addParameter','file') == 2)
    addParameter(p,'numexp', default_nexp, @isnumeric);
else
    addParamValue(p,'numexp',default_nexp, @isnumeric); %#ok<NVREPL>
end

p.parse(orientation,varargin{:});

nexp = p.Results.numexp;

if ~isempty(fieldnames(p.Unmatched))
    disp('Extra inputs:')
    disp(p.Unmatched)
end

if strncmpi(orientation,'tall',4)
    orientation = 'tall';
    position = [0.15 0.12 0.7 0.75];
    
    if (nexp > 1)   % to replicate the 'two_experiments' behaviour
        ybot        = 0.06 + nexp*0.035;  % room for dates/files
        ytop        = 0.125;              % room for title (always 2 lines)
        dy          = 1.0 - ytop - ybot;
        position    = [0.15 ybot 0.7 dy];
    end
    
else
    orientation = 'landscape';
    position = [0.10 0.15 0.8 0.7];
    
    if (nexp > 1)   % to replicate the 'two_experiments' behaviour
        ybot        = 0.06 + nexp*0.075;  % room for dates/files
        ytop        = 0.125;              % room for title (always 2 lines)
        dy          = 1.0 - ytop - ybot;
        position    = [0.10 ybot 0.8 dy];
    end
    
end

fontsize      = 16;
linewidth     = 2.5;
obs_color     = [215/255  10/255  83/255]; % obs_red
ges_color     = [  0/255 128/255   0/255]; % prior_green
anl_color     = [  0/255   0/255 255/255]; % poste_blue
rmse_color    = [  0/255   0/255   0/255]; % black
copy_color    = [  0/255 128/255 128/255]; % teal
purple        = [ 153,51,255 ]/255;
orange        = [ 255,153,51 ]/255;
obs_marker    = 'o';
ges_marker    = '*';
anl_marker    = 'd';
marker1       = 'o';
marker2       = 's';
ges_linestyle = '-';
anl_linestyle = '-';
dashed        = '--';
solid         = '-';

figdata = struct( ...
    'expcolors',  {{'k','b','m','g','c','y','r'}}, ...
    'expsymbols', {{'o','s','d','p','h','s','*'}}, ...
    'prpolines',  {{'-','--'}}, ...
    'position'     , position, ...
    'fontsize'     , fontsize, ...
    'orientation'  , orientation, ...
    'linewidth'    , linewidth, ...
    'obs_color'    , obs_color, ...
    'ges_color'    , ges_color, ...
    'anl_color'    , anl_color, ...
    'rmse_color'   , rmse_color, ...
    'copy_color'   , copy_color, ...
    'obs_marker'   , obs_marker, ...
    'ges_marker'   , ges_marker, ...
    'anl_marker'   , anl_marker, ...
    'marker1'      , marker1, ...
    'marker2'      , marker2, ...
    'ges_linestyle', ges_linestyle, ...
    'anl_linestyle', anl_linestyle, ...
    'dashed'       , dashed, ...
    'solid'        , solid );

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
