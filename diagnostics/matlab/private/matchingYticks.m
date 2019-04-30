function matchingYticks(ax1, ax2)
%% This takes the existing Y ticks from ax1 (presumed nice)
% and determines nice-looking matching labels for ax2.
%
% In our implementation (for _evolution.m applications):
% the data plotted on ax1 are the quantities of interest.
% the data plotted on ax2 are the observation counts.
%
% USAGE - given two coincident axes (see plot_evolution.m for an example):
%
% matchingYticks(ax1,ax2);

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

Dlimits = get(ax1,'YLim');
DYticks = get(ax1,'YTick');
nYticks = length(DYticks);
ylimits = get(ax2,'YLim');

% make sure the #_obs ticks have whole numbers.
% The log of the y range is useful because it separates the magnitude from the number.
% The ylimits sometimes exaggerate actual range of obs values too much.

yrange_log = log10(ylimits(2) - ylimits(1));

% Separate the log of the range (magnitude.span) into the 2 parts.
% E.g. yrange = 1500 means that yrange_log = 3.176,
%   so the size is 10^3 and the span is 10^.176 = 1.5

yrange_magnitude = floor(yrange_log);
yrange_span      = 10^(yrange_log - yrange_magnitude);

% Match right axis ticks to left axis ticks
% Generate a nicely sized spacing between ticks, rounded up from the original minimum size.
% In matlab 10^[log10(x)] = x + delta where delta can be > eps (2.e-16)

delta = 1.0e-14;

% Rescale the magnitude and span if the span is to small to be useful.

if ((nYticks -1  - yrange_span) > delta )
    yrange_span      = yrange_span * 10;
    yrange_magnitude = yrange_magnitude - 1;
end

% Here's the distance between obs ticks, in units of # of observations.
% It's used to label the ticks and determine the new YLim.
% Prevent fractional tick spans (when there are only a few obs, and relatively many ticks).

yrange_tick_size = ceil((yrange_span - delta) / (nYticks -1)) * (10^yrange_magnitude);
yrange_tick_size = max([yrange_tick_size 1.0]);

% Set the tick values for the obs number axis

yticks0 = floor(ylimits(1)/yrange_tick_size -1) * yrange_tick_size;
iticks  = 1:nYticks;
yticks  = yticks0 + yrange_tick_size*iticks;

% New axis bottom and top for obs number axis.
% These are returned to the calling program.

D_tick_size = DYticks(2) - DYticks(1);
ylimits(1)  = yticks(   1   ) - yrange_tick_size*((DYticks(1) - Dlimits(   1   ))/D_tick_size);
ylimits(2)  = yticks(nYticks) + yrange_tick_size*((Dlimits(2) - DYticks(nYticks))/D_tick_size);

newticklabels = num2str(round(10*yticks')/10);

% use the new ticks and labels.

set(ax2,'YTick', yticks, 'YTickLabel', newticklabels,'YLim', ylimits)

% sometimes the ax2 position changes from the new ticks. Make sure ax1 follows.

set(ax1,'Position',get(ax2,'Position'))
grid

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
