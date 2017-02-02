function xscale = matchingXticks(ax1, ax2)
%% This takes the existing X ticks from ax1 (presumed nice)
% and determines the nice, matching labels for ax2.
%
% In our implementation (for _profile.m applications):
% the data plotted on ax1 are the quantities of interest.
% the data plotted on ax2 are the observation counts.
%
% xscale = matchingXticks(ax1, ax2)
%
% xscale   the scaling factor between ticks and labels ...
%
% USAGE - given two coincident axes (see plot_profile.m for an example):
%
% xscale = matchingXticks(ax1,ax2);
% set(get(ax2,'Xlabel'),'String',['# of obs (o=poss, +=used) x' int2str(uint32(xscale))])

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

Dlimits = get(ax1,'XLim');
DXticks = get(ax1,'XTick');
nXticks = length(DXticks);
xlimits = get(ax2,'XLim');

% make sure the #_obs ticks have whole numbers.
% The log of the x range is useful because it separates the magnitude from the number.
% The xlimits sometimes exaggerate actual range of obs values too much.

xrange_log = log10(xlimits(2) - xlimits(1));

% Separate the log of the range (magnitude.span) into the 2 parts.
% E.g. xrange = 1500 means that xrange_log = 3.176,
%   so the size is 10^3 and the span is 10^.176 = 1.5

xrange_magnitude = floor(xrange_log);
xrange_span      = 10^(xrange_log - xrange_magnitude);

% Match top axis ticks to bottom axis ticks
% Generate a nicely sized spacing between ticks, rounded up from the original minimum size.
% In matlab 10^[log10(x)] = x + delta where delta can be > eps (2.e-16)

delta = 1.0e-14;

% Rescale the magnitude and span if the span is to small to be useful.

if ((nXticks -1  - xrange_span) > delta )
   xrange_span      = xrange_span * 10 ;
   xrange_magnitude = xrange_magnitude - 1;
end

% pass the scaling factor back for the axis label

xscale = 10.^xrange_magnitude;

% Here's the distance between obs ticks, in units of SCALED # of observations.
% It's used to label the ticks and determine the new XLim.
% Prevent fractional tick spans (when there are only a few obs, and relatively many ticks).

xrange_tick_size = ceil((xrange_span - delta) / (nXticks - 1)) * xscale ;
xrange_tick_size = max([xrange_tick_size 1.0]);

% Set the tick values for the obs number axis

xticks0 = floor(xlimits(1)/xrange_tick_size - 1) * xrange_tick_size;
iticks  = 1:nXticks;
xticks  = xticks0 + xrange_tick_size*iticks;

% New axis left and right for obs number axis.
% These are returned to the calling program.

D_tick_size = DXticks(2) - DXticks(1);
xlimits(1)  = xticks(   1   ) - xrange_tick_size*((DXticks(1) - Dlimits(   1   ))/D_tick_size);
xlimits(2)  = xticks(nXticks) + xrange_tick_size*((Dlimits(2) - DXticks(nXticks))/D_tick_size);

newticklabels = num2str(round((10/xscale)*xticks') /10);

% use the new ticks and labels.

set(ax2,'XTick', xticks, 'XTicklabel', newticklabels, 'XLim', xlimits)

% sometimes the ax2 position changes from the new ticks. Make sure ax1 follows.

set(ax1,'Position',get(ax2,'Position'))
grid

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
