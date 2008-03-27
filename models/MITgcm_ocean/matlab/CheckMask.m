function CheckMask()
% CheckMask 
% 
%

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2007, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------

mitbase = '/fs/image/home/nancy/subversion/trunk/models/MITgcm_ocean/data2/';

S     = rdmds(sprintf('%s/%s.0000040992',mitbase,'S'));
Sinds = find(S == 0.0); clear S;

T     = rdmds(sprintf('%s/%s.0000040992',mitbase,'T'));
Tinds = find(T == 0.0); clear T;

U     = rdmds(sprintf('%s/%s.0000040992',mitbase,'U'));
Uinds = find(U == 0.0); clear U;

V     = rdmds(sprintf('%s/%s.0000040992',mitbase,'V'));
Vinds = find(V == 0.0); clear V;

SSH = rdmds(sprintf('%s/%s.0000040992',mitbase,'Eta'));

disp(sprintf('S has %d zeros',length(Sinds)))
disp(sprintf('T has %d zeros',length(Tinds)))
disp(sprintf('U has %d zeros',length(Uinds)))
disp(sprintf('V has %d zeros',length(Vinds)))

