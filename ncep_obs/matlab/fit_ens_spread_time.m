function fit_ens_spread_time(ddir)
% fit_ens_spread_time(ddir)
%
%
% USAGE
%
% fit_ens_spread_time('plot')

% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$

% This ensures the directory with the datafiles
% is in Matlab's search path.

if ( nargin > 0 )
   if ( exist(ddir,'dir') ~= 7)
      error(sprintf('%s is not a valid directory',ddir))
   end
else
   ddir = [];  % defaults to local directory.
end

p     = load(fullfile(ddir,'Tanl_times_level.dat'));
level = p(1);

Tges  = fullfile(ddir,sprintf('Tges_times_%04dmb.dat',level));
Tanl  = fullfile(ddir,sprintf('Tanl_times_%04dmb.dat',level));
Wges  = fullfile(ddir,sprintf('Wges_times_%04dmb.dat',level));
Wanl  = fullfile(ddir,sprintf('Wanl_times_%04dmb.dat',level));

%----------------------------------------------------------------------
figure(1); clf; % Temperature
%----------------------------------------------------------------------

   str1   = sprintf('%4d hPa T  Ens Spread',level);
   region = 'Northern Hemisphere';
   topstr = {str1, region};
   Myplot(Tges,Tanl,topstr,1)
   
   region = 'Southern Hemisphere';
   topstr = {str1, region};
   Myplot(Tges,Tanl,topstr,2)
   
   region = 'Tropics';
   topstr = {str1, region};
   Myplot(Tges,Tanl,topstr,3)
   
   region = 'North America';
   topstr = {str1, region};
   Myplot(Tges,Tanl,topstr,4)
   
%----------------------------------------------------------------------
figure(2); clf; % Now for the Winds
%----------------------------------------------------------------------
   
   str1   = sprintf('%4d hPa Wind Ens Spread',level); 
   region = 'Northern Hemisphere';
   topstr = {str1, region};
   Myplot(Wges,Wanl,topstr,1)
   
   region = 'Southern Hemisphere';
   topstr = {str1, region};
   Myplot(Wges,Wanl,topstr,2)
   
   region = 'Tropics';
   topstr = {str1, region};
   Myplot(Wges,Wanl,topstr,3)
   
   region = 'North America';
   topstr = {str1, region};
   Myplot(Wges,Wanl,topstr,4)

%----------------------------------------------------------------------
% common
%----------------------------------------------------------------------

   print(1,'-dpsc','t_ens_spread_time.ps');
   print(2,'-dpsc','w_ens_spread_time.ps');



function Myplot(file1,file2,topstr,region)
%
% After the first column, each set of 3 columns
% represents a different region.
% Depends on the format written by new_obs_diag.f90
%
p         = load(file1);
a         = load(file2);
xp        = p(:,1);
xa        = a(:,1);
count     = 3+(region-1)*3;
yp_spread = p(:,count);
ya_spread = a(:,count);
%
subplot(2,2,region)
plot(xp, yp_spread, 'k+-', xa, ya_spread, 'ro-', 'LineWidth', 1.5)
grid
xlabel('Time interval', 'fontsize', 10) ;
ylabel('RMSE ', 'fontsize', 10)
title(topstr, 'fontsize', 14, 'FontWeight', 'bold')
legend('guess', 'analysis')
h = legend;
legend(h,'boxoff')
