function obs_num_time(ddir)
% obs_num_time     Plots the total number of observations as a function of time for a given level for several different regions.
% 
% ddir     an optional argument specifying the directory containing
%          the data files as preprocessed by the support routines.
%
% USAGE: if the preprocessed data files are in a directory called 'plot'
% 
% obs_num_time('plot')

% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$

% Make sure the directory with the files is available to matlab
if (nargin > 0 )
   datafile = fullfile(ddir,'Tanl_times_level.dat');
else
   datafile = 'Tanl_times_level.dat';
   ddir = [];
end

if ( exist(datafile,'file') == 2 ) 
   p     = load(datafile);
   level = p(1);
else
   error('%s cannot be found.',datafile)
end

Tfname = fullfile(ddir,sprintf('Tges_times_%04dmb.dat',level));
Wfname = fullfile(ddir,sprintf('Wges_times_%04dmb.dat',level));

if ( exist(Tfname,'file') ~= 2 )
   error(sprintf('%s does not seem to exist.',Tfname))
end
if ( exist(Wfname,'file') ~= 2 )
   error(sprintf('%s does not seem to exist.',Wfname))
end

   %------------------------------------------------------------
   % Plot Number of Temperature Observations
   %------------------------------------------------------------

   figure(1); clf;

   ylab   = '# of Temp observations';

   region      = 'Northern Hemisphere';
   top         = sprintf('%s  %4d hPa',region,level);
   [nT_NH,xax] = Myplot(Tfname,top,ylab,1,'gs-');
   
   region      = 'Southern Hemisphere';
   top         = sprintf('%s  %4d hPa',region,level);
   [nT_SH,xax] = Myplot(Tfname,top,ylab,2,'bd-');
   
   region      = 'Tropics';
   top         = sprintf('%s  %4d hPa',region,level);
   [nT_TR,xax] = Myplot(Tfname,top,ylab,3,'ro-');
   
   region      = 'North America';
   top         = sprintf('%s  %4d hPa',region,level);
   [nT_NA,xax] = Myplot(Tfname,top,ylab,4,'k+-');

   print -dpsc         t_obs_num_time.ps
  
   % All regions on one figure

   figure(2); clf

   h = plot(xax, nT_NH, 'gs-', ...
            xax, nT_SH, 'bd-', ...
	    xax, nT_TR, 'ro-', ...
	    xax, nT_NA, 'k+-', 'LineWidth', 2.0);
   grid
   %set(gca,'YDir', 'reverse')
   title(sprintf('# of Temperature Observations at %d hPa',level), ...
                'FontSize', 14, 'FontWeight', 'bold')
   ylabel('# of Observations', 'fontsize', 10)
   xlabel('Time Interval','fontsize',10)
   h = legend('Northern Hemisphere', ...
              'Southern Hemisphere', ...
	      'Tropics', ...
	      'North America','Location','NorthWest');
   legend(h,'boxoff');

   print -dpsc -append t_obs_num_time.ps

   %------------------------------------------------------------
   % Now for the Winds
   %------------------------------------------------------------

   figure(3); clf;
   
   ylab   = '# of Wind observations';

   region      = 'Northern Hemisphere';
   top         = sprintf('%s  %4d hPa',region,level);
   [nW_NH,xax] = Myplot(Wfname,top,ylab,1,'gs-');
   
   region      = 'Southern Hemisphere';
   top         = sprintf('%s  %4d hPa',region,level);
   [nW_SH,xax] = Myplot(Wfname,top,ylab,2,'bd-');
   
   region      = 'Tropics';
   top         = sprintf('%s  %4d hPa',region,level);
   [nW_TR,xax] = Myplot(Wfname,top,ylab,3,'ro-');
   
   region      = 'North America';
   top         = sprintf('%s  %4d hPa',region,level);
   [nW_NA,xax] = Myplot(Wfname,top,ylab,4,'k+-');

   print -dpsc         w_obs_num_time.ps

   % All regions on one figure

   figure(4); clf;

   h = plot(xax, nW_NH, 'gs-', ...
            xax, nW_SH, 'bd-', ...
	    xax, nW_TR, 'ro-', ...
	    xax, nW_NA, 'k+-', 'LineWidth', 2.0);
   grid
   % set(gca,'YDir', 'reverse')
   title(sprintf('# of Wind Observations at %d hPa',level), ...
                'FontSize', 14, 'FontWeight', 'bold')
   ylabel('# of Observations', 'fontsize', 10)
   xlabel('Time Interval','fontsize',10)
   h = legend('Northern Hemisphere', ...
              'Southern Hemisphere', ...
	      'Tropics', ...
	      'North America','Location','NorthWest');
   legend(h,'boxoff');

   print -dpsc -append w_obs_num_time.ps

function [yp_num,xp] = Myplot(file1,top,ylab,region,ptype)

p      = load(file1);
xp     = p(:,1);
count  = 4+(region-1)*3;
yp_num = p(:,count);

subplot(2,2,region)
   plot(xp,yp_num,ptype,'LineWidth',2.0)
   grid
   xlabel('Time Interval', 'fontsize', 10) ;
   ylabel(ylab, 'fontsize', 10) ;
   title(top, 'fontsize', 14,'FontWeight','bold')
