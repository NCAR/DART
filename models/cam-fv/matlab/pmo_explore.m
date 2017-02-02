function inds = pmo_explore(experiment,baseday,compday)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% baseday = 5;
% compday = 6;
% cam4_d5 = pmo_explore('CAM4',baseday,compday)
%
% baseday = 1;
% compday = 2;
% cam4_d1 = pmo_explore('CAM4',baseday,compday)

inds          = [];
region        = [0 360 -90 90 -Inf Inf];
QCString      = 'Quality Control';
verbose       = 0;
ObsTypeString = {'RADIOSONDE_U_WIND_COMPONENT',
                 'RADIOSONDE_V_WIND_COMPONENT',
                 'RADIOSONDE_TEMPERATURE'};

figure(1); orient tall

for itype = 1:length(ObsTypeString)

   file1  = sprintf('obs_epoch_noon_%s_%03d.nc',experiment,baseday);
   truth1 = read_obs_netcdf(file1, ObsTypeString{itype}, region, 'truth',        QCString, verbose);
   obs1   = read_obs_netcdf(file1, ObsTypeString{itype}, region, 'observations', QCString, verbose);
   obserrors1 = obs1.obs-truth1.obs;
 
   file2  = sprintf('obs_epoch_noon_%s_%03d.nc',experiment,compday);
   truth2 = read_obs_netcdf(file2, ObsTypeString{itype}, region, 'truth',        QCString, verbose);
   obs2   = read_obs_netcdf(file2, ObsTypeString{itype}, region, 'observations', QCString, verbose);
   obserrors2 = obs2.obs-truth2.obs;

   diffs = obserrors2 - obserrors1;
   absx  = abs(diffs);
   nsync = sum(isfinite(absx) & (absx <= 0.001));

   subplot(length(ObsTypeString),1,itype)
   plot(diffs)
   h = title({ObsTypeString{itype},sprintf('day%d - day%d %s',compday,baseday,experiment)});
   set(h,'Interpreter','none')
   xlabel('observation index')
   ylabel('observation noise difference')
   h = legend(sprintf('%d +/- 0.001',nsync));
   set(h,'Box','off')

end 

print(1,'-dpdf',sprintf('Zagar_PMO_%s_day%d%d.pdf',experiment,baseday,compday))

%%
%
%

if ( 1 == 2 ) 

   inds = zeros(92,3);
   
   for itype = 1:1
   
      file1 = sprintf('obs_epoch_noon_%s_%03d.nc',experiment,baseday);
      truth = read_obs_netcdf(file1, ObsTypeString{itype}, region, 'truth',        QCString, verbose);
      obs   = read_obs_netcdf(file1, ObsTypeString{itype}, region, 'observations', QCString, verbose);
      diffs1 = obs.obs-truth.obs;
   
      for iepoch = 1:92,
   
         fname = sprintf('obs_epoch_noon_%s_%03d.nc',experiment,iepoch);
         truth = read_obs_netcdf(fname, ObsTypeString{itype}, region, 'truth',        QCString, verbose);
         obs   = read_obs_netcdf(fname, ObsTypeString{itype}, region, 'observations', QCString, verbose);
         diffs2 = obs.obs-truth.obs;
   
         inds(iepoch,itype) = firstfinite(diffs2-diffs1);
   
         fprintf('%s day %d first diff obs noise is %d \n',...
                    ObsTypeString{itype},iepoch,inds(iepoch,itype))
   
      end
   
   end
   
   figure(2)
   xax = 1:size(inds,1);
   plot(xax,inds(:,1),'*',xax,inds(:,2),'d',xax,inds(:,3),'o')
   h = title({'First observation with substantially different obs noise', ...
          sprintf('between noontime obs on day %d and other days',baseday),...
          sprintf('experiment %s',experiment)});
   set(h,'Interpreter','none')
   h = legend(ObsTypeString{1},ObsTypeString{1},ObsTypeString{3});
   set(h,'Interpreter','none')
   xlabel('day')
   ylabel('obs index of first real difference')
   orient landscape
   grid
   
   print(2,'-dpdf',sprintf('Zagar_PMO_%s_day%d_obindex.pdf',experiment,baseday))

end

function [z,nsync] = firstfinite(x)

absx  =  abs(x);
nsync =  sum(isfinite(absx) & (absx <= 0.01));
y     = find(isfinite(absx) & (absx  > 0.01));

if (isempty(y))
   z = NaN;
else
   z = y(1);
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
