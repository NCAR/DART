function h = cable_hop_test
%% DART: Compares two CABLE restart files for bit-for-bit behavior.
%
% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

file1 = 'CABLE_restart.onehop.nc';
fileN = 'CABLE_restart.Nhops.nc';

[~, ~] = system('rm Compare_netCDF_files_results.txt');
% Compare_netCDF_files(fname1,fnameN)

gridfile = 'CABLE_gridinfo.nc';

latitude  = nc_varget(file1,'latitude');
longitude = nc_varget(file1,'longitude');

mp_patch = size(latitude);

%% Some of these are 2D variables, some are 1D
%  all the 2D variables have mp_patch as the second matlab dimension

varlist = { 'tgg', 'wb', 'wbice', 'tss', 'albsoilsn', 'rtsoil', 'gammzz', 'runoff', ...
   'rnof1', 'rnof2', 'tggsn', 'ssdnn', 'ssdn', 'snowd', 'smass', 'sdepth', ...
   'osnowd', 'isflag', 'cansto', 'ghflux', 'sghflux', 'ga', 'dgdtg', 'fev', ...
   'fes', 'fhs', 'wbtot0', 'osnowd0', 'albedo', 'trad', 'ejmax', 'vcmax'};

varlist = { 'tgg' };

for ivar =1:length(varlist)
   
   varunits = nc_attget(file1,varlist{ivar},'units');
   slab1    = nc_varget(file1,varlist{ivar});
   slabN    = nc_varget(fileN,varlist{ivar});
   
   if (numel(slab1) > mp_patch)
      fprintf('Plotting the 2D difference of %s ...\n',varlist{ivar})
      h = compare_2D_variables(slab1, slabN, latitude, longitude, varlist{ivar},varunits);
   else
      fprintf('Plotting the 1D difference of %s ...\n',varlist{ivar})
      h = compare_1D_variables(slab1, slabN, latitude, longitude, varlist{ivar},varunits);
   end
   
end


function h = compare_2D_variables(datmat1, datmat2, lats, lons, varname, units)

[nlevels, nx] = size(datmat1);

for ilevel = 1:nlevels
   x = datmat2-datmat1;
   
   h = scatter(lons,lats,10,x(ilevel,:));
   set(h,'LineWidth',1.5)
   h1 = title(sprintf('%s level %d of %d',varname,ilevel,nlevels));
   set(h1,'Interpreter','none')
   h2 = colorbar;
   set(get(h2,'YLabel'),'String',sprintf('difference - %s',units))
   disp('pausing, hit any key to continue ...'); pause
end

function h = compare_1D_variables(datmat1, datmat2, lats, lons, varname, units)

   x = datmat2-datmat1;
   
   h = scatter(lons,lats,10,x);
   set(h,'LineWidth',1.5)
   h1 = title(varname);
   set(h1,'Interpreter','none')
   h2 = colorbar;
   set(get(h2,'YLabel'),'String',sprintf('difference - %s',units))
   disp('pausing, hit any key to continue ...'); pause


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

