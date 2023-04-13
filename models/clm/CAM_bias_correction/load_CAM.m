function [cam_main] = load_CAM(varname,ens,year,path,n)
% function load_CAM: loads all years of tower met forcing (PLUMBER2 format)
% Input:  'varname' is the CAM met name; 'ens' is the ensemble member    
%         'year' is year of CAM reanalysis
%         'path' is path to the CAM reanalysis; 'n' indicates Solar (1hr inst),1hr,or 3hr
% Output: 'cam_main' is the CAM array indexed by ensemble and time

switch n

  case SOLAR  % SOLAR instantaneous 1hr
  cam_main = ncread([path ens '/CAM6_NR1.cpl_' ens '.ha2x1hi.' year '.nc'],varname); 
	
  case hour3 % average 3hr
  cam_main = ncread([path ens '/CAM6_NR1.cpl_' ens '.ha2x3h.' year '.nc'],varname); 

  case hour1 % average 1hr
  cam_main = ncread([path ens '/CAM6_NR1.cpl_' ens '.ha2x1h.' year '.nc'],varname); 

end
