function [cam_main] = load_CAM(varname,ens,year,path,n)
% function load_CAM: loads ensemble-years for CAM4/6 reanalysis
% Input:  'varname' is the CAM met name; 'ens' is the ensemble member
%         'year' is year of CAM reanalysis
%         'path' is path to the CAM reanalysis; 'n' indicates Solar (1hr inst),1hr,or 3hr
% Output: 'cam_main' is the CAM array indexed by ensemble member and time

switch n
    
    case 'SOLAR' % SOLAR instantaneous 1hr
        cam_main = ncread([path ens '/CAM6_NR1.cpl_' ens '.ha2x1hi.' year '.nc'],varname);
        
    case 'hour3' % average 3hr
        cam_main = ncread([path ens '/CAM6_NR1.cpl_' ens '.ha2x3h.' year '.nc'],varname);
        
    case 'hour1' % average 1hr
        cam_main = ncread([path ens '/CAM6_NR1.cpl_' ens '.ha2x1h.' year '.nc'],varname);
        
    case 'CAM4'  % average 6hr
        cam_main = ncread([path ens '/CAM4_NR1.cpl_' ens '.ha2x1dx6h.' year '.nc'],varname);
        
end
