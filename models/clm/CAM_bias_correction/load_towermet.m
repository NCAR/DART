function [tower_main] = load_towermet(metvar,years,path)
% function load_towermet: loads all years of tower met forcing (PLUMBER2 format)
% Input:  'metvar' is PLUMBER2 met data name; 'years' is array of tower years to load
%         'path' is path to the PLUMBER2 met data
% Output: 'tower_main' is concatenated met data set
tower_main=[];

for i = 1:length(years);
    tower_dummy=squeeze(ncread([path 'Forcing_Tower_US-NR1_' years{i} '_ver.2022.11.29.nc'],metvar));
    tower_main=[tower_main;tower_dummy];
    clear tower_dummy
end

end
