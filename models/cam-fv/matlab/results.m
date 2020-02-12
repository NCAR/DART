%% *REANALYSIS OBS-SPACE DIAGNOSTICS*
% This is an html file containing diagnostics summary for the reanalysis
% project. Below are profile and time-series evolution plots for 21
% observations kinds (see Contents). The results are split into 3 regions:
%
% # Southern Hemisphere: |LAT = [-90, -20]| 
% # Tropics: |LAT = [-20, 20]|
% # Northern Hemisphere: |LAT = [20, 90]|
%
% For the vertical profiles, the results are calculated at:
%
% * Pressure layers (hPa): |999, 925, 831, 687, 525, 400, 312, 250, 200, 150,
% 100, 55, 25, 8| 
% * Height layers (m): |315, 750, 1405, 2775, 4675, 6560, 8285, 9830, 11410,
% 13470, 16435, 20890, 26525, 36245| 
%
% For each observation kind, the root-mean-squared-error (RMSE) is plotted
% together with: (i) _BIAS_ and (ii) _TOTALSPREAD_. The number of
% available, assimilated (evaluated) and rejected observations are also
% shown on the plots. 
%
% For more information, see our web site:
% <https://dart.ucar.edu/ *DART-GROUP*>

%%
% 
% This script is called by "gen_rean_diags.m" and it's not intended to be run 
% on its own. 

global path
close

% observations:
diag.types = {  'ACARS_TEMPERATURE'           , ...
                'ACARS_U_WIND_COMPONENT'      , ...
                'ACARS_V_WIND_COMPONENT'      , ...
                'ACARS_HORIZONTAL_WIND'       , ... 
                'AIRCRAFT_TEMPERATURE'        , ...
                'AIRCRAFT_U_WIND_COMPONENT'   , ...
                'AIRCRAFT_V_WIND_COMPONENT'   , ...
                'AIRCRAFT_HORIZONTAL_WIND'    , ...
                'AIRS_TEMPERATURE'            , ...
                'AIRS_SPECIFIC_HUMIDITY'      , ...
                'GPSRO_REFRACTIVITY'          , ...
                'LAND_SFC_ALTIMETER'          , ...
                'MARINE_SFC_ALTIMETER'        , ...
                'RADIOSONDE_TEMPERATURE'      , ...
                'RADIOSONDE_U_WIND_COMPONENT' , ...
                'RADIOSONDE_V_WIND_COMPONENT' , ...
                'RADIOSONDE_SURFACE_ALTIMETER', ...
                'RADIOSONDE_SPECIFIC_HUMIDITY', ...
                'SAT_U_WIND_COMPONENT'        , ...
                'SAT_V_WIND_COMPONENT'        , ...
                'SAT_HORIZONTAL_WIND'         , ...
                                                };

% styles:
diag.style = {  'profile'                     , ...
                'evolution'                     };
            
% metrics: 
diag.metric = { 'bias'                        , ...
                'totalspread'                   };


% for now, i have to list all of them manually (and not in a loop) because
% i need a different title for each section. MATLAB only supports static
% titles for sections (cells). One could explore using MATLAB REPORT
% GENERATOR to do this more efficiently. 
            
%% 1.a. *ACARS_TEMPERATURE:* _bias_
var = diag.types{1};
met = diag.metric{1};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 1.b. *ACARS_TEMPERATURE:* _totalspread_
var = diag.types{1};
met = diag.metric{2};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 2.a. *ACARS_U_WIND_COMPONENT:* _bias_
var = diag.types{2};
met = diag.metric{1};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 2.b. *ACARS_U_WIND_COMPONENT:* _totalspread_
var = diag.types{2};
met = diag.metric{2};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 3.a. *ACARS_V_WIND_COMPONENT:* _bias_
var = diag.types{3};
met = diag.metric{1};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 3.b. *ACARS_V_WIND_COMPONENT:* _totalspread_
var = diag.types{3};
met = diag.metric{2};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 4.a. *ACARS_HORIZONTAL_WIND:* _bias_
var = diag.types{4};
met = diag.metric{1};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 4.b. *ACARS_HORIZONTAL_WIND:* _totalspread_
var = diag.types{4};
met = diag.metric{2};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 5.a. *AIRCRAFT_TEMPERATURE:* _bias_
var = diag.types{5};
met = diag.metric{1};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 5.b. *AIRCRAFT_TEMPERATURE:* _totalspread_
var = diag.types{5};
met = diag.metric{2};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 6.a. *AIRCRAFT_U_WIND_COMPONENT:* _bias_
var = diag.types{6};
met = diag.metric{1};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 6.b. *AIRCRAFT_U_WIND_COMPONENT:* _totalspread_
var = diag.types{6};
met = diag.metric{2};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 7.a. *AIRCRAFT_V_WIND_COMPONENT:* _bias_
var = diag.types{7};
met = diag.metric{1};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 7.b. *AIRCRAFT_V_WIND_COMPONENT:* _totalspread_
var = diag.types{7};
met = diag.metric{2};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 8.a. *AIRCRAFT_HORIZONTAL_WIND:* _bias_
var = diag.types{8};
met = diag.metric{1};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 8.b. *AIRCRAFT_HORIZONTAL_WIND:* _totalspread_
var = diag.types{8};
met = diag.metric{2};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 9.a. *AIRS_TEMPERATURE:* _bias_
var = diag.types{9};
met = diag.metric{1};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 9.b. *AIRS_TEMPERATURE:* _totalspread_
var = diag.types{9};
met = diag.metric{2};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 10.a. *AIRS_SPECIFIC_HUMIDITY:* _bias_
var = diag.types{10};
met = diag.metric{1};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 10.b. *AIRS_SPECIFIC_HUMIDITY:* _totalspread_
var = diag.types{10};
met = diag.metric{2};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 11.a. *GPSRO_REFRACTIVITY:* _bias_
var = diag.types{11};
met = diag.metric{1};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 11.b. *GPSRO_REFRACTIVITY:* _totalspread_
var = diag.types{11};
met = diag.metric{2};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 12.a. *LAND_SFC_ALTIMETER:* _bias_
var = diag.types{12};
met = diag.metric{1};
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 12.b. *LAND_SFC_ALTIMETER:* _totalspread_
var = diag.types{12};
met = diag.metric{2};
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 13.a. *MARINE_SFC_ALTIMETER:* _bias_
var = diag.types{13};
met = diag.metric{1};
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 13.b. *MARINE_SFC_ALTIMETER:* _totalspread_
var = diag.types{13};
met = diag.metric{2};
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 14.a. *RADIOSONDE_TEMPERATURE:* _bias_
var = diag.types{14};
met = diag.metric{1};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 14.b. *RADIOSONDE_TEMPERATURE:* _totalspread_
var = diag.types{14};
met = diag.metric{2};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 15.a. *RADIOSONDE_U_WIND_COMPONENT:* _bias_
var = diag.types{15};
met = diag.metric{1};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 15.b. *RADIOSONDE_U_WIND_COMPONENT:* _totalspread_
var = diag.types{15};
met = diag.metric{2};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 16.a. *RADIOSONDE_V_WIND_COMPONENT:* _bias_
var = diag.types{16};
met = diag.metric{1};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 16.b. *RADIOSONDE_V_WIND_COMPONENT:* _totalspread_
var = diag.types{16};
met = diag.metric{2};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 17.a. *RADIOSONDE_SURFACE_ALTIMETER:* _bias_
var = diag.types{17};
met = diag.metric{1};
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 17.b. *RADIOSONDE_SURFACE_ALTIMETER:* _totalspread_
var = diag.types{17};
met = diag.metric{2};
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 18.a. *RADIOSONDE_SPECIFIC_HUMIDITY:* _bias_
var = diag.types{18};
met = diag.metric{1};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 18.b. *RADIOSONDE_SPECIFIC_HUMIDITY:* _totalspread_
var = diag.types{18};
met = diag.metric{2};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 19.a. *SAT_U_WIND_COMPONENT:* _bias_
var = diag.types{19};
met = diag.metric{1};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 19.b. *SAT_U_WIND_COMPONENT:* _totalspread_
var = diag.types{19};
met = diag.metric{2};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 20.a. *SAT_V_WIND_COMPONENT:* _bias_
var = diag.types{20};
met = diag.metric{1};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 20.b. *SAT_V_WIND_COMPONENT:* _totalspread_
var = diag.types{20};
met = diag.metric{2};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 21.a. *SAT_HORIZONTAL_WIND:* _bias_
var = diag.types{21};
met = diag.metric{1};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)

%% 21.b. *SAT_HORIZONTAL_WIND:* _totalspread_
var = diag.types{21};
met = diag.metric{2};
invoke_diag(path.obs_space_diags, diag.style{1}, met, var)
invoke_diag(path.obs_space_diags, diag.style{2}, met, var, -1)
