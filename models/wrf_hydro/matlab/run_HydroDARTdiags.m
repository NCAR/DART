clear
close all
clc

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
% 
% Example script to run the HydroDARTdiags.m function
%
% DART $Id: run_HydroDARTdiags.m $

%% Diagnostic parameters

% Which experiment(s) is/are requested: 
just_1_exp = true;

dir_exps = {'florence/manu/R_mp80_prior_post_inf_loc200'};
% dir_exps = {'florence/manu/mp80_prior_post_inf_l200', 'florence/manu/R_mp80_prior_post_inf_loc200' };

if length(dir_exps) > 1
    just_1_exp = false;
end


% Open-loop directory if available;
% if not dir_ol should be removed from HydroDARTdiags call
dir_ol = {'florence/manu/mp80_ol'}; 


% Select the gauges where hydrographs are displayed
% options: 1, 2, 3, 4, 5, 6
% 1: See obs1 set below
% 2: See obs2 set below
% 3: See obs3 set below
% 4: Gauges displayed in the HESS paper (El Gharamti et al., 2021)
% 5: Gauges that only exist in the florence cut-out domain
% 6: All available gauges

which_gauges = 5; 

obs1 = [ 208735012, 2093877  , 2102908  , 2101800  , 2094659  , 2082950  , ...
         2094775  , 2096846  , 208524090, 208675010, 210166029, 2077303  , ...
         2095181  , 209722970, 2097464  , 208732885, 2097280  , 2091000  , ...
         209782609, 209399200, 2077200  , 2099000  , 2095271  , 2130840  , ...
         209734440, 2093800  , 2084160  , 2087580  , 2094770  , 208524975, ...
         2086849  , 2095500  , 2097517  , 208758850, 208521324, 209553650  ];
     
obs2 = [ 209741955, 2128000  , 2074500, 2092500  , 2087275  , 2133500, ...
         2081500  , 2086500  , 2130900, 2085070  , 2096960  , 2131000, ...
         208773375, 2088000  , 2087324, 2094500  , 2075500  , 2105769, ...
         2090380  , 2073000  , 2130910, 2134480  , 2129000  , 2102500, ...
         2097314  , 2074000  , 2096500, 2082770  , 212378405, 2102000, ...
         2087183  , 2100500  , 2071000, 2130980  , 2089500  , 2089000  ];
     
obs3 = [ 2132320, 2104220, 208726005, 2105500  , 208250410, 2081942, ...
         2081747, 2109500, 2088500  , 2103000  , 2083000  , 2088383, ...
         2085000, 2087359, 2086624  , 208111310, 2085500  , 2101726, ...
         2082585, 2075045, 2083500  , 2134500  , 2098206  , 2106500, ...
         2102192, 2095000, 2093000  , 2077670  , 2130561  , 2126000, ...
         2108000, 2134170, 2091500  , 2087500  , 2133624  ] ;
     
local_gauges = [02075500, 02082585, 02083500, 02102000]; 
inflt_gauges = [02087500, 02089000, 0208250410];
verfy_gauges = [0209553650, 02134170, 02105769];
biased_gauge = [02126000, 02130561];

obs4 = [local_gauges, inflt_gauges, verfy_gauges, biased_gauge]; 
obs5 = [02086849, 0208521324, 02085000, 0208524975, 02085500, 02086500, 02085070];
obs6 = 0;

switch which_gauges
    case 1, obs = obs1;
    case 2, obs = obs2;
    case 3, obs = obs3;
    case 4, obs = obs4;
    case 5, obs = obs5;
    case 6, obs = obs6;
end


%% Calling HydroDARTdiags

% Toggle to display the results
disp_res = true;

% State-space stream networks
if just_1_exp
    plot_state = 1; % optional, slows things down (0 to deactivate)  
else
    plot_state = 0; % always inactive here    
end

% y : observations statistics
% Xo: open loop statistics
% Xf: Prior statistics
% Xa: Posterior statistics
% E : Ensemble/inflation data 
[y, Xo, Xf, Xa, E] = HydroDARTdiags(dir_exps, obs, dir_ol, disp_res, plot_state);


% % <next few lines under version control, do not edit>
% % $URL: $
% % $Revision: $
% % $Date: $
