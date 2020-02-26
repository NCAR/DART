% Script to generate: 
%   (1) standard reanalysis diagnostics (.ps, .txt, .pdf files)
%   (2) an html directory useful to view everything in one place 

% Uses script "results.m" and function "invoke_diag.m"
% It creates a Diags_MMM-yyyy directory in the current working directory.

% Requirement: 'obs_diag_output.nc' file

clear
close all
clc

global path

view_format             = 'html';
dont_show_code          = false ;
fig_format              = 'png' ;

path.scripts_dir        = pwd;

% Figure which month are we doing:
obs_file = 'obs_diag_output.nc';
if exist(obs_file, 'file') ~= 2
    error('DART diagnostics file "obs_diag_output.nc" does not exit in this directory.')
end

Time                    = ncread(obs_file, 'time');
unit_time               = ncreadatt(obs_file, 'time', 'units');

origin                  = regexp(unit_time, '\d*', 'Match');
origin                  = datenum(str2double(origin));
current                 = Time + origin;
period                  = datestr(current, 'mmm-yyyy');

reana_month             = period(ceil(length(period)/2), :);
path.obs_space_dir      = strcat(path.scripts_dir, '/', 'Diags_', reana_month, '/');
path.web_dir            = strcat(path.obs_space_dir, 'web_', reana_month);  

if ~ exist(path.obs_space_dir, 'dir')
    mkdir(path.obs_space_dir)
end

system ([ 'cp obs_diag_output.nc ' path.obs_space_dir ]);

path.obs_space_diags    = strcat(path.obs_space_dir, 'obs_diag_output.nc');
path.inflation          = '';
path.fireoff_script     = strcat(path.scripts_dir, '/results.m');


% Fill in the publishing options
options = struct(   'format'        , view_format       , ...
                    'showCode'      , dont_show_code    , ...
                    'imageFormat'   , fig_format        , ...
                    'outputDir'     , path.web_dir      );
                
publish(path.fireoff_script, options);

system ([ 'mv *.ps *.pdf *.txt ' path.obs_space_dir ]);
system ([ 'tar -zcvf ' path.web_dir '.tar.gz ' path.web_dir ]);
