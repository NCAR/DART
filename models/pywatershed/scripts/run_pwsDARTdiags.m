clear
close all
clc


% Where experiments live
source_dir = '/PATH/TO/PWS-DART/EXPERIMENTS/DIRECTORY/'; 
exps_dir   = strcat(source_dir, 'pywatershed/');
disp       = true;
dir_ol     = {strcat(exps_dir, 'OL')};
dir_exps   = {strcat(exps_dir, 'DA')};

% Where to save the figures
figs_dir = strcat(exps_dir, 'hydrographs/');

if ~isfolder(figs_dir), mkdir(figs_dir); end

% image format (only for separate gauges and NOT regions)
form = '.png'; %png, jpg, tif, emf, eps

% Run pwsDART diags
[y, Xo, Xf, Xa, E] = pwsDARTdiags(dir_exps, dir_ol, disp, figs_dir, form);
