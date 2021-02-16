function inflation_test
% compare update_inflate.m results to adaptive_inflate_mod.f90:update_inflation()
%
% The testing process is to :
% 1) enable adaptive_inflate_mod.f90 to write out 'inflation_input.m',
% 2) compile the L63 filter
% 3) replace the input.nml with input.inflationtest.nml
% 4) run filter
% 5) run this script
% 6) choose another inflation algorithm in input.nml and repeat 4,5
%
% This script is designed to be run in DART_LAB/matlab/private
% so that it can use the private function update_inflate.m directly,
% however, 'filter' writes out the file in the L63 'work' directory.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

addpath('../../../../models/lorenz_63/work')

if exist('inflation_input.m','file') ~= 2
   fprintf('Matlab path unable to locate "lorenz_63/work/inflation_input.m"\n')
   error('matlabpath needs the models/lorenz_63/work directory.')
end

clear all
inflation_input

% These are the only two options implemented in Matlab at the moment
if flavor == 2
   distribution = 'Gaussian';
else
   distribution = 'I-Gamma';
end

% There is no previous inflation happenning ...
ss_inflate_base = 1.0;

[matlab_inflate, matlab_inflate_sd] = update_inflate(prior_mean, ...
                                                   prior_var, ...
                                                   obs, ...
                                                   obs_var, ...
                                                   ss_inflate_base, ...
                                                   inflate, ...
                                                   inflate_sd, ...
                                                   inf_lower_bound, ...
                                                   inf_upper_bound, ...
                                                   gamma_corr, ...
                                                   sd_lower_bound, ...
                                                   ens_size, ...
                                                   distribution);

fprintf('fortran inflation value, sd is %f %f\n',fortran_inflate, fortran_inflate_sd)
fprintf('matlab  inflation value, sd is %f %f\n', matlab_inflate,  matlab_inflate_sd)

