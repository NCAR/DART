function PlotRegFactor(fname,obsnum)
%% PlotRegFactor : Plots regression factor
%
% at present, the regression factor files are a bit non-conforming.
% They don't have enough metadata to make for nice labelling.
%
% Example 1:
% fname  = 'time_mean_reg';
% obsnum = 20;
% PlotRegFactor(fname, obsnum);

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if ( exist(fname,'file') ~=2 ), error('%s does not exist on the Matlab search path',fname); end

% Open the file and read the data into blobs.

[fid,message] = fopen(fname,'rt');
if (fid <= 0 ), disp(message); end
A = fscanf(fid,'%d',2);
B = fscanf(fid,'%f');
fclose(fid);

% recover the dimensions.

num_obs_in_set = A(1);
model_size     = A(2);
c = reshape(B,[3 num_obs_in_set*model_size]);

obs_num = c(1,:);
mod_num = c(2,:);
reg_fac = c(3,:);
obs_mat = reshape(obs_num,[model_size num_obs_in_set]);
mod_mat = reshape(mod_num,[model_size num_obs_in_set]);
reg_mat = reshape(reg_fac,[model_size num_obs_in_set]);

clear obs_num mod_num reg_fac

% [ obs_mat(:, obsnum) reg_mat(:,obsnum) ]

plot(reg_mat(:,obsnum))
ax = axis;
axis([ax(1) ax(2) 0 1])
h = title({fname}); set(h,'interpreter','none','fontsize',16)
ylabel('regression factor')
xlabel({'state variable (indexical)', ...
        sprintf('observation number %d',obsnum)})

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
