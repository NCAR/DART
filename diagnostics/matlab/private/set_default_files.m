function [truth_file, diagn_file] = set_default_files()
%% set_diagn_file sets the global variable 'diagn_file' that is used in
% the following routines:  plot_total_err, diagn_file

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

bob = evalin('base','exist(''truth_file'',''var'')');
if (bob == 1)
    truth_file = evalin('base','truth_file');
else
    truth_file = 'true_state.nc';
end

bob = evalin('base','exist(''diagn_file'',''var'')');
if (bob == 1)
    diagn_file = evalin('base','diagn_file');
else
    diagn_file = 'preassim.nc';
end

assignin('base','truth_file',truth_file)
assignin('base','diagn_file',diagn_file)

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
