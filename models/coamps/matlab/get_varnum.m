function varnum=get_varnum(restart_info_file,varname,level)
%% varnum = get_varnum(restart_info_file, varname, level)
%
% Given a variable name and sigma index as well as the location of
% the restart.vars file used for the ensemble run, generates the
% variable number in the long state vector.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

  % format string for the restart.vars file
  % u2 01 NOPERTS 1.0 M QTY_U_WIND UPDATE FALSE NOPOSDEF
  %fmt_str='%s %d %s %f %c %s %s %s %s';
  fmt_str = '%s %d %*[^\n]';  
  [name,sigma_level] = textread(restart_info_file,fmt_str, ...
				'headerlines',1);

  % Convert the supplied name to a cell array so we can use the
  % string comparison function (and don't need to worry about
  % string padding)
  name_cell = cellstr(varname);
  name_match = find(strcmp(name,name_cell));
  
  % Now we need to match the levels - the equality wants the sizes
  % to match
  level_match= find(sigma_level == repmat(level,[length(sigma_level) 1]));

  % Find out where these match
  num_name_matches = length(name_match);
  num_lev_matches = length(level_match);
  Name_indices = repmat(name_match,[1 num_lev_matches]);
  Lev_indices  = repmat(level_match.',[num_name_matches 1]);
  IsMatch = Name_indices - Lev_indices;
  [r,c] = find(IsMatch == 0);
  
  % Grab the index
  varnum = Name_indices(r,c);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
