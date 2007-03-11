function var_vec = get_var_series(fname, varname, copynum, state_var)
%GET_VAR_SERIES Gets a particular copy of a state variable from netcdf file
%
% Retrieves a particular copy of a state variable from a file whose
% full or relative path is specified in the file argument.
%
% Example 1:
% fname     = '../work/Prior_Diag.nc';
% varname   = 'state';      % State Variable
% copynum   = 8;            % Ensemble Member
% state_var = 3;            % which state variable
% var_vec   = get_var_series(fname, varname, copynum, state_var);

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2007, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

f = netcdf(fname,'nowrite');
var_atts   = dim(f{varname});       % cell array of dimensions for the var
num_copies = length(var_atts{2});
num_vars   = length(var_atts{3});

if ( ~ strcmp( name(var_atts{1}), 'time') )
    disp( sprintf('%s first dimension ( %s ) is not ''time''',fname,name(var_atts{1})))
end
if ( ~ strcmp( name(var_atts{2}), 'copy') )
    disp( sprintf('%s second dimension ( %s ) is not ''copy''',fname,name(var_atts{2})))
end
if (copynum > num_copies ) 
    disp( sprintf('%s only has %d ''copies/Ensemble members of %s''',fname,num_copies,varname))
    error(sprintf('you wanted copy %d ', copynum))
end
if (state_var > num_vars) 
   disp( sprintf('%s only has %d %s variables',fname,num_vars,varname))
   error(sprintf('you wanted variable %d ', state_var))
end
close(f);

% Get only the appropriate copy of the state and return
var_vec = getnc(fname, varname, [-1, copynum, state_var], ...
                                [-1, copynum, state_var]);

