function bob = get_varnames(fname)
%% get_varnames returns NON-coordinate variables.
%
% the result is a cell array of strings ... must use {} notation to address elements.
%
% EXAMPLE:
% fname    = 'obs_seq.final.nc';
% varnames = get_varnames(fname);
% varnames{:}
% nvars = length(varnames);
% disp(sprintf('first atmospheric variable (of %d) is %s',nvars,varnames{1}))

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

fileinfo    = ncinfo(fname);
Nvarnames   = length(fileinfo.Variables);

inds = [];

for i = 1:Nvarnames
    
    varname    = fileinfo.Variables(i).Name;
    isatmosvar = 1;
    
    % Reject the obvious coordinate variables and some that are
    % specific to DART
    rank = length(fileinfo.Variables(i).Size);

    if (rank == 1)
        dimname = fileinfo.Variables(i).Dimensions(1).Name;
    else
        dimname = [];
    end
    
    if (strcmp( varname ,            dimname)), isatmosvar = 0; end
    if (strcmp( varname ,         'inputnml')), isatmosvar = 0; end
    if (strcmp( varname ,      'time_bounds')), isatmosvar = 0; end
    if (strcmp( varname ,     'region_names')), isatmosvar = 0; end
    if (strcmp( varname ,     'CopyMetaData')), isatmosvar = 0; end
    if (strcmp( varname ,   'MemberMetadata')), isatmosvar = 0; end
    if (strcmp( varname , 'ObservationTypes')), isatmosvar = 0; end
    
    % keep track of the 'good' variables
    if (isatmosvar > 0)
        inds = [inds i];
    end
end

if (isempty(inds))
   error('There are only coordinate variables in %s',fname)
end

% coerce just the names into a cell array
ngood = length(inds);
bob = cell(ngood,1);
for i = 1:ngood
   bob{i} = fileinfo.Variables(inds(i)).Name;
end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
