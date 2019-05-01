function bob = get_DARTvars(fname)
%% get_DARTvars returns just the variable names in the netCDF file that are likely to
% be DART variables.
%
% the result is a cell array of strings ... must use {} notation to address elements.
%
% EXAMPLE:
% fname    = 'preassim.nc';
% DARTvars = get_DARTvars(fname);
% DARTvars{:}
% nvars = length(DARTvars);
% disp(sprintf('first atmospheric variable (of %d) is %s',nvars,DARTvars{1}))

% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

fileinfo  = ncinfo(fname);
nvars     = length(fileinfo.Variables);
isDARTvar = zeros(nvars,1);

for i = 1:nvars

    rank = length(fileinfo.Variables(i).Dimensions);
    for idim = 1:rank
        
       % Anything with a 'member' dimension is probably a DART state vector variable.
       dimname = fileinfo.Variables(i).Dimensions(idim).Name;
       if (strcmp(dimname,'member')), isDARTvar(i) = 1; end

       % Reject the obvious metadata variable
       varname = fileinfo.Variables(i).Name;
       if (strcmp(varname, 'MemberMetadata')), isDARTvar(i) = 0; end

       % Anything with a 'copy' dimension is probably a DART state vector variable.
       dimname = fileinfo.Variables(i).Dimensions(idim).Name;
       if (strcmp(dimname,'copy')), isDARTvar(i) = 1; end

       % Reject the obvious metadata variable
       varname = fileinfo.Variables(i).Name;
       if (strcmp(varname, 'CopyMetaData')), isDARTvar(i) = 0; end

    end

    % If the variable is 1D and the name and the dimension name are the same,
    % the variable is a coordinate variable and can be rejected.

    if (rank == 1 && strcmp(varname,fileinfo.Variables(i).Dimensions(1).Name)),
       isDARTvar(i) = 0;
    end

end

if (sum(isDARTvar) == 0)
   error('No DART state variables in %s',fname)
end

% coerce just the names into a cell array
bob = cell(sum(isDARTvar),1);
nDARTvars = 0;
for i = 1:nvars
   if (isDARTvar(i) > 0)
      nDARTvars = nDARTvars + 1;
      bob{nDARTvars} = fileinfo.Variables(i).Name;
   end
end

% Each of the candidate variables may have a _mean or _sd counterpart
% TJH ... not sure why this is here ...

% varind = nDARTvars;
% extensions = {'mean','sd','priorinf_mean', 'priorinf_sd','postinf_mean','postinf_sd'};
% for i = 1:nDARTvars
%    for ext = 1:length(extensions)
%       % construct candidate name and see if it exists
%       suffix = extensions{ext};
%       candidate_name = sprintf('%s_%s',bob{i},suffix);
%       for ivar = 1:nvars
%          if ( strcmp(candidate_name, fileinfo.Variables(ivar).Name ) )
%             varind = varind + 1;
%             bob{varind} = candidate_name;
%          end
%       end
%    end
% end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
