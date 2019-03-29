function copy_index = get_copy_index(fname, copystring, varargin)
%% GET_COPY_INDEX  Gets an index corresponding to copy metadata string
% Retrieves index associated with a given string in the
% CopyMetaData netCDF variable in the given file. If the string
% does not exist - a fatal error is thrown.
%
% Example:
% fname = 'obs_diag_output.nc';
% copystring = 'N_DARTqc_5';
% copy_index = get_copy_index(fname, copystring);
%
% If you prefer the error to be non-fatal, you can do that too.
% If the index does not exist, the copy_index will be -1
% Example:
% copy_index = get_copy_index(fname, copystring, 'fatal', false);

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

defaultContext = [];
defaultFatality = true;
defaultVerbose = false;

p = inputParser;
addRequired(p,'fname',@ischar);
addRequired(p,'copystring',@ischar);

if (exist('inputParser/addParameter','file') == 2)
    addParameter(p,'context',defaultContext, @ischar);
    addParameter(p,'fatal',  defaultFatality,@islogical);
    addParameter(p,'verbose',defaultVerbose, @islogical);
else
    addParamValue(p,'context',defaultContext, @ischar);
    addParamValue(p,'fatal',  defaultFatality,@islogical);
    addParamValue(p,'verbose',defaultVerbose, @islogical);
end

p.parse(fname, copystring, varargin{:});

errorstring = sprintf('\nERROR: "%s" is not a valid CopyMetaData value for file %s\n', ...
    strtrim(p.Results.copystring), p.Results.fname);

if (isempty(p.Results.context))
    msgstring = 'valid values for CopyMetaData are';
else
    msgstring = sprintf('valid CopyMetaData values for "%s" are', p.Results.context);
end

if ( exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

% Matlab seems to always need to transpose character variables.
copy_meta_data  = ncread(fname,'CopyMetaData')';
[num_copies, ~] = nc_dim_info(fname,'copy');
[metalen, ~]    = nc_dim_info(fname,'stringlength');

if( size(copy_meta_data,1) ~= num_copies || size(copy_meta_data,2) ~= metalen)
    error('%s from %s does not have the shape expected',copystring,fname)
end

nowhitecs = dewhite(copystring);

% Figure out which copy is the matching one
copy_index = -1;
for i = 1:num_copies,

    % for matching -- we want to ignore whitespace -- find it & remove it
    nowhitemd = dewhite(copy_meta_data(i,:));

    if strcmp(nowhitemd , nowhitecs) == 1
        copy_index = i;
    end
end

% Provide modest error support

if (copy_index < 0 && p.Results.fatal)
    for i = 1:num_copies,
        msgstring = sprintf('%s\n%s',msgstring,deblank(copy_meta_data(i,:)));
    end
    error('%s\n%s',errorstring,msgstring)
end


function str2 = dewhite(str1)
%  function to remove ALL whitespace from a character string

str2 = str1(~isspace(str1));


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
