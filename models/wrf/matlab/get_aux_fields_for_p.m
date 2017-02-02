function [ mu, dnw, phi, theta, qv ] =  ...
                      get_aux_fields_for_p( filename, T0, varargin )
%% get_aux_fields_for_p - Retrieves various 2d and 3d fields needed to calculate pressure
% from the netcdf file "filename".
%
% Other inputs: T0       = wrf caries theta as deviation from this value
%               varargin = {} (if fields come from wrfinput files)
%                          {time_index, copy_index}
%                             (if fields come from DART diagnostic files)
% Outputs: mu    = full mu (2d)
%          dnw   = intervals between w levels (1d)
%          phi   = full geopotential (3d)
%          theta = full theta (3d)
%          qv    = water-vapor mixing ratio (3d)
%
% Example to read a DART netcdf file.
% filename = 'Prior_Diag.nc';
% time_index = 1;
% copy_index = 1;
% dom_id = 1;
% T0 = 270.0;   % this is an example - I have no idea what it should be.
% [ mu, dnw, phi, theta, qv ] =  ...
%       get_aux_fields_for_p( filename, T0, time_index, copy_index, dom_id  )

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (exist(filename,'file') ~= 2)
   error('%s does not exist.',filename)
end

%% Retrieve required fields

if isempty( varargin )

   %% Read all data from 'wrfinput' style files

   varexist(filename,{'MU','MUB','DNW','PH','PHB','T','QVAPOR'})

   dnw   = nc_varget(filename,   'DNW');
   qv    = nc_varget(filename,'QVAPOR');
   theta = nc_varget(filename,     'T') + T0;
   mu    = nc_varget(filename,    'MU') + nc_varget(filename,'MUB');
   phi   = nc_varget(filename,    'PH') + nc_varget(filename,'PHB');

elseif length( varargin ) == 2

   %% Read all data from DART diagnostic file - no domain
   %  note that DART diagnostic files include multiple times, copies

   varexist(filename,{'MU','MUB','DNW','PH','PHB','T','QVAPOR'})

   time_index = varargin{1};
   copy_index = varargin{2}; 

   start = [time_index, copy_index,  1,  1] - 1;
   count = [         1,          1, -1, -1];

   dnw   = nc_varget(filename,   'DNW');
   qv    = nc_varget(filename,'QVAPOR',start,count);
   theta = nc_varget(filename,     'T',start,count) + T0;
   mu    = nc_varget(filename,    'MU',start,count) + nc_varget(filename,'MUB');
   phi   = nc_varget(filename,    'PH',start,count) + nc_varget(filename,'PHB');

elseif length( varargin ) == 3

   %% note that DART diagnostic files include multiple times, copies, domains

   time_index = varargin{1};
   copy_index = varargin{2};
   id         = varargin{3};

   varexist(filename,{[    'MU_d0',int2str(id)], ['MUB_d0',int2str(id)], ...
                      [    'PH_d0',int2str(id)], ['PHB_d0',int2str(id)], ...
                      [     'T_d0',int2str(id)], ...
                      [   'DNW_d0',int2str(id)], ...
                      ['QVAPOR_d0',int2str(id)]} )

   start = [time_index, copy_index,  1,  1] - 1;
   count = [         1,          1, -1, -1];

   dnw = nc_varget(filename,['DNW_d0',int2str(id)]);
   mu  = nc_varget(filename,[ 'MU_d0',int2str(id)], start, count) + ...
         nc_varget(filename,['MUB_d0',int2str(id)]);

   start = [time_index, copy_index,  1,  1,  1] - 1;
   count = [         1,          1, -1, -1, -1];

   qv    = nc_varget(filename,['QVAPOR_d0',int2str(id)], start, count);
   theta = nc_varget(filename,[     'T_d0',int2str(id)], start, count) + T0;
   phi   = nc_varget(filename,[    'PH_d0',int2str(id)], start, count) + ...
           nc_varget(filename,[   'PHB_d0',int2str(id)]);

else
   disp('*** Incorrect number of input arguments in get_aux_fields_for_p')
end


function varexist(filename, varnames)
%% We already know the file exists by this point.
% Lets check to make sure that file contains all needed variables.

nvars = length(varnames);

for i = 1:nvars

   if ( ~ nc_isvar(filename,varnames{i}) )
      error('%s is not a variable in %s',varnames{i},filename)
   end

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
