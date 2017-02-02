function [ qr, qg, qs ] = get_aux_fields_for_ref( filename, varargin )
%% [ qr, qg, qs ] = get_aux_fields_for_ref( filename, varargin )
%
% Retrieves various 3d fields needed to calculate reflectivity
% from the netcdf file "filename".
%
% Other inputs: varargin = {}                             (if filename is a wrfinput file)
%                          {time_index, member_index}     (if filename is a DART diagnostic file)
%                          {time_index, member_index, id} (if there is more than 1 domain)
% Outputs: qr = rain water mixing ratio (3d)
%          qg =    graupel mixing ratio (3d)
%          qs =       snow mixing ratio (3d)
%
% Example to read a DART netcdf file.
% filename = 'Prior_Diag.nc';
% time_index = 1;
% copy_index = 1;
% dom_id = 1;
% [ qr, qg, qs ] = get_aux_fields_for_ref( filename, time_index, copy_index, dom_id );

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
   varexist(filename, {'QRAIN','QGRAUP','QSNOW'});

   qr = nc_varget(filename, 'QRAIN');
   qg = nc_varget(filename,'QGRAUP');
   qs = nc_varget(filename, 'QSNOW');

elseif length( varargin ) == 2

   %% Read all data from DART diagnostic file - no domain
   %  note that DART diagnostic files include multiple times, copies

   time_index = varargin{1};
   copy_index = varargin{2}; 

   varexist(filename, {'QRAIN','QGRAUP','QSNOW'});

   start = [time_index, copy_index,  1,  1,  1] - 1;
   count = [         1,          1, -1, -1, -1];
   qr    = nc_varget(filename, 'QRAIN',start,count);
   qg    = nc_varget(filename,'QGRAUP',start,count);
   qs    = nc_varget(filename, 'QSNOW',start,count);

elseif length( varargin ) == 3

   %% note that DART diagnostic files include multiple times, copies, domains

   time_index = varargin{1};
   copy_index = varargin{2};
   id         = varargin{3};

   start = [time_index, copy_index,  1,  1,  1] - 1;
   count = [         1,          1, -1, -1, -1];

   varexist(filename, {[ 'QRAIN_d0',int2str(id)], ...
                       ['QGRAUP_d0',int2str(id)], ...
                       [ 'QSNOW_d0',int2str(id)]} );

   qr = nc_varget(filename,[ 'QRAIN_d0',int2str(id)], start, count);
   qg = nc_varget(filename,['QGRAUP_d0',int2str(id)], start, count);
   qs = nc_varget(filename,[ 'QSNOW_d0',int2str(id)], start, count);

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
