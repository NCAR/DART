function [ mu, dnw, phi, theta, qv ] =  ...
                      get_aux_fields_for_p( filename, T0, varargin )
%
% Retrieves various 2d and 3d fields needed to calculate pressure
% from the netcdf file "filename".
%
% Other inputs: T0       = wrf caries theta as deviation from this value
%               varargin = {} (if fields come from wrfinput files)
%                          {time_index, member_index}
%                             (if fields come from DART diagnostic files)
% Outputs: mu    = full mu (2d)
%          dnw   = intervals between w levels (1d)
%          phi   = full geopotential (3d)
%          theta = full theta (3d)
%          qv    = water-vapor mixing ratio (3d)

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

 % Retrieve required fields
 if isempty( varargin )
 % Read all data
   nc = netcdf( filename , 'nowrite' ) ;

   mu    = nc{'MU'}(:,:) + nc{'MUB'} ;
   dnw   = nc{'DNW'}(:) ;
   phi   = nc{'PH'}(:,:,:) + nc{'PHB'} ;
   theta = nc{'T'}(:,:,:) + T0 ;
   qv    = nc{'QVAPOR'}(:,:,:) ;
   close(nc);
 elseif length( varargin ) == 2
   time_index = varargin{1} ; mem_index = varargin{2} ; 
      % note that DART diagnostic files include multiple times, members
 % Read all data
   nc = netcdf( filename , 'nowrite' ) ;

   mu    = squeeze(nc{'MU'}(time_index,mem_index,  :,:)) + nc{'MUB'} ;
   dnw   = nc{'DNW'}(:) ;
   phi   = squeeze(nc{'PH'}(time_index,mem_index,:,:,:)) + nc{'PHB'} ;
   theta = squeeze(nc{'T'}(time_index,mem_index,:,:,:)) + T0 ;
   qv    = squeeze(nc{'QVAPOR'}(time_index,mem_index,:,:,:)) ;
   close(nc);
 elseif length( varargin ) == 3
   time_index = varargin{1} ; mem_index = varargin{2} ; id = varargin{3} ;
% note that DART diagnostic files include multiple times, members, domains
  mu = getnc(filename,['MU_d0',int2str(id)],[time_index mem_index -1 -1], ...
	     [time_index mem_index -1 -1],[1 1 1 1]) + ...
       getnc(filename,['MUB_d0',int2str(id)],[-1 -1],[-1 -1],[1 1]);
  dnw = getnc(filename,['DNW_d0',int2str(id)],[-1],[-1],[-1]);
  phi = getnc(filename,['PH_d0',int2str(id)],[time_index mem_index -1 -1 -1], ...
	      [time_index mem_index -1 -1 -1],[1 1 1 1 1]) + ...
        getnc(filename,['PHB_d0',int2str(id)],[-1 -1 -1],[-1 -1 -1],[1 1 1]);
  theta = getnc(filename,['T_d0',int2str(id)],[time_index mem_index -1 -1 -1], ...
	      [time_index mem_index -1 -1 -1],[1 1 1 1 1]) + T0;
  qv = getnc(filename,['QVAPOR_d0',int2str(id)],[time_index mem_index -1 -1 -1], ...
	     [time_index mem_index -1 -1 -1],[1 1 1 1 1]);
 else
   disp('*** Incorrect number of input arguments in get_aux_fields_for_p')
 end
