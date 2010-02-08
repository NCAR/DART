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

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
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

   qr    = nc{'QRAIN'}(:,:,:) ;
   qg    = nc{'QGRAUP'}(:,:,:) ;
   qs    = nc{'QSNOW'}(:,:,:) ;
   close(nc);
 elseif length( varargin ) == 2
   time_index = varargin{1} ; mem_index = varargin{2} ; 
      % note that DART diagnostic files include multiple times, members
 % Read all data
   nc = netcdf( filename , 'nowrite' ) ;

   qr    = squeeze(nc{'QRAIN'}(time_index,mem_index,:,:,:)) ;
   qg    = squeeze(nc{'QGRAUP'}(time_index,mem_index,:,:,:)) ;
   qs    = squeeze(nc{'QSNOW'}(time_index,mem_index,:,:,:)) ;
   close(nc);
 elseif length( varargin ) == 3
   time_index = varargin{1} ; mem_index = varargin{2} ; id = varargin{3} ;
% note that DART diagnostic files include multiple times, members, domains

  qr = getnc(filename,['QRAIN_d0',int2str(id)],[time_index mem_index -1 -1 -1], ...
	     [time_index mem_index -1 -1 -1],[1 1 1 1 1]);
  qg = getnc(filename,['QGRAUP_d0',int2str(id)],[time_index mem_index -1 -1 -1], ...
	     [time_index mem_index -1 -1 -1],[1 1 1 1 1]);
  qs = getnc(filename,['QSNOW_d0',int2str(id)],[time_index mem_index -1 -1 -1], ...
	     [time_index mem_index -1 -1 -1],[1 1 1 1 1]);
 else
   disp('*** Incorrect number of input arguments in get_aux_fields_for_p')
 end
