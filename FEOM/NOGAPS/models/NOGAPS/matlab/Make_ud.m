function dart = Make_ud(dartfile)
%% Make_ud : creates a synthetic DART state vector that has bogus, 
% but different, values for each variable/level. This file can be 
% used to check the dart_to_nogaps.F90 to ensure the parsing/unparsing
% is happening as expected.  
%
% Make_ud does read the 'allones.ics' file created by check_model_mod.f90
% simply to get the fortran record information. 
%
% EXAMPLE: 
% dartfile = 'dart_ics';
% x        = Make_ud( dartfile );

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$
%
% This is from the original assembla server we used during collaboration.
% $orgURL: https://svn2.assembla.com/svn/ngdart/matlab/Make_ud.m $
% $orgId: Make_ud.m 96 2010-05-26 17:06:04Z thoar $
% $orgRevision: 96 $
% $orgDate: 2010-05-26 11:06:04 -0600 (Wed, 26 May 2010) $

nx      = 144;
ny      = 72;
nz      = 42;
dart.seconds = 21600;
dart.days    = 149446; % 3 March 2010 

modelsize = nz*ny*nx;

%% create some fake data

datmat    = peaks(nx);

dart.U  = zeros(nx,ny,nz);
dart.V  = zeros(nx,ny,nz);
dart.T  = zeros(nx,ny,nz);
dart.Q  = zeros(nx,ny,nz);
dart.PS = ones(nx,ny) * 1013.0;

for i = 1:nz
   bob = (i/nz)*datmat;
   dart.U(:,:,i) = 0 + (nz-i+1) + bob(:,1:ny);
   dart.V(:,:,i) = 1 + (nz-i+1) + bob(:,1:ny);
   dart.T(:,:,i) = 2 + (nz-i+1) + bob(:,1:ny);
   dart.Q(:,:,i) = 3 + (nz-i+1) + bob(:,1:ny);
end

% Read a dummy binary state vector to ensure we have
% right record header lengths, etc.

% Open and read timetag for state
fid     = fopen('allones.ics','rb','ieee-le');
rec1    = fread(fid, 1, 'int32') 
seconds = fread(fid, 1, 'int32');
days    = fread(fid, 1, 'int32');
rec2    = fread(fid, 1, 'int32') 
rec3    = fread(fid, 1, 'int32') 

% filesize = U,V,T,Q * (nx*ny*nz) + PS * (nx*ny)
storage  = 8;
size2d   = nx*ny;
size3d   = nx*ny*nz;
n2ditems = 1*size2d;
n3ditems = 4*size3d;
rec1size = 4+(4+4)+4;  % time stamps ... 
rec2size = 4+(n3ditems*storage + n2ditems*storage)+4;

fsize    = rec1size + rec2size;
disp(sprintf('with a modelsize of %d the file size should be %d bytes', ...
     modelsize,fsize))

% Open and read timetag for state
fid = fopen(dartfile,'wb','ieee-le');
fwrite(fid,         rec1, 'int32');
fwrite(fid, dart.seconds, 'int32');
fwrite(fid,    dart.days, 'int32');
fwrite(fid,         rec2, 'int32');

fprintf('DART days seconds %d %d\n', dart.days, dart.seconds);

% Successively write state vector variables.
fwrite( fid, rec3, 'int32');
fwrite( fid, dart.U  , 'float64');
fwrite( fid, dart.V  , 'float64');
fwrite( fid, dart.T  , 'float64');
fwrite( fid, dart.Q  , 'float64');
fwrite( fid, dart.PS , 'float64');
fwrite( fid, rec3, 'int32');
fclose(fid);

% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

