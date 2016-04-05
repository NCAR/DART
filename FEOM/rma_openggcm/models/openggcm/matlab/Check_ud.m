function Check_ud(openggcmfile,dartfile)
%% Check_ud : check openggcm_to_dart.f90 ... the conversion of a openggcm restart to a DART state vector file.
%
% openggcmfile = '../data/da0002.dart.pot' 
% dartfile     = 'openggcm.pot.dat';
% Check_ud(openggcmfile, dartfile);
%

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Read the original openggcm file values.
if (exist(openggcmfile,'file') ~= 2)
   error('openggcm file %s does not exist.',openggcmfile)
end
if (exist(dartfile,'file') ~= 2)
   error('DART file %s does not exist.',dartfile)
end

ggcm   = fopen(openggcmfile,'r');
dart   = fopen(    dartfile,'r');

byte = fread(ggcm,1,'int32');
glon = fread(ggcm,1,'int32')
glat = fread(ggcm,1,'int32')
byte = fread(ggcm,1,'int32');

byte = fread(dart,1,'int32');
dlon = fread(dart,1,'int32')
dlat = fread(dart,1,'int32')
byte = fread(dart,1,'int32');

byte  = fread(ggcm,1,'int32');
gyear = fread(ggcm,1,'int32');
gmon  = fread(ggcm,1,'int32');
gday  = fread(ggcm,1,'int32');
ghour = fread(ggcm,1,'int32');
gmin  = fread(ggcm,1,'int32');
gsec  = fread(ggcm,1,'int32');
byte  = fread(ggcm,1,'int32');

byte  = fread(dart,1,'int32');
dyear = fread(dart,1,'int32');
dmon  = fread(dart,1,'int32');
dday  = fread(dart,1,'int32');
dhour = fread(dart,1,'int32');
dmin  = fread(dart,1,'int32');
dsec  = fread(dart,1,'int32');
byte  = fread(dart,1,'int32');

[ gyear gmon gday ghour gmin gsec ;
  dyear dmon dday dhour dmin dsec ]

byte = fread(ggcm,1,'int32');
gdat = fread(ggcm,glon*glat,'float32');
byte = fread(ggcm,1,'int32');

byte = fread(dart,1,'int32');
ddat = fread(dart,dlon*dlat,'float32');
byte = fread(dart,1,'int32');

difference = gdat - ddat;
plot(difference,'*')
max(difference)

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

