function [dart modl] = Check_ncommas_to_dart(modlfile, dartfile )
%% Check_ncommas_to_dart : check ncommas_to_dart.f90 ... the conversion of a NCOMMAS restart to a DART state vector file.
%
% modlfile = 'ncommas_restart.nc';
% dartfile = 'True_State.nc';
%
% [dart modl] = Check_ncommas_to_dart(modlfile, dartfile );

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Read the original NCOMMAS file values.
if (exist(modlfile,'file') ~= 2)
   error('NCOMMAS file %s does not exist.',modlfile)
end
if (exist(dartfile,'file') ~= 2)
   error('DART file %s does not exist.',dartfile)
end

iyear   = nc_attget(modlfile,nc_global,'YEAR');
imonth  = nc_attget(modlfile,nc_global,'MONTH');
iday    = nc_attget(modlfile,nc_global,'DAY');
ihour   = nc_attget(modlfile,nc_global,'HOUR');
iminute = nc_attget(modlfile,nc_global,'MINUTE');
isecond = nc_attget(modlfile,nc_global,'SECOND');

fprintf('NCOMMAS year  month  day  hour  minute  second %d %d %d %d %d %d\n',  ...
        iyear,imonth,iday,ihour,iminute,isecond);

% The nc_varget() function returns the variables with the fastest 
% varying dimension on the right. This is opposite to the Fortran
% convention of the fastest varying dimension on the left ... so 
% one of the variables must be permuted in order to be compared.

U  = nc_varget(modlfile,  'U'); modl.U  = permute( U, [4 3 2 1]);
V  = nc_varget(modlfile,  'V'); modl.V  = permute( V, [4 3 2 1]);
W  = nc_varget(modlfile,  'W'); modl.W  = permute( W, [4 3 2 1]);
QS = nc_varget(modlfile, 'QS'); modl.QS = permute(QS, [4 3 2 1]);

disp(sprintf('modl.W min/max are %0.8g %0.8g',min(modl.W(:)),max(modl.W(:))))

[nx ny nz nt] = size(modl.U);
fprintf('vert dimension size is %d\n',nz)
fprintf('N-S  dimension size is %d\n',ny)
fprintf('E-W  dimension size is %d\n',nx)

% And now for the DART-generated netCDF files - that have the copy dimension
% nc_varget squeezes out singleton dimensions ...

dU  = nc_varget(dartfile,  'U'); dart.U  = permute( dU, [3 2 1]);
dV  = nc_varget(dartfile,  'V'); dart.V  = permute( dV, [3 2 1]);
dW  = nc_varget(dartfile,  'W'); dart.W  = permute( dW, [3 2 1]);
dQS = nc_varget(dartfile, 'QS'); dart.QS = permute(dQS, [3 2 1]);

 Udiffs =  modl.U(:,:,:,2) -  dart.U; [min( Udiffs(:)) max( Udiffs(:))]
 Vdiffs =  modl.V(:,:,:,2) -  dart.V; [min( Vdiffs(:)) max( Vdiffs(:))]
 Wdiffs =  modl.W(:,:,:,2) -  dart.W; [min( Wdiffs(:)) max( Wdiffs(:))]
QSdiffs = modl.QS(:,:,:,2) - dart.QS; [min(QSdiffs(:)) max(QSdiffs(:))]

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
