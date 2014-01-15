function x = clm_get_var(fname,varname,copystring,levelindex,timeindex,clmfname)
%% DART clm_get_var - reads a variable from a CLM restart or DART diagnostic netCDF file and reconstitute a matrix.
%
% EXAMPLE 1:
% fname      = 'Prior_Diag.2000-01-05-00000.nc';
% varname    = 'frac_sno';
% levelindex = 1;
% timeindex  = 1;
% copystring = 'ensemble member 3';
% x          = clm_get_var(fname,varname,copystring,levelindex,timeindex);
%
% EXAMPLE 2: as above, compare to comparable field in CLM history file.
% clmfname   = '/glade/scratch/thoar/enstest_0907/enstest_0907.clm2_0003.r.2000-01-05-00000.nc';
% x          = clm_get_var(fname,varname,copystring,levelindex,timeindex,clmfname);
%
% fname      = 'Prior_Diag.2000-01-06-00000.nc';
% clmfname   = '../clmruns/enstest_0906.clm2_0006.r.2000-01-06-00000.nc';
%
% % 'frac_sno',    'KIND_SNOWCOVER_FRAC',
% % 'DZSNO',       'KIND_SNOW_THICKNESS',
% % 'H2OSOI_LIQ',  'KIND_LIQUID_WATER',
% % 'H2OSOI_ICE',  'KIND_ICE',
% % 'T_SOISNO',    'KIND_SOIL_TEMPERATURE',
%

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (exist(fname,'file') ~= 2)
   error('%s does not exist.',fname)
end
x.filename = fname;
x.clmfile  = [];

if (nargin == 6) 
   if (exist(clmfname,'file') ~= 2)
      error('%s does not exist.',clmfname)
   else
      x.clmfile = clmfname;
   end
end

if (~nc_isvar(fname,varname))
   error('%s does not have a "%s" variable.',fname,varname)
end
x.varname  = varname;

%% define and get the basic coordinate variables

variables = {'area', 'lon', 'lat', 'landfrac'};

for ivar = 1:length(variables)
   if (~nc_isvar(fname,variables{ivar}))
      error('%s does not have a "%s" variable.',fname,variables{ivar})
   else
      expr = sprintf('x.%s = nc_varget(fname,variables{ivar});',variables{ivar}); 
      eval(expr)
   end
end 

%% Must find if variable is partitioned based on land unit, column, or pft.
% Lets hope for 'column' or 'pft'.

vinfo = nc_getvarinfo(fname,varname);
ndims = length(vinfo.Dimension);
weights = [];
for idim = 1:ndims
   switch lower(vinfo.Dimension{idim})
      case {'time','copy','levgrnd'}
      case {'column'}
         weights  = nc_varget(fname,'cols1d_wtxy');
         loninds  = nc_varget(fname,'cols1d_ixy');
         latinds  = nc_varget(fname,'cols1d_jxy');
         ambigdim = idim;
      case {'pft'}
         weights  = nc_varget(fname,'pfts1d_wtxy');
         loninds  = nc_varget(fname,'pfts1d_ixy');
         latinds  = nc_varget(fname,'pfts1d_jxy');
         ambigdim = idim;
      case {'levtot','levsno'}
         % error('%s has dimension "%s" - need ZSNO',varname,vinfo.Dimension{idim})
         leveldimension = idim;
      otherwise
         error('unsupported dimension "%s"',vinfo.Dimension{idim})
   end 
end

if ( exist('ambigdim','var') ~= 1) 
   error('Could not determine weights/indexes for %s',varname)
end

%% unpack the sparse representation into a matrix.

myinfo.fname      = fname;
myinfo.timeindex  = timeindex;
myinfo.copyindex  = get_copy_index(fname,copystring);
myinfo.levelindex = levelindex;
[start,count]     = GetNCindices(myinfo,'fname',varname);
x.dat             = nc_varget(fname,varname,start,count);

%% unpack the sparse representation into a matrix.
x.datmat = zeros(size(x.area));
x.wgtmat = zeros(size(x.area));

for i = 1:length(x.dat)
   ilon = loninds(i);
   jlat = latinds(i);
   if (isfinite(x.dat(i))) 
      x.datmat(jlat,ilon) = x.datmat(jlat,ilon) + weights(i)*x.dat(i);
      x.wgtmat(jlat,ilon) = x.wgtmat(jlat,ilon) + weights(i);
   end
end

x.datmat = x.datmat./x.wgtmat;
x.datmat(x.landfrac==0.0) = NaN;

%% Get the same data from the CLM restart file and see what happens.
if (~isempty(x.clmfile))

   t_org = nc_varget(x.clmfile,varname);
   t_org(t_org > 1e+35) = NaN;

   x.org = squeeze(t_org(:,levelindex));

   x.clmres = zeros(size(x.area));
   x.wgtmat(:) = 0.0;

   for i = 1:length(x.org)
      ilon = loninds(i);
      jlat = latinds(i);
      if (isfinite(x.org(i))) 
         x.clmres(jlat,ilon) = x.clmres(jlat,ilon) + weights(i)*x.org(i);
         x.wgtmat(jlat,ilon) = x.wgtmat(jlat,ilon) + weights(i);
      end
   end

   x.clmres = x.clmres./x.wgtmat;
   x.clmres(x.landfrac==0.0) = NaN;

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

