function CheckModelCompatibility(file1,file2);
% CheckModelCompatibility   tries to ensure that two netcdf files can be compared.

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

if ( exist(file1) ~= 2 )
   error(sprintf('(file1) %s does not exist.',file1))
end
if ( exist(file2) ~= 2 )
   error(sprintf('(file2) %s does not exist.',file2))
end

% Get some information from the file1
% dimensions are f1(xxxx), variables are f1{xxxx}
f1 = netcdf(file1);
tmodel      = f1.model(:);
if (isempty(tmodel)) 
   error(sprintf('%s has no ''model'' global attribute.',file1))
end

if ( VarExist(f1,'copy') ) 
   tnum_copies = length(f1('copy')); % determine # of ensemble members
else
   error(sprintf('%s has no ''copy'' dimension.',file1))
end
if ( VarExist(f1,'time') ) 
   tnum_times  = length(f1('time')); % determine # of output times
else
   error(sprintf('%s has no ''time'' dimension.',file1))
end

[tnum_vars,tdims] = ModelDimension(f1,tmodel);
if ( tnum_vars <= 0 )
   error(sprintf('Unable to determine resolution of %s.',file1))
end
close(f1); 


% Get some information from the file2
f2 = netcdf(file2);
dmodel      = f2.model(:);
if (isempty(dmodel)) 
   error(sprintf('%s has no ''model'' global attribute.',file2))
end
if (VarExist(f2,'copy')) 
   dnum_copies = length(f2('copy')); % determine # of ensemble members
else
   error(sprintf('%s has no ''copy'' dimension.',file2))
end
if (VarExist(f2,'time')) 
   dnum_times  = length(f2('time')); % determine # of output times
else
   error(sprintf('%s has no ''time'' dimension.',file2))
end
[dnum_vars,ddims] = ModelDimension(f2,dmodel);
if ( dnum_vars <= 0 )
   error(sprintf('Unable to determine resolution of %s.',file2))
end
close(f2); 

% rudimentary bulletproofing
if (strcmp(tmodel,dmodel) ~= 1)
   disp(sprintf('%s has model %s ',file1,tmodel))
   disp(sprintf('%s has model %s ',file2,dmodel))
   error('no No NO ... models must be the same')
end
if (prod(tnum_vars) ~= prod(dnum_vars))
   disp(sprintf('%s has %d state variables',file1,prod(tnum_vars)))
   disp(sprintf('%s has %d state variables',file2,prod(dnum_vars)))
   error('no No NO ... both files must have same shape of state variables.')
end
if (tnum_times ~= dnum_times)
   disp(sprintf('%s has %d timesteps',file1,tnum_times))
   disp(sprintf('%s has %d timesteps',file2,dnum_times))
   error('ugh ... both files must have same number of timesteps.')
end


function x = VarExist(ncid,varname)

x = 0;   % false ... assumed not to exist.
variables = var(ncid);
for i=1:length(variables)
   if ( strmatch(name(variables{i}), varname) == 1 )
      x = 1;   % true ... variables exists in the netcdf file.
   end
end


function x = DimExist(ncid,dimname)

x = 0;   % false ... assumed not to exist.
dimensions = dim(ncid);
for i=1:length(dimensions)
   if ( strmatch(name(dimensions{i}), dimname) == 1 )
      x = 1;   % true ... variables exists in the netcdf file.
   end
end


function [x,y] = ModelDimension(ncid,modelname)

x = 0;
y = NaN;

switch lower(modelname)

   case 'cam'
      lonexist = VarExist(ncid,'lon'  );
      latexist = VarExist(ncid,'lat'  );
      lvlexist = VarExist(ncid,'lev');
      if ( latexist && lonexist && lvlexist )
         dnum_lons = prod(size(ncid('lon')));
         dnum_lats = prod(size(ncid('lat')));
         dnum_lvls = prod(size(ncid('lev')));
         x = 3;
         y = [dnum_lons dnum_lats dnum_lvls];
      end

   case 'fms_bgrid'
      lonexist = VarExist(ncid,'TmpI'  );
      latexist = VarExist(ncid,'TmpJ'  );
      lvlexist = VarExist(ncid,'level');
      if ( latexist && lonexist && lvlexist )
         dnum_lons = prod(size(ncid('TmpI')));
         dnum_lats = prod(size(ncid('TmpJ')));
         dnum_lvls = prod(size(ncid('level')));
         x = 3;
         y = [dnum_lons dnum_lats dnum_lvls];
      end

   case 'simple_advection'
      if ( VarExist(ncid,'loc1d')) 
         y = prod(size(ncid('loc1d')));
	 x = 1;
      end

   otherwise
%     disp(sprintf('working with %s',modelname))
      if ( VarExist(ncid,'StateVariable')) 
         y = prod(size(ncid('StateVariable')));
	 x = 1;
      end

end

