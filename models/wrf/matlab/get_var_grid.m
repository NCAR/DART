function output = get_var_grid(fname, varname)
%% get_var_grid gets the grid variables for a given WRF variable.
%
% Example using a DART file: 
% fname      = 'Prior_Diag.nc';
% varname    = 'U_d01';
%
% mygrid = get_var_grid(fname, varname);
%
% Example using a WRF file: 
% fname      = 'wrfinput_d01';
% varname    = 'U';
%
% mygrid = get_var_grid(fname, varname);


%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%% Copy the global attributes of interest
filename    = fname;
varname     = varname;

if (exist(fname,'file') ~= 2)
   error('%s does not exist',fname)
end

if (nc_isvar(fname,varname) ~= 1 )
   error('%s does not exist in %s',varname,fname)
end

%% Look for a telltale DART attribute to determine DART file or WRF file.

bob = nc_info(fname);

for i = 1:length(bob.Attribute)
   if (strcmpi(bob.Attribute(i).Name, 'model'))
       chunk.DARTfile               = 1;
       chunk.Zname         = 'ZNU_d01';
       chunk.Zname_stag    = 'ZNW_d01';
       chunk.Zdimname      = 'bottom_top_d01';
       chunk.Zdimname_stag = 'bottom_top_stag_d01';
       break
   else
       chunk.DARTfile      = -1;
       chunk.Zname         = 'ZNU';
       chunk.Zname_stag    = 'ZNW';
       chunk.Zdimname      = 'bottom_top';
       chunk.Zdimname_stag = 'bottom_top_stag';
   end
end

x = nc_getdiminfo(fname,chunk.Zdimname);      chunk.bottom_top      = x.Length;
x = nc_getdiminfo(fname,chunk.Zdimname_stag); chunk.bottom_top_stag = x.Length;

%% Get the variable information so we can query the dimensions etc.
%  "old-school" DART variables do not have a memoryorder, so we must guess.

varinfo = nc_getvarinfo(filename, varname);

memoryorder = [];
for i = 1:length(varinfo.Attribute)
   attname = varinfo.Attribute(i).Name;
   switch lower(attname)
      case 'memoryorder'
         memoryorder = deblank(varinfo.Attribute(i).Value);
   end
end

% here's the guess part

if (isempty(memoryorder))

  if     ( length(varinfo.Size) == 3 ) 
     memoryorder = 'Z';
  elseif ( length(varinfo.Size) == 4 ) 
     memoryorder = 'XY';
  elseif ( length(varinfo.Size) == 5 ) 
     memoryorder = 'XYZ';
  end

end

switch lower(memoryorder)
case '0'   % Time
   error('unsupported storage order (%s) for %s',memoryorder, varname)

case 'z'   % Time, bottom_top
   error('unsupported storage order (%s) for %s',memoryorder, varname)
   
case 'xy'  % Time, south_north, west_east

   coordinates      = strread(nc_attget(fname, varname, 'coordinates'),'%s');
   output.xvarname  = char(bob(1,:));
   output.yvarname  = char(bob(2,:));
   output.xvar      = double(nc_varget(fname, output.xvarname));
   output.yvar      = double(nc_varget(fname, output.yvarname));
   
case 'xyz' % .... you get the picture

   bob = get_z_info(fname, varname, chunk);

   coordinates      = strread(nc_attget(fname, varname, 'coordinates'),'%s');
   output.xvarname  = char(coordinates(1,:));
   output.yvarname  = char(coordinates(2,:));
   output.zvarname  = bob.zvarname;

   output.zvarlabel = bob.zvarlabel;
   output.zvarunits = bob.zvarunits;

   output.xvar = double(nc_varget(fname,output.xvarname));
   output.yvar = double(nc_varget(fname,output.yvarname));
   output.zvar = double(nc_varget(fname,output.zvarname));
   
otherwise
   error('unknown storage order (%s) for %s',memoryorder, varname)
end

inds              = output.xvar < 0.0;
output.xvar(inds) = output.xvar(inds) + 360.0;


%======================================================================


function zinfo = get_z_info(fname, varname, chunk);
% There is not much help determining metadata for the Z coordinate
% Must match dimensions and make educated guesses.

   varinfo = nc_getvarinfo(fname, varname);

   bob = strfind(lower(varinfo.Dimension), 'bottom');
   zindex = -1;
   for i = 1:length(bob)
      if ( bob{i} == 1 )
         zindex = i;
         zinfo.zdimname = varinfo.Dimension{i};
         zinfo.zdimsize = varinfo.Size(i);
      end
   end

   if (zindex < 0) 
      error('dagnabbit')
   end
   if (    zinfo.zdimsize == chunk.bottom_top )
           zinfo.zvarname =  chunk.Zname;
   elseif (zinfo.zdimsize == chunk.bottom_top_stag )
           zinfo.zvarname =  chunk.Zname_stag;
   else
       error('gollygoshdagnabbit')
   end

   % get the units for the Z variable being used.
   % preserve whatever units are being used
   % preserve whatever the long_name or description is (presumably either one works)
   % use the first word of the long name as the 'label' ... usually 'eta' or 'height' or ...

   if (nc_isvar(fname,zinfo.zvarname))

      varinfo = nc_getvarinfo(fname, zinfo.zvarname);
      for i = 1:length(varinfo.Attribute)
         switch lower(varinfo.Attribute(i).Name)
            case 'units'
               zinfo.zvarunits = varinfo.Attribute(i).Value;
            case {'description','long_name'}
               zinfo.zvarlongname = varinfo.Attribute(i).Value;
               zinfo.zvarlabel = sscanf(varinfo.Attribute(i).Value,'%s',1);
         end
      end

   else
      error('unable to find vertical variable %s in %s',output.zvarname,fname)
   end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
