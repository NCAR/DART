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


%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

%% Copy the global attributes of interest
filename    = fname;
varname     = varname;
memoryorder = deblank(nc_attget(fname,varname,'MemoryOrder'));

if (exist(fname,'file') ~= 2)
   error('%s does not exist',fname)
end

if (nc_isvar(fname,varname) ~= 1 )
   error('%s does not exist in %s',varname,fname)
end

bob      = nc_info(fname);
DARTfile = -1;

%% Look for a telltale DART attribute to determine DART file or WRF file.

for i = 1:length(bob.Attribute)
   if (strcmpi(bob.Attribute(i).Name, 'model'))
       DARTfile               = 1;
       Zname         = 'ZNU_d01';
       Zname_stag    = 'ZNW_d01';
       Zdimname      = 'bottom_top_d01';
       Zdimname_stag = 'bottom_top_stag_d01';
       break
   else
       Zname         = 'ZNU';
       Zname_stag    = 'ZNW';
       Zdimname      = 'bottom_top';
       Zdimname_stag = 'bottom_top_stag';
   end
end

x = nc_getdiminfo(fname,Zdimname);      bottom_top      = x.Length;
x = nc_getdiminfo(fname,Zdimname_stag); bottom_top_stag = x.Length;

%% Get the variable information so we can query the dimensions etc.

varinfo = nc_getvarinfo(filename, varname);

switch lower(memoryorder)
case '0'   % Time
   error('unsupported storage order (%s) for %s',memoryorder, varname)

case 'z'   % Time, bottom_top
   error('unsupported storage order (%s) for %s',memoryorder, varname)
   
case 'xy'  % Time, south_north, west_east
   coordinates     = nc_attget(fname, varname, 'coordinates');
   bob             = strread(coordinates,'%s');
   output.xvarname = char(bob(1,:));
   output.yvarname = char(bob(2,:));
   output.xvar     = nc_varget(fname, output.xvarname);
   output.yvar     = nc_varget(fname, output.yvarname);
   
case 'xyz' % .... you get the picture
   coordinates     = nc_attget(fname, varname, 'coordinates');
   bob             = strread(coordinates,'%s');
   output.xvarname = char(bob(1,:));
   output.yvarname = char(bob(2,:));

   % WRF doesn't provide the same mechanism for a Z coordinate
   % find dimension index for vert
   bob = strfind(lower(varinfo.Dimension), 'bottom') ;
   zindex = -1;
   for i = 1:length(bob)
      if ( bob{i} == 1 )
         zindex = i;
         output.zdimname = varinfo.Dimension{i};
         output.zdimsize = varinfo.Size(i);
      end
   end
   if (zindex < 0) 
      error('dagnabbit')
   end
   if (    output.zdimsize == bottom_top )
           output.zvarname =  Zname;
   elseif (output.zdimsize == bottom_top_stag )
           output.zvarname =  Zname_stag;
   else
       error('gollygoshdagnabbit')
   end

   output.xvar      = nc_varget(fname,output.xvarname);
   output.yvar      = nc_varget(fname,output.yvarname);
   output.zvar      = nc_varget(fname,output.zvarname);
   
otherwise
   error('unknown storage order (%s) for %s',memoryorder, varname)
end

inds              = output.xvar < 0;
output.xvar(inds) = output.xvar(inds) + 360;
output.minlon     = min(output.xvar(:));
output.maxlon     = max(output.xvar(:));

%% Set the latitude limits - must check for southern hemisphere

absmin            = min(abs(output.yvar(:)));
output.minlat     = min(output.yvar(:));
output.maxlat     = max(output.yvar(:));
