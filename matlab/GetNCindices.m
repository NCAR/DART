function [start, count] = GetNCindices(pinfo, whichfile, varname)
%% GETNCindices returns a start,count array for use with nc_getvar.
% At present, all times, all copies for a specific level,lat,lon.
% Does not assume anything about the dimension of the variable.
%
% USAGE:
% [start, count] = GetNCindices(pinfo, whichfile, varname);
%
% The structure 'pinfo' must have ONE of the following components:
%                    pinfo.[prior,posterior,truth,diagn,fname]
%     and may have
%                    pinfo.timeindex
%                    pinfo.copyindex
%                    pinfo.levelindex
%                    pinfo.latindex
%                    pinfo.lonindex
%                    pinfo.stateindex
%                    pinfo.regionindex
%                    pinfo.cellindex
%                    pinfo.columnindex
%                    pinfo.pftindex
%                    pinfo.rankbinindex
%                    pinfo.boundsindex
%                    pinfo.obstypesindex
%
% whichfile          is a character string specifying which
%                    filename component of 'pinfo' will be used.
%                    ['prior','posterior','truth','diagn','fname']
% varname            is the netcdf variable being extracted.
%

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% GetNCindices replaces the following hardwired piece of code.
%
%      if ( strcmp(lower(vname),'ps') ==1 ) %  PS(time, copy, lat, lon)
%         start = [  1 NaN pinfo.latindex pinfo.lonindex];
%         count = [ -1 NaN pinfo.latindex pinfo.lonindex];
%      else % U(time, copy, lev, lat, lon)
%         start = [  1 NaN pinfo.levelindex pinfo.latindex pinfo.lonindex ];
%         count = [ -1 NaN pinfo.levelindex pinfo.latindex pinfo.lonindex ];
%      end
switch lower(whichfile)
    case 'prior'
        fname = pinfo.prior_file;
    case 'posterior'
        fname = pinfo.posterior_file;
    case 'truth'
        fname = pinfo.truth_file;
    case 'diagn'
        fname = pinfo.diagn_file;
    otherwise
        fname = pinfo.fname;
end

if ( exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

% If the structure has subsetting information, we will ultimately use it.
% Otherwise, use the whole extent.

lat1        = 0; latN        = -1;
lon1        = 0; lonN        = -1;
time1       = 0; timeN       = -1;
copy1       = 0; copyN       = -1;
level1      = 0; levelN      = -1;
state1      = 0; stateN      = -1;
region1     = 0; regionN     = -1;
cell1       = 0; cellN       = -1;
column1     = 0; columnN     = -1;
pft1        = 0; pftN        = -1;
rankbin1    = 0; rankbinN    = -1;
bounds1     = 0; boundsN     = -1;
obstypes1   = 0; obstypesN   = -1;

if (isfield(pinfo,'timeindex'))
   time1 = pinfo.timeindex - 1;
   timeN = 1;
end
if (isfield(pinfo,'tindex1'))
   time1 = pinfo.tindex1 - 1;
end
if (isfield(pinfo,'tcount'))
   timeN = pinfo.tcount;
end

if (isfield(pinfo,'levelindex'))
   level1 = pinfo.levelindex - 1;
   levelN = 1;
end
if (isfield(pinfo,'level1'))
   level1 = pinfo.level1 - 1;
end
if (isfield(pinfo,'levelcount'))
   levelN = pinfo.levelcount;
end

if (isfield(pinfo,'latindex'))
   lat1 = pinfo.latindex - 1;
   latN = 1;
end
if (isfield(pinfo,'lat1'))
   lat1 = pinfo.lat1 - 1;
end
if (isfield(pinfo,'latcount'))
   latN = pinfo.latcount;
end

if (isfield(pinfo,'lonindex'))
   lon1 = pinfo.lonindex - 1;
   lonN = 1;
end
if (isfield(pinfo,'lon1'))
   lon1 = pinfo.lon1 - 1;
end
if (isfield(pinfo,'loncount'))
   lonN = pinfo.loncount;
end

if (isfield(pinfo,'stateindex'))
   state1 = pinfo.stateindex - 1;
   stateN = 1;
end
if (isfield(pinfo,'state1'))
   state1 = pinfo.state1 - 1;
end
if (isfield(pinfo,'statecount'))
   stateN = pinfo.statecount;
end

if (isfield(pinfo,'copyindex'))
   copy1 = pinfo.copyindex - 1;
   copyN = 1;
end
if (isfield(pinfo,'copy1'))
   copy1 = pinfo.copy1 - 1;
end
if (isfield(pinfo,'copycount'))
   copyN = pinfo.copycount;
end

if (isfield(pinfo,'regionindex'))
   region1 = pinfo.regionindex - 1;
   regionN = 1;
end
if (isfield(pinfo,'region1'))
   region1 = pinfo.region1 - 1;
end
if (isfield(pinfo,'regioncount'))
   regionN = pinfo.regioncount;
end

if (isfield(pinfo,'cellindex'))
   cell1 = pinfo.cellindex - 1;
   cellN = 1;
end
if (isfield(pinfo,'cell1'))
   cell1 = pinfo.cell1 - 1;
end
if (isfield(pinfo,'cellcount'))
   cellN = pinfo.cellcount;
end

if (isfield(pinfo,'columnindex'))
   column1 = pinfo.columnindex - 1;
   columnN = 1;
end
if (isfield(pinfo,'column1'))
   column1 = pinfo.column1 - 1;
end
if (isfield(pinfo,'columncount'))
   columnN = pinfo.columncount;
end

if (isfield(pinfo,'pftindex'))
   pft1 = pinfo.pftindex - 1;
   pftN = 1;
end
if (isfield(pinfo,'pft1'))
   pft1 = pinfo.pft1 - 1;
end
if (isfield(pinfo,'pftcount'))
   pftN = pinfo.pftcount;
end

if (isfield( pinfo,'rankbinindex'))
   rankbin1 = pinfo.rankbinindex - 1;
   rankbinN = 1;
end
if (isfield( pinfo,'rankbin1'))
   rankbin1 = pinfo.rankbin1 - 1;
end
if (isfield( pinfo,'rankbincount'))
   rankbinN = pinfo.rankbincount;
end

if (isfield( pinfo,'boundsindex'))
   bounds1 = pinfo.boundsindex - 1;
   boundsN = 1;
end
if (isfield( pinfo,'bounds1'))
   bounds1 = pinfo.bounds1 - 1;
end
if (isfield( pinfo,'boundscount'))
   boundsN = pinfo.boundscount;
end

if (isfield( pinfo,'obstypesindex'))
   obstypes1 = pinfo.obstypesindex - 1;
   obstypesN = 1;
end
if (isfield( pinfo,'obstypes1'))
   obstypes1 = pinfo.obstypes1 - 1;
end
if (isfield( pinfo,'obstypescount'))
   obstypesN = pinfo.obstypescount;
end

% Determine shape of variable in question.

varinfo = nc_getvarinfo(fname,varname);
ndims   = length(varinfo.Dimension);
start   = zeros(1,ndims);
count   = zeros(1,ndims);

% varinfo.Dimension is a cell array of the Dimension strings
% varinfo.Size is an N-D array describing the variable shape
% varinfo.Attribute is a struct holding the variable attribues

% loop over all of the variables dimensions and
% build up the start/endpoint arrays

for i = 1:ndims

   diminfo = nc_getdiminfo(fname,varinfo.Dimension{i});
   dimname = diminfo.Name;

   % the coordinate variable is the same name as the dimension
   % some grids have multiple vertical levels so there is no one
   % predictable coordinate variable name. By convention, these
   % variables should be tagged with a 'cartesian_axis' attribute.
   % ditto for lat, lon ... (on staggered grids, for example)
   % So the XG coordinate dimension has 'cartesian_axis = X',
   % for example.

   [~, status, value] = is_dimension_cartesian(fname, diminfo.Name);

   if (status > 0)
      dimname = value;
   else
      % Then there is no 'cartesian_axis' attribute and the best we can
      % hope for is a standard dimension name [time,copy,lat,lon,lev]
   end

   if (length(dimname) >= 4)
       switch lower(dimname(1:4)) % loop over all likely coordinate variables
           case 'time'
               start(i) = time1;
               count(i) = timeN;
           case 'copy'
               start(i) = copy1;
               count(i) = copyN;
           case {'surf','unde','hlev','mlev','plev','heig','leve','bott', ...
                 'ilev','nver','levt','levs','dept'}
               start(i) = level1;
               count(i) = levelN;
           case {'tmpj','sout','slat','iyrt'}
               start(i) = lat1;
               count(i) = latN;
           case {'tmpi','west','slon','ixrt'}
               start(i) = lon1;
               count(i) = lonN;
           case {'stat','xdim','ydim','loc1'}
               % the lorenz_96_2scale has the unfortunate choice of
               % 'Xdim' and 'YDim' for their state variable names.
               start(i) = state1;
               count(i) = stateN;
           case 'regi'
               start(i) = region1;
               count(i) = regionN;
           case 'ncel'
               start(i) = cell1;
               count(i) = cellN;
           case 'colu'
               start(i) = column1;
               count(i) = columnN;
           case 'rank'
               start(i) = rankbin1;
               count(i) = rankbinN;
           case 'boun'
               start(i) = bounds1;
               count(i) = boundsN;
           case 'obst'
               start(i) = obstypes1;
               count(i) = obstypesN;
           otherwise
               fprintf('GetNCindices encountered unknown coordinate variable %s\n',dimname)
       end
   else
       switch lower(dimname) % loop over all likely coordinate variables
           case 't'
               start(i) = time1;
               count(i) = timeN;
           case {'lev','z'}
               start(i) = level1;
               count(i) = levelN;
           case {'lat','y','iy'}
               start(i) = lat1;
               count(i) = latN;
           case {'lon','x','ix'}
               start(i) = lon1;
               count(i) = lonN;
           case 'pft'
               start(i) = pft1;
               count(i) = pftN;
           otherwise
               fprintf('GetNCindices encountered unknown coordinate variable %s\n',dimname)
       end
   end
end

count(count < 1) = -1;



function [len, status, value] = is_dimension_cartesian(fname,dimname)

status   = 0;
len      = 0;
value    = [];

if ( nc_isvar(fname,dimname) )

    % Good - the coordinate variable exists.

    Cvarinfo = nc_getvarinfo(fname, dimname);

   for j = 1:length(Cvarinfo.Attribute);
      attribute = Cvarinfo.Attribute(j);
      switch lower(attribute.Name)
         case{'cartesian_axis'}
            status = 1;
            len    = Cvarinfo.Size;
            value  = attribute.Value;
            break
          otherwise
      end
   end

else % there is no coordinate variable ... use something else

end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

