function [start, count] = GetNCindices(pinfo, whichfile, varname)
%% GETNCindices returns a start,count array for use with ncread.
% At present, all times, all copies for a specific level,lat,lon.
% Does not assume anything about the dimension of the variable.
%
% ncread uses 1-based addressing for the 'start', 
% rather than the 0-based system used in C-style indexing.
%
% USAGE:
% [start, count] = GetNCindices(pinfo, whichfile, varname);
%
% The structure 'pinfo' must have ONE of the following components:
%                    pinfo.[prior,posterior,truth,diagn,fname]
%     and may have
%                    pinfo.timeindex
%                    pinfo.memberindex
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

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% GetNCindices replaces the following hardwired piece of code.
%
%      if ( strcmp(lower(vname),'ps') ==1 ) %  PS(time, member, lat, lon)
%         start = [  1 NaN pinfo.latindex pinfo.lonindex];
%         count = [Inf NaN pinfo.latindex pinfo.lonindex];
%      else % U(time, member, lev, lat, lon)
%         start = [  1 NaN pinfo.levelindex pinfo.latindex pinfo.lonindex ];
%         count = [Inf NaN pinfo.levelindex pinfo.latindex pinfo.lonindex ];
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

lat1        = 1; latN        = Inf;
lon1        = 1; lonN        = Inf;
time1       = 1; timeN       = Inf;
copy1       = 1; copyN       = Inf;
member1     = 1; memberN     = Inf;
level1      = 1; levelN      = Inf;
state1      = 1; stateN      = Inf;
region1     = 1; regionN     = Inf;
cell1       = 1; cellN       = Inf;
column1     = 1; columnN     = Inf;
pft1        = 1; pftN        = Inf;
rankbin1    = 1; rankbinN    = Inf;
bounds1     = 1; boundsN     = Inf;
obstypes1   = 1; obstypesN   = Inf;

if (isfield(pinfo,'timeindex'))
   time1 = pinfo.timeindex;
   timeN = 1;
end
if (isfield(pinfo,'tindex1'))
   time1 = pinfo.tindex1;
end
if (isfield(pinfo,'tcount'))
   timeN = pinfo.tcount;
end

if (isfield(pinfo,'levelindex'))
   level1 = pinfo.levelindex;
   levelN = 1;
end
if (isfield(pinfo,'level1'))
   level1 = pinfo.level1;
end
if (isfield(pinfo,'levelcount'))
   levelN = pinfo.levelcount;
end

if (isfield(pinfo,'latindex'))
   lat1 = pinfo.latindex;
   latN = 1;
end
if (isfield(pinfo,'lat1'))
   lat1 = pinfo.lat1;
end
if (isfield(pinfo,'latcount'))
   latN = pinfo.latcount;
end

if (isfield(pinfo,'lonindex'))
   lon1 = pinfo.lonindex;
   lonN = 1;
end
if (isfield(pinfo,'lon1'))
   lon1 = pinfo.lon1;
end
if (isfield(pinfo,'loncount'))
   lonN = pinfo.loncount;
end

if (isfield(pinfo,'stateindex'))
   state1 = pinfo.stateindex;
   stateN = 1;
end
if (isfield(pinfo,'state1'))
   state1 = pinfo.state1;
end
if (isfield(pinfo,'statecount'))
   stateN = pinfo.statecount;
end

if (isfield(pinfo,'copyindex'))
   copy1 = pinfo.copyindex;
   copyN = 1;
end
if (isfield(pinfo,'copy1'))
   copy1 = pinfo.copy1;
end
if (isfield(pinfo,'copycount'))
   copyN = pinfo.copycount;
end

if (isfield(pinfo,'memberindex'))
   member1 = pinfo.memberindex;
   memberN = 1;
end
if (isfield(pinfo,'member1'))
   member1 = pinfo.member1;
end
if (isfield(pinfo,'membercount'))
   memberN = pinfo.membercount;
end

if (isfield(pinfo,'regionindex'))
   region1 = pinfo.regionindex;
   regionN = 1;
end
if (isfield(pinfo,'region1'))
   region1 = pinfo.region1;
end
if (isfield(pinfo,'regioncount'))
   regionN = pinfo.regioncount;
end

if (isfield(pinfo,'cellindex'))
   cell1 = pinfo.cellindex;
   cellN = 1;
end
if (isfield(pinfo,'cell1'))
   cell1 = pinfo.cell1;
end
if (isfield(pinfo,'cellcount'))
   cellN = pinfo.cellcount;
end

if (isfield(pinfo,'columnindex'))
   column1 = pinfo.columnindex;
   columnN = 1;
end
if (isfield(pinfo,'column1'))
   column1 = pinfo.column1;
end
if (isfield(pinfo,'columncount'))
   columnN = pinfo.columncount;
end

if (isfield(pinfo,'pftindex'))
   pft1 = pinfo.pftindex;
   pftN = 1;
end
if (isfield(pinfo,'pft1'))
   pft1 = pinfo.pft1;
end
if (isfield(pinfo,'pftcount'))
   pftN = pinfo.pftcount;
end

if (isfield( pinfo,'rankbinindex'))
   rankbin1 = pinfo.rankbinindex;
   rankbinN = 1;
end
if (isfield( pinfo,'rankbin1'))
   rankbin1 = pinfo.rankbin1;
end
if (isfield( pinfo,'rankbincount'))
   rankbinN = pinfo.rankbincount;
end

if (isfield( pinfo,'boundsindex'))
   bounds1 = pinfo.boundsindex;
   boundsN = 1;
end
if (isfield( pinfo,'bounds1'))
   bounds1 = pinfo.bounds1;
end
if (isfield( pinfo,'boundscount'))
   boundsN = pinfo.boundscount;
end

if (isfield( pinfo,'obstypesindex'))
   obstypes1 = pinfo.obstypesindex;
   obstypesN = 1;
end
if (isfield( pinfo,'obstypes1'))
   obstypes1 = pinfo.obstypes1;
end
if (isfield( pinfo,'obstypescount'))
   obstypesN = pinfo.obstypescount;
end

% Determine shape of variable in question.

varinfo = ncinfo(fname,varname);
ndims   = length(varinfo.Dimensions);
start   = zeros(1,ndims);
count   = zeros(1,ndims);

% varinfo.Dimension is a cell array of the Dimension strings
% varinfo.Size is an N-D array describing the variable shape
% varinfo.Attribute is a struct holding the variable attribues

% loop over all of the variables dimensions and
% build up the start/endpoint arrays

for i = 1:ndims

   dimname = varinfo.Dimensions(i).Name;

   % the coordinate variable is the same name as the dimension
   % some grids have multiple vertical levels so there is no one
   % predictable coordinate variable name. By convention, these
   % variables should be tagged with a 'cartesian_axis' attribute.
   % ditto for lat, lon ... (on staggered grids, for example)
   % So the XG coordinate dimension has 'cartesian_axis = X',
   % for example.

   [~, status, value] = is_dimension_cartesian(fname, dimname);

   if (status > 0)
      dimname = value;
   else
      % Then there is no 'cartesian_axis' attribute and the best we can
      % hope for is a standard dimension name [time,member,lat,lon,lev]
   end

   if (length(dimname) >= 4)
       switch lower(dimname(1:4)) % loop over all likely coordinate variables
           case 'time'
               start(i) = time1;
               count(i) = timeN;
           case 'copy'
               start(i) = copy1;
               count(i) = copyN;
           case 'memb'
               start(i) = member1;
               count(i) = memberN;
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
           case {'stat','xloc','yloc','loca'}
               % the lorenz_96_2scale has the unfortunate choice of
               % 'XLocation' and 'YLocation' for their state variable names.
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
           case {'lev','z','k'}
               start(i) = level1;
               count(i) = levelN;
           case {'lat','y','iy','j'}
               start(i) = lat1;
               count(i) = latN;
           case {'lon','x','ix','i'}
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

%% helper function

function [len, status, value] = is_dimension_cartesian(fname,dimname)

len      = 0;
status   = 0;
value    = [];

fileinfo = ncinfo(fname);
nvars    = length(fileinfo.Variables);

for ivar = 1:nvars
    
    if strcmp(fileinfo.Variables(ivar).Name,dimname)
        
        natts = length(fileinfo.Variables(ivar).Attributes);
        
        for j = 1:natts
            attribute = fileinfo.Variables(ivar).Attributes(j);
            switch lower(attribute.Name)
                case{'cartesian_axis'}
                    status = 1;
                    len    = fileinfo.Variables(ivar).Size;
                    value  = attribute.Value;
                    return
                otherwise
            end
        end
        % if we found the coordinate variable and exhaused the atts
        % there is no point continuing
        break 
        
    else % there is no coordinate variable ... use something else
        
    end
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
