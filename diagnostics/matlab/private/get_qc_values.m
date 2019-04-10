function qcvalues = get_qc_values(fname, varname, varargin)
%%

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

defaultRegion   = [];
defaultLevel    = [];
defaultTimeStep = [];
defaultVerbose  = false;
defaultFatal    = true;

p = inputParser;
addRequired(p,'fname',@ischar);
addRequired(p,'varname',@ischar);

if (exist('inputParser/addParameter','file') == 2)
    addParameter(p,'regionindex',defaultRegion,@isnumeric);
    addParameter(p,'levelindex',defaultLevel,@isnumeric);
    addParameter(p,'timeindex',defaultTimeStep,@isnumeric);
    addParameter(p,'verbose',defaultVerbose,@islogical);
    addParameter(p,'fatal',defaultFatal,@islogical);
else
    addParamValue(p,'regionindex',defaultRegion,@isnumeric);
    addParamValue(p,'levelindex',defaultLevel,@isnumeric);
    addParamValue(p,'timeindex',defaultTimeStep,@isnumeric);
    addParamValue(p,'verbose',defaultVerbose,@islogical);
    addParamValue(p,'fatal',defaultFatal,@islogical);
end

p.parse(fname, varname, varargin{:});

myinfo.diagn_file = p.Results.fname;

regionstring = 'all regions';
levelstring  = 'all levels';
timestring   = 'all times';

if ~isempty(p.Results.regionindex)
    myinfo.regionindex = p.Results.regionindex;
    regionstring = sprintf('region %d',myinfo.regionindex);
end

if ~isempty(p.Results.levelindex)
    myinfo.levelindex = p.Results.levelindex;
    levelstring = sprintf('level %d',myinfo.levelindex);
end

if ~isempty(p.Results.timeindex)
    myinfo.timeindex = p.Results.timeindex;
    timestring = sprintf('time %d',myinfo.timeindex);
end

qcvalues = struct('fname',fname,'varname',varname);

myinfo.copyindex  = get_copy_index(fname, 'Nposs');
[start, count]    = GetNCindices(myinfo,'diagn',varname);
nposs             = ncread(fname, varname, start, count);

myinfo.copyindex  = get_copy_index(fname, 'Nused');
[start, count]    = GetNCindices(myinfo,'diagn',varname);
nused             = ncread(fname, varname, start, count);

for qc = 0:8
    myinfo.copyindex = get_copy_index(fname, sprintf('N_DARTqc_%d',qc),'fatal',p.Results.fatal);
    
    if (myinfo.copyindex > 0)
        [start, count]   = GetNCindices(myinfo,'diagn',varname);
        cmd = sprintf('qcvalues.Nqc%d = ncread(fname, varname, start, count);',qc);
        eval(cmd)
    end
end

qcvalues.nposs         = nposs - qcvalues.Nqc5 - qcvalues.Nqc6;
qcvalues.nused         = nused;
qcvalues.num_evaluated = qcvalues.Nqc1 + qcvalues.Nqc3;

%=====================================================================

if (~ p.Results.verbose), return; end

%=====================================================================

if (exist('levelstring','var')==1)
    
    fprintf('\n %s %s %s %s\n',varname, regionstring, levelstring, timestring);
    
end

fprintf('DART QC == 0, n = %d\n',sum(qcvalues.Nqc0(:)))
fprintf('DART QC == 1, n = %d\n',sum(qcvalues.Nqc1(:)))
fprintf('DART QC == 2, n = %d\n',sum(qcvalues.Nqc2(:)))
fprintf('DART QC == 3, n = %d\n',sum(qcvalues.Nqc3(:)))
fprintf('DART QC == 4, n = %d\n',sum(qcvalues.Nqc4(:)))
fprintf('DART QC == 5, n = %d\n',sum(qcvalues.Nqc5(:)))
fprintf('DART QC == 6, n = %d\n',sum(qcvalues.Nqc6(:)))
fprintf('DART QC == 7, n = %d\n',sum(qcvalues.Nqc7(:)))
if (isfield(qcvalues,'Nqc8'))
    fprintf('DART QC == 8, n = %d\n',sum(qcvalues.Nqc8(:)))
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
