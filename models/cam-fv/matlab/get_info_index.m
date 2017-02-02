function get_info_index = get_info_index(info,kind,field)
% field is a cell array (size 1) of character(s).

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

get_info_index = -999;

if (kind == 'fld') 
   dum = fprintf('get_info_index variables %s %s %s %s %s %s %s \n',info.Variables.Name);
   num = length(info.Variables) 
elseif (kind == 'dim')
   dum = fprintf('get_info_index dimensions %s %s %s %s %s %s %s \n',info.Dimensions.Name);
   num = length(info.Dimensions)
else
   dum = fprintf('input variable "kind" must equal "fld" or "dim"')
end

for i = 1:num
    if (kind == 'fld') 
       % msg = ['Does' field ' match ' info.Variables(i).Name]
       % info_Var_Name = info.Variables(i).Name
       % field
       same = strcmp(info.Variables(i).Name,char(field)) ;
    elseif (kind == 'dim')
       same = strcmp(info.Dimensions(i).Name,char(field)) ;
    end
    if (same == 1)
       get_info_index = i ;
       % NetCDF functions uses C indexing
       if (kind == 'fld') 
          get_info_index = get_info_index - 1 ;
       end
       break
    end
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
