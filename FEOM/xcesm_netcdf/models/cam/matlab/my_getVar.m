
function data = my_getVar(ncid,varid,FillValue)

% Usage data = my_getVar(ncid,Var_ind,FillValue);

data = netcdf.getVar(ncid, varid);
if ( ~ isempty(FillValue) )
   data(data == FillValue) = NaN;
end

