% Given a NetCDF file handle, the times in question, the ensemble
% member in question, and the elements to read, reads in data from
% a DART NetCDF file
function data=read_field(ncFileID,times,member,elements,variable)
  data = squeeze(ncFileID{variable}(times,member,elements));
