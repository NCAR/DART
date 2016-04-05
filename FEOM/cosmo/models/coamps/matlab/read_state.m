function data=read_state(ncFileID,times,member,elements)
%% data = read_state(ncFileID, times, member, elements)
%
% Given a NetCDF file handle, the times in question, the ensemble
% member in question, and the elements to read, reads in data from
% a DART NetCDF file

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

data=read_field(ncFileID,times,member,elements,'state');
