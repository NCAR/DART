function indices=get_state_indices(varnum,ijarea)
% Given a size of the field and a variable number, determines the
% indices in the state vector array that correspond to that
% particular variable.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

  indices = ((varnum-1)*ijarea + 1):(varnum*ijarea);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
