% Given a size of the field and a variable number, determines the
% indices in the state vector array that correspond to that
% particular variable.
function indices=get_state_indices(varnum,ijarea)
  indices = ((varnum-1)*ijarea + 1):(varnum*ijarea);