function oba=om(i,th)

% DART $Id$
% CREDIT: Alexey Morozov

if i==1
      oba=[1 0 0
           0 cos(th) sin(th)
           0 -sin(th) cos(th)];
elseif i==2
      oba=[cos(th) 0 -sin(th)
           0 1 0
           sin(th) 0 cos(th)];
elseif i==3
    oba=[cos(th) sin(th) 0
        -sin(th) cos(th) 0
        0 0 1];
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
