%

% <next four lines automatically updated by CVS, do not edit>
% $Source$ 
% $Revision$ 
% $Date$ 
% $Author$ 

load obs_def;
x = obs_def(:, 1);
y = obs_def(:, 2);
plot(x, y, '*')
axis([0 360 -90 90])
