% Data Assimilation Research Testbed -- DART
% Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
% Matlab code to do fig A for paper


% <next four lines automatically updated by CVS, do not edit>
% $Source$ 
% $Revision$ 
% $Date$ 
% $Author$ 


% This is schematic text of filter impacts

figure(1)
clf
subplot(2, 2, 1)
x = [1.5];
y = [1.5];
text(x, y, 'A')
hold on;
x = [1.9];
y = [1.9];
text(x, y, 'B')
x = [2.3];
y = [2.3];
text(x, y, 'C')
x = [3.2];
y = [3.2];
text(x, y, 'D')
x = [7.7];
y = [7.7];
text(x, y, 'E')
x = [8.1];
y = [8.1];
text(x, y, 'F')
x = [8.6];
y = [8.6];
text(x, y, 'G')
x = [9.2];
y = [9.2];
text(x, y, 'H')
xlabel('X1')
ylabel('X2')
% Put on quadrant lines
x = [5 5]
y = [0 10]
plot(x, y, '--')
x = [0 10]
y = [5 5]
plot(x, y, '--')
x = [10 10]
y = [0 10]
plot(x, y)
x = [0 10]
y = [10 10]
plot(x, y)



title('(a)')
axis([0, 10, 0, 10])

% Plot a normal observational error density, too
x = 1:0.1:10
y = -1.0 * (x - 5.0) .* (x - 5.0) ./ 4
z = 2 .* 2.71828 .^ y
plot(x, z)

% Pattern for using kernel assim from AA
subplot(2, 2, 3)
x = [6.5]
y = [3.2]
text(x, y, 'A')
hold on;
x = [7.7]
y = [7.2]
text(x, y, 'B')
x = [7.2]
y = [2.5]
text(x, y, 'C')
x = [2.8]
y = [4.0]
text(x, y, 'D')
x = [2.5]
y = [8.2]
text(x, y, 'E')
x = [4.0]
y = [6.5]
text(x, y, 'F')
x = [8.2]
y = [2.8]
text(x, y, 'G')
x = [3.2]
y = [7.7]
text(x, y, 'H')
xlabel('X1')
ylabel('X2')
% Put on quadrant lines
x = [5 5]
y = [0 10]
plot(x, y, '--')
x = [0 10]
y = [5 5]
plot(x, y, '--')
x = [10 10]
y = [0 10]
plot(x, y)
x = [0 10]
y = [10 10]
plot(x, y)
title('(b)')
axis([0, 10, 0, 10])


% Pattern from single Gaussian AA assim
subplot(2, 2, 2)
x = [5.5]
y = [7.5]
text(x, y, 'A')
hold on;
x = [3.9]
y = [5.9]
text(x, y, 'B')
x = [5.2]
y = [4.2]
text(x, y, 'C')
x = [3.2]
y = [4.4]
text(x, y, 'D')
x = [5.7]
y = [8.0]
text(x, y, 'E')
x = [4.9]
y = [5.3]
text(x, y, 'F')
x = [6.3]
y = [4.9]
text(x, y, 'G')
x = [3.9]
y = [7.1]
text(x, y, 'H')
xlabel('X1')
ylabel('X2')
% Put on quadrant lines
x = [5 5]
y = [0 10]
plot(x, y, '--')
x = [0 10]
y = [5 5]
plot(x, y, '--')
x = [10 10]
y = [0 10]
plot(x, y)
x = [0 10]
y = [10 10]
plot(x, y)
title('(c)')
axis([0, 10, 0, 10])


subplot(2, 2, 4)
x = [2.5]
y = [2.5]
text(x, y, 'A')
hold on;
x = [2.8]
y = [2.8]
text(x, y, 'B')
x = [3.2]
y = [3.2]
text(x, y, 'C')
x = [4.0]
y = [4.0]
text(x, y, 'D')
x = [6.5]
y = [6.5]
text(x, y, 'E')
x = [7.2]
y = [7.2]
text(x, y, 'F')
x = [7.7]
y = [7.7]
text(x, y, 'G')
x = [8.2]
y = [8.2]
text(x, y, 'H')
xlabel('X1')
ylabel('X2')
% Put on quadrant lines
x = [5 5]
y = [0 10]
plot(x, y, '--')
x = [0 10]
y = [5 5]
plot(x, y, '--')
x = [10 10]
y = [0 10]
plot(x, y)
x = [0 10]
y = [10 10]
plot(x, y)
title('(d)')
axis([0, 10, 0, 10])



