%

% Data Assimilation Research Testbed -- DART
% Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next four lines automatically updated by CVS, do not edit>
% $Source$ 
% $Revision$ 
% $Date$ 
% $Author$ 

figure(1);
load bin;
c = size(bin);
hold off;
nbins = c(1);

for i = 1:min(6, c(2)-1)
   x = bin(nbins*i + 1 : nbins*(i+1));
   subplot(2, 3, i), bar(x);
   title (['True Variable ', num2str(i)])
end



hold on;
grid on;



figure(2);
load bin_obs;
d = size(bin_obs);
hold off;
nbins = d(1);

for i = 1:min(6, d(2)-1)
   x = bin_obs(nbins*i + 1 : nbins*(i+1));
   subplot(2, 3, i), bar(x);
   title (['Obs. Variable ', num2str(i)])
end



n = input('What variable to plot')
outname = ['load out', num2str(n), ';'];
outcopy = ['out = out', num2str(n), ';'];
eval(outname);
eval(outcopy);

a = size(out);
figure(3);
hold off;

for i = 1:a(2) - 1
   plot(out(a(1)*i+1:a(1)*(i+1)), '+m')
   hold on;
end
plot(out(1:a(1)), '*');

obsname = ['load obs', num2str(n), ';'];
obscopy = ['obs = obs', num2str(n), ';'];
eval(obsname);
eval(obscopy);

b = size(obs);
plot(obs(1:b(1)), obs(b(1) + 1: 2*b(1)), '*b')

