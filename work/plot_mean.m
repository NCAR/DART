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
   subplot(2, 3, i), bar(x, 'b');
   title (['True Variable ', num2str(i)])
   xlabel BIN
   ylabel FREQUENCY
end



hold on;
grid on;


n = input('What variable to plot')
outname = ['load out', num2str(n), ';'];
outcopy = ['out = out', num2str(n), ';'];
eval(outname);
eval(outcopy);

a = size(out);
figure(2);
hold off;

for i = 4:a(2) - 1
   plot(out(1:a(1)), out(a(1)*(i)+1:a(1)*(i+1)), '+m')
   hold on;
end
plot(out(1:a(1)), out(a(1) + 1:a(1) + a(1)), '*g');

obsname = ['load obs', num2str(n), ';'];
obscopy = ['obs = obs', num2str(n), ';'];
eval(obsname);
eval(obscopy);

b = size(obs);
%plot(obs(1:b(1)), obs(b(1) + 1: 2*b(1)), '*b')

figure(3);
hold off;

plot(out(1:a(1)), out(a(1) + 1:a(1) + a(1)), '*g');
hold on;
plot(out(1:a(1)), out(a(1)*2 + 1:a(1)*3), '+r');
plot(obs(1:b(1)), obs(b(1) + 1: 2*b(1)), '*b')


figure(4);
hold off;

plot(out(1:a(1)), sqrt(out(a(1)*3 + 1:a(1)*4)), 'r');
hold on;
plot(out(1:a(1)), abs(out(a(1)*2 + 1:a(1)*3) - out(a(1)*1 + 1:a(1)*2)), 'b');
cc = corrcoef(sqrt(out(a(1)*3 + 1:a(1)*4)), abs(out(a(1)*2 + 1:a(1)*3) - out(a(1)*1 + 1:a(1)*2)));
title(['Correlation coefficient ', num2str(cc(2))]);

