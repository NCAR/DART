figure(1);
load bin;
c = size(bin);
hold off;
nbins = c(1);

for i = 1:1
   x = bin(nbins*i + 1 : nbins*(i+1));
   bar(x, 'b');
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
%axes('position', [0.1, 0.1, 0.8, 0.4]);
hold off;

for i = 4:a(2) - 1
% Plot the individual ensemble members
%   plot(out(1:a(1)), out(a(1)*(i)+1:a(1)*(i+1)), '+m')
   plot(out(1:a(1)), out(a(1)*(i)+1:a(1)*(i+1)), 'b')
   hold on;
end
% Plot the verification as a thicker line
%plot(out(1:a(1)), out(a(1) + 1:a(1) + a(1)), '*g');
q = plot(out(1:a(1)), out(a(1) + 1:a(1) + a(1)), 'k');
set(q, 'LineWidth', [4.0]);

% Plot the ensemble mean as a dashed red thicker line
r = plot(out(1:a(1)), out(a(1)*2 + 1:a(1)*3), '--k');
set(r, 'LineWidth', [4.0]);

obsname = ['load obs', num2str(n), ';'];
obscopy = ['obs = obs', num2str(n), ';'];
eval(obsname);
eval(obscopy);

b = size(obs);
%plot(obs(1:b(1)), obs(b(1) + 1: 2*b(1)), '*b')

xlabel('TIME')


%figure(3);
%hold off;

%plot(out(1:a(1)), out(a(1) + 1:a(1) + a(1)), '*g');
%hold on;
%plot(out(1:a(1)), out(a(1)*2 + 1:a(1)*3), '+r');
%%plot(obs(1:b(1)), obs(b(1) + 1: 2*b(1)), '*b')


figure(4);
axes('position', [0.1, 0.1, 0.8, 0.4]);
hold off;

q = plot(out(1:a(1)), sqrt(out(a(1)*3 + 1:a(1)*4)));
set(q, 'LineWidth', [2.0]);
hold on;
r = plot(out(1:a(1)), abs(out(a(1)*2 + 1:a(1)*3) - out(a(1)*1 + 1:a(1)*2)), '--');
set(r, 'LineWidth', [2.0]);
cc = corrcoef(sqrt(out(a(1)*3 + 1:a(1)*4)), abs(out(a(1)*2 + 1:a(1)*3) - out(a(1)*1 + 1:a(1)*2)));
title(['Correlation coefficient ', num2str(cc(2))]);
xlabel('TIME')

