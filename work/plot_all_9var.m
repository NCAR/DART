% Data Assimilation Research Testbed -- DART
% Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
% <next four lines automatically updated by CVS, do not edit>
% $Source$ 
% $Revision$ 
% $Date$ 
% $Author$ 

direct = pwd

orient tall;
figure(1);

load bin;
c = size(bin);
hold off;
nbins = c(1);

for i = 1:min(6, c(2)-1)
   x = bin(nbins*i + 1 : nbins*(i+1));
   subplot(3, 2, i), bar(x, 'b');
   title (['Var ', num2str(i), direct])
   xlabel BIN
   ylabel FREQUENCY
end
print -deps fig1.eps

figure(2);
orient tall;

dir = pwd
l = size(dir)
dir_name = dir(l(2)-5:l(2))

for j = 1:3
   for k = 1:2
      n = 2 * (j - 1) + k
      outname = ['load out', num2str(n), ';'];
      outcopy = ['out = out', num2str(n), ';'];
      eval(outname);
      eval(outcopy);
%      subplot(3, 2, n);
      subplot(6, 1, n);
      a = size(out);
      
      plot(out(1:a(1)), abs(out(a(1)*2 + 1:a(1)*3) - out(a(1)*1 + 1:a(1)*2)), 'b');

      v = axis;
%      v(1) = 960;
%      v(2) = 5760;
%      axis(v);
      hold on;
      plot(out(1:a(1)), sqrt(out(a(1)*3 + 1:a(1)*4)), 'r--');
      cc = corrcoef(sqrt(out(a(1)*3 + 1:a(1)*4)), abs(out(a(1)*2 + 1:a(1)*3) - out(a(1)*1 + 1:a(1)*2)));
      err_mean = sum(abs(out(a(1)*2 + 1:a(1)*3) - out(a(1)*1 + 1:a(1)*2))) / a(1);
      spread_mean = sum(sqrt(out(a(1)*3 + 1:a(1)*4))) / a(1);
      base_ratio = err_mean / spread_mean
      desired = sqrt(11/20)
      final_ratio = base_ratio / desired


      title(['CC = ', num2str(cc(2)), '     RMS = ', num2str(err_mean), '     Spread = ', num2str(spread_mean), ' ratio = ', num2str(final_ratio)]);
%      xlabel(['TIME   Var = ', num2str(n), '   ', dir_name])

   end
end

xlabel(['RMS solid and Spread dashed  ', direct]);
print -depsc fig2.eps


%figure(2);
%for j = 1:3
%   for k = 1:2
%      n = 2 * (j - 1) + k
%      outname = ['load out', num2str(n), ';'];
%      outcopy = ['out = out', num2str(n), ';'];
%      eval(outname);
%      eval(outcopy);
%      subplot(3, 2, n);
%      a = size(out);
%      for i = 4:a(2) - 1
%% Plot the individual ensemble members
%         plot(out(1:a(1)), out(a(1)*(i)+1:a(1)*(i+1)), 'm')
%         hold on;
%      end
%% Plot the verification as a thicker line
%      q = plot(out(1:a(1)), out(a(1) + 1:a(1) + a(1)), 'g');
%      set(q, 'LineWidth', [4.0]);
%
%% Plot the ensemble mean as a dashed red thicker line
%      r = plot(out(1:a(1)), out(a(1)*2 + 1:a(1)*3), '--b');
%      set(r, 'LineWidth', [4.0]);
%      xlabel('TIME')
%   end
%end




% Figure for error time series
figure(3);
orient tall;
for j = 1:3
   for k = 1:2
      n = 2 * (j - 1) + k
      outname = ['load out', num2str(n), ';'];
      outcopy = ['out = out', num2str(n), ';'];
      eval(outname);
      eval(outcopy);
      subplot(6, 1, n);
%      subplot(3, 2, n);
      a = size(out);
      for i = 4:a(2) - 1
%% Plot the individual ensemble members difference from truth
         plot(out(1:a(1)), out(a(1)*(i)+1:a(1)*(i+1)) - out(a(1) + 1:a(1) + a(1)), 'm--')
         hold on;
      end
%% Plot the ensemble mean error as a dashed red thicker line
      r = plot(out(1:a(1)), out(a(1)*2 + 1:a(1)*3) - out(a(1) + 1:a(1) + a(1)), 'b');
      set(r, 'LineWidth', [2.0]);
% Plot a zero line for comparison
      r = plot(out(1:a(1)), out(a(1)*2 + 1:a(1)*3) - out(a(1)*2 + 1:a(1)*3), 'g-.');
%      xlabel('TIME')
   end
end

xlabel(['Ensemble and mean (solid)', direct]);
print -depsc fig3.eps

% Figure for gravity wave amplitude
figure(4);
orient tall;
for j = 1:3
   for k = 1:2
      n = 2 * (j - 1) + k
      stuffname = ['load stuff', num2str(n), ';'];
      stuffcopy = ['stuff = stuff', num2str(n), ';'];
      eval(stuffname);
      eval(stuffcopy);
%      subplot(3, 2, n);
      subplot(6, 1, n);
%      a = size(stuff);
      r = plot(stuff(:, 1), stuff(:, 4), 'b');
      hold;
% Move the 1 value of rat to 0 to allow easy overlay
      rat = stuff(:, 7) - 1;
% Want max of rat to be same as max of wave amplitude
      wave_max = max(stuff(:, 4));
      rat_max = max(rat);
      rat = rat / rat_max * wave_max;
      r = plot(stuff(:, 1), rat, 'r--');
%      xlabel('TIME')

   end
end

xlabel(['Gravity wave amplitude (solid) and normalized ratio - 1 ', direct]);
print -depsc fig4.eps
