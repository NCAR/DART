%

% <next four lines automatically updated by CVS, do not edit>
% $Source$ 
% $Revision$ 
% $Date$ 
% $Author$ 

figure(1);
orient tall;

dir = pwd
l = size(dir)
dir_name = dir(l(2)-5:l(2))

for j = 1:4
   for k = 1:2
      n = 2 * (j - 1) + k
      outname = ['load out', num2str(n), ';'];
      outcopy = ['out = out', num2str(n), ';'];
      eval(outname);
      eval(outcopy);
      subplot(4, 2, n);
      a = size(out);
      
      plot(out(1:a(1)), abs(out(a(1)*2 + 1:a(1)*3) - out(a(1)*1 + 1:a(1)*2)), 'b');

      v = axis;
      v(1) = 960;
      v(2) = 5760;
      axis(v);
      hold on;
      plot(out(1:a(1)), sqrt(out(a(1)*3 + 1:a(1)*4)), 'r');
      cc = corrcoef(sqrt(out(a(1)*3 + 1:a(1)*4)), abs(out(a(1)*2 + 1:a(1)*3) - out(a(1)*1 + 1:a(1)*2)));
      err_mean = sum(abs(out(a(1)*2 + 1:a(1)*3) - out(a(1)*1 + 1:a(1)*2))) / a(1);
      spread_mean = sum(sqrt(out(a(1)*3 + 1:a(1)*4))) / a(1);
      title(['Correl. = ', num2str(cc(2)), '     RMS = ', num2str(err_mean), '     Spread = ', num2str(spread_mean)]);
      xlabel(['TIME   Var = ', num2str(n), '   ', dir_name])

   end
end


%for j = 1:4
%   for k = 1:2
%      n = 2 * (j - 1) + k
%      outname = ['load out', num2str(n), ';'];
%      outcopy = ['out = out', num2str(n), ';'];
%      eval(outname);
%      eval(outcopy);
%      subplot(4, 2, n);
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
