function dates = check_reanalysis_times(year)

directory = '/glade/p/nsc/ncis0006/Reanalyses/f.e21.FHIST_BGC.f09_025.CAM6assim.011';

bases = { 'cpl/hist/0080/f.e21.FHIST_BGC.f09_025.CAM6assim.011.cpl_0080.ha2x3h', ...
          'cpl/hist/0080/f.e21.FHIST_BGC.f09_025.CAM6assim.011.cpl_0080.ha2x1hi', ...
          'cpl/hist/0080/f.e21.FHIST_BGC.f09_025.CAM6assim.011.cpl_0080.ha2x1d', ...
          'cpl/hist/0080/f.e21.FHIST_BGC.f09_025.CAM6assim.011.cpl_0080.ha2x1h', ...
          'cpl/hist/0080/f.e21.FHIST_BGC.f09_025.CAM6assim.011.cpl_0080.hr2x', ...
          'rof/hist/0080/f.e21.FHIST_BGC.f09_025.CAM6assim.011.mosart_0080.h0', ...
          'lnd/hist/0080/f.e21.FHIST_BGC.f09_025.CAM6assim.011.clm2_0080.h1', ...
          'lnd/hist/0080/f.e21.FHIST_BGC.f09_025.CAM6assim.011.clm2_0080.h0' };


for ifile = 1:length(bases)

   fname  = sprintf('%s/%s.%d.nc',directory,bases{ifile},year);
   times  = double(nc_read_time(fname,'time'));
   ntimes = length(times);
   dates  = datestr(times);

   fprintf('\n%s.%d.nc has %d times\n',bases{ifile}, year, ntimes)
   fprintf('   first date is %s\n', dates(     1,:))
   fprintf('   last  date is %s\n', dates(ntimes,:))

   gapcheck(bases{ifile},times)

end


function gapcheck(filetype,filetimes)
% predicting how many timesteps there should be ... and any gaps
% ha2x3h    is every 3 hours ...
% ha2x1hi   is wonky ... every hour, some 75 mins, some 45 mins
% ha2x1d    is every 6 hours ... 3,9,15,21 
% ha2x1h    is every hour ... on the :30
% hr2x      is every 6 hours ... on the :00
% mosart.h0 is the first of every month
% clm.h1    is every 6 hours ... on the :00
% clm.h0    is every 6 hours ... on the :00
%
% All the time units are 'days since ...', so the units of delta are easy.

[filepath,~,ext] = fileparts(filetype);
n = length(filetimes);

switch ext
   case '.ha2x3h'           % every 3 hours on the half hour
      delta = 3.0/24.0;
   case {'.ha2x1d','.hr2x'} % every 6 hours on the hour
      delta = 6.0/24.0;
   case '.ha2x1h'           % every hour on the half-hour
      delta = 1.0/24.0;
   case '.ha2x1hi'          % ignore ... wonky
      delta = 1.0/24.0;
      disp('   this file is a conundrum - unclear if it is needed.')
   otherwise
end

% both mosart and lnd have an h0 file, mosart is monthly,
% which has no good delta ... so we are just working with
% the land history files.

switch filepath
   case 'lnd/hist/0080'     % every 6 hours on the hour
      delta = 6.0/24.0;
   otherwise
end

if (exist('delta','var') == 1)
   temptime = [filetimes(1):delta:filetimes(n)];
   expected_length = length(temptime);
   actual_length = n;
   if (actual_length ~= expected_length) 
      fprintf('   WARNING: possible gap.\n')
      fprintf('   WARNING: expected ntimes = %d\n',expected_length)
      fprintf('   WARNING: actual   ntimes = %d\n',  actual_length)
   end
end

