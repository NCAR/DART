%% cam_correl_plot
% usage: assumes cam_correl has already been run and has set a slew of global
% variables, and then cam_correl_read_and_plot has also been run to compute 
% other vars.  as long as you don't change the data files, you can replot 
% with this routine.  if you change the input files, you must rerun 
% cam_correl_read_and_plot.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Construct a colormap with white for points near 0 and colors only when
% the differences become 'significant' (use range to set that width)
% The width of the white area in the middle is set above.
extremes = jet;  % default red to blue rainbow map
zero_pt  = 32;   % 64 bins in the full map
extremes(zero_pt-colorbar_blank:zero_pt+colorbar_blank,:) = 1;   % white

% compute X and Y ticks
xspacing = (plot_lon_max - plot_lon_min) / (plot_lon_ntics - 1);
yspacing = (plot_lat_max - plot_lat_min) / (plot_lat_ntics - 1);
xticks = [plot_lon_min:xspacing:plot_lon_max];
yticks = [plot_lat_min:yspacing:plot_lat_max];

% parse out the print format for the file extension.
file_ext = print_format(3:5);   % -dXXXaaaa   pull out the XXX chars

% what pressure corresponds to level
base_level_pressure  = lev(base_level);
probe_level_pressure = lev(probe_level);


% which direction are we going?
if (direction > 0)     % forward
  probe_lag  = 1;
  probe_time = 1;
  probe_file = 1;
else                   % adjoint
  probe_lag  = num_lags;
  probe_time = steps_per_file;
  probe_file = nfiles;
end


% suggestion for version 2: save the probe location and var type from
% the previous iteration, and if it hasn't changed, don't reread the probe
% file.  this will speed things up a lot and still plot accurate numbers.
% or, plan B - read in entire probe field and then simply squeeze out
% the actual lat,lon values here - which is what the original script did.

% get probe data:
%   use field with the right time - first or last depending on direction
filename = sprintf(file_base_pattern, file_ids(probe_file));

% get probe field (1-d list, single point per ensemble member)
probe_list = Get_Point(filename, probe_var, probe_time, ...
                       probe_pt_lat, probe_pt_lon, probe_level, ens_size);


% get the combined times into valid date formats
tbase = datenum(1601,1,1);
timearray = t+tbase; 
   
for i = 1:length(t)
   disp(sprintf('timestep %d is %s',i,datestr(timearray(i),0)))
end
   
% suggestion for version 2:  right now num_lags has to be an even multiple
% of the number of files times the number of timesteps per file.  this could
% be relaxed by reading in all the timesteps in each file, but then in the
% following loop, for adjoint, start at N to num_lags, where N is the total 
% number of lags read in minus the number to plot, and for forwards, start 
% at 1 but end at num_lags.

% compute the cross correlation coefficients
for lag = 1:num_lags

   for i = 1:lat_res
      for j = 1:lon_res
         fcorr = corrcoef(probe_list, x(lag, :, i, j));
         corr(i, j) = fcorr(1, 2);
      end
   end


   % Do the plotting
   
   % Set regional area
   loninds = find((lon >= plot_lon_min) & (lon <= plot_lon_max));
   latinds = find((lat >= plot_lat_min) & (lat <= plot_lat_max));

   rcorr = corr(latinds, loninds);
   
   figure(lag);                         % plot each lag in its own window
   imagesc(lon(loninds),lat(latinds), rcorr,[-1 1]); % always plots upside-down
   set(gca,'YDir','normal');            % rectify that
   set(gca,'XTick',xticks);             % custom ticklabels & gridlines
   set(gca,'YTick',yticks);
   axis image;
   grid;
   world_land;

   text(lon(probe_pt_lon), lat(probe_pt_lat), '+', 'fontsize', 36);

   colormap(extremes);

   % unless told not to, print axes labels, a colorbar, and a title
   if (nolabels == 0)
     xlabel('degrees east' ,'FontSize',16) % enlarge x labelling
     ylabel('degrees north','FontSize',16) % enlarge y labelling

     h  = colorbar;
     h1 = get(h,'YLabel');
     set(h1,'String',base_varunits,'FontSize',18);
   end
   if (notitle == 0)
     titlestring(1) = { sprintf('%s: %s at %d hPa, at %4.1f E, %4.1f N', ...
                           model, probe_var, floor(probe_level_pressure), ...
                           lon(probe_pt_lon), lat(probe_pt_lat)) };
     titlestring(2) = { sprintf('vs. %s at %d hPa, %s', ...
                              base_var, floor(base_level_pressure), ...
                              datestr(timearray(lag),0)) };
     h = title(titlestring);
     set(h,'FontSize',18,'Interpreter','none') % enlarge title, do not use TeX
   end 

   set(gca,'FontSize',18)                % enlarge the rest of the text.
   
   if (print_to_file)
     if (direction > 0)
       pname = sprintf('%s%d.%s', 'test_forward', lag, file_ext);
     else
       pname = sprintf('%s%d.%s', 'test_adjoint', lag, file_ext);
     end
     print(print_format, pname);
   end

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
