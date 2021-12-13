
% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% This is an example matlab script that checks met forcing variable
% characteristics for CAM reanalysis.  This checks
% for missing time steps or suspicious variance across members 

dirname = '/glade/collections/rda/data/ds199.1';
dirname = '/glade/p/cisl/dares/thoar/CAM_DATM/4xdaily';

% 	float a2x6h_Faxa_swndf(time, a2x6h_ny, a2x6h_nx) ;

if (exist('variance','var') ~= 1)

   ntimes = 1460;
   a2x6h_ny = 96;
   a2x6h_nx = 144;
   members  = 80;
   tensor = zeros(a2x6h_nx,a2x6h_ny,ntimes,members,'single');

% for iyear = [1999 2008 2009 2010]
   for iyear = 2008:2008
   for imem = 1:members
   
      filename = sprintf('%s/CAM_DATM.cpl_%04d.ha2x1dx6h.%d.nc', ...
                         dirname,imem,iyear);

   if ( exist(filename,'file') ~= 2)
      error('%s does not exist',filename)
   end

      times       = ncread(filename,'time');
   %  timeunits   = ncread(filename,'time','units');
      ncid        = netcdf.open(filename);
      varid       = netcdf.inqVarID(ncid,'time');
      timeunits   = netcdf.getAtt(ncid,varid,'units');
      netcdf.close(ncid)
   timebase    = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
   timeorigin  = datenum(timebase(1),timebase(2),timebase(3));
   ttimes      = times + timeorigin;

   ntimes = length(ttimes);
   deltat = diff(ttimes);
   delta  = median(deltat);
   inds   = find(deltat > delta);

   if ((length(ttimes) == 365*4) && isempty(inds)) 
      fprintf('%s has all the timesteps.\n',filename)
   else
      fprintf('WARNING: %s has %d timesteps, not %d.\n',filename,ntimes,365*4)
      disp('WARNING: problem occurrs around')
      datestr(ttimes(inds))
   end

   %   a2x6h_Faxa_swndr = ncread(filename,'a2x6h_Faxa_swndr');
   %   mymin = min(a2x6h_Faxa_swndr(:));
   %   mymax = max(a2x6h_Faxa_swndr(:));
   %   fprintf(1,'%04d %04d swndr min %f max %f \n',iyear,imem,mymin,mymax)
   %
   %   a2x6h_Faxa_swvdr = ncread(filename,'a2x6h_Faxa_swvdr');
   %   mymin = min(a2x6h_Faxa_swvdr(:));
   %   mymax = max(a2x6h_Faxa_swvdr(:));
   %   fprintf(1,'%04d %04d swvdr min %f max %f \n',iyear,imem,mymin,mymax)
   %
   %   a2x6h_Faxa_swndf = ncread(filename,'a2x6h_Faxa_swndf');
   %   mymin = min(a2x6h_Faxa_swndf(:));
   %   mymax = max(a2x6h_Faxa_swndf(:));
   %   fprintf(1,'%04d %04d swndf min %f max %f \n',iyear,imem,mymin,mymax)
   
      a2x6h_Faxa_swvdf = ncread(filename,'a2x6h_Faxa_swvdf');
      mymin = min(a2x6h_Faxa_swvdf(:));
      mymax = max(a2x6h_Faxa_swvdf(:));
      fprintf(1,'%04d %04d swvdf min %f max %f \n',iyear,imem,mymin,mymax)
   
      tensor(:,:,:,imem) = a2x6h_Faxa_swvdf;
      clear a2x6h_Faxa_swvdf
   
   end
   end
end

t_interest = [352, 880:1000];
t_interest = [352];

cmax = zeros(size(t_interest));

mycount = 1;
%for itime=1:1
for itime=t_interest

   datmat    = squeeze(tensor(:,:,itime,:));  %   datmat(lon,lat,member)
   variance  = var(datmat,0,3,'omitnan');     % variance(lon,lat)
   quantiles = quantile(datmat, [0.25 0.5 0.75], 3);
   maxvar    = max(variance(:));

   [i,j] = find(variance == maxvar);

   subplot(4,1,1)
      imagesc(variance')
      set(gca,'YDir','normal')
      colorbar
      title(sprintf('variance %d of %d',itime,ntimes))

   subplot(4,1,2)
      iqrange   = quantiles(:,:,3) - quantiles(:,:,1);
      imagesc(iqrange')
      set(gca,'YDir','normal')
      colorbar
      title('interquartile range')

   subplot(4,1,3)
      datmatmax = max(datmat,[],3,'omitnan') - quantiles(:,:,3);
      k = datmatmax ./ iqrange;
      imagesc(k')
      set(gca,'YDir','normal')
      colorbar
      title('tukey''s fence value')

   % remove the maximum value at each location and recompute the mean, 
   % divide the max value by the new mean
    subplot(4,1,4)
       datmatmax    = max(datmat,[],3,'omitnan');       % datmatmax(nlon,nlat)
       datmatmedian = median(datmat,3,'omitnan');
       k = datmatmax ./ datmatmedian;
       imagesc(k')
       set(gca,'YDir','normal')
       colorbar
       title('max/median')

      cmax(mycount) = max(k(:));
      [ii,jj] = find(k == cmax(mycount));
      fprintf(1,'timestep %5d max variance %12.3f %3d %3d max/median %12.3f %3d %3d\n', ...
                        itime,       maxvar,  i,  j,   cmax(mycount), ii, jj)

      mycount = mycount + 1;

%     suspicious = squeeze(datmat(37,61,:))
      suspicious = squeeze(datmat(ii,jj,:))

    pause(0.1)
   % pause()
end

% timestep 352 max variance 25982343.392381
% explore lots of isolated high variance pixels

if ( 0 )

   figure(2)

   t_inds = [352 880 887 891 895 899 908 911 919 924 927 936 940 947 955 956 959 963 967 980 984 987 992];
   i_inds = [144   1  37  37  37  37   5  37  37   4  37   5   4  37  36   4  37  37  37   4   4  37   4];
   j_inds = [ 76  67  61  61  61  61  74  61  61  74  61  75  76  61  61  76  61  61  61  76  74  61   4];
   
   for i = 1:length(i_inds)
      members = squeeze(tensor(i_inds(i),j_inds(i),t_inds(i),:));
      datmat(:)
      plot(datmat,'*')
   end

end

%%
%% max swvdf is 3x bigger than the rest
% /glade/p_old/image/thoar/CAM_DATM/4xdaily/CAM_DATM.cpl_0016.ha2x1dx6h.2008.nc has all the timesteps.
% 2008 0016 swvdf min 0.000000 max 1260.572998 
%  time = 1460;
%  a2x6h_ny = 96;
%  a2x6h_nx = 144;
%  members  = 80;

% ncap2 -v a2x6h_Faxa_swvdf -d a2x6h_nx,143 -d a2x6h_ny,75 -d time,351 CAM_DATM.cpl_0064.ha2x1dx6h.2008.nc
