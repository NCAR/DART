%% map_wrf_diff_time

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

field_name = input('Input field type (U, V, W, PH, T, MU, QV, QC, QR, XLAND, VECT, HDIV): ');

map_proj   = {'lambert', 'ups', 'mercator'};
state_name = {'True_State' 'Prior_Diag' 'Posterior_Diag'};

disp('Plotting horizontal map of a linear combination of'); disp (state_name)

fact = input('Input coefficients (between [ ]):');

if fact(1) ~= 0.0
   fname = char(state_name(1));
elseif fact(2) ~= 0.0
   fname = char(state_name(2));
elseif fact(3) ~= 0.0
   fname = char(state_name(3));
else
   error('Nothing to plot.')
end

if (exist(fname,'file') ~= 2)
   error('%s does not exist.',fname)
end

truelat1  = nc_attget(fname, nc_global, 'TRUELAT1');
truelat2  = nc_attget(fname, nc_global, 'TRUELAT2');
cen_lat   = nc_attget(fname, nc_global, 'CEN_LAT');
cen_lon   = nc_attget(fname, nc_global, 'CEN_LON');
stand_lon = nc_attget(fname, nc_global, 'STAND_LON');
mp        = nc_attget(fname, nc_global, 'MAP_PROJ');
dx        = nc_attget(fname, nc_global, 'DX'      );

num_domains = size(dx,1);

if (num_domains > 1)

   disp(['Number of domains: ',int2str(num_domains)])
   id = input('Input domain id: ');

else

   id = 1;

end

xlon  = nc_varget(fname, [ 'XLON_d0',int2str(id)]); we = size( xlon, 2);
xlat  = nc_varget(fname, [ 'XLAT_d0',int2str(id)]); sn = size( xlat, 1);
level = nc_varget(fname, ['level_d0',int2str(id)]); bt = size(level, 1);

cop = zeros(1,3);
cop(1) = -1;

if (fact(2) ~= 0.0) || (fact(3) ~= 0.0)
   if fact(2) ~= 0.0
      fname = char(state_name(2));
   else
      fname = char(state_name(3));
   end
   ncopy = size(nc_varget(fname, 'copy'), 1);
   copy_name = deblank(nc_varget(fname, 'CopyMetaData'));
   for icopy = 1:ncopy
      disp(['Copy # ' int2str(icopy) ' is ' copy_name(icopy,:)])
   end
   cop(2:3) = input('Input copy for Prior and/or Posterior: ');
   string2 = sprintf('%s',copy_name(cop(2),:));
else
   string2 = sprintf('True state');
end

%% FIXME ... this is horrid! use vector _something_ ...

for iy = 1:sn
for ix = 1:we
   if(xlon(iy,ix) > 0.0)
         xlon(iy,ix) = xlon(iy,ix) - 360.0;
   end
end
end

minlat = min(xlat(:)); maxlat = max(xlat(:));
minlon = min(xlon(:)); maxlon = max(xlon(:));

true_times     = nc_varget(fname, 'time');
num_true_times = size(true_times, 1);

stime = input('Initial time : ');
ftime = input('End time : ');

%% Get level for free atmosphere fields

if strcmp(field_name,'MU') || strcmp(field_name,'XLAND')
   field_level = 1;
   vert_coord = 1;
   lev_units = '';
else
   vert_coord = 0;
   while vert_coord ~= 1 && vert_coord ~= 2 && vert_coord ~= 3
      vert_coord = input('Input vertical coordinate: model(1), pressure(2), height(3): ');
      if vert_coord == 1
         field_level = input('Input model level: ');
         lev_units   = 'model level';
      elseif vert_coord == 2
         field_level = input('Input pressure level (hPa): ');
         lev_units   = 'hPa';
      elseif vert_coord == 3
         field_level = input('Input height (meters): ');
         lev_units   = 'm';
      else
         disp('Invalid choice. Please try again')
      end
   end
end

var_name = [field_name,'_d0',int2str(id)];
uname = ['U_d0',int2str(id)];
vname = ['V_d0',int2str(id)];

if strcmp(field_name,'U') || strcmp(field_name,'V') || strcmp(field_name,'VECT')
   var_units = 'm/s';
   iso = 0.5:1:5;
end
if strcmp(field_name,'W')
   var_units = 'm/s';
   iso = 0.01:0.01:0.1;
end
if strcmp(field_name,'PH')
   var_units = 'm^2/s^2';
   iso = 50:50:300;
end
if strcmp(field_name,'T')
   var_units = 'K';
   iso = 250:5:310;
end
if strcmp(field_name,'MU')
   var_units = 'Pa';
   iso = 100:100:600;
end
if strcmp(field_name,'QV')
   var_units = 'kg/kg';
   iso = 0.0001:0.0001:0.001;
end
if strcmp(field_name,'QC')
   var_units = 'kg/kg';
   iso = 0.00001:0.00001:0.0001;
end
if strcmp(field_name,'QR')
   var_units = 'kg/kg';
   iso = 0.00001:0.00001:0.0001;
end
if strcmp(field_name,'XLAND')
   var_units = '-';
end
if strcmp(field_name,'HDIV')
   var_units = '(m/s)/km';
   iso = -0.5:0.05:0.5;
end
if strcmp(field_name,'REF')
   var_units = 'dBZ';
   iso = -10:5:80;
end

m = ceil(sqrt(ftime-stime+1));

%iso = -3:0.25:3;
%iso = 0.1:0.2:5.25;
%iso = 0.01:0.02:5.25;

pane = 1;

for itime = stime:ftime

   plot_title = sprintf('%s (%s)  %i (%s)  day = %i  sec = %i  %s', ...
         field_name, var_units, field_level, lev_units, fix(true_times(itime)), ...
         fix((true_times(itime)-fix(true_times(itime)))*86400), string2);

   %% Extract field

   field  = zeros(sn,we);
   fieldu = zeros(sn,we);
   fieldv = zeros(sn,we);

   for istate = 1:3

      if fact(istate) ~= 0.0

	 fname = char(state_name(istate));

         %% Set up, compute pressure 
         [ Cp, Rd, gamma, Rv, L_c, g, T0, p0 ] = get_constants ;

         [ mu, dnw, phi, theta, qv ] =  ...
            get_aux_fields_for_p( fname, T0, itime, cop(istate), id ) ;

         pres = compute_pressure( mu, dnw, phi, theta, qv, Rd,Rv,gamma,p0 ) ;

         var_in = zeros(bt,sn,we);

         %--Retrieve specified variable from netcdf file

         if vert_coord == 1

	    if strcmp(field_name,'MU')
	       start = [itime cop(istate)  1  1] - 1;
               count = [    1           1 -1 -1];
            elseif strcmp(field_name,'XLAND')
               start = [ 1  1];
               count = [-1 -1];
            else
               start = [itime cop(istate) field_level  1  1] - 1;
               count = [    1           1           1 -1 -1];
            end

	    if strcmp(field_name,'VECT')

               stag_field = nc_varget(fname,uname,start,count);
               for iy = 1:sn
                  for ix = 1:we
	             fieldu(iy,ix) = fieldu(iy,ix) + ...
	                             fact(istate)*(stag_field(iy,ix) + ...
                                     stag_field(iy,ix+1))/2.0;
                  end
               end
               stag_field = nc_varget(fname,vname,start,count);
               for iy = 1:sn
                  for ix = 1:we
                     fieldv(iy,ix) = fieldv(iy,ix) + ...
	       fact(istate)*(stag_field(iy,ix) + stag_field(iy+1,ix))/2.0;
                  end
               end

	    elseif strcmp(field_name,'HDIV')

               stag_field = nc_varget(fname,uname,start,count);
               for iy = 1:sn
                  for ix = 1:we
                     field(iy,ix) = field(iy,ix) - ...
	       fact(istate)*(stag_field(iy,ix+1) - stag_field(iy,ix))/dx(id);
                  end
               end
               stag_field = nc_varget(fname,vname,start,count);
               for iy = 1:sn
                  for ix = 1:we
                     field(iy,ix) = field(iy,ix) - ...
	       fact(istate)*(stag_field(iy+1,ix) - stag_field(iy,ix))/dx(id);
                  end
               end
	       field = field*1000.0;

            else

               stag_field = nc_varget(fname, var_name,start,count);
	       if strcmp(field_name,'U')
                  for iy = 1:sn
                     for ix = 1:we
                        field(iy,ix) = field(iy,ix) + ...
		  fact(istate)*(stag_field(iy,ix) + stag_field(iy,ix+1))/2.0;
                     end
                  end
	       elseif strcmp(field_name,'V')
                  for iy = 1:sn
                     for ix = 1:we
                        field(iy,ix) = field(iy,ix) + ...
		  fact(istate)*(stag_field(iy,ix) + stag_field(iy+1,ix))/2.0;
                     end
                  end
               else
                  field = field + fact(istate)*stag_field;
               end

            end

         elseif vert_coord == 2

            start = [itime cop(istate)  1  1  1] -1;
            count = [    1           1 -1 -1 -1];

	    if strcmp(field_name,'VECT')
               var_inu = zeros(bt,sn,we);
               var_inv = zeros(bt,sn,we);

               stag_field = nc_varget(fname, uname, start, count);
               for iz = 1:bt
	          for iy = 1:sn
	             for ix = 1:we
                        var_inu(iz,iy,ix) = (stag_field(iz,iy,ix) + stag_field(iz,iy,ix+1))/2.0;
                     end
                  end
               end

	       stag_field = nc_varget(fname, vname, start, count);
               for iz = 1:bt
		  for iy = 1:sn
                     for ix = 1:we
                        var_inv(iz,iy,ix) = (stag_field(iz,iy,ix) + stag_field(iz,iy+1,ix))/2.0;
                     end
                  end
	       end

	    elseif strcmp(field_name,'HDIV')

	       stag_field = nc_varget(fname, uname, start, count);
               for iz = 1:bt
		  for iy = 1:sn
	             for ix = 1:we
                        var_in(iz,iy,ix) = var_in(iz,iy,ix) - ...
		  (stag_field(iz,iy,ix+1) - stag_field(iz,iy,ix))/dx(id);
                     end
                  end
               end

	       stag_field = nc_varget(fname, vname, start, count);
               for iz = 1:bt
		  for iy = 1:sn
                     for ix = 1:we
		        var_in(iz,iy,ix) = var_in(iz,iy,ix) - ...
		  (stag_field(iz,iy+1,ix) - stag_field(iz,iy,ix))/dx(id);
                     end
	          end
               end

	       var_in = var_in*1000.0;
            else

	       if strcmp(field_name,'PH')         % Use what get_aux_fields_for_p produced
                  stag_field = phi ;
               elseif strcmp(field_name,'T')     % Same here
                  stag_field = theta ;
               else
                  stag_field = nc_varget(fname, var_name, start, count);
               end

	       if strcmp(field_name,'U')
	          for iz = 1:bt
		     for iy = 1:sn
	                for ix = 1:we
		           var_in(iz,iy,ix) = (stag_field(iz,iy,ix) + stag_field(iz,iy,ix+1))/2.0;
                        end
                     end
                  end
               elseif strcmp(field_name,'V')
	          for iz = 1:bt
	             for iy = 1:sn
	                for ix = 1:we
		           var_in(iz,iy,ix) = (stag_field(iz,iy,ix) + stag_field(iz,iy+1,ix))/2.0;
                        end
                     end
                  end
	       elseif strcmp(field_name,'W') || strcmp(field_name,'PH')
	          for iz = 1:bt
	             for iy = 1:sn
	                for ix = 1:we
		           var_in(iz,iy,ix) = (stag_field(iz,iy,ix) + stag_field(iz+1,iy,ix))/2.0;
                        end
                     end
                  end
               else
                  var_in = stag_field;
               end

            end

	    if strcmp(field_name,'VECT')

	       var_p = interp_to_pressure( var_inu, pres, field_level*100.0) ;
               fieldu = fieldu + fact(istate)*var_p;
               var_p = interp_to_pressure( var_inv, pres, field_level*100.0) ;
               fieldv = fieldv + fact(istate)*var_p;

            else

               var_p = interp_to_pressure( var_in, pres, field_level*100.0) ;
               field = field + fact(istate)*var_p;

            end

         else

            height = compute_height( phi, g ) ;

	        if strcmp(field_name,'REF')

               [ qr, qg, qs ] = get_aux_fields_for_ref( fname, itime, cop(istate), id ) ;

               rho    = compute_density( mu, dnw, phi ) ;
               temp   = compute_temperature( pres, theta, Cp, Rd, p0 ) ;
               var_in = compute_reflectivity( qr, qg, qs, rho, temp ) ;

            end

            var_lev = interp_to_height( var_in, height, field_level) ;
            field = field + fact(istate)*var_lev;

         end

      end

   end

%% Plot field

   subplot(m,m,pane);

   axesm(map_proj{mp(id)},'Origin',[0 cen_lon(id) 0],'MapParallels', ...
	 [truelat1(id) truelat2(id)],...
	 'MapLatLimit',[minlat maxlat],'MapLonLimit',[minlon maxlon]);
   framem;

   plotm(coast,'color',[0 0 0]);
   plotm(usalo('statebvec'),'color',[0 0 0]);
   plotm(usalo('conusvec'),'color',[0 0 0]);

   [xlim ylim]=mfwdtran([xlat(1,1) xlat(sn,we)],[xlon(1,1) xlon(sn,we)]);
   set(gca,'xlim',[min(xlim(:)) max(xlim(:))]);
   set(gca,'ylim',[min(ylim(:)) max(ylim(:))]);

   if strcmp(field_name,'VECT')

      scale = ceil(we/20);

      % legend

      quiverm(xlat(1:scale:sn,1:scale:we),  xlon(1:scale:sn,1:scale:we), ...
            fieldv(1:scale:sn,1:scale:we),fieldu(1:scale:sn,1:scale:we),'k')

   else

      if min(min(field)) ~= max(max(field))

%         field = field/max(max(field));

%         if (max(max(field)) > min(iso)) & (min(min(field)) < max(iso))
%            [C h] = contourm(xlat,xlon,field, iso, 'r','LineWidth',2);
%            h = clabelm(C,h,'labelspacing',288);  set(h,'Fontsize',12);
%            hold on
%         end
%         if (max(max(field)) > min(-iso)) & (min(min(field)) < max(-iso))
%            [Cm hm] = contourm(xlat,xlon,field, -iso, 'b--','LineWidth',2);
%            hm = clabelm(Cm,hm,'labelspacing',288);  set(hm,'Fontsize',12);
%         end

         if strcmp(field_name,'REF')
            eval('dbz_colors')
            h = pcolorm(xlat,xlon,field);
         else 
            [C h] = contourfm(xlat,xlon,field, iso); caxis([min(iso(:)),max(iso(:))]);
         end

         %[C,h] = contour (field, iso);
         %hold on
         %[Cm,hm] = contour (field, -iso, '--');
         cb = colorbar('vert'); set(cb,'Fontsize',12);
         %clabel(C, h);
         %clabel(Cm, hm);

      end

   end

   title(plot_title,'Fontsize',12)

   wysiwyg

   pane = pane + 1;

end

% Loop for another try
%map_wrf_diff_time;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
