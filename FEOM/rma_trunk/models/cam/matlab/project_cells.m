function project_cells(file_in,ob_selection)
%% Displays an ob location and the enclosing cell, as viewed from the center of the earth,
%  so that there's no distortion of mapping the curved cell onto a flat picture.
%  The ob and cell corner locations need to be in the form ({ob,corner_name} lon lat), e.g.
%     45    9.331000000000E+01  0.000000000000E+00
%     12302 9.382917960675E+01 -4.363499758141E-15
%     12299 9.382917960675E+01 -8.273287929204E-01
%     12292 9.300000000000E+01 -8.280434032478E-01
%     12295 9.300000000000E+01 -4.367269282491E-15

%
% ob_selection requests a subset of the obs available in file_in.
%   Note that the ob numbers in file_in are monotonic, but may not be sequential,
%   especially if redundant obs locations (horizontal) have been weeded out.

% I want to plot an individual ob location with its enclosing cell.
% The clearest way to see if the ob is in the cell is to view them
% from the center of the earth.

% borderline obs to investigate:
% ob_selection = [45 360 2700 6975 20340 20475 22185 22860 22950 25155 ]

% The input data file was generated from the original output from 
% $p/Models/SE_RUD_pmo2/cesm.stdout.273578
% which was written by model_mod:coord_ind_cs in loop over num_close, after unit_square_location
%    if (found_quad == 1 ) then
%       ! kdr location check
%       if (lon_lat_lev(3) == 500.) then
%          write(*,'(A,1p2E20.12) ')'{lon,lat} ob  = ', lon_lat_lev(1), lon_lat_lev(2)
%       end if
%    
%       exit
%    end if
% 
% and model_mod:interp_cubed_sphere    after call to coord_ind_cs
%    ! kdr location check
%    if (lon_lat_lev(3) == 500.) then
%       do i=1,4
%          write(*,'(A,I8,A,1p2E20.12) ')'{lon,lat} ',quad_corners(i),' > ', &
%                lon%vals(quad_corners(i)), lat%vals(quad_corners(i))
%       enddo
%    end if
% I grepped for {lon  or something similar to extract the lines I needed.
% That was converted to an array by project_cells_make_array.csh

% 
% I may want perspective projection using a Z-buffer:
%     CameraViewAngle determines extent of scene at CameraTarget.
% 
% h = plotm(lat,lon,PropertyName,PropertyValue,...) 
%     Displays projected line objects on the current map axes. 
%     lat and lon are the latitude and longitude coordinates, 
%     respectively, of the line object to be projected. 
%     Note that this ordering is conceptually reversed from the MATLAB® line function, 
%     because the vertical (y) coordinate comes first. 
%     However, the ordering latitude, then longitude, is standard geographic usage. 
%     lat and lon must be the same size, and in the AngleUnits of the map axes.
%     Name/property value pairs for any properties recognized by the MATLAB line function 
%     except for XData, YData, and ZData.
% 
% camproj('projection_type') 
%     sets the projection type in the current axes to the specified value. 
%     Possible values for projection_type are orthographic and perspective.
%     Comes after the plotting routine.
% 
% campos([camera_position])
%     Sets the position of the camera in the current axes to the specified value. 
%     Specify the position as a three-element vector containing the x-, y-, and z-coordinates 
%     of the desired location in the data units of the axes. 
% 
% camtarget([camera_target])  not used
%     Sets the camera target in the current axes to the specified value. 
%     Specify the target as a three-element vector containing the x-, y-, and z-coordinates 
%     of the desired location in the data units of the axes. 
% 
% camva('auto' or angle_in_degrees)  used, but not needed, since I use default?
%     the camera view angle adjusts so that the scene fills the available space in the window. 
% 
% 'globe' axes
%     In the three-dimensional sense, globe is true in scale, equal-area, conformal, 
%     minimum error, and equidistant everywhere. 
%     When displayed, however, it looks like an Orthographic azimuthal projection, 
%     provided that the MATLAB® axes Projection property is set to 'orthographic'.
%     This is the only three-dimensional representation provided for display. 
%     Unless some other display purpose requires three dimensions, 
%     the Orthographic projection's display is equivalent.
%     
%     axesm ('globe','Grid', 'on');
%     view(60,60)
%     axis off
%     % Display coastline vectors
%     load coast
%     plotm(lat, long)
% 
%     The Globe display is based on a coordinate transformation, and is not a map projection. 
%     While it has none of the distortions inherent in planar projections, it is a 
%     three-dimensional model of a planet that cannot be displayed without distortion or in its entirety. 
%     That is, in order to render the globe in a figure window, either a perspective or orthographic 
%     transformation must be applied, both of which necessarily involve setting a viewpoint, 
%     hiding the back side and distortions of shape, scale, and angles.
% 
%     Over-the-Horizon 3-D Views Using Camera Positioning Functions
% 
%     You can create dramatic 3-D views using the Globe display. 
%     The camtargm and camposm functions (Mapping Toolbox functions corresponding to camtarget and campos) 
%     enable you to position focal point and a viewpoint, respectively, in geographic coordinates, 
%     so you do not need to deal with 3-D Cartesian figure coordinates.
% 
% view
%     The position of the viewer (the viewpoint) determines the orientation of the axes. 
%     You specify the viewpoint in terms of azimuth and elevation, or by a point in three-dimensional space.
% 
%     view(az,el) and view([az,el]) set the viewing angle for a three-dimensional plot. 
%     The azimuth, az, is the horizontal rotation about the z-axis as measured 
%     in degrees from the negative y-axis. 
%     Positive values indicate counterclockwise rotation of the viewpoint. 
%     el is the vertical elevation of the viewpoint in degrees. 
%     Positive values of elevation correspond to moving above the object; 
%     negative values correspond to moving below the object.
% 
%     view([x,y,z]) sets the view direction to the Cartesian coordinates x, y, and z. 
%     The magnitude of (x,y,z) is ignored.
% 
%     legend(string1,string2,...,...
%            'Location',{'choose_string',1x4PositionVector},...
%            'boxoff',...
%     text(x,y,[z,]'string','PropertyName',PropertyValue)
% 
% Tim suggests using the matlab m_map package that he loaded onto
% yellow/caldera; examples of using the mapping toolbox.
% 
data = load(file_in);

% Generate dimensions of arrays
% Each ob has 5 lines and 3 columns associated with it.
data_size = size(data)
col_len = data_size(1);
num_obs = data_size(1)./5
obs_todo = length(ob_selection);

% Gather the ob numbers and locations from the data array.
% Line 1 of each ob section.
obs =     data(1:5:col_len,1);
obs_lon = data(1:5:col_len,2);
obs_lat = data(1:5:col_len,3);

ob_incr = obs(2) - obs(1);

% Gather the corner numbers and locations from the data array.
% Lines 2-5 of each ob section.
% Define 5 corners to close the cell.
corners     = zeros(num_obs,5);
corners_lon = zeros(num_obs,5);
corners_lat = zeros(num_obs,5);
for o = 1:num_obs
   oline1 = (o-1) * 5 + 1;
   for c = 1:4
      corners(o,c)     = data(oline1+c,1);
      corners_lon(o,c) = data(oline1+c,2);
      corners_lat(o,c) = data(oline1+c,3);
   end
   corners(o,5)     = corners(o,1);
   corners_lon(o,5) = corners_lon(o,1);
   corners_lat(o,5) = corners_lat(o,1);
end

% Loop over all selected obs to display each within its associated cell.
label_str = cell(5,1);
ob_count = 1;
for o = 1:num_obs
o_name = o * ob_incr;
% cond1 = ob_count <= obs_todo
% cond2 = o_name == ob_selection(ob_count)
% ? Why won't the first frame allow frame menu/control panel control of picture,
%   but second frame will?
if ((ob_count <= obs_todo) && (o_name == ob_selection(ob_count)) ) 

   figure(1);
%  axesm turns 'hold' to 'on' by default
   ax_h = axesm ('globe')
   % getm (ax_h)
   % get (ax_h)
   axis off

   camproj('perspective')
   campos([0 0 0])
   camva('auto');

% Plot the 4 corners of the cell and lines (= great circles in this view) between them.
   cs_cell = plotm([corners_lat(o,1:5)],[corners_lon(o,1:5)] );

% Resize the window to accomodate coordinates text.
   lat_lims = get(ax_h,'YLim') ;
   lat_factor = abs(lat_lims(2)-lat_lims(1)) ;
   lat_lims(1) = lat_lims(1) - .1*lat_factor ;
   lat_lims(2) = lat_lims(2) + .1*lat_factor ;

   lon_lims = get(ax_h,'XLim') ;
   lon_factor = abs(lon_lims(2)-lon_lims(1)) ;
   lon_lims(1) = lon_lims(1) - .1*lon_factor ;
   lon_lims(2) = lon_lims(2) + .5*lon_factor ;

   z_lims = get(ax_h,'ZLim') ;
   z_factor = abs(z_lims(2)-z_lims(1)) ;
   z_lims(1) = z_lims(1) - .1*z_factor ;
   z_lims(2) = z_lims(2) + .1*z_factor ;

   set(ax_h, 'YLim',[lat_lims],'XLim',[lon_lims],'ZLim',[z_lims]) ;
%    set(ax_h, 'YLim',[lat_lims])
%    set(ax_h, 'XLim',[lon_lims])
%    set(ax_h, 'ZLim',[z_lims])
    
%    lat_new = get(ax_h,'YLim')
%    lon_new = get(ax_h,'XLim')

% Generate and plot text labels for cell corners; lon,lat, and node number
   for c = 1:4
      label_str{c} = sprintf('(%.2f,%.2f,%d)',corners_lon(o,c),corners_lat(o,c),corners(o,c));
      textm(corners_lat(o,c),corners_lon(o,c) ,label_str{c})
   end
% Same for observation, but ob number within obs_seq file replaces node number.
   ob = plotm([obs_lat(o)], [obs_lon(o)], ':xr');
   label_str{5} = sprintf('(%.2f,%.2f,%d)',obs_lon(o),obs_lat(o),obs(o));
   textm(obs_lat(o),obs_lon(o) ,label_str{5},'Color','red')

   title(sprintf('ob %d',o_name))

%    hold off
   disp('pausing - hit any key to continue ...'); pause
   close 
   ob_count = ob_count + 1;
end
end
