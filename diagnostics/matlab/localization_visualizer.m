function localization_visualizer
% localization_visualizer helps explore the localization impact for a given cutoff value.
%           You can change projections, resolutions, location, etc.
%           Have fun with it!
%
% The m_map toolbox from Rich Pawlowicz at UBC is required.
% m_map is freely available from http://www.eos.ubc.ca/~rich/map.html
% The script will throw an error if m_map is not found in your Matlab path.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% TJH thoughts about possible extensions:
% * change colorbar ticks to match the levels of the contour plots.
% * give option to read the lat/lon arrays (and projection?) from a
%   netCDF file and plot the grid on the graphic.
% * pick the data point with the cursor.
% * allow user to directly input a known cutoff value in addition to
%   using the slider.
% * outside the zero contour should be white, not cyan (or color 1).
% * should support some cylindrical projection

if  exist('m_map', 'dir') == 0
    error('m_map:notExist', ...
         ['The m_map toolbox is required and is freely available \n' ...
          'from http://www.eos.ubc.ca/~rich/map.html']);
end

a = ver('matlab');
if sscanf(a.Version,'%f') >= 8.4
    error('matlab:graphics', ...
         ['The new graphics requirements in 2014b are not fully \n' ...
          'supported for this yet. Must use an earlier version.']);
end

% Longitude and Latitude of Boulder, Colorado, the default location.
location(1,1) = -105.27; % x is longitude
location(1,2) = 40.015;  % y is latitude

% calculate distance between each x,y gridpoint and the target location.
covarianceGui(location);

end

%----------------------------------------------------------------------
% End of the main function.
%----------------------------------------------------------------------
% Start of the GUI and its functions. There are some scoping
%    requirements for the covarianceGui that require the use of
%    internal functions.
%----------------------------------------------------------------------

function covarianceGui(data)
% GUI which shows the weighing factor for an onservation as a contour plot.
% The m_map toolbox from Rich Pawlowicz at UBC is required to run the
% function. m_map is freely available from http://www.eos.ubc.ca/~rich/map.html

% How the contour plot is created:
%  calculate distance between all the lon lat grid points and the data point
%  calculate the comp_cov_factor for each point
%  plot the contour of the 0 value for comp_cov_factor = 0

% keep initial data values
init_data = data;

% Initial values for cutoff contour (in radians)
minCutOff_rad = 0.005;
maxCutOff_rad = 0.7;
   cutOff_rad = maxCutOff_rad/2.0;

% kilometers
minCutOff_km = rad2km(minCutOff_rad);
maxCutOff_km = rad2km(maxCutOff_rad);
   cutOff_km = rad2km(   cutOff_rad);

% Resolution of contour plot
% TJH should same value be used for both lon and lat
res = 360;

% Inital values for map projection
Long = data(1);
Lat  = data(2);
Rad  = 90;

% Inital projection type
projection = 'Stereographic';
level = 0;   % plot the contour of the 0 value for comp_cov_factor = 0

%x=linspace(-180,180,res); % Longitude
%y=linspace( -90, 90,res); % Latitude

x = linspace(  0,360,res); % Longitude
y = linspace(-90, 90,res); % Latitude

[X,Y] = meshgrid(x,y); % should remesh finer when user zooms in

% calculate distance between all the grid points and each data point
% calculate the comp_cov_factor for each point
distance = calc_dist(data, X, Y);
type = 'Gaspari-Cohn';
Z = comp_cov_factor(distance, cutOff_km, type);

%% Create and then TJH 'do not' hide the GUI as it is being constructed
f = figure('Visible','on', ...
           'Units','centimeters', ...
           'Position',[8,8,28,20], ...
           'Color',[0.8 0.8 0.8]);

%% Construct the components

% data point
dataTitle  = uicontrol('Style', 'text', ...
                       'String', 'data point', ...
                       'Units', 'centimeters', ...
                       'Position', [18, 15.8, 3, 1], ...
                       'FontWeight', 'bold');
dataValue1 = uicontrol('Style', 'edit', ...
                       'String', num2str(data(1)), ...
                       'Units', 'centimeters', ...
                       'Position', [21, 16, 2.7, 1], ...
                       'Callback', @move_data1);
dataValue2 = uicontrol('Style', 'edit', ...
                       'String', num2str(data(2)), ...
                       'Units', 'centimeters', ...
                       'Position', [23.7, 16, 2.7, 1], ...
                       'Callback', @move_data2);

% cutoff
hcutOffTitle    = uicontrol('Style', 'text', ...
                            'String', 'Cut Off (halfwidth)', ...
                            'Units', 'centimeters', ...
                            'Position', [18, 14, 8.4, 1], ...
                            'FontWeight', 'bold');
hcutOffRadTitle = uicontrol('Style', 'text', ...
                            'String', 'radians', ...
                            'Units', 'centimeters', ...
                            'Position', [18, 12, 3, 1]);
hcutOffKmTitle  = uicontrol('Style', 'text', ...
                            'String', 'kilometers', ...
                            'Units', 'centimeters', ...
                            'Position', [18, 11, 3, 1]);
hcutOffRad      = uicontrol('Style', 'text', ...
                            'String', num2str(cutOff_rad), ...
                            'Units', 'centimeters', ...
                            'Position', [21, 12, 5.4, 1]);
hcutOffKm       = uicontrol('Style', 'text', ...
                            'String', num2str(cutOff_km), ...
                            'Units', 'centimeters', ...
                            'Position', [21, 11, 5.4, 1]);

hCutOffType     = uicontrol('Style', 'popupmenu', ...
                            'String', {'Gaspari-Cohn', 'Boxcar', 'Ramped Boxcar'}, ...
                            'Units', 'centimeters', ...
                            'Position', [18, 13, 8.4, 1], ...
                            'Callback', @cut_pop);
hslider         = uicontrol('Style', 'slider', ...
                            'Units', 'centimeters', ...
                            'Position', [18, 10, 8.4, 1], ...
                            'Min', minCutOff_rad, ...
                            'Max', maxCutOff_rad, ...
                            'Value', cutOff_rad, ...
                            'Callback', {@display_slider_value});

% Map
ha = axes('Units', 'centimeters', ...
          'Position',[2,3,15,15]);

hTitleCentre = uicontrol('Style', 'text', ...
                         'String', 'Map Center', ...
                         'Units', 'centimeters', ...
                         'Position', [18, 6, 8.4, 1], ...
                         'FontWeight', 'bold');
hLong        = uicontrol('Style', 'edit', ...
                         'String', num2str(Long), ...
                         'Units', 'centimeters', ...
                         'Position', [18, 5, 2.8, 1], ...
                         'Callback', @longitude);
hTitleLong   = uicontrol('Style', 'text', ...
                         'String', 'Longitude', ...
                         'Units', 'centimeters', ...
                         'Position', [18, 4, 2.8, 1]);
hLat         = uicontrol('Style', 'edit', ...
                         'String', num2str(Lat), ...
                         'Units', 'centimeters', ...
                         'Position', [20.8, 5, 2.8, 1], ...
                         'Callback', @latitude);
hTitleLat    = uicontrol('Style', 'text', ...
                         'String', 'Latitude', ...
                         'Units', 'centimeters', ...
                         'Position', [20.8, 4, 2.8, 1]);
hRad         = uicontrol('Style', 'edit', ...
                         'String', num2str(Rad), ...
                         'Units', 'centimeters', ...
                         'Position', [23.6, 5, 2.8, 1], ...
                         'Callback', {@radius});
hTitleRad    = uicontrol('Style', 'text', ...
                         'String', 'Radius', ...
                         'Units', 'centimeters', ...
                         'Position', [23.6, 4, 2.8, 1]);
hpopupTitle  = uicontrol('Style', 'text', ...
                         'String', 'Projection', ...
                         'Units', 'centimeters', ...
                         'Position', [18, 8, 8.4, 1], ...
                         'FontWeight', 'bold');

hpopup       = uicontrol('Style', 'popupmenu', ...
                         'String',{ 'Stereographic', ...
                                    'Orthographic', ...
                                    'Azimuthal Equal-area', ...
                                    'Azimuthal Equidistant', ...
                                    'Satellite', ...
                                    'Equidistant Cylindrical', ...
                                    }, ...
                         'Units','centimeters', ...
                         'Position',[18,7,8.4,1], ...
                         'Callback',{@popup_menu_Callback});

% reset button
hbutton = uicontrol('Style', 'pushbutton', ...
                    'String', 'Reset', ...
                    'Units', 'centimeters', ...
                    'Position', [23 2 3 1], ...
                    'Callback', @start_again);
htext   = uicontrol('Style', 'text', ...
                    'String', 'Localization contours. Zero contour is black.', ...
                    'Units', 'centimeters', ...
                    'Position', [2 1 18 1]);

%% Initialize the GUI
% Change units to normalized so components resize automatically.
set([f, ha, hslider, hpopup, hLong, hLat, hTitleCentre, hTitleLat, ...
     hTitleLong, hRad, hTitleRad, hCutOffType, hpopupTitle, ...
     hcutOffTitle, hcutOffKm, hcutOffRad, hcutOffKmTitle, hcutOffRadTitle, ...
     dataTitle, dataValue1, dataValue2, hbutton, htext],'Units','normalized');

% Create a plot in the axes.
% Assign the GUI a name to appear in the window title.

plot_image(projection, Long, Lat, Rad, data, level, X,Y,Z)
set(f,'Name','Covariance Factor');
% movegui(f,'center'); % Move the GUI to the center of the screen.

% set font size
set([hslider, hpopup, hLong, hLat, hTitleCentre, hTitleLat, hTitleLong, ...
     hRad, hTitleRad, hCutOffType, hpopupTitle, hcutOffTitle, hcutOffKm, ...
     hcutOffRad, hcutOffKmTitle, hcutOffRadTitle, dataTitle, dataValue1, ...
     dataValue2, htext], 'FontSize', 14, 'Backgroundcolor',[0.8 0.8 0.8]);

% Make the GUI visible.
set(f,'Visible','on');

    %------------------------------------------------------------------
    % Start of the functions scoped internal to covarianceGUI
    %------------------------------------------------------------------

    %% convert radians to km assuming the earth-like radius

    function km = rad2km(radians)
        km = radians*6371;
    end


    %% calclate distance between the data point and every grid point

    function d_out = calc_dist(point_in, gridX, gridY)
        % The m_map function m_lldist gives the distance in kilometers
        % between successive points in the input vectors. Thus array
        % is built to be data point, grid point, data point, grid point, etc.

        array = zeros(2*size(gridX(:),1), 2);

        array(:, 1) = point_in(1); % lon
        array(:, 2) = point_in(2); % lat
        array(1:2:end, :) = [gridX(:), gridY(:)]; % lon lat

        all_dists = m_lldist(array(:,1), array(:,2));

        % only want half these values TJH DEBUG ... WHY
        d_array = all_dists(1:2:end);
        d_out = reshape(d_array, res, res);

    end


    %% Grab slider value


    function display_slider_value(source, ~)
        % The cut off is grabbed from the slider value then used to
        % recalculate the contours.
        set(hcutOffRad, 'String', num2str(get(source,'Value'), '%6.4g'));
        set(hcutOffKm, 'String', num2str(rad2km(get(source,'Value')), '%6.4g'));

        cutOff_rad = get(source, 'Value');
        cutOff_km = rad2km(get(source, 'Value'));% recalculate cutoff
        Z = comp_cov_factor(distance, cutOff_km, type);
        plot_image(projection, Long, Lat, Rad, data, level, X,Y,Z);
    end


    %%  Does the plotting


    function plot_image(projection, Long, Lat, Rad, data, level, X,Y,Z)


        cla
        switch projection
            case 'Equidistant Cylindrical'
                m_proj(projection, 'lon', Long, 'lat', Lat);
            otherwise
                m_proj(projection, 'lon', Long, 'lat', Lat, 'rad', Rad);
        end
        hold on

        [~, hf] = m_contourf(X, Y, Z, 5); shading flat;

        % make the last contour translucent
        % TJH This is the part that blows up in 2014b.
        % TJH This is the part that blows up for Equidistant Cylindrical (all versions)
        allH = allchild(hf);
        set(allH(end),'FaceAlpha',0.2);

        % [level, level] to get a single line
        [~, ~] = m_contour(X, Y, Z, [level, level], 'lineColor', 'black');

        m_plot(data(:, 1), data(:, 2), '+k')
        m_coast;
        m_grid;

        % m_map declares  a 'bug' in 2013b, 2014a ...
        % this is the workaround that does not invade m_map code.
        a = ver('matlab');
        if sscanf(a.Version,'%f') > 8.1
            set(gca,'dataaspectratio',[1 1 1e16]);
        end

        colorbar;
        caxis([0 1]);
        colormap(cool);
        hold off;
    end


    %% Map center and extent


    function longitude(source, ~)
        temp = get(source, 'String'); % get the longitude value
        Long = str2double(temp);
        str  = get(hpopup, 'String'); % get the projection in use
        projection = str{get(hpopup,'Value')};
        plot_image(projection, Long, Lat, Rad, data, level, X,Y,Z);
    end

    function latitude(source, ~)
        temp = get(source, 'String'); % get the latitude value
        Lat  = str2double(temp);
        str  = get(hpopup, 'String'); % get the projection in use
        projection = str{get(hpopup,'Value')};
        plot_image(projection, Long, Lat, Rad, data, level, X,Y,Z);
    end

    function radius(source, ~)
        temp = get(source, 'String'); % get the radius value
        Rad  = str2double(temp);
        str  = get(hpopup, 'String'); % get the projection in use
        projection = str{get(hpopup,'Value')};
        plot_image(projection, Long, Lat, Rad, data, level, X,Y,Z);

    end


    %% Pop-up menu for cutoff type


    function cut_pop(source, ~)
        str = get(source, 'String');
        val = get(source, 'Value');
        type = str{val};

        switch str{val}
            case 'Gaspari-Cohn'
                Z = comp_cov_factor(distance, cutOff_km, 'Gaspari-Cohn');
                plot_image(projection, Long, Lat, Rad, data, level, X,Y,Z);
            case 'Boxcar'
                Z = comp_cov_factor(distance, cutOff_km, 'Boxcar');
                plot_image(projection, Long, Lat, Rad, data, level, X,Y,Z);
            case 'Ramped Boxcar'
                Z = comp_cov_factor(distance, cutOff_km, 'Ramped Boxcar');
                plot_image(projection, Long, Lat, Rad, data, level, X,Y,Z);
        end
    end


    %% Enter data point


    function move_data1(source, ~)

        % get the data point
        % move the map center
        % need to recalculate distance
        % recalculate the comp_cov_factor

        data(1)  = str2double(get(source, 'String'));
        Long     = data(1);
        distance = calc_dist(data, X, Y);
        Z        = comp_cov_factor(distance, cutOff_km, type);
        set(hLong, 'String', num2str(data(1), '%6.4g'));
        plot_image(projection, Long, Lat, Rad, data, level, X,Y,Z);
    end

    function move_data2(source, ~)

        % get the data point
        % move the map center
        % need to recalculate distance
        % recalculate the comp_cov_factor

        data(2)  = str2double(get(source, 'String'));
        Lat      = data(2);
        distance = calc_dist(data, X, Y);
        Z        = comp_cov_factor(distance, cutOff_km, type);
        set(hLat, 'String', num2str(data(2), '%6.4g'));
        plot_image(projection, Long, Lat, Rad, data, level, X,Y,Z);
    end


    %%  Pop-up menu callback. Read the pop-up menu Value property


    function popup_menu_Callback(source, ~)
        % Determine the selected projection.
        str = get(source, 'String');
        val = get(source,'Value');

        projection = str{val};

        switch str{val}; % Set current data to the selected data set.
            case 'Stereographic'
                m_proj('Stereographic',           'lon', Long, 'lat', Lat, 'rad', Rad);
            case 'Orthographic'
                m_proj('Orthographic',            'lon', Long, 'lat', Lat, 'rad', Rad);
            case  'Azimuthal Equal-area'
                m_proj('Azimuthal Equal-area',    'lon', Long, 'lat', Lat, 'rad', Rad);
            case 'Azimuthal Equidistant'
                m_proj('Azimuthal Equidistant',   'lon', Long, 'lat', Lat, 'rad', Rad);
            case 'Satellite'
                m_proj('Satellite',               'lon', Long, 'lat', Lat, 'rad', Rad);
            case 'Equidistant Cylindrical'
                m_proj('Equidistant Cylindrical', 'lon', Long, 'lat', Lat);
        end

        plot_image(projection, Long, Lat, Rad, data, level, X,Y,Z);

    end


    %% Reset GUI


    function start_again(~, ~)

        % map
        Long = init_data(1);
        Lat  = init_data(2);
        Rad  = 90;
        set(hLong, 'String', num2str(Long));
        set(hLat,  'String', num2str(Lat));
        set(hRad,  'String', num2str(Rad));
        projection = 'Stereographic';
        set(hpopup, 'Value', 1);

        % cutoff
        cutOff_rad = maxCutOff_rad/2.0;
        cutOff_km  = rad2km(cutOff_rad);
        set(hcutOffRad,  'String', num2str(cutOff_rad));
        set(hcutOffKm,   'String', num2str(cutOff_km));
        set(hCutOffType, 'Value', 1);
        set(hslider,     'Value', cutOff_rad);

        % data
        data = init_data;
        set(dataValue1, 'String', num2str(data(1)));
        set(dataValue2, 'String', num2str(data(2)));


        % need to recalculate distance and
        % calculate the comp_cov_factor for each point
        distance = calc_dist(data, X, Y);
        Z = comp_cov_factor(distance, cutOff_km, 'Gaspari-Cohn');

        plot_image(projection, Long, Lat, Rad, data, level, X,Y,Z)

    end


end

%----------------------------------------------------------------------
% End of the GUI and its functions.
%----------------------------------------------------------------------
% Start of the covariance calculations.
% comp_cov_factor also exists in the DART_LAB/matlab section.
%----------------------------------------------------------------------


function CF = comp_cov_factor(z_in, c, type)
%% comp_cov_factor
%  z_in is the distance
%  c    is the cutoff
%  type is the shape of the cutoff function

Zed = abs(z_in);

if strcmp(type,'Gaspari-Cohn')
    CF = arrayfun(@calc_cov_factor, Zed);
elseif strcmp(type,'Boxcar')
    CF = arrayfun(@boxcar, Zed);
elseif strcmp(type,'Ramped Boxcar')
    CF = arrayfun(@ramped, Zed);
end

CF(CF < 0) = 0; % This removes the negative numbers, but I don't know why they occur

    function cov_factor = calc_cov_factor(x)
        % Can get spurious negative cov_factor sometimes.  See fudge above.
        if( x >= c*2.0)
            cov_factor = 0;
        elseif( x <= c )
            r = x / c;
            cov_factor = ((( -0.25*r +0.5)*r +0.625)*r -5.0/3.0)*r^2 + 1.0;
        else
            r = x / c;
            cov_factor = ((((r/12 -0.5)*r +0.625)*r +5.0/3.0)*r -5.0)*r + 4.0 - 2.0 / (3.0 * r);
        end

    end

    function cov_factor = boxcar(x)
        if (x < c*2.0)
            cov_factor = 1.0;
        else
            cov_factor = 0.0;
        end
    end

    function cov_factor = ramped(x)
        if(x >= 2.0 * c)
            cov_factor = 0.0;
        elseif x >= c && x < 2.0 * c
            cov_factor = (2.0 * c - x) / c;
        else
            cov_factor = 1.0;
        end

    end

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
