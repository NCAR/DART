function plot_connections(x, ROUTE_FILE, fig_pos, show_gauges_legend, colorbar_title)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
% 
% ** x: Hydrologic State to show (could be a time-average, a snapshot, ...) 
%
% ** ROUTE_FILE: Location of the RouteLink.nc file that was used in the
%                experiment
%
% ** fig_pos: Position of the figure; could be separate or a subplot
%
% ** show_gauges_legend: Option to display the legend of the gauges
%
% ** colorbar_title: title of the colorbar (usuallu units "cms", "m", ...)
%
%   example: x                  = mean(ncread('all_preassim_mean.nc', 'qlink1'), 2);
%            ROUTE_FILE         = 'RouteLink.nc';
%            fig_pos            = get(gca, 'position');
%            show_gauges_legend = true;
%            colorbar_title     = 'cms';
%
%            plot_connections(x, ROUTE_FILE, fig_pos, show_gauges_legend, colorbar_title);
%
% DART $Id: plot_connections.m $


LAT           = double(ncread(ROUTE_FILE,'lat'))';
LON           = double(ncread(ROUTE_FILE,'lon'))';

fromIndices   = double(ncread(ROUTE_FILE, 'fromIndices'));
fromIndsStart = double(ncread(ROUTE_FILE, 'fromIndsStart'));
fromIndsEnd   = double(ncread(ROUTE_FILE, 'fromIndsEnd'));

[~, n_links]  = netcdf.inqDim(netcdf.open(ROUTE_FILE,'NC_NOWRITE'), 0);
no_up_links   = fromIndsStart == 0;
num_up_links  = fromIndsEnd - fromIndsStart + 1; num_up_links(no_up_links) = 0;

gauges_from_routelink = strtrim(ncread(ROUTE_FILE, 'gages')');
gauges_manual         = gauges_2_indices(gauges_from_routelink);

xL = round(min(LON), 2); xR = round(max(LON), 2);
yB = round(min(LAT), 2); yT = round(max(LAT), 2);

marker_set = {'o', '<', '*', 'd', '>', '^', 'p', 'v', 'h', 's', '.', '+'};

%%
xlabel('Longitude', 'FontSize', 13)
ylabel('Latitude', 'FontSize', 13)

set(gca, 'FontSize', 16, 'XLim', [xL, xR], 'YLim', [yB, yT], ...
         'YTick', linspace(yB, yT, 5), ...
         'YTickLabelRotation', 90); grid on; hold on
    
% make sure we are not in a "single" situation
x  = double(x);
Nx = length(x);

bL = [  30, 144, 255 ]/255;
rD = [ 255,  51,  51 ]/255;
gR = [   0, 153,   0 ]/255;
pR = [ 153,  51, 255 ]/255;
oR = [ 255, 153,  51 ]/255;
lB = [ 153, 255, 255 ]/255;
lR = [ 255, 153, 153 ]/255;
lG = [ 153, 255, 204 ]/255;
lP = [ 204, 153, 255 ]/255;
lO = [ 255, 204, 153 ]/255;

% choose a colormap
if min(x) < -0.1 % increment situation
    cmap = [pR; lP; bL; lB; gR; lG; oR; lO; rD; lR]; 
else
    cmap = jet;
end
size_cmap = size(cmap, 1);

cn = (x - min(x))/(max(x)-min(x));
cn = ceil(cn * size_cmap);
cn = max(cn, 1);

if min(x) < -0.1
    m_cn = mode(cn);
    Np = length(cn(cn >= m_cn));
    Nn = Nx - Np;
    lw = [linspace(1, 3, Np), linspace(1, 3, Nn)];
else
    lw = cn/max(cn) * 3 + 1;
end
 
cmap(size_cmap+1, :) = [150, 150, 150]/255; % Grey color for zero reaches!

zero_links = x > -.001 & x < .001;
cn(zero_links) = size(cmap, 1); 
lw(zero_links) = 1.0;

for i = 1:n_links
    % ith link
    lon_i = LON(i);
    lat_i = LAT(i);
    col_i = cmap(cn(i), :); 
    wid_i = lw(i);
    num_i = num_up_links(i);       

    if num_i > 0
        z  = fromIndices(fromIndsStart(i):fromIndsEnd(i));
        x1 = lon_i*ones(1, num_i); y1 = LON(z); p1 = [x1; y1];
        x2 = lat_i*ones(1, num_i); y2 = LAT(z); p2 = [x2; y2];
        
        S = line(p1, p2); set(S, 'color', col_i, 'LineWidth', wid_i);
    end  
end
 
if length(gauges_manual) <= 10

    Legend_str = cell(1); 
    for k = 1:length(gauges_manual)
        x_loc_gauge = LON( gauges_manual(k,1) );
        y_loc_gauge = LAT( gauges_manual(k,1) );

        leg.H(k) = plot( x_loc_gauge, y_loc_gauge, marker_set{k}, ...
              'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 10); hold on

        Legend_str{k} = [num2str(gauges_manual(k,1)), '; ', num2str(gauges_manual(k,2))];

        %text(x_loc_gauge, y_loc_gauge-0.01, num2str(gauges_manual(k,1)), 'FontSize', 9, 'Color', 'k');
        %text(x_loc_gauge, y_loc_gauge+0.01, num2str(gauges_manual(k,2)), 'FontSize', 9, 'Color', 'k');
    end
    if show_gauges_legend, legend(leg.H, Legend_str, 'Location', 'NorthEast'); end
else

    % Too many gauges: not enough symbols
    for k = 1:length(gauges_manual)
        x_loc_gauge = LON( gauges_manual(k,1) );
        y_loc_gauge = LAT( gauges_manual(k,1) );

        plot( x_loc_gauge, y_loc_gauge, marker_set{1}, ...
              'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 10); hold on
    end
    warning('Number of gauges is large. The image will be cluttered!!')
end

[~, uni_cols] = unique(cn); 
num_colors    = length(uni_cols);
col_levels    = cmap(cn(uni_cols), :);

[~, r] = sort(x(uni_cols));

colormap(gca, col_levels(r,:));

hc = colorbar; 

hc_cubes  = min(6, num_colors); 
hc_range  = linspace(0, num_colors, hc_cubes);
hc_ticks  = hc_range/num_colors;

hc_range(1)  = 1; 
hc_range     = ceil(hc_range); 
hc_full_cols = sort(round(x(uni_cols), 3));

ind_tmp = find(hc_full_cols == 0);
bob     = hc_range - ind_tmp;

pos_tmp = abs(bob) == min( abs(bob) ) ;

hc_range(pos_tmp) = ind_tmp;
hc_ticks(pos_tmp) = ind_tmp/num_colors; 

hc_labels_tmp = hc_full_cols(hc_range);

if max(x) < 10
    if find(pos_tmp) ~= 1 % inc. scenario
        hc_labels = '';
        for k = 1:hc_cubes
            hc_labels(k, :) = sprintf('%9.4f', hc_labels_tmp(k));
        end
        hc_labels(pos_tmp, :) = '< |0.001|';
    else
        hc_labels = '';
        for k = 1:hc_cubes
            hc_labels(k, :) = sprintf('%7.3f', hc_labels_tmp(k));
        end
        if min(x) == 0 && max(x) == 0
            hc_labels(pos_tmp, :) = 'Exact 0';
        else
            hc_labels(pos_tmp, :) = '< 0.001';
        end
    end
else
    hc_labels = hc_labels_tmp;
end

hc_pos(1) = fig_pos(1) + .96*fig_pos(3);
hc_pos(2) = fig_pos(2) + .03*fig_pos(4);
hc_pos(3) = fig_pos(3) * .04;
hc_pos(4) = fig_pos(4) * .25;

set(hc, 'YTick', hc_ticks, 'YTickLabel', hc_labels, 'Position', hc_pos, 'FontSize', 12);    
set(get(hc, 'title'), 'String', colorbar_title)



function obs_ind_id = gauges_2_indices(gauges)

n_links = length(gauges);

k = 0;
for i = 1:n_links
    ob_id = gauges(i, :);
    
    if sum( isspace(ob_id) ) ~= 10
        k = k + 1;
        obs_ind_id(k, :) = [i, str2double(ob_id)]; %#ok
    end
end


% % <next few lines under version control, do not edit>
% % $URL: $
% % $Revision: $
% % $Date: $
