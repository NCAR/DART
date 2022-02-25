function plot_connections(x, ROUTE_FILE, fig_pos, colorbar_title, tiny_flow)

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
% ** colorbar_title: title of the colorbar (usually units "cms", "m", ...)
%
% ** tiny_flow: low flow conditions (e.g., 10cm for streamflow, 1mm for bucket)
%
% NOTE: THIS FUNCTION REQUIRES 'cbrewer'. THIS CAN BE DOWNLOADED FROM: 
% https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab 
%
%   example: x                  = mean(ncread('all_preassim_mean.nc', 'qlink1'), 2);
%            ROUTE_FILE         = 'RouteLink.nc';
%            fig_pos            = get(gca, 'position');
%            colorbar_title     = 'cms';
%            tiny_flow          = 10;
%
%            plot_connections(x, ROUTE_FILE, fig_pos, colorbar_title, tiny_flow);
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

xL = -80.7; %round(min(LON), 2); 
xR = -76.1; %round(max(LON), 2);
yB = 33.60; %round(min(LAT), 2); 
yT = 36.90; %round(max(LAT), 2);


%%
xlabel('Longitude', 'FontSize', 13)
ylabel('Latitude', 'FontSize', 13)

set(gca, 'FontSize', 16, 'XLim', [xL, xR], 'YLim', [yB, yT], ...
         'YTick', [34.4250, 35.2500, 36.0750], ...
         'YTickLabelRotation', 90); grid on; hold on
    
% make sure we are not in a "single" situation
x  = double(x);

pos_links = sort(x(x>tiny_flow));
zer_links = x>=-tiny_flow & x<=tiny_flow;
neg_links = sort(x(x<-tiny_flow));

pos_len = length(pos_links);
neg_len = length(neg_links);

if numel(find(x<0)) / n_links * 100 > 5  % negative values less than 5% of all links

    cmap = [ ...
             flipud(cbrewer('seq', 'Reds', length(neg_links), 'PCHIP')); ...
             [.8, .8, .8]; ...
             cbrewer('seq', 'Blues', length(pos_links), 'PCHIP') ...
           ]; 
    size_cmap = size(cmap, 1);
else
    mval = abs(x);
    cmap = hsv(n_links); %sstpal(n_links) %cbrewer('div', 'Spectral', n_links, 'PCHIP');
    
    size_cmap = size(cmap, 1);
    
    cmap(1, :) = [.8, .8, .8]; % Grey color for zero reaches!    
end

z = x;
y = abs(x);
b = (z - min(z))/(max(z)-min(z));
d = (y - min(y))/(max(y)-min(y));

cn = ceil(b * size_cmap);
cn = max(cn, 1);
lw = d * 8 + .1;
 
cn(zer_links) = neg_len + 1; 

for i = 55000:n_links %1:n_links (to speed things up, only show major streams)
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

if pos_len > 0 && neg_len > 0 

    colormap(gca, cmap)
    hc = colorbar; 
    
    hc1 = neg_len / size_cmap;

    hc_major = [0, neg_len, neg_len+1, neg_len+1+pos_len];

    zer_point = hc_major(3);
    if hc1 > 0.5
        neg_points = [1, ceil(neg_len/3), ceil(2*neg_len/3)];
        neg_labels = neg_links(neg_points)';

        pos_points = [ceil(pos_len/2), pos_len];
        pos_labels = pos_links(pos_points)';

    else

        neg_points = [1, ceil(neg_len/2)];
        neg_labels = neg_links(neg_points)';

        pos_points = [ceil(pos_len/3), ceil(2*pos_len/3), pos_len];
        pos_labels = pos_links(pos_points)';

    end
    hc_labels  = [neg_labels, 0, pos_labels];
    hc_labels  = round(hc_labels, 3, 'significant');

    hc_locate  = [hc_major(1), neg_points(2:end), ...
                    zer_point, pos_points+neg_len+1] / size_cmap;
                
else
    
    [~, uni_cols] = unique(cn); 
    num_colors    = length(uni_cols);
    col_levels    = cmap(cn(uni_cols), :);

    [~, r] = sort(mval(uni_cols));

    colormap(gca, col_levels(r,:));

    hc = colorbar; 

    hc_cubes   = min(7, num_colors); 
    hc_range   = linspace(0, num_colors, hc_cubes);
    hc_locate  = hc_range/num_colors;

    hc_range(1)  = 1; 
    hc_range     = ceil(hc_range); 
    hc_full_cols = sort(round(mval(uni_cols), 5));

    ind_tmp = find(hc_full_cols == 0);
    bob     = hc_range - ind_tmp;

    if size(ind_tmp, 1) > 0
        pos_tmp = find(abs(bob) == min(abs(bob))) ;
        pos_tmp = pos_tmp(end);

        hc_range(pos_tmp) = ind_tmp;
    end

    hc_labels_tmp = hc_full_cols(hc_range);

    hc_labels = round(hc_labels_tmp, 3, 'significant');
    
end

hc_pos(1) = fig_pos(1) + .96*fig_pos(3);
hc_pos(2) = fig_pos(2) + .03*fig_pos(4);
hc_pos(3) = fig_pos(3) * .04;
hc_pos(4) = fig_pos(4) * .30;

set(hc, 'YTick', hc_locate, 'YTickLabel', hc_labels, 'Position', hc_pos, 'FontSize', 12)

set(get(hc, 'title'), 'String', colorbar_title, 'FontSize', 14)

% text(-80.4, 34.0, 'SC', 'FontSize', 20, 'FontWeight', 'Bold')
% text(-80.4, 36.0, 'NC', 'FontSize', 20, 'FontWeight', 'Bold')
% text(-80.4, 36.7, 'VA', 'FontSize', 20, 'FontWeight', 'Bold')



% % <next few lines under version control, do not edit>
% % $URL: $
% % $Revision: $
% % $Date: $

