function disp_ATS_localization

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id: disp_ATS_localization.m $

    ROUTE_FILE    = 'Route_Link_fromVars.nc';

    LAT           = double(ncread(ROUTE_FILE,'lat'))';
    LON           = double(ncread(ROUTE_FILE,'lon'))';
    gauges        = ncread(ROUTE_FILE, 'gages')';
    
    fromIndices   = double(ncread(ROUTE_FILE, 'fromIndices'));
    fromIndsStart = double(ncread(ROUTE_FILE, 'fromIndsStart'));
    fromIndsEnd   = double(ncread(ROUTE_FILE, 'fromIndsEnd'));

    [~, n_links]  = netcdf.inqDim(netcdf.open(ROUTE_FILE,'NC_NOWRITE'), 0);
    no_up_links   = fromIndsStart == 0;
    num_up_links  = fromIndsEnd - fromIndsStart + 1; num_up_links(no_up_links) = 0;

    gauges_ind    = [65481, 57281, 65689, 66373, 61903, 62849]';
    num_gauges    = length(gauges_ind);
    gauges_loc    = [LON(gauges_ind); LAT(gauges_ind)] ;

    d = 100;
    c = 0.5*d/6371;
    f = zeros(1, n_links);
            
    for o = 1:num_gauges
        ifile = sprintf('ind_%d.txt' , gauges_ind(o));
        dfile = sprintf('dist_%d.txt', gauges_ind(o));  
        
        load(ifile);
        load(dfile);
        
        l = eval(sprintf('ind_%d' , gauges_ind(o)));
        z = eval(sprintf('dist_%d', gauges_ind(o)));
        for i = 1:length(z)
            f(l(i)) = get_cov_factor(z(i), c);
        end
    end
    
    xL =  -81; xR = -76;
    yB = 33.5; yT =  37;
    
    %% 
    figure('uni','pi','pos', [200, 300, 830, 700])
    
    xlabel('Longitude', 'FontSize', 16)
    ylabel('Latitude', 'FontSize', 16)
    
    set(gca, 'FontSize', 16, 'XLim', [xL, xR], 'YLim', [yB, yT], ...
             'XTick', (xL:.5:xR), 'YTick', (yB:.5:yT)); grid on; hold on

    cmap       = jet;
    cmap(1, :) = [150, 150, 150]/255;
    
    cn = ceil(f * size(cmap, 1));
    cn = max(cn, 1);

    for i = 1:n_links
        % ith link
        lon_i = LON(i);
        lat_i = LAT(i);
        col_i = cmap(cn(i), :);
        num_i = num_up_links(i);

        if num_i > 0
            z  = fromIndices(fromIndsStart(i):fromIndsEnd(i));
            x1 = lon_i*ones(1, num_i); y1 = LON(z); p1 = [x1; y1];
            x2 = lat_i*ones(1, num_i); y2 = LAT(z); p2 = [x2; y2];
            S  = line(p1, p2); set(S, 'color', col_i)
        end  
    end
    
    O1 = plot(gauges_loc(1, 1), gauges_loc(2, 1), '>', ...
             'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    O2 = plot(gauges_loc(1, 2), gauges_loc(2, 2), 'o', ...
             'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    O3 = plot(gauges_loc(1, 3), gauges_loc(2, 3), '^', ...
             'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    O4 = plot(gauges_loc(1, 4), gauges_loc(2, 4), 's', ...
             'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    O5 = plot(gauges_loc(1, 5), gauges_loc(2, 5), '<', ...
             'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    O6 = plot(gauges_loc(1, 6), gauges_loc(2, 6), 'd', ...
             'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 8);   

    [~, uni_cols] = unique(cn); 
    num_colors    = length(uni_cols);
    col_levels    = cmap(cn(uni_cols), :);
    
    colormap(col_levels);

    hc = colorbar; 

    hc_cubes  = 6; 
    hc_range  = linspace(0, num_colors, hc_cubes);
    hc_ticks  = hc_range/num_colors;

    hc_full_cols = f(uni_cols); hc_range(1) = 1; hc_range = ceil(hc_range);
    hc_labels    = round(hc_full_cols(hc_range), 1);

    fig_pos   = get(gca, 'position');
    hc_pos(1) = fig_pos(1)+0.926*fig_pos(4);
    hc_pos(2) = fig_pos(2)+0.250*fig_pos(3);
    hc_pos(3) = fig_pos(3)*0.04;
    hc_pos(4) = fig_pos(4)*0.25;

    set(hc, 'YTick', hc_ticks, 'YTickLabel', hc_labels, ...
            'Position', hc_pos);
    set(get(hc, 'title'), 'String', {'Covariance','Factor'})
    
    LEG = legend([O1, O2, O3, O4, O5, O6], ...
           ['obs: gauge ID ' strtrim( gauges(gauges_ind(1), :) ) ], ...
           ['obs: gauge ID ' strtrim( gauges(gauges_ind(2), :) ) ], ...
           ['obs: gauge ID ' strtrim( gauges(gauges_ind(3), :) ) ], ...
           ['obs: gauge ID ' strtrim( gauges(gauges_ind(4), :) ) ], ...
           ['obs: gauge ID ' strtrim( gauges(gauges_ind(5), :) ) ], ...
           ['obs: gauge ID ' strtrim( gauges(gauges_ind(6), :) ) ], ...
           'Location', 'SouthEast');

    set(LEG, 'EdgeColor', 'w')  
    title(['Florence Domain: Along-the-Stream Localization (Distance = ' num2str(d) 'km)'], 'FontSize', 18, 'FontWeight', 'normal')
    
    %saveas(gcf, 'localization_visualizer.png') 
end

% <next few lines under version control, do not edit>
% $URL: $
% $Revision: $
% $Date: $
