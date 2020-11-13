function hpol = polar_dares(THETA, RHO, SPEC)
% This is a modified version of MATLAB's polar function. 
% The modifications are done in such a way that the visualization of the
% DART tutorial is faster. 
%       - The X- and Y- labels are turned off. 
%       - The function allows for zooming 

global MEAN_DIST

    if nargin > 2
        varargin{1} = THETA;
        varargin{2} = RHO;
        varargin{3} = SPEC;
    else
        varargin{1} = THETA;
        varargin{2} = RHO;
    end
    
    % Parse possible Axes input
    [cax, args, nargs] = axescheck(varargin{:});
    
    if nargs < 1
        error(message('MATLAB:narginchk:notEnoughInputs'));
    elseif nargs > 3
        error(message('MATLAB:narginchk:tooManyInputs'));
    end
    
    if nargs < 1 || nargs > 3
        error(message('MATLAB:polar:InvalidDataInputs'));
    elseif nargs == 2
        theta = args{1};
        rho = args{2};
        if ischar(rho)
            line_style = rho;
            rho = theta;
            [mr, nr] = size(rho);
            if mr == 1
                theta = 1 : nr;
            else
                th = (1 : mr)';
                theta = th(:, ones(1, nr));
            end
        else
            line_style = 'auto';
        end
    elseif nargs == 1
        theta = args{1};
        line_style = 'auto';
        rho = theta;
        [mr, nr] = size(rho);
        if mr == 1
            theta = 1 : nr;
        else
            th = (1 : mr)';
            theta = th(:, ones(1, nr));
        end
    else % nargs == 3
        [theta, rho, line_style] = deal(args{1 : 3});
    end
    if ischar(theta) || ischar(rho)
        error(message('MATLAB:polar:InvalidInputType'));
    end
    if ~isequal(size(theta), size(rho))
        error(message('MATLAB:polar:InvalidInputDimensions'));
    end
    try
        theta = full(double(theta));
        rho = full(double(rho));
    catch
        error(message('MATLAB:specgraph:private:specgraph:nonNumericInput'));
    end
    
    % get hold state
    cax = newplot(cax);
    
    next = lower(get(cax, 'NextPlot'));
    hold_state = ishold(cax);

    if isa(handle(cax),'matlab.graphics.axis.PolarAxes')
        error(message('MATLAB:polar:PolarAxes'));
    end
    
    % make the concentric circles a bit darker gray
    % get x-axis text properties
    tc = [0.7, 0.7, 0.7]; 
    ls = get(cax, 'GridLineStyle');

    % Hold on to current Text defaults, reset them to the
    % Axes' font attributes so tick marks use them.
    fAngle = get(cax, 'DefaultTextFontAngle');
    fName = get(cax, 'DefaultTextFontName');
    fSize = get(cax, 'DefaultTextFontSize');
    fWeight = get(cax, 'DefaultTextFontWeight');
    fUnits = get(cax, 'DefaultTextUnits');
    set(cax, ...
        'DefaultTextFontAngle', get(cax, 'FontAngle'), ...
        'DefaultTextFontName', get(cax, 'FontName'), ...
        'DefaultTextFontSize', get(cax, 'FontSize'), ...
        'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
        'DefaultTextUnits', 'data');
    
    % only do grids if hold is off
    if ~hold_state
        
        % make a radial grid
        hold(cax, 'on');
        % ensure that Inf values don't enter into the limit calculation.
        arho = abs(rho(:));
        maxrho = max(arho(arho ~= Inf));
        hhh = line([-maxrho, -maxrho, maxrho,  maxrho], ...
                   [-maxrho,  maxrho, maxrho, -maxrho], 'Parent', cax);
        set(cax, 'DataAspectRatio', [1, 1, 1], 'PlotBoxAspectRatioMode', 'auto');
        v = [get(cax, 'XLim') get(cax, 'YLim')];
        ticks = sum(get(cax, 'YTick') >= 0);
        delete(hhh);
        % check radial limits and ticks
        rmin = 0;
        rmax = v(4);
        rticks = max(ticks - 1, 2);
        if rticks > 5   % see if we can reduce the number
            if rem(rticks, 2) == 0
                rticks = rticks / 2;
            elseif rem(rticks, 3) == 0
                rticks = rticks / 3;
            end
        end
        
        % define a circle
        th = 0 : pi / 50 : 2 * pi;
        xunit = cos(th);
        yunit = sin(th);
        % now really force points on x/y axes to lie on them exactly
        inds = 1 : (length(th) - 1) / 4 : length(th);
        xunit(inds(2 : 2 : 4)) = zeros(2, 1);
        yunit(inds(1 : 2 : 5)) = zeros(3, 1);
        % plot background if necessary
        if ~ischar(get(cax, 'Color'))
            patch('XData', xunit * rmax, 'YData', yunit * rmax, ...
                'EdgeColor', tc, 'FaceColor', get(cax, 'Color'), ...
                'HandleVisibility', 'off', 'Parent', cax);
        end
        
        % draw radial circles (remove annotation; not useful here!)
        c82 = cos(80 * pi / 180);
        s82 = sin(80 * pi / 180);
        rinc = (rmax - rmin) / rticks;
        for i = (rmin + rinc) : rinc : rmax - 1
            hhh = line(xunit * i, yunit * i, 'LineStyle', ls, 'Color', tc, ...
                'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', cax);
            text((i - 2 + rinc / 20) * c82, (i - 2 + rinc / 20) * s82, ...
                ['  ' num2str(i-MEAN_DIST)], 'VerticalAlignment', 'bottom', ...
                'HandleVisibility', 'off', 'Parent', cax, 'FontSize', 10);
        end
        set(hhh, 'LineStyle', '-'); % Make outer circle solid
        
        % plot spokes
        th = (1 : 4) * 2 * pi / 8;
        cst = cos(th);
        snt = sin(th);
        cs = [-cst; cst];
        sn = [-snt; snt];
        line(rmax * cs, rmax * sn, 'LineStyle', ls, 'Color', tc, 'LineWidth', 1, ...
            'HandleVisibility', 'off', 'Parent', cax);
        
        % annotate spokes in degrees
        rt = 1.1 * rmax;
        grid1 = [ 5, 10, 15, 20];  
        grid2 = [25, 30, 35, 40];
        for i = 1 : length(th)
            text( rt * cst(i),  rt * snt(i), int2str(grid1(i)), ...
                'HorizontalAlignment', 'center', ...
                'HandleVisibility', 'off', 'Parent', cax, 'FontSize', 14); 
            text(-rt * cst(i), -rt * snt(i), int2str(grid2(i)), ...
                'HorizontalAlignment', 'center', ...
                'HandleVisibility', 'off', 'Parent', cax, 'FontSize', 14); 
        end
        
        % set view to 2-D
        view(cax, 2);
        % set axis limits
        axis(cax, rmax * [-1, 1, -1.15, 1.15]);
    end
    
    % Reset defaults.
    set(cax, ...
        'DefaultTextFontAngle', fAngle , ...
        'DefaultTextFontName', fName , ...
        'DefaultTextFontSize', fSize, ...
        'DefaultTextFontWeight', fWeight, ...
        'DefaultTextUnits', fUnits );
    
    % transform data to Cartesian coordinates.
    xx = rho .* cos(theta);
    yy = rho .* sin(theta);
    
    % plot data on top of grid
    if strcmp(line_style, 'auto')
        q = plot(xx, yy, 'Parent', cax);
    else
        q = plot(xx, yy, line_style, 'Parent', cax);
    end
    
    if nargout == 1
        hpol = q;
    end
    
    if ~hold_state
        set(cax, 'DataAspectRatio', [1, 1, 1]), axis(cax, 'off');
        set(cax, 'NextPlot', next);
    end
    
    % Disable pan and zoom
    p = hggetbehavior(cax, 'Pan');
    p.Enable = true;
    z = hggetbehavior(cax, 'Zoom');
    z.Enable = true;
    
    if ~isempty(q) && ~isdeployed
        makemcode('RegisterHandle', cax, 'IgnoreHandle', q, 'FunctionName', 'polar');
    end
end

