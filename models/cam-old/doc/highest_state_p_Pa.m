function highest_state_p_Pa(fname, diffusion_flag, cutoff, vert_norm, varargin)
% plots the vertical damping profile for CAM, WACCM
%
%  highest_state_p_Pa(fname, diffusion_flag, cutoff, vert_norm)
%  highest_state_p_Pa(fname, diffusion_flag, cutoff, vert_norm, 'nu_top', nu_top)
%
%  fname = pathname of initial files containing vertical grid information.
%
%  diffusion_flag = CAM namelist div24del2flag value, except for 0, which is not used by CAM;
%                  0  for CAM-SE,
%                  2  for CAM-FV  div2 damping
%                  4  for CAM-FV  div4 damping
%                  24 for CAM-FV  del2+div4 damping
%
%  cutoff = DART namelist variable 'cutoff'
%
%  vert_norm = DART namelist variable 'vert_normalization_pressure'
%
%  Optional arguments:
%  nu_top = needed only for CAM-SE (diffusion_flag = 0).
%              Get from CAM's atm_in namelist.
%
% Plot CAM's (WACCM's) extra diffusion in the top layers
% These are functions of CAM's ptop and -FV's div24del2flag.
% Plot cam/model_mod.f90's highest_state_pressure_Pa and distance increment added to
% the nominal distances between an ob and a state variable location,
% which are functions of DART's cutoff and vert_normalization_pressure,
% as well as CAM's ptop and div24del2flag.
%
% EXAMPLE:
% fname         = 'caminput.nc';
% diffusion_flag = 2;
% cutoff        = 0.2;      % DART input.nml cutoff
% vert_norm     = 20000.0;  % DART input.nml vert_normalization_pressure
% highest_state_p_Pa(fname, diffusion_flag, cutoff, vert_norm)
%
% EXAMPLE:
% fname         = 'caminput.nc';
% diffusion_flag = 0;
% cutoff        = 0.2;      % DART input.nml cutoff
% vert_norm     = 20000.0;  % DART input.nml vert_normalization_pressure
% highest_state_p_Pa(fname, diffusion_flag, cutoff, vert_norm,'nu_top',0.00005)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

default_nu_top = 2.0e5;
args = inputParser;

addRequired(args,'fname',@ischar);
addRequired(args,'diffusion_flag',@isnumeric);
addRequired(args,'cutoff',@isnumeric);
addRequired(args,'vert_norm',@isnumeric);
addParamValue(args,'nu_top',default_nu_top,@isnumeric);
parse(args,fname,diffusion_flag,cutoff,vert_norm,varargin{:});

if ~isempty(fieldnames(args.Unmatched))
    disp('Extra inputs:')
    disp(args.Unmatched)
end

% if you want to echo the input
fprintf('fname                        : %s\n', args.Results.fname)
fprintf('diffusion_flag               : %i\n', args.Results.diffusion_flag)
fprintf('cutoff                       : %f radians\n', args.Results.cutoff)
fprintf('vert_normalization_pressure  : %f Pa/radian\n', args.Results.vert_norm)
if (diffusion_flag == 0)
    fprintf('nu_top                       : %f\n', args.Results.nu_top)
end

%   Read in nlevs, hybrid coords (As and Bs) from the initial file.
if (exist(fname,'file') ~= 2)
    error('file/fname <%s> does not exist',fname)
end

if (diffusion_flag == 4)
    error('diffusion flavor 4 ("div4") is independent of height.  Can set highest_state_pressure_Pa = 0')
end

As    = local_nc_varget(fname, 'hyam');
Bs    = local_nc_varget(fname, 'hybm');
fulls = local_nc_varget(fname, 'hyai');

ncid = netcdf.open(fname,'NOWRITE');
dimid = netcdf.inqDimID(ncid,'lev');
[~, nlevs] = netcdf.inqDim(ncid, dimid);
netcdf.close(ncid);

% Convert cutoff from radians to Pa.
% For illustration; DART does all of this in radian space.

c_v  = cutoff*vert_norm;
fprintf('cutoff*vert_norm is %f\n',c_v)

% Calculate p(k=1,nlevs) from As and Bs,

ptop  = 1.0e5 * fulls(1);
p     = 1.0e5 * (As + Bs);


if (diffusion_flag == 0)
    % Calculate the "regular diffusion" in the top 3 layers for CAM-SE
    diff_bot_k = 3;
    % Calculate the diffusion coefficient without height restrictions.
    tau = zeros(nlevs,1);
    tau(1) = args.Results.nu_top * 4;
    tau(2) = args.Results.nu_top * 2;
    tau(3) = args.Results.nu_top;
    
elseif (diffusion_flag == 2)
    % Calculate the diffusion coefficient without height restrictions.
    tau = 8.0*(1.0 + tanh(log(ptop./p)));
    % Calculate the div2 diffusion coefficient from p
    % and the pressure level above which distance increments will be added.
    tau_div2 = ones(nlevs,1);
    for k = 1:nlevs
        if (tau(k) < 1.0)
            diff_bot_k = k-1;
            break
        else
            tau_div2(k) = tau(k);
        end
    end
elseif (diffusion_flag == 24)
    % Calculate the diffusion coefficient without height restrictions.
    tau = 8.0*(1.0 + tanh(log(ptop./p)));
    % Calculate the del2 diffusion coefficient from p.
    tau_del2 = zeros(nlevs,1);
    for k = 1:nlevs
        if (tau(k) < 0.3)
            diff_bot_k = k-1;
            break
        else
            tau_del2(k) = tau(k);
        end
    end
else
    error('diffusion_flag not recognized.  Choose from {0,2,24,4}')
end

highest_spPa = p(diff_bot_k) + c_v*(1.0 + sqrt(1.0 + 2*(p(diff_bot_k) - ptop)./(c_v)));

% Distance increment in units of Pa (in model_mod it's radians)
delta_d = zeros(nlevs,1);
for k = 1:nlevs
    if (p(k) < highest_spPa)
        delta_d(k) = (highest_spPa - p(k))^2 ./ (highest_spPa - ptop);
    else
        % Save a k-range for plotting only the interesting pressure range.
        %highest_spPa_k = k+3;
        highest_spPa_k = nlevs;
        break
    end
end

clf; orient tall; wysiwyg

% little summary of the input at the bottom
axes('position',[0.1, 0.03, 0.8, 0.1125]);
axis off
axis([0 1 0 1.5])
h1 = text(0.1, 0.00, sprintf('file providing model levels = %s',args.Results.fname));
h2 = text(0.1, 0.25, sprintf('diffusion_flag              = %i',args.Results.diffusion_flag));
h3 = text(0.1, 0.50, sprintf('cutoff                      = %f radians',args.Results.cutoff));
h4 = text(0.1, 0.75, sprintf('vert_normalization_pressure = %f Pa/radian',args.Results.vert_norm));
if (diffusion_flag == 0)
    h5 = text(0.1, 1.00, sprintf('nu_top                     = %f',args.Results.nu_top));
    set(h5,'Interpreter','none','FontName','Courier New','FontSize',12)
end
set(h1,'Interpreter','none','FontName','Courier New','FontSize',12)
set(h2,'Interpreter','none','FontName','Courier New','FontSize',12)
set(h3,'Interpreter','none','FontName','Courier New','FontSize',12)
set(h4,'Interpreter','none','FontName','Courier New','FontSize',12)

% Plot on a reversed log10 scale
ax1 = axes('position',[0.125, 0.2, 0.75, 0.7]);
set(ax1, 'Ydir','reverse', ...
    'Yscale','log', ...
    'Layer','top', ...
    'FontSize',20)

hold on;
ylabel('Pressure (Pa)')

if (diffusion_flag == 0)
    % Plot (2) (div2_tau vs p(k)
    h = plot(tau(1:highest_spPa_k), p(1:highest_spPa_k), 'k.-', ...
        'LineWidth',2.0, 'MarkerSize',20);
    ptop_x = [0.0, max(tau)*1.05];
elseif (diffusion_flag == 2)
    % Plot (2) (div2_tau vs p(k)
    h = plot(tau_div2(1:highest_spPa_k), p(1:highest_spPa_k), 'k.-', ...
        'LineWidth',2.0, 'MarkerSize',20);
    ptop_x = [0.0, max(tau_div2)*1.05];
elseif (diffusion_flag == 24)
    h = plot(tau_del2(1:highest_spPa_k), p(1:highest_spPa_k), 'k.-', ...
        'LineWidth',2.0, 'MarkerSize',20);
    ptop_x = [0.0, max(tau_del2)*1.05];
end

% A line at p(top)

ptop_y = [ptop,          ptop];
ptop_line = line(ptop_x,ptop_y,'Color','k');
set(ptop_line,'LineStyle','--','Marker','none', 'LineWidth',2.0);

% A line at highest_state_pressure_Pa for default cutoff*vert_norm
% Other cutoff*vert_norm  in a lighter shade?  NO; handle in plot of H_Pa vs c*v, below

highest_spPa_y = [highest_spPa, highest_spPa];
highest_spPa_line = line(ptop_x,highest_spPa_y,'Color','r');
set(highest_spPa_line,'LineStyle',':','Marker','none', 'LineWidth',2.0);

xlabel('diffusion coefficient \tau')

% Second x-axis (on top) for plotting distance increments.

ax2_XLim = [-.05*max(delta_d), 1.05*max(delta_d)];
ax2 = axes('position',get(ax1,'Position'), ...
    'FontSize',20, ...
    'XAxisLocation','top', ...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','r','YColor','k',...
    'XLim',ax2_XLim, ...
    'YLim',get(ax1,'YLim'), ...
    'YDir',get(ax1,'YDir'), ...
    'Layer','top', ...
    'Yscale','log');
% set(ax2, 'Ydir','reverse', 'Yscale','log', 'Layer','top')

ylabel('Pressure (Pa)')
xlabel('Extra distance used in cutoff calculation (Pa)')

% A curve of delta-d above highest_state_pressure_Pa for default cutoff*vert_norm (only)
d_line = line(delta_d(1:highest_spPa_k),p(1:highest_spPa_k),'Color','r');
set(d_line,'Marker','x', 'LineWidth',2.0);

% A line at 2*cutoff
cutoff_x      = [2*c_v, 2*c_v];
y_axis_height = get(ax2,'YLim');
cutoff_y      = [y_axis_height(1), p(nlevs-5)];
cutoff_line   = line(cutoff_x, cutoff_y, 'Color','r');
set(cutoff_line,'LineStyle','--','Marker','none', 'LineWidth',2.0);

% Plot (min) highest_state_pressure_Pa vs cutoff*vert_norm.

legend_handles = [cutoff_line, d_line, highest_spPa_line, ptop_line, h];
lh = legend(legend_handles,'obs impact limit','extra "distance"', ...
    'highest\_state\_pressure\_Pa','pressure at model top','CAM \tau');
set(lh,'Interpreter','tex','Location','SouthEast','box','off')

% create a png file for import to powerpoint

outfname = sprintf('tau_dist_c%3.2f_v%i.png', cutoff,vert_norm);
print(gcf,'-dpng',outfname);

%=====================================================================


function value = local_nc_varget(fname,varname)
%% If the variable exists in the file, return the contents of the variable.
% if the variable does not exist, return empty value instead of error-ing
% out.

[variable_present, varid] = nc_var_exists(fname,varname);
if (variable_present)
    ncid  = netcdf.open(fname,'NOWRITE');
    value = netcdf.getVar(ncid, varid);
    netcdf.close(ncid)
else
    value = [];
end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
