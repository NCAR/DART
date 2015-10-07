function PlotSawtooth( pinfo )
%% PlotSawtooth: Plots time series of ensemble members, mean and truth.
%
%
% PlotSawtooth is intended to be called by 'plot_sawtooth'
% The only input argument is a structure with model-dependent
% components.
%
% Both the prior and posterior estimates are plotted as a single
% trajectory. If there is no change in the model state, this should
% appear as a series of steps. This necessitates plotting the 'posterior'
% first ... think about it ...
%
% USAGE: PlotSawtooth( pinfo );
%
% STRUCTURE COMPONENTS FOR low-order models
% truth_file      name of netCDF DART file with copy tagged 'true state'
% prior_file      name of forecast netCDF DART file with desired copy'
% posterior_file  name of nudged netCDF DART file with desired copy'
% var             name of netCDF variable of interest
% var_inds        indices of variables of interest
%
% Example 1 ( forced_lorenz_96 model )
%%--------------------------------------------------------
% pinfo.truth_file     = 'True_State.nc';
% pinfo.prior_file     = 'Prior_Diag.nc';
% pinfo.posterior_file = 'Posterior_Diag.nc';
% pinfo.var            = 'state';
% pinfo.var_inds       = [ 23 36 42 ];
% PlotSawtooth( pinfo );
%
% Example 2 ( fms_bgrid_model )
%%--------------------------------------------------------
% pinfo.truth_file     = 'True_State.nc';
% pinfo.prior_file     = 'Prior_Diag.nc';
% pinfo.posterior_file = 'Posterior_Diag.nc';
% pinfo.var        = 'u';
% pinfo.level      = 3;
% pinfo.latitude   = 23.5;
% pinfo.longitude  = 45.67;
% PlotSawtooth( pinfo )

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Get some information from the truth_file, if it exists.
if ( exist(pinfo.truth_file,'file') == 2 )
    truth = CheckModelCompatibility(pinfo.truth_file, pinfo.posterior_file);
    truth.truth_index = get_copy_index(pinfo.truth_file, 'true state' );
else
    truth = [];
end

%% Get some information from the prior_file
%  The metadata is queried to determine which "copy" is appropriate
%  and a 'doubled up' x axis plotting array is created.

prior.ens_mean_index = get_copy_index(pinfo.prior_file,     'ensemble mean');
post.ens_mean_index = get_copy_index(pinfo.posterior_file, 'ensemble mean');

x          = zeros(2,pinfo.time_series_length);
x(1,:)     = pinfo.time;
x(2,:)     = pinfo.time;
pinfo.xax  = x(:);
metadata   = nc_varget(pinfo.prior_file,'CopyMetaData');

%% The (usually simple) model states that are not stored in prognostic
%  variables are plotted with thei indices in to the StateVector, if the
%  model state has been parsed into prognostic variables, we use the second
%  routine.

if isfield(pinfo,'var_inds')
    PlotGivenIndices( pinfo, prior, post, truth, metadata);
elseif isfield(pinfo,'var_names')
    PlotGivenVariable(pinfo, truth, metadata);
end



function PlotGivenIndices(pinfo, prior, post, truth, metadata)

%% Plot given an index into the (low-order-model) state-space vector.
%  Each variable gets its own figure.

iplot = 0;

for ivar = pinfo.var_inds,

    % Get the data from the netcdf files

    po_ens_mean = get_hyperslab('fname',pinfo.posterior_file, 'varname',pinfo.var, ...
        'copyindex',post.ens_mean_index, 'stateindex',ivar, ...
        'tindex1',pinfo.posterior_time(1), 'tcount',pinfo.posterior_time(2));
    pr_ens_mean = get_hyperslab('fname',pinfo.prior_file, 'varname',pinfo.var, ...
        'copyindex',prior.ens_mean_index, 'stateindex',ivar, ...
        'tindex1',pinfo.prior_time(1), 'tcount',pinfo.prior_time(2));

    % Now we paste them together in a clever way to show
    % the effect of the assimilation
    ens_mean(1,:) = pr_ens_mean;
    ens_mean(2,:) = po_ens_mean;
    a             = ens_mean(:);

    % Plot the true trajectory if it exists; the ens mean; annotate

    iplot = iplot + 1;
    figure(iplot); clf;

    if ( exist(pinfo.truth_file,'file') == 2 )

        true_trajectory = get_hyperslab('fname',pinfo.truth_file, 'varname',pinfo.var, ...
            'copyindex',truth.truth_index, 'stateindex',ivar);

        plot(truth.time, true_trajectory, 'k-','linewidth',1.0); hold on;
        plot(pinfo.xax, a, 'k-','linewidth',2.0);
        legend('truth','ensemble mean')
    else
        plot(pinfo.xax, a, 'k-','linewidth',2.0); hold on;
        legend('ensemble mean')
    end

    ylabel(sprintf('''%s'' index %d', pinfo.var, ivar ))
    xlabel(sprintf('model "days" (%d timesteps)',pinfo.time_series_length))
    title(sprintf('%s Trajectories',pinfo.model), ...
        'interpreter','none','fontweight','bold')

    % Now check to see if we are overlaying any individual ensemble members.
    % if pinfo.copyindices = [], nothing happens.

    ens_colors = get(gca,'ColorOrder');   % trying to cycle through colors
    ncolors = size(ens_colors,1) - 1;     % last one is black, already used.

    nmem = 0;
    for imem = pinfo.copyindices
        nmem = nmem + 1;

        str1 = deblank(metadata(imem,:));
        copy_index = get_copy_index(pinfo.prior_file, str1);

        po_series  = get_hyperslab('fname',pinfo.posterior_file, 'varname',pinfo.var, ...
            'copyindex',copy_index, 'stateindex',ivar, ...
            'tindex1',pinfo.posterior_time(1), 'tcount',pinfo.posterior_time(2));
        pr_series  = get_hyperslab('fname',pinfo.prior_file, 'varname',pinfo.var, ...
            'copyindex',copy_index, 'stateindex',ivar, ...
            'tindex1',pinfo.prior_time(1), 'tcount',pinfo.prior_time(2));

        ens_member(1,:) = pr_series;
        ens_member(2,:) = po_series;
        b               = ens_member(:);

        hold on;
        memcolor = 1 + mod(nmem-1,ncolors); % cycles through colors [1,6]
        h = plot(pinfo.xax, b,'linewidth',0.5,'Color',ens_colors(memcolor,:));

        [~, ~, outh, outm] = legend;
        nlines             = length(outm);
        outm{nlines+1}     = str1;
        legend([outh; h],outm,0);
    end
    legend boxoff

end



function PlotGivenVariable(pinfo, truth, metadata)
%% Plot given an variable name.

var_names = parse(pinfo.var_names);

iplot = 0;

for ivar = 1:length(var_names)

    iplot = iplot + 1;
    figure(iplot); clf;
    ens_colors = get(gca,'ColorOrder');   % trying to cycle through colors
    ncolors    = size(ens_colors,1) - 1;  % last one is black, already used.

    vname = var_names{ivar};
    handles = zeros(1,pinfo.copies);
    strings = cell(1,pinfo.copies);
    for i = 1:pinfo.copies

        %% multiple copies can get overlain on same axis

        imem = pinfo.copyindices(i);

        if isfield(pinfo,'cellindex')
            pr_series  = get_hyperslab('fname',pinfo.prior_file, 'varname',vname, ...
                'tindex1',pinfo.prior_time(1), 'tcount',pinfo.prior_time(2), ...
                'cellindex',pinfo.cellindex, ...
                'levelindex',pinfo.levelindex, 'copyindex',imem);
            po_series  = get_hyperslab('fname',pinfo.posterior_file, 'varname',vname, ...
                'tindex1',pinfo.posterior_time(1), 'tcount',pinfo.posterior_time(2), ...
                'cellindex',pinfo.cellindex, ...
                'levelindex',pinfo.levelindex, 'copyindex',imem);
        else

            pr_series  = get_hyperslab('fname',pinfo.prior_file, 'varname',vname, ...
                'tindex1',pinfo.prior_time(1), 'tcount',pinfo.prior_time(2), ...
                'lonindex',pinfo.lonindex,'latindex',pinfo.latindex, ...
                'levelindex',pinfo.levelindex, 'copyindex',imem);
            po_series  = get_hyperslab('fname',pinfo.posterior_file, 'varname',vname, ...
                'tindex1',pinfo.posterior_time(1), 'tcount',pinfo.posterior_time(2), ...
                'lonindex',pinfo.lonindex,'latindex',pinfo.latindex, ...
                'levelindex',pinfo.levelindex, 'copyindex',imem);
        end

        % Paste Prior/Posterior into one series

        ens_member(1,:) = pr_series;
        ens_member(2,:) = po_series;
        b               = ens_member(:);
        strings{i}      = deblank(metadata(imem,:));

        hold on;
        memcolor = 1 + mod(i-1,ncolors); % cycles through colors [1,6]

        if (strcmp(strings{i},'ensemble mean') == 1)
            handles(i) = plot(pinfo.xax, b,'k-','linewidth',2.0);
        else
            handles(i) = plot(pinfo.xax, b,'linewidth',0.5,'Color',ens_colors(memcolor,:));
        end

    end

    %% Plot the true trajectory if it exists

    if ( exist(pinfo.truth_file,'file') == 2 )      % plot it

        if isfield(pinfo,'cellindex')
            true_trajectory  = get_hyperslab('fname',pinfo.truth_file, 'varname',vname, ...
                'tindex1',truth.truth_time(1), 'tcount',truth.truth_time(2), ...
                'cellindex',pinfo.cellindex, ...
                'levelindex',pinfo.levelindex, 'copyindex',truth.truth_index);
        else
            true_trajectory  = get_hyperslab('fname',pinfo.truth_file, 'varname',vname, ...
                'tindex1',truth.truth_time(1), 'tcount',truth.truth_time(2), ...
                'lonindex',pinfo.lonindex,'latindex',pinfo.latindex, ...
                'levelindex',pinfo.levelindex, 'copyindex',truth.truth_index);
        end

        h = plot(truth.time, true_trajectory, 'k-*','linewidth',1.0); hold on;

        strings{length(strings)+1} = 'truth';
        handles = [handles h];
    end

    % annotate

    h = ylabel(sprintf('''%s''  (%.3f,%.3fE) level index %d', vname, ...
        pinfo.latitude, pinfo.longitude, pinfo.levelindex ));
    set(h,'Interpreter','none')
    xdates(pinfo.time)
    title(sprintf('%s Trajectories',pinfo.model), ...
        'Interpreter','none','fontweight','bold')
    legend(handles,strings);
    legend boxoff

end



function xdates(dates)
if (length(dates) < 5)
    set(gca,'XTick',dates);
    datetick('x',31,'keepticks','keeplimits');
    xlabel('Model date (YYYY-MM-DD HH:MM:SS)')
else
    datetick('x','mm.dd.HH','keeplimits'); % 'mm/dd'
    monstr = datestr(dates(1),31);
    xlabel(sprintf('month.day.HH - %s start',monstr))
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

