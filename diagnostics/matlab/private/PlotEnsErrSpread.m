function PlotEnsErrSpread( pinfo )
%% PlotEnsErrSpread     Creates summary plots of error and spread
%
% PlotEnsErrSpread is intended to be called by 'plot_ens_err_spread'.
% The only input argument is a structure with model-dependent
% components.
%
% USAGE: EnsErrSpread( pinfo )
%
% STRUCTURE COMPONENTS FOR low-order models
% truth_file      name of netCDF DART file with copy tagged 'true state'
% diagn_file      name of netCDF DART file with copies tagged 'ensemble mean'
%                 and 'ensemble spread'
% var             name of netCDF variable of interest
% var_inds        indices of variables of interest
%
% Example 0   (9var  model)
%%--------------------------------------------------------
% pinfo.truth_file = 'true_state.nc';
% pinfo.diagn_file = 'preassim.nc';
% pinfo.var        = 'state';
% pinfo.var_inds   = [ 1 2 3 4 5 6 7 8 9 ];
% PlotEnsErrSpread(pinfo)
%
% Example 1   (Lorenz_96  model)
%%--------------------------------------------------------
% pinfo.truth_file = 'true_state.nc';
% pinfo.diagn_file = 'preassim.nc';
% pinfo.var        = 'state';
% pinfo.var_inds   = [ 3 4 36 39 22 ];
% PlotEnsErrSpread(pinfo)
%
% Example 2   (Lorenz_96_2scale  model)
%%--------------------------------------------------------
% pinfo.truth_file = 'true_state.nc';
% pinfo.diagn_file = 'preassim.nc';
% pinfo.var        = 'X';
% pinfo.var_inds   = [ 3 18 27 ];
% PlotEnsErrSpread(pinfo)
%
% Example 3 (FMS BGrid model)
%%--------------------------------------------------------
% pinfo.truth_file = 'true_state.nc';
% pinfo.diagn_file = 'preassim.nc';
% pinfo.var        = 'u';
% pinfo.level      = 3;
% pinfo.latitude   = 23.5;
% pinfo.longitude  = 45.67;
% PlotEnsErrSpread(pinfo)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

switch lower(pinfo.model)
    
    case '9var'
        
        truth      = get_hyperslab('fname',pinfo.truth_file, ...
            'varname','state', ...
            'permute','T', ...
            'tindex1',pinfo.truth_time(1), ...
            'tcount',pinfo.truth_time(2)) ;
        
        ens_mean   = get_hyperslab('fname',pinfo.diagn_file, ...
            'varname','state_mean', ...
            'permute','T', ...
            'tindex1',pinfo.diagn_time(1), ...
            'tcount',pinfo.diagn_time(2)) ;
        
        ens_spread = get_hyperslab('fname',pinfo.diagn_file, ...
            'varname','state_sd', ...
            'permute','T', ...
            'tindex1',pinfo.diagn_time(1), ...
            'tcount',pinfo.diagn_time(2)) ;
        
        % Use three different figures with three subplots each
        for i = 1:3
            figure(i); clf
            for j = 1:3
                
                ivar = (i - 1)*3 + j;
                
                err         = total_err(ens_mean(:,ivar) , truth(:,ivar));
                errTotal    = sum(err)                / pinfo.time_series_length;
                spreadTotal = sum(ens_spread(:,ivar)) / pinfo.time_series_length;
                string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
                string2 = ['time-mean Ensemble Spread = ' num2str(spreadTotal)];
                
                fprintf('%s model Variable %d\n',pinfo.model,ivar)
                
                subplot(3, 1, j);
                plot(pinfo.time,err, 'b', ...
                    pinfo.time,ens_spread(:, ivar), 'r');
                s1 = sprintf('%s model Var %d Ensemble Error Spread', pinfo.model, ivar);
                title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')
                legend(string1,string2,'Location','NorthEast')
                legend boxoff
                xlabel(sprintf('model "days" (%d timesteps)',pinfo.time_series_length))
                ylabel('distance')
            end
        end
        
    case {'lorenz_63','lorenz_84','lorenz_96', 'lorenz_96_2scale', ...
          'lorenz_04','forced_lorenz_96','ikeda','simple_advection','lorenz_96_tracer_advection','null'}
        
        truth      = get_hyperslab('fname',pinfo.truth_file, ...
            'varname',pinfo.var, ...
            'permute','T', ...
            'tindex1',pinfo.truth_time(1), ...
            'tcount',pinfo.truth_time(2)) ;
        
        varname = sprintf('%s_mean',pinfo.var);
        ens_mean   = get_hyperslab('fname',pinfo.diagn_file, ...
            'varname',varname, ...
            'permute','T', ...
            'tindex1',pinfo.diagn_time(1), ...
            'tcount',pinfo.diagn_time(2)) ;
        
        varname = sprintf('%s_sd',pinfo.var);
        ens_spread = get_hyperslab('fname',pinfo.diagn_file, ...
            'varname',varname, ...
            'permute','T', ...
            'tindex1',pinfo.diagn_time(1), ...
            'tcount',pinfo.diagn_time(2)) ;
        
        clf; iplot = 0;
        for ivar = pinfo.var_inds,
            iplot = iplot + 1;
            err         = total_err(ens_mean(:,ivar) , truth(:,ivar));
            errTotal    = sum(err)                / pinfo.time_series_length;
            spreadTotal = sum(ens_spread(:,ivar)) / pinfo.time_series_length;
            string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
            string2 = ['time-mean Ensemble Spread = ' num2str(spreadTotal)];
            
            subplot(length(pinfo.var_inds), 1, iplot);
            plot(pinfo.time,err, 'b', ...
                pinfo.time,ens_spread(:,ivar), 'r');
            s1 = sprintf('%s model Var %d Ensemble Error Spread', pinfo.model, ivar);
            title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')
            legend(string1,string2,'Location','NorthEast')
            legend boxoff
            xlabel(sprintf('model "days" (%d timesteps)',pinfo.time_series_length))
            ylabel('distance')
        end
        
    case {'fms_bgrid','pe2lyr','mitgcm_ocean','cam','wrf','sqg','pop'}
        
        clf;
        
        truth      = get_hyperslab('fname',pinfo.truth_file, ...
            'varname', pinfo.var, ...
            'permute', 'T', ...
            'levelindex', pinfo.levelindex, ...
            'latindex', pinfo.latindex, ...
            'lonindex', pinfo.lonindex, ...
            'tindex1', pinfo.truth_time(1), ...
            'tcount',pinfo.truth_time(2)) ;
        
        varname = sprintf('%s_mean',pinfo.var);
        ens_mean   = get_hyperslab('fname',pinfo.diagn_file, ...
            'varname', varname, ...
            'permute', 'T', ...
            'levelindex', pinfo.levelindex, ...
            'latindex', pinfo.latindex, ...
            'lonindex', pinfo.lonindex, ...
            'tindex1', pinfo.diagn_time(1), ...
            'tcount',pinfo.diagn_time(2)) ;
        
        varname = sprintf('%s_sd',pinfo.var);
        ens_spread = get_hyperslab('fname',pinfo.diagn_file, ...
            'varname', varname, ...
            'permute', 'T', ...
            'levelindex', pinfo.levelindex, ...
            'latindex', pinfo.latindex, ...
            'lonindex', pinfo.lonindex, ...
            'tindex1', pinfo.diagn_time(1), ...
            'tcount',pinfo.diagn_time(2)) ;
        
        subplot(2,1,1)
        PlotLocator(pinfo);
        
        subplot(2,1,2)
        err         = total_err(ens_mean, truth);
        errTotal    = sum(err)        / pinfo.time_series_length;
        spreadTotal = sum(ens_spread) / pinfo.time_series_length;
        string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
        string2 = ['time-mean Ensemble Spread = ' num2str(spreadTotal)];
        
        plot(pinfo.time,err, 'b', pinfo.time,ens_spread, 'r');
        
        s1 = sprintf('Ensemble Mean Error, Ensemble Spread %s ''%s''',pinfo.model,pinfo.var);
        s2 = sprintf('level %d lat %.2f lon %.2f', ...
            pinfo.level, pinfo.latitude, pinfo.longitude);
        title({s1, s2, pinfo.diagn_file},'interpreter','none','fontweight','bold');
        
        legend(string1,string2,'Location','NorthEast');
        legend boxoff
        xdates(pinfo.time);
        ylabel('distance');
        
    case {'mpas_atm'}
        error('not supported yet')
        clf;
        
        truth      = get_hyperslab('fname',pinfo.truth_file, ...
            'varname', pinfo.var, ...
            'permute', 'T', ...
            'levelindex', pinfo.levelindex, ...
            'cellindex', pinfo.cellindex, ...
            'tindex1', pinfo.truth_time(1), ...
            'tcount',pinfo.truth_time(2)) ;
        
        varname    = sprintf('%s_mean',pinfo.var);
        ens_mean   = get_hyperslab('fname',pinfo.diagn_file, ...
            'varname', varname, ...
            'permute', 'T', ...
            'levelindex', pinfo.levelindex, ...
            'cellindex', pinfo.cellindex, ...
            'tindex1', pinfo.diagn_time(1), ...
            'tcount',pinfo.diagn_time(2)) ;
        
        varname    = sprintf('%s_sd',pinfo.var);
        ens_spread = get_hyperslab('fname',pinfo.diagn_file, ...
            'varname', varname, ...
            'permute', 'T', ...
            'levelindex', pinfo.levelindex, ...
            'cellindex', pinfo.cellindex, ...
            'tindex1', pinfo.diagn_time(1), ...
            'tcount',pinfo.diagn_time(2)) ;
        
        subplot(2,1,1)
        PlotLocator(pinfo);
        
        subplot(2,1,2)
        err         = total_err(ens_mean, truth);
        errTotal    = sum(err)        / pinfo.time_series_length;
        spreadTotal = sum(ens_spread) / pinfo.time_series_length;
        string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
        string2 = ['time-mean Ensemble Spread = ' num2str(spreadTotal)];
        
        plot(pinfo.time,err, 'b', pinfo.time,ens_spread, 'r');
        
        s1 = sprintf('Ensemble Mean Error, Ensemble Spread %s ''%s''',pinfo.model,pinfo.var);
        s2 = sprintf('level number %d lat %.2f lon %.2f', ...
            pinfo.level, pinfo.latCell(pinfo.cellindex), pinfo.lonCell(pinfo.cellindex));
        title({s1, s2, pinfo.diagn_file},'interpreter','none','fontweight','bold');
        
        legend(string1,string2,'Location','NorthEast');
        legend boxoff
        xdates(pinfo.time);
        ylabel('distance');
        
    otherwise
        error('model %s unknown.',pinfo.model)
end

%======================================================================
% Subfunctions
%======================================================================


function PlotLocator(pinfo)
plot(pinfo.longitude,pinfo.latitude,'pb','MarkerSize',12,'MarkerFaceColor','b');
axlims = axis;
axlims = axlims + [-20 20 -20 20];
grid on
axis image
axis(axlims)
if (axlims(2) < 0)
    continents('hollow','dateline');
else
    continents('hollow','greenwich');
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
