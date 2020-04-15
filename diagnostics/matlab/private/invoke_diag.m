% This function is called by "results.m" script and is not intended to be called 
% by anyone else. 

function invoke_diag(dart_output_file, style, metric, var, lev)

    if strcmp(style, 'profile') == 1
        plot_rmse_xxx_profile_web(dart_output_file, metric, ...
            'obsname', var, 'verbose', false, 'MarkerSize', 6);
        
    else if strcmp(style, 'norm_profile') == 1
        plot_rmse_xxx_norm_profile_web(dart_output_file, metric, ...
            'obsname', var, 'verbose', false, 'MarkerSize', 6);

    else
        plot_rmse_xxx_evolution_web(dart_output_file, metric, ...
            'obsname', var, 'level', lev, 'verbose', false, 'MarkerSize', 6);
    
    end
    snapnow
    close all
end
