% This function is called by "results.m" script and is not intended to be called 
% by anyone else. 

function invoke_diag(dart_output_file, style, metric, var, method, lev)

    if strcmp(style, 'profile') == 1
        plot_rmse_xxx_profile(dart_output_file, metric, ...
            'obsname', var, 'verbose', false, 'MarkerSize', 6, 'method', method);
        
    else
        plot_rmse_xxx_evolution(dart_output_file, metric, ...
            'obsname', var, 'level', lev, 'verbose', false, 'MarkerSize', 6, 'method', method);
    
    end
    snapnow
    close all
end
