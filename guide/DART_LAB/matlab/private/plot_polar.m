function h = plot_polar(y, x, mean_dist, string, model_size)
%% plot_polar

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% Y includes a wraparound point, x does not
x_t(model_size + 1) = x(1);
x_t(1:model_size) = x;
h = polar_dares(y, mean_dist + x_t, string);

end

