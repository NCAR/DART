function clim = add_gray_color
%% Augment the colormap and the CLim so that the lowest color index can be
% forced to a light gray without compromising the data range.
% Some figures plot the QC value. cmin may equal cmax, causing problems.

bob      = colormap;
ncolors  = length(bob);
bob(1,:) = 0.7; % set lowest color to be light gray.
colormap(bob);

clim = get(gca,'CLim');
cmin = clim(1);
cmax = clim(2);

if (cmin ~= cmax)
    dz      = linspace(double(cmin),double(cmax),ncolors-1); % desired dynamic range
    dclim   = dz(2) - dz(1);
    newcmin = double(cmin) - dclim;
    clim    = [newcmin cmax]; % add extra bin at bottom that no data will use.
else
    disp('limits are identical')
    newcmin = double(cmin);  % newcmin must be same type as 'elev'.
    clim    = [double(newcmin-1.5) double(newcmin+1.5)];
end
