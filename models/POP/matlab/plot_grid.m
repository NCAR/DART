function plot_grid(fname)
%% plot_grid ... plots the ULAT,ULON and TLAT,TLON variables from a netcdf file.
% 
% fname = 'h.A1.10.nc';
% plot_grid(fname)
%

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

ulat = nc_varget(fname,'ULAT') * 180/pi;;
ulon = nc_varget(fname,'ULON') * 180/pi;;
tlat = nc_varget(fname,'TLAT') * 180/pi;;
tlon = nc_varget(fname,'TLON') * 180/pi;;

Ulat = ulat(:);
Ulon = ulon(:);
Tlat = tlat(:);
Tlon = tlon(:);

inds = find(tlon <= 0);
tlon(inds) = tlon(inds) + 360.0;
inds = find(ulon <= 0);
ulon(inds) = ulon(inds) + 360.0;

np = 5;
nx = size(ulat,1);
ny = size(ulat,2);

figure(1); clf; orient landscape
   plot(Ulon,Ulat,'ko',Tlon,Tlat,'rx')
   legend('U,V Grid','S,T Grid')
   legend boxoff


figure(2); clf; orient landscape

   i1 =  1;
   iN = np;
   j1 = 1;
   jN = np;

   myplot(i1,iN,j1,jN,ulon,ulat,tlon,tlat)

figure(3); clf; orient landscape

   i1 =  1;
   iN = np;
   j1 = ny-np+1;
   jN = ny;

   myplot(i1,iN,j1,jN,ulon,ulat,tlon,tlat)

figure(4); clf; orient landscape

   i1 = nx-np+1;
   iN = nx;
   j1 = ny-np+1;
   jN = ny;

   myplot(i1,iN,j1,jN,ulon,ulat,tlon,tlat)

figure(5); clf; orient landscape

   i1 = nx-np+1;
   iN = nx;
   j1 = 1;
   jN = np;

   myplot(i1,iN,j1,jN,ulon,ulat,tlon,tlat)


function myplot(i1,iN,j1,jN,ulon,ulat,tlon,tlat)

   h = plot(ulon(i1:iN,j1:jN), ulat(i1:iN,j1:jN), 'ko', ...
            tlon(i1:iN,j1:jN), tlat(i1:iN,j1:jN), 'rx');
   set(h,'MarkerSize',1)

   for i = i1:iN
   for j = j1:jN
   %  h = text(ulon(i,j),ulat(i,j), sprintf('%d,%d',i,j));
      str   = sprintf('''%d,%d''',i,j);
      mystr = sprintf('h1 = text(%f,%f,%s);',ulon(i,j),ulat(i,j),str);
      eval(mystr)
      set(h1,'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontSize',10,'Color','k')

      mystr = sprintf('h2 = text(%f,%f,%s);',tlon(i,j),tlat(i,j),str);
      eval(mystr)
      set(h2,'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontSize',10,'Color','r')
   end
   end
   title('POP grid layout')
   xlabel('U,V grid in black; T,S grid in red')
   set(gca,'box','off')

%  set(h,'Visible','off'); % Get the domain right.
%   h = [h1 h2];
%  legend('U,V Grid','S,T Grid')
%  legend boxoff

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
