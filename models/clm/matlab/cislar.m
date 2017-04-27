function cislar(timeindex)
%%
%

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

priorfname = '/glade/user/thoar/CLM_leafc/preassim.nc';
postefname = '/glade/user/thoar/CLM_leafc/analysis.nc';
varname    = 'leafc';
levelindex = 1;
timeindex  = 6;
copystring = 'ensemble mean';
prior      = clm_get_var(priorfname,varname,copystring,levelindex,timeindex);
poste      = clm_get_var(postefname,varname,copystring,levelindex,timeindex);

innov        = prior;
innov.datmat = poste.datmat - prior.datmat;

% figure(1); clf
% clm_var_plot(prior);
% title('Prior')

% figure(3); clf
% clm_var_plot(innov);
% title('Innovations')

figure(2); clf
orient landscape
wysiwyg
y = myplot(poste);


function h3 = myplot(x)

   h1 = imagesc(x.lon, x.lat, x.datmat);
   set(h1,'AlphaData',~isnan(x.datmat))
   set(gca,'YDir','normal')
   set(gca,'FontSize',12,'FontWeight','bold')
   h2 = title('Posterior estimate of Leaf Carbon for 7 May 2000');
   set(h2,'Interpreter','none','FontSize',14,'FontWeight','bold')
   worldmap;
   h3 = colorbar;
   %set(h3,'FontSize',12,'FontWeight','bold')
   kids = get(h3,'YLabel');
   get(kids)
   set(kids,'String','Carbon (g/m^2)','Interpreter','TeX','FontSize',12)
   cbarpos = get(h3,'Position') 
   cbarpos(4) = cbarpos(4) - 0.1;
   cbarpos(2) = cbarpos(2) + 0.05 
   set(gca,'XTick',[0:60:360])
   set(gca,'YTick',[-90:30:90])
   hx = xlabel('Longitude (degrees East)');
   hy = ylabel('Latitude (degrees North)');
%  set(hx,'FontSize',14)
%  set(hy,'FontSize',14)
%  set(h3,'Position',cbarpos)

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
