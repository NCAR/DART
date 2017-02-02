   function ax=scalebar(xmin,ymin,wid,height,zmin,zmax)
%%  h = scalebar(xmin,ymin,wid,height,[zmin],[zmax])
%
%  Makes a scale of the size and positions specified. The values
%  zmin, zmax are annotated at the left(bottom) and right(top) ends
%  of the scalebar for reference.
%
%  If zmin,zmax are not supplied, the UserData attribute of the
%  current axis is used.
%
%  returns a handle ...
%
%  As of May 2009, scalebar now uses all the colors in the colormap;
%  previously only the first 64 colors were used.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

   if nargin <= 4,
      zmin = get(gca,'UserData');
      zmin = get(gca,'Clim');
      if length(zmin) == 0,
	 disp('WARNING: unknown data bounds, using [1,64]')
         zmax =  1;
         zmin = 64;
      else
         zmax = zmin(2);
         zmin = zmin(1);
      end
   end

   if strcmp(get(gcf,'NextPlot'),'replace'),
     set(gcf,'NextPlot','add')
   end

   ax = axes('position',[ xmin ymin wid height ]);
   frog = [1:length(colormap)];
   m = (zmax-zmin)/(length(colormap)-1);
   b = zmin - m;
   x_axis =  frog*m + b;

   if (wid > height)
      image(x_axis,[1 2],frog)
      axis([ zmin zmax 1 2 ])
      set(ax,'YTickLabel',[])
   else
      image([1 2],x_axis,frog'); set(ax,'YDir','normal');
      axis([ 1 2 zmin zmax ]);
      set(ax,'XTickLabel',[]);
      set(ax,'YAxisLocation','right')
   end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
