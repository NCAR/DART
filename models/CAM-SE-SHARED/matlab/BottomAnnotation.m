function BottomAnnotation(main)
% annotates the directory containing the data being plotted

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

subplot('position',[0.48 0.01 0.04 0.04])
axis off
fullname = which(main);   % Could be in MatlabPath
if( isempty(fullname) )
   if ( main(1) == '/' )  % must be a absolute pathname
      string1 = sprintf('data file: %s',main);
   else                   % must be a relative pathname
      mydir = pwd;
      string1 = sprintf('data file: %s/%s',mydir,main);
   end
else
   string1 = sprintf('data file: %s',fullname);
end

h = text(0.0, 0.5, string1);
set(h,'HorizontalAlignment','center', ...
      'VerticalAlignment','middle',...
      'Interpreter','none',...
      'FontSize',10)

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
