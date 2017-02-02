%% dbz_colors.m  Color scale from CIDD http://www.rap.ucar.edu/colorscales/dbz_40.colors

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%-15     -10     dark green
%-10     -6      dark olive green
%-6      -3      forest green
%-3      0       SpringGreen3
%0       3       medium sea green
%3       6       medium aquamarine
%6       9       medium slate blue
%9       12      blue
%12      15      blue4
%15      18      DarkOrchid4
%18      21      HotPink4
%21      24      maroon
%24      27      VioletRed4
%27      31      sienna
%31      35      chocolate
%35      40      goldenrod
%40      45      yellow
%45      50      dark salmon
%50      55      salmon
%55      60      firebrick2
%60      65      DeepPink1
%65      70      light grey
%70      80      snow

grey20            = [51,51,51];
dark_green        = [0,100, 0];
dark_olive_green  = [ 85,107, 47];
forest_green      = [ 34,139, 34];
SpringGreen3      = [ 0,205,102];
medium_sea_green  = [ 60,179,113];
medium_aquamarine = [102,205,170];
medium_slate_blue = [123,104,238];
blue              = [ 0, 0,255];
blue4             = [ 0, 0,139];
DarkOrchid4       = [104, 34,139];
HotPink4          = [139, 58, 98];
maroon            = [128, 0, 0];
VioletRed4        = [139, 71, 93];
sienna            = [160, 82, 45];
chocolate         = [210,105, 30];
goldenrod         = [218,165, 32];
yellow            = [255,255, 0];
dark_salmon       = [233,150,122];
salmon            = [250,128,114];
firebrick2        = [238, 44, 44];
DeepPink1         = [255, 20,147];
light_grey        = [211,211,211];
snow              = [255,250,250];

dbz_color = [grey20; ...                   % -18 --> -15 (no data)
             dark_green; dark_green; ...   % -15 --> -9
             dark_olive_green; ...         %  -9 --> -6
             forest_green; ...             %  -6 --> -3
             SpringGreen3; ...             %  -3 -->  0
	     medium_sea_green; ...         %   0 -->  3
             medium_aquamarine; ...        %   3 -->  6
             medium_slate_blue; ...        %   6 -->  9
	     blue; ...                     %   9 --> 12
             blue4; ...                    %  12 --> 15
             DarkOrchid4; ...              %  15 --> 18
             HotPink4; ...                 %  18 --> 21
             maroon; ...                   %  21 --> 24
             VioletRed4; ...               %  24 --> 27
   	     sienna; ...                   %  27 --> 30
             chocolate; chocolate; ...     %  30 --> 36
             goldenrod; ...                %  36 --> 39
             yellow; yellow; ...           %  39 --> 45
             dark_salmon; dark_salmon; ... %  45 --> 51
             salmon; ...                   %  51 --> 54
	     firebrick2; firebrick2; ...   %  54 --> 60
             DeepPink1; DeepPink1; ...     %  60 --> 66
             light_grey; ...               %  66 --> 69
             snow; snow; snow; snow];      %  69 --> 81

dbz_color = dbz_color/255;

colormap(dbz_color)

caxis([-18, 81]);


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
