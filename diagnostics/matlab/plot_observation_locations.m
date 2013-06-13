%% PLOT_OBSERVATION_LOCATIONS : Plots the locations of the input observations
%
% By default this command creates 2d plots of observation locations,
% one per time epoch, from data output from the obs_diag program if
% the 'print_obs_locations' namelist item in the &obs_diag list is .true.
%
% WARNING: obs_diag no longer creates the observation_locations.NNN.dat file.
% WARNING: Convert your observation sequence files to netCDF using 'obs_seq_to_netcdf'
% WARNING: and then use 'plot_obs_netcdf.m'
%
% There are lots of user settable options.  This script prompts you
% interactively for the most common ones.  Then it calls PlotObsLocs()
% with the proper argument list to pass in your selections.
%
% See the documentation for PlotObsLocs() -- it has a lot of arguments in the
% calling sequence.

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% setup all args to be the string 'default', which will be interpreted by 
% the PlotObsLocs routine to use the default values.   

disp(' ')
disp('WARNING: obs_diag no longer creates the observation_locations.NNN.dat file.')
disp('WARNING: Convert your observation sequence files to netCDF using ''obs_seq_to_netcdf''')
disp('WARNING: and then use ''plot_obs_netcdf.m''')
disp(' ')
 
ncfname     = 'default';
plotd       = 'default';
dartqc      = 'default';
typelist    = 'default';
box         = 'default';
epochs      = 'default';
subset      = 'default';
world       = 'default';
writeplot   = 'default';
loc2dstring = 'default';
loc3dstring = 'default';
orientation = 'default';
viewlist    = 'default';
invertz     = 'default';
plotname    = 'default';

% what the arg list looks like:
%PlotObsLocs(in_dartqc, in_box, in_typelist, in_epochlist, in_subset, in_plotd, 
%             in_world, in_invertz, in_writeplot, in_legend2dloc, in_legend_3dloc, 
%             in_viewlist, in_ncfname, in_orientation, in_plotname)

done = 0;
disp('Plot observations at their proper locations.  Many subsetting options exist.');
disp('Hitting <cr> to answer the questions will use the default value,');
disp('or - once you have made a selection - reuse the previous value.');
disp(' '); 
disp('The default plotting options are:');
disp('  2D plot, full world map, all obs types, all times, ');
disp('  no file output, Z axis increases up.');
disp(' '); 

   % What file has the metadata 
   reply = input('Enter the netCDF file name with the metadata (default: ''obs_diag_output.nc''):  ');
   if (~isempty(reply))
      ncfname = reply;
   end

% loop and keep the previous default until the user says to quit
while done == 0

   % 2D or 3D plot?
   reply = input('Input 2 for 2D plot, 3 for 3D plot (default: 2D):  ');
   if (~isempty(reply))
      plotd = reply;
   end

    
   % plot selecting on dart QC value
   disp('')
   disp('DART QC Values ... 0 == all OK')
   disp('DART QC Values ... 1 == Evaluated only')
   disp('DART QC Values ... 2 == OK but posterior forward operator failed')
   disp('DART QC Values ... 3 == Evaluated only, BUT posterior forward operator failed')
   disp('DART QC Values ... 4 == prior forward operator failed')
   disp('DART QC Values ... 5 == not used because of namelist control')
   disp('DART QC Values ... 6 == prior qc rejected')
   disp('DART QC Values ... 7 == outlier rejected')
   disp('   a negative value means everything ''up to'' that value, i.e.')
   disp('    -3 == 0, 1, 2, and 3          -99 == everything');
   reply = input('Input DART QC val (default: 0):  ');
   if (~isempty(reply))
      dartqc = reply;
   end
   
   % restrict observations to a particular observation type?
   disp('')
   reply = input('Input [obs type list] to plot only some obs types, ''default'' to reset:  ');
   if (~isempty(reply))
      typelist = reply;
   end
    
   % restrict observations to a particular subregion?
   reply = input('Input [xmin xmax ymin ymax] for bounding box, ''default'' to reset:  ');
   if (~isempty(reply))
      box = reply;
   end
    
   % restrict input to particular time epochs?
   reply = input('Input [epoch list] for particular times, ''default'' to reset:  ');
   if (~isempty(reply))
      epochs = reply;
   end
    
   % sample data to reduce counts?
   reply = input('Input count for random subset of each obs type, ''default'' to reset:  ');
   if (~isempty(reply))
      subset = reply;
   end
    
   % plot world map beneath?
   reply = input('Input 0 to remove world map, 1 to restore it (default: 1):  ');
   if (~isempty(reply))
      world = reply;
   end
    
   % write out .ps files?
   reply = input('Input 1 to write .ps files for each plot, 0 to reset (default: 0):  ');
   if (~isempty(reply))
      writeplot = reply;
   end
    
   % output plot filename
   if (writeplot == 1)
      reply = input('Set output plot filename, ''default'' to reset:  ');
      if (~isempty(reply))
         plotname = reply;
      end
   end
    
   % legendloc
   disp('Examples of valid legend strings:');
   disp(' ''North'' : inside plot box near top');
   disp(' ''South'',''East'',''West'',''NorthEast'', etc');
   disp(' ''NorthOutside'' : outside plot box near top');
   disp(' ''SouthOutside'', ''NorthEastOutside'', etc');
   disp(' ''Best'' : least conflict with data in plot');
   disp(' ''BestOutside'' : least unused space outside plot');
   disp(' defaults to: ''SouthWest'' for 2D plots; ''NorthWest'' for 3D plots');
   reply = input('Input Matlab string for legend location, ''default'' to reset:  ');
   if (~isempty(reply))
      if (plotd == 3)
         loc3dstring = reply;
      else
         loc2dstring = reply;
      end
   end
    
   % page orientation
   disp('Possible page orientations: ''landscape'', ''portrait'', ''tall'' ');
   reply = input('Input Matlab string for page orientation, ''default'' to reset (default: ''landscape''):  ');
   if (~isempty(reply))
      orientation = reply;
   end
    
   % viewlist and invert z axis
   invertz = 1;
   if (plotd == 3) 
      reply = input('Input [az el; az el] list for 3d views, ''default'' to reset:  ');
      if (~isempty(reply))
         viewlist = reply;
      end
   
      reply = input('Input 1 to invert Z axis (e.g. for pressure); 0 otherwise:  ');
      if (~isempty(reply))
         invertz = reply;
      end
   end
    
   PlotObsLocs(dartqc, box, typelist, epochs, subset, plotd, world, invertz, writeplot, ...
               loc2dstring, loc3dstring, viewlist, ncfname, orientation, plotname);
   
   reply = input('Quit now or plot again? 1 = quit, <cr> plots again:  ');
   if (~isempty(reply))
      done = reply;
   end
  
end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

