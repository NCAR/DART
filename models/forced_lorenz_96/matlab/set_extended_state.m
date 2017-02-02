function set_extended_state(fname, varargin)
%% set_extended_state modifies the extended state portion of the perfect_ics file.
%  
% The forced_lorenz_96 model has 40 state variables and 40 forcing variables.
% This function simply makes it easy to set the 40 forcing variables rather
% than edit the filter_ics file(s) by hand.
% 
% fname      is the name of the file to read AND WRITE.
% forcing    is the forcing for the extended state - may either be a scalar
%            or an array with exactly 40 elements.
%
% (bad) Example (produces a plot):
% fname = 'perfect_ics';
% forcing = randn(1,40) + 8.0;
% set_extended_state(fname,'forcing',forcing)
%
% Example - no graphic output:
% set_extended_state(fname,'forcing',forcing,'graphic',0)
%
% Example - unusual forcing:
% forcing = [ones(1,30)*6 ones(1,10)*12] + randn(1,40);
% set_extended_state(fname,'forcing',forcing)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%%------------------------------------------------------------------------------
% This bit just parses the input and uses default values if need be.

default_forcing = 40;
default_graphic = 1;
p = inputParser;

addRequired(p,'fname',@ischar);
addParamValue(p,'forcing',default_forcing,@isnumeric);
addParamValue(p,'graphic',default_graphic,@isnumeric);
parse(p, fname, varargin{:});

if ~isempty(fieldnames(p.Unmatched))
   disp('Extra inputs:')
   disp(p.Unmatched)
end

forcing = p.Results.forcing;
graphic = p.Results.graphic;

if (exist(fname,'file') ~= 2), error('file <%s> does not exist',fname), end

%%------------------------------------------------------------------------------
% Read the existing data
% Check to make sure the read was successful.

fid = fopen(fname,'r');
timeline = fgetl(fid);          % first record is the valid time of the model state
[modelstate, nlines] = fscanf(fid,'%f');  % vectorized read of the rest of the file.
fclose(fid);

if (nlines ~= 80)
   fprintf('Read %d elements for the model state.\n',nlines) 
   fprintf('Expected to read 80 elements.\n') 
   error('%s not the correct size for forced_lorenz_96.',fname)
end

%%------------------------------------------------------------------------------
% Check user input for sanity

if numel(forcing) == 1
   modelstate(41:80) = forcing;
elseif numel(forcing) == 40
   modelstate(41:80) = forcing;
else
   fprintf('forcing array can have exactly 1 or exactly 40 values.\n')
   error('forcing array had %d values.',numel(forcing))
end

%%------------------------------------------------------------------------------
% open the file
% write the timestamp
% write the model state. 

fid = fopen(fname,'w');
fprintf(fid,'%s\n',timeline);
fprintf(fid,'  %.14g\n',modelstate);
fclose(fid);

%%------------------------------------------------------------------------------
% Optional visualizer of the observation locations.

if  graphic ~= 0 

   xax = 1:40;

   figure; clf

   subplot(3,1,1)
   h1 = plot(xax,modelstate(1:40),'k*-');
   set(h1,'MarkerSize',10,'LineWidth',2)
   title('Original model state','FontSize',14)
   xlabel('state index')

   subplot(3,1,2)
   h2 = plot(xax+40,modelstate(41:80),'b*-');
   set(h2,'MarkerSize',20,'LineWidth',2)
   title('Forcing','FontSize',14);
   xlabel('"extended" state index')
   
   subplot(3,1,3)
   h3 = plot(modelstate,'b*-');
   set(h3,'MarkerSize',10,'LineWidth',1.5)
   title('Entire model state','FontSize',14);
   xlabel('state index')
   
end

end % main function

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
