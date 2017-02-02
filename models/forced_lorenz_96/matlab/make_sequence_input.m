function make_sequence_input(varargin)
%% make_sequence_input creates input for create_obs_sequence
%  
% Just as an example, this will create the text input required
% to create N evenly-spaced observation prototypes. There are no
% observation values in these. The output file name is 'input_to_create_obs_sequence.txt'
% create_obs_sequence can easily be run by hand for a small number
% of observations, but it gets tedious for large numbers. This
% script just saves on some typing - nothing more.
% 
% num_obs       is the number of observations on the unit circle
% obs_err_var   is the observation error variance
% loc_offset    is a static offset for each observation location.
%               loc_offset = 0.0 will make the initial observation
%               coincident with the model grid.
%
% Three Examples:
% make_sequence_input('num_obs',40)
% make_sequence_input('num_obs',40,'obs_err_var',1.0)
% make_sequence_input('num_obs',40,'obs_err_var',1.0,'loc_offset',0.02)
%
% Example on how to use the output file ...
%
% [unixprompt] ./create_obs_sequence < input_to_create_obs_sequence.txt
%
% The output file from create_obs_sequence will be called 'set_def.out'
% and can be used as input to create_fixed_network_sequence. 

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%%------------------------------------------------------------------------------
% This bit just parses the input and uses default values if need be.

default_num_obs     = 40;
default_obs_err_var = 1.0;
default_loc_offset  = 0.0;
default_graphic     = 0;
p = inputParser;

addParamValue(p,'num_obs',    default_num_obs,    @isnumeric);
addParamValue(p,'obs_err_var',default_obs_err_var,@isnumeric);
addParamValue(p,'loc_offset', default_loc_offset, @isnumeric);
addParamValue(p,'graphic',    default_graphic,    @isnumeric);
parse(p,varargin{:});

if ~isempty(fieldnames(p.Unmatched))
   disp('Extra inputs:')
   disp(p.Unmatched)
end

Nobs = p.Results.num_obs;
observation_error_variance = p.Results.obs_err_var;
loc_offset = p.Results.loc_offset;
graphic = p.Results.graphic;

locations = (0:Nobs-1)/Nobs + loc_offset;  % set the location to be uniform.
newlocs = mod(locations,1.0);
filename = 'input_to_create_obs_sequence.txt';

%%------------------------------------------------------------------------------
% Some stuff gets written to the file once at the beginning. 

fid = fopen(filename,'w');

fprintf(fid,'%d \n',Nobs); % upper bound on number of obs in sequence
fprintf(fid,'0 \n');          % number of copies of data
fprintf(fid,'0 \n');          % number of quality control values per field
fprintf(fid,'0 \n');          % 0 is a flag to continue 

%%------------------------------------------------------------------------------
% Some stuff gets written for every observation

for i = 1:Nobs

   fprintf(fid,'%d \n',1);             % observing RAW_STATE_VARIABLE
   fprintf(fid,'%f \n',newlocs(i));  % the unit sphere location 
   fprintf(fid,'0 0 \n');              % days and seconds of observation time
   fprintf(fid,'%f \n',observation_error_variance);

   if  i < Nobs
      fprintf(fid,'0 \n'); % 0 is a flag to continue 
   end

end

%%------------------------------------------------------------------------------
% Some stuff gets written to the file at the end. 

fprintf(fid,'set_def.out\n'); % name of observation sequence file prototype
fclose(fid);
fprintf('Created %s\n',filename)

%%------------------------------------------------------------------------------
% Optional visualizer of the observation locations.

if  graphic ~= 0 

   model_size = 40;
   gridlocs = 2*pi*(0:model_size)/model_size;
   x = ones(size(gridlocs));
   h1 = polar(gridlocs,x,'o');
   set(h1,'MarkerSize',10,'LineWidth',2)
   relabel(gca)
   
   hold on;
   y = ones(size(newlocs));
   h2 = polar(2*pi*newlocs,y,'r*');
   set(h2,'MarkerSize',20,'LineWidth',2)
   
   h = title('forced_lorenz_96');
   set(h,'Interpreter','none','FontSize',20,'FontWeight','bold');
   h = xlabel({'unit circle',filename});
   set(h,'Interpreter','none');
   
   h = legend('grid locations','observation locations');
   set(h,'Location','Best','FontSize',14);
   hold off;
   
end

end % main function


%% Helper function below ...

function relabel(handle)
    
% remove dotted radii and tickmarks and lables (actuall, all text)

t = findall(handle,'type','text');          %get handles for all text in polar plot
% set(t,'FontSize',20,'FontWeight','bold')
set(t,'FontSize',20)

for i = 1:numel(t)
   switch get(t(i),'String')
       case {'0'}
           set(t(i),'String',' 0.0', ...
                    'HorizontalAlignment','left')
       case {'90'}
           set(t(i),'String','0.25')
       case {'180'}
           set(t(i),'String','0.50', ...
                    'HorizontalAlignment','right')
       case {'270'}
           set(t(i),'String','0.75')
       otherwise
           set(t(i),'String',' ')
   end        
end

end % function relabel

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
