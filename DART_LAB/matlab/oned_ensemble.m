function varargout = oned_ensemble(varargin)
% ONED_ENSEMBLE explore the details of ensemble data assimilation for a scalar.
%
%      Click on the 'Create New Ensemble' button to activate the interactive 
%      observation generation mechanism and lay down a set of 'observations'
%      representative of your ensemble. (Think: Some H() operator has 
%      converted the model state to an expected observation.)
%
%      After you have an ensemble and an observation, choose an assimilation 
%      algorithm and click 'Update Ensemble'. The algorithm is applied and the 
%      Posterior (blue) is plotted below the Prior (green).  Choose 'EAKF' 
%      and click 'Update' ...  multiple times. Do the same for 'EnKF' and
%      'RHF'. Hmnnnn ....
%
%      change the Observation Error SD, lay down an ensemble pretty far away
%      from the observation - have fun with it.
%
% See also: gaussian_product, oned_model, twod_ensemble, run_lorenz_63, 
%           run_lorenz_96

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2009, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

% Last Modified by GUIDE v2.5 25-Mar-2009 14:34:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @oned_ensemble_OpeningFcn, ...
                   'gui_OutputFcn',  @oned_ensemble_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before oned_ensemble is made visible.
function oned_ensemble_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to oned_ensemble (see VARARGIN)

help oned_ensemble

% Choose default command line output for oned_ensemble
handles.output = hObject;

% Insert the ensemble structure into this
handles.ens_size = 0;
handles.ens_members = 0;
handles.h_obs_plot = 0;
handles.h_update_ens = 0;
handles.h_ens_member = 0;
handles.h_obs_ast = 0;

% Update handles structure
guidata(hObject, handles);

% Go ahead and plot the initial observational error distribution
h_observation = get(handles.edit1);
h_obs_error_sd = get(handles.edit2);
observation = str2double(h_observation.String);
obs_error_sd = str2double(h_obs_error_sd.String);
handles.h_obs_plot = plot_gaussian(observation, obs_error_sd, 1);
set(handles.h_obs_plot, 'Color', 'r', 'Linestyle', '--', 'Linewidth', 2);
hold on

% Plot an asterisk 
handles.h_obs_ast = plot(observation, 0, 'r*', 'MarkerSize', 16);

% Set a basic plotting domain range that includes mean +/- 3 obs SDs
lower = observation - 3*obs_error_sd;
upper = observation + 3*obs_error_sd;
axis([lower upper -0.2 1]);


set(gca, 'YTick', [0 0.2 0.4 0.6 0.8]);
set(gca, 'YTickLabel', [0 0.2 0.4 0.6 0.8]);

hold on
plot([lower upper], [0 0], 'k', 'Linewidth', 2);

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes oned_ensemble wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%----------------------------------------------------------------------


% --- Outputs from this function are returned to the command line.
function varargout = oned_ensemble_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%----------------------------------------------------------------------


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disable the update ensemble button and all other active buttons
set(handles.pushbutton1, 'Enable', 'Off');
set(handles.pushbutton2, 'Enable', 'Off');
set(handles.edit1, 'Enable', 'Off');
set(handles.edit2, 'Enable', 'Off');

% Clear out any old ensemble members if they exist
for i = 1:handles.ens_size
   set(handles.h_ens_member(i), 'Visible', 'off');
end

% Remove mean and sd of old ensemble
set(handles.text2, 'String', 'Prior Mean = ');
set(handles.text3, 'String', 'Prior SD = ');

hold on
% Set a basic plotting domain range that includes mean +/- 3 obs SDs
h_observation = get(handles.edit1);
h_obs_error_sd = get(handles.edit2);
observation = str2double(h_observation.String);
obs_error_sd = str2double(h_obs_error_sd.String);
lower = observation - 3*obs_error_sd;
upper = observation + 3*obs_error_sd;
axis([lower upper -0.2 1]);

set(gca, 'YTick', [0 0.2 0.4 0.6 0.8]);
set(gca, 'YTickLabel', [0 0.2 0.4 0.6 0.8]);


% Messages should start 1/10 of the way across the screen
x_message = lower + 0.1 * (upper - lower);
h_click = text(x_message, 0.7, 'Click on x-axis to create member', 'FontSize', 16);

% Need to guarantee at least 2 ensemble members
ens_size = 0;

h_err_text = text(x_message, 0.9, 'Click inside graphics box to select member', ...
   'Color', 'r', 'FontSize', 16, 'Visible', 'off');

while ens_size < 2
   [xt, zt] = ginput(1);

   % Clear out the error message if it's been made visible
   set(h_err_text, 'Visible', 'off');

   if(xt < lower | xt > upper) 
      set(h_err_text, 'Visible', 'on');
   else
      ens_size = ens_size + 1;
      x(ens_size) = xt;
      y(ens_size) = 0;
      handles.h_ens_member(ens_size) = ...
         plot(x(ens_size), y(ens_size), '*', 'MarkerSize', 16, 'Color', [0 0.73 0]);

      % Display the prior mean and sd
      prior_mean = mean(x);
      prior_sd = std(x);
      set(handles.text2, 'String', ['Prior Mean = ', num2str(prior_mean)]);
      set(handles.text3, 'String', ['Prior SD = ', num2str(prior_sd)]);

   end
end


h_finish = text(x_message, 0.5, 'Click outside of plot to finish', 'Fontsize', 16);

while ens_size < 1000
   [xt, zt] = ginput(1);
   % Terminate by clicking outside of graph range
   if(xt > upper | xt < lower)
      break;
   else
      ens_size = ens_size + 1;
      x(ens_size) = xt;
      y(ens_size) = 0;
      handles.h_ens_member(ens_size) = ...
         plot(x(ens_size), y(ens_size), '*', 'MarkerSize', 16, 'Color', [0 0.73 0]);

      % Display the prior mean and sd
      prior_mean = mean(x);
      prior_sd = std(x);
      set(handles.text2, 'String', ['Prior Mean = ', num2str(prior_mean)]);
      set(handles.text3, 'String', ['Prior SD = ', num2str(prior_sd)]);

   end

end

% Ensemble created, comupte mean and sd, clean up and return
% Set the global gui storage
handles.ens_size = ens_size;
handles.ens_members = x;

% Update handles structure
guidata(hObject, handles);


% Turn off the data entry messages
set(h_click, 'Visible', 'off');
set(h_finish, 'Visible', 'off');

% Enable the update ensemble button
set(handles.pushbutton1, 'Enable', 'On');
set(handles.pushbutton2, 'Enable', 'On');
set(handles.edit1, 'Enable', 'On');
set(handles.edit2, 'Enable', 'On');


%----------------------------------------------------------------------



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

% Enable things that an error might have turned off
set(handles.edit2, 'Enable', 'on')
set(handles.pushbutton1, 'Enable', 'on')

% Only enable the update ensemble pushbutton if an ensemble has been created
if(handles.ens_size > 0)
   set(handles.pushbutton2, 'Enable', 'on');
end

% Get the value of the observation
if(isfinite(str2double(get(hObject, 'String'))))
   observation = str2double(get(hObject, 'String'));
else
   set(handles.edit1, 'String', '???');

   % Disable other input to guarantee only one error at a time!
   set(handles.edit2, 'Enable', 'off')
   set(handles.pushbutton1, 'Enable', 'off')
   set(handles.pushbutton2, 'Enable', 'off')
   return
end

% Get the value of the observation error sd
h_obs_error_sd = get(handles.edit2);
obs_error_sd = str2double(h_obs_error_sd.String);

% Plot the updated distribution
set(handles.h_obs_plot, 'Visible', 'Off');
handles.h_obs_plot = plot_gaussian(observation, obs_error_sd, 1);
set(handles.h_obs_plot, 'Color', 'r', 'Linestyle', '--', 'Linewidth', 2);

% Move the observation asterisk
set(handles.h_obs_ast, 'Visible', 'Off');
handles.h_obs_ast = plot(observation, 0, 'r*', 'MarkerSize', 16);

% Set a basic plotting domain range that includes mean +/- 3 obs SDs
h_observation = get(handles.edit1);
h_obs_error_sd = get(handles.edit2);
observation = str2double(h_observation.String);
obs_error_sd = str2double(h_obs_error_sd.String);
lower = observation - 3*obs_error_sd;
upper = observation + 3*obs_error_sd;
axis([lower upper -0.2 1]);

set(gca, 'YTick', [0 0.2 0.4 0.6 0.8]);
set(gca, 'YTickLabel', [0 0.2 0.4 0.6 0.8]);

hold on
plot([lower upper], [0 0], 'k', 'Linewidth', 2);

% Update handles structure
guidata(hObject, handles);


%----------------------------------------------------------------------


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%----------------------------------------------------------------------


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

% Enable things that an error might have turned off
set(handles.edit1, 'Enable', 'on')
set(handles.pushbutton1, 'Enable', 'on')

% Only enable the update ensemble pushbutton if an ensemble has been created
if(handles.ens_size > 0)
   set(handles.pushbutton2, 'Enable', 'on');
end

% Get the value of the observation
if(isfinite(str2double(get(hObject, 'String'))) && ...
   str2double(get(hObject, 'String')) > 0)
   obs_error_sd = str2double(get(hObject, 'String'));
else
   set(handles.edit2, 'String', '???');

   % Disable other input to guarantee only one error at a time!
   set(handles.edit1, 'Enable', 'off')
   set(handles.pushbutton1, 'Enable', 'off')
   set(handles.pushbutton2, 'Enable', 'off')
   return
end

% Get the value of the observation
h_observation = get(handles.edit1);
observation = str2double(h_observation.String);

% Plot the updated distribution
set(handles.h_obs_plot, 'Visible', 'off');
handles.h_obs_plot = plot_gaussian(observation, obs_error_sd, 1);
set(handles.h_obs_plot, 'Color', 'r', 'Linestyle', '--', 'Linewidth', 2);

% Set a basic plotting domain range that includes mean +/- 3 obs SDs
h_observation = get(handles.edit1);
h_obs_error_sd = get(handles.edit2);
observation = str2double(h_observation.String);
obs_error_sd = str2double(h_obs_error_sd.String);
lower = observation - 3*obs_error_sd;
upper = observation + 3*obs_error_sd;
axis([lower upper -0.2 1]);

set(handles.h_obs_plot, 'Color', 'r', 'Linestyle', '--', 'Linewidth', 2);

set(gca, 'YTick', [0 0.2 0.4 0.6 0.8]);
set(gca, 'YTickLabel', [0 0.2 0.4 0.6 0.8]);

hold on
plot([lower upper], [0 0], 'k', 'Linewidth', 2);

% Update handles structure
guidata(hObject, handles);


%----------------------------------------------------------------------


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%----------------------------------------------------------------------


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


%----------------------------------------------------------------------


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%----------------------------------------------------------------------


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Turn off any old points
set(handles.h_update_ens, 'Visible', 'off');

ensemble = handles.ens_members;
h_observation = get(handles.edit1);
h_obs_error_sd = get(handles.edit2);
observation = str2double(h_observation.String);
obs_error_sd = str2double(h_obs_error_sd.String);


% Figure out which filter option is currently selected
h_filter_kind = get(handles.popupmenu1);

filter_type = char(h_filter_kind.String(h_filter_kind.Value));

switch filter_type
   case 'EAKF'
      [obs_increments, err] = ...
         obs_increment_eakf(ensemble, observation, obs_error_sd^2);
   case 'EnKF'
      [obs_increments, err] = ...
         obs_increment_enkf(ensemble, observation, obs_error_sd^2);
   case 'RHF'
      [obs_increments, err] = ...
         obs_increment_rhf(ensemble, observation, obs_error_sd^2);
end

% Add on increments to get new ensemble
new_ensemble = ensemble + obs_increments;

y(1:size(ensemble)) = -0.1;
handles.h_update_ens = plot(new_ensemble, y, '*', 'MarkerSize', 16, 'Color', 'Blue');
guidata(hObject, handles)
