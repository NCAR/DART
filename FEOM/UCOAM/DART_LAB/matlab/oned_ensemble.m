function varargout = oned_ensemble(varargin)
%% ONED_ENSEMBLE explore the details of ensemble data assimilation for a scalar.
%
%      Click on the 'Create New Ensemble' button to activate the interactive
%      observation generation mechanism and lay down a set of 'observations'
%      representative of your ensemble. (Think: Some H() operator has
%      converted the model state to an expected observation.) This is done by
%      placing the cursor near the axis in the plot and clicking. When you
%      have all the ensemble members you want, click in the grey area of
%      the window outside of the white axis plot.
%
%      After you have an ensemble and an observation, click 'Update Ensemble'.
%      The algorithm is applied and the Posterior (blue) is plotted below the
%      Prior (green). The mean and standard deviation of the posterior are
%      also printed on the plot.
%
%      The type of ensemble Kalman filter update can be chosen using the
%      pulldown menu at the bottom.
%
%      Checking the 'Show Inflation' box will also apply inflation to the
%      prior before doing the update and will print the mean and standard
%      deviation of the inflated prior and the resulting posterior. The
%      inflated prior and posterior are plotted on an axis below the
%      axis for the uninflated ensemble.
%
%      The 'EAKF' is a stochastic algorithm so repeated updates can be done
%      for the same prior and observation.
%
%      change the Observation Error SD, lay down an ensemble pretty far away
%      from the observation - have fun with it.
%
% See also: gaussian_product, oned_model, twod_ensemble, run_lorenz_63,
%           run_lorenz_96

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Last Modified by GUIDE v2.5 28-Aug-2009 16:29:57

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
function oned_ensemble_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to oned_ensemble (see VARARGIN)

help oned_ensemble

% set random number seed to same value to generate known sequences
rng('default')

% Choose default command line output for oned_ensemble
handles.output = hObject;

% Insert the ensemble structure into this
handles.ens_size         = 0;
handles.ens_members      = 0;
handles.h_obs_plot       = 0;
handles.h_update_ens     = 0;
handles.h_ens_member     = 0;
handles.h_obs_ast        = 0;
handles.h_update_lines   = 0;
handles.observation      = 0;
handles.obs_error_sd     = 0;
handles.inflation        = 1.5;
handles.plot_inflation   = false;
handles.h_inf_ens_member = 0;
handles.h_inf_up_ens     = 0;
handles.h_inf_lines      = 0;
handles.h_inf_axis       = 0;

% Update handles structure
guidata(hObject, handles);

% Get the initial observation, obs_error_sd and inflation from the gui
handles.observation  = str2double(get(handles.edit_observation,  'String'));
handles.obs_error_sd = str2double(get(handles.edit_obs_error_sd, 'String'));
handles.inflation    = str2double(get(handles.edit_inflation,    'String'));

% Go ahead and plot the initial observational error distribution
handles.h_obs_plot = plot_gaussian(handles.observation, handles.obs_error_sd, 1);
set(handles.h_obs_plot, 'Color', 'r', 'Linestyle', '--', 'Linewidth', 2);
hold on

% Plot an asterisk
handles.h_obs_ast = plot(handles.observation, 0, 'r*', 'MarkerSize', 16,'LineWidth',2.0);

% Set a basic plotting domain range that includes mean +/- 3 obs SDs
xlower = handles.observation - 3*handles.obs_error_sd;
xupper = handles.observation + 3*handles.obs_error_sd;
ylower = -0.4;
yupper = 1.0;
axis([xlower xupper ylower yupper]);

set(gca, 'YTick', [0 0.2 0.4 0.6 0.8]);
set(gca, 'YTickLabel', [0 0.2 0.4 0.6 0.8]);

hold on
plot([xlower xupper], [0 0], 'k', 'Linewidth', 2);

% Update handles structure
guidata(hObject, handles);

% Reset focus to the menu gui window
% Setting the axes clears the legend, gcbo restores focus
axes(handles.axes1);


% UIWAIT makes oned_ensemble wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%----------------------------------------------------------------------


% --- Outputs from this function are returned to the command line.
function varargout = oned_ensemble_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%----------------------------------------------------------------------


% --- Executes on button press in pushbutton_create_new.
function pushbutton_create_new_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_create_new (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disable the update ensemble button and all other active buttons
set(handles.pushbutton_update_ens, 'Enable', 'Off');
set(handles.edit_observation,      'Enable', 'Off');
set(handles.edit_obs_error_sd,     'Enable', 'Off');
set(handles.edit_inflation,        'Enable', 'Off');

% Clear out any old ensemble members if they exist
set(handles.h_ens_member,          'Visible', 'off');
set(handles.h_inf_ens_member,      'Visible', 'off');

set(handles.h_update_lines,        'Visible', 'off');
set(handles.h_inf_lines,           'Visible', 'off');
set(handles.h_inf_axis,            'Visible', 'off');

% Turn off any old update points
set(handles.h_update_ens,          'Visible', 'off');
set(handles.h_inf_up_ens,          'Visible', 'off');
set(handles.h_inf_ens_member,      'Visible', 'off');

clear_labels(handles);

hold on

% Set a basic plotting domain range that includes mean +/- 3 obs SDs
xlower = min(handles.observation - 3*handles.obs_error_sd, min(handles.ens_members));
xupper = max(handles.observation + 3*handles.obs_error_sd, max(handles.ens_members));
ylower = -0.4;
yupper = 1.0;
axis([xlower xupper ylower yupper]);

set(gca, 'YTick',      [0 0.2 0.4 0.6 0.8]);
set(gca, 'YTickLabel', [0 0.2 0.4 0.6 0.8]);

% Messages are centered in the middle.
xmid = (xupper + xlower) / 2.0;
h_click    = text(xmid,  0.8, {'Click inside graphics box to create member',...
    '(only X value is used)'}, 'FontSize', 16, 'HorizontalAlignment', 'center');

h_err_text = text(xmid, -0.15, 'An ensemble has to have at least 2 members.', ...
    'FontSize', 16, 'Visible', 'on', 'HorizontalAlignment', 'center','Color', 'r');

h_finish   = text(xmid, -0.15, 'Click outside of plot to finish', ...
    'Fontsize', 16, 'Visible', 'off', 'HorizontalAlignment', 'center');

ens_size = 0;
while ens_size < 100
    [xt, yt] = ginput(1);

    if(xt >= xlower && xt <= xupper && yt >= ylower && yt <= yupper)
        ens_size = ens_size + 1;
        x(ens_size) = xt;
        y(ens_size) = 0;
        handles.h_ens_member(ens_size) = ...
            plot(x(ens_size), y(ens_size), '*', 'MarkerSize', 16, 'Color', [0 0.73 0],'LineWidth',2.0);

        % Display the prior mean and sd
        prior_mean = mean(x);
        prior_sd = std(x);
        set(handles.text2, 'String', ['Prior Mean = ', num2str(prior_mean)]);
        set(handles.text3, 'String', ['Prior SD = ', num2str(prior_sd)]);
    elseif (ens_size < 2)
        set(h_err_text,'FontWeight','bold')
    else
        break;
    end

    % swap messages once you have a minimal ensemble.
    if (ens_size == 2)
        set(h_err_text, 'Visible', 'off');
        set(h_finish, 'Visible', 'on')
    end
end

% Ensemble created, compute mean and sd, clean up and return
% Set the global gui storage
handles.ens_size    = ens_size;
handles.ens_members = x;

% Update handles structure
guidata(hObject, handles);

% Turn off the data entry messages
set(h_click,  'Visible', 'off');
set(h_finish, 'Visible', 'off');

% Enable the update ensemble button
set(handles.pushbutton_update_ens, 'Enable', 'On');
set(handles.edit_observation,      'Enable', 'On');
set(handles.edit_obs_error_sd,     'Enable', 'On');
set(handles.edit_inflation,        'Enable', 'On');


%----------------------------------------------------------------------



function edit_observation_Callback(hObject, ~, handles)
% hObject    handle to edit_observation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_observation as text
%        str2double(get(hObject,'String')) returns contents of edit_observation as a double

% Turn off any old updated points
set(handles.h_update_ens,     'Visible', 'off');
set(handles.h_inf_up_ens,     'Visible', 'off');
set(handles.h_inf_ens_member, 'Visible', 'off');

% Remove mean and sd of old posterior
clear_labels(handles);

% And the lines in between
set(handles.h_update_lines,        'Visible', 'off');
set(handles.h_inf_lines,           'Visible', 'off');
set(handles.h_inf_axis,            'Visible', 'off');

% Enable things that an error might have turned off
set(handles.edit_obs_error_sd,     'Enable', 'on')
set(handles.edit_inflation,        'Enable', 'on')
set(handles.pushbutton_create_new, 'Enable', 'on')

% Only enable the update ensemble pushbutton if an ensemble has been created
if(handles.ens_size > 0)
    set(handles.pushbutton_update_ens, 'Enable', 'on');
end

% Get the value of the observation
if(isfinite(     str2double(get(hObject, 'String'))))
    observation = str2double(get(hObject, 'String'));
else
    set(handles.edit_observation, 'String', '???');

    % Disable other input to guarantee only one error at a time!
    set(handles.edit_obs_error_sd,     'Enable', 'off')
    set(handles.edit_inflation,        'Enable', 'off')
    set(handles.pushbutton_create_new, 'Enable', 'off')
    set(handles.pushbutton_update_ens, 'Enable', 'off')
    return
end

% Update the global storage
handles.observation = observation;

% Plot the updated distribution
set(handles.h_obs_plot, 'Visible', 'Off');
handles.h_obs_plot = plot_gaussian(handles.observation, handles.obs_error_sd, 1);
set(handles.h_obs_plot, 'Color', 'r', 'Linestyle', '--', 'Linewidth', 2);

% Move the observation asterisk
set(handles.h_obs_ast, 'Visible', 'Off');
handles.h_obs_ast = plot(handles.observation, 0, 'r*', 'MarkerSize', 16,'LineWidth',2.0);

% Set a basic plotting domain range that includes mean +/- 3 obs SDs
xlower = min(handles.observation - 3*handles.obs_error_sd, min(handles.ens_members));
xupper = max(handles.observation + 3*handles.obs_error_sd, max(handles.ens_members));
ylower = -0.4;
yupper = 1.0;
axis([xlower xupper ylower yupper]);

set(gca, 'YTick',      [0 0.2 0.4 0.6 0.8]);
set(gca, 'YTickLabel', [0 0.2 0.4 0.6 0.8]);

hold on
plot([xlower xupper], [0 0], 'k', 'Linewidth', 2);

% Update handles structure
guidata(hObject, handles);


%----------------------------------------------------------------------


% --- Executes during object creation, after setting all properties.
function edit_observation_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_observation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%----------------------------------------------------------------------


function edit_obs_error_sd_Callback(hObject, ~, handles)
% hObject    handle to edit_obs_error_sd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_obs_error_sd as text
%        str2double(get(hObject,'String')) returns contents of edit_obs_error_sd as a double

% Turn off any old updated points
set(handles.h_update_ens,          'Visible', 'off');
set(handles.h_inf_up_ens,          'Visible', 'off');
set(handles.h_inf_ens_member,      'Visible', 'off');

% Remove mean and sd of old posterior
clear_labels(handles);

% And the lines in between
set(handles.h_update_lines,        'Visible', 'off');
set(handles.h_inf_lines,           'Visible', 'off');
set(handles.h_inf_axis,            'Visible', 'off');

% Enable things that an error might have turned off
set(handles.edit_observation,      'Enable', 'on')
set(handles.edit_inflation,        'Enable', 'on')
set(handles.pushbutton_create_new, 'Enable', 'on')

% Only enable the update ensemble pushbutton if an ensemble has been created
if(handles.ens_size > 0)
    set(handles.pushbutton_update_ens, 'Enable', 'on');
end

% Get the value of the observation
if(isfinite(      str2double(get(hObject, 'String'))) && ...
        str2double(get(hObject, 'String')) > 0)
    obs_error_sd = str2double(get(hObject, 'String'));
else
    set(handles.edit_obs_error_sd, 'String', '???');

    % Disable other input to guarantee only one error at a time!
    set(handles.edit_observation,      'Enable', 'off')
    set(handles.edit_inflation,        'Enable', 'off')
    set(handles.pushbutton_create_new, 'Enable', 'off')
    set(handles.pushbutton_update_ens, 'Enable', 'off')
    return
end

% Update the value in global storage
handles.obs_error_sd = obs_error_sd;

% Plot the updated distribution
set(handles.h_obs_plot, 'Visible', 'off');
handles.h_obs_plot = plot_gaussian(handles.observation, handles.obs_error_sd, 1);
set(handles.h_obs_plot, 'Color', 'r', 'Linestyle', '--', 'Linewidth', 2);

% Set a basic plotting domain range that includes mean +/- 3 obs SDs
xlower = min(handles.observation - 3*handles.obs_error_sd, min(handles.ens_members));
xupper = max(handles.observation + 3*handles.obs_error_sd, max(handles.ens_members));
ylower = -0.4;
yupper = 1.0;
axis([xlower xupper ylower yupper]);

set(handles.h_obs_plot, 'Color', 'r', 'Linestyle', '--', 'Linewidth', 2);

set(gca, 'YTick',      [0 0.2 0.4 0.6 0.8]);
set(gca, 'YTickLabel', [0 0.2 0.4 0.6 0.8]);

hold on
plot([xlower xupper], [0 0], 'k', 'Linewidth', 2);

% Update handles structure
guidata(hObject, handles);


%----------------------------------------------------------------------


% --- Executes during object creation, after setting all properties.
function edit_obs_error_sd_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_obs_error_sd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%----------------------------------------------------------------------


% --- Executes on selection change in popupmenu_filter_kind.
function popupmenu_filter_kind_Callback(~, ~, ~)
% hObject    handle to popupmenu_filter_kind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_filter_kind contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_filter_kind


%----------------------------------------------------------------------


% --- Executes during object creation, after setting all properties.
function popupmenu_filter_kind_CreateFcn(hObject, ~, ~)
% hObject    handle to popupmenu_filter_kind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%----------------------------------------------------------------------


% --- Executes on button press in pushbutton_update_ens.
function pushbutton_update_ens_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_update_ens (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Turn off any old points
set(handles.h_update_ens,     'Visible', 'off');
set(handles.h_inf_up_ens,     'Visible', 'off');
set(handles.h_inf_ens_member, 'Visible', 'off');

% Remove mean and sd of old posterior
clear_labels(handles);

% And the lines in between
set(handles.h_update_lines, 'Visible', 'off');
set(handles.h_inf_lines,    'Visible', 'off');
set(handles.h_inf_axis,     'Visible', 'off');

ensemble = handles.ens_members;

% Figure out which filter option is currently selected
h_filter_kind = get(handles.popupmenu_filter_kind);

filter_type = char(h_filter_kind.String(h_filter_kind.Value));

switch filter_type
    case 'EAKF'
        [obs_increments, ~] = ...
            obs_increment_eakf(ensemble, handles.observation, handles.obs_error_sd^2);
    case 'EnKF'
        [obs_increments, ~] = ...
            obs_increment_enkf(ensemble, handles.observation, handles.obs_error_sd^2);
    case 'RHF'
        [obs_increments, ~] = ...
            obs_increment_rhf(ensemble, handles.observation, handles.obs_error_sd^2);
end

% Add on increments to get new ensemble
new_ensemble = ensemble + obs_increments;

y(1:size(ensemble)) = -0.1;
handles.h_update_ens = plot(new_ensemble, y, '*', 'MarkerSize', 16, 'Color', 'Blue');

% Plot lines connecting the prior and posterior ensemble members
for i = 1:size(ensemble, 2)
    x_line = [handles.ens_members(i), new_ensemble(i)];
    y_line = [0, -0.1];
    handles.h_update_lines(i) = plot(x_line, y_line, 'k');
end

% Add in a label of the updated mean and sd
new_mean = mean(new_ensemble);
new_sd = std(new_ensemble);

% Update mean and sd of old posterior
set(handles.text8, 'String', ['Posterior Mean = ', num2str(new_mean)]);
set(handles.text8, 'Visible', 'on');
set(handles.text7, 'String', ['Posterior SD = ', num2str(new_sd)]);
set(handles.text7, 'Visible', 'on');

% If the checkbox isn't set, return now
if(not(get(handles.checkbox_inflation, 'Value')))
    guidata(hObject, handles)
    return
end

% Plot the inflated prior ensemble
y = -0.2;
prior_mean = mean(handles.ens_members(1:handles.ens_size));

inf_ens = zeros(1,handles.ens_size);
for i = 1: handles.ens_size
    inf_ens(i) = (handles.ens_members(i) - prior_mean) * sqrt(handles.inflation) + ...
        prior_mean;
    handles.h_inf_ens_member(i) = plot(inf_ens(i), y, '*', 'MarkerSize', 16, 'Color', [0 0.73 0],'LineWidth',2.0);
end

% Update mean and sd of old posterior
inf_prior_sd = std(inf_ens(1:handles.ens_size));
set(handles.text9,  'String', ['Inflated = ', num2str(prior_mean)]);
set(handles.text9,  'Visible', 'on');
set(handles.text10, 'String', ['Inflated = ', num2str(inf_prior_sd)]);
set(handles.text10, 'Visible', 'on');


% Get the update for the inflated ensemble
switch filter_type
    case 'EAKF'
        [obs_increments, ~] = ...
            obs_increment_eakf(inf_ens, handles.observation, handles.obs_error_sd^2);
    case 'EnKF'
        [obs_increments, ~] = ...
            obs_increment_enkf(inf_ens, handles.observation, handles.obs_error_sd^2);
    case 'RHF'
        [obs_increments, ~] = ...
            obs_increment_rhf(inf_ens, handles.observation, handles.obs_error_sd^2);
end

% Add on increments to get new ensemble
new_ensemble = inf_ens + obs_increments;

y(1:size(ensemble)) = -0.3;
handles.h_inf_up_ens = plot(new_ensemble, y, '*', 'MarkerSize', 16, 'Color', 'Blue');

% Plot lines connecting the prior and posterior ensemble members
for i = 1:size(ensemble, 2)
    x_line = [inf_ens(i), new_ensemble(i)];
    y_line = [-0.2, -0.3];
    handles.h_inf_lines(i) = plot(x_line, y_line, 'k');
end

% Set a basic plotting domain range that includes mean +/- 3 obs SDs
% Plus all inflated members
xlower = min(handles.observation - 3*handles.obs_error_sd, min(inf_ens));
xupper = max(handles.observation + 3*handles.obs_error_sd, max(inf_ens));
ylower = -0.4;
yupper = 1.0;
axis([xlower xupper ylower yupper]);

% Plot the axes for the two priors
plot([xlower xupper], [0 0], 'k', 'Linewidth', 2);
handles.h_inf_axis = plot([xlower xupper], [-0.2 -0.2], 'k', 'Linewidth', 2);

% Update mean and sd of old posterior
update_inf_mean = mean(new_ensemble(1:handles.ens_size));
update_inf_sd = std(new_ensemble(1:handles.ens_size));
set(handles.text12, 'String', ['Inflated = ', num2str(update_inf_mean)]);
set(handles.text12, 'Visible', 'on');
set(handles.text11, 'String', ['Inflated = ', num2str(update_inf_sd)]);
set(handles.text11, 'Visible', 'on');

guidata(hObject, handles)





% --- Executes on button press in checkbox_inflation.
function checkbox_inflation_Callback(~, ~, ~)
% hObject    handle to checkbox_inflation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_inflation




function clear_labels(handles)

% Turns off all labels except for the prior mean and SD
set(handles.text7,  'Visible', 'off');
set(handles.text8,  'Visible', 'off');
set(handles.text9,  'Visible', 'off');
set(handles.text10, 'Visible', 'off');
set(handles.text11, 'Visible', 'off');
set(handles.text12, 'Visible', 'off');




function edit_inflation_Callback(hObject, ~, handles)
% hObject    handle to edit_inflation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_inflation as text
%        str2double(get(hObject,'String')) returns contents of edit_inflation as a double

% Turn off any old updated points
set(handles.h_update_ens,          'Visible', 'off');
set(handles.h_inf_up_ens,          'Visible', 'off');
set(handles.h_inf_ens_member,      'Visible', 'off');

% Remove mean and sd of old posterior
clear_labels(handles);

% And the lines in between
set(handles.h_update_lines,        'Visible', 'off');
set(handles.h_inf_lines,           'Visible', 'off');
set(handles.h_inf_axis,            'Visible', 'off');

% Enable things that an error might have turned off
set(handles.edit_observation,      'Enable', 'on')
set(handles.edit_obs_error_sd,     'Enable', 'on')
set(handles.pushbutton_create_new, 'Enable', 'on')

% Only enable the update ensemble pushbutton if an ensemble has been created
if(handles.ens_size > 0)
    set(handles.pushbutton_update_ens, 'Enable', 'on');
end

% Get the value of the observation
if(isfinite(   str2double(get(hObject, 'String'))) && ...
        str2double(get(hObject, 'String')) > 0)
    inflation = str2double(get(hObject, 'String'));
else
    set(handles.edit_inflation, 'String', '???');

    % Disable other input to guarantee only one error at a time!
    set(handles.edit_observation,      'Enable', 'off')
    set(handles.edit_obs_error_sd,     'Enable', 'off')
    set(handles.pushbutton_create_new, 'Enable', 'off')
    set(handles.pushbutton_update_ens, 'Enable', 'off')
    return
end

% Update the value in global storage
handles.inflation = inflation;

% Plot the updated distribution
set(handles.h_obs_plot, 'Visible', 'off');
handles.h_obs_plot = plot_gaussian(handles.observation, handles.obs_error_sd, 1);
set(handles.h_obs_plot, 'Color', 'r', 'Linestyle', '--', 'Linewidth', 2);

% Set a basic plotting domain range that includes mean +/- 3 obs SDs
xlower = min(handles.observation - 3*handles.obs_error_sd, min(handles.ens_members));
xupper = max(handles.observation + 3*handles.obs_error_sd, max(handles.ens_members));
ylower = -0.4;
yupper = 1.0;
axis([xlower xupper ylower yupper]);

set(handles.h_obs_plot, 'Color', 'r', 'Linestyle', '--', 'Linewidth', 2);

set(gca, 'YTick',      [0 0.2 0.4 0.6 0.8]);
set(gca, 'YTickLabel', [0 0.2 0.4 0.6 0.8]);

hold on
plot([xlower xupper], [0 0], 'k', 'Linewidth', 2);

% Update handles structure
guidata(hObject, handles);







% --- Executes during object creation, after setting all properties.
function edit_inflation_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_inflation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

