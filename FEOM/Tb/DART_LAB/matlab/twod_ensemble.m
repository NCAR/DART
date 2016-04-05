function varargout = twod_ensemble(varargin)
%% TWOD_ENSEMBLE demonstrates a powerful aspect of ensemble data assimilation.
%      The ensemble makes it possible to estimate the impact of observations
%      on unobserved state variables.
%
%      Click on the 'Create New Ensemble' button to activate the interactive
%      observation generation mechanism and lay down a set of ensemble
%      samples of an unobserved variable (vertical axis) and an observed
%      variable (horizontal axis). The ensemble members are created by
%      left clicking in the central portion of the figure window.
%      Start out small, say 6 or so.
%      In this case, some H() operator would generate the Observed Quantity.
%      The Unobserved State Variable could simply be some portion of the
%      model state.
%
%      After creating the ensemble, the correlation between the Observed
%      Quantity and the Unobserved State Variable is calculated.
%      Select an assimilation algorithm and click 'Update Ensemble'.
%      The increments are shown as red lines, the new Posterior estimates
%      are in blue.
%
% See also: gaussian_product, oned_model, oned_ensemble, run_lorenz_63,
%           run_lorenz_96

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @twod_ensemble_OpeningFcn, ...
    'gui_OutputFcn',  @twod_ensemble_OutputFcn, ...
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



%----------------------------------------------------------------------



%% --- Executes just before twod_ensemble is made visible.
function twod_ensemble_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to twod_ensemble (see VARARGIN)

help twod_ensemble

% set random number seed to same value to generate known sequences
rng('default')

% Choose default command line output for twod_ensemble
handles.output = hObject;

% Insert the ensemble structure into this
handles.ens_size        = 0;
handles.ens_members     = 0;
handles.h_update_ens    = 0;
handles.h_ens_member    = 0;
handles.h_best_fit      = 0;
handles.h_marg_obs_plot = 0;
handles.h_obs_ast       = 0;
handles.h_obs_marg      = 0;
handles.h_gui_marg      = 0;
handles.h_unobs         = 0;
handles.h_marg          = 0;
handles.h_marg_update   = 0;
handles.h_marg_inc      = 0;
handles.h_marg_state    = 0;
handles.h_state_inc     = 0;
handles.h_joint_update  = 0;
handles.h_joint_inc     = 0;
handles.h_correl        = 0;
handles.first_correl    = true;

% Also include the subplot handles r1, r2, r3
handles.r1 = 0;
handles.r2 = 0;
handles.r3 = 0;

% Update handles structure
guidata(hObject, handles);

% Go ahead and plot the initial observational error distribution
h_observation  = get(handles.edit1);
h_obs_error_sd = get(handles.edit2);
observation    = str2double(h_observation.String);
obs_error_sd   = str2double(h_obs_error_sd.String);

% Plot this on the marginal plot on the gui figure
handles.h_marg_obs_plot = plot_gaussian(observation, obs_error_sd, 1);
axis([0 10 -0.2 1]);

% Set the ticks
set(gca, 'YTick', [0 0.2 0.4 0.6 0.8]);
set(gca, 'YTickLabel', [0 0.2 0.4 0.6 0.8]);
set(handles.h_marg_obs_plot, 'Color', 'r', 'Linestyle', '--', ...
    'LineWidth', 2);
xlabel('Observed Quantity', 'Fontsize', 14);
ylabel('Observation Likelihood', 'FontSize', 14);
hold on

% Plot an asterisk for the observed value
handles.h_obs_ast = plot(observation, 0, 'r*', 'MarkerSize', 16,'LineWidth',2.0);

% Plot an axis; display is fixed from x = 0 to 10
plot([0 10], [0 0], 'k', 'LineWidth', 2);

% Setup the joint distribution plot plus the two marginals
figure(1); clf

% Start with unobserved state variable marginal
% TJH POSSIBLE IMPROVEMENT ... annotate initial mean, sd
handles.r1 = subplot(2, 2, 1);
axis([0 1 0 10]);
ylabel('Unobserved State Variable', 'Fontsize', 14);

hold on

% Get a subplot for the unobserved state variable marginal
handles.r2 = subplot(2, 2, 2);
axis([0 10 0 10]);
grid on
hold on

% Get a subplot for the observed variable marginal
handles.r3 = subplot(2, 2, 4);
% May want to mess with vertical axis for prior density
axis([0 10 0 1]);
hold on
xlabel('Observed Quantity', 'Fontsize', 14);
% Plot the observation value as an asterisk
handles.h_obs_marg = plot(observation, 0, 'r*', 'MarkerSize', 16,'LineWidth',2.0);

% Position the marginal and joint plot boxes
set(handles.r1, 'Position', [0.14, 0.23, 0.10, 0.6550]);
set(handles.r2, 'Position', [0.25, 0.23, 0.6550, 0.6550]);
set(handles.r3, 'Position', [0.25, 0.1100, 0.6550, 0.1000]);

% Turn off the unwanted tick labels
set(handles.r2, 'Xticklabel', []);
set(handles.r2, 'Yticklabel', []);
set(handles.r1, 'Xticklabel', []);
set(handles.r3, 'Yticklabel', []);


% Update handles structure
guidata(hObject, handles);

% Reset focus to the menu gui window
% Setting the axes clears the legend, gcbo restores focus
axes(handles.axes1);

% UIWAIT makes twod_ensemble wait for user response (see UIRESUME)
% uiwait(handles.figure1);



%----------------------------------------------------------------------



%% --- Outputs from this function are returned to the command line.
function varargout = twod_ensemble_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



%----------------------------------------------------------------------



%% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, ~, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disable the update ensemble button and all other active buttons
set(handles.pushbutton1, 'Enable', 'Off');
set(handles.pushbutton2, 'Enable', 'Off');
set(handles.edit1, 'Enable', 'Off');
set(handles.edit2, 'Enable', 'Off');

% Clear out any old ensemble members if they exist
figure(1)
for i = 1:handles.ens_size
    set(handles.h_ens_member(i), 'Visible', 'off');
    set(handles.h_gui_marg(i), 'Visible', 'off');
    set(handles.h_unobs(i), 'Visible', 'off');
    set(handles.h_marg(i), 'Visible', 'off');
end

% Turn off any posterior old plotting
set(handles.h_update_ens, 'Visible', 'off');
set(handles.h_marg_update, 'Visible', 'off');
set(handles.h_marg_inc, 'Visible', 'off');
set(handles.h_marg_state, 'Visible', 'off');
set(handles.h_state_inc, 'Visible', 'off');
set(handles.h_joint_update, 'Visible', 'off');
set(handles.h_joint_inc, 'Visible', 'off');

% Clear out the old best fit line
set(handles.h_best_fit, 'Visible', 'off');

% Remove correlation of old ensemble
%%%set(handles.text2, 'String', 'Correlation = ');

% If this is not the first creation of an ensemble, remove the old label
if(handles.first_correl)
    handles.first_correl = false;
else
    set(handles.h_correl, 'Visible', 'off');
end

% Work in the joint distribution plot
subplot(handles.r2);
hold on

h_click = text(1, 2, 'Click to create member', 'FontSize', 16);

% Need to guarantee at least 2 ensemble members
ens_size = 0;

h_err_text = text(1, 3, 'Click inside graphics box to create member', ...
    'Color', 'r', 'FontSize', 16, 'Visible', 'off');

while ens_size < 2
    [xt, yt] = ginput(1);

    % Clear out the error message if it's been made visible
    set(h_err_text, 'Visible', 'off');

    % Need to make sure that we were clicked in the right subplot
    % and with proper limits
    if(xt < 0 || xt > 10 || yt < 0 || yt > 10 || gca ~= handles.r2)
        set(h_err_text, 'Visible', 'on');
        subplot(handles.r2);
    else
        ens_size = ens_size + 1;
        x(1, ens_size) = xt;
        x(2, ens_size) = yt;
        subplot(handles.r2);
        handles.h_ens_member(ens_size) = ...
            plot(x(1, ens_size), x(2, ens_size), '*', 'MarkerSize', 16, 'Color', [0 0.73 0],'LineWidth',2.0);

        % Plot the marginal for the unobserved state variable
        % TJH POSSIBLE IMPROVEMENT ... annotate new marginal mean, sd
        subplot(handles.r1);
        handles.h_unobs(ens_size) = ...
            plot(0, x(2, ens_size), '*', 'MarkerSize', 16, 'Color', [0 0.73 0],'LineWidth',2.0);

        % Plot the marginal for the observed quantity
        subplot(handles.r3);
        handles.h_marg(ens_size) = ...
            plot(x(1, ens_size), 0, '*', 'MarkerSize', 16, 'Color', [0 0.73 0],'LineWidth',2.0);

        % Plot the marginal in the gui frame
        axes(handles.axes1);
        handles.h_gui_marg(ens_size) = ...
            plot(x(1, ens_size), 0, '*', 'MarkerSize', 16, 'Color', [0 0.73 0],'LineWidth',2.0);
        % Then switch back to figure(1)
        figure(1);

    end
end

subplot(handles.r2);
h_finish = text(1, 1, 'Click outside of plot to finish', 'Fontsize', 16);

while ens_size < 1000
    [xt, yt] = ginput(1);
    % Make sure that the click was in the correct set of axes
    gca;
    % Terminate by clicking outside of graph range
    if(xt > 10 || xt < 0 || yt > 10 || yt < 0 || gca ~= handles.r2)
        subplot(handles.r2);
        break;
    else
        ens_size = ens_size + 1;
        x(1, ens_size) = xt;
        x(2, ens_size) = yt;

        subplot(handles.r2);
        handles.h_ens_member(ens_size) = ...
            plot(x(1, ens_size), x(2, ens_size), '*', 'MarkerSize', 16, 'Color', [0 0.73 0],'LineWidth',2.0);

        % Plot the marginal for the unobserved state variable
        % TJH POSSIBLE IMPROVEMENT ... annotate new marginal mean, sd
        subplot(handles.r1);
        handles.h_unobs(ens_size) = ...
            plot(0, x(2, ens_size), '*', 'MarkerSize', 16, 'Color', [0 0.73 0],'LineWidth',2.0);

        % Plot the marginal for the observed quantity
        subplot(handles.r3);
        handles.h_marg(ens_size) = ...
            plot(x(1, ens_size), 0, '*', 'MarkerSize', 16, 'Color', [0 0.73 0],'LineWidth',2.0);

        % Plot the marginal in the gui frame
        axes(handles.axes1);
        handles.h_gui_marg(ens_size) = ...
            plot(x(1, ens_size), 0, '*', 'MarkerSize', 16, 'Color', [0 0.73 0],'LineWidth',2.0);

        % Then switch back to figure(1)
        figure(1);

        % Display the prior correlation
        prior_correl = corrcoef(x(1, :), x(2, :));
        %%%set(handles.text2, 'String', ...
        %%%['Correlation = ', num2str(prior_correl(1, 2))]);

        subplot(handles.r2)
        if(ens_size > 1)
            set(handles.h_correl, 'Visible', 'Off');
        end
        handles.h_correl = text(5, 9, ['Correlation = ', num2str(prior_correl(1, 2))], ...
            'FontSize', 16, 'Color', [0 0.73 0], 'FontWeight', 'Bold');

    end

end

%% Ensemble created, compute mean and sd, clean up and return
% Set the global gui storage
handles.ens_size = ens_size;
handles.ens_members = x;

% Plot the best fit line on the ensemble
prior_mean = mean(x, 2);
prior_cov = cov(x(1, :), x(2, :));
slope = prior_cov(1, 2) / var(x(1, :));

best_x = [0 10];
best_y(1) = prior_mean(2) - (prior_mean(1)) * slope;
best_y(2) = best_y(1) + 10 * slope;
handles.h_best_fit = plot(best_x, best_y, 'g', 'LineWidth', 2.0);
set(handles.h_best_fit, 'Color', [0 0.73 0]);

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

% Reset focus to the menu gui window
% Setting the axes clears the legend, gcbo restores focus
%%%axes(handles.axes1);
[~, gcbo_fig] = gcbo;
figure(gcbo_fig);



%----------------------------------------------------------------------



function edit1_Callback(hObject, ~, handles)
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
set(handles.h_marg_obs_plot, 'Visible', 'off');
handles.h_marg_obs_plot = plot_gaussian(observation, obs_error_sd, 1);
set(handles.h_marg_obs_plot, 'Color', 'r', 'Linestyle', '--', 'Linewidth', 2);

% Update the observation asterisk
set(handles.h_obs_ast, 'Visible', 'off');
handles.h_obs_ast = plot(observation, 0, 'r*', 'MarkerSize', 16,'LineWidth',2.0);

% Plot the updated obs distribution on the marginal subplot
figure(1);
subplot(handles.r3);
% Plot the updated observation in the marginal
set(handles.h_obs_marg, 'Visible', 'off');
handles.h_obs_marg = plot(observation, 0, 'r*', 'MarkerSize', 16,'LineWidth',2.0);

% Update handles structure
guidata(hObject, handles);

% Reset focus to the menu gui window
% Setting the axes clears the legend, gcbo restores focus
%%%axes(handles.axes1);
[~, gcbo_fig] = gcbo;
figure(gcbo_fig);



%----------------------------------------------------------------------



%% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, ~, ~)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%----------------------------------------------------------------------



function edit2_Callback(hObject, ~, handles)
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

% Plot the updated distribution on the menu plot
set(handles.h_marg_obs_plot, 'Visible', 'off');
handles.h_marg_obs_plot = plot_gaussian(observation, obs_error_sd, 1);
set(handles.h_marg_obs_plot, 'Color', 'r', 'Linestyle', '--', 'Linewidth', 2);

% Update the observation asterisk
set(handles.h_obs_ast, 'Visible', 'off');
handles.h_obs_ast = plot(observation, 0, 'r*', 'MarkerSize', 16,'LineWidth',2.0);

% Plot the updated distribution
figure(1);
subplot(handles.r3);
% Plot the updated observation in the marginal
set(handles.h_obs_marg, 'Visible', 'off');
handles.h_obs_marg = plot(observation, 0, 'r*', 'MarkerSize', 16,'LineWidth',2.0);

% Update handles structure
guidata(hObject, handles);

% Reset focus to the menu gui window
% Setting the axes clears the legend, gcbo restores focus
%%%axes(handles.axes1);
[~, gcbo_fig] = gcbo;
figure(gcbo_fig);



%----------------------------------------------------------------------



%% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, ~, ~)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%----------------------------------------------------------------------



%% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(~, ~, ~)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1



%----------------------------------------------------------------------



%% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, ~, ~)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%----------------------------------------------------------------------



%% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, ~, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Turn off any old points
set(handles.h_update_ens, 'Visible', 'off');
set(handles.h_marg_update, 'Visible', 'off');
set(handles.h_marg_inc, 'Visible', 'off');
set(handles.h_marg_state, 'Visible', 'off');
set(handles.h_state_inc, 'Visible', 'off');
set(handles.h_joint_update, 'Visible', 'off');
set(handles.h_joint_inc, 'Visible', 'off');

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
        [obs_increments, ~] = ...
            obs_increment_eakf(ensemble(1, :), observation, obs_error_sd^2);
    case 'EnKF'
        [obs_increments, ~] = ...
            obs_increment_enkf(ensemble(1, :), observation, obs_error_sd^2);
    case 'RHF'
        [obs_increments, ~] = ...
            obs_increment_rhf(ensemble(1, :), observation, obs_error_sd^2);
end

% Add on increments to get new ensemble
new_ensemble = ensemble(1, :) + obs_increments;

y(1:handles.ens_size) = -0.1;
handles.h_update_ens = plot(new_ensemble, y, '*', 'MarkerSize', 16, 'Color', 'Blue');

% Plot the increments in the state marginal plot
figure(1);
subplot(handles.r3);

% Need to sort ensemble to get nice ordering for increments
[~, sort_obs_ind] = sort(ensemble(1, :));
for i = 1:handles.ens_size
    y(i) = i / (handles.ens_size + 1);
    handles.h_marg_update(i) = ...
        plot(new_ensemble(sort_obs_ind(i)), y(i), '*', 'MarkerSize', 16, 'Color', 'Blue');
    % Also plot a segment in blue
    handles.h_marg_inc(i) = ...
        plot([ensemble(1, sort_obs_ind(i)), new_ensemble(1, sort_obs_ind(i))], ...
        [y(i), y(i)], 'r');
end

% Figure out the increments for the unobserved variable
covar = cov(ensemble');
state_inc = obs_increments * covar(1, 2) / covar(1, 1);
new_state = ensemble(2, :) + state_inc;
% TJH POSSIBLE IMPROVEMENT ... annotate new marginal mean, sd


% Now need to sort the state variable ensemble to get nice ordering
subplot(handles.r1);
[~, sort_ind] = sort(ensemble(2, :));
for i = 1:handles.ens_size
    handles.h_marg_state(i) = ...
        plot(y(i), new_state(sort_ind(i)), '*', 'MarkerSize', 16, 'Color', 'Blue');
    % Also plot a segment in blue
    handles.h_state_inc(i) = plot([y(i), y(i)], ...
        [ensemble(2, sort_ind(i)), new_state(sort_ind(i))], 'r');
end

% Plot the updated joint distribution points
subplot(handles.r2);
for i = 1:handles.ens_size
    handles.h_joint_update(i) = plot(new_ensemble(i), new_state(i), ...
        '*', 'MarkerSize', 16, 'Color', 'Blue');
    handles.h_joint_inc(i) = plot([ensemble(1, i), new_ensemble(1, i)], ...
        [ensemble(2, i), new_state(i)], 'r');
end

% Return the focus to the window with pushbuttons
axes(handles.axes1);

guidata(hObject, handles);

% Reset focus to the menu gui window
% Setting the axes clears the legend, gcbo restores focus
%%%axes(handles.axes1);
[~, gcbo_fig] = gcbo;
figure(gcbo_fig);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

