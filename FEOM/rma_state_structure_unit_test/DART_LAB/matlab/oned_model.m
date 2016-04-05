function varargout = oned_model(varargin)
%% ONED_MODEL simple ensemble data assimilation example.
%
%      There are no input arguments. Simply invoke by typing the name.
%
%      ONED_MODEL demonstrates the simplest possible case of ensemble data
%      assimilation. It is possible to explore assimilation algorithms,
%      ensemble sizes, model biases, etc. on-the-fly. The posterior
%      of the state is indicated by blue asterisks, the states evolve along
%      a trajectory indicated by the green lines to wind up at a prior state
%      for the assimilation - indicated by the green asterisks. After the
%      assimilation, the (posterior) state is indicated in blue and the
%      process is ready to repeat.
%
%      ONED_MODEL opens two windows. A gui control window that also plots
%      the most recent prior, posterior, and observation, and a figure
%      window that plots time sequences of the assimilation, the RMS error,
%      spread and kurtosis, and prior and posterior rank histograms.
%
%      The top button alternates between "Advance Model" and "Assimilate" to
%      single-step the model. The "Start Free Run" button is useful to watch
%      the system evolve and generate estimates from many assimilation cycles.
%
%      Since this is a 'perfect model' experiment, we know the true state,
%      the amount of noise added to the observations, etc.; so it is possible to
%      calculate the error of the ensemble in addition to the spread. The
%      Truth is not (in general) the same as the observation!
%
%      This also introduces the concept of the 'rank histogram'. With N
%      ensemble members, there will always be N+1 'bins' that encompass the
%      state (to the left of the lowest ensemble member, to the right of the
%      highest, and all the ones in-between). The rank histogram tracks the
%      number of times the truth lies in each bin. With a good ensemble,
%      (i.e. a good assimilation system) the true state will be
%      indistinguishable from any of the ensemble members. This results
%      in a flat rank histogram, given enough samples.
%
% See also: gaussian_product, oned_ensemble, twod_ensemble, run_lorenz_63,
%           run_lorenz_96

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Last Modified by GUIDE v2.5 27-Aug-2009 08:54:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @oned_model_OpeningFcn, ...
    'gui_OutputFcn',  @oned_model_OutputFcn, ...
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




% --- Executes just before oned_model is made visible.
function oned_model_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to oned_model (see VARARGIN)

help oned_model

% set random number seed to same value to generate known sequences
rng('default')

% Choose default command line output for oned_model
handles.output = hObject;

% Set up global storage with initial values
handles.ens_size = 4;
handles.time_step = 1;
handles.ready_to_advance = true;
handles.alpha = 0;
handles.obs_error_sd = 1;
handles.observation = 0;
handles.error = 0;
handles.spread = 0;
handles.kurtosis = 0;
handles.model_bias = 0;
handles.inflation = 1.0;

% An array to keep track of rank histograms
handles.prior_rank(1 : handles.ens_size + 1) = 0;
handles.posterior_rank(1 : handles.ens_size + 1) = 0;

% Handles to subregions of plots
handles.r1 = 0;
handles.r2 = 0;
handles.r3 = 0;
handles.r4 = 0;
handles.r5 = 0;

% Set up initial axes for the gui window
y_max = 1 / (sqrt(2 * pi) * handles.obs_error_sd);
axis([-4 4 -0.2 y_max + 0.02]);

% Turn off the negative labels
set(gca, 'YTick', [0 0.1 0.2 0.3 0.4]);
set(gca, 'YTickLabel', [0 0.1 0.2 0.3 0.4]);

% Work in a figure window for three timeseries plots
figure(1); clf;

% Subplot for ensemble time series
handles.r1 = subplot(4, 1, 1);

% Draw initial ensemble values from Normal(0, 1)
handles.ens = randn([1 4]);
x(1:handles.ens_size) = handles.time_step + 0.1;
plot(x, handles.ens, 'b*', 'MarkerSize', 6);
axis([1 10 -4 4]);
hold on;

% Include the 0 line as the truth for all times
plot([1 100000], [0 0], 'k--');

% Want the y axis limits to take care of themselves
set(gca, 'YLimMode', 'Auto');

ylabel('State', 'FontSize', 14);

% Next subplot is for mean and spread
handles.r2 = subplot(4, 1, 2);
axis([1 10 0 4]);
hold on;
ylabel('Error, Spread', 'FontSize', 14);

%  Compute the initial error (truth is 0.0) and spread (standard deviation)
handles.error  = calculate_rmse(handles.ens, 0.0);
handles.spread = std(handles.ens);

% calculate rmse and spread ... expectation over a long time is
% that they would be the same.

% Third subplot is for kurtosis
handles.r3 = subplot(4, 1, 3);
axis([1 10 0 4]);
hold on
xlabel('Timestep', 'FontSize', 14);
ylabel('Kurtosis', 'FontSize', 14);

% Compute initial kurtosis
handles.kurtosis = kurt(handles.ens);

% Fourth and fifth subplot are for rank histogram
handles.r4 = subplot(4, 2, 7);
ylabel('Frequency');
xlabel('Rank');
title 'Prior Rank Histogram';

handles.r5 = subplot(4, 2, 8);
ylabel('Frequency');
xlabel('Rank');
title 'Posterior Rank Histogram';

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes oned_model wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = oned_model_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, ~, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Move the model ahead a step or assimilate observations as appropriate
step_ahead(hObject, handles)




function edit1_Callback(hObject, ~, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

% Get the value of the model_bias
handles.model_bias = str2double(get(hObject, 'String'));
if(not(isfinite(handles.model_bias)))
    % Indicate input error in text box
    set(handles.edit1, 'String', '???');
    
    % After this, only this edit box will work
    turn_off_controls(handles);
    set(handles.edit1, 'Enable', 'On');
    
    return
end

% Turn on all controls if successful
turn_on_controls(handles)
% Update handles structure
guidata(hObject, handles);





% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, ~, ~)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(~, ~, ~)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double




% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, ~, ~)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(~, ~, ~)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double



% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, ~, ~)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(~, ~, ~)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, ~, ~)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, ~, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double

% Get the value of the model nonlinearity parameter
handles.alpha= str2double(get(hObject, 'String'));
if(not(isfinite(handles.alpha)) ||  handles.alpha < 0)
    % Indicate input error in text box
    set(handles.edit5, 'String', '???');
    
    % After this, only this edit box will work
    turn_off_controls(handles);
    set(handles.edit5, 'Enable', 'On');
    
    return
end

% Enable all controls
turn_on_controls(handles);
% Update handles structure
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, ~, ~)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, ~, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% Change the ensemble size if input value is numeric.
if(not(isfinite(str2double(get(hObject, 'String')))))
    set(handles.edit6, 'String', '???');
    
    % Disable other input to guarantee only one error at a time!
    turn_off_controls(handles); % After this, only this edit box will work
    set(handles.edit6, 'Enable', 'On'); 
    return
end

% Get the value of the ensemble size
new_ens_size = round(str2double(get(hObject, 'String')));
if(new_ens_size < 2)
    set(handles.edit6, 'String', '???');
    
    % Disable other input to guarantee only one error at a time!
    turn_off_controls(handles); % After this, only this edit box will work
    set(handles.edit6, 'Enable', 'On');
    return
else
    set(handles.edit6,'String',sprintf('%d',new_ens_size))
end

% Legal value for ensemble size; enable all controls
turn_on_controls(handles);

% Reset the histograms
handles = reset_histograms(new_ens_size, hObject, handles);

% Generate a new ensemble by truncating old ensemble OR adding new
if(new_ens_size == handles.ens_size)
    return
elseif(new_ens_size < handles.ens_size)
    % Get rid of extra ensemble members, recompute mean, spread and kurtosis
    handles.ens = handles.ens(1:new_ens_size);
    handles.ens_size = new_ens_size;
else
    % Add new ensemble members drawn from present distribution
    handles.ens(handles.ens_size + 1 : new_ens_size) = ...
        randn([1 new_ens_size - handles.ens_size]);
    handles.ens_size = new_ens_size;
end

% Update moments
handles.error    = calculate_rmse(handles.ens, 0.0);
handles.spread   = std(handles.ens);
handles.kurtosis = kurt(handles.ens);

% Update handles structure
guidata(hObject, handles);





% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, ~, ~)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%----------------------------------------
function[handles] = reset_histograms(ens_size, hObject, handles)

% Clear out the histograms and reset storage

% Need to reset the rank histograms when ensemble size changes
handles.prior_rank(1 : ens_size + 1) = 0;
handles.posterior_rank(1 : ens_size + 1) = 0;

% Clear the histograms
figure(1);
subplot(handles.r4);
bar(handles.prior_rank(1:ens_size + 1));
ylabel('Frequency');
xlabel('Rank');
title 'Prior Rank Histogram';
axis tight;

subplot(handles.r5);
bar(handles.posterior_rank(1:ens_size + 1));
ylabel('Frequency');
xlabel('Rank');
title 'Posterior Rank Histogram';
axis tight;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in reset_pushbutton.
function reset_pushbutton_Callback(hObject, ~, handles)
% hObject    handle to reset_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Reset all the figures and the data structures
% Keep the current filter type, ensemble size and obs characteristics

% Reset the time to 1 and be ready to advance
handles.time_step = 1;
handles.ready_to_advance = true;

% Set the pushbutton to say Assimilate Obs
set(handles.pushbutton1, 'String', 'Advance Model');

% Reset the histograms
handles = reset_histograms(handles.ens_size, hObject, handles);

% Reset the ensemble time series
figure(1)
subplot(handles.r1);
hold off

% Draw reset ensemble values from Normal(0, 1)
handles.ens = randn([1 handles.ens_size]);
x(1:handles.ens_size) = handles.time_step + 0.1;
plot(x, handles.ens, 'b*', 'MarkerSize', 6);
hold on;
axis([1 10 -4 4]);

% Include the 0 line as the truth for all times
plot([1 100000], [0 0], 'k--');

% Want the y axis limits to take care of themselves
set(gca, 'YLimMode', 'Auto');

ylabel('State', 'FontSize', 14);


% Reset the error and spread plot
handles.error    = 0;
handles.spread   = 0;
handles.kurtosis = 0;

subplot(handles.r2);
hold off
% Reset the plot and set up the colors for the legend
plot(1, 1, 'b', 'visible', 'off');
hold on
plot(1, 1, 'r', 'visible', 'off');
axis([1 10 0 4]);
ylabel('Error, Spread', 'FontSize', 14);

%  Compute the initial error (truth is 0.0) and spread (standard deviation)

handles.error  = calculate_rmse(handles.ens, 0.0);
handles.spread = std(handles.ens);

handles.r3 = subplot(4, 1, 3);
hold off
plot(1, 1, 'visible', 'off');
axis([1 10 0 4]);
hold on
xlabel('Timestep', 'FontSize', 14);
ylabel('Kurtosis', 'FontSize', 14);

% Compute initial kurtosis
handles.kurtosis = kurt(handles.ens);

% Reset focus to the menu gui window
[~, gcbo_fig] = gcbo;
figure(gcbo_fig);


% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Turn off all the other controls to avoid a mess
turn_off_controls(handles)
set(handles.pushbutton_run, 'Enable', 'On');

if(strcmp(get(hObject, 'String'), 'Stop Free Run'))
    % Being told to stop; switch to not running status
    set(hObject, 'String', 'Start Free Run');
    % Update handles structure
    guidata(hObject, handles);
else
    % Being told to start run
    % Change the button to 'Stop Free Run')
    set(hObject, 'String', 'Stop Free Run');
    % Update the handles structure
    guidata(hObject, handles);
    
    % Start looping
    for i = 1:1000
        % Check to see if stop ; get latest copy of global data
        my_data = guidata(gcbo);
        status_string = get(my_data.pushbutton_run, 'String');
        if(strcmp(status_string, 'Start Free Run'))
            % Turn all the other controls back on
            turn_on_controls(handles);
            return
        end
        % Do the next advance or assimilation step
        step_ahead(hObject, my_data)
        pause(0.2)
    end
end

% Turn all the other controls back on
turn_on_controls(handles);


%----------- Moves the model ahead or assimilates next observations ------
function step_ahead(hObject, handles)

% Start out working in ensemble time series plot
figure(1);
subplot(handles.r1);

% If this is an advance, get and plot new model state, advance time
if(handles.ready_to_advance)
    % Set to do an assimilation next time
    handles.ready_to_advance = false;
    % Advance the model
    ens_new = advance_oned(handles.ens, handles.alpha, handles.model_bias);
    
    % Inflate the model
    ens_new_mean = mean(ens_new);
    ens_new = (ens_new - ens_new_mean) * sqrt(handles.inflation) + ens_new_mean;
    
    
    handles.time_step = handles.time_step + 1;
    h = plot(handles.time_step - 0.1, ens_new, '*', 'MarkerSize', 6);
    set(h, 'Color', [0 0.73 0]);
    for i = 1:handles.ens_size
        h = plot([handles.time_step - 1 + + 0.1, handles.time_step - 0.1], ...
            [handles.ens(i), ens_new(i)]);
        set(h, 'Color', [0 0.73 0]);
    end
    
    % Get the range of x on the plot to use for the gui plot
    ens_axis = axis;
    
    % Plot the segment for the prior error
    subplot(handles.r2);
    prior_error = calculate_rmse(ens_new, 0.0);
    plot([handles.time_step - 1 + 0.1, handles.time_step - 0.1], ...
        [handles.error, prior_error],'LineWidth',2.0);
    hold on;
    
    handles.error = prior_error;
    
    % Plot the segement for the prior spread
    prior_spread = std(ens_new);
    plot([handles.time_step - 1 + 0.1, handles.time_step - 0.1], ...
        [handles.spread, prior_spread], 'r', 'LineWidth', 2.0);
    handles.spread = prior_spread;
    
    % Put on a legend
    legend('Error', 'Spread', 'Location', 'NorthEast');
    legend boxoff
    
    % Want the lower y limit to stay 0 for error spread
    set(gca, 'YLimMode', 'Auto');
    y_lim = get(gca, 'YLim');
    y_lim(1) = 0;
    set(gca, 'YLim', y_lim);
    
    % Plot the segement for the prior kurtosis
    subplot(handles.r3);
    prior_kurtosis= kurt(ens_new);
    plot([handles.time_step - 1 + 0.1, handles.time_step - 0.1], ...
        [handles.kurtosis, prior_kurtosis], 'r','LineWidth',2);
    handles.kurtosis= prior_kurtosis;
    
    % Want the lower y limit to stay 0 for error spread
    set(gca, 'YLimMode', 'Auto');
    y_lim = get(gca, 'YLim');
    y_lim(1) = 0;
    set(gca, 'YLim', y_lim);
    
    % Update the rank data
    subplot(handles.r4);
    ens_rank = get_ens_rank(ens_new, 0);
    %fprintf([sprintf('\ntimestep %d bin edges are ',handles.time_step), ...
    %         num2str(sort([-Inf ens_new Inf]),'%10.4f'),'\n'])
    %fprintf('timestep %d bin/"rank" is %d\n',handles.time_step, ens_rank)
    
    % Plot the latest rank entry as a different color
    temp_rank(:, 1) = handles.prior_rank(1:handles.ens_size + 1);
    temp_rank(:, 2) = 0;
    temp_rank(ens_rank, 2) = 1;
    
    bar(temp_rank, 'stacked');
    ylabel('Frequency');
    xlabel('Rank');
    title 'Prior Rank Histogram';
    axis tight;
    
    % Plot the gui window for this update
    axes(handles.axes1);
    hold off;
    % Put on a black axis line
    plot([ens_axis(3), ens_axis(4)], [0, 0], 'k', 'Linewidth', 2);
    hold on;
    
    % Plot the prior ensemble members in green
    % Plotting ticks instead of asterisks makes bins clearer
    tick_half = 0.015;
    for n_tick = 1:handles.ens_size
        hold on;
        plot([ens_new(n_tick), ens_new(n_tick)], ...
            [-tick_half, tick_half], 'Color', [0 0.73 0], ...
            'LineWidth', 2);
    end
    
    % The height of the obs likelihood controls the vertical axis
    y_max = 1 / (sqrt(2 * pi) * handles.obs_error_sd);
    % Set the axis for the gui window
    axis([ens_axis(3), ens_axis(4), -0.2 y_max + 0.02]);
    
    % Turn off the negative labels
    set(gca, 'YTick', [0 0.1 0.2 0.3 0.4]);
    set(gca, 'YTickLabel', [0 0.1 0.2 0.3 0.4]);
    
    
    % Plot the truth (at 0) as a tick
    plot([0 0], [-0.02 0.02], 'k', 'LineWidth', 2);
    
    % Put in some information about the bin, x position tricky
    base_x = max(0, ens_axis(3));
    text_width = (ens_axis(4) - ens_axis(3)) / 3;
    if((base_x + text_width) > ens_axis(4))
        base_x = ens_axis(4) - text_width;
    end
    text(base_x, -0.1, ['Truth in prior bin ', num2str(ens_rank)],...
        'FontSize', 16, 'FontWeight', 'bold');
    % Draw an arrow from the label string to the truth
    % The existing DART arrow utility won't work here
    %%%arrow([base_x + text_width / 8, -0.08], [0, -0.03]);
    % Draw a line for now but come back for an arrow
    plot([base_x + text_width / 8, 0], [-0.08, -0.03], 'k');
    
    % Label this plot
    xlabel('State', 'Fontsize', 14);
    title('Latest ensemble prior', 'FontSize', 14);
    
    
    % Switch back to figure 1 for focus
    figure(1);
    
    
    % Update the permanent storage of the rank values
    handles.prior_rank(ens_rank) = handles.prior_rank(ens_rank) + 1;
    
    % Set the pushbutton to say Assimilate Obs
    set(handles.pushbutton1, 'String', 'Assimilate Obs');
    
    % Update the global storage of the ensemble
    handles.ens = ens_new;
    
else
    % Ready to do an assimilation
    % Next step should be an advance
    handles.ready_to_advance = true;
    % Generate the observation as a draw Normal(0, 1)
    obs_error_sd = handles.obs_error_sd;
    observation = randn(obs_error_sd);
    % Plot the observation
    plot(handles.time_step, observation, 'r*', 'MarkerSize', 10);
    
    % Set a legend
    %%%legend('Prior', 'Posterior', 'Observation');
    
    % Set the pushbutton to say Advance Model
    set(handles.pushbutton1, 'String', 'Advance Model');
    
    % Adjust the horizontal range of the plot window as needed
    if( mod(handles.time_step, 5) == 0)
        subplot(handles.r1);
        % Adjust the x axis limits to shift periodically
        set(gca, 'XLim', [handles.time_step - 4, handles.time_step + 6]);
        
        subplot(handles.r2);
        % Adjust the x axis limits to shift periodically
        set(gca, 'XLim', [handles.time_step - 4, handles.time_step + 6]);
        % Want the lower y limit to stay 0 for error spread
        set(gca, 'YLimMode', 'Auto');
        y_lim = get(gca, 'YLim');
        y_lim(1) = 0;
        set(gca, 'YLim', y_lim);
        
        subplot(handles.r3);
        % Adjust the x axis limits to shift periodically
        set(gca, 'XLim', [handles.time_step - 4, handles.time_step + 6]);
        % Want the lower y limit to stay 0 for error spread
        set(gca, 'YLimMode', 'Auto');
        y_lim = get(gca, 'YLim');
        y_lim(1) = 0;
        set(gca, 'YLim', y_lim);
    end
    
    % Do the assimilation
    ens = handles.ens;
    obs_error_sd = handles.obs_error_sd;
    
    % Figure out which filter option is currently selected
    h_filter_kind = get(handles.popupmenu1);
    
    filter_type = char(h_filter_kind.String(h_filter_kind.Value));
    
    switch filter_type
        case 'EAKF'
            [obs_increments, ~] = ...
                obs_increment_eakf(ens, observation, obs_error_sd^2);
        case 'EnKF'
            [obs_increments, ~] = ...
                obs_increment_enkf(ens, observation, obs_error_sd^2);
        case 'RHF'
            [obs_increments, ~] = ...
                obs_increment_rhf(ens, observation, obs_error_sd^2);
    end
    
    new_ens = ens + obs_increments;
    subplot(handles.r1);
    plot(handles.time_step + 0.1, new_ens, 'b*', 'MarkerSize', 6);
    handles.ens = new_ens;
    % Need axis info for the gui plot
    ens_axis = axis;
    
    % Update the rank data
    subplot(handles.r5);
    ens_rank = get_ens_rank(handles.ens, 0);
    
    % Plot the latest rank entry as a different color
    temp_rank(:, 1) = handles.posterior_rank(1:handles.ens_size + 1);
    temp_rank(:, 2) = 0;
    temp_rank(ens_rank, 2) = 1;
    
    bar(temp_rank, 'stacked');
    
    ylabel('Frequency');
    xlabel('Rank');
    title 'Posterior Rank Histogram';
    axis tight;
    
    % Update the permanent storage of the rank values
    handles.posterior_rank(ens_rank) = handles.posterior_rank(ens_rank) + 1;
    
    % Plot the segment for the updated error
    post_error = calculate_rmse(new_ens, 0.0);
    subplot(handles.r2);
    plot([handles.time_step - 0.1, handles.time_step + 0.1], ...
        [handles.error, post_error],'LineWidth',2.0);
    handles.error = post_error;
    
    % Plot the segment for the updated spread
    post_spread = std(new_ens);
    subplot(handles.r2);
    plot([handles.time_step - 0.1, handles.time_step + 0.1], ...
        [handles.spread, post_spread], 'r','LineWidth',2.0);
    
    % Want the lower y limit to stay 0 for error spread
    set(gca, 'YLimMode', 'Auto');
    y_lim = get(gca, 'YLim');
    y_lim(1) = 0;
    set(gca, 'YLim', y_lim);
    
    handles.spread= post_spread;
    
    % Plot the segment for the updated kurtosis
    post_kurtosis= kurt(new_ens);
    subplot(handles.r3);
    plot([handles.time_step - 0.1, handles.time_step + 0.1], ...
        [handles.kurtosis, post_kurtosis], 'r','LineWidth',2.0);
    
    % Want the lower y limit to stay 0 for error spread
    set(gca, 'YLimMode', 'Auto');
    y_lim = get(gca, 'YLim');
    y_lim(1) = 0;
    set(gca, 'YLim', y_lim);
    
    handles.kurtosis= post_kurtosis;
    
    % Plot the gui window for this update
    axes(handles.axes1);
    hold off
    
    % Plot an axis
    plot(ens_axis(3:4), [0 0], 'k', 'LineWidth', 2);
    hold on
    
    % Plot the prior ensemble members in green
    % Plotting ticks instead of asterisks makes bins clearer
    tick_half = 0.015;
    for n_tick = 1:handles.ens_size
        hold on;
        hg_prior = plot([ens(n_tick), ens(n_tick)], ...
            [-tick_half, tick_half], 'Color', [0 0.73 0], ...
            'LineWidth', 2);
    end
    
    % Plot the posterior ensemble members in blue
    for n_tick = 1:handles.ens_size
        hg_post = plot([new_ens(n_tick), new_ens(n_tick)],...
            [-0.1 - tick_half, -0.1 + tick_half], 'b', 'LineWidth', 2);
    end
    
    % Plot the observation likelihood
    hg_like = plot_gaussian(observation, obs_error_sd, 1.0);
    set(hg_like, 'color', 'red', 'LineWidth', 2,'LineStyle','--');
    
    % Plot the observation (at 0) as an asterisk
    plot(observation, 0, 'r*', 'MarkerSize', 14, 'LineWidth', 2.0);
    
    % Plot the truth (at -0.1 and 0) as a tick
    tick_half = 0.02;
    plot([0 0], [-0.1 - tick_half, -0.1 + tick_half], ...
        'k', 'LineWidth', 2.0);
    plot([0 0], [- tick_half, tick_half], ...
        'k', 'LineWidth', 2.0);
    
    % Put in some information about the bin, x position tricky
    base_x = max(0, ens_axis(3));
    text_width = (ens_axis(4) - ens_axis(3)) / 3;
    if((base_x + text_width) > ens_axis(4))
        base_x = ens_axis(4) - text_width;
    end
    text(base_x, -0.18, ['Truth in posterior bin ', num2str(ens_rank)],...
        'FontSize', 16, 'FontWeight', 'bold');
    % Draw an arrow from the label string to the truth
    % The existing DART arrow utility won't work here
    %%%arrow([base_x + text_width / 8, -0.16], [0, -0.13]);
    % Draw a line for now but come back for an arrow
    plot([base_x + text_width / 8, 0], [-0.16, -0.13], 'k');
    
    % Fix up the final axis
    def_axis = axis;
    def_axis(1:2) = ens_axis(3:4);
    def_axis(3) = -0.2;
    def_axis(4) = 1 / (sqrt(2 * pi) * handles.obs_error_sd) + 0.02;
    axis(def_axis);
    
    % Turn off the negative labels
    set(gca, 'YTick', [0 0.1 0.2 0.3 0.4]);
    set(gca, 'YTickLabel', [0 0.1 0.2 0.3 0.4]);
    grid on;
    
    % Plot an additional axis
    plot(ens_axis(3:4), [-0.1 -0.1], 'k', 'LineWidth', 2);
    hold on
    
    % Label this plot
    xlabel('State', 'Fontsize', 14);
    title('Latest ensemble prior, likelihood, posterior', 'FontSize', 14);
    
    % Put on legend
    legend([hg_prior, hg_post, hg_like], 'Prior', 'Posterior', 'Likelihood');
    
    % Return to the general figure window
    figure(1);
    
end

% Reset focus to the menu gui window
% Setting the axes clears the legend, gcbo restores focus
%%%axes(handles.axes1);
[~, gcbo_fig] = gcbo;
figure(gcbo_fig);

% Update handles structure
guidata(hObject, handles);



function edit7_Callback(hObject, ~, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double

% Get the value of the model_bias
handles.inflation= str2double(get(hObject, 'String'));
if(not(isfinite(handles.inflation)) ||  handles.inflation < 1)
    % Indicate input error in text box
    set(handles.edit7, 'String', '???');
    
    % After this, only this edit box will work
    turn_off_controls(handles);
    set(handles.edit7, 'Enable', 'On');
    
    return
end

% Turn on all the controls
turn_on_controls(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, ~, ~)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%-----------Turns off all the controls----------
function turn_off_controls(handles)

% Turn off all the other controls to avoid a mess
set(handles.pushbutton1,      'Enable', 'Off');
set(handles.pushbutton_run,   'Enable', 'Off');
set(handles.reset_pushbutton, 'Enable', 'Off');
set(handles.edit1,            'Enable', 'Off');
set(handles.edit5,            'Enable', 'Off');
set(handles.edit6,            'Enable', 'Off');
set(handles.edit7,            'Enable', 'Off');
set(handles.popupmenu1,       'Enable', 'Off');




%-----------Turns on all the controls----------
function turn_on_controls(handles)

% Turn on all the other controls to avoid a mess
set(handles.pushbutton1,      'Enable', 'On');
set(handles.pushbutton_run,   'Enable', 'On');
set(handles.reset_pushbutton, 'Enable', 'On');
set(handles.edit1,            'Enable', 'On');
set(handles.edit5,            'Enable', 'On');
set(handles.edit6,            'Enable', 'On');
set(handles.edit7,            'Enable', 'On');
set(handles.popupmenu1,       'Enable', 'On');


%----------- calculates rmse ----------
function y = calculate_rmse(x, truth)

squared_error = (x - truth).^2;
y = sqrt(mean(squared_error));

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

