function varargout = run_lorenz_96(varargin)
%% RUN_LORENZ_96 ensemble data assimilation with a 40-variable implementation
%      of the Lorenz '96 dynamical model.
%
%      To demonstrate the analogue to the atmosphere, the model is a cyclic 
%      1D domain with equally-spaced nodal points. There are 20 ensemble 
%      members initially in this example.
%
%      The model can be single-stepped through model advance and assimilation
%      steps using the top pushbutton, or allowed to run free using the
%      'Start Free Run' button. A variety of assimilation algorithms can
%      be selected from the first pulldown. Model error in the assimilating
%      model (an imperfect model assimilation experiment) can be selected
%      with the second pulldown. The localization, inflation and ensemble 
%      size can be changed with the three dialogue boxes. Changing the 
%      ensemble size resets the diagnostic displays. The figure window
%      displays time sequences of the prior and posterior error and prior
%      and posterior (if assimilation is on) rank histograms.
%      
%      It takes about twenty timesteps for the intially small ensemble
%      perturbations to grow large enough to be seen using the default
%      settings.
% 
% See also: gaussian_product, oned_model, oned_ensemble, twod_ensemble, 
%           run_lorenz_63

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @run_lorenz_96_OpeningFcn, ...
                   'gui_OutputFcn',  @run_lorenz_96_OutputFcn, ...
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


% --- Executes just before run_lorenz_96 is made visible.
function run_lorenz_96_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to run_lorenz_96 (see VARARGIN)

help run_lorenz_96

% Choose default command line output for run_lorenz_96
handles.output = hObject;

% Global semaphore; ready to advance or assimilate?
handles.ready_to_advance = true;

% Initialize the L96 model
lorenz_96_static_init_model;

% Initialize global storage for control trajectory and ensembles

global FORCING;
global MODEL_SIZE;

handles.ens_size = str2double(get(handles.edit_ens_size, 'String'));
handles.inflation= str2double(get(handles.edit_inflation, 'String'));
handles.localization= str2double(get(handles.edit_localization, 'String'));
handles.model_size = MODEL_SIZE;
handles.true_state(1, 1:MODEL_SIZE) = FORCING;
handles.true_state(1, 1) = 1.001 * FORCING;
handles.time = 1;
handles.prior_rms(1) = 0;
handles.posterior_rms(1) = 0;
% Used to normalize the polar plotting
handles.mean_dist = 35;
handles.h_ens = 0;
handles.h_truth = 0;

% Generate set of ensemble perturbations
for n = 1:handles.ens_size
   handles.post(1, 1:MODEL_SIZE, n) = handles.true_state(1, :);
end
handles.post(1, 1:MODEL_SIZE, 1:handles.ens_size) = ...
   handles.post(1, 1:MODEL_SIZE, 1:handles.ens_size) + ...
   0.001 * randn(1, MODEL_SIZE, handles.ens_size);

% For convenience make the first prior identical to the first posterior
handles.prior(1, 1:MODEL_SIZE, 1:handles.ens_size) = handles.post;

% An array to keep track of rank histograms
handles.prior_rank(1 : handles.ens_size + 1) = 0;
handles.posterior_rank(1 : handles.ens_size + 1) = 0;

% Handles to subregions of plot
handles.r1 = 0;
handles.r2 = 0;
handles.r3 = 0;

% Set up subdomains in figure window for timeseries and rank historgrams
figure(1);
handles.r1 = subplot(2, 1, 1);
% Want the y axis limits to take care of themselves
set(gca, 'YLimMode', 'Auto');

% Get a legend on here from the beginning
plot(handles.prior_rms, 'g');
hold on
plot(handles.posterior_rms, 'b');
legend('Prior', 'Posterior', 'Location', 'NorthWest')
ylabel('RMS Error', 'FontSize', 14);
xlabel('Time ', 'FontSize', 14);


handles.r2 = subplot(2, 2, 3);
ylabel('Frequency');
xlabel('Rank');
title ('Prior Rank Histogram', 'FontSize', 14);
hold on

handles.r3 = subplot(2, 2, 4);
ylabel('Frequency');
xlabel('Rank');
title ('Posterior Rank Histogram', 'FontSize', 14);
hold on

% Select the first plotting box; make the thing look polar
axes(handles.axes1);

hold off
% Plot a single point to make sure the axis limits are okay
% Unclear if these can be set more cleanly with polar
% This also gets the observations into the legend
x(1) = 14.9;
y_2 = [0, 2*pi];
h = plot_polar(y_2, x, handles.mean_dist, 'r*', 1);
hold on
set(h, 'Visible', 'Off');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes run_lorenz_96 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = run_lorenz_96_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



%%% SETTINGS REQUIRED FOR PUSHBUTTONS IN USE WITH THIS SCRIPT
% The single_step pushbutton requires:
% Busy action:   queue
% Enable:        on
% Interruptible: off
% Units:         normalized  (required for resizing)

% The free_run pushbutton requires:
% Busy action:   queue
% Enable:        on
% Interruptible: on
% Units:         normalized  (required for resizing)


% --- Executes on button press in pushbutton_single_step.
function pushbutton_single_step_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_single_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Move the model ahead a step or assimilate observations as appropriate
step_ahead(hObject, handles);



% --- Executes on button press in pushbutton_free_run.
function pushbutton_free_run_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_free_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Turn off all the other model status controls to avoid a mess
turn_off_controls(handles);
set(handles.pushbutton_free_run, 'Enable', 'On');

% Check the button label to see if we are starting or stopping a free run
if(strcmp(get(hObject, 'String'), 'Stop Free Run'))

   % Turn off the free run pushbutton until everything has completely stopped
   set(hObject, 'Enable', 'Off');

   % Being told to stop; switch to not running status
   set(hObject, 'String', 'Start Free Run');
 
   % Update the handles global structure
   guidata(hObject, handles);

else
   % Being told to start free run
   % Change the pushbutton to stop
   set(hObject, 'String', 'Stop Free Run');
   % Update the handles global structure
   guidata(hObject, handles);

   % Loop through advance and assimilate steps until stopped
   while(true)
      % Check to see if stop has been pushed; get latest copy of global data
      my_data = guidata(gcbo);

      status_string = get(my_data.pushbutton_free_run, 'String');
      if(strcmp(status_string, 'Start Free Run'))

         % Turn all the other model status controls back on
         turn_on_controls(handles);

         % Very last, turn on the start free run button
         %%%set(hObject, 'Enable', 'On'); 
 
         return
      end
      % Do the next advance or assimilation step
      step_ahead(hObject, my_data)
   end
end

% Only way to get here is at the end of a stop; No need to clean up
% Since the next thing that will happen is a test for stop in the free
% run loop.



%----------- Moves the model ahead or assimilates next observations ------
function step_ahead(hObject, handles)

% Test on semaphore, either advance or assimilate
if(handles.ready_to_advance)
   % Set semaphore to indicate that next step is an assimilation
   handles.ready_to_advance = false;

   % Set the pushbutton text to indicate that next step is an assimilate
   set(handles.pushbutton_single_step, 'String', 'Assimilate Obs');

   % Code for advancing model comes next (delete two exisiting lines)
   time = handles.time;
   [new_truth, new_time] = lorenz_96_adv_1step(handles.true_state(time, :), time);
   handles.time = new_time;
   handles.true_state(new_time, :) = new_truth;

   % See if ensembles should have model forcing error
   h_model_error= get(handles.popupmenu_model_error);
   h_err_index = h_model_error.Value;

   global FORCING
   if(h_err_index == 2) FORCING = 6; end

   %  Advance the ensemble members; posterior -> new prior
   for n = 1:handles.ens_size
      [new_ens, new_time] = lorenz_96_adv_1step(handles.post(time, :, n), time);
      handles.prior(new_time, :, n) = new_ens;
   end

   % Reset the forcing to 8
   FORCING = 8;

   % Inflate ensemble
   global MODEL_SIZE;
   for i = 1:MODEL_SIZE
      ens_mean = mean(handles.prior(new_time, i, :));
      handles.prior(new_time, i, :) = ens_mean + ...
         sqrt(handles.inflation) * (handles.prior(new_time, i, :) - ens_mean); 
   end




   % Plot the true state and the ensembles on polar plot
   y = (0:MODEL_SIZE) / MODEL_SIZE * 2 * pi;

   % Select the first plotting box
   axes(handles.axes1);

   hold off
   % Plot a single point to make sure the axis limits are okay
   % Unclear if these can be set more cleanly with polar
   % This also gets the observations into the legend
   x(1) = 14.9;
   y_2 = [0, 2*pi];
   h = plot_polar(y_2, x, handles.mean_dist, 'r*', 1);
   hold on
   set(h, 'Visible', 'Off');


   for n = 1:handles.ens_size
      handles.h_ens = plot_polar(y, handles.prior(new_time, :, n), ...
         handles.mean_dist, 'g', MODEL_SIZE);
      hold on
   end
   handles.h_truth = plot_polar(y, new_truth, handles.mean_dist, 'k', MODEL_SIZE);
   % Truth is in black

   % Get a legend shifted outside the plot
   h_leg = legend([handles.h_truth handles.h_ens, h], ...
      'True State', 'Ensemble', 'Observations', 'Location', 'NorthEast');
   pos = get(h_leg, 'Position');
   pos(1) = pos(1) + 0.1;
   set(h_leg, 'Position', pos);

   % Update the time label
   set(handles.text_time, 'String', ['Time = ', num2str(new_time)]);

   % Compute the prior RMS error of ensemble mean
   prior_rms = rms_error(new_truth, handles.prior(new_time, :, :));
   handles.prior_rms(new_time) = prior_rms;

   % Save the information about the histograms from before
   temp_rank(:, 1) = handles.prior_rank(1:handles.ens_size + 1);
   temp_rank(:, 2) = 0;

   % Compute the prior rank histograms
   for i = 1:handles.ens_size
      ens_rank = get_ens_rank(squeeze(handles.prior(new_time, i, :)), squeeze(new_truth(i)));
      handles.prior_rank(ens_rank) = handles.prior_rank(ens_rank) + 1;
      temp_rank(ens_rank, 2) = temp_rank(ens_rank, 2) + 1;
   end
   
% Plot the prior_rms error time series
   figure(1);
   subplot(handles.r1); 
   hold on
   plot(handles.prior_rms, 'g');

   % Plot the rank histogram for the prior
   subplot(handles.r2);
   bar(temp_rank, 'stacked');
   axis tight;

   % Select the first plotting box
   axes(handles.axes1);

else
   % Set semaphore to indicate that next step is a model advance
   handles.ready_to_advance = true;

   % Set the pushbutton text to indicate that the next step is a model advance
   set(handles.pushbutton_single_step, 'String', 'Advance Model');

   % Get current time step
   time = handles.time;

   % Determine what type of assimilation is being done (none, EAKF, EnKF, RHF)
   h_filter_kind = get(handles.popupmenu_assim_type);
   filter_type = char(h_filter_kind.String(h_filter_kind.Value));

   if(strcmp(filter_type, 'No Assimilation'))
      % Just copy prior to posterior
      handles.post(time, :, :) = handles.prior(time, :, :);
   else
      % Code for doing the assimilation comes here


%-----------------
      % Generate noisy observations of the truth
      obs_sd = 4;
      obs_error_var = obs_sd^2;

      % Do fully sequential assimilation algorithm
      temp_ens = squeeze(handles.prior(time, :, :));
   
      % Select the first plotting box
      axes(handles.axes1);

      % Observe each state variable independently
      global MODEL_SIZE
      for i = 1:MODEL_SIZE
         obs_prior = temp_ens(i, :);
         obs(i) = handles.true_state(time, i) + obs_sd * randn;

         % Compute the increments for observed variable
         switch filter_type
            case 'EAKF'
               [obs_increments, err] = ...
                  obs_increment_eakf(obs_prior, obs(i), obs_error_var);
            case 'EnKF'
               [obs_increments, err] = ...
                  obs_increment_enkf(obs_prior, obs(i), obs_error_var);
            case 'RHF'
               [obs_increments, err] = ...
                  obs_increment_rhf(obs_prior, obs(i), obs_error_var);
         end

         % Regress the increments onto each of the state variables
         for j = 1:MODEL_SIZE
            state_incs = get_state_increments(temp_ens(j, :), ...
               obs_prior, obs_increments);
            
            % Compute distance between obs and state for localization
            dist = abs(i - j) / 40;
            if(dist > 0.5) dist = 1 - dist; end

            % Compute the localization factor
            cov_factor = comp_cov_factor(dist, handles.localization);

            temp_ens(j, :) = temp_ens(j, :) + state_incs * cov_factor;
         end
      end

      % Plot the observations
      y = (0:MODEL_SIZE) / MODEL_SIZE * 2 * pi;
      h_obs = plot_polar(y, obs, handles.mean_dist, 'r*', MODEL_SIZE);

      % Update the posterior
      handles.post(time, :, :) = temp_ens;

      % Compute the posterior rms
      posterior_rms = rms_error(handles.true_state(time, :), handles.post(time, :, :));
      handles.posterior_rms(time) = posterior_rms;

      % Save the information about the histograms from before
      temp_rank(:, 1) = handles.posterior_rank(1:handles.ens_size + 1);
      temp_rank(:, 2) = 0;
    
      % Compute the posterior rank histograms
      for i = 1:handles.ens_size
         ens_rank = get_ens_rank(squeeze(handles.post(time, i, :)), ...
            squeeze(handles.true_state(time, i)));
         handles.posterior_rank(ens_rank) = handles.posterior_rank(ens_rank) + 1;
         temp_rank(ens_rank, 2) = temp_rank(ens_rank, 2) + 1;
      end
   
      % Plot the posterior_rms error time series
      figure(1);
      subplot(handles.r1);
      hold on
      plot(handles.posterior_rms, 'b');

      % Plot the rank histogram for the prior
      subplot(handles.r3);
      bar(temp_rank, 'stacked');
      axis tight;

      % Select the first plotting box
      axes(handles.axes1);

      pause(0.2)

      

%-----------------
   end

end

% If using multiple windows might need to reset focus to the gui window here
[gcbo_h, gcbo_fig] = gcbo;
figure(gcbo_fig);

% Update the global storage and return
guidata(hObject, handles);


% --- Executes on selection change in popupmenu_assim_type.
function popupmenu_assim_type_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_assim_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_assim_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_assim_type


% --- Executes during object creation, after setting all properties.
function popupmenu_assim_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_assim_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in popupmenu_model_error.
function popupmenu_model_error_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_model_error (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_model_error contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_model_error


% --- Executes during object creation, after setting all properties.
function popupmenu_model_error_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_model_error (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_ens_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ens_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ens_size as text
%        str2double(get(hObject,'String')) returns contents of edit_ens_size as a double

% Set the ensemble size global value to the update
handles.ens_size = str2double(get(hObject, 'String'));
if(not(isfinite(handles.ens_size)) | handles.ens_size < 2) 
   set(handles.edit_ens_size, 'String', '???');
   
   % After this, only this edit box will work
   turn_off_controls(handles);
   set(handles.edit_ens_size, 'Enable', 'On');

   return
end

% Enable all controls
turn_on_controls(handles);

% Need to reset the ensemble and the time
global FORCING
handles.true_state(1, 1:handles.model_size) = FORCING;
handles.true_state(1, 1) = 1.001 * FORCING;
handles.time = 1;

% Generate set of ensemble perturbations
handles.post = [0];
for n = 1:handles.ens_size
   handles.post(1, 1:handles.model_size, n) = handles.true_state(1, :);
end
handles.post(1, 1:handles.model_size, 1:handles.ens_size) = ...
   handles.post(1, 1:handles.model_size, 1:handles.ens_size) + ...
   0.001 * randn(1, handles.model_size, handles.ens_size);

% For convenience make the first prior identical to the first posterior
handles.prior = [0];
handles.prior(1, 1:handles.model_size, 1:handles.ens_size) = ...
   handles.post(1, 1:handles.model_size, 1:handles.ens_size);

% Clear the rms plots
handles.prior_rms = [0];
handles.posterior_rms = [0];

% Reset the array to keep track of rank histograms
handles.prior_rank(1 : handles.ens_size + 1) = 0;
handles.posterior_rank(1 : handles.ens_size + 1) = 0;

% Select the plot
figure(1);
handles.r1 = subplot(2, 1, 1);

% Clear plot and start over
hold off
handles.prior_rms
plot(handles.prior_rms, 'g');
hold on
plot(handles.posterior_rms, 'b');
ylabel('RMS Error', 'FontSize', 14);
xlabel('Time ', 'FontSize', 14);

% Get a legend on here from the beginning
legend('Prior', 'Posterior', 'Location', 'NorthWest')

% Reset the rank histograms to the new size
handles.prior_rank = [1 : handles.ens_size + 1] * 0;
handles.posterior_rank = [1 : handles.ens_size + 1] * 0;
subplot(handles.r2);
hold off
bar(handles.prior_rank);
ylabel('Frequency', 'FontSize', 14);
xlabel('Rank', 'FontSize', 14);
title('Prior Rank Histogram', 'FontSize', 14);
hold on
axis tight;

subplot(handles.r3);
hold off
bar(handles.posterior_rank);
hold on
ylabel('Frequency');
xlabel('Rank');
title('Posterior Rank Histogram', 'FontSize', 14);
axis tight;


% Switch back to the gui window
axes(handles.axes1);


% Update handles structure
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function edit_ens_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ens_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when uipanel2 is resized.
function uipanel2_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to uipanel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




function turn_off_controls(handles)

set(handles.pushbutton_single_step,     'Enable', 'Off');
set(handles.pushbutton_free_run,        'Enable', 'Off');
set(handles.popupmenu_assim_type,       'Enable', 'Off');
set(handles.popupmenu_model_error,      'Enable', 'Off');
set(handles.edit_localization,          'Enable', 'Off');
set(handles.edit_inflation,             'Enable', 'Off');
set(handles.edit_ens_size,              'Enable', 'Off');




function turn_on_controls(handles)

set(handles.pushbutton_single_step,     'Enable', 'On');
set(handles.pushbutton_free_run,        'Enable', 'On');
set(handles.popupmenu_assim_type,       'Enable', 'On');
set(handles.popupmenu_model_error,      'Enable', 'On');
set(handles.edit_localization,          'Enable', 'On');
set(handles.edit_inflation,             'Enable', 'On');
set(handles.edit_ens_size,              'Enable', 'On');



function edit_inflation_Callback(hObject, eventdata, handles)
% hObject    handle to edit_inflation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_inflation as text
%        str2double(get(hObject,'String')) returns contents of edit_inflation as a double

% Set the inflation value to the update
handles.inflation= str2double(get(hObject, 'String'));
if(not(isfinite(handles.inflation)) | handles.inflation< 0) 
   set(handles.edit_inflation, 'String', '???');
   
   % After this, only this edit box will work
   turn_off_controls(handles);
   set(handles.edit_inflation, 'Enable', 'On');

   return
end

% Enable all controls
turn_on_controls(handles);

% Update handles structure
guidata(hObject, handles);





% --- Executes during object creation, after setting all properties.
function edit_inflation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_inflation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_localization_Callback(hObject, eventdata, handles)
% hObject    handle to edit_localization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_localization as text
%        str2double(get(hObject,'String')) returns contents of edit_localization as a double

% Set the localization value to the update
handles.localization= str2double(get(hObject, 'String'));
if(not(isfinite(handles.localization)) | handles.localization< 0) 
   set(handles.edit_localization, 'String', '???');
   
   % After this, only this edit box will work
   turn_off_controls(handles);
   set(handles.edit_localization, 'Enable', 'On');

   return
end

% Enable all controls
turn_on_controls(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function edit_localization_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_localization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function ens_mean_rms = rms_error(truth, ens)

ens_mean = mean(squeeze(ens)');
ens_mean_rms = sqrt(sum((truth - ens_mean).^2) / size(truth, 2));

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

