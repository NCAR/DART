function varargout = run_lorenz_96(varargin)
% RUN_LORENZ_96 ensemble data assimilation with a 40-variable implementation
%      of the Lorenz '96 dynamical model.
%
%      To demonstrate the analogue to the atmosphere, the model is a cyclic 
%      1D domain with equally-spaced nodal points. There are 20 ensemble 
%      members in this example, each with 0.1% noise from a random normal 
%      distribution.
%
%      There are several experiments to perform - including generating the
%      true state with one forcing and assimilaing with a model that has 
%      the wrong forcing. If there's such a thing as a perfect model 
%      experiment with an imperfect model, this is IT!
%
%      This utility also explores the effect of Localization - the ability 
%      to restrict the impact of an observation to a subset of the state 
%      vector.
% 
% See also: gaussian_product, oned_model, oned_ensemble, twod_ensemble, 
%           run_lorenz_63

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

% Last Modified by GUIDE v2.5 27-Aug-2009 16:35:11

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

ens_size = 20;
handles.ens_size = ens_size;
handles.model_size = MODEL_SIZE;
handles.true_state(1, 1:MODEL_SIZE) = FORCING;
handles.true_state(1, 1) = 1.001 * FORCING;
handles.time = 1;
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
% MAKE SURE TO INCLUDE OTHER CONTROLS HERE
set(handles.pushbutton_single_step, 'Enable', 'Off');
set(handles.popupmenu_assim_type,   'Enable', 'Off');
set(handles.popupmenu_localization, 'Enable', 'Off');
set(handles.popupmenu_model_error,  'Enable', 'Off');
set(handles.popupmenu_inflation,    'Enable', 'Off');

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
         % MAKE SURE TO INCLUDE OTHER CONTROLS HERE
         set(handles.pushbutton_single_step, 'Enable', 'On');
         set(handles.popupmenu_assim_type,   'Enable', 'On');
         set(handles.popupmenu_localization, 'Enable', 'On');
         set(handles.popupmenu_model_error,  'Enable', 'On');
         set(handles.popupmenu_inflation,    'Enable', 'On');

         % Very last, turn on the start free run button
         set(hObject, 'Enable', 'On'); 
 
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

   % See if inflation is to be applied
   h_inflation = get(handles.popupmenu_inflation);
   h_inflation_index = h_inflation.Value;

   % If index is 2, inflate ensemble
   global MODEL_SIZE;
   if(h_inflation_index == 2)
      inflation = sqrt(1.10);
      for i = 1:MODEL_SIZE
         ens_mean = mean(handles.prior(new_time, i, :));
         handles.prior(new_time, i, :) = ens_mean + ...
            inflation * (handles.prior(new_time, i, :) - ens_mean); 
      end
   end




   % Plot the true state and the ensembles on polar plot
   y = (0:MODEL_SIZE) / MODEL_SIZE * 2 * pi;

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

      % Get localization from pull-down (large or 0.3 for now)
      h_localization = get(handles.popupmenu_localization);
      h_loc_index = h_localization.Value;
      if(h_loc_index == 1)
         localization_half_width = 10000.0;
      else
         localization_half_width = 0.3;
      end

%-----------------
      % Generate noisy observations of the truth
      obs_sd = 4;
      obs_error_var = obs_sd^2;

      % Do fully sequential assimilation algorithm
      temp_ens = squeeze(handles.prior(time, :, :));

      % Observe each state variable independently
      global MODEL_SIZE
      for i = 1:MODEL_SIZE
         obs_prior = temp_ens(i, :);
         obs(i) = handles.true_state(time, i) + obs_sd * randn;
         % Compute the increments for observed variable
         [obs_increments, err] = obs_increment_eakf(obs_prior, obs(i), obs_error_var);

         % Regress the increments onto each of the state variables
         for j = 1:MODEL_SIZE
            state_incs = get_state_increments(temp_ens(j, :), ...
               obs_prior, obs_increments);
            
            % Compute distance between obs and state for localization
            dist = abs(i - j) / 40;
            if(dist > 0.5) dist = 1 - dist; end

            % Compute the localization factor
            cov_factor = comp_cov_factor(dist, localization_half_width);

            temp_ens(j, :) = temp_ens(j, :) + state_incs * cov_factor;
         end
      end

      % Plot the observations
      y = (0:MODEL_SIZE) / MODEL_SIZE * 2 * pi;
      h_obs = plot_polar(y, obs, handles.mean_dist, 'r*', MODEL_SIZE);

      % Update the posterior
      handles.post(time, :, :) = temp_ens;
      pause(0.5)

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


% --- Executes on selection change in popupmenu_localization.
function popupmenu_localization_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_localization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_localization contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_localization


% --- Executes during object creation, after setting all properties.
function popupmenu_localization_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_localization (see GCBO)
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


% --- Executes on selection change in popupmenu_inflation.
function popupmenu_inflation_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_inflation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_inflation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_inflation


% --- Executes during object creation, after setting all properties.
function popupmenu_inflation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_inflation (see GCBO)
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
   %%%turn_off_controls(handles);
   set(handles.edit_ens_size, 'Enable', 'On');

   return
end

% Enable all controls
%%%turn_on_controls(handles);

% Need to reset the ensemble and the time
slkjdf



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


