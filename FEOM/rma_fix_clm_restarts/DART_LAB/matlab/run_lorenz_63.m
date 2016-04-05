function varargout = run_lorenz_63(varargin)
%% RUN_LORENZ_63 ensemble data assimilation with the 3-variable
%      Lorenz '63 dynamical model - the "butterfly" model.
%
%      There are 20 ensemble members in this example. The initial
%      conditions are chosen to provide an interesting trajectory.
%
%      This is another 'perfect_model' experiment. As the model is
%      advanced, observations are taken from the true state and are
%      assimilated by the ensemble members. The true trajectory, the
%      observations, the observation increments, and the Prior and Posterior
%      ensemble states are displayed.
%
%      To provide context and highlight the details of the assimilation,
%      two views are presented. The larger view is a detailed view of
%      the model space in the immediate vicinity of the True State. The
%      smaller view provides the context of the entire model space.
%      Both views are fundamentally views 'from above' ... looking down
%      on the Z-axis although they can be rotated when the assimilation is
%      stopped.
%
%      After you get the feel for a few single steps through the process
%      (by repeatedly pressing the 'Advance/Assimilate' button), select
%      'No Assimilation' and start a free run. After a while some of the
%      ensemble members wind up on the opposite lobe of the attractor.
%
%      Restart the experiment with some sort of assimilation and watch
%      the ensemble members diverge from the True State, only to get
%      nudged back to the True State by the observations. This is
%      particularly difficult - especially in the highly nonlinear
%      'saddle' region.
%
% See also: gaussian_product, oned_model, oned_ensemble, twod_ensemble,
%           run_lorenz_96

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @run_lorenz_63_OpeningFcn, ...
                   'gui_OutputFcn',  @run_lorenz_63_OutputFcn, ...
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



function run_lorenz_63_OpeningFcn(hObject, ~, handles, varargin)
%% --- Executes just before run_lorenz_63 is made visible.
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to run_lorenz_63 (see VARARGIN)

help run_lorenz_63

% set random number seed to same value to generate known sequences
rng('default')

% Choose default command line output for run_lorenz_63
handles.output = hObject;

% Global semaphore; ready to advance or assimilate?
handles.ready_to_advance = true;

% Initialize the l63 model
lorenz_63_static_init_model;

% Initialize global storage for control trajectory and ensemble
ens_size = 20;
handles.ens_size = ens_size;
handles.model_size = 3;
handles.true_state(1, 1:3) = [10.2471, 2.8443, 36.2666];
handles.time = 1;
handles.hstar = 0;
handles.obs(1:3) = 0;
handles.obs_sd = 1;
handles.obs_error_var = handles.obs_sd^2;
handles.h_global_obs = 0;
handles.h_global_init = false;

% Generate a set of perturbations for the first posterior ensembles
for n = 1:handles.ens_size
   handles.post(1, :, n) = handles.true_state(1, :);
end
handles.post(1, 1:3, 1:ens_size) = handles.post(1, 1:3, 1:ens_size) + ...
   0.1 * randn(1, 3, ens_size);

% For convenience, make the first prior identical to the first posterior
handles.prior(1, 1:3, 1:ens_size) = handles.post(1, 1:3, 1:ens_size);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes run_lorenz_63 wait for user response (see UIRESUME)
% uiwait(handles.figure1);



function varargout = run_lorenz_63_OutputFcn(~, ~, handles)
%% --- Outputs from this function are returned to the command line.
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



function pushbutton_single_step_Callback(hObject, ~, handles)
%% --- Executes on button press in pushbutton_single_step.
% hObject    handle to pushbutton_single_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Move the model ahead a step or assimilate observations as appropriate
step_ahead(hObject, handles);



function pushbutton_free_run_Callback(hObject, ~, handles)
%% --- Executes on button press in pushbutton_free_run.
% hObject    handle to pushbutton_free_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Turn off all the other model status controls to avoid a mess
% MAKE SURE TO INCLUDE OTHER CONTROLS HERE
set(handles.pushbutton_single_step, 'Enable', 'Off');
set(handles.popupmenu_assim_type,   'Enable', 'Off');

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



function step_ahead(hObject, handles)
%%----------- Moves the model ahead or assimilates next observations ------


axes(handles.axes1)

% Test on semaphore, either advance or assimilate
if(handles.ready_to_advance)
   % Set semaphore to indicate that next step is an assimilation
   handles.ready_to_advance = false;

   % Set the pushbutton text to indicate that next step is an assimilate
   set(handles.pushbutton_single_step, 'String', 'Assimilate Obs');

   % Turn off the most recent global space observation plot if it exists
   if(handles.h_global_init)
      set(handles.h_global_obs, 'Visible', 'off');
   end

   % Advance a number of steps in between each assimilation
   num_steps_to_advance = 20;
   for i = 1:num_steps_to_advance

      % Code for advancing model comes next (delete two exisiting lines)
      time = handles.time;
      [new_truth, new_time] = lorenz_63_adv_1step(handles.true_state(time, :), time);
      handles.time = new_time;
      handles.true_state(new_time, :) = new_truth;

      % Advance the ensemble members; posterior -> new prior
      for n = 1:handles.ens_size
         [new_ens, new_time] = lorenz_63_adv_1step(handles.post(time, :, n), time);
         handles.prior(new_time, :, n) = new_ens;
      end

      % Plot a long trajectory of truth in small window for reference
      axes(handles.axes2);
      hold on;
      plot3(handles.true_state(new_time-1:new_time, 1), ...
            handles.true_state(new_time-1:new_time, 2), ...
            handles.true_state(new_time-1:new_time, 3), 'k');

      % Also plot an asterisk on the leading edge
      if(new_time > 2)
         set(handles.h_star, 'Visible', 'Off');
      end
      clear handles.h_star;
      handles.h_star = plot3(handles.true_state(new_time, 1), ...
                     handles.true_state(new_time, 2), ...
                     handles.true_state(new_time, 3), ...
                     'k*', 'MarkerSize', 16, 'LineWidth', 2);
      view([2 -1 1]);
      axis([-20 20 -20 20 5 45]);

      % Plot the close-up view of the ensmble
      axes(handles.axes1)

      % Plot the truth trajectory for the last 8 steps
      hold off;

      btime = new_time - 7;
      if(btime < 1), btime = 1; end
      plot3(handles.true_state(btime:new_time, 1), ...
            handles.true_state(btime:new_time, 2), ...
            handles.true_state(btime:new_time, 3), 'k', 'linewidth', 2);

      hold on
      % Set an appropriate consistent view angle
      view([2, -1 1]);

      % Plot an asterisk at the head of the trajectory
      plot3(handles.true_state(new_time, 1), ...
            handles.true_state(new_time, 2), ...
            handles.true_state(new_time, 3), 'k*', 'MarkerSize', 16, 'LineWidth', 2);

      % Adjust the axes to follow the truth
      xb = handles.true_state(new_time, 1);
      yb = handles.true_state(new_time, 2);
      zb = handles.true_state(new_time, 3);
      limits = [xb - 3,  xb + 3,  yb - 3,  yb + 3, zb - 3, zb + 3];
      axis(limits);


      % Plot the ensemble members advance trajectories, too
      for n = 1:handles.ens_size
         h = plot3(handles.prior(btime:new_time, 1, n), ...
               handles.prior(btime:new_time, 2, n), ...
               handles.prior(btime:new_time, 3, n), 'g');
         set(h,'Color',[0.0, 0.73, 0.0])
      end

      % Put on a legend
      %%%legend('True State', 'Ensembles');

      % Update the time label
      set(handles.text_time, 'String', ['Time = ', num2str(new_time)]);

      % Force the buffers to flush and plot the advance
      drawnow

      % Last prior update will get overwritten when assimilation is done
      handles.post(new_time, :, :) = handles.prior(new_time, :, :);
   end

   % Compute the observations for this time and save
   % NOTE: this cannot be inside the if statement that follows because
   % assimilation could be turned on after this advance!
   for i = 1:3
      handles.obs(i) = handles.true_state(new_time, i) + handles.obs_sd * randn;
   end

   % Plot the observation as a red asterisk if assimilating
   h_filter_kind = get(handles.popupmenu_assim_type);
   filter_type = char(h_filter_kind.String(h_filter_kind.Value));
   if(~strcmp(filter_type, 'No Assimilation'))

      plot3(handles.obs(1), handles.obs(2), handles.obs(3), 'r*', 'MarkerSize', 20);
      % Also plot observation in the global attractor plot
      axes(handles.axes2);
      handles.h_global_obs = ...
         plot3(handles.obs(1), handles.obs(2), handles.obs(3), 'r*', 'MarkerSize', 20);
      handles.h_global_init = true;
      axes(handles.axes1);

      pause(0.5)
   end

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

      % Generate noisy observations of the truth
      %%%obs_sd = 1;
      %%%obs_error_var = obs_sd^2;

      % Do fully sequential assimilation algorithm
      temp_ens = squeeze(handles.prior(time, :, :));

      % Observe each state variable independently
      for i = 1:3
         obs_prior = temp_ens(i, :);
         %%%obs(i) = handles.true_state(time, i) + obs_sd * randn;
         obs(i) = handles.obs(i);

         % Compute the increments for observed variable
         switch filter_type
            case 'EAKF'
               [obs_increments, ~] = ...
                  obs_increment_eakf(obs_prior, obs(i), handles.obs_error_var);
            case 'EnKF'
               [obs_increments, ~] = ...
                  obs_increment_enkf(obs_prior, obs(i), handles.obs_error_var);
            case 'RHF'
               [obs_increments, ~] = ...
                  obs_increment_rhf(obs_prior, obs(i), handles.obs_error_var);
         end


         %%%% Fully localized test
         %%%temp_ens(i, :) = temp_ens(i, :) + obs_increments;

         % Regress the increments onto each of the three state variables
         for j = 1:3
            state_incs = get_state_increments(temp_ens(j, :), ...
               obs_prior, obs_increments);
            temp_ens(j, :) = temp_ens(j, :) + state_incs;
         end
      end


      % Update the posterior
      handles.post(time, :, :) = temp_ens;

      % Plot a segment showing the impact of the observation
      for n = 1:handles.ens_size
         xup = [handles.prior(time, 1, n), handles.post(time, 1, n)];
         yup = [handles.prior(time, 2, n), handles.post(time, 2, n)];
         zup = [handles.prior(time, 3, n), handles.post(time, 3, n)];
         plot3(xup, yup, zup, 'r');
      end

      % Pause to allow a view of the observation impact
      pause(0.5)
   end
end

% If using multiple windows might need to reset focus to the gui window here
[~, gcbo_fig] = gcbo;
figure(gcbo_fig);

% Update the global storage and return
guidata(hObject, handles);



% --- Executes on selection change in popupmenu_assim_type.
function popupmenu_assim_type_Callback(~, ~, ~)
% hObject    handle to popupmenu_assim_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_assim_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_assim_type



% --- Executes during object creation, after setting all properties.
function popupmenu_assim_type_CreateFcn(hObject, ~, ~)
% hObject    handle to popupmenu_assim_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

