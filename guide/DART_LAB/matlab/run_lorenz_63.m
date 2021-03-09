function run_lorenz_63
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
% See also: gaussian_product.m oned_model.m oned_model_inf.m oned_ensemble.m
%           twod_ensemble.m run_lorenz_96.m run_lorenz_96_inf.m

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%>@ TODO FIXME if there is a way to preserve the view angle ...

help run_lorenz_63

global DELTAT
global SIGMA
global R
global B
global MODEL_SIZE

atts = stylesheet;  % get the default fonts and colors

%% Creates the figure layout

handles.figure1 = figure('Position', [450 250 720 630], ...
    'Name', 'run_lorenz_63', ...
    'Color', atts.background);

% Position has units of normalized and therefore the components must be
% Fractions of figure. Also, the actual text  has FontUnits which are
% normalized, so the Font Size is a fraction of the text box/ edit box that
% contains it. Making the font size and box size normalized allows the text
% size and box size to change proportionately if the user resizes the
% figure window. This is true for all text/edit boxes and buttons.

%% Create the axis with the local view. This is the main point of the
%  function, so it should dominate the figure.
%  All other object positions are based on this axis position.

axes1Xmin   = 0.05;
axes1Ymin   = 0.40;
axes1Width  = 0.65;
axes1Height = 0.57;
nudge       = 0.025;

handles.local_view = axes('Position', [axes1Xmin axes1Ymin axes1Width axes1Height]);
set(handles.local_view,'FontSize',atts.fontsize)

%% Creates text in the top right corner with the time elapsed in the model
%  Start at the top and work our way down.
% The problem with ui text object is there is no way to vertically align
% the text in the object.
% These are used for every 'button'

controlXmin   = axes1Xmin + axes1Width + nudge; % just off the right edge
controlWidth  = 1 - controlXmin - nudge;        % symmetric in whats left
controlHeight = axes1Height/6;   % button height for the two bigger buttons

Yposition   = 1 - controlHeight/5 - 2.0*nudge;
handles.ui_text_time = uicontrol('Style', 'Text', ...
    'Units', 'Normalized', ...
    'Position', [controlXmin Yposition controlWidth controlHeight/2] , ...
    'String', 'Time = 0' , ...
    'BackgroundColor', atts.background, ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.5);

%% Creates a button that does a single step ahead when clicked

Yposition = Yposition - controlHeight/2 - 2.0*nudge;
handles.ui_button_Single_Step = uicontrol('Style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Position', [controlXmin Yposition controlWidth controlHeight] , ...
    'String', 'Advance Model', ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.35, ...
    'Callback', @SingleStep);

%% Creates a button that continously does a series of single steps until
%  the Pause Auto Run button is hit

Yposition = Yposition - controlHeight - nudge;
handles.ui_button_Auto_Run = uicontrol('Style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Position', [controlXmin Yposition controlWidth controlHeight] , ...
    'String', 'Start Auto Run', ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.35, ...
    'Callback', @AutoRun);

%% Create a Radio Button with choices for the types of assimilation

buttonGroupHeight = controlHeight * 1.7;

Yposition = Yposition - nudge - buttonGroupHeight;
handles.ui_button_group_assimilation = uibuttongroup('Visible','off', ...
    'Position',[controlXmin Yposition controlWidth buttonGroupHeight], ...
    'BackgroundColor' , atts.background, ...
    'BorderType', 'none', ...
    'SelectionChangeFcn', @Assimilation_selection);

% Positions are relative to the button Group
radioXmin   = 0.05;
radioWidth  = 0.9;
radioHeight = 0.2;

radioYMin = 1 - radioHeight - nudge;
% Create four radio buttons in the button group.
handles.ui_radio_noAssimilation = uicontrol(handles.ui_button_group_assimilation, ...
    'Style','radiobutton', ...
    'String','No Assimilation', ...
    'Units' , 'Normalized', ...
    'Position',[radioXmin radioYMin radioWidth radioHeight], ...
    'BackgroundColor' , atts.background, ...
    'FontName', atts.fontname, ...
    'FontUnits', 'Normalized' , ...
    'FontSize' , 0.85, ...
    'HandleVisibility','off');

radioYMin = radioYMin - radioHeight - nudge;
handles.ui_radio_EAKF = uicontrol(handles.ui_button_group_assimilation, ...
    'Style','radiobutton', ...
    'String','EAKF', ...
    'Units' , 'Normalized', ...
    'Position',[radioXmin radioYMin radioWidth radioHeight], ...
    'BackgroundColor' , atts.background, ...
    'FontName', atts.fontname, ...
    'FontUnits', 'Normalized' , ...
    'FontSize' , 0.85, ...
    'HandleVisibility','off');

radioYMin = radioYMin - radioHeight - nudge;
handles.ui_radio_EnKF = uicontrol(handles.ui_button_group_assimilation, ...
    'Style','radiobutton', ...
    'String','EnKF', ...
    'Units' , 'Normalized', ...
    'Position',[radioXmin radioYMin radioWidth radioHeight], ...
    'BackgroundColor' , atts.background, ...
    'FontName', atts.fontname, ...
    'FontUnits', 'Normalized' , ...
    'FontSize' , 0.85, ...
    'HandleVisibility','off');

radioYMin = radioYMin - radioHeight - nudge;
handles.ui_radio_RHF = uicontrol(handles.ui_button_group_assimilation, ...
    'Style','radiobutton', ...
    'String','RHF', ...
    'Units' , 'Normalized', ...
    'Position',[radioXmin radioYMin radioWidth radioHeight], ...
    'BackgroundColor' , atts.background, ...
    'FontName', atts.fontname, ...
    'FontUnits', 'Normalized' , ...
    'FontSize' , 0.85, ...
    'HandleVisibility','off');

% Make the uibuttongroup visible after creating child objects.
set(handles.ui_button_group_assimilation, 'Visible','On');

%% Create a Reset button to clear the graphs and set the time back to 0
handles.ui_button_Reset = uicontrol('Style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Position', [0.744 0.429 0.095 0.063] , ...
    'String', 'Reset', ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontWeight', 'bold', ...
    'FontSize', 0.4, ...
    'Callback', @reset);

%% Create a 'local view' label for the graph in the top left
textViewWidth  = axes1Width/5;
textViewHeight = axes1Height/20;
textViewYmin   = axes1Ymin + axes1Height - textViewHeight;
handles.ui_text_Local_view = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [axes1Xmin textViewYmin textViewWidth textViewHeight] , ...
    'String', 'Local View' , ...
    'BackgroundColor', atts.background, ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.8);

%% Create a panel in the bottom left with a key as to what the
% different colors mean in the two graphs

keyPanelYmin   = nudge;
keyPanelWidth  = 0.65*axes1Width;
keyPanelHeight = axes1Ymin - keyPanelYmin - 2*nudge;

handles.Key_Panel = uipanel( ...
    'BackgroundColor' , atts.background, ...
    'Position',[0.029 0.028 0.423 0.325]);

% Since there are 5 lines of text to put in the keyPanel, getting the
% spacing is tricky for these ui text objects that don't align vertically.
% All these positions are relative to handles.Key_Panel

borderNudge   = 0.10;
legendXmin    = 0.05;
legendXwidth  = 1.0 - 2.0 * legendXmin;
legendYheight = 1.0/6;

% Again, start at the top and work our way down.

legendYmin    = 1.0 - legendYheight - borderNudge;
handles.ui_text_black_key = uicontrol(handles.Key_Panel, ...
    'Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [legendXmin legendYmin legendXwidth legendYheight] , ...
    'String', 'Black is True Trajectory' , ...
    'BackgroundColor' , atts.background, ...
    'HorizontalAlignment', 'center', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.5);

legendYmin = legendYmin - legendYheight;
handles.ui_text_green_key = uicontrol(handles.Key_Panel, ...
    'Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [legendXmin legendYmin legendXwidth legendYheight] ,  ...
    'String', 'Green is Ensemble' , ...
    'BackgroundColor' , atts.background, ...
    'HorizontalAlignment', 'center', ...
    'ForegroundColor', atts.green, ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.5);

legendYmin = legendYmin - legendYheight;
handles.ui_text_red_key = uicontrol(handles.Key_Panel, ...
    'Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [legendXmin legendYmin legendXwidth legendYheight] , ...
    'String', 'Red is Obs. Increments' , ...
    'BackgroundColor' , atts.background, ...
    'HorizontalAlignment', 'center', ...
    'ForegroundColor', atts.red, ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.5);

legendYmin = legendYmin - legendYheight;
handles.ui_text_red_star_key = uicontrol(handles.Key_Panel, ...
    'Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [legendXmin legendYmin legendXwidth legendYheight] , ...
    'String', 'Red * is Observation' , ...
    'BackgroundColor' , atts.background, ...
    'HorizontalAlignment', 'center', ...
    'ForegroundColor', atts.red, ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.5);

%>@ TODO FIXME implement a slider for obs error variance.
legendYmin = legendYmin - legendYheight;
handles.ui_text_obs_err_key = uicontrol(handles.Key_Panel, ...
    'Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [legendXmin legendYmin legendXwidth legendYheight] , ...
    'String', 'Obs. error SD is N(0,1)' , ...
    'BackgroundColor' , atts.background, ...
    'HorizontalAlignment', 'center', ...
    'ForegroundColor', 'Black', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.5);

%% Creates one graph in the bottom right
axes2Xmin = axes1Xmin + keyPanelWidth + 2*nudge; %0.05 away from panel
axes2Ymin = keyPanelYmin;
axes2Width = 1 - axes2Xmin - nudge; %0.03 away from edge, looks balanced because of axis label taking up space
axes2Height = keyPanelHeight + 2*nudge;

handles.global_view = axes('Position', [axes2Xmin axes2Ymin axes2Width axes2Height]);
set(handles.global_view,'FontSize',atts.fontsize)

%% Create a text which is a label for global_view

handles.ui_text_Global_View = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.844 0.394 0.13 0.028] , ...
    'String', 'Global View' , ...
    'HorizontalAlignment', 'right' , ...
    'BackgroundColor', atts.background, ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.8);

% Initialize the l63 model
L63 = lorenz_63_static_init_model;

DELTAT     = L63.deltat;
SIGMA      = L63.sigma;
R          = L63.r;
B          = L63.b;
MODEL_SIZE = L63.model_size;

% Call the reset button to initialize all the handles variables
reset;

%% ----------------------------------------------------------------------
%  All function below can use the variables defined above.
%% ----------------------------------------------------------------------

    function SingleStep(~,~)
        %This Function is called whenever the button_Single_Step is
        %pressed. It disables all the buttons, calls step ahead, then re-enables
        %all the buttons
        
        %Disable all the buttons
        set(handles.ui_button_Single_Step,          'Enable', 'Off');
        set(handles.ui_radio_noAssimilation,        'Enable', 'Off');
        set(handles.ui_radio_EAKF,                  'Enable', 'Off');
        set(handles.ui_radio_EnKF,                  'Enable', 'Off');
        set(handles.ui_radio_RHF,                   'Enable', 'Off');
        set(handles.ui_button_Auto_Run,             'Enable', 'Off');
        set(handles.ui_button_Reset,                'Enable', 'Off');
        
        %Advance the model one time.
        step_ahead;
        
        %If the menu has a value of 'No Assimilation' then go ahead and
        %call step ahead one more time as there are no observations to
        %assimilate
        if (strcmp(handles.filter_kind, 'No Assimilation') && ...
                strcmp(get(handles.ui_button_Single_Step, 'String'), 'Assimilate Obs'))
            step_ahead();
        end
        
        %Re-Enable All the buttons
        set(handles.ui_button_Auto_Run,         'Enable', 'On');
        set(handles.ui_button_Single_Step,      'Enable', 'On');
        set(handles.ui_radio_noAssimilation,    'Enable', 'On');
        set(handles.ui_radio_EAKF,              'Enable', 'On');
        set(handles.ui_radio_EnKF,              'Enable', 'On');
        set(handles.ui_radio_RHF,               'Enable', 'On');
        set(handles.ui_button_Reset,            'Enable', 'On');
        
    end

%% ----------------------------------------------------------------------

    function AutoRun(~,~)
        
        % This function is called when the ui_button_Auto_Run is pressed. It
        % continuously calls the step_ahead function until the ui_button_Auto_Run is
        % pressed again.
        
        % Turn off all the other model status controls to avoid a mess
        % MAKE SURE TO INCLUDE OTHER CONTROLS HERE
        set(handles.ui_button_Single_Step,          'Enable', 'Off');
        % In 2015, there is a way to disable the button group, but it is not
        % compatible with 2014, so we must enable/disable each radio button
        % seperately
        set(handles.ui_radio_noAssimilation,        'Enable', 'Off');
        set(handles.ui_radio_EAKF,                  'Enable', 'Off');
        set(handles.ui_radio_EnKF,                  'Enable', 'Off');
        set(handles.ui_radio_RHF,                   'Enable', 'Off');
        set(handles.ui_button_Reset,                'Enable', 'Off');
        
        % Check the label to see if we are starting or stopping a free run
        if(strcmp(get(handles.ui_button_Auto_Run, 'String'), 'Pause Auto Run'))
            
            % Being told to stop; switch to not running status
            set(handles.ui_button_Auto_Run, 'Enable', 'Off');
            set(handles.ui_button_Auto_Run, 'String', 'Start Auto Run');
            
        else
            % Being told to start free run
            
            set(handles.ui_button_Auto_Run, 'String', 'Pause Auto Run');
            
            % Loop through advance and assimilate steps until stopped
            while(true)
                
                % Check to see if stop has been pushed
                status_string = get(handles.ui_button_Auto_Run, 'String');
                
                if(strcmp(status_string, 'Start Auto Run'))
                    
                    % Turn all the other model status controls back on
                    % MAKE SURE TO INCLUDE OTHER CONTROLS HERE
                    set(handles.ui_button_Single_Step,     'Enable', 'On');
                    set(handles.ui_radio_noAssimilation,   'Enable', 'On');
                    set(handles.ui_radio_EAKF,             'Enable', 'On');
                    set(handles.ui_radio_EnKF,             'Enable', 'On');
                    set(handles.ui_radio_RHF,              'Enable', 'On');
                    set(handles.ui_button_Reset,           'Enable', 'On');
                    % Very last, turn on the start free run button
                    set(handles.ui_button_Auto_Run,        'Enable', 'On');
                    
                    return
                end
                % Do the next advance or assimilation step
                step_ahead;
                pause(0.1)
            end
        end
    end

%% ----------------------------------------------------------------------

    function step_ahead()
        % Moves the model ahead or assimilates next observations
        
        % Test on semaphore, either advance or assimilate
        if(handles.ready_to_advance)
            % Set semaphore to indicate that next step is an assimilation
            handles.ready_to_advance = false;
            
            % Set the text to indicate that next step is an assimilate
            set(handles.ui_button_Single_Step, 'String', 'Assimilate Obs');
            
            % Turn off the recent observation plot if it exists
            if(handles.global_init)
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
                for imem = 1:handles.ens_size
                    [new_ens, new_time] = lorenz_63_adv_1step(handles.post(time, :, imem), time);
                    handles.prior(new_time, :, imem) = new_ens;
                end
                
                % Plot a long trajectory of truth in small window for reference
                axes(handles.global_view); %#ok<LAXES>
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
                axis([-25 25 -25 25 5 45]);
                
                % Plot the close-up view of the ensemble
                axes(handles.local_view) %#ok<LAXES>
                
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
                for imem = 1:handles.ens_size
                    %Axes continually changes, so reset grid to on
                    grid on;
                    plot3(handles.prior(btime:new_time, 1, imem), ...
                        handles.prior(btime:new_time, 2, imem), ...
                        handles.prior(btime:new_time, 3, imem), '-', 'Color', atts.green)
                end
                
                % Update the time label
                set(handles.ui_text_time, 'String', ['Time = ', num2str(new_time)]);
                
                % Force the buffers to flush and plot the advance
                drawnow
                
                % Last prior update will get overwritten when assimilation is done
                handles.post(new_time, :, :) = handles.prior(new_time, :, :);
            end
            
            % Compute the observations for this time and save
            for i = 1:3
                handles.obs(i) = handles.true_state(new_time, i) + ...
                    handles.obs_sd * randn;
            end
            
            % Plot the observation as a red asterisk in both axes
            h = plot3(handles.obs(1), handles.obs(2), handles.obs(3), ...
                'r*', 'MarkerSize', 20);
            set(h,'Color',atts.red);
            
            axes(handles.global_view);
            handles.h_global_obs = ...
                plot3(handles.obs(1), handles.obs(2), handles.obs(3), ...
                'r*', 'MarkerSize', 20);
            set(handles.h_global_obs,'Color',atts.red);
            handles.global_init = true;
            axes(handles.local_view);
            
        else
            % Set semaphore to indicate that next step is a model advance
            handles.ready_to_advance = true;
            
            % Set the pushbutton text to indicate that the next step is a model advance
            set(handles.ui_button_Single_Step, 'String', 'Advance Model');
            
            % Get current time step
            time = handles.time;
            
            % Determine what type of assimilation is being done (none, EAKF, EnKF, RHF)
            if(strcmp(handles.filter_kind, 'No Assimilation'))
                % Just copy prior to posterior
                handles.post(time, :, :) = handles.prior(time, :, :);
            else
                % Code for doing the assimilation comes here
                
                % Do fully sequential assimilation algorithm
                temp_ens = squeeze(handles.prior(time, :, :));
                
                % Observe each state variable independently
                obs = zeros(1,3);
                for i = 1:3
                    obs_prior = temp_ens(i, :);
                    obs(i) = handles.obs(i);
                    
                    % Compute the increments for observed variable
                    switch handles.filter_kind
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
                for imem = 1:handles.ens_size
                    xup = [handles.prior(time, 1, imem), handles.post(time, 1, imem)];
                    yup = [handles.prior(time, 2, imem), handles.post(time, 2, imem)];
                    zup = [handles.prior(time, 3, imem), handles.post(time, 3, imem)];
                    h = plot3(xup, yup, zup, 'r');
                    set(h,'Color',atts.red)
                end
            end
        end
        
    end

%% ----------------------------------------------------------------------

    function reset(~,~)
        %This function resets handles to it's original values and clears the graphs
        % Also, called at the beginning to initialize the variables
        
        % set random number seed to same value to generate known sequences
        % rng('default')  is the Mersenne Twister with seed 0
        rng(0,'twister')
        
        % Global semaphore; ready to advance or assimilate?
        handles.ready_to_advance = true;
        
        %Set handles to original values
        ens_size                   = 20;
        handles.ens_size           = ens_size;
        handles.model_size         = MODEL_SIZE;
        handles.true_state(1, 1:3) = [10.2471, 2.8443, 36.2666];
        handles.time               = 1;
        handles.hstar              = 0;
        handles.obs(1:3)           = 0;
        handles.obs_sd             = 1;
        handles.obs_error_var      = handles.obs_sd^2;
        handles.h_global_obs       = [];
        handles.global_init        = false;
        handles.filter_kind        = 'No Assimilation';
        
        handles.post = zeros(1, MODEL_SIZE, ens_size);
        for n = 1:handles.ens_size
            handles.post(1, :, n) = handles.true_state(1, :) + ...
                0.1 * randn(1, MODEL_SIZE);
        end
        
        % Make first prior identical to the first posterior
        handles.prior = handles.post;
        
        %clears the two graphs
        fontsize = get(handles.local_view,'FontSize');
        cla(handles.local_view); set(handles.local_view,'FontSize',fontsize)
        cla(handles.global_view); set(handles.global_view,'FontSize',fontsize)
        
        set(handles.ui_text_time         , 'String', 'Time = 0');
        set(handles.ui_button_Single_Step, 'String', 'Advance Model');
        set(handles.ui_button_group_assimilation,'SelectedObject',handles.ui_radio_noAssimilation);
        
        % Plot the initial state
        initial_plot(handles.local_view);
        initial_plot(handles.global_view);
        
    end

%% ----------------------------------------------------------------------

    function initial_plot(hax)
        % Function is called in reset function, plots the initial values for
        % the 2 axes
        axes(hax)
        hold off;
        FontSize = get(hax,'FontSize');
        
        % Plot an asterisk
        handles.initial_ob = plot3(handles.true_state(1, 1), ...
            handles.true_state(1, 2), ...
            handles.true_state(1, 3), 'k*', 'MarkerSize', 16, 'LineWidth', 2);
        
        view([2, -1 1]); % Set an appropriate consistent view angle
        
        % Adjust the axes to follow the truth
        xb = handles.true_state(1, 1);
        yb = handles.true_state(1, 2);
        zb = handles.true_state(1, 3);
        limits = [xb - 3,  xb + 3,  yb - 3,  yb + 3, zb - 3, zb + 3];
        axis(limits);
        hold on
        
        % Plot the ensemble members
        for imem = 1:handles.ens_size
            plot3(handles.prior(1, 1, imem), ...
                handles.prior(1, 2, imem), ...
                handles.prior(1, 3, imem), 'k.', 'Color', atts.green);
        end
        grid on;
        set(gca, 'FontSize', FontSize)
    end

%% ----------------------------------------------------------------------

    function Assimilation_selection(~, eventdata)
        % Function is called whenever a radio button has been selected, it sets
        % the global filter variable
        
        % eventdata refers to the data in the GUI when a radio button in the
        % group is changed
        
        % Set the filter_type_string to newest radiobutton Value
        handles.filter_kind = get(eventdata.NewValue,'String');
    end

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
