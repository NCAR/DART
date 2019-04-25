function twod_ensemble
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
% See also: gaussian_product.m oned_model.m oned_model_inf.m oned_ensemble.m
%           run_lorenz_63.m run_lorenz_96.m run_lorenz_96_inf.m

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

help twod_ensemble;

atts = stylesheet;  % get the default fonts and color

figureXmin   = 100;
figureYmin   = 250;
figureWidth  = 900;
figureHeight = 480;

%% Create a figure layout
handles.figure1 = figure('Name', 'twod_ensemble', ...
    'Position', [figureXmin figureYmin figureWidth figureHeight], ...
    'Color', atts.background);

% Position has units of normalized and therefore the components must be
% Fractions of figure. Also, the actual text  has FontUnits which are
% normalized, so the Font Size is a fraction of the text box/ edit box that
% contains it. Making the font size and box size normalized allows the text
% size and box size to change proportionately if the user resizes the
% figure window. This is true for all text/edit boxes and buttons.

%% Create a button that lets you create a new ensemble.
handles.ui_button_create_ensemble = uicontrol('Style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Position', [480/figureWidth 420/figureHeight 175/figureWidth 40/figureHeight] , ...
    'String', 'Create New Ensemble', ...
    'BackgroundColor', 'White', ...
    'Callback', @create_ensemble_Callback, ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.4);

%% Create a Button that updates the ensemble.
% This button is disabled at first because there is no ensemble to update.
handles.ui_button_update_ensemble = uicontrol('Style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Position', [480/figureWidth 360/figureHeight 175/figureWidth 40/figureHeight] , ...
    'String', 'Update Ensemble', ...
    'BackgroundColor', 'White', ...
    'Callback', @update_ensemble_Callback, ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.4, ...
    'Enable', 'Off');

%% Create Radio Button group with choices for the types of assimilation

handles.ui_assimilation_buttons = uibuttongroup('Visible','off', ...
    'Position', [480/figureWidth 260/figureHeight 150/figureWidth 80/figureHeight] , ...
    'BorderType', 'none', ...
    'BackgroundColor', atts.background, ...
    'SelectionChangeFcn', @Assimilation_selection);

% Positions are relative to the Button Group
radioXmin   = 0.05;
radioWidth  = 0.9;
nudge       = 0.025;
radioHeight = (1 - 4*nudge)/3;

radioYmin   = 1.0 - radioHeight - nudge;
handles.ui_radio_EAKF = uicontrol(handles.ui_assimilation_buttons, ...
    'Style','radiobutton', ...
    'Units' , 'Normalized', ...
    'Position',[radioXmin radioYmin radioWidth radioHeight], ...
    'String','EAKF', ...
    'BackgroundColor' , atts.background, ...
    'Foreground' , 'Black', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'Normalized' , ...
    'FontWeight', 'Bold' , ...
    'FontSize', 0.8, ...
    'HandleVisibility','off');

radioYmin   = radioYmin - radioHeight - nudge;
handles.ui_radio_EnKF = uicontrol(handles.ui_assimilation_buttons, ...
    'Style','radiobutton', ...
    'String','EnKF', ...
    'Units' , 'Normalized', ...
    'Position',[radioXmin radioYmin radioWidth radioHeight], ...
    'BackgroundColor' , atts.background, ...
    'FontName', atts.fontname, ...
    'FontUnits', 'Normalized' , ...
    'FontWeight', 'Bold' , ...
    'FontSize' , 0.8, ...
    'HandleVisibility','off');

radioYmin   = radioYmin - radioHeight - nudge;
handles.ui_radio_RHF = uicontrol(handles.ui_assimilation_buttons, ...
    'Style','radiobutton', ...
    'String','RHF', ...
    'Units' , 'Normalized', ...
    'Position',[radioXmin radioYmin radioWidth radioHeight], ...
    'BackgroundColor' , atts.background, ...
    'FontName', atts.fontname, ...
    'FontUnits', 'Normalized' , ...
    'FontWeight', 'Bold' , ...
    'FontSize' , 0.8, ...
    'HandleVisibility','off');

% Make the uibuttongroup visible after creating child objects.
set(handles.ui_assimilation_buttons, 'Visible','On');
selected = get(handles.ui_assimilation_buttons,'SelectedObject');
handles.filter_type = get(selected,'String');

%% Create a red panel with two text boxes and two edit boxes next to them

handles.ObservationPanel = uipanel('BackgroundColor', atts.red, ...
    'BorderType','none', ...
    'Units', 'Normalized', ...
    'Position',[670/figureWidth 370/figureHeight 210/figureWidth 90/figureHeight]);

% Positions are relative to the panel
nudge  = 0.05;
Xmin   = 0.05;
dX1    = 0.65;
dX2    = 1.0 - dX1 - 2*nudge;
Height = (1.0 - 3*nudge)/2;

Ymin   = 1.0 - Height - nudge;
handles.ui_text_observation = uicontrol(handles.ObservationPanel, ...
    'Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [Xmin Ymin-nudge dX1 Height] , ...
    'String', 'Observation' , ...
    'BackgroundColor', atts.red, ...
    'ForegroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontWeight', 'Bold', ...
    'FontSize', 0.50);

handles.ui_edit_observation = uicontrol(handles.ObservationPanel, ...
    'Style', 'edit', ...
    'Units', 'Normalized', ...
    'Position', [Xmin+dX1 Ymin dX2 Height] , ...
    'String', '5' , ...
    'Callback', @edit_observation_Callback, ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.50);

Ymin   = Ymin - Height - nudge;
handles.ui_text_obs_error_sd = uicontrol(handles.ObservationPanel, ...
    'Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [Xmin Ymin-nudge dX1 Height] , ...
    'String', 'Obs. Error SD' , ...
    'BackgroundColor', atts.red, ...
    'ForegroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontWeight', 'Bold', ...
    'FontSize', 0.50);

handles.ui_edit_obs_error_sd = uicontrol(handles.ObservationPanel, ...
    'Style', 'edit', ...
    'Units', 'Normalized', ....
    'Position', [Xmin+dX1 Ymin dX2 Height] , ...
    'String', '1' , ...
    'Callback', @edit_obs_error_sd_Callback, ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.50);

%% Creates an error text inside the graph that is invisible at first but
% appears if a user enters an incorrect number for the observation or the
% observation error sd
handles.ui_text_error = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.220 0.550 0.20 0.20] , ...
    'String', 'Error' , ...
    'BackgroundColor', 'white', ...
    'ForegroundColor', atts.red, ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontWeight', 'Bold', ...
    'FontSize', 0.2, ...
    'Visible', 'Off');

%Set up all the handles variables and create the graphs
initialize();

%% Creates the graph of the observation and the observation error sd
%  well as the ensemble estimates of the observation.

h_observation  = get(handles.ui_edit_observation);
h_obs_error_sd = get(handles.ui_edit_obs_error_sd);
observation    = str2double(h_observation.String);
obs_error_sd   = str2double(h_obs_error_sd.String);

plotposition = [ 40/figureWidth 110/figureHeight  60/figureWidth 350/figureHeight; ...
    115/figureWidth 110/figureHeight 340/figureWidth 350/figureHeight; ...
    115/figureWidth  40/figureHeight 340/figureWidth  55/figureHeight];

%% Create the axis for the unobserved marginal graphic

handles.h_unobMarginal = axes('Position', plotposition(1,:), ...
    'FontName', atts.fontname, ...
    'FontSize', atts.fontsize);
axis([0 1 0 10]);
ylabel('Unobserved State Variable', 'Fontsize', atts.fontsize);
hold on

%% Get a subplot for the joint

handles.h_joint = axes('Position', plotposition(2,:), ...
    'FontName', atts.fontname, ...
    'FontSize', atts.fontsize);
axis([0 10 0 10]);
grid on
hold on
title('Joint Distribution');

h_click    = text(5, 9, {'Click inside Joint Distribution plot','to create member.'}, ...
    'FontSize', 16, 'HorizontalAlignment','center');
h_finish   = text(5, 8, 'Click outside of plot to finish.', ...
    'FontSize', 16, 'HorizontalAlignment','center');
h_err_text = text(5, 4, 'An ensemble must have at least 2 members.', 'Color', atts.red, ...
    'FontSize', 16, 'HorizontalAlignment','center');

handles.h_correl = text(5, 9, ' ', 'Color', atts.green, 'FontWeight', 'Bold', ...
    'FontSize', 16, 'HorizontalAlignment','center');

%% Create a subplot for the observed variable marginal

handles.h_obMarginal = axes('Position', plotposition(3,:), ...
    'FontName', atts.fontname, ...
    'FontSize', atts.fontsize);

% May want to mess with vertical axis for prior density
axis([0 10 0 1]);
hold on
xlabel('Observed Quantity');
handles.h_obs_marg = plot(observation, 0, '*', 'MarkerSize', 16, 'LineWidth', 2.0);
set(handles.h_obs_marg,'Color',atts.red)

% Turn off the unwanted tick labels
set(handles.h_unobMarginal, 'Xticklabel', []);
set(handles.h_joint       , 'Xticklabel', [], 'Yticklabel', []);
set(handles.h_obMarginal  , 'Yticklabel', []);

%This graph is the graph of the observation
handles.h_obs_likelihood = axes( ...
    'Position', [500/figureWidth 40/figureHeight 390/figureWidth 200/figureHeight], ...
    'FontName', atts.fontname, ...
    'FontSize', atts.fontsize);

handles.h_marg_obs_plot = plot_gaussian(observation, obs_error_sd, 1);

axis([0 10 -handles.y_max/5 1]);

% Set the ticks
set(gca, 'YTick',      [0 0.2 0.4 0.6 0.8]);
set(handles.h_marg_obs_plot, 'Color', atts.red, ...
    'Linestyle', '--', ...
    'LineWidth', 2);
xlabel('Observed Quantity',      'FontSize', atts.fontsize);
ylabel('Observation Likelihood', 'FontSize', atts.fontsize);
title('Marginal Distribution of Observation', 'FontSize', atts.fontsize);
hold on

% Plot an asterisk for the observed value
handles.h_obs_ast = plot(observation, 0, 'r*', 'MarkerSize', 16, 'LineWidth', 2.0);
set(handles.h_obs_ast,'Color',atts.red)

% Plot an axis; display is fixed from x = 0 to 10
plot([0 10], [0 0], 'k', 'LineWidth', 2);

%% ----------------------------------------------------------------------------

    function create_ensemble_Callback(~,~)
        % Allows the user to create a new ensemble in the left axes
        
        % Disable the update ensemble button and all other active buttons
        set(handles.ui_button_create_ensemble, 'Enable', 'Off');
        set(handles.ui_button_update_ensemble, 'Enable', 'Off');
        set(handles.ui_edit_observation,       'Enable', 'Off');
        set(handles.ui_edit_obs_error_sd,      'Enable', 'Off');
        
        % Clear out any old ensemble members if they exist
        axes(handles.h_joint)
        for i = 1:handles.ens_size
            set(handles.h_ens_member(i),    'Visible', 'off');
            set(handles.h_gui_marg(i),      'Visible', 'off');
            set(handles.h_unobs(i),         'Visible', 'off');
            set(handles.h_marg(i),          'Visible', 'off');
        end
        
        % Turn off any posterior old plotting
        set(handles.h_update_ens,   'Visible', 'off');
        set(handles.h_marg_update,  'Visible', 'off');
        set(handles.h_marg_inc,     'Visible', 'off');
        set(handles.h_marg_state,   'Visible', 'off');
        set(handles.h_state_inc,    'Visible', 'off');
        set(handles.h_joint_update, 'Visible', 'off');
        set(handles.h_joint_inc,    'Visible', 'off');
        
        % Clear out the old best fit line
        set(handles.h_best_fit, 'Visible', 'off');
        set(handles.h_correl,   'Visible', 'off');
        
        % Work in the joint distribution plot
        axes(handles.h_joint);
        hold on
        
        % Need to guarantee at least 2 ensemble members
        ens_size = 0;
        
        while ens_size < 1000
            
            [xt, yt] = ginput(1);
            gca;
            % Make sure that the click was in the correct set of axes
            % Terminate by clicking outside of graph range
            if(xt < 0 || xt > 10 || yt < 0 || yt > 10 || gca ~= handles.h_joint)
                axes(handles.h_joint); %#ok<LAXES>
                break;
            else
                
                ens_size = ens_size + 1;
                x(1, ens_size) = xt; %#ok<AGROW>
                x(2, ens_size) = yt; %#ok<AGROW>
                
                axes(handles.h_joint); %#ok<LAXES>
                handles.h_ens_member(ens_size) = ...
                    plot(x(1, ens_size), x(2, ens_size), '*', ...
                    'MarkerSize', 16, 'Color', atts.green, 'LineWidth',2.0);
                
                % Plot the marginal for the unobserved state variable
                %>@ TODO  POSSIBLE IMPROVEMENT ... annotate new marginal mean, sd
                axes(handles.h_unobMarginal); %#ok<LAXES>
                handles.h_unobs(ens_size) = ...
                    plot(0, x(2, ens_size), '*', 'MarkerSize', 16, 'Color', atts.green, 'LineWidth',2.0);
                
                % Plot the marginal for the observed quantity
                axes(handles.h_obMarginal); %#ok<LAXES>
                handles.h_marg(ens_size) = ...
                    plot(x(1, ens_size), 0, '*', 'MarkerSize', 16, 'Color', atts.green, 'LineWidth',2.0);
                
                % Plot the marginal in the gui frame
                axes(handles.h_obs_likelihood); %#ok<LAXES>
                handles.h_gui_marg(ens_size) = ...
                    plot(x(1, ens_size), 0, '*', 'MarkerSize', 16, 'Color', atts.green, 'LineWidth',2.0);
                
                % Then switch back to axes(handles.h_joint)
                axes(handles.h_joint); %#ok<LAXES>
                
                if (ens_size < 2)
                    continue
                end
                
                % Clear out the error message if it's been made visible
                set(h_err_text, 'Visible', 'off');
                set(h_click,    'Visible', 'off');
                
                prior_correl = corrcoef(x(1, :), x(2, :));
                str1         = sprintf('Correlation = %f', prior_correl(1,2));
                set(handles.h_correl,'String', str1, 'Visible', 'on')
                
            end
        end
        
        % it is possible that they click outside the box before completing a viable
        % ensemble ... in this case, just return and let them start over.
        if (ens_size > 0)
            
            % Turn off the data entry messages
            set(h_finish, 'Visible', 'off');
            
            %% Ensemble created, compute mean and sd, clean up and return
            % Set the global gui storage
            handles.ens_size    = ens_size;
            handles.ens_members = x;
            
            % Plot the best fit line on the ensemble
            prior_mean = mean(x, 2);
            prior_cov  = cov(x(1, :), x(2, :));
            slope      = prior_cov(1, 2) / var(x(1, :));
            
            best_x = [0 10];
            best_y(1) = prior_mean(2) - (prior_mean(1)) * slope;
            best_y(2) = best_y(1) + 10 * slope;
            handles.h_best_fit = plot(best_x, best_y, 'g', 'LineWidth', 2.0);
            set(handles.h_best_fit, 'Color', atts.green);
            
        end
        
        % Enable the update ensemble button
        set(handles.ui_button_create_ensemble, 'Enable', 'On');
        set(handles.ui_button_update_ensemble, 'Enable', 'On');
        set(handles.ui_edit_observation,       'Enable', 'On');
        set(handles.ui_edit_obs_error_sd,      'Enable', 'On');
        
        % Reset focus to the menu gui window
        axes(handles.h_obs_likelihood);
        
    end

%% -------------------------------------------------------------------------

    function update_ensemble_Callback(~,~)
        % Uses the assimilation to update the ensemble according to the
        % observation, and then plots it on the main graph on the left, the
        % two marginals, and the right observation graph
        
        axes(handles.h_obs_likelihood);
        % Turn off any old points
        set(handles.h_update_ens,   'Visible', 'off');
        set(handles.h_marg_update,  'Visible', 'off');
        set(handles.h_marg_inc,     'Visible', 'off');
        set(handles.h_marg_state,   'Visible', 'off');
        set(handles.h_state_inc,    'Visible', 'off');
        set(handles.h_joint_update, 'Visible', 'off');
        set(handles.h_joint_inc,    'Visible', 'off');
        
        ensemble = handles.ens_members;
        h_observation  = get(handles.ui_edit_observation);
        h_obs_error_sd = get(handles.ui_edit_obs_error_sd);
        observation    = str2double(h_observation.String);
        obs_error_sd   = str2double(h_obs_error_sd.String);
        
        %If ensemble is not empty
        if (size(ensemble,2) > 0)
            switch handles.filter_type
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
            
            %Set the y-coordinate of the ensembles, to be halfway between 0 and
            %the bottom of the graph;
            y(1:handles.ens_size) = -handles.y_max/10;
            
            handles.h_update_ens = plot(new_ensemble, y, '*', 'MarkerSize', 16, 'Color', atts.blue);
            
            % Plot the increments in the state marginal plot
            axes(handles.h_obMarginal);
            
            % Need to sort ensemble to get nice ordering for increments
            [~, sort_obs_ind] = sort(ensemble(1, :));
            for i = 1:handles.ens_size
                y(i) = i / (handles.ens_size + 1);
                handles.h_marg_update(i) = ...
                    plot(new_ensemble(sort_obs_ind(i)), y(i), '*', 'MarkerSize', 16, 'Color', atts.blue);
                % Also plot a segment in blue
                handles.h_marg_inc(i) = ...
                    plot([ensemble(1, sort_obs_ind(i)), new_ensemble(1, sort_obs_ind(i))], ...
                    [y(i), y(i)], 'c');
            end
            
            % Figure out the increments for the unobserved variable
            
            axes(handles.h_unobMarginal);
            
            covar     = cov(ensemble');
            state_inc = obs_increments * covar(1, 2) / covar(1, 1);
            new_state = ensemble(2, :) + state_inc;
            %>@ TODO POSSIBLE IMPROVEMENT ... annotate new marginal mean, sd
            
            % Now need to sort the state variable ensemble to get nice ordering
            [~, sort_ind] = sort(ensemble(2, :));
            for i = 1:handles.ens_size
                handles.h_marg_state(i) = ...
                    plot(y(i), new_state(sort_ind(i)), '*', 'MarkerSize', 16, 'Color', atts.blue);
                % Also plot a segment in blue
                handles.h_state_inc(i) = plot([y(i), y(i)], ...
                    [ensemble(2, sort_ind(i)), new_state(sort_ind(i))], 'c');
            end
            
            % Plot the updated joint distribution points
            axes(handles.h_joint);
            for i = 1:handles.ens_size
                handles.h_joint_update(i) = plot(new_ensemble(i), new_state(i), ...
                    '*', 'MarkerSize', 16, 'Color', atts.blue);
                handles.h_joint_inc(i) = plot([ensemble(1, i), new_ensemble(1, i)], ...
                    [ensemble(2, i), new_state(i)], 'c');
            end
            
            % Return the focus to the window with pushbuttons
            axes(handles.h_obs_likelihood);
            
        end
    end

%% -------------------------------------------------------------------------

    function edit_observation_Callback(~, ~)
        
        % Enable things that an error might have turned off
        set(handles.ui_edit_obs_error_sd,      'Enable', 'on');
        set(handles.ui_button_create_ensemble, 'Enable', 'on');
        
        % Only enable the update ensemble pushbutton if an ensemble has been created
        if(handles.ens_size > 0)
            set(handles.ui_button_update_ensemble, 'Enable', 'on');
        end
        
        % Get the value of the observation
        if( isfinite( str2double(    get(handles.ui_edit_observation, 'String'))))
            observation = str2double(get(handles.ui_edit_observation, 'String'));
            
            if (observation > 10)
                set(handles.ui_edit_observation, 'String', '<10!');
                input_error('observation');
                return;
            elseif (observation < 0)
                set(handles.ui_edit_observation, 'String', '>0!');
                input_error('observation');
                return;
            end
            
            %Set background color to normal and error text off
            set(handles.ui_edit_observation, 'BackgroundColor', 'white');
            set(handles.ui_text_error,  'Visible', 'Off');
            set(handles.h_ens_member,   'Visible', 'On');
            set(handles.h_best_fit,     'Visible', 'On');
            set(handles.h_joint_update, 'Visible', 'On');
            set(handles.h_joint_inc,    'Visible', 'On');
        else
            set(handles.ui_edit_observation, 'String', '?');
            input_error('observation');
            return
        end
        
        % Get the value of the observation error sd
        h_obs_error_sd = get(handles.ui_edit_obs_error_sd);
        obs_error_sd   = str2double(h_obs_error_sd.String);
        
        % Plot the updated distribution
        set(handles.h_marg_obs_plot, 'Visible', 'off');
        handles.h_marg_obs_plot = plot_gaussian(observation, obs_error_sd, 1);
        set(handles.h_marg_obs_plot, 'Color', atts.red, 'Linestyle', '--', 'Linewidth', 2);
        
        % Update the observation asterisk
        set(handles.h_obs_ast, 'Visible', 'off');
        handles.h_obs_ast = plot(observation, 0, 'r*', 'MarkerSize', 16,'LineWidth',2.0);
        set(handles.h_obs_ast,'Color',atts.red)
        
        % Plot the updated obs distribution on the marginal subplot
        axes(handles.h_obMarginal);
        
        % Plot the updated observation in the marginal
        set(handles.h_obs_marg, 'Visible', 'off');
        handles.h_obs_marg = plot(observation, 0, 'r*', 'MarkerSize', 16,'LineWidth',2.0);
        set(handles.h_obs_marg,'Color',atts.red)
        
        % Replot the update ensemble members so the correlate to new
        % observation
        update_ensemble_Callback();
        
        axes(handles.h_obs_likelihood);
        
    end

%% --------------------------------------------------------------------

    function edit_obs_error_sd_Callback(~, ~)
        
        % Enable things that an error might have turned off
        set(handles.ui_edit_observation,       'Enable', 'on')
        set(handles.ui_button_create_ensemble, 'Enable', 'on')
        
        % Only enable the update ensemble pushbutton if an ensemble has been created
        if(handles.ens_size > 0)
            set(handles.ui_button_update_ensemble, 'Enable', 'on');
        end
        
        % Get the value of the observation error standard deviation
        if(isfinite(str2double(get(handles.ui_edit_obs_error_sd, 'String'))) && ...
                str2double(get(handles.ui_edit_obs_error_sd, 'String')) > 0)
            obs_error_sd = str2double(get(handles.ui_edit_obs_error_sd, 'String'));
            
            %Set background color to normal and error text off
            set(handles.ui_edit_obs_error_sd, 'BackgroundColor', 'white');
            set(handles.ui_text_error,  'Visible', 'Off');
            set(handles.h_ens_member,   'Visible', 'On');
            set(handles.h_best_fit,     'Visible', 'On');
            set(handles.h_joint_update, 'Visible', 'On');
            set(handles.h_joint_inc,    'Visible', 'On');
            
        else
            set(handles.ui_edit_obs_error_sd, 'String', '?');
            input_error('standard deviation');
            return
        end
        
        % Get the value of the observation
        h_observation = get(handles.ui_edit_observation);
        observation   = str2double(h_observation.String);
        handles.y_max = norm_pdf(observation, observation, obs_error_sd);
        
        %Give 0.2 cushion to y_max
        handles.y_max = handles.y_max + 0.2;
        
        %Update the axis based on the new y_max
        axis([0 10 -handles.y_max/5 handles.y_max]);
        
        set(gca,'YTickMode','auto')
        ticks    = get(gca,'YTick');
        inds     = (ticks >= 0); % Only show ticks for values greater than 0
        newticks = ticks(inds);
        set(gca,'YTick',newticks)
        
        % Replot the update ensemble members so the correlate to new obs_sd
        update_ensemble_Callback();
        
        % Plot the updated distribution on the menu plot
        
        set(handles.h_marg_obs_plot, 'Visible', 'off');
        handles.h_marg_obs_plot = plot_gaussian(observation, obs_error_sd, 1);
        set(handles.h_marg_obs_plot, 'Color', atts.red, 'Linestyle', '--', 'Linewidth', 2);
        
        % Update the observation asterisk
        
        set(handles.h_obs_ast, 'Visible', 'off');
        handles.h_obs_ast = plot(observation, 0, 'r*', 'MarkerSize', 16,'LineWidth',2.0);
        set(handles.h_obs_ast, 'Color', atts.red)
        
        % Plot the updated observation in the marginal
        axes(handles.h_obMarginal);
        
        set(handles.h_obs_marg, 'Visible', 'off');
        handles.h_obs_marg = plot(observation, 0, 'r*', 'MarkerSize', 16,'LineWidth',2.0);
        set(handles.h_obs_marg, 'Color', atts.red)
        
        % Reset focus to the menu gui window
        axes(handles.h_obs_likelihood);
        
    end

%% ----------------------------------------------------------------------------
    function initialize()
        % Insert the ensemble structure into this
        handles.ens_size        = 0;
        handles.ens_members     = [];
        handles.h_update_ens    = [];
        handles.h_ens_member    = [];
        handles.h_best_fit      = [];
        handles.h_marg_obs_plot = [];
        handles.h_obs_ast       = [];
        handles.h_obs_marg      = [];
        handles.h_gui_marg      = [];
        handles.h_unobs         = [];
        handles.h_marg          = [];
        handles.h_marg_update   = [];
        handles.h_marg_inc      = [];
        handles.h_marg_state    = [];
        handles.h_state_inc     = [];
        handles.h_joint_update  = [];
        handles.h_joint_inc     = [];
        handles.h_correl        = [];
        handles.first_correl    = true;
        handles.y_max           = 1;
    end

%% ----------------------------------------------------------------------------

    function Assimilation_selection(~, eventdata)
        % Function is called whenever a radio button has been selected, it sets
        % the global filter variable
        
        % eventdata refers to the data in the GUI when a radio button in the
        % group is changed
        
        % Set the filter_type string to newest radiobutton Value
        
        handles.filter_type = get(eventdata.NewValue,'String');
    end

%% ----------------------------------------------------------------------------

    function input_error(mystring)
        switch(lower(mystring))
            case 'observation'
                set(handles.ui_edit_observation,      'BackgroundColor', atts.red);
                set(handles.h_ens_member,             'Visible', 'Off');
                set(handles.h_best_fit,               'Visible', 'Off');
                set(handles.h_joint_update,           'Visible', 'Off');
                set(handles.h_joint_inc,              'Visible', 'Off');
                set(handles.ui_text_error,            'String' , 'Observation must be a number between 0 and 10');
                set(handles.ui_text_error,            'Visible', 'On');
                
                % Disable other input to guarantee only one error at a time!
                set(handles.ui_edit_obs_error_sd,      'Enable', 'off');
                set(handles.ui_button_create_ensemble, 'Enable', 'off');
                set(handles.ui_button_update_ensemble, 'Enable', 'off');
                
            otherwise
                
                set(handles.ui_edit_obs_error_sd, 'BackgroundColor', atts.red);
                set(handles.h_ens_member,   'Visible', 'Off');
                set(handles.h_best_fit,     'Visible', 'Off');
                set(handles.h_joint_update, 'Visible', 'Off');
                set(handles.h_joint_inc,    'Visible', 'Off');
                set(handles.ui_text_error,  'String' , 'Observation Error SD must be a number greater than 0');
                set(handles.ui_text_error,  'Visible', 'On');
                
                % Disable other input to guarantee only one error at a time!
                set(handles.ui_edit_observation,       'Enable', 'off')
                set(handles.ui_button_create_ensemble, 'Enable', 'off')
                set(handles.ui_button_update_ensemble, 'Enable', 'off')
        end   % of switch
    end  % function input_error

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
