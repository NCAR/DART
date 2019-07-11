function oned_ensemble
%% ONED_ENSEMBLE explore the details of ensemble data assimilation for a scalar.
%
%      Push on the 'Create New Ensemble' button to activate the interactive
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
% See also: gaussian_product.m oned_model.m oned_model_inf.m
%           twod_ensemble.m run_lorenz_63.m run_lorenz_96.m run_lorenz_96_inf.m

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

help oned_ensemble

atts = stylesheet; % get the default fonts and colors

% Set random number seed to same value to generate known sequences
rng('default')

% Initialize the basic structure we will be using to hold everything.
handles.ens_size         = 0;
handles.ens_members      = 0;
handles.h_obs_plot       = [];
handles.h_update_ens     = [];
handles.h_ens_member     = [];
handles.h_obs_ast        = [];
handles.h_update_lines   = [];
handles.plot_inflation   = false;
handles.h_inf_ens_member = [];
handles.h_inf_up_ens     = [];
handles.h_inf_lines      = [];
handles.h_inf_axis       = [];

%% -----------------------------------------------------------------------------

% Specify the figure size in pixels. After that, all positions are
% specified as fractions (units=Normalized). That way, the objects
% scale proportionally as the figure gets resized.
figXmin   = 450; % The horizontal position of the entire figure, in pixels
figYmin   = 250; % The vertical   position of the entire figure, in pixels
figWidth  = 670; % The width  of the entire figure, in pixels
figHeight = 500; % The height of the entire figure, in pixels

handles.figure1 = figure('Position', [figXmin figYmin figWidth figHeight], ...
    'Units', 'Pixels', ...
    'Name', 'oned_ensemble', ...
    'Color', atts.background);

%% -----------------------------------------------------------------------------
%Position has units of normalized and therefore the components must be
%Fractions of figure. Also, the actual text  has FontUnits which are
%normalized, so the Font Size is a fraction of the text box/ edit box that
%contains it. Making the font size and box size normalized allows the text
%size and box size to change proportionately if the user resizes the
%figure window. This is true for all text/edit boxes and buttons.

%% -----------------------------------------------------------------------------
%  Set up a parent container so we can move the one container around instead of
%  trying to manipulate the positions of all the components.

handles.observation_panel = uipanel('BackgroundColor',atts.red, ...
    'BorderType','none', ...
    'Units', 'Normalized', ...
    'Position',[0.66  0.744 0.325 0.2]);

handles.ui_text_observation = uicontrol(handles.observation_panel, ...
    'Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.035 0.600 0.600 0.275], ...
    'String', 'Observation' , ...
    'BackgroundColor', atts.red, ...
    'ForegroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.6, ...
    'FontWeight', 'Bold', ...
    'HorizontalAlignment', 'center');

handles.ui_edit_observation = uicontrol(handles.observation_panel, ...
    'Style', 'edit', ...
    'Units', 'Normalized', ...
    'Position', [0.675 0.562 0.300 0.373], ...
    'String', '1', ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.6, ...
    'Callback', @edit_observation_Callback);
handles.observation = str2double(get(handles.ui_edit_observation,'String'));

handles.ui_text_obs_error_sd = uicontrol(handles.observation_panel, ...
    'Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.035 0.100 0.600 0.275], ...
    'String', 'Obs. Error SD', ...
    'BackgroundColor', atts.red,...
    'ForegroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.6, ...
    'FontWeight','Bold', ...
    'HorizontalAlignment', 'center');

handles.ui_edit_obs_error_sd = uicontrol(handles.observation_panel, ...
    'Style', 'edit', ...
    'Units', 'Normalized', ...
    'Position', [0.675 0.062 0.300 0.373], ...
    'String', '1', ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.6, ...
    'Callback', @edit_obs_error_sd_Callback);
handles.obs_error_sd = str2double(get(handles.ui_edit_obs_error_sd,'String'));

%% -----------------------------------------------------------------------------
%  Try to center the boxes with what is on top.

box = get(handles.observation_panel,'Position');
center = box(1) + box(3)/2.0;
mywid = 0.250;
myleft = center - mywid/2.0;

handles.ui_button_create_new_ens = uicontrol('Style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Position', [myleft 0.645 mywid 0.075], ...
    'String', 'Create New Ensemble', ...
    'BackgroundColor', 'White', ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontSize', 0.4, ...
    'Callback', @button_create_new_ens_Callback);

mywid = 0.200;
myleft = center - mywid/2.0;
handles.ui_button_update_ens = uicontrol('style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Enable','Off', ...
    'Position', [myleft 0.549 mywid 0.075], ...
    'String', 'Update Ensemble' , ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.4, ...
    'Callback', @button_update_ens_Callback);

%% -----------------------------------------------------------------------------

handles.inflation_panel = uipanel('BackgroundColor',atts.lightblue, ...
    'BorderType','none', ...
    'Position',[0.661 0.282 0.32 0.248]);

handles.ui_checkbox_inflation = uicontrol(handles.inflation_panel, ...
    'Style', 'checkbox', ...
    'Units', 'Normalized', ...
    'Position', [0.037 0.624 0.900 0.335], ...
    'String', 'Apply Inflation', ...
    'BackgroundColor', atts.lightblue,...
    'ForegroundColor', 'k', ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontSize', 0.4, ...
    'FontWeight','Bold', ...
    'Callback', @inflation_toggle_Callback);

handles.ui_slider_inflation = uicontrol(handles.inflation_panel, ...
    'Style', 'slider', ...
    'Units', 'Normalized', ...
    'Position', [0.058 0.411 0.893 0.18], ...
    'Value', 1,...
    'Max', 5,...
    'Min', 1,...
    'Sliderstep',[0.05 0.2], ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontSize', 0.5, ...
    'Enable', 'Off', ...
    'Callback', @slider_Callback);

handles.ui_text_inflation = uicontrol(handles.inflation_panel, ...
    'Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.018 0.006 0.671 0.29], ...
    'String', 'Inflation Amount', ...
    'BackgroundColor', atts.lightblue,...
    'ForegroundColor', 'k', ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontSize', 0.4, ...
    'Enable', 'Off', ...
    'FontWeight','Bold');

handles.ui_edit_inflation_label = uicontrol(handles.inflation_panel, ...
    'Style', 'edit', ...
    'Units', 'Normalized', ...
    'Position', [0.718 0.052 0.263 0.258], ...
    'String', get(handles.ui_slider_inflation,'Value'), ...
    'BackgroundColor', 'White', ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontSize', 0.5, ...
    'Enable', 'Off', ...
    'Callback', @edit_inflation_Callback);
handles.inflation = str2double(get(handles.ui_edit_inflation_label,'String'));

%% -----------------------------------------------------------------------------
%  doesn't seem to do a great job of centering ... but the distribute is nice.

hlist = [handles.observation_panel,   ...
    handles.ui_button_create_new_ens,  ...
    handles.ui_button_update_ens,  ...
    handles.inflation_panel];

align(hlist,'Center','Distribute')

%% -----------------------------------------------------------------------------

handles.ui_radio_button_group = uibuttongroup ('BackgroundColor', atts.background, ...
    'BorderType', 'none', ...
    'Position',[0.776 0.013 0.200 0.256]);

handles.ui_radio_button_eakf = uicontrol(handles.ui_radio_button_group, ...
    'Style', 'radio button', ...
    'Units', 'Normalized', ...
    'Position', [0.007 0.667 670 0.253], ...
    'String', 'EAKF', ...
    'BackgroundColor', atts.background, ...
    'Foreground', 'Black', ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontWeight','Bold', ...
    'FontSize', 0.5, ...
    'HandleVisibility', 'On');

handles.ui_radio_button_enkf = uicontrol(handles.ui_radio_button_group, ...
    'Style', 'radio button', ...
    'Units', 'Normalized', ...
    'Position', [0.007 0.347 670 0.253], ...
    'String', 'EnKF', ...
    'BackgroundColor', atts.background, ...
    'Foreground', 'Black', ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontWeight','Bold', ...
    'FontSize', 0.5, ...
    'HandleVisibility', 'Off');

handles.ui_radio_button_rhf = uicontrol(handles.ui_radio_button_group, ...
    'Style', 'radio button', ...
    'Units', 'Normalized', ...
    'Position', [0.007 0.06 670 0.253], ...
    'String', 'RHF', ...
    'BackgroundColor', atts.background, ...
    'Foreground', 'Black', ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontWeight','Bold', ...
    'FontSize', 0.5, ...
    'HandleVisibility', 'Off');

%% -----------------------------------------------------------------------------

handles.axes = axes ('Units', 'Normalized', ...
    'Position', [30/figWidth 30/figWidth 0.6000 0.9000], ...
    'FontName', atts.fontname, ...
    'FontSize', atts.fontsize, ...
    'Color','White');

% This section specifies the annotation for the values on the main graph.
% By specifying them this way, we can turn them On/Off at will and specify
% that they scale when the object is resized.
% FIXME ... these should be part of some uipanel that is on the graphic.

handles.ui_text_prior_mean = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.052 0.852 0.233 0.065], ...
    'String', 'Prior Mean = ', ...
    'BackgroundColor', 'White', ...
    'ForegroundColor', atts.green, ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontSize', 0.4, ...
    'FontWeight','Bold', ...
    'HorizontalAlignment', 'right');

handles.ui_text_inflated_prior_mean = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.052 0.806 0.233 0.065], ...
    'String', 'Inflated =      ', ...
    'Visible', 'Off', ...
    'BackgroundColor', 'White', ...
    'ForegroundColor', atts.green, ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontSize', 0.4, ...
    'FontWeight','Bold', ...
    'HorizontalAlignment','right');

handles.ui_text_prior_sd = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.052 0.762 0.233 0.065], ...
    'String', 'Prior SD = ', ...
    'BackgroundColor', 'White', ...
    'ForegroundColor', atts.green, ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontSize', 0.4, ...
    'FontWeight','Bold', ...
    'HorizontalAlignment', 'right');

handles.ui_text_inflated_prior_sd = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.052 0.713 0.233 0.065], ...
    'String', 'Inflated = ', ...
    'Visible', 'Off', ...
    'BackgroundColor', 'White', ...
    'ForegroundColor',atts.green,...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontSize', 0.4, ...
    'FontWeight','Bold', ...
    'HorizontalAlignment','right');

handles.ui_text_post_mean = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.340 0.852 0.270 0.065], ...
    'String', 'Posterior Mean = ', ...
    'BackgroundColor', 'White', ...
    'ForegroundColor', atts.blue, ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontSize', 0.4, ...
    'FontWeight','Bold', ...
    'HorizontalAlignment','right');

handles.ui_text_inflated_post_mean = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.340 0.806 0.270 0.065], ...
    'String', 'Inflated = ', ...
    'Visible', 'Off', ...
    'BackgroundColor', 'White', ...
    'ForegroundColor', atts.blue, ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontSize', 0.4, ...
    'FontWeight','Bold', ...
    'HorizontalAlignment','right');

handles.ui_text_post_sd = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.340 0.762 0.270 0.065], ...
    'String', 'Posterior SD = ', ...
    'BackgroundColor', 'White', ...
    'ForegroundColor', atts.blue, ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontSize', 0.4, ...
    'FontWeight','Bold', ...
    'HorizontalAlignment','right');

handles.ui_text_inflated_post_sd = uicontrol('Style', 'text', ...
    'Units'              , 'Normalized', ...
    'Position'           , [0.340 0.713 0.270 0.065], ...
    'String'             , 'Inflated = ', ...
    'Visible'            , 'Off', ...
    'BackgroundColor'    , 'White', ...
    'ForegroundColor'    , atts.blue, ...
    'FontUnits'          , 'normalized', ...
    'FontName'           , atts.fontname, ...
    'FontSize'           , 0.4, ...
    'FontWeight'         , 'Bold', ...
    'HorizontalAlignment', 'right');

% justify all the strings:

set(handles.ui_text_prior_mean          , 'String', 'Prior Mean =       ');
set(handles.ui_text_prior_sd            , 'String',   'Prior SD =       ');
set(handles.ui_text_inflated_prior_mean , 'String',   'Inflated =       ');
set(handles.ui_text_inflated_prior_sd   , 'String',   'Inflated =       ');

set(handles.ui_text_post_mean           , 'String', 'Posterior Mean =       ');
set(handles.ui_text_post_sd             , 'String',   'Posterior SD =       ');
set(handles.ui_text_inflated_post_mean  , 'String',       'Inflated =       ');
set(handles.ui_text_inflated_post_sd    , 'String',       'Inflated =       ');

% FIXME These are error messages that may or may not be present.
% We create them and then turn them off. When we need to, IF we need to,
% they can be turned back on.

handles.ui_text_inf_err_print = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [40/figWidth 260/figHeight 400/figWidth 40/figHeight], ...
    'String', 'ERROR: Inflation value must be between 1 and 5.', ...
    'BackgroundColor', 'White', ...
    'ForegroundColor', atts.red, ...
    'FontSize', atts.fontsize, ...
    'FontWeight', 'Bold', ...
    'FontName', atts.fontname, ...
    'Visible', 'Off');

handles.ui_text_obs_err_print = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [40/figWidth 260/figHeight 400/figWidth 40/figHeight], ...
    'String', 'ERROR: Observation value must be numeric.', ...
    'BackgroundColor', 'White', ...
    'ForegroundColor', atts.red, ...
    'FontSize', atts.fontsize, ...
    'FontWeight', 'Bold', ...
    'FontName', atts.fontname, ...
    'Visible', 'Off');

handles.ui_text_obs_sd_err_print = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [40/figWidth 260/figHeight 400/figWidth 40/figHeight], ...
    'String', 'ERROR: Obs. Error SD value must be numeric.', ...
    'BackgroundColor', 'White', ...
    'ForegroundColor', atts.red, ...
    'FontSize', atts.fontsize, ...
    'FontWeight', 'Bold', ...
    'FontName', atts.fontname, ...
    'Visible','Off');

% Go ahead and plot the initial observational error distribution

handles.h_obs_plot = plot_gaussian(handles.observation, handles.obs_error_sd, 1);
set(handles.h_obs_plot, 'Color', atts.red, 'Linestyle', '--', 'Linewidth', 2.0);
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
set(gca, 'FontSize', atts.fontsize)
set(gca, 'XGrid', 'on')

plot([xlower xupper], [0 0], 'k', 'Linewidth', 1.7);
title('oned_ensemble','Interpreter','none')

%% -----------------------------------------------------------------------------

    function button_create_new_ens_Callback(~,~)
        
        % Disable the update ensemble button and all other active buttons
        set(handles.ui_button_update_ens,     'Enable', 'Off');
        set(handles.ui_edit_observation,      'Enable', 'Off');
        set(handles.ui_edit_obs_error_sd,     'Enable', 'Off');
        set(handles.ui_edit_inflation_label,  'Enable', 'Off');
        
        % Clear out any old ensemble members if they exist
        set(handles.h_ens_member,          'Visible', 'Off');
        set(handles.h_inf_ens_member,      'Visible', 'Off');
        
        set(handles.h_update_lines,        'Visible', 'Off');
        set(handles.h_inf_lines,           'Visible', 'Off');
        set(handles.h_inf_axis,            'Visible', 'Off');
        
        % Turn Off any old update points
        set(handles.h_update_ens,          'Visible', 'Off');
        set(handles.h_inf_up_ens,          'Visible', 'Off');
        set(handles.h_inf_ens_member,      'Visible', 'Off');
        
        clear_ui_labels;
        
        hold on
        
        % Set a basic plotting domain range that includes mean +/- 3 obs SDs
        xlower = min(handles.observation - 3*handles.obs_error_sd, min(handles.ens_members));
        xupper = max(handles.observation + 3*handles.obs_error_sd, max(handles.ens_members));
        ylower = -0.4;
        yupper = 1.0;
        axis([xlower xupper ylower yupper]);
        
        set(gca, 'YTick',      [0 0.2 0.4 0.6 0.8]);
        
        % Messages are centered in the middle.
        xmid = (xupper + xlower) / 2.0;
        h_click    = text(xmid,  0.6, {'Click inside graphics box to create member', ...
            '(only X value is used)'}, 'FontSize', atts.fontsize, 'HorizontalAlignment', 'center');
        
        h_err_text = text(xmid, -0.15, 'An ensemble has to have at least 2 members.', ...
            'FontSize', atts.fontsize, 'Visible', 'on', 'HorizontalAlignment', 'center','Color', atts.red);
        
        h_finish   = text(xmid, -0.15, 'Click outside of plot to finish', ...
            'Fontsize', atts.fontsize, 'Visible', 'Off', 'HorizontalAlignment', 'center');
        
        ens_size = 0;
        
        while ens_size < 100
            [xt, yt] = ginput(1);
            
            if(xt >= xlower && xt <= xupper && yt >= ylower && yt <= yupper)
                ens_size = ens_size + 1;
                x(ens_size) = xt; %#ok<AGROW>
                y(ens_size) = 0; %#ok<AGROW>
                handles.h_ens_member(ens_size) = ...
                    plot(x(ens_size), y(ens_size), '*', 'MarkerSize', 16, 'Color', atts.green,'LineWidth',2.0);
                
                % Display the prior mean and sd
                prior_mean = mean(x);
                prior_sd   = std(x);
                str1 = sprintf('Prior Mean = %.4f', prior_mean);
                set(handles.ui_text_prior_mean, 'String', str1);
                str1 = sprintf('Prior SD = %.4f', prior_sd);
                set(handles.ui_text_prior_sd,   'String', str1);
                
            elseif (ens_size < 2)
                set(h_err_text,'FontWeight','bold')
                
            else
                break;
                
            end
            
            % Swap messages once you have a minimal ensemble.
            if (ens_size == 2)
                set(h_err_text, 'Visible', 'Off');
                set(h_finish,   'Visible', 'on');
                
            end
            
        end
        
        % Ensemble created, compute mean and sd, clean up and return
        % Set the global gui storage
        handles.ens_size    = ens_size;
        handles.ens_members = x;
        
        % Turn Off the data entry messages
        set(h_click,  'Visible', 'Off');
        set(h_finish, 'Visible', 'Off');
        
        % Enable the update ensemble button
        set(handles.ui_button_update_ens,     'Enable', 'On');
        set(handles.ui_edit_observation,      'Enable', 'On');
        set(handles.ui_edit_obs_error_sd,     'Enable', 'On');
        set(handles.ui_edit_inflation_label,  'Enable', 'On');
        
    end


%% -----------------------------------------------------------------------------

    function inflation_toggle_Callback (~, ~)
        
        enabled = get(handles.ui_checkbox_inflation, 'Value');
        if (enabled)
            set(handles.ui_slider_inflation,         'Enable',  'On');
            set(handles.ui_text_inflation,           'Enable',  'On');
            set(handles.ui_edit_inflation_label,     'Enable',  'On');
            set(handles.ui_text_inflated_prior_mean, 'Visible', 'On');
            set(handles.ui_text_inflated_post_mean,  'Visible', 'On');
            set(handles.ui_text_inflated_prior_sd,   'Visible', 'On');
            set(handles.ui_text_inflated_post_sd,    'Visible', 'On');
            
        else
            set(handles.ui_slider_inflation,         'Enable',  'Off');
            set(handles.ui_text_inflation,           'Enable',  'Off');
            set(handles.ui_edit_inflation_label,     'Enable',  'Off');
            set(handles.ui_text_inflated_prior_mean, 'Visible', 'Off');
            set(handles.ui_text_inflated_post_mean,  'Visible', 'Off');
            set(handles.ui_text_inflated_prior_sd,   'Visible', 'Off');
            set(handles.ui_text_inflated_post_sd,    'Visible', 'Off');
        end
    end

%% -----------------------------------------------------------------------------

    function slider_Callback (~, ~)
        
        handles.inflation = get(handles.ui_slider_inflation, 'Value');
        
        str1 = sprintf('%.4f',handles.inflation);
        set(handles.ui_edit_inflation_label, 'String', str1);
        
        % Just in case the inflation label was in the error state, reset
        set(handles.ui_edit_inflation_label, 'BackgroundColor', 'White', 'FontWeight', 'Normal');
        set(handles.ui_text_inf_err_print, 'Visible', 'Off')
        
        % Disable other input to guarantee only one error at a time!
        set(handles.ui_edit_observation,      'Enable', 'On')
        set(handles.ui_edit_obs_error_sd,     'Enable', 'On')
        set(handles.ui_button_create_new_ens, 'Enable', 'On')
        set(handles.ui_button_update_ens,     'Enable', 'On')
        
    end

%% -----------------------------------------------------------------------------

    function button_update_ens_Callback (~, ~)
        
        % Turn Off any old points
        set(handles.h_update_ens,     'Visible', 'Off');
        set(handles.h_inf_up_ens,     'Visible', 'Off');
        set(handles.h_inf_ens_member, 'Visible', 'Off');
        
        % Remove mean and sd of old posterior
        clear_ui_labels;
        
        % And the lines in between
        set(handles.h_update_lines, 'Visible', 'Off');
        set(handles.h_inf_lines,    'Visible', 'Off');
        set(handles.h_inf_axis,     'Visible', 'Off');
        
        ensemble = handles.ens_members;
        
        % Figure out which filter option is currently selected
        val = get(handles.ui_radio_button_group,'SelectedObject');
        filter_type = get(val,'String');
        
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
        handles.h_update_ens = plot(new_ensemble, y, '*', 'MarkerSize', 16, 'Color', atts.blue);
        
        % Plot lines connecting the prior and posterior ensemble members
        for i = 1:size(ensemble, 2)
            x_line = [handles.ens_members(i), new_ensemble(i)];
            y_line = [0, -0.1];
            handles.h_update_lines(i) = plot(x_line, y_line, 'k');
        end
        
        % Add in a label of the updated mean and sd
        new_mean = mean(new_ensemble);
        new_sd   = std(new_ensemble);
        
        % Update mean and sd of old posterior
        str1 = sprintf('Posterior Mean = %.4f',new_mean);
        set(handles.ui_text_post_mean, 'String', str1, 'Visible', 'on');
        
        str1 = sprintf('Posterior SD = %.4f',new_sd);
        set(handles.ui_text_post_sd,   'String', str1, 'Visible', 'on');
        
        % If the checkbox isn't set, return now
        if(not(get(handles.ui_checkbox_inflation, 'Value')))
            return
        end
        
        % Plot the inflated prior ensemble
        y = -0.2;
        handles.prior_mean = mean(handles.ens_members(1:handles.ens_size));
        
        inf_ens = zeros(1,handles.ens_size);
        
        for i = 1: handles.ens_size
            inf_ens(i) = (handles.ens_members(i) - handles.prior_mean) * sqrt(handles.inflation) + ...
                handles.prior_mean;
            handles.h_inf_ens_member(i) = plot(inf_ens(i), y, '*', 'MarkerSize', 16, 'Color', atts.green,'LineWidth',2.0);
            
        end
        
        % Update mean and sd of old posterior
        handles.inf_prior_sd = std(inf_ens(1:handles.ens_size));
        
        str1 = sprintf('Inflated = %.4f',handles.prior_mean);
        set(handles.ui_text_inflated_prior_mean,'String',str1,'Visible','on');
        str1 = sprintf('Inflated = %.4f',handles.inf_prior_sd);
        set(handles.ui_text_inflated_prior_sd,  'String',str1,'Visible','on');
        
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
        handles.h_inf_up_ens = plot(new_ensemble, y, '*', 'MarkerSize', 16, 'Color', atts.blue);
        
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
        plot([xlower xupper], [0 0], 'k', 'Linewidth', 1.7);
        handles.h_inf_axis = plot([xlower xupper], [-0.2 -0.2], 'k', 'Linewidth', 1.7);
        
        % Update mean and sd of old posterior
        handles.update_inf_mean = mean(new_ensemble(1:handles.ens_size));
        handles.update_inf_sd =   std (new_ensemble(1:handles.ens_size));
        
        str1 = sprintf('Inflated = %.4f',handles.update_inf_mean);
        set(handles.ui_text_inflated_post_mean, 'String', str1, 'Visible','on');
        
        str1 = sprintf('Inflated = %.4f',handles.update_inf_sd);
        set(handles.ui_text_inflated_post_sd, 'String', str1, 'Visible', 'on');
        
    end

%% -----------------------------------------------------------------------------

    function clear_ui_labels()
        
        % Turns Off all labels except for the prior mean and SD
        set(handles.ui_text_post_sd,             'Visible', 'Off');
        set(handles.ui_text_post_mean,           'Visible', 'Off');
        set(handles.ui_text_inflated_prior_mean, 'Visible', 'Off');
        set(handles.ui_text_inflated_prior_sd,   'Visible', 'Off');
        set(handles.ui_text_inflated_post_sd,    'Visible', 'Off');
        set(handles.ui_text_inflated_post_mean,  'Visible', 'Off');
        
    end

%% -----------------------------------------------------------------------------

    function edit_inflation_Callback(~, ~)
        
        % Turn Off any old updated points
        set(handles.h_update_ens,          'Visible', 'Off');
        set(handles.h_inf_up_ens,          'Visible', 'Off');
        set(handles.h_inf_ens_member,      'Visible', 'Off');
        
        % Remove mean and sd of old posterior
        clear_ui_labels;
        
        % And the lines in between
        set(handles.h_update_lines,        'Visible', 'Off');
        set(handles.h_inf_lines,           'Visible', 'Off');
        set(handles.h_inf_axis,            'Visible', 'Off');
        
        % Enable things that an error might have turned Off
        set(handles.ui_edit_observation,      'Enable', 'on')
        set(handles.ui_edit_obs_error_sd,     'Enable', 'on')
        set(handles.ui_button_create_new_ens, 'Enable', 'on')
        
        % Only enable the update ensemble pushbutton if an ensemble has been created
        if(handles.ens_size > 0)
            set(handles.ui_button_update_ens, 'Enable', 'on');
            
        end
        
        % Get the value of the inflation
        inf_value = str2double(get(handles.ui_edit_inflation_label, 'String'));
        
        if( isfinite(inf_value) && (inf_value >= 1) && (inf_value <= 5))
            inflation = inf_value;
            
        else
            set(handles.ui_edit_inflation_label, 'String', '??','FontWeight','Bold', ...
                'BackgroundColor', atts.red);
            set(handles.ui_text_inf_err_print,'Visible','On')
            
            fprintf('ERROR: Inflation value must be between 1 and 5.\n')
            fprintf('ERROR: Inflation value must be between 1 and 5.\n')
            
            % Disable other input to guarantee only one error at a time!
            set(handles.ui_edit_observation,      'Enable', 'Off')
            set(handles.ui_edit_obs_error_sd,     'Enable', 'Off')
            set(handles.ui_button_create_new_ens, 'Enable', 'Off')
            set(handles.ui_button_update_ens,     'Enable', 'Off')
            
            return
            
        end
        
        % Update the value in global storage
        handles.inflation = inflation;
        set(handles.ui_edit_inflation_label, 'BackgroundColor', 'White', 'FontWeight', 'Normal');
        set(handles.ui_slider_inflation,'Value', handles.inflation);
        set(handles.ui_text_inf_err_print,'Visible','Off')
        
        
        % Plot the updated distribution
        set(handles.h_obs_plot, 'Visible', 'Off');
        handles.h_obs_plot = plot_gaussian(handles.observation, handles.obs_error_sd, 1);
        set(handles.h_obs_plot, 'Color', atts.red, 'Linestyle', '--', 'Linewidth', 1.7);
        
        % Set a basic plotting domain range that includes mean +/- 3 obs SDs
        xlower = min(handles.observation - 3*handles.obs_error_sd, min(handles.ens_members));
        xupper = max(handles.observation + 3*handles.obs_error_sd, max(handles.ens_members));
        ylower = -0.4;
        yupper = 1.0;
        axis([xlower xupper ylower yupper]);
        
        set(handles.h_obs_plot, 'Color', atts.red, 'Linestyle', '--', 'Linewidth', 1.7);
        
        set(gca, 'YTick', [0 0.2 0.4 0.6 0.8]);
        
        hold on
        
        plot([xlower xupper], [0 0], 'k', 'Linewidth', 1.7);
        
    end

%% -----------------------------------------------------------------------------

    function edit_observation_Callback(~, ~)
        
        % Turn Off any old updated points
        set(handles.h_update_ens,     'Visible', 'Off');
        set(handles.h_inf_up_ens,     'Visible', 'Off');
        set(handles.h_inf_ens_member, 'Visible', 'Off');
        
        % Remove mean and sd of old posterior
        clear_ui_labels;
        
        % And the lines in between
        set(handles.h_update_lines,        'Visible', 'Off');
        set(handles.h_inf_lines,           'Visible', 'Off');
        set(handles.h_inf_axis,            'Visible', 'Off');
        
        % Enable things that an error might have turned Off
        set(handles.ui_edit_obs_error_sd,     'Enable', 'on')
        set(handles.ui_edit_inflation_label,  'Enable', 'on')
        set(handles.ui_button_create_new_ens, 'Enable', 'on')
        
        % Only enable the update ensemble pushbutton if an ensemble has been created
        if(handles.ens_size > 0)
            set(handles.ui_button_update_ens, 'Enable', 'on');
            
        end
        
        % Get the value of the observation
        if(isfinite(      str2double(get(handles.ui_edit_observation, 'String'))))
            observation = str2double(get(handles.ui_edit_observation, 'String'));
            
        else
            set(handles.ui_edit_observation, 'String', '??','FontWeight','Bold', ...
                'BackgroundColor', atts.red);
            set(handles.ui_text_obs_err_print,'Visible','On')
            
            fprintf('ERROR: Observation value must be numeric.\n')
            fprintf('ERROR: Observation value must be numeric.\n')
            
            
            % Disable other input to guarantee only one error at a time!
            set(handles.ui_edit_obs_error_sd,     'Enable', 'Off')
            set(handles.ui_edit_inflation_label,  'Enable', 'Off')
            set(handles.ui_button_create_new_ens, 'Enable', 'Off')
            set(handles.ui_button_update_ens,     'Enable', 'Off')
            
            return
            
        end
        
        % Update the global storage
        handles.observation = observation;
        set(handles.ui_edit_observation, 'BackgroundColor', 'White','FontWeight', 'Normal');
        set(handles.ui_text_obs_err_print,'Visible','Off')
        
        % Plot the updated distribution
        set(handles.h_obs_plot, 'Visible', 'Off');
        handles.h_obs_plot = plot_gaussian(handles.observation, handles.obs_error_sd, 1);
        set(handles.h_obs_plot, 'Color', atts.red, 'Linestyle', '--', 'Linewidth', 1.7);
        
        % Move the observation asterisk
        set(handles.h_obs_ast, 'Visible', 'Off');
        handles.h_obs_ast = plot(handles.observation, 0, 'r*', 'MarkerSize', 16,'LineWidth',2.0);
        
        % Set a basic plotting domain range that includes mean +/- 3 obs SDs
        xlower = min(handles.observation - 3*handles.obs_error_sd, min(handles.ens_members));
        xupper = max(handles.observation + 3*handles.obs_error_sd, max(handles.ens_members));
        ylower = -0.4;
        yupper = 1.0;
        axis([xlower xupper ylower yupper]);
        
        set(gca, 'YTick', [0 0.2 0.4 0.6 0.8]);
        
        hold on
        plot([xlower xupper], [0 0], 'k', 'Linewidth', 1.7);
        
    end

%% -----------------------------------------------------------------------------

    function edit_obs_error_sd_Callback(~, ~)
        
        % Turn Off any old updated points
        set(handles.h_update_ens,          'Visible', 'Off');
        set(handles.h_inf_up_ens,          'Visible', 'Off');
        set(handles.h_inf_ens_member,      'Visible', 'Off');
        
        % Remove mean and sd of old posterior
        clear_ui_labels;
        
        % And the lines in between
        set(handles.h_update_lines,        'Visible', 'Off');
        set(handles.h_inf_lines,           'Visible', 'Off');
        set(handles.h_inf_axis,            'Visible', 'Off');
        
        % Enable things that an error might have turned Off
        set(handles.ui_edit_observation,      'Enable', 'on')
        set(handles.ui_edit_inflation_label,  'Enable', 'on')
        set(handles.ui_button_create_new_ens, 'Enable', 'on')
        
        % Only enable the update ensemble pushbutton if an ensemble has been created
        if(handles.ens_size > 0)
            set(handles.ui_button_update_ens, 'Enable', 'on');
        end
        
        % Get the value of the observation error sd
        obs_error_value = str2double(get(handles.ui_edit_obs_error_sd, 'String'));
        
        if(isfinite(obs_error_value) && (obs_error_value > 0))
            obs_error_sd = obs_error_value;
            
        else
            
            set(handles.ui_edit_obs_error_sd, 'String', '??','FontWeight','Bold', ...
                'BackgroundColor', atts.red);
            set(handles.ui_text_obs_sd_err_print,'Visible','On')
            
            fprintf('ERROR: Obs. Error SD value must be numeric.\n')
            fprintf('ERROR: Obs. Error SD value must be numeric.\n')
            
            
            % Disable other input to guarantee only one error at a time!
            set(handles.ui_edit_observation,      'Enable', 'Off')
            set(handles.ui_edit_inflation_label,  'Enable', 'Off')
            set(handles.ui_button_create_new_ens, 'Enable', 'Off')
            set(handles.ui_button_update_ens,     'Enable', 'Off')
            
            return
            
        end
        
        % Update the value in global storage
        handles.obs_error_sd = obs_error_sd;
        set(handles.ui_edit_obs_error_sd, 'BackgroundColor', 'White',  'FontWeight', 'Normal');
        set(handles.ui_text_obs_sd_err_print,'Visible','Off')
        
        
        % Plot the updated distribution
        set(handles.h_obs_plot, 'Visible', 'Off');
        handles.h_obs_plot = plot_gaussian(handles.observation, handles.obs_error_sd, 1);
        set(handles.h_obs_plot, 'Color', atts.red, 'Linestyle', '--', 'Linewidth', 1.7);
        
        % Set a basic plotting domain range that includes mean +/- 3 obs SDs
        xlower = min(handles.observation - 3*handles.obs_error_sd, min(handles.ens_members));
        xupper = max(handles.observation + 3*handles.obs_error_sd, max(handles.ens_members));
        ylower = -0.4;
        yupper = 1.0;
        axis([xlower xupper ylower yupper]);
        
        set(handles.h_obs_plot, 'Color', atts.red, 'Linestyle', '--', 'Linewidth', 1.7);
        
        set(gca, 'YTick', [0 0.2 0.4 0.6 0.8]);
        
        hold on
        
        plot([xlower xupper], [0 0], 'k', 'Linewidth', 1.7);
        
    end

%% -----------------------------------------------------------------------------

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
