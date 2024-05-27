function oned_cycle
%% ONED_CYCLEthe details of ensemble data assimilation for a scalar.
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
%      also printed on the plot. For the EAKF and RHF, the continuous prior distribution
%      that best fits the prior ensemble (green) and the posterior continuous 
%      ensemble from which the QCEFF algorithm determines the posterior ensemble
%      members (blue) are also plotted. The EnKF is not a QCEFF filter and does
%      not make use of a continuous prior distribution fit.  
%
%      The type of ensemble Kalman filter update can be chosen using the
%      pulldown menu at the bottom.
%
%      The 'EnKF' is a stochastic algorithm so repeated updates can be done
%      for the same prior and observation.
%
%      The mean and standard deviation of the likelihood
%      can be changed in the red box. Change the Observation Error SD, lay down 
%      an ensemble pretty far away from the observation - have fun with it.
%
% See also: bounded_oned_ensemble.m gaussian_product.m oned_model.m oned_model_inf.m
%           twod_ensemble.m run_lorenz_63.m run_lorenz_96.m run_lorenz_96_inf.m

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

help oned_cycle

atts = stylesheet; % get the default fonts and colors

% Set random number seed to same value to generate known sequences
rng('default')

% Initialize the basic structure we will be using to hold everything.
handles.ens_size         = 0;
handles.ens_members      = 0;
handles.h_obs_plot       = [];
handles.h_update_ens     = [];
handles.h_prior_pdf      = [];
handles.h_post_pdf       = [];
handles.h_prior_cont_pdf = [];
handles.h_post_cont_pdf  = [];
handles.h_ens_member     = [];
handles.h_obs_ast        = [];
handles.h_update_lines   = [];

handles.prior_cont_mean = 1;
handles.prior_cont_sd   = 2;

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
    'Name', 'oned_cycle', ...
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

handles.ui_text_growth_rate = uicontrol(handles.observation_panel, ...
    'Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.035 0.600 0.600 0.275], ...
    'String', 'Growth Rate' , ...
    'BackgroundColor', atts.red, ...
    'ForegroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.6, ...
    'FontWeight', 'Bold', ...
    'HorizontalAlignment', 'center');

handles.ui_edit_growth_rate= uicontrol(handles.observation_panel, ...
    'Style', 'edit', ...
    'Units', 'Normalized', ...
    'Position', [0.675 0.562 0.300 0.373], ...
    'String', '1', ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.6, ...
    'Callback', @edit_growth_rate_Callback);
handles.growth_rate = str2double(get(handles.ui_edit_growth_rate,'String'));

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

% The true observation value can be initialized to 0 for plotting
handles.observation = 0;

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
handles.ui_button_cycle_da = uicontrol('style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Enable','On', ...
    'Position', [myleft 0.549 mywid 0.075], ...
    'String', 'Cycle DA' , ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.4, ...
    'Callback', @button_cycle_da_Callback);

%% -----------------------------------------------------------------------------

%  doesn't seem to do a great job of centering ... but the distribute is nice.

hlist = [handles.observation_panel,   ...
    handles.ui_button_create_new_ens,  ...
    handles.ui_button_cycle_da];

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

handles.ui_text_prior_ens_mean = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.052 0.152 0.273 0.065], ...
    'String', 'Prior Ens Mean = ', ...
    'BackgroundColor', 'White', ...
    'ForegroundColor', atts.green, ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontSize', 0.4, ...
    'FontWeight','Bold', ...
    'HorizontalAlignment', 'right');

handles.ui_text_prior_cont_mean = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.052 0.122 0.273 0.065], ...
    'String', 'Prior Cont Mean = ', ...
    'BackgroundColor', 'White', ...
    'ForegroundColor', atts.green, ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontSize', 0.4, ...
    'FontWeight','Bold', ...
    'HorizontalAlignment', 'right');

handles.ui_text_prior_ens_sd = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.052 0.062 0.273 0.065], ...
    'String', 'Prior Ens SD = ', ...
    'BackgroundColor', 'White', ...
    'ForegroundColor', atts.green, ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontSize', 0.4, ...
    'FontWeight','Bold', ...
    'HorizontalAlignment', 'right');

handles.ui_text_prior_cont_sd = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.052 0.032 0.273 0.065], ...
    'String', 'Prior Cont SD = ', ...
    'BackgroundColor', 'White', ...
    'ForegroundColor', atts.green, ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontSize', 0.4, ...
    'FontWeight','Bold', ...
    'HorizontalAlignment', 'right');

handles.ui_text_post_ens_mean = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.340 0.152 0.300 0.065], ...
    'String', 'Posterior Ens Mean = ', ...
    'BackgroundColor', 'White', ...
    'ForegroundColor', atts.blue, ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontSize', 0.4, ...
    'FontWeight','Bold', ...
    'HorizontalAlignment','right');

handles.ui_text_post_cont_mean = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.340 0.122 0.300 0.065], ...
    'String', 'Posterior Cont Mean = ', ...
    'BackgroundColor', 'White', ...
    'ForegroundColor', atts.blue, ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontSize', 0.4, ...
    'FontWeight','Bold', ...
    'HorizontalAlignment','right');

handles.ui_text_post_ens_sd = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.340 0.062 0.300 0.065], ...
    'String', 'Posterior Ens SD = ', ...
    'BackgroundColor', 'White', ...
    'ForegroundColor', atts.blue, ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontSize', 0.4, ...
    'FontWeight','Bold', ...
    'HorizontalAlignment','right');

handles.ui_text_post_cont_sd = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.340 0.032 0.300 0.065], ...
    'String', 'Posterior Cont SD = ', ...
    'BackgroundColor', 'White', ...
    'ForegroundColor', atts.blue, ...
    'FontUnits', 'normalized', ...
    'FontName', atts.fontname, ...
    'FontSize', 0.4, ...
    'FontWeight','Bold', ...
    'HorizontalAlignment','right');

% justify all the strings, but don't turn them on yet:
set(handles.ui_text_prior_ens_mean, 'String', 'Prior Ens Mean =       ',      'visible', 'off');
set(handles.ui_text_prior_ens_sd,   'String', 'Prior Ens SD =       ',        'visible', 'off');

set(handles.ui_text_post_ens_mean,  'String', 'Posterior Ens Mean =       ',  'visible', 'off');
set(handles.ui_text_post_cont_mean, 'String', 'Posterior Cont Mean =       ', 'visible', 'off');
set(handles.ui_text_post_ens_sd,    'String', 'Posterior Ens SD =       ',    'visible', 'off');
set(handles.ui_text_post_cont_sd,   'String', 'Posterior Cont SD =       ',   'visible', 'off');

% Initial values of the continuous prior distribution
str1 = sprintf('Prior Cont Mean = %.4f', handles.prior_cont_mean);
set(handles.ui_text_prior_cont_mean, 'String', str1, 'Visible', 'on');
str1 = sprintf('Prior Cont SD = %.4f', handles.prior_cont_sd);
set(handles.ui_text_prior_cont_sd,   'String', str1, 'Visible', 'on');

% These are error messages that may or may not be present.
% We create them and then turn them off. When we need to, IF we need to,
% they can be turned back on.

handles.ui_text_growth_rate_err_print = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [40/figWidth 260/figHeight 400/figWidth 40/figHeight], ...
    'String', 'ERROR: Growth Rate must be numeric.', ...
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

% Set a basic plotting domain range that includes mean +/- 3 obs SDs
xlower = - 3*handles.obs_error_sd;
xupper = + 3*handles.obs_error_sd;
ylower = -0.4;
yupper = 1.0;
axis([xlower xupper ylower yupper]);
hold on

%set(gca, 'YTick', [0 0.2 0.4 0.6 0.8]);
set(gca, 'FontSize', atts.fontsize)
set(gca, 'XGrid', 'on')

plot([xlower xupper], [0 0], 'k', 'Linewidth', 1.7);
title('oned_cycle','Interpreter','none')

%% -----------------------------------------------------------------------------

    function button_create_new_ens_Callback(~,~)
        
        % Disable the update ensemble button and all other active buttons
        set(handles.ui_button_cycle_da,     'Enable', 'Off');
        set(handles.ui_edit_growth_rate,      'Enable', 'Off');
        set(handles.ui_edit_obs_error_sd,     'Enable', 'Off');
        
        % Clear out any old ensemble members if they exist
        set(handles.h_ens_member,          'Visible', 'Off');
        
        set(handles.h_update_lines,        'Visible', 'Off');
        
        % Turn Off any old update points
        set(handles.h_update_ens,          'Visible', 'Off');
        set(handles.h_prior_pdf,           'Visible', 'Off');
        set(handles.h_post_pdf,            'Visible', 'Off');
        set(handles.h_prior_cont_pdf,      'Visible', 'Off');
        set(handles.h_post_cont_pdf,       'Visible', 'Off');
        
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
        
        h_err_text = text(xmid, 0.35, 'An ensemble has to have at least 2 members.', ...
            'FontSize', atts.fontsize, 'Visible', 'on', 'HorizontalAlignment', 'center','Color', atts.red);
        
        h_finish   = text(xmid, 0.35, 'Click outside of plot to finish', ...
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
                
                % Display the prior ens mean and sd
                prior_ens_mean = mean(x);
                prior_ens_sd   = std(x);
                str1 = sprintf('Prior Ens Mean = %.4f', prior_ens_mean);
                set(handles.ui_text_prior_ens_mean, 'String', str1, 'Visible', 'on');
                str1 = sprintf('Prior Ens SD = %.4f', prior_ens_sd);
                set(handles.ui_text_prior_ens_sd,   'String', str1, 'Visible', 'on');
                
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
        
        % Enable the cycle DA button
        set(handles.ui_button_cycle_da,     'Enable', 'On');
        set(handles.ui_edit_growth_rate,      'Enable', 'On');
        set(handles.ui_edit_obs_error_sd,     'Enable', 'On');
        
    end


%% -----------------------------------------------------------------------------

    function button_cycle_da_Callback (~, ~)

        % Turn Off any old points
        set(handles.h_ens_member,     'Visible', 'Off');
        set(handles.h_update_ens,     'Visible', 'Off');
        set(handles.h_prior_pdf,      'Visible', 'Off');
        set(handles.h_post_pdf,       'Visible', 'Off');
        set(handles.h_prior_cont_pdf, 'Visible', 'Off');
        set(handles.h_post_cont_pdf,  'Visible', 'Off');
        
        % Remove mean and sd of old posterior
        clear_ui_labels;
        
        % And the lines in between
        set(handles.h_update_lines, 'Visible', 'Off');

        % The true value is zero in this example; obs error sd comes from the box
        % Generate a random observation
        handles.observation = randn * handles.obs_error_sd;
        
        % Move the observation asterisk
        set(handles.h_obs_ast, 'Visible', 'Off');
        handles.h_obs_ast = plot(handles.observation, 0, 'r*', 'MarkerSize', 16,'LineWidth',2.0);
        hold on;

        % Do the assimilation for the continuous normal distribution
        post_cont_var = 1 / (1 / handles.prior_cont_sd^2 + 1 / handles.obs_error_sd^2);
        post_cont_sd = sqrt(post_cont_var);
        post_cont_mean = post_cont_var * (handles.prior_cont_mean / handles.prior_cont_sd^2 + ...
           handles.observation / handles.obs_error_sd^2);

        % Plot the prior and posterior continuous normal (KF) PDFs
        handles.h_prior_cont_pdf = plot_gaussian(handles.prior_cont_mean, handles.prior_cont_sd, 1);
        set(handles.h_prior_cont_pdf, 'Color', atts.green, 'Linestyle', '-', 'Linewidth', 1.7, 'Visible', 'on');
        handles.h_post_cont_pdf  = plot_gaussian(post_cont_mean, post_cont_sd, 1);
        set(handles.h_post_cont_pdf, 'Color', atts.blue, 'Linestyle', '-', 'Linewidth', 1.7);
        
        % Plot the updated distribution
        set(handles.h_obs_plot, 'Visible', 'Off');
        handles.h_obs_plot = plot_gaussian(handles.observation, handles.obs_error_sd, 1);
        set(handles.h_obs_plot, 'Color', atts.red, 'Linestyle', '-', 'Linewidth', 1.7);

        % Display the prior cont mean and sd
        str1 = sprintf('Prior Cont Mean = %.4f', handles.prior_cont_mean);
        set(handles.ui_text_prior_cont_mean, 'String', str1, 'Visible', 'on');
        str1 = sprintf('Prior Cont SD = %.4f', handles.prior_cont_sd);
        set(handles.ui_text_prior_cont_sd,   'String', str1, 'Visible', 'on');
       
        % Update mean and sd of old cont posterior
        str1 = sprintf('Posterior Cont Mean = %.4f', post_cont_mean);
        set(handles.ui_text_post_cont_mean, 'String', str1, 'Visible', 'on');
        str1 = sprintf('Posterior Cont SD = %.4f', post_cont_sd);
        set(handles.ui_text_post_cont_sd,   'String', str1, 'Visible', 'on');

        % Update the continuous for the next cycle
        handles.prior_cont_mean = post_cont_mean * handles.growth_rate;
        handles.prior_cont_sd = post_cont_sd * handles.growth_rate;
        
        % If ensemble has been created, update an plot it also
        ens_exists = handles.ens_size > 1;
        if(ens_exists) 
           % Get the prior ensmeble members 
           ensemble = handles.ens_members;
        
           % Figure out which filter option is currently selected
           val = get(handles.ui_radio_button_group,'SelectedObject');
           filter_type = get(val,'String');
        
           switch filter_type
            
               case 'EAKF'
                   [obs_increments, ~] = ...
                       obs_increment_eakf(ensemble, handles.observation, handles.obs_error_sd^2);

                   handles.h_prior_pdf = plot_gaussian(mean(ensemble), std(ensemble), 1);
                   set(handles.h_prior_pdf, 'linewidth', 2, 'color', atts.green, 'Linestyle', '--');
               case 'EnKF'
                   [obs_increments, ~] = ...
                       obs_increment_enkf(ensemble, handles.observation, handles.obs_error_sd^2);
                   % There is no prior distribution to plot for the EnKF, it's not a QCEF
                   % However, here we can plot the fit to the ensemble anyway to show how it compares
                   handles.h_prior_pdf = plot_gaussian(mean(ensemble), std(ensemble), 1);
                   set(handles.h_prior_pdf, 'linewidth', 2, 'color', atts.green, 'Linestyle', '--');
               case 'RHF'
                   bounded_left = false;
                   [obs_increments, err, xp_prior, yp_prior, xp_post, yp_post] = ...
                       obs_increment_rhf(ensemble, handles.observation, handles.obs_error_sd^2, bounded_left);
                   handles.h_prior_pdf = plot(xp_prior, yp_prior, 'linewidth', 2, 'color', atts.green, 'Linestyle', '--');
                   handles.h_post_pdf = plot(xp_post, yp_post, 'linewidth', 2, 'color', atts.blue, 'Linestyle', '--');
           end
        
           % Add on increments to get new ensemble
           new_ensemble = ensemble + obs_increments;
       
           % Plot the prior and posterior ensemble members 
           ens_size = size(ensemble, 2);
           for i = 1:ens_size
              handles.h_ens_member(i) = ...
                       plot(handles.ens_members(i), 0, '*', 'MarkerSize', 16, 'Color', atts.green);
           end
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
                
           % Display the prior ens mean and sd
           prior_ens_mean = mean(handles.ens_members);
           prior_ens_sd   = std(handles.ens_members);
           str1 = sprintf('Prior Ens Mean = %.4f', prior_ens_mean);
           set(handles.ui_text_prior_ens_mean, 'String', str1, 'Visible', 'on');
           str1 = sprintf('Prior Ens SD = %.4f', prior_ens_sd);
           set(handles.ui_text_prior_ens_sd,   'String', str1, 'Visible', 'on');
        
           % Update mean and sd of old ens posterior
           str1 = sprintf('Posterior Ens Mean = %.4f',new_mean);
           set(handles.ui_text_post_ens_mean, 'String', str1, 'Visible', 'on');
           str1 = sprintf('Posterior Ens SD = %.4f',new_sd);
           set(handles.ui_text_post_ens_sd,   'String', str1, 'Visible', 'on');
       
           % Plot the posterior sample continuous for the EAKF or ENKF
           if(strcmp(filter_type, 'EAKF') | strcmp(filter_type, 'EnKF'))
                handles.h_post_pdf = plot_gaussian(new_mean, new_sd, 1);
                set(handles.h_post_pdf, 'linewidth', 2, 'color', atts.blue, 'Linestyle', '--');
           end

           % Update the ensemble for next cycle
           handles.ens_members = new_ensemble * handles.growth_rate;
        end

    end

%% -----------------------------------------------------------------------------

    function clear_ui_labels()
        
        % Turns Off all labels except for the prior ens mean and SD
        set(handles.ui_text_post_ens_sd,             'Visible', 'Off');
        set(handles.ui_text_post_ens_mean,           'Visible', 'Off');
        
        % Turns Off all labels except for the prior ens mean and SD
        set(handles.ui_text_post_cont_sd,             'Visible', 'Off');
        set(handles.ui_text_post_cont_mean,           'Visible', 'Off');
        
    end

%% -----------------------------------------------------------------------------


    function edit_growth_rate_Callback(~, ~)
        
        % Turn Off any old updated points
        set(handles.h_update_ens,     'Visible', 'Off');
        set(handles.h_prior_pdf,      'Visible', 'Off');
        set(handles.h_post_pdf,       'Visible', 'Off');
        set(handles.h_prior_cont_pdf, 'Visible', 'Off');
        set(handles.h_post_cont_pdf,  'Visible', 'Off');
        
        % Remove mean and sd of old posterior
        clear_ui_labels;
        
        % And the lines in between
        set(handles.h_update_lines,        'Visible', 'Off');
        
        % Enable things that an error might have turned Off
        set(handles.ui_edit_obs_error_sd,     'Enable', 'on')
        set(handles.ui_button_create_new_ens, 'Enable', 'on')
        
        % Only enable the cycle DA pushbutton if an ensemble has been created
        if(handles.ens_size > 0)
            set(handles.ui_button_cycle_da, 'Enable', 'on');
            
        end
        
        % Get the value of the growth_rate
        if(isfinite(      str2double(get(handles.ui_edit_growth_rate, 'String'))))
            growth_rate = str2double(get(handles.ui_edit_growth_rate, 'String'));
            
        else
            set(handles.ui_edit_growth_rate, 'String', '??','FontWeight','Bold', ...
                'BackgroundColor', atts.red);
            set(handles.ui_text_growth_rate_err_print,'Visible','On')
            
            fprintf('ERROR: Growth Rate value must be numeric.\n')
            fprintf('ERROR: Growth Rate value must be numeric.\n')
            
            
            % Disable other input to guarantee only one error at a time!
            set(handles.ui_edit_obs_error_sd,     'Enable', 'Off')
            set(handles.ui_button_create_new_ens, 'Enable', 'Off')
            set(handles.ui_button_cycle_da,     'Enable', 'Off')
            
            return
            
        end
        
        % Update the global storage
        handles.growth_rate = growth_rate;
        set(handles.ui_edit_growth_rate, 'BackgroundColor', 'White','FontWeight', 'Normal');
        set(handles.ui_text_growth_rate_err_print,'Visible','Off')
        
    end

%% -----------------------------------------------------------------------------

    function edit_obs_error_sd_Callback(~, ~)
        
        % Turn Off any old updated points
        set(handles.h_update_ens,          'Visible', 'Off');
        set(handles.h_prior_pdf,           'Visible', 'Off');
        set(handles.h_post_pdf,            'Visible', 'Off');
        set(handles.h_prior_cont_pdf,      'Visible', 'Off');
        set(handles.h_post_cont_pdf,       'Visible', 'Off');
        
        % Remove mean and sd of old posterior
        clear_ui_labels;
        
        % And the lines in between
        set(handles.h_update_lines,        'Visible', 'Off');
        
        % Enable things that an error might have turned Off
        set(handles.ui_edit_growth_rate,      'Enable', 'on')
        set(handles.ui_button_create_new_ens, 'Enable', 'on')
        
        % Only enable the update ensemble pushbutton if an ensemble has been created
        if(handles.ens_size > 0)
            set(handles.ui_button_cycle_da, 'Enable', 'on');
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
            set(handles.ui_edit_growth_rate,      'Enable', 'Off')
            set(handles.ui_button_create_new_ens, 'Enable', 'Off')
            set(handles.ui_button_cycle_da,     'Enable', 'Off')
            
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
