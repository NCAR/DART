function run_lorenz_96_inf
%% RUN_LORENZ_96_INF ensemble data assimilation with a 40-variable implementation of the Lorenz '96 dynamical model.
%
%      To demonstrate the analogue to the atmosphere, the model is a cyclic
%      1D domain with equally-spaced nodal points. There are 20 ensemble
%      members initially in this example.
%
%      The model can be single-stepped through model advance and assimilation
%      steps using the top pushbutton, or allowed to run Auto using the
%      'Start Auto Run' button. A variety of assimilation algorithms can
%      be selected from the first pulldown. Model error in the assimilating
%      model (an imperfect model assimilation experiment) can be selected
%      with the second pulldown. The localization, (prior) adaptive-inflation and ensemble
%      size can be changed with the three dialogue boxes. Different observation
%      networks can be selected to test he performance. 
%
%      For prior adaptive inflation algorithm; there are 2 choices: 
%      1- Inflation is assumed to follow a Gaussian pdf (from Anderson 2009)
%      2- Inflation is assumed to follow an Inverse Gamma pdf (from El Gharamti 2018)
%
%      Changing the ensemble size resets the diagnostic displays. The figure window
%      displays time sequences of the prior and posterior error and prior
%      and posterior (if assimilation is on) rank histograms.
%
%      It takes about 20 timesteps for the intially small ensemble
%      perturbations to grow large enough to be seen using the default
%      settings.
%
% See also: gaussian_product.m oned_model.m oned_model_inf.m oned_ensemble.m
%           twod_ensemble.m run_lorenz_63.m run_lorenz_96.m

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

help run_lorenz_96_inf

% Initialize global storage for control trajectory and ensembles
% Ordinarily, each MATLAB function has its own local variables, which are
% separate from those of other functions and from those of the base workspace.
% However, if several functions all declare a particular variable name as global,
% then they all share a single copy of that variable. Any change of value to
% that variable, in any function, is visible to all the functions that declare
% it as global.

global TRUE_FORCING
global FORCING
global MODEL_SIZE
global DELTA_T
global MEAN_DIST

LOG_FILE = strcat(mfilename, '.log');

first_call_to_reset = true;

atts = stylesheet; % get the default fonts and colors

figWidth  = 1024;    % in pixels
figHeight = 990;    % in pixels

%% Create figure Layout
figure('position', [400 100 figWidth figHeight], ...
    'Units', 'pixels', ...
    'Name', 'run_lorenz_96_inf', ...
    'Color', atts.background);

%% Create text in the top right corner with the elapsed model time
handles.time = 1;
handles.ui_text_time = uicontrol('Style', 'text', ...
    'Units' , 'Normalized', ...
    'Position', [0.48 0.48 0.20 0.07], ...
    'BackgroundColor', atts.background, ...
    'HorizontalAlignment', 'left', ...
    'String', sprintf('Time = %d',handles.time), ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.45);

%% Create an [Advance Model/Assimilate Obs] button that lets you do a single_step
% The 'Callback' forces the SingleStep function to be called when clicked.
handles.ui_button_Single_Step = uicontrol('Style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Position', [0.02 0.900 0.24 0.083], ...
    'BackgroundColor', 'White', ...
    'String', 'Advance Model', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontWeight', 'bold', ...
    'FontSize', 0.37, ...
    'Callback', @SingleStep_Callback);

%% Create a [Start/Pause Auto Run] button that lets you do a Auto run
handles.ui_button_Auto_Run = uicontrol('Style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Position', [0.02 0.800 0.24 0.083], ...
    'BackgroundColor', 'White', ...
    'String', 'Start Auto Run', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontWeight', 'bold', ...
    'FontSize', 0.37, ...
    'Callback', @AutoRun_Callback);

%% Create a label for the Radio Group
handles.ui_text_group_label = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.287 0.952 0.283 0.033], ...
    'BackgroundColor', atts.background, ...
    'String', 'Assimilation Type:', ...
    'HorizontalAlignment', 'left', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontWeight', 'bold', ...
    'FontSize', 0.8);

%% Create a Radio Button with choices for the types of assimilation

handles.ui_button_group_assimilation = uibuttongroup('Visible','off', ...
    'Position',[0.267 0.800 0.211 0.150], ...
    'BorderType', 'none', ...
    'BackgroundColor' , atts.background, ...
    'SelectionChangeFcn', @Assimilation_selection);

% Create the four radio buttons in the button group.
handles.ui_radio_noAssimilation = uicontrol(handles.ui_button_group_assimilation, ...
    'Style','radiobutton', ...
    'Units' , 'Normalized', ...
    'Position',[0.105 0.767 0.842 0.222], ...
    'BackgroundColor' , atts.background, ...
    'String','No Assimilation', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'Normalized' , ...
    'FontSize' , 0.6, ...
    'HandleVisibility','off');

handles.ui_radio_EAKF = uicontrol(handles.ui_button_group_assimilation, ...
    'Style','radiobutton', ...
    'Units' , 'Normalized', ...
    'Position',[0.105 0.522 0.842 0.222], ...
    'BackgroundColor' , atts.background, ...
    'String', 'EAKF', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'Normalized' , ...
    'FontSize' , 0.6, ...
    'HandleVisibility','off');

handles.ui_radio_EnKF = uicontrol(handles.ui_button_group_assimilation, ...
    'Style', 'radiobutton', ...
    'Units' , 'Normalized', ...
    'Position',[0.105 0.278 0.842 0.222], ...
    'BackgroundColor' , atts.background, ...
    'String', 'EnKF', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'Normalized' , ...
    'FontSize' , 0.6, ...
    'HandleVisibility','off');

handles.ui_radio_RHF = uicontrol(handles.ui_button_group_assimilation, ...
    'Style', 'radiobutton', ...
    'Units' , 'Normalized', ...
    'Position',[0.105 0.033 0.842 0.222], ...
    'BackgroundColor' , atts.background, ...
    'String', 'RHF', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'Normalized' , ...
    'FontSize' , 0.6, ...
    'HandleVisibility','off');

% Make the uibuttongroup visible after creating child objects.
set(handles.ui_button_group_assimilation, 'Visible','On');

handles.filter_type_string = 'No Assimilation';

%% Create a label for the obs Radio Group
handles.ui_text_group_label = uicontrol('Style', 'text', ...
    'Units'                 , 'Normalized', ...
    'Position'              , [0.69 0.64 0.30 0.033], ...
    'BackgroundColor'       , atts.background, ...
    'String'                , 'Observation Network', ...
    'HorizontalAlignment'   , 'left', ...
    'FontName'              , atts.fontname, ...
    'FontWeight'            , 'bold', ...
    'FontUnits'             , 'normalized', ...
    'FontSize'              , 0.8);

%% Create a Radio Button with choices for the observation network
handles.ui_button_group_obsnetwork = uibuttongroup('Visible','off', ...
    'Position'          ,[0.63 0.55 0.35 0.08], ...
    'BorderType', 'none', ...
    'BackgroundColor'   , atts.background, ...
    'FontWeight'        ,'bold', ...
    'SelectionChangeFcn', @ObsNetwork_selection);

% Create the six radio buttons in the button group.
handles.ui_radio_full = uicontrol(handles.ui_button_group_obsnetwork, ...
    'Style'             , 'radiobutton', ...
    'Units'             , 'Normalized', ...
    'Position'          , [0.12 0.6 0.4 0.4], ...
    'BackgroundColor'   , atts.background, ...
    'String'            , '1:40:1', ...
    'FontName'          , atts.fontname, ...
    'FontUnits'         , 'Normalized' , ...
    'FontSize'          , 0.75, ...
    'HandleVisibility'  ,'off');

handles.ui_radio_half = uicontrol(handles.ui_button_group_obsnetwork, ...
    'Style'             , 'radiobutton', ...
    'Units'             , 'Normalized', ...
    'Position'          , [0.38 0.6 0.4 0.4], ...
    'BackgroundColor'   , atts.background, ...
    'String'            , '1:40:2', ...
    'FontName'          , atts.fontname, ...
    'FontUnits'         , 'Normalized' , ...
    'FontSize'          , 0.75, ...
    'HandleVisibility'  ,'off');

handles.ui_radio_quarter = uicontrol(handles.ui_button_group_obsnetwork, ...
    'Style'             , 'radiobutton', ...
    'Units'             , 'Normalized', ...
    'Position'          , [0.67 0.6 0.4 0.4], ...
    'BackgroundColor'   , atts.background, ...
    'String'            , '1:40:4', ...
    'FontName'          , atts.fontname, ...
    'FontUnits'         , 'Normalized' , ...
    'FontSize'          , 0.75, ...
    'HandleVisibility'  ,'off');

handles.ui_radio_first20 = uicontrol(handles.ui_button_group_obsnetwork, ...
    'Style'             , 'radiobutton', ...
    'Units'             , 'Normalized', ...
    'Position'          , [0.12 0.1 0.4 0.4], ...
    'BackgroundColor'   , atts.background, ...
    'String'            , '1:20', ...
    'FontName'          , atts.fontname, ...
    'FontUnits'         , 'Normalized' , ...
    'FontSize'          , 0.75, ...
    'HandleVisibility'  ,'off');

handles.ui_radio_mid20 = uicontrol(handles.ui_button_group_obsnetwork, ...
    'Style'             , 'radiobutton', ...
    'Units'             , 'Normalized', ...
    'Position'          , [0.38 0.1 0.4 0.4], ...
    'BackgroundColor'   , atts.background, ...
    'String'            , '10:30', ...
    'FontName'          , atts.fontname, ...
    'FontUnits'         , 'Normalized' , ...
    'FontSize'          , 0.75, ...
    'HandleVisibility'  ,'off');

handles.ui_radio_firstlast10 = uicontrol(handles.ui_button_group_obsnetwork, ...
    'Style'             , 'radiobutton', ...
    'Units'             , 'Normalized', ...
    'Position'          , [0.67 0.1 0.4 0.4], ...
    'BackgroundColor'   , atts.background, ...
    'String'            , '1:10; 30:40', ...
    'FontName'          , atts.fontname, ...
    'FontUnits'         , 'Normalized' , ...
    'FontSize'          , 0.5, ...
    'HandleVisibility'  ,'off');

% Make the uibuttongroup visible after creating child objects.
set(handles.ui_button_group_obsnetwork, 'Visible','On');

handles.obs_network_string = '1:40:1';

%% Create a Text Box in above the slider that defines the true model state value
handles.ui_text_slider_caption = uicontrol('style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.430 0.84 0.211 0.067], ...
    'BackgroundColor', atts.background, ...
    'String', 'True State has F = 8', ...
    'HorizontalAlignment', 'center', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.35);

%% Create a Slider for Model Forcing/Error
% Sliderstep: 8 values from 4-12.
%    1/8 forces slider to move 1 value (1/8 of 8) when arrow is clicked.
%    1/2 forces slider to move 4 values when the actual slider is pressed
handles.ui_slider_error = uicontrol('style', 'slider', ...
    'Units', 'Normalized', ...
    'Position', [0.43 0.750 0.222 0.062], ...
    'BackgroundColor', [0.2 0.2 0.2], ...
    'ForegroundColor', atts.background, ...
    'Min', 4, ...
    'Max' ,12, ...
    'Value', 8 , ...
    'Sliderstep',[1/8 1/2], ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.4, ...
    'Callback', @Forcing_Callback);

%% Create an edit field if the user wants a more specific forcing value
handles.ui_text_forcing = uicontrol('style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.43 0.817 0.144 0.045], ...
    'BackgroundColor', atts.background, ...
    'String', 'Forcing:', ...
    'HorizontalAlignment', 'center', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.5);

handles.ui_edit_forcing = uicontrol('style', 'edit', ...
    'Units', 'Normalized', ....
    'Position', [0.56 0.82 0.060 0.043], ...
    'String', '8', ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.6, ...
    'Callback', @edit_forcing_Callback);

%% Create a panel consisting of two boxes and two edit boxes next to them.
% Only the edit boxes have a callback.

handles.RedPanel = uipanel('Units','Normalized', ...
    'Position',[0.68 0.88 0.278 0.104], ...
    'BorderType', 'none', ...
    'BackgroundColor', atts.background);

handles.ui_text_localization = uicontrol(handles.RedPanel, ...
    'Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.022 0.583 0.427 0.291], ...
    'String', 'Localization', ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', atts.background, ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontWeight','normal', ...
    'FontSize', 0.7);

handles.ui_edit_localization = uicontrol(handles.RedPanel, ...
    'Style', 'edit', ...
    'Units', 'Normalized', ...
    'Position', [0.055 0.013 0.361 0.539], ...
    'String', '1.0', ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.7, ...
    'Callback', @edit_localization_Callback);

handles.ui_text_ens_size = uicontrol(handles.RedPanel, ...
    'Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.618 0.583 0.382 0.291], ...
    'String', 'Ens. Size', ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', atts.background, ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontWeight','normal', ...
    'FontSize', 0.8);

handles.ui_edit_ens_size = uicontrol(handles.RedPanel, ...
    'Style', 'edit', ...
    'Units', 'Normalized', ....
    'Position', [0.618 0.013 0.361 0.539], ...
    'String', '20', ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.7, ...
    'Callback', @edit_ens_size_Callback);


%% Create a panel consisting of three boxes and three edit boxes next to them.
% Only the edit boxes have a callback.

handles.InfPanel = uipanel('Units','Normalized', ...
    'Position',[0.667 0.68 0.311 0.17], ...
    'BorderType', 'none', ...
    'BackgroundColor', atts.background);

handles.ui_text_inflation = uicontrol(handles.InfPanel, ...
    'Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.1 0.7 0.800 0.30], ...
    'String', 'Adaptive Inflation', ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', atts.background, ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontWeight','bold', ...
    'FontSize', 0.52);

% Select inflation algorithm
handles.ui_button_group_infalgo = uibuttongroup('Visible','off', ...
    'Position'          ,[0.67 0.765 0.35 0.05], ...
    'BorderType', 'none', ...
    'BackgroundColor'   , atts.background, ...
    'FontWeight'        ,'bold', ...
    'SelectionChangeFcn', @InfAlgo_selection);

% Create the two radio buttons in the button group.
handles.ui_radio_gauss = uicontrol(handles.ui_button_group_infalgo, ...
    'Style'             , 'radiobutton', ...
    'Units'             , 'Normalized', ...
    'Position'          , [0.10 0.3 0.4 0.4], ...
    'BackgroundColor'   , atts.background, ...
    'String'            , 'Gaussian', ...
    'FontName'          , atts.fontname, ...
    'FontUnits'         , 'Normalized' , ...
    'FontSize'          , 1, ...
    'HandleVisibility'  ,'off');

handles.ui_radio_gamma = uicontrol(handles.ui_button_group_infalgo, ...
    'Style'             , 'radiobutton', ...
    'Units'             , 'Normalized', ...
    'Position'          , [0.46 0.3 0.4 0.4], ...
    'BackgroundColor'   , atts.background, ...
    'String'            , 'I-Gamma', ...
    'FontName'          , atts.fontname, ...
    'FontUnits'         , 'Normalized' , ...
    'FontSize'          , 1, ...
    'HandleVisibility'  ,'off');

% Make the uibuttongroup visible after creating child objects.
set(handles.ui_button_group_infalgo, 'Visible','On');

handles.inf_algorithm_string = 'Gaussian';

handles.ui_text_inflation_Std = uicontrol(handles.InfPanel, ...
    'Style', 'text', ...
    'Units', 'Normalized', ....
    'Position', [0.05 0.27 0.200 0.250], ...
    'String', 'Inf. Std', ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', atts.background, ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontWeight','normal', ...
    'FontSize', 0.45);

handles.ui_edit_inflation_Std = uicontrol(handles.InfPanel, ...
    'Style', 'edit', ...
    'Units', 'Normalized', ....
    'Position', [0.05 0.07 0.210 0.270], ...
    'String', '0.3', ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.6, ...
    'Callback', @edit_inflation_Std_Callback);

handles.ui_text_inflation_Damp = uicontrol(handles.InfPanel, ...
    'Style', 'text', ...
    'Units', 'Normalized', ....
    'Position', [0.36 0.27 0.300 0.250], ...
    'String', 'Inf. Damp', ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', atts.background, ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontWeight','normal', ...
    'FontSize', 0.45);

handles.ui_edit_inflation_Damp = uicontrol(handles.InfPanel, ...
    'Style', 'edit', ...
    'Units', 'Normalized', ....
    'Position', [0.40 0.07 0.210 0.270], ...
    'String', '1.0', ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.6, ...
    'Callback', @edit_inflation_Damp_Callback);

handles.ui_text_inflation_Min = uicontrol(handles.InfPanel, ...
    'Style', 'text', ...
    'Units', 'Normalized', ....
    'Position', [0.70 0.27 0.300 0.250], ...
    'String', 'Inf. Min', ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', atts.background, ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontWeight','normal', ...
    'FontSize', 0.45);

handles.ui_edit_inflation_Min = uicontrol(handles.InfPanel, ...
    'Style', 'edit', ...
    'Units', 'Normalized', ....
    'Position', [0.75 0.07 0.210 0.270], ...
    'String', '0.0', ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.6, ...
    'Callback', @edit_inflation_Min_Callback);

%% Reset button - clear the whole thing

handles.ResetButton = uicontrol('Style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Position', [0.74 0.02 0.18 0.070], ...
    'String', 'Reset', ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.6, ...
    'FontWeight', 'bold', ...
    'Callback', @reset_Callback);

%% Button to reset the histograms

handles.ClearStats = uicontrol('Style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Position', [0.50 0.02 0.18 0.070], ...
    'String', 'Clear Statistics', ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'normalized', ...
    'FontSize', 0.3, ...
    'FontWeight', 'bold', ...
    'Callback', @ClearStats_Callback);

%% -----------------------------------------------------------------------------
%  These appear to be error messages that can be turned on or off.
%  They start turned off.

handles.ui_text_ens_size_err_print = uicontrol('Style', 'text', ...
    'Units'             , 'Normalized', ...
    'Position'          , [0.1 0.6000 0.3000 0.0800], ...
    'String'            , 'ERROR: Ens. Size value must be greater than or equal to 2.', ...
    'BackgroundColor'   , 'White', ...
    'ForegroundColor'   , atts.red, ...
    'FontName'          , atts.fontname, ...
    'FontSize'          , atts.fontsize, ...
    'FontWeight'        , 'Bold', ...
    'Visible'           ,'Off');

handles.ui_text_localization_err_print = uicontrol('Style', 'text', ...
    'Units'             , 'Normalized', ...
    'Position'          , [0.1 0.6000 0.3000 0.0800], ...
    'String'            , 'ERROR: localization must be greater than 0.', ...
    'BackgroundColor'   , 'White', ...
    'ForegroundColor'   , atts.red, ...
    'FontName'          , atts.fontname, ...
    'FontSize'          , atts.fontsize, ...
    'FontWeight'        , 'Bold', ...
    'Visible'           ,'Off');

handles.ui_text_inf_std_err_print = uicontrol('Style', 'text', ...
    'Units'             , 'Normalized', ...
    'Position'          , [0.1 0.6000 0.3000 0.0800], ...
    'String'            , 'ERROR: Inf. Std value must be greater than 0.', ...
    'BackgroundColor'   , 'White', ...
    'ForegroundColor'   , atts.red, ...
    'FontName'          , atts.fontname, ...
    'FontSize'          , atts.fontsize, ...
    'FontWeight'        , 'Bold', ...
    'Visible'           ,'Off');

handles.ui_text_inf_damp_err_print = uicontrol('Style', 'text', ...
    'Units'             , 'Normalized', ...
    'Position'          , [0.1 0.6000 0.3000 0.0800], ...
    'String'            , 'ERROR: Inf. Damp value must be between 0.1 and 1.', ...
    'BackgroundColor'   , 'White', ...
    'ForegroundColor'   , atts.red, ...
    'FontName'          , atts.fontname, ...
    'FontSize'          , atts.fontsize, ...
    'FontWeight'        , 'Bold', ...
    'Visible'           ,'Off');

handles.ui_text_inf_damp_warn_print = uicontrol('Style', 'text', ...
    'Units'             , 'Normalized', ...
    'Position'          , [0.07 0.755 0.400 0.0300], ...
    'String'            , 'Warning: This is a no inflation scenario!', ...
    'BackgroundColor'   , atts.background, ...
    'ForegroundColor'   , atts.orange, ...
    'FontName'          , atts.fontname, ...
    'FontSize'          , 15, ...
    'FontWeight'        , 'Bold', ...
    'Visible'           ,'Off');

handles.ui_text_inf_min_err_print = uicontrol('Style', 'text', ...
    'Units'             , 'Normalized', ...
    'Position'          , [0.1 0.6000 0.3000 0.0800], ...
    'String'            , 'ERROR: Inf. Min value must be greater or equal to 0.', ...
    'BackgroundColor'   , 'White', ...
    'ForegroundColor'   , atts.red, ...
    'FontName'          , atts.fontname, ...
    'FontSize'          , atts.fontsize, ...
    'FontWeight'        , 'Bold', ...
    'Visible'           ,'Off');

%% Set the positions to the plot components. These handles must exist before
% the call to reset_Callback()

positions = [ ...
    50/figWidth 560/figHeight 380/figWidth 150/figHeight; ...
    50/figWidth 340/figHeight 380/figWidth 150/figHeight; ...
    50/figWidth  80/figHeight 185/figWidth 170/figHeight; ...
    245/figWidth  80/figHeight 185/figWidth 170/figHeight; ...
    510/figWidth 100/figHeight 400/figWidth 400/figHeight];

handles.timeseries           = subplot('Position', positions(1,:));
handles.infseries            = subplot('Position', positions(2,:));
handles.prior_rank_histogram = subplot('Position', positions(3,:));
handles.post_rank_histogram  = subplot('Position', positions(4,:));
handles.polar_plot           = subplot('Position', positions(5,:));

% Set the initial values for all the components.
% reset_Callback() sets the default plotting arrays, so it has to be called
% before the rest of the setup.

reset_Callback();

% Plot default values for rmse and spread so we can set the line types etc
% for the legend. After that, make the default objects invisible.
% Really just setting the stuff in the legend.

subplot(handles.timeseries);
set(handles.timeseries, 'FontSize', atts.fontsize);

h_prior_rms        = plot(handles.prior_rms       , '-' , 'LineWidth',2.0, ...
    'Color', atts.green); hold on
h_posterior_rms    = plot(handles.posterior_rms   , '-' , 'LineWidth',2.0, ...
    'Color', atts.blue );
h_prior_spread     = plot(handles.prior_spread    , '-.', 'LineWidth',2.0, ...
    'Color', atts.green);
h_posterior_spread = plot(handles.posterior_spread, '-.', 'LineWidth',2.0, ...
    'Color', atts.blue );

% Sadly, these dont seem to scale - even when normalized.
h = legend('Prior RMSE', 'Posterior RMSE', 'Prior Spread', 'Posterior Spread');
set(h, 'FontSize', atts.fontsize, 'Position',[0.46 0.62 0.118 0.148], 'EdgeColor', 'w');

if verLessThan('matlab','R2017a')
    % Convince Matlab to not autoupdate the legend with each new line.
    % Before 2017a, this was the default behavior, so do nothing.
    % We do not want to add the bias line to the legend, for example.
else
    h.AutoUpdate = 'off';
end

ylabel('RMSE & Spread', 'FontSize', atts.fontsize);
xlabel('Time',          'FontSize', atts.fontsize);

rms_time = 1;

set(h_prior_rms,        'Visible', 'on')
set(h_posterior_rms,    'Visible', 'on')
set(h_prior_spread,     'Visible', 'on')
set(h_posterior_spread, 'Visible', 'on')

subplot(handles.infseries);
set(handles.infseries, 'FontSize', atts.fontsize);

h_inflation = plot(handles.inflate, '-x', 'Color', atts.blue, ...
    'MarkerFaceColor', atts.red);
set(gca, 'XLim', [1, MODEL_SIZE], 'XTick', 5:5:MODEL_SIZE ); box on

ylabel('Deflation | Inflation', 'FontSize', atts.fontsize);
xlabel('Location', 'FontSize', atts.fontsize);

set(h_inflation, 'Visible', 'off')

subplot(handles.prior_rank_histogram);
set(handles.prior_rank_histogram, 'FontSize', atts.fontsize);
ylabel('Frequency');
xlabel('Rank');
title ('Prior Rank Histogram', 'FontSize', atts.fontsize);
hold on

subplot(handles.post_rank_histogram);
set(handles.post_rank_histogram, 'FontSize', atts.fontsize, 'YAxisLocation','right')
ylabel('Frequency');
xlabel('Rank');
title ('Posterior Rank Histogram', 'FontSize', atts.fontsize);
hold on

subplot(handles.polar_plot);
polar_y = (0:MODEL_SIZE) / MODEL_SIZE * 2 * pi;

% Plot a single point to make sure the axis limits are okay
% Unclear if these can be set more cleanly with polar
% This also gets the observations into the legend

h_obs   = plot_polar([0, 2*pi], 14.9, MEAN_DIST, 'r*', 1); hold on
handles.h_ens   = plot_polar([0, 2*pi], 1000, MEAN_DIST, '-g', 1);
handles.h_truth = plot_polar([0, 2*pi], 1000, MEAN_DIST, '-k', 1);

set(handles.h_truth, 'linewidth', 3);
set(handles.h_ens  , 'linewidth', 1, 'Color', atts.green);
set(        h_obs  , 'linewidth', 1, 'Color', atts.red, 'Visible', 'off');

L = legend( [handles.h_truth, handles.h_ens, h_obs], 'True State', ...
    'Ensemble', 'Observations', 'Location', 'NorthEast');
pos   = get(L, 'Position') + [0.09 0.01 0.021 0.012];
set(L, 'Position', pos, 'FontSize', atts.fontsize, 'EdgeColor', 'w');

if verLessThan('matlab','R2017a')
    % Convince Matlab to not autoupdate the legend with each new line.
    % Before 2017a, this was the default behavior, so do nothing.
    % We do not want to add the bias line to the legend, for example.
else
    L.AutoUpdate = 'off';
end

% Plot the localization width
h_loc = plot_localization;


%% -----------------------------------------------------------------------------
%  Initiate log file
if exist(LOG_FILE, 'file') == 2
    logfileid = fopen(LOG_FILE, 'a');
else
    logfileid = fopen(LOG_FILE, 'w');
    fprintf(logfileid, '-----------------------------------------------------------\n');
    fprintf(logfileid, '------------------------ DART_LAB -------------------------\n');
    fprintf(logfileid, '-----------------------------------------------------------\n\n');
    fprintf(logfileid, '********************* %s **********************\n\n', mfilename);
end

fprintf(logfileid, '\n\nNEW RUN: Starting date and time %s\n', datestr(datetime));
fprintf(logfileid, '========\n\n');

fprintf(logfileid, '# Time step: %.2f (Initial configuration)\n', DELTA_T);
fprintf(logfileid, '  - Forcing `F` parameter = %.2f\n', TRUE_FORCING);
fprintf(logfileid, '  - Assimilation type is `%s`\n', handles.filter_type_string);
fprintf(logfileid, '  - Observation Network is `%s`\n', handles.obs_network_string);
fprintf(logfileid, '  - Ensemble size = %d\n', handles.ens_size);
fprintf(logfileid, '  - Localization = %.2f\n', handles.localization);
fprintf(logfileid, '  - Inflation Algorithm is `%s`\n', handles.inf_algorithm_string);
fprintf(logfileid, '  - Adaptive inflation lower bound = %.2f\n', ...
    handles.inf_lower_bound);
fprintf(logfileid, '  - Adaptive inflation upper bound = %.2f\n', ...
    handles.inf_upper_bound);
fprintf(logfileid, '  - Adaptive inflation damping factor = %.2f\n', ...
    handles.inf_damping);
fprintf(logfileid, '  - Adaptive inflation standard deviation = %.2f\n', ...
    handles.inflate_sd);
fprintf(logfileid, '  - Adaptive inflation standard deviation lower bound = %.2f\n\n', ...
    handles.inf_sd_lower_bound);

fclose(logfileid);


%% ----------------------------------------------------------------------
%  All function below can use the variables defined above.
%% ----------------------------------------------------------------------

    function SingleStep_Callback(~, ~)
        % Called whenever [Advance Model/Assimilate Obs] is pressed.
        % If Assimilation is turned off ... we are always ready to advance.
        
        % Signal that something has happened.
        set(handles.ui_button_Single_Step, 'Enable', 'Off');
        
        if(strcmp(handles.filter_type_string, 'No Assimilation'))
            if (strcmp(get(handles.ui_button_Single_Step, 'String') , 'Assimilate Obs'))
                %If it says Assimilate Obs, that means another assimilation
                %was being used, so step ahead.
                step_ahead();
            end
            handles.ready_to_advance = true;
        end
        
        step_ahead();
        
        set(handles.ui_button_Single_Step, 'Enable', 'On');
    end

%% ----------------------------------------------------------------------

    function AutoRun_Callback(~, ~)
        % Specifies action for the 'Start/Pause Auto Run' button.
        
        % Turn off all the other model status controls to avoid a mess
        turn_off_controls();
        set(handles.ui_button_Auto_Run, 'Enable', 'On');
        
        % Check the button label to see if we are starting or stopping a Auto run
        if(strcmp(get(handles.ui_button_Auto_Run, 'String'), 'Pause Auto Run'))
            
            % Turn off the Auto run pushbutton until everything has completely stopped
            set(handles.ui_button_Auto_Run, 'Enable', 'Off');
            
            % Being told to stop; switch to not running status
            set(handles.ui_button_Auto_Run, 'String', 'Start Auto Run');
            
        else
            % Being told to start Auto run
            % Change the pushbutton to stop
            set(handles.ui_button_Auto_Run, 'String', 'Pause Auto Run');
            
            % Loop through advance and assimilate steps until stopped
            while(true)
                
                % Check to see if stop has been pushed
                status_string = get(handles.ui_button_Auto_Run, 'String');
                if(strcmp(status_string, 'Start Auto Run'))
                    turn_on_controls();
                    return
                end
                
                % Do the next advance or assimilation step
                step_ahead();
                drawnow
            end
        end
    end

%% ----------------------------------------------------------------------

    function Forcing_Callback(~, ~)
        %Called when the slider has been changed
        
        err         = get(handles.ui_slider_error, 'Value');
        old_forcing = FORCING;
        
        % Round the Value of the slider to the nearest integer
        % Set the Value of the slider to the rounded number. This will
        % Create a snap into place effect
        
        FORCING = round(err);
        set(handles.ui_slider_error, 'Value' , FORCING);
        set(handles.ui_edit_forcing, 'String' , sprintf('%d',FORCING));
        
        Update_log_file(rms_time, handles.time, handles.prior_rms, ...
            handles.prior_spread, 'Forcing (F)', old_forcing, FORCING);
    end

%% ----------------------------------------------------------------------

    function edit_inflation_Std_Callback(~, ~)
        % Is called when the edit_inflation field is changed
        
        inflate_sd_new = str2double(get(handles.ui_edit_inflation_Std, 'String'));
        old_inflate_sd = handles.inflate_sd;
        
        % Set the inflation value to the update
        handles.inflate_sd     = inflate_sd_new;
        handles.inf_sd_lower_bound = handles.inflate_sd;
        
        if(not(isfinite(handles.inflate_sd)) || handles.inflate_sd <= 0)
            % After this, only this edit box will work
            turn_off_controls();
            set(handles.ui_edit_inflation_Std, 'Enable', 'On');
            
            set(handles.ui_edit_inflation_Std, 'String', '?', ...
                'FontWeight','Bold','BackgroundColor', atts.red);
            set(handles.ui_text_inf_std_err_print, 'Visible','On')
            
            fprintf('ERROR: inflation std must be greater than 0.\n')
            fprintf('ERROR: unable to interpret inflation std value, please try again.\n')
            
            return
        end
        
        Update_log_file(rms_time, handles.time, handles.prior_rms, ...
            handles.prior_spread, 'Inflation Std', old_inflate_sd, handles.inflate_sd);
        
        % Enable all controls
        turn_on_controls();
        set(handles.ui_edit_inflation_Std, 'BackgroundColor', 'White', ...
            'FontWeight', 'Normal');
    end

%% ----------------------------------------------------------------------

    function edit_inflation_Damp_Callback(~, ~)
        % Is called when the edit_inflation field is changed
        
        set(handles.ui_text_inf_damp_warn_print, 'Visible','Off')
        
        inflate_damp_new = str2double(get(handles.ui_edit_inflation_Damp, 'String'));
        old_damping      = handles.inf_damping;
        
        % Set the inflation value to the update
        handles.inf_damping = inflate_damp_new;
        
        if( not(isfinite(handles.inf_damping)) || handles.inf_damping > 1. )
            % After this, only this edit box will work
            turn_off_controls();
            set(handles.ui_edit_inflation_Damp, 'Enable', 'On');
            
            set(handles.ui_edit_inflation_Damp, 'String', '?', ...
                'FontWeight', 'Bold', 'BackgroundColor', atts.red);
            set(handles.ui_text_inf_damp_err_print, 'Visible','On')
            
            fprintf('ERROR: inflation damping must be between 0.1 and 1.\n')
            fprintf('ERROR: unable to interpret inflation damp value, please try again.\n')
            
            return
            
        end
        
        Update_log_file(rms_time, handles.time, handles.prior_rms, ...
            handles.prior_spread, 'Inflation Damping', old_damping, handles.inf_damping);
        
        % Enable all controls
        turn_on_controls();
        set(handles.ui_edit_inflation_Damp, 'BackgroundColor', 'White', ...
            'FontWeight','Normal');
    end

%% ----------------------------------------------------------------------

    function edit_inflation_Min_Callback(~, ~)
        % Is called when the edit_inflation field is changed
        
        inflate_Min_new = str2double(get(handles.ui_edit_inflation_Min, 'String'));
        old_inf_min     = handles.inf_lower_bound;
        
        % Set the inflation value to the update
        handles.inf_lower_bound = inflate_Min_new;
        
        if(not(isfinite(handles.inf_lower_bound)) || handles.inf_lower_bound < 0. )
            % After this, only this edit box will work
            turn_off_controls();
            set(handles.ui_edit_inflation_Min, 'Enable', 'On');
            
            set(handles.ui_edit_inflation_Min, 'String', '?', ...
                'FontWeight', 'Bold', 'BackgroundColor', atts.red);
            set(handles.ui_text_inf_min_err_print, 'Visible','On')
            
            fprintf('ERROR: inflation lower bound must be greater than or equal 0.\n')
            fprintf('ERROR: unable to interpret inflation min value, please try again.\n')
            
            return
        end
        
        Update_log_file(rms_time, handles.time, handles.prior_rms, ...
            handles.prior_spread, 'Inflation Min', old_inf_min, handles.inf_lower_bound);
        
        % Enable all controls
        turn_on_controls();
        set(handles.ui_edit_inflation_Min, 'BackgroundColor', 'White', ...
            'FontWeight','Normal');
    end

%% ----------------------------------------------------------------------

    function edit_localization_Callback(~, ~)
        % Specifies the action for the 'Localization' text box
        
        localization_new = str2double(get(handles.ui_edit_localization, 'String'));
        old_localization = handles.localization;
        
        % Set the localization value to the update
        handles.localization= localization_new;
        
        if(not(isfinite(handles.localization)) || handles.localization <= 0)
            
            % After this, only this edit box will work
            turn_off_controls();
            set(handles.ui_edit_localization, 'Enable', 'On');
            set(handles.ui_edit_localization, 'String', '?', ...
                'FontWeight', 'Bold', 'BackgroundColor', atts.red );
            set(handles.ui_text_localization_err_print, 'Visible','On')
            
            fprintf('\nERROR: localization must be greater than 0.\n')
            fprintf('ERROR: localization must be greater than 0.\n')
            
            return
        end
        
        Update_log_file(rms_time, handles.time, handles.prior_rms, ...
            handles.prior_spread, 'Localization', old_localization, handles.localization);
        
        % Enable all controls
        turn_on_controls();
        set(handles.ui_edit_localization, 'BackgroundColor', 'White', ...
            'FontWeight','Normal');
        set(handles.ui_text_localization_err_print, 'Visible','Off')
        
        % Update the localization plot
        delete(h_loc)
        h_loc = plot_localization;
    end

%% ----------------------------------------------------------------------

    function edit_ens_size_Callback(~, ~)
        
        % Check to see if the new ensemble size is valid
        new_ens_size = str2double(get(handles.ui_edit_ens_size, 'String'));
        old_ens_size = handles.ens_size;
        
        if(not(isfinite(new_ens_size)) || new_ens_size < 2 || new_ens_size > 40)
            
            % After this, only this edit box will work
            turn_off_controls();
            set(handles.ui_edit_ens_size, 'Enable', 'On');
            
            set(handles.ui_edit_ens_size, 'String', '?', ...
                'FontWeight', 'Bold', ...
                'BackgroundColor', atts.red );
            set(handles.ui_text_ens_size_err_print, 'Visible','On')
            
            fprintf('\nERROR: Must input an integer Ens. Size greater than 1 and less than 40\n');
            fprintf('ERROR: Must input an integer Ens. Size greater than 1 and less than 40\n');
            
            return
            
        end
        
        % Enable all controls
        turn_on_controls();
        
        % clear out the old graphics
        cla(handles.polar_plot)
        cla(handles.timeseries)
        cla(handles.infseries)
        cla(handles.prior_rank_histogram)
        cla(handles.post_rank_histogram)
        
        set(handles.ui_edit_ens_size, 'BackgroundColor', 'White', 'FontWeight', 'Normal');
        
        % Set the ensemble size global value to the update
        handles.ens_size = new_ens_size;
        
        Update_log_file(rms_time, handles.time, handles.prior_rms, ...
            handles.prior_spread, 'Ensemble size', old_ens_size, handles.ens_size);
        
        % Need to reset the ensemble and the time
        clear handles.true_state
        handles.true_state(1, 1:handles.model_size) = TRUE_FORCING;
        handles.true_state(1, 1) = 1.001 * TRUE_FORCING;
        handles.time = 1;
        
        % Generate set of ensemble perturbations
        handles.posterior = zeros(1, handles.model_size, handles.ens_size);
        for imem = 1:handles.ens_size
            handles.posterior(1, 1:handles.model_size, imem) = ...
                handles.true_state(1, :) + 0.001 * randn(1, handles.model_size);
        end
        
        % For convenience make the first prior identical to the first posterior
        handles.prior            = handles.posterior;
        handles.inflate        = ones(1, MODEL_SIZE);
        handles.prior_rms        = 0;
        handles.prior_spread     = 0;
        handles.posterior_rms    = 0;
        handles.posterior_spread = 0;
        
        % An array to keep track of rank histograms
        handles.prior_rank(    1 : handles.ens_size + 1) = 0;
        handles.posterior_rank(1 : handles.ens_size + 1) = 0;
        
        %Reset button to 'Advance Model'
        set(handles.ui_button_Single_Step, 'String' , 'Advance Model');
        handles.ready_to_advance = true;
        
        %Reset the time text
        set(handles.ui_text_time,'String', 'Time = 1');
        
        %Reset the text annotation for RMS
        rms_time = 1;
        
    end

%% ----------------------------------------------------------------------

    function reset_Callback(~, ~)
        % Sets the graphs and the values to original values
        
        % set random number seed to same value to generate known sequences
        % rng('default')  is the Mersenne Twister with seed 0
        rng(0,'twister')
        
        % Initialize the L96 model
        L96          = lorenz_96_static_init_model;
        TRUE_FORCING = L96.forcing;
        FORCING      = L96.forcing;
        MODEL_SIZE   = L96.model_size;
        DELTA_T      = L96.delta_t;
        
        % Set the edit fields
        
        set(handles.ui_edit_localization   , 'String', '0.3');
        set(handles.ui_edit_ens_size       , 'String', '20');
        set(handles.ui_edit_inflation_Std  , 'String', '0.6');
        set(handles.ui_edit_inflation_Damp , 'String', '0.9');
        set(handles.ui_edit_inflation_Min  , 'String', '1.0');
        
        set(handles.ui_button_group_assimilation, ...
            'SelectedObject', handles.ui_radio_noAssimilation);
        handles.filter_type_string = 'No Assimilation';
        set(handles.ui_button_Single_Step, 'String' , 'Advance Model');
        
        set(handles.ui_button_group_obsnetwork,'SelectedObject',handles.ui_radio_full);
        handles.obs_network_string = '1:40:1';
        
        set(handles.ui_button_group_infalgo,'SelectedObject',handles.ui_radio_gauss);
        handles.inf_algorithm_string = 'Gaussian';
        
        set(handles.ui_slider_error       , 'Value' , 10);
        set(handles.ui_edit_forcing       , 'String', 10);
        
        handles.localization = str2double(get(handles.ui_edit_localization,   'String'));
        handles.ens_size     = str2double(get(handles.ui_edit_ens_size,       'String'));
        handles.inflate_sd   = str2double(get(handles.ui_edit_inflation_Std,  'String'));
        handles.inf_damping  = str2double(get(handles.ui_edit_inflation_Damp, 'String'));
        handles.inf_lower_bound = str2double(get(handles.ui_edit_inflation_Min,'String'));
        handles.inf_upper_bound    = 5;
        handles.inf_sd_lower_bound = handles.inflate_sd;
        
        clear handles.true_state
        
        handles.model_size                  = MODEL_SIZE;
        handles.true_state(1, 1:MODEL_SIZE) = TRUE_FORCING;
        handles.true_state(1, 1)            = 1.001 * TRUE_FORCING;
        handles.time                        = 1;
        handles.prior                       = 0;
        handles.prior_rms                   = 0;
        handles.prior_spread                = 0;
        handles.inflate(1, 1:MODEL_SIZE)    = 1;
        handles.posterior                   = 0;
        handles.posterior_rms               = 0;
        handles.posterior_spread            = 0;
        
        set(handles.ui_text_time,'String', sprintf('Time = %d',handles.time))
        
        % Used to normalize the polar plotting
        MEAN_DIST = 35;
        
        handles.h_ens   = [];
        handles.h_truth = [];
        
        % Global semaphore; ready to advance or assimilate?
        handles.ready_to_advance = true;
        
        % Generate set of ensemble perturbations
        handles.posterior = zeros(1, handles.model_size, handles.ens_size);
        for imem = 1:handles.ens_size
            handles.posterior(1, 1:handles.model_size, imem) = ...
                handles.true_state(1, :) + 0.001 * randn(1, handles.model_size);
        end
        
        % For convenience make the first prior identical to the first posterior
        handles.prior = handles.posterior;
        
        % An array to keep track of rank histograms
        handles.prior_rank(    1 : handles.ens_size + 1) = 0;
        handles.posterior_rank(1 : handles.ens_size + 1) = 0;
        
        % Clear out the old graphics. The legends remain, which is nice.
        cla(handles.polar_plot)
        cla(handles.timeseries)
        cla(handles.infseries)
        cla(handles.prior_rank_histogram)
        cla(handles.post_rank_histogram)
        
        % Put back the localization if this is not the initial setup of the graphics
        if(first_call_to_reset)
            first_call_to_reset = false;
        else
            h_loc = plot_localization;
        end
        
        % Reset the RMS averaging:
        rms_time = 1;
        subplot(handles.timeseries); title( ' ' );
        
        if handles.time > 1
            Update_log_file(rms_time, handles.time, handles.prior_rms, ...
                handles.prior_spread, 'RESET');
        end
    end

%% ----------------------------------------------------------------------

    function ClearStats_Callback(~, ~)
        
        % Update log file
        Update_log_file(rms_time, handles.time, handles.prior_rms, ...
            handles.prior_spread, 'Statistics cleared');
        
        % Reset to current time step
        rms_time = handles.time;
        
        % Clear title where RMS and Spread are displayed
        subplot(handles.timeseries)
        str1 = ['Averaging over steps: (' num2str(rms_time) ':n)'];
        str2 = 'Prior: error = ---, spread = ---';
        
        title({str1, str2});
        
        % An array to keep track of rank histograms
        handles.prior_rank(    1 : handles.ens_size + 1) = 0;
        handles.posterior_rank(1 : handles.ens_size + 1) = 0;
        
        % Clear out the old graphics. The legends remain, which is nice.
        cla(handles.prior_rank_histogram)
        cla(handles.post_rank_histogram)
        
    end


%% -----------------------------------------------------------------------------

    function Update_log_file(t1, t2, RMS, AES, info, p1, p2)
        
        logfileid = fopen(LOG_FILE, 'a');
        
        fprintf(logfileid, '# Time step: %d\n', t2);
        fprintf(logfileid, '  >> Statistics over period (%d:%d): avg. Prior RMSE = %.2f, avg. Prior Spread = %.2f\n', ...
            t1, t2, mean(RMS(t1:t2)), mean(AES(t1:t2)));
        
        if strcmp(info, 'Ensemble size') == 1
            fprintf(logfileid, '  $$ User input: %s has been changed from %d to %d\n\n', info, p1, p2);
        elseif strcmp(info, 'Statistics cleared') == 1
            fprintf(logfileid, '  $$ User input: %s; Histograms, RMS and Spread values have been reset\n\n', info);
        elseif strcmp(info, 'RESET') == 1
            fprintf(logfileid, '  $$ User input: %s; Everything has been reset to initial configuration\n\n', info);
        elseif strcmp(info, 'Assimilation Type') == 1
            fprintf(logfileid, '  $$ User input: %s has been changed from `%s` to `%s`\n\n', info, p1, p2);
        elseif strcmp(info, 'Observation Network') == 1
            fprintf(logfileid, '  $$ User input: %s has been changed from `%s` to `%s`\n\n', info, p1, p2);
        elseif strcmp(info, 'Inflation Algorithm') == 1
            fprintf(logfileid, '  $$ User input: %s has been changed from `%s` to `%s`\n\n', info, p1, p2);
        else
            fprintf(logfileid, '  $$ User input: %s has been changed from %.2f to %.2f\n\n', info, p1, p2);
        end
        
        fprintf(logfileid, '     Current configuration:\n');
        fprintf(logfileid, '     - Forcing `F` parameter = %.2f\n', TRUE_FORCING);
        fprintf(logfileid, '     - Assimilation type is `%s`\n', handles.filter_type_string);
        fprintf(logfileid, '     - Observation Network is `%s`\n', handles.obs_network_string);
        fprintf(logfileid, '     - Ensemble size = %d\n', handles.ens_size);
        fprintf(logfileid, '     - Localization = %.2f\n', handles.localization);
        fprintf(logfileid, '     - Inflation Algorithm is `%s`\n', handles.inf_algorithm_string);
        fprintf(logfileid, '     - Adaptive inflation lower bound = %.2f\n', handles.inf_lower_bound);
        fprintf(logfileid, '     - Adaptive inflation upper bound = %.2f\n', handles.inf_upper_bound);
        fprintf(logfileid, '     - Adaptive inflation damping factor = %.2f\n', handles.inf_damping);
        fprintf(logfileid, '     - Adaptive inflation standard deviation = %.2f\n', handles.inflate_sd);
        fprintf(logfileid, '     - Adaptive inflation standard deviation lower bound = %.2f\n\n', handles.inf_sd_lower_bound);
        
        fclose(logfileid);
        
    end

%% ----------------------------------------------------------------------

    function step_ahead()
        % Specifies the action for the [Assimilate Obs/Advance Model] button.
        
        % Test on semaphore, either advance or assimilate
        if(handles.ready_to_advance)
            
            % Set semaphore to indicate that next step may be an assimilation
            % Set the pushbutton text to indicate that next step is an assimilate
            % only if we have selected a filter algorithm
            if( strcmp(handles.filter_type_string, 'No Assimilation'))
                handles.ready_to_advance = true;
                set(handles.ui_button_Single_Step, 'String', 'Advance Model');
            else
                handles.ready_to_advance = false;
                set(handles.ui_button_Single_Step, 'String', 'Assimilate Obs');
            end
            
            % Code for advancing model comes next
            time = handles.time;
            [new_truth, new_time] = lorenz_96_adv_1step(handles.true_state(time, :), ...
                time, TRUE_FORCING);
            handles.time = new_time;
            handles.true_state(new_time, :) = new_truth;
            
            %  Advance the ensemble members; posterior -> new prior
            for imem = 1:handles.ens_size
                [new_ens, new_time] = lorenz_96_adv_1step(handles.posterior(time, :, imem), ...
                    time, FORCING);
                handles.prior(new_time, :, imem) = new_ens;
            end
            
            % Inflate ensemble
            handles.inflate = 1.0 + handles.inf_damping * ( handles.inflate - 1.0 );
            for i = 1:MODEL_SIZE
                ens_mean = mean(handles.prior(new_time, i, :));
                handles.prior(new_time, i, :) = ens_mean + ...
                    sqrt(handles.inflate(1, i)) * (handles.prior(new_time, i, :) - ens_mean);
            end
            
            if(strcmp(handles.filter_type_string, 'No Assimilation'))
                % we are not assimilating
                % just copy prior to posterior
                handles.posterior(new_time, :, :) = handles.prior(new_time, :, :);
            end
            
            % Plot a single invisible point to wipe out the previous plot
            % and maintain axis limits of polar plot.
            axes(handles.polar_plot); cla;
            h_obs = plot_polar([0, 2*pi], 14.9, MEAN_DIST, 'r*', 1);
            set(h_obs, 'Visible', 'Off', 'Color', atts.red)
            
            % Plot the ensemble members (green) and the truth (black)
            for imem = 1:handles.ens_size
                handles.h_ens = plot_polar(polar_y, handles.prior(new_time, :, imem), ...
                    MEAN_DIST, 'g', MODEL_SIZE);
                set(handles.h_ens,'Color',atts.green)
            end
            handles.h_truth = plot_polar(polar_y, new_truth, MEAN_DIST, 'k', MODEL_SIZE);
            % Make truth wider so it is easier to distinguish
            set(handles.h_truth, 'linewidth', 3);
            
            % Plot a graphical indication of the localization halfwidth;
            % is expense of this a problem?
            h_loc = plot_localization;
            
            % Update the time label
            set(handles.ui_text_time, 'String', sprintf('Time = %d', handles.time));
            
            % Compute the prior RMS error of ensemble mean
            handles.prior_rms(new_time) = rms_error(new_truth, handles.prior(new_time, :, :));
            handles.prior_spread(new_time) = ens_spread(handles.prior(new_time, :, :));
            
            % Save the information about the histograms from before
            temp_rank(:, 1) = handles.prior_rank(1:handles.ens_size + 1);
            temp_rank(:, 2) = 0;
            
            % Compute the prior rank histograms
            for i = 1:handles.model_size
                ens_rank = get_ens_rank(squeeze(handles.prior(new_time, i, :)), ...
                    squeeze(new_truth(i)));
                handles.prior_rank(ens_rank) = handles.prior_rank(ens_rank) + 1;
                temp_rank(ens_rank, 2) = temp_rank(ens_rank, 2) + 1;
            end
            
            % Plot the prior_rms error time series
            
            % Change Focus to time evolution of rmse & spread
            subplot(handles.timeseries)
            plot(handles.prior_rms,    '-' ,'Color',atts.green,'LineWidth',2.0);
            plot(handles.prior_spread, '-.','Color',atts.green,'LineWidth',2.0);
            set (handles.timeseries,'YGrid','on')
            
            % Display average RMS of previous run!
            show_rms_on_plot(handles.prior_rms, handles.prior_spread, ...
                [ rms_time, handles.time ]);
            
            % Plot inflation
            subplot(handles.infseries)
            h_inflation = plot(handles.inflate, '-x', 'Color', atts.red ); hold on
            set(handles.infseries, 'YGrid', 'on', 'YLim', [0, 3], 'YTick', [1, 2])
            ylabel('Inflation | Deflation', 'FontSize', atts.fontsize);
            xlabel('Location', 'FontSize', atts.fontsize);
            title( [ 'Overall mean = ' sprintf( '%.4f', mean(handles.inflate) ) ], ...
                'FontSize', 16);
            hold off
            
            % Plot the rank histogram for the prior
            subplot(handles.prior_rank_histogram);
            B = bar(temp_rank, 0.7, 'stacked');
            B(1).FaceColor= atts.blue   ; B(1).EdgeColor= 'k';
            B(2).FaceColor= atts.yellow ; B(2).EdgeColor= 'k';
            axis tight
            
        else % We need to do an assimilation.
            
            % Get current time step
            time = handles.time;
            
            % Generate noisy observations of the truth
            obs_sd = 4;
            obs_error_var = obs_sd^2;
            
            % Get current (prior) inflation
            ss_inflate_base = handles.inflate;
            
            % Do fully sequential assimilation algorithm
            temp_ens = squeeze(handles.prior(time, :, :));
            
            % Select the first plotting box
            axes(handles.polar_plot);
            
            % Figure out the obs network and start assimilation
            [Ny, obs_loc] = get_obs_network(handles.obs_network_string);
            
            % Observe each state variable independently
            obs = NaN * ones(1, MODEL_SIZE);
            for o = 1:Ny
                
                % index of observed variable
                i = obs_loc(o);
                
                obs_prior = temp_ens(i, :);
                obs(i) = handles.true_state(time, i) + obs_sd * randn;
                
                % Compute the increments for observed variable
                switch handles.filter_type_string
                    case 'EAKF'
                        obs_increments = obs_increment_eakf(obs_prior, obs(i), obs_error_var);
                    case 'EnKF'
                        obs_increments = obs_increment_enkf(obs_prior, obs(i), obs_error_var);
                    case 'RHF'
                        obs_increments = obs_increment_rhf (obs_prior, obs(i), obs_error_var);
                    case 'No Assimilation'
                        %No Incrementation
                        obs_increments = 0;
                end
                
                % Regress the increments onto each of the state variables +
                % update inflation values
                for j = 1:MODEL_SIZE
                    [state_incs, r_xy]  = get_state_increments(temp_ens(j, :), ...
                        obs_prior, obs_increments);
                    
                    % Compute distance between obs and state for localization
                    dist = abs(i - j) / MODEL_SIZE;
                    if(dist > 0.5), dist = 1 - dist; end
                    
                    % Compute the localization factor
                    cov_factor = comp_cov_factor(dist, handles.localization);
                    
                    temp_ens(j, :) = temp_ens(j, :) + state_incs * cov_factor;
                    
                    % Get the correlation factor between the observation
                    % and the state variables:
                    gamma = cov_factor * abs(r_xy);
                    
                    % Bayesian update of the inflation
                    ens_obs_mean = mean(obs_prior);
                    ens_obs_var  = var(obs_prior);
                    
                    handles.inflate(j) = update_inflate(ens_obs_mean, ...
                        ens_obs_var, ...
                        obs(i), ...
                        obs_error_var, ...
                        ss_inflate_base(j), ...
                        handles.inflate(j), ...
                        handles.inflate_sd, ...
                        handles.inf_lower_bound, ...
                        handles.inf_upper_bound, ...
                        gamma, ...
                        handles.inf_sd_lower_bound, ...
                        handles.ens_size, ...
                        handles.inf_algorithm_string);
                end
            end
            
            % Plot the observations
            subplot(handles.polar_plot)
            h_obs = plot_polar(polar_y, obs, MEAN_DIST, 'r*', MODEL_SIZE);
            set(h_obs, 'Color', atts.red, 'MarkerSize', 8)
            
            % Update the posterior
            handles.posterior(time, :, :) = temp_ens;
            
            % Compute the posterior rms, spread
            handles.posterior_rms(time) = rms_error(handles.true_state(time, :), ...
                handles.posterior(time, :, :));
            handles.posterior_spread(time) = ens_spread(handles.posterior(time, :, :));
            
            % Save the information about the histograms from before
            temp_rank(:, 1) = handles.posterior_rank(1:handles.ens_size + 1);
            temp_rank(:, 2) = 0;
            
            % Compute the posterior rank histograms
            for i = 1:handles.model_size
                ens_rank = get_ens_rank(squeeze(handles.posterior(time, i, :)), ...
                    squeeze(handles.true_state(time, i)));
                handles.posterior_rank(ens_rank) = handles.posterior_rank(ens_rank) + 1;
                temp_rank(ens_rank, 2) = temp_rank(ens_rank, 2) + 1;
            end
            
            % Plot the posterior_rms error time series
            subplot(handles.timeseries);
            plot(handles.posterior_rms,    '-' , 'Color', atts.blue, 'LineWidth', 2.0);
            plot(handles.posterior_spread, '-.', 'Color', atts.blue, 'LineWidth', 2.0);
            
            % Plot the rank histogram for the prior
            subplot(handles.post_rank_histogram);
            B = bar(temp_rank, 0.7, 'stacked');
            B(1).FaceColor= atts.blue   ; B(1).EdgeColor= 'k';
            B(2).FaceColor= atts.yellow ; B(2).EdgeColor= 'k';
            axis tight
            
            % Set semaphore to indicate that next step is a model advance
            handles.ready_to_advance = true;
            
            % Set the pushbutton text to indicate that the next step is a model advance
            set(handles.ui_button_Single_Step, 'String', 'Advance Model');
            
        end
        
        
    end

%% ----------------------------------------------------------------------

    function ens_mean_rms = rms_error(truth, ens)
        % Calculates the rms_error
        
        ens_mean = mean(squeeze(ens),2)';
        ens_mean_rms = sqrt(sum((truth - ens_mean).^2) / size(truth, 2));
    end

%% ----------------------------------------------------------------------

    function spread = ens_spread(ens)
        % Calculates the ens_spread
        % Remove the mean of each of the 40 model variables (40 locations).
        % resulting matrix is 40x20 ... each row/location is centered (zero mean).
        
        [~, model_size, ens_size] = size(ens);
        datmat = detrend(squeeze(ens)','constant'); % remove the mean of each location.
        denom  = (model_size - 1)*ens_size;
        spread = sqrt(sum(datmat(:).^2) / denom);
    end

%% ----------------------------------------------------------------------

    function show_rms_on_plot(prior_rms_vals, prior_aes_vals, ranges)
        
        subplot(handles.timeseries)
        
        prior_rms_vals_new = prior_rms_vals(ranges(1): ranges(2));
        prior_aes_vals_new = prior_aes_vals(ranges(2): ranges(2));
        
        str1 = ['Averaging over steps: (' num2str(ranges(1)) ':' num2str(ranges(2)) ')'];
        str2 = ['Prior: error = ' ...
            sprintf('%.2f', mean(prior_rms_vals_new) ) ', spread = ' ...
            sprintf('%.2f', mean(prior_aes_vals_new) )];
        
        show_old_rms = title( {str1, str2});
        set(show_old_rms, 'FontSize', 16)
    end

%% ----------------------------------------------------------------------

    function turn_off_controls()
        % Disables all the buttons,menus, and edit fields
        
        set(handles.ui_button_Single_Step,      'Enable', 'Off');
        set(handles.ui_button_Auto_Run,         'Enable', 'Off');
        
        % In 2015, there is a way to disable the button group, but it is not
        % compatible with 2014, so we must enable/disable each radio button
        % seperately
        set(handles.ui_radio_noAssimilation,    'Enable', 'Off');
        set(handles.ui_radio_EAKF,              'Enable', 'Off');
        set(handles.ui_radio_EnKF,              'Enable', 'Off');
        set(handles.ui_radio_RHF,               'Enable', 'Off');
        
        set(handles.ui_radio_full,              'Enable', 'Off');
        set(handles.ui_radio_half,              'Enable', 'Off');
        set(handles.ui_radio_quarter,           'Enable', 'Off');
        set(handles.ui_radio_first20,           'Enable', 'Off');
        set(handles.ui_radio_mid20,             'Enable', 'Off');
        set(handles.ui_radio_firstlast10,       'Enable', 'Off');
        
        set(handles.ui_slider_error,            'Enable', 'Off');
        set(handles.ui_edit_forcing,            'Enable', 'Off');
        set(handles.ui_edit_localization,       'Enable', 'Off');
        set(handles.ui_edit_ens_size,           'Enable', 'Off');
        
        set(handles.ui_radio_gauss,             'Enable', 'Off');
        set(handles.ui_radio_gamma,             'Enable', 'Off');
        set(handles.ui_edit_inflation_Std,      'Enable', 'Off');
        set(handles.ui_edit_inflation_Damp,     'Enable', 'Off');
        set(handles.ui_edit_inflation_Min,      'Enable', 'Off');
        
        set(handles.ResetButton,                'Enable', 'Off');
        set(handles.ClearStats,                 'Enable', 'Off');
    end

%% ----------------------------------------------------------------------

    function turn_on_controls()
        % Enables all the buttons,menus, and edit fields
        
        set(handles.ui_button_Single_Step,      'Enable', 'On');
        set(handles.ui_button_Auto_Run,         'Enable', 'On');
        
        % In 2015, there is a way to disable the button group,
        % but it is not compatible with 2014, so we must enable/disable
        % each radio button seperately
        set(handles.ui_radio_noAssimilation,    'Enable', 'On');
        set(handles.ui_radio_EAKF,              'Enable', 'On');
        set(handles.ui_radio_EnKF,              'Enable', 'On');
        set(handles.ui_radio_RHF,               'Enable', 'On');
        
        set(handles.ui_radio_full,              'Enable', 'On');
        set(handles.ui_radio_half,              'Enable', 'On');
        set(handles.ui_radio_quarter,           'Enable', 'On');
        set(handles.ui_radio_first20,           'Enable', 'On');
        set(handles.ui_radio_mid20,             'Enable', 'On');
        set(handles.ui_radio_firstlast10,       'Enable', 'On');
        
        set(handles.ui_slider_error,            'Enable', 'On');
        set(handles.ui_edit_forcing,            'Enable', 'On');
        set(handles.ui_edit_localization,       'Enable', 'On');
        set(handles.ui_edit_ens_size,           'Enable', 'On');
        
        set(handles.ui_radio_gauss,             'Enable', 'On');
        set(handles.ui_radio_gamma,             'Enable', 'On');
        set(handles.ui_edit_inflation_Std,      'Enable', 'On');
        set(handles.ui_edit_inflation_Damp,     'Enable', 'On');
        set(handles.ui_edit_inflation_Min,      'Enable', 'On');
        
        set(handles.ResetButton,                'Enable', 'On');
        set(handles.ClearStats,                 'Enable', 'On');
    end

%% -----------------------------------------------------------------------

    function Assimilation_selection(~, eventdata)
        %eventdata refers to the data in the GUI when a radio button in the
        %group is changed
        %Set the filter_type_string to newest radiobutton Value
        
        filter_new = get(eventdata.NewValue,'String');
        old_filter = handles.filter_type_string;
        
        handles.filter_type_string = filter_new;
        
        Update_log_file(rms_time, handles.time, handles.prior_rms, handles.prior_spread, ...
            'Assimilation Type', old_filter, handles.filter_type_string);
    end

%% -----------------------------------------------------------------------

    function ObsNetwork_selection(~, eventdata)
        %eventdata refers to the data in the GUI when a radio button in the
        %group is changed
        %Set the obs_network_string to newest radiobutton Value
        
        network_new = get(eventdata.NewValue,'String');
        old_network = handles.obs_network_string;
        
        handles.obs_network_string = network_new;
        
        Update_log_file(rms_time, handles.time, handles.prior_rms, handles.prior_spread, ...
            'Observation Network', old_network, handles.obs_network_string);
    end

%% -----------------------------------------------------------------------

    function InfAlgo_selection(~, eventdata)
        %eventdata refers to the data in the GUI when a radio button in the
        %group is changed
        %Set the inf_algorithm_string to newest radiobutton Value
        
        infalgo_new = get(eventdata.NewValue,'String');
        old_infalgo = handles.inf_algorithm_string;
        
        handles.inf_algorithm_string = infalgo_new;
        
        Update_log_file(rms_time, handles.time, handles.prior_rms, handles.prior_spread, ...
            'Inflation Algorithm', old_infalgo, handles.inf_algorithm_string);
    end

%% -----------------------------------------------------------------------

    function [Ny, obs_loc] = get_obs_network(user_defined_network)
        
        switch user_defined_network
            case '1:40:1'
                obs_loc = 1:MODEL_SIZE;
                
            case '1:40:2'
                obs_loc = 2:2:MODEL_SIZE;
                
            case '1:40:4'
                obs_loc = 2:4:MODEL_SIZE;
                
            case '1:20'
                obs_loc = 2:21;
                
            case '10:30'
                obs_loc = 11:31;
                
            case '1:10; 30:40'
                obs_loc = [1:11, 31:MODEL_SIZE];
                
        end
        Ny = length(obs_loc);
        
    end

%% -----------------------------------------------------------------------

    function edit_forcing_Callback(~,~)
        % This function is called whenever the edit field beneath the
        % slider is changed, the slider and the edit field are connected,
        % the edit field is simply used to allow for more precise forcing
        % values
        
        % Undo any changes that could have been made by erros
        turn_on_controls();
        set(handles.ui_edit_forcing, 'BackgroundColor' , 'White');
        set(handles.ui_edit_forcing, 'FontWeight' , 'Normal');
        
        if(isfinite(str2double(get(handles.ui_edit_forcing, 'String'))))
            if (str2double(get(handles.ui_edit_forcing, 'String')) >= 4 && ...
                    str2double(get(handles.ui_edit_forcing, 'String')) <= 12)
                FORCING = str2double(get(handles.ui_edit_forcing, 'String'));
                set(handles.ui_slider_error, 'Value' , FORCING);
                
                % Fix everything created by a potential previous error
                turn_on_controls();
                set(handles.ui_edit_forcing, 'BackgroundColor' , 'White');
                return;
            elseif (str2double(get(handles.ui_edit_forcing, 'String')) <= 4)
                % There is an error
                % Prevent 2 errors at once
                turn_off_controls();
                set(handles.ui_edit_forcing,'Enable', 'On');
                set(handles.ui_edit_forcing, 'String' , '>4');
                set(handles.ui_edit_forcing, 'BackgroundColor' , atts.red);
                set(handles.ui_edit_forcing, 'FontWeight' , 'Bold');
                fprintf('ERROR: Number must be greater than or equal to 4\n');
                return;
            elseif (str2double(get(handles.ui_edit_forcing, 'String')) >= 12)
                % There is an error
                % Prevent 2 errors at once
                turn_off_controls();
                set(handles.ui_edit_forcing,'Enable', 'On');
                set(handles.ui_edit_forcing, 'String' , '<12');
                set(handles.ui_edit_forcing, 'BackgroundColor' , atts.red);
                set(handles.ui_edit_forcing, 'FontWeight' , 'Bold');
                fprintf('ERROR: Number must be less than or equal to 12\n');
                return
            end
        else
            % There is an error, not a number
            % Prevent 2 errors at once
            turn_off_controls();
            set(handles.ui_edit_forcing,'Enable', 'On');
            set(handles.ui_edit_forcing, 'String' , '?');
            set(handles.ui_edit_forcing, 'BackgroundColor' , atts.red);
            set(handles.ui_edit_forcing, 'FontWeight' , 'Bold');
            
            fprintf('ERROR: Must enter a number between 4 and 12\n');
            return;
        end
    end

%% -----------------------------------------------------------------------

    function my_h_loc = plot_localization
        % Plot a graphical indication of the localization halfwidth
        subplot(handles.polar_plot);
        dist = 16;
        
        % Localization is in halfwidth, fraction of domain (NOT RADIANS AS IN 3D MODELS).
        % Convert to halfwidth in radians for plotting
        half_radians = handles.localization * 2 * pi;
        
        % Plot 4 ranges
        my_h_loc   = zeros(1, 4);
        my_col_loc = atts.colors4loc;
        for ipl = 1:4
            ymax = min([half_radians * (5.-ipl) / 2., pi]);
            ymin = -ymax;
            % Use 40 points for each range
            y = ymin:ymax/20:ymax;
            
            % index the plots, so when you delete it (you delete all
            % localization lines) not only the last one!
            my_h_loc(ipl) = polar_dares(y, dist*ones(size(y)));     
            
            hold on
            
            % Lines get wider for larger localization
            set(my_h_loc(ipl), 'linewidth', 2*ipl, 'Color', my_col_loc(ipl, :));
        end
        
        % Plot a label for the localization graphic
        h_loc_text = text(-12, 0, 'Localization');
        set(h_loc_text, 'color', 'k', 'fontsize', 15, 'fontweight', 'bold');
        
        % Plot an observation asterisk
        plot(dist, 0, '*', ...
            'MarkerSize', 12, ...
            'MarkerFaceColor', atts.red, ...
            'MarkerEdgeColor', atts.red);
    end

end

