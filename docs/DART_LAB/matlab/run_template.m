function varargout = run_template(varargin)
%% RUN_TEMPLATE M-file for run_template.fig
%      RUN_TEMPLATE, by itself, creates a new RUN_TEMPLATE or raises the existing
%      singleton*.
%
%      H = RUN_TEMPLATE returns the handle to a new RUN_TEMPLATE or the handle to
%      the existing singleton*.
%
%      RUN_TEMPLATE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RUN_TEMPLATE.M with the given input arguments.
%
%      RUN_TEMPLATE('Property','Value',...) creates a new RUN_TEMPLATE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before run_template_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to run_template_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help run_template

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @run_template_OpeningFcn, ...
    'gui_OutputFcn',  @run_template_OutputFcn, ...
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



%% --- Executes just before run_template is made visible.
function run_template_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to run_template (see VARARGIN)

% set random number seed to same value to generate known sequences
rng('default')

% Choose default command line output for run_template
handles.output = hObject;

% Global semaphore; ready to advance or assimilate?
handles.ready_to_advance = true;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes run_template wait for user response (see UIRESUME)
% uiwait(handles.figure1);



%% --- Outputs from this function are returned to the command line.
function varargout = run_template_OutputFcn(~, ~, handles)
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



%% --- Executes on button press in pushbutton_single_step.
function pushbutton_single_step_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_single_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Move the model ahead a step or assimilate observations as appropriate
step_ahead(hObject, handles);



%% --- Executes on button press in pushbutton_free_run.
function pushbutton_free_run_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_free_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Turn off all the other model status controls to avoid a mess
% MAKE SURE TO INCLUDE OTHER CONTROLS HERE
set(handles.pushbutton_single_step, 'Enable', 'Off');

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



%%----------- Moves the model ahead or assimilates next observations ------
function step_ahead(hObject, handles)

% Test on semaphore, either advance or assimilate
if(handles.ready_to_advance)
    % Set semaphore to indicate that next step is an assimilation
    handles.ready_to_advance = false;
    
    % Set the pushbutton text to indicate that next step is an assimilate
    set(handles.pushbutton_single_step, 'String', 'Assimilate Obs');
    
    % Code for advancing model comes next (delete two exisiting lines)
    pause(2)
    plot([1 2], [2 1], 'r');
    
else
    % Set semaphore to indicate that next step is a model advance
    handles.ready_to_advance = true;
    
    % Set the pushbutton text to indicate that the next step is a model advance
    set(handles.pushbutton_single_step, 'String', 'Advance Model');
    
    % Code for doing the assimilation comes here (delete two exisiting lines)
    pause(2)
    plot([1 2], [2 1], 'b');
    
end

% If using multiple windows might need to reset focus to the gui window here
[~, gcbo_fig] = gcbo;
figure(gcbo_fig);

% Update the global storage and return
guidata(hObject, handles);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
