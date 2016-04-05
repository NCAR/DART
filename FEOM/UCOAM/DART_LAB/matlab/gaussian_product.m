function varargout = gaussian_product(varargin)
%% GAUSSIAN_PRODUCT demonstrates the product of two gaussian distributions.
%
%    This is fundamental to Kalman filters and to ensemble
%    data assimilation. Change the parameters of the
%    gaussian for the Prior (green) and the Observation (red)
%    and click on 'Plot Posterior'.
%
%    The product (in this case, the 'Posterior') of two gaussians is a gaussian.
%    If the parameters of the two gaussians are known, the parameters of the 
%    resulting gaussian can be calculated.
%
% See also: oned_model, oned_ensemble, twod_ensemble,
%           run_lorenz_63, run_lorenz_96

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gaussian_product_OpeningFcn, ...
                   'gui_OutputFcn',  @gaussian_product_OutputFcn, ...
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


% --- Executes just before gaussian_product is made visible.
function gaussian_product_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gaussian_product (see VARARGIN)

help gaussian_product

% set random number seed to same value to generate known sequences
rng('default')

% Choose default command line output for gaussian_product
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Plot the initial prior and observation likelihood pdf's
h = guihandles;
g_prod_plot(h);

% UIWAIT makes gaussian_product wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gaussian_product_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(~, ~, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

g_prod_plot(handles);


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

function edit2_Callback(~, ~, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

g_prod_plot(handles);


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



function edit3_Callback(~, ~, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double

g_prod_plot(handles);


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



function edit4_Callback(~, ~, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double

g_prod_plot(handles);


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, ~, ~)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(~, ~, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Need to replot prior and obs, then compute posterior and plot
[prior_mean, prior_sd, obs_mean, obs_err_sd, is_err] = g_prod_plot(handles);

% If there is an error, zero out the posterior text values
% don't try to do posterior computation
if(is_err)
    set(handles.text7, 'String', strcat('Posterior Mean = '));
    set(handles.text8, 'String', strcat('Posterior SD = '));
    set(handles.text9, 'String', strcat('Weight = '));
    return;
end

% Compute the posterior mean, sd and weight
[post_mean, post_sd, weight] = ...
    product_of_gaussians(prior_mean, prior_sd, obs_mean, obs_err_sd);
post_handle = plot_gaussian(post_mean, post_sd, 1);
set(post_handle, 'Color', 'b', 'LineWidth', 2);

% Print values
set(handles.text7, 'String', ['Posterior Mean = ', num2str(post_mean)]);
set(handles.text8, 'String', ['Posterior SD = ', num2str(post_sd)]);

% Also plot the weighted posterior as dashed
post_handle = plot_gaussian(post_mean, post_sd, weight);
set(post_handle, 'Color', 'b', 'Linestyle', '--');
set(handles.text9, 'String', ['Weight = ', num2str(weight)]);

legend('Prior', 'Obs. Likelihood', 'Posterior', 'Weighted Posterior');


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

