function varargout = follow_flux(varargin)
% follow_flux MATLAB code for follow_flux.fig
%      follow_flux, by itself, creates a new follow_flux or raises the existing
%      singleton*.
%
%      H = follow_flux returns the handle to a new follow_flux or the handle to
%      the existing singleton*.
%
%      follow_flux('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in follow_flux.M with the given input arguments.
%
%      follow_flux('Property','Value',...) creates a new follow_flux or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before follow_flux_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to follow_flux_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help follow_flux

% Last Modified by GUIDE v2.5 17-Apr-2013 16:22:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @follow_flux_OpeningFcn, ...
    'gui_OutputFcn',  @follow_flux_OutputFcn, ...
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

global cf8;
global cf4;
global cf2;
global cf1;
global boxx;
global mesh_phi
global mesh_theta
global step 



boxx=8;

step=3;


% End initialization code - DO NOT EDIT


% --- Executes just before follow_flux is made visible.
function follow_flux_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to follow_flux (see VARARGIN)

% Choose default command line output for follow_flux
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes follow_flux wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = follow_flux_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function txtParamCL_Callback(hObject, eventdata, handles)
% hObject    handle to txtParamCL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtParamCL as text
%        str2double(get(hObject,'String')) returns contents of txtParamCL as a double


% --- Executes during object creation, after setting all properties.
function txtParamCL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtParamCL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtParamAexp_Callback(hObject, eventdata, handles)
% hObject    handle to txtParamAexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtParamAexp as text
%        str2double(get(hObject,'String')) returns contents of txtParamAexp as a double


% --- Executes during object creation, after setting all properties.
function txtParamAexp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtParamAexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnSetParam.
function btnSetParam_Callback(hObject, eventdata, handles)
% hObject    handle to btnSetParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global cf1
global cf2
global cf4
global cf8
global iIndmax


ParamCL = get(handles.txtParamCL,'String');
ParamAexp = get(handles.txtParamAexp,'String');

new_env(ParamCL,'csf',ParamAexp,'win');
cf1=flux_sphere(1);
cf2=flux_sphere(2);
cf4=flux_sphere(4);
cf8=flux_sphere(8);

global NCELL
global mesh_phi
global mesh_theta
global boxx
iIndmax=NCELL*5/2;
[~,mesh_phi,mesh_theta] = sphere_grid(boxx);




% --- Executes on button press in btnUp.
function btnUp_Callback(hObject, eventdata, handles)
% hObject    handle to btnUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global step

iCurrIndex = str2num(get(handles.txtIndex,'String'));
iCurrIndex = iCurrIndex + step;
UpadteIndex(iCurrIndex, hObject, eventdata, handles);
% if (iCurrIndex < 100)
%     iCurrIndex = iCurrIndex + 1;
%     set(handles.txtIndex,'String', num2str(iCurrIndex));
%     set(handles.txtRadius,'String', num2str(iCurrIndex));
% end

% --- Executes on button press in btnDown.
function btnDown_Callback(hObject, eventdata, handles)
% hObject    handle to btnDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global step
iCurrIndex = str2num(get(handles.txtIndex,'String'));
iCurrIndex = iCurrIndex - step;
UpadteIndex(iCurrIndex, hObject, eventdata, handles);
% if (iCurrIndex > 1)
%     iCurrIndex = iCurrIndex - 1;
%     set(handles.txtIndex,'String', num2str(iCurrIndex));
%     set(handles.txtRadius,'String', num2str(iCurrIndex));
% end


function txtIndex_Callback(hObject, eventdata, handles)
% hObject    handle to txtIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtIndex as text
%        str2double(get(hObject,'String')) returns contents of txtIndex as a double
iCurrIndex = str2num(get(handles.txtIndex,'String'));

UpadteIndex(iCurrIndex, hObject, eventdata, handles);
% if (iCurrIndex >= 1 & iCurrIndex <= 100)
%     set(handles.txtRadius,'String', num2str(iCurrIndex));
% end

% --- Executes during object creation, after setting all properties.
function txtIndex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtRadius_Callback(hObject, eventdata, handles)
% hObject    handle to txtRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtRadius as text
%        str2double(get(hObject,'String')) returns contents of txtRadius as a double


% --- Executes during object creation, after setting all properties.
function txtRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function UpadteIndex(iIndex, hObject, eventdata, handles)
global iIndmax
if (iIndex >= 1 & iIndex <= iIndmax)
    set(handles.txtIndex,'String', num2str(iIndex));
    % set(handles.txtRadius,'String', num2str(iIndex));
   
    global cf1
    global cf2
    global cf4
    global cf8
    global NCELL
    global zred
    global hub
    
    
    if iIndex<NCELL
        ind=iIndex;
        boxx=1;
        shel=squeeze(sum(cf1(ind-1:ind+1,:,:)));
    elseif iIndex <NCELL*3/2
        ind=iIndex-NCELL/2;
        boxx=2;
        shel=squeeze(sum(cf2(ind-1:ind+1,:,:)));
    elseif iIndex < NCELL*2
        ind=iIndex-NCELL;
        boxx=4;
        shel=squeeze(sum(cf4(ind-1:ind+1,:,:)));
    else
        ind=iIndex-NCELL*3/2;
        boxx=8;
        shel=squeeze(sum(cf8(ind-1:ind+1,:,:),1));
    end
        
    
    shel(shel>0)=1e-30;
    rad=ind/NCELL*0.5*boxx/(1+zred);
    
    set(handles.txtRadius,'String', rad);
    set(handles.txtRadRv,'String', rad/hub/get_rvir);
    
    global clim
    clim(1) = str2num(get(handles.txtClim1,'String'));
    clim(2) = str2num(get(handles.txtClim2,'String'));
   
    global mesh_phi
    global mesh_theta
    mphi=squeeze(mesh_phi(ind,:,:));
    mthet=squeeze(mesh_theta(ind,:,:));
    axes(handles.axes2);
    plot_hammer_shell(log10(abs(shel)),mthet,mphi,'nc',500);
    caxis(clim)
    
    axes(handles.axes1);
    plot_sphere_shell(log10(abs(shel)),mthet,mphi,'nc',300,'colorbar');
    axis equal;
    xlabel('X');ylabel('Y');zlabel('Z')
     caxis(clim)
     
     pack;
    
end


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

if (strcmp(eventdata.Key,'downarrow'))
    btnDown_Callback(hObject, eventdata, handles);
elseif (strcmp(eventdata.Key,'uparrow'))
    btnUp_Callback(hObject, eventdata, handles)
end


% --- Executes on button press in btnRotate.
function btnRotate_Callback(hObject, eventdata, handles)
% hObject    handle to btnRotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 rotate3d(handles.axes1);
% Hint: get(hObject,'Value') returns toggle state of btnRotate



function txtRadRv_Callback(hObject, eventdata, handles)
% hObject    handle to txtRadRv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtRadRv as text
%        str2double(get(hObject,'String')) returns contents of txtRadRv as a double


% --- Executes during object creation, after setting all properties.
function txtRadRv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtRadRv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtClim1_Callback(hObject, eventdata, handles)
% hObject    handle to txtClim1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtClim1 as text
%        str2double(get(hObject,'String')) returns contents of txtClim1 as a double


% --- Executes during object creation, after setting all properties.
function txtClim1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtClim1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtClim2_Callback(hObject, eventdata, handles)
% hObject    handle to txtClim2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtClim2 as text
%        str2double(get(hObject,'String')) returns contents of txtClim2 as a double


% --- Executes during object creation, after setting all properties.
function txtClim2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtClim2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
