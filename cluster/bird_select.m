function varargout = bird_select(varargin)
% BIRD_SELECT M-file for bird_select.fig
%      BIRD_SELECT, by itself, creates a new BIRD_SELECT or raises the existing
%      singleton*.
%
%      H = BIRD_SELECT returns the handle to a new BIRD_SELECT or the handle to
%      the existing singleton*.
%
%      BIRD_SELECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BIRD_SELECT.M with the given input arguments.
%
%      BIRD_SELECT('Property','Value',...) creates a new BIRD_SELECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bird_select_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bird_select_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bird_select

% Last Modified by GUIDE v2.5 06-Jan-2008 18:03:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bird_select_OpeningFcn, ...
                   'gui_OutputFcn',  @bird_select_OutputFcn, ...
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


% --- Executes just before bird_select is made visible.
function bird_select_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bird_select (see VARARGIN)

% Choose default command line output for bird_select
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes bird_select wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = bird_select_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
env_local

value = get( handles.cube_n,'value' ) ;
str = get( handles.cube_n,'string' ) ;

file = str2num( str{value} ) ;
tm=T(file);
rog=RHOG(file);
% vx=Vx(file);
% vy=Vy(file);
% vz=Vz(file);

data.tm = tm ;
data.rog = rog ;
% data.vx=single(vx);data.vy=single(vy);data.vz=single(vz);

set( handles.axes1, 'userdata', data ) ;
disp('Done!' ) ;


% --- Executes on selection change in cube_n.
function cube_n_Callback(hObject, eventdata, handles)
% hObject    handle to cube_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns cube_n contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cube_n
%load_Callback

% --- Executes during object creation, after setting all properties.
function cube_n_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cube_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = get( handles.axes1, 'userdata' ) ;
BB = data.tm.*data.rog;
AA = squeeze( mean( bb(:,100:150,:),2) ) ;
imagesc( AA, 'parent', handles.axes1 );
set( handles.axes1, 'userdata', data ) ;

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


%g = get( handles.slider1, 'value' ) ;

%data = get( handles.axes1, 'userdata' ) ;
%lim1 = ceil(g*246) ;
%tmw = data.tm.*data.rog;
%AA = squeeze( mean( tmw(:,lim1:lim1+10,:),2) ) ;
%imagesc( AA, 'parent', handles.axes1 );
%data.tm = single(data.tm) ;
%data.rog = single(data.rog) ;

%set( handles.axes1, 'userdata', data ) ;
make_cmap(handles)


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
make_cmap(handles)

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in orient.
function orient_Callback(hObject, eventdata, handles)
% hObject    handle to orient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns orient contents as cell array
%        contents{get(hObject,'Value')} returns selected item from orient


% --- Executes during object creation, after setting all properties.
function orient_CreateFcn(hObject, eventdata, handles)
% hObject    handle to orient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function make_cmap(handles)
%recieves all parameters and makes color map
data = get( handles.axes1, 'userdata' ) ;

len=size(data.A,3);
wdth = ceil(len.*get( handles.slider2, 'value' )) ; %width of slice
strt = ceil((len-wdth).*get( handles.slider1, 'value' )) ; %starting point from bottom

AA = squeeze( mean( data.A(:,:,strt:strt+wdth),3) ) ;
imagesc( AA, 'parent', handles.axes1 );
set( handles.axes1, 'userdata', data ) ;


% --- Executes on selection change in field.
function field_Callback(hObject, eventdata, handles)
% hObject    handle to field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns field contents as cell array
%        contents{get(hObject,'Value')} returns selected item from field

value = get( handles.field,'value' ) ;
str = get( handles.field,'string' ) ;
data = get( handles.axes1, 'userdata' ) ;

switch str{value}
    case 'T' 
        data.A=data.tm.*data.rog;
%    case 's' data.A=data.tm./(data.rog.^(2./3));
%    case 'v_r' data.A=
%    case 'rho_g' data.A=data.rog;
%    case 'rho_dm' data.A=data.rodm;   
%    case 'rho_g' data.A=data.rotot;
    otherwise
        disp('no field selected');
end  
set( handles.axes1, 'userdata', data ) ;

% --- Executes during object creation, after setting all properties.
function field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in vfield.
function vfield_Callback(hObject, eventdata, handles)
% hObject    handle to vfield (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vfield


