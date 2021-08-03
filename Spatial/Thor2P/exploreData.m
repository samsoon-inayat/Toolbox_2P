function varargout = exploreData(varargin)
% EXPLOREDATA MATLAB code for exploreData.fig
%      EXPLOREDATA, by itself, creates a new EXPLOREDATA or raises the existing
%      singleton*.
%
%      H = EXPLOREDATA returns the handle to a new EXPLOREDATA or the handle to
%      the existing singleton*.
%
%      EXPLOREDATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EXPLOREDATA.M with the given input arguments.
%
%      EXPLOREDATA('Property','Value',...) creates a new EXPLOREDATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before exploreData_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to exploreData_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help exploreData

% Last Modified by GUIDE v2.5 27-Jun-2017 16:31:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @exploreData_OpeningFcn, ...
                   'gui_OutputFcn',  @exploreData_OutputFcn, ...
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


% --- Executes just before exploreData is made visible.
function exploreData_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to exploreData (see VARARGIN)

% Choose default command line output for exploreData
handles.output = hObject;
handles.ei = varargin{1};
ii = varargin{2};
if nargin == 6
    binning = varargin{3};
else
    binning = 1;
end
handles.totalFiles = handles.ei.totalFrames;
% set(handles.slider1,'userdata',handles.fid);
handles.binning = binning;
img = getImg(handles,1);
handles.fileNames = varargin{1};
% set(handles.slider1,'userdata',varargin{1});
imagesc(handles.axes1,img);
colormap(gray);
txt = sprintf('%s - %s - %s',handles.ei.db(ii).mouse_name,handles.ei.db(ii).date,handles.ei.db(ii).expText);
set(handles.folderName,'String',txt);
fct = sprintf('%d/%d',1,handles.totalFiles);
set(handles.frameCounter,'String',fct);
min_step = 1/handles.totalFiles;
max_step = 10/handles.totalFiles;
set(handles.slider1,'Max',handles.totalFiles);
set(handles.slider1,'Min',1);
set(handles.slider1,'SliderStep',[min_step max_step]);
set(handles.slider1,'value',1);
addlistener(handles.slider1,'ContinuousValueChange',@slider_frames_Callback1);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes exploreData wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = exploreData_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function frame = getImg(handles,frameNumber)
binning = handles.binning;
fid = fopen(handles.ei.rawFile);
row = handles.ei.pixelY;
col = handles.ei.pixelX;
fseek(fid, (frameNumber-1)*row*col*2, 'bof');
I1 = fread(fid,row*col,'uint16=>double',0,'l'); 
Z = reshape(I1,row,col); 
if binning > 1
    frame = imresize(Z',1/binning);
else
    frame = Z';
end
fclose(fid);

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
fileNumber = round(get(hObject,'Value'));
set(hObject,'Value',fileNumber);
img = getImg(handles,fileNumber);
handles.axes1;
imagesc(img);
fct = sprintf('%d/%d',fileNumber,handles.totalFiles);
set(handles.frameCounter,'String',fct);

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
function slider_frames_Callback1(hObject, eventdata, handles)
% hObject    handle to slider_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles = guidata(hObject);
fileNumber = round(get(hObject,'Value'));
set(hObject,'Value',fileNumber);
img = getImg(handles,fileNumber);
imagesc(handles.axes1,img);
fct = sprintf('%d/%d',fileNumber,handles.totalFiles);
set(handles.frameCounter,'String',fct);
refresh(handles.figure1);


% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
type = get(hObject,'String')
colormap(type);
