function varargout = markROIs(varargin)
% MARKROIS MATLAB code for markROIs.fig
%      MARKROIS, by itself, creates a new MARKROIS or raises the existing
%      singleton*.
%
%      H = MARKROIS returns the handle to a new MARKROIS or the handle to
%      the existing singleton*.
%
%      MARKROIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARKROIS.M with the given input arguments.
%
%      MARKROIS('Property','Value',...) creates a new MARKROIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before markROIs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to markROIs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help markROIs

% Last Modified by GUIDE v2.5 30-May-2016 22:10:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @markROIs_OpeningFcn, ...
                   'gui_OutputFcn',  @markROIs_OutputFcn, ...
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


% --- Executes just before markROIs is made visible.
function markROIs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to markROIs (see VARARGIN)

% Choose default command line output for markROIs
handles.output = hObject;

handles.expInfo = varargin{1};
handles = loadData(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes markROIs wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function handles = loadData(handles)
totalFrames = handles.expInfo.totalFrames;
set(handles.slider_frames,'Max',totalFrames);
set(handles.slider_frames,'Min',1);
minStep = 1/totalFrames;
maxStep = 10/totalFrames;
set(handles.slider_frames,'SliderStep',[minStep maxStep]);
set(handles.slider_frames,'value',1);
metaData.playFlag = 0;
set(handles.pushbutton_play,'userData',metaData);
set(handles.pushbutton_stop,'visible','off');
loadFrame(handles,1);


function loadFrame(handles,frameNumber)
if frameNumber == 0 % for loading average image
    load(makeName('averageImageMC.mat',handles.expInfo.pDataFolder));
    set(handles.figure1,'currentAxes',handles.axes_main);
    thisFrame = averageImage;
end
if frameNumber > 0
    fid = fopen(handles.expInfo.rawFile,'r');
    fseek(fid, (frameNumber-1)*handles.expInfo.pixelX*handles.expInfo.pixelY*2, 'bof');
    thisFrame = fread(fid,handles.expInfo.pixelX*handles.expInfo.pixelY,'uint16=>double',0,'l'); 
    fclose(fid);
    thisFrame = reshape(thisFrame,handles.expInfo.pixelY,handles.expInfo.pixelX);
    thisFrame = thisFrame';    
end
imagesc(thisFrame);
axis off; axis equal;
colormap gray;
set(handles.text_frameNumber,'String',sprintf('%d of %d',frameNumber,handles.expInfo.totalFrames));


% --- Outputs from this function are returned to the command line.
function varargout = markROIs_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider_frames_Callback(hObject, eventdata, handles)
% hObject    handle to slider_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameNumber = round(get(hObject,'Value'));
set(hObject,'Value',frameNumber);
loadFrame(handles,frameNumber);


% --- Executes during object creation, after setting all properties.
function slider_frames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton_play.
function pushbutton_play_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'visible','off');
set(handles.pushbutton_stop,'visible','on');
metaData = get(handles.axes_main,'userData');
metaData.playFlag = 1;
set(hObject,'userData',metaData);
currentFrame = round(get(handles.slider_frames,'Value'));
maxFrame = get(handles.slider_frames,'Max');

while 1
   metaData = get(hObject,'userData');
    if metaData.playFlag == 0
        break;
    end
    if currentFrame <maxFrame
        currentFrame = currentFrame + 1;
        loadFrame(handles,currentFrame);
    else
        pushbutton_stop_Callback(handles.pushbutton_stop, eventdata, handles)
        break;
    end
    pause(0.01);
end

% --- Executes on button press in pushbutton_stop.
function pushbutton_stop_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
metaData = get(handles.pushbutton_play,'userData');
metaData.playFlag = 0;
set(handles.pushbutton_play,'userData',metaData);
set(hObject,'visible','off');
set(handles.pushbutton_play,'visible','on');
pause(0.3);