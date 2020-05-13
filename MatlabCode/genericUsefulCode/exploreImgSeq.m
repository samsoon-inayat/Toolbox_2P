function varargout = exploreImgSeq(varargin)
% EXPLOREIMGSEQ MATLAB code for exploreImgSeq.fig
%      EXPLOREIMGSEQ, by itself, creates a new EXPLOREIMGSEQ or raises the existing
%      singleton*.
%
%      H = EXPLOREIMGSEQ returns the handle to a new EXPLOREIMGSEQ or the handle to
%      the existing singleton*.
%
%      EXPLOREIMGSEQ('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EXPLOREIMGSEQ.M with the given input arguments.
%
%      EXPLOREIMGSEQ('Property','Value',...) creates a new EXPLOREIMGSEQ or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before exploreImgSeq_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to exploreImgSeq_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help exploreImgSeq

% Last Modified by GUIDE v2.5 07-Jul-2016 11:43:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @exploreImgSeq_OpeningFcn, ...
                   'gui_OutputFcn',  @exploreImgSeq_OutputFcn, ...
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


% --- Executes just before exploreImgSeq is made visible.
function exploreImgSeq_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to exploreImgSeq (see VARARGIN)

% Choose default command line output for exploreImgSeq
handles.output = hObject;

handles.ei = varargin{1};
handles = loadData(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes exploreImgSeq wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function handles = loadData(handles)
fileName = makeName('ImgSeq',handles.ei.pathName);
load(fileName);
handles.ImgSeq = ImgSeq;
totalFrames = size(ImgSeq,3);
handles.ei.totalFrames = totalFrames;
set(handles.slider_frames,'Max',totalFrames);
set(handles.slider_frames,'Min',1);
minStep = 1/totalFrames;
maxStep = 10/totalFrames;
set(handles.slider_frames,'SliderStep',[minStep maxStep]);
set(handles.slider_frames,'value',1);
metaData.playFlag = 0;
set(handles.pushbutton_play,'userData',metaData);
set(handles.pushbutton_stop,'visible','off');
addlistener(handles.slider_frames,'ContinuousValueChange',@slider_frames_Callback1);
loadFrame(handles,1);


function loadFrame(handles,frameNumber)
if frameNumber == 0 % for loading average image
end
if frameNumber > 0
    thisFrame = handles.ImgSeq(:,:,frameNumber);
end
% set(handles.figure1,'currentAxes',handles.axes_main);
imagesc(thisFrame,'Parent',handles.axes_main);
% axis off; axis equal;
colormap jet;
set(handles.text_frameNumber,'String',sprintf('%d of %d',frameNumber,handles.ei.totalFrames));
set(handles.slider_frames,'value',frameNumber);


% --- Outputs from this function are returned to the command line.
function varargout = exploreImgSeq_OutputFcn(hObject, eventdata, handles) 
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

% --- Executes on slider movement.
function slider_frames_Callback1(hObject, eventdata, handles)
% hObject    handle to slider_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% handles = guihandles;
% hObject = handles.slider_frames;
frameNumber = round(get(hObject,'Value'));
% set(hObject,'Value',frameNumber);
handles = guidata(gcbo);
if frameNumber > 0
    thisFrame = handles.ImgSeq(:,:,frameNumber);
end
imagesc(thisFrame,'Parent',handles.axes_main);
colormap jet;
set(handles.text_frameNumber,'String',sprintf('%d of %d',frameNumber,handles.ei.totalFrames));
set(handles.slider_frames,'value',frameNumber);


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
