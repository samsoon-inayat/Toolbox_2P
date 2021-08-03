function varargout = identifyROIs(varargin)
%IDENTIFYROIS M-file for identifyROIs.fig
%      IDENTIFYROIS, by itself, creates a new IDENTIFYROIS or raises the existing
%      singleton*.
%
%      H = IDENTIFYROIS returns the handle to a new IDENTIFYROIS or the handle to
%      the existing singleton*.
%
%      IDENTIFYROIS('Property','Value',...) creates a new IDENTIFYROIS using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to identifyROIs_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      IDENTIFYROIS('CALLBACK') and IDENTIFYROIS('CALLBACK',hObject,...) call the
%      local function named CALLBACK in IDENTIFYROIS.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help identifyROIs

% Last Modified by GUIDE v2.5 29-Sep-2014 12:24:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
myGUI_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  myGUI_Singleton, ...
                   'gui_OpeningFcn', @identifyROIs_OpeningFcn, ...
                   'gui_OutputFcn',  @identifyROIs_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before identifyROIs is made visible.
function identifyROIs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for identifyROIs
handles.expInfo = varargin{1};
handles = loadData(handles);
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes identifyROIs wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- My Function to load Data
function handles = loadData(handles)
eI = handles.expInfo;
rawFile = uigetfile(sprintf('%s/*.raw',eI.pDataFolder));
rawFile = makeName(rawFile,eI.pDataFolder);
handles.frameSize = [eI.pixelX eI.pixelY];
handles.ROIsFile = eI.ROIsFile;
if exist(eI.ROIsFile);
    load(eI.ROIsFile);
else
    ROIs = [];
end
handles.ROIs = ROIs;
tempString = cellstr(num2str((1:length(ROIs))'));
set(handles.listbox_ROIs,'String',tempString);

handles.totalFrames = eI.totalFrames;
set(handles.slider_recordingPlayer,'Max',handles.totalFrames);
set(handles.slider_recordingPlayer,'Min',1);
minStep = 1/handles.totalFrames;
maxStep = 10/handles.totalFrames;
set(handles.slider_recordingPlayer,'SliderStep',[minStep maxStep]);
set(handles.slider_recordingPlayer,'value',1);
set(handles.playingFlag,'value',0);
set(handles.pushbutton_stop,'visible','off');
set(handles.slider_ROIVisibilityValue,'Value',0.5);
% 
loadFrame(handles,0,0);

% --- My function to load Frame
function loadFrame (handles,listingNumber,frameNumber,varargin)
if frameNumber == 0
    load(makeName('averageImageMC.mat',handles.expInfo.pDataFolder));
    set(handles.figure1,'currentAxes',handles.axesImage);
    imagesc(averageImage);
    axis off; axis equal;
    colormap gray;
    return;
end
channelNumber = str2num(get(handles.edit_channelNumber,'String'));
[rn tn] = findNumbers(handles,listingNumber);
tsn = handles.animal.recording{rn}.TSeriesFolder{tn}.name;
tsnn = handles.animal.recording{rn}.TSeriesFolder{tn}.eye;
if frameNumber == -1000
    roiNumber = get(handles.listbox_ROIs,'Value');
    if length(roiNumber) > 1
        roiNumber = roiNumber(1);
        set(handles.listbox_ROIs,'Value',roiNumber);
    end
%     iNumber = handles.ROIs{roiNumber}.masterInstance;
    ttsn = handles.animal.recording{rn}.TSeriesFolder{tn};
    iNumber = getInstanceNumberFromROI(handles.ROIs{roiNumber},ttsn);
    frameNumber = handles.ROIs{roiNumber}.instances{iNumber}.pixelFrame;
    rn = handles.ROIs{roiNumber}.instances{iNumber}.recordingNumber;
    tsn = handles.ROIs{roiNumber}.instances{iNumber}.folderName;
    if iscellstr(tsn)
        tsn = tsn{1};
    end
    recordingType = ttsn.eye;
%     [recordingType TSN] = getRecordingType(handles.animal,rn,tsn);
%     set(handles.listbox_allTSeries,'Value',sub2ind([3 length(handles.recordingNumbers)],TSN,rn));
end
frameImage = getFrame(handles.animal,rn,tsn,frameNumber,channelNumber);
if frameImage == 0
    return;
end
if frameNumber>=1
    [conditionNumber,stimType,timeOfFrame] = getStimType(handles.animal,rn,tsn,frameNumber);
else
    conditionNumber = -1;
    stimType = -1;
end
if frameNumber ~= -10001
    minagc = min(min(frameImage));
    agc = frameImage - minagc;
    agc = agc/max(max(agc));
    averageImage = agc;
    if get(handles.checkbox_enhanceContrast,'value')
        averageImage = imadjust(averageImage);
    end
    fileName = handles.animal.recording{rn}.TSeriesFrameShiftFile{tn};
    if nameExists(fileName)
        load(fileName);
        if shiftValues.rotation ~=0
            averageImage = imrotate(averageImage,shiftValues.rotation,'crop');
        end
        if shiftValues.xShift ~=0 || shiftValues.yShift ~=0 
            averageImage = circshift(averageImage,[shiftValues.yShift shiftValues.xShift]);
        end
        set(handles.text_frameXShift,'String',sprintf('X Shift = %d',shiftValues.xShift));
        set(handles.text_frameYShift,'String',sprintf('Y Shift = %d',shiftValues.yShift));
        set(handles.text_frameRotation,'String',sprintf('Rotation = %.1f',shiftValues.rotation));
    else
        set(handles.text_frameXShift,'String',sprintf('X Shift = %d',0));
        set(handles.text_frameYShift,'String',sprintf('Y Shift = %d',0));
        set(handles.text_frameRotation,'String',sprintf('Rotation = %.1f',0));
    end


    if get(handles.radiobutton_imageWithIdentifiedROIs,'Value')
        ROIsFrame = makeROIsFrame(handles,rn,tsn);
        nPixels = find(ROIsFrame);
        aIntensity = mean(frameImage(nPixels));
        set(handles.text_averageIntensity,'String',num2str(aIntensity));
        ROIVisibilityValue = get(handles.slider_ROIVisibilityValue,'Value');
        if size(averageImage,1) == size(ROIsFrame,1)
            averageImage = (255 * 0.9 * averageImage + (255 * ROIVisibilityValue * ROIsFrame));
        end
    end

    minLUTSlider = get(handles.slider_minLUT,'Value');
else
    averageImage = frameImage;
end
set(handles.figure1,'currentAxes',handles.axesImage);
imagesc(averageImage);
if get(handles.radiobutton_colorMapGray,'value')
    colormap gray;
end
if get(handles.radiobutton_colorMapJet,'value')
    colormap Jet;
end

axis off;
% if get(handles.checkbox_useAverageImagePixels,'Value') && get(handles.radiobutton_imageWithIdentifiedROIs,'Value')
%     displayROINumbers(handles,rn,tsn);
% end
if get(handles.checkbox_displayBoxes,'Value') && nargin < 4
    displayBoxNumbers(handles,rn,tsn);
end
if frameNumber>=1
    set(handles.slider_recordingPlayer,'Value',frameNumber);
    set(handles.textCurrentFrame,'String',sprintf('%d/%d',frameNumber,handles.totalFrames));
    set(handles.textCurrentFrame,'UserData',frameNumber);
end
if frameNumber <=0
    set(handles.textCurrentFrame,'String',sprintf('%d/%d',frameNumber,handles.totalFrames));
    set(handles.textCurrentFrame,'UserData',frameNumber);
end
pnd = handles.animal.recording{rn}.pnd;
th = text(5,5,sprintf('PND%d-%s  (CN,O) = %d , %.1f',pnd,tsnn,conditionNumber,stimType),'color','w','FontSize',12,'FontWeight','bold');
if exist('timeOfFrame')
    th = text(5,15,sprintf('t = %.2f secs',timeOfFrame),'color','w','FontSize',12,'FontWeight','bold');
else
    th = text(5,15,sprintf('t = %s secs','-'),'color','w','FontSize',12,'FontWeight','bold');
end
set(handles.axesImage,'userdata',listingNumber);
set(handles.textCurrentFrame,'UserData',frameNumber);
if frameNumber == -10001
    nothing = 0;
    fImage = getimage(handles.axesImage);
end
% if get(handles.checkbox_displayTimeLine,'Value')
%     subplotNumber = get(handles.checkbox_displayTimeLine,'userdata');
%     figure(80);subplot(3,1,subplotNumber);
%     plot([frameNumber (frameNumber+1)],[0 2000],'g');
% end

% --- Outputs from this function are returned to the command line.
function varargout = identifyROIs_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in buttonDone.
function buttonDone_Callback(hObject, eventdata, handles)
% hObject    handle to buttonDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.figure1);

% --- Executes on mouse press over axes background.
function axesImage_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.uipanel_ROIbuttonGroup,'visible'),'off')
    return;
else
    pushbutton_selectROIToView_Callback(handles.pushbutton_selectROIToView, eventdata, handles);
end

% --- Executes when selected object is changed in uipanel_imageDisplaySelection.
function uipanel_imageDisplaySelection_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_imageDisplaySelection 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.axesImage,'userdata');
if ~get(handles.playingFlag,'value')
    frameNumber = get(handles.textCurrentFrame,'UserData');
    loadFrame(handles,value,frameNumber);
end
if get(handles.radiobutton_image,'Value')
    set(handles.checkbox_displayROIs,'Value',0);
else
    set(handles.checkbox_displayROIs,'Value',1);
end


% --- Executes on button press in pushbutton_play.
function pushbutton_play_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'visible','off');
set(handles.pushbutton_stop,'visible','on');
set(handles.playingFlag,'value',1);
currentFrame = round(get(handles.slider_recordingPlayer,'Value'));
maxFrame = get(handles.slider_recordingPlayer,'Max');

while get(handles.playingFlag,'value')
    if currentFrame < maxFrame
        currentFrame = currentFrame + 1;
    else
        currentFrame = 1;
        for lii = 1:3
            beep;
        end
        pushbutton_stop_Callback(handles.pushbutton_stop, eventdata, handles)
        break;
    end
    value = get(handles.axesImage,'userdata');
    loadFrame(handles,value,currentFrame);
    delayValue = get(handles.slider_playSpeed,'value');
    pause(delayValue);
end


% --- Executes on slider movement.
function slider_recordingPlayer_Callback(hObject, eventdata, handles)
% hObject    handle to slider_recordingPlayer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameNumber = round(get(hObject,'Value'));
set(hObject,'Value',frameNumber);
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,frameNumber);

% --- Executes during object creation, after setting all properties.
function slider_recordingPlayer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_recordingPlayer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton_stop.
function pushbutton_stop_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.playingFlag,'value',0);
set(hObject,'visible','off');
set(handles.pushbutton_play,'visible','on');
pause(0.3);

% --- Executes on button press in playingFlag.
function playingFlag_Callback(hObject, eventdata, handles)
% hObject    handle to playingFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of playingFlag


% --- Executes on slider movement.
function slider_numberOfStd_Callback(hObject, eventdata, handles)
% hObject    handle to slider_numberOfStd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
value = get(hObject,'Value');
set(handles.text_numberOfStd,'String',num2str(value));

frameNumber = get(handles.textCurrentFrame,'UserData');
valueAX = get(handles.axesImage,'userdata');
loadFrame(handles,valueAX,frameNumber);
valueRNTN = get(handles.axesImage,'userdata');
[rn tsn] = findNumbers(handles,valueRNTN);
tsnStr = handles.animal.recording{rn}.TSeriesFolder{tsn};
currentFrame = getFrame(handles.animal,rn,tsn,frameNumber);
if ~isfield(handles,'ROIx1')
    return;
end
x1 = handles.ROIx1;
x2 = handles.ROIx2;
y1 = handles.ROIy1;
y2 = handles.ROIy2;

ROI = currentFrame(y1:y2,x1:x2);

ROI = ROI - min(min(ROI));
ROI = ROI/max(max(ROI));
meanROI = mean(mean(ROI(~isnan(ROI))));
stdROI = std(ROI(~isnan(ROI)));
%points of interest
numberOfstd = value;
poi = find(ROI > (meanROI+(numberOfstd*stdROI)));
temped = zeros(size(ROI));
temped(poi) = ones(1,length(poi));
set(handles.figure1,'currentAxes',handles.axesZoomed);
imagesc(temped);
axis off;
handles.temped = temped;
guidata(handles.figure1,handles);


% --- Executes during object creation, after setting all properties.
function slider_numberOfStd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_numberOfStd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton_deletePixels.
function pushbutton_deletePixels_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_deletePixels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'visible','off');
pause(0.5);
ed = handles.temped;
set(handles.figure1,'currentAxes',handles.axesZoomed);
k = 0;
while k == 0
    k = waitforbuttonpress;
    if strcmp(get(handles.figure1,'SelectionType'),'open')
        break;
    end
    point1 = get(gca,'CurrentPoint');    % button down detected
    finalRect = rbbox;                   % return figure units
    point2 = get(gca,'CurrentPoint');    % button up detected
    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    x1 = min(floor(point1(1)),floor(point2(1)));
    x2 = max(floor(point1(1)),floor(point2(1)));
    y1 = min(floor(point1(2)),floor(point2(2)));
    y2 = max(floor(point1(2)),floor(point2(2)));

    if x1 < 1
        x1 = 1;
    end
    if y1 < 1
        y1 = 1;
    end
    if x2 > size(ed,2)
        x2 = size(ed,2);
    end
    if y2 > size(ed,1)
        y2 = size(ed,1);
    end
    ed(y1:y2,x1:x2) = zeros(length(y1:y2),length(x1:x2));
    handles.temped = ed;
    guidata(gcbo,handles);
    %set(handles.figure1,'currentAxes',handles.axisOverLayed);
    set(handles.figure1,'currentAxes',handles.axesZoomed);
    imagesc(handles.temped);
    axis off;
end
set(hObject,'visible','on');


% --- Executes on button press in pushbutton_storePixels.
function pushbutton_storePixels_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_storePixels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

thisBox = [handles.ROIx1 handles.ROIx2 handles.ROIy1 handles.ROIy2];
if isequal(handles.lastBox,thisBox)
    return;
end

ROI_selectedList = get(handles.listbox_ROIs,'Value');
if length(ROI_selectedList) == length(handles.ROIs)
    allWereSelected = 1;
else
    allWereSelected = 0;
end
frameNumber = get(handles.textCurrentFrame,'UserData');
x1 = handles.ROIx1;
x2 = handles.ROIx2;
y1 = handles.ROIy1;
y2 = handles.ROIy2;
handles.lastBox = [handles.ROIx1 handles.ROIx2 handles.ROIy1 handles.ROIy2];
numberOfPointsX = x2-x1+1;
numberOfPointsY = y2-y1+1;
% currentPixels = find(handles.temped);
% [row col] = ind2sub([numberOfPointsY numberOfPointsX],currentPixels);rowG = y1+row-1;colG = x1+col-1;
% ROIpixels = sub2ind(handles.frameSize,rowG,colG);

% existingROIs = findExistingROIs(ROIpixels,handles);
% if isempty(existingROIs)
%     newROI = 1;
% else
%     newROI = 0;
% end

[rn tn] = findNumbers(handles,get(handles.listbox_allTSeries,'Value'));
folderName = handles.animal.recording{rn}.TSeriesFolder{tn}.name;
if eventdata(1) >= 2
    newROI = 1;
else
    newROI = 0;
end
if newROI == 0 && get(handles.checkbox_makeSubBox,'Value') == 1
    smallBox = [x1 y1 x2 y2];
    [roiNumber instanceNumber] = getROINumber(folderName,handles.ROIs,smallBox);
    ROI = handles.ROIs{roiNumber};
    ROI.instances{instanceNumber}.smallBox = [x1 y1 x2 y2];
    handles.ROIs{roiNumber} = ROI;
    guidata(gcbo,handles);
    saveROIs(handles);
    return;
end
roiList = get(handles.listbox_ROIs,'String');
if length(roiList) == 1 && length(roiList{1}) == 0
    firstROI = 1;
else
    firstROI = 0;
    roiListNumbers = get(handles.listbox_ROIs,'Value');
    if length(roiListNumbers) > 1 && newROI == 0
        display('Choose one ROI');
        return;
    end
end

if eventdata(1) == 1 || eventdata(1) == 3
    ROITypeString = 'Soma';
elseif eventdata(1) == 4
    ROITypeString = 'Blood Vessel';
else
    ROITypeString = selectROIType;   
end

if strcmp(ROITypeString,'')
    return;
end

if firstROI || newROI
    ROI = [];
    instanceNumber = 1;
    ROI.masterInstance = 1;
    ROI.id = datenum(clock);
end

if firstROI
    roiNumber = 1;
end

ROIs = handles.ROIs;
if newROI
    if firstROI
        roiNumber = 1;
    else
        roiNumber = length(roiList) + 1;
    end
else
    roiNumber = roiListNumbers;
    ROI = ROIs{roiNumber};
    instanceNumber = getInstanceNumberFromROI(ROI,folderName);
    if instanceNumber == 0
        instanceNumber = length(ROI.instances) + 1;
    end
end
   
ROI.instances{instanceNumber}.recordingNumber = rn;
% ROI.instances{instanceNumber}.pixels = ROIpixels;
ROI.instances{instanceNumber}.box = [x1 y1 x2 y2];
ROI.instances{instanceNumber}.roiPoly = handles.roiPoly;
ROI.instances{instanceNumber}.folderName = folderName;
ROI.instances{instanceNumber}.identificationFrame = frameNumber;
ROI.type = ROITypeString;
ROIs{roiNumber} = ROI;
handles.ROIs = ROIs;
guidata(gcbo,handles);
saveROIs(handles);
uipanel_imageDisplaySelection_SelectionChangeFcn(handles.uipanel_imageDisplaySelection, eventdata, handles);
tempString = cellstr(num2str((1:length(ROIs))'));
set(handles.listbox_ROIs,'String',tempString);
if allWereSelected
    set(handles.listbox_ROIs,'Value',1:length(get(handles.listbox_ROIs,'String')));
else
    set(handles.listbox_ROIs,'Value',roiNumber);
end
listbox_ROIs_Callback(handles.listbox_ROIs, [], handles);
% if ~get(handles.checkbox_displayROIs,'Value')
%     set(handles.checkbox_displayROIs,'Value',1);
%     checkbox_displayROIs_Callback(handles.checkbox_displayROIs,[],handles);
% end
% if strcmp(lower(ROITypeString),'blood vessel');
%     fileName = makeName('bloodVesselROIs.mat',handles.animal.recording{rn}.TSeriesFolder{tn}.processedDataFolderPath);
%     if nameExists(fileName)
%         load(fileName);
%         bloodVesselROIs = [bloodVesselROIs roiNumber];
%     else
%         bloodVesselROIs = roiNumber;
%     end
%     save(fileName,'bloodVesselROIs');
% end

% --- Executes on slider movement.
function slider_playSpeed_Callback(hObject, eventdata, handles)
% hObject    handle to slider_playSpeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_playSpeed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_playSpeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton_loadAverageImage.
function pushbutton_loadAverageImage_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadAverageImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,0);

% --- Executes on button press in pushbutton_sliderLeft.
function pushbutton_sliderLeft_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_sliderLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
frameNumber = round(get(handles.slider_recordingPlayer,'Value'));
frameNumber = frameNumber - 1;
if frameNumber < 1
    frameNumber = 1;
    beep;
end
set(handles.slider_recordingPlayer,'Value',frameNumber);
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,frameNumber);

% --- Executes on button press in pushbutton_sliderRight.
function pushbutton_sliderRight_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_sliderRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
frameNumber = round(get(handles.slider_recordingPlayer,'Value'));
frameNumber = frameNumber + 1;
if frameNumber > handles.totalFrames
    frameNumber = handles.totalFrames;
    beep;
end
set(handles.slider_recordingPlayer,'Value',frameNumber);
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,frameNumber);

% --- Executes on button press in pushbutton_loadMaxImage.
function pushbutton_loadMaxImage_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadMaxImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,-1);

% --- Executes on button press in pushbutton_storeSomaPixels.
function pushbutton_storeSomaPixels_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_storeSomaPixels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pushbutton_storePixels_Callback(handles.pushbutton_storePixels, 1, handles);
value = get(handles.axesImage,'userdata');
handles = guidata(handles.figure1);
loadFrame(handles,value,0);

% --- Executes on button press in pushbutton_deleteROIInstance.
function pushbutton_deleteROIInstance_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_deleteROIInstance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[rn tsnv] = findNumbers(handles,get(handles.listbox_allTSeries,'Value'));
tsns = handles.animal.recording{rn}.TSeriesFolder{tsnv};
folderName = tsns.name;
roiListNumbers = get(handles.listbox_ROIs,'Value');
% if length(roiListNumbers) > 1
%     return;
% end
ROIs = handles.ROIs;
for ii = 1:length(roiListNumbers)
    ROI = ROIs{roiListNumbers(ii)};
    instanceNumber = getInstanceNumberFromROI(ROI,folderName);
    if instanceNumber == 0
        continue;
    end
    if ROI.masterInstance == instanceNumber
        display('Cannot delete master instance');
        continue;
    end
    ROI.instances(instanceNumber) = [];
    ROIs{roiListNumbers(ii)} = ROI;
end
handles.ROIs = ROIs;
guidata(gcbo,handles);
saveROIs(handles);
uipanel_imageDisplaySelection_SelectionChangeFcn(handles.uipanel_imageDisplaySelection, eventdata, handles);
tempString = cellstr(num2str((1:length(ROIs))'));
if length(ROIs) == 0
    set(handles.listbox_ROIs,'String','');
else
    set(handles.listbox_ROIs,'String',tempString);
end


% --- Executes on button press in pushbutton_selectROIToView.
function pushbutton_selectROIToView_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_selectROIToView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   
% set(handles.uipanel_ROIbuttonGroup,'visible','off');

set(handles.figure1,'currentAxes',handles.axesImage);
h = impoly;
position = wait(h);
position = floor(position);
h.delete;
handles.roiPoly = position;
handles.ROIx1 = min(position(:,1));
handles.ROIx2 = max(position(:,1));
handles.ROIy1 = min(position(:,2));
handles.ROIy2 = max(position(:,2));
guidata(hObject, handles);
slider_numberOfStd_Callback(handles.slider_numberOfStd, eventdata, handles);
if isempty(eventdata)
    pushbutton_storePixels_Callback(handles.pushbutton_storePixels, 2, handles);
elseif eventdata == 1
    pushbutton_storePixels_Callback(handles.pushbutton_storePixels, [1 1], handles);
end

% set(handles.uipanel_ROIbuttonGroup,'visible','on');

% BW = poly2mask(position(:,1),position(:,2),handles.frameSize(1),handles.frameSize(2));

% k = 0;

% k = 0;
% while k == 0
%     k = waitforbuttonpress;
%     point1 = get(gca,'CurrentPoint');    % button down detected
%     finalRect = rbbox;                   % return figure units
%     point2 = get(gca,'CurrentPoint');    % button up detected
%     point1 = point1(1,1:2);point2 = point2(1,1:2);              % extract x and y
%     x1 = min(floor(point1(1)),floor(point2(1)));x2 = max(floor(point1(1)),floor(point2(1)));
%     y1 = min(floor(point1(2)),floor(point2(2)));y2 = max(floor(point1(2)),floor(point2(2)));
%     if x1 < 1
%         x1 = 1;
%     end
%     if y1 < 1
%         y1 = 1;
%     end
%     if x2 > handles.frameSize(2)
%         x2 = handles.frameSize(2);
%     end
%     if y2 > handles.frameSize(1)
%         y2 = handles.frameSize(1);
%     end
%     numberOfPointsX = x2-x1+1;
%     numberOfPointsY = y2-y1+1;
%     totalPoints = numberOfPointsX * numberOfPointsY;
%     k = 1;
% end
% handles.ROIx1 = x1;
% handles.ROIx2 = x2;
% handles.ROIy1 = y1;
% handles.ROIy2 = y2;
% guidata(hObject, handles);
% set(handles.uipanel_ROIbuttonGroup,'visible','on');
% slider_numberOfStd_Callback(handles.slider_numberOfStd, eventdata, handles);

% --- Executes on button press in checkbox_useAverageImagePixels.
function checkbox_useAverageImagePixels_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_useAverageImagePixels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_useAverageImagePixels
uipanel_imageDisplaySelection_SelectionChangeFcn(handles.uipanel_imageDisplaySelection, eventdata, handles);


% --- Executes on slider movement.
function slider_ROIVisibilityValue_Callback(hObject, eventdata, handles)
% hObject    handle to slider_ROIVisibilityValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
uipanel_imageDisplaySelection_SelectionChangeFcn(handles.uipanel_imageDisplaySelection, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function slider_ROIVisibilityValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_ROIVisibilityValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton_moveLeft.
function pushbutton_moveLeft_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_moveLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
moveROIs(handles,'L');

% --- Executes on button press in pushbutton_moveRight.
function pushbutton_moveRight_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_moveRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
moveROIs(handles,'R');

% --- Executes on button press in pushbutton_moveUp.
function pushbutton_moveUp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_moveUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
moveROIs(handles,'U');

% --- Executes on button press in pushbutton_moveDown.
function pushbutton_moveDown_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_moveDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
moveROIs(handles,'D');

% --- Executes on button press in pushbutton_rotateCW.
function pushbutton_rotateCW_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_rotateCW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
moveROIs(handles,'CW');

% --- Executes on button press in pushbutton_rotateCCW.
function pushbutton_rotateCCW_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_rotateCCW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
moveROIs(handles,'CCW');

% --- My function
function moveROIs(handles,movementParameter)

if get(handles.checkbox_makeSubBox,'value')
    smallBox = 1;
else
    smallBox = 0;
end
ROIsToDisplace = get(handles.listbox_ROIs,'Value');
[rn tn] = findNumbers(handles,get(handles.axesImage,'userdata'));
tsn = handles.animal.recording{rn}.TSeriesFolder{tn};
rotationAngle = 5;
ROIs = handles.ROIs;
if strcmp(lower(movementParameter),'expand') | strcmp(lower(movementParameter),'shrink') ...
        | strcmp(lower(movementParameter),'cw') | strcmp(lower(movementParameter),'ccw')
    if length(ROIsToDisplace) > 1
        ROIsToDisplace = ROIsToDisplace(end);
    end
end
for ii = 1:length(ROIsToDisplace)
    roiNumber = ROIsToDisplace(ii);
    ROI = ROIs{roiNumber};
    instanceNumber = getInstanceNumberFromROI(ROI,tsn);
    if instanceNumber == 0
        continue;
    end
%     pixels = ROI.instances{instanceNumber}.pixels;
    pixels = 100:300;
    if smallBox
        box = ROI.instances{instanceNumber}.smallBox;
    else
        box = ROI.instances{instanceNumber}.box;
    end
    roiPoly = ROI.instances{instanceNumber}.roiPoly;
    
    [I,J] = ind2sub(handles.frameSize,pixels); 
    minI = min(I); maxI = max(I); 
    minJ = min(J); maxJ = max(J);
    if strcmp(lower(movementParameter),'l')
        J = J - 1;
        box(1) = box(1) - 1;box(3) = box(3) - 1;
        if box(1) == 0
            box(1) == 1;
        end
        roiPoly(:,1) = roiPoly(:,1) - 1;
    end
    if strcmp(lower(movementParameter),'r')
        J = J + 1;
        box(1) = box(1) + 1;box(3) = box(3) + 1;
        if box(3) > handles.frameSize(2)
            box(3) = handles.frameSize(2);
        end
        roiPoly(:,1) = roiPoly(:,1) + 1;
    end
    if strcmp(lower(movementParameter),'u')
        I = I - 1;
        box(2) = box(2) - 1;box(4) = box(4) - 1;
        if box(2) == 0
            box(2) = 1;
        end
        roiPoly(:,2) = roiPoly(:,2) - 1;
    end
    if strcmp(lower(movementParameter),'d')
        I = I + 1;
        box(2) = box(2) + 1;box(4) = box(4) + 1;
        if box(4) > handles.frameSize(1)
            box(4) = handles.frameSize(1);
        end
        roiPoly(:,2) = roiPoly(:,2) + 1;
    end
   
    if strcmp(lower(movementParameter),'cw')
        [ geom, iner, cpmo ] = polygeom( roiPoly(:,1), roiPoly(:,2) );
        x_cent = floor(geom(2)); y_cent = floor(geom(3));
        BW = poly2mask(roiPoly(:,1),roiPoly(:,2),handles.frameSize(1),handles.frameSize(2));
        imx_center = floor(handles.frameSize(2)/2);
        imy_center = floor(handles.frameSize(1)/2);
        xShift = x_cent - imx_center;
        yShift = y_cent - imy_center;
        sBW = circshift(BW,[-yShift -xShift]);
        rsBW = imrotate(sBW,-rotationAngle,'crop');
        srsBW = circshift(rsBW,[yShift xShift]);
        bnd = bwboundaries(srsBW,'noholes');
        roiPoly = [];
        roiPoly(:,1) = bnd{1}(:,2);
        roiPoly(:,2) = bnd{1}(:,1);
        
%         subRegion = zeros(maxI-minI+1,maxJ-minJ+1);
%         sPixels = sub2ind(size(subRegion),I-minI+1,J-minJ+1);
%         subRegion(sPixels) = ones(size(sPixels));
%         rotatedSubRegion = imrotate(subRegion,-rotationAngle);
%         rotatedSubRegion(find(rotatedSubRegion)) = ceil(rotatedSubRegion(find(rotatedSubRegion)));
%         [ro,co] = size(subRegion);
%         [rr,cr] = size(rotatedSubRegion);
%         dr = rr - ro; dc = cr - co;
%         newMinI = minI - round(dr/2);
%         newMinJ = minJ - round(dc/2);
%         currentPixels = find(rotatedSubRegion);
%         [row col] = ind2sub(size(rotatedSubRegion),currentPixels); 
%         I = newMinI+row-1;
%         J = newMinJ+col-1;
    end
    if strcmp(lower(movementParameter),'ccw')
        [ geom, iner, cpmo ] = polygeom( roiPoly(:,1), roiPoly(:,2) );
        x_cent = floor(geom(2)); y_cent = floor(geom(3));
        BW = poly2mask(roiPoly(:,1),roiPoly(:,2),handles.frameSize(1),handles.frameSize(2));
        imx_center = floor(handles.frameSize(2)/2);
        imy_center = floor(handles.frameSize(1)/2);
        xShift = x_cent - imx_center;
        yShift = y_cent - imy_center;
        sBW = circshift(BW,[-yShift -xShift]);
        rsBW = imrotate(sBW,rotationAngle,'crop');
        srsBW = circshift(rsBW,[yShift xShift]);
        bnd = bwboundaries(srsBW,'noholes');
        roiPoly = [];
        roiPoly(:,1) = bnd{1}(:,2);
        roiPoly(:,2) = bnd{1}(:,1);
%         subRegion = zeros(maxI-minI+1,maxJ-minJ+1);
%         sPixels = sub2ind(size(subRegion),I-minI+1,J-minJ+1);
%         subRegion(sPixels) = ones(size(sPixels));
%         rotatedSubRegion = imrotate(subRegion,rotationAngle);
%         rotatedSubRegion(find(rotatedSubRegion)) = ceil(rotatedSubRegion(find(rotatedSubRegion)));
%         [ro,co] = size(subRegion);
%         [rr,cr] = size(rotatedSubRegion);
%         dr = rr - ro; dc = cr - co;
%         newMinI = minI - round(dr/2);
%         newMinJ = minJ - round(dc/2);
%         currentPixels = find(rotatedSubRegion);
%         [row col] = ind2sub(size(rotatedSubRegion),currentPixels); 
%         I = newMinI+row-1;
%         J = newMinJ+col-1;
    end
    if strcmp(lower(movementParameter),'shrink')
        [ geom, iner, cpmo ] = polygeom( roiPoly(:,1), roiPoly(:,2) );
        x_cent = floor(geom(2)); y_cent = floor(geom(3));
        ind = roiPoly(:,1) < x_cent; roiPoly(ind,1) = roiPoly(ind,1) + 1;
        ind = roiPoly(:,1) >= x_cent; roiPoly(ind,1) = roiPoly(ind,1) - 1;
        ind = roiPoly(:,2) < y_cent; roiPoly(ind,2) = roiPoly(ind,2) + 1;
        ind = roiPoly(:,2) >= y_cent; roiPoly(ind,2) = roiPoly(ind,2) - 1;
%         BW = poly2mask(roiPoly(:,1),roiPoly(:,2),handles.frameSize(1),handles.frameSize(2));
%         bnd = bwboundaries(BW,'noholes');
%         roiPoly = [];
%         roiPoly(:,1) = bnd{1}(:,2);
%         roiPoly(:,2) = bnd{1}(:,1);
    end
%     ROI.instances{instanceNumber}.pixels = sub2ind(handles.frameSize,I,J);
    if strcmp(lower(movementParameter),'expand')
        [ geom, iner, cpmo ] = polygeom( roiPoly(:,1), roiPoly(:,2) );
        x_cent = floor(geom(2)); y_cent = floor(geom(3));
        ind = roiPoly(:,1) < x_cent; roiPoly(ind,1) = roiPoly(ind,1) - 1;
        ind = roiPoly(:,1) >= x_cent; roiPoly(ind,1) = roiPoly(ind,1) + 1;
        ind = roiPoly(:,2) < y_cent; roiPoly(ind,2) = roiPoly(ind,2) - 1;
        ind = roiPoly(:,2) >= y_cent; roiPoly(ind,2) = roiPoly(ind,2) + 1;
        BW = poly2mask(roiPoly(:,1),roiPoly(:,2),handles.frameSize(1),handles.frameSize(2));
%         bnd = bwboundaries(BW,'noholes');
%         roiPoly = [];
%         roiPoly(:,1) = bnd{1}(:,2);
%         roiPoly(:,2) = bnd{1}(:,1);
    end
    if smallBox
        ROI.instances{instanceNumber}.smallBox = box;
    else
        ROI.instances{instanceNumber}.box = box;
    end
    ROI.instances{instanceNumber}.roiPoly = roiPoly;
    theBox = [min(roiPoly(:,1)) min(roiPoly(:,2)) max(roiPoly(:,1)) max(roiPoly(:,2))];
    ROI.instances{instanceNumber}.box = theBox;
    ROIs{roiNumber} = ROI;
end
handles.ROIs = ROIs;
guidata(gcbo,handles);
saveROIs(handles);
value = get(handles.axesImage,'userdata');
frameNumber = get(handles.textCurrentFrame,'UserData');
loadFrame(handles,value,frameNumber);
uicontrol(handles.listbox_allTSeries);
% uipanel_imageDisplaySelection_SelectionChangeFcn(handles.uipanel_imageDisplaySelection, [], handles);

% --- My Function to make ROI Frame
function ROIFrame = makeROIsFrame(handles,rn,tsn)
ROIFrame = zeros(handles.frameSize);
if strfind(lower(tsn),'tser')
    processedDataFolder = makeName(tsn,handles.animal.recording{rn}.processedDataFolderPath);
else
    commandText = sprintf('processedDataFolder = handles.animal.recording{rn}.processedDataFolderPath%s;',tsn);
    eval(commandText);
end
stdValue = get(handles.slider_numberOfStd,'Value');
try
    load(makeName('pnROIPixels.mat',processedDataFolder));
    load(makeName('pROIPixels.mat',processedDataFolder));
    load(makeName('apROIPixels.mat',processedDataFolder));
catch
    nothing = 0;
end
indexOfStd = (stdValue/0.25)+1;
ROIs = handles.ROIs;
if isempty(ROIs)
    return;
end
ROINumbersToDisplay = get(handles.listbox_ROIs,'Value');
ROIFrame = zeros(handles.frameSize);
if isempty(strfind(tsn,'TSer'))
    folderName = eval(sprintf('handles.animal.recording{rn}.%sFolderName;',lower(tsn)));
else
    folderName = tsn;
end
for ii = ROINumbersToDisplay
    oneROI = ROIs{ii};
    instanceFound = 0;
    for jj = 1:length(oneROI.instances)
        if isempty(oneROI.instances{jj})
            continue;
        end
        if strcmp(oneROI.instances{jj}.folderName,folderName)
            instanceFound = 1;
            break;
        end
    end
    if instanceFound
        try
%             box = oneROI.instances{jj}.box;
%             sROIFrame = ROIPixelFrame(box(2):box(4),box(1):box(3));
            if get(handles.checkbox_displayNeuropil,'Value')
                neuropilPixels = nROIPixels{ii}{indexOfStd};
            else
                neuropilPixels = [];
            end
            if get(handles.checkbox_displayROIs,'Value')
                if get(handles.checkbox_useAverageImagePixels,'Value')
                    roiPixels = aROIPixels{ii}{indexOfStd};
                else
                    roiPixels = ROIPixels{ii}{indexOfStd};
                end
            else
                roiPixels = [];
            end
            pixels = union(roiPixels,neuropilPixels);
        catch
            pixels = 100:300;
        end
    else
        continue;
    end
%     ROIFrame(box(2):box(4),box(1):box(3)) = sROIFrame;
    ROIFrame(pixels) = ones(size(pixels));
end


% --- My Function To find existing ROIs
function existingROIs = findExistingROIs(ROIpixels,handles)
existingROIIDs = [];
instanceNumbers = [];
ROIs = handles.ROIs;
try
    for hh = 1:length(ROIs)
        oneROI = ROIs{hh};
        for ii = 1:length(oneROI.instances)
            instance = oneROI.instances{ii};
            if isempty(instance)
                continue;
            else
                pixels = instance.pixels;
            end
            for jj = 1:length(ROIpixels)
                indices = find(pixels == ROIpixels(jj));
                if ~isempty(indices)
                    existingROIIDs = [existingROIIDs hh];
                    instanceNumbers = [instanceNumbers ii];
                    break;
                end
            end
        end
    end
    existingROIs = [existingROIIDs' instanceNumbers'];
    existingROIs = unique(existingROIs,'rows');
catch
    existingROIs = [existingROIIDs' instanceNumbers'];
    existingROIs = unique(existingROIs,'rows');
    return;
end

% --- My Function
function displayBoxNumbers(handles,rn,tsn)
ROIs = handles.ROIs;
ROINumbersToDisplay = get(handles.listbox_ROIs,'Value');
ROIFrame = zeros(handles.frameSize);
if isempty(strfind(tsn,'TSer'))
    folderName = eval(sprintf('handles.animal.recording{rn}.%sFolderName;',lower(tsn)));
else
    folderName = tsn;
end
for ii = ROINumbersToDisplay
    try
        oneROI = ROIs{ii};
    catch
        break;
    end
    instanceFound = 0;
    for jj = 1:length(oneROI.instances)
        if isempty(oneROI.instances{jj})
            continue;
        end
        if strcmp(oneROI.instances{jj}.folderName,folderName)
            instanceFound = 1;
            break;
        end
    end
    if instanceFound
            box = oneROI.instances{jj}.box;
%             box = oneROI.instances{jj}.box;
            if isfield(oneROI.instances{jj},'smallBox');
                smallBox = oneROI.instances{jj}.smallBox;
            else
                smallBox = [];
            end
            if isfield(oneROI.instances{jj},'roiPoly')
                set(handles.figure1,'currentAxes',handles.axesImage);
                line(oneROI.instances{jj}.roiPoly(:,1),oneROI.instances{jj}.roiPoly(:,2));
                lastOne = size(oneROI.instances{jj}.roiPoly,1);
                line(oneROI.instances{jj}.roiPoly([1 lastOne],1),oneROI.instances{jj}.roiPoly([1 lastOne],2));
                displayBox1(handles,oneROI.instances{jj}.box,ii);
            end
            if isfield(oneROI.instances{jj},'box') & ~isfield(oneROI.instances{jj},'roiPoly')
                displayBox(handles,oneROI.instances{jj}.box,ii);
            end
    else
        continue;
    end
%     displayBox(handles,box,ii);
    if ~isempty(smallBox)
        displayBox(handles,smallBox,ii);
    end
end

function rh = displayBox (handles,box,number)
minI = box(2); maxI = box(4); 
minJ = box(1); maxJ = box(3);
set(handles.figure1,'currentAxes',handles.axesImage);
if number == 0
    rh = rectangle('Position',[minJ,minI,maxJ-minJ,maxI-minI],'Curvature',0.1,...
     'EdgeColor','k','LineWidth',1,'LineStyle','--');
else
    rh = rectangle('Position',[minJ,minI,maxJ-minJ,maxI-minI],'Curvature',0.1,...
         'EdgeColor','k','LineWidth',1,'LineStyle','--');
end
if number > 0
    th = text(minJ-1,minI-1,num2str(number),'color','r','FontSize',10,'FontWeight','bold');
end

function rh = displayBox1 (handles,box,number)
minI = box(2); maxI = box(4); 
minJ = box(1); maxJ = box(3);
set(handles.figure1,'currentAxes',handles.axesImage);
if number > 0
    th = text(minJ-1,minI-1,num2str(number),'color','m','FontSize',8,'FontWeight','bold');
end

% --- My Function
function displayROINumbers(handles,rn,tsn)
ROIs = handles.ROIs;
ROINumbersToDisplay = get(handles.listbox_ROIs,'Value');
ROIFrame = zeros(handles.frameSize);
if isempty(strfind(tsn,'TSer'))
    folderName = eval(sprintf('handles.animal.recording{rn}.%sFolderName;',lower(tsn)));
else
    folderName = tsn;
end
for ii = ROINumbersToDisplay
    oneROI = ROIs{ii};
    instanceFound = 0;
    for jj = 1:length(oneROI.instances)
        if isempty(oneROI.instances{jj})
            continue;
        end
        if strcmp(oneROI.instances{jj}.folderName,folderName)
            instanceFound = 1;
            break;
        end
    end
    if instanceFound
        pixels = oneROI.instances{jj}.pixels;
    else
        continue;
    end
%     displayBoxNumbers(handles,pixels,ii);
end

% function displayBoxNumbers (handles,pixels,number)
% [I,J] = ind2sub(handles.frameSize,pixels); 
% minI = min(I); maxI = max(I); 
% minJ = min(J); maxJ = max(J);
% set(handles.figure1,'currentAxes',handles.axesImage);
% rh = rectangle('Position',[minJ,minI,maxJ-minJ,maxI-minI],'Curvature',0.4,...
%      'EdgeColor','w','LineWidth',1,'LineStyle','--');
% th = text(minJ-1,minI-1,num2str(number),'color','k','FontSize',11,'FontWeight','bold');

% --- Executes on button press in pushbutton_createROIInstance.
function pushbutton_createROIInstance_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_createROIInstance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[rn tsnv] = findNumbers(handles,get(handles.listbox_allTSeries,'Value'));
% tsns = handles.animal.recording{rn}.TSeriesFolderName;
tsns = handles.animal.recording{rn}.TSeriesFolder{tsnv}.name;
roiListNumbers = get(handles.listbox_ROIs,'Value');
ROIs = handles.ROIs;

if get(handles.checkbox_makeSubBox,'value')
    for ii = 1:length(roiListNumbers)
        roiNumber = roiListNumbers(ii);
        ROI = ROIs{roiNumber};
        for jj = 1:length(ROI.instances)
            if isfield(ROI.instances{jj},'smallBox')
                smallBox = ROI.instances{jj}.smallBox;
                break;
            end
        end
        if jj == length(ROI.instances) && ~isfield(ROI.instances{jj},'smallBox')
            continue;
        end
        instanceWithSmallBox = jj;
        for jj = 1:length(ROI.instances)
            if ~isfield(ROI.instances{jj},'smallBox')
                bigBoxofSmallBox = ROI.instances{instanceWithSmallBox}.box;
                thisBigBox = ROI.instances{jj}.box;
                translation = thisBigBox - bigBoxofSmallBox;
                ROI.instances{jj}.smallBox = smallBox + translation;
            end
        end
        ROIs{roiNumber} = ROI;
    end
else
    for ii = 1:length(roiListNumbers)
        roiNumber = roiListNumbers(ii);
        ROI = ROIs{roiNumber};
        if getInstanceNumberFromROI(ROI,tsns) == 0
            masterInstance = ROI.masterInstance;
            instanceNumber = length(ROI.instances)+1;
            ROI.instances{instanceNumber} = ROI.instances{masterInstance};
            ROI.instances{instanceNumber}.recordingNumber = rn;
            ROI.instances{instanceNumber}.folderName = tsns;
            ROI.instances{instanceNumber}.identificationFrame = [];
            ROIs{roiNumber} = ROI;
        end
    end
end
handles.ROIs = ROIs;
guidata(gcbo,handles);
saveROIs(handles);
frameNumber = get(handles.textCurrentFrame,'UserData');
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,frameNumber);

function saveROIs(handles)
ROIs = handles.ROIs;
save(handles.animal.recording{handles.recordingNumbers(1)}.ROIsFilePath,'ROIs');


% --- Executes on selection change in listbox_ROIs.
function listbox_ROIs_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_ROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_ROIs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_ROIs
if get(handles.checkbox_masterInstance,'Value')
    checkbox_masterInstance_Callback(handles.checkbox_masterInstance, eventdata, handles);
    return;
end
if ~get(handles.playingFlag,'value')
    value = get(handles.axesImage,'userdata');
    frameNumber = get(handles.textCurrentFrame,'UserData');
    loadFrame(handles,value,frameNumber);
    if get(handles.checkbox_plot,'value')
        pushbutton_plotResponses_Callback(handles.pushbutton_storeROIIdentities, [], handles)
    end
end

    
% --- Executes during object creation, after setting all properties.
function listbox_ROIs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_ROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_masterInstance.
function checkbox_masterInstance_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_masterInstance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_masterInstance

if ~get(handles.playingFlag,'value') && get(hObject,'Value')
    value = get(handles.axesImage,'userdata');
    frameNumber = -1000;
    loadFrame(handles,value,frameNumber);
end


% --- Executes on button press in pushbutton_nextTSeries.
function pushbutton_nextTSeries_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_nextTSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentValue = get(handles.axesImage,'userdata');
value = get(handles.listbox_allTSeries,'Value');

if length(value) == 1
    value = value + 1;
    if value > (length(handles.recordingNumbers) * 3)
        value = 1;
    end
    frameNumber = get(handles.textCurrentFrame,'userdata');
    loadFrame(handles,value,frameNumber);
    set(handles.listbox_allTSeries,'Value',value);
    return;
end

ii = find(value == currentValue);
if ii == length(value)
    ii = 1;
else
    ii = ii + 1;
end
frameNumber = get(handles.textCurrentFrame,'userdata');
loadFrame(handles,value(ii),frameNumber);

% --- Executes on button press in pushbutton_previousTSeries.
function pushbutton_previousTSeries_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_previousTSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentValue = get(handles.axesImage,'userdata');
value = get(handles.listbox_allTSeries,'Value');

if length(value) == 1
    value = value - 1;
    if value < 1 
        value = length(handles.recordingNumbers) * 3;
    end
    frameNumber = get(handles.textCurrentFrame,'userdata');
    loadFrame(handles,value,frameNumber);
    set(handles.listbox_allTSeries,'Value',value);
    return;
end

ii = find(value == currentValue);
if ii == 1
    ii = length(value);
else
    ii = ii - 1;
end
frameNumber = get(handles.textCurrentFrame,'userdata');
loadFrame(handles,value(ii),frameNumber);


% --- Executes on button press in pushbutton_playTSeries.
function pushbutton_playTSeries_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_playTSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'visible','off');
set(handles.pushbutton_stopTSeries,'visible','on');
set(handles.checkbox_playingFlagTSeries,'value',1);
currentFrame = round(get(handles.slider_recordingPlayer,'Value'));
while get(handles.checkbox_playingFlagTSeries,'value')
    pushbutton_nextTSeries_Callback(handles.pushbutton_nextTSeries, eventdata, handles);
    delayValue = get(handles.slider_playSpeedTSeries,'value');
    pause(delayValue);
end
set(hObject,'visible','on');
set(handles.pushbutton_stopTSeries,'visible','off');

% --- Executes on button press in pushbutton_stopTSeries.
function pushbutton_stopTSeries_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_stopTSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.checkbox_playingFlagTSeries,'value',0);
set(hObject,'visible','off');

% --- Executes on button press in checkbox_playingFlagTSeries.
function checkbox_playingFlagTSeries_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_playingFlagTSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_playingFlagTSeries


% --- Executes on slider movement.
function slider_playSpeedTSeries_Callback(hObject, eventdata, handles)
% hObject    handle to slider_playSpeedTSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_playSpeedTSeries_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_playSpeedTSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in checkbox_displayROIs.
function checkbox_displayROIs_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_displayROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_displayROIs
if get(hObject,'Value')
    set(handles.radiobutton_imageWithIdentifiedROIs,'Value',1);
    set(handles.radiobutton_image,'Value',0);
else
    set(handles.radiobutton_imageWithIdentifiedROIs,'Value',0);
    set(handles.radiobutton_image,'Value',1);
end
uipanel_imageDisplaySelection_SelectionChangeFcn(hObject, eventdata, handles);


% --- Executes on button press in pushbutton_createROIInstanceLocal.
function pushbutton_createROIInstanceLocal_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_createROIInstanceLocal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[rn tn] = findNumbers(handles,get(handles.axesImage,'userdata'));
% tsns{1} = handles.animal.recording{rn}.contraFolderName;
% tsns{2} = handles.animal.recording{rn}.ipsiFolderName;
% tsns{3} = handles.animal.recording{rn}.binocularFolderName;
roiNumbers = get(handles.listbox_ROIs,'Value');
% 
% if length(roiNumber)>1
%     display('Select only one ROI');
%     return;
% end
tSeriesList = get(handles.listbox_allTSeries,'String');
tempString = tSeriesList;
returnValue = generalGUIForSelection(tempString,1);
if returnValue == 0
    return;
end
ROIs = handles.ROIs;
for ii = 1:length(roiNumbers)
    roiNumber = roiNumbers(ii);
    ROI = ROIs{roiNumber};
    instanceToBeCopied = getInstanceNumberFromROI(ROI,handles.animal.recording{rn}.TSeriesFolder{returnValue});
    if ~instanceToBeCopied
        continue;
    end
    instanceNumber = getInstanceNumberFromROI(ROI,handles.animal.recording{rn}.TSeriesFolder{tn});
    if instanceNumber == 0
        instanceNumber = length(ROI.instances) + 1;
    end
    ROI.instances{instanceNumber} = ROI.instances{instanceToBeCopied};
    ROI.instances{instanceNumber}.recordingNumber = rn;
    ROI.instances{instanceNumber}.folderName = handles.animal.recording{rn}.TSeriesFolder{tn}.name;
    ROI.instances{instanceNumber}.identificationFrame = [];
    ROIs{roiNumber} = ROI;
end
handles.ROIs = ROIs;
guidata(gcbo,handles);
saveROIs(handles);
value = get(handles.axesImage,'userdata');
frameNumber = get(handles.textCurrentFrame,'UserData');
loadFrame(handles,value,frameNumber);


% --- Executes on button press in pushbutton_deleteROI.
function pushbutton_deleteROI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_deleteROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~confirmToProceed
    return;
end
roiListNumbers = get(handles.listbox_ROIs,'Value');
ROIs = handles.ROIs;
ROIs(roiListNumbers) = [];
handles.ROIs = ROIs;
guidata(gcbo,handles);
saveROIs(handles);
% uipanel_imageDisplaySelection_SelectionChangeFcn(handles.uipanel_imageDisplaySelection, eventdata, handles);
tempString = cellstr(num2str((1:length(ROIs))'));
set(handles.listbox_ROIs,'Value',1);
if length(ROIs) == 0
    set(handles.listbox_ROIs,'String','');
else
    set(handles.listbox_ROIs,'String',tempString);
end


% --- Executes on selection change in listbox_allTSeries.
function listbox_allTSeries_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_allTSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_allTSeries contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_allTSeries
if get(handles.checkbox_plot,'value')
    pushbutton_plotResponses_Callback(handles.pushbutton_storeROIIdentities, [], handles);
    return;
end
value = get(hObject,'Value');
[rn tn] = findNumbers(handles,value);
mouseClickType = get(handles.figure1,'SelectionType');
if strcmp(mouseClickType,'open')
    set(handles.figure1,'SelectionType','normal');
    ROIs = handles.ROIs;
    ROIValues = [];
    for ii = 1:length(ROIs)
        ROI = ROIs{ii};
        tsn = handles.animal.recording{rn}.TSeriesFolder{tn}.name;
        if strcmp(ROI.instances{ROI.masterInstance}.folderName,tsn)
            ROIValues = [ROIValues ii];
            continue;
        end
    end
    if ~isempty(ROIValues)
        set(handles.listbox_ROIs,'Value',ROIValues);
        listbox_ROIs_Callback(handles.listbox_ROIs,[],handles);
    end
    return;
end

if length(value) > 1
    frameNumber = get(handles.textCurrentFrame,'userdata');
    loadFrame(handles,value(1),frameNumber);
    return;
end
frameNumber = get(handles.textCurrentFrame,'userdata');
loadFrame(handles,value,frameNumber);


% --- Executes during object creation, after setting all properties.
function listbox_allTSeries_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_allTSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_storeROIIdentities.
function pushbutton_storeROIIdentities_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_storeROIIdentities (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.listbox_allTSeries,'value');
[rn tn] = findNumbers(handles,value);
ROI_Numbers = get(handles.listbox_ROIs,'Value');
tempString = {'VIP' 'Gad2' 'PV' 'SOM' 'Astrocytes'};
returnValue = generalGUIForSelection(tempString,1);
if returnValue == 0
    return;
end
type = tempString{returnValue};
fileName = makeName('specialROIs.mat',handles.animal.recording{rn}.processedDataFolderPath);
save(fileName,'type','ROI_Numbers');

function [rn tsn] = findNumbers(handles,value)
[t r] = ind2sub([length(handles.animal.recording{handles.recordingNumbers(1)}.TSeriesFolder) length(handles.recordingNumbers)],value);
rn = [handles.recordingNumbers(r)]'; tsn = t';

% --- Executes on key press with focus on listbox_ROIs and none of its controls.
function listbox_ROIs_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to listbox_ROIs (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

if isempty(eventdata.Modifier)
    return;
end
if strcmp(eventdata.Modifier,'control') &&  strcmp(eventdata.Key,'delete')
    pushbutton_deleteROI_Callback(handles.pushbutton_deleteROI, 1, handles);
end


% --- Executes on button press in pushbutton_selectBloodVessels.
function pushbutton_selectBloodVessels_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_selectBloodVessels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.listbox_allTSeries,'Value');
[rn tn] = findNumbers(handles,value);
ROIs = handles.ROIs;
ROIValues = [];
for ii = 1:length(ROIs)
    ROI = ROIs{ii};
    tsn = handles.animal.recording{rn}.TSeriesFolder{tn}.name;
    if strcmp(lower(ROI.type),'blood vessel')
        ROIValues = [ROIValues ii];
        continue;
    end
end
if ~isempty(ROIValues)
    set(handles.listbox_ROIs,'Value',ROIValues);
    listbox_ROIs_Callback(handles.listbox_ROIs,[],handles);
end

% --- Executes on button press in checkbox_plot.
function checkbox_plot_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot


% --- Executes on button press in pushbutton_createNewROI.
function pushbutton_createNewROI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_createNewROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pushbutton_storePixels_Callback(handles.pushbutton_storePixels, 2, handles);


% --- Executes on button press in pushbutton_loadMaxMinusAverageImage.
function pushbutton_loadMaxMinusAverageImage_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadMaxMinusAverageImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,-2);


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmp(eventdata.Key,'rightarrow')
    pushbutton_nextTSeries_Callback(handles.pushbutton_nextTSeries, [], handles);
    return;
end
if strcmp(eventdata.Key,'leftarrow')
    pushbutton_previousTSeries_Callback(handles.pushbutton_previousTSeries, [], handles);
    return;
end
if strcmp(eventdata.Key,'control')
    if get(handles.checkbox_makeSubBox,'Value') == 1
        pushbutton_storePixels_Callback(handles.pushbutton_storePixels, 1, handles);
        handles = guidata(handles.figure1);
        value = get(handles.axesImage,'userdata');
        loadFrame(handles,value,0);
    else
        pushbutton_createNewSomaROI_Callback(hObject, eventdata, handles);
    end
    return;
end
if strcmp(eventdata.Key,'alt')
    pushbutton_storePixels_Callback(handles.pushbutton_storePixels, 4, handles);
    return;
end
if strcmp(eventdata.Key,'shift')
    pushbutton_selectROIToView_Callback(hObject, eventdata, handles);
end



% --- Executes on button press in pushbutton_createNewSomaROI.
function pushbutton_createNewSomaROI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_createNewSomaROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pushbutton_storePixels_Callback(handles.pushbutton_storePixels, 3, handles);


% --- Executes on button press in pushbutton_flashROIs.
function pushbutton_flashROIs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_flashROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'visible','off');
set(handles.pushbutton_stopFlash,'visible','on');
set(handles.checkbox_flashingROIs,'value',1);
while get(handles.checkbox_flashingROIs,'value')
    if get(handles.checkbox_displayROIs,'Value')
        set(handles.checkbox_displayROIs,'Value',0);
    else
        set(handles.checkbox_displayROIs,'Value',1);
    end
    checkbox_displayROIs_Callback(handles.checkbox_displayROIs, eventdata, handles);
    delayValue = get(handles.slider_playSpeedTSeries,'value');
    pause(delayValue);
end
set(hObject,'visible','on');
set(handles.pushbutton_stopFlash,'visible','off');

% --- Executes on button press in pushbutton_stopFlash.
function pushbutton_stopFlash_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_stopFlash (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.checkbox_flashingROIs,'value',0);
set(hObject,'visible','off');

% --- Executes on button press in checkbox_flashingROIs.
function checkbox_flashingROIs_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_flashingROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_flashingROIs


% --- Executes on button press in pushbutton_loadRefImage.
function pushbutton_loadRefImage_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadRefImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
frameNumber = -10000;
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,frameNumber);



% --- Executes on button press in pushbutton_processROI.
function pushbutton_processROI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_processROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'greenChannel');
    return;
end
value = get(handles.listbox_allTSeries,'value');
[RN TN] = findNumbers(handles,value);
if length(RN) > 1
    display('Choose one Recording');
    return;
end
% ROINumber = get(handles.listbox_ROIs,'value');
% if length(ROINumber) > 1
%     display('Choose one ROI');
%     return;
% end
recordingNumber = RN; recordingEyes = {'Contra' 'Ipsi' 'Binocular'};
hh = TN;
animal = handles.animal;
commandText = sprintf('processedDataFolder = animal.recording{recordingNumber}.processedDataFolderPath%s;',recordingEyes{hh});
eval(commandText);
commandText = sprintf('dataFolder = animal.recording{recordingNumber}.%sFolderPath;',lower(recordingEyes{hh}));
eval(commandText);
commandText = sprintf('tsn = animal.recording{recordingNumber}.%sFolderName;',lower(recordingEyes{hh}));
eval(commandText);
load(makeName('stimulusData.mat',processedDataFolder));
load(makeName('relativeTime.mat',processedDataFolder));
totalFrames = size(handles.greenChannel,3);
allGreenChannel = handles.greenChannel;
allStimulusChannel = handles.stimulusChannel;
ROIs = handles.ROIs;
numberOfStd = 3;%get(handles.slider_numberOfStd,'Value');
wh = waitbar(0,'Finding pixel values');
roiList = get(handles.listbox_ROIs,'Value');
for iii = 1:length(roiList)
    ii = roiList(iii);
    waitbar(iii/length(roiList),wh,sprintf('ROI Number %d/%d',iii,length(roiList)));
    ROI = ROIs{ii};
    instanceNumber = getInstanceNumberFromROI(ROI,tsn);
    if instanceNumber == 0
        continue;
    end
    box = ROI.instances{instanceNumber}.box;
    boxError = 0;
    try
        meanROI = mean(mean(allGreenChannel(box(2):box(4),box(1):box(3),:)));
        ROIFrames = allGreenChannel(box(2):box(4),box(1):box(3),:);
        sROIFrames = allStimulusChannel(box(2):box(4),box(1):box(3),:);
    catch
        boxError = 1;
        break;
    end
    if boxError
        continue;
    end
    meanROI = reshape(meanROI,1,numel(meanROI));
    globalMean = mean(reshape(ROIFrames,1,numel(ROIFrames)));
    globalStd = std(reshape(ROIFrames,1,numel(ROIFrames)));
    meanMeanROI = mean(meanROI);
    stdMeanROI = std(meanROI);
    meanROI = meanROI - meanMeanROI;
    index = find(meanROI == max(meanROI));
    frame = getFrame(animal,recordingNumber,tsn,index);
    pixelValues = frame(box(2):box(4),box(1):box(3));
    mPV = mean(mean(pixelValues));
    sPV = std(std(pixelValues));
%     poi = find(pixelValues > (mPV + numberOfStd * sPV));
%     poi = find(pixelValues > (meanMeanROI + numberOfStd * stdMeanROI));
    poi = find(pixelValues > (globalMean + numberOfStd * globalStd));
    [row col] = ind2sub(size(pixelValues),poi); 
    rowG = box(2) + row - 1; colG = box(1) + col - 1;
%     ROI.instances{instanceNumber}.pixels = sub2ind(size(frame),rowG,colG);
    ROIs{ii} = ROI;
end
while ~close(wh)
end
handles.ROIs = ROIs;
guidata(gcbo,handles);
uipanel_imageDisplaySelection_SelectionChangeFcn(handles.uipanel_imageDisplaySelection, eventdata, handles);
nothing = 0;
% 
% 
% 
% 
% ROI = handles.ROIs{ROINumber};
% instanceNumber = getInstanceNumberFromROI(ROI,tsn);
% if instanceNumber == 0
%     display('No ROI instance found');
%     return;
% end
% ROIpixels = ROI.instances{instanceNumber}.pixels;
% tic
% wh = waitbar(0,'Processing, please wait ...');
% for jj = 1:totalFrames
%     waitbar(jj/totalFrames,wh,sprintf('Frame Number %d/%d',jj,totalFrames));
%     greenChannel = handles.greenChannel(:,:,jj);
%     stimulusChannel = handles.stimulusChannel(:,:,jj);
%     averageIntensity(jj) = mean(mean(greenChannel(ROIpixels)));
%     averageIntensityStimulus(jj) = mean(mean(stimulusChannel(ROIpixels)));
% end
% while ~close(wh)
% end
% toc
% roi = processROI (averageIntensity,averageIntensityStimulus,stimulusData,relativeTime);
% roi.ROIpixels = ROIpixels;
% fileName = makeName('TestingForROIpixelsVariation\processedROI.mat',processedDataFolder);
% save(fileName,'roi');


% --- Executes on button press in checkbox_displayBoxes.
function checkbox_displayBoxes_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_displayBoxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_displayBoxes
uipanel_imageDisplaySelection_SelectionChangeFcn(handles.uipanel_imageDisplaySelection, eventdata, handles);


% --- Executes on button press in pushbutton_selectROIsInARegion.
function pushbutton_selectROIsInARegion_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_selectROIsInARegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'visible','off');
k = 0;
while k == 0
    k = waitforbuttonpress;
    point1 = get(gca,'CurrentPoint');    % button down detected
    finalRect = rbbox;                   % return figure units
    point2 = get(gca,'CurrentPoint');    % button up detected
    point1 = point1(1,1:2);point2 = point2(1,1:2);              % extract x and y
    x1 = min(floor(point1(1)),floor(point2(1)));x2 = max(floor(point1(1)),floor(point2(1)));
    y1 = min(floor(point1(2)),floor(point2(2)));y2 = max(floor(point1(2)),floor(point2(2)));
    if x1 < 1
        x1 = 1;
    end
    if y1 < 1
        y1 = 1;
    end
    if x2 > handles.frameSize(2)
        x2 = handles.frameSize(2);
    end
    if y2 > handles.frameSize(1)
        y2 = handles.frameSize(1);
    end
    numberOfPointsX = x2-x1+1;
    numberOfPointsY = y2-y1+1;
    totalPoints = numberOfPointsX * numberOfPointsY;
    k = 1;
end
bigBox = [x1 y1 x2 y2];
ROIs = handles.ROIs;
ROIValues = [];
for ii = 1:length(ROIs)
    ROI = ROIs{ii};
    for jj = 1:length(ROI.instances)
        box = ROI.instances{jj}.box;
        if box(1) > bigBox(1) && box(1) < bigBox(3) && box(2) > bigBox(2) && box(2) < bigBox(4)
            ROIValues = [ROIValues ii];
            break;
        end
    end
end
if ~isempty(ROIValues)
    currentROIValues = get(handles.listbox_ROIs,'Value');
    ROIValues = intersect(ROIValues,currentROIValues);
    set(handles.listbox_ROIs,'Value',ROIValues);
    listbox_ROIs_Callback(handles.listbox_ROIs,[],handles);
end
set(hObject,'visible','on');

% --- Executes on button press in pushbutton_copyROIs.
function pushbutton_copyROIs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_copyROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[rn tsnv] = findNumbers(handles,get(handles.listbox_allTSeries,'Value'));
roiListNumbers = get(handles.listbox_ROIs,'Value');
if length(rn) ~= 2
    display('Select two recordings');
    return;
end
allTSeriesValue = get(handles.listbox_allTSeries,'Value');
allTSeriesString = get(handles.listbox_allTSeries,'String');
tempString{1} = allTSeriesString{allTSeriesValue(1)};
tempString{2} = allTSeriesString{allTSeriesValue(2)};
returnValue = generalGUIForSelection(tempString,1,'Select Source TSeries');
if returnValue == 0
    return;
end
tsns = handles.animal.recording{rn(returnValue)}.TSeriesFolderName;
tsn = tsns(tsnv(returnValue));

if returnValue == 1
    destinationRecordingNumber = 2;
end
if returnValue == 2
    destinationRecordingNumber = 1;
end
ROIs = handles.ROIs;
for ii = 1:length(roiListNumbers)
    roiNumber = roiListNumbers(ii);
    ROI = ROIs{roiNumber};
    sourceInstanceNumber = getInstanceNumberFromROI(ROI,tsn);
    if sourceInstanceNumber == 0
        continue;
    end
    instanceNumber = length(ROI.instances)+1;
    ROI.instances{instanceNumber} = ROI.instances{sourceInstanceNumber};
    ROI.instances{instanceNumber}.recordingNumber = rn(destinationRecordingNumber);
    dtsns = handles.animal.recording{rn(destinationRecordingNumber)}.TSeriesFolderName;
    dtsn = dtsns{tsnv(destinationRecordingNumber)};
    ROI.instances{instanceNumber}.folderName = dtsn;
    ROI.instances{instanceNumber}.identificationFrame = [];
    ROIs{roiNumber} = ROI;
end
handles.ROIs = ROIs;
guidata(gcbo,handles);
saveROIs(handles);
frameNumber = get(handles.textCurrentFrame,'UserData');
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,frameNumber);


% --- Executes on button press in pushbutton_moveBox.
function pushbutton_moveBox_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_moveBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roiNumber = get(handles.listbox_ROIs,'Value');
if length(roiNumber) > 1
    return;
end
ROIs = handles.ROIs;
ROI = ROIs{roiNumber};
[rn tsnv] = findNumbers(handles,get(handles.listbox_allTSeries,'Value'));
if length(rn) > 1
    return;
end
tsns = handles.animal.recording{rn}.TSeriesFolderName;
tsn = tsns(tsnv);
instanceNumber = getInstanceNumberFromROI(ROI,tsn);
if instanceNumber == 0
    return;
end
roiToChange.ROI = ROI;
roiToChange.roiNumber = roiNumber;
roiToChange.instanceNumber = instanceNumber;
set(handles.pushbutton_moveBox,'userdata',roiToChange);
set(handles.figure1,'WindowButtonDownFcn',@wbdcb);


function wbdcb(src,evnt) % left button down function
handles = guidata(src);
if strcmp(get(src,'SelectionType'),'normal')   
    set(src,'WindowButtonDownFcn','');
    cp = get(handles.axesImage,'CurrentPoint');
    x1 = ceil(cp(1,1)); y1 = ceil(cp(1,2));
    if x1 < 1 | x1 > handles.frameSize(2)
        return;
    end
    if y1 < 1 | y1 > handles.frameSize(1)
        return;
    end
    clickedPoint = [x1 y1];
    roiToChange = get(handles.pushbutton_moveBox,'userdata');
    roiToChange.clickedPoint = clickedPoint;
    ROI = roiToChange.ROI;
    box = ROI.instances{roiToChange.instanceNumber}.box;
    roiToChange.rh = displayBox(handles,box,0);
    set(handles.pushbutton_moveBox,'userdata',roiToChange);
    set(src,'WindowButtonMotionFcn',@wbmcb);
    set(src,'WindowButtonUpFcn',@wbucb);
end

function wbmcb(src,evnt) % button motion function ... executes when the mouse is in motion
handles = guidata(src);
cp = get(handles.axesImage,'CurrentPoint');
x1 = ceil(cp(1,1)); y1 = ceil(cp(1,2));
if x1 < 1 | x1 > handles.frameSize(2)
    return;
end
if y1 < 1 | y1 > handles.frameSize(1)
    return;
end
currentPoint = [x1 y1];
roiToChange = get(handles.pushbutton_moveBox,'userdata');
oldrh = roiToChange.rh;
delete(oldrh);
ROI = roiToChange.ROI;
box = ROI.instances{roiToChange.instanceNumber}.box;
clickedPoint = roiToChange.clickedPoint;
displacement = currentPoint - clickedPoint;
box(1) = box(1) + displacement(1);
box(3) = box(3) + displacement(1);
box(2) = box(2) + displacement(2);
box(4) = box(4) + displacement(2);
roiToChange.rh = displayBox(handles,box,0);
roiToChange.box = box;
set(handles.pushbutton_moveBox,'userdata',roiToChange);

function wbucb(src,evnt) % button up function
handles = guidata(src);
set(src,'WindowButtonMotionFcn','');
set(src,'WindowButtonUpFcn','');
roiToChange = get(handles.pushbutton_moveBox,'userdata');
oldrh = roiToChange.rh;
delete(oldrh);
ROI = roiToChange.ROI;
ROI.instances{roiToChange.instanceNumber}.box = roiToChange.box;
ROIs = handles.ROIs;
ROIs{roiToChange.roiNumber} = ROI;
handles.ROIs = ROIs;
guidata(src,handles);
frameNumber = get(handles.textCurrentFrame,'UserData');
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,frameNumber);


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user 


function anotherwbdcb(src,evnt) % left button down function
handles = guidata(src);
if strcmp(get(src,'SelectionType'),'open')
    set(src,'SelectionType','normal');
    pushbutton_moveBox_Callback(handles.pushbutton_moveBox, [], handles);
end


% --- Executes on button press in pushbutton_findPixels.
function pushbutton_findPixels_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_findPixels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[rn tsnv] = findNumbers(handles,get(handles.listbox_allTSeries,'Value'));


% --- Executes on button press in pushbutton_frameShiftLeft.
function pushbutton_frameShiftLeft_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_frameShiftLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listingNumber = get(handles.axesImage,'userdata');
[rn tsnv] = findNumbers(handles,listingNumber);
fileName = handles.animal.recording{rn}.TSeriesFrameShiftFile{tsnv};
if nameExists(fileName)
    load(fileName);
else
    shiftValues = createNewFrameShiftFile(fileName);
end
shiftValues.xShift = shiftValues.xShift - 1;
save(fileName,'shiftValues');
frameNumber = get(handles.textCurrentFrame,'UserData');
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,frameNumber);


% --- Executes on button press in pushbutton_frameShiftRight.
function pushbutton_frameShiftRight_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_frameShiftRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listingNumber = get(handles.axesImage,'userdata');
[rn tsnv] = findNumbers(handles,listingNumber);
fileName = handles.animal.recording{rn}.TSeriesFrameShiftFile{tsnv};
if nameExists(fileName)
    load(fileName);
else
    shiftValues = createNewFrameShiftFile(fileName);
end
shiftValues.xShift = shiftValues.xShift + 1;
save(fileName,'shiftValues');
frameNumber = get(handles.textCurrentFrame,'UserData');
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,frameNumber);

% --- Executes on button press in pushbutton_frameShiftUp.
function pushbutton_frameShiftUp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_frameShiftUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listingNumber = get(handles.axesImage,'userdata');
[rn tsnv] = findNumbers(handles,listingNumber);
fileName = handles.animal.recording{rn}.TSeriesFrameShiftFile{tsnv};
if nameExists(fileName)
    load(fileName);
else
    shiftValues = createNewFrameShiftFile(fileName);
end
shiftValues.yShift = shiftValues.yShift - 1;
save(fileName,'shiftValues');
frameNumber = get(handles.textCurrentFrame,'UserData');
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,frameNumber);

% --- Executes on button press in pushbutton_frameShiftDown.
function pushbutton_frameShiftDown_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_frameShiftDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listingNumber = get(handles.axesImage,'userdata');
[rn tsnv] = findNumbers(handles,listingNumber);
fileName = handles.animal.recording{rn}.TSeriesFrameShiftFile{tsnv};
if nameExists(fileName)
    load(fileName);
else
    shiftValues = createNewFrameShiftFile(fileName);
end
shiftValues.yShift = shiftValues.yShift + 1;
save(fileName,'shiftValues');
frameNumber = get(handles.textCurrentFrame,'UserData');
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,frameNumber);

% --- Executes on button press in pushbutton_frameRotateCW.
function pushbutton_frameRotateCW_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_frameRotateCW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listingNumber = get(handles.axesImage,'userdata');
[rn tsnv] = findNumbers(handles,listingNumber);
fileName = handles.animal.recording{rn}.TSeriesFrameShiftFile{tsnv};
if nameExists(fileName)
    load(fileName);
else
    shiftValues = createNewFrameShiftFile(fileName);
end
shiftValues.rotation = shiftValues.rotation - 0.5;
save(fileName,'shiftValues');
frameNumber = get(handles.textCurrentFrame,'UserData');
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,frameNumber);

% --- Executes on button press in pushbutton_frameRotateCCW.
function pushbutton_frameRotateCCW_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_frameRotateCCW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listingNumber = get(handles.axesImage,'userdata');
[rn tsnv] = findNumbers(handles,listingNumber);
fileName = handles.animal.recording{rn}.TSeriesFrameShiftFile{tsnv};
if nameExists(fileName)
    load(fileName);
else
    shiftValues = createNewFrameShiftFile(fileName);
end
shiftValues.rotation = shiftValues.rotation + 0.5;
save(fileName,'shiftValues');
frameNumber = get(handles.textCurrentFrame,'UserData');
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,frameNumber);


function shiftValues = createNewFrameShiftFile(fileName)
shiftValues.xShift = 0;
shiftValues.yShift = 0;
shiftValues.rotation = 0;
save(fileName,'shiftValues');


% --- Executes on button press in checkbox_displayNeuropil.
function checkbox_displayNeuropil_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_displayNeuropil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_displayNeuropil


% --- Executes on button press in checkbox_displayTimeLine.
function checkbox_displayTimeLine_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_displayTimeLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_displayTimeLine


% --- Executes on button press in pushbutton_redChannelImage.
function pushbutton_redChannelImage_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_redChannelImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
frameNumber = get(handles.textCurrentFrame,'UserData');
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,-3);


% --- Executes on button press in pushbutton_loadCorrectedAverageImage.
function pushbutton_loadCorrectedAverageImage_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadCorrectedAverageImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,0.1);


% --- Executes on button press in pushbutton_loadStimulusAverageImage.
function pushbutton_loadStimulusAverageImage_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadStimulusAverageImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,0.2);

% --- Executes on button press in pushbutton_loadBaseLineAverageImage.
function pushbutton_loadBaseLineAverageImage_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadBaseLineAverageImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,0.3);

% --- Executes on button press in pushbutton_loadStimulusMinusBaseLineAverageImage.
function pushbutton_loadStimulusMinusBaseLineAverageImage_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadStimulusMinusBaseLineAverageImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,0.4);


% --- Executes on slider movement.
function slider_minLUT_Callback(hObject, eventdata, handles)
% hObject    handle to slider_minLUT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_minLUT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_minLUT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton_createBloodVesselROI.
function pushbutton_createBloodVesselROI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_createBloodVesselROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pushbutton_storePixels_Callback(handles.pushbutton_storePixels, 4, handles);



function edit_channelNumber_Callback(hObject, eventdata, handles)
% hObject    handle to edit_channelNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_channelNumber as text
%        str2double(get(hObject,'String')) returns contents of edit_channelNumber as a double


% --- Executes during object creation, after setting all properties.
function edit_channelNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_channelNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_makeSubBox.
function checkbox_makeSubBox_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_makeSubBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_makeSubBox


% --- Executes on button press in pushbutton_fuseImages.
function pushbutton_fuseImages_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_fuseImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
channelNumber = str2num(get(handles.edit_channelNumber,'String'));
value = get(handles.axesImage,'userdata');
[rn tn] = findNumbers(handles,value);
tsn = handles.animal.recording{rn}.TSeriesFolder{tn}.name;
tsnn = handles.animal.recording{rn}.TSeriesFolder{tn}.eye;

rFrameImage = getFrame(handles.animal,rn,tsn,-10000,channelNumber);
rFrameImageRGB = repmat(rFrameImage/max(rFrameImage(:)), [1 1 3]);
rFrameImageRGB(:,:,3) = zeros(size(rFrameImageRGB(:,:,1)));
rFrameImageRGB(:,:,2) = zeros(size(rFrameImageRGB(:,:,1)));
aFrameImage = getFrame(handles.animal,rn,tsn,-10002,channelNumber);
aFrameImageRGB = repmat(aFrameImage/max(aFrameImage(:)), [1 1 3]);
aFrameImageRGB(:,:,1) = zeros(size(aFrameImageRGB(:,:,1)));
aFrameImageRGB(:,:,3) = zeros(size(aFrameImageRGB(:,:,1)));

fFrameImage(:,:,1) = rFrameImageRGB(:,:,1);
fFrameImage(:,:,2) = aFrameImageRGB(:,:,2);
fFrameImage(:,:,3) = zeros(size(rFrameImageRGB(:,:,1)));

% fFrameImage = double(imfuse(aFrameImage,rFrameImage));
% fFrameImage = fFrameImage - min(min(min(fFrameImage)));
% fFrameImage = fFrameImage/max(max(max(fFrameImage)));

fileName = makeName('fusedImage.mat',handles.animal.recording{rn}.processedDataFolderPath);
save(fileName,'fFrameImage');
value = get(handles.axesImage,'userdata');
loadFrame(handles,value,-10001);


% --- Executes on button press in pushbutton_openExplorer.
function pushbutton_openExplorer_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_openExplorer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[rn tn] = findNumbers(handles,1);
commandText = sprintf('explorer %s',handles.animal.recording{rn}.folderPath);
dos(commandText);


% --- Executes on button press in checkbox_enhanceContrast.
function checkbox_enhanceContrast_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_enhanceContrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_enhanceContrast
value = get(handles.axesImage,'userdata');
frameNumber = get(handles.textCurrentFrame,'UserData');
loadFrame(handles,value,frameNumber);

% --- Executes on button press in pushbutton_markElectrodeROI.
function pushbutton_markElectrodeROI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_markElectrodeROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roiNumber = get(handles.listbox_ROIs,'value');
fileName = makeName('electrodeROI.mat',handles.animal.recording{handles.recordingNumbers(1)}.processedDataFolderPath);
save(fileName,'roiNumber');


% --- Executes on button press in pushbutton_modifyROI.
function pushbutton_modifyROI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_modifyROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.axesImage,'userdata');
frameNumber = get(handles.textCurrentFrame,'UserData');
loadFrame(handles,value,frameNumber,-1);
pushbutton_selectROIToView_Callback(handles.pushbutton_selectROIToView, 1, handles);


% --- Executes when selected object is changed in uipanel10.
function uipanel10_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel10 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.axesImage,'userdata');
frameNumber = get(handles.textCurrentFrame,'UserData');
loadFrame(handles,value,frameNumber);


% --- Executes on button press in pushbutton_OpenProcFolder.
function pushbutton_OpenProcFolder_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OpenProcFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[rn tn] = findNumbers(handles,1);
commandText = sprintf('explorer %s',handles.animal.recording{rn}.processedDataFolderPath);
dos(commandText);


% --- Executes on button press in pushbutton8_ROIShrink.
function pushbutton8_ROIShrink_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8_ROIShrink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
moveROIs(handles,'shrink');

% --- Executes on button press in pushbutton_ROIExpand.
function pushbutton_ROIExpand_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ROIExpand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
moveROIs(handles,'expand');


% --- Executes on button press in pushbutton_selectROIs.
function pushbutton_selectROIs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_selectROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.listbox_allTSeries,'value');
[rn tn] = findNumbers(handles,value);
% tempString = {'VIP' 'Gad2' 'PV' 'SOM' 'Astrocytes'};
% returnValue = generalGUIForSelection(tempString,1);
% if returnValue == 0
%     return;
% end
% type = tempString{returnValue};
fileName = makeName('specialROIs.mat',handles.animal.recording{rn}.processedDataFolderPath);
load(fileName);
set(handles.listbox_ROIs,'Value',ROI_Numbers);
value = get(handles.axesImage,'userdata');
frameNumber = get(handles.textCurrentFrame,'UserData');
loadFrame(handles,value,frameNumber);
