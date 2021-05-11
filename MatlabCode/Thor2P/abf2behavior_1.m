function o = abf2behavior_1(ei,saveDataFolder,varargin)
%
p = inputParser;
default_overwrite_behavior = 0;
default_downsample = 500;
addRequired(p,'ei',@isstruct);
addRequired(p,'saveDataFolder',@ischar);
addOptional(p,'overwrite_behavior',default_overwrite_behavior,@isnumeric);
addOptional(p,'downsample',default_downsample,@isnumeric);
parse(p,ei,saveDataFolder,varargin{:});
overwrite = p.Results.overwrite_behavior;
downSampleFactor = p.Results.downsample;

fileName = makeName(sprintf('behavior%d.mat',ei.db.selectedPlane),saveDataFolder);
if exist(fileName,'file') && overwrite == 0
    o = load(fileName);
else
db = ei.db;

filename = ei.abf_file;
[d,si] = abf2load(filename);

if size(d,2) ~= length(db.channel)
    display('Check your meta data file ... number of channels not consistent with number of channel names provided');
    error;
end
o.channelsChanged = 0;
i_channels = identify_abf_channels(d,si);
if ~isequal(db.channel,i_channels)
    db.channel = i_channels;
    disp('Channel designation in xml file is not correct, please check ... using identified channels');
    o.channelsChanged = 1;
end

temp = ~cellfun(@isempty,(strfind(lower(db.channel),'lfp')));
[val,idx] = find(temp);
if ~isempty(val)
    lfp = d(:,idx);
    thisfn = makeName('lfp.mat',saveDataFolder);
    save(thisfn,'lfp','-v7.3');
    idxs = 1:size(d,2);
    idxs(idx) = [];
    d = d(:,idxs);
%     d = downsample(d,5);
end

for ii = 1:size(d,2)
    channel = db.channel{ii};
    cmdText = sprintf('o.%s_f = find_falling_edge(d(:,ii),-0.5,2);',channel);
    eval(cmdText);
    cmdText = sprintf('o.%s_r = find_rising_edge(d(:,ii),0.5,2);',channel);
    eval(cmdText);
    if strcmp(channel,'air_puff')
        stim = d(:,ii)';
        stim = stim - min(stim);
        o.air_puff_raw = stim/max(stim);
    end
    if strcmp(channel,'photo_sensor')
        stim = d(:,ii)';
        stim = stim - min(stim);
        o.photo_sensor_raw = stim/max(stim);
    end
    if strcmp(channel,'stim')
        stim = d(:,ii)';
        stim = stim - min(stim);
        o.stim_raw = stim/max(stim);
    end
end
o.number_of_samples = size(d(:,1));
o.si = si;

%%%%%%%
% new code for finding dist January 18, 2018
% based on both encoder signals to see which direction the movement was
ind = find(strcmp(db.channel,'ch_a'));
cha = d(:,ind);

ind = find(strcmp(db.channel,'ch_b'));
chb = d(:,ind);
o.encoderCount = processEncodeSignals(cha,chb);
%%%%%%% 
data_rate = 1/(si*1e-6);
if ei.zFastEnable
    frames_diff_threshold = floor(data_rate/((ei.frameRate)*(ei.zSteps+1))+25);
else
    frames_diff_threshold = floor(data_rate/((ei.frameRate))+25);
end
% disp(sprintf('Using frame diff threshold %d',frames_diff_threshold));
dFf = diff(o.frames_f);
temp = find(dFf > frames_diff_threshold);
o.frames_f(temp+1) = [];

dFr = diff(o.frames_r);
temp = find(dFr > frames_diff_threshold);
o.frames_r(temp+1) = [];
if isfield(o,'air_puff_f')
    [o.air_puff_f,apitor] = getRidOfCloseRepetitions(o.air_puff_f,5000);
    if ~isempty(apitor)
        o.air_puff_r(apitor) = [];
    end
    [o.air_puff_r,apitor] = getRidOfCloseRepetitions(o.air_puff_r,5000);
    if ~isempty(apitor)
        o.air_puff_f(apitor) = [];
    end
end
if ei.zFastEnable
    try
        ii = ei.db.selectedPlane;
        nFrames = ei.totalFrames;
        framesToSelect = ii:(ei.zSteps+1):(nFrames*(ei.zSteps+1));
        o.frames_f = o.frames_f(framesToSelect);
        o.frames_r = o.frames_r(framesToSelect);
    catch
        if length(o.frames_f) < ei.streaming_frames
            diffLengths = ei.streaming_frames - length(o.frames_f);
            remain = mod(diffLengths,ei.zSteps+1);
            ii = remain + ei.db.selectedPlane;
            nFrames = length(o.frames_f);
            framesToSelect = ii:(ei.zSteps+1):nFrames;
            o.frames_f = o.frames_f(framesToSelect);
            o.frames_r = o.frames_r(framesToSelect);
        else
            error;
        end
    end
end

if isfield(o,'air_puff_f')
    if o.air_puff_f(1) < o.air_puff_r(1)
        o.air_puff_f(1) = [];
    end

    if length(o.air_puff_f) ~= length(o.air_puff_r)
%         error('air_puff_f is not equal to air_puff_r');
    end
end
if abs(length(o.frames_f) - length(o.frames_r)) > 2
    error('frames_f is not equal to frames_r');
end

ts = (0:(size(d(:,1))-1)) * si * 1e-6;
hf = figure(77777);clf;
frames = zeros(size(ts));
frames(o.frames_f) = 0.5;
subplot 211;
plot(ts,frames,'linewidth',0.25);hold on;
if isfield(o,'air_puff_f')
    plot(ts,o.air_puff_raw,'linewidth',4);
    for ii = 1:length(o.air_puff_f)
        text(ts(o.air_puff_f(ii))+5,0.75,num2str(ii));
    end
end
subplot 212;
frames = zeros(size(ts));
frames(o.ch_a_r) = 0.5;
plot(ts,frames,'linewidth',0.25);hold on;
if isfield(o,'air_puff_f')
    plot(ts,o.air_puff_raw,'linewidth',4);
end

set(gcf,'Position',get(0,'ScreenSize'));

% prompt = {'Enter trials to use'};
% dlg_title = 'Input';
% num_lines = 1;
% defaultans = {''};
% answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
% if strfind(answer{1},':')
%     cmdTxt = sprintf('o.trials = [%s];',answer{1});
%     eval(cmdTxt);
% else
%     o.trials = str2num(answer{1});
% end
if isfield(o,'air_puff_f')
o.trials = 1:length(o.air_puff_r);
end

try
    save(fileName,'-struct','o','-v7.3');
catch
    mkdir(saveDataFolder);
    try
        save(fileName,'-struct','o','-v7.3');
    catch
        disp('unknown error');
        lasterror;
    end
end

close(hf);
end

function [edge,temp1] = getRidOfCloseRepetitions(edge,minimum_diff)
temp = find(diff(edge)<=minimum_diff);
edge(temp+1) = [];
temp1 = temp + 1;

function edge = find_rising_edge(signal,threshold,minimum_diff)
edge = find(diff(signal) >= threshold);
temp = find(diff(edge)<=minimum_diff);
edge(temp+1) = [];

function edge = find_falling_edge(signal,threshold,minimum_diff)
edge = find(diff(signal) <= threshold);
temp = find(diff(edge)<=minimum_diff);
edge(temp+1) = [];

function dist = processEncodeSignals(cha,chb)
n = 0;
chat = cha > 2.5;
chbt = chb > 2.5;
encoderCount = 0;
valP = [chat(1) chbt(1)];
dist(1) = encoderCount;
for ii = 2:length(cha)
    valC = [chat(ii) chbt(ii)];
    if valC == valP
        valP = valC;
        dist(ii) = encoderCount;
        continue;
    end
    if isequal(valP,[0 0]) && isequal(valC,[1 0])
        encoderCount = encoderCount + 1;
        valP = valC;
        dist(ii) = encoderCount;
        continue;
    end
    if isequal(valP,[0 0]) && isequal(valC,[0 1])
        encoderCount = encoderCount - 1;
        valP = valC;
        dist(ii) = encoderCount;
        continue;
    end
    if isequal(valP,[1 0]) && isequal(valC,[1 1])
        encoderCount = encoderCount + 1;
        valP = valC;
        dist(ii) = encoderCount;
        continue;
    end
    if isequal(valP,[1 0]) && isequal(valC,[0 0])
        encoderCount = encoderCount - 1;
        valP = valC;
        dist(ii) = encoderCount;
        continue;
    end
    if isequal(valP,[0 1]) && isequal(valC,[0 0])
        encoderCount = encoderCount + 1;
        valP = valC;
        dist(ii) = encoderCount;
        continue;
    end
    if isequal(valP,[0 1]) && isequal(valC,[1 1])
        encoderCount = encoderCount - 1;
        valP = valC;
        dist(ii) = encoderCount;
        continue;
    end
    if isequal(valP,[1 1]) && isequal(valC,[0 1])
        encoderCount = encoderCount + 1;
        valP = valC;
        dist(ii) = encoderCount;
        continue;
    end
    if isequal(valP,[1 1]) && isequal(valC,[1 0])
        encoderCount = encoderCount - 1;
        valP = valC;
        dist(ii) = encoderCount;
        continue;
    end
end

