function o = process_abf_data(d,si,db)

if size(d,2) ~= length(db.channel)
    display('Check your meta data file ... number of channels not consistent with number of channel names provided');
    error;
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
[o.air_puff_f,apitor] = getRidOfCloseRepetitions(o.air_puff_f,5000);
if ~isempty(apitor)
    o.air_puff_r(apitor) = [];
end
[o.air_puff_r,apitor] = getRidOfCloseRepetitions(o.air_puff_r,5000);
if ~isempty(apitor)
    o.air_puff_f(apitor) = [];
end

if o.air_puff_f(1) < o.air_puff_r(1)
    o.air_puff_f(1) = [];
end

if length(o.air_puff_f) ~= length(o.air_puff_r)
    n = 0;
end

ts = (0:(size(d(:,1))-1)) * si * 1e-6;
o.trials = 1:length(o.air_puff_r);


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

