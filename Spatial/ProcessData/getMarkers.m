function [markersOn,markersOff] = getMarkers(ei,context,trials,markerType,pMarkers)

markersOn = [];
markersOff = [];

if strcmp(lower(markerType),'air')
%     if strcmp(lower(context),'air - brake') %size(trials,1) > 1 | 
%         trials = trials(1,:);
%         onsets = ei.b.air_puff_r(trials);
%         offsets = ei.b.air_puff_f(trials);
%         toffsets = offsets + round(1e6 * 1/ei.b.si);
%         tonsets = onsets - round(1e6 * 1/ei.b.si);
%         markersOn = tonsets;
%         markersOff = toffsets;
%     else
        onsets = ei.b.air_puff_r(trials);
        offsets = ei.b.air_puff_f(trials);
        markersOn = onsets;
        markersOff = offsets;
%     end
end


if strcmp(lower(markerType),'airr')
%     if strcmp(lower(context),'air - brake') %size(trials,1) > 1 | 
%         trials = trials(1,:);
%         onsets = ei.b.air_puff_r(trials);
%         offsets = ei.b.air_puff_f(trials);
%         toffsets = offsets + round(1e6 * 1/ei.b.si);
%         tonsets = onsets - round(1e6 * 1/ei.b.si);
%         markersOn = tonsets;
%         markersOff = toffsets;
%     else
        onsets = ei.b.air_puff_r(trials);
        offsets = ei.b.air_puff_f(trials);
        markersOn = onsets;
        markersOff = offsets;
%     end
end

if strcmp(lower(markerType),'air44')
    onsets = ei.b.air_puff_r(trials);
    offsets = ei.b.air_puff_f(trials);
    toffsets = offsets + round(1e6 * 4/ei.b.si);
    tonsets = onsets - round(1e6 * 4/ei.b.si);
    markersOn = tonsets;
    markersOff = toffsets;
end

if strcmp(lower(markerType),'air33')
    onsets = ei.b.air_puff_r(trials);
    offsets = ei.b.air_puff_f(trials);
    toffsets = offsets + round(1e6 * 3/ei.b.si);
    tonsets = onsets - round(1e6 * 3/ei.b.si);
    markersOn = tonsets;
    markersOff = toffsets;
end

if strcmp(lower(markerType),'air77')
    onsets = ei.b.air_puff_r(trials);
    offsets = ei.b.air_puff_f(trials);
    toffsets = offsets + round(1e6 * 7/ei.b.si);
    tonsets = onsets - round(1e6 * 7/ei.b.si);
    markersOn = tonsets;
    markersOff = toffsets;
end

if strcmp(lower(markerType),'air55')
    onsets = ei.b.air_puff_r(trials);
    offsets = ei.b.air_puff_f(trials);
    toffsets = offsets + round(1e6 * 5/ei.b.si);
    tonsets = onsets - round(1e6 * 5/ei.b.si);
    markersOn = tonsets;
    markersOff = toffsets;
end

if strcmp(lower(markerType),'light')
%     if size(trials,1) > 1
%         trials = trials(2,:);
%     end
%     if strcmp(lower(context),'light - brake')
        onsets = ei.b.stim_r(trials);
        offsets = ei.b.stim_f(trials);
        toffsets = offsets;% + round(1e6 * 1/ei.b.si);
        tonsets = onsets;% - round(1e6 * 1/ei.b.si);
        markersOn = tonsets;
        markersOff = toffsets;
%     end
end

if strcmp(lower(markerType),'light11')
%     if size(trials,1) > 1
%         trials = trials(2,:);
%     end
%     if strcmp(lower(context),'light - brake')
        onsets = ei.b.stim_r(trials);
        offsets = ei.b.stim_f(trials);
        toffsets = offsets + round(1e6 * 1/ei.b.si);
        tonsets = onsets - round(1e6 * 1/ei.b.si);
        markersOn = tonsets;
        markersOff = toffsets;
%     end
end

if strcmp(lower(markerType),'light22')
%     if size(trials,1) > 1
%         trials = trials(2,:);
%     end
%     if strcmp(lower(context),'light - brake')
        onsets = ei.b.stim_r(trials);
        offsets = ei.b.stim_f(trials);
        toffsets = offsets + round(1e6 * 2/ei.b.si);
        tonsets = onsets - round(1e6 * 2/ei.b.si);
        markersOn = tonsets;
        markersOff = toffsets;
%     end
end

if strcmp(lower(markerType),'light0p30p3')
%     if size(trials,1) > 1
%         trials = trials(2,:);
%     end
%     if strcmp(lower(context),'light - brake')
        onsets = ei.b.stim_r(trials);
        offsets = ei.b.stim_f(trials);
        toffsets = offsets + round(1e6 * 0.3/ei.b.si);
        tonsets = onsets - round(1e6 * 0.3/ei.b.si);
        markersOn = tonsets;
        markersOff = toffsets;
%     end
end

if strcmp(lower(markerType),'tone')
%     if size(trials,1) > 1
%         trials = trials(2,:);
%     end
%     if strcmp(lower(context),'tone - brake')
        onsets = ei.b.stim_r(trials);
        offsets = ei.b.stim_f(trials);
        toffsets = offsets;% + round(1e6 * 1/ei.b.si);
        tonsets = onsets;% - round(1e6 * 1/ei.b.si);
        markersOn = tonsets;
        markersOff = toffsets;
%     end
end

if strcmp(lower(markerType),'tone11')
%     if size(trials,1) > 1
%         trials = trials(2,:);
%     end
%     if strcmp(lower(context),'tone - brake')
        onsets = ei.b.stim_r(trials);
        offsets = ei.b.stim_f(trials);
        toffsets = offsets + round(1e6 * 1/ei.b.si);
        tonsets = onsets - round(1e6 * 1/ei.b.si);
        markersOn = tonsets;
        markersOff = toffsets;
%     end
end

if strcmp(lower(markerType),'tone22')
%     if size(trials,1) > 1
%         trials = trials(2,:);
%     end
%     if strcmp(lower(context),'tone - brake')
        onsets = ei.b.stim_r(trials);
        offsets = ei.b.stim_f(trials);
        toffsets = offsets + round(1e6 * 2/ei.b.si);
        tonsets = onsets - round(1e6 * 2/ei.b.si);
        markersOn = tonsets;
        markersOff = toffsets;
%     end
end

if strcmp(lower(markerType),'tone0p30p3')
%     if size(trials,1) > 1
%         trials = trials(2,:);
%     end
%     if strcmp(lower(context),'tone - brake')
        onsets = ei.b.stim_r(trials);
        offsets = ei.b.stim_f(trials);
        toffsets = offsets + round(1e6 * 0.3/ei.b.si);
        tonsets = onsets - round(1e6 * 0.3/ei.b.si);
        markersOn = tonsets;
        markersOff = toffsets;
%     end
end

if strcmp(lower(markerType),'airi')
    onsets = ei.b.air_puff_r(trials);
    offsets = ei.b.air_puff_f(trials);
    tonsets = offsets;%(1:(end-1));% + round(1e6 * 3/ei.b.si);
    toffsets = onsets(2:end);% - round(1e6 * 0.5/ei.b.si);
    dur = mean(ei.b.ts(toffsets(1:5))-ei.b.ts(tonsets(1:5)));
    toffsets(length(toffsets)+1) = tonsets(end) + round(1e6 * dur/ei.b.si);
    if length(onsets) == 30 % to handle combined conditions 3 4 5
        toffsets(10) = tonsets(10) + round(1e6 * dur/ei.b.si);
        toffsets(20) = tonsets(20) + round(1e6 * dur/ei.b.si);
    end
    markersOn = tonsets;
    markersOff = toffsets;
end

if strcmp(lower(markerType),'airir')
    onsets = ei.b.air_puff_r(trials);
    offsets = ei.b.air_puff_f(trials);
    tonsets = offsets;%(1:(end-1));% + round(1e6 * 3/ei.b.si);
    toffsets = onsets(2:end);% - round(1e6 * 0.5/ei.b.si);
    dur = mean(ei.b.ts(toffsets(1:5))-ei.b.ts(tonsets(1:5)));
    toffsets(length(toffsets)+1) = tonsets(end) + round(1e6 * dur/ei.b.si);
    markersOn = tonsets;
    markersOff = toffsets;
end

if strcmp(lower(markerType),'belt')
    if isfield(ei.b,'air_puff_r')
        onsets = ei.b.air_puff_r(trials);
        offsets = ei.b.air_puff_f(trials);
    else
        offsets = ei.b.photo_sensor_f(end) + round(1e6 * 0.3/ei.b.si);
        onsets = ei.b.photo_sensor_f(1) - round(1e6 * 0.3/ei.b.si);
    end
    photo_sensor = ei.b.photo_sensor_f(ei.b.photo_sensor_f>onsets(1) & ei.b.photo_sensor_f<offsets(end));
    dists = ei.b.dist(photo_sensor);
    diff_dists = diff(dists);
    inds = find(diff_dists < 100);
    temp_photo_sensor = photo_sensor;
    temp_photo_sensor(inds+1) = [];
    markersOn = temp_photo_sensor(1:(length(temp_photo_sensor)-1));
    markersOff = temp_photo_sensor(2:end);
end

if strcmp(lower(markerType),'belti')
    offsets = ei.b.photo_sensor_f(end) + round(1e6 * 0.3/ei.b.si);
    onsets = ei.b.photo_sensor_f(1) - round(1e6 * 0.3/ei.b.si);
    photo_sensor = ei.b.photo_sensor_f(ei.b.photo_sensor_f>onsets(1) & ei.b.photo_sensor_f<offsets(end));
    dists = ei.b.dist(photo_sensor);
    diff_dists = diff(dists);
    inds = find(diff_dists < 100);
    temp_photo_sensor = photo_sensor;
    temp_photo_sensor(inds+1) = [];
    markersOn = temp_photo_sensor(1:(length(temp_photo_sensor)-1));
    markersOff = temp_photo_sensor(2:end);
end

if strcmp(lower(markerType),'aironsets22')
    onsets = ei.b.air_puff_r(trials);
    offsets = ei.b.air_puff_f(trials);
    timeBefore = 2; % sec
    timeAfter = 2; % sec
    markersOn = onsets - round(1e6 * timeBefore/ei.b.si);
    markersOff = onsets + round(1e6 * timeAfter/ei.b.si);
    inds = find(markersOn < 0);
    markersOn(inds) = [];
    markersOff(inds) = [];
end

if strcmp(lower(markerType),'airoffsets22')
    onsets = ei.b.air_puff_r(trials);
    offsets = ei.b.air_puff_f(trials);
    timeBefore = 2; % sec
    timeAfter = 2; % sec
    markersOn = offsets - round(1e6 * timeBefore/ei.b.si);
    markersOff = offsets + round(1e6 * timeAfter/ei.b.si);
    inds = find(markersOn < 0);
    markersOn(inds) = [];
    markersOff(inds) = [];
end


if strcmp(lower(markerType),'aironsets27')
    onsets = ei.b.air_puff_r(trials);
    offsets = ei.b.air_puff_f(trials);
    timeBefore = 2; % sec
    timeAfter = 7; % sec
    markersOn = onsets - round(1e6 * timeBefore/ei.b.si);
    markersOff = onsets + round(1e6 * timeAfter/ei.b.si);
    inds = find(markersOn < 0);
    markersOn(inds) = [];
    markersOff(inds) = [];
end

if strcmp(lower(markerType),'airoffsets27')
    onsets = ei.b.air_puff_r(trials);
    offsets = ei.b.air_puff_f(trials);
    timeBefore = 2; % sec
    timeAfter = 7; % sec
    markersOn = offsets - round(1e6 * timeBefore/ei.b.si);
    markersOff = offsets + round(1e6 * timeAfter/ei.b.si);
    inds = find(markersOn < 0);
    markersOn(inds) = [];
    markersOff(inds) = [];
end


if strcmp(lower(markerType),'aironsets11')
    onsets = ei.b.air_puff_r(trials);
    offsets = ei.b.air_puff_f(trials);
    timeBefore = 1; % sec
    timeAfter = 1; % sec
    markersOn = onsets - round(1e6 * timeBefore/ei.b.si);
    markersOff = onsets + round(1e6 * timeAfter/ei.b.si);
    inds = find(markersOn < 0);
    markersOn(inds) = [];
    markersOff(inds) = [];
end

if strcmp(lower(markerType),'airoffsets11')
    onsets = ei.b.air_puff_r(trials);
    offsets = ei.b.air_puff_f(trials);
    timeBefore = 1; % sec
    timeAfter = 1; % sec
    markersOn = offsets - round(1e6 * timeBefore/ei.b.si);
    markersOff = offsets + round(1e6 * timeAfter/ei.b.si);
    inds = find(markersOn < 0);
    markersOn(inds) = [];
    markersOff(inds) = [];
end

if strcmp(lower(markerType),'aironsets01')
    onsets = ei.b.air_puff_r(trials);
    offsets = ei.b.air_puff_f(trials);
    timeBefore = 0; % sec
    timeAfter = 1; % sec
    markersOn = onsets - round(1e6 * timeBefore/ei.b.si);
    markersOff = onsets + round(1e6 * timeAfter/ei.b.si);
    inds = find(markersOn < 0);
    markersOn(inds) = [];
    markersOff(inds) = [];
end

if strcmp(lower(markerType),'airoffsets01')
    onsets = ei.b.air_puff_r(trials);
    offsets = ei.b.air_puff_f(trials);
    timeBefore = 0; % sec
    timeAfter = 1; % sec
    markersOn = offsets - round(1e6 * timeBefore/ei.b.si);
    markersOff = offsets + round(1e6 * timeAfter/ei.b.si);
    inds = find(markersOn < 0);
    markersOn(inds) = [];
    markersOff(inds) = [];
end

if strcmp(lower(markerType),'aironsets010')
    onsets = ei.b.air_puff_r(trials);
    offsets = ei.b.air_puff_f(trials);
    timeBefore = 0; % sec
    timeAfter = 10; % sec
    markersOn = onsets - round(1e6 * timeBefore/ei.b.si);
    markersOff = onsets + round(1e6 * timeAfter/ei.b.si);
    inds = find(markersOn < 0);
    markersOn(inds) = [];
    markersOff(inds) = [];
end

if strcmp(lower(markerType),'airoffsets010')
    onsets = ei.b.air_puff_r(trials);
    offsets = ei.b.air_puff_f(trials);
    timeBefore = 0; % sec
    timeAfter = 10; % sec
    markersOn = offsets - round(1e6 * timeBefore/ei.b.si);
    markersOff = offsets + round(1e6 * timeAfter/ei.b.si);
    inds = find(markersOn < 0);
    markersOn(inds) = [];
    markersOff(inds) = [];
end


if strcmp(lower(markerType),'motiononsets22') & ~isnan(trials)
    spSig = ei.b.fSpeed;
    onsets = ei.b.air_puff_r - round(1e6 * 1/ei.b.si);
    offsets = ei.b.air_puff_f + round(1e6 * 1/ei.b.si);
    mOn = pMarkers.markersOn;
    mOff = pMarkers.markersOff;
%     diffs = mOn - mOff;
%     inds = find(diffs < 10000);
%     mOn(inds) = [];
%     mOff(inds) = [];
    all{1} = 1:(onsets(1)-1);
    for ii = 1:(length(onsets)-1)
        all{ii+1} = (offsets(ii)+1):(onsets(ii+1)-1);
    end
    all{ii+2} = (offsets(ii+1)+1):length(spSig);
    valsOn = [];
    valsOff = [];
    for ii = 1:length(all)
        vals = all{ii};
        motionOnsets = intersect(vals,mOn);
        motionOffsets = intersect(vals,mOff);
        if isempty(motionOnsets) | isempty(motionOffsets)
            continue;
        end
        valsOn = [valsOn motionOnsets];
        valsOff = [valsOff motionOffsets];
    end
    timeBefore = 2; % sec
    timeAfter = 2; % sec
    markersOn = valsOn - round(1e6 * timeBefore/ei.b.si);
    markersOff = valsOn + round(1e6 * timeAfter/ei.b.si);
    inds = find(markersOn < 0);
    markersOn(inds) = [];
    markersOff(inds) = [];
    n=0;
end

if strcmp(lower(markerType),'motionoffsets22') & ~isnan(trials)
    spSig = ei.b.fSpeed;
    onsets = ei.b.air_puff_r - round(1e6 * 1/ei.b.si);
    offsets = ei.b.air_puff_f + round(1e6 * 1/ei.b.si);;
    mOn = pMarkers.markersOn;
    mOff = pMarkers.markersOff;
%     diffs = mOn - mOff;
%     inds = find(diffs < 10000);
%     mOn(inds) = [];
%     mOff(inds) = [];
    all{1} = 1:(onsets(1)-1);
    for ii = 1:(length(onsets)-1)
        all{ii+1} = (offsets(ii)+1):(onsets(ii+1)-1);
    end
    all{ii+2} = (offsets(ii+1)+1):length(spSig);
    valsOn = [];
    valsOff = [];
    for ii = 1:length(all)
        vals = all{ii};
        motionOnsets = intersect(vals,mOn);
        motionOffsets = intersect(vals,mOff);
        if isempty(motionOnsets) | isempty(motionOffsets)
            continue;
        end
        valsOn = [valsOn motionOnsets];
        valsOff = [valsOff motionOffsets];
    end
    timeBefore = 2; % sec
    timeAfter = 2; % sec
    markersOn = valsOff - round(1e6 * timeBefore/ei.b.si);
    markersOff = valsOff + round(1e6 * timeAfter/ei.b.si);
    inds = find(markersOn < 0);
    markersOn(inds) = [];
    markersOff(inds) = [];
    n=0;
end

if strcmp(lower(markerType),'motionoffsetaironset')
    spSig = ei.b.fSpeed;
    onsets = ei.b.air_puff_r(trials);
    mOff = pMarkers.markersOff;
    markersOn = mOff(1:(end-1)) + round(1e6 * 1/ei.b.si);
    markersOff = onsets(2:end) - round(1e6 * 0.25/ei.b.si);
    [integrity pl] = checkIntegrity(markersOn,markersOff,[]);
    if integrity < 0
        error
    end
end
if strcmp(lower(markerType),'motiononsetsoffsets')
    spSig = ei.b.fSpeed;
    onsets = ei.b.air_puff_r(trials) - round(1e6 * 0.25/ei.b.si);
    offsets = ei.b.air_puff_f(trials);
    mOn = pMarkers.markersOn;
    mOff = pMarkers.markersOff;
    inds1 = [];
    for ii = 1:length(onsets)
        if ii == 10
            n = 0;
        end
        inds = find(mOn > onsets(ii),1,'first');
        if isempty(inds)
            break;
        end
        markersOn(ii) = mOn(inds);
        inds = find(mOff > offsets(ii),1,'first');
        if isempty(inds)
            markersOn(end) = [];
            break;
        end
        markersOff(ii) = mOff(inds);
        if markersOn(ii) > markersOff(ii)
            inds1 = [inds1 ii];
        end
    end
    markersOn(inds1) = [];
    markersOff(inds1) = [];
    [integrity pl] = checkIntegrity(markersOn,markersOff,spSig);
    if integrity < 0
        error
    end
end

if strcmp(lower(markerType),'motioni') & isnan(trials)
    spSig = ei.b.fSpeed;
    try
        onsets = ei.b.air_puff_r;
        offsets = ei.b.air_puff_f;
    catch
        psL = length(ei.b.photo_sensor_f);
        onsets = ei.b.photo_sensor_f(1:(psL-1));
        offsets = ei.b.photo_sensor_f(2:psL);
    end
    mOn = pMarkers.markersOn;
    mOff = pMarkers.markersOff;
    all{1} = 1:(onsets(1)-1);
    for ii = 1:(length(onsets)-1)
        all{ii+1} = (offsets(ii)+1):(onsets(ii+1)-1);
    end
    all{ii+2} = (offsets(ii+1)+1):length(spSig);
    valsOn = [];
    valsOff = [];
    for ii = 1:length(all)
        vals = all{ii};
        motionOnsets = intersect(vals,mOn);
        motionOffsets = intersect(vals,mOff);
        if ~isempty(motionOnsets) & ~isempty(motionOffsets)
            if motionOnsets(1) > motionOffsets(1)
                motionOffsets(1) = [];
            end
        end
        if isempty(motionOnsets) | isempty(motionOffsets)
            continue;
        end
        if (length(motionOnsets) ~= length(motionOffsets))
            [motionOnsets,motionOffsets] = resolveMotionOnsetsOffsets(ei,motionOnsets,motionOffsets,spSig);
        end
        valsOn = [valsOn motionOnsets];
        valsOff = [valsOff motionOffsets];
        n = 0;
    end
    motionOnsets = valsOn; motionOffsets = valsOff;
    [integrity pl] = checkIntegrity(motionOnsets,motionOffsets,spSig);
    if integrity < 0
        error
    end
    
    diffs = motionOnsets(2:end) - motionOffsets(1:(end-1));
    inds = find(diffs < 10000);
    motionOnsets(inds+1) = [];
    motionOffsets(inds) = [];
    
    diffs = motionOffsets - motionOnsets;
    inds = find(diffs < 10000);
    motionOnsets(inds) = [];
    motionOffsets(inds) = [];
    markersOn = motionOnsets;
    markersOff = motionOffsets;
    n = 0;
end

if strcmp(lower(markerType),'motion') & isnan(trials)
    spSig = ei.b.fSpeed;
    motionOnsets = find_rising_edge(spSig > 0.25,0.1,500); % find places where speed is above threshold
    motionOffsets = find_falling_edge(spSig > 0.25,-0.1,500); % find places where speed is below threshold
    if isempty(motionOnsets)
        markersOn = [];
        markerOff = [];
    else
        [motionOnsets,motionOffsets] = resolveMotionOnsetsOffsets(ei,motionOnsets,motionOffsets,spSig);
        diffs = motionOnsets(2:end) - motionOffsets(1:(end-1));
        inds = find(diffs < 500);
        motionOnsets(inds+1) = [];
        motionOffsets(inds) = [];
        markersOn = motionOnsets;
        markersOff = motionOffsets;
    end
end


if sum(markersOn < 0) > 0 | sum(markersOff < 0) > 0
    n = 0;
end
b = ei.b;
try
inds = [];
for ii = 1:length(markersOn)
    st = markersOn(ii);
    se = markersOff(ii);
    frames = (find(b.frames_f >= st & b.frames_f <= se));
    if isempty(frames)
        inds = [inds ii];
    end
end
markersOn(inds) = [];
markersOff(inds) = [];
catch
end
if length(markersOn) ~= length(markersOff)
    if length(markersOn) > length(markersOff)
        markersOn(end) = [];
    else
        markersOff(end) = [];
    end
    [integrity pl] = checkIntegrity(markersOn,markersOff,spSig);
    if integrity < 0
        error
    end
    if length(markersOn) ~= length(markersOff)
        error;
    end
    n = 0;
end
        


function [motionOnsets,motionOffsets] = resolveMotionOnsetsOffsets(ei,motionOnsets,motionOffsets,spSig)
while 1 % see if there is a first motion offset which is less than the first motion onset
    if motionOffsets(1) < motionOnsets(1)
        motionOffsets(1) = [];
    else
        break;
    end
end
while 1 % see if there is an isolated last motion onset which is greater than the last motion onset
    if motionOnsets(end) > motionOffsets(end)
        motionOnsets(end) = [];
    else
        break;
    end
end
motionOffsets1 = [];
motionOnsets1 = [];
ind = 1;
for ii = 1:length(motionOnsets) % process each motion onset and find its paired motion offset
    st = motionOnsets(ii);
    if ii > 1
        if st < motionOffsets1(end)
            n = 0;
        end
    end
    inds = find(motionOffsets > st);
    indsOn = find(motionOnsets > st);
    if isempty(inds)
        break;
    end
    if isempty(indsOn)
        break;
    end
    if motionOnsets(indsOn(1)) < motionOffsets(inds(1))
        continue;
    end
    if length(inds) > 5
        inds = inds(1:5);
    end
    for jj = 1:length(inds)
        speeds = spSig(st:motionOffsets(inds(jj)));
        temp(jj) = sum(speeds<0.1);
    end
    mInd = find(temp == min(temp));
    motionOffsets1(ind) = motionOffsets(inds(mInd(1)));
    motionOnsets1(ind) = st;
    ind = ind + 1;
end
motionOffsets = motionOffsets1;
motionOnsets = motionOnsets1;
[integrity pl] = checkIntegrity(motionOnsets,motionOffsets,spSig);
if integrity < 0
    error
end

% integriy check ... that is see if every motion onset has a corresponding motion offset
function [integrity ii] = checkIntegrity(motionOnsets,motionOffsets,spSig)
integrity = 1;
for ii = 1:(length(motionOnsets)-1)
    st = motionOnsets(ii); se = motionOnsets(ii+1);
    inds = find(motionOffsets > st & motionOffsets < se);
    if length(inds) > 1
        integrity = -1;
        return;
    end
end

for ii = 1:(length(motionOffsets)-1)
    st = motionOffsets(ii); se = motionOffsets(ii+1);
    inds = find(motionOnsets > st & motionOnsets < se);
    if length(inds) > 1
        integrity = -2;
        return;
    end
end

function edge = find_rising_edge(signal,threshold,minimum_diff)
edge = find(diff(signal) >= threshold);
temp = find(diff(edge)<=minimum_diff);
edge(temp+1) = [];

function edge = find_falling_edge(signal,threshold,minimum_diff)
edge = find(diff(signal) <= threshold);
temp = find(diff(edge)<=minimum_diff);
edge(temp+1) = [];

function plotMotionAll(ei,motionOnsets,motionOffsets)
% check
for ii = 1:length(motionOnsets) 
    figure(101);clf
    plot(ei.b.ts,ei.b.fSpeed);hold on;
    tMOn = zeros(size(ei.b.fSpeed));
    tMOn(motionOnsets(ii)) = 15;
    plot(ei.b.ts,tMOn,'m');
    tMOn = zeros(size(ei.b.fSpeed));
    tMOn(motionOffsets(ii)) = 15;
    plot(ei.b.ts,tMOn,'c');
    plot(ei.b.ts,ei.b.air_puff_raw,'k');
    xlim([ei.b.ts(motionOnsets(ii))-15 ei.b.ts(motionOffsets(ii))+15]);
    pause(0.3);
end