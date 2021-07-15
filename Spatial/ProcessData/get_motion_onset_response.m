function ei_C = get_motion_onset_response(ei_C,owr)
for ii = 1:length(ei_C)
    ei_C{ii} = process_one_animal(ei_C{ii},1,owr);
end

function ei = process_one_animal(ei,thr,owr)
[mOnst,mOnse,mOffst,mOffse] = get_monoffs(ei,thr);
% disp(sprintf('%d,%d - %d,%d',length(mOnst),length(mOnse),length(mOffst),length(mOffse)));
binwidths = evalin('base','binwidths');
for pp = 1:length(ei.plane)
    ei.tP = ei.plane{pp}.tP;
    thispFolder = ei.plane{pp}.folder;
    ei.deconv = ei.plane{pp}.tP.deconv;
    rasters = make_rasters(ei,pp,mOnst,mOnse,'time',binwidths);
    trials = 1:size(rasters.sp_rasters1,1);
    rasters = findRasterProperties_1(thispFolder,0,'motionOnsets11T',rasters,'time',trials,owr);
    ei.plane{pp}.motionOnset_rasters = rasters;
    ei.plane{pp}.motionOnset_rasters.activity_speed_corr = find_activity_speed_corr(ei.deconv.spSigAll,ei.b.fSpeed(ei.plane{pp}.b.frames_f));

    rasters = make_rasters(ei,pp,mOffst,mOffse,'time',binwidths);
    trials = 1:size(rasters.sp_rasters1,1);
    rasters = findRasterProperties_1(thispFolder,0,'motionOffsets11T',rasters,'time',trials,owr);
    ei.plane{pp}.motionOffset_rasters = rasters;
    ei.plane{pp}.motionOffset_rasters.activity_speed_corr = ei.plane{pp}.motionOnset_rasters.activity_speed_corr;
end

function out = find_activity_speed_corr(activity,speed)
speed = speed(1:size(activity,2));
out = corr(activity',speed');


function out = find_speed_response(activity,speed)


function [mOnst,mOnse,mOffst,mOffse] = get_monoffs(ei,thr)
spSig = ei.b.fSpeed;
motionOnsets = find_rising_edge(spSig > thr,0.1,500);
motionOffsets = find_falling_edge(spSig > thr,-0.1,500); % find places where speed is below threshold
[mOnst,mOnse] = get_markers_motion(ei,motionOnsets,motionOffsets);
[mOffst,mOffse] = get_markers_motion(ei,motionOffsets,motionOnsets);
disp(sprintf('%d,%d - %d,%d',length(mOnst),length(mOnse),length(mOffst),length(mOffse)));
% diplay_with_air_puff(ei.b,mOnst,mOnse);
% n = 0;
% diplay_with_air_puff(ei.b,mOffst,mOffse);
% n = 0;

 

function [st2,se2] = get_markers_motion(ei,motionOnsets,motionOffsets)
thr_t = 0.5;
tb = 1.5;
b = ei.b;
st = motionOnsets - round(1e6 * tb/b.si);
se = motionOnsets + round(1e6 * tb/b.si);
inds = find(st<0);
if ~isempty(inds)
    st(inds) = [];
    se(inds) = [];
end
[st1,se1] = removeOverlapWithAirPuff(b,st,se);
[se2,st2] = removeOverlapWithAirPuff(b,se1,st1);
[st2,se2] = eliminate_close_ones(b,b.air_puff_r,st2,se2,thr_t); 
[st2,se2] = eliminate_close_ones(b,b.air_puff_f,st2,se2,thr_t); 
[st2,se2] = eliminate_close_ones(b,motionOffsets,st2,se2,thr_t); 

[st2,se2] = check_for_frames(ei,st2,se2);

function [mon,moff] = check_for_frames(ei,mon,moff)
for pp = 1:length(ei.plane)
    firstFrame(pp) = ei.plane{pp}.b.frames_f(1);
    lastFrame(pp) = ei.plane{pp}.b.frames_f(end);
end
mlastFrame = min(lastFrame);
MfirstFrame = max(firstFrame);
inds = find(mon < MfirstFrame); mon(inds) = []; moff(inds) = [];
inds = find(moff > mlastFrame); mon(inds) = []; moff(inds) = [];


function [st,se] = eliminate_close_ones(b,aon,st,se,thr)
[st,se] = check_close_to_air_puff(b,aon,st,se,thr); 
[se,st] = check_close_to_air_puff(b,aon,se,st,thr);

function [st,se] = check_close_to_air_puff(b,aon,st,se,thr)
inds = [];
for ii = 1:length(st)
    d = b.ts(aon) - b.ts(st(ii));
    dlz = abs(d(d<0));
    dgz = d(d>0);
    if sum(min(dlz) < thr)>0 || sum(min(dgz) < thr)>0
        inds = [inds ii];
    end
end
if ~isempty(inds)
    st(inds) = []; se(inds) = [];
end


function [st1,se1] = removeOverlapWithAirPuff(b,st,se)
dar = b.air_puff_raw < 0.5;
inds = find(st > length(dar));
st(inds) = [];
se(inds) = [];
zm = zeros(size(b.air_puff_raw));
zm(st) = 1; 
if length(dar) ~= length(zm)
    n = 0;
end
zm1 = dar.*zm;
st1 = find(zm1);
inds = [];
for ii = 1:length(st1)
    inds(ii) = find(st == st1(ii));
end
se1 = se(inds);

