function place_field_analysis_2_for_172110_Feb_07_2018_Recording

ei = evalin('base','ei{1}');
GF = gausswin(3);
b = ei.b;
signals = ei.tP.signals;

% S = ei.S;
% ccsi = find_place_cells(S,'adjRSquare_threshold',0.8);
ccsi = 1:length(ei.areCells);
ccs = ei.areCells(ccsi);

min_pc_width = 2;
max_pc_width = 120;
cm_per_bin = 150/100;

spSigAll = ei.deconv.spSigAll;
caSigAll = ei.deconv.caSigAll;
coi = [];
coiSI = [];
allSIs =[];
% load('temp.mat');
% labelCells(ei,coiSI,allSIs);
% OF_frames = OF_caSig(ei,ccsi(coi));
% for cc = 1:length(ccsi)

cc = 1;
while cc<=length(ccsi)
%     if ~ismember(coi,cc)
%         continue;
%     end
%     if ~ismember(coiSI,cc)
%         continue;
%     end
    tsp = spSigAll{ccsi(cc)}';
    caSig = signals(ccs(cc),:)';
    sCaSig = caSigAll{ccsi(cc)}';
    A_cueless = b.trials(2:12); 
    A_Trials = b.trials(14:24);
    AP_Trials = b.trials(28:38);
    
    onsets = b.air_puff_r(A_cueless);
    offsets = b.air_puff_f(A_cueless);
    ac_onsets = onsets;
    ac_offsets = offsets;
%     tSignals = plotTrials(b,caSig,tsp,onsets,offsets);
%     tSignals = plotInterTrials(ei,caSig);
    AC = getDistRaster_1(b,caSig,tsp,onsets,offsets,0);
    plotDistRaster(105,AC);
%     
%     onsets = b.air_puff_r(A_Trials);
%     offsets = b.air_puff_f(A_Trials);
%     a_onsets = onsets;
%     a_offsets = offsets;
% %     tSignals = plotTrials(b,caSig,tsp,onsets,offsets);
% %     tSignals = plotInterTrials(ei,caSig);
%     A = getDistRaster_1(b,caSig,tsp,onsets,offsets,0);
%     
%     onsets = b.air_puff_r([A_Trials,AP_Trials]);
%     offsets = b.air_puff_f([A_Trials,AP_Trials]);
% %     plotMarkers(ei,onsets,offsets,100);
%     Ar = getTimeRaster_simple(b,caSig,tsp,onsets-round(1/(ei.b.si*1e-6)),onsets+round(1/(ei.b.si*1e-6)),0);
%     Af = getTimeRaster_simple(b,caSig,tsp,offsets-round(1/(ei.b.si*1e-6)),offsets+round(1/(ei.b.si*1e-6)),0);
%     
%     photo_sensor = b.photo_sensor_f(b.photo_sensor_f>a_onsets(1) & b.photo_sensor_f<a_offsets(end));
%     diff_photo_sensor = diff(photo_sensor);
%     inds = find(diff_photo_sensor < 10000);
%     temp_photo_sensor = photo_sensor;
%     temp_photo_sensor(inds) = [];
%     
%     onsets = temp_photo_sensor(1:(length(temp_photo_sensor)-1));
%     offsets = temp_photo_sensor(2:end);
% %      plotMarkers(ei,onsets,offsets,100);
%     AP = getDistRaster_1(b,caSig,tsp,onsets,offsets,0);
%     
%     
% %     AP_Trials = b.trials(21:40);
%     ap_onsets = b.air_puff_r(AP_Trials);
%     ap_offsets = b.air_puff_f(AP_Trials);
%     
%     photo_sensor = b.photo_sensor_f(b.photo_sensor_f>ap_onsets(1) & b.photo_sensor_f<ap_offsets(end));
%     diff_photo_sensor = diff(photo_sensor);
%     inds = find(diff_photo_sensor < 100000);
%     temp_photo_sensor = photo_sensor;
%     temp_photo_sensor(inds) = [];
%     plotMarkers(ei,ap_onsets,ap_offsets,100);
%     onsets = temp_photo_sensor(1:(length(temp_photo_sensor)-1));
%     offsets = temp_photo_sensor(2:end);
%     P = getDistRaster_1(b,caSig,tsp,ap_onsets,ap_offsets,0);
%     
% %     Pr = getTimeRaster_simple(b,caSig,tsp,onsets-round(2/(ei.b.si*1e-6)),onsets+round(2/(ei.b.si*1e-6)),0);
%     
%     mFR = mean(tsp);
%     ccSignalA = A.distSigRaster./A.distDurRaster;
%     xValsA = A.dists;
%     
%     onsets = b.air_puff_f(A_Trials(1:(end-1)));
%     offsets = b.air_puff_r(A_Trials(2:end));
%     AI = getTimeRaster_1(b,caSig,tsp,onsets,offsets,0);
%     
% %     starting = find(b.air_puff_f>temp_photo_sensor(1),1,'first');
% %     onsets = b.air_puff_f(b.sTrials(starting:(end-1)));
% %     offsets = b.air_puff_r(b.sTrials((starting+1):end));
%     onsets = b.air_puff_f(AP_Trials(1:(end-1)));
%     offsets = b.air_puff_r(AP_Trials(2:end));
%     PI = getTimeRaster_1(b,caSig,tsp,onsets,offsets,0);
% 
%     plotTwoFrameRefs(A,AP,P,AI,PI);
%     figure(104);clf;
%     subplot 211;     imagesc(Ar.sigRaster);title('Air start');
%     subplot 212;     imagesc(Af.sigRaster);title('Air stop');
    cc
    cc = keyboardInput(cc,[1 length(ccsi)],[1 5],'');
    if cc < 0
        break;
    end
end
coi


function SI = spatial_information(lambdaX,deltaT)
lambdaX = lambdaX/max(lambdaX);
pxN = deltaT/sum(deltaT);
mR = sum(lambdaX.*pxN);
SI = nansum(lambdaX.*log2(lambdaX/mR).*pxN)/mR;


function PC = find_place_cells_Rules (mSig,ccSignalA,varargin)

p = inputParser;
default_cm_per_bin = 157/100;
default_min_pc_width = 5;
default_max_pc_width = 120;

addRequired(p,'mSig',@isnumeric);
addRequired(p,'ccSignalA',@isnumeric);
addOptional(p,'cm_per_bin',default_cm_per_bin,@isnumeric);
addOptional(p,'max_pc_width',default_max_pc_width,@isnumeric);
addOptional(p,'min_pc_width',default_min_pc_width,@isnumeric);
parse(p,mSig,ccSignalA,varargin{:});
cm_per_bin = p.Results.cm_per_bin;
min_pc_width = p.Results.min_pc_width;
max_pc_width = p.Results.max_pc_width;

PC = 0;
initial_threshold = 0.3*(max(mSig) - min(mSig));%+min(mSig);
idxs = find(mSig > initial_threshold);
iidxs = diff(idxs);
gOnes = find(iidxs>1);
if ~isempty(gOnes)
    if gOnes(1) ~= 1
        gOnes = [1 gOnes];
    end
    if gOnes(end) ~= length(iidxs)
        gOnes = [gOnes length(iidxs)];
    end

    d_gOnes = diff(gOnes)*cm_per_bin;
    if sum(d_gOnes>min_pc_width & d_gOnes<max_pc_width) == 0
        return;
    end
    longest_chunk = find(d_gOnes>min_pc_width & d_gOnes<max_pc_width);
    if length(longest_chunk) > 1
        longest_chunk = find(d_gOnes == max(d_gOnes(longest_chunk)));
    end
    idxs_of_longest_chunk = idxs((gOnes(longest_chunk)+2):gOnes(longest_chunk+1));
else
    idxs_of_longest_chunk = idxs;
end

if length(idxs_of_longest_chunk)*cm_per_bin > max_pc_width
    return;
end
if length(idxs_of_longest_chunk)*cm_per_bin < min_pc_width
    return;
end
mean_activity_inside_place_field = mean(mSig(idxs_of_longest_chunk));
mean_activity_outside_place_field = mean(mSig(mSig < initial_threshold));
if mean_activity_inside_place_field < (3*mean_activity_outside_place_field)
    return;
end
[~,idx_of_peak_of_trials] = max(ccSignalA,[],2);
count = 0;
for kk = 1:length(idx_of_peak_of_trials)
    count = count + sum(idxs_of_longest_chunk == idx_of_peak_of_trials(kk));
end
if size(ccSignalA,2)/2 < count
    return;
end
PC = 1;



