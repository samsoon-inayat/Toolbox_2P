function place_field_analysis

ei = evalin('base','ei{1}');
GF = gausswin(6); 
b = ei.b;
signals = ei.tP.signals;

S = ei.S;
ccsi = find_place_cells(S,'adjRSquare_threshold',0.7);
ccs = ei.areCells(ccsi);

% ccs = ei.areCells;
min_pc_width = 15;
max_pc_width = 120;
cm_per_bin = 157/100;
coi = [];
coiSI = [];
rng(0,'twister');
display_figs = 1;
display_figs_raw = 1;
caSigAll = ei.deconv.caSigAll;
spSigAll = ei.deconv.spSigAll;

for cc = 1:length(ccs)
    thisSignal = signals(ccs(cc),:);
    ca_sig = caSigAll{ccsi(cc)}';
    tsp = spSigAll{ccsi(cc)}';
    titleText = sprintf('Animal - %s - %s - Rec_%d, Cell %d of %d',ei.animal_id,ei.exp_date,ei.recording_number,cc,length(ccs));
    if display_figs_raw
        figure(100);clf;
        ts = ei.b.ts(ei.b.frames_f);
        stem(ts,tsp,'.');hold on;
        plot(ts,ca_sig,'r');
        plot(ei.b.ts,(ei.b.air_puff_raw)*0.5*max(tsp(:)),'color','k','linewidth',1.5);
        plot(ei.b.ts,ei.b.speed+max(ca_sig(:))+10,'m');
        plot(ei.b.ts,(ei.b.air_puff_raw)*0.5*max(tsp(:))+min(min(ei.b.speed+max(ca_sig(:))+10)),'color','k','linewidth',1.5);
        xlim([ts(1) ts(end)]);
        for ii = ei.b.trials
            text(ei.b.ts(ei.b.air_puff_r(ii)),max(ca_sig(:))+5,num2str(ii));
        end
    end
    [xValsA,ccSignalA,trials,mInfo] = getTrialSignalsDist(tsp,b,b.air_puff_r(b.trials),b.air_puff_f(b.trials)+100);
    minSig = min(ccSignalA(:)); maxSig = max(ccSignalA(:));
    figure(101);clf;
    imagesc(ccSignalA,[minSig maxSig]);colorbar;
    figure(102);clf;
    plot(mean(ccSignalA));
%     set(gca,'YTicks',trials);
    n = 0;
    continue;
    mSig = mean(ccSignalA);
    mDur = circshift(mean(mInfo.trialDurations),5);
    
    mSig = filter(GF,1,mSig);
    f = fit(xValsA',mSig','gauss2');
%     mSig = (feval(f,xValsA))';
    
    initial_threshold = 0.3*(max(mSig) - min(mSig));%+min(mSig);
    if display_figs
        ff = makeFigureRowsCols(110,[NaN NaN 13 5],'RowsCols',[2 1],'spaceRowsCols',[0.11 0],'rightUpShifts',[0.05 0.1],'widthHeightAdjustment',[-100 -150]);
        
        [xValsA,ccSignalA] = getTrialSignalsDist(ca_sig,b,b.air_puff_r(b.trials),b.air_puff_f(b.trials)+100);
        minSig = min(ccSignalA(:)); maxSig = max(ccSignalA(:));
        meanSig = mean(ccSignalA);
        axes(ff.h_axes(1,1));
        imagesc(ccSignalA,[minSig maxSig]);colorbar;
        
        [xValsAS,ccSignalAS] = getTrialSignalsDist(tsp,b,b.air_puff_r(b.trials),b.air_puff_f(b.trials)+100);
        minSigS = min(ccSignalAS(:)); maxSigS = max(ccSignalAS(:));
        meanSigS = mean(ccSignalAS);
        axes(ff.h_axes(2,1));
        imagesc(ccSignalAS,[minSigS maxSigS]);colorbar;

        
%         onsetsT = b.air_puff_r(b.trials)-floor((25)*1e6/b.si);
%         offsetsT = b.air_puff_r(b.trials);
%         [xValsAT,ccSignalAT] = getTrialSignalsTime(thisSignal,b,onsetsT,offsetsT);
%         minSigT = min(ccSignalAT(:)); maxSigT = max(ccSignalAT(:));
%         meanSigT = mean(ccSignalAT);
%         axes(ff.h_axes(2,1));
%         imagesc(ccSignalAT,[minSigT maxSigT]);colorbar;
        
        figure(300);clf;
%         subplot 211
        plot(xValsA,meanSigS);hold on;
        xlims = xlim;
        initial_threshold = 0.3*(max(meanSigS) - min(meanSigS));%+min(mSig);
        plot(xlims,[initial_threshold initial_threshold],'r');
%         plot(f,xValsA,mSig);
%         title(titleText);
%         subplot 212
%         plot(xValsA,mDur/sum(mDur));
    end
    n = 0;
    
    SI = spatial_information(mSig,mDur);
%     [H, SIm] = entropy_2D(mSig,mDur);
%     [SIm entropy fd_bins]= mutualinformationx(mSig,mean(mInfo.trialDurations),100)
%     SIc = condentropy(mSig,mean(mInfo.trialDurations));
    for ii = 1:1000
        mSigS = mSig(randperm(length(mSig)));
%         mDurS = mDur(randperm(length(mDur)));
%         [H, SIms(ii)] = entropy_2D(mSigS,mDurS);
%         [SIms(ii) entropy fd_bins]= mutualinformationx(mSigS,mean(mInfo.trialDurations),100);
        SIs(ii) = spatial_information(mSigS,mDur);
    end
    pSIs = prctile(SIs,95);
%     pSIms = prctile(SIms,95);
    [SI pSIs]
%     [SIm pSIms]
    if SI > pSIs
        coiSI = [coiSI cc];
    end
%     mSig = filter(GF,1,mSig);
    
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
            continue;
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
        continue;
    end
    if length(idxs_of_longest_chunk)*cm_per_bin < min_pc_width
        continue;
    end
    mean_activity_inside_place_field = mean(mSig(idxs_of_longest_chunk));
%     mean_activity_inside_place_field = mean(mSig(mSig > initial_threshold));
    mean_activity_outside_place_field = mean(mSig(mSig < initial_threshold));
    if mean_activity_inside_place_field < (3*mean_activity_outside_place_field)
        continue;
    end
    [~,idx_of_peak_of_trials] = max(ccSignalA,[],2);
    count = 0;
    for kk = 1:length(idx_of_peak_of_trials)
        count = count + sum(idxs_of_longest_chunk == idx_of_peak_of_trials(kk));
    end
    if size(ccSignalA,2)/3 < count
        continue;
    end
    coi = [coi cc];
end
coi
n = 0;

% coi =
% 
%      5     7    12    13    25    32    71    84   130   144   153   158   231   236   336   344   865


function SI = spatial_information(mSig,binDurations)
sumBD = sum(binDurations);
pi = binDurations./sumBD;
sumf = sum(mSig);
fif = mSig./sumf;
log2fif = log2(fif);
prod = pi.*fif.*log2fif;
SI = -nansum(prod,2);


function hf = lPlotTrials1(thisSignal,b,figNum,position,titleText)
pp.numberOfRows = 2;
pp.spaceBetweenRows = [0.05];
pp.rowHeights = [0.15 0.45];

pp.numberOfCols = [1 1];
pp.colWidths = {[0.9]; [0.9]};
pp.spaceBetweenCols = {[0];[0]};

pp.leftOffset = 0.07;
pp.bottomOffset = 0.15;

pp = getPanelsProps(pp);

[xValsA,ccSignalA] = getTrialSignalsDist(thisSignal,b,b.air_puff_r(b.trials),b.air_puff_f(b.trials)+100);

minSig = min(ccSignalA(:)); maxSig = max(ccSignalA(:));
hf = makeFigureWindow(figNum,position,1);

meanSig = mean(ccSignalA);

% air trials
pp = makeAxes(hf,pp,2,1,[0 0 0 0]);
imagesc(ccSignalA,[minSig maxSig]);colorbar;
axis off;
set(gca,'YTickLabel',[]);
% xTicksVals = get(gca,'XTick');
% set(gca,'XTickLabel',round(xValsA(1,xTicksVals)));
text(15,-1,titleText,'FontSize',10);
% title(titleText);

% mean of air trials
pp = makeAxes(hf,pp,1,1,[0 0 0 0]);
imagesc(meanSig/max(meanSig),[0 1]);colorbar;
set(gca,'YTickLabel',[]);
% xTicksVals = get(gca,'XTick');
% set(gca,'XTickLabel',round(xValsA(1,xTicksVals)));
xlabel('Distance (cm)');
colormap parula

