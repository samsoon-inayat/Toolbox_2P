function single_day_explorer_belt_photo_sensor

ei = evalin('base','ei{6}');
b = ei.b;
signals = ei.tP.signals;

% S = ei.S;
% ccsi = find_place_cells(S,'adjRSquare_threshold',0.7);
ccsi = 1:length(ei.areCells);
ccs = ei.areCells(ccsi);
caSigAll = ei.deconv.caSigAll;
spSigAll = ei.deconv.spSigAll;
psCells = [];
for cc = 1:length(ccsi)
    tsp = spSigAll{ccsi(cc)}';
    onsetT = (b.photo_sensor_f-floor((5)*1e6/b.si));
    offsetT = (b.photo_sensor_f+floor((5)*1e6/b.si));
    [xValsA,ccSignalA] = getTrialSignalsTime_photo_sensor(tsp,b,onsetT,offsetT);
    mSig = mean(ccSignalA);
%     f = fit(xValsA',mSig','gauss2');
%     mSig = (feval(f,xValsA))';
    dvdr = size(ccSignalA,2)/2;
%     [p,tbl,stats] = anova1(ccSignalA,[ones(1,dvdr) 2*ones(1,dvdr)],'off');
    mSig1 = mean(mSig(1:dvdr));std_mSig1 = std(mSig(1:dvdr));
    mSig2 = mean(mSig((dvdr+1):end));
    thresh = mSig1 + 1*std_mSig1;
    if mSig2 < thresh
        continue;
    end
    psCells = [psCells cc];
end

ccsi = 1:length(psCells);
ccs = psCells;

nCols = 3;
nRows = 3;
totalPlots = nRows * nCols;
spind = reshape(1:totalPlots,nCols,nRows)';

startI = 1:totalPlots:length(ccs);
endI = [totalPlots:totalPlots:length(ccs),length(ccs)];


for ccc = 1:length(startI)
    thisInds = startI(ccc):endI(ccc);
    figure(101);clf;
    for row =1:nRows
        for  col = 1:nCols
            dd = sub2ind([nRows nCols],row,col);
            cc = thisInds(dd);
            tsp = spSigAll{ccs(cc)}';
%             tsp = thisSignal;
%             titleText = sprintf('Animal - %s - %s - Rec_%d, Cell %d of %d',ei.animal_id,ei.exp_date,ei.recording_number,cc,length(ccs));
            onsetT = (b.photo_sensor_f-floor((5)*1e6/b.si));
            offsetT = (b.photo_sensor_f+floor((5)*1e6/b.si));
            [xValsA,ccSignalA] = getTrialSignalsTime_photo_sensor(tsp,b,onsetT,offsetT);
            mSig = mean(ccSignalA);
            subplot(nRows,nCols,spind(row,col));
            f = fit(xValsA',mSig','gauss2');
            mSigFF = (feval(f,xValsA))';
            cla
            plot(xValsA-5,mean(ccSignalA));hold on;
            plot(xValsA-5,mSigFF);
            plot([0 0],ylim);
%             imagesc(ccSignalA,[minSig maxSig]);colorbar;
            title(ccs(cc));
    %         figure(102);clf;
    %         plot(mean(ccSignalA));
        end
    end
    n = 0;
end

