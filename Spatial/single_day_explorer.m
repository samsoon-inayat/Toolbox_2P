function place_field_analysis_1

ei = evalin('base','ei{5}');
b = ei.b;
b.sTrials = 1:43
pp = 1
signals = ei.plane{pp}.tP.signals;
b.frames_f = ei.plane{pp}.b.frames_f;

% S = ei.S;
% ccsi = find_place_cells(S,'adjRSquare_threshold',0.7);
ccsi = 1:length(ei.plane{pp}.tP.areCells);
ccs = find(ei.plane{pp}.tP.areCells(ccsi));

caSigAll = ei.plane{pp}.tP.deconv.caSigAll;
spSigAll = ei.plane{pp}.tP.deconv.spSigAll;

nCols = 3;
nRows = 3;
totalPlots = nRows * nCols;
spind = reshape(1:totalPlots,nCols,nRows)';

startI = 1:totalPlots:length(ccs);
endI = [totalPlots:totalPlots:length(ccs),length(ccs)];

figure(101);clf;
for ccc = 1:length(startI)
    thisInds = startI(ccc):endI(ccc);
    for row =1:nRows
        for  col = 1:nCols
            dd = sub2ind([nRows nCols],row,col);
            cc = thisInds(dd);
            thisSignal = signals(ccs(cc),:);
            ca_sig = caSigAll{ccsi(cc)}';
            tsp = spSigAll{ccsi(cc)}';
%             tsp = thisSignal;
            titleText = sprintf('Animal - %s - %s - Rec_%d, Cell %d of %d',ei.db(pp).mouse_name,ei.db(pp).date,pp,cc,length(ccs));
            [xValsA,ccSignalA,trials,mInfo] = getTrialSignalsDist(tsp,b,b.air_puff_r(b.sTrials),b.air_puff_f(b.sTrials)+100);
            minSig = min(ccSignalA(:)); maxSig = max(ccSignalA(:));
            subplot(nRows,nCols,spind(row,col));
            mSig = mean(ccSignalA);
            f = fit(xValsA',mSig','gauss2');
            mSigFF = (feval(f,xValsA))';
            cla
            plot(xValsA,mean(ccSignalA));hold on;
            plot(xValsA,mSigFF);
            xlim([min(xValsA) max(xValsA)]);
%             plot([0 0],ylim);
%             imagesc(ccSignalA,[minSig maxSig]);colorbar;
            title(cc);
    %         figure(102);clf;
    %         plot(mean(ccSignalA));
        end
    end
    n = 0;
end

