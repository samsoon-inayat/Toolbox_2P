function speedFigure

aei = evalin('base','ei10');
ei = aei([1:4]);
mData = evalin('base','mData');

n = 0

speed = [];
iSpeed = [];
for ii = 1:4
    onsets = ei{ii}.plane{1}.contexts(1).markers.air_onsets;
    offsets = ei{ii}.plane{1}.contexts(1).markers.air_offsets;
    iOnsets = offsets(1:(end-1));
    iOffsets = onsets(2:end);
    for jj = 1:length(onsets)
        st = onsets(jj);
        se = offsets(jj);
        speed(ii,jj) = mean(ei{ii}.b.fSpeed(st:se));
    end
    for jj = 1:length(iOnsets)
        st = iOnsets(jj);
        se = iOffsets(jj);
        iSpeed(ii,jj) = mean(ei{ii}.b.fSpeed(st:se));
    end
end

[mVals(1) semVals(1)] = findMeanAndStandardError(mean(speed,2));
[mVals(2) semVals(2)] = findMeanAndStandardError(mean(iSpeed,2));


% max(speed,[],2)
% max(iSpeed,[],2)


[h,p] = ttest2(mean(speed,2),mean(iSpeed,2),'tail','right')

figure(101);clf;
set(gcf,'units','inches');
set(gcf,'Position',[10 4 1 1]);
plotBarsWithSigLines(mVals',semVals',[1 2],[h p],'colors',{'r','b'},'ySpacing',10,'maxY',45,'sigTestName','','sigAsteriskFontSize',10, ...
                    'barWidth',0.6,'sigLinesStartYFactor',0.05);
changePosition(gca,[0.05 0.09 -0.4 -0.06]);
xlim([0.5 2.5]);
h = ylabel('Speed (cm/sec)'); changePosition(h,[0.6 -5 0]);
set(gca,'Ydir','normal','XTickLabel',{'Trials','InterTrials'},'FontSize',6,'FontWeight','Bold','TickDir','out');
xtickangle(25);
text(0.75,35,{'t-test','(N = 4 Mice)'},'FontSize',5);
save_pdf(gcf,mData.pdf_folder,sprintf('speedAvg.pdf'),600);




