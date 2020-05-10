function populationSpatialInformationPlots_3C(tei)
% close all
thisCols_all = {[0 0 0],'b','r',[0 0.7 0.3],'m','c'};
if ~exist('tei','var')
    ei = evalin('base','ei');
    tei = ei([5 6]);
end

owsi = 0;

fileName = makeName('trialsRastersMetaData.mat',tei{1}.folders.thispFolder);
trialsRastersMD = load(fileName);
trials = trialsRastersMD.trials;
% trials = {[],[1:15],[16:31],[32:47]};
trl1 = trials{1}; trl2 = trials{2}; trl3 = trials{3}; trl4 = trials{4};
% trl1 = 1:10; trl3 = 11:26; trl4 = 27:42;
if ~isempty(trl1)
    nbb_air = getRasters_air(tei,trl1);
    temp = find_mutual_information(tei,nbb_air,owsi); nbb_air.SI = temp.zMI;
    nbb_belt = getRasters_belt(tei,trl1);
    temp = find_mutual_information(tei,nbb_belt,owsi); nbb_belt.SI = temp.zMI;
end
nbc_air = getRasters_air(tei,trl2);
temp = find_mutual_information(tei,nbc_air,owsi); nbc_air.SI = temp.zMI;
nbc_belt = getRasters_belt(tei,trl2);
temp = find_mutual_information(tei,nbc_belt,owsi); nbc_belt.SI = temp.zMI;
bc_air_belt = getRasters_air_belt(tei,trl3);
temp = find_mutual_information(tei,bc_air_belt,owsi); bc_air_belt.SI = temp.zMI;
bb_air_belt = getRasters_air_belt(tei,trl4);
temp = find_mutual_information(tei,bb_air_belt,owsi); bb_air_belt.SI = temp.zMI;

maxSI = max([nbc_air.SI nbc_belt.SI bc_air_belt.SI bb_air_belt.SI]);
minSI = min([nbc_air.SI nbc_belt.SI bc_air_belt.SI bb_air_belt.SI]);
bins = minSI:2:maxSI;
[bar1 xs] = hist(nbc_air.SI,bins); bar1 = 100*bar1/sum(bar1);
[bar2 xs] = hist(nbc_belt.SI,bins); bar2 = 100*bar2/sum(bar2);
[bar3 xs] = hist(bc_air_belt.SI,bins); bar3 = 100*bar3/sum(bar3);
[bar4 xs] = hist(bb_air_belt.SI,bins); bar4 = 100*bar4/sum(bar4);
% allBars = [100*bar1/sum(bar1);100*bar2/sum(bar2);100*bar3/sum(bar3);100*bar4/sum(bar4)];
allBars = [bar1;bar2;bar3;bar4];


ff = makeFigureWindow__one_axes_only(5,[6 4 5 2.5],[0.13 0.27 0.85 0.7]);
axes(ff.ha);hold on;
hb = bar(xs,allBars');
for ii = 1:length(hb)
    set(hb(ii),'facecolor',thisCols_all{ii});
end
xlim([bins(1) bins(end)]);
set(gca,'TickDir','out','FontSize',14,'FontWeight','Bold');
xlabel('Mutual Information (z-score)');
ylabel('Percentage');

x1 = 5; x2 = x1+3; y1 = (60:-4.5:0); y1 = y1(1:4); y2 = y1;
legendFontSize = 11;
legs = {'C1-Air','C1-Belt','C2-Air/Belt','C3-Air/Belt'};
for ii = 1:length(legs)
    plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols_all{ii},'linewidth',2);
    text(x2+0.5,y1(ii),sprintf('%s',legs{ii}),'Color',thisCols_all{ii},'FontSize',legendFontSize);
end

cs1 = cumsum(bar1);cs2 = cumsum(bar2);
cs3 = cumsum(bar3);cs4 = cumsum(bar4);
axesPos = ff.pos + [0.6 0.35 0 0];
axesPos(3:4) = [0.25 0.3];
axes('Position',axesPos);hold on;
plot(xs,cs1,'color',thisCols_all{1},'linewidth',1.5);
plot(xs,cs2,'color',thisCols_all{2},'linewidth',1.5);
plot(xs,cs3,'color',thisCols_all{3},'linewidth',1.5);
plot(xs,cs4,'color',thisCols_all{4},'linewidth',1.5);
xlim([bins(1) bins(end)]);
ylim([0 100]);
set(gca,'TickDir','out','FontSize',8,'FontWeight','Bold');
% xlabel('MI(z-score)');
% ylabel('%');

save2pdf('distSIs.pdf',ff.hf,600);

annovaVar = [nbc_air.SI' nbc_belt.SI' bc_air_belt.SI' bb_air_belt.SI'];
[p,tbl,stats] = anova1(annovaVar,legs);
% [p,tbl,stats] = kruskalwallis(y2wa,subst,'on');
figure(2001);
[c,~,~,gnames] = multcompare(stats,'CType','hsd');
pdf = c(:,6);
hdf = pdf<0.05;
selGroups = [1 2 3 4];
varNames = {'nbc_air.SI','nbc_belt.SI','bc_air_belt.SI','bb_air_belt.SI'};
for iii = 1:length(selGroups)
    ii = selGroups(iii);
    cmdText = sprintf('thisVar = %s(:);',varNames{ii});
    eval(cmdText);
    cmdText = sprintf('mVals(iii) = nanmean(thisVar);');
    eval(cmdText);
    cmdText = sprintf('semVals(iii) = std(thisVar)/sqrt(length(thisVar));');
    eval(cmdText);
end

combs = nchoosek(1:length(selGroups),2);

for ii = 1:size(combs,1)
    inds = ismember(c(:,[1 2]),[selGroups(combs(ii,1)) selGroups(combs(ii,2))],'rows');
    prT(ii,1) = c(inds,6);
end
hrT = prT<0.05;

% bar graph from anova
% ff = makeFigureWindow__one_axes_only(2,[1 4 1.75 2],[0.2 0.35 0.7 0.6]);
ff = makeFigureWindow__one_axes_only(3,[2 4 2 2],[0.2 0.17 0.75 0.81]);
axes(ff.ha);
plotBarsWithSigLines(mVals,semVals,combs,[hrT prT],'colors',thisCols_all(selGroups),'ySpacingFactor',10);
xlim([0.4 0.6+length(selGroups)]);
ylim([0 10]);
hyl = ylabel('Average MI (z-score)');
pos = get(hyl,'Position');pos = pos + [+0.1 0 0];set(hyl,'Position',pos);
set(ff.ha,'linewidth',1.25);
set(ff.ha,'TickDir','out','FontSize',9);
set(ff.ha,'XTickLabel',legs);
xtickangle(20);
save2pdf('bargraphSImeans.pdf',ff.hf,600);