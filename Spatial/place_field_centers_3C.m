function place_field_centers_3C(tei)
% close all
thisCols_all = {[0 0 0],'b','r',[0 0.7 0.3],[0 0.7 0.3],'r',[0 0.7 0.3],[0 0 0],'b','m','r'};

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

% selCells{1} = 1:size(nbb_air.rasters,3);
zScoreTh = 5;
selCells{1} = find(nbc_belt.SI > zScoreTh);% & nbb_belt.SI < zScoreTh & bc_air_belt.SI > zScoreTh);
varNames = {'nbc_air','nbc_belt','bc_air_belt','bb_air_belt'};
for ii = 1:length(varNames)
    cmdTxt = sprintf('[pcw,pcc] = getPlaceCellProps(%s,selCells{1}); %s.pcw = pcw; %s.pcc = pcc;',varNames{ii},varNames{ii},varNames{ii});
    eval(cmdTxt);
end
% 
% varNames = {'bc_air_belt','bb_air_belt'};
% for ii = 1:length(varNames)
%     cmdTxt = sprintf('pfs = getPlaceCellPropsGauss(%s,selCells); %s.pfs = pfs;',varNames{ii},varNames{ii});
%     eval(cmdTxt);
% end

% 
% plotRastersMulti({nbc_air,nbc_belt,bc_air_belt,bb_air_belt},selCells{1},0,0);
% return;

% for distribution of place field centers
maxSI = 142;
minSI = 1;
bins = minSI:10:maxSI;
[bar1 xs] = hist(nbc_air.pcc(~isnan(nbc_air.pcc)),bins); bar1 = 100*bar1/sum(bar1);
[bar2 xs] = hist(nbc_belt.pcc(~isnan(nbc_belt.pcc)),bins); bar2 = 100*bar2/sum(bar2);
[bar3 xs] = hist(bc_air_belt.pcc(~isnan(bc_air_belt.pcc)),bins); bar3 = 100*bar3/sum(bar3);
[bar4 xs] = hist(bb_air_belt.pcc(~isnan(bb_air_belt.pcc)),bins); bar4 = 100*bar4/sum(bar4);
% allBars = [100*bar1/sum(bar1);100*bar2/sum(bar2);100*bar3/sum(bar3);100*bar4/sum(bar4)];
allBars = [bar1;bar2;bar3;bar4];

ff = makeFigureWindow__one_axes_only(6,[6 4 7 2.5],[0.13 0.27 0.85 0.7]);
axes(ff.ha);hold on;
hb = bar(xs,allBars');
for ii = 1:length(hb)
    set(hb(ii),'facecolor',thisCols_all{ii});
end
xlim([bins(1) bins(end)]);
ylim([0 25]);
set(gca,'TickDir','out','FontSize',14,'FontWeight','Bold');
xlabel('Position (cm)');
ylabel('Percentage');

x1 = 45; x2 = x1+3; y1 = (25:-1.5:0); y1 = y1(1:4); y2 = y1;
legendFontSize = 11;
legs = {'Air-NoBrake','Belt-NoBrake','Belt-Brake','Blank-Brake'};
for ii = 1:length(legs)
    plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols_all{ii},'linewidth',2);
    text(x2+0.5,y1(ii),sprintf('%s',legs{ii}),'Color',thisCols_all{ii},'FontSize',legendFontSize);
end

% cs1 = cumsum(bar1);cs2 = cumsum(bar2);
% cs3 = cumsum(bar3);cs4 = cumsum(bar4);
% axesPos = ff.pos + [0.6 0.35 0 0];
% axesPos(3:4) = [0.25 0.3];
% axes('Position',axesPos);hold on;
% plot(xs,cs1,'color',thisCols_all{1},'linewidth',1.5);
% plot(xs,cs2,'color',thisCols_all{2},'linewidth',1.5);
% plot(xs,cs3,'color',thisCols_all{3},'linewidth',1.5);
% plot(xs,cs4,'color',thisCols_all{4},'linewidth',1.5);
% xlim([bins(1) bins(end)]);
% ylim([0 100]);
% set(gca,'TickDir','out','FontSize',8,'FontWeight','Bold');


% figure(300);clf;
% SIs = [nbb_air.SI;nbb_belt.SI;bc_air_belt.SI;bb_air_belt.SI];
% imagesc(SIs);
% colorbar