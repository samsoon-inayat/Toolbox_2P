function figure4_populationVectors_dynamics(fn,allRs,ccs)

aei = evalin('base','ei');
ei = aei([1 2 3 4 5 6]);
cei = combineRecordings(ei);
% cei.cellsN = size(ei{1}.signals,1);
% cei.cellsInds = 1:cei.cellsN;

mData.cellsN = cei.cellsN;
mData.cellsInds(1,:) = cei.cellsInds;
mData.cellsInds(2,:) = cei.cellsIndsRecs;
mData.sith = [10 10 10 10 5 5 5 5];
% mData.dataList = dataList;

selContexts = [1 2 3 4 5 6 7 8 9];
[a_ddata,a_tdata] = getDataContexts(cei,selContexts,'air');
[ai_ddata,ai_tdata] = getDataContexts(cei,selContexts,'airI');
[moffaon_ddata,moffaon_tdata] = getDataContexts(cei,selContexts,'motionOffsetAirOnset');
[b_ddata,b_tdata] = getDataContexts(cei,selContexts,'belt');
[moo_ddata,moo_tdata] = getDataContexts(cei,selContexts,'motionOnsetsOffsets');
[m_ddata,m_tdata] = getDataContexts(cei,selContexts,'motionI');
[mon_ddata,mon_tdata] = getDataContexts(cei,selContexts,'motionOnsets22');
[moff_ddata,moff_tdata] = getDataContexts(cei,selContexts,'motionOffsets22');
[aon_ddata,aon_tdata] = getDataContexts(cei,selContexts,'airOnsets22');
[aoff_ddata,aoff_tdata] = getDataContexts(cei,selContexts,'airOffsets22');
% cells which respond to belt distance in context 1 and 2 --> landmark
% vector cells
 % these are cells which could be considered pure path integrators responding to reference as air onset
% ids = [1 5]; notPurePIs = selectCells(data,mData,ids,[0 1],[sith sith]); % these are cells which are not pure path integrators necessarily
n = 0;

%%

ff = makeFigureRowsCols(105,[25 0.5 4 1],'RowsCols',[2 4],...
    'spaceRowsCols',[-0.06 -0.065],'rightUpShifts',[0.1 0.04],'widthHeightAdjustment',...
    [30 -15]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[25 4 3.5 2]);

data_1 = a_ddata(1:5);
data_2 = b_ddata(1:5);
data = {data_1{1:2} data_2{3:4}};
numberOfRows = 2;
sith = 7.5;
trials = 3:10;
cells1 = selectCells(data,mData,1,[1],[sith]);
cells2 = selectCells(data,mData,2,[1],[sith]);
eCells2 = setdiff(cells2,cells1);
cells3 = selectCells(data,mData,3,[1],[sith]);
eCells3 = setdiff(cells3,cells1);
cells4 = selectCells(data,mData,4,[1],[sith]);
eCells4 = setdiff(cells4,cells1);
purePIs = eCells2;

allPeakPos = []; allPR = []; allP = []; 
for sii = 1:4
    % trials = 1:4;
    allRs = data(sii);
    allRsA = a_ddata(1);
    ccs = purePIs;
    ptc = findMeanRasters(data{2});
    ptc = ptc(ccs,:);
    ptc = normalizeSignal(ptc,2);
    [~,peakPos] = max(ptc,[],2);
    [~,cellNums] = sort(peakPos);
    % [~,cellNums] = sort(allRsA{1}.SI(purePIs));
%     cellNums(1:10)
    FS = 5;

%     rows = numberOfRows;%length(allRs);

%     indices = 1:(2*rows);
%     indices = reshape(indices,2,rows);
%     indices = indices';
    % figure(fn);clf;
    for ii = 1%:rows
        Rs = allRs{ii};
        P = []; C = [];
        [P,C] = findPopulationVectorPlot(Rs,ccs,trials,cellNums);
        [~,temp] = max(P,[],2);
        allPeakPos(:,sii) = Rs.dist(temp);
        axes(ff.h_axes(1,sii));
        changePosition(gca,[0 0.05 -0.091 -0.1]);
        imagesc(P);
        temp = findMeanRasters(Rs,trials);
        temp = temp(ccs,:);
        allPR(:,:,sii) = temp(cellNums,:);
        allP(:,:,sii) = P;
    %     colorbar
        axis off
        
        if sii == 1
            text(-5,1,'1','FontSize',FS,'FontWeight','Normal');
            text(-13,size(P,1)-3,sprintf('%d',size(P,1)),'FontSize',FS,'FontWeight','Normal');
            text(-21,25,sprintf('Cells'),'FontSize',FS+3,'FontWeight','Bold','rotation',90);
        end

        text(3,size(P,1)+round(size(P,1)/10),sprintf('Context %d',sii),'FontSize',FS,'FontWeight','Normal');
            % caxis([minC maxC]);
        set(gca,'Ydir','Normal','linewidth',1.5,'FontSize',FS,'FontWeight','Bold','YTick',[1 length(purePIs)]);
%         text(1,-11,'0','FontSize',FS,'FontWeight','Normal');
%         text(44,-11,'142','FontSize',FS,'FontWeight','Normal');
%         if sii == 2
%             text(3,-30,sprintf('Position (cm)'),'FontSize',FS+3,'FontWeight','Bold','rotation',0);
%         end
        cols = size(Rs.rasters(:,:,1),2);
        colsHalf = ceil(cols/2);
        ts = round(Rs.dist);
        set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
        set(gca,'XTick',[]);
    %     h = xlabel('Position (cm)');
    %     changePosition(h,[0 0 0]);
        axes(ff.h_axes(2,sii));
        dec = -0.09;
        changePosition(gca,[0.0 0.05 dec dec]);
        allCs(:,:,sii) = C;
        hCCi(sii) = imagesc(C,[-1 1]);
        miI(sii) = min(C(:));
%         maI(sii) = max(C(:));
    %     xlim([0 142]);
    %     ylim([0 142]);
        box off;
    %     hc = colorbar
    %     changePosition(hc,[0.0 0.05 0 0.33]);
        colormap jet;
        % caxis([minC maxC]);
        axis equal
        axis off
        if sii == 1        
            text(-8,3,'0','FontSize',FS,'FontWeight','Normal');
            text(-15,50,'142','FontSize',FS,'FontWeight','Normal');
        end
        text(-1,-3,'0','FontSize',FS,'FontWeight','Normal');
        text(44,-3,'142','FontSize',FS,'FontWeight','Normal');
        if sii == 2
            text(35,-13,sprintf('Position (cm)'),'FontSize',FS+3,'FontWeight','Bold','rotation',0);
        end
        if sii == 1
            text(-21,-3,sprintf('Position (cm)'),'FontSize',FS+3,'FontWeight','Bold','rotation',90);
        end
        set(gca,'Ydir','Normal','linewidth',1,'FontSize',FS,'FontWeight','Bold');
        h = xlabel('Position (cm)');
        changePosition(h,[0 0 0]);
        h = ylabel('Position (cm)');
        changePosition(h,[1 0 0]);
        cols = size(Rs.rasters(:,:,1),2);
        colsHalf = ceil(cols/2);
        ts = round(Rs.dist);
        set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
        set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1) ts(colsHalf) ts(cols)],'Ydir','Normal');
    end
end
colormap parula
mI = min(miI);
for ii = 1:3
    axes(ff.h_axes(2,sii));
    caxis([mI 1]);
%     set(hCCi(ii),
end

axes(ff.h_axes(2,4));
pos = get(gca,'Position');
h = axes('Position',pos);
changePosition(h,[0.14 0.06 0.0 -0.1]);
hc = colorbar; axis off
set(hc,'YTick',[],'linewidth',0.01);
changePosition(hc,[0 0 0 -0.01]);
ylims = ylim;
xlims = xlim;
text(xlims(1)+0.33,ylims(1)-0.05,sprintf('%.2f',mI),'FontSize',4.5);
text(xlims(1)+0.39,ylims(2)+0.045,'1','FontSize',4.5);

axes(ff.h_axes(1,4));
pos = get(gca,'Position');
h = axes('Position',pos);
changePosition(h,[0.14 0.06 0.0 -0.1]);
hc = colorbar; axis off
set(hc,'YTick',[],'linewidth',0.01);
changePosition(hc,[0 0 0 -0.01]);
ylims = ylim;
xlims = xlim;
text(xlims(1)+0.4,ylims(1)-0.05,sprintf('0'),'FontSize',4.5);
text(xlims(1)+0.39,ylims(2)+0.045,'1','FontSize',4.5);

save2pdf('figure_2_pop_vec_dynamics.pdf',ff.hf,600);

%% Significance tests
% temp = ones(size(allCs(:,:,1)));
% temp = triu(temp,1);
% inds = find(temp);
% temp = allCs(:,:,1);
% annovaVar(:,1) = temp(inds);
% temp = allCs(:,:,2);
% annovaVar(:,2) = temp(inds);
% temp = allCs(:,:,3);
% annovaVar(:,3) = temp(inds);
% 
annovaVar = reshape(allP,size(allP,1),200);
groups = [ones(1,50) ones(1,50)*2 ones(1,50)*3 ones(1,50)*4];

[p,tbl,stats] = anova1(annovaVar,groups);%,names);
    % [p,tbl,stats] = kruskalwallis(y2wa,subst,'on');
figure(2001);
[c,~,~,gnames] = multcompare(stats,'CType','hsd');


%% widths centers

ff = makeFigureWindow__one_axes_only(106,[6 4 5 2.5],[0.19 0.2 0.76 0.75]);
set(gcf,'color','w');
set(gcf,'Position',[25 4 2 2]);
axes(ff.ha);hold on;

thisCols_all = {[0 0 0],'b','r',[0 0.7 0.3],'m','c'};

incr = 20;
distD(1,:) = allPeakPos(:,2)'-allPeakPos(:,1)';
distD(2,:) = allPeakPos(:,3)'-allPeakPos(:,2)';
maxB = ceil(max(distD(:)));
minB = ceil(min(distD(:)));
bins = minB:incr:maxB;
allBars = [];
for ii = 1:size(distD,1)
    bd = distD(ii,:);
%     bd(isnan(bd)) = [];
    [bar1 xs] = hist(bd,bins); bar1 = 100*bar1/sum(bar1);
    allBars = [allBars;bar1];
end

hb = bar(xs,allBars');
for ii = 1:length(hb)
    set(hb(ii),'facecolor',thisCols_all{ii},'barwidth',0.7,'EdgeColor','none');
end
xlim([bins(1)-(incr/2) bins(end)+(incr/0.5)]);
ylim([0 max(allBars(:))+15]);
set(gca,'TickDir','out','FontSize',8,'FontWeight','Bold');
h = xlabel('Difference Peak Position (cm)');
changePosition(h,[0 0.5 0]);
h = ylabel('Percentage');
changePosition(h,[0.2 0 0]);

% legends
x1 = -110; x2 = x1+1; y1 = (80:-7:0); y1 = y1(1:4); y2 = y1;
legendFontSize = 8;
legs = {'Trials 3-6 - Trials 1-2','Trials 7-11 - Trials 3-6'};
for ii = 1:length(legs)
%     plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols_all{ii},'linewidth',0.5);
    text(x2+0.5,y1(ii),sprintf('%s',legs{ii}),'Color',thisCols_all{ii},'FontSize',legendFontSize,'FontWeight','Bold');
end

axesPos = ff.pos + [0.5 0.25 0 0];
axesPos(3:4) = [0.3 0.25];
axes('Position',axesPos);hold on;
xlims = xlim;
for ii = 1:size(distD,1)
    cs = cumsum(allBars(ii,:));
    plot(xs,cs,'color',thisCols_all{ii},'linewidth',0.5);
end
% ytick(1:2) = [];
xlim([bins(1)-(incr/2) bins(end)+(incr/2)]);
ylim([0 100]);
set(gca,'TickDir','out','FontSize',6,'FontWeight','Bold');
% h = xlabel('Difference');
% changePosition(h,[0 0.5 0]);
% h = ylabel('%');
% changePosition(h,[40 0 0]);
title('Cumulative');

save2pdf('peak_shifts_dist_trials.pdf',ff.hf,600);


function [ptc,CRc] = findPopulationVectorPlot(Rs,ccs,trials,cellNums)
if ~isempty(trials)
    ptc = findMeanRasters(Rs,trials);
else
    ptc = findMeanRasters(Rs);
end
ptc = ptc(ccs,:);
ptc = normalizeSignal(ptc,2);
if ~exist('cellNums','var')
    [~,peakPos] = max(ptc,[],2);
    [~,cellNums] = sort(peakPos);
    ptc = ptc(cellNums,:);
else
    ptc = ptc(cellNums,:);
end
CRc = corrcoef(ptc);





