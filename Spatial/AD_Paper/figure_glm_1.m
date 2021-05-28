function figure_glm_1(fn,allRs,ccs)

ei = evalin('base','ei10');
mData = evalin('base','mData');
glms = evalin('base','glms');
glmsI = evalin('base','glmsI');
colors = mData.colors;
sigColor = mData.sigColor;
selAnimals = 1:8;
% selAnimals = 5:8;
% selAnimals = 3;
mData.belt_length = ei{selAnimals(1)}.b.belt_length;
n = 0;

%%
selCells = 'All';
planeNumbers = 1;
maxDistTime = [140 5];
contextNumbers = [1 1];
% stimMarkers = {'air','air','air','air'};
% rasterTypes = {'dist','time','dist','dist'};
trials = 3:10;
trials10 = 3:9;
cns = [];
gAllVals = [];
for jj = 1:length(selAnimals)
    for ii = 1:length(contextNumbers)
        contextNumber = contextNumbers(ii);
        if ii == 1
            [tdevs,tcns] = get_param_glm('deviances',ei,glms,selAnimals(jj),contextNumber);
        else
            [tdevs,tcns] = get_param_glm('deviances',ei,glmsI,selAnimals(jj),contextNumber);
        end
        areCells = logical(tcns(:,4));
        pTemp = tdevs.pvals>0;
        pSel = logical(ones(size(pTemp,1),1));
        for pp = 1:size(pTemp,2)
            pSel = pSel & pTemp(:,pp);
        end
        selCells = pSel & areCells;
        dataT{jj,ii} = tdevs.STD_T(selCells);
        dataD{jj,ii} = tdevs.STD_D(selCells);
        dataC{jj,ii} = tcns;
    end
end
thr = 100;
for rr = 1:size(dataT,1)
    inds = logical(zeros(size(dataT{rr,1})));
    for cc = 1:size(dataT,2)
        d_deviances = dataT{rr,cc} - dataD{rr,cc};
        data{rr,cc} = d_deviances;
    end
end
n = 0;
%%
runthis = 0;
if runthis
    minBin = -100;%min(gAllVals);
    maxBin = 100;%max(gAllVals);
    incr = (maxBin-minBin)/20;
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 6.9 2.5],'color','w');
    hold on;
%     [ha,hb,hca,sigR] = plotDistributions(data,'colors',colors,'maxY',70,'cumPos',[0.5 0.55 0.25 0.25],'min',minBin,'incr',incr,'max',maxBin);
    [ha,hb,hca,sigR] = plotDistributions(data,'colors',colors,'maxY',90,'min',minBin,'incr',incr,'max',maxBin);
    hold on;
    legs = [];
    for ii = 1:length(contextNumbers)
        legs{ii} = sprintf('Condition - %d (N = %d)',contextNumbers(ii),size(data,1));
    end
%     ylim([0 120]);
%     xlims = xlim; dx = xlims(2) - xlims(1); ylims = ylim; dy = ylims(2) - ylims(1);
%     legs{ii+1} = [xlims(1)+dx/1.5 dx/50 ylims(1)+dy/2 dy/15];
%     if size(data,2) > 2
%         putLegend(ha,legs,'colors',colors,'sigR',{sigR,'anova',sigColor,8});
%     else
%         putLegend(ha,legs,'colors',colors,'sigR',{sigR,'ttest',sigColor,10});
%     end
%     axes(ha);axes(hca);
%     h = xlabel('Mutual Information (bits)');%changePosition(h,[0 -dy/3 0]);
%     h = ylabel('Percentage');changePosition(h,[0 0 0]);
%     set(gca,'FontSize',axes_font_size,'FontWeight','Bold');changePosition(ha,[-0.04 0.09 0.13 -0.05]);
%     save_pdf(hf,mData.pdf_folder,'Distribution Of zMI',600);
    sigR.anova
    return;
end
%%

cns = dataC{1,1};
areCells = find(logical(dataC{1,1}(:,4)));
STD_T = dataT{1,1}; STD_D = dataD{1,1};

%%
d_Dev = STD_T - STD_D;
figure(40000);clf;
scatter(STD_T,STD_D);
set(gca,'XScale','log','YScale','log');
return;
d_dev_thr = 50;
LD = find(d_Dev < -d_dev_thr);
LD_D = LD;
for cc = 1:length(LD)
%     disp([cc length(L)]);
    nn = areCells(LD(cc));
    ann = cns(nn,1);
    pp = cns(nn,2);
    cn = cns(nn,3);
    [dataT,~,~] = getParamValues('',ei(ann),pp,contextNumber,'air','time','All',maxDistTime);
    [dataD,~,~] = getParamValues('',ei(ann),pp,contextNumber,'air','dist','All',maxDistTime);
%     figure(1000);clf;
%     subplot 121;imagesc(dataT.sp_rasters(:,:,cn));title(dataT.info_metrics.ShannonMI_Zsh(cn));
%     subplot 122;imagesc(dataD.sp_rasters(:,:,cn));title(dataD.info_metrics.ShannonMI_Zsh(cn));
%     pause;
    diff_T_D_D(cc) = (dataT.info_metrics.ShannonMI_Zsh(cn)) - (dataD.info_metrics.ShannonMI_Zsh(cn));
end
LD = find(d_Dev > d_dev_thr);
LD_T = LD;
for cc = 1:length(LD)
%     disp([cc length(L)]);
    nn = areCells(LD(cc));
    ann = cns(nn,1);
    pp = cns(nn,2);
    cn = cns(nn,3);
    [dataT,~,~] = getParamValues('',ei(ann),pp,contextNumber,'air','time','All',maxDistTime);
    [dataD,~,~] = getParamValues('',ei(ann),pp,contextNumber,'air','dist','All',maxDistTime);
%     figure(1000);clf;
%     subplot 121;imagesc(dataT.sp_rasters(:,:,cn));title(dataT.info_metrics.ShannonMI_Zsh(cn));
%     subplot 122;imagesc(dataD.sp_rasters(:,:,cn));title(dataD.info_metrics.ShannonMI_Zsh(cn));
%     pause;
    diff_T_D_T(cc) = (dataT.info_metrics.ShannonMI_Zsh(cn)) - (dataD.info_metrics.ShannonMI_Zsh(cn));
end
bins = -25:0.5:25;
barsD = hist(diff_T_D_D,bins);
barsT = hist(diff_T_D_T,bins);
bars = [barsD' barsT'];
figure(30000);clf;bar(bins,bars);
[h,p,ci,t_stat] = ttest2(diff_T_D_D,diff_T_D_T)
disp([nanmean(diff_T_D_D) nanmean(diff_T_D_T)]);
% return;
L = LD_T;

for cc = 1:length(L)
%     disp([cc length(L)]);
    nn = areCells(L(cc));
    ann = cns(nn,1);
    pp = cns(nn,2);
    cn = cns(nn,3);
    [dataT,~,~] = getParamValues('',ei(ann),pp,contextNumber,'air','time','All',maxDistTime);
    [dataD,~,~] = getParamValues('',ei(ann),pp,contextNumber,'air','dist','All',maxDistTime);
    firstHalfTrials = nanmean(dataT.sp_rasters_nan_corrected(4:7,:,cn));
    secondHalfTrials = nanmean(dataT.sp_rasters_nan_corrected(8:10,:,cn));
    timeCCt = corrcoef(firstHalfTrials,secondHalfTrials); timeCC(cc) = timeCCt(1,2);
    firstHalfTrials = nanmean(dataD.sp_rasters_nan_corrected(4:7,:,cn));
    secondHalfTrials = nanmean(dataD.sp_rasters_nan_corrected(8:10,:,cn));
    distCCt = corrcoef(firstHalfTrials,secondHalfTrials); distCC(cc) = distCCt(1,2);
end
figure(1100);clf;scatter(timeCC,distCC);
return;
for cc = 1:length(L)
    disp([cc length(L)]);
    nn = areCells(L(cc));
    ann = cns(nn,1);
    pp = cns(nn,2);
    cn = cns(nn,3);
    [dataT,~,~] = getParamValues('',ei(ann),pp,contextNumber,'air','time','All',maxDistTime);
    [dataD,~,~] = getParamValues('',ei(ann),pp,contextNumber,'air','dist','All',maxDistTime);
    figure(2000);clf;
    subplot 121;
    imagesc(dataT.sp_rasters_nan_corrected(:,:,cn));title(sprintf('Time - %.3f',dataT.info_metrics.ShannonMI_Zsh(cn)));
    subplot 122;
    imagesc(dataD.sp_rasters_nan_corrected(:,:,cn));title(sprintf('Dist -  %.3f',dataD.info_metrics.ShannonMI_Zsh(cn)));
    
    figure(2001);clf;
    subplot 121;
    firstHalfTrials = nanmean(dataT.sp_rasters_nan_corrected(4:7,:,cn));
    secondHalfTrials = nanmean(dataT.sp_rasters_nan_corrected(8:10,:,cn));
    plot(firstHalfTrials,'b');hold on; plot(secondHalfTrials,'r');title(sprintf('Time - %.3f',dataT.info_metrics.ShannonMI_Zsh(cn)));
    subplot 122;
    firstHalfTrials = nanmean(dataD.sp_rasters_nan_corrected(4:7,:,cn));
    secondHalfTrials = nanmean(dataD.sp_rasters_nan_corrected(8:10,:,cn));
    plot(firstHalfTrials,'b');hold on; plot(secondHalfTrials,'r');
    title(sprintf('Dist -  %.3f',dataD.info_metrics.ShannonMI_Zsh(cn)));
    pause;
    
end
return;

%%
ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 4],...
    'spaceRowsCols',[-0.01 -0.05],'rightUpShifts',[0.1 0.06],'widthHeightAdjustment',...
    [15 -30]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[1 6 3.5 2]);
FS = 4;
for sii = 1:length(stimMarkers)
    P = allP{sii};
    axes(ff.h_axes(1,sii));changePosition(gca,[0 0.05 -0.091 -0.1]);
    imagesc(P);
    box off;
%     axis off
%     if sii == 1
%         text(-5,1,'1','FontSize',FS,'FontWeight','Normal');
%         if size(P,1) < 100
%             text(-7,size(P,1),sprintf('%d',size(P,1)),'FontSize',FS,'FontWeight','Normal');
%         else
%             text(-10,size(P,1),sprintf('%d',size(P,1)),'FontSize',FS,'FontWeight','Normal');
%         end
        if sii == 1
            text(-21,25,sprintf('Cells'),'FontSize',FS+3,'FontWeight','Bold','rotation',90);
        end
%     end
    text(3,size(P,1)+round(size(P,1)/10),sprintf('Condition %d',sii),'FontSize',FS,'FontWeight','Normal');
    set(gca,'Ydir','Normal','linewidth',0.25,'FontSize',FS,'FontWeight','Bold','YTick',[1 size(P,1)]);
    cols = size(P,2);
    colsHalf = ceil(cols/2);
    ts = round(time_xs{sii}(1:cols));
%     set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    set(gca,'XTick',[]);
    axes(ff.h_axes(2,sii));
    dec = -0.09;
    changePosition(gca,[0.0 0.05 dec dec]);
    imagesc(allC{sii},[-1 1]);
    minC(sii) = min(allC{sii}(:));
    box off;
%     axis equal
%     axis off
%     if sii == 1        
%         text(-8,3,'0','FontSize',FS,'FontWeight','Normal');
%         text(-10,50,num2str(mData.belt_length),'FontSize',FS,'FontWeight','Normal');
%     end
%     text(-1,-3,'0','FontSize',FS,'FontWeight','Normal');
%     text(44,-3,num2str(mData.belt_length),'FontSize',FS,'FontWeight','Normal');
%     if sii == 2
%         text(35,-13,sprintf('Position (cm)'),'FontSize',FS+3,'FontWeight','Bold','rotation',0);
%     end
%     if sii == 1
%         text(-21,-3,sprintf('Position (cm)'),'FontSize',FS+3,'FontWeight','Bold','rotation',90);
%     end
    set(gca,'Ydir','Normal','linewidth',1,'FontSize',FS,'FontWeight','Bold');
    if strcmp(rasterTypes{sii},'dist')
        h = xlabel('Position (cm)');    changePosition(h,[0 0 0]);
    else
        h = xlabel('Time (sec)');    changePosition(h,[0 0 0]);
    end
    if sii == 1
%         h = ylabel('Position (cm)');    changePosition(h,[1 0 0]);
    end
    cols = size(P,2);
    colsHalf = ceil(cols/2);
    ts = round(time_xs{sii}(1:cols));
    set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
end

colormap parula
mI = min(minC);
for ii = 1:4
    axes(ff.h_axes(2,ii));
    caxis([mI 1]);
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

% delete(ff.h_axes(1,4));
% delete(ff.h_axes(2,4));

save_pdf(ff.hf,mData.pdf_folder,'figure_place_cells_py_10.pdf',600);

return;
%%
% bar graph from anova
sigR = significanceTesting(allP);
ff = makeFigureWindow__one_axes_only(3,[2 4 2 2],[0.3 0.22 0.68 0.72]);
set(gcf,'color','w');
set(gcf,'Position',[3 4 1.2 1.5]);
axes(ff.ha); hs = sigR.anova.multcompare.h; ps = sigR.anova.multcompare.p;
plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[hs ps],'colors',colors,'sigColor',sigColor,'maxY',0.3,'ySpacing',0.03,'sigTestName','ANOVA');
% plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[hs ps],'colors',colors,'sigColor',sigColor,'maxY',0.57,'ySpacing',0.04,'sigTestName','ANOVA');
xlim([0.4 0.6+length(sigR.means)]);
%     ylim([0 max(mVals)]);
hyl = ylabel('Average Normalized FR');
pos = get(hyl,'Position');pos = pos + [+0.4 0 0];set(hyl,'Position',pos);
% set(ff.ha,'linewidth',1);
set(ff.ha,'TickDir','out','FontSize',7,'FontWeight','Normal');
set(ff.ha,'XTick',[1 2 3 4],'XTickLabel',{'Context 1','Context 2','Context 3','Context4'});
xtickangle(25);
save2pdf('figure_place_cells_bars.pdf',ff.hf,600);