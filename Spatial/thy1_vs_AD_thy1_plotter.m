function thy1_vs_AD_thy1_plotter

% the following files are used to plot graphs on the poster
% thy1_vs_AD_thy1.m to load data
% then
% this file
% tracePlots.m
% behaviorPlot.m
% place_cell_properties.m




ei = evalin('base','ei');


%%
cMIs = [];
for ii = 1:4
    cMIs = [cMIs ei{ii}.rasters.zMI];
end

dMIs = [];
for ii = 5:8
    dMIs = [dMIs ei{ii}.rasters.zMI];
end

vals1 = cMIs;
vals2 = dMIs;

minVal = min([vals1 vals2]);
maxVal = max([vals1 vals2]);
% maxVal = 2;
bins = linspace(minVal,maxVal,20);
% bins = -30:5:70;
[bar1 xs] = hist(vals1,bins);
[bar2 xs] = hist(vals2,bins);
allBars = [100*bar1/sum(bar1);100*bar2/sum(bar2)];

figure(100);clf;
hbars = bar(xs,allBars');
set(hbars(1),'facecolor','b');
set(hbars(2),'facecolor','r');
xlim([bins(1) bins(end)]);
set(gca,'TickDir','out','FontSize',14,'FontWeight','Bold','linewidth',1.5);box off;
xlabel('Normalized Mutual Information (Z-Score)');
ylabel('Percentage');
hold on;
legendText = {'APP(NL-G-F)-Thy1-GCaMP','Thy1-GCaMP'};
legsN = [2 1];
thisCols = {'r','b'};
x1 = 6.5; x2 = x1+1; y1 = (27:1.5:100); y1 = y1(1:2); y2 = y1;
legendFontSize = 11;
legs = legendText;
for ii = 1:length(legs)
    plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols{ii},'linewidth',2);
    eval(sprintf('len = length(vals%d);',legsN(ii)));
    text(x2+0.5,y1(ii),sprintf('%s, (n=%d cells, 2 animals)',legs{ii},len),'Color',thisCols{ii},'FontSize',legendFontSize);
end

axesPos = get(gca,'Position');
cAllBars = cumsum(allBars,2);
subAxesPos(1) = axesPos(1) + 0.45;
subAxesPos(2) = axesPos(2) + 0.3;
subAxesPos(3) = 0.25;
subAxesPos(4) = 0.25;
axes('Position',subAxesPos);
plot(xs,cAllBars(1,:),'linewidth',1.5,'color','b');hold on;
plot(xs,cAllBars(2,:),'linewidth',1.5,'color','r');
set(gca,'TickDir','out','FontSize',12,'FontWeight','Bold','linewidth',1.5);box off
ht = title('Cumulative Distribution');
set(ht,'FontSize',12);
tPos = get(ht,'Position');
set(ht,'Position',(tPos + [0 3 0]));
xlim([0 max(xs)]);

% save2pdf('SI_diff_M2_RSC.pdf',gcf,600);
[p,h] = ranksum(vals1,vals2,'Tail','right')
[ht,pt] = ttest2(vals1,vals2,'Tail','right')
[cM cE cS] = findMeanAndStandardError(cMIs);
[dM dE dS] = findMeanAndStandardError(dMIs);

% text(0,-40,sprintf('***, p = %.3f, Student''s t-test',p),'FontSize',12,'FontWeight','bold');
save2pdf('Dist_Norm_MI_Thy1VsThy1AD',gcf,600);

%%
% CC_populationVector_c.pdf
recs = {[1 2],[3 4],[5 6],[7 8]};
for ii = 1:4
    normSel = 1;
    rec = recs{ii};
    ptc = getPositionTuning(ei(rec),'sp');
    [~,peakPos] = max(ptc,[],2);[~,cellNums] = sort(peakPos);
    ptc = ptc(cellNums,:);
    ptc = normalizeSignal(ptc,2);
    CRc = corrcoef(ptc);
    hf = figure(1000);clf;tempPos = get(gcf,'Position');tempPos(3:4) = [8.75 2.5];set(gcf,'Units','Inches');set(gcf,'Position',tempPos);
    subplot 121
    imagesc(ptc);colorbar
    colormap jet;
    % caxis([minC maxC]);
    set(gca,'Ydir','Normal','linewidth',1.5,'FontSize',12,'FontWeight','Bold');
    xlabel('Position (cm)');ylabel('Cells');
    set(gca,'XTick',[1 25 50],'XTickLabel',[0 71 142]);
    text(73,100,'Normalized Spike Rate','Rotation',90,'FontWeight','Bold');

    subplot 122
    imagesc(CRc);colorbar
    colormap jet;
    % caxis([minC maxC]);
    set(gca,'Ydir','Normal','linewidth',1.5,'FontSize',12,'FontWeight','Bold');
    axis equal
    set(gca,'XTick',[1 25 50],'YTick',[1 25 50],'XTickLabel',[0 71 142],'YTickLabel',[0 71 142]);
    xlabel('Position (cm)');ylabel('Position (cm)');
    text(76,-5,'Pearson Correlation Coefficient','Rotation',90,'FontWeight','Bold');
    save2pdf(sprintf('CC_populationVector_c%d.pdf',ii),hf,600);
    CRAll(:,:,ii) = CRc;
end

CRc = mean(CRAll(:,:,[1 2]),3);
CRd = mean(CRAll(:,:,[3 4]),3);

minC = min([CRc(:)' CRd(:)']);
maxC = max([CRc(:)' CRd(:)']);

hf = figure(1000);clf;tempPos = get(gcf,'Position');tempPos(3:4) = [5 4];set(gcf,'Units','Inches');set(gcf,'Position',tempPos);
imagesc(CRc);colorbar
colormap jet;
caxis([minC maxC]);
set(gca,'Ydir','Normal','linewidth',1.5,'FontSize',12,'FontWeight','Bold');
axis equal
set(gca,'XTick',[1 25 50],'YTick',[1 25 50],'XTickLabel',[0 71 142],'YTickLabel',[0 71 142]);
xlabel('Position (cm)');ylabel('Position (cm)');
save2pdf('CC_populationVector_c.pdf',hf,600);

hf = figure(1000);clf;tempPos = get(gcf,'Position');tempPos(3:4) = [5 4];set(gcf,'Units','Inches');set(gcf,'Position',tempPos);
imagesc(CRd);colorbar
colormap jet;
caxis([minC maxC]);
set(gca,'Ydir','Normal','linewidth',1.5,'FontSize',12,'FontWeight','Bold');
axis equal
set(gca,'XTick',[1 25 50],'YTick',[1 25 50],'XTickLabel',[0 71 142],'YTickLabel',[0 71 142]);
xlabel('Position (cm)');ylabel('Position (cm)');
save2pdf('CC_populationVector_d.pdf',hf,600);
n = 0;

% tempO = ones(size(CRc));
% ltr = tril(tempO,-1);
% inds = find(ltr);
% vals1 = CRc(inds)'; vals2 = CRd(inds)';
% 
% minVal = min([vals1 vals2]);
% maxVal = max([vals1 vals2]);
% % maxVal = 2;
% bins = linspace(minVal,maxVal,20);
% % bins = -30:5:70;
% [bar1 xs] = hist(vals1,bins);
% [bar2 xs] = hist(vals2,bins);
% allBars = [100*bar1/sum(bar1);100*bar2/sum(bar2)];
% 
% figure(101);clf;
% hbars = bar(xs,allBars');
% set(hbars(1),'facecolor','b');
% set(hbars(2),'facecolor','r');
% xlim([bins(1) bins(end)]);
% set(gca,'TickDir','out','FontSize',14,'FontWeight','Bold','linewidth',1.5);box off;
% xlabel('Pearson Correlation Coefficients');
% ylabel('Percentage');
% hold on;
% legendText = {'APP(NL-G-F)-Thy1-GCaMP','Thy1-GCaMP'};
% legsN = [2 1];
% thisCols = {'r','b'};
% x1 = 0.1; x2 = x1+0.05; y1 = (23:1:100); y1 = y1(1:2); y2 = y1;
% legendFontSize = 11;
% legs = legendText;
% for ii = 1:length(legs)
%     plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols{ii},'linewidth',2);
%     eval(sprintf('len = length(vals%d);',legsN(ii)));
%     text(x2+0.01,y1(ii),sprintf('%s, (n=%d)',legs{ii},len),'Color',thisCols{ii},'FontSize',legendFontSize);
% end
% 
% axesPos = get(gca,'Position');
% cAllBars = cumsum(allBars,2);
% subAxesPos(1) = axesPos(1) + 0.45;
% subAxesPos(2) = axesPos(2) + 0.3;
% subAxesPos(3) = 0.25;
% subAxesPos(4) = 0.25;
% axes('Position',subAxesPos);
% plot(cAllBars(1,:),'linewidth',1.5,'color','b');hold on;
% plot(cAllBars(2,:),'linewidth',1.5,'color','r');
% set(gca,'TickDir','out','FontSize',12,'FontWeight','Bold','linewidth',1.5);box off
% ht = title('Cumulative Distribution');
% set(ht,'FontSize',12);
% tPos = get(ht,'Position');
% set(ht,'Position',(tPos + [0 3 0]));

% 
% % save2pdf('SI_diff_M2_RSC.pdf',gcf,600);
% [p,h] = ranksum(vals1,vals2,'Tail','left')
% [ht,pt] = ttest2(vals1,vals2,'Tail','left')
% 
% text(0,-40,sprintf('***, p = %.3f, Student''s t-test',p),'FontSize',12,'FontWeight','bold');
% save2pdf('PCC_Thy1VsThy1AD.pdf',gcf,600);
% n=0;

place_field_widths;