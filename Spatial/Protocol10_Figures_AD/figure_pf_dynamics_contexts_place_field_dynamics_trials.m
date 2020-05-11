function figure_pf_dynamics_coontexts_place_field_dynamics_trials(fn,allRs,ccs)

adata = evalin('base','data');
mData = evalin('base','mData');
colors = mData.colors;
sigColor = mData.sigColor;
mData.axes_font_size = mData.axes_font_size-1;
selAnimals = 1:3;
selAnimals = 6:11;
mData.belt_length = adata{selAnimals(1)}{1}{1}.belt_length;
n = 0;

%%
selContexts = [2 3];
selCells = selectCells10(selAnimals,'Common23');
% selCells = selectCells10(selAnimals,'Disrupted_C1');
% selCells = selectCells10(selAnimals,'Remained_C1');

%%

trials = 3:10;
trials10 = 3:9;
trialsC = {[1:2],[3:4],[5:6],[7:10]};
trialsC10 = {[1:2],[3:4],[5:6],[7:9]};
for kk = 1:length(trialsC)
for ii = 1:4
    distDi = [];
    mRsi = [];
    for jj = 1:length(selAnimals)
        [ii jj selAnimals(jj)]
         [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'placeCells5',ii);
         distDi = [distDi tempD];
         [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'rasters',ii);
         try
             trials = trialsC{kk};
             mR = findMeanRasters(tempD,trials);
         catch
             trials = trialsC10{kk};
             mR = findMeanRasters(tempD,trials10);
         end
         mRsi = [mRsi;mR];
    end
    distD{kk,ii} = distDi;
    allRs{kk,ii} = mRsi;
    [allP{kk,ii},allC{kk,ii}] = findPopulationVectorPlot(mRsi,selCells);
end
end
pcs5 = distD; distD = [];
n=0;
%%
runThis = 0;
if runThis
context = 4;
ff = makeFigureRowsCols(105,[25 0.5 4 1],'RowsCols',[2 4],...
    'spaceRowsCols',[-0.06 -0.05],'rightUpShifts',[0.1 0.04],'widthHeightAdjustment',...
    [20 -15]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[15 6 3.5 2]);
FS = 7;
for sii = 1:4
    axes(ff.h_axes(1,sii));changePosition(gca,[0 0.05 -0.091 -0.1]);
    imagesc(allP{sii,context});
    P = allP{sii,context};
    axis off
    if sii == 1
        text(-5,1,'1','FontSize',FS,'FontWeight','Normal');
        if size(P,1) < 100
            text(-10,size(P,1),sprintf('%d',size(P,1)),'FontSize',FS,'FontWeight','Normal');
        else
            text(-13,size(P,1),sprintf('%d',size(P,1)),'FontSize',FS,'FontWeight','Normal');
        end
        text(-21,floor(size(P,1)/4),sprintf('Cells'),'FontSize',FS+2,'FontWeight','Bold','rotation',90);
    end
    trials = trialsC{sii};
    text(3,size(P,1)+round(size(P,1)/10),sprintf('Trials %d-%d',trials(1),trials(end)),'FontSize',FS,'FontWeight','Normal');
    set(gca,'Ydir','Normal','linewidth',1.5,'FontSize',FS,'FontWeight','Bold','YTick',[1 size(allP{sii,context},1)]);
    cols = 50;%size(Rs.rasters(:,:,1),2);
    colsHalf = ceil(cols/2);
    ts = round(adata{selAnimals(1)}{1}{1}.dist);
    set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    set(gca,'XTick',[]);
    axes(ff.h_axes(2,sii));
    dec = -0.09;
    changePosition(gca,[0.0 0.05 dec dec]);
    imagesc(allC{sii,context},[-1 1]);
    minC(sii) = min(allC{sii,context}(:));
    box off;
    axis equal
    axis off
    if sii == 1        
        text(-8,3,'0','FontSize',FS,'FontWeight','Normal');
        text(-13,50,sprintf('%d',mData.belt_length),'FontSize',FS,'FontWeight','Normal');
    end
    text(-1,-3,'0','FontSize',FS,'FontWeight','Normal');
    text(44,-3,sprintf('%d',mData.belt_length),'FontSize',FS,'FontWeight','Normal');
    if sii == 2
        text(35,-13,sprintf('Position (cm)'),'FontSize',FS+2,'FontWeight','Bold','rotation',0);
    end
    if sii == 1
        text(-21,-3,sprintf('Position (cm)'),'FontSize',FS+2,'FontWeight','Bold','rotation',90);
    end
    set(gca,'Ydir','Normal','linewidth',1,'FontSize',FS,'FontWeight','Bold');
    h = xlabel('Position (cm)');
    changePosition(h,[0 0 0]);
    h = ylabel('Position (cm)');
    changePosition(h,[1 0 0]);
%     cols = 50;%size(Rs.rasters(:,:,1),2);
    colsHalf = ceil(cols/2);
%     ts = round(Rs.dist);
    set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1) ts(colsHalf) ts(cols)],'Ydir','Normal');
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
save_pdf(ff.hf,mData.pdf_folder,sprintf('figure_pf_dynamics_trials_%d',context),600);
return;
end
n = 0;
%% peak shifts 1-2 
runThis = 1;
if runThis
peakpositions = [1 2];
for iii = 1:length(selContexts)
    ii = selContexts(iii);
    allPeakPos = [];
    for sii = 1:length(peakpositions)
        siii = peakpositions(sii);
        P = allP{siii,ii};
        [~,temp] = max(P,[],2);
        allPeakPos(:,sii) = adata{selAnimals(1)}{1}{1}.dist(temp);
    end
    dAllPeakPos{iii} = (diff(allPeakPos,[],2))';%allPeakPos(:,1)'-allPeakPos(:,2)';
end
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[15 6 3.5 2],'color','w');hold on;
distD = []; distD = dAllPeakPos;
[ha,hb,hca,sigR] = plotDistributions(distD,'colors',colors(selContexts),'maxY',80,'cumPos',[0.6 0.5 0.2 0.2],'min',-50,'incr',5,'max',40);
legs = {'C 1','C 2','C 3','C 4'}; legs = legs(selContexts); 
legs{length(legs)+1} = [-20 3 60 8];
putLegend(ha,legs,'colors',colors(selContexts),'sigR',{sigR,'anova',sigColor,8});
h = xlabel('Difference Peak Position (cm)');changePosition(h,[0 -7.0 0]);
h = ylabel('Percentage');changePosition(h,[-0.5 -1 0]);
set(gca,'FontSize',mData.axes_font_size,'FontWeight','Bold');changePosition(ha,[0 0.11 0.08 -0.09]);
axes(hca);set(gca,'FontSize',9);
fileName = sprintf('figure_pf_dynamics_contexts_peak_shifts_%d_%d.pdf',peakpositions(1),peakpositions(2));
save_pdf(hf,mData.pdf_folder,fileName,600);
return;
end

%% peak shifts 1-2 
runThis = 1;
if runThis
peakpositions = [3 4];
for ii = 1:4
    allPeakPos = [];
    for sii = 1:length(peakpositions)
        siii = peakpositions(sii);
        P = allP{siii,ii};
        [~,temp] = max(P,[],2);
        allPeakPos(:,sii) = adata{selAnimals(1)}{1}{1}.dist(temp);
    end
    dAllPeakPos{ii} = (diff(allPeakPos,[],2))';%allPeakPos(:,1)'-allPeakPos(:,2)';
end
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[15 6 3.5 2],'color','w');hold on;
distD = []; distD = dAllPeakPos;
[ha,hb,hca,sigR] = plotDistributions(distD,'colors',colors,'maxY',80,'cumPos',[0.6 0.5 0.2 0.2],'min',-50,'incr',5,'max',40);
legs = {'C 1','C 2','C 3','C 4',[-20 3 80 8]};
putLegend(ha,legs,'colors',colors,'sigR',{sigR,'anova',sigColor,8});
h = xlabel('Difference Peak Position (cm)');changePosition(h,[0 -7.0 0]);
h = ylabel('Percentage');changePosition(h,[-0.5 -1 0]);
set(gca,'FontSize',mData.axes_font_size,'FontWeight','Bold');changePosition(ha,[0 0.1 0.08 -0.09]);
axes(hca);set(gca,'FontSize',8);
fileName = sprintf('figure_pf_dynamics_peak_shifts_%d_%d.pdf',peakpositions(1),peakpositions(2));
save_pdf(hf,mData.pdf_folder,fileName,600);
return;
end