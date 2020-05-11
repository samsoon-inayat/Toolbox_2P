function figure_place_field_dynamics_contexts(fn,allRs,ccs)

adata = evalin('base','data');
mData = evalin('base','mData');
tdata = evalin('base','dataT');
colors = mData.colors;
sigColor = mData.sigColor;
selAnimals = 1:11;
mData.belt_length = adata{selAnimals(2)}{1}{1}.belt_length;
n = 0;

%%
for ii = 1:2%length(data)
    excRi = []; inhRi = [];
    for jj = 1:length(selAnimals)
        [tempD cnsjj] = getVariableValues(tdata{selAnimals(jj)},'excR',ii);
        excRi = [excRi;tempD'];
        [tempD cnsjj] = getVariableValues(tdata{selAnimals(jj)},'inhR',ii);
        inhRi = [inhRi;tempD'];
    end
    excR{ii} = excRi; inhR{ii} = inhRi;
end
n = 0;
%%
selCells = selectCells15(selAnimals,'Remained_C2');
% selCells = selectCells15(selAnimals,'Only_C1');
prepost = 2;
% selCells = logical(inhR{prepost}(:,1));
% selCells = logical(SIs{2});
% selCells = selectCells15(selAnimals,'New_C2');
% selCells = selectCells15(selAnimals,'Disrupted_C1');
% selCells = selectCells15(selAnimals,'Remained_C1');
trials = 3:10;
trials10 = 3:9;
for ii = 1:3%length(data)
    mRsi = [];
    for jj = 1:length(selAnimals)
        [ii jj selAnimals(jj)]
        [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'rasters',ii);
        try
         mR = findMeanRasters(tempD,trials);
        catch
         mR = findMeanRasters(tempD,trials10);
        end
        mRsi = [mRsi;mR];
    end
    allRs{ii} = mRsi;
    [allP{ii},allC{ii}] = findPopulationVectorPlot(mRsi,selCells);
end
trials = 1:10;
for ii = 1:2%length(data)
    mRsi = [];
    SIsi = [];
    for jj = 1:length(selAnimals)
        [ii jj selAnimals(jj)]
        [tempD cnsjj] = getVariableValues(tdata{selAnimals(jj)},'rasters',ii);
         mR = findMeanRasters(tempD,trials);
        mRsi = [mRsi;mR];
        [tempD cnsjj] = getVariableValues(tdata{selAnimals(jj)},'SI',ii);
        SIsi = [SIsi;tempD'];
    end
    allRsT{ii} = mRsi;
    [allPT{ii},allCT{ii}] = findPopulationVectorPlot(mRsi,selCells);
    SIs{ii} = SIsi(selCells);
end
% allP{4} = allPT{prepost};
% allC{4} = allCT{prepost};
n = 0;
%
runThis = 1;
if runThis
ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 length(allP)],...
    'spaceRowsCols',[-0.06 -0.05],'rightUpShifts',[0.1 0.04],'widthHeightAdjustment',...
    [15 -15]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[1 6 3.5 2]);
FS = 10;
for sii = 1:length(allP)
    P = allP{sii};
    axes(ff.h_axes(1,sii));changePosition(gca,[0 0.05 -0.091 -0.1]);
    imagesc(P);
    axis off
%     if sii == 1
        text(-5,1,'1','FontSize',FS,'FontWeight','Normal');
        if size(P,1) < 100
            text(-7,size(P,1),sprintf('%d',size(P,1)),'FontSize',FS,'FontWeight','Normal');
        else
            text(-10,size(P,1),sprintf('%d',size(P,1)),'FontSize',FS,'FontWeight','Normal');
        end
        if sii == 1
            text(-21,25,sprintf('Cells'),'FontSize',FS+3,'FontWeight','Bold','rotation',90);
        end
%     end
    text(3,size(P,1)+round(size(P,1)/10),sprintf('Context %d',sii),'FontSize',FS,'FontWeight','Normal');
    set(gca,'Ydir','Normal','linewidth',1.5,'FontSize',FS,'FontWeight','Bold','YTick',[1 sum(selCells)]);
    cols = 50;
    colsHalf = ceil(cols/2);
    ts = round(adata{selAnimals(1)}{1}{1}.dist);
    set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    set(gca,'XTick',[]);
    axes(ff.h_axes(2,sii));
    dec = -0.09;
    changePosition(gca,[0.0 0.05 dec dec]);
    imagesc(allC{sii},[-1 1]);
    minC(sii) = min(allC{sii}(:));
    box off;
    axis equal
    axis off
    if sii == 1        
        text(-8,3,'0','FontSize',FS,'FontWeight','Normal');
        text(-10,50,num2str(mData.belt_length),'FontSize',FS,'FontWeight','Normal');
    end
    text(-1,-3,'0','FontSize',FS,'FontWeight','Normal');
    text(44,-3,num2str(mData.belt_length),'FontSize',FS,'FontWeight','Normal');
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
    cols = 50;%size(Rs.rasters(:,:,1),2);
    colsHalf = ceil(cols/2);
%     ts = round(Rs.dist);
    set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1) ts(colsHalf) ts(cols)],'Ydir','Normal');
end

colormap parula
mI = min(minC);
for ii = 1:3
    axes(ff.h_axes(2,ii));
    caxis([mI 1]);
end

axes(ff.h_axes(2,3));
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

axes(ff.h_axes(1,3));
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
% delete(ff.h_axes(1,4));delete(ff.h_axes(2,4));
save_pdf(ff.hf,mData.pdf_folder,'figure_pf_dynamics_contexts_population_vector_15.pdf',600);
% return;
end
n = 0;


