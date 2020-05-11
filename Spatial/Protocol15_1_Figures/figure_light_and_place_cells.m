function figure_light_and_place_cells(fn,allRs,ccs)

adata = evalin('base','data');
mData = evalin('base','mData');
tdata = evalin('base','dataT');
colors = mData.colors;
sigColor = mData.sigColor;
selAnimals = 1:11;
mData.belt_length = adata{selAnimals(2)}{1}{1}.belt_length;
n = 0;

%%
runthis = 1
if runthis
for kk = 1:3
    for ii = 1:2%length(data)
        for jj = 1:length(selAnimals)
            selCellsC = selectCells15(selAnimals(jj),sprintf('Only_C%d',kk));
            [tempD cnsjj] = getVariableValues(tdata{selAnimals(jj)},'excR',ii);
            exc(ii,kk,jj) = 100*sum(selCellsC & tempD(1,:)')/length(selCellsC);
            [tempD cnsjj] = getVariableValues(tdata{selAnimals(jj)},'inhR',ii);
            inh(ii,kk,jj) = 100*sum(selCellsC & tempD(1,:)')/length(selCellsC);
        end
    end
end

for ii = 1:2
    for jj = 1:3
        tempE = squeeze(exc(ii,jj,:));
        tempI = squeeze(inh(ii,jj,:));
        eAll{jj} = tempE; iAll{jj} = tempI;
    end
    sigE{ii} = significanceTesting(eAll);
    sigI{ii} = significanceTesting(iAll);
end

n = 0;

hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 5 4.5 2.5],'color','w');
hold on;
xdatas = {[1 2 3],[5 6 7]};
for ii = 1:2
    sigR = sigE{ii};
    hs = sigR.anova.multcompare.h; ps = sigR.anova.multcompare.p;
    plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[hs ps],'colors',colors,'sigColor',sigColor,...
        'maxY',0.27,'ySpacing',0.1,'sigTestName','ANOVA','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdatas{ii},'sigFontSize',7);
end

xdatas = {[9 10 11],[13 14 15]};
for ii = 1:2
    sigR = sigE{ii};
    hs = sigR.anova.multcompare.h; ps = sigR.anova.multcompare.p;
    plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[hs ps],'colors',colors,'sigColor',sigColor,...
        'maxY',0.27,'ySpacing',0.1,'sigTestName','ANOVA','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdatas{ii},'sigFontSize',7);
end
xlabel('Exc/Inh/Context'); ylabel('Percentage of Cells');xlim([0 16])
legs = {'Context 1','Context 2','Context 3',[11 0.25 0.25 0.025]};
putLegend(gca,legs,'colors',colors,'sigR',{[],'ks',sigColor,10});
set(gca,'XTick',[2 6 10 14],'XTickLabel',{'Pre-Exc','Post-Exc','Pre-Inh','Post-Inh'});
set(gca,'FontSize',mData.axes_font_size,'FontWeight','Bold','TickDir','out');changePosition(gca,[-0.02 0.01 0.1 -0.05]);
text(1,0.29,'Light Responsive and Place Cells','FontSize',mData.axes_font_size,'FontWeight','Bold');
save_pdf(hf,mData.pdf_folder,'Perc Of light and Place Cells_15',600);
return;
end

%%
runthis = 1
if runthis
for ii = 1:2%length(data)
    for jj = 1:length(selAnimals)
        [tempD cnsjj] = getVariableValues(tdata{selAnimals(jj)},'excR',ii);
        selCellsC = tempD(1,:)';
        exc(jj,ii) = 100*sum(selCellsC)/length(selCellsC);
        [tempD cnsjj] = getVariableValues(tdata{selAnimals(jj)},'inhR',ii);
        selCellsC = tempD(1,:)';
        inh(jj,ii) = 100*sum(selCellsC)/length(selCellsC);
    end
end
% total = exc + inh;
theData{1} = exc(:,1); theData{2} = inh(:,1);
theData{3} = exc(:,2); theData{4} = inh(:,2);
sigR = significanceTesting(theData);

hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 5 1.5 2.5],'color','w');
hold on;
xdatas = {[1 2 3 4],[5 6 7]};

hs = sigR.anova.multcompare.h; ps = sigR.anova.multcompare.p;
plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[hs ps],'colors',{'r','b','r','b'},'sigColor',sigColor,...
    'maxY',8,'ySpacing',0.1,'sigTestName','ANOVA','sigLineWidth',0.25,'BaseValue',0.001,...
    'xdata',xdatas{1},'sigFontSize',9);

% xlabel('Exc/Inh'); 
ylabel('Percentage of Cells');
legs = {'Pre-Exc','Pre-Inh','Post-Exc','Post-Inh',[1 0.25 10 1]};
% putLegend(gca,legs,'colors',colors,'sigR',{[],'ks',sigColor,10});
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'Pre-Exc','Pre-Inh','Post-Exc','Post-Inh'});
xtickangle(45);
set(gca,'FontSize',mData.axes_font_size,'FontWeight','Bold','TickDir','out');changePosition(gca,[-0.02 0.01 0.1 -0.05]);
% text(1,0.29,'Light Responsive and Place Cells','FontSize',mData.axes_font_size,'FontWeight','Bold');
save_pdf(hf,mData.pdf_folder,'Perc Of light Cells_15',600);
return;
end
