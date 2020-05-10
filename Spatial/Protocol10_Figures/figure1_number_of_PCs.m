function figure1_Distributions

adata = evalin('base','data');
mData = evalin('base','mData');
colors = mData.colors;
sigColor = mData.sigColor;
% allCells = mData.allCells;
% selAnimals = 6:11;
% selAnimals = 1:3;
selAnimals = 1:8;
% selAnimals = 5:8;
n = 0;
%%
for jj = 1:length(selAnimals)
    for ii = 1:4%length(data)
         [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'placeCells3',ii);
         distD{jj,ii} = tempD;
         cns{jj,ii} = cnsjj;
         pcs(jj,ii) = sum(tempD);
         [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'placeCells4',ii);
         distD4{jj,ii} = tempD;
         cns4{jj,ii} = cnsjj;
         pcs4(jj,ii) = sum(tempD);
         [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'placeCells5',ii);
         distD5{jj,ii} = tempD;
         cns5{jj,ii} = cnsjj;
         pcs5(jj,ii) = sum(tempD);
    end
    numCells(jj,1) = length(tempD);
end

for jj = 1:length(selAnimals)
    allpcsU = distD{jj,1}; allpcsU4 = distD4{jj,1}; allpcsU5 = distD5{jj,1};
    for ii = 2:4%length(data)
         allpcsU = allpcsU | distD{jj,ii}; allpcsU4 = allpcsU4 | distD4{jj,ii}; allpcsU5 = allpcsU5 | distD5{jj,ii};
         lastPCs = distD{jj,ii-1}; currentPCs = distD{jj,ii};
         remained(jj,ii) = 100*sum(lastPCs & currentPCs)/sum(lastPCs);
         disrupted(jj,ii) = 100*sum(lastPCs & ~currentPCs)/sum(lastPCs);
         newones(jj,ii) = 100*sum(~lastPCs & currentPCs)/sum(currentPCs);
    end
    percent_pcs(jj,1) = 100*sum(allpcsU)/length(distD{jj,1});
    percent_pcsi(jj,:) = 100*pcs(jj,:)./length(distD{jj,1});
    upcs(jj,1) = sum(allpcsU); upcs4(jj,1) = sum(allpcsU4); upcs5(jj,1) = sum(allpcsU5);
end

change_percent_pcsi = percent_pcsi - percent_pcsi(:,1);

%%
runthis = 0;
if runthis
hf = figure(1002);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 5 6 2.5],'color','w');
hold on;
barVar = [numCells upcs upcs4 upcs5];
hb = bar(barVar);
xlabel('Animal Number'); ylabel('Number of Cells');
hleg = legend(...
    sprintf('All Cells (N = %d)',sum(numCells)),...
    sprintf('Place Cells (zMI>3, N = %d - %.0f%%)',sum(upcs),100*sum(upcs)/sum(numCells)),...
    sprintf('Place Cells (zMI>4, N = %d - %.0f%%)',sum(upcs4),100*sum(upcs4)/sum(numCells)),...
    sprintf('Place Cells (zMI>5, N = %d - %.0f%%)',sum(upcs5),100*sum(upcs5)/sum(numCells)));
    changePosition(hleg,[0.07 0.08 0 0]);
legend boxoff
for ii = 1:4
    set(hb(ii),'FaceColor',colors{ii},'EdgeColor',colors{ii});
end
for ii = 1:size(numCells,1)
    text(ii,upcs(ii)+50,sprintf('%.0f%%',100*upcs5(ii)/numCells(ii)),'Color',colors{4});
end
set(gca,'XTick',[1:8],'XTickLabel',{'1','2','3','4','5','6','7','8'});
set(gca,'FontSize',mData.axes_font_size,'FontWeight','Bold','TickDir','out');changePosition(gca,[-0.02 0.01 0.1 -0.05]);
save_pdf(hf,mData.pdf_folder,'Number Of Place Cells',600);
return
end
%%
runthis = 1;
if runthis
pcsC = [];
for kk = 1:4
    for jj = 1:length(selAnimals)
        selCellsC = selectCells10(selAnimals(jj),sprintf('Only_C%d',kk),adata);
        pcsC(jj,kk) = 100*sum(selCellsC)/length(selCellsC);
    end
end

for jj = 1:4
    theData{jj} = pcsC(:,jj);
end

sigR = significanceTesting(theData);

hf = figure(10030);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 5 5 3.5],'color','w');
hold on;
xdatas = {[1 2 3 4],[5 6 7]};
hs = sigR.anova.multcompare.h; ps = sigR.anova.multcompare.p;
plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[hs ps],'colors',colors,'sigColor',sigColor,...
    'maxY',20,'ySpacing',1.5,'sigTestName','ANOVA','sigLineWidth',0.25,'BaseValue',0.001,...
    'xdata',xdatas{1},'sigFontSize',12,'sigAsteriskFontSize',17);
xlabel('Context'); ylabel('Percentage of Cells');
% legs = {'Context 1','Context 2','Context 3',[11 0.25 0.27 0.025]};
% putLegend(gca,legs,'colors',colors,'sigR',{[],'ks',sigColor,10});
set(gca,'XTick',[1:4]);
set(gca,'FontSize',mData.axes_font_size+4,'FontWeight','Bold','TickDir','out');changePosition(gca,[-0.01 0.01 0.06 -0.01]);
% text(1,0.29,'Light Responsive and Place Cells','FontSize',mData.axes_font_size,'FontWeight','Bold');
save_pdf(hf,mData.pdf_folder,'Perc Of Place Cells_contexts_10',600);
return;
end

%%
runthis = 1;
if runthis
theData = [];
for jj = 1:3
    theData{jj} = remained(:,jj+1);
end
sigR = significanceTesting(theData);
hf = figure(1005);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 5 1.5 2.5],'color','w');
hold on;
xdatas = {[1 2 3],[5 6 7]};
hs = sigR.anova.multcompare.h; ps = sigR.anova.multcompare.p;
plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[hs ps],'colors',colors,'sigColor',sigColor,...
    'maxY',120,'ySpacing',1.5,'sigTestName','ANOVA','sigLineWidth',0.25,'BaseValue',0.001,...
    'xdata',xdatas{1},'sigFontSize',9,'sigAsteriskFontSize',15);
xlabel('Contexts'); ylabel('Percentage of Cells');
% legs = {'Context 1','Context 2','Context 3',[11 0.25 0.27 0.025]};
% putLegend(gca,legs,'colors',colors,'sigR',{[],'ks',sigColor,10});
set(gca,'XTick',[1:4]);xlim([0 5]);
set(gca,'FontSize',mData.axes_font_size,'FontWeight','Bold','TickDir','out');changePosition(gca,[-0.01 0.01 0.06 -0.01]);
% text(1,0.29,'Light Responsive and Place Cells','FontSize',mData.axes_font_size,'FontWeight','Bold');
save_pdf(hf,mData.pdf_folder,'Perc Of Place Cells_Remained_contexts_10',600);
return;
end

%%
runthis = 1;
if runthis
theData = [];
for jj = 1:3
    theData{jj} = disrupted(:,jj+1);
end
sigR = significanceTesting(theData);
hf = figure(1006);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 5 1.5 2.5],'color','w');
hold on;
xdatas = {[1 2 3],[5 6 7]};
hs = sigR.anova.multcompare.h; ps = sigR.anova.multcompare.p;
plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[hs ps],'colors',colors,'sigColor',sigColor,...
    'maxY',120,'ySpacing',1.5,'sigTestName','ANOVA','sigLineWidth',0.25,'BaseValue',0.001,...
    'xdata',xdatas{1},'sigFontSize',9,'sigAsteriskFontSize',15);
xlabel('Contexts'); ylabel('Percentage of Cells');
% legs = {'Context 1','Context 2','Context 3',[11 0.25 0.27 0.025]};
% putLegend(gca,legs,'colors',colors,'sigR',{[],'ks',sigColor,10});
set(gca,'XTick',[1:4]);xlim([0 5])
set(gca,'FontSize',mData.axes_font_size,'FontWeight','Bold','TickDir','out');changePosition(gca,[-0.01 0.01 0.06 -0.01]);
% text(1,0.29,'Light Responsive and Place Cells','FontSize',mData.axes_font_size,'FontWeight','Bold');
save_pdf(hf,mData.pdf_folder,'Perc Of Place Cells_Disrupted_contexts_10',600);
return;
end

%%
runthis = 1;
if runthis
theData = [];
for jj = 1:3
    theData{jj} = newones(:,jj+1);
end
sigR = significanceTesting(theData);
hf = figure(1007);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 5 1.5 2.5],'color','w');
hold on;
xdatas = {[2 3 4],[5 6 7]};
hs = sigR.anova.multcompare.h; ps = sigR.anova.multcompare.p;
plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[hs ps],'colors',colors(2:end),'sigColor',sigColor,...
    'maxY',120,'ySpacing',1.5,'sigTestName','ANOVA','sigLineWidth',0.25,'BaseValue',0.001,...
    'xdata',xdatas{1},'sigFontSize',9,'sigAsteriskFontSize',15);
xlabel('Contexts'); ylabel('Percentage of Cells');
% legs = {'Context 1','Context 2','Context 3',[11 0.25 0.27 0.025]};
% putLegend(gca,legs,'colors',colors,'sigR',{[],'ks',sigColor,10});
set(gca,'XTick',[1:4],'XTickLabel',num2cell(1:4'));
set(gca,'FontSize',mData.axes_font_size,'FontWeight','Bold','TickDir','out');changePosition(gca,[-0.01 0.01 0.06 -0.01]);
% text(1,0.29,'Light Responsive and Place Cells','FontSize',mData.axes_font_size,'FontWeight','Bold');
save_pdf(hf,mData.pdf_folder,'Perc Of Place Cells_NewOnes_contexts_10',600);
return;
end