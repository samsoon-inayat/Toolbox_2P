function figure1_number_of_PCs

protocol_C = '15_C';
protocol_A = '15_A';
ei_C = evalin('base','ei15_C');
ei_A = evalin('base','ei15_A');
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ET_C = evalin('base',sprintf('ET%s',protocol_C));
ET_A = evalin('base',sprintf('ET%s',protocol_A));
selAnimals_C = 1:length(ei_C)
selAnimals_A = 1:length(ei_A)

% in the following variable all the measurements are in the matrices form
% for each variable colums indicate raster and stim marker types specified 
% the rows indicate condition numbers.
paramMs_C = parameter_matrices('get','15_C');
paramMs_A = parameter_matrices('get','15_A');
% after getting all matrics, we can apply selection criteria to select a
% subgroup of cells
% here is the selection criteria in make_selC_structure function
cellsOrNot = 1; planeNumber = NaN; zMI_Th = 3; fwids = [1 120]; fcens = [0 140]; rs_th = 0.4;
% cellsOrNot = 1; planeNumber = NaN; zMI_Th = 3; fwids = [1 120]; fcens = [0 140]; rs_th = 0.4;
% cellsOrNot = NaN; planeNumber = NaN; zMI_Th = 3; fwids = NaN; fcens = NaN; rs_th = 0.4;
conditionsAndRasterTypes = [31 41 51];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN);
[cpMs_C,pMs_C] = parameter_matrices('select','15_C',{paramMs_C,selC});
[cpMs_A,pMs_A] = parameter_matrices('select','15_A',{paramMs_A,selC});
perc_cells_C = parameter_matrices('print','15_C',{cpMs_C,pMs_C,ET_C,selAnimals_C});
perc_cells_A = parameter_matrices('print','15_A',{cpMs_A,pMs_A,ET_A,selAnimals_A});

all_conds = []; all_rts = [];
for rr = 1:size(pMs_C,1)
    for cc = 1:size(pMs_C,2)
        tcond = conditionsAndRasterTypes(rr,cc);
        nds = dec2base(tcond,10) - '0';
        varNames{rr,cc} = sprintf('C%dR%d',nds(1),nds(2));
        all_conds = [all_conds nds(1)]; all_rts = [all_rts nds(2)];
        xticklabels{cc,rr} = sprintf('%s-%s',paramMs_C.stimMarkers{nds(2)},paramMs_C.rasterTypes{nds(2)}(1));
    end
end
all_conds = unique(all_conds); all_rts = unique(all_rts);
%%
runthis = 1;
if runthis
    numCols = length(all_rts);
    data = perc_cells_C;
    cmdTxt = sprintf('dataT_C = table(');
    for ii = 1:(size(data,2)-1)
        cmdTxt = sprintf('%sdata(:,%d),',cmdTxt,ii);
    end
    cmdTxt = sprintf('%sdata(:,size(data,2)));',cmdTxt);
    eval(cmdTxt);
    data = perc_cells_A;
    cmdTxt = sprintf('dataT_A = table(');
    for ii = 1:(size(data,2)-1)
        cmdTxt = sprintf('%sdata(:,%d),',cmdTxt,ii);
    end
    cmdTxt = sprintf('%sdata(:,size(data,2)));',cmdTxt);
    eval(cmdTxt);
    dataT = [dataT_C;dataT_A]
    dataT.Properties.VariableNames = varNames;
    dataT = [table([ones(length(ei_C),1);2*ones(length(ei_A),1)]) dataT];
    dataT.Properties.VariableNames{1} = 'Group';
    dataT.Group = categorical(dataT.Group)
    colVar1 = [ones(1,numCols) 2*ones(1,numCols) 3*ones(1,numCols)];    colVar2 = [1:numCols 1:numCols 1:numCols 1:numCols];
    within = table(colVar1');
    within.Properties.VariableNames = {'Condition'};
    within.Condition = categorical(within.Condition);
%     within.Raster_Type = categorical(within.Raster_Type);

    rm = fitrm(dataT,sprintf('%s,%s,%s~Group',varNames{1},varNames{2},varNames{3}),'WithinDesign',within,'WithinModel','Condition');
    rtable = ranova(rm,'WithinModel',rm.WithinModel);
%     writetable(dataT,'Percentage_of_PCs_15.xlsx','WriteRowNames',true)
%     writetable(rtable,'RANOVA_Number_of_PCs_15.xlsx','WriteRowNames',true)
    mauchlytbl = mauchly(rm);
    % multcompare(rm,'Day','ComparisonType','bonferroni')
    mcGroup = find_sig_mctbl(multcompare(rm,'Group','By','Condition','ComparisonType','bonferroni'),6);
    mcGroup = multcompare(rm,'Group','By','Condition','ComparisonType','bonferroni');
    writetable(mcGroup,'RANOVA_Number_of_PCs_multcompare_15.xlsx')
%     mcTI = find_sig_mctbl(multcompare(rm,'TI','By','Day','ComparisonType','bonferroni'),6);
%     mcDays = find_sig_mctbl(multcompare(rm,'Day','By','TI','ComparisonType','bonferroni'),6);
    

    [mVarT,semVarT] = findMeanAndStandardError(perc_cells_C);
    [mVarT_A,semVarT_A] = findMeanAndStandardError(perc_cells_A);

    mVar = NaN(1,2*(length(mVarT)));
    semVar = mVar;
    mVar(1:2:6) = mVarT;semVar(1:2:6) = semVarT;
    mVar(2:2:6) = mVarT_A;semVar(2:2:6) = semVarT_A;
    combs = nchoosek(1:6,2); p = ones(size(combs,1),1); h = logical(zeros(size(combs,1),1));
%     row = [1 2]; ii = ismember(combs,row,'rows'); p(ii) = mcGroup{1,6}; h(ii) = 1; 
    % row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
    % row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 

    xdata = [1 2 4 5 7 8]; maxY = 10;
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 2.25 1],'color','w');
    hold on;
    tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3}};
    hbs = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'maxY',maxY,'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',-0.3);
    for ii = 2:2:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+1],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:2:end)+0.5; xticklabels = {'C1','C2','C3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.02 0.03 0.02 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{'Place Cells (%)',[0 0 0]});
    rectangle(gca,'Position',[0.75 9 1 2],'edgecolor','k','facecolor','k');
    text(1.85,10,'Control','FontSize',5);
    rectangle(gca,'Position',[6 9 1 2],'edgecolor','k');
    text(7.2,10,'APP','FontSize',5);
%     text(3.5,18,'Control','FontSize',7);
%     text(18.5,18,'APP','FontSize',7);
    % applyhatch_plusC(gcf
    save_pdf(hf,mData.pdf_folder,sprintf('Percentage of PCs_15'),600);
return;
end


%%
runthis = 0;
if runthis
pcsC = [];
for kk = 1:4
    for jj = 1:length(selAnimals)
        pcsC(jj,kk) = 100*sum(distD{jj,kk})/length(distD{jj,kk});
    end
end

sigR = significanceTesting(pcsC);

hf = figure(10030);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 5 5 3.5],'color','w');
hold on;
xdatas = {[1 2 3 4],[5 6 7]};
hs = sigR.anova.multcompare.h; ps = sigR.anova.multcompare.p;
plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[hs ps],'colors',colors,'sigColor',sigColor,...
    'maxY',30,'ySpacing',1.5,'sigTestName','ANOVA','sigLineWidth',0.25,'BaseValue',0.001,...
    'xdata',xdatas{1},'sigFontSize',12,'sigAsteriskFontSize',17);
xlabel('Condition'); ylabel('Percentage of Cells');
% legs = {'Context 1','Context 2','Context 3',[11 0.25 0.27 0.025]};
% putLegend(gca,legs,'colors',colors,'sigR',{[],'ks',sigColor,10});
set(gca,'XTick',[1:4]);
set(gca,'FontSize',mData.axes_font_size+4,'FontWeight','Bold','TickDir','out');changePosition(gca,[-0.01 0.01 0.06 -0.01]);
% text(1,0.29,'Light Responsive and Place Cells','FontSize',mData.axes_font_size,'FontWeight','Bold');
save_pdf(hf,mData.pdf_folder,'Perc Of Place Cells_contexts_10',600);
return;
end

%%
runthis = 0;
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
runthis = 0;
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