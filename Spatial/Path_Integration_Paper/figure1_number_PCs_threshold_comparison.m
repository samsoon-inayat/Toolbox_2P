function figure1_number_PCs_threshold_comparison

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
% all cells
conditionsAndRasterTypes = [11;21;31;41];
tcolors = colors;
tcolors = {'k','r'};
%%
if 0
    cellsOrNot = 1; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN; FR = NaN;%[0.1 5000];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN,FR);
out_all = read_data_from_base_workspace_AD(selC)

out = out_all;
ei_C = out.eis{1}; ei_A = out.eis{2};
pMs_C = out.pMs{1}; pMs_A = out.pMs{2};
paramMs_C = out.paramMs{1}; paramMs_A = out.paramMs{2};
selAnimals_C = out.selAnimals{1}; selAnimals_A = out.selAnimals{2};
perc_cells_C = out.perc_cells{1}; perc_cells_A = out.perc_cells{2};

    
numset1 = pMs_C{1}.areCells;
numset2 = pMs_A{1}.areCells;
[h,p,ci,stat] = ttest2(numset1,numset2)
[mC,semC] = findMeanAndStandardError(numset1);
[mA,semA] = findMeanAndStandardError(numset2);
    mVar = [mC mA];
    semVar = [semC semA];
    combs = [1 2]; 
    xdata = 1:2;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.5,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    set(gca,'xlim',[0.25 2.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = [1 2]; xticklabels = {'Control','APP'}; 
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    for ii = 1:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
%     xtickangle(30);
%     rectangle(gca,'Position',[0.75 30.5 1 3.5],'edgecolor','k','facecolor','k');
%     text(1.85,30.5,'Trials','FontSize',5);
%     rectangle(gca,'Position',[6 30.5 1 3.5],'edgecolor','k');
%     text(7.2,30.5,'Inter-Trials','FontSize',5);
    changePosition(gca,[0.1 0.13 -0.2 -0.1])
    put_axes_labels(gca,{[],[0 0 0]},{{'Number of Cells'},[0 0 0]});
    
    save_pdf(hf,mData.pdf_folder,'Figure_Number_Of_cells.pdf',600);
return;
end

%%
if 1
cellsOrNot = 1; planeNumber = NaN; zMI_Th = 1.96; fwids = [1 140]; fcens = [1 150]; rs_th = 0.3; FR = [0.1 5000];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN,FR);
out_196 = read_data_from_base_workspace_AD(selC)
out = out_196;

ei_C = out.eis{1}; ei_A = out.eis{2};
pMs_C = out.pMs{1}; pMs_A = out.pMs{2};
paramMs_C = out.paramMs{1}; paramMs_A = out.paramMs{2};
selAnimals_C = out.selAnimals{1}; selAnimals_A = out.selAnimals{2};
perc_cells_C = out.perc_cells{1}; perc_cells_A = out.perc_cells{2};

    
numset1 = pMs_C{1}.areCells;
numset2 = pMs_A{1}.areCells;
[h,p,ci,stat] = ttest2(numset1,numset2)
[mC,semC] = findMeanAndStandardError(numset1);
[mA,semA] = findMeanAndStandardError(numset2);
    mVar = [mC mA];
    semVar = [semC semA];
    combs = [1 2]; 
    xdata = 1:2;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.5,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    set(gca,'xlim',[0.25 2.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = [1 2]; xticklabels = {'Control','APP'}; 
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    for ii = 1:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
%     xtickangle(30);
%     rectangle(gca,'Position',[0.75 30.5 1 3.5],'edgecolor','k','facecolor','k');
%     text(1.85,30.5,'Trials','FontSize',5);
%     rectangle(gca,'Position',[6 30.5 1 3.5],'edgecolor','k');
%     text(7.2,30.5,'Inter-Trials','FontSize',5);
    changePosition(gca,[0.1 0.13 -0.2 -0.1])
    put_axes_labels(gca,{[],[0 0 0]},{{'Number of Cells'},[0 0 0]});
    
    save_pdf(hf,mData.pdf_folder,'Figure_Number_Of_cells.pdf',600);
return;
end


%%
%%
if 1
cellsOrNot = 1; planeNumber = NaN; zMI_Th = 3; fwids = [1 140]; fcens = [1 150]; rs_th = 0.3; FR = [0.1 5000];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN,FR);
out_3 = read_data_from_base_workspace_AD(selC)
out = out_3;

ei_C = out.eis{1}; ei_A = out.eis{2};
pMs_C = out.pMs{1}; pMs_A = out.pMs{2};
paramMs_C = out.paramMs{1}; paramMs_A = out.paramMs{2};
selAnimals_C = out.selAnimals{1}; selAnimals_A = out.selAnimals{2};
perc_cells_C = out.perc_cells{1}; perc_cells_A = out.perc_cells{2};

    
numset1 = pMs_C{1}.areCells;
numset2 = pMs_A{1}.areCells;
[h,p,ci,stat] = ttest2(numset1,numset2)
[mC,semC] = findMeanAndStandardError(numset1);
[mA,semA] = findMeanAndStandardError(numset2);
    mVar = [mC mA];
    semVar = [semC semA];
    combs = [1 2]; 
    xdata = 1:2;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.5,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    set(gca,'xlim',[0.25 2.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = [1 2]; xticklabels = {'Control','APP'}; 
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    for ii = 1:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
%     xtickangle(30);
%     rectangle(gca,'Position',[0.75 30.5 1 3.5],'edgecolor','k','facecolor','k');
%     text(1.85,30.5,'Trials','FontSize',5);
%     rectangle(gca,'Position',[6 30.5 1 3.5],'edgecolor','k');
%     text(7.2,30.5,'Inter-Trials','FontSize',5);
    changePosition(gca,[0.1 0.13 -0.2 -0.1])
    put_axes_labels(gca,{[],[0 0 0]},{{'Number of Cells'},[0 0 0]});
    
    save_pdf(hf,mData.pdf_folder,'Figure_Number_Of_cells.pdf',600);
return;
end




% conditionsAndRasterTypes = [11 21 31 41];

