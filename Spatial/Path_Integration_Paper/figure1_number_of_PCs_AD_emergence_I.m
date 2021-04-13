function figure1_number_of_PCs_AD_emergence_I

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
cellsOrNot = 1; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN;
cellsOrNot = 1; planeNumber = NaN; zMI_Th = 1.96; fwids = [1 150]; fcens = [1 150]; rs_th = 0.3;
conditionsAndRasterTypes = [11 21 31 41];
% conditionsAndRasterTypes = [11 21 31 41];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN);
out = read_data_from_base_workspace_AD(selC)

ei_C = out.eis{1}; ei_A = out.eis{2};
pMs_C = out.pMs{1}; pMs_A = out.pMs{2};
paramMs_C = out.paramMs{1}; paramMs_A = out.paramMs{2};
cpMs_C = out.cpMs{1}; cpMs_A = out.cpMs{2};
selAnimals_C = out.selAnimals{1}; selAnimals_A = out.selAnimals{2};
perc_cells_C = out.perc_cells{1}; perc_cells_A = out.perc_cells{2};

for ii = 1:4
    num_cells_C(:,ii) = pMs_C{ii}.numCells;  num_cells_A(:,ii) = pMs_A{ii}.numCells;
end

for rr = 1:4
    % common cells
    cmdTxt = sprintf('conditionsAndRasterTypes = [%d5];',rr); eval(cmdTxt);
    cellsOrNot = 1; planeNumber = NaN; zMI_Th = 1.96; fwids = NaN; fcens = NaN; rs_th = NaN;
    selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN);
    out = read_data_from_base_workspace_AD(selC);
    temp5_C = out.pMs{1}{1}.cellSel; temp5_A = out.pMs{2}{1}.cellSel;
    
    cmdTxt = sprintf('conditionsAndRasterTypes = [%d1];',rr); eval(cmdTxt);
    cellsOrNot = 1; planeNumber = NaN; zMI_Th = 1.96; fwids = [1 150]; fcens = [1 150]; rs_th = 0.3;
    selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN);
    out = read_data_from_base_workspace_AD(selC);
    temp1_C = out.pMs{1}{1}.cellSel; temp1_A = out.pMs{2}{1}.cellSel;
    
    for an = 1:length(temp1_C)
        c_C(an,rr) = 100*sum(temp5_C{an} & temp1_C{an})/num_cells_C(an,rr);
    end
    
    for an = 1:length(temp1_A)
        c_A(an,rr) = 100*sum(temp5_A{an} & temp1_A{an})/num_cells_A(an,rr);
    end
end
n = 0;

%%
selType = 3;
type_cells = {'common','new','disrupted'};
if selType == 1
    perc_cells_C = c_C; perc_cells_A = c_A;
end


if 1
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
    dataT.Properties.VariableNames = {'C1','C2','C3','C4'};
    dataT = [table([ones(length(ei_C),1);2*ones(length(ei_A),1)]) dataT];
    dataT.Properties.VariableNames{1} = 'Group';
    dataT.Group = categorical(dataT.Group)
    
    colVar1 = [1 2 3 4];
    within = table(colVar1');
    within.Properties.VariableNames = {'Condition'};
    within.Condition = categorical(within.Condition);
    ra = repeatedMeasuresAnova(dataT,within);

    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
     xdata = [1 2 3 4 6:9]; 
%     xdata = [1 2 3 4];
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 2 1],'color','w');
    hold on;
    tcolors ={colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 4:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C12','C23','C34','C12','C23','C34'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0.075 0.0 0 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Percentage of',sprintf('%s cells',type_cells{selType})},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Percentage of %s PCs Inter',type_cells{selType}),600);
return;
end


