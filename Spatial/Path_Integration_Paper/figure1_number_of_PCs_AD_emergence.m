function figure1_number_of_PCs_AD_emergence

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
cellsOrNot = 1; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN;
cellsOrNot = 1; planeNumber = NaN; zMI_Th = 1.96; fwids = [1 150]; fcens = [1 150]; rs_th = 0.3;
cellsOrNot = 1; planeNumber = NaN; zMI_Th = 3; fwids = [1 120]; fcens = [0 140]; rs_th = 0.4; HaFD_th = NaN; HiFD_th = NaN;
conditionsAndRasterTypes = [11 21 31 41];
% conditionsAndRasterTypes = [11 -41];
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
    for cc = 1:4
        % common cells
        cmdTxt = sprintf('conditionsAndRasterTypes = [%d1 %d1];',rr,cc); eval(cmdTxt);
        selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN);
        out = read_data_from_base_workspace_AD(selC);
        common_cells_C(rr,cc,:) = out.cpMs{1}.numCells; common_cells_A(rr,cc,:) = out.cpMs{2}.numCells;
        
        % new cells
        cmdTxt = sprintf('conditionsAndRasterTypes = [-%d1 %d1];',rr,cc); eval(cmdTxt);
        selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN);
        out = read_data_from_base_workspace_AD(selC);
        new_cells_C(rr,cc,:) = out.cpMs{1}.numCells; new_cells_A(rr,cc,:) = out.cpMs{2}.numCells;
        
        % disrupted cells
        cmdTxt = sprintf('conditionsAndRasterTypes = [%d1 -%d1];',rr,cc); eval(cmdTxt);
        selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN);
        out = read_data_from_base_workspace_AD(selC);
        disrupted_cells_C(rr,cc,:) = out.cpMs{1}.numCells; disrupted_cells_A(rr,cc,:) = out.cpMs{2}.numCells;
    end
end
n = 0;

for an = 1:size(common_cells_C,3)
    for rr = 1:4
        for cc = 1:4
            perc_common_cells_C(rr,cc,an) = 100*common_cells_C(rr,cc,an)/num_cells_C(an,rr);
            perc_new_cells_C(rr,cc,an) = 100*new_cells_C(rr,cc,an)/num_cells_C(an,cc);
            perc_disrupted_cells_C(rr,cc,an) = 100*disrupted_cells_C(rr,cc,an)/num_cells_C(an,rr);
        end
    end
end

mask = ones(4,4); mask = triu(mask,1) & ~triu(mask,2);

for an = 1:size(common_cells_C,3)
    temp = perc_common_cells_C(:,:,an); p_c_C(an,:) = temp(mask)';
    temp = perc_new_cells_C(:,:,an); p_n_C(an,:) = temp(mask)';
    temp = perc_disrupted_cells_C(:,:,an); p_d_C(an,:) = temp(mask)';
end

for an = 1:size(common_cells_A,3)
    for rr = 1:4
        for cc = 1:4
            perc_common_cells_A(rr,cc,an) = 100*common_cells_A(rr,cc,an)/num_cells_A(an,rr);
            perc_new_cells_A(rr,cc,an) = 100*new_cells_A(rr,cc,an)/num_cells_A(an,cc);
            perc_disrupted_cells_A(rr,cc,an) = 100*disrupted_cells_A(rr,cc,an)/num_cells_A(an,rr);
        end
    end
end

mask = ones(4,4); mask = triu(mask,1) & ~triu(mask,2);

for an = 1:size(common_cells_A,3)
    temp = perc_common_cells_A(:,:,an); p_c_A(an,:) = temp(mask)';
    temp = perc_new_cells_A(:,:,an); p_n_A(an,:) = temp(mask)';
    temp = perc_disrupted_cells_A(:,:,an); p_d_A(an,:) = temp(mask)';
end
n=0
%%
selType = 3;
type_cells = {'common','new','disrupted'};
if selType == 1
    perc_cells_C = p_c_C; perc_cells_A = p_c_A;
end
if selType == 2
    perc_cells_C = p_n_C; perc_cells_A = p_n_A;
end
if selType == 3
    perc_cells_C = p_d_C; perc_cells_A = p_d_A;
end


if 1
    numCols = 3;%length(all_rts);
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
    dataT.Properties.VariableNames = {'CC1','CC2','CC3'};
    dataT = [table([ones(length(ei_C),1);2*ones(length(ei_A),1)]) dataT];
    dataT.Properties.VariableNames{1} = 'Group';
    dataT.Group = categorical(dataT.Group)
    
    colVar1 = [1 2 3];
    within = table(colVar1');
    within.Properties.VariableNames = {'Condition'};
    within.Condition = categorical(within.Condition);
    ra = repeatedMeasuresAnova(dataT,within);

    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
     xdata = [1 2 3 5:7]; 
%     xdata = [1 2 3 4];
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 2 1],'color','w');
    hold on;
    tcolors ={colors{1};colors{2};colors{3};colors{1};colors{2};colors{3};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',20,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
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
    save_pdf(hf,mData.pdf_folder,sprintf('Percentage of %s PCs',type_cells{selType}),600);
return;
end


