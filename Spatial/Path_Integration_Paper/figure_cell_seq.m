function figure_place_cells_vs_other_cells_1(fn,allRs,ccs)

ei = evalin('base','ei10');
mData = evalin('base','mData');

colors = mData.colors;
sigColor = mData.sigColor;
axes_font_size = mData.axes_font_size;

T = evalin('base','T10.T(selRecs,:)');

selAnimals = [1:13];
% in the following variable all the measurements are in the matrices form
% for each variable colums indicate raster and stim marker types specified 
% the rows indicate condition numbers.
paramMs = parameter_matrices('get');
% after getting all matrics, we can apply selection criteria to select a
% subgroup of cells
% here is the selection criteria in make_selC_structure function
% cellsOrNot = NaN; planeNumber = NaN; zMI_Th = 3; fwids = [0 140]; fcens = [0 140]; rs_th = 0.4;
cellsOrNot = NaN; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN;
conditionsAndRasterTypes = [12 22 32 42]; selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th);
[cpMs,pMs] = parameter_matrices('select',{paramMs,selC});
% parameter_matrices('print percentages',{cpMs,pMs,T,selAnimals});

%%
all_trials = {1:2,3:4,5:7,8:10};%3:10;
all_trials = {1,2,3,4,5,6,7,8,9};%3:10;
sel_trials = [1 2 3 4];
% trials10 = 1;%3:9;
% align cells
stimMarkers = paramMs.stimMarkers;
rasterTypes = paramMs.rasterTypes;
CNi = 3; 
rasterTypeN = 1;

for ani = 1:length(selAnimals)
    an = selAnimals(ani);
    tei = ei(an);
    selCells = pMs{1}.cellSel{an};
    cns = paramMs.all_cns{an};
    a_cell_pos = [];;
for crt = 1:length(conditionsAndRasterTypes)
    all_cellSeq = [];
for si = 1:length(all_trials)
    [ani crt si]
    tcond = conditionsAndRasterTypes(crt);
    Ndigits = dec2base(tcond,10) - '0';
    mRsi = [];
%     for ani = 1:length(selAnimals)
% %         [si ani]
%         an = selAnimals(ani);
%         tei = ei(an);
%         selCells = pMs{1}.cellSel{an};
%         selCells = cpMs.cellSel{an};
%         cns = paramMs.all_cns{an};
        maxDistTime = paramMs.maxDistTime;
        if Ndigits(2) == 2
        [tempD cnso] = getParamValues('fromFrames.sp_rasters',tei,selC.plane_number,Ndigits(1),stimMarkers{Ndigits(2)},rasterTypes{Ndigits(2)},...
        cns(selCells,2:3),maxDistTime);
        else
        [tempD cnso] = getParamValues('rasters',tei,selC.plane_number,Ndigits(1),stimMarkers{Ndigits(2)},rasterTypes{Ndigits(2)},...
            cns(selCells,2:3),maxDistTime);
        end
        if length(tempD) == 0
            continue;
        end
        trials = all_trials{si};
         mR = findMeanRasters(tempD,trials);
        mRsi = [mRsi;mR];
%     end
    [temp,~,~] = getParamValues('',ei(1),1,1,stimMarkers{Ndigits(2)},rasterTypes{Ndigits(2)},'areCells',[Inf Inf]);
    dxs = diff(temp.xs); bin_width = dxs(1); xs = 0:bin_width:1000;
    allRs{si} = mRsi;
    time_xs{si} = xs(1:size(mRsi,2));
    raster_labels{si} = sprintf('Cond - %d, Rast - %d',Ndigits(1),Ndigits(2));
    [~,~,all_cellSeq(si,:)] = findPopulationVectorPlot(allRs{si},[]);
end

cell_pos = [];
for cc = 1:size(all_cellSeq,2)
    for rr = 1:size(all_cellSeq,1)
        cell_pos(cc,rr) = find(all_cellSeq(rr,:) == cc);
    end
end

a_cell_pos = [a_cell_pos cell_pos];
% delta_cell_pos(crt,:) = mean(100*diff(cell_pos,1,2)/size(all_cellSeq,2),2);
% data{crt} = delta_cell_pos(crt,:);

cell_pos = 100*diff(cell_pos,1,2)/size(cell_pos,1);
mean_over_cells(crt,:) = mean(abs(cell_pos));

cell_pos = cell_pos - min(cell_pos(:));
cell_pos = cell_pos/max(cell_pos(:));
all_cell_pos{crt} = cell_pos;
e(crt) = entropy(cell_pos);
f(crt) = BoxCountfracDim(cell_pos);

end
all_mean_over_cells(:,:,ani) = mean_over_cells;
a_cell_pos_an(ani) = {a_cell_pos};
d_cell_pos = mean(abs(100*diff(a_cell_pos,1,2)/size(a_cell_pos,1)));
inds = size(mean_over_cells,2)+1; inds = inds:inds:100; inds = inds(1:(size(mean_over_cells,1)-1));
da_cell_pos_an(ani,:) = d_cell_pos(inds);
a_d_cell_pos(ani,:) = d_cell_pos;
end
result = descriptiveStatistics(all_mean_over_cells,3);
int_env = descriptiveStatistics(da_cell_pos_an);
n = 0;
%%
figure(101);clf;hold on;
for ii = 1:crt
%     plot(1:9,result.avg(ii,:),'color',colors{ii});
    errorbar(1:8,result.avg(ii,:),result.sem(ii,:),'color',colors{ii});
%     subplot(1,crt,ii);
%     imagesc(all_cell_pos{ii});
end
% return;
% sigRM = significanceTesting(mean_over_cells');

n = 0;

%%
% perform repeated measures anova
    numConds = length(conditionsAndRasterTypes);
    numTrialDiffs = length(all_trials)-1;
    colNums = []; acolNums = 1:numTrialDiffs;
    for ii = 1:numConds
        if ii == 1
            colNums = acolNums;
        else
            colNums = [colNums (100+ii-1) acolNums];
        end
    end
    varNames = [];
    for jj = 1:length(colNums)
        varNames{jj} = sprintf('TD%d',jj);
    end
    data = a_d_cell_pos;
    colVar1 = [ones(1,numTrialDiffs) 12 2*ones(1,numTrialDiffs) 23 3*ones(1,numTrialDiffs) 34 4*ones(1,numTrialDiffs)];
    colVar2 = [1:numTrialDiffs 10 1:numTrialDiffs 20 1:numTrialDiffs 30 1:numTrialDiffs];
     within = table(colVar1',colVar2'); within.Properties.VariableNames = {'Condition','TrialDiff'};
    ra = repeatedMeasuresAnova(data,varNames,within);
    rm = ra.rm;
    mcTI = find_sig_mctbl(multcompare(rm,'TrialDiff','By','Condition','ComparisonType','bonferroni'),6);
    mcConds = find_sig_mctbl(multcompare(rm,'Condition','By','TrialDiff','ComparisonType','bonferroni'),6);
    [combs,h,p] = populate_multcomp_h_p(data,within,mcTI,mcConds);
    [mVar semVar] = findMeanAndStandardError(data);
    ds = descriptiveStatistics(data)
    
n = 0;
%%
    axdata = [1:1.5:(10*size(data,2))]; axdata = axdata(1:size(data,2)); maxY = 45;
    xdata = axdata;
    hf = figure(15);clf;set(gcf,'Units','Inches');set(gcf,'Position',[3 3 4 1.5],'color','w');
    hold on;
    ind = 1;
    for ii = 1:numConds
        for jj = 1:numTrialDiffs
            tcolors{ind} = colors{ii};
            ind = ind + 1;
        end
        tcolors{ind} = 'm';
        ind = ind + 1;
    end
    hbs = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'maxY',maxY,'ySpacing',1.25,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.7,'sigLinesStartYFactor',0);
    set(gca,'xlim',[0.25 max(xdata)+.75],'ylim',[20 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata; 
    xticklabels = {'T1','T2','T3','T4','T5','T6','T7','T8','T9'};xticklabels = repmat(xticklabels,1,4);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0.1 0.02 -0.03 -0.011])
    put_axes_labels(gca,{[],[0 0 0]},{{'Percent cell_seq _shift','(z-score)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('cell_seq'),600);
    
    
%%
% perform repeated measures anova
    numCols = size(all_mean_over_cells,2);    numRows = size(all_mean_over_cells,1);
    ind = 1;
    varNames = [];
    for ii = 1:numRows
        for jj = 1:numCols
            varNames{ind} = sprintf('C%dTD%d',ii,jj);
            ind = ind + 1;
        end
    end
    data = [];
    for ii = 1:size(all_mean_over_cells,3)
        thisd = all_mean_over_cells(:,:,ii);
        data(ii,:) = reshape(thisd',1,numRows*numCols);
    end
    colVar1 = [ones(1,numCols) 2*ones(1,numCols) 3*ones(1,numCols) 4*ones(1,numCols)];    colVar2 = [1:numCols 1:numCols 1:numCols 1:numCols];
    within = table(colVar1',colVar2'); within.Properties.VariableNames = {'Condition','TrialDiff'};
    ra = repeatedMeasuresAnova(data,varNames,within);
    rm = ra.rm;
    mcTI = find_sig_mctbl(multcompare(rm,'TrialDiff','By','Condition','ComparisonType','bonferroni'),6);
    mcConds = find_sig_mctbl(multcompare(rm,'Condition','By','TrialDiff','ComparisonType','bonferroni'),6);
    [combs,h,p] = populate_multcomp_h_p(data,within,mcTI,mcConds);
    [mVar semVar] = findMeanAndStandardError(data);
    ds = descriptiveStatistics(data)
%%
    n = 0;
    axdata = [1:1.5:(10*size(data,2))]; axdata = axdata(1:numCols); maxY = 45;
    xdata = []; ixdata = []; offset = 2;
    for ii = 1:numRows
        if ii == 1
            xdata = axdata;
        else
            ixdata(ii-1) = xdata(end) + ((axdata(1) + xdata(end) + offset) - xdata(end))/2;
            xdata = [xdata (axdata + xdata(end) + offset)];
        end
    end
    hf = figure(15);clf;set(gcf,'Units','Inches');set(gcf,'Position',[3 3 4 1.5],'color','w');
    hold on;
    ind = 1;
    for ii = 1:numRows
        for jj = 1:numCols
            tcolors{ind} = colors{ii};
            ind = ind + 1;
        end
    end
    hbs = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'maxY',maxY,'ySpacing',1.25,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.7,'sigLinesStartYFactor',0);
    set(gca,'xlim',[0.25 max(xdata)+.75],'ylim',[20 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata; 
    xticklabels = {'T1','T2','T3','T4','T5','T6','T7','T8'};xticklabels = repmat(xticklabels,1,4);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    hbis = bar(ixdata,int_env.avg,'barWidth',0.05,'BaseValue',0.1,'ShowBaseline','off');
    set(hbis,'FaceColor','m','EdgeColor','m');
    errorbar(ixdata,int_env.avg,int_env.sem,'linestyle', 'none','CapSize',3);
    changePosition(gca,[0.1 0.02 -0.03 -0.011])
    put_axes_labels(gca,{[],[0 0 0]},{{'Percent cell_seq _shift','(z-score)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('cell_seq'),600);
    


%%
if 0
gAllVals = delta_cell_pos(:);
minBin = min(gAllVals);
maxBin = max(gAllVals);
incr = (maxBin-minBin)/20;
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 6 4 3],'color','w');
hold on;
[ha,hb,hca,sigR] = plotDistributions(data,'colors',colors,'maxY',15,'cumPos',[0.5 0.26 0.25 0.5],...
    'min',minBin,'incr',incr,'max',maxBin,'do_mean','No');
hold on;
legs = [];
for ii = 1:length(conditionsAndRasterTypes)
    tcond = conditionsAndRasterTypes(ii);
    Ndigits = dec2base(tcond,10) - '0';
    legs{ii} = sprintf('%d-%d',Ndigits(1),Ndigits(2));
end
ylim([0 15]);
xlims = xlim; dx = xlims(2) - xlims(1); ylims = ylim; dy = ylims(2) - ylims(1);
legs{ii+1} = [xlims(1)+dx/4 dx/30 ylims(1)+dy/1 dy/15];
putLegend(ha,legs,'colors',colors,'sigR',{sigR,'anova',sigColor,10});
axes(ha);
h = xlabel('Mean Diff Pos');%changePosition(h,[0 -dy/3 0]);
h = ylabel('Percentage');changePosition(h,[-0 0 0]);
set(gca,'FontSize',axes_font_size+4,'FontWeight','Bold');changePosition(ha,[0.07 0.01 -0.05 -0.05]);
save_pdf(hf,mData.pdf_folder,sprintf('Distribution Of pos cell seq'),600);
end

%%
if 1
    
end