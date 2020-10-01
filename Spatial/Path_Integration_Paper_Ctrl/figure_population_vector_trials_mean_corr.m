function figure_place_cells_vs_other_cells_1(fn,allRs,ccs)

ei = evalin('base','ei10');
mData = evalin('base','mData');
T = evalin('base','T10.T(selRecs,:)');
colors = mData.colors

selAnimals = [1:5 7 9 11:13];
paramMs = parameter_matrices('get');
%%
a_conditionsAndRasterTypes = [12 22 32 42];
for crt = 1:length(a_conditionsAndRasterTypes)
    meanC = [];
    for ii = 1:length(selAnimals)
        [crt ii]
        meanC(ii,:) = find_mean_corr(a_conditionsAndRasterTypes(crt),selAnimals(ii),paramMs,ei,mData,T);
    end
    a_meanC(crt,:,:) = meanC';
end
result = descriptiveStatistics(meanC);
n = 0;
%%
% perform repeated measures anova
    all_mean_over_cells = a_meanC;
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
    axdata = [1:1.5:(10*size(data,2))]; axdata = axdata(1:numCols); maxY = 0.24;
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
        'maxY',maxY,'ySpacing',0.01,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.7,'sigLinesStartYFactor',0);
    set(gca,'xlim',[0.25 max(xdata)+.75],'ylim',[0.02 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata; 
    xticklabels = {'T1','T2','T3','T4','T5','T6','T7','T8','T9','T10'};xticklabels = repmat(xticklabels,1,4);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
   
%     changePosition(gca,[0.1 0.02 -0.03 -0.011])
    put_axes_labels(gca,{[],[0 0 0]},{{'Percent cell_seq _shift','(z-score)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('cell_seq'),600);




%%
figure(101);clf;hold on;
errorbar(1:9,result.avg(1,:),result.sem(1,:),'color',colors{1});
sigR = significanceTesting(meanC);

n = 0;

function meanC = find_mean_corr(conditionsAndRasterTypes,selAnimals,paramMs,ei,mData,T)

cellsOrNot = NaN; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN;
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th);
[cpMs,pMs] = parameter_matrices('select',{paramMs,selC});

all_trials = {1:2,3:4,5:7,8:10};%3:10;
tcond = conditionsAndRasterTypes(1);
Ndigits = dec2base(tcond,10) - '0';
if Ndigits(2) == 1 || Ndigits(2) == 3 || Ndigits(2) == 2
    all_trials = {1,2,3,4,5,6,7,8,9,10 };%3:10;
    sel_trials = [1 2 3 4 5 6 7 8 9,10];
else
    all_trials = {1,2,3,4,5,6,7,8,9};%3:10;
    sel_trials = [1 2 3 4 5 6 7 8 9];
end
% trials10 = 1;%3:9;
% align cells
stimMarkers = paramMs.stimMarkers;
rasterTypes = paramMs.rasterTypes;
CNi = 3; 
rasterTypeN = 1;


for si = 1:length(all_trials)
    tcond = conditionsAndRasterTypes(1);
    Ndigits = dec2base(tcond,10) - '0';
    mRsi = [];
    for ani = 1:length(selAnimals)
%         [si ani]
        an = selAnimals(ani);
        tei = ei(an);
        selCells = pMs{1}.cellSel{an};
%         selCells = cpMs.cellSel{an};
        cns = paramMs.all_cns{an};
        maxDistTime = paramMs.maxDistTime;
        if Ndigits(2) == 1 || Ndigits(2) == 3
            [tempD cnso] = getParamValues('rasters',tei,selC.plane_number,Ndigits(1),stimMarkers{Ndigits(2)},rasterTypes{Ndigits(2)},...
            cns(selCells,2:3),maxDistTime);
        else
            [tempD cnso] = getParamValues('fromFrames.sp_rasters',tei,selC.plane_number,Ndigits(1),stimMarkers{Ndigits(2)},rasterTypes{Ndigits(2)},...
            cns(selCells,2:3),maxDistTime);
            [temp_rst] = getParamValues('gauss_fit_on_mean.rst',tei,selC.plane_number,Ndigits(1),stimMarkers{Ndigits(2)},rasterTypes{Ndigits(2)},...
            cns(selCells,2:3),maxDistTime);
        end
        if length(tempD) == 0
            continue;
        end
        trials = all_trials{si};
        mR = findMeanRasters(tempD,trials);
        mRsi = [mRsi;mR];
    end
    [temp,~,~] = getParamValues('',ei(1),1,1,stimMarkers{Ndigits(2)},rasterTypes{Ndigits(2)},'areCells',[Inf Inf]);
    dxs = diff(temp.xs); bin_width = dxs(1); xs = 0:bin_width:1000;
    allRs{si} = mRsi;
    time_xs{si} = xs(1:size(mRsi,2));
    raster_labels{si} = sprintf('Cond - %d, Rast - %d',Ndigits(1),Ndigits(2));
    [~,~,all_cellSeq{si}] = findPopulationVectorPlot(allRs{si},[]);
end

[~,~,cellNums] = findPopulationVectorPlot(allRs{end},[]);
for ii = 1:length(all_trials(sel_trials))
    mRsi = allRs{sel_trials(ii)};
    [allP{ii},allC{ii}] = findPopulationVectorPlot(mRsi,[],cellNums);
    thisC = allC{ii}; cinds = triu(true(size(thisC)),1);% & tril(true(size(thisC)),5);
    meanC(ii) = nanmean(thisC(cinds)); 
end

n = 0;
