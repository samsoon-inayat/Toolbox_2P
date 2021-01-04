function figure1_number_of_PCs

protocol_C = '10_C';
% protocol_A = '10_A';
ei_C = evalin('base','ei10_C');
ei_A = evalin('base','ei10_A');
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ET_C = evalin('base',sprintf('ET%s',protocol_C));
% ET_A = evalin('base',sprintf('ET%s',protocol_A));
selAnimals_C = 1:length(ei_C)
selAnimals_A = 1:length(ei_A)

% in the following variable all the measurements are in the matrices form
% for each variable colums indicate raster and stim marker types specified 
% the rows indicate condition numbers.
paramMs_C = parameter_matrices('get','10_CD_Ctrl');
paramMs_A = parameter_matrices('get','10_CC_Ctrl');
% after getting all matrics, we can apply selection criteria to select a
% subgroup of cells
% here is the selection criteria in make_selC_structure function
% cellsOrNot = 1; planeNumber = NaN; zMI_Th = 3; fwids = [1 120]; fcens = [0 140]; rs_th = 0.4;
cellsOrNot = 1; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN;
conditionsAndRasterTypes = [11 21 31 41];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN);
[cpMs_C,pMs_C] = parameter_matrices('select','10_C',{paramMs_C,selC});
[cpMs_A,pMs_A] = parameter_matrices('select','10_A',{paramMs_A,selC});
% perc_cells_C = parameter_matrices('print','10_C',{cpMs_C,pMs_C,ET_C,selAnimals_C});
% perc_cells_A = parameter_matrices('print','10_A',{cpMs_A,pMs_A,ET_A,selAnimals_A});

all_conds = []; all_rts = [];
gAllVals_C = []; gAllVals_A = [];
for rr = 1:size(pMs_C,1)
    for cc = 1:size(pMs_C,2)
        tcond = conditionsAndRasterTypes(rr,cc);
        nds = dec2base(tcond,10) - '0';
        varNames{rr,cc} = sprintf('C%dR%d',nds(1),nds(2));
        all_conds = [all_conds nds(1)]; all_rts = [all_rts nds(2)];
        xticklabels{cc,rr} = sprintf('%s-%s',paramMs_C.stimMarkers{nds(2)},paramMs_C.rasterTypes{nds(2)}(1));
        for an = 1:length(selAnimals_C)
            zMIs_C(an,rr,cc) = nanmean(squeeze(pMs_C{rr,cc}.all_zMIs{selAnimals_C(an)}(nds(1),nds(2),:)));
            a_zMIs_C{an,rr,cc} = squeeze(pMs_C{rr,cc}.all_zMIs{selAnimals_C(an)}(nds(1),nds(2),:));
            gAllVals_C = [gAllVals_C;a_zMIs_C{an,rr,cc}];
        end
    end
end
for rr = 1:size(pMs_A,1)
    for cc = 1:size(pMs_A,2)
        tcond = conditionsAndRasterTypes(rr,cc);
        nds = dec2base(tcond,10) - '0';
        for an = 1:length(selAnimals_A)
            zMIs_A(an,rr,cc) = nanmean(squeeze(pMs_A{rr,cc}.all_zMIs{selAnimals_A(an)}(nds(1),nds(2),:)));
            a_zMIs_A{an,rr,cc} = squeeze(pMs_A{rr,cc}.all_zMIs{selAnimals_A(an)}(nds(1),nds(2),:));
            gAllVals_A = [gAllVals_A;a_zMIs_A{an,rr,cc}];
        end
    end
end
all_conds = unique(all_conds); all_rts = unique(all_rts);
var_oi_A = squeeze(zMIs_A);
var_oi_C = squeeze(zMIs_C);
n = 0;
%%
runthis = 1;
if runthis
    numCols = length(all_rts);
    data = var_oi_C;
    cmdTxt = sprintf('dataT_C = table(');
    for ii = 1:(size(data,2)-1)
        cmdTxt = sprintf('%sdata(:,%d),',cmdTxt,ii);
    end
    cmdTxt = sprintf('%sdata(:,size(data,2)));',cmdTxt);
    eval(cmdTxt);
    data = var_oi_A;
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
    
    colVar1 = [ones(1,numCols) 2*ones(1,numCols) 3*ones(1,numCols) 4*ones(1,numCols)];    colVar2 = [1:numCols 1:numCols 1:numCols 1:numCols];
    within = table(colVar1');
    within.Properties.VariableNames = {'Condition'};
    within.Condition = categorical(within.Condition);
    ra = repeatedMeasuresAnova(dataT,within);
%%
    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
%     row = [1 2]; ii = ismember(combs,row,'rows'); p(ii) = mcGroup{1,6}; h(ii) = 1; 
    % row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
    % row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 

    xdata = [1 2 3 4 6:9]; maxY = 10;
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
%     tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4}};
    tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'maxY',maxY,'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 5:length(hbs)
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+1],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = {'C1','C2','C3','C4'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.03 0.03 0.04 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{'Average zMI',[0 0 0]});
    
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bargraph',mfilename),600);
return;
end

%%
runthis = 1;
if runthis
    data_C = squeeze(a_zMIs_C);
    data_A = squeeze(a_zMIs_A);
    ii = 4; data = [data_C(:,ii) data_A(:,ii)];
    minBin = min([gAllVals_C;gAllVals_A]);
    maxBin = 15;%max(gAllVals);
    incr = (maxBin-minBin)/100;
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 6 2 1.5],'color','w');
    hold on;
    [ha,hb,hca,sigR] = plotDistributions(data,'colors',colors,'maxY',90,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    hold on;
    legs = {sprintf('%s-C',varNames{ii}),sprintf('%s-A',varNames{ii})};
%     paramMs = paramMs_A;
%     for ii = 1:length(varNames)
%         if ii == 1
%             legs{ii} = sprintf('%s-%s (N = %d)',paramMs.stimMarkers{ii},paramMs.rasterTypes{ii}(1),size(data,1));
%         else
%             legs{ii} = sprintf('%s-%s',paramMs.stimMarkers{ii},paramMs.rasterTypes{ii}(1));
%         end
%     end
    ylim([0 110]);
    xlim([-3 12])
    xlims = xlim; dx = xlims(2) - xlims(1); ylims = ylim; dy = ylims(2) - ylims(1);
    legs{length(legs)+1} = [xlims(1)+dx/1.5 dx/30 ylims(1)+dy/3 dy/15];
    putLegend(ha,legs,'colors',colors);
    axes(ha);
    h = xlabel('Mutual Information Z-Score');%changePosition(h,[0 -dy/3 0]);
    h = ylabel('Percentage');changePosition(h,[-0.1 0 0]);
    set(gca,'FontSize',axes_font_size,'FontWeight','Bold');changePosition(ha,[0.07 0.01 -0.05 -0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution Of zMI %d',ii),600);
return;
end