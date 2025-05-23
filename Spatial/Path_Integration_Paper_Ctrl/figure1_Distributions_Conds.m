function figure1_Distributions

rsel = 1
ei_names = {'ei10_C','ei10_A','ei_comb'};
param_names = {'10_CD_Ctrl','10_CC_Ctrl','10_C_Comb'};

ei_C = evalin('base',sprintf('%s',ei_names{rsel}));
% ei_A = evalin('base','ei10_A');
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
selAnimals_C = 1:length(ei_C)
% selAnimals_A = 1:length(ei_A)

% in the following variable all the measurements are in the matrices form
% for each variable colums indicate raster and stim marker types specified 
% the rows indicate condition numbers.
paramMs_C = parameter_matrices('get',sprintf('%s',param_names{rsel}));
% paramMs_A = parameter_matrices('get','10_A');
% after getting all matrics, we can apply selection criteria to select a
% subgroup of cells
% here is the selection criteria in make_selC_structure function
% cellsOrNot = 1; planeNumber = NaN; zMI_Th = 3; fwids = [1 120]; fcens = [0 140]; rs_th = 0.4;
cellsOrNot = 1; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN;
conditionsAndRasterTypes = [11 12];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN);
[cpMs_C,pMs_C] = parameter_matrices('select','10_C',{paramMs_C,selC});
% [cpMs_A,pMs_A] = parameter_matrices('select','10_A',{paramMs_A,selC});
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
% for rr = 1:size(pMs_A,1)
%     for cc = 1:size(pMs_A,2)
%         tcond = conditionsAndRasterTypes(rr,cc);
%         nds = dec2base(tcond,10) - '0';
%         for an = 1:length(selAnimals_A)
%             zMIs_A(an,rr,cc) = nanmean(squeeze(pMs_A{rr,cc}.all_zMIs{selAnimals_A(an)}(nds(1),nds(2),:)));
%             a_zMIs_A{an,rr,cc} = squeeze(pMs_A{rr,cc}.all_zMIs{selAnimals_A(an)}(nds(1),nds(2),:));
%             gAllVals_A = [gAllVals_A;a_zMIs_A{an,rr,cc}];
%         end
%     end
% end
all_conds = unique(all_conds); all_rts = unique(all_rts);
% var_oi_A = squeeze(zMIs_A);
zMIs = squeeze(a_zMIs_C);
var_oi_C = squeeze(zMIs_C);
gAllVals = gAllVals_C;
paramMs = paramMs_C;
n = 0;
%%
runthis = 1;
if runthis
    data = zMIs;
    minBin = min(gAllVals);
    maxBin = 15;%max(gAllVals);
    incr = (maxBin-minBin)/100;
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[4 4 1.25 1],'color','w');
    hold on;
    [ha,hb,hca,sigR] = plotDistributions(data,'colors',colors,'maxY',90,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    hold on;
    legs = [];
    ylim([0 100]);xlim([-3 8]);
    if rsel==1
        xlim([-3 15]);
    end
    xlims = xlim; dx = xlims(2) - xlims(1); ylims = ylim; dy = ylims(2) - ylims(1);
    legs = {'Dist (C1-AD)','Time (C1-AT)'};%,'C3','C4'};
    legs{5} = [xlims(1)+dx/1.5 dx/30 ylims(1)+dy/1.75 dy/7];
    putLegend(ha,legs,'colors',colors,'sigR',{[],'',sigColor,6});
    axes(ha);
    h = xlabel({'Mutual Information','Z-Score (zMI)'});%changePosition(h,[0 -dy/3 0]);
    h = ylabel('Percentage');changePosition(h,[0.3 0 0]);
    set(gca,'FontSize',axes_font_size,'FontWeight','Bold');
    changePosition(ha,[0.07 0.01 -0.016 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution Of zMI 1 and 2',all_conds(1)),600);
    
    %%
    numCols = length(all_rts);
    data = var_oi_C;
    cmdTxt = sprintf('dataT_C = table(');
    for ii = 1:(size(data,2)-1)
        cmdTxt = sprintf('%sdata(:,%d),',cmdTxt,ii);
    end
    cmdTxt = sprintf('%sdata(:,size(data,2)));',cmdTxt);
    eval(cmdTxt);
    dataT = [dataT_C]
    dataT.Properties.VariableNames = varNames;
    colVar1 = [1 2];%[ones(1,numCols) 2*ones(1,numCols) 3*ones(1,numCols) 4*ones(1,numCols)];    colVar2 = [1:numCols 1:numCols 1:numCols 1:numCols];
    within = table(colVar1');
    within.Properties.VariableNames = {'Condition'};
    within.Condition = categorical(within.Condition);
    ra = repeatedMeasuresAnova(dataT,within);
    
%%
%     [h,p] = ttest(dataT{:,1},dataT{:,2});
    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
    xdata = [1 2]; maxY = 10;
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.2 1],'color','w');
    hold on;
    tcolors = {colors{1};colors{2};colors{3};colors{4};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'maxY',maxY,'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.3,'sigLinesStartYFactor',0.1);
    for ii = 1:length(hbs)
        set(hbs(ii),'facecolor',tcolors{ii},'edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.35],'ylim',[00 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = {'C1','C2','C3','C4'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.04 0.03 -0.15 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{'Average zMI',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bargraph',mfilename),600);
end
