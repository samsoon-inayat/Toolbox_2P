function figure1_Distributions_Conds_Dist_Time

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;


cellsOrNot = 1; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN;
conditionsAndRasterTypes = [11 12 21 22 31 32 41 42];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN);

out = read_data_from_base_workspace(selC)

pMs_C = out.pMs{1}; pMs_A = out.pMs{2};
paramMs_C = out.paramMs{1}; paramMs_A = out.paramMs{2};
selAnimals_C = out.selAnimals{1}; selAnimals_A = out.selAnimals{2};

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
            temp = squeeze(pMs_C{rr,cc}.all_zMIs{selAnimals_C(an)}(nds(1),nds(2),:));
            if sum(temp>100) > 0
                temp(temp>100) = NaN;
                large_values = large_values + 1;
                disp('there were large values');
                n = 0;
            end
            zMIs_C(an,rr,cc) = nanmean(squeeze(pMs_C{rr,cc}.all_zMIs{selAnimals_C(an)}(nds(1),nds(2),:)));
            a_zMIs_C{an,rr,cc} = squeeze(pMs_C{rr,cc}.all_zMIs{selAnimals_C(an)}(nds(1),nds(2),:));
            gAllVals_C = [gAllVals_C;a_zMIs_C{an,rr,cc}];
        end
    end
end
large_values = 1;
for rr = 1:size(pMs_A,1)
    for cc = 1:size(pMs_A,2)
        tcond = conditionsAndRasterTypes(rr,cc);
        nds = dec2base(tcond,10) - '0';
        for an = 1:length(selAnimals_A)
            temp = squeeze(pMs_A{rr,cc}.all_zMIs{selAnimals_A(an)}(nds(1),nds(2),:));
            if sum(temp>100) > 0
                temp(temp>100) = NaN;
                large_values = large_values + 1;
                disp('there were large values');
                n = 0;
            end
            zMIs_A(an,rr,cc) = nanmean(temp);%nanmean(squeeze(pMs_A{rr,cc}.all_zMIs{selAnimals_A(an)}(nds(1),nds(2),:)));
            a_zMIs_A{an,rr,cc} = squeeze(pMs_A{rr,cc}.all_zMIs{selAnimals_A(an)}(nds(1),nds(2),:));
            gAllVals_A = [gAllVals_A;a_zMIs_A{an,rr,cc}];
        end
    end
end
all_conds = unique(all_conds); all_rts = unique(all_rts);
% var_oi_A = squeeze(zMIs_A);
zMIs = squeeze(a_zMIs_C);
var_oi_C = squeeze(zMIs_C);
gAllVals = gAllVals_C;
paramMs = paramMs_C;

zMIs = squeeze(a_zMIs_A);
var_oi_C = squeeze(zMIs_A);
gAllVals = gAllVals_A;
paramMs = paramMs_A;
% large_values
n = 0;
%%
runthis = 1;
if runthis
    data = zMIs;
    minBin = min(gAllVals);
    maxBin = 15;%max(gAllVals);
    incr = (maxBin-minBin)/100;
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[4 4 1.5 1],'color','w');
    hold on;
    [ha,hb,hca,sigR] = plotDistributions(data,'colors',colors,'maxY',90,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    hold on;
    legs = [];
    ylim([0 100]);xlim([-3 7]);
    xlims = xlim; dx = xlims(2) - xlims(1); ylims = ylim; dy = ylims(2) - ylims(1);
    legs = {'C1','C2'};%,'C3','C4'};
    legs{5} = [xlims(1)+dx/1.5 dx/30 ylims(1)+dy/1.75 dy/7];
    putLegend(ha,legs,'colors',colors,'sigR',{[],'',sigColor,6});
    axes(ha);
    h = xlabel({'Mutual Information','Z-Score (zMI)'});%changePosition(h,[0 -dy/3 0]);
    h = ylabel('Percentage');changePosition(h,[0.3 0 0]);
    set(gca,'FontSize',axes_font_size,'FontWeight','Bold');
    changePosition(ha,[0.07 0.01 -0.016 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution Of zMI %d',all_conds(1)),600);
    
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
    colVar1 = [ones(1,numCols) 2*ones(1,numCols) 3*ones(1,numCols) 4*ones(1,numCols)];    colVar2 = [1:numCols 1:numCols 1:numCols 1:numCols];
%     colVar1 = [1 2 3 4];
    within = table(colVar1');
    within = table(colVar1',colVar2');
%     within.Properties.VariableNames = {'Condition'};
    within.Properties.VariableNames = {'Condition','Raster'};
    within.Condition = categorical(within.Condition);
    within.Raster = categorical(within.Raster);
    ra = repeatedMeasuresAnova(dataT,within);
    
%%
%     [h,p] = ttest(dataT{:,1},dataT{:,2});
    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
    xdata = [1 2 3 4 5 6 7 8]; maxY = 10;
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 2 1],'color','w');
    hold on;
    tcolors ={colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.6,'sigLinesStartYFactor',0.1);
    for ii = 1:length(hbs)
        set(hbs(ii),'facecolor',tcolors{ii},'edgecolor',tcolors{ii});
    end
    for ii = 2:2:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.35],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = {'C1-AD','C1-AT','C2-AD','C2-AT','C3-AD','C3-AT','C4-AD','C4-AT'};
    xtickangle(30)
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     changePosition(gca,[0.04 0.03 -0.15 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{'zMI',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bargraph',mfilename),600);
end
