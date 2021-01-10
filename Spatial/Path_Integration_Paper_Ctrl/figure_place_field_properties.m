function figure_place_field_properties(fn,allRs,ccs)

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
% cellsOrNot = 1; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN;
cellsOrNot = 1; planeNumber = NaN; zMI_Th = 1.96; fwids = [1 150]; fcens = [0 150]; rs_th = 0.3;
conditionsAndRasterTypes = [11;21;31;41];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN);
out = read_data_from_base_workspace(selC)

ei_C = out.eis{1}; ei_A = out.eis{2};
pMs_C = out.pMs{1}; pMs_A = out.pMs{2};
paramMs_C = out.paramMs{1}; paramMs_A = out.paramMs{2};
selAnimals_C = out.selAnimals{1}; selAnimals_A = out.selAnimals{2};

param_names = {'all_zMIs','all_fFR','all_fwidths'};
ylabels = {'zMI','Firing Rate (Hz)','Field Width (cm)'};

sp = 1;

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
            cmdTxt = sprintf('val_mat = pMs_C{rr,cc}.%s{selAnimals_C(an)}(nds(1),nds(2),:);',param_names{sp});
            eval(cmdTxt);
            zMIs_C(an,rr,cc) = nanmean(squeeze(val_mat));
            a_zMIs_C{an,rr,cc} = squeeze(val_mat);
            gAllVals_C = [gAllVals_C;a_zMIs_C{an,rr,cc}];
        end
    end
end
for rr = 1:size(pMs_A,1)
    for cc = 1:size(pMs_A,2)
        tcond = conditionsAndRasterTypes(rr,cc);
        nds = dec2base(tcond,10) - '0';
        for an = 1:length(selAnimals_A)
            cmdTxt = sprintf('val_mat = pMs_A{rr,cc}.%s{selAnimals_A(an)}(nds(1),nds(2),:);',param_names{sp});
            eval(cmdTxt);
            zMIs_A(an,rr,cc) = nanmean(squeeze(val_mat));
            a_zMIs_A{an,rr,cc} = squeeze(val_mat);
            gAllVals_A = [gAllVals_A;a_zMIs_A{an,rr,cc}];
        end
    end
end
all_conds = unique(all_conds); all_rts = unique(all_rts);
var_oi_A = squeeze(zMIs_A);
var_oi_C = squeeze(zMIs_C);
n = 0;
%%
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
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 2 1],'color','w');
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
    maxY = maxY + 5 + 2;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+1],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = {'C1','C2','C3','C4'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     xtickangle(30)
    changePosition(gca,[0.03 0.03 0.04 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{ylabels{sp},[0 0 0]});
    rectangle(gca,'Position',[0.75 maxY-3 1 1.25],'edgecolor','k','facecolor','k');     text(1.85,maxY-3+1,'CRTG','FontSize',6);
    rectangle(gca,'Position',[6 maxY-3 1 1.25],'edgecolor','k');     text(7.2,maxY-3+1,'CPTG','FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('props_bargraph_%s',ylabels{sp}(1:3)),600);
return;
end



%% scatter place field widths vs centers
runThis = 0;
if runThis
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 6 2.5],'color','w');hold on;
ya = []; xa = []; gr = [];
for ii = 1:2%length(all_pws_c)
    pwidths = all_pws_c{ii};
    pcenters = all_cs_c{ii};
    scatter(pcenters,pwidths,20,'.','MarkerEdgeColor',colors{ii},'MarkerFaceColor','none');
    ft = fittype('poly1');
    [ftc,gof,output] = fit(pcenters,pwidths,ft);    rsq(ii) = gof.rsquare;
    co = coeffvalues(ftc);
    pwf = ft(co(1),co(2),pcenters);
    hold on;plot(pcenters,pwf,'color',colors{ii},'LineWidth',1)
    xa = [xa;pcenters]; ya = [ya;pwidths]; gr = [gr;ones(size(pcenters))*ii];
end
sigR = do_ancova(xa,ya,gr);
ylim([0 100]);
set(gca,'TickDir','out','FontSize',mData.axes_font_size,'FontWeight','Bold');
hx = xlabel('Place Field Center (cm)');changePosition(hx,[0 0.0 0]);
hy = ylabel('Place Field Width (cm)');changePosition(hy,[2.2 -3.0 0]);
for ii = 1:length(all_pws_c)
    legs{ii} = sprintf('C %d',ii);%{'Context 1','Context 2','Context 3','Context 4',[30 3 35 3]};
end
legs{ii+1} = [35 5 93 7];putLegend(gca,legs,'colors',colors,'sigR',{sigR,'ancovas',sigColor,9}); text(30,100,'Slope','FontWeight','Bold');
legs{ii+1} = [100 5 93 7];putLegend(gca,legs,'colors',colors,'sigR',{sigR,'ancovai',sigColor,9}); text(95,100,'Intercept','FontWeight','Bold');
changePosition(gca,[-0.03 0.09 0.09 -0.06]);
save_pdf(hf,mData.pdf_folder,'figure_scatter_centers_widths.pdf',600);
return;
end

%% scatter zMI vs centers
runThis = 0;
if runThis
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 6 2.5],'color','w');hold on;
ya = []; xa = []; gr = [];
for ii = 1:length(pw)
    SIs = SI{ii}(find(pcs5{ii}))';
    pcenters = pc{ii}(find(pcs5{ii}))';
    scatter(pcenters,SIs,20,'.','MarkerEdgeColor',colors{ii},'MarkerFaceColor','none');
    ft = fittype('poly1');
    [ftc,gof,output] = fit(pcenters,SIs,ft);  rsq(ii) = gof.rsquare;
    co = coeffvalues(ftc);
    pwf = ft(co(1),co(2),pcenters);
    hold on;plot(pcenters,pwf,'color',colors{ii},'LineWidth',1)
    xa = [xa;pcenters]; ya = [ya;SIs]; gr = [gr;ones(size(pcenters))*ii];
end
sigR = do_ancova(xa,ya,gr);
ylim([0 35]);
set(gca,'TickDir','out','FontSize',mData.axes_font_size,'FontWeight','Bold');
hx = xlabel('Place Field Center (cm)');changePosition(hx,[0 0.0 0]);
hy = ylabel('Mutual Info. (Z Score)');changePosition(hy,[0 0.0 0]);
for ii = 1:length(pw)
    legs{ii} = sprintf('C %d',ii);%{'Context 1','Context 2','Context 3','Context 4',[30 3 35 3]};
end
legs{ii+1} = [130 5 33 2];putLegend(gca,legs,'colors',colors,'sigR',{sigR,'ancovai',sigColor,9}); text(125,36,'Intercept','FontWeight','Bold'); 
legs{end} = [70 5 33 2];putLegend(gca,legs,'colors',colors,'sigR',{sigR,'ancovas',sigColor,9}); text(65,36,'Slope','FontWeight','Bold');

changePosition(gca,[-0.03 0.09 0.09 -0.06]);
save_pdf(hf,mData.pdf_folder,'figure_scatter_centers_zMI.pdf',600);
return;
end

