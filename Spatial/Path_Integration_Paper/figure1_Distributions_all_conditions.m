function figure1_Distributions_all_conditions


protocol = '10';
% protocol = '15';
ei = evalin('base',sprintf('ei%s',protocol));
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ET = evalin('base',sprintf('ET%s',protocol));
selAnimals = eval(sprintf('mData.selAnimals%s',protocol));

% in the following variable all the measurements are in the matrices form
% for each variable colums indicate raster and stim marker types specified 
% the rows indicate condition numbers.
paramMs = parameter_matrices('get',protocol);
% after getting all matrics, we can apply selection criteria to select a
% subgroup of cells
% here is the selection criteria in make_selC_structure function
% cellsOrNot = NaN; planeNumber = NaN; zMI_Th = 2; fwids = [0 140]; fcens = [0 140]; rs_th = NaN;
cellsOrNot = NaN; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN;
conditionsAndRasterTypes = [11 12 13 14 15 21 22 23 24 25 31 32 33 34 35 41 42 43 44 45]';
% conditionsAndRasterTypes = [11 13 21 23 31 33 41 43]';
% conditionsAndRasterTypes = [11 21 31 41]';
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th);
[cpMs,pMs] = parameter_matrices('select',protocol,{paramMs,selC});
perc_cells = parameter_matrices('print',protocol,{cpMs,pMs,ET,selAnimals});

for rr = 1:size(pMs,1)
    for cc = 1:size(pMs,2)
        tcond = conditionsAndRasterTypes(rr,cc);
        nds = dec2base(tcond,10) - '0';
        varNames{rr,cc} = sprintf('C%dR%d',nds(1),nds(2));
        for an = 1:length(selAnimals)
            zMIs(an,rr,cc) = nanmean(squeeze(pMs{rr,cc}.all_zMIs{selAnimals(an)}(nds(1),nds(2),:)));
        end
    end
end

all_conds = []; all_rts = [];
for rr = 1:size(pMs,1)
    for cc = 1:size(pMs,2)
        tcond = conditionsAndRasterTypes(rr,cc);
        nds = dec2base(tcond,10) - '0';
        varNames{rr,cc} = sprintf('C%dR%d',nds(1),nds(2));
        xticklabels{cc,rr} = sprintf('%s-%s',paramMs.stimMarkers{nds(2)},paramMs.rasterTypes{nds(2)}(1));
        all_conds = [all_conds nds(1)]; all_rts = [all_rts nds(2)];
    end
end
all_conds = unique(all_conds); all_rts = unique(all_rts);
n = 0;
cN = 1;
%%
runthis = 1;
if runthis
%     cN = 1;
    numCols = length(all_rts); 
    data = zMIs;
%     cmdTxt = sprintf('dataT = table(');
%     for ii = 1:(size(data,2)-1)
%         cmdTxt = sprintf('%sdata(:,%d),',cmdTxt,ii);
%     end
%     cmdTxt = sprintf('%sdata(:,size(data,2)));',cmdTxt);
%     eval(cmdTxt);
%     dataT.Properties.VariableNames = varNames;
    colVar1 = [ones(1,numCols) 2*ones(1,numCols) 3*ones(1,numCols) 4*ones(1,numCols)];    colVar2 = [1:numCols 1:numCols 1:numCols 1:numCols];
    within = table(colVar1',colVar2');
    within.Properties.VariableNames = {'Condition','Raster_Type'};
    within.Condition = categorical(within.Condition);
    within.Raster_Type = categorical(within.Raster_Type);
    ra = repeatedMeasuresAnova(data,varNames,within);
    rm = ra.rm;
    if length(all_rts) > 1
        mcTI = find_sig_mctbl(multcompare(rm,'Raster_Type','By','Condition','ComparisonType','bonferroni'),6);
        mcDays = find_sig_mctbl(multcompare(rm,'Condition','By','Raster_Type','ComparisonType','bonferroni'),6);
    else
        mcTI = [];
        mcDays = find_sig_mctbl(multcompare(rm,'Condition','ComparisonType','bonferroni'),5);
    end
    [mVar semVar] = findMeanAndStandardError(data);
    [combs,h,p] = populate_multcomp_h_p(data,within,mcTI,[]);
    
    xdata = [1:1.5:(10*size(data,2))]; xdata = xdata(1:size(data,2)); 
    ind = 1;
    for ii = 1:length(all_conds)
        for jj = 1:length(all_rts)
            tcolors{ind} = colors{ii};
            ind = ind + 1;
        end
    end
    hf = figure(15);clf;set(gcf,'Units','Inches');set(gcf,'Position',[3 3 6 1.5],'color','w');
    hold on;
    [hbs,maxYr] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 max(xdata)+.75],'ylim',[0 maxYr],'FontSize',7,'FontWeight','Bold','TickDir','out');
    xticks = xdata; 
%     xticklabels = repmat(xticklabels,length(all_rts),length(all_conds));
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0 0.02 0.03 -0.011])
    put_axes_labels(gca,{[],[0 0 0]},{{'Mutual Information','(z-score)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Mean zMI all Conditions'),600);
return;
end


runthis = 1;
if runthis
%     cN = 1;
    planeNumbers = 'All';
    maxDistTime = [142 15];
   
    stimMarkers = {'air','air','belt','airI'};
    stimLabel = {'Air','Air','Belt','InterTrial'};
    xstimLabel = {'AirD','AirT','BeltD','AirIT'};
    rasterTypes = {'dist','time','dist','time'};
    rasterLabel = {'Distance','Time','Distance','Time'};
%     stimMarkers = {'air','air','air','air'};
%     rasterTypes = {'dist','dist','dist','dist'};
    selCells = 'areCells';
    allVals = []; gAllVals = [];
    all_data = [];
    sel_ss = 4;
    for cnsi = 1:4
        for ss = 1:length(stimMarkers)
            distD = [];
            for jj = 1:length(selAnimals)
                [pcs cns areCells] = getParamValues('placeCells3',ei(selAnimals(jj)),planeNumbers,cnsi,stimMarkers{sel_ss},rasterTypes{sel_ss},selCells,maxDistTime);
%                 pcs = logical(ones(1,length(pcs)));
                npcs(jj,cnsi) = sum(pcs)/length(pcs);
                [tempVals cns ACs] = getParamValues('info_metrics.ShannonMI_Zsh',ei(selAnimals(jj)),planeNumbers,cnsi,stimMarkers{ss},rasterTypes{ss},selCells,maxDistTime);
                distD = [distD;tempVals(pcs)];
                zMI_mean(jj,ss) = nanmean(tempVals(pcs));
                data{jj,ss} = tempVals(pcs);
                datam(jj,ss) = nanmean(tempVals(pcs));
                all_data_zMIs{jj,cnsi,ss} = [tempVals];
                all_data_pcs{jj,cnsi,ss} = [pcs];
            end
            allVals{ss} = distD;
            gAllVals = [gAllVals;distD];
        end
        all_data = [all_data datam];
    end
    n=0;
    %%
    
    ss = sel_ss;
    azMI = []; apcs = [];
    for cnsi = 1:4
        zMI_C =[]; pcs_C = [];
        for an = 1:size(all_data_pcs,1)
            this = all_data_zMIs{an,cnsi,ss};
            zMI_C = [zMI_C;this(:,1)];
            this = all_data_pcs{an,cnsi,ss};
            pcs_C = [pcs_C;this(:,1)];
        end
        azMI(:,cnsi) = zMI_C;
        apcs(:,cnsi) = pcs_C;
    end
    combs = [1 2;2 3;3 4];
    n=0;
    apcs = logical(apcs);
    the_pcs = apcs(:,1) | apcs(:,2) | apcs(:,3) | apcs(:,4); 
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[8 6 6 2.5],'color','w');hold on;
    ya = []; xa = []; gr = [];
    for ii = 1:size(combs,1)
        xvals = azMI(the_pcs,combs(ii,1));
        yvals = azMI(the_pcs,combs(ii,2));
        indsx = isnan(xvals); indsy = isnan(yvals);
        xvals(indsx | indsy) = []; yvals(indsx | indsy) = [];
        scatter(xvals,yvals,20,'.','MarkerEdgeColor',colors{ii},'MarkerFaceColor','none');
        ft = fittype('poly1');
        [ftc,gof,output] = fit(xvals,yvals,ft);    rsq(ii) = gof.rsquare;
        co = coeffvalues(ftc);
        yhvals = ft(co(1),co(2),xvals);
        hold on;plot(xvals,yhvals,'color',colors{ii},'LineWidth',1)
        text(20,15-3*ii,sprintf('Rsq = %.3f',rsq(ii)),'color',colors{ii});
        xa = [xa;xvals]; ya = [ya;yvals]; gr = [gr;ones(size(xvals))*ii];
    end
    sigR = do_ancova(xa,ya,gr);
    ylim([0 30]);
    set(gca,'TickDir','out','FontSize',mData.axes_font_size,'FontWeight','Bold');
    hx = xlabel('Mutual Information (Z score)');%changePosition(hx,[0 0.0 0]);
    hy = ylabel('Mutual Information (Z score)');%changePosition(hy,[2.2 -3.0 0]);
    for ii = 1:size(combs,1)
        legs{ii} = sprintf('C %d-%d',combs(ii,1),combs(ii,2));%{'Context 1','Context 2','Context 3','Context 4',[30 3 35 3]};
    end
    legs{ii+1} = [5 1 25 3];putLegend(gca,legs,'colors',colors,'sigR',{sigR,'ancovas',sigColor,9}); text(3,28,'Slope','FontWeight','normal');
    legs{ii+1} = [20 1 25 3];putLegend(gca,legs,'colors',colors,'sigR',{sigR,'ancovai',sigColor,9}); text(20,28,'Intercept','FontWeight','normal');
    changePosition(gca,[-0.03 0.09 0.09 -0.06]);
    
    save_pdf(hf,mData.pdf_folder,sprintf('scatter zMI'),600);
    %%
return;
end
