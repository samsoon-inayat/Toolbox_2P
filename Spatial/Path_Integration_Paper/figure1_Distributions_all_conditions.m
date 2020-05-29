function figure1_Distributions_all_conditions

% dataAir = evalin('base','data');
% dataBelt = evalin('base','datab');
% dataAOn = evalin('base','dataAOn010');
% dataAOff = evalin('base','dataAOff010');
ei = evalin('base','ei10');
mData = evalin('base','mData');
colors = mData.colors;
sigColor = mData.sigColor;
axes_font_size = mData.axes_font_size;
% allCells = mData.allCells;
selAnimals = [1:9];
% selAnimals = 5:8;
% selAnimals = 1:8;
n = 0;
cN = 1;
%%
runthis = 0;
if runthis
%     cN = 1;
    planeNumbers = 'All';
    maxDistTime = [142 15];
    contextNumber = ones(1,4).*cN;
    
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
    for cnsi = 1:4
        for ss = 1:length(stimMarkers)
            distD = [];
            for jj = 1:length(selAnimals)
                [pcs cns areCells] = getParamValues('placeCells3',ei(selAnimals(jj)),planeNumbers,cnsi,'belt','dist',selCells,maxDistTime);
%                 [clus] = getParamValues('cluster3',ei(selAnimals(jj)),planeNumbers,cnsi,'air','dist',selCells,maxDistTime);
%                 [pcs cns areCells] = getParamValues('placeCells1',ei(selAnimals(jj)),planeNumbers,cnsi,'air','time',selCells,maxDistTime);
                pcs = logical(ones(size(pcs)));
%                 pcs = logical(clus(areCells));
                npcs(jj,cnsi) = sum(pcs)/length(pcs);
                [tempVals cns ACs] = getParamValues('info_metrics.ShannonMI_Zsh',ei(selAnimals(jj)),planeNumbers,cnsi,stimMarkers{ss},rasterTypes{ss},selCells,maxDistTime);
    %             [tempVals cns ACs] = getParamValues('gauss_fit_on_mean.rs',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),stimMarkers{ss},rasterTypes{ss},selCells,maxDistTime);
    %             [tempVals cns ACs] = getParamValues('data.rs',ei(selAnimals(jj)),planeNumbers,contextNumber,stimMarkers{ss},rasterTypes{ss},selCells,maxDistTime);
    %             [temppws cns ACs] = getParamValues('data.mean_trial_corr',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),stimMarkers{ss},rasterTypes{ss},selCells,maxDistTime);
    %             tempVals = temppws;
    %             tempVals = tempVals(logical(pcs));
                distD = [distD;tempVals(pcs)];
                zMI_mean(jj,ss) = nanmean(tempVals(pcs));
                data{jj,ss} = tempVals(pcs);
                datam(jj,ss) = nanmean(tempVals(pcs));
            end
            allVals{ss} = distD;
            gAllVals = [gAllVals;distD];
        end
        all_data = [all_data datam];
    end    
    ind = 1;
    for ii = 1:4
        for jj = 1:4
            varNames{ind} = sprintf('%s_Condition%d',xstimLabel{jj},ii);
            ind = ind + 1;
        end
    end

    data = all_data;
    cmdTxt = sprintf('dataT = table(');
    for ii = 1:(size(data,2)-1)
        cmdTxt = sprintf('%sdata(:,%d),',cmdTxt,ii);
    end
    cmdTxt = sprintf('%sdata(:,16));',cmdTxt);
    eval(cmdTxt);
    dataT.Properties.VariableNames = varNames;
    within = table([1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4]',[1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4]');
    cnRT = [[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4]',[1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4]'];
    within.Properties.VariableNames = {'Condition','Raster_Type'};
    within.Condition = categorical(within.Condition);
    within.Raster_Type = categorical(within.Raster_Type);

    % writetable(between,'Training_Data.xls');
    cmdTxt = sprintf('rm = fitrm(dataT,''');
    for ii = 1:(length(varNames)-1)
        cmdTxt = sprintf('%s%s,',cmdTxt,varNames{ii});
    end
    cmdTxt = sprintf('%s%s~1'');',cmdTxt,varNames{16});
    eval(cmdTxt);
    rm.WithinDesign = within;
    rm.WithinModel = 'Condition+Raster_Type';
    rtable = ranova(rm,'WithinModel',rm.WithinModel);
    mauchlytbl = mauchly(rm);
    % multcompare(rm,'Day','ComparisonType','bonferroni')
    mcTI = find_sig_mctbl(multcompare(rm,'Raster_Type','By','Condition','ComparisonType','bonferroni'),6);
    mcDays = find_sig_mctbl(multcompare(rm,'Condition','By','Raster_Type','ComparisonType','bonferroni'),6);
    [mVar semVar] = findMeanAndStandardError(all_data);
    
    combs = nchoosek(1:size(all_data,2),2); p = ones(size(combs,1),1); h = logical(zeros(size(combs,1),1));
    for rr = 1:size(mcTI,1)
        thisRow = mcTI(rr,:);
        conditionN =  thisRow{1,1}; Rtype1 = thisRow{1,2}; Rtype2 = thisRow{1,3};
        Num1 = find(ismember(within{:,:},[conditionN Rtype1],'rows'));
        Num2 = find(ismember(within{:,:},[conditionN Rtype2],'rows'));
        row = [Num1 Num2]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{1,6}; h(ii) = 1;
    end
    for rr = 1:size(mcDays,1)
        thisRow = mcDays(rr,:);
        Rtype =  thisRow{1,1}; Condition1 = thisRow{1,2}; Condition2 = thisRow{1,3};
        Num1 = find(ismember(within{:,:},[Condition1 Rtype],'rows'));
        Num2 = find(ismember(within{:,:},[Condition2 Rtype],'rows'));
        row = [Num1 Num2]; ii = ismember(combs,row,'rows'); p(ii) = mcDays{1,6}; h(ii) = 1;
    end
%%
    n = 0;
    xdata = [1:1.5:(10*size(all_data,2))]; xdata = xdata(1:size(all_data,2)); maxY = 30;

    hf = figure(15);clf;set(gcf,'Units','Inches');set(gcf,'Position',[3 3 4 1.5],'color','w');
    hold on;
    ind = 1;
    for ii = 1:4
        for jj = 1:4
            tcolors{ind} = colors{ii};
            ind = ind + 1;
        end
    end
    hbs = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'maxY',maxY,'ySpacing',1.25,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.7,'sigLinesStartYFactor',0);
    set(gca,'xlim',[0.25 max(xdata)+.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata; 
    xticklabels = {'AirD','AirT','BeltD','AirIT'};xticklabels = repmat(xticklabels,1,4);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0.1 0.02 -0.03 -0.011])
    put_axes_labels(gca,{[],[0 0 0]},{{'Mutual Information','(z-score)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Mean zMI all Conditions'),600);
    
    npcs = 100*npcs;
    data = npcs;
    cmdTxt = sprintf('dataT = table(');
    for ii = 1:(size(data,2)-1)
        cmdTxt = sprintf('%sdata(:,%d),',cmdTxt,ii);
    end
    cmdTxt = sprintf('%sdata(:,%d));',cmdTxt,size(data,2));
    eval(cmdTxt);
    varNames = {'C1','C2','C3','C4'};
    dataT.Properties.VariableNames = varNames;
    within = table([1 2 3 4]');
    cnRT = [[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4]',[1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4]'];
    within.Properties.VariableNames = {'Condition'};
    within.Condition = categorical(within.Condition);

    % writetable(between,'Training_Data.xls');
    cmdTxt = sprintf('rm = fitrm(dataT,''');
    for ii = 1:(length(varNames)-1)
        cmdTxt = sprintf('%s%s,',cmdTxt,varNames{ii});
    end
    cmdTxt = sprintf('%s%s~1'');',cmdTxt,varNames{length(varNames)});
    eval(cmdTxt);
    rm.WithinDesign = within;
    rm.WithinModel = 'Condition';
    rtable = ranova(rm,'WithinModel',rm.WithinModel);
    mauchlytbl = mauchly(rm);
    % multcompare(rm,'Day','ComparisonType','bonferroni')
    mcDays = find_sig_mctbl(multcompare(rm,'Condition','ComparisonType','bonferroni'))
    
    [mVar semVar] = findMeanAndStandardError(npcs);
    all_data = npcs;
    combs = nchoosek(1:size(all_data,2),2); p = ones(size(combs,1),1); h = logical(zeros(size(combs,1),1));
    for rr = 1:size(mcDays,1)
        thisRow = mcTI(rr,:);
        conditionN =  thisRow{1,1}; Rtype1 = thisRow{1,2}; Rtype2 = thisRow{1,3};
        Num1 = find(ismember(within{:,:},[conditionN Rtype1],'rows'));
        Num2 = find(ismember(within{:,:},[conditionN Rtype2],'rows'));
        row = [Num1 Num2]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{1,6}; h(ii) = 1;
    end
%%
    n = 0;
    xdata = [1:1.5:(10*size(all_data,2))]; xdata = xdata(1:size(all_data,2)); maxY = 50;

    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1 1],'color','w');
    hold on;
    tcolors = mData.colors;
    hbs = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'maxY',maxY,'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0);
    set(gca,'xlim',[0.25 max(xdata)+.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata; 
    xticklabels = {'C1','C2','C3','C4'};xticklabels = repmat(xticklabels,1,4);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0.1 0.02 -0.03 -0.011])
    put_axes_labels(gca,{[],[0 0 0]},{{'Perc Cells'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Percent Cells zMI threshold'),600);
    
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
