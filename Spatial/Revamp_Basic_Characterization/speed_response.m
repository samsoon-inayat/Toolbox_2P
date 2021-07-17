function speed_response

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei = evalin('base','ei'); 
var_names = {'linear','sigmoid','gauss'};
for ii = 1:length(ei)
    tei = ei{ii};
    psp = [];
    psp1 = tei.plane{1}.speed_response;
    if length(tei.plane) == 2
        psp2 = tei.plane{2}.speed_response;
        psp.corr = [psp1.corr;psp2.corr];
        psp.FR_vs_speed = [psp1.FR_vs_speed;psp2.FR_vs_speed];
         for vv = 1:length(var_names)
            cmdTxt = sprintf('psp.fits.%s.fitted = [psp1.fits.%s.fitted;psp2.fits.%s.fitted];',var_names{vv},var_names{vv},var_names{vv});
            eval(cmdTxt);
            cmdTxt = sprintf('psp.fits.%s.coeffsrs = [psp1.fits.%s.coeffsrs;psp2.fits.%s.coeffsrs];',var_names{vv},var_names{vv},var_names{vv});
            eval(cmdTxt);
         end
        psp.bin_centers = psp1.bin_centers;
        speedRs{ii} = psp;
    else
        speedRs{ii} = psp1;
    end
end
n = 0;
%% visualize the data
if 1
    an = 4;
    d.bcs = speedRs{an}.bin_centers;
    d.FR = speedRs{an}.FR_vs_speed;
    fitg = speedRs{an}.fits.gauss; fits = speedRs{an}.fits.sigmoid; fitl = speedRs{an}.fits.linear;
    d.fFRl = fitl.fitted; d.fFRs = fits.fitted; d.fFRg = fitg.fitted;
    d.cl = fitl.coeffsrs(:,3); d.cs = fits.coeffsrs(:,3); d.cg = fitg.coeffsrs(:,3);
    [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,d.bcs(2)-d.bcs(1));
    inds = centers < 1 | centers > 39 | rs < 0.25 | PWs < 10;% | PWs > 20 | PWs < 10;
    inds = ~inds;
    100*sum(inds)/length(inds)
    d.FR = d.FR(inds,:); d.fFRl = d.fFRl(inds,:); d.fFRs = d.fFRs(inds,:); d.fFRg = d.fFRg(inds,:);
    d.cl = d.cl(inds); d.cs = d.cs(inds); d.cg = d.cg(inds);
    generic_multi_plot(1000,[3,4,size(d.FR,1)],'plotSpeedTuning',d)
    return;
end
%% see the distribution of rs for linear, sigmoid, and gaussian fitting.
if 0
   tcolors = {'k','m','b'};
   distD = [];
   ind = [3 3 4];
   for ii = 1:length(ei)
       ii
       for vv = 1:length(var_names)
           cmdTxt = sprintf('distD{ii,vv} = speedRs{ii}.fits.%s.coeffsrs(:,ind(vv));',var_names{vv});
           eval(cmdTxt);
           infinds = find(distD{ii,vv}==-inf | distD{ii,vv}==inf);
           distD{ii,vv}(infinds) = NaN;
       end
   end
   [distDo,allVals,allValsG] = plotDistributions(distD);
   minBin = -2;
   maxBin = max(allVals);
   incr = 0.01;
   hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.85 1],'color','w');
   hold on;
   [ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
   set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
   legs = {'Linear','Sigmoid','Gaussian',[-1.75 0.3 85 5]};
   putLegend(gca,legs,tcolors);
   ylim([0 110]);
   changePosition(gca,[0.09 0.13 -0.05 -0.13]);
    put_axes_labels(gca,{'RS',[0 0 0]},{{'Percentage','of Neurons'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_speed_rs'),600);
end

%% mean rs values
if 0
    meanDistD = arrayfun(@(x)nanmean(x{1}),distD);
    within = make_within_table({'Type'},3);
    dataT = make_between_table({meanDistD},{'Rs_L','Rs_S','Rs_G'});
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w'); hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.25,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.05);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[-0.75 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = [xdata(1:end)]; xticklabels = {'Linear','Sigmoid','Gaussian'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30)
    changePosition(gca,[0.2 0.03 -0.1 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Ordinary','R-Square'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('mean_Rsquared'),600);
end

%% distribution preferred speed
if 1
    for ii = 1:length(ei)
        bcs = speedRs{ii}.bin_centers;
        fitg = speedRs{ii}.fits.gauss;
        [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,bcs(2)-bcs(1));
        inds = centers < 1 | centers > 39; centers(inds) = [];
       distD{ii,1} = centers;
    end
    [distDo,allVals,allValsG] = plotDistributions(distD);
   minBin = min(allVals);
   maxBin = max(allVals);
   incr = 5;
   hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.85 1],'color','w');
   hold on;
   [ha,hb,hca] = plotAverageDistributions(distD,'colors',{'k'},'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
   set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
   legs = {'Linear','Sigmoid','Gaussian',[-1.75 0.3 85 5]};
   putLegend(gca,legs,tcolors);
   ylim([0 110]);
   changePosition(gca,[0.09 0.13 -0.05 -0.13]);
    put_axes_labels(gca,{'RS',[0 0 0]},{{'Percentage','of Neurons'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_speed_preferred'),600);
    return;
end

%% Percentage of Responsive Cells
if 1
    within = make_within_table({'Cond'},3);
    dataT = make_between_table({resp_fractionC*100},{'C21','C22'});
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w'); hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',colors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.05);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 30],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = [xdata(1:end)]; xticklabels = {'C2','C2'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30)
    any_mean = mean(100*resp_OR_fractionC);    any_sem = std(100*resp_OR_fractionC)/sqrt(5);
    pmchar=char(177); any_text = sprintf('%.0f%c%.0f%%',any_mean,pmchar,any_sem); text(0.75,25,any_text,'FontSize',6);
    changePosition(gca,[0.2 0.03 -0.5 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Air Responsive','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('percentage_air_responsive'),600);
end