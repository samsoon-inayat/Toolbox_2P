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
    an = 2;
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
%%
if 1
    an = 4;
    d.bcs = speedRs{an}.bin_centers;
    d.FR = speedRs{an}.FR_vs_speed;
    fitg = speedRs{an}.fits.gauss; fits = speedRs{an}.fits.sigmoid; fitl = speedRs{an}.fits.linear;
    d.fFRl = fitl.fitted; d.fFRs = fits.fitted; d.fFRg = fitg.fitted;
    d.cl = fitl.coeffsrs(:,3); d.cs = fits.coeffsrs(:,3); d.cg = fitg.coeffsrs(:,4);
    [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,d.bcs(2)-d.bcs(1));
    inds = centers < 1 | centers > 39 | rs < 0.25 | PWs < 10;% | PWs > 20 | PWs < 10;
    inds = ~inds;
    100*sum(inds)/length(inds)
    d.FR = d.FR(inds,:); d.fFRl = d.fFRl(inds,:); d.fFRs = d.fFRs(inds,:); d.fFRg = d.fFRg(inds,:);
    d.cl = d.cl(inds); d.cs = d.cs(inds); d.cg = d.cg(inds);
    cell_inds = [24 48 54];
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 3],...
        'spaceRowsCols',[0.15 0.06],'rightUpShifts',[0.14 0.26],'widthHeightAdjustment',...
        [-100 -415]);
    set(gcf,'color','w'); set(gcf,'Position',[10 4 3 1]);
    for iii = 1:3
        ii = cell_inds(iii);
        rs = [d.cl(ii) d.cs(ii) d.cg(ii)];
        [~,mind] = max(rs);
        axes(ff.h_axes(1,iii));
        plot(d.bcs,d.FR(ii,:),'c');hold on;
        pl(1) = plot(d.bcs,d.fFRl(ii,:),'k','linewidth',1);
        pl(2) = plot(d.bcs,d.fFRs(ii,:),'m','linewidth',0.5);
        pl(3) = plot(d.bcs,d.fFRg(ii,:),'b','linewidth',1);
%         set(pl(mind),'linewidth',2);
%         title(sprintf('Cell %d',ii));
        if iii > 1
            set(gca,'YTick',[]);
        else
            ylabel({'Firing','Rate (AU)'});
        end
        if iii == 2
            xlabel('Speed (cm/sec)')
        end
        if iii == 1
            legs = {'Linear','Sigmoid','Gaussian',[25 2 0.4 0.01]};
           putLegend(gca,legs,{'k','m','b'});
        end
        box off;
        set(gca,'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('sample_speed_tuning'),600);
end
%% see the distribution of rs for linear, sigmoid, and gaussian fitting.
if 1
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
   hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
   hold on;
   [ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
   set(gca,'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
   legs = {'Linear','Sigmoid','Gaussian',[-1.75 0.15 100 5]};
   putLegend(gca,legs,tcolors);
   ylim([0 110]);
   changePosition(gca,[0.2 0.13 -0.15 -0.13]);
    put_axes_labels(gca,{'R-Squared',[0 0 0]},{{'Percentage','of Neurons'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_speed_rs'),600);
end

%% mean rs values
if 1
    meanDistD = arrayfun(@(x)nanmean(x{1}),distD);
    within = make_within_table({'Type'},3);
    dataT = make_between_table({meanDistD},{'Rs_L','Rs_S','Rs_G'});
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w'); hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.25,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[-0.75 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = [xdata(1:end)]; xticklabels = {'Linear','Sigmoid','Gaussian'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30)
    changePosition(gca,[0.11 0.03 -0.3 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'R-Squared'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('mean_Rsquared'),600);
end

%% distribution preferred speed
if 1
    distD = [];
    for ii = 1:length(ei)
        bcs = speedRs{ii}.bin_centers;
        fitg = speedRs{ii}.fits.gauss;
        [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,bcs(2)-bcs(1));
        inds = centers < 1 | centers > 39 | PWs < 1 | PWs > 40; centers(inds) = [];
       distD{ii,1} = centers;
       meanPWs(ii) = (mean(centers));
    end
    [distDo,allVals] = getAveragesAndAllValues(distD);
   minBin = min(allVals);
   maxBin = max(allVals);
   incr = 2;
   hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
   hold on;
   [ha,hb,hca] = plotAverageDistributions(distD,'colors',{'k'},'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
   set(gca,'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
%    legs = {'Linear','Sigmoid','Gaussian',[-1.75 0.3 85 5]};
%    putLegend(gca,legs,tcolors);
   ylim([0 110]);
   pmchar=char(177); any_text = sprintf('%.0f%c%.0f cm/sec',mean(meanPWs),pmchar,std(meanPWs)/sqrt(5)); 
   text(10,15,any_text,'FontSize',6);
   changePosition(gca,[0.12 0.13 -0.3 -0.13]);
    put_axes_labels(gca,{'Pref. Speed (cm/sec)',[-4 0 0]},{{'Neurons(%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_speed_preferred'),600);
    return;
end

%% distribution tuning width
if 1
    distD = [];
    for ii = 1:length(ei)
        bcs = speedRs{ii}.bin_centers;
        fitg = speedRs{ii}.fits.gauss;
        [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,bcs(2)-bcs(1));
        inds = centers < 1 | centers > 39 | PWs < 1 | PWs > 40; centers(inds) = []; PWs(inds) = [];
        distD{ii,1} = PWs;
        meanPWs(ii) = (mean(PWs));
    end
    [distDo,allVals] = getAveragesAndAllValues(distD);
   minBin = min(allVals);
   maxBin = max(allVals);
   incr = 2;
   hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
   hold on;
   [ha,hb,hca] = plotAverageDistributions(distD,'colors',{'k'},'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
   set(gca,'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
%    legs = {'Linear','Sigmoid','Gaussian',[-1.75 0.3 85 5]};
%    putLegend(gca,legs,tcolors);
   ylim([0 110]);
   pmchar=char(177); any_text = sprintf('%.0f%c%.0f cm/sec',mean(meanPWs),pmchar,std(meanPWs)/sqrt(5)); 
   text(10,15,any_text,'FontSize',6);
   changePosition(gca,[0.12 0.13 -0.3 -0.13]);
    put_axes_labels(gca,{'Tuning Width (cm/sec)',[-10 0 0]},{{'Neurons (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_speed_tuning_width'),600);
    return;
end
%% relationship between centers and PWs
if 1
    hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
   hold on;
   colors = {'k','r','b','c','m'};
   for ii = 5%1:length(ei)
        bcs = speedRs{ii}.bin_centers;
        fitg = speedRs{ii}.fits.gauss;
        [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,bcs(2)-bcs(1));
        inds = centers < 1 | centers > 39 | PWs < 1 | PWs > 40; centers(inds) = []; PWs(inds) = [];
        scatter(centers,PWs,'.');
        [p,S,mu] = polyfit(centers,PWs,1);
        sl(ii) = p(1);
        lineRS(ii) = S.R(1,1);
        PWsY = polyval(p,centers);
        plot(centers,PWsY,'color',colors{ii});
        corrs(ii) = corr(centers',PWs');
    end
   set(gca,'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
   changePosition(gca,[0.2 0.13 -0.2 -0.13]);
    put_axes_labels(gca,{'Pref. Speed (cm/sec)',[-4 0 0]},{{'Tuning Width (cm/sec)'},[0 -15 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Scatter_Centers_TW'),600);
end