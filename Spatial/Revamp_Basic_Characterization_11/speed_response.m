function speed_response

%% visualize the data
while 1
    an = 4;
    d.bcs = speedRs{an}.bin_centers;
    d.FR = speedRs{an}.FR_vs_speed;
    fitg = speedRs{an}.fits.gauss; fits = speedRs{an}.fits.sigmoid; fitl = speedRs{an}.fits.linear;
    d.fFRl = fitl.fitted; d.fFRs = fits.fitted; d.fFRg = fitg.fitted;
    d.cl = fitl.coeffsrs(:,3); d.cs = fits.coeffsrs(:,3); d.cg = fitg.coeffsrs(:,3);
    [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,d.bcs(2)-d.bcs(1));
    inds = centers < 1 | centers > 39 | rs < 0.25 | PWs < 10;% | PWs > 20 | PWs < 10;
    inds = logical(resp_speed{an,4});
%     t_resp = cell_list_op(resp,[],'or');
%     inds = ~inds;
%     inds = inds & ~t_resp{an}';
    100*sum(inds)/length(inds)
    d.FR = d.FR(inds,:); d.fFRl = d.fFRl(inds,:); d.fFRs = d.fFRs(inds,:); d.fFRg = d.fFRg(inds,:);
    d.cl = d.cl(inds); d.cs = d.cs(inds); d.cg = d.cg(inds);
    generic_multi_plot(1000,[3,4,size(d.FR,1)],'plotSpeedTuning',d)
    break;
end

%% visualize the data
while 1
    an = 4;
    clear d;
    d.bcs = speedRs{an,3}.bin_centers;
    d.FR = speedRs{an,3}.FR_vs_speed;
    fitg = speedRs{an,3}.fits.gaussR; 
    d.fFRg = fitg.fitted; d.cg = fitg.coeffsrs(:,3);
%     [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,d.bcs(2)-d.bcs(1));
    inds = resp_speedAcc{an,1};
    100*sum(inds)/length(inds)
    d.FR = d.FR(inds,:); d.fFRg = d.fFRg(inds,:); %d.PWs = PWs(inds);
%     d.cg = post_do_gauss_fit(d.fFRg);     [rs,MFR,centers,PWs] = get_gauss_fit_parameters(d.cg,d.bcs(2)-d.bcs(1));
    [~,~,centersO,PWsO] = get_gauss_fit_parameters(fitg.coeffsrs(inds,:),1);
    d.PWs = centersO;
%     d.cg = d.cg(inds);
    generic_multi_plot(1000,[3,4,size(d.FR,1)],'plotSpeedTuningG',d)
    break;
end

%%
while 1
    an = 4;
%     resp_sp = resp_speed(:,[1 3]);
    
    sMcN = speedRs{an,2};
    resp_M = resp_speed(an,3);
    d.bcs = speedRs{an}.bin_centers;
    d.FR = speedRs{an}.FR_vs_speed;
    fitg = speedRs{an}.fits.gauss; fits = speedRs{an}.fits.sigmoid; fitl = speedRs{an}.fits.linear;
    d.fFRl = fitl.fitted; d.fFRs = fits.fitted; d.fFRg = fitg.fitted;
    d.cl = fitl.coeffsrs(:,3); d.cs = fits.coeffsrs(:,3); d.cg = fitg.coeffsrs(:,4);
    [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,d.bcs(2)-d.bcs(1));
    inds = centers < 1 | centers > 39 | rs < 0.25 | PWs < 10;% | PWs > 20 | PWs < 10;
    inds = ~inds;
    resp_sp_and = cell_list_op([resp_M,{inds'}],[],'and'); %resp_sp_and = resp_sp_and(:,1);
    inds = resp_speed{an,4};
    inds_cn = find(inds);
    100*sum(inds)/length(inds)
    d.FR = d.FR(inds,:); d.fFRl = d.fFRl(inds,:); d.fFRs = d.fFRs(inds,:); d.fFRg = d.fFRg(inds,:);
    d.cl = d.cl(inds); d.cs = d.cs(inds); d.cg = d.cg(inds);
    cell_inds = [8 13 10];
    cell_inds = [32 41 24];
%     cell_inds = [1 2 3]+(5*3);
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 3],...
        'spaceRowsCols',[0.15 0.065],'rightUpShifts',[0.15 0.28],'widthHeightAdjustment',...
        [-100 -415]);    set(gcf,'color','w'); set(gcf,'Position',[10 4 3 1]);
%     ff1 = makeFigureRowsCols(2021,[0.5 0.5 4 1],'RowsCols',[1 3],...
%         'spaceRowsCols',[0.15 0.06],'rightUpShifts',[0.14 0.26],'widthHeightAdjustment',...
%         [-100 -415]);    set(gcf,'color','w'); set(gcf,'Position',[10 6 3 1]);
    for iii = 1:3
        ii = cell_inds(iii);
        rs = [d.cl(ii) d.cs(ii) d.cg(ii)];
        [~,mind] = max(rs);
        axes(ff.h_axes(1,iii));
        plot(d.bcs,d.FR(ii,:),'c');hold on;
        pl(1) = plot(d.bcs,d.fFRl(ii,:),'k','linewidth',1);
        pl(2) = plot(d.bcs,d.fFRs(ii,:),'m','linewidth',0.5);
        pl(3) = plot(d.bcs,d.fFRg(ii,:),'b','linewidth',1);
        ylims = ylim;
        ylim([0 ylims(2)]);
%         set(pl(mind),'linewidth',2);
%         title(sprintf('Cell %d',ii));
        if iii > 1
%             set(gca,'YTick',[]);
        else
            ylabel({'Firing','Rate (AU)'});
        end
        if iii == 2
            xlabel('Speed (cm/sec)')
        end
        if iii == 1
            legs = {'Linear','Sigmoid','Gaussian',[18 2 0.27 0.02]};
           putLegend(gca,legs,{'k','m','b'});
        end
        box off;
        set(gca,'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
        format_axes(gca);
%         iic = inds_cn(ii);
%         axes(ff1.h_axes(1,iii));
%         if sMcN.speed_resp(iic) == 0
%             continue;
%         end
%        if sMcN.speed_resp(iic) == 1
%             velocity_tuning = sMcN.speed_tuning_inc(iic,:);
%         end
%         if sMcN.speed_resp(iic) == -1
%             velocity_tuning = sMcN.speed_tuning_dec(iic,:);
%         end
%         [p,S,mu] = polyfit(sMcN.bins,velocity_tuning,1);
%         [f_vt,delta] = polyval(p,sMcN.bins,S,mu);
%         plot(sMcN.bins,velocity_tuning,'.','markersize',4);hold on;
%         plot(sMcN.bins,f_vt);
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('sample_speed_tuning'),600);
    break;
end
%% see the distribution of rs for linear, sigmoid, and gaussian fitting.
while 1
    var_names = {'linear','sigmoid','gauss'};
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
           infinds = find(distD{ii,vv} < -20);
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
    break;
end

%% mean rs values
while 1
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
    format_axes(gca);
    save_pdf(hf,mData.pdf_folder,sprintf('mean_Rsquared'),600);
    break;
end

%% distribution preferred speed
while 1
    distD = [];
    for ii = 1:length(ei)
        bcs = speedRs{ii,3}.bin_centers;
        fitg = speedRs{ii,3}.fits.gaussR; 
        inds = (resp_speedAcc{ii,1}); 
        [~,~,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs(inds,:),1);
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
   text(1,103,any_text,'FontSize',6);
   changePosition(gca,[0.12 0.13 -0.3 -0.13]);
    put_axes_labels(gca,{'Pref. Speed (cm/sec)',[-2 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_speed_preferred'),600);
    format_axes(gca);
    break;
end

%% distribution tuning width
while 1
    distD = [];
    for ii = 1:length(ei)
        bcs = speedRs{ii,3}.bin_centers;
        fitg = speedRs{ii,3}.fits.gaussR;
        inds = (resp_speedAcc{ii,1}); 
        [~,~,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs(inds,:),1);
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
   text(5,10,any_text,'FontSize',6);
   changePosition(gca,[0.12 0.13 -0.3 -0.13]);
    put_axes_labels(gca,{'Tuning Width (cm/sec)',[-2 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_speed_tuning_width'),600);
    format_axes(gca);
    break;
end

%% exploring McN speed curves
while 1
    %%
    an = 4; 
    sMcN = speedRs{an,2};
    old = speedRs{an,1}.fits;
    bins = sMcN.bins;
    for ii = 1:size(sMcN.speed_resp,1)
        if sMcN.speed_resp(ii) == 0
            continue;
        end
       if sMcN.speed_resp(ii) == 1
            velocity_tuning = sMcN.speed_tuning_inc(ii,:);
        end
        if sMcN.speed_resp(ii) == -1
            velocity_tuning = sMcN.speed_tuning_dec(ii,:);
        end
        [p,S,mu] = polyfit(bins,velocity_tuning,1);
        [f_vt,delta] = polyval(p,bins,S,mu);
        figure(100000);clf;plot(bins,velocity_tuning,'.');hold on;
        plot(bins,f_vt);
        pause(0.1);
    end
    %%
    break;
end

%% Percentage of excitatory inhibitory responsive cells
while 1
    si = [MOn_T MOff_T];
    Rs = o.Rs(:,si);ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    
    respME = cell_list_op(props1.exc,[],'or');     respMI = cell_list_op(props1.inh,[],'or');
    respM = cell_list_op(respME,respMI,'or');
    
    resp1 = [respM(:,1) resp_speedAcc];
    presp1 = 100 * exec_fun_on_cell_mat(resp1,'sum')./exec_fun_on_cell_mat(resp1,'length');
    
    [within,dvn,xlabels] = make_within_table({'var'},[3]);
    dataT = make_between_table({presp1},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'var','hsd'},[1 1 1]);
    xdata = make_xdata([3],[1 2]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);

    tcolors = mData.dcolors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Motion','Speed','Accel'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.06 0.0 -0.4 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Tuned Cells (%)',[0 0 0]});
    
    save_pdf(hf,mData.pdf_folder,'motion_responsive_exc_inh',600);
    %%
    break;
end
