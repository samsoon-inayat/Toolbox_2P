function accel_response

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

selContexts = [1 4 6 2 7 3 3 4 4 5 5 3 3 4 4 5 5 0 0];
rasterNames = {'light22T','light22T','light22T','air55T','air55T','airD','airID','airD','airID','airD','airID','airT','airIT','airT','airIT','airT','airIT','motionOnsets','motionOffsets'};
o = get_data(ei,selContexts,rasterNames);

si_air_dist_trials = [6 8 10]; si_no_brake = [6 13 8 15 10 17];
props1 = get_props_Rs(o.Rs,50); 
si = si_air_dist_trials;
resp = props1.good_FR(:,si);

n = 0;
%% Percentage speed responsive
resp = resp_speed;
while 1
    t_resp = cell_list_op(resp,[],'or');
    for an = 1:5
        d.bcs = speedRs{an,4}.bin_centers;
        d.FR = speedRs{an,4}.FR_vs_speed;
        fitg = speedRs{an,4}.fits.gauss; fits = speedRs{an,4}.fits.sigmoid; fitl = speedRs{an,4}.fits.linear;
        d.fFRl = fitl.fitted; d.fFRs = fits.fitted; d.fFRg = fitg.fitted;
        d.cl = fitl.coeffsrs(:,3); d.cs = fits.coeffsrs(:,3); d.cg = fitg.coeffsrs(:,3);
        [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,d.bcs(2)-d.bcs(1));
        inds = centers < 1 | centers > 39 | rs < 0.25 | PWs < 10;% | PWs > 20 | PWs < 10;
        inds = ~inds;
        speed_tuned_cells{an,1} = inds;
%         inds = inds & ~t_resp{an}';
        pR(an) = 100*sum(inds)/length(inds);
        resp_speed{an} = inds';
        [mVals,semVals] = findMeanAndStandardError(pR);
    end
    break;
end
fileName = fullfile(mData.pd_folder,sprintf('%s_tuned_cells',mfilename));
save(fileName,'speed_tuned_cells');
%% visualize the data
if 1
    an = 4;
    d.bcs = speedRs{an,4}.bin_centers;
    d.FR = speedRs{an,4}.FR_vs_accel;
    fitg = speedRs{an,4}.fits.gauss; fits = speedRs{an,4}.fits.sigmoid; fitl = speedRs{an,4}.fits.linear;
    d.fFRl = fitl.fitted; d.fFRs = fits.fitted; d.fFRg = fitg.fitted;
    d.cl = fitl.coeffsrs(:,3); d.cs = fits.coeffsrs(:,3); d.cg = fitg.coeffsrs(:,3);
    [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,d.bcs(2)-d.bcs(1));
    inds = rs < 0.5;% | PWs > 20 | PWs < 10;
    inds = ~inds;
%     inds = resp_speed{an,6};
%     t_resp = cell_list_op(resp,[],'or');
%     inds = ~inds;
%     inds = inds & ~t_resp{an}';
    100*sum(inds)/length(inds)
    d.FR = d.FR(inds,:); d.fFRl = d.fFRl(inds,:); d.fFRs = d.fFRs(inds,:); d.fFRg = d.fFRg(inds,:);
    d.cl = d.cl(inds); d.cs = d.cs(inds); d.cg = d.cg(inds);
    generic_multi_plot(1000,[3,4,size(d.FR,1)],'plotSpeedTuning',d)
    return;
end
%%
while 1
    an = 4;
    clear d
    d.bcs = speedRs{an,4}.bin_centers;
    d.FR = speedRs{an,4}.FR_vs_accel;
    fitg = speedRs{an,4}.fits.gauss; fits = speedRs{an,4}.fits.sigmoid; fitl = speedRs{an,4}.fits.linear;
    d.fFRl = fitl.fitted; d.fFRs = fits.fitted; d.fFRg = fitg.fitted;
    d.cl = fitl.coeffsrs(:,3); d.cs = fits.coeffsrs(:,3); d.cg = fitg.coeffsrs(:,4);
    [rs2,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,d.bcs(2)-d.bcs(1));
    inds = rs2 < 0.5;
    inds = ~inds;
    inds_cn = find(inds);
    100*sum(inds)/length(inds)
    d.FR = d.FR(inds,:); d.fFRl = d.fFRl(inds,:); d.fFRs = d.fFRs(inds,:); d.fFRg = d.fFRg(inds,:);
    d.cl = d.cl(inds); d.cs = d.cs(inds); d.cg = d.cg(inds);
    cell_inds = [49 7 49];
%     cell_inds = [1 2 3]+(5*3);
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 3],...
        'spaceRowsCols',[0.15 0.06],'rightUpShifts',[0.14 0.28],'widthHeightAdjustment',...
        [-100 -415]);    set(gcf,'color','w'); set(gcf,'Position',[10 4 3 1]);
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
            xlabel('Acceleration (cm/sec^2)')
        end
        if iii == 2
            legs = {'Linear','Sigmoid','Gaussian',[25 .0 0.06 0.01]};
           putLegend(gca,legs,{'k','m','b'});
        end
        box off;
        set(gca,'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('sample_speed_tuning'),600);
    break;
end
%% see the distribution of rs for linear, sigmoid, and gaussian fitting.
if 1
    var_names = {'linear','sigmoid','gauss'};
   tcolors = {'k','m','b'};
   distD = [];
   ind = [3 3 4];
   for ii = 1:length(ei)
       ii
       for vv = 1:length(var_names)
           cmdTxt = sprintf('distD{ii,vv} = speedRs{ii,4}.fits.%s.coeffsrs(:,ind(vv));',var_names{vv});
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
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_accel_rs'),600);
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
        'ySpacing',0.15,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    ylims = ylim;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[ylims(1) maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = [xdata(1:end)]; xticklabels = {'Linear','Sigmoid','Gaussian'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30)
    changePosition(gca,[0.11 0.03 -0.3 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'R-Squared'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('mean_Rsquared_accel'),600);
end

%% distribution preferred speed
if 1
    distD = [];
    for ii = 1:length(ei)
        bcs = speedRs{ii,4}.bin_centers;
        fitg = speedRs{ii,4}.fits.gauss;
        [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,bcs(2)-bcs(1));
%         inds = centers < 1 | centers > 39 | PWs < 1 | PWs > 40; 
        inds = ~resp_speed{ii,5};
        centers(inds) = [];
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
   pmchar=char(177); any_text = sprintf('%.0f%c%.0f cm/sec^2',mean(meanPWs),pmchar,std(meanPWs)/sqrt(5)); 
   text(1,100,any_text,'FontSize',6);
   changePosition(gca,[0.12 0.13 -0.3 -0.13]);
    put_axes_labels(gca,{'Pref. Acceleration (cm/sec^2)',[0 0 0]},{{'Neurons(%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_accel_preferred'),600);
    return;
end

%% distribution tuning width
if 1
    distD = [];
    for ii = 1:length(ei)
        bcs = speedRs{ii}.bin_centers;
        fitg = speedRs{ii}.fits.gauss;
        [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,bcs(2)-bcs(1));
        inds = centers < 1 | centers > 39 | PWs < 1 | PWs > 40; 
        inds = ~resp_speed{ii,5};
        centers(inds) = []; PWs(inds) = [];
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
   pmchar=char(177); any_text = sprintf('%.0f%c%.0f cm/sec^2',mean(meanPWs),pmchar,std(meanPWs)/sqrt(5)); 
   text(10,15,any_text,'FontSize',6);
   changePosition(gca,[0.12 0.13 -0.3 -0.13]);
    put_axes_labels(gca,{'Tuning Width (cm/sec^2)',[-5 0 0]},{{'Neurons (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_accel_tuning_width'),600);
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
        inds = centers < 1 | centers > 39 | PWs < 1 | PWs > 40; 
        inds = ~resp_speed{ii,1};
        centers(inds) = []; PWs(inds) = [];
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


%% Overlap Indices ImageSC 
while 1
    ntrials = 50;
    si = [Lb_T ArL_L_T Lbs_T Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
    props1A = get_props_Rs(o.Rs(:,si),ntrials);
    respO = [props1A.good_FR];% resp_speed];
    for rr = 1:size(resp_speed,1)
        for cc = 1:size(resp_speed,2)
            resp_speed_abs{rr,cc} = abs(resp_speed{rr,cc});
        end
    end
    resp = [resp_speed_abs respO];% resp_speed];
%     resp = resp_speed_abs;
    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = [{'Sp','SpM'} rasterNamesTxt(si)]; 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 1.5 1.5]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    text(0.5,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.03 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean_motion.pdf',ntrials),600);
    %%
    break;
end

%% agglomerative hierarchical clustering zMIT>zMID
while 1
    mOI1 = mOI;
    mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1);
    tree = linkage(Di);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 2.5 1]);
    set(H,'linewidth',1);
    set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel({'Eucledian','Distance'});%changePosition(hx,[-0.051 0 0]);
    changePosition(gca,[0 0.0 0.07 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster_motion.pdf'),600);
    %%
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

%% speed gauss shuffle
while 1
    speedRG = speedRs(:,3);
    for ii = 1:length(speedRG)
        tsrg = speedRG{ii};
        gaussR = tsrg.fits.gaussR;
        bcs = tsrg.bin_centers;
        [rs,MFR,centers,PWs] = get_gauss_fit_parameters(gaussR.coeffsrs,bcs(2)-bcs(1));
        nshuffle = size(gaussR.coeffsrsS,3);
        for jj = 1:nshuffle
            tcoeff = gaussR.coeffsrsS(:,:,jj);
            [rsS(jj,:),~,~,~] = get_gauss_fit_parameters(tcoeff,bcs(2)-bcs(1));
        end
        rsST = rsS';
        rs_r = repmat(rs',1,nshuffle);
        p_vals = sum(rs_r > rsST,2)/nshuffle;
        sp_cells = p_vals > 0.95;
        for cn = 1:length(sp_cells)
            if sp_cells(cn) == 0
                continue;
            end
            figure(10000);clf;
            plot(bcs,tsrg.FR_vs_speed(cn,:));hold on;
            plot(bcs,gaussR.fitted(cn,:));
            pause(0.1);
        end
        
    end
    %%
    break;
end