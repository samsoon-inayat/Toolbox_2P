function firing_rate_motion_vs_rest
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei_C = evalin('base','ei'); 

    selContexts = [3 3 4 4 5 5];
    rasterNames = {'airD','airIT','airD','airIT','airD','airIT'};

    Rs_C = get_rasters_data(ei_C,selContexts,rasterNames);
    typeP = 'Population';
    thr = -1;
    fileName = fullfile(mData.pd_folder,sprintf('%s_cond_345',mfilename));
    if 0
        out_CC345 = get_spike_rate_conditions_345(ei_C,Rs_C,thr);
        out_CC = get_spike_rate_conditions(ei_C,thr);
        save(fileName,'out_CC345','out_CC');
    else
        temp = load(fileName);
        out_CC345 = temp.out_CC345;
        out_CC = temp.out_CC;
    end
    %%
    fileName = fullfile(mData.pd_folder,sprintf('%s_%s',mfilename,typeP));
    if 0
        out_C = get_spike_rate(ei_C,thr);
       save(fileName,'out_C','thr');
    %    save(fileName,'out_C','out_CC','thr');
    else
        temp = load(fileName);
        out_C = temp.out_C;
    %     out_CC = temp.out_CC;
        thr = temp.thr;
    end
%         out_CC = get_spike_rate_conditions(ei_C,thr);

    tcolors = {'k','r','k','r'};
    n=0;
    break
end
n = 0;
%%
while 1
    for ii = 1:length(ei_C)
        ei = ei_C{ii};
        tspSigAll = [];
        for pp = 1:length(ei.plane)
            this_cell_list = logical(ei.plane{pp}.tP.iscell(:,1));
            tspSigAll = [tspSigAll;ei.plane{pp}.tP.deconv.spSigAll(this_cell_list,:)];
        end
        mean_tsp = mean(tspSigAll,2);
        perc_silent(ii) = 100*sum((mean_tsp < 0.1))/length(mean_tsp);
    end
    break;
end
%% Spike Rate Conditions 345 Trials Inter-trials
while 1
    [within,dvn,xlabels] = make_within_table({'Conditions','TI'},[3,2]);
    dataT = make_between_table({out_CC345.mean_msr},dvn);
    ra = RMA(dataT,within);
    ra.ranova
%     [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Conditions_by_TI','hsd'},[1 1 1]);
    
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 9 1.5 1],'color','w');

    hold on;
    tcolors = [mData.colors(3); mData.colors(3); mData.colors(4); mData.colors(4); mData.colors(5); mData.colors(5);];
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.0001);
    make_bars_hollow(hbs(2:2:end));
    rectangle(gca,'Position',[0.75 maxY+0.01 1 0.025],'edgecolor','k','facecolor','k');
    text(1.85,maxY+0.02,'Trials','FontSize',5);
    rectangle(gca,'Position',[4 maxY+0.01 1 0.025],'edgecolor','k');
    text(5.2,maxY+0.02,'Inter-Trials','FontSize',5);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+0.051],'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
    xticks = [1.5 4.5 7.5]; xticklabels = {'Ar','ArL','Ar*'};
    xtickangle(30)
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.1 0.05 -0.08 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Avg. FR (AU)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Firing_Rate_conditions345'),600);
break;
end

%% Ca Conditions 345 Trials Inter-trials
while 1
    [within,dvn,xlabels] = make_within_table({'Conditions','TI'},[3,2]);
    dataT = make_between_table({out_CC345.mean_mcr*100},dvn);
    ra = RMA(dataT,within);
    ra.ranova
%     ra = repeatedMeasuresAnova(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Conditions_by_TI','hsd'},[1 1 1]);
%     [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
    
    
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 9 1.35 1],'color','w');

    hold on;
    tcolors = [mData.colors(3); mData.colors(3); mData.colors(4); mData.colors(4); mData.colors(5); mData.colors(5);];
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',100,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.0001);
    make_bars_hollow(hbs(2:2:end));
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
    xticks = xdata; xticklabels = {'C3-T','C3-IT','C4-T','C4-IT','C3''-T','C3''-IT'};
    xtickangle(30)
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.1 0.05 -0.08 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{sprintf('Avg. %cF/Fo (%%)',916)},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Ca_Signal_conditions345'),600);
break;
end
%%
while 1
    [within,dvn,xlabels] = make_within_table({'Conditions'},7);
    dataT = make_between_table({out_CC.mean_mcr*100},dvn);
    ra = RMA(dataT,within);
%     ra = repeatedMeasuresAnova(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Conditions','hsd'},[1 1 1]);
%     [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
    
    
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 9 1.35 1],'color','w');

    hold on;
%     tcolors ={colors{1};colors{2};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',colors,'sigColor','k',...
        'ySpacing',100,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.0001);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
    xticks = xdata; xticklabels = {'Lb','Ab','Ar','ArL','Ar*','Lb*','Ab*'};
    xtickangle(30)
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.1 0.05 -0.08 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{sprintf('Avg. %cF/Fo (%%)',916)},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Ca_Signal_conditions'),600);
break;
end
%%
while 1
    [within,dvn,xlabels] = make_within_table({'Conditions'},7);
    dataT = make_between_table({out_CC.mean_msr},dvn);
    ra = RMA(dataT,within);
%     [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Conditions','hsd'},[1 1 1]);
    
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 9 1.5 1],'color','w');

    hold on;
%     tcolors ={colors{1};colors{2};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',colors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.0001);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
    xticks = xdata; xticklabels = {'Lb','Ab','Ar','ArL','Ar*','Lb*','Ab*'};
    xtickangle(30)
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.1 0.05 -0.08 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Avg. FR (AU)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Firing_Rate_conditions'),600);
break;
end

%%
while 1
    
    s = generate_shades(4);
    tcolors = s.g;
    [within,dvn,xlabels] = make_within_table({'Type'},2);
    dataT = make_between_table({out_C.m_sp_animal_level_rest',out_C.m_sp_animal_level_motion'},dvn);
    ra = RMA(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Type','hsd'},[1 1 1]);
    
     xdata = [1 1.75]; 
    colors = mData.colors;
        hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 9 1.25 1],'color','w');

    hold on;
%     tcolors ={colors{1};colors{2};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.4,'sigLinesStartYFactor',0.0001);
    set(gca,'xlim',[0.5 xdata(end)+0.5],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
    xticks = xdata; xticklabels = {'Rest','Run'};
    hatch(hbs(2),30,'k','-',4,0.1); %angle,color,style,step,width
%     hatch(hbs(1),150,'k','-',2,0.1); %angle,color,style,step,width
    xtickangle(30)
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.17 0.05 -0.4 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Avg. FR (AU)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Firing_Rate_Motion_vs_Rest'),600);
break;
end

%%
while 1
    
    s = generate_shades(4);
    tcolors = s.g;
    [within,dvn,xlabels] = make_within_table({'Type'},2);
    dataT = make_between_table({out_C.m_sp_animal_level_rest',out_C.m_sp_animal_level_motion'},dvn);
    ra = repeatedMeasuresAnova(dataT,within);

    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
     xdata = [1 1.75]; 
    colors = mData.colors;
        hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 9 1.25 1],'color','w');

    hold on;
%     tcolors ={colors{1};colors{2};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.4,'sigLinesStartYFactor',0.0001);
    set(gca,'xlim',[0.5 xdata(end)+0.5],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
    xticks = xdata; xticklabels = {'Rest','Run/Walk'};
    hatch(hbs(2),30,'k','-',4,0.1); %angle,color,style,step,width
%     hatch(hbs(1),150,'k','-',2,0.1); %angle,color,style,step,width
    xtickangle(30)
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.17 0.05 -0.4 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Avg. FR (AU)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Firing_Rate_Motion_vs_Rest'),600);
break;
end


%%
while 1
    tcolors = {'k','r'};
   distD(:,1) = out_C.m_sp_animal_motion;
   distD(:,2) = out_C.m_sp_animal_rest;
   [distDo,allVals,allValsG] = plotDistributions(distD);
   minBin = min(allVals);
   maxBin = max(allVals);
   [h,p,ks2stat] = kstest2(allValsG{1},allValsG{2});
   [h,p,cd,ks2stat] = ttest2(allValsG{1},allValsG{2});
    tcolors = mData.colors;
   incr = 0.001; %maxBin =
   hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.85 1],'color','w');
   hold on;
%    [ha,hb,hca] = plotDistributions(distD,'colors',tcolors,'maxY',maxBin,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
   [ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
   set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
   changePosition(gca,[0.09 0.13 -0.05 -0.13]);
    put_axes_labels(gca,{'Average Firing Rate (AU)',[0 0 0]},{{'Percentage','of Neurons'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_firing_rate'),600);
    break;
end

function out = get_spike_rate_conditions_345(ei_C,Rs,thr)
msr_cc = []; mcr_cc = [];
for rr = 1:size(Rs,1)
    ei = ei_C{rr};
    b = ei.b;
    nsamp = round(1e6 * 1/b.si);
    for cc = 1:size(Rs,2)
        tRs = Rs{rr,cc};
        C_s = tRs.onsets;
        C_e = tRs.offsets;
        msr = []; mcr = [];
        for cn = 1:length(C_s)
            tts = C_s(cn); tte = C_e(cn);
            tmsr = []; tmcr = [];
            for pp = 1:length(ei.plane)
                tspSigAll = []; caSig = [];
                this_cell_list = logical(ei.plane{pp}.tP.iscell(:,1));
                tspSigAll = ei.plane{pp}.tP.deconv.spSigAll(this_cell_list,:);
                caSig = get_calcium_data(ei,pp);
                mask = tspSigAll > thr;
                tspSigAll(~mask) = NaN;
                caSig(~mask ) = NaN;
                frames_f = ei.plane{pp}.b.frames_f;
                frameNums = find(frames_f > tts & frames_f < tte);
                tmsr = [tmsr;nanmean(tspSigAll(:,frameNums),2)];
                tmcr = [tmcr;nanmean(caSig(:,frameNums),2)];
            end
            msr(cn) = mean(tmsr);
            mcr(cn) = mean(tmcr);
        end
        msr_cc(rr,cc) = mean(msr);
        mcr_cc(rr,cc) = mean(mcr);
    end
end
out.mean_msr = msr_cc;
out.mean_mcr = mcr_cc;

function out = get_spike_rate_conditions(ei_C,thr)
allVals = []; allVals_motion = []; allVals_rest = []; allVals_th = [];
allVals_an = [];
for ii = 1:length(ei_C)
    ei = ei_C{ii};
    b = ei.b;
    nsamp = round(1e6 * 1/b.si);
    C_s = [b.stim_r(1);b.air_puff_r(1);b.air_puff_r(11);...
        b.air_puff_r(21);b.air_puff_r(31);b.stim_r(21);b.air_puff_r(41)]-nsamp;
    C_e = [b.stim_r(10);b.air_puff_r(10);b.air_puff_r(20);...
        b.air_puff_r(30);b.air_puff_r(40);b.stim_r(30);b.air_puff_r(50)]+nsamp;
    for cn = 1:length(C_s)
        tts = C_s(cn); tte = C_e(cn);
        tmsr = []; tmcr = [];
        for pp = 1:length(ei.plane)
            tspSigAll = []; caSig = [];
            this_cell_list = logical(ei.plane{pp}.tP.iscell(:,1));
            tspSigAll = ei.plane{pp}.tP.deconv.spSigAll(this_cell_list,:);
            caSig = get_calcium_data(ei,pp);
            mask = tspSigAll > thr;
            tspSigAll(~mask) = NaN;
            caSig(~mask ) = NaN;
            frames_f = ei.plane{pp}.b.frames_f;
            frameNums = find(frames_f > tts & frames_f < tte);
            tmsr = [tmsr;nanmean(tspSigAll(:,frameNums),2)];
            tmcr = [tmcr;nanmean(caSig(:,frameNums),2)];
        end
        msr{ii,cn} = tmsr;
        mcr{ii,cn} = tmcr;
    end
end
out.mean_msr = arrayfun(@(x) nanmean(x{1}),msr);
out.msr = msr
out.mean_mcr = arrayfun(@(x) nanmean(x{1}),mcr);
out.mcr = mcr;

function out = get_spike_rate(ei_C,thr)
allVals = []; allVals_motion = []; allVals_rest = []; allVals_th = [];
allVals_an = [];
for ii = 1:length(ei_C)
    ei = ei_C{ii};
    
%     speed = ei.b.air_puff_raw;
%     inds_motion = find(speed > 0.5);
%     inds_rest = find(speed < 0.5);
    
    speed = ei.b.fSpeed;
%     inds_motion = find(speed > 0);
    inds_motion = find(speed ~= 0);
    inds_rest = find(speed == 0);
    onset_first_trial = ei.b.air_puff_r(1)-1000;
    offset_last_trial = ei.b.air_puff_f(end)+1000;
    inds_motion1 = inds_motion(find(inds_motion > onset_first_trial & inds_motion < offset_last_trial));
    inds_rest1 = inds_rest(find(inds_rest > onset_first_trial & inds_rest < offset_last_trial));
    percent_motion(ii) = 100*length(inds_motion1)/(offset_last_trial - onset_first_trial);
    percent_rest(ii) = 100*length(inds_rest1)/(offset_last_trial - onset_first_trial);
    
    spSigAll = [];
    all_spSigAll = [];
    for pp = 1:length(ei.plane)
        tspSigAll = [];
        this_cell_list = logical(ei.plane{pp}.tP.iscell(:,1));
        tspSigAll = ei.plane{pp}.tP.deconv.spSigAll(this_cell_list,:);
%         for cn = 1:length(this_cell_list)
%             tspSigAll(cn,:) = ei.plane{pp}.tP.deconv.spSigAll{this_cell_list(cn)}';
%         end
        mask = tspSigAll > thr;
        tspSigAll(~mask) = NaN;
        spSigAll{pp} = tspSigAll;
        all_spSigAll = [all_spSigAll;tspSigAll];
        frames_f = ei.plane{pp}.b.frames_f;
        [~,ind_frames_motion{pp},~] = intersect(frames_f,inds_motion);
        [~,ind_frames_rest{pp},~] = intersect(frames_f,inds_rest);
    end
    tm_sp_animal_motion = []; tm_sp_animal_rest =[];
    for pp = 1:length(ei.plane)
        tspSigAll = spSigAll{pp};
        temp_frames = ind_frames_motion{pp}; temp_frames = temp_frames(temp_frames <= size(tspSigAll,2));
        tm_sp_animal_motion = [tm_sp_animal_motion;nanmean(tspSigAll(:,temp_frames),2)];
        temp_frames = ind_frames_rest{pp}; temp_frames = temp_frames(temp_frames <= size(tspSigAll,2));
        tm_sp_animal_rest = [tm_sp_animal_rest;nanmean(tspSigAll(:,temp_frames),2)];
    end
    m_sp_animal{ii} = nanmean(all_spSigAll,2);
    m_sp_animal_motion{ii} = tm_sp_animal_motion;
    m_sp_animal_rest{ii} = tm_sp_animal_rest ;

    allVals = [allVals;m_sp_animal{ii}];
    allVals_an{ii} = all_spSigAll;
    allVals_motion = [allVals_motion;m_sp_animal_motion{ii}];
    allVals_rest = [allVals_rest;m_sp_animal_rest{ii}];
    m_sp_animal_level(ii) = nanmean(m_sp_animal{ii});
    m_sp_animal_level_motion(ii) = nanmean(m_sp_animal_motion{ii});
    m_sp_animal_level_rest(ii) = nanmean(m_sp_animal_rest{ii});
%     thr = nanmean(spSigAll,2) + 3*nanstd(spSigAll,[],2);
%     for cn = 1:size(spSigAll,1)
%         m_sp_animal_th{ii}(cn,1) = nanmean(spSigAll(cn,spSigAll(cn,:) > thr(cn)));
%     end
%     m_sp_animal_level_th(ii) = nanmean(m_sp_animal_th{ii});
%     allVals_th = [allVals;m_sp_animal_th{ii}];
end
out.m_sp_animal = m_sp_animal;
out.allVals = allVals;
out.m_sp_animal_level = m_sp_animal_level;
out.m_sp_animal_motion = m_sp_animal_motion;
out.allVals_motion = allVals_motion;
out.m_sp_animal_level_motion = m_sp_animal_level_motion;
out.m_sp_animal_rest = m_sp_animal_rest;
out.allVals_rest = allVals_rest;
out.m_sp_animal_level_rest = m_sp_animal_level_rest;
out.percent_motion = percent_motion;
out.allVals_an = allVals_an;

% out.m_sp_animal_th = m_sp_animal_th;
% out.m_sp_animal_level_th = m_sp_animal_level_th;
% out.allVals_th = allVals_th;