function compare_spontaneous_activity

protocol_C = '10_C';
protocol_A = '10_A';
ei_C = evalin('base','ei10_C');
ei_A = evalin('base','ei10_A');
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ET_C = evalin('base',sprintf('ET%s',protocol_C));
ET_A = evalin('base',sprintf('ET%s',protocol_A));
selAnimals_C = 1:length(ei_C);
selAnimals_A = 1:length(ei_A);
n = 0;
%%
out_C = get_spike_rate(ei_C);
out_A = get_spike_rate(ei_A);

n=0;
%%
if 0
    data_C = out_C.m_sp_animal_th;
    data_A = out_A.m_sp_animal_th;
    data = [data_C' data_A'];
    minBin = min([out_C.allVals_th;out_A.allVals_th]);
    maxBin = max([out_C.allVals_th;out_A.allVals_th])+20;
    incr = (maxBin-minBin)/50;
    hf = figure(1002);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 6 2 1.5],'color','w');
    hold on;
    [ha,hb,hca,sigR] = plotDistributions(data,'colors',colors,'maxY',200,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    hold on;
    legs = {'Ctrl','AD'};
    ylim([0 40]);
    xlim([0 maxBin+10]);
%     ylim([0 15]);
    xlim([0 maxBin+10])
    xlims = xlim; dx = xlims(2) - xlims(1); ylims = ylim; dy = ylims(2) - ylims(1);
    legs{length(legs)+1} = [xlims(1)+dx/1.5 dx/30 ylims(1)+dy/3 dy/15];
    putLegend(ha,legs,'colors',colors);
    axes(ha);
    h = xlabel('Spike Rate (Hz)');%changePosition(h,[0 -dy/3 0]);
    h = ylabel('Percentage');changePosition(h,[-0.1 0 0]);
    set(gca,'FontSize',axes_font_size,'FontWeight','Bold');changePosition(ha,[0.07 0.01 -0.05 -0.05]);
    file_name = sprintf('%s_distribution spike rate',mfilename);
    save_pdf(hf,mData.pdf_folder,file_name,600);
end

%%
if 0
    [mVar_C,semVar_C] = findMeanAndStandardError(out_C.m_sp_animal_level_th);
    [mVar_A,semVar_A] = findMeanAndStandardError(out_A.m_sp_animal_level_th);
    mVar = [mVar_C,mVar_A];
    semVar = [semVar_C,semVar_A];
    [h,p,ci,stats] = ttest2(out_C.m_sp_animal_level_th,out_A.m_sp_animal_level_th)
    combs = [1 2]; %p = ones(size(combs,1),1); h = logical(zeros(size(combs,1),1));
    xdata = [1 2]; maxY = 40;
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 2.25 1],'color','w');
    hold on;
    tcolors = {colors{1};colors{1}};
    hbs = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'maxY',maxY,'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 2:2:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+1],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = {'Ctrl','AD'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.02 0.03 0.02 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{'FR (Hz)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bargraph',mfilename),600);
end

%%
if 1
    data_C_motion = out_C.m_sp_animal_level_motion;
    data_A_motion = out_A.m_sp_animal_level_motion;
    data_C_rest = out_C.m_sp_animal_level_rest;
    data_A_rest = out_A.m_sp_animal_level_rest;
    
    dataT_C = array2table(fliplr([data_C_motion' data_C_rest'])); dataT_C.Properties.VariableNames = {'Rest','Motion'};
    dataT_A = array2table(fliplr([data_A_motion' data_A_rest'])); dataT_A.Properties.VariableNames = {'Rest','Motion'};
    groupT = table([ones(size(dataT_C,1),1);2*ones(size(dataT_A,1),1)]); groupT.Properties.VariableNames = {'Group'};
    between = [groupT [dataT_C;dataT_A]]; between.Group = categorical(between.Group);
    within = table([1 2]');
    within.Properties.VariableNames = {'State'};
    within.State = categorical(within.State);
    rmaR = repeatedMeasuresAnova(between,within);
    writetable(between,fullfile(mData.pdf_folder,sprintf('%s_dataT_Spike_Rate_motion_rest.xlsx',mfilename)));
    writetable(rmaR.ranova,fullfile(mData.pdf_folder,sprintf('%s_dataT_Spike_Rate_motion_rest_stats_output.xlsx',mfilename)),'WriteRowNames',true);
    
%     findMeanAndStandardError
    mVar = rmaR.est_marginal_means.Mean;
    semVar = rmaR.est_marginal_means.Formula_StdErr;
    combs = rmaR.combs;
    p = rmaR.p; h = p<0.05;
    xdata = [1 2 3 4]; maxY = 1;
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 2.25 1],'color','w');
    hold on;
    tcolors = {colors{1};colors{2};colors{1};colors{2}};
    hbs = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'maxY',maxY,'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 2:2:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+0],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = {'Ctrl','AD'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.02 0.03 0.02 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{'FR (Hz)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bargraph_overall',mfilename),600);
end

n = 0;
%%
if 0
    data_C = out_C.m_sp_animal_level_th;
    data_A = out_A.m_sp_animal_level_th;
    [h,p,ci,stats] = ttest2(data_C,data_A)
    [mVarT,semVarT] = findMeanAndStandardError(data_C);
    [mVarT_A,semVarT_A] = findMeanAndStandardError(data_A);

    mVar = [mVarT mVarT_A];
    semVar = [semVarT semVarT_A];
    combs = nchoosek(1:2,2);
    
    xdata = [1 2]; maxY = 40;
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3}};
    hbs = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'maxY',maxY,'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 2:2:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+5],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata; xticklabels = {'Ctrl','AD'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.06 0.03 0.02 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{'Spike Rate (Hz)',[0 0 0]});
    file_name = sprintf('%s_bargraph spike rate',mfilename);
    save_pdf(hf,mData.pdf_folder,file_name,600);
end
n=0;
%%
if 1
    data_C_motion = out_C.m_sp_animal_level_motion;
    data_A_motion = out_A.m_sp_animal_level_motion;
    data_C_rest = out_C.m_sp_animal_level_rest;
    data_A_rest = out_A.m_sp_animal_level_rest;
    
    dataT_C = array2table(fliplr([data_C_motion' data_C_rest'])); dataT_C.Properties.VariableNames = {'Rest','Motion'};
    dataT_A = array2table(fliplr([data_A_motion' data_A_rest'])); dataT_A.Properties.VariableNames = {'Rest','Motion'};
    groupT = table([ones(size(dataT_C,1),1);2*ones(size(dataT_A,1),1)]); groupT.Properties.VariableNames = {'Group'};
    between = [groupT [dataT_C;dataT_A]]; between.Group = categorical(between.Group);
    within = table([1 2]');
    within.Properties.VariableNames = {'State'};
    within.State = categorical(within.State);
    rmaR = repeatedMeasuresAnova(between,within);
    writetable(between,fullfile(mData.pdf_folder,sprintf('%s_dataT_Spike_Rate_motion_rest.xlsx',mfilename)));
    writetable(rmaR.ranova,fullfile(mData.pdf_folder,sprintf('%s_dataT_Spike_Rate_motion_rest_stats_output.xlsx',mfilename)),'WriteRowNames',true);
    
    mVar = rmaR.est_marginal_means{:,3}';
    semVar = rmaR.est_marginal_means{:,end}';
    combs = rmaR.combs;
    p = rmaR.p; h = p < 0.05;
    
    xdata = [1 2 3 4]; maxY = 1;
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3}};
    hbs = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'maxY',maxY,'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 3:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+1],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata; xticklabels = {'Ctrl-R','Ctrl-M','AD-R','AD-M'};xtickangle(45)
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.06 0.03 0.02 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{'Spike Rate (Hz)',[0 0 0]});
    file_name = sprintf('%s_bargraph spike rate_motion_rest',mfilename);
    save_pdf(hf,mData.pdf_folder,file_name,600);
end

function out = get_spike_rate(ei_C)
allVals = []; allVals_motion = []; allVals_rest = []; allVals_th = [];
for ii = 1:length(ei_C)
    ei = ei_C{ii};
    spSigAll = nan(length(ei.plane{1}.tP.deconv.caSigAll),length(ei.plane{1}.tP.deconv.caSigAll{1}));
    for cn = 1:length(ei.plane{1}.tP.deconv.caSigAll)
        spSigAll(cn,:) = ei.plane{1}.tP.deconv.spSigAll{cn}';
    end
    iscell = logical(ei.plane{1}.tP.iscell(:,1));
    spSigAll = spSigAll(iscell,:);
    speed = ei.b.fSpeed;
    inds_motion = find(speed > 0); % find indices of array where speed is greater than zero
    inds_rest = find(speed == 0); % rest period
%     signal = spSigAll(1,:);
%     figure(100);clf;plot(signal);hold on;plot(ones(size(signal))*(mean(signal)+3*std(signal)))
    frames_f = ei.plane{1}.b.frames_f;
    [~,ind_frames_motion,~] = intersect(frames_f,inds_motion); % frames where animal was moving
    [~,ind_frames_rest,~] = intersect(frames_f,inds_rest); % frames where animal was at rest
    m_sp_animal{ii} = nanmean(spSigAll,2); % mean spike rate 
    m_sp_animal_motion{ii} = nanmean(spSigAll(:,ind_frames_motion),2); % mean spike rate when animal moving
    m_sp_animal_rest{ii} = nanmean(spSigAll(:,ind_frames_rest),2); % mean spiek rate when animal at rest
    allVals = [allVals;m_sp_animal{ii}]; % accumulating values over animals but for all cells
    allVals_motion = [allVals_motion;m_sp_animal_motion{ii}]; % accumulating values over animals for all cells
    allVals_rest = [allVals_rest;m_sp_animal_rest{ii}];
    m_sp_animal_level(ii) = nanmean(m_sp_animal{ii}); % average over cells
    m_sp_animal_level_motion(ii) = nanmean(m_sp_animal_motion{ii});
    m_sp_animal_level_rest(ii) = nanmean(m_sp_animal_rest{ii});
%     thr = nanmean(spSigAll,2) + 3*nanstd(spSigAll,[],2);
    thr = 3*nanstd(spSigAll,[],2);
    for cn = 1:size(spSigAll,1)
        m_sp_animal_th{ii}(cn,1) = nanmean(spSigAll(cn,spSigAll(cn,:) > thr(cn)));
%         figure(100);clf;plot(spSigAll(cn,:));hold on;plot(ones(size(spSigAll(cn,:)))*thr(cn))
    end
    m_sp_animal_level_th(ii) = nanmean(m_sp_animal_th{ii}); %finding mean over cells
    allVals_th = [allVals;m_sp_animal_th{ii}];
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

out.m_sp_animal_th = m_sp_animal_th;
out.m_sp_animal_level_th = m_sp_animal_level_th;
out.allVals_th = allVals_th;