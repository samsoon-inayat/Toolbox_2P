function firing_rate_motion_vs_rest

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;

typeP = 'Population';cellsOrNot = 1; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN; FR = NaN;
typeP = 'Place';cellsOrNot = 1; planeNumber = NaN; zMI_Th = 1.96; fwids = [1 150]; fcens = [1 150]; rs_th = 0.3; FR = [0.1 5000];

conditionsAndRasterTypes = [11;21;31;41];
% conditionsAndRasterTypes = [11 21 31 41];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN,FR);
out = read_data_from_base_workspace(selC)

ei_C = out.eis{1}; ei_A = out.eis{2};
pMs_C = out.pMs{1}; pMs_A = out.pMs{2};
paramMs_C = out.paramMs{1}; paramMs_A = out.paramMs{2};
selAnimals_C = out.selAnimals{1}; selAnimals_A = out.selAnimals{2};
perc_cells_C = out.perc_cells{1}; perc_cells_A = out.perc_cells{2};
thr = -1;
n = 0;
%%
fileName = fullfile(mData.pd_folder,sprintf('%s_%s',mfilename,typeP));
if 0
    out_C = get_spike_rate(ei_C,pMs_C,thr);
    out_A = get_spike_rate(ei_A,pMs_A,thr);
    save(fileName,'out_C','out_A','thr');
else
    temp = load(fileName);
    out_C = temp.out_C;
    out_A = temp.out_A;
    thr = temp.thr;
end

n=0;
%%
if 0
    perc_cells_or_C = out_C.m_sp_animal_level;
    perc_cells_or_A = out_A.m_sp_animal_level;
    [h,p,ci,stats] = ttest2(perc_cells_or_C,perc_cells_or_A)
    
    mVar = [mean(perc_cells_or_C) mean(perc_cells_or_A)]; semVar = [std(perc_cells_or_C)/sqrt(3) std(perc_cells_or_A)/sqrt(5)];
    combs = [1 2]; %p = ra.mcs.p; h = ra.mcs.p < 0.05;
    xdata = [1:2];
%     xdata = [1 2 3 4];
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors ={colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 5:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'RSEG','PSEG'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0.2 0 -0.4 -0.09]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Overall Firing','Rate (Hz)'},[0 0 0]});
    
    save_pdf(hf,mData.pdf_folder,sprintf('Firing_Rate_Motion_vs_Rest_overall'),600);
return;
end
%%
if 0
    perc_cells_or_C = out_C.percent_motion;
    perc_cells_or_A = out_A.percent_motion;
    [h,p,ci,stats] = ttest2(perc_cells_or_C,perc_cells_or_A)
    
    mVar = [mean(perc_cells_or_C) mean(perc_cells_or_A)]; semVar = [std(perc_cells_or_C)/sqrt(3) std(perc_cells_or_A)/sqrt(5)];
    combs = [1 2]; %p = ra.mcs.p; h = ra.mcs.p < 0.05;
    xdata = [1:2];
%     xdata = [1 2 3 4];
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors ={colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 5:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'RSEG','PSEG'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0.2 0 -0.4 -0.09]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Percent','Motion (%)'},[0 0 0]});
    
    save_pdf(hf,mData.pdf_folder,sprintf('Firing_Rate_Percent_Motion'),600);
return;
end
%%
if 1
    data_C = [out_C.m_sp_animal_level_motion' out_C.m_sp_animal_level_rest'];
    data_A = [out_A.m_sp_animal_level_motion' out_A.m_sp_animal_level_rest'];
    data = [data_C;data_A];
    dataT = [table([ones(size(data_C,1),1);(2*ones(size(data_A,1),1))]) array2table(data)];
    dataT.Properties.VariableNames = {'Group','Air','NoAir'};
    dataT.Group = categorical(dataT.Group)

    within = table([1;2]);
    within.Properties.VariableNames = {'Type'};
    within.Type = categorical(within.Type);
    ra = repeatedMeasuresAnova(dataT,within);

    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
     xdata = [1 2 4 5]; 
    colors = mData.colors;
    if isnan(zMI_Th)
       hf = figure(7);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    else
        hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 9 1.5 1],'color','w');
    end
    hold on;
    tcolors ={colors{1};colors{2};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    for ii = 3:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'Motion','Rest','Motion','Rest'};
    xtickangle(15)
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.075 0.0 0 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Firing Rate (Hz)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Firing_Rate_Motion_vs_Rest'),600);
return;
end


function out = get_spike_rate(ei_C,pMs_C,thr)
allVals = []; allVals_motion = []; allVals_rest = []; allVals_th = [];
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
    
    cell_list = pMs_C{1}.all_cns{ii};
    spSigAll = [];
    all_spSigAll = [];
    for pp = 1:length(ei.plane)
        tspSigAll = [];
        this_cell_list = cell_list(cell_list(:,2) == pp,3);
        for cn = 1:length(this_cell_list)
            tspSigAll(cn,:) = ei.plane{pp}.tP.deconv.spSigAll{this_cell_list(cn)}';
        end
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

% out.m_sp_animal_th = m_sp_animal_th;
% out.m_sp_animal_level_th = m_sp_animal_level_th;
% out.allVals_th = allVals_th;