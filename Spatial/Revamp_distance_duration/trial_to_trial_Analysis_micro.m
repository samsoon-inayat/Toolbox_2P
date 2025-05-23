function trial_to_trial_Analysis

%% load data

ntrials = 40;
si = [Lb_T Ab_t_T Ab_i_T Ar_t_T Ar_i_T ArL_t_T ArL_i_T Ars_t_T Ars_i_T Lbs_T Abs_t_T Abs_i_T];
si = [Ab_t_T Ab_i_T Ar_t_T Ar_i_T ArL_t_T ArL_i_T Ars_t_T Ars_i_T Abs_t_T Abs_i_T];
si = [Ab_t_T Ab_i_T Ar_t_D Ar_i_T ArL_t_D ArL_i_T Ars_t_D Ars_i_T Abs_t_T Abs_i_T];
% si = [Ar_t_T Ar_i_T ArL_t_T ArL_i_T Ars_t_T Ars_i_T];
% si = [Ar_t_D Ar_i_T ArL_t_D ArL_i_T Ars_t_D Ars_i_T];
% si = [Ab_t_T Ab_i_T Abs_t_T Abs_i_T];
si_names = rasterNamesTxt(si);

siG = si; RsG = o.Rs(:,si); propsG = get_props_Rs(RsG,ntrials); resp_All = propsG.all;
mRsG = calc_mean_rasters(RsG,1:10);
trials = mat2cell([1:10]',ones(size([1:10]')));
[allRsC,allmRsT,allmRsT_Raw] = get_trial_Rs(o,si,1:10);
lnsi = length(si);
%     respDT = combine_distance_time_rasters(o.Rs(:,si(1:3)),o.Rs(:,si(4:6)),ntrials);
disp('Done');

% Overlap Indices ImageSC all
avgProps = get_props_Rs(RsG,[50,100]); 
respG = avgProps.vals;
an  = 1:5; eic = 1; sp = 0; intersect_with_global = 0; only_global = 0;
allresp = []; ind = 1;
all_peakL = []; allresp_trials = []; allpeakL_trials = []; allpeakV_trials = []; all_peakV = [];
for cn = 1:length(si)
    mRsCTi = allmRsT{cn}; mRsCTi_Raw = allmRsT_Raw{cn};
    for rrm = 1:size(mRsCTi,1)
        for ccm = 1:size(mRsCTi,2)
            mtemp = mRsCTi{rrm,ccm};
            mRsCT{rrm,ccm} = mtemp(:,1:15);
        end
    end
    mRsCT = mRsCTi; mRsCT_Raw = mRsCTi_Raw;
    size_raster_2(cn) = size(mRsCT{1,1},2);
    resp = []; peak_locations = [];
    for rr = 1:size(mRsCT,1)
        respTrials = []; pLsTrials = [];  pVsTrials = []; 
        for cc = 1:size(mRsCT,2)
            this_mat = mRsCT{rr,cc}; this_mat_Raw = mRsCT_Raw{rr,cc};
            [~,peakL] = max(this_mat,[],2);
            [peakV,peakL_Raw] = max(this_mat_Raw,[],2);
%                 size_tmat(rr,cc) = size(this_mat,2);
            resp{rr,cc} = sum(this_mat,2) > 0;
            if cc == 1
                respTrials = resp{rr,cc};
            else
                respTrials = [respTrials resp{rr,cc}];
            end
            if intersect_with_global
                resp{rr,cc} = resp{rr,cc} & respG{rr,cn};
            end
            if only_global
                resp{rr,cc} = respG{rr,cn};
            end
            if sp == 1
                if cn == 1
                    respSe = respSeL{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,1};
                end
                if cn == 2
                    respSe = respSeL{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,2};
                end
                if cn == 3
                    respSe = respSeL{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,3};
                end
                if cn == 4
                    respSe = respSeA{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,1};
                end
                if cn == 5
                    respSe = respSeA{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,2};
                end
            end
            peakL(~resp{rr,cc}) = NaN; peakV(~resp{rr,cc}) = NaN; peakL_Raw(~resp{rr,cc}) = NaN;
            peak_locations{rr,cc} = peakL; peak_values{rr,cc} = peakV;
            if cc == 1
                pLsTrials = peakL;  pVsTrials = peakV;
            else
                pLsTrials = [pLsTrials peakL];  pVsTrials = [pVsTrials peakV];
            end
            if rr == 1
                txl{ind} = sprintf('C%dT%d',cn,cc);
                ind = ind + 1;
            end
        end
        allresp_trials{rr,cn} = respTrials;
        allresp_trialsC{rr,cn} = mat2cell(respTrials,[size(respTrials,1)],[ones(1,size(respTrials,2))]);
        allpeakL_trials{rr,cn} = pLsTrials; allpeakV_trials{rr,cn} = pVsTrials;
%             oc(rr,cn) = find_cells_based_on_cluster(cell2mat(resp(rr,:)));
    end
    allresp = [allresp resp]; all_peakL = [all_peakL peak_locations]; all_peakV = [all_peakV peak_values];

end
i_allresp = cell_list_op(allresp,[],'not');

allrespOR = cell_list_op(allresp,[],'or',1);
allrespAND = cell_list_op(allresp,[],'and',1);

pallrespOR = 100*exec_fun_on_cell_mat(allrespOR,'sum')./exec_fun_on_cell_mat(allrespOR,'length');
pallrespAND = 100*exec_fun_on_cell_mat(allrespAND,'sum')./exec_fun_on_cell_mat(allrespAND,'length');

[mparOR,semparOR] = findMeanAndStandardError(pallrespOR);

disp('Done');


[OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(allresp,0.5,0.05);
%     [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(all_resp_T,0.5,0.05);
 disp('Done');

%% plot average distributions of the peak locations themselves
    all_mVals = []; all_semVals = []; pL_vals_trials_all_C = []; pL_vals_trials_all_CLin = []; pL_vals_trials_all_CLinD = [];
    pL_vals_trials_all_CD = []; mean_pVs = [];
    for cn = 1:length(si)
        all_dist_vals = [];
        pL_vals_trials_all = []; pL_vals_trials_allD = []; pV_vals_trials_all = [];
        pL_vals_trials_allLin = []; pL_vals_trials_allLinD = [];
        for an = 1:5
            pLs = allpeakL_trials{an,cn}; dpLs = diff(pLs,[],2)/size_raster_2(cn);  pVs = allpeakV_trials{an,cn};
            resp = allresp_trials{an,cn};
            perc_resp{an,cn} = 100*sum(resp,1)/size(resp,1);
            mean_pVs(an,cn) = mode(pVs(:));
            pL_vals_all = [];
            minBin = 0;maxBin = 1;BinWidth = 0.2;  binEs = minBin:BinWidth:maxBin; binCs = binEs(1:(end-1)) + BinWidth/2;
            minBinD = -1;maxBin = 1;BinWidth = 1;  binEsD = minBinD:BinWidth:maxBin; binCsD = binEsD(1:(end-1)) + BinWidth/2;
            minBinFR = 0;maxBinFR = 1;BinWidthFR = 0.2;  binEsFR = minBinFR:BinWidthFR:maxBinFR; binCsFR = binEsFR(1:(end-1)) + BinWidthFR/2;
            pL_vals_trials = []; pL_vals_trialsD = []; pV_vals_trials = [];
            pL_vals_trialsLin = []; pL_vals_trialsLinD = [];
            for trN = 1:10
                pL_vals = pLs(:,trN);  pV_vals = pVs(resp(:,trN),trN);
                pL_vals = pL_vals(resp(:,trN))/size_raster_2(cn);
                pL_vals_all{trN,1} = pL_vals;
                [bar1,binEs,binVals] = histcounts(pL_vals,binEs,'Normalization','probability');
%                 [bar1,binEs,binVals] = histcounts(pL_vals,binEs);
%                 bar1 = bar1/size(pLs,1);
                bar2 = NaN(length(binCs),length(binCsFR));
                bar2 = NaN(length(binCs),1);
                for bii = 1:length(binCs)
                    indsbins = binVals == bii;
                    tpV_vals = pV_vals(indsbins);
                    [bar2FR,binEsFR,~] = histcounts(tpV_vals,binEsFR,'Normalization','probability');
%                      bar2(:,bii) = bar2FR';
                    if ~isempty(tpV_vals)
                        bar2(bii) = mode(tpV_vals);
                    end
                end

                pL_vals_trials = [pL_vals_trials;bar1]; pV_vals_trials = [pV_vals_trials;(bar2)'];
                pL_vals_trialsLin = [pL_vals_trialsLin bar1];
                if trN < 10
                    tdpLs = dpLs((resp(:,trN) & resp(:,trN+1)),trN);
                    [bar1D,binEsD] = histcounts(tdpLs,binEsD,'Normalization','probability');
%                     [bar1D,binEsD] = histcounts(tdpLs,binEsD);
%                     bar1D = bar1D/size(pLs,1);
                    pL_vals_trialsD = [pL_vals_trialsD;bar1D];
                    pL_vals_trialsLinD = [pL_vals_trialsLinD bar1D];
                end
            end
            pL_vals_trials_all(:,:,an) = pL_vals_trials;  pL_vals_trials_allD(:,:,an) = pL_vals_trialsD;
            pV_vals_trials_all(:,:,an) = pV_vals_trials;
            pL_vals_trials_allLin(an,:) = pL_vals_trialsLin;
            pL_vals_trials_allLinD(an,:) = pL_vals_trialsLinD;
%             dist_vals = pL_vals_all;
%             [distDo,allVals] = getAveragesAndAllValues(dist_vals);
%             minBin = 0;%min(allVals);
%             maxBin = 1;%max(allVals);
%             incr = 0.1;
%             tcolors = mData.colors;
%             hf = get_figure(8,[5 7 2.25 1.5]);hold on;
%             [ha,hb,~,bins,mVals,semVals] = plotAverageDistributions(dist_vals,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'pdf_or_cdf','pdf');
%             all_dist_vals = [all_dist_vals;mVals];
%             title(sprintf('%s - %d',rasterNamesTxt{si(cn)},an));
%             pause(0.1);
        end
        pL_vals_trials_all_C{cn} = pL_vals_trials_all;  pL_vals_trials_all_CD{cn} = pL_vals_trials_allD;
        pV_vals_trials_all_C{cn} = pV_vals_trials_all;  
        m_pL_vals_trials_all{cn} = mean(pL_vals_trials_all,3);  m_pL_vals_trials_allD{cn} = mean(pL_vals_trials_allD,3);
        m_pV_vals_trials_all{cn} = mean(pV_vals_trials_all,3);
        max_maps(cn) = max(m_pL_vals_trials_all{cn}(:));  max_mapsD(cn) = max(m_pL_vals_trials_allD{cn}(:));
        max_mapsV(cn) = max(m_pV_vals_trials_all{cn}(:));
        min_maps(cn) = min(m_pL_vals_trials_all{cn}(:));  min_mapsD(cn) = min(m_pL_vals_trials_allD{cn}(:));
        min_mapsV(cn) = min(m_pV_vals_trials_all{cn}(:));
        
        sem_pL_vals_trials_all{cn} = std(pL_vals_trials_all,[],3)/sqrt(5);
        pL_vals_trials_all_CLin = [pL_vals_trials_all_CLin pL_vals_trials_allLin];
        pL_vals_trials_all_CLinD = [pL_vals_trials_all_CLinD pL_vals_trials_allLinD];
%         [mVals,semVals] = findMeanAndStandardError(all_dist_vals);
%         hf = get_figure(8,[5 7 2.25 1.5]);hold on;
%         plot(bins,mVals);
%         shadedErrorBar(bins,mVals,semVals);
%         all_mVals = [all_mVals;mVals];
%         all_semVals = [all_semVals;semVals];
    end
mM = max(max_maps);  mMD = max(max_mapsD); mMV = max(max_mapsV);
mm = min(min_maps);  mmD = min(min_mapsD); mmV = min(min_mapsV);

%% big ANOVA
[within,dvn,xlabels,awithinD] = make_within_table({'Conf','Ph','Tr','pL'},[5,2,10,length(binCs)]); %these are trials here
dataT = make_between_table({pL_vals_trials_all_CLin},dvn);
%     ra = RMA(dataT,within,{0.05,{'hsd'}});
raBA = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(raBA)

%% find average over trials and then run ANOVA
pLVals3 = [];
for gni = 1:length(si)
    tpLvals = pL_vals_trials_all_C{gni};
    mtpLvals = (squeeze(mean(tpLvals,1)))';
    pLVals3 = [pLVals3 mtpLvals];
end
[within,dvn,xlabels,awithinD] = make_within_table({'Conf','Ph','pL'},[5,2,length(binCs)]); %these are trials here
dataT = make_between_table({pLVals3},dvn);
%     ra = RMA(dataT,within,{0.05,{'hsd'}});
raBA3 = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(raBA3)

%% run post hoc ANOVA for each configuration
clc
alpha = 0.05/5;
for confi = 1:5
    redF = [1]; redV = {[confi]};
    [dataTR3,withinR] = reduce_within_between(dataT,within,redF,redV);
    raBAR3{confi} = RMA(dataTR3,withinR,{alpha,{'hsd'}});
%     raR.ranova
    print_for_manuscript(raBAR3{confi})
end

%% C2 compare with C7 and C3 and similarly others

redF = [1]; redV = {[1 5]};
[dataTR3b,withinR] = reduce_within_between(dataT,within,redF,redV);
raBAR3b = RMA(dataTR3b,withinR,{alpha,{'hsd'}});
%     raR.ranova
print_for_manuscript(raBAR3b)

%%
clc
alpha = (0.05/5)/2;
alpha = 0.05;
for confi = 1:5
    for phii = 1:2
        redF = [1,2]; redV = {[confi],[phii]};
        [dataTR31,withinR] = reduce_within_between(dataT,within,redF,redV);
        raBAR31{confi,phii} = RMA(dataTR31,withinR,{alpha,{'hsd'}});
    %     raR.ranova
        print_for_manuscript(raBAR31{confi,phii})
    end
end

%% figure bar graphs percent cells trial-wise dists

magfac = mData.magfac;
ff = makeFigureRowsCols(107,[3 5 6.9 1.3],'RowsCols',[1 10],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.35],...
    'widthHeightAdjustment',[10 -530]);
MY = 0.46; ysp = 0.0275; mY = 0; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.28*magfac; widths = [ones(1,10)*0.55]*magfac; gap = 0.115*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = repmat(mData.dcolors(1:5),1,1); 
confi = [1 1 2 2 3 3 4 4 5 5];
phi = [1 2 1 2 1 2 1 2 1 2];
for gni = 1:10
%     tra = raBART{confi(gni)};
    tra = raBAR31{confi(gni),phi(gni)};
    axes(ff.h_axes(1,gni));

    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,tra,{'pL','hsd'},[1.5 1 1]);
    xdata = make_xdata([5],[1 3]);
%     shadedErrorBar(xdata-1,mVar,semVar,{'color',tcolors{gni}});
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
    tcolors1 = tcolors;
if tra.ranova{3,tra.selected_pval_col} > alpha
    combs = [];
end
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors1,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,'capsize',1,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.25 xdata(end)+0.75],[mY MY]); format_axes(gca); xinds = [1 10]; xticks = xdata; 
if gni > 1
    set(gca,'YTick',[]);
end
if ~mod(gni,2)
    make_bars_hollow(hbs);
end
% xticklabels = cellstr(num2str((0:10)')); xticklabels = [xticklabels;xticklabels]; xticklabels = xticklabels(xinds);
xticklabels = {'L1','L2','31','L4','L5'};
set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
if gni == 1
    ylabel('Probability');
end

if gni == 5
     hxl = xlabel('Location of peak firing (binned) as a percentage of total phase length'); changePosition(hxl,[0 -0.12 0]);
end

if ismember(gni,[1 3 5 7 9])
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'AOn'},{[0.001 0.0051]});
else
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'AOff'},{[0.001 0.0051]});
end

box off;
format_axes(gca);

end

confnames = {'C2','C3','C4','C5','C7'};
for ii = 1:5
    titletxts{ii} = sprintf('%s',confnames{ii});
end
set_sub_graph_text(ff,2,titletxts,[0 0.5 0 0],[0.02 0.07 0 0]);


save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);
disp('Done')

%% figure bar graphs percent cells trial-wise dists

magfac = mData.magfac;
ff = makeFigureRowsCols(107,[3 5 6.9 1.3],'RowsCols',[1 5],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.35],...
    'widthHeightAdjustment',[10 -530]);
MY = 0.6; ysp = 0.0375; mY = 0; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.28*magfac; widths = [ones(1,5)*1.245]*magfac; gap = 0.0859*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = repmat(mData.dcolors(1:5),1,2); 
for gni = 1:5
    tra = raBAR3{gni};
axes(ff.h_axes(1,gni));

[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,tra,{'Ph_by_pL','hsd'},[1.5 1 1]);
xdata = make_xdata([5 5],[1 3]);

[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,'capsize',1,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.25 xdata(end)+0.75],[mY MY]); format_axes(gca); xinds = [1 10]; xticks = xdata; 
if gni > 1
    set(gca,'YTick',[]);
end
make_bars_hollow(hbs(6:end));
% xticklabels = cellstr(num2str((0:10)')); xticklabels = [xticklabels;xticklabels]; xticklabels = xticklabels(xinds);
xticklabels = {'L1','L2','L3','L4','L5'};
set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
% if mod(gni,2) == 1
%     titletxt = 'AOn';
% else
%     titletxt = 'AOff';
% end
% ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0.01 -0.2 0.1 0]);set(ht,'FontWeight','NOrmal');
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'AOn','AOff'},{[0.001 0.0051]});
if gni == 3
    hxl = xlabel('Location of peak firing (binned) as a percentage of total phase length'); changePosition(hxl,[0 -0.15 0]);
end

if gni == 1
    ylabel('Cell Probability Density');
end
titletxt = (sprintf('%s',rasterNamesTxt{si(gni)}));
% ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.07 0 0]);set(ht,'FontWeight','NOrmal');
box off;
format_axes(gca);
end

confnames = {'C2','C3','C4','C5','C7'};
for ii = 1:5
    titletxts{ii} = sprintf('%s',confnames{ii});
end
set_sub_graph_text(ff,1,titletxts,[0 0.5 0 0],[0.02 0.07 0 0]);


save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);
disp('Done')

%% figure bar graphs PLs C2

magfac = mData.magfac;
ff = makeFigureRowsCols(107,[10 5 1.25 1.25],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.3],...
    'widthHeightAdjustment',[10 -470]);
MY = 0.2; ysp = 0.075; mY = 0; titletxt = 'C2'; ylabeltxt = {'Cells (%)'}; % for all cells (vals) MY = 80
stp = 0.28*magfac; widths = ([0.75 0.65 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = repmat(mData.colors(1:5),1,6);
axes(ff.h_axes(1,1));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,raBAR3{1},{'Ph','hsd'},[1.5 1 1]);
    xdata = make_xdata([2],[1 1.5]);   
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca); xticks = xdata; 
xticklabels = {'AOn','AOff'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30);
% make_bars_hollow(hbs(6:end));
% [~,hyl] = put_axes_labels(gca,{'',[]},{ylabeltxt,[]}); set(hyl,'FontWeight','bold');
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'C3','C4','C5','C3','C4','C5'},{[0 0.03]});
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,6,{'Air','No-Air'},{[-0.1 -0.012]});
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[-0.01 -0.05 0.1 0]);set(ht,'FontWeight','NOrmal');
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'Pooled'},{[0 0]});
ylabel(ylabeltxt);
format_axes(gca);


save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);


%% discarded figure bar graphs PLs C7

magfac = mData.magfac;
ff = makeFigureRowsCols(107,[10 5 1.65 1.25],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.3],...
    'widthHeightAdjustment',[10 -470]);
MY = 0.8; ysp = 0.075; mY = 0; titletxt = 'C7'; ylabeltxt = {'Cells (%)'}; % for all cells (vals) MY = 80
stp = 0.28*magfac; widths = ([1.4 0.65 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = repmat(mData.colors(1:5),1,6);
axes(ff.h_axes(1,1));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,raBAR3{1},{'Ph_by_pL','hsd'},[1.5 1 1]);
    xdata = make_xdata([5 5],[1 1.5]);   
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'0.1','0.3','0.5','0.7','0.9'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30);
make_bars_hollow(hbs(6:end));
% [~,hyl] = put_axes_labels(gca,{'',[]},{ylabeltxt,[]}); set(hyl,'FontWeight','bold');
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'C3','C4','C5','C3','C4','C5'},{[0 0.03]});
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,6,{'Air','No-Air'},{[-0.1 -0.012]});
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.05 0 0]);set(ht,'FontWeight','NOrmal');
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'AOn','AOff'},{[0 0]});
ylabel(ylabeltxt);
format_axes(gca);


save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%% discarded figure bar graphs PLs C3 C4 C5

magfac = mData.magfac;
ff = makeFigureRowsCols(107,[10 5 1.25 1.25],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.3],...
    'widthHeightAdjustment',[10 -470]);
MY = 0.5; ysp = 0.075; mY = 0; titletxt = 'C5-AOn/AOff'; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.28*magfac; widths = ([0.75 0.65 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = repmat(mData.colors(1:5),1,6);
axes(ff.h_axes(1,1));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,raBAR3{2},{'pL','hsd'},[1.5 1 1]);
    xdata = make_xdata([5],[1 1.5]);   
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'0.1','0.3','0.5','0.7','0.9'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30);
make_bars_hollow(hbs(6:end));
% [~,hyl] = put_axes_labels(gca,{'',[]},{ylabeltxt,[]}); set(hyl,'FontWeight','bold');
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'C3','C4','C5','C3','C4','C5'},{[0 0.03]});
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,6,{'Air','No-Air'},{[-0.1 -0.012]});
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[-0.01 -0.05 0.1 0]);set(ht,'FontWeight','NOrmal');
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'Pooled'},{[0 0]});
ylabel(ylabeltxt);
format_axes(gca);


save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%% Find difference in peak locations trial wise 1st order then 2nd order so on and so forth
    all_mVals = []; all_semVals = []; pL_vals_trials_all_C = []; pL_vals_trials_all_CLin = []; pL_vals_trials_all_CLinD = [];
    pL_vals_trials_all_CD = []; mean_pVs = [];
    for cn = 1:length(si)
        all_dist_vals = [];
        pL_vals_trials_all = []; pL_vals_trials_allD = []; pV_vals_trials_all = [];
        pL_vals_trials_allLin = []; pL_vals_trials_allLinD = [];
        combsTR = nchoosek(1:10,2);
        selCombsTR = combsTR (diff(combsTR,[],2) == 1,:);
        
        for an = 1:5
            pLs = allpeakL_trials{an,cn};% dpLs = diff(pLs,[],2)/size_raster_2(cn);  pVs = allpeakV_trials{an,cn};
            resp = allresp_trials{an,cn};
            minBinD = -1;maxBin = 1;BinWidth = 0.5;  binEsD = minBinD:BinWidth:maxBin; binCsD = binEsD(1:(end-1)) + BinWidth/2;
            pL_vals_trialsD = []; pL_vals_trialsLinD = [];
            for trN = 1:size(selCombsTR,1)
                tr1 = selCombsTR(trN,1);
                tr2 = selCombsTR(trN,2);
                respdTR = resp(:,tr1) & resp(:,tr2);
                tdpLs = (pLs(respdTR,tr2) - pLs(respdTR,tr1))/size_raster_2(cn);
                [bar1D,binEsD] = histcounts(tdpLs,binEsD,'Normalization','probability');
%                     [bar1D,binEsD] = histcounts(tdpLs,binEsD);
%                     bar1D = bar1D/size(pLs,1);
                pL_vals_trialsD = [pL_vals_trialsD;bar1D];
                pL_vals_trialsLinD = [pL_vals_trialsLinD bar1D];
            end
            pL_vals_trials_allD(:,:,an) = pL_vals_trialsD;
            pL_vals_trials_allLinD(an,:) = pL_vals_trialsLinD;
        end
        pL_vals_trials_all_CD{cn} = pL_vals_trials_allD;
        max_mapsD(cn) = max(m_pL_vals_trials_allD{cn}(:));
        min_mapsD(cn) = min(m_pL_vals_trials_allD{cn}(:));
        pL_vals_trials_all_CLinD = [pL_vals_trials_all_CLinD pL_vals_trials_allLinD];

    end
mMD = max(max_mapsD); mMV = max(max_mapsV);
mmD = min(min_mapsD); mmV = min(min_mapsV);

% big ANOVA Diff
[within,dvn,xlabels,awithinD] = make_within_table({'Conf','Ph','Tr','pL'},[5,2,size(selCombsTR,1),length(binCsD)]); %these are trials here
dataT = make_between_table({pL_vals_trials_all_CLinD},dvn);
%     ra = RMA(dataT,within,{0.05,{'hsd'}});
raBAD = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(raBAD)

%%
leng = 2 * 9 * length(binCsD);

st = 1:leng:size(pL_vals_trials_all_CLinD,2);
et = leng:leng:size(pL_vals_trials_all_CLinD,2);
pL_valsD = [];
for ii = 1:length(st)
    pL_valsD(:,:,ii) = pL_vals_trials_all_CLinD(:,st(ii):et(ii));
end
apL_valsD = mean(pL_valsD,3);

[within,dvn,xlabels,awithinD] = make_within_table({'Ph','Tr','pL'},[2,9,length(binCsD)]); %these are trials here
dataT = make_between_table({apL_valsD},dvn);
%     ra = RMA(dataT,within,{0.05,{'hsd'}});
raBAD3 = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(raBAD3)

%%
clc
alpha = 0.05/2;
for phi = 1:2
    redF = [1]; redV = {[phi]};
    [dataTR3D,withinR] = reduce_within_between(dataT,within,redF,redV);
    raBAR3D{phi} = RMA(dataTR3D,withinR,{alpha,{'hsd'}});
%     raR.ranova
    print_for_manuscript(raBAR3D{phi})
end
print_for_manuscript(raBAR3D{2})

%%
clc
alpha = (0.05/2)/10;
for phi = 1:9
    redF = [1 2]; redV = {[1],[phi]};
    [dataTR3D1,withinR] = reduce_within_between(dataT,within,redF,redV);
    raBAR3D1{phi} = RMA(dataTR3D1,withinR,{alpha,{'hsd'}});
%     raR.ranova
    print_for_manuscript(raBAR3D1{phi})
end

%% figure bar graphs dPLs Air ON

magfac = mData.magfac;
ff = makeFigureRowsCols(107,[3 5 6.7 1.25],'RowsCols',[1 9],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.3],...
    'widthHeightAdjustment',[10 -470]);
MY = 0.9; ysp = 0.075; mY = 0; titletxt = 'Air-On'; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.28*magfac; widths = ([ones(1,10)*0.65])*magfac; gap = 0.067*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = repmat(mData.colors(1:5),1,6);
for ii = 1:9
axes(ff.h_axes(1,ii));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,raBAR3D1{ii},{'pL','hsd'},[1.5 1 1]);
    xdata = make_xdata([5],[1 1.5]);   
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'-0.8','-0.4','0','0.4','0.8'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30);
% make_bars_hollow(hbs(6:end));
% [~,hyl] = put_axes_labels(gca,{'',[]},{ylabeltxt,[]}); set(hyl,'FontWeight','bold');
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'C3','C4','C5','C3','C4','C5'},{[0 0.03]});
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,6,{'Air','No-Air'},{[-0.1 -0.012]});
if ii == 1
    titletxt = sprintf('%s T%d-T%d',titletxt,ii,ii+1);
else
    titletxt = sprintf('T%d-T%d',ii,ii+1);
end
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.05 0 0]);set(ht,'FontWeight','NOrmal');
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'Pooled'},{[0 0]});
ylabel(ylabeltxt);
format_axes(gca);
end

save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%% figure bar graphs dPLs Air OFF

magfac = mData.magfac;
ff = makeFigureRowsCols(107,[10 5 1.65 1.25],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.3],...
    'widthHeightAdjustment',[10 -470]);
MY = 0.9; ysp = 0.075; mY = 0; titletxt = 'Air-Off'; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.28*magfac; widths = ([0.9 0.65 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = repmat(mData.colors(1:5),1,6);
axes(ff.h_axes(1,1));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,raBAR3D{2},{'pL','hsd'},[1.5 1 1]);
    xdata = make_xdata([5],[1 1.5]);   
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'-0.8','-0.4','0','0.4','0.8'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30);
% make_bars_hollow(hbs(6:end));
% [~,hyl] = put_axes_labels(gca,{'',[]},{ylabeltxt,[]}); set(hyl,'FontWeight','bold');
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'C3','C4','C5','C3','C4','C5'},{[0 0.03]});
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,6,{'Air','No-Air'},{[-0.1 -0.012]});
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.05 0 0]);set(ht,'FontWeight','NOrmal');
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'Pooled'},{[0 0]});
ylabel(ylabeltxt);
format_axes(gca);


save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);


%%
clc
alpha = 0.05/5;
for confi = 1:5
    redF = [1]; redV = {[confi]};
    [dataTR,withinR] = reduce_within_between(dataT,within,redF,redV);
    raBAR{confi} = RMA(dataTR,withinR,{alpha,{''}});
%     raR.ranova
    print_for_manuscript(raBAR{confi})
end
%%
clc
alpha = (0.05/5)/2;
for confi = 1:5
    for phi = 1:2
        redF = [1,2]; redV = {[confi],[phi]};
        [dataTR,withinR] = reduce_within_between(dataT,within,redF,redV);
        raBAR1{confi,phi} = RMA(dataTR,withinR,{alpha,{'hsd'}});
    %     raR.ranova
        print_for_manuscript(raBAR1{confi,phi})
    end
end


%% figure bar graphs percent cells trial-wise dists

magfac = mData.magfac;
ff = makeFigureRowsCols(107,[3 5 6.9 1.3],'RowsCols',[1 10],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.35],...
    'widthHeightAdjustment',[10 -530]);
MY = 0.8; ysp = 0.00775; mY = 0; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.28*magfac; widths = [ones(1,10)*0.55]*magfac; gap = 0.115*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = repmat(mData.colors(1:10),1,1); 
confi = [1 1 2 2 3 3 4 4 5 5];
phi = [1 2 1 2 1 2 1 2 1 2];
for gni = 1:10
%     tra = raBART{confi(gni)};
    tra = raBAR1{confi(gni),phi(gni)};
    axes(ff.h_axes(1,gni));

    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,tra,{'pL','hsd'},[1.5 1 1]);
    xdata = make_xdata([length(binCsD)],[1 3]);
%     shadedErrorBar(xdata-1,mVar,semVar,{'color',tcolors{gni}});
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
    tcolors1 = tcolors;
if tra.ranova{5,tra.selected_pval_col} > alpha
    combs = [];
end
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors1,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,'capsize',1,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.25 xdata(end)+0.75],[mY MY]); format_axes(gca); xinds = [1 10]; xticks = xdata; 
if gni > 1
    set(gca,'YTick',[]);
end
if ~mod(gni,2)
    make_bars_hollow(hbs);
end
% xticklabels = cellstr(num2str((0:10)')); xticklabels = [xticklabels;xticklabels]; xticklabels = xticklabels(xinds);
xticklabels = {'-L5','-L4','-L3','-L2','-L1','L1','L2','31','L4','L5'};
set(gca,'xtick',xdata(1:3:end),'xticklabels',xticklabels(1:3:end)); xtickangle(15);
if gni == 1
    ylabel('Probability');
end

if gni == 5
     hxl = xlabel('Difference in location of peak firing (binned) between adjacent trials'); changePosition(hxl,[0 -0.15 0]);
end

if ismember(gni,[1 3 5 7 9])
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,length(binCsD),{'AOn'},{[0.001 0.0051]});
else
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,length(binCsD),{'AOff'},{[0.001 0.0051]});
end

box off;
format_axes(gca);   

end

confnames = {'C2','C3','C4','C5','C7'};
for ii = 1:5
    titletxts{ii} = sprintf('%s',confnames{ii});
end
set_sub_graph_text(ff,2,titletxts,[0 0.5 0 0],[0.02 0.07 0 0]);


save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);
disp('Done')


%% run ANOVArms PLs
clc
all_ras = [];
[within,dvn,xlabels,awithinD] = make_within_table({'Tr','pL'},[10,length(binCs)]); %these are trials here
for gn = 1:10
    tvals = pL_vals_trials_all_C{gn};
    tvalsL = [];
    for an = 1:5
        tvalsL = [tvalsL;reshape(tvals(:,:,an)',1,numel(tvals(:,:,an)))];
    end
    dataT = make_between_table({tvalsL},dvn);
%     ra = RMA(dataT,within,{0.05,{'hsd'}});
    ra = RMA(dataT,within,{0.005,{'hsd'}});
    all_ras{gn} = ra;
%     ra = RMA(dataT,within,{0.05,''});
% %     ra.ranova
    all_pvals(:,gn) = ra.ranova{[3 5 7],ra.selected_pval_col};
    all_pvalsFs(:,gn) = ra.ranova{[3 5 7],4};
    print_for_manuscript(ra)
end

%% figure heat maps of probabilities trial wise
magfac = mData.magfac;
ff = makeFigureRowsCols(107,[2 3 6.9 1.75],'RowsCols',[2 10],'spaceRowsCols',[0.04 0.13],'rightUpShifts',[0.03 0.2],...
    'widthHeightAdjustment',[-100 -180]);
MY = 70; ysp = 5; mY = 0; titletxt = 'Activated Cells'; ylabeltxt = {'Cells (%)'}; % for all cells (vals) MY = 80
stp = 0.35*magfac; widths = ([ones(1,10)*0.55]-0.05)*magfac; gap = 0.15*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
% changePosition(ff.h_axes(1,:),[0 0.1 0 -0.05]);
% changePosition(ff.h_axes(2,:),[0 0.3 0 -0.05]);
% changePosition(ff.h_axes(3,:),[0 0 0 0.15]);

[within,dvn,xlabels,awithinD] = make_within_table({'pL'},[length(binCs)]);
for gn = 1:10
%     axes(ff.h_axes(2,gn));
%     imagesc(1:length(binCs),1:10,m_pL_vals_trials_all{gn},[mm mM]);
    
    tvals = pL_vals_trials_all_C{gn};
    mtvals = (squeeze(mean(tvals,1)))';
%     mtvals = (squeeze(tvals(10,:,:)))';
    dataT = make_between_table({mtvals},dvn);
%     ra = RMA(dataT,within,{0.05,{'bonferroni','hsd'}});
%     ra = RMA(dataT,within,{0.05,''});
%     ra.ranova
%     print_for_manuscript(ra)
%     ra = all_ras{gn};
%     pval(gn) = ra.ranova{5,ra.selected_pval_col};
%     pvalph(gn) = ra.MC.hsd.pL{12,5};
%     pvalphB(gn) = ra.MC.bonferroni.pL{12,5};

    axes(ff.h_axes(1,gn));
    imagesc(binCs,1:10,tvals(:,:,3),[mm mM]);
%     titletxt = (sprintf('%s - %s',rasterNamesTxt{si(gn)},getNumberOfAsterisks(pval(gn))));
    titletxt = (sprintf('%s',rasterNamesTxt{si(gn)}));
    ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.1 0 0]);set(ht,'FontWeight','NOrmal');
    set(gca,'Xticklabels',[],'YTick',[1 5 10],'YDir','Normal');
    if gn > 1
        set(gca,'Yticklabels',[]);
    else
        ylabel('Trials');
    end
    format_axes(gca);
    axes(ff.h_axes(2,gn));
    imagesc(binCs,1:10,m_pL_vals_trials_all{gn},[mm mM]);
    set(gca,'YTick',[1 5 10],'YDir','Normal','XTick',binCs);xtickangle(30);
    if gn > 1
        set(gca,'Yticklabels',[]);
    else
        ylabel('Trials');
    end
    if gn == 5
        xlabel('Location of Peak Firing as a Fraction of total Phase Length');
        hxl = xlabel('Location of peak firing (binned) as a percentage of total phase length'); changePosition(hxl,[0 -0.15 0]);
    end
    xticklabels = {'L1','L2','L3','L4','L5'};
    set(gca,'xtick',binCs,'xticklabels',xticklabels); xtickangle(30);

% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'AOn','AOff'},{[0.001 0.0051]});


    format_axes(gca);
    if gn == 10
        hc = putColorBar(ff.h_axes(1,gn),[0.0 0.03 0 -0.05],[mm mM],6,'eastoutside',[0.07 0.07 0.1 0.1]);
%         hc = putColorBar(ff.h_axes(2,gn),[0.0 0.03 0 -0.05],[mm mM],6,'eastoutside',[0.07 0.07 0.1 0.1]);
    end
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('peak_firingdists.pdf'),600);

%% run ANOVArms PLs-Diff
clc
all_ras = [];
[within,dvn,xlabels,awithinD] = make_within_table({'Tr','pL'},[9,length(binCsD)]); %these are trials here
for gn = 1:10
    tvals = pL_vals_trials_all_CD{gn};
    tvalsL = [];
    for an = 1:5
        tvalsL = [tvalsL;reshape(tvals(:,:,an)',1,numel(tvals(:,:,an)))];
    end
    dataT = make_between_table({tvalsL},dvn);
%     ra = RMA(dataT,within,{0.05,{'hsd'}});
    ra = RMA(dataT,within,{0.05,{'hsd'}});
    all_ras{gn} = ra;
%     ra = RMA(dataT,within,{0.05,''});
% %     ra.ranova
    all_pvals(:,gn) = ra.ranova{[3 5 7],ra.selected_pval_col};
    all_pvalsFs(:,gn) = ra.ranova{[3 5 7],4};
    print_for_manuscript(ra)
end

%% figure diff PLs
magfac = mData.magfac;
ff = makeFigureRowsCols(107,[2 3 6.9 3],'RowsCols',[3 10],'spaceRowsCols',[0.09 0.13],'rightUpShifts',[0.03 0.13],...
    'widthHeightAdjustment',[-100 -150]);
MY = 70; ysp = 5; mY = 0; titletxt = 'Activated Cells'; ylabeltxt = {'Cells (%)'}; % for all cells (vals) MY = 80
stp = 0.35*magfac; widths = ([ones(1,10)*0.55]-0.05)*magfac; gap = 0.15*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
changePosition(ff.h_axes(1,:),[0 0.1 0 -0.05]);
changePosition(ff.h_axes(2,:),[0 0.3 0 -0.05]);
changePosition(ff.h_axes(3,:),[0 0 0 0.15]);

[within,dvn,xlabels,awithinD] = make_within_table({'pL'},[length(binCsD)]);
for gn = 1:10
%     axes(ff.h_axes(2,gn));
%     imagesc(1:length(binCsD),1:10,m_pL_vals_trials_all{gn},[mm mM]);
%     
    tvals = pL_vals_trials_all_CD{gn};
    mtvals = (squeeze(mean(tvals,1)))';
%     mtvals = (squeeze(tvals(10,:,:)))';
    dataT = make_between_table({mtvals},dvn);

    ra = all_ras{gn};
    pval(gn) = ra.ranova{5,ra.selected_pval_col};
%     pvalph(gn) = ra.MC.hsd.pL{12,5};
%     pvalphB(gn) = ra.MC.bonferroni.pL{12,5};

    axes(ff.h_axes(1,gn));
    imagesc(binCsD,1:9,tvals(:,:,3),[mmD mMD]);
    titletxt = (sprintf('%s - %s',rasterNamesTxt{si(gn)},getNumberOfAsterisks(pval(gn))));
    titletxt = (sprintf('%s',rasterNamesTxt{si(gn)}));
    ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.1 0 0]);set(ht,'FontWeight','NOrmal');
    set(gca,'Xticklabels',[],'YTick',[1 5 9],'YDir','Normal');
    if gn > 1
        set(gca,'Yticklabels',[]);
    else
        ylabel('Trials');
    end
    
    axes(ff.h_axes(2,gn));
    imagesc(binCsD,1:9,m_pL_vals_trials_allD{gn},[mmD mMD]);
    set(gca,'YTick',[1 5 9],'YDir','Normal');
    if gn > 1
        set(gca,'Yticklabels',[]);
    else
        ylabel('Trials');
    end
    axes(ff.h_axes(3,gn));
    ysp = 0.03;
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'pL','hsd'},[1.5 1 1]);
    xdata = make_xdata([length(binCsD)],[1 1.5]); 
    tcolors = mData.dcolors;
    if pval(gn) >= 0.05
        combs = [];
    end 
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    ylim([0 0.75]); xlim([0.5 xdata(end)+0.75]);
%     set(gca,'XTick',[1 3 5],'XTickLabels',{'0.1','0.3','0.5','0.7','0.9'});
    set(gca,'XTick',[1 3 5],'XTickLabels',{'-0.9','0','0.9'});
    if gn > 1
        set(gca,'Yticklabels',[]);
    else
        ylabel('Cells (%)');
    end
    if gn == 6
        xlabel('Trial-wise Differernce in Location of Peak Firing as Fraction of total Phase Length');
    end
    if pval(gn) < 0.05
%     titletxt = (sprintf('%s',getNumberOfAsterisks(pval(gn))));
%     ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.01 0 0]);set(ht,'FontWeight','NOrmal','FontSize',5);
%     titletxt = (sprintf('F(4,16) = %.3f',all_pvalsFs(2,gn)));
%     ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.031 0.2 0]);set(ht,'FontWeight','NOrmal','FontSize',5);
    titletxt = (sprintf('%s (%.3f)',getNumberOfAsterisks(pval(gn)),all_pvalsFs(2,gn)));
    ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.031 0.2 0]);set(ht,'FontWeight','NOrmal','FontSize',6.5);
    end
    if gn == 10
        hc = putColorBar(ff.h_axes(1,gn),[0.0 0.03 0 -0.05],[mmD mMD],6,'eastoutside',[0.07 0.07 0.1 0.1]);
%         hc = putColorBar(ff.h_axes(2,gn),[0.0 0.03 0 -0.05],[mm mM],6,'eastoutside',[0.07 0.07 0.1 0.1]);
    end
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('peak_firingdists.pdf'),600);


 %% heatmap conj
hf = get_figure(6,[8 3 3.25 3.25]);
%     allresp_sp = [allresp resp_speed(:,4)];
allresp_sp = allresp(1:5,:) ;
[mOI,semOI] = heatmap_conj_comp(gca,allresp_sp,1,{si,rasterNamesTxt});
save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_mean_tu_spatial.pdf'),600);

figdim = 1;
hf = get_figure(7,[5 2 figdim+0.1 figdim]);
minI = min(semOI(:)); maxI = max(semOI(:));
sz = size(mOI,1);
oM = ones(size(mOI));mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
imAlpha(mask1 == 1) = 0;
im1 = imagesc(semOI,[minI,maxI]);    im1.AlphaData = imAlpha;
set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
format_axes(gca);
set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',[],'yticklabels',[],'Ydir','reverse'); xtickangle(75);
changePosition(gca,[0.0 0 -0.05 0]);
set(gca,'Ydir','normal');
box off
colormap jet
save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_sem.pdf'),600);

%% trial by trial recruitment and position in sequence
    an = 1;
    cn = 1;
    resps = [];
    for trN = 1:100
        resps(:,trN) = allresp{an,trN};
        pks(:,trN) = all_peakL{an,trN};
    end
    sresps = sum(resps,2);
    
    psresps = 100*sresps/size(resps,2);
    
%     [clus_inds,clus_centers] = kmeans(psresps,5);
    
    [sr,sri] = sort(sresps);
    sort_resps = resps(sri,:);
%     [uvs,yvs] = trial_collections(resps);
    %% align by peaks then view
    s_cell_list = [];
    for ii = 1:100
    peak1 = pks(:,ii);
    [sr,sri] = sort(peak1);
    s_cell_list(:,ii) = sri;
    sort_resps = resps(sri,:);
    figure(1000);clf;imagesc(sort_resps)
    pause(0.01);
    end
    figure(1000);clf;imagesc(s_cell_list);colorbar
    %%
    sort_resps = resps(clus_inds == 3,:);
    figure(1000);clf;imagesc(sort_resps)
    %% finding all cells that were recruited in different configurations or their phases
    st = 1:10:(length(si)*10);
    et = 10:10:(length(si)*10);
    for ii = 1:length(si)
        t_resp = resps(:,st(ii):et(ii));
        c_ids = [];
        for jj = 1:10
            tt_resp = t_resp(:,jj);
            c_ids = [c_ids;find(tt_resp)];
        end
        ucids(ii) = length(unique(c_ids));
    end
    
    %%
    pLs = all_peakL{an,cn};
    resp_cells_pLs = pLs(resps);
    
    presp = find_percent(allresp);
    varC = presp;
    [within,dvn,xlabels,awithinD] = make_within_table({'Cond','TI','Tr'},[2,2,10]);
    dataT = make_between_table({varC},dvn);
    ra = RMA(dataT,within,{0.05,{'bonferroni','hsd'}});
    ra.ranova
    print_for_manuscript(ra)
    %% Heatmap of active cells
    st = 1:10:(length(si)*10);
    et = 10:10:(length(si)*10);
    for ii = 1:length(si)
        t_resp = allresp(:,st(ii):et(ii));
        allresp_OR(:,ii) = cell_list_op(t_resp,[],'or',1);
    end
%     si = [
    hf = get_figure(6,[8 3 3.25 3.25]);
    [mOI,semOI] = heatmap_conj_comp(gca,allresp_OR,1,{si,rasterNamesTxt,1});
    
    %%
    all_pos_seq = [];
    for an = 1:5
        for cn = 1:length(si)
            pLs = allpeakL_trials{an,cn};
            dpLs = diff(pLs,[],2)/size_raster_2(cn);
            resp = allresp_trials{an,cn};
            resp = sum(resp,2)>2;
            pLs = pLs(resp,:);
            tr_seq = [];
            for trn = 1:10
                [~,inds] = sort(pLs(:,trn));
                [~,inds1] = sort(inds);
                nan_temp = double(isnan(pLs(:,trn)));
                nan_temp(nan_temp == 1) = NaN; nan_temp(nan_temp == 0) = 1;
                tr_seq(:,trn) = inds1 .* nan_temp;
            end
            all_pos_seq{an,cn} = tr_seq;
            seq_shift = diff(tr_seq)/size(tr_seq,1);
            m_seq_shift = nanmean(seq_shift,2); % this is mean over trials;
            f_seq_shift_trials_LZ = sum(seq_shift<0,2)/10;
            mm_seq_shift(an,cn) = nanmean(m_seq_shift); % this is mean over cells
            f_seq_shift_LZ(an,cn) = sum(m_seq_shift<-0.2)/size(tr_seq,1);
            
            mdpLs{an,cn} = nanmean(dpLs,2);
        end
    end
    
    %% run stats on average shifts of peak locations
    inds = [3 4 5 6 7 8];
%     inds = [1 2 3 4];% inds = [1 2 5 6]; inds = [1 2 9 10];
%     inds = [1 2 9 10];
    varC = f_seq_shift_LZ(:,inds);
    varC = mm_seq_shift(:,inds);
    [within,dvn,xlabels,awithinD] = make_within_table({'Cnds','TI'},[length(inds)/2,2]);
    dataT = make_between_table({varC},dvn);
    ra = RMA(dataT,within,{0.05,{'bonferroni','hsd'}});
    ra.ranova
    print_for_manuscript(ra)
    %% for all active cells in a trial, I want to see the peak locations and average shift of peak locations for animals
    binwidths = evalin('base','binwidths');
    an = 1; cn = 1;
    mshifts = [];mshiftsB = []; mshiftsF = []; fshiftB = []; fshiftF = []; fshiftZ = [];
    all_shifts = [];
    for an = 1:5
        for cn = 1:length(si)
            pLs = allpeakL_trials{an,cn};
            resp = allresp_trials{an,cn};
            perc_resp{an,cn} = 100*sum(resp,1)/size(resp,1); m_pLs_shifts = [];
            for trN = 1:9
                t_resp = resp(:,trN) & resp(:,trN+1);
                pLs_shifts = (pLs(t_resp,trN+1) - pLs(t_resp,trN))/size_raster_2(cn);
                m_pLs_shifts(trN) = mean(pLs_shifts); % mean over cells
            end
            all_m_pLs_shifts(an,cn) = mean(m_pLs_shifts); % mean over trial pairs
            all_m_pLs_shiftsT{an,cn} = m_pLs_shifts; % mean over trial pairs
            resp = sum(resp,2)>2; % to find cells that replied in at least 3 trials
            pLs = pLs(resp,:);
            [~,inds] = sort(pLs(:,1));
        %     figure(1000);clf;imagesc(pLs(inds,:));colorbar;
%             figure(1000);clf;imagesc(pLs);colorbar;
            100*sum(~isnan(pLs))/size(pLs,1);
            prob_participation = [];
            mshift = []; 
            for ii = 1:size(pLs,1)
                tcell = pLs(ii,:);
                prob_participation(ii,1) = sum(~isnan(tcell))/length(tcell);
                ntcell = tcell(~isnan(tcell));
                mshift(ii,1) = mean(diff(ntcell)/size_raster_2(cn));
            end
            all_shifts{an,cn} = mshift;
            mshifts(an,cn) = mean(mshift); % average over cells
            mshiftsB(an,cn) = mean(mshift(mshift<0)); % average backward shift
            mshiftsF(an,cn) = mean(mshift(mshift>0)); % average forward shift
            fshiftB(an,cn) = sum(mshift<0)/length(mshift); % fraction of cells that had backward mean shift
            fshiftF(an,cn) = sum(mshift>0)/length(mshift); % fraction of cells that had forward mean shift
            fshiftZ(an,cn) = sum(mshift==0)/length(mshift); % fraction of cells that had zero shift
        end
    end
    
    %% run stats on average shifts of peak locations
    inds = [3 4 5 6 7 8];
%     inds = [1 2 3 4];
%     inds = 1:6;
    varC = [abs(mshiftsB(:,inds)) mshiftsF(:,inds)];
    varC = all_m_pLs_shifts(:,inds);
    [within,dvn,xlabels,awithinD] = make_within_table({'Cnds','TI'},[length(inds)/2,2]);
    dataT = make_between_table({varC},dvn);
    ra = RMA(dataT,within,{0.05,{'bonferroni'}});
    ra.ranova
    print_for_manuscript(ra)
    
    %%
    inds = [1 2 3 4]; inds = [1 2 9 10];
    inds = [3 4 5 6 7 8];
%     varC = [mshiftsF(:,inds) abs(mshiftsB(:,inds))];
    varC = [fshiftF(:,inds) fshiftB(:,inds) fshiftZ(:,inds)];
    varC = [fshiftF(:,inds) fshiftB(:,inds)];
    [within,dvn,xlabels,awithinD] = make_within_table({'Po','Cnds','TI'},[2,3,2]);
    dataT = make_between_table({varC},dvn);
    ra = RMA(dataT,within,{0.05,{'bonferroni','hsd'}});
    ra.ranova
    print_for_manuscript(ra)
    %%
    an = 5;
    big_mR = [];
    for ii = 1:size(mRsG,2)
        big_mR = [big_mR mRsG{an,ii}];
    end
    
    CRc = corr(big_mR);
    figure(1000);clf;
    subplot(2,1,1);imagesc(big_mR);colorbar;set(gca,'Ydir','Normal');
    subplot(2,1,2);imagesc(CRc);colorbar;set(gca,'Ydir','Normal');
%% scrambling of cellular activations

activated_cells_OR = cell_list_op(allresp,[],'or',1);

p_activated_cells_OR = find_percent(activated_cells_OR);

sci = [];
for sii = 1:length(si)
    for ani = 1:5
       tresp_trials = allresp_trials{ani,sii};
       tdists_centroids = dists_centroids{ani};
       tcentroids = centroids{ani};
       tac = activated_cells_OR{ani};
       for tri = 1:9
           ttrac1 = tresp_trials(:,tri);
           ttrac2 = tresp_trials(:,tri+1);
           conj1 = ttrac1 & tac;
           thcentroids1 = tcentroids(conj1,:);
           conj2 = ttrac2 & tac;
           thcentroids2 = tcentroids(conj2,:);
           tcombs1 = nchoosek(1:size(thcentroids1),2);
           tcombs2 = nchoosek(1:size(thcentroids2),2);
           parfor cmbi = 1:size(tcombs1,1)
               dist_cent1(cmbi) = sqrt(sum((thcentroids1(tcombs1(cmbi,1),:) - thcentroids1(tcombs1(cmbi,2),:)).^2));
           end
           parfor cmbi = 1:size(tcombs2,1)
               dist_cent2(cmbi) = sqrt(sum((thcentroids2(tcombs2(cmbi,1),:) - thcentroids2(tcombs2(cmbi,2),:)).^2));
           end
           [h,p,ks2stat] = kstest2(dist_cent1,dist_cent2);
           all_p(sii,tri,ani) = p;
           all_h(sii,tri,ani) = h;
       end
    end
end



%%
[within,dvn,xlabels,awithinD] = make_within_table({'cns'},[length(si)]);
[within,dvn,xlabels,awithinD] = make_within_table({'Br','Ph'},[2,2]);
dataT = make_between_table({mean_pVs(:,[1 2 3 4])},dvn);
ra = RMA(dataT,within,{0.05,{'bonferroni','hsd'}});
print_for_manuscript(ra)
%%
ff = makeFigureRowsCols(108,[1 1 6.9 3],'RowsCols',[2 10],...
'spaceRowsCols',[0.09 0.025],'rightUpShifts',[0.03 0.13],'widthHeightAdjustment',...
[-30 -200]);
[within,dvn,xlabels,awithinD] = make_within_table({'pL'},[length(binCsD)]);
for gn = 1:10
    tvals = pL_vals_trials_all_CD{gn};
    mtvals = (squeeze(mean(tvals,1)))';
%     mtvals = (squeeze(tvals(10,:,:)))';
    dataT = make_between_table({mtvals},dvn);
    ra = RMA(dataT,within,{0.05,{'bonferroni','hsd'}});
%     ra = RMA(dataT,within,{0.05,''});
%     ra.ranova
%     print_for_manuscript(ra)
%     ra = all_rasD{gn};
    pval(gn) = ra.ranova{3,ra.selected_pval_col};
%     pvalph(gn) = ra.MC.hsd.pL{12,5};
%     pvalphB(gn) = ra.MC.bonferroni.pL{12,5};
    axes(ff.h_axes(1,gn));
    imagesc(1:5,1:9,m_pL_vals_trials_allD{gn},[mmD mMD]);
    title(sprintf('%s - %s',rasterNamesTxt{si(gn)},getNumberOfAsterisks(pval(gn))));
%     plot(mean(mtvals));
%     shadedErrorBar(1:5,mean(mtvals),std(mtvals)/sqrt(5));
%     title(pval);
    axes(ff.h_axes(2,gn));
    ysp = 0.03;
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'pL','hsd'},[1.5 1 1]);
    xdata = make_xdata([length(binCsD)],[1 1.5]); 
    tcolors = mData.dcolors;
    if pval(gn) >= 0.05
        combs = [];
    end 
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    ylim([0 0.5]); xlim([0.5 length(binCsD)+0.75]);
    if gn > 1
        set(gca,'Yticklabels',[]);
    end
    if gn == 10
        hc = putColorBar(ff.h_axes(1,gn),[0.0 0.03 0 -0.05],[mm mM],6,'eastoutside',[0.07 0.07 0.1 0.1]);
    end
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('peak_firingdists.pdf'),600);


%%
ff = makeFigureRowsCols(108,[1 1 6.9 1.5],'RowsCols',[1 10],...
'spaceRowsCols',[0.05 0.025],'rightUpShifts',[0.03 0.13],'widthHeightAdjustment',...
[-30 -250]);
[within,dvn,xlabels,awithinD] = make_within_table({'pL'},[10]); %these are trials here
for gn = 1:10
    tvals = pL_vals_trials_all_C{gn};
    mtvals = (squeeze(mean(tvals,1)))';
    mtvals = (squeeze(tvals(:,5,:)))';
    dataT = make_between_table({mtvals},dvn);
    ra = RMA(dataT,within,{0.05,{'bonferroni','hsd'}});
%     ra = RMA(dataT,within,{0.05,''});
%     ra.ranova
    print_for_manuscript(ra)
%     pval(gn) = ra.ranova{3,ra.selected_pval_col};
%     pvalph(gn) = ra.MC.hsd.pL{12,5};
%     pvalphB(gn) = ra.MC.bonferroni.pL{12,5};
    axes(ff.h_axes(1,gn));
%     imagesc(m_pL_vals_trials_all{gn});
%     plot(mean(mtvals));
%     shadedErrorBar(1:5,mean(mtvals),std(mtvals)/sqrt(5));
%     title(pval);
    ysp = 0.03;
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'pL','hsd'},[1.5 1 1]);
    xdata = make_xdata([10],[1 1.5]); 
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    ylim([0 0.5]);
    title(sprintf('%s',rasterNamesTxt{si(gn)}));
end

%%
% ff = makeFigureRowsCols(108,[1 1 6.9 1.5],'RowsCols',[1 10],...
% 'spaceRowsCols',[0.05 0.025],'rightUpShifts',[0.03 0.13],'widthHeightAdjustment',...
% [-30 -250]);
[within,dvn,xlabels,awithinD] = make_within_table({'Tr','pL'},[10,5]); %these are trials here
all_vals = [];
for gn = 1:10
    tvals = pL_vals_trials_all_C{gn};
    tvals = pV_vals_trials_all_C{gn};
    tvalsL = [];
    for an = 1:5
        tvalsL = [tvalsL;reshape(tvals(:,:,an)',1,50)];
    end
    if gn < 5
        all_vals = [all_vals tvalsL];
    end
    dataT = make_between_table({tvalsL},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd'}});
    all_ras{gn} = ra;
end

%%
[within,dvn,xlabels,awithinD] = make_within_table({'Br','Ph','Tr','pL'},[2,2,10,5]); %these are trials here
dataT = make_between_table({all_vals},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra);


%%
all_vals = [];
for gn = [1 2 9 10]
    tvals = pL_vals_trials_all_C{gn};
    tvalsL = [];
    for an = 1:5
        tvalsL = [tvalsL;reshape(tvals(:,:,an)',1,50)];
    end
%     if gn < 5
        all_vals = [all_vals tvalsL];
%     end
end
[within,dvn,xlabels,awithinD] = make_within_table({'Br','Ph','Tr','pL'},[2,2,10,5]); %these are trials here
dataT = make_between_table({all_vals},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra);
%%
all_vals = [];
for gn = [1 3]
    tvals = pV_vals_trials_all_C{gn};
    tvalsL = [];
    for an = 1:5
        tvalsL = [tvalsL;reshape(tvals(:,:,an)',1,50)];
    end
%     if gn < 5
        all_vals = [all_vals tvalsL];
%     end
end
[within,dvn,xlabels,awithinD] = make_within_table({'Br','Tr','pL'},[2,10,5]); %these are trials here
dataT = make_between_table({all_vals},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra);
%%
[within,dvn,xlabels,awithinD] = make_within_table({'Br','Ph','Tr','pL'},[2,2,10,5]); %these are trials here
dataT = make_between_table({all_vals},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra);


%%
hf = get_figure(8,[5 5 6.9 3]);hold on;
mean_pL_vals_trials_all_CLin = mean(pL_vals_trials_all_CLin);
sem_pL_vals_trials_all_CLin = std(pL_vals_trials_all_CLin)/sqrt(5);

plot(1:500,mean_pL_vals_trials_all_CLin);
shadedErrorBar(1:500,mean_pL_vals_trials_all_CLin,sem_pL_vals_trials_all_CLin);
    %%
    hf = get_figure(8,[5 5 2.25 5]);hold on;
    xs = binEs;
    cn = 1; ys = all_mVals(cn,:); eys = all_semVals(cn,:); 
    for cn = 1:length(si)
        if cn > 1
            ys = max(ys) + all_mVals(cn,:); eys = all_semVals(cn,:);
        end
        plot(xs,ys);
        shadedErrorBar(xs,ys,eys);
    end
    hf = get_figure(8,[5 5 2.25 5]);hold on;
    plot(xs,all_mVals');
    
    %%
    for cn = 1:length(si)
        hf = get_figure(8,[5 5 5 5]);hold on;
        imagesc(binCs,1:10,m_pL_vals_trials_all{cn},[0 0.5]);colorbar;
        pause(0.1);
    end
    
    %%
    inds = [3 4 5 6 7 8];
    varC = pL_vals_trials_all_CLin(:,101:400);
    [within,dvn,xlabels,awithinD] = make_within_table({'Co','Ph','Tr','Sp'},[3,2,10,5]);
    dataT = make_between_table({varC},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd'}});
%     ra = RMA(dataT,within,{0.05,''});
    ra.ranova
    print_for_manuscript(ra)
    %%
    redF = [1,2]; redV = {[1,2],1};
    [dataTR,withinR] = reduce_within_between(dataT,within,redF,redV);
    raR = RMA(dataTR,withinR,{0.025,{'hsd'}});
    raR.ranova
    print_for_manuscript(raR)
    
    %%
    for cn = 1:length(si)
    varC = cell2mat(all_m_pLs_shiftsT(:,cn));
    [within,dvn,xlabels,awithinD] = make_within_table({'T'},[9]);
    dataT = make_between_table({varC},dvn);
    ra = RMA(dataT,within,{0.05,{'bonferroni','hsd'}});
    ra.ranova
    print_for_manuscript(ra)
    
    [mVals,semVals] = findMeanAndStandardError(varC);
    hf = get_figure(8,[5 7 2.25 1.5]);hold on;
    plot(1:9,mVals);
    shadedErrorBar(1:9,mVals,semVals);
    set(gca,'Ylim',[-0.5 0.5]);
    pause(0.3);
    end
    %%
    %% TI cond
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Sp_by_Co','hsd'},[1.5 1 1]);
    xdata = make_xdata([5 5 5],[1 1.5]);
    hf = get_figure(5,[8 7 6.9 1]);
    tcolors = repmat(mData.colors(1:5),1,3);
    MmVar = max(mVar);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',MmVar/3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'D','T'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.06 0.01 -0.05 0]); 
%     put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
%     save_pdf(hf,mData.pdf_folder,sprintf('var_TI_by_Cond.pdf'),600);

