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
[mRsG,~,mRsG_raw] = calc_mean_rasters(RsG,1:10);
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
%     for rrm = 1:size(mRsCTi,1)
%         for ccm = 1:size(mRsCTi,2)
%             mtemp = mRsCTi{rrm,ccm};
%             mRsCT{rrm,ccm} = mtemp(:,1:15);
%         end
%     end
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

%% find the location of a cell in an ensemble when all cells are aligned based on peaks
eLs = []; % store ensemble locations as a percentage of the total number of activated cells in that trial
minBin = 0;maxBin = 1;BinWidth = 0.2;  binEs = minBin:BinWidth:maxBin; binCs = binEs(1:(end-1)) + BinWidth/2;
d_eLs = []; d_eLs_Lin = []; m_d_eLs_Lin = [];
d_eLs_C = [];
for cn = 1:length(si)
    mRsCTi = allmRsT{cn}; mRsCTi_Raw = allmRsT_Raw{cn};
%     for rrm = 1:size(mRsCTi,1)
%         for ccm = 1:size(mRsCTi,2)
%             mtemp = mRsCTi{rrm,ccm};
%             mRsCT{rrm,ccm} = mtemp(:,1:15);
%         end
%     end
    mRsCT = mRsCTi; mRsCT_Raw = mRsCTi_Raw;
%     size_raster_2(cn) = size(mRsCT{1,1},2);
    an_eLs = []; m_an_eLs = [];
    an_eLs_C = []; 
    for rr = 1:size(mRsCT,1)
        td_eLs = [];
        for cc = 1:size(mRsCT,2)
            tmRsCT = mRsCT{rr,cc};
            tresp = allresp_trials{rr,cn}(:,cc);
            tpLs = allpeakL_trials{rr,cn}(:,cc);
            
            rtmRsCT = tmRsCT(tresp,:);
            rpLs = tpLs(tresp,:);
            
%             [ia,ib] = sort(rpLs); % ib has the indices of cells when they are aligned by peaks which are in rpLs
%             [iaa,ibb] = sort(ib);
%             teLs = tpLs; % temporarirly assign the same
%             
            [ia,ib] = sort(tpLs); % ib has the indices of cells when they are aligned by peaks which are in rpLs
            [iaa,ibb] = sort(ib);
            teLs = tpLs; % temporarirly assign the same
            
            perc_eL = ibb/length(ibb);
            perc_eL(~tresp) = NaN;
%             teLs(tresp) = perc_eL;
            eLs{rr,cn}(:,cc) = teLs;
            
            [bar1,binEs,binVals] = histcounts(perc_eL,binEs,'Normalization','probability');
            td_eLs = [td_eLs bar1];
%             artmRsCT = rtmRsCT(ib,:);
%             figure(100);clf;
%             imagesc(artmRsCT);set(gca,'Ydir','normal');
        end
        an_eLs = [an_eLs;td_eLs];
        m_an_eLs = [m_an_eLs;mean((reshape(td_eLs,5,10))')];
        an_eLs_C(:,:,rr) = (reshape(td_eLs,5,10))';
    end
    d_eLs{cn} = an_eLs;     d_eLs_Lin = [d_eLs_Lin an_eLs];     m_d_eLs_Lin = [m_d_eLs_Lin m_an_eLs];
    d_eLs_C{cn} = an_eLs_C;
    mmm(cn) = min(an_eLs_C(:)); MMM(cn) = max(an_eLs_C(:));
end
mm = min(mmm); mM = max(MMM);
%% Figure 4C figure heat maps of probabilities trial wise
magfac = mData.magfac;
ff = makeFigureRowsCols(107,[2 3 6.9 1.5],'RowsCols',[2 10],'spaceRowsCols',[0.04 0.13],'rightUpShifts',[0.03 0.23],...
    'widthHeightAdjustment',[-100 -180]);
MY = 70; ysp = 5; mY = 0; titletxt = 'Activated Cells'; ylabeltxt = {'Cells (%)'}; % for all cells (vals) MY = 80
stp = 0.35*magfac; widths = ([ones(1,10)*0.55]-0.05)*magfac; gap = 0.15*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
[within,dvn,xlabels,awithinD] = make_within_table({'pL'},[length(binCs)]);
for gn = 1:10

    tvals = d_eLs_C{gn};
    mtvals = mean(tvals,3);
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
    imagesc(binCs,1:10,mtvals,[mm mM]);
    set(gca,'YTick',[1 5 10],'YDir','Normal','XTick',binCs);xtickangle(30);
    if gn > 1
        set(gca,'Yticklabels',[]);
    else
        ylabel('Trials');
    end
    if gn == 5
%         xlabel('Location of Peak Firing as a Fraction of total Phase Length');
        hxl = xlabel('Location of peak firing (binned) as a percentage of total phase length'); changePosition(hxl,[0 0.51 0]);
    end
    xticklabels = {'r1','r2','r3','r4','r5'};
    set(gca,'xtick',binCs,'xticklabels',xticklabels); xtickangle(45);

% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'AOn','AOff'},{[0.001 0.0051]});


    format_axes(gca);
    if gn == 10
        hc = putColorBar(ff.h_axes(1,gn),[0.0 0.03 0 -0.05],[mm mM],6,'eastoutside',[0.07 0.07 0.1 0.1]);
%         hc = putColorBar(ff.h_axes(2,gn),[0.0 0.03 0 -0.05],[mm mM],6,'eastoutside',[0.07 0.07 0.1 0.1]);
    end
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('peak_firingdists.pdf'),600);

%% big ANOVA
[within,dvn,xlabels,awithinD] = make_within_table({'Conf','Ph','Tr','pL'},[5,2,10,length(binCs)]); %these are trials here
dataT = make_between_table({d_eLs_Lin},dvn);
%     ra = RMA(dataT,within,{0.05,{'hsd'}});
raBA = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(raBA)

%% big ANOVA on average over trials
[within,dvn,xlabels,awithinD] = make_within_table({'Conf','Ph','pL'},[5,2,length(binCs)]); %these are trials here
dataT = make_between_table({m_d_eLs_Lin},dvn);
%     ra = RMA(dataT,within,{0.05,{'hsd'}});
raBAT = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(raBAT)
%%

clc
alpha = 0.05/5;
for confi = 1:5
    redF = [1]; redV = {[confi]};
    [dataTR3,withinR] = reduce_within_between(dataT,within,redF,redV);
    raBAR3{confi} = RMA(dataTR3,withinR,{alpha,{''}});
%     raR.ranova
    print_for_manuscript(raBAR3{confi})
end

%% figure bar graphs percent cells trial-wise dists

magfac = mData.magfac;
ff = makeFigureRowsCols(107,[3 5 6.9 1.5],'RowsCols',[1 5],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.3],...
    'widthHeightAdjustment',[10 -500]);
MY = 0.24; ysp = 0.0040375; mY = 0.0856; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.28*magfac; widths = [ones(1,5)*1.1]*magfac; gap = 0.0859*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = repmat(mData.dcolors(1:5),1,2); 
for gni = 1:5
    tra = raBAR3{gni};
    axes(ff.h_axes(1,gni));

%     if ismember(gni,[1 5])
%         [xdata,mVar,semVar,combs,p,h,nB] = get_vals_RMA(mData,tra,{'Ph:pL','hsd',0.05},[1 2]);
%     else
        [xdata,mVar,semVar,combs,p,h,nB] = get_vals_RMA(mData,tra,{'pL','hsd',0.05},[1 2]);
%     end
%     ccis = [];
%     for cci = 1:size(combs,1)
%         if ismember(combs(cci,1),[1 2 3 4 5]) & ismember(combs(cci,2),[6 7 8 9 10])
%             ccis = [ccis;cci];
%         end
%     end
%     h(ccis) = 0;
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
% if ismember(gni,[1 5])
%     set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'AOn','AOff'},{[0.001 0.0051]});
% else
    for hbi = 1:length(hbs)
        set(hbs(hbi),'facealpha',0.5,'linestyle',':','edgecolor',tcolors{hbi});
    end
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'AOn/AOff'},{[0.001 0.0051]});
% end
if gni == 3
    hxl = xlabel('Location of peak firing (binned) as a percentage of total phase length'); changePosition(hxl,[0 -0.12 0]);
end

if gni == 1
    ylabel('Probability');
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
set_sub_graph_text(ff,1,titletxts,[0 0.6 0 0],[0.02 0.07 0 0]);


save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);
disp('Done')

%% Find difference in peak locations trial wise 1st order then 2nd order so on and so forth
clc
all_perc_vals = [];
% for 
sel_trialpairs = 1:9;
    all_mVals = []; all_semVals = []; pL_vals_trials_all_C = []; pL_vals_trials_all_CLin = []; pL_vals_trials_all_CLinD = [];
    pL_vals_trials_all_CD = []; mean_pVs = []; all_dpLs = []; all_perc_dpLs = [];
    for cn = 1:length(si)
        all_dist_vals = [];
        pL_vals_trials_all = []; pL_vals_trials_allD = []; pV_vals_trials_all = [];
        pL_vals_trials_allLin = []; pL_vals_trials_allLinD = []; an_perc_dPLs = [];
        combsTR = nchoosek(1:10,2);
        selCombsTR = combsTR(diff(combsTR,[],2) == 1,:);
        
        for an = 1:5
            pLs = eLs{an,cn};% dpLs = diff(pLs,[],2)/size_raster_2(cn);  pVs = allpeakV_trials{an,cn};
            resp = allresp_trials{an,cn};
            minBinD = -1;maxBin = 1;BinWidth = 0.2;  binEsD = minBinD:BinWidth:maxBin; binCsD = binEsD(1:(end-1)) + BinWidth/2;
            pL_vals_trialsD = []; pL_vals_trialsLinD = []; tall_dpLs = []; perc_dpLs = [];
            
            for trNi = 1:length(sel_trialpairs)
                trN = sel_trialpairs(trNi);
                tr1 = selCombsTR(trN,1);
                tr2 = selCombsTR(trN,2);
                respdTR = resp(:,tr1) & resp(:,tr2);
                tdpLs = (pLs(respdTR,tr2) - pLs(respdTR,tr1))/sum(respdTR);
                [bar1D,binEsD] = histcounts(tdpLs,binEsD,'Normalization','probability');
%                     [bar1D,binEsD] = histcounts(tdpLs,binEsD);
%                     bar1D = bar1D/size(pLs,1);
                pL_vals_trialsD = [pL_vals_trialsD;bar1D];
                pL_vals_trialsLinD = [pL_vals_trialsLinD bar1D];
                tall_dpLs = [tall_dpLs;tdpLs];
                
                percNt = sum(tdpLs<0)/length(tdpLs);
                percZt = sum(tdpLs==0)/length(tdpLs);
                percPt = sum(tdpLs>0)/length(tdpLs);
                perc_dpLs = [perc_dpLs percNt percZt percPt];
                all_perc_vals(cn,an,trNi,1) = percNt; all_perc_vals(cn,an,trNi,2) = percZt; all_perc_vals(cn,an,trNi,3) = percPt;
            end
            pL_vals_trials_allD(:,:,an) = pL_vals_trialsD;
            pL_vals_trials_allLinD(an,:) = pL_vals_trialsLinD;
            all_dpLs{an,cn} = tall_dpLs;
            an_perc_dPLs(an,:) = perc_dpLs;
        end
        pL_vals_trials_all_CD{cn} = pL_vals_trials_allD;
        max_mapsD(cn) = max(pL_vals_trials_all_CD{cn}(:));
        min_mapsD(cn) = min(pL_vals_trials_all_CD{cn}(:));
        pL_vals_trials_all_CLinD = [pL_vals_trials_all_CLinD pL_vals_trials_allLinD];
        all_perc_dpLs = [all_perc_dpLs an_perc_dPLs];
    end
mMD = max(max_mapsD); 
mmD = min(min_mapsD); 
% [within,dvn,xlabels,awithinD] = make_within_table({'Conf','Ph','Pol'},[5,2,3]); %these are trials here
% dataT = make_between_table({all_perc_dpLs},dvn);
%     ra = RMA(dataT,within,{0.05,{''}});
%     ra_nzp{sel_trialpairs} = ra;
% % raDpol = RMA(dataT,within,{0.05,{''}});
% print_for_manuscript(ra)
% allps(:,sel_trialpairs) = ra.ranova{[1 3 5 7 9 11 13 15],ra.selected_pval_col};
% end

%% Figure 4A figure heat maps of probabilities trial wise diff pL
magfac = mData.magfac;
ff = makeFigureRowsCols(107,[2 3 6.9 1.5],'RowsCols',[2 10],'spaceRowsCols',[0.04 0.13],'rightUpShifts',[0.03 0.24],...
    'widthHeightAdjustment',[-100 -180]);
MY = 70; ysp = 5; mY = 0; titletxt = 'Activated Cells'; ylabeltxt = {'Cells (%)'}; % for all cells (vals) MY = 80
stp = 0.35*magfac; widths = ([ones(1,10)*0.55]-0.05)*magfac; gap = 0.15*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
% changePosition(ff.h_axes(1,:),[0 0.1 0 -0.05]);
% changePosition(ff.h_axes(2,:),[0 0.3 0 -0.05]);
% changePosition(ff.h_axes(3,:),[0 0 0 0.15]);

[within,dvn,xlabels,awithinD] = make_within_table({'pL'},[length(binCsD)]);
for gn = 1:10
    tvals = pL_vals_trials_all_CD{gn};

    axes(ff.h_axes(1,gn));
    imagesc(binCsD,1:9,tvals(:,:,3),[mmD mMD]);
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
    imagesc(binCsD,1:9,mean(tvals,3),[mmD mMD]);
    set(gca,'YTick',[1 5 10],'YDir','Normal','XTick',binCs);xtickangle(30);
    if gn > 1
        set(gca,'Yticklabels',[]);
    else
        ylabel('Trials');
    end
    if gn == 5
        hxl = xlabel('Difference in location of peak firing between adjacent trials'); changePosition(hxl,[0 0.15 0]);
    end
    xticklabels = {'-L5','-L4','-L3','-L2','-L1','L1','L2','31','L4','L5'};
    set(gca,'xtick',binCsD(1:3:end),'xticklabels',xticklabels(1:3:end)); xtickangle(15);
    
%     set(gca,'xtick',binCsD,'xticklabels',xticklabels); xtickangle(30);

% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'AOn','AOff'},{[0.001 0.0051]});


    format_axes(gca);
    if gn == 10
        hc = putColorBar(ff.h_axes(1,gn),[0.0 0.03 0 -0.05],[mmD mMD],6,'eastoutside',[0.07 0.07 0.1 0.1]);
%         hc = putColorBar(ff.h_axes(2,gn),[0.0 0.03 0 -0.05],[mm mM],6,'eastoutside',[0.07 0.07 0.1 0.1]);
    end
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('peak_firingdists.pdf'),600);

%% big ANOVA Diff
[within,dvn,xlabels,awithinD] = make_within_table({'Conf','Ph','Tr','pL'},[5,2,size(selCombsTR,1),length(binCsD)]); %these are trials here
dataT = make_between_table({pL_vals_trials_all_CLinD},dvn);
%     ra = RMA(dataT,within,{0.05,{'hsd'}});
raBAD = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(raBAD)
%%
ra = raBAD;
tcolors = repmat(mData.colors,1,10);
figure(100);clf; ha = gca;
view_results_rmanova(ha,ra,'Conf:pL','hsd',[1 2],tcolors,[0 1 0.001],mData)
%% 5 ANOVAs for configurations bonferroni adjusted alpha
clc
alpha = 0.05/5;
for confi = 1:5
    redF = [1]; redV = {[confi]};
    [dataTR,withinR] = reduce_within_between(dataT,within,redF,redV);
    raBAR{confi} = RMA(dataTR,withinR,{alpha,{''}});
%     raR.ranova
    print_for_manuscript(raBAR{confi})
end



%% Figure 4D figure bar graphs percent cells trial-wise dists
magfac = mData.magfac;
ff = makeFigureRowsCols(107,[3 5 6.9 1.5],'RowsCols',[1 5],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.3],...
    'widthHeightAdjustment',[10 -450]);
MY = 0.6; ysp = 0.01945125; mY = 0; titletxt = ''; ylabeltxt = {'PDF'}; iplt = 1.5% for all cells (vals) MY = 80
stp = 0.28*magfac; widths = [iplt ones(1,2)*0.8 iplt iplt]*magfac; gap = 0.0859*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = repmat(mData.dcolors(1:10),1,2);
for gni = 1:5
    tra = raBAR{gni};
    axes(ff.h_axes(1,gni));

    if ismember(gni,[1 4 5])
        [xdata,mVar,semVar,combs,p,h,nB] = get_vals_RMA(mData,tra,{'Ph:pL','hsd',0.05},[1 2]);
    else
        [xdata,mVar,semVar,combs,p,h,nB] = get_vals_RMA(mData,tra,{'pL','hsd',0.05},[1 2]);
    end
    ccis = [];
    for cci = 1:size(combs,1)
        if ismember(combs(cci,1),1:10) & ismember(combs(cci,2),11:20)
            ccis = [ccis;cci];
        end
        if ~ismember(combs(cci,1),[5 15 6 16])% & ismember(combs(cci,2),11:20)
            ccis = [ccis;cci];
        end
    end
    h(ccis) = 0;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,'capsize',1,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.25 xdata(end)+0.75],[mY MY]); format_axes(gca); xinds = [1 10]; xticks = xdata; 
if gni > 1
    set(gca,'YTick',[]);
end
make_bars_hollow(hbs(11:end));
% xticklabels = cellstr(num2str((0:10)')); xticklabels = [xticklabels;xticklabels]; xticklabels = xticklabels(xinds);
xticklabels = {'-L5','-L4','-L3','-L2','-L1','L1','L2','31','L4','L5'};
set(gca,'xtick',xdata(1:3:end),'xticklabels',xticklabels(1:3:end)); xtickangle(15);
if ismember(gni,[1 4 5])
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,nB(1),{'AOn','AOff'},{[0.001 0.0051]});
else
    for hbi = 1:length(hbs)
        set(hbs(hbi),'facealpha',0.5,'linestyle',':','edgecolor',tcolors{hbi});
    end
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,nB,{'AOn/AOff'},{[0.001 0.0051]});
end
if gni == 3
    hxl = xlabel('Difference in location of peak firing between adjacent trials'); changePosition(hxl,[0 -0.092 0]);
end

if gni == 1
    ylabel('Probability');
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
set_sub_graph_text(ff,1,titletxts,[0 0.6 0 0],[0.02 0.07 0 0]);


save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);
disp('Done')


