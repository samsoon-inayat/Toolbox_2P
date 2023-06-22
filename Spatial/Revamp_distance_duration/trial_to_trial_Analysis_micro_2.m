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

%% get dists of peaks
minBin = 0;maxBin = 1;BinWidth = 0.2;  binEs = minBin:BinWidth:maxBin; binCs = binEs(1:(end-1)) + BinWidth/2;
cellpops = {'resp','conj','comp1','comp2'};
for ii = 1:length(cellpops)
    [pLdists{ii},pLdistsL{ii},pLdists_T{ii},pLdists_TL{ii},allMI{ii},allMIL{ii},allMITr{ii},allMITrL{ii}] = get_dist_of_peak_locations(cellpops{ii},si,mRsG,allpeakL_trials,allresp_trials,binEs,propsG); % si,an,tr,pL
    % cn an tr pL
end
disp('Done');
%% big ANOVA cell types MI

pLdistsL_cells = [allMIL{2} allMIL{3} allMIL{4}];
pLdistsL_cells = fillmissing(pLdistsL_cells,'linear',2,'EndValues','nearest');
[within,dvn,xlabels,awithinD] = make_within_table({'CT','Conf','Ph','Tr','pL'},[3,5,2,9,length(binCs)]); %these are trials here
dataT = make_between_table({pLdistsL_cells},dvn);
%     ra = RMA(dataT,within,{0.05,{'hsd'}});
raBA = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(raBA)

%% big ANOVA cell types MI Trial Averaged

pLdistsL_cells = [allMITrL{2} allMITrL{3} allMITrL{4}];
pLdistsL_cells = fillmissing(pLdistsL_cells,'linear',2,'EndValues','nearest');
[within,dvn,xlabels,awithinD] = make_within_table({'CT','Conf','Ph','pL'},[3,5,2,length(binCs)]); %these are trials here
dataT = make_between_table({pLdistsL_cells},dvn);
%     ra = RMA(dataT,within,{0.05,{'hsd'}});
raBA = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(raBA)

%% sub anova
clc
raBA_R = RMA_R(raBA,{'Conf','Ph'});
print_for_manuscript(raBA_R)
%%
ra = raBA_R.ras{10};
print_for_manuscript(ra);
tcolors = repmat(mData.colors,1,10);
figure(100);clf; ha = gca;
view_results_rmanova(ha,ra,'CT','hsd',[1 2],tcolors,[0 2 0.071],mData);
%% Figure
clc
magfac = mData.magfac;
ff = makeFigureRowsCols(107,[3 5 6.9 1.3],'RowsCols',[1 1++2+1+2+1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.3],...
    'widthHeightAdjustment',[10 -500]);
MY = 1.75; ysp = 0.1285; mY = 0; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.28*magfac; widths = [ones(1,7)*0.8]*magfac; gap = 0.115*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 2];

tcolors = repmat(mData.dcolors(6:end),1,5); 
[xdata,hbs] = view_results_rmanova(ff.h_axes(1,1),raBA_R.ras{1},'CT','hsd',xs_gaps,tcolors,[mY MY ysp],mData);
xticklabels = {'Conj','Comp 1','Comp 2'}; set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'PL Pooled'},{[0.001 0.0051]});
ylabel('Bits');
axes_title(ff,{1},{'C2 - AOn'},axes_title_shifts_line,axes_title_shifts_text);

tcolors = repmat(mData.dcolors(6:end),1,5); 
[xdata,hbs] = view_results_rmanova(ff.h_axes(1,2),raBA_R.ras{2},'CT','hsd',xs_gaps,tcolors,[mY MY ysp],mData);
xticklabels = {'Conj','Comp 1','Comp 2'}; set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'PL Pooled'},{[0.001 0.0051]});
axes_title(ff,{2},{'C2 - AOff'},axes_title_shifts_line,axes_title_shifts_text);

tcolors = repmat(mData.dcolors(6:end),1,5); 
[xdata,hbs] = view_results_rmanova(ff.h_axes(1,3),raBA_R.ras{4},'CT','hsd',xs_gaps,tcolors,[mY MY ysp],mData);
xticklabels = {'Conj','Comp 1','Comp 2'}; set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'PL Pooled'},{[0.001 0.0051]});
axes_title(ff,{3},{'C3 - AOff'},axes_title_shifts_line,axes_title_shifts_text);

tcolors = repmat(mData.dcolors(6:end),1,5); 
[xdata,hbs] = view_results_rmanova(ff.h_axes(1,4),raBA_R.ras{5},'CT','hsd',xs_gaps,tcolors,[mY MY ysp],mData);
xticklabels = {'Conj','Comp 1','Comp 2'}; set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'PL Pooled'},{[0.001 0.0051]});
axes_title(ff,{4},{'C4 - AOn'},axes_title_shifts_line,axes_title_shifts_text);

tcolors = repmat(mData.dcolors(6:end),1,5); 
[xdata,hbs] = view_results_rmanova(ff.h_axes(1,5),raBA_R.ras{6},'CT','hsd',xs_gaps,tcolors,[mY MY ysp],mData);
xticklabels = {'Conj','Comp 1','Comp 2'}; set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'PL Pooled'},{[0.001 0.0051]});
axes_title(ff,{5},{'C4 - AOff'},axes_title_shifts_line,axes_title_shifts_text);

tcolors = repmat(mData.dcolors(6:end),1,5); 
[xdata,hbs] = view_results_rmanova(ff.h_axes(1,6),raBA_R.ras{8},'CT','hsd',xs_gaps,tcolors,[mY MY ysp],mData);
xticklabels = {'Conj','Comp 1','Comp 2'}; set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'PL Pooled'},{[0.001 0.0051]});
axes_title(ff,{6},{'C5 - AOff'},axes_title_shifts_line,axes_title_shifts_text);

tcolors = repmat(mData.dcolors(6:end),1,5); 
[xdata,hbs] = view_results_rmanova(ff.h_axes(1,7),raBA_R.ras{10},'CT','hsd',xs_gaps,tcolors,[mY MY ysp],mData);
xticklabels = {'Conj','Comp 1','Comp 2'}; set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'PL Pooled'},{[0.001 0.0051]});
axes_title(ff,{7},{'C7 - AOff'},axes_title_shifts_line,axes_title_shifts_text);

save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graph.pdf'),600);


%% sub anova
clc
raBA_R1 = RMA_R(raBA_R.ras{1},{'Ph'});
print_for_manuscript(raBA_R1)
%% big ANOVA cell types

pLdistsL_cells = [pLdistsL{2} pLdistsL{3} pLdistsL{4}];
[within,dvn,xlabels,awithinD] = make_within_table({'CT','Conf','Ph','Tr','pL'},[3,5,2,9,length(binCs)]); %these are trials here
dataT = make_between_table({pLdistsL_cells},dvn);
%     ra = RMA(dataT,within,{0.05,{'hsd'}});
raBA = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(raBA)

%% sub anova
clc
raBA_R = RMA_R(raBA,{'CT'});
print_for_manuscript(raBA_R)
%%
ra = raBA_R.ras{5};
print_for_manuscript(ra);
tcolors = repmat(mData.colors,1,10);
figure(100);clf; ha = gca;
view_results_rmanova(ha,ra,'pL','hsd',[1 2],tcolors,[0 2 0.1],mData);

%% big anova after trial averaging
pLdistsL_cells = [pLdists_TL{2} pLdists_TL{3} pLdists_TL{4}];
[within,dvn,xlabels,awithinD] = make_within_table({'CT','Conf','Ph','pL'},[3,5,2,length(binCs)]); %these are trials here
dataT = make_between_table({pLdistsL_cells},dvn);
%     ra = RMA(dataT,within,{0.05,{'hsd'}});
raBA = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(raBA)
%%
%%
clc
raBA_R = RMA_R(raBA,{'Conf','Ph'});
print_for_manuscript(raBA_R)
%% figure bar graphs percent cells trial-wise dists

magfac = mData.magfac;
ff = makeFigureRowsCols(107,[3 5 6.9 1.25],'RowsCols',[1 10],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.35],...
    'widthHeightAdjustment',[10 -500]);
MY = 1; ysp = 0.05; mY = 0; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.28*magfac; widths = [ones(1,10)*0.55]*magfac; gap = 0.115*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 2];
tcolors = repmat(mData.dcolors(1:5),1,3); 
for gni = 1:10
    tra = raBA_R.ras{gni};
    [xdata,hbs] = view_results_rmanova(ff.h_axes(1,gni),raBA_R.ras{gni},{'CT:pL','hsd',0.05},xs_gaps,tcolors,[mY MY ysp],mData);

    if gni > 1
        set(gca,'YTick',[]);
    end
    if ~mod(gni,2)
        make_bars_hollow(hbs);
    end
    % xticklabels = cellstr(num2str((0:10)')); xticklabels = [xticklabels;xticklabels]; xticklabels = xticklabels(xinds);
    xticklabels = {'L1','L2','L3','L4','L5'};
    set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
    if gni == 1
        ylabel('Bits');
    end

    if gni == 5
         hxl = xlabel('Location of peak firing (binned) as a percentage of total phase length'); changePosition(hxl,[0 -0.37 0]);
    end
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'Conj','Comp 1','Comp 2'},{[0.001 0.0051]});
    box off;
    format_axes(gca);
end

confnames = {'C2','C3','C4','C5','C7'};
for ii = 1:5
    titletxts{ii} = sprintf('%s',confnames{ii});
end
set_sub_graph_text(ff,2,titletxts,[0 0.56 0 0],[0.02 0.07 0 0]);


save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);
disp('Done')

%% Figure
magfac = mData.magfac;
ff = makeFigureRowsCols(107,[3 5 6.9 1.3],'RowsCols',[1 1+1+1+1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.3],...
    'widthHeightAdjustment',[10 -500]);
MY = 0.43; ysp = 0.035285; mY = 0; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.28*magfac; widths = [1.35 1 2.85 1]*magfac; gap = 0.115*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 2];

tcolors = repmat(mData.colors(1:2),1,3); 
[xdata,hbs] = view_results_rmanova(ff.h_axes(1,1),raBA_R1.ras{5},'CT:Ph','hsd',xs_gaps,tcolors,[mY MY ysp],mData);
xticklabels = {'AOn','AOff'}; set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Conj','Comp1','Comp2'},{[0.001 0.0051]});
ylabel('Probability');
axes_title(ff,{1},{'C2 - L5'},axes_title_shifts_line,axes_title_shifts_text);

tcolors = repmat(mData.dcolors(1:5),1,3); 
[xdata,hbs] = view_results_rmanova(ff.h_axes(1,2),raBA_R.ras{2},'pL','hsd',xs_gaps,tcolors,[mY MY ysp],mData);
xticklabels = {'L1','L2','L3','L4','L5'};
set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'Conj','Comp1','Comp2'},{[0.001 0.0051]});
axes_title(ff,{2},{'C3'},axes_title_shifts_line,axes_title_shifts_text);


tcolors = repmat(mData.dcolors(1:5),1,3); 
[xdata,hbs] = view_results_rmanova(ff.h_axes(1,3),raBA_R2.ras{1},'CT:pL','hsd',xs_gaps,tcolors,[mY MY ysp],mData);
xticklabels = {'L1','L2','L3','L4','L5'};
set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'Conj','Comp1','Comp2'},{[0.001 0.0051]});
axes_title(ff,{3},{'C4 - AOn'},axes_title_shifts_line,axes_title_shifts_text);

tcolors = repmat(mData.dcolors(1:5),1,3); 
[xdata,hbs] = view_results_rmanova(ff.h_axes(1,4),raBA_R.ras{4},'pL','hsd',xs_gaps,tcolors,[mY MY ysp],mData);
xticklabels = {'L1','L2','L3','L4','L5'};
set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'Conj','Comp1','Comp2'},{[0.001 0.0051]});
axes_title(ff,{4},{'C5'},axes_title_shifts_line,axes_title_shifts_text);

save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graph.pdf'),600);

%% anova w.r.t pL
clc
raBA_R = RMA_R(raBA,{'pL'});
print_for_manuscript(raBA_R)

%% sub anova
clc
raBA_R1 = RMA_R(raBA_R.ras{1},{'Ph'});
print_for_manuscript(raBA_R1)

%% Figure
magfac = mData.magfac;
ff = makeFigureRowsCols(107,[3 5 6.9 1.3],'RowsCols',[1 1+1+2+1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.3],...
    'widthHeightAdjustment',[10 -500]);
MY = 0.46; ysp = 0.0275; mY = 0; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.28*magfac; widths = [2.325 1.25 0.85 0.425 1.25]*magfac; gap = 0.115*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 2];

tcolors = repmat(mData.colors(6:10),1,3); 
[xdata,hbs] = view_results_rmanova(ff.h_axes(1,1),raBA_R1.ras{1},'CT:Conf','hsd',xs_gaps,tcolors,[mY MY ysp],mData);
ylabel('Probability');
xticklabels = {'C2','C3','C4','C5','C7'};
set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'Conj','Comp1','Comp2'},{[0.001 0.0051]});
axes_title(ff,{1},{'L1'},axes_title_shifts_line,axes_title_shifts_text);

tcolors = repmat(mData.colors(1:2),1,3); 
[xdata,hbs] = view_results_rmanova(ff.h_axes(1,2),raBA_R.ras{2},'CT:Ph','hsd',xs_gaps,tcolors,[mY MY ysp],mData);
xticklabels = {'AOn','AOff'}; set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Conj','Comp1','Comp2'},{[0.001 0.0051]});
axes_title(ff,{2},{'L2'},axes_title_shifts_line,axes_title_shifts_text);

tcolors = repmat(mData.colors(6:10),1,3); 
[xdata,hbs] = view_results_rmanova(ff.h_axes(1,3),raBA_R.ras{4},'Conf','hsd',xs_gaps,tcolors,[mY MY ysp],mData);
xticklabels = {'C2','C3','C4','C5','C7'}; set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);

tcolors = repmat(mData.colors(1:2),1,3); 
[xdata,hbs] = view_results_rmanova(ff.h_axes(1,4),raBA_R.ras{4},'Ph','hsd',xs_gaps,tcolors,[mY MY ysp],mData);
xticklabels = {'AOn','AOff'}; set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
axes_title(ff,{3:4},{'L4'},axes_title_shifts_line,axes_title_shifts_text);

tcolors = repmat(mData.colors(1:2),1,3); 
[xdata,hbs] = view_results_rmanova(ff.h_axes(1,5),raBA_R.ras{5},'CT:Ph','hsd',xs_gaps,tcolors,[mY MY ysp],mData);
xticklabels = {'AOn','AOff'}; set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Conj','Comp1','Comp2'},{[0.001 0.0051]});
axes_title(ff,{5},{'L5'},axes_title_shifts_line,axes_title_shifts_text);

save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graph.pdf'),600);


%% anova w.r.t. CT
clc
raBA_R = RMA_R(raBA,{'CT'});
print_for_manuscript(raBA_R)

%% sub anova
clc
raBA_R1 = RMA_R(raBA_R.ras{1},{'CT'});
print_for_manuscript(raBA_R1)

%% Figure
clc
magfac = mData.magfac;
ff = makeFigureRowsCols(107,[3 5 6.9 2],'RowsCols',[1 1+1+1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.2],...
    'widthHeightAdjustment',[10 -350]);
titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.28*magfac; widths = [1 2.6 2.6]*magfac; gap = 0.115*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
MY = 0.65; ysp = 0.0275; mY = 0; 
axes_title_shifts_line = [0 0.675 0 0]; axes_title_shifts_text = [0.02 0.071 0 0]; xs_gaps = [1 2];

tcolors = repmat(mData.dcolors(1:5),1,3); 
[xdata,hbs] = view_results_rmanova(ff.h_axes(1,1),raBA_R.ras{1},'Ph:pL','hsd',xs_gaps,tcolors,[mY MY ysp],mData);
ylabel('Probability');
xticklabels = {'L1','L2','L3','L4','L5'};
set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'AOn','AOff'},{[0.001 0.0051]});
axes_title(ff,{1},{'Conj'},axes_title_shifts_line,axes_title_shifts_text);

tcolors = repmat(mData.dcolors(1:5),1,5); 
[xdata,hbs] = view_results_rmanova(ff.h_axes(1,2),raBA_R.ras{2},'Conf:pL','hsd',xs_gaps,tcolors,[mY MY ysp],mData);
xticklabels = {'L1','L2','L3','L4','L5'}; set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'C2','C3','C4','C5','C7'},{[0.001 0.0051]});
axes_title(ff,{2},{'Comp1'},axes_title_shifts_line,axes_title_shifts_text);

tcolors = repmat(mData.dcolors(1:5),1,5); 
[xdata,hbs] = view_results_rmanova(ff.h_axes(1,3),raBA_R.ras{3},'Conf:pL','hsd',xs_gaps,tcolors,[mY MY ysp],mData);
xticklabels = {'L1','L2','L3','L4','L5'}; set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'C2','C3','C4','C5','C7'},{[0.001 0.0051]});
axes_title(ff,{3},{'Comp2'},axes_title_shifts_line,axes_title_shifts_text);


save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graph.pdf'),600);

%% random visualize results
ra = raBA_R.ras{2};
tcolors = repmat(mData.colors,1,10);
figure(300);clf; ha = gca;
view_results_rmanova(ha,ra,'Conf:pL','hsd',[1 2],tcolors,[0 0.3 0.05],mData);

%%
ra = raBA_R1.ras{5};
tcolors = repmat(mData.colors,1,10);
figure(300);clf; ha = gca;
view_results_rmanova(ha,ra,'CT:Ph','hsd',[1 2],tcolors,[0 0.3 0.05],mData)
%%
clc
close(figure(100));
tra = raBA_R.ras{1};
ra_pL1 = RMA_R(tra,{'Ph'});
print_for_manuscript(ra_pL1)
view_results_rmanova([],ra_pL1.ras{1},'CT:Conf','hsd',[1 2],tcolors,[0 1 0.05],mData)

%% Figure 5A figure heat maps of probabilities trial wise 3 rows Conj, Comp1, Comp2
magfac = mData.magfac;
ff = makeFigureRowsCols(107,[2 3 6.9 2],'RowsCols',[3 10],'spaceRowsCols',[0.04 0.13],'rightUpShifts',[0.03 0.178],...
    'widthHeightAdjustment',[-100 -110]);
MY = 70; ysp = 5; mY = 0; titletxt = 'Activated Cells'; ylabeltxt = {'Cells (%)'}; % for all cells (vals) MY = 80
stp = 0.35*magfac; widths = ([ones(1,10)*0.55]-0.05)*magfac; gap = 0.15*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
rep_animal = 0; % if 0 then plot average
rep_an = 3;
for cti = 2:4
    ttvals = pLdists{cti}; % cn an tr pL
    ttvals = permute(ttvals,[3 4 1 2]);
    size(ttvals);
    if rep_animal
        mttvals = ttvals(:,:,:,rep_an);
    else
        mttvals = squeeze(mean(ttvals,4));
    end
    amm(cti) = min(mttvals(:));
    amM(cti) = max(mttvals(:));
end
mM = max(amM)/1.25; mm = min(amm);
for cti = 2:4
    ttvals = pLdists{cti}; % cn an tr pL
    ttvals = permute(ttvals,[3 4 1 2]);
    size(ttvals);
    if rep_animal
        mttvals = ttvals(:,:,:,rep_an);
    else
        mttvals = squeeze(mean(ttvals,4));
    end
%     mm = amm(cti); mM = amM(cti);
for gn = 1:10
    axes(ff.h_axes(cti-1,gn));
    if cti == 1
        imagesc(binCs,1:10,mttvals(:,:,gn),[mm mM]);
        set(gca,'YTick',[1 5 10],'YDir','Normal','XTick',binCs);xtickangle(30);
    else
        imagesc(binCs,1:9,mttvals(:,:,gn),[mm mM]);
        set(gca,'YTick',[1 9],'YDir','Normal','XTick',binCs);xtickangle(30);
    end
    if gn > 1
        set(gca,'Yticklabels',[]);
    else
        switch cti
            case 1
                ylabel({'Activ','Trials'});
            case 2
                ylabel({'Conj','Trials'});
            case 3
                ylabel({'Comp 1','Trials'});
            case 4
                ylabel({'Comp 2','Trials'});
        end
    end
    if gn == 5 && cti-1 == 3
%         xlabel('Location of Peak Firing as a Fraction of total Phase Length');
        hxl = xlabel('Location of peak firing (binned) as a percentage of total phase length'); changePosition(hxl,[0 0.51 0]);
    end
    xticklabels = {'L1','L2','L3','L4','L5'};
    if cti-1 == 3
        set(gca,'xtick',binCs,'xticklabels',xticklabels); xtickangle(30);
    else
        set(gca,'xtick',binCs,'xticklabels',[]); xtickangle(30);
    end
    if cti == 2
    titletxt = (sprintf('%s',rasterNamesTxt{si(gn)}));
    ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.1 0 0]);set(ht,'FontWeight','NOrmal');
    end
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'AOn','AOff'},{[0.001 0.0051]});


    format_axes(gca);
    if gn == 10 && cti == 2
        hc = putColorBar(ff.h_axes(cti-1,gn),[0.0 0.03 0 -0.05],[mm mM],6,'eastoutside',[0.1 0.07 0.1 0.1]);
%         hc = putColorBar(ff.h_axes(2,gn),[0.0 0.03 0 -0.05],[mm mM],6,'eastoutside',[0.07 0.07 0.1 0.1]);
    end
end
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('peak_firingdists.pdf'),600);

%% Figure 5A figure heat maps of MI trial wise 3 rows Conj, Comp1, Comp2
magfac = mData.magfac;
ff = makeFigureRowsCols(107,[2 3 6.9 2],'RowsCols',[3 10],'spaceRowsCols',[0.04 0.13],'rightUpShifts',[0.03 0.178],...
    'widthHeightAdjustment',[-100 -110]);
MY = 70; ysp = 5; mY = 0; titletxt = 'Activated Cells'; ylabeltxt = {'Cells (%)'}; % for all cells (vals) MY = 80
stp = 0.35*magfac; widths = ([ones(1,10)*0.55]-0.05)*magfac; gap = 0.15*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
rep_animal = 0; % if 0 then plot average
rep_an = 3;
for cti = 2:4
    ttvals = allMI{cti}; % cn an tr pL
    ttvals = permute(ttvals,[3 4 1 2]);
    size(ttvals);
    if rep_animal
        mttvals = ttvals(:,:,:,rep_an);
    else
        mttvals = squeeze(nanmean(ttvals,4));
    end
    amm(cti) = min(mttvals(:));
    amM(cti) = max(mttvals(:));
end
mM = max(amM); mm = min(amm);
for cti = 2:4
    ttvals = allMI{cti}; % cn an tr pL
    ttvals = permute(ttvals,[3 4 1 2]);
    size(ttvals);
    if rep_animal
        mttvals = ttvals(:,:,:,rep_an);
    else
        mttvals = squeeze(mean(ttvals,4));
    end
    mttvals = fillmissing(mttvals,'linear',2,'EndValues','nearest');
%     mm = amm(cti); mM = amM(cti);
for gn = 1:10
    axes(ff.h_axes(cti-1,gn));
    if cti == 1
        imagesc(binCs,1:10,mttvals(:,:,gn),[mm mM]);
        set(gca,'YTick',[1 5 10],'YDir','Normal','XTick',binCs);xtickangle(30);
    else
        imagesc(binCs,1:9,mttvals(:,:,gn),[mm mM]);
        set(gca,'YTick',[1 9],'YDir','Normal','XTick',binCs);xtickangle(30);
    end
    if gn > 1
        set(gca,'Yticklabels',[]);
    else
        switch cti
            case 1
                ylabel({'Activ','Trials'});
            case 2
                ylabel({'Conj','Trials'});
            case 3
                ylabel({'Comp 1','Trials'});
            case 4
                ylabel({'Comp 2','Trials'});
        end
    end
    if gn == 5 && cti-1 == 3
%         xlabel('Location of Peak Firing as a Fraction of total Phase Length');
        hxl = xlabel('Location of peak firing (binned) as a percentage of total phase length'); changePosition(hxl,[0 0.51 0]);
    end
    xticklabels = {'L1','L2','L3','L4','L5'};
    if cti-1 == 3
        set(gca,'xtick',binCs,'xticklabels',xticklabels); xtickangle(30);
    else
        set(gca,'xtick',binCs,'xticklabels',[]); xtickangle(30);
    end
    if cti == 2
    titletxt = (sprintf('%s',rasterNamesTxt{si(gn)}));
    ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.1 0 0]);set(ht,'FontWeight','NOrmal');
    end
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'AOn','AOff'},{[0.001 0.0051]});


    format_axes(gca);
    if gn == 10 && cti == 2
        hc = putColorBar(ff.h_axes(cti-1,gn),[0.0 0.03 0 -0.05],[mm mM],6,'eastoutside',[0.1 0.07 0.1 0.1]);
%         hc = putColorBar(ff.h_axes(2,gn),[0.0 0.03 0 -0.05],[mm mM],6,'eastoutside',[0.07 0.07 0.1 0.1]);
    end
end
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('peak_firingdists.pdf'),600);

%% Figure 5A figure heat maps of MI trial wise
magfac = mData.magfac;
ff = makeFigureRowsCols(107,[2 3 6.9 2.25],'RowsCols',[4 10],'spaceRowsCols',[0.04 0.13],'rightUpShifts',[0.03 0.16],...
    'widthHeightAdjustment',[-100 -90]);
MY = 70; ysp = 5; mY = 0; titletxt = 'Activated Cells'; ylabeltxt = {'Cells (%)'}; % for all cells (vals) MY = 80
stp = 0.35*magfac; widths = ([ones(1,10)*0.55]-0.05)*magfac; gap = 0.15*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
rep_animal = 1; % if 0 then plot average
rep_an = 3;
for cti = 1:4
    ttvals = allMI{cti}; % cn an tr pL
    ttvals = permute(ttvals,[3 4 1 2]);
    size(ttvals);
    if rep_animal
        mttvals = ttvals(:,:,:,rep_an);
    else
        mttvals = squeeze(mean(ttvals,4));
    end
    amm(cti) = min(mttvals(:));
    amM(cti) = max(mttvals(:));
end
mM = max(amM); mm = min(amm);
for cti = 1:4
    ttvals = allMI{cti}; % cn an tr pL
    ttvals = permute(ttvals,[3 4 1 2]);
    size(ttvals);
    if rep_animal
        mttvals = ttvals(:,:,:,rep_an);
    else
        mttvals = squeeze(mean(ttvals,4));
    end
%     mm = amm(cti); mM = amM(cti);
for gn = 1:10
    axes(ff.h_axes(cti,gn));
    if cti == 1
        imagesc(binCs,1:10,mttvals(:,:,gn),[mm mM]);
        set(gca,'YTick',[1 5 10],'YDir','Normal','XTick',binCs);xtickangle(30);
    else
        imagesc(binCs,1:9,mttvals(:,:,gn),[mm mM]);
        set(gca,'YTick',[1 9],'YDir','Normal','XTick',binCs);xtickangle(30);
    end
    if gn > 1
        set(gca,'Yticklabels',[]);
    else
        switch cti
            case 1
                ylabel({'Activ','Trials'});
            case 2
                ylabel({'Conj','Trials'});
            case 3
                ylabel({'Comp 1','Trials'});
            case 4
                ylabel({'Comp 2','Trials'});
        end
    end
    if gn == 5 && cti == 4
%         xlabel('Location of Peak Firing as a Fraction of total Phase Length');
        hxl = xlabel('Location of peak firing (binned) as a percentage of total phase length'); changePosition(hxl,[0 0.51 0]);
    end
    xticklabels = {'L1','L2','L3','L4','L5'};
    if cti == 4
        set(gca,'xtick',binCs,'xticklabels',xticklabels); xtickangle(30);
    else
        set(gca,'xtick',binCs,'xticklabels',[]); xtickangle(30);
    end
    if cti == 1
    titletxt = (sprintf('%s',rasterNamesTxt{si(gn)}));
    ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.1 0 0]);set(ht,'FontWeight','NOrmal');
    end
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'AOn','AOff'},{[0.001 0.0051]});


    format_axes(gca);
    if gn == 10 && cti == 2
        hc = putColorBar(ff.h_axes(cti-1,gn),[0.0 0.03 0 -0.05],[mm mM],6,'eastoutside',[0.1 0.07 0.1 0.1]);
%         hc = putColorBar(ff.h_axes(2,gn),[0.0 0.03 0 -0.05],[mm mM],6,'eastoutside',[0.07 0.07 0.1 0.1]);
    end
end
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('peak_firingdists.pdf'),600);

%% Figure 5A figure heat maps of probabilities trial wise (only resp top row)
magfac = mData.magfac;
ff = makeFigureRowsCols(107,[2 3 6.9 1],'RowsCols',[1 10],'spaceRowsCols',[0.04 0.13],'rightUpShifts',[0.03 0.336],...
    'widthHeightAdjustment',[-100 -500]);
MY = 70; ysp = 5; mY = 0; titletxt = 'Activated Cells'; ylabeltxt = {'Cells (%)'}; % for all cells (vals) MY = 80
stp = 0.35*magfac; widths = ([ones(1,10)*0.55]-0.05)*magfac; gap = 0.15*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
rep_animal = 1; % if 0 then plot average
rep_an = 3;
for cti = 1%:4
    ttvals = pLdists{cti}; % cn an tr pL
    ttvals = permute(ttvals,[3 4 1 2]);
    size(ttvals);
    if rep_animal
        mttvals = ttvals(:,:,:,rep_an);
    else
        mttvals = squeeze(mean(ttvals,4));
    end
    amm(cti) = min(mttvals(:));
    amM(cti) = max(mttvals(:));
end
mM = max(amM)/1.25; mm = min(amm);
for cti = 1%:4
    ttvals = pLdists{cti}; % cn an tr pL
    ttvals = permute(ttvals,[3 4 1 2]);
    size(ttvals);
    if rep_animal
        mttvals = ttvals(:,:,:,rep_an);
    else
        mttvals = squeeze(mean(ttvals,4));
    end
%     mm = amm(cti); mM = amM(cti);
for gn = 1:10
    axes(ff.h_axes(cti,gn));
    if cti == 1
        imagesc(binCs,1:10,mttvals(:,:,gn),[mm mM]);
        set(gca,'YTick',[1 5 10],'YDir','Normal','XTick',binCs);xtickangle(30);
    else
        imagesc(binCs,1:9,mttvals(:,:,gn),[mm mM]);
        set(gca,'YTick',[1 9],'YDir','Normal','XTick',binCs);xtickangle(30);
    end
    if gn > 1
        set(gca,'Yticklabels',[]);
    else
        switch cti
            case 1
                ylabel({'Trials'});
            case 2
                ylabel({'Conj','Trials'});
            case 3
                ylabel({'Comp 1','Trials'});
            case 4
                ylabel({'Comp 2','Trials'});
        end
    end
    if gn == 5 && cti == 1
%         xlabel('Location of Peak Firing as a Fraction of total Phase Length');
        hxl = xlabel('Location of peak firing (binned) as a percentage of total phase length'); changePosition(hxl,[0 0.51 0]);
    end
    xticklabels = {'L1','L2','L3','L4','L5'};
    if cti == 1
        set(gca,'xtick',binCs,'xticklabels',xticklabels); xtickangle(30);
    else
        set(gca,'xtick',binCs,'xticklabels',[]); xtickangle(30);
    end
    if cti == 1
    titletxt = (sprintf('%s',rasterNamesTxt{si(gn)}));
    ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.051 0 0]);set(ht,'FontWeight','NOrmal');
    end
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'AOn','AOff'},{[0.001 0.0051]});


    format_axes(gca);
    if gn == 10 && cti == 1
        hc = putColorBar(ff.h_axes(cti,gn),[0.0 0.03 0 -0.05],[mm mM],6,'eastoutside',[0.1 0.07 0.1 0.1]);
%         hc = putColorBar(ff.h_axes(2,gn),[0.0 0.03 0 -0.05],[mm mM],6,'eastoutside',[0.07 0.07 0.1 0.1]);
    end
end
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('peak_firingdists.pdf'),600);

%% Figure 5A figure heat maps of MI (only resp top row)
magfac = mData.magfac;
ff = makeFigureRowsCols(107,[2 3 6.9 1],'RowsCols',[1 10],'spaceRowsCols',[0.04 0.13],'rightUpShifts',[0.03 0.336],...
    'widthHeightAdjustment',[-100 -500]);
MY = 70; ysp = 5; mY = 0; titletxt = 'Activated Cells'; ylabeltxt = {'Cells (%)'}; % for all cells (vals) MY = 80
stp = 0.35*magfac; widths = ([ones(1,10)*0.55]-0.05)*magfac; gap = 0.15*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
rep_animal = 1; % if 0 then plot average
rep_an = 3;
for cti = 1%:4
    ttvals = allMI{cti}; % cn an tr pL
    ttvals = permute(ttvals,[3 4 1 2]);
    size(ttvals);
    if rep_animal
        mttvals = ttvals(:,:,:,rep_an);
    else
        mttvals = squeeze(mean(ttvals,4));
    end
    amm(cti) = min(mttvals(:));
    amM(cti) = max(mttvals(:));
end
mM = max(amM)/1.25; mm = min(amm);
for cti = 1%:4
    ttvals = allMI{cti}; % cn an tr pL
    ttvals = permute(ttvals,[3 4 1 2]);
    size(ttvals);
    if rep_animal
        mttvals = ttvals(:,:,:,rep_an);
    else
        mttvals = squeeze(mean(ttvals,4));
    end
%     mm = amm(cti); mM = amM(cti);
for gn = 1:10
    axes(ff.h_axes(cti,gn));
    if cti == 1
        imagesc(binCs,1:10,mttvals(:,:,gn),[mm mM]);
        set(gca,'YTick',[1 5 10],'YDir','Normal','XTick',binCs);xtickangle(30);
    else
        imagesc(binCs,1:9,mttvals(:,:,gn),[mm mM]);
        set(gca,'YTick',[1 9],'YDir','Normal','XTick',binCs);xtickangle(30);
    end
    if gn > 1
        set(gca,'Yticklabels',[]);
    else
        switch cti
            case 1
                ylabel({'Trials'});
            case 2
                ylabel({'Conj','Trials'});
            case 3
                ylabel({'Comp 1','Trials'});
            case 4
                ylabel({'Comp 2','Trials'});
        end
    end
    if gn == 5 && cti == 1
%         xlabel('Location of Peak Firing as a Fraction of total Phase Length');
        hxl = xlabel('Location of peak firing (binned) as a percentage of total phase length'); changePosition(hxl,[0 0.51 0]);
    end
    xticklabels = {'L1','L2','L3','L4','L5'};
    if cti == 1
        set(gca,'xtick',binCs,'xticklabels',xticklabels); xtickangle(30);
    else
        set(gca,'xtick',binCs,'xticklabels',[]); xtickangle(30);
    end
    if cti == 1
    titletxt = (sprintf('%s',rasterNamesTxt{si(gn)}));
    ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.051 0 0]);set(ht,'FontWeight','NOrmal');
    end
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'AOn','AOff'},{[0.001 0.0051]});


    format_axes(gca);
    if gn == 10 && cti == 1
        hc = putColorBar(ff.h_axes(cti,gn),[0.0 0.03 0 -0.05],[mm mM],6,'eastoutside',[0.1 0.07 0.1 0.1]);
%         hc = putColorBar(ff.h_axes(2,gn),[0.0 0.03 0 -0.05],[mm mM],6,'eastoutside',[0.07 0.07 0.1 0.1]);
    end
end
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('peak_firingdists.pdf'),600);

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
%     axes(ff.h_axes(2,gn));
%     imagesc(1:length(binCs),1:10,m_pL_vals_trials_all{gn},[mm mM]);
    
    tvals = pL_vals_trials_all_CD{gn};
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
    imagesc(binCsD,1:10,tvals(:,:,3),[mmD mMD]);
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
    imagesc(binCsD,1:10,m_pL_vals_trials_allD{gn},[mmD mMD]);
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
    
    %% visualizing overall correlation
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


