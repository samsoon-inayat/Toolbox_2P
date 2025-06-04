function time_trial_distance_speed
%% Load Data
udata = evalin('base','udata1');
ei = evalin('base','ei');

outT = get_the_binned_data(udata,'time',0.3);
outD = get_the_binned_data(udata,'dist',3);

filename = 'FR_tuning.mat';
if exist(filename,'file')
    load(filename);
else
    % Get the metrics and run stats
    variable_combs = {'FR_time','FR_dist','FR_speed'};
    trialsOrConcat = 'concatenate'; %'trials' or "concatenate"
    nshuffles = 1000;
    clc
    for vn = 1:length(variable_combs)
        vn
        met_valsT{vn} = get_metrics_FR(outT,variable_combs{vn},trialsOrConcat,nshuffles);
        met_valsD{vn} = get_metrics_FR(outD,variable_combs{vn},trialsOrConcat,nshuffles);
    end
    
    save(filename,"met_valsT","met_valsD");
end

% cell_typesT = get_cell_types(outT,met_valsT);
% cell_typesD = get_cell_types(outD,met_valsD);

n = 0;

%% find the lengths of the bins in air on and off phases
alensig = [];
for an = 1:5
    len_sig = [];
    for cn = 1:3
        for ap = 1:2
            if ap == 2
                out = outT;
            else
                out = outD;
            end
            len_sig = [len_sig length(out.atimecc{an,cn,ap})];
        end
    end
    alensig = [alensig;len_sig];
end
%% Get the metrics and run stats (Time-Speed and Distance-Speed
variable_combs = {'time_speed','dist_speed'};
met = 'MI'; trialsOrConcat = 'concatenate'; %'trials' or "concatenate"
avar = [];
for vn = 1:length(variable_combs)
    avar = [avar get_metrics(outT,variable_combs{vn},met,trialsOrConcat)];
end

for vn = 1:length(variable_combs)
    avar = [avar get_metrics(outD,variable_combs{vn},met,trialsOrConcat)];
end

clc
descriptiveStatistics(avar(:)); %(mean ± sem: 0.281 ± 0.022, range: -0.957, 0.999, median: 0.485) concatenate MI (mean ± sem: 0.746 ± 0.051, range: 0.187, 1.998, median: 0.593)
% (Intercept):PT [F(2,8) = 116.50, p < 0.001, η2 = .84] <--
fac_names = {'BT','PT','CN','AP','TN'}; fac_levels = [2,2,3,2,10];
% fac_names = {'PT','CN','AP'}; fac_levels = [3,3,2];
% fac_names = {'CN','PT'}; fac_levels = [3,3];
fac_names = {'BT','PT','CN','AP'}; fac_levels = [2,2,3,2];
[within,dvn,xlabels,awithinD] = make_within_table(fac_names,fac_levels);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)

%% Run Reduced Subset RM-ANOVA
clc
raR = RMA_subset(ra,'AP');
% Run Bonferroni adjusted alpha RM-ANOVA
clc
raRB = RMA_bonferroni(raR{1},'BT');

%% Plot rasters
out = outT;
for an = 1:5
    for cn = 1:3
        for ap = 2
            rasters = out.arasters{an,cn,ap}; rx = out.arx{an,cn,ap}; ry = out.ary{an,cn,ap};
            num_neurons = size(rasters,3);
            for nn = 1:num_neurons
                figure(100);clf;
                imagesc(rx,ry,rasters(:,:,nn));colorbar;
                pause(0.3);
            end
        end
    end
end

%% Plot active cells heatmap
tmvT = met_valsT{1};
mctFR = zeros(500,500);
for an = 1:5
    for cn = 1
        for ap = 1
            if ap == 1
                out = outD;
            end
            if ap == 2
                out = outT;
            end
            tPCT = ~isnan(tmvT{an,cn,ap}.PC(:,1));
            tFR = (out.aFRcc{an,cn,ap})'; tTrials = out.atrialcc{an,cn,ap};
            tFR = tFR(tPCT,:);
            [rasters,x,y] = build_rasters(tTrials,tFR);
            mRasters = (squeeze(nanmean(rasters,1)))';
            [~,midx] = max(mRasters,[],2);
            [~,smidx] = sort(midx);
            mmRasters = mRasters(smidx,:);
            % figure(100);clf;imagesc(x,y,mmRasters);set(gca,'Ydir','normal');colorbar;
            % ctFR = corr(mmRasters);
            [~,midx] = max(tFR,[],2);
            [~,smidx] = sort(midx);
            mtFR = tFR(smidx,:);
            % figure(100);clf;imagesc(tTrials,1:size(tFR,1),mtFR,[0 80]);set(gca,'Ydir','normal');colorbar;
            ctFR = corr(tFR);
            % figure(200);clf;imagesc(tTrials,tTrials,ctFR); set(gca,'Ydir','normal');colorbar;

        end
    end
    sz_ctFR = size(ctFR,1);
    if size(ctFR,1) < 500
        ctFR = padarray(ctFR, [500-size(ctFR,1), 500-size(ctFR,1)], 'replicate', 'post'); 
    end
    mctFR = mctFR + ctFR;
end
mctFR = fillmissing(mctFR,'nearest','EndValues','extrap');
mctFR = fillmissing(mctFR','nearest','EndValues','extrap');
mctFR = mctFR'/5;
figure(200);clf;imagesc(tTrials,tTrials,mctFR,[0 0.5]); set(gca,'Ydir','normal');colorbar;
%% Active Cells non-NaN (For Firing Rate and finding the tuning of different neurons with respect to time, distance, and speed)
variable_combs = {'FR_time','FR_dist','FR_speed'};
% Factors, CN, 
clc
active_cells = {};
out = outT;
ap = 1;
avar = [];
for an = 1:5
    anvar = [];
    for cn = 1:3
        for ap = 1:2
            AC = out.active_cells{an,cn,ap};
            for vn = 1%:length(variable_combs)
                tmvT = met_valsT{vn};
                tmvD = met_valsD{vn};
                tPCT = ~isnan(tmvT{an,cn,ap}.PC(:,1));
                tPCD = ~isnan(tmvD{an,cn,ap}.PC(:,1));
                % tPCT = AC;
                active_cells{an,cn,ap} = tPCT;
                NR = [tPCT tPCD];
                if sum(tPCT) ~= sum(tPCD)
                    disp('Differences Found')
                else
                    fResp = (sum(tPCT)/length(tPCT))*100;
                end
                anvar = [anvar fResp];
            end
        end
    end
    avar = [avar;anvar];
end
descriptiveStatistics(avar(:)); %(mean ± sem: 0.281 ± 0.022, range: -0.957, 0.999, median: 0.485) concatenate MI (mean ± sem: 0.746 ± 0.051, range: 0.187, 1.998, median: 0.593)
% (Intercept):PT [F(2,8) = 116.50, p < 0.001, η2 = .84] <--
fac_names = {'CN','AP'}; fac_levels = [3,2];
[within,dvn,xlabels,awithinD] = make_within_table(fac_names,fac_levels);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)
% Run Reduced Subset RM-ANOVA
raR = RMA_subset(ra,'AP');
descriptiveStatistics(raR{1}.ds)
descriptiveStatistics(raR{2}.ds)

% Sanity Check
% for an = 1:5
%     for cn = 1:3
%         for ap = 1:2
%             AC = out.active_cells{an,cn,ap};
%             if AC ~= active_cells{an,cn,ap}
%                 disp('Differences_Found');
%             end
%         end
%     end
% end
%% Trial-Wise --> Active Cells FR > 0

variable_combs = {'FR_time','FR_dist','FR_speed'};
% Factors, CN, 
clc
active_cells = {};
out = outT;
ap = 1;
avar = [];
for an = 1:5
    anvar = [];
    for cn = 1:3
        for ap = 1:2
            AC = out.active_cells_trials{an,cn,ap};
            fResp = (sum(AC)/size(AC,1))*100;
            anvar = [anvar fResp];
        end
    end
    avar = [avar;anvar];
end
descriptiveStatistics(avar(:)); %(mean ± sem: 0.281 ± 0.022, range: -0.957, 0.999, median: 0.485) concatenate MI (mean ± sem: 0.746 ± 0.051, range: 0.187, 1.998, median: 0.593)
% (Intercept):PT [F(2,8) = 116.50, p < 0.001, η2 = .84] <--
fac_names = {'CN','AP','TN'}; fac_levels = [3,2,10];
[within,dvn,xlabels,awithinD] = make_within_table(fac_names,fac_levels);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)
% Run Reduced Subset RM-ANOVA
raR = RMA_subset(ra,'AP');
descriptiveStatistics(raR{1}.ds)
descriptiveStatistics(raR{2}.ds)

%% Trial-Wise --> Active Cells FR > 0 ... comparison of percentage distributions
%

variable_combs = {'FR_time','FR_dist','FR_speed'};
% Factors, CN, 
clc
active_cells = {};
out = outT;
ap = 1;
avar = [];
for an = 1:5
    anvar = [];
    for cn = 1:3
        for ap = 1:2
            AC = out.active_cells_trials{an,cn,ap};
            fResp = hist(sum(AC,2),[0:10]);
            fResp = 100*fResp/sum(fResp);
            if length(fResp) < 11
                err
            end
            anvar = [anvar fResp];
        end
    end
    avar = [avar;anvar];
end
descriptiveStatistics(avar(:)); %(mean ± sem: 0.281 ± 0.022, range: -0.957, 0.999, median: 0.485) concatenate MI (mean ± sem: 0.746 ± 0.051, range: 0.187, 1.998, median: 0.593)
% (Intercept):PT [F(2,8) = 116.50, p < 0.001, η2 = .84] <--
fac_names = {'CN','AP','NT'}; fac_levels = [3,2,11];
[within,dvn,xlabels,awithinD] = make_within_table(fac_names,fac_levels);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)
% Run Reduced Subset RM-ANOVA
raR = RMA_subset(ra,'AP');
descriptiveStatistics(raR{1}.ds)
descriptiveStatistics(raR{2}.ds)

%% Response Fidelity --> Active Cells FR > 0
variable_combs = {'FR_time','FR_dist','FR_speed'};
% Factors, CN, 
clc
active_cells = {};
out = outT;
ap = 1;
avar = [];
for an = 1:5
    anvar = [];
    for cn = 1:3
        for ap = 1:2
            AC = out.active_cells_trials{an,cn,ap};
            fResp = (sum(AC,2)/10)*100;
            anvar = [anvar fResp];
        end
    end
    avar = [avar;anvar];
end
descriptiveStatistics(avar(:)); %(mean ± sem: 0.281 ± 0.022, range: -0.957, 0.999, median: 0.485) concatenate MI (mean ± sem: 0.746 ± 0.051, range: 0.187, 1.998, median: 0.593)
% (Intercept):PT [F(2,8) = 116.50, p < 0.001, η2 = .84] <--
fac_names = {'CN','AP'}; fac_levels = [3,2];
[within,dvn,xlabels,awithinD] = make_within_table(fac_names,fac_levels);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)
% Run Reduced Subset RM-ANOVA
raR = RMA_subset(ra,'AP');
descriptiveStatistics(raR{1}.ds)
descriptiveStatistics(raR{2}.ds)

%% low, medium, high Response Fidelity --> Active Cells FR > 0 but categorized as low, medium, high
variable_combs = {'FR_time','FR_dist','FR_speed'};
% Factors, CN, 
clc
active_cells = {};
out = outT;
ap = 1;
avar = [];
for an = 1:5
    anvar = [];
    for cn = 1:3
        for ap = 1:2
            AC = out.active_cells_trials{an,cn,ap};
            AAC1 = sum(AC,2) > 0 & sum(AC,2) <= 3;
            AAC2 = sum(AC,2) > 3 & sum(AC,2) <= 6;
            AAC3 = sum(AC,2) > 6 & sum(AC,2) < 10;
            fResp = [((sum(AAC1)/length(AAC1))*100) ((sum(AAC2)/length(AAC2))*100) ((sum(AAC3)/length(AAC3))*100)] ;
            anvar = [anvar fResp];
        end
    end
    avar = [avar;anvar];
end
descriptiveStatistics(avar(:)); %(mean ± sem: 0.281 ± 0.022, range: -0.957, 0.999, median: 0.485) concatenate MI (mean ± sem: 0.746 ± 0.051, range: 0.187, 1.998, median: 0.593)
% (Intercept):PT [F(2,8) = 116.50, p < 0.001, η2 = .84] <--
fac_names = {'CN','AP','RF'}; fac_levels = [3,2,3];
[within,dvn,xlabels,awithinD] = make_within_table(fac_names,fac_levels);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)
% Run Reduced Subset RM-ANOVA
raR = RMA_subset(ra,'AP');
descriptiveStatistics(raR{1}.ds)
descriptiveStatistics(raR{2}.ds)

%% Conjunctive and Complementary populations
variable_combs = {'FR_time','FR_dist','FR_speed'};
% Factors, CN, 
clc
ap = 1;
avar = [];
tmvT = met_valsT{1};
tmvD = met_valsD{vn};
for an = 1:5
    for cn = 1:3
        tPCT1 = ~isnan(tmvT{an,cn,1}.PC(:,1));
        tPCT2 = ~isnan(tmvT{an,cn,2}.PC(:,1));
        conj(an,cn) = 100*(sum((tPCT1 & tPCT2))/length(tPCT1));
        comp1(an,cn) = 100*(sum((tPCT1 & ~tPCT2))/length(tPCT1));
        comp2(an,cn) = 100*(sum((~tPCT1 & tPCT2))/length(tPCT1));
    end
end
avar = [conj comp1 comp2];
fac_names = {'CT','CN'}; fac_levels = [3,3];
[within,dvn,xlabels,awithinD] = make_within_table(fac_names,fac_levels);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)

%% Air-On --> Get the metrics and run stats FR --> time, distance, speed
variable_combs = {'FR_time','FR_dist','FR_speed'};
% metric = 'MI';
metric = 'PC';
RV = 'perc';
% RV = 'rfid';
% RV = 'mval';


avar = [];
for an = 1:5
    anvar = [];
    for cn = 1:3
        for ap = 1:2
            AC = outT.active_cells_trials{an,cn,ap}; mAC = outT.active_cells{an,cn,ap}; RF = sum(AC,2);
            for bn = 1:2
                if bn == 1
                    out = outT; MV = met_valsT; 
                else
                    out = outD; MV = met_valsD; 
                end
                selcells = mAC;% & (RF>0);
                if strcmp(metric,'PC')
                    idx = 3; pvals = [MV{1}{an,cn,ap}.PC(:,idx) MV{2}{an,cn,ap}.PC(:,idx) MV{3}{an,cn,ap}.PC(:,idx)];
                    idx = 2; zvals = [MV{1}{an,cn,ap}.PC(:,idx) MV{2}{an,cn,ap}.PC(:,idx) MV{3}{an,cn,ap}.PC(:,idx)];
                else
                    idx = 3; pvals = [MV{1}{an,cn,ap}.MI(:,idx) MV{2}{an,cn,ap}.MI(:,idx) MV{3}{an,cn,ap}.MI(:,idx)];
                    idx = 2; zvals = [MV{1}{an,cn,ap}.MI(:,idx) MV{2}{an,cn,ap}.MI(:,idx) MV{3}{an,cn,ap}.MI(:,idx)];
                end
                % criteria for finding the types of cells
                tCells = ((pvals < 0.05) * [4; 2; 1]) == 4; dCells = ((pvals < 0.05) * [4; 2; 1]) == 2; sCells = ((pvals < 0.05) * [4; 2; 1]) == 1;
                time_cells{an,cn,ap,bn} = tCells; distance_cells{an,cn,ap,bn} = dCells; speed_cells{an,cn,ap,bn} = sCells;
                if strcmp(RV,'perc')
                    anvar = [anvar sum(tCells & selcells)/length(AC) sum(dCells & selcells)/length(AC) sum(sCells & selcells)/length(AC)];
                end
                if strcmp(RV,'mval')
                    anvar = [anvar nanmean(zvals(tCells&selcells,1)) nanmean(zvals(dCells&selcells,2)) nanmean(zvals(sCells&selcells,3))];
                end
                if strcmp(RV,'rfid')
                    anvar = [anvar mean(RF(tCells&selcells,1)) mean(RF(dCells&selcells,1)) mean(RF(sCells&selcells,1))];
                end
                if sum(isnan(anvar)) > 0
                    [an cn ap bn rfi]
                end

            end
        end
    end
    if strcmp(RV,'perc')
        avar = [avar;100*anvar];
    else
        avar = [avar;anvar];
    end
end

fac_names = {'CN','AP','BT','TT'}; fac_levels = [3,2,2,3];
% fac_names = {'CN','AP','BT','RF','CT'}; fac_levels = [3,2,2,3,4];
% fac_names = {'CN','AP','BT','CT'}; fac_levels = [3,2,2,7];
[within,dvn,xlabels,awithinD] = make_within_table(fac_names,fac_levels);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
clc
print_for_manuscript(ra)
% 'NoT' 'speed'  'dist'  'dist-speed'  'time'  'time-speed'  'time-dist'  'TDS'
% separation by AP - air phase
clc
ra_AP = RMA_subset(ra,'AP');
%% visualizing the results in the previous section
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova([],ra_AP{2},{'BT:TT','hsd',0.05},[1 1.75],tcolors,[mY MY ysp ystf ysigf],mData);
%% separation AP by MT - air phase by metric type
clc
ra_AP1_BT = RMA_bonferroni(ra_AP{2},'TT');
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova([],ra_AP1_BT{3},{'BT','hsd',0.05},[1 1.75],tcolors,[mY MY ysp ystf ysigf],mData);
% ra_AP1_BT2_CN = RMA_subset(ra_AP1_BT{2},'CN');
%% Overlap indices
ti_cells = []; di_cells = []; sp_cells = [];
bn = 1;
for an = 1:5
    cc = 1;
    for cn = 1:3
        for ap = 1:2
            ti_cells{an,cc} = time_cells{an,cn,ap,bn};
            di_cells{an,cc} = distance_cells{an,cn,ap,bn};
            sp_cells{an,cc} = speed_cells{an,cn,ap,bn};
            cc = cc + 1;
        end
    end
end

allresp = [ti_cells di_cells sp_cells];

% [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(allresp,0.5,0.05);\
hf = figure(100);clf;clc
txl = {'C3-AOn-T','C3-AOff-T','C4-AOn-T','C4-AOff-T','C5-AOn-T','C5-AOff-T','C3-AOn-D','C3-AOff-D','C4-AOn-D','C4-AOff-D','C5-AOn-D','C5-AOff-D','C3-AOn-S','C3-AOff-S','C4-AOn-S','C4-AOff-S','C5-AOn-S','C5-AOff-S'};
[mOI,semOI] = heatmap_conj_comp(gca,allresp,0,{1:18,txl,1});
mOI1 = mOI;
    mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1);
    tree = linkage(Di);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','right','ColorThreshold',1.305);
    hf = gcf;
    % set(hf,'Position',[7 3 1.25 2]);
    set(H,'linewidth',1);
    set(gca,'yticklabels',txl(TC));ytickangle(30);
    format_axes(gca);
    hx = xlabel('Eucledian Distance');%changePosition(hx,[-0.051 0 0]);
    changePosition(gca,[0 0.0 0.05 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);
%% transition of time, distance, and speed cells across air-on and air-off phases
ti_cellsC = []; di_cellsC = []; sp_cellsC = [];
bn = 1;
for an = 1:5
    cc = 1; ti_cells = []; di_cells = []; sp_cells = [];
    for cn = 1%:3
        for ap = 1:2
            ti_cells(:,cc) = time_cells{an,cn,ap,bn};
            di_cells(:,cc) = distance_cells{an,cn,ap,bn};
            sp_cells(:,cc) = speed_cells{an,cn,ap,bn};
            cc = cc + 1;
        end
    end
    ti_cellsC{an,1} = ti_cells;
    di_cellsC{an,1} = di_cells;
    sp_cellsC{an,1} = sp_cells;
end


%% average speed and other speed characteristics
% for total time, distance use the variable_combs = {'time'} or distance
% and then use the max_fun. you may also change ap to 1 or 2 depending upon
% whether we want air on phase or air off phase
out = outT;
mean_fun = @(x) mean(x);
std_fun = @(x) std(x, 0);  % Standard deviation (unbiased)
cv_fun = @(x) std(x, 0) ./ mean(x);  % Coefficient of Variation (CV)
skew_fun = @(x) skewness(x);  % Skewness
kurt_fun = @(x) kurtosis(x);  % Kurtosis
max_fun = @(x) max(x); 
min_fun = @(x) min(x);
% latency_fun = @(x, t) t(find(x > 0, 1, 'first'));  
% rlatency_fun = @(x, t) t(find(x > 0, 1, 'last'));  

variable_combs = {'speed'};
% variable_combs = {'time'};
avar = [];
for an = 1:5
    anvar = [];
    for cn = 1:3
        for ap = 1:2
            timecc = out.otime{an,cn,ap}; distcc = out.odist{an,cn,ap}; speedcc = out.ospeed{an,cn,ap}; trialcc = out.otrial{an,cn,ap};
            for vn = 1:length(variable_combs)
                [an cn ap vn]
                var_name = variable_combs{vn};
                idx_us = strfind(var_name,'_');
                if length(variable_combs) == 1
                    cmdTxt = sprintf('var1v = %scc;',var_name);eval(cmdTxt);
                end
                thisvar = (accumarray(trialcc,var1v,[],max_fun))';
                anvar = [anvar thisvar];% outD.trial_metrics(:,idx)'];
                n = 0;
            end
        end
    end
    avar = [avar;anvar];
end
%
n = 0;
clc
[within,dvn,xlabels,awithinD] = make_within_table({'CN','AP','TN'},[3,2,10]);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)
%%
% the following or the previous one to this will give error depending on
% the value of ap if it is only air on phase or air off phase
[within,dvn,xlabels,awithinD] = make_within_table({'CN','TN'},[3,10]);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)
%% total time, distance and other things 

mean_fun = @(x) mean(x);
std_fun = @(x) std(x, 0);  % Standard deviation (unbiased)
cv_fun = @(x) std(x, 0) ./ mean(x);  % Coefficient of Variation (CV)
skew_fun = @(x) skewness(x);  % Skewness
kurt_fun = @(x) kurtosis(x);  % Kurtosis
max_fun = @(x) max(x); 
min_fun = @(x) min(x);
latency_fun = @(x, t) t(find(x > 0, 1, 'first'));  % movement_latency = accumarray(trialcc, speed, [], @(x) latency_fun(x, time));
rlatency_fun = @(x, t) t(find(x > 0, 1, 'last'));  % rest_latency = accumarray(trialcc, speed, [], @(x) latency_fun(x, time));
sp_thr = 0;
mon_fun = @(x) length(find_rising_edge(x > sp_thr,0.1,500));
moff_fun = @(x) length(find_falling_edge(x > sp_thr,-0.1,500));
latency_fun1 = @(x, t) t(find_first_rising_edge(x > sp_thr,0.1,0));  % movement_latency = accumarray(trialcc, speed, [], @(x) latency_fun(x, time));
rlatency_fun1 = @(x, t) t(find_first_falling_edge(x > sp_thr,-0.1,0));

variable_combs = {'speed'};
avar = [];
for an = 1:5
    anvar = [];
    for cn = 1:3
        for ap = 2
            timecc = o_atimecc{an,cn,ap}; distcc = o_adistcc{an,cn,ap}; speedcc = o_aspeedcc{an,cn,ap}; trialcc = o_atrialcc{an,cn,ap};
            for vn = 1:length(variable_combs)
                [an cn ap vn]
                var_name = variable_combs{vn};
                idx_us = strfind(var_name,'_');
                if length(variable_combs) == 1
                    cmdTxt = sprintf('var1v = %scc;',var_name);eval(cmdTxt);
                end
                grouped_var1v = accumarray(trialcc, var1v, [], @(x) {x}); grouped_timecc = accumarray(trialcc, timecc, [], @(x) {x});
                thisvar = arrayfun(@(i) rlatency_fun1(grouped_var1v{i}, grouped_timecc{i}), 1:length(grouped_var1v));
                % thisvar = (accumarray(trialcc,var1v,[],moff_fun))';
                anvar = [anvar thisvar];
                n = 0;
            end
        end
    end
    avar = [avar;anvar];
end
descriptiveStatistics(avar(:));

%%
clc
raR = RMA_bonferroni(ra,1);
%%
clc
[within,dvn,xlabels,awithinD] = make_within_table({'CN','TN'},[3,10]);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)


%%
n = 0;
clc
fac_names = {'CN','TN'}; fac_levels = [3,10];
fac_names = {'CN','PT','TN'}; fac_levels = [3,3,10];
% fac_names = {'CN','AP','TN'}; fac_levels = [3,2,10];
fac_names = {'CN','AP','PT','TN'}; fac_levels = [3,2,3,10];
% fac_names = {'CN','AP','PT'}; fac_levels = [3,2,3];
% fac_names = {'CN','PT'}; fac_levels = [3,3];
[within,dvn,xlabels,awithinD] = make_within_table(fac_names,fac_levels);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)
%%
ra1 = ra;
ra2 = ra;
%%
clc
raR = RMA_bonferroni(ra,2);
%%
clc
raR = RMA_subset(ra,'AP');
%%
raRR = RMA_bonferroni(raR{1},'PT');
% [xdata,mVar,semVar,combs,p,h,nB] = RMA_get_multcompare(raR{1},{'TN','hsd',0.05},[1 1]);
%%
raRB = RMA_bonferroni(raR.ras{2},3);
[xdata,mVar,semVar,combs,p,h,nB] = RMA_get_multcompare(raRB{1},{'TN','hsd',0.05},[1 1]);
%%
n = 0;
clc
[within,dvn,xlabels,awithinD] = make_within_table({'CN','TN','PT'},[3,10,3]);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)
%% speed time and dist glm
out = outD;
avar = [];
for an = 1:5
    anvar = [];
    for cn = 1:3
        for ap = 1:2
            timecc = out.atimecc{an,cn,ap}; distcc = out.adistcc{an,cn,ap}; speedcc = out.aspeedcc{an,cn,ap}; trialcc = out.atrialcc{an,cn,ap};
            Y = speedcc; % Firing rate for current neuron (timepoints x 1)
            X = [timecc,distcc]; % Predictors (time, distance, speed)
            
            % Fit the GLM (here we use a linear model, but you can modify it for other GLM families)
            mdl = fitglm(X, Y, 'linear');
            aglm{an,cn,ap} = mdl;
            anvar = [anvar mdl.Coefficients.pValue(1) mdl.Coefficients.pValue(2)];
        end
    end
    avar = [avar;anvar];
end
clc
[within,dvn,xlabels,awithinD] = make_within_table({'CN','AP','DT'},[3,2,2]);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)
%% GLM
aglm = {};
for an = 1:5
    for cn = 1:3
        for ap = 1:2
            [an cn ap]
            timecc = atimecc{an,cn,ap}; distcc = adistcc{an,cn,ap}; speedcc = aspeedcc{an,cn,ap}; FRcc = aFRcc{an,cn,ap};
            % Assuming `timecc`, `distcc`, `speedcc`, and `FRcc` are already defined

            numNeurons = size(FRcc, 2); % Number of neurons
            results = cell(numNeurons, 1); % Store GLM results for each neuron
            parfor neuron_idx = 1:numNeurons
                % Prepare the data for the GLM
                Y = FRcc(:, neuron_idx); % Firing rate for current neuron (timepoints x 1)
                X = [timecc, distcc, speedcc]; % Predictors (time, distance, speed)
                
                % Fit the GLM (here we use a linear model, but you can modify it for other GLM families)
                mdl = fitglm(X, Y, 'linear'); % You can change 'linear' to other distributions as needed
                
                % Store the model results for the current neuron
                results{neuron_idx} = mdl;
            end
            aglm{an,cn,ap} = results;
        end
    end
end
%%
% Assuming aglm{an, cn, ap} contains the GLM results for each animal, configuration, and phase
numAnimals = size(aglm, 1); % Number of animals
numConfigurations = size(aglm, 2); % Number of configurations
numAirPhases = size(aglm, 3); % Number of air-phases
% numNeurons = size(FRcc{1,1,1}, 2); % Number of neurons (assuming FRcc has neurons in columns)

% Initialize a cell to store the tuning results (for each animal, configuration, and phase)
tuningResults = cell(numAnimals, numConfigurations, numAirPhases);

% Loop through each animal, configuration, and air-phase
for an = 1:numAnimals
    for cn = 1:numConfigurations
        for ap = 1:numAirPhases
            % Get the GLM results for this animal, configuration, and air-phase
            results = aglm{an, cn, ap};
            numNeurons = length(results);
            % Initialize an array to store the tuning category for each neuron
            neuronTuning = cell(numNeurons, 1);
            
            % Loop through each neuron
            for neuron_idx = 1:numNeurons
                % Get the GLM model for this neuron
                mdl = results{neuron_idx};
                
                % Extract the coefficients for time, distance, and speed
                coefficients = mdl.Coefficients.Estimate(2:4); % Time, Distance, Speed (assuming first coefficient is intercept)
                pvals = mdl.Coefficients.pValue(2:4); % Time, Distance, Speed (assuming first coefficient is intercept)
                sigpval = pvals' < 0.05;
                sigpvaldec = (sigpval(1) * 4 + sigpval(2) * 2 + sigpval(3) * 1);
                tuningTypes = {'Untuned','Speed','Distance','Dist-Speed','Time','Time-Speed','Time-Dist','Time-Dist-Speed'};
                thistuning = tuningTypes{sigpvaldec + 1};
                neuronTuning{neuron_idx} = thistuning;
            end
            % Store the tuning results for this animal, configuration, and air-phase
            tuningResults{an, cn, ap} = neuronTuning;
        end
    end
end
%%
% Initialize arrays to store percentages (or fractions) for time, distance, and speed tuning
timeTunedPercent = zeros(numAnimals, numConfigurations, numAirPhases);
distanceTunedPercent = zeros(numAnimals, numConfigurations, numAirPhases);
speedTunedPercent = zeros(numAnimals, numConfigurations, numAirPhases);

% Loop through the tuning results and calculate the percentages (or fractions)
avar = [];
for an = 1:numAnimals
    anvar = [];
    for cn = 1:numConfigurations
        for ap = 1:numAirPhases
            % Get the tuning results for this animal, configuration, and phase
            neuronTuning = tuningResults{an, cn, ap};
            numNeurons = length(neuronTuning);
            % tuningTypes = {'Untuned','Speed','Distance','Dist-Speed','Time','Time-Speed','Time-Dist','Time-Dist-Speed'};
            numtuned = [];
            for tii = 1:length(tuningTypes)
                fractuned(1,tii) = sum(strcmp(neuronTuning, tuningTypes{tii}))/numNeurons;
            end
            anvar = [anvar fractuned([2 3 5])];
        end
    end
    avar = [avar;anvar];
end

