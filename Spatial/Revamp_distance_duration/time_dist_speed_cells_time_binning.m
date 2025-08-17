%% first find the variables time_cells, distance_cells, and speed_cells in the final_time_dist_bin...m file
variable_combs = {'FR_time','FR_dist','FR_speed'};
metric = 'MI';
metric = 'PC';
clear time_cells distance_cells speed_cells
clear t_cells_T d_cells_T s_cells_T td_cells_T ds_cells_T ts_cells_T tds_cells_T ntds_cells_T
clear t_cells_I d_cells_I s_cells_I td_cells_I ds_cells_I ts_cells_I tds_cells_I ntds_cells_I
all_cells = {};
si = [Ar_t_T Ar_i_T Ar_t_D Ar_i_D ArL_t_T ArL_i_T ArL_t_D ArL_i_D Ars_t_T Ars_i_T Ars_t_D Ars_i_D]; propsPL = get_props_Rs(o.Rs(:,si),30);
si = [Ar_t_T Ar_i_T ArL_t_T ArL_i_T Ars_t_T Ars_i_T]; si_cn_ap = [[1 1 2 2 3 3];[1 2 1 2 1 2]];
% si = [Ar_t_D Ar_i_D ArL_t_D ArL_i_D Ars_t_D Ars_i_D]; 
props = get_props_Rs(o.Rs(:,si),30);
propsT = get_props_new(outT,met_valsT,props,si_cn_ap);
propsD = get_props_new(outD,met_valsD,props,si_cn_ap);
% [SinT,MixT,AllT] = get_the_pops(propsT,propsD);
%% untuned or non-responsive cells
pop_names = {'R','NR'};
all_cellsnew = []; all_cellsnew_vals = [];
for ii = 1:length(pop_names)
    cmdTxt = sprintf('all_cellsnew = [all_cellsnew propsT.newPC.cells_%s propsD.newPC.cells_%s propsT.newMI.cells_%s propsD.newMI.cells_%s];',pop_names{ii},pop_names{ii},pop_names{ii},pop_names{ii}); eval(cmdTxt);
    % cmdTxt = sprintf('all_cellsnew_vals = [all_cellsnew propsT.newPC.cells_%s_zvals propsD.newPC.cells_%s propsT.newMI.cells_%s propsD.newMI.cells_%s];',pop_names{ii},pop_names{ii},pop_names{ii},pop_names{ii}); eval(cmdTxt);
end

avar = exec_fun_on_cell_mat(all_cellsnew,'percent');
fac_names = {'PoT','MT','BT','CN','AP'}; fac_levels = [2,2,2,3,2];
[within,dvn,xlabels,awithinD] = make_within_table(fac_names,fac_levels);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
clc
print_for_manuscript(ra)
%%
clc
raR = RMA_bonferroni(ra,'PoT');

%%
clc
raRR1 = RMA_bonferroni(raR{1},'MT');
%%
clc
raRR2 = RMA_bonferroni(raR{2},'MT');


%% this is basically for singularly_tuned mixed_tuned etc., but I am just checking for other to use the same code
pop_names = {'singularly_tuned','mixed_tuned','all_tuned'};
% % pop_names = {'NR','speed','dist','dist_speed','time','time_speed','time_dist','time_dist_speed','singularly_tuned','mixed_tuned','all_tuned'};
% pop_names = {'time','dist','speed','time_speed','time_dist','dist_speed'};
all_cellsnew = [];
for ii = 1:length(pop_names)
    cmdTxt = sprintf('all_cellsnew = [all_cellsnew propsT.newPC.cells_%s propsD.newPC.cells_%s propsT.newMI.cells_%s propsD.newMI.cells_%s];',pop_names{ii},pop_names{ii},pop_names{ii},pop_names{ii}); eval(cmdTxt);
end

avar = exec_fun_on_cell_mat(all_cellsnew,'percent');
fac_names = {'PoT','MT','BT','CN','AP'}; fac_levels = [length(pop_names),2,2,3,2];
[within,dvn,xlabels,awithinD] = make_within_table(fac_names,fac_levels);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
clc
print_for_manuscript(ra)
%%
clc
raR = RMA_bonferroni(ra,'PoT');
% RMA_visualize(raR{1},{'MT:BT','hsd',0.05});
%%
clc
raRR1 = RMA_bonferroni(raR{1},'MT');
%%
clc
raRR2 = RMA_bonferroni(raR{2},'MT');
%%
tcolors = repmat(mData.dcolors(1:10),1,3); MY = 100; ysp = 5; mY = 0; ystf = 5; ysigf = 0.025;titletxt = ''; ylabeltxt = {'Cells (%)'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova([],raR{1},{'MT:BT','hsd',0.05},[1 1.75],tcolors,[mY MY ysp ystf ysigf],mData);


%%
avar = cell_list_op_percent(props.good_FR,props.good_FR);

fac_names = {'CN','AP'}; fac_levels = [3,2];
[within,dvn,xlabels,awithinD] = make_within_table(fac_names,fac_levels);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
clc
print_for_manuscript(ra)

%%
% {'NR','speed','dist','dist_speed','time','time_speed','time_dist','time_dist_speed','singularly_tuned','mixed_tuned','all_tuned'};
pop_names = {'time','dist','speed'};
% pop_names = {'time','dist','speed','time_speed','time_dist','dist_speed'};
all_cellsnew = [];  all_cellsnew_hm = []; txl_all = {}; cnti = 1;
all_cellsnew_zvals = []; all_cellsnew_rs = []; MTs = {'PC','MI'};
all_cellsnew_peakvalsT = []; all_cellsnew_peakvalsD = [];
all_cellsnew_GWT = []; all_cellsnew_GWD = [];
for ii = 1:length(pop_names)
    if ii == 1 || ii == 3
    cmdTxt = sprintf('all_cellsnew = [all_cellsnew propsT.newPC.cells_%s propsT.newMI.cells_%s];',pop_names{ii},pop_names{ii}); eval(cmdTxt);
    cmdTxt = sprintf('all_cellsnew_hm = [all_cellsnew_hm propsT.newMI.cells_%s];',pop_names{ii}); eval(cmdTxt);
    cmdTxt = sprintf('all_cellsnew_zvals = [all_cellsnew_zvals get_zvals(propsT.newPC.cells_%s_zvals,ii) get_zvals(propsT.newMI.cells_%s_zvals,ii)];',pop_names{ii},pop_names{ii}); eval(cmdTxt);
    else
    cmdTxt = sprintf('all_cellsnew = [all_cellsnew propsD.newPC.cells_%s propsD.newMI.cells_%s];',pop_names{ii},pop_names{ii}); eval(cmdTxt);
    cmdTxt = sprintf('all_cellsnew_hm = [all_cellsnew_hm propsD.newMI.cells_%s];',pop_names{ii}); eval(cmdTxt);
    cmdTxt = sprintf('all_cellsnew_zvals = [all_cellsnew_zvals get_zvals(propsD.newPC.cells_%s_zvals,ii) get_zvals(propsD.newMI.cells_%s_zvals,ii)];',pop_names{ii},pop_names{ii}); eval(cmdTxt);
    end
    si = [Ar_t_T Ar_i_T ArL_t_T ArL_i_T Ars_t_T Ars_i_T]; propsPL = get_props_Rs(o.Rs(:,si),30);
    all_cellsnew_peakvalsT = [all_cellsnew_peakvalsT propsPL.peak_locations propsPL.peak_locations];
    all_cellsnew_GWT = [all_cellsnew_GWT propsPL.PWs propsPL.PWs];
    si = [Ar_t_D Ar_i_D ArL_t_D ArL_i_D Ars_t_D Ars_i_D]; propsPL = get_props_Rs(o.Rs(:,si),30);
    all_cellsnew_peakvalsD = [all_cellsnew_peakvalsD propsPL.peak_locations propsPL.peak_locations];
    all_cellsnew_GWD = [all_cellsnew_GWD propsPL.PWs propsPL.PWs];
    for cn = 1:3
        for ap = 1:2
            for mti = 2
                txl_all{cnti} = sprintf('C%d-A%d-%s-%s',cn+2,ap,MTs{mti},upper(pop_names{ii}(1)));
                cnti = cnti + 1;
            end
        end
    end
end

for ii = 1:length(pop_names)
    if ii == 1 || ii == 3
    cmdTxt = sprintf('all_cellsnew = [all_cellsnew propsD.newPC.cells_%s propsD.newMI.cells_%s];',pop_names{ii},pop_names{ii}); eval(cmdTxt);
    cmdTxt = sprintf('all_cellsnew_hm = [all_cellsnew_hm propsD.newPC.cells_%s];',pop_names{ii}); eval(cmdTxt);
    cmdTxt = sprintf('all_cellsnew_zvals = [all_cellsnew_zvals get_zvals(propsD.newPC.cells_%s_zvals,ii) get_zvals(propsD.newMI.cells_%s_zvals,ii)];',pop_names{ii},pop_names{ii}); eval(cmdTxt);
    else
    cmdTxt = sprintf('all_cellsnew = [all_cellsnew propsT.newPC.cells_%s propsT.newMI.cells_%s];',pop_names{ii},pop_names{ii}); eval(cmdTxt);
    cmdTxt = sprintf('all_cellsnew_hm = [all_cellsnew_hm propsT.newPC.cells_%s];',pop_names{ii}); eval(cmdTxt);
    cmdTxt = sprintf('all_cellsnew_zvals = [all_cellsnew_zvals get_zvals(propsT.newPC.cells_%s_zvals,ii) get_zvals(propsT.newMI.cells_%s_zvals,ii)];',pop_names{ii},pop_names{ii}); eval(cmdTxt);
    end
    si = [Ar_t_T Ar_i_T ArL_t_T ArL_i_T Ars_t_T Ars_i_T]; propsPL = get_props_Rs(o.Rs(:,si),30);
    all_cellsnew_peakvalsT = [all_cellsnew_peakvalsT propsPL.peak_locations propsPL.peak_locations];
    all_cellsnew_GWT = [all_cellsnew_GWT propsPL.PWs propsPL.PWs];
    si = [Ar_t_D Ar_i_D ArL_t_D ArL_i_D Ars_t_D Ars_i_D]; propsPL = get_props_Rs(o.Rs(:,si),30);
    all_cellsnew_peakvalsD = [all_cellsnew_peakvalsD propsPL.peak_locations propsPL.peak_locations];
    all_cellsnew_GWD = [all_cellsnew_GWD propsPL.PWs propsPL.PWs];
    for cn = 1:3
        for ap = 1:2
            for mti = 1
                txl_all{cnti} = sprintf('C%d-A%d-%s-%s',cn+2,ap,MTs{mti},upper(pop_names{ii}(1)));
                cnti = cnti + 1;
            end
        end
    end
end

avar = exec_fun_on_cell_mat(all_cellsnew,'percent');
avar = exec_fun_on_cell_mat(all_cellsnew_zvals,'nanmean',all_cellsnew);
% avar = exec_fun_on_cell_mat(all_cellsnew_GWT,'nanmean',all_cellsnew);
% avar = exec_fun_on_cell_mat(all_cellsnew_peakvalsT,'nanmean',all_cellsnew);
fac_names = {'FM','PoT','MT','CN','AP'}; fac_levels = [2,length(pop_names),2,3,2];
% fac_names = {'MT','CN','AP'}; fac_levels = [2,3,2];
[within,dvn,xlabels,awithinD] = make_within_table(fac_names,fac_levels);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
clc
print_for_manuscript(ra)
%%
for cid = 1:3
rows = awithinD(:,3) == 1 & awithinD(:,1) == 1 & awithinD(:,2) == cid;
all_cellsnew_PC = all_cellsnew(:,find(rows));
for an = 1:5
    all_cellsnew_PC_mat_time{an,cid} = cell2mat(all_cellsnew_PC(an,:));
end
end
%%
clc
raS = RMA_subset(ra,'MT');
% raS = RMA_bonferroni(ra,'FM');
%%
clc
raR = RMA_bonferroni(raS{1},'FM');

%%
clc
raR1 = RMA_bonferroni(raR{1},'MT');

%%
clc
raR = RMA_bonferroni(raS{2},'PoT');

%% visualizing the results in the previous section
tcolors = repmat(mData.dcolors(1:10),1,3); MY = 5; ysp = 1; mY = 0; ystf = 2; ysigf = 0.025;titletxt = ''; ylabeltxt = {'Cells (%)'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova([],raS{1},{'MT:AP','hsd',0.05},[1 1.75],tcolors,[mY MY ysp ystf ysigf],mData);

%%
clc
raRR = RMA_bonferroni(raR{3},'AP');

%% visualizing the results in the previous section
tcolors = repmat(mData.dcolors(1:10),1,3); MY = 10; ysp = 1; mY = 0; ystf = 2; ysigf = 0.025;titletxt = ''; ylabeltxt = {'Cells (%)'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova([],ra,{'PoT:BT','hsd',0.05},[1 1.75],tcolors,[mY MY ysp ystf ysigf],mData);
[EMs,GS] = RMA_get_EM_GS(raR{2}.rm,{'MT','AP'});

%%
% {'NR','speed','dist','dist_speed','time','time_speed','time_dist','time_dist_speed','singularly_tuned','mixed_tuned','all_tuned'};
pop_names = {'time','dist','speed'};
all_cellsnew = []; all_cellsnew_zvals = [];
for ii = 1:length(pop_names)
    cmdTxt = sprintf('all_cellsnew = [all_cellsnew propsT.newPC.cells_%s propsT.newMI.cells_%s propsD.newPC.cells_%s propsD.newMI.cells_%s];',pop_names{ii},pop_names{ii},pop_names{ii},pop_names{ii}); eval(cmdTxt);
    cmdTxt = sprintf(['all_cellsnew_zvals = [all_cellsnew_zvals get_zvals(propsT.newPC.cells_%s_zvals,ii) get_zvals(propsT.newMI.cells_%s_zvals,ii) ' ...
        'get_zvals(propsD.newPC.cells_%s_zvals,ii) get_zvals(propsD.newMI.cells_%s_zvals,ii)];'],pop_names{ii},pop_names{ii},pop_names{ii},pop_names{ii}); eval(cmdTxt);
end

avar = exec_fun_on_cell_mat(all_cellsnew,'percent');
% avar = exec_fun_on_cell_mat(all_cellsnew_zvals,'mean',all_cellsnew);
fac_names = {'PoT','MT','BT','CN','AP'}; fac_levels = [length(pop_names),2,2,3,2];
% fac_names = {'MT','CN','AP'}; fac_levels = [2,3,2];
[within,dvn,xlabels,awithinD] = make_within_table(fac_names,fac_levels);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
clc
print_for_manuscript(ra)
%%
clc
raS = RMA_subset(ra,'PoT');

%% Plot rasters for one animal top-down concatenated time-distance-and speed cells
% out = outT; bn = 1;
out = outD; bn = 2;
for an = 1%:5
    for cn = 3%:3
        for ap = 1
            rasters = out.arasters{an,cn,ap}; rx = out.arx{an,cn,ap}; ry = out.ary{an,cn,ap};
            tCells = time_cells{an,cn,ap,bn};
            dCells = distance_cells{an,cn,ap,bn};
            sCells = speed_cells{an,cn,ap,bn};
            num_neurons = size(rasters,3);
            for nn = 1:num_neurons
                if tCells(nn) == 0
                    continue;
                end
                figure(100);clf;
                imagesc(rx,ry,rasters(:,:,nn));colorbar;
                title(nn)
                pause(0.3);
            end
        end
    end
end

%%
%% Find time cells with PC in air_on phase and plot graphs and look at population stats
% this is with time-binning
% RV is response variable
tt = 1;
clear pvals zvals mi_zvals tCells
RV = 'perc';
% RV = 'rfid';
% RV = 'mval';
% RV = 'MI';
% RV = 'PC';
avar = [];
for an = 1:5
    anvar = [];
    for cn = 1:3
        ap = 1;
        AC = outT.active_cells_trials{an,cn,ap}; mAC1 = outT.active_cells{an,cn,ap}; RF = sum(AC,2); mAC = RF>5;
        out = outT; MV = met_valsT;         % out = outD; MV = met_valsD; 
        selcells = mAC1;% & (RF>0);
        idx = 3; pvals = [MV{tt}{an,cn,ap}.PC(:,idx)]; idx = 2; zvals = [MV{tt}{an,cn,ap}.PC(:,idx)];
        idx = 2; mi_zvals = [MV{tt}{an,cn,ap}.MI(:,idx)];
        % idx = 3; pvals = [MV{1}{an,cn,ap}.MI(:,idx)];
        % idx = 2; zvals = [MV{1}{an,cn,ap}.MI(:,idx)];
        % criteria for finding the types of cells
        tCells = pvals < 0.05; %time_cells{an,cn,ap,bn} = tCells< 0.05;
        if strcmp(RV,'perc')
            anvar = [anvar sum(tCells & selcells)/length(AC)];
        end
        if strcmp(RV,'mval')
            anvar = [anvar nanmean(zvals(tCells&selcells,1))];
        end
        if strcmp(RV,'rfid')
            anvar = [anvar mean(RF(tCells&selcells,1))];
        end
        if strcmp(RV,'MI')
            anvar = [anvar mean(mi_zvals(tCells&selcells,1))];
        end
        if sum(isnan(anvar)) > 0
            [an cn]
        end
    end
    if strcmp(RV,'perc')
        avar = [avar;100*anvar];
    else
        avar = [avar;anvar];
    end
end

fac_names = {'CN'}; fac_levels = [3];
% fac_names = {'CN','AP','BT','RF','CT'}; fac_levels = [3,2,2,3,4];
% fac_names = {'CN','AP','BT','CT'}; fac_levels = [3,2,2,7];
[within,dvn,xlabels,awithinD] = make_within_table(fac_names,fac_levels);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
clc
print_for_manuscript(ra)
% (Intercept) [F(1,4) = 93.20, p < 0.001, η2 = .90] <--
%% visualizing the results in the previous section
tcolors = repmat(mData.dcolors(8:10),1,3); MY = 30; ysp = 1; mY = 0; ystf = 2; ysigf = 0.025;titletxt = ''; ylabeltxt = {'Cells (%)'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova([],ra,{'CN','hsd',0.05},[1 1.75],tcolors,[mY MY ysp ystf ysigf],mData);