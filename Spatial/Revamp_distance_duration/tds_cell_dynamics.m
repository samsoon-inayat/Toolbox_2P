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

%%
MVT = met_valsT; MVD = met_valsD;
for an = 1:5
    for cni = 1:size(si_cn_ap,2)
        cn = si_cn_ap(1,cni); ap = si_cn_ap(2,cni);
        idx = 3; pvalsPC = [MVT{1}{an,cn,ap}.PC(:,idx) MVD{2}{an,cn,ap}.PC(:,idx) MVT{3}{an,cn,ap}.PC(:,idx)];
        idx = 3; pvalsMI = [MVT{1}{an,cn,ap}.MI(:,idx) MVD{2}{an,cn,ap}.MI(:,idx) MVT{3}{an,cn,ap}.MI(:,idx)];
        idvPC{an,cni} = pvalsPC < 0.05;
        idvMI{an,cni} = pvalsMI < 0.05;
        pvalsPC_a{an,cni} = pvalsPC;
        pvalsMI_a{an,cni} = pvalsMI;
        idx = 2; zvalsPC = [MVT{1}{an,cn,ap}.PC(:,idx) MVD{2}{an,cn,ap}.PC(:,idx) MVT{3}{an,cn,ap}.PC(:,idx)];
        idx = 2; zvalsMI = [MVT{1}{an,cn,ap}.MI(:,idx) MVD{2}{an,cn,ap}.MI(:,idx) MVT{3}{an,cn,ap}.MI(:,idx)];
        zvalsPC_a{an,cni} = zvalsPC;
        zvalsMI_a{an,cni} = zvalsMI;
    end
end

idv = idvPC;
pvals = pvalsPC_a; zvals = zvalsPC_a;
disp('Done')
%%
%% Set input (you already have these in workspace)
% idv = idvPC;              % or: idv = idvMI;

%% --- Pick one animal for detailed view
a = 1;
[D_cell, D_mean, D_sem, N_pair, P_same] = jaccard_animal(idv, a);
adj_idx = sub2ind([6 6], 1:5, 2:6);
mean_adj = mean(D_mean(adj_idx), 'omitnan');
fprintf('Animal %d: mean adjacent identity distance = %.3f\n', a, mean_adj);

% Plots: per-animal mean & SEM
figure('Color','w');
subplot(1,3,1); imagesc(D_mean,[0 1]); axis square; colormap(parula); colorbar
title(sprintf('Jaccard distance (mean over cells) — A%d',a));
subplot(1,3,2); imagesc(P_same,[0 1]); axis square; colormap(parula); colorbar
title('Agreement = 1 - distance');
subplot(1,3,3); imagesc(D_sem,[0 0.3]); axis square; colormap(parula); colorbar
title('SEM over cells');

%% --- Group Jaccard stats
G = jaccard_group(idv);
figure('Color','w');
subplot(1,2,1); imagesc(G.D_group_mean,[0 1]); axis square; colormap(parula); colorbar
title('Group Jaccard distance (mean over animals)');
subplot(1,2,2); imagesc(G.P_group_mean,[0 1]); axis square; colormap(parula); colorbar
title('Group agreement (=1 - distance)');

%% --- Stability percentages (Stable / Dynamic / Middle per label)
P = stability_percentages(idv);
labels = {'Time','Distance','Speed'};
fprintf('\nPer-animal %%Stable / %%Dynamic / %%Middle (denominators = cells that ever had the label)\n');
for a = 1:size(idv,1)
    fprintf('Animal %d:\n', a);
    for b = 1:3
        fprintf('  %-8s  S=%.1f%%  D=%.1f%%  M=%.1f%%   (n=%d)\n', ...
            labels{b}, P.stable(a,b), P.dynamic(a,b), P.middle(a,b), P.denom(a,b));
    end
end

fprintf('\nGroup means ± SEM (%%):\n');
for b = 1:3
    fprintf('  %-8s  Stable  %5.1f ± %4.1f   Dynamic %5.1f ± %4.1f   Middle %5.1f ± %4.1f\n', ...
        labels{b}, ...
        P.mean_sem.stable_mean(b),  P.mean_sem.stable_sem(b), ...
        P.mean_sem.dynamic_mean(b), P.mean_sem.dynamic_sem(b), ...
        P.mean_sem.middle_mean(b),  P.mean_sem.middle_sem(b));
end

%% (Optional) save results
% save results_idv.mat D_cell D_mean D_sem N_pair P_same G P

