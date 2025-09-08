function G = jaccard_group(idv)
% Collect per-animal Jaccard maps and compute group summaries (diag=NaN).

A = size(idv,1);  C = size(idv,2);

D_mean_all = nan(C,C,A);
D_sem_all  = nan(C,C,A);
N_pair_all = nan(C,C,A);

for a = 1:A
    [~, Dm, Ds, Nij] = jaccard_animal(idv,a);
    D_mean_all(:,:,a) = Dm;        % already has NaN on diagonal
    D_sem_all(:,:,a)  = Ds;        % already has NaN on diagonal
    N_pair_all(:,:,a) = Nij;
end

% Group mean/SEM across animals (omit NaNs)
D_group_mean = mean(D_mean_all, 3, 'omitnan');
na = sum(~isnan(D_mean_all),3);
D_group_sem  = std(D_mean_all, 0, 3, 'omitnan') ./ sqrt(max(na,1));

% Agreement summaries (same diagonal NaNs)
P_group_mean = 1 - D_group_mean;
P_group_sem  = D_group_sem;
% P_group_mean(1:7:end) = NaN;
% D_group_sem(1:7:end) = NaN;

G.D_mean_all   = D_mean_all;
G.D_sem_all    = D_sem_all;
G.N_pair_all   = N_pair_all;
G.D_group_mean = D_group_mean;
G.D_group_sem  = D_group_sem;
G.P_group_mean = P_group_mean;
G.P_group_sem  = P_group_sem;
end
