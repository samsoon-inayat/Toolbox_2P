function G = jaccard_group(idv)
% Loop animals, collect per-animal Jaccard maps, and compute group mean/SEM.
% Returns struct G with fields:
%   .D_mean_all  : 6x6xA per-animal means
%   .D_sem_all   : 6x6xA per-animal SEMs (over cells)
%   .D_group_mean: 6x6 mean over animals
%   .D_group_sem : 6x6 SEM over animals
%   .P_group_mean: 6x6 group agreement (=1 - D_group_mean)
%   .N_pair_all  : 6x6xA per-animal contributing counts

A = size(idv,1);  C = size(idv,2);  assert(C==6);

D_mean_all = nan(6,6,A);
D_sem_all  = nan(6,6,A);
N_pair_all = nan(6,6,A);

for a = 1:A
    [~, Dm, Ds, Nij] = jaccard_animal(idv,a);
    D_mean_all(:,:,a) = Dm;
    D_sem_all(:,:,a)  = Ds;
    N_pair_all(:,:,a) = Nij;
end

D_group_mean = mean(D_mean_all, 3, 'omitnan');
na = sum(~isnan(D_mean_all),3);
D_group_sem  = std(D_mean_all, 0, 3, 'omitnan') ./ sqrt(max(na,1));

G.D_mean_all   = D_mean_all;
G.D_sem_all    = D_sem_all;
G.D_group_mean = D_group_mean;
G.D_group_sem  = D_group_sem;
G.P_group_mean = 1 - D_group_mean;
G.N_pair_all   = N_pair_all;
end
