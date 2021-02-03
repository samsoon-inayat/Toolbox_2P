function n_clus_inds = getClusters(data)
propMatrix = [data.MI' data.Ent' data.HaFD' data.HiFD'];
num_clus = 20;
[clus_inds,clus_centers] = kmeans(propMatrix,num_clus);
%     data.cluster = clus_inds;
for ii = 1:num_clus
    cluster_median_MI(ii) = median(data.MI(clus_inds == ii));
    cluster_mean_MI(ii) = mean(data.MI(clus_inds == ii));
    cluster_min_MI(ii) = min(data.MI(clus_inds == ii));
end
MI_th = nanmean(data.MI) + 2*nanstd(data.MI);
clus_More = find(cluster_min_MI >= MI_th);
n_clus_inds = zeros(size(clus_inds));
for ii = 1:length(clus_More)
    n_clus_inds = n_clus_inds | clus_inds == clus_More(ii);
end

