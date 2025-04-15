function [rasters,x,y] = build_rasters(trials,FR)


num_neurons = size(FR,1);
trialidxnums = accumarray(trials,1,[],@sum);
rasters = NaN(length(trialidxnums),max(trialidxnums),num_neurons);

for nn = 1:size(FR,1)
    tFR = FR(nn,:);
    for ii = 1:length(trialidxnums)
        lenT = sum(trials == ii);
        rasters(ii,1:lenT,nn) = tFR(trials == ii);
    end
end

x = 1:max(trialidxnums);
y = 1:length(trialidxnums);
