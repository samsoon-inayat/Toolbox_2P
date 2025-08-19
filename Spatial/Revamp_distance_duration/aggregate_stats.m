function M = aggregate_stats(statsCell, field)
% statsCell: 1×N_anim cell, each has S.field (18×18)
C = size(statsCell{1}.(field),1);
N = numel(statsCell);
stack = NaN(C,C,N);
for i=1:N, stack(:,:,i) = statsCell{i}.(field); end
M.mean = nanmean(stack,3);
M.sem  = nanstd(stack,[],3) ./ sqrt(N);
end
