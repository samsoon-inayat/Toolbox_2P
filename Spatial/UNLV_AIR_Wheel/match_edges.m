function M = match_edges(t_ref, t_mov, maxAbsDiff)
% Match each ref time to nearest moving time, within tolerance
% Returns struct with pairs + diffs

t_ref = t_ref(:);
t_mov = t_mov(:);

idx_mov = nan(size(t_ref));
dt      = nan(size(t_ref));

for i = 1:numel(t_ref)
    [dmin, j] = min(abs(t_mov - t_ref(i)));
    if dmin <= maxAbsDiff
        idx_mov(i) = j;
        dt(i) = t_mov(j) - t_ref(i);
    end
end

M.idx_mov = idx_mov;
M.dt      = dt;
M.n_matched = sum(~isnan(dt));
end
