function [is_stable, is_dynamic, thrQ] = stability_classify(S_bit)
% Quartile-based classification per label (columns 1..3).
% Returns logical masks (nCells x 3) and thresholds thrQ = [Q1; Q3] (2x3).

nCells = size(S_bit,1);
is_stable  = false(nCells,3);
is_dynamic = false(nCells,3);
thrQ = nan(2,3);

for b = 1:3
    v = S_bit(:,b);
    mask = ~isnan(v);
    if ~any(mask), continue; end
    q = quantile(v(mask), [0.25 0.75]);
    thrQ(:,b) = q(:);
    is_stable(mask,b)  = v(mask) >= q(2);
    is_dynamic(mask,b) = v(mask) <= q(1);
end
end
