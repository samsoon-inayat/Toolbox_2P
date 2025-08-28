function S_bit = stability_per_label(idv, animal_idx)
% Compute per-cell, per-label retention: P(retain next | present now).
% Returns S_bit: [nCells x 3], columns = [T D S].

C = size(idv,2);  assert(C==6);
nCells = size(idv{animal_idx,1},1);

% Stack bits
bits = false(nCells,3,C);
for c = 1:C
    x = logical(idv{animal_idx,c});
    if size(x,1)~=nCells || size(x,2)~=3, error('Size mismatch'); end
    bits(:,:,c) = x;
end

S_bit = nan(nCells,3);
for r = 1:nCells
    for b = 1:3
        starts = 0; stays = 0;
        for t = 1:C-1
            if bits(r,b,t)==1
                starts = starts + 1;
                if bits(r,b,t+1)==1, stays = stays + 1; end
            end
        end
        if starts>0
            S_bit(r,b) = stays/starts;
        end
    end
end
end
