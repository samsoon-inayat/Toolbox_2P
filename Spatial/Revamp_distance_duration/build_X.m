function X = build_X(all_cells, animal_idx, cols18)
% BUILD_X  Return cells×18 logical matrix for one animal and lens.
% Columns are ordered as 6 configs × [T D S].
nCells = numel(all_cells{animal_idx, cols18(1)});
X = false(nCells, numel(cols18));
for j = 1:numel(cols18)
    v = all_cells{animal_idx, cols18(j)};
    X(:,j) = logical(v(:));
    if numel(v) ~= nCells
        error('Inconsistent cell count for animal %d', animal_idx);
    end
end
end
