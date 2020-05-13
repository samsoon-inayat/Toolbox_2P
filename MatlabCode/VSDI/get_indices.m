function indices = get_indices (stim_folders_info,ss,gg,aa,sc)
indxs = logical(ones(size(stim_folders_info.list,1),1));

if ~isempty(ss)
    indxs = indxs & stim_folders_info.list(:,1) == ss;
end

if ~isempty(gg)
    indxs = indxs & stim_folders_info.list(:,2) == gg;
end

if ~isempty(aa)
    indxs = indxs & stim_folders_info.list(:,3) == aa;
end



indices = find(indxs);

if ~isempty(sc)
    list = [];
    for ii = 1:length(indices)
        if ~isempty(strfind(stim_folders_info.root_folders{indices(ii)},sc))
            list = [list ii];
        end
    end
    indices(list) = [];
end

