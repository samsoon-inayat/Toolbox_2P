function bidishifts = get_bidishift(T)
fileName = 'bidishift.mat';
for ii = 1:size(T,1)
    pd_path = T{ii,7};
    fn = fullfile(pd_path,fileName);
    temp = load(fn{1});
    bidishifts(ii) = temp.bidishift;
end