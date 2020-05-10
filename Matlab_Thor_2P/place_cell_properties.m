function props = place_cell_properties (mSig,ccSignalA,varargin)

p = inputParser;
default_cm_per_bin = 1.5;
default_min_pc_width = 5;
default_max_pc_width = 120;

addRequired(p,'mSig',@isnumeric);
addRequired(p,'ccSignalA',@isnumeric);
addOptional(p,'cm_per_bin',default_cm_per_bin,@isnumeric);
addOptional(p,'max_pc_width',default_max_pc_width,@isnumeric);
addOptional(p,'min_pc_width',default_min_pc_width,@isnumeric);
parse(p,mSig,ccSignalA,varargin{:});
cm_per_bin = p.Results.cm_per_bin;
min_pc_width = p.Results.min_pc_width;
max_pc_width = p.Results.max_pc_width;

props.width = NaN;
props.center = NaN;

initial_threshold = 0.2*(max(mSig) - min(mSig));%+min(mSig);
idxs = find(mSig > initial_threshold);
iidxs = diff(idxs);
gOnes = find(iidxs>1);
if ~isempty(gOnes)
    if gOnes(1) ~= 1
        gOnes = [1 gOnes];
    end
    if gOnes(end) ~= length(iidxs)
        gOnes = [gOnes length(iidxs)];
    end

    d_gOnes = diff(gOnes)*cm_per_bin;
    if sum(d_gOnes>min_pc_width & d_gOnes<max_pc_width) == 0
        return;
    end
    longest_chunk = find(d_gOnes>min_pc_width & d_gOnes<max_pc_width);
    if length(longest_chunk) > 1
        longest_chunk = find(d_gOnes == max(d_gOnes(longest_chunk)));
    end
    idxs_of_longest_chunk = idxs((gOnes(longest_chunk)+2):gOnes(longest_chunk+1));
else
    idxs_of_longest_chunk = idxs;
end

% PC = length(idxs_of_longest_chunk)*cm_per_bin;
% 
if length(idxs_of_longest_chunk)*cm_per_bin > max_pc_width
    return;
end
if length(idxs_of_longest_chunk)*cm_per_bin < min_pc_width
    return;
end
mean_activity_inside_place_field = mean(mSig(idxs_of_longest_chunk));
mean_activity_outside_place_field = mean(mSig(mSig < initial_threshold));
if mean_activity_inside_place_field < (3*mean_activity_outside_place_field)
    return;
end
[~,idx_of_peak_of_trials] = max(ccSignalA,[],2);
count = 0;
for kk = 1:length(idx_of_peak_of_trials)
    count = count + sum(idxs_of_longest_chunk == idx_of_peak_of_trials(kk));
end
if size(ccSignalA,2)/2 < count
    return;
end
% figure(10002);plot(1:50,mSig);
props.width = length(idxs_of_longest_chunk)*cm_per_bin;
[~,temp] = max(mSig(idxs_of_longest_chunk));
props.center = idxs_of_longest_chunk(temp)* cm_per_bin;