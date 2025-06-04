function [uvs,yvs] = trial_collections(resps)

uvs = 1:size(resps,2);
for rr = 1:size(resps,1)
    yv = find_trial_collections(resps(rr,:));
    yvs(rr,:) = yv;
end
n = 0;


function yv = find_trial_collections(resp)
uvs = 1:size(resp,2);
yv = ones(size(uvs));

n = 0;
dresp = diff(resp);
for ii = 1:length(uvs)
    
end
