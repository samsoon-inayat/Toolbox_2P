function oedge = find_first_falling_edge(signal,threshold,minimum_diff)
edge = find(diff(signal) <= threshold);
temp = find(diff(edge)<=minimum_diff);
edge(temp+1) = [];

if ~isempty(edge)
    oedge = edge(1);
else
    oedge = length(signal);
end