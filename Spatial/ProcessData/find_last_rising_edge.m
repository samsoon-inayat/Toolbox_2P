function edge = find_first_rising_edge(signal,threshold,minimum_diff)
edge = find(diff(signal) >= threshold);
temp = find(diff(edge)<=minimum_diff);
edge(temp+1) = [];

if ~isempty(edge)
    edge = egde(end);
end
