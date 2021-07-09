function h = get_obj(ff,obj)
h = [];
if isstruct(ff)
    has = ff.h_axes;
else
    has = ff;
end

for ii = 1:size(has,1)
    for jj = 1:size(has,2)
        h(ii,jj) = get(has(ii,jj),obj);
    end
end