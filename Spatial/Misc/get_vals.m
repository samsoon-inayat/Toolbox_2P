function valso = get_vals(vals,resp)

for rr = 1:size(vals,1)
    for cc = 1:size(vals,2)
        valso{rr,cc} = vals{rr,cc}(resp{rr,cc},:);
    end
end