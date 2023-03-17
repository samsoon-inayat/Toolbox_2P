function valso = get_vals(vals,resp)
if size(vals,2) == size(resp,2)
    for rr = 1:size(vals,1)
        for cc = 1:size(vals,2)
            valso{rr,cc} = vals{rr,cc}(logical(resp{rr,cc}),:);
        end
    end
else
    b = size(resp,2)/size(vals,2);
    for rr = 1:size(vals,1)
        ind = 1;
        for cc = 1:size(vals,2)
            for bb = 1:b
                valso{rr,ind} = vals{rr,cc}(resp{rr,ind},:);
                ind = ind + 1;
            end
        end
    end
end