function resp = sep_cell_list(resp_valsC,resp_valsC1)

for rr = 1:size(resp_valsC,1)
    for cc = 1:size(resp_valsC,2)
        c = resp_valsC{rr,cc}; c1 = resp_valsC1{rr,cc};
        resp{rr,cc} = c & ~c1;
    end
end
