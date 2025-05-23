function resp = sep_cell_list(resp_valsCi,resp_valsC1i)



if isstruct(resp_valsCi)
    resp_valsC = resp_valsCi.vals;
    resp_valsC1 = resp_valsC1i.vals;
else
    resp_valsC = resp_valsCi;
    resp_valsC1 = resp_valsC1i;
end

for rr = 1:size(resp_valsC,1)
    for cc = 1:size(resp_valsC,2)
        c = resp_valsC{rr,cc}; c1 = resp_valsC1{rr,cc};
        resp{rr,cc} = c & ~c1;
    end
end
