function [all_OI] = get_overlap_index(resp_valsi)

% resp_valsCi = resp_valsC;
for rr = 1:size(resp_valsi,1)
    ccs = [];
    for cc = 1:size(resp_valsi,2)
        ccs(:,cc) = resp_valsi{rr,cc};
    end
    resp_vals{rr} = ccs;
end


for ii = 1:length(resp_vals)
    ccs = resp_vals{ii};
    OI = NaN(size(ccs,2),size(ccs,2));
    mask = triu(ones(size(ccs,2)),1);
    for rr = 1:size(ccs,2)
        ccs1 = ccs(:,rr);
        for cc = 1:size(ccs,2)
            if mask(rr,cc)
                ccs2 = ccs(:,cc);
                shared = ccs1 & ccs2;
                OI(rr,cc) = sum(shared)/(sum(ccs1)+sum(ccs2)-sum(shared));
            end
        end
    end
    all_OI{ii} = OI;
end