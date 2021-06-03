function [resp_fraction,resp_vals,all_OI,mean_OI] = get_responsive_fraction(Rs)


for rr = 1:size(Rs,1)
    ccs = [];
    for cc = 1:size(Rs,2)
        R = Rs{rr,cc};
        ccs(:,cc) = R.resp.vals';
        resp_fraction(rr,cc) = R.resp.fraction;
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

all_OIm = repmat(OI,1,1,size(ccs,2));
for ii = 1:length(resp_vals)
    OI = all_OI{ii};
    all_OIm(:,:,ii) = OI;
end

mean_OI = nanmean(all_OIm,3);